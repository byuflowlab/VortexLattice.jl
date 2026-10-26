# Sectional polars are read through `sectional_coefficients(polar, alpha, Re, Mach)`, the
# interface FLOWUnsteadyCore defines, so the polars LiftingLines takes work here too. Until
# VortexLattice depends on the core, the `CCBlade.AlphaAF` method is duplicated here.

"""
    sectional_coefficients(polar, alpha, Re, Mach) -> (cl, cd)

Lift and drag coefficient of a section polar at angle of attack `alpha` [rad]. For a
`CCBlade.AlphaAF` (what `CCBlade.AlphaAF(file)` or `FLOWUnsteadyCore.read_polar`
return) this forwards to `CCBlade.afeval`; another polar type adds a method.
"""
@inline sectional_coefficients(af::CCBlade.AlphaAF, alpha, Re, Mach) =
    CCBlade.afeval(af, alpha, Re, Mach)

"""
    viscous!(props, Γ, surface, grid, polars, basis, ref)

Correct one surface's panel forces and circulation with sectional polars, after the
inviscid solve: each spanwise strip's lift is moved to the polar's cl at the strip's
effective angle (additively, see below) and its cd is added along the local flow.

## Arguments

- `props`: the surface's panel properties (nc x ns), updated in place
- `Γ`: the surface's panel circulations (column-major, nc x ns), updated in place
- `surface`: the surface panels
- `grid`: the surface grid (chord lengths)
- `polars`: one polar per spanwise strip (anything with a [`sectional_coefficients`](@ref)
  method, e.g. a `CCBlade.AlphaAF`). Reynolds and Mach number are not computed and are
  passed as zero, so polars that correct for them are not supported yet
- `basis`: the surface's axes in global coordinates (columns: chordwise, spanwise,
  normal); the strips' lift and drag directions lie in its x-z plane
- `ref`: reference quantities used for dimensional scaling
"""
function viscous!(props::AbstractMatrix{PanelProperties{TF}}, Γ, surface::AbstractMatrix{SurfacePanel{TF}}, grid,
        polar_array::AbstractVector, Rp::SMatrix{3,3}, ref::Reference) where TF

    # index for circulation strengths
    iΓ = 1

    # dynamic pressure for dimensionalizing forces
    q = 0.5 * RHO * ref.V * ref.V

    # rotation from the global frame to the surface's axes
    R = transpose(Rp)

    # loop over spanwise sections in this surface
    for j in axes(surface, 2)

        # get the circulation per unit length of this section
        γ = zero(TF)

        # get induced velocity
        v_induced = zero(SVector{3,TF})

        # get net aerodynamic force on this section
        cf = zero(SVector{3,TF})

        # extract the polar for this section
        polar = polar_array[j]

        # initialize dynamic pressure
        q_local = zero(TF)

        # number of chordwise panels for this surface
        nc_section = size(surface, 1)

        # loop over chordwise panels in this section
        for i in axes(surface, 1)

            # project bound vortex onto y axis of this frame
            ds = R * (surface[i,j].rtr - surface[i,j].rtl)
            dy = abs(ds[2]) / norm(ds)

            # Γ[iΓ] is the CUMULATIVE bound circulation from the leading edge
            # through panel i (see the cumulative-sum reconstruction below), so the
            # section's 2D-equivalent total circulation is just the trailing-edge
            # (last chordwise) panel's Γ -- not a sum over all nc panels, which
            # over-counts at nc>1 (Γ[1]+Γ[2]+...+Γ[nc] double/triple/quadruple-counts
            # since each Γ[i] already includes every prior segment). Harmless at
            # nc=1 (single panel, trivially the trailing-edge value) but corrupted
            # cl_vlm/Δcl at nc>1, producing nonphysical CT (2026-08-12 fix).
            if i == nc_section
                γ = Γ[iΓ] * dy
            end
            iΓ += 1

            # accumulate induced velocity contribution at this bound vortex
            v_induced += props[i,j].velocity * ref.V # convert from non-dimensionalized velocity

            # accumulate aerodynamic force contribution from this bound vortex
            cf += props[i,j].cfb
        end

        # average dynamic pressure over the section
        v_induced /= size(surface, 1) # average induced velocity over the section
        v_local = norm(v_induced) # local velocity magnitude at this section
        q_local = 0.5 * RHO * v_local * v_local

        # get chord length
        le = SVector(grid[1,1,j], grid[2,1,j], grid[3,1,j])
        te = SVector(grid[1,end,j], grid[2,end,j], grid[3,end,j])
        c = norm(le - te)

        # average of chord lengths of either side
        le = SVector(grid[1,1,j+1], grid[2,1,j+1], grid[3,1,j+1])
        te = SVector(grid[1,end,j+1], grid[2,end,j+1], grid[3,end,j+1])
        c += norm(le - te)
        c *= 0.5

        # convert summed force coefficient back to dimensional force
        cf *= q * ref.S

        # project aerodynamic force into xz plane
        cf = R * cf # rotate into this frame

        # get 2-D lift magnitude of this section
        l_2d_norm = sqrt(cf[1]*cf[1] + cf[3]*cf[3])

        # force per length
        l_2d_norm /= abs(( R * (surface[1,j].rtl - surface[1,j].rtr) )[2])

        # calculate effective cl predicted by the VLM, = 2π * α_eff
        cl_vlm = -2 * RHO * γ * (γ / (l_2d_norm + eps(l_2d_norm))) / c * sign(γ)

        # get effective α
        α_eff = cl_vlm / (2 * pi) # in radians

        # additive viscous correction: the polar's cl at the effective α minus the
        # inviscid (thin-airfoil) cl there, which is cl_vlm itself
        cl_visc, cd = sectional_coefficients(polar, α_eff, zero(TF), zero(TF))
        Δcl = cl_visc - cl_vlm

        # lift direction in the strip xz plane: perpendicular to local
        # flow projected into xz. well-defined whenever there is any
        # freestream over the section.
        v_induced_strip = R * v_induced

        dhat_strip = SVector(v_induced_strip[1], 0.0, v_induced_strip[3])
        dhat_strip /= norm(dhat_strip)
        lhat_strip = SVector(-dhat_strip[3], 0.0, dhat_strip[1])

        # lift direction in the global frame (used to project the
        # globally-stored cfb / Δs / V vectors below)
        lhat = Rp * lhat_strip

        # spanwise extent of this strip (used to convert Δcl → ΔL)
        Δs_y = ( R * (surface[1,j].rtl - surface[1,j].rtr) )[2]

        # strip-level prescribed dimensional lift increment
        nc = size(surface, 1)
        Δl_strip = Δcl * q_local * c * Δs_y

        # viscous drag: dimensionalized the same way as Δl_strip above (cd * q_local *
        # c * Δs_y), NOT via the γ-based l_2d_norm/γ reconstruction. That reconstruction
        # divides by γ² with no epsilon guard, so it produced Inf/NaN whenever γ crossed
        # zero (e.g. Ct~990 instead of ~0.77) -- see prior disable note. q_local comes
        # from the total local velocity (props.velocity), independent of γ, so it stays
        # well-behaved through zero-circulation stations. Re-enabled 2026-08-12.
        # The strip width is |Δs_y|: Δs_y (left minus right) is negative on a grid laid out
        # along +y. Δl_strip above may keep the sign, because cl_vlm is signed the same
        # way and the two cancel; the drag has no second sign, and with Δs_y it pushed
        # upstream (fixed 2026-10-26).
        D_visc_strip = cd * q_local * c * abs(Δs_y)
        # cfb (added to below) is a force COEFFICIENT, non-dimensionalized by q*ref.S
        # (see the (RHO*Δγ_j)*cross_jl/(q*ref.S) term a few lines down) -- D_visc_strip
        # is a dimensional force, so it must go through the same normalization before
        # being added, or it inflates the reported force by a factor of q*ref.S (~1e5-1e6
        # for this rotor). Missing normalization caused CT~1770 instead of ~0.9 when this
        # was first re-enabled 2026-08-12.
        d_viscous = (D_visc_strip / nc) * (Rp * dhat_strip) / (q * ref.S)

        # ---- per-segment additive correction ------------------------
        # The strip's prescribed dimensional lift increment Δl_strip is
        # distributed equally across the nc chordwise bound segments.
        # Each segment's bound circulation changes by the Kutta-Joukowski
        # amount
        #
        #     Δγ_j = (Δl_strip / nc) / (ρ (V_j × Δs_j) · l̂)
        #
        # using the local effective velocity V_j and bound-vortex segment
        # Δs_j that nearfield.jl uses, and its force by ρ Δγ_j (V_j × Δs_j).
        # Additive, so the unsteady dΓ/dt force in cfb stays force and is
        # not converted into circulation (rebuilding γ from the total
        # segment force did that, 2026-10-26); in steady flow it equals
        # solving l_j + Δl/nc = ρ proj_j γ_j_new.
        #
        # Δγ_j changes the *segment* circulation Γ_b,i = Γ[i] - Γ[i-1], so
        # panel i's circulation changes by the cumulative sum Σ_{k≤i} Δγ_k.
        Δγ_accum = zero(TF)
        for i in axes(surface, 1)
            # unpack props
            (; gamma, velocity, cfb, cfl, cfr) = props[i,j]

            # local effective velocity at the bound vortex midpoint
            # (props.velocity is stored as Vi/ref.V; see nearfield.jl:554)
            Vj = velocity * ref.V

            # bound vortex segment vector (global frame)
            Δsj = top_vector(surface[i,j])

            # KJ projection onto the lift direction
            cross_jl = cross(Vj, Δsj)
            proj_j = dot(cross_jl, lhat)

            # change of this segment's circulation from local KJ
            Δγ_j = (Δl_strip / nc) / (RHO * proj_j)

            # its bound-vortex force (nearfield.jl's F_b = ρ Γ_b (V × Δs));
            # cfl/cfr untouched
            cfb_new = cfb + (RHO * Δγ_j) * cross_jl / (q * ref.S) + d_viscous

            props[i,j] = PanelProperties(
                gamma,
                velocity,
                cfb_new,  # apply lift correction factor to bound circulation contribution and add viscous drag
                cfl, # * f_cl,  # apply lift correction factor to left edge contribution
                cfr, # * f_cl,  # apply lift correction factor to right edge contribution
            )

            # panel circulation: Γ[i] += Σ_{k≤i} Δγ_k
            Δγ_accum += Δγ_j
            Γ[iΓ - nc + (i - 1)] += Δγ_accum
        end
    end
    return nothing
end
