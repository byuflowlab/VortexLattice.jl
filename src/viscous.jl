# --- Diagnostic instrumentation: capture per-section effective AoA (deg) computed by the ---
# --- viscous correction below, keyed by surface index, one Vector per call to viscous!.  ---
# --- Added 2026-07-29 to compare local AoA against a CCBlade reference; does not affect   ---
# --- forces/circulation. Enable with `ALPHA_LOG_ENABLED[] = true`; read out `ALPHA_LOG`.  ---
const VISC_ITER_DEBUG = Ref(false)
const ALPHA_LOG_ENABLED = Ref(false)
const ALPHA_LOG = Dict{Int, Vector{Vector{Float64}}}()
# Local strip-frame velocity components (v[1], v[3] -- the xz in-plane components used to
# build dhat_strip below), one (v1,v3) pair per section per call. Added 2026-08-10 alongside
# the alpha bias investigation, to decompose VL's local relative wind into axial/tangential-like
# components for comparison against CCBlade's induction factors (a, ap). Same enable flag/reset
# as ALPHA_LOG since both are populated together in the same loop.
const VELOCITY_LOG = Dict{Int, Vector{Vector{Tuple{Float64,Float64}}}}()
# Corrected lift coefficient (cl_vlm + Δcl, i.e. the polar-informed cl actually used to
# build Δl_strip) and the RAW polar-lookup drag coefficient BEFORE the `*0.0` disable
# (see the viscous-drag comment below) -- added 2026-08-11 to compare VL's own Cl/Cd
# against a CCBlade polar-lookup reference, same enable flag/history convention as
# ALPHA_LOG/VELOCITY_LOG. Logging the raw (undisabled) Cd does not re-enable it --
# forces/circulation are unaffected, this is read-only instrumentation.
const CL_LOG = Dict{Int, Vector{Vector{Float64}}}()
const CD_LOG = Dict{Int, Vector{Vector{Float64}}}()
const CLDIRECT_LOG = Dict{Int, Vector{Vector{Float64}}}()
function reset_alpha_log!()
    empty!(ALPHA_LOG)
    empty!(VELOCITY_LOG)
    empty!(CL_LOG)
    empty!(CD_LOG)
    empty!(CLDIRECT_LOG)
    return nothing
end

struct Polar{TF}
    alphas::Vector{TF}
    cls_inv::Vector{TF}
    cls_delta::Vector{TF}
    cls_visc::Vector{TF}
    cds_visc::Vector{TF}
end

"""
    Polar(alphas, cls_visc, cds_visc)

Construct a `Polar` from a viscous polar table. Inviscid lift is taken as the
thin-airfoil approximation `cl_inv = 2π·α` (with `alphas` in degrees), and the
additive correction table is precomputed as `cls_delta = cls_visc - cls_inv`.
This allows the viscous correction to be applied as a smooth additive lookup
indexed by inviscid `cl`, rather than as a ratio `cl_visc / cl_inv`.
"""
function Polar(alphas, cls_visc, cds_visc)
    TF = eltype(cls_visc)
    cls_inv = TF(2π) .* (alphas .* TF(π/180))  # thin-airfoil: cl_inv = 2π·α[rad]
    cls_delta = cls_visc .- cls_inv
    return Polar{TF}(alphas, cls_inv, cls_delta, cls_visc, cds_visc)
end

function Polar{TF}(alphas, cls_visc, cds_visc) where TF
    cls_inv = TF(2π) .* (alphas .* TF(π/180))  # thin-airfoil: cl_inv = 2π·α[rad]
    cls_delta = cls_visc .- cls_inv
    return Polar{TF}(TF.(alphas), cls_inv, cls_delta, TF.(cls_visc), TF.(cds_visc))
end

"""
    get_polars2(section_rs, rotor_file, TF=Float64; data_path)

Generates a vector of vectors of `::Polar` objects for use in the `viscous!` function.

**Arguments**
- `section_rs::Vector{Vector{Float64}}`: `section_rs[i]` contains a vector of radial coordinates of each airfoil in the `i`th surface, normalized by the radius. Note that a 0-length vector indicates no airfoil correction is to be used.
- `rotor_file::String`: name of the file containing the radial coordinates and polar files

**Optional Arguments**
- `data_path::String`: path where `rotors/` and `airfoils/` directories live

**Returns**
- `polars::Vector{Vector{Polar{Float64}}}`: polar object for each section in each surface of the corresponding `::System`
"""
function get_polars2(section_rs::Vector{<:Vector}, rotor_files, TF=Float64; data_path="VortexLattice_rotor_data")
    
    # check vector lengths
    @assert length(section_rs) == length(rotor_files)
    
    # create vector of polars to push to
    polars = Vector{Vector{Polar{TF}}}(undef, length(section_rs))

    # loop over surfaces
    for i_surface in eachindex(section_rs)

        section_r = section_rs[i_surface]
        if length(section_r) > 0

            # instantiate polar vector
            sections = Vector{Polar{TF}}(undef, length(section_r))

            # get this rotor file
            rotor_file = rotor_files[i_surface]

            # read rotor data
            data = readdlm(joinpath(data_path, "rotors", rotor_file), ',', skipstart=1)

            # generate polar objects
            polars_list = Vector{Polar{TF}}(undef, size(data,1))
            af_files = String.(data[:,3])
            for i_polar in axes(data,1)
                # get file name
                af_file = af_files[i_polar]

                # read file
                af_data = readdlm(joinpath(data_path, "airfoils", af_file), ',', skipstart=1)

                # generate polar
                polars_list[i_polar] = Polar{TF}(af_data[:,1], af_data[:,2], af_data[:,3])
            end

            # calculate ranges for each airfoil
            rRs = Float64.(data[:,1])
            rRmids = (rRs[1:end-1] .+ rRs[2:end]) .* 0.5
            rRranges = zeros(length(rRs))
            rRranges[1:end-1] .= rRmids
            rRranges[end] = 1.0

            # loop over sections
            for i_r in eachindex(section_r)
                
                # radial coordinate of this section
                r = section_r[i_r]

                # find which airfoil corresponds
                i_polar = findfirst((x) -> (x>=r), rRranges)

                # populate sections
                sections[i_r] = polars_list[i_polar]

            end

            polars[i_surface] = sections
        
        else

            # push 0-length vector
            polars[i_surface] = Vector{Polar{TF}}(undef,0)
        end
    end

    return polars
end

"""
    viscous!(properties, Γ, dΓdt, surfaces, grids, frames, frames_index, polars, ref, dt)

Apply viscous corrections to the aerodynamic forces and circulation strengths
using sectional airfoil polar data. The corrections are applied to each
surface's panel `properties`, to the circulation strengths `Γ`, and to their
time derivatives `dΓdt`.

If `polars === nothing`, this function returns without modifying its inputs.

Note: `dΓdt` should contain `-Γ` from the previous timestep.

## Arguments

- `properties`: panel aerodynamic properties to be updated in place
- `Γ`: panel circulation strengths, updated in place
- `dΓdt`: circulation time derivative workspace, updated in place
- `surfaces`: surface panel geometry
- `grids`: grid coordinates for each surface
- `frames`: reference frames for each surface
- `frames_index`: map from surfaces to entries in `frames`
- `polars`: sectional viscous polar data for each surface, or `nothing` to skip
  viscous corrections
- `ref`: aerodynamic reference quantities used for dimensional scaling
- `dt`: timestep used to convert the updated `Γ` values into `dΓdt`

## Returns

Mutates `properties`, `Γ`, and `dΓdt` in place and returns `nothing`.
"""
function viscous!(properties::Vector{Matrix{PanelProperties{TF}}}, Γ, surfaces::Vector{Matrix{SurfacePanel{TF}}}, grids, frames::Vector{<:ReferenceFrame}, frames_index::Vector{Int}, polars::Vector{<:Vector{<:Polar}}, ref::Reference, dt) where TF
    # properties contains:
    # * cf
    # surface contains:
    # * position
    # grid contains:
    # * chord information
    # Γ contains:
    # * circulation
    # frame contains:
    # * rotation matrix into this frame

    # index for circulation strengths
    iΓ = 1

    # dynamic pressure for dimensionalizing forces
    q = 0.5 * RHO * ref.V * ref.V

    # diagnostic AoA capture for this call (see ALPHA_LOG above)
    _alpha_this_call = ALPHA_LOG_ENABLED[] ? Dict{Int, Vector{Float64}}() : nothing
    _velocity_this_call = ALPHA_LOG_ENABLED[] ? Dict{Int, Vector{Tuple{Float64,Float64}}}() : nothing
    _cl_this_call = ALPHA_LOG_ENABLED[] ? Dict{Int, Vector{Float64}}() : nothing
    _cd_this_call = ALPHA_LOG_ENABLED[] ? Dict{Int, Vector{Float64}}() : nothing
    # Diagnostic (2026-08-17): l_2d_norm/(q_local*c), i.e. the sectional cl you'd get from the
    # REAL local dynamic pressure (q_local, built from the actual |v_induced|) instead of from
    # cl_vlm's self-referential KJ formula (-2*RHO*gamma^2/(c*l_2d_norm), which implicitly uses
    # an inferred "V_KJ" that need not equal the real |v_induced| -- see the 2026-08-14
    # v_local-vs-V_KJ finding). Comparing this against cl_vlm tests whether the slow multi-
    # iteration convergence is caused by that mismatch (relevant here since this is a 45-deg
    # swept wing, where the bound-vortex segment isn't orthogonal to the local velocity).
    _cldirect_this_call = ALPHA_LOG_ENABLED[] ? Dict{Int, Vector{Float64}}() : nothing

    # loop over surfaces
    for isurf in eachindex(surfaces)

        if frames_index[isurf] > 0

            # extract containers corresponding to this surface
            surface = surfaces[isurf]
            props = properties[isurf]
            grid = grids[isurf]
            polar_array = polars[isurf]

            # get reference frame for this surface
            frame = frames[frames_index[isurf]]

            # get rotation matrix from global frame to this frame
            Rp = frame.Rp2g * frame.R # rotation from this frame to global frame
            R = transpose(Rp)

            # @show R # verified: basis vectors of the frame expressed in global coordinates

            if !isnothing(_alpha_this_call)
                _alpha_this_call[isurf] = Float64[]
                _velocity_this_call[isurf] = Tuple{Float64,Float64}[]
                _cl_this_call[isurf] = Float64[]
                _cd_this_call[isurf] = Float64[]
                _cldirect_this_call[isurf] = Float64[]
            end

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

                # Effective sectional cl, from the REAL local dynamic pressure (q_local, built
                # from |v_induced|) -- not from the self-referential KJ formula
                # (-2*RHO*gamma^2/(c*l_2d_norm)) this used to use. That formula implicitly
                # assumes the bound-vortex segment is perpendicular to the local velocity
                # (valid for an unswept wing) via the scalar relation L'=RHO*V*Gamma; the real
                # relation is the vector cross product F=RHO*Gamma*(V x Ds), which picks up a
                # sin(theta) factor whenever the segment isn't orthogonal to V (i.e. any swept
                # wing). Confirmed via the Weber & Brebner 45-deg swept wing (2026-08-17): the
                # old cl_vlm and this cl_direct disagree by up to ~0.35 and iterating the old
                # scheme crawls (13,000+ iterations) to a degenerate near-zero-Cl fixed point
                # trying to reconcile the two; cl_direct has no such inconsistency to resolve.
                cl_vlm = sign(γ) * l_2d_norm / (q_local * c + eps(l_2d_norm))

                # get effective α (thin-airfoil inverse, same approximation as before, now
                # applied to the physically-correct cl instead of the swept-biased one)
                α_eff = cl_vlm / (2 * pi) * 180 / pi # in degrees

                if !isnothing(_alpha_this_call)
                    push!(_alpha_this_call[isurf], α_eff)
                end

                if !isnothing(_cldirect_this_call)
                    push!(_cldirect_this_call[isurf], cl_vlm)
                end

                # additive viscous correction: Δcl(cl_inv), indexed by the
                # inviscid cl so the table stays smooth near zero lift
                Δcl = FLOWMath.linear(polar.cls_inv, polar.cls_delta, cl_vlm)

                if !isnothing(_cl_this_call)
                    push!(_cl_this_call[isurf], cl_vlm + Δcl)
                end

                # lift direction in the strip xz plane: perpendicular to local
                # flow projected into xz. well-defined whenever there is any
                # freestream over the section.
                v_induced_strip = R * v_induced

                if !isnothing(_velocity_this_call)
                    push!(_velocity_this_call[isurf], (v_induced_strip[1], v_induced_strip[3]))
                end

                dhat_strip = SVector(v_induced_strip[1], 0.0, v_induced_strip[3])
                dhat_strip /= norm(dhat_strip)
                lhat_strip = SVector(-dhat_strip[3], 0.0, dhat_strip[1])

                # lift direction in the global frame (used to project the
                # globally-stored cfb / Δs / V vectors below)
                lhat = Rp * lhat_strip

                # spanwise extent of this strip (used to convert Δcl → ΔL). This is a SIGNED
                # coordinate difference (kept signed here -- see the abs()'d `Δs_y_mag` below for
                # why the sign must NOT be dropped for the lift term specifically).
                Δs_y = ( R * (surface[1,j].rtl - surface[1,j].rtr) )[2]

                # magnitude-only version for the DRAG term below: on a mirrored wing Δs_y comes
                # out negative across the ENTIRE span (confirmed: -0.06223 at every j on the Weber
                # & Brebner 45deg swept-wing case, -0.3 at every j on the unswept NACA4415 case --
                # never flips sign either way), so using the signed Δs_y for D_visc_strip made a
                # positive polar cd produce a NEGATIVE dimensional drag force -- the "viscous
                # drag" term SUBTRACTED drag instead of adding it (confirmed: Weber-wing total CD
                # went from +0.00344 inviscid to -0.00232 with the bug). Confirmed via the real
                # package call that using abs() ONLY here (not also on Δl_strip below) fixes that
                # (CD -> +0.0092, physically correct: viscous = inviscid + profile drag, never
                # less than inviscid) while leaving the NACA4415/Cocco lift correction unchanged
                # (CL=1.1279, matching the pre-fix/memory-validated value) -- using abs() on BOTH
                # terms (tried first) collapsed NACA4415's CL to 0.28, so the lift term's sign
                # convention is NOT the same bug and must be left alone. Confirmed 2026-08-19,
                # viscous_drag_sign_fix_prototype.jl (VPM-Validation repo).
                Δs_y_mag = abs(Δs_y)

                # strip-level prescribed dimensional lift increment
                nc = size(surface, 1)
                Δl_strip = Δcl * q_local * c * Δs_y

                # viscous drag: dimensionalized the same way as Δl_strip above (cd * q_local *
                # c * Δs_y), NOT via the γ-based l_2d_norm/γ reconstruction. That reconstruction
                # divides by γ² with no epsilon guard, so it produced Inf/NaN whenever γ crossed
                # zero (e.g. Ct~990 instead of ~0.77) -- see prior disable note. q_local comes
                # from the total local velocity (props.velocity), independent of γ, so it stays
                # well-behaved through zero-circulation stations. Re-enabled 2026-08-12.
                cd_raw = FLOWMath.linear(polar.alphas, polar.cds_visc, α_eff)
                if !isnothing(_cd_this_call)
                    push!(_cd_this_call[isurf], cd_raw)
                end
                cd = cd_raw
                D_visc_strip = cd * q_local * c * Δs_y_mag
                # cfb (added to below) is a force COEFFICIENT, non-dimensionalized by q*ref.S
                # (see the (RHO*γ_j_new)*cross_jl/(q*ref.S) term a few lines down) -- D_visc_strip
                # is a dimensional force, so it must go through the same normalization before
                # being added, or it inflates the reported force by a factor of q*ref.S (~1e5-1e6
                # for this rotor). Missing normalization caused CT~1770 instead of ~0.9 when this
                # was first re-enabled 2026-08-12.
                d_viscous = (D_visc_strip / nc) * (Rp * dhat_strip) / (q * ref.S)

                # ---- per-segment additive correction ------------------------
                # The strip's prescribed dimensional lift increment Δl_strip is
                # distributed equally across the nc chordwise bound segments.
                # For each segment we solve the local Kutta-Joukowski relation
                #
                #     l_j_new = ρ ((V_j × Δs_j) · l̂) γ_j_new
                #
                # for the new bound-vortex circulation γ_j_new, using the local
                # effective velocity V_j and bound-vortex segment Δs_j that
                # nearfield.jl uses to compute panel forces. This guarantees
                # that the downstream cfb is exactly consistent with γ_j_new.
                #
                # γ_j here is the *segment* (bound-vortex) circulation
                # Γ_b,i = Γ[i] - Γ[i-1] (with Γ_b,1 = Γ[1]). After solving for
                # all γ_j_new, panel circulations are recovered by cumulative
                # sum and written back into Γ.
                Γ_b_new_accum = zero(TF)
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

                    # existing inviscid lift on this segment (dimensional)
                    l_j = dot(cfb, lhat) * q * ref.S

                    # new segment circulation from local KJ
                    γ_j_new = (l_j + Δl_strip / nc) / (RHO * proj_j)

                    # recompute the bound-vortex contribution to cfb directly
                    # from γ_j_new so cfb stays exactly consistent with the
                    # corrected circulation (matches nearfield.jl form
                    # F_b = ρ Γ_b (V × Δs)). preserves cfl/cfr untouched.
                    cfb_new = (RHO * γ_j_new) * cross_jl / (q * ref.S) + d_viscous

                    props[i,j] = PanelProperties(
                        gamma,
                        velocity,
                        cfb_new,  # apply lift correction factor to bound circulation contribution and add viscous drag
                        cfl, # * f_cl,  # apply lift correction factor to left edge contribution
                        cfr, # * f_cl,  # apply lift correction factor to right edge contribution
                    )

                    # reconstruct panel circulation by cumulative sum of
                    # segment circulations: Γ[i] = Σ_{k≤i} γ_k_new
                    Γ_b_new_accum += γ_j_new
                    Γ[iΓ - nc + (i - 1)] = Γ_b_new_accum
                end
            end
        else
            iΓ += size(surfaces[isurf], 1) * size(surfaces[isurf], 2) # skip circulation strengths for this surface
        end
    end

    if !isnothing(_alpha_this_call)
        for (isurf, row) in _alpha_this_call
            push!(get!(ALPHA_LOG, isurf, Vector{Float64}[]), row)
        end
        for (isurf, row) in _velocity_this_call
            push!(get!(VELOCITY_LOG, isurf, Vector{Tuple{Float64,Float64}}[]), row)
        end
        for (isurf, row) in _cl_this_call
            push!(get!(CL_LOG, isurf, Vector{Float64}[]), row)
        end
        for (isurf, row) in _cd_this_call
            push!(get!(CD_LOG, isurf, Vector{Float64}[]), row)
        end
        for (isurf, row) in _cldirect_this_call
            push!(get!(CLDIRECT_LOG, isurf, Vector{Float64}[]), row)
        end
    end
end

function viscous!(properties::Vector{Matrix{PanelProperties{TF}}}, Γ, dΓdt, surfaces::Vector{Matrix{SurfacePanel{TF}}}, grids, frames::Vector{<:ReferenceFrame}, frames_index::Vector{Int}, polar::Nothing, ref, dt) where TF
    return nothing
end

"""
    viscous_iterative!(properties, Γ, surfaces, wakes, grids, frames_index,
        polars, ref, fs; symmetric, nwake, surface_id, wake_finite_core,
        wake_shedding_locations, trailing_vortices, xhat, Vh=nothing, Vv=nothing,
        dΓdt=nothing, maxiter=200, tol=1e-6, m=4, beta=0.3)

Alpha-target Anderson-accelerated alternative to `viscous!`. Instead of `viscous!`'s single
additive `Δcl(cl_inv)` correction applied once to the frozen inviscid `Γ` (a blind accumulator
with no genuine fixed point when iterated -- see
`vortexlattice_viscous_correction_diverges_confirmed_2026_08_17` /
`vortexlattice_native_viscous_accumulator_bug_diagnosed_20260818` memory), this reassigns each
viscous section's circulation directly to its target `Γ_target = 0.5 * cl_polar(α_local) * c *
|V_local|` and iterates that map to an actual fixed point via Anderson(m) acceleration,
refreshing induced velocities with a real `near_field_forces!` call (`calculate_vlm_induced =
true`, i.e. the O(N²) direct VLM Biot-Savart pass, NOT whatever FMM-precomputed `Vh`/`Vv` term
the outer unsteady loop may be using -- each trial `Γ` during the internal iteration needs its
own self-consistent induced velocity, which a cached FMM evaluation from the ORIGINAL Γ cannot
provide) once per internal iteration -- unlike the additive scheme, this map is provably
convergent under damped/accelerated iteration away from stall (see
`vortexlattice_native_viscous_target_fix_prototyped_tested_20260819` memory).

Hard-capped at `maxiter` iterations (default 200, no fallback needed beyond that): if the cap is
hit, the partially-converged `Γ` from the last Anderson step is used as-is. This is expected only
within a few degrees of a local extremum in the polar's `dCl/dα` (e.g. right at post-stall
rollover, confirmed structurally as a locally-expanding fixed-point map there, not fixable by any
damping constant -- see `vortexlattice_intelligent_iteration_methods_tested_20260819` memory) and
is a reasonable fallback since the outer unsteady time-marching re-seeds from a fresh (and
typically nearby) inviscid `Γ` next step regardless (see
`vortexlattice_warmstart_assumption_corrected_20260819` memory -- there is no persisted viscous
`Γ` across timesteps in production to lose by falling back here).

**`beta` default lowered to 0.3 (not full-step `beta=1.0`) after testing on a SWEPT wing**: the
unswept NACA4415 validation converges cleanly at `beta=1.0`, but the 45°-swept Weber & Brebner
case develops a spanwise checkerboard (alternating-sign cl error from station to station) at
`beta=1.0`/`m=4` -- the same family of mirror-pair/cross-station-coupling instability documented
for Cocco's per-station accelerators on swept geometry (see
`vortexlattice_naca4415_smoothing_removed_secant_scheme_20260818` memory), here showing up in a
plain Anderson-accelerated joint fixed-point map rather than a per-station one. Confirmed the fix
is just damping harder: `beta=0.3` (with `m=4`) or plain damped Picard at `beta<=0.1` both
converge cleanly to the same answer as `viscous!`/inviscid on the swept case, and `beta=0.3` costs
no accuracy on the unswept case either (same CL to <0.02 away from the fold, `m=4` and 200
iterations is generous headroom for the swept case's slower convergence at this damping).

**Restricted to `nc == 1`** (one chordwise panel per section) for every surface with polar data;
raises an `ArgumentError` otherwise. This mirrors the only configuration validated so far.
Multiple surfaces ARE supported and iterated jointly (a single flat Anderson vector spanning
every viscous station across every surface, so cross-surface induction is captured exactly, not
approximated).

The local (chordwise, normal) frame used to compute each station's effective angle of attack is
derived directly from the panel/grid GEOMETRY (chordwise = leading-to-trailing-edge direction,
normal = chordwise × spanwise), not from `frames`/`ReferenceFrame` rotations -- `ReferenceFrame`
encodes a vehicle body-axis CONVENTION (e.g. `Rp2g` can include an axis flip such as
`BackRightUp` vs `ForwardRightDown`) that is unrelated to the airfoil's actual chord direction;
using it to split velocity into "chordwise"/"normal" components silently corrupted the effective
angle of attack (confirmed via a debug case: `Rp2g = diag(-1,1,1)` flipped an 8° local alpha into
~173°). The additive `viscous!` never hit this because it derives `cl` from force magnitude +
`sign(γ)`, never from a raw velocity angle -- deriving the angle directly from geometry instead
of any external frame convention is the general-case-correct fix, valid for any sweep/rotation.

**Performance note**: because each internal iteration needs a full induced-velocity refresh at
`calculate_vlm_induced = true`, the per-call cost scales like the direct O(N²) VLM Biot-Savart
sum times the number of iterations (up to `maxiter`) -- validated cheap (a few ms/step) at
~20-90 viscous stations; for very large panel counts (where the outer solve relies on FMM
specifically to avoid an O(N²) cost) this internal iteration could become the dominant per-
timestep cost and should be profiled before use at that scale.

**`calculate_vlm_induced`/`additional_velocity` overrides (added 2026-08-19)**: default
`calculate_vlm_induced=true` recomputes each station's local velocity directly from
`surfaces`/`wakes` only -- correct for the steady `WakePanel` case, where `wakes` genuinely holds
the entire wake. For the unsteady `PanelParticleWake` path, `wakes` is just the thin near-wake
panel buffer; the bulk of the wake is vortex particles (`wake.pfield`) that never enter that direct
sum, so a post-hoc call against an unsteady system silently sees only a ~1-chord-wake local flow
regardless of how long the simulated wake actually is (confirmed: flat CL across a 1-16 chord
sweep -- see `vortexlattice_viscous_iterative_steady_vs_unsteady_gap_found_20260819` memory).
`system.Vh`/`system.Vv` do NOT help here -- they hold velocity due to rigid-body *vehicle motion*
only (zero for a non-maneuvering vehicle in a fixed freestream), not wake/particle induction;
setting `calculate_vlm_induced=false` with those fields just drops induction entirely (confirmed
worse: CL becomes bit-identical across every wake length, since the local AoA collapses to the
bare freestream angle). The correct fix is `additional_velocity` (already an unconditionally-added
term in `near_field_forces!`, independent of `calculate_vlm_induced`): keep
`calculate_vlm_induced=true` (still recomputes bound self + near-wake-row induction, which DOES
depend on the trial Γ each Anderson iteration) and pass an `additional_velocity(rc)` closure
supplying the particle wake's induced velocity at point `rc`. Since the particles are a frozen
background field during this post-hoc call (only bound Γ changes iteration to iteration), that
contribution can be evaluated ONCE per call (e.g. via one `FastMultipole.fmm!` pass of a
`ProbeSystem` against `wake.pfield` alone, keeping bound self-induction out of the source set to
avoid double-counting what `calculate_vlm_induced=true` already supplies) and reused every
iteration.
"""
function viscous_iterative!(properties, Γ, surfaces, wakes, grids, frames_index::Vector{Int},
        polars, ref::Reference, fs; symmetric, nwake, surface_id,
        wake_finite_core, wake_shedding_locations, trailing_vortices, xhat, Vh=nothing, Vv=nothing,
        dΓdt=nothing, maxiter=200, tol=1e-6, m=4, beta=0.3, calculate_vlm_induced=true,
        additional_velocity=nothing)

    TF = eltype(Γ)
    nsurf = length(surfaces)

    # ---- flat list of viscous stations (one per section) across every surface. `station_iΓ`
    # holds the TRAILING-EDGE (cumulative-total) panel index for each station -- for nc==1 this
    # is trivially the (only) panel's own index, preserving all prior nc==1 behavior exactly;
    # for nc>1, `station_iΓ_first`/`station_nc` give the rest of that station's chordwise column
    # (see station_target's Γ-weighted velocity average and write_station!'s proportional
    # rescale below, both ported from naca4415_cocco_validation_test.jl's validated nc>1 handling
    # -- vortexlattice_viscous_nc_gt1_support_added_20260820 memory).
    iΓ = 1
    station_iΓ = Int[]
    station_iΓ_first = Int[]
    station_nc = Int[]
    station_isurf = Int[]
    station_j = Int[]
    for isurf in 1:nsurf
        surface = surfaces[isurf]
        n_i, n_j = size(surface)
        if frames_index[isurf] > 0 && !isempty(polars[isurf])
            for j in 1:n_j
                push!(station_iΓ_first, iΓ)
                push!(station_iΓ, iΓ + n_i - 1)
                push!(station_nc, n_i)
                push!(station_isurf, isurf)
                push!(station_j, j)
                iΓ += n_i
            end
        else
            iΓ += n_i * n_j
        end
    end
    n = length(station_iΓ)
    n == 0 && return nothing

    # section chord lengths + local (chordwise, normal) unit vectors, straight from panel/grid
    # geometry in the SAME global frame `properties[...].velocity` is expressed in (see docstring
    # note on why this must NOT go through `frames`/`ReferenceFrame` rotations)
    c_arr = zeros(TF, n)
    that_arr = Vector{SVector{3,TF}}(undef, n)
    nhat_arr = Vector{SVector{3,TF}}(undef, n)
    for k in 1:n
        isurf = station_isurf[k]; j = station_j[k]
        grid = grids[isurf]
        le1 = SVector(grid[1,1,j], grid[2,1,j], grid[3,1,j])
        te1 = SVector(grid[1,end,j], grid[2,end,j], grid[3,end,j])
        le2 = SVector(grid[1,1,j+1], grid[2,1,j+1], grid[3,1,j+1])
        te2 = SVector(grid[1,end,j+1], grid[2,end,j+1], grid[3,end,j+1])
        c_arr[k] = 0.5 * (norm(le1 - te1) + norm(le2 - te2))
        t_hat = normalize((te1 - le1) + (te2 - le2))  # chordwise, leading- to trailing-edge
        s_hat = normalize(le2 - le1)                  # spanwise
        that_arr[k] = t_hat
        nhat_arr[k] = normalize(cross(t_hat, s_hat))  # airfoil "up" (suction-side) normal
    end

    function refresh_forces!()
        near_field_forces!(properties, surfaces, wakes, ref, fs, Γ; dΓdt=dΓdt,
            additional_velocity=additional_velocity, Vh=Vh, Vv=Vv, symmetric, nwake, surface_id,
            wake_finite_core, wake_shedding_locations, trailing_vortices, xhat,
            calculate_vlm_induced=calculate_vlm_induced)
        return nothing
    end

    # local (chordwise, normal) flow angle at this station, and the polar's cl-target
    # circulation implied by it: Γ_target = 0.5 * cl(α) * c * |V_local|
    function station_target(k)
        isurf = station_isurf[k]; j = station_j[k]
        nc_k = station_nc[k]
        if nc_k == 1
            v_induced = properties[isurf][1, j].velocity * ref.V
        else
            # Γ-difference-weighted average velocity across chordwise rows: each row's
            # bound-vortex FORCE contribution is weighted by its own net circulation
            # Γ[i]-Γ[i-1], not a uniform 1/nc share -- see naca4415_cocco_validation_test.jl's
            # cocco_targets for the original derivation/validation of this exact formula.
            first = station_iΓ_first[k]
            v_induced = zero(SVector{3,TF})
            wsum = zero(TF)
            Γprev = zero(TF)
            for i in 1:nc_k
                Γcur = Γ[first + i - 1]
                Γrow_net = Γcur - Γprev
                v_induced += Γrow_net * properties[isurf][i, j].velocity * ref.V
                wsum += Γrow_net
                Γprev = Γcur
            end
            if wsum == 0
                # zero-circulation seed: Γ-weighted sum degenerates to 0/0 -- fall back to a
                # plain unweighted average so the first iteration can bootstrap away from Γ=0.
                v_induced = zero(SVector{3,TF})
                for i in 1:nc_k
                    v_induced += properties[isurf][i, j].velocity * ref.V
                end
                v_induced /= nc_k
            else
                v_induced /= wsum
            end
        end
        u2D_t = dot(v_induced, that_arr[k])
        u2D_n = dot(v_induced, nhat_arr[k])
        u2D_mag = sqrt(u2D_t^2 + u2D_n^2)
        α = atan(u2D_n, u2D_t) * 180 / pi
        polar = polars[isurf][j]
        cl_target = FLOWMath.linear(polar.alphas, polar.cls_visc, α)
        return 0.5 * cl_target * c_arr[k] * u2D_mag, α, u2D_mag
    end

    # writes a new TOTAL (trailing-edge/cumulative) station circulation, redistributing the
    # change across the station's chordwise column by uniformly rescaling every row's existing
    # NET circulation -- preserves whatever chordwise loading shape the inviscid/prior-iteration
    # solution already has while correcting the total to match the viscous target, without
    # introducing a frame-dependent per-segment KJ solve (which would reintroduce the same class
    # of ReferenceFrame-convention risk found in viscous!()'s per-segment reconstruction -- see
    # vortexlattice_viscous_frame_mismatch_root_cause_20260820 memory). Exact for nc==1.
    function write_station!(k, new_total)
        nc_k = station_nc[k]
        if nc_k == 1
            Γ[station_iΓ[k]] = new_total
        else
            first = station_iΓ_first[k]
            old_total = Γ[station_iΓ[k]]
            if old_total == 0
                for i in 1:nc_k
                    Γ[first + i - 1] = new_total * i / nc_k
                end
            else
                r = new_total / old_total
                for i in 1:nc_k
                    Γ[first + i - 1] *= r
                end
            end
        end
        return nothing
    end

    # ---- Anderson(m) acceleration of the Picard map g(x) = x + R(x), seeded from the
    # incoming Γ (the fresh inviscid Γ this call actually receives -- see
    # vortexlattice_warmstart_assumption_corrected_20260819 memory) ----
    x = [Γ[station_iΓ[k]] for k in 1:n]
    Xs = Vector{Vector{TF}}()
    Fs = Vector{Vector{TF}}()
    refresh_forces!()
    for it in 0:maxiter-1
        f = zeros(TF, n)
        for k in 1:n
            Γ_target, _, _ = station_target(k)
            f[k] = Γ_target - x[k]
        end
        resid = maximum(abs.(f))
        Γscale = maximum(abs.(x)) + eps(TF)
        if VISC_ITER_DEBUG[]
            println("  [visc_iter] it=$it resid/scale=$(resid/Γscale) Γ_target[1]=$(station_target(1)[1]) alpha[1]=$(station_target(1)[2]) x[1]=$(x[1])")
        end
        if resid / Γscale < tol
            break
        end
        push!(Xs, copy(x)); push!(Fs, copy(f))
        mk = min(m, length(Xs) - 1)
        if mk == 0
            x = x .+ beta .* f
        else
            ΔF = hcat([Fs[end-i+1] .- Fs[end-i] for i in 1:mk]...)
            ΔX = hcat([Xs[end-i+1] .- Xs[end-i] for i in 1:mk]...)
            γ = ΔF \ f
            x = x .+ beta .* f .- (ΔX .+ beta .* ΔF) * γ
        end
        for k in 1:n
            write_station!(k, x[k])
        end
        refresh_forces!()
    end
    for k in 1:n
        write_station!(k, x[k])
    end
    refresh_forces!()  # final refresh so `properties` is consistent with the returned Γ

    # ---- viscous drag: same convention as viscous!() (cd*q_local*c*Δs_y added to cfb), all in
    # the global frame (matches how cfb/Δs are stored -- see docstring note on `frames`).
    # Distributed as D_visc_strip/nc across every chordwise row (matches viscous!()'s per-segment
    # d_viscous convention) rather than dumped entirely onto row 1.
    q = 0.5 * RHO * ref.V * ref.V
    for k in 1:n
        isurf = station_isurf[k]; j = station_j[k]; nc_k = station_nc[k]
        surface = surfaces[isurf]
        _, α, u2D_mag = station_target(k)
        polar = polars[isurf][j]
        cd = FLOWMath.linear(polar.alphas, polar.cds_visc, α)
        q_local = 0.5 * RHO * u2D_mag^2
        c = c_arr[k]
        Δs_y = norm(surface[1,j].rtl - surface[1,j].rtr)
        v_induced = properties[isurf][1, j].velocity * ref.V
        u2D_t = dot(v_induced, that_arr[k])
        u2D_n = dot(v_induced, nhat_arr[k])
        dhat_global = (u2D_t * that_arr[k] + u2D_n * nhat_arr[k]) / u2D_mag
        D_visc_strip = cd * q_local * c * Δs_y
        d_viscous = (D_visc_strip / nc_k) * dhat_global / (q * ref.S)
        for i in 1:nc_k
            (; gamma, velocity, cfb, cfl, cfr) = properties[isurf][i, j]
            properties[isurf][i, j] = PanelProperties(gamma, velocity, cfb + d_viscous, cfl, cfr)
        end
    end

    return nothing
end

"""
    viscous_iterative_shed!(system, wake::PanelParticleWake, trailing_edge_filaments,
        Γ_wake, polars, ref, fs, frames_index; kwargs...)

Per-step, fully Anderson-converged viscous coupling for the `PanelParticleWake` unsteady path,
intended to REPLACE the single-pass `viscous!()` call in `_simulate_step!` when the shed wake's
own strength needs to track the viscous-corrected circulation (enable via
`simulate!(...; viscous_iterative_shed=true)`).

**Why this exists** (2026-08-20): `viscous_iterative!()` applied post-hoc (after `simulate!`
finishes, on an otherwise purely-inviscid-shed wake) plateaus far from the steady target,
independent of wake length -- root-caused to the near-wake buffer/particle wake having been
shed with the INVISCID Γ history the whole run, so its own downwash never reflects the
viscous-corrected circulation (see
vortexlattice_viscous_iterative_frozen_wake_selfconsistency_ROOT_CAUSE_20260820 memory). The
fix is to converge the viscous circulation BEFORE shedding, every step -- not just once at the
end -- exactly like single-pass `viscous!()` already does (it mutates `Γ_wake` in place before
`shed_wake!` uses it), except with a full self-consistent Anderson solve instead of one additive
increment (which is known to diverge when repeated every step, see
vortexlattice_viscous_correction_diverges_confirmed_2026_08_17).

**How the trial-Γ-dependent velocity is derived** (avoiding the ~2x error from the direct
`calculate_vlm_induced=true` Biot-Savart recompute, see
vortexlattice_viscous_iterative_induced_velocity_mismatch_20260820): the WAKE-on-vehicle
contribution to `system.Vh`/`Vv` (populated earlier this step by `wake_on_all!`, before the AIC
solve) is Γ-independent -- frozen once as a baseline. Only the VEHICLE-on-vehicle contribution
(bound self + TE-to-wake_shedding_location interface ring) depends on the trial Γ, so
`vehicle_on_all!` is re-run each Anderson iteration with the trial Γ written into `system.Γ`,
added on top of the frozen baseline, then `near_field_forces!` runs with
`calculate_vlm_induced=false` -- the IDENTICAL code path production uses, just re-evaluated at
each trial Γ instead of only once.

On return, `Γ_wake` holds the converged viscous circulation (for `shed_wake!`) and
`system.properties`/`system.Vh`/`system.Vv` are left consistent with it (mirroring how
`viscous!()` also leaves `system.properties` viscous-corrected). `system.Γ` is restored to its
original (inviscid) value, since the next step's AIC solve overwrites it anyway but leaving it
viscous-corrected mid-step would be surprising to any code reading it before then.
"""
function viscous_iterative_shed!(system, wake, trailing_edge_filaments, Γ_wake, polars,
        ref, fs, frames_index; maxiter=200, tol=1e-6, m=4, beta=0.3,
        fmm_wake_args=NamedTuple(), fmm_vehicle_args=NamedTuple())

    nsurf = length(system.surfaces)
    TF = eltype(system.Γ)

    # ---- flat list of viscous stations (see viscous_iterative!'s equivalent block for the
    # nc>1 generalization -- station_iΓ is the TRAILING-EDGE/cumulative-total panel index) ----
    iΓ = 1
    station_iΓ = Int[]; station_iΓ_first = Int[]; station_nc = Int[]
    station_isurf = Int[]; station_j = Int[]
    for isurf in 1:nsurf
        n_i, n_j = size(system.surfaces[isurf])
        if frames_index[isurf] > 0 && !isempty(polars[isurf])
            for j in 1:n_j
                push!(station_iΓ_first, iΓ)
                push!(station_iΓ, iΓ + n_i - 1)
                push!(station_nc, n_i)
                push!(station_isurf, isurf); push!(station_j, j)
                iΓ += n_i
            end
        else
            iΓ += n_i * n_j
        end
    end
    n = length(station_iΓ)
    n == 0 && return nothing

    grids_ = system.grids
    c_arr = zeros(TF, n)
    that_arr = Vector{SVector{3,TF}}(undef, n)
    nhat_arr = Vector{SVector{3,TF}}(undef, n)
    for k in 1:n
        isurf = station_isurf[k]; j = station_j[k]; grid = grids_[isurf]
        le1 = SVector(grid[1,1,j], grid[2,1,j], grid[3,1,j]); te1 = SVector(grid[1,end,j], grid[2,end,j], grid[3,end,j])
        le2 = SVector(grid[1,1,j+1], grid[2,1,j+1], grid[3,1,j+1]); te2 = SVector(grid[1,end,j+1], grid[2,end,j+1], grid[3,end,j+1])
        c_arr[k] = 0.5 * (norm(le1 - te1) + norm(le2 - te2))
        t_hat = normalize((te1 - le1) + (te2 - le2))
        s_hat = normalize(le2 - le1)
        that_arr[k] = t_hat; nhat_arr[k] = normalize(cross(t_hat, s_hat))
    end

    Γ_seed = copy(system.Γ) # the inviscid solve this step already produced; restored at the end

    # ---- frozen wake-on-vehicle baseline (Γ-independent) ----
    # system.Vh/Vv at the call site (_simulate_step!, right after the INVISCID near_field_forces!
    # call) already hold kinematic + wake_on_all!'s contribution + vehicle_on_all!'s contribution
    # AT THE INVISCID Γ (vehicle_on_all! already ran earlier this step, before this function is
    # called) -- NOT just the wake-only piece. Reusing that directly would double-count bound
    # self-induction (stale inviscid + fresh trial, added again by vehicle_on_all! below every
    # iteration). Recompute a clean wake-only baseline explicitly instead.
    for isurf in 1:nsurf
        system.Vh[isurf] .= Ref(zero(eltype(system.Vh[isurf])))
        system.Vv[isurf] .= Ref(zero(eltype(system.Vv[isurf])))
        system.Vcp[isurf] .= Ref(zero(eltype(system.Vcp[isurf])))
        system.Vte[isurf] .= Ref(zero(eltype(system.Vte[isurf])))
    end
    wake_on_all!(system, wake, trailing_edge_filaments; fmm_wake_args...)
    Vh_base = deepcopy(system.Vh); Vv_base = deepcopy(system.Vv)

    nsurf_range = 1:nsurf
    symmetric_f = fill(false, nsurf)
    trailing_vortices_f = fill(false, nsurf)

    function refresh_forces!()
        for isurf in 1:nsurf
            system.Vh[isurf] .= Vh_base[isurf]
            system.Vv[isurf] .= Vv_base[isurf]
        end
        vehicle_on_all!(system, wake, trailing_edge_filaments; fmm_vehicle_args...)
        near_field_forces!(system.properties, system.surfaces, system.wakes, ref, fs, system.Γ;
            dΓdt=nothing, additional_velocity=nothing, Vh=system.Vh, Vv=system.Vv,
            symmetric=symmetric_f, nwake=system.nwake, surface_id=nsurf_range,
            wake_finite_core=fill(true, nsurf), wake_shedding_locations=nothing,
            trailing_vortices=trailing_vortices_f, xhat=system.xhat[], calculate_vlm_induced=false)
    end

    function station_target(k)
        isurf = station_isurf[k]; j = station_j[k]
        nc_k = station_nc[k]
        if nc_k == 1
            v_induced = system.properties[isurf][1, j].velocity * ref.V
        else
            first = station_iΓ_first[k]
            v_induced = zero(SVector{3,TF})
            wsum = zero(TF)
            Γprev = zero(TF)
            for i in 1:nc_k
                Γcur = system.Γ[first + i - 1]
                Γrow_net = Γcur - Γprev
                v_induced += Γrow_net * system.properties[isurf][i, j].velocity * ref.V
                wsum += Γrow_net
                Γprev = Γcur
            end
            if wsum == 0
                v_induced = zero(SVector{3,TF})
                for i in 1:nc_k
                    v_induced += system.properties[isurf][i, j].velocity * ref.V
                end
                v_induced /= nc_k
            else
                v_induced /= wsum
            end
        end
        u2D_t = dot(v_induced, that_arr[k]); u2D_n = dot(v_induced, nhat_arr[k])
        u2D_mag = sqrt(u2D_t^2 + u2D_n^2)
        α = atan(u2D_n, u2D_t) * 180 / pi
        pol = polars[isurf][j]
        cl_target = FLOWMath.linear(pol.alphas, pol.cls_visc, α)
        return 0.5 * cl_target * c_arr[k] * u2D_mag
    end

    function write_station!(k, new_total)
        nc_k = station_nc[k]
        if nc_k == 1
            system.Γ[station_iΓ[k]] = new_total
        else
            first = station_iΓ_first[k]
            old_total = system.Γ[station_iΓ[k]]
            if old_total == 0
                for i in 1:nc_k
                    system.Γ[first + i - 1] = new_total * i / nc_k
                end
            else
                r = new_total / old_total
                for i in 1:nc_k
                    system.Γ[first + i - 1] *= r
                end
            end
        end
        return nothing
    end

    x = [system.Γ[station_iΓ[k]] for k in 1:n]
    Xs = Vector{Vector{TF}}(); Fs = Vector{Vector{TF}}()
    refresh_forces!()
    for it in 0:maxiter-1
        f = zeros(TF, n)
        for k in 1:n
            f[k] = station_target(k) - x[k]
        end
        resid = maximum(abs.(f)); Γscale = maximum(abs.(x)) + eps(TF)
        if VISC_ITER_DEBUG[]
            println("  [visc_iter_shed] it=$it resid/scale=$(resid/Γscale)")
        end
        resid / Γscale < tol && break
        push!(Xs, copy(x)); push!(Fs, copy(f))
        mk = min(m, length(Xs) - 1)
        if mk == 0
            x = x .+ beta .* f
        else
            ΔF = hcat([Fs[end-i+1] .- Fs[end-i] for i in 1:mk]...)
            ΔX = hcat([Xs[end-i+1] .- Xs[end-i] for i in 1:mk]...)
            γ = ΔF \ f
            x = x .+ beta .* f .- (ΔX .+ beta .* ΔF) * γ
        end
        for k in 1:n
            write_station!(k, x[k])
        end
        refresh_forces!()
    end
    for k in 1:n
        write_station!(k, x[k])
    end
    refresh_forces!()

    # ---- viscous drag: missing block, found while chasing the residual CD gap after the
    # row-1-staleness fix (vortexlattice_row1_staleness_bug_ROOT_CAUSE_FIXED_20260820) --
    # this function only ever converged Γ toward the polar's cl_target (lift), it never added
    # the polar's cd_target profile-drag contribution to cfb the way viscous_iterative! (steady)
    # and viscous!() (single-pass unsteady) both do. Same convention as viscous_iterative!'s
    # drag block below (cd*q_local*c*Δs_y added to cfb, all in the global frame), distributed
    # across every chordwise row when nc>1 (see viscous_iterative!'s equivalent block).
    q = 0.5 * RHO * ref.V * ref.V
    for k in 1:n
        isurf = station_isurf[k]; j = station_j[k]; nc_k = station_nc[k]
        surface = system.surfaces[isurf]
        v_induced = system.properties[isurf][1, j].velocity * ref.V
        u2D_t = dot(v_induced, that_arr[k]); u2D_n = dot(v_induced, nhat_arr[k])
        u2D_mag = sqrt(u2D_t^2 + u2D_n^2)
        α = atan(u2D_n, u2D_t) * 180 / pi
        pol = polars[isurf][j]
        cd = FLOWMath.linear(pol.alphas, pol.cds_visc, α)
        q_local = 0.5 * RHO * u2D_mag^2
        c = c_arr[k]
        Δs_y = norm(surface[1,j].rtl - surface[1,j].rtr)
        dhat_global = (u2D_t * that_arr[k] + u2D_n * nhat_arr[k]) / u2D_mag
        D_visc_strip = cd * q_local * c * Δs_y
        d_viscous = (D_visc_strip / nc_k) * dhat_global / (q * ref.S)
        for i in 1:nc_k
            (; gamma, velocity, cfb, cfl, cfr) = system.properties[isurf][i, j]
            system.properties[isurf][i, j] = PanelProperties(gamma, velocity, cfb + d_viscous, cfl, cfr)
        end
    end

    Γ_wake .= system.Γ
    # DO NOT restore system.Γ to Γ_seed here (previously done, see
    # vortexlattice_row1_staleness_bug_ROOT_CAUSE_FIXED_20260820): the next step's
    # `update_trailing_edge_filaments!(trailing_edge_filaments, system.surfaces, system.Γ)`
    # (unsteady.jl, called BEFORE that step's own AIC solve overwrites system.Γ) reads
    # system.Γ to refresh the near-wake ring closest to the TE -- restoring to the inviscid
    # value here silently discarded the viscous correction on that ring every single step,
    # right before wake_on_all! used it. Leaving system.Γ at the converged viscous value keeps
    # that ring consistent with what shed_wake! (below, via Γ_wake) actually just shed; the
    # AIC solve next step still overwrites system.Γ completely regardless (ldiv! has no
    # dependence on Γ's incoming value), so nothing else observes this changed convention.

    return nothing
end
