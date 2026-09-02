# --- Diagnostic instrumentation: capture per-section effective AoA (deg) computed by the ---
# --- viscous correction below, keyed by surface index, one Vector per call to viscous!.  ---
# --- Added 2026-07-29 to compare local AoA against a CCBlade reference; does not affect   ---
# --- forces/circulation. Enable with `ALPHA_LOG_ENABLED[] = true`; read out `ALPHA_LOG`.  ---
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
function reset_alpha_log!()
    empty!(ALPHA_LOG)
    empty!(VELOCITY_LOG)
    empty!(CL_LOG)
    empty!(CD_LOG)
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
    _dyncamber_calls[] += 1

    # dynamic pressure for dimensionalizing forces
    q = 0.5 * RHO * ref.V * ref.V

    # diagnostic AoA capture for this call (see ALPHA_LOG above)
    _alpha_this_call = ALPHA_LOG_ENABLED[] ? Dict{Int, Vector{Float64}}() : nothing
    _velocity_this_call = ALPHA_LOG_ENABLED[] ? Dict{Int, Vector{Tuple{Float64,Float64}}}() : nothing
    _cl_this_call = ALPHA_LOG_ENABLED[] ? Dict{Int, Vector{Float64}}() : nothing
    _cd_this_call = ALPHA_LOG_ENABLED[] ? Dict{Int, Vector{Float64}}() : nothing

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

                # calculate effective cl predicted by the VLM, = 2π * α_eff
                cl_vlm = -2 * RHO * γ * (γ / (l_2d_norm + eps(l_2d_norm))) / c * sign(γ)

                # get effective α
                # dynamic camber (2026-09-02): the lattice already carries 2π·Δα of the polar's extra lift
                dα_dc = (DYNAMIC_CAMBER[] && !isempty(CAMBER_DALPHA)) ? CAMBER_DALPHA[isurf][j] : 0.0 # deg
                cl_vlm_true = cl_vlm - 2 * pi * dα_dc * pi / 180
                α_eff = cl_vlm_true / (2 * pi) * 180 / pi # in degrees
                if DYNAMIC_CAMBER[] && DYNAMIC_CAMBER_ALPHA_V[] && !isempty(CAMBER_ALPHA0)
                    # mode B: geometric flow angle at the ¼-chord from the strip velocity, measured from
                    # the zero-lift line (α_L0 from CAMBER_ALPHA0); independent of the Δα rotation
                    pn1 = surface[1, j]; pnN = surface[end, j]
                    cvec = 0.5 * (pnN.rbl + pnN.rbr) - 0.5 * (pn1.rtl + pn1.rtr); cvec /= norm(cvec)
                    svec = pn1.rtr - pn1.rtl
                    ngeo = cross(cvec, svec); ngeo /= norm(ngeo)
                    dot(ngeo, pn1.ncp) < 0 && (ngeo = -ngeo)
                    vsec = v_induced - (dot(v_induced, svec) / dot(svec, svec)) * svec
                    α_geo = -atan(dot(vsec, ngeo), dot(vsec, cvec)) * 180 / pi   # ncp points to the pressure side
                    α_eff = α_geo - CAMBER_ALPHA0[isurf][j]
                end

                if !isnothing(_alpha_this_call)
                    push!(_alpha_this_call[isurf], α_eff)
                end

                # additive viscous correction: Δcl(cl_inv), indexed by the
                # inviscid cl so the table stays smooth near zero lift
                Δcl = FLOWMath.linear(polar.cls_inv, polar.cls_delta, cl_vlm_true)
                if DYNAMIC_CAMBER[] && !isempty(CAMBER_DALPHA)
                    ω_dc = DYNAMIC_CAMBER_OMEGA[]
                    if DYNAMIC_CAMBER_ALPHA_V[]
                        # mode C: feed-forward from the polar at the velocity-based α (no lattice lift in the
                        # loop); Δα → (cl_visc(α) − 2πα)/2π, relaxed, held at 0 during the start transient.
                        # The force adds cl_visc(α) − cl_vlm so the total is polar-consistent either way.
                        cl_target = FLOWMath.linear(polar.alphas, polar.cls_visc, α_eff)
                        Δcl = cl_target - cl_vlm
                        if _dyncamber_calls[] >= DYNAMIC_CAMBER_START[]
                            dα_ff = (cl_target - 2 * pi * α_eff * pi / 180) / (2 * pi) * 180 / pi
                            CAMBER_DALPHA[isurf][j] = (1 - ω_dc) * dα_dc + ω_dc * dα_ff
                        end
                    else
                        CAMBER_DALPHA[isurf][j] = (1 - ω_dc) * dα_dc + ω_dc * Δcl / (2 * pi) * 180 / pi
                        Δcl -= 2 * pi * dα_dc * pi / 180   # only the part the lattice does not yet produce
                    end
                    DYNAMIC_CAMBER_LOG[] && isurf == 1 && (j in (8, 13, 20)) &&
                        println("DYNCAMBER j=$j α_eff=$(round(α_eff, digits=2)) cl_vlm=$(round(cl_vlm, digits=3)) Δcl=$(round(Δcl, digits=3)) Δα=$(round(CAMBER_DALPHA[isurf][j], digits=3)) Γ=$(round(γ, digits=2))")
                end

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
                cd_raw = FLOWMath.linear(polar.alphas, polar.cds_visc, α_eff)
                if !isnothing(_cd_this_call)
                    push!(_cd_this_call[isurf], cd_raw)
                end
                cd = cd_raw
                D_visc_strip = cd * q_local * c * Δs_y
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
    end
end

function viscous!(properties::Vector{Matrix{PanelProperties{TF}}}, Γ, dΓdt, surfaces::Vector{Matrix{SurfacePanel{TF}}}, grids, frames::Vector{<:ReferenceFrame}, frames_index::Vector{Int}, polar::Nothing, ref, dt) where TF
    return nothing
end
