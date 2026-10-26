"""
    steady_analysis(surfaces, reference, freestream, frames; kwargs...)

Perform a steady vortex lattice method analysis.  Return an object of type
[`VortexLatticeSystem`](@ref) containing the system state.

# Arguments
 - `surfaces`:
   - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
   or
   - Vector of matrices of shape (nc, ns) containing surface panels (see
   [`SurfacePanel`](@ref))
   where `nc` is the number of chordwise panels and `ns` is the number of
   spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: the core's `Freestream` (the air velocity in the surfaces' axes)
 - `frames`: the SixDOF frame tree that moves the surfaces; every surface must share one
   rigid motion, which sets the relative wind and the body rates (see [`flight_condition`](@ref)
   to build both from angles)

# Keyword Arguments
 - `symmetric`: Flag for each surface indicating whether a mirror image across
    the X-Z plane should be used when calculating induced velocities. Defaults to
    `false` for each surface
 - `wakes`: Matrix of wake panels (see [`WakePanel`](@ref)) for each surface.  Each
    matrix has shape (nw, ns) where `nw` is the number of chordwise wake panels
    and `ns` is the number of spanwise panels for each surface, defaults to no
    wake panels for each surface
 - `nwake`: Number of chordwise wake panels to use from each wake in `wakes`,
    defaults to all wake panels for each surface
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
    By default, all surfaces are assigned their own IDs
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating a wake's influence on itself and
    surfaces/wakes with the same surface ID.  Defaults to `true` for each surface.
 - `trailing_vortices`: Flags to enable/disable trailing vortices for each surface,
    defaults to `true` for each surface
 - `xhat`: Direction in which to shed trailing vortices, defaults to `[1, 0, 0]`
 - `additional_velocity`: Function which defines additional velocity as a
    function of location.
 - `fcore`: function which sets the finite core size for each surface based on
    the chord length and/or the panel width. Defaults to `(c, Δs) -> 1e-3`.
    Only used for grid inputs.
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix has already been calculated.  Re-using the same AIC matrix
    will reduce calculation times when the underlying geometry has not changed.
    Defaults to `true`.  Note that this argument is only valid for the pre-allocated
    version of this function.
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces. Defaults to `true`.
 - `derivatives`: Flag indicating whether the derivatives with respect
    to the freestream variables should be calculated. Defaults to `true`.
"""
function steady_analysis(grids, reference, freestream::Freestream, frames; kwargs...)

    # pre-allocate system storage
    system = VortexLatticeSystem(grids)

    return steady_analysis!(system, reference, freestream, frames; kwargs...,
        calculate_influence_matrix = true)
end

"""
    steady_analysis!(system, reference, freestream, frames; kwargs...)

Pre-allocated version of `steady_analysis`.
"""
steady_analysis!(system, ref, freestream::Freestream, frames; kwargs...) =
    steady_analysis!(system, ref, body_freestream(freestream, frames, ref, length(system.surfaces)); kwargs...)

# A mirrored half (symmetric = true) is exact only in symmetric flight: sideslip or a roll
# or yaw rate is mirrored too, which silently zeroes the lateral loads (CL 3.4x off at a
# roll rate of 0.05, checked 2026-10-26 against the full mirrored wing).
function _check_symmetric_flight(fs)
    (iszero(fs.beta) && iszero(fs.Omega[1]) && iszero(fs.Omega[3])) || throw(ArgumentError(
        "symmetric = true needs symmetric flight (beta = 0, no roll or yaw rate); model the full surface (mirror = true in wing_to_grid) instead"))
    return nothing
end

function steady_analysis!(system, ref, fs::BodyFreestream;
    symmetric = fill(false, length(system.surfaces)),
    wakes = [Matrix{WakePanel{Float64}}(undef, 0, size(system.surfaces[i],2)) for i = 1:length(system.surfaces)],
    nwake = size.(wakes, 1),
    surface_id = 1:length(system.surfaces),
    wake_finite_core = fill(true, length(system.surfaces)),
    trailing_vortices = fill(true, length(system.surfaces)),
    xhat = system.xhat[],
    additional_velocity = nothing,
    fcore = (c, Δs) -> 1e-3,
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # number of surfaces
    nsurf = length(system.surfaces)

    # if only one value is provided, use it for all surfaces
    symmetric = isa(symmetric, Number) ? fill(symmetric, nsurf) : symmetric
    nwake = isa(nwake, Number) ? fill(nwake, nsurf) : nwake
    surface_id = isa(surface_id, Number) ? fill(surface_id, nsurf) : surface_id
    any(symmetric) && _check_symmetric_flight(fs)
    wake_finite_core = isa(wake_finite_core, Number) ? fill(wake_finite_core, nsurf) : wake_finite_core
    trailing_vortices = isa(trailing_vortices, Number) ? fill(trailing_vortices, nsurf) : trailing_vortices

    for isurf = 1:nsurf
        update_surface_panels!(system.surfaces[isurf], system.grids[isurf]; ratios = system.ratios[isurf], fcore)
    end

    # update other parameters stored in `system`
    system.reference[] = ref
    system.freestream[] = fs
    system.symmetric .= symmetric
    system.wakes .= wakes
    system.nwake .= nwake
    system.surface_id .= surface_id
    system.wake_finite_core .= wake_finite_core
    system.trailing_vortices .= trailing_vortices
    system.xhat[] = xhat

    # unpack variables stored in `system`
    surfaces = system.surfaces
    AIC = system.AIC
    w = system.w
    Γ = system.Γ
    dw = system.dw
    dΓ = system.dΓ
    properties = system.properties
    dproperties = system.dproperties

    # see if wake panels are being used
    wake_panels = nwake .> 0

    # calculate/re-calculate AIC matrix (if necessary)
    if calculate_influence_matrix
        influence_coefficients!(AIC, surfaces;
            symmetric = symmetric,
            surface_id = surface_id,
            trailing_vortices = trailing_vortices .& .!wake_panels,
            xhat = xhat)
    end

    # calculate RHS
    if derivatives
        normal_velocity_derivatives!(w, dw, surfaces, wakes, ref, fs;
            additional_velocity = additional_velocity,
            Vcp = nothing, # no velocity at control points due to surface motion
            symmetric = symmetric,
            surface_id = surface_id,
            nwake = nwake,
            wake_finite_core = wake_finite_core,
            trailing_vortices = trailing_vortices,
            xhat = xhat)
    else
        normal_velocity!(w, surfaces, wakes, ref, fs;
            additional_velocity = additional_velocity,
            Vcp = nothing, # no velocity at control points due to surface motion
            symmetric = symmetric,
            surface_id = surface_id,
            nwake = nwake,
            wake_finite_core = wake_finite_core,
            trailing_vortices = trailing_vortices,
            xhat = xhat)
    end

    # solve for the circulation distribution
    if derivatives
        circulation_derivatives!(Γ, dΓ, AIC, w, dw)
    else
        circulation!(Γ, AIC, w)
    end

    if near_field_analysis
        # perform a near field analysis to obtain panel properties
        if derivatives
            near_field_forces_derivatives!(properties, dproperties, surfaces, wakes,
                ref, fs, Γ, dΓ;
                dΓdt = nothing, # no unsteady forces
                additional_velocity = additional_velocity,
                Vh = nothing, # no velocity due to surface motion
                Vv = nothing, # no velocity due to surface motion
                symmetric = symmetric,
                nwake = nwake,
                surface_id = surface_id,
                wake_finite_core = wake_finite_core,
                wake_shedding_locations = nothing, # shedding location at trailing edge
                trailing_vortices = trailing_vortices,
                xhat = xhat)
        else
            near_field_forces!(properties, surfaces, wakes, ref, fs, Γ;
                dΓdt = nothing, # no unsteady forces
                additional_velocity = additional_velocity,
                Vh = nothing, # no velocity due to surface motion
                Vv = nothing, # no velocity due to surface motion
                symmetric = symmetric,
                nwake = nwake,
                surface_id = surface_id,
                wake_finite_core = wake_finite_core,
                wake_shedding_locations = nothing, # shedding location at trailing edge
                trailing_vortices = trailing_vortices,
                xhat = xhat)
        end
    end

    # save flags indicating whether certain analyses have been performed
    system.near_field_analysis[] = near_field_analysis
    system.derivatives[] = derivatives

    # return the modified system
    return system
end
