"""
    steady_analysis(surfaces, reference, freestream; kwargs...)

Perform a steady vortex lattice method analysis.  Return an object of type
[`System`](@ref) containing the system state.

# Arguments
 - `surfaces`:
   - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
   or
   - Vector of matrices of shape (nc, ns) containing surface panels (see
   [`SurfacePanel`](@ref))
   where `nc` is the number of chordwise panels and `ns` is the number of
   spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))

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
 - `additional_velocity`: Function which defines additional velocity (divided by
    the freestream velocity) as a function of location.
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
function steady_analysis(surfaces, reference, freestream; kwargs...)

    # pre-allocate system storage
    system = System(surfaces)

    return steady_analysis!(system, surfaces, reference, freestream; kwargs...,
        calculate_influence_matrix = true)
end

"""
    steady_analysis!(system, surfaces, reference, freestream; kwargs...)

Pre-allocated version of `steady_analysis`.
"""
function steady_analysis!(system, surfaces, ref, fs;
    symmetric = fill(false, length(surfaces)),
    wakes = [Matrix{WakePanel{Float64}}(undef, 0, size(surfaces[i],2)) for i = 1:length(surfaces)],
    nwake = size.(wakes, 1),
    surface_id = 1:length(surfaces),
    wake_finite_core = fill(true, length(surfaces)),
    trailing_vortices = fill(true, length(surfaces)),
    xhat = SVector(1, 0, 0),
    additional_velocity = nothing,
    fcore = (c, Δs) -> 1e-3,
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # number of surfaces
    nsurf = length(surfaces)

    # check if input is a grid
    grid_input = typeof(surfaces) <: AbstractVector{<:AbstractArray{<:Any, 3}}

    # if only one value is provided, use it for all surfaces
    symmetric = isa(symmetric, Number) ? fill(symmetric, nsurf) : symmetric
    nwake = isa(nwake, Number) ? fill(nwake, nsurf) : nwake
    surface_id = isa(surface_id, Number) ? fill(surface_id, nsurf) : surface_id
    wake_finite_core = isa(wake_finite_core, Number) ? fill(wake_finite_core, nsurf) : wake_finite_core
    trailing_vortices = isa(trailing_vortices, Number) ? fill(trailing_vortices, nsurf) : trailing_vortices

    # update surface panels stored in `system`
    if grid_input
        # generate and store surface panels
        for isurf = 1:nsurf
            update_surface_panels!(system.surfaces[isurf], surfaces[isurf]; fcore)
        end
    else
        # store provided surface panels
        for isurf = 1:nsurf
            system.surfaces[isurf] .= surfaces[isurf]
        end
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

"""
    unsteady_analysis(surfaces, reference, freestream, dt; kwargs...)

Perform a unsteady vortex lattice method analysis.  Return an object of type
[`System`](@ref) containing the final system state, a matrix of surface panels
(see [`SurfacePanel`](@ref) for each surface at each time step, a matrix of surface
panel properties (see [`PanelProperties`](@ref)) for each surface at each time step,
and a matrix of wake panels (see [`WakePanel`](@ref)) for each surface at each time
step.

# Arguments
 - `surfaces`:
   - Grids of shape (3, nc+1, ns+1) which represent lifting surfaces
     or
   - Matrices of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns)
    where `nc` is the number of chordwise panels and `ns` is the number of
    spanwise panels. Alternatively, a vector containing surface shapes/positions
    at each time step (including at `t=0`) may be provided to model
    moving/deforming lifting surfaces.
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters for each time step (see [`Freestream`](@ref))
 - `dt`: Time step vector

# Keyword Arguments
 - `symmetric`: Flags indicating whether a mirror image (across the X-Z plane) should
    be used when calculating induced velocities, defaults to `false` for each surface
 - `initial_wakes`: Vector of initial wakes corresponding to each surface, represented
    by matrices of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where
    `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels. Defaults to no wake panels for each surface
 - `initial_circulation`: Vector containing the initial circulation of all surface
    panels in the system.  Defaults to `zeros(N)` where `N` is the total number
    of surface panels in `surfaces`.
 - `nwake`: Maximum number of wake panels in the chordwise direction for each
    surface.  Defaults to `length(dx)` for all surfaces.
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.  Defaults to `true` for each surface.
 - `save`: Time indices at which to save the time history, defaults to `1:length(dx)`
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix needs to be calculated.  Re-using the same AIC matrix
    will (slightly) reduce calculation times when the underlying geometry has not changed.
    Defaults to `true`.  Note that this argument only affects the pre-allocated
    version of this function.
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces for the final
    time step. Defaults to `true`.
 - `derivatives`: Flag indicating whether the derivatives with respect to the
    freestream variables should be calculated for the final time step. Defaults
    to `true`.
"""
unsteady_analysis

# same grids/surfaces at each time step
function unsteady_analysis(surfaces::AbstractVector{<:AbstractArray}, ref, fs, dt;
    nwake = fill(length(dt), length(surfaces)), kwargs...)

    # pre-allocate system storage
    system = System(surfaces; nw = nwake)

    return unsteady_analysis!(system, surfaces, ref, fs, dt;
        kwargs..., nwake, calculate_influence_matrix = true)
end

# different grids/surfaces at each time step
function unsteady_analysis(surfaces::AbstractVector{<:AbstractVector{<:AbstractArray}},
    ref, fs, dt; nwake = fill(length(dt), length(surfaces[1])), kwargs...)

    # pre-allocate system storage
    system = System(surfaces[1]; nw = nwake)

    return unsteady_analysis!(system, surfaces, ref, fs, dt;
        kwargs..., nwake, calculate_influence_matrix = true)
end

"""
    unsteady_analysis!(system, surfaces, reference, freestream, dt; kwargs...)

Pre-allocated version of `unsteady_analysis`.
"""
function unsteady_analysis!(system, surfaces, ref, fs, dt;
    symmetric = fill(false, length(system.surfaces)),
    initial_circulation = zero(system.Γ),
    initial_wakes = [Matrix{WakePanel{Float64}}(undef, 0, size(surfaces[i], 2)) for i = 1:length(system.surfaces)],
    nwake = fill(length(dt), length(system.surfaces)),
    surface_id = 1:length(system.surfaces),
    wake_finite_core = fill(true, length(system.surfaces)),
    additional_velocity = nothing,
    fcore = (c, Δs) -> 1e-3,
    save = 1:length(dt),
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # float number type
    TF = eltype(system)

    # number of surfaces
    nsurf = length(system.surfaces)

    # surface motion?
    surface_motion = typeof(surfaces) <: AbstractVector{<:AbstractVector{<:AbstractArray}}

    # grid input?
    if surface_motion
        grid_input = typeof(surfaces) <: AbstractVector{<:AbstractVector{<:AbstractArray{<:Any, 3}}}
    else
        grid_input = typeof(surfaces) <: AbstractVector{<:AbstractArray{<:Any, 3}}
    end

    # convert scalar inputs to vectors
    symmetric = isa(symmetric, Number) ? fill(symmetric, nsurf) : symmetric
    nwake = isa(nwake, Number) ? fill(nwake, nsurf) : nwake
    surface_id = isa(surface_id, Number) ? fill(surface_id, nsurf) : surface_id
    wake_finite_core = isa(wake_finite_core, Number) ? fill(wake_finite_core, nsurf) : wake_finite_core

    # convert single freestream input to vector (if applicable)
    fs = isa(fs, Freestream) ? fill(fs, length(dt)) : fs

    # update surface panels stored in `system`
    if grid_input
        # generate and store surface panels
        if surface_motion
            for isurf = 1:nsurf
                update_surface_panels!(system.surfaces[isurf], grids[1][isurf]; fcore)
            end
        else
            for isurf = 1:nsurf
                update_surface_panels!(system.surfaces[isurf], grids[isurf]; fcore)
            end
        end
    else
        # store provided surface panels
        if surface_motion
            for isurf = 1:nsurf
                system.surfaces[isurf] .= surfaces[1][isurf]
            end
        else
            for isurf = 1:nsurf
                system.surfaces[isurf] .= surfaces[isurf]
            end
        end
    end

    # check if existing wake panel storage is sufficient, replace if necessary
    for isurf = 1:nsurf
        if size(system.wakes[isurf], 1) < nwake[isurf]
            # get surface/wake dimensions
            nc, ns = size(surfaces[isurf])
            nw = nwake[isurf]

            # update wake panel storage
            system.wakes[isurf] = Matrix{WakePanel{TF}}(undef, nwake[isurf], size(surfaces[isurf], 2))
        end
    end

    # copy initial wake panels to pre-allocated storage
    for isurf = 1:nsurf
        for I in CartesianIndices(initial_wakes[isurf])
            system.wakes[isurf][I] = initial_wakes[isurf][I]
        end
    end

    # zero out surface motion if there is none
    if !surface_motion
        for i = 1:nsurf
            system.Vcp[i] .= Ref(SVector(0.0, 0.0, 0.0))
            system.Vh[i] .= Ref(SVector(0.0, 0.0, 0.0))
            system.Vv[i] .= Ref(SVector(0.0, 0.0, 0.0))
            system.Vte[i] .= Ref(SVector(0.0, 0.0, 0.0))
        end
    end

    # update other parameters stored in `system`
    system.reference[] = ref
    system.freestream[] = fs[1]
    system.symmetric .= symmetric
    system.nwake .= nwake
    system.surface_id .= surface_id
    system.wake_finite_core .= wake_finite_core
    system.trailing_vortices .= false
    system.Γ .= initial_circulation

    # unpack variables stored in `system`
    AIC = system.AIC
    w = system.w
    Γ = system.Γ
    current_surfaces = system.surfaces
    previous_surfaces = system.previous_surfaces
    properties = system.properties
    dproperties = system.dproperties
    wakes = system.wakes
    wake_velocities = system.V
    dw = system.dw
    dΓ = system.dΓ
    wake_shedding_locations = system.wake_shedding_locations
    Vcp = system.Vcp
    Vh = system.Vh
    Vv = system.Vv
    Vte = system.Vte
    dΓdt = system.dΓdt

    # find intersecting surfaces
    repeated_points = repeated_trailing_edge_points(current_surfaces)

    # set the initial number of wake panels for each surface
    iwake = [min(size(initial_wakes[isurf], 1), nwake[isurf]) for isurf = 1:nsurf]

    # trailing vortices are disabled
    trailing_vortices = fill(false, nsurf)
    xhat = SVector(1, 0, 0)

    # distance to wake shedding location (most users probably shouldn't be changing this)
    eta = 0.25

    # initialize solution history for each time step
    surface_history = Vector{Vector{Matrix{SurfacePanel{TF}}}}(undef, length(save))
    property_history = Vector{Vector{Matrix{PanelProperties{TF}}}}(undef, length(save))
    wake_history = Vector{Vector{Matrix{WakePanel{TF}}}}(undef, length(save))
    isave = 1

    # # loop through all time steps
    for it = 1 : length(dt)

        # NOTE: Each step models the transition from `t = t[it]` to `t = [it+1]`
        # (e.g. the first step models from `t = 0` to `t = dt`).  Properties are
        # assumed to be constant during each time step and resulting transient
        # forces/moments are provided at `t[it+1]` (e.g. the properties returned
        # for the first step correspond to `t = dt`).

        first_step = it == 1
        last_step = it == length(dt)

        if surface_motion
            # move geometry and calculate velocities for this time step

            for isurf = 1:nsurf

                # current surfaces are now previous surfaces
                previous_surfaces[isurf] .= current_surfaces[isurf]

                # set new surface shape...
                if grid_input
                    # ...based on grid inputs
                    update_surface_panels!(current_surfaces[isurf], surfaces[1+it][isurf]; fcore)
                else
                    # ...based on surface panels
                    current_surfaces[isurf] .= surfaces[1+it][isurf]
                end

                # calculate surface velocities
                get_surface_velocities!(Vcp[isurf], Vh[isurf], Vv[isurf], Vte[isurf],
                    current_surfaces[isurf], previous_surfaces[isurf], dt[it])
            end

        end

        # update stored freestream parameters for this time step
        system.freestream[] = fs[it]

        # update number of wake panels for this time step
        for isurf = 1:nsurf
            system.nwake[isurf] = iwake[isurf]
        end

        # update the wake shedding location for this time step
        update_wake_shedding_locations!(wakes, wake_shedding_locations,
            current_surfaces, ref, fs[it], dt[it], additional_velocity, Vte,
            iwake, eta)

        # calculate/re-calculate AIC matrix (if necessary)
        if first_step && calculate_influence_matrix
            influence_coefficients!(AIC, current_surfaces;
                symmetric = symmetric,
                wake_shedding_locations = wake_shedding_locations,
                surface_id = surface_id,
                trailing_vortices = trailing_vortices,
                xhat = xhat)
        end

        # update the AIC matrix to use the new wake shedding locations
        update_trailing_edge_coefficients!(AIC, current_surfaces;
            symmetric = symmetric,
            wake_shedding_locations = wake_shedding_locations,
            trailing_vortices = trailing_vortices)

        # calculate RHS
        if last_step && derivatives
            normal_velocity_derivatives!(w, dw, current_surfaces, wakes,
                ref, fs[it]; additional_velocity, Vcp, symmetric, nwake=iwake,
                surface_id, wake_finite_core, trailing_vortices, xhat)
        else
            normal_velocity!(w, current_surfaces, wakes, ref, fs[it];
                additional_velocity, Vcp, symmetric, nwake=iwake, surface_id,
                wake_finite_core, trailing_vortices, xhat)
        end

        # save (negative) previous circulation in dΓdt
        dΓdt .= .-Γ

        # solve for the new circulation
        if last_step && derivatives
            circulation_derivatives!(Γ, dΓ, AIC, w, dw)
        else
            circulation!(Γ, AIC, w)
        end

        # solve for dΓdt using finite difference `dΓdt = (Γ - Γp)/dt`
        dΓdt .+= Γ # add newly computed circulation
        dΓdt ./= dt[it] # divide by corresponding time step

        # compute transient forces on each panel (if necessary)
        if it in save || (last_step && near_field_analysis)
            if last_step && derivatives
                near_field_forces_derivatives!(properties, dproperties,
                    current_surfaces, wakes, ref, fs[it], Γ, dΓ; dΓdt,
                    additional_velocity, Vh, Vv, symmetric, nwake=iwake,
                    surface_id, wake_finite_core, wake_shedding_locations,
                    trailing_vortices, xhat)
            else
                near_field_forces!(properties, current_surfaces, wakes,
                    ref, fs[it], Γ; dΓdt, additional_velocity, Vh, Vv,
                    symmetric, nwake=iwake, surface_id, wake_finite_core,
                    wake_shedding_locations, trailing_vortices, xhat)
            end
        end

        # save the panel history
        if it in save
            surface_history[isave] = [copy(system.surfaces[isurf]) for isurf = 1:nsurf]
            property_history[isave] = [copy(system.properties[isurf]) for isurf = 1:nsurf]
            wake_history[isave] = [wakes[isurf][1:iwake[isurf], :] for isurf = 1:nsurf]
            isave += 1
        end

        # update wake velocities
        get_wake_velocities!(wake_velocities, current_surfaces,
            wakes, ref, fs[it], Γ, additional_velocity, Vte, symmetric,
            repeated_points, iwake, surface_id, wake_finite_core,
            wake_shedding_locations, trailing_vortices, xhat)

        # shed additional wake panel (and translate existing wake panels)
        shed_wake!(wakes, wake_shedding_locations, wake_velocities,
            dt[it], current_surfaces, Γ, iwake)

        # increment wake panel counter for each surface
        for isurf = 1:nsurf
            if iwake[isurf] < nwake[isurf]
                iwake[isurf] += 1
            end
        end

    end

    # save flags indicating whether certain analyses have been performed
    system.near_field_analysis[] = near_field_analysis
    system.derivatives[] = derivatives

    # return the modified system and time history
    return system, surface_history, property_history, wake_history
end
