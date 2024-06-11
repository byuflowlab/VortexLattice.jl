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
    derivatives = true,
    vertical_segments = false)

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
    filaments = system.filaments
    AIC = system.AIC
    w = system.w
    Γ = system.Γ
    dw = system.dw
    dΓ = system.dΓ
    properties = system.properties
    filament_properties = system.filament_properties
    dproperties = system.dproperties
    fmm_velocity_probes = system.fmm_velocity_probes

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
                xhat = xhat, vertical_segments = vertical_segments)
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
                xhat = xhat, vertical_segments = vertical_segments)

            # # testing filaments here

            # # update filaments with current surface locations and strengths
            # wake_shedding_locations = nothing
            # dΓdt = nothing
            # update_filaments!(filaments, surfaces, wake_shedding_locations, Γ, dΓdt; trailing_edge=true, trailing_vortices)

            # # update probes with current surface locations
            # update_probes!(fmm_velocity_probes, surfaces, 0)

            # # compute the surface influence on the probes
            # FastMultipole.direct!(fmm_velocity_probes, system; scalar_potential=false, vector_potential=false, velocity=true, velocity_gradient=false)

            # # compute the semi-infinite wake influence on the probes
            # trailing_induced_velocity(fmm_velocity_probes, surfaces, Γ, trailing_vortices; xhat, symmetric)

            # # compute forces
            # near_field_forces!(filament_properties, filaments, surfaces, ref, fs, fmm_velocity_probes;
            #     additional_velocity, Vh=nothing, Vv=nothing, vertical_segments)

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
    surface.  Defaults to `length(dt)` for all surfaces.
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
    nwake = fill(length(dt), length(surfaces)),
    fmm_toggle=false, fmm_p=4, fmm_ncrit=20, fmm_theta=0.4,
    kwargs...)

    # pre-allocate system storage
    system = System(surfaces; nw = nwake, fmm_toggle, fmm_p, fmm_ncrit, fmm_theta)

    return unsteady_analysis!(system, surfaces, ref, fs, dt;
        kwargs..., nwake, calculate_influence_matrix = true)
end

# different grids/surfaces at each time step
function unsteady_analysis(surfaces::AbstractVector{<:AbstractVector{<:AbstractArray}},
    ref, fs, dt; nwake = fill(length(dt), length(surfaces[1])),
    fmm_toggle=false, fmm_p=4, fmm_ncrit=20, fmm_theta=0.4,
    kwargs...)

    # pre-allocate system storage
    system = System(surfaces[1]; nw = nwake, fmm_toggle, fmm_p, fmm_ncrit, fmm_theta)

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
    derivatives = true,
    vertical_segments = false)

    # float number type
    TF = eltype(system)

    # number of surfaces
    nsurf = length(system.surfaces)

    # surface motion?
    surface_motion = eltype(surfaces) <: AbstractVector

    # --- Input Pre-Processing --- #

    # convert scalar inputs to appropriately sized vectors
    symmetric = isa(symmetric, Number) ? fill(symmetric, nsurf) : symmetric
    nwake = isa(nwake, Number) ? fill(nwake, nsurf) : nwake
    surface_id = isa(surface_id, Number) ? fill(surface_id, nsurf) : surface_id
    wake_finite_core = isa(wake_finite_core, Number) ? fill(wake_finite_core, nsurf) : wake_finite_core
    fs = isa(fs, Freestream) ? fill(fs, length(dt)) : fs

    # extract initial surface panels from input
    if surface_motion
        # surface moves, store initial surface panel locations
        initial_surfaces = surfaces[1]
    else
        # surface doesn't move
        initial_surfaces = surfaces
    end

    # --- Update System Parameters --- #

    system.reference[] = ref
    system.symmetric .= symmetric
    system.surface_id .= surface_id
    system.wake_finite_core .= wake_finite_core
    system.trailing_vortices .= false

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

    # find intersecting surfaces
    repeated_points = repeated_trailing_edge_points(initial_surfaces)

    # --- Set Initial Simulation Variables --- #

    # store initial surface panels in `system`
    if eltype(initial_surfaces) <: AbstractArray{<:Any, 3}
        # initial surfaces are input as a grid, convert to surface panels
        for isurf = 1:nsurf
            update_surface_panels!(system.surfaces[isurf], initial_surfaces[isurf]; fcore)
        end
    else
        # initial surfaces are input as matrices of surface panels
        for isurf = 1:nsurf
            system.surfaces[isurf] .= initial_surfaces[isurf]
        end
    end

    # store initial wake panels in `system`
    for isurf = 1:nsurf
        for I in CartesianIndices(initial_wakes[isurf])
            system.wakes[isurf][I] = initial_wakes[isurf][I]
        end
    end

    # store initial freestream parameters in `system`
    system.freestream[] = fs[1]

    # store initial circulation parameters in `system`
    system.Γ .= initial_circulation

    # set the initial number of wake panels for each surface
    iwake = [min(size(initial_wakes[isurf], 1), nwake[isurf]) for isurf = 1:nsurf]

    # --- Begin Simulation --- #

    # initialize solution history for each time step
    surface_history = Vector{Vector{Matrix{SurfacePanel{TF}}}}(undef, length(save))
    property_history = Vector{Vector{Matrix{PanelProperties{TF}}}}(undef, length(save))
    wake_history = Vector{Vector{Matrix{WakePanel{TF}}}}(undef, length(save))
    isave = 1

    # # loop through all time steps
    for it = 1 : length(dt)

        first_step = it == 1
        last_step = it == length(dt)

        if surface_motion
            propagate_system!(system, surfaces[1+it], fs[it], dt[it];
                additional_velocity, reuse_kinematic_velocity=false,
                repeated_points,
                nwake = iwake,
                eta = 0.25,
                calculate_influence_matrix = first_step && calculate_influence_matrix,
                near_field_analysis = it in save || (last_step && near_field_analysis),
                derivatives = last_step && derivatives, vertical_segments = vertical_segments,
                fmm_velocity_probes = system.fmm_toggle ? system.fmm_velocity_probes : nothing)
        else
            propagate_system!(system, fs[it], dt[it];
                additional_velocity, reuse_kinematic_velocity=false,
                repeated_points,
                nwake = iwake,
                eta = 0.25,
                calculate_influence_matrix = first_step && calculate_influence_matrix,
                near_field_analysis = it in save || (last_step && near_field_analysis),
                derivatives = last_step && derivatives, vertical_segments = vertical_segments,
                fmm_velocity_probes = system.fmm_toggle ? system.fmm_velocity_probes : nothing)
        end

        # increment wake panel counter for each surface
        for isurf = 1:nsurf
            if iwake[isurf] < nwake[isurf]
                iwake[isurf] += 1
            end
        end

        # save the surface shape, properties, and resulting shed wake
        if it in save
            surface_history[isave] = [copy(system.surfaces[isurf]) for isurf = 1:nsurf]
            property_history[isave] = [copy(system.properties[isurf]) for isurf = 1:nsurf]
            wake_history[isave] = [system.wakes[isurf][1:iwake[isurf], :] for isurf = 1:nsurf]
            isave += 1
        end

    end

    # return the modified system and time history
    return system, surface_history, property_history, wake_history
end

"""
    propagate_system!(system, [surfaces, ] freestream, dt; kwargs...)

Propagate the state variables in `system` forward one time step using the
unsteady vortex lattice method system of equations.

# Arguments
 - `system`: Object of type `system` which contains the current system state
 - `surfaces`: Surface locations at the end of this time step. If omitted,
   surfaces are assumed to be stationary.
 - `freestream`: Freestream parameters corresponding to this time step.
 - `dt`: Time increment

# Keyword Arguments
 - `additional_velocity`: Function which defines additional velocity as a
    function of location.
 - `repeated_points`: Dictionary of the form `Dict((isurf, i) => [(jsurf1, j1),
    (jsurf2, j2)...]` which defines repeated trailing edge points.  Trailing edge
    point `i` on surface `isurf` is repeated on surface `jsurf1` at point `j1`,
    `jsurf2` at point `j2`, and so forth. See [`repeated_trailing_edge_points`](@ref)
 - `nwake`: Number of wake panels in the chordwise direction for each surface.
 - `eta`: Time step fraction used to define separation between trailing
    edge and wake shedding location.  Typical values range from 0.2-0.3.
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix needs to be calculated. If argument `surfaces` is provided
    the influence matrix will always be recalculated.
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces.
 - `derivatives`: Flag indicating whether the derivatives with respect to the
    freestream variables should be calculated.
"""
propagate_system!

# stationary surfaces
propagate_system!(system, fs, dt; kwargs...) = propagate_system!(system,
    nothing, fs, dt; kwargs...)

const DEBUG = Array{Bool,0}(undef)
DEBUG[] = false

function wake_on_all!(system::System)
    # compile sources
    update_filaments!(system.filaments, system.wakes, system.nwake)
    if DEBUG[] && isnan(system)
        throw("found nans: 6")
    end

    # preallocate probes
    n_control_points = 0
    n_surface_filaments = 0
    for surface in system.surfaces
        nc, ns = size(surface)
        n_control_points += nc * ns # control points
        n_surface_filaments += nc * ns + nc * ns + nc # top (nc*ns) + left (nc*ns) + right (nc)
    end

    # n_wake_shedding_locations = 0
    # for wake_shedding_location in wake_shedding_locations
    #     n_wake_shedding_locations += length(wake_shedding_location)
    # end

    n_wake_corners = 0
    for (nc, wake) in zip(system.nwake, system.wakes)
        if nc > 0
            ns = size(wake, 2)
            n_wake_corners += (nc+1)*(ns+1)
        end
    end

    update_n_probes!(system.fmm_velocity_probes, n_control_points + n_surface_filaments + n_wake_corners)

    # compile targets
    update_probes!(system.fmm_velocity_probes, system.surfaces, 0) # control points and filament centers
    update_probes!(system.fmm_velocity_probes, system.wakes, system.nwake, n_control_points + n_surface_filaments) # wake corners

    # reset to zero
    reset!(system.fmm_velocity_probes)
    if DEBUG[] && isnan(system.fmm_velocity_probes)
        throw("found nans! 6.5")
    end

    # run FMM
    if length(system.filaments) > 0
        if system.fmm_toggle[]
            FastMultipole.fmm!(system.fmm_velocity_probes, system; expansion_order=system.fmm_p, leaf_size_source=system.fmm_ncrit, leaf_size_target=system.fmm_ncrit, multipole_threshold=system.fmm_theta)
        else
            FastMultipole.direct!(system.fmm_velocity_probes, system)
        end
        if DEBUG[] && isnan(system)
            throw("found nans: 7")
        end
    end
end

# moving/deforming surfaces
function propagate_system!(system, surfaces, fs, dt;
    additional_velocity, reuse_kinematic_velocity,
    repeated_points,
    nwake,
    eta,
    calculate_influence_matrix,
    near_field_analysis,
    derivatives, vertical_segments,
    fmm_velocity_probes)

    # NOTE: Each step models the transition from `t = t[it]` to `t = [it+1]`
    # (e.g. the first step models from `t = 0` to `t = dt`).  Properties are
    # assumed to be constant during each time step and resulting transient
    # forces/moments are provided at `t[it+1]` (e.g. the properties returned
    # for the first step correspond to `t = dt`).

    nsurf = length(system.surfaces)

    # unpack constant system parameters
    ref = system.reference[]
    symmetric = system.symmetric
    surface_id = system.surface_id
    wake_finite_core = system.wake_finite_core
    trailing_vortices = system.trailing_vortices
    xhat = system.xhat[]

    # unpack system storage (including state variables)
    previous_surfaces = system.previous_surfaces
    current_surfaces = system.surfaces
    properties = system.properties
    dproperties = system.dproperties
    wakes = system.wakes
    wake_velocities = system.V
    wake_shedding_locations = system.wake_shedding_locations
    AIC = system.AIC
    w = system.w
    dw = system.dw
    Γ = system.Γ
    dΓ = system.dΓ
    dΓdt = system.dΓdt
    Vcp = system.Vcp
    Vh = system.Vh
    Vv = system.Vv
    Vte = system.Vte
    fmm_toggle = system.fmm_toggle
    filaments = system.filaments
    fmm_p = system.fmm_p
    fmm_ncrit = system.fmm_ncrit
    fmm_theta = system.fmm_theta

    # check if the surfaces are moving
    surface_motion = !isnothing(surfaces)

    if DEBUG[] && isnan(system)
        throw("found nans: 1")
    end
    # move geometry and calculate velocities for this time step
    if !reuse_kinematic_velocity
        if surface_motion

            # check if grids are used to represent the new surfaces
            grid_input = eltype(surfaces) <: AbstractArray{<:Any, 3}

            for isurf = 1:nsurf

                # save current surfaces as previous surfaces
                previous_surfaces[isurf] .= current_surfaces[isurf]

                # set new surface shape...
                if grid_input
                    # ...based on grid inputs
                    update_surface_panels!(current_surfaces[isurf], surfaces[isurf]; fcore)
                else
                    # ...based on surface panels
                    current_surfaces[isurf] .= surfaces[isurf]
                end

                # calculate surface velocities
                get_surface_velocities!(Vcp[isurf], Vh[isurf], Vv[isurf], Vte[isurf],
                    current_surfaces[isurf], previous_surfaces[isurf], dt)
            end
        else
            # zero out surface motion
            for isurf = 1:nsurf
                system.Vcp[isurf] .= Ref(SVector(0.0, 0.0, 0.0))
                system.Vh[isurf] .= Ref(SVector(0.0, 0.0, 0.0))
                system.Vv[isurf] .= Ref(SVector(0.0, 0.0, 0.0))
                system.Vte[isurf] .= Ref(SVector(0.0, 0.0, 0.0))
            end
        end
    end

    if DEBUG[] && isnan(system)
        throw("found nans: 2")
    end

    # update stored freestream parameters for this time step
    system.freestream[] = fs

    # update number of wake panels for each surface for this time step
    system.nwake .= nwake

    # update the wake shedding location for this time step
    # (based on freestream/kinematic/other velocity only)
    # become the top corners of the next generation of wake
    update_wake_shedding_locations!(wakes, wake_shedding_locations,
        current_surfaces, ref, fs, dt, additional_velocity, Vte,
        nwake, eta)

    if DEBUG[] && isnan(system)
        throw("found nans: 3")
    end
    # calculate/re-calculate AIC matrix (if necessary)
    if surface_motion || calculate_influence_matrix
        influence_coefficients!(AIC, current_surfaces;
            symmetric = symmetric,
            wake_shedding_locations = wake_shedding_locations,
            surface_id = surface_id,
            trailing_vortices = trailing_vortices,
            xhat = xhat)
    end

    if DEBUG[] && isnan(system)
        throw("found nans: 4")
    end
    # update the AIC matrix to use the new wake shedding locations
    update_trailing_edge_coefficients!(AIC, current_surfaces;
        symmetric = symmetric,
        wake_shedding_locations = wake_shedding_locations,
        trailing_vortices = trailing_vortices)

    if DEBUG[] && isnan(system)
        throw("found nans: 5")
    end
    # wake-on-all
    if fmm_toggle # wake on all influence
        wake_on_all!(system)
    end

    # calculate RHS (doesn't store velocity induced by wake)
    if derivatives
        normal_velocity_derivatives!(w, dw, current_surfaces, wakes,
        ref, fs, fmm_velocity_probes; additional_velocity, Vcp, symmetric, nwake,
        surface_id, wake_finite_core, trailing_vortices, xhat)
    else
        normal_velocity!(w, current_surfaces, wakes, ref, fs, fmm_velocity_probes;
            additional_velocity, Vcp, symmetric, nwake, surface_id,
            wake_finite_core, trailing_vortices, xhat)
    end
    if DEBUG[] && isnan(system)
        throw("found nans: 8")
    end

    # save (negative) previous circulation in dΓdt
    dΓdt .= .-Γ

    # solve for the new circulation
    if derivatives
        circulation_derivatives!(Γ, dΓ, AIC, w, dw)
    else
        circulation!(Γ, AIC, w)
    end
    if DEBUG[] && isnan(system)
        throw("found nans: 9")
    end

    # solve for dΓdt using finite difference `dΓdt = (Γ - Γp)/dt`
    dΓdt .+= Γ # add newly computed circulation
    dΓdt ./= dt # divide by corresponding time step

    # surface-on-all
    if fmm_toggle

        # update sources
        update_filaments!(filaments, current_surfaces, wake_shedding_locations, Γ)
        if DEBUG[] && isnan(system)
            throw("found nans: 10")
        end

        # run FMM
		FastMultipole.fmm!(fmm_velocity_probes, system; expansion_order=fmm_p, leaf_size_source=fmm_ncrit, leaf_size_target=fmm_ncrit, multipole_threshold=fmm_theta)
        if DEBUG[] && isnan(system)
            throw("found nans: 11")
        end
    end

    # compute transient forces on each panel (if necessary)
    if near_field_analysis
        if derivatives
            near_field_forces_derivatives!(properties, dproperties,
                current_surfaces, wakes, ref, fs, Γ, dΓ, fmm_velocity_probes; dΓdt,
                additional_velocity, Vh, Vv, symmetric, nwake,
                surface_id, wake_finite_core, wake_shedding_locations,
                trailing_vortices, xhat, vertical_segments)
        else
            near_field_forces!(properties, current_surfaces, wakes,
                ref, fs, Γ, fmm_velocity_probes; dΓdt, additional_velocity, Vh, Vv,
                symmetric, nwake, surface_id, wake_finite_core,
                wake_shedding_locations, trailing_vortices, xhat, vertical_segments)
        end
        if DEBUG[] && isnan(system)
            throw("found nans: 12")
        end

        # save flag indicating that a near-field analysis has been performed
        system.near_field_analysis[] = near_field_analysis
    end

    # save flag indicating that derivatives wrt freestream variables have been obtained
    system.derivatives[] = derivatives

    # update wake velocities
    get_wake_velocities!(wake_velocities, current_surfaces,
        wakes, ref, fs, Γ, additional_velocity, Vte, symmetric,
        repeated_points, nwake, surface_id, wake_finite_core,
        wake_shedding_locations, trailing_vortices, xhat,
        fmm_velocity_probes)
    if DEBUG[] && isnan(system)
        throw("found nans: 13")
    end

    # shed additional wake panel (and translate existing wake panels)
    shed_wake!(wakes, wake_shedding_locations, wake_velocities,
        dt, current_surfaces, Γ, nwake)
    if DEBUG[] && isnan(system)
        throw("found nans: 14")
    end

    return system
end

function Base.isnan(array::AbstractArray{<:Real})
    for v in array
        if isnan(v)
            return true
        end
    end
    return false
end

function Base.isnan(array::AbstractArray{<:AbstractArray})
    for a in array
        if isnan(a)
            return true
        end
    end
    return false
end

function Base.isnan(object, fields::NTuple{<:Any,Symbol})
    for field in fields
        if isnan(getfield(object,field))
            throw("nan in $field")
            return true
        end
    end
    return false
end

function Base.isnan(system::System)
    for field in (:AIC,:w,:Γ,:V,:xhat,:Vcp,:Vh,:Vv,:Vte,:dΓdt)
        if isnan(getfield(system,field))
            throw("nan in system.$field")
            return true
        end
    end
    for dw in system.dw
        if isnan(dw)
            throw("nan in system.dw")
            return true
        end
    end
    for dΓ in system.dΓ
        if isnan(dΓ)
            throw("nan in system.dΓ")
            return true
        end
    end
    for s in system.surfaces
        for p in s
            if isnan(p,(:rtl, :rtc, :rtr, :rbl, :rbc, :rbr, :rcp, :ncp, :core_size, :chord))
                throw("nan in system.surfaces")
                return true
            end
        end
    end
    for prop in system.properties
        for p in prop
            if isnan(p,(:gamma,:velocity,:cfb,:cfl,:cfr))
                throw("nan in system.properties")
                return true
            end
        end
    end
    for (wake,nc) in zip(system.wakes, system.nwake)
        for j in 1:size(wake,2)
            for i in 1:nc
                p = wake[i,j]
                if isnan(p,(:rtl,:rtr,:rbl,:rbr,:core_size,:gamma))
                    throw("nan in system.wakes")
                    return true
                end
            end
        end
    end
    if isnan(system.reference[],(:S,:c,:b,:r,:V))
        return true
    end
    if isnan(system.freestream[],(:Vinf,:alpha,:beta,:Omega,:rho))
        return true
    end
    for props in system.dproperties
        for prop in props
            for p in prop
                if isnan(p,(:gamma,:velocity,:cfb,:cfl,:cfr))
                    return true
                end
            end
        end
    end
    if !isnothing(system.wake_shedding_locations)
        for wsl in system.wake_shedding_locations
            if isnan(wsl)
                @show wsl
                throw("nan in system.wake_shedding_locations")
                return true
            end
        end
    end
    for p in system.filaments
        if isnan(p,(:x1,:x2,:xc,:γ,:dγdt,:core_size))
            throw("nan in system.filaments")
            return true
        end
    end
    if isnan(system.fmm_velocity_probes)
        throw("nan in system.fmm_velocity_probes")
        return true
    end
    return false
end

