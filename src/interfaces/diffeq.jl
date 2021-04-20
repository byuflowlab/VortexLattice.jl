"""
    ODEProblem(surfaces, reference, freestream, tshed, tspan; kwargs...)

Models the unsteady vortex lattice method as an ODE problem (with a singular mass
matrix) that can be solved using the DifferentialEquations package.

# Arguments
 - `surfaces`: Surface geometry at the beginning of the simulation.
    Represented as either a:
     - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
    or
     - Vector of matrices of shape (nc, ns) containing surface panels (see [`SurfacePanel`](@ref))
       where `nc` is the number of chordwise panels and `ns` is the number of
       spanwise panels.
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters at the beginning of the simulation (see [`Freestream`](@ref))
 - `tshed`: Vector containing times at which to shed additional wake panels
 - `tspan`: Time span over which to perform the unsteady vortex lattice method analysis

# Keyword Arguments
 - `wakes`: Vector of initial wakes corresponding to each surface, represented
    by matrices of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where
    `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels. Defaults to no wake panels for each surface
 - `circulation`: Vector containing the initial circulation of all surface
    panels in the system.  Defaults to `zeros(N)` where `N` is the total number
    of surface panels in `surfaces`.
 - `surface_geometry_history`: Function which returns surface definitions
    as a function of time, where surfaces are defined as described by the argument
    `surfaces`. If not provided, no surface motion/deformation is assumed.
 - `surface_geometry_prototype`: Prototype output from `surface_geometry_history`.
    If provided, `surface_geometry_history` is assumed to be in-place.
 - `surface_velocity_history`: Function which returns surface velocities
    as a function of time, where surface velocities are defined as described by
    the argument `surfaces`, except velocities are used in place of positions.
 - `surface_velocity_prototype`: Prototype output from `surface_velocity_history`.
    If provided, `surface_velocity_history` is assumed to be in-place.
 - `freestream_history`: Function which returns freestream parameters as a function
    of time.  If not provided, the freestream parameters are assumed to be
    constant throughout the simulation.
 - `symmetric`: Flags indicating whether the geometry for each surface is
    symmetric about the X-Z plane.  Note that applying symmetry to surfaces is
    only valid when the freestream conditions are symmetric as well.
    Defaults to `false` for each surface.
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.  Defaults to `true` for each surface.
 - `nwake`: Maximum number of wake panels in the chordwise direction for each
    surface.  Defaults to `length(tshed)` for all surfaces.
 - `fcore`: function for setting the finite core size when generating surface
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width. Only used if surfaces
    are provided as grids.
 - `preserve_core_size`: Flag indicating whether the finite core size should be
    preserved from the previous set of panels when updating surface geometries.
    Defaults to `true`.
"""
function DifferentialEquations.ODEProblem(surfaces, reference::Reference, freestream::Freestream,
    tshed, tspan; nwake = fill(length(tshed), length(surfaces)), kwargs...)
    # number of steps
    nstep = length(tshed)
    # pre-allocate system storage
    system = UnsteadySystem(surfaces, reference, freestream, nstep; nwake)
    # create the ODE Problem
    return DifferentialEquations.ODEProblem(system, tshed, tspan; nwake, kwargs...)
end

"""
    ODEProblem(system::UnsteadySystem, tshed, tspan; kwargs...)

Models the unsteady vortex lattice method as an ODE problem (with a singular mass
matrix) that can be solved using the DifferentialEquations package.

# Arguments
 - `system`: Object which holds vortex lattice method inputs and current state
 - `tshed`: Vector containing times at which to shed additional wake panels
 - `tspan`: Time span over which to perform the unsteady vortex lattice method analysis

# Keyword Arguments
 - `surface_geometry_history`: Function which returns surface definitions
    as a function of time, where surfaces are defined as described by the argument
    `surfaces`. If not provided, no surface motion/deformation is assumed.
 - `surface_geometry_prototype`: Prototype output from `surface_geometry_history`.
    If provided, `surface_geometry_history` is assumed to be in-place.
 - `surface_velocity_history`: Function which returns surface velocities
    as a function of time, where surface velocities are defined as described by
    the argument `surfaces`, except velocities are used in place of positions.
 - `surface_velocity_prototype`: Prototype output from `surface_velocity_history`.
    If provided, `surface_velocity_history` is assumed to be in-place.
 - `freestream_history`: Function which returns freestream parameters as a function
    of time.  If not provided, the freestream parameters are assumed to be
    constant throughout the simulation.
 - `fcore`: function for setting the finite core size when generating surface
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width. Only used when
    surfaces are provided as grids.
 - `preserve_core_size`: Flag indicating whether the finite core size should be
    preserved from the previous set of panels when updating surface geometries.
    Defaults to `true`.

The following keyword arguments may be provided to update the inputs and states
of the vortex lattice system prior to performing the unsteady analysis.

# Additional Keyword Arguments
 - `surfaces`: Surface geometry at the beginning of the simulation.
    Represented as either a:
     - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
    or
     - Vector of matrices of shape (nc, ns) containing surface panels (see [`SurfacePanel`](@ref))
       where `nc` is the number of chordwise panels and `ns` is the number of
       spanwise panels.
 - `wakes`: Vector of initial wakes corresponding to each surface, represented
    by matrices of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where
    `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels.
 - `circulation`: Vector containing the initial circulation of all surface
    panels in the system.
 - `circulation_rates`: Vector containing the initial circulation rates of all
    surface panels in the system.
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Initial freestream parameters (see [`Freestream`](@ref))
 - `symmetric`: Flags indicating whether the geometry for each surface is
    symmetric about the X-Z plane.  Note that applying symmetry to surfaces is
    only valid when the freestream conditions are symmetric as well.
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.
 - `nwake`: Maximum number of wake panels in the chordwise direction for each
    surface.
"""
function DifferentialEquations.ODEProblem(system::UnsteadySystem{TF}, tshed, tspan;
    surface_geometry_history = nothing,
    surface_geometry_prototype = nothing,
    surface_velocity_history = nothing,
    surface_velocity_prototype = nothing,
    freestream_history = nothing,
    calculate_influence_matrix = true,
    fcore = (c, Δs) -> 1e-3,
    preserve_core_size = true,
    kwargs...) where TF

    # number of surfaces
    nsurf = length(system.surfaces)

    # --- Set Initial System Variables --- #

    # update system inputs
    update_inputs!(system, kwargs...)

    # check for surface motion
    surface_motion = !isnothing(surface_geometry_history)

    if surface_motion
        # make surface_geometry_history in-place
        if isnothing(surface_geometry_prototype)
            surface_geometry_prototype = surface_geometry_history(zero(TF))
            surface_geometry_history = make_inplace(surface_geometry_history)
        end
        # make surface_velocity_history in-place
        if isnothing(surface_velocity_prototype)
            surface_velocity_prototype = surface_velocity_history(zero(TF))
            surface_velocity_history = make_inplace(surface_velocity_history)
        end
    else
        # zero out surface velocities
        initialize_surface_velocities!(system)
        # calculate initial influence coefficient matrix
        if calculate_influence_matrix
            update_influence_coefficients!(system)
        end
    end

    # create freestream history function (if not defined)
    if isnothing(freestream_history)
        freestream = deepcopy(system.freestream)
        freestream_history = (t) -> freestream
    end

    # initialize undefined wake panels
    initialize_wake_panels!(system)

    # get initial wake circulation product
    Γwl0 = get_wake_circulation(system) .* get_wake_filament_length(system)

    # get initial wake core sizes
    wake_core_size0 = [get_core_size.(system.wakes[isurf]) for isurf = 1:nsurf]

    # --- Define State Variable Shape and Initial Values --- #

    # initial wake vertex locations
    ζw = get_wake_vertices(system)

    # vector of wake vertex locations
    wake_vertex_points = vcat(reshape.(ζw, :)...)

    # matrix of wake vertex locations
    wake_vertex_matrix = copy(reinterpret(reshape, eltype(eltype(wake_vertex_points)), wake_vertex_points))

    # initial state variable vector
    u0 = wake_vertex_matrix

    # indices of vertices corresponding to each surface
    iu2 = cumsum(length.(ζw))
    iu1 = vcat(1, iu2[1:end-1] .- 1)
    iu = [iu1[isurf]:iu2[isurf] for isurf = 1:nsurf]

    # --- Define Parameters --- #
    # store parameters as a named tuple
    p = (
        # number of state variables corresponding to each surface
        iu = iu,
        # time dependent parameters
        surface_geometry_history = surface_geometry_history,
        surface_geometry_prototype = surface_geometry_prototype,
        surface_velocity_history = surface_velocity_history,
        surface_velocity_prototype = surface_velocity_prototype,
        freestream_history = freestream_history,
        # temporary variables from `system`
        AIC = deepcopy(system.AIC),
        w = deepcopy(system.w),
        Γ = deepcopy(system.Gamma),
        surfaces = deepcopy(system.surfaces),
        repeated_points = deepcopy(system.repeated_points),
        wake_shedding_locations = deepcopy(system.wake_shedding_locations),
        wakes = deepcopy(system.wakes),
        Vcp = deepcopy(system.Vcp),
        Vh = deepcopy(system.Vh),
        Vv = deepcopy(system.Vv),
        Vte = deepcopy(system.Vte),
        reference = deepcopy(system.reference),
        symmetric = deepcopy(system.symmetric),
        surface_id = deepcopy(system.surface_id),
        wake_finite_core = deepcopy(system.wake_finite_core),
        iwake0 = deepcopy(system.iwake),
        iwake = deepcopy(system.nwake),
        nwake = deepcopy(system.nwake),
        trailing_vortices = deepcopy(system.trailing_vortices),
        xhat = deepcopy(system.xhat),
        # additional temporary variables
        Γw = similar(Γwl0, TF),
        Γwl = deepcopy(Γwl0),
        Γwl0 = deepcopy(Γwl0),
        wake_core_size0 = deepcopy(wake_core_size0)
    )

    # --- Set Up ODE Solution

    # Wake Shedding Time Pre-Processing
    tshed = unique(tshed) # remove duplicate values
    tshed = sort(unique(tshed)) # sort wake shedding times
    # prevent zero area wake panels from being shed
    if tshed[1] == tspan[1]
        for isurf = 1:nsurf, j = 1:size(system.surfaces[isurf], 2)
            # is the first wake vertices collocated with trailing edge?
            if all(ζw[isurf][1,j] .== bottom_left(system.surfaces[isurf][end,j])) .&
               all(ζw[isurf][1,j+1] .== bottom_right(system.surfaces[isurf][end,j]))
               # if so, make sure that we don't shed a wake panel at tspan[1]
               tshed = tshed[2:end]
               # tshed is now fixed, we can now exit the loop(s)
               break
           end
       end
   end

    # UVLM Wake Progression Ordinary Differential Equation
    ode = uvlm_ode!

    # UVLM Initialization Callback
    init_cb = PresetTimeCallback([tspan[1]], uvlm_init!)

    # UVLM Wake Shedding Callback
    shed_cb = PresetTimeCallback(tshed, uvlm_shed!)

    # # UVLM Result Saving Callback
    # saved_values = SavedValues(TF, NamedTuple{
    #         (:surfaces, :properties, :wakes),
    #         Tuple{
    #             Vector{Matrix{SurfacePanel{TF}}},
    #             Vector{Matrix{PanelProperties{TF}}},
    #             Vector{Matrix{WakePanel{TF}}}
    #         }
    #     }
    # )
    # save_cb = SavingCallback(uvlm_save, saved_values)

    # construct the ODE problem
    prob = ODEProblem(ode, u0, tspan, p; callback = CallbackSet(init_cb, shed_cb))

    return prob
end

# UVLM Wake Progression Ordinary Differential Equation
function uvlm_ode!(du, u, p, t)

    # --- Extract Parameters --- #

    # time-dependent inputs
    surface_geometry_history = p.surface_geometry_history
    surface_velocity_history = p.surface_velocity_history
    freestream_history = p.freestream_history

    # flag for surface motion
    surface_motion = !isnothing(surface_geometry_history)

    # storage for outputs from time-dependent inputs
    TT = eltype(t)
    if TT <: ForwardDiff.Dual && surface_motion
        if eltype(eltype(p.surface_geometry_prototype)) <: SurfacePanel
            # surface panel input
            surface_geometry = [SurfacePanel{TT}.(surf) for surf in p.surface_geometry_prototype]
        else
            # grid input
            surface_geometry = [TT.(surf) for surf in p.surface_geometry_prototype]
        end
        if eltype(eltype(p.surface_velocity_prototype)) <: SurfacePanel
            # surface panel input
            surface_velocities = [SurfacePanel{TT}.(vsurf) for vsurf in p.surface_velocity_prototype]
        else
            # grid input
            surface_velocities = [TT.(vsurf) for vsurf in p.surface_velocity_prototype]
        end
        surfaces = [SurfacePanel{TT}.(surf) for surf in p.surfaces]
        Vcp = [SVector{3,TT}.(V) for V in p.Vcp]
        Vh = [SVector{3,TT}.(V) for V in p.Vh]
        Vv = [SVector{3,TT}.(V) for V in p.Vv]
        Vte = [SVector{3,TT}.(V) for V in p.Vte]
    else
        surface_geometry = p.surface_geometry_prototype
        surface_velocities = p.surface_velocity_prototype
        surfaces = p.surfaces
        Vcp = p.Vcp
        Vh = p.Vh
        Vv = p.Vv
        Vte = p.Vte
    end

    # intermediate variables
    TU = eltype(u)
    if TU <: ForwardDiff.Dual
        AIC = TU.(p.AIC)
        w = TU.(p.w)
        Γ = TU.(p.Γ)
        Γw = TU.(p.Γw)
        Γwl = TU.(p.Γwl)
        wake_shedding_locations = [SVector{3,TU}.(points) for points in p.wake_shedding_locations]
        wakes = [WakePanel{TU}.(wake) for wake in p.wakes]
    else
        AIC = p.AIC
        w = p.w
        Γ = p.Γ
        Γw = p.Γw
        Γwl = p.Γwl
        wake_shedding_locations = p.wake_shedding_locations
        wakes = p.wakes
    end

    # state variable indices corresponding to each surface
    iu = p.iu

    # control parameters
    symmetric = p.symmetric
    surface_id = p.surface_id
    repeated_points = p.repeated_points
    wake_finite_core = p.wake_finite_core
    iwake = p.iwake
    nwake = p.nwake
    trailing_vortices = p.trailing_vortices

    # other parameters
    reference = p.reference
    xhat = p.xhat

    # --- Extract States and State Rates --- #

    u_ζw = [reinterpret(reshape, SVector{3, eltype(u)}, reshape(view(u, :, iu[isurf]), 3, size(wakes[isurf], 1) + 1, size(wakes[isurf], 2) + 1)) for isurf = 1:length(surfaces)]
    du_ζw = [reinterpret(reshape, SVector{3, eltype(du)}, reshape(view(du, :, iu[isurf]), 3, size(wakes[isurf], 1) + 1, size(wakes[isurf], 2) + 1)) for isurf = 1:length(surfaces)]

    # --- Update Time-Dependent Parameters --- #

    if surface_motion
        # update surface geometry
        surface_geometry_history(surface_geometry, t)
        update_surface_panels!(surfaces, surface_geometry)
        # update surface velocities
        surface_velocity_history(surface_velocities, t)
        update_surface_velocities!(Vcp, Vh, Vv, Vte, surface_velocities)
    end

    freestream = freestream_history(t)

    # --- Update Wake Panels --- #

    # set wake shedding locations
    get_wake_shedding_locations!(wake_shedding_locations, u_ζw)

    # set wake vertices
    set_wake_vertices!(wakes, u_ζw)

    # set wake circulation
    vortex_filament_length!(Γw, u_ζw)
    for i = 1:length(Γw)
        if iszero(Γw[i]) || iszero(Γwl[i])
            Γw[i] = 0.0
        else
            Γw[i] = Γwl[i] / Γw[i]
        end
    end
    set_circulation!(wakes, Γw)

    # --- Update Surface Circulation --- #

    if surface_motion
        influence_coefficients!(AIC, surfaces; symmetric, surface_id,
            wake_shedding_locations, trailing_vortices, xhat)
    else
        trailing_edge_coefficients!(AIC, surfaces; symmetric, surface_id,
            wake_shedding_locations, trailing_vortices, xhat)
    end

    normal_velocities!(w, surfaces, wakes, reference, freestream;
        Vcp, symmetric, surface_id, wake_finite_core, iwake, trailing_vortices, xhat)

    circulation!(Γ, AIC, w)

    # --- Calculate Wake Velocities --- #

    wake_velocities!(du_ζw, surfaces, repeated_points,
        wake_shedding_locations, wakes, reference, freestream, Γ; Vte, symmetric,
        surface_id, wake_finite_core, iwake, trailing_vortices, xhat)

    return du
end

# UVLM Discrete State Variable Initialization
function uvlm_init!(integrator)
    # set wake core size to initial wake core size
    wakes = integrator.p.wakes
    wake_core_size = integrator.p.wake_core_size0
    for (isurf, wake) in enumerate(wakes)
        for I in eachindex(wake)
            core_size = integrator.p.wake_core_size0[isurf][I]
            wake[I] = set_core_size(wake[I], core_size)
        end
    end
    # set wake circulation product to initial wake circulation product
    integrator.p.Γwl .= integrator.p.Γwl0
    # set number of wake panels to initial number of wake panels
    integrator.p.iwake .= integrator.p.iwake0
    # indicate that no state variables have been modified
    u_modified!(integrator, false)
end

# UVLM Wake Shedding Callback
function uvlm_shed!(integrator)

    # get state vector, parameters, and current time
    u, p, t = integrator.u, integrator.p, integrator.t

    # --- Extract Parameters --- #

    # time-dependent inputs
    surface_geometry_history = p.surface_geometry_history
    surface_velocity_history = p.surface_velocity_history
    freestream_history = p.freestream_history

    # flag for surface motion
    surface_motion = !isnothing(surface_geometry_history)

    # storage for outputs from time-dependent inputs
    TT = eltype(t)
    if TT <: ForwardDiff.Dual && surface_motion
        if eltype(eltype(p.surface_geometry_prototype)) <: SurfacePanel
            # surface panel input
            surface_geometry = [SurfacePanel{TT}.(surf) for surf in p.surface_geometry_prototype]
        else
            # grid input
            surface_geometry = [TT.(surf) for surf in p.surface_geometry_prototype]
        end
        if eltype(eltype(p.surface_velocity_prototype)) <: SurfacePanel
            # surface panel input
            surface_velocities = [SurfacePanel{TT}.(vsurf) for vsurf in p.surface_velocity_prototype]
        else
            # grid input
            surface_velocities = [TT.(vsurf) for vsurf in p.surface_velocity_prototype]
        end
        surfaces = [SurfacePanel{TT}.(surf) for surf in p.surfaces]
        Vcp = [SVector{3,TT}.(V) for V in p.Vcp]
        Vh = [SVector{3,TT}.(V) for V in p.Vh]
        Vv = [SVector{3,TT}.(V) for V in p.Vv]
        Vte = [SVector{3,TT}.(V) for V in p.Vte]
    else
        surface_geometry = p.surface_geometry_prototype
        surface_velocities = p.surface_velocity_prototype
        surfaces = p.surfaces
        Vcp = p.Vcp
        Vh = p.Vh
        Vv = p.Vv
        Vte = p.Vte
    end

    # intermediate variables
    TU = eltype(u)
    if TU <: ForwardDiff.Dual
        AIC = TU.(p.AIC)
        w = TU.(p.w)
        Γ = TU.(p.Γ)
        Γw = TU.(p.Γw)
        Γwl = TU.(p.Γwl)
        wake_shedding_locations = [SVector{3,TU}.(points) for points in p.wake_shedding_locations]
        wakes = [WakePanel{TU}.(wake) for wake in p.wakes]
    else
        AIC = p.AIC
        w = p.w
        Γ = p.Γ
        Γw = p.Γw
        Γwl = p.Γwl
        wake_shedding_locations = p.wake_shedding_locations
        wakes = p.wakes
    end

    # state variable indices corresponding to each surface
    iu = p.iu

    # control parameters
    symmetric = p.symmetric
    surface_id = p.surface_id
    repeated_points = p.repeated_points
    wake_finite_core = p.wake_finite_core
    iwake = p.iwake
    nwake = p.nwake
    trailing_vortices = p.trailing_vortices

    # other parameters
    reference = p.reference
    xhat = p.xhat

    # --- Extract States --- #

    u_ζw = [reinterpret(reshape, SVector{3, eltype(u)}, reshape(view(u, :, iu[isurf]), 3, size(wakes[isurf], 1) + 1, size(wakes[isurf], 2) + 1)) for isurf = 1:length(surfaces)]

    # --- Update Time-Dependent Parameters --- #

    if surface_motion
        # update surface geometry
        surface_geometry_history(surface_geometry_prototype, t)
        update_surface_panels!(surfaces, surface_geometry_prototype)
        # update surface velocities
        surface_velocity_history(surface_velocity_prototype, t)
        update_surface_velocities!(Vcp, Vh, Vv, Vte, surface_velocity_prototype)
    end

    freestream = freestream_history(t)

    # --- Update Wake Panels --- #

    # set wake shedding locations
    get_wake_shedding_locations!(wake_shedding_locations, u_ζw)

    # set wake vertices
    set_wake_vertices!(wakes, u_ζw)

    # set wake circulation
    vortex_filament_length!(Γw, u_ζw)
    for i = 1:length(Γw)
        if iszero(Γw[i]) || iszero(Γwl[i])
            Γw[i] = 0.0
        else
            Γw[i] = Γwl[i] / Γw[i]
        end
    end
    set_circulation!(wakes, Γw)

    # --- Update Surface Circulation --- #

    if surface_motion
        influence_coefficients!(AIC, surfaces; symmetric, surface_id,
            wake_shedding_locations, trailing_vortices, xhat)
    else
        trailing_edge_coefficients!(AIC, surfaces; symmetric, surface_id,
            wake_shedding_locations, trailing_vortices, xhat)
    end

    normal_velocities!(w, surfaces, wakes, reference, freestream;
        Vcp, symmetric, surface_id, wake_finite_core, iwake, trailing_vortices, xhat)

    circulation!(Γ, AIC, w)

    # --- Shed Wake Panel --- #

    iΓs = 0 # index for accessing surface circulation
    iΓwl = 0 # index for accessing wake circulation product
    for isurf = 1:length(surfaces)
        # number of surface and wake panels for this surface and wake
        Ns = length(surfaces[isurf])
        Nw = length(wakes[isurf])
        # number of spanwise and chordwise panels for this surface and wake
        nc, ns = size(surfaces[isurf])
        nw, ns = size(wakes[isurf])

        # --- modify wake vertices --- #
        # replace oldest row of coordinates with new row of coordinates
        for j = 1:ns
            rte = bottom_left(surfaces[isurf][end, j])
            u_ζw[isurf][end, j] = rte
        end
        rte = bottom_right(surfaces[isurf][end, end])
        u_ζw[isurf][end, end] = rte
        # shift data to put newly shed geometry first
        rowshift!(u_ζw[isurf])

        # --- modify wake circulation product --- #
        # get surface circulation and wake circulation product corresponding to this surface
        vΓs = reshape(view(Γ, iΓs+1:iΓs+Ns), nc, ns)
        vΓwl = reshape(view(Γwl, iΓwl+1:iΓwl+Nw), nw, ns)
        # replace oldest wake circulation product with new wake circulation product
        for j = 1:ns
            # get wake panel corners
            rtl = u_ζw[isurf][1,j]
            rtr = u_ζw[isurf][1,j+1]
            rbl = u_ζw[isurf][2,j]
            rbr = u_ζw[isurf][2,j+1]
            # get wake filament lengths
            lt = norm(rtl - rtr)
            lb = norm(rbl - rbr)
            ll = norm(rtl - rbl)
            lr = norm(rtr - rbr)
            # sum wake filament lengths
            l = (lt + lb + ll + lr)
            # replace oldest wake circulation product with new wake circulation product
            vΓwl[end,j] = vΓs[end,j]*l
        end
        # shift data to put newly shed wake circulation product first
        rowshift!(vΓwl)

        # replace oldest core size with newest core size
        for j = 1:ns
            wakes[isurf][end,j] = set_core_size(wakes[isurf][end,j], get_core_size(surfaces[isurf][end,j]))
        end
        # shift to put newly shed wake panels first
        rowshift!(wakes[isurf])

        # --- move to next surface --- #
        iΓs += Ns
        iΓwl += Nw
    end

    # --- Update Control Parameters --- #

    # increment number of wake panels corresponding to each surface
    for isurf = 1:length(surfaces)
        iwake[isurf] = min(iwake[isurf] + 1, nwake[isurf])
    end
end

# UVLM Circulation Calculation
function uvlm_circulation(u, p, t)
    # --- Extract Parameters --- #

    # time-dependent inputs
    surface_geometry_history = p.surface_geometry_history
    surface_velocity_history = p.surface_velocity_history
    freestream_history = p.freestream_history

    # flag for surface motion
    surface_motion = !isnothing(surface_geometry_history)

    # storage for outputs from time-dependent inputs
    TT = eltype(t)
    if TT <: ForwardDiff.Dual && surface_motion
        if eltype(eltype(p.surface_geometry_prototype)) <: SurfacePanel
            # surface panel input
            surface_geometry = [SurfacePanel{TT}.(surf) for surf in p.surface_geometry_prototype]
        else
            # grid input
            surface_geometry = [TT.(surf) for surf in p.surface_geometry_prototype]
        end
        if eltype(eltype(p.surface_velocity_prototype)) <: SurfacePanel
            # surface panel input
            surface_velocities = [SurfacePanel{TT}.(vsurf) for vsurf in p.surface_velocity_prototype]
        else
            # grid input
            surface_velocities = [TT.(vsurf) for vsurf in p.surface_velocity_prototype]
        end
        surfaces = [SurfacePanel{TT}.(surf) for surf in p.surfaces]
        Vcp = [SVector{3,TT}.(V) for V in p.Vcp]
        Vh = [SVector{3,TT}.(V) for V in p.Vh]
        Vv = [SVector{3,TT}.(V) for V in p.Vv]
        Vte = [SVector{3,TT}.(V) for V in p.Vte]
    else
        surface_geometry = p.surface_geometry_prototype
        surface_velocities = p.surface_velocity_prototype
        surfaces = p.surfaces
        Vcp = p.Vcp
        Vh = p.Vh
        Vv = p.Vv
        Vte = p.Vte
    end

    # intermediate variables
    TU = eltype(u)
    if TU <: ForwardDiff.Dual
        AIC = TU.(p.AIC)
        w = TU.(p.w)
        Γ = TU.(p.Γ)
        Γw = TU.(p.Γw)
        Γwl = TU.(p.Γwl)
        wake_shedding_locations = [SVector{3,TU}.(points) for points in p.wake_shedding_locations]
        wakes = [WakePanel{TU}.(wake) for wake in p.wakes]
    else
        AIC = p.AIC
        w = p.w
        Γ = p.Γ
        Γw = p.Γw
        Γwl = p.Γwl
        wake_shedding_locations = p.wake_shedding_locations
        wakes = p.wakes
    end

    # state variable indices corresponding to each surface
    iu = p.iu

    # control parameters
    symmetric = p.symmetric
    surface_id = p.surface_id
    repeated_points = p.repeated_points
    wake_finite_core = p.wake_finite_core
    iwake = p.iwake
    nwake = p.nwake
    trailing_vortices = p.trailing_vortices

    # other parameters
    reference = p.reference
    xhat = p.xhat

    # --- Extract States --- #

    u_ζw = [reinterpret(reshape, SVector{3, eltype(u)}, reshape(view(u, :, iu[isurf]), 3, size(wakes[isurf], 1) + 1, size(wakes[isurf], 2) + 1)) for isurf = 1:length(surfaces)]

    # --- Update Time-Dependent Parameters --- #

    if surface_motion
        # update surface geometry
        surface_geometry_history(surface_geometry_prototype, t)
        update_surface_panels!(surfaces, surface_geometry_prototype)
        # update surface velocities
        surface_velocity_history(surface_velocity_prototype, t)
        update_surface_velocities!(Vcp, Vh, Vv, Vte, surface_velocity_prototype)
    end

    freestream = freestream_history(t)

    # --- Update Wake Panels --- #

    # set wake shedding locations
    get_wake_shedding_locations!(wake_shedding_locations, u_ζw)

    # set wake vertices
    set_wake_vertices!(wakes, u_ζw)

    # set wake circulation
    vortex_filament_length!(Γw, u_ζw)
    for i = 1:length(Γw)
        if iszero(Γw[i]) || iszero(Γwl[i])
            Γw[i] = 0.0
        else
            Γw[i] = Γwl[i] / Γw[i]
        end
    end
    set_circulation!(wakes, Γw)

    # --- Update Surface Circulation --- #

    if surface_motion
        influence_coefficients!(AIC, surfaces; symmetric, surface_id,
            wake_shedding_locations, trailing_vortices, xhat)
    else
        trailing_edge_coefficients!(AIC, surfaces; symmetric, surface_id,
            wake_shedding_locations, trailing_vortices, xhat)
    end

    normal_velocities!(w, surfaces, wakes, reference, freestream;
        Vcp, symmetric, surface_id, wake_finite_core, iwake, trailing_vortices, xhat)

    circulation!(Γ, AIC, w)

    return Γ
end

# UVLM Force/Moment Calculation Callback
function uvlm_save(u, t, integrator)

    p = integrator.p

    # --- Extract Parameters --- #

    # time-dependent inputs
    surface_geometry_history = p.surface_geometry_history
    surface_velocity_history = p.surface_velocity_history
    freestream_history = p.freestream_history

    # flag for surface motion
    surface_motion = !isnothing(surface_geometry_history)

    # storage for outputs from time-dependent inputs
    TT = eltype(t)
    if TT <: ForwardDiff.Dual && surface_motion
        if eltype(eltype(p.surface_geometry_prototype)) <: SurfacePanel
            # surface panel input
            surface_geometry = [SurfacePanel{TT}.(surf) for surf in p.surface_geometry_prototype]
        else
            # grid input
            surface_geometry = [TT.(surf) for surf in p.surface_geometry_prototype]
        end
        if eltype(eltype(p.surface_velocity_prototype)) <: SurfacePanel
            # surface panel input
            surface_velocities = [SurfacePanel{TT}.(vsurf) for vsurf in p.surface_velocity_prototype]
        else
            # grid input
            surface_velocities = [TT.(vsurf) for vsurf in p.surface_velocity_prototype]
        end
        surfaces = [SurfacePanel{TT}.(surf) for surf in p.surfaces]
        Vcp = [SVector{3,TT}.(V) for V in p.Vcp]
        Vh = [SVector{3,TT}.(V) for V in p.Vh]
        Vv = [SVector{3,TT}.(V) for V in p.Vv]
        Vte = [SVector{3,TT}.(V) for V in p.Vte]
    else
        surface_geometry = p.surface_geometry_prototype
        surface_velocities = p.surface_velocity_prototype
        surfaces = p.surfaces
        Vcp = p.Vcp
        Vh = p.Vh
        Vv = p.Vv
        Vte = p.Vte
    end

    # intermediate variables
    TU = eltype(u)
    if TU <: ForwardDiff.Dual
        AIC = TU.(p.AIC)
        w = TU.(p.w)
        Γ = TU.(p.Γ)
        Γw = TU.(p.Γw)
        Γwl = TU.(p.Γwl)
        wake_shedding_locations = [SVector{3,TU}.(points) for points in p.wake_shedding_locations]
        wakes = [WakePanel{TU}.(wake) for wake in p.wakes]
    else
        AIC = p.AIC
        w = p.w
        Γ = p.Γ
        Γw = p.Γw
        Γwl = p.Γwl
        wake_shedding_locations = p.wake_shedding_locations
        wakes = p.wakes
    end

    # state variable indices corresponding to each surface
    iu = p.iu

    # control parameters
    symmetric = p.symmetric
    surface_id = p.surface_id
    repeated_points = p.repeated_points
    wake_finite_core = p.wake_finite_core
    iwake = p.iwake
    nwake = p.nwake
    trailing_vortices = p.trailing_vortices

    # other parameters
    reference = p.reference
    xhat = p.xhat

    # --- Extract States --- #

    u_ζw = [reinterpret(reshape, SVector{3, eltype(u)}, reshape(view(u, :, iu[isurf]), 3, size(wakes[isurf], 1) + 1, size(wakes[isurf], 2) + 1)) for isurf = 1:length(surfaces)]

    # --- Update Time-Dependent Parameters --- #

    if surface_motion
        # update surface geometry
        surface_geometry_history(surface_geometry_prototype, t)
        update_surface_panels!(surfaces, surface_geometry_prototype)
        # update surface velocities
        surface_velocity_history(surface_velocity_prototype, t)
        update_surface_velocities!(Vcp, Vh, Vv, Vte, surface_velocity_prototype)
    end

    freestream = freestream_history(t)

    # --- Update Wake Panels --- #

    # set wake shedding locations
    get_wake_shedding_locations!(wake_shedding_locations, u_ζw)

    # set wake vertices
    set_wake_vertices!(wakes, u_ζw)

    # set wake circulation
    vortex_filament_length!(Γw, u_ζw)
    for i = 1:length(Γw)
        if iszero(Γw[i]) || iszero(Γwl[i])
            Γw[i] = 0.0
        else
            Γw[i] = Γwl[i] / Γw[i]
        end
    end
    set_circulation!(wakes, Γw)

    # --- Update Surface Circulation --- #

    if surface_motion
        influence_coefficients!(AIC, surfaces; symmetric, surface_id,
            wake_shedding_locations, trailing_vortices, xhat)
    else
        trailing_edge_coefficients!(AIC, surfaces; symmetric, surface_id,
            wake_shedding_locations, trailing_vortices, xhat)
    end

    normal_velocities!(w, surfaces, wakes, reference, freestream;
        Vcp, symmetric, surface_id, wake_finite_core, iwake, trailing_vortices, xhat)

    circulation!(Γ, AIC, w)

    # --- Update Surface Circulation Rate --- #

    # Γdot = ∂Γ∂u * dudt + ∂Γ∂p * dpdt
    ∂Γ∂u = ForwardDiff.jacobian((u) -> uvlm_circulation(u, p, t), reshape(u, :))
    dudt = reshape(integrator(t, Val{1}), :)
    ∂Γ∂p_dpdt = ForwardDiff.derivative((t) -> uvlm_circulation(u, p, t), t)
    Γdot = ∂Γ∂u*dudt + ∂Γ∂p_dpdt

    # --- Calculate Near Field Properties --- #
    properties = near_field_properties(surfaces, wakes, reference, freestream, Γ,
        Γdot; wake_shedding_locations, Vh, Vv, symmetric, surface_id,
        wake_finite_core, iwake, trailing_vortices, xhat)

    # create copy of surfaces for outputs
    surfaces = [copy(surface) for surface in surfaces]

    # create copy of wakes for outputs, with only defined wake panels
    wakes = [wakes[isurf][1:iwake[isurf], :] for isurf = 1:length(wakes)]

    # return named tuple of outputs
    return (surfaces=surfaces, properties=properties, wakes=wakes)
end

function initialize_surface_velocities!(system)
    for isurf = 1:length(system.surfaces)
        system.Vcp[isurf] .= Ref(@SVector zeros(3))
        system.Vh[isurf] .= Ref(@SVector zeros(3))
        system.Vv[isurf] .= Ref(@SVector zeros(3))
        system.Vte[isurf] .= Ref(@SVector zeros(3))
    end
    return system
end

"""
    make_inplace(out, f)

Make a
"""
function make_inplace(out, f)
    if isnothing(out)
        # function is out-of-place
        f! = function(out, t)
            out .= f(t)
        end
    else
        # function is already in-place
        f! = f
    end
    return f!
end

# ----------------------------------- #


# ------------------------------------ #

# NOTE: These functions were taken from CoupledSystems, and may be used directly
# from there in the future

"""
    vector_length(var)

Return the length of a variable when expressed as a vector
"""
vector_length(x) = add_vector_length(x, 0)

@inline add_vector_length(x::Number, nx) = nx + 1
@inline add_vector_length(x::NTuple{N,T}, nx) where {N,T<:Number} = nx + length(x)
@inline add_vector_length(x::AbstractArray{T,N}, nx) where {T<:Number,N} = nx + length(x)
@inline function add_vector_length(x::AbstractArray, nx)
    for xi in x
        nx = add_vector_length(xi, nx)
    end
    return nx
end
@inline function add_vector_length(x::Tuple, nx)
    nx = add_vector_length(first(x), nx)
    return add_vector_length(Base.tail(x), nx)
end
@inline add_vector_length(x::Tuple{}, nx) = nx


"""
    combine(x)

Combine the data in `x` into a vector.
"""
combine(x) = [x]
function combine(x::Union{<:Tuple, <:AbstractArray})
    x = collect(Iterators.flatten(x))
    while any(x->typeof(x) <: Union{<:Tuple, <:AbstractArray}, x)
        x = collect(Iterators.flatten(x))
    end
    return x
end
combine(x::Tuple{}) = Float64[]

"""
    combine!(v, x)

Combine the data in `x` into the vector `v`.
"""
@inline function combine!(v, x)
    vo, v = offset_combine!(v, x, firstindex(v))
    return v
end
@inline combine!(v, x::Tuple{}) = v

"""
    offset_combine!(v, x, vo)

Copies the data in `x` into the vector `v` starting at offset `vo`.  Return the
new offset `vo` and the modified vector `v`.
"""
offset_combine!

@inline function offset_combine!(v, x::Number, vo)
    nx = length(x)
    copyto!(v, vo, x, 1, nx)
    vo += nx
    return vo, v
end
@inline function offset_combine!(v, x::NTuple{N,T}, vo) where {N, T<:Number}
    nx = length(x)
    copyto!(v, vo, x, 1, nx)
    vo += nx
    return vo, v
end
@inline function offset_combine!(v, x::AbstractArray{T,N}, vo) where {T<:Number, N}
    nx = length(x)
    copyto!(v, vo, x, 1, nx)
    vo += nx
    return vo, v
end
@inline function offset_combine!(v, x::AbstractArray, vo)
    for xi in x
        vo, v = offset_combine!(v, xi, vo)
    end
    return vo, v
end
@inline function offset_combine!(v, x::Tuple, vo)
    vo, v = offset_combine!(v, first(x), vo)
    return offset_combine!(v, Base.tail(x), vo)
end
@inline offset_combine!(v, x::Tuple{}, vo) = vo, v

"""
    separate(x, v)

Return the data in the vector `v` with the shape of `x`.
"""
separate(x, v) = separate!(deepcopy(x), v)
separate(x::Tuple{}, v) = v

"""
    separate!(x, v)

Return the data in the vector `v` with the shape of `x`.  Modifies `x` in order
to avoid allocations, if possible.
"""
@inline function separate!(x, v)
    vo, x = offset_separate!(x, v, firstindex(v))
    return x
end
separate!(x::Tuple{}, v) = v

"""
    offset_separate!(x, v, vo)

Return the data in the vector `v` with the shape of `x` starting at offset `vo`.
Modify `x` in order to avoid allocations, if possible.
"""
offset_separate!

@inline function offset_separate!(x::Number, v, vo)
    x = v[vo]
    vo += 1
    return vo, x
end
@inline function offset_separate!(x::NTuple{N,T}, v, vo) where {N, T<:Number, TV}
    x = NTuple{N,eltype(v)}(view(v, vo : vo + N - 1))
    vo += N
    return vo, x
end
@inline function offset_separate!(x::AbstractArray{T,N}, v, vo) where {T<:Number, N, TV}
    nx = length(x)
    x = reshape(view(v, vo : vo + nx - 1), size(x))
    vo += nx
    return vo, x
end
@inline function offset_separate!(x::StaticArray{S, T, N}, v, vo) where {S, T<:Number, N}
    nx = length(x)
    x = SArray{S}(view(v, vo : vo + nx - 1))
    vo += nx
    return vo, x
end
function offset_separate!(x::AbstractArray, v, vo)
    if ismutable(x)
        # modify in place (if possible)
        for i = 1:length(x)
            # extract new element type
            vo, xi = offset_separate!(x[i], v, vo)
            # update element type of x if necessary
            TE = promote_type(typeof(xi), eltype(x))
            if !(TE <: eltype(x))
                x = TE.(x)
            end
            # store result
            x[i] = xi
        end
    else
        # convert inputs to tuple
        xt = Tuple(x)
        # get outputs as a tuple
        vo, xt = offset_separate!(xt, v, vo)
        # update element type of x if necessary
        TE = promote_type(eltype(x), typeof.(xt)...)
        if !(TE <: eltype(x))
            x = TE.(x)
        end
        # convert outputs to correct type
        x = typeof(x)(xt)
    end
    return vo, x
end
function offset_separate!(x::Tuple, v, vo)
    vo, xi = offset_separate!(first(x), v, vo)
    vo, x = offset_separate!(Base.tail(x), v, vo)
    return vo, (xi, x...)
end
offset_separate!(x::Tuple{}, v, vo) = vo, ()

# ----------------------------------------- #

"""
    rowshift!(A)

Circularly shifts the rows of a matrix down one row.
"""
function rowshift!(A)

    ni, nj = size(A)

    for j = 1:nj
        tmp = A[ni,j]
        for i = ni:-1:2
            A[i,j] = A[i-1,j]
        end
        A[1,j] = tmp
    end

    return A
end
