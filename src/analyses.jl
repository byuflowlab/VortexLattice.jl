"""
    steady_analysis(surfaces, reference, freestream; kwargs...)

Perform a steady vortex lattice method analysis.  Return an object of type
[`System`](@ref) which contains the system inputs and states.

# Arguments
 - `surfaces`: Surface geometry. Represented as either a:
     - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
    or
     - Vector of matrices of shape (nc, ns) containing surface panels (see [`SurfacePanel`](@ref))
    where `nc` is the number of chordwise panels and `ns` is the number of
    spanwise panels.
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))

# Keyword Arguments
 - `wakes`: Wakes. Represented as either a:
     - Vector of grids of shape (3, nw+1, ns+1) which represent wake geometry
    or
     - Vector of matrices of shape (nw, ns) containing wake panels
       (see [`WakePanel`](@ref))
    where `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels.  Defaults to no wake panels.
 - `wake_circulation`: Vector containing the initial wake circulation for all
    wake panels. Overwrites wake circulation provided by keyword argument `wakes`.
    Defaults to no wake circulation.
 - `additional_velocity`: Function of the form `V = f(r)` which defines the
    additional velocity `V` at location `r` in the global coordinate frame.
    Defaults to no additional velocity field.
 - `symmetric`: Flags indicating whether the geometry for each surface is
    symmetric.  Note that applying symmetry to surfaces is only theoretically
    valid when the freestream conditions are symmetric as well.  Defaults to
    `false` for each surface.
 - `surface_id`: Surface ID for each surface.  Defaults to 1:length(surfaces).
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.  Defaults to `true` for each wake.
 - `iwake`: Number of wake panels in the chordwise direction to use for each
    surface. Defaults to all provided wake panels
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from
    each surface.  Defaults to `true` for steady simulations.
 - `xhat`: Direction in which trailing vortices are shed.  Defaults to [1,0,0].
    Only used when using the fixed-wake model.  (`free_wake = false`)
 - `near_field_analysis`: Flag indicating whether to perform a near field analysis
    to obtain surface panel properties.  Defaults to `true`.
 - `derivatives`: Flag indicating whether to calculate derivatives of the near
    field properties with respect to the freestream variables.  Defaults to `true`
 - `fcore`: Function for setting the finite core size when generating surface/wake
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width.
 - `preserve_core_size`: Flag indicating whether the finite core size should be
    preserved from the previous set of surface/wake panels (if applicable).
    Defaults to `false`.
 - `free_wake`: Flag indicating whether the wake geometry (and trailing vortices)
    should be free to move. Defaults to `false`.
 - `dt`: Time step.  Used for defining the separation between wake vertices when
    using the free-wake model (`free_wake = true`)
"""
function steady_analysis(surfaces, reference, freestream; kwargs...)

    # pre-allocate system storage
    system = System(surfaces, reference, freestream)

    return steady_analysis!(system; kwargs..., calculate_influence_matrix = true)
end

"""
    steady_analysis!(system::System; kwargs...)

Perform a steady analysis on the vortex lattice system in `system`.

# Keyword Arguments
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix should be calculated. Defaults to `true`.  Re-using
    the same AIC matrix can reduce calculation times, but is only valid when the
    underlying geometry has not changed.
 - `near_field_analysis`: Flag indicating whether to perform a near field analysis
    to obtain surface panel properties.  Defaults to `true`.
 - `derivatives`: Flag indicating whether to calculate derivatives of the near
    field properties with respect to the freestream variables.  Defaults to `true`
 - `fcore`: Function for setting the finite core size when generating surface/wake
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width.
 - `preserve_core_size`: Flag indicating whether the finite core size should be
    preserved from the previous set of surface/wake panels (if applicable).
    Defaults to `false`.
 - `free_wake`: Flag indicating whether the wake geometry (and trailing vortices)
    should be free to move. Defaults to `false`.
 - `dt`: Time step.  Used for defining the separation between wake vertices when
    using the free-wake model (`free_wake = true`)

The following keyword arguments may be provided to update the inputs to the
vortex lattice system prior to performing the steady analysis.

# Additional Keyword Arguments
 - `surfaces`: Surface geometry. Represented as either a:
     - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
    or
     - Vector of matrices of shape (nc, ns) containing surface panels (see [`SurfacePanel`](@ref))
    where `nc` is the number of chordwise panels and `ns` is the number of
    spanwise panels.
 - `wakes`: Wake panels. Represented as either a:
     - Vector of grids of shape (3, nw+1, ns+1) which represent wake geometry
    or
     - Vector of matrices of shape (nw, ns) containing wake panels (see [`WakePanel`](@ref))
    where `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels.
 - `wake_circulation`: Vector containing the initial wake circulation for all
    wake panels. Overwrites wake circulation provided by keyword argument `wakes`
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))
 - `additional_velocity`: Function of the form `V = f(r)` which defines the
    additional velocity `V` at location `r` in the global coordinate frame. By
    default, no additional velocity field is assumed.
 - `symmetric`: Flags indicating whether the geometry for each surface is
    symmetric.  Note that applying symmetry to surfaces is only theoretically
    valid when the freestream conditions are symmetric as well.
 - `surface_id`: Surface ID for each surface.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.
 - `iwake`: Number of wake panels in the chordwise direction to use for each surface.
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from
    each surface.
 - `xhat`: Direction in which trailing vortices are shed.  Only used if the
    wakes are fixed (`fixed_wake = true`).
"""
function steady_analysis!(system::System; additional_velocity =
    system.additional_velocity, kwargs...)

    # store new additional velocity function
    update_inputs!(system; additional_velocity)

    return steady_analysis!(system, additional_velocity; kwargs...)

end

function steady_analysis!(system::System, additional_velocity;
    calculate_influence_matrix = true, near_field_analysis = true, derivatives = true,
    free_wake = false, dt = 1.0, kwargs...)

    # update system inputs using keyword arguments
    update_inputs!(system; kwargs...)

    # initialize system variables for a steady simulation
    initialize_steady_simulation!(system)

    # set freestream direction (if free wake)
    if free_wake
        tmp = freestream_velocity(system.freestream)
        system.xhat = tmp ./ norm(tmp)
    end

    # unpack system inputs
    surfaces = system.surfaces
    wakes = system.wakes
    reference = system.reference
    freestream = system.freestream
    Vcp = system.Vcp
    Vh = system.Vh
    Vv = system.Vv

    # unpack system internal variables
    AIC = system.AIC
    w = system.w
    dw = system.dw
    Γ = system.Gamma
    dΓ = system.dGamma
    Γdot = system.Gammadot
    properties = system.properties
    dproperties = system.dproperties
    wake_shedding_locations = system.wake_shedding_locations
    Vw = system.Vw

    # unpack system control variables
    repeated_points = system.repeated_points
    symmetric = system.symmetric
    surface_id = system.surface_id
    wake_finite_core = system.wake_finite_core
    iwake = system.iwake
    nwake = system.nwake
    eta = system.eta
    trailing_vortices = system.trailing_vortices
    xhat = system.xhat

    # determine if any wake panels exist
    wake_panels = any(iw -> iw > 0, iwake)

    # calculate/re-calculate AIC matrix
    if calculate_influence_matrix
        influence_coefficients!(AIC, surfaces, symmetric, surface_id,
            trailing_vortices, xhat, wake_shedding_locations, wakes,
            wake_finite_core, iwake)
    end

    # calculate RHS
    if derivatives
        steady_normal_velocities_and_derivatives!(w, dw, surfaces, reference,
            freestream, additional_velocity; symmetric, surface_id, trailing_vortices, xhat)
    else
        steady_normal_velocities!(w, surfaces, reference, freestream,
            additional_velocity; symmetric, surface_id, trailing_vortices, xhat)
    end

    # solve for the free wake vertex locations
    if free_wake
        # calculate wake shedding locations
        wake_shedding_locations!(wake_shedding_locations, surfaces,
            reference, freestream, additional_velocity, 0, eta*dt)
        # calculate free wake vertices
        if wake_panels
            # update wake panels to account for wake shedding locations
            update_wake_panels!(wakes, wake_shedding_locations, iwake)
            # solve residual expression for wake vertex locations
            free_wake_vertices!(wakes, surfaces, repeated_points,
                wake_shedding_locations, reference, freestream, AIC, w, Γ; symmetric,
                surface_id, wake_finite_core, iwake, trailing_vortices, xhat, dt)
        end
        # update trailing edge coefficients of AIC matrix
        trailing_coefficients!(AIC, surfaces, symmetric, surface_id,
            trailing_vortices, xhat, wake_shedding_locations, wakes,
            wake_finite_core, iwake)
    end

    # solve for the circulation distribution
    if derivatives
        circulation_and_derivatives!(Γ, dΓ, AIC, w, dw)
    else
        circulation!(Γ, AIC, w)
    end

    # set wake panel circulation equal to the trailing edge circulation
    if wake_panels
        set_steady_circulation!(wakes, surfaces, Γ; iwake)
    end

    # perform a near field analysis
    if near_field_analysis
        if derivatives
            near_field_properties_and_derivatives!(properties, dproperties,
                surfaces, wakes, reference, freestream, additional_velocity,
                Γ, dΓ, Γdot;
                wake_shedding_locations, Vh, Vv, symmetric, surface_id,
                wake_finite_core, iwake, trailing_vortices, xhat)
        else
            near_field_properties!(properties, surfaces, wakes, reference,
                freestream, additional_velocity, Γ, Γdot; wake_shedding_locations, Vh, Vv, symmetric,
                surface_id, wake_finite_core, iwake, trailing_vortices, xhat)
        end
    end

    # return the modified system
    return system
end

"""
    unsteady_analysis(surfaces, reference, freestream, tvec; kwargs...)

Perform a unsteady vortex lattice method analysis.  Return an object of type
[`System`](@ref) containing the final system state, a vector of times at which
properties are saved, a matrix of surface panels (see [`SurfacePanel`](@ref) for
each surface at each saved time step, a matrix of surface panel properties
(see [`PanelProperties`](@ref)) for each surface at each saved time step, and a
matrix of wake panels (see [`WakePanel`](@ref)) for each surface at each saved
time step.

# Arguments
 - `surfaces`: Surface geometry at the beginning of the simulation.
    Represented as either a:
     - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
    or
     - Vector of matrices of shape (nc, ns) containing surface panels (see [`SurfacePanel`](@ref))
       where `nc` is the number of chordwise panels and `ns` is the number of
       spanwise panels.
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters for each time step (see [`Freestream`](@ref))
 - `tvec`: Time vector

# Keyword Arguments
 - `circulation`: Vector containing the initial circulation of all surface
    panels in the system.  Defaults to `zeros(N)` where `N` is the total number
    of surface panels in `surfaces`.
 - `wakes`: Vector of initial wakes corresponding to each surface, represented
    by matrices of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where
    `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels. Defaults to no wake panels for each surface
 - `wake_circulation`: Vector containing the initial wake circulation for all
    wake panels. Overwrites wake circulation provided by keyword argument `wakes`.
    Defaults to no wake circulation.
 - `additional_velocity`: Function of the form `V = f(r, t)` which defines the
    initial additional velocity `V` at location `r` in the global coordinate
    frame. Defaults to no additional velocity field.
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
    surface.  Defaults to `length(tvec)-1` for all surfaces.
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
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces for the final
    time step. Defaults to `true`.
 - `derivatives`: Flag indicating whether the derivatives with respect to the
    freestream variables should be calculated for the final time step. Defaults
    to `true`.
 - `fcore`: function for setting the finite core size when generating surface
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width. Only used if surfaces
    are provided as grids.
 - `preserve_core_size`: Flag indicating whether the finite core size should be
    preserved from the previous set of panels when updating surface geometries.
    Defaults to `true`.
 - `save_start`: Indicates whether to save the initial conditions. Defaults to `true`.
 - `save_steps`: Indices of the intermediate time steps at which to save the
    time history, defaults to `1:length(tvec)-1`.
 - `save_end`: Indicates whether to save the final time step.  Defaults to `true`.
"""
unsteady_analysis

function unsteady_analysis(surfaces, reference, freestream, tvec; kwargs...)

    # pre-allocate system storage
    system = System(surfaces, reference, freestream, length(tvec)-1)

    return unsteady_analysis!(system, tvec; kwargs..., calculate_influence_matrix = true)
end

"""
    unsteady_analysis!(system, tvec; kwargs...)

Perform an unsteady analysis on the vortex lattice system in `system` given the
vector of times in `tvec`.

# Keyword Arguments
 - `surfaces`: Surface geometry at the beginning of the simulation.
    Represented as either a:
     - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
    or
     - Vector of matrices of shape (nc, ns) containing surface panels (see [`SurfacePanel`](@ref))
       where `nc` is the number of chordwise panels and `ns` is the number of
       spanwise panels.
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters for each time step (see [`Freestream`](@ref))
 - `tvec`: Time vector

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
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces for the final
    time step. Defaults to `true`.
 - `derivatives`: Flag indicating whether the derivatives with respect to the
    freestream variables should be calculated for the final time step. Defaults
    to `true`.
 - `fcore`: function for setting the finite core size when generating surface
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width. Only used if surfaces
    are provided as grids.
 - `preserve_core_size`: Flag indicating whether the finite core size should be
    preserved from the previous set of panels when updating surface geometries.
    Defaults to `true`.
 - `save_start`: Indicates whether to save the initial conditions. Defaults to `true`.
 - `save_steps`: Indices of the intermediate time steps at which to save the
    time history, defaults to `1:length(tvec)-1`.
 - `save_end`: Indicates whether to save the final time step.  Defaults to `true`.

 The following keyword arguments may be provided to update the inputs to the
 vortex lattice system prior to performing the steady analysis.

# Additional Keyword Arguments
 - `surfaces`: Initial surface geometry. Represented as either a:
     - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
    or
     - Vector of matrices of shape (nc, ns) containing surface panels (see [`SurfacePanel`](@ref))
    where `nc` is the number of chordwise panels and `ns` is the number of
    spanwise panels.
 - `wakes`: Initial wake panels. Represented as either a:
     - Vector of grids of shape (3, nw+1, ns+1) which represent wake geometry
    or
     - Vector of matrices of shape (nw, ns) containing wake panels (see [`WakePanel`](@ref))
    where `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels.
 - `wake_circulation`: Vector containing the initial wake circulation for all
    wake panels. Overwrites wake circulation provided by keyword argument `wakes`
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))
 - `additional_velocity`: Function of the form `V = f(r, t)` which defines the
    additional velocity `V` at location `r` in the global coordinate frame. By
    default, no additional velocity field is assumed.
 - `symmetric`: Flags indicating whether the geometry for each surface is
    symmetric.  Note that applying symmetry to surfaces is only theoretically
    valid when the freestream conditions are symmetric as well.
 - `surface_id`: Surface ID for each surface.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.
 - `iwake`: Initial number of wake panels in the chordwise direction to use for
    each surface.
 - `nwake`: Maximum number of wake panels in the chordwise direction to use for
    each surface.
 - `eta`: Time step fraction used to define separation between trailing
    edge and wake shedding location.  Typical values range from 0.2 to 0.3.
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from
    each surface.
 - `xhat`: Direction in which trailing vortices are shed (if shed).
"""
function unsteady_analysis!(system, tvec; additional_velocity =
    system.additional_velocity, kwargs...)

    update_inputs!(system; additional_velocity)

    return unsteady_analysis!(system, additional_velocity, tvec; kwargs...)
end

function unsteady_analysis!(system, additional_velocity, tvec;
    surface_geometry_history = nothing,
    surface_geometry_prototype = nothing,
    surface_velocity_history = nothing,
    surface_velocity_prototype = nothing,
    freestream_history = nothing,
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true,
    fcore = (c, Δs) -> 1e-3,
    preserve_core_size = true,
    save_start = true,
    save_steps = 1:length(tvec)-1,
    save_end = true,
    kwargs...)

    # get floating point type
    TF = eltype(system)

    # update system inputs using keyword arguments
    update_inputs!(system; kwargs...)

    # unpack system inputs
    surfaces = system.surfaces
    wakes = system.wakes
    reference = system.reference
    freestream = system.freestream
    Γ = system.Gamma
    Γdot = system.Gammadot

    # unpack system internal variables
    AIC = system.AIC
    w = system.w
    dw = system.dw
    dΓ = system.dGamma
    properties = system.properties
    dproperties = system.dproperties
    wake_shedding_locations = system.wake_shedding_locations
    Vcp = system.Vcp
    Vh = system.Vh
    Vv = system.Vv
    Vw = system.Vw

    # unpack system control variables
    repeated_points = system.repeated_points
    symmetric = system.symmetric
    surface_id = system.surface_id
    wake_finite_core = system.wake_finite_core
    iwake = system.iwake
    nwake = system.nwake
    eta = system.eta
    trailing_vortices = system.trailing_vortices
    xhat = system.xhat

    # number of surfaces
    nsurf = length(surfaces)

    # --- Initialize Motion Dependent Inputs/Parameters --- #

    # check for surface motion
    surface_motion = !isnothing(surface_geometry_history)

    # Motion Dependent Pre-Processing
    if surface_motion
        # make in-place geometry function
        if isnothing(surface_geometry_prototype)
            surface_geometry_prototype = surface_geometry_history(zero(TF))
            surface_geometry_history = make_inplace(surface_geometry_history)
        end

        # make in-place velocity function
        if isnothing(surface_velocity_prototype)
            surface_velocity_prototype = surface_velocity_history(zero(TF))
            surface_velocity_history = make_inplace(surface_velocity_history)
        end

        # set initial surface geometries
        surface_geometry_history(surface_geometry_prototype, tvec[1])
        update_surface_panels!(surfaces, surface_geometry_prototype)

        # set initial surface velocities
        surface_velocity_history(surface_velocity_prototype, tvec[1])
        update_surface_velocities!(Vcp, Vh, Vv, surface_velocity_prototype)

    else

        # zero out surface velocities
        for isurf = 1:nsurf
            Vcp[isurf] .= Ref(@SVector zeros(3))
            Vh[isurf] .= Ref(@SVector zeros(3))
            Vv[isurf] .= Ref(@SVector zeros(3))
        end

        # pre-calculate initial influence coefficient matrix
        if calculate_influence_matrix
            influence_coefficients!(AIC, surfaces, symmetric, surface_id,
                trailing_vortices, xhat)
        end

    end

    # create freestream history function, if not defined
    if isnothing(freestream_history)
        freestream_history = (t) -> freestream
    end

    # --- Begin Simulation --- #

    # initialize solution history for each time step
    nsave = save_start + length(save_steps) + save_end
    time_history = Vector{TF}(undef, nsave)
    surface_history = Vector{Vector{Matrix{SurfacePanel{TF}}}}(undef, nsave)
    property_history = Vector{Vector{Matrix{PanelProperties{TF}}}}(undef, nsave)
    wake_history = Vector{Vector{Matrix{WakePanel{TF}}}}(undef, nsave)

    # current index for saving
    isave = 1

    # save initial timepoint, if specified
    if save_start
        # calculate near field properties
        near_field_properties!(properties, surfaces, wakes, reference,
            freestream, additional_velocity, Γ, Γdot; wake_shedding_locations,
            Vh, Vv, symmetric, surface_id, wake_finite_core, iwake,
            trailing_vortices, xhat)

        # save state
        time_history[isave] = tvec[1]
        surface_history[isave] = [copy(system.surfaces[isurf]) for isurf = 1:nsurf]
        property_history[isave] = [copy(system.properties[isurf]) for isurf = 1:nsurf]
        wake_history[isave] = [system.wakes[isurf][1:iwake[isurf], :] for isurf = 1:nsurf]

        # increment index for saving
        isave += 1
    end

    # set most recent time for which the surface circulation (and its rate) is known
    tprev = tvec[1]

    # loop through all time steps
    for it = 1 : length(tvec) - 1

        # times corresponding to this time step
        t1 = tvec[it] # start of time step
        t2 = tvec[it + 1] # end of time step
        dt = t2 - t1 # time step size
        teta = t1 + eta*dt # time at which properties are evaluated

        # wake shedding locations at `t + eta*dt`
        wake_shedding_locations!(wake_shedding_locations, surfaces,
            reference, freestream, additional_velocity, t1, teta)

        # update wake panels to account for new wake shedding locations
        update_wake_panels!(wakes, wake_shedding_locations, iwake)

        # surface geometries/velocities at `t + eta*dt`
        if surface_motion
            # update surface geometry
            surface_geometry_history(surface_geometry_prototype, teta)
            update_surface_panels!(surfaces, surface_geometry_prototype)
            # update surface velocities
            surface_velocity_history(surface_velocity_prototype, teta)
            update_surface_velocities!(Vcp, Vh, Vv, surface_velocity_prototype)
        end

        # influence coefficient matrix at `t + eta*dt`
        if surface_motion
            influence_coefficients!(AIC, surfaces, symmetric, surface_id,
                trailing_vortices, xhat, wake_shedding_locations)
        else
            trailing_coefficients!(AIC, surfaces, symmetric, surface_id,
                trailing_vortices, xhat, wake_shedding_locations)
        end

        # normal velocities at `t + eta*dt`
        unsteady_normal_velocities!(w, surfaces, wakes, reference, freestream,
            additional_velocity, Vcp; symmetric, surface_id, wake_finite_core,
            iwake, trailing_vortices, xhat)

        # initialize circulation rate calculation (Γdot_init = 2/dt * Γ + Γdot)
        Γdot .= 2 ./ (teta - tprev) .* Γ .+ Γdot

        # circulation at `t + eta*dt`
        circulation!(Γ, AIC, w)

        # calculate circulation rate at `t + eta*dt`
        # Γdot = 2/(t - t_prev) * (Γ - Γ_prev) - Γdot_prev
        Γdot .= 2 ./ (teta - tprev) .* Γ .- Γdot

        # set most recent time for which the surface circulation (and its rate) is known
        tprev = teta

        # wake velocities at `t + eta*dt`
        wake_velocities!(Vw, surfaces, repeated_points,
            wake_shedding_locations, wakes, reference, freestream,
            additional_velocity, Γ; symmetric, surface_id, wake_finite_core,
            iwake, trailing_vortices, xhat)

        # println("velocities")
        # display(Vw[1][1:iwake[1]+1, :])

        if it in save_steps
            # near field properties at `t + eta*dt`
            near_field_properties!(properties, surfaces, wakes, reference,
                freestream, additional_velocity, Γ, Γdot; wake_shedding_locations,
                Vh, Vv, symmetric, surface_id, wake_finite_core, iwake,
                trailing_vortices, xhat)

            # save state at `t + eta*dt`
            time_history[isave] = teta
            surface_history[isave] = [copy(system.surfaces[isurf]) for isurf = 1:nsurf]
            property_history[isave] = [copy(system.properties[isurf]) for isurf = 1:nsurf]
            wake_history[isave] = [system.wakes[isurf][1:iwake[isurf], :] for isurf = 1:nsurf]

            # increment index for saving
            isave += 1
        end

        # translate wake using the explicit Euler method
        translate_wake!(wakes, wake_shedding_locations, Vw, dt; iwake)

        # NOTE: This time stepping scheme is modified from the classical explicit
        # Euler method because properties from `t1 + eta*(t2 - t1)` are used to
        # calculate wake velocities rather than properties at `t1`.  This is
        # standard practice for unsteady vortex lattice methods.

        # println("vertices")
        # display(get_vertices(wake_shedding_locations[1], wakes[1][1:iwake[1],:]))

        # shed an additional wake panel
        shed_wake!(wakes, wake_shedding_locations, surfaces, Γ; iwake)

        # increment wake panel counter for each surface
        for isurf = 1:nsurf
            if iwake[isurf] < nwake[isurf]
                iwake[isurf] += 1
            end
        end

        # println("vertices")
        # display(get_vertices(wake_shedding_locations[1], wakes[1][1:iwake[1],:]))

        # set surface geometry and velocities at `t + dt`
        if surface_motion
            # update surface geometry
            surface_geometry_history(surface_geometry_prototype, t2)
            update_surface_panels!(surfaces, surface_geometry_prototype)
            # update surface velocities
            surface_velocity_history(surface_velocity_prototype, t2)
            update_surface_velocities!(Vcp, Vh, Vv, surface_velocity_prototype)
        end

    end

    # --- Calculate Final Properties --- #

    # calculate/re-calculate AIC matrix at the end of the time interval
    if surface_motion
        influence_coefficients!(AIC, surfaces, symmetric, surface_id,
            trailing_vortices, xhat)
    else
        trailing_coefficients!(AIC, surfaces, symmetric, surface_id,
            trailing_vortices, xhat)
    end

    # calculate/re-calculate normal velocities at the end of the time interval
    if derivatives
        unsteady_normal_velocities_and_derivatives!(w, dw, surfaces, wakes, reference,
            freestream, additional_velocity, Vcp; symmetric, surface_id, wake_finite_core,
            iwake, trailing_vortices, xhat)
    else
        unsteady_normal_velocities!(w, surfaces, wakes, reference, freestream,
            additional_velocity, Vcp; symmetric, surface_id, wake_finite_core,
            iwake, trailing_vortices, xhat)
    end

    # initialize circulation rate calculation (Γdot_init = 2/dt * Γ + Γdot)
    Γdot .= 2 ./ (tvec[end] - tprev) .* Γ .+ Γdot

    # calculate circulation at the end of the time interval
    if derivatives
        circulation_and_derivatives!(Γ, dΓ, AIC, w, dw)
    else
        circulation!(Γ, AIC, w)
    end

    # calculate circulation rate at the end of the time interval
    # Γdot = 2/(t - t_prev) * (Γ - Γ_prev) - Γdot_prev
    Γdot .= 2 ./ (tvec[end] - tprev) .* Γ .- Γdot

    # set most recent time for which the surface circulation (and its rate) is known
    tprev = tvec[end]

    # perform a near field analysis at the end of the time interval
    if near_field_analysis || save_end
        if derivatives
            near_field_properties_and_derivatives!(properties, dproperties,
                surfaces, wakes, reference, freestream, additional_velocity,
                Γ, dΓ, Γdot; wake_shedding_locations, Vh, Vv, symmetric, surface_id,
                wake_finite_core, iwake, trailing_vortices, xhat)
        else
            near_field_properties!(properties, surfaces, wakes, reference,
                freestream, additional_velocity, Γ, Γdot; wake_shedding_locations,
                Vh, Vv, symmetric, surface_id, wake_finite_core, iwake,
                trailing_vortices, xhat)
        end

        if save_end
            time_history[isave] = tvec[end]
            surface_history[isave] = [copy(system.surfaces[isurf]) for isurf = 1:nsurf]
            property_history[isave] = [copy(system.properties[isurf]) for isurf = 1:nsurf]
            wake_history[isave] = [system.wakes[isurf][1:iwake[isurf], :] for isurf = 1:nsurf]

            # increment index for saving
            isave += 1
        end
    end

    # return the modified system and time history
    return system, time_history, surface_history, property_history, wake_history
end
