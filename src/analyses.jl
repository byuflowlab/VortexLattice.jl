"""
    steady_analysis(surface, reference, freestream; kwargs...)

Perform a steady vortex lattice method analysis.  Return an object of type
`System` containing the system state.

# Arguments
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see `Reference`)
 - `freestream`: Freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image of the panels
    in `surface`, should be used when calculating induced velocities
 - `trailing_vortices`: Flag to enable/disable trailing vortices, defaults to `true`
 - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]
 - `wake`: Matrix of wake panels of shape (nw, ns) where `nw` is the number of
    chordwise wake panels and `ns` is the number of spanwise panels, defaults to
    no wake panels
 - `nwake`: Number of chordwise wake panels to use from `wake`, defaults to all
    specified wake panels
 - `surface_id`: Surface ID.  The finite core model is disabled when calculating
    the influence of surfaces/wakes that share the same ID.  Additionally, if a
    surface/wake's ID is negative, the finite core model will always be enabled,
    even when calculating the influence of the surface/wake on itself. By default
    each surface has its own ID.
 - `wake_id`: Wake ID.  The finite core model is disabled when calculating
    the influence of surfaces/wakes that share the same ID. Additionally, if a
    surface/wake's ID is negative, the finite core model will always be enabled,
    even when calculating the influence of the surface/wake on itself. By default
    the finite core model for the wakes is always enabled.
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix has already been calculated.  Re-using the same AIC matrix
    will reduce calculation times when the underlying geometry has not changed.
    Defaults to `true`.  Note that this argument is only valid for the pre-allocated
    version of this function.
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces. Defaults to `true`
 - `derivatives`: Flag indicating whether the derivatives with respect
    to the freestream variables should be calculated. Defaults to `true`
"""
function steady_analysis(surface::AbstractMatrix, reference, freestream; kwargs...)

    system = System(surface)

    return steady_analysis!(system, surface, reference, freestream; kwargs...,
        calculate_influence_matrix = true)
end

"""
    steady_analysis(surfaces, reference, freestream; kwargs...)

Perform a steady vortex lattice method analysis.  Return an object of type
`System` containing the system state.

# Arguments
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: Reference parameters (see `Reference`)
 - `freestream`: Freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: (required) Flags indicating whether a mirror image (across the X-Z plane)
    should be used when calculating induced velocities
 - `wakes`: Vector of wakes corresponding to each surface, represented by matrices
    of panels of shape (nw, ns) where `nw` is the number of chordwise wake panels
    and `ns` is the number of spanwise panels
 - `nwake`: Number of chordwise wake panels to use from each wake in `wakes`, defaults
    to all specified wake panels
 - `trailing_vortices`: Flags to enable/disable trailing vortices, defaults to `true`
    for each surface
 - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
    Additionally, if a surface/wake's ID is negative the finite core model will
    always be enabled, even when calculating the influence of the surface/wake
    on itself. By default each surface has its own ID.
 - `wake_id`: Wake ID for each wake.  The finite core model is disabled when
    calculating the influence of surfaces/wakes that share the same ID.
    Additionally, if a surface/wake's ID is negative, the finite core model will
    always be enabled, even when calculating the influence of the surface/wake
    on itself. By default the finite core model for the wakes is always enabled.
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix has already been calculated.  Re-using the same AIC matrix
    will reduce calculation times when the underlying geometry has not changed.
    Defaults to `true`.  Note that this argument is only valid for the pre-allocated
    version of this function.
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces. Defaults to `true`
 - `derivatives`: Flag indicating whether the derivatives with respect
    to the freestream variables should be calculated. Defaults to `true`
"""
function steady_analysis(surfaces::AbstractVector{<:AbstractMatrix}, reference,
    freestream; kwargs...)

    system = System(surfaces)

    return steady_analysis!(system, surfaces, reference, freestream; kwargs...,
        calculate_influence_matrix = true)
end

"""
    steady_analysis!(system, surface, reference, freestream; kwargs...)

Pre-allocated version of `steady_analysis`.
"""
function steady_analysis!(system, surface::AbstractMatrix, ref, fs;
    symmetric,
    wake = Matrix{Wake{Float64}}(undef, 0, size(surface,2)),
    nwake = size(wake, 1),
    trailing_vortices = true,
    xhat = SVector(1, 0, 0),
    surface_id = 1,
    wake_id = -1,
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # unpack system variables
    AIC = system.AIC
    b = system.b
    Γ = system.gamma
    db = system.db
    dΓ = system.dgamma

    # see if wake panels are being used
    wake_panels = nwake > 0

    # store new wake panels in system
    system.wakes[1] = wake

    # calculate AIC matrix
    if calculate_influence_matrix
        influence_coefficients!(AIC, surface;
            symmetric = symmetric,
            trailing_vortices = trailing_vortices && !wake_panels,
            xhat = xhat,
            surface_id = surface_id)
    end

    # calculate RHS
    if derivatives
        normal_velocity_derivatives!(b, db, surface, ref, fs)
    else
        normal_velocity!(b, surface, ref, fs)
    end

    # add wake's contribution to RHS
    if wake_panels
        add_wake_normal_velocity!(b, surface, wake;
            symmetric = symmetric,
            nwake = nwake,
            trailing_vortices = trailing_vortices,
            xhat = xhat,
            surface_id = surface_id,
            wake_id = wake_id)
    end

    # solve for the circulation distribution
    if derivatives
        circulation_derivatives!(Γ, dΓ, AIC, b, db)
    else
        circulation!(Γ, AIC, b)
    end

    if near_field_analysis
        # perform a near field analysis to obtain panel properties
        if derivatives
            near_field_forces_derivatives!(system, surface, ref, fs;
                symmetric = symmetric,
                nwake = nwake,
                trailing_vortices = trailing_vortices,
                xhat = xhat,
                surface_id = surface_id,
                wake_id = wake_id)
        else
            near_field_forces!(system, surface, ref, fs;
                symmetric = symmetric,
                nwake = nwake,
                trailing_vortices = trailing_vortices,
                xhat = xhat,
                surface_id = surface_id,
                wake_id = wake_id)
        end
    end

    # return the modified system
    return system
end

"""
    steady_analysis!(system, surfaces, reference, freestream; kwargs...)

Pre-allocated version of `steady_analysis`.
"""
function steady_analysis!(system, surfaces::AbstractVector{<:AbstractMatrix},
    ref, fs;
    symmetric,
    wakes = [Matrix{Wake{Float64}}(undef, 0, size(surfaces[i],2)) for i = 1:length(surfaces)],
    nwake = size.(wakes, 1),
    trailing_vortices = fill(true, length(surfaces)),
    xhat = SVector(1, 0, 0),
    surface_id = 1:length(surfaces),
    wake_id = -1:-1:-length(surfaces),
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # unpack system variables
    AIC = system.AIC
    b = system.b
    Γ = system.gamma
    db = system.db
    dΓ = system.dgamma

    # see if wake panels are being used
    wake_panels = nwake .> 0

    # store new wake panels in system
    for i = 1:length(wakes)
        system.wakes[i] = wakes[i]
    end

    # calculate AIC matrix
    if calculate_influence_matrix
        influence_coefficients!(AIC, surfaces;
            symmetric = symmetric,
            trailing_vortices = trailing_vortices .& .!wake_panels,
            xhat = xhat,
            surface_id = surface_id)
    end

    # calculate RHS
    if derivatives
        normal_velocity_derivatives!(b, db, surfaces, ref, fs)
    else
        normal_velocity!(b, surfaces, ref, fs)
    end

    # add wake's contribution to RHS
    if any(wake_panels)
        add_wake_normal_velocity!(b, surfaces, wakes;
            symmetric = symmetric,
            nwake = nwake,
            trailing_vortices = trailing_vortices,
            xhat = xhat,
            surface_id = surface_id,
            wake_id = wake_id)
    end

    # solve for the circulation distribution
    if derivatives
        circulation_derivatives!(Γ, dΓ, AIC, b, db)
    else
        circulation!(Γ, AIC, b)
    end

    if near_field_analysis
        # perform a near field analysis to obtain panel properties
        if derivatives
            near_field_forces_derivatives!(system, surfaces, ref, fs;
                symmetric = symmetric,
                nwake = nwake,
                trailing_vortices = trailing_vortices,
                xhat = xhat,
                surface_id = surface_id,
                wake_id = wake_id)
        else
            near_field_forces!(system, surfaces, ref, fs;
                symmetric = symmetric,
                nwake = nwake,
                trailing_vortices = trailing_vortices,
                xhat = xhat,
                surface_id = surface_id,
                wake_id = wake_id)
        end
    end

    # return the modified system
    return system
end

"""
    unsteady_analysis(surface, reference, freestream, time; kwargs...)

Perform a unsteady vortex lattice method analysis.  Return an object of type
`System` containing the system state.

# Arguments
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see `Reference`)
 - `freestream`: Freestream parameters (see `Freestream`)
 - `time`: Time vector

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image of the panels in `surface`,
    should be used when calculating induced velocities, defaults to `false`
 - `trailing_vortices`: Flag to enable/disable trailing vortices from the wake,
    defaults to `false`
 - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]
 - `nwake`: Maximum number of wake panels in the chordwise direction.  Defaults
    to `length(time)-1`.
 - `initial_wake`: Initial wake panels of the system, defaults to no wake panels
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
    Additionally, if a surface/wake's ID is negative the finite core model will
    always be enabled, even when calculating the influence of the surface/wake
    on itself. By default each surface has its own ID.
 - `wake_id`: Wake ID for each wake.  The finite core model is disabled when
    calculating the influence of surfaces/wakes that share the same ID.
    Additionally, if a surface/wake's ID is negative, the finite core model will
    always be enabled, even when calculating the influence of the surface/wake
    on itself. By default the finite core model for the wakes is always enabled.
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix has already been calculated.  Re-using the same AIC matrix
    will reduce calculation times when the underlying geometry has not changed.
    Defaults to `true`.  Note that this argument is only valid for the pre-allocated
    version of this function.
 - `save`: Time indices at which to save the time history, defaults to `1:length(time)`
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces for the last
    time step. Defaults to `true`.
 - `derivatives`: Flag indicating whether the derivatives with respect
    to the freestream variables should be calculated for the last time step.
    Defaults to `true`.
"""
function unsteady_analysis(surface::AbstractMatrix, reference, freestream, time;
    nwake = length(time)-1, kwargs...)

    system = System(surface, nwake=nwake)

    return unsteady_analysis!(system, surface, reference, freestream, time;
        kwargs..., nwake = nwake, calculate_influence_matrix = true)
end

"""
    unsteady_analysis(surfaces, reference, freestream, time; kwargs...)

Perform a unsteady vortex lattice method analysis to find the unsteady circulation
and wake shape of a group of vortex lattice panels.

# Arguments
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: Reference parameters (see `Reference`)
 - `freestream`: Freestream parameters (see `Freestream`)
 - `time`: time vector

# Keyword Arguments
 - `nwake`: Maximum number of wake panels in the chordwise direction.  Defaults
    to `length(time)-1`.
 - `initial_wake`: Initial wake panels of the system
 - `symmetric`: Flags indicating whether a mirror image (across the X-Z plane) should
    be used when calculating induced velocities, defaults to `false` for each surface
 - `trailing_vortices`: Flags to enable/disable trailing vortices, defaults to `false`
    for each surface
 - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
    Additionally, if a surface/wake's ID is negative the finite core model will
    always be enabled, even when calculating the influence of the surface/wake
    on itself. By default each surface has its own ID.
 - `wake_id`: Wake ID for each wake.  The finite core model is disabled when
    calculating the influence of surfaces/wakes that share the same ID.
    Additionally, if a surface/wake's ID is negative, the finite core model will
    always be enabled, even when calculating the influence of the surface/wake
    on itself. By default the finite core model for the wakes is always enabled
    for wakes.
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix has already been calculated.  Re-using the same AIC matrix
    will reduce calculation times when the underlying geometry has not changed.
    Defaults to `true`.  Note that this argument is only valid for the pre-allocated
    version of this function.
 - `save`: Time indices at which to save the time history, defaults to `1:length(time)`
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces for the last
    time step. Defaults to `true`.
 - `derivatives`: Flag indicating whether the derivatives with respect
    to the freestream variables should be calculated for the last time step.
    Defaults to `true`.
"""
function unsteady_analysis(surfaces::AbstractVector{<:AbstractMatrix}, reference,
    freestream, time; nwake = fill(length(time) - 1, length(surfaces)), kwargs...)

    system = System(surfaces; nwake = nwake)

    return unsteady_analysis!(system, surfaces, reference, freestream, time;
        kwargs..., nwake = nwake, calculate_influence_matrix = true)
end

"""
    unsteady_analysis!(system, surface, reference, freestream, time; kwargs...)

Pre-allocated version of `unsteady_analysis`.
"""
function unsteady_analysis!(system, surface::AbstractMatrix, ref, fs, time;
    symmetric,
    trailing_vortices = false,
    xhat = SVector(1, 0, 0),
    initial_wake = Matrix{Wake{Float64}}(undef, 0, size(surface, 2)),
    nwake = length(time) - 1,
    surface_id = 1,
    wake_id = -1,
    calculate_influence_matrix = true,
    save = 1:length(time),
    near_field_analysis = true,
    derivatives = true)

    TF = eltype(system)

    # check if existing wake panel storage is sufficient, replace if necessary
    if size(system.wakes[1], 2) < nwake
        # construct new wake panel storage
        wake = Matrix{Wake{TF}}(undef, nwake, size(surface, 2))
        # replace old wake panels
        system.wakes[1] = wake
    end

    # copy initial wake panels to pre-allocated storage
    for I in CartesianIndices(initial_wake)
        wake[I] = initial_wake[I]
    end

    # check for trailing edge points that are repeated
    repeated_points = repeated_trailing_edge_points([surface])

    # unpack pre-allocated storage
    wake = system.wakes[1]
    wake_velocities = system.wake_velocities[1]
    Γ = system.gamma

    # current number of wake panels
    iwake = min(size(initial_wake, 1), nwake)

    # initialize solution history for each time step
    isave = 1
    surface_history = Vector{Matrix{PanelProperties{TF}}}(undef, length(save))
    wake_history = Vector{Matrix{Wake{TF}}}(undef, length(save))

    # # loop through all time steps
    for it = 1:length(time)

        # perform a steady analysis
        steady_analysis!(system, surface, ref, fs;
            symmetric = symmetric,
            surface_id = surface_id,
            wake_id = wake_id,
            wake = wake,
            nwake = iwake,
            trailing_vortices = trailing_vortices,
            xhat = xhat,
            calculate_influence_matrix = it == 1 && calculate_influence_matrix,
            near_field_analysis = false,
            derivatives = false)

        # save the panel history
        if it in save
            surface_history[isave] = deepcopy(system.panels[1])
            wake_history[isave] = wake[1:iwake, :]
            isave += 1
        end

        if it < length(time)

            # calculate wake corner point velocities
            VortexLattice.get_wake_velocities!(wake_velocities, surface, wake, ref, fs, Γ;
                symmetric = symmetric,
                surface_id = surface_id,
                wake_id = wake_id,
                trailing_vortices = trailing_vortices,
                xhat = xhat,
                nwake = iwake,
                repeated_points = repeated_points)

            # increment wake panel counter
            if iwake < nwake
                iwake += 1
            end

            # translate the wake
            dt = time[it+1] - time[it]
            VortexLattice.shed_wake!(wake, wake_velocities, dt, surface, Γ; nwake=iwake)
        end

    end

    # return the modified system and the time history
    return system, surface_history, wake_history
end

"""
    unsteady_analysis!(system, surfaces, reference, freestream, time; kwargs...)

Pre-allocated version of `unsteady_analysis`.
"""
function unsteady_analysis!(system, surfaces::AbstractVector{<:AbstractMatrix}, ref, fs, time;
    symmetric,
    trailing_vortices = fill(false, length(surfaces)),
    xhat = SVector(1, 0, 0),
    surface_id = 1:length(surfaces),
    wake_id = -1:-1:-length(surfaces),
    initial_wakes = [Matrix{Wake{Float64}}(undef, 0, size(surfaces[i], 2)) for i = 1:length(surfaces)],
    nwake = fill(length(time) - 1, length(surfaces)),
    calculate_influence_matrix = true,
    save = 1:length(time),
    near_field_analysis = true,
    derivatives = true)

    TF = eltype(system)
    nsurf = length(surfaces)

    # check if existing wake panel storage is sufficient, replace if necessary
    for isurf = 1:nsurf
        if size(system.wakes[isurf], 2) < nwake[isurf]
            # construct new wake panel storage
            system.wakes[isurf] = Matrix{Wake{TF}}(undef, nwake[isurf], size(surfaces[isurf], 2))
        end
    end

    # copy initial wake panels to pre-allocated storage
    for isurf = 1:nsurf
        for I in CartesianIndices(initial_wakes[isurf])
            system.wakes[isurf][I] = initial_wakes[isurf][I]
        end
    end

    # check for trailing edge points that are repeated
    repeated_points = repeated_trailing_edge_points(surfaces)

    # unpack pre-allocated storage
    wakes = system.wakes
    wake_velocities = system.wake_velocities
    Γ = system.gamma

    # current number of wake panels
    iwake = [min(size(initial_wakes[isurf], 1), nwake[isurf]) for isurf = 1:nsurf]

    # initialize solution history for each time step
    isave = 1
    surface_history = Vector{Vector{Matrix{PanelProperties{TF}}}}(undef, length(save))
    wake_history = Vector{Vector{Matrix{Wake{TF}}}}(undef, length(save))

    # # loop through all time steps
    for it = 1:length(time)

        # perform a steady analysis
        steady_analysis!(system, surfaces, ref, fs;
            symmetric = symmetric,
            surface_id = surface_id,
            wake_id = wake_id,
            wakes = wakes,
            nwake = iwake,
            trailing_vortices = trailing_vortices,
            xhat = xhat,
            calculate_influence_matrix = it == 1 && calculate_influence_matrix,
            near_field_analysis = false,
            derivatives = false)

        # save the panel history
        if it in save
            surface_history[isave] = deepcopy(system.panels)
            wake_history[isave] = [wakes[isurf][1:iwake[isurf], :] for isurf = 1:nsurf]
            isave += 1
        end

        if it < length(time)

            # calculate wake corner point velocities
            VortexLattice.get_wake_velocities!(wake_velocities, surfaces, wakes, ref, fs, Γ;
                symmetric = symmetric,
                surface_id = surface_id,
                wake_id = wake_id,
                trailing_vortices = trailing_vortices,
                xhat = xhat,
                nwake = iwake,
                repeated_points = repeated_points)

            # increment wake panel counter
            for isurf = 1:nsurf
                if iwake[isurf] < nwake[isurf]
                    iwake[isurf] += 1
                end
            end

            # shed an additional wake panel from each surface
            dt = time[it+1] - time[it]
            VortexLattice.shed_wake!(wakes, wake_velocities, dt, surfaces, Γ; nwake=iwake)
        end

    end

    # return the modified system and time history
    return system, surface_history, wake_history
end
