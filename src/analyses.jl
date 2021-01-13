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
 - `wake_id`: WakePanel ID.  The finite core model is disabled when calculating
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
 - `wake_id`: WakePanel ID for each wake.  The finite core model is disabled when
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
    wake = Matrix{WakePanel{Float64}}(undef, 0, size(surface,2)),
    nwake = size(wake, 1),
    trailing_vortices = true,
    xhat = SVector(1, 0, 0),
    surface_id = 1,
    wake_id = -1,
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # see if wake panels are being used
    wake_panels = nwake > 0

    # store new wake panels in system
    system.wakes[1] = wake

    # unpack system variables
    AIC = system.AIC
    w = system.w
    Γ = system.Γ
    dw = system.dw
    dΓ = system.dΓ
    dΓdt = system.dΓdt

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
        normal_velocity_derivatives!(w, dw, surface, ref, fs)
    else
        normal_velocity!(w, surface, ref, fs)
    end

    # add wake's contribution to RHS
    if wake_panels
        subtract_wake_normal_velocity!(w, surface, wake;
            symmetric = symmetric,
            nwake = nwake,
            trailing_vortices = trailing_vortices,
            xhat = xhat,
            surface_id = surface_id,
            wake_id = wake_id)
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
    wakes = [Matrix{WakePanel{Float64}}(undef, 0, size(surfaces[i],2)) for i = 1:length(surfaces)],
    nwake = size.(wakes, 1),
    trailing_vortices = fill(true, length(surfaces)),
    xhat = SVector(1, 0, 0),
    surface_id = 1:length(surfaces),
    wake_id = -1:-1:-length(surfaces),
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # see if wake panels are being used
    wake_panels = nwake .> 0

    # store new wake panels in system
    for i = 1:length(wakes)
        system.wakes[i] = wakes[i]
    end

    # unpack system variables
    AIC = system.AIC
    w = system.w
    Γ = system.Γ
    dw = system.dw
    dΓ = system.dΓ

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
        normal_velocity_derivatives!(w, dw, surfaces, ref, fs)
    else
        normal_velocity!(w, surfaces, ref, fs)
    end

    # add wake's contribution to RHS
    if any(wake_panels)
        subtract_wake_normal_velocity!(w, surfaces, wakes;
            symmetric = symmetric,
            nwake = nwake,
            trailing_vortices = trailing_vortices,
            xhat = xhat,
            surface_id = surface_id,
            wake_id = wake_id)
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
    unsteady_analysis(surface, reference, freestream, dt; kwargs...)

Perform a unsteady vortex lattice method analysis.  Return `system, surface_history, wake_history`
where `system` is the modified system state, `surface_history` is the time
history of the surface, and `wake_history` is the time history of the wake.

# Arguments
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see `Reference`)
 - `freestream`: Freestream parameters for each time step (see `Freestream`)
 - `dt`: Non-dimensional time step size for each time step (defined as `Vinf*Δt/cref`)

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image of the panels in `surface`,
    should be used when calculating induced velocities, defaults to `false`
 - `initial_wake`: Initial wake panels of the system, defaults to no wake panels
 - `nwake`: Maximum number of wake panels in the chordwise direction.  Defaults
    to `length(dt)`.
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
    Additionally, if a surface/wake's ID is negative the finite core model will
    always be enabled, even when calculating the influence of the surface/wake
    on itself. By default each surface has its own ID.
 - `wake_id`: WakePanel ID for each wake.  The finite core model is disabled when
    calculating the influence of surfaces/wakes that share the same ID.
    Additionally, if a surface/wake's ID is negative, the finite core model will
    always be enabled, even when calculating the influence of the surface/wake
    on itself. By default the finite core model for the wakes is always enabled.
 - `save`: Time steps at which to save the time history, defaults to `1:1+length(dt)`
    where index `1` corresponds to solution for the initial conditions.
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix has already been calculated.  Re-using the same AIC matrix
    will reduce calculation times when the underlying geometry has not changed.
    Defaults to `true`.  Note that this argument is only valid for the pre-allocated
    version of this function.
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces for the last
    time step. Defaults to `true`.
 - `derivatives`: Flag indicating whether the derivatives with respect
    to the freestream variables should be calculated for the last time step.
    Defaults to `true`.
"""
function unsteady_analysis(surface::AbstractMatrix, reference, freestream, dt;
    nwake = length(dt), kwargs...)

    system = System(surface, nwake=nwake)

    return unsteady_analysis!(system, surface, reference, freestream, dt;
        kwargs..., nwake = nwake, calculate_influence_matrix = true)
end

"""
    unsteady_analysis(surfaces, reference, freestream, dt; kwargs...)

Perform a unsteady vortex lattice method analysis.  Return `system, surface_history, wake_history`
where `system` is the modified system state, `surface_history` is the time history
of the surfaces, and `wake_history` is the time history of the wakes.

# Arguments
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: Reference parameters (see `Reference`)
 - `freestream`: Freestream parameters for each time step (see `Freestream`)
 - `Vinf`: Freestream velocity magnitude for each time step
 - `dt`: Non-dimensional time step size for each time step (defined as `Vinf*Δt/cref`)

# Keyword Arguments
 - `symmetric`: Flags indicating whether a mirror image (across the X-Z plane) should
    be used when calculating induced velocities, defaults to `false` for each surface
 - `nwake`: Maximum number of wake panels in the chordwise direction.  Defaults
    to `length(dt)`.
 - `initial_wake`: Initial wake panels of the system
 - `circulation_derivatives`:  Vector containing time derivatives of the surface
    panel circulation strengths at the initial time step. Defaults to the current
    values stored in `system`.
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
    Additionally, if a surface/wake's ID is negative the finite core model will
    always be enabled, even when calculating the influence of the surface/wake
    on itself. By default each surface has its own ID.
 - `wake_id`: WakePanel ID for each wake.  The finite core model is disabled when
    calculating the influence of surfaces/wakes that share the same ID.
    Additionally, if a surface/wake's ID is negative, the finite core model will
    always be enabled, even when calculating the influence of the surface/wake
    on itself. By default the finite core model for the wakes is always enabled
    for wakes.
 - `save`: Time indices at which to save the time history, defaults to `1:1+length(dt)`
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix has already been calculated.  Re-using the same AIC matrix
    will reduce calculation times when the underlying geometry has not changed.
    Defaults to `true`.  Note that this argument is only valid for the pre-allocated
    version of this function.
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces for the last
    time step. Defaults to `true`.
 - `derivatives`: Flag indicating whether the derivatives with respect
    to the freestream variables should be calculated for the last time step.
    Defaults to `true`.
"""
function unsteady_analysis(surfaces::AbstractVector{<:AbstractMatrix}, reference,
    freestream, dt; nwake = fill(length(dt), length(surfaces)), kwargs...)

    system = System(surfaces; nwake = nwake)

    return unsteady_analysis!(system, surfaces, reference, freestream, dt;
        kwargs..., nwake = nwake, calculate_influence_matrix = true)
end

"""
    unsteady_analysis!(system, surface, reference, freestream, dt; kwargs...)

Pre-allocated version of `unsteady_analysis`.
"""
function unsteady_analysis!(system, surface::AbstractMatrix, ref, fs, dt;
    symmetric,
    initial_circulation = zeros(length(surface)),
    initial_wake = Matrix{WakePanel{Float64}}(undef, 0, size(surface, 2)),
    nwake = length(dt),
    surface_id = 1,
    wake_id = -1,
    save = 1:1+length(dt),
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    TF = eltype(system)

    N = length(surface)
    trailing_vortices = false
    xhat = SVector(1, 0, 0)

    # check if existing wake panel storage is sufficient, replace if necessary
    if size(system.wakes[1], 2) < nwake
        # construct new wake panel storage
        wake = Matrix{WakePanel{TF}}(undef, nwake, size(surface, 2))
        # replace old wake panels
        system.wakes[1] = wake
    end

    # copy initial wake panels to pre-allocated storage
    for I in CartesianIndices(initial_wake)
        system.wakes[1][I] = initial_wake[I]
    end

    # unpack pre-allocated storage
    AIC = system.AIC
    w = system.w
    Γ = system.Γ
    wake = system.wakes[1]
    wake_velocities = system.V[1]
    dw = system.dw
    dΓ = system.dΓ
    dΓdt = system.dΓdt

    # check for trailing edge points that are repeated
    repeated_points = repeated_trailing_edge_points(surface)

    # current number of wake panels
    iwake = min(size(initial_wake, 1), nwake)

    # initialize solution history for each time step
    isave = 1
    surface_history = Vector{Matrix{PanelProperties{TF}}}(undef, length(save))
    wake_history = Vector{Matrix{WakePanel{TF}}}(undef, length(save))

    # loop through all time steps (including t = 0)
    for it = 1 : 1 + length(dt)

        first_step = it == 1
        last_step = it == 1 + length(dt)
        wake_panels = iwake > 0

        # calculate AIC matrix
        if first_step && calculate_influence_matrix
            influence_coefficients!(AIC, surface;
                symmetric = symmetric,
                trailing_vortices = trailing_vortices && !wake_panels,
                xhat = xhat,
                surface_id = surface_id)
        end

        # calculate RHS
        if last_step && derivatives
            normal_velocity_derivatives!(w, dw, surface, ref, fs)
        else
            normal_velocity!(w, surface, ref, fs)
        end

        # add wake's contribution to RHS
        if wake_panels
            subtract_wake_normal_velocity!(w, surface, wake;
                symmetric = symmetric,
                nwake = iwake,
                trailing_vortices = trailing_vortices,
                xhat = xhat,
                surface_id = surface_id,
                wake_id = wake_id)
        end

        # begin constructing dΓdt using finite difference `dΓdt = (Γ - Γp)/dt`
        if first_step
            for i = 1:N
                dΓdt[i] = -initial_circulation[i]
            end
        else
            dΓdt .= -Γ
        end

        # solve for the circulation distribution
        if last_step && derivatives
            circulation_derivatives!(Γ, dΓ, AIC, w, dw)
        else
            circulation!(Γ, AIC, w)
        end

        # finish constructing dΓdt using finite difference `dΓdt = (Γ - Γp)/dt`
        if first_step
            dΓdt .+= Γ
            dΓdt ./= 1e-6 # arbitrarily small value (change is instantaneous)
        else
            dΓdt .+= Γ
            dΓdt ./= dt[it-1]
        end

        if it in save || (last_step && near_field_analysis)
            # perform a near field analysis to obtain panel properties
            if last_step && derivatives
                near_field_forces_derivatives!(system, surface, ref, fs;
                    symmetric = symmetric,
                    nwake = iwake,
                    trailing_vortices = trailing_vortices,
                    xhat = xhat,
                    surface_id = surface_id,
                    wake_id = wake_id)
            else
                near_field_forces!(system, surface, ref, fs;
                    symmetric = symmetric,
                    nwake = iwake,
                    trailing_vortices = trailing_vortices,
                    xhat = xhat,
                    surface_id = surface_id,
                    wake_id = wake_id)
            end
        end

        # save the panel history
        if it in save
            surface_history[isave] = deepcopy(system.panels[1])
            wake_history[isave] = wake[1:iwake, :]
            isave += 1
        end

        # move the wake
        if it < 1 + length(dt)

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
            dx = dt[it]*ref.c
            VortexLattice.shed_wake!(wake, wake_velocities, dx, surface, Γ; nwake=iwake)
        end

    end

    # return the modified system and the time history
    return system, surface_history, wake_history
end

"""
    unsteady_analysis!(system, surfaces, reference, freestream, dt; kwargs...)

Pre-allocated version of `unsteady_analysis`.
"""
function unsteady_analysis!(system, surfaces::AbstractVector{<:AbstractMatrix}, ref, fs, dt;
    symmetric,
    initial_wakes = [Matrix{WakePanel{Float64}}(undef, 0, size(surfaces[i], 2)) for i = 1:length(surfaces)],
    nwake = fill(length(dt), length(surfaces)),
    surface_id = 1:length(surfaces),
    wake_id = -1:-1:-length(surfaces),
    save = 1:1+length(dt),
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    TF = eltype(system)
    nsurf = length(surfaces)

    trailing_vortices = fill(false, length(surfaces))
    xhat = SVector(1, 0, 0)

    # check if existing wake panel storage is sufficient, replace if necessary
    for isurf = 1:nsurf
        if size(system.wakes[isurf], 2) < nwake[isurf]
            # construct new wake panel storage
            system.wakes[isurf] = Matrix{WakePanel{TF}}(undef, nwake[isurf], size(surfaces[isurf], 2))
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
    AIC = system.AIC
    w = system.w
    Γ = system.Γ
    wakes = system.wakes
    wake_velocities = system.V
    dw = system.dw
    dΓ = system.dΓ
    dΓdt = system.dΓdt

    # current number of wake panels
    iwake = [min(size(initial_wakes[isurf], 1), nwake[isurf]) for isurf = 1:nsurf]

    # initialize solution history for each time step
    isave = 1
    surface_history = Vector{Vector{Matrix{PanelProperties{TF}}}}(undef, length(save))
    wake_history = Vector{Vector{Matrix{WakePanel{TF}}}}(undef, length(save))

    # # loop through all time steps
    for it = 1:1+length(dt)

        first_step = it == 1
        last_step = it == 1 + length(dt)
        wake_panels = iwake .> 0

        # calculate AIC matrix
        if first_step && calculate_influence_matrix
            influence_coefficients!(AIC, surfaces;
                symmetric = symmetric,
                trailing_vortices = trailing_vortices .& .!wake_panels,
                xhat = xhat,
                surface_id = surface_id)
        end

        # calculate RHS
        if last_step && derivatives
            normal_velocity_derivatives!(w, dw, surfaces, ref, fs)
        else
            normal_velocity!(w, surfaces, ref, fs)
        end

        # add wake's contribution to RHS
        if any(wake_panels)
            subtract_wake_normal_velocity!(w, surfaces, wakes;
                symmetric = symmetric,
                nwake = iwake,
                trailing_vortices = trailing_vortices,
                xhat = xhat,
                surface_id = surface_id,
                wake_id = wake_id)
        end

        # begin constructing dΓdt using finite difference `dΓdt = (Γ - Γp)/dt`
        if !first_step # dΓdt at first step is given
            dΓdt .= -Γ
        end

        # solve for the circulation distribution
        if last_step && derivatives
            circulation_derivatives!(Γ, dΓ, AIC, w, dw)
        else
            circulation!(Γ, AIC, w)
        end

        # finish constructing dΓdt using finite difference `dΓdt = (Γ - Γp)/dt`
        if !first_step # dΓdt at first step is given
            dΓdt .+= Γ
            dΓdt ./= dt[it-1]
        end

        if it in save || (last_step && near_field_analysis)
            # perform a near field analysis to obtain panel properties
            if last_step && derivatives
                near_field_forces_derivatives!(system, surfaces, ref, fs;
                    symmetric = symmetric,
                    nwake = iwake,
                    trailing_vortices = trailing_vortices,
                    xhat = xhat,
                    surface_id = surface_id,
                    wake_id = wake_id,
                    Gamma_t = dΓdt)
            else
                near_field_forces!(system, surfaces, ref, fs;
                    symmetric = symmetric,
                    nwake = iwake,
                    trailing_vortices = trailing_vortices,
                    xhat = xhat,
                    surface_id = surface_id,
                    wake_id = wake_id,
                    Gamma_t = dΓdt)
            end
        end

        # save the panel history
        if it in save
            surface_history[isave] = deepcopy(system.panels)
            wake_history[isave] = [wakes[isurf][1:iwake[isurf], :] for isurf = 1:nsurf]
            isave += 1
        end

        if it < 1 + length(dt)

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
            dx = dt[it]*ref.c
            VortexLattice.shed_wake!(wakes, wake_velocities, dx, surfaces, Γ; nwake=iwake)
        end

    end

    # return the modified system and time history
    return system, surface_history, wake_history
end
