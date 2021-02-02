"""
    steady_analysis(surface, reference, freestream; kwargs...)

Perform a steady vortex lattice method analysis.  Return an object of type
[`System`](@ref) containing the system state.

# Arguments
 - `surface`: Matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image (across the
    X-Z plane) of the panels in `surface` should be used when calculating induced velocities
 - `wake`: Matrix of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns)
    where `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels, defaults to no wake panels
 - `nwake`: Number of chordwise wake panels to use from `wake`, defaults to all
    provided wake panels
 - `wake_finite_core`: Flag indicating whether the finite core model should be
    enabled when calculating a wake's influence on itself and its corresponding
    surface. Defaults to `true`
 - `trailing_vortices`: Flag to enable/disable trailing vortices, defaults to `true`
 - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]
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
function steady_analysis(surface::AbstractMatrix, reference, freestream; kwargs...)

    system = System(surface)

    return steady_analysis!(system, surface, reference, freestream; kwargs...,
        calculate_influence_matrix = true)
end

"""
    steady_analysis(surfaces, reference, freestream; kwargs...)

Perform a steady vortex lattice method analysis.  Return an object of type
[`System`](@ref) containing the system state.

# Arguments
 - `surfaces`: Vector of surfaces, represented by matrices of surface panels
    (see [`SurfacePanel`](@ref) of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`]@ref)

# Keyword Arguments
 - `symmetric`: (required) Flag for each surface indicating whether a mirror
    image across the X-Z plane should be used when calculating induced velocities
 - `wakes`: Vector of wakes corresponding to each surface, represented by matrices
    of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where `nw` is the
    number of chordwise wake panels and `ns` is the number of spanwise panels.
    Defaults to no wake panels.
 - `nwake`: Number of chordwise wake panels to use from each wake in `wakes`,
    defaults to all provided wake panels
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.  Defaults to `true` for each surface.
 - `trailing_vortices`: Flags to enable/disable trailing vortices, defaults to
    `true` for each surface
 - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]
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
    wake_finite_core = true,
    trailing_vortices = true,
    xhat = SVector(1, 0, 0),
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

    # calculate AIC matrix
    if calculate_influence_matrix
        influence_coefficients!(AIC, surface;
            symmetric = symmetric,
            trailing_vortices = trailing_vortices && !wake_panels,
            xhat = xhat)
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
            near_field_forces_derivatives!(system, surface, ref, fs;
                symmetric = symmetric,
                nwake = nwake,
                wake_finite_core = wake_finite_core,
                trailing_vortices = trailing_vortices,
                xhat = xhat)
        else
            near_field_forces!(system, surface, ref, fs;
                symmetric = symmetric,
                nwake = nwake,
                wake_finite_core = wake_finite_core,
                trailing_vortices = trailing_vortices,
                xhat = xhat)
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
    surface_id = 1:length(surfaces),
    wake_finite_core = fill(true, length(surfaces)),
    trailing_vortices = fill(true, length(surfaces)),
    xhat = SVector(1, 0, 0),
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # see if wake panels are being used
    wake_panels = nwake .> 0

    # store new wake panels in the system
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
            surface_id = surface_id,
            trailing_vortices = trailing_vortices .& .!wake_panels,
            xhat = xhat)
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
            surface_id = surface_id,
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
            near_field_forces_derivatives!(system, surfaces, ref, fs;
                symmetric = symmetric,
                nwake = nwake,
                surface_id = surface_id,
                wake_finite_core = wake_finite_core,
                trailing_vortices = trailing_vortices,
                xhat = xhat)
        else
            near_field_forces!(system, surfaces, ref, fs;
                symmetric = symmetric,
                nwake = nwake,
                surface_id = surface_id,
                wake_finite_core = wake_finite_core,
                trailing_vortices = trailing_vortices,
                xhat = xhat)
        end
    end

    # return the modified system
    return system
end

"""
    prescribed_motion(dt; kwargs...)

Return (`Vinf`, `freestream`) given an aircraft/wing's trajectory where `Vinf`
is the freestream velocity magnitude for each time step and `freestream` is an
object of type [`Freestream`](@ref).

# Arguments:
 - `dt`: Time step vector (seconds)

# Keyword Arguments:
 - `Xdot = zeros(length(dt))`: Global frame x-velocity for each time step
 - `Ydot = zeros(length(dt))`: Global frame y-velocity for each time step
 - `Zdot = zeros(length(dt))`: Global frame z-velocity for each time step
 - `p = zeros(length(dt))`: Angular velocity about x-axis for each time step
 - `q = zeros(length(dt))`: Angular velocity about y-axis for each time step
 - `r = zeros(length(dt))`: Angular velocity about z-axis for each time step
 - `phi0 = 0`: Roll angle for initial time step
 - `theta0 = 0`: Pitch angle for initial time step
 - `psi0 = 0`: Yaw angle for initial time step
"""
function prescribed_motion(dt;
    Xdot = zeros(length(dt)), Ydot = zeros(length(dt)), Zdot = zeros(length(dt)),
    p = zeros(length(dt)), q = zeros(length(dt)), r = zeros(length(dt)),
    phi0 = 0, theta0 = 0, psi0 = 0)

    # change scalar inputs to vectors
    if isa(Xdot, Number)
        Xdot = fill(Xdot, length(dt))
    end

    if isa(Ydot, Number)
        Ydot = fill(Ydot, length(dt))
    end

    if isa(Zdot, Number)
        Zdot = fill(Zdot, length(dt))
    end

    if isa(p, Number)
        p = fill(p, length(dt))
    end

    if isa(q, Number)
        q = fill(q, length(dt))
    end

    if isa(r, Number)
        r = fill(r, length(dt))
    end

    # get common floating point type
    TF = promote_type(eltype(t), typeof(c), eltype(Xdot), eltype(Ydot), eltype(Zdot),
        eltype(p), eltype(q), eltype(r), typeof(phi0), typeof(theta0), typeof(psi0))

    # number of discrete time steps
    nt = length(dt)

    # velocity magnitude corresponding to each time step
    Vinf = [sqrt(Xdot[it]^2 + Ydot[it]^2 + Zdot[it]^2) for it = 1:nt]

    # freestream parameters corresponding to each time step
    fs = Vector{Freestream{TF}}(undef, nt)

    # set initial orientation
    ϕ = phi0
    θ = theta0
    ψ = psi0

    # populate freestream parameters for each time step
    for it = 1:length(dt)

        sϕ, cϕ = sincos(ϕ)
        sθ, cθ = sincos(θ)
        sψ, cψ = sincos(ψ)

        Rϕ = @SMatrix [1 0 0; 0 cϕ sϕ; 0 -sϕ cϕ]
        Rθ = @SMatrix [cθ 0 -sθ; 0 1 0; sθ 0 cθ]
        Rψ = @SMatrix [cψ sψ 0; -sψ cψ 0; 0 0 1]

        # freestream velocity vector
        V = Rϕ*Rθ*Rψ*SVector(-Xdot[it], -Ydot[it], -Zdot[it])

        # convert to angular representation
        α = atan(V[3]/V[1])
        β = -asin(V[2]/norm(V))
        Ω = SVector(p[it]/Vinf[it], q[it]/Vinf[it], r[it]/Vinf[it])

        # assemble freestream parameters
        fs[it] = Freestream(α, β, Ω)

        if it < length(t)
            # update orientation
            ϕ += p[it]*dt[it]
            θ += q[it]*dt[it]
            ψ += r[it]*dt[it]
        end
    end

    return Vinf, fs
end

"""
    unsteady_analysis(surface, reference, freestream, dx; kwargs...)

Perform an unsteady vortex lattice method analysis.  Return an object of type
[`System`](@ref) containing the final system state, a matrix of surface panel
properties (see [`PanelProperties`](@ref)) at each time step, and a matrix of
wake panels(see [`WakePanel`](@ref)) at each time step.

# Arguments
 - `surface`: Matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the
    number of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters for each time step (see [`Freestream`](@ref))
 - `dx`: Time step vector, multiplied by the freestream velocity for each time step.

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image of the panels
    in `surface` should be used when calculating induced velocities
 - `initial_wake`: Matrix of initial wake panels (see [`WakePanel`](@ref)) of shape
    (nw, ns) where `nw` is the number of chordwise wake panels and `ns` is the
    number of spanwise panels, defaults to no wake panels
 - `initial_circulation`: Vector containing the initial circulation of all surface
    panels in the system.  Defaults to `zeros(N)` where `N` is the total number
    of surface panels in the system.
 - `nwake`: Maximum number of wake panels in the chordwise direction.  Defaults
    to `length(dx)`.
 - `wake_finite_core`: Flag indicating whether the finite core model should be
    enabled when calculating the wake's influence on itself and its corresponding
    surface. Defaults to `true`
 - `save`: Time steps at which to save the time history, defaults to `1:length(dx)`.
 - `calculate_influence_matrix`: Flag indicating whether the aerodynamic influence
    coefficient matrix has already been calculated.  Re-using the same AIC matrix
    will (slightly) reduce calculation times when the underlying geometry has
    not changed. Defaults to `true`.  Note that this argument is only valid for
    the pre-allocated version of this function.
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces for the final
    time step. Defaults to `true`.
 - `derivatives`: Flag indicating whether the derivatives with respect to the
    freestream variables should be calculated for the final time step. Defaults
    to `true`.
"""
function unsteady_analysis(surface::AbstractMatrix, reference, freestream, dx;
    nwake = length(dx), kwargs...)

    system = System(surface; nwake)

    return unsteady_analysis!(system, surface, reference, freestream, dx;
        kwargs..., nwake, calculate_influence_matrix = true)
end

"""
    unsteady_analysis(surfaces, reference, freestream, dx; kwargs...)

Perform a unsteady vortex lattice method analysis.  Return an object of type
[`System`](@ref) containing the final system state, a matrix of surface panel
properties (see [`PanelProperties`](@ref)) for each surface at each time step,
and a matrix of wake panels (see [`WakePanel`](@ref)) for each wake at each time
step.

# Arguments
 - `surfaces`: Vector of surfaces, represented by matrices of surface panels
    (see [`SurfacePanel`](@ref) of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters for each time step (see [`Freestream`](@ref))
 - `dx`: Time step vector, multiplied by the freestream velocity at each time step.

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
    coefficient matrix has already been calculated.  Re-using the same AIC matrix
    will (slightly) reduce calculation times when the underlying geometry has not changed.
    Defaults to `true`.  Note that this argument is only valid for the pre-allocated
    version of this function.
 - `near_field_analysis`: Flag indicating whether a near field analysis should be
    performed to obtain panel velocities, circulation, and forces for the final
    time step. Defaults to `true`.
 - `derivatives`: Flag indicating whether the derivatives with respect to the
    freestream variables should be calculated for the final time step. Defaults
    to `true`.
"""
function unsteady_analysis(surfaces::AbstractVector{<:AbstractMatrix}, reference,
    freestream, dx; nwake = fill(length(dx), length(surfaces)), kwargs...)

    system = System(surfaces; nwake)

    return unsteady_analysis!(system, surfaces, reference, freestream, dx;
        kwargs..., nwake, calculate_influence_matrix = true)
end

"""
    unsteady_analysis!(system, surface, reference, freestream, dx; kwargs...)

Pre-allocated version of `unsteady_analysis`.
"""
function unsteady_analysis!(system, surface::AbstractMatrix, ref, fs, dx;
    symmetric,
    initial_circulation = zeros(length(surface)),
    initial_wake = Matrix{WakePanel{Float64}}(undef, 0, size(surface, 2)),
    nwake = length(dx),
    wake_finite_core = true,
    save = 1:length(dx),
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # float number type
    TF = eltype(system)

    # dimensions
    nc, ns = size(surface)
    nw = nwake

    # convert single freestream input to vector (if applicable)
    if isa(fs, Freestream)
        fs = fill(fs, length(dx))
    end

    # check if existing wake panel storage is sufficient, replace if necessary
    if size(system.wakes[1], 1) < nw
        system.wakes[1] = Matrix{WakePanel{TF}}(undef, nw, ns)
    end

    # copy initial wake panels to pre-allocated storage
    for I in CartesianIndices(initial_wake)
        system.wakes[1][I] = initial_wake[I]
    end

    # copy initial circulation to pre-allocated storage
    system.Γ .= initial_circulation

    # find intersecting surfaces
    repeated_points = repeated_trailing_edge_points(surface)

    # set the initial number of wake panels
    iwake = min(size(initial_wake, 1), nw)

    # disable trailing vortices throughout the analysis
    trailing_vortices = false

    # analysis is unsteady
    unsteady = true

    # hard code distance to wake shedding location (most users probably shouldn't be changing this)
    eta = 0.25

    # unpack pre-allocated storage
    AIC = system.AIC # AIC matrix
    w = system.w # downwash
    Γ = system.Γ # circulation
    wake = system.wakes[1] # wake panels
    wake_velocities = system.V[1] # wake panel vertex velocities
    dw = system.dw # derivatives of downwash wrt freestream variables
    dΓ = system.dΓ # derivatives of circulation wrt freestream variables
    wake_shedding_locations = system.wake_shedding_locations[1] # wake shedding location
    dΓdt = system.dΓdt # derivative of circulation wrt time

    # initialize solution history for each time step
    surface_history = Vector{Matrix{PanelProperties{TF}}}(undef, length(save))
    wake_history = Vector{Matrix{WakePanel{TF}}}(undef, length(save))
    isave = 1

    # loop through all time steps
    for it = 1 : length(dx)

        first_step = it == 1
        last_step = it == length(dx)

        # flag indicating presence of wake panels
        wake_panels = iwake > 0

        # update time step size
        dt = VINF*dx[it]

        # update the wake shedding location
        update_wake_shedding_locations!(wake, wake_shedding_locations,
            surface, fs[it], ref, dt; nwake = iwake, eta = eta)

        # calculate AIC matrix
        if first_step && calculate_influence_matrix
                influence_coefficients!(AIC, surface;
                    symmetric = symmetric,
                    wake_shedding_locations = wake_shedding_locations,
                    trailing_vortices = trailing_vortices)
        end

        # update the AIC matrix to use the new wake shedding locations
        update_trailing_edge_coefficients!(AIC, surface;
            symmetric = symmetric,
            wake_shedding_locations = wake_shedding_locations,
            trailing_vortices = trailing_vortices)

        # calculate downwash due to freestream velocity and aircraft motion
        if last_step && derivatives
            normal_velocity_derivatives!(w, dw, surface, ref, fs[it])
        else
            normal_velocity!(w, surface, ref, fs[it])
        end

        # subtract downwash due to the wake panels
        if wake_panels
            subtract_wake_normal_velocity!(w, surface, wake;
                symmetric = symmetric,
                nwake = iwake,
                wake_finite_core = wake_finite_core,
                trailing_vortices = trailing_vortices)
        end

        # save (-) previous circulation in dΓdt
        dΓdt .= -Γ

        # solve for the new circulation
        if last_step && derivatives
            circulation_derivatives!(Γ, dΓ, AIC, w, dw)
        else
            circulation!(Γ, AIC, w)
        end

        # solve for dΓdt using finite difference `dΓdt = (Γ - Γp)/dt`
        dΓdt .+= Γ # add newly computed circulation
        dΓdt ./= dt # divide by corresponding time step

        # perform a near field analysis to obtain panel properties (if necessary)
        if it in save || (last_step && near_field_analysis)
            if last_step && derivatives
                near_field_forces_derivatives!(system, surface, ref, fs[it];
                    symmetric = symmetric,
                    unsteady = unsteady,
                    wake_shedding_locations = wake_shedding_locations,
                    nwake = iwake,
                    wake_finite_core = wake_finite_core,
                    trailing_vortices = trailing_vortices)
            else
                near_field_forces!(system, surface, ref, fs[it];
                    symmetric = symmetric,
                    unsteady = unsteady,
                    wake_shedding_locations = wake_shedding_locations,
                    nwake = iwake,
                    wake_finite_core = wake_finite_core,
                    trailing_vortices = trailing_vortices)
            end
        end

        # update wake velocities
        VortexLattice.get_wake_velocities!(wake_velocities, surface,
            wake_shedding_locations, wake, ref, fs[it], Γ;
            symmetric = symmetric,
            repeated_points = repeated_points,
            nwake = iwake,
            wake_finite_core = wake_finite_core,
            trailing_vortices = trailing_vortices)

        # shed additional wake panel (and translate existing wake panels)
        VortexLattice.shed_wake!(wake, wake_shedding_locations, wake_velocities,
            dt, surface, Γ; nwake = iwake)

        # increment the wake panel counter
        if iwake < nwake
            iwake += 1
        end

        # save the panel history (if applicable)
        if it in save
            surface_history[isave] = deepcopy(system.panels[1])
            wake_history[isave] = wake[1:iwake, :]
            isave += 1
        end

    end

    # return the modified system and the time history
    return system, surface_history, wake_history
end

"""
    unsteady_analysis!(system, surfaces, reference, freestream, dx; kwargs...)

Pre-allocated version of `unsteady_analysis`.
"""
function unsteady_analysis!(system, surfaces::AbstractVector{<:AbstractMatrix}, ref, fs, dx;
    symmetric,
    initial_circulation = [zeros(size(surfaces[i])) for i = 1:length(surfaces)],
    initial_wakes = [Matrix{WakePanel{Float64}}(undef, 0, size(surfaces[i], 2)) for i = 1:length(surfaces)],
    nwake = fill(length(dx), length(surfaces)),
    surface_id = 1:length(surfaces),
    wake_finite_core = fill(true, length(surfaces)),
    save = 1:length(dx),
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # float number type
    TF = eltype(system)

    # number of surfaces
    nsurf = length(surfaces)

    # convert single freestream input to vector (if applicable)
    if isa(fs, Freestream)
        fs = fill(fs, length(dt))
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

    # copy initial circulation to pre-allocated storage
    system.Γ .= initial_circulation

    # find intersecting surfaces
    repeated_points = repeated_trailing_edge_points(surface)

    # set the initial number of wake panels for each surface
    iwake = [min(size(initial_wakes[isurf], 1), nwake[isurf]) for isurf = 1:nsurf]

    # disable trailing vortices throughout the analysis
    trailing_vortices = fill(false, length(surfaces))

    # analysis is unsteady
    unsteady = true

    # hard code distance to wake shedding location (most users probably shouldn't be changing this)
    eta = 0.25

    # unpack pre-allocated storage
    AIC = system.AIC
    w = system.w
    Γ = system.Γ
    wakes = system.wakes
    wake_velocities = system.V
    dw = system.dw
    dΓ = system.dΓ
    wake_shedding_locations = system.wake_shedding_locations
    dΓdt = system.dΓdt

    # initialize solution history for each time step
    surface_history = Vector{Vector{Matrix{PanelProperties{TF}}}}(undef, length(save))
    wake_history = Vector{Vector{Matrix{WakePanel{TF}}}}(undef, length(save))
    isave = 1

    # # loop through all time steps
    for it = 1:length(dt)

        first_step = it == 1
        last_step = it == 1 + length(dt)

        # flags indicating presence of wake panels
        wake_panels = iwake .> 0

        # update time step size
        dt = VINF*dx[it]

        # update the wake shedding location
        update_wake_shedding_locations!(wakes, wake_shedding_locations,
            surfaces, fs[it], ref, dt; nwake = iwake, eta = eta)

        # calculate AIC matrix
        if first_step && calculate_influence_matrix
            influence_coefficients!(AIC, surfaces;
                symmetric = symmetric,
                wake_shedding_locations = wake_shedding_locations,
                surface_id = surface_id,
                trailing_vortices = trailing_vortices)
        end

        # update the AIC matrix to use the new wake shedding locations
        update_trailing_edge_coefficients!(AIC, surfaces;
            symmetric = symmetric,
            wake_shedding_locations = wake_shedding_locations,
            trailing_vortices = trailing_vortices)

        # calculate downwash due to freestream velocity and aircraft motion
        if last_step && derivatives
            normal_velocity_derivatives!(w, dw, surfaces, ref, fs[it])
        else
            normal_velocity!(w, surfaces, ref, fs[it])
        end

        # subtract downwash due to the wake panels
        if any(wake_panels)
            subtract_wake_normal_velocity!(w, surfaces, wakes;
                symmetric = symmetric,
                nwake = iwake,
                surface_id = surface_id,
                wake_finite_core = wake_finite_core,
                trailing_vortices = trailing_vortices,
                xhat = xhat)
        end

        # save (-) previous circulation in dΓdt
        dΓdt .= -Γ

        # solve for the new circulation
        if last_step && derivatives
            circulation_derivatives!(Γ, dΓ, AIC, w, dw)
        else
            circulation!(Γ, AIC, w)
        end

        # solve for dΓdt using finite difference `dΓdt = (Γ - Γp)/dt`
        dΓdt .+= Γ # add newly computed circulation
        dΓdt ./= dt # divide by corresponding time step

        if it in save || (last_step && near_field_analysis)
            # perform a near field analysis to obtain panel properties
            if last_step && derivatives
                near_field_forces_derivatives!(system, surfaces, ref, fs[it];
                    symmetric = symmetric,
                    wake_shedding_locations = wake_shedding_locations,
                    nwake = iwake,
                    surface_id = surface_id,
                    wake_finite_core = wake_finite_core,
                    trailing_vortices = trailing_vortices,
                    unsteady = unsteady)
            else
                near_field_forces!(system, surfaces, ref, fs[it];
                    symmetric = symmetric,
                    wake_shedding_locations = wake_shedding_locations,
                    nwake = iwake,
                    surface_id = surface_id,
                    wake_finite_core = wake_finite_core,
                    trailing_vortices = trailing_vortices,
                    unsteady = unsteady)
            end
        end

        # update wake velocities
        VortexLattice.get_wake_velocities!(wake_velocities, surfaces,
            wake_shedding_locations, wakes, ref, fs[it], Γ;
            symmetric = symmetric,
            repeated_points = repeated_points,
            nwake = iwake,
            surface_id = surface_id,
            wake_finite_core = wake_finite_core,
            trailing_vortices = trailing_vortices)

        # shed additional wake panel (and translate existing wake panels)
        VortexLattice.shed_wake!(wakes, wake_shedding_locations, wake_velocities,
            dt, surfaces, Γ; nwake = iwake)

        # increment wake panel counter for each surface
        for isurf = 1:nsurf
            if iwake[isurf] < nwake[isurf]
                iwake[isurf] += 1
            end
        end

        # save the panel history
        if it in save
            surface_history[isave] = deepcopy(system.panels)
            wake_history[isave] = [wakes[isurf][1:iwake[isurf], :] for isurf = 1:nsurf]
            isave += 1
        end

    end

    # return the modified system and time history
    return system, surface_history, wake_history
end
