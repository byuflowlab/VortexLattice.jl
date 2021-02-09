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
steady_analysis

function steady_analysis(grids::AbstractVector{<:AbstractArray{<:Any, 3}},
    reference, freestream; fcore = (c, Δs) -> 1e-3, kwargs...)

    # pre-allocate system storage
    system = System(grids)

    # generate surface panels
    for i = 1:length(grids)
        update_surface_panels!(system.surfaces[i], grids[i]; fcore)
    end

    return steady_analysis!(system, system.surfaces, reference, freestream; kwargs...,
        calculate_influence_matrix = true)
end

function steady_analysis(surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}},
    reference, freestream; kwargs...)

    # pre-allocate system storage
    system = System(surfaces)

    return steady_analysis!(system, surfaces, reference, freestream; kwargs...,
        calculate_influence_matrix = true)
end

"""
    steady_analysis!(system, surfaces, reference, freestream; kwargs...)

Pre-allocated version of `steady_analysis`.
"""
function steady_analysis!(system, surfaces::AbstractVector{<:AbstractMatrix},
    ref, fs;
    symmetric = fill(false, length(surfaces)),
    wakes = [Matrix{WakePanel{Float64}}(undef, 0, size(surfaces[i],2)) for i = 1:length(surfaces)],
    nwake = size.(wakes, 1),
    surface_id = 1:length(surfaces),
    wake_finite_core = fill(true, length(surfaces)),
    trailing_vortices = fill(true, length(surfaces)),
    xhat = SVector(1, 0, 0),
    additional_velocity = nothing,
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # number of surfaces
    nsurf = length(surfaces)

    # convert scalar inputs to vectors
    symmetric = isa(symmetric, Number) ? fill(symmetric, nsurf) : symmetric
    nwake = isa(nwake, Number) ? fill(nwake, nsurf) : nwake
    surface_id = isa(surface_id, Number) ? fill(surface_id, nsurf) : surface_id
    wake_finite_core = isa(wake_finite_core, Number) ? fill(wake_finite_core, nsurf) : wake_finite_core
    trailing_vortices = isa(trailing_vortices, Number) ? fill(trailing_vortices, nsurf) : trailing_vortices

    # update parameters stored in `system`
    system.surfaces .= surfaces
    system.reference[] = ref
    system.freestream[] = fs
    system.symmetric .= symmetric
    system.wakes .= wakes
    system.nwake .= nwake
    system.surface_id .= surface_id
    system.wake_finite_core .= wake_finite_core
    system.trailing_vortices .= trailing_vortices
    system.xhat[] = xhat

    # unpack system variables
    AIC = system.AIC
    w = system.w
    Γ = system.Γ
    dw = system.dw
    dΓ = system.dΓ

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
        normal_velocity_derivatives!(w, dw, surfaces, ref, fs, additional_velocity)
    else
        normal_velocity!(w, surfaces, ref, fs, additional_velocity)
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
            near_field_forces_derivatives!(system; additional_velocity)
        else
            near_field_forces!(system; additional_velocity)
        end
    end

    # save flags indicating whether certain analyses have been performed
    system.near_field_analysis[] = near_field_analysis
    system.derivatives[] = derivatives

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
    TF = promote_type(eltype(dt), eltype(Xdot), eltype(Ydot), eltype(Zdot),
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

        if it < length(dt)
            # update orientation
            ϕ += p[it]*dt[it]
            θ += q[it]*dt[it]
            ψ += r[it]*dt[it]
        end
    end

    return Vinf, fs
end

"""
    unsteady_analysis(surfaces, reference, freestream, dx; kwargs...)

Perform a unsteady vortex lattice method analysis.  Return an object of type
[`System`](@ref) containing the final system state, a matrix of surface panel
properties (see [`PanelProperties`](@ref)) for each surface at each time step,
and a matrix of wake panels (see [`WakePanel`](@ref)) for each wake at each time
step.

# Arguments
 - `surfaces`:
   - Grids of shape (3, nc+1, ns+1) which represent lifting surfaces
     or
   - Matrices of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns)
    where `nc` is the number of chordwise panels and `ns` is the number of
    spanwise panels.  Moving lifting surfaces may be modeled by passing in a
    vector containing `surfaces` at each time step.
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
unsteady_analysis

# grid inputs
function unsteady_analysis(grids::AbstractVector{<:AbstractArray{<:Any, 3}},
    reference, freestream;
    fcore = (c, Δs) -> 1e-3,
    nwake = fill(length(dx), length(surfaces)),
    kwargs...)

    # pre-allocate system storage
    system = System(grids; nw = nwake)

    # generate surface panels
    for i = 1:length(grids)
        update_surface_panels!(system.surfaces[i], grids[i]; fcore)
    end

    return unsteady_analysis!(system, system.surfaces, reference, freestream, dx;
        kwargs..., nwake, calculate_influence_matrix = true)
end

# surface panel inputs
function unsteady_analysis(surfaces::AbstractVector{<:AbstractMatrix}, reference,
    freestream, dx; nwake = fill(length(dx), length(surfaces)), kwargs...)

    system = System(surfaces; nw = nwake)

    return unsteady_analysis!(system, surfaces, reference, freestream, dx;
        kwargs..., nwake, calculate_influence_matrix = true)
end

"""
    unsteady_analysis!(system, surfaces, reference, freestream, dx; kwargs...)

Pre-allocated version of `unsteady_analysis`.
"""
function unsteady_analysis!(system, surfaces::AbstractVector{<:AbstractMatrix}, ref, fs, dx;
    symmetric = fill(false, length(surfaces)),
    initial_circulation = zero(system.Γ),
    initial_wakes = [Matrix{WakePanel{Float64}}(undef, 0, size(surfaces[i], 2)) for i = 1:length(surfaces)],
    nwake = fill(length(dx), length(surfaces)),
    surface_id = 1:length(surfaces),
    wake_finite_core = fill(true, length(surfaces)),
    additional_velocity = nothing,
    save = 1:length(dx),
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # float number type
    TF = eltype(system)

    # number of surfaces
    nsurf = length(surfaces)

    # convert scalar inputs to vectors
    symmetric = isa(symmetric, Number) ? fill(symmetric, nsurf) : symmetric
    nwake = isa(nwake, Number) ? fill(nwake, nsurf) : nwake
    surface_id = isa(surface_id, Number) ? fill(surface_id, nsurf) : surface_id
    wake_finite_core = isa(wake_finite_core, Number) ? fill(wake_finite_core, nsurf) : wake_finite_core

    # convert single freestream input to vector (if applicable)
    if isa(fs, Freestream)
        fs = fill(fs, length(dx))
    end

    # update parameters stored in `system`

    # copy initial surface panels to pre-allocated storage
    for isurf = 1:nsurf
        system.surfaces[isurf] .= surfaces[isurf]
    end

    system.reference[] = ref
    system.freestream[] = fs[1]
    system.symmetric .= symmetric

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

    system.nwake .= nwake
    system.surface_id .= surface_id
    system.wake_finite_core .= wake_finite_core
    system.trailing_vortices .= false
    system.Γ .= initial_circulation

    # unpack pre-allocated storage
    AIC = system.AIC
    w = system.w
    Γ = system.Γ
    surfaces = system.surfaces
    wakes = system.wakes
    wake_velocities = system.V
    dw = system.dw
    dΓ = system.dΓ
    wake_shedding_locations = system.wake_shedding_locations
    dΓdt = system.dΓdt

    # find intersecting surfaces
    repeated_points = repeated_trailing_edge_points(surfaces)

    # set the initial number of wake panels for each surface
    iwake = [min(size(initial_wakes[isurf], 1), nwake[isurf]) for isurf = 1:nsurf]

    # analysis is unsteady
    unsteady = true

    # trailing vortices are disabled
    trailing_vortices = fill(false, nsurf)

    # hard code distance to wake shedding location (most users probably shouldn't be changing this)
    eta = 0.25

    # initialize solution history for each time step
    surface_history = Vector{Vector{Matrix{SurfacePanel{TF}}}}(undef, length(save))
    property_history = Vector{Vector{Matrix{PanelProperties{TF}}}}(undef, length(save))
    wake_history = Vector{Vector{Matrix{WakePanel{TF}}}}(undef, length(save))
    isave = 1

    # # loop through all time steps
    for it = 1 : length(dx)

        first_step = it == 1
        last_step = it == length(dx)

        # flags indicating presence of wake panels
        wake_panels = iwake .> 0

        # update time step size
        dt = VINF*dx[it]

        # update freestream parameters
        system.freestream[] = fs[it]

        # update current number of wake panels
        system.nwake .= iwake

        # update the wake shedding location
        update_wake_shedding_locations!(wakes, wake_shedding_locations,
            surfaces, ref, fs[it], additional_velocity, dt; nwake = iwake, eta = eta)

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
            normal_velocity_derivatives!(w, dw, surfaces, ref, fs[it], additional_velocity)
        else
            normal_velocity!(w, surfaces, ref, fs[it], additional_velocity)
        end

        # subtract downwash due to the wake panels
        if any(wake_panels)
            subtract_wake_normal_velocity!(w, surfaces, wakes;
                symmetric = symmetric,
                nwake = iwake,
                surface_id = surface_id,
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

        if it in save || (last_step && near_field_analysis)
            # perform a near field analysis to obtain panel properties
            if last_step && derivatives
                near_field_forces_derivatives!(system; additional_velocity, unsteady)
            else
                near_field_forces!(system; additional_velocity, unsteady)
            end
        end

        # save the panel history
        if it in save
            surface_history[isave] = deepcopy(system.surfaces)
            property_history[isave] = deepcopy(system.properties)
            wake_history[isave] = [wakes[isurf][1:iwake[isurf], :] for isurf = 1:nsurf]
            isave += 1
        end

        # update wake velocities
        VortexLattice.get_wake_velocities!(wake_velocities, surfaces,
            wake_shedding_locations, wakes, ref, fs[it], additional_velocity, Γ;
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

    end

    # save flags indicating whether certain analyses have been performed
    system.near_field_analysis[] = near_field_analysis
    system.derivatives[] = derivatives

    # return the modified system and time history
    return system, surface_history, property_history, wake_history
end
