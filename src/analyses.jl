"""
    steady_analysis(surface, reference, freestream; kwargs...)

Perform a steady vortex lattice method analysis to find the circulation
given a specific combination of vortex lattice panels.

# Arguments
 - `system`: Pre-allocated system properties
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see `Reference`)
 - `freestream`: Freestream parameters (see `Freestream`)

# Keyword Arguments
 - `wake`: Matrix of wake panels of shape (nw, ns) where `nw` is the number of
    chordwise wake panels and `ns` is the number of spanwise panels, defaults to
    no wake panels
 - `symmetric`: Flag indicating whether a mirror image of the panels in `surface`,
    should be used when calculating induced velocities, defaults to `false`
 - `trailing_vortices`: Flag to enable/disable trailing vortices, defaults to `true`
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
function steady_analysis(surface::AbstractMatrix, reference, freestream; kwargs...)

    system = System(surface)

    return steady_analysis!(system, surface, reference, freestream; kwargs...,
        calculate_influence_matrix = true)
end

"""
    steady_analysis(surfaces, reference, freestream; kwargs...)

Perform a steady vortex lattice method analysis to find the circulation
given a specific combination of vortex lattice panels.

# Arguments
 - `system`: Pre-allocated system properties
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: Reference parameters (see `Reference`)
 - `freestream`: Freestream parameters (see `Freestream`)

# Keyword Arguments
 - `wakes`: Vector of wakes corresponding to each surface, represented by matrices
    of panels of shape (nw, ns) where `nw` is the number of chordwise wake panels
    and `ns` is the number of spanwise panels
 - `symmetric`: Flags indicating whether a mirror image (across the y-axis) should
    be used when calculating induced velocities, defaults to `false` for each surface
 - `trailing_vortices`: Flags to enable/disable trailing vortices, defaults to `true`
    for each surface
 - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]
 - `surface_id`: ID for each surface.  May be used to deactivate the finite core
    model by setting all surface ID's to the same value. By default all surfaces
    have their own IDs
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
    wake = Matrix{Wake{Float64}}(undef, 0, size(surface,2)),
    symmetric = false,
    trailing_vortices = true,
    xhat = SVector(1, 0, 0),
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # unpack system variables
    AIC = system.AIC
    b = system.b
    Γ = system.gamma
    db = system.db
    dΓ = system.dgamma

    # number of wake panels
    nwake = size(wake, 1)

    # see if wake panels are being used
    wake_panels = nwake > 0

    # store new wake panels in system
    system.wakes[1] = wake

    # calculate AIC matrix
    if calculate_influence_matrix
        influence_coefficients!(AIC, surface, symmetric;
            trailing_vortices=trailing_vortices && !wake_panels, xhat=xhat)
    end

    # calculate RHS
    if derivatives
        normal_velocity_derivatives!(b, db, surface, ref, fs)
    else
        normal_velocity!(b, surface, ref, fs)
    end

    # add wake's contribution to RHS
    if wake_panels
        add_wake_normal_velocity!(b, surface, wake, symmetric;
            nwake = nwake,
            trailing_vortices = trailing_vortices,
            xhat = xhat)
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
                xhat = xhat)
        else
            near_field_forces!(system, surface, ref, fs;
                symmetric = symmetric,
                nwake = nwake,
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
function steady_analysis!(system, surfaces::AbstractVector{<:AbstractMatrix}, ref, fs;
    wakes = [Matrix{Wake{Float64}}(undef, 0, size(surfaces[i],2)) for i = 1:length(surfaces)],
    symmetric = fill(false, length(surfaces)),
    trailing_vortices = fill(true, length(surfaces)),
    xhat = SVector(1, 0, 0),
    surface_id = 1:length(surfaces),
    calculate_influence_matrix = true,
    near_field_analysis = true,
    derivatives = true)

    # unpack system variables
    AIC = system.AIC
    b = system.b
    Γ = system.gamma
    db = system.db
    dΓ = system.dgamma

    # number of wake panels
    nwake = size.(wakes, 1)

    # see if wake panels are being used
    wake_panels = nwake .> 0

    # store new wake panels in system
    for i = 1:length(wakes)
        system.wakes[i] = wakes[i]
    end

    # calculate AIC matrix
    if calculate_influence_matrix
        influence_coefficients!(AIC, surfaces, surface_id, symmetric;
            trailing_vortices = trailing_vortices .& .!wake_panels,
            xhat = xhat)
    end

    # calculate RHS
    if derivatives
        normal_velocity_derivatives!(b, db, surfaces, ref, fs)
    else
        normal_velocity!(b, surfaces, ref, fs)
    end

    # add wake's contribution to RHS
    if any(wake_panels)
        add_wake_normal_velocity!(b, surfaces, wakes, symmetric;
            trailing_vortices = trailing_vortices,
            nwake = nwake,
            xhat = xhat)
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
                surface_id = surface_id)
        else
            near_field_forces!(system, surfaces, ref, fs;
                symmetric = symmetric,
                nwake = nwake,
                trailing_vortices = trailing_vortices,
                xhat = xhat,
                surface_id = surface_id)
        end
    end

    # return the modified system
    return system
end
#
# """
#     unsteady_analysis(surface, reference, freestream; kwargs...)
#
# Perform a unsteady vortex lattice method analysis to find the unsteady
# circulation and wake shape of a group of vortex lattice panels.
#
# # Arguments
#  - `system`: Pre-allocated system properties
#  - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
#     chordwise panels and `ns` is the number of spanwise panels
#  - `reference`: Reference parameters (see `Reference`)
#  - `freestream`: Freestream parameters (see `Freestream`)
#
# # Keyword Arguments
#
#  - `dt`: Step size for each time step, defaults to 0.1
#  - `nstep`: Total number of time steps, defaults to 100
#  - `nwake`: Maximum number of wake panels in the chordwise direction.  Defaults
#     to 100.
#  - `iwake`: Initial number of wake panels in the chordwise direction. Defaults to
#     zero initial wake panels
#  - `wake`: Pre-initialized matrix of wake panels of shape (nwake, ns) where
#     `nwake` is the maximum number of chordwise wake panels and `ns` is the number
#     of spanwise panels, defaults to
#  - `symmetric`: Flag indicating whether a mirror image of the panels in `surface`,
#     should be used when calculating induced velocities, defaults to `false`
#  - `trailing_vortices`: Flag to enable/disable trailing vortices, defaults to `true`
#  - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]
# """
# function unsteady_analysis(surface::AbstractMatrix, reference, freestream; kwargs...)
#
#     system = System(surface)
#
#     return steady_analysis!(surface, reference, freestream; kwargs...)
# end
#
# function unsteady_analysis!(system, surface, wake, symmetric, ref, fs; dt=0.1, nsteps=1,
#     trailing_vortices=false, xhat=SVector(1, 0, 0))
#
#     nw, ns = size(wake)
#
#     Vwake = zeros(eltype(eltype(wake)), 3, nw+1, ns+1)
#
#     for istep = 1:nsteps
#         nwake = min(istep, size(wake,1))
#
#         AIC = influence_coefficients(surface, symmetric; xhat=xhat)
#         b = normal_velocity(surface, ref, fs)
#         b = add_wake_normal_velocity!(b, surface, wake, trailing_vortices, symmetric;
#             nwake=nwake, xhat=xhat)
#         Γ = circulation(AIC, b)
#         wake_velocities!(V, surface, wake, false, symmetric, ref, fs)
#         translate_wake!(wake, V, dt; nwake=nwake)
#         shed_wake!(wake, V, dt; nwake=nwake)
#     end
#
#     return surface, wake
# end
