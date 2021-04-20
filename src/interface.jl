"""
    update_surface_panels!(system, surfaces; kwargs...)

Update the surface panels stored in `system` to correspond to the surface panels/
grids in `surfaces`.

# Keyword Arguments
 - `fcore`: function for setting the finite core size when generating surface
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width. Only used if keyword
    argument `surfaces` is provided and is a collection of grids.
 - `preserve_core_size`: Flag indicating whether the finite core size should be
    preserved from the previous set of panels. Defaults to `true`. Only used if
    keyword argument `surfaces` is provided.
"""
function update_surface_panels!(system::AbstractSystem, surfaces; kwargs...)
    return update_surface_panels!(system.surfaces, surfaces; kwargs...)
end

"""
    initialize_wake_panels!(system)

Initialize the undefined wake panels in `wake_panels`.
"""
function initialize_wake_panels!(system::AbstractSystem;
    wake_shedding_locations = system.wake_shedding_locations,
    iwake = system.iwake)
    return initialize_wake_panels!(system.wakes, wake_shedding_locations; iwake)
end

"""
    update_wake_panels!(system, wakes; kwargs...)

Update the wake panels in `system` to correspond to the wake panels/vertices
in `wakes`.

# Keyword Arguments:
 - `helmholtz: Flag indicating whether circulation should be updated based on
    Helmholtz' fourth vortex theorem. (The strength of a vortex tube is constant).
    Defaults to `true`
"""
function update_wake_panels!(system::AbstractSystem, wakes; kwargs...)
    return update_wake_panels!(system.wakes, wakes; kwargs...)
end

"""
    initialize_wake_shedding_locations!(wake_shedding_locations, surfaces)

Initialize the wake shedding locations for one or more surfaces based on the
provided surfaces.
"""
function initialize_wake_shedding_locations!(system::AbstractSystem, surfaces)
    return initialize_wake_shedding_locations!(system.wake_shedding_locations, surfaces)
end

"""
    set_wake_shedding_locations!(wake_shedding_locations, wakes)

Update the wake shedding locations for one or more surfaces based on the
provided wakes.
"""
function set_wake_shedding_locations!(system::AbstractSystem, wakes)
    return get_wake_shedding_locations!(system.wake_shedding_locations, wakes)
end

"""
    update_influence_coefficients!([AIC,] system; kwargs...)

Update the influence coefficients for `system`.

Store the result in `AIC` if provided, otherwise store in `system`.
"""
# steady/unsteady, system input and output
function update_influence_coefficients!(system::AbstractSystem)
    return update_influence_coefficients!(system.AIC, system)
end

# steady, default inputs from system
function update_influence_coefficients!(AIC, system::SteadySystem;
    surfaces = system.surfaces,
    symmetric = system.symmetric,
    surface_id = system.surface_id,
    trailing_vortices = system.trailing_vortices,
    xhat = system.xhat)

    return influence_coefficients!(AIC, surfaces; symmetric, surface_id,
        trailing_vortices, xhat)
end

# unsteady, default inputs from system
function update_influence_coefficients!(AIC, system::UnsteadySystem;
    surfaces = system.surfaces,
    symmetric = system.symmetric,
    surface_id = system.surface_id,
    wake_shedding_locations = system.wake_shedding_locations,
    trailing_vortices = system.trailing_vortices,
    xhat = system.xhat)

    return influence_coefficients!(AIC, surfaces; symmetric, surface_id,
        wake_shedding_locations, trailing_vortices, xhat)
end

"""
    update_trailing_edge_coefficients!([AIC,] system; kwargs...)

Update the trailing edge influence coefficients for `system`.

Store the result in `AIC` if provided, otherwise store in `system`.
"""
update_trailing_edge_coefficients!

# steady/unsteady, system input and output
function update_trailing_edge_coefficients!(system::AbstractSystem)
    return update_trailing_edge_coefficients!(system.AIC, system)
end

# steady, default inputs from system
function update_trailing_edge_coefficients!(AIC, system::SteadySystem;
    surfaces = system.surfaces,
    symmetric = system.symmetric,
    surface_id = system.surface_id,
    trailing_vortices = system.trailing_vortices,
    xhat = system.xhat)

    return trailing_edge_coefficients!(AIC, surfaces; symmetric, surface_id,
        trailing_vortices, xhat)
end

# unsteady, default inputs from system
function update_trailing_edge_coefficients!(AIC, system::UnsteadySystem;
    surfaces = system.surfaces,
    symmetric = system.symmetric,
    surface_id = system.surface_id,
    wake_shedding_locations = system.wake_shedding_locations,
    trailing_vortices = system.trailing_vortices,
    xhat = system.xhat)

    return trailing_edge_coefficients!(AIC, surfaces; symmetric, surface_id,
        wake_shedding_locations, trailing_vortices, xhat)
end

"""
    update_normal_velocities!([w,] system)

Update the normal velocities for `system`.

Store the result in `w` if provided, otherwise store in `system`.
"""
update_normal_velocities!(system::AbstractSystem)

# steady/unsteady, system input and output
function update_normal_velocities!(system::AbstractSystem)
    return update_normal_velocities!(system.w, system)
end

# steady, default inputs from system
function update_normal_velocities!(w, system::SteadySystem;
    surfaces = system.surfaces,
    reference = system.reference,
    freestream = system.freestream)

    return normal_velocities!(w, surfaces, reference, freestream)
end

# unsteady, default inputs from system
function update_normal_velocities!(w, system::UnsteadySystem;
    surfaces = system.surfaces,
    wakes = system.wakes,
    reference = system.reference,
    freestream = system.freestream,
    Vcp = system.Vcp,
    symmetric = system.symmetric,
    surface_id = system.surface_id,
    wake_finite_core = system.wake_finite_core,
    iwake = system.iwake,
    trailing_vortices = system.trailing_vortices,
    xhat = system.xhat)

    return normal_velocities!(w, surfaces, wakes, reference, freestream;
        Vcp, symmetric, surface_id, wake_finite_core, iwake, trailing_vortices, xhat)
end

"""
    update_normal_velocities_and_derivatives!([w, dw,] system)

Update the normal velocities for `system` and their derivatives with respect to
the freestream parameters.

Store the result in `w` and `dw` if provided, otherwise store in `system`.
"""
update_normal_velocities!

# steady/unsteady, system input and output
function update_normal_velocities_and_derivatives!(system::AbstractSystem)
    return update_normal_velocities_and_derivatives!(system.w, system.dw, system)
end

# steady, system input
function update_normal_velocities_and_derivatives!(w, dw, system::SteadySystem;
    surfaces = system.surfaces,
    reference = system.reference,
    freestream = system.freestream)

    return normal_velocities_and_derivatives!(w, dw, surfaces, reference, freestream)
end

# unsteady, system input
function update_normal_velocities_and_derivatives!(w, dw, system::UnsteadySystem;
    surfaces = system.surfaces,
    wakes = system.wakes,
    reference = system.reference,
    freestream = system.freestream,
    Vcp = system.Vcp,
    symmetric = system.symmetric,
    surface_id = system.surface_id,
    wake_finite_core = system.wake_finite_core,
    iwake = system.iwake,
    trailing_vortices = system.trailing_vortices,
    xhat = system.xhat)

    return normal_velocities_and_derivatives!(w, dw, surfaces, wakes,
        ref, fs; Vcp, symmetric, surface_id, wake_finite_core, iwake,
        trailing_vortices, xhat)
end

"""
    update_circulation!([Γ,] system)

Update the circulation strengths of the surface panels in `system`.

Store the result in `Γ` if provided, otherwise store in `system`.

This function assumes that the influence coefficients and normal velocities
have been updated (and are stored in `system`) using [`update_influence_coefficients!`](@ref)
and [`update_normal_velocities!`](@ref), respectively.
"""
update_circulation!

# system inputs and outputs
update_circulation!(system::AbstractSystem) = update_circulation!(system.Gamma, system)

# system inputs
function update_circulation!(Γ, system::AbstractSystem;
    AIC = system.AIC,
    w = system.w)

    return circulation!(Γ, AIC, w)
end

"""
    update_circulation_and_derivatives!([Γ, dΓ,] system)

Update the circulation strengths of the surface panels for `system` and their
derivatives with respect to the freestream variables.

Store the result in `Γ` and `dΓ` if provided, otherwise store in `system`.

This function assumes that the influence coefficients and normal velocities (and
their derivatives) have been updated (and are stored in `system`) using
[`update_influence_coefficients`](@ref) and [`update_normal_velocities_and_derivatives`](@ref),
respectively.
"""
update_circulation_and_derivatives!

# system inputs and outputs
function update_circulation_and_derivatives!(system::AbstractSystem)
    Γ = system.Gamma
    dΓ = system.dGamma
    return update_circulation_and_derivatives!(Γ, dΓ, system)
end

# system inputs
function update_circulation_and_derivatives!(Γ, dΓ, system::AbstractSystem;
    AIC = system.AIC,
    w = system.w,
    dw = system.dw)

    return circulation_and_derivatives!(Γ, dΓ, AIC, w, dw)
end

"""
    update_surface_velocities!([Vcp, Vh, Vv, Vte,] system, surface_velocities)

Update the surface velocities in `system` to correspond to the velocities given
in `surface_velocities`.

Store the result in `Vcp`, `Vh`, `Vv`, and `Vte` if provided, otherwise store in
the corresponding fields in `system`
"""
update_surface_velocities!

# multiple surfaces, surface panels in `system`
function update_surface_velocities!(system::AbstractSystem, surface_velocities)
    return update_surface_velocities!(system.Vcp, system.Vh, system.Vv, system.Vte, system, surface_velocities)
end

# multiple surfaces, provided surface panels
function update_surface_velocities!(Vcp, Vh, Vv, Vte, system::AbstractSystem, surface_velocities)
    return surface_velocities!(Vcp[isurf], Vh[isurf], Vv[isurf], Vte[isurf], surface_velocities[isurf])
end

"""
    update_near_field_properties!([properties,] system)

Update the near field properties for `system`.

Store the result in `properties` if provided, otherwise store in `system`.

This function assumes that the panel circulation (and circulation rates) have
already been computed (and stored in `system`) using [`update_circulation!`](@ref).
"""
update_near_field_properties!

# steady/unsteady, system input and output
function update_near_field_properties!(system::AbstractSystem; kwargs...)
    return update_near_field_properties!(system.properties, system; kwargs...)
end

# steady, system input
function update_near_field_properties!(properties, system::SteadySystem;
    surfaces = system.surfaces,
    reference = system.reference,
    freestream = system.freestream,
    Gamma = system.Gamma,
    symmetric = system.symmetric,
    surface_id = system.surface_id,
    trailing_vortices = system.trailing_vortices,
    xhat = system.xhat)

    return near_field_properties!(properties, surfaces, reference,
        freestream, Gamma; symmetric, surface_id, trailing_vortices, xhat)
end

# unsteady, system input
function update_near_field_properties!(properties, system::UnsteadySystem;
    surfaces = system.surfaces,
    wakes = system.wakes,
    reference = system.reference,
    freestream = system.freestream,
    Gamma = system.Gamma,
    Gammadot = system.Gammadot,
    wake_shedding_locations = system.wake_shedding_locations,
    Vh = system.Vh,
    Vv = system.Vv,
    symmetric = system.symmetric,
    surface_id = system.surface_id,
    trailing_vortices = system.trailing_vortices,
    xhat = system.xhat)

    return near_field_properties!(properties, surfaces, wakes, ref, fs,
        Gamma, Gammadot; wake_shedding_locations, Vh, Vv, symmetric, surface_id,
        wake_finite_core, iwake, trailing_vortices, xhat)
end

"""
    update_near_field_properties_and_derivatives!([properties, dproperties,]
        system)
    update_near_field_properties_and_derivatives!([properties, dproperties,]
        system::UnsteadySystem)

Update the near field properties for `system` and their derivatives with respect
to the freestream variables.

Store the result in `properties` and `dproperties` if provided, otherwise store
in `system`.

This function assumes that the panel circulation, circulation rates, and
circulation derivatives with respect to the freestream variables have been
calculated (and stored in `system`) using [`update_circulation_and_derivatives!`](@ref).
"""
update_near_field_properties_and_derivatives!

# steady/unsteady, system input and output
function update_near_field_properties_and_derivatives!(system::AbstractSystem; kwargs...)
    return update_near_field_properties_and_derivatives!(system.properties, system.dproperties, system; kwargs...)
end

# steady, system input
function update_near_field_properties_and_derivatives!(properties, dproperties, system::SteadySystem;
    surfaces = system.surfaces,
    reference = system.reference,
    freestream = system.freestream,
    Gamma = system.Gamma,
    dGamma = system.dGamma,
    symmetric = system.symmetric,
    surface_id = system.surface_id,
    trailing_vortices = system.trailing_vortices,
    xhat = system.xhat)

    return near_field_properties_and_derivatives!(properties, dproperties,
        surfaces, reference, freestream, Gamma, dGamma; symmetric, surface_id,
        trailing_vortices, xhat)
end

# unsteady, system input
function update_near_field_properties_and_derivatives!(properties, dproperties, system::UnsteadySystem;
    surfaces = system.surfaces,
    wakes = system.wakes,
    reference = system.reference,
    freestream = system.freestream,
    Gamma = system.Gamma,
    Gammadot = system.Gammadot,
    dGamma = system.dGamma,
    wake_shedding_locations = system.wake_shedding_locations,
    Vh = system.Vh,
    Vv = system.Vv,
    symmetric = system.symmetric,
    surface_id = system.surface_id,
    trailing_vortices = system.trailing_vortices,
    xhat = system.xhat)

    return near_field_properties_and_derivatives!(properties, dproperties,
        surfaces, wakes, ref, fs, Gamma, dGamma, Gammadot;
        wake_shedding_locations, Vh, Vv, symmetric, surface_id,
        wake_finite_core, iwake, trailing_vortices, xhat)
end

"""
    body_forces(system::AbstractSystem; frame=Body())

Return the body force coefficients `CF` and `CM` in the frame specified by the
keyword argument `frame`.

Options for keyword argument `frame` are [`Body()`](@ref) (which is the default),
[`Stability()`](@ref), and [`Wind()`](@ref)`

This function assumes that a near-field analysis has been performed (and its
results have been stored in `system`), either by setting the `near_field_analysis`
flag to true (which is the default) when using [`steady_analysis`](@ref) or
[`unsteady_analysis`](@ref) or by using the [`update_near_field_properties`](@ref)
function.
"""
function body_forces(system;
    surfaces = system.surfaces,
    properties = system.properties,
    reference = system.reference,
    freestream = system.freestream,
    symmetric = system.symmetric,
    frame = Body())

    return body_forces(surfaces, properties, reference, freestream, symmetric, frame)
end

"""
    body_forces_derivatives(system::AbstractSystem)

Return the body force coefficients `CF` and `CM` and their derivatives with
respect to the freestream variables in the body frame.

This function assumes that a near-field analysis has been performed with freestream
parameter derivatives (and has been stored in `system`), either by setting the
`near_field_analysis` and `derivatives` flags to true (which is the default) when
using [`steady_analysis`](@ref) or [`unsteady_analysis`](@ref) or by using the
[`update_near_field_properties_and_derivatives`](@ref) function.
"""
function body_forces_derivatives(system::AbstractSystem;
    surfaces = system.surfaces,
    properties = system.properties,
    dproperties = system.dproperties,
    reference = system.reference,
    freestream = system.freestream,
    symmetric = system.symmetric
    )

    return body_forces_derivatives(surfaces, properties, dproperties, reference, freestream, symmetric)
end

"""
    body_forces_history(system, surface_history, property_history; frame=Body())

Return the body force coefficients `CF`, `CM` at each time step in `property_history`.

# Arguments:
 - `system`: Object of type [`SteadySystem`](@ref) which holds system properties
 - `surface_history`: Vector of surfaces at each time step, where each surface is
    represented by a matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `property_history`: Vector of surface properties for each surface at each
    time step, where surface properties are represented by a matrix of panel
    properties (see [`PanelProperties`](@ref)) of shape (nc, ns) where `nc` is
    the number of chordwise panels and `ns` is the number of spanwise panels

# Keyword Arguments
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)`
"""
function body_forces_history(system::AbstractSystem, surface_history, property_history;
    symmetric = system.symmetric,
    reference = system.reference,
    freestream = system.freestream,
    frame=Body())

    return body_forces_history(surface_history, property_history, reference,
        freestream, symmetric, frame)
end

"""
    lifting_line_coefficients(system, r, c)

Return the force and moment coefficients (per unit span) for each spanwise segment
of a lifting line representation of the geometry.

# Arguments
 - `system`: Object of type [`SteadySystem`](@ref) that holds precalculated
    system properties.
 - `r`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns+1) which contains the x, y, and z coordinates
    of the resulting lifting line coordinates
 - `c`: Vector with length equal to the number of surfaces, with each element
    being a vector of length `ns+1` which contains the chord lengths at each
    lifting line coordinate.

# Return Arguments:
 - `cf`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns) which contains the x, y, and z direction
    force coefficients (per unit span) for each spanwise segment.
 - `cm`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns) which contains the x, y, and z direction
    moment coefficients (per unit span) for each spanwise segment.

This function assumes that a near-field analysis has been performed (and its
results have been stored in `system`), either by setting the `near_field_analysis`
flag to true (which is the default) when using [`steady_analysis`](@ref) or
[`unsteady_analysis`](@ref) or by using the [`update_near_field_properties`](@ref)
function.
"""
function lifting_line_coefficients(system, r, c;
    surfaces = system.surfaces,
    properties = system.properties,
    reference = system.reference)

    return lifting_line_coefficients(surfaces, properties, reference, r, c)
end

"""
    lifting_line_coefficients!(cf, cm, system, r, c)

In-place version of [`lifting_line_coefficients`](@ref)
"""
function lifting_line_coefficients!(cf, cm, system, r, c;
    surfaces = system.surfaces,
    properties = system.properties,
    reference = system.reference)

    return lifting_line_coefficients!(cf, cm, surfaces, properties, reference, r, c)
end

"""
    far_field_drag(system)

Computes induced drag using the Trefftz plane (far field method).

This function assumes that panel circulation strengths have been calculated
(and are stored in `system`) by using [`steady_analysis`](@ref),
[`unsteady_analysis`](@ref), [`update_circulation!`](@ref), or
[`update_circulation_and_derivatives!`](@ref).
"""
function far_field_drag(system;
    surfaces = system.surfaces,
    trefftz = system.trefftz,
    reference = system.reference,
    freestream = system.freestream,
    symmetric = system.symmetric,
    Gamma = system.Gamma)

    return far_field_drag(trefftz, surfaces, reference, freestream, symmetric, Gamma)
end

"""
    body_derivatives(system)

Returns the derivatives of the body forces and moments with respect to the
freestream velocity components `(u, v, w)` and the angular velocity components
`(p, q, r)` in the body frame.

The derivatives are returned as two named tuples: `dCF, dCM`

This function assumes that a near-field analysis has been performed with freestream
parameter derivatives (and has been stored in `system`), either by setting the
`near_field_analysis` and `derivatives` flags to true (which is the default) when
using [`steady_analysis`](@ref) or [`unsteady_analysis`](@ref) or by using the
[`update_near_field_properties_and_derivatives`](@ref) function.
"""
function body_derivatives(system;
    surfaces = system.surfaces,
    properties = system.properties,
    dproperties = system.dproperties,
    reference = system.reference,
    freestream = system.freestream,
    symmetric = system.symmetric)

    return body_derivatives(surfaces, properties, dproperties, reference, freestream, symmetric)
end

"""
    stability_derivatives(system)

Returns the derivatives of the body forces and moments in the stability frame
with respect to the freestream velocity components `(alpha, beta)` and the angular
velocity components `(p, q, r)` in the stability frame.

The derivatives are returned as two named tuples: `dCF, dCM`

This function assumes that a near-field analysis has been performed with freestream
parameter derivatives (and has been stored in `system`), either by setting the
`near_field_analysis` and `derivatives` flags to true (which is the default) when
using [`steady_analysis`](@ref) or [`unsteady_analysis`](@ref) or by using the
[`update_near_field_properties_and_derivatives`](@ref) function.
"""
function stability_derivatives(system;
    surfaces = system.surfaces,
    properties = system.properties,
    dproperties = system.dproperties,
    reference = system.reference,
    freestream = system.freestream,
    symmetric = system.symmetric)

    return stability_derivatives(surfaces, properties, dproperties, reference, freestream, symmetric)
end

"""
    update_wake_velocities!([wake_velocities,] system)

Update the wake velocities for `system`.

Store the result in `wake_velocities` if provided, otherwise store in `system`.
"""
update_wake_velocities!

# unsteady, system input and output
update_wake_velocities!(system::AbstractSystem) = update_wake_velocities!(system.Vw, system)

# unsteady, system input
function update_wake_velocities!(Vw, system;
    surfaces = system.surfaces,
    repeated_points = system.repeated_points,
    wake_shedding_locations = system.wake_shedding_locations,
    wakes = system.wakes,
    reference = system.reference,
    freestream = system.freestream,
    Gamma = system.Gamma,
    Vte = system.Vte,
    symmetric = system.symmetric,
    surface_id = system.surface_id,
    wake_finite_core = system.wake_finite_core,
    iwake = system.iwake,
    trailing_vortices = system.trailing_vortices,
    xhat = system.xhat,
    )

    return wake_velocities!(Vw, surfaces, repeated_points,
        wake_shedding_locations, wakes, reference, freestream, Gamma; Vte,
        symmetric, surface_id, wake_finite_core, iwake, trailing_vortices, xhat)
end
