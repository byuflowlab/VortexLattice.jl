"""
    near_field_forces!(system, surface, ref, fs; kwargs...)

Perform a near-field analysis to calculate panel forces in the body frame.  Return
the modified system containing the calculated panel forces.

Note that this function assumes that the circulation distribution has already
been calculated and stored in `system`.

# Arguments:
 - `system`: Object of type [`System`](@ref) which holds system properties
 - `surface`: Matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: reference parameters (see [`Reference`](@ref))
 - `freestream`: freestream parameters (see [`Freestream`](@ref))

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image (across the X-Z plane) of
    the panels in `surface` should be used when calculating induced velocities.
 - `unsteady`: Flag indicating whether to include unsteady forces. Defaults to `false`.
 - `wake_shedding_locations`: Wake shedding locations for each surface.  Defaults
    to the trailing edge.
 - `nwake`: number of chordwise wake panels to use from the wake stored in `system`.
    Defaults to all wake panels in `system`.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.  Defaults to `true` for each surface.
 - `trailing_vortices`: Flag indicating whether trailing vortices should be shed
    from the last chordwise panel in `wake`
 - `xhat`: direction in which trailing vortices are shed
"""
@inline function near_field_forces!(system, surface::AbstractMatrix, ref, fs;
    symmetric = false,
    wake_shedding_locations = nothing,
    nwake = size(system.wakes[1], 1),
    wake_finite_core = true,
    trailing_vortices = true,
    kwargs...)

    return near_field_forces!(system, [surface], ref, fs;
        symmetric = [symmetric],
        wake_shedding_locations = [wake_shedding_locations],
        nwake = [nwake],
        wake_finite_core = [wake_finite_core],
        trailing_vortices = [trailing_vortices],
        kwargs...,
        surface_id = 1:1)
end

"""
    near_field_forces!(system, surfaces, reference, freestream; kwargs...)

Perform a near-field analysis to calculate panel forces in the body frame.  Return
the modified `system` containing the calculated panel forces.

Note that this function assumes that the circulation distribution has already
been calculated and stored in `system`.

# Arguments
 - `system`: Object of type [`System`](@ref) which holds system properties
 - `surfaces`: Vector of surfaces, represented by matrices of surface panels
    (see [`SurfacePanel`](@ref) of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`]@ref)

# Keyword Arguments
 - `symmetric`: (required) Flag for each surface indicating whether a mirror image
    (across the X-Z plane) should be used when calculating induced velocities
 - `unsteady`: Flag indicating whether to include unsteady forces. Defaults to `false`.
 - `wake_shedding_locations`: Wake shedding locations for each surface.  Defaults
    to the trailing edge.
 - `nwake`: Number of chordwise wake panels to use from each wake in `wakes`,
    Defaults to all wake panels in `system`
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.  Defaults to `true` for each surface.
 - `trailing_vortices`: Flags to enable/disable trailing vortices, defaults to
    `true` for each surface
 - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]

"""
@inline function near_field_forces!(system, surfaces::AbstractVector{<:AbstractMatrix},
    ref, fs;
    symmetric = fill(false, length(surfaces)),
    unsteady = false,
    wake_shedding_locations = fill(nothing, length(surfaces)),
    nwake = size.(system.wakes, 1),
    surface_id = 1:length(surfaces),
    wake_finite_core = fill(true, length(surfaces)),
    trailing_vortices = fill(true, length(surfaces)),
    xhat = SVector(1, 0, 0))

    # unpack system storage
    Γ = system.Γ
    dΓdt = system.dΓdt

    props = system.panels

    # number of surfaces
    nsurf = length(surfaces)

    # see if wake panels are being used
    wake_panels = nwake .> 0

    # loop through receiving surfaces
    iΓ = 0 # index for accessing Γ
    for isurf = 1:nsurf

        receiving = surfaces[isurf]
        nr = length(receiving)
        nr1, nr2 = size(receiving)
        cr = CartesianIndices(receiving)

        # loop through receiving panels
        for i = 1:length(receiving)

            I = cr[i]

            # --- Calculate forces on the panel bound vortex -- #
            rc = top_center(receiving[I])

            # start with external velocity at the bound vortex center
            Vi = external_velocity(fs, rc, ref.r)

            # now add the velocity induced by each sending panel
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                sending = surfaces[jsurf]
                ns = length(sending)

                # extract circulation values corresonding to the sending surface
                vΓ = view(Γ, jΓ+1:jΓ+ns)

                # add induced velocity from the sending surface
                if isurf == jsurf
                    # special induced velocity calculation to avoid bound vortex
                    Vi += induced_velocity(I, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = wake_shedding_locations[jsurf],
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels[jsurf],
                        xhat = xhat)
                else
                    # regular induced velocity calculation
                    Vi += induced_velocity(rc, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = wake_shedding_locations[jsurf],
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels[jsurf],
                        xhat = xhat)
                end

                # add induced velocity from wake panels
                if wake_panels[jsurf]
                    Vi += induced_velocity(rc, system.wakes[jsurf];
                        symmetric = symmetric[jsurf],
                        finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end

                jΓ += ns # increment Γ index for sending panels
            end

            # circulation of current panel (and its time derivative)
            if isone(I[1])
                Γi = Γ[iΓ+i]
                dΓdti = dΓdt[iΓ+i]
            else
                Γi = Γ[iΓ+i] - Γ[iΓ+i-1]
                dΓdti = (dΓdt[iΓ+i] + dΓdt[iΓ+i-1])/2
            end

            # circulation vector
            Δs = top_vector(receiving[I])

            # now use the Kutta-Joukowski theorem to solve for the panel force
            tmp = cross(Vi, Δs)

            # steady part of Kutta-Joukowski theorem
            Fbi = RHO*Γi*tmp

            # unsteady part of Kutta-Joukowski theorem
            if unsteady
                # panel chord length
                c = receiving[I].chord

                #TODO: decide whether to divide by perpindicular velocity like
                # Drela does in ASWING?

                # unsteady part of Kutta-Joukowski theorem
                Fbi += RHO*dΓdti*c*tmp
            end

            # --- Calculate forces on the left bound vortex --- #

            rc = left_center(receiving[I])

            # get effective velocity at the bound vortex center
            Veff = external_velocity(fs, rc, ref.r)

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # now use the Kutta-Joukowski theorem to calculate the forces
            Γli = Γ[iΓ+i]
            Δs = left_vector(receiving[I])
            Fbli = RHO*Γli*cross(Veff, Δs)

            # --- Calculate forces on the right bound vortex --- #

            rc = right_center(receiving[I])

            # get effective velocity at the bound vortex center
            Veff = external_velocity(fs, rc, ref.r)

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # now use the Kutta-Joukowski theorem to calculate the forces
            Γri = Γ[iΓ+i]
            Δs = right_vector(receiving[I])
            Fbri = RHO*Γri*cross(Veff, Δs)

            # store panel circulation, velocity, and forces
            props[isurf][I] = PanelProperties(Γ[iΓ+i]/VINF, Vi/VINF,
                Fbi/(QINF*ref.S), Fbli/(QINF*ref.S), Fbri/(QINF*ref.S))
        end

        # increment Γ index for receiving panels
        iΓ += nr
    end

    return system
end

"""
    near_field_forces_derivatives!(system, surface[s], ref, fs; kwargs...)

Version of [`near_field_forces!`](@ref) that also calculates the derivatives of
the panel forces with respect to the freestream variables.
"""
near_field_forces_derivatives!

# single surface
@inline function near_field_forces_derivatives!(system, surface::AbstractMatrix, ref, fs;
    symmetric = false,
    wake_shedding_locations = nothing,
    nwake = size(system.wakes[1], 1),
    wake_finite_core = true,
    trailing_vortices = true,
    kwargs...)

    return near_field_forces_derivatives!(system, [surface], ref, fs;
        symmetric = [symmetric],
        wake_shedding_locations = [wake_shedding_locations],
        nwake = [nwake],
        wake_finite_core = [wake_finite_core],
        trailing_vortices = [trailing_vortices],
        kwargs...,
        surface_id = 1:1)
end

# multiple surfaces
@inline function near_field_forces_derivatives!(system, surfaces::AbstractVector{<:AbstractMatrix},
    ref, fs;
    symmetric = fill(false, length(surfaces)),
    unsteady = false,
    wake_shedding_locations = fill(nothing, length(surfaces)),
    nwake = size.(system.wakes, 1),
    surface_id = 1:length(surfaces),
    wake_finite_core = fill(true, length(surfaces)),
    trailing_vortices = fill(true, length(surfaces)),
    xhat = SVector(1, 0, 0))

    # unpack system storage
    Γ = system.Γ
    Γ_a, Γ_b, Γ_p, Γ_q, Γ_r = system.dΓ
    dΓdt = system.dΓdt

    props = system.panels
    props_a, props_b, props_p, props_q, props_r = system.dpanels

    # number of surfaces
    nsurf = length(surfaces)

    # see if wake panels are being used
    wake_panels = nwake .> 0

    # loop through receiving surfaces
    iΓ = 0 # index for accessing Γ
    for isurf = 1:nsurf

        receiving = surfaces[isurf]
        nr = length(receiving)
        nr1, nr2 = size(receiving)
        cr = CartesianIndices(receiving)

        # loop through receiving panels
        for i = 1:length(receiving)

            I = cr[i]

            # --- Calculate forces on the panel bound vortex -- #
            rc = top_center(receiving[I])

            # start with external velocity at the bound vortex center
            Vi, dVi = external_velocity_derivatives(fs, rc, ref.r)

            # unpack derivatives
            Vi_a, Vi_b, Vi_p, Vi_q, Vi_r = dVi

            # now add the velocity induced by each sending panel
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                sending = surfaces[jsurf]
                ns = length(sending)

                # extract circulation values corresonding to the sending surface
                vΓ = view(Γ, jΓ+1:jΓ+ns)

                vΓ_a = view(Γ_a, jΓ+1:jΓ+ns)
                vΓ_b = view(Γ_b, jΓ+1:jΓ+ns)
                vΓ_p = view(Γ_p, jΓ+1:jΓ+ns)
                vΓ_q = view(Γ_q, jΓ+1:jΓ+ns)
                vΓ_r = view(Γ_r, jΓ+1:jΓ+ns)

                vdΓ = (vΓ_a, vΓ_b, vΓ_p, vΓ_q, vΓ_r)

                # add induced velocity from the sending surface
                if isurf == jsurf
                    # special implementation that ignores chosen bound vortex
                    Vind, dVind = induced_velocity_derivatives(I, surfaces[jsurf], vΓ, vdΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
                        wake_shedding_locations = wake_shedding_locations[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels[jsurf],
                        xhat = xhat)
                else
                    # regular induced velocity calculation
                    Vind, dVind = induced_velocity_derivatives(rc, surfaces[jsurf], vΓ, vdΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
                        wake_shedding_locations = wake_shedding_locations[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels[jsurf],
                        xhat = xhat)
                end

                Vind_a, Vind_b, Vind_p, Vind_q, Vind_r = dVind

                Vi += Vind

                Vi_a += Vind_a
                Vi_b += Vind_b
                Vi_p += Vind_p
                Vi_q += Vind_q
                Vi_r += Vind_r

                # add induced velocity from wake panels (assume wake is frozen)
                if wake_panels[jsurf]
                    Vi += induced_velocity(rc, system.wakes[jsurf];
                        symmetric = symmetric[jsurf],
                        finite_core = wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end

                jΓ += ns
            end

            # circulation of current panel (and its time derivative)
            if isone(I[1])
                Γi = Γ[iΓ+i]
                dΓdti = dΓdt[iΓ+i]

                Γi_a = Γ_a[iΓ+i]
                Γi_b = Γ_b[iΓ+i]
                Γi_p = Γ_p[iΓ+i]
                Γi_q = Γ_q[iΓ+i]
                Γi_r = Γ_r[iΓ+i]
            else
                Γi = Γ[iΓ+i] - Γ[iΓ+i-1]
                dΓdti = (dΓdt[iΓ+i] + dΓdt[iΓ+i-1])/2

                Γi_a = Γ_a[iΓ+i] - Γ_a[iΓ+i-1]
                Γi_b = Γ_b[iΓ+i] - Γ_b[iΓ+i-1]
                Γi_p = Γ_p[iΓ+i] - Γ_p[iΓ+i-1]
                Γi_q = Γ_q[iΓ+i] - Γ_q[iΓ+i-1]
                Γi_r = Γ_r[iΓ+i] - Γ_r[iΓ+i-1]
            end

            # circulation vector
            Δs = top_vector(receiving[I])

            # now use the Kutta-Joukowski theorem to solve for the panel force
            tmp = cross(Vi, Δs)

            # steady part of Kutta-Joukowski theorem
            Fbi = RHO*Γi*tmp

            Fbi_a = RHO*(Γi_a*tmp + Γi*cross(Vi_a, Δs))
            Fbi_b = RHO*(Γi_b*tmp + Γi*cross(Vi_b, Δs))
            Fbi_p = RHO*(Γi_p*tmp + Γi*cross(Vi_p, Δs))
            Fbi_q = RHO*(Γi_q*tmp + Γi*cross(Vi_q, Δs))
            Fbi_r = RHO*(Γi_r*tmp + Γi*cross(Vi_r, Δs))

            # unsteady part of Kutta-Joukowski theorem
            if unsteady
                # panel chord length
                c = receiving[I].chord

                #TODO: decide whether to divide by perpindicular velocity like
                # Drela does in ASWING?

                # unsteady part of Kutta-Joukowski theorem
                Fbi += RHO*dΓdti*c*tmp
            end

            # --- Calculate forces for the left bound vortex --- #

            rc = left_center(receiving[I])

            # get effective velocity at the bound vortex center
            Veff, dVeff = external_velocity_derivatives(fs, rc, ref.r)

            # unpack derivatives
            Veff_a, Veff_b, Veff_p, Veff_q, Veff_r = dVeff

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # now use the Kutta-Joukowski theorem to calculate the forces
            Γli = Γ[iΓ+i]

            Γli_a = Γ_a[iΓ+i]
            Γli_b = Γ_b[iΓ+i]
            Γli_p = Γ_p[iΓ+i]
            Γli_q = Γ_q[iΓ+i]
            Γli_r = Γ_r[iΓ+i]

            Δs = left_vector(receiving[I])

            tmp = cross(Veff, Δs)

            Fbli = RHO*Γli*tmp

            Fbli_a = RHO*(Γli_a*tmp + Γli*cross(Veff_a, Δs))
            Fbli_b = RHO*(Γli_b*tmp + Γli*cross(Veff_b, Δs))
            Fbli_p = RHO*(Γli_p*tmp + Γli*cross(Veff_p, Δs))
            Fbli_q = RHO*(Γli_q*tmp + Γli*cross(Veff_q, Δs))
            Fbli_r = RHO*(Γli_r*tmp + Γli*cross(Veff_r, Δs))

            # --- Calculate forces on the right bound vortex --- #

            rc = right_center(receiving[I])

            # get effective velocity at the bound vortex center
            Veff, dVeff = external_velocity_derivatives(fs, rc, ref.r)

            # unpack derivatives
            Veff_a, Veff_b, Veff_p, Veff_q, Veff_r = dVeff

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # now use the Kutta-Joukowski theorem to calculate the forces
            Γri = Γ[iΓ+i]

            Γri_a = Γ_a[iΓ+i]
            Γri_b = Γ_b[iΓ+i]
            Γri_p = Γ_p[iΓ+i]
            Γri_q = Γ_q[iΓ+i]
            Γri_r = Γ_r[iΓ+i]

            Δs = right_vector(receiving[I])

            tmp = cross(Veff, Δs)

            Fbri = RHO*Γri*tmp

            Fbri_a = RHO*(Γri_a*tmp + Γri*cross(Veff_a, Δs))
            Fbri_b = RHO*(Γri_b*tmp + Γri*cross(Veff_b, Δs))
            Fbri_p = RHO*(Γri_p*tmp + Γri*cross(Veff_p, Δs))
            Fbri_q = RHO*(Γri_q*tmp + Γri*cross(Veff_q, Δs))
            Fbri_r = RHO*(Γri_r*tmp + Γri*cross(Veff_r, Δs))

            # store panel circulation, velocity, and forces
            props[isurf][I] = PanelProperties(Γ[iΓ+i]/VINF, Vi/VINF, Fbi/(QINF*ref.S),
                Fbli/(QINF*ref.S), Fbri/(QINF*ref.S))

            props_a[isurf][I] = PanelProperties(Γ_a[iΓ+i]/VINF, Vi_a/VINF, Fbi_a/(QINF*ref.S),
                Fbli_a/(QINF*ref.S), Fbri_a/(QINF*ref.S))
            props_b[isurf][I] = PanelProperties(Γ_b[iΓ+i]/VINF, Vi_b/VINF, Fbi_b/(QINF*ref.S),
                Fbli_b/(QINF*ref.S), Fbri_b/(QINF*ref.S))
            props_p[isurf][I] = PanelProperties(Γ_p[iΓ+i]/VINF, Vi_p/VINF, Fbi_p/(QINF*ref.S),
                Fbli_p/(QINF*ref.S), Fbri_p/(QINF*ref.S))
            props_q[isurf][I] = PanelProperties(Γ_q[iΓ+i]/VINF, Vi_q/VINF, Fbi_q/(QINF*ref.S),
                Fbli_q/(QINF*ref.S), Fbri_q/(QINF*ref.S))
            props_r[isurf][I] = PanelProperties(Γ_r[iΓ+i]/VINF, Vi_r/VINF, Fbi_r/(QINF*ref.S),
                Fbli_r/(QINF*ref.S), Fbri_r/(QINF*ref.S))
        end

        # increment Γ index for receiving panels
        iΓ += nr
    end

    return system
end

"""
    body_forces(system, surface, reference, freestream; kwargs...)

Return the body force coefficients for `system`

**Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.**

# Arguments:
 - `system`: Object of type [`System`](@ref) which holds system properties
 - `surface`: Matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: reference parameters (see [`Reference`](@ref))
 - `freestream`: freestream parameters (see [`Freestream`](@ref))

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image (across the X-Z plane) of
    the panels in `surface` should be used when calculating induced velocities.
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)`
"""
function body_forces(system::System, surface::AbstractMatrix, ref, fs; symmetric,
    frame = Body())

    return body_forces([surface], system.panels, ref, fs; symmetric=[symmetric], frame)
end

"""
    body_forces(surface, surface_properties, reference, freestream; kwargs...)

Return the body force coefficients given the panel properties of `surface`

**Note that this function assumes that a near-field analysis has already been
performed to obtain `surface_properties`.**

# Arguments:
 - `surface_propertiess`: Surface properties, represented by matrices of panel
    properties (see [`PanelProperties`](@ref)) of shape (nc, ns) where `nc` is the number
    of chordwise panels and `ns` is the number of spanwise panels
 - `surface`: Matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: reference parameters (see [`Reference`](@ref))
 - `freestream`: freestream parameters (see [`Freestream`](@ref))

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image (across the X-Z plane) of
    the panels in `surface` should be used when calculating induced velocities.
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)`
"""
function body_forces(surface::AbstractMatrix, surface_properties, ref, fs; symmetric,
    frame = Body())

    return body_forces([surface], [surface_properties], ref, fs; symmetric=[symmetric], frame)
end

"""
    body_forces(system, surfaces, reference, freestream; kwargs...)

Return the body force coefficients for `system`

**Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.**

# Arguments
 - `system`: Object of type [`System`](@ref) which holds system properties
 - `surfaces`: Vector of surfaces, represented by matrices of surface panels
    (see [`SurfacePanel`](@ref) of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`]@ref)

# Keyword Arguments
 - `symmetric`: (required) Flag for each surface indicating whether a mirror image
   (across the X-Z plane) should be used when calculating induced velocities
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)`
"""
function body_forces(system::System, surfaces::AbstractVector{<:AbstractMatrix}, ref, fs;
    symmetric = fill(false, length(surfaces)),
    frame=Body())

    return body_forces(surfaces, system.panels,  ref, fs; symmetric, frame)
end

"""
    body_forces(surfaces, surface_properties, reference, freestream; kwargs...)

Return the body force coefficients given the panel properties for `surfaces`

Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.

# Arguments:
 - `surfaces`: Vector of surfaces, represented by matrices of surface panels
    (see [`SurfacePanel`](@ref) of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `surface_properties`: Vector of properties for each surface, represented by
    matrices of panel properties (see [`PanelProperties`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`]@ref)

# Keyword Arguments
 - `symmetric`: (required) Flag for each surface indicating whether a mirror image
   (across the X-Z plane) should be used when calculating induced velocities
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)`
"""
function body_forces(surfaces::AbstractVector, surface_properties, ref, fs; symmetric,
    frame = Body())

    # float number type
    TF = eltype(eltype(eltype(surface_properties)))

    # initialize body force coefficients
    CF = @SVector zeros(TF, 3)
    CM = @SVector zeros(TF, 3)

    # loop through all surfaces
    for isurf = 1:length(surfaces)

        # initialize surface contribution to body force coefficients
        CFi = @SVector zeros(TF, 3)
        CMi = @SVector zeros(TF, 3)

        # loop through all panels on this surface
        for i = 1:length(surfaces[isurf])

            # top bound vortex
            rc = top_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = surface_properties[isurf][i].cfb
            CFi += cf
            CMi += cross(Δr, cf)

            # left bound vortex
            rc = left_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = surface_properties[isurf][i].cfl
            CFi += cf
            CMi += cross(Δr, cf)

            # right bound vortex
            rc = right_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = surface_properties[isurf][i].cfr
            CFi += cf
            CMi += cross(Δr, cf)
        end

        # adjust forces from this surface to account for symmetry
        if symmetric[isurf]
            CFi = SVector(2*CFi[1], 0.0, 2*CFi[3])
            CMi = SVector(0.0, 2*CMi[2], 0.0)
        end

        # add to body forces
        CF += CFi
        CM += CMi

    end

    # add reference length in moment normalization
    reference_length = SVector(ref.b, ref.c, ref.b)
    CM = CM ./ reference_length

    # switch to specified frame
    CF, CM = body_to_frame(CF, CM, ref, fs, frame)

    return CF, CM
end

"""
    body_forces_derivatives(system, surface, reference, freestream; kwargs...)

Return the body force coefficients for the `system` and their derivatives with
respect to the freestream variables

Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.

# Arguments:
 - `system`: Object of type `System` which holds system properties
 - `surfaces`: Surface, represented by matrices of panels of shape (nc, ns) where
    `nc` is the number of chordwise panels and `ns` is the number of spanwise panels
 - `reference`: reference parameters (see `Reference`)
 - `freestream`: freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: Flag indicating whether a mirror image of the panels in `surface`
    should be used when calculating induced velocities.
"""
@inline function body_forces_derivatives(system::System, surface::AbstractMatrix, ref, fs;
    symmetric=false)

    return body_forces_derivatives(system, [surface], ref, fs; symmetric=[symmetric])
end

"""
    body_forces_derivatives(system, surfaces, reference, freestream; kwargs...)

Return the body force coefficients for the `system` and their derivatives with
respect to the freestream variables

Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.

# Arguments:
 - `system`: Object of type `System` which holds system properties
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: reference parameters (see `Reference`)
 - `freestream`: freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: Flags indicating whether a mirror image of each panel in `surfaces`
    should be used when calculating induced velocities.
"""
@inline function body_forces_derivatives(system::System, surfaces::AbstractVector, ref, fs;
    symmetric = fill(false, length(surfaces)))

    # float number type
    TF = eltype(system)

    # extract surface properties and their derivatives
    surface_properties = system.panels
    props_a, props_b, props_p, props_q, props_r = system.dpanels

    # initialize body force coefficients
    CF = @SVector zeros(TF, 3)
    CM = @SVector zeros(TF, 3)

    CF_a = @SVector zeros(TF, 3)
    CF_b = @SVector zeros(TF, 3)
    CF_p = @SVector zeros(TF, 3)
    CF_q = @SVector zeros(TF, 3)
    CF_r = @SVector zeros(TF, 3)

    CM_a = @SVector zeros(TF, 3)
    CM_b = @SVector zeros(TF, 3)
    CM_p = @SVector zeros(TF, 3)
    CM_q = @SVector zeros(TF, 3)
    CM_r = @SVector zeros(TF, 3)

    # loop through all surfaces
    for isurf = 1:length(surfaces)

        # initialize surface contribution to body force coefficients
        CFi = @SVector zeros(TF, 3)
        CMi = @SVector zeros(TF, 3)

        CFi_a = @SVector zeros(TF, 3)
        CFi_b = @SVector zeros(TF, 3)
        CFi_p = @SVector zeros(TF, 3)
        CFi_q = @SVector zeros(TF, 3)
        CFi_r = @SVector zeros(TF, 3)

        CMi_a = @SVector zeros(TF, 3)
        CMi_b = @SVector zeros(TF, 3)
        CMi_p = @SVector zeros(TF, 3)
        CMi_q = @SVector zeros(TF, 3)
        CMi_r = @SVector zeros(TF, 3)

        for i = 1:length(surfaces[isurf])

            # top bound vortex
            rc = top_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = surface_properties[isurf][i].cfb
            CFi += cf
            CMi += cross(Δr, cf)

            cf_a = props_a[isurf][i].cfb
            CFi_a += cf_a
            CMi_a += cross(Δr, cf_a)

            cf_b = props_b[isurf][i].cfb
            CFi_b += cf_b
            CMi_b += cross(Δr, cf_b)

            cf_p = props_p[isurf][i].cfb
            CFi_p += cf_p
            CMi_p += cross(Δr, cf_p)

            cf_q = props_q[isurf][i].cfb
            CFi_q += cf_q
            CMi_q += cross(Δr, cf_q)

            cf_r = props_r[isurf][i].cfb
            CFi_r += cf_r
            CMi_r += cross(Δr, cf_r)

            # left bound vortex
            rc = left_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = surface_properties[isurf][i].cfl
            CFi += cf
            CMi += cross(Δr, cf)

            cf_a = props_a[isurf][i].cfl
            CFi_a += cf_a
            CMi_a += cross(Δr, cf_a)

            cf_b = props_b[isurf][i].cfl
            CFi_b += cf_b
            CMi_b += cross(Δr, cf_b)

            cf_p = props_p[isurf][i].cfl
            CFi_p += cf_p
            CMi_p += cross(Δr, cf_p)

            cf_q = props_q[isurf][i].cfl
            CFi_q += cf_q
            CMi_q += cross(Δr, cf_q)

            cf_r = props_r[isurf][i].cfl
            CFi_r += cf_r
            CMi_r += cross(Δr, cf_r)

            # right bound vortex
            rc = right_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = surface_properties[isurf][i].cfr
            CFi += cf
            CMi += cross(Δr, cf)

            cf_a = props_a[isurf][i].cfr
            CFi_a += cf_a
            CMi_a += cross(Δr, cf_a)

            cf_b = props_b[isurf][i].cfr
            CFi_b += cf_b
            CMi_b += cross(Δr, cf_b)

            cf_p = props_p[isurf][i].cfr
            CFi_p += cf_p
            CMi_p += cross(Δr, cf_p)

            cf_q = props_q[isurf][i].cfr
            CFi_q += cf_q
            CMi_q += cross(Δr, cf_q)

            cf_r = props_r[isurf][i].cfr
            CFi_r += cf_r
            CMi_r += cross(Δr, cf_r)
        end

        # adjust forces from this surface to account for symmetry
        if symmetric[isurf]
            CFi = SVector(2*CFi[1], 0.0, 2*CFi[3])
            CMi = SVector(0.0, 2*CMi[2], 0.0)

            CFi_a = SVector(2*CFi_a[1], 0.0, 2*CFi_a[3])
            CFi_b = SVector(2*CFi_b[1], 0.0, 2*CFi_b[3])
            CFi_p = SVector(2*CFi_p[1], 0.0, 2*CFi_p[3])
            CFi_q = SVector(2*CFi_q[1], 0.0, 2*CFi_q[3])
            CFi_r = SVector(2*CFi_r[1], 0.0, 2*CFi_r[3])

            CMi_a = SVector(0.0, 2*CMi_a[2], 0.0)
            CMi_b = SVector(0.0, 2*CMi_b[2], 0.0)
            CMi_p = SVector(0.0, 2*CMi_p[2], 0.0)
            CMi_q = SVector(0.0, 2*CMi_q[2], 0.0)
            CMi_r = SVector(0.0, 2*CMi_r[2], 0.0)
        end

        # add to body forces
        CF += CFi
        CM += CMi

        CF_a += CFi_a
        CF_b += CFi_b
        CF_p += CFi_p
        CF_q += CFi_q
        CF_r += CFi_r

        CM_a += CMi_a
        CM_b += CMi_b
        CM_p += CMi_p
        CM_q += CMi_q
        CM_r += CMi_r

    end

    # add reference length in moment normalization
    reference_length = SVector(ref.b, ref.c, ref.b)
    CM = CM ./ reference_length

    CM_a = CM_a ./ reference_length
    CM_b = CM_b ./ reference_length
    CM_p = CM_p ./ reference_length
    CM_q = CM_q ./ reference_length
    CM_r = CM_r ./ reference_length

    # pack up derivatives
    dCF = (CF_a, CF_b, CF_p, CF_q, CF_r)
    dCM = (CM_a, CM_b, CM_p, CM_q, CM_r)

    return CF, CM, dCF, dCM
end

"""
    body_forces_history(surface, surface_history, ref, fs; symmetric, frame=Body())

Return the body force coefficients `CF`, `CM` at each time step in `surface_history`.

**Note that body forces are normalized by the instantaneous freestream velocity
rather than a common freestream velocity**

# Arguments:
 - `surface`: Matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `surface_history`: Vector of surface properties at each time step, where surface
    properties are represented by a matrix of panel properties (see [`PanelProperties`](@ref))
    of shape (nc, ns) where `nc` is the number of chordwise panels and `ns` is
    the number of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`]@ref)

# Keyword Arguments
 - `symmetric`: (required) Flag for each surface indicating whether a mirror image
   (across the X-Z plane) should be used when calculating induced velocities
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)
"""
function body_forces_history(surface::AbstractMatrix, surface_history, ref, fs;
    symmetric, frame=Body())

    # float type
    TF = eltype(eltype(eltype(surface_history)))

    # number of time steps
    nt = length(surface_history)

    # convert single freestream input to vector
    if isa(fs, Freestream)
        fs = fill(fs, nt)
    end

    # initialize time history coefficients
    CF = Vector{SVector{3, TF}}(undef, nt)
    CM = Vector{SVector{3, TF}}(undef, nt)

    # populate time history coefficients
    for it = 1:nt
        CF[it], CM[it] = body_forces(surface, surface_history[it], ref, fs[it]; symmetric, frame)
    end

    return CF, CM
end

"""
    body_forces_history(surfaces, surface_history, ref, fs; symmetric, frame=Body())

Return the body force coefficients `CF`, `CM` at each time step in `surface_history`.

**Note that body forces are normalized by the instantaneous freestream velocity
rather than a common freestream velocity**

# Arguments:
 - `surface`: Matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `surface_history`: Vector of surface properties for each surface at each
    time step, where surface properties are represented by a matrix of panel
    properties (see [`PanelProperties`](@ref)) of shape (nc, ns) where `nc` is
    the number of chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`]@ref)

# Keyword Arguments
 - `symmetric`: (required) Flag for each surface indicating whether a mirror image
   (across the X-Z plane) should be used when calculating induced velocities
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)`
"""
function body_forces_history(surfaces::AbstractVector{<:AbstractMatrix},
    surface_history, ref, fs; symmetric, frame=Body())

    # float type
    TF = eltype(eltype(eltype(eltype(surface_history))))

    # number of time steps
    nt = length(surface_history)

    # convert single freestream input to vector
    if isa(fs, Freestream)
        fs = fill(fs, nt)
    end

    # initialize time history coefficients
    CF = Vector{SVector{3, TF}}(undef, nt)
    CM = Vector{SVector{3, TF}}(undef, nt)

    # populate time history coefficients
    for it = 1:nt
        CF[it], CM[it] = body_forces(surfaces, surface_history[it], ref, fs[it]; symmetric, frame)
    end

    return CF, CM
end

"""
    body_to_frame(CF, CM, reference, freestream, frame)

Transform the coefficients `CF` and `CM` from the body frame to the frame
specified in `frame`
"""
body_to_frame

@inline body_to_frame(CF, CM, ref, fs, ::Body) = CF, CM

@inline function body_to_frame(CF, CM, ref, fs, ::Stability)
    R = body_to_stability(fs)
    return R*CF, R*CM
end

@inline function body_to_frame(CF, CM, ref, fs, ::Wind)
    # remove reference lengths
    reflen = SVector(ref.b, ref.c, ref.b)
    CM = CM .* reflen

    # rotate
    R = body_to_wind(fs)
    CF = R*CF
    CM = R*CM

    # reapply reference lengths
    CM = CM ./ reflen

    return CF, CM
end
