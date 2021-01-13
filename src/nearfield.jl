"""
    near_field_forces!(system, surface, ref, fs; kwargs...)

Perform a near-field analysis to calculate panel forces in the body frame.  Return
the modified `system` containing the calculated panel forces.

Note that this function assumes that the circulation distribution has already
been calculated and stored in `system`.

# Arguments:
 - `system`: Object of type `System` which holds system properties
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: reference parameters (see `Reference`)
 - `freestream`: freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: Flag indicating whether a mirror image of the panels in `surface`
    should be used when calculating induced velocities.
 - `trailing_vortices`: Flag indicating whether trailing vortices should be shed
    from the last chordwise panel in `wake`
 - `xhat`: direction in which trailing vortices are shed
 - `nwake`: number of chordwise wake panels to use from the wake stored in `system`.
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
"""
@inline function near_field_forces!(system, surface::AbstractMatrix, ref, fs;
    symmetric = false,
    nwake = size(system.wakes[1], 1),
    trailing_vortices = true,
    surface_id = 1,
    wake_id = -1,
    kwargs...)

    return near_field_forces!(system, [surface], ref, fs; symmetric=[symmetric],
        nwake = [nwake], trailing_vortices = [trailing_vortices],
        surface_id = surface_id:surface_id, wake_id = wake_id:wake_id, kwargs...)
end

"""
    near_field_forces!(system, surfaces, reference, freestream; kwargs...)

Perform a near-field analysis to calculate panel forces in the body frame.  Return
the modified `system` containing the calculated panel forces.

Note that this function assumes that the circulation distribution has already
been calculated and stored in `system`.

# Arguments:
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: reference parameters (see `Reference`)
 - `freestream`: freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: Flags indicating whether a mirror image (across the X-Z plane) should
    be used when calculating induced velocities
 - `trailing_vortices`: Flags to enable/disable trailing vortices
    for each surface
 - `xhat`: Direction in which to shed trailing vortices
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
"""
@inline function near_field_forces!(system, surfaces::AbstractVector{<:AbstractMatrix},
    ref, fs;
    symmetric = false,
    trailing_vortices = fill(false, length(surfaces)),
    xhat = SVector(1, 0, 0),
    nwake = size.(system.wakes, 1),
    surface_id = 1:length(surfaces),
    wake_id = -1:-1:length(surfaces))

    # unpack system storage
    Γ = system.Γ
    dΓdt = system.dΓdt

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

            # --- Calculate forces on the panel bound vortex -- #
            I = cr[i]

            rc = top_center(receiving[I])

            # start with external velocity at the bound vortex center
            Vi = external_velocity(fs, rc, ref.r)

            # now add the velocity induced by each sending panel
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                sending = surfaces[jsurf]
                ns = length(sending)

                # find out whether surfaces have the same index and/or ID
                vΓ = view(Γ, jΓ+1:jΓ+ns)

                if isurf == jsurf
                    # special implementation that ignores chosen bound vortex
                    Vi += induced_velocity(I, surfaces[jsurf], vΓ;
                        finite_core = surface_id[jsurf] < 0 || surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels[jsurf],
                        xhat=xhat)
                else
                    # normal induced velocity calculation
                    Vi += induced_velocity(rc, surfaces[jsurf], vΓ;
                        finite_core = surface_id[jsurf] < 0 || surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels[jsurf],
                        xhat=xhat)
                end

                if wake_panels[jsurf]
                    Vi += induced_velocity(rc, system.wakes[jsurf];
                        symmetric = symmetric[jsurf],
                        finite_core = wake_id[jsurf] < 0 || surface_id[isurf] != wake_id[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end

                jΓ += ns # increment Γ index for sending panels
            end

            # circulation of current panel
            if isone(I[1])
                Γi = Γ[iΓ+i]
            else
                Γi = Γ[iΓ+i] - Γ[iΓ+i-1]
            end

            # now use the Kutta-Joukowski theorem to solve for the panel force
            Δs = top_vector(receiving[I])
            Fbi = RHO*Γi*cross(Vi, Δs)

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

            # extract circulation derivative
            dΓdti = dΓdt[iΓ+i]
            # dimensionalize non-dimensional time
            dΓdti *= ref.c/VINF
            # unsteady forces are in normal direction
            ncp = normal(receiving[I])
            # unsteady forces act over the panel area
            area = receiving[I].area
            # calculate unsteady forces
            Fui = RHO*dΓdti*area*ncp

            # store panel circulation, velocity, and forces
            system.panels[isurf][I] = PanelProperties(Γ[iΓ+i]/VINF, Vi/VINF,
                Fbi/(QINF*ref.S), Fbli/(QINF*ref.S), Fbri/(QINF*ref.S), Fui/(QINF*ref.S))
        end
        # increment Γ index for receiving panels
        iΓ += nr
    end

    return system
end

"""
    near_field_forces_derivatives!(system, surface, ref, fs; kwargs...)

Perform a near-field analysis to calculate panel forces and their derivatives
with respect to the freestream variables in the body frame.  Return the modified
`system` containing the calculated panel forces and their derivatives.

Note that this function assumes that the circulation distribution has already
been calculated and stored in `system`.

# Arguments:
 - `system`: Object of type `System` which holds system properties
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: reference parameters (see `Reference`)
 - `freestream`: freestream parameters (see `Freestream`)

# Keyword Arguments
 - `surface_id`: Surface ID.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
    Additionally, if a surface/wake's ID is negative the finite core model will
    always be enabled, even when calculating the influence of the surface/wake
    on itself. By default each surface has its own ID.
 - `wake_id`: Wake ID.  The finite core model is disabled when
    calculating the influence of surfaces/wakes that share the same ID.
    Additionally, if a surface/wake's ID is negative, the finite core model will
    always be enabled, even when calculating the influence of the surface/wake
    on itself. By default the finite core model for the wakes is always enabled.
 - `symmetric`: Flag indicating whether a mirror image of the panels in `surface`
    should be used when calculating induced velocities.
 - `nwake`: number of chordwise wake panels to use from the wake stored in `system`.
 - `trailing_vortices`: Flag indicating whether trailing vortices should be shed
    from the last chordwise panel in `wake`
 - `xhat`: direction in which trailing vortices are shed
"""
@inline function near_field_forces_derivatives!(system, surface::AbstractMatrix, ref, fs;
    surface_id = 1,
    wake_id = -1,
    symmetric = false,
    nwake = size(system.wakes[1], 1),
    steady_wake = false,
    trailing_vortices = true,
    kwargs...)

    return near_field_forces_derivatives!(system, [surface], ref, fs;
        symmetric = [symmetric],
        surface_id = surface_id:surface_id,
        wake_id = wake_id:wake_id,
        nwake = [nwake],
        trailing_vortices = [trailing_vortices],
        kwargs...)
end

"""
    near_field_forces_derivatives!(system, surfaces, reference, freestream; kwargs...)

Perform a near-field analysis to calculate panel forces and their derivatives
with respect to the freestream variables in the body frame.  Return the modified
`system` containing the calculated panel forces and their derivatives.

Note that this function assumes that the circulation distribution has already
been calculated and stored in `system`.

# Arguments:
 - `system`: Object of type `System` which holds system properties
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: reference parameters (see `Reference`)
 - `freestream`: freestream parameters (see `Freestream`)

# Keyword Arguments
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
 - `symmetric`: Flags indicating whether a mirror image (across the X-Z plane) should
    be used when calculating induced velocities
 - `trailing_vortices`: Flags to enable/disable trailing vortices
 - `xhat`: Direction in which to shed trailing vortices
"""
@inline function near_field_forces_derivatives!(system, surfaces::AbstractVector{<:AbstractMatrix},
    ref, fs;
    surface_id = 1:length(surfaces),
    wake_id = -1:-1:-length(surfaces),
    symmetric = fill(false, length(surfaces)),
    nwake = size.(system.wakes, 1),
    trailing_vortices = fill(true, length(surfaces)),
    xhat = SVector(1, 0, 0),
    Gamma_t = nothing)

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

            # --- Calculate forces for the bound vortices --- #
            rc = top_center(receiving[I])

            # get external velocity at the bound vortex center
            Vi, dVi = external_velocity_derivatives(fs, rc, ref.r)

            # unpack derivatives
            Vi_a, Vi_b, Vi_p, Vi_q, Vi_r = dVi

            # now add the velocity induced by each sending panel
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                sending = surfaces[jsurf]
                ns = length(sending)

                vΓ = view(Γ, jΓ+1:jΓ+ns)

                vΓ_a = view(Γ_a, jΓ+1:jΓ+ns)
                vΓ_b = view(Γ_b, jΓ+1:jΓ+ns)
                vΓ_p = view(Γ_p, jΓ+1:jΓ+ns)
                vΓ_q = view(Γ_q, jΓ+1:jΓ+ns)
                vΓ_r = view(Γ_r, jΓ+1:jΓ+ns)

                vdΓ = (vΓ_a, vΓ_b, vΓ_p, vΓ_q, vΓ_r)

                if isurf == jsurf
                    # special implementation that ignores chosen bound vortex
                    Vind, dVind = induced_velocity_derivatives(I, surfaces[jsurf], vΓ, vdΓ;
                        finite_core = surface_id[jsurf] < 0 || surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels[jsurf],
                        xhat = xhat)
                else
                    # normal induced velocity calculation
                    Vind, dVind = induced_velocity_derivatives(rc, surfaces[jsurf], vΓ, vdΓ;
                        finite_core = surface_id[jsurf] < 0 || surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
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

                # assume a frozen wake
                if wake_panels[jsurf]
                    Vi += induced_velocity(rc, system.wakes[jsurf];
                        symmetric = symmetric[jsurf],
                        finite_core = wake_id[jsurf] < 0 || surface_id[isurf] != wake_id[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end

                jΓ += ns
            end

            # circulation of current panel
            if isone(I[1])
                Γi = Γ[iΓ+i]
                Γi_a = Γ_a[iΓ+i]
                Γi_b = Γ_b[iΓ+i]
                Γi_p = Γ_p[iΓ+i]
                Γi_q = Γ_q[iΓ+i]
                Γi_r = Γ_r[iΓ+i]
            else
                Γi = Γ[iΓ+i] - Γ[iΓ+i-1]
                Γi_a = Γ_a[iΓ+i] - Γ_a[iΓ+i-1]
                Γi_b = Γ_b[iΓ+i] - Γ_b[iΓ+i-1]
                Γi_p = Γ_p[iΓ+i] - Γ_p[iΓ+i-1]
                Γi_q = Γ_q[iΓ+i] - Γ_q[iΓ+i-1]
                Γi_r = Γ_r[iΓ+i] - Γ_r[iΓ+i-1]
            end

            # now use the Kutta-Joukowski theorem to calculate the forces at this point
            Δs = top_vector(receiving[I])
            tmp = cross(Vi, Δs)
            Fbi = RHO*Γi*tmp

            Fbi_a = RHO*(Γi_a*tmp + Γi*cross(Vi_a, Δs))
            Fbi_b = RHO*(Γi_b*tmp + Γi*cross(Vi_b, Δs))
            Fbi_p = RHO*(Γi_p*tmp + Γi*cross(Vi_p, Δs))
            Fbi_q = RHO*(Γi_q*tmp + Γi*cross(Vi_q, Δs))
            Fbi_r = RHO*(Γi_r*tmp + Γi*cross(Vi_r, Δs))

            # --- Calculate forces for the left bound vortex --- #

            rc = left_center(receiving[I])

            # get effective velocity at the bound vortex center
            Veff, dVeff = external_velocity_derivatives(fs, rc, ref.r)

            # unpack derivatives
            Veff_a, Veff_b, Veff_p, Veff_q, Veff_r = dVeff

            # note that we don't include induced velocity in the effective velocity
            # because its influence is generally negligible once we take the cross product
            # with the bound vortex vector, this is also assumed in AVL

            # now use the Kutta-Joukowski theorem to calculate the forces
            Γli = Γ[iΓ+i]
            Δs = left_vector(receiving[I])
            tmp = cross(Veff, Δs)
            Fbli = RHO*Γli*tmp

            Γli_a = Γ_a[iΓ+i]
            Γli_b = Γ_b[iΓ+i]
            Γli_p = Γ_p[iΓ+i]
            Γli_q = Γ_q[iΓ+i]
            Γli_r = Γ_r[iΓ+i]

            Fbli_a = RHO*(Γli_a*tmp + Γli*cross(Veff_a, Δs))
            Fbli_b = RHO*(Γli_b*tmp + Γli*cross(Veff_b, Δs))
            Fbli_p = RHO*(Γli_p*tmp + Γli*cross(Veff_p, Δs))
            Fbli_q = RHO*(Γli_q*tmp + Γli*cross(Veff_q, Δs))
            Fbli_r = RHO*(Γli_r*tmp + Γli*cross(Veff_r, Δs))

            # --- Calculate forces for the right bound vortex --- #

            rc = right_center(receiving[I])

            # get effective velocity at the bound vortex center
            Veff, dVeff = external_velocity_derivatives(fs, rc, ref.r)

            # unpack derivatives
            Veff_a, Veff_b, Veff_p, Veff_q, Veff_r = dVeff

            # note that we don't include induced velocity in the effective velocity
            # because its influence is generally negligible once we take the cross product
            # with the bound vortex vector, this is also assumed in AVL

            # now use the Kutta-Joukowski theorem to calculate the forces
            Γri = Γ[iΓ+i]
            Δs = right_vector(receiving[I])
            tmp = cross(Veff, Δs)
            Fbri = RHO*Γri*tmp

            Γri_a = Γ_a[iΓ+i]
            Γri_b = Γ_b[iΓ+i]
            Γri_p = Γ_p[iΓ+i]
            Γri_q = Γ_q[iΓ+i]
            Γri_r = Γ_r[iΓ+i]

            Fbri_a = RHO*(Γri_a*tmp + Γri*cross(Veff_a, Δs))
            Fbri_b = RHO*(Γri_b*tmp + Γri*cross(Veff_b, Δs))
            Fbri_p = RHO*(Γri_p*tmp + Γri*cross(Veff_p, Δs))
            Fbri_q = RHO*(Γri_q*tmp + Γri*cross(Veff_q, Δs))
            Fbri_r = RHO*(Γri_r*tmp + Γri*cross(Veff_r, Δs))

            # --- Calculate Unsteady Forces on the Panel --- #
            # extract circulation derivative
            dΓdti = dΓdt[iΓ+i]
            # dimensionalize non-dimensional time
            dΓdti *= ref.c/VINF
            # unsteady forces are in normal direction
            ncp = normal(receiving[I])
            # unsteady forces act over the panel area
            area = receiving[I].area
            # calculate unsteady forces
            Fui = RHO*dΓdti*area*ncp

            # assemble normalized panel outputs
            props[isurf][I] = PanelProperties(Γ[iΓ+i]/VINF, Vi/VINF, Fbi/(QINF*ref.S),
                Fbli/(QINF*ref.S), Fbri/(QINF*ref.S), (@SVector zeros(3)))

            # also assemble normalized derivative panel outputs
            props_a[isurf][I] = PanelProperties(Γ_a[iΓ+i]/VINF, Vi_a/VINF, Fbi_a/(QINF*ref.S),
                Fbli_a/(QINF*ref.S), Fbri_a/(QINF*ref.S), (@SVector zeros(3)))
            props_b[isurf][I] = PanelProperties(Γ_b[iΓ+i]/VINF, Vi_b/VINF, Fbi_b/(QINF*ref.S),
                Fbli_b/(QINF*ref.S), Fbri_b/(QINF*ref.S), (@SVector zeros(3)))
            props_p[isurf][I] = PanelProperties(Γ_p[iΓ+i]/VINF, Vi_p/VINF, Fbi_p/(QINF*ref.S),
                Fbli_p/(QINF*ref.S), Fbri_p/(QINF*ref.S), (@SVector zeros(3)))
            props_q[isurf][I] = PanelProperties(Γ_q[iΓ+i]/VINF, Vi_q/VINF, Fbi_q/(QINF*ref.S),
                Fbli_q/(QINF*ref.S), Fbri_q/(QINF*ref.S), (@SVector zeros(3)))
            props_r[isurf][I] = PanelProperties(Γ_r[iΓ+i]/VINF, Vi_r/VINF, Fbi_r/(QINF*ref.S),
                Fbli_r/(QINF*ref.S), Fbri_r/(QINF*ref.S), (@SVector zeros(3)))
        end
        # increment Γ index
        iΓ += nr
    end

    return system
end

"""
    body_forces(system, surface, reference, freestream; kwargs...)

Return the body force coefficients for the `system`

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
 - `frame`: frame in which to return `CF` and `CM`, options are `Body()` (default),
   `Stability()`, and `Wind()`
"""
function body_forces(system::System, surface::AbstractMatrix, ref, fs;
    symmetric=false,
    frame=Body())

    return body_forces(system.panels, [surface], ref, fs; symmetric=[symmetric], frame)
end

"""
    body_forces(surface_properties, surface, reference, freestream; kwargs...)

Return the body force coefficients given the panel properties for `surface`

Note that this function assumes that a near-field analysis has already been
performed to obtain `surface_properties`.

# Arguments:
 - `surface_propertiess`: Surface properties, represented by matrices of panel
    properties (see `PanelProperties`) of shape (nc, ns) where `nc` is the number
    of chordwise panels and `ns` is the number of spanwise panels
 - `surface`: Surface, represented by matrices of panels of shape (nc, ns) where
    `nc` is the number of chordwise panels and `ns` is the number of spanwise panels
 - `reference`: reference parameters (see `Reference`)
 - `freestream`: freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: Flag indicating whether a mirror image of the panels in `surface`
    should be used when calculating induced velocities.
 - `frame`: frame in which to return `CF` and `CM`, options are `Body()` (default),
   `Stability()`, and `Wind()`
"""
function body_forces(surface_properties, surface::AbstractMatrix, ref, fs;
    symmetric=false,
    frame=Body())

    return body_forces([surface_properties], [surface], ref, fs; symmetric=[symmetric], frame)
end

"""
    body_forces(system, surfaces, reference, freestream; kwargs...)

Return the body force coefficients for the `system` in the given `frame`.

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
 - `frame`: frame in which to return `CF` and `CM`, options are `Body()` (default),
   `Stability()`, and `Wind()`
"""
function body_forces(system::System, surfaces::AbstractVector{<:AbstractMatrix}, ref, fs;
    symmetric = fill(false, length(surfaces)),
    frame=Body())

    return body_forces(system.panels, surfaces, ref, fs; symmetric=symmetric, frame)
end

"""
    body_forces(surface_properties, surfaces, reference, freestream; kwargs...)

Return the body force coefficients given the panel properties for `surfaces`

Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.

# Arguments:
 - `surface_propertiess`: Vector of properties for each surface, represented by
    matrices of panel properties (see `PanelProperties`) of shape (nc, ns) where
    `nc` is the number of chordwise panels and `ns` is the number of spanwise panels
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: reference parameters (see `Reference`)
 - `freestream`: freestream parameters (see `Freestream`)

 # Keyword Arguments
  - `symmetric`: Flags indicating whether a mirror image of each panel in `surfaces`
     should be used when calculating induced velocities.
 - `frame`: frame in which to return `CF` and `CM`, options are `Body()` (default),
   `Stability()`, and `Wind()`
"""
function body_forces(surface_properties, surfaces::AbstractVector, ref, fs;
    symmetric = fill(false, length(surfaces)),
    frame=Body())

    TF = eltype(eltype(eltype(surface_properties)))

    CF = @SVector zeros(TF, 3)
    CM = @SVector zeros(TF, 3)

    for isurf = 1:length(surfaces)

        CFi = @SVector zeros(TF, 3)
        CMi = @SVector zeros(TF, 3)

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

            # unsteady forces
            rc = controlpoint(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = surface_properties[isurf][i].cfu
            CFi += cf
            CMi += cross(Δr, cf)
        end

        # adjust surface forces to account for symmetry
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

    TF = eltype(system)

    props = system.panels
    props_a, props_b, props_p, props_q, props_r = system.dpanels

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

    for isurf = 1:length(surfaces)

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
            cf = props[isurf][i].cfb
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
            cf = props[isurf][i].cfl
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
            cf = props[isurf][i].cfr
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

            # unsteady forces
            rc = controlpoint(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = props[isurf][i].cfu
            CFi += cf
            CMi += cross(Δr, cf)

            cf_a = props_a[isurf][i].cfu
            CFi_a += cf_a
            CMi_a += cross(Δr, cf_a)

            cf_b = props_b[isurf][i].cfu
            CFi_b += cf_b
            CMi_b += cross(Δr, cf_b)

            cf_p = props_p[isurf][i].cfu
            CFi_p += cf_p
            CMi_p += cross(Δr, cf_p)

            cf_q = props_q[isurf][i].cfu
            CFi_q += cf_q
            CMi_q += cross(Δr, cf_q)

            cf_r = props_r[isurf][i].cfu
            CFi_r += cf_r
            CMi_r += cross(Δr, cf_r)
        end

        # adjust surface forces to account for symmetry
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

    dCF = (CF_a, CF_b, CF_p, CF_q, CF_r)
    dCM = (CM_a, CM_b, CM_p, CM_q, CM_r)

    return CF, CM, dCF, dCM
end

function body_forces_over_time(surface_properties, surfaces, ref, fs; symmetric, frame=Body())

    TF = eltype(eltype(eltype(surface_properties)))

    nt = length(surface_properties)

    CF = zeros(TF, 3, nt)
    CM = zeros(TF, 3, nt)

    for it = 1:length(surface_properties)
        CF[:,it], CM[:,it] = body_forces(surface_properties[it], surfaces, ref, fs; symmetric, frame)
    end

    return CF, CM
end

"""
    body_to_frame(CF, CM, ref, fs, frame)

Transforms the coefficients `CF` and `CM` which are provided in the body frame
to another frame
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
