"""
    near_field_forces!(properties, surfaces, wakes, reference, freestream, Γ;
        dΓdt, additional_velocity, Vh, Vv, symmetric, nwake, surface_id,
        wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)

Calculate local panel forces in the body frame.
"""
@inline function near_field_forces!(props, surfaces, wakes, ref, fs, Γ;
    dΓdt, additional_velocity, Vh, Vv, symmetric, nwake, surface_id,
    wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)

    # number of surfaces
    nsurf = length(surfaces)

    # loop through receiving surfaces
    iΓ = 0 # index for accessing Γ
    for isurf = 1:nsurf

        receiving = surfaces[isurf]
        nr = length(receiving)
        nr1, nr2 = size(receiving)
        cr = CartesianIndices(receiving)

        # loop through receiving panels
        for i = 1:length(receiving)

            # get panel cartesian index
            I = cr[i]

            # --- Calculate forces on the panel bound vortex --- #

            # bound vortex location
            rc = top_center(receiving[i])

            # freestream velocity
            Vi = freestream_velocity(fs)

            # rotational velocity
            Vi += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                Vi += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vh)
                Vi += Vh[isurf][i]
            end

            # induced velocity from surfaces and wakes
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                # number of panels on sending surface
                sending = surfaces[jsurf]
                Ns = length(sending)

                # see if wake panels are being used
                wake_panels = nwake[jsurf] > 0

                # check if we need to shift shedding locations
                if isnothing(wake_shedding_locations)
                    shedding_locations = nothing
                else
                    shedding_locations = wake_shedding_locations[jsurf]
                end

                # extract circulation values corresonding to the sending surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)

                # induced velocity from this surface
                if isurf == jsurf
                    # induced velocity on self
                    Vi += induced_velocity(I, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)
                else
                    # induced velocity on another surface
                    Vi += induced_velocity(rc, surfaces[jsurf], vΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)
                end

                # induced velocity from corresponding wake
                if wake_panels
                    Vi += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end

                jΓ += Ns # increment Γ index for sending panels
            end

            # steady part of Kutta-Joukowski theorem
            Γi = I[1] == 1 ? Γ[iΓ+i] : Γ[iΓ+i] - Γ[iΓ+i-1] # net circulation
            Δs = top_vector(receiving[I]) # bound vortex vector
            tmp = cross(Vi, Δs)
            Fbi = RHO*Γi*tmp

            if !isnothing(dΓdt)
                # unsteady part of Kutta-Joukowski theorem

                #TODO: decide whether to divide by perpindicular velocity like
                # Drela does in ASWING?

                dΓdti = I[1] == 1 ? dΓdt[iΓ+i] : (dΓdt[iΓ+i] + dΓdt[iΓ+i-1])/2
                c = receiving[I].chord
                Fbi += RHO*dΓdti*c*tmp

            end

            # --- Calculate forces on the left bound vortex --- #

            # bound vortex location
            rc = left_center(receiving[I])

            # freestream velocity
            Veff = freestream_velocity(fs)

            # rotational velocity
            Veff += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                Veff += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Veff += Vv[isurf][I[1], I[2]]
            end

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # steady part of Kutta-Joukowski theorem
            Γli = Γ[iΓ+i]
            Δs = left_vector(receiving[I])
            Fbli = RHO*Γli*cross(Veff, Δs)

            # --- Calculate forces on the right bound vortex --- #

            rc = right_center(receiving[I])

            # freestream velocity
            Veff = freestream_velocity(fs)

            # rotational velocity
            Veff += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                Veff += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Veff += Vv[isurf][I[1], I[2]+1]
            end

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # steady part of Kutta-Joukowski theorem
            Γri = Γ[iΓ+i]
            Δs = right_vector(receiving[I])
            Fbri = RHO*Γri*cross(Veff, Δs)

            # store panel circulation, velocity, and forces
            q = 1/2*RHO*ref.V^2

            props[isurf][i] = PanelProperties(Γ[iΓ+i]/ref.V, Vi/ref.V,
                Fbi/(q*ref.S), Fbli/(q*ref.S), Fbri/(q*ref.S))
        end

        # increment Γ index for receiving panels
        iΓ += nr
    end

    return props
end

"""
    near_field_forces_derivatives!(properties, dproperties, surfaces, reference,
        freestream, Γ, dΓ; dΓdt, additional_velocity, Vh, Vv, symmetric, nwake,
        surface_id, wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)

Version of [`near_field_forces!`](@ref) that also calculates the derivatives of
the local panel forces with respect to the freestream variables.
"""
near_field_forces_derivatives!

@inline function near_field_forces_derivatives!(props, dprops, surfaces, wakes,
    ref, fs, Γ, dΓ; dΓdt, additional_velocity, Vh, Vv, symmetric, nwake,
    surface_id, wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)

    # unpack derivatives
    props_a, props_b, props_p, props_q, props_r = dprops
    Γ_a, Γ_b, Γ_p, Γ_q, Γ_r = dΓ

    # number of surfaces
    nsurf = length(surfaces)

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
            rc = top_center(receiving[i])

            # freestream velocity
            Vi, dVi = freestream_velocity_derivatives(fs)
            Vi_a, Vi_b = dVi

            # rotational velocity
            Vrot, dVrot = rotational_velocity_derivatives(rc, fs, ref)
            Vi += Vrot
            Vi_p, Vi_q, Vi_r = dVrot

            # additional velocity field
            if !isnothing(additional_velocity)
                Vi += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vh)
                Vi += Vh[isurf][i]
            end

            # induced velocity from surfaces and wakes
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                # number of panels on sending surface
                sending = surfaces[jsurf]
                Ns = length(sending)

                # see if wake panels are being used
                wake_panels = nwake[jsurf] > 0

                # check if we need to shift shedding locations
                if isnothing(wake_shedding_locations)
                    shedding_locations = nothing
                else
                    shedding_locations = wake_shedding_locations[jsurf]
                end

                # extract circulation values corresonding to the sending surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)

                vΓ_a = view(Γ_a, jΓ+1:jΓ+Ns)
                vΓ_b = view(Γ_b, jΓ+1:jΓ+Ns)
                vΓ_p = view(Γ_p, jΓ+1:jΓ+Ns)
                vΓ_q = view(Γ_q, jΓ+1:jΓ+Ns)
                vΓ_r = view(Γ_r, jΓ+1:jΓ+Ns)

                vdΓ = (vΓ_a, vΓ_b, vΓ_p, vΓ_q, vΓ_r)

                # induced velocity from this surface
                if isurf == jsurf
                    # induced velocity on self
                    Vind, dVind = induced_velocity_derivatives(I, surfaces[jsurf], vΓ, vdΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)
                else
                    # induced velocity on another surface
                    Vind, dVind = induced_velocity_derivatives(rc, surfaces[jsurf], vΓ, vdΓ;
                        finite_core = surface_id[isurf] != surface_id[jsurf],
                        wake_shedding_locations = shedding_locations,
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                        xhat = xhat)
                end

                Vind_a, Vind_b, Vind_p, Vind_q, Vind_r = dVind

                Vi += Vind

                Vi_a += Vind_a
                Vi_b += Vind_b
                Vi_p += Vind_p
                Vi_q += Vind_q
                Vi_r += Vind_r

                # induced velocity from corresponding wake
                if wake_panels
                    Vi += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end

                jΓ += Ns # increment Γ index for sending panels
            end

            # steady part of Kutta-Joukowski theorem
            if I[1] == 1
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

            # bound vortex vector
            Δs = top_vector(receiving[I])

            tmp = cross(Vi, Δs)

            Fbi = RHO*Γi*tmp

            Fbi_a = RHO*(Γi_a*tmp + Γi*cross(Vi_a, Δs))
            Fbi_b = RHO*(Γi_b*tmp + Γi*cross(Vi_b, Δs))
            Fbi_p = RHO*(Γi_p*tmp + Γi*cross(Vi_p, Δs))
            Fbi_q = RHO*(Γi_q*tmp + Γi*cross(Vi_q, Δs))
            Fbi_r = RHO*(Γi_r*tmp + Γi*cross(Vi_r, Δs))

            if !isnothing(dΓdt)
                # unsteady part of Kutta-Joukowski theorem

                #TODO: decide whether to divide by perpindicular velocity like
                # Drela does in ASWING?

                dΓdti = I[1] == 1 ? dΓdt[iΓ+i] : (dΓdt[iΓ+i] + dΓdt[iΓ+i-1])/2
                c = receiving[I].chord
                Fbi += RHO*dΓdti*c*tmp

            end

            # --- Calculate forces for the left bound vortex --- #

            rc = left_center(receiving[I])

            # freestream velocity
            Vfs, dVfs = freestream_velocity_derivatives(fs)
            Veff = Vfs
            Veff_a, Veff_b = dVfs

            # rotational velocity
            Vrot, dVrot = rotational_velocity_derivatives(rc, fs, ref)
            Veff += Vrot
            Veff_p, Veff_q, Veff_r = dVrot

            # additional velocity field
            if !isnothing(additional_velocity)
                Veff += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Veff += Vv[isurf][I[1], I[2]]
            end

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # steady part of Kutta-Joukowski theorem
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

            # freestream velocity
            Vfs, dVfs = freestream_velocity_derivatives(fs)
            Veff = Vfs
            Veff_a, Veff_b = dVfs

            # rotational velocity
            Vrot, dVrot = rotational_velocity_derivatives(rc, fs, ref)
            Veff += Vrot
            Veff_p, Veff_q, Veff_r = dVrot

            # additional velocity field
            if !isnothing(additional_velocity)
                Veff += additional_velocity(rc)
            end

            # velocity due to surface motion
            if !isnothing(Vv)
                Veff += Vv[isurf][I[1], I[2]+1]
            end

            # NOTE: We don't include induced velocity in the effective velocity
            # for the vertical segments because its influence is likely negligible
            # once we take the cross product with the bound vortex vector. This
            # is also assumed in AVL. This could change in the future.

            # steady part of Kutta-Joukowski theorem
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
            q = 1/2*RHO*ref.V^2

            props[isurf][I] = PanelProperties(Γ[iΓ+i]/ref.V, Vi/ref.V, Fbi/(q*ref.S),
                Fbli/(q*ref.S), Fbri/(q*ref.S))

            props_a[isurf][I] = PanelProperties(Γ_a[iΓ+i]/ref.V, Vi_a/ref.V, Fbi_a/(q*ref.S),
                Fbli_a/(q*ref.S), Fbri_a/(q*ref.S))
            props_b[isurf][I] = PanelProperties(Γ_b[iΓ+i]/ref.V, Vi_b/ref.V, Fbi_b/(q*ref.S),
                Fbli_b/(q*ref.S), Fbri_b/(q*ref.S))
            props_p[isurf][I] = PanelProperties(Γ_p[iΓ+i]/ref.V, Vi_p/ref.V, Fbi_p/(q*ref.S),
                Fbli_p/(q*ref.S), Fbri_p/(q*ref.S))
            props_q[isurf][I] = PanelProperties(Γ_q[iΓ+i]/ref.V, Vi_q/ref.V, Fbi_q/(q*ref.S),
                Fbli_q/(q*ref.S), Fbri_q/(q*ref.S))
            props_r[isurf][I] = PanelProperties(Γ_r[iΓ+i]/ref.V, Vi_r/ref.V, Fbi_r/(q*ref.S),
                Fbli_r/(q*ref.S), Fbri_r/(q*ref.S))
        end

        # increment Γ index for receiving panels
        iΓ += nr
    end

    return props, dprops
end

"""
    body_forces(system; kwargs...)

Return the body force coefficients given the panel properties for `surfaces`

Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.

# Arguments
 - `system`: Object of type [`System`](@ref) which holds system properties

# Keyword Arguments
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)`
"""
function body_forces(system::System{TF}; frame = Body()) where TF

    @assert system.near_field_analysis[] "Near field analysis required"

    # unpack parameters stored in `system`
    surfaces = system.surfaces # surface panels defining each surface
    properties = system.properties # panel properties
    ref = system.reference[] # reference parameters
    fs = system.freestream[] # freestream parameters
    symmetric = system.symmetric # symmetric flag for each surface

    return body_forces(surfaces, properties, ref, fs, symmetric, frame)
end

"""
    body_forces(surfaces, properties, reference, freestream, symmetric; kwargs...)

Return the body force coefficients given the panel properties for `surfaces`

Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.

# Arguments:
 - `surfaces`: Collection of surfaces, where each surface is represented by a
    matrix of surface panels (see [`SurfacePanel`](@ref)) of shape (nc, ns)
    where `nc` is the number of chordwise panels and `ns` is the number of
    spanwise panels
 - `properties`: Surface properties for each surface, where surface
    properties for each surface are represented by a matrix of panel properties
    (see [`PanelProperties`](@ref)) of shape (nc, ns) where `nc` is the number
    of chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`]@ref)
 - `symmetric`: (required) Flag for each surface indicating whether a mirror image
   (across the X-Z plane) was used when calculating induced velocities
 - `frame`: frame in which to return `CF` and `CM`, options are [`Body()`](@ref) (default),
   [`Stability()`](@ref), and [`Wind()`](@ref)
"""
function body_forces(surfaces, properties, ref, fs, symmetric, frame)

    TF = eltype(eltype(eltype(properties)))

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
            cf = properties[isurf][i].cfb
            CFi += cf
            CMi += cross(Δr, cf)

            # left bound vortex
            rc = left_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfl
            CFi += cf
            CMi += cross(Δr, cf)

            # right bound vortex
            rc = right_center(surfaces[isurf][i])
            Δr = rc - ref.r
            cf = properties[isurf][i].cfr
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
    body_forces_derivatives(system)

Return the body force coefficients for the `system` and their derivatives with
respect to the freestream variables

Note that this function assumes that a near-field analysis has already been
performed to obtain the panel forces.

# Arguments:
 - `system`: Object of type `System` which holds system properties
"""
@inline function body_forces_derivatives(system::System)

    # float number type
    TF = eltype(system)

    @assert system.near_field_analysis[] "Near field analysis required"
    @assert system.derivatives[] "Derivative computations required"

    # unpack parameters stored in `system`
    surfaces = system.surfaces # surface panels defining each surface
    ref = system.reference[] # reference parameters
    fs = system.freestream[] # freestream parameters
    symmetric = system.symmetric # symmetric flag for each surface
    properties = system.properties
    props_a, props_b, props_p, props_q, props_r = system.dproperties

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
            cf = properties[isurf][i].cfb
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
            cf = properties[isurf][i].cfl
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
            cf = properties[isurf][i].cfr
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
    body_forces_history(system, surface_history, property_history; frame=Body())

Return the body force coefficients `CF`, `CM` at each time step in `property_history`.

# Arguments:
 - `system`: Object of type [`System`](@ref) which holds system properties
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
function body_forces_history(system, surface_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}},
    property_history; frame=Body())

    # unpack system parameters
    symmetric = system.symmetric
    ref = system.reference[]
    fs = system.freestream[]

    # float type
    TF = eltype(system)

    # number of time steps
    nt = length(property_history)

    # convert single freestream input to vector
    if isa(fs, Freestream)
        fs = fill(fs, nt)
    end

    # initialize time history coefficients
    CF = Vector{SVector{3, TF}}(undef, nt)
    CM = Vector{SVector{3, TF}}(undef, nt)

    # populate time history coefficients
    for it = 1:nt
        CF[it], CM[it] = body_forces(surface_history[it], property_history[it],
            ref, fs[it], symmetric, frame)
    end

    return CF, CM
end

"""
    lifting_line_coefficients(system, r, c; frame=Body())

Return the force and moment coefficients (per unit span) for each spanwise segment
of a lifting line representation of the geometry.

This function requires that a near-field analysis has been performed on `system`
to obtain panel forces.

# Arguments
 - `system`: Object of type [`System`](@ref) that holds precalculated
    system properties.
 - `r`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns+1) which contains the x, y, and z coordinates
    of the resulting lifting line coordinates
 - `c`: Vector with length equal to the number of surfaces, with each element
    being a vector of length `ns+1` which contains the chord lengths at each
    lifting line coordinate.

# Keyword Arguments
 - `frame`: frame in which to return `cf` and `cm`, possible options are
    [`Body()`](@ref) (default), [`Stability()`](@ref), and [`Wind()`](@ref)`

# Return Arguments:
 - `cf`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns) which contains the x, y, and z direction
    force coefficients (per unit span) for each spanwise segment.
 - `cm`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns) which contains the x, y, and z direction
    moment coefficients (per unit span) for each spanwise segment.
"""
function lifting_line_coefficients(system, r, c; frame=Body())
    TF = promote_type(eltype(system), eltype(eltype(r)), eltype(eltype(c)))
    nsurf = length(system.surfaces)
    cf = Vector{Matrix{TF}}(undef, nsurf)
    cm = Vector{Matrix{TF}}(undef, nsurf)
    for isurf = 1:nsurf
        ns = size(system.surfaces[isurf], 2)
        cf[isurf] = Matrix{TF}(undef, 3, ns)
        cm[isurf] = Matrix{TF}(undef, 3, ns)
    end
    return lifting_line_coefficients!(cf, cm, system, r, c; frame)
end

"""
    lifting_line_coefficients!(cf, cm, system, r, c; frame=Body())

In-place version of [`lifting_line_coefficients`](@ref)
"""
function lifting_line_coefficients!(cf, cm, system, r, c; frame=Body())

    # number of surfaces
    nsurf = length(system.surfaces)

    # check that a near field analysis has been performed
    @assert system.near_field_analysis[] "Near field analysis required"

    # extract reference properties
    ref = system.reference[]
    fs = system.freestream[]

    # iterate through each lifting surface
    for isurf = 1:nsurf
        nc, ns = size(system.surfaces[isurf])
        # extract current surface panels and panel properties
        panels = system.surfaces[isurf]
        properties = system.properties[isurf]
        # loop through each chordwise set of panels
        for j = 1:ns
            # calculate segment length
            rls = SVector(r[isurf][1,j], r[isurf][2,j], r[isurf][3,j])
            rrs = SVector(r[isurf][1,j+1], r[isurf][2,j+1], r[isurf][3,j+1])
            ds = norm(rrs - rls)
            # calculate reference location
            rs = (rls + rrs)/2
            # calculate reference chord
            cs = (c[isurf][j] + c[isurf][j+1])/2
            # calculate section force and moment coefficients
            cfj = @SVector zeros(eltype(cf[isurf]), 3)
            cmj = @SVector zeros(eltype(cm[isurf]), 3)
            for i = 1:nc
                # add influence of bound vortex
                rb = top_center(panels[i,j])
                cfb = properties[i,j].cfb
                cfj += cfb
                cmj += cross(rb - rs, cfb)
                # add influence of left vortex leg
                rl = left_center(panels[i,j])
                cfl = properties[i,j].cfl
                cfj += cfl
                cmj += cross(rl - rs, cfl)
                # add influence of right vortex leg
                rr = right_center(panels[i,j])
                cfr = properties[i,j].cfr
                cfj += cfr
                cmj += cross(rr - rs, cfr)
            end
            # update normalization
            cfj *= ref.S/(ds*cs)
            cmj *= ref.S/(ds*cs^2)
            # change coordinate frame
            cfj, cmj = body_to_frame(cfj, cmj, ref, fs, frame)
            # save coefficients
            cf[isurf][:,j] = cfj
            cm[isurf][:,j] = cmj
        end
    end

    return cf, cm
end

"""
    Convenience function for strip_theory_drag!()
"""
function strip_theory_drag(cf, afs, Res)
    @assert system.near_field_analysis[] "Near field analysis required"
    TF = promote_type(eltype(system), eltype(eltype(r)), eltype(eltype(c)))
    nsurf = length(system.surfaces)
    cf = Vector{Matrix{TF}}(undef, nsurf)
    cm = Vector{Matrix{TF}}(undef, nsurf)
    return strip_theory_drag!(cf, afs, Res)
end

"""
    strip_theory_drag!(cf, cm, system, r, c, cl_to_cd)

Return the updated force coefficients per unit span for each spanwise segment using 
strip theory to account for viscous affects

This function requires that a near-field analysis has been performed on `system`
to obtain panel forces.

# Arguments
 - `cf`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns) which contains the x, y, and z components of 
    force coefficients (per unit span) for each spanwise segment.
 - `cm`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns) which contains the x, y, and z components of
    moment coefficients (per unit span) for each spanwise segment.
 - `r`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns+1) which contains the x, y, and z coordinates
    of the resulting lifting line coordinates
 - `c`: Vector with length equal to the number of surfaces, with each element
    being a vector of length `ns+1` which contains the chord lengths at each
    lifting line coordinate.
 - `cl_2_cd`: A function accepting lift coefficient and reynolds number, and returning
    viscous drag coefficent.

# Return Arguments:
 - `cf`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns) which contains the x, y, and z direction
    force coefficients (per unit span) in the wind frame for each spanwise segment.  
    The x element of these matrices includes a viscous drag correction.
"""
function strip_theory_drag!(cf, Res)  # TODO: attach airfoils to surface
    # iterate over each surface
    nsurf = length(cf)
    CDv = 0.0
    CM = 0.0
    
    for i = 1:nsurf
        # extract for this surface
        forces = cf[i]
        af1 = surface.  # TODO
        af2 = surface.
        Re = Res[i]
        y = surface. # TODO
        
        # broadcast across each airfoil, cl, Re
        viscous = cl_2_cd.(Ref(af1), Ref(af2), forces[3, :], Re, y)
        cd = getindex.(viscous, 1)
        cm = getindex.(viscous, 1)

        forces[1, :] .+= cd
        forces[2, :] .+= cm

        multiplier = symmetric[i] ? 2.0 : 1.0
        CDv += multiplier * cd * c * panel_width / ref.S
        # TODO: add cm
    end

    return cf, CDv
end

struct AirfoilData{TM, TV, TI}
    cl::TM
    cd::TM
    cm::TM
    Re::TV
    clcache::TV
    cdcache::TV
    cmcache::TV
    nalpha::TI
end

function AirfoilData(cl, cd, cm, Re)
    nalpha = length(cl[:, 1])
    clcache = zeros(nalpha)
    cdcache = zeros(nalpha)
    cmcache = zeros(nalpha)

    return AirfoilData(cl, cd, cm, Re, clcache, cdcache, cmcache, nalpha)
end


function cl_2_cd(af1, af2, cl, Re, eta)

    cd1, cm1 = cl_2_cd(af1, cl, Re)
    cd2, cm2 = cl_2_cd(af2, cl, Re)
    cd = cd1 + eta*(cd2 - cd1)
    cm = cm1 + eta*(cm2 - cm1)

    return cd, cm
end

function cl_2_cd(af, cl, Re)

    for i in 1:af.nalpha
        af.clcache[i] = linear(af.Res, view(af.cl, i, :), Re)
        af.cdcache[i] = linear(af.Res, view(af.cd, i, :), Re)
        af.cmcache[i] = linear(af.Res, view(af.cm, i, :), Re)
    end
    cd = linear(af.clcache, af.cdcache, cl)
    cm = linear(af.clcache, af.cmcache, cl)
    
    return cd, cm
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
