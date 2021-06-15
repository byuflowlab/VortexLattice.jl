"""
    surface_velocities!(Vcp, Vh, Vv, surface_velocities)

Calculate the velocities experienced at the control points `Vcp`, the horizontal
bound vortex centers `Vh`, and the vertical bound vortex centers `Vv` due to
surface motion/deformation.
"""
surface_velocities!

# single surface, surface panels input
function surface_velocities!(Vcp, Vh, Vv, surface_velocities::AbstractMatrix{<:SurfacePanel})
    nc, ns = size(surface_velocities)
    # velocity at the control points
    for j = 1:ns, i = 1:nc
        Vcp[i,j] = -controlpoint(surface_velocities[i,j])
    end
    # velocity at the horizontal bound vortices
    for j = 1:ns
        for i = 1:nc
            Vh[i,j] = -top_center(surface_velocities[i,j])
        end
        Vh[end,j] = -bottom_center(surface_velocities[end,j])
    end
    # velocity at the vertical bound vortices
    for i = 1:nc
        for j = 1:ns
            Vv[i,j] = -left_center(surface_velocities[i,j])
        end
        Vv[i,end] = -right_center(surface_velocities[i,end])
    end
    return Vcp, Vh, Vv
end

# single surface, grid input
function surface_velocities!(Vcp, Vh, Vv, surface_velocities::AbstractArray{TF, 3}) where TF
    nc = size(surface_velocities, 2) - 1
    ns = size(surface_velocities, 3) - 1
    for j = 1:ns
        # leading edge panel corner velocities
        r1n = SVector(surface_velocities[1,1,j], surface_velocities[2,1,j], surface_velocities[3,1,j]) # top left velocity
        r2n = SVector(surface_velocities[1,1,j+1], surface_velocities[2,1,j+1], surface_velocities[3,1,j+1]) # top right velocity
        r3n = SVector(surface_velocities[1,2,j], surface_velocities[2,2,j], surface_velocities[3,2,j]) # bottom left velocity
        r4n = SVector(surface_velocities[1,2,j+1], surface_velocities[2,2,j+1], surface_velocities[3,2,j+1]) # bottom right velocity
        for i = 1:nc-1
            # corner velocities of current panel
            r1 = r1n # top left of panel
            r2 = r2n # top right of panel
            r3 = r3n # bottom left of panel
            r4 = r4n # bottom right of panel
            # corner velocities of next panel
            r1n = r3 # top left
            r2n = r4 # top right
            r3n = SVector(surface_velocities[1,i+2,j], surface_velocities[2,i+2,j], surface_velocities[3,i+2,j]) # bottom left
            r4n = SVector(surface_velocities[1,i+2,j+1], surface_velocities[2,i+2,j+1], surface_velocities[3,i+2,j+1]) # bottom right
            # corner velocities of ring vortex
            rtl = linearinterp(0.25, r1, r3)
            rtr = linearinterp(0.25, r2, r4)
            rbl = linearinterp(0.25, r1n, r3n)
            rbr = linearinterp(0.25, r2n, r4n)
            # velocity experienced at the control point
            rtop = linearinterp(0.5, r1, r2)
            rbot = linearinterp(0.5, r3, r4)
            Vcp[i,j] = -linearinterp(0.75, rtop, rbot)
            # velocity experienced at the bound vortices
            Vh[i,j] = -linearinterp(0.5, rtl, rtr)
            Vv[i,j] = -linearinterp(0.5, rbl, rtl)
            Vv[i,j+1] = -linearinterp(0.5, rbr, rtr)
        end
        # corner velocities of current panel
        r1 = r1n # top left
        r2 = r2n # top right
        r3 = r3n # bottom left
        r4 = r4n # bottom right
        # corner velocities of ring vortex
        rtl = linearinterp(0.25, r1, r3)
        rtr = linearinterp(0.25, r2, r4)
        rbl = r3
        rbr = r4
        # velocity experienced at the control point
        rtop = linearinterp(0.5, r1, r2)
        rbot = linearinterp(0.5, r3, r4)
        Vcp[i,j] = -linearinterp(0.75, rtop, rbot)
        # velocity experienced at the bound vortices
        Vh[nc,j] = -linearinterp(0.5, rtl, rtr)
        Vv[nc,j] = -linearinterp(0.5, rbl, rtl)
        Vv[nc,j+1] = -linearinterp(0.5, rbr, rtr)
        # velocity experienced at the trailing edge
        Vh[nc+1,j] = -rbot
    end
    return Vcp, Vh, Vv
end

# multiple surfaces
function surface_velocities!(Vcp, Vh, Vv, surface_velocities)
    for isurf = 1:length(surfaces)
        surface_velocities!(Vcp[isurf], Vh[isurf], Vv[isurf], surface_velocities[isurf])
    end
    return Vcp, Vh, Vv
end

"""
    steady_normal_velocities!(w, surfaces, reference, freestream;
        symmetric, surface_id, trailing_vortices, xhat)

Compute the downwash on the surface panels in `surfaces` due to all velocity
sources except for induced velocities from the surface and wake panels.
"""
function steady_normal_velocities!(w, surfaces, ref, fs, additional_velocity;
    symmetric, surface_id, trailing_vortices, xhat)

    # number of surfaces
    nsurf = length(surfaces)

    # index for keeping track of where we are in the w vector
    iw = 0

    # loop through receiving surfaces
    for isurf = 1:nsurf

        # current surface
        surface = surfaces[isurf]

        # iterate through all control points on this surface
        for i = 1:length(surface)

            # control point
            rcp = controlpoint(surface[i])

            # normal vector
            nhat = normal(surface[i])

            # freestream velocity
            V = freestream_velocity(fs)

            # rotational velocity
            V += rotational_velocity(rcp, fs, ref)

            # additional velocity field
            V += SVector{3}(additional_velocity(rcp))

            # get normal component
            w[iw+i] = -dot(V, nhat)
        end

        # increment position in AIC matrix
        iw += length(surface)
    end

    return w
end

"""
    unsteady_normal_velocities!(w, surfaces, wakes, reference, freestream, Vcp;
        symmetric, surface_id, wake_finite_core, iwake, trailing_vortices, xhat)

Compute the downwash on the surface panels in `surfaces` due to all velocity
sources except for induced velocities from the surface panels.
"""
function unsteady_normal_velocities!(w, surfaces, wakes, ref, fs,
    additional_velocity, Vcp; symmetric, surface_id, wake_finite_core, iwake,
    trailing_vortices, xhat)

    # number of surfaces
    nsurf = length(surfaces)

    # index for keeping track of where we are in the w vector
    iw = 0

    # calculate steady portion of normal velocities
    steady_normal_velocities!(w, surfaces, ref, fs, additional_velocity;
        symmetric, surface_id, trailing_vortices, xhat)

    # add unsteady portion of normal velocities
    for isurf = 1:nsurf

        # current surface
        surface = surfaces[isurf]

        # iterate through all control points on this surface
        for i = 1:length(surface)

            # control point
            rcp = controlpoint(surface[i])

            # normal vector
            nhat = normal(surface[i])

            # velocity due to surface motion
            V = Vcp[isurf][i]

            # velocity due to the wake panels
            for jsurf = 1:nsurf
                if iwake[jsurf] > 0
                    V += induced_velocity(rcp, wakes[jsurf];
                        nc = iwake[jsurf],
                        finite_core = wake_finite_core[isurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end
            end

            # get normal component
            w[iw+i] -= dot(V, nhat)
        end

        # increment position in AIC matrix
        iw += length(surface)
    end

    return w
end

"""
    steady_normal_velocities_and_derivatives!(w, dw, surfaces, reference, freestream;
        symmetric, surface_id, trailing_vortices, xhat)

Compute the downwash on the surface panels in `surfaces` due to all velocity
sources except for induced velocities from the surface and wake panels.  Also
compute the derivative of the downwash with respect to the freestream parameters.
"""
function steady_normal_velocities_and_derivatives!(w, dw, surfaces, ref, fs,
    additional_velocity; symmetric, surface_id, trailing_vortices, xhat)

    # number of surfaces
    nsurf = length(surfaces)

    # unpack derivatives
    (w_a, w_b, w_p, w_q, w_r) = dw

    # index for keeping track of where we are in the w vector
    iw = 0

    # loop through receiving surfaces
    for isurf = 1:nsurf

        # current surface
        surface = surfaces[isurf]

        # iterate through all control points on this surface
        for i = 1:length(surface)

            # control point
            rcp = controlpoint(surface[i])

            # normal vector
            nhat = normal(surface[i])

            # freestream velocity and its derivatives
            V, dVf = freestream_velocity_derivatives(fs)

            V_a, V_b = dVf

            # rotational velocity and its derivatives
            Vrot, dVrot = rotational_velocity_derivatives(rcp, fs, ref)

            V += Vrot

            V_p, V_q, V_r = dVrot

            # additional velocity field
            V += SVector{3}(additional_velocity(rcp))

            # get normal component
            w[iw+i] = -dot(V, nhat)

            # associated derivatives
            w_a[iw+i] = -dot(V_a, nhat)
            w_b[iw+i] = -dot(V_b, nhat)
            w_p[iw+i] = -dot(V_p, nhat)
            w_q[iw+i] = -dot(V_q, nhat)
            w_r[iw+i] = -dot(V_r, nhat)

        end

        # increment position in AIC matrix
        iw += length(surface)
    end

    # pack up derivatives
    dw = (w_a, w_b, w_p, w_q, w_r)

    return w, dw
end

"""
    unsteady_normal_velocities_and_derivatives!(w, dw, surfaces, wakes,
        reference, freestream, additional_velocity, Vcp; symmetric, surface_id,
        wake_finite_core, iwake, trailing_vortices, xhat)

Compute the downwash on the surface panels in `surfaces` due to all velocity
sources except for induced velocities from the surface panels.  Also compute the
derivative of the downwash with respect to the freestream parameters.
"""
function unsteady_normal_velocities_and_derivatives!(w, dw, surfaces, wakes,
    ref, fs, additional_velocity, Vcp; symmetric, surface_id, wake_finite_core, iwake,
    trailing_vortices, xhat)

    # number of surfaces
    nsurf = length(surfaces)

    # index for keeping track of where we are in the w vector
    iw = 0

    # calculate steady portion of normal velocities (and derivatives)
    steady_normal_velocities_and_derivatives!(w, dw, surfaces, ref, fs,
        additional_velocity; symmetric, surface_id, trailing_vortices, xhat)

    # add unsteady portion of normal velocities
    for isurf = 1:nsurf

        # current surface
        surface = surfaces[isurf]

        # iterate through all control points on this surface
        for i = 1:length(surface)

            # control point
            rcp = controlpoint(surface[i])

            # normal vector
            nhat = normal(surface[i])

            # velocity due to surface motion
            V = Vcp[isurf][i]

            # velocity due to the wake panels
            for jsurf = 1:nsurf
                if iwake[jsurf] > 0
                    V += induced_velocity(rcp, wakes[jsurf];
                        nc = iwake[jsurf],
                        finite_core = wake_finite_core[isurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end
            end

            # get normal component
            w[iw+i] -= dot(V, nhat)
        end

        # increment position in AIC matrix
        iw += length(surface)
    end

    return w, dw
end

"""
    circulation!(Γ, AIC, w)

Calculate the circulation strength of the surface panels.
"""
circulation!(Γ, AIC, w) = ldiv!(Γ, lu(AIC), w)

"""
    circulation_and_derivatives!(Γ, dΓ, AIC, w, dw)

Calculate the circulation strength of the surface panels and its derivatives
with respect to the freestream variables
"""
function circulation_and_derivatives!(Γ, dΓ, AIC, w, dw)
    # unpack derivatives
    (Γ_a, Γ_b, Γ_p, Γ_q, Γ_r) = dΓ
    (w_a, w_b, w_p, w_q, w_r) = dw
    # factorize AIC matrix
    fAIC = lu(AIC)
    # solve for circulation and its derivatives
    ldiv!(Γ, fAIC, w)
    ldiv!(Γ_a, fAIC, w_a)
    ldiv!(Γ_b, fAIC, w_b)
    ldiv!(Γ_p, fAIC, w_p)
    ldiv!(Γ_q, fAIC, w_q)
    ldiv!(Γ_r, fAIC, w_r)
    return Γ, dΓ
end
