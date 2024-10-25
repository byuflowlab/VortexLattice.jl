"""
    get_surface_velocities!(Vcp, Vh, Vv, Vte, current_surface, previous_surface, dt)

Calculate the velocities experienced by the surface at the control points, horizontal
bound vortex centers, vertical bound vortex centers and trailing edge vertices due
to surface motion.
"""
function get_surface_velocities!(Vcp, Vh, Vv, Vte, current_surface, previous_surface, dt)

    nc, ns = size(current_surface)

    # velocity at the control points
    for j = 1:ns, i = 1:nc
        Vcp[i,j] = (controlpoint(previous_surface[i,j]) - controlpoint(current_surface[i,j]))/dt
    end

    # velocity at the horizontal bound vortices
    for j = 1:ns
        for i = 1:nc
            Vh[i,j] = (top_center(previous_surface[i,j]) - top_center(current_surface[i,j]))/dt
        end
        Vh[end,j] = (bottom_center(previous_surface[end,j]) - bottom_center(current_surface[end,j]))/dt
    end

    # velocity at the vertical bound vortices
    for i = 1:nc
        for j = 1:ns
            Vv[i,j] = (left_center(previous_surface[i,j]) - left_center(current_surface[i,j]))/dt
        end
        Vv[i,end] = (right_center(previous_surface[i,end]) - right_center(current_surface[i,end]))/dt
    end

    # velocity at the trailing edge vertices
    for j = 1:ns
        Vte[j] = (bottom_left(previous_surface[end,j]) - bottom_left(current_surface[end,j]))/dt
    end
    Vte[end] = (bottom_right(previous_surface[end,end]) - bottom_right(current_surface[end,end]))/dt

end

"""
    normal_velocity!(w, surfaces, wakes, ref, fs; additional_velocity,
        Vcp, symmetric, nwake, surface_id, wake_finite_core, trailing_vortices, xhat)

Compute the downwash at the control points on `surfaces` due to the freestream
velocity, rotational velocity, additional velocity field, surface motion, and
induced velocity from the wake panels.

This forms the right hand side of the circulation linear system solve.
"""
function normal_velocity!(w, surfaces, wakes, ref, fs; additional_velocity,
    Vcp, symmetric, nwake, surface_id, wake_finite_core, trailing_vortices, xhat)

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
            if !isnothing(additional_velocity)
                V += additional_velocity(rcp)
            end

            # velocity due to surface motion
            if !isnothing(Vcp)
                V += Vcp[isurf][i]
            end

            # velocity due to the wake panels
            for jsurf = 1:length(wakes)
                if nwake[jsurf] > 0
                    V += induced_velocity(rcp, wakes[jsurf];
                        finite_core = wake_finite_core[isurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end
            end

            # get normal component
            w[iw+i] = -dot(V, nhat)
        end

        # increment position in AIC matrix
        iw += length(surface)
    end

    return w
end

"""
    normal_velocity_derivatives!(w, dw, surfaces, wakes, ref, fs;
        additional_velocity, Vcp, symmetric, nwake, surface_id,
        wake_finite_core, trailing_vortices, xhat)

Compute the downwash at the control points on `surfaces` due to the freestream
velocity, rotational velocity, additional velocity field, surface motion, and
induced velocity from the wake panels.  Also calculate its derivatives with
respect to the freestream parameters.

This forms the right hand side of the circulation linear system solve (and its derivatives).
"""
function normal_velocity_derivatives!(w, dw, surfaces, wakes, ref, fs;
    additional_velocity, Vcp, symmetric, nwake, surface_id, wake_finite_core,
    trailing_vortices, xhat)

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
            if !isnothing(additional_velocity)
                V += additional_velocity(rcp)
            end

            # velocity due to surface motion
            if !isnothing(Vcp)
                V += Vcp[isurf][i]
            end

            # velocity due to the wake panels
            for jsurf = 1:length(wakes)
                if nwake[jsurf] > 0
                    V += induced_velocity(rcp, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || (surface_id[isurf] != surface_id[jsurf]),
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end
            end

            # right hand side vector
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

# --- circulation solve --- #

"""
    circulation(AIC, w)

Solve for the circulation distribution.
"""
circulation(AIC, w) = AIC\w

"""
    circulation!(Γ, AIC, w)

Pre-allocated version of `circulation`
"""
circulation!(Γ, AIC, w) = ldiv!(Γ, lu(AIC), w)

"""
    circulation_derivatives(AIC, w, dw)

Solve for the circulation distribution and its derivatives with respect to
the freestream parameters.
"""
function circulation_derivatives(AIC, w, dw)

    # unpack derivatives
    (w_a, w_b, w_p, w_q, w_r) = dw

    # factorize AIC matrix (since we'll be reusing it)
    fAIC = factorize(AIC)

    # solve for circulation and its derivatives
    Γ = fAIC\w

    Γ_a = fAIC\w_a
    Γ_b = fAIC\w_b
    Γ_pb = fAIC\w_p
    Γ_qb = fAIC\w_q
    Γ_rb = fAIC\w_r

    # pack up derivatives
    dΓ = (Γ_a, Γ_b, Γ_pb, Γ_qb, Γ_rb)

    return Γ, dΓ
end

"""
    circulation_derivatives!(Γ, dΓ, AIC, w, dw)

Pre-allocated version of `circulation_derivatives`
"""
function circulation_derivatives!(Γ, dΓ, AIC, w, dw)

    # unpack derivatives
    (Γ_a, Γ_b, Γ_p, Γ_q, Γ_r) = dΓ
    (w_a, w_b, w_p, w_q, w_r) = dw

    # factorize AIC matrix (since we'll be reusing it)
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
