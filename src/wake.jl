"""
    update_wake_shedding_locations!(wakes, wake_shedding_locations,
        surfaces, ref, fs, dt, additional_velocity, Vte, nwake, eta)

Update the wake shedding locations.  Also update the first chordwise wake panels
to account for the new wake shedding location

# Arguments
 - `wakes`: Vector of wakes corresponding to each surface, represented by matrices
    of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where `nw` is the
    number of chordwise wake panels and `ns` is the number of spanwise panels.
 - `wake_shedding_locations`: Shedding location coordinates for each surface for
    each trailing edge vertex.
 - `surfaces`: Vector of surfaces, represented by matrices of surface panels
    (see [`SurfacePanel`](@ref) of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))
 - `dt`: Time step (seconds)
 - `additional_velocity`: Function defining additional velocity field
 - `Vte`: Velocity experienced at the trailing edge due to surface motion.
 - `nwake`: Number of chordwise wake panels to use from each wake in `wakes`
 - `eta`: Time step fraction used to define separation between trailing
    edge and wake shedding location.  Typical values range from 0.2-0.3.
"""
@inline function update_wake_shedding_locations!(wakes, wake_shedding_locations,
    surfaces, ref, fs, dt, additional_velocity, Vte, nwake, eta)

    # get number of surfaces
    nsurf = length(surfaces)

    # loop through all surfaces
    for isurf = 1:nsurf

        # number of spanwise panels
        ns = size(wakes[isurf], 2)

        # update wake shedding location
        for j = 1:ns+1

            # extract trailing edge coordinate
            if j < ns + 1
                rte = bottom_left(surfaces[isurf][end, j])
            else
                rte = bottom_right(surfaces[isurf][end, j-1])
            end

            # freestream velocity
            V = freestream_velocity(fs)

            # rotational velocity
            V += rotational_velocity(rte, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                V += additional_velocity(rte)
            end

            # velocity due to surface motion
            if !isnothing(Vte)
                V += Vte[isurf][j]
            end

            # update wake shedding location coordinates
            wake_shedding_locations[isurf][j] = rte + eta*V*dt

        end

        if nwake[isurf] > 0
            # loop through first row of wake panels
            for j = 1:ns
                # update wake panel with wake shedding location coordinates
                rtl = wake_shedding_locations[isurf][j]
                rtr = wake_shedding_locations[isurf][j+1]

                # preserve other wake panel coordinates
                rbl = bottom_left(wakes[isurf][1,j])
                rbr = bottom_right(wakes[isurf][1,j])

                # preserve core size
                core_size = get_core_size(wakes[isurf][1,j])

                # preserve circulation strength
                gamma = circulation_strength(wakes[isurf][1,j])

                # replace the old wake panel
                wakes[isurf][1,j] = WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
            end
        end

    end

    return wakes, wake_shedding_locations
end

"""
    get_wake_velocities!(wake_velocities, surfaces, wakes, ref, fs, Γ,
        additional_velocity, Vte, symmetric, repeated_points, nwake,
        surface_id, wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)

# Arguments
 - `wake_velocities`: Velocities at the corners of the wake panels in `wakes`
 - `surfaces`: Vector of surfaces, represented by matrices of surface panels
    (see [`SurfacePanel`](@ref) of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `wakes`: Vector of wakes corresponding to each surface, represented by matrices
    of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where `nw` is the
    number of chordwise wake panels and `ns` is the number of spanwise panels.
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))
 - `Γ`: Circulation of all surface panels stored in a single vector
 - `additional_velocity`: Function defining additional velocity field
 - `Vte`: Velocity at the trailing edge vertices on each surface due to surface motion
 - `symmetric`: (required) Flag for each surface indicating whether a mirror
    image across the X-Z plane should be used when calculating induced velocities
 - `repeated_points`: Dictionary of the form `Dict((isurf, i) => [(jsurf1, j1),
    (jsurf2, j2)...]` which defines repeated trailing edge points.  Trailing edge
    point `i` on surface `isurf` is repeated on surface `jsurf1` at point `j1`,
    `jsurf2` at point `j2`, and so forth. See [`repeated_trailing_edge_points`](@ref)
 - `nwake`: Number of chordwise wake panels to use from each wake in `wakes`,
    defaults to all provided wake panels
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.  Defaults to `true` for each surface.
 - `wake_shedding_locations`: Shedding location coordinates for each surface for
    each trailing edge vertex.
 - `trailing_vortices`: Flags to enable/disable trailing vortices, defaults to
    `true` for each surface
 - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]
"""
@inline function get_wake_velocities!(wake_velocities, surfaces, wakes, ref, fs, Γ,
    additional_velocity, Vte, symmetric, repeated_points, nwake,
    surface_id, wake_finite_core, wake_shedding_locations, trailing_vortices, xhat)

    # number of surfaces
    nsurf = length(surfaces)

    # loop through all surfaces
    for isurf = 1:nsurf

        # number of chordwise and spanwise panels
        nc, ns = size(surfaces[isurf])
        nw = nwake[isurf]

        # velocity at the wake shedding locations
        for is = 1:ns+1

            # check if this point is a duplicate, skip if it is
            if (isurf, is) in keys(repeated_points)
                for (jsurf, js) in repeated_points[(isurf, is)]
                    # NOTE: we assume that a point is not repeated on the same surface
                    if jsurf < isurf
                        wake_velocities[isurf][1, is] = wake_velocities[jsurf][1, js]
                        continue
                    end
                end
            end

            # get vertex location
            rc = wake_shedding_locations[isurf][is]

            # freestream velocity
            V = freestream_velocity(fs)

            # rotational velocity
            V += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                V += additional_velocity(rc)
            end

            # velocity at the trailing edge
            wake_velocities[isurf][1,is] = V
        end

        # velocity at all other wake vertices
        cr = CartesianIndices((2:nw+1, 1:ns+1))

        # loop through all vertices
        for I in cr

            # check if this point is a duplicate, skip if it is
            if (isurf, I[2]) in keys(repeated_points)
                for (jsurf, js) in repeated_points[(isurf, I[2])]
                    # NOTE: we assume that a point is not repeated on the same surface
                    if jsurf < isurf
                        wake_velocities[isurf][I] = wake_velocities[jsurf][I[1], js]
                        continue
                    end
                end
            end

            # get vertex location
            if I[1] <= nw && I[2] <= ns
                rc = top_left(wakes[isurf][I[1], I[2]])
            elseif I[1] == nw + 1 && I[2] <= ns
                rc = bottom_left(wakes[isurf][I[1]-1, I[2]])
            elseif I[1] <= nw && I[2] == ns + 1
                rc = top_right(wakes[isurf][I[1], I[2]-1])
            else # I[1] == nw + 1 && I[2] == ns + 1
                rc = bottom_right(wakes[isurf][I[1]-1, I[2]-1])
            end

            # freestream velocity
            wake_velocities[isurf][I] = freestream_velocity(fs)

            # rotational velocity
            wake_velocities[isurf][I] += rotational_velocity(rc, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                wake_velocities[isurf][I] += additional_velocity(rc)
            end

            # induced velocity from each surface and wake
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                # number of panels on sending surface
                Ns = length(surfaces[jsurf])

                # check if receiving point is repeated on the sending surface
                if isurf == jsurf
                    # the surfaces are the same
                    same_surface = true
                    # vertex spanwise coordinate
                    js = I[2]
                else
                    # the surfaces are different, but the point could still be a duplicate
                    if (isurf, I[2]) in keys(repeated_points)
                        # the vertex is duplicated

                        # check if the vertex is on the sending surface
                        idx = findfirst(x -> x[1] == jsurf, repeated_points[(isurf, I[2])])

                        if isnothing(idx)
                            # the vertex is not duplicated on the sending surface
                            same_surface = false
                        else
                            # the vertex is duplicated on the sending surface
                            same_surface = true
                            # vertex spanwise coordinate on sending surface
                            js = repeated_points[(isurf, I[2])][idx][2]
                        end
                    else
                        # the vertex is not duplicated

                        # the vertex is not on the sending surface
                        same_surface = false
                    end
                end

                # extract circulation values corresponding to the sending surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)

                # induced velocity from this surface
                wake_velocities[isurf][I] += induced_velocity(rc, surfaces[jsurf], vΓ;
                    finite_core = surface_id[isurf] != surface_id[jsurf],
                    wake_shedding_locations = wake_shedding_locations[jsurf],
                    symmetric = symmetric[jsurf],
                    trailing_vortices = false,
                    xhat = xhat)

                # add induced velocity from the wake
                if same_surface
                    # vertex location on wake
                    J = CartesianIndex(I[1], js)

                    # induced velocity from wake on its own vertex
                    wake_velocities[isurf][I] += induced_velocity(J, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                else
                    # induced velocity from wake on another wake's vertex
                    wake_velocities[isurf][I] += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end

                jΓ += Ns # increment Γ index
            end
        end
    end

    return wake_velocities
end

"""
    translate_wake(panel, wake_velocities, dt)

Return a translated copy of the wake panel `panel` given the wake corner velocities
`wake_velocities` and the time step `dt`

# Arguments
 - `panel`: Wake panel (of type [`WakePanel`](@ref))
 - `wake_velocities`: Matrix containing the velocities at each of the four corners
    of `panel`
 - `dt`: Time step (seconds)
"""
@inline function translate_wake(panel::WakePanel, wake_velocities, dt)

    # extract corners
    rtl = top_left(panel)
    rtr = top_right(panel)
    rbl = bottom_left(panel)
    rbr = bottom_right(panel)

    # get vortex filament length
    lt = norm(rtr - rtl)
    lb = norm(rbl - rbr)
    ll = norm(rtl - rbl)
    lr = norm(rbr - rtr)
    l1 = lt + lb + ll + lr

    # translate corners
    rtl += wake_velocities[1,1]*dt
    rtr += wake_velocities[1,2]*dt
    rbl += wake_velocities[2,1]*dt
    rbr += wake_velocities[2,2]*dt

    # get new vortex filament length
    lt = norm(rtr - rtl)
    lb = norm(rbl - rbr)
    ll = norm(rtl - rbl)
    lr = norm(rbr - rtr)
    l2 = lt + lb + ll + lr

    # use previous core size
    core_size = get_core_size(panel)

    # correct vorticity for vortex stretching
    gamma = circulation_strength(panel)*l1/l2

    return WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
end

"""
    translate_wake!(wake, wake_velocities, dt; nwake = size(wake, 1))

Translate the wake panels in `wake` given the corner velocities `wake_velocities`
and the time step `dt`.

# Arguments
 - `wake`: Matrix of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns)
    where `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels, defaults to no wake panels
 - `wake_velocities`: Velocities at each of the vertices corresponding to the
    wake panels in `wake`
 - `dt`: Time step

# Keyword Arguments
 - `nwake`: Number of chordwise wake panels to use from `wake`, defaults to all
    provided wake panels
"""
@inline function translate_wake!(wake, wake_velocities, dt; nwake = size(wake, 1))

    nw = nwake
    ns = size(wake, 2)
    cw = CartesianIndices((nw, ns))

    for I in cw

        panel = wake[I]

        vV = view(wake_velocities, I[1]:I[1]+1, I[2]:I[2]+1)

        wake[I] = translate_wake(panel, vV, dt)
    end

    return wake
end

"""
    shed_wake!(wake, wake_shedding_locations, wake_velocities, dt, surface, Γ, nwake)

Shed a new wake panel from the wake shedding locations and translate existing
wake panels.

# Arguments
 - `wake`: Matrix of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns)
    where `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels
 - `wake_shedding_locations`: Vector of length `ns` which stores the coordinates
    where wake panels are shed from the trailing edge of `surface`.
 - `wake_velocities`: Velocities at each of the vertices corresponding to the
    wake panels in `wake`
 - `dt`: Time step (seconds)
 - `surface`: Matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `Γ`: Circulation strength of each surface panel in `surface`
 - `nwake`: Number of chordwise wake panels to use from `wake`, defaults to all
    provided wake panels
"""
@inline function shed_wake!(wake::AbstractMatrix, wake_shedding_locations,
    wake_velocities, dt, surface, Γ, nwake)

    nc, ns = size(surface)
    nw = size(wake, 1)
    ls = LinearIndices((nc, ns))

    # replace the last chordwise panels with the newly shed wake panels
    for j = 1:ns

        # shedding location coordinates
        rtl = wake_shedding_locations[j]
        rtr = wake_shedding_locations[j+1]

        # shed coordinates
        rbl = rtl + wake_velocities[1, j]*dt
        rbr = rtr + wake_velocities[1, j+1]*dt

        # use core size from the shedding panel
        core_size = get_core_size(surface[end, j])

        # use circulation strength from the shedding panel
        gamma = Γ[ls[end,j]]

        # replace the oldest wake panel
        wake[end,j] = WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
    end

    # translate all existing wake panels except the most recently shed panels
    translate_wake!(wake, wake_velocities, dt, nwake = min(nwake, nw-1))

    # shift wake panels to make newly shed wake panel first
    rowshift!(wake)

    return wake
end

"""
    shed_wake!(wakes, wake_shedding_locations, wake_velocities, dt, surfaces, Γ, nwake)

Shed a new wake panel from the wake shedding locations and translate existing
wake panels.

# Arguments
 - `wakes`: Vector of wakes corresponding to each surface, represented by matrices
    of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where `nw` is the
    number of chordwise wake panels and `ns` is the number of spanwise panels.
 - `wake_shedding_locations`: Shedding location coordinates for each surface for
    each trailing edge vertex.
 - `wake_velocities`: Velocities at each of the vertices corresponding to the
    wake panels in `wake`
 - `dt`: Time step (seconds)
 - `surfaces`: Vector of surfaces, represented by matrices of surface panels
    (see [`SurfacePanel`](@ref) of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `Γ`: Circulation strength of each surface panel in `surfaces`
 - `nwake`: Number of chordwise wake panels to use from each wake in `wakes`,
    defaults to all provided wake panels
"""
@inline function shed_wake!(wakes::AbstractVector{<:AbstractMatrix}, wake_shedding_locations,
    wake_velocities, dt, surfaces, Γ, nwake)

    iΓ = 0
    for i = 1:length(surfaces)

        N = length(surfaces[i])

        vΓ = view(Γ, iΓ+1:iΓ+N)

        shed_wake!(wakes[i], wake_shedding_locations[i], wake_velocities[i], dt,
            surfaces[i], vΓ, nwake[i])

        iΓ += N
    end

    return wakes
end

"""
    rowshift!(A)

Circularly shifts the rows of a matrix down one row.
"""
function rowshift!(A)

    ni, nj = size(A)

    for j = 1:nj
        tmp = A[ni,j]
        for i = ni:-1:2
            A[i,j] = A[i-1,j]
        end
        A[1,j] = tmp
    end

    return A
end
