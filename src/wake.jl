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
function update_wake_shedding_locations!(wakes, wake_shedding_locations,
    surfaces, ref, fs, dt, additional_velocity, Vte, nwake, eta)

    # get number of surfaces
    nsurf = length(surfaces)

    # loop through all surfaces
    for isurf = 1:nsurf

        # number of spanwise panels
        ns = length(wake_shedding_locations[isurf]) - 1

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
                # preserve other wake panel coordinates
                rtl = wake_shedding_locations[isurf][j]
                rtr = wake_shedding_locations[isurf][j+1]
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
    Similar to `update_wake_shedding_locations`, but assumes `wake_shedding_locations` contain the previous trailing edge location.
"""
function update_wake_shedding_locations_unsteady!(wakes, wake_shedding_locations,
    surfaces, ref, fs, dt, additional_velocity, Vte, nwake, eta; sync_panels=true)

    # get number of surfaces
    nsurf = length(surfaces)

    # loop through all surfaces
    for isurf = 1:nsurf

        # number of spanwise panels
        ns = length(wake_shedding_locations[isurf]) - 1
        wsl = wake_shedding_locations[isurf]
        surface = surfaces[isurf]
        wake = wakes[isurf]

        # update wake shedding location
        for j = 1:ns+1

            # extract trailing edge coordinate
            if j < ns + 1
                rte = bottom_left(surface[end, j])
            else
                rte = bottom_right(surface[end, j-1])
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
            wsl[j] = rte + eta*V*dt

        end

        if sync_panels && nwake[isurf] > 0
            # loop through first row of wake panels
            for j = 1:ns
                # update wake panel with wake shedding location coordinates
                rtl = wsl[j]
                rtr = wsl[j+1]

                # preserve other wake panel coordinates
                rbl = bottom_left(wakes[isurf][1,j])
                rbr = bottom_right(wakes[isurf][1,j])

                # preserve core size
                core_size = get_core_size(wakes[isurf][1,j])

                # preserve circulation strength
                gamma = circulation_strength(wakes[isurf][1,j])

                # replace the old wake panel
                wake[1,j] = WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
            end
        end

    end

    return wakes, wake_shedding_locations
end

function initial_wake_panels!(wakes, wake_shedding_locations, surfaces, Γ, eta)
    for isurf in eachindex(surfaces)
        wsl = wake_shedding_locations[isurf]
        surface = surfaces[isurf]
        wake = wakes[isurf]
        nc, ns = size(surface)
        ls = LinearIndices((nc, ns))

        # get trailing edge point
        rte = bottom_left(surface[nc,1])

        # wsl point
        rwsl = wsl[1]

        # extend from the trailing edge to the end of the wake panel
        dx = rwsl - rte
        rbl_wake = rte + dx / eta

        for j in 1:ns
            
            # get right trailing edge point
            rte = bottom_right(surface[nc,j])

            # wsl point
            rwsl = wsl[j+1]

            # extend from the trailing edge to the end of the wake panel
            dx = rwsl - rte
            rbr_wake = rte + dx / eta

            # update wake panel with wake shedding location coordinates
            rtl = wsl[j]
            rtr = wsl[j+1]
            core_size = get_core_size(surface[nc, j])
            gamma = Γ[ls[nc, j]]

            # replace the old wake panel
            wake[1,j] = WakePanel(rtl, rtr, rbl_wake, rbr_wake, core_size, gamma)

            # recurse rbl
            rbl_wake = rbr_wake

        end
    end
end

function store_trailing_edge!(wake_shedding_locations, surfaces)

    # get number of surfaces
    nsurf = length(surfaces)

    # loop through all surfaces
    for isurf = 1:nsurf

        # number of spanwise panels
        ns = length(wake_shedding_locations[isurf]) - 1

        # update wake shedding location
        for j = 1:ns+1

            # extract trailing edge coordinate
            if j < ns + 1
                rte = bottom_left(surfaces[isurf][end, j])
            else
                rte = bottom_right(surfaces[isurf][end, j-1])
            end

            # update wake shedding location coordinates
            wake_shedding_locations[isurf][j] = rte

        end

    end

    return wake_shedding_locations
end

function update_vpm_shedding_locations!(wakes, ref, fs, dt, additional_velocity, Vwake)

    # get number of surfaces
    nsurf = length(wakes)

    # loop through all surfaces
    for isurf = 1:nsurf

        # number of spanwise panels
        ns = size(wakes[isurf], 2)

        #--- get left velocity ---#

        # extract trailing edge coordinate
        rte = wakes[isurf][1,1].rtl

        # freestream velocity
        V = freestream_velocity(fs)

        # rotational velocity
        V += rotational_velocity(rte, fs, ref)

        # additional velocity field
        if !isnothing(additional_velocity)
            V += additional_velocity(rte)
        end

        # velocity due to surface motion
        if !isnothing(Vwake)
            V += Vwake[isurf][1,1]
        end

        # update vpm shedding location
        new_rbl = rte + V*dt

        # update vpm shedding location
        for j = 1:ns

            # extract trailing edge coordinate
            rte = wakes[isurf][1,j].rtr

            # freestream velocity
            V = freestream_velocity(fs)

            # rotational velocity
            V += rotational_velocity(rte, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                V += additional_velocity(rte)
            end

            # velocity due to surface motion
            if !isnothing(Vwake)
                V += Vwake[isurf][1,j+1]
            end

            # update wake shedding location coordinates
            new_rbr = rte + V*dt

            # preserve other wake panel coordinates
            rtl = top_left(wakes[isurf][1,j])
            rtr = top_right(wakes[isurf][1,j])

            # preserve core size
            core_size = get_core_size(wakes[isurf][1,j])

            # preserve circulation strength
            gamma = circulation_strength(wakes[isurf][1,j])

            # replace the old wake panel
            wakes[isurf][1,j] = WakePanel(rtl, rtr, new_rbl, new_rbr, core_size, gamma)

            # recurse new_rbl for next panel
            new_rbl = new_rbr

        end

    end

    return wakes
end

function update_vpm_shedding_TE!(wakes, ref, fs, dt, additional_velocity, Vwake)

    # get number of surfaces
    nsurf = length(wakes)

    # loop through all surfaces
    for isurf = 1:nsurf

        # number of spanwise panels
        ns = size(wakes[isurf], 2)

        #--- get left velocity ---#

        # extract trailing edge coordinate
        rte = wakes[isurf][1,1].rtl

        # freestream velocity
        V = freestream_velocity(fs)

        # rotational velocity
        V += rotational_velocity(rte, fs, ref)

        # additional velocity field
        if !isnothing(additional_velocity)
            V += additional_velocity(rte)
        end

        # velocity due to surface motion
        if !isnothing(Vwake)
            V += Vwake[isurf][1,1]
        end

        # update vpm shedding location
        new_rbl = rte + V*dt

        # update vpm shedding location
        for j = 1:ns

            # extract trailing edge coordinate
            rte = wakes[isurf][1,j].rtr

            # freestream velocity
            V = freestream_velocity(fs)

            # rotational velocity
            V += rotational_velocity(rte, fs, ref)

            # additional velocity field
            if !isnothing(additional_velocity)
                V += additional_velocity(rte)
            end

            # velocity due to surface motion
            if !isnothing(Vwake)
                V += Vwake[isurf][1,j+1]
            end

            # update wake shedding location coordinates
            new_rbr = rte + V*dt

            # preserve other wake panel coordinates
            rtl = top_left(wakes[isurf][1,j])
            rtr = top_right(wakes[isurf][1,j])

            # preserve core size
            core_size = get_core_size(wakes[isurf][1,j])

            # preserve circulation strength
            gamma = circulation_strength(wakes[isurf][1,j])

            # replace the old wake panel
            wakes[isurf][1,j] = WakePanel(rtl, rtr, new_rbl, new_rbr, core_size, gamma)

            # recurse new_rbl for next panel
            new_rbl = new_rbr

        end

    end

    return wakes
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
function get_wake_velocities!(wake_velocities, surfaces, wakes, ref, fs, Γ,
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
                v_surface = induced_velocity(rc, surfaces[jsurf], vΓ;
                    finite_core = surface_id[isurf] != surface_id[jsurf],
                    wake_shedding_locations = wake_shedding_locations[jsurf],
                    symmetric = symmetric[jsurf],
                    trailing_vortices = false,
                    xhat = xhat)
                wake_velocities[isurf][I] += v_surface

                # add induced velocity from the wake
                if same_surface
                    nc_wake = max(nwake[jsurf] - 1, 0)
                    # I[1] ranges up to nw+1 (the newly-shed/attached row), which is
                    # excluded by nc_wake = nw-1; induced_velocity's CartesianIndex
                    # branches assume I[1] <= nc_wake+1, so skip the last row here
                    # rather than falling through to an incorrect branch (this used
                    # to throw a BoundsError, e.g. wake[nw, 0], for I[1] == nw+1).
                    if nc_wake > 1 && I[1] <= nc_wake + 1
                        # vertex location on wake
                        J = CartesianIndex(I[1], js)

                        # induced velocity from wake on its own vertex
                        v_wake = induced_velocity(J, wakes[jsurf];
                            finite_core = wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
                            symmetric = symmetric[jsurf],
                            nc = nc_wake,
                            trailing_vortices = trailing_vortices[jsurf],
                            xhat = xhat)
                        wake_velocities[isurf][I] += v_wake
                        if isurf == 1 && I == CartesianIndex(2, 2)
                            @info "wake velocity contributions" jsurf same_surface v_surface v_wake total=wake_velocities[isurf][I]
                        end
                    end
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
function translate_wake(panel::WakePanel, wake_velocities, dt)

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
function translate_wake!(wake, wake_velocities, dt; nwake = size(wake, 1))

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
function shed_wake!(wake::AbstractMatrix, wake_shedding_locations,
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
function shed_wake!(wakes::AbstractVector{<:AbstractMatrix}, wake_shedding_locations,
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

#--- Wake Shedding Methods ---#

"""
    WakeSheddingMethod

Abstract supertype for particle-shedding strategies used by
[`PanelParticleWake`](@ref).
"""
abstract type WakeSheddingMethod end

"""No particle shedding."""
struct NoShed <: WakeSheddingMethod end

"""Gaussian particle shedding with prescribed smoothing width `sigma`."""
struct SigmaPPS{TF} <: WakeSheddingMethod
    sigma::TF
    p_per_step::Int
end

"""Overlap-based particle shedding with target overlap ratio."""
struct OverlapPPS{TF} <: WakeSheddingMethod
    overlap::TF
    p_per_step::Int
end

"""Gaussian particle shedding with fixed `sigma` and target `overlap`.
   The number of particles per step is computed as
   `ceil(overlap * dist / sigma)`.
"""
struct SigmaOverlap{TF} <: WakeSheddingMethod
    sigma::TF
    overlap::TF
end

function _shed_particles!(pfield, r1, r2, Γ, method::OverlapPPS)
    dist = norm(r2 - r1)
    dist < eps(typeof(dist)) && return zero(typeof(Γ))
    sigma = dist * method.overlap / method.p_per_step
    return _shed_particles!(pfield, r1, r2, Γ, SigmaPPS(sigma, method.p_per_step))
end

function _shed_particles!(pfield, r1, r2, Γ, method::SigmaPPS)
    sigma = method.sigma
    p_per_step = method.p_per_step
    distance_vector = (r2 - r1) / p_per_step
    Xp = r1 + distance_vector * 0.5
    Γp = Γ * distance_vector
    circ_scalar = Γ * norm(distance_vector)
    total_added = zero(typeof(circ_scalar))
    for _ in 1:p_per_step
        FLOWVPM.add_particle(pfield, Xp, Γp, sigma; circulation=circ_scalar)
        total_added += circ_scalar
        Xp += distance_vector
    end
    return total_added
end

function _shed_particles!(pfield, r1, r2, Γ, method::SigmaOverlap)
    dist = norm(r2 - r1)
    dist < eps(typeof(dist)) && return zero(typeof(Γ))

    pps = max(1, ceil(Int, method.overlap * dist / method.sigma))
    return _shed_particles!(pfield, r1, r2, Γ, SigmaPPS(method.sigma, pps))
end

function _shed_particles!(pfield, r1, r2, Γ, ::NoShed)
    return zero(typeof(Γ))
end

#--- Panel + Particle Wake ---#

"""
    PanelParticleWake(system; optargs...)

Hybrid wake model that combines VortexLattice's panel wake buffer with a
FLOWVPM particle field. The panel buffer holds the most recent `nwakerows` of
shed vorticity; once the buffer fills, the oldest row is converted to vortex
particles on each subsequent shed. The panel arrays are aliased from `system`
so existing VLM kernels continue to operate on the same storage.
"""
struct PanelParticleWake{TF, MT<:WakeSheddingMethod, MU<:WakeSheddingMethod, TPF, TFW, TBFW, TFWakeFMM, TVehicleFMM}
    wakes::Vector{Matrix{WakePanel{TF}}}
    wake_shedding_locations::Vector{Vector{SVector{3,TF}}}
    wake_velocities::Vector{Matrix{SVector{3,TF}}}
    nwake::Vector{Int}
    nwakerows::Int
    overflowed::Base.RefValue{Bool}
    pfield::TPF
    trailing_edge_filaments::TFW
    boundary_filaments::TBFW
    fmm_wake::TFWakeFMM
    fmm_vehicle::TVehicleFMM
    method_trailing::Vector{MT}
    method_unsteady::Vector{MU}
    prev_bottom_gamma::Vector{Vector{TF}}
    pending_overflow::Vector{Vector{WakePanel{TF}}}
    has_pending::Vector{Bool}
    eta::TF
    probes_active::ProbeSystem{TF}
end

function PanelParticleWake(system;
        nwakerows::Int=2,
        max_particles::Int=10_000,
        eta::Real=0.3,
        fmm::FLOWVPM.FMM=FLOWVPM.FMM(),
        fmm_wake::Union{Nothing, FLOWVPM.FMM}=nothing,
        fmm_vehicle::Union{Nothing, FLOWVPM.FMM}=nothing,
        method_trailing::WakeSheddingMethod=OverlapPPS(1.3, 2),
        method_unsteady::WakeSheddingMethod=OverlapPPS(1.3, 2),
        vpm_kwargs::NamedTuple=NamedTuple(),
    )

    TF = eltype(system.wake_shedding_locations[1][1])
    nsurf = length(system.surfaces)

    for i in 1:nsurf
        @assert size(system.wakes[i], 1) == nwakerows "System.wakes[$i] has " *
            "$(size(system.wakes[i], 1)) rows but PanelParticleWake expects " *
            "nwakerows=$nwakerows. Rebuild the System with nw=fill($nwakerows, nsurf)."
    end

    # Alias panel-wake storage from System (shared arrays — not copies).
    wakes = system.wakes
    wake_shedding_locations = system.wake_shedding_locations
    wake_velocities = system.V             # persistent wake-node velocity buffer
    nwake = zeros(Int, nsurf)              # active rows (grows 0 → nwakerows)
    overflowed = Ref(false)

    # Particle field (FLOWVPM)
    pfield = FLOWVPM.ParticleField(max_particles, TF;
        fmm=fmm, vpm_kwargs...)

    fmm_wake = something(fmm_wake, fmm)
    fmm_vehicle = something(fmm_vehicle, fmm)

    # Trailing-edge filament wrapper for FMM coupling
    trailing_edge_filaments = FilamentWrapper(system.wakes)

    # Boundary filament: top edge of the most-recently converted wake row.
    boundary_filaments = BoundaryFilamentWrapper(system.wakes)

    # Independent runtime FMM states for wake and vehicle coupling.
    fmm_wake = Base.RefValue{FLOWVPM.FMM}(fmm_wake)
    fmm_vehicle = Base.RefValue{FLOWVPM.FMM}(fmm_vehicle)

    # Per-surface shedding methods (plan spec)
    method_trailing_vec = [method_trailing for _ in 1:nsurf]
    method_unsteady_vec = [method_unsteady for _ in 1:nsurf]

    # Per-surface memory of the previously-converted oldest row's circulation
    # (one entry per spanwise panel). Used by _convert_to_particles! as Γ_tm1.
    prev_bottom_gamma = [zeros(TF, size(system.wakes[i], 2)) for i in 1:nsurf]

    # Overflow buffer: holds the oldest row saved before each shift so that
    # _convert_to_particles! can read it without needing an extra buffer row in
    # the system allocation. has_pending[i] is true iff pending_overflow[i] is
    # populated and ready to convert.
    pending_overflow = [Vector{WakePanel{TF}}(undef, size(system.wakes[i], 2)) for i in 1:nsurf]
    has_pending = fill(false, nsurf)

    # Reusable scratch buffer for the "active probes" subset passed to each FMM
    # call. Sized to the full persistent probe count up front so every later
    # `resize_active!` call (always to <= this capacity) reuses this storage
    # instead of allocating a fresh ProbeSystem.
    probes_active = ProbeSystem(length(system.probes.position), TF)

    return PanelParticleWake{TF, typeof(method_trailing), typeof(method_unsteady), typeof(pfield), typeof(trailing_edge_filaments), typeof(boundary_filaments), typeof(fmm_wake), typeof(fmm_vehicle)}(
        wakes, wake_shedding_locations, wake_velocities, nwake, nwakerows, overflowed, pfield,
        trailing_edge_filaments, boundary_filaments, fmm_wake, fmm_vehicle, method_trailing_vec, method_unsteady_vec,
        prev_bottom_gamma, pending_overflow, has_pending, TF(eta), probes_active,
    )
end

"""
    reset!(w::PanelParticleWake)

Zero the wake-node velocity buffer and reset the particle field velocity/Jacobian
and SFS properties. Preserves particle positions, strengths, and panel geometry.
"""
function reset!(w::PanelParticleWake)
    # zero wake-node velocities (panel wake)
    for V in w.wake_velocities
        fill!(V, zero(eltype(V)))
    end

    # reset particle velocity and Jacobian fields (preserve position and strength)
    FLOWVPM._reset_particles(w.pfield)

    # reset particle SFS properties
    FLOWVPM._reset_particles_sfs(w.pfield)

    return w
end

"""
    update_TE!(w::PanelParticleWake, system::System)

Snap the first row of wake panels to the current `wake_shedding_locations`,
preserving the bottom edge, core size, and circulation of each panel.

The freestream-based update of `wake_shedding_locations` itself lives in the
propagate step (see `update_wake_shedding_locations_unsteady!`).
"""
function update_TE!(w::PanelParticleWake, system)
    for isurf in eachindex(w.wakes)
        w.nwake[isurf] == 0 && continue   # no active wake rows yet
        wake = w.wakes[isurf]
        wsl = w.wake_shedding_locations[isurf]
        ns = length(wsl) - 1
        for j in 1:ns
            rtl = wsl[j]
            rtr = wsl[j+1]
            rbl = bottom_left(wake[1, j])
            rbr = bottom_right(wake[1, j])
            core_size = get_core_size(wake[1, j])
            gamma = circulation_strength(wake[1, j])
            wake[1, j] = WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
        end
    end
    return w
end

"""
    propagate!(w::PanelParticleWake, dt; scheme=EulerScheme(), Vinf=nothing, relax=true)

Convect the wake forward by `dt`: translate active wake panels by their stored
node velocities (`w.wake_velocities`) and advance the particle field using the
selected `IntegrationScheme`. If `Vinf` is provided it is seeded onto every
active particle via `apply_freestream!` before the integration step.
"""
function propagate!(w::PanelParticleWake, dt;
        scheme=EulerScheme(), Vinf=nothing, relax=true)
    # panel wake — translate nodes via stored wake_velocities (only active rows)
    for i_surf in eachindex(w.wakes)
        w.nwake[i_surf] == 0 && continue
        translate_wake!(w.wakes[i_surf], w.wake_velocities[i_surf], dt;
            nwake = w.nwake[i_surf])
    end

    # apply freestream to particles (once, pre-step) if requested
    Vinf === nothing || apply_freestream!(w, Vinf)

    # convect particles via selected scheme
    _integrate_particles!(w, dt, scheme; relax)

    return w
end

"""
    _convert_to_particles!(w::PanelParticleWake)

Convert the oldest active row of each surface's panel wake into vortex
particles in `w.pfield`. After conversion, `w.prev_bottom_gamma[isurf]` is
updated with the row's circulations so the next call can form `Γ - Γ_tm1`.
"""
function _convert_to_particles!(w::PanelParticleWake{TF}) where TF
    for isurf in eachindex(w.wakes)
        w.has_pending[isurf] || continue
        pending   = w.pending_overflow[isurf]
        ns        = length(pending)
        method_t  = w.method_trailing[isurf]
        method_u  = w.method_unsteady[isurf]
        prev_gamma = w.prev_bottom_gamma[isurf]

        r1_le_first = top_left(pending[1])
        rend_le     = top_right(pending[ns])
        wraps = norm(r1_le_first - rend_le) < 5 * eps(TF)
        Γ_last = wraps ? circulation_strength(pending[ns]) : zero(TF)

        # track how much scalar circulation particles contributed per panel
        added_by_particles = zeros(TF, ns)

        for j in 1:ns
            panel = pending[j]
            Γ     = circulation_strength(panel)
            r1_le = top_left(panel)
            r1_te = bottom_left(panel)
            r2_te = bottom_right(panel)

            added_t = _shed_particles!(w.pfield, r1_le, r1_te, Γ - Γ_last, method_t)

            Γ_tm1 = prev_gamma[j]
            added_u = _shed_particles!(w.pfield, r1_te, r2_te, Γ - Γ_tm1, method_u)

            # ensure numeric zero for any missing returns
            added_t === nothing && (added_t = zero(TF))
            added_u === nothing && (added_u = zero(TF))
            added_by_particles[j] = added_t + added_u

            Γ_last = Γ
            prev_gamma[j] = Γ
        end

        if !wraps
            panel = pending[ns]
            Γ     = circulation_strength(panel)
            r_le  = top_right(panel)
            r_te  = bottom_right(panel)
            _shed_particles!(w.pfield, r_le, r_te, -Γ, method_t)
        end

        # Boundary filament: top edge of the converted row.
        # The WakeBufferRings ring goes rtl→rbl→rbr→rtr (Ring B, CW from above),
        # so its top edge is rtr→rtl with +Γ = rtl→rtr with -Γ. The filament
        # stored as r1=rtl, r2=rtr therefore needs -Γ to match.
        bfw = w.boundary_filaments
        bfw.active[isurf] = true
        for j in 1:ns
            panel = pending[j]
            bfw.r1[isurf][j]        = top_left(panel)
            bfw.r2[isurf][j]        = top_right(panel)
            bfw.gamma[isurf][j]     = -circulation_strength(panel)
            bfw.core_size[isurf][j] = panel.core_size
        end

    end
    return w
end

"""
    shed_wake!(w::PanelParticleWake, system::System, dt, Gamma)

Advance the hybrid wake one step: convert the oldest row of any surface whose
buffer is full into particles, shed a new row of wake panels at the trailing
edge, and grow `nwake` until the buffer is saturated.
"""
function shed_wake!(w::PanelParticleWake, system, dt, Gamma)  # dt unused; kept for call-site compatibility
    iΓ = 0
    for isurf in eachindex(w.wakes)
        surface = system.surfaces[isurf]
        wake    = w.wakes[isurf]
        wsl     = w.wake_shedding_locations[isurf]
        nc, ns  = size(surface)
        ls      = LinearIndices((nc, ns))

        # save propagated row-1 top nodes before the shift overwrites them
        if w.nwake[isurf] > 0
            saved_rtl = [top_left(wake[1, j]) for j in 1:ns]
            saved_rtr = [top_right(wake[1, j]) for j in 1:ns]
        end

        # when the buffer is full, save the oldest row before it is overwritten
        if w.nwake[isurf] == w.nwakerows
            for j in 1:ns
                w.pending_overflow[isurf][j] = wake[w.nwakerows, j]
            end
            w.has_pending[isurf] = true
        else
            w.has_pending[isurf] = false
        end

        # shift only rows that fit — oldest row was already saved above
        nkeep = min(w.nwake[isurf], w.nwakerows - 1)
        for j = 1:ns, i = nkeep:-1:1
            wake[i+1, j] = wake[i, j]
        end

        for j = 1:ns
            rtl = wsl[j]
            rtr = wsl[j+1]
            if w.nwake[isurf] > 0
                rbl = saved_rtl[j]
                rbr = saved_rtr[j]
            else
                # first shed: mirror initial_wake_panels! — bottom corner is the
                # TE projected through wsl by the eta ratio, so kinematics in wsl
                # (from update_wake_shedding_locations_unsteady!) are included
                rte_l = bottom_left(surface[end, j])
                rte_r = bottom_right(surface[end, j])
                rbl = rte_l + (rtl - rte_l) / w.eta
                rbr = rte_r + (rtr - rte_r) / w.eta
            end
            core_size = get_core_size(surface[end, j])
            gamma = Gamma[iΓ + ls[end, j]]
            wake[1, j] = WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
        end

        iΓ += length(surface)

        if w.nwake[isurf] < w.nwakerows
            w.nwake[isurf] += 1
        end
    end

    _convert_to_particles!(w)

    if any(w.has_pending)
        w.overflowed[] = true
    end

    return w
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
