"""
    get_wake_velocities!(wake_velocities, surface, wake, ref, fs; symmetric,
        trailing_vortices, xhat, nwake)

Returns the velocities at the corners of the wake panels in `wake`
"""
@inline function get_wake_velocities!(wake_velocities, surface::AbstractMatrix,
    wake::AbstractMatrix, ref, fs, Γ; symmetric, surface_id, wake_id,
    trailing_vortices, xhat, nwake, repeated_points)

    get_wake_velocities!([wake_velocities], [surface], [wake], ref, fs, Γ;
        symmetric = [symmetric],
        surface_id = surface_id,
        wake_id = wake_id,
        trailing_vortices = [trailing_vortices],
        xhat = xhat,
        nwake = [nwake],
        repeated_points = repeated_points)

    return wake_velocities
end

"""
    get_wake_velocities!(wake_velocities, surfaces, wakes, ref, fs; symmetric,
        surface_id, wake_id, trailing_vortices, xhat, nwake, repeated_points)

Returns the velocities at the corners of the wake panels in `wakes`
"""
@inline function get_wake_velocities!(wake_velocities, surfaces::AbstractVector{<:AbstractMatrix},
    wakes::AbstractVector{<:AbstractMatrix}, ref, fs, Γ; symmetric,
    surface_id, wake_id, trailing_vortices, xhat, nwake, repeated_points)

    nsurf = length(surfaces)

    for isurf = 1:nsurf

        nc, ns = size(surfaces[isurf])
        nw = nwake[isurf]

        # velocity at the trailing edge
        for is = 1:ns+1

            # check if this point is repeated, skip if it is
            if (isurf, is) in keys(repeated_points)
                for (jsurf, js) in repeated_points[(isurf, is)]
                    # NOTE: we assume that a point is not repeated on the same surface
                    if jsurf < isurf
                        wake_velocities[isurf][1,is] = wake_velocities[jsurf][1,js]
                        continue
                    end
                end
            end

            # wake corner location
            if is <= ns
                rc = bottom_left(surfaces[isurf][end,is])
            else
                rc = bottom_right(surfaces[isurf][end,is-1])
            end

            # external velocity at the wake corner location
            wake_velocities[isurf][1,is] = external_velocity(fs, rc, ref.r)

            # induced velocity at the wake corner location
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                Ns = length(surfaces[jsurf])

                # check if receiving point is repeated on the sending surface
                if isurf == jsurf
                    same_surface = true
                    js = is
                else
                    if (isurf, is) in keys(repeated_points)
                        # check if point is on the sending surface
                        idx = findfirst(x -> x[1] == jsurf, repeated_points[(isurf, is)])

                        if isnothing(idx)
                            same_surface = false
                        else
                            same_surface = true
                            js = repeated_points[(isurf, is)][idx][2]
                        end
                    else
                        same_surface = false
                    end
                end

                # check if finite core model is enabled
                finite_core = surface_id[jsurf] < 0 || surface_id[isurf] != surface_id[jsurf]

                # add the induced velocity from the surface, ignoring co-located
                # bound and trailing vortices
                vΓ = view(Γ, jΓ+1:jΓ+Ns)

                if same_surface
                    wake_velocities[isurf][1,is] += induced_velocity(js, surfaces[jsurf], vΓ;
                        symmetric = symmetric[jsurf],
                        finite_core = finite_core,
                        trailing_vortices = false)
                else
                    wake_velocities[isurf][1,is] += induced_velocity(rc, surfaces[jsurf], vΓ;
                        symmetric = symmetric[jsurf],
                        finite_core = finite_core,
                        trailing_vortices = false)
                end

                # add induced velocity from wake, ignoring co-located bound and
                # trailing vortices
                if same_surface
                    wake_velocities[isurf][1,is] += induced_velocity(
                        CartesianIndex(1, js), wakes[jsurf];
                        symmetric = symmetric[jsurf],
                        finite_core = wake_id[jsurf] < 0 || surface_id[isurf] != wake_id[jsurf],
                        nwake = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)

                else
                    wake_velocities[isurf][1,is] += induced_velocity(rc, wakes[jsurf];
                        symmetric = symmetric[jsurf],
                        finite_core = wake_id[jsurf] < 0 || surface_id[isurf] != wake_id[jsurf],
                        nwake = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end

                jΓ += Ns # increment Γ index
            end
        end

        # velocity at all other wake corners
        cr = CartesianIndices((2:nw+1, 1:ns+1))

        for I in cr

            # check if this point is repeated, skip if it is
            if (isurf, I[2]) in keys(repeated_points)
                for (jsurf, js) in repeated_points[(isurf, I[2])]
                    # NOTE: we assume that a point is not repeated on the same surface
                    if jsurf < isurf
                        wake_velocities[isurf][I] = wake_velocities[jsurf][I[1], js]
                        continue
                    end
                end
            end

            # extract relevant corner
            if I[1] <= nw && I[2] <= ns
                rc = top_left(wakes[isurf][I[1], I[2]])
            elseif I[1] == nw + 1 && I[2] <= ns
                rc = bottom_left(wakes[isurf][I[1]-1, I[2]])
            elseif I[1] <= nw && I[2] == ns + 1
                rc = top_right(wakes[isurf][I[1], I[2]-1])
            else # I[1] == nw + 1 && I[2] == ns + 1
                rc = bottom_right(wakes[isurf][I[1]-1, I[2]-1])
            end

            # external velocity
            wake_velocities[isurf][I] = external_velocity(fs, rc, ref.r)

            # induced velocity from each surface and wake
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                Ns = length(surfaces[jsurf])

                # check if receiving point is repeated on the sending surface
                if isurf == jsurf
                    same_surface = true
                    js = I[2]
                else
                    if (isurf, I[2]) in keys(repeated_points)
                        # check if point I[2] on the sending surface
                        idx = findfirst(x -> x[1] == jsurf, repeated_points[(isurf, I[2])])

                        if isnothing(idx)
                            same_surface = false
                        else
                            same_surface = true
                            js = repeated_points[(isurf, I[2])][idx][2]
                        end
                    else
                        same_surface = false
                    end
                end

                # add induced velocity from surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)
                wake_velocities[isurf][I] += induced_velocity(rc, surfaces[jsurf], vΓ;
                    symmetric = symmetric[jsurf],
                    finite_core = surface_id[jsurf] < 0 || surface_id[isurf] != surface_id[jsurf],
                    trailing_vortices = false)

                # add induced velocity from wake
                if same_surface
                    J = CartesianIndex(I[1], js)

                    tmp = induced_velocity(J, wakes[jsurf];
                        symmetric = symmetric[jsurf],
                        finite_core = wake_id[jsurf] < 0 || surface_id[isurf] != wake_id[jsurf],
                        nwake = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)

                    if any(isnan.(tmp))
                        error()
                    end

                    wake_velocities[isurf][I] += tmp
                else
                    tmp = induced_velocity(rc, wakes[jsurf];
                        symmetric = symmetric[jsurf],
                        finite_core = wake_id[jsurf] < 0 || surface_id[isurf] != wake_id[jsurf],
                        nwake = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)

                    if any(isnan.(tmp))
                        error()
                    end

                    wake_velocities[isurf][I] += tmp
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
"""
@inline function translate_wake(panel::Wake, wake_velocities, dt)

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
    gamma = circulation_strength(panel)*l2/l1

    return Wake(rtl, rtr, rbl, rbr, core_size, gamma)
end

"""
    translate_wake!(wake, wake_velocities, dt; nwake = size(wake, 1))

Translate the wake panels in `wake` given a the corner velocities `wake_velocities`
and the time step `dt`
"""
@inline function translate_wake!(wake, wake_velocities, dt; nwake = size(wake, 1))

    nw = nwake
    ns= size(wake, 2)
    cw = CartesianIndices((nw, ns))

    for I in cw

        panel = wake[I]

        vV = view(wake_velocities, I[1]:I[1]+1, I[2]:I[2]+1)

        wake[I] = translate_wake(panel, vV, dt)
    end

    return wake
end


"""
    shed_wake!(wake, wake_velocities, dt, surface, Γ; nwake)

Sheds the wake from a surface, given the corner velocities of the wake panels `wake_velocities`,
the time step `dt`, and the circulation strength at the trailing edge `Γ_te`.

The wake panels in `wake` are shifted chordwise to make room for the newly shed
wake panel.
"""
shed_wake!

# single surface and wake
@inline function shed_wake!(wake::AbstractMatrix, wake_velocities, dt, surface, Γ; nwake)

    nc, ns = size(surface)
    nw = size(wake, 1)
    ls = LinearIndices((nc, ns))

    # replace the last chordwise panels with the newly shed wake panels
    for j = 1:ns
        # trailing edge coordinates
        rtl = bottom_left(surface[end, j])
        rtr = bottom_right(surface[end, j])

        # shed coordinates
        rbl = rtl + wake_velocities[1, j]*dt
        rbr = rtr + wake_velocities[1, j+1]*dt

        # use core size from the trailing edge
        core_size = get_core_size(surface[end, j])

        # use circulation strength from the trailing edge
        gamma = Γ[ls[end,j]]

        # replace the oldest wake panel
        wake[end,j] = Wake(rtl, rtr, rbl, rbr, core_size, gamma)
    end

    # translate all wake panels except the newly shed panels
    translate_wake!(wake, wake_velocities, dt; nwake = nwake - 1)

    # shift wake panels to make the newly shed panels first
    rowshift!(wake)

    return wake
end

# multiple surfaces and wakes
@inline function shed_wake!(wake::AbstractVector{<:AbstractMatrix}, wake_velocities,
    dt, surfaces, Γ; nwake)

    iΓ = 0
    for i = 1:length(surfaces)

        N = length(surfaces[i])

        vΓ = view(Γ, iΓ+1:iΓ+N)

        shed_wake!(wake[i], wake_velocities[i], dt, surfaces[i], vΓ; nwake=nwake[i])
    end

    return wake
end

"""
    rowshift!(A)

Circularly shifts data down one row.
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
