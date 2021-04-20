"""
    initialize_wake_shedding_locations!(wake_shedding_locations, surfaces)

Initialize the wake shedding locations for one or more surfaces based on the
provided surfaces.
"""
initialize_wake_shedding_locations!(wake_shedding_locations, surfaces)

# multiple surfaces
function initialize_wake_shedding_locations!(wake_shedding_locations, surfaces)
    for isurf = 1:length(surfaces)
        initialize_wake_shedding_locations!(wake_shedding_locations[isurf], surfaces[isurf])
    end
end

# single wake, wake panels input
function initialize_wake_shedding_locations!(wake_shedding_locations, surface::AbstractMatrix)
    ns = size(surface, 2)
    for j = 1:ns
        wake_shedding_locations[j] = bottom_left(surface[end,j])
    end
    wake_shedding_locations[end] = bottom_right(surface[end,end])
    return wake_shedding_locations
end

# single wake, grid input
function initialize_wake_shedding_locations!(wake_shedding_locations, surface::AbstractArray{<:Number, 3})
    for j = 1:size(surface, 3)
        wake_shedding_locations[j] = SVector(surface[1,end,j], surface[2,end,j], surface[3,end,j])
    end
    return wake_shedding_locations
end

"""
    get_wake_shedding_locations!(wake_shedding_locations, wakes)

Get the wake shedding locations for one or more surfaces based on the
provided wakes or wake vertices.
"""
get_wake_shedding_locations!(wake_shedding_locations, wakes)

# multiple surfaces
function get_wake_shedding_locations!(wake_shedding_locations, wakes)
    for isurf = 1:length(wakes)
        get_wake_shedding_locations!(wake_shedding_locations[isurf], wakes[isurf])
    end
end

# single wake, wake panels input
function get_wake_shedding_locations!(wake_shedding_locations, wake::AbstractMatrix{<:WakePanel})
    ns = size(wake, 2)
    for j = 1:ns
        wake_shedding_locations[j] = top_left(wake[1,j])
    end
    wake_shedding_location[j] = top_right(wake[1,end])
    return wake_shedding_location
end

# single wake, vertex inputs
function get_wake_shedding_locations!(wake_shedding_locations, vertices::AbstractMatrix{<:SVector{3, <:Any}})
    for j = 1:size(vertices, 2)
        wake_shedding_locations[j] = vertices[1,j]
    end
    return wake_shedding_locations
end

"""
    wake_velocities!(Vw, surfaces, repeated_points,
        wake_shedding_locations, wakes, reference, freestream, Γ; Vte, symmetric,
        surface_id, wake_finite_core, iwake, trailing_vortices, xhat)

Calculate the velocities at the wake verticies.

Store the result in `Vw` if provided, otherwise store in `system`.
"""
function wake_velocities!(Vw, surfaces, repeated_points,
    wake_shedding_locations, wakes, ref, fs, Γ; Vte, symmetric, surface_id,
    wake_finite_core, iwake, trailing_vortices, xhat)
    # number of surfaces
    nsurf = length(surfaces)
    # additional velocity function, wrapped to produce static vector of known type
    additional_velocity = (r) -> SVector{3, eltype(eltype(eltype(Vw)))}(fs.additional_velocity(r[1], r[2], r[3]))
    # loop through all surfaces
    for isurf = 1:nsurf
        # number of chordwise and spanwise panels
        nc, ns = size(surfaces[isurf])
        nw = size(wakes[isurf], 1)
        # velocities at wake shedding locations
        for is = 1:ns+1
            # check if this point is a duplicate, skip if it is
            if (isurf, is) in keys(repeated_points)
                for (jsurf, js) in repeated_points[(isurf, is)]
                    # NOTE: we assume that a point is not repeated on the same surface
                    if jsurf < isurf
                        Vw[isurf][1,is] = Vw[jsurf][1, js]
                        continue
                    end
                end
            end
            # get vertex location
            rc = wake_shedding_locations[isurf][is]
            # freestream velocity
            Vw[isurf][1,is] = freestream_velocity(fs)
            # rotational velocity
            Vw[isurf][1,is] += rotational_velocity(rc, fs, ref)
            # additional velocity field
            Vw[isurf][1,is] += additional_velocity(rc)
        end
        # velocities at the vertices of wake panels that have been shed
        cr = CartesianIndices((2:iwake[isurf]+1, 1:ns+1))
        for I in cr
            # check if this point is a duplicate, skip if it is
            if (isurf, I[2]) in keys(repeated_points)
                for (jsurf, js) in repeated_points[(isurf, I[2])]
                    # NOTE: we assume that a point is not repeated on the same surface
                    if jsurf < isurf
                        Vw[isurf][I] = Vw[jsurf][I[1], js]
                        continue
                    end
                end
            end
            # get vertex location
            if I[1] <= iwake[isurf] && I[2] <= ns
                rc = top_left(wakes[isurf][I[1], I[2]])
            elseif I[1] == iwake[isurf] + 1 && I[2] <= ns
                rc = bottom_left(wakes[isurf][I[1]-1, I[2]])
            elseif I[1] <= iwake[isurf] && I[2] == ns + 1
                rc = top_right(wakes[isurf][I[1], I[2]-1])
            else # I[1] == nw + 1 && I[2] == ns + 1
                rc = bottom_right(wakes[isurf][I[1]-1, I[2]-1])
            end
            # freestream velocity
            Vw[isurf][I] = freestream_velocity(fs)
            # rotational velocity
            Vw[isurf][I] += rotational_velocity(rc, fs, ref)
            # additional velocity field
            Vw[isurf][I] += additional_velocity(rc)
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
                Vw[isurf][I] += induced_velocity(rc, surfaces[jsurf], vΓ;
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
                    Vw[isurf][I] += induced_velocity(J, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
                        nc = iwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                else
                    # induced velocity from wake on another wake's vertex
                    Vw[isurf][I] += induced_velocity(rc, wakes[jsurf];
                        finite_core = wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
                        symmetric = symmetric[jsurf],
                        nc = iwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                end

                jΓ += Ns # increment Γ index
            end
        end
        # velocity at the vertices of wake panels that have not yet been shed
        cr = CartesianIndices((iwake[isurf]+2:nw+1, 1:ns+1))
        for I in cr
            Vw[isurf][I] = -Vte[isurf][I[2]]
        end
    end

    return Vw
end

function shed_wake!(wakes, surfaces)
    # set circulation of the newly shed panel
    iΓs = 0 # index for accessing surface circulation
    iΓw = 0 # index for accessing wake circulation
    for isurf = 1:nsurf
        Ns = length(system.surfaces[isurf])
        Nw = length(system.wakes[isurf])
        # get view of circulations for this surface/wake
        vΓs = reshape(view(u_Γ, iΓ+1:iΓ+Ns[isurf]), nc[isurf], ns[isurf])
        vΓw = reshape(view(u_Γw, iΓ+1:iΓ+Nw[isurf]), nw[isurf], ns[isurf])
        # replace storage for last panel with new panel
        copyto!(view(vΓw, nw[isurf], :), view(vΓs, nc[isurf], :))
        # shift wake panel data downstream
        rowshift!(vΓw)
        # move indices to next surface
        iΓs += Ns[isurf]
        iΓw += Nw[isurf]
    end
    # reset the wake shedding location
    for isurf = 1:nsurf
        for j = 1:ns[isurf]+1
            # extract trailing edge coordinate
            if j < ns + 1
                rte = bottom_left(system.surfaces[isurf][end, j])
            else
                rte = bottom_right(system.surfaces[isurf][end, j-1])
            end
            # replace storage for coordinates of last wake panel
            copyto!(view(u_ζw, 1:3, nw[isurf], :), rte)
        end
        # shift wake panel coordinates downstream
        for i = 1:3
            rowshift!(view(u_ζw, i, 1:nw[isurf], :))
        end
    end
end

# """
#     translate_wake(panel, Vw, dt)
#
# Return a translated copy of the wake panel `panel` given the wake corner velocities
# `Vw` and the time step `dt`
#
# # Arguments
#  - `panel`: Wake panel (of type [`WakePanel`](@ref))
#  - `Vw`: Matrix containing the velocities at each of the four corners
#     of `panel`
#  - `dt`: Time step (seconds)
# """
# @inline function translate_wake(panel::WakePanel, Vw, dt)
#
#     # extract corners
#     rtl = top_left(panel)
#     rtr = top_right(panel)
#     rbl = bottom_left(panel)
#     rbr = bottom_right(panel)
#
#     # get vortex filament length
#     lt = norm(rtr - rtl)
#     lb = norm(rbl - rbr)
#     ll = norm(rtl - rbl)
#     lr = norm(rbr - rtr)
#     l1 = lt + lb + ll + lr
#
#     # translate corners
#     rtl += Vw[1,1]*dt
#     rtr += Vw[1,2]*dt
#     rbl += Vw[2,1]*dt
#     rbr += Vw[2,2]*dt
#
#     # get new vortex filament length
#     lt = norm(rtr - rtl)
#     lb = norm(rbl - rbr)
#     ll = norm(rtl - rbl)
#     lr = norm(rbr - rtr)
#     l2 = lt + lb + ll + lr
#
#     # use previous core size
#     core_size = get_core_size(panel)
#
#     # correct vorticity for vortex stretching
#     gamma = circulation_strength(panel)*l1/l2
#
#     return WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
# end
#
# """
#     translate_wake!(wake, Vw, dt; nwake = size(wake, 1))
#
# Translate the wake panels in `wake` given the corner velocities `Vw`
# and the time step `dt`.
#
# # Arguments
#  - `wake`: Matrix of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns)
#     where `nw` is the number of chordwise wake panels and `ns` is the number of
#     spanwise panels, defaults to no wake panels
#  - `Vw`: Velocities at each of the vertices corresponding to the
#     wake panels in `wake`
#  - `dt`: Time step
#
# # Keyword Arguments
#  - `nwake`: Number of chordwise wake panels to use from `wake`, defaults to all
#     provided wake panels
# """
# @inline function translate_wake!(wake, Vw, dt; nwake = size(wake, 1))
#
#     nw = nwake
#     ns = size(wake, 2)
#     cw = CartesianIndices((nw, ns))
#
#     for I in cw
#
#         panel = wake[I]
#
#         vV = view(Vw, I[1]:I[1]+1, I[2]:I[2]+1)
#
#         wake[I] = translate_wake(panel, vV, dt)
#     end
#
#     return wake
# end
#
# """
#     shed_wake!([wakes,] system)
#
# Shed a new wake panel from the wake shedding locations and translate existing
# wake panels.
#
# Store the result in `wakes` if provided, otherwise store in `system`.
# """
# shed_wake!
#
# # unsteady, system input and output
# shed_wake!(system::AbstractSystem) = shed_wake!(system.wakes, system)
#
# # unsteady, system input
# function shed_wake!(wakes, system;
#     wake_shedding_locations = system.wake_shedding_locations,
#     Vw = system.Vw,
#     surfaces = system.surfaces,
#     Gamma = system.Gamma,
#     nwake = system.nwake,
#     dt = system.dt)
#
#     return shed_wake!(wakes, wake_shedding_locations, Vw,
#         surfaces, Gamma, nwake, dt)
# end
#
# # unsteady, multiple surfaces
# @inline function shed_wake!(wakes::AbstractVector{<:AbstractMatrix}, wake_shedding_locations,
#     Vw, surfaces, Γ, nwake, dt)
#
#     iΓ = 0
#     for i = 1:length(surfaces)
#
#         N = length(surfaces[i])
#
#         vΓ = view(Γ, iΓ+1:iΓ+N)
#
#         shed_wake!(wakes[i], wake_shedding_locations[i], Vw[i],
#             surfaces[i], vΓ, nwake[i], dt)
#
#         iΓ += N
#     end
#
#     return wakes
# end
#
# # unsteady, single surface
# function shed_wake!(wake::AbstractMatrix, wake_shedding_locations,
#     Vw, surface, Γ, nwake, dt)
#
#     nc, ns = size(surface)
#     nw = size(wake, 1)
#     ls = LinearIndices((nc, ns))
#
#     # replace the last chordwise panels with the newly shed wake panels
#     for j = 1:ns
#
#         # shedding location coordinates
#         rtl = wake_shedding_locations[j]
#         rtr = wake_shedding_locations[j+1]
#
#         # shed coordinates
#         rbl = rtl + Vw[1, j]*dt
#         rbr = rtr + Vw[1, j+1]*dt
#
#         # use core size from the shedding panel
#         core_size = get_core_size(surface[end, j])
#
#         # use circulation strength from the shedding panel
#         gamma = Γ[ls[end,j]]
#
#         # replace the oldest wake panel
#         wake[end,j] = WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
#     end
#
#     # translate all existing wake panels except the most recently shed panels
#     translate_wake!(wake, Vw, dt, nwake = min(nwake, nw-1))
#
#     # shift wake panels to make newly shed wake panel first
#     rowshift!(wake)
#
#     return wake
# end
#
