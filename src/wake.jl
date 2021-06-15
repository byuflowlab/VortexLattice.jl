"""
    trailing_edge_velocities!(V, surfaces, reference, freestream, additional_velocity)

Calculate freestream velocities at the trailing edge vertices.  Neglect induced
velocities.
"""
@inline function trailing_edge_velocities!(V, surfaces, reference, freestream,
    additional_velocity)

    # number of surfaces
    nsurf = length(surfaces)

    # loop through all surfaces
    for isurf = 1:nsurf

        # number of spanwise panels
        ns = size(surfaces, 2)

        # loop through each trailing edge vertex
        for j = 1:ns+1

            # extract trailing edge coordinate
            if j < ns + 1
                rte = bottom_left(surfaces[isurf][end, j])
            else
                rte = bottom_right(surfaces[isurf][end, j-1])
            end

            # freestream velocity
            V[isurf][j] = freestream_velocity(freestream)

            # rotational velocity
            V[isurf][j] += rotational_velocity(rte, freestream, reference)

            # additional velocity
            V[isurf][j] += SVector{3}(additional_velocity(rte))

            # neglect induced velocity

        end

    end

    return V
end

"""
    wake_shedding_locations!(wake_shedding_locations, surfaces, reference,
        freestream, additional_velocity, t1, t2)

Calculate a new set of wake shedding locations.
"""
@inline function wake_shedding_locations!(wake_shedding_locations, surfaces,
    reference, freestream, additional_velocity, t1, t2)

    # calculate trailing edge freestream velocities, store in `wake_shedding_locations`
    trailing_edge_velocities!(wake_shedding_locations, surfaces, reference,
        freestream, additional_velocity)

    # time step size from t1 to teta
    dt = t2 - t1

    # calculate new wake shedding locations
    nsurf = length(surfaces)
    for isurf = 1:nsurf
        ns = size(surfaces[isurf], 2)
        for j = 1:ns+1
            # extract trailing edge coordinate
            if j < ns + 1
                rte = bottom_left(surfaces[isurf][end, j])
            else
                rte = bottom_right(surfaces[isurf][end, j-1])
            end
            # extract freestream velocity at this vertex
            V = wake_shedding_locations[isurf][j]
            # calculate new shedding location
            wake_shedding_locations[isurf][j] = rte + V*dt
        end
    end

    return wake_shedding_locations
end

"""
    free_wake_vertices!(wakes, wake_shedding_locations, surfaces,
        repeated_points, ref, fs, AIC, w, Γ; symmetric, surface_id,
        wake_finite_core, iwake, trailing_vortices, xhat)

Calculate the force-free wake vertex locations for steady-state operating
conditions and store in `wakes`.
"""
function free_wake_vertices!(wakes, wake_shedding_locations, surfaces,
    repeated_points, ref, fs, AIC, w, Γ; symmetric, surface_id,
    wake_finite_core, iwake, trailing_vortices, xhat, dt)

    # return early if no wake panels exist
    if !any(iw -> iw > 0, system.iwake)
        return wakes
    end

    # get total number of vertices on each surface
    nv = (iwake .+ 1) .* (size.(wakes, 2) .+ 1)

    # indices for accessing state variable vector
    iu2 = cumsum(nv)
    iu1 = vcat(1, iu2[1:end-1] .- 1)
    iu = [iu1[isurf]:iu2[isurf] for isurf = 1:nsurf]

    # freestream velocity
    Vinf = freestream_velocity(fs)

    # construct initial guess for wake vertices
    wake_vertices = zeros(eltype(eltype(eltype(wakes))), 3, sum(nv))
    for isurf = 1:nsurf
        # wake dimensions
        iw = iwake[isurf] # current number of chordwise wake panels
        ns = size(wakes[isurf], 2) # number of spanwise wake panels
        # view of current surface wake vertices
        vertices = reshape(view(wake_vertices, iu[isurf]), 3, iw + 1, ns + 1)
        # set initial guess for wake vertices
        for j = 1:ns+1
            # wake shedding location
            rshed = wake_shedding_locations[isurf][j]
            # vertex location is wake shedding location + freestream velocity * dt
            for i = 1:nw+1
                vertices[:, i, j] = rshed + Vinf * (i-1)*dt
            end
        end
    end

    # initial state variables are wake vertices
    u0 = wake_vertices

    # construct residual function
    f! = (r, u) -> free_wake_residuals!(r, u, iu; AIC, w, Γ, surfaces,
        repeated_points, wake_shedding_locations, wakes, reference, freestream,
        symmetric, surface_id, wake_finite_core, iwake, trailing_vortices, xhat,
        dt)

    # solve residual function for wake vertices
    result = nlsolve(f!, u0; autodiff = :forward)

    # save calculated wake vertices
    set_wake_vertices!(wakes, result.zero)

    return wakes
end

function free_wake_residuals!(r::AbstractMatrix{TF}, u::AbstractMatrix{TF}, iu;
    AIC, w, Γ, surfaces, repeated_points, wake_shedding_locations, wakes,
    reference, freestream, symmetric, surface_id, wake_finite_core, iwake,
    trailing_vortices, xhat, dt) where TF

    # TODO: Use collocation method instead of explicit euler.

    # initialize new versions of internal variables (if using autodiff)
    if TF <: ForwardDiff.Dual
        AIC = TF.(AIC)
        Γ = TF.(Γ)
        wakes = [WakePanel{TF}.(wake_panels) for wake_panels in wakes]
    end

    # reinterpret state and residual vector into better format
    u_ζw = [reinterpret(
        reshape, # reshape to create matrix for each surface
        SVector{3, eltype(u)}, # each matrix entry is a static vector
        # wake vertices for this surfaces as array of shape (3, nc, ns)
        reshape(
            view(u, :, iu[isurf]), # portion of residual vector for this surface
            3, # number of dimensions
            size(wakes[isurf], 1) + 1, # number of chordwise points
            size(wakes[isurf], 2) + 1 # number of spanwise points
            )
        ) for isurf = 1:length(surfaces)]
    r_ζw = [reinterpret(
        reshape, # reshape to create matrix for each surface
        SVector{3, eltype(r)}, # each matrix entry is a static vector
        # wake vertices for this surfaces as array of shape (3, nc, ns)
        reshape(
            view(r, :, iu[isurf]), # portion of residual vector for this surface
            3, # number of dimensions
            size(wakes[isurf], 1) + 1, # number of chordwise points
            size(wakes[isurf], 2) + 1 # number of spanwise points
            )
        ) for isurf = 1:length(surfaces)]

    # set wake panel vertices to provided coordinates
    set_wake_vertices!(wakes, u_ζw)

    # update trailing edge coefficients to account for the new wake locations
    trailing_coefficients!(AIC, surfaces, symmetric, surface_id,
        trailing_vortices, xhat, wake_shedding_locations, wakes,
        wake_finite_core, iwake)

    # update surface circulation
    circulation!(Γ, AIC, w)

    # calculate wake velocities, store in residual vector
    wake_velocities!(r_ζw, surfaces, repeated_points, wake_shedding_locations,
        wakes, reference, freestream, additional_velocity, Γ; symmetric,
        surface_id, wake_finite_core, iwake, trailing_vortices, xhat)

    # calculate residuals
    for isurf = 1:nsurf
        # wake dimensions
        iw = iwake[isurf]
        ns = size(surfaces[isurf], 2)
        # set current surface wake vertices
        for j = 1:ns+1
            # first set of wake vertices should be at wake shedding location
            r_ζw[isurf][1,j] = u_ζw[isurf][1,j] .- system.wake_shedding_locations[isurf][j]
            # remaining wake vertices depend on position and velocity of previous vertex
            for i = 2:iw+1
                r_ζw[isurf][i,j] = u_ζw[isurf][i-1,j] + r_ζw[isurf][i-1,j]*dt
            end
        end
    end

    return r
end

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
    wake_shedding_locations[end] = top_right(wake[1,end])
    return wake_shedding_locations
end

# single wake, vertex inputs
function get_wake_shedding_locations!(wake_shedding_locations, vertices::AbstractMatrix{<:SVector{3, <:Any}})
    for j = 1:size(vertices, 2)
        wake_shedding_locations[j] = vertices[1,j]
    end
    return wake_shedding_locations
end

"""
    set_wake_shedding_locations!(wake_shedding_locations, wakes)

Set the wake shedding locations for one or more surfaces based on the
provided wakes or wake vertices.
"""
set_wake_shedding_locations!(wake_shedding_locations, wakes) =
    get_wake_shedding_locations!(wake_shedding_locations, wakes)

"""
    wake_velocities!(Vw, surfaces, repeated_points, wake_shedding_locations,
        wakes, reference, freestream, additional_velocity, Γ; symmetric,
        surface_id, wake_finite_core, iwake, trailing_vortices, xhat)

Calculate the velocities at the wake verticies.
"""
function wake_velocities!(Vw, surfaces, repeated_points, wake_shedding_locations,
    wakes, reference, freestream, additional_velocity, Γ; symmetric, surface_id,
    wake_finite_core, iwake, trailing_vortices, xhat)

    # number of surfaces
    nsurf = length(surfaces)

    # loop through all receiving wakes
    for isurf = 1:nsurf

        # current wake velocities
        receiving = Vw[isurf]

        # number of chordwise and spanwise wake panels
        ns = size(surfaces[isurf], 2)
        nw = iwake[isurf]

        # number of vertices on this wake
        Nr = (nw+1)*(ns+1)

        # cartesian indices of this wake's vertices
        cr = CartesianIndices((1:nw+1, 1:ns+1))

        # loop through each receiving vertex
        for I in cr

            # check if this point is a duplicate, skip if it is
            is = I[2]
            if (isurf, is) in keys(repeated_points)
                repeats = repeated_points[(isurf, is)]
                for (jsurf, js) in repeats
                    # NOTE: we assume that a point is not repeated on the same surface
                    if jsurf < isurf
                        Vw[isurf][I] = Vw[jsurf][I[1],js]
                        continue
                    end
                end
            end

            # vertex location
            if iwake[isurf] <= 0
                rc = wake_shedding_locations[isurf][I[2]]
            else
                if I[1] <= iwake[isurf] && I[2] <= ns
                    rc = top_left(wakes[isurf][I[1], I[2]])
                elseif I[1] == iwake[isurf] + 1 && I[2] <= ns
                    rc = bottom_left(wakes[isurf][I[1]-1, I[2]])
                elseif I[1] <= iwake[isurf] && I[2] == ns + 1
                    rc = top_right(wakes[isurf][I[1], I[2]-1])
                else # I[1] == nw + 1 && I[2] == ns + 1
                    rc = bottom_right(wakes[isurf][I[1]-1, I[2]-1])
                end
            end

            # freestream velocity
            Vw[isurf][I] = freestream_velocity(freestream)

            # rotational velocity
            Vw[isurf][I] += rotational_velocity(rc, freestream, reference)

            # additional velocity field
            Vw[isurf][I] += SVector{3}(additional_velocity(rc))

            # index for accessing surface circulation strength
            jΓ = 0

            # add induced velocity from each surface and wake
            for jsurf = 1:nsurf

                # current sending surface
                sending = surfaces[jsurf]

                # number of panels on this surface
                Ns = length(sending)

                # check if receiving point is repeated on the sending surface
                if isurf == jsurf
                    same_surface = true # the surfaces are the same
                    js = I[2] # vertex spanwise coordinate
                else
                    # check if this point is a duplicate
                    is = I[2]
                    if (isurf, is) in keys(repeated_points)
                        # the vertex is duplicated
                        repeats = repeated_points[(isurf, is)]
                        idx = findfirst(repeat -> repeat[1] == jsurf, repeats)
                        if isnothing(idx)
                            # the vertex is not duplicated on the sending surface
                            same_surface = false
                        else
                            # the vertex is duplicated on the sending surface
                            same_surface = true
                            js = repeated_points[(isurf, is)][idx][2]
                        end
                    else # the vertex is not duplicated
                        same_surface = false
                    end
                end

                # see if wake panels are being used
                wake_panels = iwake[jsurf] > 0

                # extract circulation values corresonding to the sending surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)

                # induced velocity from the sending surface
                # if same_surface && I[1] == 1
                #     Vw[isurf][I] += induced_velocity(js, surfaces[jsurf], vΓ;
                #         finite_core = surface_id[isurf] != surface_id[jsurf],
                #         wake_shedding_locations = wake_shedding_locations[jsurf],
                #         symmetric = symmetric[jsurf],
                #         trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                #         xhat = xhat)
                # else
                #     Vw[isurf][I] += induced_velocity(rc, surfaces[jsurf], vΓ;
                #         finite_core = surface_id[isurf] != surface_id[jsurf],
                #         wake_shedding_locations = wake_shedding_locations[jsurf],
                #         symmetric = symmetric[jsurf],
                #         trailing_vortices = trailing_vortices[jsurf] && !wake_panels,
                #         xhat = xhat)
                # end

                # induced velocity from the corresponding wake
                if wake_panels
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
                end

                jΓ += Ns # increment Γ index
            end
        end
    end

    return Vw
end

"""
    translate_wake(panel, Vw, dt)

Return a translated copy of the wake panel `panel` given the wake corner velocities
`Vw` and the time step `dt`
"""
@inline function translate_wake(panel::WakePanel, Vw, dt)

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
    rtl += Vw[1,1]*dt
    rtr += Vw[1,2]*dt
    rbl += Vw[2,1]*dt
    rbr += Vw[2,2]*dt

    # get new vortex filament length
    lt = norm(rtr - rtl)
    lb = norm(rbl - rbr)
    ll = norm(rtl - rbl)
    lr = norm(rbr - rtr)
    l2 = lt + lb + ll + lr

    # use previous core size
    core_size = get_core_size(panel)

    # correct vorticity for vortex stretching
    gamma = get_circulation(panel)*l1/l2

    return WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
end

"""
    translate_wake!(wake, Vw, dt; iwake)

Translate the wake panels in `wake` given the corner velocities `Vw`
and the time step `dt`.
"""
@inline function translate_wake!(wake::AbstractMatrix, wake_shedding_locations,
    Vw, dt; iwake)

    nw = iwake
    ns = size(wake, 2)
    cw = CartesianIndices((nw, ns))

    # translate wake shedding locations
    for j = 1:ns+1
        wake_shedding_locations[j] = wake_shedding_locations[j] + Vw[1,j]*dt
    end

    # translate wake panels
    for I in cw

        panel = wake[I]

        vV = view(Vw, I[1]:I[1]+1, I[2]:I[2]+1)

        wake[I] = translate_wake(panel, vV, dt)
    end

    return wake
end

@inline function translate_wake!(wakes::AbstractVector{<:AbstractMatrix},
    wake_shedding_locations, Vw, dt; iwake)

    nsurf = length(wakes)
    for isurf = 1:nsurf
        translate_wake!(wakes[isurf], wake_shedding_locations[isurf], Vw[isurf],
            dt; iwake = iwake[isurf])
    end

    return wakes
end

# function shed_wake!(wakes, surfaces)
#     # set circulation of the newly shed panel
#     iΓs = 0 # index for accessing surface circulation
#     iΓw = 0 # index for accessing wake circulation
#     for isurf = 1:nsurf
#         Ns = length(system.surfaces[isurf])
#         Nw = length(system.wakes[isurf])
#         # get view of circulations for this surface/wake
#         vΓs = reshape(view(u_Γ, iΓ+1:iΓ+Ns[isurf]), nc[isurf], ns[isurf])
#         vΓw = reshape(view(u_Γw, iΓ+1:iΓ+Nw[isurf]), nw[isurf], ns[isurf])
#         # replace storage for last panel with new panel
#         copyto!(view(vΓw, nw[isurf], :), view(vΓs, nc[isurf], :))
#         # shift wake panel data downstream
#         rowshift!(vΓw)
#         # move indices to next surface
#         iΓs += Ns[isurf]
#         iΓw += Nw[isurf]
#     end
#     # reset the wake shedding location
#     for isurf = 1:nsurf
#         for j = 1:ns[isurf]+1
#             # extract trailing edge coordinate
#             if j < ns + 1
#                 rte = bottom_left(system.surfaces[isurf][end, j])
#             else
#                 rte = bottom_right(system.surfaces[isurf][end, j-1])
#             end
#             # replace storage for coordinates of last wake panel
#             copyto!(view(u_ζw, 1:3, nw[isurf], :), rte)
#         end
#         # shift wake panel coordinates downstream
#         for i = 1:3
#             rowshift!(view(u_ζw, i, 1:nw[isurf], :))
#         end
#     end
# end



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
# shed_wake!(system::System) = shed_wake!(system.wakes, system)
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

@inline function shed_wake!(wakes::AbstractVector{<:AbstractMatrix},
    wake_shedding_locations, surfaces, Γ; iwake)

    iΓ = 0
    for isurf = 1:length(surfaces)

        N = length(surfaces[isurf])

        vΓ = view(Γ, iΓ+1:iΓ+N)

        shed_wake!(wakes[isurf], wake_shedding_locations[isurf], surfaces[isurf],
            vΓ; iwake = iwake[isurf])

        iΓ += N
    end

    return wakes
end

function shed_wake!(wake::AbstractMatrix, wake_shedding_locations, surface, Γ; iwake)

    nc, ns = size(surface)
    nw = size(wake, 1)
    ls = LinearIndices((nc, ns))

    # replace the last chordwise panels with the newly shed wake panels
    for j = 1:ns

        # trailing edge is top side
        rtl = bottom_left(surface[end, j])
        rtr = bottom_right(surface[end, j])

        # current wake shedding location is bottom side
        rbl = wake_shedding_locations[j]
        rbr = wake_shedding_locations[j+1]

        # use core size from the shedding panel
        core_size = get_core_size(surface[end, j])

        # use circulation strength from the shedding panel
        gamma = Γ[ls[end,j]]

        # replace the oldest wake panel
        wake[end,j] = WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
    end

    # shift wake panels to make newly shed wake panel first
    rowshift!(wake)

    # set new wake shedding locations based on new wake panels
    set_wake_shedding_locations!(wake_shedding_locations, wake)

    return wake, wake_shedding_locations
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
