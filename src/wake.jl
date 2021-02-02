"""
    update_wake_shedding_locations!(wake, wake_shedding_locations, surface,
        freestream, reference, dt; kwargs...)

Update the wake shedding locations.  Also update the first chordwise wake panel
to account for the new wake shedding location

# Arguments
 - `wake`: Matrix of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns)
    where `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels, defaults to no wake panels
 - `wake_shedding_locations`: Vector of length `ns` which stores the coordinates
    where wake panels are shed from the trailing edge of `surface`.
 - `surface`: Matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))
 - `dt`: Time step (seconds)

# Keyword Arguments
 - `nwake = size(wake, 1)`: Number of chordwise wake panels.
 - `eta = 0.25`: Time step fraction used to define separation between trailing
    edge and wake shedding location.  Typical values range from 0.2-0.3.
"""
@inline function update_wake_shedding_locations!(wake::AbstractMatrix,
    wake_shedding_locations, surface::AbstractMatrix, fs, ref, dt;
    nwake = size(wake, 1), eta = 0.25)

    # number of spanwise panels
    ns = size(wake, 2)

    # update wake shedding location
    for j = 1:ns+1

        # extract trailing edge coordinate
        if j < ns + 1
            rte = bottom_left(surface[end, j])
        else
            rte = bottom_right(surface[end, j-1])
        end

        # extract corresponding velocity
        Vte = external_velocity(fs, rte, ref.r)

        # update wake shedding location coordinates
        wake_shedding_locations[j] = rte + eta*Vte*dt

    end

    if nwake > 0
        # loop through first row of wake panels
        for j = 1:ns
            # update wake panel with wake shedding location coordinates
            rtl = wake_shedding_locations[j]
            rtr = wake_shedding_locations[j+1]

            # preserve other wake panel coordinates
            rbl = bottom_left(wake[1,j])
            rbr = bottom_right(wake[1,j])

            # preserve core size
            core_size = get_core_size(wake[1,j])

            # preserve circulation strength
            gamma = circulation_strength(wake[1,j])

            # replace the old wake panel
            wake[1,j] = WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
        end
    end

    return wake, wake_shedding_locations
end

"""
    update_wake_shedding_locations!(wakes, wake_shedding_locations, surfaces,
        freestream, reference, dt; kwargs...)

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

# Keyword Arguments
 - `nwake`: Number of chordwise wake panels to use from each wake in `wakes`
 - `eta = 0.25`: Time step fraction used to define separation between trailing
    edge and wake shedding location.  Typical values range from 0.2-0.3.
"""
@inline function update_wake_shedding_locations!(wakes::AbstractVector{<:AbstractMatrix},
    wake_shedding_locations, surfaces::AbstractVector{<:AbstractMatrix}, fs, ref, dt;
    nwake = size.(wakes, 1),
    eta = 0.25)

    # get number of surfaces
    nsurf = length(surfaces)

    # loop through all surfaces
    for i = 1:nsurf

        # update the wake shedding location for this surface
        update_wake_shedding_locations!(wakes[i], wake_shedding_locations[i],
            surfaces[i], fs, ref, dt; nwake = nwake[i], eta)

    end

    return wakes, wake_shedding_locations
end

"""
    wake_normal_velocity(surface, wake; kwargs...)

Compute the normal component of the velocity induced on a surface by its own
wake panels

This forms part of the right hand side of the circulation linear system solve.

# Arguments
 - `surface`: Matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `wake`: Matrix of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns)
    where `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image (across the
    X-Z plane) of the panels in `surface` should be used when calculating induced velocities
 - `nwake`: Number of chordwise wake panels to use from `wake`, defaults to all
    wake panels
 - `trailing_vortices`: Flag to enable/disable trailing vortices, defaults to `true`
 - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]
"""
function wake_normal_velocity(surface::AbstractMatrix, wake::AbstractMatrix; kwargs...)

    TF = promote_type(eltype(eltype(surface)), eltype(eltype(wake)))
    N = length(surface)
    b = zeros(TF, N)

    return -subtract_wake_normal_velocity!(b, surface, wake; kwargs...)
end

"""
    wake_normal_velocity(surfaces, wakes; kwargs...)

Compute the normal component of the velocity induced on multiple surfaces by
their wake panels

This forms part of the right hand side of the circulation linear system solve.

# Arguments
 - `surfaces`: Vector of surfaces, represented by matrices of surface panels
    (see [`SurfacePanel`](@ref) of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `wakes`: Vector of wakes corresponding to each surface, represented by matrices
    of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where `nw` is the
    number of chordwise wake panels and `ns` is the number of spanwise panels.

# Keyword Arguments
 - `symmetric`: Flag indicating whether a mirror image of the panels in `wake`
    should be used when calculating induced velocities
 - `nwake`: Number of chordwise wake panels to use from each wake in `wakes`,
    defaults to all provided wake panels
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.  Defaults to `true` for each surface.
 - `trailing_vortices`: Flags to enable/disable trailing vortices, defaults to
    `true` for each surface
 - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]
"""
function wake_normal_velocity(surfaces::AbstractVector{<:AbstractMatrix},
    wakes::AbstractVector{<:AbstractMatrix}; kwargs...)

    TF = promote_type(eltype(eltype(eltype(surfaces))), eltype(eltype(eltype(wakes))))
    N = sum(length.(surfaces))
    b = zeros(TF, N)

    return -subtract_wake_normal_velocity!(b, surfaces, wakes; kwargs...)
end

"""
    subtract_wake_normal_velocity!(b, surface, wake; kwargs...)

Pre-allocated version of [`wake_normal_velocity`](@ref) which subtracts the
normal induced velocity on the surface panels in `surface` created by the wake
panels in `wake` from the existing vector `b`
"""
@inline function subtract_wake_normal_velocity!(b, surface::AbstractMatrix, wake::AbstractMatrix;
    wake_finite_core = false,
    nwake = size(wake, 1),
    kwargs...)

    finite_core = wake_finite_core

    # loop over receiving panels
    for i = 1:length(surface)

        # control point location
        rcp = controlpoint(surface[i])

        # normal vector body axis
        nhat = normal(surface[i])

        # get induced velocity at rcp from the wake
        Vind = induced_velocity(rcp, wake;
            finite_core = wake_finite_core,
            nc = nwake,
            kwargs...)

        # subtract normal velocity
        b[i] -= dot(Vind, nhat)
    end

    return b
end

"""
    subtract_wake_normal_velocity!(b, surfaces, wakes; kwargs...)

Pre-allocated version of [`wake_normal_velocity`](@ref) which subtracts the normal
induced velocity on the surface panels in `surfaces` created by the wake panels
in `wakes` from the existing vector `b`
"""
function subtract_wake_normal_velocity!(b, surfaces::AbstractVector{<:AbstractMatrix},
    wakes::AbstractVector{<:AbstractMatrix};
    symmetric = fill(false, length(surfaces)),
    nwake = size.(wakes, 1),
    surface_id = 1:length(surfaces),
    wake_finite_core = fill(true, length(surfaces)),
    trailing_vortices = fill(true, length(surfaces)),
    xhat = SVector(1, 0, 0))

    nsurf = length(surfaces)

    # index for keeping track of where we are in the b vector
    ib = 0

    # loop through receiving surfaces
    for i = 1:nsurf

        # number of panels on this surface
        n = length(surfaces[i])

        # create view into RHS vector
        vb = view(b, ib+1:ib+n)

        # fill in RHS vector
        for j = 1:nsurf

            subtract_wake_normal_velocity!(vb, surfaces[i], wakes[j];
                finite_core = wake_finite_core[j] || (surface_id[i] != surface_id[j]),
                symmetric = symmetric[j],
                nwake = nwake[j],
                trailing_vortices = trailing_vortices[j],
                xhat = xhat)
        end

        # increment position in AIC matrix
        ib += n
    end

    return b
end

"""
    get_wake_velocities!(wake_velocities, surface, wake, wake_shedding_locations,
    reference, freestream; kwargs...)

Return the velocities at the corners of the wake panels in `wake`

# Arguments
 - `wake_velocities`: Velocities at the corners of the wake panels in `wake`
 - `surface`: Matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels'
 - `wake`: Matrix of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns)
    where `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels, defaults to no wake panels
 - `wake_shedding_locations`: Vector of length `ns` which stores the coordinates
    where wake panels are shed from the trailing edge of `surface`.
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image of the panels
    in `surface` should be used when calculating induced velocities
 - `repeated_points`: Dictionary of the form `Dict((isurf, i) => [(jsurf1, j1),
    (jsurf2, j2)...]` which defines repeated trailing edge points.  Trailing edge
    point `i` on surface `isurf` is repeated on surface `jsurf1` at point `j1`,
    `jsurf2` at point `j2`, and so forth. See [`repeated_trailing_edge_points`](@ref)
 - `nwake`: Maximum number of wake panels in the chordwise direction.  Defaults
    to `length(dx)`.
 - `wake_finite_core`: Flag indicating whether the finite core model should be
    enabled when calculating the wake's influence on itself and its corresponding
    surface. Defaults to `true`
 - `trailing_vortices`: Flag to enable/disable trailing vortices, defaults to `true`
 - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]
"""
@inline function get_wake_velocities!(wake_velocities, surface::AbstractMatrix,
    wake_shedding_locations::AbstractVector,
    wake::AbstractMatrix, ref, fs, Γ;
    symmetric,
    repeated_points = Dict{NTuple{2, Int}, Vector{NTuple{2, Int}}}(),
    nwake = size(wake, 1),
    wake_finite_core = true,
    trailing_vortices = false,
    xhat = SVector(1, 0, 0))

    get_wake_velocities!([wake_velocities], [surface], [wake_shedding_locations],
        [wake], ref, fs, Γ;
        symmetric = [symmetric],
        surface_id = 1:1,
        wake_finite_core = wake_finite_core,
        trailing_vortices = [trailing_vortices],
        xhat = xhat,
        nwake = [nwake],
        repeated_points = repeated_points)

    return wake_velocities
end

"""
    get_wake_velocities!(wake_velocities, surfaces, wakes, wake_shedding_locations,
        ref, fs; symmetric, surface_id, wake_finite_core,
        trailing_vortices, xhat, nwake, repeated_points)

# Arguments
 - `wake_velocities`: Velocities at the corners of the wake panels in `wakes`
 - `surfaces`: Vector of surfaces, represented by matrices of surface panels
    (see [`SurfacePanel`](@ref) of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `wakes`: Vector of wakes corresponding to each surface, represented by matrices
    of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where `nw` is the
    number of chordwise wake panels and `ns` is the number of spanwise panels.
 - `wake_shedding_locations`: Shedding location coordinates for each surface for
    each trailing edge vertex.
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))

# Keyword Arguments
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
 - `trailing_vortices`: Flags to enable/disable trailing vortices, defaults to
    `true` for each surface
 - `xhat`: Direction in which to shed trailing vortices, defaults to [1, 0, 0]
"""
@inline function get_wake_velocities!(wake_velocities, surfaces::AbstractVector{<:AbstractMatrix},
    wake_shedding_locations::AbstractVector{<:AbstractVector},
    wakes::AbstractVector{<:AbstractMatrix}, ref, fs, Γ;
    symmetric,
    repeated_points = Dict{NTuple{2, Int}, Vector{NTuple{2, Int}}}(),
    nwake = size.(wakes, 1),
    surface_id = 1:length(surfaces),
    wake_finite_core = fill(true, length(surfaces)),
    trailing_vortices = fill(false, length(surfaces)),
    xhat = SVector(1, 0, 0))

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

            # velocity of trailing edge
            wake_velocities[isurf][1,is] = external_velocity(fs, rc, ref.r)
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

            # start with external velocity
            wake_velocities[isurf][I] = external_velocity(fs, rc, ref.r)

            # add induced velocity from each surface and wake
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

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

                # add induced velocity from the surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)
                wake_velocities[isurf][I] += induced_velocity(rc, surfaces[jsurf], vΓ;
                    symmetric = symmetric[jsurf],
                    finite_core = surface_id[isurf] != surface_id[jsurf],
                    wake_shedding_locations = wake_shedding_locations[jsurf],
                    trailing_vortices = false)

                # add induced velocity from the wake
                if same_surface
                    # vertex location on wake
                    J = CartesianIndex(I[1], js)

                    # induced velocity from wake on its own vertex
                    wake_velocities[isurf][I] += induced_velocity(J, wakes[jsurf];
                        symmetric = symmetric[jsurf],
                        finite_core = wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
                        nc = nwake[jsurf],
                        trailing_vortices = trailing_vortices[jsurf],
                        xhat = xhat)
                else
                    # induced velocity from wake on another wake's vertex
                    wake_velocities[isurf][I] += induced_velocity(rc, wakes[jsurf];
                        symmetric = symmetric[jsurf],
                        finite_core = wake_finite_core[jsurf] || surface_id[isurf] != surface_id[jsurf],
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
    shed_wake!(wake, wake_shedding_locations, wake_velocities, dt, surface, Γ; nwake)

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

# Keyword Arguments
 - `nwake`: Number of chordwise wake panels to use from `wake`, defaults to all
    provided wake panels
"""
@inline function shed_wake!(wake::AbstractMatrix, wake_shedding_locations,
    wake_velocities, dt, surface, Γ; nwake)

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
    translate_wake!(wake, wake_velocities, dt; nwake = min(nwake, nw-1))

    # shift wake panels to make newly shed wake panel first
    rowshift!(wake)

    return wake
end

"""
    shed_wake!(wakes, wake_shedding_locations, wake_velocities, dt, surfaces, Γ; nwake)

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

# Keyword Arguments
 - `nwake`: Number of chordwise wake panels to use from each wake in `wakes`,
    defaults to all provided wake panels
"""
@inline function shed_wake!(wakes::AbstractVector{<:AbstractMatrix}, wake_shedding_locations,
    wake_velocities, dt, surfaces, Γ; nwake)

    iΓ = 0
    for i = 1:length(surfaces)

        N = length(surfaces[i])

        vΓ = view(Γ, iΓ+1:iΓ+N)

        shed_wake!(wakes[i], wake_shedding_locations[i], wake_velocities[i], dt,
            surfaces[i], vΓ; nwake=nwake[i])

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
