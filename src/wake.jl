"""
    Wake{TF} <: AbstractPanel

Wake panel element.

**Fields**
 - `rtl`: position of the left side of the top bound vortex
 - `rtr`: position of the right side of the top bound vortex
 - `rbl`: position of the left side of the bottom bound vortex
 - `rbr`: position of the right side of the bottom bound vortex
 - `core_size`: finite core size
 - `gamma`: circulation strength of the panel
"""
struct Wake{TF} <: AbstractPanel{TF}
    rtl::SVector{3, TF}
    rtr::SVector{3, TF}
    rbl::SVector{3, TF}
    rbr::SVector{3, TF}
    core_size::TF
    gamma::TF
end

"""
    Wake(rtl, rtr, rbl, rbr, core_size, gamma)

Construct and return a wake panel element.

**Arguments**
 - `rtl`: position of the left side of the top bound vortex
 - `rtr`: position of the right side of the top bound vortex
 - `rbl`: position of the left side of the bottom bound vortex
 - `rbr`: position of the right side of the bottom bound vortex
 - `core_size`: finite core size
 - `gamma`: circulation strength of the panel
"""
function Wake(rtl, rtr, rbl, rbr, core_size, gamma)

    TF = promote_type(eltype(rtl), eltype(rtr), eltype(rbl), eltype(rbr), typeof(core_size), typeof(gamma))

    return Wake{TF}(rtl, rtr, rbl, rbr, core_size, gamma)
end

@inline Base.eltype(::Type{Wake{TF}}) where TF = TF
@inline Base.eltype(::Wake{TF}) where TF = TF

@inline top_left(panel::Wake) = panel.rtl

@inline top_right(panel::Wake) = panel.rtr

@inline bottom_left(panel::Wake) = panel.rbl

@inline bottom_right(panel::Wake) = panel.rbr

@inline get_core_size(panel::Wake) = panel.core_size

@inline circulation_strength(panel::Wake) = panel.gamma

@inline panel_induced_velocity(rcp, panel::Wake, trailing; kwargs...) = ring_induced_velocity(
    rcp, top_left(panel), top_right(panel), bottom_left(panel), bottom_right(panel),
    trailing; kwargs..., core_size=get_core_size(panel))

@inline function translate(panel::Wake, r)

    rtl = panel.rtl + r
    rtr = panel.rtr + r
    rbl = panel.rbl + r
    rbr = panel.rbr + r
    core_size = panel.core_size
    gamma = panel.gamma

    return Wake(rtl, rtr, rbl, rbr, core_size, gamma)
end

@inline function reflect(panel::Wake)

    rtl = flipy(panel.rtr)
    rtr = flipy(panel.rtl)
    rbl = flipy(panel.rbr)
    rbr = flipy(panel.rbl)
    core_size = panel.core_size
    gamma = panel.gamma

    return Ring(rtl, rtr, rbl, rbr, core_size, gamma)
end

"""
    wake_induced_velocity(rcp, wake; kwargs...)

Returns the induced velocity at `rcp` from the wake panels in `wake`

# Arguments
 - `rcp`: location where induced velocity is calculated
 - `wake`: Matrix of wake panels of shape (nw, ns) where `nw` is the number of
    chordwise wake panels and `ns` is the number of spanwise panels

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image of the panels in `wake`
    should be used when calculating induced velocities.
 - `same_surface`: Flag indicating whether `rcp` corresponds to a panel corner
    on `wake`
 - `same_id`: Flag indicating whether `rcp` is on a surface with the same ID as
    `wake`
 - `trailing_vortices`: Flag that may be used to enable/disable trailing vortices
    shed from `wake`
 - `xhat`: direction in which trailing vortices are shed
 - `nwake`: number of chordwise wake panels to use from `wake`
 - `I`: cartesian index corresponding to the location of `rcp` on `wake`. (1,1)
    corresponds to the top left corner and (nw+1, ns+1) corresponds to the bottom
    right corner.  By default, `rcp` is not assumed to correspond to a panel
    corner on `wake`.
"""
@inline function wake_induced_velocity(rcp, wake; symmetric,
    same_surface, same_id, trailing_vortices, xhat, nwake, I=nothing)

    TF = promote_type(eltype(rcp), eltype(eltype(wake)))

    nw = nwake
    ns = size(wake, 2)
    Ns = nw*ns
    cs = CartesianIndices((nwake, ns))

    finite_core = !same_id

    Vind = @SVector zeros(TF, 3)
    for j = 1:Ns
        J = cs[j]

        leading_edge = J[1] == 1
        trailing_edge = J[1] == nw
        left_side = J[2] == 1
        right_side = J[2] == ns

        # exclude bound vortices that correspond to the corner where `rcp` is
        if isnothing(I)
            include_top = true
            include_bottom = true
            include_left = true
            include_right = true
            include_left_trailing = true
            include_right_trailing = true

            if symmetric
                include_top_mirrored = true
                include_bottom_mirrored = true
                include_left_mirrored = true
                include_right_mirrored = true
                include_left_trailing_mirrored = true
                include_right_trailing_mirrored = true
            end
        else
            vertically_adjacent = same_surface && I[1] == J[1]+1 || I[1] == J[1]
            horizontally_adjacent = same_surface && I[2] == J[2]+1 || I[2] == J[2]

            include_top = !(horizontally_adjacent && I[1] == J[1])
            include_bottom = !(horizontally_adjacent && I[1] == J[1]+1)
            include_left = !(vertically_adjacent && I[2] == J[2])
            include_right = !(vertically_adjacent && I[2] == J[2]+1)
            include_left_trailing = include_left || include_bottom
            include_right_trailing = include_right || include_bottom

            if symmetric
                # ignore the reflected vortex filaments if they also intersect
                # with the control point, otherwise include them
                include_top_mirrored = include_top ||
                    (not_on_symmetry_plane(top_left(wake[J])) &&
                    not_on_symmetry_plane(top_right(wake[J])))
                include_bottom_mirrored = include_bottom ||
                    (not_on_symmetry_plane(bottom_left(wake[J])) &&
                    not_on_symmetry_plane(bottom_right(wake[J])))
                include_left_mirrored = include_left ||
                    (not_on_symmetry_plane(top_left(wake[J])) &&
                    not_on_symmetry_plane(bottom_left(wake[J])))
                include_right_mirrored = include_right ||
                    (not_on_symmetry_plane(top_right(wake[J])) &&
                    not_on_symmetry_plane(bottom_right(wake[J])))
                include_left_trailing_mirrored = include_left_trailing ||
                    not_on_symmetry_plane(bottom_left(wake[J]))
                include_right_trailing_mirrored = include_right_trailing ||
                    not_on_symmetry_plane(bottom_right(wake[J]))
            end
        end

        # trailing vortices instead of bottom bound vortex?
        trailing = trailing_edge && trailing_vortices

        if true#finite_core

            # induced velocity from the panel
            vt, vb, vl, vr, vlt, vrt = panel_induced_velocity(rcp, wake[J], trailing;
                finite_core = finite_core,
                reflect = false,
                xhat = xhat,
                include_top = include_top,
                include_bottom = include_bottom,
                include_left = include_left,
                include_right = include_right,
                include_left_trailing = include_left_trailing,
                include_right_trailing = include_right_trailing)

            # add induced velocity from reflected panel if symmetric
            if symmetric

                vt_s, vb_s, vl_s, vr_s, vlt_s, vrt_s = panel_induced_velocity(rcp,
                    wake[J], trailing;
                    finite_core = finite_core,
                    reflect = true,
                    xhat = xhat,
                    include_top = include_top_mirrored,
                    include_bottom = include_bottom_mirrored,
                    include_left = include_left_mirrored,
                    include_right = include_right_mirrored,
                    include_left_trailing = include_left_trailing_mirrored,
                    include_right_trailing = include_right_trailing_mirrored)

                vt += vt_s
                vb += vb_s
                vl += vl_s
                vr += vr_s
                vlt += vlt_s
                vrt += vrt_s
            end


            # add velocity from this panel
            Vind += (vt + vb + vl + vr + vlt + vrt) * circulation_strength(wake[J])

        else
            # use more efficient formulation when finite core is disabled

            # induced velocity from the panel, excluding already computed sides
            vt, vb, vl, vr, vlt, vrt = panel_induced_velocity(rcp, wake[J], trailing;
                finite_core = finite_core,
                reflect = false,
                xhat = xhat,
                include_top = leading_edge && include_top,
                include_bottom = include_bottom,
                include_left = left_side && include_left,
                include_right = include_right,
                include_left_trailing = left_side && include_left_trailing,
                include_right_trailing = include_right_trailing)

            # add induced velocity from reflected panel if symmetric
            if symmetric

                vt_s, vb_s, vl_s, vr_s, vlt_s, vrt_s = panel_induced_velocity(rcp, wake[J], trailing;
                    finite_core = finite_core,
                    reflect = true,
                    xhat = xhat,
                    include_top = leading_edge && include_top_mirrored,
                    include_bottom = include_bottom_mirrored,
                    include_left = left_side && include_left_mirrored,
                    include_right = include_right_mirrored,
                    include_left_trailing = left_side && include_left_trailing_mirrored,
                    include_right_trailing = include_right_trailing_mirrored)

                vt += vt_s
                vb += vb_s
                vl += vl_s
                vr += vr_s
                vlt += vlt_s
                vrt += vrt_s
            end

            # add velocity from this panel (excluding already computed components)
            Vind += (vt + vb + vl + vr + vlt + vrt) * circulation_strength(wake[J])

            # additional contribution from bottom shared edge
            if !trailing_edge
                Vind -= vb * circulation_strength(wake[j+1])
            end

            # additional contribution from right shared edge
            if !right_side
                Vind -= vr * circulation_strength(wake[j+nw])
            end

            # additional contribution from right shared trailing vortex
            if !right_side
                Vind -= vrt * circulation_strength(wake[j+nw])
            end
        end
    end

    return Vind
end

"""
    get_wake_velocities!(wake_velocities, surface::AbstractMatrix{<:AbstractPanel},
        wake::AbstractMatrix{<:Wake}, ref, fs; symmetric, trailing_vortices, xhat, nwake)

Returns the velocities at the corners of the wake panels in `wake`
"""
@inline function get_wake_velocities!(wake_velocities, surface::AbstractMatrix,
    wake::AbstractMatrix, ref, fs, Γ; symmetric, trailing_vortices, xhat, nwake)

    get_wake_velocities!([wake_velocities], [surface], [wake], ref, fs, Γ;
        symmetric = [symmetric],
        surface_id = 1:1,
        trailing_vortices = [trailing_vortices],
        xhat = xhat,
        nwake = [nwake])

    return wake_velocities
end

"""
    get_wake_velocities!(wake_velocities, surfaces::AbstractVector{<:AbstractMatrix{AbstractPanel}},
        wakes::AbstractVector{<:AbstractMatrix{<:Wake}}, ref, fs; symmetric,
        surface_id, trailing_vortices, xhat, nwake)

Returns the velocities at the corners of the wake panels in `wakes`
"""
@inline function get_wake_velocities!(wake_velocities, surfaces::AbstractVector{<:AbstractMatrix},
    wakes::AbstractVector{<:AbstractMatrix}, ref, fs, Γ; symmetric,
    surface_id, trailing_vortices, xhat, nwake)

    nsurf = length(surfaces)

    for isurf = 1:nsurf

        nc, ns = size(surfaces[isurf])
        nw = nwake[isurf]

        # velocity at the trailing edge
        for is = 1:ns+1

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

                # check if receiving point is on the sending surface
                same_surface = isurf == jsurf

                # check if surfaces have the same ID
                same_id = surface_id[isurf] == surface_id[jsurf]

                # add the induced velocity from the surface, ignoring co-located
                # bound and trailing vortices
                vΓ = view(Γ, jΓ+1:jΓ+Ns)
                wake_velocities[isurf][1,is] += surface_induced_velocity(rc, surfaces[jsurf], vΓ;
                    symmetric = symmetric[jsurf],
                    same_surface = same_surface,
                    same_id = same_id,
                    trailing_vortices = false,
                    xhat = xhat, # not used
                    Ic = CartesianIndex(nc+1, is))

                # add induced velocity from wake, ignoring co-located bound and
                # trailing vortices
                wake_velocities[isurf][1,is] += wake_induced_velocity(rc, wakes[jsurf];
                    symmetric = symmetric[jsurf],
                    same_surface = same_surface,
                    same_id = same_id,
                    trailing_vortices = trailing_vortices[jsurf],
                    xhat = xhat,
                    nwake = nwake[jsurf],
                    I = CartesianIndex(1, is))

                jΓ += Ns # increment Γ index
            end
        end

        # velocity at all other wake corners
        cr = CartesianIndices((2:nw+1, 1:ns+1))

        for I in cr

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

                # also check if it has the same ID
                same_id = surface_id[isurf] == surface_id[jsurf]

                # add induced velocity from surface
                vΓ = view(Γ, jΓ+1:jΓ+Ns)
                wake_velocities[isurf][I] += surface_induced_velocity(rc, surfaces[jsurf], vΓ;
                    symmetric = symmetric[jsurf],
                    same_surface = false,
                    same_id = same_id,
                    trailing_vortices = false,
                    xhat = xhat)

                # add induced velocity from wake
                wake_velocities[isurf][I] += wake_induced_velocity(rc, wakes[jsurf];
                    symmetric = symmetric[jsurf],
                    same_surface = isurf == jsurf,
                    same_id = same_id,
                    trailing_vortices = trailing_vortices[jsurf],
                    xhat = xhat,
                    nwake = nwake[jsurf],
                    I = I)

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
        vΓ = view(Γ, nc*(j-1)+1:(nc)*(j-1)+nc)
        gamma = trailing_edge_circulation(surface, vΓ)

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

    @inbounds for j = 1:nj
        tmp = A[ni,j]
        @inbounds for i = ni:-1:2
            A[i,j] = A[i-1,j]
        end
        A[1,j] = tmp
    end

    return A
end
