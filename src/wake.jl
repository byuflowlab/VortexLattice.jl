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

@inline induced_velocity(rcp, panel::Wake, trailing; kwargs...) = ring_induced_velocity(
    rcp, top_left(panel), top_right(panel), bottom_left(panel), bottom_right(panel),
    trailing; kwargs..., core_size=get_core_size(panel))

# @inline panel_induced_velocity(receiving, Ir, sending::Wake, Is, same_surface, trailing;
#     kwargs...) = induced_velocity(top_center(receiving), sending, trailing; kwargs...,
#     include_top = !(same_surface && Ir == Is),
#     include_bottom = !(same_surface && (Ir[1] == Is[1]+1 && Ir[2] == Is[2])))

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
    wake_induced_velocity(rcp, wake, symmetric, same_id, trailing_vortices; kwargs...)

Returns the induced velocity at `rcp` from the wake panels in `wake`

# Arguments
 - `rcp`: location where induced velocity is calculated
 - `wake`: Matrix of wake panels of shape (nw, ns) where `nw` is the number of
    chordwise wake panels and `ns` is the number of spanwise panels
 - `symmetric`: Flag indicating whether a mirror image of the panels in `wake`
    should be used when calculating induced velocities.
 - `trailing_vortices`: Flag that may be used to enable/disable trailing vortices
    shed from `wake`
 - `same_surface`: Flag indicating whether `rcp` corresponds to a panel corner
    on `wake`
 - `same_id`: Flag indicating whether `rcp` is on a surface with the same ID as
    `wake`

# Keyword Arguments
 - `xhat`: direction in which trailing vortices are shed, defaults to [1, 0, 0]
 - `nwake`: number of chordwise wake panels to use from `wake`, defaults to all
    wake panels
 - `I`: cartesian index corresponding to the location of `rcp` on `wake`. (1,1)
    corresponds to the top left corner and (nw+1, ns+1) corresponds to the bottom
    right corner.  By default, `rcp` is not assumed to correspond to a panel
    corner on `wake`.
"""
function wake_induced_velocity(rcp, wake, same_surface, same_id, trailing_vortices;
    xhat=SVector(1, 0, 0), nwake=size(wake, 1), I=CartesianIndex(-1, -1))

    TF = promote_type(eltype(rcp), eltype(eltype(wake)))

    nc, ns = size(wake)
    cs = CartesianIndices((nwake, ns))

    finite_core = !same_id

    Vind = @SVector zeros(TF, 3)
    for j = 1:Ns
        J = cs[j]

        leading_edge = J[1] == 1
        trailing_edge = J[1] == nwake
        left_side = J[2] == 1
        right_side = J[2] == ns

        trailing = trailing_edge && trailing_vortices

        # exclude bound vortices that correspond to the corner where `rcp` is
        vertically_adjacent = I[1] == J[1]-1 || I[1] == J[1]
        horizontally_adjacent = I[2] == J[2]-1 || I[2] == J[2]

        include_top = !(same_surface && I[1] == J[1] && horizontally_adjacent)
        include_bottom = !(same_surface && I[1] == J[1]-1 && horizontally_adjacent)
        include_left = !(same_surface || I[2] == J[2] && vertically_adjacent)
        include_right = !(same_surface || I[2] != J[2]-1 && vertically_adjacent)
        include_left_trailing = include_left || include_bottom
        include_right_trailing = include_right || include_bottom

        if finite_core

            # induced velocity from the panel
            vt, vb, vl, vr, vlt, vrt = induced_velocity(rcp, wake[j], trailing;
                finite_core = finite_core, reflect = false, xhat = xhat,
                include_top=include_top, include_bottom=include_bottom,
                include_left=include_left, include_right=include_right,
                include_left_trailing=include_left_trailing,
                include_right_trailing=include_right_trailing)

            # add induced velocity from reflected panel if symmetric
            if symmetric
                vt_s, vb_s, vl_s, vr_s, vlt_s, vrt_s = induced_velocity(rcp,
                    wake[j], trailing; finite_core = finite_core, reflect = true,
                    xhat = xhat,
                    include_top=include_top, include_bottom=include_bottom,
                    include_left=include_left, include_right=include_right,
                    include_left_trailing=include_left_trailing,
                    include_right_trailing=include_right_trailing)

                vt += vt_s
                vb += vb_s
                vl += vl_s
                vr += vr_s
                vlt += vlt_s
                vrt += vrt_s
            end

            # add velocity from this panel
            Vind += (vt + vb + vl + vr + vlt + vrt) * circulation_strength(wake[j])

        else
            # use more efficient formulation when finite core is disabled

            # only comopute top bound vortex contribution for the leading edge,
            # reuse the bottom bound vortex contribution everywhere else
            include_top = include_top && leading_edge

            # only comopute left bound vortex contribution for the left side,
            # reuse the right bound vortex contribution everywhere else
            include_left = include_left && left_side

            # only comopute left trailing vortex contribution for the left side,
            # reuse the right bound vortex contribution everywhere else
            include_left_trailing = include_left_trailing && left_side

            # induced velocity from the panel, excluding already computed sides
            vt, vb, vl, vr, vlt, vrt = induced_velocity(rcp, wake[j], trailing;
                finite_core = finite_core,  reflect = false, xhat = xhat,
                include_top=include_top, include_bottom=include_bottom,
                include_left=include_left, include_right=include_right,
                include_left_trailing=include_left_trailing,
                include_right_trailing=include_right_trailing)

            # add induced velocity from reflected panel if symmetric
            if symmetric
                vt_s, vb_s, vl_s, vr_s, vlt_s, vrt_s = induced_velocity(rcp, wake[j], trailing;
                    finite_core = finite_core,  reflect = true, xhat = xhat,
                    include_top=include_top, include_bottom=include_bottom,
                    include_left=include_left, include_right=include_right,
                    include_left_trailing=include_left_trailing,
                    include_right_trailing=include_right_trailing)

                vt += vt_s
                vb += vb_s
                vl += vl_s
                vr += vr_s
                vlt += vlt_s
                vrt += vrt_s
            end

            # add velocity from this panel (excluding already computed components)
            Vind += (vt + vb + vl + vr + vlt + vrt) * circulation_strength(wake[j])

            # additional contribution from bottom shared edge
            if !trailing_edge
                Vind -= vb * circulation_strength(wake[j+1])
            end

            # additional contribution from right shared edge
            if !right_side
                Vind -= vr * circulation_strength(wake[j+nchordwise])
            end

            # additional contribution from right shared trailing vortex
            if !right_side
                Vind -= vrt * circulation_strength(wake[j+nchordwise])
            end
        end
    end

    return Vind
end

"""
    wake_velocities!(V, surface::AbstractMatrix{<:AbstractPanel},
        wake::AbstractMatrix{<:Wake}, surface_id, trailing_vortices,
        symmetric, ref, fs; nwake=size(wakes, 1), xhat=SVector(1, 0, 0))

Returns the velocities at the corners of the wake panels in `wake`
"""
function wake_velocities!(V, surface::AbstractMatrix{<:AbstractPanel},
    wake::AbstractMatrix{<:Wake}, surface_id, trailing_vortices,
    symmetric, ref, fs; nwake=size(wakes, 1), xhat=SVector(1, 0, 0))

    nsurf = length(surfaces)


    for isurf = 1:nsurf
        receiving = wake[i]

        nw = nwake[isurf]
        ns = size(receiving, 2)

        Nr = nw*ns
        cr = CartesianIndices((nw, ns))

        # loop through all wake panels
        for i = 1:Nr
            I = cr[i]

            V[I] = 0.0

            # add contribution from surfaces
            for jsurf = 1:nsurf
                sending = surfaces[jsurf]

                # also check if it has the same ID
                same_id = surface_id[isurf] == surface_id[jsurf]

                # add induced velocity from surface
                same_surface = false
                surface_trailing_vortices = false
                V[isurf][I] += surface_induced_velocity(rcp, wake, same_surface,
                    same_id, surface_trailing_vortices)

                # add induced velocity from wake
                same_surface = isurf == jsurf
                wake_trailing_vortices = trailing_vortices
                V[isurf][I] += wake_induced_velocity(rcp, wake, same_surface,
                    same_id, trailing_vortices; xhat=xhat, nwake=nwake[jsurf], I=I)

            end
        end
    end

    return V
end

function wake_velocities!(V, surfaces::AbstractVector{<:AbstractMatrix{AbstractPanel}},
    wakes::AbstractVector{<:AbstractMatrix{<:Wake}}, surface_id, trailing_vortices,
    symmetric, ref, fs; nwake=size.(wakes, 1), xhat=SVector(1, 0, 0))

    nsurf = length(surfaces)


    for isurf = 1:nsurf
        receiving = wake[i]

        nw = nwake[isurf]
        ns = size(receiving, 2)

        Nr = nw*ns
        cr = CartesianIndices((nw+1, ns+1))

        # loop through all wake panels
        for i = 1:Nr
            I = cr[i]

            # extract relevant corner
            if I[1] <= nw && I[2] <= ns
                rc = top_left(receiving[I[1], I[2]])
            elseif I[1] == nw + 1 && I[2] <= ns
                rc = bottom_left(receiving[I[1]-1, I[2]])
            elseif I[1] <= nw && I[2] == ns + 1
                rc = top_right(receiving[I[1], I[2]-1])
            else # I[1] == nw + 1 && I[2] == ns + 1
                rc = bottom_right(receiving[I[1]-1, I[2]-1])
            end

            # external velocity
            V[isurf][I] = external_velocity(fs, rc, ref.r)

            # induced velocity from each surface and wake
            for jsurf = 1:nsurf
                sending = surfaces[jsurf]

                # also check if it has the same ID
                same_id = surface_id[isurf] == surface_id[jsurf]

                # add induced velocity from surface
                same_surface = false
                surface_trailing_vortices = false
                V[isurf][I] += surface_induced_velocity(rc, wake, same_surface,
                    same_id, surface_trailing_vortices)

                # add induced velocity from wake
                same_surface = isurf == jsurf
                wake_trailing_vortices = trailing_vortices
                V[isurf][I] += wake_induced_velocity(rc, wake, same_surface,
                    same_id, trailing_vortices; xhat=xhat, nwake=nwake[jsurf], I=I)

            end
        end
    end

    return V
end
