# --- Vortex Filament Induced Velocity Functions --- #

"""
    bound_induced_velocity(r1, r2, finite_core, core_size)

Compute the induced velocity (per unit circulation) for a bound vortex, at a
control point located at `r1` relative to the start of the bound vortex and `r2`
relative to the end of the bound vortex
"""
@inline function bound_induced_velocity(r1, r2, finite_core, core_size)

    nr1 = norm(r1)
    nr2 = norm(r2)

    if finite_core
        rdot = dot(r1, r2)
        r1s, r2s, εs = nr1^2, nr2^2, core_size^2
        f1 = cross(r1, r2)/(r1s*r2s - rdot^2 + εs*(r1s + r2s - 2*nr1*nr2))
        f2 = (r1s - rdot)/sqrt(r1s + εs) + (r2s - rdot)/sqrt(r2s + εs)
    else
        f1 = cross(r1, r2)/(nr1*nr2 + dot(r1, r2))
        f2 = (1/nr1 + 1/nr2)
    end

    Vhat = (f1*f2)/(4*pi)

    return Vhat
end

"""
    trailing_induced_velocity(r1, r2, xhat, finite_core, core_size)

Compute the induced velocity (per unit circulation) for a vortex trailing in
the `xhat` direction, at a control point located at `r` relative to the start of the
trailing vortex.
"""
@inline function trailing_induced_velocity(r, xhat, finite_core, core_size)

    nr = norm(r)

    rdot = dot(r, xhat)

    if finite_core
        εs = core_size^2
        tmp = εs/(nr + rdot)
        f = cross(r, xhat)/(nr*(nr - rdot + tmp))
    else
        f = cross(r, xhat)/(nr*(nr - rdot))
    end

    Vhat = -f/(4*pi)

    return Vhat
end

"""
    ring_induced_velocity(rcp, r11, r12, r21, r22; finite_core = false,
        core_size = 0.0, symmetric = false, xhat = [1,0,0], top = true, bottom = true,
        left = true, right = true, left_trailing = false, right_trailing = false)

Compute the induced velocity (per unit circulation) for a vortex ring defined by
the corners `r11`, `r12`, `r21`, and `r22` at a control point located at `rcp`

Also returns the induced velocity resulting from shared edges from adjacent panels
on the top, bottom, left, and right sides of panel.
"""
@inline function ring_induced_velocity(rcp, r11, r12, r21, r22;
    finite_core = false, core_size = 0.0, symmetric = false, xhat = SVector(1, 0, 0),
    top = true, bottom = true, left = true, right = true, left_trailing = false,
    right_trailing = false, reflected_top=true, reflected_bottom = true,
    reflected_left = true, reflected_right = true, reflected_left_trailing = false,
    reflected_right_trailing = false)

    # move origin to control point
    r1 = rcp - r11
    r2 = rcp - r12
    r3 = rcp - r22
    r4 = rcp - r21

    Vhat = zero(rcp) # current panel (from specified edges)
    Vhat_t = zero(rcp) # panel above this panel (from shared edges)
    Vhat_l = zero(rcp) # panel to the left of this panel (from shared edges)
    Vhat_r = zero(rcp) # panel to the right of this panel (from shared edges)
    Vhat_b = zero(rcp) # panel below this panel (from shared edges)

    # top edge
    if top
        vt = bound_induced_velocity(r1, r2, finite_core, core_size)
        # shared edge contribution from current panel
        Vhat += vt
        # shared edge contribution from panel above this panel
        Vhat_t -= vt
    end

    # right edge
    if right
        vr = bound_induced_velocity(r2, r3, finite_core, core_size)
        # shared edge contribution from current panel
        Vhat += vr
        # shared edge contribution from panel to the right of this panel
        Vhat_r -= vr
    end

    # right trailing vortex
    if right_trailing
        vrt = trailing_induced_velocity(r3, xhat, finite_core, core_size)
        # shared edge contribution from current panel
        Vhat += vrt
        # shared edge contribution from panel to the right of this panel
        Vhat_r -= vrt
    end

    # bottom edge
    if bottom
        vb = bound_induced_velocity(r3, r4, finite_core, core_size)
        # shared edge contribution from current panel
        Vhat += vb
        # shared edge contribution from panel below this panel
        Vhat_b -= vb
    end

    # left trailing vortex
    if left_trailing
        vlt = trailing_induced_velocity(r4, xhat, finite_core, core_size)
        # shared edge contribution from current panel
        Vhat -= vlt
        # shared edge contribution from panel to the left of this panel
        Vhat_l += vlt
    end

    # left edge
    if left
        vl = bound_induced_velocity(r4, r1, finite_core, core_size)
        # shared edge contribution from current panel
        Vhat += vl
        # shared edge contribution from panel to the left of this panel
        Vhat_l -= vl
    end

    if symmetric

        # move origin to control point
        r1 = rcp - flipy(r11)
        r2 = rcp - flipy(r12)
        r3 = rcp - flipy(r22)
        r4 = rcp - flipy(r21)

        # top edge
        if reflected_top
            vt = bound_induced_velocity(r2, r1, finite_core, core_size)
            # shared edge contribution from current panel
            Vhat += vt
            # shared edge contribution from panel above this panel
            Vhat_t -= vt
        end

        # right edge
        if reflected_right
            vr = bound_induced_velocity(r3, r2, finite_core, core_size)
            # shared edge contribution from current panel
            Vhat += vr
            # shared edge contribution from panel to the right of this panel
            Vhat_r -= vr
        end

        # right trailing vortex
        if reflected_right_trailing
            vrt = trailing_induced_velocity(r3, xhat, finite_core, core_size)
            # shared edge contribution from current panel
            Vhat -= vrt
            # shared edge contribution from panel to the right of this panel
            Vhat_r += vrt
        end

        # bottom edge
        if reflected_bottom
            vb = bound_induced_velocity(r4, r3, finite_core, core_size)
            # shared edge contribution from current panel
            Vhat += vb
            # shared edge contribution from panel below this panel
            Vhat_b -= vb
        end

        # left trailing vortex
        if reflected_left_trailing
            vlt = trailing_induced_velocity(r4, xhat, finite_core, core_size)
            # shared edge contribution from current panel
            Vhat += vlt
            # shared edge contribution from panel to the left of this panel
            Vhat_l -= vlt
        end

        # left edge
        if reflected_left
            vl = bound_induced_velocity(r1, r4, finite_core, core_size)
            # shared edge contribution from current panel
            Vhat += vl
            # shared edge contribution from panel to the left of this panel
            Vhat_l -= vl
        end
    end

    return Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r
end

"""
    ring_induced_velocity(rcp, panel; kwargs...)

Compute the induced velocity (per unit circulation) for the vortex ring panel
`panel` at a control point located at `rcp`

Also returns the induced velocity resulting from shared edges from adjacent panels
on the top, bottom, left, and right sides of panel.
"""
@inline function ring_induced_velocity(rcp, panel; kwargs...)

    r11 = top_left(panel)
    r12 = top_right(panel)
    r21 = bottom_left(panel)
    r22 = bottom_right(panel)
    core_size = get_core_size(panel)

    return ring_induced_velocity(rcp, r11, r12, r21, r22; core_size=core_size, kwargs...)
end

"""
    influence_coefficients!(AIC, surface; kwargs...)

Construct the aerodynamic influence coefficient matrix for a single surface.

# Arguments:
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels

# Keyword Arguments
 - `symmetric`: Flag indicating whether a mirror image of the panels in `surface`
    should be used when calculating induced velocities.
 - `trailing_vortices`: Flag to enable/disable trailing vortices.
 - `xhat`: Direction in which trailing vortices are shed if `trailing_vortices = true`.
    Defaults to [1, 0, 0]
"""
@inline influence_coefficients!(AIC, surface::AbstractMatrix; surface_id=1, kwargs...) =
    influence_coefficients!(AIC, surface, surface; kwargs..., finite_core = surface_id < 0)

"""
    influence_coefficients!(AIC, surfaces; kwargs...)

Construct the aerodynamic influence coefficient matrix for multiple surfaces.

# Arguments:
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels

# Keyword Arguments:
 - `surface_id`: ID for each surface.  May be used to deactivate the finite core
    model by setting all surface ID's to the same value.  Defaults to a unique ID
    for each surface
 - `symmetric`: Flags indicating whether a mirror image (across the X-Z plane) should
    be used when calculating induced velocities. Defaults to `false` for each surface.
 - `trailing_vortices`: Flags to indicate whether trailing vortices are used for
    each surface. Defaults to `true` for each surface.
 - `xhat`: Direction in which trailing vortices are shed if `trailing_vortices = true`.
    Defaults to [1, 0, 0]
"""
@inline function influence_coefficients!(AIC, surfaces::AbstractVector{<:AbstractMatrix};
    surface_id = 1:length(surfaces), symmetric = fill(false, length(surfaces)),
    trailing_vortices = fill(true, length(surfaces)), xhat = SVector(1, 0, 0))

    nsurf = length(surfaces)

    # indices for keeping track of where we are in the AIC matrix
    iAIC = 0
    jAIC = 0

    # loop through receving surfaces
    for i = 1:nsurf
        receiving = surfaces[i]

        # extract number of panels on this receiving surface
        nr = length(receiving) # number of panels on this receiving surface

        # loop through sending surfaces
        jAIC = 0
        for j = 1:nsurf
            sending = surfaces[j]

            # extract number of panels on this sending surface
            ns = length(sending)

            # extract portion of AIC matrix for the two surfaces
            vAIC = view(AIC, iAIC+1:iAIC+nr, jAIC+1:jAIC+ns)

            # check if it's the same surface
            finite_core = surface_id[i] != surface_id[j]

            # populate entries in the AIC matrix
            influence_coefficients!(vAIC, receiving, sending;
                finite_core = finite_core,
                symmetric = symmetric[j],
                trailing_vortices = trailing_vortices[j],
                xhat = xhat)

            # increment position in AIC matrix
            jAIC += ns
        end
        # increment position in AIC matrix
        iAIC += nr
    end

    return AIC
end

"""
    influence_coefficients!(AIC, receiving, sending; kwargs...)

Compute the AIC coefficients corresponding to the influence of the panels in
`sending` on the panels in `receiving`.

# Keyword Arguments
 - `finite_core`: Flag indicating whether the finite core model is enabled. Defaults
    to `true`
 - `symmetric`: Flag indicating whether sending panels should be mirrored across the X-Z plane
 - `trailing_vortices`: Indicates whether trailing vortices are used. Defaults to `true`.
 - `xhat`: Direction in which trailing vortices are shed if `trailing_vortices = true`.
    Defaults to [1, 0, 0]
"""
@inline function influence_coefficients!(AIC, receiving, sending;
    finite_core = true,
    symmetric = false,
    trailing_vortices = true,
    xhat = SVector(1, 0, 0))

    # zero out AIC matrix
    AIC .= 0.0

    # get sending surface dimensions
    nc, ns = size(sending)
    ls = LinearIndices(sending)

    # determine control flow parameters based on inputs
    reuse_edges = !finite_core

    # loop over receiving panels
    for i in eachindex(receiving)

        # control point location
        rcp = controlpoint(receiving[i])

        # normal vector body axis
        nhat = normal(receiving[i])

        if reuse_edges
            # use more efficient loop when we can reuse edges

            # calculate influence of all panels except right and trailing edges
            for j2 in 1:ns-1, j1 in 1:nc-1

                # compute induced velocity
                Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, sending[j1, j2];
                    finite_core = finite_core,
                    symmetric = symmetric,
                    bottom = false, reflected_bottom = false,
                    right = false, reflected_right = false)

                # add partial contribution from current panel
                AIC[i,ls[j1, j2]] += dot(Vhat, nhat)

                # add partial contribution from panel above this one (if applicable)
                if j1 > 1
                    AIC[i,ls[j1-1, j2]] += dot(Vhat_t, nhat)
                end

                # add partial contribution from panel to the left (if applicable)
                if j2 > 1
                    AIC[i,ls[j1, j2-1]] += dot(Vhat_l, nhat)
                end
            end

            # panels on the right edge (excluding the bottom right corner)
            j2 = ns
            for j1 in 1:nc-1

                # compute induced velocity
                Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, sending[j1, j2];
                    finite_core = finite_core,
                    symmetric = symmetric,
                    bottom = false, reflected_bottom = false)

                # add partial contribution from current panel
                AIC[i,ls[j1, j2]] += dot(Vhat, nhat)

                # add partial contribution from panel above this one
                if j1 > 1
                    AIC[i,ls[j1-1, j2]] += dot(Vhat_t, nhat)
                end

                # add partial contribution from panel to the left
                if j2 > 1
                    AIC[i,ls[j1, j2-1]] += dot(Vhat_l, nhat)
                end
            end

            # panels on the trailing edge (excluding bottom right corner)
            j1 = nc
            for j2 in 1:ns-1

                # compute induced velocity
                Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, sending[j1, j2];
                    finite_core = finite_core,
                    symmetric = symmetric,
                    bottom = !trailing_vortices,
                    reflected_bottom = !trailing_vortices,
                    right = false,
                    reflected_right = false,
                    left_trailing = trailing_vortices,
                    reflected_left_trailing = trailing_vortices,
                    right_trailing = false,
                    reflected_right_trailing = false)

                # add partial contribution from current panel
                AIC[i,ls[j1, j2]] += dot(Vhat, nhat)

                # add partial contribution from panel above this one
                if j1 > 1
                    AIC[i,ls[j1-1, j2]] += dot(Vhat_t, nhat)
                end

                # add partial contribution from panel to the left
                if j2 > 1
                    AIC[i,ls[j1, j2-1]] += dot(Vhat_l, nhat)
                end

            end

            # bottom right corner
            j1 = nc
            j2 = ns

            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, sending[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = !trailing_vortices,
                reflected_bottom = !trailing_vortices,
                left_trailing = trailing_vortices,
                reflected_left_trailing = trailing_vortices,
                right_trailing = trailing_vortices,
                reflected_right_trailing = trailing_vortices)

            # add partial contribution from current panel
            AIC[i,ls[j1, j2]] += dot(Vhat, nhat)

            # add partial contribution from panel above this one
            if j1 > 1
                AIC[i,ls[j1-1, j2]] += dot(Vhat_t, nhat)
            end

            # add partial contribution from panel to the left
            if j2 > 1
                AIC[i,ls[j1, j2-1]] += dot(Vhat_l, nhat)
            end

        else
            # we can't reuse edges, probably because the finite-core model is active

            # calculate influence of all panels except trailing edge panels
            for j2 in 1:ns, j1 in 1:nc-1
                # compute induced velocity
                Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, sending[j1, j2];
                    finite_core = finite_core,
                    symmetric = symmetric)

                # add to AIC matrix
                j = ls[j1, j2]
                AIC[i,j] += dot(Vhat, nhat)
            end

            # calculate influence of all trailing edge panels
            j1 = nc
            for j2 in 1:ns
                Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, sending[j1, j2];
                    finite_core = finite_core,
                    symmetric = symmetric,
                    bottom = !trailing_vortices,
                    reflected_bottom = !trailing_vortices,
                    left_trailing = trailing_vortices,
                    reflected_left_trailing = trailing_vortices,
                    right_trailing = trailing_vortices,
                    reflected_right_trailing = trailing_vortices)

                # add to AIC matrix
                j = ls[j1, j2]
                AIC[i,j] += dot(Vhat, nhat)
            end
        end
    end

    return AIC
end

"""
    induced_velocity(rcp, surface::AbstractMatrix{<:VortexRing}, Γ; kwargs...)

Compute the velocity induced by the grid of vortex ring panels in `surface` at
control point `rcp`

# Keyword Arguments
 - `finite_core`: Flag indicating whether the finite core model should be used
 - `symmetric`: Flag indicating whether sending panels should be mirrored across the X-Z plane
 - `trailing_vortices`: Indicates whether trailing vortices are used. Defaults to `true`
 - `xhat`: direction in which trailing vortices are shed if `trailing_vortices = true`.
    Defaults to [1, 0, 0]
"""
@inline function induced_velocity(rcp, surface, Γ; finite_core = true,
    symmetric = false, trailing_vortices = true, xhat = SVector(1, 0, 0))

    @assert length(surface) == length(Γ)

    # initialize output
    Vind = zero(rcp)

    # get surface surface dimensions
    nc, ns = size(surface)
    ls = LinearIndices(surface)

    # determine control flow parameters based on inputs
    reuse_edges = !finite_core

    if reuse_edges
        # use more efficient loop when we can reuse edges

        # calculate influence of all panels except right and trailing edges
        for j2 in 1:ns-1, j1 in 1:nc-1

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = false,
                reflected_bottom = false,
                right = false,
                reflected_right = false)

            # add partial contribution from current panel
            Vind += Vhat * Γ[ls[j1, j2]]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * Γ[ls[j1-1, j2]]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * Γ[ls[j1, j2-1]]
            end
        end

        # panels on the right edge (excluding the bottom right corner)
        j2 = ns
        for j1 in 1:nc-1

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = false,
                reflected_bottom = false)

            # add partial contribution from current panel
            Vind += Vhat * Γ[ls[j1, j2]]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * Γ[ls[j1-1, j2]]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * Γ[ls[j1, j2-1]]
            end
        end

        # panels on the trailing edge (exculding the bottom right corner)
        j1 = nc
        for j2 in 1:ns-1

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                right = false,
                reflected_right = false,
                bottom = !trailing_vortices,
                reflected_bottom = !trailing_vortices,
                left_trailing = trailing_vortices,
                reflected_left_trailing = trailing_vortices,
                right_trailing = false,
                reflected_right_trailing = false)

            # add partial contribution from current panel
            Vind += Vhat * Γ[ls[j1, j2]]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * Γ[ls[j1-1, j2]]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * Γ[ls[j1, j2-1]]
            end
        end

        # bottom right corner
        j1 = nc
        j2 = ns

        Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
            finite_core = finite_core,
            symmetric = symmetric,
            bottom = !trailing_vortices,
            reflected_bottom = !trailing_vortices,
            left_trailing = trailing_vortices,
            reflected_left_trailing = trailing_vortices,
            right_trailing = trailing_vortices,
            reflected_right_trailing = trailing_vortices)

        # add partial contribution from current panel
        Vind += Vhat * Γ[ls[j1, j2]]

        # add partial contribution from panel above this one (if applicable)
        if j1 > 1
            Vind += Vhat_t * Γ[ls[j1-1, j2]]
        end

        # add partial contribution from panel to the left (if applicable)
        if j2 > 1
            Vind += Vhat_l * Γ[ls[j1, j2-1]]
        end

    else
        # we can't reuse edges, probably because the finite-core model is active

        # calculate influence of all panels except trailing edge panels
        for j2 in 1:ns, j1 in 1:nc-1
            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric)

            # add contribution to induced velocity
            Vind += Vhat * Γ[ls[j1, j2]]
        end

        # calculate influence of all trailing edge panels
        j1 = nc
        for j2 in 1:ns
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = !trailing_vortices,
                reflected_bottom = !trailing_vortices,
                left_trailing = trailing_vortices,
                reflected_left_trailing = trailing_vortices,
                right_trailing = trailing_vortices,
                reflected_right_trailing = trailing_vortices)

            # add contribution to induced velocity
            Vind += Vhat * Γ[ls[j1, j2]]
        end
    end

    return Vind
end

"""
    induced_velocity(I::CartesianIndex, surface::AbstractMatrix{<:VortexRing}, Γ; kwargs...)

Compute the velocity induced by the grid of vortex ring panels in `surface` on
the top bound vortex of panel `I` in `surface`.

# Keyword Arguments
 - `finite_core`: Flag indicating whether the finite core model should be used
 - `symmetric`: Flag indicating whether sending panels should be mirrored across the X-Z plane
 - `trailing_vortices`: Indicates whether trailing vortices are used. Defaults to `true`.
 - `xhat`: Direction in which trailing vortices are shed if `trailing_vortices = true`.
    Defaults to [1, 0, 0]
"""
@inline function induced_velocity(I::CartesianIndex, surface, Γ; finite_core = true,
    symmetric = false, trailing_vortices = true, xhat = SVector(1, 0, 0))

    @assert length(surface) == length(Γ)

    # extract the control point of interest
    rcp = top_center(surface[I])

    # initialize output
    Vind = zero(rcp)

    # get surface dimensions
    nc, ns = size(surface)
    ls = LinearIndices(surface)

    # determine control flow parameters based on inputs
    reuse_edges = !finite_core

    if reuse_edges
        # use more efficient loop when we can reuse edges

        # calculate influence of all panels except right and trailing edges
        for j2 in 1:ns-1, j1 in 1:nc-1

            J = CartesianIndex(j1, j2)
            include_top = J != I

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                bottom = false,
                reflected_bottom = false,
                right = false,
                reflected_right = false)

            # add partial contribution from current panel
            Vind += Vhat * Γ[ls[j1, j2]]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * Γ[ls[j1-1, j2]]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * Γ[ls[j1, j2-1]]
            end
        end

        # panels on the right edge (excluding the bottom right corner)
        j2 = ns
        for j1 in 1:nc-1

            J = CartesianIndex(j1, j2)
            include_top = J != I

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                bottom = false,
                reflected_bottom = false)

            # add partial contribution from current panel
            Vind += Vhat * Γ[ls[j1, j2]]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * Γ[ls[j1-1, j2]]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * Γ[ls[j1, j2-1]]
            end
        end

        # panels on the trailing edge (excluding the bottom right corner)
        j1 = nc
        for j2 in 1:ns-1

            J = CartesianIndex(j1, j2)
            include_top = J != I

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                bottom = !trailing_vortices,
                reflected_bottom = !trailing_vortices,
                right = false,
                reflected_right = false,
                left_trailing = trailing_vortices,
                reflected_left_trailing = trailing_vortices,
                right_trailing = false,
                reflected_right_trailing = false)

            # add partial contribution from current panel
            Vind += Vhat * Γ[ls[j1, j2]]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * Γ[ls[j1-1, j2]]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * Γ[ls[j1, j2-1]]
            end
        end

        # bottom right corner
        j1 = nc
        j2 = ns

        J = CartesianIndex(j1, j2)
        include_top = J != I

        Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
            finite_core = finite_core,
            symmetric = symmetric,
            top = include_top,
            bottom = !trailing_vortices,
            reflected_bottom = !trailing_vortices,
            left_trailing = trailing_vortices,
            reflected_left_trailing = trailing_vortices,
            right_trailing = trailing_vortices,
            reflected_right_trailing = trailing_vortices)

        # add partial contribution from current panel
        Vind += Vhat * Γ[ls[j1, j2]]

        # add partial contribution from panel above this one (if applicable)
        if j1 > 1
            Vind += Vhat_t * Γ[ls[j1-1, j2]]
        end

        # add partial contribution from panel to the left (if applicable)
        if j2 > 1
            Vind += Vhat_l * Γ[ls[j1, j2-1]]
        end
    else
        # we can't reuse edges, probably because the finite-core model is active

        # calculate influence of all panels except trailing edge panels
        for j2 in 1:ns, j1 in 1:nc-1
            J = CartesianIndex(j1, j2)
            include_top = J != I
            include_bottom = J != CartesianIndex(I[1]+1, I[2])

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                bottom = include_bottom)

            # add contribution to induced velocity
            Vind += Vhat * Γ[ls[j1, j2]]
        end

        # calculate influence of all trailing edge panels
        j1 = nc
        for j2 in 1:ns
            J = CartesianIndex(j1, j2)
            include_top = J != I

            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                bottom = !trailing_vortices,
                reflected_bottom = !trailing_vortices,
                left_trailing = trailing_vortices,
                reflected_left_trailing = trailing_vortices,
                right_trailing = trailing_vortices,
                reflected_right_trailing = trailing_vortices)

            # add contribution to induced velocity
            Vind += Vhat * Γ[ls[j1, j2]]
        end
    end

    return Vind
end

"""
    induced_velocity(is::Integer, surface, Γ; kwargs...)

Compute the induced velocity from the grid of vortex ring panels in `surface` at
the trailing edge vertex corresponding to index `is`

# Keyword Arguments
 - `finite_core`: Flag indicating whether the finite core model should be used
 - `symmetric`: Flag indicating whether sending panels should be mirrored across the X-Z plane
 - `trailing_vortices`: Indicates whether trailing vortices are used. Defaults to `true`.
 - `xhat`: Direction in which trailing vortices are shed if `trailing_vortices = true`.
    Defaults to [1, 0, 0]
"""
@inline function induced_velocity(is::Integer, surface, Γ; finite_core = true,
    symmetric = false, trailing_vortices = true, xhat = SVector(1, 0, 0))

    @assert length(surface) == length(Γ)

    # get surface surface dimensions
    nc, ns = size(surface)
    ls = LinearIndices(surface)

    # extract the control point of interest
    if is <= ns
        rcp = bottom_left(surface[nc, is])
    else
        rcp = bottom_right(surface[nc, is-1])
    end

    # initialize output
    Vind = zero(rcp)

    # determine control flow parameters based on inputs
    reuse_edges = !finite_core
    keep_mirrored = not_on_symmetry_plane(rcp)

    if reuse_edges
        # use more efficient loop when we can reuse edges

        # calculate influence of all panels except right and trailing edges
        for j2 in 1:ns-1, j1 in 1:nc-1

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = false,
                reflected_bottom = false,
                right = false,
                reflected_right = false)

            # add partial contribution from current panel
            Vind += Vhat * Γ[ls[j1, j2]]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * Γ[ls[j1-1, j2]]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * Γ[ls[j1, j2-1]]
            end
        end

        # panels on the right edge (excluding the bottom right corner
        j2 = ns
        for j1 in 1:nc-1

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = false,
                reflected_bottom = false)

            # add partial contribution from current panel
            Vind += Vhat * Γ[ls[j1, j2]]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * Γ[ls[j1-1, j2]]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * Γ[ls[j1, j2-1]]
            end
        end

        # panels on the trailing edge (excluding the bottom right corner)
        j1 = nc
        for j2 in 1:ns-1

            include_left = j2 != is
            include_right = j2 != is - 1
            include_bottom = include_left && include_right

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = !trailing_vortices && include_bottom,
                reflected_bottom = !trailing_vortices && (include_bottom || keep_mirrored),
                left = include_left,
                reflected_left = include_left || keep_mirrored,
                right = false,
                reflected_right = false,
                left_trailing = trailing_vortices && include_left,
                reflected_left_trailing = trailing_vortices && (include_left || keep_mirrored),
                right_trailing = false,
                reflected_right_trailing = false)

            # add partial contribution from current panel
            Vind += Vhat * Γ[ls[j1, j2]]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * Γ[ls[j1-1, j2]]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * Γ[ls[j1, j2-1]]
            end
        end

        # bottom right corner
        j1 = nc
        j2 = ns

        include_left = j2 != is
        include_right = j2 != is - 1
        include_bottom = include_left && include_right

        # compute induced velocity
        Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
            finite_core = finite_core,
            symmetric = symmetric,
            bottom = !trailing_vortices && include_bottom,
            reflected_bottom = !trailing_vortices && (include_bottom || keep_mirrored),
            left = include_left,
            reflected_left = include_left || keep_mirrored,
            right = include_right,
            reflected_right = include_right || keep_mirrored,
            left_trailing = trailing_vortices && include_left,
            reflected_left_trailing = trailing_vortices && (include_left || keep_mirrored),
            right_trailing = trailing_vortices && include_left,
            reflected_right_trailing = trailing_vortices && (include_left || keep_mirrored))

        # add partial contribution from current panel
        Vind += Vhat * Γ[ls[j1, j2]]

        # add partial contribution from panel above this one (if applicable)
        if j1 > 1
            Vind += Vhat_t * Γ[ls[j1-1, j2]]
        end

        # add partial contribution from panel to the left (if applicable)
        if j2 > 1
            Vind += Vhat_l * Γ[ls[j1, j2-1]]
        end

    else
        # we can't reuse edges, probably because the finite-core model is active

        # calculate influence of all panels except trailing edge panels
        for j2 in 1:ns, j1 in 1:nc-1
            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric)

            # add contribution to induced velocity
            Vind += Vhat * Γ[ls[j1, j2]]
        end

        # calculate influence of all trailing edge panels
        j1 = nc
        for j2 in 1:ns
            include_left = j2 != is
            include_right = j2 != is - 1
            include_bottom = include_left && include_right

            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = !trailing_vortices && include_bottom,
                reflected_bottom = !trailing_vortices && (include_bottom || keep_mirrored),
                left = include_left,
                reflected_left = include_left || keep_mirrored,
                right = include_right,
                reflected_right = include_right || keep_mirrored,
                left_trailing = trailing_vortices && include_left,
                reflected_left_trailing = trailing_vortices && (include_left || keep_mirrored),
                right_trailing = trailing_vortices && include_right,
                reflected_right_trailing = trailing_vortices && (include_right || keep_mirrored))

            # add contribution to induced velocity
            Vind += Vhat * Γ[ls[j1, j2]]
        end
    end

    return Vind
end

"""
    induced_velocity_derivatives(rcp, surface::AbstractMatrix{<:VortexRing}, Γ, dΓ; kwargs...)

Compute the velocity induced by the grid of vortex ring panels in `surface` at
control point `rcp` and its derivatives with respect to the freestream variables

# Keyword Arguments
 - `finite_core`: Flag indicating whether the finite core model should be used
 - `symmetric`: Flag indicating whether sending panels should be mirrored across the X-Z plane
 - `trailing_vortices`: Indicates whether trailing vortices are used. Defaults to `true`.
 - `xhat`: Direction in which trailing vortices are shed if `trailing_vortices = true`.
    Defaults to [1, 0, 0]
"""
@inline function induced_velocity_derivatives(rcp, surface, Γ, dΓ; finite_core = true,
    symmetric = false, trailing_vortices = true, xhat = SVector(1, 0, 0))

    # unpack derivatives
    Γ_a, Γ_b, Γ_p, Γ_q, Γ_r = dΓ

    # check input size
    @assert length(surface) == length(Γ) == length(Γ_a) == length(Γ_b) ==
        length(Γ_p) == length(Γ_q) == length(Γ_r)

    # initialize outputs
    Vind = zero(rcp)

    Vind_a = zero(rcp)
    Vind_b = zero(rcp)
    Vind_p = zero(rcp)
    Vind_q = zero(rcp)
    Vind_r = zero(rcp)

    # get surface surface dimensions
    nc, ns = size(surface)
    ls = LinearIndices(surface)

    # determine control flow parameters based on inputs
    reuse_edges = !finite_core

    if reuse_edges
        # use more efficient loop when we can reuse edges

        # calculate influence of all panels except right and trailing edges
        for j2 in 1:ns-1, j1 in 1:nc-1

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = false,
                reflected_bottom = false,
                right = false,
                reflected_right = false)

            # add partial contribution from current panel
            j = ls[j1, j2]

            Vind += Vhat * Γ[j]

            Vind_a += Vhat*Γ_a[j]
            Vind_b += Vhat*Γ_b[j]
            Vind_p += Vhat*Γ_p[j]
            Vind_q += Vhat*Γ_q[j]
            Vind_r += Vhat*Γ_r[j]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                j = ls[j1-1, j2]

                Vind += Vhat_t * Γ[j]

                Vind_a += Vhat_t*Γ_a[j]
                Vind_b += Vhat_t*Γ_b[j]
                Vind_p += Vhat_t*Γ_p[j]
                Vind_q += Vhat_t*Γ_q[j]
                Vind_r += Vhat_t*Γ_r[j]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                j = ls[j1, j2-1]

                Vind += Vhat_l * Γ[j]

                Vind_a += Vhat_l*Γ_a[j]
                Vind_b += Vhat_l*Γ_b[j]
                Vind_p += Vhat_l*Γ_p[j]
                Vind_q += Vhat_l*Γ_q[j]
                Vind_r += Vhat_l*Γ_r[j]
            end
        end

        # panels on the right edge (excluding the bottom right corner)
        j2 = ns
        for j1 in 1:nc-1

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = false,
                reflected_bottom = false)

            # add partial contribution from current panel
            j = ls[j1, j2]

            Vind += Vhat * Γ[j]

            Vind_a += Vhat*Γ_a[j]
            Vind_b += Vhat*Γ_b[j]
            Vind_p += Vhat*Γ_p[j]
            Vind_q += Vhat*Γ_q[j]
            Vind_r += Vhat*Γ_r[j]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                j = ls[j1-1, j2]

                Vind += Vhat_t * Γ[j]

                Vind_a += Vhat_t*Γ_a[j]
                Vind_b += Vhat_t*Γ_b[j]
                Vind_p += Vhat_t*Γ_p[j]
                Vind_q += Vhat_t*Γ_q[j]
                Vind_r += Vhat_t*Γ_r[j]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                j = ls[j1, j2-1]

                Vind += Vhat_l * Γ[j]

                Vind_a += Vhat_l*Γ_a[j]
                Vind_b += Vhat_l*Γ_b[j]
                Vind_p += Vhat_l*Γ_p[j]
                Vind_q += Vhat_l*Γ_q[j]
                Vind_r += Vhat_l*Γ_r[j]
            end
        end

        # panels on the trailing edge (excluding the bottom right corner)
        j1 = nc
        for j2 in 1:ns-1

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                right = false,
                reflected_right = false,
                bottom = !trailing_vortices,
                reflected_bottom = !trailing_vortices,
                left_trailing = trailing_vortices,
                reflected_left_trailing = trailing_vortices,
                right_trailing = false,
                reflected_right_trailing = false)

            # add partial contribution from current panel
            j = ls[j1, j2]

            Vind += Vhat * Γ[j]

            Vind_a += Vhat*Γ_a[j]
            Vind_b += Vhat*Γ_b[j]
            Vind_p += Vhat*Γ_p[j]
            Vind_q += Vhat*Γ_q[j]
            Vind_r += Vhat*Γ_r[j]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                j = ls[j1-1, j2]

                Vind += Vhat_t * Γ[ls[j1-1, j2]]

                Vind_a += Vhat_t*Γ_a[j]
                Vind_b += Vhat_t*Γ_b[j]
                Vind_p += Vhat_t*Γ_p[j]
                Vind_q += Vhat_t*Γ_q[j]
                Vind_r += Vhat_t*Γ_r[j]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                j = ls[j1, j2-1]

                Vind += Vhat_l * Γ[j]

                Vind_a += Vhat_l*Γ_a[j]
                Vind_b += Vhat_l*Γ_b[j]
                Vind_p += Vhat_l*Γ_p[j]
                Vind_q += Vhat_l*Γ_q[j]
                Vind_r += Vhat_l*Γ_r[j]
            end
        end

        # bottom right corner
        j1 = nc
        j2 = ns

        # compute induced velocity
        Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
            finite_core = finite_core,
            symmetric = symmetric,
            bottom = !trailing_vortices,
            reflected_bottom = !trailing_vortices,
            left_trailing = trailing_vortices,
            reflected_left_trailing = trailing_vortices,
            right_trailing = trailing_vortices,
            reflected_right_trailing = trailing_vortices)

        # add partial contribution from current panel
        j = ls[j1, j2]

        Vind += Vhat * Γ[j]

        Vind_a += Vhat*Γ_a[j]
        Vind_b += Vhat*Γ_b[j]
        Vind_p += Vhat*Γ_p[j]
        Vind_q += Vhat*Γ_q[j]
        Vind_r += Vhat*Γ_r[j]

        # add partial contribution from panel above this one (if applicable)
        if j1 > 1
            j = ls[j1-1, j2]

            Vind += Vhat_t * Γ[j]

            Vind_a += Vhat_t*Γ_a[j]
            Vind_b += Vhat_t*Γ_b[j]
            Vind_p += Vhat_t*Γ_p[j]
            Vind_q += Vhat_t*Γ_q[j]
            Vind_r += Vhat_t*Γ_r[j]
        end

        # add partial contribution from panel to the left (if applicable)
        if j2 > 1
            j = ls[j1, j2-1]

            Vind += Vhat_l * Γ[j]

            Vind_a += Vhat_l*Γ_a[j]
            Vind_b += Vhat_l*Γ_b[j]
            Vind_p += Vhat_l*Γ_p[j]
            Vind_q += Vhat_l*Γ_q[j]
            Vind_r += Vhat_l*Γ_r[j]
        end

    else
        # we can't reuse edges, probably because the finite-core model is active

        # calculate influence of all panels except trailing edge panels
        for j2 in 1:ns, j1 in 1:nc-1
            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric)

            # add contribution to induced velocity
            j = ls[j1, j2]

            Vind += Vhat * Γ[j]

            Vind_a += Vhat*Γ_a[j]
            Vind_b += Vhat*Γ_b[j]
            Vind_p += Vhat*Γ_p[j]
            Vind_q += Vhat*Γ_q[j]
            Vind_r += Vhat*Γ_r[j]
        end

        # calculate influence of all trailing edge panels
        j1 = nc
        for j2 in 1:ns
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = !trailing_vortices,
                reflected_bottom = !trailing_vortices,
                left_trailing = trailing_vortices,
                reflected_left_trailing = trailing_vortices,
                right_trailing = trailing_vortices,
                reflected_right_trailing = trailing_vortices)

            # add contribution to induced velocity
            j = ls[j1, j2]

            Vind += Vhat * Γ[j]

            Vind_a += Vhat*Γ_a[j]
            Vind_b += Vhat*Γ_b[j]
            Vind_p += Vhat*Γ_p[j]
            Vind_q += Vhat*Γ_q[j]
            Vind_r += Vhat*Γ_r[j]
        end
    end

    dVind = (Vind_a, Vind_b, Vind_p, Vind_q, Vind_r)

    return Vind, dVind
end

"""
    induced_velocity_derivatives(I::CartesianIndex, surface::AbstractMatrix{<:Wake}; kwargs...)

Compute the velocity induced by the grid of vortex ring panels in `surface` on
the top bound vortex of panel `I` in `surface` and its derivatives with respect
to the freestream variables

# Keyword Arguments
 - `finite_core`: Flag indicating whether the finite core model should be used
 - `symmetric`: Flag indicating whether sending panels should be mirrored across the X-Z plane
 - `trailing_vortices`: Indicates whether trailing vortices are used. Defaults to `true`.
 - `xhat`: Direction in which trailing vortices are shed if `trailing_vortices = true`.
    Defaults to [1, 0, 0]
"""
@inline function induced_velocity_derivatives(I::CartesianIndex, surface, Γ, dΓ;
    finite_core = true, symmetric = false, trailing_vortices = true, xhat = SVector(1, 0, 0))

    # unpack derivatives
    Γ_a, Γ_b, Γ_p, Γ_q, Γ_r = dΓ

    # check input size
    @assert length(surface) == length(Γ) == length(Γ_a) == length(Γ_b) ==
        length(Γ_p) == length(Γ_q) == length(Γ_r)

    # extract the control point of interest
    rcp = top_center(surface[I])

    # initialize outputs
    Vind = zero(rcp)

    Vind_a = zero(rcp)
    Vind_b = zero(rcp)
    Vind_p = zero(rcp)
    Vind_q = zero(rcp)
    Vind_r = zero(rcp)

    # get surface surface dimensions
    nc, ns = size(surface)
    ls = LinearIndices(surface)

    # determine control flow parameters based on inputs
    reuse_edges = !finite_core

    if reuse_edges
        # use more efficient loop when we can reuse edges

        # calculate influence of all panels except right and trailing edges
        for j2 in 1:ns-1, j1 in 1:nc-1

            J = CartesianIndex(j1, j2)
            include_top = J != I

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                bottom = false,
                reflected_bottom = false,
                right = false,
                reflected_right = false)

            # add partial contribution from current panel
            j = ls[j1, j2]

            Vind += Vhat * Γ[j]

            Vind_a += Vhat*Γ_a[j]
            Vind_b += Vhat*Γ_b[j]
            Vind_p += Vhat*Γ_p[j]
            Vind_q += Vhat*Γ_q[j]
            Vind_r += Vhat*Γ_r[j]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                j = ls[j1-1, j2]

                Vind += Vhat_t * Γ[j]

                Vind_a += Vhat_t*Γ_a[j]
                Vind_b += Vhat_t*Γ_b[j]
                Vind_p += Vhat_t*Γ_p[j]
                Vind_q += Vhat_t*Γ_q[j]
                Vind_r += Vhat_t*Γ_r[j]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                j = ls[j1, j2-1]

                Vind += Vhat_l * Γ[j]

                Vind_a += Vhat_l*Γ_a[j]
                Vind_b += Vhat_l*Γ_b[j]
                Vind_p += Vhat_l*Γ_p[j]
                Vind_q += Vhat_l*Γ_q[j]
                Vind_r += Vhat_l*Γ_r[j]
            end
        end

        # panels on the right edge (excluding the bottom right corner)
        j2 = ns
        for j1 in 1:nc-1

            J = CartesianIndex(j1, j2)
            include_top = J != I

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                bottom = false,
                reflected_bottom = false)

            # add partial contribution from current panel
            j = ls[j1, j2]

            Vind += Vhat * Γ[j]

            Vind_a += Vhat*Γ_a[j]
            Vind_b += Vhat*Γ_b[j]
            Vind_p += Vhat*Γ_p[j]
            Vind_q += Vhat*Γ_q[j]
            Vind_r += Vhat*Γ_r[j]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                j = ls[j1-1, j2]

                Vind += Vhat_t * Γ[j]

                Vind_a += Vhat_t*Γ_a[j]
                Vind_b += Vhat_t*Γ_b[j]
                Vind_p += Vhat_t*Γ_p[j]
                Vind_q += Vhat_t*Γ_q[j]
                Vind_r += Vhat_t*Γ_r[j]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                j = ls[j1, j2-1]

                Vind += Vhat_l * Γ[j]

                Vind_a += Vhat_l*Γ_a[j]
                Vind_b += Vhat_l*Γ_b[j]
                Vind_p += Vhat_l*Γ_p[j]
                Vind_q += Vhat_l*Γ_q[j]
                Vind_r += Vhat_l*Γ_r[j]
            end
        end

        # panels on the trailing edge (excluding bottom right corner)
        j1 = nc
        for j2 in 1:ns-1

            J = CartesianIndex(j1, j2)
            include_top = J != I

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                bottom = !trailing_vortices,
                reflected_bottom = !trailing_vortices,
                right = false,
                reflected_right = false,
                left_trailing = trailing_vortices,
                reflected_left_trailing = trailing_vortices,
                right_trailing = false,
                reflected_right_trailing = false)

            # add partial contribution from current panel
            j = ls[j1, j2]

            Vind += Vhat * Γ[j]

            Vind_a += Vhat*Γ_a[j]
            Vind_b += Vhat*Γ_b[j]
            Vind_p += Vhat*Γ_p[j]
            Vind_q += Vhat*Γ_q[j]
            Vind_r += Vhat*Γ_r[j]

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                j = ls[j1-1, j2]

                Vind += Vhat_t * Γ[j]

                Vind_a += Vhat_t*Γ_a[j]
                Vind_b += Vhat_t*Γ_b[j]
                Vind_p += Vhat_t*Γ_p[j]
                Vind_q += Vhat_t*Γ_q[j]
                Vind_r += Vhat_t*Γ_r[j]
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                j = ls[j1, j2-1]

                Vind += Vhat_l * Γ[j]

                Vind_a += Vhat_l*Γ_a[j]
                Vind_b += Vhat_l*Γ_b[j]
                Vind_p += Vhat_l*Γ_p[j]
                Vind_q += Vhat_l*Γ_q[j]
                Vind_r += Vhat_l*Γ_r[j]
            end
        end

        # bottom right corner
        j1 = nc
        j2 = ns

        J = CartesianIndex(j1, j2)
        include_top = J != I

        Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
            finite_core = finite_core,
            symmetric = symmetric,
            top = include_top,
            bottom = !trailing_vortices,
            reflected_bottom = !trailing_vortices,
            left_trailing = trailing_vortices,
            reflected_left_trailing = trailing_vortices,
            right_trailing = trailing_vortices,
            reflected_right_trailing = trailing_vortices)

        # add partial contribution from current panel
        j = ls[j1, j2]

        Vind += Vhat * Γ[j]

        Vind_a += Vhat*Γ_a[j]
        Vind_b += Vhat*Γ_b[j]
        Vind_p += Vhat*Γ_p[j]
        Vind_q += Vhat*Γ_q[j]
        Vind_r += Vhat*Γ_r[j]

        # add partial contribution from panel above this one (if applicable)
        if j1 > 1
            j = ls[j1-1, j2]

            Vind += Vhat_t * Γ[ls[j1-1, j2]]

            Vind_a += Vhat_t*Γ_a[j]
            Vind_b += Vhat_t*Γ_b[j]
            Vind_p += Vhat_t*Γ_p[j]
            Vind_q += Vhat_t*Γ_q[j]
            Vind_r += Vhat_t*Γ_r[j]
        end

        # add partial contribution from panel to the left (if applicable)
        if j2 > 1
            j = ls[j1, j2-1]

            Vind += Vhat_l * Γ[j]

            Vind_a += Vhat_l*Γ_a[j]
            Vind_b += Vhat_l*Γ_b[j]
            Vind_p += Vhat_l*Γ_p[j]
            Vind_q += Vhat_l*Γ_q[j]
            Vind_r += Vhat_l*Γ_r[j]
        end

    else
        # we can't reuse edges, probably because the finite-core model is active

        # calculate influence of all panels except trailing edge panels
        for j2 in 1:ns, j1 in 1:nc-1
            J = CartesianIndex(j1, j2)
            include_top = J != I
            include_bottom = J != CartesianIndex(I[1]-1, I[2])

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                bottom = include_bottom)

            # add contribution to induced velocity
            j = ls[j1, j2]

            Vind += Vhat * Γ[j]

            Vind_a += Vhat*Γ_a[j]
            Vind_b += Vhat*Γ_b[j]
            Vind_p += Vhat*Γ_p[j]
            Vind_q += Vhat*Γ_q[j]
            Vind_r += Vhat*Γ_r[j]
        end

        # calculate influence of all trailing edge panels
        j1 = nc
        for j2 in 1:ns
            J = CartesianIndex(j1, j2)
            include_top = J != I

            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                bottom = !trailing_vortices,
                reflected_bottom = !trailing_vortices,
                left_trailing = trailing_vortices,
                reflected_left_trailing = trailing_vortices,
                right_trailing = trailing_vortices,
                reflected_right_trailing = trailing_vortices)

            # add contribution to induced velocity
            j = ls[j1, j2]

            Vind += Vhat * Γ[j]

            Vind_a += Vhat*Γ_a[j]
            Vind_b += Vhat*Γ_b[j]
            Vind_p += Vhat*Γ_p[j]
            Vind_q += Vhat*Γ_q[j]
            Vind_r += Vhat*Γ_r[j]
        end
    end

    dVind = (Vind_a, Vind_b, Vind_p, Vind_q, Vind_r)

    return Vind, dVind
end

"""
    induced_velocity(rcp, surface::AbstractMatrix{<:Wake}; kwargs...)

Compute the velocity induced by the grid of vortex ring panels in `surface` at
control point `rcp`

# Keyword Arguments
 - `finite_core`: Flag indicating whether the finite core model should be used
 - `symmetric`: Flag indicating whether sending panels should be mirrored across the X-Z plane
 - `nwake`: number of chordwise wake panels to use from each `wake`
 - `trailing_vortices`: Indicates whether trailing vortices are used. Defaults to `true`.
 - `xhat`: Direction in which trailing vortices are shed if `trailing_vortices = true`.
    Defaults to [1, 0, 0]
"""
@inline function induced_velocity(rcp, surface; finite_core = true, symmetric = false,
    nwake = size(surface, 1), trailing_vortices = true, xhat = SVector(1, 0, 0))

    # initialize output
    Vind = zero(rcp)

    # automatically return if there is nothing to compute
    if iszero(nwake)
        return Vind
    end

    # get surface surface dimensions
    nw = nwake
    ns = size(surface, 2)

    # determine control flow parameters based on inputs
    reuse_edges = !finite_core

    if reuse_edges
        # use more efficient loop when we can reuse edges

        # calculate influence of all panels except right and trailing edges
        for j2 in 1:ns-1, j1 in 1:nw-1

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = false,
                reflected_bottom = false,
                right = false,
                reflected_right = false)

            # add partial contribution from current panel
            Vind += Vhat * surface[j1, j2].gamma

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * surface[j1-1, j2].gamma
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * surface[j1, j2-1].gamma
            end
        end

        # panels on the right edge (excluding the trailing edge)
        j2 = ns
        for j1 in 1:nw-1

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = false,
                reflected_bottom = false)

            # add partial contribution from current panel
            Vind += Vhat * surface[j1, j2].gamma

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * surface[j1-1, j2].gamma
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * surface[j1, j2-1].gamma
            end
        end

        # panels on the trailing edge
        j1 = nw
        for j2 in 1:ns

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = !trailing_vortices,
                reflected_bottom = !trailing_vortices,
                left_trailing = trailing_vortices,
                reflected_left_trailing = trailing_vortices,
                right_trailing = trailing_vortices,
                reflected_right_trailing = trailing_vortices)

            # add partial contribution from current panel
            Vind += Vhat * surface[j1, j2].gamma

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * surface[j1-1, j2].gamma
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * surface[j1, j2-1].gamma
            end
        end

    else
        # we can't reuse edges, probably because the finite-core model is active

        # calculate influence of all panels except trailing edge panels
        for j2 in 1:ns, j1 in 1:nw-1
            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric)

            # add contribution to induced velocity
            Vind += Vhat * surface[j1, j2].gamma
        end

        # calculate influence of all trailing edge panels
        j1 = nw
        for j2 in 1:ns
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                bottom = !trailing_vortices,
                reflected_bottom = !trailing_vortices,
                left_trailing = trailing_vortices,
                reflected_left_trailing = trailing_vortices,
                right_trailing = trailing_vortices,
                reflected_right_trailing = trailing_vortices)

            # add contribution to induced velocity
            Vind += Vhat * surface[j1, j2].gamma

        end
    end

    return Vind
end

"""
    induced_velocity(I::CartesianIndex, surface; kwargs...)

Compute the induced velocity from the grid of vortex ring wake panels in `surface`
at the vertex corresponding to index `I`

# Keyword Arguments
 - `finite_core`: Flag indicating whether the finite core model should be used
 - `symmetric`: Flag indicating whether sending panels should be mirrored across the X-Z plane
 - `nwake`: number of chordwise wake panels to use from each `wake`
 - `trailing_vortices`: Indicates whether trailing vortices are used. Defaults to `true`.
 - `xhat`: Direction in which trailing vortices are shed if `trailing_vortices = true`.
    Defaults to [1, 0, 0]
"""
@inline function induced_velocity(I::CartesianIndex, surface; finite_core = true,
    symmetric = false, nwake = size(surface, 1), trailing_vortices = true, xhat = SVector(1, 0, 0))

    # get surface surface dimensions
    nw = nwake
    ns = size(surface, 2)

    # initialize output
    Vind = @SVector zeros(eltype(eltype(surface)), 3)

    # automatically return if there is nothing to compute
    if iszero(nwake)
        return Vind
    end

    # extract the control point of interest
    if I[1] <= nw && I[2] <= ns
        rcp = top_left(surface[I[1], I[2]])
    elseif I[1] <= nw && I[2] == ns + 1
        rcp = top_right(surface[I[1], I[2]-1])
    elseif I[1] == nw + 1 && I[2] <= ns
        rcp = bottom_left(surface[I[1]-1, I[2]])
    else # I[1] == nw + 1 && I[2] == ns + 1
        rcp = bottom_right(surface[I[1]-1, I[2]-1])
    end

    # determine control flow parameters based on inputs
    reuse_edges = !finite_core
    keep_mirrored = not_on_symmetry_plane(rcp)

    # set panel coordinates to skip
    skip_top = (I, CartesianIndex(I[1], I[2]-1))
    skip_bottom = (CartesianIndex(I[1]-1, I[2]), CartesianIndex(I[1]-1, I[2]-1))
    skip_left = (I, CartesianIndex(I[1]-1, I[2]))
    skip_right = (CartesianIndex(I[1], I[2]-1), CartesianIndex(I[1]-1, I[2]-1))
    skip_left_trailing = (CartesianIndex(I[1]-1, I[2]),)
    skip_right_trailing = (CartesianIndex(I[1]-1, I[2]-1),)

    if reuse_edges
        # use more efficient loop when we can reuse edges

        # calculate influence of all panels except right and trailing edges
        for j2 in 1:ns-1, j1 in 1:nw-1

            J = CartesianIndex(j1, j2)
            include_top = !(J in skip_top)
            include_left = !(J in skip_left)

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                reflected_top = include_top || keep_mirrored,
                left = include_left,
                reflected_left = include_left || keep_mirrored,
                bottom = false,
                reflected_bottom = false,
                right = false,
                reflected_right = false)

            # add partial contribution from current panel
            Vind += Vhat * surface[j1, j2].gamma

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * surface[j1-1, j2].gamma
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * surface[j1, j2-1].gamma
            end
        end

        # panels on the right edge (excluding the bottom right corner)
        j2 = ns
        for j1 in 1:nw-1

            J = CartesianIndex(j1, j2)
            include_top = !(J in skip_top)
            include_left = !(J in skip_left)
            include_right = !(J in skip_right)

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                reflected_top = include_top || keep_mirrored,
                bottom = false,
                reflected_bottom = false,
                left = include_left,
                reflected_left = include_left || keep_mirrored,
                reflected_right = include_right || keep_mirrored)

            # add partial contribution from current panel
            Vind += Vhat * surface[j1, j2].gamma

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * surface[j1-1, j2].gamma
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * surface[j1, j2-1].gamma
            end
        end

        # panels on the trailing edge (excluding the bottom right corner)
        j1 = nw
        for j2 in 1:ns-1

            J = CartesianIndex(j1, j2)
            include_top = !(J in skip_top)
            include_bottom = !(J in skip_bottom)
            include_left = !(J in skip_left)
            include_right = !(J in skip_right)
            include_left_trailing = !(J in skip_left_trailing)
            include_right_trailing = !(J in skip_right_trailing)

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                reflected_top = include_top || keep_mirrored,
                bottom = !trailing_vortices && include_bottom,
                reflected_bottom = !trailing_vortices && (include_bottom || keep_mirrored),
                left = include_left,
                reflected_left = include_left || keep_mirrored,
                right = include_right,
                reflected_right = include_right || keep_mirrored,
                left_trailing = trailing_vortices && include_left_trailing,
                reflected_left_trailing = trailing_vortices && (include_left_trailing || keep_mirrored),
                right_trailing = trailing_vortices && include_right_trailing,
                reflected_right_trailing = trailing_vortices && (include_right_trailing || keep_mirrored))

            # add partial contribution from current panel
            Vind += Vhat * surface[j1, j2].gamma

            # add partial contribution from panel above this one (if applicable)
            if j1 > 1
                Vind += Vhat_t * surface[j1-1, j2].gamma
            end

            # add partial contribution from panel to the left (if applicable)
            if j2 > 1
                Vind += Vhat_l * surface[j1, j2-1].gamma
            end
        end
    else
        # we can't reuse edges, probably because the finite-core model is active

        # calculate influence of all panels except trailing edge panels
        for j2 in 1:ns, j1 in 1:nw-1

            J = CartesianIndex(j1, j2)
            include_top = !(J in skip_top)
            include_bottom = !(J in skip_bottom)
            include_left = !(J in skip_left)
            include_right = !(J in skip_right)

            # compute induced velocity
            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                reflected_top = include_top || keep_mirrored,
                bottom = include_bottom,
                reflected_bottom = include_bottom || keep_mirrored,
                left = include_left,
                reflected_left = include_left || keep_mirrored,
                right = include_right,
                reflected_right = include_right || keep_mirrored)

            # add contribution to induced velocity
            Vind += Vhat * surface[j1, j2].gamma
        end

        # calculate influence of all trailing edge panels
        j1 = nw
        for j2 in 1:ns
            J = CartesianIndex(j1, j2)
            include_top = !(J in skip_top)
            include_bottom = !(J in skip_bottom)
            include_left = !(J in skip_left)
            include_right = !(J in skip_right)
            include_left_trailing = !(J in skip_left_trailing)
            include_right_trailing = !(J in skip_right_trailing)

            Vhat, Vhat_t, Vhat_b, Vhat_l, Vhat_r = ring_induced_velocity(rcp, surface[j1, j2];
                finite_core = finite_core,
                symmetric = symmetric,
                top = include_top,
                reflected_top = include_top || keep_mirrored,
                bottom = !trailing_vortices && include_bottom,
                reflected_bottom = !trailing_vortices && (include_bottom || keep_mirrored),
                left = include_left,
                reflected_left = include_left || keep_mirrored,
                right = include_right,
                reflected_right = include_right || keep_mirrored,
                left_trailing = trailing_vortices && include_left_trailing,
                reflected_left_trailing = trailing_vortices && (include_left_trailing || keep_mirrored),
                right_trailing = trailing_vortices && include_right_trailing,
                reflected_right_trailing = trailing_vortices && (include_right_trailing || keep_mirrored))

            # add contribution to induced velocity
            Vind += Vhat * surface[j1, j2].gamma
        end
    end

    return Vind
end
