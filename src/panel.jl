"""
    AbstractPanel{TF}

Supertype of vortex lattice method panel types.
"""
abstract type AbstractPanel{TF} end

# --- implementations --- #

"""
    Horseshoe{TF}

Horseshoe shaped panel element with trailing vortices that extend to [Inf, 0, 0]

**Fields**
 - `rl`: position of the left side of the horseshoe vortex
 - `rr`: position of the right side of the horseshoe vortex
 - `rcp`: position of the panel control point
 - `theta`: twist angle of the panel (radians)
"""
struct Horseshoe{TF} <: AbstractPanel{TF}
    rl::SVector{3, TF}
    rr::SVector{3, TF}
    rcp::SVector{3, TF}
    theta::TF
end

"""
    Horseshoe(rl, rr, rcp, theta)

Construct and return a horseshoe vortex panel

**Arguments**
- `rl`: position of the left side of the horseshoe vortex
- `rr`: position of the right side of the horseshoe vortex
- `rcp`: position of the panel control point
- `theta`: twist angle of the panel (radians)
"""
function Horseshoe(rl, rr, rcp, theta)
    TF = promote_type(eltype(rl), eltype(rr), eltype(rcp), typeof(theta))
    return Horseshoe{TF}(SVector{3}(rl), SVector{3}(rr), SVector{3}(rcp), theta)
end

@inline Base.eltype(::Type{Horseshoe{TF}}) where TF = TF
@inline Base.eltype(::Horseshoe{TF}) where TF = TF

# -------------------------------

"""
    Ring{TF}

Vortex ring shaped panel element.

**Fields**
 - `rl`: position of the left side of the bound vortex
 - `rr`: position of the right side of the bound vortex
 - `rln`: position of the left side of the next chordwise bound vortex
 - `rrn`: position of the right side of the next chordwise bound vortex
 - `rcp`: position of the panel control point
 - `normal`: normal vector of the panel at the panel control point
 - `trailing`: indicates whether the panel is on the trailing edge
"""
struct Ring{TF} <: AbstractPanel{TF}
    rtl::SVector{3, TF}
    rtr::SVector{3, TF}
    rbl::SVector{3, TF}
    rbr::SVector{3, TF}
    rcp::SVector{3, TF}
    normal::SVector{3, TF}
    trailing::Bool
end

"""
    Ring(rtl, rtr, rbl, rbr, rcp, normal, trailing)

Construct and return a vortex ring panel

**Arguments**
- `rtl`: position of the top left side of the horseshoe vortex
- `rtr`: position of the top right side of the horseshoe vortex
- `rbl`: position of the bottom left side of the horseshoe vortex
- `rbr`: position of the bottom right side of the horseshoe vortex
- `rcp`: position of the panel control point
- `normal`: normal vector of the panel at the panel control point
- `trailing`: indicates whether the panel is on the trailing edge
"""
function Ring(rtl, rtr, rbl, rbr, rcp, normal, trailing)
    TF = promote_type(eltype(rtl), eltype(rtr), eltype(rbl), eltype(rbr), eltype(rcp), eltype(normal))
    return Ring{TF}(rtl, rtr, rbl, rbr, rcp, normal, trailing)
end

@inline Base.eltype(::Type{Ring{TF}}) where TF = TF
@inline Base.eltype(::Ring{TF}) where TF = TF

# --- required methods --- #

"""
    controlpoint(panel::AbstractPanel)

Return the control point of `panel` (typically located at the 3/4 chord)
"""
controlpoint

@inline controlpoint(panel::Horseshoe) = panel.rcp

@inline controlpoint(panel::Ring) = panel.rcp

"""
    midpoint(panel::AbstractPanel)

Compute the bound vortex midpoint of `panel` (typically the quarter chord)
"""
midpoint

@inline midpoint(panel::Horseshoe) = (panel.rl + panel.rr)/2

@inline midpoint(panel::Ring) = (panel.rtl + panel.rtr)/2

"""
    panel_width(panel::AbstractPanel)

Return the panel width
"""
panel_width

@inline panel_width(panel::Horseshoe) = panel.rr - panel.rl

@inline panel_width(panel::Ring) = panel.rtr - panel.rtl

"""
    normal(panel::AbstractPanel)

Compute the normal vector of `panel`
"""
normal

@inline function normal(panel::Horseshoe)

    delta = panel.rr - panel.rl
    dy = delta[2]
    dz = delta[3]
    ds = sqrt(dy^2 + dz^2)

    cp, sp = dy/ds, dz/ds
    ct, st = cos(panel.theta), sin(panel.theta)
    nhat = SVector(st, -ct*sp, ct*cp)

    return nhat
end

@inline normal(panel::Ring) = panel.normal

"""
    trefftz_plane_normal(panel, xhat)

Compute the normal vector for `panel` when projected onto the Trefftz plane
(including magnitude), which is defined perpindicular to `xhat`
"""
trefftz_plane_normal

@inline function trefftz_plane_normal(panel::Horseshoe)

    delta = panel.rr - panel.rl
    dy = delta[2]
    dz = delta[3]

    nhat = SVector(0, -dz, dy)

    return nhat
end

@inline function trefftz_plane_normal(panel::Ring)
    delta = panel.rbr - panel.rbl
    dy = delta[2]
    dz = delta[3]

    nhat = SVector(0, -dz, dy)

    return nhat
end

"""
    translate(panel::AbstractPanel, r)

Return a copy of `panel` translated the distance specified by vector `r`
"""
translate

function translate(panel::Horseshoe, r)

    rl = panels[i].rl + r
    rr = panels[i].rr + r
    rcp = panels[i].rcp + r
    theta = panels[i].theta

    return Horseshoe(rl, rr, rcp, theta)
end

function translate(panel::Ring, r)

    rtl = panels.rtl + r
    rtr = panels.rtr + r
    rbl = panels.rbl + r
    rbr = panels.rbr + r
    rcp = panels.rcp + r
    normal = panels.normal
    trailing = panels.trailing

    return Ring(rtl, rtr, rbl, rbr, rcp, normal, trailing)
end

"""
    induced_velocity(rcp, panel, symmetric, args...)

Computes the normalized velocity induced at control point `rcp` from `panel`.
"""
induced_velocity

@inline function induced_velocity(rcp, panel::Horseshoe, symmetric, include_bound=true)

    r1 = rcp - panel.rl
    r2 = rcp - panel.rr

    # trailing vortices
    Vhat = trailing_induced_velocity(r1, r2)

    # bound vortex
    if include_bound
        Vhat += bound_induced_velocity(r1, r2)
    end

    # add contribution from other side (if applicable)
    if symmetric && not_on_symmetry_plane(panel.rl, panel.rr)
        # flip sign for y, but now left is right and right is left.
        r1 = rcp - flipy(panel.rr)
        r2 = rcp - flipy(panel.rl)
        Vhat += trailing_induced_velocity(r1, r2)
        Vhat += bound_induced_velocity(r1, r2)
    end

    return Vhat

end

@inline function induced_velocity(rcp, panel::Ring, symmetric, xhat=SVector(1,0,0), include_top=true,
    include_bottom=true)

    # position of control point relative to top of vortex ring
    r11 = rcp - panel.rtl
    r12 = rcp - panel.rtr

    # position of control point relative to bottom of vortex ring
    r21 = rcp - panel.rbl
    r22 = rcp - panel.rbr

    # left bound vortex
    Vhat = bound_induced_velocity(r21, r11)

    # right bound vortex
    Vhat += bound_induced_velocity(r12, r22)

    # top bound vortex
    if include_top
        Vhat += bound_induced_velocity(r11, r12)
    end

    if panel.trailing
        # append horseshoe vortex to the back of the vortex ring
        Vhat += trailing_induced_velocity(r21, r22, xhat)
        # we omit the bottom bound vortex since it is cancelled by the horseshoe vortex
    else
        # bottom bound vortex
        if include_bottom
            Vhat += bound_induced_velocity(r22, r21)
        end
    end

    # add contribution from other side (if applicable)
    if symmetric && not_on_symmetry_plane(panel.rtl, panel.rtr)
        # flip sign for y, but now left is right and right is left.
        r11 = rcp - flipy(panel.rtr) # top left
        r12 = rcp - flipy(panel.rtl) # top right
        r21 = rcp - flipy(panel.rbr) # bottom left
        r22 = rcp - flipy(panel.rbl) # bottom right
        Vhat += bound_induced_velocity(r21, r11) # left
        Vhat += bound_induced_velocity(r12, r22) # right
        Vhat += bound_induced_velocity(r11, r12) # top
        if panel.trailing
            Vhat += trailing_induced_velocity(r21, r22, xhat) # trailing
        else
            Vhat += bound_induced_velocity(r22, r21) # bottom
        end
    end

    return Vhat
end

# --- derived methods --- #

"""
    translate!(panels, r)

Translate the panels contained in `panels` the distance specified by vector `r`
"""
function translate!(panels, r)

    for i in eachindex(panels)
        panels[i] = translate(panels[i], r)
    end

    return panels
end

# --- internal methods --- #

"""
    flipy(r)

Flip sign of y-component of vector `r` (used for symmetry)
"""
@inline flipy(r) = SVector{3}(r[1], -r[2], r[3])

"""
    not_on_symmetry_plane(r1, r2, tol=eps())

Test whether points `r1` and `r2` are on the symmetry plane (y = 0)
"""
@inline function not_on_symmetry_plane(r1, r2, tol=eps())
    return !(isapprox(r1[2], 0.0, atol=tol) && isapprox(r2[2], 0.0, atol=tol))
end

"""
    trailing_induced_velocity(r1, r2, xhat=[1,0,0])

Compute the induced velocity (per unit circulation) for two vortices trailing in
the `xhat` direction, at a control point located at `r1` relative to the start of the
left trailing vortex and `r2` relative to the start of the right trailing vortex.
"""
@inline function trailing_induced_velocity(r1, r2, xhat=SVector(1, 0, 0))

    nr1 = norm(r1)
    nr2 = norm(r2)

    f1 = cross(r1, xhat)/(nr1 - r1[1])/nr1
    f2 = cross(r2, xhat)/(nr2 - r2[1])/nr2

    Vhat = (f1 - f2)/(4*pi)

    return Vhat
end

"""
    bound_induced_velocity(r1, r2)

Compute the induced velocity (per unit circulation) for a bound vortex, at a
control point located at `r1` relative to the start of the bound vortex and `r2`
relative to the end of the bound vortex
"""
@inline function bound_induced_velocity(r1, r2)

    nr1 = norm(r1)
    nr2 = norm(r2)
    f1 = cross(r1, r2)/(nr1*nr2 + dot(r1, r2))
    f2 = (1/nr1 + 1/nr2)

    Vhat = (f1*f2)/(4*pi)

    return Vhat
end
