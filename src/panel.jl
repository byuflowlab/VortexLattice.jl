"""
    AbstractPanel{TF}

Supertype of vortex lattice method panel types.
"""
abstract type AbstractPanel{TF} end

"""
    Horseshoe{TF}

Horseshoe shaped panel element with trailing vortices that extend in the +x
direction to the trailing edge, and then into the farfield.

**Fields**
 - `rl`: position of the left side of the bound vortex
 - `rc`: position of the center of the bound vortex
 - `rr`: position of the right side of the bound vortex
 - `rcp`: position of the panel control point
 - `ncp`: normal vector at the panel control point
 - `xl_te`: x-distance from the left side of the bound vortex to the trailing edge
 - `xc_te`: x-distance from the center of the bound vortex to the trailing edge
 - `xr_te`: x-distance from the right side of the bound vortex to the trailing edge
 - `core_size`: vortex core size
"""
struct Horseshoe{TF} <: AbstractPanel{TF}
    rl::SVector{3, TF}
    rc::SVector{3, TF}
    rr::SVector{3, TF}
    rcp::SVector{3, TF}
    ncp::SVector{3, TF}
    xl_te::TF
    xc_te::TF
    xr_te::TF
    core_size::TF
end

"""
    Horseshoe(rl, rr, rcp, xl_te, xc_te, xr_te, theta; rc)

Construct and return a horseshoe vortex panel.

**Arguments**
- `rl`: position of the left side of the bound vortex
- `rr`: position of the right side of the bound vortex
- `rcp`: position of the panel control point
- `ncp`: normal vector at the panel control point
- `xl_te`: x-distance from the left side of the bound vortex to the trailing edge
- `xr_te`: x-distance from the right side of the bound vortex to the trailing edge
- `rc`: (optional) position of the center of the bound vortex, defaults to `(rl+rr)/2`
- `xc_te`: (optional) x-distance from the center of the bound vortex to the
    trailing edge, defaults to `(xl_te+xr_te)/2`
"""
Horseshoe(rl, rr, rcp, ncp, xl_te, xr_te; rc=(rl+rr)/2, xc_te=(xl_te+xr_te)/2) = Horseshoe(rl, rc, rr, rcp, ncp, xl_te, xc_te, xr_te)

function Horseshoe(rl, rc, rr, rcp, ncp, xl_te, xc_te, xr_te, core_size)
    TF = promote_type(eltype(rl), eltype(rc), eltype(rr), eltype(rcp), eltype(ncp), typeof(xl_te), typeof(xc_te), typeof(xr_te), typeof(core_size))
    return Horseshoe{TF}(rl, rc, rr, rcp, ncp, xl_te, xc_te, xr_te, core_size)
end

@inline Base.eltype(::Type{Horseshoe{TF}}) where TF = TF
@inline Base.eltype(::Horseshoe{TF}) where TF = TF

# -------------------------------

"""
    Ring{TF}

Vortex ring shaped panel element.

**Fields**
 - `rtl`: position of the left side of the top bound vortex
 - `rtc`: position of the center of the top bound vortex
 - `rtr`: position of the right side of the top bound vortex
 - `rbl`: position of the left side of the bottom bound vortex
 - `rbc`: position of the center of the bottom bound vortex
 - `rbr`: position of the right side of the bottom bound vortex
 - `rcp`: position of the panel control point
 - `ncp`: normal vector at the panel control point
 - `core`: vortex core size
"""
struct Ring{TF} <: AbstractPanel{TF}
    rtl::SVector{3, TF}
    rtc::SVector{3, TF}
    rtr::SVector{3, TF}
    rbl::SVector{3, TF}
    rbc::SVector{3, TF}
    rbr::SVector{3, TF}
    rcp::SVector{3, TF}
    ncp::SVector{3, TF}
    core_size::TF
end

"""
    Ring(rtl, rtr, rbl, rbr, rcp, ncp, core_size; rtc, rbc)

Construct and return a vortex ring panel.

**Arguments**
 - `rtl`: position of the left side of the top bound vortex
 - `rtr`: position of the right side of the top bound vortex
 - `rbl`: position of the left side of the bottom bound vortex
 - `rbr`: position of the right side of the bottom bound vortex
 - `rcp`: position of the panel control point
 - `ncp`: normal vector at the panel control point
 - `core_size`: vortex core size
 - `rtc`: (optional) position of the center of the top bound vortex, defaults to `(rtl+rtr)/2`
 - `rbc`: (optional) position of the center of the bottom bound vortex, defaults to `(rbl+rbr)/2`
"""
Ring(rtl, rtr, rbl, rbr, rcp, ncp, core_size; rtc=(rtl+rtr)/2, rbc=(rbl+rbr)/2) =
    Ring(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)

function Ring(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)
    TF = promote_type(eltype(rtl), eltype(rtc), eltype(rtr), eltype(rbl), eltype(rbc), eltype(rbr), eltype(rcp), eltype(ncp), typeof(core_size))
    return Ring{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)
end

@inline Base.eltype(::Type{Ring{TF}}) where TF = TF
@inline Base.eltype(::Ring{TF}) where TF = TF

# --- geometric functions --- #

"""
    translate(panel::AbstractPanel, r)

Return a copy of `panel` translated the distance specified by vector `r`
"""
translate

@inline function translate(panel::Horseshoe, r)

    rl = panel.rl + r
    rc = panel.rc + r
    rr = panel.rr + r
    rcp = panel.rcp + r
    ncp = panel.ncp
    xl_te = panel.xl_te
    xc_te = panel.xc_te
    xr_te = panel.xr_te
    core_size = panel.core_size

    return Horseshoe(rl, rc, rr, rcp, ncp, xl_te, xc_te, xr_te, core_size)
end

@inline function translate(panel::Ring, r)

    rtl = panel.rtl + r
    rtc = panel.rtc + r
    rtr = panel.rtr + r
    rbl = panel.rbl + r
    rbc = panel.rbc + r
    rbr = panel.rbr + r
    rcp = panel.rcp + r
    ncp = panel.ncp
    core_size = panel.core_size

    return Ring(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)
end

"""
    translate(panels::AbstractVector{<:AbstractPanel}, r)

Return a copy of `panels` translated the distance specified by vector `r`
"""
translate(panels::AbstractVector{<:AbstractPanel}, r) = translate.(panels, Ref(r))

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

"""
    reflect(panel::AbstractPanel)

Reflects a panel across the y-axis.
"""
reflect

function reflect(panel::Horseshoe)

    rl = flipy(panel.rr)
    rc = flipy(panel.rc)
    rr = flipy(panel.rl)
    rcp = flipy(panel.rcp)
    ncp = flipy(panel.ncp)
    xl_te = panel.xr_te
    xc_te = panel.xc_te
    xr_te = panel.xl_te
    core_size = panel.core_size

    return Horseshoe(rl, rc, rr, rcp, ncp, xl_te, xc_te, xr_te, core_size)
end

function reflect(panel::Ring)

    rtl = flipy(panel.rtr)
    rtc = flipy(panel.rtc)
    rtr = flipy(panel.rtl)
    rbl = flipy(panel.rbr)
    rbc = flipy(panel.rbc)
    rbr = flipy(panel.rbl)
    rcp = flipy(panel.rcp)
    ncp = flipy(panel.ncp)
    core_size = panel.core_size

    return Ring(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)
end

"""
    reflect(panels::AbstractMatrix{<:AbstractPanel})

Reflects panels about the y-axis, preserving panel grid-based ordering
"""
reflect(panels::AbstractMatrix{<:AbstractPanel}) = reflect.(reverse(panels, dims=2))

"""
    controlpoint(panel::AbstractPanel)

Return the control point of `panel` (typically located at the 3/4 chord)
"""
controlpoint

@inline controlpoint(panel::Horseshoe) = panel.rcp

@inline controlpoint(panel::Ring) = panel.rcp

"""
    normal(panel::AbstractPanel)

Return the normal vector of `panel`
"""
normal

@inline normal(panel::Horseshoe) = panel.ncp

@inline normal(panel::Ring) = panel.ncp

"""
    midpoint(panel::AbstractPanel)

Return the midpoint of the bound vortex on `panel`
"""
midpoint

@inline midpoint(panel::Horseshoe) = panel.rc

@inline midpoint(panel::Ring) = panel.rtc

"""
    left_midpoint(panel::AbstractPanel)

Return the midpoint of the left bound vortex on `panel`
"""
left_midpoint

@inline left_midpoint(panel::Horseshoe) = panel.rl + SVector(panel.xl_te/2, 0, 0)

@inline left_midpoint(panel::Ring) = (panel.rtl + panel.rbl)/2

"""
    right_midpoint(panel::AbstractPanel)

Return the midpoint of the right bound vortex on `panel`
"""
right_midpoint

@inline right_midpoint(panel::Horseshoe) = panel.rr + SVector(panel.xr_te/2, 0, 0)
@inline right_midpoint(panel::Ring) = (panel.rtr + panel.rbr)/2

"""
    top_vector(panel)

Returns the path of the top bound vortex for `panel`
"""
top_vector
@inline top_vector(panel::Horseshoe) = panel.rr - panel.rl
@inline top_vector(panel::Ring) = panel.rtr - panel.rtl

"""
    left_vector(panel)

Return the path of the left bound vortex for `panel`
"""
left_vector
@inline left_vector(panel::Horseshoe) = -SVector(panel.xl_te, 0, 0)
@inline left_vector(panel::Ring) = panel.rtl - panel.rbl

"""
    right_vector(panel)

Return the path of the right bound vortex for `panel`
"""
right_vector
@inline right_vector(panel::Horseshoe) = SVector(panel.xr_te, 0, 0)
@inline right_vector(panel::Ring) = panel.rbr - panel.rtr

"""
    endpoints(panel)

Return the endpoints of the bound vortices in `panel`.
"""
endpoints
@inline function endpoints(panel::Horseshoe)
    r11 = panel.rl
    r12 = panel.rr
    r21 = panel.rl + SVector(panel.xl_te, 0, 0)
    r22 = panel.rr + SVector(panel.xr_te, 0, 0)
    return r11, r12, r21, r22
end
@inline endpoints(panel::Ring) = panel.rtl, panel.rtr, panel.rbl, panel.rbr

"""
    reflected_endpoints(panel)

Return the endpoints of the bound vortice in `panel`, reflected across the y-axis.
"""
reflected_endpoints
@inline function reflected_endpoints(panel::Horseshoe)
    rtl = flipy(panel.rr)
    rtr = flipy(panel.rl)
    rbl = flipy(panel.rr) + SVector(panel.xr_te, 0, 0)
    rbr = flipy(panel.rl) + SVector(panel.xl_te, 0, 0)
    return rtl, rtr, rbl, rbr
end
@inline function reflected_endpoints(panel::Ring)
    rtl = flipy(panel.rtr)
    rtr = flipy(panel.rtl)
    rbl = flipy(panel.rbr)
    rbr = flipy(panel.rbl)
    return rtl, rtr, rbl, rbr
end

"""
    trefftz_endpoints(panel)

Return the two points where trailing vortices originating from `panel` are shed
into the farfield.
"""
@inline function trefftz_endpoints(panel::Horseshoe)
    rl = panel.rl + SVector(panel.xl_te, 0, 0)
    rr = panel.rr + SVector(panel.xr_te, 0, 0)
    return rl, rr
end

@inline trefftz_endpoints(panel::Ring) = panel.rbl, panel.rbr

"""
    trefftz_center(panel)

Return the center of the panel in the farfield created by vortices shed by `panel`
"""
trefftz_center
@inline trefftz_center(panel::Horseshoe) = panel.rc + SVector(panel.xc_te, 0, 0)
@inline trefftz_center(panel::Ring) = panel.rbc

# --- induced velocity functions --- #

"""
    induced_velocity(rcp, panel, xhat, symmetric, same_surface, trailing, args...)

Computes the normalized induced velocity at control point `rcp` from `panel`.

# Arguments
 - `rcp`: Control point for receiving panel
 - `panel`: Sending panel
 - `symmetric`: Flag indicating whether sending surface is symmetric
 - `same_surface`: Flag indicating whether the sending and receiving surfaces are the same
 - `trailing`: Flag indicating whether sending panel sheds trailing vortices
 - `xhat`: Direction in which to shed trailing vortices (if applicable)

# Optional Arguments for Horseshoe Vortices
 - `include_bound`: Indicates whether to include the bound vortex in induced velocity
    calculations, defaults to true

# Optional Arguments for Vortex Rings
 - `include_top`: Indicates whether to include the top bound vortex in induced velocity
    calculations, defaults to true
 - `include_bottom`: Indicates whether to include the bottom bound vortex in induced velocity
    calculations, defaults to true
"""
induced_velocity

# horseshoe vortex
@inline function induced_velocity(rcp, panel::Horseshoe, xhat, symmetric,
    same_surface, include_bound=true)

    rtl, rtr, rbl, rbr = endpoints(panel)
    r11, r12, r21, r22 = rcp - rtl, rcp - rtr, rcp - rbl, rcp - rbr
    core_size = panel.core_size
    Vhat = trailing_induced_velocity(r21, r22, xhat, core_size, same_surface) # trailing vortices
    Vhat += bound_induced_velocity(r21, r11, core_size, same_surface) # left side (to trailing edge)
    Vhat += bound_induced_velocity(r12, r22, core_size, same_surface) # right side (to trailing edge)
    if include_bound
        Vhat += bound_induced_velocity(r11, r12, core_size, same_surface) # bound vortex
    end

    # add contribution from other side (if applicable)
    if symmetric && not_on_symmetry_plane(panel.rl, panel.rr)
        rtl, rtr, rbl, rbr = reflected_endpoints(panel)
        r11, r12, r21, r22 = rcp - rtl, rcp - rtr, rcp - rbl, rcp - rbr
        Vhat += trailing_induced_velocity(r21, r22, xhat, core_size, same_surface) # trailing vortices
        Vhat += bound_induced_velocity(r21, r11, core_size, same_surface) # left side (to trailing edge)
        Vhat += bound_induced_velocity(r12, r22, core_size, same_surface) # right side (to trailing edge)
        Vhat += bound_induced_velocity(r11, r12, core_size, same_surface) # bound vortex
    end

    return Vhat
end

# vortex ring
@inline function induced_velocity(rcp, panel::Ring, xhat, symmetric, same_surface,
    trailing, include_top=true, include_bottom=true)

    rtl, rtr, rbl, rbr = endpoints(panel)
    r11, r12, r21, r22 = rcp - rtl, rcp - rtr, rcp - rbl, rcp - rbr
    core_size = panel.core_size
    Vhat = bound_induced_velocity(r21, r11, core_size, same_surface) # left bound vortex
    Vhat += bound_induced_velocity(r12, r22, core_size, same_surface) # right bound vortex

    if include_top
        Vhat += bound_induced_velocity(r11, r12, core_size, same_surface) # top bound vortex
    end

    if trailing
        Vhat += trailing_induced_velocity(r21, r22, xhat, core_size, same_surface) # trailing vortices
    else
        if include_bottom
            Vhat += bound_induced_velocity(r22, r21, core_size, same_surface) # bottom bound vortex
        end
    end

    # add contribution from other side (if applicable)
    if symmetric && not_on_symmetry_plane(panel.rtl, panel.rtr)
        rtl, rtr, rbl, rbr = reflected_endpoints(panel)
        r11, r12, r21, r22 = rcp - rtl, rcp - rtr, rcp - rbl, rcp - rbr
        Vhat += bound_induced_velocity(r21, r11, core_size, same_surface) # left bound vortex
        Vhat += bound_induced_velocity(r12, r22, core_size, same_surface) # right bound vortex
        Vhat += bound_induced_velocity(r11, r12, core_size, same_surface) # top bound vortex
        if trailing
            Vhat += trailing_induced_velocity(r21, r22, xhat, core_size, same_surface) # trailing vortices
        else
            Vhat += bound_induced_velocity(r22, r21, core_size, same_surface) # bottom bound vortex
        end
    end

    return Vhat
end


# --- helper functions --- #

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
    trailing_induced_velocity(r1, r2, xhat, core_size, same_surface)

Compute the induced velocity (per unit circulation) for two vortices trailing in
the `xhat` direction, at a control point located at `r1` relative to the start of the
left trailing vortex and `r2` relative to the start of the right trailing vortex.
"""
@inline function trailing_induced_velocity(r1, r2, xhat, core_size, same_surface)

    nr1 = norm(r1)
    nr2 = norm(r2)

    if same_surface || core_size == 0
        # no finite core
        f1 = cross(r1, xhat)/(nr1 - dot(r1, xhat))/nr1
        f2 = cross(r2, xhat)/(nr2 - dot(r2, xhat))/nr2
    else
        # finite core
        f1 = cross(r1, xhat)/(nr1^2 - dot(r1, xhat)^2 + core_size^2)
        f2 = cross(r2, xhat)/(nr2^2 - dot(r2, xhat)^2 + core_size^2)
    end

    Vhat = (f1 - f2)/(4*pi)

    return Vhat
end

"""
    bound_induced_velocity(r1, r2, core_size, same_surface)

Compute the induced velocity (per unit circulation) for a bound vortex, at a
control point located at `r1` relative to the start of the bound vortex and `r2`
relative to the end of the bound vortex
"""
@inline function bound_induced_velocity(r1, r2, core_size, same_surface)

    nr1 = norm(r1)
    nr2 = norm(r2)

    if same_surface || core_size == 0
        # no finite core
        f1 = cross(r1, r2)/(nr1*nr2 + dot(r1, r2))
        f2 = (1/nr1 + 1/nr2)
    else
        # finite core
        rdot = dot(r1, r2)
        nr1s, nr2s = nr1^2, nr2^2
        f1 = cross(r1, r2)/(nr1s*nr2s - rdot^2 + (nr1s + nr2s - 2*nr1*nr2)*core_size^2)
        f2 = (nr1s - rdot/sqrt(nr1s + core_size^2)) + (nr2s - rdot/sqrt(nr2s + core_size^2))
    end

    Vhat = (f1*f2)/(4*pi)

    return Vhat
end

"""
    vortex_induced_drag(rj, Γj, ri, Γi, nhati)

Return induced drag from vortex `j` induced on panel `i`
"""
@inline function vortex_induced_drag(rj, Γj, ri, Γi, nhati, core_size, same_surface)

    rij = ri - rj
    ε = ifelse(same_surface, zero(core_size), core_size)
    Vthetai = SVector(0, -Γj*rij[3], Γj*rij[2]) / (2*pi*sqrt(rij[2]^2 + rij[3]^2 + core_size^2))
    Vn = -dot(Vthetai, nhati)

    Di = RHO/2.0*Γi*Vn

    return Di
end
