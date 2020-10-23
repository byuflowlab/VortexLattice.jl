"""
    AbstractPanel{TF}

Supertype of vortex lattice method panel types.
"""
abstract type AbstractPanel{TF} end

"""
    top_left(panel::AbstractPanel)

Return the top left vertex of `panel`
"""
top_left

"""
    top_center(panel::AbstractPanel)

Return the top center vertex of `panel`
"""
top_center

"""
    top_right(panel::AbstractPanel)

Return the top right vertex of `panel`
"""
top_right

"""
    bottom_left(panel::AbstractPanel)

Return the bottom left vertex of `panel`
"""
bottom_left

"""
    bottom_center(panel::AbstractPanel)

Return the bottom center vertex of `panel`
"""
bottom_center

"""
    bottom_right(panel::AbstractPanel)

Return the bottom right vertex of `panel`
"""
bottom_right

"""
    controlpoint(panel::AbstractPanel)

Return the control point of `panel` (typically located at the 3/4 chord)
"""
controlpoint

"""
    normal(panel::AbstractPanel)

Return the normal vector of `panel`
"""
normal

"""
    get_core_size(panel::AbstractPanel)

Return the panel core size
"""
get_core_size

"""
    translate(panel::AbstractPanel, r)

Return a copy of `panel` translated the distance specified by vector `r`
"""
translate(panel::AbstractPanel, r)

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
reflect(panel::AbstractPanel)

"""
    reflect(panels::AbstractMatrix{<:AbstractPanel})

Reflects panels about the y-axis, preserving panel grid-based ordering
"""
reflect(panels::AbstractMatrix{<:AbstractPanel}) = reflect.(reverse(panels, dims=2))

"""
    induced_velocity(rcp, panel, same_id, symmetric, xhat, args...)

Computes the normalized induced velocity at control point `rcp` from `panel`.

# Arguments
 - `rcp`: Control point for receiving panel
 - `panel`: Sending panel
 - `same_id`: Flag indicating whether the sending and receiving surfaces have the same surface ID
 - `symmetric`: Flag indicating whether sending panel should be mirrored across the y-axis
 - `xhat`: Direction in which to shed trailing vortices

# Optional Arguments for Horseshoe Vortices
 - `include_bound`: Indicates whether to include the bound vortex in induced velocity
    calculations, defaults to true

# Optional Arguments for Vortex Rings
 - `trailing`: Flag indicating whether sending panel sheds trailing vortices
 - `include_top`: Indicates whether to include the top bound vortex in induced velocity
    calculations, defaults to true
 - `include_bottom`: Indicates whether to include the bottom bound vortex in induced velocity
    calculations, defaults to true
"""
induced_velocity

"""
    panel_induced_velocity(receiving, Ir, sending, Is, same_surface, same_id, symmetric, xhat, trailing)

Computes the normalized induced velocity on panel `receiving` from panel `sending`.

# Arguments
 - `receiving`: Receiving panel
 - `Ir`: Indices of `receiving` panel
 - `sending`: Sending panel
 - `Is`: Indices of `sending` panel
 - `same_surface`: Flag indicating whether the sending and receiving surfaces are the same
 - `same_id`: Flag indicating whether the sending and receiving surfaces have the same surface ID
 - `symmetric`: Flag indicating whether sending panel should be mirrored across the y-axis
 - `xhat`: Direction in which to shed trailing vortices
 - `trailing`: Indicates whether the panel is on the trailing edge
"""
panel_induced_velocity

"""
    panel_circulation(panel, Γ1, Γ2)

Return the circulation on `panel` given the circulation strength of the previous
bound vortex Γ1 and the current bound vortex Γ2
"""
panel_circulation

"""
    influence_coefficients!(AIC, receiving::AbstractMatrix{<:AbstractPanel},
        sending::AbstractMatrix{<:AbstractPanel}, same_id, symmetric, xhat)

Compute the AIC coefficients corresponding to the influence of the panels in
`sending` on the panels in `receiving`.
"""
influence_coefficients!(AIC, receiving::AbstractMatrix{<:AbstractPanel},
    sending::AbstractMatrix{<:AbstractPanel}, same_id, symmetric, xhat)

"""
    left_center(panel::AbstractPanel)

Return the center of the left bound vortex on `panel`
"""
@inline left_center(panel::AbstractPanel) = (top_left(panel) + bottom_left(panel))/2

"""
    right_center(panel::AbstractPanel)

Return the center of the right bound vortex on `panel`
"""
@inline right_center(panel::AbstractPanel) = (top_right(panel) + bottom_right(panel))/2

"""
    top_vector(panel::AbstractPanel)

Returns the path of the top bound vortex for `panel`
"""
@inline top_vector(panel::AbstractPanel) = top_right(panel) - top_left(panel)

"""
    left_vector(panel)

Return the path of the left bound vortex for `panel`
"""
@inline left_vector(panel::AbstractPanel) = top_left(panel) - bottom_left(panel)

"""
    right_vector(panel)

Return the path of the right bound vortex for `panel`
"""
@inline right_vector(panel::AbstractPanel) = bottom_right(panel) - top_right(panel)

"""
    bottom_vector(panel)

Return the path of the bottom bound vortex for `panel`
"""
@inline bottom_vector(panel::AbstractPanel) = bottom_left(panel) - bottom_right(panel)

# --- Horseshoe Vortex Implementation --- #

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
 - `chord`: section chord length for finite core calculations
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
    chord::TF
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
- `chord`: section chord length for finite core calculations
- `rc`: (optional) position of the center of the bound vortex, defaults to `(rl+rr)/2`
- `xc_te`: (optional) x-distance from the center of the bound vortex to the
    trailing edge, defaults to `(xl_te+xr_te)/2`
"""
Horseshoe(rl, rr, rcp, ncp, xl_te, xr_te, chord; rc=(rl+rr)/2, xc_te=(xl_te+xr_te)/2) = Horseshoe(rl, rc, rr, rcp, ncp, xl_te, xc_te, xr_te, chord)

function Horseshoe(rl, rc, rr, rcp, ncp, xl_te, xc_te, xr_te, chord)
    TF = promote_type(eltype(rl), eltype(rc), eltype(rr), eltype(rcp), eltype(ncp), typeof(xl_te), typeof(xc_te), typeof(xr_te), typeof(chord))
    return Horseshoe{TF}(rl, rc, rr, rcp, ncp, xl_te, xc_te, xr_te, chord)
end

@inline Base.eltype(::Type{Horseshoe{TF}}) where TF = TF
@inline Base.eltype(::Horseshoe{TF}) where TF = TF

@inline top_left(panel::Horseshoe) = panel.rl

@inline top_center(panel::Horseshoe) = panel.rc

@inline top_right(panel::Horseshoe) = panel.rr

@inline bottom_left(panel::Horseshoe) = panel.rl + SVector(panel.xl_te, 0, 0)

@inline bottom_center(panel::Horseshoe) = panel.rc + SVector(panel.xc_te, 0, 0)

@inline bottom_right(panel::Horseshoe) = panel.rr +  SVector(panel.xr_te, 0, 0)

@inline controlpoint(panel::Horseshoe) = panel.rcp

@inline normal(panel::Horseshoe) = panel.ncp

@inline get_core_size(panel::Horseshoe) = panel.chord

@inline function induced_velocity(rcp, panel::Horseshoe, symmetric, same_id,
    xhat, include_bound)

    # get distance to control point from each vortex filament vertex
    rtl = top_left(panel)
    rtr = top_right(panel)
    r11 = rcp - rtl
    r12 = rcp - rtr
    r21 = rcp - bottom_left(panel)
    r22 = rcp - bottom_right(panel)

    core_size = get_core_size(panel)

    # add contribution from trailing vortices
    Vhat = trailing_induced_velocity(r21, r22, xhat, same_id, core_size)

    # add contribution from left bound vortex
    Vhat += bound_induced_velocity(r21, r11, same_id, core_size)

    # add contribution from right bound vortex
    Vhat += bound_induced_velocity(r12, r22, same_id, core_size)

    # add contribution from top bound vortex (unless otherwise specified)
    if include_bound
        Vhat += bound_induced_velocity(r11, r12, same_id, core_size)
    end

    # add contribution from other side (if applicable)
    if symmetric && not_on_symmetry_plane(rtl, rtr)
        Vhat += induced_velocity(rcp, reflect(panel), false, same_id, xhat, true)
    end

    return Vhat
end

@inline function panel_induced_velocity(receiving, Ir, sending::Horseshoe, Is,
    same_surface, symmetric, same_id, xhat, trailing)

    rc = top_center(receiving)

    include_bound = !(same_surface && Ir == Is)

    return induced_velocity(rc, sending, symmetric, same_id, xhat, include_bound)
end

@inline panel_circulation(panel::Horseshoe, Γ1, Γ2) = Γ2

@inline function influence_coefficients!(AIC, receiving::AbstractMatrix{<:AbstractPanel},
    sending::AbstractMatrix{<:Horseshoe}, symmetric, same_id, xhat)

    Nr = length(receiving)
    Ns = length(sending)

    cr = CartesianIndices(receiving)
    cs = CartesianIndices(sending)
    nchordwise = size(sending,1)

    # always include bound vortices
    include_bound = true

    # loop over receiving panels
    for i = 1:Nr
        I = cr[i]

        # control point location
        rcp = controlpoint(receiving[I])

        # normal vector body axis
        nhat = normal(receiving[I])

        # loop over sending panels
        for j = 1:Ns
            J = cs[j]

            #TODO: take advantage of shared edges when calculating influence coefficients
            # (this should speed up the calculations)

            Vhat = induced_velocity(rcp, sending[J], symmetric, same_id,
                xhat, include_bound)

            AIC[i, j] = dot(Vhat, nhat)
        end
    end

    return AIC
end

@inline function translate(panel::Horseshoe, r)

    rl = panel.rl + r
    rc = panel.rc + r
    rr = panel.rr + r
    rcp = panel.rcp + r
    ncp = panel.ncp
    xl_te = panel.xl_te
    xc_te = panel.xc_te
    xr_te = panel.xr_te
    chord = panel.chord

    return Horseshoe(rl, rc, rr, rcp, ncp, xl_te, xc_te, xr_te, chord)
end

@inline function reflect(panel::Horseshoe)

    rl = flipy(panel.rr)
    rc = flipy(panel.rc)
    rr = flipy(panel.rl)
    rcp = flipy(panel.rcp)
    ncp = flipy(panel.ncp)
    xl_te = panel.xr_te
    xc_te = panel.xc_te
    xr_te = panel.xl_te
    chord = panel.chord

    return Horseshoe(rl, rc, rr, rcp, ncp, xl_te, xc_te, xr_te, chord)
end

# --- Vortex Ring Implementation --- #

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
 - `chord`: chord length for finite core calculations
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
    chord::TF
end

"""
    Ring(rtl, rtr, rbl, rbr, rcp, ncp, chord; rtc, rbc)

Construct and return a vortex ring panel.

**Arguments**
 - `rtl`: position of the left side of the top bound vortex
 - `rtr`: position of the right side of the top bound vortex
 - `rbl`: position of the left side of the bottom bound vortex
 - `rbr`: position of the right side of the bottom bound vortex
 - `rcp`: position of the panel control point
 - `ncp`: normal vector at the panel control point
 - `chord`: chord length for finite core calculations
 - `rtc`: (optional) position of the center of the top bound vortex, defaults to `(rtl+rtr)/2`
 - `rbc`: (optional) position of the center of the bottom bound vortex, defaults to `(rbl+rbr)/2`
"""
Ring(rtl, rtr, rbl, rbr, rcp, ncp, chord; rtc=(rtl+rtr)/2, rbc=(rbl+rbr)/2) =
    Ring(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, chord)

function Ring(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, chord)
    TF = promote_type(eltype(rtl), eltype(rtc), eltype(rtr), eltype(rbl), eltype(rbc), eltype(rbr), eltype(rcp), eltype(ncp), typeof(chord))
    return Ring{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, chord)
end

@inline Base.eltype(::Type{Ring{TF}}) where TF = TF
@inline Base.eltype(::Ring{TF}) where TF = TF

@inline controlpoint(panel::Ring) = panel.rcp

@inline normal(panel::Ring) = panel.ncp

@inline midpoint(panel::Ring) = panel.rtc

@inline top_left(panel::Ring) = panel.rtl

@inline top_center(panel::Ring) = panel.rtc

@inline top_right(panel::Ring) = panel.rtr

@inline bottom_left(panel::Ring) = panel.rbl

@inline bottom_center(panel::Ring) = panel.rbc

@inline bottom_right(panel::Ring) = panel.rbr

@inline get_core_size(panel::Ring) = panel.chord

@inline function induced_velocity(rcp, panel::Ring, symmetric, same_id,
    xhat, trailing, include_top, include_bottom)

    # get distance to control point from each vortex filament vertex
    rtl = top_left(panel)
    rtr = top_right(panel)
    r11 = rcp - rtl
    r12 = rcp - rtr
    r21 = rcp - bottom_left(panel)
    r22 = rcp - bottom_right(panel)

    core_size = get_core_size(panel)

    # add contribution from left bound vortex
    Vhat = bound_induced_velocity(r21, r11, same_id, core_size)

    # add contribution from right bound vortex
    Vhat += bound_induced_velocity(r12, r22, same_id, core_size)

    # add contribution from top bound vortex (unless otherwise specified)
    if include_top
        Vhat += bound_induced_velocity(r11, r12, same_id, core_size)
    end

    if trailing
        # add contribution from trailing vortices
        Vhat += trailing_induced_velocity(r21, r22, xhat, same_id, core_size)
    else
        # add contribution from bottom bound vortex (unless otherwise specified)
        if include_bottom
            Vhat += bound_induced_velocity(r22, r21, same_id, core_size)
        end
    end

    # add contribution from other side (if applicable)
    if symmetric && not_on_symmetry_plane(rtl, rtr)
        Vhat += induced_velocity(rcp, reflect(panel), false, same_id, xhat,
            trailing, true, true)
    end

    return Vhat
end

@inline function panel_induced_velocity(receiving, Ir, sending::Ring, Is,
    same_surface, symmetric, same_id, xhat, trailing)

    rc = top_center(receiving)

    include_top = !(same_surface && Ir == Is)
    include_bottom = !(same_surface && (Ir[1] == Is[1]+1 && Ir[2] == Is[2]))

    return induced_velocity(rc, sending, symmetric, same_id, xhat, trailing,
        include_top, include_bottom)
end

@inline panel_circulation(panel::Ring, Γ1, Γ2) = Γ2 - Γ1

@inline function influence_coefficients!(AIC, receiving::AbstractMatrix{<:AbstractPanel},
    sending::AbstractMatrix{<:Ring}, symmetric, same_id, xhat)

    Nr = length(receiving)
    Ns = length(sending)

    cr = CartesianIndices(receiving)
    cs = CartesianIndices(sending)
    nchordwise = size(sending,1)

    # loop over receiving panels
    for i = 1:Nr
        I = cr[i]

        # control point location
        rcp = controlpoint(receiving[I])

        # normal vector body axis
        nhat = normal(receiving[I])

        # loop over sending panels
        for j = 1:Ns
            J = cs[j]

            trailing = J[1] == nchordwise
            include_top = true
            include_bottom = true
            Vhat = induced_velocity(rcp, sending[J], symmetric, same_id,
                xhat, trailing, include_top, include_bottom)

            AIC[i, j] = dot(Vhat, nhat)
        end
    end

    return AIC
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
    chord = panel.chord

    return Ring(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, chord)
end

@inline function reflect(panel::Ring)

    rtl = flipy(panel.rtr)
    rtc = flipy(panel.rtc)
    rtr = flipy(panel.rtl)
    rbl = flipy(panel.rbr)
    rbc = flipy(panel.rbc)
    rbr = flipy(panel.rbl)
    rcp = flipy(panel.rcp)
    ncp = flipy(panel.ncp)
    chord = panel.chord

    return Ring(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, chord)
end

# --- Struct to store panel properties --- #

"""
    PanelProperties

Panel specific properties calculated during the vortex lattice method analysis.

**Fields**
 - `gamma`: Panel circulation strength (normalized by the freestream velocity)
 - `v`: Local velocity at the panel's center (typically the quarter-chord), normalized
 - `cf`: Bound vortex force per unit length, normalized by `QINF*S` where `QINF`
    is the dynamic pressure and `S` is the user-provided reference area.
 - `cfl`: Left vortex force per unit length, normalized by `QINF*S`
 - `cfr`: Right vortex force per unit length, normalized by `QINF*S`
"""
struct PanelProperties{TF}
    gamma::TF
    v::SVector{3, TF}
    cf::SVector{3, TF}
    cfl::SVector{3, TF}
    cfr::SVector{3, TF}
end

function PanelProperties(gamma, v, cf, cfl, cfr)

    TF = promote_type(typeof(gamma), typeof(v), eltype(cf), eltype(cfl), eltype(cfr))

    return PanelProperties{TF}(gamma, v, cf, cfl, cfr)
end

Base.eltype(::Type{PanelProperties{TF}}) where TF = TF
Base.eltype(::PanelProperties{TF}) where TF = TF


# --- internal functions --- #

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
    trailing_induced_velocity(r1, r2, xhat, same_id, core_size)

Compute the induced velocity (per unit circulation) for two vortices trailing in
the `xhat` direction, at a control point located at `r1` relative to the start of the
left trailing vortex and `r2` relative to the start of the right trailing vortex.
"""
@inline function trailing_induced_velocity(r1, r2, xhat, same_id, core_size)

    nr1 = norm(r1)
    nr2 = norm(r2)

    r1dot = dot(r1, xhat)
    r2dot = dot(r2, xhat)

    if same_id || iszero(core_size)
        # no finite core
        f1 = cross(r1, xhat)/(nr1 - r1dot)/nr1
        f2 = cross(r2, xhat)/(nr2 - r2dot)/nr2
    else
        # finite core
        εs = core_size^2
        tmp1 = εs/(nr1 + r1dot)
        tmp2 = εs/(nr2 + r2dot)
        f1 = cross(r1, xhat)/(nr1 - r1dot + tmp1)/nr1
        f2 = cross(r2, xhat)/(nr2 - r2dot + tmp2)/nr2
    end

    Vhat = (f1 - f2)/(4*pi)

    return Vhat
end

"""
    bound_induced_velocity(r1, r2, same_id, core_size)

Compute the induced velocity (per unit circulation) for a bound vortex, at a
control point located at `r1` relative to the start of the bound vortex and `r2`
relative to the end of the bound vortex
"""
@inline function bound_induced_velocity(r1, r2, same_id, core_size)

    nr1 = norm(r1)
    nr2 = norm(r2)

    if same_id || iszero(core_size)
        # no finite core
        f1 = cross(r1, r2)/(nr1*nr2 + dot(r1, r2))
        f2 = (1/nr1 + 1/nr2)
    else
        # finite core
        rdot = dot(r1, r2)
        r1s, r2s, εs = nr1^2, nr2^2, core_size^2
        f1 = cross(r1, r2)/(r1s*r2s - rdot^2 + εs*(r1s + r2s - 2*nr1*nr2))
        f2 = (r1s - rdot)/sqrt(r1s + εs) + (r2s - rdot)/sqrt(r2s + εs)
    end

    Vhat = (f1*f2)/(4*pi)

    return Vhat
end
