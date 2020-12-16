"""
    VortexRing{TF}

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
 - `core_size`: finite core size
"""
struct VortexRing{TF}
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
    VortexRing(rtl, rtr, rbl, rbr, rcp, ncp, core_size; rtc, rbc)

Construct and return a vortex ring panel.

**Arguments**
 - `rtl`: position of the left side of the top bound vortex
 - `rtr`: position of the right side of the top bound vortex
 - `rbl`: position of the left side of the bottom bound vortex
 - `rbr`: position of the right side of the bottom bound vortex
 - `rcp`: position of the panel control point
 - `ncp`: normal vector at the panel control point
 - `core_size`: core_size length for finite core calculations
 - `rtc`: (optional) position of the center of the top bound vortex, defaults to `(rtl+rtr)/2`
 - `rbc`: (optional) position of the center of the bottom bound vortex, defaults to `(rbl+rbr)/2`
"""
VortexRing(args...; kwargs...)

VortexRing(rtl, rtr, rbl, rbr, rcp, ncp, core_size; rtc=(rtl+rtr)/2, rbc=(rbl+rbr)/2) =
    VortexRing(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)

function VortexRing(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)
    TF = promote_type(eltype(rtl), eltype(rtc), eltype(rtr), eltype(rbl), eltype(rbc), eltype(rbr), eltype(rcp), eltype(ncp), typeof(core_size))
    return VortexRing{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)
end

@inline Base.eltype(::Type{VortexRing{TF}}) where TF = TF
@inline Base.eltype(::VortexRing{TF}) where TF = TF

"""
    top_left(panel::VortexRing)

Return the top left vertex of `panel`
"""
@inline top_left(panel::VortexRing) = panel.rtl

"""
    top_center(panel::VortexRing)

Return the top center vertex of `panel`
"""
@inline top_center(panel::VortexRing) = panel.rtc

"""
    top_right(panel::VortexRing)

Return the top right vertex of `panel`
"""
@inline top_right(panel::VortexRing) = panel.rtr

"""
    bottom_left(panel::VortexRing)

Return the bottom left vertex of `panel`
"""
@inline bottom_left(panel::VortexRing) = panel.rbl

"""
    bottom_center(panel::VortexRing)

Return the bottom center vertex of `panel`
"""
@inline bottom_center(panel::VortexRing) = panel.rbc

"""
    bottom_right(panel::VortexRing)

Return the bottom right vertex of `panel`
"""
@inline bottom_right(panel::VortexRing) = panel.rbr

"""
    controlpoint(panel::VortexRing)

Return the control point of `panel` (typically located at the 3/4 chord)
"""
@inline controlpoint(panel::VortexRing) = panel.rcp

"""
    normal(panel::VortexRing)

Return the normal vector of `panel`
"""
@inline normal(panel::VortexRing) = panel.ncp

"""
    get_core_size(panel::VortexRing)

Return the panel core size
"""
@inline get_core_size(panel::VortexRing) = panel.core_size

"""
    translate(panel::VortexRing, r)

Return a copy of `panel` translated the distance specified by vector `r`
"""
@inline function translate(panel::VortexRing, r)

    rtl = panel.rtl + r
    rtc = panel.rtc + r
    rtr = panel.rtr + r
    rbl = panel.rbl + r
    rbc = panel.rbc + r
    rbr = panel.rbr + r
    rcp = panel.rcp + r
    ncp = panel.ncp
    core_size = panel.core_size

    return VortexRing(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)
end

"""
    translate(surface, r)

Return a copy of `surface` translated the distance specified by vector `r`
"""
translate(surface::AbstractMatrix, r) = translate.(surface, Ref(r))

"""
    translate!(surface, r)

Translate the panels in `surface` the distance specified by vector `r`
"""
function translate!(surface, r)

    for i in eachindex(surface)
        surface[i] = translate(surface[i], r)
    end

    return surface
end

"""
    reflect(panel::VortexRing)

Reflects a panel across the X-Z plane.
"""
@inline function reflect(panel::VortexRing)

    rtl = flipy(panel.rtr)
    rtc = flipy(panel.rtc)
    rtr = flipy(panel.rtl)
    rbl = flipy(panel.rbr)
    rbc = flipy(panel.rbc)
    rbr = flipy(panel.rbl)
    rcp = flipy(panel.rcp)
    ncp = flipy(panel.ncp)
    core_size = panel.core_size

    return VortexRing(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)
end

"""
    reflect(surface)

Reflects `surface` about the X-Z plane
"""
reflect(surface::AbstractMatrix) = reflect.(reverse(surface, dims=2))

"""
    left_center(panel::VortexRing)

Return the center of the left bound vortex on `panel`
"""
@inline left_center(panel::VortexRing) = (top_left(panel) + bottom_left(panel))/2

"""
    right_center(panel::VortexRing)

Return the center of the right bound vortex on `panel`
"""
@inline right_center(panel::VortexRing) = (top_right(panel) + bottom_right(panel))/2

"""
    top_vector(panel::VortexRing)

Returns the path of the top bound vortex for `panel`
"""
@inline top_vector(panel::VortexRing) = top_right(panel) - top_left(panel)

"""
    left_vector(panel)

Return the path of the left bound vortex for `panel`
"""
@inline left_vector(panel::VortexRing) = top_left(panel) - bottom_left(panel)

"""
    right_vector(panel)

Return the path of the right bound vortex for `panel`
"""
@inline right_vector(panel::VortexRing) = bottom_right(panel) - top_right(panel)

"""
    bottom_vector(panel)

Return the path of the bottom bound vortex for `panel`
"""
@inline bottom_vector(panel::VortexRing) = bottom_left(panel) - bottom_right(panel)

"""
    Wake{TF}

Wake panel element.

**Fields**
 - `rtl`: position of the left side of the top bound vortex
 - `rtr`: position of the right side of the top bound vortex
 - `rbl`: position of the left side of the bottom bound vortex
 - `rbr`: position of the right side of the bottom bound vortex
 - `core_size`: finite core size
 - `gamma`: circulation strength of the panel
"""
struct Wake{TF}
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
Wake(args...; kwargs...)

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

"""
    circulation_strength(panel::Wake)

Return circulation strength of the wake panel.
"""
@inline circulation_strength(panel::Wake) = panel.gamma

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
    TrefftzPanel{TF}

Panel in the Trefftz plane.
"""
struct TrefftzPanel{TF}
    rl::SVector{3,TF}
    rc::SVector{3,TF}
    rr::SVector{3,TF}
    Γ::TF
end

@inline Base.eltype(::Type{TrefftzPanel{TF}}) where TF = TF
@inline Base.eltype(::TrefftzPanel{TF}) where TF = TF

"""
    normal(panel::TrefftzPanel)

Return the normal vector of `panel`, including magnitude
"""
@inline function normal(panel::TrefftzPanel)

    rl = panel.rl
    rr = panel.rr

    dy = rr[2] - rl[2]
    dz = rr[3] - rl[3]

    nhat = SVector(0, -dz, dy)

    return nhat
end

"""
    trefftz_panels(surface[s], freestream, Γ)

Constructs a set of panels for Trefftz plane calculations
"""
trefftz_panels

function trefftz_panels(surface::AbstractMatrix, fs, Γ)

    TF = eltype(eltype(panels))

    panels = Vector{TrefftzPanel{TF}}(undef, size(surface, 2))

    return trefftz_panels!(panels, surface, fs, Γ)
end

function trefftz_panels(surfaces::AbstractVector{<:AbstractMatrix}, fs, Γ)

    TF = eltype(eltype(eltype(surfaces)))

    panels = [Vector{TrefftzPanel{TF}}(undef, size(surfaces[i], 2)) for i = 1:length(surfaces)]

    return trefftz_panels!(panels, surfaces, fs, Γ)
end

# single surface
function trefftz_panels!(panels, surface::AbstractMatrix, fs, Γ)

    n1, n2 = size(surface)

    rΓ = reshape(Γ, n1, n2)

    R = body_to_wind(fs)

    for j = 1:n2
        # use trailing edge panel for geometry
        rl = bottom_left(surface[end,j])
        rc = bottom_center(surface[end,j])
        rr = bottom_right(surface[end,j])

        # panel circulation is trailing edge panel circulation
        Γt = rΓ[end,j]

        # rotate into Trefftz plane frame of reference
        rl = R*rl
        rc = R*rc
        rr = R*rr

        # zero out x-components
        rl = SVector(0, rl[2], rl[3])
        rc = SVector(0, rc[2], rc[3])
        rr = SVector(0, rr[2], rr[3])

        panels[j] = TrefftzPanel(rl, rc, rr, Γt)
    end

    return panels
end

# multiple surfaces
function trefftz_panels!(panels, surfaces::AbstractVector{<:AbstractMatrix}, fs, Γ)

    nsurf = length(surfaces)

    # indices for keeping track of circulation
    iΓ = 0

    # loop through surfaces
    for i = 1:nsurf

        # number of panels on this surface
        n = length(surfaces[i])

        # create view into `panels` and `Γ` vector
        vΓ = view(Γ, iΓ+1:iΓ+n)

        # populate vector of Trefftz panels
        trefftz_panels!(panels[i], surfaces[i], fs, vΓ)

        # increment position in circulation index
        iΓ += n
    end

    return panels
end

"""
    trefftz_panel_induced_drag(receiving::TrefftzPanel, sending::TrefftzPanel; kwargs...)

Induced drag on `receiving` panel induced by `sending` panel.

# Keyword Arguments
 - `symmetric`: Flag indicating whether a mirror image of `sending` should be
    used when calculating the induced drag
"""
@inline function trefftz_panel_induced_drag(receiving::TrefftzPanel,
    sending::TrefftzPanel; symmetric)

    rl = sending.rl
    rr = sending.rr
    Γs = sending.Γ

    rc = receiving.rc
    nc = normal(receiving)
    Γr = receiving.Γ

    Di = vortex_induced_drag(rl, -Γs, rc, Γr, nc)
    Di += vortex_induced_drag(rr, Γs, rc, Γr, nc)
    if symmetric
        Di += vortex_induced_drag(flipy(rr), -Γs, rc, Γr, nc)
        Di += vortex_induced_drag(flipy(rl), Γs, rc, Γr, nc)
    end

    return Di
end

"""
    vortex_induced_drag(rj, Γj, ri, Γi, ni)

Return induced drag from vortex `j` induced on panel `i`
"""
@inline function vortex_induced_drag(rj, Γj, ri, Γi, ni)

    rij = ri - rj
    Vthetai = SVector(0, -Γj*rij[3], Γj*rij[2]) / (2*pi*(rij[2]^2 + rij[3]^2))
    Vn = -dot(Vthetai, ni)

    Di = RHO/2.0*Γi*Vn

    return Di
end
