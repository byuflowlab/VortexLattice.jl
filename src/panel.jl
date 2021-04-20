"""
    SurfacePanel{TF}

Lifting surface panel with attached vortex ring

**Fields**
 - `rtl`: position of the left side of the top bound vortex
 - `rtc`: position of the center of the top bound vortex
 - `rtr`: position of the right side of the top bound vortex
 - `rbl`: position of the left side of the bottom bound vortex
 - `rbc`: position of the center of the bottom bound vortex
 - `rbr`: position of the right side of the bottom bound vortex
 - `rcp`: position of the panel control point
 - `ncp`: normal vector at the panel control point
 - `core_size`: finite core size (for use when the finite core smoothing model is enabled)
 - `chord`: panel chordwise length (for determining unsteady forces)
"""
struct SurfacePanel{TF}
    rtl::SVector{3, TF}
    rtc::SVector{3, TF}
    rtr::SVector{3, TF}
    rbl::SVector{3, TF}
    rbc::SVector{3, TF}
    rbr::SVector{3, TF}
    rcp::SVector{3, TF}
    ncp::SVector{3, TF}
    core_size::TF
    chord::TF
end

"""
    SurfacePanel(rtl, rtr, rbl, rbr, rcp, ncp, core_size, chord; kwargs...)

Construct and return a vortex ring panel.

# Arguments
 - `rtl`: position of the left side of the top bound vortex
 - `rtr`: position of the right side of the top bound vortex
 - `rbl`: position of the left side of the bottom bound vortex
 - `rbr`: position of the right side of the bottom bound vortex
 - `rcp`: position of the panel control point
 - `ncp`: normal vector at the panel control point
 - `core_size`: finite core size (for use when the finite core smoothing model is enabled)
 - `chord`: panel chord length (for determining unsteady forces)

# Keyword Arguments
 - `rtc`: position of the center of the top bound vortex, defaults to `(rtl+rtr)/2`
 - `rbc`: position of the center of the bottom bound vortex, defaults to `(rbl+rbr)/2`
"""
SurfacePanel(args...; kwargs...)

function SurfacePanel(rtl, rtr, rbl, rbr, rcp, ncp, core_size, chord;
    rtc = (rtl+rtr)/2, # default to average of corners
    rbc = (rbl+rbr)/2) # default to average of corners

    return SurfacePanel(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size, chord)
end

function SurfacePanel(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size, chord)

    TF = promote_type(eltype(rtl), eltype(rtc), eltype(rtr), eltype(rbl), eltype(rbc),
        eltype(rbr), eltype(rcp), eltype(ncp), typeof(core_size), typeof(chord))

    return SurfacePanel{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size, chord)
end

@inline Base.eltype(::Type{SurfacePanel{TF}}) where TF = TF
@inline Base.eltype(::SurfacePanel{TF}) where TF = TF

"""
    top_left(panel::SurfacePanel)

Return the top left vertex of the vortex ring associated with `panel`
"""
@inline top_left(panel::SurfacePanel) = panel.rtl

"""
    top_center(panel::SurfacePanel)

Return the top center vertex of the vortex ring associated with `panel`
"""
@inline top_center(panel::SurfacePanel) = panel.rtc

"""
    top_right(panel::SurfacePanel)

Return the top right vertex of the vortex ring associated with `panel`
"""
@inline top_right(panel::SurfacePanel) = panel.rtr

"""
    bottom_left(panel::SurfacePanel)

Return the bottom left vertex of the vortex ring associated with `panel`
"""
@inline bottom_left(panel::SurfacePanel) = panel.rbl

"""
    bottom_center(panel::SurfacePanel)

Return the bottom center vertex of the vortex ring associated with `panel`
"""
@inline bottom_center(panel::SurfacePanel) = panel.rbc

"""
    bottom_right(panel::SurfacePanel)

Return the bottom right vertex of the vortex ring associated with `panel`
"""
@inline bottom_right(panel::SurfacePanel) = panel.rbr

"""
    controlpoint(panel::SurfacePanel)

Return the control point of `panel` (typically located at the 3/4 chord)
"""
@inline controlpoint(panel::SurfacePanel) = panel.rcp

"""
    normal(panel::SurfacePanel)

Return the normal vector of `panel` at the panel control point
"""
@inline normal(panel::SurfacePanel) = panel.ncp

"""
    get_core_size(panel::SurfacePanel)

Return the core size (smoothing parameter) corresponding to the vortex ring
associated with `panel`
"""
@inline get_core_size(panel::SurfacePanel) = panel.core_size

"""
    set_core_size(panel::SurfacePanel, core_size)

Return a copy of `panel` with core size equal to `core_size`.
"""
@inline set_core_size(panel::SurfacePanel{TF}, core_size) where TF = SurfacePanel{TF}(
    panel.rtl, panel.rtc, panel.rtr, panel.rbl, panel.rbc, panel.rbr, panel.rcp,
    panel.ncp, core_size, panel.chord)

"""
    translate(panel::SurfacePanel, r)

Return a copy of `panel` translated the distance specified by vector `r`
"""
@inline function translate(panel::SurfacePanel, r)

    rtl = panel.rtl + r
    rtc = panel.rtc + r
    rtr = panel.rtr + r
    rbl = panel.rbl + r
    rbc = panel.rbc + r
    rbr = panel.rbr + r
    rcp = panel.rcp + r
    ncp = panel.ncp
    core_size = panel.core_size
    chord = panel.chord

    return SurfacePanel(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size, chord)
end

"""
    translate(surface, r)

Return a copy of the panels in `surface` translated the distance specified by vector `r`
"""
translate(surface::AbstractMatrix, r) = translate.(surface, Ref(r))

"""
    translate!(surface, r)

Translate the panels in `surface` the distance specified by vector `r`
"""
function translate!(surface::AbstractMatrix, r)

    for i in eachindex(surface)
        surface[i] = translate(surface[i], r)
    end

    return surface
end

"""
    translate(surfaces, r)

Return a copy of the surfaces in `surfaces` translated the distance specified by vector `r`
"""
translate(surfaces::AbstractVector{<:AbstractMatrix}, r) = translate.(surfaces, Ref(r))

"""
    translate!(surfaces, r)

Translate the surfaces in `surfaces` the distance specified by vector `r`
"""
function translate!(surfaces::AbstractVector{AbstractMatrix}, r)

    for i in eachindex(surfaces)
        surfaces[i] = translate!(surfaces[i], r)
    end

    return surfaces
end

"""
    rotate(panel::SurfacePanel, R, r = [0,0,0])

Return a copy of `panel` rotated about point `r` using the rotation matrix `R`
"""
@inline function rotate(panel::SurfacePanel, R, r = (@SVector zeros(3)))

    rtl = R*(panel.rtl - r) + r
    rtc = R*(panel.rtc - r) + r
    rtr = R*(panel.rtr - r) + r
    rbl = R*(panel.rbl - r) + r
    rbc = R*(panel.rbc - r) + r
    rbr = R*(panel.rbr - r) + r
    rcp = R*(panel.rcp - r) + r
    ncp = R*panel.ncp
    core_size = panel.core_size
    chord = panel.chord

    return SurfacePanel(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size, chord)
end

"""
    rotate(surface, R, r = [0,0,0])

Return a copy of the panels in `surface` rotated about point `r` using the rotation matrix `R`
"""
rotate(surface::AbstractMatrix, R, r = (@SVector zeros(3))) = rotate.(surface, Ref(R), Ref(r))

"""
    rotate!(surface, R, r = [0,0,0])

Rotate the panels in `surface` about point `r` using the rotation matrix `R`
"""
function rotate!(surface::AbstractMatrix, R, r = (@SVector zeros(3)))

    for i in eachindex(surface)
        surface[i] = rotate(surface[i], R, r)
    end

    return surface
end

"""
    rotate(surfaces, R, r = [0,0,0])

Return a copy of the surfaces in `surfaces` rotated about point `r` using the
rotation matrix `R`
"""
rotate(surfaces::AbstractVector{<:AbstractMatrix}, R, r = (@SVector zeros(3))) =
    rotate.(surfaces, Ref(R), Ref(r))

"""
    rotate!(surfaces, R, r = [0,0,0])

Rotate the surfaces in `surfaces` about point `r` using the rotation matrix `R`
"""
function rotate!(surfaces::AbstractVector{AbstractMatrix}, R, r = (@SVector zeros(3)))

    for i in eachindex(surfaces)
        surfaces[i] = rotate!(surfaces[i], R, r)
    end

    return surfaces
end

"""
    reflect(panel::SurfacePanel)

Reflect `panel` across the X-Z plane.
"""
@inline function reflect(panel::SurfacePanel)

    rtl = flipy(panel.rtr)
    rtc = flipy(panel.rtc)
    rtr = flipy(panel.rtl)
    rbl = flipy(panel.rbr)
    rbc = flipy(panel.rbc)
    rbr = flipy(panel.rbl)
    rcp = flipy(panel.rcp)
    ncp = flipy(panel.ncp)
    core_size = panel.core_size
    chord = panel.chord

    return SurfacePanel(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size, chord)
end

"""
    reflect(surface)

Reflects the panels in `surface` across the X-Z plane
"""
reflect(surface::AbstractMatrix) = reflect.(reverse(surface, dims=2))

"""
    set_normal(panel::SurfacePanel, ncp)

Return a copy of `panel` with the new normal vector `ncp`
"""
@inline set_normal(panel::SurfacePanel, ncp) = SurfacePanel(panel.rtl, panel.rtc,
    panel.rtr, panel.rbl, panel.rbc, panel.rbr, panel.rcp, ncp, panel.core_size,
    panel.chord)

"""
    left_center(panel::SurfacePanel)

Return the center of the left bound vortex of the vortex ring associated with `panel`
"""
@inline left_center(panel::SurfacePanel) = (top_left(panel) + bottom_left(panel))/2

"""
    right_center(panel::SurfacePanel)

Return the center of the right bound vortex of the vortex ring associated with `panel`
"""
@inline right_center(panel::SurfacePanel) = (top_right(panel) + bottom_right(panel))/2

"""
    top_vector(panel::SurfacePanel)

Return the path of the top bound vortex of the vortex ring associated with `panel`
"""
@inline top_vector(panel::SurfacePanel) = top_right(panel) - top_left(panel)

"""
    left_vector(panel)

Return the path of the left bound vortex of the vortex ring associated with `panel`
"""
@inline left_vector(panel::SurfacePanel) = top_left(panel) - bottom_left(panel)

"""
    right_vector(panel)

Return the path of the right bound vortex of the vortex ring associated with `panel`
"""
@inline right_vector(panel::SurfacePanel) = bottom_right(panel) - top_right(panel)

"""
    bottom_vector(panel)

Return the path of the bottom bound vortex of the vortex ring associated with `panel`
"""
@inline bottom_vector(panel::SurfacePanel) = bottom_left(panel) - bottom_right(panel)

"""
    trailing_edge_points(surface[s])

Return points on the trailing edge of each surface.
"""
trailing_edge_points(surfaces) = [bottom_left.(surfaces[end, :])..., bottom_right(surfaces[end, end])]

"""
    repeated_trailing_edge_points(surface[s])

Generates a dictionary of the form `Dict((isurf, i) => [(jsurf1, j1), (jsurf2, j2)...]` which
defines repeated trailing edge points.  Trailing edge point `i` on surface
`isurf` is repeated on surface `jsurf1` at point `j1`, `jsurf2` at point `j2`,
and so forth.
"""
repeated_trailing_edge_points

repeated_trailing_edge_points(surface::AbstractMatrix) = repeated_trailing_edge_points([surface])

function repeated_trailing_edge_points(surfaces::AbstractVector{<:AbstractMatrix})

    nsurf = length(surfaces)
    points = trailing_edge_points.(surfaces)
    repeated_points = Dict{NTuple{2, Int}, Vector{NTuple{2, Int}}}()

    # loop through all points
    for isurf = 1:nsurf, i = 1:length(points[isurf])
        # start with empty vector
        repeats = NTuple{2, Int}[]
        # add all repeated points to the vector
        for jsurf = 1:nsurf, j = 1:length(points[jsurf])
            if (isurf, i) != (jsurf, j) && isapprox(points[isurf][i], points[jsurf][j])
                push!(repeats, (jsurf, j))
            end
        end
        # only store vector if nonempty
        if !isempty(repeats)
            repeated_points[(isurf, i)] = repeats
        end
    end

    return repeated_points
end

"""
    WakePanel{TF}

SurfacePanel used for modeling wakes.

**Fields**
 - `rtl`: position of the left side of the top bound vortex
 - `rtr`: position of the right side of the top bound vortex
 - `rbl`: position of the left side of the bottom bound vortex
 - `rbr`: position of the right side of the bottom bound vortex
 - `core_size`: finite core size (for use when the finite core smoothing model is enabled)
 - `gamma`: circulation strength of the panel
"""
struct WakePanel{TF}
    rtl::SVector{3, TF}
    rtr::SVector{3, TF}
    rbl::SVector{3, TF}
    rbr::SVector{3, TF}
    core_size::TF
    gamma::TF
end

WakePanel{TF}(panel::WakePanel) where TF = WakePanel{TF}(panel.rtl, panel.rtr,
    panel.rbl, panel.rbr, panel.core_size, panel.gamma)

"""
    WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)

Construct and return a wake panel.

**Arguments**
 - `rtl`: position of the left side of the top bound vortex
 - `rtr`: position of the right side of the top bound vortex
 - `rbl`: position of the left side of the bottom bound vortex
 - `rbr`: position of the right side of the bottom bound vortex
 - `core_size`: finite core size
 - `gamma`: circulation strength of the panel
"""
WakePanel(args...; kwargs...)

function WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)

    TF = promote_type(eltype(rtl), eltype(rtr), eltype(rbl), eltype(rbr), typeof(core_size), typeof(gamma))

    return WakePanel{TF}(rtl, rtr, rbl, rbr, core_size, gamma)
end

@inline Base.eltype(::Type{WakePanel{TF}}) where TF = TF
@inline Base.eltype(::WakePanel{TF}) where TF = TF

@inline top_left(panel::WakePanel) = panel.rtl

@inline top_right(panel::WakePanel) = panel.rtr

@inline bottom_left(panel::WakePanel) = panel.rbl

@inline bottom_right(panel::WakePanel) = panel.rbr

@inline get_core_size(panel::WakePanel) = panel.core_size

@inline set_core_size(panel::WakePanel{TF}, core_size) where TF = WakePanel{TF}(
    panel.rtl, panel.rtr, panel.rbl, panel.rbr, core_size, panel.gamma)

"""
    get_circulation(wakes; kwargs...)

Return the circulation strength of the wake panels in `wakes` as a vector
"""
get_circulation(wakes)

# single panel
@inline get_circulation(panel::WakePanel) = panel.gamma

# single wake
function get_circulation(wake::AbstractMatrix{<:WakePanel{TF}}) where TF
    Γw = Vector{TF}(undef, length(wake))
    return get_circulation!(Γw, wake)
end

# multiple wakes
function get_circulation(wakes::AbstractVector{<:AbstractMatrix{<:WakePanel{TF}}}) where TF
    Γw = Vector{TF}(undef, sum(length.(wakes)))
    return get_circulation!(Γw, wakes)
end

"""
    get_circulation!(wakes; kwargs...)

Pre-allocated version of [`get_circulation`](@ref)
"""
get_circulation!(wakes)

# single wake
function get_circulation!(Γw, wake::AbstractMatrix{<:WakePanel}; iwake = size(wake, 1))
    # wake (and wake storage) dimensions
    iw = iwake
    nw, ns = size(wake)
    # reshape circulation strength to more convenient form
    Γw = reshape(Γw, nw, ns)
    # return circulation strength for defined wake panels
    for i = 1:iw, j = 1:ns
        Γw[i,j] = get_circulation(wake[i,j])
    end
    # return zero for non-existant wake panels
    Γw[iw+1:nw,:] .= 0.0
    return Γw
end

# multiple wakes
function get_circulation!(Γw, wakes::AbstractVector{<:AbstractMatrix{<:WakePanel}})
    iΓ = 0
    for wake in wakes
        N = length(wake)
        vΓw = view(Γw, iΓ+1:iΓ+N)
        get_circulation!(vΓw, wake)
        iΓ += N
    end
    return Γw
end

"""
    get_vertices(wake_shedding_locations, wakes)

Return the vertices of the wake panels, given the wake shedding locations and
the wake panels.
"""
get_vertices

# single wake
function get_vertices(wake_shedding_locations, wake::AbstractMatrix{WakePanel{TF}};
    iwake = size(wake, 1)) where TF
    nw, ns = size(wake)
    ζw = Matrix{SVector{3, TF}}(undef, iwake+1, ns+1)
    return get_vertices!(ζw, wake_shedding_locations, wake; iwake)
end

# multiple wakes
function get_vertices(wake_shedding_locations, wakes::AbstractVector{<:AbstractMatrix{WakePanel{TF}}};
    iwake = size.(wake, 1)) where TF
    ζw = Vector{Matrix{SVector{3, TF}}}(undef, length(wakes))
    for i = 1:length(wakes)
        ζw[i] = get_vertices(wake_shedding_locations[i], wakes[i]; iwake=iwake[i])
    end
    return ζw
end

"""
    get_vertices!(wake_shedding_locations, wakes)

Pre-allocated version of [`get_vertices`](@ref)
"""
get_vertices!

# single wake
function get_vertices!(ζw, wake_shedding_locations, wake::AbstractMatrix{<:WakePanel};
    iwake = size(wake, 1))
    nw, ns = size(wake)
    # first row corresponds to wake shedding locations
    for j = 1:ns+1
        ζw[1,j] = wake_shedding_locations[j]
    end
    # remaining rows correspond to wake panel coordinates
    for j = 1:ns+1, i = 2:iwake+1
        if i == nw+1 && j == ns+1
            ζw[i,j] = bottom_right(wake[end,end])
        elseif i == nw+1
            ζw[i,j] = bottom_left(wake[end,j])
        elseif j == ns+1
            ζw[i,j] = top_right(wake[i,end])
        else
            ζw[i,j] = top_left(wake[i,j])
        end
    end
    return ζw
end

# multiple wakes
function get_vertices!(ζw, wake_shedding_locations, wakes::AbstractVector{<:AbstractMatrix{<:WakePanel}};
    iwake = size.(wakes, 1))
    for i = 1:length(wakes)
        get_vertices!(ζw[i], wake_shedding_locations[i], wakes[i]; iwake=iwake[i])
    end
    return Γw
end

"""
    set_circulation(wakes, Γw)

Return a copy of the panels in `wakes` with the circulation strength of each panel
set to the circulation strengths in `Γw`.
"""
set_circulation

@inline function set_circulation(panel::WakePanel, gamma)

    rtl = panel.rtl
    rtr = panel.rtr
    rbl = panel.rbl
    rbr = panel.rbr
    core_size = panel.core_size

    return WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
end

# single surface
function set_circulation(wake::AbstractMatrix{<:WakePanel}, Γw; kwargs...)
    return set_circulation!(deepcopy(wake), Γw; kwargs...)
end

# multiple surfaces
function set_circulation(wake::AbstractVector{<:AbstractMatrix{<:WakePanel}}, Γw; kwargs...)
    return set_circulation!(deepcopy(wake), Γw; kwargs...)
end

"""
    set_circulation!(wakes, Γw)

Set the circulation strength of the panels in `wakes` to the circulations strengths
in `Γw`
"""
set_circulation!

# single surface
function set_circulation!(wake::AbstractMatrix{<:WakePanel}, Γw; iwake = size(wakes, 1))
    nw, ns = size(wake)
    vΓw = reshape(Γw, nw, ns)
    for i = 1:iwake
        for j = 1:ns
            wake[i,j] = set_circulation(wake[i,j], vΓw[i,j])
        end
    end
    return wake
end

# multiple surfaces
function set_circulation!(wakes, Γw; iwake = size.(wakes, 1))
    iΓ = 0
    for isurf = 1:length(wakes)
        N = length(wakes[isurf])
        vΓw = view(Γw, iΓ+1:iΓ+N)
        set_circulation!(wakes[isurf], vΓw; iwake = iwake[isurf])
        iΓ += N
    end
    return wakes
end

@inline function translate(panel::WakePanel, r)

    rtl = panel.rtl + r
    rtr = panel.rtr + r
    rbl = panel.rbl + r
    rbr = panel.rbr + r
    core_size = panel.core_size
    gamma = panel.gamma

    return WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
end

@inline function reflect(panel::WakePanel)

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

"""
    PanelProperties

Holds surface panel properties calculated during the vortex lattice method analysis.

**Fields**
 - `gamma`: Vortex ring circulation strength, normalized by the reference velocity
 - `velocity`: Local velocity at the panel's bound vortex center, normalized by
    the reference velocity
 - `cfb`: Net force on the panel's bound vortex, as calculated using the
    Kutta-Joukowski theorem, normalized by the reference dynamic pressure and area
 - `cfl`: Force on the left bound vortex from this panel's vortex ring, as
    calculated by the Kutta-Joukowski theorem, normalized by the reference
    dynamic pressure and area
 - `cfr`: Force on the right bound vortex from this panel's vortex ring, as
    calculated by the Kutta-Joukowski theorem, normalized by the reference
    dynamic pressure and area
"""
struct PanelProperties{TF}
    gamma::TF
    velocity::SVector{3, TF}
    cfb::SVector{3, TF}
    cfl::SVector{3, TF}
    cfr::SVector{3, TF}
end

function PanelProperties(gamma, velocity, cfb, cfl, cfr)

    TF = promote_type(typeof(gamma), typeof(velocity), eltype(cfb), eltype(cfl),
        eltype(cfr))

    return PanelProperties{TF}(gamma, velocity, cfb, cfl, cfr)
end

Base.eltype(::Type{PanelProperties{TF}}) where TF = TF
Base.eltype(::PanelProperties{TF}) where TF = TF

"""
    update_surface_panels!(surface_panels, surfaces; kwargs...)

Update the surface panels in `surface_panels` to correspond to the surface
panels/grids in `surfaces`.

# Keyword Arguments
 - `fcore`: function for setting the finite core size when generating surface
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width. Only used if keyword
    argument `surfaces` is provided and is a collection of grids.
 - `preserve_core_size`: Flag indicating whether the finite core size should be
    preserved from the previous set of panels. Defaults to `true`. Only used if
    keyword argument `surfaces` is provided.
"""
update_surface_panels!

# multiple surfaces
function update_surface_panels!(surface_panels, surfaces; kwargs...)
    nsurf = length(surfaces)
    for isurf = 1:nsurf
        update_surface_panels!(surface_panels[isurf], surfaces[isurf]; kwargs...)
    end
    return surface_panels
end

# single surface, surface panels input
function update_surface_panels!(surface_panels, surface::AbstractMatrix;
    fcore = (c, Δs) -> 1e-3, preserve_core_size = true)
    for (i, panel) in enumerate(surface)
        if preserve_core_size && isassigned(surface_panels, i)
            surface_panels[i] = set_core_size(panel, get_core_size(surface_panels[i]))
        else
            surface_panels[i] = surface[i]
        end
    end
    return surface_panels
end

# single surface, grid input
function update_surface_panels!(surface_panels, surface::AbstractArray{<:Number, 3};
    fcore = (c, Δs) -> 1e-3, preserve_core_size = true)

    TF = eltype(eltype(surface_panels))

    nc, ns = size(surface_panels) # number of chordwise and spanwise panels

    # populate each panel
    for j = 1:ns

        # get leading edge panel corners
        r1n = SVector(surface[1,1,j], surface[2,1,j], surface[3,1,j]) # top left
        r2n = SVector(surface[1,1,j+1], surface[2,1,j+1], surface[3,1,j+1]) # top right
        r3n = SVector(surface[1,2,j], surface[2,2,j], surface[3,2,j]) # bottom left
        r4n = SVector(surface[1,2,j+1], surface[2,2,j+1], surface[3,2,j+1]) # bottom right

        # also get chord length for setting finite core size
        cl = norm(surface[:,end,j] - surface[:,1,j])
        cr = norm(surface[:,end,j+1] - surface[:,1,j+1])
        c = (cl + cr)/2

        for i = 1:nc-1

            # grid corners of current panel
            r1 = r1n # top left of panel
            r2 = r2n # top right of panel
            r3 = r3n # bottom left of panel
            r4 = r4n # bottom right of panel

            # grid corners of next panel
            r1n = r3 # top left
            r2n = r4 # top right
            r3n = SVector(surface[1,i+2,j], surface[2,i+2,j], surface[3,i+2,j]) # bottom left
            r4n = SVector(surface[1,i+2,j+1], surface[2,i+2,j+1], surface[3,i+2,j+1]) # bottom right

            # top left corner of ring vortex
            rtl = linearinterp(0.25, r1, r3)

            # top right corner of ring vortex
            rtr = linearinterp(0.25, r2, r4)

            # center of top bound vortex
            rtc = (rtl + rtr)/2

            # bottom left corner of ring vortex
            rbl = linearinterp(0.25, r1n, r3n)

            # bottom right corner of ring vortex
            rbr = linearinterp(0.25, r2n, r4n)

            # center of bottom bound vortex
            rbc = (rbl + rbr)/2

            # control point
            rtop = linearinterp(0.5, r1, r2)
            rbot = linearinterp(0.5, r3, r4)
            rcp = linearinterp(0.75, rtop, rbot)

            # surface normal
            ncp = cross(rcp - rtr, rcp - rtl)
            ncp /= norm(ncp)

            # set finite core size
            if preserve_core_size && isassigned(surface_panels, i, j)
                core_size = get_core_size(surface_panels[i,j])
            else
                Δs = sqrt((rtr[2]-rtl[2])^2 + (rtr[3]-rtl[3])^2)
                core_size = fcore(c, Δs)
            end

            # get chord length of current panel
            chord = norm((r1 + r2)/2 - (r3 + r4)/2)

            surface_panels[i,j] = SurfacePanel{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp,
                ncp, core_size, chord)
        end

        # grid corners of current panel
        r1 = r1n # top left
        r2 = r2n # top right
        r3 = r3n # bottom left
        r4 = r4n # bottom right

        # top left corner of ring vortex
        rtl = linearinterp(0.25, r1, r3)

        # top right corner of ring vortex
        rtr = linearinterp(0.25, r2, r4)

        # center of top bound vortex
        rtc = (rtl + rtr)/2

        # bottom left corner of ring vortex
        rbl = r3

        # bottom right corner of ring vortex
        rbr = r4

        # center of bottom bound vortex
        rbc = (rbl + rbr)/2

        # control point
        rtop = linearinterp(0.5, r1, r2)
        rbot = linearinterp(0.5, r3, r4)
        rcp = linearinterp(0.75, rtop, rbot)

        # surface normal
        ncp = cross(rcp - rtr, rcp - rtl)
        ncp /= norm(ncp)

        # set finite core size
        if preserve_core_size && isassigned(surface_panels, nc, j)
            core_size = get_core_size(surface_panels[end,j])
        else
            Δs = sqrt((rtr[2]-rtl[2])^2 + (rtr[3]-rtl[3])^2)
            core_size = fcore(c, Δs)
        end

        # get chord length of current panel
        chord = norm((r1 + r2)/2 - (r3 + r4)/2)

        surface_panels[end, j] = SurfacePanel{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp,
            core_size, chord)
    end

    return surface_panels
end


"""
    initialize_wake_panels!(wake_panels, wake_shedding_locations; kwargs...)

Initialize the undefined wake panels in `wake_panels`.

# Keyword Arguments:
 - `iwake`: Number of chordwise wake panels on each surface which have already
    been initialized.  Defaults to zero wake panels
"""
initialize_wake_panels!

# multiple surfaces
function initialize_wake_panels!(wake_panels, wake_shedding_locations;
    iwake = zeros(length(wake_shedding_locations)), kwargs...)

    for isurf = 1:length(wake_shedding_locations)
        initialize_wake_panels!(wake_panels[isurf], wake_shedding_locations[isurf];
            iwake=iwake[isurf], kwargs...)
    end

    return wake_panels
end

# single wake, wake panels input
function initialize_wake_panels!(wake_panels::AbstractMatrix, wake_shedding_locations; iwake = 0)
    nw, ns = size(wake_panels)
    for j = 1:ns
        rtl = rbl = wake_shedding_locations[j]
        rtr = rbr = wake_shedding_locations[j+1]
        for i = iwake+1:nw
            wake_panels[i,j] = WakePanel(rtl, rtr, rbl, rbr, 0.0, 0.0)
        end
    end
    return wake_panels
end

"""
    update_wake_panels!(wake_panels, wakes; kwargs...)

Update the wake panels in `wake_panels` to correspond to the wake panels/vertices
in `wakes`.

# Keyword Arguments:
 - `helmholtz: Flag indicating whether circulation should be updated based on
    Helmholtz' fourth vortex theorem. (The strength of a vortex tube is constant).
    Defaults to `true`
"""
update_wake_panels!

# multiple surfaces, wake panel input
function update_wake_panels!(wake_panels, wakes::AbstractVector{<:AbstractMatrix{<:WakePanel}};
    iwake = size.(wakes, 1), kwargs...)
    for isurf = 1:length(wakes)
        update_wake_panels!(wake_panels[isurf], wakes[isurf]; iwake = iwake[isurf], kwargs...)
    end
end

function update_wake_panels!(wake_panels, wakes::AbstractVector{<:AbstractMatrix{<:SVector}};
    iwake = size.(wakes, 1) .- 1, kwargs...)
    for isurf = 1:length(wakes)
        update_wake_panels!(wake_panels[isurf], wakes[isurf]; iwake = iwake[isurf], kwargs...)
    end
end

# multiple surfaces, grid input
function update_wake_panels!(wake_panels, wakes::AbstractVector{<:AbstractArray{<:Number, 3}};
    iwake = size.(wakes, 2) .- 1, kwargs...)
    for isurf = 1:length(wakes)
        update_wake_panels!(wake_panels[isurf], wakes[isurf]; iwake = iwake[isurf], kwargs...)
    end
end

# single wake, wake panels input
function update_wake_panels!(wake_panels, wake::AbstractMatrix{<:WakePanel}; iwake = size(wake, 1),
    helmholtz = true)

    nw, ns = size(wake)

    for i = 1:iwake, j = 1:ns
        if isassigned(wake_panels, i, j)
            if helmholtz
                # get original vortex filament length
                rtl = top_left(wake_panels[i,j])
                rtr = top_right(wake_panels[i,j])
                rbl = bottom_left(wake_panels[i,j])
                rbr = bottom_right(wake_panels[i,j])
                lt = norm(rtr - rtl)
                ll = norm(rtl - rbl)
                lr = norm(rtr - rbr)
                lb = norm(rbr - rbl)
                l = lt + ll + lr + lb
                # set new wake panel coordinates
                rtl = top_left(wake[i,j])
                rtr = top_right(wake[i,j])
                rbl = bottom_left(wake[i,j])
                rbr = bottom_right(wake[i,j])
                # get new vortex filament length
                lt = norm(rtr - rtl)
                ll = norm(rtl - rbl)
                lr = norm(rtr - rbr)
                lb = norm(rbr - rbl)
                new_l = new_lt + new_ll + new_lr + new_lb
                # new core size is old core size
                core_size = get_core_size(wake_panels[i,j])
                # set new vortex filament strength
                gamma = get_circulation(wake_panels[i,j]) * l / new_l
            else
                # set new wake panel coordinates
                rtl = top_left(wake[i,j])
                rtr = top_right(wake[i,j])
                rbl = bottom_left(wake[i,j])
                rbr = bottom_right(wake[i,j])
                # new core size is old core size
                core_size = get_core_size(wake_panels[i,j])
                # new circulation strength is old circulation strength
                gamma = get_circulation(wake_panels[i,j])
            end
        else
            # set new wake panel coordinates
            rtl = top_left(wake[i,j])
            rtr = top_right(wake[i,j])
            rbl = bottom_left(wake[i,j])
            rbr = bottom_right(wake[i,j])
            # set new core size to zero
            core_size = 0.0
            # set new vortex filament strength
            gamma = 0.0
        end
        # update wake panel
        wake_panels[i,j] = WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
    end

    return wake_panels
end

# single wake, grid input
function update_wake_panels!(wake_panels, wake::AbstractMatrix{<:SVector};
    iwake = size(wake, 1) - 1, helmholtz = true)

    nw = size(wake, 1) - 1
    ns = size(wake, 2) - 1

    for i = 1:iwake, j = 1:ns
        if isassigned(wake_panels, i, j)
            if helmholtz
                # get original vortex filament length
                rtl = top_left(wake_panels[i,j])
                rtr = top_right(wake_panels[i,j])
                rbl = bottom_left(wake_panels[i,j])
                rbr = bottom_right(wake_panels[i,j])
                lt = norm(rtr - rtl)
                ll = norm(rtl - rbl)
                lr = norm(rtr - rbr)
                lb = norm(rbr - rbl)
                l = lt + ll + lr + lb
                # set new wake panel coordinates
                rtl = wake[i,j]
                rtr = wake[i,j+1]
                rbl = wake[i+1,j]
                rbr = wake[i+1,j+1]
                # get new vortex filament length
                lt = norm(rtr - rtl)
                ll = norm(rtl - rbl)
                lr = norm(rtr - rbr)
                lb = norm(rbr - rbl)
                new_l = lt + ll + lr + lb
                # new core size is old core size
                core_size = get_core_size(wake_panels[i,j])
                # set new vortex filament strength
                gamma = get_circulation(wake_panels[i,j]) * l / new_l
            else
                # set new wake panel coordinates
                rtl = wake[i,j]
                rtr = wake[i,j+1]
                rbl = wake[i+1,j]
                rbr = wake[i+1,j+1]
                # new core size is old core size
                core_size = get_core_size(wake_panels[i,j])
                # new circulation strength is old circulation strength
                gamma = get_circulation(wake_panels[i,j])
            end
        else
            # set new wake panel coordinates
            rtl = wake[i,j]
            rtr = wake[i,j+1]
            rbl = wake[i+1,j]
            rbr = wake[i+1,j+1]
            # set new core size to zero
            core_size = 0.0
            # set new vortex filament strength
            gamma = 0.0
        end
        # update wake panel
        wake_panels[i,j] = WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
    end

    return wake_panels
end

# single wake, grid input
function update_wake_panels!(wake_panels, wake::AbstractArray{<:Number, 3};
    iwake = size(wake, 2) - 1, helmholtz = true)

    nw = size(wake, 2) - 1
    ns = size(wake, 3) - 1

    for i = 1:iwake, j = 1:ns
        if isassigned(wake_panels, i, j)
            if helmholtz
                # get original vortex filament length
                rtl = top_left(wake_panels[i,j])
                rtr = top_right(wake_panels[i,j])
                rbl = bottom_left(wake_panels[i,j])
                rbr = bottom_right(wake_panels[i,j])
                lt = norm(rtr - rtl)
                ll = norm(rtl - rbl)
                lr = norm(rtr - rbr)
                lb = norm(rbr - rbl)
                l = lt + ll + lr + lb
                # set new wake panel coordinates
                rtl = SVector(wake[1,i,j], wake[2,i,j], wake[3,i,j])
                rtr = SVector(wake[1,i,j+1], wake[2,i,j+1],wake[3,i,j+1])
                rbl = SVector(wake[1,i+1,j], wake[2,i+1,j], wake[3,i+1,j])
                rbr = SVector(wake[1,i+1,j+1], wake[2,i+1,j+1], wake[3,i+1,j+1])
                # get new vortex filament length
                lt = norm(rtr - rtl)
                ll = norm(rtl - rbl)
                lr = norm(rtr - rbr)
                lb = norm(rbr - rbl)
                new_l = lt + ll + lr + lb
                # new core size is old core size
                core_size = get_core_size(wake_panels[i,j])
                # set new vortex filament strength
                gamma = get_circulation(wake_panels[i,j]) * l / new_l
            else
                # set new wake panel coordinates
                rtl = SVector(wake[1,i,j], wake[2,i,j], wake[3,i,j])
                rtr = SVector(wake[1,i,j+1], wake[2,i,j+1],wake[3,i,j+1])
                rbl = SVector(wake[1,i+1,j], wake[2,i+1,j], wake[3,i+1,j])
                rbr = SVector(wake[1,i+1,j+1], wake[2,i+1,j+1], wake[3,i+1,j+1])
                # new core size is old core size
                core_size = get_core_size(wake_panels[i,j])
                # new circulation strength is old circulation strength
                gamma = get_circulation(wake_panels[i,j])
            end
        else
            # set new wake panel coordinates
            rtl = SVector(wake[1,i,j], wake[2,i,j], wake[3,i,j])
            rtr = SVector(wake[1,i,j+1], wake[2,i,j+1],wake[3,i,j+1])
            rbl = SVector(wake[1,i+1,j], wake[2,i+1,j], wake[3,i+1,j])
            rbr = SVector(wake[1,i+1,j+1], wake[2,i+1,j+1], wake[3,i+1,j+1])
            # set new core size to zero
            core_size = 0.0
            # set new vortex filament strength
            gamma = 0.0
        end
        # update wake panel
        wake_panels[i,j] = WakePanel(rtl, rtr, rbl, rbr, core_size, gamma)
    end

    return wake_panels
end

"""
    vortex_filament_length(panels)

Calculate the vortex filament length for each surface/wake panel in `panels`.
"""
vortex_filament_length

@inline function vortex_filament_length(rtl, rtr, rbl, rbr)
    lt = norm(rtl - rtr)
    lb = norm(rbl - rbr)
    ll = norm(rtl - rbl)
    lr = norm(rtr - rbr)
    return lt + lb + ll + lr
end

@inline function vortex_filament_length(panel)
    rtl = top_left(panel)
    rtr = top_right(panel)
    rbl = bottom_left(panel)
    rbr = bottom_right(panel)
    lt = norm(rtl - rtr)
    lb = norm(rbl - rbr)
    ll = norm(rtl - rbl)
    lr = norm(rtr - rbr)
    return lt + lb + ll + lr
end

@inline function vortex_filament_length(panels::AbstractMatrix)
    return vortex_filament_length.(reshape(panels, :))
end

@inline function vortex_filament_length(vertices::AbstractArray{TF, 3}) where TF
    nc = size(vertices, 2) - 1
    ns = size(vertices, 3) - 1
    N = nc*ns
    l = Vector{TF}(undef, N)
    return vortex_filament_length!(l, vertices)
end

@inline function vortex_filament_length(panels::AbstractVector{<:AbstractMatrix})
    TF = eltype(eltype(eltype(panels)))
    N = sum(length.(panels))
    l = Vector{TF}(undef, N)
    return vortex_filament_length!(l, panels)
end

@inline function vortex_filament_length(vertices::AbstractVector{<:AbstractArray{TF, 3}}) where TF
    nc = size.(vertices, 2) .- 1
    ns = size.(vertices, 3) .- 1
    N = sum(nc .* ns)
    l = Vector{TF}(undef, N)
    return vortex_filament_length!(l, vertices)
end

"""
    vortex_filament_length!(l, panels)

In-place version of `vortex_filament_length`
"""
vortex_filament_length!

function vortex_filament_length!(l, panels::AbstractMatrix)
    for (i, panel) in enumerate(panels)
        l[i] = vortex_filament_length(panel)
    end
    return l
end

function vortex_filament_length!(l, vertices::AbstractMatrix{<:SVector})
    nc = size(vertices, 1) - 1
    ns = size(vertices, 2) - 1
    rl = reshape(l, nc, ns)
    for i = 1:nc, j = 1:ns
        rtl = vertices[i,j]
        rtr = vertices[i,j+1]
        rbl = vertices[i+1,j]
        rbr = vertices[i+1,j+1]
        rl[i,j] = vortex_filament_length(rtl, rtr, rbl, rbr)
    end
    return l
end

function vortex_filament_length!(l, vertices::AbstractArray{<:Number, 3})
    nc = size(vertices, 2) - 1
    ns = size(vertices, 3) - 1
    rl = reshape(l, nc, ns)
    for i = 1:nc, j = 1:ns
        rtl = SVector(vertices[1,i,j], vertices[2,i,j], vertices[3,i,j])
        rtr = SVector(vertices[1,i,j+1], vertices[2,i,j+1], vertices[3,i,j+1])
        rbl = SVector(vertices[1,i+1,j], vertices[2,i+1,j], vertices[3,i+1,j])
        rbr = SVector(vertices[1,i+1,j+1], vertices[2,i+1,j+1], vertices[3,i+1,j+1])
        rl[i,j] = vortex_filament_length(rtl, rtr, rbl, rbr)
    end
    return l
end

function vortex_filament_length!(l, panels::AbstractVector{<:AbstractMatrix})
    nsurf = length(panels)
    il = 0
    for isurf = 1:nsurf
        N = length(panels[isurf])
        vl = view(l, il+1:il+N)
        vortex_filament_length!(vl, panels[isurf])
        il += N
    end
    return l
end

function vortex_filament_length!(l, vertices::AbstractVector{<:AbstractMatrix{<:SVector}})
    nsurf = length(vertices)
    il = 0
    for isurf = 1:nsurf
        nc = size(vertices[isurf], 1) - 1
        ns = size(vertices[isurf], 2) - 1
        N = nc*ns
        vl = view(l, il+1:il+N)
        vortex_filament_length!(vl, vertices[isurf])
        il += N
    end
    return l
end

function vortex_filament_length!(l, vertices::AbstractVector{<:AbstractArray{<:Number, 3}})
    nsurf = length(vertices)
    il = 0
    for isurf = 1:nsurf
        nc = size(vertices[isurf], 2) - 1
        ns = size(vertices[isurf], 3) - 1
        N = nc*ns
        vl = view(l, il+1:il+N)
        vortex_filament_length!(vl, vertices[isurf])
        il += N
    end
    return l
end
