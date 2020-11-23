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

# single surface, horseshoe
function trefftz_panels!(panels, surface::AbstractMatrix{<:Horseshoe}, fs, Γ)

    n1, n2 = size(surface)

    rΓ = reshape(Γ, n1, n2)

    R = body_to_wind(fs)

    for j = 1:n2
        # use trailing edge panel for geometry
        rl = bottom_left(surface[end,j])
        rc = bottom_center(surface[end,j])
        rr = bottom_right(surface[end,j])

        # sum circulation across chordwise strip to get panel circulation
        Γt = zero(eltype(Γ))
        for i = 1:n1
            Γt += rΓ[i,j]
        end

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

# single surface, vortex ring
function trefftz_panels!(panels, surface::AbstractMatrix{<:Ring}, fs, Γ)

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
    trefftz_panel_induced_drag(receiving::TrefftzPanel, sending::TrefftzPanel, same_id, symmetric)

Induced drag on `receiving` panel induced by `sending` panel.
"""
@inline function trefftz_panel_induced_drag(receiving::TrefftzPanel, sending::TrefftzPanel,
    symmetric)

    rl = sending.rl
    rr = sending.rr
    Γs = sending.Γ

    rc = receiving.rc
    nc = normal(receiving)
    Γr = receiving.Γ

    Di = vortex_induced_drag(rl, -Γs, rc, Γr, nc)
    Di += vortex_induced_drag(rr, Γs, rc, Γr, nc)
    if symmetric && not_on_symmetry_plane(rl, rr)
        Di += vortex_induced_drag(flipy(rr), -Γs, rc, Γr, nc)
        Di += vortex_induced_drag(flipy(rl), Γs, rc, Γr, nc)
    end

    return Di
end

"""
    far_field_drag(system, surface, reference, freestream; kwargs...)

Computes induced drag using the Trefftz plane (far field method).

Note that this function assumes that the circulation distribution has already
been computed and is present in `system`

# Arguments
 - `system`: Pre-allocated system properties
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see `Reference`)
 - `freestream`: Freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: Flag indicating whether a mirror image of the panels in `surface`,
    should be used when calculating induced velocities, defaults to `false`
"""
function far_field_drag(system, surface::AbstractMatrix, ref, fs;
    symmetric = false)

    # unpack system
    Γ = system.gamma
    trefftz = system.trefftz[1]

    trefftz_panels!(trefftz, surface, fs, Γ)

    return far_field_drag(trefftz, trefftz, ref, fs, symmetric)
end

"""
    far_field_drag(system, surface, reference, freestream; kwargs...)

Computes induced drag using the Trefftz plane (far field method).

Note that this function assumes that the circulation distribution has already
been computed and is present in `system`

# Arguments
 - `system`: Pre-allocated system properties
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see `Reference`)
 - `freestream`: Freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: Flag indicating whether a mirror image of the panels in `surface`,
    should be used when calculating induced velocities, defaults to `false`
"""
function far_field_drag(system, surfaces::AbstractVector{<:AbstractMatrix}, ref, fs;
    symmetric = fill(false, length(surfaces)))

    # unpack system
    Γ = system.gamma
    trefftz = system.trefftz

    # construct trefftz panels
    trefftz_panels!(trefftz, surfaces, fs, Γ)

    # perform far field analysis
    nsurf = length(surfaces)
    CD = zero(eltype(eltype(eltype(surfaces))))
    for i = 1:nsurf, j = 1:nsurf
        CD += far_field_drag(trefftz[i], trefftz[j], ref, fs, symmetric[j])
    end

    return CD
end

# one surface on another surface
function far_field_drag(receiving, sending, ref, fs, symmetric)

    TF = promote_type(eltype(eltype(receiving)), eltype(eltype(sending)), eltype(ref), eltype(fs))

    Nr = length(receiving)
    Ns = length(sending)

    # add up drag
    Di = zero(TF)
    for j = 1:Ns
        for i = 1:Nr
            Di += trefftz_panel_induced_drag(receiving[i], sending[j], symmetric)
        end
    end

    # apply symmetry
    if symmetric
        Di *= 2
    end

    # normalize
    CDi = Di / (QINF*ref.S)

    return CDi
end

# --- internal functions --- #

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
