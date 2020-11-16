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

# single surface, horseshow
function trefftz_panels(surface::AbstractMatrix{<:Horseshoe}, fs, Γ)

    TF = promote_type(eltype(eltype(surface)), eltype(fs), eltype(Γ))

    n1, n2 = size(surface)

    rΓ = reshape(Γ, n1, n2)

    R = body_to_wind(fs)

    panels = Vector{TrefftzPanel{TF}}(undef, n2)

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
function trefftz_panels(surface::AbstractMatrix{<:Ring}, fs, Γ)

    TF = promote_type(eltype(eltype(surface)), eltype(fs), eltype(Γ))

    n1, n2 = size(surface)

    rΓ = reshape(Γ, n1, n2)

    R = body_to_wind(fs)

    panels = Vector{TrefftzPanel{TF}}(undef, n2)

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
function trefftz_panels(surfaces::AbstractVector{<:AbstractMatrix}, fs, Γ)

    TF = promote_type(eltype(eltype(eltype(surfaces))), eltype(fs), eltype(Γ))

    nsurf = length(surfaces)

    # indices for keeping track of circulation
    iΓ = 0

    panels = Vector{Vector{TrefftzPanel{TF}}}(undef, nsurf)

    # loop through surfaces
    for i = 1:nsurf

        # number of panels on this surface
        n = length(surfaces[i])

        # create view into `panels` and `Γ` vector
        vΓ = view(Γ, iΓ+1:iΓ+n)

        # populate vector of Trefftz panels
        panels[i] = trefftz_panels(surfaces[i], fs, vΓ)

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
    far_field_drag(surface[s], ref, fs, symmetric, Γ)

Computes induced drag using the Trefftz plane (far field method).
"""
far_field_drag

# one surface
function far_field_drag(surface::AbstractMatrix, ref, fs, symmetric, Γ)

    panels = trefftz_panels(surface, fs, Γ)

    return far_field_drag(panels, panels, ref, fs, symmetric)
end

# multiple surfaces
function far_field_drag(surfaces::AbstractVector{<:AbstractMatrix}, ref, fs,
    symmetric, Γ)

    nsurf = length(surfaces)

    panels = trefftz_panels(surfaces, fs, Γ)

    CD = zero(eltype(eltype(eltype(surfaces))))
    for i = 1:nsurf, j = 1:nsurf
        CD += far_field_drag(panels[i], panels[j], ref, fs, symmetric)
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
