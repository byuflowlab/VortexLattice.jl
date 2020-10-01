"""
    PanelOutputs

Panel specific outputs from the VLM analysis.

**Fields**
 - `gamma`: Panel circulation strength (normalized by the freestream velocity)
 - `v`: Local velocity at the panel's bound vortex midpoint (typically the quarter-chord, normalized)
 - `cf`: Forces on this panel, normalized by `QINF*S` where `QINF=0.5` is the dynamic
    pressure used in this package and `S` is the user-provided reference area.
"""
struct PanelOutputs{TF}
    gamma::TF
    v::SVector{3, TF}
    cf::SVector{3, TF}
end

function PanelOutputs(gamma, v, cf)

    TF = promote_type(typeof(gamma), typeof(v), eltype(cf))

    return PanelOutputs{TF}(cf, v, gamma)
end

Base.eltype(::Type{PanelOutputs{TF}}) where TF = TF
Base.eltype(::PanelOutputs{TF}) where TF = TF

# --------------------------

"""
    Outputs

Outputs from the VLM analysis

**Arguments**
- `CF`: force coeffients (CD, CY, CL) normalized by reference area
- `CM`: moment coeffients (CMX, CMY, CMZ) normalized by reference area
- `panels`: vector of outputs for each panel
"""
struct Outputs{TF}
    CF::SVector{3, TF}
    CM::SVector{3, TF}
    CDiff::Float64
    panels::Vector{PanelOutputs{TF}}
end

function Outputs(CF, CM, CDiff, panels)
    TF = promote_type(eltype(CF), eltype(CM), typeof(CDiff), eltype(eltype(panels)))
    return Outputs{TF}(CF, CM, CDiff, panels)
end

Base.eltype(::Type{Outputs{TF}}) where TF = TF
Base.eltype(::Outputs{TF}) where TF = TF

# ----------- Near Field Solution for Forces ---------------

"""
    forces_moments(panels, reference, freestream, Γ, symmetric, xhat)

Computes the forces and moments acting on the aircraft given the circulation
distribution `Γ`
"""
forces_moments

@inline function forces_moments(panels::AbstractVector{<:Horseshoe}, reference, freestream, Γ, symmetric)

    # float number type
    TF = promote_type(eltype(eltype(panels)), eltype(reference), eltype(freestream), eltype(Γ))

    # number of panels
    N = length(panels)

    # rotation vector from body to wind frame
    R = body_to_wind(freestream)

    # initialize panel outputs
    paneloutputs = Vector{PanelOutputs{TF}}(undef, N)

    # initialize body outputs
    Fb = @SVector zeros(TF, 3)  # forces
    Mb = @SVector zeros(TF, 3)  # moments

    # reference quantities
    Sref = reference.S
    bref = reference.b
    cref = reference.c

    for i = 1:N # control points

        # compute induced velocity at quarter-chord midpoints (rather than at control points)
        Vind = @SVector zeros(3)
        rmid = midpoint(panels[i])

        for j = 1:N  # vortices
            Vij = induced_velocity(rmid, panels[j], symmetric, i != j)
            Vind += Vij*Γ[j]
        end

        # add external velocity
        Vext = external_velocity(freestream, rmid, reference.rcg)
        Vi = Vind + Vext

        # forces and moments
        Δs = panel_width(panels[i])
        Fbi = RHO*Γ[i]*cross(Vi, Δs)
        Mbi = cross(rmid - reference.rcg, Fbi)

        # add panel forces and moments to body forces and moments
        Fb += Fbi
        Mb += Mbi

        # force per unit length along wing (y and z)
        ds = sqrt(Δs[2]^2 + Δs[3]^2)
        Fp = Fbi/ds

        # rotate distributed forces from body to wind frame
        Fp = R*Fp

        # assemble normalized panel outputs
        paneloutputs[i] = PanelOutputs(Γ[i]/VINF, Vi/VINF, Fp/(QINF*Sref))
    end

    # adjust body forces to account for symmetry
    if symmetric
        Fb = SVector(2*Fb[1], 0.0, 2*Fb[3])
        Mb = SVector(0.0, 2*Mb[2], 0.0)
    end

    Fw = R*Fb
    Mw = R*Mb

    # normalize body forces
    CF = Fw/(QINF*Sref)
    CM = Mw./(QINF*Sref*SVector(bref, cref, bref))

    return CF, CM, paneloutputs
end

@inline function forces_moments(panels::AbstractVector{<:Ring}, reference, freestream, Γ, symmetric)

    # float number type
    TF = promote_type(eltype(eltype(panels)), eltype(reference), eltype(freestream), eltype(Γ))

    # number of panels
    N = length(panels)

    # rotation vector from body to wind frame
    R = body_to_wind(freestream)

    # initialize panel outputs
    paneloutputs = Vector{PanelOutputs{TF}}(undef, N)

    # initialize body outputs
    Fb = @SVector zeros(TF, 3)  # forces
    Mb = @SVector zeros(TF, 3)  # moments

    # reference quantities
    Sref = reference.S
    bref = reference.b
    cref = reference.c

    for i = 1:N # control points

        # compute induced velocity at quarter-chord midpoints (rather than at control points)
        Vind = @SVector zeros(3)
        rmid = midpoint(panels[i])

        for j = 1:N  # vortices

            # note that we assume panels are grouped into chordwise strips
            # that are ordered from leading edge to trailing edge
            # e.g.
            #     1 4 7 10         1 10 7 4
            #     2 5 8 11   and   2 11 8 5   are both valid orderings
            #     3 6 9 12         3 12 9 6
            include_top = i != j
            include_bottom = i != j+1
            Vij = induced_velocity(rmid, panels[j], symmetric, include_top, include_bottom)
            Vind += Vij*Γ[j]
        end


        # add external velocity
        Vext = external_velocity(freestream, rmid, reference.rcg)
        Vi = Vind + Vext

        # extract circulation for this panel
        if i == 1 || panels[i-1].trailing
            # panel is on the leading edge
            Γi = Γ[i]
        else
            Γi = Γ[i] - Γ[i-1]
        end

        # forces and moments
        Δs = panel_width(panels[i])
        Fbi = RHO*Γi*cross(Vi, Δs)
        Mbi = cross(rmid - reference.rcg, Fbi)

        # add panel forces and moments to body forces and moments
        Fb += Fbi
        Mb += Mbi

        # force per unit length along wing (y and z)
        ds = sqrt(Δs[2]^2 + Δs[3]^2)
        Fp = Fbi/ds

        # rotate distributed forces from body to wind frame
        Fp = R*Fp

        # assemble normalized panel outputs
        paneloutputs[i] = PanelOutputs(Γi/VINF, Vi/VINF, Fp/(QINF*Sref))
    end

    # adjust body forces to account for symmetry
    if symmetric
        Fb = SVector(2*Fb[1], 0.0, 2*Fb[3])
        Mb = SVector(0.0, 2*Mb[2], 0.0)
    end

    Fw = R*Fb
    Mw = R*Mb

    # normalize body forces
    CF = Fw/(QINF*Sref)
    CM = Mw./(QINF*Sref*SVector(bref, cref, bref))

    return CF, CM, paneloutputs
end

# ----------- Far Field Solution for Forces ---------------

"""
    project_panels(panels, freestream)

Non-allocating version of `project_panels`. Overwrites `panels`.
"""
project_panels

@inline function project_panels(panels::AbstractVector{<:Horseshoe}, freestream)

    R = body_to_wind(freestream)

    panels = rotate.(panels, Ref(R))

    return panels
end

@inline function project_panels(panels::AbstractVector{<:Ring}, freestream)

    # get rotation matrix for rotating to freestream
    R = body_to_wind(freestream)

    panels = rotate.(panels, Ref(R))

    return panels
end

"""
    vortex_induced_drag(rj, Γj, ri, Γi, nhat_ds_i)

Induced drag from vortex `j` induced on panel `i`
"""
@inline function vortex_induced_drag(rj, Γj, ri, Γi, nhat_ds_i, xhat=SVector(1, 0, 0))

    rij = ri - rj
    rij -= dot(xhat, rij)*xhat  # 2D plane (no xhat-component)
    Vthetai = cross(Γj*xhat, rij) / (2*pi*norm(rij)^2)
    Vn = -dot(Vthetai, nhat_ds_i)

    Di = RHO/2.0*Γi*Vn

    return Di
end


"""
    panel_induced_drag(panel_j, Γj, panel_i, Γi, symmetric)

Induced drag from `panel_j` induced on `panel_i`
"""
panel_induced_drag

@inline function panel_induced_drag(panel_j::Horseshoe, Γj, panel_i::Horseshoe, Γi, symmetric)

    nhat_i = trefftz_plane_normal(panel_i)
    ri = midpoint(panel_i)

    rl_j = panel_j.rl
    rr_j = panel_j.rr

    Di = vortex_induced_drag(rl_j, -Γj, ri, Γi, nhat_i)
    Di += vortex_induced_drag(rr_j, Γj, ri, Γi, nhat_i)

    if symmetric && not_on_symmetry_plane(rl_j, rr_j)
        Di += vortex_induced_drag(flipy(rr_j), -Γj, ri, Γi, nhat_i)
        Di += vortex_induced_drag(flipy(rl_j), Γj, ri, Γi, nhat_i)
    end

    return Di
end

@inline function panel_induced_drag(panel_j::Ring, Γj, panel_i::Ring, Γi, symmetric)

    nhat_i = trefftz_plane_normal(panel_i)
    ri = midpoint(panel_i)

    rl_j = panel_j.rbl
    rr_j = panel_j.rbr

    Di = vortex_induced_drag(rl_j, -Γj, ri, Γi, nhat_i)
    Di += vortex_induced_drag(rr_j, Γj, ri, Γi, nhat_i)

    if symmetric && not_on_symmetry_plane(rl_j, rr_j)
        Di += vortex_induced_drag(flipy(rr_j), -Γj, ri, Γi, nhat_i)
        Di += vortex_induced_drag(flipy(rl_j), Γj, ri, Γi, nhat_i)
    end

    return Di
end

"""
    trefftz_induced_drag(panels, reference, freestream, Γ, symmetric)

Computes induced drag using the Trefftz plane (far field method).
"""
trefftz_induced_drag

@inline function trefftz_induced_drag(panels::AbstractVector{<:Horseshoe}, reference, freestream, Γ, symmetric)

    # rotate panels into wind coordinate system
    panels = project_panels(panels, freestream)

    N = length(panels)
    TF = promote_type(eltype(eltype(panels)), eltype(Γ))

    # add up drag
    Di = zero(TF)
    for j = 1:N
        for i = 1:N
            Di += panel_induced_drag(panels[j], Γ[j], panels[i], Γ[i], symmetric)
        end
    end

    # apply symmetry
    if symmetric
        Di *= 2
    end

    # normalize
    CDi = Di / (QINF*reference.S)

    return CDi, panels
end

@inline function trefftz_induced_drag(panels::AbstractVector{<:Ring}, reference, freestream, Γ, symmetric)

    # get panels that shed trailing vortices
    ip = findall((p)->p.trailing, panels)

    # rotate panels into wind coordinate system
    wake_panels = project_panels(panels[ip], freestream)

    N = length(wake_panels)
    TF = promote_type(eltype(eltype(wake_panels)), eltype(Γ))

    Di = zero(TF)
    for j = 1:N
        for i = 1:N
            Di += panel_induced_drag(wake_panels[j], Γ[ip[j]], wake_panels[i], Γ[ip[i]], symmetric)
        end
    end

    if symmetric
        Di *= 2
    end

    # normalize
    CDi = Di / (QINF*reference.S)

    return CDi, panels
end

"""
    vlm(panels, reference, freestream, symmetric)

Runs the vortex lattice method.
"""
function vlm(panels, reference, freestream, symmetric)

    Γ = circulation(panels, reference, freestream, symmetric)
    CF, CM, paneloutputs = forces_moments(panels, reference, freestream, Γ, symmetric)
    CDiff, trefftz_panels = trefftz_induced_drag(panels, reference, freestream, Γ, symmetric)

    return Outputs(CF, CM, CDiff, paneloutputs)
end
