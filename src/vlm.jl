# ----------- Structs ---------------

"""
    Panel(rl, rr, chord, theta)

Construct and return an object of type `Panel`

**Arguments**
- `rl`: position of the left side of the horseshoe vortex
- `rr`: position of the right side of the horseshoe vortex
- `rcp`: position of the panel control point
- `theta`: twist angle of the panel (radians)
"""
struct Panel{TF}
    rl::SVector{3, TF}
    rr::SVector{3, TF}
    rcp::SVector{3, TF}
    theta::TF
end

function Panel(rl, rr, rcp, theta)
    TF = promote_type(eltype(rl), eltype(rr), eltype(rcp))
    return Panel{TF}(SVector{3}(rl), SVector{3}(rr), SVector{3}(rcp), theta)
end

Base.eltype(::Type{Panel{TF}}) where TF = TF
Base.eltype(::Panel{TF}) where TF = TF

"""
    midpoint(panel)

Compute the bound vortex midpoint of `panel`
"""
@inline midpoint(panel) = (panel.rl + panel.rr)/2

"""
    normal(panel)

Compute the normal vector of `panel`
"""
@inline function normal(panel)

    delta = panel.rr - panel.rl
    dy = delta[2]
    dz = delta[3]
    ds = sqrt(dy^2 + dz^2)
    sphi = dz/ds
    cphi = dy/ds

    ct, st = cos(panel.theta), sin(panel.theta)
    nhat = SVector(st, -ct*sphi, ct*cphi)

    return nhat
end

"""
    trefftz_plane_normal(panel)

Compute the normal vector for `panel` when projected onto the Trefftz plane (including magnitude)
"""
@inline function trefftz_plane_normal(panel)

    # includes magnitude
    delta = panel.rr - panel.rl
    dy = delta[2]
    dz = delta[3]

    nhat = SVector(0.0, -dz, dy)

    return nhat
end

"""
    translate!(panels, r)

Translates the panels contained in `panels` the distance specified by vector `r`
"""
function translate!(panels, r)

    TF = eltype(eltype(panels))

    npanels = length(panels)
    for i = 1:npanels
        rl = panels[i].rl + r
        rr = panels[i].rr + r
        rcp = panels[i].rcp + r
        theta = panels[i].theta
        panels[i] = Panel{TF}(rl, rr, rcp, theta)
    end

    return panels
end

# -------------------------------

"""
    Freestream(alpha, beta, Omega, additional_velocity=nothing)

Define the freestream properties.

**Arguments**
- `alpha`: angle of attack (rad)
- `beta`: sideslip angle (rad)
- `Omega`: rotation vector (p, q, r) of the body frame about the center of
    gravity, normalized by Vinf
- `additional_velocity`: a function of the form: V = additional_velocity(r) which returns
    the additional velocity `V` (normalized by the freestream velocity) at
    position `r`.  Defaults to `nothing`.
"""
struct Freestream{TF, TV}
    alpha::TF
    beta::TF
    Omega::SVector{3, TF}
    additional_velocity::TV
end

function Freestream(alpha, beta, Omega, additional_velocity = nothing)
    TF = promote_type(typeof(alpha), typeof(beta), eltype(Omega))
    return Freestream{TF, typeof(additional_velocity)}(alpha, beta, Omega, additional_velocity)
end

Base.eltype(::Type{Freestream{TF, TV}}) where {TF,TV} = TF
Base.eltype(::Freestream{TF, TV}) where {TF,TV} = TF

"""
    body_to_wind(fs::Freestream)

Construct a rotation matrix from the body axis to the wind axis
"""
@inline function body_to_wind(fs::Freestream)

    alpha = fs.alpha
    beta = fs.beta

    cb, sb = cos(beta), sin(beta)
    ca, sa = cos(alpha), sin(alpha)

    Rb = @SMatrix [cb -sb 0; sb cb 0; 0 0 1]
    Ra = @SMatrix [ca 0 sa; 0 1 0; -sa 0 ca]
    R = Rb*Ra

    return R
end

# -------------------------------

"""
    Reference(S, c, b, rcg)

Reference quantities.

**Arguments**
- `S`: reference area
- `c`: reference chord
- `b`: reference span
- `rcg`: location of center of gravity
"""
struct Reference{TF}
    S::TF
    c::TF
    b::TF
    rcg::SVector{3, TF}
end

function Reference(S, c, b, rcg)
    TF = promote_type(typeof(S), typeof(c), typeof(b), eltype(rcg))
    return Reference{TF}(S, c, b, rcg)
end

Base.eltype(::Type{Reference{TF}}) where TF = TF
Base.eltype(::Reference{TF}) where TF = TF

# -------------------------------

"""
    PanelOutputs

Outputs from the VLM analysis for a specific panel
"""
struct PanelOutputs{TF}
    ymid::TF
    zmid::TF
    cf::SVector{3, TF}
    ds::TF
    V::SVector{3, TF}
    Gamma::TF
end

function PanelOutputs(ymid, zmid, cf, ds, V, Gamma)
    TF = promote_type(typeof(ymid), typeof(zmid), eltype(cf), typeof(ds), eltype(V), typeof(Gamma))
    return PanelOutputs{TF}(ymid, zmid, cf, ds, V, Gamma)
end

Base.eltype(::Type{PanelOutputs{TF}}) where TF = TF
Base.eltype(::PanelOutputs{TF}) where TF = TF

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

"""
    StabilityDerivatives

Holds stability derivatives

# Fields:
 - `alpha`: angle of attack
 - `beta`: sideslip
 - `p`: roll rate
 - `q`: pitch rate
 - `r`: yaw rate
"""
struct StabilityDerivatives{TF}
    alpha::TF
    beta::TF
    p::TF
    q::TF
    r::TF
end

function StabilityDerivatives(alpha, beta, p, q, r)
    TF = promote_type(typeof(alpha), typeof(beta), typeof(p), typeof(q), typeof(r))
    return StabilityDerivatives{TF}(alpha, beta, p, q, r)
end

Base.eltype(::Type{StabilityDerivatives{TF}}) where TF = TF
Base.eltype(::StabilityDerivatives{TF}) where TF = TF

# ----------- Convenience Methods ---------------

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

# ----------- Induced Velocity ---------------

"""
    trailing_induced_velocity(r1, r2)

Compute the induced velocity (per unit circulation) for two vortices trailing in
the +x direction, at a control point located at `r1` relative to the start of the
left trailing vortex and `r2` relative to the start of the right trailing vortex.
"""
@inline function trailing_induced_velocity(r1, r2)

    nr1 = norm(r1)
    nr2 = norm(r2)
    xhat = SVector(1, 0, 0)

    f1 = cross(r1, xhat)/(nr1 - r1[1])/nr1
    f2 = cross(r2, xhat)/(nr2 - r2[1])/nr2

    Vhat = (f1 - f2)/(4*pi)

    return Vhat
end

"""
    horseshoe_induced_velocity(r1, r2)

Compute the induced velocity (per unit circulation) for a horseshoe vortex
trailing in the +x direction, at a control point located at `r1` relative to the
start of the left trailing vortex and `r2` relative to the start of the right
trailing vortex.
"""
@inline function horseshoe_induced_velocity(r1, r2)

    # contribution from trailing vortices
    Vhat = trailing_induced_velocity(r1, r2)

    # contribution from bound vortex
    nr1 = norm(r1)
    nr2 = norm(r2)
    f1 = cross(r1, r2)/(nr1*nr2 + dot(r1, r2))
    f2 = (1/nr1 + 1/nr2)

    Vhat += (f1*f2)/(4*pi)

    return Vhat
end

"""
    induced_velocity(rcp, rl, rr, symmetric, include_bound)

Compute the induced velocity (per unit circulation) for a horseshoe vortex, with
trailing vortices in the +x direction.  The velocity is computed at point `rcp`
from a panel defined by position `rl` and `rr`.
"""
@inline function induced_velocity(rcp, rl, rr, symmetric, include_bound)

    if include_bound
        vfunc = horseshoe_induced_velocity
    else
        vfunc = trailing_induced_velocity
    end

    r1 = rcp - rl
    r2 = rcp - rr
    Vhat = vfunc(r1, r2)

    # add contribution from other side (if applicable)
    if symmetric && not_on_symmetry_plane(rl, rr)
        # flip sign for y, but r1 is on left which now corresponds to rr and vice vesa

        vfunc = horseshoe_induced_velocity  # always include bound vortex from other side

        r1 = rcp - flipy(rr)
        r2 = rcp - flipy(rl)
        Vhat += vfunc(r1, r2)
    end

    return Vhat
end

# ----------- left hand side (AIC) -----------

"""
    influence_coefficients(panels, symmetric)

Construct the aerodynamic influence coefficient matrix
"""
@inline function influence_coefficients(panels, symmetric)

    N = length(panels)
    TF = eltype(eltype(panels))
    AIC = Matrix{TF}(undef, N, N)

    influence_coefficients!(AIC, panels, symmetric)

    return influence_coefficients!(AIC, panels, symmetric)
end

"""
    influence_coefficients!(AIC, panels, symmetric)

Pre-allocated version of `influence_coefficients`
"""
@inline function influence_coefficients!(AIC, panels, symmetric)

    bound_vortices = true

    N = length(panels)

    # loop over control points
    for i = 1:N  # CP

        # normal vector body axis
        nhati = normal(panels[i])

        # loop over bound vortices
        for j = 1:N
            Vij = induced_velocity(panels[i].rcp, panels[j].rl, panels[j].rr, symmetric, bound_vortices)
            AIC[i, j] = dot(Vij, nhati)
        end
    end

    return AIC
end

# ---------- right hand side (boundary conditions) --------------

"""
    external_velocity(freestream, r, rcg)

Compute the external velocity at location `r`
"""
@inline function external_velocity(freestream, r, rcg)

    # Freestream velocity in body coordinate
    ca, sa = cos(freestream.alpha), sin(freestream.alpha)
    cb, sb = cos(freestream.beta), sin(freestream.beta)

    Vext = VINF*SVector(ca*cb, -sb, sa*cb)

    # add velocity due to body rotation
    Vext -= VINF*cross(freestream.Omega, r - rcg)  # unnormalize

    # add contribution due to additional velocity field
    if !isnothing(freestream.additional_velocity)
        Vext += VINF*fs.additional_velocity(r)  # unnormalize
    end

    return Vext
end

"""
    normal_velocity(panels, reference, freestream)

Compute the normal component of the external velocity along the geometry.
This forms the right hand side of the circulation linear system solve.
"""
@inline function normal_velocity(panels, reference, freestream)

    N = length(panels)
    TF = promote_type(eltype(eltype(panels)), eltype(reference), eltype(freestream))
    b = Vector{TF}(undef, N)

    return normal_velocity!(b, panels, reference, freestream)
end

"""
    normal_velocity!(b, panels, reference, freestream)

Non-allocating version of `normal_velocity`
"""
@inline function normal_velocity!(b, panels, reference, freestream)

    N = length(panels)

    # iterate through panels
    for i = 1:N

        # normal vector
        nhat = normal(panels[i])

        # external velocity
        Vext = external_velocity(freestream, panels[i].rcp, reference.rcg)

        # right hand side vector
        b[i] = -dot(Vext, nhat)

    end

    return b
end

# ------ circulation solve ---------

"""
    circulation(panels, reference, freestream, symmetric)

Solve for the circulation distribution.
"""
@inline function circulation(panels, reference, freestream, symmetric)
    AIC = influence_coefficients(panels, symmetric)
    b = normal_velocity(panels, reference, freestream)
    return AIC\b
end

"""
    circulation!(Γ, AIC, b, panels, reference, freestream, symmetric)

Pre-allocated version of `circulation`
"""
@inline function circulation!(Γ, AIC, b, panels, reference, freestream, symmetric)
    AIC = influence_coefficients!(AIC, panels, symmetric)
    b = vn!(b, panels, reference, freestream)
    return ldiv!(Γ, factorize(AIC), b)
end

# ----------- forces/moments --------------

"""
    forces_moments(panels, reference, freestream, Γ, symmetric)

Computes the forces and moments acting on the aircraft given the circulation
distribution `Γ`
"""
@inline function forces_moments(panels, reference, freestream, Γ, symmetric)

    N = length(panels)
    TF = promote_type(eltype(eltype(panels)), eltype(reference), eltype(freestream))

    # reference quantities
    qinf = 0.5*RHO*VINF^2
    Sref = reference.S
    bref = reference.b
    cref = reference.c

    # rotation vector from body to wind frame
    Rot = body_to_wind(freestream)

    # initialize body outputs
    Fb = @SVector zeros(TF, 3)  # forces
    Mb = @SVector zeros(TF, 3)  # moments

    # initialize panel outputs
    paneloutputs = Vector{PanelOutputs{TF}}(undef, N)

    for i = 1:N # control points

        # compute induced velocity at quarter-chord midpoints (rather than at control points)
        Vind = @SVector zeros(3)
        rmid = midpoint(panels[i])
        for j = 1:N  # vortices
            Vij = induced_velocity(rmid, panels[j].rl, panels[j].rr, symmetric, i != j)  # include bound vortices if i != j
            Vind += Vij*Γ[j]
        end

        # add external velocity
        Vext = external_velocity(freestream, rmid, reference.rcg)
        Vi = Vind + Vext

        # forces and moments
        Δs = panels[i].rr - panels[i].rl
        Fbi = RHO*Γ[i]*cross(Vi, Δs)
        Mbi = cross(rmid - reference.rcg, Fbi)

        # add panel forces and moments to body forces and moments
        Fb += Fbi
        Mb += Mbi

        # force per unit length along wing (y and z)
        ds = sqrt(Δs[2]^2 + Δs[3]^2)
        Fp = Fbi/ds  #*0.5*RHO*norm(Vi)^2*panels[i].chord)  # normalize by local velocity not freestream

        # rotate distributed forces from body to wind frame
        Fp = Rot*Fp

        # assemble normalized panel outputs
        paneloutputs[i] = PanelOutputs(rmid[2], rmid[3], Vi/VINF, Γ[i]/VINF, Fp/(qinf*Sref), ds)
    end

    # adjust body forces to account for symmetry
    if symmetric
        Fb = SVector(2*Fb[1], 0.0, 2*Fb[3])
        Mb = SVector(0.0, 2*Mb[2], 0.0)
    end

    Fw = Rot*Fb
    Mw = Rot*Mb

    # normalize body forces
    CF = Fw/(qinf*Sref)
    CM = Mw./(qinf*Sref*SVector(bref, cref, bref))

    return CF, CM, paneloutputs
end

# --------- Trefftz Plane --------

"""
    project_panels(panels, freestream)

Project panels onto Trefftz plane (rotate into wind coordinate system)
"""
@inline project_panels(panels, freestream) = project_panels!(deepcopy(panels), freestream)

"""
    project_panels!(panels, freestream)

Non-allocating version of `project_panels`. Overwrites `panels`.
"""
@inline function project_panels!(panels, freestream)

    Rot = body_to_wind(freestream)

    N = length(panels)

    for i = 1:N
        rl_wind = Rot*panels[i].rl
        rr_wind = Rot*panels[i].rr
        rcp_wind = Rot*panels[i].rcp

        panels[i] = Panel(rl_wind, rr_wind, rcp_wind, panels[i].theta)
    end

    return panels
end

"""
    vortex_induced_drag(rj, Γj, ri, Γi, nhat_ds_i)

Induced drag from vortex `j` induced on panel `i`
"""
@inline function vortex_induced_drag(rj, Γj, ri, Γi, nhat_i)

    rij = ri - rj
    rij = SVector(0.0, rij[2], rij[3])  # 2D plane (no x-component)
    Vthetai = cross(SVector(Γj, 0.0, 0.0), rij) / (2*pi*norm(rij)^2)
    Vn = -dot(Vthetai, nhat_i)

    Di = RHO/2.0*Γi*Vn

    return Di
end


"""
    panel_induced_drag(panel_j, Γj, panel_i, Γi, symmetric)

Induced drag from `panel_j` induced on `panel_i`
"""
@inline function panel_induced_drag(panel_j, Γj, panel_i, Γi, symmetric)

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

"""
    trefftz_induced_drag(panels, reference, freestream, Γ, symmetric)

Compute induced drag using the Trefftz plane (far field method).
"""
@inline function trefftz_induced_drag(panels, reference, freestream, Γ, symmetric)

    return trefftz_induced_drag!(copy(panels), reference, freestream, Γ, symmetric)
end

"""
    trefftz_induced_drag!(panels, freestream, Γ, symmetric)

Pre-allocated version of trefftz_induced_drag!.  Overwrites `panels` with the
panels in the Trefftz plane.
"""
@inline function trefftz_induced_drag!(panels, reference, freestream, Γ, symmetric)

    # rotate panels into wind coordinate system
    project_panels!(panels, freestream)

    N = length(panels)
    TF = promote_type(eltype(eltype(panels)), eltype(Γ))

    Di = zero(TF)
    for j = 1:N
        for i = 1:N
            Di += panel_induced_drag(panels[j], Γ[j], panels[i], Γ[i], symmetric)
        end
    end

    if symmetric
        Di *= 2
    end

    # normalize
    qinf = 0.5*RHO*VINF^2
    Sref = reference.S
    CDi = Di / (qinf*Sref)

    return CDi, panels
end

# ------------ run method --------------------

"""
    vlm(panels, reference, freestream, symmetric)

Runs the vortex lattice method.
"""
function vlm(panels, reference, freestream, symmetric)

    # make sure panels is of concrete type
    TF = promote_type(eltype.(panels)...)
    panels = Vector{Panel{TF}}(panels)

    # reference quantities
    qinf = 0.5*RHO*VINF^2
    Sref = reference.S
    bref = reference.b
    cref = reference.c

    Γ = circulation(panels, reference, freestream, symmetric)
    CF, CM, paneloutputs = forces_moments(panels, reference, freestream, Γ, symmetric)
    CDiff, trefftz_panels = trefftz_induced_drag(panels, reference, freestream, Γ, symmetric)

    return Outputs(CF, CM, CDiff, paneloutputs)
end

function stability_analysis(panels, reference, freestream, symmetric)

    # make sure panels is of concrete type
    TF = promote_type(eltype.(panels)...)
    panels = Vector{Panel{TF}}(panels)

    x = vcat(freestream.alpha, freestream.beta, freestream.Omega)

    f = function(x)
        freestream = Freestream(x[1], x[2], SVector(x[3], x[4], x[5]))
        Γ = circulation(panels, reference, freestream, symmetric)
        CF, CM, _ = forces_moments(panels, reference, freestream, Γ, symmetric)
        return vcat(CF, CM)
    end

    dfdx = ForwardDiff.jacobian(f, x)

    iCF = SVector(1,2,3)
    CF_α = dfdx[iCF,1]
    CF_β = dfdx[iCF,2]
    CF_p = dfdx[iCF,3]*2*VINF/reference.b
    CF_q = dfdx[iCF,4]*2*VINF/reference.c
    CF_r = dfdx[iCF,5]*2*VINF/reference.b
    dCF = StabilityDerivatives(CF_α, CF_β, CF_p, CF_q, CF_r)

    iCM = SVector(4,5,6)
    CM_α = dfdx[iCM,1]
    CM_β = dfdx[iCM,2]
    CM_p = dfdx[iCM,3]*2*VINF/reference.b
    CM_q = dfdx[iCM,4]*2*VINF/reference.c
    CM_r = dfdx[iCM,5]*2*VINF/reference.b
    dCM = StabilityDerivatives(CM_α, CM_β, CM_p, CM_q, CM_r)

    return dCF, dCM

end
