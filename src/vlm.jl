"""
    AbstractPanel{TF}

Supertype of vortex lattice method panel types.
"""
abstract type AbstractPanel{TF} end

"""
    Horseshoe{TF}

Struct which holds the properties of a horseshoe vortex panel.

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

"""
    Ring{TF}

Struct which holds the properties of a vortex ring panel.

**Fields**
 - `rtl`: position of the top left side of the horseshoe vortex
 - `rtr`: position of the top right side of the horseshoe vortex
 - `rbl`: position of the bottom left side of the horseshoe vortex
 - `rbr`: position of the bottom right side of the horseshoe vortex
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

"""
    midpoint(panel::AbstractPanel)

Compute the bound vortex midpoint of `panel`
"""
midpoint

@inline midpoint(panel::Horseshoe) = (panel.rl + panel.rr)/2

@inline midpoint(panel::Ring) = (panel.rtl + panel.rtr)/2

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

@inline trefftz_plane_normal(panel::Ring, xhat) = cross(xhat, cross(panel.normal, xhat))

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

# -------------------------------

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
    trailing_induced_velocity(r, xhat=[1, 0, 0])

Compute the induced velocity (per unit circulation) induced by a vortex trailing
in the `xhat` direction, at a control point located at `r` relative to the start
of the trailing vortex
"""
@inline function trailing_induced_velocity(r, xhat=SVector(1, 0, 0))

    nr = norm(r)

    Vhat = cross(r, xhat)/(nr - dot(r, xhat))/nr/(4*pi)

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

"""
    horseshoe_induced_velocity(rcp, rl, rr, xhat, symmetric, include_bound)

Compute the induced velocity (per unit circulation) for a horseshoe vortex, with
trailing vortices in the `xhat` direction.  The velocity is computed at point `rcp`
from a panel defined by position `rl` and `rr`.
"""
@inline function horseshoe_induced_velocity(rcp, rl, rr, symmetric, include_bound)

    r1 = rcp - rl
    r2 = rcp - rr

    # left trailing vortex
    Vhat = trailing_induced_velocity(r1)

    # right trailing vortex
    Vhat -= trailing_induced_velocity(r2)

    # bound vortex
    if include_bound
        Vhat += bound_induced_velocity(r1, r2)
    end

    # add contribution from other side (if applicable)
    if symmetric && not_on_symmetry_plane(rl, rr)
        # flip sign for y, but r1 is on left which now corresponds to rr and vice vesa
        rl = flipy(rr)
        rr = flipy(rl)

        symmetric = false
        include_bound = true
        Vhat += horseshoe_induced_velocity(rcp, rl, rr, symmetric, include_bound)
    end

    return Vhat
end

"""
    ring_induced_velocity(rcp, rl, rr, symmetric, include_bound, trailing, xhat)

Compute the induced velocity (per unit circulation) for a ring vortex. The velocity
is computed at point `rcp` from a panel defined by positions `rtl`, `rtr`, `rbl`,
and `rbr`.
"""
@inline function horseshoe_induced_velocity(rcp, rtl, rtr, rbl, rbr, symmetric,
    include_top_bound, include_bottom_bound, trailing, xhat)

    # define position of control point relative to each panel corner
    r1 = rcp - rtl
    r2 = rcp - rtr
    r3 = rcp - rbr
    r4 = rcp - rbl

    # left bound vortex
    Vhat = bound_induced_velocity(r4, r1)

    # right bound vortex
    Vhat += bound_induced_velocity(r2, r3)

    # top bound vortex
    if include_top_bound
        Vhat += bound_induced_velocity(r1, r2)
    end

    if trailing
        # left trailing vortex
        Vhat += trailing_induced_velocity(r4, xhat)
        # right trailing vortex
        Vhat -= trailing_induced_velocity(r3, xhat)
    else
        if include_bottom_bound
            Vhat += bound_induced_velocity(r3, r4)
        end
    end

    # add contribution from other side (if applicable)
    if symmetric && not_on_symmetry_plane(rl, rr)
        # flip sign for y, but now left is right and right is left.
        rtl = flipy(rtr)
        rtr = flipy(rtl)
        rbl = flipy(rbr)
        rbr = flipy(rbl)

        symmetric = false
        include_top_bound = true
        include_bottom_bound = true
        Vhat += horseshoe_induced_velocity(rcp, rtl, rtr, rbl, rbr, symmetric,
            include_top_bound, include_bottom_bound, trailing, xhat)
    end

    return Vhat
end

# ----------- left hand side (AIC) -----------

"""
    influence_coefficients(panels::AbstractVector{<:Horseshoe}, symmetric)
    influence_coefficients(panels::AbstractMatrix{<:Ring}, symmetric, xhat)
    influence_coefficients(panels::AbstractVector{<:AbstractMatrix{<:Ring}}, symmetric, xhat)

Construct the aerodynamic influence coefficient matrix.
"""
influence_coefficients

@inline function influence_coefficients(panels::AbstractVector{<:Horseshoe}, symmetric)

    N = length(panels)
    TF = eltype(eltype(panels))
    AIC = Matrix{TF}(undef, N, N)

    return influence_coefficients!(AIC, panels, symmetric)
end

@inline function influence_coefficients(panels::AbstractMatrix{<:Ring}, symmetric, xhat)

    N = length(panels)
    TF = eltype(eltype(panels))
    AIC = Matrix{TF}(undef, N, N)

    return influence_coefficients!(AIC, panels, symmetric, xhat)
end

@inline function influence_coefficients(panels::AbstractVector{AbstractMatrix{<:Ring}}, symmetric, xhat)

    N = prod(length.(panels))
    TF = eltype(eltype(eltype(panels)))
    AIC = Matrix{TF}(undef, N, N)

    return influence_coefficients!(AIC, panels, symmetric)
end

"""
    influence_coefficients!(AIC, panels::AbstractVector{<:Horseshoe}, symmetric)
    influence_coefficients!(AIC, panels::AbstractMatrix{<:Ring}, symmetric, xhat)
    influence_coefficients!(AIC, panels::AbstractVector{<:AbstractMatrix{<:Ring}}, symmetric, xhat)
    influence_coefficients!(AIC, receiving_panels::AbstractMatrix{<:Ring}, sending_panels::AbstractMatrix{<:Ring}, symmetric, xhat)

Pre-allocated version of `influence_coefficients`
"""
influence_coefficients!

@inline function influence_coefficients!(AIC, panels::AbstractVector{<:Horseshoe}, symmetric)

    bound_vortices = true

    N = length(panels)

    # loop over control points
    for i = 1:N

        # normal vector body axis
        nhat = normal(panels[i])

        # loop over bound vortices
        for j = 1:N
            V = horseshoe_induced_velocity(panels[i].rcp, panels[j].rl, panels[j].rr, symmetric, bound_vortices)
            AIC[i, j] = dot(V, nhat)
        end
    end

    return AIC
end

@inline function influence_coefficients!(AIC, panels::AbstractMatrix{<:Ring}, symmetric, xhat)

    AIC = influence_coefficients!(AIC, panels, panels, symmetric, xhat)

    return AIC
end

@inline function influence_coefficients!(AIC, panels::AbstractVector{AbstractMatrix{<:Ring}}, symmetric, xhat)

    N = size(AIC, 1)
    nsurf = length(panels)

    # initialize panel counters
    ip = 0
    jp = 0

    # loop over receiving surfaces
    for i = 1:nsurf
        # number of receiving panels
        npi = length(panels[i])
        # corresponding range in AIC matrix
        iAIC = ip+1:ip+npi
        # loop over sending surfaces
        for j = 1:nsurf
            # number of sending panels
            npj = length(panels[j])
            # corresponding range in AIC matrix
            jAIC = jp+1:jp+npj
            # construct view of relevant portion of AIC matrix
            vAIC = view(AIC, iAIC, jAIC)
            # fill in entries in AIC matrix
            influence_coefficients!(vAIC, panels[i], panels[j], symmetric, xhat)
            # increment sending panel counter
            jp += npj
        end
        # increment receiving panel counter
        ip += npi
    end

    return AIC
end


@inline function influence_coefficients!(AIC, receiving_panels::AbstractMatrix{<:Ring}, sending_panels::AbstractMatrix{<:Ring}, symmetric, xhat)

    # get number of panels in i and j directions
    ni_r = size(receiving_panels, 1)
    nj_r = size(receiving_panels, 2)
    ni_s = size(sending_panels, 1)
    nj_s = size(sending_panels, 2)

    # reshape AIC matrix for easier indexing
    rAIC = reshape(AIC, ni_r, nj_r, ni_s, nj_s)

    for j_r = 1:nj_r
        for i_r = 1:ni_r
            # control point location
            rcp = receiving_panels[i_r, j_r].rcp

            # influence coefficients from bound vortices at the leading edge
            for j_s = 1:nj_s
                rl = sending_panels[1, j_s].rl
                rr = sending_panels[1, j_s].rr

                r1 = rcp - rl
                r2 = rcp - rr
                Vhat = bound_induced_velocity(r1, r2, xhat)
                rAIC[i_r, j_r, i_s, j_s] = Vhat

                if symmetric && not_on_symmetry_plane(rl, rr)
                    # flip sign for y, but r1 is on left which now corresponds to rr and vice vesa
                    r1 = rcp - flipy(rr)
                    r2 = rcp - flipy(rl)

                    Vhat += bound_induced_velocity(r1, r2)
                end
            end

            # influence coefficients from bound vortices at
            for i_s = 2:ni_s-1
                # bottom bound vortex of current (sending) panel
                r1 = receiving_panels[i_r, j_r].rcp - sending_panels[i_s+1, j_s].rl
                r2 = receiving_panels[i_r, j_r].rcp - sending_panels[i_s+1, j_s].rr
                Vhat = bound_induced_velocity(r1, r2, xhat)
                rAIC[i_r, j_r, i_s, j_s] -= Vhat
                # top bound vortex of next (sending) panel
                rAIC[i_r, j_r, i_s+1, j_s] = Vhat

            end

            for j_s = 2:nj_s-1
                # top bound vortex of current (sending) panel

                for i_s = 2:ni_s-1
                    # right/bound vortex of current (sending) panel
                    r1 = receiving_panels[i_r, j_r].rcp - sending_panels[i_s, j_s].rr
                    r2 = receiving_panels[i_r, j_r].rcp - sending_panels[i_s+1, j_s].rr
                    Vhat = bound_induced_velocity(r1, r2, xhat)
                    rAIC[i_r, j_r, i_s, j_s] -= Vhat
                    # left bound vortex of next (sending) panel
                    rAIC[i_r, j_r, i_s, j_s+1] = Vhat
                    # bottom bound vortex of current (sending) panel
                    r1 = receiving_panels[i_r, j_r].rcp - sending_panels[i_s+1, j_s].rl
                    r2 = receiving_panels[i_r, j_r].rcp - sending_panels[i_s+1, j_s].rr
                    Vhat = bound_induced_velocity(r1, r2, xhat)
                    rAIC[i_r, j_r, i_s, j_s] -= Vhat
                    # top bound vortex of next (sending) panel
                    rAIC[i_r, j_r, i_s+1, j_s] = Vhat
                end


            end

        end
    end


    # loop over control points
    for i = 1:N

        # normal vector body axis
        nhat = normal(panels[i])

        # loop over bound vortices
        for j = 1:N
            V = horseshoe_induced_velocity(panels[i].rcp, panels[j].rl, panels[j].rr, symmetric, bound_vortices)
            AIC[i, j] = dot(V, nhat)
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
        Vext += VINF*freestream.additional_velocity(r)  # unnormalize
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
