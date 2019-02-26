#=
Author: Andrew Ning

Vortex Lattice Method

- spanwise and chordwise panels
- additional Trefftz plane analysis for induced drag
- stability derivatives
- general inflow

=#

module VortexLatticeMethod

import LinearAlgebra: norm, cross, dot

export Panel, Freestream, Reference, Outputs
export solve


include("geometry.jl")  # defines some convenience functions for generating geometry
include("sderiv.jl")  # defines stability derivative type and some associated methods


# --------- structs -------------

"""
    Panel(rl, rr, chord, theta)

Define one panel

**Arguments**
- `rl::Array{Float64}(3)`: vector position for the left side of the horseshoe vortex
- `rr::Array{Float64}(3)`: vector position for the right side of the horseshoe vortex
- `rcp::Array{Float64}(3)`: vector position for the control point
- `theta::Float64`: twist angle (radians)
"""
struct Panel{TF}
    rl::Array{TF}
    rr::Array{TF}
    rcp::Array{TF}
    theta::TF
end

# chord::TF

"""
    Freestream(alpha, beta, Omega, vother)

Define the freestream properties.

**Arguments**
- `alpha::Float64`: angle of attack (rad)
- `beta::Float64`: sideslip attack (rad)
- `Omega::Array{Float64}(3)`: rotation (p, q, r) about the c.g., normalized by Vinf
- `vother`: a function of the form: Vext = vother(r).  returns velocity at position r.  can also use nothing.
"""
struct Freestream{TF, TV}
    alpha::TF
    beta::TF
    Omega::Array{TF}
    vother::TV  
end


"""
    Reference(S, c, b, rcg)

Reference quantities.

**Arguments**
- `S::Float64`: reference area
- `c::Float64`: reference chord
- `b::Float64`: reference span
- `rcg::Array{Float64}(3)`: location of center of gravity
"""
struct Reference{TF}
    S::TF
    c::TF
    b::TF
    rcg::Array{TF}
end


"""
    Outputs()

Outputs from the VLM analysis

**Arguments**
- `CF::Float64`: force coeffients (CD, CY, CL) normalized by reference area
TODO
"""
struct Outputs{TF, TSD}
    CF::Array{TF}
    CM::Array{TF}
    dCF::TSD
    dCM::TSD
    CDiff::TF
    ymid::Array{TF}
    zmid::Array{TF}
    cf::Array{TF}
    ds::Array{TF}
    V::Array{TF}
    Gamma::Array{TF}
end


# normalized (so they don't matter, but are included just for clarity in the algorithms)
const RHO = 1.0
const VINF = 1.0


# -------------------------------

# ----------- Convenience Methods ---------------


"""
    flipy(r)

flip sign of y-component of vector (used for symmetry)
"""
function flipy(r)
    return [r[1]; -r[2]; r[3]]
end


"""
    not_on_symmetry_plane(r1, r2)

Checks whether r1 and r2 are on the symmetry plane (y = 0)
"""
function not_on_symmetry_plane(r1, r2)
    return !(isapprox(r1[2], 0.0, atol=1e-12) && isapprox(r2[2], 0.0, atol=1e-12))
end



"""
    mid_point(p::Panel)

Compute the positions of the middle of each bound vortex
"""
function mid_point(p::Panel)
    return 0.5*(p.rl + p.rr)
end



"""
rotation matrix from body axis to wind axis
"""
function body_to_wind(fs::Freestream)

    alpha = fs.alpha
    beta = fs.beta

    Rb = [cos(beta) -sin(beta) 0;
        sin(beta) cos(beta) 0.0;
        0 0 1]
    Ra = [cos(alpha) 0.0 sin(alpha);
        0 1 0;
        -sin(alpha) 0 cos(alpha)]
    R = Rb*Ra

    # stability derivatives
    dRa = SDeriv([-sin(alpha) 0.0 cos(alpha);
            0 0 0;
            -cos(alpha) 0 -sin(alpha)], 
            zeros(3, 3), zeros(3, 3), zeros(3, 3), zeros(3, 3))
    dRb = SDeriv(zeros(3, 3), 
            [-sin(beta) -cos(beta) 0;
            cos(beta) -sin(beta) 0.0;
            0 0 0.0], 
            zeros(3, 3), zeros(3, 3), zeros(3, 3))

    dR = Rb*dRa + dRb*Ra
    
    return R, dR
end


# -------------------------------


# ----------- Induced Velocity ---------------

"""
    vhat_trailing(r1, r2)

Computes the induced velocity (per unit circulation) for two vortices trailing in the +x
direction, at a control point located at r1 relative to the start of the left trailing vortex
and r2 relative to the start of the right trailing vortex.
"""
function vhat_trailing(r1, r2)
    
    nr1 = norm(r1)
    nr2 = norm(r2)
    xhat = [1.0; 0; 0]
    
    f1 = cross(r1, xhat)/(nr1 - r1[1])/nr1
    f2 = cross(r2, xhat)/(nr2 - r2[1])/nr2

    Vhat = (f1 - f2)/(4*pi)

    return Vhat
end


"""
    vhat_horseshoe(r1, r2)

Computes the induced velocity (per unit circulation) for a horseshoe vortex trailing in the +x
direction, at a control point located at r1 relative to the start of the left trailing vortex
and r2 relative to the start of the right trailing vortex.
"""
function vhat_horseshoe(r1, r2)

    # contribution from trailing vortices
    Vhat = vhat_trailing(r1, r2)
    
    # contribution from bound vortex
    nr1 = norm(r1)
    nr2 = norm(r2)
    f1 = cross(r1, r2)/(nr1*nr2 + dot(r1, r2))
    f2 = (1.0/nr1 + 1.0/nr2)
    
    Vhat += (f1*f2)/(4*pi)

    return Vhat
end


"""
    vhat(rcp, rl, rr, symmetric, include_bound)

Computes the induced velocity (per unit circulation) for a horseshoe vortex, with trailing
vortices in the +x direction.  The velocity is computed at point rcp
from a panel defined by position rl and rr.
"""
function vhat(rcp, rl, rr, symmetric, include_bound)
    
    if include_bound
        vfunc = vhat_horseshoe
    else
        vfunc = vhat_trailing
    end

    r1 = rcp - rl
    r2 = rcp - rr        
    Vhat = vfunc(r1, r2)

    if symmetric && not_on_symmetry_plane(rl, rr) # add contribution from other side
        # flip sign for y, but r1 is on left which now corresponds to rr and vice vesa

        vfunc = vhat_horseshoe  # always include bound vortex from other side
        
        r1 = rcp - flipy(rr)
        r2 = rcp - flipy(rl)
        Vhat += vfunc(r1, r2)
    end

    return Vhat
end



# -------------------------------


# ----------- AIC -----------

"""
    normal_vector(p::Panel)

Compute the normal vector for the panel
"""
function normal_vector(p::Panel)

    delta = p.rr - p.rl
    dy = delta[2]
    dz = delta[3]
    ds = sqrt(dy^2 + dz^2)
    sphi = dz/ds
    cphi = dy/ds

    nhat = [sin(p.theta); 
            -cos(p.theta)*sphi;
            cos(p.theta)*cphi]

    return nhat
end



"""
    aic(v::Vehicle, symmetric)

aerodynamic influence coefficients
"""
function aic(panels::Array{Panel, 1}, symmetric)

    # rename
    N = length(panels)

    include_bound = true  # include bound vortices
    
    AIC = zeros(N, N)
    for i = 1:N  # CP

        # normal vector body axis
        nhati = normal_vector(panels[i])

        for j = 1:N  # QC
            Vij = vhat(panels[i].rcp, panels[j].rl, panels[j].rr, symmetric, include_bound)
            AIC[i, j] = dot(Vij, nhati)
        end
    end

    return AIC
end


# -------------------------------------------


# ---------- external velocity for b.c. --------------


"""
    ext_velocity(fs::Freestream, r, rcg)

Compute the external velocity at location r
and the corresonding (partial) stability derivatives.
"""
function ext_velocity(fs::Freestream, r, rcg)

    # Freestream velocity in body coordinate
    Vext = VINF*[cos(fs.alpha)*cos(fs.beta);
            -sin(fs.beta);
            sin(fs.alpha)*cos(fs.beta)]

    Vext -= VINF*cross(fs.Omega, r - rcg)  # unnormalize

    if fs.vother != nothing
        Vext += VINF*fs.vother(r)  # unnormalize
    end

    # for stability derivatives
    dVda = VINF*[-sin(fs.alpha)*cos(fs.beta);
            0.0;
            cos(fs.alpha)*cos(fs.beta)]
    dVdb = VINF*[-cos(fs.alpha)*sin(fs.beta);
        -cos(fs.beta);
        -sin(fs.alpha)*sin(fs.beta)]
    rvec = r - rcg
    dVdp = VINF*[0.0; rvec[3]; -rvec[2]]
    dVdq = VINF*[-rvec[3]; 0.0; rvec[1]]
    dVdr = VINF*[rvec[2]; -rvec[1]; 0.0]

    dVext = SDeriv(dVda, dVdb, dVdp, dVdq, dVdr)

    return Vext, dVext
end


"""
    vn_ext(panels, ref, fs)

Compute the normal component of the external velocity along the geometry.
This forms the right hand side of the circulation linear system solve.
"""
function vn_ext(panels::Array{Panel, 1}, ref::Reference, fs::Freestream)
    
    # initialize
    N = length(panels)
    b = zeros(N)
    db = SDeriv(N)

    # iterate through panels
    for i = 1:N

        # normal vector
        nhat = normal_vector(panels[i])

        # external velocity
        Vext, dVext = ext_velocity(fs, panels[i].rcp, ref.rcg)

        # right hand side vector
        b[i] = -dot(Vext, nhat)

        # (partial) stability derivatives
        db[i] = -dot(dVext, nhat)
    end

    return b, db
end

# -------------------------------------


# ------ circulation solve ---------

"""
    circulation(panels, ref, fs, symmetric)

Solve for circulation distribution.
"""
function circulation(panels::Array{Panel, 1}, ref::Reference, fs::Freestream, symmetric)

    AIC = aic(panels, symmetric)
    b, db = vn_ext(panels, ref, fs)

    Gamma = AIC\b

    dGamma = SDeriv(AIC\db.alpha, AIC\db.beta, AIC\db.p, AIC\db.q, AIC\db.r)

    return Gamma, dGamma
end


# -----------------------------


# ----------- forces/moments --------------

"""
    forces_moments(panels, ref, fs, Gamma, symmetric)

Computes the forces and moments acting on the aircraft using the given circulation.
"""
function forces_moments(panels::Array{Panel, 1}, ref::Reference, fs::Freestream, Gamma, dGamma, symmetric)

    N = length(panels)

    # initialize
    Fb = zeros(3)  # forces
    Mb = zeros(3)  # moments
    Fpvec = zeros(3, N)  # distributed forces
    Vtotal = zeros(3, N)  # total velocities
    ds = zeros(N)  # panel size
    
    dFb = SDeriv(3)
    dMb = SDeriv(3)

    for i = 1:N  # control points

        Vindi = zeros(3)
        dVindi = SDeriv(3)

        rmid = mid_point(panels[i])  # compute induced velocity at quarter-quard midpoints (rather than at control points)

        for j = 1:N  # vortices  
            Vij = vhat(rmid, panels[j].rl, panels[j].rr, symmetric, i != j)  # include bound vortices if i != j
            Vindi += Vij*Gamma[j]
        
            dVindi += Vij*dGamma[j]
        end

        # add external velocity
        Vext, dVext = ext_velocity(fs, rmid, ref.rcg)
        Vi = Vindi + Vext
        
        dVi = dVindi + dVext
        
        # forces and moments
        Delta_s = panels[i].rr - panels[i].rl
        Fbi = RHO*Gamma[i]*cross(Vi, Delta_s)
        Mbi = cross(rmid - ref.rcg, Fbi)
        Fb += Fbi
        Mb += Mbi

        dFbi = RHO*cross(Gamma[i]*dVi + dGamma[i]*Vi, Delta_s)
        dFb += dFbi
        dMb += cross(rmid - ref.rcg, dFbi)
    
        # force per unit length along wing (y and z)
        ds[i] = sqrt(Delta_s[2]^2 + Delta_s[3]^2)
        Fpvec[:, i] = Fbi/ds[i]  #*0.5*RHO*norm(Vi)^2*panels[i].chord)  # normalize by local velocity not freestream

        # save in array
        Vtotal[:, i] = Vi
    end


    if symmetric
        Fb *= 2
        Mb *= 2
        Fb[2] = 0.0
        Mb[1] = 0.0
        Mb[3] = 0.0
        
        dFb *= 2
        dMb *= 2
        dFb[2] = 0.0
        dMb[1] = 0.0
        dMb[3] = 0.0
    end

    Rot, dRot = body_to_wind(fs)
    Fw = Rot*Fb
    Mw = Rot*Mb

    dFw = Rot*dFb + dRot*Fb
    dMw = Rot*dMb + dRot*Mb

    # rotate distributed forces from body to wind frame
    Fwpvec = zeros(3, N)
    for i = 1:N
        Fwpvec[:, i] = Rot*Fpvec[:, i]
    end

    return Fw, Mw, dFw, dMw, Fwpvec, ds, Vtotal
end

# ---------------------------------


# --------- Trefftz Plane --------

"""
project panels onto Trefftz plane (rotate into wind coordinate system)
"""
function project_panels(panels::Array{Panel, 1}, fs::Freestream)

    Rot, _ = body_to_wind(fs)

    N = length(panels)
    newpanels = Array{Panel}(undef, N)

    for i = 1:N
        rl_wind = Rot*panels[i].rl
        rr_wind = Rot*panels[i].rr
        rcp_wind = Rot*panels[i].rcp

        newpanels[i] = Panel(rl_wind, rr_wind, rcp_wind, panels[i].theta)
    end

    return newpanels
end


"""
Compute the normal vector for the panel when projected in to the Trefftz plane (including magnitude)
"""
function normal_vector_magnitude_2d(p::Panel)

    # includes magnitude
    delta = p.rr - p.rl
    dy = delta[2]
    dz = delta[3]
    
    nhat_ds = [0.0; -dz; dy]

    return nhat_ds
end


"""
Induced drag from velocity from vortex j induced onto panel i
"""
function Disubsub(rj, Gammaj, ri, Gammai, nhat_ds_i)

    rij = ri - rj
    rij[1] = 0.0  # 2D plane (no x-component)
    Vthetai = cross([Gammaj; 0.0; 0.0], rij) / (2*pi*norm(rij)^2)
    Vn_ds = -dot(Vthetai, nhat_ds_i)
    
    Di = RHO/2.0*Gammai*Vn_ds

    return Di
end


"""
Induced drag from velocity from panel j induced onto panel i
"""
function Disub(panel_j, Gamma_j, panel_i, Gamma_i, symmetric)

    nhat_ds_i = normal_vector_magnitude_2d(panel_i)
    ri = mid_point(panel_i)

    rl_j = panel_j.rl
    rr_j = panel_j.rr

    Di = Disubsub(rl_j, -Gamma_j, ri, Gamma_i, nhat_ds_i)
    Di += Disubsub(rr_j, Gamma_j, ri, Gamma_i, nhat_ds_i)

    if symmetric && not_on_symmetry_plane(rl_j, rr_j)
        Di += Disubsub(flipy(rr_j), -Gamma_j, ri, Gamma_i, nhat_ds_i)
        Di += Disubsub(flipy(rl_j), Gamma_j, ri, Gamma_i, nhat_ds_i)
    end

    return Di
end
    


"""
Far-field method to compute induced drag.
"""
function Di_trefftz(panels::Array{Panel, 1}, fs::Freestream, Gamma, symmetric)
    
    # rotate into wind coordinate system
    newpanels = project_panels(panels, fs)
    
    N = length(newpanels)
    Di = 0.0

    for j = 1:N
        for i = 1:N
            Di += Disub(newpanels[j], Gamma[j], newpanels[i], Gamma[i], symmetric)
        end
    end
    
    if symmetric
        Di *= 2
    end

    return Di
end


# -----------------------------------------


# ------------ run method --------------------

"""
    solve(panels, ref, fs, symmetric)

Run the vortex lattice method.

# Returns 
- CF: force coefficients, wind axes
- CM: moment coefficients, wind axes, about provided c.g.
- ymid, zmid: middle points along quarter chord bound vortex
- l, cl: lift distribution (Lp/(q cref)), lift coefficient distribution (Lp/(q c))
- dCF, dCM: stability derivatives, wind axes
"""
function solve(panels::Array{Panel, 1}, ref::Reference, fs::Freestream, symmetric)

    Gamma, dGamma = circulation(panels, ref, fs, symmetric)
    F, M, dF, dM, Fp, ds, Vvec = forces_moments(panels, ref, fs, Gamma, dGamma, symmetric)


    # force and moment coefficients
    qinf = 0.5*RHO*VINF^2
    Sref = ref.S
    bref = ref.b
    cref = ref.c

    # normalize forces/moments
    CF = F/(qinf*Sref)
    CM = M./(qinf*Sref*[bref; cref; bref])
    
    dCF = dF/(qinf*Sref)
    dCM = SDeriv(3)    
    dCM[1] = dM[1]/(qinf*Sref*bref)  # I'm sure there's a better way to overload this.
    dCM[2] = dM[2]/(qinf*Sref*cref)
    dCM[3] = dM[3]/(qinf*Sref*bref)


    # trefftz plane analysis for drag
    CDiff = Di_trefftz(panels, fs, Gamma, symmetric) / (qinf*Sref)  

    # # normalize p, q, r    
    dCF = SDeriv(dCF.alpha, dCF.beta, dCF.p*2*VINF/bref, dCF.q*2*VINF/cref, dCF.r*2*VINF/bref)
    dCM = SDeriv(dCM.alpha, dCM.beta, dCM.p*2*VINF/bref, dCM.q*2*VINF/cref, dCM.r*2*VINF/bref)

    # rotate p, q, r, from wind to stability axes

    # lift and cl dist
    ymid = [mid_point(p)[2] for p in panels]
    zmid = [mid_point(p)[3] for p in panels]
    # chord = [p.chord for p in panels]
    # Np = sqrt.(Fp[2, :].^2 + Fp[3, :].^2)  # normal force to panel
    # cl = Np./(qinf*chord)  # so it's not exactly cl for nonplanar wings
    # l = Np/(qinf*cref)

    # l = 2*Gamma.*Vmag./(VINF^2.*cref)
    # cl = 2*Gamma.*Vmag./(Vmag.^2.*chord)
    
    # return CF, CM, ymid, zmid, l, cl, dCF, dCM
    return Outputs(CF, CM, dCF, dCM, CDiff, ymid, zmid, Fp./(qinf*Sref), ds, Vvec/VINF, Gamma/VINF)
end


end # module
