module VLM

struct Panel
    rl
    rr    
    chord
    theta
end

struct Freestream
    Vinf
    alpha
    beta
    Omega  # length: 3 (p, q, r)
    vother  # Vother = vother(x, y, z).  can also use nothing
end

struct Reference
    S
    c
    b
    rcg
end


include("geometry.jl")  # defines some convenience functions for generating geometry
include("sderiv.jl")  # defines stability derivative type and some associated methods


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
end


"""
    flipy(r)

flip sign of y-component of vector (used for symmetry)
"""
function flipy(r)
    return [r[1]; -r[2]; r[3]]
end


function not_on_symmetry_plane(r1, r2)
    return !(isapprox(r1[2], 0.0, atol=1e-12) && isapprox(r2[2], 0.0, atol=1e-12))
end

function not_on_symmetry_plane(p::Panel)
    return not_on_symmetry_plane(p.rl, p.rr)
end

"""
    vhat(rcp, rl, rr, symmetric, include_bound)

Computes the induced velocity (per unit circulation) for a horseshoe vortex, with trailing
vortices in the +x direction.  The velocity is computed at control point rcp
with from a panel defined by positino rl and rr
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
        
        r1 = rcp - flipy(rr)
        r2 = rcp - flipy(rl)
        Vhat += vfunc(r1, r2)
    end

    return Vhat
end


"""
    mid_point(p::Panel)

Compute the positions of the middle of each bound vortex
"""
function mid_point(p::Panel)
    return 0.5*(p.rl + p.rr)
end


"""
    control_point(p::Panel)

Compute the positions of the control points for each panel
"""
function control_point(p::Panel)

    rcp = mid_point(p)
    rcp[1] += p.chord/2.0  # add 3/4 chord point

    return rcp
end


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
    normal_vector_magnitude_2d(p::Panel)

Compute the normal vector for the panel when projected in to the Trefftz plane (including magnitude)
"""
function normal_vector_magnitude_2d(p::Panel)
    # includes magnitude

    delta = p.rr - p.rl
    dy = delta[2]
    dz = delta[3]
    
    nhat = [0.0; -dz; dy]

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
        nhat = normal_vector(panels[i])
        rcp = control_point(panels[i])

        for j = 1:N  # QC
            Vij = vhat(rcp, panels[j].rl, panels[j].rr, symmetric, include_bound)
            AIC[i, j] = dot(Vij, nhat)
        end
    end

    return AIC
end

"""
    ext_velocity(fs::Freestream, x, y, z)

Compute the external velocity at location r
and the corresonding (partial) stability derivatives
"""
function ext_velocity(fs::Freestream, r, rcg)

    # Freestream velocity in body coordinate
    Vinf = fs.Vinf*[cos(fs.alpha)*cos(fs.beta);
        -sin(fs.beta);
        sin(fs.alpha)*cos(fs.beta)]

    Vext = Vinf - cross(fs.Omega, r - rcg)

    if fs.vother != nothing
        Vext += fs.vother(r)
    end

    # for stability derivatives
    dVda = fs.Vinf*[-sin(fs.alpha)*cos(fs.beta);
        0.0;
        cos(fs.alpha)*cos(fs.beta)]
    dVdb = fs.Vinf*[-cos(fs.alpha)*sin(fs.beta);
        -cos(fs.beta);
        -sin(fs.alpha)*sin(fs.beta)]
    rvec = r - rcg
    dVdp = [0.0; rvec[3]; -rvec[2]]
    dVdq = [-rvec[3]; 0.0; rvec[1]]
    dVdr = [rvec[2]; -rvec[1]; 0.0]

    dVext = SDeriv(dVda, dVdb, dVdp, dVdq, dVdr)

    return Vext, dVext
end


"""
    vehicle::Vehicle, fs::Freestream

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
        Vext, dVext = ext_velocity(fs, control_point(panels[i]), ref.rcg)

        # right hand side vector
        b[i] = -dot(Vext, nhat)

        # (partial) stability derivatives
        db[i] = -dot(dVext, nhat)
    end

    return b, db
end


"""
    circulation(vehicle::Vehicle, fs::Freestream, symmetric)

Solve for circulation distribution.
"""
function circulation(panels::Array{Panel, 1}, ref::Reference, fs::Freestream, symmetric)

    AIC = aic(panels, symmetric)
    b, db = vn_ext(panels, ref, fs)

    Gamma = AIC\b

    dGamma = SDeriv(AIC\db.alpha, AIC\db.beta, AIC\db.p, AIC\db.q, AIC\db.r)

    return Gamma, dGamma
end


"""
    forces_moments(geom::Geometry, fs::Freestream, Gamma, dGamma, symmetric)

Computes the forces and moments acting on the aircraft using the given circulation.
"""
function forces_moments(panels::Array{Panel, 1}, ref::Reference, fs::Freestream, Gamma, dGamma, symmetric)

    rho = 1.0  # cancels out from normalization
    N = length(panels)

    # initialize
    Fb = zeros(3)  # forces
    Mb = zeros(3)  # moments
    Fpvec = zeros(3, N)  # distributed forces
    Vvec = zeros(3, N)  # distributed velocity
    
    dFb = SDeriv(3)
    dMb = SDeriv(3)

    include_bound = false  # bound vortices don't contribute

    for i = 1:N  # control points

        Vindi = 0.0
        dVindi = SDeriv(0.0, 0.0, 0.0, 0.0, 0.0)

        rmid = mid_point(panels[i])  # compute induced velocity at quarter-quard midpoints (rather than at control points)

        for j = 1:N  # vortices  
            Vij = vhat(rmid, panels[j].rl, panels[j].rr, symmetric, include_bound)
            Vindi += Vij*Gamma[j]
        
            dVindi += Vij*dGamma[j]
        end

        # add external velocity
        Vext, dVext = ext_velocity(fs, rmid, ref.rcg)
        Vi = Vindi + Vext
        dVi = dVindi + dVext
        Vvec[:, i] = Vi  # save velocity distribution in a vector
        
        # forces and moments
        Delta_s = panels[i].rr - panels[i].rl
        Fbi = rho*Gamma[i]*cross(Vi, Delta_s)
        Fb += Fbi
        Mb += cross(rmid - ref.rcg, Fbi)

        dFbi = rho*cross(Gamma[i]*dVi + dGamma[i]*Vi, Delta_s)
        dFb += dFbi
        dMb += cross(rmid - ref.rcg, dFbi)
    
        # force per unit length
        Fpvec[:, i] = Fbi/norm(Delta_s)
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

    # rotation matrix from body to wind coordinate system
    alpha = fs.alpha
    beta = fs.beta

    Rb = [cos(beta) -sin(beta) 0;
        sin(beta) cos(beta) 0.0;
        0 0 1]
    Ra = [cos(alpha) 0.0 sin(alpha);
        0 1 0;
        -sin(alpha) 0 cos(alpha)]
    Fw = Rb*Ra*Fb
    Mw = Rb*Ra*Mb

    # stability derivatives
    dRa = SDeriv([-sin(alpha) 0.0 cos(alpha);
            0 0 0;
            -cos(alpha) 0 -sin(alpha)], 0.0, 0.0, 0.0, 0.0)
    dRb = SDeriv(0.0, [-sin(beta) -cos(beta) 0;
            cos(beta) -sin(beta) 0.0;
            0 0 0.0], 0.0, 0.0, 0.0)

    dFw = Rb*Ra*dFb + Rb*dRa*Fb + dRb*Ra*Fb
    dMw = Rb*Ra*dMb + Rb*dRa*Mb + dRb*Ra*Mb

    # rotate distributed forces
    Fwpvec = zeros(3, N)
    for i = 1:N
        Fwpvec[:, i] = Rb*Ra*Fpvec[:, i]
    end

    return Fw, Mw, Vvec, Fwpvec, dFw, dMw
end



"""
a subcomponent of Trefftz plane induced drag calculation.
"""
function Disub(rleft, rright, Gammaj, ri, Gammai, ndsi)
    rho = 1.0  # normalizes away
    
    # left
    rij = ri - rleft
    Vthetai = -cross([-Gammaj; 0.0; 0.0], rij) / (2*pi*norm(rij)^2)
    Di = rho/2.0*Gammai*dot(Vthetai, ndsi)

    # right
    rij = ri - rright
    Vthetai = -cross([Gammaj; 0.0; 0.0], rij) / (2*pi*norm(rij)^2)
    
    Di += rho/2.0*Gammai*dot(Vthetai, ndsi)

    return Di
end


"""
    dic(y, z, rho, symmetric)

Far-field method.  Induced drag coefficient matrix.
"""
function dic(panels::Array{Panel, 1}, Gamma, symmetric)
    
    N = length(panels)
    
    Di = 0.0
    for j = 1:N
        for i = 1:N
            ndsi = normal_vector_magnitude_2d(panels[i])
            ri = mid_point(panels[i])

            # sum up vortex contributions
            Di += Disub(panels[j].rl, panels[j].rr, Gamma[j], ri, Gamma[i], ndsi)
            
            if symmetric && not_on_symmetry_plane(panels[j])
                Di += Disub(flipy(panels[j].rr), flipy(panels[j].rl), Gamma[j], ri, Gamma[i], ndsi)
            end
        end
    end
    
    if symmetric
        Di *= 2
    end

    return Di
end


# function get_lic(rho, Vinf, y, symmetric)
#     # y is QC

#     N = length(y) - 1
#     LIC = zeros(N)
#     factor = symmetric ? 2.0 : 1.0

#     for i = 1:N
#         LIC[i] = factor*rho*Vinf*(y[i+1] - y[i])
#     end

#     return LIC
# end

"""
    function run(vehicle::Vehicle, fs::Freestream, symmetric)

run the vortex lattice method

# Returns 
- CF: force coefficients, wind axes
- CM: moment coefficients, wind axes, about provided c.g.
- ymid, zmid: middle points along quarter chord bound vortex
- l, cl: lift distribution (Lp/(q cref)), lift coefficient distribution (Lp/(q c))
- dCF, dCM: stability derivatives, wind axes
"""
function run(panels::Array{Panel, 1}, ref::Reference, fs::Freestream, symmetric)

    Gamma, dGamma = circulation(panels, ref, fs, symmetric)
    F, M, Vmag, Fp, dF, dM = forces_moments(panels, ref, fs, Gamma, dGamma, symmetric)
    Lp = Fp[3, :]

    # trefftz plane analysis for drag
    Di = dic(panels, Gamma, symmetric)
    F[1] = Di

    # force and moment coefficients
    rho = 1.0  # cancels out from normalization
    Vinf = fs.Vinf
    q = 0.5*rho*Vinf^2
    Sref = ref.S
    bref = ref.b
    cref = ref.c

    CF = F./(q*Sref)
    CM = M./(q*Sref*[bref; cref; bref])
    
    dCF = dF./(q*Sref)
    dCM = dM./(q*Sref*[bref; cref; bref])

    # # normalize p, q, r    
    dCF = SDeriv(dCF.alpha, dCF.beta, dCF.p*2*Vinf/bref, dCF.q*2*Vinf/cref, dCF.r*2*Vinf/bref)
    dCM = SDeriv(dCM.alpha, dCM.beta, dCM.p*2*Vinf/bref, dCM.q*2*Vinf/cref, dCM.r*2*Vinf/bref)

    # rotate p, q, r, to stability axes
    dCF = SDeriv(dCF.alpha, dCF.beta, dCF.p*cos(fs.alpha) + dCF.r*sin(fs.alpha), dCF.q, -dCF.p*sin(fs.alpha) + dCF.r*cos(fs.alpha))
    dCM = SDeriv(dCM.alpha, dCM.beta, dCM.p*cos(fs.alpha) + dCM.r*sin(fs.alpha), dCM.q, -dCM.p*sin(fs.alpha) + dCM.r*cos(fs.alpha))

    # lift and cl dist
    ymid = [mid_point(p)[2] for p in panels]
    zmid = [mid_point(p)[3] for p in panels]
    chord = [p.chord for p in panels]
    cl = Lp./(q*chord)
    l = Lp/(q*cref)
    
    return CF, CM, ymid, zmid, l, cl, dCF, dCM
end

end