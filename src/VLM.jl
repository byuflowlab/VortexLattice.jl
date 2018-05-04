module VLM

struct Freestream
    rho
    Vinf
    alpha
    beta
    Omega  # length: 3 (p, q, r)
    rcg  # length: 3
    vother  # Vother = vother(x, y, z).  can also use nothing
end

struct Geometry
    x  # length: npanels + 1
    y
    z
    c  # length: npanels
    theta
    npanels
    Sref
    cref
    bref
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
    vhat_horseshoe(xcp, ycp, zcp, i, x, y, z, j, symmetric)

Computes the induced velocity (per unit circulation) for a horseshoe vortex, with trailing
vortices in the +x direction.  The velocity is computed at control point i
with vortices j and j+1
"""
function vhat(xcp, ycp, zcp, i, x, y, z, j, symmetric, include_bound)
    
    if include_bound
        vfunc = vhat_horseshoe
    else
        vfunc = vhat_trailing
    end

    r1 = [xcp[i] - x[j]; 
            ycp[i] - y[j]; 
            zcp[i] - z[j]]
    r2 = [xcp[i] - x[j+1]; 
            ycp[i] - y[j+1]; 
            zcp[i] - z[j+1]]
        
    Vhat = vfunc(r1, r2)

    if symmetric  # add contribution from other side
        # flip sign for y, but r1 is on left which now corresponds to j+1 and vice vesa
        r1 = [xcp[i] - x[j+1]; 
              ycp[i] + y[j+1]; 
              zcp[i] - z[j+1]]
        r2 = [xcp[i] - x[j]; 
              ycp[i] + y[j]; 
              zcp[i] - z[j]]

        Vhat += vfunc(r1, r2)
    end

    return Vhat
end


"""
    mid_points(geom::Geometry)

Compute the positions of the middle of each horseshoe vortex.
"""
function mid_points(geom::Geometry)
    
    N = geom.npanels
    xmid = 0.5*(geom.x[1:N] + geom.x[2:N+1])
    ymid = 0.5*(geom.y[1:N] + geom.y[2:N+1])
    zmid = 0.5*(geom.z[1:N] + geom.z[2:N+1])

    return xmid, ymid, zmid
end


"""
    control_points(geom::Geometry)

Compute the positions of the control points for each horseshoe vortex.
"""
function control_points(geom::Geometry)

    xcp, ycp, zcp = mid_points(geom)
    xcp += geom.c/2.0  # add 3/4 chord point

    return xcp, ycp, zcp
end


"""
    normalvector(geom::Geometry, i)

Compute the normal vector for section i of the wing.
"""
function normalvector(geom::Geometry, i)

    dy = geom.y[i+1] - geom.y[i]
    dz = geom.z[i+1] - geom.z[i]
    ds = sqrt(dy^2 + dz^2)
    sphi = dz/ds
    cphi = dy/ds

    nhat = [sin(geom.theta[i]); 
            -cos(geom.theta[i])*sphi;
            cos(geom.theta[i])*cphi]

    return nhat
end



"""
    aic(geom::Geometry, symmetric)

aerodynamic influence coefficients
"""
function aic(geom::Geometry, symmetric)

    N = geom.npanels
    x = geom.x
    y = geom.y
    z = geom.z
    
    # control point locations
    xcp, ycp, zcp = control_points(geom)

    include_bound = true  # include bound vortices
    
    AIC = zeros(N, N)
    for i = 1:N  # CP

        # normal vector body axis
        nhat = normalvector(geom, i)

        for j = 1:N  # QC
            Vij = vhat(xcp, ycp, zcp, i, x, y, z, j, symmetric, include_bound)
            AIC[i, j] = dot(Vij, nhat)
        end

    end

    return AIC
end

"""
    ext_velocity(fs::Freestream, x, y, z)

Compute the external velocity at location x, y, z,
and the corresonding (partial) stability derivatives
"""
function ext_velocity(fs::Freestream, x, y, z)

    # Freestream velocity in body coordinate
    Vinf = fs.Vinf*[cos(fs.alpha)*cos(fs.beta);
        -sin(fs.beta);
        sin(fs.alpha)*cos(fs.beta)]

    Vext = Vinf - cross(fs.Omega, [x; y; z] - fs.rcg)

    if fs.vother != nothing
        Vext += fs.vother(x, y, z)
    end

    # for stability derivatives
    dVda = fs.Vinf*[-sin(fs.alpha)*cos(fs.beta);
        0.0;
        cos(fs.alpha)*cos(fs.beta)]
    dVdb = fs.Vinf*[-cos(fs.alpha)*sin(fs.beta);
        -cos(fs.beta);
        -sin(fs.alpha)*sin(fs.beta)]
    rvec = [x; y; z] - fs.rcg
    dVdp = [0.0; rvec[3]; -rvec[2]]
    dVdq = [-rvec[3]; 0.0; rvec[1]]
    dVdr = [rvec[2]; -rvec[1]; 0.0]

    dVext = SDeriv(dVda, dVdb, dVdp, dVdq, dVdr)

    return Vext, dVext
end


"""
    vn_ext(geom::Geometry, fs::Freestream)

Compute the normal component of the external velocity along the geometry.
This forms the right hand side of the circulation linear system solve.
"""
function vn_ext(geom::Geometry, fs::Freestream)
    
    # control point locations
    xcp, ycp, zcp = control_points(geom)
    
    # initialize
    N = geom.npanels
    b = zeros(N)
    db = SDeriv(N)

    # iterate through panels
    for i = 1:N

        # normal vector
        nhat = normalvector(geom, i)

        # external velocity
        Vext, dVext = ext_velocity(fs, xcp[i], ycp[i], zcp[i])

        # right hand side vector
        b[i] = -dot(Vext, nhat)

        # (partial) stability derivatives
        db[i] = -dot(dVext, nhat)
    end

    return b, db
end


"""
    forces_moments(geom::Geometry, fs::Freestream, Gamma, dGamma, symmetric)

Computes the forces and moments acting on the aircraft using the given circulation.
"""
function forces_moments(geom::Geometry, fs::Freestream, Gamma, dGamma, symmetric)
    N = geom.npanels
    x = geom.x
    y = geom.y
    z = geom.z

    # forces evaluated at midpoint of panels (on bound vortex)
    xmid, ymid, zmid = mid_points(geom)

    # initialize
    Fb = zeros(3)  # forces
    Mb = zeros(3)  # moments
    Fpvec = zeros(3, N)  # distributed forces
    Vvec = zeros(3, N)  # distributed velocity
    
    dFb = SDeriv(3)
    dMb = SDeriv(3)

    include_bound = false  # bound vortices don't contribute

    # compute induced velocity at quarter-quard midpoints (rather than at control points)
    for i = 1:N  # control points

        Vindi = 0.0
        dVindi = SDeriv(0.0, 0.0, 0.0, 0.0, 0.0)

        for j = 1:N  # vortices  
            Vij = vhat(xmid, ymid, zmid, i, x, y, z, j, symmetric, include_bound)
            Vindi += Vij*Gamma[j]
        
            dVindi += Vij*dGamma[j]
        end

        # add external velocity
        Vext, dVext = ext_velocity(fs, xmid[i], ymid[i], zmid[i])
        Vi = Vindi + Vext
        dVi = dVindi + dVext
        Vvec[:, i] = Vi  # save velocity distribution in a vector
        
        # forces and moments
        Delta_s = [x[i+1] - x[i]; 
                   y[i+1] - y[i]; 
                   z[i+1] - z[i]]
        Fbi = fs.rho*Gamma[i]*cross(Vi, Delta_s)
        Fb += Fbi
        Mb += cross([xmid[i]; ymid[i]; zmid[i]] - fs.rcg, Fbi)

        dFbi = fs.rho*cross(Gamma[i]*dVi + dGamma[i]*Vi, Delta_s)
        dFb += dFbi
        dMb += cross([xmid[i]; ymid[i]; zmid[i]] - fs.rcg, dFbi)
    
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
    Rb = [cos(fs.beta) -sin(fs.beta) 0;
        sin(fs.beta) cos(fs.beta) 0.0;
        0 0 1]
    Ra = [cos(fs.alpha) 0.0 sin(fs.alpha);
        0 1 0;
        -sin(fs.alpha) 0 cos(fs.alpha)]
    Fw = Rb*Ra*Fb
    Mw = Rb*Ra*Mb

    # stability derivatives
    dRa = SDeriv([-sin(fs.alpha) 0.0 cos(fs.alpha);
            0 0 0;
            -cos(fs.alpha) 0 -sin(fs.alpha)], 0.0, 0.0, 0.0, 0.0)
    dRb = SDeriv(0.0, [-sin(fs.beta) -cos(fs.beta) 0;
            cos(fs.beta) -sin(fs.beta) 0.0;
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
    circulation(geom::Geometry, fs::Freestream, symmetric)

Solve for circulation distribution.
"""
function circulation(geom::Geometry, fs::Freestream, symmetric)

    AIC = aic(geom, symmetric)
    b, db = vn_ext(geom, fs)

    Gamma = AIC\b

    dGamma = SDeriv(AIC\db.alpha, AIC\db.beta, AIC\db.p, AIC\db.q, AIC\db.r)

    return Gamma, dGamma
end


"""
    dic(y, z, rho, symmetric)

Far-field method.  Induced drag coefficient matrix.
"""
function dic(Gamma, y, z, rho, symmetric)
    # y, z are QC
    
    # #panels is one less than #vortices
    N = length(y) - 1

    # center points of panels
    ybar = 0.5*(y[1:N] + y[2:N+1])
    zbar = 0.5*(z[1:N] + z[2:N+1])

    function getk(i, j)
        signj = sign(j)
        j = abs(j)

        ry = ybar[i] - signj*y[j]
        rz = zbar[i] - z[j]
        r = ry^2 + rz^2
        kij = (rz*(z[i+1] - z[i]) + ry*(y[i+1] - y[i]))/r

        return kij
    end

    Di = 0.0

    if symmetric
        for j = 1:N
            for i = 1:N
                Di += rho*Gamma[i]*Gamma[j]/(2*pi)*(getk(i, j) - getk(i, j+1) - getk(i, -j) + getk(i, -(j+1)))
            end
        end
    else  # not symmetric
        for j = 1:N
            for i = 1:N
                Di += rho*Gamma[i]*Gamma[j]/(4*pi)*(getk(i, j) - getk(i, j+1))
            end
        end
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
    function run(geom::Geometry, fs::Freestream, symmetric)

run the vortex lattice method

# Returns 
- CF: force coefficients, wind axes
- CM: moment coefficients, wind axes, about provided c.g.
- ymid, zmid: middle points along quarter chord bound vortex
- l, cl: lift distribution (Lp/(q cref)), lift coefficient distribution (Lp/(q c))
- dCF, dCM: stability derivatives, wind axes
"""
function run(geom::Geometry, fs::Freestream, symmetric)

    Gamma, dGamma = circulation(geom, fs, symmetric)
    F, M, Vmag, Fp, dF, dM = forces_moments(geom, fs, Gamma, dGamma, symmetric)
    Lp = Fp[3, :]

    # trefftz plane analysis for drag
    Di = dic(Gamma, geom.y, geom.z, fs.rho, symmetric)
    F[1] = Di

    # force and moment coefficients
    q = 0.5*fs.rho*fs.Vinf^2

    CF = F./(q*geom.Sref)
    CM = M./(q*geom.Sref*[geom.bref; geom.cref; geom.bref])
    
    dCF = dF./(q*geom.Sref)
    dCM = M./(q*geom.Sref*[geom.bref; geom.cref; geom.bref])

    # lift and cl dist
    xmid, ymid, zmid = mid_points(geom)
    cl = Lp./(q*geom.c)
    l = Lp/(q*geom.cref)
    
    return CF, CM, ymid, zmid, l, cl, dCF, dCM
end

end

# symmetric = true

# b = 5.0
# AR = 8.0
# λ = 0.6
# Λ = 20.0*pi/180
# ϕ = 0.0
# θt = -2.0*pi/180
# npanels = 100
# geom = simpletapered(b, AR, λ, Λ, ϕ, θt, npanels, symmetric)

# rho = 1.0
# Vinf = 2.0
# alpha = 5.0*pi/180
# beta = 0.0
# Omega = [0.0; 0.0; 0.0]
# rcg = [0.0; 0.0; 0.0]
# vother = nothing
# fs = Freestream(rho, Vinf, alpha, beta, Omega, rcg, vother)

# CF, CM, ymid, zmid, l, cl, dCF = vlm(geom, fs, symmetric)

# using PyPlot
# figure()
# plot(ymid, l)
# plot(ymid, cl)
# gcf()