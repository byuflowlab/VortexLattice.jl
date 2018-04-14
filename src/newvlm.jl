struct Freestream
    rho
    Vinf
    alpha
    beta
    Omega  # length: 3 (p, q, r)
    rcg  # length: 3
    vother  # Vother = vother(r).  can also use nothing
end

struct Geometry
    x  # length: npanels + 1
    y
    z
    c  # length: npanels
    theta
    phi
    npanels
    Sref
    cref
    bref
end

struct SDeriv
    alpha
    beta
    p
    q
    r
end

# import Base: +


SDeriv() = SDeriv(zeros(1), zeros(1), zeros(1), zeros(1), zeros(1)) 
SDeriv(N) = SDeriv(zeros(N), zeros(N), zeros(N), zeros(N), zeros(N)) 
SDeriv(M, N) = SDeriv(zeros(M, N), zeros(M, N), zeros(M, N), zeros(M, N), zeros(M, N)) 
Base.:+(x::SDeriv, y::SDeriv) = SDeriv(x.alpha + y.alpha, x.beta + y.beta, x.p + y.p, x.q + y.q, x.r + y.r)
Base.:*(x::Any, y::SDeriv) = SDeriv(x * y.alpha, x * y.beta, x * y.p, x * y.q, x * y.r)
Base.:*(x::SDeriv, y::Any) = SDeriv(x.alpha * y, x.beta * y, x.p * y, x.q * y, x.r * y)
Base.:-(x::SDeriv) = SDeriv(-x.alpha, -x.beta, -x.p, -x.q, -x.r)
Base.:/(x::SDeriv, y::Any) = SDeriv(x.alpha / y, x.beta / y, x.p / y, x.q / y, x.r / y)
Base.dot(x::SDeriv, y::Array{T, 1}) where T<:Any = SDeriv(dot(x.alpha, y), dot(x.beta, y), dot(x.p, y), dot(x.q, y), dot(x.r, y))
Base.cross(x::SDeriv, y::Array{T, 1}) where T<:Any = SDeriv(cross(x.alpha, y), cross(x.beta, y), cross(x.p, y), cross(x.q, y), cross(x.r, y))
Base.sum(x::SDeriv, n) = SDeriv(sum(x.alpha, n), sum(x.beta, n), sum(x.p, n), sum(x.q, n), sum(x.r, n))

function Base.getindex(x::SDeriv, i)
    return SDeriv(x.alpha[i], x.beta[i], x.p[i], x.q[i], x.r[i])
end
function Base.getindex(x::SDeriv, i, j)
    return SDeriv(x.alpha[i, j], x.beta[i, j], x.p[i, j], x.q[i, j], x.r[i, j])
end
function Base.setindex!(x::SDeriv, y::SDeriv, i)
    x.alpha[i] = y.alpha
    x.beta[i] = y.beta
    x.p[i] = y.p
    x.q[i] = y.q
    x.r[i] = y.r
end
function Base.setindex!(x::SDeriv, y::Any, i)
    x.alpha[i] = y
    x.beta[i] = y
    x.p[i] = y
    x.q[i] = y
    x.r[i] = y
end
function Base.setindex!(x::SDeriv, y::SDeriv, i, j)
    x.alpha[i, j] = y.alpha
    x.beta[i, j] = y.beta
    x.p[i, j] = y.p
    x.q[i, j] = y.q
    x.r[i, j] = y.r
end
Base.endof(x) = length(x.alpha)


function vhat_trailing(r1, r2)
    
    nr1 = norm(r1)
    nr2 = norm(r2)
    
    f3 = cross(r1, [1.0; 0; 0])/(nr1 - r1[1])
    f4 = 1.0/nr1
    f5 = cross(r2, [1.0; 0; 0])/(nr2 - r2[1])
    f6 = 1.0/nr2

    Vhat = (f3*f4 - f5*f6)/(4*pi)

    return Vhat
end

function vhat(r1, r2)
    
    nr1 = norm(r1)
    nr2 = norm(r2)

    Vhat = vhat_trailing(r1, r2)
    
    f1 = cross(r1, r2)/(nr1*nr2 + dot(r1, r2))
    f2 = (1.0/nr1 + 1.0/nr2)
    
    Vhat += (f1*f2)/(4*pi)

    return Vhat
end


function mid_points(geom)
    
    N = geom.npanels
    xmid = 0.5*(geom.x[1:N] + geom.x[2:N+1])
    ymid = 0.5*(geom.y[1:N] + geom.y[2:N+1])
    zmid = 0.5*(geom.z[1:N] + geom.z[2:N+1])

    return xmid, ymid, zmid
end

function control_points(geom)

    xcp, ycp, zcp = mid_points(geom)
    xcp += geom.c/2.0

    return xcp, ycp, zcp
end

function normalvector(geom, i)

    nhat = [sin(geom.theta[i]); 
            -cos(geom.theta[i])*sin(geom.phi[i]);
            cos(geom.theta[i])*cos(geom.phi[i])]

    return nhat
end


function aic(geom::Geometry, symmetric)

    N = geom.npanels
    x = geom.x
    y = geom.y
    z = geom.z
    
    # center points of panels (control points)
    xcp, ycp, zcp = control_points(geom)
    
    AIC = zeros(N, N)
    for i = 1:N  # CP
        # normal vector body axis
        nhat = normalvector(geom, i)

        for j = 1:N  # QC
            r1 = [xcp[i] - x[j]; ycp[i] - y[j]; zcp[i] - z[j]]
            r2 = [xcp[i] - x[j+1]; ycp[i] - y[j+1]; zcp[i] - z[j+1]]
            Vij = vhat(r1, r2)
            AIC[i, j] = dot(Vij, nhat)

            if symmetric
                # flip sign for y, but r1 is on left which now corresponds to j+1 and vice vesa
                r1 = [xcp[i] - x[j+1]; ycp[i] + y[j+1]; zcp[i] - z[j+1]]
                r2 = [xcp[i] - x[j]; ycp[i] + y[j]; zcp[i] - z[j]]
                Vij = vhat(r1, r2)
                AIC[i, j] += dot(Vij, nhat)
            end
        end

    end

    return AIC
end


function ext_velocity(fs::Freestream, r)

    # Freestream velocity in body coordinate
    Vinf = fs.Vinf*[cos(fs.alpha)*cos(fs.beta);
        -sin(fs.beta);
        sin(fs.alpha)*cos(fs.beta)]

    Vext = Vinf - cross(fs.Omega, r - fs.rcg)

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
    rvec = r - fs.rcg
    dVdp = [0.0; rvec[3]; -rvec[2]]
    dVdq = [-rvec[3]; 0.0; rvec[1]]
    dVdr = [rvec[2]; -rvec[1]; 0.0]

    dVext = SDeriv(dVda, dVdb, dVdp, dVdq, dVdr)

    return Vext, dVext
end


function vn_ext(geom::Geometry, fs::Freestream)
    
    # initialize
    N = geom.npanels
    theta = geom.theta
    phi = geom.phi
    
    # center points of panels (control points)
    xcp, ycp, zcp = control_points(geom)
    
    b = zeros(N)
    db = SDeriv(N)

    for i = 1:N

        nhat = normalvector(geom, i)

        Vext, dVext = ext_velocity(fs, [xcp[i]; ycp[i]; zcp[i]])

        b[i] = -dot(Vext, nhat)

        # for stability derivatives
        db[i] = -dot(dVext, nhat)
    end

    return b, db
end


function forces_moments(geom::Geometry, fs::Freestream, Gamma, dGamma, symmetric)
    N = geom.npanels
    x = geom.x
    y = geom.y
    z = geom.z

    xmid, ymid, zmid = mid_points(geom)

    Fvec = zeros(3, N)
    Mvec = zeros(3, N)
    Vvec = zeros(3, N)
    Fpvec = zeros(3, N)
    
    dVvec = SDeriv(3, N)
    dFvec = SDeriv(3, N)
    dMvec = SDeriv(3, N)

    for i = 1:N  # control points

        Vindi = 0.0
        dVi = SDeriv(0.0, 0.0, 0.0, 0.0, 0.0)
        for j = 1:N  # vortices  
            r1 = [xmid[i] - x[j]; ymid[i] - y[j]; zmid[i] - z[j]]
            r2 = [xmid[i] - x[j+1]; ymid[i] - y[j+1]; zmid[i] - z[j+1]]
            Vij = vhat_trailing(r1, r2)
            Vindi += Vij*Gamma[j]
        
            dVi += Vij*dGamma[j]

            if symmetric
                # flip sign for y, but r1 is on left which now corresponds to j+1 and vice vesa
                r1 = [xmid[i] - x[j+1]; ymid[i] + y[j+1]; zmid[i] - z[j+1]]
                r2 = [xmid[i] - x[j]; ymid[i] + y[j]; zmid[i] - z[j]]
                Vij = vhat_trailing(r1, r2)
                Vindi += Vij*Gamma[j]

                dVi += Vij*dGamma[j]
            end
        end

        rmidi = [xmid[i]; ymid[i]; zmid[i]]
        Vext, dVext = ext_velocity(fs, rmidi)
        Vvec[:, i] = Vindi + Vext
        
        dVvec[:, i] = dVi + dVext
        
        Delta_s = [x[i+1] - x[i]; y[i+1] - y[i]; z[i+1] - z[i]]
        Fvec[:, i] = rho*Gamma[i]*cross(Vvec[:, i], Delta_s)
        Mvec[:, i] = cross(rmidi - fs.rcg, Fvec[:, i])

        dFvec[:, i] = rho*cross(Gamma[i]*dVvec[:, i] + dGamma[i]*Vvec[:, i], Delta_s)
    
        # force per unit length
        Fpvec[:, i] = Fvec[:, i]/norm(Delta_s)
    end

    Fb = sum(Fvec, 2)
    Mb = sum(Mvec, 2)
    Vmag = sqrt.(Vvec[1, :].^2 + Vvec[2, :].^2 + Vvec[3, :].^2)

    dFb = sum(dFvec, 2)

    if symmetric
        Fb *= 2
        Mb *= 2
        Fb[2] = 0.0
        Mb[1] = 0.0
        Mb[3] = 0.0
        
        dFb *= 2
        dFb[2] = 0.0
    end

    # rotation matrix
    Rb = [cos(fs.beta) -sin(fs.beta) 0;
        sin(fs.beta) cos(fs.beta) 0.0;
        0 0 1]
    Ra = [cos(fs.alpha) 0.0 sin(fs.alpha);
        0 1 0;
        -sin(fs.alpha) 0 cos(fs.alpha)]
    Fw = Rb*Ra*Fb
    Mw = Rb*Ra*Mb

    dRa = SDeriv([-sin(fs.alpha) 0.0 cos(fs.alpha);
            0 0 0;
            -cos(fs.alpha) 0 -sin(fs.alpha)], 0.0, 0.0, 0.0, 0.0)
    dRb = SDeriv(0.0, [-sin(fs.beta) -cos(fs.beta) 0;
            cos(fs.beta) -sin(fs.beta) 0.0;
            0 0 0.0], 0.0, 0.0, 0.0)

    dFw = Rb*Ra*dFb + Rb*dRa*Fb + dRb*Ra*Fb

    Fwpvec = zeros(3, N)
    for i = 1:N
        Fwpvec[:, i] = Rb*Ra*Fpvec[:, i]
    end
    Lp = Fwpvec[3, :]

    return Fw, Mw, Vmag, Lp, dFw
end

function simpletapered(b, AR, λ, Λ, ϕ, θt, npanels, symmetric)
    
    # check number of panels
    if !symmetric && (npanels % 2 != 0)
        warn("must have an even number of panels if nonsymmetric. adding one panel")
        npanels += 1
    end
    
    # geometry parsing
    S = b^2/AR
    cr = 2*S/(b*(1 + λ))
    ct = cr*λ

    # number of vortex lines
    N = npanels + 1

    # number of vortex lines on half of the wing
    Nhalf = symmetric ? N : (N+1)/2
    
    # create geometru
    y = linspace(0, b/2, Nhalf)
    x = y*tan(Λ)
    z = zeros(Nhalf)
    
    ybar = 0.5*(y[1:end-1] + y[2:end])
    c = cr + (ct - cr)*ybar*2/b
    theta = θt*ybar*2/b
    phi = ϕ*ones(npanels)
        
    if !symmetric
        # mirror over y axis
        y = [-y[end:-1:1]; y[2:end]]
        x = [x[end:-1:1]; x[2:end]]
        z = [z[end:-1:1]; z[2:end]]
        
        c = [c[end:-1:1]; c]
        theta = [theta[end:-1:1]; theta]
        phi = [-phi[end:-1:1]; phi]
    end

    geom = Geometry(x, y, z, c, theta, phi, npanels, S, S/b, b)

    return geom
end


function circulation(geom::Geometry, fs::Freestream, symmetric)

    AIC = aic(geom, symmetric)
    b, db = vn_ext(geom, fs)

    Gamma = AIC\b

    dGamma = SDeriv(AIC\db.alpha, AIC\db.beta, AIC\db.p, AIC\db.q, AIC\db.r)

    return Gamma, dGamma
end






function dic(y, z, rho, symmetric)
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

    # DIC matrix
    DIC = zeros(N, N)

    if symmetric
        for j = 1:N
            for i = 1:N
                DIC[i, j] = rho/(2*pi)*(getk(i, j) - getk(i, j+1) - getk(i, -j) + getk(i, -(j+1)))
            end
        end
    else  # not symmetric
        for j = 1:N
            for i = 1:N
                DIC[i, j] = rho/(4*pi)*(getk(i, j) - getk(i, j+1))
            end
        end
    end

    return DIC
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


function vlm(geom::Geometry, fs::Freestream, symmetric)

    Gamma, dGamma = circulation(geom, fs, symmetric)
    F, M, Vmag, Lp, dF = forces_moments(geom, fs, Gamma, dGamma, symmetric)

    # trefftz plane analysis for drag
    DIC = dic(geom.y, geom.z, fs.rho, symmetric)
    Di = dot(Gamma, DIC*Gamma)
    F[1] = Di

    # force and moment coefficients
    CF = F/(0.5*rho*Vinf^2*geom.Sref)
    CM = M./(0.5*rho*Vinf^2*[geom.bref; geom.cref; geom.bref])

    dCF = dF/(0.5*rho*Vinf^2*geom.Sref)

    # lift and cl dist
    xmid, ymid, zmid = mid_points(geom)
    cl = Lp./(0.5*rho*Vinf^2*geom.c)
    l = Lp/(0.5*rho*Vinf^2*geom.cref)
    
    return CF, CM, ymid, zmid, l, cl, dCF
end


symmetric = true

b = 5.0
AR = 8.0
λ = 0.6
Λ = 20.0*pi/180
ϕ = 0.0
θt = -2.0*pi/180
npanels = 100
geom = simpletapered(b, AR, λ, Λ, ϕ, θt, npanels, symmetric)

rho = 1.0
Vinf = 2.0
alpha = 5.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
rcg = [0.0; 0.0; 0.0]
vother = nothing
fs = Freestream(rho, Vinf, alpha, beta, Omega, rcg, vother)

CF, CM, ymid, zmid, l, cl, dCF = vlm(geom, fs, symmetric)

using PyPlot
figure()
plot(ymid, l)
plot(ymid, cl)
gcf()