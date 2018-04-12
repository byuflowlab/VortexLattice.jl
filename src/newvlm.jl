struct freestream
    rho
    Vinf
    alpha
    beta
    Omega  # length: 3
    rcg  # length: 3
    vother  # Vother = vother(r).  can also use nothing
end

struct panelgeometry
    x  # length: npanels + 1
    y
    z
    c  # length: npanels
    theta
    phi
    npanels
end


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


function aic(geom::panelgeometry, symmetric)

    N = geom.npanels
    x = geom.x
    y = geom.y
    z = geom.z
    
    # center points of panels (control points)
    xcp, ycp, zcp = control_points(geom)
    
    AIC = zeros(N, N)
    for i = 1:N  # CP
        # normal vector body axis
        nhat = [sin(geom.theta[i]); 
            -cos(geom.theta[i])*sin(geom.phi[i]);
            cos(geom.theta[i])*cos(geom.phi[i])]

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


function ext_velocity(fs::freestream, r)

    # freestream velocity in body coordinate
    Vinf = fs.Vinf*[cos(fs.alpha)*cos(fs.beta);
        -sin(fs.beta);
        sin(fs.alpha)*cos(fs.beta)]

    Vext = Vinf - cross(fs.Omega, r - fs.rcg)

    if fs.vother != nothing
        Vext += fs.vother(r)
    end

    return Vext
end


function vn_ext(geom::panelgeometry, fs::freestream)
    
    # initialize
    N = geom.npanels
    theta = geom.theta
    phi = geom.phi
    
    # center points of panels (control points)
    xcp, ycp, zcp = control_points(geom)
    
    b = zeros(N)
    for i = 1:N

        nhat = [sin(theta[i]); 
            -cos(theta[i])*sin(phi[i]);
            cos(theta[i])*cos(phi[i])]

        Vext = ext_velocity(fs, [xcp[i]; ycp[i]; zcp[i]])

        b[i] = -dot(Vext, nhat)
    end

    return b
end


function forces_moments(geom::panelgeometry, fs::freestream, Gamma, symmetric)
    N = geom.npanels
    x = geom.x
    y = geom.y
    z = geom.z

    xmid, ymid, zmid = mid_points(geom)

    Fvec = zeros(3, N)
    Mvec = zeros(3, N)

    for i = 1:N  # control points

        Vindi = 0.0
        for j = 1:N  # vortices  
            r1 = [xmid[i] - x[j]; ymid[i] - y[j]; zmid[i] - z[j]]
            r2 = [xmid[i] - x[j+1]; ymid[i] - y[j+1]; zmid[i] - z[j+1]]
            Vij = vhat_trailing(r1, r2)
            Vindi += Vij*Gamma[j]
        end

        rmidi = [xmid[i]; ymid[i]; zmid[i]]
        Vext = ext_velocity(fs, rmidi)
        Vi = Vindi + Vext
        
        Deltas = [x[i+1] - x[i]; y[i+1] - y[i]; z[i+1] - z[i]]
        Fvec[:, i] = rho*Gamma[i]*cross(Vi, Deltas)
        Mvec[:, i] = cross(rmidi - fs.rcg, Fvec[:, i])
    end

    Fb = sum(Fvec, 2)
    Mb = sum(Mvec, 2)

    # rotation matrix
    Rb = [cos(fs.beta) -sin(fs.beta) 0;
        sin(fs.beta) cos(fs.beta) 0.0;
        0 0 1]
    Ra = [cos(fs.alpha) 0.0 sin(fs.alpha);
        0 1 0;
        -sin(fs.alpha) 0 cos(fs.alpha)]
    Fw = Rb*Ra*Fb
    Mw = Rb*Ra*Mb
    
    return Fw, Mw
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
        
    if symmetric
        geom = panelgeometry(x, y, z, c, theta, phi, npanels)
    else
        # mirror over y axis
        y = [-y[end:-1:1]; y[2:end]]
        x = [x[end:-1:1]; x[2:end]]
        z = [z[end:-1:1]; z[2:end]]
        
        c = [c[end:-1:1]; c]
        theta = [theta[end:-1:1]; theta]
        phi = [-phi[end:-1:1]; phi]

        geom = panelgeometry(x, y, z, c, theta, phi, npanels)
    end

    return geom
end


function circulation(geom::panelgeometry, fs::freestream, symmetric)

    AIC = aic(geom, symmetric)
    b = vn_ext(geom, fs)

    Gamma = AIC\b

    return Gamma
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


function vlm(geom::panelgeometry, fs::freestream, symmetric)

    Gamma = circulation(geom, fs, symmetric)
    F, M = forces_moments(geom, fs, Gamma, symmetric)

    # trefftz plane analysis for drag
    DIC = dic(geom.y, geom.z, fs.rho, symmetric)
    Di = dot(Gamma, DIC*Gamma)
    F[1] = Di

    # lift and cl dist
    xmid, ymid, zmid = mid_points(geom)
    cl = 2*Gamma./(fs.Vinf*geom.c)
    l = 2*Gamma./(fs.Vinf*mean(geom.c))
    
    return F, M, ymid, zmid, l, cl
end


symmetric = true

b = 5.0
AR = 8.0
λ = 0.6
Λ = 0.0
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
fs = freestream(rho, Vinf, alpha, beta, Omega, rcg, vother)

F, M, ymid, zmid, l, cl = vlm(geom, fs, symmetric)

using PyPlot
figure()
plot(ymid, l)
plot(ymid, cl)
gcf()