struct freestream
    rho
    Vinf
    alpha
    beta
    vother
end

struct panelgeometry
    x
    y
    z
    c
    theta
    phi
    npanels
end


function getVhat(r1, r2)
    
    nr1 = norm(r1)
    nr2 = norm(r2)
    
    f1 = cross(r1, r2)/(nr1*nr2 + dot(r1, r2))
    f2 = (1.0/nr1 + 1.0/nr2)
    f3 = cross(r1, [1.0; 0; 0])/(nr1 - r1[1])
    f4 = 1.0/nr1
    f5 = cross(r2, [1.0; 0; 0])/(nr2 - r2[1])
    f6 = 1.0/nr2

    Vhat = f1*f2 + f3*f4 - f5*f6

    return Vhat
end


function getAIC(geom::panelgeometry, fs::freestream, symmetric)

    N = geom.npanels
    x = geom.x
    y = geom.y
    z = geom.z

    # center points of panels (control points)
    xbar = 0.5*(x[1:N] + x[2:N+1]) + geom.c/2.0
    ybar = 0.5*(y[1:N] + y[2:N+1])
    zbar = 0.5*(z[1:N] + z[2:N+1])

    # rotation matrices
    RA = [cos(fs.alpha) 0.0 sin(fs.alpha);
        0 1 0;
        -sin(fs.alpha) 0 cos(fs.alpha)]
    RB = [cos(fs.beta) -sin(fs.beta) 0.0;
        sin(fs.beta) cos(fs.beta) 0.0;
        0 0 1]
    
    AIC = zeros(N, N)
    for i = 1:N  # CP
        # normal vector body axis
        nb = [sin(geom.theta[i]); 
            -cos(geom.theta[i])*sin(geom.phi[i]);
            cos(geom.theta[i])*cos(geom.phi[i])]
        # normal vector wind axis
        nw = RB*RA*nb

        for j = 1:N  # QC
            r1 = [xbar[i] - x[j], ybar[i] - y[j], zbar[i] - z[j]]
            r2 = [xbar[i] - x[j+1], ybar[i] - y[j+1], zbar[i] - z[j+1]]
            Vij = getVhat(r1, r2)
            AIC[i, j] = dot(Vij, nw)
        end

        if symmetric

            # flip sign for y, but also r1 (left is j+1) now and vice vesa
            for j = 1:N  # QC
                r1 = [xbar[i] - x[j+1], ybar[i] + y[j+1], zbar[i] - z[j+1]]
                r2 = [xbar[i] - x[j], ybar[i] + y[j], zbar[i] - z[j]]
                Vij = getVhat(r1, r2)
                AIC[i, j] += dot(Vij, nw)
            end
        end
    end

    return AIC
end


function Vn_ext(geom::panelgeometry, fs::freestream)
    """fs.vother is a function that should take in xb, yb, zb
    and return [Vxb; Vyb; Vzb].  Can also use nothing."""

    theta = geom.theta
    phi = geom.phi
    
    # freestream velocity in body coordinate
    Vinfb = fs.Vinf*[cos(fs.alpha)*cos(fs.beta);
        -sin(fs.beta);
        sin(fs.alpha)*cos(fs.beta)]
    
    # initialize
    N = length(theta)
    b = zeros(N)

    for i = 1:N

        nb = [sin(theta[i]); 
            -cos(theta[i])*sin(phi[i]);
            cos(theta[i])*cos(phi[i])]

        Vextb = Vinfb
        if fs.vother != nothing
            Vextb += fs.vother(geom.x[i], geom.y[i], geom.z[i])
        end

        b[i] = -dot(Vextb, nb)
    end

    return b
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
    z = zeros(N)
    
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


function getGamma(geom::panelgeometry, fs::freestream, symmetric)

    AIC = getAIC(geom, fs, symmetric)
    b = Vn_ext(geom, fs)

    Gamma = AIC\b

    return Gamma
end




function getDIC(y, z, rho, symmetric)
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


function getLIC(rho, Vinf, y, symmetric)
    # y is QC

    N = length(y) - 1
    LIC = zeros(N)
    factor = symmetric ? 2.0 : 1.0

    for i = 1:N
        LIC[i] = factor*rho*Vinf*(y[i+1] - y[i])
    end

    return LIC
end


function vlm(geom::panelgeometry, fs::freestream, symmetric)

    Gamma = getGamma(geom, fs, symmetric)
    LIC = getLIC(fs.rho, fs.Vinf, geom.y, symmetric)
    DIC = getDIC(geom.y, geom.z, fs.rho, symmetric)

    L = dot(LIC, Gamma)
    Di = dot(Gamma, DIC*Gamma)

    ybar = 0.5*(geom.y[1:end-1] + geom.y[2:end])
    return ybar, Gamma, L, Di
end


symmetric = false

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
vother = nothing
fs = freestream(rho, Vinf, alpha, beta, vother)

y, Gamma, L, Di = vlm(geom, fs, symmetric)

using PyPlot
figure()
plot(y, Gamma)
gcf()