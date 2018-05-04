
"""
A simple tapered wing
"""
function simpletapered(b, AR, λ, Λ, ϕ, θr, θt, Sref, cref, bref, npanels, symmetric, spacing)
    
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
    
    # create basic geometry
    if spacing == "uniform"
        y = linspace(0, b/2, Nhalf)
    elseif spacing == "cosine"
        y = b/2*cos.(linspace(pi/2, 0, Nhalf))
    end
    x = y*tan(Λ)

    # linear chord/twist
    ybar = 0.5*(y[1:end-1] + y[2:end])
    c = cr + (ct - cr)*ybar*2/b
    theta = θr + (θt - θr)*ybar*2/b

    # rotate for dihedral
    z = y*sin(ϕ)
    y = y*cos(ϕ)
        
    if !symmetric
        # mirror over y axis
        y = [-y[end:-1:1]; y[2:end]]
        x = [x[end:-1:1]; x[2:end]]
        z = [z[end:-1:1]; z[2:end]]
        
        c = [c[end:-1:1]; c]
        theta = [theta[end:-1:1]; theta]
    end

    geom = Geometry(x, y, z, c, theta, npanels, Sref, cref, bref)

    return geom
end


