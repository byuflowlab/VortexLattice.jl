

function linearsections(xle, yle, zle, chord, theta, npanels, duplicate, spacing)

    # initialize
    nsections = length(xle) - 1  # number of segments on lifting surface
    nptotal = sum(npanels)  # total number of panels
    if duplicate
        panels = Array{Panel}(2*nptotal)  # twice the number of panels
    else
        panels = Array{Panel}(nptotal)
    end
    start = 0

    for i = 1:nsections

        # create basic geometry
        if spacing == "uniform"
            eta = linspace(0, 1, npanels[i]+1)
        elseif spacing == "cosine"
            eta = cos.(linspace(pi/2, 0, npanels[i]+1))
        end

        # midpoints
        etabar = 0.5*(eta[1:end-1] + eta[2:end])
        
        # linearly interpolate
        xvec = (1 - eta)*(xle[i] + chord[i]/4.0) + eta*(xle[i+1] + chord[i+1]/4.0)
        yvec = (1 - eta)*yle[i] + eta*yle[i+1]
        zvec = (1 - eta)*zle[i] + eta*zle[i+1]
        cvec = (1 - etabar)*chord[i] + etabar*chord[i+1]
        tvec = (1 - etabar)*theta[i] + etabar*theta[i+1]
        

        # create panels
        for j = 1:npanels[i]
            rleft = [xvec[j]; yvec[j]; zvec[j]]
            rright = [xvec[j+1]; yvec[j+1]; zvec[j+1]]
            panels[start + j] = Panel(rleft, rright, cvec[j], tvec[j])
            
            if duplicate
                rleft = [xvec[j+1]; -yvec[j+1]; zvec[j+1]]
                rright = [xvec[j]; -yvec[j]; zvec[j]]
                panels[nptotal + start + j] = Panel(rleft, rright, cvec[j], tvec[j])
            end
        end

        # update starting index
        start += npanels[i]
    end

    return panels
end

function translate!(panels, r)

    npanels = length(panels)

    for i = 1:npanels
        panels[i].rl[:] += r[:]
        panels[i].rr[:] += r[:]
    end
end


function simplewing(b, AR, λ, Λ, ϕ, θr, θt, npanels, duplicate, spacing)
    
    # geometry parsing
    S = b^2/AR
    cr = 2*S/(b*(1 + λ))
    ct = cr*λ

    xle = [0.0; cr/4.0 + b/2.0*tan(Λ) - ct/4.0]
    yle = [0.0; b/2*cos(ϕ)]
    zle = [0.0; b/2*sin(ϕ)]
    chord = [cr; ct]
    theta = [θr; θt]

    return linearsections(xle, yle, zle, chord, theta, [npanels], duplicate, spacing)

end


