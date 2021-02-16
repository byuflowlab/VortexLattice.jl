export linearsections



"""linearly interpolate between two points, eta is fraction between 0 (rstart) and 1 (rend)"""
function linearinterp(eta, rstart, rend)

    return (1-eta)*rstart + eta*rend
end


function spanwise_spacing(n, stype)

    if stype == "u"  # uniform

        eta = range(0, 1.0, length=n)
        eta_mid = linearinterp(0.5, eta[1:end-1], eta[2:end])

    elseif stype == "c"  # cosine

        theta = range(0, pi, length=n)
        eta = (1.0 .- cos.(theta))/2.0

        # note that control points are also placed with cosine spacing as this improves accuracy tremendously
        theta_mid = linearinterp(0.5, theta[1:end-1], theta[2:end])
        eta_mid = (1.0 .- cos.(theta_mid))/2.0

    else
        error("invalid spacing type: ")
    end

    return eta, eta_mid

end


function chordwise_spacing(n, stype)

    if stype == "u"  # uniform

        eta = range(0, 1.0, length=n)

    elseif stype == "c"  # cosine

        theta = range(0, pi, length=n)
        eta = (1.0 .- cos.(theta))/2.0

    else
        error("invalid spacing type: ")
    end

    eta_qtr = linearinterp(0.25, eta[1:end-1], eta[2:end])
    eta_thrqtr = linearinterp(0.75, eta[1:end-1], eta[2:end])

    return eta_qtr, eta_thrqtr

end


#=
Panel order numbering:

1 ------- 2
|         |
|         |
3-------- 4

=#


"""discretize a linear segement of a wing"""
function creategrid_onesegment(r1, r2, r3, r4, ns, spacing_s, nc, spacing_c, thetaL, thetaR)


    etas, etabar = spanwise_spacing(ns+1, spacing_s)
    eta_qtr, eta_thrqtr = chordwise_spacing(nc+1, spacing_c)
    chord1 = r3[1] - r1[1]
    chord2 = r4[1] - r2[1]

    panels = Vector{Panel}(undef, ns*nc)
    for j = 1:ns
        for i = 1:nc

            rtop = linearinterp(etas[j], r1, r2)
            rbot = linearinterp(etas[j], r3, r4)
            rl = linearinterp(eta_qtr[i], rtop, rbot)

            rtop = linearinterp(etas[j+1], r1, r2)
            rbot = linearinterp(etas[j+1], r3, r4)
            rr = linearinterp(eta_qtr[i], rtop, rbot)

            rtop = linearinterp(etabar[j], r1, r2)
            rbot = linearinterp(etabar[j], r3, r4)
            rcp = linearinterp(eta_thrqtr[i], rtop, rbot)

            theta = linearinterp(etabar[j], thetaL, thetaR)
            chord = linearinterp(etabar[j], chord1, chord2)

            panels[(j-1)*nc + i] = Panel(rl, rr, rcp, theta, chord)

        end
    end


    return panels
end


function linearsections(xle, yle, zle, chord, theta, ns, spacing_s, nc, spacing_c, duplicate)

    panels = Vector{Panel}(undef, 0)

    n = length(ns)
    for i = 1:n
        r1 = [xle[i]; yle[i]; zle[i]]
        r2 = [xle[i+1]; yle[i+1]; zle[i+1]]
        r3 = r1 + [chord[i]; 0; 0]
        r4 = r2 + [chord[i+1]; 0; 0]

        panels = [panels; creategrid_onesegment(r1, r2, r3, r4, ns[i], spacing_s[i], nc[i], spacing_c[i], theta[i], theta[i+1])]
    end

    if duplicate

        for i = 1:n
            r1 = [xle[i+1]; -yle[i+1]; zle[i+1]]
            r2 = [xle[i]; -yle[i]; zle[i]]
            r3 = r1 + [chord[i+1]; 0; 0]
            r4 = r2 + [chord[i]; 0; 0]

            panels = [panels; creategrid_onesegment(r1, r2, r3, r4, ns[i], spacing_s[i], nc[i], spacing_c[i], theta[i], theta[i+1])]
        end

    end

    return panels
end


function translate!(panels, r)

    npanels = length(panels)

    for i = 1:npanels
        panels[i].rl[:] += r[:]
        panels[i].rr[:] += r[:]
        panels[i].rcp[:] += r[:]
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

    return linearsections(xle, yle, zle, chord, theta, [npanels], [spacing], [1], ["u"], duplicate)

end


# # TODO: With Pkg3 I can make this a separate module within same repo so base VLM doesn't depending on plotting.
# function visualizegeometry(panels)

#     for i = 1:length(panels)

#         x = [panels[i].rl[1]; panels[i].rr[1]]
#         y = [panels[i].rl[2]; panels[i].rr[2]]
#         z = [panels[i].rl[3]; panels[i].rr[3]]
#         plot3D(x, y, z, color="b")
#         plot3D(x, -y, z, color="b")
#         x = [panels[i].rl[1]; panels[i].rl[1] + 3.0/4*panels[i].chord]
#         y = [panels[i].rl[2]; panels[i].rl[2]]
#         z = [panels[i].rl[3]; panels[i].rl[3]]
#         plot3D(x, y, z, color="b")
#         plot3D(x, -y, z, color="b")
#         x = [panels[i].rr[1]; panels[i].rr[1] + 3.0/4*panels[i].chord]
#         y = [panels[i].rr[2]; panels[i].rr[2]]
#         z = [panels[i].rr[3]; panels[i].rr[3]]
#         plot3D(x, y, z, color="b")
#         plot3D(x, -y, z, color="b")
#         x = [panels[i].rl[1] - 1.0/4*panels[i].chord; panels[i].rr[1] - 1.0/4*panels[i].chord]
#         y = [panels[i].rl[2]; panels[i].rr[2]]
#         z = [panels[i].rl[3]; panels[i].rr[3]]
#         plot3D(x, y, z, color="0.5")
#         plot3D(x, -y, z, color="0.5")
#         x = [panels[i].rl[1] + 3.0/4*panels[i].chord; panels[i].rr[1] + 3.0/4*panels[i].chord]
#         y = [panels[i].rl[2]; panels[i].rr[2]]
#         z = [panels[i].rl[3]; panels[i].rr[3]]
#         plot3D(x, y, z, color="0.5")
#         plot3D(x, -y, z, color="0.5")


#     end

#     grid("off")
#     # gca()[:view_init](90.0, 0.0)
#     gca()[:view_init](20, -135)
#     axis("equal")
#     xlabel("x")
#     ylabel("y")
#     zlabel("z")
# end



# function visualizeoutline2d(panels)

#     for i = 1:length(panels)

#         x = [panels[i].rl[1] - 1.0/4*panels[i].chord; panels[i].rr[1] - 1.0/4*panels[i].chord]
#         y = [panels[i].rl[2]; panels[i].rr[2]]
#         plot(y, -x, color="k")
#         plot(-y, -x, color="k")
#         x = [panels[i].rl[1] + 3.0/4*panels[i].chord; panels[i].rr[1] + 3.0/4*panels[i].chord]
#         y = [panels[i].rl[2]; panels[i].rr[2]]
#         plot(y, -x, color="k")
#         plot(-y, -x, color="k")
#     end

#     x = [panels[end].rr[1] - 1.0/4*panels[end].chord; panels[end].rr[1] + 3.0/4*panels[end].chord]
#     y = [panels[end].rr[2]; panels[end].rr[2]]
#     plot(y, -x, color="k")
#     p, = plot(-y, -x, color="k")

#     axis("equal")
#     axis("off")

#     return p
# end



# function twoview(panels)

#     p = visualizeoutline2d(panels)
#     xleft = panels[1].rl[1] + 3.0/4*panels[1].chord
#     xright = panels[end].rr[1] + 3.0/4*panels[end].chord
#     xmax = max(xleft, xright)
#     xoffset = xmax + panels[1].chord

#     p2 = p
#     for i = 1:length(panels)

#         y = [panels[i].rl[2]; panels[i].rr[2]]
#         z = [panels[i].rl[3]; panels[i].rr[3]]

#         plot(y, z - xoffset, color="r")
#         p2, = plot(-y, z - xoffset, color="r")
#     end

#     legend([p, p2], ["top view", "back view"])

# end
