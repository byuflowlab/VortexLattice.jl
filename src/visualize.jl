using PyPlot

function visualizegeometry(panels)
    for i = 1:length(panels)
        
        x = [panels[i].rl[1]; panels[i].rr[1]]
        y = [panels[i].rl[2]; panels[i].rr[2]]
        z = [panels[i].rl[3]; panels[i].rr[3]]
        plot3D(x, y, z, color="b")
        plot3D(x, -y, z, color="b")
        x = [panels[i].rl[1]; panels[i].rl[1] + 3.0/4*panels[i].chord]
        y = [panels[i].rl[2]; panels[i].rl[2]]
        z = [panels[i].rl[3]; panels[i].rl[3]]
        plot3D(x, y, z, color="b")
        plot3D(x, -y, z, color="b")
        x = [panels[i].rr[1]; panels[i].rr[1] + 3.0/4*panels[i].chord]
        y = [panels[i].rr[2]; panels[i].rr[2]]
        z = [panels[i].rr[3]; panels[i].rr[3]]
        plot3D(x, y, z, color="b")
        plot3D(x, -y, z, color="b")
        x = [panels[i].rl[1] - 1.0/4*panels[i].chord; panels[i].rr[1] - 1.0/4*panels[i].chord]
        y = [panels[i].rl[2]; panels[i].rr[2]]
        z = [panels[i].rl[3]; panels[i].rr[3]]
        plot3D(x, y, z, color="0.5")
        plot3D(x, -y, z, color="0.5")
        x = [panels[i].rl[1] + 3.0/4*panels[i].chord; panels[i].rr[1] + 3.0/4*panels[i].chord]
        y = [panels[i].rl[2]; panels[i].rr[2]]
        z = [panels[i].rl[3]; panels[i].rr[3]]
        plot3D(x, y, z, color="0.5")
        plot3D(x, -y, z, color="0.5")

        
    end
    
    grid("off")
    # gca()[:view_init](90.0, 0.0)
    gca()[:view_init](20, -135)
    axis("equal")
    xlabel("x")
    ylabel("y")
    zlabel("z")
end
