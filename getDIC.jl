function getDIC(yFF,zFF,rho,yFF2,zFF2)
    if 3 == 3
        yFF2 = yFF
        zFF2 = zFF
    end

    n = length(yFF)-1
    n2 = length(yFF2)-1

     # ------------ define useful variables ------------------
    yc = 1/2*(yFF[2:n+1] + yFF[1:n])   # center of panels
    zc = 1/2*(zFF[2:n+1] + zFF[1:n])
    ny = zFF[2:n+1] - zFF[1:n]
    nz = yFF[2:n+1] - yFF[1:n]
     # ----------------------------------------------------

    DIC = zeros(n,n2)
     #------- normal wash calculation ----------------
    for i = 1:n
        for j = 1:n2
            ry = yFF2[j] - yc[i]
            rz = zFF2[j] - zc[i]
            r2 = ry^2+rz^2
            DIC[i,j] = rho/4/pi/r2*(-1)*(rz*ny[i]+ry*nz[i])

             # add other side of panel
            ry = yFF2[j+1] - yc[i]
            rz = zFF2[j+1] - zc[i]
            r2 = ry^2+rz^2
            DIC[i,j] = DIC[i,j] + rho/4/pi/r2*(1)*(rz*ny[i]+ry*nz[i])
        end
    end
    return DIC
end

