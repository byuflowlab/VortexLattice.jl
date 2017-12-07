module NVLM

struct CPdata
    chord::Array{Float64, 1}
    twist::Array{Float64, 1}
    tc::Array{Float64, 1}
    x::Array{Float64, 1}
    y::Array{Float64, 1}
    z::Array{Float64, 1}
    dihedral::Array{Float64, 1}
    sweep::Array{Float64, 1}
    ds::Array{Float64, 1}
end

struct position
    x::Array{Float64, 1}
    y::Array{Float64, 1}
    z::Array{Float64, 1}
end

struct wingsection
    span
    chord
    twist
    tc
    sweep
    dihedral
    N
end

struct fs_def
    alpha
    CL
    method
end
struct ref_def
    S
    c
    CLmax
end
struct pdrag_def
    polar
    alt
    xt
    method
end
struct mvr_def
    qN
    n
    kbar
end

## -------- geometry ---------------

function geometry(wing::wingsection)

    # ----------------- rename for convenience --------------------
    b = wing.span
    Lambda = wing.sweep
    chord = wing.chord
    twist = wing.twist
    tc = wing.tc
    phi = wing.dihedral
    N = wing.N

    # --------------------------------------------------------------

    P = round.(Int64, b/sum(b)*N) # divide up panels

    # -------------  Quarter Chord Locations --------------------------
    M = 1 + sum(P)

    c = zeros(M)
    t = zeros(M)
    thickness = zeros(M)
    # QC = position(zeros(M), zeros(M), zeros(M))
    # LE = position(zeros(M), zeros(M), zeros(M))
    # TE = position(zeros(M), zeros(M), zeros(M))
    # CP = CPdata(zeros(N), zeros(N), zeros(N), zeros(N), zeros(N), zeros(N),
    #     zeros(N), zeros(N), zeros(N))
    QCx = zeros(M)
    QCy = zeros(M)
    QCz = zeros(M)
    LEx = zeros(M)
    LEy = zeros(M)
    LEz = zeros(M)
    TEx = zeros(M)
    TEy = zeros(M)
    TEz = zeros(M)
    CPx = zeros(N)
    CPy = zeros(N)
    CPz = zeros(N)

    last = 1

    for i = 1:length(b)
        first = last
        last = first + P[i]
        eta = linspace(0, b[i], P[i]+1)

        QCx[first:last] = QCx[first] + eta*tan(Lambda[i])
        QCy[first:last] = QCy[first] + eta*cos(phi[i])
        QCz[first:last] = QCz[first] + eta*sin(phi[i])
        c[first:last] = chord[i] + eta*(chord[i+1]-chord[i])/b[i] # chord
        t[first:last] = twist[i] + eta*(twist[i+1]-twist[i])/b[i] # twist
        thickness[first:last] = tc[i]*chord[i] + eta*(tc[i+1]*chord[i+1]-tc[i]*chord[i])/b[i] # thickness
    end
    QC = position(QCx, QCy, QCz)
    # --------------------------------------------------------------

    # ------------ Trailing Edge Locations --------------------------
    TEx = QCx + 3.0/4*c
    TEy = QCy
    TEz = QCz
    TE = position(TEx, TEy, TEz)
    # ----------------------------------------------------------------

    # -------------- Leading Edge locations ------------------------
    LEx = QCx - c/4
    LEy = QCy
    LEz = QCz
    LE = position(LEx, LEy, LEz)
    # -----------------------------------------------------------------

    # ------------- Control Point Locations --------------------------
    chord = 0.5*(c[1:N] + c[2:N+1])
    twist = 0.5*(t[1:N] + t[2:N+1])
    tc = 0.5*(thickness[1:N] + thickness[2:N+1])./chord
    CPx = 0.5*(QC.x[1:N] + QC.x[2:N+1]) + 0.5*chord
    CPy = 0.5*(QC.y[1:N] + QC.y[2:N+1])
    CPz = 0.5*(QC.z[1:N] + QC.z[2:N+1])

    dihedral = zeros(N)
    sweep = zeros(N)

    last = 0
    for i = 1:length(b)
      first = last + 1
      last = first + P[i] - 1
      dihedral[first:last] = phi[i]
      sweep[first:last] = Lambda[i]
    end
    ds = sqrt.((QC.y[2:N+1]-QC.y[1:N]).^2 + (QC.z[2:N+1]-QC.z[1:N]).^2)

    CP = CPdata(chord, twist, tc, CPx, CPy, CPz, dihedral, sweep, ds)

    return QC, TE, CP, LE

end


## -------- influence coefficients ---------------


#=
--------------------
lift influence coefficients

Author: S. Andrew Ning
---------------------
=#
function getLIC(CP::CPdata, rho::Float64, Vinf::Float64)
    LIC = 2*rho*Vinf*cos.(CP.dihedral).*CP.ds
end

#=
-----------------------------------
Pitching moment coefficients about (x) c.g. location

Author: S. Andrew Ning
Updates: 1/30/09 - removed dependence on QC since this is known
-------------------------------------
=#
function getMIC(CP::CPdata, rho::Float64, Vinf::Float64, xcg::Float64)

    xQC = CP.x - 0.5*CP.chord # quarter chord locations

    MIC = -(xQC - xcg)*2*rho*Vinf.*cos.(CP.dihedral).*CP.ds

end


#=
------------------------------------------------
Induced drag influence coefficients with a Trefftz plane calculation
Di = gamma'*DIC*gamma


Author: S. Andrew Ning

Updated feb. 13 2008
1 is control points, 2 is vorticies.  If you leave out last two
arguments the control points and vorticies are on one wing.  Include them
if you need effect of wing on tail for example.
-------------------------------------------------
=#

function getDIC(yFF::Array{Float64, 1}, zFF::Array{Float64, 1}, rho::Float64)
    getDIC(yFF, zFF, rho, yFF, zFF)
end

function getDIC(yFF::Array{Float64, 1}, zFF::Array{Float64, 1}, rho::Float64, yFF2::Array{Float64, 1}, zFF2::Array{Float64, 1})

    n = length(yFF)-1
    n2 = length(yFF2)-1

     # ------------ define useful variables ------------------
    yc = 0.5*(yFF[2:n+1] + yFF[1:n])   # center of panels
    zc = 0.5*(zFF[2:n+1] + zFF[1:n])
    ny = zFF[2:n+1] - zFF[1:n]
    nz = yFF[2:n+1] - yFF[1:n]
     # ----------------------------------------------------

    DIC = zeros(n,n2)
     #------- normal wash calculation ----------------
    for i = 1:n
        for j = 1:n2
            ry = yFF2[j] - yc[i]
            rz = zFF2[j] - zc[i]
            r2 = ry^2 + rz^2
            DIC[i, j] = -rho/(2*pi*r2)*(rz*ny[i]+ry*nz[i])

             # add other side of panel
            ry = yFF2[j+1] - yc[i]
            rz = zFF2[j+1] - zc[i]
            r2 = ry^2 + rz^2
            DIC[i, j] += rho/(2*pi*r2)*(rz*ny[i]+ry*nz[i])
        end
    end

    yFF2 = -yFF2 # add left side

    for i = 1:n
        for j = 1:n2
            ry = yFF2[j] - yc[i]
            rz = zFF2[j] - zc[i]
            r2 = ry^2 + rz^2
            DIC[i,j] += rho/(2*pi*r2)*(rz*ny[i]+ry*nz[i])

            # add other side of panel
            ry = yFF2[j+1] - yc[i]
            rz = zFF2[j+1] - zc[i]
            r2 = ry^2 + rz^2
            DIC[i,j] -= rho/(2*pi*r2)*(rz*ny[i]+ry*nz[i])
        end
    end

    return DIC
end



#=
---------------------
assumes a parabolic variation in viscous drag with section lift coefficient
Dp = D0 + D1*gamma + D2*gamma.^2

Author: S. Andrew Ning
-----------------------
=#

function getViscousDrag(cd1::Float64, cd2::Float64, CP::CPdata, rho::Float64, Vinf::Float64)

  D1 = cd1*rho*Vinf*CP.ds
  D2 = cd2*2*rho./CP.chord.*CP.ds

  return D1, D2
end


#=
----------------------------------------
 Calculates velocity induced by vortex (with unit vorticity)
 starting at rA ending at rB at control point rC. (see Bertin/Smith)

 Author: S. Andrew Ning
------------------------------------------
=#

function vortex(xA::Array{Float64, 1}, yA::Array{Float64, 1}, zA::Array{Float64, 1},
        xB::Array{Float64, 1}, yB::Array{Float64, 1}, zB::Array{Float64, 1},
        xC::Array{Float64, 1}, yC::Array{Float64, 1}, zC::Array{Float64, 1})

    m = length(yC)
    n = length(yA)
    v = zeros(m, n)
    w = zeros(m, n)

    for i = 1:m
        for j = 1:n
            # define repeatedly used values
            x1 = xC[i] - xA[j]
            x2 = xC[i] - xB[j]
            x21 = xB[j] - xA[j]
            y1 = yC[i] - yA[j]
            y2 = yC[i] - yB[j]
            y21 = yB[j] - yA[j]
            z1 = zC[i] - zA[j]
            z2 = zC[i] - zB[j]
            z21 = zB[j] - zA[j]

            denom1 = (y1*z2-y2*z1)^2 + (x1*z2-x2*z1)^2 + (x1*y2-x2*y1)^2
            # i1 = y1*z2 - y2*z1
            j1 = -x1*z2 + x2*z1
            k1 = x1*y2 - x2*y1
            frac1_1 = (x21*x1 + y21*y1 + z21*z1)/sqrt(x1^2 + y1^2 + z1^2)
            frac1_2 = (x21*x2 + y21*y2 + z21*z2)/sqrt(x2^2 + y2^2 + z2^2)
            frac1 = frac1_1 - frac1_2
            # u[i, j] = 1.0/(4*pi)*i1/denom1*frac1
            v[i, j] = 1.0/(4*pi)*j1/denom1*frac1
            w[i, j] = 1.0/(4*pi)*k1/denom1*frac1
        end
    end

    return v, w
end


#=
------------------------------------------------
This function calculates the induced velocies at control points
xCP,yCP,zCP from the vorticies which are bound at the quarter chord,
follow the wing to the trailing edge, and leave the wing parallel to the
freestream (drag-free wake).  All values are computed for unit vorticity

Author: Andrew Ning
------------------------------------------------
=#
function velocity3d(QC::position, TE::position, CP::CPdata)

    # first calculate induced velocities from vorticies on right side of wing.
    m = length(CP.y)
    n = length(QC.y)-1

    # u = zeros(m, n)
    v = zeros(m, n)
    w = zeros(m, n)

    # Induced velocities due to bound vortex (BC)
    dv, dw = vortex(QC.x[1:n], QC.y[1:n], QC.z[1:n], QC.x[2:n+1], QC.y[2:n+1], QC.z[2:n+1],
        CP.x, CP.y, CP.z)
    # u += du
    v += dv
    w += dw

    # Induced velocities due to AB - left vortex bound to wing
    dv, dw = vortex(TE.x[1:n], TE.y[1:n], TE.z[1:n], QC.x[1:n], QC.y[1:n], QC.z[1:n],
        CP.x, CP.y, CP.z)
    # u += du
    v += dv
    w += dw

    # Induced velocities due to CD - right vortex bound to wing
    dv, dw = vortex(QC.x[2:n+1], QC.y[2:n+1], QC.z[2:n+1], TE.x[2:n+1], TE.y[2:n+1], TE.z[2:n+1],
        CP.x, CP.y, CP.z)
    # u += du
    v += dv
    w += dw

    # Induced velocities due to trailing vortex from A
    l = 1e9  # some large number to represent infinity
    # xINF = TE.x + l*cos(alpha)*ones(size(TE.x))
    # yINF = TE.y
    # zINF = TE.z + l*sin(alpha)*ones(size(TE.z))
    xINF = TE.x + l
    yINF = TE.y
    zINF = TE.z

    dv, dw = vortex(xINF[1:n], yINF[1:n], zINF[1:n], TE.x[1:n], TE.y[1:n], TE.z[1:n],
        CP.x, CP.y, CP.z)
    # u += du
    v += dv
    w += dw

    # Induced velocities due to trailing vortex from D
    dv, dw = vortex(TE.x[2:n+1], TE.y[2:n+1], TE.z[2:n+1], xINF[2:n+1], yINF[2:n+1], zINF[2:n+1],
        CP.x, CP.y, CP.z)
    # u += du
    v += dv
    w += dw

    # ---------- Now repeat with contribution from left side of the wing -----
    # y switched signs (except control points of course)
    # QCy = -QC.y
    # TEy = -TE.y
    # yINF = -yINF
    # Also all induced velocities have opposite sign because vorticies travel
    # in opposite direction.

    # Induced velocities due to bound vortex (BC)
    dv, dw = vortex(QC.x[1:n], -QC.y[1:n], QC.z[1:n], QC.x[2:n+1], -QC.y[2:n+1], QC.z[2:n+1],
        CP.x, CP.y, CP.z)
    # u -= du
    v -= dv
    w -= dw

    # Induced velocities due to AB - left vortex bound to wing
    dv, dw = vortex(TE.x[1:n], -TE.y[1:n], TE.z[1:n], QC.x[1:n], -QC.y[1:n], QC.z[1:n],
        CP.x, CP.y, CP.z)
    # u -= du
    v -= dv
    w -= dw

    # Induced velocities due to CD - right vortex bound to wing
    dv, dw = vortex(QC.x[2:n+1], -QC.y[2:n+1], QC.z[2:n+1], TE.x[2:n+1], -TE.y[2:n+1], TE.z[2:n+1],
        CP.x, CP.y, CP.z)
    # u -= du
    v -= dv
    w -= dw

    # Induced velocities due to trailing vortex from A

    dv, dw = vortex(xINF[1:n], -yINF[1:n], zINF[1:n], TE.x[1:n], -TE.y[1:n], TE.z[1:n],
        CP.x, CP.y, CP.z)
    # u -= du
    v -= dv
    w -= dw

    # Induced velocities due to trailing vortex from D
    dv, dw = vortex(TE.x[2:n+1], -TE.y[2:n+1], TE.z[2:n+1], xINF[2:n+1], -yINF[2:n+1], zINF[2:n+1],
        CP.x, CP.y, CP.z)
    # u -= du
    v -= dv
    w -= dw

    return v, w
end



#=
------------------------------------
Aerodynamic influence coefficients
Vn = AIC*gamma

Author: S. Andrew Ning
Updates: summer 08 - simplified so that twist/camber are only applied in boundary conditions
                     consistent with assumptions of VLM.
--------------------------------------
=#

function getAIC(QC::position, TE::position, CP::CPdata)

    v, w =  velocity3d(QC, TE, CP)

    m = length(CP.y)
    n = length(QC.y)-1

    AIC = zeros(m, n)

    for j = 1:n
        AIC[:, j] = -v[:, j].*sin.(CP.dihedral) + w[:, j].*cos.(CP.dihedral)
    end

    return AIC
end



#=
-------------------------------------------------
this function computes integrated bending moment over thickness and the
bending moment distribution
W = WIC*gamma;
bending moment distribution = BMM*gamma

Author: S. Andrew Ning
Updates: 1/23/09 - corrected integration along spar and vectorized for speed
         1/30/09 - allows for thickness distribution
---------------------------------------------------
=#


function getWIC(CP::CPdata, rho::Float64, Vinf::Float64)

    x = CP.x - 0.5*CP.chord
    y = CP.y
    z = CP.z

    # rename for convenience
    N = length(y)
    phi = CP.dihedral
    Lambda = CP.sweep
    chord = CP.chord
    tc = CP.tc
    ds = CP.ds

    R = zeros(N, N)

    for i = 1:N
        for j = i+1:N
            cij = cos(CP.dihedral[j])*cos(CP.dihedral[i]) + sin(CP.dihedral[j])*sin(CP.dihedral[i])
            R[i, j] = (x[j] - x[i])*cij*sin(CP.sweep[i]) +
                     (y[j] - y[i])*cos(CP.dihedral[j])*cos(CP.sweep[i]) +
                     (z[j] - z[i])*sin(CP.dihedral[j])*cos(CP.sweep[i])

            # if j <= i, R(i,j) = 0
        end
    end


    # integrate along structural span
    ds_str = CP.ds./cos.(CP.sweep)

    BMM = rho*Vinf*R.*repmat(ds.', N, 1)  #(ones(N, 1)*ds)
    WIC = (2*(ds_str./(tc.*chord)).'*BMM).'

    return WIC, BMM
end

# TODO: haven't tested which is faster in Julia

# # compute R matrix
# PhiM = ((cos(phi)*cos(phi)') + (sin(phi)*sin(phi)'))
# Rx = (sin(Lambda)*x' - (x.*sin(Lambda))*ones(1,N)).*PhiM
# Ry = cos(Lambda)*(y.*cos(phi))' - (y.*cos(Lambda))*cos(phi)'
# Rz = cos(Lambda)*(z.*sin(phi))' - (z.*cos(Lambda))*sin(phi)'
# R = Rx + Ry + Rz
# R = R - tril(R) # R(i,j) = 0 for j <= i


# % below version is exactly the same (slower since its in loops, but left
# % here since it's easier to read then vectorized version above)
#
# % R = zeros(N);
#
# % for i = 1:N
# %     for j = 1:N
# %         cij = cos(CP.dihedral(j))*cos(CP.dihedral(i)) + sin(CP.dihedral(j))*sin(CP.dihedral(i));
# %         R(i,j) = (x(j) - x(i))*cij*sin(CP.sweep(i)) + ...
# %                  (y(j) - y(i))*cos(CP.dihedral(j))*cos(CP.sweep(i)) + ...
# %                  (z(j) - z(i))*sin(CP.dihedral(j))*cos(CP.sweep(i));
# %         if j <= i
# %             R(i,j) = 0;
# %         end
# %     end
# % end

# % ds_str = CP.ds./cos(CP.sweep);
#
# % WIC = zeros(1,N);
# %
# % for j = 1:N
# %     WIC(j) = 2*rho*U/tc*CP.ds(j).*ds_str./CP.chord*R(:,j);
# % end

#=
-----------------------------------------------------
This function computes angle of attack necessary to match a specified
lift

Author: S. Andrew Ning
Updates: 1/28/09 - created
         9/2/09 - changed the way alpha was computed - seems more robust
                 to larger angles of attack by using arcsin rather than arccos
---------------------------------------------------
=#

function getAlpha(Ltarget::Float64, CP::CPdata, LIC::Array{Float64, 1}, AIC::Array{Float64, 2}, Vinf::Float64)

    # solve quadratic equation for alpha
    qq = LIC/AIC
    a = -U*qq*sin.(CP.twist)'
    b = -U*qq*(cos.(CP.twist).*cos.(CP.dihedral))'
    c = Ltarget
    # alpha = acos((a*c + b*sqrt(a^2+b^2-c^2))/(a^2+b^2));
    alpha = asin(c/b - a/b*(a*c + b*sqrt(a^2+b^2-c^2))/(a^2 + b^2))

    return alpha
end


## -------- vortex lattice method ---------------


#=
----------------------------------------------------------------------
This version assumes (one) symmetric wing.

Author: S. Andrew Ning
Updates: repackaged with features from several different versions - 1/30/09
        2/2/09 - fixed computation of t/c_avg to be properly area weighted
                (used in compressibility and Re dependent parasite drag)
               - clmax is now a function of t/c (polynomial fit from 241 data)
        2/3/09 - added viscous dependent induced drag if PASS method is used to compute
                 parasite drag
        2/17/09 - added area dependent weight, correct area in quadratic parasite
                  drag calculation, parameter to choose pdrag method
----------------------------------------------------------------------
=#

function VLM(wing, fs, ref, pdrag, mvr, plots)

    # freestream properties
    if fs.method == "CL" # specify CL
        CLref = fs.CL
        alpha = 0.0
    else
        alpha = fs.alpha
        CLref = 0.0
    end

    # mach = fs.mach
    rho = 1.0
    Vinf = 1.0
    q = 0.5*rho*Vinf^2

    # reference quantities
    Sref = ref.S
    cref = ref.c

    # viscous drag (2 possible methods)
    if pdrag.method == "pass" # Reynolds number dependent method
        alt = pdrag.alt
        xt = pdrag.xt

        # not used
        cd0 = 0.0
        cd1 = 0.0
        cd2 = 0.0

    else # quadratic variation with section lift coefficient
        cd0 = pdrag.polar[1]
        cd1 = pdrag.polar[2]
        cd2 = pdrag.polar[3]

        # not used
        alt = 0.0
        xt = 0.0
    end

    # structures
    qmvrN = mvr.qN # ratio of maneuver dynamic pressure to cruise dynamic pressure
    n = mvr.n # load factor
    kbar = mvr.kbar # coefficient used in computing area dependent weight (3.57e4 or 0)

    # stall
    CLmax = ref.CLmax
    # -------------------------------------

    # ------------- geometry ------------------
    QC, TE, CP, LE = geometry(wing)

    # compute wing area
    S = 2*sum(CP.chord.*CP.ds)
    # ------------------------------------------

    # ---------- force influence coefficients --------------------
    # aerodynamic influence coefficients
    AIC = getAIC(QC, TE, CP)

    # lift
    LIC = getLIC(CP, rho, Vinf)

    # induced drag
    DIC = getDIC(TE.y, TE.z, rho)

    # viscous drag
    D1, D2 = getViscousDrag(cd1, cd2, CP, rho, Vinf)

    # weight
    WIC, BMM = getWIC(CP, rho, Vinf)
    # -----------------------------------------------

    # ---------- compute angle of attack necessary to match CL -----
    if fs.method == "CL"
        alpha = getAlpha(CLref*q*Sref, CP, LIC, AIC, Vinf)
    end
    # ------------------------------------------------------------------

    # -------- compute circulation -----------------
    Vn = -Vinf*(cos(alpha)*sin.(CP.twist) + sin(alpha)*cos.(CP.twist).*cos.(CP.dihedral))
    gamma = AIC\Vn
    # ----------------------------------------------

    # --------- aerodynamic forces ------------------
    L = dot(LIC, gamma)
    Di = gamma'*DIC*gamma
    Dp = cd0*q*S + D1'*gamma + D2'*gamma.^2

    CL = L/q/Sref
    CDi = Di/q/Sref
    CDp = Dp/q/Sref
    # ------------------------------------------------

    # --------- weight (integrated bending moment over thickness --------
    # circulation at maneuver load
    bbb = AIC\cos.(CP.dihedral)
    LL = dot(LIC, bbb)
    gamma_mvr = gamma + ((n/qmvrN-1)*(LIC'*gamma)/LL*bbb')'

    # compute weight
    W = qmvrN*WIC'*gamma_mvr

    # add area dependent weight
    W = W + kbar*S

    # weight coefficient
    CW = W/q/Sref/cref
    # -----------------------------------------------------------------

    # # ----- estimate compressibility and parasite drag (strip theory + PASS method) -----
    # supercrit = 1
    # P = round(wing.span/sum(wing.span)*wing.N) # number of panels in each section
    # cdc = zeros(1,length(wing.span))
    # cdp = zeros(1,length(wing.span))
    # area = zeros(1,length(wing.span))
    # start = 1
    #
    # for i = 1:length(wing.span)
    #   finish = start + P[i] - 1
    #
    #   # rename for convenience
    #   cr = wing.chord[i]
    #   ct = wing.chord[i+1]
    #   cbar = 1/2*(cr + ct)
    #   tcbar = (wing.tc[i]*wing.chord[i] + wing.tc[i+1]*wing.chord[i+1])/(wing.chord[i]+wing.chord[i+1])
    #   area[i] = cbar*wing.span[i]
    #   mac = 2/3*(cr + ct - cr*ct/(cr+ct))
    #   CL_local = 0.1*sum(gamma[start:finish]'.*CP.ds[start:finish])*2/U/area[i]
    #
    #   # compressibility drag
    #   cdc[i] = Cdrag(CL_local,wing.sweep[i],tcbar,mach,supercrit)
    #
    #   # parasite drag
    #   cdp[i] = Pdrag(alt,mach,xt,mac,wing.sweep[i],tcbar)
    #
    #   start = finish + 1
    # end
    #
    # # compressibility drag - area weighted average
    # CDc = 2*sum(cdc.*area)/Sref[1]
    #
    # if pass_mthd
    #   # parasite drag - area weighted average
    #   CDp = 2*sum(cdp.*area)/Sref[1]
    #
    #   # add viscous dependent induced drag
    #   Lambda_bar = area*wing.sweep'/sum(area)
    #   CDi = CDi + 0.38*CDp*CL^2/cos(Lambda_bar)^2
    # end
    # # ----------------------------------------------------------------------

    # ----------- cl distribution at CLmax ----------------
    cl = 2/Vinf*gamma./CP.chord
    clmax_dist = cl + rho*Vinf*(CLmax*Sref-L/q)/LL*bbb./CP.chord

    # clmax as a function of thickness - polynomial fit
    tc = CP.tc*100
    clmax = -1.748 + 0.8013*tc - 0.06567*tc.^2 + 0.0022307*tc.^3 - 2.7634e-5*tc.^4

    cl_margin = clmax - clmax_dist
    # ----------------------------------------

    # --------- pitching moment about a.c. --------------------
    # pitching moment about quarter chord
    MIC = getMIC(CP, rho, Vinf, QC.x[1])

    # find aerodynamic center
    dRHS = -sin(alpha)*sin.(CP.twist) + cos(alpha)*cos.(CP.twist).*cos.(CP.dihedral)
    dbc = AIC\dRHS
    dMda = MIC'*dbc
    dLda = LIC'*dbc
    xac = -dMda/dLda + QC.x[1]

    # find Moment about a.c.
    MICac = getMIC(CP, rho, Vinf, xac[1])
    Mac = MICac'*gamma

    Cmac = Mac/q/Sref/cref
    # ----------------------------------------------------------

    # --------------- structures --------------------------
    # bending moment distribution
    Mb = qmvrN*BMM*gamma

    # distance along structural span
    ds_str = CP.ds./cos.(CP.sweep)
    eta_str = [0; cumsum(ds_str, 2)]
    eta_str = 0.5*(eta_str[1:end-1] + eta_str[2:end])
    # --------------------------------------------------------------

    # if (plots)
    #   # ------------- plots --------------------
    #   # plot wings
    # #       plot_wing(LE,QC,TE,CP)
    #   N = length(QC.x)
    #   #       axis equal
    #   PyPlot.figure()
    # #       if TE[3] < 0:
    # #         TE.y = -TE.y
    # #       end
    #   for i = 1:N-1
    #     PyPlot.plot([LE.y[i], LE.y[i+1]], -[LE.x[i], LE.x[i+1]], "b")
    #     PyPlot.plot([TE.y[i], TE.y[i+1]],-[TE.x[i], TE.x[i+1]], "b")
    #   end
    #   PyPlot.plot([LE.y[1], TE.y[1]],-[LE.x[1], TE.x[1]], "b")
    # #       PyPlot.plot([LE.y[end], TE.y[end]],-[LE.x[end], TE.x[end]], "b")
    #   for i = 1:N-1
    #       PyPlot.plot([-LE.y[i], -LE.y[i+1]],-[LE.x[i], LE.x[i+1]],"b")
    #       PyPlot.plot([-TE.y[i], -TE.y[i+1]],-[TE.x[i], TE.x[i+1]],"b")
    #   end
    #   PyPlot.plot([-LE.y[1], -TE.y[1]],-[LE.x[1], TE.x[1]], "b")
    # #       PyPlot.plot([-LE.y[end], -TE.y[end]],-[LE.x[end], TE.x[end]],"b")
    #   PyPlot.xlabel("x")
    #   PyPlot.ylabel("y")
    #   PyPlot.title("Plot of wing")
    #
    #   # plot lift
    #   #       figure(50) hold on
    #   PyPlot.figure()
    #   l = 2*gamma/U/(ref.c)
    #   l_mvr = 2*gamma_mvr/U/(ref.c)
    #   eta = linspace(0,0.5,length(l))
    #   PyPlot.plot(eta,l,"b")
    #   PyPlot.plot(eta,l_mvr,"r")
    #   PyPlot.xlabel("xi / b")
    #   PyPlot.ylabel("c_l c / c_ref")
    # #   PyPlot.title("")
    #
    #   # plot cl
    #   PyPlot.figure()
    #   PyPlot.plot(eta,cl,"b")
    #   PyPlot.plot(eta,clmax_dist,"r")
    #   PyPlot.plot(eta,clmax,"r--")
    #   PyPlot.xlabel("xi / b")
    #   PyPlot.ylabel("c_l")
    #   PyPlot.title("Plot of c_l")
    #
    #   # plot bending over thickness
    #   PyPlot.figure()
    #   PyPlot.plot(eta_str/eta_str[end]*0.5,Mb./(CP.tc.*CP.chord))
    #   PyPlot.xlabel("xi_str/b_str")
    #   PyPlot.ylabel("M_b/t")
    #   PyPlot.title("Plot of bending over thickness")
    #   PyPlot.show()
    #   # -------------------------------------------------
    # end

    return CL, CDi, CDp, CW, Cmac, cl_margin, gamma, CP, cl
end

end
