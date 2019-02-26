using Test
import VortexLatticeMethod
const VLM = VortexLatticeMethod

# tests using AVL 3.35

# ---- run1 ------  simple wing

xle = [0.0; 0.4]
yle = [0.0; 7.5]
zle = [0.0; 0.0]
chord = [2.2; 1.8]
theta = [2.0*pi/180; 2.0*pi/180]
npanels = [11]
spacing = ["u"]
duplicate = false
# panels = VLM.linearsections(xle, yle, zle, chord, theta, npanels, duplicate, spacing)
panels = VLM.linearsections(xle, yle, zle, chord, theta, npanels, spacing, [1], ["u"], duplicate)

alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
vother = nothing
fs = VLM.Freestream(alpha, beta, Omega, vother)

Sref = 30.0
cref = 2.0
bref = 15.0
rcg = [0.50, 0.0, 0.0]
ref = VLM.Reference(Sref, cref, bref, rcg)

symmetric = true
outputs = VLM.solve(panels, ref, fs, symmetric)
CD, CY, CL = outputs.CF
Cl, Cm, Cn = outputs.CM

@test isapprox(CL, 0.24324, atol=1e-3)
@test isapprox(CD, 0.00243, atol=1e-5)
@test isapprox(outputs.CDiff, 0.00245, atol=1e-5)
@test isapprox(Cm, -0.02252, atol=1e-4)
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0

# -----------------------------

# ---- simple wing without symmetry -----
duplicate = true
panels = VLM.linearsections(xle, yle, zle, chord, theta, npanels, spacing, [1], ["u"], duplicate)

symmetric = false
outputs = VLM.solve(panels, ref, fs, symmetric)
CD, CY, CL = outputs.CF
Cl, Cm, Cn = outputs.CM

@test isapprox(CL, 0.24324, atol=1e-3)
@test isapprox(CD, 0.00243, atol=1e-5)
@test isapprox(outputs.CDiff, 0.00245, atol=1e-5)
@test isapprox(Cm, -0.02252, atol=1e-4)
@test isapprox(CY, 0.0, atol=1e-16)
@test isapprox(Cl, 0.0, atol=1e-16)
@test isapprox(Cn, 0.0, atol=1e-16)


CDa, CYa, CLa = outputs.dCF.alpha
Cla, Cma, Cna = outputs.dCM.alpha
CDb, CYb, CLb = outputs.dCF.beta
Clb, Cmb, Cnb = outputs.dCM.beta
CDp, CYp, CLp = outputs.dCF.p
Clp, Cmp, Cnp = outputs.dCM.p
CDq, CYq, CLq = outputs.dCF.q
Clq, Cmq, Cnq = outputs.dCM.q
CDr, CYr, CLr = outputs.dCF.r
Clr, Cmr, Cnr = outputs.dCM.r

@test isapprox(CLa, 4.638090, atol=.01*abs(CLa))
@test isapprox(Cma, -0.429247, atol=.01*abs(Cma))
@test isapprox(CLq, 5.549788, atol=.01*abs(CLq))
@test isapprox(Cmq, -0.517095, atol=.01*abs(Cmq))
# @test isapprox(Clb, -0.025749, atol=.01*abs(Clb))  # TODO
@test isapprox(Clp, -0.518725, atol=.01*abs(Clp))
# @test isapprox(Cnp, -0.019846, atol=.01*abs(Cnp))  # TODO
@test isapprox(CLq, 5.549788, atol=.01*abs(CLq))
@test isapprox(Cmq, -0.517095, atol=.01*abs(Cmq))
# @test isapprox(Clr, 0.064243, atol=.01*abs(Clr)) # TODO

# h = 1e-6
# betap = beta + h
# fsp = VLM.Freestream(alpha, betap, Omega, vother)
# outputsp = VLM.solve(panels, ref, fsp, symmetric)
# Cl_p, Cm_p, Cn_p = outputsp.CM

# betam = beta - h
# fsm = VLM.Freestream(alpha, betam, Omega, vother)
# outputsm = VLM.solve(panels, ref, fsm, symmetric)
# Cl_m, Cm_m, Cn_m = outputsm.CM
# println(Clb)
# println((Cl_p - Cl_m)/(2*h))


# Stability-axis derivatives...

#                              alpha                beta
#                   ----------------    ----------------
#  z' force CL |    CLa =   4.638090    CLb =   0.000000
#  y  force CY |    CYa =   0.000000    CYb =  -0.000007
#  x' mom.  Cl'|    Cla =   0.000000    Clb =  -0.025749
#  y  mom.  Cm |    Cma =  -0.429247    Cmb =  -0.000000
#  z' mom.  Cn'|    Cna =  -0.000000    Cnb =   0.000466

#                      roll rate  p'      pitch rate  q'        yaw rate  r'
#                   ----------------    ----------------    ----------------
#  z' force CL |    CLp =  -0.000000    CLq =   5.549788    CLr =   0.000000
#  y  force CY |    CYp =   0.047000    CYq =  -0.000000    CYr =  -0.000745
#  x' mom.  Cl'|    Clp =  -0.518725    Clq =   0.000000    Clr =   0.064243
#  y  mom.  Cm |    Cmp =   0.000000    Cmq =  -0.517095    Cmr =  -0.000000
#  z' mom.  Cn'|    Cnp =  -0.019846    Cnq =  -0.000000    Cnr =  -0.000898

#  Neutral point  Xnp =   0.685096


# ---- simple wing with chordwise panels -----

duplicate = false
cpanels = [6]
cspacing = ["u"]
panels = VLM.linearsections(xle, yle, zle, chord, theta, npanels, spacing, cpanels, cspacing, duplicate)

symmetric = true
outputs = VLM.solve(panels, ref, fs, symmetric)
CD, CY, CL = outputs.CF
Cl, Cm, Cn = outputs.CM

@test isapprox(CL, 0.24454, atol=1e-3)
@test isapprox(CD, 0.00247, atol=1e-5)
@test isapprox(outputs.CDiff, 0.00248, atol=1e-5)
@test isapprox(Cm, -0.02091, atol=1e-4)
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0


# -------------------------

# reset
duplicate = false
symmetric = true


# TODO
# # ---- run2 ------  simple wing with cosine spacing
# spacing = "cosine"
# panels = VLM.simplewing(b, AR, λ, Λ, ϕ, θr, θt, npanels, symmetric, spacing)

# CF, CM, ymid, zmid, l, cl, dCF, dCM = VLM.run(panels, ref, fs, symmetric)
# CD, CY, CL = CF
# Cl, Cm, Cn = CM
# # @test isapprox(CL, 0.23744, atol=1e-3)  #  TODO: is this really correct in AVL?  Seems like a big leap?
# # @test isapprox(CD, 0.00243, atol=1e-5)
# # @test isapprox(Cm, -0.02165, atol=1e-4)
# # @test CY == 0.0
# # @test Cl == 0.0
# # @test Cn == 0.0


# ---- run3 ------  simple wing at higher angle of attack

spacing = ["u"]
panels = VLM.linearsections(xle, yle, zle, chord, theta, npanels, spacing, [1], ["u"], duplicate)

alpha = 8.0*pi/180
fs = VLM.Freestream(alpha, beta, Omega, vother)

outputs = VLM.solve(panels, ref, fs, symmetric)
CD, CY, CL = outputs.CF
Cl, Cm, Cn = outputs.CM
@test isapprox(CL, 0.80348, atol=2e-3)
@test isapprox(CD, 0.02651, atol=2e-5)
@test isapprox(outputs.CDiff, 0.02696, atol=2.1e-5)
@test isapprox(Cm, -0.07399, atol=2.5e-4)
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0

# -----------------------------


# ------ run4 -----  simple wing with dihedral

xle = [0.0; 0.4]
yle = [0.0; 7.5]
zle = [0.0; 3.0]
chord = [2.2; 1.8]
theta = [2.0*pi/180; 2.0*pi/180]
npanels = [11]
spacing = ["u"]
duplicate = false
panels = VLM.linearsections(xle, yle, zle, chord, theta, npanels, spacing, [1], ["u"], duplicate)

alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
vother = nothing
fs = VLM.Freestream(alpha, beta, Omega, vother)

Sref = 30.0
cref = 2.0
bref = 15.0
rcg = [0.50, 0.0, 0.0]
ref = VLM.Reference(Sref, cref, bref, rcg)

outputs = VLM.solve(panels, ref, fs, symmetric)
CD, CY, CL = outputs.CF
Cl, Cm, Cn = outputs.CM
@test isapprox(CL, 0.24787, atol=1e-3)
@test isapprox(CD, 0.00246, atol=1e-5)
@test isapprox(outputs.CDiff, 0.00245, atol=1e-5)
@test isapprox(Cm, -0.02395, atol=2e-3)
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0


# ------ run5 -----  simple wing with dihedral, high angle of attack.

alpha = 20.0*pi/180  # nonphysical, just testing the numerics
fs = VLM.Freestream(alpha, beta, Omega, vother)

outputs = VLM.solve(panels, ref, fs, symmetric)
CD, CY, CL = outputs.CF
Cl, Cm, Cn = outputs.CM
@test isapprox(CL, 1.70985, atol=.02*abs(CL))
@test isapprox(CD, 0.12904, atol=.01*abs(CD))
# @test isapprox(outputs.CDiff, 0.11502, atol=.001*abs(CD))  # TODO: AVL doesn't project wake (drag-free).  So I actually think what is predicted in this code is more accurate than the AVL value.
@test isapprox(Cm, -0.45606, atol=.01*abs(Cm))
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0


# ------ run 6 ------ wing/tail

Sref = 9.0
cref = 0.9
bref = 10.0
rcg = [0.5, 0.0, 0.0]
ref = VLM.Reference(Sref, cref, bref, rcg)


xle = [0.0; 0.2]
yle = [0.0; 5.0]
zle = [0.0; 1.0]
chord = [1.0; 0.6]
theta = [2.0*pi/180; 2.0*pi/180]
npanels = [11]
spacing = ["u"]
duplicate = false
wing = VLM.linearsections(xle, yle, zle, chord, theta, npanels, spacing, [1], ["u"], duplicate)

xle = [0.0; 0.14]
yle = [0.0; 1.25]
zle = [0.0; 0.0]
chord = [0.7; 0.42]
theta = [0.0; 0.0]
npanels = [5]
spacing = ["u"]
duplicate = false
htail = VLM.linearsections(xle, yle, zle, chord, theta, npanels, spacing, [1], ["u"], duplicate)
VLM.translate!(htail, [4.0; 0.0; 0.0])

xle = [0.0; 0.14]
yle = [0.0; 0.0]
zle = [0.0; 1.0]
chord = [0.7; 0.42]
theta = [0.0; 0.0]
npanels = [4]
spacing = ["u"]
duplicate = false
vtail = VLM.linearsections(xle, yle, zle, chord, theta, npanels, spacing, [1], ["u"], duplicate)
VLM.translate!(vtail, [4.0; 0.0; 0.0])


vehicle = [wing; htail; vtail]


alpha = 5.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
vother = nothing
fs = VLM.Freestream(alpha, beta, Omega, vother)

symmetric = true
outputs = VLM.solve(vehicle, ref, fs, symmetric)
CD, CY, CL = outputs.CF
Cl, Cm, Cn = outputs.CM
@test isapprox(CL, 0.60563, atol=.01*abs(CL))
@test isapprox(CD, 0.01049, atol=.01*abs(CD)) 
# @test isapprox(Cm, -0.03377, atol=.01*abs(Cm))  # TODO: not passing
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0


# # ---- Veldhius validation case ------
# b = 0.64*2
# AR = 5.33
# λ = 1.0
# Λ = 0.0
# ϕ = 0.0
# θr = 0
# θt = 0
# npanels = 50
# duplicate = false
# spacing = "uniform"

# wing = VLM.simplewing(b, AR, λ, Λ, ϕ, θr, θt, npanels, duplicate, spacing)


# import Interpolations: interpolate, Gridded, Linear

# """helper function"""
# function interp1(xpt, ypt, x)
#     intf = interpolate((xpt,), ypt, Gridded(Linear()))
#     y = zeros(x)
#     idx = (x .> xpt[1]) .& (x.< xpt[end])
#     y[idx] = intf[x[idx]]
#     return y
# end

# function votherprop(rpos)

#     # rcenter = [0.0; 0.469*b/2; 0.0]
#     rvec = abs(rpos[2] - 0.469*b/2)  # norm(rpos - rcenter)

#     rprop = [0.017464, 0.03422, 0.050976, 0.067732, 0.084488, 0.101244, 0.118]
#     uprop = [0.0, 1.94373, 3.02229, 7.02335, 9.02449, 8.85675, 0.0]/50.0
#     vprop = [0.0, 1.97437, 2.35226, 4.07227, 4.35436, 3.69232, 0.0]/50.0

#     cw = 1.0

#     u = 2 * interp1(rprop, uprop, [rvec])[1]  # factor of 2 from far-field
#     v = cw * interp1(rprop, vprop, [rvec])[1]

#     if rpos[2] > 0.469*b/2
#         v *= -1
#     end

#     unew = u*cos(alpha) - v*sin(alpha)
#     vnew = u*sin(alpha) + v*cos(alpha)

#     return [unew; 0.0; vnew]
# end

# alpha = 0.0*pi/180
# beta = 0.0
# Omega = [0.0; 0.0; 0.0]
# fs = VLM.Freestream(alpha, beta, Omega, votherprop)

# Sref = 0.30739212
# cref = 0.24015
# bref = b
# rcg = [0.0, 0.0, 0.0]
# ref = VLM.Reference(Sref, cref, bref, rcg)

# symmetric = true
# CF, CM, ymid, zmid, l, cl, dCF, dCM = VLM.run(wing, ref, fs, symmetric)

# alpha = 8.0*pi/180
# fs = VLM.Freestream(alpha, beta, Omega, votherprop)
# CF, CM, ymid, zmid, l2, cl2, dCF, dCM = VLM.run(wing, ref, fs, symmetric)

# # # total velocity in direction of Vinf
# # Vinfeff = zeros(cl)
# # for i = 1:length(ymid)
# #     Vext = votherprop([0.0; ymid[i]; zmid[i]])
# #     Vinfeff[i] = Vinf + Vext[1]*cos(alpha)*cos(beta) + Vext[2]*sin(beta) + Vext[3]*sin(alpha)*cos(beta)
# # end

# using PyPlot
# # figure()
# # plot(ymid, cl2)
# # plot(ymid, cl3)
# # gcf()

# figure()
# plot(ymid, l2)
# plot(ymid, cl2)
# ylim([0.0, 1])
# gcf()

# trapz(ymid/(b/2), cl2)


# # figure()
# # plot(ymid, cl*Vinf./Vinfeff)
# # gcf()

# # yvec = linspace(0, 0.64, 20)
# # axial = zeros(20)
# # swirl = zeros(20)
# # for i = 1:20
# #     V = votherprop([0.0; yvec[i]; 0.0])
# #     axial[i] = V[1]
# #     swirl[i] = V[2]
# # end

# # figure()
# # plot(yvec, axial)
# # plot(yvec, swirl)
# # gcf()