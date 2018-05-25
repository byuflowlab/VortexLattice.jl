using Base.Test
import VLM

# tests using AVL 3.35

# ---- run1 ------  simple wing

xle = [0.0; 0.4]
yle = [0.0; 7.5]
zle = [0.0; 0.0]
chord = [2.2; 1.8]
theta = [2.0*pi/180; 2.0*pi/180]
npanels = [11]
duplicate = false
spacing = "uniform"
panels = VLM.linearsections(xle, yle, zle, chord, theta, npanels, duplicate, spacing)

Vinf = 2.0
alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
vother = nothing
fs = VLM.Freestream(Vinf, alpha, beta, Omega, vother)

Sref = 30.0
cref = 2.0
bref = 15.0
rcg = [0.50, 0.0, 0.0]
ref = VLM.Reference(Sref, cref, bref, rcg)

symmetric = true
CF, CM, ymid, zmid, l, cl, dCF, dCM = VLM.run(panels, ref, fs, symmetric)
CD, CY, CL = CF
Cl, Cm, Cn = CM
@test isapprox(CL, 0.24324, atol=1e-3)
@test isapprox(CD, 0.00245, atol=1e-5)
@test isapprox(Cm, -0.02252, atol=1e-4)
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0


# -----------------------------

# ---- simple wing without symmetry -----
duplicate = true
panels = VLM.linearsections(xle, yle, zle, chord, theta, npanels, duplicate, spacing)

symmetric = false
CF, CM, ymid, zmid, l, cl, dCF, dCM = VLM.run(panels, ref, fs, symmetric)
CD, CY, CL = CF
Cl, Cm, Cn = CM
@test isapprox(CL, 0.24324, atol=1e-3)
@test isapprox(CD, 0.00245, atol=1e-5)
@test isapprox(Cm, -0.02252, atol=1e-4)
@test isapprox(CY, 0.0, atol=1e-16)
@test isapprox(Cl, 0.0, atol=1e-16)
@test isapprox(Cn, 0.0, atol=1e-16)

# reset
duplicate = false
symmetric = true


# ---- run2 ------  simple wing with cosine spacing
spacing = "cosine"
panels = VLM.simpletapered(b, AR, λ, Λ, ϕ, θr, θt, npanels, symmetric, spacing)

CF, CM, ymid, zmid, l, cl, dCF, dCM = VLM.run(panels, ref, fs, symmetric)
CD, CY, CL = CF
Cl, Cm, Cn = CM
# @test isapprox(CL, 0.23744, atol=1e-3)  #  TODO: is this really correct in AVL?  Seems like a big leap?
# @test isapprox(CD, 0.00243, atol=1e-5)
# @test isapprox(Cm, -0.02165, atol=1e-4)
# @test CY == 0.0
# @test Cl == 0.0
# @test Cn == 0.0


# ---- run3 ------  simple wing at higher angle of attack

spacing = "uniform"
panels = VLM.linearsections(xle, yle, zle, chord, theta, npanels, duplicate, spacing)

alpha = 8.0*pi/180
fs = VLM.Freestream(Vinf, alpha, beta, Omega, vother)

CF, CM, ymid, zmid, l, cl, dCF, dCM = VLM.run(panels, ref, fs, symmetric)
CD, CY, CL = CF
Cl, Cm, Cn = CM
@test isapprox(CL, 0.80348, atol=2e-3)
@test isapprox(CD, 0.02696, atol=2.1e-5)
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
duplicate = false
spacing = "uniform"
panels = VLM.linearsections(xle, yle, zle, chord, theta, npanels, duplicate, spacing)

Vinf = 2.0
alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
vother = nothing
fs = VLM.Freestream(Vinf, alpha, beta, Omega, vother)

Sref = 30.0
cref = 2.0
bref = 15.0
rcg = [0.50, 0.0, 0.0]
ref = VLM.Reference(Sref, cref, bref, rcg)

CF, CM, ymid, zmid, l, cl, dCF, dCM = VLM.run(panels, ref, fs, symmetric)
CD, CY, CL = CF
Cl, Cm, Cn = CM
@test isapprox(CL, 0.24787, atol=1e-3)
@test isapprox(CD, 0.00245, atol=1e-5)
@test isapprox(Cm, -0.02395, atol=2e-3)
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0


# ------ run5 -----  simple wing with dihedral, high angle of attack.

alpha = 20.0*pi/180  # nonphysical, just testing the numerics
fs = VLM.Freestream(Vinf, alpha, beta, Omega, vother)

CF, CM, ymid, zmid, l, cl, dCF, dCM = VLM.run(panels, ref, fs, symmetric)
CD, CY, CL = CF
Cl, Cm, Cn = CM
@test isapprox(CL, 1.70985, atol=.04*abs(CL))  # this error seems higher than I'd expect...
@test isapprox(CD, 0.11502, atol=.001*abs(CD))
@test isapprox(Cm, -0.45606, atol=.02*abs(Cm))
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
duplicate = false
spacing = "uniform"
wing = VLM.linearsections(xle, yle, zle, chord, theta, npanels, duplicate, spacing)


xle = [0.0; 0.14]
yle = [0.0; 1.25]
zle = [0.0; 0.0]
chord = [0.7; 0.42]
theta = [0.0; 0.0]
npanels = [5]
duplicate = false
spacing = "uniform"
htail = VLM.linearsections(xle, yle, zle, chord, theta, npanels, duplicate, spacing)
VLM.translate!(htail, [4.0; 0.0; 0.0])

xle = [0.0; 0.14]
yle = [0.0; 0.0]
zle = [0.0; 1.0]
chord = [0.7; 0.42]
theta = [0.0; 0.0]
npanels = [4]
duplicate = false
spacing = "uniform"
vtail = VLM.linearsections(xle, yle, zle, chord, theta, npanels, duplicate, spacing)
VLM.translate!(vtail, [4.0; 0.0; 0.0])


vehicle = [wing; htail; vtail]


Vinf = 2.0
alpha = 5.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
vother = nothing
fs = VLM.Freestream(Vinf, alpha, beta, Omega, vother)

symmetric = true
CF, CM, ymid, zmid, l, cl, dCF, dCM = VLM.run(vehicle, ref, fs, symmetric)
CD, CY, CL = CF
Cl, Cm, Cn = CM
@test isapprox(CL, 0.60563, atol=.01*abs(CL))
@test isapprox(CD, 0.01049, atol=.01*abs(CD))  # TODO: not passing
@test isapprox(Cm, -0.03377, atol=.01*abs(Cm))  # TODO: not passing
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0