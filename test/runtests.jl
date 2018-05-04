using Base.Test
import VLM

# tests using AVL 3.35

# cat commands - | ./avl3.35 simplewing1.avl 

# ---- run1 ------

symmetric = true
halfspan = 7.5
b = 2*halfspan
cr = 2.2
ct = 1.8
AR = b^2/(b*0.5*(cr + ct))
λ = ct/cr
Λ = atan((0.4 + ct/4 - cr/4)/halfspan)
ϕ = 0.0
θr = 2.0*pi/180  
θt = 2.0*pi/180  
npanels = 11
spacing = "uniform"
Sref = 30.0
cref = 2.0
bref = 15.0
geom = VLM.simpletapered(b, AR, λ, Λ, ϕ, θr, θt, Sref, cref, bref, npanels, symmetric, spacing)

rho = 1.0
Vinf = 2.0
alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
rcg = [0.50 - cr/4, 0.0, 0.0]
vother = nothing
fs = VLM.Freestream(rho, Vinf, alpha, beta, Omega, rcg, vother)

CF, CM, ymid, zmid, l, cl, dCF = VLM.run(geom, fs, symmetric)
CD, CY, CL = CF
Cl, Cm, Cn = CM
@test isapprox(CL, 0.24324, atol=1e-3)
@test isapprox(CD, 0.00245, atol=1e-5)
@test isapprox(Cm, -0.02252, atol=1e-4)
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0

# -----------------------------

# ---- run2 ------
spacing = "cosine"
geom = VLM.simpletapered(b, AR, λ, Λ, ϕ, θr, θt, npanels, symmetric, spacing)

CF, CM, ymid, zmid, l, cl, dCF = VLM.run(geom, fs, symmetric)
CD, CY, CL = CF
Cl, Cm, Cn = CM
@test isapprox(CL, 0.23744, atol=1e-3)
@test isapprox(CD, 0.00243, atol=1e-5)
@test isapprox(Cm, -0.02165, atol=1e-4)
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0

# TODO: multi-component aircraft
# TODO: cosine spacing
# TODO: set reference values


# ---- run3 ------

spacing = "uniform"
geom = VLM.simpletapered(b, AR, λ, Λ, ϕ, θr, θt, Sref, cref, bref, npanels, symmetric, spacing)

alpha = 8.0*pi/180
fs = VLM.Freestream(rho, Vinf, alpha, beta, Omega, rcg, vother)

CF, CM, ymid, zmid, l, cl, dCF = VLM.run(geom, fs, symmetric)
CD, CY, CL = CF
Cl, Cm, Cn = CM
@test isapprox(CL, 0.80348, atol=2e-3)
@test isapprox(CD, 0.02696, atol=2e-5)
@test isapprox(Cm, -0.07399, atol=2.5e-4)
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0

# -----------------------------


# ------ run4 -----
symmetric = true
halfspan = 7.5
height = 3.0
b = 2*(sqrt(halfspan^2 + height^2))
cr = 2.2
ct = 1.8
AR = b^2/(b*0.5*(cr + ct))
λ = ct/cr
Λ = atan((0.4 + ct/4 - cr/4)/halfspan)
ϕ = atan(height/halfspan)
θr = 2.0*pi/180  
θt = 2.0*pi/180  
npanels = 11
spacing = "uniform"
Sref = 30.0
cref = 2.0
bref = 15.0
geom = VLM.simpletapered(b, AR, λ, Λ, ϕ, θr, θt, Sref, cref, bref, npanels, symmetric, spacing)

rho = 1.0
Vinf = 2.0
alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
rcg = [0.50 - cr/4, 0.0, 0.0]
vother = nothing
fs = VLM.Freestream(rho, Vinf, alpha, beta, Omega, rcg, vother)

CF, CM, ymid, zmid, l, cl, dCF = VLM.run(geom, fs, symmetric)
CD, CY, CL = CF
Cl, Cm, Cn = CM
@test isapprox(CL, 0.24787, atol=1e-3)
@test isapprox(CD, 0.00245, atol=1e-5)
@test isapprox(Cm, -0.02395, atol=2e-3)
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0


# ------ run5 -----

alpha = 20.0*pi/180  # nonphysical, just testing the numerics
fs = VLM.Freestream(rho, Vinf, alpha, beta, Omega, rcg, vother)

CF, CM, ymid, zmid, l, cl, dCF = VLM.run(geom, fs, symmetric)
CD, CY, CL = CF
Cl, Cm, Cn = CM
@test isapprox(CL, 1.70985, atol=.04*abs(CL))  # this error seems higher than I'd expect...
@test isapprox(CD, 0.11502, atol=.001*abs(CD))
@test isapprox(Cm, -0.45606, atol=.02*abs(Cm))
@test CY == 0.0
@test Cl == 0.0
@test Cn == 0.0

