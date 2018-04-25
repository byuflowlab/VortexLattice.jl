using Base.Test
import VLM


symmetric = true
halfspan = 7.5
b = 2*halfspan
cr = 2.2
ct = 1.8
AR = b^2/(b*0.5*(cr + ct))
λ = ct/cr
Λ = atan((0.4 + ct/4 - cr/4)/halfspan)
ϕ = 0.0
θt = 0.0  
npanels = 11
geom = VLM.simpletapered(b, AR, λ, Λ, ϕ, θt, npanels, symmetric)
# geom.Sref = 30.0
# geom.cref = 2.0
# geom.bref = 15.0

rho = 1.0
Vinf = 2.0
alpha = 3.0*pi/180
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

#  Alpha =   1.00000     pb/2V =   0.00000     p'b/2V =   0.00000
#   Beta  =   0.00000     qc/2V =   0.00000
#   Mach  =     0.000     rb/2V =   0.00000     r'b/2V =   0.00000

#   CXtot =   0.00181     Cltot =   0.00000     Cl'tot =   0.00000
#   CYtot =   0.00000     Cmtot =  -0.02252
#   CZtot =  -0.24325     Cntot =   0.00000     Cn'tot =   0.00000

#   CLtot =   0.24324
#   CDtot =   0.00243
#   CDvis =   0.00000     CDind =   0.00243
#   CLff  =   0.24328     CDff  =   0.00245    | Trefftz
#   CYff  =   0.00000         e =    1.0257    | Plane  

