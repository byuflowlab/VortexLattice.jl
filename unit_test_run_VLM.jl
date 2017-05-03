#=
----------------------------------------------------------------------
 This version assumes (one) symmetric wing.

 Sample Inputs Below

 Author: S. Andrew Ning
----------------------------------------------------------------------
=#

include("VLM.jl")
using NVLM
close("all") # KRM include this in case you haven't in your profile

#TEST WITH ORIGINAL SETUP - non-variable twist distribution, and quad viscous drag etc
# wing geometry  - all degrees in radians
wing = NVLM.wingsection([29, 65], [43, 26, 11], [0, 0, 0]*pi/180, [0.13, 0.12, 0.11], [25, 30]*pi/180, [0, 0], 50)
# wing.span = [29 65]
# wing.chord = [43 26 11]
# wing.twist = [0 0 0]*pi/180
# wing.tc = [0.13 0.12 0.11]
# wing.sweep = [25 30]*pi/180
# wing.dihedral = [0 0]
# wing.N = 50

# println(wing.span)
# println(wing.chord)
# println(wing.chord[1:end-1]+wing.chord[2:end])
Sref = dot(wing.chord[1:end-1]+wing.chord[2:end], wing.span)
bref = 2*sum(wing.span.*cos(wing.dihedral))

# freestream parameters - 2 methods either specify CL or angle of attack
fs = NVLM.fs_def(0.01, 5*pi/180, 0.5, "alpha")
# fs.mach = 0.78
# fs.alpha = 5*pi/180 # only used for alpha method
# fs.CL = 0.5 # only used for CL method
# fs.method = "alpha"

# reference and other parameters (used for force/moment coefficients)
ref = NVLM.ref_def(Sref, Sref/bref, 1.4)
# ref.S = Sref
# ref.c = Sref/bref
# ref.CLmax = 1.4

# parasite drag - 2 methods either PASS method or strip theory with a quadratic varition in cl
pdrag = NVLM.pdrag_def([0.007 0 0.003], 35000.0, 0.05, "quad")
# pdrag.polar = [0.007 0 0.003] # only used in quad method
# pdrag.alt = 35000 # only used in pass method
# pdrag.xt = 0.05 # only used in pass method
# pdrag.method = "pass"

# structures
mvr = NVLM.mvr_def(2.5, 2.5, 0)
# mvr.qN = 2.5 # ratio of mvr dynamic pressure to cruise dynamic pressure
# mvr.n = 2.5 # limit load factor
# mvr.kbar = 0 # coefficient for area-dependent weight

# QC = QC_def([0],[0],[0])
# TE = TE_def([0],[0],[0])
# LE = LE_def([0],[0],[0])
# CP = CP_def([0],[0],[0],[0],[0],[0],[0],[0],[0])

plots = false
CL1, CDi1, CDp1, CDc1, CW1, Cmac1, cl_margin1, gamma1, CP1 = NVLM.VLM(wing, fs, ref, pdrag, mvr, plots)
# println(CL)
# println(CDi)
# println(CDp)
# println(CDc)
# println(Cmac)
d1 = ["TEST WITH ORIGINAL SETUP - non-variable twist distribution, and quad viscous drag etc \n \n CL \n $CL1 \n  CDi \n $CDi1 \n  CDp \n $CDp1 \n  CDc \n $CDc1 \n  CW \n $CW1 \n  Cmac \n $Cmac1 \n  cl_margin \n $cl_margin1 \n  gamma \n $gamma1 \n \n \n \n"]


#TEST WITH MODIFIED SETUP - variable twist distribution, and pass viscous drag etc

# wing geometry  - all degrees in radians
wing = NVLM.wingsection([29, 65], [43, 26, 11], [0, 0, 0]*pi/180, [0.13, 0.12, 0.11], [25, 30]*pi/180, [0, 0], 50)
# wing.span = [29 65]
# wing.chord = [43 26 11]
# wing.twist = [0 0 0]*pi/180
# wing.tc = [0.13 0.12 0.11]
# wing.sweep = [25 30]*pi/180
# wing.dihedral = [0 0]
# wing.N = 50
wing.twist = zeros(wing.N)

# println(wing.span)
# println(wing.chord)
# println(wing.chord[1:end-1]+wing.chord[2:end])
Sref = dot(wing.chord[1:end-1]+wing.chord[2:end], wing.span)
bref = 2*sum(wing.span.*cos(wing.dihedral))

# freestream parameters - 2 methods either specify CL or angle of attack
fs = NVLM.fs_def(0.01, 5*pi/180, 0.5, "alpha")
# fs.mach = 0.78
# fs.alpha = 5*pi/180 # only used for alpha method
# fs.CL = 0.5 # only used for CL method
# fs.method = "alpha"
fs.mach = ones(wing.twist)*.01

# reference and other parameters (used for force/moment coefficients)
ref = NVLM.ref_def(Sref, Sref/bref, 1.4)
# ref.S = Sref
# ref.c = Sref/bref
# ref.CLmax = 1.4

# parasite drag - 2 methods either PASS method or strip theory with a quadratic varition in cl
pdrag = NVLM.pdrag_def([0.007 0 0.003], 35000.0, 0.05, "pass")
# pdrag.polar = [0.007 0 0.003] # only used in quad method
# pdrag.alt = 35000 # only used in pass method
# pdrag.xt = 0.05 # only used in pass method
# pdrag.method = "pass"

# structures
mvr = NVLM.mvr_def(2.5, 2.5, 0)
# mvr.qN = 2.5 # ratio of mvr dynamic pressure to cruise dynamic pressure
# mvr.n = 2.5 # limit load factor
# mvr.kbar = 0 # coefficient for area-dependent weight

# QC = QC_def([0],[0],[0])
# TE = TE_def([0],[0],[0])
# LE = LE_def([0],[0],[0])
# CP = CP_def([0],[0],[0],[0],[0],[0],[0],[0],[0])

plots = false
CL, CDi, CDp, CDc, CW, Cmac, cl_margin, gamma, CP = NVLM.VLM(wing, fs, ref, pdrag, mvr, plots)

d2 = ["TEST WITH MODIFIED SETUP - variable twist distribution, pass viscous drag, compressibility etc \n \n CL \n $CL \n  CDi \n $CDi \n  CDp \n $CDp \n  CDc \n $CDc \n  CW \n $CW \n  Cmac \n $Cmac \n  cl_margin \n $cl_margin \n  gamma \n $gamma"]
d = [d1,d2]
writedlm("check_output.csv", d)
