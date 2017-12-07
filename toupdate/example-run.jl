#=
----------------------------------------------------------------------
 This version assumes (one) symmetric wing.

 Sample Inputs Below

 Author: S. Andrew Ning
----------------------------------------------------------------------
=#

include("VLM.jl")


# wing geometry  - all degrees in radians
span = [29, 65]
chord = [43, 26, 11]
twist = [0, 0, 0]*pi/180
tc = [0.13, 0.12, 0.11]
sweep = [25, 30]*pi/180
dihedral = [0, 0]
N = 50
wing = wingsection(span, chord, twist, tc, sweep, dihedral, N)

# println(wing.span)
# println(wing.chord)
# println(wing.chord[1:end-1]+wing.chord[2:end])

# freestream parameters - 2 methods either specify CL or angle of attack
mach = 0.78
alpha = 5*pi/180 # only used for alpha method
CL = 0.5 # only used for CL method
method = "alpha"
fs = fs_def(mach, alpha, CL, method)

# reference and other parameters (used for force/moment coefficients)
Sref = dot(wing.chord[1:end-1]+wing.chord[2:end], wing.span)
bref = 2*sum(wing.span.*cos(wing.dihedral))
cref = Sref/bref
CLmax = 1.4

ref = ref_def(Sref, cref, CLmax)

# parasite drag - 2 methods either PASS method or strip theory with a quadratic varition in cl
polar = [0.007 0 0.003] # only used in quad method
alt = 35000 # only used in pass method
xt = 0.05 # only used in pass method
method = "pass"
pdrag = pdrag_def(polar, alt, xt, method)

# structures
qN = 2.5 # ratio of mvr dynamic pressure to cruise dynamic pressure
n = 2.5 # limit load factor
kbar = 0 # coefficient for area-dependent weight
mvr = mvr_def(qN, n, kbar)

plots = true
CL, CDi, CDp, CW, Cmac, cl_margin = VLM(wing, fs, ref, pdrag, mvr, plots)
println(CL)
println(CDi)
println(CDp)
println(Cmac)
