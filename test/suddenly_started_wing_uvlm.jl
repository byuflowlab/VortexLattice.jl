using VortexLattice
using DelimitedFiles
using PythonPlot

# See Katz and Plotkin: Figure 13.37
# AR = ∞
# Vinf*Δt/c = 1/16
# α = 5°

# essentially infinite aspect ratio
AR = 1e2

# chord length
c = 1

# span length
b = AR*c

# planform area
S = b*c

# geometry
xle = [0.0, 0.0]
yle = [-b/2, b/2]
zle = [0.0, 0.0]
chord = [c, c]
theta = [0.0, 0.0]*pi/180
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section
ns = 1
nc = 4
spacing_s = Uniform()
spacing_c = Uniform()
mirror = false
symmetric = false

# reference parameters
cref = c
bref = b
Sref = S
rref = [0.0, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref, cref, bref, rref, Vinf)

# freestream parameters
alpha = 5.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha, beta, Omega)

# non-dimensional time (t*Vinf/c)
t = range(0.0, 7.0, step=1/8)

# time step
dt = [(t[i+1]-t[i]) for i = 1:length(t)-1]

# create vortex rings
grid, _ = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

_, _, surface = grid_to_surface_panels(grid) # Uniform spacing means ratios are not needed

# create vector containing all surfaces and grids
surfaces = [surface]
grids = [grid]

# run steady analysis
system = steady_analysis(grids, ref, fs; symmetric)

# extract steady forces
CFs, CMs = body_forces(system; frame=Wind())

# run transient analysis
system, surface_history, property_history, wake_history = unsteady_analysis(
    surfaces, ref, fs, dt; symmetric=symmetric)

# extract transient forces
CF, CM = body_forces_history(system, surface_history, property_history; frame=Wind())

# Computational Results
CL = getindex.(CF, 3)
CD = getindex.(CF, 1)
CLs = getindex(CFs, 3)
CDs = getindex(CFs, 1)

fig = figure("uvlm_ssw")
fig.clear()
fig.add_subplot(121, xlabel=L"t^*", ylabel=L"C_L")
fig.add_subplot(122, xlabel=L"t^*", ylabel=L"C_D")
ax = fig.get_axes()[0]
ax2 = fig.get_axes()[1]

ax.plot(t[2:end], CL./CLs, label="UVLM")

# analytic solution (Wagner's function)
Φ(t) = 1 - 0.165*exp(-0.045*t) - 0.335*exp(-0.3*t)
ax.plot(t, Φ.(2 .* t), ":", label="analytic")

# CDs
Ds_uvlm = [CF[i][1] for i in eachindex(CF)]
ax2.plot(t[2:end], CD ./ CDs, label="UVLM")

writedlm("uvlm_ssw.csv", hcat(t[2:end], CL, CD), ',')