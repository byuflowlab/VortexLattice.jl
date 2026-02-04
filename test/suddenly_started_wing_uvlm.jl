using VortexLattice
using PythonPlot

# Katz and Plotkin: Figures 13.34 and 13.35
# AR = [4, 8, 12, 20, ∞]
# Vinf*Δt/c = 1/16
# α = 5°

# AR = 20 # last aspect ratio is essentially infinite
AR = 6 # last aspect ratio is essentially infinite

# non-dimensional time (t*Vinf/c)
t = range(0.0, 10.0, step=1/16)
# t = range(0.0, 2.0, step=1/16)

# chord length
c = 1

# time step
dt = [t[i+1]-t[i] for i = 1:length(t)-1]

# span length
b = AR*c

# planform area
S = b*c

# geometry
xle = [0.0, 0.0]
yle = [-b/2, b/2]
zle = [0.0, 0.0]
chord = [c, c]
theta = [0.0, 0.0]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section
ns = 13
nc = 1

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

# create vortex rings
grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

# create vector containing grids
grids = [grid]
ratios = [ratio]

grid, ratio, surface = grid_to_surface_panels(grid; ratios=ratio)
surfaces = [surface]

# run analysis
system, surface_history, property_history, wake_history =
    unsteady_analysis(surfaces, ref, fs, dt; symmetric, wake_finite_core = false)

# extract forces at each time step
CF, CM = body_forces_history(system, surface_history,
    property_history; frame=Wind())

# save vtk files
write_vtk("uvlm_ssw", surface_history, property_history,
    wake_history, dt; symmetric=false)

fig = figure("uvlm_ssw")
fig.clear()
fig.add_subplot(121, xlabel=L"t^*", ylabel=L"C_L")
fig.add_subplot(122, xlabel=L"t^*", ylabel=L"C_D")
ax = fig.get_axes()[0]
ax2 = fig.get_axes()[1]

CLs_uvlm = [CF[i][3] for i in eachindex(CF)]
ax.plot(t[1:end-1]*Vinf/cref, CLs_uvlm, label="UVLM")
Ds_uvlm = [CF[i][1] for i in eachindex(CF)]
ax2.plot(t[1:end-1]*Vinf/cref, Ds_uvlm, label="UVLM")