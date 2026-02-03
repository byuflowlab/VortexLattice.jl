using VortexLattice
using StaticArrays
using PythonPlot

constant_maneuver!(frames, system, wake, t) = nothing

# generate wing
# AR = 20.0
AR = 6.0
c = 1.0
b = AR * c
xle = [0.0, 0.0]
yle = [-b/2, b/2]
zle = [0.0, 0.0]
chord = [c, c]
theta = [0.0, 0.0]
phi = [0.0, 0.0]
ns = 1
# ns = 1
# nc = 4
nc = 1
fc = fill((xc) -> 0, length(yle)) # camberline function for each section
spacing_s = Uniform()
spacing_c = Uniform()
mirror = false
grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

grids = [grid]
ratios = [ratio]

system_ssw = System(grids; ratios)

Sref = c*b
cref = c
bref = b
rref = [0.0, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref, cref, bref, rref, Vinf)
system_ssw.reference[] = ref

# freestream parameters
alpha = 5.0 * pi/180
beta = 0.0
Omega = [0; 0.0; 0.0]
fs = Freestream(Vinf, alpha, beta, Omega)
system_ssw.freestream[] = fs

println("steady analysis:")
steady_analysis!(system_ssw, system_ssw.reference[], system_ssw.freestream[]; symmetric=false, trailing_vortices=true)

steady_system_ssw = deepcopy(system_ssw)

CF_steady, _ = body_forces(steady_system_ssw; frame=Wind())
CL_steady = CF_steady[3]
CD_steady = CF_steady[1]

write_vtk("ssw_steady", steady_system_ssw; trailing_vortices=true)
frames = ReferenceFrame(system_ssw;
        origin = SVector{3}(0.0, 0.0, 0.0),
        v = SVector{3}(0.0, 0.0, 0.0),
        ω_axis = SVector{3}(1.0, 0.0, 0.0),
        ω = 0.0,
        R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name = "vehicle",
        child_index = Int[],
        dependent_index = collect(1:length(system_ssw.surfaces))
    )

Vinf_func(t) = SVector{3,Float64}(cos(alpha), 0.0, sin(alpha)) * Vinf
dt = c / Vinf / 16
max_Vinf_t_c = 10.0
tmax = max_Vinf_t_c * cref / Vinf
nsteps = Int(ceil(tmax / dt))
t_range = range(start=0.0, stop=nsteps*dt, length=nsteps+1)
monitors = (VortexLattice.ForcesMonitor(length(t_range); frame=Wind()),)
benchmark = @elapsed wake = simulate!(system_ssw, frames, constant_maneuver!, Vinf_func, t_range; 
            monitors,
            # particle_trailing_methods=fill(VortexLattice.NoShed(), length(system_ssw.surfaces)),
            particle_trailing_methods=fill(VortexLattice.OverlapPPS(1.3, 1), length(system_ssw.surfaces)),
            # particle_unsteady_methods=fill(VortexLattice.SigmaPPS(5.5, 1), length(system_ssw.surfaces)),
            # particle_unsteady_methods=fill(VortexLattice.OverlapPPS(1.3, 80), length(system_ssw.surfaces)),
            particle_unsteady_methods=fill(VortexLattice.NoShed(), length(system_ssw.surfaces)),
            eta = 0.25,
            derivatives = false,
            vtk_args=(trailing_vortices=false,),
            fmm_wake_args=(leaf_size_source=1000,),
            # nonlinear_analysis=true,
            # nonlinear_args=(polar_correction=false,),
            # calculate_influence_matrix=true,
            # path=nothing,
            name="suddenly_started_wing",
            # wake_args=(relaxation=VortexLattice.FLOWVPM.relaxation_none,)
        )

CLs = [monitors[1].CF[i][3] for i in eachindex(monitors[1].CF)]
# CLs = Ls ./ (0.5 * Vinf^2 * Sref)
# Ds = [monitors[1].CF[i][1] for i in eachindex(monitors[1].CF)]
CDs = [monitors[1].CF[i][1] for i in eachindex(monitors[1].CF)]
# CDs = Ds ./ (0.5 * Vinf^2 * Sref)
# cls = Ls ./ (0.5 * Vinf^2 * c * b)
tstar = Vinf * collect(t_range) / cref

fig = figure("vpm")
fig.clear()
fig.add_subplot(121, xlabel=L"t^*", ylabel=L"C_L")
fig.add_subplot(122, xlabel=L"t^*", ylabel=L"C_D")
ax = fig.get_axes()[0]
ax.plot(tstar, CLs, label="VPM")
ax.plot(tstar, fill(CL_steady, length(tstar)), "--", label="steady VLM")
ax.set_ylim(0.0, 1.0)
ax.legend()
ax2 = fig.get_axes()[1]
ax2.plot(tstar, CDs, label="VPM")
ax2.plot(tstar, fill(CD_steady, length(tstar)), "--", label="steady VLM")
# ax2.set_ylim(0.0, 0.008)
ax2.legend()
tight_layout()

# check rotation/frame