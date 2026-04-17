using VortexLattice
using CCBlade
using StaticArrays

save_path = "vortex_lattice_simulation"
# Empty out the data_path directory
isdir(save_path) && rm(save_path, recursive=true, force=true)
mkdir(save_path)

const rho = 1.225

inch2meter(x) = x * 0.0254
ft2meter(x) = x * 0.3048

ns = 20
nc = 1

constant_maneuver!(frames, system, wake, t) = nothing

aoa=4.2 
Vinf=ft2meter(163.0)

polars = CCBlade.AlphaAF("VortexLattice_rotor_data/airfoils/rae101_12.dat"; radians=false)

polars = VortexLattice.Polar(polars.alpha, polars.cl, polars.cd)
polars = fill(polars, ns) # wrap in array for compatibility with simulate!
polars = [polars]

halfspan = inch2meter(98.0)/2
root_chord = inch2meter(20.0)

# 45° sweep (swept-back tip)
sweep_deg = 45.0
sweep = deg2rad(sweep_deg)
xle = [0.0, tan(sweep) * halfspan]
yle = [0.0, halfspan]
zle = zeros(2)
chord = [root_chord, root_chord]  # adjust tip chord ratio
theta = zeros(2)
phi = zeros(2)

# reference parameters
Sref = halfspan * 2 * root_chord
cref = root_chord
bref = halfspan * 2
rref = [0.0, 0.0, 0.0]
ref = Reference(Sref, cref, bref, rref, Vinf)

spacing_s = Uniform() # spanwise discretization scheme
spacing_c = Uniform() # chordwise discretization scheme
grid, ratios = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc; spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)

grids = [grid]
ratios = [ratios]

system = System(grids; ratios)
empty_system = deepcopy(system)
system.reference[] = ref
fs = Freestream(Vinf, deg2rad(aoa), 0.0, [0.0; 0.0; 0.0])
steady_analysis!(system, system.reference[], fs; symmetric=false)

cf_steady, cm_steady = lifting_line_coefficients(system; normalized=true)

frames = ReferenceFrame(system;
    origin = SVector{3}(0.0, 0.0, -10.0),
    v = SVector{3}(0.0, 0.0, 0.0),
    ω_axis = SVector{3}(0.0, 1.0, 0.0),
    ω = 0.0 * 2 * pi,
    R = SMatrix{3,3}(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0),
    name = "vehicle",
    child_index = Int[],
    dependent_index = collect(1:length(system.surfaces))
)

alpha = aoa * (pi/180)
Uinf(t) = SVector{3,Float64}(Vinf * cos(alpha), 0.0, Vinf * sin(alpha))
t_range = range(start=0.0, stop=0.125, length=51)
monitors = (VortexLattice.ForcesMonitor(length(t_range)),)

wake = simulate!(system, frames, constant_maneuver!, Uinf, t_range; 
            monitors,
            # particle_trailing_methods=fill(VortexLattice.NoShed(), length(system.surfaces)),
            particle_trailing_methods=fill(VortexLattice.OverlapPPS(1.3,1), length(system.surfaces)),
            # particle_unsteady_methods=fill(VortexLattice.OverlapPPS(1.3,5), length(system.surfaces)),
            particle_unsteady_methods=fill(VortexLattice.NoShed(), length(system.surfaces)),
            eta = 1.0,
            derivatives = false,
            vtk_args=(trailing_vortices=false,),
            polars,
            nwakerows=1
            )

# alpha = 4.2
true_x = vec([0 0.041 0.082 0.163 0.245 0.367 0.510 0.653 0.898 0.949])
true_cl = vec([0.235 0.241 0.248 0.253 0.251 0.251 0.251 0.246 0.192 0.171])
true_cd = vec([0.059 0.025 0.016 0.009 0.007 0.006 0.006 0.004 -0.002 -0.007])
p = scatter(true_x, true_cl, label="Experimental Data", markershape=:diamond, ylabel="Cl", xlabel="x/c", title="Lift Coefficient Distribution at α=4.2°")
p_cd = scatter(true_x, true_cd, label="Experimental Data", markershape=:diamond, ylabel="Cd", xlabel="x/c", title="Drag Coefficient Distribution at α=4.2°")

cf, cm = lifting_line_coefficients(system; normalized=true);