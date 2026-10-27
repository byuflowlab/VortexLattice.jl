#=
Weber & Brebner's 45° swept wing (RAE 101 12 % section, AR 4.9) at 4.2°: sectional lift along
the span against the experiment, from the steady inviscid analysis and from an unsteady run
on FLOWUnsteadyCore with the RAE 101 polar's viscous correction. Not part of the unit suite.
Run from the package root:
    julia --project=<env with VortexLattice, CCBlade, StaticArrays, Plots> test/weber_wing.jl
=#
using VortexLattice
using CCBlade
using StaticArrays
using Plots

inch2meter(x) = x * 0.0254
ft2meter(x) = x * 0.3048

ns = 20                 # spanwise panels per half
nc = 1
aoa = 4.2 * pi / 180
Vinf = ft2meter(163.0)

halfspan = inch2meter(98.0) / 2
root_chord = inch2meter(20.0)
sweep = deg2rad(45.0)
xle = [0.0, tan(sweep) * halfspan]
yle = [0.0, halfspan]
chord = [root_chord, root_chord]

ref = Reference(2halfspan * root_chord, root_chord, 2halfspan, [0.0, 0.0, 0.0], Vinf)
mkgrid() = wing_to_grid(xle, yle, zeros(2), chord, zeros(2), zeros(2), ns, nc;
    mirror = true, spacing_s = Uniform(), spacing_c = Uniform())

# sectional lift coefficient along the right half (strip centres at 2y/b)
function sectional_cl(system)
    cf, _ = lifting_line_coefficients(system; normalized = true)
    return [cf[1][3, j] for j in ns+1:2ns]
end
eta = [(j - 0.5) / ns for j in 1:ns]

# --- steady, inviscid ---
grid, ratios = mkgrid()
steady = VortexLatticeSystem([grid]; ratios = [ratios])
steady_analysis!(steady, ref, flight_condition(steady, ref, Vinf; alpha = aoa)...)
cl_steady = sectional_cl(steady)

# --- unsteady on FLOWUnsteadyCore, with the RAE 101 polar ---
polar = CCBlade.AlphaAF(joinpath(@__DIR__, "..", "VortexLattice_rotor_data", "airfoils", "rae101_12.dat"); radians = false)
grid, ratios = mkgrid()
system = VortexLatticeSystem([grid]; ratios = [ratios])
system.reference[] = ref
bodies, states = surface_bodies(system; polars = [fill(polar, 2ns)])
frames = ReferenceFrame(system; origin = SVector(0.0, 0.0, 0.0))
freestream = Freestream(SVector(Vinf * cos(aoa), 0.0, Vinf * sin(aoa)))
nsteps = 200
dt = root_chord / Vinf / 4                  # a quarter chord per step: 50 chords
wake = PanelParticleWake(Float64, bodies; nwakerows = 2, max_particles = 50_000)
monitor = ForcesMonitor(system, nsteps; frame = Wind())
simulate!(states, bodies, frames, wake, freestream; nsteps, dt, monitor)
cl_unsteady = sectional_cl(system)

# --- the experiment (Weber & Brebner, alpha = 4.2°): sectional cl at 2y/b ---
true_eta = [0, 0.041, 0.082, 0.163, 0.245, 0.367, 0.510, 0.653, 0.898, 0.949]
true_cl = [0.235, 0.241, 0.248, 0.253, 0.251, 0.251, 0.251, 0.246, 0.192, 0.171]

println("wing CL: steady inviscid ", round(body_forces(steady; frame = Wind())[1][3], sigdigits = 4),
        ", unsteady with polar ", round(monitor.CF[end][3], sigdigits = 4))
println("2y/b    experiment  steady  unsteady+polar")
for (e, c) in zip(true_eta, true_cl)
    k = clamp(round(Int, e * ns + 0.5), 1, ns)
    println(rpad(e, 8), rpad(c, 12), rpad(round(cl_steady[k], digits = 3), 8), round(cl_unsteady[k], digits = 3))
end

plot(xlabel = "2y/b", ylabel = "sectional cl", title = "Weber & Brebner 45° swept wing, α = 4.2°", legend = :bottomleft)
scatter!(true_eta, true_cl, label = "Experiment", markershape = :diamond)
plot!(eta, cl_steady, label = "VortexLattice, steady inviscid")
plot!(eta, cl_unsteady, label = "VortexLattice on FLOWUnsteadyCore, RAE 101 polar")
savefig("weber_wing_cl.png")
