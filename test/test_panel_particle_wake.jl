using StaticArrays
using LinearAlgebra
using VortexLattice
using FLOWVPM

# --- minimal rectangular wing ---
c  = 1.0
b  = 8.0
xle   = [0.0, 0.0]
yle   = [-b/2, b/2]
zle   = [0.0, 0.0]
chord = [c, c]
theta = [0.0, 0.0]
phi   = [0.0, 0.0]
fc    = fill((xc) -> 0, 2)
ns    = 4
nc    = 2

grid, ratios = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=false, fc=fc, spacing_s=Uniform(), spacing_c=Uniform())

nwakerows = 3
system = System([grid]; nw=[nwakerows], ratios=[ratios])

# reference / initial freestream
Vinf_mag = 10.0
ref = Reference(b*c, c, b, [0.0, 0.0, 0.0], Vinf_mag)
system.reference[] = ref
system.freestream[] = Freestream(Vinf_mag, 0.0, 0.0, [0.0, 0.0, 0.0])

# initialize surface panels from the grid
for isurf in eachindex(system.surfaces)
    VortexLattice.update_surface_panels!(system.surfaces[isurf], system.grids[isurf];
        ratios=system.ratios[isurf],
        fcore=(c, Δs) -> system.core_size)
end

# frames (single top-level vehicle frame)
frames = ReferenceFrame(system; origin=SVector{3,Float64}(0.0, 0.0, 0.0))

# scalar callbacks for simulate!
Vinf_vec = SVector{3,Float64}(Vinf_mag, 0.0, 0.0)
Vinf_func(t) = Vinf_vec
Ωinf_func(t) = SVector{3,Float64}(0.0, 0.0, 0.0)
maneuver!(frames, system, wake, t) = false

# time range
dt = 0.05
n_steps = 12
t_range = collect(0.0:dt:(n_steps*dt))

outdir = joinpath(@__DIR__, "panel_particle_wake_output")
isdir(outdir) && rm(outdir; recursive=true, force=true)

wake = simulate!(system, frames, maneuver!, Vinf_func, t_range, Ωinf_func;
    wake_type=PanelParticleWake,
    nwakerows=nwakerows, max_particles=1000,
    method_trailing=OverlapPPS(1.3, 2),
    method_unsteady=OverlapPPS(1.3, 2),
    name="wing", path=outdir)

println("nwake             = ", wake.nwake)
println("overflowed        = ", wake.overflowed[])
println("num particles     = ", FLOWVPM.get_np(wake.pfield))
println("prev_bottom_gamma = ", wake.prev_bottom_gamma)
println("output written to: ", outdir)
