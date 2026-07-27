# Allocation / time profiling harness for the PanelParticleWake simulate! path.
#
# Usage (interactive REPL, with ProfileView loaded):
#
#   using ProfileView
#   include("test/profile_panel_particle_wake.jl")
#
#   run_case()                    # warm up / precompile
#   @profview run_case(n_steps=200)        # CPU flame graph
#   @profview_allocs run_case(n_steps=200) sample_rate=1.0   # allocation profile
#
# Or without ProfileView, for raw @time numbers:
#   include("test/profile_panel_particle_wake.jl")
#   run_case()                     # warm up
#   @time run_case(n_steps=200)    # timed, steady-state allocation count

using StaticArrays
using LinearAlgebra
using VortexLattice
using FLOWVPM

function build_case(; n_steps=200, alpha=5.0*pi/180)
    c  = 1.0
    b  = 8.0
    xle   = [0.0, 0.0]
    yle   = [-b/2, b/2]
    zle   = [0.0, 0.0]
    chord = [c, c]
    theta = [0.0, 0.0]
    phi   = [0.0, 0.0]
    fc    = fill((xc) -> 0, 2)
    ns    = 13
    nc    = 1

    grid, ratios = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=false, fc=fc, spacing_s=Uniform(), spacing_c=Uniform())

    nwakerows = 1
    system = System([grid]; nw=[nwakerows], ratios=[ratios])

    Vinf_mag = 10.0
    ref = Reference(b*c, c, b, [0.0, 0.0, 0.0], Vinf_mag)
    system.reference[] = ref
    system.freestream[] = Freestream(Vinf_mag, alpha, 0.0, [0.0, 0.0, 0.0])

    for isurf in eachindex(system.surfaces)
        VortexLattice.update_surface_panels!(system.surfaces[isurf], system.grids[isurf];
            ratios=system.ratios[isurf],
            fcore=(c, Δs) -> system.core_size)
    end

    frames = ReferenceFrame(system; origin=SVector{3,Float64}(0.0, 0.0, 0.0))

    Vinf_vec = SVector{3,Float64}(Vinf_mag*cos(alpha), 0.0, Vinf_mag*sin(alpha))
    Vinf_func(t) = Vinf_vec
    Ωinf_func(t) = SVector{3,Float64}(0.0, 0.0, 0.0)
    maneuver!(frames, system, wake, t) = false

    dt = 0.05
    t_range = collect(0.0:dt:(n_steps*dt))

    return system, frames, maneuver!, Vinf_func, Ωinf_func, t_range, nwakerows
end

function run_case(; n_steps=36, max_particles=20_000, verbose=false)
    system, frames, maneuver!, Vinf_func, Ωinf_func, t_range, nwakerows = build_case(; n_steps=n_steps)

    wake = simulate!(system, frames, maneuver!, Vinf_func, t_range, Ωinf_func;
        wake_type=PanelParticleWake,
        nwakerows=nwakerows, max_particles=max_particles,
        method_trailing=OverlapPPS(1.3, 2),
        method_unsteady=OverlapPPS(1.3, 2),
        name=nothing, path=nothing, verbose=verbose)

    return system, wake
end

run_case()