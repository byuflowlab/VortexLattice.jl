# Allocation / time profiling harness for the PanelParticleWake simulate! path.
#
# ALWAYS warm up before profiling. The first call to run_case compiles the entire
# simulate! call tree, and profiling that call reports compiler allocations
# (Core.Compiler.typeinf and friends) rather than the allocations you care about.
# profile_case below does the warm-up for you.
#
# Usage (interactive REPL, with ProfileView loaded):
#
#   using ProfileView
#   include("test/profile_panel_particle_wake.jl")
#
#   profile_case(n_steps=200)                # warms up, then allocation profile
#   run_case(n_steps=4)                      # warm up by hand...
#   @profview run_case(n_steps=200)          # ...then CPU flame graph
#
# Or without ProfileView, for raw @time numbers:
#   include("test/profile_panel_particle_wake.jl")
#   run_case(n_steps=4)            # warm up
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

"""
    profile_case(; n_steps=36, warmup_steps=6, max_particles=20_000)

Allocation-profile `run_case`, warming up first so the profile reflects steady
state rather than compilation.

`warmup_steps` defaults to 6 rather than 1 because the particle-shedding and
buffer-overflow paths are not reached on the first step: with `nwakerows=1` the
wake buffer does not overflow into the particle field until step 2, and code that
first runs at step 3 would otherwise still be compiling inside the profiled call.
"""
function profile_case(; n_steps=72, warmup_steps=6, max_particles=20_000)
    run_case(; n_steps=warmup_steps, max_particles=max_particles)
    @profview_allocs run_case(; n_steps=n_steps, max_particles=max_particles) sample_rate=1E-3
    @time run_case(; n_steps=n_steps, max_particles=max_particles)
    return nothing
end

profile_case()