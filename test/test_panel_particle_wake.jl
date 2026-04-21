using StaticArrays
using LinearAlgebra
using Test
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

function build_restart_case()
    c = 1.0
    b = 8.0
    xle = [0.0, 0.0]
    yle = [-b/2, b/2]
    zle = [0.0, 0.0]
    chord = [c, c]
    theta = [0.0, 0.0]
    phi = [0.0, 0.0]
    fc = fill((xc) -> 0, 2)

    grid, ratios = wing_to_grid(xle, yle, zle, chord, theta, phi, 2, 2;
        mirror=false, fc=fc, spacing_s=Uniform(), spacing_c=Uniform())

    system = System([grid]; nw=[2], ratios=[ratios])
    ref = Reference(b*c, c, b, [0.0, 0.0, 0.0], 10.0)
    system.reference[] = ref
    system.freestream[] = Freestream(10.0, 0.0, 0.0, [0.0, 0.0, 0.0])

    for isurf in eachindex(system.surfaces)
        VortexLattice.update_surface_panels!(system.surfaces[isurf], system.grids[isurf];
            ratios=system.ratios[isurf], fcore=(c, Δs) -> system.core_size)
    end

    frames = ReferenceFrame(system; origin=SVector{3,Float64}(0.0, 0.0, 0.0))
    maneuver!(frames, system, wake, t) = nothing
    Vinf_func(t) = SVector{3,Float64}(10.0, 0.0, 0.0)
    Ωinf_func(t) = SVector{3,Float64}(0.0, 0.0, 0.0)
    t_range = collect(0.0:0.05:0.10)

    return system, frames, maneuver!, Vinf_func, Ωinf_func, t_range
end

@testset "PanelParticleWake restart checkpoint" begin
    baseline_dir = mktempdir()
    restart_dir = mktempdir()

    system_a, frames_a, maneuver_a, Vinf_a, Ωinf_a, t_range = build_restart_case()
    baseline_wake = simulate!(system_a, frames_a, maneuver_a, Vinf_a, t_range, Ωinf_a;
        wake_type=PanelParticleWake,
        nwakerows=2,
        max_particles=500,
        eta=0.3,
        method_trailing=OverlapPPS(1.3, 2),
        method_unsteady=OverlapPPS(1.3, 2),
        name="full_run",
        path=baseline_dir,
        derivatives=false,
        verbose=false,
    )

    system_b, frames_b, maneuver_b, Vinf_b, Ωinf_b, _ = build_restart_case()
    partial_wake = simulate!(system_b, frames_b, maneuver_b, Vinf_b, t_range[1:2], Ωinf_b;
        wake_type=PanelParticleWake,
        nwakerows=2,
        max_particles=500,
        eta=0.3,
        method_trailing=OverlapPPS(1.3, 2),
        method_unsteady=OverlapPPS(1.3, 2),
        name="restart_run",
        path=restart_dir,
        derivatives=false,
        verbose=false,
    )

    system_c, frames_c, maneuver_c, Vinf_c, Ωinf_c, _ = build_restart_case()
    resumed_wake = simulate!(system_c, frames_c, maneuver_c, Vinf_c, t_range, Ωinf_c;
        wake_type=PanelParticleWake,
        nwakerows=2,
        max_particles=500,
        eta=0.3,
        method_trailing=OverlapPPS(1.3, 2),
        method_unsteady=OverlapPPS(1.3, 2),
        name="restart_run",
        path=restart_dir,
        restart_from=joinpath(restart_dir, "restart_run"),
        restart_idx=1,
        derivatives=false,
        verbose=false,
    )

    @test baseline_wake.nwake == resumed_wake.nwake
    @test baseline_wake.overflowed[] == resumed_wake.overflowed[]
    @test baseline_wake.pfield.np == resumed_wake.pfield.np
    @test baseline_wake.prev_bottom_gamma == resumed_wake.prev_bottom_gamma
    @test isapprox(system_a.Γ, system_c.Γ; atol=0, rtol=0)
    @test isapprox(system_a.freestream[].Vinf, system_c.freestream[].Vinf; atol=0, rtol=0)
end

function build_rotor_restart_case()
    data_path = joinpath(@__DIR__, "..", "VortexLattice_rotor_data")
    ns = 12
    nc = 1
    RPM = 4200.0

    grids, ratios, _, _ = VortexLattice.generate_rotor("DJI9443.csv", data_path;
        turbine_flag=false,
        clockwise=true,
        ns,
        nc,
        spacing_s=Uniform(),
        interpolate_airfoils=true,
        RPM=RPM)

    nwakerows = 2
    core_size = 1e-3
    system = System(grids; ratios, core_size, nw=fill(nwakerows, length(grids)))

    ref = Reference(1.0, 1.0, 1.0, [0.0, 0.0, 0.0], 1.0)
    system.reference[] = ref
    system.freestream[] = Freestream(0.0, 0.0, 0.0, [0.0, 0.0, 0.0])

    for isurf in eachindex(system.surfaces)
        VortexLattice.update_surface_panels!(system.surfaces[isurf], system.grids[isurf];
            ratios=system.ratios[isurf],
            fcore=(c, Δs) -> system.core_size)
    end

    frames = ReferenceFrame(system;
        origin=SVector{3}(0.0, 0.0, 0.0),
        v=SVector{3}(0.0, 0.0, 0.0),
        ω_axis=SVector{3}(1.0, 0.0, 0.0),
        ω=-RPM * 2 * pi / 60,
        R=SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0,
                                 0.0, 1.0, 0.0,
                                 0.0, 0.0, 1.0),
        name="vehicle",
        child_index=Int[],
        dependent_index=collect(1:length(system.surfaces)))

    maneuver!(frames, system, wake, t) = nothing
    Uinf(t) = SVector{3,Float64}(0.0, 0.0, 0.0)
    Ωinf(t) = SVector{3,Float64}(0.0, 0.0, 0.0)

    n_revs = 0.4
    timestep_per_rev = 24
    ttot = n_revs / (RPM / 60)
    t_range = collect(range(start=0.0, stop=ttot, length=Int(round(n_revs * timestep_per_rev)) + 1))

    return system, frames, maneuver!, Uinf, Ωinf, t_range
end

@testset "PanelParticleWake restart visual comparison outputs" begin
    visual_root = joinpath(@__DIR__, "restart_visual_compare_output")
    full_dir = joinpath(visual_root, "full_run")
    restart_dir = joinpath(visual_root, "restart_run")
    isdir(visual_root) && rm(visual_root; recursive=true, force=true)
    mkpath(full_dir)
    mkpath(restart_dir)

    system_full, frames_full, maneuver_full, Vinf_full, Ωinf_full, t_range_vis = build_rotor_restart_case()
    restart_idx = max(2, Int(floor((length(t_range_vis)-1) / 2)))

    # Full continuous rotor run with rotating-frame kinematics.
    wake_full = simulate!(system_full, frames_full, maneuver_full, Vinf_full, t_range_vis, Ωinf_full;
        wake_type=PanelParticleWake,
        nwakerows=2,
        max_particles=5000,
        eta=0.3,
        method_trailing=OverlapPPS(1.3, 2),
        method_unsteady=OverlapPPS(1.3, 2),
        name="full_run",
        path=full_dir,
        derivatives=false,
        verbose=false,
    )

    # Run to restart point only.
    system_part, frames_part, maneuver_part, Vinf_part, Ωinf_part, _ = build_rotor_restart_case()
    simulate!(system_part, frames_part, maneuver_part, Vinf_part,
        t_range_vis[1:restart_idx+1], Ωinf_part;
        wake_type=PanelParticleWake,
        nwakerows=2,
        max_particles=5000,
        eta=0.3,
        method_trailing=OverlapPPS(1.3, 2),
        method_unsteady=OverlapPPS(1.3, 2),
        name="restart_run",
        path=restart_dir,
        derivatives=false,
        verbose=false,
    )

    # Resume from checkpoint and finish the same timeline.
    system_resume, frames_resume, maneuver_resume, Vinf_resume, Ωinf_resume, _ = build_rotor_restart_case()
    wake_resume = simulate!(system_resume, frames_resume, maneuver_resume, Vinf_resume, t_range_vis, Ωinf_resume;
        wake_type=PanelParticleWake,
        nwakerows=2,
        max_particles=5000,
        eta=0.3,
        method_trailing=OverlapPPS(1.3, 2),
        method_unsteady=OverlapPPS(1.3, 2),
        name="restart_run",
        path=restart_dir,
        restart_from=joinpath(restart_dir, "restart_run"),
        restart_idx=restart_idx,
        derivatives=false,
        verbose=false,
    )

    @test wake_full.nwake == wake_resume.nwake
    @test wake_full.overflowed[] == wake_resume.overflowed[]
    @test wake_full.pfield.np == wake_resume.pfield.np
    @test isapprox(system_full.Γ, system_resume.Γ; atol=1e-3, rtol=1e-3)

    full_bodies = joinpath(full_dir, "full_run_bodies.pvd")
    restart_bodies = joinpath(restart_dir, "restart_run_bodies.pvd")
    full_wake = joinpath(full_dir, "full_run_wake.pvd")
    restart_wake = joinpath(restart_dir, "restart_run_wake.pvd")
    @test isfile(full_bodies)
    @test isfile(restart_bodies)
    @test isfile(full_wake)
    @test isfile(restart_wake)

    println("visual full-run output: ", full_dir)
    println("visual restarted-run output: ", restart_dir)
    println("open full_run_bodies.pvd/full_run_wake.pvd and restart_run_bodies.pvd/restart_run_wake.pvd in ParaView")
end
