using Test
using StaticArrays
using LinearAlgebra
using VortexLattice
using FLOWVPM
const VL = VortexLattice

@testset "Part 3: write_vtk(System) and write_vtk(PanelParticleWake)" begin

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

    grid, ratios = VL.wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=false, fc=fc, spacing_s=VL.Uniform(), spacing_c=VL.Uniform())

    nwakerows = 3
    system = VL.System([grid]; nw=[nwakerows], ratios=[ratios])

    # reference + freestream so steady_analysis! populates properties / Γ
    Vinf_mag = 10.0
    Sref = b * c
    ref  = VL.Reference(Sref, c, b, [0.0, 0.0, 0.0], Vinf_mag)
    fs   = VL.Freestream(Vinf_mag, 0.0, 0.0, [0.0, 0.0, 0.0])
    system.reference[] = ref
    # pass pre-sized wake buffers so steady_analysis! doesn't clobber them
    init_wakes = [Matrix{VL.WakePanel{Float64}}(undef, nwakerows, ns)]
    VL.steady_analysis!(system, ref, fs; symmetric=false, wakes=init_wakes,
        nwake=[0], derivatives=false)

    # --- seed wake_shedding_locations from the TE of the surface ---
    surface = system.surfaces[1]
    wsl = system.wake_shedding_locations[1]
    for j in 1:ns
        wsl[j] = VL.bottom_left(surface[end, j])
    end
    wsl[ns+1] = VL.bottom_right(surface[end, ns])

    # --- seed wake_velocities with freestream so panels translate on propagate! ---
    Vinf_vec = SVector{3,Float64}(Vinf_mag, 0.0, 0.0)
    for V in system.V
        fill!(V, Vinf_vec)
    end

    dt = 0.05
    n_steps = nwakerows + 3   # a few overflow steps so we get particles

    wake = VL.PanelParticleWake(system;
        nwakerows=nwakerows, max_particles=1000,
        method_trailing=VL.OverlapPPS(1.3, 2),
        method_unsteady=VL.OverlapPPS(1.3, 2),
    )

    # --- dedicated output dir ---
    outdir = mktempdir(; prefix="vl_part3_viz_")
    body_name = joinpath(outdir, "wing_bodies")
    wake_name = joinpath(outdir, "wing_wake")

    # --- initial write (step 0): empty wake + zero particles ---
    VL.write_vtk(body_name, system, 0, 0.0; overwrite=true)
    VL.write_vtk(wake_name, wake,   0, 0.0; overwrite=true)

    @test isfile(body_name * ".pvd")
    @test isfile(wake_name * ".pvd")
    @test isfile(wake_name * "_particles.pvd")
    @test isfile(joinpath(outdir, "wing_bodies", "wing_bodies_0.vtm"))
    # wake VTM exists but has no surface blocks yet (nwake==0)
    @test isfile(joinpath(outdir, "wing_wake", "wing_wake_0.vtm"))
    # particle VTP exists (empty)
    @test isfile(joinpath(outdir, "wing_wake_particles", "wing_wake_particles_0.vtp"))

    # --- time-step loop: shed, propagate, write ---
    # One Γ per panel (flat vector over all surfaces)
    ntot = sum(length.(system.surfaces))
    Gamma = ones(Float64, ntot)

    for i in 1:n_steps
        t = i * dt

        VL.update_TE!(wake, system)
        VL.shed_wake!(wake, system, dt, Gamma)
        VL.propagate!(wake, dt; relax=false)

        VL.write_vtk(body_name, system, i, t; overwrite=false)
        VL.write_vtk(wake_name, wake,   i, t; overwrite=false)

        # per-step files exist
        @test isfile(joinpath(outdir, "wing_bodies", "wing_bodies_$i.vtm"))
        @test isfile(joinpath(outdir, "wing_wake",   "wing_wake_$i.vtm"))
        @test isfile(joinpath(outdir, "wing_wake_particles",
                              "wing_wake_particles_$i.vtp"))
    end

    # --- after n_steps > nwakerows, particles must have been shed ---
    @test wake.overflowed[] == true
    @test FLOWVPM.get_np(wake.pfield) > 0
    @test wake.nwake[1] == nwakerows

    # --- PVD index grew with every write (n_steps + 1 entries) ---
    body_pvd = read(body_name * ".pvd", String)
    wake_pvd = read(wake_name * ".pvd", String)
    parts_pvd = read(wake_name * "_particles.pvd", String)
    # count DataSet entries — one per write_vtk call
    @test count("<DataSet", body_pvd) == n_steps + 1
    @test count("<DataSet", wake_pvd) == n_steps + 1
    @test count("<DataSet", parts_pvd) == n_steps + 1

    @info "Part 3 VTK output written to $outdir"
end
