using VortexLattice
using StaticArrays
using DelimitedFiles
using FLOWMath
using Plots

function main()

    data_path="./VortexLattice_rotor_data"

    save_path = abspath(joinpath(@__DIR__, "..", "vortex_lattice_simulation"))
    # Empty out the output directory
    isdir(save_path) && rm(save_path, recursive=true, force=true)
    mkpath(save_path)

    constant_maneuver!(frames, system, wake, t) = nothing

    R = 63.0
    RPM = 9.1552
    rho             = 1.071778                  # (kg/m^3) air density
    mu              = 1.85508e-5                # (kg/ms) air dynamic viscosity
    speedofsound    = 342.35                    # (m/s) speed of sound
    magVinf         = -8.0
    Uinf(t) = SVector{3,Float64}(-1.0, 0.0, 0.0) * magVinf
    J               = magVinf/(RPM/60 * 2*R)

    ns = 20
    nc = 1

    grids, ratios, polars, frames = VortexLattice.generate_rotor("NREL5MW.csv", data_path;
                                                                turbine_flag=false,
                                                                clockwise=true,
                                                                ns,
                                                                nc,
                                                                spacing_s=Uniform(),
                                                                interpolate_airfoils=true);


    core_size = 1e-3
    nwakerows = 1
    system = System(grids; ratios, core_size, nw=fill(nwakerows, length(grids)));

    Sref = 1.0
    cref = 1.0
    bref = 1.0
    rref = [0.0, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)
    system.reference[] = ref

    # freestream parameters
    alpha = 0.0
    beta = 0.0
    Omega = [RPM * 2*pi/60; 0.0; 0.0]
    fs = Freestream(magVinf, alpha, beta, Omega)
    system.freestream[] = fs

    for isurf in eachindex(system.surfaces)
        VortexLattice.update_surface_panels!(system.surfaces[isurf], system.grids[isurf];
            ratios=system.ratios[isurf],
            fcore=(c, Δs) -> system.core_size)
    end

    frames = ReferenceFrame(system;
            origin = SVector{3}(0.0, 0.0, 0.0),
            v = SVector{3}(0.0, 0.0, 0.0),
            ω_axis = SVector{3}(1.0, 0.0, 0.0),
            ω = -RPM * 2 * pi / 60,
            R = SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
            name = "vehicle",
            child_index = Int[],
            dependent_index = collect(1:length(system.surfaces))
        )

    n_revs = 1
    ttot = n_revs / (RPM / 60)
    timestep_per_rev = 36
    t_range = range(start=0.0, stop=ttot, length=n_revs * timestep_per_rev + 1)
    overlap = 1.3
    p_per_step = 2
    nsteps_per_rev = length(t_range) / n_revs
    sigma = overlap * 2*pi*R / (nsteps_per_rev*p_per_step)

    monitor = VortexLattice.PanelForcesMonitor(length(t_range), system)
    monitor1 = VortexLattice.LiftingLineCoefficientsMonitor(length(t_range), system; normalized=false)

    fd_monitor = FluidDomainMonitor(
        range(-R, 3R, step=R/2),   # x: upstream to 3 diameters downstream
        range(-R, R,  step=R/2),   # y: ±1 radius
        range(-R, R,  step=R/2);   # z: ±1 radius
        vtk_interval = 5,
        name = "fluid_domain",
        path = joinpath(save_path, "fluid_domain"),
    )

    monitors = (monitor, monitor1, fd_monitor)
    monitors = (monitor, monitor1)
    Ωinf(_) = SVector{3,Float64}(0.0, 0.0, 0.0)

    fmm_wake = VortexLattice.fmm(;
        p = 20,
        ncrit = 3,
        autotune_p = true,
        autotune_ncrit = true,
        autotune_reg_error = true
    )

    fmm_vehicle = VortexLattice.fmm(;
        p = 20,
        ncrit = 3,
        autotune_p = true,
        autotune_ncrit = true,
        autotune_reg_error = true
    )

    wake = simulate!(system, frames, constant_maneuver!, Uinf, t_range, Ωinf;
                wake_type=PanelParticleWake,
                method_trailing=SigmaPPS(sigma, p_per_step),
                method_unsteady=SigmaPPS(sigma, p_per_step),
                eta=0.3,
                monitors,
                name = "NREL5MW",
                # path = nothing,
                write_restart = false,
                derivatives=false,
                polars,
                frames_index=fill(1, length(system.surfaces)),
                verbose=true,
                max_particles=50000,
                fmm_wake=fmm_wake,
                fmm_vehicle=fmm_vehicle,
            )

    # Normal force along the blade
    RHO = 1.0
    r_hub = 1.5
    dr = (R - r_hub) / ns
    x = r_hub .+ dr * (1:ns)

    # Panel forces: CF is (3, nc, ns, nt) — component 1 (axial/thrust), chordwise 1, all spans, last step
    F_panel = monitors[1].CF[1, 1, :, end] .* (0.5 * RHO * ref.V^2 * ref.S / dr)
    # Lifting line: normalized=false gives dimensional N/m
    F_lift = monitors[2].CF[1][1, :, end]

    p = plot(x, F_panel; xlabel="r (m)", ylabel="Normal force (N/m)", label="Panel forces")
    plot!(p, x, F_lift; label="Lifting line forces")
    display(p)

    # Coefficient of thrust vs time
    B = 3
    q = 0.5 * rho * magVinf^2
    A_disk = pi * R^2
    nt = length(t_range)
    CT_time = [abs(B * FLOWMath.trapz(x, monitors[2].CF[1][1, :, it])) / (q * A_disk) for it in 1:nt]
    println("Final CT: $(CT_time[end])")

    p2 = plot(t_range, CT_time; xlabel="t (s)", ylabel="CT", label="CT")
    display(p2)
end
main();
