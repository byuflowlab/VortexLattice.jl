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
    nwakerows = 2
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
    Omega = [0; 0.0; 0.0]
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
    # t_range = range(start=0.0, stop=ttot/36, length=2)
    overlap = 1.3
    p_per_step = 2
    nsteps_per_rev = length(t_range) / n_revs
    sigma = overlap * 2*pi*R / (nsteps_per_rev*p_per_step)

    monitor = VortexLattice.PanelForcesMonitor(length(t_range), system)
    monitor1 = VortexLattice.LiftingLineCoefficientsMonitor(length(t_range), system; normalized=false)
    monitors = (monitor, monitor1)
    Ωinf(_) = SVector{3,Float64}(0.0, 0.0, 0.0)

    fmm_wake = VortexLattice.fmm(;
        p = 12,
        ncrit = 3,
        autotune_p = true,
        autotune_ncrit = true,
        autotune_reg_error = true
    )

    fmm_vehicle = VortexLattice.fmm(;
        p = 12,
        ncrit = 3,
        autotune_p = true,
        autotune_ncrit = true,
        autotune_reg_error = true
    )

    @profview_allocs wake = simulate!(system, frames, constant_maneuver!, Uinf, t_range, Ωinf;
                wake_type=PanelParticleWake,
                method_trailing=SigmaPPS(sigma, p_per_step),
                method_unsteady=SigmaPPS(sigma, p_per_step),
                eta=0.3,
                # monitors,
                name = "NREL5MW",
                path = save_path,
                derivatives=false,
                polars,
                frames_index=fill(1, length(system.surfaces)),
                verbose=true,
                max_particles=50000,
                fmm_wake=fmm_wake,
                fmm_vehicle=fmm_vehicle,
            )
    # RHO = 1
    # R = 63.0
    # r = 11.75
    # dr = (R - r) / ns
    # x = r .+ dr * (1:ns)
    # p = plot()
    # F_panel = monitors[1].CF[1,1,:,end-1] .* 0.5*RHO*Vinf^2 * ref.S ./ dr #Panel forces monitor
    # F_lift = monitors[2].CF[1][1,:,end] #Lifting line monitor
    # p = plot(x,F_panel, legend=true, xlabel="r (m)", ylabel="Force (N)",label="Panel forces")
    # p = plot(p, x,F_lift, legend=true, xlabel="r (m)", ylabel="Force (N/m)",label="Lifting line forces")
    # display(p)

    # Calculate coefficient of thrust using trapz integration
    # B = 3
    # T_blade = trapz(x, F)
    # T_total = B * T_blade
    # q = 0.5 * rho * magVinf^2
    # A_disk = pi * R^2
    # CT = abs(T_total) / (q * A_disk)
    # println("Turbine thrust coefficient (CT): $CT")
end
main();