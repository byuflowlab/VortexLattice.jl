using VortexLattice
using StaticArrays
using DelimitedFiles

data_path="./VortexLattice_rotor_data"

constant_maneuver!(frames, system, wake, t) = nothing

R = 0.12
RPM = 5400.0
J = 0.0001
rho             = 1.071778                  # (kg/m^3) air density
mu              = 1.85508e-5                # (kg/ms) air dynamic viscosity
speedofsound    = 342.35                    # (m/s) speed of sound
magVinf         = J*RPM/60*(2*R) * 0.0
Uinf(t) = SVector{3,Float64}(-1.0, 0.0, 0.0) * magVinf
Ωinf(t) = SVector{3,Float64}(0.0, 0.0, 0.0)

ns = 40
nc = 1

grids, ratios, polars, frames = VortexLattice.generate_rotor("DJI9443.csv", data_path;
                                                            turbine_flag=false,
                                                            clockwise=true,
                                                            ns,
                                                            nc,
                                                            spacing_s=Uniform(),
                                                            interpolate_airfoils=true,
                                                            RPM=RPM,);

nwakerows = 3
core_size = 1e-3
system = System(grids; ratios, core_size, nw=fill(nwakerows, length(grids)));

for isurf in eachindex(system.surfaces)
    VortexLattice.update_surface_panels!(system.surfaces[isurf], system.grids[isurf];
        ratios=system.ratios[isurf],
        fcore=(c, Δs) -> system.core_size)
end

Sref = 1.0
cref = 1.0
bref = 1.0
rref = [0.0, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref, cref, bref, rref, Vinf)
system.reference[] = ref

alpha = 0.0
beta = 0.0
Omega = [0; 0.0; 0.0]
fs = Freestream(magVinf, alpha, beta, Omega)
system.freestream[] = fs

save_path = abspath(joinpath(@__DIR__, "..", "rotor_hover_simulation"))
isdir(save_path) && rm(save_path, recursive=true, force=true)
mkpath(save_path)

n_revs = 5
ttot = n_revs / (RPM / 60)
timestep_per_rev = 36
t_range = range(start=0.0, stop=ttot, length=n_revs * timestep_per_rev + 1)
overlap = 1.3
p_per_step = 2
nsteps_per_rev = length(t_range) / n_revs
sigma = overlap * 2*pi*R / (nsteps_per_rev*p_per_step)

fmm_wake = VortexLattice.fmm(;
    ncrit = 3,
    p = 25,
    autotune_p = true,
    autotune_ncrit = true,
    autotune_reg_error = false,
)

fmm_vehicle = VortexLattice.fmm(;
    ncrit = 3,
    autotune_p = true,
    autotune_ncrit = true,
    autotune_reg_error = false,
)

monitors = (VortexLattice.ForcesMonitor(length(t_range)),)
wake = simulate!(system, frames, constant_maneuver!, Uinf, t_range, Ωinf;
            wake_type=PanelParticleWake,
            method_trailing=SigmaPPS(sigma, p_per_step),
            method_unsteady=NoShed(),
            eta=0.3,
            nwakerows=nwakerows,
            max_particles=1_000_000,
            fmm_wake=fmm_wake,
            fmm_vehicle=fmm_vehicle,
            monitors,
            name="rotorhover",
            path=save_path,
            derivatives=false,
            polars,
            frames_index=fill(1, length(system.surfaces)),
            verbose=true,
        )

# post-process
Ts = [monitors[1].CF[i][1] for i in 1:length(t_range)]
Ts .*= (0.5 * VortexLattice.RHO * ref.V^2 * ref.S)
CTs = Ts ./ (rho * (RPM/60)^2 * (2*R)^4)

# comparison
CT_exp = 0.072
CT_URANS = 0.071

di = timestep_per_rev * 1
CT_vpm = sum(CTs[end-di+1 : end]) / length(CTs[end-di+1 : end])
percent_error = abs((CT_vpm - CT_exp) / CT_exp) * 100
println("VPM CT: $CT_vpm\nExperiment CT: $CT_exp\nURANS CT: $CT_URANS\nPercent Error (VPM vs Experiment): $percent_error %")

# save csv with CT vs time
name = "rotor_hover_eta0.3_ns20_nc1_nt36_pps4_overlap1.3"
data = hcat(collect(t_range), CTs)
writedlm(name*".csv", data, ',')
