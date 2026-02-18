using VortexLattice
using StaticArrays
# using PythonPlot
using DelimitedFiles

data_path="./VortexLattice_rotor_data"

constant_maneuver!(frames, system, wake, t) = nothing

R = 63.0
RPM = 9.1552
rho             = 1.071778                  # (kg/m^3) air density
mu              = 1.85508e-5                # (kg/ms) air dynamic viscosity
speedofsound    = 342.35                    # (m/s) speed of sound
magVinf         = 8.0
Uinf(t) = SVector{3,Float64}(-1.0, 0.0, 0.0) * magVinf
J               = magVinf/(RPM/60 * 2*R)

ns = 20
nc = 1

grids, ratios = VortexLattice.generate_rotor("NREL5MW.csv", data_path; 
                                                            turbine_flag=false, 
                                                            clockwise=true, 
                                                            ns, 
                                                            nc, 
                                                            spacing_s=Uniform(), 
                                                            interpolate_airfoils=true);


core_size = 1e-3
system = System(grids; ratios, core_size);
# system = System(grids; ratios, sections);

Sref = 2.0
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

steady_analysis!(system, system.reference[], system.freestream[]; symmetric=false);

write_vtk("rotor_hover_initial", system; write_wakes=false, trailing_edge_list=fill(false, length(system.surfaces)));

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

# generate correction functions
#### include("polar_correction.jl") # generates cl_correction and cd_correction functions based on XFOIL data for the airfoil sections at 70% span

# cl_correction, cd_correction = get_viscous_corrections("corrections.csv")
# cl_alpha0, delta_cl_fun, cd_visc_fun = get_viscous_corrections2("corrections.csv")

# filename = "corrections.csv"
# data = readdlm(filename, ',', skipstart=0)
# @show data
# cls_inv = data[:,1]
# cls_visc = data[:,2]
# cds_visc = data[:,3]
# alphas = data[:,4]

# polar = VortexLattice.Polar(alphas, cls_visc, cds_visc .* 0.0)

# get section_rs
# section_rs = (yle_p1[1:end-1] .+ yle_p1[2:end]) .* 0.5 ./ yle_p1[end]
# blade_files = fill("dji_9443_airfoils.csv", 2)
# polars = VortexLattice.get_polars2([section_rs, section_rs], blade_files)
# @show length(polars[1])
# error()

# or just use the same polar for all sections
# polars = fill(polar, size(system.surfaces, 2))
# polars = [polars, polars]
# polars = nothing

# function plot_polars(polars::Vector{VortexLattice.Polar{TF}}, labels) where TF
#     fig = figure("airfoils")
#     fig.clear()
#     fig.add_subplot(121, xlabel=L"\alpha (^\circ)", ylabel=L"c_l")
#     fig.add_subplot(122, xlabel=L"\alpha (^\circ)", ylabel=L"c_d")
#     axs = fig.get_axes()

#     # loop over polars
#     for (ip,polar) in enumerate(polars)
#         alpha = polar.alphas
#         cl = polar.cls_visc
#         cd = polar.cds_visc
#         @show length(alpha), length(cl), length(cd) alpha cl cd
#         axs[0].plot(alpha, cl, label=labels[ip])
#         axs[1].plot(alpha, cd, label=labels[ip])
#     end
#     axs[0].legend()
#     axs[1].legend()
# end

# plot_polars(polars[1], ["sec$i" for i in 1:length(polars[1])])

monitors = (VortexLattice.ForcesMonitor(length(t_range)),)
benchmark = @elapsed wake = simulate!(system, frames, constant_maneuver!, Uinf, t_range; 
            monitors, name = "rotorhover", 
            # particle_trailing_methods=fill(VortexLattice.NoShed(), length(system.surfaces)),
            # particle_trailing_methods=fill(VortexLattice.OverlapPPS(overlap, p_per_step), length(system.surfaces)),
            particle_trailing_methods=fill(VortexLattice.SigmaOverlap(sigma, overlap), length(system.surfaces)),
            # particle_unsteady_methods=fill(VortexLattice.SigmaOverlap(sigma, overlap), length(system.surfaces)),
            particle_unsteady_methods=fill(VortexLattice.NoShed(), length(system.surfaces)),
            eta = 0.3,
            derivatives = false,
            # vtk_args=(trailing_vortices=false,),
            # wake_args=(SFS=VortexLattice.FLOWVPM.SFS_Cd_twolevel_nobackscatter,),
            # nonlinear_analysis=true,
            # nonlinear_args=(polar_correction=false,),
            # calculate_influence_matrix=true,
            # path=nothing,
            # wake_args=(relaxation=VortexLattice.FLOWVPM.relaxation_none,),
            polars, frames_index = fill(1, length(system.surfaces))
        )

# post-process
Ts = [monitors[1].CF[i][1] for i in 1:length(t_range)]
CTs = Ts ./ (rho * (RPM/60)^2 * (2*R)^4)
# fig = figure("CT")
# fig.clear()
# fig.add_subplot(111, xlabel=L"t", ylabel=L"C_T")
# ax = fig.get_axes()[0]
# ax.plot(collect(t_range), CTs, label="VPM")
# ax.set_ylim(-1.0, 1.0)

# comparison
CT_exp = 0.072
CT_URANS = 0.071
# ax.plot(collect(t_range), fill(CT_exp, length(t_range)), "--", label="experiment")

di = timestep_per_rev * 1
CT_vpm = sum(CTs[end-di+1 : end]) / length(CTs[end-di+1 : end])
percent_error = abs((CT_vpm - CT_exp) / CT_exp) * 100
println("VPM CT: $CT_vpm\nExperiment CT: $CT_exp\nURANS CT: $CT_URANS\nPercent Error (VPM vs Experiment): $percent_error %")

# save csv with CT vs time
name = "rotor_hover_eta0.3_ns20_nc1_nt36_pps4_overlap1.3"
data = hcat(collect(t_range), CTs)
writedlm(name*".csv", data, ',')
