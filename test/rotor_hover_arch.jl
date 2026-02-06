using VortexLattice
using StaticArrays
using PythonPlot
using DelimitedFiles

data_path="./VortexLattice_rotor_data"

constant_maneuver!(frames, system, wake, t) = nothing

R = 0.12
RPM = 5400.0
J = 0.0001
rho             = 1.071778                  # (kg/m^3) air density
mu              = 1.85508e-5                # (kg/ms) air dynamic viscosity
speedofsound    = 342.35                    # (m/s) speed of sound
magVinf         = J*RPM/60*(2*R)
Uinf(t) = SVector{3,Float64}(-1.0, 0.0, 0.0) * magVinf

ns = 13
nc = 1

grids, ratios, sections, invert_normals = VortexLattice.generate_rotor("DJI9443.csv", data_path; 
                                                            turbine_flag=false, 
                                                            clockwise=true, 
                                                            ns, 
                                                            nc, 
                                                            spacing_s=Uniform(), 
                                                            interpolate_airfoils=true);

# generate surface panels for propeller
xle_p1 = [-0.007760952, -0.00912684020509993, -0.01054884338296846, -0.011250349849440629, -0.011772865802853139, -0.012119569829317039, -0.012290082732736044, -0.012268795505995866, -0.012043987625913162, -0.011895522315499299, -0.011461081868185891, -0.010904219636113025, -0.010404176972364522, -0.009799134225632307, -0.00927010352088215, -0.009051792223013484, -0.008616634170407296, -0.008196262641487507, -0.007893374592914745, -0.007768353692954446, -0.007640989960027567, -0.0073486327398695215, -0.007153683499931082, -0.006766626723735409, -0.006228888689944864, -0.00288816]
yle_p1 = [0.004874435999999999, 0.01093368, 0.01699296, 0.0204, 0.02305212, 0.0264, 0.0288, 0.0291114, 0.0324, 0.035170679999999996, 0.04122984, 0.04728912, 0.0533484, 0.05940756, 0.06546684, 0.07152612, 0.07758527999999999, 0.08364456, 0.08970383999999999, 0.095763, 0.10182228, 0.10788155999999999, 0.11394072, 0.11639999999999999, 0.1176, 0.12]
zle_p1 = [0.0017399304062728094, 0.001033148694103248, 5.515520297142833e-5, -0.00039695657142857163, -0.0007488910010571427, -0.0011093118561081603, -0.0013439497122163205, -0.0013743939740463542, -0.0013387395267986484, -0.0012032992967165625, -0.0005335055964156584, 4.7425651801566416e-5, 0.0006065310797702348, 0.0009497390876361914, 0.0012103887940491129, 0.0013988792124036567, 0.0016551317107757106, 0.0019260359731138059, 0.0021145239782349853, 0.0023030049522443858, 0.0024916919017199017, 0.0026830811171171167, 0.002821339876272293, 0.002859589115700257, 0.002878252743800171, 0.0029155799999999996]
chord_p1 = [0.0144, 0.020876759999999998, 0.02597172, 0.028706999999999996, 0.030367079999999998, 0.031616399999999996, 0.03180876, 0.03178992, 0.031150079999999997, 0.030362519999999997, 0.02841924, 0.02627268, 0.02442492, 0.02249088, 0.0207618, 0.01945428, 0.0179046, 0.01653588, 0.015304680000000001, 0.01429284, 0.013367519999999999, 0.01233264, 0.011391324, 0.01069332, 0.009962628, 0.005859876]
theta_p1 = [-0.26234070166665546, -0.31991308290106085, -0.3419082713692357, -0.3454966520785375, -0.3454966520785375, -0.3395274176652338, -0.3327597060800743, -0.3319272068845572, -0.3232074235076193, -0.3148776834780317, -0.29334247725428353, -0.2746396431003036, -0.25713950069538405, -0.23847750818494157, -0.22117737756768568, -0.20365437888236765, -0.18562387512780293, -0.16783277004609531, -0.15107104792086348, -0.1357267208862114, -0.12474499791294327, -0.11647593523729516, -0.11181989297684268, -0.10931562044809875, -0.105924650504665, -0.09440346297697169]
phi_p1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
ns_p1 = ns
nc_p1 = nc
fc_p1 = fill((xc) -> 0, length(yle_p1)) # camberline function for each section
spacing_s_p1 = Sine()
spacing_c_p1 = Uniform()
mirror_p1 = false
p1grid1, p1ratio1 = wing_to_grid(xle_p1, yle_p1, zle_p1, chord_p1, theta_p1, phi_p1, ns_p1, nc_p1;
    mirror=mirror_p1, fc=fc_p1, spacing_s=spacing_s_p1, spacing_c=spacing_c_p1)
p1grid2, p1ratio2 = wing_to_grid(xle_p1, yle_p1, zle_p1, chord_p1, theta_p1, phi_p1, ns_p1, nc_p1;
mirror=mirror_p1, fc=fc_p1, spacing_s=spacing_s_p1, spacing_c=spacing_c_p1)
translate!(p1grid1, SVector{3}(-chord_p1[1]*0.5, 0.0, 0.0))
translate!(p1grid2, SVector{3}(-chord_p1[1]*0.5, 0.0, 0.0))
R1_b1 = VortexLattice.Rodrigues(SVector{3}(0.0, 1.0, 0.0), -pi*0.5)
rotate!(p1grid1, R1_b1)
R1_b2 = VortexLattice.Rodrigues(SVector{3}(1.0,0,0), pi*1.0) * R1_b1
rotate!(p1grid2, R1_b2)

grids = [p1grid1, p1grid2]
ratios = [p1ratio1, p1ratio2]

system = System(grids; ratios, sections);

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

steady_analysis!(system, system.reference[], system.freestream[]; symmetric=false)

write_vtk("rotor_hover_initial", system; write_wakes=false, trailing_edge_list=fill(false, length(system.surfaces)))

frames = ReferenceFrame(system;
        origin = SVector{3}(0.0, 0.0, 0.0),
        v = SVector{3}(0.0, 0.0, 0.0),
        ω_axis = SVector{3}(1.0, 0.0, 0.0),
        ω = -RPM * 2 * pi / 60,
        R = SMatrix{3,3}(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0),
        name = "vehicle",
        child_index = Int[],
        dependent_index = collect(1:length(system.surfaces))
    )


n_revs = 10
ttot = n_revs / (RPM / 60)
timestep_per_rev = 36
t_range = range(start=0.0, stop=ttot, length=n_revs * timestep_per_rev + 1)
overlap = 1.3
p_per_step = 4
nsteps_per_rev = length(t_range) / n_revs
sigma = overlap * 2*pi*R / (nsteps_per_rev*p_per_step)

monitors = (VortexLattice.ForcesMonitor(length(t_range)),)
benchmark = @elapsed wake = simulate!(system, frames, constant_maneuver!, Uinf, t_range; 
            monitors, name = "rotorhover", 
            # particle_trailing_methods=fill(VortexLattice.NoShed(), length(system.surfaces)),
            # particle_trailing_methods=fill(VortexLattice.OverlapPPS(overlap, p_per_step), length(system.surfaces)),
            particle_trailing_methods=fill(VortexLattice.SigmaOverlap(sigma, overlap), length(system.surfaces)),
            # particle_unsteady_methods=fill(VortexLattice.SigmaOverlap(sigma, overlap), length(system.surfaces)),
            particle_unsteady_methods=fill(VortexLattice.NoShed(), length(system.surfaces)),
            eta = 0.4,
            derivatives = false,
            vtk_args=(trailing_vortices=false,),
            # nonlinear_analysis=true,
            # nonlinear_args=(polar_correction=false,),
            # calculate_influence_matrix=true,
            # path=nothing,
            # wake_args=(relaxation=VortexLattice.FLOWVPM.relaxation_none,)
        )

# post-process
Ts = [-monitors[1].CF[i][1] for i in 1:length(t_range)]
CTs = Ts ./ (rho * (RPM/60)^2 * (2*R)^4)
fig = figure("CT")
fig.clear()
fig.add_subplot(111, xlabel=L"t", ylabel=L"C_T")
ax = fig.get_axes()[0]
ax.plot(collect(t_range), CTs, label="VPM")
# ax.set_ylim(-1.0, 1.0)

# comparison
CT_exp = 0.072
CT_URANS = 0.071
ax.plot(collect(t_range), fill(CT_exp, length(t_range)), "--", label="experiment")

di = timestep_per_rev >> 1
CT_vpm = sum(CTs[end-di: end]) / length(CTs[end-di: end])
percent_error = abs((CT_vpm - CT_exp) / CT_exp) * 100
println("VPM CT: $CT_vpm\nExperiment CT: $CT_exp\nURANS CT: $CT_URANS\nPercent Error (VPM vs Experiment): $percent_error %")

# save csv with CT vs time
name = "rotor_hover_eta0.4_ns13cos_nc1"
data = hcat(collect(t_range), CTs)
writedlm(name*".csv", data, ',')