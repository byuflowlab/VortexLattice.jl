using DelimitedFiles
using Xfoil
using PythonPlot
using FLOWMath

#=
cst_u, cst_l = cst_foil_fit(x, yu, x, yl, n_cst=10, xn1=0.3, xn2=1.0)
x, yu, yl, tmax, rLE = cst_foil(1001, cst_u, cst_l, x=None, t=None, tail=0.0)   # black
x, yu, yl, tmax, rLE = cst_foil(1001, cst_u, cst_l, x=None, t=None, tail=0.02)  # red
x, yu, yl, tmax, rLE = cst_foil(1001, cst_u, cst_l, x=None, t=0.05, tail=0.0)   # blue
x, yu, yl, tmax, rLE = cst_foil(1001, cst_u, cst_l, x=None, t=0.05, tail=0.01)  # green
=#

function airfoil_geometry(contour_file; airfoil_path=joinpath(@__DIR__, "..", "VortexLattice_rotor_data", "airfoils"))
    data = readdlm(joinpath(airfoil_path, contour_file), ',', skipstart=1)
    x = data[:,1]
    y = data[:,2]
    return x, y
end

function get_camber(x, y)
    i_upper = findfirst(v -> isapprox(v, 0.0), x)
    x_upper = reverse(x[1:i_upper])
    y_upper = reverse(y[1:i_upper])
    x_lower = x[i_upper:end]
    y_lower = y[i_upper:end]
    y_l = FLOWMath.linear(x_lower, y_lower, x_upper)
    camber = 0.5 .* (y_upper + y_l)
    return x_upper, camber
end

# extract contour
rfl_contour = "DJI9443-airfoilsec2.csv"
x, y = airfoil_geometry(rfl_contour)

fig = figure("contour")
fig.clear()
ax = fig.add_subplot(111, xlabel=L"x/c", ylabel=L"y/c")
ax.plot(x, y, label="contour")
ax.set_aspect("equal", adjustable="box")
# ax.legend()

set_coordinates(x, y)
pane()

# quick sanity check of viscous vs. inviscid:
alpha = 5.0
cl_inv, cm_inv = solve_alpha(alpha; mach=0.0)

# Re = 13131.0

# check Re
R = 0.12
RPM = 5400.0

rho = 1.071778                  # (kg/m^3) air density
mu = 1.85508e-5                # (kg/ms) air dynamic viscosity
xle_p1 = [-0.007760952, -0.00912684020509993, -0.01054884338296846, -0.011250349849440629, -0.011772865802853139, -0.012119569829317039, -0.012290082732736044, -0.012268795505995866, -0.012043987625913162, -0.011895522315499299, -0.011461081868185891, -0.010904219636113025, -0.010404176972364522, -0.009799134225632307, -0.00927010352088215, -0.009051792223013484, -0.008616634170407296, -0.008196262641487507, -0.007893374592914745, -0.007768353692954446, -0.007640989960027567, -0.0073486327398695215, -0.007153683499931082, -0.006766626723735409, -0.006228888689944864, -0.00288816]
yle_p1 = [0.004874435999999999, 0.01093368, 0.01699296, 0.0204, 0.02305212, 0.0264, 0.0288, 0.0291114, 0.0324, 0.035170679999999996, 0.04122984, 0.04728912, 0.0533484, 0.05940756, 0.06546684, 0.07152612, 0.07758527999999999, 0.08364456, 0.08970383999999999, 0.095763, 0.10182228, 0.10788155999999999, 0.11394072, 0.11639999999999999, 0.1176, 0.12]
zle_p1 = [0.0017399304062728094, 0.001033148694103248, 5.515520297142833e-5, -0.00039695657142857163, -0.0007488910010571427, -0.0011093118561081603, -0.0013439497122163205, -0.0013743939740463542, -0.0013387395267986484, -0.0012032992967165625, -0.0005335055964156584, 4.7425651801566416e-5, 0.0006065310797702348, 0.0009497390876361914, 0.0012103887940491129, 0.0013988792124036567, 0.0016551317107757106, 0.0019260359731138059, 0.0021145239782349853, 0.0023030049522443858, 0.0024916919017199017, 0.0026830811171171167, 0.002821339876272293, 0.002859589115700257, 0.002878252743800171, 0.0029155799999999996]
chord_p1 = [0.0144, 0.020876759999999998, 0.02597172, 0.028706999999999996, 0.030367079999999998, 0.031616399999999996, 0.03180876, 0.03178992, 0.031150079999999997, 0.030362519999999997, 0.02841924, 0.02627268, 0.02442492, 0.02249088, 0.0207618, 0.01945428, 0.0179046, 0.01653588, 0.015304680000000001, 0.01429284, 0.013367519999999999, 0.01233264, 0.011391324, 0.01069332, 0.009962628, 0.005859876]
theta_p1 = [-0.26234070166665546, -0.31991308290106085, -0.3419082713692357, -0.3454966520785375, -0.3454966520785375, -0.3395274176652338, -0.3327597060800743, -0.3319272068845572, -0.3232074235076193, -0.3148776834780317, -0.29334247725428353, -0.2746396431003036, -0.25713950069538405, -0.23847750818494157, -0.22117737756768568, -0.20365437888236765, -0.18562387512780293, -0.16783277004609531, -0.15107104792086348, -0.1357267208862114, -0.12474499791294327, -0.11647593523729516, -0.11181989297684268, -0.10931562044809875, -0.105924650504665, -0.09440346297697169]

chord_07rR = FLOWMath.linear(yle_p1, chord_p1, 0.7*R)
v_07rR = 2*pi*RPM/60*0.7*R
Re_07rR = rho * chord_07rR * v_07rR / mu

cl_visc, cd_visc, cdp_visc, cm_visc, conv = solve_alpha(alpha, Re_07rR; mach=0.0, iter=50, ncrit=1, reinit=false, xtrip=(1.0,1.0))

# @show cl_inv, cl_visc, cd_visc

# generate full polars
alpha_range = range(-10.0, stop=15.0, length=41)
cls_inv, cms_inv = alpha_sweep(x, y, alpha_range; mach=0.0, npan=140, printdata=false, filename=nothing, zeroinit=true)
cls_visc, cds_visc, cdps_visc, cms_visc, convs_visc = alpha_sweep(x, y, alpha_range, Re_07rR; mach=0.0, iter=50, ncrit=1, reinit=false, xtrip=(1.0,1.0), printdata=false, filename=nothing, clmaxstop=false, clminstop=false)

fig = figure("polar")
fig.clear()
ax = fig.add_subplot(121, xlabel=L"\alpha", ylabel=L"c_l")
ax.plot(alpha_range, cls_inv, label="inviscid")
ax.plot(alpha_range, cls_visc, label="viscous")
ax.legend()
ax2 = fig.add_subplot(122, xlabel=L"\alpha", ylabel="correction factor")
ax2.plot(alpha_range, cls_visc ./ cls_inv)

# generate correction functions
function get_viscous_corrections(filename)
    data = readdlm(filename, ',', skipstart=0)
    # @show data
    cls_inv = data[:,1]
    cls_visc = data[:,2]
    cds_visc = data[:,3]
    cl_ratios = cls_visc ./ cls_inv
    cd_ratios = cds_visc ./ cls_inv
    cl_correction = (cl) -> begin
        if cl > 0.125
            val = FLOWMath.linear(cls_inv, cl_ratios, cl)
            return clamp(val, 0.0, 1.0)
        else
            return 1.0
        end
    end

    cd_correction = (cl) -> begin
        return 0.0
        # val = FLOWMath.linear(cls_inv, cd_ratios, cl)
        # return clamp(val, -1.0, 1.0)
    end
    return cl_correction, cd_correction
end

function get_viscous_corrections2(filename)
    data = readdlm(filename, ',', skipstart=0)
    # @show data
    cls_inv = data[:,1]
    cls_visc = data[:,2]
    cds_visc = data[:,3]
    alphas = data[:,4]

    # get alpha=0 cl
    cl_inv_func(alpha) = FLOWMath.linear(alphas, cls_inv, alpha)
    cl_alpha0 = cl_inv_func(0.0)

    # get cl correction factor as function of cl
    cl_visc_func = (alpha) -> FLOWMath.linear(alphas, cls_visc, alpha)
    delta_cl_func = (alpha) -> cl_visc_func(alpha) - cl_inv_func(alpha)

    # cd as a function of cl
    cd_visc_func = (cl_v) -> 0.0 # FLOWMath.linear(cls_visc, cds_visc, cl_v)

    return cl_alpha0, delta_cl_func, cd_visc_func
end

function write_polar(af_file, contour_file, Re; alpha_range=range(-10,stop=15,length=51),
        airfoil_path=joinpath(@__DIR__, "..", "VortexLattice_rotor_data", "airfoils"),
        ncrit=1, xfoil_args...
    )
    # read contour
    x, y = airfoil_geometry(contour_file; airfoil_path)

    # run xfoil
    cls_visc, cds_visc, cdps_visc, cms_visc, convs_visc = 
        alpha_sweep(x, y, alpha_range, Re; ncrit, xfoil_args...)

    if !prod(convs_visc)
        @warn "Some Xfoil runs did not converge for $af_file"
    end
    
    # write to file
    data = Matrix{Any}(undef, length(cls_visc)+1, 4)
    data[1,:] .= ["alpha", "cl", "cd", "cm"]
    data[2:end, 1] .= alpha_range
    data[2:end, 2] .= cls_visc
    data[2:end, 3] .= cds_visc
    data[2:end, 4] .= cms_visc
    writedlm(joinpath(airfoil_path, af_file), data, ',')
end

function write_polars(rs, chords, RPM, af_files, contour_files; 
        rho = 1.071778, mu = 1.85508e-5,
        alpha_range=range(-10,stop=15,length=51), 
        airfoil_path=joinpath(@__DIR__, "..", "VortexLattice_rotor_data", "airfoils"),
        ncrit=1, xfoil_args...
    )
    @assert length(rs) == length(chords)
    @assert length(chords) == length(af_files)
    @assert length(af_files) == length(contour_files)
    for i in eachindex(rs)
        r = rs[i]
        c = chords[i]
        v = 2 * pi * RPM / 60 * r
        Re = rho * c * v / mu
        write_polar(af_files[i], contour_files[i], Re; alpha_range, airfoil_path, ncrit, xfoil_args...)
    end
end

write_polar("testPolar.txt", rfl_contour, Re_07rR)

# get rs for known sections
Rhub = 0.00624
Rtip = 0.12
rs = [0.0, 0.0857143, 0.185714, 0.371429, 0.714286, 0.942857, 1.0] .* (Rtip - Rhub) .+ Rhub
chords = FLOWMath.linear(yle_p1, chord_p1, rs)
RPM = 5400.0
af_files = ["dji_9443_$i.csv" for i in 1:length(rs)]
contour_files = ["DJI9443-airfoilsec6.csv", "DJI9443-airfoilsec6.csv", "DJI9443-airfoilsec4.csv", "DJI9443-airfoilsec3.csv", "DJI9443-airfoilsec2.csv", "DJI9443-airfoilsec1.csv", "DJI9443-airfoilsec1.csv"]
write_polars(rs, chords, RPM, af_files, contour_files)

cl_correction, cd_correction = get_viscous_corrections("corrections.csv")

fig2 = figure("splines")
fig2.clear()
ax = fig2.add_subplot(121, xlabel=L"c_l", ylabel="lift correction")
cl_range = range(minimum(cls_inv), stop=maximum(cls_inv), length=100)
ax.plot(cl_range, cl_correction.(cl_range), label="cl correction")
ax = fig2.add_subplot(122, xlabel=L"c_l", ylabel="drag correction")
ax.plot(cl_range, cd_correction.(cl_range), label="cd correction")
ax.legend()

writedlm("corrections.csv", hcat(cls_inv, cls_visc, cds_visc, collect(alpha_range)), ',')