
"""
    generate_rotor(rotor_file::String, data_path="."; optargs...)

Generate a grids, ratios from a rotor file. Explained in the docs:

 # Arguments
 - `rotor_file`: File with rotor parameters
 - `data_path`: Path to the rotor_file folder

 # Keyword Arguments
 - `interpolate_airfoils`: If true, the airfoil polars are interpolated linearly between the airfoils.
    Defaults to false.
 - `zero_at_root`: If true, the airfoil positions are defined from the root,
    otherwise from the center of the hub. Defaults to false.
 - `polar_in_radians`: If true, the airfoil polars are in radians, otherwise in degrees.
    Defaults to false.
 - `turbine_flag`: If true, the rotor is a turbine, otherwise a propeller. Defaults to false.
    If the true then the angle of attack is inverted in the nonlinear solver.
 - `clockwise`: If true, the rotor is intended to rotate clockwise, otherwise counter-clockwise. Defaults to false.
 - `ns`: Number of spanwise panels. Defaults to 10.
 - `nc`: Number of chordwise panels. Defaults to 1.
 - `spacing_s`: Spanwise spacing. Defaults to `VortexLattice.Sine()`.
 - `spacing_c`: Chordwise spacing. Defaults to `VortexLattice.Uniform()`.
 - `initial_azimuthal_angle`: Initial azimuthal angle of the rotor. Defaults to 0.0.
"""
function generate_rotor(rotor_file::String, data_path; optargs...)
    Rtip, Rhub, B, blade_file = _read_rotor(rotor_file, data_path)

    return _generate_rotor(Rtip, Rhub, B, blade_file,
                            data_path; optargs...)
end

function _read_rotor(rotor_file::String, data_path)

    # Path to rotor files
    rotor_path = joinpath(data_path, "rotors")

    data = readdlm(joinpath(rotor_path, rotor_file),',';skipstart=1)[:,2]
    Rtip = data[1]
    Rhub = data[2]
    B = Int64(data[3])
    blade_file = String(data[4])

    return Rtip, Rhub, B, blade_file
end

function _read_blade(blade_file::String, data_path)

    # Path to rotor files
    rotor_path = joinpath(data_path, "rotors")

    # Read blade
    files = readdlm(joinpath(rotor_path, blade_file),',';skipstart=1)

    num_files = size(files, 1)

    chorddist = readdlm(joinpath(rotor_path, files[1, 2]),',';skipstart=1)
    pitchdist = readdlm(joinpath(rotor_path, files[2, 2]),',';skipstart=1)
    sweepdist = readdlm(joinpath(rotor_path, files[3, 2]),',';skipstart=1)
    heightdist = readdlm(joinpath(rotor_path, files[4, 2]),',';skipstart=1)
    airfoil_files = readdlm(joinpath(rotor_path, files[5, 2]),',';skipstart=1)
    if num_files > 5
        if !isa(files[6,2], Number)
            airfoil_reference = readdlm(joinpath(rotor_path, files[6, 2]),',';skipstart=1)
        else
            airfoil_reference = zeros(1,1)
        end
    else
        airfoil_reference = zeros(1,1)
    end

    af = airfoil_files
    airfoil_files = [(Float64(af[i, 1]), String(af[i, 2]), String(af[i, 3]))
                        for i in axes(af, 1)]

    return chorddist, pitchdist, sweepdist, heightdist, airfoil_files, airfoil_reference
end

function _generate_rotor(Rtip::Real, Rhub::Real, B::Int, blade_file::String,
    data_path; optargs...)

    (chorddist, pitchdist, sweepdist, heightdist,
    airfoil_files, airfoil_reference) = _read_blade(blade_file, data_path)

    return _generate_rotor(Rtip, Rhub, B, chorddist, pitchdist, sweepdist,
                            heightdist, airfoil_files, airfoil_reference,
                            data_path; optargs...)
end

function _generate_rotor(Rtip, Rhub, B::Int,
    chorddist,
    pitchdist,
    sweepdist,
    heightdist,
    airfoil_files::Array{Tuple{TF,String,String},1},
    airfoil_reference,
    data_path;
    # INPUT OPTIONS
    optargs...) where TF

    # Read airfoil contours
    # Airfoils along the blade as
    # airfoil_contours=[ (pos1, contour1, polar1), (pos2, contour2, pol2), ...]
    # with contour=(x,y) and pos the position from root to tip between 0 and 1.
    # pos1 must equal 0 (root airfoil) and the last must be 1 (tip airfoil)
    airfoil_contours = Tuple{TF,Array{TF, 2},String}[]
    airfoil_path = joinpath(data_path, "airfoils")
    for (r, rfl_file, clcurve_file) in airfoil_files
        contour = readdlm(joinpath(airfoil_path,rfl_file),',';skipstart=1)
        rfl = Array{TF,2}(undef,size(contour,1),2)
        rfl[:,1] .= contour[:,1]
        rfl[:,2] .= contour[:,2]

        push!(airfoil_contours, (r, rfl, clcurve_file))
    end

    return _generate_rotor(Rtip, Rhub, B,
            chorddist, pitchdist, sweepdist, heightdist,
            airfoil_contours, airfoil_reference,
            data_path; optargs...)
end

function _generate_rotor(Rtip, Rhub, B::Int,
    chorddist,
    pitchdist,
    sweepdist,
    heightdist,
    airfoil_contours,
    airfoil_reference,
    data_path;
    interpolate_airfoils=false,
    # INPUT OPTIONS
    zero_at_root=false,
    polar_in_radians=false,
    # PROCESSING OPTIONS
    turbine_flag=false,
    clockwise=false,
    ns=10, nc=1,
    spacing_s=VortexLattice.Sine(),
    spacing_c=VortexLattice.Uniform(),
    RPM=0.0,
    rotor_name="rotor",
    frames = Vector{VortexLattice.ReferenceFrame{Float64}}(undef, 0),
    surface_index = collect(1:B),
    parent_index = -1) # radians

    yle = chorddist[:,1] .* Rtip
    if zero_at_root
        yle = (Rhub .+ yle .* (Rtip-Rhub))
    end

    chord = FLOWMath.linear(chorddist[:,1] .* Rtip, chorddist[:,2] .* Rtip, yle)
    theta = .-FLOWMath.linear(pitchdist[:,1] .* Rtip, deg2rad.(pitchdist[:,2]), yle)
    xle = .-FLOWMath.linear(sweepdist[:,1] .* Rtip, sweepdist[:,2] .* Rtip, yle)
    zle = .-FLOWMath.linear(heightdist[:,1] .* Rtip, heightdist[:,2] .* Rtip, yle)

    if size(airfoil_reference,1) < length(yle)
        if length(airfoil_reference) != 1
            @warn "Airfoil reference line has different length than spanwise stations. Ignoring reference line."
        end
        airfoil_reference = zeros(length(yle),2)
    end

    grid, ratio = wing_to_grid(xle,yle,zle,chord,theta,zeros(length(yle)),
                    ns,nc;reference_line=airfoil_reference,
                    spacing_s=spacing_s, spacing_c=spacing_c)

    translate!(grid, SVector{3}( -chord[1]*0.5,0.0, 0.0))
    R = VortexLattice.Rodrigues(SVector{3}(0.0, 1.0, 0.0), -pi*0.5)
    VortexLattice.rotate!(grid, R)

    grids = Vector{typeof(grid)}(undef,B)
    ratios = Vector{typeof(ratio)}(undef,B)

    grids[1] = deepcopy(grid)
    ratios[1] = deepcopy(ratio)

    diff_angle = 2π/B
    R1 = VortexLattice.Rodrigues(SVector{3}(1.0, 0.0, 0.0), diff_angle)
    for i = 2:B
        grids[i] = deepcopy(grids[i-1])
        VortexLattice.rotate!(grids[i], R1)
    end

    surfaces = Vector{Matrix{SurfacePanel{Float64}}}(undef,B)
    for i = eachindex(grids)
        grids[i], ratios[i], surfaces[i] = grid_to_surface_panels(grids[i]; ratios = ratio)
    end

    airfoils = Vector{Tuple{Float64, CCBlade.AlphaAF{Float64, String, Akima{Vector{Float64}, Vector{Float64}, Float64}}}}(undef,length(airfoil_contours))
    contours = Vector{Array{Float64,2}}(undef,length(airfoil_contours))
    for (rfli, (pos, contour, file_name)) in enumerate(airfoil_contours)
        polar = get_polars(joinpath(data_path, "airfoils", file_name); radians=polar_in_radians)

        if zero_at_root
            pos = (Rhub + pos*(Rtip-Rhub))/Rtip
        end

        airfoils[rfli] = (pos*Rtip, polar)
        contours[rfli] = contour
    end
    airfoils, contours = redo_airfoils(airfoils,contours,surfaces[1]; interpolate=interpolate_airfoils)
    polar = Vector{Polar{eltype(airfoils[1].alpha)}}(undef, length(airfoils))
    for i in eachindex(airfoils)
        polar[i] = Polar(rad2deg.(airfoils[i].alpha), airfoils[i].cl, airfoils[i].cd)
    end
    polars = Vector{typeof(polar)}(undef,B)
    for i in eachindex(polars)
        polars[i] = deepcopy(polar)
    end

    frames = add_rotor_frames!(frames, rotor_name, parent_index, surface_index, RPM, B)

    return grids, ratios, polars, frames
end

function add_rotor_frames!(frames::Vector{ReferenceFrame{TF}}, rotor_name, parent_index, surface_index, RPM, B) where TF
    origin = SVector{3}(0.0, 0.0, 0.0)

    frame = ReferenceFrame(origin,
            zero(SVector{3,TF}), # v
            SVector{3,TF}(1.0, 0.0, 0.0), # ω_axis
            -RPM * 2 * pi / 60, # ω
            SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0), # R
            SMatrix{3,3}(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0), # Rp2g
            rotor_name, # name
            parent_index, # parent_index
            Int[], # child_indices
            surface_index, # surface_index
    )
    push!(frames, frame)
    rotor_index = length(frames)

    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    Rot = VortexLattice.Rodrigues(SVector{3}(1.0, 0.0, 0.0), 2π/B)

    for i in 1:B
        add_frame!(frames, "$rotor_name-blade-$i", rotor_index, origin, [surface_index[i]];
            R = R,
        )
        R = R * Rot
    end
    return frames
end

function redo_airfoils(airfoils, contours, surface; interpolate=false)
    nc, ns = size(surface)
    new_airfoils = Vector{CCBlade.AlphaAF{Float64, String, Akima{Vector{Float64}, Vector{Float64}, Float64}}}(undef,ns)
    new_contours = Vector{eltype(contours)}(undef,ns)
    r = zeros(3)
    for i in 1:ns
        r .= 0.0
        for j in 1:nc
            r .+= surface[j,i].rcp
        end
        r = r ./ nc
        radius = norm(r)

        index = select_airfoil_index(radius, airfoils)

        if !interpolate
            new_airfoils[i] = airfoils[index][2]
            new_contours[i] = contours[index]
        else
            interpolate_airfoil!(new_airfoils, new_contours, airfoils, contours, radius, index, i)
        end
    end
    return new_airfoils, new_contours
end

function select_airfoil_index(radius, airfoils)
    for i in eachindex(airfoils)
        if radius < airfoils[i][1]
            return i-1
        end
    end
    return length(airfoils)
end

function interpolate_airfoil!(new_airfoils, new_contours, airfoils, contours, r, index, i)
    if index >= length(airfoils)
        new_airfoils[i] = airfoils[index][2]
        new_contours[i] = contours[index]
        return
    end
    rs = [airfoils[index][1], airfoils[index+1][1]]
    airfoil1 = airfoils[index][2]
    airfoil2 = airfoils[index+1][2]
    alpha_l = max(minimum(airfoil1.alpha), minimum(airfoil2.alpha))
    alpha_h = min(maximum(airfoil1.alpha), maximum(airfoil2.alpha))
    alpha1 = airfoil1.alpha[alpha_h .>= airfoil1.alpha .>= alpha_l]
    alpha2 = airfoil2.alpha[alpha_h .>= airfoil2.alpha .>= alpha_l]
    alpha = sort(unique(union(alpha1,alpha2)))
    cl = zeros(length(alpha))
    cd = zeros(length(alpha))
    coefs = zeros(2)

    for i in eachindex(alpha)
        coefs[1] = airfoil1.clspline(alpha[i])
        coefs[2] = airfoil2.clspline(alpha[i])
        cl[i] = FLOWMath.akima(rs,coefs,r)

        coefs[1] = airfoil1.cdspline(alpha[i])
        coefs[2] = airfoil2.cdspline(alpha[i])
        cd[i] = FLOWMath.akima(rs,coefs,r)
    end

    new_airfoils[i] = CCBlade.AlphaAF(alpha, cl, cd, "Interpolated Airfoil from airfoils $(index) and $(index+1)")
    new_contours[i] = contours[index]
end

function get_polars(filename; radians=true)
    if endswith(lowercase(filename), ".csv")
        data = readdlm(filename,',';skipstart=1)
        alpha = data[:,1]
        cl = data[:,2]
        cd = data[:,3]
        if !radians
            alpha *= pi/180
        end
        return CCBlade.AlphaAF(alpha, cl, cd, filename)
    else
        return CCBlade.AlphaAF(filename; radians=radians)
    end
end