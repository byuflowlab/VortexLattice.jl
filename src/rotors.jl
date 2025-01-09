
"""
    generate_rotor(rotor_file::String; data_path=".",optargs...)

Generate a rotor from a rotor file. The rotor file follows the same format as FLOWUnsteady.jl:
"""
function generate_rotor(rotor_file::String; data_path=".",optargs...)
    Rtip, Rhub, B, blade_file = read_rotor(rotor_file; data_path=data_path)

    return generate_rotor(Rtip, Rhub, B, blade_file;
                            data_path=data_path, optargs...)
end

function read_rotor(rotor_file::String; data_path=".")

    # Path to rotor files
    rotor_path = joinpath(data_path, "rotors")

    data = readdlm(joinpath(rotor_path, rotor_file),',';skipstart=1)[:,2]
    Rtip = data[1]
    Rhub = data[2]
    B = Int64(data[3])
    blade_file = String(data[4])

    return Rtip, Rhub, B, blade_file
end

function read_blade(blade_file::String; data_path=def_data_path)

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
    airfoil_reference = zeros(length(airfoil_files),2)
        if num_files > 5
            if tryparse(Float64, files[6, 2]) !== nothing # check if files[6,2] is a Float
            else
                airfoil_reference = readdlm(joinpath(rotor_path, files[6, 2]),',';skipstart=1)
            end
        end

    af = airfoil_files
    airfoil_files = [(Float64(af[i, 1]), String(af[i, 2]), String(af[i, 3]))
                        for i in axes(af, 1)]

    return chorddist, pitchdist, sweepdist, heightdist, airfoil_files, airfoil_reference
end

function generate_rotor(Rtip::Real, Rhub::Real, B::Int, blade_file::String;
    data_path=def_data_path, optargs...)

    (chorddist, pitchdist, sweepdist, heightdist,
    airfoil_files, airfoil_reference) = read_blade(blade_file; data_path=data_path)

    return generate_rotor(Rtip, Rhub, B, chorddist, pitchdist, sweepdist,
                            heightdist, airfoil_files, airfoil_reference;
                            data_path=data_path, optargs...)
end

function generate_rotor(Rtip, Rhub, B::Int,
    chorddist,
    pitchdist,
    sweepdist,
    heightdist,
    airfoil_files::Array{Tuple{TF,String,String},1},
    airfoil_reference;
    # INPUT OPTIONS
    data_path=def_data_path, optargs...) where TF

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

    return generate_rotor(Rtip, Rhub, B,
            chorddist, pitchdist, sweepdist, heightdist,
            airfoil_contours, airfoil_reference;
            data_path=data_path, optargs...)
end

function generate_rotor(Rtip, Rhub, B::Int,
    chorddist,
    pitchdist,
    sweepdist,
    heightdist,
    airfoil_contours,
    airfoil_reference;
    # INPUT OPTIONS
    zero_at_root=false,
    polar_in_radians=false,
    data_path=default_database,
    # PROCESSING OPTIONS
    turbine_flag=false,
    clockwise=false,
    ns=10, nc=1,
    spacing_s=VortexLattice.Sine(),
    spacing_c=VortexLattice.Uniform(),
    initial_azimuthal_angle = 0.0, # radians
    # OUTPUT OPTIONS
    verbose=false, v_lvl=1)

    if verbose; println("\t"^v_lvl*"Generating geometry..."); end;

    turbine_mod = (-1)^turbine_flag # if turbine_flag then -1
    if turbine_flag
        if clockwise
            clockwise = false
        else
            clockwise = true
        end
    end

    yle = chorddist[:,1] .* Rtip
    chord = chorddist[:,2] .* Rtip
    theta = deg2rad.(pitchdist[:,2]) * turbine_mod
    xle = sweepdist[:,2] .* Rtip
    zle = heightdist[:,2] .* Rtip

    grid, ratio = wing_to_grid(xle,yle,zle,chord,theta,zeros(length(yle)),
                    ns,nc;reference_line=airfoil_reference,
                    spacing_s=spacing_s, spacing_c=spacing_c)

    grids = Vector{typeof(grid)}(undef,B)
    ratios = Vector{typeof(ratio)}(undef,B)
    Rotate_grid!(grid, -pi/2 * turbine_mod, 2)

    if clockwise
        reverse!(grid,dims=3)
        grid[2,:,:] = -grid[2,:,:]
        Rotate_grid!(grid, pi, 1)
    end

    Rotate_grid!(grid, initial_azimuthal_angle, 1)
    grids[1] = deepcopy(grid)
    ratios[1] = deepcopy(ratio)

    diff_angle = 2π/B
    for i = 2:B
        Rotate_grid!(grid, diff_angle, 1)
        grids[i] = deepcopy(grid)
    end

    surfaces = Vector{Matrix{SurfacePanel{Float64}}}(undef,B)
    for i = eachindex(grids)
        grids[i], ratios[i], surfaces[i] = grid_to_surface_panels(grids[i]; ratios = ratio)
    end

    if verbose; println("\t"^v_lvl*"Generating airfoil data..."); end;

    airfoils = Vector{Tuple{Float64, CCBlade.AlphaAF{Float64, String, Akima{Vector{Float64}, Vector{Float64}, Float64}}}}(undef,length(airfoil_contours))
    contours = Vector{Array{Float64,2}}(undef,length(airfoil_contours))
    for (rfli, (pos, contour, file_name)) in enumerate(airfoil_contours)

        if verbose; println("\t"^(v_lvl+1)*"$file_name"); end;
        polar = CCBlade.AlphaAF(joinpath(data_path, "airfoils", file_name); radians=polar_in_radians)

        if zero_at_root
            pos = (Rhub + pos*(Rtip-Rhub))/Rtip
        end

        airfoils[rfli] = (pos*Rtip, polar)
        contours[rfli] = contour
    end
    airfoils, contours = redo_airfoils(airfoils,contours,surfaces[1])

    section = grid_to_sections(grids[1], airfoils; ratios=ratios[1], contours)
    sections = Vector{typeof(section)}(undef,B)
    sections[1] = deepcopy(section)
    for i = 2:B
        sections[i] = deepcopy(section)
    end
    
    return grids, ratios, sections
end

function Rotate_grid!(grid::Array{Float64,3}, angle, axis::Int)
    if angle == 0.0
        return grid
    end
    R = RotationMatrix(angle, axis)
    for k = axes(grid,3)
        for j = axes(grid,2)
            grid[:,j,k] = R*grid[:,j,k]
        end
    end
    return grid
end

function RotationMatrix(angle, axis::Int)
    st, ct = sincos(angle)
    if axis == 1
        return [1 0 0; 0 ct -st; 0 st ct]
    elseif axis == 2
        return [ct 0 st; 0 1 0; -st 0 ct]
    elseif axis == 3
        return [ct -st 0; st ct 0; 0 0 1]
    else
        error("Invalid axis")
    end
end

function MirrorGrid!(grid::Array{Float64,3}, axis::Int)
    if axis == 1
        grid[1,:,:] = -grid[1,:,:]
    elseif axis == 2
        grid[2,:,:] = -grid[2,:,:]
    elseif axis == 3
        grid[3,:,:] = -grid[3,:,:]
    else
        error("Invalid axis")
    end
    return grid
end

function redo_airfoils(airfoils, contours, surface)
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

        new_airfoils[i] = airfoils[index][2]
        new_contours[i] = contours[index]
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