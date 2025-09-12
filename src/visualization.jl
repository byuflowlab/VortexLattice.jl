# TODO: Add grid based visualization
"""
    write_vtk(name, system; write_surfaces = true, write_wakes = false, kwargs...)

Write geometry from surfaces and/or wakes to Paraview files for visualization.

# Arguments
 - `name`: Base name for the generated files
 - `system`: System object containing surfaces and/or wakes

# Keyword Arguments:
    - `write_surfaces = true`: Flag indicating whether to write surface geometry
    - `write_wakes = false`: Flag indicating whether to write wake geometry
    - `symmetric`: (required if `surface_properties` is provided) Flags indicating whether a
        mirror image (across the X-Z plane) was used when calculating induced velocities
        for each surface.
    - `trailing_vortices`: Flag indicating whether the model uses trailing vortices.
        Defaults to `true` when wake panels are absent, `false` otherwise
    - `xhat`: Direction in which trailing vortices extend if used. Defaults to [1, 0, 0].
    - `wake_length`: Distance to extend trailing vortices. Defaults to 10
    - `metadata`: Dictionary of metadata to include in generated files
"""


function write_vtk(name::String, system::System; write_surfaces = true, write_wakes = false, xhat = system.xhat[], 
    trailing_edge_list=fill(true, length(system.surfaces)), kwargs...)

    if write_surfaces && write_wakes
        write_vtk(name, system.surfaces, system.wakes, system.properties; trailing_edge_list, symmetric=system.symmetric, kwargs...)
    elseif write_surfaces
        write_vtk(name, system.surfaces, system.properties; trailing_edge_list, symmetric=system.symmetric, xhat, kwargs...)
    elseif write_wakes
        write_vtk(name, system.wakes; symmetric=system.symmetric, kwargs...)
    else
        error("At least one of `write_surfaces` or `write_wakes` must be true")
    end

    return nothing
end

"""
    write_vtk(name, surfaces, [surface_properties]; kwargs...)
    write_vtk(name, wakes; kwargs...)
    write_vtk(name, surfaces, wakes, [surface_properties]; kwargs...)

Write geometry from surfaces and/or wakes to Paraview files for visualization.

# Arguments
 - `name`: Base name for the generated files
 - `surfaces`:
   - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
   or
   - Vector of matrices of shape (nc, ns) containing surface panels (see
    [`SurfacePanel`](@ref))
   where `nc` is the number of chordwise panels and `ns` is the number of
   spanwise panels
 - `wakes`: (optional) Vector of wakes corresponding to each surface, represented
    by matrices of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where
    `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels.
 - `surface_properties`: (optional) Vector of surface panel properties for each
    surface, stored as matrices of panel properties (see [`PanelProperties`](@ref))
    of shape (nc, ns) where `nc` is the number of chordwise panels and `ns` is
    the number of spanwise panels

# Keyword Arguments:
 - `symmetric`: (required if `surface_properties` is provided) Flags indicating whether a
    mirror image (across the X-Z plane) was used when calculating induced velocities
    for each surface.
 - `trailing_vortices`: Flag indicating whether the model uses trailing vortices.
    Defaults to `true` when wake panels are absent, `false` otherwise
 - `xhat`: Direction in which trailing vortices extend if used. Defaults to [1, 0, 0].
 - `wake_length`: Distance to extend trailing vortices. Defaults to 10
 - `metadata`: Dictionary of metadata to include in generated files
"""
write_vtk(name, surfaces::AbstractVector{<:AbstractMatrix}, args...; kwargs...)

function write_vtk(name, surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}},
    properties=fill(nothing, length(surfaces)); symmetric=fill(nothing, length(surfaces)), 
    trailing_edge_list=fill(true, length(surfaces)), kwargs...)

    # create paraview multiblock file
    vtk_multiblock(name) do vtmfile
        # loop through all surfaces
        for i = 1:length(surfaces)
            # add paraview files corresponding to the surface to the multiblock file
            write_vtk!(vtmfile, surfaces[i], properties[i]; trailing_edge=trailing_edge_list[i], symmetric=symmetric[i], kwargs...)
        end
    end

    return nothing
end

function write_vtk(name, wakes::AbstractVector{<:AbstractMatrix{<:WakePanel}}; kwargs...)

    # create paraview multiblock file
    vtk_multiblock(name) do vtmfile
        # loop through all wakes
        for i = 1:length(wakes)
            # add paraview files corresponding to the wake to the multiblock file
            write_vtk!(vtmfile, wakes[i]; kwargs..., surface_index = i)
        end
    end

    return nothing
end

function write_vtk(name, surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}},
    wakes::AbstractVector{<:AbstractMatrix{<:WakePanel}}, properties=nothing; 
        trailing_edge_list=fill(true, length(surfaces)), kwargs...)

    # create multiblock file
    vtk_multiblock(name) do vtmfile

        # loop through all surfaces
        for i = 1:length(surfaces)

            # extract circulation at the trailing edge of `surface`
            if isnothing(properties)
                surface_circulation = zeros(size(surfaces[i], 2))
            else
                surface_circulation = getproperty.(properties[i][end,:], :gamma)
            end

            # extract circulation at the leading edge of `wake`
            wake_circulation = getproperty.(wakes[i][1,:], :gamma)

            # add paraview files corresponding to the surface to the multiblock file
            write_vtk!(vtmfile, surfaces[i], properties[i]; wake_circulation,
                trailing_edge = isempty(wakes[i]) && trailing_edge_list[i], kwargs..., trailing_vortices = false)

            # add paraview files corresponding to the wake to the multiblock file
            write_vtk!(vtmfile, wakes[i]; surface_circulation, kwargs...)
        end
    end

    return nothing
end

"""
    write_vtk(name, surface_history, property_history, wake_history; kwargs...)

Writes unsteady simulation geometry to Paraview files for visualization.

# Arguments
 - `name`: Base name for the generated files
 - `surface_history`: Vector of surfaces at each time step, where each surface is
    represented by a matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `property_history`: Vector of surface properties for each surface at each
    time step, where surface properties are represented by a matrix of panel
    properties (see [`PanelProperties`](@ref)) of shape (nc, ns) where `nc` is
    the number of chordwise panels and `ns` is the number of spanwise panels
 - `wake_history`: Vector of wakes corresponding to each surface at each time step,
    where each wake is represented by a matrix of wake panels (see [`WakePanel`](@ref))
    of shape (nw, ns) where `nw` is the number of chordwise wake panels and
    `ns` is the number of spanwise panels.
 - `dt`: Time step vector

# Keyword Arguments:
 - `symmetric`: (required if `properties` is provided) Flags indicating whether a
    mirror image (across the X-Z plane) was used when calculating induced velocities
    for each surface.
 - `wake_length`: Distance to extend trailing vortices. Defaults to 10
 - `metadata`: Dictionary of metadata to include in generated files
"""
function write_vtk(name, surface_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}},
    property_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}},
    wake_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}}, dt;
    symmetric = fill(nothing, length(surface_history[1])), kwargs...)

    symmetric = isa(symmetric, Number) ? fill(symmetric, length(surface_history[1])) : symmetric

    # create paraview collection file
    paraview_collection(name) do pvdfile

        # construct time vector
        time = cumsum(dt)

        # loop through each time step
        for it = 1:length(time)

            # construct multiblock file for each time step
            vtk_multiblock(name*"-step$it") do vtmfile

                # loop through all surfaces
                for i = 1:length(surface_history[it])

                    # extract circulation at the trailing edge of `surface`
                    surface_circulation = getproperty.(property_history[it][i][end,:], :gamma)

                    # extract circulation at the leading edge of `wake`
                    if isempty(wake_history[it][i])
                        wake_circulation = zeros(size(surface_history[it][i], 2))
                    else
                        wake_circulation = getproperty.(wake_history[it][i][1,:], :gamma)
                    end

                    # add paraview files corresponding to the surface to the multiblock file
                    write_vtk!(vtmfile, surface_history[it][i], property_history[it][i];
                        wake_circulation,
                        trailing_edge = isempty(wake_history[it][i]),
                        symmetric = symmetric[i],
                        kwargs...,
                        trailing_vortices = false)

                    # add paraview files corresponding to the wake to the multiblock file
                    write_vtk!(vtmfile, wake_history[it][i];
                        surface_circulation,
                        symmetric = symmetric[i],
                        kwargs...)
                end

                # add multiblock file to the paraview collection file
                pvdfile[time[it]] = vtmfile
            end
        end
    end

    return nothing
end

"""
    write_vtk!(vtmfile, surface, [surface_properties]; kwargs...)

Writes geometry to Paraview files for visualization.

# Arguments
 - `vtmfile`: Multiblock file handle
 - `surface`: Matrix of surface panels (see [`SurfacePanel`](@ref)) of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `surface_properties`: (optional) Matrix of panel properties for each non-wake panel
    where each element of the matrix is of type [`PanelProperties`](@ref).

# Keyword Arguments:
 - `symmetric`: (required if `properties` is provided) Flag indicating whether a
    mirror image (across the X-Z plane) was used when calculating induced velocities.
 - `trailing_vortices = true`: Flag indicating whether the model uses trailing vortices
 - `xhat = [1, 0, 0]`: Direction in which trailing vortices extend if used
 - `wake_length = 10`: Distance to extend trailing vortices
 - `wake_circulation = zeros(size(surfaces, 2))`: Contribution to the trailing
    edge circulation from the wake attached to this surface
 - `metadata = Dict()`: Dictionary of metadata to include in generated files
"""
function write_vtk!(vtmfile, surface::AbstractMatrix{<:SurfacePanel}, properties=nothing;
    symmetric = nothing,
    trailing_vortices = true,
    trailing_edge = true,
    xhat = SVector(1, 0, 0),
    wake_length = 10,
    wake_circulation = zeros(size(surface, 2)),
    metadata = Dict())

    # get float type
    TF = eltype(eltype(surface))

    # check to make sure `symmetric` is provided if `properties` is provided
    @assert !(!isnothing(properties) && isnothing(symmetric)) "Keyword argument `symmetric` is required when optional argument `properties` is provided"

    # grid dimensions
    nc, ns = size(surface)
    N = length(surface)

    # extract geometry as a grid
    xyz = Array{TF, 4}(undef, 3, nc+1, ns+1, 1)
    for (i, I) in enumerate(CartesianIndices((nc+1, ns+1)))
        if I[1] <= nc && I[2] <= ns
            xyz[:, I, 1] = top_left(surface[I[1], I[2]])
        elseif I[1] <= nc && I[2] == ns + 1
            xyz[:, I, 1] = top_right(surface[I[1], I[2]-1])
        elseif I[1] == nc + 1 && I[2] <= ns
            xyz[:, I, 1] = bottom_left(surface[I[1]-1, I[2]])
        else # I[1] == nc + 1 && I[2] == ns + 1
            xyz[:, I, 1] = bottom_right(surface[I[1]-1, I[2]-1])
        end
    end

    # convert grid to points
    points = reshape(xyz, 3, :)

    # now extract bound vortex geometries
    li = LinearIndices((size(xyz, 2), size(xyz, 3)))
    # horizontal bound vortices
    lines_h = [MeshCell(PolyData.Lines(), [
        li[i,j],
        li[i,j+1]
        ]) for j = 1:ns for i = 1:nc]
    # vertical bound vortices
    lines_v = [MeshCell(PolyData.Lines(), [
        li[i,j],
        li[i+1,j]
        ]) for j = 1:ns+1 for i = 1:nc]

    # now extract data (if applicable)
    if !isnothing(properties)
        # horizontal bound vortex circulation strength
        gamma_h = Matrix{TF}(undef, nc, ns)
        for i = 1:nc, j = 1:ns
            previous_gamma = i == 1 ? 0.0 : properties[i-1, j].gamma
            current_gamma = properties[i,j].gamma
            # check if we need to account for symmetry
            if symmetric && on_symmetry_plane(top_left(surface[i,j]), top_right(surface[i,j]))
                gamma_h[i,j] = 0.0
            else
                gamma_h[i,j] = current_gamma - previous_gamma
            end
        end

        # horizontal bound vortex force coefficient
        cf_h = Array{TF}(undef, 3, nc, ns)
        for i = 1:nc, j = 1:ns

            # check if we need to account for symmetry
            if symmetric && on_symmetry_plane(top_left(surface[i,j]), top_right(surface[i,j]))
                cf_h[:,i,j] .= 0.0
            else
                cf_h[:,i,j] = properties[i,j].cfb
            end
        end

        # horizontal bound vortex local velocity
        v_h = Array{TF}(undef, 3, nc, ns)
        for i = 1:nc, j = 1:ns
            v_h[:,i,j] = properties[i,j].velocity
        end

        # vertical bound vortex circulation strength
        gamma_v = Matrix{TF}(undef, nc, ns+1)
        for i = 1:nc
            current_gamma = properties[i,1].gamma
            # check if we need to account for symmetry
            if symmetric && on_symmetry_plane(bottom_left(surface[i,1]), top_left(surface[i,1]))
                gamma_v[i,1] = 0.0
            else
                gamma_v[i,1] = -current_gamma
            end
            for j = 2:ns
                previous_gamma = current_gamma
                current_gamma = properties[i,j].gamma
                # check if we need to account for symmetry
                if symmetric && on_symmetry_plane(bottom_left(surface[i,j]), top_left(surface[i,j]))
                    gamma_v[i,j] = 0.0
                else
                    gamma_v[i,j] = previous_gamma - current_gamma
                end
            end
            previous_gamma = current_gamma
            # check if we need to account for symmetry
            if symmetric && on_symmetry_plane(bottom_right(surface[i,end]), top_right(surface[i,end]))
                gamma_v[i,end] = 0.0
            else
                gamma_v[i,end] = previous_gamma
            end
        end

        # vertical bound vortex force coefficient
        cf_v = Array{TF}(undef, 3, nc, ns+1)
        for i = 1:nc
            #  check if we need to account for symmetry
            if symmetric && on_symmetry_plane(bottom_left(surface[i,1]), top_left(surface[i,1]))
                cf_v[:,i,1] .= 0.0
            else
                cf_v[:,i,1] = properties[i,1].cfl
            end
            for j = 2:ns
                previous_cf = properties[i,j-1].cfr
                current_cf = properties[i,j].cfl
                #  check if we need to account for symmetry
                if symmetric && on_symmetry_plane(bottom_left(surface[i,j]), top_left(surface[i,j]))
                    cf_v[i,j] = 0.0
                else
                    cf_v[:,i,j] = previous_cf + current_cf
                end
            end
            #  check if we need to account for symmetry
            if symmetric && on_symmetry_plane(bottom_right(surface[i,end]), top_right(surface[i,end]))
                cf_v[:,i,end] .= 0.0
            else
                cf_v[:,i,end] = properties[i,end].cfr
            end
        end
    end

    # horizontal bound vortices
    vtk_grid(vtmfile, points, lines_h) do vtkfile

        # add metadata
        for (key, value) in pairs(metadata)
            vtkfile[string(key)] = value
        end


        if !isnothing(properties)
            # circulation strength
            vtkfile["circulation"] = reshape(gamma_h, :)

            # local velocity
            vtkfile["velocity"] = reshape(v_h, 3, :)

            # force coefficient
            vtkfile["force"] = reshape(cf_h, 3, :)
        end
    end

    # vertical bound vortices
    vtk_grid(vtmfile, points, lines_v) do vtkfile

        # add metadata
        for (key, value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        if !isnothing(properties)
            # circulation strength
            vtkfile["circulation"] = reshape(gamma_v, :)

            # force coefficient
            vtkfile["force"] = reshape(cf_v, 3, :)
        end
    end

    # --- control points ---

    points_cp = Matrix{TF}(undef, 3, N)
    for i = 1:N
        ipoint = i
        points_cp[:,ipoint] = controlpoint(surface[i])
    end

    cells_cp = [MeshCell(PolyData.Verts(), [i]) for i = 1:N]

    vtk_grid(vtmfile, points_cp, cells_cp) do vtkfile

        # add metadata
        for (key, value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        # add normal
        data = Matrix{TF}(undef, 3, N)
        for i = 1:length(surface)
            data[:,i] = normal(surface[i])
        end
        vtkfile["normal"] = data

    end

    # trailing edge and/or trailing vortices
    if trailing_vortices

        # trailing vortices points as a grid
        xyz_t = Array{TF}(undef, 3, 2, ns+1)
        for j = 1:ns
            xyz_t[:,1,j] = bottom_left(surface[end,j])
            xyz_t[:,2,j] = xyz_t[:,1,j] + wake_length*xhat
        end
        xyz_t[:,1,end] = bottom_right(surface[end,end])
        xyz_t[:,2,end] = xyz_t[:,1,end] + wake_length*xhat

        # generate points
        points_t = reshape(xyz_t, 3, :)

        li = LinearIndices((2, ns+1))

        # horizontal bound vortices
        lines_h = [MeshCell(PolyData.Lines(), [li[1,j], li[1,j+1]]) for j = 1:ns]

        # vertical bound vortices
        lines_v = [MeshCell(PolyData.Lines(), [li[1,j], li[2,j]]) for j = 1:ns+1]

        # combine
        lines_t = vcat(lines_h, lines_v)

        # now extract data (if applicable)
        if !isnothing(properties)
            # horizontal bound vortex circulation strength
            gamma_h = zeros(ns)

            # vertical bound vortex circulation strength
            gamma_v = Vector{TF}(undef, ns+1)
            current_gamma = properties[end,1].gamma

            # check if we need to account for symmetry
            if symmetric && on_symmetry_plane(bottom_left(surface[end,1]))
                gamma_v[1] = 0.0
            else
                gamma_v[1] = -current_gamma
            end

            for j = 2:ns
                previous_gamma = current_gamma
                current_gamma = properties[end,j].gamma

                # check if we need to account for symmetry
                if symmetric && on_symmetry_plane(bottom_left(surface[end,j]))
                    gamma_v[j] = 0.0
                else
                    gamma_v[j] = previous_gamma - current_gamma
                end

            end
            previous_gamma = current_gamma

            # check if we need to account for symmetry
            if symmetric && on_symmetry_plane(bottom_right(surface[end,end]))
                gamma_v[end] = 0.0
            else
                gamma_v[end] = previous_gamma
            end

            # combine
            gamma_t = vcat(gamma_h, gamma_v)
        end

    else
        # only generate if we plan to include the trailing edge
        if trailing_edge

            # generate points
            points_t = Matrix{TF}(undef, 3, ns+1)
            for j = 1:ns
                points_t[:,j] = bottom_left(surface[end,j])
            end
            points_t[:,end] = bottom_right(surface[end,end])

            # horizontal bound vortices
            lines_t = [MeshCell(PolyData.Lines(), j:j+1) for j = 1:ns]

            # now extract data (if applicable)
            if !isnothing(properties)
                # horizontal bound vortex circulation strength
                gamma_t = Vector{TF}(undef, ns)
                for j = 1:ns
                    # check if we need to account for symmetry
                    if symmetric && on_symmetry_plane(bottom_left(surface[end,j]), bottom_right(surface[end,j]))
                        gamma_t[j] = 0.0
                    else
                        gamma_t[j] = wake_circulation[j] - properties[end,j].gamma
                    end
                end
            end
        end

    end

    if trailing_vortices || trailing_edge

        # trailing edge and trailing vortices
        vtk_grid(vtmfile, points_t, lines_t) do vtkfile

            # add metadata
            for (key, value) in pairs(metadata)
                vtkfile[string(key)] = value
            end

            if !isnothing(properties)
                # circulation strength
                vtkfile["circulation"] = gamma_t
            end
        end
    end

    return nothing
end

"""
    write_vtk!(vtmfile, wake; kwargs...)

Writes geometry to Paraview files for visualization.

# Arguments
 - `vtmfile`: Paraview file handle
 - `wake`: Matrix of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns)
    where `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels

# Keyword Arguments:
 - `symmetric`: (required) Flag indicating whether a mirror image (across the
    X-Z plane) was used when calculating induced velocities.
 - `trailing_vortices = false`: Flag indicating whether the model uses trailing vortices
 - `xhat = [1, 0, 0]`: Direction in which trailing vortices extend if used
 - `wake_length = 10`: Distance to extend trailing vortices
 - `surface_circulation = zeros(size(wake, 2))`: Contribution to the leading edge
    circulation from the surface attached to this wake.
 - `metadata = Dict()`: Dictionary of metadata to include in generated files
"""
function write_vtk!(vtmfile, wake::AbstractMatrix{<:WakePanel};
    symmetric,
    trailing_vortices = false,
    xhat = SVector(1, 0, 0),
    wake_length = 10,
    surface_circulation = zeros(size(wake, 2)),
    metadata=Dict())

    # do nothing if no wake panels are present
    if isempty(wake)
        return vtmfile
    end

    # extract float type
    TF = eltype(eltype(wake))

    # get wake dimensions
    nc, ns = size(wake)
    N = length(wake)

    # extract geometry as a grid
    xyz = Array{TF, 4}(undef, 3, nc+1, ns+1, 1)
    for (i, I) in enumerate(CartesianIndices((nc+1, ns+1)))
        if I[1] <= nc && I[2] <= ns
            xyz[:, I, 1] = top_left(wake[I[1], I[2]])
        elseif I[1] <= nc && I[2] == ns + 1
            xyz[:, I, 1] = top_right(wake[I[1], I[2]-1])
        elseif I[1] == nc + 1 && I[2] <= ns
            xyz[:, I, 1] = bottom_left(wake[I[1]-1, I[2]])
        else # I[1] == nc + 1 && I[2] == ns + 1
            xyz[:, I, 1] = bottom_right(wake[I[1]-1, I[2]-1])
        end
    end

    # convert to points
    points = reshape(xyz, 3, :)

    # now extract bound vortex geometries
    li = LinearIndices(size(xyz)[2:3])
    # horizontal bound vortices
    lines_h = [MeshCell(PolyData.Lines(), [
        li[i,j],
        li[i,j+1]
        ]) for j = 1:ns for i = 1:nc]
    # vertical bound vortices
    lines_v = [MeshCell(PolyData.Lines(), [
        li[i,j],
        li[i+1,j]
        ]) for j = 1:ns+1 for i = 1:nc]

    # horizontal bound vortex circulation strength
    gamma_h = Matrix{TF}(undef, nc, ns)
    for i in 1:nc, j = 1:ns
        previous_gamma = i == 1 ? surface_circulation[j] : wake[i-1, j].gamma
        current_gamma = wake[i,j].gamma
        if symmetric && on_symmetry_plane(top_left(wake[nc, ns]), top_right(wake[nc, ns]))
            gamma_h[i,j] = 0.0
        else
            gamma_h[i,j] = current_gamma - previous_gamma
        end
    end

    # vertical bound vortex circulation strength
    gamma_v = Matrix{TF}(undef, nc, ns+1)
    for i = 1:nc
        current_gamma = wake[i,1].gamma
        if symmetric && on_symmetry_plane(bottom_left(wake[i,1]), top_left(wake[i,1]))
            gamma_v[i,1] = 0.0
        else
            gamma_v[i,1] = -current_gamma
        end
        for j = 2:ns
            previous_gamma = current_gamma
            current_gamma = wake[i,j].gamma
            if symmetric && on_symmetry_plane(bottom_left(wake[i,j]), top_left(wake[i,j]))
                gamma_v[i,j] = 0.0
            else
                gamma_v[i,j] = previous_gamma - current_gamma
            end
        end
        previous_gamma = current_gamma
        if symmetric && on_symmetry_plane(bottom_right(wake[i,end]), top_right(wake[i,end]))
            gamma_v[i,end] = 0.0
        else
            gamma_v[i,end] = previous_gamma
        end
    end

    # horizontal bound vortices
    vtk_grid(vtmfile, points, lines_h) do vtkfile

        # add metadata
        for (key, value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        # circulation strength
        vtkfile["circulation"] = reshape(gamma_h, :)
    end

    # vertical bound vortices
    vtk_grid(vtmfile, points, lines_v) do vtkfile

        # add metadata
        for (key, value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        # circulation strength
        vtkfile["circulation"] = reshape(gamma_v, :)
    end

    if trailing_vortices

        # trailing vortices points as a grid
        xyz_t = Array{TF}(undef, 3, 2, ns+1)
        for j = 1:ns
            xyz_t[:,1,j] = bottom_left(wake[end,j])
            xyz_t[:,2,j] = xyz_t[:,1,j] + wake_length*xhat
        end
        xyz_t[:,1,end] = bottom_right(wake[end,end])
        xyz_t[:,2,end] = xyz_t[:,1,end] + wake_length*xhat

        # generate points
        points_t = reshape(xyz_t, 3, :)

        li = LinearIndices((2, ns+1))

        # horizontal bound vortices
        lines_h = [MeshCell(PolyData.Lines(), [li[1,j], li[1,j+1]]) for j = 1:ns]

        # vertical bound vortices
        lines_v = [MeshCell(PolyData.Lines(), [li[1,j], li[2,j]]) for j = 1:ns+1]

        # combine
        lines_t = vcat(lines_h, lines_v)

        # horizontal bound vortex circulation strength
        gamma_h = zeros(ns)

        # vertical bound vortex circulation strength
        gamma_v = Vector{TF}(undef, ns+1)
        current_gamma = wake[end,1].gamma
        if symmetric && on_symmetry_plane(bottom_left(wake[end,1]))
            gamma_v[1] = 0.0
        else
            gamma_v[1] = -current_gamma
        end
        for j = 2:ns
            previous_gamma = current_gamma
            current_gamma = wake[end,j].gamma
            if symmetric && on_symmetry_plane(bottom_left(wake[end,j]))
                gamma_v[j] = 0.0
            else
                gamma_v[j] = previous_gamma - current_gamma
            end
        end
        previous_gamma = current_gamma
        if symmetric && on_symmetry_plane(bottom_right(wake[end,end]))
            gamma_v[j] = 0.0
        else
            gamma_v[end] = previous_gamma
        end

        # combine
        gamma_t = vcat(gamma_h, gamma_v)

    else

        # generate points
        points_t = Matrix{TF}(undef, 3, ns+1)
        for j = 1:ns
            points_t[:,j] = bottom_left(wake[end,j])
        end
        points_t[:,end] = bottom_right(wake[end,end])

        # horizontal bound vortices
        lines_t = [MeshCell(PolyData.Lines(), j:j+1) for j = 1:ns]

        # horizontal bound vortex circulation strength
        gamma_t = Vector{TF}(undef, ns)
        for j = 1:ns
            if symmetric && on_symmetry_plane(bottom_left(wake[end,j]), bottom_right(wake[end,j]))
                gamma_t[j] = 0.0
            else
                gamma_t[j] = -wake[end,j].gamma
            end
        end

    end

    # trailing edge and trailing vortices
    vtk_grid(vtmfile, points_t, lines_t) do vtkfile

        # add metadata
        for (key, value) in pairs(metadata)
            vtkfile[string(key)] = value
        end

        # circulation strength
        vtkfile["circulation"] = reshape(gamma_t, :)
    end

    return nothing
end

#--- reference frames ---#

function write_vtk(name::String, frames::Vector{<:ReferenceFrame}; kwargs...)

    # create paraview file
    vtk(name) do vtkfile

        # loop through all frames
        for frame in frames

            # get transformation to global frame


            # add paraview file corresponding to the frame
            write_vtk!(vtkfile, frame; kwargs..., frame_index = i)
        end
    end

    return nothing
end

# #--- vortex filaments ---#

# function write_vtk(fname, vortex_filaments::VortexFilaments)
#     # create points
#     pts = zeros(SVector{3,eltype(vortex_filaments)}, length(vortex_filaments.filaments) * 2)
#     ic = 1
#     for i in 1:length(vortex_filaments.filaments)
#         pts[ic] = vortex_filaments.filaments[i].r1
#         ic += 1
#         pts[ic] = vortex_filaments.filaments[i].r2
#         ic += 1
#     end

#     # create lines
#     lines = [MeshCell(PolyData.Lines(), (2*i-1, 2*i)) for i in 1:length(vortex_filaments.filaments)]

#     # save strengths
#     strengths = Vector{SVector{3,eltype(vortex_filaments)}}(undef, length(vortex_filaments.filaments))
#     for i in 1:length(vortex_filaments.filaments)
#         filament = vortex_filaments.filaments[i]
#         s = filament.r2 - filament.r1
#         if norm(s) == 0.0
#             strengths[i] = SVector{3,eltype(vortex_filaments)}(0.0, 0, 0)
#         else
#             strengths[i] = filament.strength * s / norm(s)
#         end
#     end

#     # save as VTK
#     vtk_grid(fname, pts, lines) do vtk
#         vtk["strength"] = strengths
#         vtk["velocity"] = vortex_filaments.velocity
#     end
# end