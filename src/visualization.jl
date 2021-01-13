"""
    write_vtk(name, surface, [properties]; kwargs...)
    write_vtk(name, wake; kwargs...)
    write_vtk(name, surface, wake, [properties]; kwargs...)

Writes geometry to Paraview files for visualization.

# Arguments
 - `name`: Name for the generated files
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `wake`: (optional) Matrix of wake panels of shape (nw, ns) where `nw` is the
    number of chordwise wake panels and `ns` is the number of spanwise panels,
    defaults to no wake panels
 - `properties`: (optional) Matrix of panel properties for each non-wake panel
    where each element of the matrix is of type `PanelProperties`.

# Keyword Arguments:
 - `mirror`: Flag indicating whether to include a mirror image of the panels in
    `surface`. Defaults to false.
 - `trailing_vortices`: Flag indicating whether the model uses trailing vortices.
    Defaults to `true` when wake panels are absent, `false` otherwise
 - `xhat`: Direction in which trailing vortices extend if used. Defaults to [1, 0, 0].
 - `wake_length`: Distance to extend trailing vortices. Defaults to 10
 - `metadata`: Dictionary of metadata to include in generated files
 """
write_vtk(name, surface::AbstractMatrix, args...; kwargs...)

# surface implementation
function write_vtk(name, surface::AbstractMatrix{<:SurfacePanel}, properties=nothing; kwargs...)

    vtk_multiblock(name) do vtmfile
        write_vtk!(vtmfile, surface, properties; kwargs...)
    end

    return nothing
end

# wake implementation
function write_vtk(name, surface::AbstractMatrix{<:WakePanel}; kwargs...)

    vtk_multiblock(name) do vtmfile
        write_vtk!(vtmfile, surface; kwargs...)
    end

    return nothing
end

# surface + wake implementation
function write_vtk(name, surface::AbstractMatrix{<:SurfacePanel},
    wake::AbstractMatrix{<:WakePanel}, properties=nothing; kwargs...)

    vtk_multiblock(name) do vtmfile

        # circulation at trailing edge of surface
        if isnothing(properties)
            surface_circulation = zeros(size(surface, 2))
        else
            surface_circulation = getproperty.(properties[end,:], :gamma)
        end

        # circulation at leading edge of wake
        wake_circulation = getproperty.(wake[1,:], :gamma)

        # surface
        write_vtk!(vtmfile, surface, properties; wake_circulation, trailing_edge = false,
            kwargs..., trailing_vortices = false)

        # wake
        write_vtk!(vtmfile, wake; surface_circulation, kwargs...)
    end

    return nothing
end

"""
    write_vtk(name, surface, panel_history, wake_history; kwargs...)

Writes unsteady simulation geometry to Paraview files for visualization.

# Arguments
 - `name`: Name for generated files
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `panel_history`: Matrices corresponding to the panel properties for each
    non-wake panel at each time step.
 - `wake_history`: Matrices corresponding to the wake panels at each time step

# Keyword Arguments:
 - `mirror`: Flag indicating whether to include a mirror image of the panels in
    `surface`. Defaults to false.
 - `trailing_vortices`: Flag indicating whether the model uses trailing vortices.
    Defaults to `true` when wake panels are absent, `false` otherwise
 - `xhat`: Direction in which trailing vortices extend if used. Defaults to [1, 0, 0].
 - `wake_length`: Distance to extend trailing vortices. Defaults to 10
 - `metadata`: Dictionary of metadata to include in generated files
"""
function write_vtk(name, surface::AbstractMatrix,
    panel_history::AbstractVector{<:AbstractMatrix},
    wake_history::AbstractVector{<:AbstractMatrix}, time; kwargs...)

    paraview_collection(name) do pvdfile
        for it = 1:length(time)
            vtk_multiblock(name*"-step$it") do vtmfile

                surface_circulation = getproperty.(panel_history[it][end,:], :gamma)

                if isempty(wake_history[it])
                    wake_circulation = zeros(size(surface, 2))
                else
                    wake_circulation = getproperty.(wake_history[it][1,:], :gamma)
                end

                # surface
                write_vtk!(vtmfile, surface, panel_history[it]; wake_circulation,
                    trailing_edge = isempty(wake_history[it]), kwargs...,
                    trailing_vortices = false)

                # wake
                write_vtk!(vtmfile, wake_history[it]; surface_circulation,
                    kwargs...)

                # save
                pvdfile[time[it]] = vtmfile
            end
        end
    end

    return nothing
end

"""
    write_vtk(name, surfaces, [properties]; kwargs...)
    write_vtk(name, wakes; kwargs...)
    write_vtk(name, surfaces, wakes, [properties]; kwargs...)

Writes geometry to Paraview files for visualization.

# Arguments
 - `name`: Name for generated files
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `wakes`: (optional) Vector of wakes corresponding to each surface, represented
    by matrices of panels of shape (nw, ns) where `nw` is the number of chordwise
    wake panels and `ns` is the number of spanwise panels
 - `properties`: (optional) Vector of matrices corresponding to the panel
    properties for each non-wake panel.

# Keyword Arguments:
 - `mirror`: Flag indicating whether to include a mirror image of the panels in
    `surface`. Defaults to false.
 - `trailing_vortices`: Flag indicating whether the model uses trailing vortices.
    Defaults to `true` when wake panels are absent, `false` otherwise
 - `xhat`: Direction in which trailing vortices extend if used. Defaults to [1, 0, 0].
 - `wake_length`: Distance to extend trailing vortices. Defaults to 10
 - `metadata`: Dictionary of metadata to include in generated files
"""
write_vtk(name, surfaces::AbstractVector{<:AbstractMatrix}, args...; kwargs...)

function write_vtk(name, surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}},
    properties=nothing; kwargs...) where T

    vtk_multiblock(name) do vtmfile
        for i = 1:length(surfaces)
            write_vtk!(vtmfile, surfaces[i], properties[i]; kwargs...)
        end
    end

    return nothing
end

function write_vtk(name, surfaces::AbstractVector{<:AbstractMatrix{<:WakePanel}}; kwargs...) where T

    vtk_multiblock(name) do vtmfile
        for i = 1:length(surfaces)
            write_vtk!(vtmfile, surfaces[i]; kwargs..., surface_index = i)
        end
    end

    return nothing
end

function write_vtk(name, surfaces::AbstractVector{<:AbstractMatrix},
    wakes::AbstractVector{<:AbstractMatrix}, properties=nothing; kwargs...)

    vtk_multiblock(name) do vtmfile
        for i = 1:length(surfaces)

            # circulation at trailing edge of surface
            if isnothing(properties)
                surface_circulation = zeros(size(surfaces[i], 2))
            else
                surface_circulation = getproperty.(properties[i][end,:], :gamma)
            end

            # circulation at leading edge of wake
            wake_circulation = getproperty.(wakes[i][1,:], :gamma)

            # surface
            write_vtk!(vtmfile, surfaces[i], properties[i]; wake_circulation,
                trailing_edge = isempty(wakes[i]), kwargs..., trailing_vortices = false)

            # wake
            write_vtk!(vtmfile, wakes[i]; surface_circulation, kwargs...)
        end
    end

    return nothing
end

"""
    write_vtk(name, surfaces, panel_history, wake_history; kwargs...)

Writes unsteady simulation geometry to Paraview files for visualization.

# Arguments
 - `name`: Name for generated files
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `panel_history`: Matrices corresponding to the panel properties for each
    non-wake panel for each surface at each time step.
 - `wake_history`: Matrices corresponding to the wake panels for each surface for
    each time step

# Keyword Arguments:
 - `mirror`: Flag indicating whether to include a mirror image of the panels in
    `surface`. Defaults to false.
 - `trailing_vortices`: Flag indicating whether the model uses trailing vortices.
    Defaults to `true` when wake panels are absent, `false` otherwise
 - `xhat`: Direction in which trailing vortices extend if used. Defaults to [1, 0, 0].
 - `wake_length`: Distance to extend trailing vortices. Defaults to 10
 - `metadata`: Dictionary of metadata to include in generated files
"""
function write_vtk(name, surfaces::AbstractVector{<:AbstractMatrix},
    panel_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}},
    wake_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}}, time; kwargs...)

    paraview_collection(name) do pvdfile
        for it = 1:length(time)
            vtk_multiblock(name*"-step$it") do vtmfile
                # generate VTK files
                for i = 1:length(surfaces)

                    surface_circulation = getproperty.(panel_history[it][i][end,:], :gamma)

                    if isempty(wake_history[it][i])
                        wake_circulation = zeros(size(surfaces[i], 2))
                    else
                        wake_circulation = getproperty.(wake_history[it][i][1,:], :gamma)
                    end

                    # surface
                    write_vtk!(vtmfile, surfaces[i], panel_history[it][i];
                        wake_circulation, trailing_edge = isempty(wake_history[it][i]),
                        kwargs..., trailing_vortices = false)

                    # wake
                    write_vtk!(vtmfile, wake_history[it][i]; surface_circulation,
                        kwargs...)
                end

                # save
                pvdfile[time[it]] = vtmfile
            end
        end
    end

    return nothing
end

"""
    write_vtk!(vtmfile, surface, [properties]; kwargs...)

Writes geometry to Paraview files for visualization.

# Arguments
 - `vtmfile`: Paraview file handle
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `properties`: (optional) Matrix of panel properties for each non-wake panel
    where each element of the matrix is of type `PanelProperties`.

# Keyword Arguments:
 - `mirror = false`: Flag indicating whether to include a mirror image of the
    panels in `surface`
 - `trailing_vortices = true`: Flag indicating whether the model uses trailing vortices
 - `xhat = [1, 0, 0]`: Direction in which trailing vortices extend if used
 - `wake_length = 10`: Distance to extend trailing vortices
 - `wake_circulation = zeros(size(surfaces, 2))`: Contribution to the trailing
    edge circulation from the wake attached to this surface
 - `metadata = Dict()`: Dictionary of metadata to include in generated files
"""
function write_vtk!(vtmfile,
    surface::AbstractMatrix{<:SurfacePanel}, properties=nothing;
    mirror = false, trailing_vortices = true, trailing_edge = true, xhat = SVector(1, 0, 0),
    wake_length = 10, wake_circulation = zeros(size(surface, 2)), metadata=Dict())

    TF = eltype(eltype(surface))

    if mirror
        surface = vcat(reflect(surface), surface)
    end

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

    # now extract data (if applicable)
    if !isnothing(properties)
        # horizontal bound vortex circulation strength
        gamma_h = Matrix{TF}(undef, nc, ns)
        for i = 1:nc, j = 1:ns
            previous_gamma = i == 1 ? 0.0 : properties[i-1, j].gamma
            current_gamma = properties[i,j].gamma
            gamma_h[i,j] = current_gamma - previous_gamma
        end

        # horizontal bound vortex local velocity
        v_h = Array{TF}(undef, 3, nc, ns)
        for i = 1:nc, j = 1:ns
            v_h[:,i,j] = properties[i,j].velocity
        end

        # horizontal bound vortex force coefficient
        cf_h = Array{TF}(undef, 3, nc, ns)
        for i = 1:nc, j = 1:ns
            cf_h[:,i,j] = properties[i,j].cfb
        end

        # vertical bound vortex circulation strength
        gamma_v = Matrix{TF}(undef, nc, ns+1)
        for i = 1:nc
            current_gamma = properties[i,1].gamma
            gamma_v[i,1] = -current_gamma
            for j = 2:ns
                previous_gamma = current_gamma
                current_gamma = properties[i,j].gamma
                gamma_v[i,j] = previous_gamma - current_gamma
            end
            previous_gamma = current_gamma
            gamma_v[i,end] = previous_gamma
        end

        # vertical bound vortex force coefficient
        cf_v = Array{TF}(undef, 3, nc, ns+1)
        for i = 1:nc
            cf_v[:,i,1] = properties[i,1].cfl
            for j = 2:ns
                previous_cf = properties[i,j-1].cfr
                current_cf = properties[i,j].cfl
                cf_v[:,i,j] = previous_cf + current_cf
            end
            cf_v[:,i,end] = properties[i,end].cfr
        end

        # mirror across X-Z plane
        if mirror
            gamma_h = cat(reverse(gamma_h, dims=2), gamma_v, dims=2)
            v_h = cat(reverse(v_h, dims=3), v_h, dims=3)
            cf_h = cat(reverse(cf_h, dims=3), cf_h, dims=3)
            gamma_v = cat(reverse(gamma_v, dims=2), gamma_v, dims=2)
            cf_v = cat(reverse(cf_v, dims=3), cf_v, dims=3)
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
            gamma_v[1] = -current_gamma
            for j = 2:ns
                previous_gamma = current_gamma
                current_gamma = properties[end,j].gamma
                gamma_v[j] = previous_gamma - current_gamma
            end
            previous_gamma = current_gamma
            gamma_v[end] = previous_gamma

            # mirror across X-Z plane
            if mirror
                gamma_h = vcat(reverse(gamma_h), gamma_h)
                gamma_v = vcat(reverse(gamma_v), gamma_v)
            end

            # combine
            gamma_t = vcat(gamma_h, gamma_v)
        end

    else

        # only generate if we plane to include the trailing edge
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
                    gamma_t[j] = wake_circulation[j] - properties[end,j].gamma
                end

                # mirror across X-Z plane
                if mirror
                    gamma_t = vcat(reverse(gamma_t), gamma_t)
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
 - `wake`: Matrix of wake panels of shape (nw, ns) where `nw` is the number of
    chordwise wake panels and `ns` is the number of spanwise panels,
    defaults to no wake panels

# Keyword Arguments:
 - `mirror = false`: Flag indicating whether to include a mirror image of the
    panels in `surface`
 - `trailing_vortices = false`: Flag indicating whether the model uses trailing vortices
 - `xhat = [1, 0, 0]`: Direction in which trailing vortices extend if used
 - `wake_length = 10`: Distance to extend trailing vortices
 - `surface_circulation = zeros(size(wake, 2))`: Contribution to the leading edge
    circulation from the surface attached to this wake.
 - `metadata = Dict()`: Dictionary of metadata to include in generated files
"""
function write_vtk!(vtmfile, wake::AbstractMatrix{<:WakePanel}; mirror = false,
    trailing_vortices = false, xhat = SVector(1, 0, 0), wake_length = 10,
    surface_circulation = zeros(size(wake, 2)), metadata=Dict())

    if isempty(wake)
        return vtmfile
    end

    TF = eltype(eltype(wake))

    if mirror
        wake = vcat(reflect(wake), wake)
    end

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
        gamma_h[i,j] = current_gamma - previous_gamma
    end

    # vertical bound vortex circulation strength
    gamma_v = Matrix{TF}(undef, nc, ns+1)
    for i = 1:nc
        current_gamma = wake[i,1].gamma
        gamma_v[i,1] = -current_gamma
        for j = 2:ns
            previous_gamma = current_gamma
            current_gamma = wake[i,j].gamma
            gamma_v[i,j] = previous_gamma - current_gamma
        end
        previous_gamma = current_gamma
        gamma_v[i,end] = previous_gamma
    end

    # mirror across X-Z plane
    if mirror
        gamma_h = cat(reverse(gamma_h, dims=2), gamma_h, dims=2)
        gamma_v = cat(reverse(gamma_v, dims=2), gamma_v, dims=2)
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
        gamma_v[1] = -current_gamma
        for j = 2:ns
            previous_gamma = current_gamma
            current_gamma = wake[end,j].gamma
            gamma_v[j] = previous_gamma - current_gamma
        end
        previous_gamma = current_gamma
        gamma_v[end] = previous_gamma

        # mirror across X-Z plane
        if mirror
            gamma_h = vcat(reverse(gamma_h), gamma_h)
            gamma_v = vcat(reverse(gamma_v), gamma_v)
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
            gamma_t[j] = -wake[end,j].gamma
        end

        # mirror across X-Z plane
        if mirror
            gamma_t = vcat(reverse(gamma_t), gamma_t)
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
