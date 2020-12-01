"""
    write_vtk(name, surface, [system]; kwargs...)

Writes geometry to Paraview files for visualization.

# Keyword Arguments:
 - `xhat=[1, 0, 0]`: direction in which trailing vortices are shed
 - `mirror=false`: creates a mirror image of the geometry across the X-Z axis
 - `wake_length=10`: distance to extend the trailing vortices
 - `metadata=Dict()`: dictionary of metadata to include in generated files
"""
function write_vtk(name, surface, system=nothing; xhat=SVector(1, 0, 0),
    wake_length=10, mirror=false, surface_index = 1, metadata=Dict())

    TF = eltype(eltype(surface))

    if mirror
        surface = vcat(reflect(surface), surface)
    end

    nc, ns = size(surface)
    N = length(surface)

    vtk_multiblock(name) do vtmfile

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
        # top bound vortices
        lines_t = [MeshCell(PolyData.Lines(), [
            li[i,j],
            li[i,j+1]
            ]) for i = 1:nc for j = 1:ns]
        # left bound vortices
        lines_l = [MeshCell(PolyData.Lines(), [
            li[i,j],
            ifelse(on_trailing_edge(surface[i,j], false),
                li[end,j], li[i+1,j])
            ]) for i = 1:nc for j = 1:ns]
        # right bound vortices
        lines_r = [MeshCell(PolyData.Lines(), [
            ifelse(on_trailing_edge(surface[i,j], false),
                li[end,j+1], li[i+1,j+1]),
            li[i,j+1]
            ]) for i = 1:nc for j = 1:ns]

        # now extract data (if applicable)
        if !isnothing(system)
            # circulation strength
            gamma = Vector{TF}(undef, N)
            for (i, I) in enumerate(CartesianIndices(system.panels[surface_index]))
                gamma[i] = system.panels[surface_index][i].gamma
            end

            # local velocity
            v = Matrix{TF}(undef, 3, N)
            for (i, I) in enumerate(CartesianIndices(system.panels[surface_index]))
                v[:,i] = system.panels[surface_index][i].v
            end

            # force coefficient
            cf = Matrix{TF}(undef, 3, N)
            for (i, I) in enumerate(CartesianIndices(system.panels[surface_index]))
                cf[:,i] = system.panels[surface_index][i].cf
            end

            # force coefficient
            cfl = Matrix{TF}(undef, 3, N)
            for (i, I) in enumerate(CartesianIndices(system.panels[surface_index]))
                cfl[:,i] = system.panels[surface_index][i].cfl
            end

            # force coefficient
            cfr = Matrix{TF}(undef, 3, N)
            for (i, I) in enumerate(CartesianIndices(system.panels[surface_index]))
                cfr[:,i] = system.panels[surface_index][i].cfr
            end
        end

        # top bound vortices
        vtk_grid(vtmfile, points, lines_t) do vtkfile

            # add metadata
            for (key, value) in pairs(metadata)
                vtkfile[string(key)] = value
            end


            if !isnothing(system)
                # circulation strength
                vtkfile["circulation"] = gamma

                # local velocity
                vtkfile["velocity"] = v

                # force coefficient
                vtkfile["force"] = cf
            end
        end

        # left bound vortices
        vtk_grid(vtmfile, points, lines_l) do vtkfile

            # add metadata
            for (key, value) in pairs(metadata)
                vtkfile[string(key)] = value
            end

            # TODO: add circulation strength

            # force coefficient
            vtkfile["force"] = cfl
        end

        # left bound vortices
        vtk_grid(vtmfile, points, lines_r) do vtkfile

            # add metadata
            for (key, value) in pairs(metadata)
                vtkfile[string(key)] = value
            end

            # TODO: add circulation strength

            # force coefficient
            vtkfile["force"] = cfr
        end

        # trailing vortices
        points_t = Matrix{TF}(undef, 3, 4*ns)
        ipoint = 0
        for is = 1:ns
            points_t[:, ipoint+1] = bottom_left(surface[is]) + wake_length*xhat
            points_t[:, ipoint+2] = bottom_left(surface[is])
            points_t[:, ipoint+3] = bottom_right(surface[is])
            points_t[:, ipoint+4] = bottom_right(surface[is]) + wake_length*xhat
            ipoint += 4
        end

        cells_t = [MeshCell(PolyData.Lines(), 4*(is-1)+1 : 4*(is-1)+4) for is = 1:ns]

        vtk_grid(vtmfile, points_t, cells_t) do vtkfile

            # add metadata
            for (key, value) in pairs(metadata)
                vtkfile[string(key)] = value
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
    end

    return nothing
end

function write_vtk(name, panels::AbstractMatrix{<:Wake}; xhat=SVector(1, 0, 0), wake_length=0, mirror=false, metadata=Dict())

    TF = eltype(eltype(panels))

    if mirror
        panels = vcat(reflect(panels), panels)
    end

    npanels = length(panels)

    vtk_multiblock(name) do vtmfile

        # --- vortex rings ---

        points_r = Matrix{TF}(undef, 3, 4*npanels)
        for i = 1:length(panels)
            ipoint = 4*(i - 1)
            points_r[:,ipoint+1] = panels[i].rbl
            points_r[:,ipoint+2] = panels[i].rtl
            points_r[:,ipoint+3] = panels[i].rtr
            points_r[:,ipoint+4] = panels[i].rbr
        end

        cells_r = [MeshCell(PolyData.Polys(), 4*(i-1)+1 : 4*(i-1)+4) for i = 1:npanels]

        vtk_grid(vtmfile, points_r, cells_r) do vtkfile

            # add metadata
            for (key, value) in pairs(metadata)
                vtkfile[string(key)] = value
            end

        end
    end

    return nothing
end
