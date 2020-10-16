"""
    write_vtk(name, panels::AbstractVector{<:Horseshoe}; kwargs...)

Writes the geometry to Paraview files for visualization.

# Keyword Arguments:
 - `xhat=[1, 0, 0]`: direction in which trailing vortices are shed
 - `mirror=false`: creates a mirror image of the geometry across the X-Z axis
 - `wake_length=10`: distance to extend the trailing vortices
 - `metadata=Dict()`: dictionary of metadata to include in generated files
"""
write_vtk

function write_vtk(name, panels::AbstractVector{<:Horseshoe}; xhat=SVector(1, 0, 0), mirror=false, wake_length=10, metadata=Dict())

    TF = eltype(eltype(panels))

    if mirror
        panels = vcat(panels, reflect_panel.(panels))
    end

    npanels = length(panels)

    vtk_multiblock(name) do vtmfile

        # --- horseshoe vortices ---

        points_h = Matrix{TF}(undef, 3, 6*npanels)
        for i = 1:npanels
            ipoint = 6*(i - 1)
            points_h[:,ipoint+1] = panels[i].rl + SVector(panels[i].xl_te, 0, 0) + wake_length*xhat
            points_h[:,ipoint+2] = panels[i].rl + SVector(panels[i].xl_te, 0, 0)
            points_h[:,ipoint+3] = panels[i].rl
            points_h[:,ipoint+4] = panels[i].rr
            points_h[:,ipoint+5] = panels[i].rr + SVector(panels[i].xr_te, 0, 0)
            points_h[:,ipoint+6] = panels[i].rr + SVector(panels[i].xr_te, 0, 0) + wake_length*xhat
        end

        cells_h = [MeshCell(PolyData.Lines(), 6*(i-1)+1 : 6*(i-1)+6) for i = 1:npanels]

        vtk_grid(vtmfile, points_h, cells_h) do vtkfile

            # add metadata
            for (key, value) in pairs(metadata)
                vtkfile[string(key)] = value
            end

        end

        # --- control points ---

        points_cp = Matrix{TF}(undef, 3, npanels)
        for i = 1:npanels
            ipoint = i
            points_cp[:,ipoint] = panels[i].rcp
        end

        cells_cp = [MeshCell(PolyData.Verts(), [i]) for i = 1:npanels]

        vtk_grid(vtmfile, points_cp, cells_cp) do vtkfile

            # add metadata
            for (key, value) in pairs(metadata)
                vtkfile[string(key)] = value
            end

            # add normal
            data = Matrix{TF}(undef, 3, npanels)
            for i = 1:length(panels)
                data[:,i] = normal(panels[i])
            end
            vtkfile["normal"] = data

        end
    end

    return nothing
end

function write_vtk(name, panels::AbstractVector{<:Ring}; xhat=SVector(1, 0, 0), wake_length=10, mirror=false, metadata=Dict())

    TF = eltype(eltype(panels))

    if mirror
        panels = vcat(panels, reflect_panel.(panels))
    end

    npanels = length(panels)
    ntrailing = count([panel.trailing for panel in panels])

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

        # --- horseshoe vortices ---

        points_h = Matrix{TF}(undef, 3, 4*ntrailing)
        ipoint = 0
        for i = 1:length(panels)
            if panels[i].trailing
                points_h[:,ipoint+1] = panels[i].rbl + wake_length*xhat
                points_h[:,ipoint+2] = panels[i].rbl
                points_h[:,ipoint+3] = panels[i].rbr
                points_h[:,ipoint+4] = panels[i].rbr + wake_length*xhat
                ipoint += 4
            end
        end

        cells_h = [MeshCell(PolyData.Lines(), 4*(i-1)+1 : 4*(i-1)+4) for i = 1:ntrailing]

        vtk_grid(vtmfile, points_h, cells_h) do vtkfile

            # add metadata
            for (key, value) in pairs(metadata)
                vtkfile[string(key)] = value
            end

        end

        # --- control points ---

        points_cp = Matrix{TF}(undef, 3, npanels)
        for i = 1:npanels
            ipoint = i
            points_cp[:,ipoint] = panels[i].rcp
        end

        cells_cp = [MeshCell(PolyData.Verts(), [i]) for i = 1:npanels]

        vtk_grid(vtmfile, points_cp, cells_cp) do vtkfile

            # add metadata
            for (key, value) in pairs(metadata)
                vtkfile[string(key)] = value
            end

            # add normal
            data = Matrix{TF}(undef, 3, npanels)
            for i = 1:length(panels)
                data[:,i] = normal(panels[i])
            end
            vtkfile["normal"] = data

        end
    end

    return nothing
end
