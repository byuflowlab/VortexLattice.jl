# TODO: Add grid based visualization
"""
    write_vtk(name, system; write_surfaces = true, write_wakes = false, kwargs...)

Write geometry from surfaces and/or wakes to Paraview files for visualization.

# Arguments
 - `name`: Base name for the generated files
 - `system`: VortexLatticeSystem object containing surfaces and/or wakes

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
"""


const _LINE_CELLS_H_CACHE = Dict{Tuple{Int,Int}, Any}()
const _LINE_CELLS_V_CACHE = Dict{Tuple{Int,Int}, Any}()
const _LINE_CELLS_HV_CACHE = Dict{Tuple{Int,Int}, Any}()
const _VERTS_CELLS_CACHE = Dict{Int, Any}()
const _TRAILING_LINE_CELLS_CACHE = Dict{Tuple{Int,Bool}, Any}()

function _line_cells_h(nc::Int, ns::Int)
    get!(_LINE_CELLS_H_CACHE, (nc, ns)) do
        li = LinearIndices((nc + 1, ns + 1))
        [MeshCell(PolyData.Lines(), SVector{2,Int}(li[i, j], li[i, j + 1])) for j = 1:ns for i = 1:nc]
    end
end

function _line_cells_v(nc::Int, ns::Int)
    get!(_LINE_CELLS_V_CACHE, (nc, ns)) do
        li = LinearIndices((nc + 1, ns + 1))
        [MeshCell(PolyData.Lines(), SVector{2,Int}(li[i, j], li[i + 1, j])) for j = 1:ns + 1 for i = 1:nc]
    end
end

function _line_cells_hv(nc::Int, ns::Int)
    get!(_LINE_CELLS_HV_CACHE, (nc, ns)) do
        lines_h = _line_cells_h(nc, ns)
        lines_v = _line_cells_v(nc, ns)
        n_h = length(lines_h)
        n_v = length(lines_v)
        lines_hv = Vector{eltype(lines_h)}(undef, n_h + n_v)
        copyto!(lines_hv, 1, lines_h, 1, n_h)
        copyto!(lines_hv, n_h + 1, lines_v, 1, n_v)
        lines_hv
    end
end

function _verts_cells(n::Int)
    get!(_VERTS_CELLS_CACHE, n) do
        [MeshCell(PolyData.Verts(), SVector{1,Int}(i)) for i = 1:n]
    end
end

function _trailing_line_cells(ns::Int; trailing_vortices::Bool)
    get!(_TRAILING_LINE_CELLS_CACHE, (ns, trailing_vortices)) do
        if trailing_vortices
            li = LinearIndices((2, ns + 1))
            lines_t = Vector{MeshCell{PolyData.Lines, SVector{2,Int}}}(undef, ns + (ns + 1))
            for j = 1:ns
                lines_t[j] = MeshCell(PolyData.Lines(), SVector{2,Int}(li[1, j], li[1, j + 1]))
            end
            for j = 1:ns + 1
                lines_t[ns + j] = MeshCell(PolyData.Lines(), SVector{2,Int}(li[1, j], li[2, j]))
            end
            lines_t
        else
            [MeshCell(PolyData.Lines(), SVector{2,Int}(j, j + 1)) for j = 1:ns]
        end
    end
end

function write_vtk(name::String, system::VortexLatticeSystem; write_surfaces = true, write_wakes = false, xhat = system.xhat[], 
    trailing_edge_list=nothing, kwargs...)

    trailing_edge_list = isnothing(trailing_edge_list) ? fill(true, length(system.surfaces)) : trailing_edge_list

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
"""
write_vtk(name, surfaces::AbstractVector{<:AbstractMatrix}, args...; kwargs...)

function write_vtk(name, surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}},
    properties=nothing; symmetric=nothing,
    trailing_edge_list=nothing, kwargs...)

    symmetric = isnothing(symmetric) ? fill(nothing, length(surfaces)) : symmetric
    trailing_edge_list = isnothing(trailing_edge_list) ? fill(true, length(surfaces)) : trailing_edge_list

    # create paraview multiblock file
    vtk_multiblock(name) do vtmfile
        # loop through all surfaces
        for i = 1:length(surfaces)
            # add paraview files corresponding to the surface to the multiblock file
            prop = isnothing(properties) ? nothing : properties[i]
            write_vtk!(vtmfile, surfaces[i], prop; trailing_edge=trailing_edge_list[i], symmetric=symmetric[i], kwargs...)
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
        symmetric=nothing,
        trailing_edge_list=nothing, kwargs...)

    symmetric = isnothing(symmetric) ? fill(nothing, length(surfaces)) : symmetric
    trailing_edge_list = isnothing(trailing_edge_list) ? fill(true, length(surfaces)) : trailing_edge_list

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
                trailing_edge = isempty(wakes[i]) && trailing_edge_list[i], symmetric=symmetric[i], kwargs..., trailing_vortices = false)

            # add paraview files corresponding to the wake to the multiblock file
            write_vtk!(vtmfile, wakes[i]; surface_circulation, symmetric=symmetric[i], kwargs...)
        end
    end

    return nothing
end

struct _SurfaceVTKScratch{TF}
    xyz::Array{TF,4}        # 3×(nc+1)×(ns+1)×1
    gamma_h::Matrix{TF}     # nc×ns
    cf_h::Array{TF,3}       # 3×nc×ns
    v_h::Array{TF,3}        # 3×nc×ns
    gamma_v::Matrix{TF}     # nc×(ns+1)
    cf_v::Array{TF,3}       # 3×nc×(ns+1)
    points_cp::Matrix{TF}   # 3×(nc*ns)
    normals::Matrix{TF}     # 3×(nc*ns)
    xyz_t::Matrix{TF}       # 3×2(ns+1), covers both trailing modes
    gamma_t::Vector{TF}     # 2ns+1, covers both trailing modes
    gamma_hv::Vector{TF}    # nc*ns + nc*(ns+1)
    cf_hv::Matrix{TF}       # 3 × (nc*ns + nc*(ns+1))
    v_hv::Matrix{TF}        # 3 × (nc*ns + nc*(ns+1))
end

function _SurfaceVTKScratch(TF::Type, nc::Int, ns::Int)
    N = nc * ns
    N_hv = N + nc * (ns + 1)
    _SurfaceVTKScratch{TF}(
        Array{TF,4}(undef, 3, nc+1, ns+1, 1),
        Matrix{TF}(undef, nc, ns),
        Array{TF,3}(undef, 3, nc, ns),
        Array{TF,3}(undef, 3, nc, ns),
        Matrix{TF}(undef, nc, ns+1),
        Array{TF,3}(undef, 3, nc, ns+1),
        Matrix{TF}(undef, 3, N),
        Matrix{TF}(undef, 3, N),
        Matrix{TF}(undef, 3, 2*(ns+1)),
        Vector{TF}(undef, 2*ns+1),
        Vector{TF}(undef, N_hv),
        Matrix{TF}(undef, 3, N_hv),
        Matrix{TF}(undef, 3, N_hv),
    )
end

struct _WakeVTKScratch{TF}
    xyz::Array{TF,4}           # 3×(nc+1)×(ns+1)×1
    gamma_h::Matrix{TF}        # nc×ns
    gamma_v::Matrix{TF}        # nc×(ns+1)
    circulation_hv::Vector{TF} # nc*ns + nc*(ns+1)
    xyz_t::Matrix{TF}          # 3×2(ns+1), covers both trailing modes
    gamma_t::Vector{TF}        # 2ns+1, covers both trailing modes
end

function _WakeVTKScratch(TF::Type, nc::Int, ns::Int)
    _WakeVTKScratch{TF}(
        Array{TF,4}(undef, 3, nc+1, ns+1, 1),
        Matrix{TF}(undef, nc, ns),
        Matrix{TF}(undef, nc, ns+1),
        Vector{TF}(undef, nc*ns + nc*(ns+1)),
        Matrix{TF}(undef, 3, 2*(ns+1)),
        Vector{TF}(undef, 2*ns+1),
    )
end

mutable struct _SystemVTKWriterState
    pvd
    block_name::String
    surface_scratch::Vector{_SurfaceVTKScratch}
    log_io::IOStream
end

function _init_system_vtk_writer(name::String, system::VortexLatticeSystem; overwrite::Bool=false)
    _parent, _base = splitdir(name)
    subdir = joinpath(_parent, _base)
    mkpath(subdir)
    block_name = joinpath(subdir, _base)
    pvd = paraview_collection(name; append=!overwrite)
    TF = eltype(eltype(system.surfaces[1]))
    scratch = [_SurfaceVTKScratch(TF, size(s, 1), size(s, 2)) for s in system.surfaces]

    # Per-step "*_log.csv" write disabled (2026-07-15) — nothing in this package or its tests
    # reads it back; it was pure write-only diagnostic output. `log_io` stays an IOStream (the
    # type _SystemVTKWriterState expects) — `open(devnull, "w")` doesn't work here (devnull is
    # a special IO sink object, not a path, and has no matching 2-arg open method) — use the
    # literal "/dev/null" path string instead, which does return a real IOStream.
    log_io = open("/dev/null", "w")

    return _SystemVTKWriterState(pvd, block_name, scratch, log_io)
end

function _append_system_vtk!(writer::_SystemVTKWriterState, system::VortexLatticeSystem, idx::Int, t::Real)
    vtm = vtk_multiblock(writer.block_name * "_$idx.vtm")
    for i = 1:length(system.surfaces)
        write_vtk!(vtm, system.surfaces[i], system.properties[i];
            trailing_edge = true,
            trailing_vortices = false,
            symmetric = system.symmetric[i],
            scratch = writer.surface_scratch[i])
    end
    writer.pvd[t] = vtm
    WriteVTK.LightXML.save_file(writer.pvd.xdoc, writer.pvd.path)
    return nothing
end

function _save_system_vtk_writer!(writer::_SystemVTKWriterState)
    WriteVTK.vtk_save(writer.pvd)
    close(writer.log_io)
    return nothing
end

"""
    write_vtk(name, system::VortexLatticeSystem, idx, t; overwrite=false)

Append one time step of `system`'s body geometry to a ParaView PVD collection.
The `.pvd` file is written at `name.pvd` and per-step multiblock files are
routed to a `name/` subdirectory as `name.{idx}.vtm`, with one block per
surface.
"""
function write_vtk(name::String, system::VortexLatticeSystem, idx::Int, t::Real; overwrite::Bool=false)
    writer = _init_system_vtk_writer(name, system; overwrite)
    _append_system_vtk!(writer, system, idx, t)
    _save_system_vtk_writer!(writer)

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
"""
function write_vtk!(vtmfile, surface::AbstractMatrix{<:SurfacePanel}, properties=nothing;
    symmetric = nothing,
    trailing_vortices = true,
    trailing_edge = true,
    xhat = SVector(1, 0, 0),
    wake_length = 10,
    wake_circulation = nothing,
    scratch::Union{Nothing,_SurfaceVTKScratch} = nothing)

    # get float type
    TF = eltype(eltype(surface))

    # check to make sure `symmetric` is provided if `properties` is provided
    @assert !(!isnothing(properties) && isnothing(symmetric)) "Keyword argument `symmetric` is required when optional argument `properties` is provided"

    # grid dimensions
    nc, ns = size(surface)
    N = length(surface)

    # extract geometry as a grid — reuse scratch buffer when available
    xyz = isnothing(scratch) ? Array{TF,4}(undef, 3, nc+1, ns+1, 1) : scratch.xyz
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

    # bound vortex geometries (h+v merged)
    lines_hv = _line_cells_hv(nc, ns)
    n_h = nc * ns
    n_v = nc * (ns + 1)

    gamma_hv_flat = nothing
    cf_hv_flat = nothing
    v_hv_flat = nothing

    if !isnothing(properties)
        gamma_h = isnothing(scratch) ? Matrix{TF}(undef, nc, ns) : scratch.gamma_h
        for i = 1:nc, j = 1:ns
            previous_gamma = i == 1 ? 0.0 : properties[i-1, j].gamma
            current_gamma = properties[i,j].gamma
            if symmetric && on_symmetry_plane(top_left(surface[i,j]), top_right(surface[i,j]))
                gamma_h[i,j] = 0.0
            else
                gamma_h[i,j] = current_gamma - previous_gamma
            end
        end

        cf_h = isnothing(scratch) ? Array{TF}(undef, 3, nc, ns) : scratch.cf_h
        for i = 1:nc, j = 1:ns
            if symmetric && on_symmetry_plane(top_left(surface[i,j]), top_right(surface[i,j]))
                cf_h[:,i,j] .= 0.0
            else
                cf_h[:,i,j] = properties[i,j].cfb
            end
        end

        v_h = isnothing(scratch) ? Array{TF}(undef, 3, nc, ns) : scratch.v_h
        for i = 1:nc, j = 1:ns
            v_h[:,i,j] = properties[i,j].velocity
        end

        gamma_v = isnothing(scratch) ? Matrix{TF}(undef, nc, ns+1) : scratch.gamma_v
        for i = 1:nc
            current_gamma = properties[i,1].gamma
            if symmetric && on_symmetry_plane(bottom_left(surface[i,1]), top_left(surface[i,1]))
                gamma_v[i,1] = 0.0
            else
                gamma_v[i,1] = -current_gamma
            end
            for j = 2:ns
                previous_gamma = current_gamma
                current_gamma = properties[i,j].gamma
                if symmetric && on_symmetry_plane(bottom_left(surface[i,j]), top_left(surface[i,j]))
                    gamma_v[i,j] = 0.0
                else
                    gamma_v[i,j] = previous_gamma - current_gamma
                end
            end
            previous_gamma = current_gamma
            if symmetric && on_symmetry_plane(bottom_right(surface[i,end]), top_right(surface[i,end]))
                gamma_v[i,end] = 0.0
            else
                gamma_v[i,end] = previous_gamma
            end
        end

        cf_v = isnothing(scratch) ? Array{TF}(undef, 3, nc, ns+1) : scratch.cf_v
        for i = 1:nc
            if symmetric && on_symmetry_plane(bottom_left(surface[i,1]), top_left(surface[i,1]))
                cf_v[:,i,1] .= 0.0
            else
                cf_v[:,i,1] = properties[i,1].cfl
            end
            for j = 2:ns
                previous_cf = properties[i,j-1].cfr
                current_cf = properties[i,j].cfl
                if symmetric && on_symmetry_plane(bottom_left(surface[i,j]), top_left(surface[i,j]))
                    cf_v[i,j] = 0.0
                else
                    cf_v[:,i,j] = previous_cf + current_cf
                end
            end
            if symmetric && on_symmetry_plane(bottom_right(surface[i,end]), top_right(surface[i,end]))
                cf_v[:,i,end] .= 0.0
            else
                cf_v[:,i,end] = properties[i,end].cfr
            end
        end

        gamma_hv = isnothing(scratch) ? Vector{TF}(undef, n_h + n_v) : scratch.gamma_hv
        cf_hv = isnothing(scratch) ? Matrix{TF}(undef, 3, n_h + n_v) : scratch.cf_hv
        v_hv = isnothing(scratch) ? Matrix{TF}(undef, 3, n_h + n_v) : scratch.v_hv
        copyto!(gamma_hv, 1, gamma_h, 1, n_h)
        copyto!(gamma_hv, n_h + 1, gamma_v, 1, n_v)
        copyto!(cf_hv, 1, cf_h, 1, 3 * n_h)
        copyto!(cf_hv, 3 * n_h + 1, cf_v, 1, 3 * n_v)
        copyto!(v_hv, 1, v_h, 1, 3 * n_h)
        fill!(view(v_hv, :, n_h + 1:n_h + n_v), zero(TF))
        gamma_hv_flat = gamma_hv
        cf_hv_flat = cf_hv
        v_hv_flat = v_hv
    end

    # bound vortices (h+v merged)
    vtk_grid(vtmfile, points, lines_hv) do vtkfile
        if !isnothing(properties)
            vtkfile["circulation"] = gamma_hv_flat
            vtkfile["velocity"] = v_hv_flat
            vtkfile["force"] = cf_hv_flat
        end
    end

    # --- control points ---

    points_cp = isnothing(scratch) ? Matrix{TF}(undef, 3, N) : scratch.points_cp
    for i = 1:N
        points_cp[:,i] = controlpoint(surface[i])
    end

    normals = isnothing(scratch) ? Matrix{TF}(undef, 3, N) : scratch.normals
    for i = 1:N
        normals[:,i] = normal(surface[i])
    end

    cells_cp = _verts_cells(N)

    vtk_grid(vtmfile, points_cp, cells_cp) do vtkfile
        vtkfile["normal"] = normals
    end

    # trailing edge and/or trailing vortices
    if trailing_vortices

        xyz_t_3d = isnothing(scratch) ? Array{TF}(undef, 3, 2, ns+1) : reshape(scratch.xyz_t, 3, 2, ns+1)
        for j = 1:ns
            xyz_t_3d[:,1,j] = bottom_left(surface[end,j])
            xyz_t_3d[:,2,j] = xyz_t_3d[:,1,j] + wake_length*xhat
        end
        xyz_t_3d[:,1,end] = bottom_right(surface[end,end])
        xyz_t_3d[:,2,end] = xyz_t_3d[:,1,end] + wake_length*xhat
        points_t = reshape(xyz_t_3d, 3, :)

        lines_t = _trailing_line_cells(ns; trailing_vortices=true)

        if !isnothing(properties)
            gamma_tv = isnothing(scratch) ? Vector{TF}(undef, ns+1) : view(scratch.gamma_t, ns+2:2ns+2)
            gamma_t  = isnothing(scratch) ? Vector{TF}(undef, 2ns+1) : view(scratch.gamma_t, 1:2ns+1)
            current_gamma = properties[end,1].gamma
            if symmetric && on_symmetry_plane(bottom_left(surface[end,1]))
                gamma_tv[1] = 0.0
            else
                gamma_tv[1] = -current_gamma
            end
            for j = 2:ns
                previous_gamma = current_gamma
                current_gamma = properties[end,j].gamma
                if symmetric && on_symmetry_plane(bottom_left(surface[end,j]))
                    gamma_tv[j] = 0.0
                else
                    gamma_tv[j] = previous_gamma - current_gamma
                end
            end
            previous_gamma = current_gamma
            if symmetric && on_symmetry_plane(bottom_right(surface[end,end]))
                gamma_tv[end] = 0.0
            else
                gamma_tv[end] = previous_gamma
            end
            fill!(view(gamma_t, 1:ns), zero(TF))
            gamma_t[ns+1:end] .= gamma_tv
        end

    else
        if trailing_edge
            points_t = isnothing(scratch) ? Matrix{TF}(undef, 3, ns+1) : view(scratch.xyz_t, :, 1:ns+1)
            for j = 1:ns
                points_t[:,j] = bottom_left(surface[end,j])
            end
            points_t[:,end] = bottom_right(surface[end,end])

            lines_t = _trailing_line_cells(ns; trailing_vortices=false)

            if !isnothing(properties)
                gamma_t = isnothing(scratch) ? Vector{TF}(undef, ns) : view(scratch.gamma_t, 1:ns)
                for j = 1:ns
                    if symmetric && on_symmetry_plane(bottom_left(surface[end,j]), bottom_right(surface[end,j]))
                        gamma_t[j] = 0.0
                    else
                        wake_gamma = isnothing(wake_circulation) ? zero(TF) : wake_circulation[j]
                        gamma_t[j] = wake_gamma - properties[end,j].gamma
                    end
                end
            end
        end
    end

    if trailing_vortices || trailing_edge
        vtk_grid(vtmfile, points_t, lines_t) do vtkfile
            if !isnothing(properties)
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
"""
function write_vtk!(vtmfile, wake::AbstractMatrix{<:WakePanel};
    symmetric,
    trailing_vortices = false,
    trailing_edge = true,
    xhat = SVector(1, 0, 0),
    wake_length = 10,
    surface_circulation = nothing,
    scratch::Union{Nothing,_WakeVTKScratch} = nothing)

    if isempty(wake)
        return vtmfile
    end

    TF = eltype(eltype(wake))
    nc, ns = size(wake)

    xyz = isnothing(scratch) ? Array{TF,4}(undef, 3, nc+1, ns+1, 1) : scratch.xyz
    for (i, I) in enumerate(CartesianIndices((nc+1, ns+1)))
        if I[1] <= nc && I[2] <= ns
            xyz[:, I, 1] = top_left(wake[I[1], I[2]])
        elseif I[1] <= nc && I[2] == ns + 1
            xyz[:, I, 1] = top_right(wake[I[1], I[2]-1])
        elseif I[1] == nc + 1 && I[2] <= ns
            xyz[:, I, 1] = bottom_left(wake[I[1]-1, I[2]])
        else
            xyz[:, I, 1] = bottom_right(wake[I[1]-1, I[2]-1])
        end
    end
    points = reshape(xyz, 3, :)

    lines_hv = _line_cells_hv(nc, ns)

    gamma_h = isnothing(scratch) ? Matrix{TF}(undef, nc, ns) : scratch.gamma_h
    for i in 1:nc, j = 1:ns
        previous_gamma = i == 1 ? (isnothing(surface_circulation) ? zero(TF) : surface_circulation[j]) : wake[i-1, j].gamma
        current_gamma = wake[i,j].gamma
        if symmetric && on_symmetry_plane(top_left(wake[nc, ns]), top_right(wake[nc, ns]))
            gamma_h[i,j] = 0.0
        else
            gamma_h[i,j] = current_gamma - previous_gamma
        end
    end
    gamma_h_flat = vec(gamma_h)

    gamma_v = isnothing(scratch) ? Matrix{TF}(undef, nc, ns+1) : scratch.gamma_v
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
    gamma_v_flat = vec(gamma_v)

    n_h = length(gamma_h_flat)
    n_v = length(gamma_v_flat)
    circulation_hv = isnothing(scratch) ? Vector{TF}(undef, n_h + n_v) : scratch.circulation_hv
    copyto!(circulation_hv, 1, gamma_h_flat, 1, n_h)
    copyto!(circulation_hv, n_h + 1, gamma_v_flat, 1, n_v)

    vtk_grid(vtmfile, points, lines_hv) do vtkfile
        vtkfile["circulation"] = circulation_hv
    end

    if trailing_vortices || trailing_edge
    if trailing_vortices

        xyz_t_3d = isnothing(scratch) ? Array{TF}(undef, 3, 2, ns+1) : reshape(scratch.xyz_t, 3, 2, ns+1)
        for j = 1:ns
            xyz_t_3d[:,1,j] = bottom_left(wake[end,j])
            xyz_t_3d[:,2,j] = xyz_t_3d[:,1,j] + wake_length*xhat
        end
        xyz_t_3d[:,1,end] = bottom_right(wake[end,end])
        xyz_t_3d[:,2,end] = xyz_t_3d[:,1,end] + wake_length*xhat
        points_t = reshape(xyz_t_3d, 3, :)

        lines_t = _trailing_line_cells(ns; trailing_vortices=true)

        gamma_tv = isnothing(scratch) ? Vector{TF}(undef, ns+1) : view(scratch.gamma_t, ns+2:2ns+2)
        gamma_t  = isnothing(scratch) ? Vector{TF}(undef, 2ns+1) : view(scratch.gamma_t, 1:2ns+1)
        current_gamma = wake[end,1].gamma
        if symmetric && on_symmetry_plane(bottom_left(wake[end,1]))
            gamma_tv[1] = 0.0
        else
            gamma_tv[1] = -current_gamma
        end
        for j = 2:ns
            previous_gamma = current_gamma
            current_gamma = wake[end,j].gamma
            if symmetric && on_symmetry_plane(bottom_left(wake[end,j]))
                gamma_tv[j] = 0.0
            else
                gamma_tv[j] = previous_gamma - current_gamma
            end
        end
        previous_gamma = current_gamma
        if symmetric && on_symmetry_plane(bottom_right(wake[end,end]))
            gamma_tv[end] = 0.0
        else
            gamma_tv[end] = previous_gamma
        end
        fill!(view(gamma_t, 1:ns), zero(TF))
        gamma_t[ns+1:end] .= gamma_tv

    else

        points_t = isnothing(scratch) ? Matrix{TF}(undef, 3, ns+1) : view(scratch.xyz_t, :, 1:ns+1)
        for j = 1:ns
            points_t[:,j] = bottom_left(wake[end,j])
        end
        points_t[:,end] = bottom_right(wake[end,end])

        lines_t = _trailing_line_cells(ns; trailing_vortices=false)

        gamma_t = isnothing(scratch) ? Vector{TF}(undef, ns) : view(scratch.gamma_t, 1:ns)
        for j = 1:ns
            if symmetric && on_symmetry_plane(bottom_left(wake[end,j]), bottom_right(wake[end,j]))
                gamma_t[j] = 0.0
            else
                gamma_t[j] = -wake[end,j].gamma
            end
        end

    end

    vtk_grid(vtmfile, points_t, lines_t) do vtkfile
        vtkfile["circulation"] = gamma_t
    end
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