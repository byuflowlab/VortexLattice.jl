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

function write_vtk(name::String, system::System; write_surfaces = true, write_wakes = false, xhat = system.xhat[], 
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
"""
function write_vtk(name, surface_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}},
    property_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}},
    wake_history::AbstractVector{<:AbstractVector{<:AbstractMatrix}}, dt;
    symmetric = nothing, kwargs...)

    symmetric = isnothing(symmetric) ? fill(nothing, length(surface_history[1])) : symmetric
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
    write_vtk(name, system::System, idx, t; overwrite=false)

Append one time step of `system`'s body geometry to a ParaView PVD collection.
The `.pvd` file is written at `name.pvd` and per-step multiblock files are
routed to a `name/` subdirectory as `name.{idx}.vtm`, with one block per
surface.
"""
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

mutable struct _RestartCheckpointWriterState
    dir::String
    index_file::String
end

function _restart_checkpoint_prefix(name::String)
    base = endswith(name, ".pvd") ? name[1:end-4] : name
    if endswith(base, "_bodies")
        return base[1:end-7]
    end
    return base
end

function _init_restart_checkpoint_writer(name::String; overwrite::Bool=false)
    prefix = _restart_checkpoint_prefix(name)
    dir = prefix * "_restart_vtk"
    mkpath(dir)
    index_file = joinpath(dir, "index.tsv")
    if overwrite || !isfile(index_file)
        open(index_file, "w") do io
            println(io, "idx\ttime\tfile")
        end
    end
    return _RestartCheckpointWriterState(dir, index_file)
end

function _flatten_reference_frames(frames::AbstractVector{<:ReferenceFrame}, TF)
    nframes = length(frames)
    frame_x = Vector{TF}(undef, 3nframes)
    frame_v = Vector{TF}(undef, 3nframes)
    frame_axis = Vector{TF}(undef, 3nframes)
    frame_omega = Vector{TF}(undef, nframes)
    frame_R = Vector{TF}(undef, 9nframes)
    frame_Rp2g = Vector{TF}(undef, 9nframes)
    for i in 1:nframes
        i3 = 3(i - 1)
        i9 = 9(i - 1)
        frame = frames[i]
        frame_x[i3+1:i3+3] .= frame.x
        frame_v[i3+1:i3+3] .= frame.v
        frame_axis[i3+1:i3+3] .= frame.ω_axis
        frame_omega[i] = frame.ω
        frame_R[i9+1:i9+9] .= vec(frame.R)
        frame_Rp2g[i9+1:i9+9] .= vec(frame.Rp2g)
    end
    return frame_x, frame_v, frame_axis, frame_omega, frame_R, frame_Rp2g
end

function _flatten_wake_panels(wakes, TF)
    dims = Vector{TF}(undef, 2length(wakes))
    n_panels = sum(length, wakes)
    rtl = Vector{TF}(undef, 3n_panels)
    rtr = Vector{TF}(undef, 3n_panels)
    rbl = Vector{TF}(undef, 3n_panels)
    rbr = Vector{TF}(undef, 3n_panels)
    core = Vector{TF}(undef, n_panels)
    gamma = Vector{TF}(undef, n_panels)
    ip = 0
    for isurf in eachindex(wakes)
        wake = wakes[isurf]
        nr, ns = size(wake)
        dims[2isurf-1] = nr
        dims[2isurf] = ns
        for j in 1:ns, i in 1:nr
            ip += 1
            i3 = 3(ip - 1)
            panel = wake[i, j]
            rtl[i3+1:i3+3] .= panel.rtl
            rtr[i3+1:i3+3] .= panel.rtr
            rbl[i3+1:i3+3] .= panel.rbl
            rbr[i3+1:i3+3] .= panel.rbr
            core[ip] = panel.core_size
            gamma[ip] = panel.gamma
        end
    end
    return dims, rtl, rtr, rbl, rbr, core, gamma
end

function _flatten_wsl(wsl, TF)
    dims = Vector{TF}(undef, length(wsl))
    npts = sum(length, wsl)
    data = Vector{TF}(undef, 3npts)
    ip = 0
    for isurf in eachindex(wsl)
        dims[isurf] = length(wsl[isurf])
        for p in wsl[isurf]
            ip += 1
            i3 = 3(ip - 1)
            data[i3+1:i3+3] .= p
        end
    end
    return dims, data
end

function _flatten_wake_velocities(wake_velocities, TF)
    dims = Vector{TF}(undef, 2length(wake_velocities))
    npts = 0
    for V in wake_velocities
        npts += length(V)
    end
    data = Vector{TF}(undef, 3npts)
    ip = 0
    for isurf in eachindex(wake_velocities)
        V = wake_velocities[isurf]
        nr, ns = size(V)
        dims[2isurf-1] = nr
        dims[2isurf] = ns
        for j in 1:ns, i in 1:nr
            ip += 1
            i3 = 3(ip - 1)
            data[i3+1:i3+3] .= V[i, j]
        end
    end
    return dims, data
end

function _flatten_grids(grids, TF)
    dims = Vector{TF}(undef, 3length(grids))
    nvals = 0
    for grid in grids
        nvals += length(grid)
    end
    data = Vector{TF}(undef, nvals)
    i0 = 0
    for isurf in eachindex(grids)
        grid = grids[isurf]
        n1, n2, n3 = size(grid)
        dims[3isurf-2] = n1
        dims[3isurf-1] = n2
        dims[3isurf] = n3
        n = length(grid)
        data[i0+1:i0+n] .= vec(grid)
        i0 += n
    end
    return dims, data
end

function _append_restart_checkpoint!(writer::_RestartCheckpointWriterState, idx::Int, t::Real,
        system::System, wake::PanelParticleWake, frames::AbstractVector)
    TF = eltype(system.Γ)
    fname = "restart_" * lpad(string(idx), 8, '0') * ".vtp"
    fpath = joinpath(writer.dir, fname)

    fs = system.freestream[]
    frame_x, frame_v, frame_axis, frame_omega, frame_R, frame_Rp2g =
        _flatten_reference_frames(frames, TF)
    wake_dims, wake_rtl, wake_rtr, wake_rbl, wake_rbr, wake_core, wake_gamma =
        _flatten_wake_panels(wake.wakes, TF)
    wsl_dims, wsl_data = _flatten_wsl(wake.wake_shedding_locations, TF)
    wake_vel_dims, wake_vel_data = _flatten_wake_velocities(wake.wake_velocities, TF)
    grid_dims, grid_data = _flatten_grids(system.grids, TF)
    prev_bottom_gamma = vcat((copy(g) for g in wake.prev_bottom_gamma)...)

    points = zeros(TF, 3, 1)
    cells = [WriteVTK.MeshCell(WriteVTK.PolyData.Verts(), 1:1)]
    vtk_grid(fpath[1:end-4], points, cells) do vtkfile
        vtkfile["restart_idx", WriteVTK.VTKFieldData()] = TF[idx]
        vtkfile["restart_time", WriteVTK.VTKFieldData()] = TF[t]
        vtkfile["freestream", WriteVTK.VTKFieldData()] = TF[fs.Vinf, fs.alpha, fs.beta, fs.Omega[1], fs.Omega[2], fs.Omega[3]]
        vtkfile["Gamma", WriteVTK.VTKFieldData()] = copy(system.Γ)
        vtkfile["dGamma_dt", WriteVTK.VTKFieldData()] = copy(system.dΓdt)
        vtkfile["nwake_active", WriteVTK.VTKFieldData()] = TF.(wake.nwake)
        vtkfile["wake_overflowed", WriteVTK.VTKFieldData()] = TF[wake.overflowed[] ? 1 : 0]
        vtkfile["prev_bottom_gamma", WriteVTK.VTKFieldData()] = prev_bottom_gamma
        vtkfile["frame_count", WriteVTK.VTKFieldData()] = TF[length(frames)]
        vtkfile["frame_x", WriteVTK.VTKFieldData()] = frame_x
        vtkfile["frame_v", WriteVTK.VTKFieldData()] = frame_v
        vtkfile["frame_omega_axis", WriteVTK.VTKFieldData()] = frame_axis
        vtkfile["frame_omega", WriteVTK.VTKFieldData()] = frame_omega
        vtkfile["frame_R", WriteVTK.VTKFieldData()] = frame_R
        vtkfile["frame_Rp2g", WriteVTK.VTKFieldData()] = frame_Rp2g
        vtkfile["wake_dims", WriteVTK.VTKFieldData()] = wake_dims
        vtkfile["wake_rtl", WriteVTK.VTKFieldData()] = wake_rtl
        vtkfile["wake_rtr", WriteVTK.VTKFieldData()] = wake_rtr
        vtkfile["wake_rbl", WriteVTK.VTKFieldData()] = wake_rbl
        vtkfile["wake_rbr", WriteVTK.VTKFieldData()] = wake_rbr
        vtkfile["wake_core", WriteVTK.VTKFieldData()] = wake_core
        vtkfile["wake_gamma", WriteVTK.VTKFieldData()] = wake_gamma
        vtkfile["wsl_dims", WriteVTK.VTKFieldData()] = wsl_dims
        vtkfile["wake_shedding_locations", WriteVTK.VTKFieldData()] = wsl_data
        vtkfile["wake_vel_dims", WriteVTK.VTKFieldData()] = wake_vel_dims
        vtkfile["wake_velocities", WriteVTK.VTKFieldData()] = wake_vel_data
        vtkfile["grid_dims", WriteVTK.VTKFieldData()] = grid_dims
        vtkfile["grid_data", WriteVTK.VTKFieldData()] = grid_data
    end

    open(writer.index_file, "a") do io
        println(io, "$(idx)\t$(t)\t$(fname)")
    end
    return nothing
end

"""
    write_restart_checkpoint(name, idx, t, system, wake, frames; overwrite=false)

Write a full VTK restart checkpoint that can be resumed later with
[`restore_restart!`](@ref).
"""
function write_restart_checkpoint(name::String, idx::Int, t::Real,
        system::System, wake::PanelParticleWake, frames::AbstractVector;
        overwrite::Bool=false)
    writer = _init_restart_checkpoint_writer(name; overwrite)
    _append_restart_checkpoint!(writer, idx, t, system, wake, frames)
    return nothing
end

function _latest_restart_checkpoint_entry(index_file::String)
    lines = readlines(index_file)
    length(lines) <= 1 && error("No restart checkpoint entries found in $(index_file)")
    entry = split(lines[end], '\t')
    length(entry) == 3 || error("Malformed restart checkpoint index entry: $(lines[end])")
    return (idx=parse(Int, entry[1]), t=parse(Float64, entry[2]), file=entry[3])
end

function _find_restart_checkpoint_entry(index_file::String, idx::Int)
    for line in Iterators.drop(eachline(index_file), 1)
        entry = split(line, '\t')
        length(entry) == 3 || continue
        if parse(Int, entry[1]) == idx
            return (idx=idx, t=parse(Float64, entry[2]), file=entry[3])
        end
    end
    error("Restart checkpoint entry not found for idx=$(idx) in $(index_file)")
end

function _restart_getdata(field_data, key::String)
    haskey(field_data, key) || error("Missing restart checkpoint field: $(key)")
    value = field_data[key]
    return value isa ReadVTK.VTKDataArray ? ReadVTK.get_data(value) : value
end

"""
    read_restart_checkpoint_info(name; idx=nothing)

Read a VTK restart checkpoint and return `(idx, t, data)`.
"""
function read_restart_checkpoint_info(name::String; idx::Union{Nothing,Int}=nothing)
    prefix = _restart_checkpoint_prefix(name)
    dir = prefix * "_restart_vtk"
    index_file = joinpath(dir, "index.tsv")
    isfile(index_file) || error("Restart checkpoint index not found: $(index_file)")

    entry = isnothing(idx) ? _latest_restart_checkpoint_entry(index_file) :
        _find_restart_checkpoint_entry(index_file, idx)
    fpath = joinpath(dir, entry.file)
    isfile(fpath) || error("Restart checkpoint file not found: $(fpath)")

    vtk = ReadVTK.VTKFile(fpath)
    field_data = ReadVTK.get_field_data(vtk)
    data = Dict{String,Any}()
    for (key, value) in field_data
        data[key] = ReadVTK.get_data(value)
    end

    return (idx=entry.idx, t=entry.t, data=data)
end

function _restore_pfield_from_vtp!(pfield, vtp_path::String, TF)
    vtk = ReadVTK.VTKFile(vtp_path)
    X = TF.(ReadVTK.get_points(vtk))  # 3 × np
    np = size(X, 2)
    np > 0 || return
    np <= size(pfield.particles, 2) || error("Particle buffer overflow during restart restore")
    pd = ReadVTK.get_point_data(vtk)
    pfield.particles[FLOWVPM.X_INDEX, 1:np] .= X
    pfield.particles[FLOWVPM.GAMMA_INDEX, 1:np] .= reshape(TF.(ReadVTK.get_data(pd["gamma"])), 3, np)
    pfield.particles[FLOWVPM.SIGMA_INDEX, 1:np] .= vec(TF.(ReadVTK.get_data(pd["sigma"])))
    pfield.particles[FLOWVPM.VOL_INDEX, 1:np] .= vec(TF.(ReadVTK.get_data(pd["vol"])))
    pfield.particles[FLOWVPM.CIRCULATION_INDEX, 1:np] .= vec(TF.(ReadVTK.get_data(pd["circulation"])))
    pfield.particles[FLOWVPM.U_INDEX, 1:np] .= reshape(TF.(ReadVTK.get_data(pd["velocity"])), 3, np)
    pfield.particles[FLOWVPM.VORTICITY_INDEX, 1:np] .= reshape(TF.(ReadVTK.get_data(pd["vorticity"])), 3, np)
    pfield.particles[FLOWVPM.J_INDEX, 1:np] .= reshape(TF.(ReadVTK.get_data(pd["velocity_gradient"])), 9, np)
    pfield.np = np
end

"""
    restore_restart!(system, wake, frames, name; idx=nothing)

Load a VTK restart checkpoint and apply it in one call. Returns `(idx, t)` for
the resolved restart step.
"""
function restore_restart!(system::System{TF}, wake::PanelParticleWake{TF},
        frames::AbstractVector{<:ReferenceFrame{TF}}, name::String;
        idx::Union{Nothing,Int}=nothing) where TF
    info = read_restart_checkpoint_info(name; idx)
    data = info.data

    fs = TF.(_restart_getdata(data, "freestream"))
    system.freestream[] = Freestream(fs[1], fs[2], fs[3], SVector{3,TF}(fs[4], fs[5], fs[6]))
    system.Γ .= TF.(_restart_getdata(data, "Gamma"))
    system.dΓdt .= TF.(_restart_getdata(data, "dGamma_dt"))

    wake.nwake .= Int.(round.(TF.(_restart_getdata(data, "nwake_active"))))
    wake.overflowed[] = Int(round(TF(_restart_getdata(data, "wake_overflowed")[1]))) != 0

    prev = TF.(_restart_getdata(data, "prev_bottom_gamma"))
    i0 = 0
    for isurf in eachindex(wake.prev_bottom_gamma)
        n = length(wake.prev_bottom_gamma[isurf])
        wake.prev_bottom_gamma[isurf] .= view(prev, i0+1:i0+n)
        i0 += n
    end

    nframes = Int(round(TF(_restart_getdata(data, "frame_count")[1])))
    nframes == length(frames) || error("Frame count mismatch during restart restore")
    frame_x = TF.(_restart_getdata(data, "frame_x"))
    frame_v = TF.(_restart_getdata(data, "frame_v"))
    frame_axis = TF.(_restart_getdata(data, "frame_omega_axis"))
    frame_omega = TF.(_restart_getdata(data, "frame_omega"))
    frame_R = TF.(_restart_getdata(data, "frame_R"))
    frame_Rp2g = TF.(_restart_getdata(data, "frame_Rp2g"))
    for i in 1:nframes
        i3 = 3(i - 1)
        i9 = 9(i - 1)
        frame_old = frames[i]
        x = SVector{3,TF}(view(frame_x, i3+1:i3+3))
        v = SVector{3,TF}(view(frame_v, i3+1:i3+3))
        omega_axis = SVector{3,TF}(view(frame_axis, i3+1:i3+3))
        omega = frame_omega[i]
        R = SMatrix{3,3,TF,9}(Tuple(view(frame_R, i9+1:i9+9)))
        Rp2g = SMatrix{3,3,TF,9}(Tuple(view(frame_Rp2g, i9+1:i9+9)))
        frames[i] = ReferenceFrame(x, v, omega_axis, omega, R, Rp2g,
            frame_old.name, frame_old.parent_index, frame_old.child_index, frame_old.dependent_index)
    end

    wake_dims = Int.(round.(TF.(_restart_getdata(data, "wake_dims"))))
    rtl = TF.(_restart_getdata(data, "wake_rtl"))
    rtr = TF.(_restart_getdata(data, "wake_rtr"))
    rbl = TF.(_restart_getdata(data, "wake_rbl"))
    rbr = TF.(_restart_getdata(data, "wake_rbr"))
    core = TF.(_restart_getdata(data, "wake_core"))
    gamma = TF.(_restart_getdata(data, "wake_gamma"))
    ip = 0
    for isurf in eachindex(wake.wakes)
        nr = wake_dims[2isurf-1]
        ns = wake_dims[2isurf]
        size(wake.wakes[isurf]) == (nr, ns) || error("Wake shape mismatch during restart restore")
        for j in 1:ns, i in 1:nr
            ip += 1
            i3 = 3(ip - 1)
            wake.wakes[isurf][i, j] = WakePanel{TF}(
                SVector{3,TF}(view(rtl, i3+1:i3+3)),
                SVector{3,TF}(view(rtr, i3+1:i3+3)),
                SVector{3,TF}(view(rbl, i3+1:i3+3)),
                SVector{3,TF}(view(rbr, i3+1:i3+3)),
                core[ip], gamma[ip])
        end
    end

    wsl_dims = Int.(round.(TF.(_restart_getdata(data, "wsl_dims"))))
    wsl_data = TF.(_restart_getdata(data, "wake_shedding_locations"))
    ip = 0
    for isurf in eachindex(wake.wake_shedding_locations)
        n = wsl_dims[isurf]
        length(wake.wake_shedding_locations[isurf]) == n || error("Wake shedding shape mismatch during restart restore")
        for i in 1:n
            ip += 1
            i3 = 3(ip - 1)
            wake.wake_shedding_locations[isurf][i] = SVector{3,TF}(view(wsl_data, i3+1:i3+3))
        end
    end

    wake_vel_dims = Int.(round.(TF.(_restart_getdata(data, "wake_vel_dims"))))
    wake_vel_data = TF.(_restart_getdata(data, "wake_velocities"))
    ip = 0
    for isurf in eachindex(wake.wake_velocities)
        nr = wake_vel_dims[2isurf-1]
        ns = wake_vel_dims[2isurf]
        size(wake.wake_velocities[isurf]) == (nr, ns) || error("Wake velocity shape mismatch during restart restore")
        for j in 1:ns, i in 1:nr
            ip += 1
            i3 = 3(ip - 1)
            wake.wake_velocities[isurf][i, j] = SVector{3,TF}(view(wake_vel_data, i3+1:i3+3))
        end
    end

    grid_dims = Int.(round.(TF.(_restart_getdata(data, "grid_dims"))))
    grid_data = TF.(_restart_getdata(data, "grid_data"))
    i0 = 0
    for isurf in eachindex(system.grids)
        n1 = grid_dims[3isurf-2]
        n2 = grid_dims[3isurf-1]
        n3 = grid_dims[3isurf]
        size(system.grids[isurf]) == (n1, n2, n3) || error("Grid shape mismatch during restart restore")
        n = n1 * n2 * n3
        system.grids[isurf] .= reshape(view(grid_data, i0+1:i0+n), n1, n2, n3)
        i0 += n
        update_surface_panels!(system.surfaces[isurf], system.grids[isurf];
            ratios=system.ratios[isurf],
            fcore=(c, Δs) -> system.core_size)
    end

    wake.pfield.particles[:, :] .= zero(TF)
    wake.pfield.np = 0
    particles_vtp = joinpath(name * "_wake_particles",
        basename(name) * "_wake_particles_$(info.idx - 1).vtp")
    if isfile(particles_vtp)
        _restore_pfield_from_vtp!(wake.pfield, particles_vtp, TF)
    end
    wake.pfield.t = info.t
    wake.pfield.nt = info.idx

    return (idx=info.idx, t=info.t)
end

function _init_system_vtk_writer(name::String, system::System; overwrite::Bool=false)
    _parent, _base = splitdir(name)
    subdir = joinpath(_parent, _base)
    mkpath(subdir)
    block_name = joinpath(subdir, _base)
    pvd = paraview_collection(name; append=!overwrite)
    TF = eltype(eltype(system.surfaces[1]))
    scratch = [_SurfaceVTKScratch(TF, size(s, 1), size(s, 2)) for s in system.surfaces]

    log_path = name * "_log.csv"
    write_header = overwrite || !isfile(log_path)
    log_io = open(log_path, overwrite ? "w" : "a")
    if write_header
        nsurf = length(system.surfaces)
        nwake_cols = join(("nwake_$i" for i in 1:nsurf), ",")
        println(log_io, "step,time,Vinf,alpha,beta,Omega_x,Omega_y,Omega_z,wake_overflowed,$nwake_cols")
        flush(log_io)
    end

    return _SystemVTKWriterState(pvd, block_name, scratch, log_io)
end

function _append_system_vtk!(writer::_SystemVTKWriterState, system::System, idx::Int, t::Real)
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

function _append_step_log!(writer::_SystemVTKWriterState, system::System, wake::PanelParticleWake, i_step::Int, t::Real)
    fs = system.freestream[]
    Ω = fs.Omega
    overflowed = wake.overflowed[] ? 1 : 0
    nwake = join(wake.nwake, ",")
    println(writer.log_io, "$i_step,$t,$(fs.Vinf),$(fs.alpha),$(fs.beta),$(Ω[1]),$(Ω[2]),$(Ω[3]),$overflowed,$nwake")
    flush(writer.log_io)
    return nothing
end

function _save_system_vtk_writer!(writer::_SystemVTKWriterState)
    WriteVTK.vtk_save(writer.pvd)
    close(writer.log_io)
    return nothing
end

function write_vtk(name::String, system::System, idx::Int, t::Real; overwrite::Bool=false)
    writer = _init_system_vtk_writer(name, system; overwrite)
    _append_system_vtk!(writer, system, idx, t)
    _save_system_vtk_writer!(writer)

    return nothing
end

"""
    write_vtk(name, wake::PanelParticleWake, idx, t; overwrite=false)

Append one time step of the buffer-overflow wake to ParaView PVD collections.
Writes the panel portion to `name.pvd` (per-step `.vtm` under `name/`, one VTS
block per surface) and the particle portion to `name_particles.pvd` (per-step
`.vtp` under `name_particles/`).
"""
mutable struct _WakeVTKWriterState
    panel_pvd
    panel_block_name::String
    particles_pvd
    particles_block::String
    particle_cells
    particle_empty_points
    particle_empty_cells
    wake_scratch::Vector{_WakeVTKScratch}
end

function _init_wake_vtk_writer(name::String, wake::PanelParticleWake; overwrite::Bool=false)
    _parent, _base = splitdir(name)

    panel_subdir = joinpath(_parent, _base)
    mkpath(panel_subdir)
    panel_block_name = joinpath(panel_subdir, _base)
    panel_pvd = paraview_collection(name; append=!overwrite)

    particles_pvd_name = joinpath(_parent, _base * "_particles")
    mkpath(particles_pvd_name)
    particles_block = joinpath(particles_pvd_name, _base * "_particles")
    particles_pvd = paraview_collection(particles_pvd_name; append=!overwrite)

    particle_cells = [WriteVTK.MeshCell(WriteVTK.PolyData.Verts(), 1:1)]
    particle_empty_points = zeros(Float64, 3, 0)
    particle_empty_cells = Vector{WriteVTK.MeshCell{WriteVTK.PolyData.Verts, UnitRange{Int}}}()

    TF = eltype(eltype(wake.wakes[1]))
    wake_scratch = [_WakeVTKScratch(TF, wake.nwakerows, size(w, 2)) for w in wake.wakes]

    return _WakeVTKWriterState(panel_pvd, panel_block_name, particles_pvd, particles_block,
        particle_cells, particle_empty_points, particle_empty_cells, wake_scratch)
end

function _append_wake_vtk!(writer::_WakeVTKWriterState, wake::PanelParticleWake, idx::Int, t::Real)
    # panel wake
    vtm = vtk_multiblock(writer.panel_block_name * "_$idx.vtm")
    for i = 1:length(wake.wakes)
        n = wake.nwake[i]
        n == 0 && continue
        wake_view = view(wake.wakes[i], 1:n, :)
        write_vtk!(vtm, wake_view; symmetric=false, trailing_vortices=false,
            scratch = writer.wake_scratch[i])
    end
    writer.panel_pvd[t] = vtm

    # particle wake
    np = wake.pfield.np
    X = view(wake.pfield.particles, FLOWVPM.X_INDEX, 1:np)
    cells = writer.particle_cells
    cells[1] = WriteVTK.MeshCell(WriteVTK.PolyData.Verts(), 1:max(np, 1))

    vtp_filename = writer.particles_block * "_$idx.vtp"
    if np > 0
        vtp = WriteVTK.vtk_grid(vtp_filename, X, cells)
        vtp["gamma", WriteVTK.VTKPointData()] = view(wake.pfield.particles, FLOWVPM.GAMMA_INDEX, 1:np)
        vtp["sigma", WriteVTK.VTKPointData()] = view(wake.pfield.particles, FLOWVPM.SIGMA_INDEX, 1:np)
        vtp["vol", WriteVTK.VTKPointData()] = view(wake.pfield.particles, FLOWVPM.VOL_INDEX, 1:np)
        vtp["circulation", WriteVTK.VTKPointData()] = view(wake.pfield.particles, FLOWVPM.CIRCULATION_INDEX, 1:np)
        vtp["velocity", WriteVTK.VTKPointData()] = view(wake.pfield.particles, FLOWVPM.U_INDEX, 1:np)
        vtp["vorticity", WriteVTK.VTKPointData()] = view(wake.pfield.particles, FLOWVPM.VORTICITY_INDEX, 1:np)
        vtp["velocity_gradient", WriteVTK.VTKPointData()] = view(wake.pfield.particles, FLOWVPM.J_INDEX, 1:np)
    else
        vtp = WriteVTK.vtk_grid(vtp_filename, writer.particle_empty_points, writer.particle_empty_cells)
    end

    writer.particles_pvd[t] = vtp
    WriteVTK.LightXML.save_file(writer.panel_pvd.xdoc, writer.panel_pvd.path)
    WriteVTK.LightXML.save_file(writer.particles_pvd.xdoc, writer.particles_pvd.path)

    return nothing
end

function _save_wake_vtk_writer!(writer::_WakeVTKWriterState)
    WriteVTK.vtk_save(writer.panel_pvd)
    WriteVTK.vtk_save(writer.particles_pvd)
    return nothing
end

function write_vtk(name::String, wake::PanelParticleWake, idx::Int, t::Real; overwrite::Bool=false)
    writer = _init_wake_vtk_writer(name, wake; overwrite)
    _append_wake_vtk!(writer, wake, idx, t)
    _save_wake_vtk_writer!(writer)

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