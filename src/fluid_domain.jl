"""
    FluidDomainMonitor

Pre-allocated monitor for evaluating velocity and vorticity on a structured
rectilinear grid. Designed for zero per-step allocation after construction.

Two usage modes:
  1. During simulation — pass as an element of `monitors` in `simulate!`.
  2. Post-processing  — call `evaluate_fluid_domain_from_restarts!` to iterate
     over saved restart checkpoints without re-running the simulation.

`velocity` and `vorticity` always reflect the most recently evaluated step;
they are overwritten each call. Use `vtk_interval` to write `.vtr` files each
N steps, or supply a `callback` to `evaluate_fluid_domain_from_restarts!` to
capture results before they are overwritten.
"""
struct FluidDomainMonitor{TF,
        TX <: AbstractVector{TF},
        TY <: AbstractVector{TF},
        TZ <: AbstractVector{TF}}
    x::TX
    y::TY
    z::TZ
    velocity::Array{SVector{3,TF}, 3}    # (nx, ny, nz); overwritten each call
    vorticity::Array{SVector{3,TF}, 3}   # (nx, ny, nz); overwritten each call
    probes::FastMultipole.ProbeSystemStatic{TF}  # flat, column-major (x fastest)
    vtk_interval::Int                    # write every N steps (0 = never)
    name::String
    path::String
    pvd::WriteVTK.CollectionFile
end

"""
    FluidDomainMonitor(x, y, z; vtk_interval, name, path, TF)

Construct a `FluidDomainMonitor` on the rectilinear grid defined by coordinate
vectors (or ranges) `x`, `y`, `z`.

# Keyword arguments
- `vtk_interval::Int=0`: write a `.vtr` file every this many steps; 0 disables VTK output.
- `name::String="fluid_domain"`: base name for output files and PVD collection.
- `path::String="."`: output directory; created if it does not exist.
- `TF::Type=Float64`: floating-point type for all pre-allocated arrays.
"""
function FluidDomainMonitor(x::AbstractVector, y::AbstractVector, z::AbstractVector;
        vtk_interval::Int = 0,
        name::String = "fluid_domain",
        path::String = ".",
        TF::Type = Float64)

    xv = Vector{TF}(x)
    yv = Vector{TF}(y)
    zv = Vector{TF}(z)
    nx, ny, nz = length(xv), length(yv), length(zv)
    n = nx * ny * nz

    velocity  = Array{SVector{3,TF}}(undef, nx, ny, nz)
    vorticity = Array{SVector{3,TF}}(undef, nx, ny, nz)
    probes    = FastMultipole.ProbeSystem(n, TF)

    # populate probe positions once; column-major order (x index varies fastest)
    i = 0
    @inbounds for iz in 1:nz, iy in 1:ny, ix in 1:nx
        i += 1
        probes.position[i] = SVector{3,TF}(xv[ix], yv[iy], zv[iz])
    end

    mkpath(path)
    pvd = WriteVTK.paraview_collection(joinpath(path, name))

    return FluidDomainMonitor(xv, yv, zv, velocity, vorticity, probes,
        vtk_interval, name, path, pvd)
end

"""
    (m::FluidDomainMonitor)(system, wake, i_step; fmm_wake_args, fmm_vehicle_args)

Monitor callable. Evaluates total velocity (induced + freestream) and vorticity
(curl of velocity via velocity-gradient hessian) at every grid point and stores
results in `m.velocity` and `m.vorticity`. Writes a VTK file when
`vtk_interval > 0` and `i_step % vtk_interval == 0`.
"""
function (m::FluidDomainMonitor{TF})(system::System, wake::PanelParticleWake,
        i_step::Int; fmm_wake_args=(), fmm_vehicle_args=()) where TF

    FastMultipole.reset!(m.probes)

    np   = FLOWVPM.get_np(wake.pfield)
    tef  = WakeBufferRings(wake)
    nfil = FastMultipole.get_n_bodies(tef)

    wake_kw    = _fmm_kwargs(wake.fmm_wake[],    wake.pfield.useGPU)
    vehicle_kw = _fmm_kwargs(wake.fmm_vehicle[], wake.pfield.useGPU)

    # wake particles + panel-buffer filaments → fluid domain probes
    if np > 0 && nfil > 0
        fmm!((m.probes,), (wake.pfield, tef);
            wake_kw..., hessian=SVector{1}(true), fmm_wake_args...)
    elseif np > 0
        fmm!((m.probes,), (wake.pfield,);
            wake_kw..., hessian=SVector{1}(true), fmm_wake_args...)
    elseif nfil > 0
        fmm!((m.probes,), (tef,);
            wake_kw..., hessian=SVector{1}(true), fmm_wake_args...)
    end

    # body surfaces → fluid domain probes
    fmm!((m.probes,), (system,);
        vehicle_kw..., hessian=SVector{1}(true), fmm_vehicle_args...)

    Vinf = freestream_velocity(system.freestream[])
    nx   = length(m.x)
    nxy  = nx * length(m.y)

    @inbounds for k in eachindex(m.probes.gradient)
        V = m.probes.gradient[k] + Vinf
        H = m.probes.hessian[k]
        ω = SVector{3,TF}(H[3,2] - H[2,3], H[1,3] - H[3,1], H[2,1] - H[1,2])
        km1 = k - 1
        ix  = mod(km1, nx)  + 1
        iy  = mod(km1 ÷ nx, length(m.y)) + 1
        iz  = km1 ÷ nxy + 1
        m.velocity[ix, iy, iz]  = V
        m.vorticity[ix, iy, iz] = ω
    end

    if m.vtk_interval > 0 && i_step % m.vtk_interval == 0
        _write_fluid_domain_vtk(m, i_step, wake.pfield.t, Vinf)
    end

    return nothing
end

"""
    evaluate_fluid_domain!(monitor, system, wake, i_step=0; kwargs...)

Evaluate the fluid domain at the current system/wake state. Convenience wrapper
around the monitor callable; use after a manual [`restore_restart!`](@ref) call.
"""
function evaluate_fluid_domain!(monitor::FluidDomainMonitor, system::System,
        wake::PanelParticleWake, i_step::Int=0; kwargs...)
    monitor(system, wake, i_step; kwargs...)
    return monitor
end

"""
    evaluate_fluid_domain_from_restarts!(monitor, system, wake, frames, restart_name;
        indices, callback, kwargs...)

Iterate over saved restart checkpoints and evaluate the fluid domain at each
step, reusing the pre-allocated `monitor` arrays to keep peak memory at
`O(n_grid_points)` regardless of how many steps are processed.

# Arguments
- `indices`: subset of checkpoint step indices to process. Accepts any
  `AbstractVector{Int}` (e.g. `1:5:200`) or `nothing` to process all
  available checkpoints. Steps not found in the checkpoint index are silently
  skipped.
- `callback(monitor, idx, t)`: optional function called after each step while
  `monitor.velocity` and `monitor.vorticity` still hold that step's values.
  Use this to copy results or accumulate statistics before they are overwritten.
"""
function evaluate_fluid_domain_from_restarts!(
        monitor::FluidDomainMonitor,
        system::System,
        wake::PanelParticleWake,
        frames::AbstractVector{<:ReferenceFrame},
        restart_name::String;
        indices::Union{Nothing, AbstractVector{Int}} = nothing,
        callback::Union{Nothing, Function} = nothing,
        kwargs...)

    all_idx, all_t = _read_restart_checkpoint_index(restart_name)
    isempty(all_idx) && return monitor

    if isnothing(indices)
        step_idx = all_idx
        step_t   = all_t
    else
        idx_set  = Set(indices)
        mask     = [i ∈ idx_set for i in all_idx]
        step_idx = all_idx[mask]
        step_t   = all_t[mask]
    end

    for (idx, t) in zip(step_idx, step_t)
        restore_restart!(system, wake, frames, restart_name; idx=idx)
        monitor(system, wake, idx; kwargs...)
        isnothing(callback) || callback(monitor, idx, t)
    end

    return monitor
end

# ---- internal helpers -------------------------------------------------------

function _write_fluid_domain_vtk(m::FluidDomainMonitor{TF}, i_step::Int,
        t::Real, Vinf::SVector{3,TF}) where TF
    filename = joinpath(m.path, m.name * "_$(lpad(i_step, 8, '0'))")
    nx, ny, nz = length(m.x), length(m.y), length(m.z)

    vmag  = [norm(m.velocity[i,j,k])  for i in 1:nx, j in 1:ny, k in 1:nz]
    omag  = [norm(m.vorticity[i,j,k]) for i in 1:nx, j in 1:ny, k in 1:nz]
    v_ind = m.velocity .- Ref(Vinf)
    vmag_ind = [norm(v_ind[i,j,k]) for i in 1:nx, j in 1:ny, k in 1:nz]

    vts = WriteVTK.vtk_grid(filename, m.x, m.y, m.z)
    vts["velocity",       WriteVTK.VTKPointData()] = reinterpret(reshape, TF, m.velocity)
    vts["velocity_mag",   WriteVTK.VTKPointData()] = vmag
    vts["v_induced",      WriteVTK.VTKPointData()] = reinterpret(reshape, TF, v_ind)
    vts["v_induced_mag",  WriteVTK.VTKPointData()] = vmag_ind
    vts["vorticity",      WriteVTK.VTKPointData()] = reinterpret(reshape, TF, m.vorticity)
    vts["vorticity_mag",  WriteVTK.VTKPointData()] = omag
    m.pvd[t] = vts
    WriteVTK.LightXML.save_file(m.pvd.xdoc, m.pvd.path)
    return nothing
end

function _read_restart_checkpoint_index(name::String)
    base = endswith(name, ".pvd") ? name[1:end-4] : name
    endswith(base, "_bodies") && (base = base[1:end-7])
    index_file = joinpath(base * "_restart_vtk", "index.tsv")
    isfile(index_file) || error("Restart checkpoint index not found: $(index_file)")

    indices = Int[]
    times   = Float64[]
    for line in Iterators.drop(eachline(index_file), 1)
        entry = split(line, '\t')
        length(entry) == 3 || continue
        push!(indices, parse(Int, entry[1]))
        push!(times,   parse(Float64, entry[2]))
    end
    return indices, times
end
