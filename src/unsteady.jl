struct DerivativesMonitor{TF}
    CFalpha::Vector{SVector{3,TF}}
    CFbeta::Vector{SVector{3,TF}}
    CFp::Vector{SVector{3,TF}}
    CFq::Vector{SVector{3,TF}}
    CFr::Vector{SVector{3,TF}}
    CMalpha::Vector{SVector{3,TF}}
    CMbeta::Vector{SVector{3,TF}}
    CMp::Vector{SVector{3,TF}}
    CMq::Vector{SVector{3,TF}}
    CMr::Vector{SVector{3,TF}}
end

function DerivativesMonitor(nt::Int, TF=Float64)
    CFalpha = zeros(SVector{3,TF}, nt)
    CFbeta = zeros(SVector{3,TF}, nt)
    CFp = zeros(SVector{3,TF}, nt)
    CFq = zeros(SVector{3,TF}, nt)
    CFr = zeros(SVector{3,TF}, nt)
    CMalpha = zeros(SVector{3,TF}, nt)
    CMbeta = zeros(SVector{3,TF}, nt)
    CMp = zeros(SVector{3,TF}, nt)
    CMq = zeros(SVector{3,TF}, nt)
    CMr = zeros(SVector{3,TF}, nt)

    return DerivativesMonitor{TF}(CFalpha, CFbeta, 
        CFp, CFq, CFr, CMalpha, CMbeta, CMp, CMq, CMr)
end

function (monitor::DerivativesMonitor)(system::System, wake, i_step::Int)
    CF, CM = stability_derivatives(system)
    monitor.CFalpha[i_step + 1] = CF.alpha
    monitor.CFbeta[i_step + 1] = CF.beta
    monitor.CFp[i_step + 1] = CF.p
    monitor.CFq[i_step + 1] = CF.q
    monitor.CFr[i_step + 1] = CF.r
    monitor.CMalpha[i_step + 1] = CM.alpha
    monitor.CMbeta[i_step + 1] = CM.beta
    monitor.CMp[i_step + 1] = CM.p
    monitor.CMq[i_step + 1] = CM.q
    monitor.CMr[i_step + 1] = CM.r
end

struct ForcesMonitor{TF,F}
    CF::Vector{SVector{3,TF}}
    CM::Vector{SVector{3,TF}}
    frame::F
end

function ForcesMonitor(nt::Int, TF=Float64; frame=Body())
    CF = zeros(SVector{3,TF}, nt)
    CM = zeros(SVector{3,TF}, nt)

    return ForcesMonitor{TF,typeof(frame)}(CF, CM, frame)
end

function (monitor::ForcesMonitor)(system::System, wake, i_step::Int)
    CF, CM = body_forces(system.surfaces, system.properties,
                            system.reference[], system.freestream[], 
                            system.symmetric, monitor.frame)
    monitor.CF[i_step + 1] = CF
    monitor.CM[i_step + 1] = CM
end

struct PanelForcesMonitor{TF}
    CF::Array{TF, 4}
    surface_index::Int
    ns::Int
    nc::Int
end

function PanelForcesMonitor(nt::Int, system::System, TF=Float64; surface_index=1)
    nc, ns = size(system.surfaces[surface_index])
    CF = zeros(TF, 3, nc, ns, nt)
    return PanelForcesMonitor{TF}(CF, surface_index, ns, nc)
end

function (monitor::PanelForcesMonitor)(system::System, wake, i_step::Int)
    CF = view(monitor.CF, 1:3, 1:monitor.nc, 1:monitor.ns, i_step + 1)
    ns = monitor.ns
    nc = monitor.nc
    properties = system.properties[monitor.surface_index]
    for j in 1:ns
        for i in 1:nc
            CF[:, i, j] .= properties[i, j].cfb
        end
    end
end

struct LiftingLineCoefficientsMonitor{TF, F}
    CF::Vector{Array{TF, 3}}
    CM::Vector{Array{TF, 3}}
    ns::Vector{Int}
    frame::F
    normalized::Bool
end

function LiftingLineCoefficientsMonitor(nt::Int, system::System, TF=Float64; frame=Body(), normalized=true)
    nsurf = length(system.surfaces)
    ns = zeros(Int, nsurf)
    CF = Vector{Array{TF, 3}}(undef, nsurf)
    CM = Vector{Array{TF, 3}}(undef, nsurf)
    for isurf in 1:nsurf
        ns[isurf] = size(system.surfaces[isurf], 2)
        CF[isurf] = zeros(TF, 3, ns[isurf], nt)
        CM[isurf] = zeros(TF, 3, ns[isurf], nt)
    end
    return LiftingLineCoefficientsMonitor{TF,typeof(frame)}(CF, CM, ns, frame, normalized)
end

function (monitor::LiftingLineCoefficientsMonitor)(system::System, wake, i_step::Int)
    cf, cm = lifting_line_coefficients(system; frame=monitor.frame, normalized=monitor.normalized)
    for isurf in 1:length(system.surfaces)
        CF = view(monitor.CF[isurf], 1:3, 1:monitor.ns[isurf], i_step + 1)
        CM = view(monitor.CM[isurf], 1:3, 1:monitor.ns[isurf], i_step + 1)
        CF .= cf[isurf]
        CM .= cm[isurf]
    end
end

struct FrameForcesMonitor{TF,F}
    CF::Vector{SVector{3,TF}}
    CM::Vector{SVector{3,TF}}
    frame::F
end

function simulate!(system::System{TF}, frames::AbstractVector{<:ReferenceFrame}, maneuver!::Function, Vinf::Function, t_range;
        particle_trailing_methods=fill(OverlapPPS(1.3, 2), length(system.surfaces)),
        particle_unsteady_methods=fill(OverlapPPS(1.3, 2), length(system.surfaces)),
        wake_args=(), kwargs...) where TF

    # construct particle field
    n_particles_per_step = get_max_particles(system, particle_trailing_methods, particle_unsteady_methods)
    wake = ParticleField(n_particles_per_step * length(t_range), TF; Uinf=Vinf, wake_args...)

    # begin simulation
    simulate!(system, wake, frames, maneuver!, Vinf, t_range;
        particle_trailing_methods, particle_unsteady_methods, kwargs...)

    return wake
end

"""
    simulate!(system, frames, maneuver!, Vinf, t_range, Ωinf;
              wake_type=PanelParticleWake, nwakerows, max_particles, kwargs...)

Auto-construct a `PanelParticleWake` for `system` and drive the unsteady
simulation. The `system` must have been allocated with `nw=fill(nwakerows, ...)`
so that the panel buffer has the expected row count.
"""
function simulate!(system::System, frames::AbstractVector{<:ReferenceFrame},
        maneuver!::Function, Vinf::Function, t_range,
        Ωinf::Function; wake_type::Type{PanelParticleWake},
        nwakerows::Int=size(system.wakes[1], 1),
        max_particles::Int=10_000,
        eta::Real=0.3,
        fmm::FLOWVPM.FMM=FLOWVPM.FMM(; p=20),
        fmm_wake::Union{Nothing, FLOWVPM.FMM}=nothing,
        fmm_vehicle::Union{Nothing, FLOWVPM.FMM}=nothing,
        method_trailing::WakeSheddingMethod=OverlapPPS(1.3, 2),
        method_unsteady::WakeSheddingMethod=OverlapPPS(1.3, 2),
        vpm_kwargs::NamedTuple=NamedTuple(),
        kwargs...)

    wake = PanelParticleWake(system;
        nwakerows=nwakerows,
        max_particles=max_particles,
        eta=eta,
        fmm=fmm,
        fmm_wake=fmm_wake,
        fmm_vehicle=fmm_vehicle,
        method_trailing=method_trailing,
        method_unsteady=method_unsteady,
        vpm_kwargs=vpm_kwargs)

    simulate!(system, wake, frames, maneuver!, Vinf, t_range, Ωinf; kwargs...)

    return wake
end

function force!(val, force_val::Nothing)
    return val
end

function force!(val, force_val)
    val .= force_val
end

function check_for_nans(system::System)
    if any(isnan.(system.Γ))
        error("NaN detected in circulation")
    end
    if any(isnan.(system.w))
        error("NaN detected in normal velocity")
    end
    for isurf in 1:length(system.surfaces)
        if any(isnan.(system.V[isurf]))
            error("NaN detected in velocity on surface $isurf")
        end
        if any(isnan.(system.Vcp[isurf]))
            error("NaN detected in velocity due to surface motion on surface $isurf")
        end
        if any(isnan.(system.Vh[isurf]))
            error("NaN detected in velocity due to heave on surface $isurf")
        end
        if any(isnan.(system.Vv[isurf]))
            error("NaN detected in velocity due to pitch on surface $isurf")
        end
        if any(isnan.(system.Vte[isurf]))
            error("NaN detected in velocity due to trailing edge motion on surface $isurf")
        end
        wake = system.wakes[isurf]
        nw_active = system.nwake[isurf]
        for j in 1:size(wake, 2), i in 1:nw_active
            panel = wake[i, j]
            if isnan(panel.core_size)
                error("NaN detected in wake panel core size on surface $isurf, panel ($i,$j)")
            end
            if isnan(panel.gamma)
                error("NaN detected in wake panel circulation on surface $isurf, panel ($i,$j)")
            end
            if any(isnan.(panel.rtl))
                error("NaN detected in wake panel rtl on surface $isurf, panel ($i,$j)")
            end
            if any(isnan.(panel.rbl))
                error("NaN detected in wake panel rbl on surface $isurf, panel ($i,$j)")
            end
            if any(isnan.(panel.rtr))
                error("NaN detected in wake panel rtr on surface $isurf, panel ($i,$j)")
            end
            if any(isnan.(panel.rbr))
                error("NaN detected in wake panel rbr on surface $isurf, panel ($i,$j)")
            end
        end
    end
end

function check_for_nans(wake::ParticleField)
    if any(isnan.(wake.particles[1:3, :]))
        @show wake.particles[:, 1:10]
        error("NaN detected in wake particle positions")
    end
    if any(isnan.(wake.particles[4:6, :]))
        @show wake.particles[:, 1:10]
        error("NaN detected in wake particle gamma")
    end
    if any(isnan.(wake.particles[7,:]))
        @show wake.particles[:, 1:10]
        error("NaN detected in wake particle core size")
    end
    if any(isnan.(wake.particles[10:12,:]))
        @show wake.particles[:, 1:10]
        error("NaN detected in wake particle velocity")
    end
    if any(isnan.(wake.particles[16:24,:]))
        @show wake.particles[:, 1:10]
        error("NaN detected in wake particle J")
    end
end

function Base.isnan(val::SVector{N,TF}) where {N,TF}
    return any(isnan.(val))
end


function _init_simulate!(system, wake::PanelParticleWake, frames, name, path,
        restart_from, restart_idx, write_restart, vtk_interval, t_range)
    @assert wake.nwakerows >= 1 "PanelParticleWake requires nwakerows >= 1, got $(wake.nwakerows)"

    if !isnothing(path) && !isdir(path)
        mkpath(path)
    end
    if isnothing(path)
        write_restart = false
        vtk_interval  = 0
    end

    restart_state = nothing
    if !isnothing(restart_from)
        restart_state = restore_restart!(system, wake, frames, restart_from; idx=restart_idx)
    end

    body_name       = isnothing(path) ? nothing : joinpath(path, name * "_bodies")
    checkpoint_base = isnothing(path) ? name    : joinpath(path, name)
    overwrite_first = isnothing(restart_state)
    body_writer = isnothing(body_name) ? nothing :
        _init_system_vtk_writer(body_name, system; overwrite=overwrite_first)
    wake_writer = isnothing(path) ? nothing :
        _init_wake_vtk_writer(joinpath(path, name * "_wake"), wake; overwrite=overwrite_first)

    if write_restart && !isnothing(body_name) && isnothing(restart_state)
        write_restart_checkpoint(checkpoint_base, 0, t_range[1], system, wake, frames; overwrite=true)
    end

    return (;restart_state, body_name, checkpoint_base, body_writer, wake_writer,
             write_restart, vtk_interval)
end

function _warm_start!(system, wake::PanelParticleWake, frames, maneuver!, Vinf, Ωinf,
        t_range, trailing_edge_filaments, ref, symmetric, surface_id, wake_finite_core,
        xhat, derivatives, fmm_wake_args)
    t0  = t_range[1]
    dt0 = t_range[2] - t_range[1]

    Vcp = system.Vcp; Vh = system.Vh; Vv = system.Vv; Vte = system.Vte
    for isurf in 1:length(system.surfaces)
        Vcp[isurf] .= Ref(zero(eltype(Vcp[isurf])))
        Vh[isurf]  .= Ref(zero(eltype(Vh[isurf])))
        Vv[isurf]  .= Ref(zero(eltype(Vv[isurf])))
        Vte[isurf] .= Ref(zero(eltype(Vte[isurf])))
    end
    system.w .= zero(eltype(system.w))
    system.nwake .= 0
    maneuver!(frames, system, wake, t0)
    system.freestream[] = velocity_to_freestream(Vinf(t0), Ωinf(t0))
    kinematic_velocity!(Vcp, Vh, Vv, Vte, system.surfaces, frames; skip_top_level=false)
    update_wake_shedding_locations!(system.wakes, system.wake_shedding_locations,
        system.surfaces, ref, system.freestream[], dt0, nothing, Vte, system.nwake, wake.eta)
    influence_coefficients!(system.AIC, system.surfaces;
        symmetric, wake_shedding_locations=system.wake_shedding_locations,
        surface_id, trailing_vortices=system.trailing_vortices, xhat,
        force_finite_core=fill(true, length(system.surfaces)))
    update_trailing_edge_coefficients!(system.AIC, system.surfaces;
        symmetric, wake_shedding_locations=system.wake_shedding_locations,
        trailing_vortices=system.trailing_vortices)
    system.fAIC[] = lu(system.AIC)

    # solve twice: once without wake influence, once with the first shed row
    for pass in 1:2
        system.w .= zero(eltype(system.w))
        if derivatives
            normal_velocity_derivatives!(system.w, system.dw, system.surfaces, system.wakes,
                ref, system.freestream[]; additional_velocity=nothing, Vcp, symmetric,
                nwake=system.nwake, surface_id, wake_finite_core,
                trailing_vortices=system.trailing_vortices, xhat, include_wakes=false)
            circulation_derivatives!(system.Γ, system.dΓ, system.AIC, system.w, system.dw)
        else
            normal_velocity!(system.w, system.surfaces, system.wakes, ref, system.freestream[];
                additional_velocity=nothing, Vcp, symmetric, nwake=system.nwake,
                surface_id, wake_finite_core, trailing_vortices=system.trailing_vortices,
                xhat, include_wakes=false)
            ldiv!(system.Γ, system.fAIC[], system.w)
        end
        pass == 2 && break
        # pre-seed the first wake row so step 0 sees wake influence (makes dΓdt ≈ 0)
        shed_wake!(wake, system, dt0, system.Γ)
        system.nwake .= wake.nwake
        reset!(wake)
        update_trailing_edge_filaments!(trailing_edge_filaments, system.surfaces, system.Γ)
        wake.pfield.SFS(wake.pfield, FLOWVPM.BeforeUJ())
        wake_on_all!(system, wake, trailing_edge_filaments; fmm_wake_args...)
    end
end

function _simulate_step!(system, wake::PanelParticleWake, frames, trailing_edge_filaments,
        config, i_step, start_step, t, t_range)
    (;ref, symmetric, surface_id, wake_finite_core, xhat,
      recalculate_influence_matrix, derivatives, fmm_wake_args, fmm_vehicle_args,
      polars, frames_index, monitors, vtk_interval, body_name, write_restart,
      checkpoint_base, body_writer, wake_writer, Vinf, Ωinf, maneuver!, verbose,
      Γ_wake, dΓdt_wake) = config

    verbose && println("\tstep $(i_step)/$(length(t_range)-1) at time $(t) | particles: $(FLOWVPM.get_np(wake.pfield))")

    #------- reset system -------#

    for isurf in 1:length(system.surfaces)
        system.Vcp[isurf] .= Ref(zero(eltype(system.Vcp[isurf])))
        system.Vh[isurf]  .= Ref(zero(eltype(system.Vh[isurf])))
        system.Vv[isurf]  .= Ref(zero(eltype(system.Vv[isurf])))
        system.Vte[isurf] .= Ref(zero(eltype(system.Vte[isurf])))
    end
    system.w .= zero(eltype(system.w))
    for i in eachindex(system.dw)
        system.dw[i] .= zero(eltype(system.dw[i]))
    end
    reset!(wake)

    #------- controls -------#

    maneuver!(frames, system, wake, t)

    #------- external flow -------#

    vinf = Vinf(t)
    fs   = velocity_to_freestream(vinf, Ωinf(t))
    system.freestream[] = fs
    kinematic_velocity!(system.Vcp, system.Vh, system.Vv, system.Vte, system.surfaces, frames; skip_top_level=false)

    dt = i_step == length(t_range) - 1 ? t_range[end] - t_range[end-1] : t_range[i_step + 2] - t_range[i_step + 1]

    if i_step == 0
        update_wake_shedding_locations!(system.wakes, system.wake_shedding_locations,
            system.surfaces, ref, fs, dt, nothing, system.Vte, system.nwake, wake.eta)
    end

    #------- wake coupling -------#

    system.nwake .= wake.nwake
    update_TE!(wake, system)
    update_trailing_edge_filaments!(trailing_edge_filaments, system.surfaces, system.Γ)
    if FLOWVPM.get_np(wake.pfield) > 0 || any(>(0), wake.nwake)
        wake.pfield.SFS(wake.pfield, FLOWVPM.BeforeUJ())
        wake_on_all!(system, wake, trailing_edge_filaments; fmm_wake_args...)
    end

    #------- AIC + solve -------#

    AIC = system.AIC
    Γ   = system.Γ
    trailing_vortices = system.trailing_vortices

    if recalculate_influence_matrix || i_step == start_step
        influence_coefficients!(AIC, system.surfaces;
            symmetric, wake_shedding_locations=system.wake_shedding_locations,
            surface_id, trailing_vortices, xhat,
            force_finite_core=fill(true, length(system.surfaces)))
        update_trailing_edge_coefficients!(AIC, system.surfaces;
            symmetric, wake_shedding_locations=system.wake_shedding_locations,
            trailing_vortices)
        system.fAIC[] = lu(AIC)
    else
        update_wake_shedding_locations!(system.wakes, system.wake_shedding_locations,
            system.surfaces, ref, fs, dt, nothing, system.Vte, system.nwake, wake.eta)
    end

    if derivatives
        normal_velocity_derivatives!(system.w, system.dw, system.surfaces, system.wakes,
            ref, fs; additional_velocity=nothing, Vcp=system.Vcp, symmetric,
            nwake=system.nwake, surface_id, wake_finite_core,
            trailing_vortices, xhat, include_wakes=false)
        system.dΓdt .= .-Γ
        circulation_derivatives!(Γ, system.dΓ, AIC, system.w, system.dw)
    else
        normal_velocity!(system.w, system.surfaces, system.wakes, ref, fs;
            additional_velocity=nothing, Vcp=system.Vcp, symmetric, nwake=system.nwake,
            surface_id, wake_finite_core, trailing_vortices, xhat, include_wakes=false)
        system.dΓdt .= .-Γ
        ldiv!(Γ, system.fAIC[], system.w)
    end
    system.dΓdt .+= Γ
    system.dΓdt ./= dt

    vehicle_on_all!(system, wake, trailing_edge_filaments; fmm_vehicle_args...)

    #------- near-field forces + viscous -------#

    if derivatives
        near_field_forces_derivatives!(system.properties, system.dproperties,
            system.surfaces, system.wakes, ref, fs, Γ, system.dΓ;
            dΓdt=system.dΓdt, additional_velocity=nothing, Vh=system.Vh, Vv=system.Vv,
            symmetric, nwake=system.nwake, surface_id, wake_finite_core,
            wake_shedding_locations=system.wake_shedding_locations,
            trailing_vortices, xhat, calculate_vlm_induced=false)
    else
        near_field_forces!(system.properties, system.surfaces, system.wakes,
            ref, fs, Γ; dΓdt=system.dΓdt, additional_velocity=nothing,
            Vh=system.Vh, Vv=system.Vv, symmetric, nwake=system.nwake, surface_id,
            wake_finite_core, wake_shedding_locations=system.wake_shedding_locations,
            trailing_vortices, xhat, calculate_vlm_induced=false)
    end

    Γ_wake .= Γ
    if !isnothing(polars)
        viscous!(system.properties, Γ_wake, system.surfaces, system.grids, frames,
            frames_index, polars, ref, dt)
    end
    dΓdt_wake .+= Γ_wake
    dΓdt_wake ./= dt

    #------- monitors -------#

    system.near_field_analysis[] = true
    for monitor in monitors
        monitor(system, wake, i_step)
    end

    #------- propagate system -------#

    _seed_wake_velocity!(wake, fs)
    propagate!(wake, dt; Vinf=vinf, relax=true)
    store_trailing_edge!(system.wake_shedding_locations, system.surfaces)
    propagate_kinematics!(system, frames, dt)

    idx = i_step == length(t_range) - 1 ? i_step + 1 : i_step + 2
    system.freestream[] = velocity_to_freestream(Vinf(t_range[idx]), Ωinf(t_range[idx]))

    for isurf in 1:length(system.surfaces)
        system.Vcp[isurf] .= Ref(zero(eltype(system.Vcp[isurf])))
        system.Vh[isurf]  .= Ref(zero(eltype(system.Vh[isurf])))
        system.Vv[isurf]  .= Ref(zero(eltype(system.Vv[isurf])))
        system.Vte[isurf] .= Ref(zero(eltype(system.Vte[isurf])))
    end
    kinematic_velocity!(system.Vcp, system.Vh, system.Vv, system.Vte, system.surfaces, frames; skip_top_level=false)

    update_wake_shedding_locations_unsteady!(system.wakes, system.wake_shedding_locations,
        system.surfaces, ref, system.freestream[], dt, nothing,
        system.Vte, system.nwake, wake.eta; sync_panels=false)

    shed_wake!(wake, system, dt, Γ_wake)

    if vtk_interval > 0 && !isnothing(body_writer) && mod(i_step, vtk_interval) == 0
        _append_system_vtk!(body_writer, system, i_step, t)
        _append_wake_vtk!(wake_writer, wake, i_step, t)
        _append_step_log!(body_writer, system, wake, i_step, t)
        if write_restart && !isnothing(body_name) && i_step < length(t_range) - 1
            write_restart_checkpoint(checkpoint_base, i_step + 1, t, system, wake, frames;
                overwrite=false)
        end
    end
end

function _simulate_step!(system, wake::ParticleField, frames, trailing_edge_filaments,
        config, i_step, t, t_range)
    (;eta, trailing_vortices, derivatives, recalculate_influence_matrix,
      fmm_wake_args, fmm_vehicle_args, particle_trailing_methods, particle_unsteady_methods,
      polars, frames_index, monitors, path, name, vtk_args, vtk_postshed,
      Vinf, Ωinf, maneuver!, verbose, Γ_wake, dΓdt_wake) = config

    verbose && println("\tstep $(i_step)/$(length(t_range)-1) at time $(t) | particles: $(FLOWVPM.get_np(wake))")

    #------- reset system -------#

    for isurf in 1:length(system.surfaces)
        system.Vcp[isurf] .= Ref(zero(eltype(system.Vcp[isurf])))
        system.Vh[isurf]  .= Ref(zero(eltype(system.Vh[isurf])))
        system.Vv[isurf]  .= Ref(zero(eltype(system.Vv[isurf])))
        system.Vte[isurf] .= Ref(zero(eltype(system.Vte[isurf])))
        system.V[isurf]   .= Ref(zero(eltype(system.V[isurf])))
    end
    system.w .= zero(eltype(system.w))
    for i in eachindex(system.dw)
        system.dw[i] .= zero(eltype(system.dw[i]))
    end
    FLOWVPM._reset_particles(wake)
    FLOWVPM._reset_particles_sfs(wake)

    #------- controls -------#

    maneuver!(frames, system, wake, t)
    kinematic_velocity!(system.Vcp, system.Vh, system.Vv, system.Vte, system.surfaces, frames; skip_top_level=false)

    #------- aerodynamics -------#

    ref              = system.reference[]
    symmetric        = system.symmetric
    surface_id       = system.surface_id
    wake_finite_core = system.wake_finite_core
    xhat             = system.xhat[]
    symmetric       .= false

    wakes                  = system.wakes
    wake_shedding_locations = system.wake_shedding_locations
    nwake                  = system.nwake
    Γ                      = system.Γ
    fs                     = system.freestream[]

    dt = i_step == length(t_range) - 1 ? t_range[end] - t_range[end-1] : t_range[i_step + 2] - t_range[i_step + 1]
    if i_step == 0
        dt = t_range[2] - t_range[1]
        update_wake_shedding_locations!(wakes, wake_shedding_locations,
            system.surfaces, ref, fs, dt, nothing, system.Vte, nwake, eta)
        initial_wake_panels!(wakes, wake_shedding_locations, system.surfaces, Γ, eta)
    end

    update_trailing_edge_filaments!(trailing_edge_filaments, system.surfaces, Γ)
    wake.SFS(wake, FLOWVPM.BeforeUJ())
    wake_on_all!(system, wake, trailing_edge_filaments; fmm_wake_args...)

    #------- AIC + solve -------#

    AIC = system.AIC

    if recalculate_influence_matrix || i_step == 0
        influence_coefficients!(AIC, system.surfaces;
            symmetric, wake_shedding_locations,
            surface_id, trailing_vortices, xhat,
            force_finite_core=fill(true, length(system.surfaces)))
        update_trailing_edge_coefficients!(AIC, system.surfaces;
            symmetric, wake_shedding_locations, trailing_vortices)
        system.fAIC[] = lu(AIC)
    else
        update_wake_shedding_locations!(wakes, wake_shedding_locations,
            system.surfaces, ref, fs, dt, nothing, system.Vte, nwake, eta)
    end

    if derivatives
        normal_velocity_derivatives!(system.w, system.dw, system.surfaces, wakes,
            ref, fs; additional_velocity=nothing, Vcp=system.Vcp, symmetric, nwake,
            surface_id, wake_finite_core, trailing_vortices, xhat, include_wakes=false)
        system.dΓdt .= .-Γ
        circulation_derivatives!(Γ, system.dΓ, AIC, system.w, system.dw)
    else
        normal_velocity!(system.w, system.surfaces, wakes, ref, fs;
            additional_velocity=nothing, Vcp=system.Vcp, symmetric, nwake,
            surface_id, wake_finite_core, trailing_vortices, xhat, include_wakes=false)
        system.dΓdt .= .-Γ
        ldiv!(Γ, system.fAIC[], system.w)
    end
    system.dΓdt .+= Γ
    system.dΓdt ./= dt

    vehicle_on_all!(system, wake, trailing_edge_filaments; fmm_vehicle_args...)

    #------- near-field forces + viscous -------#

    if derivatives
        near_field_forces_derivatives!(system.properties, system.dproperties,
            system.surfaces, wakes, ref, fs, Γ, system.dΓ; dΓdt=system.dΓdt,
            additional_velocity=nothing, Vh=system.Vh, Vv=system.Vv, symmetric, nwake,
            surface_id, wake_finite_core, wake_shedding_locations,
            trailing_vortices, xhat, calculate_vlm_induced=false)
    else
        near_field_forces!(system.properties, system.surfaces, wakes,
            ref, fs, Γ; dΓdt=system.dΓdt, additional_velocity=nothing,
            Vh=system.Vh, Vv=system.Vv, symmetric, nwake, surface_id, wake_finite_core,
            wake_shedding_locations, trailing_vortices, xhat, calculate_vlm_induced=false)
    end

    Γ_wake .= Γ
    viscous!(system.properties, Γ_wake, system.surfaces, system.grids, frames,
        frames_index, polars, ref, dt)
    dΓdt_wake .+= Γ_wake
    dΓdt_wake ./= dt

    #------- monitors -------#

    system.near_field_analysis[] = true
    for monitor in monitors
        monitor(system, wake, i_step)
    end

    #------- save state -------#

    if !isnothing(path)
        write_vtk(joinpath(path, name * "_step_$i_step"), system; vtk_args...)
        FLOWVPM.save(wake, name * "_wake"; add_num=true, num=i_step, path, overwrite_time=i_step)
    end

    #------- propagate system -------#

    if i_step < length(t_range)
        FLOWVPM._euler(wake, dt; relax=true)
        update_vpm_shedding_TE!(wakes, ref, fs, dt, nothing, nothing)
        store_trailing_edge!(wake_shedding_locations, system.surfaces)
        propagate_kinematics!(system, frames, dt)

        idx = i_step == length(t_range) - 1 ? i_step + 1 : i_step + 2
        fs_next = velocity_to_freestream(Vinf(t_range[idx]), Ωinf(t_range[idx]))
        system.freestream[] = fs_next

        update_wake_shedding_locations_unsteady!(wakes, wake_shedding_locations,
            system.surfaces, ref, fs_next, dt, nothing, nothing, nwake, eta)

        shed_wake!(wake, system, dt, Γ_wake, system.dΓdt,
            particle_trailing_methods, particle_unsteady_methods)
        dΓdt_wake .= -Γ_wake

        if !isnothing(path) && vtk_postshed
            write_vtk(joinpath(path, name * "_postshed_step_$i_step"), system; write_wakes=true, vtk_args...)
            FLOWVPM.save(wake, name * "_postshed_wake"; add_num=true, num=i_step, path, overwrite_time=i_step)
        end
    end
end

"""
    simulate!(system::System, wake::PanelParticleWake, frames, maneuver!, Vinf, t_range,
              Ωinf=(t)->SVector{3}(0,0,0); kwargs...)

Drive an unsteady simulation against a buffer-overflow `PanelParticleWake`.
"""
function simulate!(system::System, wake::PanelParticleWake,
        frames::AbstractVector{<:ReferenceFrame}, maneuver!::Function,
        Vinf::Function, t_range, Ωinf=(t)->SVector{3}(0.0, 0.0, 0.0);
        name="vortex_lattice_simulation", path="./vortex_lattice_simulation",
        fmm_wake_args=(), fmm_vehicle_args=(),
        derivatives=false,
        monitors=(),
        recalculate_influence_matrix=true,
        polars=nothing, frames_index=fill(-1, length(system.surfaces)),
        restart_from::Union{Nothing,String}=nothing,
        restart_idx::Union{Nothing,Int}=nothing,
        vtk_interval::Int=1,
        write_restart::Bool=true,
        verbose=true,
    )
    (;restart_state, body_name, checkpoint_base, body_writer, wake_writer,
      write_restart, vtk_interval) = _init_simulate!(system, wake, frames, name, path,
        restart_from, restart_idx, write_restart, vtk_interval, t_range)

    ref               = system.reference[]
    Γ_wake            = zeros(length(system.Γ))
    dΓdt_wake         = zeros(length(system.Γ))
    trailing_edge_filaments = WakeBufferRings(wake)
    symmetric         = system.symmetric
    surface_id        = system.surface_id
    wake_finite_core  = system.wake_finite_core
    xhat              = system.xhat[]
    symmetric        .= false

    if isnothing(restart_state)
        _warm_start!(system, wake, frames, maneuver!, Vinf, Ωinf, t_range,
            trailing_edge_filaments, ref, symmetric, surface_id, wake_finite_core,
            xhat, derivatives, fmm_wake_args)
    end

    # bundle constant simulation parameters for the per-step function
    config = (;ref, symmetric, surface_id, wake_finite_core, xhat,
               recalculate_influence_matrix, derivatives, fmm_wake_args, fmm_vehicle_args,
               polars, frames_index, monitors, vtk_interval, body_name, write_restart,
               checkpoint_base, body_writer, wake_writer, Vinf, Ωinf, maneuver!, verbose,
               Γ_wake, dΓdt_wake)

    println()
    start_step = isnothing(restart_state) ? 0 : restart_state.idx
    for (it, t) in enumerate(t_range)
        i_step = it - 1
        i_step < start_step && continue
        _simulate_step!(system, wake, frames, trailing_edge_filaments, config, i_step, start_step, t, t_range)
    end

    if !isnothing(body_writer)
        _save_system_vtk_writer!(body_writer)
        _save_wake_vtk_writer!(wake_writer)
    end

    return wake
end

function simulate!(system::System, wake::ParticleField, frames::AbstractVector{<:ReferenceFrame}, maneuver!::Function, Vinf::Function, t_range, Ωinf=(t)->SVector{3}(0.0, 0.0, 0.0);
        name="vortex_lattice_simulation", path="./vortex_lattice_simulation",
        vtk_args=(trailing_vortices=false, write_wakes=false), vtk_postshed=false,
        fmm_wake_args=(), fmm_vehicle_args=(),
        derivatives=false,
        eta=0.3, 
        particle_trailing_methods=fill(OverlapPPS(1.3, 2), length(system.surfaces)),
        particle_unsteady_methods=fill(OverlapPPS(1.3, 2), length(system.surfaces)),
        trailing_vortices=fill(false, length(system.surfaces)),
        shedding_surfaces=fill(true, length(system.surfaces)),
        monitors=(),
        recalculate_influence_matrix=true,
        polars=nothing, frames_index=fill(-1, length(system.surfaces)), # viscous correction
        verbose=true
    )
    if !isnothing(path) && !isdir(path)
        mkpath(path)
    end

    # one row of wake panels per surface
    system.nwake .= 1
    for isurf in 1:length(system.surfaces)
        if shedding_surfaces[isurf]
            TF = eltype(system.Γ)
            panels = Matrix{WakePanel{TF}}(undef, 1, size(system.surfaces[isurf], 2))
            z3 = SVector{3,TF}(0, 0, 0)
            for i in 1:size(system.surfaces[isurf], 2)
                panels[1, i] = WakePanel{TF}(z3, z3, z3, z3, zero(TF), zero(TF))
            end
            system.wakes[isurf] = panels
        end
    end

    trailing_edge_filaments = FilamentWrapper(system.wakes)

    for isurf in 1:length(system.surfaces)
        nc, ns = size(system.wakes[isurf])
        system.V[isurf] = zeros(SVector{3,eltype(system.Γ)}, nc+1, ns+1)
    end

    system.trailing_vortices .= trailing_vortices
    system.freestream[] = velocity_to_freestream(Vinf(t_range[1]), Ωinf(t_range[1]))

    Γ_wake    = zeros(length(system.Γ))
    dΓdt_wake = zeros(length(system.Γ))

    config = (;eta, trailing_vortices, derivatives, recalculate_influence_matrix,
               fmm_wake_args, fmm_vehicle_args, particle_trailing_methods,
               particle_unsteady_methods, polars, frames_index, monitors,
               path, name, vtk_args, vtk_postshed, Vinf, Ωinf, maneuver!, verbose,
               Γ_wake, dΓdt_wake)

    println()
    i_step = 0
    for t in t_range
        _simulate_step!(system, wake, frames, trailing_edge_filaments, config, i_step, t, t_range)
        i_step += 1
    end
end


function update_trailing_edge_filaments!(trailing_edge_filaments::FilamentWrapper, current_surfaces, Γ::Vector{TF}) where TF
    # loop over surfaces
    wakes = trailing_edge_filaments.wakes
    iΓ = 0
    for isurf = eachindex(current_surfaces)
        surface = current_surfaces[isurf]
        wake = wakes[isurf]
        nc, ns = size(surface)
        for j in 1:ns
            # strength
            iΓ += nc

            # core size
            core_size = surface[nc, j].core_size

            # update wake panel circulation
            wp = wake[1, j]
            wake[1, j] = WakePanel{TF}(wp.rtl, wp.rtr, wp.rbl, wp.rbr, core_size, Γ[iΓ]) # use previous timestep's circulation
        end
    end
end

function wake_on_all!(system::System, wake::ParticleField, trailing_edge_filaments::FilamentWrapper; fmm_wake_args...)
    # reset probes
    FastMultipole.reset!(system.probes)

    # update probe positions
    update_probes!(system)
    
    # solve n-body problem
    fmm!((wake, system.probes), (wake, trailing_edge_filaments); hessian=SVector{2}(true,false), fmm_wake_args...) # solve N-body problem
    probes_to_surfaces!(system) # update Vcp, Vv, Vh, Vte, V based on probes
end

#------- wake shedding -------#

# WakeSheddingMethod, NoShed, SigmaPPS, OverlapPPS, SigmaOverlap are defined in wake.jl.

function shed_wake!(pfield::FLOWVPM.ParticleField, system, dt, Γ, dΓdt,
        shedding_trailing::AbstractVector{<:WakeSheddingMethod}, shedding_unsteady::AbstractVector{<:WakeSheddingMethod})
    # shed trailing edge particles
    shed_trailing_edge!(pfield, system.surfaces, system.wakes, Γ, shedding_trailing)

    # shed unsteady particles
    shed_unsteady!(pfield, system.surfaces, system.wakes, dΓdt, dt, shedding_unsteady)
end

function shed_trailing_edge!(pfield::FLOWVPM.ParticleField, surfaces, wakes, Γ, shedding_methods)
    # loop over surfaces
    iΓ = 0
    for isurf = eachindex(surfaces)
        surface = surfaces[isurf]
        wake = wakes[isurf]
        method = shedding_methods[isurf]
        nc, ns = size(surface)
        Γlast = zero(eltype(Γ))
        for j in 1:ns
            # strength
            iΓ += nc

            # get vertices
            panel = wake[1, j]
            r2 = top_left(panel)
            r1 = bottom_left(panel)

            # shed left particles
            Γthis = Γ[iΓ]
            shed_particles!(pfield, r1, r2, Γthis - Γlast, method)

            # recurse
            Γlast = Γthis
        end

        # get vertices
        panel = wake[1, end]
        r1 = top_right(panel)
        r2 = bottom_right(panel)

        # shed right particles
        shed_particles!(pfield, r1, r2, Γlast, method)
    end
end

function shed_unsteady!(pfield::FLOWVPM.ParticleField, surfaces, wakes, dΓdt, dt, shedding_methods)
    # loop over surfaces
    iΓ = 0
    for isurf = eachindex(surfaces)
        surface = surfaces[isurf]
        wake = wakes[isurf]
        method = shedding_methods[isurf]
        nc, ns = size(surface)
        for j in 1:ns
            # strength
            iΓ += nc

            # get vertices
            panel = wake[1, j]
            r2 = bottom_left(panel)
            r1 = bottom_right(panel)
            Γ = dΓdt[iΓ] * dt

            # shed unsteady particles
            shed_particles!(pfield, r1, r2, Γ, method)
        end
    end
end

function shed_particles!(pfield, r1, r2, Γ, method::OverlapPPS)
    # shed particles with overlap and p_per_step
    overlap = method.overlap
    p_per_step = method.p_per_step
    sigma = norm(r2 - r1) * overlap / p_per_step
    return shed_particles!(pfield, r1, r2, Γ, SigmaPPS(sigma, p_per_step))
end

function shed_particles!(pfield, r1, r2, Γ, method::SigmaOverlap)
    # shed particles with sigma and overlap
    sigma = method.sigma
    overlap = method.overlap
    p_per_step = ceil(Int, overlap * norm(r2 - r1) / sigma)
    return shed_particles!(pfield, r1, r2, Γ, SigmaPPS(sigma, p_per_step))
end

function shed_particles!(pfield, r1, r2, Γ, method::SigmaPPS)
    # shed particles with sigma and p_per_step
    sigma = method.sigma
    p_per_step = method.p_per_step

    # add particles
    distance_vector = (r2 - r1) / p_per_step
    Xp = r1 + distance_vector * 0.5
    Γp = Γ * distance_vector
end

function shed_particles!(pfield, r1, r2, Γ, method::NoShed)
    # do not shed particles
    return nothing
end

function get_max_particles(surface::AbstractMatrix{<:SurfacePanel}, method::Union{<:SigmaPPS, <:OverlapPPS})
    pps = method.p_per_step
    _, ns = size(surface)
    np = (ns + 1) * pps # pps particles at each trailing edge vertex
    return np
end

function get_max_particles(surface::AbstractMatrix{<:SurfacePanel}, method::SigmaOverlap)
    # estimate p_per_step
    pps = 8

    return get_max_particles(surface, SigmaPPS(method.sigma, pps))
end
    
function get_max_particles(surface::AbstractMatrix{<:SurfacePanel}, method::NoShed)
    # no particles shed
    return 0
end

function get_max_particles(system::System, particle_trailing_methods)
    np = 0
    for (isurf, surface) in enumerate(system.surfaces)
        np += get_max_particles(surface, particle_trailing_methods[isurf])
    end
    return np
end

function get_max_particles(system, particle_trailing_methods, particle_unsteady_methods)
    np_trailing = get_max_particles(system, particle_trailing_methods)
    np_unsteady = get_max_particles(system, particle_unsteady_methods)
    return np_trailing + np_unsteady
end

function vehicle_on_all!(system::System, wake::ParticleField, trailing_edge_filaments::FilamentWrapper; fmm_vehicle_args...)
    # reset probes
    FastMultipole.reset!(system.probes)

    # n-body problem
    fmm!((wake, system.probes), (system, ); hessian=SVector{2}(true,false), fmm_vehicle_args...)

    # update Vcp, Vv, Vh, and Vte based on probes
    probes_to_surfaces!(system)
end
