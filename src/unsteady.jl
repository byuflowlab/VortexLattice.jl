# WAKE_RHS_LOG (diagnostic, 2026-08-14): compares the wake's influence as seen by the
# Gamma-solve's boundary condition (Vcp, populated by wake_on_all! BEFORE the AIC solve,
# since normal_velocity! is called with include_wakes=false) against the wake's influence
# as seen by near_field_forces!/viscous.jl's induced-velocity readback (Vh, populated by
# vehicle_on_all! AFTER the solve) -- same wake state within the same timestep, two
# independent computations. If they disagree, the Gamma-solve's boundary condition is
# getting a different (and possibly wrong) picture of the wake than the verified-correct
# force/velocity readback uses, which would explain a self-consistent-but-wrong Gamma.
const WAKE_RHS_LOG_ENABLED = Ref(false)
const WAKE_RHS_LOG = Vector{NTuple{6,Float64}}()  # (Vcp_x,y,z, Vh_x,y,z) per logged step, post-solve (wake+kinematic)
const WAKE_RHS_KINEMATIC_LOG = Vector{NTuple{6,Float64}}()  # same, but right after kinematic_velocity! (pre-wake, kinematic only)
const WAKE_RHS_POSTWAKE_LOG = Vector{NTuple{6,Float64}}()  # same, but right after wake_on_all! (post-wake, pre-solve, pre-vehicle_on_all!) -- isolates whether vehicle_on_all! (which runs later, post-solve) modifies Vh at all

# VH_DECOMP_LOG (diagnostic, 2026-08-17): full-span, strip-frame-rotated (vx,vz) at blade 1
# (surface 1), captured at the same three _simulate_step! stages as WAKE_RHS_LOG above --
# right after kinematic_velocity!, right after wake_on_all!, and right after
# vehicle_on_all! -- to test whether the wake_on_all! pass (particles + wake panels/
# filaments) or the vehicle_on_all! pass (bound blade self-induction) is the source of the
# tangential/swirl over-prediction found in vortexlattice_axial_tangential_asymmetry_found_
# 2026_08_14. Rotation matches viscous.jl's `R = transpose(frame.Rp2g * frame.R)` convention
# exactly, and freestream_velocity(fs) is added so these are directly comparable to
# Velocity.txt's (vx,vz). Differencing WAKE_LOG-KINEMATIC_LOG isolates wake_on_all!'s
# contribution; VEHICLE_LOG-WAKE_LOG isolates vehicle_on_all!'s. Enable with
# VH_DECOMP_LOG_ENABLED[] = true; read out VH_DECOMP_KINEMATIC_LOG/_WAKE_LOG/_VEHICLE_LOG.
const VH_DECOMP_LOG_ENABLED = Ref(false)
const VH_DECOMP_KINEMATIC_LOG = Vector{Vector{Tuple{Float64,Float64}}}()
const VH_DECOMP_WAKE_LOG      = Vector{Vector{Tuple{Float64,Float64}}}()
const VH_DECOMP_VEHICLE_LOG   = Vector{Vector{Tuple{Float64,Float64}}}()

# RHS_WAKE_DUMP (diagnostic, 2026-08-17): ground-truth cross-check of the wake's
# contribution to the Gamma-solve's RHS. Dumps the exact wake state (particles +
# active buffer-ring panels, plus the boundary-filament body count) and the Vcp
# value at one control point immediately after wake_on_all! -- i.e. exactly what
# normal_velocity! will read as the wake's contribution to the RHS (Vcp is not
# touched again before the solve). An independent direct-sum Biot-Savart
# evaluation over this dumped state can then be diffed against VL's own Vcp,
# isolating the FMM/probe-writeback pipeline that feeds the RHS specifically --
# previously only the post-solve Vh/force-readback path was checked this way
# (see vortexlattice_probe_indexing_ruled_out_2026_08_14).
const RHS_WAKE_DUMP_ENABLED = Ref(false)
const RHS_WAKE_DUMP_STEP = Ref(-1)
const RHS_WAKE_DUMP_DIR = Ref("")

function _dump_rhs_wake_state(system, wake, i_step)
    isurf = 1
    ns = size(system.Vcp[isurf], 2)
    j = ns ÷ 2  # matches WAKE_RHS_LOG's jmid convention, for cross-referencing
    i = 1
    rcp = controlpoint(system.surfaces[isurf][i, j])
    vcp = system.Vcp[isurf][i, j]

    dir = RHS_WAKE_DUMP_DIR[]
    mkpath(dir)

    open(joinpath(dir, "cp.csv"), "w") do io
        println(io, "x,y,z,vcp_x,vcp_y,vcp_z")
        println(io, "$(rcp[1]),$(rcp[2]),$(rcp[3]),$(vcp[1]),$(vcp[2]),$(vcp[3])")
    end

    np = FLOWVPM.get_np(wake.pfield)
    open(joinpath(dir, "particles.csv"), "w") do io
        println(io, "x,y,z,gx,gy,gz,sigma")
        for p in 1:np
            X = FLOWVPM.get_X(wake.pfield, p)
            G = FLOWVPM.get_Gamma(wake.pfield, p)
            s = FLOWVPM.get_sigma(wake.pfield, p)[]
            println(io, "$(X[1]),$(X[2]),$(X[3]),$(G[1]),$(G[2]),$(G[3]),$s")
        end
    end

    open(joinpath(dir, "rings.csv"), "w") do io
        println(io, "rtl_x,rtl_y,rtl_z,rtr_x,rtr_y,rtr_z,rbl_x,rbl_y,rbl_z,rbr_x,rbr_y,rbr_z,core_size,gamma")
        for k in eachindex(wake.wakes)
            wk = wake.wakes[k]
            nwk = wake.nwake[k]
            nc, nsk = size(wk)
            for jj in 1:nsk, ii in 1:min(nwk, nc)
                p = wk[ii, jj]
                println(io, "$(p.rtl[1]),$(p.rtl[2]),$(p.rtl[3]),$(p.rtr[1]),$(p.rtr[2]),$(p.rtr[3]),$(p.rbl[1]),$(p.rbl[2]),$(p.rbl[3]),$(p.rbr[1]),$(p.rbr[2]),$(p.rbr[3]),$(p.core_size),$(p.gamma)")
            end
        end
    end

    bfw = wake.boundary_filaments
    open(joinpath(dir, "boundary_filaments.csv"), "w") do io
        println(io, "r1_x,r1_y,r1_z,r2_x,r2_y,r2_z,core_size,gamma")
        for isurf in eachindex(bfw.active)
            bfw.active[isurf] || continue
            for j in eachindex(bfw.r1[isurf])
                r1 = bfw.r1[isurf][j]
                r2 = bfw.r2[isurf][j]
                println(io, "$(r1[1]),$(r1[2]),$(r1[3]),$(r2[1]),$(r2[2]),$(r2[3]),$(bfw.core_size[isurf][j]),$(bfw.gamma[isurf][j])")
            end
        end
    end

    nbf = FastMultipole.get_n_bodies(wake.boundary_filaments)

    open(joinpath(dir, "meta.txt"), "w") do io
        println(io, "i_step=$i_step np=$np nbf=$nbf")
    end

    println("RHS_WAKE_DUMP written to $dir at step $i_step (np=$np, nbf=$nbf)")
end

# AIC_CHECK_DUMP (diagnostic, 2026-08-17): dump surface-1's bound-panel geometry,
# wake_shedding_locations, trailing_vortices/symmetric flags, solved Gamma, and the
# surface-1 diagonal block of the production AIC matrix -- so a standalone script can
# independently rebuild that AIC block via VortexLattice's OWN influence_coefficients!,
# both through the production (finite_core=false, edge-reuse/bookkeeping) branch and
# the simple (finite_core=true, no-bookkeeping, per-column) branch, and check that the
# two code paths agree off the near-diagonal (where finite-core specifics shouldn't
# matter at ~1.5 m core vs ~50 m span) -- isolating a possible indexing/wiring bug in
# the edge-reuse bookkeeping from the (already independently verified, see
# RHS_WAKE_DUMP) Biot-Savart kernel itself.
function _dump_aic_check_state(system, config, i_step)
    isurf = 1
    surface = system.surfaces[isurf]
    nc, ns = size(surface)
    dir = RHS_WAKE_DUMP_DIR[]
    mkpath(dir)

    open(joinpath(dir, "bound_panels.csv"), "w") do io
        println(io, "rtl_x,rtl_y,rtl_z,rtr_x,rtr_y,rtr_z,rbl_x,rbl_y,rbl_z,rbr_x,rbr_y,rbr_z,rcp_x,rcp_y,rcp_z,ncp_x,ncp_y,ncp_z,core_size")
        for j in 1:ns, i in 1:nc
            p = surface[i, j]
            println(io, "$(p.rtl[1]),$(p.rtl[2]),$(p.rtl[3]),$(p.rtr[1]),$(p.rtr[2]),$(p.rtr[3]),$(p.rbl[1]),$(p.rbl[2]),$(p.rbl[3]),$(p.rbr[1]),$(p.rbr[2]),$(p.rbr[3]),$(p.rcp[1]),$(p.rcp[2]),$(p.rcp[3]),$(p.ncp[1]),$(p.ncp[2]),$(p.ncp[3]),$(p.core_size)")
        end
    end

    wsl = system.wake_shedding_locations[isurf]
    open(joinpath(dir, "wake_shedding_locations.csv"), "w") do io
        println(io, "x,y,z")
        for r in wsl
            println(io, "$(r[1]),$(r[2]),$(r[3])")
        end
    end

    n1 = nc * ns
    Gamma1 = system.Γ[1:n1]
    open(joinpath(dir, "gamma_surf1.csv"), "w") do io
        println(io, "gamma")
        for g in Gamma1
            println(io, "$g")
        end
    end

    AIC1 = system.AIC[1:n1, 1:n1]
    open(joinpath(dir, "aic_surf1_block.csv"), "w") do io
        for i in 1:n1
            println(io, join(AIC1[i, :], ","))
        end
    end

    open(joinpath(dir, "aic_meta.txt"), "w") do io
        println(io, "nc=$nc ns=$ns symmetric=$(config.symmetric[isurf]) trailing_vortices=$(system.trailing_vortices[isurf]) xhat=$(config.xhat)")
    end

    println("AIC_CHECK_DUMP written to $dir at step $i_step (nc=$nc, ns=$ns)")
end

function _vh_decomp_row(system, frames, frames_index, fs)
    isurf = 1
    frame = frames[frames_index[isurf]]
    R = transpose(frame.Rp2g * frame.R)
    Vfs = freestream_velocity(fs)
    Vh = system.Vh[isurf]
    ns = size(Vh, 2)
    row = Vector{Tuple{Float64,Float64}}(undef, ns)
    for j in 1:ns
        v = R * (Vfs + Vh[1, j])
        row[j] = (v[1], v[3])
    end
    return row
end

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
    @assert wake.nwakerows >= 0 "PanelParticleWake requires nwakerows >= 0, got $(wake.nwakerows)"

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
        xhat, derivatives, fmm_wake_args, force_finite_core)
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
        force_finite_core)
    update_trailing_edge_coefficients!(system.AIC, system.surfaces;
        symmetric, wake_shedding_locations=system.wake_shedding_locations,
        trailing_vortices=system.trailing_vortices, force_finite_core)
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
        # (nwakerows==0: only the boundary filament is seeded — no particles are
        # emitted here, so the impulsive-start transient isn't double-counted)
        shed_wake!(wake, system, dt0, system.Γ; emit_particles=false)
        system.nwake .= wake.nwake
        reset!(wake)
        update_trailing_edge_filaments!(trailing_edge_filaments, system.surfaces, system.Γ)
        wake.pfield.SFS(wake.pfield, FLOWVPM.BeforeUJ())
        wake_on_all!(system, wake, trailing_edge_filaments; fmm_wake_args...)
    end
end

function _simulate_step!(system, wake::PanelParticleWake, frames, trailing_edge_filaments,
        config, i_step, start_step, t, t_range)
    (;ref, symmetric, surface_id, wake_finite_core, xhat, force_finite_core,
      recalculate_influence_matrix, derivatives, fmm_wake_args, fmm_vehicle_args,
      polars, frames_index, monitors, vtk_interval, body_name, write_restart,
      checkpoint_base, body_writer, wake_writer, Vinf, Ωinf, maneuver!, verbose,
      Γ_wake, dΓdt_wake, viscous_iterative_shed, viscous_frames) = config

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

    if WAKE_RHS_LOG_ENABLED[]
        jmid_diag2 = size(system.Vcp[1], 2) ÷ 2
        vcp0 = system.Vcp[1][1, jmid_diag2]
        vh0 = system.Vh[1][1, jmid_diag2]
        push!(WAKE_RHS_KINEMATIC_LOG, (vcp0[1], vcp0[2], vcp0[3], vh0[1], vh0[2], vh0[3]))
    end

    if VH_DECOMP_LOG_ENABLED[]
        push!(VH_DECOMP_KINEMATIC_LOG, _vh_decomp_row(system, frames, frames_index, fs))
    end

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

    if RHS_WAKE_DUMP_ENABLED[] && i_step == RHS_WAKE_DUMP_STEP[]
        _dump_rhs_wake_state(system, wake, i_step)
    end

    if WAKE_RHS_LOG_ENABLED[]
        jmid_diag3 = size(system.Vcp[1], 2) ÷ 2
        vcp1 = system.Vcp[1][1, jmid_diag3]
        vh1 = system.Vh[1][1, jmid_diag3]
        push!(WAKE_RHS_POSTWAKE_LOG, (vcp1[1], vcp1[2], vcp1[3], vh1[1], vh1[2], vh1[3]))
    end

    if VH_DECOMP_LOG_ENABLED[]
        push!(VH_DECOMP_WAKE_LOG, _vh_decomp_row(system, frames, frames_index, fs))
    end

    #------- AIC + solve -------#

    AIC = system.AIC
    Γ   = system.Γ
    trailing_vortices = system.trailing_vortices

    if recalculate_influence_matrix || i_step == start_step
        influence_coefficients!(AIC, system.surfaces;
            symmetric, wake_shedding_locations=system.wake_shedding_locations,
            surface_id, trailing_vortices, xhat,
            force_finite_core)
        update_trailing_edge_coefficients!(AIC, system.surfaces;
            symmetric, wake_shedding_locations=system.wake_shedding_locations,
            trailing_vortices, force_finite_core)
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

    if RHS_WAKE_DUMP_ENABLED[] && i_step == RHS_WAKE_DUMP_STEP[]
        _dump_aic_check_state(system, config, i_step)
    end

    if WAKE_RHS_LOG_ENABLED[]
        jmid_diag = size(system.Vcp[1], 2) ÷ 2
        vcp = system.Vcp[1][1, jmid_diag]
        vh = system.Vh[1][1, jmid_diag]
        push!(WAKE_RHS_LOG, (vcp[1], vcp[2], vcp[3], vh[1], vh[2], vh[3]))
    end

    if VH_DECOMP_LOG_ENABLED[]
        push!(VH_DECOMP_VEHICLE_LOG, _vh_decomp_row(system, frames, frames_index, fs))
    end

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
        if viscous_iterative_shed
            # fully-converged (Anderson) per-step viscous coupling, fed to shed_wake! below,
            # instead of the single-pass additive viscous!() correction (which diverges when
            # applied every step -- see vortexlattice_viscous_correction_diverges_confirmed_2026_08_17).
            # Reuses the SAME FMM machinery (wake_on_all!/vehicle_on_all!) this step already
            # computed, re-run against a trial Gamma each iteration, so the shed wake's own
            # strength tracks the viscous-corrected circulation instead of staying frozen at the
            # inviscid history (see vortexlattice_viscous_iterative_frozen_wake_selfconsistency_ROOT_CAUSE_20260820).
            viscous_iterative_shed!(system, wake, trailing_edge_filaments, Γ_wake, polars,
                ref, fs, frames_index)
        else
            # NOTE: intentionally uses `viscous_frames` (a static, non-rotated axis reference),
            # NOT the kinematic `frames` this step otherwise uses for Vh/Vv/motion bookkeeping.
            # `viscous!()` derives its local (chordwise, normal) decomposition from
            # transpose(frame.Rp2g * frame.R); the kinematic vehicle frame's R (chosen to orient
            # Uinf(t)/Ωinf(t) correctly) does not represent the grid's own body-axis orientation
            # and silently feeds viscous!() the wrong rotation (confirmed: nets to a Z-flip
            # instead of the X-flip the validated steady single-pass usage relies on, corrupting
            # every station's local alpha and making the per-step correction have ~zero net
            # effect on the reported CL/CD regardless of wake resolution -- see
            # vortexlattice_viscous_frame_mismatch_ROOT_CAUSE_FOUND_20260820 memory).
            viscous!(system.properties, Γ_wake, system.surfaces, system.grids, viscous_frames,
                frames_index, polars, ref, dt)
        end
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
    (;eta, trailing_vortices, derivatives, recalculate_influence_matrix, force_finite_core,
      fmm_wake_args, fmm_vehicle_args, particle_trailing_methods, particle_unsteady_methods,
      polars, frames_index, monitors, path, name, vtk_args, vtk_postshed,
      Vinf, Ωinf, maneuver!, verbose, Γ_wake, dΓdt_wake, viscous_frames) = config

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
            force_finite_core)
        update_trailing_edge_coefficients!(AIC, system.surfaces;
            symmetric, wake_shedding_locations, trailing_vortices, force_finite_core)
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
    if !isnothing(polars)
        # see the NOTE at the PanelParticleWake call site above: must use the static
        # `viscous_frames`, not the kinematic `frames`.
        viscous!(system.properties, Γ_wake, system.surfaces, system.grids, viscous_frames,
            frames_index, polars, ref, dt)
    end
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
        viscous_iterative_shed::Bool=false,
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
    # real near/far wake panels + shed particles are the sole representation of trailing-vortex
    # downwash for this wake type -- leaving system.trailing_vortices at its seeding-call default
    # (true) double-counts that downwash (once via the AIC's semi-infinite filament assumption,
    # again via the real wake), which compounds every step into runaway divergence (confirmed:
    # rectangular NACA4415 AR=12 wing, aoa=4.2deg, inviscid, goes from steady CL=0.371 to
    # unsteady CL=9453 at 8 chords without this, converging to ~0.367 at 16 chords with it).
    # steady_analysis! avoids this by masking trailing_vortices off whenever wake panels are
    # present (analyses.jl); do the same thing here for the PanelParticleWake stepping loop.
    system.trailing_vortices .= false
    force_finite_core = system.surface_finite_core

    if isnothing(restart_state)
        _warm_start!(system, wake, frames, maneuver!, Vinf, Ωinf, t_range,
            trailing_edge_filaments, ref, symmetric, surface_id, wake_finite_core,
            xhat, derivatives, fmm_wake_args, force_finite_core)
    end

    # static, non-rotated axis reference for viscous!()'s local (chordwise, normal) velocity
    # decomposition -- see the NOTE at its call site in _simulate_step! for why this must NOT be
    # the kinematic `frames` used for Vh/Vv/motion elsewhere in the step.
    viscous_frames = ReferenceFrame(system;
        origin=SVector{3}(0.0, 0.0, 0.0), v=SVector{3}(0.0, 0.0, 0.0),
        ω_axis=SVector{3}(1.0, 0.0, 0.0), ω=0.0,
        R=SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="viscous_static", child_index=Int[], dependent_index=collect(1:length(system.surfaces)))

    # bundle constant simulation parameters for the per-step function
    config = (;ref, symmetric, surface_id, wake_finite_core, xhat, force_finite_core,
               recalculate_influence_matrix, derivatives, fmm_wake_args, fmm_vehicle_args,
               polars, frames_index, monitors, vtk_interval, body_name, write_restart,
               checkpoint_base, body_writer, wake_writer, Vinf, Ωinf, maneuver!, verbose,
               Γ_wake, dΓdt_wake, viscous_iterative_shed, viscous_frames)

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
    force_finite_core = system.surface_finite_core

    # static, non-rotated axis reference for viscous!() -- see the NOTE at its call site.
    viscous_frames = ReferenceFrame(system;
        origin=SVector{3}(0.0, 0.0, 0.0), v=SVector{3}(0.0, 0.0, 0.0),
        ω_axis=SVector{3}(1.0, 0.0, 0.0), ω=0.0,
        R=SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="viscous_static", child_index=Int[], dependent_index=collect(1:length(system.surfaces)))

    config = (;eta, trailing_vortices, derivatives, recalculate_influence_matrix, force_finite_core,
               fmm_wake_args, fmm_vehicle_args, particle_trailing_methods,
               particle_unsteady_methods, polars, frames_index, monitors,
               path, name, vtk_args, vtk_postshed, Vinf, Ωinf, maneuver!, verbose,
               Γ_wake, dΓdt_wake, viscous_frames)

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
            _shed_particles!(pfield, r1, r2, Γthis - Γlast, method)

            # recurse
            Γlast = Γthis
        end

        # get vertices
        panel = wake[1, end]
        r1 = top_right(panel)
        r2 = bottom_right(panel)

        # shed right particles
        _shed_particles!(pfield, r1, r2, Γlast, method)
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
            _shed_particles!(pfield, r1, r2, Γ, method)
        end
    end
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
