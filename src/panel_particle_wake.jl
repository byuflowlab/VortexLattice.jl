# ============================================================================
# PanelParticleWake helpers (Part 5 of the PanelParticleWake port).
#
# The `PanelParticleWake` struct itself and its core lifecycle methods
# (`reset!`, `update_TE!`, `propagate!`, `_convert_to_particles!`, `shed_wake!`)
# live in `wake.jl`. This file holds the surrounding simulation helpers:
# integration schemes and freestream application for now, with room for
# `_reset_system!`, `update_wake_shedding_locations_unsteady!`, etc. as the
# `simulate!` rewrite progresses.
# ============================================================================

"""
    IntegrationScheme

Abstract supertype for particle-field time integrators used by
`propagate!(::PanelParticleWake, dt; scheme)`. Concrete subtypes select the
underlying FLOWVPM stepper.
"""
abstract type IntegrationScheme end

"""
    EulerScheme()

First-order explicit Euler. Dispatches to `FLOWVPM._euler`.
"""
struct EulerScheme <: IntegrationScheme end

"""
    RK3Scheme()

Low-storage third-order Runge–Kutta. Stub: particle-only RK3 is a small
addition (dispatch to `FLOWVPM._rungekutta3`), while body-influence-aware RK3
requires restructuring `simulate!`'s FMM calls to re-evaluate at each substage.
Not yet wired.
"""
struct RK3Scheme <: IntegrationScheme end

"""
    _integrate_particles!(w::PanelParticleWake, dt, scheme; relax)

Internal dispatch target for `propagate!`. Selects the FLOWVPM stepper that
advances the particle field by `dt`.
"""
function _integrate_particles!(w::PanelParticleWake, dt, ::EulerScheme; relax)
    FLOWVPM._euler(w.pfield, dt; relax)
    return w
end

function _integrate_particles!(::PanelParticleWake, dt, ::RK3Scheme; relax=true)
    error("RK3Scheme not yet wired for PanelParticleWake (dt=$dt, relax=$relax). " *
          "Particle-only RK3 is a drop-in call to FLOWVPM._rungekutta3; " *
          "body-influence-aware RK3 requires restructuring simulate!'s FMM " *
          "calls per substage.")
end

"""
    PanelBufferFilaments{TF}

FMM-compatible source wrapper over the *active* portion of a
`PanelParticleWake`'s panel buffer. Unlike `FilamentWrapper` — which exposes
every slot in `Vector{Matrix{WakePanel}}` — this wrapper honors each surface's
`nwake[isurf]` count so that inactive (not-yet-shed) buffer rows are hidden
from FMM entirely. Inactive rows therefore contribute no influence and are
never read, so the buffer may be left uninitialized until `shed_wake!`
populates it.
"""
struct PanelBufferFilaments{TF}
    wakes::Vector{Matrix{WakePanel{TF}}}
    nwake::Vector{Int}
end

PanelBufferFilaments(w::PanelParticleWake) = PanelBufferFilaments(w.wakes, w.nwake)

function _fmm_kwargs(fmm::FLOWVPM.FMM, useGPU::Int)
    return (
        expansion_order = max(fmm.p - 1, 0),
        leaf_size_source = max(fmm.ncrit, fmm.min_ncrit),
        multipole_acceptance = fmm.theta,
        error_tolerance = FastMultipole.PowerRelativeGradient{fmm.relative_tolerance, fmm.absolute_tolerance, true}(),
        tune = true,
        nearfield_device = (useGPU > 0),
    )
end

function _update_fmm_autotune!(w::PanelParticleWake, which::Symbol, fmm_args)
    fmm_state = which === :wake ? w.fmm_wake : w.fmm_vehicle
    fmm = fmm_state[]
    optargs = fmm_args[1]

    new_p = fmm.autotune_p ? optargs.expansion_order + 1 : fmm.p
    if new_p < fmm.p
        new_p = fmm.p
    end

    new_ncrit = fmm.autotune_ncrit ? optargs.leaf_size_source[1] : fmm.ncrit

    new_fmm = FLOWVPM.FMM(
        p = new_p,
        ncrit = new_ncrit,
        theta = fmm.theta,
        shrink_recenter = fmm.shrink_recenter,
        relative_tolerance = fmm.relative_tolerance,
        absolute_tolerance = fmm.absolute_tolerance,
        autotune_p = fmm.autotune_p,
        autotune_ncrit = fmm.autotune_ncrit,
        autotune_reg_error = fmm.autotune_reg_error,
        default_rho_over_sigma = fmm.default_rho_over_sigma,
        min_ncrit = fmm.min_ncrit,
    )

    if which === :wake
        w.fmm_wake[] = new_fmm
    else
        w.fmm_vehicle[] = new_fmm
    end

    return nothing
end

Base.eltype(::PanelBufferFilaments{TF}) where TF = TF

function _active_to_matrix_index(pbf::PanelBufferFilaments, n)
    n_counter = 0
    for k in eachindex(pbf.wakes)
        nwk = pbf.nwake[k]
        nwk == 0 && continue
        ns = size(pbf.wakes[k], 2)
        block = nwk * ns
        if n_counter + block >= n
            i_plus = n - n_counter
            i = mod(i_plus - 1, nwk) + 1
            j = div(i_plus - 1, nwk) + 1
            return k, i, j
        end
        n_counter += block
    end
    error("PanelBufferFilaments index $n out of active range")
end

function FastMultipole.get_n_bodies(pbf::PanelBufferFilaments)
    n = 0
    for k in eachindex(pbf.wakes)
        n += pbf.nwake[k] * size(pbf.wakes[k], 2)
    end
    return n
end

FastMultipole.data_per_body(::PanelBufferFilaments)    = 12
FastMultipole.strength_dims(::PanelBufferFilaments)    = 1
FastMultipole.has_vector_potential(::PanelBufferFilaments) = true

function FastMultipole.get_position(pbf::PanelBufferFilaments, i)
    k, ir, jc = _active_to_matrix_index(pbf, i)
    panel = pbf.wakes[k][ir, jc]
    return 0.5 * (panel.rtl + panel.rbr)
end

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, pbf::PanelBufferFilaments, i_body)
    k, ir, jc = _active_to_matrix_index(pbf, i_body)
    panel = pbf.wakes[k][ir, jc]
    buffer[1:3, i_buffer] .= 0.5 * (panel.rtl + panel.rbr)
    buffer[4, i_buffer]    = 0.5 * norm(panel.rbl - panel.rtr) + panel.core_size
    buffer[5, i_buffer]    = panel.gamma
    buffer[6:8,  i_buffer] .= panel.rtl
    buffer[9:11, i_buffer] .= panel.rtr
    buffer[12, i_buffer]   = panel.core_size
end

function FastMultipole.body_to_multipole!(pbf::PanelBufferFilaments, multipole_coefficients, buffer::Matrix, center, bodies_index, harmonics, expansion_order)
    for i_body in bodies_index
        rtl = FastMultipole.get_vertex(buffer, pbf, i_body, 1)
        rtr = FastMultipole.get_vertex(buffer, pbf, i_body, 2)
        gamma = FastMultipole.get_strength(buffer, pbf, i_body)[1]
        body_to_multipole_vl!(multipole_coefficients, harmonics, rtl, rtr, center, gamma, expansion_order)
    end
end

function FastMultipole.direct!(target_system, target_index, ::DerivativesSwitch{PS,VS,GS}, source_system::PanelBufferFilaments{TF}, source_buffer, source_index) where {PS,VS,GS,TF}
    @inbounds for j_target in target_index
        target = FastMultipole.get_position(target_system, j_target)
        v = SVector{3,TF}(0.0, 0.0, 0.0)
        @inbounds for i_source in source_index
            v1 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 1)
            v2 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 2)
            gamma = FastMultipole.get_strength(source_buffer, source_system, i_source)[1]
            cs = source_buffer[12, i_source]
            if VS
                v += bound_induced_velocity(target - v1, target - v2, true, cs) * gamma
            end
        end
        FastMultipole.set_gradient!(target_system, j_target, v)
    end
end

function FastMultipole.buffer_to_target_system!(target_system::PanelBufferFilaments, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    @warn "A PanelBufferFilaments should not be used as a target in an FMM call."
end

"""
    update_trailing_edge_filaments!(pbf::PanelBufferFilaments, current_surfaces, Γ)

Populate row-1 strengths of the active panel-buffer filaments from the
current circulation vector `Γ`. Row-1 only — deeper buffer rows keep the
strength they were given when originally shed. Surfaces with `nwake==0` are
skipped (their Γ indices are still advanced).
"""
function update_trailing_edge_filaments!(pbf::PanelBufferFilaments, current_surfaces, Γ::Vector{TF}) where TF
    wakes = pbf.wakes
    iΓ = 0
    for isurf in eachindex(current_surfaces)
        surface = current_surfaces[isurf]
        nc, ns = size(surface)
        if pbf.nwake[isurf] == 0
            iΓ += nc * ns
            continue
        end
        wake = wakes[isurf]
        for j in 1:ns
            iΓ += nc
            core_size = surface[nc, j].core_size
            wp = wake[1, j]
            wake[1, j] = WakePanel{TF}(wp.rtl, wp.rtr, wp.rbl, wp.rbr, core_size, Γ[iΓ])
        end
    end
end

"""
    wake_on_all!(system, wake::PanelParticleWake, trailing_edge_filaments; fmm_wake_args...)

FMM wake-influence pass for `PanelParticleWake`. Sources are the overflow
particle field (`wake.pfield`) plus the active panel-buffer filaments;
targets are the particle field and the body probes.
"""
function wake_on_all!(system, wake::PanelParticleWake,
        trailing_edge_filaments::PanelBufferFilaments; fmm_wake_args...)
    FastMultipole.reset!(system.probes)
    update_probes!(system)
    np = FLOWVPM.get_np(wake.pfield)
    nfil = FastMultipole.get_n_bodies(trailing_edge_filaments)
    if np > 0 && nfil > 0
        fmm_args = fmm!((wake.pfield, system.probes), (wake.pfield, trailing_edge_filaments);
            _fmm_kwargs(wake.fmm_wake[], wake.pfield.useGPU)...,
            hessian=SVector{2}(true, false), fmm_wake_args...)
        _update_fmm_autotune!(wake, :wake, fmm_args)
    elseif np > 0
        fmm_args = fmm!((wake.pfield, system.probes), (wake.pfield,);
            _fmm_kwargs(wake.fmm_wake[], wake.pfield.useGPU)...,
            hessian=SVector{2}(true, false), fmm_wake_args...)
        _update_fmm_autotune!(wake, :wake, fmm_args)
    elseif nfil > 0
        for V in system.V
            fill!(V, zero(eltype(V)))
        end
    end
    probes_to_surfaces!(system)
    return wake
end

"""
    vehicle_on_all!(system, wake::PanelParticleWake, trailing_edge_filaments; fmm_vehicle_args...)

FMM vehicle-influence pass for `PanelParticleWake`. Sources are the body
surfaces; targets are the particle field and the body probes.
"""
function vehicle_on_all!(system, wake::PanelParticleWake,
        trailing_edge_filaments::PanelBufferFilaments; fmm_vehicle_args...)
    FastMultipole.reset!(system.probes)
    update_probes!(system)
    np = FLOWVPM.get_np(wake.pfield)
    if np > 0
        fmm_args = fmm!((wake.pfield, system.probes), (system,);
            _fmm_kwargs(wake.fmm_vehicle[], wake.pfield.useGPU)...,
            hessian=SVector{2}(true, false), fmm_vehicle_args...)
        _update_fmm_autotune!(wake, :vehicle, fmm_args)
    else
        fmm!((system.probes,), (system,);
            hessian=SVector{1}(false), fmm_vehicle_args...)
    end
    probes_to_surfaces!(system)
    return wake
end

"""
    apply_freestream!(w::PanelParticleWake, Vinf)

Seed the freestream velocity onto all active particles in `w.pfield` so that
the next integration step advects them with the external flow. `Vinf` may be
an `SVector{3}`, tuple, or any 3-indexable.
"""
function apply_freestream!(w::PanelParticleWake, Vinf)
    np = FLOWVPM.get_np(w.pfield)
    @inbounds for p in 1:np, d in 1:3
        w.pfield.particles[FLOWVPM.U_INDEX[d], p] += Vinf[d]
    end
    return w
end

function _seed_wake_velocity!(w::PanelParticleWake, fs)
    Vfs = freestream_velocity(fs)
    for isurf in eachindex(w.wakes)
        nwk = w.nwake[isurf]
        nwk == 0 && continue
        V = w.wake_velocities[isurf]
        ns = size(w.wakes[isurf], 2)
        for j in 1:ns+1, i in 1:nwk+1
            V[i, j] += Vfs
        end
    end
end

function _update_wsl!(wsl_all, surfaces, fs, dt, eta)
    Vfs = freestream_velocity(fs)
    for isurf in eachindex(surfaces)
        wsl = wsl_all[isurf]
        for j in eachindex(wsl)
            wsl[j] += eta * Vfs * dt
        end
    end
end
