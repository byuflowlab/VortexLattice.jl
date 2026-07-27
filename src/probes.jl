# --- FMM probe target system --- #

"""
    ProbeSystem(n_bodies, TF=Float64)

Target system holding `n_bodies` evaluation points along with the scalar
potential, gradient (velocity), and hessian accumulated onto them by an FMM
call.

This is `FastMultipole.ProbeSystem` plus a retained estimate of the influence
computed by the preceding FMM pass (`previous_potential`, `previous_gradient`),
exposed to FastMultipole through its previous-influence metadata interface.
FastMultipole tunes a relative error tolerance against that estimate
(`ε = max(branch.min_gradient*RET, AET)`); a target system that does not
provide it reports zero, and since a branch takes the minimum over all of its
target systems, a single non-participating system drags the whole branch back
to absolute-tolerance behavior — including for the particle field sharing that
branch. `FastMultipole.ProbeSystemStatic` has nowhere to store such an estimate,
so probes are given their own type here rather than reusing it.

Field names match `FastMultipole.ProbeSystemStatic` so the two are
interchangeable at every use site.
"""
struct ProbeSystem{TF}
    position::Vector{SVector{3,TF}}
    scalar_potential::Vector{TF}
    gradient::Vector{SVector{3,TF}}
    hessian::Vector{SMatrix{3,3,TF,9}}
    previous_potential::Vector{TF}
    previous_gradient::Vector{TF} # magnitude only -- all the metadata interface consumes
end

function ProbeSystem(n_bodies, TF=Float64)
    return ProbeSystem{TF}(
        zeros(SVector{3,TF}, n_bodies),
        zeros(TF, n_bodies),
        zeros(SVector{3,TF}, n_bodies),
        zeros(SMatrix{3,3,TF,9}, n_bodies),
        zeros(TF, n_bodies),
        zeros(TF, n_bodies),
    )
end

Base.eltype(::ProbeSystem{TF}) where TF = TF
FastMultipole.numtype(::ProbeSystem{TF}) where TF = TF

"""
    reset!(probes::ProbeSystem)

Zero the accumulated influence, first rolling it into `previous_potential` /
`previous_gradient` so the next FMM pass can use it for relative error control.
"""
function FastMultipole.reset!(probes::ProbeSystem{TF}) where TF
    for i in eachindex(probes.gradient)
        probes.previous_potential[i] = abs(probes.scalar_potential[i])
        probes.previous_gradient[i] = norm(probes.gradient[i])
        probes.scalar_potential[i] = zero(TF)
        probes.gradient[i] = zero(SVector{3,TF})
        probes.hessian[i] = zero(SMatrix{3,3,TF,9})
    end
end

"""
    seed_previous_influence!(dest::ProbeSystem, src::ProbeSystem, n)

Copy the first `n` retained previous-influence estimates from `src` to `dest`.
Used when a subset of a persistent probe set is packed into a temporary
`ProbeSystem` for an FMM call.
"""
function seed_previous_influence!(dest::ProbeSystem, src::ProbeSystem, n=length(dest.position))
    @views dest.previous_potential[1:n] .= src.previous_potential[1:n]
    @views dest.previous_gradient[1:n] .= src.previous_gradient[1:n]
    return dest
end

"""
    resize_active!(probes::ProbeSystem, n)

Resize `probes` to `n` active bodies, reusing its existing backing arrays
rather than allocating new ones (`resize!` down then back up within the same
capacity is a no-op allocation-wise). `n` must not exceed the capacity
`probes` was originally constructed with. Zeros the FMM-accumulated fields
(`scalar_potential`, `gradient`, `hessian`); `position` and the
previous-influence fields are left for the caller to fill.
"""
function resize_active!(probes::ProbeSystem{TF}, n) where TF
    resize!(probes.position, n)
    resize!(probes.scalar_potential, n)
    resize!(probes.gradient, n)
    resize!(probes.hessian, n)
    resize!(probes.previous_potential, n)
    resize!(probes.previous_gradient, n)
    fill!(probes.scalar_potential, zero(TF))
    fill!(probes.gradient, zero(SVector{3,TF}))
    fill!(probes.hessian, zero(SMatrix{3,3,TF,9}))
    return probes
end

#--- FastMultipole target interface ---#

FastMultipole.data_per_body(::ProbeSystem) = 3
FastMultipole.strength_dims(::ProbeSystem) = 0
FastMultipole.get_n_bodies(probes::ProbeSystem) = length(probes.position)
FastMultipole.get_position(probes::ProbeSystem, i) = probes.position[i]

FastMultipole.metadata_per_body(::ProbeSystem) = 2
FastMultipole.previous_potential_metadata_index(::ProbeSystem) = 1
FastMultipole.previous_gradient_metadata_index(::ProbeSystem) = 2

function FastMultipole.metadata_to_buffer!(buffer, switch, i_buffer, probes::ProbeSystem, i_body)
    buffer[FastMultipole.metadata_index(switch, 1), i_buffer] = probes.previous_potential[i_body]
    buffer[FastMultipole.metadata_index(switch, 2), i_buffer] = probes.previous_gradient[i_body]
end

function FastMultipole.buffer_to_target_system!(probes::ProbeSystem, i_target, switch::FastMultipole.DerivativesSwitch{PS,GS,HS}, target_buffer, i_buffer) where {PS,GS,HS}
    if PS
        probes.scalar_potential[i_target] += FastMultipole.get_scalar_potential(target_buffer, switch, i_buffer)
    end
    if GS
        probes.gradient[i_target] += FastMultipole.get_gradient(target_buffer, switch, i_buffer)
    end
    if HS
        probes.hessian[i_target] += FastMultipole.get_hessian(target_buffer, switch, i_buffer)
    end
end

#--- unused source interface (probes are targets only) ---#

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, probes::ProbeSystem, i_body)
    throw("a ProbeSystem cannot be used as a source system")
end

FastMultipole.body_to_multipole!(probes::ProbeSystem, args...) = nothing

FastMultipole.direct!(target_system, target_index, derivatives_switch, probes::ProbeSystem, source_buffer, source_index) = nothing
