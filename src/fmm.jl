#------- FastMultipole compatibility functions -------#

function vlm_to_fmm_index(system::System, i_surf, i, j)
    surfaces = system.surfaces
    n = 0
    for k in 1:i_surf-1
        nc, ns = size(surfaces[k])
        i += nc * ns
    end
    
    nc, ns = size(surfaces[i_surf])
    n += (j-1) * nc + i
    return n
end

function fmm_to_vlm_index(system::System, n)
    n_counter = 0
    for k in 1:length(system.surfaces)
        nc, ns = size(system.surfaces[k])
        if n_counter + nc * ns >= n # found the surface
            i_plus = n - n_counter
            i = mod(i_plus-1, nc) + 1
            j = div(i_plus-1, nc) + 1
            i_surf = k
            return i_surf, i, j
        end
        n_counter += nc * ns
    end
end

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, system::System, i_body)
    
    # vlm index
    i_surf, i, j = fmm_to_vlm_index(system, i_body)

    # get panel
    panel = system.surfaces[i_surf][i, j]

    # update buffer
    buffer[1:3, i_buffer] .= panel.rcp
    buffer[4, i_buffer] = norm(panel.rbc - panel.rtr)
    buffer[5, i_buffer] = system.Γ[i_body]
    buffer[6:8,i_buffer] .= panel.rtr
    buffer[9:11,i_buffer] .= panel.rtl
    buffer[12:14,i_buffer] .= panel.rbl
    buffer[15:17,i_buffer] .= panel.rbr
    buffer[18,i_buffer] = panel.core_size
end

function FastMultipole.data_per_body(system::System)
    return 18
end

# function reset!(system::System{TF}) where TF
#     system.potential .= zero(TF)
# end

function FastMultipole.get_position(system::System, i)
    i_surf, i, j = fmm_to_vlm_index(system, i)
    panel = system.surfaces[i_surf][i, j]
    return panel.rcp
end

function FastMultipole.strength_dims(system::System)
    return 1
end

FastMultipole.get_n_bodies(system::System) = length(system.Γ)

FastMultipole.body_to_multipole!(system::System, args...) = FastMultipole.body_to_multipole_quad!(Panel{Dipole}, system, args...)

function FastMultipole.direct!(target_system, target_index, ::DerivativesSwitch{PS,VS,GS}, source_system::System, source_buffer, source_index) where {PS,VS,GS}
    @inbounds for i_source in source_index
        v1 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 1)
        v2 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 2)
        v3 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 3)
        v4 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 4)
        gamma = FastMultipole.get_strength(source_buffer, source_system, i_source)[1]
        cs = source_buffer[18, i_source]
        @inbounds for j_target in target_index
            target = FastMultipole.get_position(target_system, j_target)
            if VS
                v = bound_induced_velocity(target-v1, target-v2, true, cs)
                v += bound_induced_velocity(target-v2, target-v3, true, cs)
                v += bound_induced_velocity(target-v3, target-v4, true, cs)
                v += bound_induced_velocity(target-v4, target-v1, true, cs)
                FastMultipole.set_velocity!(target_system, j_target, v)
            end
        end
    end
end

function FastMultipole.buffer_to_target_system!(target_system::System, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    # get values
    # TF = eltype(target_buffer)
    # scalar_potential = PS ? FastMultipole.get_scalar_potential(target_buffer, i_buffer) : zero(TF)
    # velocity = VS ? FastMultipole.get_velocity(target_buffer, i_buffer) : zero(SVector{3,TF})
    # velocity_gradient = GS ? FastMultipole.get_velocity_gradient(target_buffer, i_buffer) : zero(SMatrix{3,3,TF,9})

    # update system
    # target_system.potential[i_POTENTIAL[1], i_target] = scalar_potential

    # target_system.potential[i_VELOCITY, i_target] .= velocity
    @warn "A VortexLattice.System object should not be used as a target in an FMM call."
end
