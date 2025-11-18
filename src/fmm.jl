#------- FastMultipole compatibility functions for System -------#

# function vlm_to_fmm_index(system::System, i_surf, i, j)
#     surfaces = system.surfaces
#     n = 0
#     for k in 1:i_surf-1
#         nc, ns = size(surfaces[k])
#         i += nc * ns
#     end
    
#     nc, ns = size(surfaces[i_surf])
#     n += (j-1) * nc + i
#     return n
# end

function fmm_to_vlm_index(system::System, n)
    n_counter = 0
    for i_surf in 1:length(system.surfaces)
        nc, ns = size(system.surfaces[i_surf])
        if n_counter + nc * ns + ns >= n # found the surface
            i_plus = n - n_counter
            if i_plus <= nc * ns
                i = mod(i_plus-1, nc) + 1
                j = div(i_plus-1, nc) + 1
                return i_surf, i, j
            else
                return i_surf, -1, i_plus - nc * ns
            end
        end
        n_counter += nc * ns + ns
    end
end

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, system::System, i_body)
    
    # vlm index
    i_surf, i, j = fmm_to_vlm_index(system, i_body)

    # gamma index
    i_Γ = 0
    for i_surf2 in 1:i_surf-1
        nc2, ns2 = size(system.surfaces[i_surf2])
        i_Γ += ns2 * nc2
    end
    # i_Γ is now at end of previous surface

    # dimensions
    nc, ns = size(system.surfaces[i_surf])

    if i > 0 # regular panel
        # get panel
        panel = system.surfaces[i_surf][i, j]

        # vertices
        rtl = panel.rtl
        rtr = panel.rtr
        rbr = panel.rbr
        rbl = panel.rbl
        rcp = panel.rcp

        # adjust gamma index
        i_Γ += (j-1) * nc + i

        # whether or not to include bottom bound vortex
        # include_bottom_vortex = i != nc 
        # ^^^ this really isn't necessary since it is always used

    else # wake interface panel

        # get panel info
        panel = system.surfaces[i_surf][nc, j]
        rtl = panel.rbl
        rtr = panel.rbr
        wsl = system.wake_shedding_locations[i_surf]
        rbl = wsl[j]
        rbr = wsl[j+1]
        rcp = 0.25 * (rtl + rtr + rbr + rbl)

        # adjust gamma index
        i_Γ += nc * j

    end
        
    # core size
    core_size = panel.core_size

    # update buffer
    buffer[1:3, i_buffer] .= rcp
    buffer[4, i_buffer] = 0.5 * max(norm(rbl - rtr), norm(rbr - rtl)) + core_size
    buffer[5, i_buffer] = system.Γ[i_Γ]
    buffer[6:8,i_buffer] .= rtl
    buffer[9:11,i_buffer] .= rtr
    buffer[12:14,i_buffer] .= rbr
    buffer[15:17,i_buffer] .= rbl
    buffer[18,i_buffer] = core_size
    # buffer[19,i_buffer] = one(eltype(buffer)) # always include bottom vortex

end

function FastMultipole.data_per_body(system::System)
    return 18
end

# function reset!(system::System{TF}) where TF
#     system.potential .= zero(TF)
# end

function FastMultipole.get_position(system::System, i)
    i_surf, i, j = fmm_to_vlm_index(system, i)
    if i > 0
        panel = system.surfaces[i_surf][i, j]
        return panel.rcp
    else
        # wake interface panel
        panel = system.surfaces[i_surf][end, j]
        wsl = system.wake_shedding_locations[i_surf]
        return 0.25 * (panel.rbl + panel.rbr + wsl[j] + wsl[j+1])
    end
end

function FastMultipole.strength_dims(system::System)
    return 1
end

FastMultipole.has_vector_potential(system::System) = true

function FastMultipole.get_n_bodies(system::System) 
    n_bodies = length(system.Γ) # regular panels
    # add wake interface panels
    for isurf in eachindex(system.surfaces)
        n_bodies += size(system.surfaces[isurf], 2)
    end
    return n_bodies
end

# FastMultipole.body_to_multipole!(system::System, args...) = FastMultipole.body_to_multipole_quad!(Panel{Dipole}, system, args...)

function body_to_multipole_vl!(multipole_coefficients, harmonics, x1, x2, center, gamma, expansion_order)
    # delta
    xu = x2 - x1
    x0 = x1 - center

    # vector strength
    Gamma = gamma * xu / norm(xu)

    # update values
    FastMultipole.body_to_multipole_filament!(FastMultipole.Filament{FastMultipole.Vortex}, multipole_coefficients, harmonics, x0, xu, Gamma, expansion_order)
end

function FastMultipole.body_to_multipole!(system::System, multipole_coefficients, buffer::Matrix, center, bodies_index, harmonics, expansion_order)
    # loop over bodies
    for i_body in bodies_index
       
        # extract vertices from buffer
        rtl = FastMultipole.get_vertex(buffer, system, i_body, 1)
        rtr = FastMultipole.get_vertex(buffer, system, i_body, 2)
        rbr = FastMultipole.get_vertex(buffer, system, i_body, 3)
        rbl = FastMultipole.get_vertex(buffer, system, i_body, 4)

        # extract strength from buffer
        gamma = FastMultipole.get_strength(buffer, system, i_body)[1]

        # top bound vortex
        body_to_multipole_vl!(multipole_coefficients, harmonics, rtl, rtr, center, gamma, expansion_order)

        # right vortex
        body_to_multipole_vl!(multipole_coefficients, harmonics, rtr, rbr, center, gamma, expansion_order)

        # left vortex
        body_to_multipole_vl!(multipole_coefficients, harmonics, rbl, rtl, center, gamma, expansion_order)

        # bottom bound vortex
        # if buffer[19, i_body] > 0.0 # omit trailing edge vortices # always include bottom vortex
        body_to_multipole_vl!(multipole_coefficients, harmonics, rbr, rbl, center, gamma, expansion_order)
        # end

    end
end

function FastMultipole.direct!(target_system, target_index, ::DerivativesSwitch{PS,VS,GS}, source_system::System, source_buffer, source_index) where {PS,VS,GS}
    @inbounds for j_target in target_index
        target = FastMultipole.get_position(target_system, j_target)
        v = zero(SVector{3,eltype(target_system)})
        @inbounds for i_source in source_index
            v1 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 1)
            v2 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 2)
            v3 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 3)
            v4 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 4)
            gamma = FastMultipole.get_strength(source_buffer, source_system, i_source)[1]
            cs = source_buffer[18, i_source]
            # include_bottom = source_buffer[19, i_source] > 0.0

            if VS
                this_v = bound_induced_velocity(target-v1, target-v2, true, cs)
                this_v += bound_induced_velocity(target-v2, target-v3, true, cs)
                # if include_bottom
                    this_v += bound_induced_velocity(target-v3, target-v4, true, cs)
                # end
                this_v += bound_induced_velocity(target-v4, target-v1, true, cs)
                v += this_v * gamma
            end
        end
        if VS
            FastMultipole.set_gradient!(target_system, j_target, v)
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

#------- FastMultipole compatibility functions for Vector{<:WakePanel} -------#

struct FilamentWrapper{TF}
    wakes::Vector{Matrix{WakePanel{TF}}}
end

Base.eltype(filaments::FilamentWrapper{TF}) where TF = TF

function fmm_to_vlm_index(filaments::FilamentWrapper, n)
    n_counter = 0
    wakes = filaments.wakes
    for k in eachindex(wakes)
        nc, ns = size(wakes[k])
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

function FastMultipole.source_system_to_buffer!(buffer, i_buffer, filaments::FilamentWrapper, i_body)
    
    # vlm index
    wakes = filaments.wakes
    i_surf, i, j = fmm_to_vlm_index(filaments, i_body)

    # get panel
    panel = wakes[i_surf][i, j]

    # update buffer
    buffer[1:3, i_buffer] .= 0.5 * (panel.rtl + panel.rbr)
    buffer[4, i_buffer] = 0.5 * norm(panel.rbl - panel.rtr) + panel.core_size
    buffer[5, i_buffer] = panel.gamma
    buffer[6:8,i_buffer] .= panel.rtl
    buffer[9:11,i_buffer] .= panel.rtr
    buffer[12,i_buffer] = panel.core_size
end

function FastMultipole.data_per_body(wakes::FilamentWrapper)
    return 12
end

function FastMultipole.get_position(filaments::FilamentWrapper, i)
    wakes = filaments.wakes
    i_surf, i, j = fmm_to_vlm_index(filaments, i)
    panel = wakes[i_surf][i, j]
    return 0.5 * (panel.rtl + panel.rbr)
end

function FastMultipole.strength_dims(filaments::FilamentWrapper)
    return 1
end

FastMultipole.has_vector_potential(filaments::FilamentWrapper) = true

function FastMultipole.get_n_bodies(filaments::FilamentWrapper)
    nwakes = 0
    for wake in filaments.wakes
        nc, ns = size(wake)
        nwakes += nc * ns
    end
    return nwakes
end

function FastMultipole.body_to_multipole!(filaments::FilamentWrapper, multipole_coefficients, buffer::Matrix, center, bodies_index, harmonics, expansion_order)
    # loop over bodies
    for i_body in bodies_index
       
        # extract vertices from buffer
        rtl = FastMultipole.get_vertex(buffer, filaments, i_body, 1)
        rtr = FastMultipole.get_vertex(buffer, filaments, i_body, 2)

        # extract strength from buffer
        gamma = FastMultipole.get_strength(buffer, filaments, i_body)[1]

        # top bound vortex
        body_to_multipole_vl!(multipole_coefficients, harmonics, rtl, rtr, center, gamma, expansion_order)
    end
end

function FastMultipole.direct!(target_system, target_index, ::DerivativesSwitch{PS,VS,GS}, source_system::FilamentWrapper{TF}, source_buffer, source_index) where {PS,VS,GS,TF}
    @inbounds for j_target in target_index
        target = FastMultipole.get_position(target_system, j_target)
        v = SVector{3,TF}(0.0, 0.0, 0.0)
        @inbounds for i_source in source_index
            v1 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 1)
            v2 = FastMultipole.get_vertex(source_buffer, source_system, i_source, 2)
            gamma = FastMultipole.get_strength(source_buffer, source_system, i_source)[1]
            cs = source_buffer[12, i_source]
            if VS
                v += bound_induced_velocity(target-v1, target-v2, true, cs) * gamma
            end
        end
        FastMultipole.set_gradient!(target_system, j_target, v)
    end
end

function FastMultipole.buffer_to_target_system!(target_system::FilamentWrapper, i_target, ::FastMultipole.DerivativesSwitch{PS,VS,GS}, target_buffer, i_buffer) where {PS,VS,GS}
    @warn "A VortexLattice.System.wakes object should not be used as a target in an FMM call."
end
