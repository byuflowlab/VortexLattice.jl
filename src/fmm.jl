#####
##### overload fmm compatibility functions for `system::System.filaments`
#####

Base.getindex(system::System, i, ::FastMultipole.Position) = system.filaments[i].xc
Base.getindex(system::System, i, ::FastMultipole.Radius) = norm(system.filaments[i].x2 - system.filaments[i].xc)
# Base.getindex(system::System{TF}, i, ::FastMultipole.VectorPotential) where TF = zero(SVector{3,TF})
# Base.getindex(system::System{TF}, i, ::FastMultipole.ScalarPotential) where TF = zero(TF)
Base.getindex(system::System{TF}, i, ::FastMultipole.Velocity) where TF = zero(SVector{3,TF})
# Base.getindex(system::System{TF}, i, ::FastMultipole.VelocityGradient) where TF = zero(SMatrix{3,3,TF,9})
Base.getindex(system::System, i, ::FastMultipole.ScalarStrength) = system.filaments[i].gamma
Base.getindex(system::System, i, ::FastMultipole.Body) = system.filaments[i]#, system.fmm_Vcp[i]

function Base.getindex(system::System, i, ::FastMultipole.Vertex, i_vertex)
    if i_vertex == 1
        return system.filaments[i].x1
    elseif i_vertex == 2
        return system.filaments[i].x2
    else # i_vertex == 4
        throw("only 2 vertices exist for a vortex filament; requested $i_vertex")
    end
end

function Base.setindex!(system::System, val, i, ::FastMultipole.Body)
    system.filaments[i] = val
end

# function Base.setindex!(system::System, val, i, ::FastMultipole.ScalarPotential)
#     return nothing
# end
# function Base.setindex!(system::System, val, i, ::FastMultipole.VectorPotential)
#     g.potential[i_POTENTIAL[2:4],i] .= val
# end
# function Base.setindex!(system::System, val, i, ::FastMultipole.Velocity)
#     system.fmm_Vcp[i] = val
# end
# function Base.setindex!(system::System, val, i, ::FastMultipole.VelocityGradient)
#     reshape(g.potential[i_VELOCITY_GRADIENT,i],3,3) .= val
# end

FastMultipole.get_n_bodies(system::System) = length(system.filaments)

# Base.eltype(::System{TF}) where TF = TF

FastMultipole.buffer_element(system::System) = deepcopy(system.filaments[1])#, deepcopy(system.fmm_Vcp[1])

FastMultipole.B2M!(system::System, branch, bodies_index, harmonics, expansion_order) =
    nothing
#    FastMultipole.B2M!_filament(system, branch, bodies_index, harmonics, expansion_order, FastMultipole.Vortex)

function FastMultipole.direct!(target_system, target_index, derivatives_switch::FastMultipole.DerivativesSwitch{<:Any,<:Any,VS,GS}, source_system::System{TF}, source_index) where {TF,VS,GS}
    # loop over targets
    for j_target in target_index
        rcp = target_system[j_target,FastMultipole.POSITION]
        V = zero(SVector{3,TF})
        Vgrad = zero(SMatrix{3,3,TF,9})

        # loop over source filaments
        for i_source in source_index
            filament = source_system.filaments[i_source]

            # get panel vertices, strength, and core size
            x1 = filament.x1
            x2 = filament.x2
            gamma = filament.gamma
            core_size = filament.core_size

            # move origin to control point
            r1 = rcp - x1
            r2 = rcp - x2

            # add induced velocity
            finite_core = true
			if VS
            	V += bound_induced_velocity_fmm(r1, r2, finite_core, core_size) * gamma
			end

            # velocity gradient
			if GS
            	Vgrad += bound_velocity_gradient(r1, r2, finite_core, core_size) * gamma
			end
        end

        # add to target
		VS && (target_system[j_target,FastMultipole.VELOCITY] = target_system[j_target,FastMultipole.VELOCITY] + V)
		GS && (target_system[j_target,FastMultipole.VELOCITY_GRADIENT] = target_system[j_target,FastMultipole.VELOCITY_GRADIENT] + Vgrad)
    end
end

#####
##### probes
#####
function update_n_probes!(probes::FastMultipole.ProbeSystem{TF,Nothing,Nothing,TV,Nothing}, n) where {TF,TV}
    resize!(probes.position, n)
    resize!(probes.velocity, n)
end

function update_n_probes!(probes::FastMultipole.ProbeSystem{TF,Nothing,Nothing,TV,TVG}, n) where {TF,TV,TVG}
    resize!(probes.position, n)
    resize!(probes.velocity, n)
    resize!(probes.velocity_gradient, n)
end

function reset!(probes::FastMultipole.ProbeSystem{TF,Nothing,Nothing,TV,Nothing}) where {TF,TV}
    for i in eachindex(probes.velocity)
        probes.velocity[i] = zero(eltype(probes.velocity))
    end
end

function reset!(probes::FastMultipole.ProbeSystem{TF,Nothing,Nothing,TV,TVG}) where {TF,TV,TVG}
    for i in eachindex(probes.velocity)
        probes.velocity[i] = zero(eltype(probes.velocity))
        probes.velocity_gradient[i] = zero(eltype(probes.velocity_gradient))
    end
end

"Places probes at the control point of all surface panels, and additional probes at the center of all surface filaments"
function update_probes!(fmm_velocity_probes::FastMultipole.ProbeSystem, surfaces::Vector{Matrix{SurfacePanel{TF}}}, i_start) where TF # wake corners
    i_probe = 1

    # place probes at control points
    for surface in surfaces
        for panel in surface
            fmm_velocity_probes.position[i_probe + i_start] = panel.rcp
            i_probe += 1
        end
    end

    # place probes at vortex filament centers (for force computation)
    for surface in surfaces
        nc, ns = size(surface)
        for i in 1:ns
            for j in 1:nc

                # top
                fmm_velocity_probes.position[i_probe + i_start] = surface[j,i].rtc
                i_probe += 1

                # left
                fmm_velocity_probes.position[i_probe + i_start] = left_center(surface[j,i])
                i_probe += 1
            end

            # NOTE: As the trailing edge panels have the same strength as the wake transition panel, their bottom vortex
            #       has zero strength and can be neglected when computing the force. Similar behavior to a Kutta condition.
            #       For now, the bottom vortex is ignored.
            # # bottom
            # fmm_velocity_probes.position[i_probe + i_start] = panel.rbc
            # i_probe += 1
        end

        for j in 1:nc
            # right
            fmm_velocity_probes.position[i_probe + i_start] = right_center(surface[j,ns])
            i_probe += 1
        end
    end
end

"Places probes at the corners of all wake panels"
function update_probes!(fmm_velocity_probes::FastMultipole.ProbeSystem{<:Any,<:Any,<:Any,<:Any,<:Any}, wakes::Vector{Matrix{WakePanel{TF}}}, nwake, i_start) where TF # wake corners
    i_corner = i_start + 1
    for (nc, wake) in zip(nwake, wakes)
        if nc > 0
            ns = size(wake,2) # number of spanwise panels

            # get wake_shedding_locations (first row of wake points)
            for i in 1:ns
                fmm_velocity_probes.position[i_corner] = wake[1,i].rtl
                i_corner += 1
            end
            fmm_velocity_probes.position[i_corner] = wake[1,ns].rtr
            i_corner += 1

            # bottom left of each panel
            for i in 1:ns
                for j in 1:nc
                    fmm_velocity_probes.position[i_corner] = wake[j,i].rbl
                    i_corner += 1
                end
            end

            # bottom right of the last column
            for j in 1:nc
                fmm_velocity_probes.position[i_corner] = wake[j,ns].rbr
                i_corner += 1
            end
        end
    end
end

# "Places probes at all wake shedding locations"
# function update_probes!(fmm_velocity_probes::FastMultipole.ProbeSystem{<:Any,<:Any,<:Any,<:Any,<:Any}, wake_shedding_locations::Vector{Vector{SVector{3,TF}}}, i_start) where TF # wake transition points
#     i_shed = 1
#     for wsl in wake_shedding_locations
#         for location in wsl
#             fmm_velocity_probes.position[i_shed + i_start] = location
#             i_shed += 1
#         end
#     end
# end
