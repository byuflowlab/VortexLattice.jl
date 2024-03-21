#####
##### overload fmm compatibility functions for `system::System.fmm_panels`
#####

Base.getindex(system::System, i, ::FLOWFMM.Position) = system.fmm_panels[i].rcp
Base.getindex(system::System, i, ::FLOWFMM.Radius) = system.fmm_panels[i].radius
# Base.getindex(system::System{TF}, i, ::FLOWFMM.VectorPotential) where TF = zero(SVector{3,TF})
# Base.getindex(system::System{TF}, i, ::FLOWFMM.ScalarPotential) where TF = zero(TF)
Base.getindex(system::System{TF}, i, ::FLOWFMM.Velocity) where TF = zero(SVector{3,TF})
# Base.getindex(system::System{TF}, i, ::FLOWFMM.VelocityGradient) where TF = zero(SMatrix{3,3,TF,9})
Base.getindex(system::System, i, ::FLOWFMM.ScalarStrength) = system.fmm_panels[i].gamma
Base.getindex(system::System, i, ::FLOWFMM.Body) = system.fmm_panels[i]#, system.fmm_Vcp[i]
function Base.getindex(system::System, i, ::FLOWFMM.Vertex, i_vertex)
    if i_vertex == 1
        return system.fmm_panels[i].rtl
    elseif i_vertex == 2
        return system.fmm_panels[i].rbl
    elseif i_vertex == 3
        return system.fmm_panels[i].rbr
    else # i_vertex == 4
        return system.fmm_panels[i].rtr
    end
end
Base.getindex(system::System, i, ::FLOWFMM.Normal) = system.fmm_panels[i].ncp
function Base.setindex!(system::System, val, i, ::FLOWFMM.Body)
    # panel, velocity = val
    system.fmm_panels[i] = val
    # system.fmm_Vcp[i] = velocity
    return nothing
end
# function Base.setindex!(system::System, val, i, ::FLOWFMM.ScalarPotential)
#     return nothing
# end
# function Base.setindex!(system::System, val, i, ::FLOWFMM.VectorPotential)
#     g.potential[i_POTENTIAL[2:4],i] .= val
# end
# function Base.setindex!(system::System, val, i, ::FLOWFMM.Velocity)
#     system.fmm_Vcp[i] = val
# end
# function Base.setindex!(system::System, val, i, ::FLOWFMM.VelocityGradient)
#     reshape(g.potential[i_VELOCITY_GRADIENT,i],3,3) .= val
# end
FLOWFMM.get_n_bodies(system::System) = length(system.fmm_panels)
# Base.eltype(::System{TF}) where TF = TF

FLOWFMM.buffer_element(system::System) = deepcopy(system.fmm_panels[1])#, deepcopy(system.fmm_Vcp[1])

FLOWFMM.B2M!(system::System, branch, bodies_index, harmonics, expansion_order) = 
    FLOWFMM.B2M!_quadpanel(system, branch, bodies_index, harmonics, Val(system.fmm_p), FLOWFMM.UniformNormalDipolePanel())

# no velocity gradient needed if the target is a `::VortexLattice.System` or a FLOWFMM.ProbeSystem without velocity gradient
function FLOWFMM.direct!(target_system::FLOWFMM.ProbeSystem{<:Any,<:Any,<:Any,<:Any,Nothing}, target_index, source_system::System{TF}, source_index) where TF
    # loop over targets
    for j_target in target_index
        rcp = target_system[j_target,FLOWFMM.POSITION]
        V = zero(SVector{3,TF})

        # loop over source panels
        for i_source in source_index
            panel = source_system.fmm_panels[i_source]

            # get panel vertices, strength, and core size
            r11 = panel.rtl
            r12 = panel.rtr
            r21 = panel.rbl
            r22 = panel.rbr
            gamma = panel.gamma
            core_size = panel.core_size
            
            # move origin to control point
            r1 = rcp - r11
            r2 = rcp - r12
            r3 = rcp - r22
            r4 = rcp - r21

            # add induced velocity, ignoring finite core size for now
            finite_core = true
            this_V = zero(SVector{3,TF})
            this_V += bound_induced_velocity_fmm(r1, r2, finite_core, core_size)
            this_V += bound_induced_velocity_fmm(r2, r3, finite_core, core_size)
            this_V += bound_induced_velocity_fmm(r3, r4, finite_core, core_size)
            this_V += bound_induced_velocity_fmm(r4, r1, finite_core, core_size)
            V += this_V * gamma
        end

        # add to target
        target_system[j_target,FLOWFMM.VELOCITY] += V
    end
end

function FLOWFMM.direct!(target_system, target_index, source_system::System{TF}, source_index) where TF
    # loop over targets
    for j_target in target_index
        rcp = target_system[j_target,FLOWFMM.POSITION]
        V = zero(SVector{3,TF})
        Vgrad = zero(SMatrix{3,3,TF,9})

        # loop over source panels
        for i_source in source_index
            panel = source_system.fmm_panels[i_source]

            # get panel vertices, strength, and core size
            r11 = panel.rtl
            r12 = panel.rtr
            r21 = panel.rbl
            r22 = panel.rbr
            gamma = panel.gamma
            core_size = panel.core_size
            
            # move origin to control point
            r1 = rcp - r11
            r2 = rcp - r12
            r3 = rcp - r22
            r4 = rcp - r21

            # add induced velocity, ignoring finite core size for now
            finite_core = false
            this_V = zero(SVector{3,TF})
            this_V += bound_induced_velocity(r1, r2, finite_core, core_size)
            this_V += bound_induced_velocity(r2, r3, finite_core, core_size)
            this_V += bound_induced_velocity(r3, r4, finite_core, core_size)
            this_V += bound_induced_velocity(r4, r1, finite_core, core_size)
            V += this_V * gamma

            # velocity gradient
            this_Vgrad = zero(SMatrix{3,3,TF,9})
            this_Vgrad += bound_velocity_gradient(r1, r2)
            this_Vgrad += bound_velocity_gradient(r2, r3)
            this_Vgrad += bound_velocity_gradient(r3, r4)
            this_Vgrad += bound_velocity_gradient(r4, r1)
            Vgrad += this_Vgrad
        end

        # add to target
        target_system[j_target,FLOWFMM.VELOCITY] = target_system[j_target,FLOWFMM.VELOCITY] + V
        target_system[j_target,FLOWFMM.VELOCITY_GRADIENT] = target_system[j_target,FLOWFMM.VELOCITY_GRADIENT] + Vgrad
    end
end

#####
##### probes
#####
function update_n_probes!(probes::FLOWFMM.ProbeSystem{TF,Nothing,Nothing,TV,Nothing}, n) where {TF,TV}
    resize!(probes.position, n)
    resize!(probes.velocity, n)
end

# function update_n_probes!(probes::FLOWFMM.ProbeSystem{TF,Nothing,Nothing,TV,TVG}, n) where {TF,TV,TVG}
#     resize!(probes.position, n)
#     resize!(probes.velocity, n)
#     resize!(probes.velocity_gradient, n)
# end

function reset!(probes::FLOWFMM.ProbeSystem{TF,Nothing,Nothing,TV,Nothing}) where {TF,TV}
    for i in eachindex(probes.velocity)
        probes.velocity[i] = zero(eltype(probes.velocity))
    end
end

# function reset!(probes::FLOWFMM.ProbeSystem{TF,Nothing,Nothing,TV,TVG}) where {TF,TV,TVG}
#     for i in eachindex(probes.velocity)
#         probes.velocity[i] = zero(eltype(probes.velocity))
#         probes.velocity_gradient[i] = zero(eltype(probes.velocity_gradient))
#     end
# end

"Places probes at the control point of all surface panels, and additional probes at the center of all surface filaments"
function update_probes!(fmm_velocity_probes::FLOWFMM.ProbeSystem{<:Any,<:Any,<:Any,<:Any,<:Any}, surfaces::Vector{Matrix{SurfacePanel{TF}}}, i_start) where TF # wake corners
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
function update_probes!(fmm_velocity_probes::FLOWFMM.ProbeSystem{<:Any,<:Any,<:Any,<:Any,<:Any}, wakes::Vector{Matrix{WakePanel{TF}}}, nwake, i_start) where TF # wake corners
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
# function update_probes!(fmm_velocity_probes::FLOWFMM.ProbeSystem{<:Any,<:Any,<:Any,<:Any,<:Any}, wake_shedding_locations::Vector{Vector{SVector{3,TF}}}, i_start) where TF # wake transition points
#     i_shed = 1
#     for wsl in wake_shedding_locations
#         for location in wsl
#             fmm_velocity_probes.position[i_shed + i_start] = location
#             i_shed += 1
#         end
#     end
# end
