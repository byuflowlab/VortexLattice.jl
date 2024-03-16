#####
##### overload fmm compatibility functions for `system::System.fmm_panels`
#####

Base.getindex(system::System, i, ::FLOWFMM.Position) = system.fmm_panels[i].rcp
Base.getindex(system::System, i, ::FLOWFMM.Radius) = system.fmm_panels[i].radius
# Base.getindex(system::System, i, ::FLOWFMM.VectorPotential) = view(g.potential,2:4,i)
# Base.getindex(system::System, i, ::FLOWFMM.ScalarPotential) = g.potential[1,i]
Base.getindex(system::System, i, ::FLOWFMM.Velocity) = system.fmm_Vcp[i]
# Base.getindex(system::System, i, ::FLOWFMM.VelocityGradient) = reshape(view(g.potential,i_VELOCITY_GRADIENT,i),3,3)
Base.getindex(system::System, i, ::FLOWFMM.ScalarStrength) = system.fmm_panels[i].gamma
Base.getindex(system::System, i, ::FLOWFMM.Body) = system.fmm_panels[i], system.fmm_Vcp[i]
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
    panel, velocity = val
    system.fmm_panels[i] = panel
    system.fmm_Vcp[i] = velocity
    return nothing
end
# function Base.setindex!(system::System, val, i, ::FLOWFMM.ScalarPotential)
#     g.potential[i_POTENTIAL[1],i] = val
# end
# function Base.setindex!(system::System, val, i, ::FLOWFMM.VectorPotential)
#     g.potential[i_POTENTIAL[2:4],i] .= val
# end
function Base.setindex!(system::System, val, i, ::FLOWFMM.Velocity)
    system.fmm_Vcp[i] = val
end
# function Base.setindex!(system::System, val, i, ::FLOWFMM.VelocityGradient)
#     reshape(g.potential[i_VELOCITY_GRADIENT,i],3,3) .= val
# end
FLOWFMM.get_n_bodies(system::System) = length(system.fmm_panels)
Base.eltype(::System{TF}) where TF = TF

FLOWFMM.buffer_element(system::System) = deepcopy(system.fmm_panels[1]), deepcopy(system.fmm_Vcp[1])

FLOWFMM.B2M!(system::System, branch, bodies_index, harmonics, expansion_order) = 
    FLOWFMM.B2M!_quadpanel(system, branch, bodies_index, harmonics, system.fmm_p, FLOWFMM.UniformNormalDipolePanel())

# no velocity gradient needed if the target is a `::VortexLattice.System` or a FLOWFMM.ProbeSystem without velocity gradient
function FLOWFMM.direct!(target_system::Union{System,FLOWFMM.ProbeSystem{<:Any,<:Any,<:Any,Vector{SVector{3,<:Any}},Nothing}}, target_index, source_system::System{TF}, source_index) where TF
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
            finite_core = false
            this_V = zero(SVector{3,TF})
            this_V += bound_induced_velocity(r1, r2, finite_core, core_size)
            this_V += bound_induced_velocity(r2, r3, finite_core, core_size)
            this_V += bound_induced_velocity(r3, r4, finite_core, core_size)
            this_V += bound_induced_velocity(r4, r1, finite_core, core_size)
            V += this_V * gamma
        end

        # add to target
        target_system[j_target,FLOWFMM.VELOCITY] = target_system[j_target,FLOWFMM.VELOCITY] + V
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