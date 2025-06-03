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

function simulate!(system::System{TF}, frames::AbstractVector{<:ReferenceFrame}, maneuver!::Function, Vinf::Function, t_range;
        wake_args=(), kwargs...) where TF
   
    # number of shed locations
    n = 0
    for surface in system.surfaces
        nc, ns = size(surface)
        n += (2*ns + 1) # * p_per_step
    end

    # construct particle field
    wake = ParticleField((n + div(n,8)) * length(t_range), TF; Uinf=Vinf, wake_args...)
    trailing_edge_particles = ParticleField(n, TF; wake_args...)
    
    # begin simulation
    simulate!(system, wake, trailing_edge_particles, frames, maneuver!, Vinf, t_range; kwargs...)

    return wake
end

function force!(val, force_val::Nothing)
    return val
end

function force!(val, force_val)
    val .= force_val
end

function simulate!(system::System, wake::ParticleField, trailing_edge_particles::ParticleField, frames::AbstractVector{<:ReferenceFrame}, maneuver!::Function, Vinf::Function, t_range;
        name="vortex_lattice_simulation", path="./vortex_lattice_simulation",
        vtk_args=(), fmm_wake_args=(), fmm_vehicle_args=(),
        nonlinear_analysis=false, nonlinear_args=(),
        eta=0.3, derivatives=false,
        trailing_vortices=fill(false, length(system.surfaces)),
        shedding_surfaces=fill(true, length(system.surfaces)),
        monitors,
    )
    # create save path if it does not exist
    if !isnothing(path) && !isdir(path)
        mkpath(path)
    end

    # begin simulation
    i_step = 0
    for t in t_range[1:end-1]
        dt = t_range[i_step + 2] - t
        println("\tstep $(i_step+1)/$(length(t_range)) at time $(t)")
        
        #------- reset system -------#

        Vcp = system.Vcp
        Vh = system.Vh
        Vv = system.Vv
        Vte = system.Vte
        for isurf in 1:length(system.surfaces)
            Vcp[isurf] .= Ref(zero(eltype(Vcp[isurf])))
            Vh[isurf] .= Ref(zero(eltype(Vh[isurf])))
            Vv[isurf] .= Ref(zero(eltype(Vv[isurf])))
            Vte[isurf] .= Ref(zero(eltype(Vte[isurf])))
        end
        system.w .= zero(eltype(system.w))
        for i in eachindex(system.dw)
            system.dw[i] .= zero(eltype(system.dw[i]))
        end
        FastMultipole.reset!(system.probes)

        # particle field
        FLOWVPM._reset_particles(wake)
        FLOWVPM._reset_particles_sfs(wake)

        #------- controls -------#

        # update frames based on maneuver
        # (RPMs, tilting systems, prescribed trajectory, etc.)
        maneuver!(frames, system, wake, t)

        # update kinematic velocity due to rigid body motion
        # (structural deflections should be remembered from the previous step)
        # NOTE: this skips the top level frame, which is captured in system.fs
        current_surfaces = system.surfaces
        kinematic_velocity!(Vcp, Vh, Vv, Vte, current_surfaces, frames; skip_top_level=true)

        #------- aerodynamics -------#

        # update freestream velocity based on the top level frame
        # (includes top level frame's velocity and rotation)
        ref = system.reference[]
        fs = Freestream(frames[1], ref, Vinf(t))
        system.freestream[] = fs

        # unpack constant system parameters
        symmetric = system.symmetric
        surface_id = system.surface_id
        wake_finite_core = system.wake_finite_core
        trailing_vortices = system.trailing_vortices
        xhat = system.xhat[]

        # unpack system storage (including state variables)
        previous_surfaces = system.previous_surfaces
        properties = system.properties
        dproperties = system.dproperties
        wakes = system.wakes
        wake_velocities = system.V
        wake_shedding_locations = system.wake_shedding_locations
        nwake = system.nwake
        AIC = system.AIC
        w = system.w
        dw = system.dw
        Γ = system.Γ
        dΓ = system.dΓ
        dΓdt = system.dΓdt
        
        # default arguments
        additional_velocity = nothing
        system.trailing_vortices .= trailing_vortices
        symmetric .= false
        
        # update the wake shedding location for this time step
        update_wake_shedding_locations!(wakes, wake_shedding_locations,
            current_surfaces, ref, fs, dt, additional_velocity, Vte,
            nwake, eta)

        # shed new particles representing unsteady effects
        n_sheds = shed_wake!(wake, current_surfaces, wake_shedding_locations, Γ)

        # wake-on-all
        wake.SFS(wake, FLOWVPM.BeforeUJ())
        wake_on_all!(wake, system; fmm_wake_args...)

        #--- solve the system ---#

        # calculate/re-calculate AIC matrix (if necessary)
        # if surface_motion || calculate_influence_matrix
            influence_coefficients!(AIC, current_surfaces;
                symmetric, wake_shedding_locations, particles=shedding_surfaces,
                surface_id, trailing_vortices, xhat)
        # end

        # # update the AIC matrix to use the new wake shedding locations
        # update_trailing_edge_coefficients!(AIC, current_surfaces;
        #     symmetric, wake_shedding_locations, trailing_vortices)

        # calculate RHS
        if derivatives
            normal_velocity_derivatives!(w, dw, current_surfaces, wakes,
                ref, fs; additional_velocity, Vcp, symmetric, nwake,
                surface_id, wake_finite_core, trailing_vortices, xhat)
        else
            normal_velocity!(w, current_surfaces, wakes, ref, fs;
                additional_velocity, Vcp, symmetric, nwake, surface_id,
                wake_finite_core, trailing_vortices, xhat)
        end

        # save (negative) previous circulation in dΓdt
        dΓdt .= .-Γ

        # solve for the new circulation
        if derivatives
            circulation_derivatives!(Γ, dΓ, AIC, w, dw)
        else
            circulation!(Γ, AIC, w)
        end

        # solve for dΓdt using finite difference `dΓdt = (Γ - Γp)/dt`
        dΓdt .+= Γ # add newly computed circulation
        dΓdt ./= dt # divide by corresponding time step

        #--- vehicle-on-all ---#

        # add new particles to trailing edge pfield
        update_locations!(trailing_edge_particles, current_surfaces, wake_shedding_locations)
        update_strengths!(trailing_edge_particles, current_surfaces, wake_shedding_locations, Γ; i_start=1, reset_strengths=true)
        vehicle_on_all!(system, trailing_edge_particles, wake; fmm_vehicle_args...)

        # update new particle strengths
        update_strengths!(wake, current_surfaces, wake_shedding_locations, Γ; i_start=wake.np - n_sheds + 1, reset_strengths=false)
        
        #--- forces and moments ---#

        # compute transient forces on each panel (if necessary)
        if derivatives
            near_field_forces_derivatives!(properties, dproperties,
                current_surfaces, wakes, ref, fs, Γ, dΓ; dΓdt,
                additional_velocity, Vh, Vv, symmetric, nwake,
                surface_id, wake_finite_core, wake_shedding_locations,
                trailing_vortices, xhat)
        else
            near_field_forces!(properties, current_surfaces, wakes,
                ref, fs, Γ; dΓdt, additional_velocity, Vh, Vv,
                symmetric, nwake, surface_id, wake_finite_core,
                wake_shedding_locations, trailing_vortices, xhat)
        end

        # save flag indicating that a near-field analysis has been performed
        system.near_field_analysis[] = true

        # nonlinear airfoil analysis
        if nonlinear_analysis
            nonlinear_analysis!(system, ref, fs; nonlinear_args...)
        end
        
        #------- other solvers -------#
        
        # e.g. structures, acoustics, dynamics, etc.
        
        #------- update state -------#
        
        
        #------- save state -------#

        if !isnothing(path)
            # VortexLattice system
            write_vtk(joinpath(path, name * "_step_$i_step"), system; vtk_args...)

            # FLOWVLM particle field
            FLOWVPM.save(wake, name * "_wake"; add_num=true, num=i_step, path, overwrite_time=i_step)
        end

        for monitor in monitors
            monitor(system, wake, i_step)
        end

        #------- propagate system -------#
        
        # propagate rigid-body kinematics
        propagate_kinematics!(system, frames, dt)
        
        # propagate wake
        FLOWVPM._euler(wake, dt; relax=true)
        
        # increment step
        i_step += 1
    end
            
    println("\tstep $(i_step+1)/$(length(t_range)) at time $(t_range[end])")
    
    #------- save state -------#

    if !isnothing(path)
        # VortexLattice system
        write_vtk(joinpath(path, name * "_step_$i_step"), system; vtk_args...)

        # FLOWVLM particle field
        FLOWVPM.save(wake, joinpath(path, name * "_wake_step_$i_step"); add_num=false)
    end
    
end

function wake_on_all!(wake::ParticleField, system::System; fmm_wake_args...)
    update_probes!(system) # update probe positions and reset
    fmm!((wake, system.probes), wake; velocity_gradient=SVector{2}(true,false), fmm_wake_args...) # solve N-body problem
    probes_to_surfaces!(system) # update Vv, Vh, Vcp, etc. based on probes
end

function shed_wake!(pfield::FLOWVPM.ParticleField, surfaces, wake_shedding_locations, Γ)
    # loop over surfaces
    iΓ = 0
    np = 0
    for isurf = 1:length(surfaces)
        surface = surfaces[isurf]
        wsl = wake_shedding_locations[isurf]
        nc, ns = size(surface)
        for j in 1:ns
            # strength
            iΓ += nc

            # get vertices
            panel = surface[end, j]
            r11 = bottom_left(panel)
            r12 = bottom_right(panel)
            r21 = wsl[j]
            r22 = wsl[j+1]

            # left particle
            r1, r2 = r21, r11
            Xp, Γp, sigmap = filament_to_particle(r1, r2, 0.0)
            FLOWVPM.add_particle(pfield, Xp, Γp, sigmap)
            np += 1
            
            # bottom particle
            r1, r2 = r22, r21
            Xp, Γp, sigmap = filament_to_particle(r1, r2, -Γ[iΓ])
            FLOWVPM.add_particle(pfield, Xp, Γp, sigmap)
            np += 1
            
        end
        
        # get vertices
        panel = surface[end, end]
        r11 = bottom_left(panel)
        r12 = bottom_right(panel)
        r21 = wsl[end-1]
        r22 = wsl[end]
        
        # right
        r1, r2 = r12, r22
        Xp, Γp, sigmap = filament_to_particle(r1, r2, 0.0)
        FLOWVPM.add_particle(pfield, Xp, Γp, sigmap)
        np += 1
    end

    return np
end

function update_locations!(trailing_edge_particles::FLOWVPM.ParticleField, surfaces, wake_shedding_locations)
    # empty particle field
    trailing_edge_particles.np = 0
    
    # loop over surfaces
    iΓ = 0
    for isurf = 1:length(surfaces)
        surface = surfaces[isurf]
        wsl = wake_shedding_locations[isurf]
        nc, ns = size(surface)
        for j in 1:ns
            # strength index
            iΓ += nc

            # get vertices
            panel = surface[end, j]
            r11 = bottom_left(panel)
            r12 = bottom_right(panel)
            r21 = wsl[j]
            r22 = wsl[j+1]

            # left particle
            r1, r2 = r21, r11
            Xp, Γp, sigmap = filament_to_particle(r1, r2, 1.0)
            FLOWVPM.add_particle(trailing_edge_particles, Xp, Γp, sigmap)

            # bottom particle
            r1, r2 = r22, r21
            Xp, Γp, sigmap = filament_to_particle(r1, r2, 1.0)
            FLOWVPM.add_particle(trailing_edge_particles, Xp, Γp, sigmap)

        end
        
        # get vertices
        panel = surface[end, end]
        r11 = bottom_left(panel)
        r12 = bottom_right(panel)
        r21 = wsl[end-1]
        r22 = wsl[end]

        # right
        r1, r2 = r12, r22
        Xp, Γp, sigmap = filament_to_particle(r1, r2, 1.0)
        FLOWVPM.add_particle(trailing_edge_particles, Xp, Γp, sigmap)
    end
end

function update_strengths!(pfield::FLOWVPM.ParticleField, surfaces, wake_shedding_locations, Γ; i_start=1, reset_strengths=true, unsteady=true)
    
    # total number of particles
    np = 0
    for isurf in 1:length(surfaces)
        nc, ns = size(surfaces[isurf])
        np += 2 * ns + 1
    end

    if reset_strengths
        for ip in i_start:i_start+np-1
            Gamma = FLOWVPM.get_Gamma(pfield, ip)
            Gamma .= zero(eltype(Gamma))
        end
    end

    # loop over surfaces
    iΓ = 0
    ip = i_start
    for isurf = 1:length(surfaces)
        surface = surfaces[isurf]
        wsl = wake_shedding_locations[isurf]
        nc, ns = size(surface)

        Γlast = zero(eltype(Γ))
        for j in 1:ns
            # strength index
            iΓ += nc
            # get vertices
            panel = surface[end, j]
            r11 = bottom_left(panel)
            r12 = bottom_right(panel)
            r21 = wsl[j]
            r22 = wsl[j+1]

            # left particle
            r1, r2 = r21, r11
            thisΓ = Γ[iΓ]
            Xp, Γp, sigmap = filament_to_particle(r1, r2, thisΓ-Γlast)
            Gamma = FLOWVPM.get_Gamma(pfield, ip)
            Gamma .+= Γp
            ip += 1
            
            # bottom particle
            r1, r2 = r22, r21
            Xp, Γp, sigmap = filament_to_particle(r1, r2, thisΓ)
            Gamma = FLOWVPM.get_Gamma(pfield, ip)
            if unsteady
                if VortexLattice.DEBUG[]
                    @show j, Gamma Γp Gamma .+ Γp
                end
                Gamma .+= Γp
            else
                Gamma .= zero(eltype(Gamma))
            end
            ip += 1

            # recurse Γ
            Γlast = thisΓ

        end
        
        # get vertices
        panel = surface[end, end]
        r11 = bottom_left(panel)
        r12 = bottom_right(panel)
        r21 = wsl[end-1]
        r22 = wsl[end]

        # right
        r1, r2 = r12, r22
        Xp, Γp, sigmap = filament_to_particle(r1, r2, Γ[iΓ])
        Gamma = FLOWVPM.get_Gamma(pfield, ip)
        Gamma .+= Γp
        ip += 1
    end
end

function vehicle_on_all!(system::System, trailing_edge_pfield::ParticleField, wake::ParticleField; fmm_vehicle_args...)
    FastMultipole.reset!(system.probes)
    fmm!((wake, system.probes), (system, trailing_edge_pfield); velocity_gradient=SVector{2}(true,false), fmm_vehicle_args...)
end

function unsteady_tmp!(wake::ParticleField, current_surfaces, Γ, wake_shedding_locations)
    
    iΓ = 0
    np = 0
    for isurf in 1:length(current_surfaces)
        surface = current_surfaces[isurf]
        wsl = wake_shedding_locations[isurf]
        nc, ns = size(surface)
        for i_tmp in 1:ns
            iΓ += nc
            v1 = wsl[i_tmp+1]
            v2 = wsl[i_tmp]
            Xp, Γp, sigmap = filament_to_particle(v1, v2, Γ[iΓ])
            FLOWVLM.add_particle(wake, Xp, Γp, sigmap)
            np += 1
        end
    end

    return np
end
