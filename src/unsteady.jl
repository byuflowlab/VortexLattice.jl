function simulate!(system::System, wake::ParticleField, frames::AbstractVector{<:ReferenceFrame}, maneuver!::Function, Vinf::Function, t_range;
        name="vortex_lattice_simulation", path="./vortex_lattice_simulation",
        vtk_args=(), fmm_wake_args=(), fmm_vehicle_args=(),
        nonlinear_args=(),
    )
    # create save path if it does not exist
    if !isnothing(path) && !isdir(path)
        mkpath(path)
    end

    # update VPM freestream
    wake.Uinf = Vinf
    
    # begin simulation
    i_step = 0
    for t in t_range
        dt = t_range[i_step + 2] - t
        println("\tstep $(i_step)/$(length(t_range)-1) at time $(t)")
        
        #------- controls -------#

        # update frames based on maneuver
        # (RPMs, tilting systems, prescribed trajectory, etc.)
        maneuver!(frames, system, wake, t)

        # update kinematic velocity due to rigid body motion
        # (structural deflections should be remembered from the previous step)
        # NOTE: this skips the top level frame, which is captured in system.fs
        kinematic_velocity!(system, frames; skip_top_level=true)
        
        #------- aerodynamics -------#

        # update freestream velocity based on the top level frame
        # (includes top level frame's velocity and rotation)
        fs = Freestream(frames[1], ref, Vinf(t))

        # wake-on-all
        wake.SFS(wake, FLOWVPM.BeforeUJ())
        wake_on_all!(wake, system; fmm_wake_args...)

        #--- solve the system ---#

        # unpack constant system parameters
        ref = system.reference[]
        symmetric = system.symmetric
        surface_id = system.surface_id
        wake_finite_core = system.wake_finite_core
        trailing_vortices = system.trailing_vortices
        xhat = system.xhat[]

        # unpack system storage (including state variables)
        previous_surfaces = system.previous_surfaces
        current_surfaces = system.surfaces
        properties = system.properties
        dproperties = system.dproperties
        wakes = system.wakes
        wake_velocities = system.V
        wake_shedding_locations = system.wake_shedding_locations
        AIC = system.AIC
        w = system.w
        dw = system.dw
        Γ = system.Γ
        dΓ = system.dΓ
        dΓdt = system.dΓdt
        Vcp = system.Vcp
        Vh = system.Vh
        Vv = system.Vv
        Vte = system.Vte

        # default arguments
        additional_velocity = nothing
        trailing_vortices .= false
        symmetric .= false
        
        # update the wake shedding location for this time step
        update_wake_shedding_locations!(wakes, wake_shedding_locations,
            current_surfaces, ref, fs, dt, additional_velocity, Vte,
            nwake, eta)

        # calculate/re-calculate AIC matrix (if necessary)
        # if surface_motion || calculate_influence_matrix
            influence_coefficients!(AIC, current_surfaces;
                symmetric, wake_shedding_locations,
                surface_id, trailing_vortices, xhat)
        # end

        # update the AIC matrix to use the new wake shedding locations
        update_trailing_edge_coefficients!(AIC, current_surfaces;
            symmetric, wake_shedding_location, trailing_vortices)

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

        vehicle_on_all!(system, wake; fmm_vehicle_args...)
        
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
        system.near_field_analysis[] = near_field_analysis

        # nonlinear airfoil analysis
        if nonlinear_analysis
            nonlinear_analysis!(system, ref, fs; nonlinear_args...)
        end
        
        #------- save state -------#

        if !isnothing(path)
            # VortexLattice system
            write_vtk(joinpath(path, name * "_step_$i_step"), system; vtk_args...)

            # FLOWVLM particle field
            save(joinpath(path, name * "_wake_step_$i_step"), wake)
        end

        #------- other solvers -------#

        # e.g. structures, acoustics, dynamics, etc.
        
        
        if i_step < length(t_range) - 1
            
            # propagate rigid-body kinematics
            propagate_kinematics!(system, frames, dt)

            # propagate wake
            FLOWVPM._euler(wake, dt; relax=true)

            # shed new particles
            shed_wake!(wake, system, dt) # TODO

            # increment step
            i_step += 1
            
        end
    end
    
end

function wake_on_all!(wake::ParticleField, system::System; fmm_wake_args...)
    update_probes!(system) # update probe positions and reset
    fmm!((wake, system.probes), wake; fmm_wake_args...) # solve N-body problem
    probes_to_surfaces!(system) # update Vv, Vh, Vcp, etc. based on probes
end

function vehicle_on_all!(system::System, wake::ParticleField; fmm_vehicle_args...)
    FastMultipole.reset_probes!(system)
    fmm!((wake, system.probes), system; fmm_vehicle_args...)
end
