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

function ForcesMonitor(nt::Int, TF=Float64; frame=Wind())
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
        for i in eachindex(wake)
            panel = wake[i]
            if any(isnan.(panel.core_size))
                error("NaN detected in wake panel core size on surface $isurf, panel $i")
            end
            if any(isnan.(panel.gamma))
                error("NaN detected in wake panel circulation on surface $isurf, panel $i")
            end
            if any(isnan.(panel.rtl))
                error("NaN detected in wake panel rtl on surface $isurf, panel $i")
            end
            if any(isnan.(panel.rbl))
                error("NaN detected in wake panel rbl on surface $isurf, panel $i")
            end
            if any(isnan.(panel.rtr))
                error("NaN detected in wake panel rtr on surface $isurf, panel $i")
            end
            if any(isnan.(panel.rbr))
                error("NaN detected in wake panel rbr on surface $isurf, panel $i")
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

function simulate!(system::System, wake::ParticleField, frames::AbstractVector{<:ReferenceFrame}, maneuver!::Function, Vinf::Function, t_range, Ωinf=(t)->SVector{3}(0.0, 0.0, 0.0);
        name="vortex_lattice_simulation", path="./vortex_lattice_simulation",
        vtk_args=(trailing_vortices=false, write_wakes=false), fmm_wake_args=(), fmm_vehicle_args=(),
        derivatives=false, nonlinear_analysis=false, nonlinear_args=(),
        eta=0.3, 
        particle_trailing_methods=fill(OverlapPPS(1.3, 2), length(system.surfaces)),
        particle_unsteady_methods=fill(OverlapPPS(1.3, 2), length(system.surfaces)),
        trailing_vortices=fill(false, length(system.surfaces)),
        shedding_surfaces=fill(true, length(system.surfaces)),
        monitors=(),
        calculate_influence_matrix=true
    )
    # create save path if it does not exist
    if !isnothing(path) && !isdir(path)
        mkpath(path)
    end

    # empty wake shedding locations
    # empty_wake_shedding_locations = fill(nothing, length(system.surfaces))

    # one row of wake panels per surface
    system.nwake .= 1
    for isurf in 1:length(system.surfaces)
        if shedding_surfaces[isurf]
            panels = Matrix{WakePanel{eltype(system.Γ)}}(undef, 1, size(system.surfaces[isurf], 2))
            for i in 1:size(system.surfaces[isurf], 2)
                panels[1, i] = WakePanel{eltype(system.Γ)}(SVector{3,eltype(system.Γ)}(0.0, 0.0, 0.0),
                                                        SVector{3,eltype(system.Γ)}(0.0, 0.0, 0.0),
                                                        SVector{3,eltype(system.Γ)}(0.0, 0.0, 0.0),
                                                        SVector{3,eltype(system.Γ)}(0.0, 0.0, 0.0),
                                                        zero(eltype(system.Γ)),
                                                        zero(eltype(system.Γ)))
            end
            system.wakes[isurf] = panels
        end
    end
    
    # wrap in wake object
    trailing_edge_filaments = FilamentWrapper(system.wakes)

    # velocity at wake vertices
    for isurf in 1:length(system.surfaces)
        nc, ns = size(system.wakes[isurf])
        system.V[isurf] = zeros(SVector{3,eltype(system.Γ)}, nc+1, ns+1)
    end

    # unpack system properties
    symmetric = system.symmetric
    surface_id = system.surface_id
    system.trailing_vortices .= trailing_vortices

    # freestream for initial step
    ref = system.reference[]
    vinf = Vinf(t_range[1])
    Ω = Ωinf(t_range[1])
    fs = velocity_to_freestream(vinf, Ω)
    system.freestream[] = fs

    # begin simulation
    i_step = 0
    for t in t_range
        println("\tstep $(i_step)/$(length(t_range)-1) at time $(t)")
        
        #------- reset system -------#

        Vcp = system.Vcp
        Vh = system.Vh
        Vv = system.Vv
        Vte = system.Vte
        V = system.V
        for isurf in 1:length(system.surfaces)
            Vcp[isurf] .= Ref(zero(eltype(Vcp[isurf])))
            Vh[isurf] .= Ref(zero(eltype(Vh[isurf])))
            Vv[isurf] .= Ref(zero(eltype(Vv[isurf])))
            Vte[isurf] .= Ref(zero(eltype(Vte[isurf])))
            V[isurf] .= Ref(zero(eltype(V[isurf])))
        end
        system.w .= zero(eltype(system.w))
        for i in eachindex(system.dw)
            system.dw[i] .= zero(eltype(system.dw[i]))
        end

        # particle field
        FLOWVPM._reset_particles(wake)
        FLOWVPM._reset_particles_sfs(wake)

        #------- controls -------#

        # update frames based on maneuver
        # (RPMs, tilting systems, prescribed trajectory, etc.)
        dynamics_toggle = maneuver!(frames, system, wake, t)

        # update kinematic velocity due to rigid body motion
        # (structural deflections should be remembered from the previous step)
        # NOTE: this skips the top level frame, which is captured in system.fs
        # CORRECTION: just changed this to not skip the top level frame
        current_surfaces = system.surfaces
        kinematic_velocity!(Vcp, Vh, Vv, Vte, current_surfaces, frames; skip_top_level=false)
        
        #------- aerodynamics -------#

        # unpack constant system parameters
        wake_finite_core = system.wake_finite_core
        xhat = system.xhat[]

        # unpack system storage (including state variables)
        # previous_surfaces = system.previous_surfaces
        properties = system.properties
        dproperties = system.dproperties
        wakes = system.wakes
        # wake_velocities = system.V
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
        symmetric .= false

        # if the first timestep, set up the first wake panels for particle shedding later
        dt = i_step == length(t_range) - 1 ? t_range[end] - t_range[end-1] : t_range[i_step + 2] - t_range[i_step + 1]
        if i_step == 0
            # align "wake shedding locations" with the trailing edge
            dt = t_range[2] - t_range[1]
            additional_velocity = nothing
            update_wake_shedding_locations!(wakes, wake_shedding_locations,
                current_surfaces, ref, fs, dt, additional_velocity, Vte, nwake, eta)
            
            # initial wake panels are an extension of the Kutta panel for the first timestep
            initial_wake_panels!(wakes, wake_shedding_locations, current_surfaces, eta)
        end
        
        # update trailing edge filaments with the previous circulation solution
        update_trailing_edge_filaments!(trailing_edge_filaments, current_surfaces, Γ)

        # wake-on-all
        wake.SFS(wake, FLOWVPM.BeforeUJ())
        # println("\n\n~~~ BEFORE WAKE ON ALL ~~~\n")
        # check_for_nans(system)
        # check_for_nans(wake)
        wake_on_all!(system, wake, trailing_edge_filaments; fmm_wake_args...)
        # println("\n\n~~~ AFTER WAKE ON ALL ~~~\n")
        # check_for_nans(system)
        # check_for_nans(wake)

        #--- solve the system ---#

        # calculate/re-calculate AIC matrix (if necessary)
        if calculate_influence_matrix
            influence_coefficients!(AIC, current_surfaces;
                symmetric, wake_shedding_locations,
                # ignore_trailing_edges = shedding_surfaces,
                surface_id, trailing_vortices, xhat,
                force_finite_core = fill(true, length(current_surfaces)))
        end

        # # update the AIC matrix to use the new wake shedding locations
        # update_trailing_edge_coefficients!(AIC, current_surfaces;
        #     symmetric, wake_shedding_locations, trailing_vortices)

        # calculate RHS
        if derivatives
            normal_velocity_derivatives!(w, dw, current_surfaces, wakes,
                ref, fs; additional_velocity, Vcp, symmetric, nwake,
                surface_id, wake_finite_core, trailing_vortices, xhat, include_wakes=false)
        else
            normal_velocity!(w, current_surfaces, wakes, ref, fs;
                additional_velocity, Vcp, symmetric, nwake, surface_id,
                wake_finite_core, trailing_vortices, xhat, include_wakes=false)
        end
        # @show w
        # throw("here2")

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

        #--- nonlinear airfoil analysis ---#
        if nonlinear_analysis
            call_near_field_forces!(system)
            nonlinear_analysis!(system, ref, fs; nonlinear_args...)
        end

        #--- vehicle-on-all ---#

        # solve n-body problem
        # @show V[1][1,1] V[1][1,end]
        # DEBUG[] = true
        vehicle_on_all!(system, wake, trailing_edge_filaments; fmm_vehicle_args...)
        # DEBUG[] = false
        # @show V[1][1,1] V[1][1,end]
        # throw(ErrorException("STOP HERE"))

        #--- forces and moments ---#

        if nonlinear_analysis
            update_section_forces!(system)
        end

        # compute transient forces on each panel (if necessary)
        if derivatives
            near_field_forces_derivatives!(properties, dproperties,
                current_surfaces, wakes, ref, fs, Γ, dΓ; dΓdt,
                additional_velocity, Vh, Vv, symmetric, nwake,
                surface_id, wake_finite_core, wake_shedding_locations,
                trailing_vortices, xhat,
                calculate_vlm_induced=false,
                skip_nonlinear_surfaces=nonlinear_analysis) # we've already calculated the induced velocity
                                             # in vehicle_on_all!
        else
            near_field_forces!(properties, current_surfaces, wakes,
                ref, fs, Γ; dΓdt, additional_velocity, Vh, Vv,
                symmetric, nwake, surface_id, wake_finite_core,
                wake_shedding_locations, trailing_vortices, xhat,
                calculate_vlm_induced=false,
                skip_nonlinear_surfaces=nonlinear_analysis) # we've already calculated the induced velocity
                                             # in vehicle_on_all!
        end
        
        #------- other solvers -------#
        
        # e.g. structures, acoustics, dynamics, etc.
        
        #------- update state -------#
        
        
        #------- save state -------#

        if !isnothing(path)
            # VortexLattice system
            write_vtk(joinpath(path, name * "_step_$i_step"), system; vtk_args...) # trailing_edge_list=.!shedding_surfaces, vtk_args...)

            # FLOWVLM particle field

            # check wake for NaNs
            FLOWVPM.save(wake, name * "_wake"; add_num=true, num=i_step, path, overwrite_time=i_step)
            
            # save trailing edge filaments
            # write_vtk(joinpath(path, name * "_filaments_step_$i_step"), trailing_edge_filaments)
        end

        for monitor in monitors
            monitor(system, wake, i_step)
        end

        #------- propagate system -------#

        if i_step < length(t_range)

            #--- state evolution ---#
            
            # propagate wake
            FLOWVPM._euler(wake, dt; relax=true)

            # dynamics function
            # if dynamics_toggle
            #     apply_dynamics!(system, frames)
            # end

            # calculate next step's wake trailing edge
            this_V = nothing # ignore wake- and vehicle-induced velocity for wake shedding location update
            update_vpm_shedding_TE!(wakes, ref, fs, dt, additional_velocity, this_V) # uses current step's freestream

            # store trailing edge location for next step's wsl
            store_trailing_edge!(wake_shedding_locations, current_surfaces)
            
            # propagate rigid-body kinematics
            propagate_kinematics!(system, frames, dt)

            # next step's freestream
            idx = i_step == length(t_range) - 1 ? i_step + 1 : i_step + 2
            vinf = Vinf(t_range[idx])
            Ω = Ωinf(t_range[idx])
            # fs = Freestream(frames[1], ref, vinf)
            fs = velocity_to_freestream(vinf, Ω)
            system.freestream[] = fs

            # update wake shedding locations / wake leading edge
            this_Vte = nothing # ignore wake- and vehicle-induced velocity for wake shedding location update
            update_wake_shedding_locations_unsteady!(wakes, wake_shedding_locations,
                current_surfaces, ref, fs, dt, additional_velocity, this_Vte, nwake, eta) # uses next step's freestream

            #--- shed new wake particles ---#

            shed_wake!(wake, system,  dt, 
                particle_trailing_methods, particle_unsteady_methods)

            # update wake shedding locations based on wake and vehicle
            # accounts for vehicle-induced, wake-induced, freestream,
            # and kinematic velocities
            # update_vpm_shedding_LE!(wakes, ref, fs, dt, additional_velocity, V)

            #--- update locations for the next step ---#

            if !isnothing(path)
                # VortexLattice system
                write_vtk(joinpath(path, name * "_postshed_step_$i_step"), system; write_wakes=true, vtk_args...) # trailing_edge_list=.!shedding_surfaces, vtk_args...)

                # FLOWVLM particle field

                # check wake for NaNs
                FLOWVPM.save(wake, name * "_postshed_wake"; add_num=true, num=i_step, path, overwrite_time=i_step)
                
                # save trailing edge filaments
                # write_vtk(joinpath(path, name * "_filaments_step_$i_step"), trailing_edge_filaments)
            end
        end

        # increment step
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

abstract type WakeSheddingMethod end

struct NoShed <: WakeSheddingMethod end

struct SigmaPPS{TF} <: WakeSheddingMethod
    sigma::TF
    p_per_step::Int
end

struct SigmaOverlap{TF} <: WakeSheddingMethod
    sigma::TF
    overlap::TF
end

struct OverlapPPS{TF} <: WakeSheddingMethod
    overlap::TF
    p_per_step::Int
end

function shed_wake!(pfield::FLOWVPM.ParticleField, system, dt, 
        shedding_trailing::AbstractVector{<:WakeSheddingMethod}, shedding_unsteady::AbstractVector{<:WakeSheddingMethod})
    # shed trailing edge particles
    shed_trailing_edge!(pfield, system.surfaces, system.wakes, system.Γ, shedding_trailing)

    # shed unsteady particles
    shed_unsteady!(pfield, system.surfaces, system.wakes, system.dΓdt, dt, shedding_unsteady)
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
            shed_particles!(pfield, r1, r2, Γthis - Γlast, method)

            # recurse
            Γlast = Γthis
        end

        # get vertices
        panel = wake[1, end]
        r1 = top_right(panel)
        r2 = bottom_right(panel)

        # shed right particles
        shed_particles!(pfield, r1, r2, Γlast, method)
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
            shed_particles!(pfield, r1, r2, Γ, method)
        end
    end
end

function shed_particles!(pfield, r1, r2, Γ, method::OverlapPPS)
    # shed particles with overlap and p_per_step
    overlap = method.overlap
    p_per_step = method.p_per_step
    sigma = norm(r2 - r1) * overlap / p_per_step
    return shed_particles!(pfield, r1, r2, Γ, SigmaPPS(sigma, p_per_step))
end

function shed_particles!(pfield, r1, r2, Γ, method::SigmaOverlap)
    # shed particles with sigma and overlap
    sigma = method.sigma
    overlap = method.overlap
    p_per_step = ceil(Int, overlap * norm(r2 - r1) / sigma)
    return shed_particles!(pfield, r1, r2, Γ, SigmaPPS(sigma, p_per_step))
end

function shed_particles!(pfield, r1, r2, Γ, method::SigmaPPS)
    # shed particles with sigma and p_per_step
    sigma = method.sigma
    p_per_step = method.p_per_step

    # add particles
    distance_vector = (r2 - r1) / p_per_step
    Xp = r1 + distance_vector * 0.5
    Γp = Γ * distance_vector
    for i in 1:p_per_step
        FLOWVPM.add_particle(pfield, Xp, Γp, sigma; circulation=Γ)
        Xp += distance_vector
    end
end

function shed_particles!(pfield, r1, r2, Γ, method::NoShed)
    # do not shed particles
    return nothing
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
