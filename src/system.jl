"""
    System{TF}

Pre-allocated container for vortex lattice method inputs and states.

# Fields:
 - `AIC`: Aerodynamic influence coefficient matrix from the surface panels
 - `w`: Normal velocity at the control points from external sources and wakes
 - `Gamma`: Circulation strength of the surface panels
 - `Gammadot`: Derivative of the circulation strength with respect to time
 - `surfaces`: Surfaces, represented by matrices of surface panels (see [`SurfacePanel`](@ref))
 - `repeated_points`: Dictionary of the form `Dict((isurf, i) => [(jsurf1, j1),
    (jsurf2, j2)...]` which defines repeated trailing edge points.  Trailing edge
    point `i` on surface `isurf` is repeated on surface `jsurf1` at point `j1`,
    `jsurf2` at point `j2`, and so forth.
 - `wake_shedding_locations`: Wake shedding locations for each surface
 - `wakes`: Vector of wakes corresponding to each surface, represented
    by matrices of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where
    `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels.
 - `properties`: Surface panel properties for each surface
 - `trefftz`: Trefftz panels associated with each surface
 - `Vcp`: Velocities at the control points
 - `Vh`: Velocities at the center of the horizontal bound vortices on each surface
 - `Vv`: Velocities at the center of the vertical bound vortices on each surface
 - `Vw`: Velocities at the wake panel vertices for each wake
 - `dw`: Derivatives of the R.H.S. with respect to the freestream variables
 - `dGamma`: Derivatives of the circulation strength with respect to the freestream variables
 - `dproperties`: Derivatives of the panel properties with respect to the freestream variables
 - `reference`: Current reference parameters
 - `freestream`: Current freestream parameters
 - `additional_velocity`: Additional velocity function
 - `symmetric`: Flag indicating whether surface (and freestream parameters) are symmetric
 - `surface_id`: Surface ID for each surface.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.
 - `iwake`: Current number of chordwise wake panels to use for each surface.
 - `nwake`: Maximum number of chordwise wake panels to use for each surface.
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from
    each surface.
 - `xhat`: Direction in which trailing vortices are shed.
"""
mutable struct System{TF}
    AIC::Matrix{TF}
    w::Vector{TF}
    Gamma::Vector{TF}
    Gammadot::Vector{TF}
    surfaces::Vector{Matrix{SurfacePanel{TF}}}
    repeated_points::Dict{NTuple{2,Int},Vector{NTuple{2,Int}}}
    wake_shedding_locations::Vector{Vector{SVector{3,TF}}}
    wakes::Vector{Matrix{WakePanel{TF}}}
    properties::Vector{Matrix{PanelProperties{TF}}}
    trefftz::Vector{Vector{TrefftzPanel{TF}}}
    Vcp::Vector{Matrix{SVector{3, TF}}}
    Vh::Vector{Matrix{SVector{3, TF}}}
    Vv::Vector{Matrix{SVector{3, TF}}}
    Vw::Vector{Matrix{SVector{3,TF}}}
    dw::NTuple{5, Vector{TF}}
    dGamma::NTuple{5, Vector{TF}}
    dproperties::NTuple{5, Vector{Matrix{PanelProperties{TF}}}}
    reference::Reference{TF}
    freestream::Freestream{TF}
    additional_velocity # Function
    symmetric::Vector{Bool}
    surface_id::Vector{Int}
    wake_finite_core::Vector{Bool}
    iwake::Vector{Int}
    nwake::Vector{Int}
    eta::TF
    trailing_vortices::Vector{Bool}
    xhat::SVector{3,TF}
end

@inline Base.eltype(::Type{System{TF}}) where TF = TF
@inline Base.eltype(::System{TF}) where TF = TF

"""
    System([TF], surfaces, reference, freestream, nshed = 0; kwargs...)

Return an object of type `System` which holds current vortex lattice method
inputs and states.

# Arguments:
 - `TF`: Floating point type, defaults to the floating point type used by `surface`
 - `surfaces`: Initial surfaces
        Either:
         - One or more grids of shape (3, nc+1, ns+1) which represents lifting surfaces,
         or
         - One or more matrices of surface panels (see [`SurfacePanel`](@ref)) of
           shape (nc, ns)
        where `nc` is the number of chordwise panels and `ns` is the number of
        spanwise panels.
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Initial freestream parameters (see [`Freestream`](@ref))
 - `nshed`: Number of wake shedding events.

# Keyword Arguments:
 - `wakes`: Initial wake panels
        Either:
         - One or more grids of shape (3, nw+1, ns+1) which define wake vertices,
         or
         - One or more matrices of wake panels (see [`WakePanel`](@ref)) of
           shape (nw, ns)
        where `nw` is the number of chordwise wake panels and `ns` is the number
        of spanwise panels.  Defaults to no wake panels for each surface.
 - `surface_velocities`: Surface velocities at the beginning of the simulation,
    defined as described by the argument `surfaces`, except velocities are used
    in place of positions.  Defaults to zero surface motion/deformation.
 - `circulation`: Vector containing the initial circulation of all surface
    panels in the system.  Defaults to `zeros(N)` where `N` is the total number
    of surface panels in `surfaces`.
 - `circulation_rates`: Vector containing the initial circulation rates of all
    surface panels in the system.  Defaults to `zeros(N)` where `N` is the total number
    of surface panels in `surfaces`.
 - `wake_circulation`: Vector containing the initial wake circulation for all
    wake panels in the system.  Defaults to `zeros(N)` where `N` is the total
    number of wake panels in `wakes`.
 - `additional_velocity`: Function of the form `V = f(r)` which defines the
    additional velocity `V` at location `r` in the global coordinate frame.
    Defaults to no additional velocity field.
 - `symmetric`: Flags indicating whether the geometry for each surface is
    symmetric about the X-Z plane.  Note that applying symmetry to surfaces is
    only valid when the freestream conditions are symmetric as well.
    Defaults to `false` for each surface.
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.  Defaults to `true` for each surface.
 - `iwake`: Initial number of chordwise wake panels to use for each surface.
    Defaults to `size.(wakes, 1)`.
 - `nwake`: Maximum number of chordwise wake panels to use for each surface.
    Defaults to `fill(nshed, length(surfaces))`
 - `eta`: Time step fraction used to define separation between trailing
    edge and wake shedding location.  Typical values range from 0.2 to 0.3.
    Defaults to 0.2.
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from
    each surface (or wake if present).  Defaults to `true` for each surface.
 - `xhat`: Direction in which trailing vortices are shed.  Defaults
    to `[1, 0, 0]`
"""
System(args...; kwargs...)

# grid inputs, no provided type
function System(surfaces::AbstractVector{<:AbstractArray{TF, 3}}, reference,
    freestream, nshed=0; kwargs...) where TF
    return System(TF, surfaces, reference, freestream, nshed; kwargs...)
end

# grid inputs, provided type
function System(TF::Type, surfaces::AbstractVector{<:AbstractArray{<:Any, 3}},
    reference, freestream, nshed=0; nwake = fill(nshed, length(surfaces)), kwargs...)
    nc = size.(surfaces, 2) .- 1
    ns = size.(surfaces, 3) .- 1
    return System(TF, nc, ns, nwake; surfaces, reference, freestream, kwargs...)
end

# surface inputs, no provided type
function System(surfaces::AbstractVector{<:AbstractMatrix{SurfacePanel{TF}}},
    reference, freestream, nshed=0; kwargs...) where TF
    return System(TF, surfaces, reference, freestream, nshed; kwargs...)
end

# surface inputs, provided type
function System(TF::Type, surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}},
    reference, freestream, nshed=0; nwake = fill(nshed, length(surfaces)), kwargs...)
    nc = size.(surfaces, 1)
    ns = size.(surfaces, 2)
    return System(TF, nc, ns, nwake; surfaces, reference, freestream, kwargs...)
end

"""
    System(TF, nc, ns, nw; kwargs...)

Construct an object of type `System` which holds current vortex lattice
method inputs and states.

# Arguments:
 - `TF`: Floating point type.
 - `nc`: Number of chordwise panels on each surface.
 - `ns`: Number of spanwise panels on each surface.
 - `nw`: Number of chordwise wake panels for each surface.

# Keyword Arguments
 - `surfaces`: Initial surfaces
   - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
   or
   - Vector of matrices of shape (nc, ns) containing surface panels (see
   [`SurfacePanel`](@ref))
   where `nc` is the number of chordwise panels and `ns` is the number of
   spanwise panels.
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters for each time step (see [`Freestream`](@ref))
 - `wakes`: Vector of initial wakes corresponding to each surface, represented
    by matrices of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where
    `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels. Defaults to no wake panels for each surface
 - `circulation`: Vector containing the initial circulation of all surface
    panels in the system.  Defaults to `zeros(N)` where `N` is the total number
    of surface panels in `surfaces`.
 - `symmetric`: Flags indicating whether the geometry for each surface is
    symmetric about the X-Z plane.  Note that applying symmetry to surfaces is
    only valid when the freestream conditions are symmetric as well.
    Defaults to `false` for each surface.
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.  Defaults to `true` for each surface.
 - `nwake`: Maximum number of wake panels in the chordwise direction for each
    surface.  Defaults to `nshed` for all surfaces.
 - `eta`: Time step fraction used to define separation between trailing
    edge and wake shedding location.  Typical values range from 0.2 to 0.3.
    Defaults to 0.2.
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from
    each surface.  Defaults to `true` for each surface.
 - `xhat`: Direction in which trailing vortices are shed.  Defaults
    to `[1, 0, 0]`
"""
function System(TF::Type, nc, ns, nw; fcore=(c, Δs) -> 1e-3, kwargs...)
    @assert length(nc) == length(ns) == length(nw)
    # number of surfaces
    nsurf = length(nc)
    # total number of surface panels
    N = sum(nc .* ns)
    # initialize system storage
    AIC = Matrix{TF}(undef, N, N)
    w = Vector{TF}(undef, N)
    Γ = Vector{TF}(undef, N)
    Γdot = Vector{TF}(undef, N)
    surfaces = [Matrix{SurfacePanel{TF}}(undef, nc[i], ns[i]) for i = 1:nsurf]
    repeated_points = Dict{NTuple{2,Int},Vector{NTuple{2,Int}}}()
    wake_shedding_locations = [Vector{SVector{3,TF}}(undef, ns[i]+1) for i = 1:nsurf]
    wakes = [Matrix{WakePanel{TF}}(undef, nw[i], ns[i]) for i = 1:nsurf]
    properties = [Matrix{PanelProperties{TF}}(undef, nc[i], ns[i]) for i = 1:nsurf]
    trefftz = [Vector{TrefftzPanel{TF}}(undef, ns[i]) for i = 1:nsurf]
    Vcp = [Matrix{SVector{3,TF}}(undef, nc[i], ns[i]) for i = 1:nsurf]
    Vh = [Matrix{SVector{3,TF}}(undef, nc[i]+1, ns[i]) for i = 1:nsurf]
    Vv = [Matrix{SVector{3,TF}}(undef, nc[i], ns[i]+1) for i = 1:nsurf]
    Vw = [Matrix{SVector{3,TF}}(undef, nw[i]+1, ns[i]+1) for i = 1:nsurf]
    dw = Tuple(Vector{TF}(undef, N) for i = 1:5)
    dΓ = Tuple(Vector{TF}(undef, N) for i = 1:5)
    dproperties = Tuple([Matrix{PanelProperties{TF}}(undef, nc[i], ns[i]) for i = 1:length(surfaces)] for j = 1:5)
    reference = Reference{TF}(1.0, 1.0, 1.0, [0.0, 0.0, 0.0], 1.0)
    freestream = Freestream{TF}(1.0, 0.0, 0.0, [0.0, 0.0, 0.0])
    additional_velocity = r -> @SVector zeros(3)
    # set default variables/parameters
    Γ .= 0.0
    Γdot .= 0.0
    symmetric = fill(false, nsurf)
    surface_id = collect(1:nsurf)
    wake_finite_core = fill(true, nsurf)
    iwake = fill(0, nsurf)
    nwake = nw
    eta = 0.2
    trailing_vortices = fill(true, nsurf)
    xhat = SVector(1.0, 0.0, 0.0)
    # construct system
    system = System(AIC, w, Γ, Γdot, surfaces, repeated_points,
        wake_shedding_locations, wakes, properties, trefftz, Vcp, Vh, Vv,
        Vw, dw, dΓ, dproperties, reference, freestream, additional_velocity,
        symmetric, surface_id, wake_finite_core, iwake, nwake, eta,
        trailing_vortices, xhat)
    # update inputs to correspond to provided inputs
    update_inputs!(system, fcore, false; kwargs...)
    # return result
    return system
end

"""
    update_inputs!(system::System, fcore=(c, Δs) -> 1e-3, preserve_core_size=true; kwargs...)

Updates the inputs in `system` to correspond to the provided keyword arguments.

# Arguments
 - `system`: Object that holds unsteady vortex lattice method inputs and state
 - `fcore`: Function for setting the finite core size when generating surface/wake
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width.
 - `preserve_core_size`: Flag indicating whether the finite core size should be
    preserved from the previous set of surface/wake panels (if applicable).
    Defaults to `true`.

# Keyword Arguments
 - `steady`: Flag indicating whether to update inputs for a steady simulation.
 - `surfaces`: Surface geometry at the beginning of the simulation.
    Represented as either a:
     - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
    or
     - Vector of matrices of shape (nc, ns) containing surface panels (see [`SurfacePanel`](@ref))
    where `nc` is the number of chordwise panels and `ns` is the number of
    spanwise panels.
 - `surface_velocities`: Surface velocities at the beginning of the simulation,
    defined as described by the argument `surfaces`, except velocities are used
    in place of positions.
 - `wakes`: Wake panels at the beginning of the simulation.
    Represented as either a:
     - Vector of grids of shape (3, nw+1, ns+1) which represent wake geometry
    or
     - Vector of matrices of shape (nw, ns) containing wake panels
       (see [`WakePanel`](@ref))
    where `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels.
 - `circulation`: Vector containing the initial circulation of all surface
    panels in the system.
 - `circulation_rates`: Vector containing the initial circulation rates of all
    surface panels in the system.
 - `wake_circulation`: Vector containing the initial wake circulation for all
    wake panels in the system.
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))
 - `additional_velocity`: Function of the form `V = f(r)` which defines the
    additional velocity `V` at location `r` in the global coordinate frame.
 - `symmetric`: Flags indicating whether the geometry for each surface is
    symmetric.  Note that applying symmetry to surfaces is only theoretically
    valid when the freestream conditions are symmetric as well.
 - `surface_id`: Surface ID for each surface.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.
 - `iwake`: Initial number of chordwise wake panels to use for each surface.
 - `nwake`: Maximum number of chordwise wake panels to use for each surface.
 - `eta`: Time step fraction used to define separation between trailing
    edge and wake shedding location.  Typical values range from 0.2 to 0.3.
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from
    each surface.
 - `xhat`: Direction in which trailing vortices are shed if wakes are fixed.
"""
@inline function update_inputs!(system::System, fcore = (c, Δs) -> 1e-3, preserve_core_size=true;
    surfaces = nothing, surface_velocities = nothing, wakes = nothing,
    circulation = nothing, circulation_rates = nothing, wake_circulation = nothing,
    reference = nothing, freestream = nothing, additional_velocity = nothing,
    symmetric = nothing, surface_id = nothing, wake_finite_core = nothing,
    iwake = nothing, nwake = nothing, eta = nothing, trailing_vortices = nothing,
    xhat = nothing)

    # number of surfaces
    nsurf = length(system.surfaces)

    # update surfaces (if specified)
    if !isnothing(surfaces)
        # update surfaces
        update_surface_panels!(system, surfaces; fcore, preserve_core_size)
        # re-initialize repeated trailing edge points
        system.repeated_points = repeated_trailing_edge_points(system.surfaces)
        # re-initialize wake shedding locations
        initialize_wake_shedding_locations!(system, surfaces)
        # re-initialize undefined wake panels
        initialize_wake_panels!(system)
    end

    # update surface velocities (if specified)
    if !isnothing(surface_velocities)
        update_surface_velocities!(system, surface_velocities)
    end

    # check if existing wake panel storage is sufficient, replace if necessary
    if !isnothing(nwake)
        for isurf = 1:nsurf
            if size(system.wakes[isurf], 1) < nwake[isurf]
                # get surface/wake dimensions
                nc, ns = size(system.surfaces[isurf])
                nw = nwake[isurf]
                # initialize extra panels at trailing edge
                extra_panels = Matrix{WakePanel{TF}}(undef, tmp, nw - size(system.wakes[isurf], 1), ns)
                # update wake panel storage
                system.wakes[isurf] = vcat(system.wakes[isurf], extra_panels)
            end
        end
        copyto!(system.nwake, nwake)
    end

    # update wakes (if specified)
    if !isnothing(wakes)
        for isurf = 1:nsurf
            for I in CartesianIndices(wakes[isurf])
                system.wakes[isurf][I] = initial_wakes[isurf][I]
            end
            # update current number of chordwise wake panels
            system.iwake[isurf] = size(initial_wakes[isurf],1)
        end
    end

    # update circulation (if specified)
    if !isnothing(circulation)
        system.Gamma .= circulation
    end

    # update circulation rate (if specified)
    if !isnothing(circulation_rates)
        system.Gammadot .= circulation_rates
    end

    # update wake circulation (if specified)
    if !isnothing(wake_circulation)
        set_wake_circulation!(system, wake_circulation)
    end

    # update reference parameters (if specified)
    if !isnothing(reference)
        system.reference = reference
    end

    # update freestream parameters (if specified)
    if !isnothing(freestream)
        system.freestream = freestream
    end

    # update additional velocity function (if specified)
    if !isnothing(additional_velocity)
        system.additional_velocity = additional_velocity
    end

    # update run control parameters (if specified)
    if !isnothing(symmetric)
        copyto!(system.symmetric, symmetric)
    end
    if !isnothing(surface_id)
        copyto!(system.surface_id, surface_id)
    end
    if !isnothing(wake_finite_core)
        copyto!(system.wake_finite_core, wake_finite_core)
    end
    if !isnothing(iwake)
        # set new number of wake panels
        copyto!(system.iwake, iwake)
        # re-initialize undefined wake panels
        initialize_wake_panels!(system)
    end
    if !isnothing(eta)
        system.eta = eta
    end
    if !isnothing(trailing_vortices)
        copyto!(system.trailing_vortices, trailing_vortices)
    end
    if !isnothing(xhat)
        system.xhat = SVector{3}(xhat)
    end

    return system
end

function initialize_steady_simulation!(system)

    # set circulation rate to zero
    system.Gammadot .= 0.0

    # set surface velocities to zero
    for isurf = 1:length(system.surfaces)
        system.Vcp[isurf] .= Ref(@SVector zeros(3))
        system.Vh[isurf] .= Ref(@SVector zeros(3))
        system.Vv[isurf] .= Ref(@SVector zeros(3))
    end
end
