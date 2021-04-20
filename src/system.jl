"""
    AbstractSystem{TF}

Supertype for [`SteadySystem`](@ref) and [`UnsteadySystem`](@ref).
"""
abstract type AbstractSystem{TF} end

"""
    SteadySystem{TF} <: AbstractSystem{TF}

Pre-allocated container for steady vortex lattice method inputs and states.

# Fields:
 - `AIC`: Aerodynamic influence coefficient matrix from the surface panels
 - `w`: Normal velocity at the control points from external sources
 - `Gamma`: Circulation strength of the surface panels
 - `surfaces`: Surfaces, represented by matrices of surface panels
 - `properties`: Surface panel properties for each surface
 - `trefftz`: Trefftz panels associated with each surface
 - `dw`: Derivatives of the R.H.S. with respect to the freestream variables
 - `dGamma`: Derivatives of the circulation strength with respect to the freestream variables
 - `dproperties`: Derivatives of the panel properties with respect to the freestream variables
 - `reference`: Current reference parameters
 - `freestream`: Current freestream parameters
 - `symmetric`: Flag indicating whether surface (and freestream parameters) are symmetric
 - `surface_id`: Surface ID for each surface.
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from each surface
 - `xhat`: Direction for trailing vortices
"""
mutable struct SteadySystem{TF} <: AbstractSystem{TF}
    AIC::Matrix{TF}
    w::Vector{TF}
    Gamma::Vector{TF}
    surfaces::Vector{Matrix{SurfacePanel{TF}}}
    properties::Vector{Matrix{PanelProperties{TF}}}
    trefftz::Vector{Vector{TrefftzPanel{TF}}}
    dw::NTuple{5, Vector{TF}}
    dGamma::NTuple{5, Vector{TF}}
    dproperties::NTuple{5, Vector{Matrix{PanelProperties{TF}}}}
    reference::Reference{TF}
    freestream::Freestream{TF}
    symmetric::Vector{Bool}
    surface_id::Vector{Int}
    trailing_vortices::Vector{Bool}
    xhat::Vector{Float64}
end

@inline Base.eltype(::Type{SteadySystem{TF}}) where TF = TF
@inline Base.eltype(::SteadySystem{TF}) where TF = TF

# --- Constructors --- #

"""
    SteadySystem([TF], surfaces, reference, freestream; kwargs...)

Construct an object of type `SteadySystem` which holds current vortex lattice
method inputs and states.

# Arguments:
 - `TF`: Floating point type, defaults to the floating point type used by `surface`
 - `surfaces`:
   - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
   or
   - Vector of matrices of shape (nc, ns) containing surface panels (see
   [`SurfacePanel`](@ref))
   where `nc` is the number of chordwise panels and `ns` is the number of
   spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))

# Keyword Arguments
 - `symmetric`: Flags indicating whether the geometry for each surface is
    symmetric.  Note that applying symmetry to surfaces is only theoretically
    valid when the freestream conditions are symmetric as well.
    Defaults to `false` for each surface.
 - `surface_id`: Surface ID for each surface. By default, each surface has its
    own ID.
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from
    each surface.  Defaults to `true` for each surface.
 - `xhat`: Direction in which trailing vortices are shed.  Defaults
    to `[1, 0, 0]`
 - `near_field_analysis`: Flag indicating whether to perform a near field analysis
    to obtain surface panel properties.  Defaults to `true`.
 - `derivatives`: Flag indicating whether to calculate derivatives with respect
    to the freestream variables.  Defaults to `true`
 - `fcore`: Function for setting the finite core size when generating surface
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width.
"""
SteadySystem()

# grid inputs, no provided type
function SteadySystem(surfaces::AbstractVector{<:AbstractArray{TF, 3}}, reference,
    freestream; kwargs...) where TF
    return SteadySystem(TF, surfaces, reference, freestream; kwargs...)
end

# grid inputs, provided type
function SteadySystem(TF::Type, surfaces::AbstractVector{<:AbstractArray{<:Any, 3}},
    reference, freestream; kwargs...)
    nc = size.(surfaces, 2) .- 1
    ns = size.(surfaces, 3) .- 1
    return SteadySystem(TF, nc, ns; surfaces, reference, freestream, kwargs...)
end

# surface inputs, no provided type
function SteadySystem(surfaces::AbstractVector{<:AbstractMatrix{SurfacePanel{TF}}},
    reference, freestream; kwargs...) where TF
    return SteadySystem(TF, surfaces, reference, freestream; kwargs...)
end

# surface inputs, provided type
function SteadySystem(TF::Type, surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}},
    reference, freestream; kwargs...)
    nc = size.(surfaces, 1)
    ns = size.(surfaces, 2)
    return SteadySystem(TF, nc, ns; surfaces, reference, freestream, kwargs...)
end

"""
    SteadySystem(TF, nc, ns; kwargs...)

Construct an object of type `SteadySystem` which holds current vortex lattice
method inputs and states.

# Arguments:
 - `TF`: Floating point type.
 - `nc`: Number of chordwise panels on each surface.
 - `ns`: Number of spanwise panels on each surface.

# Additional Keyword Arguments
 - `surfaces`:
   - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
   or
   - Vector of matrices of shape (nc, ns) containing surface panels (see
   [`SurfacePanel`](@ref))
   where `nc` is the number of chordwise panels and `ns` is the number of
   spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))
 - `symmetric`: Flags indicating whether the geometry for each surface is
    symmetric.  Note that applying symmetry to surfaces is only theoretically
    valid when the freestream conditions are symmetric as well.
    Defaults to `false` for each surface.
 - `surface_id`: Surface ID for each surface. By default, each surface has its
    own ID.
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from
    each surface.  Defaults to `true` for each surface.
 - `xhat`: Direction in which trailing vortices are shed.  Defaults
    to `[1, 0, 0]`
 - `fcore`: Function for setting the finite core size when generating surface
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width.
"""
function SteadySystem(TF::Type, nc, ns; fcore=(c, Δs) -> 1e-3, kwargs...)
    # check inputs for consistency
    @assert length(nc) == length(ns)
    # number of surfaces
    nsurf = length(nc)
    # total number of panels
    N = sum(nc .* ns)
    # initialize system storage
    AIC = Matrix{TF}(undef, N, N)
    w = Vector{TF}(undef, N)
    Γ = Vector{TF}(undef, N)
    surfaces = [Matrix{SurfacePanel{TF}}(undef, nc[i], ns[i]) for i = 1:nsurf]
    properties = [Matrix{PanelProperties{TF}}(undef, nc[i], ns[i]) for i = 1:nsurf]
    trefftz = [Vector{TrefftzPanel{TF}}(undef, ns[i]) for i = 1:nsurf]
    dw = Tuple(Vector{TF}(undef, N) for i = 1:5)
    dΓ = Tuple(Vector{TF}(undef, N) for i = 1:5)
    dproperties = Tuple([Matrix{PanelProperties{TF}}(undef, nc[i], ns[i]) for i = 1:length(surfaces)] for j = 1:5)
    # set reference and freestream parameters
    reference = Reference{TF}(1.0, 1.0, 1.0, [0.0, 0.0, 0.0], 1.0)
    freestream = Freestream{TF}(1.0, 0.0, 0.0, [0.0, 0.0, 0.0])
    # set default control parameters
    symmetric = fill(false, nsurf)
    surface_id = collect(1:nsurf)
    trailing_vortices = fill(true, nsurf)
    xhat = [1.0, 0.0, 0.0]
    # construct system
    system = SteadySystem(AIC, w, Γ, surfaces, properties, trefftz, dw, dΓ,
        dproperties, reference, freestream, symmetric, surface_id,
        trailing_vortices, xhat)
    # update inputs to correspond to provided inputs
    update_inputs!(system, fcore, false; kwargs..., )
    return system
end

"""
    update_inputs!(system::SteadySystem, fcore=(c, Δs) -> 1e-3, preserve_core_size=true; kwargs...)

Updates the inputs in `system` to correspond to the provided keyword arguments.

# Arguments
 - `system`: Object that holds steady vortex lattice method inputs and state
 - `fcore`: Function for setting the finite core size when generating surface
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width.
 - `preserve_core_size`: Flag indicating whether the finite core size should be
    preserved from the previous set of panels (if applicable). Defaults to `true`.

# Keyword Arguments
 - `surfaces`:
   - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
   or
   - Vector of matrices of shape (nc, ns) containing surface panels (see
   [`SurfacePanel`](@ref))
   where `nc` is the number of chordwise panels and `ns` is the number of
   spanwise panels
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))
 - `symmetric`: Flags indicating whether the geometry for each surface is
    symmetric.  Note that applying symmetry to surfaces is only theoretically
    valid when the freestream conditions are symmetric as well.
 - `surface_id`: Surface ID for each surface.
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from
    each surface.
 - `xhat`: Direction in which trailing vortices are shed.
"""
function update_inputs!(system::SteadySystem, fcore = (c, Δs) -> 1e-3, preserve_core_size=true;
    surfaces = nothing, reference = nothing, freestream = nothing,
    symmetric = nothing, surface_id = nothing, trailing_vortices = nothing,
    xhat = nothing)

    # update surfaces (if specified)
    if !isnothing(surfaces)
        update_surface_panels!(system, surfaces; fcore, preserve_core_size)
    end

    # update reference parameters (if specified)
    if !isnothing(reference)
        system.reference = reference
    end

    # update freestream parameters (if specified)
    if !isnothing(freestream)
        system.freestream = freestream
    end

    # update run control parameters (if specified)
    if !isnothing(symmetric)
        copyto!(system.symmetric, symmetric)
    end
    if !isnothing(surface_id)
        copyto!(system.surface_id, surface_id)
    end
    if !isnothing(trailing_vortices)
        copyto!(system.trailing_vortices, trailing_vortices)
    end
    if !isnothing(xhat)
        copyto!(system.xhat, xhat)
    end

    return system
end

"""
    UnsteadySystem{TF}

Pre-allocated container for unsteady vortex lattice method inputs and states.

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
 - `Vte`: Velocities at the trailing edge vertices for each surface
 - `Vw`: Velocities at the wake panel vertices for each wake
 - `dw`: Derivatives of the R.H.S. with respect to the freestream variables
 - `dGamma`: Derivatives of the circulation strength with respect to the freestream variables
 - `dproperties`: Derivatives of the panel properties with respect to the freestream variables
 - `reference`: Current reference parameters
 - `freestream`: Current freestream parameters
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
mutable struct UnsteadySystem{TF} <: AbstractSystem{TF}
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
    Vte::Vector{Vector{SVector{3, TF}}}
    Vw::Vector{Matrix{SVector{3,TF}}}
    dw::NTuple{5, Vector{TF}}
    dGamma::NTuple{5, Vector{TF}}
    dproperties::NTuple{5, Vector{Matrix{PanelProperties{TF}}}}
    reference::Reference{TF}
    freestream::Freestream{TF}
    symmetric::Vector{Bool}
    surface_id::Vector{Int}
    wake_finite_core::Vector{Bool}
    iwake::Vector{Int}
    nwake::Vector{Int}
    trailing_vortices::Vector{Bool}
    xhat::Vector{TF}
end

@inline Base.eltype(::Type{UnsteadySystem{TF}}) where TF = TF
@inline Base.eltype(::UnsteadySystem{TF}) where TF = TF

"""
    UnsteadySystem([TF], surfaces, reference, freestream, nshed; kwargs...)

Return an object of type `UnsteadySystem` which holds current unsteady vortex
lattice method inputs and states.

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
 - `nshed`: Number of time steps at which wake panels are shed. Used solely to
    determine the default value of keyword argument `nwake`.

# Keyword Arguments:
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
 - `nwake`: Maximum number of chordwise wake panels to initialize for each
    surface. Defaults to `ntime` for each surface.
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from
    each surface.  Defaults to `true` for each surface.
 - `xhat`: Direction in which trailing vortices are shed.  Defaults
    to `[1, 0, 0]`
 - `fcore`: Function for setting the finite core size when generating surface
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width.
"""
UnsteadySystem(args...; kwargs...)

# grid inputs, no provided type
function UnsteadySystem(surfaces::AbstractVector{<:AbstractArray{TF, 3}}, ref, fs,
    nshed; kwargs...) where TF
    return UnsteadySystem(TF, surfaces, reference, freestream, nshed; kwargs...)
end

# grid inputs, provided type
function UnsteadySystem(TF::Type, surfaces::AbstractVector{<:AbstractArray{<:Any, 3}},
    reference, freestream, nshed; nwake = fill(nshed, length(surfaces)), kwargs...)
    nc = size.(surfaces, 2) .- 1
    ns = size.(surfaces, 3) .- 1
    return UnsteadySystem(TF, nc, ns, nwake; surfaces, reference, freestream, kwargs...)
end

# surface inputs, no provided type
function UnsteadySystem(surfaces::AbstractVector{<:AbstractMatrix{SurfacePanel{TF}}},
    reference, freestream, nshed; kwargs...) where TF
    return UnsteadySystem(TF, surfaces, reference, freestream, nshed; kwargs...)
end

# surface inputs, provided type
function UnsteadySystem(TF::Type, surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}},
    reference, freestream, nshed; nwake = fill(nshed, length(surfaces)), kwargs...)
    nc = size.(surfaces, 1)
    ns = size.(surfaces, 2)
    return UnsteadySystem(TF, nc, ns, nwake; surfaces, reference, freestream, kwargs...)
end

"""
    UnsteadySystem(TF, nc, ns, nw; kwargs...)

Construct an object of type `UnsteadySystem` which holds current vortex lattice
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
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from
    each surface.  Defaults to `true` for each surface.
 - `xhat`: Direction in which trailing vortices are shed.  Defaults
    to `[1, 0, 0]`
 - `fcore`: Function for setting the finite core size when generating surface
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width.
"""
function UnsteadySystem(TF::Type, nc, ns, nw; fcore=(c, Δs) -> 1e-3, kwargs...)
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
    Vte = [Vector{SVector{3,TF}}(undef, ns[i]+1) for i = 1:nsurf]
    Vw = [Matrix{SVector{3,TF}}(undef, nw[i]+1, ns[i]+1) for i = 1:nsurf]
    dw = Tuple(Vector{TF}(undef, N) for i = 1:5)
    dΓ = Tuple(Vector{TF}(undef, N) for i = 1:5)
    dproperties = Tuple([Matrix{PanelProperties{TF}}(undef, nc[i], ns[i]) for i = 1:length(surfaces)] for j = 1:5)
    reference = Reference{TF}(1.0, 1.0, 1.0, [0.0, 0.0, 0.0], 1.0)
    freestream = Freestream{TF}(1.0, 0.0, 0.0, [0.0, 0.0, 0.0])
    # set default variables/parameters
    Γ .= 0.0
    Γdot .= 0.0
    symmetric = fill(false, nsurf)
    surface_id = collect(1:nsurf)
    wake_finite_core = fill(true, nsurf)
    iwake = fill(0, nsurf)
    nwake = nw
    trailing_vortices = fill(true, nsurf)
    xhat = [1.0, 0.0, 0.0]
    # construct system
    system = UnsteadySystem(AIC, w, Γ, Γdot, surfaces, repeated_points,
        wake_shedding_locations, wakes, properties, trefftz, Vcp, Vh, Vv, Vte,
        Vw, dw, dΓ, dproperties, reference, freestream, symmetric, surface_id,
        wake_finite_core, iwake, nwake, trailing_vortices, xhat)
    # update inputs to correspond to provided inputs
    update_inputs!(system, fcore, false; kwargs...)
    # return result
    return system
end

"""
    update_inputs!(system::UnsteadySystem, fcore=(c, Δs) -> 1e-3, preserve_core_size=true; kwargs...)

Updates the inputs in `system` to correspond to the provided keyword arguments.

# Arguments
 - `system`: Object that holds unsteady vortex lattice method inputs and state
 - `fcore`: Function for setting the finite core size when generating surface
    panels from grid inputs. Defaults to `(c, Δs) -> 1e-3`, where `c` is the
    cross-section chord length and `Δs` is the panel width.
 - `preserve_core_size`: Flag indicating whether the finite core size should be
    preserved from the previous set of panels (if applicable). Defaults to `true`.

# Keyword Arguments
 - `surfaces`: Surface geometry at the beginning of the simulation.
    Represented as either a:
     - Vector of grids of shape (3, nc+1, ns+1) which represent lifting surfaces
    or
     - Vector of matrices of shape (nc, ns) containing surface panels (see [`SurfacePanel`](@ref))
       where `nc` is the number of chordwise panels and `ns` is the number of
       spanwise panels.
 - `wakes`: Vector of initial wakes corresponding to each surface, represented
    by matrices of wake panels (see [`WakePanel`](@ref)) of shape (nw, ns) where
    `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels.
 - `circulation`: Vector containing the initial circulation of all surface
    panels in the system.
 - `circulation_rates`: Vector containing the initial circulation rates of all
    surface panels in the system.
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `freestream`: Freestream parameters (see [`Freestream`](@ref))
 - `symmetric`: Flags indicating whether the geometry for each surface is
    symmetric.  Note that applying symmetry to surfaces is only theoretically
    valid when the freestream conditions are symmetric as well.
 - `surface_id`: Surface ID for each surface.
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.
 - `iwake`: Initial number of chordwise wake panels to use for each surface.
 - `nwake`: Current maximum number of chordwise wake panels to use for each surface.
 - `trailing_vortices`: Flags indicating whether trailing vortices are shed from
    each surface.
 - `xhat`: Direction in which trailing vortices are shed.
"""
function update_inputs!(system::UnsteadySystem, fcore = (c, Δs) -> 1e-3, preserve_core_size=true;
    surfaces = nothing, wakes = nothing, circulation = nothing, reference = nothing,
    freestream = nothing, symmetric = nothing, surface_id = nothing,
    wake_finite_core = nothing, iwake = nothing, nwake = nothing,
    trailing_vortices = nothing, xhat = nothing)

    # number of surfaces
    nsurf = length(system.surfaces)

    if !isnothing(surfaces)
        # update surfaces
        update_surface_panels!(system, surfaces; fcore, preserve_core_size)
        # re-initialize repeated trailing edge points
        system.repeated_points = repeated_trailing_edge_points(surfaces)
        # re-initialize wake shedding locations
        initialize_wake_shedding_locations!(system, surfaces)
        # re-initialize undefined wake panels
        initialize_wake_panels!(system)
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

    # update reference parameters (if specified)
    if !isnothing(reference)
        system.reference = reference
    end

    # update freestream parameters (if specified)
    if !isnothing(freestream)
        system.freestream = freestream
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
    if !isnothing(nwake)
        copyto!(system.nwake, nwake)
    end
    if !isnothing(trailing_vortices)
        copyto!(system.trailing_vortices, trailing_vortices)
    end
    if !isnothing(xhat)
        copyto!(system.xhat, xhat)
    end

    return system
end

get_surface_circulation(system) = system.Gamma
get_surface_circulation_rate(system) = system.Gammadot
get_wake_shedding_locations(system) = system.wake_shedding_locations
get_wake_circulation(system) = get_circulation(system.wakes)
get_wake_circulation!(Γw, system) = get_circulation!(Γw, system.wakes)
get_wake_filament_length(system) = vortex_filament_length(system.wakes)
get_wake_filament_length!(l, system) = vortex_filament_length(l, system.wakes)
get_wake_vertices(system; iwake=size.(system.wakes, 1)) = get_vertices(system.wake_shedding_locations, system.wakes; iwake)
get_wake_vertices!(ζw, system; iwake=size.(system.wakes, 1)) = get_vertices!(ζw, system.wake_shedding_locations, system.wakes; iwake)
get_wake_velocities(system) = system.Vw

set_surface_circulation!(system, Γ) = system.Gamma .= Γ
set_wake_shedding_locations!(system, ζw) = get_wake_shedding_locations!(system.wake_shedding_locations, ζw)
set_wake_vertices!(system::AbstractSystem, ζw) = update_wake_panels!(system.wakes, ζw; helmholtz=false)
set_wake_circulation!(system, Γw) = set_circulation!(system.wakes, Γw)

set_wake_vertices!(wakes, ζw) = update_wake_panels!(wakes, ζw; helmholtz=false)


"""
    get_surface_properties(system::AbstractSystem)

Return a vector of surface panel properties for each surface, stored as matrices
of panel properties (see [`PanelProperties`](@ref)) of shape (nc, ns) where `nc`
is the number of chordwise panels and `ns` is the number of spanwise panels
"""
get_surface_properties(system::AbstractSystem)
get_surface_properties(system::SteadySystem) = system.properties
get_surface_properties(system::UnsteadySystem) = system.properties
