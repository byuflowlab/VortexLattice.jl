"""
    PanelProperties

Panel specific properties calculated during the vortex lattice method analysis.

**Fields**
 - `gamma`: Vortex ring circulation strength, normalized by VINF
 - `velocity`: Local velocity at the panel's bound vortex center, normalized by VINF
 - `cfb`: Net force on the panel's bound vortex, as calculated using the
    Kutta-Joukowski theorem, normalized by `QINF*S`
 - `cfl`: Force on the left bound vortex from this panel's vortex ring, as
    calculated by the Kutta-Joukowski theorem, normalized by `QINF*S`
 - `cfr`: Force on the right bound vortex from this panel's vortex ring, as
    calculated by the Kutta-Joukowski theorem, normalized by `QINF*S`
"""
struct PanelProperties{TF}
    gamma::TF
    velocity::SVector{3, TF}
    cfb::SVector{3, TF}
    cfl::SVector{3, TF}
    cfr::SVector{3, TF}
end

# constructor
function PanelProperties(gamma, velocity, cfb, cfl, cfr)

    TF = promote_type(typeof(gamma), typeof(velocity), eltype(cfb), eltype(cfl),
        eltype(cfr))

    return PanelProperties{TF}(gamma, velocity, cfb, cfl, cfr)
end

Base.eltype(::Type{PanelProperties{TF}}) where TF = TF
Base.eltype(::PanelProperties{TF}) where TF = TF

"""
    System{TF}

Contains pre-allocated storage for internal system variables.

# Fields:
 - `AIC`: Aerodynamic influence coefficient matrix from the surface panels
 - `w`: Normal velocity at the control points from external sources and wakes
 - `Γ`: Circulation strength of the surface panels
 - `V`: Velocity at the wake vertices for each surface
 - `surfaces`: Surfaces, represented by matrices of surface panels
 - `properties`: Surface panel properties for each surface
 - `wakes`: Wake panel properties for each surface
 - `trefftz`: Trefftz panels associated with each surface
 - `dw`: Derivatives of the R.H.S. with respect to the freestream variables
 - `dΓ`: Derivatives of the circulation strength with respect to the freestream variables
 - `dproperties`: Derivatives of the panel properties with respect to the freestream variables
 - `wake_shedding_locations`: Wake shedding locations for each surface
 - `dΓdt`: Derivative of the circulation strength with respect to non-dimensional time
"""
struct System{TF}
    AIC::Matrix{TF}
    w::Vector{TF}
    Γ::Vector{TF}
    V::Vector{Matrix{SVector{3,TF}}}
    surfaces::Vector{Matrix{SurfacePanel{TF}}}
    properties::Vector{Matrix{PanelProperties{TF}}}
    wakes::Vector{Matrix{WakePanel{TF}}}
    trefftz::Vector{Vector{TrefftzPanel{TF}}}
    dw::NTuple{5, Vector{TF}}
    dΓ::NTuple{5, Vector{TF}}
    dproperties::NTuple{5, Vector{Matrix{PanelProperties{TF}}}}
    wake_shedding_locations::Vector{Vector{SVector{3,TF}}}
    dΓdt::Vector{TF}
end

@inline Base.eltype(::Type{System{TF}}) where TF = TF
@inline Base.eltype(::System{TF}) where TF = TF

"""
    System([TF], surfaces; kwargs...)

Return an object of type `System` with pre-allocated storage for internal system
variables

# Arguments:
 - `TF`: Floating point type, defaults to the floating point type used by `surface`
 - `surfaces`:
        Either:
         - One or more grids of shape (3, nc+1, ns+1) which represents lifting surfaces,
         or
         - One or more matrices of surface panels (see [`SurfacePanel`](@ref)) of
           shape (nc, ns)
        where `nc` is the number of chordwise panels and `ns` is the number of
        spanwise panels

# Keyword Arguments:
 - `nw`: Number of chordwise wake panels to initialize for each surface. Defaults to
    zero wake panels for each surface.
"""
function System() end

# one grid, no provided type
function System(grid::AbstractMatrix{TF}; kwargs...) where TF

    return System(TF, grid; kwargs...)
end

# one grid, provided type
function System(TF::Type, grid::AbstractMatrix; kwargs...)

    nc = size(grid, 2)
    ns = size(grid, 3)

    return System(TF, nc, ns; kwargs...)
end

# one surface, no provided type
function System(surface::AbstractMatrix{SurfacePanel{TF}}; kwargs...) where TF

    return System(TF, surface; kwargs...)
end

# one surface, provided type
function System(TF::Type, surface::AbstractMatrix{<:SurfacePanel}; kwargs...)

    nc = size(surface, 1)
    ns = size(surface, 2)

    return System(TF, nc, ns; kwargs...)
end

# multiple grids, no provided type
function System(grids::AbstractVector{<:AbstractMatrix{TF}}; kwargs...) where TF

    return System(TF, grids; kwargs...)
end

# multiple grids, provided type
function System(TF::Type, grids::AbstractVector{<:AbstractMatrix}; kwargs...)

    nc = size(grids, 2)
    ns = size(grids, 3)

    return System(TF, nc, ns; kwargs...)
end

# multiple surfaces, no provided type
function System(surfaces::AbstractVector{<:AbstractMatrix{SurfacePanel{TF}}};
    kwargs...) where TF

    return System(TF, surfaces; kwargs...)
end

# multiple surfaces, provided type
function System(TF::Type, surfaces::AbstractVector{<:AbstractMatrix{<:SurfacePanel}};
    kwargs...)

    nc = size.(surfaces, 1)
    ns = size.(surfaces, 2)

    return System(TF, nc, ns; kwargs...)
end

"""
    System(TF, nc, ns; nw = zero(nc))

Return an object of type `System` with pre-allocated storage for internal system
variables

# Arguments:
 - `TF`: Floating point type.
 - `nc`: Number of chordwise panels on each surface.
 - `ns`: Number of spanwise panels on each surface.

# Keyword Arguments
 - `nw`: Number of chordwise wake panels for each surface. Defaults to zero wake
    panels on each surface
"""
function System(TF::Type, nc, ns; nw = zero(nc))

    @assert length(nc) == length(ns) == length(nw)

    # number of surfaces
    nsurf = length(nc)

    # number of surface panels
    N = sum(nc .* ns)

    AIC = zeros(TF, N, N)
    w = zeros(TF, N)
    Γ = zeros(TF, N)
    V = [Matrix{SVector{3, TF}}(undef, nw[i]+1, ns[i]+1) for i = 1:nsurf]
    surfaces = [Matrix{SurfacePanel{TF}}(undef, nc[i], ns[i]) for i = 1:nsurf]
    properties = [Matrix{PanelProperties{TF}}(undef, nc[i], ns[i]) for i = 1:nsurf]
    wakes = [Matrix{WakePanel{TF}}(undef, nw[i], ns[i]) for i = 1:nsurf]
    trefftz = [Vector{TrefftzPanel{TF}}(undef, ns[i]) for i = 1:nsurf]
    dw = Tuple(zeros(TF, N) for i = 1:5)
    dΓ = Tuple(zeros(TF, N) for i = 1:5)
    dproperties = Tuple([Matrix{PanelProperties{TF}}(undef, nc[i], ns[i]) for i = 1:length(surfaces)] for j = 1:5)
    wake_shedding_locations = [Vector{SVector{3, TF}}(undef, ns[i]+1) for i = 1:nsurf]
    dΓdt = zeros(TF, N)

    return System{TF}(AIC, w, Γ, V, surfaces, properties, wakes, trefftz, dw, dΓ,
        dproperties, wake_shedding_locations, dΓdt)
end

"""
    get_panel_properties(system, surface)

Return the panel properties stored in `system`
"""
get_panel_properties(system, surface::AbstractMatrix) = system.properties[1]

"""
    get_panel_properties(system, surfaces)

Return the panel properties stored in `system`
"""
get_panel_properties(system, surfaces::AbstractVector{<:AbstractMatrix}) = system.properties
