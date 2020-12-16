# --- structure for holding panel forces --- #

"""
    PanelProperties

Panel specific properties calculated during the vortex lattice method analysis.

**Fields**
 - `gamma`: Net circulation strength at the bound vortex center, normalized by VINF
 - `v`: Local velocity at the bound vortex center, normalized by VINF
 - `cf`: Bound vortex force (acting at its center), normalized by `QINF*S`
 - `cfl`: Left vortex force (acting at its midpoint on the surface), normalized by `QINF*S`
 - `cfr`: Right vortex force (acting at its midpoint on the surface), normalized by `QINF*S`
"""
struct PanelProperties{TF}
    gamma::TF
    v::SVector{3, TF}
    cf::SVector{3, TF}
    cfl::SVector{3, TF}
    cfr::SVector{3, TF}
end

function PanelProperties(gamma, v, cf, cfl, cfr)

    TF = promote_type(typeof(gamma), typeof(v), eltype(cf), eltype(cfl), eltype(cfr))

    return PanelProperties{TF}(gamma, v, cf, cfl, cfr)
end

Base.eltype(::Type{PanelProperties{TF}}) where TF = TF
Base.eltype(::PanelProperties{TF}) where TF = TF

"""
    System{TF}

Contains the system AIC matrix, R.H.S., circulation distribution, and wake shape.

# Fields:
 - `Acb`: Aerodynamic influence coefficient matrix on the control points from the body
 - `Acw`: Aerodynamic influence coefficient matrix on the control points from the wake
 - `Avb`: Aerodynamic influence coefficient matrix on the wake vertices from the body
 - `Avw`: Aerodynamic influence coefficient matrix on the wake vertices from the wake
 - `Γb`: Circulation strength of the body panels
 - `Γw`: Circulation strength of the wake panels
 - `w`: Normal velocity at the control points (except induced velocities)
 - `V`: Velocity at the wake vertices

 - `b`: System R.H.S.
 - `gamma`: Circulation strength of each panel
 - `panels`: Panel properties for each surface
 - `wakes`: Wake panels for each surface
 - `trefftz`: Trefftz panels associated with each surface
 - `dgamma`: Derivatives of the R.H.S. with respect to the freestream variables
 - `dgamma`: Derivatives of the circulation strength with respect to the freestream variables
 - `dpanels`: Derivatives of the panel properties with respect to the freestream variables
 - `wake_velocities`: Velocities at the corners of the wake panels for each surface
"""
struct System{TF}
    AIC::Matrix{TF}
    b::Vector{TF}
    gamma::Vector{TF}
    panels::Vector{Matrix{PanelProperties{TF}}}
    wakes::Vector{Matrix{Wake{TF}}}
    trefftz::Vector{Vector{TrefftzPanel{TF}}}
    db::NTuple{5, Vector{TF}}
    dgamma::NTuple{5, Vector{TF}}
    dpanels::NTuple{5, Vector{Matrix{PanelProperties{TF}}}}
    wake_velocities::Vector{Matrix{SVector{3,TF}}}
end

@inline Base.eltype(::Type{System{TF}}) where TF = TF
@inline Base.eltype(::System{TF}) where TF = TF

"""
    System([TF], surface; nwake = 0)

Returns an object of type `System` with pre-allocated storage for vortex lattice
calculations.

# Arguments:
 - `TF`: Floating point type, defaults to the floating point type used by `surface`
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels

# Keyword Arguments:
 - `nwake = 0`: Number of chordwise wake panels in pre-allocated storage
"""
System(args...; kwargs...)

System(surface::AbstractMatrix; kwargs...) = System(eltype(eltype(surface)), surface; kwargs...)

function System(TF::Type{<:AbstractFloat}, surface::AbstractMatrix; nwake=0)

    N = length(surface)
    AIC = zeros(TF, N, N)
    b = zeros(TF, N)
    gamma = zeros(TF, N)
    panels = [Matrix{PanelProperties{TF}}(undef, size(surface))]
    wakes = [Matrix{Wake{TF}}(undef, nwake, size(surface, 2))]
    trefftz = [Vector{TrefftzPanel{TF}}(undef, size(surface, 2))]
    db = Tuple(zeros(TF, N) for id = 1:5)
    dgamma = Tuple(zeros(TF, N) for id = 1:5)
    dpanels = Tuple([Matrix{PanelProperties{TF}}(undef, size(surface))] for id = 1:5)
    wake_velocities = [Matrix{SVector{3, TF}}(undef, nwake+1, size(surface, 2)+1)]

    return System{TF}(AIC, b, gamma, panels, wakes, trefftz, db, dgamma, dpanels,
        wake_velocities)
end

"""
    System([TF], surfaces; nwake = 0)

Returns an object of type `System` with pre-allocated storage for vortex lattice
calculations.

# Arguments:
 - `TF`: Floating point type, defaults to the floating point type used by `surface`
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels

# Keyword Arguments:
 - `nwake = 0`: Number of chordwise wake panels in pre-allocated storage
"""
System(surfaces::AbstractVector{<:AbstractMatrix}; kwargs...) =
    System(eltype(eltype(eltype(surfaces))), surfaces; kwargs...)

function System(TF::Type{<:AbstractFloat}, surfaces::AbstractVector{<:AbstractMatrix};
    nwake = fill(0, length(surfaces)))

    N = sum(length.(surfaces))
    AIC = zeros(TF, N, N)
    b = zeros(TF, N)
    gamma = zeros(TF, N)
    panels = [Matrix{PanelProperties{TF}}(undef, size(surfaces[i])) for i = 1:length(surfaces)]
    wakes = [Matrix{Wake{TF}}(undef, nwake[i], size(surfaces[i], 2)) for i = 1:length(surfaces)]
    trefftz = [Vector{TrefftzPanel{TF}}(undef, size(surfaces[i], 2)) for i = 1:length(surfaces)]
    db = Tuple(zeros(TF, N) for id = 1:5)
    dgamma = Tuple(zeros(TF, N) for id = 1:5)
    dpanels = Tuple([Matrix{PanelProperties{TF}}(undef, size(surfaces[i])) for i = 1:length(surfaces)] for id = 1:5)
    wake_velocities = [Matrix{SVector{3, TF}}(undef, nwake[i]+1, size(surfaces[i], 2)+1) for i = 1:length(surfaces)]

    return System{TF}(AIC, b, gamma, panels, wakes, trefftz, db, dgamma, dpanels,
        wake_velocities)
end

"""
    panel_properties(system)

Returns panel properties resulting from a near field analysis
"""
panel_properties(system) = system.panels
