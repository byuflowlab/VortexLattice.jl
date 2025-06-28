"""
    PanelProperties

Panel specific properties calculated during the vortex lattice method analysis.

**Fields**
 - `gamma`: Vortex ring circulation strength, normalized by the reference velocity
 - `velocity`: Local velocity at the panel's bound vortex center, normalized by
    the reference velocity
 - `cfb`: Net force on the panel's bound vortex, as calculated using the
    Kutta-Joukowski theorem, normalized by the reference dynamic pressure and area
 - `cfl`: Force on the left bound vortex from this panel's vortex ring, as
    calculated by the Kutta-Joukowski theorem, normalized by the reference
    dynamic pressure and area
 - `cfr`: Force on the right bound vortex from this panel's vortex ring, as
    calculated by the Kutta-Joukowski theorem, normalized by the reference
    dynamic pressure and area
 - `velocity_from_streamwise`: Local velocity at the panel's bound vortex center
    excluding the induced velocity from the bound vorticies, not normalized
"""
struct PanelProperties{TF}
    gamma::TF
    velocity::SVector{3, TF}
    cfb::SVector{3, TF}
    cfl::SVector{3, TF}
    cfr::SVector{3, TF}
    velocity_from_streamwise::SVector{3, TF}
end

# constructor
function PanelProperties(gamma, velocity, cfb, cfl, cfr, velocity_from_streamwise)

    TF = promote_type(typeof(gamma), typeof(velocity), eltype(cfb), eltype(cfl),
        eltype(cfr), eltype(velocity_from_streamwise))

    return PanelProperties{TF}(gamma, velocity, cfb, cfl, cfr, velocity_from_streamwise)
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
 - `grids`: Grids of the surfaces, represented by matrices of vertices
 - `ratios`: Ratios of the locations of each control point on each panel
 - `surfaces`: Surfaces, represented by matrices of surface panels
 - `invert_normals`: Flags indicating whether the normals of each surface should
        be inverted (used for the nonlinear VLM)
 - `sections`: Section properties for each surface (Used as part of nonlinear VLM)
 - `properties`: Surface panel properties for each surface
 - `wakes`: Wake panel properties for each surface
 - `trefftz`: Trefftz panels associated with each surface
 - `reference`: Pointer to reference parameters associated with the system (see [`Reference`](@ref))
 - `freestream`: Pointer to current freestream parameters associated with the system (see [`Freestream`](@ref))
 - `symmetric`: Flags indicating whether each surface is symmetric across the X-Z plane
 - `nwake`: Number of chordwise wake panels to use from each wake in `wakes`,
 - `surface_id`: Surface ID for each surface.  The finite core model is disabled
    when calculating the influence of surfaces/wakes that share the same ID.
 - `trailing_vortices`: Flags to enable/disable trailing vortices
 - `wake_finite_core`: Flag for each wake indicating whether the finite core
    model should be enabled when calculating the wake's influence on itself and
    surfaces/wakes with the same surface ID.
 - `near_field_analysis`: Flag indicating whether a near field analysis has been
    performed for the current system state
 - `derivatives`: Flag indicating whether the derivatives with respect to the
    freestream variables have been calculated
 - `dw`: Derivatives of the R.H.S. with respect to the freestream variables
 - `dΓ`: Derivatives of the circulation strength with respect to the freestream variables
 - `dproperties`: Derivatives of the panel properties with respect to the freestream variables
 - `wake_shedding_locations`: Wake shedding locations for each surface
 - `Vcp`: Velocity due to surface motion at the control points
 - `Vh`: Velocity due to surface motion at the horizontal bound vortex centers
 - `Vv`: Velocity due to surface motion at the vertical bound vortex centers
 - `Vte`: Velocity due to surface motion at the trailing edge vertices
 - `dΓdt`: Derivative of the circulation strength with respect to non-dimensional time
"""
struct System{TF}
    AIC::Matrix{TF}
    w::Vector{TF}
    Γ::Vector{TF}
    V::Vector{Matrix{SVector{3,TF}}}
    grids::Vector{<:AbstractArray{TF, 3}}
    ratios::Vector{Array{TF,3}}
    surfaces::Vector{Matrix{SurfacePanel{TF}}}
    invert_normals::Vector{Bool}
    sections::Vector{Vector{SectionProperties{TF}}}
    properties::Vector{Matrix{PanelProperties{TF}}}
    wakes::Vector{Matrix{WakePanel{TF}}}
    trefftz::Vector{Vector{TrefftzPanel{TF}}}
    reference::Array{Reference{TF}, 0}
    freestream::Array{Freestream{TF}, 0}
    symmetric::Vector{Bool}
    nwake::Vector{Int}
    surface_id::Vector{Int}
    wake_finite_core::Vector{Bool}
    trailing_vortices::Vector{Bool}
    xhat::Array{SVector{3, TF}, 0}
    near_field_analysis::Array{Bool, 0}
    derivatives::Array{Bool, 0}
    dw::NTuple{5, Vector{TF}}
    dΓ::NTuple{5, Vector{TF}}
    dproperties::NTuple{5, Vector{Matrix{PanelProperties{TF}}}}
    wake_shedding_locations::Vector{Vector{SVector{3,TF}}}
    previous_surfaces::Vector{Matrix{SurfacePanel{TF}}}
    Vcp::Vector{Matrix{SVector{3, TF}}}
    Vh::Vector{Matrix{SVector{3, TF}}}
    Vv::Vector{Matrix{SVector{3, TF}}}
    Vte::Vector{Vector{SVector{3, TF}}}
    dΓdt::Vector{TF}
    lifting_line_r::Vector{Matrix{TF}}
    lifting_line_c::Vector{Vector{TF}}
end

Base.eltype(::Type{System{TF}}) where TF = TF
Base.eltype(::System{TF}) where TF = TF

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

# Keyword Arguments
 - `nw`: Number of chordwise wake panels for each surface. Defaults to zero wake
    panels on each surface
 - `grids`: Grids of the surfaces, represented by matrices of vertices
 - `ratios`: Ratios of the locations of each control point on each panel
 - `sections`: Section properties for each surface (Used as part of nonlinear VLM)
 - `invert_normals`: Flags indicating whether the normals of each surface should
    be inverted
"""
System(args...; kwargs...)

# grid inputs, no provided type
function System(grids::AbstractVector{<:AbstractArray{TF, 3}}; kwargs...) where TF

    return System(TF, grids; kwargs...)
end

# grid inputs, provided type
function System(TF::Type, grids::AbstractVector{<:AbstractArray{<:Any, 3}}; kwargs...)

    nc = size.(grids, 2) .- 1
    ns = size.(grids, 3) .- 1

    return System(TF, nc, ns; grids=grids, kwargs...)
end

# surface inputs, no provided type
function System(surfaces::AbstractVector{<:AbstractMatrix{SurfacePanel{TF}}};
    kwargs...) where TF

    return System(TF, surfaces; kwargs...)
end

# surface inputs, provided type
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
 - `grids`: Grids of the surfaces, represented by matrices of vertices
 - `ratios`: Ratios of the locations of each control point on each panel
 - `sections`: Section properties for each surface (Used as part of nonlinear VLM)
 - `invert_normals`: Flags indicating whether the normals of each surface should
    be inverted
"""
function System(TF::Type, nc, ns; nw = zero(nc), grids = nothing, ratios = nothing, sections = nothing, invert_normals = nothing)

    @assert length(nc) == length(ns) == length(nw)

    # number of surfaces
    nsurf = length(nc)

    # number of surface panels
    N = sum(nc .* ns)

    if isnothing(grids)
        grids = [Array{TF,3}(undef, 3, nc[i]+1, ns[i]+1) for i = 1:nsurf]
    end
    if isnothing(ratios)
        ratios = [Array{TF}(undef, 2, nc[i], ns[i]) for i = 1:nsurf]
        for i = 1:nsurf
            ratios[i] = ratios[i] .+ [0.5;0.75]
        end
    end

    if isnothing(sections)
        sections = [Vector{SectionProperties{TF}}() for i = 1:nsurf]
    else
        redefine_gamma_index!(sections, ns, nc)
    end

    if isnothing(invert_normals)
        invert_normals = fill(false, nsurf)
    end

    AIC = zeros(TF, N, N)
    w = zeros(TF, N)
    Γ = zeros(TF, N)
    V = [fill((@SVector zeros(TF, 3)), nw[i]+1, ns[i]+1) for i = 1:nsurf]
    surfaces = [Matrix{SurfacePanel{TF}}(undef, nc[i], ns[i]) for i = 1:nsurf]
    properties = [Matrix{PanelProperties{TF}}(undef, nc[i], ns[i]) for i = 1:nsurf]
    wakes = [Matrix{WakePanel{TF}}(undef, nw[i], ns[i]) for i = 1:nsurf]
    trefftz = [Vector{TrefftzPanel{TF}}(undef, ns[i]) for i = 1:nsurf]
    reference = Array{Reference{TF}}(undef)
    freestream = Array{Freestream{TF}}(undef)
    symmetric = [false for i = 1:nsurf]
    nwake = [0 for i = 1:nsurf]
    surface_id = [i for i = 1:nsurf]
    wake_finite_core = [true for i = 1:nsurf]
    trailing_vortices = [false for i = 1:nsurf]
    xhat = fill(SVector{3,TF}(1, 0, 0))
    near_field_analysis = fill(false)
    derivatives = fill(false)
    dw = Tuple(zeros(TF, N) for i = 1:5)
    dΓ = Tuple(zeros(TF, N) for i = 1:5)
    dproperties = Tuple([Matrix{PanelProperties{TF}}(undef, nc[i], ns[i]) for i = 1:length(surfaces)] for j = 1:5)
    wake_shedding_locations = [fill((@SVector zeros(TF, 3)), ns[i]+1) for i = 1:nsurf]
    previous_surfaces = [Matrix{SurfacePanel{TF}}(undef, nc[i], ns[i]) for i = 1:nsurf]
    Vcp = [fill((@SVector zeros(TF, 3)), nc[i], ns[i]) for i = 1:nsurf]
    Vh = [fill((@SVector zeros(TF, 3)), nc[i]+1, ns[i]) for i = 1:nsurf]
    Vv = [fill((@SVector zeros(TF, 3)), nc[i], ns[i]+1) for i = 1:nsurf]
    Vte = [fill((@SVector zeros(TF, 3)), ns[i]+1) for i = 1:nsurf]
    dΓdt = zeros(TF, N)

    return System{TF}(AIC, w, Γ, V, grids, ratios, surfaces, invert_normals, sections, 
        properties, wakes, trefftz, reference, freestream, symmetric, nwake, surface_id, 
        wake_finite_core, trailing_vortices, xhat, near_field_analysis, derivatives,
        dw, dΓ, dproperties, wake_shedding_locations, previous_surfaces, Vcp, Vh,
        Vv, Vte, dΓdt)
end

"""
    get_surface_properties(system)

Return a vector of surface panel properties for each surface, stored as matrices
of panel properties (see [`PanelProperties`](@ref)) of shape (nc, ns) where `nc`
is the number of chordwise panels and `ns` is the number of spanwise panels
"""
get_surface_properties(system) = system.properties


"""
    save_system_to_bson(system, filename)
Save the system to a BSON file.
# Arguments:
 - `system`: The system to save
 - `filename`: The name of the file to save the system to
# Returns:
 - `nothing`: The function does not return anything, it saves the system to a file
"""
function save_system_to_bson(system::System, filename::AbstractString)
    JLD2.save_object(filename, system)
    return nothing
end


"""
    load_system_from_bson(filename::AbstractString)
Load a system from a BSON file.
# Arguments:
 - `filename`: The name of the file to load the system from
# Returns:
 - `system`: The system loaded from the file
"""
function load_system_from_bson(filename::AbstractString)
    system = JLD2.load_object(filename)
    return system
end