module VortexLattice

using LinearAlgebra
using StaticArrays
using Interpolations
using WriteVTK

# values for dimensionalizing, included just for clarity in the algorithms
const RHO = 1.0
const VINF = 1.0
const QINF = 0.5*RHO*VINF^2

# useful functions for dealing with tuples
tuple_len(::NTuple{N,Any}) where N = N
tuple_len(::Type{<:NTuple{N,Any}}) where N = N
tuple_type_param(t::Tuple, ::Val{i}) where {T <: Tuple, i} = typeof(t[i])
tuple_type_param(::Type{T}, ::Val{i}) where {T <: Tuple, i} = T.parameters[i]

include("panel.jl")
export VortexRing
export translate, translate!, reflect

include("wake.jl")
export Wake

include("geometry.jl")
export AbstractSpacing, Uniform, Sine, Cosine
export grid_to_vortex_rings, wing_to_vortex_rings
export translate, translate!, reflect

include("reference.jl")
export Reference
export AbstractFrame, Body, Stability, Wind

include("freestream.jl")
export Freestream

include("induced.jl")

include("circulation.jl")

include("system.jl")
export System
export PanelProperties, panel_properties
export panel_properties

include("analyses.jl")
export steady_analysis, steady_analysis!
export unsteady_analysis, unsteady_analysis!

include("nearfield.jl")
export body_forces

include("farfield.jl")
export far_field_drag

include("stability.jl")
export body_derivatives, stability_derivatives

include("visualization.jl")
export write_vtk

end # module
