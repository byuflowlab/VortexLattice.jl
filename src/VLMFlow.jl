module VLMFlow

using LinearAlgebra
using StaticArrays
using Interpolations
using WriteVTK

# values for dimensionalizing, included just for clarity in the algorithms
const RHO = 1.0
const VINF = 1.0
const QINF = 0.5*RHO*VINF^2

include("panel.jl")
export AbstractPanel, Horseshoe, Ring
export translate, translate!, reflect, rotate

include("geometry.jl")
export AbstractSpacing, Uniform, Sine, Cosine
export grid_to_horseshoe_vortices, grid_to_vortex_rings
export wing_to_horseshoe_vortices, wing_to_vortex_rings

include("reference.jl")
export Reference

include("freestream.jl")
export Freestream
export body_to_stability, stability_to_body
export stability_to_wind, wind_to_stability
export wind_to_body, body_to_wind

include("circulation.jl")
export influence_coefficients, influence_coefficients!
export normal_velocity, normal_velocity!
export circulation

include("nearfield.jl")
export AbstractFrame, Body, Stability, Wind
export PanelProperties
export near_field_forces

include("farfield.jl")
export far_field_drag

include("stability.jl")
export body_derivatives, stability_derivatives

include("visualization.jl")
export write_vtk

end # module
