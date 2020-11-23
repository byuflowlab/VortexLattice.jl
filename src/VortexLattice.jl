module VortexLattice

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
export translate, translate!, reflect

include("geometry.jl")
export AbstractSpacing, Uniform, Sine, Cosine
export grid_to_horseshoe_vortices, grid_to_vortex_rings
export wing_to_horseshoe_vortices, wing_to_vortex_rings

include("reference.jl")
export Reference, AbstractFrame, Body, Stability, Wind

include("freestream.jl")
export Freestream
export body_to_stability, stability_to_body
export stability_to_wind, wind_to_stability
export wind_to_body, body_to_wind

include("wake.jl")
export Wake

include("farfield.jl")
export TrefftzPanel
export far_field_drag

include("circulation.jl")

include("nearfield.jl")
export body_forces

include("stability.jl")
export body_derivatives, stability_derivatives

include("system.jl")
export PanelProperties
export System

include("analyses.jl")
export steady_analysis, steady_analysis!

include("visualization.jl")
export write_vtk

end # module
