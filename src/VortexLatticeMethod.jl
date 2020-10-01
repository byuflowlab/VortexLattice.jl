module VortexLatticeMethod

using LinearAlgebra
using StaticArrays
using Interpolations
using ForwardDiff

# values for dimensionalizing, included just for clarity in the algorithms
const RHO = 1.0
const VINF = 1.0
const QINF = 0.5*RHO*VINF^2

include("panel.jl")
export AbstractPanel, Horseshoe, Ring

include("geometry.jl")
export Uniform, Cosine
export grid_to_horseshoe_vortices, grid_to_vortex_rings
export wing_to_horseshoe_vortices, wing_to_vortex_rings

include("reference.jl")
export Reference

include("freestream.jl")
export Freestream

include("circulation.jl")
export influence_coefficients, influence_coefficients!
export normal_velocity, normal_velocity!
export circulation, circulation!

include("forces.jl")
export Outputs
export forces_moments, trefftz_induced_drag, vlm

include("stability.jl")
export StabilityDerivatives, stability_analysis

end # module
