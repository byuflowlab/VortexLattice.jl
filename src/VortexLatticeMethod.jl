#=
Author: Andrew Ning

Vortex Lattice Method

- spanwise and chordwise panels
- additional Trefftz plane analysis for induced drag
- stability derivatives
- general inflow

=#

module VortexLatticeMethod

using ForwardDiff, StaticArrays, LinearAlgebra

# core functionality
export Panel, Freestream, Reference, Outputs
export vlm, stability_analysis

# geometry convenience functions
export Uniform, Cosine
export create_grid, linear_sections, simple_wing, translate!

# included just for clarity in the algorithms
const RHO = 1.0
const VINF = 1.0

include("vlm.jl") # core functionality of this package
include("geometry.jl")  # convenience functions for generating panels

end # module
