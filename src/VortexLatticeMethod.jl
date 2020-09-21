#=
Author: Andrew Ning

Vortex Lattice Method

- spanwise and chordwise panels
- additional Trefftz plane analysis for induced drag
- stability derivatives
- general inflow

=#

module VortexLatticeMethod

using StaticArrays, LinearAlgebra

export Panel, Freestream, Reference, Outputs
export vlm

# normalized (so they don't matter, but are included just for clarity in the algorithms)
const RHO = 1.0
const VINF = 1.0

include("vlm.jl")
include("geometry.jl")  # defines some convenience functions for generating geometry
include("sderiv.jl")  # defines stability derivative type and some associated methods




end # module
