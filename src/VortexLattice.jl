module VortexLattice

using LinearAlgebra
using StaticArrays
using SparseArrays
using Interpolations
using DifferentialEquations
using ForwardDiff
using WriteVTK

# value for dimensionalizing, included just for clarity in the algorithms
const RHO = 1.0

include("panel.jl")
export SurfacePanel, WakePanel, TrefftzPanel, PanelProperties
export translate, translate!, rotate, rotate!, reflect, set_normal

include("geometry.jl")
export AbstractSpacing, Uniform, Sine, Cosine
export grid_to_surface_panels, wing_to_surface_panels
export lifting_line_geometry, lifting_line_geometry!

include("reference.jl")
export Reference
export AbstractFrame, Body, Stability, Wind

include("freestream.jl")
export Freestream, trajectory_to_freestream

include("induced.jl")

include("circulation.jl")

include("nearfield.jl")

include("farfield.jl")

include("stability.jl")

include("wake.jl")

include("visualization.jl")
export write_vtk

include("system.jl")
export SteadySystem, UnsteadySystem
export get_surface_properties

include("interface.jl")
export body_forces, body_forces_history
export lifting_line_coefficients, lifting_line_forces
export far_field_drag
export body_derivatives, stability_derivatives

include("analyses.jl")
export steady_analysis, steady_analysis!
export unsteady_analysis, unsteady_analysis!

# DifferentialEquations.jl interface
include("interfaces/diffeq.jl")

end # module
