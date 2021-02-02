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
export SurfacePanel, WakePanel, TrefftzPanel
export translate, translate!, reflect, set_normal

include("wake.jl")
export Wake

include("geometry.jl")
export AbstractSpacing, Uniform, Sine, Cosine
export grid_to_surface_panels, wing_to_surface_panels
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
export PanelProperties, get_panel_properties

include("analyses.jl")
export steady_analysis, steady_analysis!
export unsteady_analysis, unsteady_analysis!
export prescribed_motion

include("nearfield.jl")
export body_forces
export body_forces_history

include("farfield.jl")
export far_field_drag

include("stability.jl")
export body_derivatives, stability_derivatives

include("visualization.jl")
export write_vtk

end # module
