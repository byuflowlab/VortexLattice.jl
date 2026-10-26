module VortexLattice

using LinearAlgebra
using StaticArrays
using FLOWMath
using WriteVTK
using VSPGeom
using CCBlade
import SixDOF
import FLOWUnsteadyCore
using FLOWUnsteadyCore: Freestream
import SixDOF: ReferenceFrame, propagate_kinematics!, add_frame!, Rodrigues, inverse_Rodrigues

# value for dimensionalizing, included just for clarity in the algorithms
const RHO = 1.0

include("panel.jl")
export SurfacePanel, WakePanel, TrefftzPanel
export reflect, set_normal


include("geometry.jl")
export AbstractSpacing, Uniform, Sine, Cosine
export grid_to_surface_panels, wing_to_grid
export lifting_line_geometry, lifting_line_geometry!
export translate, translate!, rotate, rotate!, reflect

include("vspgeom.jl")
export read_degengeom, import_vsp

include("reference.jl")
export Reference
export AbstractFrame, Body, Stability, Wind

include("freestream.jl")
export Freestream   # the core's, re-exported

include("induced.jl")

include("circulation.jl")

include("system.jl")
export VortexLatticeSystem
export PanelProperties, get_surface_properties

include("analyses.jl")
export steady_analysis, steady_analysis!
export spanwise_force_coefficients

include("nearfield.jl")
export body_forces
export lifting_line_coefficients, lifting_line_coefficients!

include("farfield.jl")
export far_field_drag

include("stability.jl")
export body_derivatives, stability_derivatives

include("frames.jl")
export ReferenceFrame, propagate_kinematics!, change_convention!, BackRightUp, ForwardRightDown
export flight_condition

include("visualization.jl")
export write_vtk

include("viscous.jl")

end # module
