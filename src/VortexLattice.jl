module VortexLattice

using LinearAlgebra
using StaticArrays
using FLOWMath
using WriteVTK
using VSPGeom
using CCBlade
using FastMultipole
using FLOWVPM
using DelimitedFiles
using ReadVTK

# value for dimensionalizing, included just for clarity in the algorithms
const RHO = 1.0

include("panel.jl")
export SurfacePanel, WakePanel, TrefftzPanel
export reflect, set_normal

include("wake.jl")
export Wake
export PanelParticleWake
export WakeSheddingMethod, NoShed, SigmaPPS, OverlapPPS, SigmaOverlap
export IntegrationScheme, EulerScheme, RK3Scheme

include("panel_particle_wake.jl")

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
export Freestream, trajectory_to_freestream

include("induced.jl")

include("circulation.jl")

include("system.jl")
export System
export PanelProperties, get_surface_properties

include("analyses.jl")
export steady_analysis, steady_analysis!
export unsteady_analysis, unsteady_analysis!
export spanwise_force_coefficients

include("nearfield.jl")
export body_forces
export body_forces_history
export lifting_line_coefficients, lifting_line_coefficients!

include("farfield.jl")
export far_field_drag

include("stability.jl")
export body_derivatives, stability_derivatives

include("frames.jl")
export ReferenceFrame, propagate_kinematics!, change_convention!, BackRightUp, ForwardRightDown

include("rotors.jl")
export generate_rotor

include("visualization.jl")
export write_vtk
export write_restart_checkpoint, read_restart_checkpoint_info, restore_restart!

include("fmm.jl")
export fmm

include("fluid_domain.jl")
export FluidDomainMonitor
export evaluate_fluid_domain!, evaluate_fluid_domain_from_restarts!

include("viscous.jl")

include("unsteady.jl")
export simulate!

const DEBUG = Array{Bool,0}(undef)
DEBUG[] = false

end # module
