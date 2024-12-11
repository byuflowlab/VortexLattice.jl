# Library

```@contents
Pages = ["library.md"]
Depth = 3
```

## Public API

### Generating Lifting Surfaces

```@docs
AbstractSpacing
Uniform
Sine
Cosine
SurfacePanel
SurfacePanel()
WakePanel
WakePanel()
grid_to_surface_panels
wing_to_grid
lifting_line_geometry
lifting_line_geometry!
read_degengeom
import_vsp
set_normal
translate
translate!
rotate
rotate!
reflect(surface::AbstractMatrix)
```

### Reference Parameters and Frames
```@docs
Reference
AbstractFrame
Body
Stability
Wind
```

### Freestream Parameters
```@docs
Freestream
trajectory_to_freestream
```

### Performing an Analysis
```@docs
System
System()
steady_analysis
steady_analysis!
unsteady_analysis
unsteady_analysis!
```

### Near Field Forces and Moments
```@docs
PanelProperties
get_surface_properties
body_forces(system::System{TF}; frame = Body()) where TF
body_forces_history
lifting_line_coefficients
lifting_line_coefficients!
```

### Far Field Drag
```@docs
far_field_drag
```

### Body and Stability Derivatives
```@docs
body_derivatives
stability_derivatives
```

### Visualization
```@docs
write_vtk
```

## Private API

### Geometry
```@docs
VortexLattice.linearinterp
VortexLattice.spanwise_spacing
VortexLattice.chordwise_spacing
VortexLattice.interpolate_grid
VortexLattice.update_surface_panels!
VortexLattice.trailing_edge_points
VortexLattice.repeated_trailing_edge_points
VortexLattice.flipy
VortexLattice.on_symmetry_plane
VortexLattice.not_on_symmetry_plane
```

### Surface Panels
```@docs
VortexLattice.top_left
VortexLattice.top_center
VortexLattice.top_right
VortexLattice.bottom_left
VortexLattice.bottom_center
VortexLattice.bottom_right
VortexLattice.controlpoint
VortexLattice.normal(panel::SurfacePanel)
VortexLattice.get_core_size
VortexLattice.reflect(panel::SurfacePanel)
VortexLattice.left_center
VortexLattice.right_center
VortexLattice.top_vector
VortexLattice.left_vector
VortexLattice.right_vector
VortexLattice.bottom_vector
```

### Wake Panels
```@docs
VortexLattice.update_wake_shedding_locations!
VortexLattice.set_circulation_strength
VortexLattice.circulation_strength
VortexLattice.get_wake_velocities!
VortexLattice.translate_wake
VortexLattice.translate_wake!
VortexLattice.shed_wake!
VortexLattice.rowshift!
```

### Induced Velocity
```@docs
VortexLattice.bound_induced_velocity
VortexLattice.trailing_induced_velocity
VortexLattice.ring_induced_velocity
VortexLattice.influence_coefficients!
VortexLattice.update_trailing_edge_coefficients!
VortexLattice.induced_velocity
VortexLattice.induced_velocity_derivatives
```

### Freestream
```@docs
VortexLattice.body_to_stability
VortexLattice.body_to_wind
VortexLattice.stability_to_body
VortexLattice.stability_to_wind
VortexLattice.wind_to_body
VortexLattice.wind_to_stability
VortexLattice.body_to_stability_alpha
VortexLattice.body_to_wind_derivatives
VortexLattice.stability_to_body_alpha
VortexLattice.stability_to_wind_beta
VortexLattice.wind_to_body_derivatives
VortexLattice.wind_to_stability_beta
VortexLattice.freestream_velocity
VortexLattice.freestream_velocity_derivatives
VortexLattice.rotational_velocity
VortexLattice.rotational_velocity_derivatives
VortexLattice.get_surface_velocities!
```

### Circulation
```@docs
VortexLattice.normal_velocity!
VortexLattice.normal_velocity_derivatives!
VortexLattice.circulation
VortexLattice.circulation!
VortexLattice.circulation_derivatives
VortexLattice.circulation_derivatives!
```

### Time-Domain Analysis
```@docs
VortexLattice.propagate_system!
```

### Near-Field Analysis
```@docs
VortexLattice.near_field_forces!
VortexLattice.near_field_forces_derivatives!
VortexLattice.body_forces(surfaces, properties, ref, fs, symmetric, frame)
VortexLattice.body_forces_derivatives
VortexLattice.body_to_frame
```

### Far-Field
```@docs
VortexLattice.TrefftzPanel
VortexLattice.normal(panel::VortexLattice.TrefftzPanel)
VortexLattice.trefftz_panels
VortexLattice.trefftz_panel_induced_drag
VortexLattice.vortex_induced_drag
```

### Visualization
```@docs
VortexLattice.write_vtk!
```

## Index

```@index
```
