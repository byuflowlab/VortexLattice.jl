# Library

```@contents
Pages = ["library.md"]
Depth = 3
```

## Public API

### Generating Geometry

```@docs
AbstractSpacing
Uniform
Sine
Cosine
grid_to_horseshoe_vortices
grid_to_vortex_rings
wing_to_horseshoe_vortices
wing_to_vortex_rings
```

### Manipulating Vortex Lattice Panels
```@docs
AbstractPanel
Horseshoe
Ring
translate
translate!
reflect
rotate
```

### Reference Parameters
```@docs
Reference
```

### Freestream Parameters
```@docs
Freestream
body_to_stability
body_to_wind
stability_to_body
stability_to_wind
wind_to_body
wind_to_stability
```

### Solving for Circulation
```@docs
influence_coefficients
influence_coefficients!
normal_velocity
normal_velocity!
circulation
circulation!
```

### Near Field Forces and Moments
```@docs
AbstractFrame
Body
Stability
Wind
PanelProperties
near_field_forces
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
VortexFlow.linearinterp
VortexFlow.spanwise_spacing
VortexFlow.chordwise_spacing
VortexFlow.interpolate_grid
```

### Vortex Lattice Panels
```@docs
VortexFlow.controlpoint
VortexFlow.normal
VortexFlow.trefftz_normal
VortexFlow.midpoint
VortexFlow.left_midpoint
VortexFlow.right_midpoint
VortexFlow.top_vector
VortexFlow.left_vector
VortexFlow.right_vector
VortexFlow.endpoints
VortexFlow.reflected_endpoints
VortexFlow.trefftz_endpoints
VortexFlow.trefftz_center
VortexFlow.panel_induced_drag
VortexFlow.flipy
VortexFlow.not_on_symmetry_plane
VortexFlow.trailing_induced_velocity
VortexFlow.bound_induced_velocity
VortexFlow.induced_velocity
VortexFlow.vortex_induced_drag
```

### Freestream
```@docs
VortexFlow.body_to_stability_alpha
VortexFlow.body_to_wind_derivatives
VortexFlow.stability_to_body_alpha
VortexFlow.stability_to_wind_beta
VortexFlow.wind_to_body_derivatives
VortexFlow.wind_to_stability_beta
VortexFlow.freestream_velocity
VortexFlow.freestream_velocity_derivatives
VortexFlow.external_velocity
VortexFlow.external_velocity_derivatives
```

### Circulation
```@docs
VortexFlow.normal_velocity_derivatives
VortexFlow.normal_velocity_derivatives!
VortexFlow.circulation_derivatives
```

### Near-Field
```@docs
VortexFlow.near_field_forces_derivatives
VortexFlow.body_to_frame
```

### Far-Field
```@docs
VortexFlow.project_panels
```

## Index

```@index
```
