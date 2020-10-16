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
VLMFlow.linearinterp
VLMFlow.spanwise_spacing
VLMFlow.chordwise_spacing
VLMFlow.interpolate_grid
```

### Vortex Lattice Panels
```@docs
VLMFlow.controlpoint
VLMFlow.normal
VLMFlow.trefftz_normal
VLMFlow.midpoint
VLMFlow.left_midpoint
VLMFlow.right_midpoint
VLMFlow.top_vector
VLMFlow.left_vector
VLMFlow.right_vector
VLMFlow.endpoints
VLMFlow.reflected_endpoints
VLMFlow.trefftz_endpoints
VLMFlow.trefftz_center
VLMFlow.panel_induced_drag
VLMFlow.flipy
VLMFlow.not_on_symmetry_plane
VLMFlow.trailing_induced_velocity
VLMFlow.bound_induced_velocity
VLMFlow.induced_velocity
VLMFlow.vortex_induced_drag
```

### Freestream
```@docs
VLMFlow.body_to_stability_alpha
VLMFlow.body_to_wind_derivatives
VLMFlow.stability_to_body_alpha
VLMFlow.stability_to_wind_beta
VLMFlow.wind_to_body_derivatives
VLMFlow.wind_to_stability_beta
VLMFlow.freestream_velocity
VLMFlow.freestream_velocity_derivatives
VLMFlow.external_velocity
VLMFlow.external_velocity_derivatives
```

### Circulation
```@docs
VLMFlow.normal_velocity_derivatives
VLMFlow.normal_velocity_derivatives!
VLMFlow.circulation_derivatives
```

### Near-Field
```@docs
VLMFlow.near_field_forces_derivatives
VLMFlow.body_to_frame
```

### Far-Field
```@docs
VLMFlow.project_panels
```

## Index

```@index
```
