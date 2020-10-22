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
grid_to_horseshoe_vortices(xyz)
grid_to_vortex_rings(xyz)
grid_to_horseshoe_vortices(xyz, ns, nc)
grid_to_vortex_rings(xyz, ns, nc)
wing_to_horseshoe_vortices
wing_to_vortex_rings
```

### Manipulating Lifting Surfaces
```@docs
AbstractPanel
Horseshoe
Ring
translate(panels::AbstractVector{<:AbstractPanel}, r)
translate!(panels, r)
reflect(panels::AbstractMatrix{<:AbstractPanel})
```

### Reference Parameters and Frames
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
influence_coefficients(surface::AbstractMatrix, symmetric; xhat=SVector(1, 0, 0))
influence_coefficients(surfaces::AbstractVector{<:AbstractMatrix},
    surface_id, symmetric; xhat=SVector(1, 0, 0))
influence_coefficients!(AIC, surface::AbstractMatrix, symmetric; xhat=SVector(1, 0, 0))
influence_coefficients!(AIC, surfaces::AbstractVector{<:AbstractMatrix},
    surface_id, symmetric; xhat=SVector(1, 0, 0))
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
near_field_forces(surface::AbstractMatrix, ref, fs, symmetric, Γ)
near_field_forces(surfaces::AbstractVector{<:AbstractMatrix}, surface_id, ref, fs, symmetric, Γ)
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
```

### Vortex Lattice Panels
```@docs
VortexLattice.top_left
VortexLattice.top_center
VortexLattice.top_right
VortexLattice.bottom_left
VortexLattice.bottom_center
VortexLattice.bottom_right
VortexLattice.controlpoint
VortexLattice.normal
VortexLattice.get_core_size
VortexLattice.translate(panel::AbstractPanel, r)
VortexLattice.reflect(panel::AbstractPanel)
VortexLattice.induced_velocity
VortexLattice.panel_induced_velocity
VortexLattice.panel_circulation
VortexLattice.influence_coefficients!(AIC, receiving::AbstractMatrix{<:AbstractPanel}, sending::AbstractMatrix{<:AbstractPanel}, same_id, symmetric, xhat)
VortexLattice.left_center
VortexLattice.right_center
VortexLattice.top_vector
VortexLattice.left_vector
VortexLattice.right_vector
VortexLattice.bottom_vector
VortexLattice.flipy
VortexLattice.not_on_symmetry_plane
VortexLattice.trailing_induced_velocity
VortexLattice.bound_induced_velocity
```

### Freestream
```@docs
VortexLattice.body_to_stability_alpha
VortexLattice.body_to_wind_derivatives
VortexLattice.stability_to_body_alpha
VortexLattice.stability_to_wind_beta
VortexLattice.wind_to_body_derivatives
VortexLattice.wind_to_stability_beta
VortexLattice.freestream_velocity
VortexLattice.freestream_velocity_derivatives
VortexLattice.external_velocity
VortexLattice.external_velocity_derivatives
```

### Circulation
```@docs
VortexLattice.normal_velocity_derivatives
VortexLattice.normal_velocity_derivatives!
VortexLattice.circulation_derivatives
```

### Near-Field
```@docs
VortexLattice.body_to_frame
VortexLattice.near_field_forces_derivatives(surface::AbstractMatrix, ref, fs, symmetric, Γ, dΓ)
VortexLattice.near_field_forces_derivatives(surfaces::AbstractVector{<:AbstractMatrix}, surface_id, ref, fs, symmetric, Γ, dΓ)
```

### Far-Field
```@docs
VortexLattice.TrefftzPanel
VortexLattice.normal(panel::VortexLattice.TrefftzPanel)
VortexLattice.trefftz_panels
VortexLattice.panel_induced_drag
VortexLattice.vortex_induced_drag
```

## Index

```@index
```
