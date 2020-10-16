# VortexFlow

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://flow.byu.edu/VortexFlow.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/VortexFlow.jl/dev)
![](https://github.com/byuflowlab/VortexFlow.jl/workflows/Run%20tests/badge.svg)

*An Extensively Verified Vortex Lattice Method Flow Solver Written in Julia*

Authors: Taylor McDonnell and Andrew Ning

**VortexFlow** is a vortex lattice method package written in pure Julia.  It is designed to be fast, accurate (within theoretical limitations), and easy to use.  It has been extensively verified against results generated using Mark Drela's AVL, but has also been enhanced to incorporate features and modeling capabilties not present in AVL.

## Package Features
 - Custom Vortex Lattice Panels
  - Horseshoe Vortices (see Flight Vehicle Aerodynamics by Mark Drela[[1]](#1))
  - Vortex Rings (see Low-Speed Aerodynamics by Katz and Plotkin[[2]](#2))
 - Convenient Geometry Generation
  - From Pre-Existing Grid
  - From Lifting Surface Parameters
  - Symmmetric Geometries
 - Multiple Discretization Schemes
  - Uniform
  - Sine
  - Cosine
 - General Freestream Description
  - Freestream flow angles
  - Aircraft rotation components
  - Custom additional velocity field
 - Multiple Analyses
  - Near Field Forces in Body, Stability, or Wind Axes
  - Far Field Drag
  - Body and Stability Derivatives
 - Geometry Visualization using [WriteVTK](https://github.com/jipolanco/WriteVTK.jl)
 - Extensively verified against computational results generated using AVL.

## Installation

Enter the package manager by typing `]` and then run the following:

```julia
pkg> add https://github.com/byuflowlab/VortexFlow.jl
```

## Performance

This code has been optimized to be highly performant, primarily by maintaining type stability and minimizing allocations.  It should easily outperform other vortex lattice method codes written in other higher level languages.

## Usage

See the [documentation](https://flow.byu.edu/VortexFlow.jl/dev)

## References
<a id="1">[1]</a>
Drela, M. Flight Vehicle Aerodynamics. MIT Press, 2014.

<a id="2">[2]</a>
Katz, J., and Plotkin A. Low-Speed Aerodynamics. Vol. 13. Cambridge University Press, 2001.
