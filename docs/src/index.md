# VortexFlow

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://flow.byu.edu/VortexFlow.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/VortexFlow.jl/dev)
![](https://github.com/byuflowlab/VortexFlow.jl/workflows/Run%20tests/badge.svg)

*An Extensively Verified Vortex Lattice Method Flow Solver Written in Julia*

Authors: Taylor McDonnell and Andrew Ning

**VortexFlow** is a vortex lattice method package written in pure Julia.  It is designed to be fast, accurate (within theoretical limitations), and easy to use.  It has been extensively verified against results generated using Mark Drela's AVL, but has also been enhanced to incorporate features and modeling capabilties not present in AVL.

## Package Features
 - Custom vortex lattice panels
  - Horseshoe vortices (see Flight Vehicle Aerodynamics by Mark Drela[[1]](@ref References))
  - Vortex rings (see Low-Speed Aerodynamics by Katz and Plotkin[[2]](@ref References))
 - Convenient geometry generation
  - From pre-existing grid
  - From lifting surface parameters
  - Symmetric geometries
 - Multiple discretization schemes
  - Uniform
  - Sine
  - Cosine
 - General freestream description
  - Freestream flow angles
  - Aircraft rotation components
  - Additional velocity
 - Multiple analyses
  - Near field forces in body, stability, or wind Axes
  - Far field drag
  - Body and stability derivatives
 - Geometry visualization using [WriteVTK](https://github.com/jipolanco/WriteVTK.jl)
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
[1] Drela, M. Flight Vehicle Aerodynamics. MIT Press, 2014.

[2] Katz, J., and Plotkin A. Low-Speed Aerodynamics. Cambridge University Press, 2001.
