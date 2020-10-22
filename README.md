# VortexLattice

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://flow.byu.edu/VortexLattice.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://flow.byu.edu/VortexLattice.jl/dev)
![](https://github.com/byuflowlab/VortexLattice.jl/workflows/Run%20tests/badge.svg)

*A Comprehensive Julia implementation of the Vortex Lattice Method*

Authors: Taylor McDonnell and Andrew Ning

**VortexLattice** is a comprehensive pure-Julia implementation of the vortex lattice method.  It is designed to be fast, accurate (within theoretical limitations), easy to use, and applicable to arbitrary geometries and velocity fields.  It has been extensively verified against results generated using Mark Drela's AVL, but has also been enhanced to incorporate features and modeling capabilties not present in AVL.

![](docs/src/showoff.png)

## Package Features
 - Custom vortex lattice panels
   - Horseshoe vortices (see Flight Vehicle Aerodynamics by Mark Drela[[1]](#1))
   - Vortex rings (see Low-Speed Aerodynamics by Katz and Plotkin[[2]](#2))
   - Optional finite-core model
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
pkg> add https://github.com/byuflowlab/VortexLattice.jl
```

## Performance

This code has been optimized to be highly performant, primarily by maintaining type stability and minimizing allocations.  It should easily outperform other vortex lattice method codes written in other higher level languages.

## Usage

See the [documentation](https://flow.byu.edu/VortexLattice.jl/dev)

## References
<a id="1">[1]</a>
Drela, M. Flight Vehicle Aerodynamics. MIT Press, 2014.

<a id="2">[2]</a>
Katz, J., and Plotkin A. Low-Speed Aerodynamics. Cambridge University Press, 2001.
