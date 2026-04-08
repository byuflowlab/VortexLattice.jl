# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

VortexLattice.jl is a pure-Julia implementation of the vortex lattice method (VLM) for steady and unsteady aerodynamic analysis of arbitrary lifting surfaces. Steady results are verified against AVL; unsteady against Katz & Plotkin.

The current `frames` branch extends the package with reference-frame handling, nonlinear/viscous corrections, rotor support (CCBlade), FMM acceleration (FastMultipole), and a particle wake option (FLOWVPM). Note `Project.toml` has uncommitted dependency changes adding these.

## Commands

Run from the package root:

```bash
# Run full test suite
julia --project=. -e 'using Pkg; Pkg.test()'

# Run a single test file directly (faster iteration; uses test/Project.toml)
julia --project=test test/test_frames.jl
julia --project=test test/rotor_hover.jl

# REPL with the package loaded
julia --project=. -e 'using Revise, VortexLattice'
```

`test/runtests.jl` is organized as `@testset` blocks (mostly AVL verification cases). To run just one set, copy its block into a scratch file or wrap with `if false ... end`.

## Architecture

The package is a flat set of files included from `src/VortexLattice.jl` in a specific order — that file is the index of what exists and what is exported. Key layering:

1. **Geometry & panels** (`panel.jl`, `geometry.jl`, `wake.jl`, `vspgeom.jl`, `rotors.jl`)
   `SurfacePanel` / `WakePanel` / `TrefftzPanel` are the core data types. Geometry is built either from a raw grid (`grid_to_surface_panels`) or from wing parameters with a discretization scheme (`Uniform`, `Sine`, `Cosine`). VSP `degengeom` import and rotor generation (via CCBlade) live alongside.

2. **System assembly** (`system.jl`, `reference.jl`, `freestream.jl`, `frames.jl`)
   `System` is the central mutable container holding panels, wakes, AIC matrix, RHS, circulation, derivatives, and per-panel `PanelProperties`. `Reference` defines the non-dimensionalization. `Freestream` defines flow conditions. `frames.jl` (new on this branch) introduces `ReferenceFrame` with conventions like `BackRightUp` / `ForwardRightDown` and `propagate_kinematics!` / `change_convention!` for converting between body/stability/wind frames.

3. **Solvers** (`induced.jl`, `circulation.jl`, `analyses.jl`, `unsteady.jl`, `nonlinear.jl`, `viscous.jl`, `fmm.jl`)
   `induced.jl` computes vortex-induced velocities (Biot–Savart with optional finite core); `circulation.jl` builds the AIC system; `analyses.jl` exposes the user-facing `steady_analysis[!]` and `unsteady_analysis[!]`. `nonlinear.jl` provides `SectionProperties` + `nonlinear_analysis!` for viscous/stall corrections via 2D polars; `viscous.jl` and `fmm.jl` extend this with viscous coupling and FastMultipole acceleration. `unsteady.jl` adds the higher-level `simulate!` driver (used with FLOWVPM particle wakes).

4. **Post-processing** (`nearfield.jl`, `farfield.jl`, `stability.jl`, `visualization.jl`)
   Near-field forces from panel circulations, Trefftz-plane far-field drag, body/stability derivatives, and VTK output via WriteVTK (`write_vtk`).

### Conventions to be aware of

- The `!`-suffixed analysis functions mutate an existing `System` and are the fast path; the non-mutating versions allocate one. Reuse `System` across calls when sweeping parameters.
- `const RHO = 1.0` in `VortexLattice.jl` — densities are folded out; forces returned are non-dimensional unless a `Reference` is supplied.
- Performance relies on type stability and `StaticArrays`; avoid introducing `Any`-typed fields or abstractly-typed containers in hot paths.
- Many root-level files (`*.csv`, `*.vtp`, `*.vts`, `courier/`, `vortex_lattice_simulation/`, `VortexLattice_rotor_data/`) are working artifacts from rotor/unsteady runs, not part of the package.
