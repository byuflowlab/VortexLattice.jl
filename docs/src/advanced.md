# Advanced Use Cases

## Nonlinear Vortex Lattice (Still under development)
To incorporate airfoil polars a nonlinear vortex lattice method is used. In order to use this method system.sections must be correctly populated.

```
system = System(grids; sections)
```

sections is a vector of vector of SectionProperties objects. The length of grids and sections must be the same, for each grid there is a vector of SectionProperties. These vectors of section properties are created as follows:

```
section = grid_to_sections(grid, airfoils)
```

If grid is a 3 x n x m array, then airfoils is a vector of length m-1 with each object containing a CCBlade AlphaAF object that contains the polars for the airfoil. For more information see CCBlade.jl.

```
system = System(grids; sections)
nonlinear_analysis!(system)
```

**nonlinear_analysis!** populates the sections with angles of attack, and coefficient of lift and drag.

## Unsteady Simulation with PanelParticleWake

`PanelParticleWake` combines a small rolling buffer of panel rows (for the near wake) with a vortex particle field (for the far wake). Use it via the `simulate!` dispatcher:

```julia
Uinf(t) = SVector{3,Float64}(10.0, 0.0, 0.0)
Ωinf(t) = SVector{3,Float64}(0.0, 0.0, 0.0)
maneuver!(frames, system, wake, t) = nothing
t_range = 0.0:0.05:1.0

monitor = ForcesMonitor(length(t_range))

wake = simulate!(system, frames, maneuver!, Uinf, t_range, Ωinf;
    wake_type = PanelParticleWake,
    nwakerows = 2,
    max_particles = 10_000,
    method_trailing = OverlapPPS(1.3, 2),
    method_unsteady = OverlapPPS(1.3, 2),
    eta = 0.3,
    monitors = (monitor,),
    name = "my_sim",
    path = "output/",
)
```

Key parameters:
- `nwakerows`: number of panel rows to keep in the near-wake buffer before converting to particles.
- `max_particles`: maximum particle field size.
- `method_trailing`/`method_unsteady`: controls how panel rows are converted to particles (`NoShed`, `SigmaPPS`, or `OverlapPPS`).
- `eta`: relaxation parameter for the particle time integration (0–1).

The `System` must be allocated with `nw=fill(nwakerows, length(grids))` so the panel buffer has the correct row count.

## Restart / Checkpoints

By default, `simulate!` writes a restart checkpoint at every step when `path` is provided. To restart from a previous run:

```julia
# First run — writes checkpoints to output/
wake = simulate!(system, frames, maneuver!, Uinf, t_range, Ωinf;
    wake_type = PanelParticleWake,
    name = "my_sim", path = "output/", write_restart = true, ...)

# Partial re-run — start fresh, stop at step 5
simulate!(system, frames, maneuver!, Uinf, t_range[1:6], Ωinf;
    wake_type = PanelParticleWake,
    name = "my_sim", path = "output/", ...)

# Resume from checkpoint 5
wake = simulate!(system, frames, maneuver!, Uinf, t_range, Ωinf;
    wake_type = PanelParticleWake,
    name = "my_sim", path = "output/",
    restart_from = "output/my_sim",
    restart_idx = 5, ...)
```

You can also restore a checkpoint manually for post-processing:
```julia
restore_restart!(system, wake, frames, "output/my_sim"; idx=5)
```

## Fluid Domain

`FluidDomainMonitor` evaluates total velocity and vorticity on a rectilinear grid. Pass it in `monitors` to sample at every step, or use `evaluate_fluid_domain_from_restarts!` to post-process saved checkpoints without re-running the simulation.

```julia
fd = FluidDomainMonitor(
    range(-10.0, 30.0, step=5.0),   # x
    range(-10.0, 10.0, step=5.0),   # y
    range(-10.0, 10.0, step=5.0);   # z
    vtk_interval = 1,
    name = "fluid_domain",
    path = "output/fluid_domain",
)

# Option A: during simulation
wake = simulate!(system, frames, maneuver!, Uinf, t_range, Ωinf;
    wake_type = PanelParticleWake, monitors = (fd,), ...)

# Option B: post-processing from saved checkpoints
evaluate_fluid_domain_from_restarts!(fd, system, wake, frames, "output/my_sim")
```

After each evaluation `fd.velocity` and `fd.vorticity` hold the field values as `Array{SVector{3}, 3}` indexed `[ix, iy, iz]`.