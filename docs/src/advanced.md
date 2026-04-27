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

## Rotors
For convenience some functionality is provided for generating rotors that are compatible with the nonlinear vortex lattice analysis. **generate_rotor** provides grids, ratios, sections, and invert_normals that can be used to create the rotor in the system.

```
grids, ratios, sections, invert_normals = generate_rotor(rotor_file, data_path)
system = System(grids; ratios, sections, invert_normals)
```

### Rotor file structure

This functionality requires a specific file system. Here is the file system for a simple rotor.

```bash
.
└── rotor_data
    ├── airfoils
    │   ├── my_airfoil.csv
    │   └── my_airfoil.dat
    └── rotors
        ├── my_rotor.csv
        ├── airfoils.csv
        ├── blade.csv
        ├── chorddist.csv
        ├── heightdist.csv
        ├── pitchdist.csv
        ├── sweepdist.csv
        └── centers.csv
```

data_path points to the rotor_data folder, rotor_file is the name of the rotor file inside rotor_data/rotors, in this case filename is "my_rotor.csv". Here is a list of the contents of each file with examples:

my_rotor.csv contains the radius of the hub and blade tip, the number of blades on the rotor, and the file that contains the rotor blade information.
``` 
property,file,description
Rtip,0.75, (m) Radius of blade tip
Rhub,0.0375, (m) Radius of hub
B,3, Number of blades
blade,blade.csv, Blade file
```

blade.csv contains files names for the chord, pitch, sweep, height, and airfoils for the blade and an optional file that defines the reference point on the airfoil for the measurements of pitch, sweep, and height. If centers.csv is not provided it defaults to the leading edge of the airfoil.
```
property,file,description
chorddist,chorddist.csv, Chord distribution
pitchdist,pitchdist.csv, Pitch distribution
sweepdist,sweepdist.csv, LE sweep distribution
heightdist,heightdist.csv, LE height distribution
airfoil_files,airfoils.csv, Airfoil distribution
airfoil_reference,centers.csv,Airfoil references
```

chorddist.csv, pitchdist.csv, sweepdist.csv, and heightdist.csv, each define the chord, pitch, sweep, and height of the airfoil normalized by the tip radius of the blade. Each file follows the same setup. Here is an example of the chorddist.csv file.

```
r/R,c/R
0.0,0.134
0.086,0.137106
0.16,0.144606
...
1.0,0.0375
```

As a note, all files (except airfoil_files.csv) should use the same r/R column (where R is the tip radius) while the second column will vary. The r/R column must always go from 0.0 to 1.0. The headings for the each column are arbitrary but headings must be provided.

airfoils.csv defines which airfoils are used along the length of the blade. The first column r/R does not need to match the other files. The contour file is the normalized contour of the airfoil, this is cuurently not used by VortexLattice.jl but the column is required but the files can be all zeros. The third column points to the files that define the polars of the airfoil.

```
r/R,Contour file,Aero file
0.0,naca4412.csv,naca4412.dat
0.368421052631579,naca4412.csv,naca4412.dat
0.6842105263157894,naca4412.csv,naca4412.dat
0.8947368421052632,naca4412.csv,naca4412.dat
1.0,naca4412.csv,naca4412.dat
```

Here is an example contour file:
```
x/c,y/c
1.0,0.0012489471548600977
0.9983816051162109,0.0017436209477767937
0.996268484921455,0.0023275052942059735
...
1.0,-0.0012489471548600977
```

Here is an example airfoil file. This file must follow the convention used in CCBlade.jl, though the Reynolds and Mach number are not used during calculation. The first line is an arbitrary information line, second line is Reynold number, the third is Mach number. The fourth line to the end contains the angle of attack, the coefficent of lift, and the coefficient of drag in that order.
```
Polars info
0.0
0.0
-180.00000   0.00000   0.60000
-175.00000   0.00000   0.60000
-170.00000   0.00000   0.60000
```

The airfoil_reference.csv file is an optional file with the intent of allowing the user to define height, sweep, and pitch at any point on the airfoil that is not the leading edge (as is common with wind turbines). The file format matches the height, sweep, and pitch files and so will not be shown here.

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