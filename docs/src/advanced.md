# Advanced Use Cases

## Nonlinear Vortex Lattice
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
nonlinear_analysis!(system; max_iter=1, tol=1E-6, damping=0.01, print_iters=false)
```

**nonlinear_analysis!** populates the sections with angles of attack, coefficients of lift and drag, and force-per-unit length normalized by the density (defaulted to 1) and the reference velocity.

### Warnings
Using a significant number of spanwise panels (more than 75) with non-uniform spacing may result in non-physical behavior at the ends of the surface. In this case it is recommended to use uniform spacing.

Increasing **max_iter** will result in a more accurate solution but may also introduce noise into the solution depending on the surface.

Default **tol** is designed to not be achieved so that **max_iter** determines the length of the nonlinear analysis. If noise is present in the solution the solution may diverge and convergence may be impossible without reducing the **damping**.

Increasing **damping** causes larger step sizes at each iteration. Larger step sizes may accelerate convergence but may also lead to noise solutions. 

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