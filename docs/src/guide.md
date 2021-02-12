# Guide

This guide demonstrates the basic steady analysis capabilities of VortexLattice. See the examples for more advanced uses of VortexLattice, including unsteady simulations.

We start by loading the package.

```@example guide
using VortexLattice
nothing #hide
```

Then we need to create our geometry.  While VortexLattice can handle multiple lifting surfaces, for this guide we will be analyzing a planar wing with the following geometric properties.

```@example guide
xle = [0.0, 0.4] # leading edge x-position
yle = [0.0, 7.5] # leading edge y-position
zle = [0.0, 0.0] # leading edge z-position
chord = [2.2, 1.8] # chord length
theta = [2.0*pi/180, 2.0*pi/180] # twist (in radians)
phi = [0.0, 0.0] # section rotation about the x-axis
nothing #hide
```

Note that we are only defining half the wing since the wing is symmetric about the X-Z plane.

We also need to define the number of panels and the discretization scheme in the spanwise and chordwise directions.  There are currently three discretization
scheme options: `Uniform()`, `Sine()`, and `Cosine()`.  To maximize the accuracy of our analysis we would like to use cosine spacing in the spanwise direction.  To do this, we need to use sine spacing on the right half of the wing (since once reflected across the y-z plane, sine spacing become cosine spacing).  

```@example guide
ns = 12 # number of spanwise panels
nc = 6  # number of chordwise panels
spacing_s = Sine() # spanwise discretization scheme
spacing_c = Uniform() # chordwise discretization scheme
nothing #hide
```

We generate our lifting surface using `wing_to_surface_panels`.  We use the keyword argument `mirror` to mirror our geometry across the X-Y plane.  A grid with the
panel corners and a matrix of vortex lattice panels representing the surface
of the wing is returned from this function.  Only the latter is needed for the analysis. The former is provided simply for the user's convenience.

```@example guide
grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)
nothing #hide
```

We could have also generated our lifting surface from a pre-existing grid using
`grid_to_surface_panels`.

The last step in defining our geometry is to combine all surfaces in a single vector.  Since we only have one surface, we create a vector with a single element.

```@example guide
surfaces = [surface]
nothing #hide
```

Now that we have generated our geometry we need to define our reference parameters and freestream properties. We use the following reference parameters

```@example guide
Sref = 30.0 # reference area
cref = 2.0  # reference chord
bref = 15.0 # reference span
rref = [0.50, 0.0, 0.0] # reference location for rotations/moments (typically the c.g.)
Vinf = 1.0 # reference velocity (magnitude)
ref = Reference(Sref, cref, bref, rref, Vinf)
nothing #hide
```

We use the following freestream properties.
```@example guide
alpha = 1.0*pi/180 # angle of attack
beta = 0.0 # sideslip angle
Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
fs = Freestream(Vinf, alpha, beta, Omega)
nothing #hide
```

Since the flow conditions are symmetric, we could have modeled one half of our wing and used symmetry to model the other half.  This, however, would give incorrect results for the lateral stability derivatives so we have instead mirrored our geometry across the X-Z plane.

```@example guide
symmetric = false
nothing #hide
```

We are now ready to perform a steady state analysis. We do so by calling the `steady_analysis` function. This function:
 - Finds the circulation distribution for a given set of panels and flow conditions
 - Performs a near-field analysis to find the forces on each panel, unless
   otherwise specified through the keyword argument `near_field_analysis`
 - Determines the derivatives of the near-field analysis forces with respect to
   the freestream variables, unless otherwise specified through the keyword
   argument `derivatives`

```@example guide
system = steady_analysis(surfaces, ref, fs; symmetric)
nothing #hide
```

The result of our analysis is an object of type `system` which holds the system state.  Note that the keyword argument `symmetric` is not strictly necessary, since by default it is set to false for each surface.

Once we have performed our steady state analysis (and associated near field analysis) we can extract the body force/moment coefficients using the function `body_forces`. These forces are returned in the reference frame specified by the keyword argument `frame`, which defaults to the body reference frame.

Note that a near field analysis must have been performed on `system` for this function to return sensible results (which is the default behavior when running an analysis).

```@example guide
CF, CM = body_forces(system; frame=Wind())

# extract aerodynamic forces
CD, CY, CL = CF
Cl, Cm, Cn = CM
nothing #hide
```

Numerical noise often corrupts drag estimates from near-field analyses, therefore, it is often more accurate to compute drag in the farfield on the Trefftz plane.

```@example guide
CDiff = far_field_drag(system)
nothing #hide
```

We can also extract the body and/or stability derivatives for the aircraft easily using the functions `body_derivatives` and/or `stability_derivatives`.  

Once again, note that the derivatives of the near-field analysis forces with respect to the freestream variables must have been previously calculated (which is the default behavior when running an analysis) for these functions to yield sensible results.

```@example guide
dCFb, dCMb = body_derivatives(system)

# traditional names for each body derivative
CXu, CYu, CZu = dCFb.u
CXv, CYv, CZv = dCFb.v
CXw, CYw, CZw = dCFb.w
CXp, CYp, CZp = dCFb.p
CXq, CYq, CZq = dCFb.q
CXr, CYr, CZr = dCFb.r
Clu, Cmu, Cnu = dCMb.u
Clv, Cmv, Cnv = dCMb.v
Clw, Cmw, Cnw = dCMb.w
Clp, Cmp, Cnp = dCMb.p
Clq, Cmq, Cnq = dCMb.q
Clr, Cmr, Cnr = dCMb.r

nothing #hide
```

```@example guide
dCFs, dCMs = stability_derivatives(system)

# traditional names for each stability derivative
CDa, CYa, CLa = dCFs.alpha
Cla, Cma, Cna = dCMs.alpha
CDb, CYb, CLb = dCFs.beta
Clb, Cmb, Cnb = dCMs.beta
CDp, CYp, CLp = dCFs.p
Clp, Cmp, Cnp = dCMs.p
CDq, CYq, CLq = dCFs.q
Clq, Cmq, Cnq = dCMs.q
CDr, CYr, CLr = dCFs.r
Clr, Cmr, Cnr = dCMs.r

nothing #hide
```

Visualizing the geometry (and results) may be done in Paraview after writing the associated visualization files using `write_vtk`.

```julia
properties = get_surface_properties(system)

write_vtk("simplewing", surfaces, properties)
```

![](simple-guide.png)

For visualization purposes, positive circulation is defined in the +i and +j directions.
