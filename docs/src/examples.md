# Examples

These examples show how to use VortexLattice for various geometries, flow conditions, and analyses.  Many of these examples also provide a verification for the implementation of the vortex lattice method in this package.

```@contents
Pages = ["examples.md"]
Depth = 3
```

```@setup examples
# this is placed here to pre-install matplotlib so the documentation doesn't get cluttered with the installation print statements.
using Plots
pyplot()
```

## Steady State Analysis of a Wing

This example shows how to calculate aerodynamic coefficients and stability derivatives for a symmetric planar wing.

```@example planar-wing
using VortexLattice

# geometry (right half of the wing)
xle = [0.0, 0.4]
yle = [0.0, 7.5]
zle = [0.0, 0.0]
chord = [2.2, 1.8]
theta = [2.0*pi/180, 2.0*pi/180]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section

# discretization parameters
ns = 12
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()

# reference parameters
Sref = 30.0
cref = 2.0
bref = 15.0
rref = [0.50, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref, cref, bref, rref, Vinf)

# freestream parameters
alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha, beta, Omega)

# construct surface
grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

# create vector containing all surfaces
surfaces = [surface]

# we can use symmetry since the geometry and flow conditions are symmetric about the X-Z axis
symmetric = true

# perform steady state analysis
system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

# retrieve near-field forces
CF, CM = body_forces(system; frame=Wind())

# perform far-field analysis
CDiff = far_field_drag(system)

CD, CY, CL = CF
Cl, Cm, Cn = CM

nothing # hide
```

The aerodynamic coefficients predicted by VortexLattice are nearly identical to those predicted by AVL.

```@example planar-wing
using PrettyTables # hide
using Markdown # hide

CD_avl = 0.00247 # hide
CL_avl = 0.24454 # hide
Cm_avl = -0.02091 # hide
CDiff_avl = 0.00248 # hide

table = [  # hide
"``C_L``" CL CL_avl CL-CL_avl;  # hide
"``C_{Di}`` (nearfield)" CD CD_avl CD-CD_avl;  # hide
"``C_{Di}`` (farfield)" CDiff CDiff_avl CDiff-CDiff_avl;  # hide
"``C_M``" Cm Cm_avl Cm-Cm_avl  # hide
]  # hide
header = ["Coefficient", "VortexLattice", "AVL", "Difference"]  # hide

str = pretty_table(String, table, header; # hide
    backend=:text, # hide
    tf = tf_markdown, # hide
    alignment=[:l, :r, :r, :r], # hide
    formatters = (ft_printf("%0.5f", [2,3]), ft_printf("%0.1e", 4))) # hide

Markdown.parse(str) # hide
```

We can also generate files to visualize the results in Paraview using the function `write_vtk`.

```julia
properties = get_surface_properties(system)

write_vtk("symmetric-planar-wing", surfaces, properties; symmetric)
```

![](symmetric-planar-wing.png)

For asymmetric flow conditions and/or to obtain accurate asymmetric stability derivatives we can use the keyword argument `mirror` when constructing the geometry to reflect the geometry across the X-Z plane prior to the analysis.  We also set the `symmetric` flag to `false` since we are no longer using symmetry in the analysis.

```@example planar-wing

# construct geometry with mirror image
grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    fc=fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)

# symmetry is not used in the analysis
symmetric = false

# create vector containing all surfaces
surfaces = [surface]

# perform steady state analysis
system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

# retrieve near-field forces
CF, CM = body_forces(system; frame=Wind())

# perform far-field analysis
CDiff = far_field_drag(system)

CD, CY, CL = CF
Cl, Cm, Cn = CM

nothing # hide
```

Once again, the aerodynamic coefficients predicted by VortexLattice are nearly identical to those predicted by AVL.

```@example planar-wing
using PrettyTables # hide
using Markdown # hide

CD_avl = 0.00247 # hide
CL_avl = 0.24454 # hide
Cm_avl = -0.02091 # hide
CDiff_avl = 0.00248 # hide

table = [ # hide
"``C_L``" CL CL_avl CL-CL_avl; # hide
"``C_{Di}`` (nearfield)" CD CD_avl CD-CD_avl; # hide
"``C_{Di}`` (farfield)" CDiff CDiff_avl CDiff-CDiff_avl; # hide
"``C_M``" Cm Cm_avl Cm-Cm_avl # hide
] # hide
header = ["Coefficient", "VortexLattice", "AVL", "Difference"] # hide

str = pretty_table(String, table, header; # hide
    backend=:text, # hide
    tf = tf_markdown, # hide
    alignment=[:l, :r, :r, :r], # hide
    formatters = (ft_printf("%0.5f", [2,3]), ft_printf("%0.1e", 4))) # hide

Markdown.parse(str) # hide
```

The stability derivatives are also very close to those predicted by AVL.

```@example planar-wing

dCF, dCM = stability_derivatives(system)

CDa, CYa, CLa = dCF.alpha
Cla, Cma, Cna = dCM.alpha
CDb, CYb, CLb = dCF.beta
Clb, Cmb, Cnb = dCM.beta
CDp, CYp, CLp = dCF.p
Clp, Cmp, Cnp = dCM.p
CDq, CYq, CLq = dCF.q
Clq, Cmq, Cnq = dCM.q
CDr, CYr, CLr = dCF.r
Clr, Cmr, Cnr = dCM.r

nothing # hide
```

```@example planar-wing
using PrettyTables # hide
using Markdown # hide

CLa_avl =   4.663214 # hide
CLb_avl =   0.0 # hide
CYa_avl =   0.0 # hide
CYb_avl =  -0.000002 # hide
Cla_avl =   0.0 # hide
Clb_avl =  -0.025435 # hide
Cma_avl =  -0.397758 # hide
Cmb_avl =   0.0 # hide
Cna_avl =   0.0 # hide
Cnb_avl =   0.000452 # hide
CLp_avl =   0.0 # hide
CLq_avl =   5.649411 # hide
CLr_avl =   0.0 # hide
CYp_avl =   0.049063 # hide
CYq_avl =   0.0 # hide
CYr_avl =  -0.000828 # hide
Clp_avl =  -0.524750 # hide
Clq_avl =   0.0 # hide
Clr_avl =   0.064456 # hide
Cmp_avl =   0.0 # hide
Cmq_avl =  -1.270212 # hide
Cmr_avl =   0.0 # hide
Cnp_avl =  -0.019175 # hide
Cnq_avl =   0.0 # hide
Cnr_avl =  -0.000931 # hide

table = [ # hide
"``C_{La}``" CLa CLa_avl CLa-CLa_avl; # hide
"``C_{Lb}``" CLb CLb_avl CLb-CLb_avl; # hide
"``C_{Ya}``" CYa CYa_avl CYa-CYa_avl; # hide
"``C_{Yb}``" CYb CYb_avl CYb-CYb_avl; # hide
"``C_{la}``" Cla Cla_avl Cla-Cla_avl; # hide
"``C_{lb}``" Clb Clb_avl Clb-Clb_avl; # hide
"``C_{ma}``" Cma Cma_avl Cma-Cma_avl; # hide
"``C_{mb}``" Cmb Cmb_avl Cmb-Cmb_avl; # hide
"``C_{na}``" Cna Cna_avl Cna-Cna_avl; # hide
"``C_{nb}``" Cnb Cnb_avl Cnb-Cnb_avl; # hide
"``C_{Lp}``" CLp CLp_avl CLp-CLp_avl; # hide
"``C_{Lq}``" CLq CLq_avl CLq-CLq_avl; # hide
"``C_{Lr}``" CLr CLr_avl CLr-CLr_avl; # hide
"``C_{Yp}``" CYp CYp_avl CYp-CYp_avl; # hide
"``C_{Yq}``" CYq CYq_avl CYq-CYq_avl; # hide
"``C_{Yr}``" CYr CYr_avl CYr-CYr_avl; # hide
"``C_{lp}``" Clp Clp_avl Clp-Clp_avl; # hide
"``C_{lq}``" Clq Clq_avl Clq-Clq_avl; # hide
"``C_{lr}``" Clr Clr_avl Clr-Clr_avl; # hide
"``C_{mp}``" Cmp Cmp_avl Cmp-Cmp_avl; # hide
"``C_{mq}``" Cmq Cmq_avl Cmq-Cmq_avl; # hide
"``C_{mr}``" Cmr Cmr_avl Cmr-Cmr_avl; # hide
"``C_{np}``" Cnp Cnp_avl Cnp-Cnp_avl; # hide
"``C_{nq}``" Cnq Cnq_avl Cnq-Cnq_avl; # hide
"``C_{nr}``" Cnr Cnr_avl Cnr-Cnr_avl; # hide
] # hide
header = ["Coefficient", "VortexLattice", "AVL", "Difference"] # hide

str = pretty_table(String, table, header; # hide
    backend=:text, # hide
    tf = tf_markdown, # hide
    alignment=[:l, :r, :r, :r], # hide
    formatters = (ft_printf("%0.5f", [2,3]), ft_printf("%0.1e", 4))) # hide

Markdown.parse(str) # hide
```

Visualizing the geometry now shows the circulation distribution across the entire wing.

```julia
properties = get_surface_properties(system)

write_vtk("mirrored-planar-wing", surfaces, properties; symmetric)
```

![](mirrored-planar-wing.png)

## Steady State Analysis of a Wing with Dihedral

This example shows how to calculate aerodynamic coefficients and stability derivatives for a simple wing with dihedral.

```@example wing-with-dihedral
using VortexLattice

xle = [0.0, 0.4]
yle = [0.0, 7.5]
zle = [0.0, 3.0]
chord = [2.2, 1.8]
theta = [2.0*pi/180, 2.0*pi/180]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) #camberline function for each section

ns = 12
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()
mirror = false
symmetric = true

Sref = 30.0
cref = 2.0
bref = 15.0
rref = [0.50, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref, cref, bref, rref, Vinf)

alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha, beta, Omega)

# declare symmetry
symmetric = true

# construct surface
grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

# create vector containing all surfaces
surfaces = [surface]

# perform steady state analysis
system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

# retrieve near-field forces
CF, CM = body_forces(system; frame=Wind())

# perform far-field analysis
CDiff = far_field_drag(system)

CD, CY, CL = CF
Cl, Cm, Cn = CM

nothing # hide
```

The results predicted by VortexLattice are close to those predicted by AVL, with the difference primarily explained by the manner in which the normal vector is defined in VortexLattice and AVL, respectively.

```@example wing-with-dihedral
using PrettyTables # hide
using Markdown # hide

CD_avl = 0.00248 # hide
CL_avl = 0.24808 # hide
Cm_avl = -0.02250 # hide
CDiff_avl = 0.0024671 # hide

table = [ # hide
"``C_L``" CL CL_avl CL-CL_avl; # hide
"``C_{Di}`` (nearfield)" CD CD_avl CD-CD_avl; # hide
"``C_{Di}`` (farfield)" CDiff CDiff_avl CDiff-CDiff_avl; # hide
"``C_M``" Cm Cm_avl Cm-Cm_avl # hide
] # hide
header = ["Coefficient", "VortexLattice", "AVL", "Difference"] # hide

str = pretty_table(String, table, header; # hide
    backend=:text, # hide
    tf = tf_markdown, # hide
    alignment=[:l, :r, :r, :r], # hide
    formatters = (ft_printf("%0.5f", [2,3]), ft_printf("%0.1e", 4))) # hide

Markdown.parse(str) # hide
```

If we set the normal vectors in VortexLattice equal to those used in AVL, the results are even closer, though not necessarily more accurate.

```@example wing-with-dihedral

using LinearAlgebra

# function to construct a normal vector the way AVL does
#  - `ds` is a line representing the leading edge
#  - `theta` is the incidence angle, taken as a rotation (+ by RH rule) about
#        the surface's spanwise axis projected onto the Y-Z plane.
function avl_normal_vector(ds, theta)

    st, ct = sincos(theta)

    # bound vortex vector
    bhat = ds/norm(ds)

    # chordwise strip normal vector
    shat = [0, -ds[3], ds[2]]/sqrt(ds[2]^2+ds[3]^2)

    # camberline vector
    chat = [ct, -st*shat[2], -st*shat[3]]

    # normal vector perpindicular to camberline and bound vortex for entire chordwise strip
    ncp = cross(chat, ds)
    return ncp / norm(ncp) # normal vector used by AVL
end

# new normal vector
ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

# overwrite normal vector for each panel
for i = 1:length(surface)
    surface[i] = set_normal(surface[i], ncp)
end

# create vector containing all surfaces
surfaces = [surface]

# perform steady state analysis
system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

# retrieve near-field forces
CF, CM = body_forces(system; frame=Wind())

# perform far-field analysis
CDiff = far_field_drag(system)

CD, CY, CL = CF
Cl, Cm, Cn = CM

nothing # hide
```

```@example wing-with-dihedral
using PrettyTables # hide
using Markdown # hide

CD_avl = 0.00248 # hide
CL_avl = 0.24808 # hide
Cm_avl = -0.02250 # hide
CDiff_avl = 0.0024671 # hide

table = [ # hide
"``C_L``" CL CL_avl CL-CL_avl; # hide
"``C_{Di}`` (nearfield)" CD CD_avl CD-CD_avl; # hide
"``C_{Di}`` (farfield)" CDiff CDiff_avl CDiff-CDiff_avl; # hide
"``C_M``" Cm Cm_avl Cm-Cm_avl # hide
] # hide
header = ["Coefficient", "VortexLattice", "AVL", "Difference"] # hide

str = pretty_table(String, table, header; # hide
    backend=:text, # hide
    tf = tf_markdown, # hide
    alignment=[:l, :r, :r, :r], # hide
    formatters = (ft_printf("%0.5f", [2,3]), ft_printf("%0.1e", 4))) # hide

Markdown.parse(str) # hide
```

```julia
properties = get_surface_properties(system)

write_vtk("wing-with-dihedral", surfaces, properties; symmetric)
```

![](wing-with-dihedral.png)

## Steady State Analysis of a Wing and Tail

This example shows how to calculate aerodynamic coefficients and stability derivatives for multiple lifting surfaces.

```@example wing-tail
using VortexLattice

# wing
xle = [0.0, 0.2]
yle = [0.0, 5.0]
zle = [0.0, 1.0]
chord = [1.0, 0.6]
theta = [2.0*pi/180, 2.0*pi/180]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section
ns = 12
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()
mirror = false

# horizontal stabilizer
xle_h = [0.0, 0.14]
yle_h = [0.0, 1.25]
zle_h = [0.0, 0.0]
chord_h = [0.7, 0.42]
theta_h = [0.0, 0.0]
phi_h = [0.0, 0.0]
fc_h = fill((xc) -> 0, 2) #camberline function for each section
ns_h = 6
nc_h = 3
spacing_s_h = Uniform()
spacing_c_h = Uniform()
mirror_h = false

# vertical stabilizer
xle_v = [0.0, 0.14]
yle_v = [0.0, 0.0]
zle_v = [0.0, 1.0]
chord_v = [0.7, 0.42]
theta_v = [0.0, 0.0]
phi_v = [0.0, 0.0]
fc_v = fill((xc) -> 0, 2) #camberline function for each section
ns_v = 5
nc_v = 3
spacing_s_v = Uniform()
spacing_c_v = Uniform()
mirror_v = false

Sref = 9.0
cref = 0.9
bref = 10.0
rref = [0.5, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref, cref, bref, rref, Vinf)

alpha = 5.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha, beta, Omega)

symmetric = [true, true, false]

# generate surface panels for wing
wgrid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

# generate surface panels for horizontal tail
hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
    mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
translate!(hgrid, [4.0, 0.0, 0.0])
translate!(htail, [4.0, 0.0, 0.0])

# generate surface panels for vertical tail
vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
    mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
translate!(vgrid, [4.0, 0.0, 0.0])
translate!(vtail, [4.0, 0.0, 0.0])

grids = [wgrid, hgrid, vgrid]
surfaces = [wing, htail, vtail]
surface_id = [1, 2, 3]

system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)

CF, CM = body_forces(system; frame=Wind())

CDiff = far_field_drag(system)

CD, CY, CL = CF
Cl, Cm, Cn = CM

nothing # hide
```

The results predicted by VortexLattice are close to those predicted by AVL (with the finite core model disabled in AVL), with the difference primarily explained by the manner in which the normal vector is defined in VortexLattice and AVL, respectively.

```@example wing-tail
using PrettyTables # hide
using Markdown # hide

CD_avl = 0.01060 # hide
CL_avl = 0.60478 # hide
Cm_avl = -0.02700 # hide
CDiff_avl = 0.0104282 # hide

table = [ # hide
"``C_L``" CL CL_avl CL-CL_avl; # hide
"``C_{Di}`` (nearfield)" CD CD_avl CD-CD_avl; # hide
"``C_{Di}`` (farfield)" CDiff CDiff_avl CDiff-CDiff_avl; # hide
"``C_M``" Cm Cm_avl Cm-Cm_avl # hide
] # hide
header = ["Coefficient", "VortexLattice", "AVL", "Difference"] # hide

str = pretty_table(String, table, header; # hide
    backend=:text, # hide
    tf = tf_markdown, # hide
    alignment=[:l, :r, :r, :r], # hide
    formatters = (ft_printf("%0.5f", [2,3]), ft_printf("%0.1e", 4))) # hide

Markdown.parse(str) # hide
```

If we set the normal vectors in VortexLattice equal to those used in AVL, the results are closer, though not necessarily more accurate.

```@example wing-tail

using LinearAlgebra

# function to construct a normal vector the way AVL does
#  - `ds` is a line representing the leading edge
#  - `theta` is the incidence angle, taken as a rotation (+ by RH rule) about
#        the surface's spanwise axis projected onto the Y-Z plane.
function avl_normal_vector(ds, theta)

    st, ct = sincos(theta)

    # bound vortex vector
    bhat = ds/norm(ds)

    # chordwise strip normal vector
    shat = [0, -ds[3], ds[2]]/sqrt(ds[2]^2+ds[3]^2)

    # camberline vector
    chat = [ct, -st*shat[2], -st*shat[3]]

    # normal vector perpindicular to camberline and bound vortex for entire chordwise strip
    ncp = cross(chat, ds)
    return ncp / norm(ncp) # normal vector used by AVL
end

# new normal vector for the wing
ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

# overwrite normal vector for each wing panel
for i = 1:length(wing)
    wing[i] = set_normal(wing[i], ncp)
end
surfaces[1] = wing

# perform steady state analysis
system = steady_analysis(surfaces, ref, fs; symmetric=symmetric)

# retrieve near-field forces
CF, CM = body_forces(system; frame=Wind())

# perform far-field analysis
CDiff = far_field_drag(system)

CD, CY, CL = CF
Cl, Cm, Cn = CM

nothing # hide
```

```@example wing-tail
using PrettyTables # hide
using Markdown # hide

CD_avl = 0.01060 # hide
CL_avl = 0.60478 # hide
Cm_avl = -0.02700 # hide
CDiff_avl = 0.0104282 # hide

table = [ # hide
"``C_L``" CL CL_avl CL-CL_avl; # hide
"``C_{Di}`` (nearfield)" CD CD_avl CD-CD_avl; # hide
"``C_{Di}`` (farfield)" CDiff CDiff_avl CDiff-CDiff_avl; # hide
"``C_M``" Cm Cm_avl Cm-Cm_avl # hide
] # hide
header = ["Coefficient", "VortexLattice", "AVL", "Difference"] # hide

str = pretty_table(String, table, header; # hide
    backend=:text, # hide
    tf = tf_markdown, # hide
    alignment=[:l, :r, :r, :r], # hide
    formatters = (ft_printf("%0.5f", [2,3]), ft_printf("%0.1e", 4))) # hide

Markdown.parse(str) # hide
```

To achieve a theoretically identical setup as AVL we can place all our panels in the X-Y plane and then set the normal vector manually to match the actual lifting geometry.  In our case this involves removing the small amount of twist on the wing when creating the wing surface panels.

```@example wing-tail
using VortexLattice

# wing
xle = [0.0, 0.2]
yle = [0.0, 5.0]
zle = [0.0, 1.0]
chord = [1.0, 0.6]
theta = [0.0, 0.0]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section
ns = 12
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()
mirror = false

# horizontal stabilizer
xle_h = [0.0, 0.14]
yle_h = [0.0, 1.25]
zle_h = [0.0, 0.0]
chord_h = [0.7, 0.42]
theta_h = [0.0, 0.0]
phi_h = [0.0, 0.0]
fc_h = fill((xc) -> 0, 2) # camberline function for each section
ns_h = 6
nc_h = 3
spacing_s_h = Uniform()
spacing_c_h = Uniform()
mirror_h = false

# vertical stabilizer
xle_v = [0.0, 0.14]
yle_v = [0.0, 0.0]
zle_v = [0.0, 1.0]
chord_v = [0.7, 0.42]
theta_v = [0.0, 0.0]
phi_v = [0.0, 0.0]
fc_v = fill((xc) -> 0, 2) # camberline function for each section
ns_v = 5
nc_v = 3
spacing_s_v = Uniform()
spacing_c_v = Uniform()
mirror_v = false

Sref = 9.0
cref = 0.9
bref = 10.0
rref = [0.5, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref, cref, bref, rref, Vinf)

alpha = 5.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha, beta, Omega)

symmetric = [true, true, false]

# generate surface panels for wing
wgrid, wing = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

# generate surface panels for horizontal tail
hgrid, htail = wing_to_surface_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
    mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
translate!(hgrid, [4.0, 0.0, 0.0])
translate!(htail, [4.0, 0.0, 0.0])

# generate surface panels for vertical tail
vgrid, vtail = wing_to_surface_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
    mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
translate!(vgrid, [4.0, 0.0, 0.0])
translate!(vtail, [4.0, 0.0, 0.0])

# now set normal vectors manually
ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

# overwrite normal vector for each wing panel
for i = 1:length(wing)
    wing[i] = set_normal(wing[i], ncp)
end

grids = [wgrid, hgrid, vgrid]
surfaces = [wing, htail, vtail]
surface_id = [1, 2, 3]

system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)

CF, CM = body_forces(system; frame=Stability())

CDiff = far_field_drag(system)

CD, CY, CL = CF
Cl, Cm, Cn = CM

nothing # hide
```

The resulting aerodynamic coefficients now match very closely with AVL.

```@example wing-tail
using PrettyTables # hide
using Markdown # hide

CD_avl = 0.01060 # hide
CL_avl = 0.60478 # hide
Cm_avl = -0.02700 # hide
CDiff_avl = 0.0104282 # hide

table = [ # hide
"``C_L``" CL CL_avl CL-CL_avl; # hide
"``C_{Di}`` (nearfield)" CD CD_avl CD-CD_avl; # hide
"``C_{Di}`` (farfield)" CDiff CDiff_avl CDiff-CDiff_avl; # hide
"``C_M``" Cm Cm_avl Cm-Cm_avl # hide
] # hide
header = ["Coefficient", "VortexLattice", "AVL", "Difference"] # hide

str = pretty_table(String, table, header; # hide
    backend=:text, # hide
    tf = tf_markdown, # hide
    alignment=[:l, :r, :r, :r], # hide
    formatters = (ft_printf("%0.5f", [2,3]), ft_printf("%0.1e", 4))) # hide

Markdown.parse(str) # hide
```

By comparing these results with previous results we can see exactly how much restricting surface panels in the X-Y plane changes the results from the vortex lattice method.

```julia
properties = get_surface_properties(system)

write_vtk("wing-tail", surfaces, properties; symmetric)
```

![](wing-tail.png)

## Viscous Drag Correction

This example shows how to make corerections to account for Viscous Drag

```@example viscous-drag
#geometry of right half of wing
xle = [0.0, 0.4]
yle = [0.0, 7.5]
zle = [0.0, 0.0]
chords = [2.2, 1.8]
theta = [0.0*pi/180, 0.0*pi/180]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section

# discretization parameters
ns = 12
nc = 1
spacing_s = Uniform()
spacing_c = Uniform()

# reference parameters
Sref = 30.0
cref = 2.0
bref = 15.0
rref = [0.50, 0.0, 0.0]
rho = 1.225
mu = 1.81e-5
Vinf = 1
ref = Reference(Sref, cref, bref, rref, Vinf)

# freestream parameters
alpha = 0.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha, beta, Omega)

# construct surface
grid, surface = wing_to_surface_panels(xle, yle, zle, chords, theta, phi, ns, nc;
    fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

# combine all grid representations of surfaces into a single vector
grids = [grid]

# calculate lifting line geometry
r, c = lifting_line_geometry(grids)

```
We will now read in the airfoil data, provide the reynolds number, and perform a sweep over angle of attack.

```@example viscous-drag

#read in xfoil data
res = [150000   320000  490000  660000  830000  1000000]
cls = [
    -0.647901   -0.865237   -0.897876   -0.921655   -0.924562   -0.912016
    -0.688353   -0.839616   -0.868863   -0.886064   -0.865452   -0.848249
    -0.719029   -0.815164   -0.835761   -0.819934   -0.798952   -0.785489
    -0.724569   -0.780939   -0.790784   -0.75616    -0.734558   -0.722138
    -0.709499   -0.745029   -0.722006   -0.691479   -0.670796   -0.653588
    -0.684209   -0.706653   -0.656801   -0.622444   -0.599316   -0.582615
    -0.649806   -0.650326   -0.586158   -0.548678   -0.530815   -0.522789
    -0.611421   -0.575085   -0.509326   -0.483783   -0.466995   -0.460582
    -0.569993   -0.501801   -0.444908   -0.414056   -0.407301   -0.407514
    -0.526773   -0.4288     -0.370251   -0.353972   -0.349998   -0.354136
    -0.481819   -0.354335   -0.306358   -0.29439    -0.29748    -0.300658
    -0.437417   -0.283642   -0.239317   -0.239731   -0.243597   -0.24657
    -0.355987   -0.206903   -0.183994   -0.18603    -0.189385   -0.191869
    -0.272751   -0.142592   -0.129707   -0.132465   -0.135185   -0.137191
    -0.19158    -0.0804409  -0.0758313  -0.0786082  -0.0808187  -0.0825197
    -0.113186   -0.0229882  -0.0221598  -0.0245126  -0.0263089  -0.0275698
    -0.0424867   0.02954     0.0306951   0.0293735   0.0280789   0.026986
    0.0246814   0.0807037   0.0829588   0.0828598   0.0825106   0.082028
    0.106665    0.131261    0.135121    0.136025    0.13604     0.136118
    0.231398    0.18        0.186528    0.188861    0.18975     0.190195
    0.332249    0.236334    0.235567    0.24062     0.243155    0.244114
    0.379965    0.309605    0.28689     0.289299    0.293824    0.296698
    0.426923    0.395662    0.349986    0.341236    0.343296    0.346701
    0.47316     0.470653    0.427291    0.4041      0.396623    0.397817
    0.520035    0.517771    0.508826    0.477587    0.461357    0.454366
    0.567195    0.564904    0.566138    0.555519    0.53455     0.521906
    0.614409    0.612375    0.613351    0.615857    0.610176    0.594746
    0.661897    0.659498    0.660533    0.663514    0.665885    0.668431
    0.709736    0.706256    0.708341    0.710808    0.713209    0.715556
    0.75633     0.753501    0.755084    0.758106    0.760814    0.76251
    0.802545    0.800147    0.80302     0.805617    0.807281    0.80878
    0.848128    0.846923    0.849982    0.85223     0.853764    0.855129
    0.891559    0.892859    0.895978    0.898479    0.900326    0.90197
    0.933165    0.936628    0.940798    0.944134    0.946547    0.948488
    0.971578    0.978011    0.984136    0.988182    0.991548    0.994097
    1.00277     1.01615     1.02564     1.03238     1.03724     1.04163
    1.02813     1.05281     1.06639     1.07653     1.08365     1.08834
    1.056       1.08886     1.10721     1.12135     1.12849     1.13565
    1.08774     1.1238      1.14938     1.16289     1.1746      1.17986
    1.12409     1.15418     1.184       1.2062      1.21487     1.2263
    1.16485     1.18227     1.22046     1.24182     1.25905     1.27023
    1.20165     1.21108     1.25208     1.27826     1.29944     1.30683
    1.24163     1.23119     1.27154     1.31085     1.32602     1.34906
    1.26658     1.25802     1.29801     1.33667     1.36192     1.38469
    1.29076     1.28117     1.32457     1.35789     1.39391     1.41224
    1.29915     1.29855     1.343       1.38634     1.41698     1.44162
  ]
cds = [
    0.0654349  0.0296189   0.0228493   0.0190337   0.0179236   0.0159423
    0.0553959  0.0277749   0.0208117   0.018136    0.015916    0.0148696
    0.0459334  0.023951    0.0193385   0.0172538   0.0150951   0.0137197
    0.0385589  0.0226159   0.018273    0.0154623   0.0143091   0.0129814
    0.0338295  0.0210759   0.0172401   0.0144747   0.0130531   0.0123372
    0.0296124  0.0200311   0.0153633   0.0136195   0.0121329   0.0112241
    0.0263312  0.0176739   0.0143444   0.0125092   0.0114031   0.0105385
    0.0242647  0.0165613   0.0135285   0.0116141   0.0106384   0.00998252
    0.0225039  0.0151111   0.0123357   0.0108986   0.00996728  0.00937767
    0.0204981  0.0141967   0.0114408   0.0101346   0.00934501  0.00887346
    0.0190007  0.0131301   0.0106836   0.00952052  0.0087898   0.00836871
    0.0179213  0.0122845   0.0100597   0.00899888  0.00840774  0.00796738
    0.0168996  0.0114676   0.00949084  0.00858268  0.00803717  0.00766273
    0.0161244  0.0107389   0.00892733  0.00811376  0.00766265  0.00735471
    0.0152867  0.0100302   0.00841999  0.00771093  0.00730886  0.00704231
    0.0143646  0.00934181  0.00795968  0.00736196  0.00699525  0.00677437
    0.0135524  0.00865968  0.00750625  0.00701264  0.00672623  0.0065265
    0.0130552  0.0081212   0.00706139  0.00664781  0.00643271  0.00628345
    0.0129214  0.00770648  0.00674272  0.00634507  0.00613936  0.00601957
    0.0126324  0.00750968  0.00649519  0.00611761  0.00591863  0.00580701
    0.0121033  0.00757967  0.00633231  0.00595413  0.00578714  0.00567214
    0.0117236  0.00775288  0.00641919  0.00585671  0.00561447  0.00554661
    0.0115409  0.00794312  0.00665121  0.00599588  0.00563508  0.00546568
    0.01155    0.00811528  0.00692937  0.00625737  0.00583928  0.00558312
    0.0116789  0.00830769  0.00722311  0.00656245  0.00610511  0.00579274
    0.0119019  0.00855542  0.0074681   0.0068753   0.0064166   0.00607888
    0.0121988  0.00884388  0.00774843  0.00713832  0.00673198  0.00638581
    0.012564   0.00918564  0.00806183  0.00740902  0.00702941  0.00674493
    0.0129882  0.00956289  0.00837725  0.00775605  0.00736807  0.00708523
    0.0134068  0.00996367  0.00879841  0.00816252  0.00775411  0.00750858
    0.0138772  0.0104175   0.00923359  0.0086287   0.00829339  0.00805655
    0.0144083  0.0109608   0.00982798  0.00927505  0.00894606  0.00870128
    0.0149293  0.0116238   0.0105898   0.0100486   0.00969578  0.0094194
    0.0156143  0.0125228   0.0115174   0.0109321   0.0105371   0.0102365
    0.0165971  0.0136909   0.0126092   0.011985    0.0115137   0.0111661
    0.0182657  0.0151833   0.0138636   0.0130141   0.0124229   0.0119165
    0.0206356  0.0167678   0.0151394   0.0140169   0.0132418   0.0127244
    0.022962   0.0182899   0.016323    0.0149031   0.0141601   0.013441
    0.0253224  0.0197913   0.0172989   0.0159966   0.0148978   0.0143806
    0.0274142  0.021502    0.0187484   0.0168503   0.0160522   0.0150538
    0.0303662  0.0231384   0.0199176   0.0181825   0.0167942   0.0158637
    0.0327772  0.024691    0.0210983   0.0192633   0.0177298   0.0171334
    0.0366253  0.0273164   0.023116    0.0203368   0.0192977   0.0178622
    0.0395347  0.0291617   0.0247909   0.021866    0.0201695   0.0187629
    0.0425926  0.0312784   0.026501    0.023791    0.0213097   0.020106
    0.0479709  0.0343795   0.0289269   0.0253349   0.0230689   0.0213966
    ]

#generate cl to cd function (this is generating the template funciton provided)
cl_2_cd = VortexLattice.generate_cl_2_cd(cls,cds,res)

#select operating reynolds number and set freestream
re = 1000000
Vinf = re*mu/(c[1][1]*rho)
ref = Reference(Sref, cref, bref, rref, Vinf)

#calculate the width of each panel (this works only for uniform spacing)
w = (yle[2] - yle[1])/ns

#initialization for alpha sweep  
alphas = -10:.5:20                              #discrete angles of attack
viscous = zeros(n_alphas)                       #Cd_total for each value of alpha
inviscid = zeros(n_alphas)                      #Cd_induced for each value of alpha
local_viscous = zeros(3,ns,n_alphas)            #cfs_local_total for each spanwise location
local_inviscid =  zeros(3,ns,n_alphas)          #cfs_local_induced for each spanwise location
n_alphas = length(cls[:,1])                     #number of discrete angles of attack


#sweep of angle of attack
for j in 1:n_alphas   
    local_freestream = Freestream(Vinf, alphas[j]*pi/180, beta, Omega);
    system_local = steady_analysis!(system, surfaces, ref, local_freestream; symmetric);
    cfs, cms = lifting_line_coefficients(system_local, r, c; frame=Wind());
    local_inviscid[:,:,j] = cfs[1]
    inviscid[j] = body_forces(system_local; frame=Wind())[1][1]
    local_viscous[:,:,j] = VortexLattice.strip_theory_drag!(cfs, cms, system_local, r, c, cl_2_cd, re)[1]
    viscous[j] = sum((local_viscous[1,:,j]) .*c_local)*w*2/Sref
end

```

We can see that the results match those from Xflr5 identically.

![](strip-theory-verification.png)

## Sudden Acceleration of a Rectangular Wing into a Constant-Speed Forward Flight

This example shows how to predict the transient forces and moments on a rectangular wing when suddenly accelerated into forward flight at a five degree angle.

```@example rectangular-wing-sudden-acceleration
# Katz and Plotkin: Figures 13.34 and 13.35
# AR = [4, 8, 12, 20, ∞]
# Vinf*Δt/c = 1/16
# α = 5°

using VortexLattice

AR = [4, 8, 12, 20, 1e3] # last aspect ratio is essentially infinite

system = Vector{Any}(undef, length(AR))
surface_history = Vector{Any}(undef, length(AR))
property_history = Vector{Any}(undef, length(AR))
wake_history = Vector{Any}(undef, length(AR))
CF = Vector{Vector{Vector{Float64}}}(undef, length(AR))
CM = Vector{Vector{Vector{Float64}}}(undef, length(AR))

# non-dimensional time (t*Vinf/c)
t = range(0.0, 10.0, step=1/16)

# chord length
c = 1

# time step
dt = [t[i+1]-t[i] for i = 1:length(t)-1]

for i = 1:length(AR)

    # span length
    b = AR[i]*c

    # planform area
    S = b*c

    # geometry
    xle = [0.0, 0.0]
    yle = [-b/2, b/2]
    zle = [0.0, 0.0]
    chord = [c, c]
    theta = [0.0, 0.0]
    phi = [0.0, 0.0]
    fc = fill((xc) -> 0, 2) # camberline function for each section
    ns = 13
    nc = 4
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false
    symmetric = false

    # reference parameters
    cref = c
    bref = b
    Sref = S
    rref = [0.0, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # freestream parameters
    alpha = 5.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    # create vortex rings
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # create vector containing surfaces
    surfaces = [surface]

    # run analysis
    system[i], surface_history[i], property_history[i], wake_history[i] =
        unsteady_analysis(surfaces, ref, fs, dt; symmetric, wake_finite_core = false)

    # extract forces at each time step
    CF[i], CM[i] = body_forces_history(system[i], surface_history[i],
        property_history[i]; frame=Wind())

end

nothing # hide
```

We can visualize the solution using the `write_vtk` function.
```julia
write_vtk("acceleration-AR4", surface_history[1], property_history[1],
    wake_history[1], dt; symmetric=false)
```

![](acceleration-AR4.gif)

The transient lift and drag coefficients are similar to those shown in Figures 13.34 and 13.35 of Low-Speed Aerodynamics by Katz and Plotkin.

```@example rectangular-wing-sudden-acceleration
using Plots
pyplot()

# lift coefficient plot
plot(
    xlim = (0.0, 10.0),
    xticks = 0.0:1.0:10.0,
    xlabel = "\$ \\frac{U_\\infty t}{c} \$",
    ylim = (0.0, 0.55),
    yticks = 0.0:0.1:0.5,
    ylabel = "\$ C_{L} \$",
    grid = false,
    overwrite_figure=false
    )

for i = 1:length(AR)
    CL = [CF[i][j][3] for j = 1:length(CF[i])]
    plot!(t[2:end], CL, label="AR = $(AR[i])")
end

plot!(show=true)

savefig("rectangular-wing-sudden-acceleration-cl.svg") # hide

nothing # hide
```

![](rectangular-wing-sudden-acceleration-cl.svg)

```@example rectangular-wing-sudden-acceleration
# drag coefficient plot
plot(
    xlim = (0.0, 10.0),
    xticks = 0.0:1.0:10.0,
    xlabel = "\$ \\frac{U_\\infty t}{c} \$",
    ylim = (0.0, 0.030),
    yticks = 0.0:0.005:0.03,
    ylabel = "\$ C_{D} \$",
    grid = false,
    overwrite_figure=false
    )

for i = 1:length(AR)
    CD = [CF[i][j][1] for j = 1:length(CF[i])]
    plot!(t[2:end], CD, label="AR = $(AR[i])")
end

plot!(show=true)

savefig("rectangular-wing-sudden-acceleration-cd.svg") # hide

nothing # hide
```

![](rectangular-wing-sudden-acceleration-cd.svg)

We modeled the problem in the body-fixed reference frame (which for this problem is more straightforward), but we could have also modeled the problem in the global reference frame.

```@example rectangular-wing-sudden-acceleration
# Katz and Plotkin: Figures 13.34 and 13.35
# AR = [4, 8, 12, 20, ∞]
# Vinf*Δt/c = 1/16
# α = 5°

using VortexLattice

AR = [4, 8, 12, 20, 1e3] # last aspect ratio is essentially infinite

system_t = Vector{Any}(undef, length(AR))
surface_history_t = Vector{Any}(undef, length(AR))
property_history_t = Vector{Any}(undef, length(AR))
wake_history_t = Vector{Any}(undef, length(AR))
CF_t = Vector{Vector{Vector{Float64}}}(undef, length(AR))
CM_t = Vector{Vector{Vector{Float64}}}(undef, length(AR))

# non-dimensional time (t*Vinf/c)
t = range(0.0, 10.0, step=1/16)

# chord length
c = 1

# time step
dt = [t[i+1]-t[i] for i = 1:length(t)-1]

for i = 1:length(AR)

    # span length
    b = AR[i]*c

    # planform area
    S = b*c

    # geometry
    xle = [0.0, 0.0]
    yle = [-b/2, b/2]
    zle = [0.0, 0.0]
    chord = [c, c]
    theta = [0.0, 0.0]
    phi = [0.0, 0.0]
    fc = fill((xc) -> 0, 2) # camberline function for each section
    ns = 13
    nc = 4
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false
    symmetric = false

    # reference parameters
    cref = c
    bref = b
    Sref = S
    rref = [0.0, 0.0, 0.0]
    Vinf = 1.0 # reference velocity is 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # freestream parameters
    Vinf = 0.0 # freestream velocity is 0.0
    alpha = 5.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    # create vortex rings
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # create vector containing surfaces at each time step
    surfaces = [[VortexLattice.translate(surface,
        -t[it]*[cos(alpha), 0, sin(alpha)])] for it = 1:length(t)]

    # run analysis
    system_t[i], surface_history_t[i], property_history_t[i], wake_history_t[i] =
        unsteady_analysis(surfaces, ref, fs, dt; symmetric, wake_finite_core = false)

    # extract forces at each time step
    CF_t[i], CM_t[i] = body_forces_history(system_t[i], surface_history_t[i],
        property_history_t[i]; frame=Wind())

end

nothing # hide
```

As can be seen, the transient lift and drag coefficients for the two setups are identical.

```@example rectangular-wing-sudden-acceleration
using Plots
pyplot()

# lift coefficient plot
plot(
    xlim = (0.0, 10.0),
    xticks = 0.0:1.0:10.0,
    xlabel = "\$ \\frac{U_\\infty t}{c} \$",
    ylim = (0.0, 0.55),
    yticks = 0.0:0.1:0.5,
    ylabel = "\$ C_{L} \$",
    grid = false,
    overwrite_figure=false
    )

for i = 1:length(AR)
    CL = [CF_t[i][j][3] for j = 1:length(CF_t[i])]
    plot!(t[2:end], CL, label="AR = $(AR[i])")
end

plot!(show=true)

savefig("moving-rectangular-wing-sudden-acceleration-cl.svg") # hide

nothing # hide
```

![](moving-rectangular-wing-sudden-acceleration-cl.svg)

```@example rectangular-wing-sudden-acceleration
# drag coefficient plot
plot(
    xlim = (0.0, 10.0),
    xticks = 0.0:1.0:10.0,
    xlabel = "\$ \\frac{U_\\infty t}{c} \$",
    ylim = (0.0, 0.030),
    yticks = 0.0:0.005:0.03,
    ylabel = "\$ C_{D} \$",
    grid = false,
    overwrite_figure=false
    )

for i = 1:length(AR)
    CD = [CF_t[i][j][1] for j = 1:length(CF_t[i])]
    plot!(t[2:end], CD, label="AR = $(AR[i])")
end

plot!(show=true)

savefig("moving-rectangular-wing-sudden-acceleration-cd.svg") # hide

nothing # hide
```

![](moving-rectangular-wing-sudden-acceleration-cd.svg)

Visualizing the solution shows the movement of the body in the global reference frame.
```julia
write_vtk("acceleration-AR4-moving", surface_history_t[1], property_history_t[1],
    wake_history_t[1], dt; symmetric=false)
```

![](acceleration-AR4-moving.gif)

For infinite aspect ratios, the problem degenerates into the analysis of the sudden acceleration of a 2D flat plate, for which we have an analytical solution through the work of Herbert Wagner.

```@example rectangular-wing-sudden-acceleration
# See Katz and Plotkin: Figure 13.37
# AR = ∞
# Vinf*Δt/c = 1/16
# α = 5°

# essentially infinite aspect ratio
AR = 1e3

# chord length
c = 1

# span length
b = AR*c

# planform area
S = b*c

# geometry
xle = [0.0, 0.0]
yle = [-b/2, b/2]
zle = [0.0, 0.0]
chord = [c, c]
theta = [0.0, 0.0]*pi/180
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section
ns = 1
nc = 4
spacing_s = Uniform()
spacing_c = Uniform()
mirror = false
symmetric = false

# reference parameters
cref = c
bref = b
Sref = S
rref = [0.0, 0.0, 0.0]
Vinf = 1.0
ref = Reference(Sref, cref, bref, rref, Vinf)

# freestream parameters
alpha = 5.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha, beta, Omega)

# non-dimensional time (t*Vinf/c)
t = range(0.0, 7.0, step=1/8)

# time step
dt = [(t[i+1]-t[i]) for i = 1:length(t)-1]

# create vortex rings
grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

# create vector containing all surfaces
surfaces = [surface]

# run steady analysis
system = steady_analysis(surfaces, ref, fs; symmetric)

# extract steady forces
CFs, CMs = body_forces(system; frame=Wind())

# run transient analysis
system, surface_history, property_history, wake_history = unsteady_analysis(
    surfaces, ref, fs, dt; symmetric=symmetric)

# extract transient forces
CF, CM = body_forces_history(system, surface_history, property_history; frame=Wind())

nothing # hide
```

The results from VortexLattice compare very well with the analytical solution provided by Wagner.  As discussed in Low Speed Aerodynamics by Katz and Plotkin, the difference between the curves can be attributed to the finite acceleration rate during the first time step, which increases the lift sharply during the acceleration and then increases it moderately later.

```@example rectangular-wing-sudden-acceleration

# lift coefficient plot
plot(
    xlim = (0.0, 7.0),
    xticks = 0.0:1.0:7.0,
    xlabel = "\$ \\frac{U_\\infty t}{c} \$",
    ylim = (0.0, 1.0),
    yticks = 0.0:0.1:1.0,
    ylabel = "\$ C_{L} \$",
    grid = false,
    overwrite_figure=false
    )

# Computational Results
CL = getindex.(CF, 3)
CLs = getindex(CFs, 3)
plot!(t[2:end], CL./CLs, label="VortexLattice")

# Wagner's Function (using approximation of R. T. Jones)
Φ(t) = 1 - 0.165*exp(-0.045*t) - 0.335*exp(-0.3*t)

plot!(t, Φ.(2*t), label = "Wagner's Function")

plot!(show=true)

savefig("rectangular-wing-sudden-acceleration-wagner.svg") # hide

nothing # hide
```

![](rectangular-wing-sudden-acceleration-wagner.svg)

## Heaving Oscillations of a Rectangular Wing

This example shows how to predict the transient forces and moments for a heaving rectangular wing.

```@example heaving-rectangular-wing
# Katz and Plotkin: Figures 13.38a
# AR = 4
# k = ω*c/(2*Vinf) = [0.5, 0.3, 0.1]
# c = [1.0, 0.6, 0.2]
# α = -5°

using VortexLattice

# forward velocity
Vinf = 1

# angle of attack
alpha = -5*pi/180

# aspect ratio
AR = 4

# chord lengths
c = [1.0, 0.6, 0.2]

# reduced frequency
k = [0.5, 0.3, 0.1]

t = Vector{Vector{Float64}}(undef, length(k))
CF = Vector{Vector{Vector{Float64}}}(undef, length(k))
CM = Vector{Vector{Vector{Float64}}}(undef, length(k))

for i = 1:length(k)

    # span length
    b = AR*c[i]

    # geometry
    xle = [0.0, 0.0]
    yle = [0.0, b/2]
    zle = [0.0, 0.0]
    chord = [c[i], c[i]]
    theta = [0.0, 0.0]
    phi = [0.0, 0.0]
    fc = fill((xc) -> 0, 2) # camberline function for each section
    ns = 13
    nc = 4
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false
    symmetric = true

    # reference parameters
    cref = c[i]
    bref = b
    Sref = b*c[i]
    rref = [0.0, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # angular frequency
    ω = 2*Vinf*k[i]/c[i]

    # time
    t[i] = range(0.0, 9*pi/ω, length = 100)
    dt = t[i][2:end] - t[i][1:end-1]
    dt = Vinf*dt

    # heaving amplitude
    h = 0.1*c[i]

    # use forward and vertical velocity at beginning of each time step
    Xdot = Vinf*cos(alpha)
    Zdot = Vinf*sin(alpha) .- h*cos.(ω*t[i][1:end-1])

    # freestream parameters for each time step
    fs = trajectory_to_freestream(dt; Xdot, Zdot)

    # surface panels
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # create vector containing all surfaces
    surfaces = [surface]

    # run analysis
    system, surface_history, property_history, wake_history = unsteady_analysis(
        surfaces, ref, fs, dt; symmetric=symmetric, nwake = 50)

    # extract forces at each time step (uses instantaneous velocity as reference)
    CF[i], CM[i] = body_forces_history(system, surface_history, property_history; frame=Wind())

end

nothing # hide
```

Plotting the results reveals that the results are similar to the results in Figure 13.34 of Low-Speed Aerodynamic by Katz and Plotkin, which verifies the unsteady vortex lattice method implementation in VortexLattice.

```@example heaving-rectangular-wing
using Plots
pyplot()

# lift coefficient plot
plot(
    xlim = (6*pi, 8*pi),
    xticks = ([6*pi, 13*pi/2, 7*pi, 15*pi/2, 8*pi], ["\$ 0 \$",
        "\$ \\frac{\\pi}{2} \$", "\$ \\pi \$", "\$ \\frac{3\\pi}{2} \$",
        "\$ 2\\pi \$"]),
    xlabel = "\$ ω \\cdot t \$",
    ylim = (-1.0, 0.1),
    yticks = -1.0:0.2:0.0,
    yflip = true,
    ylabel = "\$ C_{L} \$",
    grid = false,
    )

for i = 1:length(k)
    # extract ω
    ω = 2*Vinf*k[i]/c[i]

    # extract ω*t (use time at the beginning of the time step)
    ωt = ω*t[i][1:end-1]

    # extract CL
    CL = [CF[i][it][3] for it = 1:length(t[i])-1]

    plot!(ωt, CL, label="\$ k = \\frac{\\omega c}{2 U_\\infty} = $(k[i]) \$")
end

plot!(show=true)

savefig("heaving-rectangular-wing.svg") # hide

nothing # hide
```

![](heaving-rectangular-wing.svg)

Visualizing the `k=0.5` case in ParaView yields the following animation.

![](heaving-rectangular-wing.gif)
