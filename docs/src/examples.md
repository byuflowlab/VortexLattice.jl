# Examples

These examples show how to use VortexLattice for various geometries, flow conditions, and analyses.

```@contents
Pages = ["examples.md"]
Depth = 3
```

## Symmetric Simple Wing

```@example simple-symmetric
using VortexLattice

# geometry
xle = [0.0, 0.4]
yle = [0.0, 7.5]
zle = [0.0, 0.0]
chord = [2.2, 1.8]
theta = [2.0*pi/180, 2.0*pi/180]
phi = [0.0, 0.0]

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
ref = Reference(Sref, cref, bref, rref)

# freestream parameters
alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
vother = nothing
fs = Freestream(alpha, beta, Omega, vother)

# construct surface
surface = wing_to_panels(xle, yle, zle, chord, theta, phi, ns, nc; spacing_s=spacing_s, spacing_c=spacing_c)

# declare symmetry
symmetric = true

# perform steady state analysis
system = steady_analysis(surface, ref, fs; symmetric=symmetric)

# retrieve near-field forces
CF, CM = body_forces(system, surface, ref, fs; symmetric=symmetric, frame=Wind())

# perform far-field analysis
CDiff = far_field_drag(system, surface, ref, fs; symmetric=symmetric)

CD, CY, CL = CF
Cl, Cm, Cn = CM

nothing #hide
```

![](simple-example.png)

## Mirrored Simple Wing

```@example simple-mirrored
using VortexLattice

# geometry
xle = [0.0, 0.4]
yle = [0.0, 7.5]
zle = [0.0, 0.0]
chord = [2.2, 1.8]
theta = [2.0*pi/180, 2.0*pi/180]
phi = [0.0, 0.0]

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
ref = Reference(Sref, cref, bref, rref)

# freestream parameters
alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
vother = nothing
fs = Freestream(alpha, beta, Omega, vother)

# construct surface (and mirror geometry)
surface = wing_to_panels(xle, yle, zle, chord, theta, phi, ns, nc; spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)

# declare symmetry
symmetric = false

# perform steady state analysis
system = steady_analysis(surface, ref, fs; symmetric=symmetric)

# retrieve near-field forces
CF, CM = body_forces(system, surface, ref, fs; symmetric=symmetric, frame=Wind())

# perform far-field analysis
CDiff = far_field_drag(system, surface, ref, fs; symmetric=symmetric)

CD, CY, CL = CF
Cl, Cm, Cn = CM

nothing #hide
```

![](simple-example.png)

## Simple Wing with Dihedral

```@example simple-dihedral
using VortexLattice

xle = [0.0, 0.4]
yle = [0.0, 7.5]
zle = [0.0, 3.0]
chord = [2.2, 1.8]
theta = [2.0*pi/180, 2.0*pi/180]
phi = [0.0, 0.0]
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
ref = Reference(Sref, cref, bref, rref)

alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
vother = nothing
fs = Freestream(alpha, beta, Omega, vother)

# declare symmetry
symmetric = true

# construct surface
surface = wing_to_panels(xle, yle, zle, chord, theta, phi, ns, nc; spacing_s=spacing_s, spacing_c=spacing_c)

# perform steady state analysis
system = steady_analysis(surface, ref, fs; symmetric=symmetric)

# retrieve near-field forces
CF, CM = body_forces(system, surface, ref, fs; symmetric=symmetric, frame=Wind())

# perform far-field analysis
CDiff = far_field_drag(system, surface, ref, fs; symmetric=symmetric)

CD, CY, CL = CF
Cl, Cm, Cn = CM

nothing #hide
```

![](simple-dihedral.png)

## Wing and Tail

```@example wing-tail
using VortexLattice

# wing
xle = [0.0, 0.2]
yle = [0.0, 5.0]
zle = [0.0, 1.0]
chord = [1.0, 0.6]
theta = [2.0*pi/180, 2.0*pi/180]
phi = [0.0, 0.0]
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
ns_v = 5
nc_v = 3
spacing_s_v = Uniform()
spacing_c_v = Uniform()
mirror_v = false

Sref = 9.0
cref = 0.9
bref = 10.0
rref = [0.5, 0.0, 0.0]
ref = Reference(Sref, cref, bref, rref)

alpha = 5.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
vother = nothing
fs = Freestream(alpha, beta, Omega, vother)

symmetric = [true, true, false]

# generate vortex rings
wing = wing_to_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

htail = wing_to_panels(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
    mirror=mirror_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
translate!(htail, [4.0, 0.0, 0.0])

vtail = wing_to_panels(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
    mirror=mirror_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
translate!(vtail, [4.0, 0.0, 0.0])

surfaces = [wing, htail, vtail]
surface_id = [1, 2, 3]

system = steady_analysis(surfaces, ref, fs; symmetric=symmetric, surface_id=surface_id)

CF, CM = body_forces(system, surfaces, ref, fs; symmetric=symmetric, frame=Stability())

CDiff = far_field_drag(system, surfaces, ref, fs; symmetric=symmetric)

CD, CY, CL = CF
Cl, Cm, Cn = CM

nothing #hide
```

![](wing-tail.png)

## Body/Stability Derivatives

```@example stability
using VortexLattice

xle = [0.0, 0.4]
yle = [0.0, 7.5]
zle = [0.0, 0.0]
chord = [2.2, 1.8]
theta = [2.0*pi/180, 2.0*pi/180]
phi = [0.0, 0.0]

ns = 12
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()

mirror = true
symmetric = false

Sref = 30.0
cref = 2.0
bref = 15.0
rref = [0.50, 0.0, 0.0]
ref = Reference(Sref, cref, bref, rref)

alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
vother = nothing
fs = Freestream(alpha, beta, Omega, vother)

# generate vortex rings
surface = wing_to_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

system = steady_analysis(surface, ref, fs; symmetric=symmetric)

dCFb, dCMb = body_derivatives(system, surface, ref, fs; symmetric = symmetric)

CXu, CYu, CZu = dCFb.u
CXv, CYv, CZv = dCFb.v
CXw, CYw, CZw = dCFb.w
CXp, CYp, CZp = dCFb.p
CXq, CYq, CZq = dCFb.q
CXr, CYr, CZr = dCFb.r

Clu, Cmu, Cnu = dCMb.u
Clv, Cmv, Cnv = dCMb.v
Clw, Cmw, Cnw = dCMb.w
Clp_b, Cmp_b, Cnp_b = dCMb.p
Clq_b, Cmq_b, Cnq_b = dCMb.q
Clr_b, Cmr_b, Cnr_b = dCMb.r

dCFs, dCMs = stability_derivatives(system, surface, ref, fs; symmetric = symmetric)

CDa, CYa, CLa = dCFs.alpha
CDb, CYb, CLb = dCFs.beta
CDp, CYp, CLp = dCFs.p
CDq, CYq, CLq = dCFs.q
CDr, CYr, CLr = dCFs.r

Cla, Cma, Cna = dCMs.alpha
Clb, Cmb, Cnb = dCMs.beta
Clp_s, Cmp_s, Cnp_s = dCMs.p
Clq_s, Cmq_s, Cnq_s = dCMs.q
Clr_s, Cmr_s, Cnr_s = dCMs.r

nothing #hide
```

![](simple-example.png)
