using VortexLattice, DifferentialEquations

# aspect ratio
AR = 4
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
reference = Reference(Sref, cref, bref, rref, Vinf)

# freestream parameters
alpha = 5.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
freestream = Freestream(Vinf, alpha, beta, Omega)

# create vortex rings
grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

# create vector containing surfaces
surfaces = [surface]

# times at which to shed wake panels
tshed = 0.0:0.25:10.0

# time span over which to solve
tspan = (0.0, 10.0)

# create ODE (with singular mass matrix)
prob = ODEProblem(surfaces, reference, freestream, tshed, tspan)

# solve ODE
sol = solve(prob)










    # extract forces at each time step
    CF[i], CM[i] = body_forces_history(system[i], surface_history[i],
        property_history[i]; frame=Wind())

end

nothing # hide
