# Advanced Use Cases

The time domain runs on [FLOWUnsteadyCore](https://github.com/byuflowlab/FLOWUnsteadyCore.jl): the wake, the time loop, the coupling with other bodies, checkpoints and ParaView output are the core's, and VortexLattice provides its surfaces as the core's bodies. The examples below show the pieces a vortex-lattice user needs; the core's documentation covers every keyword of `simulate!`.

## The surfaces as the core's bodies

`surface_bodies(system)` returns one body and one state per surface of a `VortexLatticeSystem`. They are handles into the system: the geometry, the circulation and the panel forces stay in it, and moving a body moves the system's grid. All surfaces of one system are solved together by the system's own influence matrix, and the core never splits them. `system.reference[]` must be set before the run (its velocity nondimensionalizes the forces).

```julia
using VortexLattice, StaticArrays

system = VortexLatticeSystem(grids; ratios)
system.reference[] = Reference(Sref, cref, bref, rref, Vinf)
bodies, states = surface_bodies(system)

frames = ReferenceFrame(system; origin = SVector(0.0, 0.0, 0.0))
freestream = Freestream(SVector(Vinf*cos(alpha), 0.0, Vinf*sin(alpha)))
wake = PanelParticleWake(Float64, bodies; nwakerows = 2, max_particles = 20_000)
monitor = ForcesMonitor(system, nsteps; frame = Wind())

simulate!(states, bodies, frames, wake, freestream; nsteps, dt, monitor)
```

The influence matrix depends only on the surfaces' relative geometry. It is built once and reused while the surfaces move together (the trailing-edge columns follow the wake's shedding line every step, and the factorization is redone only when they change), so a rigid vehicle in steady flight factors it once per run.

## The wake: rings and particles

The wake is a buffer of `nwakerows` rows of vortex rings behind each trailing edge; older rows become vortex particles. With `nwakerows` at least the number of steps nothing is converted and the method is the classical unsteady vortex lattice method (see the examples). With fewer rows the particles carry the far wake, which keeps long runs affordable.

```julia
wake = PanelParticleWake(Float64, bodies; nwakerows = 2, max_particles = 50_000,
    method_trailing = OverlapPPS(1.3, 2), method_unsteady = OverlapPPS(1.3, 2), shed_eta = 0.25)
```

- `method_trailing`, `method_unsteady`: how the spanwise and the temporal circulation differences of a converted row become particles (`NoShed`, `SigmaPPS`, `OverlapPPS`, `SigmaOverlap`).
- `OverlapPPS(overlap, n)` sizes each particle's core from the length of the edge it comes from. On a long spanwise edge (a wide strip) the cores of the shed vortex become large and smear it, which makes the lift respond too quickly to a change; cap them with `OverlapPPS(1.3, 2; max_sigma = Vinf*dt)` or keep more rows as rings.
- `shed_eta`: the fraction of one convection step behind the trailing edge at which the wake starts.

## Viscous correction with sectional polars

Pass one polar per spanwise strip of each surface (anything with a [`VortexLattice.sectional_coefficients`](@ref) method, such as a `CCBlade.AlphaAF`). After each solve the strips' lift is moved to the polar's lift at the strip's effective angle and the polar's drag is added along the local flow. Reynolds and Mach numbers are not computed yet (they are passed as zero).

```julia
using CCBlade
af = CCBlade.AlphaAF("rae101_12.dat"; radians = false)
bodies, states = surface_bodies(system; polars = [fill(af, ns)])
```

`basis` gives each surface's chordwise, spanwise and normal axes in global coordinates (matrix columns); the default identity suits a wing laid out along y. The axes turn with the surface.

## Checkpoints and restart

```julia
# write a checkpoint every 10 steps
simulate!(states, bodies, frames, wake, freestream; nsteps = 100, dt, checkpoint = "output/run", checkpoint_every = 10)

# later: build the same objects again and continue from step 50
simulate!(states, bodies, frames, wake, freestream; nsteps = 100, dt, restart = "output/run_00000050.jls")
```

A restart reproduces the uninterrupted run.

## Fluid domain

The velocity and vorticity on a rectilinear grid, sampled during the run (and written to ParaView when `path` is given):

```julia
using VortexLattice.FLOWUnsteadyCore: FluidDomain
fd = FluidDomain(range(-1.0, 5.0; step = 0.1), range(-2.0, 2.0; step = 0.1), [0.0]; every = 10, path = "output/fluid")
simulate!(states, bodies, frames, wake, freestream; nsteps, dt, fluid_domain = fd)
```

## Surfaces together with other bodies

Any body that implements the core's contract can share the simulation: for example a LiftingLines rotor ahead of a wing. The core solves each package's bodies with that package's solver and iterates between them until the circulation of every body converges.

```julia
using LiftingLines
rotor = LiftingLines.Rotor(data; ns = 12, twist_sign = 1, clockwise = true)

geometries = Any[rotor, wing_bodies...]
states = VortexLattice.FLOWUnsteadyCore.AbstractBodyState[LiftingLines.LiftingLinesSystem(Float64, LiftingLines.n_sections(rotor)), wing_states...]
frames = VortexLattice.FLOWUnsteadyCore.rotor_frames(geometries; omega = -omega)   # the rotor spins, the wing rides the root frame
wake = PanelParticleWake(Float64, geometries; max_particles = 200_000)

simulate!(states, geometries, frames, wake, freestream; nsteps, dt, monitor)
```
