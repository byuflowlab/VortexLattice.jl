# Heaving/Pitching Wing Validation

using VortexLattice
import VortexLattice.StaticArrays: SMatrix, SVector

function prepare_simulation(St, heave_amplitude, pitch_amplitude, phase, n_cycles;
        c_per_dt=0.02, # percent chord per timestep to resolve
        b_pivot=1/3,
        c=3.81e-2, Re=1100.0, ν=1e-6,
        AR=20, ns=13, nc=4,
    )

    # infer simulation parameters based on input
    uinf = Re * ν / c
    h0 = heave_amplitude * c
    A = 2 * h0
    θ0 = pitch_amplitude
    ψ = phase # assume radians
    ω = St * uinf / A * 2 * pi # angular frequency
    T = 2*pi / ω
    b = c * AR
    TF = promote_type(typeof(uinf), typeof(h0), typeof(ω), typeof(c))

    # generate wing
    xle = [0.0, 0.0]
    yle = [-b/2, b/2]
    zle = [0.0, 0.0]
    chord = [c, c]
    theta = [0.0, 0.0]
    phi = [0.0, 0.0]
    fc = fill((xc) -> 0, length(yle)) # camberline function for each section
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false
    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # initial orientation
    origin = SVector{3,TF}(b_pivot * c, 0.0, 0.0)
    θ_init = θ0 * sin(ψ)
    Rinit = VortexLattice.Rodrigues(SVector{3}(0,1.0,0), θ_init)
    VortexLattice.rotate!(grid, Rinit, origin)

    grids = [grid]
    ratios = [ratio]
    system = System(grids; ratios)

    # wing reference
    Sref = c*b
    cref = c
    bref = b
    rref = [0.0, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rref, uinf)
    system.reference[] = ref

    # freestream parameters
    alpha = 0.0 * pi/180
    beta = 0.0
    Omega = [0; 0.0; 0.0]
    fs = Freestream(uinf, alpha, beta, Omega)
    system.freestream[] = fs
    uinf_func(t) = SVector{3,Float64}(uinf, 0.0, 0.0)

    # reference frames
    frames = ReferenceFrame(system;
        origin = origin,
        v = SVector{3}(0.0, 0.0, h0 * ω),
        ω_axis = SVector{3}(0.0, 1.0, 0.0),
        ω = θ0 * ω * cos(ψ),
        R = Rinit,
        name = "wing",
        child_index = Int[],
        dependent_index = collect(1:length(system.surfaces))
    )

    # maneuver
    function heaving_maneuver!(frames::Vector{ReferenceFrame{TF}}, system, wake, t) where TF
        # extract current reference frame state
        (; x, ω_axis, R, name, parent_index, child_index, dependent_index) = frames[1]

        # compute rates of change
        this_ω = θ0 * ω * cos(ω*t + ψ)
        this_v = h0 * ω * cos(ω*t)

        # update the frame
        frames[1] = typeof(frames[1])(x, SVector{3,TF}(0.0, 0.0, this_v), ω_axis, this_ω, R, name, parent_index, child_index, dependent_index)
    end

    # timestep parameterization
    dt = c / uinf * c_per_dt
    t_final = T * n_cycles
    t_range = range(0, stop=t_final, step=dt)

    # generate monitors
    monitors = (VortexLattice.ForcesMonitor(length(t_range); frame=Wind()),)

    # force system evaluation to avoid NaNs in initialization and to speed up convergence
    steady_analysis!(system, system.reference[], system.freestream[]; symmetric=false, trailing_vortices=true)

    return system, frames, heaving_maneuver!, uinf_func, t_range, monitors
end

# establish parameters
St = 0.3
heave_amplitude = 0.25 # chord fraction
pitch_amplitude = 15 * pi/180 # in radians
phase = 90 * pi/180 # in radians
n_cycles = 5
eta = 0.3

# generate simulation inputs
system, frames, maneuver, uinf, t_range, monitors = prepare_simulation(St, heave_amplitude, pitch_amplitude, phase, n_cycles;
    c_per_dt=0.02,
    AR=6,
    nc=1
)

# run simulation
wake = simulate!(system, frames, maneuver, uinf, t_range;
                    name = "heaving", eta,
                    fmm_wake_args=(leaf_size_source=20,),
                    fmm_vehicle_args=(leaf_size_source=20,),
                    # particle_trailing_methods=fill(VortexLattice.NoShed(), length(system.surfaces)),
                    particle_trailing_methods=fill(VortexLattice.OverlapPPS(1.3,1), length(system.surfaces)),
                    # particle_unsteady_methods=fill(VortexLattice.NoShed(), length(system.surfaces)),
                    particle_unsteady_methods=fill(VortexLattice.OverlapPPS(1.3,5), length(system.surfaces)),
                    monitors
                )

