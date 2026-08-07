using VortexLattice
using StaticArrays
using Plots
using DelimitedFiles

"""
Suddenly-started high-AR wing (impulsive start into steady level flight).

Ryan's diagnostic: `eta` should only affect the transient (unsteady) part of
the solution, never the converged steady CL. This script sweeps `eta` for two
wake models — the hybrid `PanelParticleWake` (nwakerows>0, panels near the TE,
particles beyond) and the restored pure-particle wake (nwakerows=0, every
shed row becomes particles immediately, no panel buffer at all) — and checks
whether CL(t)/CL_steady converges to the same long-time value regardless of
eta for each wake model.
"""

constant_maneuver!(frames, system, wake, t) = false

# --- high-AR flat wing (approximates the 2D Wagner problem) ---
AR = 100
c = 1.0
b = AR * c
xle = [0.0, 0.0]
yle = [-b/2, b/2]
zle = [0.0, 0.0]
chord = [c, c]
theta = [0.0, 0.0]
phi = [0.0, 0.0]
ns = 1
nc = 4
fc = fill((xc) -> 0, length(yle))
spacing_s = Uniform()
spacing_c = Uniform()
mirror = false

grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

Sref = c*b
cref = c
bref = b
rref = [0.0, 0.0, 0.0]
Vinf = 1.0
alpha = 5.0 * pi/180
beta = 0.0
Omega = [0.0, 0.0, 0.0]

function build_system(nwakerows=0)
    system = System([grid]; ratios=[ratio], nw=[nwakerows])
    for isurf in eachindex(system.surfaces)
        VortexLattice.update_surface_panels!(system.surfaces[isurf], system.grids[isurf];
            ratios=system.ratios[isurf], fcore=(c, Δs) -> system.core_size)
    end
    system.reference[] = Reference(Sref, cref, bref, rref, Vinf)
    system.freestream[] = Freestream(Vinf, alpha, beta, Omega)
    return system
end

# steady reference lift (same for every eta / wake model)
steady_system = build_system()
steady_analysis!(steady_system, steady_system.reference[], steady_system.freestream[];
    symmetric=false, trailing_vortices=true)
CF_steady, _ = body_forces(steady_system; frame=Wind())
CL_steady = CF_steady[3]
println("CL_steady = ", CL_steady)

# --- unsteady setup shared across runs ---
Vinf_func(t) = SVector{3,Float64}(cos(alpha), 0.0, sin(alpha)) * Vinf
Ωinf_func(t) = SVector{3,Float64}(0.0, 0.0, 0.0)
dt = c / Vinf / 16
max_Vinf_t_c = 10.0
tmax = max_Vinf_t_c * cref / Vinf
nsteps = Int(ceil(tmax / dt))
t_range = range(start=0.0, stop=nsteps*dt, length=nsteps+1)
tstar = Vinf * collect(t_range) / cref

etas = [0.1, 0.2, 0.3, 0.5, 0.8]

function run_case(nwakerows, eta)
    system = build_system(nwakerows)

    frames = ReferenceFrame(system;
        origin=SVector{3}(0.0, 0.0, 0.0),
        v=SVector{3}(0.0, 0.0, 0.0),
        ω_axis=SVector{3}(1.0, 0.0, 0.0),
        ω=0.0,
        R=SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="vehicle",
        child_index=Int[],
        dependent_index=collect(1:length(system.surfaces)),
    )

    monitors = (VortexLattice.ForcesMonitor(length(t_range); frame=Wind()),)

    wake = simulate!(system, frames, constant_maneuver!, Vinf_func, t_range, Ωinf_func;
        wake_type=PanelParticleWake,
        nwakerows=nwakerows,
        max_particles=50_000,
        eta=eta,
        method_trailing=VortexLattice.OverlapPPS(1.3, 2),
        method_unsteady=VortexLattice.OverlapPPS(1.3, 2),
        monitors=monitors,
        derivatives=false,
        name="ssw_nwakerows$(nwakerows)_eta$(eta)",
        path=nothing,
        verbose=false,
    )

    CLs = [monitors[1].CF[i][3] for i in eachindex(monitors[1].CF)]
    return CLs
end

results = Dict{Tuple{Int,Float64}, Vector{Float64}}()
for nwakerows in (0, 4)
    for eta in etas
        println("running nwakerows=$nwakerows, eta=$eta ...")
        results[(nwakerows, eta)] = run_case(nwakerows, eta)
    end
end

# --- convergence check: does CL(t->large)/CL_steady agree across eta? ---
tail_frac = 0.1  # last 10% of the time history
for nwakerows in (0, 4)
    println("\n--- nwakerows = $nwakerows ---")
    for eta in etas
        CLs = results[(nwakerows, eta)]
        n_tail = max(1, Int(round(tail_frac * length(CLs))))
        CL_tail_mean = sum(CLs[end-n_tail+1:end]) / n_tail
        println("  eta=$eta : CL_tail/CL_steady = $(CL_tail_mean / CL_steady)")
    end
end

# analytic solution (Wagner's function)
Φ(t) = 1 - 0.165*exp(-0.045*t) - 0.335*exp(-0.3*t)

p0 = plot(xlabel="t*", ylabel="CL / CL_steady", title="pure particle wake (nwakerows=0)",
    ylim=(0.0, 1.2), xlim=(0.0, max_Vinf_t_c), legend=:bottomright, legendfontsize=6)
p1 = plot(xlabel="t*", ylabel="CL / CL_steady", title="panel-particle wake (nwakerows=4)",
    ylim=(0.0, 1.2), xlim=(0.0, max_Vinf_t_c), legend=:bottomright, legendfontsize=6)
for eta in etas
    plot!(p0, tstar, results[(0, eta)] ./ CL_steady, label="eta=$eta")
    plot!(p1, tstar, results[(4, eta)] ./ CL_steady, label="eta=$eta")
end
for p in (p0, p1)
    plot!(p, tstar, Φ.(tstar .* 2), linestyle=:dot, color=:black, label="Wagner")
end
fig = plot(p0, p1, layout=(1,2), size=(1000,450))

outdir = joinpath(@__DIR__, "suddenly_started_wing_output")
mkpath(outdir)
savefig(fig, joinpath(outdir, "eta_sweep.png"))

open(joinpath(outdir, "eta_sweep.csv"), "w") do io
    header = ["tstar"]
    for nwakerows in (0, 4), eta in etas
        push!(header, "CL_nwakerows$(nwakerows)_eta$(eta)")
    end
    writedlm(io, [header], ',')
    for i in eachindex(tstar)
        row = [tstar[i]]
        for nwakerows in (0, 4), eta in etas
            push!(row, results[(nwakerows, eta)][i])
        end
        writedlm(io, [row], ',')
    end
end

println("\nfigure written to: ", joinpath(outdir, "eta_sweep.png"))
println("data written to:   ", joinpath(outdir, "eta_sweep.csv"))
