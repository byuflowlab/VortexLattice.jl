@testset "FMM induced velocity" begin

using VortexLattice.StaticArrays

function mean(v)
    return sum(v) / length(v)
end

# wing
xle = [0.0, 0.2]
yle = [0.0, 1.0]
zle = [0.1, 0.0]
chord = [1.0, 0.9]
theta = [0.0, 0.0]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section
ns = 2
nc = 1
spacing_s = Uniform()
spacing_c = Uniform()
mirror = true

# reference
S, c, b = 1.0, mean(chord), yle[end]-yle[1]
r = SVector{3}(0.5*(xle[1]+xle[end]), 0, 0.5*(zle[1]+zle[end]))
V = 1.0
ref = Reference(S, c, b, r, V)

# create system
w1grid, w1ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc; fc = fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=true)
grids = [w1grid]
ratios = [w1ratio]
surface_id = [1]
system = System(grids; ratios)
system.reference[] = ref
symmetric = [false for _ in grids]  # no implied symmetry

# freestream
α, β = 5.0, 0.0
W = SVector{3}(0.0, 0.0, 0.0)
fs = Freestream(V, α, β, W)

# run steady analysis then initialize wake shedding locations
# (FMM wake interface panels require wake_shedding_locations to be set)
steady_analysis!(system, ref, fs; symmetric, surface_id=surface_id)
VortexLattice.update_wake_shedding_locations!(system.wakes, system.wake_shedding_locations,
    system.surfaces, ref, fs, 1.0, nothing, nothing, system.nwake, 0.2)

# induced velocity at some arbitrary point
xt = SVector{3}(0.4, 0.2, 3.4)
wsl = system.wake_shedding_locations[1]

# VortexLattice induced velocity (with wake shedding locations to match FMM)
v_VL = VortexLattice.induced_velocity(xt, system.surfaces[1], system.Γ;
    trailing_vortices=false, finite_core=true, skip_trailing_edge=false,
    wake_shedding_locations=wsl)

# VortexLattice filaments induced velocity (manual loop matching FMM structure)
v_VL_filaments = zero(v_VL)
i_panel = 0
nc, ns = size(system.surfaces[1])
for i_s in 1:ns
    for i_c in 1:nc
        i_panel += 1
        panel = system.surfaces[1][i_c, i_s]

        # all sides of ring vortex (bottom included; wake interface top cancels it)
        v_VL_filaments += VortexLattice.bound_induced_velocity(xt - panel.rtl, xt - panel.rtr, true, panel.core_size) * system.Γ[i_panel]
        v_VL_filaments += VortexLattice.bound_induced_velocity(xt - panel.rtr, xt - panel.rbr, true, panel.core_size) * system.Γ[i_panel]
        v_VL_filaments += VortexLattice.bound_induced_velocity(xt - panel.rbr, xt - panel.rbl, true, panel.core_size) * system.Γ[i_panel]
        v_VL_filaments += VortexLattice.bound_induced_velocity(xt - panel.rbl, xt - panel.rtl, true, panel.core_size) * system.Γ[i_panel]
    end
end
# wake interface panels: top cancels ring bottom, right+bottom+left are trailing vortex segments
for i_s in 1:ns
    panel = system.surfaces[1][nc, i_s]
    Γ_wake = system.Γ[nc * i_s]
    cs = panel.core_size
    v_VL_filaments += VortexLattice.bound_induced_velocity(xt - panel.rbl,    xt - panel.rbr,    true, cs) * Γ_wake  # top (cancels ring bottom)
    v_VL_filaments += VortexLattice.bound_induced_velocity(xt - panel.rbr,    xt - wsl[i_s+1],  true, cs) * Γ_wake  # right
    v_VL_filaments += VortexLattice.bound_induced_velocity(xt - wsl[i_s+1],   xt - wsl[i_s],    true, cs) * Γ_wake  # bottom
    v_VL_filaments += VortexLattice.bound_induced_velocity(xt - wsl[i_s],     xt - panel.rbl,   true, cs) * Γ_wake  # left
end

@test isapprox(v_VL, v_VL_filaments; atol=1e-10)

# FastMultipole.direct! induced velocity
probes = VortexLattice.FastMultipole.ProbeSystem(1)
probes.position[1] = xt
VortexLattice.FastMultipole.direct!(probes, system)
v_direct = probes.gradient[1]
@test isapprox(v_VL, v_direct; atol=1e-12)

# FastMultipole.fmm! induced velocity, finite core.
#
# body_to_multipole_vl! expands a SINGULAR vortex filament while direct! applies
# the finite-core regularization, so with core_size > 0 the two describe slightly
# different physics and fmm! cannot converge to the direct result at any expansion
# order -- the error plateaus at roughly |v(finite core) - v(singular)| (~1e-9 for
# core_size=1e-3 at this evaluation point) from about p=12 onward. Test against
# that floor here, and test the multipole path properly at core_size=0 below.
probes.gradient[1] = zero(eltype(probes.gradient))
_, _, _, _, m2l_list, direct_list, _ = VortexLattice.FastMultipole.fmm!(probes, system; expansion_order=20, leaf_size_source=1)
v_fmm = probes.gradient[1]
@assert length(m2l_list) > 0
@test isapprox(v_VL, v_fmm; atol=1e-8)

# FastMultipole.fmm! induced velocity, zero core: now the multipole expansion and
# the direct kernel represent the same physics, so the multipole error must fall
# exponentially with expansion order all the way down to roundoff. This is what
# actually exercises the m2l path -- a fixed tolerance at finite core cannot.
for i in eachindex(system.surfaces)
    VortexLattice.update_surface_panels!(system.surfaces[i], system.grids[i];
        ratios=system.ratios[i], fcore=(c, Δs) -> 0.0)
end
v_VL_singular = VortexLattice.induced_velocity(xt, system.surfaces[1], system.Γ;
    trailing_vortices=false, finite_core=true, skip_trailing_edge=false,
    wake_shedding_locations=wsl)

fmm_errors = Float64[]
for p in (4, 8, 12, 20)
    probes.gradient[1] = zero(eltype(probes.gradient))
    _, _, _, _, m2l_list_p, _, _ = VortexLattice.FastMultipole.fmm!(probes, system;
        expansion_order=p, leaf_size_source=1)
    @assert length(m2l_list_p) > 0
    push!(fmm_errors, VortexLattice.norm(probes.gradient[1] - v_VL_singular))
end

@test issorted(fmm_errors; rev=true) # monotone convergence in expansion order
@test fmm_errors[1] < 1e-5          # p=4
@test fmm_errors[2] < 1e-8          # p=8
@test fmm_errors[3] < 1e-11         # p=12
@test fmm_errors[4] < 1e-14         # p=20, machine precision

end
