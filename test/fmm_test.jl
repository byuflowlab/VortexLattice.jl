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
ns = 1
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
grids = [w1grid]#, w2grid]#, hgrid, vgrid]#, p1grid1, p1grid2, p2grid1, p2grid2]
ratios = [w1ratio]#, w2ratio]#, hratio, vratio]#, p1ratio1, p1ratio2, p2ratio1, p2ratio2]
surface_id = [1]#, 1]#, 3, 4]#, 5, 6, 7, 8]
system = System(grids; ratios)
system.reference[] = ref
symmetric = [false for _ in grids]  # no implied symmetry

# freestream
α, β = 5.0, 0.0
W = SVector{3}(0.0, 0.0, 0.0)
fs = Freestream(V, α, β, W)

# force generation of surface panels
steady_analysis!(system, ref, fs; symmetric, surface_id=surface_id)

# induced velocity at some arbitrary points
xt = SVector{3}(0.4, 0.2, 3.4)

# VortexLattice induced velocity
v_VL = VortexLattice.induced_velocity(xt, system.surfaces[1], system.Γ; trailing_vortices=false, finite_core=true, skip_trailing_edge = true)

# VortexLattice filaments induced velocity
v_VL_filaments = zero(v_VL)
i_panel = 0
nc, ns = size(system.surfaces[1])
for i_s in 1:ns
    for i_c in 1:nc
        # increment panel index
        i_panel += 1

        # get panel
        panel = system.surfaces[1][i_c, i_s]

        # induced velocity from bound vortices
        v_VL_filaments += VortexLattice.bound_induced_velocity(xt - panel.rtl, xt - panel.rtr, true, panel.core_size) * system.Γ[i_panel]
        v_VL_filaments += VortexLattice.bound_induced_velocity(xt - panel.rtr, xt - panel.rbr, true, panel.core_size) * system.Γ[i_panel]
        if i_c < nc
            v_VL_filaments += VortexLattice.bound_induced_velocity(xt - panel.rbr, xt - panel.rbl, true, panel.core_size) * system.Γ[i_panel]
        end
        v_VL_filaments += VortexLattice.bound_induced_velocity(xt - panel.rbl, xt - panel.rtl, true, panel.core_size) * system.Γ[i_panel]
    end
end

@test isapprox(v_VL, v_VL_filaments; atol=1e-10)

# FastMultipole.direct! induced velocity
probes = VortexLattice.FastMultipole.ProbeSystem(1)
probes.position[1] = xt
VortexLattice.FastMultipole.direct!(probes, system)
v_direct = probes.gradient[1]
@test isapprox(v_VL, v_direct; atol=1e-12)

# FastMultipole.fmm! induced velocity
probes.gradient[1] = zero(eltype(probes.gradient))
_, _, _, _, m2l_list, direct_list, _ = VortexLattice.FastMultipole.fmm!(probes, system; expansion_order=20, leaf_size_source=1)#, farfield=false)
v_fmm = probes.gradient[1]
@assert length(m2l_list) > 0
@test isapprox(v_VL, v_fmm; atol=1e-12)

end
