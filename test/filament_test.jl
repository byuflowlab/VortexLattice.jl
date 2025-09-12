@testset "VortexFilaments induced velocity" begin

# endpoints of the filament
x1 = SVector{3}(-0.4,0.1,0.2)
x2 = SVector{3}(1.0,0.3,-0.1)

# evaluation point
xt = SVector{3}(2.4, 1.2, 0.4)

# filament strength
Gamma = 0.15

# VortexLattice evaluation
finite_core, core_size = true, 0.001
v_VL = VortexLattice.bound_induced_velocity(xt-x1, xt-x2, finite_core, core_size) * Gamma

# VortexFilaments struct
filaments = VortexLattice.VortexFilaments(1)
filaments.velocity[1] = zero(SVector{3,Float64})
filaments.filaments[1] = VortexLattice.Filament(x1, x2, Gamma, core_size)

# VortexFilament buffer
buffer = zeros(VortexLattice.FastMultipole.data_per_body(filaments), 1)
VortexLattice.FastMultipole.source_system_to_buffer!(buffer, 1, filaments, 1)

# VortexFilament evaluation
probes = VortexLattice.FastMultipole.allocate_target_buffer(Float64, filaments)
probes[1:3, 1] .= xt
VortexLattice.FastMultipole.direct!(probes, 1:1, VortexLattice.FastMultipole.DerivativesSwitch(false,true,false), filaments, buffer, 1:1)
v_VF = VortexLattice.FastMultipole.get_gradient(probes, 1)

# @test isapprox(v_VL, v_VF; atol=1e-10)
@test isapprox(v_VL, v_VF; atol=1e-10)

end