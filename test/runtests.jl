using Test
using VortexLatticeMethod

ztol = sqrt(eps())

@testset "AVL - Run 1" begin

    # Simple Wing with Uniform Spacing

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    fc = fill(x->0, 2)
    ns = 12
    nc = 1
    spacing_s = Uniform()
    spacing_c = Uniform()

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rcg = [0.50, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rcg)

    alpha = 1.0*pi/180
    beta = 5.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    # horseshoe vortices with symmetry
    mirror = false
    symmetric = true

    panels = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.24324, atol=1e-3)
    @test isapprox(CD, 0.00243, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02252, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # horseshoe vortices with mirrored geometry
    mirror = true
    symmetric = false

    panels = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s, spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.24324, atol=1e-3)
    @test isapprox(CD, 0.00243, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02252, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings with symmetry, untwisted geometry
    mirror = false
    symmetric = true

    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.24324, atol=1e-3)
    @test isapprox(CD, 0.00243, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02252, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings with mirrored geometry, untwisted geometry
    mirror = true
    symmetric = false

    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.24324, atol=1e-3)
    @test isapprox(CD, 0.00243, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02252, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)


    # vortex rings with symmetry, twisted geometry
    mirror = false
    symmetric = true

    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, fc, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.24324, atol=1e-3)
    @test isapprox(CD, 0.00243, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02252, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings with mirrored geometry, twisted geometry
    mirror = true
    symmetric = false

    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, fc, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.24324, atol=1e-3)
    @test isapprox(CD, 0.00243, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02252, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

end

@testset "AVL - Run 2" begin

    # Run 2: Simple Wing with Cosine Spacing

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    fc = fill(x->0, 2)
    ns = 12
    nc = 1
    spacing_s = Cosine()
    spacing_c = Uniform()
    mirror = true
    symmetric = false

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rcg = [0.50, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rcg)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    # horseshoe vortices
    panels = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.23744, atol=1e-3)
    @test isapprox(CD, 0.00254, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00243, atol=1e-5)
    @test isapprox(Cm, -0.02165, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings, untwisted geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.23744, atol=1e-3)
    @test isapprox(CD, 0.00254, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00243, atol=1e-5)
    @test isapprox(Cm, -0.02165, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings, twisted geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, fc, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.23744, atol=1e-3)
    @test isapprox(CD, 0.00254, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00243, atol=1e-5)
    @test isapprox(Cm, -0.02165, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

end

@testset "AVL - Run 3" begin

    # Simple Wing at High Angle of Attack

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    fc = fill(x->0, 2)
    ns = 12
    nc = 1
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false
    symmetric = true

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rcg = [0.50, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rcg)

    alpha = 8.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    # horseshoe vortices
    panels = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.80348, atol=1e-5)
    @test isapprox(CD, 0.02651, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.02696, atol=1e-5)
    @test isapprox(Cm, -0.07399, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings, untwisted geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.80348, atol=1e-5)
    @test isapprox(CD, 0.02651, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.02696, atol=1e-5)
    @test isapprox(Cm, -0.07399, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings, twisted geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, fc, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.80348, rtol=0.01)
    @test isapprox(CD, 0.02651, rtol=0.01)
    @test isapprox(outputs.CDiff, 0.02696, rtol=0.01)
    @test isapprox(Cm, -0.07399, rtol=0.02)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

end

@testset "AVL - Run 4" begin

    # Simple Wing with Dihedral

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 3.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, atan(zle[2]-zle[1], yle[2]-yle[1])]
    fc = fill(x->0, 2)
    ns = 12
    nc = 1
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false
    symmetric = true

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rcg = [0.50, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rcg)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    # horseshoe vortices
    panels = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.24787, atol=1e-3)
    @test isapprox(CD, 0.00246, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02395, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings, untwisted geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.24787, atol=1e-3)
    @test isapprox(CD, 0.00246, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02395, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings, twisted geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, fc, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    # NOTE: There is some interaction between where dihedral is applied (in this
    # case at the leading edge) and twist which cause the dihedral at the quarter
    # chord in this final case to vary from that in the previous two cases.  This,
    # among other geometrical differences, explain the differences in the results
    # of this test case relative to the previous test cases with dihedral.  Arguably,
    # the results from this test case are more accurate.

    # As a result we use the current results of this package for this test case

    @test isapprox(CL, 0.24102, atol=1e-3)
    @test isapprox(CD, 0.00234, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00232, atol=1e-5)
    @test isapprox(Cm, -0.02334, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

end

@testset "AVL - Run 5" begin

    # Simple Wing with Dihedral at Very High Angle of Attack

    # NOTE: this test case is nonphysical, so it just tests the numerics

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 3.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, atan(zle[2]-zle[1], yle[2]-yle[1])]
    fc = fill(x->0, 2)
    ns = 12
    nc = 1
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false
    symmetric = true

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rcg = [0.50, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rcg)

    alpha = 20.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    # horseshoe vortices
    panels = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 1.70982, atol=0.01)
    @test isapprox(CD, 0.12904, rtol=0.01)
    @test isapprox(outputs.CDiff, 0.11502, rtol=0.01)
    @test isapprox(Cm, -0.45606, rtol=0.01)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings, untwisted geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 1.70982, rtol=0.01)
    @test isapprox(CD, 0.12904, rtol=0.01)
    @test isapprox(outputs.CDiff, 0.11502, rtol=0.01)
    @test isapprox(Cm, -0.45606, rtol=0.01)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings, twisted geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, fc, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    # NOTE: There is some interaction between where dihedral is applied (in this
    # case at the leading edge) and twist which cause the dihedral at the quarter
    # chord in this final case to vary from that in the previous two cases.  This,
    # among other geometrical differences, explain the differences in the results
    # of this test case relative to the previous test cases with dihedral.  Arguably,
    # the results from this test case are more accurate.

    # As a result we use the current results of this package for this test case

    @test isapprox(CL, 1.70982, rtol=0.01)
    @test isapprox(CD, 0.12904, rtol=0.01)
    @test isapprox(outputs.CDiff, 0.11502, rtol=0.01)
    @test isapprox(Cm, -0.45606, rtol=0.015)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

end

@testset "AVL - Run 6" begin

    # Wing and Tail

    # NOTE: AVL's finite-core model is turned off for these tests

    # wing
    xle = [0.0, 0.2]
    yle = [0.0, 5.0]
    zle = [0.0, 1.0]
    chord = [1.0, 0.6]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, atan(zle[2]-zle[1], yle[2]-yle[1])]
    fc = fill(x->0, 2)
    ns = 12
    nc = 1
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
    fc_h = fill(x->0, 2)
    ns_h = 6
    nc_h = 1
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
    fc_v = fill(x->0, 2)
    ns_v = 5
    nc_v = 1
    spacing_s_v = Uniform()
    spacing_c_v = Uniform()
    mirror_v = false

    Sref = 9.0
    cref = 0.9
    bref = 10.0
    rcg = [0.5, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rcg)

    alpha = 5.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    symmetric = true

    # horseshoe vortices
    wing = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    htail = wing_to_horseshoe_vortices(xle_h, yle_h, zle_h, chord_h, theta_h, ns_h, nc_h;
        mirror=mirror_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(htail, [4.0, 0.0, 0.0])

    vtail = wing_to_horseshoe_vortices(xle_v, yle_v, zle_v, chord_v, theta_v, ns_v, nc_v;
        mirror=mirror_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    translate!(vtail, [4.0, 0.0, 0.0])

    vehicle = vcat(wing, htail, vtail)

    outputs = vlm(vehicle, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.60408, atol=1e-4)
    @test isapprox(CD, 0.01058, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.010378, atol=1e-4)
    @test isapprox(Cm, -0.02778, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings, untwisted geometry
    wing = wing_to_vortex_rings(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    htail = wing_to_vortex_rings(xle_h, yle_h, zle_h, chord_h, theta_h, ns_h, nc_h;
        mirror=mirror_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(htail, [4.0, 0.0, 0.0])

    vtail = wing_to_vortex_rings(xle_v, yle_v, zle_v, chord_v, theta_v, ns_v, nc_v;
        mirror=mirror_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    translate!(vtail, [4.0, 0.0, 0.0])

    vehicle = vcat(wing, htail, vtail)

    outputs = vlm(vehicle, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.60408, atol=1e-4)
    @test isapprox(CD, 0.01058, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.010378, atol=1e-4)
    @test isapprox(Cm, -0.02778, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings, twisted geometry
    wing = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, fc, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    htail = wing_to_vortex_rings(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, fc_h, ns_h, nc_h;
        mirror=mirror_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(htail, [4.0, 0.0, 0.0])

    vtail = wing_to_vortex_rings(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, fc_v, ns_v, nc_v;
        mirror=mirror_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    translate!(vtail, [4.0, 0.0, 0.0])

    vehicle = vcat(wing, htail, vtail)

    outputs = vlm(vehicle, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    # NOTE: Properties are lower in this case because the dihedral is higher at
    # the quarter chord based on how we modeled the geometry.

    @test isapprox(CL, 0.60408, rtol=0.01)
    @test isapprox(CD, 0.01058, rtol=0.01)
    @test isapprox(outputs.CDiff, 0.010378, rtol=0.015)
    @test isapprox(Cm, -0.02778, atol=1e-3)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

end

@testset "AVL - Run 7" begin

    # Simple Wing with Chordwise Panels

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    fc = fill(x->0, 2)
    ns = 12
    nc = 6
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rcg = [0.50, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rcg)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    symmetric = true

    # horseshoe vortices
    panels = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.24454, atol=1e-3)
    @test isapprox(CD, 0.00247, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00248, atol=1e-5)
    @test isapprox(Cm, -0.02091, atol=1e-4)
    @test isapprox(CY, 0.0, atol=1e-16)
    @test isapprox(Cl, 0.0, atol=1e-16)
    @test isapprox(Cn, 0.0, atol=1e-16)

    # vortex rings, untwisted geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.24454, atol=1e-3)
    @test isapprox(CD, 0.00247, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00248, atol=1e-5)
    @test isapprox(Cm, -0.02091, atol=1e-4)
    @test isapprox(CY, 0.0, atol=1e-16)
    @test isapprox(Cl, 0.0, atol=1e-16)
    @test isapprox(Cn, 0.0, atol=1e-16)

    # vortex rings, twisted_geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, fc, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.24454, atol=1e-3)
    @test isapprox(CD, 0.00247, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.00248, atol=1e-5)
    @test isapprox(Cm, -0.02091, atol=1e-4)
    @test isapprox(CY, 0.0, atol=1e-16)
    @test isapprox(Cl, 0.0, atol=1e-16)
    @test isapprox(Cn, 0.0, atol=1e-16)

end

@testset "AVL - Run 8" begin

    # Simple Wing with Cosine-Spaced Spanwise and Chordwise Panels

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    fc = fill(x->0, 2)
    ns = 12
    nc = 6
    spacing_s = Cosine()
    spacing_c = Cosine()
    mirror = false

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rcg = [0.50, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rcg)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    symmetric = true

    # horseshoe vortices
    panels = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.23879, atol=1e-3)
    @test isapprox(CD, 0.00249, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.0024626, atol=1e-5)
    @test isapprox(Cm, -0.01995, atol=1e-4)
    @test isapprox(CY, 0.0, atol=1e-16)
    @test isapprox(Cl, 0.0, atol=1e-16)
    @test isapprox(Cn, 0.0, atol=1e-16)

    # vortex rings, untwisted geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.23879, atol=1e-3)
    @test isapprox(CD, 0.00249, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.0024626, atol=1e-5)
    @test isapprox(Cm, -0.01995, atol=1e-4)
    @test isapprox(CY, 0.0, atol=1e-16)
    @test isapprox(Cl, 0.0, atol=1e-16)
    @test isapprox(Cn, 0.0, atol=1e-16)

    # vortex rings, twisted_geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, fc, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    outputs = vlm(panels, ref, fs, symmetric)
    CD, CY, CL = outputs.CF
    Cl, Cm, Cn = outputs.CM

    @test isapprox(CL, 0.23879, atol=1e-3)
    @test isapprox(CD, 0.00249, atol=1e-5)
    @test isapprox(outputs.CDiff, 0.0024626, atol=1e-5)
    @test isapprox(Cm, -0.01995, atol=1e-4)
    @test isapprox(CY, 0.0, atol=1e-16)
    @test isapprox(Cl, 0.0, atol=1e-16)
    @test isapprox(Cn, 0.0, atol=1e-16)

end

@testset "AVL - Run 9" begin

    # Stability-axis derivatives...
    #
    #                             alpha                beta
    #                  ----------------    ----------------
    # y  force CY |    CYa =  -0.000011    CYb =  -0.000006
    # x' mom.  Cl'|    Cla =  -0.040657    Clb =  -0.025358
    # z' mom.  Cn'|    Cna =   0.002904    Cnb =   0.000459
    #
    #                     roll rate  p'      pitch rate  q'        yaw rate  r'
    #                  ----------------    ----------------    ----------------
    # z' force CL |    CLp =  -0.000002                        CLr =  -0.000206
    # y  force CY |    CYp =   0.046821    CYq =  -0.000013    CYr =  -0.000742
    # x' mom.  Cl'|                        Clq =  -0.048812    Clr =   0.063998
    # y  mom.  Cm |    Cmp =   0.000000                        Cmr =   0.000026
    # z' mom.  Cn'|    Cnp =  -0.019770    Cnq =   0.000901    Cnr =  -0.000895
    #

    # Run 9: Simple Wing Stability Derivatives

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    fc = fill(x->0, 2)
    ns = 12
    nc = 1
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = true
    symmetric = false

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rcg = [0.50, 0.0, 0.0]
    ref = Reference(Sref, cref, bref, rcg)

    alpha = 1.0*pi/180
    beta = 5.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    # horseshoe vortices
    panels = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    dCF, dCM = stability_derivatives(panels, ref, fs, symmetric)
    CDa, CYa, CLa = dCF.alpha
    CDb, CYb, CLb = dCF.beta
    CDp, CYp, CLp = dCF.p
    CDq, CYq, CLq = dCF.q
    CDr, CYr, CLr = dCF.r
    Cla, Cma, Cna = dCM.alpha
    Clb, Cmb, Cnb = dCM.beta
    Clp, Cmp, Cnp = dCM.p
    Clq, Cmq, Cnq = dCM.q
    Clr, Cmr, Cnr = dCM.r

    @test isapprox(CLa, 4.638088, rtol=0.01)
    @test isapprox(CLb, 0.0, atol=ztol)
    @test isapprox(CYa, 0.0, atol=ztol)
    # @test isapprox(CYb, -0.000007, rtol=0.01) #TODO 0.0
    @test isapprox(Cla, 0.0, atol=ztol)
    # @test isapprox(Clb, 0.025749, rtol=0.01) #TODO 0.0021917056419237884
    @test isapprox(Cma, -0.429247, rtol=0.01)
    @test isapprox(Cmb, 0.0, atol=ztol)
    @test isapprox(Cna, 0.0, atol=ztol)
    # @test isapprox(Cnb, 0.000466, atol=ztol) #TODO 0.0
    @test isapprox(Clp, -0.518725, rtol=0.01)
    @test isapprox(Clq, 0.0, atol=ztol)
    # @test isapprox(Clr, 0.064243, rtol=0.01) #TODO 0.051990368639157805
    @test isapprox(Cmp, 0.0, atol=ztol)
    @test isapprox(Cmq, -0.517094, rtol=0.01)
    @test isapprox(Cmr, 0.0, atol=ztol)
    # @test isapprox(Cnp, -0.019846, rtol=0.01) #TODO -0.025901883268277578
    @test isapprox(Cnq, 0.0, atol=ztol)
    @test isapprox(Cnr, -0.000898, rtol=0.01)

    # vortex rings, untwisted geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    dCF, dCM = stability_derivatives(panels, ref, fs, symmetric)
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

    # vortex rings, twisted geometry
    panels = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, fc, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    dCF, dCM = stability_derivatives(panels, ref, fs, symmetric)
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

end

# # ---- Veldhius validation case ------
# b = 0.64*2
# AR = 5.33
# λ = 1.0
# Λ = 0.0
# ϕ = 0.0
# θr = 0
# θt = 0
# npanels = 50
# duplicate = false
# spacing = "uniform"

# wing = VLM.simplewing(b, AR, λ, Λ, ϕ, θr, θt, npanels, duplicate, spacing)


# import Interpolations: interpolate, Gridded, Linear

# """helper function"""
# function interp1(xpt, ypt, x)
#     intf = interpolate((xpt,), ypt, Gridded(Linear()))
#     y = zeros(x)
#     idx = (x .> xpt[1]) .& (x.< xpt[end])
#     y[idx] = intf[x[idx]]
#     return y
# end

# function votherprop(rpos)

#     # rcenter = [0.0; 0.469*b/2; 0.0]
#     rvec = abs(rpos[2] - 0.469*b/2)  # norm(rpos - rcenter)

#     rprop = [0.017464, 0.03422, 0.050976, 0.067732, 0.084488, 0.101244, 0.118]
#     uprop = [0.0, 1.94373, 3.02229, 7.02335, 9.02449, 8.85675, 0.0]/50.0
#     vprop = [0.0, 1.97437, 2.35226, 4.07227, 4.35436, 3.69232, 0.0]/50.0

#     cw = 1.0

#     u = 2 * interp1(rprop, uprop, [rvec])[1]  # factor of 2 from far-field
#     v = cw * interp1(rprop, vprop, [rvec])[1]

#     if rpos[2] > 0.469*b/2
#         v *= -1
#     end

#     unew = u*cos(alpha) - v*sin(alpha)
#     vnew = u*sin(alpha) + v*cos(alpha)

#     return [unew; 0.0; vnew]
# end

# alpha = 0.0*pi/180
# beta = 0.0
# Omega = [0.0; 0.0; 0.0]
# fs = VLM.Freestream(alpha, beta, Omega, votherprop)

# Sref = 0.30739212
# cref = 0.24015
# bref = b
# rcg = [0.0, 0.0, 0.0]
# ref = VLM.Reference(Sref, cref, bref, rcg)

# symmetric = true
# CF, CM, ymid, zmid, l, cl, dCF, dCM = VLM.run(wing, ref, fs, symmetric)

# alpha = 8.0*pi/180
# fs = VLM.Freestream(alpha, beta, Omega, votherprop)
# CF, CM, ymid, zmid, l2, cl2, dCF, dCM = VLM.run(wing, ref, fs, symmetric)

# # # total velocity in direction of Vinf
# # Vinfeff = zeros(cl)
# # for i = 1:length(ymid)
# #     Vext = votherprop([0.0; ymid[i]; zmid[i]])
# #     Vinfeff[i] = Vinf + Vext[1]*cos(alpha)*cos(beta) + Vext[2]*sin(beta) + Vext[3]*sin(alpha)*cos(beta)
# # end

# using PyPlot
# # figure()
# # plot(ymid, cl2)
# # plot(ymid, cl3)
# # gcf()

# figure()
# plot(ymid, l2)
# plot(ymid, cl2)
# ylim([0.0, 1])
# gcf()

# trapz(ymid/(b/2), cl2)


# # figure()
# # plot(ymid, cl*Vinf./Vinfeff)
# # gcf()

# # yvec = linspace(0, 0.64, 20)
# # axial = zeros(20)
# # swirl = zeros(20)
# # for i = 1:20
# #     V = votherprop([0.0; yvec[i]; 0.0])
# #     axial[i] = V[1]
# #     swirl[i] = V[2]
# # end

# # figure()
# # plot(yvec, axial)
# # plot(yvec, swirl)
# # gcf()
