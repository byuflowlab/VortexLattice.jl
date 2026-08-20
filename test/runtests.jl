using Test
using VortexLattice
using LinearAlgebra
using StaticArrays

ztol = sqrt(eps())

# reflects a vector across the x-z plane
flipy(x) = [x[1], -x[2], x[3]]

# constructs a normal vector the way AVL does
function avl_normal_vector(dr, theta)
    st, ct = sincos(theta)
    bhat = dr/norm(dr) # bound vortex vector
    shat = [0, -dr[3], dr[2]]/sqrt(dr[2]^2+dr[3]^2) # chordwise strip normal vector
    chat = [ct, -st*shat[2], -st*shat[3]] # camberline vector
    ncp = cross(chat, dr) # normal vector perpindicular to camberline and bound vortex
    return ncp / norm(ncp) # normal vector used by AVL
end

@testset "AVL - Run 1 - Wing with Uniform Spacing" begin

    # Simple Wing with Uniform Spacing

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    ns = 12
    nc = 1
    spacing_s = Uniform()
    spacing_c = Uniform()

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

    # adjust chord length so x-chord length matches AVL
    chord = @. chord/cos(theta)

    # vortex rings with symmetry
    mirror = false
    symmetric = true

    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grids = [grid]
    ratios = [ratio]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric)

    CF, CM = body_forces(system; frame=Stability())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.24324, atol=1e-3)
    @test isapprox(CD, 0.00243, atol=1e-5)
    @test isapprox(CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02252, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings with mirrored geometry
    mirror = true
    symmetric = false

    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grids = [grid]
    ratios = [ratio]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric)

    CF, CM = body_forces(system; frame=Stability())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.24324, atol=1e-3)
    @test isapprox(CD, 0.00243, atol=1e-5)
    @test isapprox(CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02252, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)
end

@testset "AVL - Run 2 - Wing with Cosine Spacing" begin

    # Run 2: Simple Wing with Cosine Spacing

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    ns = 12
    nc = 1
    spacing_s = Cosine()
    spacing_c = Uniform()
    mirror = true
    symmetric = false

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

    # adjust chord length so x-chord length matches AVL
    chord = @. chord/cos(theta)

    # vortex rings
    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grids = [grid]
    ratios = [ratio]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric)

    CF, CM = body_forces(system; frame=Stability())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.23744, atol=1e-3)
    @test isapprox(CD, 0.00254, atol=1e-5)
    @test isapprox(CDiff, 0.00243, atol=1e-5)
    @test isapprox(Cm, -0.02165, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)
end

@testset "AVL - Run 3 - Wing at High Angle of Attack" begin

    # Simple Wing at High Angle of Attack

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    ns = 12
    nc = 1
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

    alpha = 8.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    # adjust chord length so x-chord length matches AVL
    chord = @. chord/cos(theta)

    # vortex rings
    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grids = [grid]
    ratios = [ratio]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric)

    CF, CM = body_forces(system; frame=Stability())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.80348, atol=1e-3)
    @test isapprox(CD, 0.02651, atol=1e-4)
    @test isapprox(CDiff, 0.02696, atol=1e-5)
    @test isapprox(Cm, -0.07399, atol=1e-3)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)
end

@testset "AVL - Run 4 - Wing with Dihedral" begin

    # Simple Wing with Dihedral

    # NOTE: There is some interaction between twist, dihedral, and chordwise
    # position which causes the normal vectors found by AVL to differ from those
    # computed by this package.  We therefore manually overwrite the normal
    # vectors when this occurs in order to get a better comparison.

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 3.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    ns = 24
    nc = 2
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

    # adjust chord length so x-chord length matches AVL
    chord = @. chord/cos(theta)

    # also get normal vector as AVL defines it
    ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

    # vortex rings
    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grid, ratio, surface = grid_to_surface_panels(grid; ratios=ratio)

    for (ip, p) in enumerate(surface)
        # check that our normal vector is approximately the same as AVL's
        @test isapprox(p.ncp, ncp, rtol=0.01)
        # replace our normal vector with AVL's normal vector for this test
        surface[ip] = set_normal(p, ncp)
    end

    grids = [grid]
    ratios = [ratio]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric)

    CF, CM = body_forces(system; frame=Stability())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.24787, atol=0.02)
    @test isapprox(CD, 0.00246, atol=0.02)
    @test isapprox(CDiff, 0.00245, atol=0.02)
    @test isapprox(Cm, -0.02395, atol=0.02)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)
end

@testset "AVL - Run 5 - Wing with Dihedral at Very High Angle of Attack" begin

    # Simple Wing with Dihedral at Very High Angle of Attack

    # NOTE: this test case is nonphysical, so it just tests the numerics

    # NOTE: There is some interaction between twist, dihedral, and chordwise
    # position which causes the normal vectors found by AVL to differ from those
    # computed by this package.  We therefore manually overwrite the normal
    # vectors when this occurs in order to get a better comparison.

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 3.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    ns = 12
    nc = 1
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false
    symmetric = true

    # adjust chord length to match AVL (which uses chord length in the x-direction)
    chord = @. chord/cos(theta)

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rref = [0.50, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    alpha = 20.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

    # vortex rings, untwisted geometry
    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grid, ratio, surface = grid_to_surface_panels(grid; ratios=ratio)

    for (ip, p) in enumerate(surface)
        # check that our normal vector is approximately the same as AVL's
        @test isapprox(p.ncp, ncp, rtol=0.01)
        # replace our normal vector with AVL's normal vector for this test
        surface[ip] = set_normal(p, ncp)
    end

    grids = [grid]
    ratios = [ratio]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric)

    CF, CM = body_forces(system; frame=Stability())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 1.70982, rtol=0.02)
    @test isapprox(CD, 0.12904, rtol=0.02)
    @test isapprox(CDiff, 0.11502, rtol=0.02)
    @test isapprox(Cm, -0.45606, rtol=0.02)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)
end

@testset "AVL - Run 6 - Wing and Tail without Finite Core Model" begin

    # Wing and Tail without Finite Core Model

    # NOTE: AVL's finite-core model is turned off for these tests

    # NOTE: There is some interaction between twist, dihedral, and chordwise
    # position which causes the normal vectors found by AVL to differ from those
    # computed by this package.  We therefore manually overwrite the normal
    # vectors when this occurs in order to get a better comparison.

    # wing
    xle = [0.0, 0.2]
    yle = [0.0, 5.0]
    zle = [0.0, 1.0]
    chord = [1.0, 0.6]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
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
    ns_v = 5
    nc_v = 1
    spacing_s_v = Uniform()
    spacing_c_v = Uniform()
    mirror_v = false

    # adjust chord lengths to match AVL (which uses chord length in the x-direction)
    chord = @. chord/cos(theta)
    chord_h = @. chord_h/cos(theta_h)
    chord_v = @. chord_v/cos(theta_v)

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

    ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

    # vortex rings - finite core deactivated
    wgrid, wratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grid, ratio, wing = grid_to_surface_panels(wgrid; ratios=wratio)

    for (ip, p) in enumerate(wing)
        # check that our normal vector is approximately the same as AVL's
        @test isapprox(p.ncp, ncp, rtol=0.01)
        # replace our normal vector with AVL's normal vector for this test
        wing[ip] = set_normal(p, ncp)
    end

    hgrid, hratio = wing_to_grid(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(hgrid, [4.0, 0.0, 0.0])

    vgrid, vratio = wing_to_grid(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    translate!(vgrid, [4.0, 0.0, 0.0])

    grids = [wgrid, hgrid, vgrid]
    ratios = [wratio, hratio, vratio]
    surface_id = [1, 1, 1]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric, surface_id=surface_id)

    CF, CM = body_forces(system; frame=Stability())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.60408, atol=1e-2)
    @test isapprox(CD, 0.01058, atol=1e-4)
    @test isapprox(CDiff, 0.010378, atol=1e-3)
    @test isapprox(Cm, -0.02778, atol=2e-3)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)
end

@testset "AVL - Run 7 - Wing and Tail with Finite Core Model" begin

    # Wing and Tail with Finite Core Model

    # NOTE: There is some interaction between twist, dihedral, and chordwise
    # position which causes the normal vectors found by AVL to differ from those
    # computed by this package.  We therefore manually overwrite the normal
    # vectors when this occurs in order to get a better comparison.

    # wing
    xle = [0.0, 0.2]
    yle = [0.0, 5.0]
    zle = [0.0, 1.0]
    chord = [1.0, 0.6]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
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
    ns_v = 5
    nc_v = 1
    spacing_s_v = Uniform()
    spacing_c_v = Uniform()
    mirror_v = false

    # adjust chord lengths to match AVL (which uses chord length in the x-direction)
    chord = @. chord/cos(theta)
    chord_h = @. chord_h/cos(theta_h)
    chord_v = @. chord_v/cos(theta_v)

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

    ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

    # vortex rings - finite core deactivated
    wgrid, wratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grid, ratio, wing = grid_to_surface_panels(wgrid; ratios=wratio)

    for (ip, p) in enumerate(wing)
        # check that our normal vector is approximately the same as AVL's
        @test isapprox(p.ncp, ncp, rtol=0.01)
        # replace our normal vector with AVL's normal vector for this test
        wing[ip] = set_normal(p, ncp)
    end

    hgrid, hratio = wing_to_grid(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(hgrid, [4.0, 0.0, 0.0])

    vgrid, vratio = wing_to_grid(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    translate!(vgrid, [4.0, 0.0, 0.0])

    grids = [wgrid, hgrid, vgrid]
    ratios = [wratio, hratio, vratio]
    surface_id = [1, 2, 3]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric, surface_id=surface_id, derivatives=false)

    CF, CM = body_forces(system; frame=Stability())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.60562, atol=1e-2)
    @test isapprox(CD, 0.01058, atol=1e-4)
    @test isapprox(CDiff, 0.0104855, atol=1e-3)
    # @test isapprox(Cm, -0.03377, atol=2e-3) # Why would this value be different from the previous test?
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)
end

@testset "AVL - Run 8 - Wing with Chordwise Panels" begin

    # Simple Wing with Chordwise Panels

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
    mirror = false

    # adjust chord length to match AVL (which uses chord length in the x-direction)
    chord = @. chord/cos(theta)

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

    symmetric = true

    # vortex rings
    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grids = [grid]
    ratios = [ratio]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric)

    CF, CM = body_forces(system; frame=Stability())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.24454, atol=1e-3)
    @test isapprox(CD, 0.00247, atol=1e-5)
    @test isapprox(CDiff, 0.00248, atol=1e-5)
    @test isapprox(Cm, -0.02091, atol=1e-4)
    @test isapprox(CY, 0.0, atol=1e-16)
    @test isapprox(Cl, 0.0, atol=1e-16)
    @test isapprox(Cn, 0.0, atol=1e-16)
end

@testset "AVL - Run 9 - Wing with Cosine-Spaced Spanwise and Chordwise Panels" begin

    # Simple Wing with Cosine-Spaced Spanwise and Chordwise Panels

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    ns = 12
    nc = 6
    spacing_s = Cosine()
    spacing_c = Cosine()
    mirror = false

    # adjust chord length to match AVL (which uses chord length in the x-direction)
    chord = @. chord/cos(theta)

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

    symmetric = true

    # vortex rings
    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grids = [grid]
    ratios = [ratio]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric)

    CF, CM = body_forces(system; frame=Stability())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.23879, atol=1e-3)
    @test isapprox(CD, 0.00249, atol=1e-5)
    @test isapprox(CDiff, 0.0024626, atol=1e-5)
    @test isapprox(Cm, -0.01995, atol=1e-4)
    @test isapprox(CY, 0.0, atol=1e-16)
    @test isapprox(Cl, 0.0, atol=1e-16)
    @test isapprox(Cn, 0.0, atol=1e-16)
end

@testset "AVL - Run 10 - Wing with Sideslip" begin

    # Simple Wing with Sideslip

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    ns = 12
    nc = 1
    spacing_s = Uniform()
    spacing_c = Uniform()

    # adjust chord length to match AVL (which uses chord length in the x-direction)
    chord = @. chord/cos(theta)

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rref = [0.50, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    alpha = 1.0*pi/180
    beta = 15.0*pi/180
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

    # vortex rings
    mirror = true
    symmetric = false

    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grids = [grid]
    ratios = [ratio]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric)

    CF, CM = body_forces(system; frame=Stability())

    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.22695, atol=1e-3)
    @test isapprox(CD, 0.00227, atol=1e-5)
    @test isapprox(CDiff, 0.0022852, atol=1e-5)
    @test isapprox(Cm, -0.02101, atol=1e-4)
    @test isapprox(CY, 0.0, atol=1e-5)
    @test isapprox(Cl, -0.00644, atol=1e-4)
    @test isapprox(Cn, 0.00012, atol=2e-4)
end

@testset "AVL - Run 11 - Wing Stability Derivatives" begin

    # Run 11: Simple Wing Stability Derivatives

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    ns = 12
    nc = 1
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = true
    symmetric = false

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

    # vortex rings
    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grids = [grid]
    ratios = [ratio]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric)

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

    @test isapprox(CLa, 4.638088, rtol=0.01)
    @test isapprox(CLb, 0.0, atol=ztol)
    @test isapprox(CYa, 0.0, atol=ztol)
    @test isapprox(CYb, -0.000007, atol=1e-4)
    @test isapprox(Cla, 0.0, atol=ztol)
    @test isapprox(Clb, -0.025749, atol=0.001)
    @test isapprox(Cma, -0.429247, rtol=0.01)
    @test isapprox(Cmb, 0.0, atol=ztol)
    @test isapprox(Cna, 0.0, atol=ztol)
    @test isapprox(Cnb, 0.000466, atol=1e-3)
    @test isapprox(Clp, -0.518725, rtol=0.01)
    @test isapprox(Clq, 0.0, atol=ztol)
    @test isapprox(Clr, 0.064243, rtol=0.01)
    @test isapprox(Cmp, 0.0, atol=ztol)
    @test isapprox(Cmq, -0.517094, rtol=0.01)
    @test isapprox(Cmr, 0.0, atol=ztol)
    @test isapprox(Cnp, -0.019846, rtol=0.01)
    @test isapprox(Cnq, 0.0, atol=ztol)
    @test isapprox(Cnr, -0.000898, rtol=0.01)
end

@testset "AVL - Run 12 - Rotational Velocity" begin

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    ns = 12
    nc = 1
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = true
    symmetric = false

    Sref = 30.0
    cref = 2.0
    bref = 15.0
    rref = [0.50, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [2*Vinf*0.05/bref; 0.0; 0.0]  # nondimensional pbar = 0.05
    fs = Freestream(Vinf, alpha, beta, Omega)

    # vortex rings
    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grids = [grid]
    ratios = [ratio]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric)

    CF, CM = body_forces(system; frame=Stability())

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.24323, atol=1e-3)
    @test isapprox(CD, 0.00069, atol=1e-5)
    @test isapprox(Cm, -0.02251, atol=1e-4)
    @test isapprox(CY, 0.00235, atol=2e-4)
    @test isapprox(Cl, -0.02594, atol=1e-4)
    @test isapprox(Cn, -0.00099, atol=2e-4)

end

@testset "Lifting Line Coefficients" begin
    # Simple Wing with Uniform Spacing

    xle = [0.0, 0.4]
    yle = [0.0, 7.5]
    zle = [0.0, 0.0]
    chord = [2.2, 1.8]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
    ns = 12
    nc = 1
    spacing_s = Uniform()
    spacing_c = Uniform()

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

    # adjust chord length so x-chord length matches AVL
    chord = @. chord/cos(theta)

    # vortex rings with symmetry
    mirror = false
    symmetric = true

    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    grids = [grid]
    ratios = [ratio]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric)

    r_ll, c_ll, w_ll = lifting_line_geometry(grids)

    cf, cm = lifting_line_coefficients(system, r_ll, c_ll, w_ll; frame=Stability())

    cl_avl = [0.2618, 0.2646, 0.2661, 0.2664, 0.2654, 0.2628, 0.2584, 0.2513,
        0.2404, 0.2233, 0.1952, 0.1434]
    cd_avl = [0.0029, 0.0024, 0.0023, 0.0023, 0.0023, 0.0023, 0.0024, 0.0024,
        0.0025, 0.0026, 0.0026, 0.0022]
    cm_avl = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    @test isapprox(cf[1][3,:], cl_avl, atol=1e-3, norm=(x)->norm(x, Inf))
    @test isapprox(cf[1][1,:], cd_avl, atol=1e-4, norm=(x)->norm(x, Inf))
    @test isapprox(cm[1][2,:], cm_avl, atol=1e-4, norm=(x)->norm(x, Inf))
end

@testset "Geometry Generation" begin

    # Tests of the geometry generation functions

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

    # Symmetric Geometry

    halfgrid1, ratio1 = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        spacing_s=spacing_s, spacing_c=spacing_c)

    grid, ratio, surface1 = grid_to_surface_panels(halfgrid1; ratios=ratio1)

    grid, ratio, surface2 = grid_to_surface_panels(halfgrid1, ns, nc;
        spacing_s=spacing_s, spacing_c=spacing_c)

    surface3 = similar(surface1)
    VortexLattice.update_surface_panels!(surface3, halfgrid1)

    for I in CartesianIndices(surface1)
        for field in fieldnames(SurfacePanel)
            @test isapprox(getproperty(surface1[I], field), getproperty(surface2[I], field))
        end
    end

    for I in CartesianIndices(surface1)
        for field in fieldnames(SurfacePanel)
            @test isapprox(getproperty(surface1[I], field), getproperty(surface3[I], field))
        end
    end

    # Mirrored Geometry

    mirror = true

    grid1, ratio1 = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        spacing_s=spacing_s, spacing_c=spacing_c, mirror)

    grid, ratio, surface1 = grid_to_surface_panels(halfgrid1; ratios=ratio1, mirror)

    grid, ratio, surface2 = grid_to_surface_panels(halfgrid1, ns, nc;
        spacing_s=spacing_s, spacing_c=spacing_c, mirror)

    surface3 = similar(surface1)
    VortexLattice.update_surface_panels!(surface3, grid1)

    for I in CartesianIndices(surface1)
        for field in fieldnames(SurfacePanel)
            @test isapprox(getproperty(surface1[I], field), getproperty(surface2[I], field))
        end
    end


    for I in CartesianIndices(surface1)
        for field in fieldnames(SurfacePanel)
            @test isapprox(getproperty(surface1[I], field), getproperty(surface3[I], field))
        end
    end


    # Reference line (use 0.5 * chord as reference so middle chorwise node will have x and z = 0)

    xle = zeros(2)
    yle = LinRange(0,4,length(xle))
    zle = zeros(length(xle))
    chord = [2.0, 2.0]
    theta = zeros(length(xle))# .+ 30*pi/180#LinRange(0,30*pi/180,length(xle))
    phi = zeros(length(xle))
    ns = 4
    nc = 2
    spacing_s = Uniform()
    spacing_c = Uniform()

    my_ref = zeros(length(xle),2)
    my_ref[:,1] .= 0.5

    # construct surface
    grid, ratios = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;spacing_s=spacing_s, spacing_c=spacing_c, reference_line=my_ref)

    for p = 1:ns+1
        @test grid[1,2,p] == 0.0
        @test grid[2,2,p] == p-1
        @test grid[3,2,p] == 0.0
    end
end

@testset "Update Trailing Edge Coefficients" begin

    # This test checks whether the function which updates the trailing edge
    # coefficients `update_trailing_edge_coefficients!` results in the same
    # trailing edge coefficients as `influnece_coefficients!`

    # wing
    xle = [0.0, 0.2]
    yle = [0.0, 5.0]
    zle = [0.0, 1.0]
    chord = [1.0, 0.6]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
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
    ns_v = 5
    nc_v = 1
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

    # horseshoe vortices
    wgrid, wratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    hgrid, hratio = wing_to_grid(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(hgrid, [4.0, 0.0, 0.0])

    vgrid, vratio = wing_to_grid(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    translate!(vgrid, [4.0, 0.0, 0.0])

    surfaces = [grid_to_surface_panels(wgrid; ratios=wratio)[3],
                grid_to_surface_panels(hgrid; ratios=hratio)[3],
                grid_to_surface_panels(vgrid; ratios=vratio)[3]]

    surface_id = [1, 2, 2]

    # number of panels
    N = nc*ns + nc_h*ns_h + nc_v*ns_v
    AIC1 = zeros(N, N)

    VortexLattice.influence_coefficients!(AIC1, surfaces;
        symmetric = symmetric,
        trailing_vortices = fill(false, length(surfaces)),
        surface_id = surface_id)

    AIC2 = copy(AIC1)

    VortexLattice.update_trailing_edge_coefficients!(AIC2, surfaces;
        symmetric = symmetric,
        trailing_vortices = fill(false, length(surfaces)),
        surface_id = surface_id)

    @test isapprox(AIC1, AIC2)
end

@testset "Wake Induced Velocity" begin

    # This test constructs two surfaces, one using surface panels,
    # and one using wake panels.  It then tests that the wake panel
    # implementations of `induced_velocity` yields identical results to the
    # surface panel implementations

    # generate the surface/wake geometry
    nc = 5
    ns = 10

    x = range(0, 1, length = nc+1)
    y = range(0, 2, length = ns+1)
    z = range(0, 3, length = ns+1)

    surface = Matrix{SurfacePanel{Float64}}(undef, nc, ns)
    wake = Matrix{WakePanel{Float64}}(undef, nc, ns)
    Γ = rand(nc, ns)

    for i = 1:nc, j = 1:ns
        rtl = [x[i], y[j], z[j]]
        rtr = [x[i], y[j+1], z[j+1]]
        rbl = [x[i+1], y[j], z[j]]
        rbr = [x[i+1], y[j+1], z[j+1]]
        rcp = (rtl + rtr + rbl + rbr)/4
        ncp = cross(rcp - rtr, rcp - rtl)
        core_size = 0.1
        chord = 0.0 # only used for unsteady simuulations

        surface[i,j] = SurfacePanel(rtl, rtr, rbl, rbr, rcp, ncp, core_size, chord)
        wake[i,j] = WakePanel(rtl, rtr, rbl, rbr, core_size, Γ[i, j])
    end

    # Test induced velocity calculation at an arbitrary point in space:
    rcp = [10, 11, 12]

    # no finite core model, no symmetry, no trailing vortices

    Vs = VortexLattice.induced_velocity(rcp, surface, Γ[:];
        finite_core = false,
        symmetric = false,
        trailing_vortices = false,
        xhat = [1, 0, 0])

    Vw = VortexLattice.induced_velocity(rcp, wake;
        finite_core = false,
        symmetric = false,
        trailing_vortices = false,
        xhat = [1, 0, 0])

    @test isapprox(Vs, Vw)

    # finite core model, symmetry, and trailing vortices

    Vs = VortexLattice.induced_velocity(rcp, surface, Γ[:];
        finite_core = true,
        symmetric = true,
        trailing_vortices = true,
        xhat = [1, 0, 0])

    Vw = VortexLattice.induced_velocity(rcp, wake;
        finite_core = true,
        symmetric = true,
        trailing_vortices = true,
        xhat = [1, 0, 0])

    @test isapprox(Vs, Vw)

    # Test induced velocity calculation at a trailing edge point:
    is = 3 # trailing edge index
    I = CartesianIndex(nc+1, is) # index on surface

    # no finite core model, no symmetry, no trailing vortices
    Vs = VortexLattice.induced_velocity(is, surface, Γ[:];
        finite_core = false,
        symmetric = false,
        trailing_vortices = false,
        xhat = [1, 0, 0])

    Vw = VortexLattice.induced_velocity(I, wake;
        finite_core = false,
        symmetric = false,
        trailing_vortices = false,
        xhat = [1, 0, 0])

    @test isapprox(Vs, Vw)

    # finite core model, symmetry, and trailing vortices
    Vs = VortexLattice.induced_velocity(is, surface, Γ[:];
        finite_core = true,
        symmetric = true,
        trailing_vortices = true,
        xhat = [1, 0, 0])

    Vw = VortexLattice.induced_velocity(I, wake;
        finite_core = true,
        symmetric = true,
        trailing_vortices = true,
        xhat = [1, 0, 0])

    @test isapprox(Vs, Vw)

end

@testset "Unsteady Vortex Lattice Method - Rectangular Wing" begin

    Uinf = 1.0
    AR = 4

    # reference parameters
    cref = 1.0
    bref = AR
    Sref = bref*cref
    rref = [0.0, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    # freestream parameters
    alpha = 5.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    # geometry
    xle = [0.0, 0.0]
    yle = [-bref/2, bref/2]
    zle = [0.0, 0.0]
    chord = [cref, cref]
    theta = [0.0, 0.0]
    phi = [0.0, 0.0]
    ns = 13
    nc = 4
    spacing_s = Uniform()
    spacing_c = Uniform()
    mirror = false
    symmetric = false

    # non-dimensional time
    t = range(0.0, 10.0, step=0.2)
    dt = [t[i+1]-t[i] for i = 1:length(t)-1]

    # create vortex rings
    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    _, _, surface = grid_to_surface_panels(grid; ratios=ratio)

    surfaces = [surface]

    # run analysis
    system, surface_history, property_history, wake_history = unsteady_analysis(surfaces, ref, fs, dt;
        symmetric=symmetric)

    # extract forces at each time step
    CF, CM = body_forces_history(system, surface_history, property_history; frame=Wind())
end

@testset "Unsteady Vortex Lattice Method - Wing + Tail" begin

    # Unsteady Wing and Tail

    # wing
    xle = [0.0, 0.2]
    yle = [0.0, 5.0]
    zle = [0.0, 1.0]
    chord = [1.0, 0.6]
    theta = [2.0*pi/180, 2.0*pi/180]
    phi = [0.0, 0.0]
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
    ns_v = 5
    nc_v = 1
    spacing_s_v = Uniform()
    spacing_c_v = Uniform()
    mirror_v = false

    # adjust chord lengths to match AVL (which uses chord length in the x-direction)
    chord = @. chord/cos(theta)
    chord_h = @. chord_h/cos(theta_h)
    chord_v = @. chord_v/cos(theta_v)

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

    # horseshoe vortices
    wgrid, wratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)
    _, _, wing = grid_to_surface_panels(wgrid; ratios=wratio)

    hgrid, hratio = wing_to_grid(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(hgrid, [4.0, 0.0, 0.0])
    _, _, htail = grid_to_surface_panels(hgrid; ratios=hratio)

    vgrid, vratio = wing_to_grid(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(vgrid, [4.0, 0.0, 0.0])
    _, _, vtail = grid_to_surface_panels(vgrid; ratios=vratio)

    grids = [wgrid, hgrid, vgrid]
    surfaces = [wing, htail, vtail]
    surface_id = [1, 2, 3]

    # t
    t = range(0.0, 10.0, step=0.2)
    dt = t[2:end] - t[1:end-1]

    system, surface_history, property_history, wake_history = unsteady_analysis(surfaces, ref, fs, dt; symmetric)

    # extract forces at each time step
    CF, CM = body_forces_history(system, surface_history, property_history; frame=Wind())
end

@testset "OpenVSP Geometry Import" begin
    Sref = 45.0
    cref = 2.5
    bref = 18.0
    rref = [0.625, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    fs = Freestream(Vinf, alpha, beta, Omega)

    comp = read_degengeom("samplewing.csv")
    grid, ratios, surface = import_vsp(comp[1]; mirror=true)

    symmetric = false
    grids = [grid]
    surfaces = [surface]
    ratios = [ratios]

    system = System(grids; ratios)

    steady_analysis!(system, ref, fs; symmetric=symmetric)
    CF, CM = body_forces(system; frame=Wind())

    CF_true = [2.41223539e-3, 0.0, 2.37009019e-1]
    CM_true = [0.0, -2.75091871e-1, 0.0]
    @test isapprox(CF, CF_true, atol=1e-5)
    @test isapprox(CM, CM_true, atol=1e-5)
end

#------- additional tests -------#

function _build_short_restart_case()
    grid, ratios = wing_to_grid([0.0, 0.0], [-1.0, 1.0], [0.0, 0.0],
        [1.0, 1.0], [0.0, 0.0], [0.0, 0.0], 2, 1;
        mirror=false, spacing_s=Uniform(), spacing_c=Uniform())

    system = System([grid]; nw=[2], ratios=[ratios])
    system.reference[] = Reference(2.0, 1.0, 2.0, [0.0, 0.0, 0.0], 10.0)
    system.freestream[] = Freestream(10.0, 0.0, 0.0, [0.0, 0.0, 0.0])

    for isurf in eachindex(system.surfaces)
        VortexLattice.update_surface_panels!(system.surfaces[isurf], system.grids[isurf];
            ratios=system.ratios[isurf], fcore=(c, Δs) -> system.core_size)
    end

    frames = ReferenceFrame(system;
        origin=SVector{3,Float64}(0.0, 0.0, 0.0),
        v=SVector{3,Float64}(0.0, 0.0, 0.0),
        ω_axis=SVector{3,Float64}(1.0, 0.0, 0.0),
        ω=0.0,
        R=SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0,
                                 0.0, 1.0, 0.0,
                                 0.0, 0.0, 1.0),
        name="vehicle",
        child_index=Int[],
        dependent_index=collect(1:length(system.surfaces)))
    maneuver!(frames, system, wake, t) = nothing
    Uinf(t) = SVector{3,Float64}(10.0, 0.0, 0.0)
    Ωinf(t) = SVector{3,Float64}(0.0, 0.0, 0.0)
    t_range = collect(0.0:0.05:0.10)

    return system, frames, maneuver!, Uinf, Ωinf, t_range
end

@testset "PanelParticleWake Restart (Short)" begin
    full_dir = mktempdir()
    restart_dir = mktempdir()

    system_full, frames_full, maneuver_full, Uinf_full, Ωinf_full, t_range = _build_short_restart_case()
    wake_full = simulate!(system_full, frames_full, maneuver_full, Uinf_full, t_range, Ωinf_full;
        wake_type=PanelParticleWake,
        nwakerows=2,
        max_particles=400,
        eta=0.3,
        method_trailing=OverlapPPS(1.3, 2),
        method_unsteady=OverlapPPS(1.3, 2),
        name="full_run",
        path=full_dir,
        verbose=false)

    system_part, frames_part, maneuver_part, Uinf_part, Ωinf_part, _ = _build_short_restart_case()
    simulate!(system_part, frames_part, maneuver_part, Uinf_part, t_range[1:2], Ωinf_part;
        wake_type=PanelParticleWake,
        nwakerows=2,
        max_particles=400,
        eta=0.3,
        method_trailing=OverlapPPS(1.3, 2),
        method_unsteady=OverlapPPS(1.3, 2),
        name="restart_run",
        path=restart_dir,
        verbose=false)

    system_restart, frames_restart, maneuver_restart, Uinf_restart, Ωinf_restart, _ = _build_short_restart_case()
    wake_restart = simulate!(system_restart, frames_restart, maneuver_restart, Uinf_restart, t_range, Ωinf_restart;
        wake_type=PanelParticleWake,
        nwakerows=2,
        max_particles=400,
        eta=0.3,
        method_trailing=OverlapPPS(1.3, 2),
        method_unsteady=OverlapPPS(1.3, 2),
        name="restart_run",
        path=restart_dir,
        restart_from=joinpath(restart_dir, "restart_run"),
        restart_idx=1,
        verbose=false)

    @test wake_full.nwake == wake_restart.nwake
    @test wake_full.overflowed[] == wake_restart.overflowed[]
    @test wake_full.pfield.np == wake_restart.pfield.np
    @test isapprox(system_full.Γ, system_restart.Γ; atol=0, rtol=0)
end

@testset "FluidDomainMonitor" begin
    system_fd, frames_fd, maneuver_fd, Uinf_fd, Ωinf_fd, t_range_fd = _build_short_restart_case()
    fd_dir = mktempdir()

    wake_fd = simulate!(system_fd, frames_fd, maneuver_fd, Uinf_fd, t_range_fd, Ωinf_fd;
        wake_type=PanelParticleWake,
        nwakerows=2,
        max_particles=400,
        eta=0.3,
        method_trailing=OverlapPPS(1.3, 2),
        method_unsteady=OverlapPPS(1.3, 2),
        name="fd_run",
        path=fd_dir,
        verbose=false)

    # Coarse 3×3×3 grid: well upstream/downstream and off to the side
    monitor_fd = FluidDomainMonitor(
        [-5.0, 0.0, 5.0], [-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0];
        vtk_interval=0, name="fd_test", path=fd_dir)

    evaluate_fluid_domain!(monitor_fd, system_fd, wake_fd)

    # All values must be finite
    @test all(isfinite(v[1]) && isfinite(v[2]) && isfinite(v[3])
              for v in monitor_fd.velocity)
    @test all(isfinite(v[1]) && isfinite(v[2]) && isfinite(v[3])
              for v in monitor_fd.vorticity)

    # Far upstream (x=-5, y=0, z=0) velocity x-component should be close to freestream (10 m/s)
    Vinf_ref = 10.0
    @test isapprox(monitor_fd.velocity[1, 2, 2][1], Vinf_ref; rtol=0.15)
end

@testset "Unsteady Viscous Coupling: frames_index / viscous frame regressions" begin
    # Regression coverage for two bugs found together (2026-08-20, BTV25/VPM-Validation):
    #  (1) simulate!()'s `frames_index` kwarg defaults to fill(-1, nsurf), which silently
    #      no-ops viscous!()/viscous_iterative_shed!()'s entire per-surface correction loop
    #      with no error/warning -- easy to forget, and produces plausible-looking-but-wrong
    #      results (a full pps/nwakerows sweep was run before this was caught).
    #  (2) simulate!() was internally feeding viscous!() the vehicle's KINEMATIC ReferenceFrame
    #      (built for Uinf(t)/Vh/Vv orientation) instead of a static axis-reference frame,
    #      corrupting viscous!()'s local lift-direction decomposition once frames_index was
    #      actually enabled (fixed via a dedicated `viscous_frames` built inside simulate!()).
    # Also covers nc>1 support added to viscous_iterative!()/viscous_iterative_shed!() (both
    # previously threw ArgumentError for nc>1).

    halfspan = 3.0
    root_chord = 1.0
    Vinf = 20.0
    ns = 8
    aoa_deg = 5.0

    # small synthetic cambered polar (cl0=0.3 at alpha=0, roughly thin-airfoil slope, mild
    # profile drag) -- deliberately different from the thin-airfoil cl=2*pi*alpha VLM result
    # so a no-op viscous correction is easy to detect.
    polar_alphas = collect(-20.0:2.0:20.0)
    polar_cls = 0.3 .+ 2π .* deg2rad.(polar_alphas)
    polar_cds = 0.01 .+ 0.001 .* polar_alphas.^2
    polar = VortexLattice.Polar(polar_alphas, polar_cls, polar_cds)

    function build_wing(nc)
        xle = [0.0, 0.0]; yle = [0.0, halfspan]; zle = zeros(2)
        chord = [root_chord, root_chord]; theta = zeros(2); phi = zeros(2)
        Sref = halfspan * 2 * root_chord
        cref = root_chord
        bref = halfspan * 2
        ref = Reference(Sref, cref, bref, [0.0, 0.0, 0.0], Vinf)
        grid, ratios = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
            spacing_s=Uniform(), spacing_c=Uniform(), mirror=true)
        grids = [grid]
        ratios_v = [ratios]
        ns_total = size(grid, 3) - 1
        polars = [fill(polar, ns_total)]
        return (; ref, grids, ratios_v, cref, ns_total, polars)
    end

    function body_CL(system)
        CF, _ = body_forces(system; frame=Wind())
        return CF[3]
    end

    fs = Freestream(Vinf, deg2rad(aoa_deg), 0.0, [0.0, 0.0, 0.0])

    function run_unsteady(wing; polars_arg=nothing, viscous_iterative_shed=false,
            frames_index=nothing, nwakerows=1, pps=2, nchords=4.0, steps_per_chord=3.0)
        grids_case = deepcopy(wing.grids)
        sys = System(grids_case; ratios=wing.ratios_v, nw=fill(nwakerows, length(grids_case)))
        sys.reference[] = wing.ref
        wakes = [Matrix{WakePanel{Float64}}(undef, nwakerows, size(sys.surfaces[i], 2)) for i in 1:length(sys.surfaces)]
        steady_analysis!(sys, sys.reference[], fs; symmetric=false, wakes, nwake=fill(nwakerows, length(sys.surfaces)))
        frames = ReferenceFrame(sys;
            origin=SVector{3}(0.0, 0.0, -10.0), v=SVector{3}(0.0, 0.0, 0.0),
            ω_axis=SVector{3}(0.0, 1.0, 0.0), ω=0.0,
            R=SMatrix{3,3}(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0),
            name="vehicle", child_index=Int[], dependent_index=collect(1:length(sys.surfaces)))
        alpha_rad = deg2rad(aoa_deg)
        Uinf(t) = SVector{3,Float64}(Vinf * cos(alpha_rad), 0.0, Vinf * sin(alpha_rad))
        Ωinf(t) = SVector{3,Float64}(0.0, 0.0, 0.0)
        dt = wing.cref / (Vinf * steps_per_chord)
        nsteps = round(Int, nchords * steps_per_chord)
        t_range = range(start=0.0, stop=nsteps * dt, length=nsteps + 1)
        constant_maneuver!(fr, s, wk, t) = nothing
        max_particles = (wing.ns_total + 1) * pps * nsteps + 2000
        kwargs = (; wake_type=PanelParticleWake, nwakerows=nwakerows, max_particles=max_particles,
            method_trailing=OverlapPPS(1.3, pps), method_unsteady=NoShed(),
            derivatives=false, path=nothing, verbose=false)
        if !isnothing(polars_arg)
            kwargs = (; kwargs..., polars=polars_arg, viscous_iterative_shed=viscous_iterative_shed)
        end
        if !isnothing(frames_index)
            kwargs = (; kwargs..., frames_index=frames_index)
        end
        simulate!(sys, frames, constant_maneuver!, Uinf, t_range, Ωinf; kwargs...)
        return body_CL(sys)
    end

    wing1 = build_wing(1)

    CL_inviscid_unsteady = run_unsteady(wing1)

    # (1) forgetting frames_index with polars set must equal the pure-inviscid result -- this
    # is the documented (if unfortunate) default behavior, tested explicitly so a future change
    # to this default doesn't go unnoticed either way.
    CL_forgot_frames_index = run_unsteady(wing1; polars_arg=wing1.polars, viscous_iterative_shed=true)
    @test isapprox(CL_forgot_frames_index, CL_inviscid_unsteady; atol=1e-3)

    # (2) with frames_index correctly set, both viscous coupling methods must produce a REAL,
    # substantial change from the inviscid result (guards against bug (1) reappearing) ...
    fidx = fill(1, length(wing1.grids))
    CL_plain = run_unsteady(wing1; polars_arg=wing1.polars, viscous_iterative_shed=false, frames_index=fidx)
    CL_iterative = run_unsteady(wing1; polars_arg=wing1.polars, viscous_iterative_shed=true, frames_index=fidx)
    @test abs(CL_plain - CL_inviscid_unsteady) > 0.05
    @test abs(CL_iterative - CL_inviscid_unsteady) > 0.05

    # ... and must be reasonably close to the steady single-pass viscous!() target (guards
    # against bug (2), the ReferenceFrame mismatch, reappearing -- that bug alone was capable
    # of driving CL negative).
    system_steady_visc = System(deepcopy(wing1.grids); ratios=wing1.ratios_v)
    steady_analysis!(system_steady_visc, wing1.ref, fs; symmetric=false, derivatives=false)
    frames_static = ReferenceFrame(system_steady_visc;
        origin=SVector{3}(0.0, 0.0, 0.0), v=SVector{3}(0.0, 0.0, 0.0),
        ω_axis=SVector{3}(1.0, 0.0, 0.0), ω=0.0,
        R=SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
        name="vehicle", child_index=Int[], dependent_index=collect(1:length(system_steady_visc.surfaces)))
    VortexLattice.viscous!(system_steady_visc.properties, system_steady_visc.Γ, system_steady_visc.surfaces,
        wing1.grids, frames_static, fidx, wing1.polars, wing1.ref, 0.0)
    CL_steady_target = body_CL(system_steady_visc)

    @test CL_steady_target > 0.3 # sanity: the synthetic polar's camber should raise CL well above the thin-airfoil inviscid value
    @test isapprox(CL_plain, CL_steady_target; rtol=0.3)
    @test isapprox(CL_iterative, CL_steady_target; rtol=0.3)

    # (3) nc>1 support: viscous_iterative!() (steady) must not throw for nc>1, and should be
    # close to mesh-independent in nc, matching the nc-independence already established for the
    # custom Cocco NL-VL scheme on this same style of case.
    function steady_iterative_CL(wing)
        system = System(deepcopy(wing.grids); ratios=wing.ratios_v)
        steady_analysis!(system, wing.ref, fs; symmetric=false, derivatives=false)
        fidx_local = fill(1, length(system.surfaces))
        VortexLattice.viscous_iterative!(system.properties, system.Γ, system.surfaces, system.wakes, wing.grids,
            fidx_local, wing.polars, wing.ref, fs;
            symmetric=system.symmetric, nwake=system.nwake, surface_id=system.surface_id,
            wake_finite_core=system.wake_finite_core, wake_shedding_locations=nothing,
            trailing_vortices=system.trailing_vortices, xhat=system.xhat[])
        return body_CL(system)
    end

    CL_nc1 = steady_iterative_CL(build_wing(1))
    CL_nc2 = steady_iterative_CL(build_wing(2))
    CL_nc4 = steady_iterative_CL(build_wing(4))
    @test isapprox(CL_nc2, CL_nc1; rtol=0.05)
    @test isapprox(CL_nc4, CL_nc1; rtol=0.05)

    # nc>1 must also not throw for the unsteady viscous_iterative_shed! path.
    wing2 = build_wing(2)
    CL_iterative_nc2 = run_unsteady(wing2; polars_arg=wing2.polars, viscous_iterative_shed=true,
        frames_index=fill(1, length(wing2.grids)))
    @test isfinite(CL_iterative_nc2)
end

include("fmm_test.jl")
