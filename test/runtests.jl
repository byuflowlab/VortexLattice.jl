using Test
using VLMFlow
using LinearAlgebra

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

@testset "AVL - Run 1" begin

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
    ref = Reference(Sref, cref, bref, rref)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    # adjust chord length so x-chord length matches AVL
    chord = @. chord/cos(theta)

    # horseshoe vortices with symmetry
    mirror = false
    symmetric = true

    wing = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.24324, atol=1e-3)
    @test isapprox(CD, 0.00243, atol=1e-5)
    @test isapprox(CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02252, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # horseshoe vortices with mirrored geometry
    mirror = true
    symmetric = false

    wing = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s, spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.24324, atol=1e-3)
    @test isapprox(CD, 0.00243, atol=1e-5)
    @test isapprox(CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02252, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings with symmetry
    mirror = false
    symmetric = true

    wing = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

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

    wing = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

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

@testset "AVL - Run 2" begin

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
    ref = Reference(Sref, cref, bref, rref)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    # adjust chord length so x-chord length matches AVL
    chord = @. chord/cos(theta)

    # horseshoe vortices
    wing = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.23744, atol=1e-3)
    @test isapprox(CD, 0.00254, atol=1e-5)
    @test isapprox(CDiff, 0.00243, atol=1e-5)
    @test isapprox(Cm, -0.02165, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings
    wing = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

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

@testset "AVL - Run 3" begin

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
    ref = Reference(Sref, cref, bref, rref)

    alpha = 8.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    # adjust chord length so x-chord length matches AVL
    chord = @. chord/cos(theta)

    # horseshoe vortices
    wing = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.80348, atol=1e-5)
    @test isapprox(CD, 0.02651, atol=1e-5)
    @test isapprox(CDiff, 0.02696, atol=1e-5)
    @test isapprox(Cm, -0.07399, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings
    wing = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

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

@testset "AVL - Run 4" begin

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
    ref = Reference(Sref, cref, bref, rref)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    # adjust chord length so x-chord length matches AVL
    chord = @. chord/cos(theta)

    # also get normal vector as AVL defines it
    ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

    # horseshoe vortices
    wing = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    for (ip, p) in enumerate(wing)
        # check that our normal vector is approximately the same as AVL's
        @test isapprox(p.ncp, ncp, rtol=0.01)
        # replace our normal vector with AVL's normal vector for this test
        wing[ip] = Horseshoe(p.rl, p.rc, p.rr, p.rcp, ncp, p.xl_te, p.xc_te, p.xr_te, p.core_size)
    end

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.24787, atol=1e-3)
    @test isapprox(CD, 0.00246, atol=1e-5)
    @test isapprox(CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02395, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings
    wing = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    for (ip, p) in enumerate(wing)
        # check that our normal vector is approximately the same as AVL's
        @test isapprox(p.ncp, ncp, rtol=0.01)
        # replace our normal vector with AVL's normal vector for this test
        wing[ip] = Ring(p.rtl, p.rtc, p.rtr, p.rbl, p.rbc, p.rbr, p.rcp, ncp, p.core_size)
    end

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.24787, atol=1e-3)
    @test isapprox(CD, 0.00246, atol=1e-4)
    @test isapprox(CDiff, 0.00245, atol=1e-5)
    @test isapprox(Cm, -0.02395, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)
end

@testset "AVL - Run 5" begin

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
    ref = Reference(Sref, cref, bref, rref)

    alpha = 20.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

    # horseshoe vortices
    wing = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    for (ip, p) in enumerate(wing)
        # check that our normal vector is approximately the same as AVL's
        @test isapprox(p.ncp, ncp, rtol=0.01)
        # replace our normal vector with AVL's normal vector for this test
        wing[ip] = Horseshoe(p.rl, p.rc, p.rr, p.rcp, ncp, p.xl_te, p.xc_te, p.xr_te, p.core_size)
    end

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 1.70982, rtol=0.01)
    @test isapprox(CD, 0.12904, rtol=0.01)
    @test isapprox(CDiff, 0.11502, rtol=0.01)
    @test isapprox(Cm, -0.45606, rtol=0.01)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings, untwisted geometry
    wing = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    for (ip, p) in enumerate(wing)
        # check that our normal vector is approximately the same as AVL's
        @test isapprox(p.ncp, ncp, rtol=0.01)
        # replace our normal vector with AVL's normal vector for this test
        wing[ip] = Ring(p.rtl, p.rtc, p.rtr, p.rbl, p.rbc, p.rbr, p.rcp, ncp, p.core_size)
    end

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

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

@testset "AVL - Run 6" begin

    # Wing and Tail

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
    ref = Reference(Sref, cref, bref, rref)

    alpha = 5.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    symmetric = true

    ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

    # horseshoe vortices
    wing = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    for (ip, p) in enumerate(wing)
        # check that our normal vector is approximately the same as AVL's
        @test isapprox(p.ncp, ncp, rtol=0.01)
        # replace our normal vector with AVL's normal vector for this test
        wing[ip] = Horseshoe(p.rl, p.rc, p.rr, p.rcp, ncp, p.xl_te, p.xc_te, p.xr_te, p.core_size)
    end

    htail = wing_to_horseshoe_vortices(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(htail, [4.0, 0.0, 0.0])

    vtail = wing_to_horseshoe_vortices(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    translate!(vtail, [4.0, 0.0, 0.0])

    surfaces = [wing, htail, vtail]

    AIC = influence_coefficients(surfaces, symmetric)
    b = normal_velocity(surfaces, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(surfaces, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(surfaces, ref, fs, symmetric, Γ)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.60408, atol=1e-4)
    @test isapprox(CD, 0.01058, atol=1e-5)
    @test isapprox(CDiff, 0.010378, atol=1e-4)
    @test isapprox(Cm, -0.02778, atol=1e-4)
    @test isapprox(CY, 0.0, atol=ztol)
    @test isapprox(Cl, 0.0, atol=ztol)
    @test isapprox(Cn, 0.0, atol=ztol)

    # vortex rings
    wing = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    for (ip, p) in enumerate(wing)
        # check that our normal vector is approximately the same as AVL's
        @test isapprox(p.ncp, ncp, rtol=0.01)
        # replace our normal vector with AVL's normal vector for this test
        wing[ip] = Ring(p.rtl, p.rtc, p.rtr, p.rbl, p.rbc, p.rbr, p.rcp, ncp, p.core_size)
    end

    htail = wing_to_vortex_rings(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
        mirror=mirror_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
    translate!(htail, [4.0, 0.0, 0.0])

    vtail = wing_to_vortex_rings(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
        mirror=mirror_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
    translate!(vtail, [4.0, 0.0, 0.0])

    surfaces = [wing, htail, vtail]

    AIC = influence_coefficients(surfaces, symmetric)
    b = normal_velocity(surfaces, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(surfaces, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(surfaces, ref, fs, symmetric, Γ)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.60408, atol=1e-3)
    @test isapprox(CD, 0.01058, atol=1e-4)
    @test isapprox(CDiff, 0.010378, atol=1e-4)
    @test isapprox(Cm, -0.02778, atol=2e-3)
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
    ref = Reference(Sref, cref, bref, rref)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    symmetric = true

    # horseshoe vortices
    wing = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.24454, atol=1e-3)
    @test isapprox(CD, 0.00247, atol=1e-5)
    @test isapprox(CDiff, 0.00248, atol=1e-5)
    @test isapprox(Cm, -0.02091, atol=1e-4)
    @test isapprox(CY, 0.0, atol=1e-16)
    @test isapprox(Cl, 0.0, atol=1e-16)
    @test isapprox(Cn, 0.0, atol=1e-16)

    # vortex rings
    wing = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

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

@testset "AVL - Run 8" begin

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
    ref = Reference(Sref, cref, bref, rref)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    symmetric = true

    # horseshoe vortices
    wing = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.23879, atol=1e-3)
    @test isapprox(CD, 0.00249, atol=1e-5)
    @test isapprox(CDiff, 0.0024626, atol=1e-5)
    @test isapprox(Cm, -0.01995, atol=1e-4)
    @test isapprox(CY, 0.0, atol=1e-16)
    @test isapprox(Cl, 0.0, atol=1e-16)
    @test isapprox(Cn, 0.0, atol=1e-16)

    # vortex rings
    wing = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

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

@testset "AVL - Run 9" begin

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
    ref = Reference(Sref, cref, bref, rref)

    alpha = 1.0*pi/180
    beta = 15.0*pi/180
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    ncp = avl_normal_vector([xle[2]-xle[1], yle[2]-yle[1], zle[2]-zle[1]], 2.0*pi/180)

    # horseshoe vortices
    mirror = true
    symmetric = false

    wing = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.22695, atol=1e-3)
    @test isapprox(CD, 0.00227, atol=1e-5)
    @test isapprox(CDiff, 0.0022852, atol=1e-5)
    @test isapprox(Cm, -0.02101, atol=1e-4)
    @test isapprox(CY, 0.0, atol=1e-5)
    @test isapprox(Cl, 0.00644, atol=1e-4)
    @test isapprox(Cn, -0.00012, atol=1e-4)

    # vortex rings
    mirror = true
    symmetric = false

    wing = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)
    b = normal_velocity(wing, ref, fs)
    Γ = circulation(AIC, b)
    CF, CM, props = near_field_forces(wing, ref, fs, symmetric, Γ; frame=Stability())
    CDiff = far_field_drag(wing, ref, fs, symmetric, Γ)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    @test isapprox(CL, 0.22695, atol=1e-3)
    @test isapprox(CD, 0.00227, atol=1e-5)
    @test isapprox(CDiff, 0.0022852, atol=1e-5)
    @test isapprox(Cm, -0.02101, atol=1e-4)
    @test isapprox(CY, 0.0, atol=1e-5)
    @test isapprox(Cl, 0.00644, atol=1e-4)
    @test isapprox(Cn, -0.00012, atol=2e-4)

end

@testset "AVL - Run 10" begin

    # Run 9: Simple Wing Stability Derivatives

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
    ref = Reference(Sref, cref, bref, rref)

    alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    vother = nothing
    fs = Freestream(alpha, beta, Omega, vother)

    # horseshoe vortices
    wing = wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)

    dCF, dCM = stability_derivatives(wing, ref, fs, symmetric, AIC)

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
    @test isapprox(CYb, -0.000007, atol=1e-4)
    @test isapprox(Cla, 0.0, atol=ztol)
    @test isapprox(Clb, 0.025749, rtol=0.015)
    @test isapprox(Cma, -0.429247, rtol=0.01)
    @test isapprox(Cmb, 0.0, atol=ztol)
    @test isapprox(Cna, 0.0, atol=ztol)
    @test isapprox(Cnb, -0.000466, atol=1e-4)
    @test isapprox(Clp, -0.518725, rtol=0.01)
    @test isapprox(Clq, 0.0, atol=ztol)
    @test isapprox(Clr, 0.064243, rtol=0.01)
    @test isapprox(Cmp, 0.0, atol=ztol)
    @test isapprox(Cmq, -0.517094, rtol=0.01)
    @test isapprox(Cmr, 0.0, atol=ztol)
    @test isapprox(Cnp, -0.019846, rtol=0.01)
    @test isapprox(Cnq, 0.0, atol=ztol)
    @test isapprox(Cnr, -0.000898, rtol=0.01)

    # vortex rings
    wing = wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, ns, nc;
        mirror=mirror, spacing_s=spacing_s, spacing_c=spacing_c)

    AIC = influence_coefficients(wing, symmetric)

    dCF, dCM = stability_derivatives(wing, ref, fs, symmetric, AIC)

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
    @test isapprox(Clb, 0.025749, rtol=0.015)
    @test isapprox(Cma, -0.429247, rtol=0.01)
    @test isapprox(Cmb, 0.0, atol=ztol)
    @test isapprox(Cna, 0.0, atol=ztol)
    # @test isapprox(Cnb, -0.000466, atol=1e-4)
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
# rref = [0.0, 0.0, 0.0]
# ref = VLM.Reference(Sref, cref, bref, rref)

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
