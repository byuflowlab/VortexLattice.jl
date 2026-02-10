using VortexLattice
using VortexLattice.StaticArrays
using PythonPlot
using LinearAlgebra

function get_CL(system, ref, alpha; symmetric = false) # radians

    # freestream parameters
    # alpha = 1.0*pi/180
    beta = 0.0
    Omega = [0.0; 0.0; 0.0]
    Vinf = ref.V
    fs = Freestream(Vinf, alpha, beta, Omega)

    # perform steady state analysis
    steady_analysis!(system, ref, fs; symmetric=symmetric, derivatives=false)

    # retrieve near-field forces
    CF, CM = body_forces(system; frame=Wind())

    # perform far-field analysis
    CDiff = far_field_drag(system)

    CD, CY, CL = CF
    Cl, Cm, Cn = CM

    return CL

end

function get_system(AR; ns=24, nc=6)
    # geometry (right half of the wing)
    c = 2.0
    b = AR * c
    xle = [0.0, 0.0]
    yle = [-b/2, b/2]
    zle = [0.0, 0.0]
    chord = [c, c]
    theta = [0*pi/180, 0*pi/180]
    phi = [0.0, 0.0]
    fc = fill((xc) -> 0, 2) # camberline function for each section

    # discretization parameters
    spacing_s = Uniform()
    spacing_c = Uniform()

    # reference parameters
    Sref = b * c
    cref = c
    bref = b
    rref = [0.25 * c, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)
    
    # construct grid
    grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
        fc = fc, spacing_s=spacing_s, spacing_c=spacing_c)

    # create vector containing all grids
    grids=[grid]
    ratios=[ratio]

    # Construct the system
    system = System(grids; ratios)

    return system, ref
end

function plot_stuff(ARs = [10, 20, 40, 80, 160, 320, 640])
    fig = figure("polar")
    fig.clear()
    ax = fig.add_subplot(111, xlabel=L"\alpha", ylabel=L"c_l")
    alphas = range(0.0, 10.0, length=11) .* pi/180
    cls_list = zeros(length(alphas), length(ARs))
    for (i,AR) in enumerate(ARs)
        @show AR
        ns = 2 * AR
        nc = 6
        system, ref = get_system(AR; ns, nc)
        cls = [get_CL(system, ref, alpha) for alpha in alphas]
        cls_list[:, i] .= cls

        # plot results
        ax.plot(alphas, cls, label="AR=$AR")
    end
    ax.legend()

    return cls_list
end

function get_CL_circ(system::System{TF}, alpha, vinf_vec) where TF
    println("\n==== BEGIN CL CALCULATION ====\n")

    # solve system with unit reference
    Sref = 2.0
    cref = 1.0
    bref = 1.0
    rref = [0.25 * 2, 0.0, 0.0]
    Vinf = 1.0
    ref = Reference(Sref, cref, bref, rref, Vinf)
    l = get_CL(system, ref, alpha)
    @show l

    # choose which spanwise section to analyze
    nc, ns = size(system.surfaces[1])
    surface = system.surfaces[1]
    grid = system.grids[1]
    props = system.properties[1]

    # index for circulation strengths
    
    # extract containers corresponding to this surface
    surface = system.surfaces[1]
    props = system.properties[1]
    grid = system.grids[1]
    
    # get rotation matrix from global frame to this frame
    Rp = SMatrix{3,3,Float64,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    R = transpose(Rp)
    
    # loop over spanwise sections in this surface
    cls_new = zeros(size(surface, 2))
    # for j in axes(surface, 2)
    j = ns >> 2
    iΓ = j
    # j = 1

        # get the circulation per unit length of this section
        γ = zero(TF)

        # get induced velocity
        v_induced = zero(SVector{3,TF})

        # get net aerodynamic force on this section
        cf = zero(SVector{3,TF})

        # loop over chordwise panels in this section
        for i in axes(surface, 1)

            # project bound vortex onto y axis of this frame
            ds = R * (surface[i,j].rtr - surface[i,j].rtl)
            dy = abs(ds[2]) / norm(ds)
            @show dy

            # accumulate circulation contribution from this bound vortex
            γ += system.Γ[iΓ] * dy
            iΓ += 1

            # accumulate induced velocity contribution at this bound vortex
            v_induced += props[i,j].velocity # TODO: what if system.reference[].v != 1.0?

            # accumulate aerodynamic force contribution from this bound vortex
            cf += props[i,j].cfb
        end

        # get chord length
        le = SVector(grid[1,1,j], grid[2,1,j], grid[3,1,j])
        te = SVector(grid[1,end,j], grid[2,end,j], grid[3,end,j])
        c = norm(le - te)

        # average of chord lengths of either side
        le = SVector(grid[1,1,j+1], grid[2,1,j+1], grid[3,1,j+1])
        te = SVector(grid[1,end,j+1], grid[2,end,j+1], grid[3,end,j+1])
        c += norm(le - te)
        c *= 0.5
        @show c
        
        # project aerodynamic force into xz plane
        @show cf
        cf = R * cf # rotate into this frame
        @show cf
        
        # get 2-D lift magnitude of this section
        l_2d_norm = sqrt(cf[1]*cf[1] + cf[3]*cf[3])

        # force per length
        l_2d_norm /= norm(surface[end,j].rtr - surface[1,j].rtl)
        @show l_2d_norm
        
        # calculate effective cl
        cl_vlm = 2 * VortexLattice.RHO * γ * γ / (l_2d_norm * c) * sign(γ)
        # cls_new[j] = cl_vlm

        @show cl_vlm

        # check velocity
        this_v = l_2d_norm / (VortexLattice.RHO * γ)
        @show this_v, norm(v_induced) # checks out

    # return cls_new
    return cl_vlm
end

# plot_stuff([5])
# cls_list = plot_stuff()

# slope = (cls_list[end, end] - cls_list[1, end]) / (10 * pi/180)
# println("Relative error in lift slope: ", (slope - 2*pi) / (2*pi))

#--- verify my cl prediction ---#
AR = 100
ns = 3 * AR
nc = 1
system, ref = get_system(AR; ns, nc)
@show ref.V ref.S ref.c ref.b
alpha = 10.0 * pi/180
vinf_vec = SVector(ref.V * cos(alpha), 0.0, ref.V * sin(alpha))

println("\nBefore New Stuff\n")
CL_vlm = get_CL(system, ref, alpha)
L_vlm = CL_vlm * 0.5 * VortexLattice.RHO * ref.V^2 * ref.S
@show L_vlm
@show system.properties[1][1,ns>>1].cfb

CL_circulation = get_CL_circ(system, alpha, vinf_vec)

@show CL_vlm
@show CL_circulation
@show 2 * pi * alpha
