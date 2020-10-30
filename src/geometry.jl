"""
    AbstractSpacing

Spacing discretization scheme supertype
"""
abstract type AbstractSpacing end

"""
    Uniform()

Uniform discretization scheme
"""
struct Uniform <: AbstractSpacing end

"""
    Sine()

Sine-spaced discretization scheme.
"""
struct Sine <: AbstractSpacing end

"""
    Cosine()

Cosine-spaced discretization scheme
"""
struct Cosine <: AbstractSpacing end

"""
   linearinterp(eta, rstart, rend)

Linearly interpolate between `rstart` and `rend` where eta is the fraction
between 0 (rstart) and 1 (rend)
"""
linearinterp(eta, rstart, rend) = (1-eta)*rstart + eta*rend

"""
    spanwise_spacing(n, spacing::AbstractSpacing)

Distribute `n` panel endpoints and `n-1` panel midpoints on the interval between
0 and 1 according to the discretization strategy in `spacing`.
"""
spanwise_spacing

# uniform
function spanwise_spacing(n, ::Uniform)

    eta = range(0, 1.0, length=n)
    eta_mid = linearinterp(0.5, eta[1:end-1], eta[2:end])

    return eta, eta_mid
end

# sine
function spanwise_spacing(n, ::Sine)

    theta = range(0, pi/2, length=n)
    eta = sin.(theta)

    # note that control points are also placed with cosine spacing as this improves accuracy tremendously
    theta_mid = linearinterp(0.5, theta[1:end-1], theta[2:end])
    eta_mid = sin.(theta_mid)

    return eta, eta_mid
end

# cosine
function spanwise_spacing(n, ::Cosine)
    theta = range(0, pi, length=n)
    eta = (1.0 .- cos.(theta))/2.0

    # note that control points are also placed with cosine spacing as this improves accuracy tremendously
    theta_mid = linearinterp(0.5, theta[1:end-1], theta[2:end])
    eta_mid = (1.0 .- cos.(theta_mid))/2.0

    return eta, eta_mid
end

"""
    spanwise_spacing(n, spacing::AbstractSpacing)

Distribute `n-1` vortex and `n-1` control point chordwise locations on the
interval between 0 and 1 according to the discretization strategy in `spacing`.
"""
chordwise_spacing

# uniform
function chordwise_spacing(n, ::Uniform)

    eta = range(0, 1.0, length=n)

    eta_qtr = linearinterp(0.25, eta[1:end-1], eta[2:end])
    eta_thrqtr = linearinterp(0.75, eta[1:end-1], eta[2:end])

    return eta_qtr, eta_thrqtr
end

# sine
function chordwise_spacing(n, ::Sine)

    theta = range(0, pi/2, length=n)
    eta = sin.(theta)

    eta_qtr = linearinterp(0.25, eta[1:end-1], eta[2:end])
    eta_thrqtr = linearinterp(0.75, eta[1:end-1], eta[2:end])

    return eta, eta_mid
end

# cosine
function chordwise_spacing(n, ::Cosine)

    theta = range(0, pi, length=n)
    eta = (1.0 .- cos.(theta))/2.0

    eta_qtr = linearinterp(0.25, eta[1:end-1], eta[2:end])
    eta_thrqtr = linearinterp(0.75, eta[1:end-1], eta[2:end])

    return eta_qtr, eta_thrqtr
end

"""
    interpolate_grid(xyz, eta, interp; xdir=0, ydir=1)

Interpolates the grid `xyz` along direction `dir`

# Arguments
 - `xyz`: Grid of size (3, ni, nj)
 - `eta`: New (normalized) coordinates in direction `dir` (0 <= eta <= 1)
 - `interp`: Interpolation method of form `ypt = f(x,y,xpt)`
 - `xdir`: Independent variable direction, defaults to arc length
 - `ydir`: Dependent variable direction `xyz` (`i=1`, `j=2`)
"""
function interpolate_grid(xyz, eta, interp; xdir=0, ydir=1)

    ydim = ydir + 1

    ni = size(xyz, ydim)

    xyz = mapslices(xyz; dims=[1,ydim]) do xyz_i

        # get distance between each grid location
        if xdir == 0
            ds = [norm(xyz_i[:,k] - xyz_i[:,max(1,k-1)]) for k = 1:ni]
        else
            ds = [xyz_i[xdir,k] - xyz_i[xdir,max(1,k-1)] for k = 1:ni]
        end

        # create interpolation vector
        t = cumsum(ds)

        # normalize interpolation vector
        t /= t[end]

        # interpolate x, y, and z
        x = interp(t, xyz_i[1,:], eta)
        y = interp(t, xyz_i[2,:], eta)
        z = interp(t, xyz_i[3,:], eta)

        vcat(x',y',z')
    end

    return xyz
end

"""
    grid_to_horseshoe_vortices(xyz; mirror=false, fcore=(c, Δs)->max(c/4, Δs))

Construct a set of horseshoe vortex panels given a potentially curved lifting surface
defined by a grid with dimensions (3, i, j) where `i` corresponds to the chordwise
direction (ordered from leading edge to trailing edge) and `j` corresponds to the spanwise
direction (ordered from left to right).  The bound vortex will be placed at the
1/4 chord and the control point will be placed at the 3/4 chord of each panel.

In order to be consistent with vortex lattice theory using horseshoe vortices,
the y and z-coordinates of each resulting chordwise strip of panels are set
to the same y and z-coordinates as found at the leading edge of the chordwise
strip.  The normal vectors used to satisfy the no-flow boundary conditions,
however, are set using the original geometry.

The optional argument `mirror` may be used to mirror the geometry across the
y-axis.

The argument `fcore` defines a function for setting the finite core size based
on the chord length in the x-direction and/or the panel width in the y/z directions.
The default function is the same as that used by AVL.
"""
function grid_to_horseshoe_vortices(xyz; mirror=false, fcore=(c, Δs)->max(c/4, Δs/2))

    TF = eltype(xyz)

    nc = size(grid, 2)-1 # number of chordwise panels
    ns = size(grid, 3)-1 # number of spanwise panels

    # check which side we're working with
    right_side = sum(xyz[2,:,:]) > 0

    # initialize output
    panels = Matrix{Horseshoe{TF}}(undef, nc, (1+mirror)*ns)

    # populate each panel
    for j = 1:ns

        # get y and z-coordinates for this chordwise strip
        rly = grid[2,1,j]
        rlz = grid[3,1,j]

        rry = grid[2,1,j+1]
        rrz = grid[3,1,j+1]

        rcy = (yl + yr)/2
        rcz = (zl + zr)/2

        # also get chord length for setting finite core size
        cl = grid[1,end,j] - grid[1,1,j]
        cr = grid[1,end,j+1] - grid[1,1,j+1]
        c = (cl + cr)/2

        for i = 1:nc

            # extract grid corners
            r1 = SVector(grid[1,i,j], grid[2,i,j], grid[3,i,j]) # top left
            r2 = SVector(grid[1,i,j+1], grid[2,i,j+1], grid[3,i,j+1]) # top right
            r3 = SVector(grid[1,i+1,j], grid[2,i+1,j], grid[3,i+1,j]) # bottom left
            r4 = SVector(grid[1,i+1,j+1], grid[2,i+1,j+1], grid[3,i+1,j+1]) # bottom right

            # left side of bound vortex
            rl = linearinterp(0.25, r1, r3)

            # right side of bound vortex
            rr = linearinterp(0.25, r2, r4)

            # center of bound vortex
            rc = (rl + rr)/2

            # control point
            rtop = linearinterp(0.5, r1, r2)
            rbot = linearinterp(0.5, r3, r4)
            rcp = linearinterp(0.75, rtop, rbot)

            # surface normal
            ncp = cross(rcp - rr, rcp - rl)
            ncp /= norm(ncp)

            # set y and z-coordinates to correspond to the leading edge
            rl = SVector(rl[1], rly, rlz)
            rr = SVector(rr[1], rry, rrz)
            rc = SVector(rc[1], rcy, rcz)
            rcp = SVector(rcp[1], rcy, rcz)

            # set x-distance to trailing edge
            xl_te = grid[1,end,j] - rl[1]
            xr_te = grid[1,end,j+1] - rr[1]
            xc_te = (xl_te + xr_te)/2

            # set finite core size
            Δs = sqrt((rr[2]-rl[2])^2 + (rr[3]-rl[3])^2)
            core_size = fcore(c, Δs)

            ip = i
            jp = mirror*right_side*ns + j
            panels[ip,jp] = Horseshoe{TF}(rl, rc, rr, rcp, ncp, xl_te, xc_te, xr_te, core_size)
        end
    end

    # other side
    if mirror
        for j = 1:ns
            for i = 1:nc
                ip = i
                jp1 = right_side*ns + j
                jp2 = 2*ns - right_side*ns - j + 1
                panels[ip,jp2] = reflect(panels[ip,jp1])
            end
        end
    end

    return panels
end

"""
    grid_to_vortex_rings(xyz; mirror=false, fcore=(c, Δs)->max(c/4, Δs/2))

Construct a set of vortex ring panels given a potentially curved lifting surface
defined by a grid with dimensions (3, i, j) where `i` corresponds to the chordwise
direction (ordered from leading edge to trailing edge) and `j` corresponds to the spanwise
direction (ordered from left to right).  The leading edge of each ring vortex
will be placed at the 1/4 chord and the control point will be placed at the 3/4
chord of each panel.

The optional argument `mirror` may be used to mirror the geometry across the
y-axis.

The argument `fcore` defines a function for setting the finite core size based
on the chord length in the x-direction and/or the panel width in the y/z directions.
The default function is the same as that used by AVL.
"""
function grid_to_vortex_rings(xyz; mirror=false, fcore=(c, Δs)->max(c/4, Δs/2))

    TF = eltype(xyz)

    nc = size(grid, 2)-1 # number of chordwise panels
    ns = size(grid, 3)-1 # number of spanwise panels
    N = (1+mirror)*nc*ns # total number of panels

    # check which side we're working with
    right_side = sum(xyz[2,:,:]) > 0

    # initialize output
    panels = Matrix{Ring{TF}}(undef, nc, (1+mirror)*ns)

    # populate each panel
    for j = 1:ns

        # get leading edge panel corners
        r1n = SVector(grid[1,1,j], grid[2,1,j], grid[3,1,j]) # top left
        r2n = SVector(grid[1,1,j+1], grid[2,1,j+1], grid[3,1,j+1]) # top right
        r3n = SVector(grid[1,2,j], grid[2,2,j], grid[3,2,j]) # bottom left
        r4n = SVector(grid[1,2,j+1], grid[2,2,j+1], grid[3,2,j+1]) # bottom right

        # also get chord length for setting finite core size
        cl = grid[1,end,j] - grid[1,1,j]
        cr = grid[1,end,j+1] - grid[1,1,j+1]
        c = (cl + cr)/2

        for i = 1:nc-1

            # grid corners of current panel
            r1 = r1n # top left
            r2 = r2n # top right
            r3 = r3n # bottom left
            r4 = r4n # bottom right

            # grid corners of next panel
            r1n = r3 # top left
            r2n = r4 # top right
            r3n = SVector(grid[1,i+2,j], grid[2,i+2,j], grid[3,i+2,j]) # bottom left
            r4n = SVector(grid[1,i+2,j+1], grid[2,i+2,j+1], grid[3,i+2,j+1]) # bottom right

            # top left corner of ring vortex
            rtl = linearinterp(0.25, r1, r3)

            # top right corner of ring vortex
            rtr = linearinterp(0.25, r2, r4)

            # center of top bound vortex
            rtc = (rtl + rtr)/2

            # bottom left corner of ring vortex
            rbl = linearinterp(0.25, r1n, r3n)

            # bottom right corner of ring vortex
            rbr = linearinterp(0.25, r2n, r4n)

            # center of bottom bound vortex
            rbc = (rbl + rbr)/2

            # control point
            rtop = linearinterp(0.5, r1, r2)
            rbot = linearinterp(0.5, r3, r4)
            rcp = linearinterp(0.75, rl, rr)

            # surface normal
            ncp = cross(rcp - rtr, rcp - rtl)
            ncp /= norm(ncp)

            # set finite core size
            Δs = sqrt((rtr[2]-rtl[2])^2 + (rtr[3]-rtl[3])^2)
            core_size = fcore(c, Δs)

            ipanel = mirror*right_side*nc*ns + (j-1)*nc + i
            panels[ipanel] = Ring{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)
        end

        # grid corners of current panel
        r1 = r1n # top left
        r2 = r2n # top right
        r3 = r3n # bottom left
        r4 = r4n # bottom right

        # top left corner of ring vortex
        rtl = linearinterp(0.25, r1, r3)

        # top right corner of ring vortex
        rtr = linearinterp(0.25, r2, r4)

        # center of top bound vortex
        rtc = (rtl + rtr)/2

        # bottom left corner of ring vortex
        rbl = r3

        # bottom right corner of ring vortex
        rbr = r4

        # center of bottom bound vortex
        rbc = (rbl + rbr)/2

        # control point
        rtop = linearinterp(0.5, r1, r2)
        rbot = linearinterp(0.5, r3, r4)
        rcp = linearinterp(0.75, rl, rr)

        # surface normal
        ncp = cross(rcp - rtr, rcp - rtl)
        ncp /= norm(ncp)

        # set finite core size
        Δs = sqrt((rtr[2]-rtl[2])^2 + (rtr[3]-rtl[3])^2)
        core_size = fcore(c, Δs)

        ip = i
        jp = mirror*right_side*ns + j
        panels[ip,jp] = Ring{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)
    end

    # other side
    if mirror
        for j = 1:ns
            for i = 1:nc
                ip = i
                jp1 = right_side*ns + j
                jp2 = 2*ns - right_side*ns - j + 1
                panels[ip,jp2] = reflect(panels[ip,jp1])
            end
        end
    end

    return panels
end

"""
    grid_to_horseshoe_vortices(xyz, ns, nc; mirror=false,
        fcore=(c, Δs)->max(c/4, Δs/2), spacing_s=Cosine(), spacing_c=Uniform(),
        interp_s=(x)->(x,y,xpt)->LinearInterpolation(x, y)(xpt),
        interp_c=(x,y,xpt)->LinearInterpolation(x, y)(xpt))

Discretize a potentially curved lifting surface defined by a grid with dimensions
(3, i, j) where `i` corresponds to the chordwise direction (ordered from leading
edge to trailing edge) and `j` corresponds to the spanwise direction (ordered from
left to right) into `ns` spanwise and `nc` chordwise horseshoe vortex panels
according to the spanwise discretization scheme `spacing_s` and chordwise
discretization scheme `spacing_c`.  The bound vortex will be placed at the
1/4 chord and the control point will be placed at the 3/4 chord of each panel.

In order to be consistent with vortex lattice theory using horseshoe vortices,
the y and z-coordinates of each resulting chordwise strip of panels are set
to the same y and z-coordinates as found at the leading edge of the chordwise
strip.  The normal vectors used to satisfy the no-flow boundary conditions,
however, are set using the original geometry.

# Arguments
 - `xyz`: grid of dimensions (3, i, j) where where `i` corresponds to the
    chordwise direction and `j` corresponds to the spanwise direction.
 - `ns`: number of spanwise panels
 - `nc`: number of chordwise panels
 - `mirror`:  mirror the geometry across the y-axis? defaults to `false`.
 - `fcore`: function for setting the finite core size based on the chord length
        (in the x-direction) and/or the panel width (in the y/z directions).
        Defaults to `(c, Δs)->max(c/4, Δs/2)` (which corresponds to AVL's implementation)
 - `spacing_s`: spanwise discretization scheme, defaults to `Cosine()`
 - `spacing_c`: chordwise discretization scheme, defaults to `Uniform()`
 - `interp_s`: spanwise interpolation function, defaults to linear interpolation
 - `interp_c`: chordwise interpolation function, defaults to linear interpolation
"""
function grid_to_horseshoe_vortices(xyz, ns, nc; mirror=false,
    fcore=(c, Δs)->max(c/4, Δs/2), spacing_s=Cosine(), spacing_c=Uniform(),
    interp_s=(x,y,xpt)->LinearInterpolation(x,y)(xpt),
    interp_c=(x,y,xpt)->LinearInterpolation(x,y)(xpt))

    TF = eltype(xyz)

    # get spanwise and chordwise spacing
    etas, etabar = spanwise_spacing(ns+1, spacing_s)
    eta_qtr, eta_thrqtr = chordwise_spacing(nc+1, spacing_c)

    # add leading and trailing edge to grid
    eta_qtr = vcat(0.0, eta_qtr, 1.0)

    # interpolate grid first along `i` then along `j`
    xyz_bound = interpolate_grid(xyz, eta_qtr, interp_c, xdir=1, ydir=1)
    xyz_corner = interpolate_grid(xyz_bound, etas, interp_s, ydir=2)
    xyz_center = interpolate_grid(xyz_bound, etabar, interp_s, ydir=2)

    xyz_cp = interpolate_grid(xyz, eta_thrqtr, interp_c, xdir=1, ydir=1)
    xyz_cp = interpolate_grid(xyz_cp, etabar, interp_s, ydir=2)

    # check which side we're working with
    right_side = sum(xyz[2,:,:]) > 0

    # initialize output
    panels = Matrix{Horseshoe{TF}}(undef, nc, (1+mirror)*ns)

    # populate each panel
    for j = 1:ns
        # get y and z-coordinates for this chordwise strip
        rly = xyz_corner[2,1,j]
        rlz = xyz_corner[3,1,j]

        rry = xyz_corner[2,1,j+1]
        rrz = xyz_corner[3,1,j+1]

        rcy = xyz_center[2,1,j]
        rcz = xyz_center[3,1,j]

        # also get chord length for setting the finite core size
        c = xyz_center[1,end,j] - xyz_center[1,1,j]

        for i = 1:nc
            # note that `i==1` corresponds to the leading edge for `xyz_corner` and `xyz_center`
            rl = SVector(xyz_corner[1,i+1,j], xyz_corner[2,i+1,j], xyz_corner[3,i+1,j])
            rc = SVector(xyz_center[1,i+1,j], xyz_center[2,i+1,j], xyz_center[3,i+1,j])
            rr = SVector(xyz_corner[1,i+1,j+1], xyz_corner[2,i+1,j+1], xyz_corner[3,i+1,j+1])
            rcp = SVector(xyz_cp[1,i,j], xyz_cp[2,i,j], xyz_cp[3,i,j])

            # get normal vector at the control point
            ncp = cross(rcp - rr, rcp - rl)
            ncp /= norm(ncp)

            # set y and z-coordinates to correspond to the leading edge
            rl = SVector(rl[1], rly, rlz)
            rr = SVector(rr[1], rry, rrz)
            rc = SVector(rc[1], rcy, rcz)
            rcp = SVector(rcp[1], rcy, rcz)

            # x-distance to the trailing edge
            xl_te = xyz_corner[1,end,j] - rl[1]
            xc_te = xyz_center[1,end,j] - rc[1]
            xr_te = xyz_corner[1,end,j+1] - rr[1]

            # set finite core size
            Δs = sqrt((rr[2]-rl[2])^2 + (rr[3]-rl[3])^2)
            core_size = fcore(c, Δs)

            ip = i
            jp = mirror*right_side*ns + j
            panels[ip, jp] = Horseshoe{TF}(rl, rc, rr, rcp, ncp, xl_te, xc_te, xr_te, core_size)
        end
    end

    # other side
    if mirror
        for j = 1:ns
            for i = 1:nc
                ip = i
                jp1 = right_side*ns + j
                jp2 = 2*ns - right_side*ns - j + 1
                panels[ip,jp2] = reflect(panels[ip,jp1])
            end
        end
    end

    return panels
end

"""
    grid_to_vortex_rings(xyz, ns, nc; mirror=false, fcore=(c, Δs)->max(c/4, Δs/2),
        spacing_s=Cosine(), spacing_c=Uniform(),
        interp_s=(x,y,xpt)->LinearInterpolation(x, y)(xpt),
        interp_c=(x,y,xpt)->LinearInterpolation(x, y)(xpt))

Discretize a potentially curved lifting surface defined by a grid with dimensions
(3, i, j) where `i` corresponds to the chordwise direction (ordered from leading
edge to trailing edge) and `j` corresponds to the spanwise direction (ordered from
left to right) into `ns` spanwise and `nc` chordwise vortex ring panels
according to the spanwise discretization scheme `spacing_s` and chordwise
discretization scheme `spacing_c`.  The bound vortex will be placed at the
1/4 chord and the control point will be placed at the 3/4 chord of each panel

# Arguments
 - `xyz`: grid of dimensions (3, i, j) where where `i` corresponds to the
    chordwise direction and `j` corresponds to the spanwise direction.
 - `ns`: number of spanwise panels
 - `nc`: number of chordwise panels
 - `mirror`:  mirror the geometry across the y-axis? defaults to `false`.
 - `fcore`: function for setting the finite core size based on the chord length
        (in the x-direction) and/or the panel width (in the y/z directions).
        Defaults to `(c, Δs)->max(c/4, Δs/2)` (which corresponds to AVL's implementation)
 - `spacing_s`: spanwise discretization scheme, defaults to `Cosine()`
 - `spacing_c`: chordwise discretization scheme, defaults to `Uniform()`
 - `interp_s`: spanwise interpolation function, defaults to linear interpolation
 - `interp_c`: chordwise interpolation function, defaults to linear interpolation
"""
function grid_to_vortex_rings(xyz, ns, nc; mirror=false,
    fcore=(c, Δs)->max(c/4, Δs/2), spacing_s=Cosine(), spacing_c=Uniform(),
    interp_s=(x,y,xpt)->LinearInterpolation(x, y)(xpt),
    interp_c=(x,y,xpt)->LinearInterpolation(x, y)(xpt))

    TF = eltype(xyz)

    # get spanwise and chordwise spacing
    etas, etabar = spanwise_spacing(ns+1, spacing_s)
    eta_qtr, eta_thrqtr = chordwise_spacing(nc+1, spacing_c)

    # add leading and trailing edge to grid
    eta_qtr = vcat(0.0, eta_qtr, 1.0)

    # interpolate grid first along `i` then along `j`
    xyz_bound = interpolate_grid(xyz, eta_qtr, interp_c, ydir=1)
    xyz_corner = interpolate_grid(xyz_bound, etas, interp_s, ydir=2)
    xyz_center = interpolate_grid(xyz_bound, etabar, interp_s, ydir=2)

    xyz_cp = interpolate_grid(xyz, eta_thrqtr, interp_c, ydir=1)
    xyz_cp = interpolate_grid(xyz_cp, etabar, interp_s, ydir=2)

    # check which side we're working with
    right_side = sum(xyz[2,:,:]) > 0

    # initialize output
    panels = Matrix{Ring{TF}}(undef, nc, (1+mirror)*ns)

    # populate each panel
    for j = 1:ns

        # get chord for setting finite core size
        c = xyz_center[1,end,j] - xyz_center[1,1,j]

        for i = 1:nc
            # note that `i==1` corresponds to the leading edge for `xyz_corner` and `xyz_center`
            rtl = SVector(xyz_corner[1,i+1,j], xyz_corner[2,i+1,j], xyz_corner[3,i+1,j])
            rtc = SVector(xyz_center[1,i+1,j], xyz_center[2,i+1,j], xyz_center[3,i+1,j])
            rtr = SVector(xyz_corner[1,i+1,j+1], xyz_corner[2,i+1,j+1], xyz_corner[3,i+1,j+1])
            rbl = SVector(xyz_corner[1,i+2,j], xyz_corner[2,i+2,j], xyz_corner[3,i+2,j])
            rbc = SVector(xyz_center[1,i+2,j], xyz_center[2,i+2,j], xyz_center[3,i+2,j])
            rbr = SVector(xyz_corner[1,i+2,j+1], xyz_corner[2,i+2,j+1], xyz_corner[3,i+2,j+1])
            rcp = SVector(xyz_cp[1,i,j], xyz_cp[2,i,j], xyz_cp[3,i,j])
            ncp = cross(rcp - rtr, rcp - rtl)
            ncp /= norm(ncp)

            # set finite core size
            Δs = sqrt((rtr[2]-rtl[2])^2 + (rtr[3]-rtl[3])^2)
            core_size = fcore(c, Δs)

            ip = i
            jp = mirror*right_side*ns + j
            panels[ip, jp] = Ring{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)
        end
    end

    # other side
    if mirror
        for j = 1:ns
            for i = 1:nc
                ip = i
                jp1 = right_side*ns + j
                jp2 = 2*ns - right_side*ns - j + 1
                panels[ip,jp2] = reflect(panels[ip,jp1])
            end
        end
    end

    return panels
end

"""
    wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, phi, ns, nc;
        fc=fill(x->0, length(xle)), mirror=false, fcore=(c, Δs)->max(c/4, Δs/2),
        spacing_s=Cosine(), spacing_c=Uniform(), interp_s=(x,y,xpt)->LinearInterpolation(x, y)(xpt))

Discretize a wing into `ns` spanwise and `nc` chordwise horseshoe vortex panels
according to the spanwise discretization scheme `spacing_s` and chordwise
discretization scheme `spacing_c`.

# Arguments
 - `xle`: leading edge x-coordinate of each airfoil section
 - `yle`: leading edge y-coordinate of each airfoil section
 - `zle`: leading edge z-coordinate of each airfoil section
 - `chord`: chord length of each airfoil section
 - `theta`: twist of each airfoil section
 - `phi`: dihedral angle of each airfoil section, defined by a right hand rotation about the x-axis
 - `ns`: number of spanwise panels
 - `nc`: number of chordwise panels
 - `fc`: (optional) camber line function y=f(x) of each airfoil section
 - `mirror`:  mirror the geometry across the y-axis?, defaults to `false`
 - `fcore`: function for setting the finite core size based on the chord length
        (in the x-direction) and/or the panel width (in the y/z directions).
        Defaults to `(c, Δs)->max(c/4, Δs/2)` (which corresponds to AVL's implementation)
 - `spacing_s`: spanwise discretization scheme, defaults to `Cosine()`
 - `spacing_c`: chordwise discretization scheme, defaults to `Uniform()`
 - `interp_s`: interpolation function between spanwise stations, defaults to linear interpolation
"""
wing_to_horseshoe_vortices

function wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, phi, ns, nc;
    fc=fill(x->0, length(xle)), mirror=false, fcore=(c, Δs)->max(c/4, Δs/2),
    spacing_s=Cosine(), spacing_c=Uniform(), interp_s=(x,y,xpt)->LinearInterpolation(x, y)(xpt))

    TF = promote_type(eltype(xle), eltype(yle), eltype(zle), eltype(chord), eltype(theta), eltype(phi))

    # get spanwise and chordwise spacing
    etas, etabar = spanwise_spacing(ns+1, spacing_s)
    eta_qtr, eta_thrqtr = chordwise_spacing(nc+1, spacing_c)

    # add leading and trailing edge to grid
    eta_qtr = vcat(0.0, eta_qtr, 1.0)

    # check input dimensions
    n = length(xle)
    @assert n == length(yle) == length(zle) == length(chord) == length(theta) == length(phi)

    # set bound vortex and control point chordwise locations
    xyz_bound = Array{TF, 3}(undef, 3, nc+2, n)
    xyz_cp = Array{TF, 3}(undef, 3, nc, n)
    for j = 1:n

        # rotation matrix for twist
        st, ct = sincos(theta[j])
        Rt = @SMatrix [ct 0 st; 0 1 0; -st 0 ct]

        # rotation matrix for dihedral
        sp, cp = sincos(phi[j])
        Rp = @SMatrix [1 0 0; 0 cp -sp; 0 sp cp]

        # location of leading edge
        rle = SVector(xle[j], yle[j], zle[j])

        # bound vortex chordwise locations (and leading and trailing edge)
        for i = 1:nc+2
            # bound vortex location
            xc = eta_qtr[i]
            zc = fc[j](xc)

            # location on airfoil
            r = SVector(xc, 0, zc)

            # scale by chord length
            r = chord[j] * r

            # apply twist
            r = Rt * r

            # apply dihedral
            r = Rp * r

            # add leading edge offset
            r = r + rle

            # store final location
            xyz_bound[:,i,j] = r
        end

        # control point chordwise locations
        for i = 1:nc
            # bound vortex location
            xc = eta_thrqtr[i]
            zc = fc[j](xc)

            # location on airfoil
            r = SVector(xc, 0, zc)

            # scale by chord length
            r = chord[j] * r

            # apply twist
            r = Rt * r

            # apply dihedral
            r = Rp * r

            # add leading edge offset
            r = r + rle

            # store final location
            xyz_cp[:,i,j] = r
        end
    end

    # interpolate the grids to match the specified spanwise spacing
    xyz_corner = interpolate_grid(xyz_bound, etas, interp_s; ydir=2)
    xyz_center = interpolate_grid(xyz_bound, etabar, interp_s; ydir=2)
    xyz_cp = interpolate_grid(xyz_cp, etabar, interp_s; ydir=2)

    # initialize output
    panels = Matrix{Horseshoe{TF}}(undef, nc, (1+mirror)*ns)

    # check which side we're working with
    right_side = sum(yle) > 0

    # populate each panel
    for j = 1:ns

        # get y and z-coordinates for this chordwise strip
        rly = xyz_corner[2,1,j]
        rlz = xyz_corner[3,1,j]

        rry = xyz_corner[2,1,j+1]
        rrz = xyz_corner[3,1,j+1]

        rcy = xyz_center[2,1,j]
        rcz = xyz_center[3,1,j]

        # also get chord for setting finite core size
        c = xyz_center[1,end,j] - xyz_center[1,1,j]

        for i = 1:nc
            # note that `i==1` corresponds to the leading edge for `xyz_corner` and `xyz_center`
            rl = SVector(xyz_corner[1,i+1,j], xyz_corner[2,i+1,j], xyz_corner[3,i+1,j])
            rc = SVector(xyz_center[1,i+1,j], xyz_center[2,i+1,j], xyz_center[3,i+1,j])
            rr = SVector(xyz_corner[1,i+1,j+1], xyz_corner[2,i+1,j+1], xyz_corner[3,i+1,j+1])
            rcp = SVector(xyz_cp[1,i,j], xyz_cp[2,i,j], xyz_cp[3,i,j])

            # get normal vector at the control point
            ncp = cross(rcp - rr, rcp - rl)
            ncp = ncp/norm(ncp)

            # set y and z-coordinates to correspond to the leading edge
            rl = SVector(rl[1], rly, rlz)
            rr = SVector(rr[1], rry, rrz)
            rc = SVector(rc[1], rcy, rcz)
            rcp = SVector(rcp[1], rcy, rcz)

            # x-distance to the trailing edge
            xl_te = xyz_corner[1,end,j] - rl[1]
            xc_te = xyz_center[1,end,j] - rc[1]
            xr_te = xyz_corner[1,end,j+1] - rr[1]

            # set finite core size
            Δs = sqrt((rr[2]-rl[2])^2 + (rr[3]-rl[3])^2)
            core_size = fcore(c, Δs)

            ip = i
            jp = mirror*right_side*ns + j
            panels[ip, jp] = Horseshoe{TF}(rl, rc, rr, rcp, ncp, xl_te, xc_te, xr_te, core_size)
        end
    end

    # other side
    if mirror
        for j = 1:ns
            for i = 1:nc
                ip = i
                jp1 = right_side*ns + j
                jp2 = 2*ns - right_side*ns - j + 1
                panels[ip,jp2] = reflect(panels[ip,jp1])
            end
        end
    end

    return panels
end

"""
    wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, ns, nc;
        fc=fill(x->0, length(xle)), mirror=false,
        fcore=(c, Δs)->max(c/4, Δs/2), spacing_s=Cosine(), spacing_c=Uniform(),
        interp_s=(x,y,xpt)->LinearInterpolation(x, y)(xpt))

Discretize a wing into `ns` spanwise and `nc` chordwise vortex ring panels
according to the spanwise discretization scheme `spacing_s` and chordwise
discretization scheme `spacing_c`.

# Arguments
 - `xle`: leading edge x-coordinate of each airfoil section
 - `yle`: leading edge y-coordinate of each airfoil section
 - `zle`: leading edge z-coordinate of each airfoil section
 - `chord`: chord length of each airfoil section
 - `theta`: twist of each airfoil section
 - `phi`: dihedral angle of each airfoil section, defined by a right hand rotation about the x-axis
 - `ns`: number of spanwise panels
 - `nc`: number of chordwise panels
 - `fc`: (optional) camber line function y=f(x) of each airfoil section
 - `mirror`:  mirror the geometry across the y-axis?, defaults to `false`
 - `fcore`: function for setting the finite core size based on the chord length
        (in the x-direction) and/or the panel width (in the y/z directions).
        Defaults to `(c, Δs)->max(c/4, Δs/2)` (which corresponds to AVL's implementation)
 - `spacing_s`: spanwise discretization scheme, defaults to `Cosine()`
 - `spacing_c`: chordwise discretization scheme, defaults to `Uniform()`
 - `interp_s`: interpolation function between spanwise stations, defaults to linear interpolation
"""
wing_to_vortex_rings

function wing_to_vortex_rings(xle, yle, zle, chord, theta, phi, ns, nc;
    fc=fill(x->0, length(xle)), mirror=false, fcore=(c, Δs)->max(c/4, Δs/2),
    spacing_s=Cosine(), spacing_c=Uniform(), interp_s=(x,y,xpt)->LinearInterpolation(x, y)(xpt))

    TF = promote_type(eltype(xle), eltype(yle), eltype(zle), eltype(chord), eltype(theta), eltype(phi))

    N = (1+mirror)*nc*ns

    # get spanwise and chordwise spacing
    etas, etabar = spanwise_spacing(ns+1, spacing_s)
    eta_qtr, eta_thrqtr = chordwise_spacing(nc+1, spacing_c)

    # add leading and trailing edge to grid
    eta_qtr = vcat(0.0, eta_qtr, 1.0)

    # check input dimensions
    n = length(xle)
    @assert n == length(yle) == length(zle) == length(chord) == length(theta) == length(phi)

    # set bound vortex and control point chordwise locations
    xyz_bound = Array{TF, 3}(undef, 3, nc+2, n)
    xyz_cp = Array{TF, 3}(undef, 3, nc, n)
    for j = 1:n

        # rotation matrix for twist
        st, ct = sincos(theta[j])
        Rt = @SMatrix [ct 0 st; 0 1 0; -st 0 ct]

        # rotation matrix for dihedral
        sp, cp = sincos(phi[j])
        Rp = @SMatrix [1 0 0; 0 cp -sp; 0 sp cp]

        # location of leading edge
        rle = SVector(xle[j], yle[j], zle[j])

        # bound vortex chordwise locations (and trailing edge)
        for i = 1:nc+2
            # bound vortex location
            xc = eta_qtr[i]
            zc = fc[j](xc)

            # location on airfoil
            r = SVector(xc, 0, zc)

            # scale by chord length
            r = chord[j] * r

            # apply twist
            r = Rt * r

            # apply dihedral
            r = Rp * r

            # add leading edge offset
            r = r + rle

            # store final location
            xyz_bound[:,i,j] = r
        end

        # control point chordwise locations
        for i = 1:nc
            # bound vortex location
            xc = eta_thrqtr[i]
            zc = fc[j](xc)

            # location on airfoil
            r = SVector(xc, 0, zc)

            # scale by chord length
            r = chord[j] * r

            # apply twist
            r = Rt * r

            # apply dihedral
            r = Rp * r

            # add leading edge offset
            r = r + rle

            # store final location
            xyz_cp[:,i,j] = r
        end
    end

    # interpolate the grids to match the specified spanwise spacing
    xyz_corner = interpolate_grid(xyz_bound, etas, interp_s; ydir=2)
    xyz_center = interpolate_grid(xyz_bound, etabar, interp_s; ydir=2)
    xyz_cp = interpolate_grid(xyz_cp, etabar, interp_s; ydir=2)

    # initialize output
    panels = Matrix{Ring{TF}}(undef,  nc, (1+mirror)*ns)

    # check which side we're working with
    right_side = sum(yle) > 0

    # populate each panel
    for j = 1:ns

        # get chord for setting finite core size
        c = xyz_center[1,end,j] - xyz_center[1,1,j]

        for i = 1:nc
            # note that `i==1` corresponds to the leading edge for `xyz_corner` and `xyz_center`
            rtl = SVector(xyz_corner[1,i+1,j], xyz_corner[2,i+1,j], xyz_corner[3,i+1,j])
            rtc = SVector(xyz_center[1,i+1,j], xyz_center[2,i+1,j], xyz_center[3,i+1,j])
            rtr = SVector(xyz_corner[1,i+1,j+1], xyz_corner[2,i+1,j+1], xyz_corner[3,i+1,j+1])
            rbl = SVector(xyz_corner[1,i+2,j], xyz_corner[2,i+2,j], xyz_corner[3,i+2,j])
            rbc = SVector(xyz_center[1,i+2,j], xyz_center[2,i+2,j], xyz_center[3,i+2,j])
            rbr = SVector(xyz_corner[1,i+2,j+1], xyz_corner[2,i+2,j+1], xyz_corner[3,i+2,j+1])
            rcp = SVector(xyz_cp[1,i,j], xyz_cp[2,i,j], xyz_cp[3,i,j])
            ncp = cross(rcp - rtr, rcp - rtl)
            ncp /= norm(ncp)

            # set finite core size
            Δs = sqrt((rtr[2]-rtl[2])^2 + (rtr[3]-rtl[3])^2)
            core_size = fcore(c, Δs)

            ip = i
            jp = mirror*right_side*ns + j
            panels[ip, jp] = Ring{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp, core_size)
        end
    end

    # other side
    if mirror
        for j = 1:ns
            for i = 1:nc
                ip = i
                jp1 = right_side*ns + j
                jp2 = 2*ns - right_side*ns - j + 1
                panels[ip,jp2] = reflect(panels[ip,jp1])
            end
        end
    end

    return panels
end
