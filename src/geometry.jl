"""
    AbstractSpacing

Spacing discretization scheme supertype
"""
abstract type AbstractSpacing end

"""
    Uniform()

Uniform discretization scheme.
"""
struct Uniform <: AbstractSpacing end

"""
    Sine()

Sine-spaced discretization scheme.  Using sine-spacing on the right half of a
wing effectively results in cosine spacing once symmetry is applied.
"""
struct Sine <: AbstractSpacing end

"""
    Cosine()

Cosine-spaced discretization scheme.  This is typically one of the most accurate
spacing schemes for spanwise spacing.
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
    chordwise_spacing(n, spacing::AbstractSpacing)

Distribute `n` panel edge, `n-1` vortex, and `n-1` control point chordwise
locations on the interval between 0 and 1 according to the discretization
strategy in `spacing`.
"""
chordwise_spacing

# uniform
function chordwise_spacing(n, ::Uniform)

    eta = range(0, 1.0, length=n)

    eta_qtr = linearinterp(0.25, eta[1:end-1], eta[2:end])
    eta_thrqtr = linearinterp(0.75, eta[1:end-1], eta[2:end])

    return eta, eta_qtr, eta_thrqtr
end

# sine
function chordwise_spacing(n, ::Sine)

    theta = range(0, pi/2, length=n)
    eta = sin.(theta)

    eta_qtr = linearinterp(0.25, eta[1:end-1], eta[2:end])
    eta_thrqtr = linearinterp(0.75, eta[1:end-1], eta[2:end])

    return eta, eta_qtr, etq_thrqtr
end

# cosine
function chordwise_spacing(n, ::Cosine)

    theta = range(0, pi, length=n)
    eta = (1.0 .- cos.(theta))/2.0

    eta_qtr = linearinterp(0.25, eta[1:end-1], eta[2:end])
    eta_thrqtr = linearinterp(0.75, eta[1:end-1], eta[2:end])

    return eta, eta_qtr, eta_thrqtr
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
    grid_to_surface_panels(xyz; mirror = false, fcore = (c, Δs) -> 1e-3)

Construct a set of panels with associated vortex rings given a potentially curved
lifting surface defined by a grid with dimensions (3, i, j) where `i` corresponds
to the chordwise direction (ordered from leading edge to trailing edge) and `j`
corresponds to the spanwise direction (ordered from left to right).  The leading
edge of each ring vortex will be placed at the 1/4 chord and the control point
will be placed at the 3/4 chord of each panel.

Return a grid with dimensions (3, i, j) containing the panel corners and a matrix
with dimensions (i, j) containing the generated panels.

# Keyword Arguments
- `mirror`:  mirror the geometry across the X-Z plane? defaults to `false`.
- `fcore`: function for setting the finite core size based on the chord length
       (in the x-direction) and/or the panel width (in the y/z directions).
       Defaults to `(c, Δs) -> 1e-3`
"""
function grid_to_surface_panels(xyz;
    mirror = false,
    fcore = (c, Δs) -> 1e-3)

    TF = eltype(xyz)

    nc = size(xyz, 2)-1 # number of chordwise panels
    ns = size(xyz, 3)-1 # number of spanwise panels
    N = (1+mirror)*nc*ns # total number of panels

    # check which side we're working with
    right_side = sum(xyz[2,:,:]) > 0

    # initialize output
    xyz_panels = copy(xyz)
    panels = Matrix{SurfacePanel{TF}}(undef, nc, (1+mirror)*ns)

    # populate each panel
    for j = 1:ns

        # get leading edge panel corners
        r1n = SVector(xyz[1,1,j], xyz[2,1,j], xyz[3,1,j]) # top left
        r2n = SVector(xyz[1,1,j+1], xyz[2,1,j+1], xyz[3,1,j+1]) # top right
        r3n = SVector(xyz[1,2,j], xyz[2,2,j], xyz[3,2,j]) # bottom left
        r4n = SVector(xyz[1,2,j+1], xyz[2,2,j+1], xyz[3,2,j+1]) # bottom right

        # also get chord length for setting finite core size
        cl = norm(xyz[:,end,j] - xyz[:,1,j])
        cr = norm(xyz[:,end,j+1] - xyz[:,1,j+1])
        c = (cl + cr)/2

        for i = 1:nc-1

            # xyz corners of current panel
            r1 = r1n # top left of panel
            r2 = r2n # top right of panel
            r3 = r3n # bottom left of panel
            r4 = r4n # bottom right of panel

            # xyz corners of next panel
            r1n = r3 # top left
            r2n = r4 # top right
            r3n = SVector(xyz[1,i+2,j], xyz[2,i+2,j], xyz[3,i+2,j]) # bottom left
            r4n = SVector(xyz[1,i+2,j+1], xyz[2,i+2,j+1], xyz[3,i+2,j+1]) # bottom right

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
            rcp = linearinterp(0.75, rtop, rbot)

            # surface normal
            ncp = cross(rcp - rtr, rcp - rtl)
            ncp /= norm(ncp)

            # set finite core size
            Δs = sqrt((rtr[2]-rtl[2])^2 + (rtr[3]-rtl[3])^2)
            core_size = fcore(c, Δs)

            # get chord length of current panel
            chord = norm((r1 + r2)/2 - (r3 + r4)/2)

            ip = i
            jp = mirror*right_side*ns + j
            panels[ip, jp] = SurfacePanel{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp,
                ncp, core_size, chord)
        end

        # xyz corners of current panel
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
        rcp = linearinterp(0.75, rtop, rbot)

        # surface normal
        ncp = cross(rcp - rtr, rcp - rtl)
        ncp /= norm(ncp)

        # set finite core size
        Δs = sqrt((rtr[2]-rtl[2])^2 + (rtr[3]-rtl[3])^2)
        core_size = fcore(c, Δs)

        # get chord length of current panel
        chord = norm((r1 + r2)/2 - (r3 + r4)/2)

        ip = nc
        jp = mirror*right_side*ns + j
        panels[ip,jp] = SurfacePanel{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp,
            core_size, chord)
    end

    # other side
    if mirror

        # first reflect xyz
        if right_side
            xyz_l = reverse(xyz_panels, dims=3)
            xyz_l[2,:,:] .= -xyz_l[2,:,:]
            xyz_r = xyz_panels
        else
            xyz_l = xyz_panels
            xyz_r = reverse(xyz_panels, dims=3)
            xyz_r[2,:,:] .= -xyz_r[2,:,:]
        end
        xyz_panels = cat(xyz_l, xyz_r[:,:,2:end], dims=3)

        # then populate remaining panels
        for j = 1:ns
            for i = 1:nc
                ip = i
                jp1 = right_side*ns + j
                jp2 = 2*ns - right_side*ns - j + 1
                panels[ip,jp2] = reflect(panels[ip,jp1])
            end
        end
    end

    return xyz_panels, panels
end

"""
    grid_to_surface_panels(xyz, ns, nc;
        mirror = false,
        fcore = (c, Δs) -> 1e-3,
        spacing_s = Cosine(),
        spacing_c = Uniform(),
        interp_s = (x, y, xpt) -> LinearInterpolation(x, y)(xpt),
        interp_c = (x, y, xpt) -> LinearInterpolation(x, y)(xpt))

Discretize a potentially curved lifting surface defined by a grid with dimensions
(3, i, j) where `i` corresponds to the chordwise direction (ordered from leading
edge to trailing edge) and `j` corresponds to the spanwise direction (ordered from
left to right) into `ns` spanwise and `nc` chordwise panels with associated vortex
rings according to the spanwise discretization scheme `spacing_s` and chordwise
discretization scheme `spacing_c`.  The bound vortex will be placed at the
1/4 chord and the control point will be placed at the 3/4 chord of each panel.

Return a grid with dimensions (3, i, j) containing the interpolated panel corners
and a matrix with dimensions (i, j) containing the generated panels.

# Arguments
 - `xyz`: grid of dimensions (3, i, j) where where `i` corresponds to the
    chordwise direction and `j` corresponds to the spanwise direction.
 - `ns`: number of spanwise panels
 - `nc`: number of chordwise panels

# Keyword Arguments
 - `mirror`:  mirror the geometry across the X-Z plane? defaults to `false`.
 - `fcore`: function for setting the finite core size based on the chord length
        (in the x-direction) and/or the panel width (in the y/z directions).
        Defaults to `(c, Δs) -> 1e-3`
 - `spacing_s`: spanwise discretization scheme, defaults to `Cosine()`
 - `spacing_c`: chordwise discretization scheme, defaults to `Uniform()`
 - `interp_s`: spanwise interpolation function, defaults to linear interpolation
 - `interp_c`: chordwise interpolation function, defaults to linear interpolation
"""
function grid_to_surface_panels(xyz, ns, nc;
    mirror = false,
    fcore = (c, Δs) -> 1e-3,
    spacing_s = Cosine(),
    spacing_c = Uniform(),
    interp_s = (x, y, xpt) -> LinearInterpolation(x, y)(xpt),
    interp_c = (x, y, xpt) -> LinearInterpolation(x, y)(xpt))

    TF = eltype(xyz)

    # get spanwise and chordwise spacing
    etas, etabar = spanwise_spacing(ns+1, spacing_s)
    etac, eta_qtr, eta_thrqtr = chordwise_spacing(nc+1, spacing_c)

    # add leading and trailing edge to grid
    eta_qtr = vcat(0.0, eta_qtr, 1.0)

    # interpolate grid first along `i` then along `j`
    xyz_edge = interpolate_grid(xyz, etac, interp_c, ydir=1)
    xyz_panels = interpolate_grid(xyz_edge, etas, interp_s, ydir=2)

    xyz_bound = interpolate_grid(xyz, eta_qtr, interp_c, ydir=1)
    xyz_corner = interpolate_grid(xyz_bound, etas, interp_s, ydir=2)
    xyz_center = interpolate_grid(xyz_bound, etabar, interp_s, ydir=2)

    xyz_cp = interpolate_grid(xyz, eta_thrqtr, interp_c, ydir=1)
    xyz_cp = interpolate_grid(xyz_cp, etabar, interp_s, ydir=2)

    # check which side we're working with
    right_side = sum(xyz[2,:,:]) > 0

    # initialize output
    panels = Matrix{SurfacePanel{TF}}(undef, nc, (1+mirror)*ns)

    # populate each panel
    for j = 1:ns

        # get chord for setting finite core size

        # also get chord length for setting finite core size
        c = norm(xyz_center[:,end,j] - xyz_center[:,1,j])

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

            # set chord length of current panel
            r1 = xyz_panels[:, i, j] # top left
            r2 = xyz_panels[:, i, j+1] # top right
            r3 = xyz_panels[:, i+1, j] # bottom left
            r4 = xyz_panels[:, i+1, j+1] # bottom right
            chord = norm((r1 + r2)/2 - (r3 + r4)/2)

            ip = i
            jp = mirror*right_side*ns + j
            panels[ip, jp] = SurfacePanel{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp,
                ncp, core_size, chord)
        end
    end

    # other side
    if mirror
        # first reflect grid
        if right_side
            xyz_l = reverse(xyz_panels, dims=3)
            xyz_l[2,:,:] .= -xyz_l[2,:,:]
            xyz_r = xyz_panels
        else
            xyz_l = xyz_panels
            xyz_r = reverse(xyz_panels, dims=3)
            xyz_r[2,:,:] .= -xyz_r[2,:,:]
        end
        xyz_panels = cat(xyz_l, xyz_r[:,:,2:end], dims=3)

        # then populate remaining panels
        for j = 1:ns
            for i = 1:nc
                ip = i
                jp1 = right_side*ns + j
                jp2 = 2*ns - right_side*ns - j + 1
                panels[ip,jp2] = reflect(panels[ip,jp1])
            end
        end
    end

    return xyz_panels, panels
end

"""
    wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
        fc = fill(x -> 0, length(xle)),
        mirror = false,
        fcore = (c, Δs) -> 1e-3,
        spacing_s = Cosine(),
        spacing_c = Uniform(),
        interp_s = (x, y, xpt) -> LinearInterpolation(x, y)(xpt))

Discretize a wing into `ns` spanwise and `nc` chordwise panels with associated
vortex rings according to the spanwise discretization scheme `spacing_s` and
chordwise discretization scheme `spacing_c`.

Return a grid with dimensions (3, i, j) containing the panel corners and a matrix
with dimensions (i, j) containing the generated panels.

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
 - `mirror`:  mirror the geometry across the X-Z plane?, defaults to `false`
 - `fcore`: function for setting the finite core size based on the chord length
        (in the x-direction) and/or the panel width (in the y/z directions).
        Defaults to `(c, Δs) -> 1e-3`
 - `spacing_s`: spanwise discretization scheme, defaults to `Cosine()`
 - `spacing_c`: chordwise discretization scheme, defaults to `Uniform()`
 - `interp_s`: interpolation function between spanwise stations, defaults to linear interpolation
"""
function wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;
    fc = fill(x->0, length(xle)),
    mirror = false,
    fcore = (c, Δs) -> 1e-3,
    spacing_s = Cosine(),
    spacing_c = Uniform(),
    interp_s = (x, y, xpt) -> LinearInterpolation(x, y)(xpt))

    TF = promote_type(eltype(xle), eltype(yle), eltype(zle), eltype(chord), eltype(theta), eltype(phi))

    N = (1+mirror)*nc*ns

    # get spanwise and chordwise spacing
    etas, etabar = spanwise_spacing(ns+1, spacing_s)
    etac, eta_qtr, eta_thrqtr = chordwise_spacing(nc+1, spacing_c)

    # add leading and trailing edge to grid
    eta_qtr = vcat(0.0, eta_qtr, 1.0)

    # check input dimensions
    n = length(xle)
    @assert n == length(yle) == length(zle) == length(chord) == length(theta) == length(phi)

    # set bound vortex and control point chordwise locations
    xyz_edge = Array{TF, 3}(undef, 3, nc+1, n)
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

        # panel edge chordwise locations
        for i = 1:nc+1
            # bound vortex location
            xc = etac[i]
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
            xyz_edge[:,i,j] = r
        end

        # bound vortex chordwise locations (and leading and trailing edges)
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
    xyz_panels = interpolate_grid(xyz_edge, etas, interp_s; ydir=2)
    xyz_corner = interpolate_grid(xyz_bound, etas, interp_s; ydir=2)
    xyz_center = interpolate_grid(xyz_bound, etabar, interp_s; ydir=2)
    xyz_cp = interpolate_grid(xyz_cp, etabar, interp_s; ydir=2)

    # initialize output
    panels = Matrix{SurfacePanel{TF}}(undef,  nc, (1+mirror)*ns)

    # check which side we're working with
    right_side = sum(yle) > 0

    # populate each panel
    for j = 1:ns

        # get chord for setting finite core size
        c = norm(xyz_center[:,end,j] - xyz_center[:,1,j])

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

            # set chord length of current panel
            r1 = xyz_panels[:,i,j] # top left
            r2 = xyz_panels[:,i,j+1] # top right
            r3 = xyz_panels[:,i+1, j] # bottom left
            r4 = xyz_panels[:,i+1, j+1] # bottom right
            chord = norm((r1 + r2)/2 - (r3 + r4)/2)

            ip = i
            jp = mirror*right_side*ns + j
            panels[ip, jp] = SurfacePanel{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp,
                ncp, core_size, chord)
        end
    end

    # other side
    if mirror
        # first reflect grid
        if right_side
            xyz_l = reverse(xyz_panels, dims=3)
            xyz_l[2,:,:] .= -xyz_l[2,:,:]
            xyz_r = xyz_panels
        else
            xyz_l = xyz_panels
            xyz_r = reverse(xyz_panels, dims=3)
            xyz_r[2,:,:] .= -xyz_r[2,:,:]
        end
        xyz_panels = cat(xyz_l, xyz_r[:,:,2:end], dims=3)

        # then populate remaining panels
        for j = 1:ns
            for i = 1:nc
                ip = i
                jp1 = right_side*ns + j
                jp2 = 2*ns - right_side*ns - j + 1
                panels[ip,jp2] = reflect(panels[ip,jp1])
            end
        end
    end

    return xyz_panels, panels
end

"""
    update_surface_panels!(surface, grid; fcore = (c, Δs) -> 1e-3)

Updates the surface panels in `surface` to correspond to the grid coordinates in
`grid`.
"""
function update_surface_panels!(surface, grid; fcore = (c, Δs) -> 1e-3)

    TF = eltype(eltype(surface))

    nc, ns = size(surface) # number of chordwise and spanwise panels

    # populate each panel
    for j = 1:ns

        # get leading edge panel corners
        r1n = SVector(grid[1,1,j], grid[2,1,j], grid[3,1,j]) # top left
        r2n = SVector(grid[1,1,j+1], grid[2,1,j+1], grid[3,1,j+1]) # top right
        r3n = SVector(grid[1,2,j], grid[2,2,j], grid[3,2,j]) # bottom left
        r4n = SVector(grid[1,2,j+1], grid[2,2,j+1], grid[3,2,j+1]) # bottom right

        # also get chord length for setting finite core size
        cl = norm(grid[:,end,j] - grid[:,1,j])
        cr = norm(grid[:,end,j+1] - grid[:,1,j+1])
        c = (cl + cr)/2

        for i = 1:nc-1

            # grid corners of current panel
            r1 = r1n # top left of panel
            r2 = r2n # top right of panel
            r3 = r3n # bottom left of panel
            r4 = r4n # bottom right of panel

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
            rcp = linearinterp(0.75, rtop, rbot)

            # surface normal
            ncp = cross(rcp - rtr, rcp - rtl)
            ncp /= norm(ncp)

            # set finite core size
            Δs = sqrt((rtr[2]-rtl[2])^2 + (rtr[3]-rtl[3])^2)
            core_size = fcore(c, Δs)

            # get chord length of current panel
            chord = norm((r1 + r2)/2 - (r3 + r4)/2)

            surface[i,j] = SurfacePanel{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp,
                ncp, core_size, chord)
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
        rcp = linearinterp(0.75, rtop, rbot)

        # surface normal
        ncp = cross(rcp - rtr, rcp - rtl)
        ncp /= norm(ncp)

        # set finite core size
        Δs = sqrt((rtr[2]-rtl[2])^2 + (rtr[3]-rtl[3])^2)
        core_size = fcore(c, Δs)

        # get chord length of current panel
        chord = norm((r1 + r2)/2 - (r3 + r4)/2)

        surface[end, j] = SurfacePanel{TF}(rtl, rtc, rtr, rbl, rbc, rbr, rcp, ncp,
            core_size, chord)
    end

    return surface
end

"""
    lifting_line_geometry(grids, xc=0.25)
Construct a lifting line representation of the surfaces in `grids` at the
normalized chordwise location `xc`.  Return the lifting line coordinates
and chord lengths.
# Arguments
 - `grids`: Vector with length equal to the number of surfaces.  Each element of
    the vector is a grid with shape (3, nc, ns) which defines the discretization
    of a surface into panels. `nc` is the number of chordwise panels and `ns` is
    the number of spanwise panels.
 - `xc`: Normalized chordwise location of the lifting line from the leading edge.
    Defaults to the quarter chord
# Return Arguments:
 - `r`: Vector with length equal to the number of surfaces, with each element
    being a matrix with size (3, ns+1) which contains the x, y, and z coordinates
    of the resulting lifting line coordinates
 - `c`: Vector with length equal to the number of surfaces, with each element
    being a vector of length `ns+1` which contains the chord lengths at each
    lifting line coordinate.
"""
function lifting_line_geometry(grids, xc=0.25)
    TF = eltype(eltype(grids))
    nsurf = length(grids)
    r = Vector{Matrix{TF}}(undef, nsurf)
    c = Vector{Vector{TF}}(undef, nsurf)
    for isurf = 1:nsurf
        ns = size(grids[isurf], 3) - 1
        r[isurf] = Matrix{TF}(undef, 3, ns+1)
        c[isurf] = Vector{TF}(undef, ns+1)
    end
    return lifting_line_geometry!(r, c, grids, xc)
end

"""
    lifting_line_geometry!(r, c, grids, xc=0.25)
In-place version of [`lifting_line_geometry`](@ref)
"""
function lifting_line_geometry!(r, c, grids, xc=0.25)
    nsurf = length(grids)
    # iterate through each lifting surface
    for isurf = 1:nsurf
        # extract current grid
        grid = grids[isurf]
        # dimensions of this grid
        nc = size(grid, 2)
        ns = size(grid, 3)
        # loop through each spanwise section
        for j = 1:ns
            # get leading and trailing edges
            le = SVector(grid[1,1,j], grid[2,1,j], grid[3,1,j])
            te = SVector(grid[1,end,j], grid[2,end,j], grid[3,end,j])
            # get quarter-chord
            r[isurf][:,j] = linearinterp(xc, le, te)
            # get chord length
            c[isurf][j] = norm(le - te)
        end
    end
    return r, c
end

"""
    translate(grid, r)

Return a copy of the grid points in `grid` translated the distance specified by vector `r`
"""
translate(grid, r) = grid .+ r

"""
    translate!(grid, r)

Translate the grid points in `grid` the distance specified by vector `r`
"""
function translate!(grid, r)
    grid .+= r
    return grid
end

"""
    rotate(grid, R, r = [0,0,0])

Return a copy of the grid points in `grid` rotated about point `r` using the rotation matrix `R`
"""
rotate(grid, args...) = rotate!(copy(grid), args...)

"""
    rotate!(grid, R, r = [0,0,0])

Rotate the grid points in `grid` about point `r` using the rotation matrix `R`
"""
function rotate!(grid, R, r = (@SVector zeros(3)))
    nc, ns = size(grid, 2), size(grid, 3)

    grid .-= r

    for i = 1:nc, j = 1:ns
        grid[:,i,j] = R*grid[:,i,j]
    end

    grid .+= r

    return grid
end

"""
    trailing_edge_points(surface[s])

Return points on the trailing edge of each surface.
"""
trailing_edge_points(surfaces) = [bottom_left.(surfaces[end, :])..., bottom_right(surfaces[end, end])]

"""
    repeated_trailing_edge_points(surface[s])

Generates a dictionary of the form `Dict((isurf, i) => [(jsurf1, j1), (jsurf2, j2)...]` which
defines repeated trailing edge points.  Trailing edge point `i` on surface
`isurf` is repeated on surface `jsurf1` at point `j1`, `jsurf2` at point `j2`,
and so forth.
"""
repeated_trailing_edge_points

repeated_trailing_edge_points(surface::AbstractMatrix) = repeated_trailing_edge_points([surface])

function repeated_trailing_edge_points(surfaces::AbstractVector{<:AbstractMatrix})

    nsurf = length(surfaces)
    points = trailing_edge_points.(surfaces)
    repeated_points = Dict{NTuple{2, Int}, Vector{NTuple{2, Int}}}()

    # loop through all points
    for isurf = 1:nsurf, i = 1:length(points[isurf])
        # start with empty vector
        repeats = NTuple{2, Int}[]
        # add all repeated points to the vector
        for jsurf = 1:nsurf, j = 1:length(points[jsurf])
            if (isurf, i) != (jsurf, j) && isapprox(points[isurf][i], points[jsurf][j])
                push!(repeats, (jsurf, j))
            end
        end
        # only store vector if nonempty
        if !isempty(repeats)
            repeated_points[(isurf, i)] = repeats
        end
    end

    return repeated_points
end

"""
    flipy(r)

Flip sign of y-component of vector `r` (used for symmetry)
"""
@inline flipy(r) = SVector{3}(r[1], -r[2], r[3])

"""
    on_symmetry_plane(args...; tol=eps())

Test whether points in `args` are on the symmetry plane (y = 0)
"""
on_symmetry_plane

on_symmetry_plane(r; tol=eps()) = isapprox(r[2], 0.0, atol=tol)

@inline function on_symmetry_plane(r, args...; tol=eps())
    return on_symmetry_plane(r; tol=tol) && on_symmetry_plane(args...; tol=tol)
end

"""
    not_on_symmetry_plane(args...; tol=eps())

Test whether and of the points in `args` are not on the symmetry plane (y = 0)
"""
not_on_symmetry_plane(args...; tol=eps()) = !on_symmetry_plane(args...; tol=tol)
