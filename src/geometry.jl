abstract type AbstractSpacing end

"""
    Uniform()

Uniform discretization scheme
"""
struct Uniform <: AbstractSpacing end

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
    spanwise_spacing(n, spacing_type::AbstractSpacing)

Distribute `n` panel endpoints and `n-1` panel midpoints on the interval between
0 and 1 according to the discretization strategy in `spacing_type`.
"""
function spanwise_spacing(n, ::Uniform)

    eta = range(0, 1.0, length=n)
    eta_mid = linearinterp(0.5, eta[1:end-1], eta[2:end])

    return eta, eta_mid
end

function spanwise_spacing(n, ::Cosine)
    theta = range(0, pi, length=n)
    eta = (1.0 .- cos.(theta))/2.0

    # note that control points are also placed with cosine spacing as this improves accuracy tremendously
    theta_mid = linearinterp(0.5, theta[1:end-1], theta[2:end])
    eta_mid = (1.0 .- cos.(theta_mid))/2.0

    return eta, eta_mid
end

"""
    spanwise_spacing(n, spacing_type::AbstractSpacing)

Distribute `n-1` vortex and `n-1` control point chordwise locations on the
interval between 0 and 1 according to the discretization strategy in `spacing_type`.
"""
chordwise_spacing

# uniform
function chordwise_spacing(n, ::Uniform)

    eta = range(0, 1.0, length=n)

    eta_qtr = linearinterp(0.25, eta[1:end-1], eta[2:end])
    eta_thrqtr = linearinterp(0.75, eta[1:end-1], eta[2:end])

    return eta_qtr, eta_thrqtr
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
    interpolate_grid(xyz, eta, interp; dir=1 )

Interpolates the grid `xyz` along direction `dir`

# Arguments
 - `xyz`: Grid of size (3, ni, nj)
 - `eta`: New (normalized) coordinates in direction `dir` (0 <= eta <= 1)
 - `interp`: Interpolation method of form `ypt = f(x,y,xpt)`
 - `dir`: Direction in which to interpolate `xyz` (`i=1`, `j=2`)
"""
function interpolate_grid(xyz, eta, interp; dir=1)

    dim = dir + 1

    ni = size(xyz, dim)

    xyz = mapslices(xyz; dims=[1,dim]) do xyz_i

        # get distance between each grid location
        ds = [norm(xyz_i[:,k] - xyz_i[:,max(1,k-1)]) for k = 1:ni]

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
    grid_to_horseshoe_vortices(xyz, theta; mirror=false)

Construct a set of horseshoe vortex panels given a pre-existing set of panels defined
by a grid with dimensions (3, i, j) where `i` corresponds to the chordwise
direction and `j` corresponds to the spanwise direction.  The bound vortex for
each panel will be placed at the 1/4 chord and the control point will be placed at
the 3/4 chord.

Each chordwise strip must have the same y and z-coordinates (i.e. the airfoil geometry
cannot be cambered or twisted).  Local angle of attack at each spanwise increment
is specified with the input vector `theta` which has length `size(xyz, 3)`.

The optional argument `mirror` may be used to mirror the geometry across the
y-axis.

The resulting vector of panels can be reshaped to correspond to the original
grid with `reshape(panels, ni-1, nj-1)`  (`reshape(panels, nc, 2*ns)` if the
geometry is mirrored).
"""
function grid_to_horseshoe_vortices(xyz, theta; mirror=false)

    TF = promote_type(eltype(xyz), eltype(theta))

    nc = size(grid, 2)-1
    ns = size(grid, 3)-1
    N = (1+mirror)*nc*ns

    # check which side we're working with
    right_side = sum(xyz[2,:,:]) > 0
    left_side = !right_side

    # initialize output
    panels = Vector{Horseshoe{TF}}(undef, N)

    # populate each panel
    for j = 1:ns

        # angle
        th = linearinterp(0.5, theta[j], theta[j+1])

        for i = 1:nc

            # grid corners
            r1 = SVector(grid[1,i,j], grid[2,i,j], grid[3,i,j]) # top left
            r2 = SVector(grid[1,i,j+1], grid[2,i,j+1], grid[3,i,j+1]) # top right
            r3 = SVector(grid[1,i+1,j], grid[2,i+1,j], grid[3,i+1,j]) # bottom left
            r4 = SVector(grid[1,i+1,j+1], grid[2,i+1,j+1], grid[3,i+1,j+1]) # bottom right

            # left side of bound vortex
            rl = linearinterp(0.25, r1, r3)

            # right side of bound vortex
            rr = linearinterp(0.25, r2, r4)

            # control point
            rtop = linearinterp(0.5, r1, r2)
            rbot = linearinterp(0.5, r3, r4)
            rcp = linearinterp(0.75, rtop, rbot)

            ipanel = mirror*right_side*nc*ns + (j-1)*nc + i
            panels[ipanel] = Horseshoe{TF}(rl, rr, rcp, th)
        end
    end

    # other side
    if mirror
        for j = 1:ns
            for i = 1:nc
                ipanel1 = right_side*nc*ns + (j-1)*nc + i
                ipanel2 = left_side*nc*ns + (ns-j)*nc + i
                panels[ipanel2] = reflect_y(panels[ipanel1])
            end
        end
    end

    return panels
end


"""
    grid_to_vortex_rings(xyz; mirror=false)

Construct a set of vortex ring panels given a pre-existing set of panels defined
by a grid with dimensions (3, i, j) where `i` corresponds to the chordwise
direction and `j` corresponds to the spanwise direction.  The start of each ring
vortex will be placed at the 1/4 chord and the control point will be placed at
the 3/4 chord of each panel.

The resulting vector of panels can be reshaped to correspond to the original
grid with `reshape(panels, ni-1, nj-1)`.  (`reshape(panels, ni-1, 2*(nj-1))` if
the geometry is mirrored).
"""
function grid_to_vortex_rings(xyz; mirror=false)

    TF = eltype(xyz)

    nc = size(grid, 2)-1
    ns = size(grid, 3)-1
    N = (1+mirror)*nc*ns

    # check which side we're working with
    right_side = sum(xyz[2,:,:]) > 0
    left_side = !right_side

    # initialize output
    panels = Vector{Ring{TF}}(undef, N)

    # populate each panel
    for j = 1:ns

        # angle
        th = linearinterp(0.5, theta[j], theta[j+1])

        r1n = SVector(grid[1,1,j], grid[2,1,j], grid[3,1,j]) # top left
        r2n = SVector(grid[1,1,j+1], grid[2,1,j+1], grid[3,1,j+1]) # top right
        r3n = SVector(grid[1,2,j], grid[2,2,j], grid[3,2,j]) # bottom left
        r4n = SVector(grid[1,2,j+1], grid[2,2,j+1], grid[3,2,j+1]) # bottom right

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

            # bottom left corner of ring vortex
            rbl = linearinterp(0.25, r1n, r3n)

            # bottom right corner of ring vortex
            rbr = linearinterp(0.25, r2n, r4n)

            # control point
            rtop = linearinterp(0.5, r1, r2)
            rbot = linearinterp(0.5, r3, r4)
            rcp = linearinterp(0.75, rl, rr)

            # surface normal
            normal = cross(rcp - rtr, rcp - rtl)
            normal /= norm(normal)

            # panel on trailing edge?
            trailing = false

            ipanel = mirror*right_side*nc*ns + (j-1)*nc + i
            panels[ipanel] = Ring{TF}(rtl, rtr, rbl, rbr, rcp, normal, trailing)
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

        # bottom left corner of ring vortex
        rbl = r3

        # bottom right corner of ring vortex
        rbr = r4

        # control point
        rtop = linearinterp(0.5, r1, r2)
        rbot = linearinterp(0.5, r3, r4)
        rcp = linearinterp(0.75, rl, rr)

        # panel on trailing edge?
        trailing = true

        # surface normal
        normal = cross(rcp - rtr, rcp - rtl)
        normal /= norm(normal)

        ipanel = mirror*right_side*nc*ns + (j-1)*nc + nc
        panels[ipanel] = Ring{TF}(rtl, rtr, rbl, rbr, rcp, normal, trailing)

    end

    # other side
    if mirror
        for j = 1:ns
            for i = 1:nc
                ipanel1 = right_side*nc*ns + (j-1)*nc + i
                ipanel2 = left_side*nc*ns + (ns-j)*nc + i
                panels[ipanel2] = reflect_y(panels[ipanel1])
            end
        end
    end

    return panels
end


"""
    grid_to_horseshoe_vortices(xyz, theta, ns, nc; spacing_s=Cosine(), spacing_c=Uniform(),
        interp_s=(x)->(x,y,xpt)->LinearInterpolation(x, y)(xpt), interp_c=(x,y,xpt)->LinearInterpolation(x, y)(xpt))

Discretize a grid into `ns` spanwise and `nc` chordwise horseshoe vortex panels
according to the spanwise discretization scheme `spacing_s` and chordwise
discretization scheme `spacing_c`.

Each chordwise strip must have the same y and z-coordinates (i.e. the airfoil geometry
cannot be cambered or twisted).  Local angle of attack at each spanwise increment
is specified with the input vector `theta` which has length `size(xyz, 3)`.

The resulting vector of panels can be reshaped to correspond to the original
grid with `reshape(panels, ni-1, nj-1)` (`reshape(panels, ni-1, 2*(nj-1))` if
the geometry is mirrored).
"""
function grid_to_horseshoe_vortices(xyz, theta, ns, nc; mirror=false, spacing_s=Cosine(),
    spacing_c=Uniform(), interp_s=(x,y,xpt)->LinearInterpolation(x,y)(xpt),
    interp_c=(x,y,xpt)->LinearInterpolation(x,y)(xpt))

    TF = promote_type(eltype(xyz), eltype(theta))

    N = nc*ns

    # get spanwise and chordwise spacing
    etas, etabar = spanwise_spacing(ns+1, spacing_s)
    eta_qtr, eta_thrqtr = chordwise_spacing(nc+1, spacing_c)

    # interpolate grid first along `i` then along `j`
    xyz_bound = interpolate_grid(xyz, eta_qtr, interp_c)
    xyz_bound = interpolate_grid(xyz_bound, etas, interp_s)

    xyz_cp = interpolate_grid(xyz, eta_thrqtr, interp_c)
    xyz_cp = interpolate_grid(xyz_cp, etabar, interp_s)

    # check which side we're working with
    right_side = sum(xyz[2,:,:]) > 0
    left_side = !right_side

    # interpolate local angle of attack to new locations
    ds = [norm(xyz[:,1,j] - xyz[:,1,max(1,j-1)]) for j = 1:size(xyz,3)]
    t = cumsum(dc)
    t /= t[end]
    theta = interp_s(t, theta, etabar)

    # initialize output
    panels = Vector{Horseshoe{TF}}(undef, N)

    # populate each panel
    for j = 1:ns
        for i = 1:nc
            rl = SVector(xyz_bound[1,i,j], xyz_bound[2,i,j], xyz_bound[3,i,j])
            rr = SVector(xyz_bound[1,i,j+1], xyz_bound[2,i,j+1], xyz_bound[3,i,j+1])
            rcp = SVector(xyz_cp[1,i,j], xyz_cp[2,i,j], xyz_cp[3,i,j])
            th = theta[j]

            ipanel = mirror*right_side*nc*ns + (j-1)*nc + i
            panels[ipanel] = Horseshoe{TF}(rl, rr, rcp, th)
        end
    end

    # other side
    if mirror
        for j = 1:ns
            for i = 1:nc
                ipanel1 = right_side*nc*ns + (j-1)*nc + i
                ipanel2 = left_side*nc*ns + (ns-j)*nc + i
                panels[ipanel2] = reflect_y(panels[ipanel1])
            end
        end
    end

    return panels
end

"""
    grid_to_vortex_rings(xyz, ns, nc; mirror=false, spacing_s=Cosine(), spacing_c=Uniform(),
        interp_s=(x,y,xpt)->LinearInterpolation(x, y)(xpt),
        interp_c=(x,y,xpt)->LinearInterpolation(x, y)(xpt))

Discretize a grid into `ns` spanwise and `nc` chordwise vortex ring panels
according to the spanwise discretization scheme `spacing_s` and chordwise
discretization scheme `spacing_c`.

# Arguments
 - `xyz`: grid of dimensions (3, i, j) where where `i` corresponds to the
    chordwise direction and `j` corresponds to the spanwise direction.
 - `ns`: number of spanwise panels
 - `nc`: number of chordwise panels
 - `mirror`:  mirror the geometry across the y-axis? defaults to `false`.
 - `spacing_s`: spanwise discretization scheme, defaults to `Cosine()`
 - `spacing_c`: chordwise discretization scheme, defaults to `Uniform()`
 - `interp_s`: spanwise interpolation function, defaults to linear interpolation
 - `interp_c`: chordwise interpolation function, defaults to linear interpolation

The resulting vector of panels can be reshaped into a grid format with
`reshape(panels, nc, ns)` (`reshape(panels, nc, 2*ns)` if the geometry is mirrored).
"""
function grid_to_vortex_rings(xyz, ns, nc; mirror=false, spacing_s=Cosine(),
    spacing_c=Uniform(), interp_s=(x,y,xpt)->LinearInterpolation(x, y)(xpt),
    interp_c=(x,y,xpt)->LinearInterpolation(x, y)(xpt))

    TF = eltype(xyz)

    N = nc*ns

    # get spanwise and chordwise spacing
    etas, etabar = spanwise_spacing(ns+1, spacing_s)
    eta_qtr, eta_thrqtr = chordwise_spacing(nc+1, spacing_c)

    # add trailing edge to grid
    vcat(eta_qtr, 1.0)

    # interpolate grid first along `i` then along `j`
    xyz_bound = interpolate_grid(xyz, eta_qtr, interp_c)
    xyz_bound = interpolate_grid(xyz_bound, etas, interp_s)

    xyz_cp = interpolate_grid(xyz, eta_thrqtr, interp_c)
    xyz_cp = interpolate_grid(xyz_cp, etabar, interp_s)

    # check which side we're working with
    right_side = sum(xyz[2,:,:]) > 0
    left_side = !right_side

    # initialize output
    panels = Vector{Ring{TF}}(undef, N)

    # populate each panel
    for j = 1:ns
        for i = 1:nc
            rtl = SVector(xyz_bound[1,i,j], xyz_bound[2,i,j], xyz_bound[3,i,j])
            rtr = SVector(xyz_bound[1,i,j+1], xyz_bound[2,i,j+1], xyz_bound[3,i,j+1])
            rbl = SVector(xyz_bound[1,i+1,j], xyz_bound[2,i+1,j], xyz_bound[3,i+1,j])
            rbr = SVector(xyz_bound[1,i+1,j+1], xyz_bound[2,i+1,j+1], xyz_bound[3,i+1,j+1])
            rcp = SVector(xyz_cp[1,i,j], xyz_cp[2,i,j], xyz_cp[3,i,j])

            normal = cross(rcp - rtr, rcp - rtl)
            normal /= norm(normal)

            trailing = i == nc

            ipanel = mirror*right_side*nc*ns + (j-1)*nc + i
            panels[ipanel] = Ring{TF}(rtl, rtr, rbl, rbr, rcp, normal, trailing)
        end
    end

    # other side
    if mirror
        for j = 1:ns
            for i = 1:nc
                ipanel1 = right_side*nc*ns + (j-1)*nc + i
                ipanel2 = left_side*nc*ns + (ns-j)*nc + i
                panels[ipanel2] = reflect_y(panels[ipanel1])
            end
        end
    end

    return panels
end

"""
    wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, ns, nc; mirror = false,
        spacing_s=Cosine(), spacing_c=Uniform(), interp_s=(x,y,xpt)->LinearInterpolation(x, y)(xpt))

Discretize a wing into `ns` spanwise and `nc` chordwise horseshoe vortex panels
according to the spanwise discretization scheme `spacing_s` and chordwise
discretization scheme `spacing_c`.

# Arguments
 - `xle`: leading edge x-coordinate of each station
 - `yle`: leading edge y-coordinate of each station
 - `zle`: leading edge z-coordinate of each station
 - `chord`: chord length of each station
 - `theta`: local angle of attack of each station
 - `ns`: number of spanwise panels
 - `nc`: number of chordwise panels
 - `mirror`:  mirror the geometry across the y-axis?, defaults to `false`
 - `spacing_s`: spanwise discretization scheme, defaults to `Cosine()`
 - `spacing_c`: chordwise discretization scheme, defaults to `Uniform()`
 - `interp_s`: interpolation function between spanwise stations, defaults to linear interpolation

The resulting vector of panels can be reshaped into a grid format with
`reshape(panels, nc, ns)` (`reshape(panels, nc, 2*ns)` if the geometry is mirrored).
"""
function wing_to_horseshoe_vortices(xle, yle, zle, chord, theta, ns, nc; mirror = false,
    spacing_s=Cosine(), spacing_c=Uniform(), interp_s=(x,y,xpt)->LinearInterpolation(x, y)(xpt))

    TF = promote_type(eltype(xle), eltype(yle), eltype(zle), eltype(chord), eltype(theta))

    N = (1+mirror)*nc*ns

    # get spanwise and chordwise spacing
    etas, etabar = spanwise_spacing(ns+1, spacing_s)
    eta_qtr, eta_thrqtr = chordwise_spacing(nc+1, spacing_c)

    # check input dimensions
    n = length(chord)
    @assert length(xle) == length(yle) == length(zle) == length(chord) == length(theta)

    # construct interpolation vector
    ds = [sqrt(( xle[j] - xle[max(1, j-1)])^2 +
               ( yle[j] - yle[max(1, j-1)])^2 +
               ( zle[j] - zle[max(1, j-1)])^2 ) for j = 1:n]
    t = cumsum(ds)
    t /= t[end]

    # interpolate properties to new stations
    xle_v = interp_s(t, xle, etas)
    yle_v = interp_s(t, yle, etas)
    zle_v = interp_s(t, zle, etas)

    xle_cp = interp_s(t, xle, etabar)
    yle_cp = interp_s(t, yle, etabar)
    zle_cp = interp_s(t, zle, etabar)

    chord_v = interp_s(t, chord, etas)
    chord_cp = interp_s(t, chord, etabar)

    theta = interp_s(t, theta, etabar)

    # initialize output
    panels = Vector{Horseshoe{TF}}(undef, N)

    # check which side we're working with
    right_side = sum(yle) > 0
    left_side = !right_side

    for j = 1:ns
        for i = 1:nc
            rl = SVector(xle_v[j] + eta_qtr[i]*chord_v[j], yle_v[j], zle_v[j])
            rr = SVector(xle_v[j+1] + eta_qtr[i]*chord_v[j+1], yle_v[j+1], zle_v[j+1])
            rcp = SVector(xle_cp[j]+eta_thrqtr[i]*chord_cp[j], yle_cp[j], zle_cp[j])
            th = theta[j]
            ipanel = mirror*right_side*nc*ns + (j-1)*nc + i
            panels[ipanel] = Horseshoe{TF}(rl, rr, rcp, th)
        end
    end

    # other side
    if mirror
        for j = 1:ns
            for i = 1:nc
                ipanel1 = right_side*nc*ns + (j-1)*nc + i
                ipanel2 = left_side*nc*ns + (ns-j)*nc + i
                panels[ipanel2] = reflect_y(panels[ipanel1])
            end
        end
    end

    return panels
end

"""
    wing_to_vortex_rings(xle, yle, zle, chord, theta, fc, ns, nc; mirror=false,
        spacing_s=Cosine(), spacing_c=Uniform(), interp_s=(x,y,xpt)->LinearInterpolation(x, y)(xpt))

Discretize a wing into `ns` spanwise and `nc` chordwise vortex ring panels
according to the spanwise discretization scheme `spacing_s` and chordwise
discretization scheme `spacing_c`.

# Arguments
 - `xle`: leading edge x-coordinate of each station
 - `yle`: leading edge y-coordinate of each station
 - `zle`: leading edge z-coordinate of each station
 - `chord`: chord length of each station
 - `theta`: local angle of attack at each station
 - `fc`: camber line function y=f(x) (for a normalized airfoil) at each station
 - `ns`: number of spanwise panels
 - `nc`: number of chordwise panels
 - `mirror`:  mirror the geometry across the y-axis?, defaults to `false`
 - `spacing_s`: spanwise discretization scheme, defaults to `Cosine()`
 - `spacing_c`: chordwise discretization scheme, defaults to `Uniform()`
 - `interp_s`: interpolation function between spanwise stations, defaults to linear interpolation
 - `apply_twist`: indicates whether twist should be applied to the geometry in
        addition to the normal vector, defaults to `true`

The resulting vector of panels can be reshaped into a grid format with
`reshape(panels, ni-1, nj-1)` (`reshape(panels, nc, 2*ns)` if the geometry is mirrored).
"""
function wing_to_vortex_rings(xle, yle, zle, chord, theta, fc, ns, nc; mirror=false,
    spacing_s=Cosine(), spacing_c=Uniform(), interp_s=(x,y,xpt)->LinearInterpolation(x, y)(xpt),
    apply_twist=true)

    TF = promote_type(eltype(xle), eltype(yle), eltype(zle), eltype(chord), eltype(theta))

    N = (1+mirror)*nc*ns

    # get spanwise and chordwise spacing
    etas, etabar = spanwise_spacing(ns+1, spacing_s)
    eta_qtr, eta_thrqtr = chordwise_spacing(nc+1, spacing_c)

    # add trailing edge to grid
    eta_qtr = vcat(eta_qtr, 1.0)

    # check input dimensions
    n = length(fc)
    @assert length(xle) == length(yle) == length(zle) == length(chord) == length(theta) == n

    # get bound vortex chordwise locations
    xyz_v = Array{TF, 3}(undef, 3, nc+1, n)
    for j = 1:n
        # rotation matrix for angle of attack
        ct, st = cos(theta[j]), sin(theta[j])
        Rt = @SMatrix [ct 0 st; 0 1 0; -st 0 ct]

        # location of leading edge
        rle = SVector(xle[j], yle[j], zle[j])

        for i = 1:nc+1
            # normalized chordwise location
            xc = eta_qtr[i]
            yc = fc[j](xc)

            # location on airfoil
            r = SVector(xc, yc, 0)

            # scale by chord length
            r *= chord[j]

            # apply twist (unless otherwise specified)
            if apply_twist
                r = Rt * r
            end

            # add leading edge offset
            r += rle

            # store final location
            xyz_v[:,i,j] = r
        end
    end

    # get control point chordwise locations
    xyz_cp = Array{TF, 3}(undef, 3, nc, n)
    for j = 1:n
        # rotation matrix for angle of attack
        ct, st = cos(theta[j]), sin(theta[j])
        Rt = @SMatrix [ct 0 st; 0 1 0; -st 0 ct]

        # location of leading edge
        rle = SVector(xle[j], yle[j], zle[j])

        for i = 1:nc
            # normalized chordwise location
            xc = eta_thrqtr[i]
            yc = fc[j](xc)

            # location on airfoil
            r = SVector(xc, yc, 0)

            # scale by chord length
            r *= chord[j]

            # apply twist (unless otherwise specified)
            if apply_twist
                r = Rt * r
            end

            # add leading edge offset
            r += rle

            # store final location
            xyz_cp[:,i,j] = r
        end
    end

    # now interpolate the grid to match the specified spanwise spacing
    xyz_v = interpolate_grid(xyz_v, etas, interp_s; dir=2)
    xyz_cp = interpolate_grid(xyz_cp, etabar, interp_s; dir=2)

    # interpolate twist to new spanwise stations if not included in the geometry
    if !apply_twist
        ds = [sqrt(( xle[j] - xle[max(1, j-1)])^2 +
                   ( yle[j] - yle[max(1, j-1)])^2 +
                   ( zle[j] - zle[max(1, j-1)])^2 ) for j = 1:n]
        t = cumsum(ds)/sum(ds)
        theta = interp_s(t, theta, etabar)
    end

    # initialize output
    panels = Vector{Ring{TF}}(undef, N)

    # check which side we're working with
    right_side = sum(yle) > 0
    left_side = !right_side

    # populate each panel
    for j = 1:ns
        for i = 1:nc
            rtl = SVector(xyz_v[1,i,j], xyz_v[2,i,j], xyz_v[3,i,j])
            rtr = SVector(xyz_v[1,i,j+1], xyz_v[2,i,j+1], xyz_v[3,i,j+1])
            rbl = SVector(xyz_v[1,i+1,j], xyz_v[2,i+1,j], xyz_v[3,i+1,j])
            rbr = SVector(xyz_v[1,i+1,j+1], xyz_v[2,i+1,j+1], xyz_v[3,i+1,j+1])
            rcp = SVector(xyz_cp[1,i,j], xyz_cp[2,i,j], xyz_cp[3,i,j])

            # estimate normal vector from control point and top bound vortex
            normal = cross(rcp - rtr, rcp - rtl)
            normal /= norm(normal)

            # add twist to normal vector if it wasn't included in the geometry
            if !apply_twist
                ct, st = cos(theta[j]), sin(theta[j])
                Rt = @SMatrix [ct 0 st; 0 1 0; -st 0 ct]
                normal = Rt*normal
            end

            trailing = i == nc

            ipanel = mirror*right_side*nc*ns + (j-1)*nc + i
            panels[ipanel] = Ring{TF}(rtl, rtr, rbl, rbr, rcp, normal, trailing)
        end
    end

    # other side
    if mirror
        for j = 1:ns
            for i = 1:nc
                ipanel1 = right_side*nc*ns + (j-1)*nc + i
                ipanel2 = left_side*nc*ns + (ns-j)*nc + i
                panels[ipanel2] = reflect_y(panels[ipanel1])
            end
        end
    end

    return panels
end
