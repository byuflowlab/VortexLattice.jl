export linearsections

abstract type AbstractSpacing end
struct Uniform <: AbstractSpacing end
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
0 and 1 according to the discretization strategy in `spacing_type`
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
    spanwise_spacing(n, spacing_type::AbstractSpacing; xc_cp=0.75)

Distribute `n-1` vortex x-locations and `n-1` control point x-locations on the
interval between 0 and 1 according to the discretization strategy in `spacing_type`.
The control point location can be customized by setting `xc_cp` to either a scalar
or a vector of length `n-1`
"""
function chordwise_spacing(n, ::Uniform; xc_cp=0.75)

    xc_v = 0.25

    eta = range(0, 1.0, length=n)

    eta_qtr = linearinterp(xc_v, eta[1:end-1], eta[2:end])
    eta_thrqtr = linearinterp(xc_cp, eta[1:end-1], eta[2:end])

    return eta_qtr, eta_thrqtr
end

function chordwise_spacing(n, ::Cosine; xc_cp=0.75)

    xc_v = 0.25

    theta = range(0, pi, length=n)
    eta = (1.0 .- cos.(theta))/2.0

    eta_qtr = linearinterp(xc_v, eta[1:end-1], eta[2:end])
    eta_thrqtr = linearinterp(xc_cp, eta[1:end-1], eta[2:end])

    return eta_qtr, eta_thrqtr
end

"""
    create_grid(r1, r2, r3, r4, ns, spacing_s, nc, spacing_c, thetaL, thetaR)

Discretizes a linear wing segment into `ns` spanwise panels and `nc` chordwise
panels according to the spanwise discretization scheme `spacing_s` and chordwise
discretization scheme `spacing_c`.  The points `r1`, `r2`, `r3`, `r4` are defined
as depicted:

1 ------- 2
|         |
|         |
3-------- 4

The command `reshape(panels, nc, ns)` may be used to reshape the output to match
the shape of the original grid.
"""
function create_grid(r1, r2, r3, r4, ns, spacing_s, nc, spacing_c, thetaL, thetaR)

    TF = promote_type(eltype(r1), eltype(r2), eltype(r3), eltype(r4),
        eltype(thetaL), eltype(thetaR))

    etas, etabar = spanwise_spacing(ns+1, spacing_s)
    eta_qtr, eta_thrqtr = chordwise_spacing(nc+1, spacing_c)

    panels = Vector{Panel{TF}}(undef, ns*nc)
    for j = 1:ns
        for i = 1:nc
            rtop = linearinterp(etas[j], r1, r2)
            rbot = linearinterp(etas[j], r3, r4)
            rl = linearinterp(eta_qtr[i], rtop, rbot)

            rtop = linearinterp(etas[j+1], r1, r2)
            rbot = linearinterp(etas[j+1], r3, r4)
            rr = linearinterp(eta_qtr[i], rtop, rbot)

            rtop = linearinterp(etabar[j], r1, r2)
            rbot = linearinterp(etabar[j], r3, r4)
            rcp = linearinterp(eta_thrqtr[i], rtop, rbot)

            theta = linearinterp(etabar[j], thetaL, thetaR)

            panels[(j-1)*nc + i] = Panel(rl, rr, rcp, theta)

        end
    end

    return panels
end

"""
    linear_sections(xle, yle, zle, chord, theta, ns, spacing_s,
    nc, spacing_c, mirror)

Discretizes a wing into multiple linear sections.  Each linear section is defined
using `create_grid`.

# Arguments:
 - `xle`: x-coordinate of the leading edge of each station
 - `yle`: y-coordinate of the leading edge of each station
 - `zle`: z-coordinate of the leading edge of each station
 - `chord`: chord length of each station
 - `theta`: angle of attack of each station
 - `ns`: number of spanwise panels in each section
 - `spacing_s`: spanwise discretization scheme for each section
 - `nc`: number of chordwise panels in each section
 - `spacing_c`: chordwise discretization scheme for each section
 - `mirror`: flag indicating whether the geometry should be mirrored across the y-axis
"""
function linear_sections(xle, yle, zle, chord, theta, ns, spacing_s,
    nc, spacing_c, mirror)

    TF = promote_type(eltype(xle), eltype(yle), eltype(zle), eltype(chord), eltype(theta))

    N = (1+mirror)*sum(nc)*sum(ns)

    panels = Vector{Panel{TF}}(undef, N)
    ipanel = 0

    n = length(ns)
    for i = 1:n
        r1 = SVector(xle[i], yle[i], zle[i])
        r2 = SVector(xle[i+1], yle[i+1], zle[i+1])
        r3 = r1 + SVector(chord[i], 0, 0)
        r4 = r2 + SVector(chord[i+1], 0, 0)

        panels[ipanel+1:ipanel+nc[i]*ns[i]] = create_grid(r1, r2, r3, r4, ns[i], spacing_s[i], nc[i], spacing_c[i], theta[i], theta[i+1])
        ipanel += nc[i]*ns[i]
    end

    if mirror

        for i = 1:n
            r1 = SVector(xle[i+1], -yle[i+1], zle[i+1])
            r2 = SVector(xle[i], -yle[i], zle[i])
            r3 = r1 + SVector(chord[i+1], 0, 0)
            r4 = r2 + SVector(chord[i], 0, 0)

            panels[ipanel+1:ipanel+nc[i]*ns[i]] = create_grid(r1, r2, r3, r4, ns[i], spacing_s[i], nc[i], spacing_c[i], theta[i], theta[i+1])
            ipanel += nc[i]*ns[i]
        end

    end

    return panels
end

"""
    simple_wing(b, AR, λ, Λ, ϕ, θr, θt, npanels, mirror, spacing)

Defines a simple linear wing section with one chordwise panel

# Arguments:
 - `b`: span length
 - `AR`: aspect ratio
 - `λ`: taper ratio
 - `Λ`: sweep (radians)
 - `ϕ`: dihedral (radians)
 - `θr`: root twist
 - `θt`: tip twist
 - `mirror`: flag indicating whether the geometry should be mirrored across the y-axis
 - `spacing`: spanwise discretization scheme
"""
function simplewing(b, AR, λ, Λ, ϕ, θr, θt, npanels, duplicate, spacing)

    # geometry parsing
    S = b^2/AR
    cr = 2*S/(b*(1 + λ))
    ct = cr*λ

    xle = [0.0; cr/4.0 + b/2.0*tan(Λ) - ct/4.0]
    yle = [0.0; b/2*cos(ϕ)]
    zle = [0.0; b/2*sin(ϕ)]
    chord = [cr; ct]
    theta = [θr; θt]

    return linearsections(xle, yle, zle, chord, theta, [npanels], [spacing], [1], ["u"], duplicate)

end

"""
    surface_panels(grid::AbstractArray{3,TF}; cp_xc=0.75, cp_yb=0.5)

Constructs a set of vortex lattice panels given a grid with dimensions (3, i, j).
The i-direction is assumed to be in the chordwise/streamwise direction.  The resulting
vector of panels can be reshaped to correspond to the original grid with
`reshape(panels, ni-1, nj-1)`.

# Arguments:
 - `grid`: Grid of shape (3, i, j) containing the corners of each panel
 - `cp_xc`: The normalized chordwise location of each control point provided as
        either a single value or as a matrix of size (ni-1, nj-1)
 - `cp_yb`: The normalized spanwise location of each control point provided as
        either a single value or as a matrix of size (ni-1, nj-1).
 - `twist`: (optional) the twist of each section specified as a matrix of size
        (ni-1, nj-1).  By default the twist will be estimated from the geometry.
"""
function surface_panels(grid::AbstractArray{TF, 3}; cp_xc=0.75, cp_yb=0.5,
    twist=nothing) where TF

    ni = size(grid, 2)
    nj = size(grid, 3)
    N = (ni-1)*(nj-1)

    # process optional arguments
    if cp_xc <: Number
        cp_xc = fill(cp_xc, ni, nj)
    end
    if cp_yb <: Number
        cp_yb = fill(cp_yb, ni, nj)
    end

    # initialize output
    panels = Vector{Panel{TF}}(undef, N)

    # populate each panel
    ipanel = 0
    for j = 2:nj
        for i = 2:ni
            # left vortex point at 1/4 chord
            rlx = (1-0.25)*grid[1, i-1, j-1] + 0.25*grid[1, i, j-1]
            rly = (1-0.25)*grid[2, i-1, j-1] + 0.25*grid[2, i, j-1]
            rlz = (1-0.25)*grid[3, i-1, j-1] + 0.25*grid[3, i, j-1]
            rl = SVector(rlx, rly, rlz)
            # right vortex point at 1/4 chord
            rrx = (1-0.25)*grid[1, i-1, j] + 0.25*grid[1, i, j]
            rry = (1-0.25)*grid[2, i-1, j] + 0.25*grid[2, i, j]
            rrz = (1-0.25)*grid[3, i-1, j] + 0.25*grid[3, i, j]
            rr = SVector(rrx, rry, rrz)
            # control point at 3/4 chord
            rcpx = (1-cp_xc[i,j])*((1-cp_yb[i,j])*grid[1, i-1, j-1] + cp_yb[i,j]*grid[1, i-1, j]) +
                cp_xc[i,j]*((1-cp_yb[i,j])*grid[1, i, j-1] + cp_yb[i,j]*grid[1, i, j])
            rcpy = (1-cp_xc[i,j])*((1-cp_yb[i,j])*grid[2, i-1, j-1] + cp_yb[i,j]*grid[2, i-1, j]) +
                cp_xc[i,j]*((1-cp_yb[i,j])*grid[2, i, j-1] + cp_yb[i,j]*grid[2, i, j])
            rcpz = (1-cp_xc[i,j])*((1-cp_yb[i,j])*grid[3, i-1, j-1] + cp_yb[i,j]*grid[3, i-1, j]) +
                cp_xc[i,j]*((1-cp_yb[i,j])*grid[3, i, j-1] + 0.5*grid[3, i, j])
            rcp = SVector(rcpx, rcpy, rcpz)
            # twist angle derived from the geometry at control point
            if isnothing(twist)
                dx = ((1-cp_yb[i,j])*grid[1, i, j-1] + cp_yb[i,j]*grid[1, i, j]) -
                    ((1-cp_yb[i,j])*grid[1, i-1, j-1] + cp_yb[i,j]*grid[1, i-1, j])
                dz = ((1-cp_yb[i,j])*grid[3, i-1, j-1] + cp_yb[i,j]*grid[3, i-1, j]) -
                    ((1-cp_yb[i,j])*grid[3, i, j-1] + cp_yb[i,j]*grid[3, i, j])
                theta = dz/dx
            else
                theta = twist[i,j]
            end
            # add panel to array
            ipanel += 1
            panels[ipanel] = Panel{TF}(rl, rr, rcp, theta)
        end
    end

    return panels
end
