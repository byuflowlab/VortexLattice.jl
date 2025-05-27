using VortexLattice
using VortexLattice.StaticArrays

# wing
xle = [0.0, 0.2]
yle = [0.0, 5.0]
zle = [0.0, 1.0]
chord = [1.0, 0.6]
theta = [0.0, 0.0]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) # camberline function for each section
ns = 12
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()
mirror = true

# horizontal stabilizer
xle_h = [0.0, 0.14]
yle_h = [0.0, 1.25]
zle_h = [0.0, 0.0]
chord_h = [0.7, 0.42]
theta_h = [0.0, 0.0]
phi_h = [0.0, 0.0]
fc_h = fill((xc) -> 0, 2) # camberline function for each section
ns_h = 6
nc_h = 3
spacing_s_h = Uniform()
spacing_c_h = Uniform()
mirror_h = true

# vertical stabilizer
xle_v = [0.0, 0.14]
yle_v = [0.0, 0.0]
zle_v = [0.0, 1.0]
chord_v = [0.7, 0.42]
theta_v = [0.0, 0.0]
phi_v = [0.0, 0.0]
fc_v = fill((xc) -> 0, 2) # camberline function for each section
ns_v = 5
nc_v = 3
spacing_s_v = Uniform()
spacing_c_v = Uniform()
mirror_v = false

# propeller
xle_p = [0.0, 0.1]
yle_p = [0.0, 0.8]
zle_p = [0.0, 0.0]
chord_p = [0.2, 0.1]
theta_p = [5.0, 1.0] .* pi / 180 # convert to radians
phi_p = [0.0, 0.0]
fc_p = fill((xc) -> 0, 2) # camberline function for each section
ns_p = 6
nc_p = 3
spacing_s_p = Cosine()
spacing_c_p = Uniform()
mirror_p = false

# vehicle reference
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

# generate surface panels for wing
wgrid, wratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=mirror, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)

# generate surface panels for horizontal tail
hgrid, hratio = wing_to_grid(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
    mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
translate!(hgrid, [4.0, 0.0, 0.0])

# generate surface panels for vertical tail
vgrid, vratio = wing_to_grid(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
    mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
translate!(vgrid, [4.0, 0.0, 0.0])

# generate surface panels for propeller
pgrid1, pratio1 = wing_to_grid(xle_p, yle_p, zle_p, chord_p, theta_p, phi_p, ns_p, nc_p;
    mirror=mirror_p, fc=fc_p, spacing_s=spacing_s_p, spacing_c=spacing_c_p)
pgrid2, pratio2 = wing_to_grid(xle_p, yle_p, zle_p, chord_p, theta_p, phi_p, ns_p, nc_p;
    mirror=mirror_p, fc=fc_p, spacing_s=spacing_s_p, spacing_c=spacing_c_p)
translate!(pgrid1, SVector{3}(-chord_p[1]*0.5, 0.0, 0.0))
translate!(pgrid2, SVector{3}(-chord_p[1]*0.5, 0.0, 0.0))
R_p1 = VortexLattice.Rodrigues(SVector{3}(0.0, 1.0, 0.0), pi*0.5)
rotate!(pgrid1, R_p1)
R_p2 = VortexLattice.Rodrigues(SVector{3}(1.0,0,0), pi*1.0) * R_p1
rotate!(pgrid2, R_p2)
o_p = SVector{3}(-0.4, 0.0, 0.0)
translate!(pgrid1, o_p)
translate!(pgrid2, o_p)

# create system
grids = [wgrid, hgrid, vgrid, pgrid1, pgrid2]
ratios = [wratio, hratio, vratio, pratio1, pratio2]
surface_id = [1, 2, 3, 4, 5]
system = System(grids; ratios)
symmetric = [false, false, false, false, false]  # no implied symmetry

# force generation of surface panels
steady_analysis!(system, ref, fs; symmetric, surface_id=surface_id)

# change convention to match flight dynamics
change_convention!(system, system.reference[].r, ForwardRightDown(), BackRightUp())
steady_analysis!(system, ref, fs; symmetric)

# create reference frames
frames = ReferenceFrame(system;
    origin = SVector{3}(0.0, 0.0, -10.0),
    v = SVector{3}(0.0, 0.0, 0.0),
    ω_axis = SVector{3}(0.0, 1.0, 0.0),
    ω = 1.0 * 2 * pi,
    R = SMatrix{3,3}(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0),
    name = "vehicle",
    child_index = Int[],
    dependent_index = collect(1:length(system.surfaces))
)

# add propeller and blade frames
VortexLattice.add_frame!(frames, "propeller1", "vehicle", SVector{3}(0.0, 0.0, -10.0) + o_p, [4,5]; 
    v = SVector{3}(0.0, 0.0, 0.0), ω_axis = SVector{3}(1.0, 0.0, 0.0), ω = -5.0 * 2 * pi,
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

#--- perform maneuver ---#

function test_loop!(system, frames, dt=0.1)
    write_vtk("test_frames_step_0.vtk", system; trailing_vortices=false)
    for i_step in 1:Int(round(1.0 / dt))
        println("\nStep $(i_step) at time $(i_step * dt) s")
        propagate_kinematics!(system, frames, dt)
        
        # save vtk files
        write_vtk("test_frames_step_$(i_step).vtk", system; trailing_vortices=false)
    end
end

test_loop!(system, frames, 0.01)

#=
# Define a vector to rotate
v = SVector{3}(1.0, 0.0, 0.0)

# Define the rotation axis and angle
axis = SVector{3}(0.0, 0.0, 1.0)  # Rotate around y-axis
angle = pi / 4  # 45 degrees

# Perform the rotation
R = VortexLattice.Rodrigues(axis, angle)
rotated_v = R * v

@show isapprox(rotated_v, SVector{3}(sqrt(2)/2, sqrt(2)/2, 0.0), atol=1e-6)
=#