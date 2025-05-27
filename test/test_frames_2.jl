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

# propeller 1
xle_p1 = [0.0, 0.1]
yle_p1 = [0.0, 0.8]
zle_p1 = [0.0, 0.0]
chord_p1 = [0.2, 0.1]
theta_p1 = [15.0, 3.0] .* pi / 180 # convert to radians
phi_p1 = [0.0, 0.0]
fc_p1 = fill((xc) -> 0, 2) # camberline function for each section
ns_p1 = 6
nc_p1 = 3
spacing_s_p1 = Cosine()
spacing_c_p1 = Uniform()
mirror_p1 = false

# propeller 2
xle_p2 = [0.0, 0.1]
yle_p2 = [0.0, -0.8]
zle_p2 = [0.0, 0.0]
chord_p2 = [0.2, 0.1]
theta_p2 = [15.0, 3.0] .* pi / 180 # convert to radians
phi_p2 = [0.0, 0.0]
fc_p2 = fill((xc) -> 0, 2) # camberline function for each section
ns_p2 = 6
nc_p2 = 3
spacing_s_p2 = Cosine()
spacing_c_p2 = Uniform()
mirror_p2 = false

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
w1grid, w1ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
    mirror=false, fc=fc, spacing_s=spacing_s, spacing_c=spacing_c)
w2grid = VortexLattice.flipy(w1grid)  # mirror
w2ratio = w1ratio

# generate surface panels for horizontal tail
hgrid, hratio = wing_to_grid(xle_h, yle_h, zle_h, chord_h, theta_h, phi_h, ns_h, nc_h;
    mirror=mirror_h, fc=fc_h, spacing_s=spacing_s_h, spacing_c=spacing_c_h)
translate!(hgrid, [4.0, 0.0, 0.0])

# generate surface panels for vertical tail
vgrid, vratio = wing_to_grid(xle_v, yle_v, zle_v, chord_v, theta_v, phi_v, ns_v, nc_v;
    mirror=mirror_v, fc=fc_v, spacing_s=spacing_s_v, spacing_c=spacing_c_v)
translate!(vgrid, [4.0, 0.0, 0.0])

# generate surface panels for propeller 1
p1grid1, p1ratio1 = wing_to_grid(xle_p1, yle_p1, zle_p1, chord_p1, theta_p1, phi_p1, ns_p1, nc_p1;
    mirror=mirror_p1, fc=fc_p1, spacing_s=spacing_s_p1, spacing_c=spacing_c_p1)
p1grid2, p1ratio2 = wing_to_grid(xle_p1, yle_p1, zle_p1, chord_p1, theta_p1, phi_p1, ns_p1, nc_p1;
mirror=mirror_p1, fc=fc_p1, spacing_s=spacing_s_p1, spacing_c=spacing_c_p1)

# opposite chirality for propeller 2
p2grid1 = VortexLattice.flipy(p1grid1)
p2grid2 = VortexLattice.flipy(p1grid1)

translate!(p1grid1, SVector{3}(-chord_p1[1]*0.5, 0.0, 0.0))
translate!(p1grid2, SVector{3}(-chord_p1[1]*0.5, 0.0, 0.0))
R1_b1 = VortexLattice.Rodrigues(SVector{3}(0.0, 1.0, 0.0), pi*0.5)
rotate!(p1grid1, R1_b1)
R1_b2 = VortexLattice.Rodrigues(SVector{3}(1.0,0,0), pi*1.0) * R1_b1
rotate!(p1grid2, R1_b2)
o_p1 = SVector{3}(-0.4, 2.5, 0.5)
translate!(p1grid1, o_p1)
translate!(p1grid2, o_p1)

# generate surface panels for propeller 2
p2ratio1 = p1ratio1
p2ratio2 = p1ratio2
translate!(p2grid1, SVector{3}(-chord_p2[1]*0.5, 0.0, 0.0))
translate!(p2grid2, SVector{3}(-chord_p2[1]*0.5, 0.0, 0.0))
R2_b1 = VortexLattice.Rodrigues(SVector{3}(0.0, 1.0, 0.0), -pi*0.5)
rotate!(p2grid1, R2_b1)
R2_b2 = VortexLattice.Rodrigues(SVector{3}(1.0,0,0), pi*1.0) * R2_b1
rotate!(p2grid2, R2_b2)
o_p2 = SVector{3}(-0.4, -2.5, 0.5)
translate!(p2grid1, o_p2)
translate!(p2grid2, o_p2)

# create system
grids = [w1grid, w2grid, hgrid, vgrid, p1grid1, p1grid2, p2grid1, p2grid2]
ratios = [w1ratio, w2ratio, hratio, vratio, p1ratio1, p1ratio2, p2ratio1, p2ratio2]
surface_id = [1, 2, 3, 4, 5, 6, 7, 8]
system = System(grids; ratios)
symmetric = [false for _ in grids]  # no implied symmetry

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
    ω = 0.0 * 2 * pi,
    R = SMatrix{3,3}(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0),
    name = "vehicle",
    child_index = Int[],
    dependent_index = collect(1:length(system.surfaces))
)

# add wing frames
VortexLattice.add_frame!(frames, "wing1", "vehicle", SVector{3}(0.0, 0.0, -10.0), [1]; 
    v = SVector{3}(0.0, 0.0, 0.25), ω_axis = SVector{3}(-1.0, 0.0, 0.0), ω = -0.025 * 2 * pi,
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

VortexLattice.add_frame!(frames, "wing2", "vehicle", SVector{3}(0.0, 0.0, -10.0), [2]; 
    v = SVector{3}(0.0, 0.0, 0.25), ω_axis = SVector{3}(-1.0, 0.0, 0.0), ω = 0.025 * 2 * pi,
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

# add propeller and blade frames
VortexLattice.add_frame!(frames, "propeller1", "wing1", o_p1, [5,6]; 
    v = SVector{3}(0.0, 0.0, 0.0), ω_axis = SVector{3}(1.0, 0.0, 0.0), ω = -3.0 * 2 * pi,
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

VortexLattice.add_frame!(frames, "propeller2", "wing2", o_p2, [7,8]; 
    v = SVector{3}(0.0, 0.0, 0.0), ω_axis = SVector{3}(1.0, 0.0, 0.0), ω = 3.0 * 2 * pi,
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))

#--- perform maneuver ---#

function maneuver!(frames, system, wake, t)
    if 0.0 <= t < 0.25
        (; x, R, name, parent_index, child_index, dependent_index) = frames[2]
        v = SVector{3}(0.0, 0.0, 0.25)
        ω_axis = SVector{3}(-1.0, 0.0, 0.0)
        ω = -0.025 * 2 * pi
        frames[2] = typeof(frames[2])(x, v, ω_axis, ω, R, name, parent_index, child_index, dependent_index)
        (; x, R, name, parent_index, child_index, dependent_index) = frames[3]
        frames[3] = typeof(frames[3])(x, v, ω_axis, -ω, R, name, parent_index, child_index, dependent_index)
    elseif 0.25 <= t < 0.5
        v = SVector{3}(0.0, 0.0, 0.25) * -4
        ω_axis = SVector{3}(-1.0, 0.0, 0.0)
        ω = -0.025 * 2 * pi * -8
        (; x, R, name, parent_index, child_index, dependent_index) = frames[2]
        frames[2] = typeof(frames[2])(x, v, ω_axis, ω, R, name, parent_index, child_index, dependent_index)
        (; x, R, name, parent_index, child_index, dependent_index) = frames[3]
        frames[3] = typeof(frames[3])(x, v, ω_axis, -ω, R, name, parent_index, child_index, dependent_index)
    elseif 0.5 <= t
        v = SVector{3}(0.0, 0.0, 0.25) * 3 * 0.5
        ω_axis = SVector{3}(-1.0, 0.0, 0.0)
        ω = -0.025 * 2 * pi * 7 * 0.5
        (; x, R, name, parent_index, child_index, dependent_index) = frames[2]
        frames[2] = typeof(frames[2])(x, v, ω_axis, ω, R, name, parent_index, child_index, dependent_index)
        (; x, R, name, parent_index, child_index, dependent_index) = frames[3]
        frames[3] = typeof(frames[3])(x, v, ω_axis, -ω, R, name, parent_index, child_index, dependent_index)
    end
end

function test_loop!(system, frames, nt, tf=1.0)
    write_vtk("test_frames_2_step_0.vtk", system; trailing_vortices=false)
    dt = tf / nt  # time step for the maneuver

    for (i_step, t) in enumerate(range(0.0, stop=tf, length=nt))
        println("\nStep $(i_step) at time $(t) s")
        
        # perform maneuver
        maneuver!(frames, system, system.wakes, t)
        
        # propagate kinematics
        propagate_kinematics!(system, frames, dt)
        
        # save vtk files
        write_vtk("test_frames_2_step_$(i_step).vtk", system; trailing_vortices=false)
    end

    # for i_step in 1:nt
    #     println("\nStep $(i_step) at time $(i_step * dt) s")
    #     propagate_kinematics!(system, frames, dt)
        
    #     # save vtk files
    #     write_vtk("test_frames_2_step_$(i_step).vtk", system; trailing_vortices=false)
    # end

    # # resume motion
    # (;x, v, ω_axis, ω, R, name, parent_index, child_index, dependent_index) = frames[2]
    # frames[2] = typeof(frames[2])(x, -4*v, ω_axis, -8*ω, R, name, parent_index, child_index, dependent_index)
    # (;x, v, ω_axis, ω, R, name, parent_index, child_index, dependent_index) = frames[3]
    # frames[3] = typeof(frames[3])(x, -4*v, ω_axis, -8*ω, R, name, parent_index, child_index, dependent_index)
    
    # # finish maneuver
    # for i_step in nt+1:2*nt
    #     println("\nStep $(i_step) at time $(i_step * dt) s")
    #     propagate_kinematics!(system, frames, dt)
    #     write_vtk("test_frames_2_step_$(i_step).vtk", system; trailing_vortices=false)
    # end

    # # resume motion
    # (;x, v, ω_axis, ω, R, name, parent_index, child_index, dependent_index) = frames[2]
    # frames[2] = typeof(frames[2])(x, -v/4*3*0.5, ω_axis, -ω/8*7*0.5, R, name, parent_index, child_index, dependent_index)
    # (;x, v, ω_axis, ω, R, name, parent_index, child_index, dependent_index) = frames[3]
    # frames[3] = typeof(frames[3])(x, -v/4*3*0.5, ω_axis, -ω/8*7*0.5, R, name, parent_index, child_index, dependent_index)
    
    # # finish maneuver
    # for i_step in 2*nt+1:4*nt
    #     println("\nStep $(i_step) at time $(i_step * dt) s")
    #     propagate_kinematics!(system, frames, dt)
    #     write_vtk("test_frames_2_step_$(i_step).vtk", system; trailing_vortices=false)
    # end
    
end

test_loop!(system, frames, 200)

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