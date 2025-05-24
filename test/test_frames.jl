using VortexLattice
using VortexLattice.StaticArrays

xle = [0.0, 0.4]
yle = [0.0, 7.5]
zle = [0.0, 3.0]
chord = [2.2, 1.8]
theta = [2.0*pi/180, 2.0*pi/180]
phi = [0.0, 0.0]
fc = fill((xc) -> 0, 2) #camberline function for each section

ns = 12
nc = 6
spacing_s = Uniform()
spacing_c = Uniform()
mirror = true
symmetric = false

Sref = 30.0 # reference area
cref = 2.0  # reference chord
bref = 15.0 # reference span
rref = [0.50, 0.0, 0.0] # reference location for rotations/moments (typically the c.g.)
Vinf = 1.0 # reference velocity (magnitude)
ref = Reference(Sref, cref, bref, rref, Vinf)

alpha = 1.0*pi/180
beta = 0.0
Omega = [0.0; 0.0; 0.0]
fs = Freestream(Vinf, alpha, beta, Omega)

# construct surface
grid, ratio = wing_to_grid(xle, yle, zle, chord, theta, phi, ns, nc;
    fc = fc, spacing_s, spacing_c, mirror)

# create vector containing all grids and ratios
grids = [grid]
ratios = [ratio]

# create the system object
system = System(grids; ratios)
steady_analysis!(system, ref, fs; symmetric)
change_convention!(system, system.reference[].r, ForwardRightDown(), BackRightUp())
steady_analysis!(system, ref, fs; symmetric)

# create reference frames
frames = ReferenceFrame(system;
    origin = SVector{3}(0.0, 0.0, -10.0),
    v = SVector{3}(0.0, 0.0, 0.0),
    ω_axis = SVector{3}(0.0, 1.0, 0.0),
    ω = 1.0 * 2 * pi,
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    name = "vehicle",
    child_index = Int[],
    dependent_index = collect(1:length(system.surfaces))
)

#--- perform maneuver ---#
function test_loop!(system, frames, dt=0.1)
    write_vtk("test_frames_step_0.vtk", system)
    @show frames[1].R
    println("\nBeginning test loop...\n")
    for i_step in 1:10
        println("\tStep $i_step")
        propagate_kinematics!(system, frames, dt)
        @show frames[1].R

        # save vtk files
        write_vtk("test_frames_step_$(i_step).vtk", system)
    end
end

test_loop!(system, frames, 0.1)

# Define a vector to rotate
v = SVector{3}(1.0, 0.0, 0.0)

# Define the rotation axis and angle
axis = SVector{3}(0.0, 0.0, 1.0)  # Rotate around y-axis
angle = pi / 4  # 45 degrees

# Perform the rotation
R = VortexLattice.Rodrigues(axis, angle)
rotated_v = R * v

@show isapprox(rotated_v, SVector{3}(sqrt(2)/2, sqrt(2)/2, 0.0), atol=1e-6)
