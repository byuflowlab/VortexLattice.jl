#=
The frame tree is SixDOF's (`ReferenceFrame`, `add_frame!`, `propagate_kinematics!`). Each
surface of a `VortexLatticeSystem` is one geometry of the tree: `SixDOF.move` rotates and translates its
grid in place.
=#

struct SurfaceGrid{S}
    system::S
    index::Int
end

function SixDOF.move(s::SurfaceGrid, origin, R, dx)
    rotate_translate!(s.system, s.index, origin, R, dx)
    return s
end

struct ForwardRightDown end

struct BackRightUp end

function propagate_kinematics!(system::VortexLatticeSystem, frames::Vector{<:ReferenceFrame}, dt::Real)
    propagate_kinematics!([SurfaceGrid(system, i) for i in eachindex(system.grids)], frames, dt)

    # update panels
    for isurf = 1:length(system.surfaces)
        update_surface_panels!(system.surfaces[isurf], system.grids[isurf]; ratios = system.ratios[isurf], fcore = (c, Δs) -> system.core_size)
    end
end

function rotate_translate!(system::VortexLatticeSystem, i, origin, Rω, dx)
    # get grid
    grid = system.grids[i]

    # grid size
    _, nc, ns = size(grid)

    # rotate/translate
    for i in 1:ns
        for j in 1:nc
            # relative to origin
            for k in 1:3
                grid[k,j,i] -= origin[k]
            end

            # rotate
            grid[:,j,i] .= Rω * SVector{3}(grid[1,j,i], grid[2,j,i], grid[3,j,i])

            # translate and shift origin back
            for k in 1:3
                grid[k,j,i] += dx[k] + origin[k]
            end
        end
    end
end

function rotate_translate!(system::VortexLatticeSystem{TF}, i, origin::SVector{3}, Rω::SMatrix{3,3}, dx::SVector{3}) where TF
    # get grid
    grid::Array{TF,3} = system.grids[i]

    # grid size
    _, nc, ns = size(grid)

    # rotate/translate
    for i in 1:ns
        for j in 1:nc
            # relative to origin
            for k in 1:3
                grid[k,j,i] -= origin[k]
            end

            # rotate
            grid[:,j,i] .= Rω * SVector{3}(grid[1,j,i], grid[2,j,i], grid[3,j,i])

            # translate and shift origin back
            for k in 1:3
                grid[k,j,i] += dx[k] + origin[k]
            end
        end
    end
end

function ReferenceFrame(system::VortexLatticeSystem{TF};
        # vvv all in global frame vvv
        origin = system.reference.rref,
        v = zero(SVector{3,TF}),
        omega_axis = SVector{3,TF}(0.0, 1.0, 0.0),
        omega = zero(TF),
        R = SMatrix{3,3}(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0),
        Rp2g = SMatrix{3,3}(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0), # rotation from parent frame to global frame
        name = "vehicle",  # name of this frame
        child_index = Int[],  # indices of child frames
        dependent_index = collect(1:length(system.surfaces))  # indices of dependent surfaces
        # ^^^ all in global frame ^^^
    ) where TF
    return ReferenceFrame(TF, length(system.surfaces); origin, v, omega_axis, omega, R, Rp2g,
                          name, child_index, dependent_index)
end

"""
    body_freestream(freestream, frames, ref, nsurf) -> BodyFreestream

The relative wind and rotation rates that the surfaces see in a steady analysis.
Every surface moves with one rigid motion (`SixDOF.frame_motion`). The wind is
taken at the reference point `ref.r`. The rates are the body's angular velocity
in the surface axes (back-right-up) in VortexLattice's sign convention, where
positive p rolls the right wing down, positive q pitches nose up and positive r
yaws nose right: p = -omega_x, q = omega_y, r = -omega_z.
"""
function body_freestream(freestream::Freestream, frames, ref, nsurf)
    origin, v, omega = SixDOF.frame_motion(frames, 1)
    vref = v + cross(omega, ref.r - origin)
    for i in 2:nsurf
        o, vi, wi = SixDOF.frame_motion(frames, i)
        scale = max(norm(vref), norm(omega) * norm(ref.r - origin), one(eltype(vref)))
        (isapprox(vi + cross(wi, ref.r - o), vref; atol = sqrt(eps(scale)) * scale) &&
         isapprox(wi, omega; atol = sqrt(eps(scale)) * scale)) || throw(ArgumentError(
            "a steady analysis needs every surface on one rigid motion; surface $i moves relative to surface 1"))
    end
    return velocity_to_freestream(freestream.velocity - vref, SVector(-omega[1], omega[2], -omega[3]))
end

"""
    flight_condition(system, ref, Vinf; alpha=0, beta=0, rates=zeros(3), rho=1.225, mu=1.81e-5, speed_of_sound=343.0) -> (freestream, frames)

The core's `Freestream` and a one-frame tree for a vehicle flying at speed `Vinf`
with angle of attack `alpha`, sideslip `beta` (radians) and body rates
`rates = (p, q, r)` about `ref.r`, for users who think in angles. Pass both to
[`steady_analysis!`](@ref).
"""
function flight_condition(system::VortexLatticeSystem{TF}, ref::Reference, Vinf; alpha = zero(TF), beta = zero(TF),
        rates = zero(SVector{3,TF}), rho = 1.225, mu = 1.81e-5, speed_of_sound = 343.0) where TF
    fs = BodyFreestream(Vinf, alpha, beta, SVector{3}(rates))
    omega = SVector(-rates[1], rates[2], -rates[3])
    w = norm(omega)
    frames = ReferenceFrame(system; origin = SVector{3,TF}(ref.r),
        omega_axis = iszero(w) ? SVector{3,TF}(0, 0, 1) : omega / w, omega = w)
    return Freestream(freestream_velocity(fs); rho, mu, speed_of_sound), frames
end

function change_convention!(system, origin, to::ForwardRightDown, from::BackRightUp)

    # 180 degrees about the y axis
    R180_y = SMatrix{3,3}(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0)

    for i_surf in 1:length(system.surfaces)
        rotate_translate!(system, i_surf, origin, R180_y, SVector{3}(0.0, 0.0, 0.0))
        update_surface_panels!(system.surfaces[i_surf], system.grids[i_surf]; ratios = system.ratios[i_surf], fcore = (c, Δs) -> system.core_size)
    end

    # freestream direction
    x, y, z = system.xhat[]
    system.xhat[] = SVector{3}(-x, y, -z)

    return system
end

function change_convention!(system, to::BackRightUp, from::ForwardRightDown)
    # 180 degrees about the y axis
    R180_y = SMatrix{3,3}(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0)

    for i_surf in 1:length(system.surfaces)
        rotate_translate!(system, i_surf, origin, R180_y, SVector{3}(0.0, 0.0, 0.0))
        update_surface_panels!(system.surfaces[i_surf], system.grids[i_surf]; ratios = system.ratios[i_surf], fcore = (c, Δs) -> system.core_size)
    end

    # freestream direction
    x, y, z = system.xhat[]
    system.xhat[] = SVector{3}(-x, y, -z)

    return system
end

function change_convention!(system, to::ForwardRightDown, from::ForwardRightDown)
    return system
end

function change_convention!(system, to::BackRightUp, from::BackRightUp)
    return system
end
