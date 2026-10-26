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

function Freestream(frame::ReferenceFrame{TF}, ref::Reference, vinf_ext) where TF
    # equivalent freestream about ref.r
    dx = ref.r - frame.x  # vector from frame origin to reference point
    ω = frame.omega_axis * frame.omega
    V = vinf_ext - frame.v - ω × dx  # freestream velocity in South-East-Up convention

    # create freestream object
    Omega = frame.omega_axis * frame.omega
    fs = velocity_to_freestream(V, -Omega)

    return fs
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
