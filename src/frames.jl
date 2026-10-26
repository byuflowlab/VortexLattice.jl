#=
The frame tree is SixDOF's (`ReferenceFrame`, `add_frame!`, `propagate_kinematics!`). Each
surface of a `System` is one geometry of the tree: `SixDOF.move` rotates and translates its
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

function propagate_kinematics!(system::System, frames::Vector{<:ReferenceFrame}, dt::Real)
    propagate_kinematics!([SurfaceGrid(system, i) for i in eachindex(system.grids)], frames, dt)

    # update panels
    for isurf = 1:length(system.surfaces)
        update_surface_panels!(system.surfaces[isurf], system.grids[isurf]; ratios = system.ratios[isurf], fcore = (c, Δs) -> system.core_size)
    end
end

function rotate_translate!(system::System, i, origin, Rω, dx)
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

function rotate_translate!(system::System{TF}, i, origin::SVector{3}, Rω::SMatrix{3,3}, dx::SVector{3}) where TF
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

function ReferenceFrame(system::System{TF};
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

#------- kinematic velocity -------#

function kinematic_velocity!(Vcp, Vh, Vv, Vte, surfaces, frames::AbstractVector{ReferenceFrame{TF}}; skip_top_level=false) where TF

    # capture the top level frame if requested
    if skip_top_level
        # begin recursion (skipping the top level frame, which is captured by system.fs)
        for i in frames[1].child_index
            frame = frames[1]
            kinematic_velocity!(Vcp, Vh, Vv, Vte, surfaces, frames, i, frame.x, frame.R)
        end
    else
        kinematic_velocity!(Vcp, Vh, Vv, Vte, surfaces, frames, 1, zero(SVector{3,TF}), SMatrix{3,3,TF,9}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0))
    end

end

function kinematic_velocity!(Vcp, Vh, Vv, Vte, surfaces, frames::AbstractVector{<:ReferenceFrame}, i_frame::Int, dx_parent_to_global, R_parent_to_global)
    # get frame
    frame = frames[i_frame]

    # this frame's origin in global frame
    origin_global = R_parent_to_global * frame.x + dx_parent_to_global # global origin vector

    # this frame's velocity in global frame
    v_global = R_parent_to_global * frame.v # global velocity vector
    ω_global = R_parent_to_global * frame.omega_axis * frame.omega # global angular velocity
    
    # update the kinematic velocity of each dependent surface
    for isurf in frame.dependent_index
        
        # unpack containers
        surface = surfaces[isurf]
        vcp = Vcp[isurf]
        vh = Vh[isurf]  # horizontal velocity
        vv = Vv[isurf]  # vertical velocity
        vte = Vte[isurf]  # trailing edge velocity

        nc, ns = size(surface)

        # velocity at the control points
        for j = 1:ns, i = 1:nc
            # unpack panel
            panel = surface[i,j]

            # control point velocity
            x = controlpoint(panel)
            vcp[i,j] -= v_global + cross(ω_global, (x - origin_global))
        end

        # velocity at the horizontal bound vortices
        for j = 1:ns
            for i = 1:nc
                # unpack panel
                panel = surface[i,j]

                # top center velocity
                x = top_center(panel)

                vh[i, j] -= v_global + cross(ω_global, (x - origin_global))
            end
            # unpack panel
            panel = surface[end,j]

            # bottom center velocity
            x = bottom_center(panel)

            vh[end,j] -= v_global + cross(ω_global, (x - origin_global))
        end

        # velocity at the vertical bound vortices
        for i = 1:nc
            for j = 1:ns
                panel = surface[i,j]
                x = left_center(panel)
                vv[i,j] -= v_global + cross(ω_global, (x - origin_global))
            end
            panel = surface[i,end]
            x = right_center(panel)
            vv[i,end] -= v_global + cross(ω_global, (x - origin_global))
        end

        # velocity at the trailing edge vertices
        for j = 1:ns
            panel = surface[end,j]
            x = bottom_left(panel)
            vte[j] -= v_global + cross(ω_global, (x - origin_global))
        end
        panel = surface[end,end]
        x = bottom_right(panel)
        vte[end] -= v_global + cross(ω_global, (x - origin_global))

    end

    # new dx_parent_to_global
    dx_parent_to_global = origin_global

    # new R_parent_to_global
    R_parent_to_global = R_parent_to_global * frame.R
    
    # propagate to child frames
    for i in frame.child_index
        kinematic_velocity!(Vcp, Vh, Vv, Vte, surfaces, frames, i, dx_parent_to_global, R_parent_to_global)
    end
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
