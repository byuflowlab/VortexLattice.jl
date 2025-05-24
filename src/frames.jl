struct ReferenceFrame{TF}
	x::SVector{3,TF}             # origin in parent frame
	v::SVector{3,TF}             # velocity in parent frame
	ω_axis::SVector{3,TF}        # axis of rotation in parent frame
	ω::TF                        # angular velocity in parent frame
	R::SMatrix{3,3,TF,9}         # basis vectors expressed in parent frame
	name::String                 # name of this frame
	parent_index::Int            # index of parent frame
	child_index::Vector{Int}  # child reference frames
	dependent_index::Vector{Int} # VortexLattice surfaces
end

struct ForwardRightDown end

struct BackRightUp end

function propagate_kinematics!(system::System, frames::Vector{<:ReferenceFrame}, dt::Real)
    
    # global rotation matrix from parent to global frame
    R_parent_to_global = SMatrix{3,3,Float64}(1.0,0,0,0,1.0,0,0,0,1.0)

    # begin recursion
    propagate_kinematics!(system, 1, frames, R_parent_to_global, dt)

    # update panels
    for isurf = 1:length(system.surfaces)
        update_surface_panels!(system.surfaces[isurf], system.grids[isurf]; ratios = system.ratios[isurf])
    end
end

function propagate_kinematics!(system::System, i_frame::Int, frames::Vector{<:ReferenceFrame}, R_parent_to_global::SMatrix, dt::Real)
    # get frame
    frame = frames[i_frame]

    # origin in global frame
    origin_global = R_parent_to_global * frame.x  # global origin vector

    # differential translation
    dx = frame.v * dt  # translation in parent frame
    dx_global = R_parent_to_global * dx  # global translation vector

    # differential rotation
    dω = frame.ω * dt  # angular displacement in parent frame
    Rω = Rodrigues(frame.ω_axis, dω) # rotation matrix in parent frame
    Rω_global = R_parent_to_global * Rω # global frame

    # rotate and translate dependent surfaces
    for i in frame.dependent_index
        rotate_translate!(system, i, origin_global, Rω_global, dx_global)
    end
    
    # Update the frame
    x_new = frame.x + dx
    R_new = Rω * frame.R
    frames[i_frame] = ReferenceFrame(x_new, frame.v, frame.ω_axis, frame.ω, R_new, frame.name, frame.parent_index, frame.child_index, frame.dependent_index)    

    # new R_parent_to_global
    R_parent_to_global = R_parent_to_global * R_new

    # Recursively propagate to child frames
    for i in frame.child_index
        propagate_kinematics!(system, i, frames, R_parent_to_global, dt)
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

function Rodrigues(axis, angle::TF) where TF
    s, c = sincos(-angle)
    t = 1 - c
    x, y, z = axis

    return SMatrix{3,3,TF,9}(
        t*x*x + c,     t*x*y - s*z, t*x*z + s*y,
        t*x*y + s*z,   t*y*y + c,   t*y*z - s*x,
        t*x*z - s*y,   t*y*z + s*x, t*z*z + c
    )
end

function ReferenceFrame(system::System{TF}; 
        # vvv all in global frame vvv
        origin = system.reference.rref,
        v = zero(SVector{3,TF}),
        ω_axis = SVector{3,TF}(0.0, 1.0, 0.0),
        ω = zero(TF),
        R = SMatrix{3,3}(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0),
        name = "vehicle",  # name of this frame
        child_index = Int[],  # indices of child frames
        dependent_index = collect(1:length(system.surfaces))  # indices of dependent surfaces
        # ^^^ all in global frame ^^^
    ) where TF
    parent_index = -1  # no parent frame
    vehicle_frame = ReferenceFrame{TF}(
        origin, v, ω_axis, ω, R,
        name, parent_index, child_index, dependent_index
    )
    frames = [vehicle_frame]
    return frames
end

function add_frame!(frames::Vector{ReferenceFrame{TF}}, parent_index::Int, name::String, origin, surface_indices::Vector{Int};
    v = zero(SVector{3,TF}),  # velocity in parent frame
    ω_axis = SVector{3,TF}(0.0, 1.0, 0.0),  # axis of rotation in parent frame
    ω = zero(TF),  # angular velocity in parent frame
    R = SMatrix{3,3}(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)  # basis vectors expressed in parent frame
) where TF

    # create new frame
    new_frame = ReferenceFrame{TF}(
        origin, v, ω_axis, ω, R,
        name, parent_index, Int[], surface_indices
    )
    push!(frames, new_frame)

    # inform parents
    push!(frames[parent_index].child_index, length(frames))
end

function add_frame!(frames::Vector{<:ReferenceFrame}, parent_name::String, name::String, origin, surface_indices::Vector{Int}; optargs...)
    for (i, frame) in enumerate(frames)
        if frame.name == parent_name
            return add_frame!(frames, i, name, origin, surface_indices; optargs...)
        end
    end
end

function change_convention!(system, origin, to::ForwardRightDown, from::BackRightUp)

    # 180 degrees about the y axis
    R180_y = SMatrix{3,3}(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0)

    for i_surf in 1:length(system.surfaces)
        rotate_translate!(system, i_surf, origin, R180_y, SVector{3}(0.0, 0.0, 0.0))
        update_surface_panels!(system.surfaces[i_surf], system.grids[i_surf]; ratios = system.ratios[i_surf])
    end

    return system
end

function change_convention!(system, to::BackRightUp, from::ForwardRightDown)
    # 180 degrees about the y axis
    R180_y = SMatrix{3,3}(-1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0)

    for i_surf in 1:length(system.surfaces)
        rotate_translate!(system, i_surf, origin, R180_y, SVector{3}(0.0, 0.0, 0.0))
        update_surface_panels!(system.surfaces[i_surf], system.grids[i_surf]; ratios = system.ratios[i_surf])
    end

    return system
end

function change_convention!(system, to::ForwardRightDown, from::ForwardRightDown)
    return system
end

function change_convention!(system, to::BackRightUp, from::BackRightUp)
    return system
end
