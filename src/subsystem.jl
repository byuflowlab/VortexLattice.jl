"""
    SubSystem

A struct that represents a subsystem of a larger system

# Fields:
- `indicies::Vector{CartesianIndex}`: The indicies of the grids of the subsystem in the larger system
- `r::SVector{3, TF}`: The position of the subsystem in the larger system
- `R::SMatrix{3, 3, TF}`: The rotation matrix of the subsystem in the larger system
- `velocity::SVector{3, TF}`: The velocity of the objects in the subsystem
- `angular_velocity::SVector{3, TF}`: The angular velocity of the objects in the larger system
"""

struct SubSystem{TF}
    indicies::Vector{Int}
    r::SVector{3, TF}
    R::SMatrix{3, 3, TF}
    velocity::SVector{3, TF}
    angular_velocity::SVector{3, TF}
end