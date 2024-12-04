# This file exists to allow for the nonlinear VLM solver that takes into account airfoil polars

struct SectionProperties{TF}
    α::Array{TF,0}
    cl::Array{TF,0}
    cd::Array{TF,0}
    panels::Vector{CartesianIndex}
    gammas::Vector{CartesianIndex}
    area::TF
    force::Vector{TF}
    airfoil::CCBlade.AlphaAF{TF, String, Akima{Vector{TF}, Vector{TF}, TF}}
    contour::Matrix{Float64}
end