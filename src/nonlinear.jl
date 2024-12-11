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

function SectionProperties(panels_indicies, gammas, area, airfoil, contour)
    α = cl = cd = zeros()
    force = zeros(3)
    return SectionProperties(α, cl, cd, panels_indicies, gammas, area, force, airfoil, contour)
end

function grid_to_sections(grid, airfoils; ratios, contours)
    _, _, surface = grid_to_surface_panels(grid; ratios)
    ns = size(surface, 2)
    sections = Vector{SectionProperties{TF}}(undef, ns)
    if length(airfoils != ns)
        error("Number of airfoils must match number of spanwise panels")
    end

    rows = collect(1:size(surface, 1))
    cols = zeros(Int, size(rows))
    gamma_vec = zeros(Int, size(rows))
    gamma_start = 1
    for i in 1:ns
        cols .= i
        panels = CartesianIndex.(rows, cols)
        gamma_vec .= collect(gamma_start:gamma_start + size(panels, 1) - 1)
        gammas = CartesianIndex.(gamma_vec)
        gamma_start = gamma_vec[end] + 1

        area = 
        sections[i] = SectionProperties(panels, gammas, area, airfoils[i], contour[i])
    end
    return sections
end