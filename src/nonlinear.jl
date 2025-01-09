# This file exists to allow for the nonlinear VLM solver that takes into account airfoil polars

struct SectionProperties{TF}
    α::Array{TF,0}
    cl::Array{TF,0}
    cd::Array{TF,0}
    panels::Vector{CartesianIndex{2}}
    gammas::Vector{CartesianIndex{1}}
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
    nc = size(surface, 1)
    sections = Vector{SectionProperties{typeof(surface[1].chord)}}(undef, ns)
    if length(airfoils) != ns
        error("Number of airfoils must match number of spanwise panels")
    end

    rows = collect(1:nc)
    cols = zeros(Int, size(rows))
    gamma_vec = zeros(Int, size(rows))

    gamma_start = 1
    for i in 1:ns
        cols .= i
        panels = CartesianIndex.(rows, cols)
        gamma_vec .= collect(gamma_start:gamma_start + size(panels, 1) - 1)
        gammas = CartesianIndex.(gamma_vec)
        gamma_start = gamma_vec[end] + 1

        chord = 0.0
        span = norm(surface[panels[1]].rtl - surface[panels[1]].rtr)
        for j in 1:nc 
            chord += surface[panels[j]].chord
        end
        area = chord*span

        sections[i] = SectionProperties(panels, gammas, area, airfoils[i], contours[i])
    end
    return sections
end

function redefine_gamma_index!(sections)
    gamma_start = 1
    for i in eachindex(sections)
        for j in eachindex(sections[i])
            sections[i][j].gammas .= CartesianIndex.(collect(gamma_start:gamma_start + size(sections[i][j].panels, 1) - 1))
            gamma_start += size(sections[i][j].panels, 1) - 1
        end
    end
end

function nonlinear_analysis!(system)
    vel = zeros(3)
    n_hat = zeros(3)
    c_hat = zeros(3)
    v = zeros(3)

    for i in eachindex(system.surfaces)
        surface = system.surfaces[i]
        properties = system.properties[i]
        sections = system.sections[i]
        for j in eachindex(sections)
            total_chord = 0.0

        end
    end
end