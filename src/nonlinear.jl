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
    α = zeros()
    cl = zeros()
    cd = zeros()
    force = zeros(3)
    return SectionProperties(α, cl, cd, panels_indicies, gammas, area, force, airfoil, contour)
end

function grid_to_sections(grid, airfoils; 
                    ratios=zeros(2, size(grid, 2)-1, size(grid, 3)-1) .+ [0.5;0.75], 
                    contours=zeros(1,1))
                    
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
    vel = zeros(2)
    n_hat = zeros(3)
    c_hat = zeros(3)
    v = zeros(3)

    for i in eachindex(system.surfaces)
        surface = system.surfaces[i]
        properties = system.properties[i]
        sections = system.sections[i]
        vx = zeros(size(surface,1))
        vy = zeros(size(surface,1))
        x_c = zeros(size(surface,1))
        for j in eachindex(sections)
            total_chord = 0.0
            section = sections[j]
            c_hat .= surface[section.panels[end]].rbc - surface[section.panels[1]].rtc
            c_hat ./= norm(c_hat)
            n_hat .= cross(c_hat, surface[section.panels[end]].rbr - surface[section.panels[1]].rtl)
            n_hat ./= norm(n_hat)

            for k in eachindex(section.panels)
                total_chord += surface[section.panels[k]].chord
                v .= properties[section.panels[k]].velocity_from_streamwise
                vx[k] = dot(v, n_hat)
                vy[k] = dot(v, c_hat)
            end
            for k in eachindex(x_c)
                x_c[k] = norm(surface[section.panels[1]].rtc - surface[section.panels[k]].rcp) / total_chord
            end
            if length(section.panels) > 1
                vel[1] = FLOWMath.akima(x_c, vx, 0.5)
                vel[2] = FLOWMath.akima(x_c, vy, 0.5)
            else
                vel[1] = vx[1]
                vel[2] = vy[1]
            end
            section.α[1] = atan(vel[1], vel[2]) * (-1)^system.invert_normals[i]
            section.cl .= section.airfoil.clspline(section.α[1])
            section.cd .= section.airfoil.cdspline(section.α[1])
            section.α[1] = rad2deg(section.α[1])
        end
    end
    return system
end

function get_coefficient_distribution(sections)
    α = zeros(length(sections))
    cl = zeros(length(sections))
    cd = zeros(length(sections))
    for i in eachindex(sections)
        α[i] = sections[i].α[1]
        cl[i] = sections[i].cl[1]
        cd[i] = sections[i].cd[1]
    end
    return α, cl, cd
end