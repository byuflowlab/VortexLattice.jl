# This file exists to allow for the nonlinear VLM solver that takes into account airfoil polars

"""
    SectionProperties
"""
struct SectionProperties{TF}
    α::Array{TF,0}
    cl::Array{TF,0}
    cd::Array{TF,0}
    Γs::Array{TF,1}
    λ::Vector{TF}
    panels::Vector{CartesianIndex{2}}
    gammas::Vector{CartesianIndex{1}}
    area::TF
    force::Array{SVector{3,TF},0} # replace with Array{SVector{3,TF},0} if you want to use StaticArrays
    c_hat::Array{SVector{3,TF},0} # unit vector in the direction of the chord
    n_hat::Array{SVector{3,TF},0} # unit vector normal to the section
    airfoil::CCBlade.AlphaAF{TF, String, Akima{Vector{TF}, Vector{TF}, TF}}
    contour::Matrix{Float64}
end

function SectionProperties(panels_indicies, gammas, area, airfoil, contour)
    α = zeros()
    cl = zeros()
    cd = zeros()
    force = fill(SVector{3, eltype(α)}(0, 0, 0))
    c_hat = fill(SVector{3, eltype(α)}(0, 0, 0))
    n_hat = fill(SVector{3, eltype(α)}(0, 0, 0))
    Γs = zeros(3)
    λ = ones(length(panels_indicies))
    return SectionProperties(α, cl, cd, Γs, λ, panels_indicies, gammas, area, force, c_hat, n_hat, airfoil, contour)
end

"""
    grid_to_sections(grid, airfoils; 
                    ratios=zeros(2, size(grid, 2)-1, size(grid, 3)-1) .+ [0.5;0.75], 
                    contours=zeros(1,1))
"""
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
            gamma_start += size(sections[i][j].panels, 1)
        end
    end
end


"""
    nonlinear_analysis!(system)

Perform a nonlinear analysis on the system. This assumes that reference and freestream are already contained in the system

# Arguments
- `system::System`: The system to perform the analysis on
- `max_iter`: The maximum number of iterations to perform, default 1
- `tol`: The tolerance for convergence (based on absolute difference in Γ), default 1E-6
- `damping`: The damping factor for the iteration default 0.01
- `kwargs...`: Additional keyword arguments for steady analysis

"""
function nonlinear_analysis!(system; max_iter=1, tol=1E-6, damping=0.01, kwargs...)
    ref = system.reference[]
    fs = system.freestream[]
    return nonlinear_analysis!(system, ref, fs; max_iter=max_iter, tol=tol, damping=damping, kwargs...)
end

"""
    nonlinear_analysis!(system, ref, fs)

Perform a nonlinear analysis on the system.

# Arguments
- `system::System`: The system to perform the analysis on
- `max_iter`: The maximum number of iterations to perform, default 1
- `tol`: The tolerance for convergence (based on absolute difference in Γ), default 1E-6
- `damping`: The damping factor for the iteration default 0.01
- `kwargs...`: Additional keyword arguments for steady analysis
"""

function nonlinear_analysis!(system, ref, fs; max_iter=1, tol=1E-6, damping=0.01, kwargs...)
    if !system.near_field_analysis[]
        steady_analysis!(system, ref, fs; derivatives=false, near_field_analysis=true, kwargs...)
    end

    r, _ = lifting_line_geometry(system.grids)
    if max_iter < 1
        error("max_iter must be greater than 0")
    end

    for i in eachindex(system.sections)
        surface_properties = system.properties[i]
        for j in eachindex(system.sections[i])
            section = system.sections[i][j]
            Γavg = 0.0
            for k in eachindex(section.panels)
                Γavg += surface_properties[section.panels[k]].gamma
            end
            Γavg /= length(section.panels)
            section.Γs[1] = Γavg
            for k in eachindex(section.panels)
                section.λ[k] = surface_properties[section.panels[k]].gamma / Γavg
            end
        end
    end

    T = eltype(system.Γ)
    vel = zeros(T,3)
    vx = zeros(T, maximum(size.(system.surfaces,1)))
    vy = similar(vx)
    vz = similar(vx)
    x_c = similar(vx)

    for _ in 1:max_iter
        if _nonlinear_analysis!(system, r, damping, tol, vel, vx, vy, vz, x_c) # Calculate Γs, α, cl, and cd for each section
            break
        end
    end
    update_section_forces!(system, vel, vx, vy, vz, x_c)
    system.near_field_analysis[] = true
    system.derivatives[] = false
    return system
end

function _nonlinear_analysis!(system, r, damping, tol, vel, vx, vy, vz, x_c)
    converged = true
    vel .= 0.0
    vx .= 0.0
    vy .= 0.0
    vz .= 0.0
    x_c .= 0.0
    for i in eachindex(system.surfaces)
        sections = system.sections[i]
        !isassigned(sections,1) && continue # Ensure sections are assigned, if not do not perform nonlinear analysis on this surface

        surface = system.surfaces[i]
        nc = size(surface, 1)
        properties = system.properties[i]
        vx_view = view(vx,1:nc)
        vy_view = view(vy,1:nc)
        vz_view = view(vz,1:nc)
        x_c_view = view(x_c,1:nc)

        for j in eachindex(sections)
            total_chord = 0.0
            section = sections[j]
            section.c_hat[1] = surface[section.panels[end]].rbc - surface[section.panels[1]].rtc
            section.c_hat[1] /= norm(section.c_hat[1])
            section.n_hat[1] = cross(section.c_hat[1], surface[section.panels[end]].rbr - surface[section.panels[1]].rtl)
            section.n_hat[1] /= norm(section.n_hat[1])

            for k in eachindex(section.panels)
                total_chord += surface[section.panels[k]].chord
                vx_view[k], vy_view[k], vz_view[k] = properties[section.panels[k]].velocity_from_streamwise
            end

            for k in eachindex(section.panels)
                x_c_view[k] = norm(surface[section.panels[1]].rtc - surface[section.panels[k]].rcp) / total_chord
            end

            if length(section.panels) > 1
                vel[1] = FLOWMath.linear(x_c_view, vx_view, 0.5)
                vel[2] = FLOWMath.linear(x_c_view, vy_view, 0.5)
                vel[3] = FLOWMath.linear(x_c_view, vz_view, 0.5)
            else
                vel[1] = vx[1]
                vel[2] = vy[1]
                vel[3] = vz[1]
            end

            section.α[1] = atan(dot(vel, section.n_hat[1]), dot(vel, section.c_hat[1])) * (-1)^system.invert_normals[i]
            section.cl[1] = section.airfoil.clspline(section.α[1])
            section.cd[1] = section.airfoil.cdspline(section.α[1])

            r1 = r[i][1,j+1] - r[i][1,j]
            r2 = r[i][2,j+1] - r[i][2,j]
            r3 = r[i][3,j+1] - r[i][3,j]
            r_vec = SVector{3, eltype(r1)}(r1, r2, r3)

            # cross product
            cr = SVector{3, eltype(vel)}(
                vel[2]*r_vec[3] - vel[3]*r_vec[2],
                vel[3]*r_vec[1] - vel[1]*r_vec[3],
                vel[1]*r_vec[2] - vel[2]*r_vec[1]
            )
            section.Γs[2] = ((0.5*section.cl[1]*section.area * (dot(vel, section.n_hat[1])^2 + dot(vel, section.c_hat[1])^2)) / (sqrt(dot(cr, section.n_hat[1])^2+dot(cr, section.c_hat[1])^2))) * (-1)^system.invert_normals[i]
            section.Γs[3] = section.Γs[1] + damping*(section.Γs[2] - section.Γs[1])

            if isnan(section.Γs[2]) || isnan(section.Γs[3])
                error("NaN detected in _nonlinear_analysis! at the Γs calculation")
            end

            for k in eachindex(section.gammas)
                system.Γ[section.gammas[k]] = section.Γs[3] * section.λ[k]
            end

            if converged && (abs((section.Γs[3] - section.Γs[1])) > tol)
                converged = false
            end

            section.Γs[1] = section.Γs[3]
        end
    end
    call_near_field_forces!(system)
    return converged
end

function update_section_forces!(system, vel, vx, vy, vz, x_c)
    for i in eachindex(system.surfaces)
        surface = system.surfaces[i]
        nc = size(surface, 1)
        properties = system.properties[i]
        sections = system.sections[i]
        vx_view = view(vx,1:nc)
        vy_view = view(vy,1:nc)
        vz_view = view(vz,1:nc)
        x_c_view = view(x_c,1:nc)
        for j in eachindex(sections)
            total_chord = 0.0
            section = sections[j]

            for k in eachindex(section.panels)
                total_chord += surface[section.panels[k]].chord
                vx_view[k], vy_view[k], vz_view[k] = properties[section.panels[k]].velocity * system.reference[1].V
            end

            for k in eachindex(section.panels)
                x_c_view[k] = norm(surface[section.panels[1]].rtc - surface[section.panels[k]].rcp) / total_chord
            end

            if length(section.panels) > 1
                vel[1] = FLOWMath.linear(x_c_view, vx_view, 0.5)
                vel[2] = FLOWMath.linear(x_c_view, vy_view, 0.5)
                vel[3] = FLOWMath.linear(x_c_view, vz_view, 0.5)
            else
                vel[1] = vx[1]
                vel[2] = vy[1]
                vel[3] = vz[1]
            end
            
            v_mag = norm(vel)
            lift = section.cl[1] * v_mag^2 * RHO / 2 * total_chord
            drag = section.cd[1] * v_mag^2 * RHO / 2 * total_chord
            R = transpose([section.n_hat[1] cross(section.n_hat[1], section.c_hat[1]) section.c_hat[1]])
            section.force[1] = R * SVector{3,eltype(vel)}(lift, 0.0, drag) ./ RHO ./ system.reference[1].V
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

function call_near_field_forces!(system)
    properties = system.properties
    surfaces = system.surfaces
    wakes = system.wakes
    ref = system.reference[]
    fs = system.freestream[]
    Γ = system.Γ
    symmetric = system.symmetric
    nwake = system.nwake
    surface_id = system.surface_id
    wake_finite_core = system.wake_finite_core
    trailing_vortices = system.trailing_vortices
    xhat = system.xhat[]
    
    near_field_forces!(properties, surfaces, wakes, ref, fs, Γ;
                dΓdt = system.dΓdt,
                additional_velocity = nothing,
                Vh = system.Vh,
                Vv = system.Vv,
                symmetric = symmetric,
                nwake = nwake,
                surface_id = surface_id,
                wake_finite_core = wake_finite_core,
                wake_shedding_locations = nothing,#system.wake_shedding_location,
                trailing_vortices = trailing_vortices,
                xhat = xhat)
end