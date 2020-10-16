"""
    project_panels(panels, freestream)

Rotates `panels` from the body frame into a wind-aligned coordinate system
"""
project_panels

@inline function project_panels(panels, fs)

    R = body_to_wind(fs)

    panels = rotate.(panels, Ref(R))

    return panels
end

"""
    far_field_drag(panels, ref, fs, symmetric, Γ)

Computes induced drag using the Trefftz plane (far field method).
"""
far_field_drag

@inline function far_field_drag(panels::AbstractVector{<:Horseshoe}, ref, fs, symmetric, Γ)

    # rotate panels into wind coordinate system
    panels = project_panels(panels, fs)

    N = length(panels)
    TF = promote_type(eltype(eltype(panels)), eltype(Γ))

    # add up drag
    Di = zero(TF)
    for j = 1:N
        for i = 1:N
            Di += panel_induced_drag(panels[j], Γ[j], panels[i], Γ[i], symmetric)
        end
    end

    # apply symmetry
    if symmetric
        Di *= 2
    end

    # normalize
    CDi = Di / (QINF*ref.S)

    return CDi
end


@inline function far_field_drag(panels::AbstractVector{<:Ring}, ref, fs, symmetric, Γ)

    # get panels that shed trailing vortices
    ip = findall((p)->p.trailing, panels)

    # rotate panels into wind coordinate system
    wake_panels = project_panels(panels[ip], fs)

    N = length(wake_panels)
    TF = promote_type(eltype(eltype(wake_panels)), eltype(Γ))

    Di = zero(TF)
    for j = 1:N
        for i = 1:N
            Di += panel_induced_drag(wake_panels[j], Γ[ip[j]], wake_panels[i], Γ[ip[i]], symmetric)
        end
    end

    if symmetric
        Di *= 2
    end

    # normalize
    CDi = Di / (QINF*ref.S)

    return CDi
end
