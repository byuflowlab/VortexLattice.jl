# --- left hand side - AIC matrix --- #

"""
    influence_coefficients(panels, symmetric)

Construct the aerodynamic influence coefficient matrix.
"""
@inline function influence_coefficients(panels, symmetric, args...)

    N = length(panels)
    TF = eltype(eltype(panels))
    AIC = Matrix{TF}(undef, N, N)

    return influence_coefficients!(AIC, panels, symmetric, args...)
end

"""
    influence_coefficients!(AIC, panels::AbstractVector{<:Ring}, symmetric, args...)

Pre-allocated version of `influence_coefficients`
"""
@inline function influence_coefficients!(AIC, panels, symmetric, args...)

    N = length(panels)

    # loop over control points
    for i = 1:N

        # control point location
        rcp = controlpoint(panels[i])

        # normal vector body axis
        nhat = normal(panels[i])

        # loop over bound vortices
        for j = 1:N
            Vhat = induced_velocity(rcp, panels[j], symmetric, args...)
            AIC[i, j] = dot(Vhat, nhat)
        end
    end

    return AIC
end

# --- right hand side - normal velocities at control points --- #

"""
    normal_velocity(panels, reference, freestream)

Compute the normal component of the external velocity along the geometry.
This forms the right hand side of the circulation linear system solve.
"""
@inline function normal_velocity(panels, reference, freestream)

    N = length(panels)
    TF = promote_type(eltype(eltype(panels)), eltype(reference), eltype(freestream))
    b = Vector{TF}(undef, N)

    return normal_velocity!(b, panels, reference, freestream)
end

"""
    normal_velocity!(b, panels, reference, freestream)

Non-allocating version of `normal_velocity`
"""
@inline function normal_velocity!(b, panels, reference, freestream)

    N = length(panels)

    # iterate through panels
    for i = 1:N

        # control point
        rcp = controlpoint(panels[i])

        # normal vector
        nhat = normal(panels[i])

        # external velocity
        Vext = external_velocity(freestream, rcp, reference.rcg)

        # right hand side vector
        b[i] = -dot(Vext, nhat)

    end

    return b
end

# --- circulation solve --- #

"""
    circulation(panels, reference, freestream, symmetric)

Solve for the circulation distribution.
"""
@inline function circulation(panels, reference, freestream, symmetric, args...)

    AIC = influence_coefficients(panels, symmetric, args...)
    b = normal_velocity(panels, reference, freestream)

    return AIC\b
end

"""
    circulation!(Γ, AIC, b, panels, reference, freestream, symmetric)

Pre-allocated version of `circulation`
"""
@inline function circulation!(Γ, AIC, b, panels, reference, freestream, symmetric, args...)

    AIC = influence_coefficients!(AIC, panels, symmetric, args...)
    b = normal_velocity!(b, panels, reference, freestream)

    return ldiv!(Γ, factorize(AIC), b)
end

# ------------ run method --------------------

"""
    vlm(panels, reference, freestream, symmetric)

Runs the vortex lattice method.
"""
function vlm(panels, reference, freestream, symmetric)

    Γ = circulation(panels, reference, freestream, symmetric)
    CF, CM, paneloutputs = forces_moments(panels, reference, freestream, Γ, symmetric)
    CDiff, trefftz_panels = trefftz_induced_drag(panels, reference, freestream, Γ, symmetric)

    return Outputs(CF, CM, CDiff, paneloutputs)
end
