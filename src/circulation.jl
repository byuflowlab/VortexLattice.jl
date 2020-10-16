# --- left hand side - AIC matrix --- #

"""
    influence_coefficients(panels; symmetric=false, xhat=[1,0,0])

Construct the aerodynamic influence coefficient matrix.
"""
function influence_coefficients(panels; symmetric=false, xhat=SVector(1, 0, 0))

    N = length(panels)
    TF = eltype(eltype(panels))
    AIC = Matrix{TF}(undef, N, N)

    return influence_coefficients!(AIC, panels; symmetric=symmetric, xhat=xhat)
end

"""
    influence_coefficients!(AIC, panels, symmetric)

Pre-allocated version of `influence_coefficients`
"""
function influence_coefficients!(AIC, panels; symmetric=false, xhat=SVector(1, 0, 0))

    N = length(panels)

    # loop over control points
    for i = 1:N

        # control point location
        rcp = controlpoint(panels[i])

        # normal vector body axis
        nhat = normal(panels[i])

        # loop over bound vortices
        for j = 1:N
            Vhat = induced_velocity(rcp, panels[j], symmetric, xhat)
            AIC[i, j] = dot(Vhat, nhat)
        end
    end

    return AIC
end

# --- right hand side - normal velocities at control points --- #

"""
    normal_velocity(panels, ref, fs)

Compute the normal component of the external velocity along the geometry.
This forms the right hand side of the circulation linear system solve.
"""
function normal_velocity(panels, ref, fs)

    N = length(panels)
    TF = promote_type(eltype(eltype(panels)), eltype(ref), eltype(fs))
    b = Vector{TF}(undef, N)

    return normal_velocity!(b, panels, ref, fs)
end

"""
    normal_velocity!(b, panels, ref, fs)

Non-allocating version of `normal_velocity`
"""
function normal_velocity!(b, panels, ref, fs)

    N = length(panels)

    # iterate through panels
    for i = 1:N

        # control point
        rcp = controlpoint(panels[i])

        # normal vector
        nhat = normal(panels[i])

        # external velocity
        Vext = external_velocity(fs, rcp, ref.r)

        # right hand side vector
        b[i] = -dot(Vext, nhat)

    end

    return b
end

"""
    normal_velocity_derivatives(panels, ref, fs)

Compute the normal component of the external velocity along the geometry and its
derivatives with respect to (alpha, beta, p, q, r). This forms the right hand
side of the circulation linear system solve (and its derivatives).
"""
function normal_velocity_derivatives(panels, ref, fs)

    N = length(panels)
    TF = promote_type(eltype(eltype(panels)), eltype(ref), eltype(fs))

    # RHS vector
    b = Vector{TF}(undef, N)

    # derivatives of RHS wrt freestream variables
    b_a = Vector{TF}(undef, N)
    b_b = Vector{TF}(undef, N)
    b_pb = Vector{TF}(undef, N)
    b_qb = Vector{TF}(undef, N)
    b_rb = Vector{TF}(undef, N)

    # pack up derivatives
    db = (b_a, b_b, b_pb, b_qb, b_rb)

    return normal_velocity_derivatives!(b, db, panels, ref, fs)
end

"""
    normal_velocity_derivatives!(b, db, panels, ref, fs)

Non-allocating version of `normal_velocity_derivatives`
"""
function normal_velocity_derivatives!(b, db, panels, ref, fs)

    N = length(panels)

    # unpack derivatives
    (b_a, b_b, b_pb, b_qb, b_rb) = db

    # iterate through panels
    for i = 1:N

        # control point
        rcp = controlpoint(panels[i])

        # normal vector
        nhat = normal(panels[i])

        # external velocity and its derivatives
        Vext, dVext = external_velocity_derivatives(fs, rcp, ref.r)

        # unpack derivatives
        Vext_a, Vext_b, Vext_pb, Vext_qb, Vext_rb = dVext

        # right hand side vector
        b[i] = -dot(Vext, nhat)

        # associated derivatives
        b_a[i] = -dot(Vext_a, nhat)
        b_b[i] = -dot(Vext_b, nhat)
        b_pb[i] = -dot(Vext_pb, nhat)
        b_qb[i] = -dot(Vext_qb, nhat)
        b_rb[i] = -dot(Vext_rb, nhat)

    end

    db = (b_a, b_b, b_pb, b_qb, b_rb)

    return b, db
end

# --- circulation solve --- #

"""
    circulation(AIC, b)

Solve for the circulation distribution.
"""
circulation(AIC, b) = AIC\b

"""
    circulation!(AIC, b)

Pre-allocated version of `circulation`
"""
circulation!(Γ, AIC, b) = ldiv!(Γ, factorize(AIC), b)

"""
    circulation_derivatives(AIC, b, db)

Solve for the circulation distribution and its derivatives with respect to
(alpha, beta, p, q, r)    .
"""
function circulation_derivatives(AIC, b, db)

    # unpack derivatives
    (b_a, b_b, b_pb, b_qb, b_rb) = db

    # factorize AIC matrix (since we'll be reusing it)
    fAIC = factorize(AIC)

    # solve for circulation and its derivatives
    Γ = fAIC\b

    Γ_a = fAIC\b_a
    Γ_b = fAIC\b_b
    Γ_pb = fAIC\b_pb
    Γ_qb = fAIC\b_qb
    Γ_rb = fAIC\b_rb

    # pack up derivatives
    dΓ = (Γ_a, Γ_b, Γ_pb, Γ_qb, Γ_rb)

    return Γ, dΓ
end
