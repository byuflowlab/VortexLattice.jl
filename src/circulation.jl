# --- right hand side - normal velocities at control points --- #

"""
    normal_velocity(surface[s], reference, freestream, additional_velocity)

Compute the downwash due to the freestream velocity at each control point for a
single surface or vector of surfaces.

This forms the right hand side of the circulation linear system solve.
"""
normal_velocity

# one surface
function normal_velocity(surface::AbstractMatrix, ref, fs, additional_velocity)

    N = length(surface)
    TF = promote_type(eltype(eltype(surface)), eltype(ref), eltype(fs))
    w = Vector{TF}(undef, N)

    return normal_velocity!(w, surface, ref, fs, additional_velocity)
end

# multiple surfaces
function normal_velocity(surfaces::AbstractVector{<:AbstractMatrix}, ref, fs, additional_velocity)

    N = sum(length.(surfaces))
    TF = promote_type(eltype(eltype(eltype(surfaces))), eltype(ref), eltype(fs))
    w = Vector{TF}(undef, N)

    return normal_velocity!(w, surfaces, ref, fs, additional_velocity)
end

"""
    normal_velocity!(w, surface[s], reference, freestream)

Non-allocating version of `normal_velocity`
"""
normal_velocity!

# one surface
function normal_velocity!(w, surface::AbstractMatrix, ref, fs, additional_velocity)

    N = length(surface)
    c = CartesianIndices(surface)

    # iterate through panels
    for i = 1:N

        I = c[i]

        # control point
        rcp = controlpoint(surface[I])

        # normal vector
        nhat = normal(surface[I])

        # external velocity
        Vext = external_velocity(rcp, fs, ref.r, additional_velocity)

        # right hand side vector
        w[i] = -dot(Vext, nhat)

    end

    return w
end

# multiple surfaces
function normal_velocity!(w, surfaces::AbstractVector{<:AbstractMatrix}, ref, fs, additional_velocity)

    nsurf = length(surfaces)

    # index for keeping track of where we are in the w vector
    iw = 0

    # loop through receiving surfaces
    for i = 1:nsurf

        # number of panels on this surface
        n = length(surfaces[i])

        # create view into RHS vector
        vw = view(w, iw+1:iw+n)

        # fill in RHS vector
        normal_velocity!(vw, surfaces[i], ref, fs, additional_velocity)

        # increment position in AIC matrix
        iw += n
    end

    return w
end

"""
    normal_velocity_derivatives(surface[s], ref, fs)

Compute the downwash due to the freestream velocity at each control point for a
single surface or vector of surfaces and its derivatives with respect to the
freestream parameters.

This forms the right hand side of the circulation linear system solve (and its derivatives).
"""
normal_velocity_derivatives

# one surface
function normal_velocity_derivatives(surface::AbstractMatrix, ref, fs, additional_velocity)

    N = length(surface)
    TF = promote_type(eltype(eltype(surface)), eltype(ref), eltype(fs))

    # RHS vector
    w = Vector{TF}(undef, N)

    # derivatives of RHS wrt freestream variables
    w_a = Vector{TF}(undef, N)
    w_b = Vector{TF}(undef, N)
    w_p = Vector{TF}(undef, N)
    w_q = Vector{TF}(undef, N)
    w_r = Vector{TF}(undef, N)

    # pack up derivatives
    dw = (w_a, w_b, w_p, w_q, w_r)

    return normal_velocity_derivatives!(w, dw, surface, ref, fs, additional_velocity)
end

# multiple surfaces
function normal_velocity_derivatives(surfaces::AbstractVector{<:AbstractMatrix}, ref, fs, additional_velocity)

    N = sum(length.(surfaces))
    TF = promote_type(eltype(eltype(eltype(surfaces))), eltype(ref), eltype(fs))

    # RHS vector
    w = Vector{TF}(undef, N)

    # derivatives of RHS wrt freestream variables
    w_a = Vector{TF}(undef, N)
    w_b = Vector{TF}(undef, N)
    w_p = Vector{TF}(undef, N)
    w_q = Vector{TF}(undef, N)
    w_r = Vector{TF}(undef, N)

    # pack up derivatives
    dw = (w_a, w_b, w_p, w_q, w_r)

    return normal_velocity_derivatives!(w, dw, surfaces, ref, fs, additional_velocity)
end

"""
    normal_velocity_derivatives!(b, dw, surface[s], ref, fs, additional_velocity)

Non-allocating version of `normal_velocity_derivatives`
"""
normal_velocity_derivatives!

# single surface
function normal_velocity_derivatives!(w, dw, surface, ref, fs, additional_velocity)

    N = length(surface)

    c = CartesianIndices(surface)

    # unpack derivatives
    (w_a, w_b, w_p, w_q, w_r) = dw

    # iterate through surface
    for i = 1:N

        I = c[i]

        # control point
        rcp = controlpoint(surface[I])

        # normal vector
        nhat = normal(surface[I])

        # external velocity and its derivatives
        Vext, dVext = external_velocity_derivatives(rcp, fs, ref.r, additional_velocity)

        # unpack derivatives
        Vext_a, Vext_b, Vext_pb, Vext_qb, Vext_rb = dVext

        # right hand side vector
        w[i] = -dot(Vext, nhat)

        # associated derivatives
        w_a[i] = -dot(Vext_a, nhat)
        w_b[i] = -dot(Vext_b, nhat)
        w_p[i] = -dot(Vext_pb, nhat)
        w_q[i] = -dot(Vext_qb, nhat)
        w_r[i] = -dot(Vext_rb, nhat)

    end

    dw = (w_a, w_b, w_p, w_q, w_r)

    return w, dw
end

# multiple surfaces
function normal_velocity_derivatives!(w, dw, surfaces::AbstractVector{<:AbstractMatrix}, ref, fs, additional_velocity)

    nsurf = length(surfaces)

    # unpack derivatives
    (w_a, w_b, w_p, w_q, w_r) = dw

    # index for keeping track of where we are in the w vector
    iw = 0

    # loop through receving panels
    for i = 1:nsurf

        # number of panels on this surface
        n = length(surfaces[i])

        # create view into RHS vector and its derivatives
        vb = view(w, iw+1:iw+n)
        vdb = view.(dw, Ref(iw+1:iw+n))

        # fill in RHS vector and its derivatives
        normal_velocity_derivatives!(vb, vdb, surfaces[i], ref, fs, additional_velocity)

        # increment position in AIC matrix
        iw += n
    end

    return w, dw
end

# --- circulation solve --- #

"""
    circulation(AIC, w)

Solve for the circulation distribution.
"""
circulation(AIC, w) = AIC\w

"""
    circulation!(Γ, AIC, w)

Pre-allocated version of `circulation`
"""
circulation!(Γ, AIC, w) = ldiv!(Γ, lu(AIC), w)

"""
    circulation_derivatives(AIC, w, dw)

Solve for the circulation distribution and its derivatives with respect to
the freestream parameters.
"""
function circulation_derivatives(AIC, w, dw)

    # unpack derivatives
    (w_a, w_b, w_p, w_q, w_r) = dw

    # factorize AIC matrix (since we'll be reusing it)
    fAIC = factorize(AIC)

    # solve for circulation and its derivatives
    Γ = fAIC\w

    Γ_a = fAIC\w_a
    Γ_b = fAIC\w_b
    Γ_pb = fAIC\w_p
    Γ_qb = fAIC\w_q
    Γ_rb = fAIC\w_r

    # pack up derivatives
    dΓ = (Γ_a, Γ_b, Γ_pb, Γ_qb, Γ_rb)

    return Γ, dΓ
end

"""
    circulation_derivatives!(Γ, dΓ, AIC, w, dw)

Pre-allocated version of `circulation_derivatives`
"""
function circulation_derivatives!(Γ, dΓ, AIC, w, dw)

    # unpack derivatives
    (Γ_a, Γ_b, Γ_p, Γ_q, Γ_r) = dΓ
    (w_a, w_b, w_p, w_q, w_r) = dw

    # factorize AIC matrix (since we'll be reusing it)
    fAIC = lu(AIC)

    # solve for circulation and its derivatives
    ldiv!(Γ, fAIC, w)

    ldiv!(Γ_a, fAIC, w_a)
    ldiv!(Γ_b, fAIC, w_b)
    ldiv!(Γ_p, fAIC, w_p)
    ldiv!(Γ_q, fAIC, w_q)
    ldiv!(Γ_r, fAIC, w_r)

    return Γ, dΓ
end
