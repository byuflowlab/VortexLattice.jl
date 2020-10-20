# --- left hand side - AIC matrix --- #

"""
    influence_coefficients(surface[s], symmetric; kwargs...)

Construct the aerodynamic influence coefficient matrix for a single surface or
a vector of surfaces.

# Keyword Arguments
 - `xhat`: direction in which trailing vortices are shed, defaults to [1,0,0]

# Additional Keyword Arguments for Multiple Surfaces
 - `surface_id`: ID for each surface, defaults to `1:length(surfaces)`
"""
influence_coefficients

# one surface
function influence_coefficients(surface::AbstractMatrix, symmetric; xhat=SVector(1, 0, 0))

    N = length(surface)
    TF = eltype(eltype(surface))
    AIC = Matrix{TF}(undef, N, N)

    return influence_coefficients!(AIC, surface, symmetric; xhat=xhat)
end

# multiple surfaces
function influence_coefficients(surfaces::AbstractVector{<:AbstractMatrix}, symmetric;
    xhat=SVector(1, 0, 0), surface_id=1:length(surfaces))

    N = sum(length.(surfaces))
    TF = eltype(eltype(eltype(surfaces)))
    AIC = Matrix{TF}(undef, N, N)

    return influence_coefficients!(AIC, surfaces, symmetric; xhat=xhat, surface_id=surface_id)
end

"""
    influence_coefficients!(AIC, surface[s], symmetric; kwargs...)

Pre-allocated version of `influence_coefficients`
"""
influence_coefficients!

# one surface
function influence_coefficients!(AIC, surface::AbstractMatrix, symmetric; xhat=SVector(1, 0, 0))

    receiving = surface
    sending = surface
    same_surface = true

    influence_coefficients!(AIC, surface, surface, same_surface, symmetric, xhat)

    return AIC
end

# multiple surfaces
function influence_coefficients!(AIC, surfaces::AbstractVector{<:AbstractMatrix},
    symmetric; xhat=SVector(1, 0, 0), surface_id=1:length(surfaces))

    nsurf = length(surfaces)

    # indices for keeping track of where we are in the AIC matrix
    iAIC = 0
    jAIC = 0

    # loop through receving surfaces
    for i = 1:nsurf
        receiving = surfaces[i]

        # extract number of panels on this receiving surface
        nr = length(receiving) # number of panels on this receiving surface

        # loop through sending surfaces
        jAIC = 0
        for j = 1:nsurf
            sending = surfaces[j]

            # extract number of panels on this sending surface
            ns = length(sending)

            # extract portion of AIC matrix for the two surfaces
            vAIC = view(AIC, iAIC+1:iAIC+nr, jAIC+1:jAIC+ns)

            # check if it's the same surface
            same_surface = surface_id[i] == surface_id[j]

            # populate entries in the AIC matrix
            influence_coefficients!(vAIC, receiving, sending, same_surface, symmetric, xhat)

            # increment position in AIC matrix
            jAIC += ns
        end
        # increment position in AIC matrix
        iAIC += nr
    end

    return AIC
end

# one surface on another surface
function influence_coefficients!(AIC, receiving, sending, same_surface, symmetric, xhat)

    Nr = length(receiving)
    Ns = length(sending)

    cr = CartesianIndices(receiving)
    cs = CartesianIndices(sending)
    nchordwise = size(sending,1)

    # loop over receiving panels
    for i = 1:Nr
        I = cr[i]

        # control point location
        rcp = controlpoint(receiving[I])

        # normal vector body axis
        nhat = normal(receiving[I])

        # loop over sending panels
        for j = 1:Ns
            J = cs[j]

            if typeof(sending[J]) <: Horseshoe
                include_bound = true
                Vhat = induced_velocity(rcp, sending[J], xhat, symmetric, same_surface,
                    include_bound)
            elseif typeof(sending[j]) <: Ring
                trailing = J[1] == nchordwise
                include_top = true
                include_bottom = true
                Vhat = induced_velocity(rcp, sending[J], xhat, symmetric, same_surface,
                    trailing, include_top, include_bottom)
            end

            AIC[i, j] = dot(Vhat, nhat)
        end
    end

    return AIC
end

# --- right hand side - normal velocities at control points --- #

"""
    normal_velocity(surface[s], ref, fs)

Compute the normal component of the external velocity for a single surface or
for a vector of surfaces.

This forms the right hand side of the circulation linear system solve.
"""
normal_velocity

# one surface
function normal_velocity(surface::AbstractMatrix, ref, fs)

    N = length(surface)
    TF = promote_type(eltype(eltype(surface)), eltype(ref), eltype(fs))
    b = Vector{TF}(undef, N)

    return normal_velocity!(b, surface, ref, fs)
end

# multiple surfaces
function normal_velocity(surfaces::AbstractVector{<:AbstractMatrix}, ref, fs)

    N = sum(length.(surfaces))
    TF = promote_type(eltype(eltype(eltype(surfaces))), eltype(ref), eltype(fs))
    b = Vector{TF}(undef, N)

    return normal_velocity!(b, surfaces, ref, fs)
end

"""
    normal_velocity!(b, surface[s], ref, fs)

Non-allocating version of `normal_velocity`
"""
normal_velocity!

# one surface
function normal_velocity!(b, surface::AbstractMatrix, ref, fs)

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
        Vext = external_velocity(fs, rcp, ref.r)

        # right hand side vector
        b[i] = -dot(Vext, nhat)

    end

    return b
end

# multiple surfaces
function normal_velocity!(b, surfaces::AbstractVector{<:AbstractMatrix}, ref, fs)

    nsurf = length(surfaces)

    # index for keeping track of where we are in the b vector
    ib = 0

    # loop through receiving surfaces
    for i = 1:nsurf

        # number of panels on this surface
        n = length(surfaces[i])

        # create view into RHS vector
        vb = view(b, ib+1:ib+n)

        # fill in RHS vector
        normal_velocity!(vb, surfaces[i], ref, fs)

        # increment position in AIC matrix
        ib += n
    end

    return b
end

"""
    normal_velocity_derivatives(surface[s], ref, fs)

Compute the normal component of the external velocity for a single surface or
for a vector of surfaces and its derivatives with respect to (alpha, beta, p, q, r).

This forms the right hand side of the circulation linear system solve (and its derivatives).
"""
normal_velocity_derivatives

# one surface
function normal_velocity_derivatives(surface::AbstractMatrix, ref, fs)

    N = length(surface)
    TF = promote_type(eltype(eltype(surface)), eltype(ref), eltype(fs))

    # RHS vector
    b = Vector{TF}(undef, N)

    # derivatives of RHS wrt freestream variables
    b_a = Vector{TF}(undef, N)
    b_b = Vector{TF}(undef, N)
    b_p = Vector{TF}(undef, N)
    b_q = Vector{TF}(undef, N)
    b_r = Vector{TF}(undef, N)

    # pack up derivatives
    db = (b_a, b_b, b_p, b_q, b_r)

    return normal_velocity_derivatives!(b, db, surface, ref, fs)
end

# multiple surfaces
function normal_velocity_derivatives(surfaces::AbstractVector{<:AbstractMatrix}, ref, fs)

    N = sum(length.(surfaces))
    TF = promote_type(eltype(eltype(eltype(surfaces))), eltype(ref), eltype(fs))

    # RHS vector
    b = Vector{TF}(undef, N)

    # derivatives of RHS wrt freestream variables
    b_a = Vector{TF}(undef, N)
    b_b = Vector{TF}(undef, N)
    b_p = Vector{TF}(undef, N)
    b_q = Vector{TF}(undef, N)
    b_r = Vector{TF}(undef, N)

    # pack up derivatives
    db = (b_a, b_b, b_p, b_q, b_r)

    return normal_velocity_derivatives!(b, db, surfaces, ref, fs)
end

"""
    normal_velocity_derivatives!(b, db, surface[s], ref, fs)

Non-allocating version of `normal_velocity_derivatives`
"""
normal_velocity_derivatives

# single surface
function normal_velocity_derivatives!(b, db, surface, ref, fs)

    N = length(surface)

    c = CartesianIndices(surface)

    # unpack derivatives
    (b_a, b_b, b_p, b_q, b_r) = db

    # iterate through surface
    for i = 1:N

        I = c[i]

        # control point
        rcp = controlpoint(surface[I])

        # normal vector
        nhat = normal(surface[I])

        # external velocity and its derivatives
        Vext, dVext = external_velocity_derivatives(fs, rcp, ref.r)

        # unpack derivatives
        Vext_a, Vext_b, Vext_pb, Vext_qb, Vext_rb = dVext

        # right hand side vector
        b[i] = -dot(Vext, nhat)

        # associated derivatives
        b_a[i] = -dot(Vext_a, nhat)
        b_b[i] = -dot(Vext_b, nhat)
        b_p[i] = -dot(Vext_pb, nhat)
        b_q[i] = -dot(Vext_qb, nhat)
        b_r[i] = -dot(Vext_rb, nhat)

    end

    db = (b_a, b_b, b_p, b_q, b_r)

    return b, db
end

# multiple surfaces
function normal_velocity_derivatives!(b, db, surfaces::AbstractVector{<:AbstractMatrix}, ref, fs)

    nsurf = length(surfaces)

    # unpack derivatives
    (b_a, b_b, b_p, b_q, b_r) = db

    # index for keeping track of where we are in the b vector
    ib = 0

    # loop through receving panels
    for i = 1:nsurf

        # number of panels on this surface
        n = length(surfaces[i])

        # create view into RHS vector and its derivatives
        vb = view(b, ib+1:ib+n)
        vdb = view.(b, Ref(ib+1:ib+n))

        # fill in RHS vector and its derivatives
        normal_velocity_derivatives!(vb, vdb, surfaces[i], ref, fs)

        # increment position in AIC matrix
        ib += n
    end

    return b
end

# --- circulation solve --- #

"""
    circulation(AIC, b)

Solve for the circulation distribution.
"""
circulation(AIC, b) = AIC\b

"""
    circulation_derivatives(AIC, b, db)

Solve for the circulation distribution and its derivatives with respect to
(alpha, beta, p, q, r)    .
"""
function circulation_derivatives(AIC, b, db)

    # unpack derivatives
    (b_a, b_b, b_p, b_q, b_r) = db

    # factorize AIC matrix (since we'll be reusing it)
    fAIC = factorize(AIC)

    # solve for circulation and its derivatives
    Γ = fAIC\b

    Γ_a = fAIC\b_a
    Γ_b = fAIC\b_b
    Γ_pb = fAIC\b_p
    Γ_qb = fAIC\b_q
    Γ_rb = fAIC\b_r

    # pack up derivatives
    dΓ = (Γ_a, Γ_b, Γ_pb, Γ_qb, Γ_rb)

    return Γ, dΓ
end
