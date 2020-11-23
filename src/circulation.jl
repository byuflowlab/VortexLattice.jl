"""
    influence_coefficients(surface, symmetric; kwargs...)

Construct the aerodynamic influence coefficient matrix for a single surface.

# Arguments:
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `symmetric`: Flag indicating whether a mirror image of the panels in `surface`
    should be used when calculating induced velocities.

# Keyword Arguments
 - `trailing_vortices`: Flag which may be used to enable/disable trailing vortices,
    trailing vortices are enabled by default
 - `xhat`: direction in which trailing vortices are shed if wake panels are not
    used to model the wake, defaults to [1, 0, 0]
"""
function influence_coefficients(surface::AbstractMatrix, symmetric; kwargs...)

    N = length(surface)
    TF = eltype(eltype(surface))
    AIC = Matrix{TF}(undef, N, N)

    return influence_coefficients!(AIC, surface, symmetric; kwargs...)
end

"""
    influence_coefficients(surfaces, surface_id, symmetric; kwargs...)

Construct the aerodynamic influence coefficient matrix for multiple surfaces.

# Arguments:
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `surface_id`: ID for each surface.  May be used to deactivate the finite core
    model by setting all surface ID's to the same value.
 - `symmetric`: Flags indicating whether a mirror image (across the y-axis) should
    be used when calculating induced velocities, defaults to `false` for each surface

# Keyword Arguments
 - `trailing_vortices`: Flags to enable/disable trailing vortices, defaults to `true`
    for each surface
 - `xhat`: direction in which trailing vortices are shed if wake panels are not
    used to model the wake, defaults to [1, 0, 0]
"""
function influence_coefficients(surfaces::AbstractVector{<:AbstractMatrix}, surface_id,
    symmetric; kwargs...)

    N = sum(length.(surfaces))
    TF = eltype(eltype(eltype(surfaces)))
    AIC = Matrix{TF}(undef, N, N)

    return influence_coefficients!(AIC, surfaces, surface_id, symmetric; kwargs...)
end

"""
    influence_coefficients!(AIC, surface, symmetric; kwargs...)

Pre-allocated version of `influence_coefficients`
"""
function influence_coefficients!(AIC, surface::AbstractMatrix, symmetric; kwargs...)

    receiving = surface
    sending = surface
    same_id = true

    influence_coefficients!(AIC, surface, surface, symmetric, same_id; kwargs...)

    return AIC
end

"""
    influence_coefficients!(AIC, surfaces, surface_id, symmetric; kwargs...)

Pre-allocated version of `influence_coefficients`
"""
function influence_coefficients!(AIC, surfaces::AbstractVector{<:AbstractMatrix},
    surface_id, symmetric; trailing_vortices=fill(true, length(surfaces)),
    xhat=SVector(1, 0, 0))

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
            same_id = surface_id[i] == surface_id[j]

            # populate entries in the AIC matrix
            influence_coefficients!(vAIC, receiving, sending, symmetric[j],
                same_id; trailing_vortices = trailing_vortices[j], xhat = xhat)

            # increment position in AIC matrix
            jAIC += ns
        end
        # increment position in AIC matrix
        iAIC += nr
    end

    return AIC
end

"""
    influence_coefficients!(AIC, receiving, sending, symmetric, same_id; kwargs...)

Compute the AIC coefficients corresponding to the influence of the panels in
`sending` on the panels in `receiving`.
"""
@inline function influence_coefficients!(AIC, receiving, sending, symmetric,
    same_id; trailing_vortices=true, xhat=SVector(1, 0, 0))

    Nr = length(receiving)
    Ns = length(sending)

    cr = CartesianIndices(receiving)
    cs = CartesianIndices(sending)
    nchordwise = size(sending, 1)

    finite_core = !same_id

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

            trailing_edge = J[1] == nchordwise

            trailing = trailing_vortices && has_trailing_vortices(sending[J], trailing_edge)

            # get normalized induced velocity for the panel
            Vhat = sum(panel_induced_velocity(rcp, sending[J], trailing;
                finite_core = finite_core, xhat = xhat))

            # also get normalized induced velocity for its mirror image
            if symmetric
                Vhat += sum(panel_induced_velocity(rcp, sending[J], trailing;
                    finite_core = finite_core, reflect = true, xhat = xhat))
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
normal_velocity_derivatives!

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
        vdb = view.(db, Ref(ib+1:ib+n))

        # fill in RHS vector and its derivatives
        normal_velocity_derivatives!(vb, vdb, surfaces[i], ref, fs)

        # increment position in AIC matrix
        ib += n
    end

    return b
end

# --- wake_normal_velocity --- #

"""
    wake_normal_velocity(surface, wake, trailing_vortices, symmetric; kwargs...)

Compute the normal component of the velocity induced on a surface by its own
wake panels

This forms part of the right hand side of the circulation linear system solve.

# Arguments
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `wake`: Matrix of wake panels of shape (nw, ns) where `nw` is the number of
    chordwise wake panels and `ns` is the number of spanwise panels
 - `trailing_vortices`: flags indicating whether trailing vortices should be shed from
    the last chordwise panels in `wake`.
 - `symmetric`: Flag indicating whether a mirror image of the panels in `wake`
    should be used when calculating induced velocities.

# Keyword Arguments
 - `nwake`: number of chordwise wake panels to use from `wake`, defaults to all
    wake panels
 - `xhat`: direction in which trailing vortices are shed, defaults to [1, 0, 0]
"""
function wake_normal_velocity(surface::AbstractMatrix, wake::AbstractMatrix,
    trailing_vortices, symmetric; nwake=size(wake,1), xhat=SVector(1,0,0))

    TF = promote_type(eltype(eltype(surface)), eltype(eltype(wake)))
    N = length(surface)
    b = zeros(TF, N)

    return add_wake_normal_velocity!(b, surface, wake, trailing_vortices, symmetric; nwake=nwake, xhat=xhat)
end

"""
    wake_normal_velocity(surfaces, wakes, surface_id, trailing_vortices, symmetric; kwargs...)

Compute the normal component of the velocity induced on multiple surfaces by
their wake panels

This forms part of the right hand side of the circulation linear system solve.

# Arguments
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `wakes`: Vector of wakes, represented by matrices of panels of shape (nw, ns)
    where `nw` is the number of chordwise wake panels and `ns` is the number of
    spanwise panels
 - `surface_id`: ID for each surface.  May be used to deactivate the finite core
    model by setting all surface ID's to the same value.
 - `trailing_vortices`: flags indicating whether trailing vortices should be shed from
    the last chordwise panel of each wake.
 - `symmetric`: Flag indicating whether a mirror image of the panels in `wake`
    should be used when calculating induced velocities.

# Keyword Arguments
 - `nwake`: number of chordwise wake panels to use from each `wake`, defaults to
    all wake panels
 - `xhat`: direction in which trailing vortices are shed, defaults to [1, 0, 0]
"""
function wake_normal_velocity(surfaces::AbstractVector{<:AbstractMatrix},
    wakes::AbstractVector{<:AbstractMatrix}, surface_id, trailing_vortices, symmetric;
    nwake=size.(wakes, 1), xhat = SVector(1, 0, 0))

    TF = promote_type(eltype(eltype(eltype(surface))), eltype(eltype(eltype(wake))))
    N = sum(length.(surfaces))
    b = zeros(TF, N)

    return add_wake_normal_velocity!(b, surfaces, wakes, surface_id, trailing_vortices,
        symmetric; nwake=nwake, xhat=xhat)
end

"""
    add_wake_normal_velocity!(b, surface, wake, trailing_vortices, symmetric; kwargs...)

Pre-allocated version of `wake_normal_velocity` which adds the normal induced
velocity created by the wake to the existing vector `b`
"""
function add_wake_normal_velocity!(b, surface::AbstractMatrix, wake::AbstractMatrix,
    trailing_vortices, symmetric; nwake=size(wake,1), xhat=SVector(1,0,0))

    same_id = true

    return add_wake_normal_velocity!(b, surface, wake, symmetric, same_id, nwake,
        trailing_vortices, xhat)
end

"""
    add_wake_normal_velocity!(b, surfaces, wakes, surface_id, trailing_vortices, symmetric; kwargs...)

Pre-allocated version of `wake_normal_velocity` which adds the normal induced
velocity created by the wakes to the existing vector `b`
"""
function add_wake_normal_velocity!(b, surfaces::AbstractVector{<:AbstractMatrix},
    wakes::AbstractVector{<:Wake}, surface_id, trailing_vortices, symmetric;
    nwake=size.(wakes,1), xhat=SVector(1, 0, 0))

    nsurf = length(surfaces)
    nwake = length(wakes)

    # index for keeping track of where we are in the b vector
    ib = 0

    # loop through receiving surfaces
    for i = 1:nsurf

        # number of panels on this surface
        n = length(surfaces[i])

        # create view into RHS vector
        vb = view(b, ib+1:ib+n)

        # fill in RHS vector
        for j = 1:nsurf
            same_id = surface_id[i] == surface_id[j]
            add_wake_normal_velocity!(vb, surfaces[i], wakes[j], nwake[j], symmetric, same_id, trailing_vortices, xhat)
        end

        # increment position in AIC matrix
        ib += n
    end

    return b
end

"""
    add_wake_normal_velocity!(b, surface, wake, symmetric, same_id,
        trailing_vortices; xhat=SVector(1, 0, 0), nwake=size(wake, 1))

Adds the normal induced velocity from `wake` onto the control points of `surface`
to the existing vector `b`
"""
@inline function add_wake_normal_velocity!(b, surface, wake, symmetric, same_id,
    trailing_vortices; xhat=SVector(1, 0, 0), nwake=size(wake, 1))

    # loop over receiving panels
    for i = 1:length(surface)

        # control point location
        rcp = controlpoint(receiving[i])

        # normal vector body axis
        nhat = normal(receiving[i])

        # get induced velocity at rcp from the wake
        Vind = wake_induced_velocity(rcp, wake, symmetric, same_id, trailing_vortices;
            xhat = xhat, nwake = nwake)

        # add normal velocity
        b[i] += dot(Vind, nhat)
    end

    return b
end

# --- circulation solve --- #

"""
    circulation(AIC, b)

Solve for the circulation distribution.
"""
circulation(AIC, b) = AIC\b

circulation!(Γ, AIC, b) = ldiv!(Γ, lu(AIC), b)

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

function circulation_derivatives!(Γ, dΓ, AIC, b, db)

    # unpack derivatives
    (Γ_a, Γ_b, Γ_p, Γ_q, Γ_r) = dΓ
    (b_a, b_b, b_p, b_q, b_r) = db

    # factorize AIC matrix (since we'll be reusing it)
    fAIC = lu(AIC)

    # solve for circulation and its derivatives
    ldiv!(Γ, fAIC, b)

    ldiv!(Γ_a, fAIC, b_a)
    ldiv!(Γ_b, fAIC, b_b)
    ldiv!(Γ_p, fAIC, b_p)
    ldiv!(Γ_q, fAIC, b_q)
    ldiv!(Γ_r, fAIC, b_r)

    return Γ, dΓ
end
