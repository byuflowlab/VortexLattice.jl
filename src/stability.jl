"""
    body_derivatives(system, surface, reference, freestream; kwargs...)

Returns the derivatives of the body forces and moments with respect to the
freestream velocity components `(u, v, w)` and the angular velocity components
`(p, q, r)` in the body frame.

The derivatives are returned as two named tuples: `dCF, dCM`

Note that the derivatives with respect to the freestream variables of the panel
forces must have been previously computed and stored in `system`

# Arguments:
 - `system`: Object of type `System` which holds system properties
 - `surfaces`: Surface, represented by matrices of panels of shape (nc, ns) where
    `nc` is the number of chordwise panels and `ns` is the number of spanwise panels
 - `reference`: reference parameters (see `Reference`)
 - `freestream`: freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image of the panels
    in `surface` was used when calculating induced velocities.
"""
function body_derivatives(system, surface::AbstractMatrix, ref, fs; symmetric)

    return body_derivatives(system, [surface], ref, fs; symmetric = [symmetric])
end

"""
    body_derivatives(system, surfaces, reference, freestream; kwargs...)

Returns the derivatives of the body forces and moments with respect to the
freestream velocity components `(u, v, w)` and the angular velocity components
`(p, q, r)` in the body frame.

The derivatives are returned as two named tuples: `dCF, dCM`

Note that the derivatives with respect to the freestream variables of the panel
forces must have been previously computed and stored in `system`

# Arguments:
 - `system`: Object of type `System` which holds system properties
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: reference parameters (see `Reference`)
 - `freestream`: freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image of the panels in `surface`
    was used when calculating panel induced velocities.
"""
function body_derivatives(system, surfaces::AbstractVector{<:AbstractMatrix},
    ref, fs; symmetric)

    CF, CM, dCF, dCM = body_forces_derivatives(system, surfaces, ref, fs; symmetric)

    # unpack derivatives
    (CF_a, CF_b, CF_p, CF_q, CF_r) = dCF
    (CM_a, CM_b, CM_p, CM_q, CM_r) = dCM

    # convert derivatives to be with respect to (u, v, w) instead of (alpha, beta)
    u, v, w = freestream_velocity(fs)
    u2, v2, w2 = u^2, v^2, w^2
    # a = atan(w, u)
    a_u = -w/(u2+w2)
    a_w = u/(u2+w2)
    # b = atan(v, u)
    b_u =-v/(u2+v2)
    b_v = u/(u2+v2)

    CF_u = CF_a*a_u + CF_b*b_u
    CF_v = CF_b*b_v
    CF_w = CF_a*a_w

    CM_u = CM_a*a_u + CM_b*b_u
    CM_v = CM_b*b_v
    CM_w = CM_a*a_w

    # pack up derivatives as named tuples
    dCF = (u=CF_u, v=CF_v, w=CF_w, p=CF_p, q=CF_q, r=CF_r)
    dCM = (u=CM_u, v=CM_v, w=CM_w, p=CM_p, q=CM_q, r=CM_r)

    return dCF, dCM
end

"""
    stability_derivatives(system, surfaces, reference, freestream; kwargs...)

Returns the derivatives of the body forces and moments in the stability frame
with respect to the freestream velocity components `(alpha, beta)` and the angular
velocity components `(p, q, r)` in the stability frame.

The derivatives are returned as two named tuples: `dCF, dCM`

Note that the derivatives with respect to the freestream variables of the panel
forces must have been previously computed and stored in `system`

# Arguments:
 - `system`: Object of type `System` which holds system properties
 - `surfaces`: Surface, represented by matrices of panels of shape (nc, ns) where
    `nc` is the number of chordwise panels and `ns` is the number of spanwise panels
 - `reference`: reference parameters (see `Reference`)
 - `freestream`: freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image of the panels
    in `surface` was used when calculating panel induced velocities.
"""
function stability_derivatives(system, surface::AbstractMatrix, ref, fs;
    symmetric = false)

    return stability_derivatives(system, [surface], ref, fs;
        symmetric = [symmetric])
end

"""
    stability_derivatives(system, surfaces, reference, freestream; kwargs...)

Returns the derivatives of the body forces and moments in the stability frame
with respect to the freestream velocity components `(alpha, beta)` and the angular
velocity components `(p, q, r)` in the stability frame.

The derivatives are returned as two named tuples: `dCF, dCM`

Note that the derivatives with respect to the freestream variables of the panel
forces must have been previously computed and stored in `system`

# Arguments:
 - `system`: Object of type `System` which holds system properties
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: reference parameters (see `Reference`)
 - `freestream`: freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: (required) Flag indicating whether a mirror image of the panels
    in `surface` was used when calculating panel induced velocities.
"""
function stability_derivatives(system, surfaces::AbstractVector{<:AbstractMatrix},
    ref, fs; symmetric)

    CFb, CMb, dCFb, dCMb = body_forces_derivatives(system, surfaces, ref, fs; symmetric)

    # unpack derivatives
    (CFb_a, CFb_b, CFb_pb, CFb_qb, CFb_rb) = dCFb
    (CMb_a, CMb_b, CMb_pb, CMb_qb, CMb_rb) = dCMb

    # rotate forces/moments into the stability frame
    R, R_a = body_to_stability_alpha(fs)

    CFs = R * CFb
    CFs_a = R * CFb_a + R_a * CFb
    CFs_b = R * CFb_b
    CFs_pb = R * CFb_pb
    CFs_qb = R * CFb_qb
    CFs_rb = R * CFb_rb

    CMs = R * CMb
    CMs_a = R * CMb_a + R_a * CMb
    CMs_b = R * CMb_b
    CMs_pb = R * CMb_pb
    CMs_qb = R * CMb_qb
    CMs_rb = R * CMb_rb

    # switch (p,q,r) derivatives to be wrt (p,q,r) of stability frame
    Ωb = fs.Omega # body p, q, r
    Ωs = R * Ωb # stability p, q, r

    # Ωb = R' * Ωs
    pb_ps, qb_ps, rb_ps = R[1,:]
    pb_qs, qb_qs, rb_qs = R[2,:]
    pb_rs, qb_rs, rb_rs = R[3,:]

    CFs_ps = CFs_pb*pb_ps + CFs_qb*qb_ps + CFs_rb*rb_ps
    CFs_qs = CFs_pb*pb_qs + CFs_qb*qb_qs + CFs_rb*rb_qs
    CFs_rs = CFs_pb*pb_rs + CFs_qb*qb_rs + CFs_rb*rb_rs

    CMs_ps = CMs_pb*pb_ps + CMs_qb*qb_ps + CMs_rb*rb_ps
    CMs_qs = CMs_pb*pb_qs + CMs_qb*qb_qs + CMs_rb*rb_qs
    CMs_rs = CMs_pb*pb_rs + CMs_qb*qb_rs + CMs_rb*rb_rs

    # assign outputs, and apply stability derivative normalizations
    CF_a = CFs_a
    CF_b = CFs_b
    CF_p = CFs_ps*2*VINF/ref.b
    CF_q = CFs_qs*2*VINF/ref.c
    CF_r = CFs_rs*2*VINF/ref.b

    CM_a = CMs_a
    CM_b = CMs_b
    CM_p = CMs_ps*2*VINF/ref.b
    CM_q = CMs_qs*2*VINF/ref.c
    CM_r = CMs_rs*2*VINF/ref.b

    # pack up derivatives as named tuples
    dCF = (alpha=CF_a, beta=CF_b, p=CF_p, q=CF_q, r=CF_r)
    dCM = (alpha=CM_a, beta=CM_b, p=CM_p, q=CM_q, r=CM_r)

    return dCF, dCM
end
