"""
    BodyFreestream(Vinf, alpha, beta, Omega)

Internal: the relative wind and rotation rates in body angles, the variables the
stability derivatives are taken with respect to. Built from the core's
`Freestream` and the frames by [`body_freestream`](@ref); users never construct it.

**Arguments**
- `Vinf`: freestream speed
- `alpha`: angle of attack (rad)
- `beta`: sideslip angle (rad)
- `Omega`: rotation vector (p, q, r) of the body frame about the reference center.
    Uses standard coordinate system from dynamics (positve p roll right wing down to turn right, positive q is pitch nose up, positive r is yaw nose to the right)
"""
struct BodyFreestream{TF}
    Vinf::TF
    alpha::TF
    beta::TF
    Omega::SVector{3, TF}
end

function BodyFreestream(Vinf, alpha, beta, Omega)
    TF = promote_type(typeof(Vinf), typeof(alpha), typeof(beta), eltype(Omega))
    return BodyFreestream{TF}(Vinf, alpha, beta, Omega)
end

Base.eltype(::Type{BodyFreestream{TF}}) where TF = TF
Base.eltype(::BodyFreestream{TF}) where TF = TF

"""
    body_to_stability(freestream)

Construct a rotation matrix from the body axis to the stability axis.
"""
function body_to_stability(fs::BodyFreestream)
    sa, ca = sincos(fs.alpha)
    return body_to_stability(sa, ca)
end

body_to_stability(sa, ca) = @SMatrix [ca 0 sa; 0 1 0; -sa 0 ca]

"""
    body_to_stability_alpha(freestream)

Construct a rotation matrix from the body axis to the stability axis and its
derivative with respect to `alpha`
"""
body_to_stability_alpha

function body_to_stability_alpha(fs::BodyFreestream)
    sa, ca = sincos(fs.alpha)
    return body_to_stability_alpha(sa, ca)
end

function body_to_stability_alpha(sa, ca)

    R = body_to_stability(sa, ca)

    R_a = @SMatrix [-sa 0 ca; 0 0 0; -ca 0 -sa]

    return R, R_a
end

"""
    stability_to_body(freestream)

Construct a rotation matrix from the stability axis to the body axis
"""
function stability_to_body(fs::BodyFreestream)
    sa, ca = sincos(fs.alpha)
    return stability_to_body(sa, ca)
end

stability_to_body(sa, ca) = body_to_stability(sa, ca)'

"""
    stability_to_body(freestream)

Construct a rotation matrix from the stability axis to the body axis and its
derivative with respect to `alpha`
"""
function stability_to_body_alpha(fs::BodyFreestream)
    sa, ca = sincos(fs.alpha)
    return stability_to_body_alpha(sa, ca)
end

function stability_to_body_alpha(sa, ca)

    R, R_a = body_to_stability_alpha(sa, ca)

    return R', R_a'
end

"""
    stability_to_wind(freestream)

Construct a rotation matrix from the stability axis to the wind axis
"""
stability_to_wind

function stability_to_wind(fs::BodyFreestream)
    sb, cb = sincos(fs.beta)
    return stability_to_wind(sb, cb)
end

stability_to_wind(sb, cb) = @SMatrix [cb -sb 0; sb cb 0; 0 0 1]

"""
    stability_to_wind_beta(freestream)

Construct a rotation matrix from the stability axis to the wind axis and its
derivative with respect to `beta`
"""
stability_to_wind_beta

function stability_to_wind_beta(fs::BodyFreestream)
    sb, cb = sincos(fs.beta)
    return stability_to_wind_beta(sb, cb)
end

function stability_to_wind_beta(sb, cb)

    R = stability_to_wind(sb, cb)

    R_b = @SMatrix [-sb -cb 0; cb -sb 0; 0 0 0]

    return R, R_b
end

"""
    wind_to_stability(freestream)

Construct a rotation matrix from the wind axis to the stability axis
"""
function wind_to_stability(fs::BodyFreestream)
    sb, cb = sincos(fs.beta)
    return wind_to_stability(sb, cb)
end

wind_to_stability(sb, cb) = stability_to_wind(sb, cb)'

"""
    wind_to_stability_beta(freestream)

Construct a rotation matrix from the wind axis to the stability axis and its
derivative with respect to `beta`
"""
function wind_to_stability_beta(fs::BodyFreestream)
    sb, cb = sincos(fs.beta)
    return wind_to_stability_beta(sb, cb)
end

function wind_to_stability_beta(sb, cb)

    R, R_b = stability_to_wind_beta(sb, cb)

    return R', R_b'
end

"""
    body_to_wind(freestream)

Construct a rotation matrix from the body axis to the wind axis
"""
function body_to_wind(fs::BodyFreestream)
    sa, ca = sincos(fs.alpha)
    sb, cb = sincos(fs.beta)
    return body_to_wind(sa, ca, sb, cb)
end

function body_to_wind(sa, ca, sb, cb)

    Ra = body_to_stability(sa, ca)

    Rb = stability_to_wind(sb, cb)

    return Rb*Ra
end

"""
    body_to_wind_derivatives(freestream)

Construct a rotation matrix from the body axis to the wind axis and its
derivatives with respect to `alpha` and `beta`
"""
function body_to_wind_derivatives(fs::BodyFreestream)
    sa, ca = sincos(fs.alpha)
    sb, cb = sincos(fs.beta)
    return body_to_wind_derivatives(sa, ca, sb, cb)
end

function body_to_wind_derivatives(sa, ca, sb, cb)

    Ra, Ra_a = body_to_stability_alpha(sa, ca)

    Rb, Rb_b = stability_to_wind_beta(sb, cb)

    R = Rb*Ra
    R_a = Rb*Ra_a
    R_b = Rb_b*Ra

    return R, R_a, R_b
end

"""
    wind_to_body(freestream)

Construct a rotation matrix from the wind axis to the body axis
"""
function wind_to_body(fs::BodyFreestream)
    sa, ca = sincos(fs.alpha)
    sb, cb = sincos(fs.beta)
    return wind_to_body(sa, ca, sb, cb)
end

wind_to_body(sa, ca, sb, cb) = body_to_wind(sa, ca, sb, cb)'

"""
    wind_to_body_derivatives(freestream)

Construct a rotation matrix from the wind axis to the body axis and its derivatives
with respect to `alpha` and `beta`
"""
function wind_to_body_derivatives(fs::BodyFreestream)
    sa, ca = sincos(fs.alpha)
    sb, cb = sincos(fs.beta)
    return wind_to_body_derivatives(sa, ca, sb, cb)
end

function wind_to_body_derivatives(sa, ca, sb, cb)

    R, R_a, R_b = body_to_wind_derivatives(sa, ca, sb, cb)

    return R', R_a', R_b'
end

"""
    velocity_to_freestream(V::SVector{3}, Omega::SVector{3}=zeros(3))

Given a freestream velocity vector `V` in body x-y-z coordinates (South-East-Up convention), returns a `BodyFreestream`
object with computed `Vinf`, `alpha`, and `beta`. Optionally, specify `Omega` (default zeros).
"""
function velocity_to_freestream(V::SVector{3,TF}, Omega=zero(SVector{3,TF})) where TF
    Vinf = norm(V)
    if iszero(Vinf)
        return BodyFreestream(0.0, 0.0, 0, Omega)
    end
    α = atan(V[3], V[1])
    β = -asin(V[2] / Vinf)
    
    return BodyFreestream(Vinf, α, β, Omega)
end

"""
    freestream_velocity(freestream)

Computes the freestream velocity
"""
function freestream_velocity(fs)

    Vinf = fs.Vinf
    sa, ca = sincos(fs.alpha)
    sb, cb = sincos(fs.beta)

    return freestream_velocity(Vinf, sa, ca, sb, cb)
end

function freestream_velocity(Vinf, sa, ca, sb, cb) 
    return Vinf*SVector(ca*cb, -sb, sa*cb)
end

"""
    freestream_velocity_derivatives(freestream)

Computes the freestream velocity
"""
function freestream_velocity_derivatives(fs)

    Vinf = fs.Vinf
    sa, ca = sincos(fs.alpha)
    sb, cb = sincos(fs.beta)

    return freestream_velocity_derivatives(Vinf, sa, ca, sb, cb)
end

function freestream_velocity_derivatives(Vinf, sa, ca, sb, cb)

    V = Vinf*SVector(ca*cb, -sb, sa*cb)

    V_a = Vinf*SVector(-sa*cb, 0, ca*cb)
    V_b = Vinf*SVector(-ca*sb, -cb, -sa*sb)

    dV = (V_a, V_b)


    return V, dV
end

"""
    rotational_velocity(r, freestream, reference)

Compute the velocity due to body rotations about the reference center
"""
rotational_velocity

rotational_velocity(r, fs::BodyFreestream, ref::Reference) = rotational_velocity(r, fs.Omega, ref.r)

function rotational_velocity(r, Omega, rref) 
    tmp = r - rref
    # swap signs for p and r to follow standard dynamics convention
    return SVector(tmp[2]*(-Omega[3]) - tmp[3]*Omega[2], tmp[3]*(-Omega[1]) - tmp[1]*(-Omega[3]), tmp[1]*Omega[2] - tmp[2]*(-Omega[1]))
end

"""
    rotational_velocity_derivatives(r, freestream, reference)

Compute the velocity due to body rotations about the reference center and its
derivatives with respect to (p, q, r)
"""
rotational_velocity_derivatives

rotational_velocity_derivatives(r, fs::BodyFreestream, ref::Reference) = rotational_velocity_derivatives(r, fs.Omega, ref.r)

function rotational_velocity_derivatives(r, Omega, rref)

    tmp = r - rref
    Ω = [-Omega[1], Omega[2], -Omega[3]]  # swap signs for p and r to follow standard dynamics convention

    Vrot = cross(tmp, Ω)
    Vrot_p = -SVector(0, tmp[3], -tmp[2])
    Vrot_q = SVector(-tmp[3], 0, tmp[1])
    Vrot_r = -SVector(tmp[2], -tmp[1], 0)

    dVrot = (Vrot_p, Vrot_q, Vrot_r)

    return Vrot, dVrot
end
