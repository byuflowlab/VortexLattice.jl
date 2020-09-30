"""
    Freestream(alpha, beta, Omega, additional_velocity=nothing)

Define the freestream properties.

**Arguments**
- `alpha`: angle of attack (rad)
- `beta`: sideslip angle (rad)
- `Omega`: rotation vector (p, q, r) of the body frame about the center of
    gravity, normalized by Vinf
- `additional_velocity`: a function of the form: V = additional_velocity(r) which returns
    the additional velocity `V` (normalized by the freestream velocity) at
    position `r`.  Defaults to `nothing`.
"""
struct Freestream{TF, TV}
    alpha::TF
    beta::TF
    Omega::SVector{3, TF}
    additional_velocity::TV
end

function Freestream(alpha, beta, Omega, additional_velocity = nothing)
    TF = promote_type(typeof(alpha), typeof(beta), eltype(Omega))
    return Freestream{TF, typeof(additional_velocity)}(alpha, beta, Omega, additional_velocity)
end

Base.eltype(::Type{Freestream{TF, TV}}) where {TF,TV} = TF
Base.eltype(::Freestream{TF, TV}) where {TF,TV} = TF

"""
    body_to_wind(fs::Freestream)

Construct a rotation matrix from the body axis to the wind axis
"""
@inline function body_to_wind(fs::Freestream)

    alpha = fs.alpha
    beta = fs.beta

    cb, sb = cos(beta), sin(beta)
    ca, sa = cos(alpha), sin(alpha)

    Rb = @SMatrix [cb -sb 0; sb cb 0; 0 0 1]
    Ra = @SMatrix [ca 0 sa; 0 1 0; -sa 0 ca]
    R = Rb*Ra

    return R
end

"""
    external_velocity(freestream, r, rcg)

Compute the external velocity at location `r`
"""
@inline function external_velocity(freestream, r, rcg)

    # Freestream velocity in body coordinate
    ca, sa = cos(freestream.alpha), sin(freestream.alpha)
    cb, sb = cos(freestream.beta), sin(freestream.beta)

    Vext = VINF*SVector(ca*cb, -sb, sa*cb)

    # add velocity due to body rotation
    Vext -= VINF*cross(freestream.Omega, r - rcg)  # unnormalize

    # add contribution due to additional velocity field
    if !isnothing(freestream.additional_velocity)
        Vext += VINF*freestream.additional_velocity(r)  # unnormalize
    end

    return Vext
end
