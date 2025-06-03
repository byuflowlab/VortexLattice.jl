"""
    Freestream([Vinf,] alpha, beta, Omega)

Defines the freestream and rotational velocity properties.

**Arguments**
- `Vinf`: Freestream velocity
- `alpha`: angle of attack (rad)
- `beta`: sideslip angle (rad)
- `Omega`: rotation vector (p, q, r) of the body frame about the reference center.
    Uses standard coordinate system from dynamics (positve p roll right wing down to turn right, positive q is pitch nose up, positive r is yaw nose to the right)
"""
struct Freestream{TF}
    Vinf::TF
    alpha::TF
    beta::TF
    Omega::SVector{3, TF}
end

function Freestream(Vinf, alpha, beta, Omega)
    TF = promote_type(typeof(Vinf), typeof(alpha), typeof(beta), eltype(Omega))
    return Freestream{TF}(Vinf, alpha, beta, Omega)
end

Base.eltype(::Type{Freestream{TF}}) where TF = TF
Base.eltype(::Freestream{TF}) where TF = TF

"""
    body_to_stability(freestream)

Construct a rotation matrix from the body axis to the stability axis.
"""
function body_to_stability(fs::Freestream)
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

function body_to_stability_alpha(fs::Freestream)
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
function stability_to_body(fs::Freestream)
    sa, ca = sincos(fs.alpha)
    return stability_to_body(sa, ca)
end

stability_to_body(sa, ca) = body_to_stability(sa, ca)'

"""
    stability_to_body(freestream)

Construct a rotation matrix from the stability axis to the body axis and its
derivative with respect to `alpha`
"""
function stability_to_body_alpha(fs::Freestream)
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

function stability_to_wind(fs::Freestream)
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

function stability_to_wind_beta(fs::Freestream)
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
function wind_to_stability(fs::Freestream)
    sb, cb = sincos(fs.beta)
    return wind_to_stability(sb, cb)
end

wind_to_stability(sb, cb) = stability_to_wind(sb, cb)'

"""
    wind_to_stability_beta(freestream)

Construct a rotation matrix from the wind axis to the stability axis and its
derivative with respect to `beta`
"""
function wind_to_stability_beta(fs::Freestream)
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
function body_to_wind(fs::Freestream)
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
function body_to_wind_derivatives(fs::Freestream)
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
function wind_to_body(fs::Freestream)
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
function wind_to_body_derivatives(fs::Freestream)
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

Given a freestream velocity vector `V` in body x-y-z coordinates (South-East-Up convention), returns a `Freestream`
object with computed `Vinf`, `alpha`, and `beta`. Optionally, specify `Omega` (default zeros).
"""
function velocity_to_freestream(V::SVector{3,TF}, Omega=zero(SVector{3,TF})) where TF
    Vinf = norm(V)
    if iszero(Vinf)
        return Freestream(0.0, 0.0, 0, Omega)
    end
    α = atan(V[3], V[1])
    β = -asin(V[2] / Vinf)
    
    return Freestream(Vinf, α, β, Omega)
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
    # println("NON-deriv")
    # @show Vinf*SVector(ca*cb, -sb, sa*cb)
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

    # println("Deriv:")
    # @show V

    return V, dV
end

"""
    rotational_velocity(r, freestream, reference)

Compute the velocity due to body rotations about the reference center
"""
rotational_velocity

rotational_velocity(r, fs::Freestream, ref::Reference) = rotational_velocity(r, fs.Omega, ref.r)

function rotational_velocity(r, Omega, rref) 
    Ω = [-Omega[1], Omega[2], -Omega[3]]  # swap signs for p and r to follow standard dynamics convention
    return cross(r - rref, Ω)
end

"""
    rotational_velocity_derivatives(r, freestream, reference)

Compute the velocity due to body rotations about the reference center and its
derivatives with respect to (p, q, r)
"""
rotational_velocity_derivatives

rotational_velocity_derivatives(r, fs::Freestream, ref::Reference) = rotational_velocity_derivatives(r, fs.Omega, ref.r)

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

"""
    trajectory_to_freestream(dt; kwargs...)

Convert trajectory parameters into freestream velocity parameters (see [`Freestream`](@ref))
at a collection of time steps.

# Arguments:
 - `dt`: Time step vector (seconds)

# Keyword Arguments:
 - `Xdot = zeros(length(dt))`: Global frame x-velocity for each time step
 - `Ydot = zeros(length(dt))`: Global frame y-velocity for each time step
 - `Zdot = zeros(length(dt))`: Global frame z-velocity for each time step
 - `p = zeros(length(dt))`: Angular velocity about x-axis for each time step
 - `q = zeros(length(dt))`: Angular velocity about y-axis for each time step
 - `r = zeros(length(dt))`: Angular velocity about z-axis for each time step
 - `phi0 = 0`: Roll angle for initial time step
 - `theta0 = 0`: Pitch angle for initial time step
 - `psi0 = 0`: Yaw angle for initial time step
"""
function trajectory_to_freestream(dt;
    Xdot = zeros(length(dt)), Ydot = zeros(length(dt)), Zdot = zeros(length(dt)),
    p = zeros(length(dt)), q = zeros(length(dt)), r = zeros(length(dt)),
    phi0 = 0, theta0 = 0, psi0 = 0)

    # change scalar inputs to vectors
    if isa(Xdot, Number)
        Xdot = fill(Xdot, length(dt))
    end

    if isa(Ydot, Number)
        Ydot = fill(Ydot, length(dt))
    end

    if isa(Zdot, Number)
        Zdot = fill(Zdot, length(dt))
    end

    if isa(p, Number)
        p = fill(p, length(dt))
    end

    if isa(q, Number)
        q = fill(q, length(dt))
    end

    if isa(r, Number)
        r = fill(r, length(dt))
    end

    # get common floating point type
    TF = promote_type(eltype(dt), eltype(Xdot), eltype(Ydot), eltype(Zdot),
        eltype(p), eltype(q), eltype(r), typeof(phi0), typeof(theta0), typeof(psi0))

    # number of discrete time steps
    nt = length(dt)

    # freestream parameters corresponding to each time step
    fs = Vector{Freestream{TF}}(undef, nt)

    # set initial orientation
    ϕ = phi0
    θ = theta0
    ψ = psi0

    # populate freestream parameters for each time step
    for it = 1:length(dt)

        sϕ, cϕ = sincos(ϕ)
        sθ, cθ = sincos(θ)
        sψ, cψ = sincos(ψ)

        Rϕ = @SMatrix [1 0 0; 0 cϕ sϕ; 0 -sϕ cϕ]
        Rθ = @SMatrix [cθ 0 -sθ; 0 1 0; sθ 0 cθ]
        Rψ = @SMatrix [cψ sψ 0; -sψ cψ 0; 0 0 1]

        # freestream velocity vector
        V = Rϕ*Rθ*Rψ*SVector(-Xdot[it], -Ydot[it], -Zdot[it])

        # convert to angular representation
        Vinf = norm(V)
        α = atan(V[3]/V[1])
        β = -asin(V[2]/norm(V))
        Ω = SVector(-p[it], q[it], -r[it])  # using standard dynamic convention

        # assemble freestream parameters
        fs[it] = Freestream(Vinf, α, β, Ω)

        if it < length(dt)
            # update orientation
            ϕ += p[it]*dt[it]
            θ += q[it]*dt[it]
            ψ += r[it]*dt[it]
        end
    end

    return fs
end
