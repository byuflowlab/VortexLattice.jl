"""
    StabilityDerivatives

Holds stability derivatives

# Fields:
 - `alpha`: angle of attack
 - `beta`: sideslip
 - `p`: roll rate
 - `q`: pitch rate
 - `r`: yaw rate
"""
struct StabilityDerivatives{TF}
    alpha::TF
    beta::TF
    p::TF
    q::TF
    r::TF
end

function StabilityDerivatives(alpha, beta, p, q, r)
    TF = promote_type(typeof(alpha), typeof(beta), typeof(p), typeof(q), typeof(r))
    return StabilityDerivatives{TF}(alpha, beta, p, q, r)
end

Base.eltype(::Type{StabilityDerivatives{TF}}) where TF = TF
Base.eltype(::StabilityDerivatives{TF}) where TF = TF

# -------------------------------

function stability_analysis(panels, reference, freestream, symmetric)

    # make sure panels is of concrete type
    TF = promote_type(eltype.(panels)...)
    panels = Vector{Panel{TF}}(panels)

    x = vcat(freestream.alpha, freestream.beta, freestream.Omega)

    f = function(x)
        freestream = Freestream(x[1], x[2], SVector(x[3], x[4], x[5]))
        Γ = circulation(panels, reference, freestream, symmetric)
        CF, CM, _ = forces_moments(panels, reference, freestream, Γ, symmetric)
        return vcat(CF, CM)
    end

    dfdx = ForwardDiff.jacobian(f, x)

    iCF = SVector(1,2,3)
    CF_α = dfdx[iCF,1]
    CF_β = dfdx[iCF,2]
    CF_p = dfdx[iCF,3]*2*VINF/reference.b
    CF_q = dfdx[iCF,4]*2*VINF/reference.c
    CF_r = dfdx[iCF,5]*2*VINF/reference.b
    dCF = StabilityDerivatives(CF_α, CF_β, CF_p, CF_q, CF_r)

    iCM = SVector(4,5,6)
    CM_α = dfdx[iCM,1]
    CM_β = dfdx[iCM,2]
    CM_p = dfdx[iCM,3]*2*VINF/reference.b
    CM_q = dfdx[iCM,4]*2*VINF/reference.c
    CM_r = dfdx[iCM,5]*2*VINF/reference.b
    dCM = StabilityDerivatives(CM_α, CM_β, CM_p, CM_q, CM_r)

    return dCF, dCM

end
