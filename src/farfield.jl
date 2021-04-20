"""
    far_field_drag(trefftz, surfaces, ref, fs, symmetric, Γ)

Compute induced drag in the Trefftz plane (far field method).
"""
function far_field_drag(trefftz, surfaces, ref, fs, symmetric, Γ)
    # construct trefftz panels
    trefftz_panels!(trefftz, surfaces, fs, Γ)
    # perform far field analysis
    nsurf = length(surfaces)
    CD = zero(eltype(eltype(eltype(trefftz))))
    for i = 1:nsurf, j = 1:nsurf
        CD += far_field_drag(trefftz[i], trefftz[j], ref, symmetric[j])
    end
    return CD
end

"""
    far_field_drag(receiving, sending, reference, symmetric)

Computes the induced drag on `receiving` from `sending` using the Trefftz
plane analysis.

# Arguments
 - `receiving`: Vector of receiving Trefftz panels (see [`TrefftzPanel`](@ref))
 - `sending`: Vector of sending Trefftz panels (see [`TrefftzPanel`](@ref))
 - `reference`: Reference parameters (see [`Reference`](@ref))
 - `symmetric`: Flag indicating whether a mirror image of the panels in `surface`,
    should be used when calculating induced velocities
"""
@inline function far_field_drag(receiving, sending, ref, symmetric)
    # get float type
    TF = promote_type(eltype(eltype(receiving)), eltype(eltype(sending)), eltype(ref))
    # get number of receiving and sending panels
    Nr = length(receiving)
    Ns = length(sending)
    # add up drag
    Di = zero(TF)
    for j = 1:Ns, i = 1:Nr
        Di += trefftz_panel_induced_drag(receiving[i], sending[j]; symmetric)
    end
    # apply symmetry
    if symmetric
        Di *= 2
    end
    # reference dynamic pressure
    q = 1/2*RHO*ref.V^2
    # normalize
    CDi = Di / (q*ref.S)
    # return result
    return CDi
end
