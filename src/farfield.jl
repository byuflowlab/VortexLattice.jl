"""
    far_field_drag(system)

Computes induced drag using the Trefftz plane (far field method).

Note that this function assumes that the circulation distribution has already
been computed and is present in `system`

# Arguments
 - `system`: Pre-allocated system properties
"""
function far_field_drag(system)

    # unpack system
    surfaces = system.surfaces
    trefftz = system.trefftz
    ref = system.reference[]
    fs = system.freestream[]
    symmetric = system.symmetric
    Γ = system.Γ

    # construct trefftz panels
    trefftz_panels!(trefftz, surfaces, fs, Γ)

    # perform far field analysis
    nsurf = length(surfaces)
    CD = zero(eltype(system))
    for i = 1:nsurf, j = 1:nsurf
        CD += far_field_drag(trefftz[i], trefftz[j], ref, symmetric[j])
    end

    return CD
end

# one surface on another surface
"""
    far_field_drag(receiving, sending, reference, freestream; kwargs...)

Computes the induced drag on `receiving` from `sending` using the Trefftz
plane analysis.

# Arguments
 - `receiving`: Vector of receiving Trefftz panels (see [`TrefftzPanel`](@ref))
 - `sending`: Vector of sending Trefftz panels (see [`TrefftzPanel`](@ref))
 - `reference`: Reference parameters (see [`Reference`](@ref))

# Keyword Arguments
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

    # normalize
    CDi = Di / (QINF*ref.S)

    return CDi
end
