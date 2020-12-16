

"""
    far_field_drag(system, surface, reference, freestream; kwargs...)

Computes induced drag using the Trefftz plane (far field method).

Note that this function assumes that the circulation distribution has already
been computed and is present in `system`

# Arguments
 - `system`: Pre-allocated system properties
 - `surface`: Matrix of panels of shape (nc, ns) where `nc` is the number of
    chordwise panels and `ns` is the number of spanwise panels
 - `reference`: Reference parameters (see `Reference`)
 - `freestream`: Freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: Flag indicating whether a mirror image of the panels in `surface`,
    should be used when calculating induced velocities
"""
function far_field_drag(system, surface::AbstractMatrix, ref, fs; symmetric)

    # unpack system
    Γ = system.gamma
    trefftz = system.trefftz[1]

    trefftz_panels!(trefftz, surface, fs, Γ)

    return far_field_drag(trefftz, trefftz, ref, fs; symmetric)
end

"""
    far_field_drag(system, surfaces, reference, freestream; kwargs...)

Computes induced drag using the Trefftz plane (far field method).

Note that this function assumes that the circulation distribution has already
been computed and is present in `system`

# Arguments
 - `system`: Pre-allocated system properties
 - `surfaces`: Vector of surfaces, represented by matrices of panels of shape
    (nc, ns) where `nc` is the number of chordwise panels and `ns` is the number
    of spanwise panels
 - `reference`: Reference parameters (see `Reference`)
 - `freestream`: Freestream parameters (see `Freestream`)

# Keyword Arguments
 - `symmetric`: Flag indicating whether a mirror image of the panels in `surface`,
    should be used when calculating induced velocities
"""
function far_field_drag(system, surfaces::AbstractVector{<:AbstractMatrix}, ref,
    fs; symmetric)

    # unpack system
    Γ = system.gamma
    trefftz = system.trefftz

    # construct trefftz panels
    trefftz_panels!(trefftz, surfaces, fs, Γ)

    # perform far field analysis
    nsurf = length(surfaces)
    CD = zero(eltype(eltype(eltype(surfaces))))
    for i = 1:nsurf, j = 1:nsurf
        CD += far_field_drag(trefftz[i], trefftz[j], ref, fs; symmetric = symmetric[j])
    end

    return CD
end

# one surface on another surface
function far_field_drag(receiving, sending, ref, fs; symmetric)

    TF = promote_type(eltype(eltype(receiving)), eltype(eltype(sending)), eltype(ref), eltype(fs))

    Nr = length(receiving)
    Ns = length(sending)

    # add up drag
    Di = zero(TF)
    for j = 1:Ns
        for i = 1:Nr
            Di += trefftz_panel_induced_drag(receiving[i], sending[j]; symmetric)
        end
    end

    # apply symmetry
    if symmetric
        Di *= 2
    end

    # normalize
    CDi = Di / (QINF*ref.S)

    return CDi
end
