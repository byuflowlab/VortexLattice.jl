# --- frames of reference for expressing near-field forces --- #

"""
    AbstractFrame

Supertype for the different possible reference frames used by this package.
"""
abstract type AbstractFrame end

"""
   Body <: AbstractFrame

Reference frame aligned with the global X-Y-Z axes
"""
struct Body <: AbstractFrame end

"""
    Stability <: AbstractFrame

Reference frame rotated from the body frame about the y-axis to be aligned with
the freestream `alpha`.
"""
struct Stability <: AbstractFrame end

"""
    Wind <: AbstractFrame

Reference frame rotated to be aligned with the freestream `alpha` and `beta`
"""
struct Wind <: AbstractFrame end

# --- struct to store panel properties --- #

"""
    PanelProperties

Panel specific properties calculated during the vortex lattice method analysis.

**Fields**
 - `gamma`: Panel circulation strength (normalized by the freestream velocity)
 - `v`: Local velocity at the panel's center (typically the quarter-chord), normalized
 - `cf`: Bound vortex force per unit length, normalized by `QINF*S` where `QINF`
    is the dynamic pressure and `S` is the user-provided reference area.
 - `cfl`: Left vortex force per unit length, normalized by `QINF*S`
 - `cfr`: Right vortex force per unit length, normalized by `QINF*S`
"""
struct PanelProperties{TF}
    gamma::TF
    v::SVector{3, TF}
    cf::SVector{3, TF}
    cfl::SVector{3, TF}
    cfr::SVector{3, TF}
end

function PanelProperties(gamma, v, cf)

    TF = promote_type(typeof(gamma), typeof(v), eltype(cf))

    return PanelProperties{TF}(cf, v, gamma)
end

Base.eltype(::Type{PanelProperties{TF}}) where TF = TF
Base.eltype(::PanelProperties{TF}) where TF = TF

# --- near field solution for forces and moments --- #

"""
    near_field_forces(panels, reference, freestream, symmetric, Γ; xhat=[1,0,0], frame=Body())

Compute the forces and moments acting on the aircraft given the circulation
distribution `Γ`.  Return `CF`, `CM`, and a vector of panel properties of type
`PanelProperties`.  `CF` and `CM` are returned in the frame specified by `frame`.
"""
@inline function near_field_forces(panels, ref, fs, symmetric, Γ, xhat=SVector(1, 0, 0); frame=Body())

    # float number type
    TF = promote_type(eltype(eltype(panels)), eltype(ref), eltype(fs), eltype(Γ), eltype(xhat))

    # number of panels
    N = length(panels)

    # initialize forces/moments in the body frame
    Fb = @SVector zeros(TF, 3)  # forces
    Mb = @SVector zeros(TF, 3)  # moments

    # initialize storage for panel outputs
    pprops = Vector{PanelProperties{TF}}(undef, N)

    # loop through receiving panels
    for i = 1:N

        # --- Calculate forces for the bound vortices --- #

        rc = midpoint(panels[i])

        # get external velocity at the bound vortex midpoint
        Vi = external_velocity(fs, rc, ref.r)

        # now add the velocity induced by each sending panel
        for j = 1:N

            # note that for the horseshoe ring case we assume panels are
            # grouped into chordwise strips that are ordered from leading edge
            # to trailing edge
            # e.g.
            #     1 4 7 10         1 10 7 4
            #     2 5 8 11   and   2 11 8 5   are both valid orderings
            #     3 6 9 12         3 12 9 6

            # get induced velocity for this sending panel normalized by Γ
            if eltype(panels) <: Horseshoe
                include_bound = i != j
                Vij = induced_velocity(rc, panels[j], symmetric, xhat, i != j)
            elseif eltype(panels) <: Ring
                include_top = i != j
                include_bottom = i != j+1
                Vij = induced_velocity(rc, panels[j], symmetric, xhat, include_top, include_bottom)
            end

            # add contribution to velocity
            Vi += Vij*Γ[j]
        end

        # now use the Kutta-Joukowski theorem to calculate the forces at this point
        if eltype(panels) <: Horseshoe
            Γi = Γ[i]
        elseif eltype(panels) <: Ring
            if i == 1 || panels[i-1].trailing
                # panel is on the leading edge
                Γi = Γ[i]
            else
                # panel is not on the leading edge
                Γi = Γ[i] - Γ[i-1]
            end
        end

        Δs = top_vector(panels[i])
        lenb = norm(Δs)
        Fbi = RHO*Γi*cross(Vi, Δs)

        # add the resulting forces and moments to the body forces and moments
        Δr = rc - ref.r
        Fb += Fbi
        Mb += cross(Δr, Fbi)

        # --- Calculate forces for the left bound vortex --- #

        rc = left_midpoint(panels[i])

        # get effective velocity at the bound vortex midpoint
        Veff = external_velocity(fs, rc, ref.r)

        # note that we don't include induced velocity in the effective velocity
        # because its influence is generally negligible once we take the cross product
        # with the bound vortex vector, this is also assumed in AVL

        # now use the Kutta-Joukowski theorem to calculate the forces
        Γli = Γ[i]
        Δs = left_vector(panels[i])
        lenl = norm(Δs)
        Fbli = RHO*Γli*cross(Veff, Δs)

        # add the resulting forces and moments to the body forces and moments
        Δr = rc - ref.r
        Fb += Fbli
        Mb += cross(Δr, Fbli)

        # --- Calculate forces for the right bound vortex --- #

        rc = right_midpoint(panels[i])

        # get effective velocity at the bound vortex midpoint
        Veff = external_velocity(fs, rc, ref.r)

        # note that we don't include induced velocity in the effective velocity
        # because its influence is generally negligible once we take the cross product
        # with the bound vortex vector, this is also assumed in AVL

        # now use the Kutta-Joukowski theorem to calculate the forces
        Γri = Γ[i]
        Δs = right_vector(panels[i])
        lenr = norm(Δs)
        Fbri = RHO*Γri*cross(Veff, Δs)

        # add the resulting forces and moments to the body forces and moments
        Δr = rc - ref.r
        Fb += Fbri
        Mb += cross(Δr, Fbri)

        # assemble normalized panel outputs
        pprops[i] = PanelProperties(Γi/VINF, Vi/VINF, Fbi/(QINF*ref.S*lenb),
            Fbli/(QINF*ref.S*lenl), Fbri/(QINF*ref.S*lenr))
    end

    # adjust body forces to account for symmetry
    if symmetric
        Fb = SVector(2*Fb[1], 0.0, 2*Fb[3])
        Mb = SVector(0.0, 2*Mb[2], 0.0)
    end

    # normalize body forces
    reference_length = SVector(ref.b, ref.c, ref.b)
    CF = Fb./(QINF.*ref.S)
    CM = Mb./(QINF.*ref.S.*reference_length)

    CF, CM = body_to_frame(CF, CM, ref, fs, frame)

    return CF, CM, pprops
end

"""
    near_field_forces_derivatives(panels, reference, freestream, symmetric, Γ, dΓ; xhat=[1,0,0])

Compute the forces and moments acting on the aircraft given the circulation
distribution `Γ` and their derivatives with respect to the variables in `freestream`.
Return `CF`, `CM`, `dCF`, `dCM`, and a vector of panel properties of type
`PanelProperties`.  `CF`, `CM`, `dCF`, and `dCM` are returned in the body frame.
"""
@inline function near_field_forces_derivatives(panels, ref, fs, symmetric, Γ, dΓ;
        xhat = SVector(1, 0, 0), frame=Body())

    # unpack derivatives
    Γ_a, Γ_b, Γ_p, Γ_q, Γ_r = dΓ

    # float number type
    TF = promote_type(eltype(eltype(panels)), eltype(ref), eltype(fs),
        eltype(Γ), eltype(Γ_a), eltype(Γ_b), eltype(Γ_p), eltype(Γ_q), eltype(Γ_r))

    # number of panels
    N = length(panels)

    # initialize body outputs
    Fb = @SVector zeros(TF, 3)  # forces
    Mb = @SVector zeros(TF, 3)  # moments

    # initialize derivatives of body outputs
    Fb_a = @SVector zeros(TF, 3)
    Fb_b = @SVector zeros(TF, 3)
    Fb_p = @SVector zeros(TF, 3)
    Fb_q = @SVector zeros(TF, 3)
    Fb_r = @SVector zeros(TF, 3)

    Mb_a = @SVector zeros(TF, 3)
    Mb_b = @SVector zeros(TF, 3)
    Mb_p = @SVector zeros(TF, 3)
    Mb_q = @SVector zeros(TF, 3)
    Mb_r = @SVector zeros(TF, 3)

    # initialize panel outputs
    pprops = Vector{PanelProperties{TF}}(undef, N)

    # loop through receiving panels
    for i = 1:N

        # --- Calculate forces for the bound vortices --- #

        rc = midpoint(panels[i])

        # get external velocity at the bound vortex midpoint
        Vi, dVi = external_velocity_derivatives(fs, rc, ref.r)

        # unpack derivatives
        Vi_a, Vi_b, Vi_p, Vi_q, Vi_r = dVi

        # now add the velocity induced by each sending panel
        for j = 1:N

            # note that for the horseshoe ring case we assume panels are
            # grouped into chordwise strips that are ordered from leading edge
            # to trailing edge
            # e.g.
            #     1 4 7 10         1 10 7 4
            #     2 5 8 11   and   2 11 8 5   are both valid orderings
            #     3 6 9 12         3 12 9 6

            # get induced velocity for this sending panel normalized by Γ
            if eltype(panels) <: Horseshoe
                include_bound = i != j
                Vij = induced_velocity(rc, panels[j], symmetric, xhat, include_bound)
            elseif eltype(panels) <: Ring
                include_top = i != j
                include_bottom = i != j+1
                Vij = induced_velocity(rc, panels[j], symmetric, xhat, include_top, include_bottom)
            end

            # add contribution to induced velocity
            Vi += Vij*Γ[j]

            # derivatives
            Vi_a += Vij*Γ_a[j]
            Vi_b += Vij*Γ_b[j]
            Vi_p += Vij*Γ_p[j]
            Vi_q += Vij*Γ_q[j]
            Vi_r += Vij*Γ_r[j]
        end

        # now use the Kutta-Joukowski theorem to calculate the forces at this point
        if eltype(panels) <: Horseshoe
            Γi = Γ[i]

            Γi_a = Γ_a[i]
            Γi_b = Γ_b[i]
            Γi_p = Γ_p[i]
            Γi_q = Γ_q[i]
            Γi_r = Γ_r[i]
        elseif eltype(panels) <: Ring
            if i == 1 || panels[i-1].trailing
                # panel is on the leading edge
                Γi = Γ[i]

                Γi_a = Γ_a[i]
                Γi_b = Γ_b[i]
                Γi_p = Γ_p[i]
                Γi_q = Γ_q[i]
                Γi_r = Γ_r[i]
            else
                Γi = Γ[i] - Γ[i-1]

                Γi_a = Γ_a[i] - Γ_a[i-1]
                Γi_b = Γ_b[i] - Γ_b[i-1]
                Γi_p = Γ_p[i] - Γ_p[i-1]
                Γi_q = Γ_q[i] - Γ_q[i-1]
                Γi_r = Γ_r[i] - Γ_r[i-1]
            end
        end

        Δs = top_vector(panels[i])
        lenb = norm(Δs)
        tmp = cross(Vi, Δs)
        Fbi = RHO*Γi*tmp

        Fbi_a = RHO*(Γi_a*tmp + Γi*cross(Vi_a, Δs))
        Fbi_b = RHO*(Γi_b*tmp + Γi*cross(Vi_b, Δs))
        Fbi_p = RHO*(Γi_p*tmp + Γi*cross(Vi_p, Δs))
        Fbi_q = RHO*(Γi_q*tmp + Γi*cross(Vi_q, Δs))
        Fbi_r = RHO*(Γi_r*tmp + Γi*cross(Vi_r, Δs))

        # add the resulting forces and moments to the body forces and moments
        Δr = rc - ref.r
        Fb += Fbi
        Mb += cross(Δr, Fbi)

        # also add their derivatives
        Fb_a += Fbi_a
        Fb_b += Fbi_b
        Fb_p += Fbi_p
        Fb_q += Fbi_q
        Fb_r += Fbi_r

        Mb_a += cross(Δr, Fbi_a)
        Mb_b += cross(Δr, Fbi_b)
        Mb_p += cross(Δr, Fbi_p)
        Mb_q += cross(Δr, Fbi_q)
        Mb_r += cross(Δr, Fbi_r)

        # --- Calculate forces for the left bound vortex --- #

        rc = left_midpoint(panels[i])

        # get effective velocity at the bound vortex midpoint
        Veff, dVeff = external_velocity_derivatives(fs, rc, ref.r)

        # unpack derivatives
        Veff_a, Veff_b, Veff_p, Veff_q, Veff_r = dVeff

        # note that we don't include induced velocity in the effective velocity
        # because its influence is generally negligible once we take the cross product
        # with the bound vortex vector, this is also assumed in AVL

        # now use the Kutta-Joukowski theorem to calculate the forces
        Γli = Γ[i]
        Δs = left_vector(panels[i])
        lenl = norm(Δs)
        tmp = cross(Veff, Δs)
        Fbli = RHO*Γli*tmp

        Γli_a = Γ_a[i]
        Γli_b = Γ_b[i]
        Γli_p = Γ_p[i]
        Γli_q = Γ_q[i]
        Γli_r = Γ_r[i]

        Fbli_a = RHO*(Γli_a*tmp + Γli*cross(Veff_a, Δs))
        Fbli_b = RHO*(Γli_b*tmp + Γli*cross(Veff_b, Δs))
        Fbli_p = RHO*(Γli_p*tmp + Γli*cross(Veff_p, Δs))
        Fbli_q = RHO*(Γli_q*tmp + Γli*cross(Veff_q, Δs))
        Fbli_r = RHO*(Γli_r*tmp + Γli*cross(Veff_r, Δs))

        # add the resulting forces and moments to the body forces and moments
        Δr = rc - ref.r
        Fb += Fbli
        Mb += cross(Δr, Fbli)

        # also add their derivatives
        Fb_a += Fbli_a
        Fb_b += Fbli_b
        Fb_p += Fbli_p
        Fb_q += Fbli_q
        Fb_r += Fbli_r

        Mb_a += cross(Δr, Fbli_a)
        Mb_b += cross(Δr, Fbli_b)
        Mb_p += cross(Δr, Fbli_p)
        Mb_q += cross(Δr, Fbli_q)
        Mb_r += cross(Δr, Fbli_r)

        # --- Calculate forces for the right bound vortex --- #

        rc = right_midpoint(panels[i])

        # get effective velocity at the bound vortex midpoint
        Veff, dVeff = external_velocity_derivatives(fs, rc, ref.r)

        # unpack derivatives
        Veff_a, Veff_b, Veff_p, Veff_q, Veff_r = dVeff

        # note that we don't include induced velocity in the effective velocity
        # because its influence is generally negligible once we take the cross product
        # with the bound vortex vector, this is also assumed in AVL

        # now use the Kutta-Joukowski theorem to calculate the forces
        Γri = Γ[i]
        Δs = right_vector(panels[i])
        lenr = norm(Δs)
        tmp = cross(Veff, Δs)
        Fbri = RHO*Γri*tmp

        Γri_a = Γ_a[i]
        Γri_b = Γ_b[i]
        Γri_p = Γ_p[i]
        Γri_q = Γ_q[i]
        Γri_r = Γ_r[i]

        Fbri_a = RHO*(Γri_a*tmp + Γri*cross(Veff_a, Δs))
        Fbri_b = RHO*(Γri_b*tmp + Γri*cross(Veff_b, Δs))
        Fbri_p = RHO*(Γri_p*tmp + Γri*cross(Veff_p, Δs))
        Fbri_q = RHO*(Γri_q*tmp + Γri*cross(Veff_q, Δs))
        Fbri_r = RHO*(Γri_r*tmp + Γri*cross(Veff_r, Δs))

        # add the resulting forces and moments to the body forces and moments
        Δr = rc - ref.r
        Fb += Fbri
        Mb += cross(Δr, Fbri)

        # also add their derivatives
        Fb_a += Fbri_a
        Fb_b += Fbri_b
        Fb_p += Fbri_p
        Fb_q += Fbri_q
        Fb_r += Fbri_r

        Mb_a += cross(Δr, Fbri_a)
        Mb_b += cross(Δr, Fbri_b)
        Mb_p += cross(Δr, Fbri_p)
        Mb_q += cross(Δr, Fbri_q)
        Mb_r += cross(Δr, Fbri_r)

        # assemble normalized panel outputs
        pprops[i] = PanelProperties(Γi/VINF, Vi/VINF, Fbi/(QINF*ref.S*lenb),
            Fbli/(QINF*ref.S*lenl), Fbri/(QINF*ref.S*lenr))
    end

    # adjust body forces to account for symmetry
    if symmetric
        Fb = SVector(2*Fb[1], 0.0, 2*Fb[3])
        Mb = SVector(0.0, 2*Mb[2], 0.0)

        # derivatives
        Fb_a = SVector(2*Fb_a[1], 0.0, 2*Fb_a[3])
        Fb_b = SVector(2*Fb_b[1], 0.0, 2*Fb_b[3])
        Fb_p = SVector(2*Fb_p[1], 0.0, 2*Fb_p[3])
        Fb_q = SVector(2*Fb_q[1], 0.0, 2*Fb_q[3])
        Fb_r = SVector(2*Fb_r[1], 0.0, 2*Fb_r[3])

        Mb_a = SVector(0.0, 2*Mb_a[2], 0.0)
        Mb_b = SVector(0.0, 2*Mb_b[2], 0.0)
        Mb_p = SVector(0.0, 2*Mb_p[2], 0.0)
        Mb_q = SVector(0.0, 2*Mb_q[2], 0.0)
        Mb_r = SVector(0.0, 2*Mb_r[2], 0.0)
    end

    # normalize body forces
    reference_length = SVector(ref.b, ref.c, ref.b)
    CF = Fb/(QINF*ref.S)
    CM = Mb./(QINF.*ref.S.*reference_length)

    # also normalize derivatives
    CF_a = Fb_a/(QINF*ref.S)
    CF_b = Fb_b/(QINF*ref.S)
    CF_p = Fb_p/(QINF*ref.S)
    CF_q = Fb_q/(QINF*ref.S)
    CF_r = Fb_r/(QINF*ref.S)

    CM_a = Mb_a./(QINF.*ref.S.*reference_length)
    CM_b = Mb_b./(QINF.*ref.S.*reference_length)
    CM_p = Mb_p./(QINF.*ref.S.*reference_length)
    CM_q = Mb_q./(QINF.*ref.S.*reference_length)
    CM_r = Mb_r./(QINF.*ref.S.*reference_length)

    dCF = (CF_a, CF_b, CF_p, CF_q, CF_r)
    dCM = (CM_a, CM_b, CM_p, CM_q, CM_r)

    return CF, CM, dCF, dCM, pprops
end

"""
    body_to_frame(CF, CM, ref, fs, frame)

Transforms the coefficients `CF` and `CM` which are provided in the body frame
to another frame
"""
body_to_frame

@inline body_to_frame(CF, CM, ref, fs, ::Body) = CF, CM

@inline function body_to_frame(CF, CM, ref, fs, ::Stability)
    R = body_to_stability(fs)
    return R*CF, R*CM
end

@inline function body_to_frame(CF, CM, ref, fs, ::Wind)
    # remove reference lengths
    reflen = SVector(ref.b, ref.c, ref.b)
    CM = CM .* reflen

    # rotate
    R = body_to_wind(fs)
    CF = R*CF
    CM = R*CM

    # reapply reference lengths
    CM = CM ./ reflen

    return CF, CM
end

"""
    body_to_frame_derivatives(CF, CM, dCF, dCM, ref, fs, frame)

Transforms the coefficients `CF` and `CM` (and their derivatives with respect to
the freestream variables) from the body frame into the specified `frame`
"""
body_to_frame_derivatives

@inline body_to_frame_derivatives(CF, CM, dCF, dCM, ref, fs, ::Body) = CF, CM, dCF, dCM

@inline function body_to_frame_derivatives(CF, CM, dCF, dCM, ref, fs, ::Stability)

    (CF_a, CF_b, CF_p, CF_q, CF_r) = dCF
    (CM_a, CM_b, CM_p, CM_q, CM_r) = dCM

    R, R_a = body_to_stability_alpha(fs)

    R*

    return R*CF, R*CM
end

@inline function body_to_frame_derivatives(CF, CM, ref, fs, ::Wind)
    # remove reference lengths
    reflen = SVector(ref.b, ref.c, ref.b)
    CM = CM .* reflen

    # rotate
    R = body_to_wind(fs)
    CF = R*CF
    CM = R*CM

    # reapply reference lengths
    CM = CM ./ reflen

    return CF, CM
end
