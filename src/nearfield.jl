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

function PanelProperties(gamma, v, cf, cfl, cfr)

    TF = promote_type(typeof(gamma), typeof(v), eltype(cf), eltype(cfl), eltype(cfr))

    return PanelProperties{TF}(gamma, v, cf, cfl, cfr)
end

Base.eltype(::Type{PanelProperties{TF}}) where TF = TF
Base.eltype(::PanelProperties{TF}) where TF = TF

# --- near field solution for forces and moments --- #

"""
    near_field_forces(surface[s], reference, freestream, symmetric, Γ; kwargs...)

Compute the forces and moments acting on the aircraft given the circulation
distribution `Γ`.

Return `CF`, `CM`, and a vector of panel properties of type `PanelProperties`.

# Keyword Arguments
 - `xhat`: direction in which trailing vortices are shed, defaults to [1,0,0]
 - `frame`: frame in which to return `CF` and `CM`, options are `Body()` (default),
    `Stability()`, and `Wind()`

# Additional Keyword Arguments for Multiple Surfaces
 - `surface_id`: ID for each surface, defaults to `1:length(surfaces)`
"""
near_field_forces

# single surface
function near_field_forces(surface::AbstractMatrix, ref, fs, symmetric, Γ;
    xhat=SVector(1, 0, 0), frame=Body())

    CF, CM, props = near_field_forces([surface], ref, fs, symmetric, Γ;
        xhat=xhat, frame=frame, surface_id=1:1)

    return CF, CM, props[1]
end

# multiple surfaces
function near_field_forces(surfaces::AbstractVector{<:AbstractMatrix}, ref, fs,
    symmetric, Γ; xhat=SVector(1, 0, 0), frame=Body(), surface_id=1:length(surfaces))

    # float number type
    TF = promote_type(eltype(eltype(eltype(surfaces))), eltype(ref), eltype(fs), eltype(Γ), eltype(xhat))

    # number of surfaces
    nsurf = length(surfaces)

    # initialize forces/moments in the body frame
    Fb = @SVector zeros(TF, 3)  # forces
    Mb = @SVector zeros(TF, 3)  # moments

    # initialize storage for panel properties
    props = Vector{Matrix{PanelProperties{TF}}}(undef, nsurf)

    # loop through receiving surfaces
    iΓ = 0 # index for accessing Γ
    for isurf = 1:nsurf

        receiving = surfaces[isurf]
        nr = length(receiving)
        nr1, nr2 = size(receiving)
        cr = CartesianIndices(receiving)

        # initialize panel outputs for this surface
        props[isurf] = Matrix{PanelProperties{TF}}(undef, nr1, nr2)

        # loop through receiving panels
        for i = 1:length(receiving)

            # --- Calculate forces for the bound vortices --- #
            I = cr[i]

            rc = midpoint(receiving[I])

            # start with external velocity at the bound vortex midpoint
            Vi = external_velocity(fs, rc, ref.r)

            # now add the velocity induced by each sending panel
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                sending = surfaces[jsurf]
                ns = length(sending)
                ns1, ns2 = size(sending)
                cs = CartesianIndices(sending)

                # find out whether surfaces have the same index and/or ID
                same_surface = isurf == jsurf
                same_id = surface_id[isurf] == surface_id[jsurf]

                # loop through each sending panel
                for j = 1:ns

                    J = cs[j]

                    if typeof(sending[J]) <: Horseshoe
                        include_bound = !(same_surface && I == J)
                        Vij = induced_velocity(rc, sending[J], xhat, symmetric, same_id, include_bound)
                    elseif typeof(sending[J]) <: Ring
                        include_top = !(same_surface && I == J)
                        include_bottom = !(same_surface && (I[1] == J[1]+1 && I[2] == J[2]))
                        trailing = J[1] == ns1
                        Vij = induced_velocity(rc, sending[J], xhat, symmetric,
                            same_id, trailing, include_top, include_bottom)
                    end

                    # add contribution to velocity
                    Vi += Vij*Γ[jΓ+j]
                end
                jΓ += ns
            end

            # now use the Kutta-Joukowski theorem to calculate the forces at this point
            if typeof(receiving[I]) <: Horseshoe
                Γi = Γ[iΓ+i]
            elseif typeof(receiving[I]) <: Ring
                if I[1] == 1
                    # panel is on the leading edge
                    Γi = Γ[iΓ+i]
                else
                    # panel is not on the leading edge
                    Γi = Γ[iΓ+i] - Γ[iΓ+i-1]
                end
            end

            Δs = top_vector(receiving[I])
            lenb = norm(Δs)
            Fbi = RHO*Γi*cross(Vi, Δs)

            # add the resulting forces and moments to the body forces and moments
            Δr = rc - ref.r
            Fb += Fbi
            Mb += cross(Δr, Fbi)

            # --- Calculate forces for the left bound vortex --- #

            rc = left_midpoint(receiving[I])

            # get effective velocity at the bound vortex midpoint
            Veff = external_velocity(fs, rc, ref.r)

            # note that we don't include induced velocity in the effective velocity
            # because its influence is generally negligible once we take the cross product
            # with the bound vortex vector, this is also assumed in AVL

            # now use the Kutta-Joukowski theorem to calculate the forces
            Γli = Γ[iΓ+i]
            Δs = left_vector(receiving[I])
            lenl = norm(Δs)
            Fbli = RHO*Γli*cross(Veff, Δs)

            # add the resulting forces and moments to the body forces and moments
            Δr = rc - ref.r
            Fb += Fbli
            Mb += cross(Δr, Fbli)

            # --- Calculate forces for the right bound vortex --- #

            rc = right_midpoint(receiving[I])

            # get effective velocity at the bound vortex midpoint
            Veff = external_velocity(fs, rc, ref.r)

            # note that we don't include induced velocity in the effective velocity
            # because its influence is generally negligible once we take the cross product
            # with the bound vortex vector, this is also assumed in AVL

            # now use the Kutta-Joukowski theorem to calculate the forces
            Γri = Γ[iΓ+i]
            Δs = right_vector(receiving[I])
            lenr = norm(Δs)
            Fbri = RHO*Γri*cross(Veff, Δs)

            # add the resulting forces and moments to the body forces and moments
            Δr = rc - ref.r
            Fb += Fbri
            Mb += cross(Δr, Fbri)

            # assemble normalized panel outputs
            props[isurf][I] = PanelProperties(Γi/VINF, Vi/VINF, Fbi/(QINF*ref.S*lenb),
                Fbli/(QINF*ref.S*lenl), Fbri/(QINF*ref.S*lenr))
        end
        # increment Γ index
        iΓ += nr
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

    return CF, CM, props
end

"""
    near_field_forces_derivatives(surface[s], reference, freestream, symmetric, Γ, dΓ; xhat=[1,0,0])

Compute the forces and moments acting on the aircraft given the circulation
distribution `Γ` and their derivatives with respect to the variables in `freestream`.
Return `CF`, `CM`, `dCF`, `dCM`, and a vector of panel properties of type
`PanelProperties`.  `CF`, `CM`, `dCF`, and `dCM` are returned in the body frame.
"""
near_field_forces_derivatives

# single surface
function near_field_forces_derivatives(surface::AbstractMatrix, ref, fs, symmetric, Γ, dΓ;
    xhat=SVector(1, 0, 0), frame=Body())

    CF, CM, dCF, dCM, props = near_field_forces_derivatives([surface], ref, fs, symmetric, Γ, dΓ;
        xhat=xhat, frame=frame, surface_id=1:1)

    return CF, CM, dCF, dCM, props[1]
end

# multiple surfaces
function near_field_forces_derivatives(surfaces::AbstractVector{<:AbstractMatrix}, ref, fs, symmetric, Γ, dΓ;
        xhat = SVector(1, 0, 0), frame=Body(), surface_id=1:length(surfaces))

    # unpack derivatives
    Γ_a, Γ_b, Γ_p, Γ_q, Γ_r = dΓ

    # float number type
    TF = promote_type(eltype(eltype(eltype(surfaces))), eltype(ref), eltype(fs),
        eltype(Γ), eltype(Γ_a), eltype(Γ_b), eltype(Γ_p), eltype(Γ_q), eltype(Γ_r))

    # number of surfaces
    nsurf = length(surfaces)

    # initialize forces/moments in the body frame
    Fb = @SVector zeros(TF, 3)  # forces
    Mb = @SVector zeros(TF, 3)  # moments

    # and their derivatives
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
    props = Vector{Matrix{PanelProperties{TF}}}(undef, nsurf)

    # loop through receiving surfaces
    iΓ = 0 # index for accessing Γ
    for isurf = 1:nsurf

        receiving = surfaces[isurf]
        nr = length(receiving)
        nr1, nr2 = size(receiving)
        cr = CartesianIndices(receiving)

        # initialize panel outputs for this surface
        props[isurf] = Matrix{PanelProperties{TF}}(undef, nr1, nr2)

        # loop through receiving panels
        for i = 1:length(receiving)

            I = cr[i]

            # --- Calculate forces for the bound vortices --- #
            rc = midpoint(receiving[I])

            # get external velocity at the bound vortex midpoint
            Vi, dVi = external_velocity_derivatives(fs, rc, ref.r)

            # unpack derivatives
            Vi_a, Vi_b, Vi_p, Vi_q, Vi_r = dVi

            # now add the velocity induced by each sending panel
            jΓ = 0 # index for accessing Γ
            for jsurf = 1:nsurf

                sending = surfaces[jsurf]
                ns = length(sending)
                ns1, ns2 = size(sending)
                cs = CartesianIndices(sending)

                # find out whether surfaces have the same index and/or ID
                same_surface = isurf == jsurf
                same_id = surface_id[isurf] == surface_id[jsurf]

                # loop through each sending panel
                for j = 1:ns

                    J = cs[j]

                    if typeof(sending[J]) <: Horseshoe
                        include_bound = !(same_surface && I == J)
                        Vij = induced_velocity(rc, sending[J], xhat, symmetric, same_id, include_bound)
                    elseif typeof(sending[J]) <: Ring
                        include_top = !(same_surface && I == J)
                        include_bottom = !(same_surface && (I[1] == J[1]+1 && I[2] == J[2]))
                        trailing = J[1] == ns1
                        Vij = induced_velocity(rc, sending[J], xhat, symmetric,
                            same_id, trailing, include_top, include_bottom)
                    end

                    # add contribution to velocity
                    Vi += Vij*Γ[jΓ+j]

                    # derivatives
                    Vi_a += Vij*Γ_a[j]
                    Vi_b += Vij*Γ_b[j]
                    Vi_p += Vij*Γ_p[j]
                    Vi_q += Vij*Γ_q[j]
                    Vi_r += Vij*Γ_r[j]
                end
                # increment Γ index
                jΓ += ns
            end

            # now use the Kutta-Joukowski theorem to calculate the forces at this point
            if typeof(receiving[I]) <: Horseshoe
                Γi = Γ[iΓ+i]

                Γi_a = Γ_a[iΓ+i]
                Γi_b = Γ_b[iΓ+i]
                Γi_p = Γ_p[iΓ+i]
                Γi_q = Γ_q[iΓ+i]
                Γi_r = Γ_r[iΓ+i]
            elseif typeof(receiving[I]) <: Ring
                if I[1] == 1
                    # panel is on the leading edge
                    Γi = Γ[iΓ+i]

                    Γi_a = Γ_a[iΓ+i]
                    Γi_b = Γ_b[iΓ+i]
                    Γi_p = Γ_p[iΓ+i]
                    Γi_q = Γ_q[iΓ+i]
                    Γi_r = Γ_r[iΓ+i]
                else
                    # panel is not on the leading edge
                    Γi = Γ[iΓ+i] - Γ[iΓ+i-1]

                    Γi_a = Γ_a[iΓ+i] - Γ_a[iΓ+i-1]
                    Γi_b = Γ_b[iΓ+i] - Γ_b[iΓ+i-1]
                    Γi_p = Γ_p[iΓ+i] - Γ_p[iΓ+i-1]
                    Γi_q = Γ_q[iΓ+i] - Γ_q[iΓ+i-1]
                    Γi_r = Γ_r[iΓ+i] - Γ_r[iΓ+i-1]
                end
            end

            Δs = top_vector(receiving[I])
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

            rc = left_midpoint(receiving[I])

            # get effective velocity at the bound vortex midpoint
            Veff, dVeff = external_velocity_derivatives(fs, rc, ref.r)

            # unpack derivatives
            Veff_a, Veff_b, Veff_p, Veff_q, Veff_r = dVeff

            # note that we don't include induced velocity in the effective velocity
            # because its influence is generally negligible once we take the cross product
            # with the bound vortex vector, this is also assumed in AVL

            # now use the Kutta-Joukowski theorem to calculate the forces
            Γli = Γ[iΓ+i]
            Δs = left_vector(receiving[I])
            lenl = norm(Δs)
            tmp = cross(Veff, Δs)
            Fbli = RHO*Γli*tmp

            Γli_a = Γ_a[iΓ+i]
            Γli_b = Γ_b[iΓ+i]
            Γli_p = Γ_p[iΓ+i]
            Γli_q = Γ_q[iΓ+i]
            Γli_r = Γ_r[iΓ+i]

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

            rc = right_midpoint(receiving[I])

            # get effective velocity at the bound vortex midpoint
            Veff, dVeff = external_velocity_derivatives(fs, rc, ref.r)

            # unpack derivatives
            Veff_a, Veff_b, Veff_p, Veff_q, Veff_r = dVeff

            # note that we don't include induced velocity in the effective velocity
            # because its influence is generally negligible once we take the cross product
            # with the bound vortex vector, this is also assumed in AVL

            # now use the Kutta-Joukowski theorem to calculate the forces
            Γri = Γ[iΓ+i]
            Δs = right_vector(receiving[I])
            lenr = norm(Δs)
            tmp = cross(Veff, Δs)
            Fbri = RHO*Γri*tmp

            Γri_a = Γ_a[iΓ+i]
            Γri_b = Γ_b[iΓ+i]
            Γri_p = Γ_p[iΓ+i]
            Γri_q = Γ_q[iΓ+i]
            Γri_r = Γ_r[iΓ+i]

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
            props[isurf][I] = PanelProperties(Γi/VINF, Vi/VINF, Fbi/(QINF*ref.S*lenb),
                Fbli/(QINF*ref.S*lenl), Fbri/(QINF*ref.S*lenr))
        end
        # increment Γ index
        iΓ += nr
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

    return CF, CM, dCF, dCM, props
end

"""
    body_to_frame(CF, CM, ref, fs, frame)

Transforms the coefficients `CF` and `CM` which are provided in the body frame
to another frame
"""
body_to_frame

body_to_frame(CF, CM, ref, fs, ::Body) = CF, CM

function body_to_frame(CF, CM, ref, fs, ::Stability)
    R = body_to_stability(fs)
    return R*CF, R*CM
end

function body_to_frame(CF, CM, ref, fs, ::Wind)
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
