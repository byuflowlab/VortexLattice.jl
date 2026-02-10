struct Polar{TF}
    alphas::Vector{TF}
    cls_inv::Vector{TF}
    cls_delta::Vector{TF}
    cls_visc::Vector{TF}
    cds_visc::Vector{TF}
    cl_alpha0::TF
    m_inv::TF
end

function Polar(alphas, cls_inv, cls_visc, cds_visc)
    cls_delta = cls_visc .- cls_inv
    m_inv = (cls_inv[end] - cls_inv[1]) / (alphas[end] - alphas[1])
    cl_alpha0 = FLOWMath.linear(alphas, cls_inv, 0.0)
    return Polar{eltype(cls_inv)}(alphas, cls_inv, cls_delta, cls_visc, cds_visc, cl_alpha0, m_inv)
end

"""
    viscous!(properties, Γ, dΓdt, surfaces, grids, frames, frames_index, viscous_ratio_cl, viscous_ratio_cd)

Apply viscous corrections to the aerodynamic forces and circulation strengths based on the provided viscous correction functions `viscous_ratio_cl` and `viscous_ratio_cd`.
The corrections are applied to the `properties` of each panel, as well as the circulation strengths `Γ` and their time derivatives `dΓdt`.

Note: dΓdt should contain -Γ from the PREVIOUS timestep
"""
function viscous!(properties::Vector{Matrix{PanelProperties{TF}}}, Γ, dΓdt, surfaces::Vector{Matrix{SurfacePanel{TF}}}, grids, frames::Vector{<:ReferenceFrame}, frames_index::Vector{Int}, polar::Polar, ref::Reference, dt) where TF
    # properties contains:
    # * cf
    # surface contains:
    # * position
    # grid contains:
    # * chord information
    # Γ contains:
    # * circulation
    # frame contains:
    # * rotation matrix into this frame

    # index for circulation strengths
    iΓ = 1

    # loop over surfaces
    for isurf in eachindex(surfaces)

        if frames_index[isurf] > 0

            # extract containers corresponding to this surface
            surface = surfaces[isurf]
            props = properties[isurf]
            grid = grids[isurf]

            # get reference frame for this surface
            frame = frames[frames_index[isurf]]

            # get rotation matrix from global frame to this frame
            Rp = frame.Rp2g * frame.R # rotation from this frame to global frame
            R = transpose(Rp)

            # @show R # verified: basis vectors of the frame expressed in global coordinates

            # loop over spanwise sections in this surface
            for j in axes(surface, 2)

                # get the circulation per unit length of this section
                γ = zero(TF)

                # get induced velocity
                v_induced = zero(SVector{3,TF})

                # get net aerodynamic force on this section
                cf = zero(SVector{3,TF})

                # loop over chordwise panels in this section
                for i in axes(surface, 1)

                    # project bound vortex onto y axis of this frame
                    ds = R * (surface[i,j].rtr - surface[i,j].rtl)
                    dy = abs(ds[2]) / norm(ds)

                    # accumulate circulation contribution from this bound vortex
                    γ += Γ[iΓ] * dy
                    iΓ += 1

                    # accumulate induced velocity contribution at this bound vortex
                    v_induced += props[i,j].velocity # TODO: what if system.reference[].v != 1.0?

                    # accumulate aerodynamic force contribution from this bound vortex
                    cf += props[i,j].cfb
                end

                # get chord length
                le = SVector(grid[1,1,j], grid[2,1,j], grid[3,1,j])
                te = SVector(grid[1,end,j], grid[2,end,j], grid[3,end,j])
                c = norm(le - te)

                # average of chord lengths of either side
                le = SVector(grid[1,1,j+1], grid[2,1,j+1], grid[3,1,j+1])
                te = SVector(grid[1,end,j+1], grid[2,end,j+1], grid[3,end,j+1])
                c += norm(le - te)
                c *= 0.5

                # project aerodynamic force into xz plane
                cf = R * cf # rotate into this frame

                # get 2-D lift magnitude of this section
                l_2d_norm = sqrt(cf[1]*cf[1] + cf[3]*cf[3])

                # force per length
                l_2d_norm /= norm(surface[end,j].rtr - surface[1,j].rtl)

                # calculate effective cl predicted by the VLM, = 2π * α_eff
                cl_vlm = -2 * RHO * γ * γ / (l_2d_norm * c) * sign(γ)

                # get effective α
                α_eff = cl_vlm / (2 * pi) * pi/180

                # refer to polar for viscous cl
                cl_star = FLOWMath.linear(polar.alphas, polar.cls_visc, α_eff)

                # # correct for alpha=0 cl, inviscid lift slope, and viscous correction
                # cl_star = polar.m_inv / (2*pi) * cl_vlm + polar.cl_alpha0
                # cl_star = cl_star + FLOWMath.linear(polar.cls_inv, polar.cls_delta, cl_star)

                # get viscous lift correction factor
                f_cl = cl_star / cl_vlm
                # f_cl = clamp(f_cl, 0.0, 1.0)
                # @show j, cl_star / cl_vlm, cl_star, cl_vlm, polar.m_inv, polar.cl_alpha0

                # get direction of viscous drag
                v_induced = R * v_induced # rotate into this frame
                dhat = SVector(v_induced[1], 0.0, v_induced[3])
                dhat /= norm(dhat)

                # get viscous drag coefficient
                # cd = FLOWMath.linear(polar.cls_visc, polar.cds_visc, cl_star)
                cd = FLOWMath.linear(polar.alphas, polar.cds_visc, α_eff)

                # get magnitude of viscous drag
                d_viscous_mag = cd * l_2d_norm * l_2d_norm / (2 * RHO * γ * γ * c)

                # get viscous drag vector
                d_viscous = d_viscous_mag * dhat

                # distribute evenly over chordwise panels
                d_viscous /= size(surface, 1)

                # rotate back into global frame
                d_viscous = Rp * d_viscous

                # apply viscous corrections to aerodynamic force on each panel
                for i in axes(surface, 1)
                    # unpack props
                    (; gamma, velocity, cfb, cfl, cfr, velocity_from_streamwise) = props[i,j]

                    # decompose cfb into 2-d and remaining components
                    cfb_2d = R * cfb
                    # cfb_strip = SVector{3,TF}(cfb_2d[1], 0.0, cfb_2d[3])
                    # cfb_remaining = SVector{3,TF}(0.0, cfb_2d[2], 0.0)

                    # apply lift correction factor to bound circulation contribution and add viscous drag
                    cfb_new = Rp * (SVector{3,TF}(cfb_2d[1], 0.0, cfb_2d[3]) * f_cl + SVector{3,TF}(0, cfb_2d[2], 0)) + d_viscous

                    # reassemble properties
                    # @show cfb_strip, cfb_remaining
                    # @show f_cl, cfb, cfb_new, d_viscous
                    props[i,j] = PanelProperties(
                        gamma,
                        velocity,
                        cfb_new,  # apply lift correction factor to bound circulation contribution and add viscous drag
                        cfl, # * f_cl,  # apply lift correction factor to left edge contribution
                        cfr, # * f_cl,  # apply lift correction factor to right edge contribution
                        velocity_from_streamwise
                    )
                end

                # apply viscous correction to circulation strengths and their time derivatives
                Γ[iΓ-size(surface,1):iΓ-1] .*= f_cl
            end
            dΓdt .+= Γ
            dΓdt ./= dt
        else
            iΓ += size(surface, 1) * size(surface, 2) # skip circulation strengths for this surface
        end
    end
end

function viscous!(properties::Vector{Matrix{PanelProperties{TF}}}, Γ, dΓdt, surfaces::Vector{Matrix{SurfacePanel{TF}}}, grids, frames::Vector{<:ReferenceFrame}, frames_index::Vector{Int}, polar::Nothing, ref) where TF
    return nothing
end
