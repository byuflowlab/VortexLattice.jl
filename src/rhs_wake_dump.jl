# RHS_WAKE_DUMP (diagnostic) — ground-truth cross-check of the wake's contribution
# to the RHS, recovered from the abandoned 2026-08-25 viscous-accuracy investigation
# and widened from a single control point to the full span.
#
# The original dump captured only surface 1, chordwise row 1, spanwise station ns÷2.
# The axial-induction deficit is a *spanwise* signal (VL realizes 35-50% of CCBlade's
# induction, varying with r/R), so a single station cannot discriminate "wake too
# weak" from "wake too far" from "evaluated at the wrong location". cp.csv now carries
# every (i,j) on surface 1 with its control point and Vcp, so an independent
# Biot-Savart evaluation over the dumped particle/ring/filament state can be diffed
# against VL's own Vcp station by station.

const RHS_WAKE_DUMP_ENABLED = Ref(false)
const RHS_WAKE_DUMP_STEP = Ref(-1)
const RHS_WAKE_DUMP_DIR = Ref("")
const PARTICLES_USE_GAMMA_WAKE = Ref(false)
# Camber from polars (2026-09-01, default ON): a flat nc=1 panel solves a flat-plate circulation while
# the polar-based strip forces respond to (α − α_L0); the wake then carried ~60% of the force-consistent
# circulation (see research-notes VPM-convergence/2026-09-01). Each panel normal is rotated nose-up by
# −α_L0 about its span axis (update_surface_panels!) and the polar tables are re-based so α is
# measured from the zero-lift line (cl_inv = 2π(α − α_L0), cd looked up at the same shifted α).
const CAMBER_FROM_POLARS = Ref(true)
const CAMBER_ALPHA0 = Vector{Vector{Float64}}()           # per surface, per station, degrees; empty = off
# map strip j of a surface with nj strips onto the per-input-station vectors (length n). Mirrored
# surfaces (wing_to_grid mirror=true) carry 2n strips: 1..n are the reflected half in reverse order.
@inline camber_station(j, n, nj) = nj == 2n ? (j <= n ? n + 1 - j : j - n) : j
# Option 3 (2026-09-01): impose lifting-line Γ from the full polar (FLOWVLM-style), skip the AIC result
const IMPOSE_POLAR_GAMMA = Ref(false)
const IMPOSE_POLAR_GAMMA_POLARS = Ref{Any}(nothing)
const IMPOSE_POLAR_GAMMA_PREV = Ref{Any}(nothing)
const PRESCRIBED_GAMMA = Ref{Any}(nothing)   # H9 comparator (2026-09-02): (r_m, Γ) table; overrides system.Γ every step (nc=1), no blade feedback
const _camber_isurf = Ref(1)
# DYNAMIC CAMBER (2026-09-02): per-strip extra normal rotation Δα (deg) so the lattice produces the
# polar's lift itself (Γ_pre → Γ_wake); updated in viscous! with under-relaxation. Opt-in.
const DYNAMIC_CAMBER = Ref(true)     # default on (Bruno, 2026-09-02 pm): mode C, closes the shed-Γ deficit (H12)
const DYNAMIC_CAMBER_OMEGA = Ref(0.3)
const DYNAMIC_CAMBER_ALPHA_V = Ref(true)   # α from ¼-chord velocity vs geometric chord (mode B); false = cl_vlm/2π − Δα (mode A, unstable)
const DYNAMIC_CAMBER_LOG = Ref(false)
const DRAG_ABS_SPAN = Ref(true)    # 2026-09-02 (H16): drag uses |Δs_y|; Δs_y is signed for the lift term and <0 on NREL 5MW, which flipped drag into spurious +Tp
const DRAG_SIGN_LOG = Ref(false)
const DYNAMIC_CAMBER_START = Ref(36)     # hold Δα = 0 for this many viscous! calls (impulsive-start transient)
const _dyncamber_calls = Ref(0)
const CAMBER_DALPHA = Vector{Vector{Float64}}()
# H10c (2026-09-02, opt-in): velocity used for the Kutta-Joukowski force = (1-β)·V(¼-chord bound midpoint,
# all sources) + β·V(¾-chord control point, all sources except the panel's own leading-edge segment).
# β = 0.5 approximates the half-chord onset flow without the Kutta-condition self-influence.
const KJ_HALFCHORD_BETA = Ref(0.0)   # default 0 (Bruno, 2026-09-02 pm): 0.5 was resolution-dependent at hi-fi (−39% Tp at the tip); tunable

function _dump_rhs_wake_state(system, wake, i_step, vcp_kin)
    isurf = 1
    surface = system.surfaces[isurf]
    nc, ns = size(surface)

    dir = RHS_WAKE_DUMP_DIR[]
    mkpath(dir)

    # widened: every control point on surface 1, not just (1, ns÷2)
    open(joinpath(dir, "cp.csv"), "w") do io
        println(io, "i,j,x,y,z,vcp_x,vcp_y,vcp_z,vkin_x,vkin_y,vkin_z,nx,ny,nz")
        for j in 1:ns, i in 1:nc
            rcp = controlpoint(surface[i, j])
            vcp = system.Vcp[isurf][i, j]
            vk = vcp_kin[i, j]
            n = normal(surface[i, j])
            println(io, "$i,$j,$(rcp[1]),$(rcp[2]),$(rcp[3]),$(vcp[1]),$(vcp[2]),$(vcp[3]),$(vk[1]),$(vk[2]),$(vk[3]),$(n[1]),$(n[2]),$(n[3])")
        end
    end

    np = FLOWVPM.get_np(wake.pfield)
    open(joinpath(dir, "particles.csv"), "w") do io
        println(io, "x,y,z,gx,gy,gz,sigma")
        for p in 1:np
            X = FLOWVPM.get_X(wake.pfield, p)
            G = FLOWVPM.get_Gamma(wake.pfield, p)
            s = FLOWVPM.get_sigma(wake.pfield, p)[]
            println(io, "$(X[1]),$(X[2]),$(X[3]),$(G[1]),$(G[2]),$(G[3]),$s")
        end
    end

    open(joinpath(dir, "rings.csv"), "w") do io
        println(io, "rtl_x,rtl_y,rtl_z,rtr_x,rtr_y,rtr_z,rbl_x,rbl_y,rbl_z,rbr_x,rbr_y,rbr_z,core_size,gamma")
        for k in eachindex(wake.wakes)
            wk = wake.wakes[k]
            nwk = wake.nwake[k]
            nc_w, nsk = size(wk)
            for jj in 1:nsk, ii in 1:min(nwk, nc_w)
                p = wk[ii, jj]
                println(io, "$(p.rtl[1]),$(p.rtl[2]),$(p.rtl[3]),$(p.rtr[1]),$(p.rtr[2]),$(p.rtr[3]),$(p.rbl[1]),$(p.rbl[2]),$(p.rbl[3]),$(p.rbr[1]),$(p.rbr[2]),$(p.rbr[3]),$(p.core_size),$(p.gamma)")
            end
        end
    end

    bfw = wake.boundary_filaments
    open(joinpath(dir, "boundary_filaments.csv"), "w") do io
        println(io, "r1_x,r1_y,r1_z,r2_x,r2_y,r2_z,core_size,gamma")
        for is in eachindex(bfw.active)
            bfw.active[is] || continue
            for j in eachindex(bfw.r1[is])
                r1 = bfw.r1[is][j]
                r2 = bfw.r2[is][j]
                println(io, "$(r1[1]),$(r1[2]),$(r1[3]),$(r2[1]),$(r2[2]),$(r2[3]),$(bfw.core_size[is][j]),$(bfw.gamma[is][j])")
            end
        end
    end

    open(joinpath(dir, "bound_panels_all.csv"), "w") do io
        println(io, "isurf,j,rtl_x,rtl_y,rtl_z,rtr_x,rtr_y,rtr_z,rbl_x,rbl_y,rbl_z,rbr_x,rbr_y,rbr_z,wsl_l_x,wsl_l_y,wsl_l_z,wsl_r_x,wsl_r_y,wsl_r_z,gamma,core")
        kk = 0
        for (is, sf) in enumerate(system.surfaces)
            ncs, nss = size(sf); wsl = system.wake_shedding_locations[is]
            for j in 1:nss
                pn = sf[ncs, j]; kk += ncs; g = system.Γ[kk]
                println(io, "$is,$j,$(pn.rtl[1]),$(pn.rtl[2]),$(pn.rtl[3]),$(pn.rtr[1]),$(pn.rtr[2]),$(pn.rtr[3]),$(pn.rbl[1]),$(pn.rbl[2]),$(pn.rbl[3]),$(pn.rbr[1]),$(pn.rbr[2]),$(pn.rbr[3]),$(wsl[j][1]),$(wsl[j][2]),$(wsl[j][3]),$(wsl[j+1][1]),$(wsl[j+1][2]),$(wsl[j+1][3]),$g,$(pn.core_size)")
            end
        end
    end

    nbf = FastMultipole.get_n_bodies(wake.boundary_filaments)

    open(joinpath(dir, "meta.txt"), "w") do io
        println(io, "i_step=$i_step np=$np nbf=$nbf nc=$nc ns=$ns")
    end

    println("RHS_WAKE_DUMP written to $dir at step $i_step (np=$np, nbf=$nbf, ncp=$(nc*ns))")
end

function _dump_gamma_correction(system, Γ_wake, ref, i_step)
    surface = system.surfaces[1]
    nc, ns = size(surface)
    props = system.properties[1]
    open(joinpath(RHS_WAKE_DUMP_DIR[], "gamma_correction_surf1.csv"), "w") do io
        println(io, "i,j,gamma_pre,gamma_wake,vx,vy,vz,cfb_x,cfb_y,cfb_z,ds_x,ds_y,ds_z,rtc_x,rtc_y,rtc_z,Vh_x,Vh_y,Vh_z,rbl_x,rbl_y,rbl_z,rbr_x,rbr_y,rbr_z,rtl_x,rtl_y,rtl_z,rtr_x,rtr_y,rtr_z")
        for j in 1:ns, i in 1:nc
            k = (j - 1) * nc + i
            v = props[i, j].velocity * ref.V
            cfb = props[i, j].cfb
            ds = top_vector(surface[i, j])
            pn = surface[i, j]; rtc = top_center(pn); vh = system.Vh[1][i, j]
            println(io, "$i,$j,$(system.Γ[k]),$(Γ_wake[k]),$(v[1]),$(v[2]),$(v[3]),$(cfb[1]),$(cfb[2]),$(cfb[3]),$(ds[1]),$(ds[2]),$(ds[3]),$(rtc[1]),$(rtc[2]),$(rtc[3]),$(vh[1]),$(vh[2]),$(vh[3]),$(pn.rbl[1]),$(pn.rbl[2]),$(pn.rbl[3]),$(pn.rbr[1]),$(pn.rbr[2]),$(pn.rbr[3]),$(pn.rtl[1]),$(pn.rtl[2]),$(pn.rtl[3]),$(pn.rtr[1]),$(pn.rtr[2]),$(pn.rtr[3])")
        end
    end
    println("RHS_WAKE_DUMP gamma_correction_surf1.csv written at step $i_step (ref.V=$(ref.V), ref.S=$(ref.S))")
end

# Option 3 helper: replace system.Γ (surface isurf, nc=1) by the polar lifting-line circulation
# from the bound-midpoint velocity already stored in system.properties by near_field_forces!.
function _impose_polar_gamma!(system, ref)
    polars = IMPOSE_POLAR_GAMMA_POLARS[]
    iΓ = 0
    for isurf in eachindex(system.surfaces)
        surface = system.surfaces[isurf]
        nc, ns = size(surface)
        nc == 1 || error("IMPOSE_POLAR_GAMMA implemented for nc=1 only")
        props = system.properties[isurf]
        for j in 1:ns
            panel = surface[1, j]
            V = props[1, j].velocity * ref.V
            shat = top_vector(panel); shat /= norm(shat)
            chat = 0.5 * (bottom_left(panel) + bottom_right(panel)) - 0.5 * (top_left(panel) + top_right(panel))
            chat -= dot(chat, shat) * shat; chat /= norm(chat)
            nhat = panel.ncp
            Vp = V - dot(V, shat) * shat
            W = norm(Vp)
            # VL's ncp points to the pressure side for this panel convention (the VLM's own Γ is
            # negative at positive incidence here), so incidence and lift are measured toward -ncp
            alpha = atan(-dot(Vp, nhat), dot(Vp, chat)) * 180 / pi
            polar = polars[isurf][j]
            cl = FLOWMath.linear(polar.alphas, polar.cls_visc, alpha)
            lhat = cross(Vp, shat); lhat /= norm(lhat)
            dot(lhat, nhat) > 0 && (lhat = -lhat)
            L = 0.5 * RHO * W^2 * panel.chord * cl * norm(top_vector(panel))
            proj = dot(cross(V, top_vector(panel)), lhat)
            system.Γ[iΓ + j] = L / (RHO * proj)
            if RHS_WAKE_DUMP_ENABLED[] && isurf == 1 && j in (5, 13, 21)
                println("IMPOSE_DBG j=$j alpha=$(round(alpha,digits=2)) W=$(round(W,digits=2)) cl=$(round(cl,digits=3)) proj=$(round(proj,digits=2)) Γ=$(round(system.Γ[iΓ + j],digits=2))")
            end
        end
        iΓ += nc * ns
    end
end

# H9 comparator: replace system.Γ (nc=1) by a prescribed Γ(r) table, r = distance from the x axis of
# the bound-vortex centre. The wake then carries a fixed, known circulation with no blade feedback.
function _impose_prescribed_gamma!(system)
    r_tab, Γ_tab = PRESCRIBED_GAMMA[]
    iΓ = 0
    for isurf in eachindex(system.surfaces)
        surface = system.surfaces[isurf]
        nc, ns = size(surface)
        nc == 1 || error("PRESCRIBED_GAMMA implemented for nc=1 only")
        for j in 1:ns
            rtc = top_center(surface[1, j])
            r = hypot(rtc[2], rtc[3])
            system.Γ[iΓ + j] = FLOWMath.linear(r_tab, Γ_tab, clamp(r, r_tab[1], r_tab[end]))
        end
        iΓ += nc * ns
    end
end

# zero-lift angle (deg) of a polar from its Cl = 0 crossing within ±10°; 0 for lift-less sections
function zero_lift_angle(p::Polar)
    a, cl = p.alphas, p.cls_visc
    m = (a .> -10) .& (a .< 10)
    maximum(cl[m]) <= 0 && return zero(eltype(a))
    return FLOWMath.linear(cl[m], a[m], zero(eltype(a)))
end

# fill CAMBER_ALPHA0 from `polars` and return polars re-based to the zero-lift line
function camber_rebase_polars!(polars)
    empty!(CAMBER_ALPHA0); empty!(CAMBER_DALPHA)
    out = similar(polars)
    for (isurf, ps) in enumerate(polars)
        a0s = [zero_lift_angle(p) for p in ps]
        push!(CAMBER_ALPHA0, a0s)
        push!(CAMBER_DALPHA, zeros(length(a0s)))
        out[isurf] = [Polar(p.alphas .- a0, p.cls_visc, p.cds_visc) for (p, a0) in zip(ps, a0s)]
    end
    return out
end
