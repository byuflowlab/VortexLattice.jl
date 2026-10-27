#=
VortexLattice on FLOWUnsteadyCore's time loop. A `VortexLatticeSystem` keeps its
surfaces, its AIC and its solve; the core addresses one surface at a time through
two handles that point into the system, `SurfaceBody` (the geometry) and
`SurfaceState` (the external velocity and the shedding line). All surfaces of one
system are one solve group, solved together by the system's own AIC, and the core
never splits it.

Evaluation points of a surface (nc x ns panels), in order and each block in
column-major order: the control points (nc x ns), the bound-vortex centres
(nc x ns) and the vertical-segment centres (nc x (ns+1)). The near field reads
the last two blocks as its `Vh` and `Vv`.

The trailing-edge rings end on the core's shedding line (VortexLattice's
`wake_shedding_locations`, no trailing vortices), so each ring's closing edge lies
along the line the core's first wake row starts from. The core evaluates that
row's top edge (the previous circulation) on the surface and the closing edge
carries the current one, so the pair nets to the shed change (the core's implicit
shedding, with buffered rows). Without buffered rows the closing edges come off
and the legs to the shedding line stay open.
=#

# The AIC depends only on the surfaces' relative geometry, so it is built again only
# when they moved relative to each other; the trailing-edge columns follow the shedding
# line every solve, and the LU is refactored only when they change.
mutable struct AICCache{TF}
    valid::Bool
    te::Matrix{TF}     # the trailing-edge columns the LU was factored with
end

"""
    SurfaceBody

Surface `i` of a `VortexLatticeSystem` as a body of FLOWUnsteadyCore's time loop
(built by [`surface_bodies`](@ref)). A handle: the geometry stays in the system,
and moving the body moves `system.grids[i]`. It carries the surface's polars (one
per spanwise strip, or `nothing` for no viscous correction) and its axes in global
coordinates, which turn with it.
"""
struct SurfaceBody{TF,P} <: AbstractBody
    system::VortexLatticeSystem{TF}
    index::Int
    polars::P
    basis::Base.RefValue{SMatrix{3,3,TF,9}}
    motion::Base.RefValue{Tuple{SMatrix{3,3,TF,9},SVector{3,TF}}}   # rigid motion since the AIC was built
    aic::AICCache{TF}                                              # shared by the system's surfaces
end


"""
    SurfaceState(system, i)

The per-step state of surface `i` of `system`: the velocity the core accumulates
at its evaluation points and its shedding line. Its circulation is the slice of
`system.Γ` that belongs to surface `i`.
"""
struct SurfaceState{TF} <: AbstractBodyState
    system::VortexLatticeSystem{TF}
    index::Int
    external::Vector{SVector{3,TF}}
    te_left_prev::Vector{SVector{3,TF}}
    te_right_prev::Vector{SVector{3,TF}}
    shed_left::Vector{SVector{3,TF}}
    shed_right::Vector{SVector{3,TF}}
end

Base.eltype(::SurfaceBody{TF}) where TF = TF
SurfaceBody(system::VortexLatticeSystem{TF}, i) where TF = SurfaceBody(system, i, nothing,
    Ref(SMatrix{3,3,TF,9}(I)), Ref((SMatrix{3,3,TF,9}(I), zero(SVector{3,TF}))), AICCache(false, zeros(TF, 0, 0)))
Base.eltype(::SurfaceState{TF}) where TF = TF

function SurfaceState(system::VortexLatticeSystem{TF}, i::Integer) where TF
    nc, ns = size(system.surfaces[i])
    te = _trailing_edge_nodes(system.surfaces[i])
    return SurfaceState{TF}(system, i, zeros(SVector{3,TF}, 2nc*ns + nc*(ns+1)),
        te[1:ns], te[2:ns+1], te[1:ns], te[2:ns+1])
end

"""
    surface_bodies(system; polars = nothing, basis = I) -> (bodies, states)

One [`SurfaceBody`](@ref) and one [`SurfaceState`](@ref) per surface of `system`,
for FLOWUnsteadyCore's `simulate!`. Builds the panels from the grids, and sets the
system up for the core: full surfaces (no mirror image) and no trailing vortices
(the core's wake replaces them). `system.reference[]` must be set; its velocity
nondimensionalizes the panel forces.

`polars`: per surface, one polar per spanwise strip (or `nothing`) for the viscous
correction. `basis`: per surface (or one for all), the surface's chordwise,
spanwise and normal axes in global coordinates as matrix columns; the identity
suits a wing laid out along y.
"""
function surface_bodies(system::VortexLatticeSystem{TF}; polars = nothing, basis = I) where TF
    system.symmetric .= false
    system.trailing_vortices .= false
    for i in eachindex(system.surfaces)
        update_surface_panels!(system.surfaces[i], system.grids[i]; ratios = system.ratios[i],
            fcore = (c, Δs) -> system.core_size)
    end
    n = length(system.surfaces)
    aic = AICCache(false, zeros(TF, 0, 0))
    pick(x, i) = x isa AbstractVector && !(x isa AbstractVector{<:Number}) ? x[i] : x
    bodies = [SurfaceBody(system, i, isnothing(polars) ? nothing : polars[i],
                          Ref(SMatrix{3,3,TF,9}(pick(basis, i))), Ref((SMatrix{3,3,TF,9}(I), zero(SVector{3,TF}))), aic)
              for i in 1:n]
    return bodies, [SurfaceState(system, i) for i in 1:n]
end

# the trailing-edge vertices of a surface, ns+1 of them
_trailing_edge_nodes(surface) =
    [j <= size(surface, 2) ? bottom_left(surface[end, j]) : bottom_right(surface[end, end]) for j in 1:size(surface, 2)+1]

# indices of surface i's panels in system.Γ
function _panel_range(system, i)
    k = 0
    for j in 1:i-1
        k += length(system.surfaces[j])
    end
    return k+1:k+length(system.surfaces[i])
end

# --- evaluation points and state ---

function FLOWUnsteadyCore.n_evaluation_points(b::SurfaceBody)
    nc, ns = size(b.system.surfaces[b.index])
    return 2nc*ns + nc*(ns+1)
end

function FLOWUnsteadyCore.evaluation_points!(dest, b::SurfaceBody)
    surface = b.system.surfaces[b.index]
    nc, ns = size(surface)
    n = nc*ns
    for k in 1:n
        dest[k] = controlpoint(surface[k])
        dest[n+k] = top_center(surface[k])
    end
    for j in 1:ns+1, i in 1:nc
        dest[2n + (j-1)*nc + i] = j <= ns ? left_center(surface[i, j]) : right_center(surface[i, ns])
    end
    return dest
end

FLOWUnsteadyCore.external_velocity(s::SurfaceState) = s.external

# what sheds: the trailing-edge row of the surface's circulation
function FLOWUnsteadyCore.circulation(s::SurfaceState)
    nc, ns = size(s.system.surfaces[s.index])
    r = _panel_range(s.system, s.index)
    return view(s.system.Γ, first(r)-1 .+ (nc:nc:nc*ns))
end

# what the coupling iterates on: every panel of the surface
FLOWUnsteadyCore.coupling_state(s::SurfaceState) = view(s.system.Γ, _panel_range(s.system, s.index))

function FLOWUnsteadyCore.solved_velocity(s::SurfaceState)
    v = copy(s.external)
    FLOWUnsteadyCore.self_influence!(v, s, SurfaceBody(s.system, s.index), NamedTuple(); closed = true)
    return v
end

# --- shedding ---

FLOWUnsteadyCore.n_shedding_elements(b::SurfaceBody) = size(b.system.surfaces[b.index], 2)
FLOWUnsteadyCore.shedding_locations(s::SurfaceState, ::SurfaceBody) = (s.shed_left, s.shed_right)
FLOWUnsteadyCore.previous_shedding_locations(s::SurfaceState, ::SurfaceBody) = (s.te_left_prev, s.te_right_prev)

function FLOWUnsteadyCore.store_trailing_edges!(s::SurfaceState, b::SurfaceBody)
    te = _trailing_edge_nodes(b.system.surfaces[b.index])
    ns = length(te) - 1
    copyto!(s.te_left_prev, view(te, 1:ns)); copyto!(s.te_right_prev, view(te, 2:ns+1))
    copyto!(s.shed_left, view(te, 1:ns)); copyto!(s.shed_right, view(te, 2:ns+1))
    return nothing
end

# the shedding points a fraction eta of one step behind the trailing edge, along the
# relative flow: TE + eta * (TE_prev + vinf * dt - TE)
function FLOWUnsteadyCore.set_shedding_locations!(s::SurfaceState{TF}, b::SurfaceBody, vinf, dt, eta;
        from::Symbol = :trailing_edge) where TF
    from in (:trailing_edge, :convected) || throw(ArgumentError(
        "VortexLattice sheds from the trailing edge; shed_from = $from is not supported"))
    te = _trailing_edge_nodes(b.system.surfaces[b.index])
    e = TF(eta)
    step(X) = SVector{3,TF}(vinf isa Function ? vinf(X) : vinf) * TF(dt)
    for j in eachindex(s.shed_left)
        pl = s.te_left_prev[j]; pr = s.te_right_prev[j]
        s.shed_left[j] = te[j] + e * (pl + step(pl) - te[j])
        s.shed_right[j] = te[j+1] + e * (pr + step(pr) - te[j+1])
    end
    return s
end

# --- dΓ/dt and the unsteady force ---

function FLOWUnsteadyCore.begin_dgamma_dt!(s::SurfaceState)
    r = _panel_range(s.system, s.index)
    s.system.dΓdt[r] .= .-view(s.system.Γ, r)
    return nothing
end

# dΓ/dt over the step, and its force rho * dΓ/dt * panel area along the panel normal
# (Katz & Plotkin), added to the panel forces the solve left
function FLOWUnsteadyCore.finish_dgamma_dt!(s::SurfaceState, dt)
    system = s.system
    r = _panel_range(system, s.index)
    dΓdt = view(system.dΓdt, r)
    dΓdt .= (dΓdt .+ view(system.Γ, r)) ./ dt
    ref = system.reference[]
    q = 1/2*RHO*ref.V^2
    surface = system.surfaces[s.index]
    props = system.properties[s.index]
    cr = CartesianIndices(surface)
    for k in eachindex(surface)
        I = cr[k]
        dΓdtk = I[1] == 1 ? dΓdt[k] : (dΓdt[k] + dΓdt[k-1])/2
        p = props[k]
        F = RHO*dΓdtk*surface[k].chord*norm(top_vector(surface[k]))*surface[k].ncp
        props[k] = PanelProperties(p.gamma, p.velocity, p.cfb + F/(q*ref.S), p.cfl, p.cfr)
    end
    return nothing
end

# --- the solve: every surface of the system at once ---

"""
    VortexLatticeGroup

The solve-group key of the surfaces of one `VortexLatticeSystem`: they share one
AIC and are only solved together.
"""
struct VortexLatticeGroup
    id::UInt
end

FLOWUnsteadyCore.solve_group(b::SurfaceBody) = VortexLatticeGroup(objectid(b.system))
FLOWUnsteadyCore.indivisible_group(::VortexLatticeGroup) = true

FLOWUnsteadyCore.solve_body!(s::SurfaceState, b::SurfaceBody, freestream, external, options, implicit_shed;
        workspace = nothing) =
    FLOWUnsteadyCore.solve_group!([s], [b], freestream, [external], options, [implicit_shed])

function FLOWUnsteadyCore.solve_group!(states::AbstractVector{<:SurfaceState}, bodies::AbstractVector,
        freestream, externals, options, implicit_shed; workspace = nothing, timing = nothing)
    system = first(states).system
    TF = eltype(system)
    nsurf = length(system.surfaces)
    sort([s.index for s in states]) == 1:nsurf && all(s -> s.system === system, states) ||
        throw(ArgumentError("a VortexLatticeSystem is solved whole: the group must hold every surface of one system"))
    all(b -> b.system === system, bodies) || throw(ArgumentError(
        "the bodies point at a different VortexLatticeSystem than their states (a deepcopy detached them)"))

    for (s, ext) in zip(states, externals)
        s.external === ext || copyto!(s.external, ext)
        system.wake_shedding_locations[s.index][1:end-1] .= s.shed_left
        system.wake_shedding_locations[s.index][end] = s.shed_right[end]
    end
    perm = sortperm([s.index for s in states])
    ordered = collect(states)[perm]
    closed = [implicit_shed isa Bool ? implicit_shed : Bool(implicit_shed[k]) for k in perm]

    off = fill(false, nsurf)
    _update_aic!(system, bodies, ordered, closed, off)

    # the external velocity (freestream, motion, wake, other groups) through each control point
    iw = 0
    for s in ordered
        surface = system.surfaces[s.index]
        for k in eachindex(surface)
            system.w[iw+k] = -dot(s.external[k], normal(surface[k]))
        end
        iw += length(surface)
    end
    ldiv!(system.Γ, system.fAIC[], system.w)

    # the panel forces: the external velocity at the bound and vertical centres plus
    # the surfaces' own induction (the freestream is already in the external velocity)
    Vh = [_bound_centre_velocity(s) for s in ordered]
    all(closed) || (Vh = _open_trailing_edge_velocity(Vh, system, ordered, closed))
    Vv = [_vertical_centre_velocity(s) for s in ordered]
    still = BodyFreestream(zero(TF), zero(TF), zero(TF), zero(SVector{3,TF}))
    near_field_forces!(system.properties, system.surfaces, system.wakes, system.reference[], still, system.Γ;
        dΓdt = nothing, additional_velocity = nothing, Vh, Vv, symmetric = off, nwake = zeros(Int, nsurf),
        surface_id = system.surface_id, wake_finite_core = system.wake_finite_core,
        wake_shedding_locations = system.wake_shedding_locations, trailing_vortices = off,
        xhat = system.xhat[], calculate_vlm_induced = true)

    # the sectional polars' correction, surface by surface
    for s in ordered
        b = bodies[findfirst(b -> b.index == s.index, bodies)]
        isnothing(b.polars) && continue
        viscous!(system.properties[s.index], view(system.Γ, _panel_range(system, s.index)), system.surfaces[s.index],
            system.grids[s.index], b.polars, b.basis[], system.reference[])
    end

    # the flight state body_forces turns coefficients into the stability and wind frames with
    system.freestream[] = velocity_to_freestream(SVector{3,TF}(freestream.velocity))
    system.near_field_analysis[] = true
    system.derivatives[] = false
    return nothing
end

# the external velocity at the bound-vortex centres (nc x ns) and the vertical-segment
# centres (nc x (ns+1)), as the near field's Vh and Vv
function _bound_centre_velocity(s::SurfaceState)
    n = length(s.system.surfaces[s.index])
    return reshape(view(s.external, n+1:2n), size(s.system.surfaces[s.index]))
end
function _vertical_centre_velocity(s::SurfaceState)
    nc, ns = size(s.system.surfaces[s.index])
    return reshape(view(s.external, 2nc*ns+1:2nc*ns+nc*(ns+1)), nc, ns+1)
end

# The AIC for the current shedding line: built whole when the surfaces moved relative to each
# other since the last build (or never were), else only its trailing-edge columns, and the LU
# refactored only when those changed.
function _update_aic!(system::VortexLatticeSystem{TF}, bodies, ordered, closed, off) where TF
    cache = first(bodies).aic
    R1, t1 = first(bodies).motion[]
    scale = one(TF) + maximum(maximum(abs, g) for g in system.grids)
    rigid = all(b -> (Rb = b.motion[]; maximum(abs, Rb[1] - R1) <= 1e-12 && maximum(abs, Rb[2] - t1) <= 1e-12 * scale), bodies)
    for b in bodies
        b.motion[] = (SMatrix{3,3,TF,9}(I), zero(SVector{3,TF}))
    end
    te = [first(_panel_range(system, i)) - 1 + (j-1)*size(system.surfaces[i], 1) + size(system.surfaces[i], 1)
          for i in eachindex(system.surfaces) for j in 1:size(system.surfaces[i], 2)]
    kw = (; symmetric = off, surface_id = system.surface_id, wake_shedding_locations = system.wake_shedding_locations,
            trailing_vortices = off, force_finite_core = off, xhat = system.xhat[])
    if cache.valid && rigid
        update_trailing_edge_coefficients!(system.AIC, system.surfaces; kw...)
        _open_trailing_edges!(system, ordered, closed)
        maximum(abs, view(system.AIC, :, te) - cache.te) <= 1e-12 * maximum(abs, cache.te) && return system
    else
        influence_coefficients!(system.AIC, system.surfaces; kw...)
        _open_trailing_edges!(system, ordered, closed)
        cache.valid = true
    end
    system.fAIC[] = lu(system.AIC)
    cache.te = system.AIC[:, te]
    return system
end

# The closing edges along the shedding line are part of influence_coefficients! and of the
# near field's own induction. Without implicit shedding (no buffered wake row to pair them
# with) they come off again, leaving the legs from the trailing edge to the shedding line
# open; the kernel and core are the ones influence_coefficients! used for that edge.
function _closing_edge(r, s::SurfaceState, j, finite_core)
    surface = s.system.surfaces[s.index]
    return bound_induced_velocity(r - s.shed_right[j], r - s.shed_left[j], finite_core, get_core_size(surface[end, j]))
end

function _open_trailing_edges!(system, ordered, closed)
    row = 0
    for sr in ordered
        receiving = system.surfaces[sr.index]
        for k in eachindex(receiving)
            rcp = controlpoint(receiving[k]); nhat = normal(receiving[k])
            for (ss, cl) in zip(ordered, closed)
                cl && continue
                nc, ns = size(system.surfaces[ss.index])
                col = first(_panel_range(system, ss.index)) - 1
                fc = system.surface_id[sr.index] != system.surface_id[ss.index]
                for j in 1:ns
                    system.AIC[row+k, col+(j-1)*nc+nc] -= dot(_closing_edge(rcp, ss, j, fc), nhat)
                end
            end
        end
        row += length(receiving)
    end
    return system
end

function _open_trailing_edge_velocity(Vh, system, ordered, closed)
    Vh = [copy(v) for v in Vh]
    for (i, sr) in enumerate(ordered)
        receiving = system.surfaces[sr.index]
        for k in eachindex(receiving)
            rc = top_center(receiving[k])
            for (ss, cl) in zip(ordered, closed)
                cl && continue
                G = FLOWUnsteadyCore.circulation(ss)
                fc = system.surface_id[sr.index] != system.surface_id[ss.index]
                for j in eachindex(G)
                    Vh[i][k] -= G[j] * _closing_edge(rc, ss, j, fc)
                end
            end
        end
    end
    return Vh
end

# --- monitor ---

"""
    ForcesMonitor(system, nsteps; frame = Body(), freestream = nothing)

A FLOWUnsteadyCore `simulate!` monitor that records `system`'s force and moment
coefficients ([`body_forces`](@ref) in `frame`) after each step's solve, in
`CF[step]` and `CM[step]`. The stability and wind axes follow the relative wind
at the reference point: the freestream (the run's `freestream` at the step's
time when it is given, else its reference velocity) minus the system's frame
motion there, which also sets the body rates.
"""
struct ForcesMonitor{TF,F,FS}
    system::VortexLatticeSystem{TF}
    CF::Vector{SVector{3,TF}}
    CM::Vector{SVector{3,TF}}
    frame::F
    freestream::FS
end

ForcesMonitor(system::VortexLatticeSystem{TF}, nsteps::Integer; frame = Body(), freestream = nothing) where TF =
    ForcesMonitor(system, zeros(SVector{3,TF}, nsteps), zeros(SVector{3,TF}, nsteps), frame, freestream)

function (m::ForcesMonitor)(step, t, systems, geometries, wake, frames)
    s = m.system
    TF = eltype(s)
    ref = s.reference[]
    rref = SVector{3,TF}(ref.r)
    u = isnothing(m.freestream) ? freestream_velocity(s.freestream[]) :
        SVector{3,TF}(FLOWUnsteadyCore.freestream_velocity(m.freestream, rref, t))
    i = findfirst(g -> g isa SurfaceBody && g.system === s, geometries)
    if isnothing(i) || !SixDOF.owns(frames, i)
        vref = zero(SVector{3,TF}); omega = zero(SVector{3,TF})
    else
        origin, v, omega = SixDOF.frame_motion(frames, i)
        vref = v + cross(omega, rref - origin)
    end
    s.freestream[] = velocity_to_freestream(SVector{3,TF}(u - vref), SVector{3,TF}(-omega[1], omega[2], -omega[3]))
    m.CF[step], m.CM[step] = body_forces(s; frame = m.frame)
    return nothing
end

# --- the surface as vortex filaments ---

# Surface i's rings as segments, with the net strength on each shared edge: the
# bound (horizontal) edges, the vertical edges, the legs from the trailing edge to
# the shedding line and, when `closed`, the closing edges along the shedding line.
function _surface_filaments!(f::FilamentSource, s::SurfaceState, b::SurfaceBody, core, closed::Bool)
    surface = b.system.surfaces[b.index]
    nc, ns = size(surface)
    G = reshape(view(b.system.Γ, _panel_range(b.system, b.index)), nc, ns)
    seg!(r1, r2, g) = FLOWUnsteadyCore._push_segment!(f, r1, r2, g, FLOWUnsteadyCore._segment_core(core, r1, r2))
    for j in 1:ns, k in 1:nc
        p = surface[k, j]
        seg!(p.rtl, p.rtr, G[k, j] - (k > 1 ? G[k-1, j] : zero(eltype(G))))
    end
    for jj in 1:ns+1, k in 1:nc
        top = jj <= ns ? surface[k, jj].rtl : surface[k, ns].rtr
        bot = jj <= ns ? surface[k, jj].rbl : surface[k, ns].rbr
        seg!(top, bot, (jj > 1 ? G[k, jj-1] : zero(eltype(G))) - (jj <= ns ? G[k, jj] : zero(eltype(G))))
    end
    te = _trailing_edge_nodes(surface)
    for jj in 1:ns+1
        shed = jj <= ns ? s.shed_left[jj] : s.shed_right[ns]
        seg!(te[jj], shed, (jj > 1 ? G[nc, jj-1] : zero(eltype(G))) - (jj <= ns ? G[nc, jj] : zero(eltype(G))))
    end
    if closed
        for j in 1:ns
            seg!(s.shed_right[j], s.shed_left[j], G[nc, j])
        end
    end
    return f
end

function FLOWUnsteadyCore.coupling_sources(states::Union{AbstractVector{<:SurfaceState},Tuple{Vararg{SurfaceState}}},
        bodies, options::NamedTuple; closed = false)
    TF = isempty(states) ? Float64 : eltype(first(states))
    f = FilamentSource(TF)
    core = get(options, :filament_core, 0)
    for i in eachindex(states)
        _surface_filaments!(f, states[i], bodies[i], core, closed isa Bool ? closed : closed[i])
    end
    return f
end

FLOWUnsteadyCore.append_body_filaments!(dest, s::SurfaceState, b::SurfaceBody, options; closed::Bool = false, core = 0) =
    _surface_filaments!(dest, s, b, core, closed)

FLOWUnsteadyCore.body_source(s::SurfaceState, b::SurfaceBody; closed::Bool = false, kwargs...) =
    _surface_filaments!(FilamentSource(eltype(s)), s, b, 0, closed)

function FLOWUnsteadyCore.body_induced_velocity!(dest, points, s::SurfaceState, b::SurfaceBody, options;
        closed::Bool = false)
    f = _surface_filaments!(FilamentSource(eltype(s)), s, b, get(options, :filament_core, 0), closed)
    FLOWUnsteadyCore._add_influence!(dest, points, f, DirectBackend(), false)
    return dest
end

# its own filaments at its own points, which the coupled sweep subtracts from the all-bodies pass:
# the same segments as coupling_sources, so the subtraction is exact
FLOWUnsteadyCore.self_influence!(dest, s::SurfaceState, b::SurfaceBody, options; closed::Bool = false) =
    FLOWUnsteadyCore.body_induced_velocity!(dest, evaluation_points!(similar(s.external), b), s, b, options; closed)

# --- motion ---

function SixDOF.move(b::SurfaceBody{TF}, origin, R, dx) where TF
    system = b.system
    Ra, ta = b.motion[]
    b.motion[] = (SMatrix{3,3,TF,9}(R * Ra), SVector{3,TF}(R * ta + origin + dx - R * origin))
    b.basis[] = SMatrix{3,3,TF,9}(R * b.basis[])
    rotate_translate!(system, b.index, origin, R, dx)
    update_surface_panels!(system.surfaces[b.index], system.grids[b.index]; ratios = system.ratios[b.index],
        fcore = (c, Δs) -> system.core_size)
    return b
end

# --- checkpoint and output ---

function FLOWUnsteadyCore.system_state(s::SurfaceState)
    r = _panel_range(s.system, s.index)
    return Dict{Symbol,Any}(:gamma => s.system.Γ[r], :dgamma_dt => s.system.dΓdt[r], :external => copy(s.external),
        :te_left_prev => copy(s.te_left_prev), :te_right_prev => copy(s.te_right_prev),
        :shed_left => copy(s.shed_left), :shed_right => copy(s.shed_right))
end

function FLOWUnsteadyCore.restore_system!(s::SurfaceState, d::Dict)
    r = _panel_range(s.system, s.index)
    s.system.Γ[r] .= d[:gamma]; s.system.dΓdt[r] .= d[:dgamma_dt]
    for k in (:external, :te_left_prev, :te_right_prev, :shed_left, :shed_right)
        copyto!(getfield(s, k), d[k])
    end
    return s
end

FLOWUnsteadyCore.body_state(b::SurfaceBody) = Dict{Symbol,Any}(:grid => copy(b.system.grids[b.index]), :basis => b.basis[])

function FLOWUnsteadyCore.restore_body(b::SurfaceBody, d::Dict)
    system = b.system
    system.grids[b.index] .= d[:grid]
    haskey(d, :basis) && (b.basis[] = d[:basis])
    b.aic.valid = false
    update_surface_panels!(system.surfaces[b.index], system.grids[b.index]; ratios = system.ratios[b.index],
        fcore = (c, Δs) -> system.core_size)
    return b
end

# the geometry lives in the shared system, so a snapshot is the surface's grid and axes,
# written back in place (a deepcopy would detach a copy of the whole system)
FLOWUnsteadyCore.snapshot_body(b::SurfaceBody) =
    (grid = copy(b.system.grids[b.index]), basis = b.basis[], motion = b.motion[])

function FLOWUnsteadyCore.restore_snapshot(b::SurfaceBody, snap)
    system = b.system
    system.grids[b.index] .= snap.grid
    update_surface_panels!(system.surfaces[b.index], system.grids[b.index]; ratios = system.ratios[b.index],
        fcore = (c, Δs) -> system.core_size)
    b.basis[] = snap.basis
    b.motion[] = snap.motion
    return b
end

function FLOWUnsteadyCore.vtk_body!(vtm, s::SurfaceState, b::SurfaceBody)
    write_vtk!(vtm, b.system.surfaces[b.index], s.system.properties[b.index];
        trailing_edge = true, trailing_vortices = false, symmetric = false)
    return nothing
end
