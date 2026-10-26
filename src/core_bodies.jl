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

"""
    SurfaceBody(system, i)

Surface `i` of `system` as a body of FLOWUnsteadyCore's time loop. A handle: the
geometry stays in `system`, and moving the body moves `system.grids[i]`.
"""
struct SurfaceBody{TF} <: AbstractBody
    system::VortexLatticeSystem{TF}
    index::Int
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
Base.eltype(::SurfaceState{TF}) where TF = TF

function SurfaceState(system::VortexLatticeSystem{TF}, i::Integer) where TF
    nc, ns = size(system.surfaces[i])
    te = _trailing_edge_nodes(system.surfaces[i])
    return SurfaceState{TF}(system, i, zeros(SVector{3,TF}, 2nc*ns + nc*(ns+1)),
        te[1:ns], te[2:ns+1], te[1:ns], te[2:ns+1])
end

"""
    surface_bodies(system) -> (bodies, states)

One [`SurfaceBody`](@ref) and one [`SurfaceState`](@ref) per surface of `system`,
for FLOWUnsteadyCore's `simulate!`. Builds the panels from the grids, and sets the
system up for the core: full surfaces (no mirror image) and no trailing vortices
(the core's wake replaces them). `system.reference[]` must be set; its velocity
nondimensionalizes the panel forces.
"""
function surface_bodies(system::VortexLatticeSystem)
    system.symmetric .= false
    system.trailing_vortices .= false
    for i in eachindex(system.surfaces)
        update_surface_panels!(system.surfaces[i], system.grids[i]; ratios = system.ratios[i],
            fcore = (c, Δs) -> system.core_size)
    end
    n = length(system.surfaces)
    return [SurfaceBody(system, i) for i in 1:n], [SurfaceState(system, i) for i in 1:n]
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

    for (s, ext) in zip(states, externals)
        s.external === ext || copyto!(s.external, ext)
        system.wake_shedding_locations[s.index][1:end-1] .= s.shed_left
        system.wake_shedding_locations[s.index][end] = s.shed_right[end]
    end
    perm = sortperm([s.index for s in states])
    ordered = collect(states)[perm]
    closed = [implicit_shed isa Bool ? implicit_shed : Bool(implicit_shed[k]) for k in perm]

    off = fill(false, nsurf)
    influence_coefficients!(system.AIC, system.surfaces;
        symmetric = off, surface_id = system.surface_id,
        wake_shedding_locations = system.wake_shedding_locations,
        trailing_vortices = off, force_finite_core = off, xhat = system.xhat[])
    _open_trailing_edges!(system, ordered, closed)
    system.fAIC[] = lu(system.AIC)

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

function SixDOF.move(b::SurfaceBody, origin, R, dx)
    system = b.system
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

FLOWUnsteadyCore.body_state(b::SurfaceBody) = Dict{Symbol,Any}(:grid => copy(b.system.grids[b.index]))

function FLOWUnsteadyCore.restore_body(b::SurfaceBody, d::Dict)
    system = b.system
    system.grids[b.index] .= d[:grid]
    update_surface_panels!(system.surfaces[b.index], system.grids[b.index]; ratios = system.ratios[b.index],
        fcore = (c, Δs) -> system.core_size)
    return b
end

function FLOWUnsteadyCore.vtk_body!(vtm, s::SurfaceState, b::SurfaceBody)
    write_vtk!(vtm, b.system.surfaces[b.index], s.system.properties[b.index];
        trailing_edge = true, trailing_vortices = false, symmetric = false)
    return nothing
end
