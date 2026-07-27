# Bug report: FastMultipole.jl commit `adc4f26` breaks relative-error-tolerance FMM accuracy control, exposing a latent near-field singularity in VortexLattice's `bound_induced_velocity`

## Current status (read this first if picking this up)

- **Applied and working** (see file/change summary at the bottom): FLOWVPM metadata-API port, `induced.jl` `f2`/`f3` fix, `panel_particle_wake.jl` legacy-swap retention + switch-aware `set_gradient!` fix (5 sites), `wake.jl` off-by-one fix, `runtests.jl` scoping fix. Full 36-step `PanelParticleWake` reproduction case runs clean, no `NaN`/crash. Test suite: everything passes except one tight-tolerance assertion in `"FMM induced velocity"` (section F3), believed unrelated and not chased further.
- **Open, actively being worked**: `bound_induced_velocity`'s `f1` term is a genuine, only-partially-solved numerical-methods problem — see section D3 for the full derivation. Two real, physically-distinct near-field singularities exist (corner-approach and midline-approach); the two available formulas (`current`/`legacy`) each fix one and silently break on the other; a derived, synthetically-validated candidate fix (`biv_candidate2`, see `f1_candidate_fix2.jl`) for both **hangs** when run against the real simulation, for reasons not yet root-caused. **This is where to pick up next**: harden/debug that candidate (or find why it hangs) before considering shipping a full `f1` fix. Section D3 has the exact formula, the validation results, and the failure mode.
  - **Ruled out (2026-07-27)**: hypothesized that the hang was caused by the candidate's hard `nr3sq > 1e-28` degenerate-segment cutoff blowing up unboundedly for near-but-not-exactly-degenerate segments (tiny nonzero segment length, target off the line) — tested directly (`f1_hang_hypothesis.jl` in scratchpad) by sweeping segment length `nr3` from `1e-2` down through the `1e-28` cutoff with a fixed off-line target: the candidate's magnitude decreases smoothly toward 0 as `nr3→0` (`0.0031 → 3.4e-9 → 0.0` monotonically), no blowup, no discontinuity at the cutoff. So the cliff at `nr3sq=1e-28` is not itself the hang's cause — need a different hypothesis (e.g. try `rc→0` with `nr3` NOT small — i.e. genuine near-midline approach — combined with some other degeneracy; or instrument the actual hang directly with @time/prints inside the monkey-patched function during the real run rather than guessing from synthetic geometries).
  - **New finding (2026-07-27, later same day)**: instrumented `bound_induced_velocity` (monkey-patched to the candidate) with a call counter + a ring buffer of the last few `(r1, r2, finite_core, core_size)` argument tuples, checkpointing every 500 calls to stderr, and ran the reproduction case with `n_steps` reduced all the way down to **1**. The hang reproduces within the *first* timestep — call counter reaches ~13000-13500 (in two independent runs, consistently in that range) and then stops incrementing entirely, while the Julia process continues to spin at 100% CPU (confirmed via `ps`) for 5+ more minutes with no further calls to `bound_induced_velocity` and no output. This means: (a) the earlier description of the hang manifesting "at step 3-6" was an overestimate — with this candidate it hangs on **timestep 1**, before the crash-triggering `NaN` propagation described elsewhere in this report even becomes possible; (b) since the call counter stops advancing but CPU stays pegged, execution has left `bound_induced_velocity` entirely and is stuck in some *other* function — most likely FastMultipole's own dual-tree traversal / adaptive tree-refinement logic, which would explain a silent infinite (or extremely slow) loop if a single `Inf`/`NaN`-valued (or just pathologically large) induced-velocity evaluation from the candidate feeds into FastMultipole's local error estimation and prevents its refinement criterion from ever being satisfied. This was **not chased further** — confirming it would require instrumenting FastMultipole's tree-traversal/refinement internals directly (a different package, more invasive, higher-effort investigation), which is a natural next step for whoever picks this up but wasn't pursued today given diminishing returns for a known-fragile candidate that was already flagged as unsafe to ship.
  - Also learned from this session: monkey-patching a function low in a call chain (attempted `upper_bound_abs` / `solve_ρ_over_σ` in FLOWVPM to look for a pathological `σ`/`ω`/`ε` feeding a broken `Roots.find_zero` bracket) triggers a very expensive Julia method-invalidation/recompilation cascade through the whole FMM dispatch machinery — one such attempt burned 10+ minutes of 100% CPU with zero output before being killed, with no way to tell recompilation-time apart from a real hang. This ruled out `solve_ρ_over_σ`/`upper_bound_abs` as accessible via this technique — if pursued again, edit the FLOWVPM source file directly (accepting the same one-time recompile cost, but at least isolated from repeated live-`@eval` overhead) rather than monkey-patching at runtime.
- Report sections are numbered roughly in the order they were written (A, B, C, D, D2, D3, E, F, G) — later sections refer back to earlier ones by letter.
- **Upstream issue filed (2026-07-27)**: [byuflowlab/FastMultipole.jl#37](https://github.com/byuflowlab/FastMultipole.jl/issues/37), asking FastMultipole to make the silent relative-error-to-absolute-tolerance fallback (the triggering change, see TL;DR below) louder/opt-out-required rather than a silent per-type default. This is about the *trigger* only — it does not address and is not a substitute for the `f1` fix below, which is a pre-existing VortexLattice defect that would need fixing regardless of what FastMultipole does.
- **Scope reminder**: none of the `f1`-candidate work below (`biv_candidate2` / section D3) has been applied to any tracked file in this repo. It exists only as monkey-patches inside scratch scripts (see section D3 and "Scratchpad files" below) — `git diff` on this repo reflects only the already-applied, already-verified fixes listed above. The repo is currently in a working, non-hung state; only the *next* fix (open item above) risks reintroducing the hang, and only if it's actually wired into `src/induced.jl`.

## TL;DR

`FastMultipole.jl` commit **`adc4f26`** ("support target system metadata and extra outputs in the direct list with accompanying API change", 2026-05-20) removed the `get_previous_influence(system, i)` compatibility hook that downstream packages (FLOWVPM's `ParticleField`) used to implement in order to opt in to relative-error-based FMM accuracy tuning. It was replaced with a new metadata-index API (`previous_potential_metadata_index` / `previous_gradient_metadata_index`). Neither FLOWVPM nor VortexLattice has been ported to the new API, so every FMM call in the VortexLattice `PanelParticleWake` pipeline now silently falls back to **absolute-tolerance-only** accuracy control.

This change in accuracy control is enough to expose a pre-existing, previously-dormant weakness in VortexLattice's finite-core vortex-segment velocity formula (`bound_induced_velocity`, `src/induced.jl:35`): its "desingularized" near-endpoint term does not properly bound the induced velocity as the evaluation point approaches a filament endpoint (it should saturate near `1/core_size`, but instead approaches `1/core_size²`). Under the old (relative+absolute) accuracy regime this pairing was apparently never routed to a code path where this weakness mattered; under the new (absolute-only fallback) regime it is, and the result is a velocity blow-up (~10²–10⁴× too large) at wake-panel corners, which cascades into `NaN` particle circulations/core sizes within 2–4 timesteps and an eventual crash inside `FLOWVPM`'s `solve_ρ_over_σ` (`Roots.find_zero` bracketing failure).

**Verified by bisection:** the identical VortexLattice/FLOWVPM code, run against `FastMultipole` checked out at `c8cf796` (the last commit before `adc4f26`), completes the full 36-step reproduction case with no error. The same code against current `main` (`58cf693`) fails by step 3–4.

---

## Environment

- `VortexLattice.jl` — branch `frames`, HEAD `ee127bc` ("update vpm kwargs pass in")
- `FLOWVPM.jl` — branch `master`, HEAD `e2bd487` ("Bump version to 4.0.4 for SFS/FMM bugfix release"), plus local uncommitted adaptation removing `get_previous_influence` (see below)
- `FastMultipole.jl` — branch `main`, HEAD `58cf693` ("fix multithreading bug with extra_farfield")
- Julia 1.10, macOS
- All three packages are `Pkg.develop`'d path dependencies of each other (standard FLOW lab dev setup)

## Reproduction

`test/profile_panel_particle_wake.jl` in VortexLattice (branch `frames`): a straight, untapered, unswept wing (`b=8`, `c=1`, `ns=13` spanwise panels, `nc=1` chordwise), `alpha=5°`, `Vinf=10`, `dt=0.05`, run through `simulate!` with `wake_type=PanelParticleWake`, `nwakerows=1`, `method_trailing=method_unsteady=OverlapPPS(1.3, 2)`. Calling `run_case()` throws:

```
ERROR: ArgumentError: The interval [a,b] is not a bracketing interval.
You need f(a) and f(b) to have different signs (f(a) * f(b) < 0).
  ...
 [10] solve_ρ_over_σ(...)
    @ FLOWVPM ~/.julia/dev/FLOWVPM/src/FLOWVPM_fmm.jl:42
 ...
 [20] wake_on_all!(...)
    @ VortexLattice ~/.julia/dev/VortexLattice/src/panel_particle_wake.jl:439
```

preceded, on every step, by:

```
┌ Warning: relative error prediction requested but previous influence metadata indices are unavailable for type ParticleField{...}; falling back to absolute tolerance behavior
└ @ FastMultipole ~/.julia/dev/FastMultipole/src/compatibility.jl:424
```

## Investigation timeline

1. **First symptom**: `solve_ρ_over_σ` (FLOWVPM) throws because a particle's circulation magnitude `ω` and core size `σ` are `NaN`. Traced back to a specific wake-buffer-ring panel corner acquiring a huge (~100+) coordinate value one step after being shed, then propagating to `NaN` on the next VPM integration step.

2. **Traced the huge coordinate** to the corner's *velocity* (`system.V`/`wake.wake_velocities`, aliased) being computed as ~2000–8700 (vs. ~0.02–0.4 for freestream + normal induced velocity everywhere else), evaluated at a wake-ring corner via `FastMultipole.direct!` for `VortexLattice.WakeBufferRings` (`src/panel_particle_wake.jl:270-290`), which calls `bound_induced_velocity(target-vi, target-vj, true, cs)` for each of the ring's 4 edges.

3. **Traced the velocity blow-up to a single edge term.** For the edge touching the evaluation point (`target ≈ v4`, i.e. `r1 = target - v4 ≈ (-2.7e-5, 0, 0)`, an ~2.7e-5-length near-zero vector — the probe sits essentially *at* the ring's own corner by construction), `bound_induced_velocity(r1, r2, true, cs=1e-3)` (`src/induced.jl:35-67`) evaluates:

   ```julia
   f2 = 1/(nr1 + δ1)     # δ1 = get_δ(nr1, core_size)
   get_δ(distance, core_size) = distance < core_size ? (distance - core_size)^2 : 0
   ```

   As `nr1 → 0`: `δ1 → core_size²`, so `f2 → 1/core_size²` (≈ 1e6 for `core_size=1e-3`) — **not** the expected, physically-motivated saturation at `1/core_size` (≈ 1e3) that a proper Rankine/Lamb-Oseen-style regularization (`v/(sqrt(r²+ε²))`) would give. In fact at `nr1 = 2.7e-5` (well inside the core), the "regularized" `f2` (≈ 35,943) is barely reduced (~3%) from the *unregularized* `1/nr1` (≈ 37,037) — the desingularization term is nearly inert exactly where it's needed most.

   Hand-computed value for the observed case (`nr1≈2.7e-5`, `nr2≈0.615`, `core_size=1e-3`, `gamma≈1.17`) reproduces the observed order of magnitude (~2700–3000) for this one edge's contribution alone.

4. **Compared against `bound_induced_velocity_legacy`** (same file, an older implementation still present but unused by the new `PanelParticleWake` code paths), which uses the standard `sqrt(r² + ε²)` regularization and gives a well-bounded result (~O(1)) for the identical inputs. Swapping the 6 call sites in `src/panel_particle_wake.jl` from `bound_induced_velocity` to `bound_induced_velocity_legacy` reduced the peak spurious velocity by ~40× (2108 → 55.8) and delayed the crash from step 3–4 to step 6, confirming the mechanism but not eliminating it (a residual, smaller near-singular contribution remains, not yet isolated — likely the same defect reached via `src/fmm.jl`'s separate, more heavily-used copies of the same `bound_induced_velocity` call, which we did **not** swap because of blast radius — see below).

5. **This `bound_induced_velocity` desingularization defect is old** — introduced (or at least present in its current, buggy form) back in commit `752e8e1` ("induced velocity bugfix", 2026-11-18) and effectively unchanged since — yet the reproduction script/feature apparently worked before. This pointed at something *external* changing, not the VortexLattice formula itself being newly broken.

6. **Noticed a suspicious, constant warning** on every simulation step: `"relative error prediction requested but previous influence metadata indices are unavailable ... falling back to absolute tolerance behavior"` (`FastMultipole/src/compatibility.jl:424`). Investigated why: FLOWVPM's `ParticleField` has an **uncommitted, local, in-progress** edit that *removes* its `fmm.get_previous_influence(system::ParticleField, i)` override (`FLOWVPM/src/FLOWVPM_fmm.jl`). Restoring that override and trying to run against current `FastMultipole` main fails to even *compile*:

   ```
   ERROR: LoadError: UndefVarError: `get_previous_influence` not defined
   ```

   i.e., the generic `FastMultipole.get_previous_influence` function that FLOWVPM extends **no longer exists** in current FastMultipole — confirming the local FLOWVPM edit is a (necessary, in-progress) adaptation to an upstream removal, not an independent regression.

7. **Located the removal**: `git log -S"previous_potential_metadata_index" -- .` in FastMultipole → exactly one hit, commit **`adc4f26`**, "support target system metadata and extra outputs in the direct list with accompanying API change" (2026-05-20). Diff excerpt (`src/compatibility.jl`):

   ```diff
   -function get_previous_influence(system, i)
   -    if WARNING_FLAG_MAX_INFLUENCE[]
   -        @warn "get_previous_influence not overloaded for type $(typeof(system)); relative error prediction will not be used"
   -        WARNING_FLAG_MAX_INFLUENCE[] = false
   -    end
   -    return zero(numtype(system)), zero(numtype(system))
   -end
   ```

   replaced by (current `main`, `src/compatibility.jl:395-427`):

   ```julia
   previous_potential_metadata_index(system) = 0
   previous_gradient_metadata_index(system) = 0

   target_uses_previous_influence_metadata(system) =
       previous_potential_metadata_index(system) > 0 && previous_gradient_metadata_index(system) > 0

   function warn_missing_previous_influence_metadata(target_systems::Tuple, error_tolerance)
       isnothing(error_tolerance) && return nothing
       error_tolerance isa RelativeErrorMethod || return nothing
       for system in target_systems
           if !target_uses_previous_influence_metadata(system) && WARNING_FLAG_MAX_INFLUENCE[]
               @warn "relative error prediction requested but previous influence metadata indices are unavailable for type $(typeof(system)); falling back to absolute tolerance behavior"
               ...
   ```

   Since neither `ParticleField` (FLOWVPM) nor any of VortexLattice's FMM-compatible probe/wake types (`WakeBufferRings`, `BoundaryFilamentWrapper`, `PanelBufferFilaments`, `FilamentWrapper`, `System`, `ProbeSystemStatic`) implement `previous_potential_metadata_index`/`previous_gradient_metadata_index`, `target_uses_previous_influence_metadata` returns `false` by default for all of them, and every `fmm!` call in this pipeline (which requests `PowerRelativeGradient{0.001,0.001,true}` — a `RelativeErrorMethod`) silently degrades to absolute-tolerance-only accuracy control.

8. **Confirmed by direct bisection.** With the *original*, unmodified VortexLattice and FLOWVPM source (no workarounds applied):
   - `FastMultipole` @ `c8cf796` (last commit before `adc4f26`): `run_case()` completes **all 36 steps, no error**.
   - `FastMultipole` @ `adc4f26` itself: FLOWVPM's original code (with `get_previous_influence` defined) **fails to compile** (`UndefVarError`), confirming this is exactly the breaking commit.
   - `FastMultipole` @ current `main` (`58cf693`, with FLOWVPM's necessary local adaptation removing the dead override): same `NaN`/bracketing crash as originally reported, by step 3–4.

## Causal chain (best understanding)

```
FastMultipole adc4f26 removes get_previous_influence
        │
        ▼
FLOWVPM's ParticleField (and VortexLattice's probe/wake types) no longer
opt in to relative-error accuracy tracking — not yet ported to the new
previous_*_metadata_index API
        │
        ▼
Every fmm!() call requesting PowerRelativeGradient (a RelativeErrorMethod)
silently falls back to absolute-tolerance-only accuracy control
        │
        ▼
This changes FastMultipole's autotuned accuracy/interaction-list behavior
enough that wake-ring-corner probes (which sit essentially AT a filament's
own endpoint, by construction of PanelParticleWake's WakeBufferRings) are
now evaluated in a regime where VortexLattice's bound_induced_velocity
near-endpoint desingularization (get_δ) — which has always understated the
regularization (approaches 1/core_size² instead of saturating at
1/core_size) — actually gets exercised
        │
        ▼
Resulting O(10^3–10^4) spurious induced velocity at a wake corner →
one wildly displaced panel corner → NaN circulation/core-size in the
VPM particle field within 2-4 timesteps → Roots.find_zero bracketing
failure in FLOWVPM.solve_ρ_over_σ
```

We have **not** pinned down the exact internal FastMultipole mechanism by which the absolute-tolerance fallback changes which source/target pairs get routed to direct (near-field, exact) evaluation vs. multipole-approximated — that would need FastMultipole-side investigation/instrumentation. The causal link between the API change and the observed failure is established empirically (bisection above), not yet mechanistically at the FastMultipole-internals level.

## What we changed (workaround only, not a real fix)

In `VortexLattice/src/panel_particle_wake.jl`, swapped 6 call sites from `bound_induced_velocity` to `bound_induced_velocity_legacy` (properly `sqrt(r²+ε²)`-regularized). This reduces peak spurious velocities ~40× and delays the crash (step 3→step 6) but does **not** eliminate it — we believe a residual instance of the same `bound_induced_velocity` defect, reached via `src/fmm.jl`'s separate copies of the same call (used by the much more heavily-relied-upon `System`-as-FMM-source path), is still contributing. We attempted the same swap there; it caused the run to *hang* rather than error, so we reverted it — that swap is not safe to apply blindly and needs more careful attention (`bound_induced_velocity_legacy` may not robustly handle degenerate/zero-length segment cases that occur in that broader code path, which was presumably why `bound_induced_velocity` was introduced as a replacement in the first place, per commit `6c35823`).

## Suggested next steps for the developer

1. **FastMultipole**: either (a) restore a default/fallback path that doesn't silently degrade accuracy so drastically when `previous_*_metadata_index` isn't implemented, or at minimum (b) make the accuracy-control fallback behavior change more conservative/documented, since it's currently a silent behavior change for every downstream package that hasn't migrated to the new API yet (as of this report, that includes both FLOWVPM and VortexLattice).
2. **FLOWVPM / VortexLattice**: port to the new `previous_potential_metadata_index` / `previous_gradient_metadata_index` API to restore relative-error-based tuning. (VortexLattice's `System`/probe types and `WakeBufferRings`/`BoundaryFilamentWrapper`/`PanelBufferFilaments` would also need this.)
3. **VortexLattice, independent of the above**: `get_δ` / `bound_induced_velocity` (`src/induced.jl:26-67`) should be fixed regardless — its near-endpoint regularization is not physically bounded the way `bound_induced_velocity_legacy`'s is, and this is a latent correctness issue any time an evaluation point sits within `core_size` of a filament endpoint (which `PanelParticleWake`'s wake-ring-corner probes do by construction). This should be validated against the full existing test suite before rolling out broadly, since `bound_induced_velocity` is also used in `src/fmm.jl`'s widely-used `System`-as-source kernel.

## Immediate unblock (no code changes needed)

Pin `FastMultipole` to `c8cf796` or any earlier commit until (1)/(2) above land upstream. Confirmed working via bisection above.

---

## Addendum: fixes implemented, why the workaround alone wasn't right, and what's still open

Since the original report we (a) implemented the actual FastMultipole-side accommodation in FLOWVPM, (b) properly fixed the `bound_induced_velocity` defect in VortexLattice instead of just swapping to the legacy implementation, (c) discovered *why* the naive legacy-swap is unsafe if applied more broadly, and (d) ran the full existing test suite to check for regressions. Numbers below use `maxσ`/`maxΓ` (largest particle core size / circulation magnitude across the 36-step reproduction case) as a proxy for "how contaminated is the result" — `NaN` is the extreme end of that spectrum, but a run can complete without `NaN` and still be numerically polluted.

**Clean baseline** (`FastMultipole@c8cf796`, unmodified VortexLattice/FLOWVPM): `maxσ = 0.40`, `maxΓ = 0.28`.

### A. The real fix (suggestion #2 in the original report): FLOWVPM ported to the new metadata API

FastMultipole's own test suite (`test/vortex.jl`, `test/gravitational.jl`) contains the reference pattern for the new interface — it's structurally a drop-in replacement for what `get_previous_influence` used to do:

```julia
FastMultipole.metadata_per_body(system::VortexParticles) = 2
FastMultipole.previous_potential_metadata_index(system::VortexParticles) = 1
FastMultipole.previous_gradient_metadata_index(system::VortexParticles) = 2

function FastMultipole.metadata_to_buffer!(buffer, switch, i_buffer, system::VortexParticles, i_body)
    previous_potential = zero(eltype(system))
    previous_gradient = norm(SVector{3}(system.gradient_stretching[1, i_body], ...))
    buffer[FastMultipole.metadata_index(switch, 1), i_buffer] = previous_potential
    buffer[FastMultipole.metadata_index(switch, 2), i_buffer] = previous_gradient
end
```

We ported FLOWVPM's `ParticleField` to this exactly (`FLOWVPM/src/FLOWVPM_fmm.jl`), using `get_U(system, i)` as the "previous gradient" estimate — the same value the deleted `get_previous_influence` used. **With only this change** (VortexLattice's `bound_induced_velocity` left in its original, still-defective form), the 36-step reproduction case runs to completion with no `NaN`/crash — but `maxσ = 89.08`, `maxΓ = 73.9`, i.e. still heavily contaminated relative to baseline. So restoring the intended FastMultipole accuracy behavior is necessary but not sufficient on its own; VortexLattice's own formula defect still needs fixing.

Caveat: `FastMultipole.ProbeSystemStatic` (FastMultipole's *own* built-in probe target type, used by VortexLattice for control-point/wake-node evaluation) still does **not** implement the new metadata interface — we found no implementation for it anywhere in FastMultipole's source. The warning for it persists even after our FLOWVPM fix. Empirically this doesn't seem to block the improvement above, but it's an open gap worth FastMultipole's attention: currently *no* shipped built-in target type demonstrates the new API end-to-end apart from the test-only examples.

### B. Properly fixing `bound_induced_velocity` instead of swapping to `_legacy`

The original report's suggestion #3. We isolated the defect precisely to the `f2`/`f3` terms (the near-*endpoint* regularization) — confirmed both by direct computation (calling the function with the exact `target/v1/v4/gamma/cs` captured from the live crash reproduces the exact crash value, `-3169.8`) and by a clean synthetic sweep independent of the simulation:

```
r1 = (d, 0, 0), r2 = (d, 1, 0), core_size = 0.001
d:        0.1     0.01    0.001   1e-4    1e-5     1e-6
current:  0.79    7.96    79.6    783.0   6600.5   19934.2   (unbounded — diverges as d→0)
legacy:   0.79    7.88    39.8    7.88    0.796    0.0796    (bounded — decays toward 0)
```

Fix applied in `src/induced.jl`: replace `f2 = 1/(nr1+δ1)`, `f3 = 1/(nr2+δ2)` (where `δ1,δ2 = get_δ(nr1,core_size), get_δ(nr2,core_size)` — approaches `1/core_size²` as `r→0`, i.e. barely regularizes at all) with the standard `f2 = 1/sqrt(nr1²+core_size²)`, `f3 = 1/sqrt(nr2²+core_size²)` (saturates properly at `1/core_size`). `f1` (the numerator/cross-product term) is left untouched. Verified:

- **Exact algebraic equivalence to `bound_induced_velocity_legacy` when `finite_core=false`**: 20 random `(r1, r2)` pairs, bit-identical results. Confirms these are genuinely the same underlying Biot–Savart formula in different algebraic packaging — only the finite-core regularization differs.
- **No regression on the existing test suite** (see section D).

### C. Why we did *not* just swap everything to `bound_induced_velocity_legacy` (and why the earlier attempt to do so in `src/fmm.jl` hung)

We went back to understand the hang reported in the original document (swapping `src/fmm.jl`'s copies of `bound_induced_velocity` to `_legacy` caused the run to hang rather than error). Root cause: `bound_induced_velocity_legacy`'s `f1` denominator is

```
r1s*r2s - rdot² + εs*(r1s + r2s - 2·nr1·nr2)
```

By the Lagrange identity, `r1s*r2s - rdot² = |cross(r1,r2)|²`, and `r1s + r2s - 2·nr1·nr2 = (nr1-nr2)²`. **Both terms are exactly zero whenever a filament segment is degenerate** (`v1 == v2`, i.e. `r1 == r2` for any evaluation point) — a real, unremarkable occurrence in `src/fmm.jl`'s broader `System`-as-source usage (collapsed/zero-chord panel edges, symmetric-plane duplicate points, etc.). Whether this lands on **exact** `0.0` (→ numerator also exactly `0` → `0/0 = NaN`) or a tiny nonzero floating-point rounding residual (→ safely `0`) is **non-deterministic**, dependent on instruction-level rounding order, not on `core_size` or any code path we can control. We reproduced both outcomes directly: an interactively-evaluated expression landed on exact `0.0` (`NaN`), while the compiled function itself, given the identical mathematical inputs, landed on `3.47e-18` (safe `0`). We believe an unlucky exact-cancellation case is what poisoned FastMultipole's tree/error-estimation machinery and produced the hang.

By contrast, the current (non-legacy) `bound_induced_velocity`'s `f1` denominator, `nr1·nr2 + dot(r1,r2)`, is `2·nr1²` for a degenerate segment — safely bounded away from zero whenever `nr1 > 0`, with no delicate cancellation involved. This is presumably *why* `bound_induced_velocity` was introduced in the first place (commit `6c35823`, "shed unsteady particles") — likely precisely to handle degenerate segments arising during particle shedding, at the cost of (unintentionally) weaker endpoint regularization. **This is why our fix only touches `f2`/`f3` and deliberately leaves `f1`'s structure alone.**

### D. A residual gap we chose not to chase further, and why

With (A) + (B) above (i.e. no `panel_particle_wake.jl`-specific workaround at all): `maxσ = 24.39`, `maxΓ = 15.9` — a real improvement over (A) alone (89.08/73.9), but still ~60× the clean baseline (0.40/0.28). Re-adding the `bound_induced_velocity_legacy` swap in `panel_particle_wake.jl`'s 3 wake-ring-specific call sites on top of (B) gets to `maxσ = 2.42`, `maxΓ = 1.03` — better, but still ~6× baseline.

We dug into why `_legacy` still outperforms our fixed `bound_induced_velocity` in practice, since after fix (B) their `f2`/`f3` terms are now identical:

- Random-sampled comparison (20 pairs, `finite_core=false`): 0 relative difference (proves the two formulas are the same physical quantity, as expected).
- Same comparison with `finite_core=true`, `core_size` small relative to typical separations: `f1` values diverge substantially depending on geometry (e.g. one synthetic case: `f1_new ≈ 0.925` vs. `f1_legacy ≈ 99.0`, a ~107× difference) even after our `f2`/`f3` fix. Root cause: `legacy`'s `f1` denominator is a *quadratic* form in `(nr1, nr2, dot)` (`r1s*r2s - rdot²`, i.e. `|cross|²`, plus a `core_size`-scaled quadratic correction), while the current formula's `f1` denominator is a *linear* form (`nr1·nr2 + dot`). These are algebraically equal at `core_size=0` (proven above) but respond very differently once a `core_size`-dependent regularizer is layered on top — the quadratic form damps near-collinear/near-antiparallel geometry (target near the segment's extension line) much more aggressively.
- We checked whether this was itself introduced by the November 2025 "bugfix" commit (`752e8e1`) by reverting just the `f1` regularizer to its pre-commit form (`max` of three geometric-distance-based `get_δ` terms, rather than a single value-based `get_δ(denom, core_size)` term). Numerically this pre-commit variant is **nearly identical** to the current one (differences in the 5th–6th significant digit) — so that commit's `f1` change was essentially a wash, not the source of the gap. The real gap is the older, structural `linear` vs. `quadratic` difference between the "new" and "legacy" formulas, which predates that commit.

Closing this gap properly would mean giving `f1` the quadratic form's stronger damping *without* reintroducing the degenerate-segment fragility from section C — i.e. a real (if modest) piece of numerical-methods design, not a mechanical patch, and something we did not want to rush without dedicated validation against known analytical cases. We're flagging it as open rather than guessing further. **In the interim we kept `bound_induced_velocity_legacy` specifically for the 3 `panel_particle_wake.jl` wake-ring call sites** (empirically the best-performing option there, and — unlike `src/fmm.jl`'s usage — not exposed to degenerate/zero-length wake-ring segments in practice), while every other caller (including `src/fmm.jl`) now benefits automatically from the `f2`/`f3` fix in the shared `bound_induced_velocity` function.

### D2. Prior art: FLOWVPM's `origin/flowpanel` branch already has (most of) fix (A)

Before finalizing, we checked whether this had already been solved elsewhere and hadn't yet been merged. It had: `FLOWVPM`'s `origin/flowpanel` branch (10 commits ahead of `master`, last commit `a950790`, 2026-06-30 — *after* the FastMultipole `adc4f26` break) independently implements essentially the same `metadata_per_body`/`previous_potential_metadata_index`/`previous_gradient_metadata_index`/`metadata_to_buffer!` port for `ParticleField` that we wrote from scratch, plus additional migration to FastMultipole's newer switch-aware `set_gradient!`/`set_hessian!` signatures (`DerivativesSwitch{PS,GS,HS,NO,NM}`) and vorticity-as-extra-output support that our minimal patch doesn't attempt. It targets `FastMultipole = "2.2.0"`; our checked-out `FastMultipole@main` is already at `2.3.0`, so it's compatible.

We checked out `flowpanel` in place of our hand-written patch and reran the reproduction case (VortexLattice's `induced.jl`/`panel_particle_wake.jl` fixes unchanged): full 36-step completion, `maxσ = 2.58`, `maxΓ = 1.85` — same ballpark as our patch (`2.42`/`1.03`), confirming both implementations are correct and cross-validating each other (the residual gap vs. baseline is unchanged either way, since it's the separate `f1`/section-D issue in VortexLattice, orthogonal to which FLOWVPM implementation is used). Note: switching FLOWVPM's dev branch requires re-running `Pkg.resolve()` in VortexLattice's environment, since `flowpanel` adds an explicit `LinearAlgebra` dependency that `master` doesn't declare.

We did not adopt `flowpanel` for the ongoing fix, since it isn't a clean, isolated cherry-pick (it carries unrelated commits — "relaxation filter", "working with FLOWPanel" — and switching a `Pkg.develop`'d path dependency's branch affects every project using that FLOWVPM checkout, not just this one). We kept our smaller, scoped patch instead, easy to drop once the FastMultipole-compatibility portion of `flowpanel` (or equivalent) is merged into `master` by its author. No corresponding work exists on the VortexLattice side on any branch (checked all local and remote branches for the new metadata API — zero hits) — `WakeBufferRings`/`BoundaryFilamentWrapper`/`PanelBufferFilaments`/`System` remain unported there regardless of which FLOWVPM branch is used, consistent with the `FastMultipole.ProbeSystemStatic` warning persisting in both cases.

### D3. Deriving what `f1` should actually be (both singularities are real; neither existing formula gets both right; a derived candidate fix hangs in practice and needs further hardening)

This follows up on section D's open question with an actual derivation, prompted by asking "can we just use a bigger `core_size`?" — the answer is no, and here's the reasoning.

**Setting up the exact (unregularized) formula.** With `r1 = P - A`, `r2 = P - B` (field point to each segment endpoint) and `r0 = r1 - r2` (`= B - A` up to sign), the standard closed-form Biot–Savart velocity for a straight vortex segment (e.g. Katz & Plotkin) is

```
V = (Γ/4π) · [cross(r1,r2) / |cross(r1,r2)|²] · [dot(r0,r1)/|r1| − dot(r0,r2)/|r2|]
```

We verified numerically that `bound_induced_velocity_legacy` at `core_size=0` matches this exactly, and that `bound_induced_velocity` (current) matches `bound_induced_velocity_legacy` exactly at `core_size=0` too (20 random trials, exact to machine precision) — both formulas are the *same* physical quantity in different algebraic packaging; only their `core_size>0` regularizations diverge.

**There are two distinct near-field regimes, and they need to be told apart.** Using the geometric identity `|cross(r1,r2)| = |r0| · r_c` (segment length times perpendicular distance from the field point to the *infinite* line through the segment — literally twice the area of the triangle formed by the field point and the two endpoints):

1. **Endpoint terms are not actually singular as `|r1| → 0`.** Algebraically, `dot(r0,r1)/|r1| = (r1s − rdot)/|r1| = |r1| − |r2|cos θ → −|r2|cos θ`, a *finite* limit depending only on approach direction. There is no true mathematical singularity here — only a numerically ill-conditioned way of computing a well-behaved quantity (a small-`|r1|` cancellation, not a divergence).
2. **The real singularity lives entirely in `f1`, via `r_c`.** We confirmed this with `BigFloat` (eliminating floating-point cancellation as an explanation): approaching a segment's endpoint from a direction perpendicular to it produces a genuine, unbounded `1/r_c` divergence, growing cleanly as `1/d` from `d=0.1` all the way to `d=1e-40`. This is the ordinary near-filament singularity of an idealized (zero-thickness) vortex line, and it's exactly what `core_size` is supposed to regularize — via the standard viscous-core substitution `r_c → sqrt(r_c² + core_size²)`.

**Neither existing formula regularizes `r_c` directly — each stumbles onto approximately correct behavior in only one sub-case of it.** We tested both formulas' *saturation scaling* rigorously: fixing `d = core_size/1000` (always deep inside the core) and sweeping `core_size` over 4 decades, the correct physical answer is `value · core_size → constant`.

- **Endpoint approach, perpendicular to the filament** (`r_c → 0` because the field point is near one specific endpoint, off to the side): `bound_induced_velocity_legacy`'s `value·core_size` converges cleanly to `≈7.96e-5` (`= [1/(4π)] × 10⁻³`, i.e. exactly the expected `Γ/(4π·core_size)` semi-infinite-line-like saturation). `bound_induced_velocity` (current) does **not** converge — `value·core_size` keeps *growing* (`1.4e-4 → 8.7e-4 → 7.3e-3 → 0.040 → 0.072` across the same 4 decades) — it never saturates.
- **Midline approach** (`r_c → 0` because the field point is near the segment's *middle*, roughly equidistant from both endpoints — not near either corner): `bound_induced_velocity` converges cleanly (`value·core_size ≈ 3.18e-4`, constant). `bound_induced_velocity_legacy` does **not** converge — its raw (unmultiplied) value tracks the fully *unregularized* `Γ/(2π·d)` law almost exactly (`159155.9` vs. the predicted `159154.9` at `d=1e-6`, `core_size=0.001`) — its core term simply never engages in this geometry, regardless of how small `core_size` is set.

So: each formula correctly regularizes one of the two regimes and is silently, unboundedly broken in the other. **No choice of `core_size` fixes a formula whose saturation scaling law is wrong** — it only moves *where* the wrong behavior starts, not *whether* it happens.

**A candidate fix, derived from the `r_c` identity above, and its validation.** Replacing only `f1`'s singular part with the direct, physically-motivated regularization —

```julia
f1 = cross(r1,r2) / (nr3² · r_c · sqrt(r_c² + core_size²))     # nr3 = |r1 - r2|, r_c = |cross(r1,r2)| / nr3
```

— paired with `bound_induced_velocity_legacy`'s own (self-consistent, and per the point above, not independently in need of regularization) `f2`, with a fallback to the current formula's degenerate-segment-safe form when `nr3 ≈ 0`:

- Matches the reference formula exactly at `core_size=0` (`6e-16`, machine epsilon, 50 random trials).
- Converges cleanly in **both** regimes: `value·core_size → 1/(4π) ≈ 0.0796` for the corner approach, `→ 1/(2π) ≈ 0.159` for the midline approach — and that clean factor-of-2 relationship is exactly what physical intuition predicts (an interior point of the filament "sees" contributions from both directions along the line; a corner only sees one, i.e. half).
- Stays safe (no `NaN`) for degenerate zero-length segments.
- On the real captured crash geometry, gives `-93.3` (vs. current's `-88.3` and legacy's `-2.5`) — much closer to `current` than to `legacy`, suggesting `legacy`'s small answer there was itself an *under*-correction for that particular mixed-angle geometry, not a better one.

**But: swapping this candidate into the actual reproduction case hangs.** We monkey-patched `bound_induced_velocity` with this candidate and re-ran the full 36-step `PanelParticleWake` case — it hung (5+ minutes of CPU with zero output, killed manually), the same failure mode as the earlier `bound_induced_velocity_legacy`-in-`fmm.jl` hang from section C, despite this candidate having passed every synthetic test we threw at it (including an explicit degenerate-segment check). We have not root-caused this specific hang; the likely culprit is some other near-degenerate configuration (e.g. `nr3sq` just above our `1e-28` cutoff while genuinely nearly-degenerate, amplifying floating-point noise in `r_c` unpredictably) that our synthetic tests didn't happen to construct. **This is the headline takeaway of this whole appendix: even a formula that is provably, exactly correct in the unregularized limit and demonstrably correct in its regularized asymptotic scaling in both known singular regimes can still have its own undiscovered numerical fragility once exercised against the full geometric diversity of a real simulation.** Getting this right needs dedicated hardening and testing (systematically covering near-degenerate-but-not-quite segments, mixed-regime approaches, etc.) before it's safe to ship — not something to rush in as a follow-on to this investigation.

### E. Full test suite run

Ran `test/runtests.jl` (12 AVL regression cases validated against known-correct AVL output, plus geometry/wake/unsteady/`PanelParticleWake` tests) both with and without the `induced.jl` fix. Three testsets fail — **all three confirmed pre-existing** by reverting `induced.jl` entirely and reproducing bit-identical failures:

1. **"Unsteady Vortex Lattice Method — Rectangular Wing"** (`test/runtests.jl:1272`): `BoundsError: attempt to access 50×13 Matrix{WakePanel{Float64}} at index [3, 0]`, inside `induced_velocity(::CartesianIndex, ...)` (`src/induced.jl:1102`) via the *old-style*, non-`PanelParticleWake` `unsteady_analysis!` path. Unrelated to anything touched this session.
2. **"FluidDomainMonitor"** (`test/runtests.jl:1531`): `UndefVarError: _build_short_restart_case not defined`. This helper function is defined *inside* the sibling `"PanelParticleWake Restart (Short)"` `@testset` block (line 1451) and is not visible to `"FluidDomainMonitor"` — `@testset begin...end` introduces its own scope in Julia's `Test` stdlib (confirmed with a 6-line minimal repro), so this was never going to work regardless of any of our changes. It likely went unnoticed because the suite always aborted earlier (at failure #1) before ever reaching it.
3. **"FMM induced velocity"** (`test/fmm_test.jl`): the manual-filament-loop assertion (which directly exercises our fixed `bound_induced_velocity`) **passes**; the two assertions comparing against `FastMultipole.direct!`/`fmm!` output fail with what looks like a component-shift/buffer-layout mismatch (e.g. `v_direct = [0.0, -0.00192, 0.00059]` where `v_VL = [-0.00192, 0.00059, 0.00727]` — the last two components of `v_VL` appear shifted into `v_direct`'s last two slots, with the first replaced by `0.0`). Reproduces identically with `induced.jl` reverted.

### F. Fixed all three pre-existing test failures

**F1. "FluidDomainMonitor" scoping bug.** Hoisted `_build_short_restart_case` out of the `"PanelParticleWake Restart (Short)"` `@testset` to file scope (`test/runtests.jl`), so both it and `"FluidDomainMonitor"` can call it. Both testsets pass afterward (4/4 and 3/3).

**F2. "Unsteady Vortex Lattice Method — Rectangular Wing" `BoundsError`.** Root cause: `get_wake_velocities!` (`src/wake.jl:583-599`) computes a same-surface wake self-induced-velocity contribution using a deliberately-reduced row count, `nc_wake = max(nwake[jsurf]-1, 0)`, to exclude the newest/still-attached wake row from self-influence — but the outer loop's row index `I[1]` still ranges up to `nwake[jsurf]+1`, one past what `nc_wake` accounts for. When `I[1] == nwake[jsurf]+1`, `induced_velocity(::CartesianIndex, ...)`'s branch logic (which assumes `I[1] <= nc_wake+1`) falls through to its final `else` clause under a false assumption, computing a bogus `wake[nc_wake, -1]`-style index — the observed `wake[3, 0]`. Fixed by skipping this contribution when `I[1] > nc_wake + 1`, with a comment explaining why. Both `"...Rectangular Wing"` and the same-shaped `"...Wing + Tail"` testset (which hits the identical code path but was never reached before because the suite aborted at the first failure) now run to completion.

**F3. "FMM induced velocity" component-shift mismatch — same root cause found in three other places, likely relevant to the report's still-open residual gap.** This one connects back to the FastMultipole `adc4f26` API change discussed throughout this report. `FastMultipole.get_gradient`/`set_gradient!` now have both a legacy 2-/3-arg form (hardcoded to buffer rows `5:7`, kept only for backward compatibility) and a **switch-aware** form (`..., switch::DerivativesSwitch, ...`) that computes the actual row range from `gradient_range(switch)` — which depends on `PS` (whether scalar potential output is enabled) and `NM` (number of metadata rows). VortexLattice's own `direct!` implementations write with the *legacy* 2-arg setter (rows hardcoded to `5:7`) while the tree's *reader* uses the switch-aware getter — these only agree by coincidence when `PS=true, NM=0`. `test/fmm_test.jl`'s `fmm!` call uses `PS=false`, so the switch-aware gradient range starts at row `4` instead of `5` — a one-row read/write mismatch that reproduces the observed symptom exactly: `v_direct = [0.0, v_VL[1], v_VL[2]]` (row 4 read as component 1, still zero/never written; rows 5,6 — actually holding components 1,2 — read as components 2,3; component 3, written to row 7, is never read at all).

  We found **the same bare 2-arg `set_gradient!` anti-pattern in five places** — `src/fmm.jl`'s `direct!` for `System` and `FilamentWrapper` (2 sites), and `src/panel_particle_wake.jl`'s `direct!` for `PanelBufferFilaments`, `WakeBufferRings`, and `BoundaryFilamentWrapper` (3 sites) — all silently discarding the `DerivativesSwitch` argument instead of naming and forwarding it. Fixed all five by naming the switch parameter and passing it through to the switch-aware setter.

  Effect on `fmm_test.jl`: the `FastMultipole.direct!` comparison (line 92) now **passes exactly**; the `fmm!` (multipole-approximated) comparison (line 99) improves from being wrong by orders of magnitude to matching to ~7 significant figures, but still narrowly misses the test's very tight `atol=1e-12` (actual absolute difference ≈ `1.4e-11`). We checked whether this last piece connects to the `ProbeSystemStatic`/relative-error story from earlier in this report — it doesn't: this specific `fmm!` call never passes an `error_tolerance` at all, so the relative-error/metadata fallback path isn't engaged for it. The residual looks like an ordinary multipole-truncation/floating-point artifact at `expansion_order=20`, not a bug, and we did not chase it further.

  Effect on the `PanelParticleWake` reproduction case (relevant since 3 of the 5 sites are in the hybrid-wake code): re-ran the full 36-step case with these fixes added on top of everything else in this report. Result was mixed rather than a clean further improvement — `maxΓ` essentially unchanged (`0.99` vs. `1.03`), `maxσ` moved *up* slightly (`4.02` vs. `2.42`) and `minσ` dropped well below the clean baseline's `0.22` floor (to `0.045`). This isn't a regression from these changes — it's correct information (gradients/hessians) now reaching parts of the code that were previously silently reading stale/zero/misaligned buffer rows, so the simulation's dynamics shift to reflect physics that's more correct than before, not necessarily "closer to baseline" by this rough proxy metric. We consider this a net-positive, clearly-justified correctness fix independent of its effect on that specific number; a rigorous re-assessment of the residual gap from section D would need to be redone against this corrected baseline rather than trusted from before.

### G. Reproduction case status with every fix in this report applied

Re-ran the full 36-step reproduction case one more time with everything above in place (FLOWVPM metadata port, `induced.jl` `f2`/`f3` fix, `panel_particle_wake.jl` `bound_induced_velocity_legacy` retention + switch-aware `set_gradient!` fix, `wake.jl` off-by-one fix). Completes cleanly, no crash, no `NaN`, deterministic and bit-identical across repeated runs:

| | `maxσ` | `minσ` | `maxΓ` | `mean|Γ|` |
|---|---|---|---|---|
| Clean baseline (`FastMultipole@c8cf796`, unmodified everything) | 0.40 | ~0.22 (typical particle floor) | 0.28 | — |
| All fixes in this report applied | 4.02 | 0.045 | 0.99 | 0.080 |

As noted in section F3, this is not a monotonic "closer to baseline" result, and shouldn't be read as one — the `set_gradient!` fix corrected a distinct, previously-silent buffer-misalignment bug, so the comparison basis shifted. The still-open item is section D's `f1` structural gap (`bound_induced_velocity`'s linear-form denominator vs. `bound_induced_velocity_legacy`'s quadratic-form one), which nothing in section F touches.

### Updated file/change summary (supersedes the one above)

- `FLOWVPM/src/FLOWVPM_fmm.jl`: `ParticleField` ported to the new `previous_potential_metadata_index`/`previous_gradient_metadata_index`/`metadata_to_buffer!` API (real fix for the root cause; cross-validated against FLOWVPM's own `origin/flowpanel` branch, which independently implements the same thing).
- `VortexLattice/src/induced.jl`: `bound_induced_velocity`'s `f2`/`f3` corrected to a properly-bounded `sqrt(r²+core_size²)` regularization; `f1` deliberately left untouched (see section C).
- `VortexLattice/src/panel_particle_wake.jl`: `bound_induced_velocity_legacy` kept for the 3 wake-ring-specific call sites as a targeted, empirically-justified supplement; all 5 `direct!` implementations across this file and `src/fmm.jl` fixed to use switch-aware `set_gradient!` instead of the legacy 2-arg form (section F3).
- `VortexLattice/src/wake.jl`: fixed an off-by-one in `get_wake_velocities!`'s same-surface wake self-influence bounds check (section F2).
- `VortexLattice/test/runtests.jl`: hoisted a test helper function out of testset-local scope so a sibling testset can use it (section F1).
- No FastMultipole changes.

Final full-suite status: every testset passes except one assertion in `"FMM induced velocity"` (section F3's tight-tolerance residual, believed unrelated to this report's core issue).

### H. Appendix: self-contained code for the candidate `f1` fix and the hang-debugging harness (inlined so nothing depends on ephemeral scratchpad files)

**H1. The candidate fix itself (`biv_candidate2`), exactly as tested.** This has NOT been applied to `src/induced.jl` — it exists only as a runtime monkey-patch in throwaway scripts. `get_δ` here is `VortexLattice.get_δ` (already in `src/induced.jl`, unchanged).

```julia
using StaticArrays, LinearAlgebra
get_δ(distance, core_size) = distance < core_size ? (distance-core_size)^2 : zero(distance)

# candidate2: rc-based (perpendicular-distance-to-line) regularization for f1, paired
# with legacy's own (self-consistent) f2 term -- f2 needs no separate core
# regularization since (r1s-rdot)/nr1 = nr1 - nr2*cosθ has a finite limit as nr1->0
# (verified algebraically); keep the sqrt(r1s+εs) form anyway as a numerical-conditioning
# safety net at the exact endpoint (nr1=0 identically), consistent with legacy.
function biv_candidate2(r1, r2, finite_core, cs)
    nr1, nr2 = norm(r1), norm(r2)
    rdot = dot(r1, r2)
    r1s, r2s, εs = nr1^2, nr2^2, cs^2
    num = cross(r1, r2)
    if !finite_core
        f1 = num/(nr1*nr2 + rdot)
        return (f1*(1/nr1+1/nr2))/(4pi)
    end
    r12 = r1 - r2
    nr3sq = dot(r12, r12)
    if nr3sq > 1e-28
        crossnorm = norm(num)
        rc = crossnorm / sqrt(nr3sq)
        f1 = num / (nr3sq * rc * sqrt(rc^2 + εs))
    else
        # degenerate/zero-length segment fallback: same denominator structure as
        # the current (shipped) bound_induced_velocity's f1, safe for r1==r2.
        f1 = num / (2*nr1*nr1 + get_δ(nr1*nr2+rdot, cs))
    end
    f2 = (r1s - rdot)/sqrt(r1s + εs) + (r2s - rdot)/sqrt(r2s + εs)
    return (f1*f2)/(4pi)
end
```

To wire this into VortexLattice at runtime for testing (do NOT edit `src/induced.jl` directly until the hang is understood):

```julia
using VortexLattice
@eval VortexLattice function bound_induced_velocity(r1, r2, finite_core, core_size)
    return $biv_candidate2(r1, r2, finite_core, core_size)
end
```

**H2. Validation checks that passed** (all against the reference unregularized Biot-Savart formula and against `bound_induced_velocity_legacy`):
- 50 random `(r1,r2)` trials at `core_size=0`: matches `biv_legacy` to `6e-16` (machine epsilon).
- Saturation-scaling sweep (`d = core_size/1000`, 4 decades of `core_size`): `value·core_size → 1/(4π)` for corner-perpendicular approach, `→ 1/(2π)` for midline approach — both converge cleanly, unlike either existing formula alone (see section D3 body above for full numbers).
- Exact degenerate segment (`r1 == r2`): no `NaN`, falls through to the safe branch.
- Real captured crash-case geometry (`target/v1/v4/gamma/cs` from the live original crash, see section 3 of the Investigation timeline above) reproduces a finite, sane-order-of-magnitude value.

**H3. The reproduction harness used to find the hang.** `test/profile_panel_particle_wake_impl.jl` (untracked, already in this repo) exposes `build_case`/`run_case` without an auto-run at file scope (unlike its sibling `test/profile_panel_particle_wake.jl`, which does auto-run — useful for quick interactive checks but not for scripting multiple calls in one process). Minimal repro of the hang:

```julia
using StaticArrays, LinearAlgebra, VortexLattice, FLOWVPM
# ... define get_δ and biv_candidate2 as in H1, then monkey-patch as shown ...
include("test/profile_panel_particle_wake_impl.jl")
run_case(n_steps=1)   # hangs -- does not return, spins at 100% CPU indefinitely
```

**H4. Instrumentation that localized the hang (call-tracing).** Wrap `biv_candidate2` with a call counter and periodic checkpoint print before monkey-patching:

```julia
const NCALLS = Ref(0)
function biv_candidate2_traced(r1, r2, finite_core, cs)
    NCALLS[] += 1
    if NCALLS[] % 500 == 0
        println(stderr, "call #$(NCALLS[])  r1=$r1 r2=$r2 fc=$finite_core cs=$cs")
        flush(stderr)
    end
    return biv_candidate2(r1, r2, finite_core, cs)
end
# monkey-patch bound_induced_velocity to call biv_candidate2_traced instead
```

Result (reproduced twice independently): the call counter advances steadily to ~13000-13500 total calls within `n_steps=1` (a single timestep), then **stops advancing entirely** while the OS-level process (confirm via `ps aux | grep julia` — look for `%CPU` near 100 and a large, roughly-stable RSS) continues consuming 100% CPU for 5+ more minutes with zero further output. This means:
- The hang is not inside `bound_induced_velocity`/`biv_candidate2` itself (it would keep incrementing `NCALLS` if it were looping there).
- It happens on the very first timestep, far earlier than previously assumed.
- Execution is stuck in some caller/downstream function — most likely FastMultipole's dual-tree traversal or adaptive error-based tree refinement, which could plausibly loop indefinitely if one pathological (`Inf`/`NaN`/anomalously large) induced-velocity value from `biv_candidate2` prevents a local error criterion from ever being satisfied. **Not yet confirmed** — would need direct instrumentation of FastMultipole's own traversal/refinement code (a separate package) to pin down further.

**H5. Dead-end technique to avoid repeating.** Attempting to instrument *FLOWVPM's* `upper_bound_abs`/`solve_ρ_over_σ` (to check for a pathological `σ`/`ω`/`ε` triple feeding a broken `Roots.find_zero` bracket — a plausible-sounding alternate hypothesis) via the same runtime `@eval` monkey-patch technique is **not usable**: redefining a function that's a hot, deeply-inlined dependency of FastMultipole's generic tree-traversal dispatch triggers a massive Julia method-invalidation/recompilation cascade. One such attempt burned 10+ minutes of 100% CPU with zero output before being killed manually — indistinguishable from a real hang by symptoms alone, and wasted significant time before being identified as a compile-time artifact rather than a runtime one. If this hypothesis is worth checking, either (a) add prints by editing the FLOWVPM package source file directly (accepting a one-time recompile cost, at least not repeated across multiple live-eval attempts), or (b) use a debugger/profiler attached to a run using the *unmodified* FLOWVPM, rather than redefining its internals live.

**H6. Suggested concrete next step for whoever picks this up:** rather than continuing to guess at FastMultipole's internal behavior, attach `Debugger.jl` or send `SIGINFO`/interrupt (`Ctrl-\` / `kill -QUIT`, then inspect the stacktrace Julia prints on `InterruptException` — note plain `Ctrl-C`/`SIGINT` may not interrupt a tight non-yielding numerical loop, may need `kill -TERM` and catch the resulting stacktrace, or run under `perf`/`sample` (macOS: `sample <pid> 5` while it's hung) to get a real call stack at the moment of the hang, rather than continuing to bisect by call-counting from the outside.

---

## I. Section D3's open `f1` item: RESOLVED (2026-07-27, later session)

**Status: fix applied to `src/induced.jl`, full suite run, 36-step reproduction case runs clean in ~17 s with no hang and no `NaN`.** This supersedes section D3's "open, actively being worked" item and the "current status" note at the top of this report.

### I1. First: is the `f1` defect real, or did we cause it?

Re-verified from scratch, deliberately with **nothing external in the loop** — three formula variants transcribed verbatim from git into a standalone script importing only `LinearAlgebra`/`StaticArrays` (no VortexLattice, no FastMultipole, no FLOWVPM), checked against a `BigFloat` Katz–Plotkin reference:

- `master` = `git show master:src/induced.jl` (the quadratic-denominator formula, i.e. what is now `bound_induced_velocity_legacy`)
- `pristine` = `git show HEAD:src/induced.jl` (frames branch, before this investigation's edits)
- `current` = working tree at the time (pristine + the section-B `f2`/`f3` fix)

**Provenance.** The linear-denominator `bound_induced_velocity` was introduced by commit `6c35823` ("shed unsteady particles") on the `frames` branch, which simultaneously renamed the pre-existing formula to `bound_induced_velocity_legacy` (`git log -S'bound_induced_velocity_legacy' -- src/induced.jl` → exactly one hit). `master` has only the legacy formula. Nothing in this investigation, and nothing in the FastMultipole/FLOWVPM updates, created it.

**It is a genuine defect, provable by pure arithmetic:**

1. All three variants agree with the exact Biot–Savart reference to `5e-16` at `core_size=0` — same physics, differing only in regularization.
2. **Scale invariance fails.** `|V|` per unit circulation has units `1/L`, so scaling the geometry *and* `core_size` by `s` must scale `|V|` by exactly `1/s`. Near-axis/near-endpoint probe, `s` from `1e-3` to `1e6`: `master` gives `6.934493` at every scale; `pristine`/`current` drift `1253.6 → 2356.1`. A benign far-field control point is scale-invariant for all three, isolating the failure to the `finite_core` branch of `f1`.
3. **Root cause is a units error**, visible by inspection once you look for it: `δ0 = get_δ(denom, core_size)` compares `denom = nr1*nr2 + dot(r1,r2)` (units `L²`) against `core_size` (units `L`). Whether the desingularizer engages at all therefore depends on the absolute size of the model.
4. **Worst-case overshoot** relative to the physical bound `1/(2π·core_size)`, max over 400k randomly sampled near-field points per core size (segment length 1):

   | `core_size` | master | pristine | current | fix (I3) |
   |---|---|---|---|---|
   | 1e-1 | 30266 | 1.88 | 0.94 | 0.48 |
   | 1e-2 | 51675 | 15.2 | 2.54 | 0.50 |
   | 1e-3 | 37678 | 121.7 | 7.87 | 0.50 |
   | 1e-4 | 5484 | 264.7 | 24.9 | 0.50 |
   | 1e-5 | 3281 | 696 | 77.0 | 0.50 |

   The section-B `f2`/`f3` fix helped substantially (122× → 7.9× at `cs=1e-3`) but the residual overshoot still grows without bound as `core_size` shrinks.

**Correction to section D3.** Its claim that the `f2`/`f3`-fixed formula "does not converge" in the corner regime is wrong — the sequence it reported (`1.4e-4 → … → 0.072`) is *converging* to `1/(4π) = 0.0796`, just slowly; extending the sweep to `core_size=1e-10` gives `0.0795774`. The corner regime was already fixed by section B. The actual remaining failure regime is **near-axis, near-endpoint** (`rc → 0` with `r1`, `r2` nearly anti-parallel, so `denom → 0`), which neither of D3's two probe geometries constructed.

### I2. Why `biv_candidate2` (section D3/H1) hung

It forms `rc = norm(cross(r1,r2))/sqrt(nr3sq)` explicitly and then divides by it. For any field point lying **exactly on the segment's infinite line** — which flat, unswept lattice geometry produces constantly — `cross(r1,r2)` is exactly the zero vector, so `f1 = num/(nr3sq * rc * sqrt(rc^2 + εs))` is `0/0`:

```
point on segment, at midline    [NaN, NaN, NaN]
point on line beyond endpoint   [NaN, NaN, NaN]
point exactly at endpoint A     [NaN, NaN, NaN]
```

`NaN` injected into FastMultipole's error-based refinement criterion accounts exactly for the observed symptom in H4 (call counter stops advancing, 100% CPU, no output, no progress). No FastMultipole-internals instrumentation was needed after all; H6's suggested next step is moot.

### I3. The applied fix

Regularize the one true singularity — `rc`, the perpendicular distance from the field point to the segment's infinite line — via the standard `rc → sqrt(rc² + core_size²)`, but **never form `rc`**, which is what killed `candidate2`:

```
norm(cross)^2 + core_size^2 * nr3^2  ==  nr3^2 * (rc^2 + core_size^2)
```

and rewrite `f1` with the identity `cross/(nr1*nr2 + rdot) == cross*(nr1*nr2 - rdot)/norm(cross)^2` so the regularized form still reduces exactly to the unregularized formula as `core_size → 0`:

```julia
εs = core_size * core_size
r12 = r1 - r2
nr3s = dot(r12, r12)                # squared segment length
f1_denom = dot(num, num) + εs*nr3s  # == nr3s * (rc^2 + core_size^2)
iszero(f1_denom) && return zero(num)
f1 = num*(nr1*nr2 - dot(r1, r2))/f1_denom
f2 = 1/sqrt(nr1*nr1 + εs)
f3 = 1/sqrt(nr2*nr2 + εs)
```

`f1_denom` vanishes only for a genuinely zero-length segment, where the induced velocity is zero anyway — so the guard is an exact, non-heuristic test, not `candidate2`'s absolute `nr3sq > 1e-28` cliff. `get_δ` is now referenced only in comments; left in place rather than deleted.

**Standalone validation** (same no-external-dependencies harness as I1):
- Exact vs. the `BigFloat` reference at `core_size=0`: max rel. err `1.5e-15` over 200 random points.
- Scale invariant to 7 significant figures over 9 decades of geometry scale.
- Bounded at exactly `1/(4π·core_size)` worst case for every core size tested, `0` `NaN`/`Inf` in 2M near-field samples (table in I1).
- Correct Rankine-core asymptotics in **both** regimes: `value·core_size → 1/(4π)` corner, `→ 1/(2π)` midline at the core radius, with linear decay to 0 inside the core.
- Finite (`0`, not `NaN`) on every degenerate input that broke `_legacy` and `candidate2`: zero-length segment, collapsed segment at the origin, field point exactly at an endpoint, exactly at the midline, on the line beyond an endpoint, tiny-but-nonzero segment off-line, and near-antiparallel vectors of magnitude `1e8`.

### I4. Test suite and reproduction case

`test/runtests.jl`, run **from the `test/` directory** (running it from the repo root aborts early in `"OpenVSP Geometry Import"` on a relative path to `samplewing.csv` — a harness quirk, unrelated to any of this): every testset passes except the single known tight-tolerance assertion at `test/fmm_test.jl:99` (section F3's `fmm!` multipole-truncation residual, `atol=1e-12`). Its absolute residual moved from `~1.4e-11` to `~5.3e-10` — expected, since `f1` changed and both sides of that comparison use the kernel; the `FastMultipole.direct!` comparison at line 92 still passes exactly.

36-step `PanelParticleWake` reproduction case — **completes in ~17 s, no hang, no `NaN`**:

| | `maxσ` | `minσ` | `maxΓ` | `mean|Γ|` |
|---|---|---|---|---|
| Clean baseline (`FastMultipole@c8cf796`, unmodified everything) | 0.40 | ~0.22 | 0.28 | — |
| Section G (all prior fixes, `f1` untouched) | 4.02 | 0.045 | 0.99 | 0.080 |
| **+ this `f1` fix** | **3.54** | 0.059 | **0.97** | 0.081 |
| + this `f1` fix, `_legacy` workaround also bypassed | **3.41** | 0.054 | **0.64** | 0.068 |

### I5. The `panel_particle_wake.jl` `_legacy` workaround has been removed

Last row above: routing the 3 wake-ring `bound_induced_velocity_legacy` call sites through the fixed `bound_induced_velocity` improves every metric — `maxΓ` `0.97 → 0.64`, `maxσ` `3.54 → 3.41`. That is expected: `_legacy` is the variant that diverges in the midline regime and returns `NaN` on degenerate segments, so keeping it was only ever justified while `bound_induced_velocity`'s corner behavior was worse. It no longer is.

**Applied**: all 3 sites in `src/panel_particle_wake.jl` (`PanelBufferFilaments`, `WakeBufferRings` ×4 edges, `BoundaryFilamentWrapper`) reverted to `bound_induced_velocity`, with their now-obsolete explanatory comments replaced by short notes on why the near-endpoint regime is hit there at all. Reproduction case re-run after the revert: clean, no `NaN`, matching the last row of the I4 table bit-for-bit. Full suite re-run: same single known `fmm_test.jl:99` failure, same values (that test doesn't touch this file).

`bound_induced_velocity_legacy` now has no callers anywhere in `src/` or `test/` — it is the formula `master` still ships as `bound_induced_velocity`, kept as a reference/comparison point. Left in place rather than deleted (it predates this work; deleting it is the maintainer's call).

The residual gap to the clean baseline (~8× on `maxσ`) is not necessarily a defect: as noted in F3, the baseline predates the `set_gradient!` buffer-alignment fix and runs under an entirely different FastMultipole accuracy regime, so it is no longer a like-for-like comparison.

---

## J. The `ProbeSystemStatic` warning, and the FLOWVPM buffer misalignment it exposed (2026-07-27)

**Result: the warning is gone, and the reproduction case now reproduces the clean baseline.** This closes the last remaining item from section A's caveat and supersedes section G's "residual gap".

### J1. Why the warning was more than cosmetic

`FastMultipole.ProbeSystemStatic` implements none of the previous-influence metadata interface, and no built-in FastMultipole target type does — the only reference implementations are in FastMultipole's own tests (`test/gravitational.jl:84-93`, `test/vortex.jl:143-152`).

The consequence is not confined to probes. `update_min_influence_leaf!` (`FastMultipole/src/tree.jl:1844`) takes the **minimum** over every target system in a branch, and `metadata_value(buffer, 0, i)` returns `0` for any system that hasn't opted in. So every tree branch holding both particles and probes had `min_gradient = 0`, and `ε = max(branch.min_gradient*RET, AET)` (`src/translate.jl:553`) collapsed to `AET` — i.e. the probes were silently nullifying the FLOWVPM metadata port from section A wherever the two target systems shared a branch.

### J2. Where the fix belongs

Implemented in **VortexLattice**, as a VortexLattice-owned type (`src/probes.jl`), for two reasons:

1. Every one of the ~48 FastMultipole interface methods VortexLattice defines is on a VortexLattice-owned type (`System`, `FilamentWrapper`, `PanelBufferFilaments`, `WakeBufferRings`, `BoundaryFilamentWrapper`). There is no type-piracy precedent in this repo, and the established idiom is "define your own type, implement the interface".
2. `ProbeSystemStatic` has nowhere to *store* a previous-influence estimate. Implementing the metadata methods on it — in either package — could only ever report the freshly-zeroed current values, which silences the warning without restoring anything.

`VortexLattice.ProbeSystem` therefore adds `previous_potential` / `previous_gradient` fields alongside the standard four. `FastMultipole.reset!` rolls the accumulated influence into them before zeroing, so each FMM pass automatically carries the preceding pass's influence; `seed_previous_influence!` copies them into the temporary `probes_active` packings in `panel_particle_wake.jl`. Field names match `ProbeSystemStatic` exactly, so all 46 existing `.position` / `.gradient` / `.hessian` / `.scalar_potential` accesses were left untouched — only the two struct field declarations (`src/system.jl:112`, `src/fluid_domain.jl:26`), the two constructors, and the 5 `probes_active` sites changed.

### J3. What this exposed: FLOWVPM has the section-F3 bug, unfixed

Switching the probes to `metadata_per_body = 2` initially appeared to be a spectacular win — the 36-step case dropped from `maxσ = 4.02` to `0.40`, right on the clean baseline. **It was not a win.** Two attribution runs showed why:

| configuration | `maxσ` | `maxΓ` |
|---|---|---|
| new type, previous influence seeded | 0.4034 | 0.2771 |
| new type, seeding disabled (metadata reports 0) | 0.4034 | 0.2771 |
| new type, metadata API opted out entirely (`NM=0`) | 4.826 | 1.262 |

Seeding made *no* difference; the whole effect came from `NM` changing the target buffer layout. Results must not be sensitive to that — and they were, because `FLOWVPM/src/FLOWVPM_fmm.jl` carried the exact anti-pattern section F3 fixed on the VortexLattice side and never fixed on the FLOWVPM side:

- `:143` `fmm.set_gradient!(target_buffer, j_target, val)` — bare 3-arg form, hardcoded to rows `5:7`, with `derivatives_switch` in scope and discarded.
- `:167` `fmm.set_hessian!(target_buffer, j_target, val)` — same, rows `8:16`.
- `:177` `buffer_to_target_system!` — bare 2-arg `get_gradient` / `get_hessian`, same hardcoded rows.

These FMM calls run with `scalar_potential=false` (FastMultipole's default), so `PS=false`, and FLOWVPM's own section-A metadata port set `NM=2` for `ParticleField`. The real `gradient_range` is therefore **6:8**, not `5:7`. So FLOWVPM's port had silently introduced a one-row misalignment in its own particle velocity/hessian accumulation — and VortexLattice's (correctly switch-aware, post-F3) writes into particle targets were being read back one row off. Section A's own diagnosis, that porting FLOWVPM was "necessary but not sufficient", was right for a reason nobody had identified: the port itself was half-finished.

**Fixed** in `FLOWVPM/src/FLOWVPM_fmm.jl`: all three sites now name and forward `derivatives_switch` to the switch-aware setters/getters, with the reads guarded on the switch's `VS`/`GS` parameters (the switch-aware getters `throw` when the corresponding output is disabled, where the legacy forms silently read garbage rows).

### J4. Validation after both fixes

Layout sensitivity is gone — `NM=0` and `NM=2` now agree to within the genuine accuracy-control difference rather than by an order of magnitude:

| configuration | `maxσ` | `minσ` | `maxΓ` | `mean|Γ|` |
|---|---|---|---|---|
| Clean baseline (`FastMultipole@c8cf796`, unmodified everything) | 0.40 | ~0.22 | 0.28 | — |
| **All fixes, probes opted in (`NM=2`)** | **0.4025** | **0.2234** | **0.2769** | 0.0189 |
| All fixes, probes opted out (`NM=0`, warning returns) | 0.5133 | 0.1839 | 0.2816 | 0.0239 |
| Before this section (section I state) | 3.41 | 0.0541 | 0.6366 | 0.0684 |

The reproduction case now matches the pre-regression baseline to ~1%, from ~8.5× off. Test suites: VortexLattice unchanged (same single known `fmm_test.jl:99` tight-tolerance residual, no new failures, no warnings anywhere in the log); FLOWVPM's own suite passes clean (14 testsets, exit 0) with the `set_gradient!`/`set_hessian!`/`get_gradient`/`get_hessian` change.

### J5. Honest note: the previous-influence carrying is currently inert

Verified the rollover works — after 6 steps all 95 probes hold nonzero `previous_gradient` (`0.011` to `0.82`, mean `0.25`). But seeding still changes nothing, and the reason is the tolerance settings, not a defect: `FLOWVPM.FMM` defaults to `relative_tolerance = absolute_tolerance = 1e-3` (`FLOWVPM/src/FLOWVPM_particlefield.jl:62`), so `min_gradient*RET` lands in `1.1e-5 … 8.2e-4`, always below the `AET = 1e-3` floor, and `ε = max(min_gradient*RET, AET) = AET` regardless. Relative error control is inactive in this case *by tolerance choice*.

It is still worth having implemented rather than reporting zeros: the interface is now honest, the warning is gone for a real reason, and the machinery engages automatically for any case with larger velocities or a smaller `absolute_tolerance`. But nobody should expect the metadata port alone to change results at these settings — what actually mattered here was J3's buffer alignment.

### J6. Remaining known gaps

- `FastMultipole.ProbeSystemStatic` / `ProbeSystemArray` still don't implement the metadata interface. VortexLattice no longer uses them, so the warning is gone here, but the upstream gap noted in section A stands and is worth a PR alongside issue #37.
- `test/fmm_test.jl:88` still constructs a `FastMultipole.ProbeSystem` directly. Left alone deliberately — it exercises the FastMultipole-owned target path, and that test passes no `error_tolerance`, so it emits no warning.

---

## K. The last failing assertion (`fmm_test.jl:99`) — diagnosed and the test rewritten (2026-07-27)

Sections E/F3/I/J all carried this forward as "a tight-tolerance residual, believed unrelated". It is unrelated, and now understood exactly. **The full test suite exits 0 for the first time in this investigation.**

### K1. Raising the expansion order does nothing — the error is already converged

Error of `fmm!` against the direct result, sweeping expansion order at the shipped `core_size = 1e-3`:

| `p` | 4 | 8 | 12 | 16 | 20 | 24 | 30 |
|---|---|---|---|---|---|---|---|
| error | 1.0e-6 | 4.2e-10 | 6.16e-10 | 6.1555e-10 | 6.1555e-10 | 6.1555e-10 | 6.1555e-10 |

Flat from `p=12` onward, to 11 significant figures. This is not multipole truncation.

### K2. The floor is a model difference, not a bug

`body_to_multipole_vl!` expands a **singular** vortex filament; `direct!` applies the **finite-core** regularization. At `core_size > 0` the two are not the same physical quantity, so no expansion order can reconcile them. At this evaluation point:

```
|v(finite_core=true) - v(finite_core=false)| = 1.234e-9
```

which brackets the observed plateau (6.16e-10 against the cored reference, 9.46e-10 against the singular one — `fmm!`'s answer is a blend, since the dual-tree sends some panels direct/cored and some through m2l/singular). The shipped `atol=1e-12` was unachievable by construction, and was presumably set before the core-model asymmetry existed.

Confirmed by removing the asymmetry: zero the panel core sizes *after* the analysis (`update_surface_panels!` with `fcore=(c,Δs)->0.0`; note `System(...; core_size=0.0)` alone is insufficient — `steady_analysis!` rebuilds the panels with its own `fcore`), and the multipole error converges exponentially straight to roundoff:

| `p` | 4 | 8 | 12 | 16 | 20 | 30 |
|---|---|---|---|---|---|---|
| error | 1.0e-6 | 6.2e-10 | 4.5e-13 | 8.7e-16 | 8.3e-16 | 8.3e-16 |

So the m2l path is correct, and — with the section-F3 switch-aware `set_gradient!` fix in place — accurate to machine precision.

### K3. Test rewritten

`test/fmm_test.jl` now:
- keeps the exact `FastMultipole.direct!` comparison at `atol=1e-12` (passes to 2.4e-18);
- checks the finite-core `fmm!` against an honest `atol=1e-8`, i.e. the core-model floor, with a comment explaining why a tighter tolerance is unreachable;
- adds a zero-core convergence test that asserts the multipole error is **monotonically decreasing** in expansion order and crosses `1e-5 / 1e-8 / 1e-11 / 1e-14` at `p = 4 / 8 / 12 / 20`.

That last check is strictly stronger than the assertion it replaces: a fixed tolerance at finite core could only ever measure the core-model gap, whereas this exercises the m2l path down to machine precision and would catch a buffer misalignment or a broken expansion immediately.

Suite status: **all testsets pass, exit code 0** (`FMM induced velocity`: 8/8).
