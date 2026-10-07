# VILMA2 — performance assessment & roadmap

A whole-model performance review aimed at **real transient runs** (a full glacial
cycle: thousands of time steps with a migrating coastline and a laterally varying
load), *not* the one-shot validation benchmarks. It records where the time goes,
the prioritized options, what has been implemented, and what is deferred.

For the per-step measured anchor and the lmax-scaling extrapolation see
[`performance.md`](performance.md); for the per-degree solver micro-optimisations
(band LU, degree-grouped memory, skip-negligible, OpenMP over the degree loop) see
§Performance of [`design.md`](design.md). This note sits on top of both.

## Measured results (validated this session)

The flag and warm-start changes have now been built and timed on the target
machine (Mac.fritz.box, 10-core Apple Silicon, `OMP_NUM_THREADS=8`). Anchor =
Martinec-2018 **E2**, lmax 128, 750 steps. `user` is total CPU work (the
low-noise measure); `real` is wall clock.

| build / mode | real | user | ms/step | notes |
|---|---|---|---|---|
| `-O2` (controlled baseline) | 58.8 s | 239 s | 78 | parallelizes ~4× (`user/real`≈4) |
| `-O3 -mcpu=native -funroll-loops -ffast-math` | **45.6 s** | 175 s | 61 | **1.29× wall / 1.36× CPU** vs `-O2` |
| `-O3` + warm start | 41.0 s | 172 s | 55 | warm adds ~1.5% CPU on E2 (see §2) |

**Three corrections to the prior notes fall out of this:**

1. **The flag speedup is ~1.3×, not 1.5–3×.** Real and free, but the lower bound.
2. **The model is *not* serial-bound at lmax 128.** Both `-O2` and `-O3` show
   `user/real ≈ 4` — the degree-loop OpenMP delivers ~4× here. This contradicts the
   "effectively serial, SHT-bound" diagnosis below and in [`performance.md`](performance.md).
3. **The `performance.md` anchor (162.5 s / 216 ms-step) does not reproduce.** The
   controlled `-O2` number on this machine is 58.8 s / 78 ms-step. The old run was
   likely serial (`openmp=0` / `OMP_NUM_THREADS=1`) or otherwise differently
   configured; treat `performance.md`'s extrapolation table as stale until rebased
   on these numbers.

Correctness under `-ffast-math`: full `make check` = **21/21 pass**, SLE
mass-conservation residuals unchanged at ~1e-16, E2 figure errors bit-identical.

## Cost model

For a transient run the wall time is

```
total  =  (timespan / dt)  ×  per_step
per_step  ≈  begin_step  +  SLE_fixed_point
SLE_fixed_point  ≈  n_outer × n_inner × (~3 SHTs per inner iteration)
```

- `begin_step` (`vilma_response`) does ~2 real banded solves per active `(l,m)`; the
  band LU + OpenMP-over-degree work already made this cheap (~58 ms at lmax 128).
- The **SLE fixed point's spherical-harmonic transforms dominate the *serial*
  fraction.** SHTns is linked serial and the per-step SLE cost is roughly linear in
  `n_outer × n_inner`; SHT cost scales ≈ O(lmax³). **But the step as a whole is not
  serial-bound at lmax 128** — measured `user/real ≈ 4` (the degree-loop OpenMP in
  `begin_step` parallelizes well), so wall time ≈ 61 ms/step at `-O3`, *not* the
  216 ms/step quoted in `performance.md`. The serial SHTs cap the *achievable*
  speedup (Amdahl), which is why OpenMP-SHTns (#4) still matters at higher lmax —
  but they do not make the present run single-core.

  > ⚠️ **Superseded (September 2026) — the paragraph above is wrong in three
  > ways.** It reasons from an end-to-end anchor; the model is now directly
  > instrumented, and [`performance.md`](performance.md) carries the measured
  > numbers. In short:
  >
  > 1. The SLE transforms do **not** dominate. The residual `SLE + coupling`
  >    bucket is ~61 % the adaptive stepper's rollback snapshot
  >    (`response_save_state`), and `sle_solve` itself splits ≈48 % transforms /
  >    ≈50 % grid-space work.
  > 2. Cost scales as **ℓmax^1.98** end to end (1-D) and **ℓmax^2.15** for the
  >    bucket — not O(ℓmax³), which is itself the evidence it is not
  >    transform-bound.
  > 3. **SHTns is not linked serial.** `config/common.mk` links `-lshtns_omp`
  >    whenever `openmp=1` (the default); what is missing is the
  >    `shtns_use_threads()` call, which appears nowhere in `src/`. The threaded
  >    library is linked and dormant.
  >
  > Amdahl's law is also the wrong frame here: on a Levante node the limit is
  > NUMA memory bandwidth, and thread *placement* alone is worth up to 1.6× at
  > fixed core count.

So the two highest-value levers for transient runs are **(a) the step count** (the
`dt` lever) and **(b) the SLE iteration count × per-SHT cost**. The benchmarks hide
both: they don't pay the long step count, and their clean coastlines converge the
SLE in few iterations.

## Prioritized options (3-agent review)

Ranked by return on effort for transient runs. ✅ = implemented now; ⏭ = deferred.

| # | Option | Payoff | Effort/risk | Status |
|---|--------|--------|-------------|--------|
| 1 | **Optimization compiler flags** (`-O3 -mcpu=native -funroll-loops -ffast-math`) | **1.3× measured** (free) | low | ✅ |
| 2 | **Warm-start the SLE fixed point** from the previous step's `rsl` | **~1.5% on E2**; up to ~2× only on strongly-migrating coastlines (unverified) | low | ✅ |
| 3 | ~~**ETD0** exponential memory update → larger `dt`~~ | **rejected — fails benchmarks** | — | ✗ |
| 3b | ~~**ETD1** (linear-strain φ-weights) → larger `dt`~~ | **rejected — the memory *rule*, not the strain coupling, sets the order** | — | ✗ |
| 3c | **Trapezoidal memory rule solved by coupling iteration** → 2nd-order, with an adaptive-`dt` controller (1-D + field driver 3a + step-doubling + SLE-coupled 3b + `vilma_timestep` controller — all DONE) | **order 1→2; ~1300× accuracy (1-D) / ~270× (SLE, fast load) at fixed `dt`; adaptive: ~1.6× wall on dynamic-range loads (no gain on smooth loads)** | med | ✅ |
| 4 | **OpenMP SHTns as the default** (offline *and* coupled) | several× at lmax ≥ 256 | low | ⏭ |
| 5 | **Fuse the two syntheses** to one synthesis of `N_lm − u_lm` | ~33% of inner-loop SHTs | low | ⏭ |
| 6 | **Batched multi-RHS band solve** over all `m`/re-im at fixed `l`; kill the dense degree-1 LU via nullspace projection | high at production lmax | med | ⏭ |
| 7 | **Multi-rate cadence** — re-solve the SLE less often than the cheap memory advance | ~2× | med | ⏭ |
| — | Memory footprint ~O(lmax²): ~0.5 GB @128, ~2 GB @256, ~8 GB @512 | feasibility gate | — | note |
| — | Restart I/O: write the full `tau_*` memory state only at checkpoints, not per step | avoid a multi-GB/step footgun | low | note |

Measured so far: flags + warm-start give ~1.3× on E2. ETD1 (the would-be step-count
lever) was tested and rejected (§3b), but it pointed the way: §3c then found the real
larger-`dt` lever — a **2nd-order trapezoidal memory rule solved by coupling iteration**
(order 1→2 in the 1-D stepper, validated; field-driver wiring + step-doubling pending).
The remaining large multipliers live there and in the deferred items (OpenMP-SHTns at
high lmax) and in real-coastline iteration counts the benchmark does not exercise.

## Implemented now

### 1. Optimization compiler flags — `.configme/machines/macbook.mk`

The only optimization flag was `-O2`; no architecture targeting, unrolling, or
vectorization. The hot kernels (`dissipative_rhs`, `advance_memory`, the band
solve) are tight double-precision loops that benefit directly. Now:

```make
DFLAGS_NODEBUG = -O3 -mcpu=native -funroll-loops -ffast-math
```

- `-mcpu=native` is the aarch64 / Apple-Silicon idiom (the x86 `-march=native`
  equivalent); it enables NEON vectorization. On other machine fragments use the
  platform-appropriate flag (`-march=native` on x86-64).
- **`-ffast-math` is ON by design**, for parity with the CLIMBER-X production
  build (the host enables it, and the coupled build cannot reliably be configured
  differently). It relaxes IEEE semantics, so this is the one change that **must
  be validated on the target machine**: re-run `make check` and confirm the SLE
  mass-conservation tests (~1e-16) still pass.

**Measured:** built and validated on the target machine. `make check` = 21/21 pass,
SLE mass residuals unchanged at ~1e-16 under `-ffast-math`. E2 wall time
58.8 s → 45.6 s = **1.29× (1.36× CPU)** — real, but the low end of the estimate.

### 2. Warm-start the SLE fixed point — `vilma_sle.f90`, `vilma_coupling.f90`

`sle_solve` previously did `rsl = 0.0` at the top of *every* call, discarding the
previous step's converged solution. For a transient run the coastline and RSL
barely move between adjacent steps, so that solution is a near-converged seed.

- New `sle_solver%warm_start` flag (default `.false.` → cold start, so **every
  existing benchmark/test is byte-for-byte unchanged**).
- `rsl` changed from `intent(out)` to `intent(inout)`; the cold-start reset is now
  guarded by `.not. warm_start`. With the flag off the reset still runs, so the
  semantics for existing callers are identical.
- The coupling driver (`solid_earth_init`) turns it **on**: `self%rsl` persists
  across sub-steps and coupling intervals (seeded to 0 at init), so each solve
  starts from the last converged field.

The converged fixed point is unique (the SLE is a Fredholm equation of the second
kind → a contraction) and mass is rebalanced every iteration via Δφ, so the
**answer is unchanged to tolerance** — only the inner iteration count drops
(typically to 1–2, and it tames the pathological strongly-migrating-coastline case
that can otherwise need tens of iterations).

**Measured:** E2 is now warm-start-by-default (`test_benchmark_sle`); all four cases
(C2/D3/E2/F1) give **bit-identical figures** to cold start. Last-pass inner iters
drop **1.45 → 1.00/step** (the floor). But the E2 wall-time/CPU win is only **~1.5%**,
because E2's clean spherical cap already converges in ~1.45 inner iters cold — there
is almost nothing to save, and the inner loop is a small part of per-step cost (3
outer coastline passes + `begin_step` dominate). **The "~2×" upside is therefore
*unverified*** and will only appear on a strongly-migrating *real* coastline
(ICE-6G-style), which no benchmark exercises — consistent with caveat #1 below.
Warm-start is kept on regardless: free, correct, and the right default for real
domains. `VILMA_SLE_WARM=0` forces cold start for the A/B.

## Deferred — roadmap (next sessions)

**3. ETD0 — REJECTED (reproduced failure, this session).** ETD0 is the exact
held-strain exponential update: a one-line substitution in `advance_memory`
(`vilma_viscoelastic.f90:216-235`), `(1−M) → exp(−M)` and `2μM → 2μ(1−exp(−M))`. It is
unconditionally stable and was proposed here as the highest-value `dt` lever.
**It does not work.** Applied to the current code, it worsens benchmark agreement
and fails a test:

| metric | forward-Euler | ETD0 | tol |
|---|---|---|---|
| disc VE centre, t=1 kyr (rel) | ~0.5% | **2.19%** (−68.6 m vs Spada −70.1) | — (passes, degraded) |
| Martinec-2018 A geoid | 1.94% | **4.48%** | 3% → **FAILS** |
| Martinec-2018 A uplift / horizontal | 0.52 / 0.69% | 0.83 / 0.91% | 3% |

The root cause is *accuracy, not stability*: ETD0 is exact only if strain is held
constant over the step, but in the fast early transient strain varies within a step,
so ETD0 under-relaxes (weight `1−exp(−M) < M`) — the **wrong direction** for the
known long-time under-relaxation. The Spada/Martinec benchmarks are the `dt→0`
normal-mode truth, and forward-Euler is empirically closer at our `dt`. (This had
been recorded only in project memory, never committed — hence the earlier "no trace
in the repo" note; the experiment has now been re-run and confirmed.)

**3b. ETD1 — IMPLEMENTED, MEASURED, REJECTED (this session).** ETD1 uses linear-
strain φ-function weights (`τ_{n+1} = e^{−M}τ_n − 2μM[(φ₁−φ₂)ε_n + φ₂ε_{n+1}]`,
`φ₁=(1−e^{−M})/M`, `φ₂=(M−1+e^{−M})/M²`), which account for the within-step strain
variation that breaks ETD0. It was implemented behind a scheme flag in the 1-D
stepper (forward-Euler stays default; `make check` byte-identical) and swept against
a converged reference (`test_etd1`, the standalone characterization; FE and ETD1
agree to 2e-4 at `dt`=0.1 yr, confirming a shared `dt→0` limit). **It does not help:**

| `dt` [yr] | M | FE err | ETD1 err | FE/ETD1 |
|--:|--:|--:|--:|--:|
| 25 | 0.11 | 5.8e-3 | 4.6e-2 | 0.12 |
| 100 | 0.44 | 2.3e-2 | 1.6e-1 | 0.14 |
| 1000 | 4.4 | 4.1e-1 | 7.1e-1 | 0.58 |
| 2000 | 8.8 | 1.3e0 | 8.7e-1 | 1.50 |

In the usable (resolved, `M<1`) regime ETD1 is ~8× **less** accurate than FE; both
are 1st-order in `dt`. **The root cause is the explicit strain↔memory coupling:** the
strain fed to the memory update comes from the *previous* step's memory (lagged), so
it is only 1st-order accurate — the coupling, not the memory integrator, is the order
bottleneck. ETD1's higher-order memory treatment is wasted, and its exponential
under-relaxes per step (forcing weight `2μMφ₁ < 2μM`) — the same "wrong direction"
that sank ETD0. Separately, **forward-Euler is practically unconditionally stable for
these models** (finite to `M≈35`; the elastic/self-gravity feedback damps the naive
`M<2` scalar limit), so FE's stability is not a real constraint and ETD's only edge —
boundedness past `M≈8` — lands in the under-resolved regime. ETD's unconditional
stability would matter only with a genuinely weak (low-η) layer, which no benchmark
model has.

**Where the larger-`dt` lever actually is (for the adaptive-stepping goal):** not the
memory integrator family that ETD belongs to, but a higher-order *memory rule*. The
ETD1 study's stated cause — "the explicit strain↔memory coupling is the 1st-order
bottleneck" — was half right and is corrected in §3c below.

**3c. TRAPEZOIDAL MEMORY RULE + COUPLING ITERATION — IMPLEMENTED, MEASURED, 2ND-ORDER
(1-D stepper, this session).** The decisive observation: in this affine scheme the
surface observable is an *exact algebraic function* of the memory `τ` (the balance
solve `K ε = load + D τ` is exact given `τ`), so **the observable's convergence order
equals the order of the memory time-integration alone.** Two consequences, both
measured (`test_couple_order`, held degree-2 load, homogeneous Maxwell sphere, resolved
regime `M < 0.4`):

| variant | order | err @ `dt`=5 yr |
|---|---|---|
| A — forward-Euler (explicit, historical) | 1.00 | 1.1e-3 |
| B — backward-Euler, coupling iterated to consistency | 1.00 | 1.1e-3 |
| C — trapezoidal, single predictor (no fixed point) | 0.99 | 4.2e-3 |
| **D — trapezoidal, coupling iterated** | **2.00** | **8.3e-7** |

- **Iterating the coupling is NOT itself an order lever.** Variant B makes the coupling
  fully implicit/consistent but keeps a 1st-order rule → stays 1st-order. This corrects
  the §3b framing: the original explicit scheme is *not* balance-lag-limited (it already
  uses current-time memory `τ_n` in the balance at `t_n`); it is simply forward-Euler on
  the memory ODE. The order bottleneck is the **memory rule**, not the strain coupling.
- **2nd order comes from the trapezoidal (Crank–Nicolson) rule (D)** — A-stable, and
  ~1300× more accurate than FE at fixed `dt` here. But trapezoidal is *implicit* in the
  end-of-step strain, so it must be solved by a Picard fixed point over the endpoint
  balance. **So the two original §3c "levers" are one:** an implicit 2nd-order advance
  *solved by* coupling iteration. The fixed point is cheap — `couple_tol=1e-6` reaches
  the trapezoidal truncation floor in **5–8 iterations** (the measured knee; tightening
  past it buys iterations, not accuracy).

Implementation: `SCHEME_TRAP` (+ `SCHEME_BE` as the 1st-order control) in
`advance_memory`, and a restructured `ve_step` that separates the time-aligned REPORT
(always the balance against the entering `τ_n`, exactly elastic on the first call) from
the ADVANCE (explicit single-pass for FE/ETD1 — byte-identical; iterated endpoint for
BE/TRAP). Default unchanged (`scheme=FE`, `max_couple_iter=1`); `make check` 21/21.

**Step-doubling (part ii) — DONE.** `ve_degree%step_double` gives the Richardson
local-error estimate a controller needs: advance one Δt by two Δt/2 sub-steps (the
kept result) and scale the coarse/fine memory difference by `2^p−1` (3 for
trapezoidal). M = μΔt/η is linear in Δt, so halving Δt just halves Mk — no re-init.
Validated to scale as `dt^(p+1)`: measured `dt^2.00` for FE (p=1), `dt^2.99` for
trapezoidal (p=2). Estimate-only; the accept/reject + Δt controller is the remaining
piece.

**Field driver (part 3a) — DONE, behind a flag.** Trapezoidal+iteration is wired
into the per-`(l,m)` `ve_response` (default FE = unchanged). The report (`apply`, via
`begin_step`'s frozen drift from τ_n) was already time-aligned; only `commit_step`
changes — its drift solve is factored into a reusable `solve_drift(τ)`, and the TRAP
commit Picard-iterates the endpoint with the load **frozen** at the converged σ
(re-solve drift against trial τ_{n+1}, trapezoid-advance from a τ_n snapshot, repeat
to `couple_tol`). Validated against the 1-D `ve_degree` TRAP to 3.6e-13 (displacement)
/ 1.4e-16 (geoid) in `test_ve_response`; FE byte-identical; `make check` 21/21.
**Measured E2 cost: ~6.0× per-step at fixed `dt`** (serial, lmax 128, 750 steps,
`couple_tol`=1e-6; ~8–10 drift re-solves/step) — the overhead the larger adaptive
step must amortize. E2 still passes to tolerance under the frozen-load approximation.

**SLE-coupled endpoint (part 3b) — DONE, behind `sle_solver%max_mem_iter`.** For a
**fast-evolving load** σ and τ_{n+1} are a mutual fixed point (the ocean load depends
on the relaxation state, and the trapezoidal endpoint depends on the load), so 3a's
frozen-σ endpoint is only 1st-order. 3b lifts the memory advance into `sle_solve` as a
single combined Picard loop: `prepare_endpoint` (snapshot τ_n) once, then each pass
re-converges the water load against the current τ_{n+1} estimate and calls
`advance_endpoint` — one trapezoid pass that **refreshes the report drift to the new
τ_{n+1}** so the next σ-convergence (and coastline migration) sees the advanced memory
— until the surface drift settles; `finalize_step` advances time. `commit_step` and
3a's frozen endpoint are kept intact (and DRY-shared via `trapezoid_advance_all` /
`vilma_advance`) for the held-load path that `test_ve_response` pins. FE / elastic / null
report converged after one pass, so the loop is inert for them and the FE default is
byte-identical (`make check` 21/21).

A latent bug surfaced here and was fixed: the trapezoidal rule ½(ε_n+ε_{n+1}) needs
the **start-of-step** load σ_n in ε_n, but the code used the current σ_{n+1}. Invisible
for a held load (σ_n = σ_{n+1} — why 3a, `test_couple_order` and `test_ve_response`
never caught it), but an O(`dt`) error for a varying load that caps the order at 1.
Tracking σ_n (`sigma_n`, primed after the first step so the held path stays
byte-identical) restores 2nd order. Validated by `test_sle_couple_order` (a fast ice
ramp through the full SLE under fixed ocean to isolate the integrator): a single
combined pass is **order 0.98**, co-convergence is **order 2.08** and ~270× more
accurate at the finest `dt` — the trapezoidal 2nd order carried through the driver.

**Adaptive-`dt` controller — DONE (`vilma_timestep`).** The accept/reject + step-size
loop on top of the step-doubling estimate, isolated in its own module (room for more
than one stepping strategy). `adaptive_stepper%advance(t0,t1)` crosses a coupling
interval with the ice load **linearly interpolated** between its endpoints, choosing
Δt by step-doubling on the Maxwell memory: each candidate step is taken once at Δt and
once as two Δt/2 sub-steps; the coarse/fine memory ∞-norm difference / (2^p−1) is the
local error, scaled by `atol + rtol·‖τ‖`. Accept (keep the more accurate fine state,
grow Δt by `safety·err^{-1/(p+1)}`) or reject (roll back, shrink). Rollback uses new
`ve_response` primitives `save_state`/`restore_state` (buffer A = τ_n) +
`stash_coarse`/`coarse_fine_error` (buffer B = τ_coarse); Δt changes via `set_dt`, a
cheap `Mk = (μ/η)·Δt` rescale (no operator re-factor — the band LU is Δt-independent).
`test_timestep`: the field step-doubling estimate scales as Δt³ (order 2.97), and on a
fast ice ramp the controller converges to a fine fixed-Δt reference to ~1e-4 of signal,
the error falling monotonically as `rtol` tightens. The **first step is also 3rd-order**:
the start-of-step load σ_0 is seeded by a *report-only* SLE solve at t=0 (`sle_solve`'s
`report_only`/`sigma_lm` args — converge the load against the rest memory, no advance,
so σ_0 is the elastic-consistent load), via `ve_response%prime_sigma`. Without it the
first step falls back to the σ_{n+1} proxy (O(Δt) for the relaxing SLE load); seeding σ_0
makes that step 3rd-order **and** removes the startup reject storm (the controller no
longer wastes steps shrinking through an inaccurate first step).

**When adaptivity pays — cost vs accuracy (`test_timestep` part C; cost = SLE solves,
the dominant unit; adaptive ≈ 3 solves/step for the coarse+2×fine estimate).** Two
takeaways, both honest:
- **Smooth / uniform dynamics → fixed wins.** On a load with no dynamic range (a single
  ramp), the 3×/step estimate overhead is never amortized: fixed Δt=100 yr reaches
  3.2e-5 rel-error in 20 solves, while adaptive `rtol=1e-2` gives a *worse* 2.2e-4 in 21.
  Don't pay for error control you don't need.
- **Dynamic range → adaptive wins.** With a fast 200-yr ramp then an 8-kyr hold, fixed
  is *forced* to the small ramp Δt for the entire slow hold (≥200 solves, and then
  over-accurate at 5e-7); adaptive grows Δt through the hold and hits a useful ~1e-4 in
  **57 solves / 0.68 s vs the fixed floor of 200 solves / 1.1 s** — ~3.5× fewer solves,
  ~1.6× wall-clock, the gap widening for looser tolerances. (The wall-clock win trails
  the solve-count win because a larger Δt means a larger M = μΔt/η and so more SLE↔memory
  co-convergence passes per solve — adaptive's big steps are individually dearer.)
- **Unconditional value:** error *control* to a target tolerance without hand-tuning Δt.
  Real ice-age loading (fast deglaciation + slow interglacial relaxation) is the
  dynamic-range regime, so the controller is the right default there; a steady forcing
  is better served by a fixed Δt.

The long-time under-relaxation that originally motivated §3c is a separate
`dt`/NMAX-resolution + slow-mode question, not an integrator issue. **§3c is now
complete** (1-D + 3a + step-doubling + 3b + controller); the default scheme stays FE
until a driver opts into the trapezoidal adaptive path.

`test_etd1` (and the scheme-pluggable kernel) are kept as the reproducible evidence
that ETD0/ETD1 are not re-attempted, and `test_couple_order` as the evidence for the
trapezoidal result.

**4. OpenMP SHTns as the default.** SHTns is currently linked serial
(`config/common.mk`) — a deliberate choice to avoid OpenMP nesting inside
CLIMBER-X. **Decision (this review): make the OpenMP SHTns variant the default for
both the standalone/offline driver and the coupled build.** The transforms are the
serial bottleneck and the per-SHT cost is O(lmax³), so threading them is the clean
lever that moves lmax 256 into comfortable range. The `vilma_sht` wrapper needs no
source change — only the build wiring (link the `shtns-omp` variant; it must be
built in `fesm-utils`). The model's own degree-loop OpenMP and the SHT calls run in
*different phases* of the step (never nested concurrently), so they can share one
thread pool safely. (Requires building the OpenMP SHTns variant in `fesm-utils`,
hence deferred to a build/dependency pass.)

**5. Fuse the two syntheses.** The inner loop only consumes `Sraw = N − u`;
synthesize the combined spectrum `N_lm − u_lm` once instead of `u` and `N`
separately, cutting inner-loop SHTs from 3 to 2 (~33%). Synthesize `u`/`N`
individually only on the converged iterate, for diagnostics. Pure reordering.
(`vilma_sle.f90:177-178`.)

**6. Batched multi-RHS band solve + degree-1 nullspace projection.** `begin_step`
calls the band solver one RHS at a time, twice per coefficient (re/im), but for a
fixed degree `l` the factored LU is identical across all `m` and both re/im.
Batching them into one multi-RHS triangular solve (`vilma_band.f90:114`) turns
memory-bound work into cache-friendly BLAS-3-style reuse — biggest at lmax ≥ 256.
Separately, degree 1's dense KKT border makes its solve O(nr³) and serializes one
thread; replace with a Sherman–Morrison / nullspace projection on the narrow band
(`vilma_radial_fe.f90`, bordered path). Validate against the degree-1 benchmark.

**7. Multi-rate cadence.** The viscous memory evolves on millennial timescales; the
expensive migrating-coastline SLE need not be re-converged at the same fine cadence
as the cheap Maxwell memory advance. The affine `begin_step/apply` split already
separates frozen drift from the current load, so the SLE can run on a coarser
cadence (e.g. per coupling interval) with the load frozen between. Complementary to
warm-starting (#2).

### Notes / guardrails (not bottlenecks today)

- **Memory footprint scales O(lmax²)** (six `(NLAM, ne, nk)` memory-stress arrays):
  ~0.5 GB at lmax 128, ~2 GB at 256, **~8 GB at 512** — a genuine feasibility wall
  at 512 on a laptop, and a real cost when this state is checkpointed. lmax 128 is
  the right coupled default; lmax > 256 would want a sparse/compressed memory layout
  (the skip-negligible machinery already avoids *solving* negligible coefficients
  but still *stores* them all).
- **Restart I/O.** The coupling `update` does no per-step I/O (good). But a full
  snapshot writes the entire `tau_*` memory state (multi-GB at high resolution).
  Enforce that heavy `tau_*` restart writes happen only at checkpoints; per-step
  diagnostics should use the small lon×lat field subset (`vilma_io` already supports a
  variable subset via its `nms` argument).

## Validation status of this change set

All validated on Mac.fritz.box (gfortran 15.2, `OMP_NUM_THREADS=8`):

- **Compiler flags:** `make check` 21/21 pass under `-ffast-math`; SLE
  mass-conservation residuals unchanged at ~1e-16. E2 1.29× wall / 1.36× CPU.
- **Warm-start:** now default-on in the benchmark; C2/D3/E2/F1 all give
  bit-identical figures to cold start. ~1.5% CPU win on E2 (see §2); larger benefit
  unverified pending a real migrating coastline. The coupling driver already opts in.
- **ETD0:** re-run and **rejected** (fails Martinec A geoid, 4.48% > 3%); see §3.
- **ETD1:** implemented behind a scheme flag (FE default, `make check` byte-identical)
  and **rejected** — ~8× less accurate than FE in the resolved regime; the memory
  *rule* (not the strain coupling) sets the order, and FE is practically
  unconditionally stable for these models (§3b). Kept as `test_etd1` + the scheme-
  pluggable kernel (reproducible evidence). Forward-Euler retained as the default.
- **Trapezoidal + coupling iteration (§3c):** implemented (`SCHEME_TRAP`/`SCHEME_BE`,
  iterated `ve_step`; FE default `make check` 21/21 byte-identical) and **adopted as
  the adaptive-`dt` core** — **2nd-order** in the 1-D stepper (`test_couple_order`:
  order 2.00, ~1300× more accurate than FE at fixed `dt`), the fixed point converging
  in 5–8 iterations at `couple_tol=1e-6`. The control (backward-Euler, iterated) stays
  1st-order, confirming the iteration is not itself an order lever.
- **Step-doubling (§3c ii):** `ve_degree%step_double` (1-D) Richardson estimate validated
  to scale as `dt^(p+1)` (`dt^2.00` FE, `dt^2.99` trapezoidal); the field-level estimate
  (via the controller primitives) scales as `dt^2.97` (`test_timestep`).
- **Field driver (§3c 3a):** trapezoidal+iteration wired into `ve_response` behind the
  scheme flag (FE default byte-identical, `make check` 21/21). Reproduces the 1-D
  `ve_degree` TRAP to 3.6e-13 / 1.4e-16 (`test_ve_response`). E2 cost ~6.0× per-step at
  fixed `dt` (frozen-load endpoint; passes E2 to tolerance).
- **SLE-coupled endpoint (§3c 3b):** σ↔τ co-convergence lifted into `sle_solve`
  (`sle_solver%max_mem_iter`; FE / elastic / null inert, `make check` 21/21 byte-
  identical). Fixes the start-of-step-load (σ_n) term in the trapezoidal ε_n that a
  varying load needs. `test_sle_couple_order` (fast ice ramp, fixed ocean): a single
  combined pass is order 0.98, co-convergence is order 2.08 and ~270× more accurate at
  the finest `dt` — trapezoidal 2nd order carried through the full driver.
- **Adaptive-`dt` controller (§3c):** `vilma_timestep%adaptive_stepper` — step-doubling
  accept/reject + step selection with linear load interpolation; `ve_response`
  save/restore/set_dt + `prime_sigma` primitives, `sle_solve` `report_only`/`sigma_lm`
  (`make check` 21/21 byte-identical, FE default untouched). `test_timestep`: estimate
  ~`dt^2.97` including the σ_0-primed first step; error monotone in `rtol`. Cost/accuracy
  (part C): adaptive does NOT beat tuned fixed Δt on a smooth load (3×/step overhead),
  but on a fast-ramp+long-hold load it reaches ~1e-4 accuracy in ~3.5× fewer solves /
  ~1.6× less wall-clock (fixed is pinned to the small ramp Δt through the slow hold).
