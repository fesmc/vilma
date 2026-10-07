# VILMA2 — performance note

Whole-step cost, where it goes, and what it means for real runs. For the
per-degree solver micro-optimisations (band LU, degree-grouped memory,
skip-negligible, OpenMP over the degree loop) see §Performance of
[`design.md`](design.md); for the roadmap see
[`performance-assessment.md`](performance-assessment.md).

> **This note was rewritten in September 2026** against direct instrumentation on
> DKRZ Levante. The previous version inferred the cost breakdown from an
> end-to-end anchor on a laptop and got the attribution wrong in three ways; the
> corrections are in §What changed and why, at the end, because the wrong version
> was quoted elsewhere and the record matters.

## Where the time goes — measured, not inferred

The driver reports `solid_earth_update` as **drift solve + memory advance + a
residual bucket** historically labelled `SLE + coupling`. That residual is now
decomposed by explicit `system_clock` regions (`vilma_sle`'s `t_sht`/`t_apply`/
`t_total`, `vilma_coupling`'s `t_remap`/`t_rot`, `vilma_timestep`'s `t_guard`), with an
`unattributed` remainder printed so a missed phase shows up as a number rather
than as someone else's cost. Set nothing: it is always on and costs one
`system_clock` pair per region (~25 ns).

ℓmax 64, 1-D, 16 threads, 81 coupling steps, Levante (2 × EPYC 7763):

| phase | share of the residual bucket |
|---|---:|
| `sle_solve` own work | 33 % — of which ≈48 % harmonic transforms, ≈50 % grid-space |
| stepper rollback bookkeeping | **61 %** |
| rotation (polar motion) | 4 % |
| host↔Gauss remap | <1 % |
| unattributed | 2 % |

**The bucket is mostly not the SLE.** The largest single term was
`response_save_state` — six full-array copies of the `(NLAM,ne,nk)` Maxwell
memory taken on *every* sub-step, to guard against a step rejection that did not
occur once in any run measured. It is now threaded (and the arrays are
first-touched in parallel, so their pages land on the NUMA domain that will read
them), which is worth **1.50× on the whole solver step**.

## Resolution scaling — ℓmax², not ℓmax³

Fitted over ℓmax 32…512 (1-D) and 32…256 (3-D), Levante, 128 threads:

| viscosity | fitted exponent | ℓmax 32 | 128 | 256 | 512 |
|---|---:|---:|---:|---:|---:|
| 1-D | **ℓmax^1.98** | 12.3 ms | 164.9 | 694.4 | 2985.5 |
| 3-D | **ℓmax^2.40** | 178.1 ms | 5829.7 | 26076.6 | — |

The residual bucket on its own scales as **ℓmax^2.15**. A Legendre-transform-bound
cost must approach ℓmax³, so **the model is not SHT-bound** — which is exactly
what the ≈50/50 transform/grid-space split inside `sle_solve` says directly.

## Thread scaling — the limit is memory bandwidth

One Levante node is 2 sockets × 64 cores in **eight NUMA domains of 16 cores
each** (`numactl -H`). With `OMP_PROC_BIND=close` threads pack onto consecutive
cores, so ≤16 threads share a *single* memory controller however many cores they
use. The measured curve follows the domain count, not the core count.

Same binary, same thread count, placement the only difference (ℓmax 128, 1-D):

| threads | `close` | `spread` | gain |
|---:|---:|---:|---:|
| 8 | 187.2 ms | **126.6** | 1.48× |
| 16 | 150.4 ms | **93.4** | 1.61× |
| 32 | 94.3 ms | **76.6** | 1.23× |
| 64 | 69.4 ms | 68.8 | 1.01× |
| 128 | 80.0 ms | 80.7 | 0.99× |

**16 threads spread beats 32 threads packed** — the same speed for half the cores.
The gain vanishes by 64, where four domains already supply more bandwidth than the
code can use; at 128 both placements regress equally, so that regression is not
placement but the cost of occupying every core on the node.

> **Run with `OMP_PROC_BIND=spread` and at most 64 threads.** The default `close`
> costs up to 1.6× at the thread counts an ensemble member actually uses.

Best end-to-end speed-up on one node: **9.2× (1-D)** and **13.2× (3-D)**, both at
64 threads, against a true single-thread baseline. Amdahl's law is *not* a useful
model for this curve — it assumes extra cores add only compute, when past 16
threads here they also add memory controllers; fitting `s + p/n` produces a
ceiling the measurements then exceed.

## Time to solution — the full last deglaciation

GLAC-1D (Tarasov) LGM→present, 261 coupling steps at 100 yr, ℓmax 128, RTopo-2
present-day reference, rotation on, migrating coastline. One Levante node:

| viscosity | solver | elapsed | setup | peak RSS |
|---|---:|---:|---:|---:|
| 1-D | 39.0 s | **87 s** | 1.5 s | 2.45 GB |
| 3-D (Bagge 2021) | 1584.6 s | **30 min** | 172.4 s | 10.07 GB |

A full glacial cycle at this resolution is therefore minutes (1-D) to a few hours
(3-D) — not the tens of hours the previous version of this note projected.

## Caveats

1. **SLE iteration count on real coastlines.** Per-step cost is roughly linear in
   `n_outer × n_inner`. The outer loop now **exits when the coastline stops
   moving** rather than always running `n_outer` passes; on the deglaciation that
   is 1.4 passes/solve instead of 3.0. A pathological coastline that oscillates
   rather than settles would still hit the `n_outer` cap — `sle_result%n_coast_flip`
   reports which case you are in.
2. **These are `ifx -Ofast -march=znver3` numbers on Levante.** Absolute times on
   other machines will differ; the scaling exponents and the phase attribution
   should not.
3. **I/O.** The timing runs write every coupling step, so their `vilma_write_step`
   share is an upper bound on a production run writing every fifth.

## Reproducing

```sh
make vilma                                 # ifx, OpenMP, -Ofast
OMP_NUM_THREADS=64 OMP_PROC_BIND=spread OMP_PLACES=cores \
    bin/vilma.x <your deglaciation namelist>
```

The `[PROFILE]` blocks in stdout give the setup/transient split, the
drift/memory/residual breakdown, and the residual bucket opened up. The driving
scripts and the analysis that produced the tables above live in the companion
experiments repository (`experiments/run_sweep2.sh`, `analysis/run_timing.jl`).

## What changed and why

The previous version of this note inferred the breakdown from a single
end-to-end anchor (162.5 s for Martinec E2 at ℓmax 128, on a 10-core laptop) and
drew three conclusions that direct measurement contradicts:

| previous claim | measured |
|---|---|
| "the bottleneck has moved to the SLE fixed-point's spherical-harmonic transforms" | the residual bucket is **61 % stepper rollback bookkeeping**; `sle_solve` is ≈half grid-space work |
| "SHT cost scales ≈ O(lmax³)" | **ℓmax^1.98** (1-D end to end), **ℓmax^2.15** (the bucket) |
| "the run is effectively serial… a bigger machine does not help as-is" | **9.2× / 13.2×** on one node; the limit is NUMA bandwidth, and placement alone is worth 1.6× |

One factual correction as well: this note and `performance-assessment.md` both
said **"SHTns is linked serial (deliberately, to avoid OpenMP nesting)"**. It is
not — `config/common.mk` links `-lshtns_omp` whenever `openmp=1`, which is the
default. What is true is that `shtns_use_threads()` is **never called** anywhere
in `src/`, so SHTns keeps `omp_threads=1` and selects its serial kernels. The
threaded library is linked and dormant; enabling it is one line, and remains an
untested lever rather than a measured one.
