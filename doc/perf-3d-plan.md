# 3-D viscosity performance — findings

This file originally proposed a *fused batched-transform* redesign on the premise that
the per-element spherical-harmonic transforms dominated the 3-D memory advance. That
premise turned out to be **wrong**. This is the corrected record.

## TL;DR

- **Batching the dyadic transforms was rejected** — implemented and verified correct, but
  it *regressed* the advance ~1.8× (SHTns `set_many` gives little on CPU; the assembly/
  scratch overhead dominates at lmax 64). Reverted.
- The 3-D per-step cost is **explicit sub-stepping × the memory advance**, not a slow
  single transform. Output and input-remap are negligible (<1 % each).
- We were profiling the **wrong viscosity field** (Pan 2022). The intended default is
  **Bagge 2021, floored at log10 η = 19.5** (= VILMA1's CLIMBER-X default), which is ~6×
  cheaper and physically sane (Pan's raw min ~13.7 dex is unphysical).
- The real wall-time killer was **`solid_earth_init` auto-tuning 8 SHTns configs (~136 s)**; fixed
  with `SHT_QUICK_INIT` (→ ~1 s), which unblocked lmax 128.

## Cost structure (measured)

A coupling step is `solid_earth_update`, which is `stepper_advance` (vilma_timestep). For the explicit
`scheme=fe`, the interval `span` is split into `n_sub = ceil(span / (cfl/max_rate))` equal
sub-steps; each does one SLE solve + one Maxwell memory advance. `max_rate` is set by the
**stiffest (lowest-η) Maxwell cell**, so the viscosity floor controls `n_sub`.

Full-step profile (`vilma_drive` PROFILE timers, 8 threads, Bagge floor 19.5, `SHT_QUICK_INIT`):

| | lmax 64 | lmax 128 |
|---|---|---|
| `solid_earth_init` (one-time) | 0.1 s | 1.2 s |
| `build_remap` (one-time) | 7.5 s | 13.9 s |
| **`solid_earth_update` / 100-yr step** | ~0.81 s | ~4.41 s |
| sub-steps (`n_solve`) | 8 | 8 |
| laterally-3-D elements | 51 / 217 | 66 / 217 |
| `read_ice` / `vilma_write_step` | 9 / 11 ms | 9 / 36 ms |

≈ 8 ms/simulated-yr (lmax 64), ≈ 44 ms/yr (lmax 128). VILMA1's CLIMBER-X cadence is
`n_year_geo = 10`, so ~1 sub-step per VILMA1 update (~0.5 s/update at lmax 128).

For contrast, Pan 2022 at the same floor: `solid_earth_update` ≈ 3.9 s/step at lmax 64, `n_sub` 23,
116/217 elements 3-D — i.e. the field choice was most of the apparent slowness.

## What was changed

- `vilma_sht`: `SHT_QUICK_INIT` instead of `SHT_GAUSS` — skips SHTns's ~17 s/config algorithm
  benchmark (prohibitive for the 8-config per-thread tensor-SH pool). Transforms run ~20 %
  slower than the tuned optimum; for long production runs, add `SHT_LOAD_SAVE_CFG` to
  `SHT_GAUSS` to recover tuned transforms with a cached, one-time init.
- `vilma.nml`: default 3-D field → Bagge 2021 (floored 19.5 via `visc_log10_min`).
- `vilma_drive`: coarse PROFILE timers (per-step read/update/write + `n_accept`/`n_solve`,
  one-time `build_remap` / `solid_earth_init`).

The batched primitives (`vilma_sht`/`vilma_tensor_sh` `clone_cfg_batched` / `synth_batch` /
`analysis_batch`) were fully reverted — do not re-attempt batching for CPU at these lmax.

## Remaining levers (if more speed is needed)

Ordered by leverage on `solid_earth_update`:

1. **Reduce `n_sub`** (the explicit-CFL multiplier). It is set by the stiffest Maxwell cell;
   the floor (19.5) already tames it. A per-element/per-layer local sub-cycling scheme
   (only the few stiff elements sub-step, the rest take the full Δt) would cut it further
   without changing physics. The trapezoidal/implicit scheme is *not* the answer — measured
   catastrophically slower here (step-doubling × Picard × heavier trap advance).
2. **Cheaper per-advance**: skip negligible-memory coefficients in the advance (like
   `begin_step`'s `skip_tol`).
3. **Tuned transforms for production**: `SHT_GAUSS + SHT_LOAD_SAVE_CFG` (~20 % per-step).
