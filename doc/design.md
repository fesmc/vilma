# VILMA — design

This document records the design decisions, the method comparison that led to
them, the validation plan, and the implementation pitfalls the GIA literature
flags. It is the reference for filling in the module stubs.

## 1. Goal and constraints

Replace **VILMA-v1** (the closed-source solid-Earth model coupled to CLIMBER-X)
with an open-source Fortran model that is:

- **state of the art** — full sea-level equation, rotational feedback, 3D
  (laterally varying) viscosity;
- **simple and fast** — no normal-mode root finding, no heavy 3D-FEM
  infrastructure;
- **3D-ready from the start**, validated first in 1D.

We reproduce **VILMA-v1's physics** (confirmed from Albrecht, Bagge & Klemann 2024,
*The Cryosphere* 18:4233, and Bagge et al. 2021): an **incompressible,
self-gravitating, Maxwell** viscoelastic sphere, solved by the
**spectral–finite-element, time-domain** method of **Martinec (2000)**, *GJI*
142:117, with rotational feedback (Martinec & Hagedoorn 2014) and a full
migrating-coastline sea-level equation. This is a new implementation
from the published literature.

## 2. Why this method (landscape)

| Family | Examples | 3D viscosity | Speed | Fit |
|---|---|---|---|---|
| Normal-mode + Love numbers + SLE | SELEN, TABOO, ALMA3 | no (1D only) | fast (1D) | dead end for a 3D model |
| **Time-domain spectral-FE** | **VILMA-v1** (Martinec 2000) | **yes** | **fast** | **chosen** |
| Full 3D FE/FV | CitcomSVE, ASPECT, Elmer | yes | slow, heavy | too complex |

The decisive structural fact: with the **explicit** Maxwell time scheme, the
stiffness solve has **no angular coupling**, so each spherical-harmonic degree
`l` decouples into an independent small banded radial system — embarrassingly
parallel and cheap. Lateral (3D) viscosity does **not** couple the harmonics in
the solve; it enters only through the explicitly-known **memory stress**,
evaluated pointwise on the Gauss-Legendre grid and transformed back to spectral
each step (pseudo-spectral). Therefore the **1D and 3D code paths are the same
solver** — 3D viscosity is just a spatially varying field instead of a constant.
That is what makes "3D-ready from day one, validated in 1D" cheap rather than a
rewrite.

## 3. Architecture (modules)

| Module | Role | Status |
|---|---|---|
| `vilma_precision` | working precision (= C double) | done |
| `vilma_constants` | physical / reference constants (benchmark conventions) | done |
| `vilma_params` | `vilma_param_class` + `vilma_par_load` (one `&vilma` nml group) | done + tested |
| `vilma_sht` | SHTns wrapper — the transform kernel | done + tested |
| `vilma_earth_structure` | radial layers (`build_earth`: named / custom) + optional 3D viscosity | done + tested |
| `vilma_radial_integrals` | Appendix C P1/P0 element integrals | done + tested |
| `vilma_band` | pivoted banded LU (dependency-free, re-entrant) | done + tested |
| `vilma_radial_fe` | per-degree saddle-point operator + banded-LU solve | done + tested |
| `vilma_viscoelastic` | Maxwell memory-stress time stepping (1-D, scheme-pluggable) + shared kernel | done + tested |
| `vilma_response` | surface-load → (uplift, geoid) operator: elastic + viscoelastic field driver | done + tested |
| `vilma_sle` | sea-level equation (ocean function, migration) | done + tested (elastic + VE) |
| `vilma_timestep` | adaptive-Δt controller (step-doubling on the Maxwell memory) | done + tested |
| `vilma_rotation` | rotational feedback / TPW (degree-2 Liouville + tidal VE channels, SLE-coupled) | done + tested (5a/5b/5c) |
| `vilma_coupling` | CLIMBER-X-compatible init/update/finalize API (adaptive-coupled) | done + tested |
| `vilma_io` | netCDF restart + step output (yelmo variable-table convention) | done + tested |
| `vilma_remap` | conservative lon-lat → Gauss remap (wraps fesm-utils `coords`) | done + tested |
| `vilma_drive` | standalone forced-run loop (`program vilma_main`; online remap + `i_eq`) | done + tested |
| `vilma` | umbrella re-export (incl. `vilma_param_class`, `vilma_run`) | done |

Convention in `vilma_sht`: fully-normalized real harmonics, no Condon-Shortley
phase, Gauss grid, phi-contiguous layout; spectral arrays hold `m >= 0`.

## 4. Coupling contract (from CLIMBER-X `src/geo/vilma.F90`)

- **in:** ice thickness `h_ice` [m]; **out:** relative sea level `rsl` [m] and
  bedrock `z_bed = z_bed_eq - rsl` [m].
- Host (CLIMBER-X) maps between its grid and the model's lat-lon Gauss grid with
  conservative (in) / bilinear (out) SCRIP weights; **all SH transforms stay
  inside this model**.
- Reference VILMA-v1 resolution: SH degree 170; SLE grid 1024×2048; radial FE 5 km
  (→420 km) / 10 km (→670 km) / 40–60 km (→CMB); explicit Δt 2.5 yr; coupling
  cadence `n_year_geo` (default 10 yr); viscosity floor ~1e19 Pa·s.

## 5. Validation ladder (each rung = a published benchmark)

1. **SH transform** round trip — `tests/test_sht.f90` (done).
2. **Love numbers** (elastic + viscoelastic) vs **Spada et al. (2011)**, *GJI*
   185:106, tests 2/1, 3/1, 5/1, 7/1 on Earth model **M3–L70–V01**; cross-check
   **ALMA3**.
3. **Deformation / geoid** for a disc load vs Spada (2011) tests 1/2–2/2.
4. **Sea-level equation** vs **Martinec et al. (2018)**, *GJI* 215:389, cases
   A→E (the migrating-coastline benchmark VILMA-v1 itself passed).
5. **Rotation** vs Spada (2011) test 3/2 (polar motion) — **DONE** (`test_rotation`,
   `test_rotation_sle`): degree-2 Liouville polar motion `|m(t)|` matches Table 14
   (Cw=0) to <1% for cap + disc at t=0–20 kyr; the SLE-coupled feedback (5c) matches
   Adhikari et al. (2016) eq. 8 to 1e-14 and conserves ocean mass. See §11.
6. **3D** (laterally-varying viscosity) vs **Weerdesteijn et al. (2023)**, *G-cubed*
   24:e2022GC010813, §5.2 — the ASPECT/Abaqus/TABOO low-viscosity-zone benchmark
   (disc load over a confined soft column) — **DONE** (`test_benchmark_lvz`): central
   uplift matches to 1–3% (axisymmetric). See §12.

## 6. Implementation pitfalls to design in (not patch later)

- **Degree-1 is frame-dependent** — do not hard-code `h1=l1=1, k1=0`. Work in the
  **CM frame** (Blewitt 2003). Geocenter motion *is* the degree-1 signal.
- **Rotation is degree-2** and uses **tidal** Love numbers; fix the fluid Love
  number `k_f` from the observed flattening (Mitrovica et al. 2005) to avoid the
  lithosphere-thickness paradox.
- **Lithosphere = η→∞ exactly**, never large-but-finite (a finite value injects a
  spurious slow relaxation mode that contaminates multi-cycle runs).
- **Pseudo-spectral SLE**: convolution in spectral space, ocean-function product
  on the spatial grid — this is what kills Gibbs ringing at coastlines.
- **Conserve mass, not volume**: apply `rho_ice/rho_water` exactly once each in
  the flotation criterion, the eustatic conversion, and the SLE; fix one density
  set (benchmark: ρ_ice = 931, ρ_water = 1000).
- The SLE is a **Fredholm equation of the second kind** → nested fixed-point
  iteration (~3 outer × 3 inner; Kendall, Mitrovica & Milne 2005).
- **Explicit-scheme stability**: enforce a viscosity floor; the time step is
  bounded by the shortest Maxwell time.

## 7. Build system

`configme` inserts a compiler/machine fragment into the template
`config/Makefile` (the pre-configme generator is kept in `config/legacy/`); shared
dependency wiring is in `config/common.mk`, source/rule lists in
`config/Makefile_vilma.mk`. Dependencies come from a `fesm-utils` symlink
(FFTW, SHTns, `fesmutils`) plus system netCDF (`nf-config`). SHTns was added to
fesm-utils' `build.py` as a first-class component (serial + OpenMP variants).

## 8. Key references

- Martinec (2000), *GJI* 142:117 — the spectral-FE method.
- Albrecht, Bagge & Klemann (2024), *The Cryosphere* 18:4233 — PISM–VILMA,
  confirms VILMA-v1 physics.
- Martinec & Hagedoorn (2014), *GJI* 199:1823 — time-domain rotational feedback.
- Kendall, Mitrovica & Milne (2005), *GJI* 161:679 — SLE with moving shorelines.
- Spada et al. (2011), *GJI* 185:106 — GIA benchmark.
- Martinec et al. (2018), *GJI* 215:389 — SLE benchmark.
- Weerdesteijn et al. (2023), *G-cubed* 24:e2022GC010813 — ASPECT GIA; the
  lateral-viscosity (low-viscosity-zone) benchmark vs Abaqus/TABOO (rung 6).
- Mitrovica et al. (2005), *GJI* 161:491 — revised rotation theory.
- Blewitt (2003), *JGR* 108:2103 — degree-1 / reference frames.

## 9. Radial solver formulation (rung 2) — working notes

Verified from Martinec's *Continuum Mechanics* lecture notes (Ch. 9), the Hanyk
PhD thesis, and the VEGA benchmark description (Martinec et al. 2018). Items
marked **(unverified)** await the paywalled Martinec (2000) PDF.

**Governing equations** (incompressible, self-gravitating, quasi-static):
- Incremental momentum balance `Div t^L + ρ₀ f^L + ρ^L f₀ + (pre-stress) = 0`,
  body force `f^L = −∇φ₁ + ∇u·∇φ₀` (Martinec Eq. 9.32, 9.152).
- Poisson via advected density: `∇²φ₁ = 4πG ρ₁`, `ρ₁ = −∇·(ρ₀u)` (Eq. 9.98).
- Incompressibility `Div u = 0`; isotropic stress Π is the Lagrange multiplier
  (K→∞), deviatoric part `τ^dev = 2μ ε^dev` (Maxwell in time, §D below).

**Spectral/radial:** spheroidal-only for 1D loading. Radial scalar set per degree
`{U, V, F}` nodal + pressure Π; tractions are natural BCs of the weak form.
**P1 ("tent") radial basis** (VEGA: 165 elements). Mesh = done (§3, `vilma_radial_fe`).

**Incompressibility discretization (B4) — RESOLVED.** Confirmed from the
Martinec (2000) PDF: **P1 (U,V,F nodal) / P0 (Π per element)**, `Div u = 0`
enforced weakly via the pressure block (eq 82). Assembled and validated — see
[formulation.md](formulation.md) and `vilma_radial_fe%build_dense_operator`.

**Boundary conditions — settled (whole-sphere mesh).** The fluid core IS meshed
(μ=0 region), so there is **no explicit CMB boundary condition** — free-slip
emerges. No explicit centre BC either: Martinec meshes through r=0 and the r²
weighting handles regularity (the singular `I⁷` term vanishes via `R₁=0`).
Surface: the (j+1)F(a) exterior match on the F-F diagonal + the load RHS.
Density-jump interfaces are natural conditions of the weak form. (The earlier
Wu & Peltier CMB-BC plan assumed an un-meshed core and is superseded.)

**Time scheme (rung 3) — DONE (1-D).** Explicit ω=1 Maxwell scheme (Martinec
2000 eqs 23-25), implemented in `vilma_viscoelastic%ve_degree`. Memory stress
`τ^{V,i} = (1−M)τ^{V,i-1} − 2μM ε^i`, `M = μΔt/η`; it enters the SAME elastic
operator as the dissipative RHS forcing `−∫τ^V:δε dV` (radial Gauss-2 + the
spectral double-dot over the four spheroidal tensor components, eqs 94/110).
Elastic layers (η→∞) freeze, fluid layers (μ=0) carry no memory. Stability
`Δt ≲ 2η_min/μ` ⇒ viscosity floor; VEGA Δt = 20 yr. The fixed operator is built
and LU-factored once (`vilma_band`, banded) and reused every step (~20 µs/solve).
Validated (`test_relax`): held load relaxes elastic→fluid, `t_relax ∝ η`.

**Love numbers (rung 2) — DONE.** `h_n = g₀ U(a)/φ^L`, `l_n = g₀ V(a)/φ^L`,
`k_n = −F(a)/φ^L − 1`, with `φ^L = 4πG a σ/(2n+1)` (`vilma_radial_fe%loading_love`).
The fluid (`h_n→−(2n+1)/3`, `k_n→−1`) and rigid (`→0`) limits are reproduced to
~1e-5 (`test_love`). `l_n` sign/normalization and the quantitative Spada (2011)
Test 2/1 match are the remaining calibration items.

## 10. Sea-level equation (rung 4) — working notes

**Response operator (`vilma_response`).** The SLE is built on an abstraction:
`response_apply(resp, σ_lm) → (u_lm, N_lm)` maps a spectral surface mass load
to surface uplift `u` and geoid height `N`, so the elastic and viscoelastic Earth
responses are swappable. **Geoid mapping: `N(a) = −F(a)/g`** (Bruns' formula; the
geopotential perturbation is `−F` since Martinec's `φ₁ = F → −φ^L` rigid). This
uses only `U` and `F`, *not* the unresolved horizontal Love number `l`, so the
SLE is not blocked by that open item. Three implementations:
- `elastic_response` — per-degree gains `gu(l)=U(a)`, `gn(l)=−F(a)/g` precomputed
  once from a unit-load solve;
- `null_response` — `u≡N≡0` (rigid, non-self-gravitating) → eustatic baseline;
- `ve_response` — viscoelastic field driver (below).

**SLE fixed point (`vilma_sle`) — DONE.** Pseudo-spectral (KMM 2005): the load →
response convolution is spectral, the ocean-function product `C·S` is pointwise on
the Gauss grid (kills coastline Gibbs ringing). For a fixed coastline,
`S = C·(N − u + Δφ)`, `N,u` = response to `L = ρ_i ΔI + ρ_w C·S`, and the spatial
constant `Δφ` is fixed *each iteration* by ocean-mass conservation
`ρ_w ∫C·S dA = −ρ_i ∫ΔI dA` ⇒ mass conserved to machine precision by construction.
Inner loop iterates `S`; outer loop migrates the coastline from `topo0 − S`.
Validated (`test_sle`): eustatic limit → uniform barystatic rise (mass resid
~4e-16); self-gravitating elastic M3 → mass resid 0, fixed point converges, real
spatial structure.

**Viscoelastic field driver (`ve_response`) — DONE.** The SLE needs a *field*
load with an **independent memory history per (l,m)** (each coefficient has its
own load time-series). Key insight: the response at a fixed time is **affine** in
the current load, `u_lm = gu(l)·σ_lm + drift_lm`, where `drift_lm` comes from the
*frozen* past-relaxation memory — so the SLE fixed point can call `apply()`
repeatedly without corrupting state. The step is bracketed: `begin_step` freezes
the drift (one memory-forcing solve per (l,m)), `commit_step` advances the Maxwell
memory with the converged load. Each complex (l,m) history is two real histories
(re/im) since the operator and `M = μΔt/η` are real. The per-element Maxwell
kernel (`strain_coeffs`, `ve_strain_constants`, `dissipative_rhs`,
`advance_memory`) is shared with the 1-D stepper — no duplicated algorithm.
Validated (`test_ve_response`): first-step gains == `elastic_response` exactly; a
held degree-2 load reproduces the 1-D `ve_degree` history to ~1e-13 relative.
End-to-end (`test_sle_ve`): a held 2 km ice cap → ocean mass conserved ~1e-16 per
step, eustatic mean held, ~16 m of viscoelastic relaxation over 500 yr.

**Degree-1 (geocenter) is now carried in the field driver — DONE.** Previously
skipped (`gu(1)=gn(1)=0`, `l < 2` guards) as a workaround for the *dense* j=1
operator hanging the stepper. The sparse KKT degree-1 operator (merged from
`degree1-sparse`) fixed that at the source — j=1 solves in ~2 GMRES iters in the
CM frame (wᵀd = 0 rigid-translation removal, Blewitt 2003). The `l < 2` guards
in `vilma_response` are now `l < 1`: degree 1 joins the unit-load loop in
`ve_response_init` (assembles `ops(1)`, real gains + nodal fields like any
degree), and carries memory drift through begin/apply/commit; only degree 0
(monopole geoid, no deformation, no operator) stays special. Validated
(`test_ve_response`): degree-1 first-step gains == `elastic_response` exactly
(0.0 diff), and the held (1,0) field response reproduces the 1-D `ve_degree` at
j=1 to ~1e-11 relative; `test_sle_ve` mass conservation unchanged (~6.6e-16).
The frame is inherited from the operator (same CM frame as `elastic_response`
and `ve_degree`), so the SLE's `N=−F/g` at degree 1 is self-consistent with the
displacement `u`.

**Open / next** (priority order for a fresh session):
0. ~~Fix the elastic low-degree Love-number bug.~~ **DONE.** Benchmark data
   in-repo (`data/benchmarks/love_M3-L70-V01/`; `test_benchmark_love`). Root cause:
   a transposed index in the self-gravity potential-gradient force in
   `build_dense_operator` (the U-F coupling used `I²_αβ` instead of `I²_βα`),
   which broke the U↔F symmetry the energy functional requires. Found by
   re-deriving the continuous gravity form (Martinec eq 65) term-by-term and
   cross-checking the shear block against the `vilma_viscoelastic` strain
   representation (eqs 85–88). Fix: `i2(ia,ib) → i2(ib,ia)`. Now elastic AND fluid
   M3-L70-V01 Love numbers match the benchmark to ~0.1% at every degree; the
   operator is exactly symmetric (`test_assembly`). Closes the disc offset (rungs
   2/3) at the source — a direct disc re-run to confirm <1% is a quick follow-up.
   (The earlier "lift the degree-1 skip in `ve_response`" item is now DONE on main
   — the geocenter degree-1 response is carried in the field driver.)
1. **`vilma_coupling` wiring — DONE.** The CLIMBER-X contract is wired: `ve_response`
   member, `solid_earth_init(se, p, sht, z_bed_eq, h_ice_ref)` / `solid_earth_update(se, h_ice, dt)`
   drives the adaptive stepper (SLE per Δt) across the coupling interval and returns
   `z_bed = z_bed_eq − rsl` (`test_coupling`).
2. **Grounded-ice flotation — DONE.** `ocean_function` (`vilma_sle`) is ocean only where
   `topo < 0` AND the ice floats (`ρ_i·I < −ρ_w·topo`); grounded ice over a subsided
   bed stays land. The load splits grounded ice `ρ_i·dI·(1−C)` + ocean water
   `ρ_w·C·rsl` + the subgrid sloping-coast term; floating ice contributes no load
   (buoyancy) and no melt-water source. Wired through coupling → stepper → SLE with
   absolute `ice` (flotation) and `d_ice` (load). `test_flotation`, `test_flotation_load`.
   Consistency follow-ups (this session): removed the dead `vilma_gravity` stub
   (self-gravity is solved in `vilma_radial_fe`), and the rotational polar-motion update
   now uses the SLE's converged end-of-interval load (via the stepper's `sigma_out`)
   instead of a hand re-derivation that dropped the subgrid term.
3. **Martinec et al. (2018) cases A–E quantitative match — DONE.** Case A (pure
   VE loading) in `test_benchmark_martinec` (`make check`). The four migrating-coast
   SLE cases C2/D3/E2/F1 in the standalone `test_benchmark_sle` (NOT in `make check`,
   it is a multi-minute lmax-128 run; bundled in `make check-slow`). All six fields
   (u, v_θ, v_φ, geoid, sea-surface, SLE) match the SBK reference curves (figs 10–13)
   within the inter-code scatter: D3/E2/F1 ~1–7%, C2 to ~7% on the coupled fields
   (~21–28% only on a small peak-normalized basin-region uplift/horizontal). The
   load/topo spec is the analytic giapy/Table-4 construction (caps L0–L2, exponential
   basins B0–B2, Heaviside/linear time, fixed vs migrating ocean), built in the test.
   Checked against Adhikari et al. (2016): its three validation experiments are
   elastic + qualitative (Farrell–Clark ocean-load fingerprint, a SELEN disc-melt run,
   a rotational-signature plot) with no vendored quantitative reference — subsumed by
   our Love + Martinec + Spada suite, so nothing added.
4. **Performance — DONE.** `begin_step` does 2 real solves per (l,m) per step
   (~O(nlm) solves), the cost driver at VILMA-v1 resolution. Four changes, all exact
   (results unchanged) except the threshold-controlled skip:
   - **Banded LU** (`vilma_band`) replaces the iterative solver (LIS GMRES+ILU) on the
     per-degree solve, and **LIS is removed entirely**. The operator is banded
     (half-bandwidth ~6) and the equilibrated system is effectively direct (1 GMRES
     iter), so a pivoted band LU — factor once at assemble, band solve per RHS — is
     far faster (~20 µs vs ~700 µs/solve) and cache-light (no ILU fill to evict).
     Pivoting is required (zero pressure (Π,Π) block). j=1 carries the dense KKT
     border (rigid-mode removal), so that one degree has ~full bandwidth and factors
     as a dense LU — still `vilma_band`, just wide. Dependency-free, **re-entrant** (so
     no serial-vs-OpenMP LIS variant to reconcile when linked into a larger host).
   - **Degree-grouped storage** of the per-(l,m) Maxwell memory/drift (slot `k`,
     `lm↔k` map) so the loop is contiguous and per-degree.
   - **Skip-negligible**: coefficients whose memory is < `skip_tol`×max are not
     solved (drift ≈ 0); ~2× for a localized cap.
   - **OpenMP** over the degree loop (`make openmp=1`, serial deps + `-fopenmp`);
     safe via the re-entrant band LU. ~5.4× at 8 threads.
   Net: `begin_step` at lmax 128 ≈ 7.9 s (orig LIS, lmax 64 extrapolated ~30 s) →
   ~58 ms (8 threads) — well over 100× single-thread from the band LU alone.
5. **Adaptive time stepping (§3c) — DONE.** The 2nd-order lever is the
   **trapezoidal** memory rule (Crank–Nicolson, order 2.00 vs FE order 1; the
   coupling iteration is its implicit solver) with the start-of-step load `σ_n`
   tracked. `vilma_timestep`'s `adaptive_stepper` crosses a coupling interval with the
   ice load linearly interpolated, choosing Δt by step-doubling on the Maxwell
   memory. Δt enters only as `Mk=(μ/η)Δt` (cheap rescale, no operator re-factor).
   Pays off on dynamic-range (ice-age) loads; ~1.6× wall there. See
   `doc/performance-assessment.md` §3c.
6. **Parameter type + nml + standalone driver — DONE.** One `vilma_param_class`
   loaded from a single `&vilma` namelist (`vilma_params`, yelmo `defaults_file`
   overlay; `vilma.nml` is the complete defaults; time fields in years→s).
   `build_earth(p)` (named built-in / custom layers). `solid_earth_init(se, p, …)` /
   `solid_earth_update(se, h_ice, dt)` distribute the knobs and run the adaptive controller per
   interval (fixed substeps removed). Restart persists the **full** integrator
   state (`dt_try` + `σ_n`) → bit-for-bit continuation. Umbrella module
   `vilma`; standalone `program vilma_main` (`vilma_drive`) runs a
   forced simulation from `&vilma` + an ice forcing already on the Gauss grid.
7. **Real ice forcing + lon-lat → Gauss remapping — DONE (§13).** `vilma_remap` wraps
   the fesm-utils `coords` library (great-circle conservative polygon clipping) to map
   a regular lon-lat field onto the Gauss grid; the standalone driver remaps each
   forcing slice ONLINE by default (`remap_input=.true.`; no preprocessing, raw data
   stays on disk) and an offline `vilma_remap` tool can pre-bake a Gauss forcing
   (identical engine). Reference / equilibration via `i_eq` (0: declare start slice as
   relaxed; 1: ice-free + paleotopo spin-up to LGM equilibrium under `dt_equil`).
   Validated on CLIMBER-X `geo_ice_tarasov_deglac.nc` (−26 ka→0). Depends on the
   fesm-utils `coords-dev` branch (`coords` module). See §13.

**Rung 6 (3D laterally-varying viscosity) — 6a/6b/6c DONE (§12).** Tensor-correct
pseudo-spectral memory advance (`vilma_tensor_sh`, general order `mmax≥0`), validated
against the Weerdesteijn et al. (2023) ASPECT/Abaqus/TABOO low-viscosity-zone
benchmark to 1–3%. 6c adds real lon-lat-r viscosity loading from netCDF
(`vilma_read_visc_3d`, node→element log10-mean bridge), per-thread-config element-loop
OpenMP, the implicit TRAP-3D advance, and an off-pole rotational-invariance test.
Rung 5 (rotation) is DONE — 5a/5b/5c (§11).

## 11. Rotational feedback / TPW (rung 5) — working notes

**Formulation (Spada et al. 2011 §2.1.1; time-domain à la Martinec & Hagedoorn 2014,
i.e. VILMA-v1).** Equatorial polar motion `m = m₁ + i m₂` from the GIA (quasi-static)
Liouville equation with the Chandler wobble neglected (eq. 7, justified since the
Chandler period ≪ GIA timescales):

    [1 − k^T(t)/k_s] ∗ m(t) = Ψ_L(t),   Ψ_L = I/(C−A),   I = [δ+k^L] ∗ I_rigid,

with `k^T,k^L` the degree-2 tidal/loading Love numbers and `k_s ≡ k^T_f` the secular
(fluid) tidal Love number (eq. 11). No explicit Ω — it is absorbed into `m = ω/Ω`
and `k_s`. Validated against the Cw-excluded column of Table 14 (the regime this
quasi-static form models; the Cw≠0 Chandler transient is deliberately out of scope).

**5a — tidal forcing path (`vilma_radial_fe`).** An external degree-j potential reuses
the SAME per-degree operator as a surface load; only the natural surface term differs.
In Martinec's φ₁ convention the external potential couples to F(a) with the SAME sign
as the load's own potential, `−(a/4πG)(2j+1)φ_t`, but with NO U-traction (a load
subsides; a tide-raising potential uplifts). `tidal_rhs` + `tidal_love`
(`k^T=−F/φ_t−1`). Validated (`test_tidal`) against the homogeneous-sphere Kelvin
tidal Love numbers: fluid `k^T_f=3/(2(n−1))`, `h^T_f=(2n+1)/(2(n−1))`; rigid →0;
degree-2 elastic `(3/2)/(1+μ̃)`, `μ̃=19μ/(2ρga)`.

**5b — Liouville solve (`vilma_rotation`), self-contained, degree-2 only.** Two compact
degree-2 *complex* viscoelastic channels (reusing the `vilma_viscoelastic` Maxwell
kernel; no normal modes, no convolution quadrature) carry the convolutions as memory:
a LOADING channel returns `(1+k^L)∗I_rigid = I(t)`; a TIDAL channel (`tidal_rhs`,
forced by the centrifugal potential ∝ `m`) returns `k^T∗m`. The feedback makes the
step ALGEBRAIC in `m` (the affine begin/apply/commit structure of the field driver):

    m_n = [ Ψ_L,n − dF_tidal/k_s ] / [ 1 − k^T_e/k_s ],

then both channels' memory is advanced. The rigid inertia `I₁₃+iI₂₃ =
−a⁴∫σ sinθcosθ e^{iφ}dΩ` is a DIRECT Gauss-grid quadrature of the load (3-D-ready —
no spherical-harmonic normalization assumption; verified by reproducing the paper's
published `G_cap/G_disc` to <0.5%). `k_s = k^T_f` from fluidizing the Maxwell mantle.

**Pitfall (designed-in): the lithosphere-thickness paradox (Mitrovica et al. 2005).**
The secular polar-motion slope is *pathologically* sensitive to `k_s` — a 0.5% change
in `k_s` moved the t=20 kyr `|m|` by ~6%. So `k_s` must be the consistent fluid limit:
fluidize ONLY `RHEOL_MAXWELL` layers (the viscous mantle); keeping the lithosphere
elastic is essential (fluidizing it by mistake inflated `k^T_f` 0.967→0.975 and
under-drove the late-time motion ~11%). For deep-time runs `k_s` is a parameter,
overridable to the observed-flattening value (§5c). Validated (`test_rotation`): cap + disc
`|m(t)|` match Table 14 (Cw=0) to <1% at t=0–20 kyr; `k^T_e=0.303`, `k_s=0.967`.

**5c — feedback into the SLE.** The centrifugal potential of `m` perturbs the sea
surface and deforms the solid, adding a degree-2 contribution to relative sea level
`s_rot = N_rot − u_rot` with (Adhikari et al. 2016, eq. 8) `N_rot = (1+k^T)Λ/g`,
`u_rot = h^T Λ/g`, `Λ = −Ω²a² sinθcosθ(m₁cosφ+m₂sinφ)`. The VE `(1+k^T)`,`h^T` reuse
the 5b tidal channel (which now exposes the uplift readout); the rotational fields
reuse the `m`-forced channel exactly. `s_rot` enters the SLE geometry (`Sraw`) but NOT
the surface mass load — the rotational potential forces the Earth through the tidal
channel, not as a load — and `Δφ` is recomputed so ocean mass stays conserved. The
rotation ↔ SLE coupling is a fixed point (`m` responds to the ice+ocean load);
`vilma_rotation` is split into begin_step / solve_m / s_rot / commit (affine, no memory
advance until commit) so it can be iterated. It steps with the solid Earth:
`sle_solve(…, rot, rot_dt)` runs `rotation_trial` (from the step-start state, the
explicit-FE channels sub-cycled to `dt_fe_max` with the load held) in every inner
iteration, adds its `s_rot`, and commits with `rotation_close_step` where the response
memory advances; `stepper_advance` saves/restores it alongside the response. (It used
to be applied once per coupling interval, as a predictor with `s_rot` held from the
entering `m`, which lagged the rotational sea level by one interval.) The channels
relax on the response's radial viscosity (`rotation_set_rate`). `k_s` exposes two values: the model fluid limit `k_s_fluid`
(default, reproduces Spada) and the observed-flattening closed form
`k_s_flat = 3G(C−A)/(a⁵Ω²) = 0.943` (Adhikari/Mitrovica — the recommended deep-time
value). Validated (`test_rotation_sle`, elastic): hook-off is bit-for-bit the
no-rotation SLE; `s_rot` matches Adhikari eq. 8 pointwise to 1.9e-14; ocean-mass
residual 3e-16; the fixed point converges (ocean feedback on `m` ~0.7%);
`|s_rot| ≈ 1.9 m`. End-to-end (`test_coupling`): an off-axis cap drives `|m| = 0.30°`,
mass conserved, the bed shifts up to 6.7 m vs rotation-off.

## 12. 3D laterally-varying viscosity (rung 6) — working notes

> **Superseded in one respect** (2026-09): these notes describe a spheroidal-only
> memory, λ∈{1,2,5,6}. Lateral viscosity also drives toroidal flow, which the
> model now carries (λ3,4 and the displacement W) — see `design-toroidal.md`.

**The one structural change.** Lateral viscosity makes the Maxwell factor
`M = μΔt/η` a *field* `M(θ,φ)`, so the memory update `τ⁺ = (1−M)τ − 2μM·ε` has a
pointwise lateral product that *couples* harmonics. Everything else is untouched:
the per-degree elastic solve and the dissipative RHS `−∫τ^V:δε` are linear in the
memory (no lateral product), so they stay the exact 1-D code path. Only **η** varies
laterally — **μ and ρ stay radial** — which is why the operator/LU factorisation and
the whole speed argument survive. This is the "3D-ready from day one" payoff (§2).

**The pitfall (and the false start).** The memory/strain are second-order **tensor**
fields; the four kept components λ∈{1,2,5,6} are *tensor* spherical-harmonic
coefficients (degree-dependent norms `[1, Jr/2, 2Jr², 2(Jr−1)Jr-ish]`, Jr=l(l+1) —
λ=1 scalar, λ=2 vector/gradient, λ=5,6 rank-2), **not** scalar Y_lm. The first 6a cut
scalar-synthesised each λ independently and multiplied by M(θ,φ); that is wrong for a
laterally-varying M (it neither reconstructs the physical tensor nor captures the
cross-λ/(l,m) coupling a scalar×tensor product makes). A uniform M hides it (the
factor pulls out, the round trip cancels), so the uniform test passed while the real
LVZ benchmark blew up to −8 m (Δt- and de-aliasing-independent — a formulation error).
Lesson: validate the lateral coupling with a *non-uniform* M, not just uniform.

**6a — tensor-correct pseudo-spectral memory advance — DONE** (`vilma_tensor_sh`,
`advance_memory_3d`; commits after 2026-06-25). The update `τ⁺=(1−M)τ−2μM·ε` is
pointwise in PHYSICAL space (Martinec 2000 eq 102), so per radial element and per
radial shape-coeff (A,B,C) the memory and strain TENSORS are reconstructed on the
Gauss grid via their six **dyadic components** (eqs 90/91, B10/B11), advanced pointwise
with the lateral field M(θ,φ), and projected back. `vilma_viscoelastic` stays grid-unaware;
`vilma_response` owns the transform. Set via `ve%enable_lateral_visc(sht, pert_elem)` —
per-element `(nphi,nlat,ne)` log10 η perturbation, `η_eff=η·10^pert` ⇒
`MkPerDt3=(μ/η)·10^(−pert)`; elastic/fluid elements keep `MkPerDt=0` so the
**lithosphere stays exactly elastic** and is skipped. FE scheme only; implicit TRAP-3D
is a guarded follow-up.

`vilma_tensor_sh` maps the spheroidal tensor-harmonic coefficients (λ=1,2,5,6) ↔ the six
dyadic grid fields (rr, rθ, rφ, θθ, θφ, φφ) at GENERAL order. **Synthesis is exact via
grid identities** (no recurrence, no re-analysis): rr/rθ,rφ/trace use SHTns scalar +
spheroidal-vector synth; the spin-2 channel (G,H — the only piece SHTns has no routine
for) uses ∂_φ = "multiply coeffs by im" and ∂_θθ via ∇₁²Y=−l(l+1)Y, giving
`Sg=∇₁²f−2cotθ·g_θ−2(1/sinθ)∂_φ g_φ`, `Sh=(1/sinθ)∂_φ g_θ−cotθ·g_φ`. **Analysis**:
channels 1/2/5 invert through SHTns's own scalar/vector analyses (the −l(l+1) factor
for the trace); the spin-2 channel uses the adjoint of its synthesis (NB: the vector
synth's ∫dΩ adjoint is l(l+1)·`spat_to_SHsphtor`, its *inverse* differing by the
spheroidal norm l(l+1); and the θφ projection carries `e_θφ:e_θφ=½`), normalised by a
per-degree factor calibrated once at init. Validated: `test_tensor_sh` (mmax>0) — round
trip = identity (3e-14) and the physical `∫τ:τ` (dyadic weights `[1,1,1,½,½,½]`) equals
`Σ_λ norm_λ|τ^λ|²` with the B13 norms (7e-16). `test_response_3d` (now NON-axisymmetric,
mmax=lmax) — a uniform M reduces to the 1-D advance per (l,m) (memory ~5e-13, uplift
~4e-14), masking the inert λ=6 degree-1 null (no Z⁶ there). Serial over elements: the
dyadic transforms call SHTns, which is not safe for concurrent calls on one config —
the element-loop OpenMP was reverted; per-thread configs would restore it (perf follow-up).

**6b — LVZ benchmark vs Weerdesteijn et al. (2023) — DONE** (`test_benchmark_lvz`,
standalone; not in `make check` — lmax 512 ≈ 67 s). Their §5.2 short-timescale case
(Earth M3-L70-V01 = `build_M3L70V01`; axisymmetric ice disc R=100 km, H 0→100 m over
100 yr then held to 200 yr; cylindrical LVZ under the load, radius 100 km, depth
70–170 km, η=1e19 vs 1e21). Driven through `ve_response` directly (pure ice-load
deformation, no SLE). Load given by the **exact spherical-cap SH coefficients** (a
sampled step's quadrature error makes the mass — hence the uplift — nlat-dependent).
At lmax 512 (`nlat=3·lmax` de-aliases the spin-2 G channel; the disc-edge Gibbs
oscillation settles by ~512): load-center uplift at 200 yr = **−0.731 m homogeneous
(ref −0.75, 2.6%)** and **−1.218 m with LVZ (ref −1.23, 1.0%)**, amplification
**1.67 (ref 1.64)**; the LVZ trajectory saturates (the scalar scaffold ran away).
ASPECT/Abaqus themselves differ 1–3% and ASPECT(no self-grav) vs TABOO(self-grav)
0.28% at the load center, so self-gravity/sphericity are sub-few-% here.

**6c — real-field loading, per-thread OpenMP, TRAP-3D, rotational-invariance — DONE.**
- **Real 3D (lon-lat-r) viscosity from netCDF** — `vilma_read_visc_3d` (`vilma_io`) reads a
  log10(η) field (e.g. Pan et al. 2022: `eta(r,lat,lon)`, absolute log10) and
  interpolates onto the Gauss grid × FE nodes (bilinear-in-log10 horizontally with
  periodic longitude; linear-in-radius, clamped) into `earth%visc_3d` as ABSOLUTE
  log10(η). `ve%enable_lateral_visc_from_nodes` bridges node→element by the log10-mean
  of the two bracketing nodes and forms the per-element perturbation against the
  element's radial reference η. Only genuinely Maxwell elements receive the field: the
  advance and the bridge skip elements with `MkPerDt==0`, and `ve_init` now sets that
  rate to exactly 0 for elastic/fluid layers BY RHEOLOGY (previously the elastic
  lithosphere kept η=huge ⇒ MkPerDt=μ/huge≈3e-298, nonzero), so the elastic lithosphere
  is left exactly elastic rather than overwritten by the loaded field. (`test_visc_load`:
  synthetic round trip ~2e-13; pan2022 loads finite.)
- **Element-loop OpenMP with per-thread SHTns configs** — `tensor_sh` holds a config
  pool (one clone per `omp_get_max_threads`, built serially since FFTW planning is not
  thread-safe); each element's transforms run on the calling thread's `thread_cfg`. The
  four `sht_grid` transforms gained an optional `cfg` override. (openmp=1 LVZ reproduces
  the serial −0.7306/−1.2180 m exactly.)
- **TRAP-3D** — `trapezoid_advance_all` dispatches to `advance_memory_3d_trap` when
  `lat_visc`; per shape-coeff, τ_n/ε_n/ε_{n+1} are reconstructed on the grid and advanced
  by `τ⁺=[(1−M/2)τ_n−μM(ε_n+ε_{n+1})]/(1+M/2)` with the lateral M-field. (`test_response_3d`
  TRAP case reduces to the 1-D trapezoid ~1e-12.) This unblocks the adaptive controller
  and SLE co-convergence with 3D viscosity.
- **Off-pole rotational-invariance** — `test_rotinv`: a cap+LVZ on an off-pole axis
  (all m) matches the on-pole (m=0) cap-centre uplift to ~3e-4 — a strong end-to-end
  m>0 check, run through TRAP-3D. In `make check` at lmax 16; full-res standalone via arg.

**Open / next.** Wiring `visc_3d` loading into the params/driver path (lateral
viscosity from a file in a forced run).

## 13. Real ice forcing — lon-lat → Gauss remapping + driver (working notes)

Goal: drive the standalone model with real datasets (CLIMBER-X
`geo_ice_tarasov_deglac.nc`: Tarasov PMIP4 ice+bed, 1°×0.5° lon-lat, −26 ka→0 at
100-yr slices) with **no preprocessing** — raw data on disk, remapped internally so
the remap is always consistent with the running code.

**`vilma_remap` — conservative lon-lat → Gauss.** Wraps the fesm-utils `coords` library
(`map_init_conservative` great-circle polygon clipping; `map_field` area-weighted
mean). `ll2gauss_init(map, sht, lon, lat)` builds the map once (target = the SHTns
Gauss grid: longitudes `sht%lon`, latitudes `90−sht%colat`, built south-first for
coords then flipped to the SHTns north-first row order); `ll2gauss_apply` remaps a field per
slice. `coords` derives target cell boundaries from axis midpoints (not the Gauss
weights), so `apply(conserve_mass=.true.)` rescales a mass-bearing field (ice) by one
global factor so its SHTns quadrature integral equals the source area-integral exactly
(the whole sphere is 4π sr — no planet radius needed). Geometry (bed) is remapped as-is.
`test_remap`: constant preserved ~3e-15, latitude orientation, zonal accuracy, mass
exact. **Depends on the fesm-utils `coords-dev` branch** (the `coords` module); point
the `fesm-utils` symlink at a `coords-dev` checkout and build its utils lib.

**Driver (`vilma_drive`).** General forced-transient loop, fields read per
(file, var, time-index) so experiments mix and match. `remap_input=.true.` (default)
remaps each forcing slice online; `.false.` is the legacy Gauss-grid path. The offline
`vilma_remap` tool pre-bakes a Gauss forcing with the same engine (identical
results) for workflows that prefer it. `time_init/time_end` clip the record.
`vilma.x cfg.nml [defaults.nml]` allows a sparse run config over `vilma.nml`.

**Reference / equilibration (`i_eq`).** The deglac record starts already glaciated
(26 ka), so there is no ice-free pre-glacial slice. Two references:
- `i_eq=0` — declare the first in-window slice as the relaxed reference (`z_bed_eq`=
  bed[k0], `h_ice_ref`=ice[k0], memory 0); transient load = ice change vs that slice.
- `i_eq=1` (default) — ice-free reference; a **paleotopo fixed point** (the F1 pattern,
  §10) finds the ice-free relaxed bed whose viscoelastic equilibrium under the start
  ice reproduces the data bed[k0], holding the load `dt_equil` per pass; this leaves
  the model at LGM equilibrium (bed AND viscous memory) before deglaciating with the
  absolute ice as load. Converges in ~3 passes (mean residual ~170→5→<1 m on the
  Tarasov LGM), and the residual is the (i_eq=1 vs 0) consistency check — small ⇒ the
  data bed is a good equilibrium ⇒ `i_eq=0` suffices and the spin-up can be skipped.

Both modes agree on the physical bed at the LGM; they differ in the initial memory
state, which changes the deglacial rebound — the thing `i_eq=1` gets right and the
comparison quantifies.

## 14. Forced-run roadmap — IMPLEMENTED

Captured at the close of the real-forcing session (2026-06-25); implemented on
branch `s14-forced`. Five items, in rough priority order. As built (differences
from the original plan noted):

- **(a) years output** — `vilma_io` writes the `time` axis and `dt_try` in YEARS;
  restart read converts back to SI. Done.
- **(b) BSL** — `solid_earth%bsl` (scalar, vs reference), surfaced through `vilma_io`
  as a per-step time series. The init ocean function is seeded from reference
  flotation so the seed step yields a physical value (`ocean_function` made public).
- **(c) reference frame** — adopted CLIMBER-X `i_equilibrium` semantics exactly:
  `i_eq` 0=start-slice, 1=present-day reference (default), 2=eq-file, 3=ref+rsl.
  `z_bed_ref_file`/`h_ice_ref_file` (RTopo) are remapped online with a per-file
  conservative map. With `i_eq=1` rsl≈0 at PD by construction (no `rsl_anom` param
  needed). The ice-free paleotopo fixed point is retained but non-default, gated by
  `dt_equil>0` (default 0); it supersedes the i_eq bed when on.
- **(d) 3D viscosity** — `l_visc_3d`/`visc_3d_file` wired through `solid_earth_init`;
  `vilma_read_visc_3d` moved vilma_io→vilma_earth_structure (avoids a circular dep); a
  `visc_log10_min/max` clamp imposes the viscosity floor after read (so the base
  Bagge field reproduces the `_min_19.5` variant). Bagge field is in `isostasy_data`
  (`earth_structure/viscosity/bagge2021.nc`, reformatted to the Pan layout).
- **(e) uncertainty** — `f_visc_sd` (units of σ) with a RELATIVE σ: read from
  `name_visc_sd` if set, else `f_visc_rel·log10η` (default 0.1). Clamped as in (d).

Original notes (for reference) follow.

**(a) Output time units → years (all files).** The driver/`vilma_io` currently write
the `time` axis (and restart time fields) in SECONDS; `vilma_params` already converts
years→s on load. Make YEARS the default OUTPUT unit too — the inverse conversion in
`vilma_write_step` / the restart writer — so step output, time series, and restarts are
all in years (SI stays internal). Low-risk, do first.

**(b) Barystatic sea level (BSL) diagnostic.** Add BSL as an output (scalar time
series, and optionally the uniform field). It is already in hand: the SLE computes the
eustatic offset Δφ (`res%esl`) and the melt source `ice_int`; the barystatic value is
`bary = −(ρ_i/ρ_w)∫ΔI dΩ / ∫C dΩ` (see `test_benchmark_sle`). Surface it through
`vilma_coupling`/`vilma_io` per step. Useful both as a physical diagnostic and to audit the
eustatic budget (see (c)).

**(c) PD sea level "still negative" — reference frame, not a bug (verify).** With
`i_eq=1` the reference is the ICE-FREE relaxed state, so `rsl` is measured relative to
an ice-free world. At PD the SH ocean is therefore negative by ≈ the sea-level
equivalent of present-day ice still locked up (Antarctica+Greenland, ≈ −60–65 m) plus
the local Antarctic-load depression — EXPECTED, not an error. To get "≈0 at today"
output **sea level as an anomaly relative to a chosen epoch** (PD/final slice, or LGM):
add an output `rsl_anom = rsl(t) − rsl(t_ref)`. Then verify the LGM→PD change matches
the barystatic budget from (b) (~+120–130 m global rise) — that is the real
correctness check. Decide the canonical output reference frame and document it.

**(d) Run with a loaded 3D viscosity field (driver wiring).** `vilma_read_visc_3d`
(rung 6c) already loads a lon-lat-r log10(η) field onto the Gauss grid × FE nodes;
it is NOT yet wired into `vilma_params`/`vilma_drive`. Add `l_visc_3d`, `visc_3d_file` (+ the
var/axis names) to `&vilma` and load it in `solid_earth_init`. Target two fields:
  - the current Pan et al. (2022) field (already used in `test_visc_load`);
  - the CLIMBER-X production field `~/models/climber-x/input/vilma/visc3d_Bagge2021*.nc`
    (Bagge et al. 2021; var `lgvisc(radius,lat,lon)`, 512×256×164, log10 dex). Its
    512×256 lon-lat maps cleanly onto our lmax-128 Gauss grid (512×258).

**(e) Viscosity-uncertainty sampling — `f_visc_sd`, with a RELATIVE sd.** Mirror the
CLIMBER-X VILMA-v1 scheme (`src/geo/vilma.F90`): perturb in log10 space,
`log10 η → log10 η + f_visc_sd·σ`, clamped to `[visc_log10_min, visc_log10_max]`.
`f_visc_sd` is in units of standard deviation (0 = mean field, +1 = +1σ). CLIMBER-X
uses a CONSTANT `sigma_log10_visc` (default 0.5 dex) applied uniformly — the user's
view (and ours) is that a constant floor is the wrong default. **Improvement:** define
σ as RELATIVE to the field, spatially variable — a parameter giving σ as a fraction of
the mean log10 field, `σ(x) = f_visc_rel · log10 η_mean(x)` (in log10 space), so the
perturbation tracks the viscosity structure. If the dataset ships its own uncertainty
field, prefer that; otherwise fall back to the relative `f_visc_rel`. Params:
`f_visc_sd`, `f_visc_rel`, `visc_log10_min/max`. Enables ensemble GIA runs sampling
viscosity uncertainty.
