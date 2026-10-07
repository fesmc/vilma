# Benchmark reference data — provenance

Reference datasets for validating VILMA2 against the published GIA
community benchmark. None of this is VILMA2 output; it is external
reference data, vendored here so the validation tests are self-contained.

## `love_M3-L70-V01/mod_M3-L70-V01`

Normal-mode loading Love numbers for the **M3-L70-V01** earth model
(incompressible, self-gravitating, layered Maxwell; 70 km elastic
lithosphere, three mantle layers, inviscid core), degrees 1–256.

- **Source:** the `giapy` package (S. B. Kachuck), file
  `giapy/data/earth/mod_M3-L70-V01`. https://github.com/skachuck/giapy
- **License:** MIT (giapy) — redistribution with attribution is permitted.
- **Underlying benchmark:** the Charles University GIA Benchmark
  (https://geofjv.troja.mff.cuni.cz/GIABenchmark), i.e. the model defined in
  Spada et al. (2011), *A benchmark study for glacial isostatic adjustment
  codes*, Geophys. J. Int. 185, 106–132.
- **How it was generated (not by us):** these are *normal-mode*
  (analytic / semi-analytic) viscoelastic Love numbers — for an incompressible
  layered Maxwell sphere the response is a finite sum of decaying exponentials
  (here 9 modes per degree). They are the output of a normal-mode code such as
  TABOO or ALMA (Spada's codes), not a time-domain FE code. VILMA2
  computes the *same physics* by a different route (per-degree FE solve for the
  elastic/fluid limits, time-domain memory-stress integration for the transient),
  so this table is an independent cross-check, and TABOO
  (https://github.com/danielemelini/TABOO, GPLv3) can regenerate it externally
  if needed.

### File format (per the giapy loader)

- Lines 1–5: the earth model — `layer  r[m]  rho[kg/m^3]  mu[Pa]  eta[Pa·s]`
  (η = 1e44 ≈ ∞ marks the elastic lithosphere; μ = η = 0 marks the inviscid core).
- Then, for each degree `n`:
  - header `n  nmodes  k_e  h_e  l_e`  — **elastic** (instantaneous) Love numbers,
  - `nmodes` lines of normal-mode residues + relaxation,
  - one line `... h_f l_f k_f` — the **fluid** (t→∞) Love numbers.
  - (column order in the elastic/fluid lines is k, h, l after the leading index.)

## `disc_spada2011/{u,n,dudt}_disc.txt`

Spatial response of the **disc load** test of the Spada et al. (2011) benchmark:
vertical displacement `u`, geoid `n`, and uplift rate `dudt`, for the
**M3-L70-V01** earth model. Load: a 10°-radius disc of 1000 m ice
(ρ_ice = 931 kg m⁻³) on elevated bedrock (no ocean load).

- **Source:** the `isostasy_data` repository (J. Jereczek),
  `model_outputs/Spada-2011/`. https://github.com/JanJereczek/isostasy_data
- **Underlying benchmark:** Spada et al. (2011), *A benchmark study for glacial
  isostatic adjustment codes*, Geophys. J. Int. 185, 106–132 (the Charles
  University GIA Benchmark). Same earth model as the Love table above.
- **Grid:** 201 rows, colatitude θ = 0:0.1:20°; 6 columns, times
  t = [0, 1, 2, 5, 10, 100] kyr. **Column 1 (t=0) is the elastic response.**
- **θ ordering.** All files in this repository are stored in ascending θ
  order (row 1 is θ=0°, row 201 is θ=20°) so they can be paired with the θ
  grid directly. The upstream `n_*` (geoid) files in `isostasy_data` are
  stored in reversed θ order, which the loader
  (`dataloaders.jl`) undoes with `reverse!(X, dims=1)` for the `n_` cases
  only; `u_*` and `dudt_*` are not reversed there. The `n_disc.txt` shipped
  here has been reversed once at ingest so that all three files use the same
  ascending-θ convention.
- **Sampling — important when quoting errors against these files.** The stored
  curves are **piecewise linear on a 1° θ grid**, linearly interpolated onto the
  0.1° grid they are distributed on. Only every tenth row is data; the nine
  between are straight-line fill. The signature is unambiguous: all six time
  columns of `u_disc.txt` consist of exactly 40 linear segments whose slope
  breaks fall at the *same* 35 colatitudes, at integer degrees (0.9/1.0,
  1.9/2.0, 6.0/6.1, 7.0/7.1, …) — a smooth physical field cannot have slope
  discontinuities at the same colatitudes at every epoch. `n_disc.txt` has the
  same structure with slightly different breakpoints. The files are stored at
  full double precision, which disguises this.

  Consequence: the reference's own roughness (RMS second difference along θ) is
  **2.4–2.8× the model's** at every epoch, so a point-by-point comparison on the
  0.1° grid charges the model for the reference's interpolation error. Quoted
  peak-normalized errors against these files are therefore an **upper bound** on
  the model's true error. A comparison restricted to the ~21 reference nodes
  would be the honest measure. (Whether the 1° data was digitized from the
  published figures or taken from a table is not recorded upstream.)
- **Degree-1 / geoid frame:** the displacement is in the CE-like gauge (h₁≈0,
  geocenter), the geoid in the CM frame (N₁=0). VILMA2 reproduces both: u to
  ~1% near-field, n to ~1% once the degree-1 geoid is referenced to CM (N₁=0; see
  vilma_response). The far-field (θ≳12°) forebulge in `u` is small-amplitude and
  shows larger relative differences (low-degree-truncation sensitive).

## `sle_martinec2018/*_SBK.dat`

Reference spatial response curves for the sea-level-equation benchmark,
cases A, C2, D3, E2, F1 (figs 10–13), at truncation j256.

- **Source:** `giapy`, `tests/sbk_benchmark_data.tar.gz` (MIT).
- **Underlying benchmark:** the same Charles University GIA Benchmark; the SLE
  intercomparison published as Martinec et al. (2018), *A benchmark study of
  numerical implementations of the sea level equation in GIA modelling*,
  Geophys. J. Int. 215, 389–414.
- **Two file formats** (the case-A loading file differs from the SLE-case files):
  - **`A_fig10_SBK.dat`** (case A only): 4 columns, header
    `# Colat(deg) Uplift(m) Horizontal(m) Geoid(m)`; rows are colatitude
    180°→0° at 0.25° spacing (721 rows). This case has NO ocean (B0=0): a pure
    viscoelastic cap-loading response, with degrees n=0,1 dropped (giapy jmin=2).
  - **`{C2,D3,E2,F1}_fig{10..13}_SBK.dat`** (the SLE cases): 7 data columns —
    `col1` longitude/colatitude [deg]; `col2` vertical displacement [m];
    `col3`,`col4` θ- and φ-horizontal displacement [m]; `col5` gravitational
    potential increment [SI, = −geoid·9.815]; `col6` sea-surface variation w.r.t.
    h_UF [m]; `col7` sea-level-equation result [m]. Profiles are along circles of
    constant lon or lat (figs 10–13 use lon=75, lat=25, lon=±/25, lat=35/100).
- **Benchmark spec** — the authority is **Martinec et al. (2018) Table 4**
  (`doc/refs/Martinec2018.pdf`), NOT giapy's `tests/sle_test.py` (whose working
  case labels D1/D2/E1 and basin numbering differ from the paper). The paper
  defines 5 cases (A–E), 3 ice models, 3 basins:
  - **Ice** (sqrt cap h(ψ)=h0·√[(cosψ−cosα)/(1−cosα)], α=10°): **L0** pole,1500 m;
    **L1** colat25/lon75, 1500 m; **L2** colat25/lon75, 500 m.
  - **Basin** (ζ⁽⁰⁾ = bmax − b0·exp(−ψ²/2σ²), σ=26°): **B0** none; **B1**
    colat100/lon320, 760/1200 (shallow, far from ice); **B2** colat35/lon25,
    3800/6000 (2200 m deep at centre, near ice).
  - **Time**: **T0** Heaviside (full load held, terminal time tf = t1 = 10 kyr);
    **T1** linear growth over 10 kyr then held to t2 = 15 kyr.
  - **Ocean**: **SLE0** none (case A); **SLE1** fixed ocean geometry O⁽⁰⁾ (held);
    **SLE2** time-varying (migrating coastline).
  - Densities ρ_ice/ρ_w = 931/1000 kg m⁻³, g₀ = 9.815 m s⁻².
- **SBK file ↔ paper case** (the SBK letter = paper letter + 1; verified against
  the data — cap from the fig10 uplift peak, basin from the fig12 bump, ocean mode
  from the far-field esl, since a migrating coastline drains a fixed-ocean case ~2×
  too far):
  | SBK | paper | SLE | ice | time | basin |
  |-----|-------|-----|-----|------|-------|
  | A   | A | SLE0 | L0 | T0 | B0 (none) |
  | C2  | B | SLE1 fixed | L1 (1500 m) | T0 Heaviside | B1 (760/1200, shallow) |
  | D3  | C | SLE1 fixed | L2 (500 m) | T1 | B2 (3800/6000, deep) |
  | E2  | D | SLE2 migrate | L2 (500 m) | T1 | B2 (3800/6000, deep) |
  | F1  | E | SLE2 migrate | L2 (500 m) | T1 | B2* (paleo-SET prescribed at tf) |

  Validated by `tests/test_benchmark_sle.f90` (standalone): all four reproduce the
  SBK curves to within the inter-code scatter (esl matches to ≲5 m; D3 to ~0.2 m).

## `rotation_spada2011/` — polar motion (Test 3/2) + constants

Reference for rung 5 (rotational feedback / true polar wander), transcribed
directly from the **Spada et al. (2011)** paper PDF (`doc/refs/Spada2011.pdf`).

- **`polar_motion_test3-2.txt`** — Table 14 (p. 129): modulus of polar motion
  `|m|` [deg] and rate `|ṁ|` [deg/Ma] at t = 0,1,2,5,10,20 kyr, for the three
  axisymmetric loads (cap/disc/point, Table 4), model M3-L70-V01, Heaviside
  history, load centroid θc = 25°, λc = 75°. Two contributors — **Vb**
  (Mathematica matrix-propagator, Sabadini & Vermeersen) and **Gs** (TABOO,
  Spada) — both *normal-mode* codes; **VILMA1 (Vk) did not provide m for this
  test** (only LOD via 1+k^L), so this validates a time-domain solver the same
  cross-method way the giapy table validates rung 2. Two cases: Chandler wobble
  *included* (`incl`, Cw≠0) and *excluded* (`excl`, Cw=0). A quasi-static
  time-domain method (Liouville term i·ṁ/σ_r dropped, eq. 7) targets the **`excl`
  column** — it cannot reproduce the Cw≠0 ~1 kyr Chandler transient.
- **`constants_table2.txt`** — Table 2 (p. 109): moments of inertia C, A,
  angular velocity Ω, radius a, etc. — the rotation constants the earlier rungs
  did not need (C−A and Ω enter the excitation eqs 20/31; C the LOD eq 33).

**Source-notation gotcha (decoded):** Table 14 prints a compact power-of-ten form
— a trailing **subscript** digit d = ×10⁻ᵈ, a **superscript** digit d = ×10⁺ᵈ
(e.g. `0.132₁` = 0.0132, `0.681⁵` = 6.81e4, `Gσ 0.541₄` = 0.541e-4). Verified via
the Cw-included |ṁ|(t=0)=×10⁵ Chandler spike, |m|↔|ṁ| finite-difference
consistency, and Gσ Arg = 105°W (eq. 31). The vendored `.txt` is **fully expanded
to decimal** so downstream code never re-parses the notation.

- **License / regeneration:** the values are published benchmark reference data.
  **TABOO** (https://github.com/danielemelini/TABOO, GPLv3, Spada's "Gs" code)
  can regenerate Test 3/2 externally as an independent cross-check; not required.
  This table may be contributed upstream to `isostasy_data/model_outputs/Spada-2011/`
  later (it has the spatial Spada-2011 cases but not the rotational ones).
