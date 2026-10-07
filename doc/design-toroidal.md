# VILMA2 — the toroidal degree of freedom (design)

**Status: implemented** (branch `toroidal`, 2026-09). §0 records what was
built and where it departs from the plan; §1–§6 are the original scoping
document, kept as written except where marked. It complements `design.md` (the
Martinec FE method) and `formulation.md` (the equations as implemented).

## 0. As implemented

| commit | content |
|---|---|
| `cff40ba` | `vilma_sht`: toroidal synthesis `e_r×∇₁T`, one-call sph+tor synthesis, sph+tor analysis |
| `de90472` | `vilma_tensor_sh`: Z³, Z⁴ as channels 5, 6 (λ3, λ4 appended) |
| `8d955f0` | `vilma_radial_fe`: `toroidal_operator`; shared `bordered_band` |
| `58005d0` | `vilma_viscoelastic`: toroidal strain rows, norms, RHS, memory advance |
| `745f672` | `vilma_response`: W carried as drift once a 3-D element exists |
| `c9a9269` | the switch: the 3-D advance carries λ3, λ4; `&vilma l_toroidal` |
| `a6ccc3d` | restarts carry the channel count; 4 → 6 migration |
| `cbf77bd` | `&ctl file_hor`: horizontal displacement, total and toroidal |

Every refactor-only commit was checked variable by variable against a
`dump_reference` baseline of `main`; only the switch moves results, and only
for non-axisymmetric lateral viscosity.

**Where it departs from §2–§5, and why.**

- **The equations were re-derived**, not transcribed: with S⁽⁰⁾ = e_r×∇₁Y and
  the code's ½-symmetrised dyads, u = W e_r×∇₁Y has strain (W′−W/r) Z³ +
  (W/r) Z⁴, Z³ = sym(e_r ⊗ e_r×∇₁Y), Z⁴ = sym ∇₁(e_r×∇₁Y) = G e_θφ − H(e_θθ −
  e_φφ). That reproduces §2's eq-87 rows, B13 norms and dyadic map exactly.
  SHTns' toroidal field is −e_r×∇ (pinned exactly in `test_sht`).
- **W has its own operator** instead of a fifth interleaved field (§3.1): it
  couples to nothing on the left-hand side, so a tridiagonal per-degree system
  costs nothing in the spheroidal band, and the spheroidal operator is untouched.
  §6.3's multi-border question then lives in the small W system only.
- **The μ = 0 gauge** (§6.2): dead W dofs pinned, W = 0, as proposed.
- **The degree-1 gauge** (§2, eq 83): no net rotation per solid shell,
  w = ∫ψ r³ dr — derived from ∫ x × u dV = 0 for u = W e_r×∇₁Y₁ₘ. The ∫ψ r²
  transcribed in §2 cannot be that integral. One border per solid shell, since a
  stack with a solid inner core has two rigid rotations. A gauge only: a rigid
  rotation has no strain and moves no output.
- **4 or 6 channels at run time** instead of `NLAM 4 → 6` everywhere (§3.2):
  radially symmetric and laterally uniform fields force no W, so those runs
  keep four channels and pay nothing. `NLAM` still means the spheroidal four.
- **The modal path was removed** rather than fenced off (§3.5, C8).
- **§4's zero test was stated too broadly.** Toroidal flow vanishes for
  AXISYMMETRIC configurations, not for any shared mirror plane: under a
  reflection the toroidal potential is odd, so a Y₂₀ load over a cos 2φ field
  drives W ∝ sin 2φ. V4 needs no chirality.

**Validation** (`test_sht`, `test_tensor_sh`, `test_assembly`, `test_toroidal`,
`test_rotinv`, `test_restart`): V1 W stiffness = strain-route energy to 2e-16;
V2 six channels reproduce an arbitrary symmetric tensor to 9e-15 (four miss
35 %), channel cross-talk at nlat = 2·lmax+2 ≤ 8e-13 (§6.4's aliasing is not
there); V3 W = 0 exactly for axisymmetric configurations, off-pole W falls
2.6e-4 → 1.5e-6 from lmax 16 to 32; V4 W driven, reflection selection rule to
4e-16, W ∝ δ^0.987 and its uplift feedback ∝ δ^1.986. V6/V7 (VILMA1, the
deglaciation) are the paper's, not this note's.

**Context.** Volker Klemann, VILMA1's author, is part of the effort, and
VILMA2 is developed as an open-source alternative to VILMA1. That reframes
this work. The toroidal block is not a feature to match a competitor on — it
is physics VILMA1 has and VILMA2 must not lose. It also means the oracle problem of §4 is far less severe than
it looks (§2.1).

---

## 1. What is missing, and why it may matter

VILMA2 has no toroidal degree of freedom anywhere:

| | where |
|---|---|
| `NLAM = 4`, for Martinec's λ ∈ {1,2,5,6}; λ=3,4 dropped | `src/vilma_viscoelastic.f90:41` |
| DOF map is `4nr−1` — `idx_u/idx_v/idx_f/idx_p`, no `idx_w` | `src/vilma_radial_fe.f90:173-181` |
| strain coefficients, "Martinec eq 87, W dropped" | `src/vilma_viscoelastic.f90:318-330` |
| radial stiffness, "toroidal W block dropped" | `src/vilma_radial_fe.f90:209` |
| the scalar SHT computes the toroidal part and discards it | `src/vilma_sht.f90:237` |

The justification given in the code is *"spheroidal-only 1-D loading"*. That is
correct for a radially symmetric Earth and wrong once viscosity varies
laterally. Martinec (2000) says so himself, immediately after eq 110 (p. 130) —
the code cites the first half of the sentence and the second half is the point:

> the λ=1,2,5,6 terms "are completely decoupled from those with superscripts
> λ=3,4 … **A laterally heterogeneous viscosity removes both the decoupling and
> the degeneracy of the spheroidal and toroidal displacements.**"

Mechanically: a poloidal surface load over laterally varying viscosity drives
toroidal flow, and that flow back-couples into the poloidal field through the
same pointwise `M(θ,φ)·τ` product. VILMA2 computes the forward leg
implicitly in the grid product and then annihilates it in the analysis, so the
return leg never happens.

### 1.1 The motivating evidence has changed — re-establish it first

This work was scoped when VILMA2's lateral response looked severely damped
against VILMA1. **That turned out to be a bug in the viscosity reader, not
missing physics** (fixed in `dd03b79`; `input/bagge2021.nc` stores latitude
north-first and the whole field was collapsing onto one parallel). Corrected,
at 21 ka against the block-D deglaciation:

| lateral signal (3D − 1Db rsl, area-weighted) | rms | max abs | pattern corr. vs VILMA1 |
|---|---:|---:|---:|
| VILMA2, before the reader fix | 0.470 m | 15.3 m | −0.108 |
| VILMA2, after | 2.661 m | 95.0 m | **+0.963** |
| VILMA1, clean lateral signal | 2.915 m | 149.1 m | 1 |

So the rms is now within 9 % of VILMA1 and the pattern correlates at 0.96. **The
large deficit this document was written to explain no longer exists.**

What remains is a peak-amplitude gap: 95 m against 149 m, at r = 0.96. That is a
plausible home for the toroidal block, but it is equally a plausible home for the
elastic-lid exclusion (§1.2), and the two have not been separated. Caveats on the
number: one epoch, one metric, a 5 ka run rather than the full deglaciation.

**Do not start §5 until that residual has been attributed.** The cheap
discriminator is the lithosphere, not the toroidal block: see §1.2.

### 1.2 The cheaper suspect, to rule out first

`response_enable_lateral_visc_from_nodes` skips every element whose rheology is
elastic:

```fortran
if (self%MkPerDt(e) == 0.0_wp) cycle    ! elastic/fluid: stay as-is
```
`src/vilma_response.f90:1462`

Block D's layer 1 is `rheology = 0` to 80 km, so VILMA2 runs a laterally
uniform rigid plate exactly where the Bagge field has its strongest lateral
contrast, while VILMA1's lithosphere is defined *by* the viscosity file and sees
the weak zones. Two things follow:

- there is no laterally varying lithosphere **thickness** in the model at all —
  no parameter, no field, no reader (`earth_layer` carries scalar `r_bot`/`r_top`;
  `DEPTH_LITHO` at `src/vilma_radial_fe.f90:46` is a mesh-spacing constant);
- `isostasy_data/earth_structure/lithothickness/pan2022.nc` exists, i.e. the
  input is available and unusable.

Lithospheric thickness controls the flexural wavelength, which is what a
peak-amplitude gap at high pattern correlation looks like. Rule this out before
committing to §5.

---

## 2. The equations

Sourced from `doc/refs/Martinec2000.pdf` in a scoping pass. **Verify before
use.** Symbol extraction from the PDF is lossy; these were transcribed
conservatively.

### What the existing spheroidal code implements

| Martinec eq | content | where |
|---|---|---|
| B8/B10 | the six Regge–Wheeler tensor harmonics Z^(1..6) | `vilma_tensor_sh.f90:11-15` (4 of 6) |
| B11 | E, F, G, H angular functions | `vilma_tensor_sh.f90:160-181` |
| B13 | Z^λ:Z^λ orthogonality norms | `vilma_viscoelastic.f90:335` |
| 80 | δE_shear, U/V/W stiffness | `vilma_radial_fe.f90:268-290` (W omitted) |
| 81–84 | grav / press / uniq / surface forcing | `vilma_radial_fe.f90:292-350`, `:488-491` |
| 87–88 | strain coefficients (a,b,c) from nodal U,V,W | `vilma_viscoelastic.f90:318-328` (W rows omitted) |
| 102/107 | Maxwell memory advance (3-D pointwise / 1-D spectral) | `vilma_response.f90:1730-1752`, `vilma_viscoelastic.f90:383-488` |
| 109/110 | memory shape functions + spectral double-dot | `vilma_viscoelastic.f90:342-381` |

### The missing toroidal terms, same numbering

**B13 norms**, complete set λ=1..6, with `J = j(j+1)`:

```
[ 1, J/2, J/2, ½·J(J−2), 2J², 2J(J−2) ]
```

The code's `norm = [1, J/2, 2J², 2J(J−2)]` is entries 1,2,5,6. Martinec writes
the λ=4 and λ=6 norms as `½(j−1)j(j+1)(j+2)` and `2(j−1)j(j+1)(j+2)`; since
`(j−1)(j+2) = J−2` these agree with the code and with eq 110. **New entries:
`norm₃ = J/2`, `norm₄ = J(J−2)/2`.**

**B10 dyadic components** of the two missing channels:

- `Z³ = −F e_rθ + E e_rφ` — toroidal companion of `Z² = E e_rθ + F e_rφ`
- `Z⁴ = G e_θφ − H (e_θθ − e_φφ)` — toroidal companion of `Z⁶ = G(e_θθ − e_φφ) + 4H e_θφ`

Extending `tensor_sh_synth` (`vilma_tensor_sh.f90:135-158`), the full dyadic map
becomes:

```
rr  =  T¹Y
rθ  =  T²E − T³F          rφ =  T²F + T³E
θθ  = −J·T⁵Y + T⁶G − T⁴H
φφ  = −J·T⁵Y − T⁶G + T⁴H
θφ  =  4T⁶H + T⁴G
```

Cross-checked against eq 92's `b` dyad, which carries W explicitly
(`b_θθ ∋ −W^k H`, `b_θφ ∋ W^k G`, `b_φφ ∋ +W^k H`).

**Eq 87 toroidal strain rows** — the two rows dropped at
`vilma_viscoelastic.f90:326-327`:

```
λ=3:  a³ = −W^k + W^{k+1},   b³ = −W^k,   c³ = −W^{k+1}
λ=4:  a⁴ = 0,                b⁴ =  W^k,   c⁴ =  W^{k+1}
```

These depend **only** on W; the λ∈{1,2,5,6} rows depend only on U,V. The strain
map is exactly block-diagonal — this is what makes the extension tractable.

**Eq 80 toroidal stiffness**, the whole W contribution (factor μ_k):

```
J·( I¹_ab − I³_ab − I³_ba + I⁶_ab )·W^a δW^b  +  J(J−2)·I⁶_ab·W^a δW^b
```

No U, V, F or Π appears. All integrals already exist in
`src/vilma_radial_integrals.f90`.

**Eq 83 (δE_uniq)** carries a second rank-1 penalty at j=1 —
`(4π/3)·[Σ_a K³_a W^a_{1m}]·[Σ_b K³_b δW^b_{1m}]` — the rigid-**rotation** null
space, alongside the existing rigid-translation one
(`vilma_radial_fe.f90:354-372`). Not optional: at j=1, `norm₄ = 0`, so only λ=3
survives, `ε³ = dW/dr − W/r`, which vanishes identically for `W ∝ r`. The j=1 W
block **is** exactly singular.

**Eq 84 (δF_surf) has no δW term.** A surface load never forces W directly; W is
driven only through the dissipative RHS (eq 110's λ=3,4 terms), i.e. only once
lateral viscosity has mixed the memory. Two consequences:

- the elastic gains `gu/gn/gv` (`vilma_response.f90:1112-1119`) are untouched, so
  **every elastic benchmark must stay bit-identical**;
- W has no elastic-gain counterpart to `xUn`/`xVn` — it comes entirely from the
  drift solve.

### Not sourced — treat as open

- **The SHTns↔Martinec toroidal convention.** From B1–B3,
  `Λ_Ω Y = e_r × ∇_Ω Y` has components `(−F, E)`; SHTns' `SHtor_to_spat`
  produces `((1/sinθ)∂_φT, −∂_θT) = (F, −E)`, i.e. `−Λ_Ω`. That is a derivation,
  not a citation. Pin it with the same per-degree calibration the code already
  uses for the spin-2 channel (`vilma_tensor_sh.f90:85-110`) — do not assume it.
- **Order-in-contrast.** "Toroidal flow appears at first order in the viscosity
  contrast, back-coupling at second order" is standard perturbation reasoning,
  not stated in Martinec (2000). Used in §4 only to predict a scaling exponent.
*(A third item, whether VILMA1 carries the toroidal block, is now settled — see
§2.1.)*

### 2.1 VILMA1 does carry the toroidal treatment — confirmed

Volker Klemann, VILMA1's author, has confirmed that VILMA1 treats the toroidal
component. This matters three ways:

- it removes the last doubt about the physics: the reference the model is
  measured against solves a problem VILMA2 does not;
- it makes VILMA1 a **valid oracle** for §4/V6, which is the only external check
  available for this work;
- it means the toroidal implementation can be developed against a working one
  rather than from the paper alone. Volker is joining the effort, so the sign
  and normalisation conventions of §2 — the part this document is least sure of
  — are answerable by asking rather than by re-deriving.

Nothing in `doc/vilma-v1-backend.md` or `src/vilma_v1.f90` records this; it is
written down here because the code cannot tell you.

---

## 3. Blast radius

**Bottom line: an extension, not a reformulation.** ~600–750 lines of production
change across 8 files, roughly 80 % mechanical, plus ~400–500 lines of tests.

### 3.1 DOF layout — `src/vilma_radial_fe.f90`

Adding `idx_w` makes the stride 5 and `ndof_of = 5nr − 1`. Mechanical at the
index level: every consumer goes through the `idx_*` accessors, never a literal
stride (verified across `vilma_viscoelastic`, `vilma_response`, `vilma_rotation`,
`vilma_modal`, `tests/`). `band_build` infers `kl`/`ku` from the COO pattern
(`vilma_band.f90:38-47`), so bandwidth grows automatically — roughly 2.2× per
radial solve. Assembly widens `Aloc(7,7)→(9,9)`, `gmap(7)→(9)`
(`vilma_radial_fe.f90:230-231`).

**Two real design questions, both for the author:**

1. **μ = 0 regions have no W stiffness at all.** With `build_M3L70V01` the core
   is `r < 3480 km` with μ = 0 (`vilma_earth_structure.f90:232`) — ~40 % of the
   218-node mesh. Every interior core W node gets an all-zero row and the
   operator is singular. `tests/test_assembly.f90` check (2b) will fail
   immediately, which is the right place for it to surface. Options: (a) pin
   dead W dofs with an identity row — a genuine Dirichlet condition on a dof the
   physics does not possess, not a guard clause; (b) a variable DOF map omitting
   W where no adjacent element has μ > 0, which breaks the uniform stride and
   costs much more. (a) looks right, but it is a judgement call.
2. **A second KKT border row at j=1.** `radial_operator` supports exactly one
   border (`self%w`, `self%bordered`, `ndof_solve = nd+1`,
   `vilma_radial_fe.f90:402-432`). Eq 83's W penalty needs a second, orthogonal
   one. This collides with the degree-1 machinery (`nullmode`, the `border=`
   argument on `radial_operator_solve_vec`, `response_deg1_to_cm`). Generalising
   `w(:)` → `w(:,nborder)` is ~40 lines but touches the equilibration logic,
   which folds the border row/column into the row/column maxima. **This is the
   one piece that is architecture rather than transcription.**

### 3.2 `src/vilma_viscoelastic.f90` — NLAM 4→6

Most of it is already λ-generic (`do m = 1, NLAM` at `:364, :369, :438, :445,
:459, :471, :479`). What hard-codes 4: `NLAM` (`:41`); `strain_coeffs`
(`:318-328`, gains `w1,w2` — **signature change, 6 call sites**); the `norm`
literal (`:335`); `sa/sb/sc(4,NLAM)` (`:65, :334`, and `ve_strain_constants`
`:330-340`) where the 4 is the *test dofs* and becomes 6, mirrored at
`vilma_response.f90:1092`; `dissipative_rhs` (`:342-381` — `floc(4)`, `do t = 1,4`,
and two new `idx_w` scatter targets); `advance_memory` (`:383-488`, needs
`Wn`/`Wn_prev`).

### 3.3 `src/vilma_tensor_sh.f90` — lower cost than it looks

- **SHTns already has both toroidal transforms and the code already calls one.**
  `spat_to_SHsphtor` at `vilma_sht.f90:237` computes `tlm` into a local and
  discards it; `SHtor_to_spat` exists in `shtns.f03:262` unwrapped. Z³ costs two
  thin wrappers (~35 lines) and no new mathematics.
- **Z⁴ is nearly free.** `spin2_synth` (`:160-181`) already computes both
  `sg = Σc·G` and `sh = Σc·H` and uses only `(sg, −sg, 4sh)`. Z⁴ is the other
  combination of the same two fields: `(−sh, +sh, sg)`. `spin2_adjoint`
  (`:222-257`) likewise already has both `S_g*` and `S_h*`. The per-degree
  normalisation `n4(l)` calibrates exactly as `n6` does.
- `TLAM = 4` → 6 (`:40`); synth/analysis grow ~15 lines each.
- **Verify, do not assume:** B12 guarantees `∫Z⁴:Z⁶ = 0`, so the adjoint
  projection stays diagonal *provided the grid quadrature is exact*. The spin-2
  channel already wants `nlat = 3·lmax` to de-alias
  (`tests/test_benchmark_lvz.f90:35`, and note production runs `2·lmax`). A
  second spin-2 channel on the same grid is the natural place for silent λ4↔λ6
  aliasing — test the off-diagonal explicitly.

### 3.4 `src/vilma_response.f90` — largest diff, shallow

- **`advance_shape_tensor` (`:1730-1752`) needs zero change** — it operates on
  the 6 dyadic planes and is λ-agnostic. Same for the TRAP twin (`:1858-1881`).
- `gather_tensor_coeffs` (`:1697-1728`) and its TRAP twin (`:1801-1856`) need
  nodal W. Because eq 84 gives no elastic W forcing, the affine form
  `σ·xUn + dUn_re` collapses for W to drift only: `W(e,k) = dWn_re(e,k)`. No
  `xWn` array.
- `solve_drift` (`:1176-1257`) extracts `idx_w(node)` into new `dWn_re/dWn_im`.
- All `(NLAM, ne, nk)` allocations are symbolic and grow automatically
  (`:1151-1153, :2041-2046, :2163-2165`). **Memory: +50 % on the dominant state
  array** — at lmax 170 roughly 2.4 GB → 3.6 GB for `Are…Cim` plus snapshots.
  Check against `doc/performance.md` before committing.
- `Mk3`/`MkPerDt3`, `response_enable_lateral_visc` and the stability estimate
  need **no** change: M is a scalar field applied pointwise to all six planes.

### 3.5 The modal path

*Superseded: the modal solver was removed from the code instead (§0).*

**Leave `vilma_modal` spheroidal-only with an explicit `error stop`.** Its lateral
model is a per-mode scalar rate modulation (`vilma_response.f90:1476-1566`,
`design-modal.md` §4) and is structurally incapable of poloidal→toroidal
coupling; extending it is a research project. Pin its packing
(`npk = 3·NLAM·nmax`, `vilma_modal.f90:57, :134, :163-189`) to the spheroidal λ
subset so its cost and results are literally unchanged.

### 3.6 Restart format

`src/vilma_io.f90:471` hard-errors if the file's `nlam` ≠ `NLAM`. Bumping to 6
**breaks every existing restart**; needs a migration read (old `nlam=4` → the
spheroidal channels, zero the rest). `get3d_pad` (`:601-617`) is the template.

### 3.7 Interaction with the degree-1 guard

`c36cc9b` added the `Jr <= 2` guard on `b(4)`/`c(4)` in `strain_coeffs`, because
at j=1 `norm₆ = 2J(J−2) = 0` and `Z⁶_{1m} ≡ 0` pointwise, so the 3-D dyadic path
drops a nonzero λ6 strain while the 1-D spectral advance retains it. **The
identical statement holds for the new λ=4 channel** (`norm₄ = ½J(J−2) = 0`,
`Z⁴_{1m} ≡ 0`), so the guard must be extended to cover it when λ=4 lands.

---

## 4. Validation — the biggest risk

**There is no toroidal oracle in this repo or its references.** This, not the
implementation, is what makes the work hard.

The existing suite provably cannot catch a toroidal error: toroidal flow
vanishes whenever load *and* viscosity share a mirror plane.
`tests/test_benchmark_lvz.f90:35` runs `mmax=0` (axisymmetric disc over
axisymmetric cylinder); `tests/test_rotinv.f90` is a rotation of an m=0 problem,
and the spheroidal/toroidal split commutes with rotation, so the rotated field
is still exactly toroidal-free; `tests/test_response_3d.f90` asserts reduction to
the 1-D answer under a *uniform* field. Martinec never exercised it either — of
model E he writes that its spheroidal displacements are decoupled from the
toroidal ones (p. 133).

| candidate reference | verdict |
|---|---|
| `data/benchmarks/love_M3-L70-V01`, `disc_spada2011`, `sle_martinec2018`, `rotation_spada2011` | all 1-D. Bit-identity regressions only. |
| Weerdesteijn (2023) §5.2 | axisymmetric cylinder + axisymmetric disc → toroidal ≡ 0 |
| Martinec (2018) VEGA | 1-D |
| Martinec (2000) models C, E | explicitly decoupled |
| **VILMA1 via `solver="v1"`** | **the only external oracle** — same driver, namelist, forcing, remap, output; an v2-vs-v1 comparison is a one-line namelist change. Needs ifx + a hand-installed `vega_pism.a`. |

### Proposed tests

- **V1 — W stiffness by independent re-derivation** (algebraic, no data). In
  `tests/test_assembly.f90`, rebuild the W block from the strain representation,
  `2∫μ Σ_{λ=3,4} norm_λ B^λ_i B^λ_j r² dr` (eqs 87/110), and require
  machine-precision agreement with the eq-80 assembly. This is the cross-check
  that localised the U–F self-gravity bug (`formulation.md:205-208`), and it is
  the highest-value test here: it pins the new stiffness against the new strain
  coefficients through two independent routes. Re-assert (2b) no empty rows and
  (2f) symmetry.
- **V2 — dyadic completeness** (decisive, no data). In
  `tests/test_tensor_sh.f90`, synthesise an arbitrary band-limited symmetric
  second-order tensor field (six independent dyadic planes), analyse,
  re-synthesise, require identity. **Impossible with 4 channels** — the current
  round trip only tests closure within a subspace. Add λ4↔λ6 and λ2↔λ3
  off-diagonal projections at round-off, and extend the `∫τ:τ` Parseval check to
  the full 6-entry B13 norm vector.
- **V3 — the zero test.** For any mirror-symmetric configuration the toroidal
  field must be identically zero. Add `max|W|` and `max|τ^{λ=3,4}|` ≤ round-off
  to `test_benchmark_lvz` and `test_rotinv`, and require uplift bit-identical to
  the pre-change values (−0.7306 / −1.2180 m, `design.md` §12). A sign or
  normalisation error almost certainly breaks this.
- **V4 — `test_toroidal_chiral`** (new; the core test). A configuration with **no
  mirror plane**: an axisymmetric cap load at the pole over an LVZ centred off
  the load axis with an additional azimuthal `sin(2φ)` component, so no
  reflection leaves both fields invariant. Assert (a) `max|W|` far above
  round-off; (b) **selection rule** — for a `Y_{20}` load and an `M` field with a
  single `Y_{2±2}` component, toroidal response appears only in the (l,m)
  permitted by the triple product, everything else at round-off; (c)
  **δ-scaling** — halve the log-viscosity contrast, `|W|` should halve (slope 1)
  while the toroidal-induced change in surface uplift falls ~4× (slope 2). The
  exponents come from perturbation reasoning, not Martinec, but the *ratio* of
  the two slopes is a real falsifiable prediction.
- **V5 — rotational invariance of a toroidal-carrying case.** V4's configuration
  on-pole vs off-pole, matching cap-centre uplift. **State in the test header**
  that this cannot detect a *missing* toroidal block, only a mis-normalised,
  mis-signed or mis-ordered one — the split is itself rotationally covariant.
- **V6 — VILMA1 cross-check.** V4's chiral configuration through
  `solver="v1"` and through VILMA2, same forcing and output grid. The
  acceptance gate for the whole piece of work, not an afterthought.
- **V7 — the payoff measurement.** Lateral-response amplitude with and without
  the toroidal block, against VILMA1, on the corrected baseline of §1.1. That
  number is the deliverable.

`tests/dump_reference.f90` already dumps the memory tensor and dyadic round trip
to NetCDF — use it to byte-compare the refactor-only commits.

---

## 5. Commit sequence

C1–C5 must be **bit-identical** on every existing test. C6 is the only commit
that changes any result, and only for non-mirror-symmetric viscosity.

| # | commit | changes results? | verified by |
|---|---|---|---|
| **C1** | `vilma_sht`: expose the toroidal vector transforms — return the `tlm` discarded at `vilma_sht.f90:237`, wrap `SHtor_to_spat` | no (pure addition) | `test_sht`: toroidal round trip, sph/tor orthogonality |
| **C2** | `vilma_tensor_sh`: `TLAM 4→6`, add Z³ and Z⁴ by recombining existing `spin2_synth`/`spin2_adjoint` outputs; calibrate `n3`, `n4` | no (callers still pass zeros) | **V2** |
| **C3** | `vilma_radial_fe`: `idx_w`, stride 4→5, eq-80 W block, and the μ=0 / degree-1 gauge decision (§3.1) | no — eq 84 never forces W | **V1**; `test_love`, `test_relax`, `test_benchmark_love`, `test_benchmark_disc` bit-identical |
| **C4** | `vilma_viscoelastic`: `NLAM 4→6`, eq-87 toroidal rows, B13 norms, 6 test dofs, `dissipative_rhs` W scatter; update `vilma_rotation`/`vilma_modal` call sites; extend the `Jr <= 2` guard to λ=4 (§3.7) | no — all W inputs zero | full `make check` bit-identical vs a `dump_reference` baseline |
| **C5** | `vilma_response`: carry W (`dWn_*`, `edWn_*` from `solve_drift`; `gather_tensor_coeffs` ×2) | no while η is radial | **V3** zero-assertions |
| **C6** | **the switch**: activate λ=3,4 in `advance_memory_3d` / `_trap` | **yes**, for non-mirror-symmetric lateral viscosity only | **V4**, **V5**; LVZ and rotinv still reproduce −0.7306 / −1.2180 m |
| **C7** | `vilma_io`: restart `nlam` 4→6 with a migration read instead of the hard error at `:471` | no | `test_restart` + a vendored legacy (nlam=4) restart |
| **C8** | `vilma_modal`: explicit `error stop` when toroidal is active; pin modal packing to the spheroidal subset | no | `test_modal`, `test_modal_resp`, `test_modal_visc3d` bit-identical |
| **C9** | docs: `formulation.md:35, 127-131`, `docs/discretization/spectral.qmd:29-35`, `docs/physics/rheology.qmd:70-80`, `design.md` §12, and this file | — | prose |
| **C10** | VILMA1 cross-check (**V6**) + the measurement (**V7**) — a results commit | — | side-by-side run |

If C3's gauge question needs the multi-border KKT generalisation, split it:
**C3a** = index layout + eq-80 W block + dead-dof pinning, j ≥ 2 only with an
`error stop` at j = 1; **C3b** = the second border row and the degree-1 rotation
gauge.

---

## 6. Open questions for the author

1. **Is this work justified?** §1.1 — the deficit that motivated it is gone.
   Attribute the residual peak-amplitude gap (95 m vs 149 m at r = 0.96) before
   starting. §1.2 is the cheaper suspect and should be eliminated first.
2. **The μ = 0 W gauge** — identity row for dead dofs, or a variable DOF map?
   (§3.1.1)
3. **The second KKT border row** — generalise `w(:)` → `w(:,nborder)`, or split
   C3 and defer degree 1? (§3.1.2)
4. **Grid resolution** — production runs `nlat = 2·lmax` while `design.md` §12
   says the spin-2 channel wants `3·lmax`. Adding a second spin-2 channel raises
   the stakes; settle the de-aliasing question independently of this work.
5. ~~Confirm VILMA1 carries the toroidal block before relying on it as the
   oracle.~~ **Settled** — confirmed by Volker Klemann, §2.1. The remaining
   question is a practical one: agree the sign and normalisation conventions of
   §2 with him before implementing C2, rather than re-deriving them from the
   paper.
