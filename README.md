# VILMA2

> **Under heavy development.** VILMA2 is not yet ready for use in a
> scientific production capacity. If you want to use VILMA1 for now, please
> contact Volker Klemann at GFZ.

**VILMA2** is a state-of-the-art but simple and fast **3D solid-Earth model** — a
visco-elastic deformation model coupled with the sea-level equation — intended
as an **open-source alternative to VILMA1** within the CLIMBER-X climate model.

The method is a new implementation of the **spectral–finite-element,
time-domain** approach of Martinec (2000): spherical harmonics horizontally,
finite elements radially, an incompressible Maxwell rheology integrated
explicitly in time, a self-consistent sea-level equation with migrating
coastlines, and rotational feedback. It is built **3D-ready from the start**
(laterally varying viscosity) and validated against the published GIA benchmarks.

**Documentation:** <https://fesmc.github.io/vilma2/> (physics,
discretization, benchmarks, install & run), rendered from the Quarto sources
under [`docs/`](docs/); see [doc/design.md](doc/design.md) for the design
rationale and method comparison, and
[doc/vilma-v1-backend.md](doc/vilma-v1-backend.md) for the optional VILMA1 backend
(`&vilma solver = "v1"`), which drives VILMA1 itself through this model's
driver, namelist, forcing, remap and output for a like-for-like comparison. It
is **off by default and is not a dependency**: it needs an explicit
`make vilma vilma_v1=1 VILMA_V1_ROOT=<install>` and a hand-installed VILMA1.

## Status

VILMA2 is under heavy development and not yet ready for scientific
production use (see the note at the top). Implemented and validated so far: the spectral–finite-element solver
core, viscoelastic time stepping, the self-consistent migrating-coastline
sea-level equation, rotational feedback (polar motion), and laterally varying
(3D) viscosity — plus restart, spin-up, online lon-lat→Gauss remapping, and a
host-coupling API. The toroidal flow driven by lateral viscosity is also
implemented (optional, off by default) but so far only unit-tested.

Validated against the Spada et al. (2011) and Martinec et al. (2018) community
benchmarks (radial Love numbers, disc-load response, sea-level equation) and
Spada test 3/2 (rotation). Cross-code validation of the 3D path is ongoing.

## Install

Dependencies are [fesm-utils](https://github.com/fesm-org/fesm-utils) (branch
`coords-dev`, providing FFTW, SHTns, the `fesmutils` helper and the `coords`
module) plus a system netCDF. **configme** clones/links the dependencies and
generates the machine/compiler Makefile:

```bash
configme install vilma2                      # resolve machine/compiler, clone deps
configme install vilma2 -m macbook -c gfortran
configme install vilma2 --link fesm-utils=/abs/path/to/fesm-utils   # reuse a checkout
```

## Build

```bash
make vilma        # -> bin/vilma.x        (standalone forced-run driver)
make vilma_mkref  # -> bin/vilma_mkref.x  (build a Gauss-grid reference)
make vilma_remap  # -> bin/vilma_remap.x  (offline lon-lat -> Gauss remap)
make check            # build + run the test suite
```

`make` switches: `debug=0|1|2`, `openmp=0|1` (default 1; build `openmp=1` for the
threaded degree loop at production resolutions).

## Configure & run

A run configuration has two namelist groups: `&vilma` (physics and numerics,
the `vilma_param_class` record a host fills in memory) and `&ctl` (standalone
run control: forcing, reference and output files, time window, `i_eq`,
restart). The complete, documented `&vilma` defaults are in
[`input/vilma_defaults.nml`](input/vilma_defaults.nml) and are loaded
automatically; a run config such as [`vilma.nml`](vilma.nml) (the default) or
[`examples/deglac_lgm.nml`](examples/deglac_lgm.nml) carries the `&ctl` group
plus only the `&vilma` parameters it overrides (yelmo `defaults_file`
convention). Time fields (`dt_*`, `time_*`) are given in **years** and
converted to SI seconds on load.

- **Earth structure** — `earth`: a named built-in (`"M3-L70-V01"`, `"PREM"`) or
  `"custom"` to assemble from the surface-first layer arrays.
- **Response solver** — `earth_response`: `"ve"` (full viscoelastic, default),
  `"elastic"`, `"null"`.
- **Memory scheme** — `scheme = "fe"` (1st-order explicit, default) or `"trap"`
  (2nd-order adaptive), advanced by the `vilma_timestep` controller.
- **Degree-1 frame** — `deg1_frame = "cm"` (default, geocenter motion kept) or
  `"cf"` (Spada disc-benchmark convention).
- **3D viscosity** — `l_visc_3d`, `visc_3d_file`, `visc3d_tol`; `l_toroidal`
  carries the toroidal flow (off by default).
- **Spin-up / restart** — `equil_time_max`, `pre_spinup_1d`, and
  `restart_in_file` (`&ctl`).
- **Rotation** — `rotation` (TPW feedback): on by default; `.false.` for the
  non-rotating benchmarks.
- **Backend** — `solver = "v2"` (default) or `"v1"` (see above).

Run the standalone driver on a run config:

```bash
./bin/vilma.x examples/deglac_lgm.nml    # default: vilma.nml
```

It reads a reference state and an ice-thickness forcing (`file_forcing`,
`h_ice(lon,lat,time)`; remapped from lon-lat on the fly by default), marches the
model across the forcing, and writes the diagnostic surface fields (`rsl`,
`z_bed`, …) to `file_out`, and optionally the horizontal displacement to
`file_hor`.

Or stage/submit runs and ensembles with **runme** (`-r` run, `-s` submit;
comma-lists in `-p` define ensemble dimensions):

```bash
runme -o runs/deglac -e main --omp 8 -r -p vilma.lmax=128 vilma.earth_response=ve
```

Embedding the model in a host (the CLIMBER-X coupling path) uses the same API
behind a single `use vilma`:

```fortran
use vilma
type(solid_earth) :: se
call vilma_par_load(se%par, "vilma.nml")
call solid_earth_init(se, z_bed_eq, h_ice_eq, grid=host_grid)  ! Gauss grid + remap
call solid_earth_update(se, h_ice, dt_yr)   ! advance dt_yr [years]; read se%rsl, se%z_bed
call solid_earth_finalize(se)
```
