module vilma_params
   !! Physics + numerics configuration record for the solid-Earth model, loaded from
   !! one namelist group `&vilma` (yelmo convention: a flat parameter type filled by
   !! nml_read, see fesm-utils/utils/src/nml.f90). This is the host API contract:
   !! the high-level system init (solid_earth_init) consumes the whole record and
   !! distributes the values to the sub-solvers, while the specific component inits
   !! keep their direct-argument signatures.
   !!
   !! A host model (CLIMBER-X) fills this record in memory and drives the model
   !! through the API (z_bed_eq/h_ice_eq passed to solid_earth_init, each interval to
   !! solid_earth_update). The standalone driver reads it from a file. Run-management
   !! settings the *executables* need — forcing/reference/output paths, the time
   !! window, online-remap and i_eq selectors, restart-in path — live in the separate
   !! `&ctl` group (vilma_control), which a host never reads.
   use vilma_precision, only: wp
   use vilma_constants, only: sec_per_year
   use nml
   implicit none
   private

   public :: vilma_param_class, vilma_par_load, vilma_par_print, expand_path
   public :: MAX_LAYER

   integer, parameter :: MAX_LAYER = 16   !! cap on custom earth-structure layers

   type :: vilma_param_class
      ! --- solid-earth solver backend --------------------------------------------
      character(len=16) :: solver = "v2"
         !! which solid-earth solver sits behind the vilma_coupling API:
         !!   "v2"    (default) — this model's native Gauss-grid FE/SLE solver.
         !!   "v1"    — the VILMA1 library (Martinec/Klemann), driven through the
         !!             SAME driver, namelist, forcing, remap and output. Available
         !!             ONLY in a build made with `make vilma vilma_v1=1
         !!             VILMA_V1_ROOT=<install>`; the default build compiles a stub that
         !!             aborts with an actionable message. See doc/vilma-v1-backend.md.
         !! Everything below in the &vilma record that describes the VILMA2
         !! solver (earth structure, scheme, response kind, SLE, adaptive Δt,
         !! rotation, spin-up) is IGNORED when solver="v1": VILMA1 has its own
         !! earth structure, its own sea-level equation and its own time stepping,
         !! configured through the vilma_v1_* settings below and its own input files.
         !! Only lmax/nlat/nphi (the Gauss grid the coupling and the output live on)
         !! and l_visc_3d are shared.

      ! --- spectral grid (vilma_sht) ------------------------------------------------
      integer  :: lmax  = 0       !! maximum spherical-harmonic degree (required)
      integer  :: nlat  = 0       !! Gauss latitudes  (0 => de-aliased default 2*lmax+2)
      integer  :: nphi  = 0       !! longitudes       (0 => de-aliased default 4*lmax)
      integer  :: mmax  = -1      !! maximum order    (<0 => = lmax)
      integer  :: mres  = 1       !! order stride
      real(wp) :: eps_polar = -1.0_wp   !! polar-optimization threshold (<0 => library default)

      ! --- earth structure (vilma_earth_structure) ---------------------------------
      character(len=64) :: earth = "M3-L70-V01"
         !! named built-in model, or "custom" to build from the layer arrays below
      integer  :: n_layer = 0                  !! # custom layers (surface-first)
      real(wp) :: r_earth = 6371.0e3_wp        !! custom: surface radius [m]
      real(wp) :: r_core  = 3480.0e3_wp        !! custom: core-mantle boundary radius [m]
      real(wp) :: r_bot(MAX_LAYER) = 0.0_wp    !! custom: layer inner radii [m]
      real(wp) :: r_top(MAX_LAYER) = 0.0_wp    !! custom: layer outer radii [m]
      real(wp) :: rho(MAX_LAYER)   = 0.0_wp    !! custom: layer densities [kg m^-3]
      real(wp) :: mu(MAX_LAYER)    = 0.0_wp    !! custom: layer shear moduli [Pa]
      real(wp) :: eta(MAX_LAYER)   = 0.0_wp    !! custom: layer viscosities [Pa s]
      integer  :: rheology(MAX_LAYER) = 1      !! custom: 0=elastic 1=Maxwell 2=fluid

      ! --- viscoelastic memory scheme (vilma_response / vilma_viscoelastic) ------------
      character(len=8) :: scheme = "fe"        !! fe | etd1 | trap | be
      integer  :: max_couple_iter = 20         !! SLE<->memory co-convergence cap (implicit schemes)

      ! --- response kind selector (vilma_response) ----------------------------------
      character(len=8)  :: earth_response = "ve"  !! ve | elastic | null

      ! --- sea-level equation (vilma_sle) ------------------------------------------
      integer  :: sle_n_outer      = 3
      integer  :: sle_n_inner      = 20
      real(wp) :: sle_tol          = 1.0e-7_wp
      integer  :: sle_max_mem_iter = 20
      logical  :: sle_fixed_ocean  = .false.
      logical  :: sle_subgrid      = .true.

      ! --- adaptive time stepping (vilma_timestep) ----------------------------------
      ! The Δt fields are SI [s] in the record; the nml supplies them in YEARS and
      ! vilma_par_load converts on read (so the in-memory record is uniformly SI). The
      ! coupling cadence is NOT a parameter: the host (or the standalone driver's
      ! forcing axis) passes each interval to solid_earth_update.
      real(wp) :: dt_init   = 0.0_wp         !! first trial Δt [s] (0 => try the whole interval)
      real(wp) :: dt_min    = 0.0_wp         !! Δt floor [s] (0 => none)
      real(wp) :: dt_max    = huge(1.0_wp)   !! Δt ceiling [s]
      real(wp) :: rtol      = 1.0e-4_wp      !! relative local-error tolerance (memory ∞-norm)
      real(wp) :: atol      = 1.0e-3_wp      !! absolute local-error floor
      real(wp) :: safety    = 0.9_wp         !! step-size safety factor
      real(wp) :: grow_max  = 5.0_wp         !! max Δt growth per accepted step
      real(wp) :: shrink_min = 0.2_wp        !! min Δt shrink per step
      real(wp) :: cfl       = 1.0_wp         !! explicit (fe) sub-step Maxwell-number ceiling M=μΔt/η

      ! --- rotational feedback (vilma_rotation) ------------------------------------
      logical  :: rotation = .true.          !! TPW feedback (on for real runs; off for non-rotating benchmarks)
      real(wp) :: rotation_k_s = 0.0_wp      !! secular Love number k_s in the Liouville equation;
                                             !! <= 0 uses the model's own fluid limit k^T_f
      real(wp) :: rotation_c_minus_a = 2.63e35_wp !! C − A [kg m²], normalises the load inertia Ψ_L

      ! --- LGM-memory spin-up (vilma_coupling solid_earth_spinup) -------------------
      ! A model capability the host opts into. Relax under the start-slice ice while
      ! HOLDING the reference (z_bed_eq, h_ice_eq) as the datum, so the transient
      ! enters with viscous memory. The standalone driver triggers it; a host
      ! (CLIMBER-X) leaves equil_time_max=0 / pre_spinup_1d=.false. to skip it and
      ! takes whatever state it passes in as the (zero-memory) equilibrium.
      real(wp) :: equil_time_max = 0.0_wp
         !! >0: spin up the LGM memory state by holding the start-slice ice (relaxing to
         !! isostatic equilibrium) in the FULL model BEFORE the transient, with the
         !! reference held as the datum. A cap [years]: the relaxation exits early once
         !! the bed stops moving, or at this time with a warning if it has not
         !! converged. =0 skips the full-model equilibration phase. Non-default.
      real(wp) :: equil_rate_tol = 1.0e-3_wp
         !! spin-up convergence criterion: the mean bed velocity over a relaxation pass
         !! [m/yr]. Relaxation exits once the rate drops below this.
      logical :: pre_spinup_1d = .false.
         !! .true.: before the full-model phase, run a cheap 1-D pre-equilibration to
         !! bed-stationary convergence (the 1-D radial viscosity is the lateral
         !! geometric mean of the 3-D field), then switch to the full model carrying the
         !! spun-up memory. Independent of equil_time_max.

      ! --- 3D viscosity field + uncertainty sampling (vilma_earth_structure) --------
      ! Mirrors the CLIMBER-X VILMA1 scheme (src/geo/vilma.F90) but with a RELATIVE
      ! 1-sigma instead of a constant floor: perturb log10(eta) by f_visc_sd*sigma,
      ! sigma read from the file if name_visc_sd is set, else f_visc_rel*log10(eta).
      ! Degree-1 reference frame. "cm" (default) puts both the displacement and
      ! the geoid in the centre-of-mass frame, in which the solid Earth translates --
      ! geocenter motion, a real part of the degree-1 sea-level fingerprint, and
      ! what VILMA1 computes (it reports the term in vega_deg1.dat); it cuts the
      ! disc residual against VILMA1 ~5x. "cf" keeps the historical behaviour:
      ! the degree-1 displacement gauge is the solver's own w'd = 0 (no
      ! volume-integrated translation, centre-of-figure-like) while the geoid is
      ! referenced to CM (N1 = 0), so rsl carries NO degree 1 at all.
      !
      ! This changes degree 1 ONLY; every degree >= 2 is bit-identical. Use "cf"
      ! to reproduce the block A disc benchmark, which was validated with N1
      ! dropped (Spada's tables and tests/test_benchmark_love start at degree 2).
      character(len=8) :: deg1_frame = "cm"   !! "cm" | "cf"
      logical  :: l_visc_3d   = .false.   !! load a lateral log10(eta) field
      character(len=512) :: visc_3d_file  = ""       !! lon-lat-r log10(eta) field
      character(len=64)  :: name_visc     = "eta"    !! viscosity var (log10 Pa s)
      character(len=64)  :: name_visc_lon = "lon"
      character(len=64)  :: name_visc_lat = "lat"
      character(len=64)  :: name_visc_r   = "r"
      character(len=64)  :: name_visc_sd  = ""       !! optional sigma(log10 eta) var; "" => relative
      real(wp) :: f_visc_sd      = 0.0_wp   !! perturbation in units of sigma (0 = mean field)
      real(wp) :: f_visc_rel     = 0.1_wp   !! relative sigma = f_visc_rel*log10(eta) when no sd var
      real(wp) :: visc_log10_min = 19.5_wp  !! floor on log10(eta) after read + perturbation [dex]
      real(wp) :: visc_log10_max = 30.0_wp  !! ceiling on log10(eta) after read + perturbation [dex]
      real(wp) :: visc3d_lid_depth    = 0.0_wp  !! [m] lid rule off at 0; see visc3d_lid_log10max
      real(wp) :: visc3d_lid_log10max = 22.0_wp !! Maxwell elements lying wholly above visc3d_lid_depth
         !! whose log10(eta) exceeds this are set to visc_log10_max (effectively elastic). Such
         !! elements relax over Myr, so a held-load spin-up never equilibrates them; the rule
         !! splits the lid into elements that relax within kyr and ones that do not.
      real(wp) :: visc3d_tol     = 1.0e-3_wp !! lateral log10(eta) spread [dex] above which a radial
         !! element is treated as genuinely 3-D (pays the dyadic SHT round-trip); below it the
         !! element collapses to its lateral-mean scalar rate (cheap degree-diagonal path). Raising
         !! it demotes weakly-3-D elements to 1-D and cuts the memory-advance cost (the dominant cost).
      logical  :: l_toroidal     = .false.  !! carry the toroidal degree of freedom once a 3-D element
         !! exists (Martinec 2000 after eq 110). Off by default: it is slower and its effect on
         !! the solution is not yet shown (design-toroidal.md V7). .false. is the spheroidal-only model.

      ! --- VILMA1 backend (solver="v1" only; vilma_v1) -------------------------
      ! INERT unless solver="v1". These mirror the settings CLIMBER-X's VILMA1
      ! wrapper (src/geo/vilma.F90) hard-codes or takes from geo_params, so the two
      ! wrappers can be configured to agree exactly. VILMA1's spectral resolution is a
      ! RUNTIME setting (vg%jmax written to VILMA1's stdin file), not compiled in.
      integer :: vilma_v1_jmax = 170
         !! VILMA1 spectral degree (vg%jmax). CLIMBER-X uses 170. Must be consistent
         !! with vilma_v1_grid_file, which defines the grid VILMA1's fields come back on.
      character(len=512) :: vilma_v1_input_dir = "input/vilma"
         !! directory holding VILMA1's own static inputs: densi.inp, tint.inp,
         !! SLI_data.inp and the viscosity files below.
      character(len=512) :: vilma_v1_out_dir = "vilma_v1"
         !! scratch/output directory VILMA1 writes into (io.tmp, vega.lis, rsl.nc,
         !! dflag.nc, the ice-history NetCDF, restart files, ...). Created if absent.
      character(len=512) :: vilma_v1_grid_file = "input/vilma/vilma_grid.nc"
         !! NetCDF file carrying VILMA1's own lon/lat axes (its Gauss-Legendre grid at
         !! vilma_v1_jmax). vilma_v1 builds the VILMA2-Gauss <-> VILMA1-grid map pair
         !! from these axes; see doc/vilma-v1-backend.md for which grid each field is on.
      character(len=128) :: vilma_v1_visc_1d_file = "visko.inp"
         !! 1-D radial viscosity file, relative to vilma_v1_input_dir (io_visko).
      integer :: vilma_v1_l_prem = 1
         !! VILMA1 vg%l_prem. 1 (default, and what CLIMBER-X uses): VILMA1 generates
         !! its elastic structure from a polynomial PREM, and densi.inp supplies
         !! only the layer boundaries and the radial element sizes. 0: the rho and
         !! mu columns of densi.inp are used as given. Set 0 together with a
         !! densi.inp written from VILMA2's own layer table
         !! (experiments/make_vilma_densi.jl) to make the two backends share a
         !! radial structure, which is otherwise NOT matched -- only the
         !! viscosity is.
      integer :: vilma_v1_nsub = 1
         !! Number of VILMA1 sub-steps per coupling interval. VILMA1 enforces its own
         !! Maxwell stability condition at setup and ABORTS if its time step exceeds
         !! the shortest Maxwell time in the structure -- it does not sub-step
         !! itself. With the Bagge (2021) 3-D field that limit is short: measured
         !! 3.95 yr unfloored, 11.0 yr clamped at 1e19.5 (the clamp VILMA2
         !! applies), 26.2 yr at 1e20 -- all below the 100 yr GLAC-1D coupling
         !! interval, so the 3-D case cannot run at nsub = 1 with ANY of the
         !! available floors. Set nsub so that dt_coupling/nsub is below the
         !! reported minimum Maxwell time; VILMA1 prints both numbers when it
         !! refuses, so the required value is read straight off a failed run.
         !! The ice load is held across the sub-steps of one interval, which is
         !! VILMA1's own convention (see doc/vilma-v1-backend.md).
      character(len=128) :: vilma_v1_visc_3d_file = "visc3d_Bagge2021.nc"
         !! 3-D viscosity NetCDF, relative to vilma_v1_input_dir (io_nc3in). Read only
         !! when l_visc_3d = .true. (which sets VILMA1's vg%l_mod=1).
   end type vilma_param_class

contains

   subroutine vilma_par_load(p, filename, defaults_file, group)
      !! Fill the whole parameter record from the `&vilma` group of `filename`,
      !! overlaid on a complete `defaults_file` (yelmo convention): every parameter
      !! must exist in the defaults file, but the user `filename` may set only the
      !! subset it wants to override. If `defaults_file` is omitted, `filename` IS
      !! its own defaults — so it must then be complete. Override `group` to read a
      !! differently-named namelist.
      type(vilma_param_class), intent(inout) :: p
      character(len=*),     intent(in)    :: filename
      character(len=*),     intent(in), optional :: defaults_file
      character(len=*),     intent(in), optional :: group
      character(len=64)  :: g
      character(len=512) :: df
      real(wp) :: dt_init_yr, dt_min_yr, dt_max_yr
      real(wp) :: equil_time_max_yr

      g  = "vilma";      if (present(group))         g  = group
      df = filename;    if (present(defaults_file)) df = defaults_file
      call nml_set_verbose(.false.)             ! vilma_par_print echoes a concise summary instead

      ! solver backend ("v2" | "v1"); validated in solid_earth_init
      call nml_read(filename, g, "solver",    p%solver,    defaults_file=df)

      ! grid
      call nml_read(filename, g, "lmax",      p%lmax,      defaults_file=df)
      call nml_read(filename, g, "nlat",      p%nlat,      defaults_file=df)
      call nml_read(filename, g, "nphi",      p%nphi,      defaults_file=df)
      call nml_read(filename, g, "mmax",      p%mmax,      defaults_file=df)
      call nml_read(filename, g, "mres",      p%mres,      defaults_file=df)
      call nml_read(filename, g, "eps_polar", p%eps_polar, defaults_file=df)

      ! earth structure
      call nml_read(filename, g, "earth",     p%earth,     defaults_file=df)
      call nml_read(filename, g, "n_layer",   p%n_layer,   defaults_file=df)
      call nml_read(filename, g, "r_earth",   p%r_earth,   defaults_file=df)
      call nml_read(filename, g, "r_core",    p%r_core,    defaults_file=df)
      call nml_read(filename, g, "r_bot",     p%r_bot,     defaults_file=df)
      call nml_read(filename, g, "r_top",     p%r_top,     defaults_file=df)
      call nml_read(filename, g, "rho",       p%rho,       defaults_file=df)
      call nml_read(filename, g, "mu",        p%mu,        defaults_file=df)
      call nml_read(filename, g, "eta",       p%eta,       defaults_file=df)
      call nml_read(filename, g, "rheology",  p%rheology,  defaults_file=df)

      ! memory scheme
      call nml_read(filename, g, "scheme",          p%scheme,          defaults_file=df)
      call nml_read(filename, g, "max_couple_iter", p%max_couple_iter, defaults_file=df)

      ! response kind selector
      call nml_read(filename, g, "earth_response",  p%earth_response,  defaults_file=df)

      ! sea-level equation
      call nml_read(filename, g, "sle_n_outer",      p%sle_n_outer,      defaults_file=df)
      call nml_read(filename, g, "sle_n_inner",      p%sle_n_inner,      defaults_file=df)
      call nml_read(filename, g, "sle_tol",          p%sle_tol,          defaults_file=df)
      call nml_read(filename, g, "sle_max_mem_iter", p%sle_max_mem_iter, defaults_file=df)
      call nml_read(filename, g, "sle_fixed_ocean",  p%sle_fixed_ocean,  defaults_file=df)
      call nml_read(filename, g, "sle_subgrid",      p%sle_subgrid,      defaults_file=df)

      ! adaptive time stepping. The Δt fields are given in YEARS in the nml and
      ! converted to SI seconds here (the record is uniformly SI internally).
      dt_init_yr = p%dt_init/sec_per_year
      dt_min_yr  = p%dt_min /sec_per_year
      dt_max_yr  = p%dt_max /sec_per_year
      call nml_read(filename, g, "dt_init",    dt_init_yr,   defaults_file=df)
      call nml_read(filename, g, "dt_min",     dt_min_yr,    defaults_file=df)
      call nml_read(filename, g, "dt_max",     dt_max_yr,    defaults_file=df)
      p%dt_init = dt_init_yr*sec_per_year
      p%dt_min  = dt_min_yr *sec_per_year
      p%dt_max  = dt_max_yr *sec_per_year
      call nml_read(filename, g, "rtol",       p%rtol,       defaults_file=df)
      call nml_read(filename, g, "atol",       p%atol,       defaults_file=df)
      call nml_read(filename, g, "safety",     p%safety,     defaults_file=df)
      call nml_read(filename, g, "grow_max",   p%grow_max,   defaults_file=df)
      call nml_read(filename, g, "shrink_min", p%shrink_min, defaults_file=df)
      call nml_read(filename, g, "cfl",        p%cfl,        defaults_file=df)

      ! rotation
      call nml_read(filename, g, "rotation",   p%rotation,   defaults_file=df)
      call nml_read(filename, g, "rotation_k_s",       p%rotation_k_s,       defaults_file=df)
      call nml_read(filename, g, "rotation_c_minus_a", p%rotation_c_minus_a, defaults_file=df)

      ! LGM-memory spin-up (equil_time_max given in YEARS, converted to SI below)
      equil_time_max_yr = p%equil_time_max/sec_per_year
      call nml_read(filename, g, "equil_time_max", equil_time_max_yr, defaults_file=df)
      p%equil_time_max = equil_time_max_yr*sec_per_year
      call nml_read(filename, g, "equil_rate_tol", p%equil_rate_tol, defaults_file=df)
      call nml_read(filename, g, "pre_spinup_1d",  p%pre_spinup_1d,  defaults_file=df)

      ! 3D viscosity + uncertainty
      call nml_read(filename, g, "deg1_frame",     p%deg1_frame,     defaults_file=df)
      if (trim(p%deg1_frame) /= "cf" .and. trim(p%deg1_frame) /= "cm") &
         error stop 'vilma_params: deg1_frame must be "cf" or "cm"'
      call nml_read(filename, g, "l_visc_3d",      p%l_visc_3d,      defaults_file=df)
      call nml_read(filename, g, "visc_3d_file",   p%visc_3d_file,   defaults_file=df)
      p%visc_3d_file = expand_path(p%visc_3d_file)
      call nml_read(filename, g, "visc3d_tol",     p%visc3d_tol,     defaults_file=df)
      call nml_read(filename, g, "visc3d_lid_depth",    p%visc3d_lid_depth,    defaults_file=df)
      call nml_read(filename, g, "visc3d_lid_log10max", p%visc3d_lid_log10max, defaults_file=df)
      call nml_read(filename, g, "l_toroidal",     p%l_toroidal,     defaults_file=df)
      call nml_read(filename, g, "name_visc",      p%name_visc,      defaults_file=df)
      call nml_read(filename, g, "name_visc_lon",  p%name_visc_lon,  defaults_file=df)
      call nml_read(filename, g, "name_visc_lat",  p%name_visc_lat,  defaults_file=df)
      call nml_read(filename, g, "name_visc_r",    p%name_visc_r,    defaults_file=df)
      call nml_read(filename, g, "name_visc_sd",   p%name_visc_sd,   defaults_file=df)
      call nml_read(filename, g, "f_visc_sd",      p%f_visc_sd,      defaults_file=df)
      call nml_read(filename, g, "f_visc_rel",     p%f_visc_rel,     defaults_file=df)
      call nml_read(filename, g, "visc_log10_min", p%visc_log10_min, defaults_file=df)
      call nml_read(filename, g, "visc_log10_max", p%visc_log10_max, defaults_file=df)

      ! VILMA1 backend (inert unless solver="v1")
      call nml_read(filename, g, "vilma_v1_jmax",         p%vilma_v1_jmax,         defaults_file=df)
      call nml_read(filename, g, "vilma_v1_input_dir",    p%vilma_v1_input_dir,    defaults_file=df)
      p%vilma_v1_input_dir = expand_path(p%vilma_v1_input_dir)
      call nml_read(filename, g, "vilma_v1_out_dir",      p%vilma_v1_out_dir,      defaults_file=df)
      p%vilma_v1_out_dir = expand_path(p%vilma_v1_out_dir)
      call nml_read(filename, g, "vilma_v1_grid_file",    p%vilma_v1_grid_file,    defaults_file=df)
      p%vilma_v1_grid_file = expand_path(p%vilma_v1_grid_file)
      call nml_read(filename, g, "vilma_v1_visc_1d_file", p%vilma_v1_visc_1d_file, defaults_file=df)
      call nml_read(filename, g, "vilma_v1_visc_3d_file", p%vilma_v1_visc_3d_file, defaults_file=df)
      call nml_read(filename, g, "vilma_v1_l_prem",       p%vilma_v1_l_prem,       defaults_file=df)
      if (p%vilma_v1_l_prem /= 0 .and. p%vilma_v1_l_prem /= 1) &
         error stop 'vilma_params: vilma_v1_l_prem must be 0 or 1'
      call nml_read(filename, g, "vilma_v1_nsub",         p%vilma_v1_nsub,         defaults_file=df)
      if (p%vilma_v1_nsub < 1) error stop 'vilma_params: vilma_v1_nsub must be >= 1'
   end subroutine vilma_par_load

   function expand_path(path) result(out)
      !! Expand a leading `~/` or `$HOME/` in a file path using the HOME
      !! environment variable. Fortran `open` does no shell expansion, so nml
      !! paths such as `~/data/visc.nc` would otherwise be opened literally.
      !! Any other path (relative or absolute) is returned unchanged.
      character(len=*), intent(in) :: path
      character(len=512)           :: out   ! matches the path fields' length
      character(len=512) :: home
      integer :: n, status

      out = path
      if (len_trim(path) == 0) return
      if (path(1:2) == "~/") then
         n = 1
      else if (len_trim(path) >= 6) then
         if (path(1:6) == "$HOME/") then
            n = 5
         else
            return
         end if
      else
         return
      end if
      call get_environment_variable("HOME", home, status=status)
      if (status /= 0 .or. len_trim(home) == 0) return
      out = trim(home)//path(n+1:len_trim(path))
   end function expand_path

   subroutine vilma_par_print(p, unit)
      !! Echo the active configuration (to stdout, or `unit` if given).
      type(vilma_param_class), intent(in) :: p
      integer, intent(in), optional :: unit
      integer :: u, k
      u = 6;  if (present(unit)) u = unit

      write(u,'(a)')          ' [vilma] configuration'
      write(u,'(a,a)')        '   solver: ', trim(p%solver)
      write(u,'(a,i0,a,i0,a,i0)') '   grid:   lmax=', p%lmax, '  nlat=', p%nlat, '  nphi=', p%nphi
      write(u,'(a,a)')        '   earth:  ', trim(p%earth)
      if (trim(p%earth) == "custom") then
         do k = 1, p%n_layer
            write(u,'(a,i0,a,es9.2,a,es9.2,a,f8.1,a,es9.2,a,es9.2,a,i0)') &
                 '     layer ', k, ': r=[', p%r_bot(k), ',', p%r_top(k), &
                 ']  rho=', p%rho(k), '  mu=', p%mu(k), '  eta=', p%eta(k), &
                 '  rheol=', p%rheology(k)
         end do
      end if
      write(u,'(a,a)')        '   response: ', trim(p%earth_response)
      write(u,'(a,a,a,i0)')   '   scheme: ', trim(p%scheme), '   max_couple_iter=', p%max_couple_iter
      write(u,'(a,i0,a,i0,a,es8.1,a,l1,a,l1)') &
           '   sle:    n_outer=', p%sle_n_outer, '  n_inner=', p%sle_n_inner, &
           '  tol=', p%sle_tol, '  fixed_ocean=', p%sle_fixed_ocean, '  subgrid=', p%sle_subgrid
      write(u,'(a,es9.2,a,es8.1,a,es8.1,a,f5.2)') &
           '   dt:     init=', p%dt_init, &
           '  rtol=', p%rtol, '  atol=', p%atol, '  cfl=', p%cfl
      write(u,'(a,l1,a,f7.4,a,es10.3)') '   rotation: ', p%rotation, &
           '  k_s=', p%rotation_k_s, '  C-A=', p%rotation_c_minus_a
      write(u,'(a,es9.2,a,es9.2,a,l1)') '   spinup: equil_time_max=', p%equil_time_max, &
           '  equil_rate_tol=', p%equil_rate_tol, '  pre_spinup_1d=', p%pre_spinup_1d
      if (p%l_visc_3d) then
         write(u,'(a,a)')   '   visc_3d: ', trim(p%visc_3d_file)
         write(u,'(a,f6.2,a,f6.2,a,f5.2,a,f5.2,a)') &
              '            f_visc_sd=', p%f_visc_sd, '  f_visc_rel=', p%f_visc_rel, &
              '  clamp=[', p%visc_log10_min, ',', p%visc_log10_max, ']'
         write(u,'(a,es9.2,a)') '            visc3d_tol=', p%visc3d_tol, ' dex (3-D split)'
         if (p%visc3d_lid_depth > 0.0_wp) &
            write(u,'(a,f6.1,a,f5.2,a)') '            lid rule: elements above ', &
                 p%visc3d_lid_depth*1.0e-3_wp, ' km with log10(eta) > ', p%visc3d_lid_log10max, &
                 ' -> clamp ceiling'
         write(u,'(a,l1)')      '            l_toroidal=', p%l_toroidal
      end if
      if (trim(p%solver) == "v1") then
         write(u,'(a,i0)')  '   vilma_v1: jmax=', p%vilma_v1_jmax
         write(u,'(a,a)')   '             input_dir = ', trim(p%vilma_v1_input_dir)
         write(u,'(a,a)')   '             out_dir   = ', trim(p%vilma_v1_out_dir)
         write(u,'(a,a)')   '             grid_file = ', trim(p%vilma_v1_grid_file)
         if (p%l_visc_3d) then
            write(u,'(a,a)') '             visc (3d) = ', trim(p%vilma_v1_visc_3d_file)
         else
            write(u,'(a,a)') '             visc (1d) = ', trim(p%vilma_v1_visc_1d_file)
         end if
      end if
   end subroutine vilma_par_print

end module vilma_params
