module vilma_coupling
   !! Top-level coupling API — the contract a host climate/ice model (CLIMBER-X)
   !! drives the solid-Earth model through. Mirrors the VILMA-v1 wrapper in CLIMBER-X
   !! (src/geo/vilma.F90): ice thickness goes in, relative sea level and bedrock
   !! elevation come out. The model OWNS its Gauss-Legendre transform grid (built
   !! from se%par at init) and OWNS the remap between the host grid and that Gauss
   !! grid — the host never sees the Gauss grid.
   !!
   !!   call vilma_par_load(se%par, cfg, defaults)        ! configuration -> se%par
   !!   call solid_earth_init(se, z_bed_eq, h_ice_eq, grid=host_grid)   ! reference state
   !!   ...
   !!   call solid_earth_update(se, h_ice, dt_yr)      ! advance dt_yr [years]; fills se%rsl, se%z_bed
   !!   ...
   !!   call solid_earth_finalize(se)                  ! frees the grid too
   !!
   !! Grids. If `grid` (a coords lon-lat grid) is passed at init and differs from the
   !! model's Gauss grid, the model builds a conservative (host->Gauss) + bilinear
   !! (Gauss->host) map pair (vilma_remap) and drives h_ice in / rsl out through it.
   !! If `grid` is absent (or already the Gauss grid) the fields are taken as-is
   !! (passthrough). Either way the host reads se%rsl / se%z_bed / se%bsl on the SAME
   !! grid it supplied; the Gauss-grid working state lives in se%gg (se%gg%rsl, ...).
   !! Only the smooth rsl perturbation is mapped back to the host; se%z_bed is
   !! reconstructed as se%z_bed_eq - se%rsl on the host's own (high-resolution) bed,
   !! so host bed detail is preserved.
   !!
   !! Reference (equilibrium) state. The supplied (z_bed_eq, h_ice_eq) IS the relaxed
   !! state: the Maxwell memory starts at zero, so with h_ice = h_ice_eq we get
   !! d_ice = 0, rsl = 0, z_bed = z_bed_eq. Departures from the reference ice load
   !! drive the deformation and sea-level change incrementally. The reference
   !! topography z_bed_eq doubles as the SLE's reference topo0.
   use vilma_precision,       only: wp
   use vilma_constants,       only: rho_ice, rho_water, sec_per_year, pi, kyr
   use vilma_params,          only: vilma_param_class
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, &
                                 sht_grid_synthesis, sht_grid_surface_integral
   use vilma_earth_structure, only: earth_model, build_earth, load_visc_3d
   use vilma_viscoelastic,    only: scheme_from_name
   use vilma_response,        only: response_destroy, response_enable_lateral_visc_from_nodes, response_radial_rate, response, &
                                 response_init_elastic, response_init_ve, response_init_null, RESP_VE
   use vilma_sle,             only: sle_solver, ocean_function
   use vilma_timestep,        only: stepper_advance, adaptive_stepper
   use vilma_rotation,        only: rotation_destroy, rotation_init, rotation_state, &
                                    rotation_set_rate
   use vilma_remap,           only: remap_ll_gauss, remap_init, remap_to_gauss, remap_to_ll
   use vilma_v1,           only: vilma_v1_backend, vilma_v1_init, vilma_v1_update, vilma_v1_finalize, &
                                 vilma_v1_require
   use coords,             only: grid_class
   implicit none
   private

   public :: solid_earth, gauss_state
   public :: solid_earth_init, solid_earth_update, solid_earth_finalize
   public :: solid_earth_enable_visc_3d, solid_earth_sync_host, solid_earth_spinup
   public :: solid_earth_check_solver
   public :: VILMA_UNSET

   !! Fill value for a diagnostic the ACTIVE backend does not produce. Only the
   !! VILMA-v1 backend leaves anything unset (see the "known gaps" note on the
   !! solid_earth type); with the native solver every field is populated. It is a
   !! wildly unphysical number on purpose, so an unset value cannot be mistaken for
   !! a computed zero, and it matches the missing-value convention used by the
   !! CLIMBER-X remapping layer. Writers/printers are expected to skip it.
   real(wp), parameter :: VILMA_UNSET = -9999.0_wp

   ! LGM-memory spin-up controls (the relaxation interval is internal — no user dt).
   real(wp), parameter :: SPINUP_DT0   = 200.0_wp   !! first relaxation interval [years]
   real(wp), parameter :: SPINUP_GROW  = 1.6_wp     !! interval growth per pass (the bed stiffens)
   integer,  parameter :: SPINUP_MAXP  = 60         !! internal cap on relaxation passes
   !! The convergence criterion (mean bed velocity over a pass [m/yr]) is the runtime
   !! parameter par%equil_rate_tol.

   type :: gauss_state
      !! Model fields on the model's own Gauss grid, where the physics runs. When the
      !! host supplies data already on the Gauss grid (passthrough), these mirror the
      !! host-grid fields one-to-one.
      real(wp), allocatable :: z_bed_eq(:,:)   !! relaxed bedrock [m] (nphi,nlat)
      real(wp), allocatable :: h_ice_eq(:,:)   !! reference grounded ice [m]
      real(wp), allocatable :: h_ice(:,:)      !! current grounded-ice thickness [m]
      real(wp), allocatable :: rsl(:,:)        !! relative sea level change [m] (full field)
      real(wp), allocatable :: z_bed(:,:)      !! bedrock = z_bed_eq − rsl [m]
      real(wp), allocatable :: C(:,:)          !! ocean function (1 ocean / 0 land)
      real(wp), allocatable :: C0(:,:)         !! REFERENCE ocean function, from (z_bed_eq,
                                               !! h_ice_eq). Both are fixed at init, so this
                                               !! is a RUN CONSTANT: built once in
                                               !! solid_earth_init, never refreshed.
   end type gauss_state

   type :: solid_earth
      type(vilma_param_class)     :: par       !! configuration record (one &vilma group); set by the host before init
      type(sht_grid), pointer  :: sht => null()  !! model-OWNED transform grid (built at init, freed at finalize)
      type(earth_model)        :: earth     !! radial (+ optional 3D) structure
      type(response)           :: resp      !! viscoelastic field driver (load → u, N)
      type(sle_solver)         :: sle       !! sea-level equation
      type(adaptive_stepper)   :: stepper   !! adaptive-Δt controller (vilma_timestep)
      type(rotation_state)     :: rotation  !! TPW feedback (ON by default; par%rotation)

      type(gauss_state)        :: gg        !! Gauss-grid working state (the physics)
      !! The SLE's converged spectral surface load [kg m⁻²] at the end of the last
      !! interval — the σ that rsl and z_bed were solved with, kept so a consumer
      !! (the horizontal-displacement output) can evaluate the response at it.
      complex(wp), allocatable :: sigma_lm(:)

      !! Optional VILMA-v1 backend (par%solver = "v1"), for a like-for-like
      !! comparison behind this same API. When active, `resp`/`sle`/`stepper`/
      !! `rotation`/`earth` above are NOT built — VILMA-v1 owns the physics — and the
      !! diagnostics they would otherwise fill are left at VILMA_UNSET. Specifically:
      !!   worst_mass_resid : no analogue (VILMA-v1 reports no SLE mass residual)   -> VILMA_UNSET
      !!   stepper%*        : VILMA-v1 does its own internal time stepping          -> VILMA_UNSET
      !!   resp%t_* / sle%t_* / stepper%t_guard : those phases do not exist      -> stay 0
      !! while `rsl`, `z_bed`, `C` and `bsl` ARE populated (see solid_earth_update).
      !! `t_solver` below is the one timer that is meaningful for both backends.
      logical                  :: use_vilma_v1 = .false.
      type(vilma_v1_backend)   :: vilma_v1

      ! host-grid coupling: remap between the host grid and the model Gauss grid
      logical                  :: remap = .false.  !! host grid /= Gauss grid?
      type(remap_ll_gauss)     :: map               !! conservative-in / bilinear-out map pair

      !! PROFILE: wall-clock [s] for the two pieces of solid_earth_update that are
      !! neither the solver nor the SLE, and so otherwise vanish into the driver's
      !! residual bucket. See sle_solver's accumulators for the SLE's own split.
      real(wp)                 :: t_remap = 0.0_wp  !! host<->Gauss remap, in and out
      real(wp)                 :: t_solver = 0.0_wp !! wall-clock [s] inside the backend's own
         !! solve, comparable across backends: for "v1" it is VILMA-v1's time_evolution
         !! (excluding the Gauss<->VILMA-v1 remap, which is in %t_remap); for the native
         !! solver it is left 0 because the finer-grained split above already covers it.

      ! host-grid I/O fields (== the gg fields when passthrough); the host reads these
      real(wp), allocatable    :: z_bed_eq(:,:)  !! relaxed bedrock [m], host grid (kept at full host resolution)
      real(wp), allocatable    :: h_ice_eq(:,:)  !! reference grounded ice [m], host grid
      real(wp), allocatable    :: h_ice(:,:)     !! current grounded ice [m], host grid
      real(wp), allocatable    :: rsl(:,:)       !! relative sea level change [m], host grid
      real(wp), allocatable    :: z_bed(:,:)     !! bedrock = z_bed_eq − rsl [m], host grid

      ! clock + diagnostics
      real(wp) :: time = 0.0_wp               !! model time [years]
      real(wp) :: worst_mass_resid = 0.0_wp   !! worst SLE mass residual over the last interval
      real(wp) :: bsl = 0.0_wp                !! barystatic sea level vs reference [m] (eustatic equivalent)
   end type solid_earth

contains

   subroutine solid_earth_init(self, z_bed_eq, h_ice_eq, grid, h_ice_init)
      !! Build the model from its parameter record (self%par, set by the host before
      !! this call) and set the reference state. The model builds and owns its Gauss
      !! transform grid (from self%par%lmax/nlat/nphi) and, when `grid` differs from
      !! it, the host<->Gauss remap. Reference fields are supplied on the host grid
      !! (or directly on the Gauss grid when `grid` is absent / matches). The model is
      !! built in its full configuration per self%par; the 1-D pre-spin-up
      !! (solid_earth_spinup) toggles the lateral viscosity internally, not here.
      type(solid_earth),   intent(inout)        :: self
      real(wp),             intent(in)           :: z_bed_eq(:,:)  !! relaxed bedrock [m] (host grid)
      real(wp),             intent(in)           :: h_ice_eq(:,:)  !! reference grounded ice [m] (host grid)
      type(grid_class),     intent(in), optional :: grid           !! host lon-lat grid; absent => data already on Gauss
      real(wp),             intent(in), optional :: h_ice_init(:,:)!! current ice at t0 (host grid); default h_ice_eq
      integer  :: np, nl

      call solid_earth_finalize(self)                       ! clean slate (safe on a fresh object)

      ! which solver sits behind this API (see solid_earth_check_solver; the
      ! standalone driver calls it earlier still, right after reading the config)
      call solid_earth_check_solver(self%par)
      self%use_vilma_v1 = (trim(self%par%solver) == "v1")

      ! the model owns its Gauss grid, sized from the parameter record
      allocate(self%sht)
      call build_sht(self%par, self%sht)
      np = self%sht%nphi;  nl = self%sht%nlat
      self%time = 0.0_wp

      ! host<->Gauss remap: build it when a host grid is supplied that is not the Gauss
      ! grid; otherwise the host fields are taken to be on the Gauss grid already.
      self%remap = .false.
      if (present(grid)) self%remap = .not. grid_is_gauss(grid, self%sht)
      if (self%remap) call remap_init(self%map, self%sht, real(grid%G%x, wp), real(grid%G%y, wp))

      ! reference (equilibrium) state on the host grid (kept at full host resolution),
      ! and its Gauss-grid image used by the physics.
      allocate(self%z_bed_eq, source=z_bed_eq)
      allocate(self%h_ice_eq, source=h_ice_eq)
      allocate(self%gg%z_bed_eq(np,nl), self%gg%h_ice_eq(np,nl))
      call to_gauss(self, z_bed_eq, self%gg%z_bed_eq, conserve_mass=.false.)   ! bed: geometry
      call to_gauss(self, h_ice_eq, self%gg%h_ice_eq, conserve_mass=.true.)    ! ice: mass

      if (self%use_vilma_v1) then
         ! VILMA-v1 owns the physics, so NONE of the native sub-solvers (earth
         ! structure, response, SLE, adaptive stepper, rotation) are built — that is
         ! the point of the backend swap. The response is initialised NULL so the I/O
         ! layer sees a memoryless model and writes only the common diagnostics:
         ! there is no VILMA Maxwell memory to persist, and VILMA-v1 persists its
         ! own state through its own restart files.
         call response_init_null(self%resp)
      else
         call build_solver(self, np, nl)
      end if

      ! current state — seeded at the reference (or at h_ice_init if given). At the
      ! reference nothing has moved: rsl = 0, z_bed = z_bed_eq, memory zero. The ocean
      ! function is seeded from the reference flotation so a dt=0 seed step yields a
      ! physical barystatic diagnostic (update_bsl needs C).
      allocate(self%gg%h_ice(np,nl), self%gg%rsl(np,nl), self%gg%z_bed(np,nl), &
               self%gg%C(np,nl), self%gg%C0(np,nl))
      allocate(self%sigma_lm(self%sht%nlm));  self%sigma_lm = (0.0_wp, 0.0_wp)
      if (present(h_ice_init)) then
         call to_gauss(self, h_ice_init, self%gg%h_ice, conserve_mass=.true.)
         allocate(self%h_ice, source=h_ice_init)
      else
         self%gg%h_ice = self%gg%h_ice_eq
         allocate(self%h_ice, source=h_ice_eq)
      end if
      self%gg%rsl   = 0.0_wp
      self%gg%z_bed = self%gg%z_bed_eq
      ! The reference ocean function is a run constant (z_bed_eq and h_ice_eq are
      ! fixed from here on), so build it once and keep it. The initial C is that
      ! same field, which is exactly what the old unconditional call produced.
      call ocean_function(self%gg%z_bed_eq, self%gg%h_ice_eq, self%gg%C0)
      self%gg%C = self%gg%C0

      ! host-grid current state at the reference
      allocate(self%rsl(size(z_bed_eq,1), size(z_bed_eq,2)), source=0.0_wp)
      allocate(self%z_bed, source=z_bed_eq)

      ! VILMA-v1 backend: hand it the SAME reference and start-slice ice, already on
      ! the model Gauss grid. vilma_v1 owns the second remap leg (model Gauss <->
      ! VILMA-v1's own grid at par%vilma_v1_jmax) and VILMA-v1's whole file-based setup.
      if (self%use_vilma_v1) then
         call vilma_v1_init(self%vilma_v1, self%par, self%sht, &
                            self%gg%z_bed_eq, self%gg%h_ice_eq, self%gg%h_ice)
         self%worst_mass_resid = VILMA_UNSET      ! no analogue in VILMA-v1; see the type
      end if
   end subroutine solid_earth_init

   subroutine solid_earth_check_solver(par)
      !! Validate &vilma solver and, for "v1", that this binary actually HAS the
      !! optional VILMA-v1 backend. Split out of solid_earth_init so a caller can fail
      !! the moment it has read the configuration, before doing any work: a namelist
      !! typo, or asking for VILMA-v1 in a default build, should cost nothing.
      !! vilma_v1_require aborts with an actionable rebuild message — never a link
      !! error and never a crash.
      type(vilma_param_class), intent(in) :: par
      select case (trim(par%solver))
      case ("v2")   ! native solver: nothing to check
      case ("v1"); call vilma_v1_require()
      case default
         write(*,'(a)') ' vilma: unknown solver "'//trim(par%solver)//'"'
         write(*,'(a)') '   &vilma solver must be "v2" (native, the default) or "v1"'
         error stop 'unknown &vilma solver (use v2|v1)'
      end select
   end subroutine solid_earth_check_solver

   subroutine build_solver(self, np, nl)
      !! Build the earth structure and the sub-solvers (response, SLE, adaptive
      !! stepper, optional rotation) from self%par on the model's Gauss grid.
      type(solid_earth), intent(inout) :: self
      integer,           intent(in)    :: np, nl
      real(wp) :: dt0

      self%earth = build_earth(self%par)

      ! viscoelastic driver: operators assembled + factored once, memory zeroed.
      ! Δt enters only through Mk = (μ/η)Δt, which the adaptive stepper rescales
      ! per sub-step (resp%set_dt) — so the init Δt is just a nominal seed (the host
      ! passes the real interval to solid_earth_update; this only seeds the assembly).
      dt0 = self%par%dt_init;  if (dt0 <= 0.0_wp) dt0 = kyr
      ! Degree-1 frame: read BEFORE the response is built, because the elastic and
      ! viscoelastic gains are computed at init and the frame is part of them.
      self%resp%deg1_cm = (trim(self%par%deg1_frame) == "cm")
      select case (trim(self%par%earth_response))
      case ("ve")
         call response_init_ve(self%resp, self%earth, self%sht, dt0)
         self%resp%scheme          = scheme_from_name(self%par%scheme)
         self%resp%max_couple_iter = self%par%max_couple_iter
      case ("elastic")
         call response_init_elastic(self%resp, self%earth, self%sht%lmax)
      case ("null")
         call response_init_null(self%resp)
      case default
         error stop "solid_earth_init: unknown earth_response (use ve|elastic|null)"
      end select
      self%resp%visc3d_tol = self%par%visc3d_tol   ! 3-D split threshold (read before any enable below)
      self%resp%toroidal   = self%par%l_toroidal   ! toroidal coupling (likewise)

      ! laterally-varying (3D) viscosity (rung 6c), per self%par%l_visc_3d
      if (self%par%l_visc_3d) call solid_earth_enable_visc_3d(self, self%sht)

      ! sea-level equation knobs. Warm-start the fixed point across steps: gg%rsl
      ! persists (seeded to 0), and between adjacent steps the coastline barely moves,
      ! so the previous solution is a near-converged seed — sharply cutting the inner
      ! iteration count (the dominant per-step cost is the SLE's spherical transforms).
      self%sle%n_outer      = self%par%sle_n_outer
      self%sle%n_inner      = self%par%sle_n_inner
      self%sle%tol          = self%par%sle_tol
      self%sle%max_mem_iter = self%par%sle_max_mem_iter
      self%sle%fixed_ocean  = self%par%sle_fixed_ocean
      self%sle%subgrid      = self%par%sle_subgrid
      self%sle%warm_start   = .true.

      ! adaptive-Δt controller (vilma_timestep)
      self%stepper%rtol       = self%par%rtol
      self%stepper%atol       = self%par%atol
      self%stepper%safety     = self%par%safety
      self%stepper%grow_max   = self%par%grow_max
      self%stepper%shrink_min = self%par%shrink_min
      self%stepper%dt_min     = self%par%dt_min
      self%stepper%dt_max     = self%par%dt_max
      self%stepper%cfl        = self%par%cfl           ! explicit (fe) sub-step Maxwell ceiling
      self%stepper%dt_try     = self%par%dt_init       ! 0 => first guess = whole interval

      ! rotational feedback (degree-2 Liouville polar motion → centrifugal potential
      ! fed back into the SLE; vilma_rotation).
      self%rotation%enabled = self%par%rotation
      if (self%par%rotation) then
         if (self%par%rotation_k_s > 0.0_wp) then
            call rotation_init(self%rotation, self%earth, self%sht, dt0, &
                               k_s=self%par%rotation_k_s, CminusA=self%par%rotation_c_minus_a)
         else
            call rotation_init(self%rotation, self%earth, self%sht, dt0, &
                               CminusA=self%par%rotation_c_minus_a)
         end if
         self%rotation%enabled = .true.          ! init clears it; turn back on
         ! the channels relax on the response's radial viscosity, not the layer table
         call rotation_set_rate(self%rotation, response_radial_rate(self%resp))
      end if
   end subroutine build_solver

   subroutine solid_earth_enable_visc_3d(self, sht)
      !! Read the lon-lat-r log10(eta) field (load_visc_3d) and enable the tensor-
      !! correct lateral-viscosity (3-D) memory advance. The Maxwell memory state is
      !! preserved, so this can be called either at init or AFTER a 1-D spin-up to
      !! switch the transient onto the 3-D path from the 1-D seed.
      type(solid_earth),    intent(inout)       :: self
      type(sht_grid),       intent(in), target  :: sht
      real(wp), allocatable :: visc_node(:,:)
      call load_visc_3d(self%par, sht, self%resp%r, visc_node)
      if (self%resp%kind /= RESP_VE) &
         error stop 'solid_earth_enable_visc_3d: lateral viscosity needs earth_response=ve'
      call apply_visc_nodes(self, sht, visc_node)
   end subroutine solid_earth_enable_visc_3d

   subroutine apply_visc_nodes(self, sht, visc_node)
      !! Hand a node log10(η) field to the response, and the resulting radial profile
      !! to the rotation channels, so the two always relax on the same Earth. Every
      !! viscosity change goes through here (the 3-D enable and the 1-D pre-spin).
      !! Before build_solver the rotation is not yet initialised; build_solver syncs it.
      type(solid_earth), intent(inout) :: self
      type(sht_grid),    intent(in)    :: sht
      real(wp),          intent(in)    :: visc_node(:,:)
      call response_enable_lateral_visc_from_nodes(self%resp, sht, visc_node, &
           lid_depth=self%par%visc3d_lid_depth, lid_log10max=self%par%visc3d_lid_log10max, &
           log10_cap=self%par%visc_log10_max)
      if (self%rotation%enabled) call rotation_set_rate(self%rotation, response_radial_rate(self%resp))
   end subroutine apply_visc_nodes

   subroutine solid_earth_update(self, h_ice, dt_yr)
      !! Advance the model from time to time+dt_yr [years] under the ice thickness
      !! h_ice (host grid) and store the results in the derived type: se%rsl, se%z_bed
      !! (host grid) and se%gg%* (Gauss grid). The ice load is ramped linearly from the
      !! previous h_ice to the new one across the interval; the adaptive stepper
      !! (vilma_timestep) chooses the internal Δt, solving the SLE against the current
      !! relaxation state and advancing the Maxwell memory at each sub-step.
      type(solid_earth), intent(inout) :: self
      real(wp),           intent(in)    :: h_ice(:,:)   !! grounded-ice thickness [m] (host grid)
      real(wp),           intent(in)    :: dt_yr        !! interval to advance [years]
      real(wp), allocatable :: ice_new(:,:)
      real(wp) :: dt, t0, t1
      integer  :: np, nl
      integer(kind=8) :: pc0, pc1, prate            ! PROFILE: see %t_remap

      np = self%sht%nphi;  nl = self%sht%nlat
      dt = dt_yr*sec_per_year                       ! interface is years; integrator is seconds
      t0 = self%time*sec_per_year;  t1 = t0 + dt

      ! ice load on the Gauss grid (conservative remap in, or passthrough)
      allocate(ice_new(np,nl))
      call system_clock(pc0, prate)
      call to_gauss(self, h_ice, ice_new, conserve_mass=.true.)
      call system_clock(pc1);  self%t_remap = self%t_remap + real(pc1-pc0,wp)/prate

      if (self%use_vilma_v1) then
         ! --- VILMA-v1 backend ------------------------------------------------------
         ! Same contract, different solver: ice in, relative sea level out. vilma_v1
         ! remaps ice_new (model Gauss grid) onto VILMA-v1's own grid, advances VILMA-v1
         ! over [time, time+dt_yr], and maps its rsl back onto the model Gauss grid,
         ! so everything below — and every consumer of se%gg / se%rsl / se%z_bed — is
         ! on exactly the grids it is for the native solver.
         call vilma_v1_update(self%vilma_v1, self%sht, ice_new, dt_yr, &
                              self%gg%rsl, t_remap=self%t_remap, t_solve=self%t_solver)
         self%gg%h_ice = ice_new
         self%gg%z_bed = self%gg%z_bed_eq - self%gg%rsl
         self%time     = self%time + dt_yr

         ! Diagnostics VILMA-v1 does not hand back. DERIVED HONESTLY from VILMA-v1's own
         ! output, using the SAME formulas the native solver uses, so the two
         ! backends' diagnostics are like-for-like:
         !   C   — flotation (vilma_sle ocean_function) applied to VILMA-v1's updated bed
         !         and the current ice. This is VILMA's diagnostic of VILMA-v1's
         !         state, NOT VILMA-v1's internal ocean function (which the library does
         !         not expose on this grid); it can differ from VILMA-v1's own coastline
         !         by a grid cell where the two flotation rules disagree.
         !   bsl — update_bsl's barystatic integral over that C and the ice anomaly.
         !         Again a VILMA diagnostic of VILMA-v1's state, not VILMA-v1's own
         !         ocean bookkeeping (vega_oce.dat).
         call ocean_function(self%gg%z_bed, self%gg%h_ice, self%gg%C)
         call update_bsl(self)
         ! NOT derivable: VILMA-v1 reports no sea-level-equation mass residual, so this
         ! stays at the documented fill value rather than a made-up number. The
         ! driver skips printing it (see vilma_drive).
         self%worst_mass_resid = VILMA_UNSET

         self%h_ice = h_ice
         call system_clock(pc0)
         call solid_earth_sync_host(self)
         call system_clock(pc1);  self%t_remap = self%t_remap + real(pc1-pc0,wp)/prate
         return
      end if

      if (self%rotation%enabled) then
         ! Rotational feedback, stepped with the solid Earth: every SLE solve in the
         ! stepper converges the rotation ↔ sea-level fixed point at the end of its own
         ! sub-step (vilma_rotation: rotation_open_step / trial / close_step).
         call stepper_advance(self%stepper, self%sht, self%resp, self%sle, self%gg%z_bed_eq, &
                                   self%gg%h_ice, ice_new, self%gg%h_ice_eq, t0, t1, &
                                   self%gg%rsl, self%gg%C, rot=self%rotation, sigma_out=self%sigma_lm)
      else
         call stepper_advance(self%stepper, self%sht, self%resp, self%sle, self%gg%z_bed_eq, &
                                   self%gg%h_ice, ice_new, self%gg%h_ice_eq, t0, t1, &
                                   self%gg%rsl, self%gg%C, sigma_out=self%sigma_lm)
      end if

      self%gg%h_ice         = ice_new
      self%gg%z_bed         = self%gg%z_bed_eq - self%gg%rsl
      self%worst_mass_resid = self%stepper%worst_mass_resid
      self%time             = self%time + dt_yr
      call update_bsl(self)

      self%h_ice = h_ice
      call system_clock(pc0)
      call solid_earth_sync_host(self)
      call system_clock(pc1);  self%t_remap = self%t_remap + real(pc1-pc0,wp)/prate
   end subroutine solid_earth_update

   subroutine solid_earth_sync_host(self)
      !! Refresh the host-grid output fields from the Gauss-grid state: map the smooth
      !! rsl perturbation back (bilinear, or copy when passthrough) and reconstruct
      !! z_bed on the host's own (high-resolution) bed so its detail survives. update()
      !! calls this each step; vilma_restart_read calls it after restoring the gg state.
      type(solid_earth), intent(inout) :: self
      if (self%remap) then
         call remap_to_ll(self%map, self%gg%rsl, self%rsl)
      else
         self%rsl = self%gg%rsl
      end if
      self%z_bed = self%z_bed_eq - self%rsl
   end subroutine solid_earth_sync_host

   subroutine solid_earth_spinup(self, h_ice_lgm)
      !! Bring the model to isostatic equilibrium under the start-slice (LGM) ice load,
      !! HOLDING the reference (z_bed_eq, h_ice_eq) as the datum, so the transient that
      !! follows measures rsl/bsl against the reference. Driven entirely by self%par:
      !!
      !!   pre_spinup_1d : run a cheap 1-D pre-equilibration first (the 1-D radial
      !!                   viscosity is the lateral geometric mean of the 3-D field),
      !!                   then switch to the full model carrying the spun-up memory.
      !!   equil_time_max: relax the FULL model, exiting early when the bed stops moving
      !!                   (rate < equil_rate_tol) or at this cap [years] with a warning
      !!                   if not yet converged.
      !!
      !! Both phases relax to bed-stationary convergence; the relaxation interval is
      !! chosen internally (no user dt). A no-op if both are disabled. On return the
      !! model is always in its full configuration, ready for the transient.
      type(solid_earth), intent(inout) :: self
      real(wp),           intent(in)    :: h_ice_lgm(:,:)   !! start-slice (LGM) ice [m] (host grid)
      real(wp), allocatable :: visc_node(:,:), visc_unif(:,:)
      integer :: r, nh
      logical :: pre_1d
      real(wp) :: t_max

      pre_1d = self%par%pre_spinup_1d
      t_max  = self%par%equil_time_max/sec_per_year          ! par stores SI; relax works in years
      if (.not. pre_1d .and. t_max <= 0.0_wp) return         ! nothing to do

      ! KNOWN GAP (documented, not faked): this spin-up relaxes the model's own
      ! viscous memory while HOLDING the reference as the datum. VILMA-v1 owns its
      ! memory internally and advances it only along its own time axis, so the
      ! reference-held relaxation cannot be reproduced without silently consuming
      ! VILMA-v1's clock and desynchronising it from the forcing. Refuse plainly
      ! rather than run something that looks like a spin-up but is not.
      if (self%use_vilma_v1) then
         write(*,'(a)') ' solid_earth_spinup: not available with solver="v1".'
         write(*,'(a)') '   VILMA-v1 advances its own viscous memory along its own time axis;'
         write(*,'(a)') '   the reference-held relaxation this routine performs has no VILMA-v1'
         write(*,'(a)') '   analogue. Set equil_time_max=0 and pre_spinup_1d=.false. in &vilma,'
         write(*,'(a)') '   and start the transient from the full (e.g. LGM->present) window,'
         write(*,'(a)') '   which is how CLIMBER-X drives VILMA-v1. See doc/vilma-v1-backend.md.'
         error stop 'solid_earth_spinup: unsupported with solver="v1" (see message above)'
      end if

      call solid_earth_update(self, h_ice_lgm, 0.0_wp)       ! seed the entering ice = start (LGM) ice

      if (pre_1d) then
         ! 1-D phase: replace the lateral viscosity with its lateral geometric mean
         ! (a laterally-uniform field => a 1-D model on the mean radial profile). The
         ! Maxwell memory is preserved across the enable, so it seeds the full model.
         if (self%par%l_visc_3d) then
            call load_visc_3d(self%par, self%sht, self%resp%r, visc_node)
            nh = size(visc_node, 1)
            allocate(visc_unif(nh, size(visc_node, 2)))
            do r = 1, size(visc_node, 2)
               visc_unif(:, r) = sum(visc_node(:, r)) / real(nh, wp)   ! lateral mean of log10(eta)
            end do
            call apply_visc_nodes(self, self%sht, visc_unif)
         end if
         call relax_hold(self, h_ice_lgm, -1.0_wp, "1-D ")    ! converge; internal pass cap only
         if (self%par%l_visc_3d) call solid_earth_enable_visc_3d(self, self%sht)  ! restore the real 3-D field
      end if

      if (t_max > 0.0_wp) call relax_hold(self, h_ice_lgm, t_max, "full")
   end subroutine solid_earth_spinup

   subroutine relax_hold(self, h_ice_lgm, t_cap, label)
      !! Hold h_ice_lgm and relax to bed-stationary convergence, advancing in an
      !! internally-grown interval (the adaptive stepper sub-steps within it). t_cap<0:
      !! no time cap (internal pass cap only, for the fast 1-D phase). t_cap>0: a cap
      !! [years] — exit with a WARNING + diagnostics if convergence is not reached.
      !! Logs per pass: wall time, mean internal Δt, sub-steps, bed residual, resid/tol.
      type(solid_earth), intent(inout) :: self
      real(wp),           intent(in)    :: h_ice_lgm(:,:)
      real(wp),           intent(in)    :: t_cap
      character(len=*),   intent(in)    :: label
      real(wp), allocatable :: z_prev(:,:)
      real(wp) :: pass_dt, t_done, rmean, rmax, rrate, wall0, wall1, mean_dt, tol
      integer  :: it, nsub0, nsub, np, nl
      logical  :: converged

      tol = self%par%equil_rate_tol
      np = self%sht%nphi;  nl = self%sht%nlat
      allocate(z_prev(np,nl), source=self%gg%z_bed)
      pass_dt = SPINUP_DT0;  t_done = 0.0_wp;  converged = .false.
      rmean = huge(1.0_wp);  rrate = huge(1.0_wp)
      if (t_cap > 0.0_wp) then
         write(*,'(a,a,a,f0.0,a)') ' spin-up [', trim(label), &
              ']: relax LGM memory vs reference, cap ', t_cap, ' yr'
      else
         write(*,'(a,a,a)') ' spin-up [', trim(label), &
              ']: relax LGM memory vs reference (1-D pre-spin)'
      end if
      do it = 1, SPINUP_MAXP
         if (t_cap > 0.0_wp) pass_dt = min(pass_dt, t_cap - t_done)
         if (pass_dt <= 0.0_wp) exit
         nsub0 = self%stepper%n_accept
         call cpu_time(wall0)
         call solid_earth_update(self, h_ice_lgm, pass_dt)
         call cpu_time(wall1)
         t_done  = t_done + pass_dt
         nsub    = max(1, self%stepper%n_accept - nsub0)
         mean_dt = pass_dt / real(nsub, wp)
         rmean   = sht_grid_surface_integral(self%sht, abs(self%gg%z_bed - z_prev)) / (4.0_wp*pi)
         rmax    = maxval(abs(self%gg%z_bed - z_prev))
         rrate   = rmean / pass_dt                            ! mean bed velocity over the pass [m/yr]
         write(*,'(a,i3,a,f7.2,a,f9.1,a,i0,a,f9.4,a,f9.2,a,es9.2,a,es8.1)') &
              '   pass ', it, ':  ', wall1-wall0, ' s  <dt>=', mean_dt, ' yr (', nsub, &
              ' sub)  d|z_bed|=', rmean, ' m  max=', rmax, ' m  v=', rrate, ' m/yr  v/tol=', rrate/tol
         if (rrate < tol) then;  converged = .true.;  exit;  end if
         z_prev  = self%gg%z_bed
         pass_dt = pass_dt*SPINUP_GROW
      end do
      if (t_cap > 0.0_wp .and. .not. converged) &
         write(*,'(a,a,a,f0.1,a,es9.2,a,es8.1,a)') ' WARNING: spin-up [', trim(label), &
              '] hit equil_time_max=', t_done, ' yr before convergence (v=', &
              rrate, ' m/yr, ', rrate/tol, 'x tol)'
   end subroutine relax_hold

   subroutine update_bsl(self)
      !! Diagnose the barystatic sea level from the current (Gauss-grid) state.
      !!
      !! The grounded-ice increment is masked exactly as vilma_sle masks it: the
      !! current column by the current ocean function, the REFERENCE column by
      !! the reference one. Masking the raw difference by the current C alone
      !! drops the whole column of any cell that carried grounded marine
      !! reference ice and has since flooded -- see the ΔI_g note in vilma_sle.
      !!
      !! This is the WHOLE grounded-column change, NOT the sea level vilma_sle
      !! delivers: there is no subgrid term here, while vilma_sle's mass-conservation
      !! offset carries one. Over cells that flooded since the reference, bsl
      !! therefore reads high by their below-flotation volume -- about 10 m on the
      !! LGM-datum deglaciation. It is an ice-volume-equivalent diagnostic, not the
      !! barystatic rise; for that, take the ocean-mean of C*rsl.
      type(solid_earth), intent(inout) :: self
      real(wp) :: c_int
      c_int = sht_grid_surface_integral(self%sht, self%gg%C)
      if (c_int > 0.0_wp) then
         self%bsl = -(rho_ice/rho_water) * sht_grid_surface_integral(self%sht, &
                       self%gg%h_ice   *(1.0_wp - self%gg%C) &
                     - self%gg%h_ice_eq*(1.0_wp - self%gg%C0)) / c_int
      else
         self%bsl = 0.0_wp
      end if
   end subroutine update_bsl

   subroutine solid_earth_finalize(self)
      type(solid_earth), intent(inout) :: self
      if (self%use_vilma_v1) call vilma_v1_finalize(self%vilma_v1)
      self%use_vilma_v1 = .false.
      call response_destroy(self%resp)
      call rotation_destroy(self%rotation)
      if (associated(self%sht)) then
         call sht_grid_destroy(self%sht)        ! model-owned grid
         deallocate(self%sht)
         self%sht => null()
      end if
      self%remap = .false.
      if (allocated(self%sigma_lm))    deallocate(self%sigma_lm)
      if (allocated(self%gg%z_bed_eq)) deallocate(self%gg%z_bed_eq)
      if (allocated(self%gg%h_ice_eq)) deallocate(self%gg%h_ice_eq)
      if (allocated(self%gg%h_ice))    deallocate(self%gg%h_ice)
      if (allocated(self%gg%rsl))      deallocate(self%gg%rsl)
      if (allocated(self%gg%z_bed))    deallocate(self%gg%z_bed)
      if (allocated(self%gg%C))        deallocate(self%gg%C)
      if (allocated(self%gg%C0))       deallocate(self%gg%C0)
      if (allocated(self%z_bed_eq))    deallocate(self%z_bed_eq)
      if (allocated(self%h_ice_eq))    deallocate(self%h_ice_eq)
      if (allocated(self%h_ice))       deallocate(self%h_ice)
      if (allocated(self%rsl))         deallocate(self%rsl)
      if (allocated(self%z_bed))       deallocate(self%z_bed)
      self%time = 0.0_wp
   end subroutine solid_earth_finalize

   ! --- internals --------------------------------------------------------------

   subroutine to_gauss(self, f_host, f_gauss, conserve_mass)
      !! Bring a host-grid field onto the Gauss grid: conservative remap when a host
      !! grid is in play, else a straight copy (passthrough).
      type(solid_earth), intent(in)  :: self
      real(wp),          intent(in)  :: f_host(:,:)
      real(wp),          intent(out) :: f_gauss(:,:)
      logical,           intent(in)  :: conserve_mass
      if (self%remap) then
         call remap_to_gauss(self%map, self%sht, f_host, f_gauss, conserve_mass=conserve_mass)
      else
         f_gauss = f_host
      end if
   end subroutine to_gauss

   subroutine build_sht(p, sht)
      !! Build the Gauss-Legendre transform grid from the parameter record. When
      !! nlat/nphi are unset (<=0) default to a de-aliased grid (nlat=2 lmax+2,
      !! nphi=4 lmax) sized for the SLE's quadratic ocean-function product.
      type(vilma_param_class), intent(in)    :: p
      type(sht_grid),       intent(inout) :: sht
      integer :: nlat, nphi
      nlat = p%nlat;  if (nlat <= 0) nlat = 2*p%lmax + 2
      nphi = p%nphi;  if (nphi <= 0) nphi = 4*p%lmax
      if (p%mmax >= 0 .and. p%eps_polar > 0.0_wp) then
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mmax=p%mmax, mres=p%mres, eps_polar=p%eps_polar)
      else if (p%mmax >= 0) then
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mmax=p%mmax, mres=p%mres)
      else if (p%eps_polar > 0.0_wp) then
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mres=p%mres, eps_polar=p%eps_polar)
      else
         call sht_grid_init(sht, p%lmax, nlat=nlat, nphi=nphi, mres=p%mres)
      end if
   end subroutine build_sht

   logical function grid_is_gauss(grid, sht) result(same)
      !! True when the host grid IS the model's Gauss grid (same dims and matching
      !! lon/lat axes), so no remap is needed and the host fields pass straight through.
      type(grid_class), intent(in) :: grid
      type(sht_grid),   intent(in) :: sht
      real(wp), parameter :: TOL = 1.0e-6_wp           ! degrees
      real(wp), parameter :: RAD2DEG = 57.295779513082323_wp
      integer :: j
      same = .false.
      if (grid%G%nx /= sht%nphi .or. grid%G%ny /= sht%nlat) return
      ! longitudes (SHTns lon is ascending in [0,360))
      do j = 1, sht%nphi
         if (abs(real(grid%G%x(j),wp) - sht%lon(j)*RAD2DEG) > TOL) return
      end do
      ! latitudes: SHTns colat is north-first (descending lat); the host axis is
      ! ascending, so compare against the reversed Gauss colat row.
      do j = 1, sht%nlat
         if (abs(real(grid%G%y(j),wp) - (90.0_wp - sht%colat(sht%nlat - j + 1)*RAD2DEG)) > TOL) return
      end do
      same = .true.
   end function grid_is_gauss

end module vilma_coupling
