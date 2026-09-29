module vilma_sle
   !! Sea-level equation: gravitationally self-consistent, mass-conserving
   !! redistribution of ocean water over a deforming solid Earth and geoid, with
   !! migrating coastlines (Kendall, Mitrovica & Milne 2005; Martinec et al.
   !! 2018 benchmark).
   !!
   !! Solved pseudo-spectrally: the load → (uplift, geoid) convolution is done in
   !! spectral space (via a response_operator, vilma_response), while the ocean
   !! function multiply C·S is a pointwise product on the spatial Gauss grid —
   !! this is what avoids Gibbs ringing at coastlines.
   !!
   !! The equation is a Fredholm equation of the second kind, so it is iterated.
   !! For a fixed coastline the change in relative sea level over the ocean is
   !!
   !!     S = C·( N − u + Δφ ) ,
   !!
   !! where N (geoid) and u (uplift) are the response to the TOTAL surface load
   !! L = ρ_i ΔI_g + ρ_w (C·S), and the spatial constant Δφ (a uniform shift of
   !! the equipotential) is fixed each iteration by ocean-mass conservation,
   !!
   !!     ρ_w ∫ C·S dA = −ρ_i ∫ ΔI_g dA   (melt water volume) ,
   !!     Δφ = [ −(ρ_i/ρ_w) ∫ΔI_g dΩ − ∫ C(N−u) dΩ ] / ∫ C dΩ .
   !!
   !! ΔI_g is the GROUNDED-ice increment, I·(1−C) − I⁽⁰⁾·(1−C⁽⁰⁾): each endpoint
   !! masked by its own ocean function, because only grounded ice loads the bed
   !! and only grounded ice exchanges mass with the ocean. See sle_solve.
   !!
   !! Because Δφ is built to satisfy that balance, mass is conserved to machine
   !! precision at every iteration. The inner loop iterates S (the water load
   !! feeds back through the response); the outer loop rebuilds the ocean
   !! function C from the migrated topography topo0 − S (moving shorelines).
   use vilma_precision, only: wp
   use vilma_constants, only: rho_ice, rho_water
   use vilma_sht,       only: sht_grid, sht_grid_surface_integral, sht_grid_analysis, sht_grid_synthesis
   use vilma_rotation,  only: rotation_state, rotation_open_step, rotation_trial, rotation_close_step, &
                              rotation_s_rot
   use vilma_response,  only: response_finalize_step, response_endpoint_converged, response_advance_endpoint, response_apply, response_prepare_endpoint, response_begin_step, response, response_init_elastic, response_init_ve, response_init_null
   implicit none
   private

   public :: sle_solver, sle_result, ocean_function
   public :: sle_solve

   type :: sle_result
      !! Diagnostics returned by a solve.
      integer  :: n_outer_done = 0      !! coastline iterations performed
      integer  :: n_couple_done = 0     !! SLE<->memory co-convergence passes (§3c 3b)
      integer  :: n_inner_last = 0      !! inner iterations in the last outer pass
      integer  :: n_coast_flip = 0      !! cells that changed ocean/land state in the last pass
                                        !! (0 = the coastline converged; see sle_solve)
      real(wp) :: resid        = 0.0_wp !! last inner max|ΔS| [m]
      real(wp) :: mass_resid   = 0.0_wp !! relative ocean-mass-conservation error
      real(wp) :: ocean_frac   = 0.0_wp !! ∫C dΩ / 4π
      real(wp) :: esl          = 0.0_wp !! eustatic offset Δφ [m] (uniform sea-surface shift)
      real(wp), allocatable :: u(:,:)   !! converged solid uplift [m] (nphi,nlat)
      real(wp), allocatable :: N(:,:)   !! converged geoid rise [m] (nphi,nlat)
   end type sle_result

   type :: sle_solver
      !! Cap on the paleotopography / coastline iterations, not a fixed count: the
      !! loop exits as soon as the migrated coastline stops moving (see sle_solve).
      integer  :: n_outer = 3          !! paleotopography / coastline iterations
      integer  :: n_inner = 20         !! water-load fixed-point iterations
      real(wp) :: tol     = 1.0e-7_wp  !! inner convergence on max|ΔS| / max|S| [-]
      !! SLE<->memory co-convergence cap (§3c part 3b). When the response carries
      !! viscoelastic memory under an implicit (trapezoidal) scheme and the load
      !! evolves fast within a step, the ocean load σ and the end-of-step memory
      !! τ_{n+1} are a mutual fixed point: each pass re-converges the water load
      !! against the latest τ_{n+1} estimate, then advances the memory one trapezoid
      !! pass. The response signals convergence (surface drift settled) to exit
      !! early; 1st-order / stateless responses converge in a single pass regardless,
      !! so this is inert for them (and for the FE default). > 1 enables the 2nd-order
      !! co-convergence for the trapezoidal scheme.
      integer  :: max_mem_iter = 20
      !! Ocean geometry. .false. (default) = time-varying: the coastline migrates
      !! each outer pass from the deformed surface topo0 − rsl (Martinec 2018 §2.2,
      !! the SLE2 suite). .true. = fixed: the ocean function is held at the initial
      !! O⁽⁰⁾ = (topo0 < 0) throughout (§2.1, eq 1, the SLE1 suite) — no coastline
      !! migration, so a single outer pass converges the inner water load.
      logical  :: fixed_ocean = .false.
      !! Sloping-coast ("subgrid") ocean water load. .true. (default) uses the
      !! actual water-column change (Martinec 2018 §2.2, eqs 15-19): s = C·rsl −
      !! ζ⁽⁰⁾·(C − C⁽⁰⁾), which equals rsl over permanent ocean but rsl − ζ⁽⁰⁾ over
      !! newly-flooded cells, tapering to zero at the coast where the bed meets the
      !! sea surface — the mass-correct load on a migrating, sloping coastline.
      !! .false. is the simpler binary load ρ_w·C·rsl (the full sea-level change
      !! wherever C = 1, a sharp coastline), kept for comparison; it agrees with
      !! subgrid when the coastline does not move (deep basins, fixed_ocean).
      logical  :: subgrid = .true.
      !! Warm start. .false. (default) zeroes rsl at the top of every solve (a
      !! cold start — what the benchmarks expect). .true. reuses the incoming rsl
      !! as the initial guess for the fixed point. For transient runs the coastline
      !! and RSL change very little between adjacent time steps, so the previous
      !! step's converged solution is a near-converged seed and the inner loop
      !! exits in far fewer iterations (often 1–2). The converged fixed point is
      !! unique (the SLE is a Fredholm equation of the 2nd kind → a contraction)
      !! and mass is rebalanced every iteration via Δφ, so the ANSWER is unchanged
      !! to tolerance — only the iteration count drops. The caller must keep rsl
      !! alive across calls (the coupling driver does; it turns this on at init).
      logical  :: warm_start = .false.

      !! --- PROFILE: wall-clock accumulators [s], and iteration counters ------
      !! The driver reports solid_earth_update as drift + memory + a RESIDUAL
      !! "SLE + coupling (rest)" bucket. These decompose the SLE's share of that
      !! bucket, so the split between harmonic transforms and grid-space work is
      !! measured rather than inferred:
      !!
      !!     t_total = t_sht + t_apply + t_resp + (grid-space remainder)
      !!
      !! t_resp is the nested response lifecycle (begin/prepare/advance/finalize)
      !! and is ALSO accumulated in response%t_drift / %t_mem — it is the overlap
      !! with the driver's breakdown, so subtract it before comparing. Everything
      !! else here is work no other timer sees. Cost is one system_clock pair per
      !! region (~25 ns); negligible against a transform.
      real(wp) :: t_total = 0.0_wp  !! whole sle_solve
      real(wp) :: t_sht   = 0.0_wp  !! sht_grid_analysis + sht_grid_synthesis
      real(wp) :: t_apply = 0.0_wp  !! response_apply (spectral load -> u, N)
      real(wp) :: t_resp  = 0.0_wp  !! response lifecycle (overlaps t_drift/t_mem)
      real(wp) :: t_rot   = 0.0_wp  !! rotation trial + s_rot + commit (polar motion)
      integer  :: n_solve     = 0   !! sle_solve calls
      integer  :: n_outer_tot = 0   !! coastline passes, summed over calls
      integer  :: n_inner_tot = 0   !! inner water-load iterations, summed
   end type sle_solver

contains

   subroutine sle_solve(self, sht, resp, d_ice, ice, topo0, rsl, C, res, &
                        report_only, sigma_lm, rot, rot_dt)
      !! Solve for the relative-sea-level change rsl [m] driven by a grounded-ice
      !! thickness change d_ice [m], on a reference topography topo0 [m] (solid
      !! surface relative to the reference sea surface; ocean where < 0).
      !!
      !! rsl is the FULL relative-sea-level change field, N − u + Δφ (geoid rise
      !! minus solid uplift, plus the mass-conservation offset), defined over the
      !! WHOLE grid — not just the ocean. On land it is the bedrock-vs-sea-surface
      !! change that drives bedrock motion under grounded ice; over the ocean it
      !! is the sea-level change. The ocean-masked sea level is simply C·rsl, so
      !! it is not returned separately. The bedrock relative to the sea surface is
      !! topo0 − rsl everywhere.
      !!
      !! The absolute grounded-ice thickness ice [m] is passed alongside the
      !! change d_ice: both are needed because the load is incremental but
      !! flotation is an absolute condition. ice enters the coastline test in
      !! ocean_function so that ice thick enough to ground on the bed is excluded
      !! from the ocean (it bears on the solid surface, it does not float — even
      !! where the bed has subsided below the sea surface), and ice together with
      !! the reference ice − d_ice gives the two absolute columns whose difference,
      !! each masked by its own coastline, is the grounded-ice load increment ΔI_g
      !! that actually drives the load and the melt source (see ΔI_g below). The
      !! raw d_ice is never used unmasked.
      type(sle_solver),        intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      type(response), intent(inout) :: resp
      real(wp),                 intent(in)    :: d_ice(:,:)  !! ice CHANGE [m] (load)
      real(wp),                 intent(in)    :: ice(:,:)    !! abs. ice [m] (flotation)
      real(wp),                 intent(in)    :: topo0(:,:)  !! (nphi,nlat) [m]
      !! intent(inout): zeroed here on a cold start, or reused as the initial
      !! guess when self%warm_start (see the type definition).
      real(wp),                 intent(inout) :: rsl(:,:)    !! full RSL change [m]
      real(wp),                 intent(out)   :: C(:,:)      !! (nphi,nlat) ocean fn
      type(sle_result),         intent(out)   :: res
      !! report_only (default .false.): converge the load against the FROZEN entering
      !! memory and return WITHOUT advancing the memory or time (no co-convergence) —
      !! a pure "what is the load/response at the current state" query. Used to seed the
      !! trapezoidal start-of-step load σ_0 at t=0 (where the memory is at rest, so the
      !! load is the elastic-consistent one). sigma_lm (optional out): the converged
      !! spectral surface load, in either mode.
      logical,          optional, intent(in)  :: report_only
      complex(wp),      optional, intent(out) :: sigma_lm(:)
      !! rot, rot_dt (optional, together): rotational feedback, stepped with the solid
      !! Earth over this solve's step of length rot_dt [s]. Each inner iteration's load
      !! sets the polar motion at the end of the step (rotation_trial), and its
      !! rotational sea level s_rot = N_rot − u_rot (geoid minus uplift from the
      !! centrifugal potential) enters that iteration's sea surface, so the rotation ↔
      !! SLE fixed point converges inside this solve. s_rot enters the sea-surface
      !! geometry (Sraw) but NOT the surface mass load that drives the load response:
      !! the rotational potential forces the Earth through vilma_rotation's own tidal
      !! channel, not as a surface mass. Mass is still conserved — Δφ is recomputed
      !! from Sraw including s_rot. The rotation is committed with the closing load
      !! where the response's memory advances; a report-only solve steps it by 0 and
      !! commits nothing. With rot absent (or disabled) the solve is bit-for-bit the
      !! no-rotation result.
      type(rotation_state), optional, intent(inout) :: rot
      real(wp),             optional, intent(in)    :: rot_dt

      real(wp), allocatable :: load(:,:), u(:,:), N(:,:), Sraw(:,:), rsl_new(:,:), srot(:,:)
      real(wp), allocatable :: C0(:,:), wcorr(:,:), C_next(:,:), d_ice_g(:,:)
      complex(wp), allocatable :: load_lm(:), u_lm(:), N_lm(:)
      real(wp) :: rho_ratio, ice_int, dphi, C_int, Cs_int, zeta_int, smax, dmax
      integer  :: im, io, ii, np, nl, n_mem
      logical  :: ronly, rotate
      integer(kind=8) :: pc0, pc1, pca, pcb, prate   ! PROFILE: see %t_total

      call system_clock(pc0, prate)
      self%n_solve = self%n_solve + 1

      np = sht%nphi;  nl = sht%nlat
      allocate(load(np,nl), u(np,nl), N(np,nl), Sraw(np,nl), rsl_new(np,nl))
      allocate(C0(np,nl), wcorr(np,nl), C_next(np,nl), d_ice_g(np,nl))
      allocate(load_lm(sht%nlm), u_lm(sht%nlm), N_lm(sht%nlm))

      rho_ratio = rho_ice/rho_water
      ! Cold start zeroes rsl; warm start keeps the incoming field as the initial
      ! guess (the caller's previous converged solution — see %warm_start).
      if (.not. self%warm_start) rsl = 0.0_wp
      u = 0.0_wp;  N = 0.0_wp;  dphi = 0.0_wp;  ice_int = 0.0_wp
      d_ice_g = 0.0_wp
      zeta_int = 0.0_wp;  wcorr = 0.0_wp
      res%n_inner_last = 0;  res%resid = 0.0_wp;  res%n_outer_done = 0
      res%n_coast_flip = 0

      ! Initial ocean function O⁽⁰⁾ — the reference (t0) coastline against which the
      ! subgrid term measures newly flooded / emerged cells, and the held coastline
      ! in fixed-ocean mode. It is the flotation-aware ocean function of the
      ! reference state: bathymetry ζ⁽⁰⁾ = topo0 (rsl = 0 at t0) and the reference
      ! ice ice − d_ice (= 0 for an ice-free reference ⇒ O⁽⁰⁾ = (topo0 < 0); = ice
      ! when d_ice = 0 ⇒ grounded reference ice is excluded, so a no-change solve
      ! leaves C ≡ C⁽⁰⁾ and the subgrid term vanishes).
      call ocean_function(topo0, ice - d_ice, C0)

      ronly = .false.;  if (present(report_only)) ronly = report_only
      rotate = .false.
      if (present(rot)) rotate = rot%enabled
      if (rotate) then
         if (.not. present(rot_dt)) error stop 'sle_solve: rot needs rot_dt'
         allocate(srot(np,nl))
         call system_clock(pca)
         call rotation_open_step(rot, merge(0.0_wp, rot_dt, ronly))
         call system_clock(pcb);  self%t_rot = self%t_rot + real(pcb-pca,wp)/prate
      end if

      ! Freeze the response's relaxation drift for this time step; for elastic /
      ! null responses this is a no-op.
      call system_clock(pca)
      call response_begin_step(resp, sht)
      ! Open the SLE<->memory co-convergence (§3c 3b): snapshot τ_n. The im loop
      ! re-converges the water load against the latest end-of-step memory estimate
      ! and advances the memory one trapezoid pass each time, until the response's
      ! report drift settles. For FE / elastic / null it runs exactly once (they
      ! report converged after a single advance). In report-only mode there is no
      ! memory advance, so a single load-convergence pass against τ_n suffices.
      if (.not. ronly) call response_prepare_endpoint(resp, sht)
      call system_clock(pcb);  self%t_resp = self%t_resp + real(pcb-pca,wp)/prate
      n_mem = self%max_mem_iter;  if (ronly) n_mem = 1

      do im = 1, n_mem
      do io = 1, self%n_outer
         if (self%fixed_ocean) then
            ! Fixed ocean geometry (Martinec 2018 §2.1): hold the coastline at the
            ! reference O⁽⁰⁾ for all time (no migration). Computed once.
            if (io == 1) C = C0
         else if (io == 1) then
            ! migrate the coastline using the current (full-field) sea level: ocean
            ! where the deformed solid surface topo0 − rsl is below the sea surface
            ! AND the ice there floats rather than grounds (grounded ice keeps a
            ! subsided cell as land).
            call ocean_function(topo0 - rsl, ice, C)
         else
            ! the coastline the previous pass's converged rsl implies; it was built
            ! at the foot of that pass by the convergence test, from this same rsl,
            ! so this is the identical field the old unconditional call produced.
            C = C_next
         end if
         ! GROUNDED-ice thickness increment. Only grounded ice loads the bed and
         ! only grounded ice exchanges mass with the ocean, so the increment that
         ! drives both is the difference of the two grounded columns, each masked
         ! by ITS OWN ocean function — the endpoint's by C, the reference's by
         ! C⁽⁰⁾ — not the raw increment masked by the endpoint's alone:
         !
         !     ΔI_g = I·(1 − C) − I⁽⁰⁾·(1 − C⁽⁰⁾)      [this]
         !     ΔI_g = (I − I⁽⁰⁾)·(1 − C)               [WRONG unless C ≡ C⁽⁰⁾]
         !
         ! The two agree while the cell does not change state, and for an ice-free
         ! reference (I⁽⁰⁾ = 0) they are identical — which is why every benchmark
         ! in the suite, all of which reference an essentially ice-free state,
         ! passes either way. They differ by I⁽⁰⁾·(C − C⁽⁰⁾): a cell carrying
         ! grounded reference ice (C⁽⁰⁾ = 0) that loses it and becomes ocean
         ! (C = 1). The old form masked that cell's whole column out of the melt
         ! source at the moment it flooded, so marine-grounded reference ice
         ! silently delivered no meltwater. Over the last deglaciation referenced
         ! to the LGM that is 45 m of missing barystatic rise, concentrated in
         ! Hudson Bay, the Canadian Arctic, the Baltic and West Antarctica.
         !
         ! The below-flotation part is not double counted: the water that fills
         ! such a newly flooded cell from its bed is removed again by the subgrid
         ! term ζ⁽⁰⁾(C − C⁽⁰⁾) below, leaving exactly the above-flotation volume.
         d_ice_g = ice*(1.0_wp - C) - (ice - d_ice)*(1.0_wp - C0)

         C_int = sht_grid_surface_integral(sht, C)
         if (C_int <= 0.0_wp) exit          ! no ocean: nothing to redistribute

         ! water-equivalent melt source ∝ −(ρ_i/ρ_w)∫ΔI_g dΩ over GROUNDED ice
         ! only: floating ice (C=1) is already in the ocean, so it does not change
         ! the ocean-water budget. Recomputed per coastline pass (grounded set
         ! shifts).
         ice_int = -rho_ratio * sht_grid_surface_integral(sht, d_ice_g)

         ! Subgrid sloping-coast correction (Martinec 2018 eq 17): the ocean water
         ! column change is C·rsl − ζ⁽⁰⁾·(C − C⁽⁰⁾), not C·rsl. The −ζ⁽⁰⁾(C−C⁽⁰⁾)
         ! piece accounts for the bed elevation of cells that crossed the coastline
         ! since t0 (a newly flooded cell fills from its bed, ζ⁽⁰⁾, not from the
         ! reference sea surface), so the load tapers to zero at the moving coast.
         ! Constant over the inner loop (depends only on the coastline C, C⁽⁰⁾). It
         ! also enters the mass balance via ζ̄⁽⁰⁾ = ∫ζ⁽⁰⁾(C−C⁽⁰⁾) (eqs 19-20).
         if (self%subgrid) then
            wcorr    = -rho_water * topo0 * (C - C0)
            zeta_int = sht_grid_surface_integral(sht, topo0*(C - C0))
         else
            wcorr = 0.0_wp;  zeta_int = 0.0_wp
         end if

         do ii = 1, self%n_inner
            call assemble_load()
            if (rotate) then                              ! before analysis overwrites load
               call system_clock(pca)
               call rotation_trial(rot, sht, load)
               call rotation_s_rot(rot, sht, srot)
               call system_clock(pcb);  self%t_rot = self%t_rot + real(pcb-pca,wp)/prate
            end if
            call system_clock(pca)
            call sht_grid_analysis(sht, load, load_lm)            ! analysis overwrites load
            call system_clock(pcb);  self%t_sht = self%t_sht + real(pcb-pca,wp)/prate
            call response_apply(resp, sht, load_lm, u_lm, N_lm)
            call system_clock(pca);  self%t_apply = self%t_apply + real(pca-pcb,wp)/prate
            call sht_grid_synthesis(sht, u_lm, u)
            call sht_grid_synthesis(sht, N_lm, N)
            call system_clock(pcb);  self%t_sht = self%t_sht + real(pcb-pca,wp)/prate

            Sraw = N - u
            if (rotate) Sraw = Sraw + srot              ! rotational feedback
            Cs_int = sht_grid_surface_integral(sht, C*Sraw)
            dphi   = (ice_int - Cs_int + zeta_int)/C_int ! mass-conservation offset
            rsl_new = Sraw + dphi                        ! full field, everywhere

            dmax = maxval(abs(C*(rsl_new - rsl)))        ! converge on the ocean part
            rsl  = rsl_new
            res%n_inner_last = ii;  res%resid = dmax
            smax = maxval(abs(C*rsl))
            if (dmax <= self%tol*max(smax, tiny(1.0_wp))) exit
         end do

         res%n_outer_done = io
         self%n_outer_tot = self%n_outer_tot + 1
         self%n_inner_tot = self%n_inner_tot + res%n_inner_last
         if (self%fixed_ocean) exit         ! C is fixed: one coastline pass converges

         ! Coastline convergence. Migrating C is the outer loop's ONLY job — every
         ! other quantity in a pass (C_int, ice_int, wcorr, zeta_int, and the inner
         ! fixed point itself) is a function of C and the rsl converged against it.
         ! So if the updated rsl implies the same ocean function, the next pass
         ! would reproduce this one exactly: the test below is not a tolerance, it
         ! is the fixed point. Without it the loop always ran n_outer passes, which
         ! for the default n_outer = 3 is up to two redundant ones per solve.
         !
         ! ocean_function assigns the literals 0 and 1, so the comparison is exact
         ! and n_coast_flip counts cells that changed state — a diagnostic worth
         ! having: a count that never reaches 0 means the coastline is oscillating
         ! between two states rather than settling, and n_outer is then a genuine
         ! cap rather than a formality.
         !
         ! C_next is promoted to C at the TOP of the next pass, never here. On the
         ! final pass C must stay the coastline that rsl was converged against,
         ! because the closing load (and wcorr, zeta_int) below are built from both.
         call ocean_function(topo0 - rsl, ice, C_next)
         res%n_coast_flip = count(C_next /= C)
         if (res%n_coast_flip == 0) exit
      end do

      ! Advance the relaxation memory one co-convergence pass with the converged
      ! total load (no-op for elastic / null; one Maxwell update for FE; one
      ! trapezoid endpoint pass for TRAP). Same grounded-ice masking + subgrid
      ! sloping-coast term (wcorr) as the inner load. advance_endpoint also refreshes
      ! the report drift to the new τ_{n+1}, so the next im pass's σ-convergence and
      ! coastline migration see the advanced memory.
      call assemble_load()
      call system_clock(pca)
      if (rotate) call rotation_trial(rot, sht, load)   ! m with the closing load
      call system_clock(pcb);  self%t_rot = self%t_rot + real(pcb-pca,wp)/prate
      call system_clock(pca)
      call sht_grid_analysis(sht, load, load_lm)
      call system_clock(pcb);  self%t_sht = self%t_sht + real(pcb-pca,wp)/prate
      if (ronly) exit                    ! report only: do NOT advance the memory/time
      call response_advance_endpoint(resp, sht, load_lm)
      call system_clock(pca);  self%t_resp = self%t_resp + real(pca-pcb,wp)/prate
      res%n_couple_done = im
      ! Converged when the report drift has settled (the σ<->τ fixed point); 1st-order
      ! / stateless responses report converged after a single pass.
      if (response_endpoint_converged(resp)) exit
      end do

      call system_clock(pca)
      if (.not. ronly) call response_finalize_step(resp, sht)
      if (rotate .and. .not. ronly) call rotation_close_step(rot, sht)
      call system_clock(pca);  self%t_rot = self%t_rot + real(pca-pcb,wp)/prate
      call system_clock(pcb);  self%t_resp = self%t_resp + real(pcb-pca,wp)/prate
      if (present(sigma_lm)) sigma_lm = load_lm   ! converged spectral surface load

      ! diagnostics. The conserved ocean-water volume is ∫s dΩ = ∫C·rsl − ζ̄⁽⁰⁾
      ! (the subgrid sloping-coast term; ζ̄⁽⁰⁾ = 0 in the binary case), which must
      ! balance the melt source ice_int.
      Cs_int = sht_grid_surface_integral(sht, C*rsl) - zeta_int
      res%ocean_frac = C_int/(16.0_wp*atan(1.0_wp))      ! ∫C dΩ / 4π
      if (abs(ice_int) > 0.0_wp) then
         res%mass_resid = abs(Cs_int - ice_int)/abs(ice_int)
      else
         res%mass_resid = abs(Cs_int)
      end if
      res%u = u;  res%N = N;  res%esl = dphi             ! converged fields + offset

      call system_clock(pc1);  self%t_total = self%t_total + real(pc1-pc0,wp)/prate

   contains

      subroutine assemble_load()
         !! Total surface mass load = GROUNDED ice + ocean water + subgrid term.
         !! Ice over ocean cells (C=1: open ocean or floating ice) does not press
         !! its full weight on the bed -- it is borne by buoyancy and carried by
         !! the ocean term ρ_w·C·rsl. That masking is already inside ΔI_g, built
         !! above with each endpoint against its own coastline; without it, ice
         !! overhanging a deep basin over-loads the bed. wcorr is the subgrid
         !! sloping-coast term (zero unless self%subgrid).
         !!
         !! Host-associated: reads d_ice_g, C, rsl and wcorr, writes load — all
         !! of sle_solve's own locals. Called at the head of the inner fixed
         !! point and again for the closing memory-advance load, which MUST see
         !! the same C the rsl was converged against. One definition rather than
         !! two copies kept in step by hand: the grounded-mask fix had to touch
         !! both sites, which is the failure mode this removes.
         load = rho_ice*d_ice_g + rho_water*(C*rsl) + wcorr
      end subroutine assemble_load

   end subroutine sle_solve

   subroutine ocean_function(topo, ice, C)
      !! Migrating-coastline ocean function with grounded-ice flotation. A cell
      !! is ocean (C = 1) only where BOTH
      !!   (a) the solid surface is below the sea surface, topo < 0, and
      !!   (b) the ice column is thin enough to float rather than ground:
      !!       ρ_i·I < −ρ_w·topo  (−topo > 0 is the water depth; the inequality
      !!       compares the ice draft to the column it would displace).
      !! Where ice grounds (ρ_i·I ≥ −ρ_w·topo) the cell is land (C = 0): the ice
      !! rests on and bears on the bed, so that column is not free ocean. With
      !! ice = 0 this reduces to the bare bathymetry test topo < 0.
      real(wp), intent(in)  :: topo(:,:)   !! solid surface vs. sea surface [m]
      real(wp), intent(in)  :: ice(:,:)    !! absolute grounded-ice thickness [m]
      real(wp), intent(out) :: C(:,:)
      where (topo < 0.0_wp .and. rho_ice*ice < -rho_water*topo)
         C = 1.0_wp
      elsewhere
         C = 0.0_wp
      end where
   end subroutine ocean_function

end module vilma_sle
