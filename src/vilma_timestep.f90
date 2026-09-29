module vilma_timestep
   !! Adaptive time-stepping strategies for the viscoelastic + sea-level model.
   !! Isolated from vilma_coupling so the model can carry more than one stepping strategy
   !! (fixed sub-steps, adaptive, …) behind a common driver.
   !!
   !! `adaptive_stepper` advances the VE+SLE model across a coupling interval [t0,t1]
   !! with the ice load LINEARLY INTERPOLATED between its endpoints. It carries TWO
   !! sub-stepping strategies, chosen by the memory scheme:
   !!
   !!  • EXPLICIT (forward-Euler, the default): the 1st-order memory has no embedded
   !!    error signal, so Δt is set a priori by the Maxwell STABILITY ceiling
   !!    Δt ≤ cfl/max(μ/η) and the interval is divided into equal sub-steps (this is
   !!    how VILMA-v1 marches its explicit memory). A cheap reactive guard rolls a sub-step
   !!    back (resp%save_state/restore_state) and halves it if the memory ∞-norm goes
   !!    non-finite or grows past guard_growth× — a safety net for a stiffer-than-
   !!    estimated structure that normally never fires.
   !!
   !!  • IMPLICIT (trapezoidal, SCHEME_TRAP): Δt is chosen by a step-doubling local-error
   !!    estimate on the memory (the §3c adaptive-dt controller) + the SLE↔memory
   !!    co-convergence (§3c 3b): each candidate step is taken once at Δt and once as two
   !!    Δt/2 sub-steps; the coarse/fine difference is the local error (∝ Δt^{p+1}). A
   !!    rejected step rolls back and retries smaller; an accepted step keeps the fine
   !!    state and grows Δt.
   !!
   !! Δt enters the response only through Mk = (μ/η)Δt, so changing it (resp%set_dt) is a
   !! cheap rescale — no operator re-factorization (the band LU is Δt-independent).
   use vilma_precision,    only: wp
   use vilma_sht,          only: sht_grid
   use vilma_response,     only: response_coarse_fine_error, response_stash_coarse, response_prime_sigma, response_restore_state, response_set_dt, response_save_state, response_memory_norm, response_max_rate, response
   use vilma_sle,          only: sle_solve, sle_solver, sle_result
   use vilma_rotation,     only: rotation_state, rotation_save_state, rotation_restore_state
   use vilma_viscoelastic, only: scheme_order, scheme_is_implicit
   implicit none
   private

   public :: adaptive_stepper
   public :: stepper_advance

   type :: adaptive_stepper
      !! Embedded step-doubling controller. Tolerances are on the memory ∞-norm (the
      !! surface observable is algebraic in the memory τ, so this also bounds the
      !! observable's local error). Defaults are placeholders — set per problem.
      real(wp) :: rtol      = 1.0e-4_wp   !! relative local-error tolerance (memory norm)
      real(wp) :: atol      = 1.0e-3_wp   !! absolute local-error floor (memory units)
      real(wp) :: safety    = 0.9_wp      !! step-size safety factor
      real(wp) :: grow_max  = 5.0_wp      !! max Δt growth per accepted step
      real(wp) :: shrink_min = 0.2_wp     !! min Δt shrink per step (reject or grow)
      real(wp) :: dt_min    = 0.0_wp      !! Δt floor (0 = none); a step at the floor is
                                          !! accepted even if over tolerance (no infinite
                                          !! subdivision) and tallied in n_floor
      real(wp) :: dt_max    = huge(1.0_wp)!! Δt ceiling
      !! Explicit (forward-Euler) sub-stepping. The 1st-order memory carries no
      !! embedded local-error signal, so Δt is set a priori by the Maxwell stability
      !! ceiling Δt ≤ cfl/max(μ/η): the interval is split into equal sub-steps with
      !! Maxwell number M = μΔt/η ≤ cfl. cfl = 1 matches the rotation channel's
      !! convention (M ≤ 1, a factor-2 margin below the |1−M| = 1 stability bound).
      real(wp) :: cfl       = 1.0_wp      !! explicit sub-step Maxwell-number ceiling
      !! Reactive instability guard (b): a sub-step is rolled back and halved when the
      !! memory ∞-norm is non-finite OR grows by more than guard_growth× over the step
      !! (floor atol). A generous factor ⇒ a true safety net for a stiffer-than-
      !! estimated structure, not an accuracy controller; it normally never fires.
      real(wp) :: guard_growth = 1.0e3_wp !! max memory-norm growth per sub-step
      real(wp) :: dt_try    = 0.0_wp      !! next-step Δt suggestion (carried across
                                          !! intervals; 0 ⇒ first guess = whole interval)
      integer  :: n_accept  = 0           !! accepted steps (cumulative)
      integer  :: n_reject  = 0           !! rejected attempts (cumulative)
      integer  :: n_floor   = 0           !! steps accepted at dt_min over tolerance
      integer  :: n_solve   = 0           !! SLE solves issued (the dominant cost unit)
      real(wp) :: worst_mass_resid = 0.0_wp  !! worst SLE mass residual over the LAST advance()
      !! PROFILE: wall-clock [s] of the rollback bookkeeping — response_save_state /
      !! _restore_state (a full copy of the (NLAM,ne,nk) memory arrays) and the
      !! response_memory_norm scans that feed the guards. This is paid on EVERY
      !! sub-step, outside sle_solve, so without a timer it lands in the driver's
      !! residual bucket and reads as SLE cost.
      real(wp) :: t_guard = 0.0_wp
   end type adaptive_stepper

contains

   subroutine stepper_advance(self, sht, resp, sle, topo0, ice0, ice1, ice_ref, &
                              t0, t1, rsl, C, rot, sigma_out)
      !! Advance the VE+SLE model from t0 to t1 with the ice load interpolated linearly
      !! ice(t) = ice0 + (t−t0)/(t1−t0)·(ice1−ice0), absolute; the SLE load is
      !! d_ice(t) = ice(t) − ice_ref (total change from the reference state). Δt is
      !! chosen adaptively; on return resp holds the end-of-interval memory at t1 and
      !! rsl/C hold the converged sea-level fields there.
      type(adaptive_stepper),  intent(inout) :: self
      type(sht_grid),           intent(in)    :: sht
      type(response),       intent(inout) :: resp
      type(sle_solver),         intent(inout) :: sle
      real(wp),                 intent(in)    :: topo0(:,:)
      real(wp),                 intent(in)    :: ice0(:,:), ice1(:,:)  !! abs. ice at t0,t1
      real(wp),                 intent(in)    :: ice_ref(:,:)          !! reference ice
      real(wp),                 intent(in)    :: t0, t1
      real(wp),                 intent(inout) :: rsl(:,:)
      real(wp),                 intent(out)   :: C(:,:)
      !! rot (optional): rotational feedback (vilma_rotation), stepped with the solid
      !! Earth: every SLE solve iterates it to the fixed point at the end of its own
      !! sub-step, and it is saved / restored wherever the response is.
      type(rotation_state), optional, intent(inout) :: rot
      !! sigma_out (optional): the SLE's converged spectral surface mass load [kg m⁻²]
      !! at t1 — the SAME load the response saw, including the subgrid sloping-coast
      !! term. Returned so the caller drives the end-of-interval rotational feedback
      !! from it rather than re-deriving the load (vilma_coupling). The final SLE solve of
      !! the run lands on the accepted state at t1, so its load is the t1 load.
      complex(wp),    optional, intent(out)   :: sigma_out(:)

      type(sle_result)      :: res
      real(wp), allocatable :: rsl_n(:,:), ice_now(:,:), dice_now(:,:)
      complex(wp), allocatable :: sig0(:), sig_last(:)
      real(wp) :: t, dt, err_inf, tau_inf, errsc, fac, ricfac, expo, span
      real(wp) :: rate, dt_stab, tau0, tau_run, scale
      integer  :: p, np, nl, n_sub
      logical  :: at_floor, accept, finite
      integer(kind=8) :: gc0, gc1, grate          ! PROFILE: see %t_guard

      span = t1 - t0
      ! Capture the converged load on every SLE solve; after the run sig_last holds
      ! the final (t1) solve's load. Always allocated so solve_at can fill it cheaply.
      allocate(sig_last(sht%nlm))
      self%worst_mass_resid = 0.0_wp             ! per-interval diagnostic (reset each advance)

      ! Zero-length interval: no time or memory advance, but still refresh the diagnostic
      ! fields (rsl, C) so they are consistent with the CURRENT memory + ice load. A
      ! report-only SLE solve converges the load against the frozen entering memory and
      ! returns without advancing it. The drive uses solid_earth_update(...,dt=0) to seed
      ! the entering ice and populate rsl/z_bed/C from a seeded or restored memory state;
      ! a cross-resolution restart restores only the spectral memory (not the 2-D rsl), so
      ! without this the caller's rsl would stay zero on the first output sample.
      !
      ! KNOWN ISSUE -- the diagnostic is built from ice0, the ice the model is
      ! already carrying, NOT ice1, the slice the caller just handed in. The caller
      ! (solid_earth_update) overwrites gg%h_ice with the new ice the moment this
      ! returns, so after a dt=0 call rsl/C describe the PREVIOUS ice while h_ice
      ! describes the new one, and update_bsl then combines the new h_ice with the
      ! old C. Three consequences:
      !   * solid_earth_update(se, X, 0.0) is not idempotent -- calling it twice
      !     with the same X gives different rsl/C.
      !   * the two backends differ at t0. vilma_drive seeds with dt=0 and
      !     solid_earth_init leaves h_ice = h_ice_eq, so the native seed reports
      !     the REFERENCE coastline while the VILMA-v1 path (which sets h_ice before
      !     its own diagnostics) reports the LGM one. They differ over every cell
      !     where LGM grounded ice sits on sub-sea-level bed -- Hudson Bay, the
      !     Barents shelf, West Antarctica -- in both C_ocean and the t0 bsl.
      !   * it affects only the t0 diagnostic sample, not the transient physics.
      ! Switching to ice1 is the obviously right expression and was tried; it makes
      ! solid_earth_spinup stop relaxing (relax_hold's first pass then measures
      ! ~4e-7 m of bed motion instead of 136 m and exits "converged" immediately,
      ! leaving 29 m of subsidence where the equilibrium is 239 m). That interaction
      ! is not understood, and initialization I1 -- which every production
      ! deglaciation depends on -- runs through this path, so the one-line change is
      ! NOT safe to make on its own. Fixing it means understanding why a seeded rsl
      ! suppresses the subsequent viscous advance; test_spinup is the reproducer.
      if (span <= 0.0_wp) then
         allocate(ice_now(sht%nphi, sht%nlat), dice_now(sht%nphi, sht%nlat))
         ice_now  = ice0
         dice_now = ice0 - ice_ref
         call sle_solve(sle, sht, resp, dice_now, ice_now, topo0, rsl, C, res, &
                        report_only=.true., sigma_lm=sig_last, rot=rot, rot_dt=0.0_wp)
         self%worst_mass_resid = res%mass_resid
         if (present(sigma_out)) sigma_out = sig_last
         return
      end if
      if (.not. scheme_is_implicit(resp%scheme)) then
         ! --- explicit (forward-Euler) memory: a-priori stability sub-stepping -------
         ! The 1st-order memory carries no embedded error signal, so step-doubling is
         ! meaningless. Instead the interval is divided into equal sub-steps sized by
         ! the Maxwell stability ceiling Δt ≤ cfl/max(μ/η) (M = μΔt/η ≤ cfl), and a
         ! cheap reactive guard rolls back and halves a sub-step whose memory goes
         ! non-finite OR grows past guard_growth× an ESTABLISHED memory scale (a safety
         ! net for a stiffer-than-estimated structure that normally never fires). The
         ! scale is the running max of the accepted memory norm, NOT a fixed floor —
         ! otherwise the legitimate ramp-up of memory from the relaxed (τ≈0) start, where
         ! the instantaneous norm is tiny, reads as runaway growth and trips false alarms.
         ! VILMA-v1 advances its explicit memory the same way (small fixed Δt).
         np = sht%nphi;  nl = sht%nlat
         allocate(rsl_n(np,nl), ice_now(np,nl), dice_now(np,nl))
         rate = response_max_rate(resp)
         if (rate > 0.0_wp) then
            dt_stab = self%cfl/rate                  ! M = μΔt/η ≤ cfl
         else
            dt_stab = span                           ! purely elastic: one step suffices
         end if
         dt_stab = min(dt_stab, self%dt_max)
         if (self%dt_min > 0.0_wp) dt_stab = max(dt_stab, self%dt_min)
         n_sub = max(1, ceiling(span/dt_stab - 1.0e-9_wp))
         dt    = span/real(n_sub, wp)                 ! nominal (equal) sub-step
         t = t0
         call system_clock(gc0, grate)
         tau_run = response_memory_norm(resp)                ! established memory scale at entry
         call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
         do
            if (t >= t1 - 1.0e-9_wp*span) exit
            dt   = min(dt, t1 - t)
            rsl_n = rsl
            call system_clock(gc0, grate)
            tau0 = response_memory_norm(resp)                ! entering memory ∞-norm
            call response_save_state(resp)
            if (present(rot)) call rotation_save_state(rot)
            call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
            call response_set_dt(resp, dt)
            call solve_at(t + dt, dt)                 ! advances memory by dt
            call system_clock(gc0, grate)
            err_inf = response_memory_norm(resp)
            call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
            finite  = (err_inf <= huge(1.0_wp))      ! .false. for NaN / +Inf
            scale   = max(tau0, tau_run)             ! established memory magnitude
            if (scale <= 0.0_wp) then
               accept = finite                       ! no scale yet (relaxed): finite-only
            else
               accept = finite .and. (err_inf <= self%guard_growth*scale)
            end if
            if (accept) then
               tau_run = max(tau_run, err_inf)       ! grow the running scale
               t = t + dt;  self%n_accept = self%n_accept + 1
               dt = min(span/real(n_sub, wp), 2.0_wp*dt)   ! recover toward nominal
            else
               call system_clock(gc0, grate)
               call response_restore_state(resp);  rsl = rsl_n
               if (present(rot)) call rotation_restore_state(rot)
               call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
               self%n_reject = self%n_reject + 1
               dt = 0.5_wp*dt
               if (dt <= 1.0e-12_wp*span) error stop &
                  'vilma_timestep: explicit guard sub-step collapsed (viscosity too stiff for FE)'
            end if
         end do
         if (present(sigma_out)) sigma_out = sig_last
         return
      end if

      np = sht%nphi;  nl = sht%nlat
      allocate(rsl_n(np,nl), ice_now(np,nl), dice_now(np,nl))

      ! Seed the trapezoidal start-of-step load σ_0 once, if not already tracked: at the
      ! interval start the response carries its entering memory, so a report-only SLE
      ! solve (no memory/time advance) gives the load consistent with it — at t=0 that
      ! is the elastic-consistent load. Without this the first step falls back to the
      ! σ_{n+1} proxy for σ_n (O(Δt) for the relaxing SLE load → 2nd order on that step).
      if (.not. resp%sigma_primed) then
         allocate(sig0(sht%nlm))
         ice_now = ice0;  dice_now = ice0 - ice_ref
         call sle_solve(sle, sht, resp, dice_now, ice_now, topo0, rsl, C, res, &
                        report_only=.true., sigma_lm=sig0, rot=rot, rot_dt=0.0_wp)
         call response_prime_sigma(resp, sig0)
      end if

      p      = scheme_order(resp%scheme)          ! global order (2 for trapezoidal)
      ricfac = real(2**p - 1, wp)                 ! Richardson factor for the estimate
      expo   = 1.0_wp/real(p + 1, wp)             ! step-size exponent (local err ∝ Δt^{p+1})

      t = t0
      if (self%dt_try <= 0.0_wp) self%dt_try = span
      do
         if (t >= t1 - 1.0e-9_wp*span) exit
         dt = min(self%dt_try, t1 - t)
         at_floor = (self%dt_min > 0.0_wp .and. dt <= self%dt_min*(1.0_wp + 1.0e-9_wp))

         ! --- field step-doubling around [t, t+dt] -------------------------------
         rsl_n = rsl
         call system_clock(gc0, grate)
         call response_save_state(resp)                   ! buffer A = τ_n
         if (present(rot)) call rotation_save_state(rot)
         call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
         call response_set_dt(resp, dt)
         call solve_at(t + dt, dt)                 ! coarse: one Δt
         call system_clock(gc0, grate)
         call response_stash_coarse(resp)                  ! buffer B = τ_coarse
         call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
         call system_clock(gc0, grate)
         call response_restore_state(resp);  rsl = rsl_n   ! back to τ_n (and its rsl seed)
         if (present(rot)) call rotation_restore_state(rot) ! the rotation follows the fine pair
         call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
         call response_set_dt(resp, 0.5_wp*dt)
         call solve_at(t + 0.5_wp*dt, 0.5_wp*dt)   ! fine sub-step 1
         call solve_at(t + dt, 0.5_wp*dt)          ! fine sub-step 2 → τ_fine
         call system_clock(gc0, grate)
         call response_coarse_fine_error(resp, err_inf, tau_inf)
         call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
         call response_set_dt(resp, dt)                       ! leave resp%dt at the step size
         errsc = (err_inf/ricfac) / (self%atol + self%rtol*tau_inf)
         ! Same two guards the explicit path carries. Without them a non-finite
         ! estimate makes fac NaN and dt_try NaN, so t never advances; and with the
         ! default dt_min = 0 a persistently rejected step shrinks geometrically
         ! toward zero forever. Both are an infinite hang rather than an error.
         if (.not. (errsc <= huge(1.0_wp))) error stop &
            'vilma_timestep: implicit local-error estimate went non-finite &
            &(viscosity too stiff, or the state has already diverged)'

         ! --- accept / reject ----------------------------------------------------
         accept = (errsc <= 1.0_wp) .or. at_floor
         if (accept) then
            t = t + dt;  self%n_accept = self%n_accept + 1
            if (at_floor .and. errsc > 1.0_wp) self%n_floor = self%n_floor + 1
         else
            call system_clock(gc0, grate)
            call response_restore_state(resp);  rsl = rsl_n
            if (present(rot)) call rotation_restore_state(rot)
            call system_clock(gc1);  self%t_guard = self%t_guard + real(gc1-gc0,wp)/grate
            self%n_reject = self%n_reject + 1
         end if

         ! --- step-size update ---------------------------------------------------
         if (errsc <= 0.0_wp) then
            fac = self%grow_max
         else
            fac = self%safety * errsc**(-expo)
         end if
         fac = min(self%grow_max, max(self%shrink_min, fac))
         self%dt_try = min(self%dt_max, max(self%dt_min, dt*fac))
         if (self%dt_try <= 1.0e-12_wp*span) error stop &
            'vilma_timestep: implicit adaptive step collapsed (tolerance unreachable)'
      end do

      ! sig_last now holds the final accepted step's fine sub-step at t1.
      if (present(sigma_out)) sigma_out = sig_last

   contains

      subroutine solve_at(t_eval, h)
         !! One SLE solve (co-converged memory advance) with the ice load interpolated
         !! at t_eval. Advances resp by the current resp%dt (set by the caller) and the
         !! rotation by h, the same sub-step (passed: resp%dt is not set for elastic).
         real(wp), intent(in) :: t_eval, h
         real(wp) :: frac
         frac = (t_eval - t0)/span
         ice_now  = ice0 + frac*(ice1 - ice0)
         dice_now = ice_now - ice_ref
         self%n_solve = self%n_solve + 1
         call sle_solve(sle, sht, resp, dice_now, ice_now, topo0, rsl, C, res, &
                        sigma_lm=sig_last, rot=rot, rot_dt=h)
         self%worst_mass_resid = max(self%worst_mass_resid, res%mass_resid)
      end subroutine solve_at

   end subroutine stepper_advance

end module vilma_timestep
