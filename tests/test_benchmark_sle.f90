program test_benchmark_sle
   !! Rung-4 validation against the Martinec et al. (2018) SLE benchmark: the full
   !! migrating-coastline sea-level equation with basin topography and a time-
   !! evolving off-pole ice cap. Standalone (NOT `make check`); case A (pure
   !! loading, no SLE) lives in test_benchmark_martinec.
   !!
   !! Usage:  test_benchmark_sle.x [CASE]      CASE in {C2, D3, E2, F1}; default E2.
   !!
   !! Configs from Martinec-2018 Table 4 (doc/refs/Martinec2018.pdf). The SBK data
   !! file letter = paper case letter + 1: SBK C2 = paper B, D3 = C, E2 = D, F1 = E
   !! (verified: cap from the fig10 uplift peak; basin from the fig12 bump; ocean
   !! mode from the far-field esl -- a migrating coastline drains a fixed-ocean case
   !! ~2x too far). Building blocks: a sqrt spherical cap (centre colat25/lon75,
   !! alpha=10deg) of height L1=1500 m or L2=500 m; a circular exponential basin
   !! B = bmax - b0 exp(-psi^2/2 sigma^2) (sigma=26deg): B1 (760/1200, colat100/
   !! lon320, shallow, far from ice) or B2 (3800/6000, colat35/lon25, 2200 m deep,
   !! near ice). Time: T0 = Heaviside (full ice held tf=10 kyr, 500 steps); T1 =
   !! linear growth over 10 kyr (500 steps) then held 5 kyr (250 steps), dt=20 yr.
   !! Ocean geometry: SLE1 = fixed (O held at the initial coastline; sle%fixed_ocean
   !! = .true.); SLE2 = time-varying (migrating coastline). lmax=128 (ntrunc=128).
   !!
   !!   CASE  paper  cap     basin       time  ocean        profiles fig10/11/12/13
   !!   C2    B      L1 1500 B1 shallow  T0    SLE1 fixed    lon75, c25, lon-40, c100
   !!   D3    C      L2 500  B2 deep     T1    SLE1 fixed    lon75, c25, lon25 , c35
   !!   E2    D      L2 500  B2 deep     T1    SLE2 migrate  lon75, c25, lon25 , c35
   !!   F1    E      L2 500  B2 deep*    T1    SLE2 migrate  lon75, c25, lon25 , c35
   !!
   !! F1 paleotopo spinup (paper E, basin B2*): B2 is the PRESENT-DAY sea-surface
   !! topography. The reference topo0 (held over each run) is iterated by a Heun-
   !! style fixed point so the deformed present-day surface topo0 - rsl_final lands
   !! on B2:  topo0 <- topo0 - 0.5*(topo0 - rsl_final - B2), until mean|.| < 1 m
   !! (giapy: topo -= 0.5*(sstopo[-1] - B2); sstopo = bedrock vs sea surface).
   !!
   !! Columns validated (all four cases): col2 vertical displacement u; col3/col4
   !! th/ph horizontal displacement (spheroidal V via resp%horizontal + the SHTns
   !! gradient eval); col5 = -N*g (geoid); col6 sea-surface = N+esl; col7 SLE =
   !! rsl = N - u + esl. Profiles (figs 10-13): circles of constant lon (col1=colat)
   !! or constant lat (col1=180+lon), sampled with sht_grid_eval_point[_horiz].
   use vilma_precision,       only: wp
   use vilma_constants,       only: pi, kyr, rho_ice, rho_water
   use vilma_earth_structure, only: earth_model, build_M3L70V01
   use vilma_radial_fe,       only: radial_fe_finalize
   use vilma_response,        only: response_destroy, response_horizontal, response, response_init_elastic, response_init_ve, response_init_null
   use vilma_sht,             only: sht_grid_destroy, sht_grid_eval_point_horiz, sht_grid_eval_point, sht_grid_analysis, sht_grid_surface_integral, sht_grid_init, sht_grid
   use vilma_sle,             only: sle_solve, sle_solver, sle_result, sle_result_grids
   use vilma_field,           only: spherical_cap, exp_basin
   implicit none

   character(*), parameter :: DIR = 'data/benchmarks/sle_martinec2018/'
   integer,  parameter :: LMAX = 128
   integer,  parameter :: NGROW = 500, NHOLD = 250, NSTEP = NGROW + NHOLD
   real(wp), parameter :: DEG = pi/180.0_wp
   real(wp), parameter :: G_GIAPY = 9.815_wp           ! giapy reference gravity
   real(wp), parameter :: ICE_COLAT = 25.0_wp*DEG, ICE_LON = 75.0_wp*DEG
   real(wp), parameter :: ALPHA_MAX = 10.0_wp*DEG, SIGB = 26.0_wp*DEG
   integer,  parameter :: NFIG = 4
   integer,  parameter :: MAXSPIN = 12                  ! F1 paleotopo iterations
   real(wp), parameter :: SPIN_TOL = 1.0_wp             ! F1: mean|present - B2| [m]

   ! Peak-normalized max-error tolerances. The Martinec-2018 SLE benchmark is a
   ! numerical INTER-CODE comparison (no analytical solution), so giapy is one
   ! implementation and a few-% scatter is expected. With the subgrid coast the
   ! migrating cases E2/F1 match to ~1% on the coupled fields; the loosest residuals
   ! are C2's basin-region uplift / horizontal (~23-28%, a small peak amplified) and
   ! the geoid/SLE (~7%). These bounds sit ~30-50% above the worst achieved value --
   ! tight enough to catch a real regression, loose enough for inter-code scatter.
   real(wp), parameter :: TOL_U = 3.0e-1_wp, TOL_N = 1.2e-1_wp, &
                          TOL_SS = 8.0e-2_wp, TOL_SLE = 1.2e-1_wp, TOL_H = 3.5e-1_wp

   ! --- per-case configuration (set by set_case) -------------------------------
   character(2) :: casename
   character(3) :: caseprefix                            ! e.g. 'E2_' for file names
   real(wp) :: H0, BAS_COLAT, BAS_LON, BMAX, B0
   real(wp) :: pfix(NFIG)                                ! fixed lon (Z) or colat (Y) [rad]
   logical  :: spinup                                    ! F1 paleotopo iteration
   logical  :: fixed_ocean_case                          ! SLE1 fixed ocean (C2/D3)
   logical  :: heaviside_case                            ! T0 Heaviside loading (C2)
   logical  :: binary_run                                ! 2nd arg "binary": force the binary coast
   integer  :: nsteps                                    ! steps in the time history
   character(8) :: fig(NFIG) = ['fig10', 'fig11', 'fig12', 'fig13']
   character(1) :: ptype(NFIG) = ['Z', 'Y', 'Z', 'Y']   ! Z: col1=colat, Y: col1=180+lon

   type(sht_grid)     :: sht
   type(earth_model)  :: em
   type(response)  :: resp
   type(sle_solver)   :: sle
   type(sle_result)   :: res
   real(wp), allocatable :: topo0(:,:), basin(:,:), ice_now(:,:)
   real(wp), allocatable :: rsl(:,:), C(:,:), tmp(:,:), dtopo(:,:), ug(:,:), Ng(:,:)
   complex(wp), allocatable :: u_lm(:), N_lm(:), rsl_lm(:), v_lm(:)
   real(wp) :: dt, esl, bary, shift, rho_ratio, spinerr
   integer  :: spin
   logical  :: ok
   ! warm-start instrumentation: per-history SLE iteration tallies (see VILMA_SLE_WARM)
   integer  :: tot_inner = 0, tot_outer = 0
   character(len=8) :: warmenv
   integer  :: envstat

   ok = .true.
   call set_case()

   dt = 0.02_wp*kyr
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   em = build_M3L70V01()
   call response_init_ve(resp, em, sht, dt)
   sle%n_inner = 20     ! n_outer left at default 3 (converged: 3 == 12 on E2)
   sle%fixed_ocean = fixed_ocean_case   ! SLE1 (C2/D3) vs SLE2 migrating (E2/F1)
   sle%subgrid     = .not. binary_run   ! subgrid coast by default; "binary" 2nd arg opts out
   ! Warm-start by DEFAULT (matches the coupling driver): the previous step's
   ! converged rsl seeds the SLE fixed point. rsl persists across steps here, so
   ! this just sets the flag. The benchmark answer is byte-identical to cold start
   ! (the fixed point is unique); warming exercises the production path and keeps it
   ! in view for real, strongly-migrating domains. Set VILMA_SLE_WARM=0 to force the
   ! legacy cold start (for the cold-vs-warm A/B comparison).
   call get_environment_variable("VILMA_SLE_WARM", warmenv, status=envstat)
   sle%warm_start = .not. (envstat == 0 .and. trim(warmenv) == "0")
   if (sle%warm_start) then
      write(*,'(a)') '   SLE start mode: WARM (reuse previous step rsl)'
   else
      write(*,'(a)') '   SLE start mode: COLD (rsl=0 each step)'
   end if

   allocate(topo0(sht%nphi,sht%nlat), basin(sht%nphi,sht%nlat), &
            ice_now(sht%nphi,sht%nlat), rsl(sht%nphi,sht%nlat), &
            C(sht%nphi,sht%nlat), tmp(sht%nphi,sht%nlat), dtopo(sht%nphi,sht%nlat))
   allocate(u_lm(sht%nlm), N_lm(sht%nlm), rsl_lm(sht%nlm), v_lm(sht%nlm))

   ! present-day basin topography; topo0 is the reference held over the run (= the
   ! basin for C2/D3/E2; iterated for F1 so the present-day surface lands on it).
   call exp_basin(sht, BAS_COLAT, BAS_LON, BMAX, B0, SIGB, basin)
   topo0 = basin

   write(*,'(3a,i0,a,i0,a)') ' ', casename, ' SLE benchmark: lmax=', LMAX, ', ', nsteps, ' steps @ dt=20 yr'

   if (spinup) then
      do spin = 1, MAXSPIN
         call run_history()                              ! fills rsl, C, ice_now, res
         dtopo = (topo0 - rsl) - basin                   ! present-day surface - B2
         spinerr = sht_grid_surface_integral(sht, abs(dtopo))/(16.0_wp*atan(1.0_wp))
         write(*,'(a,i2,a,f10.4,a)') '   spinup ', spin, ': mean|present - B2| = ', spinerr, ' m'
         if (spinerr < SPIN_TOL) exit
         topo0 = topo0 - 0.5_wp*dtopo                     ! Heun-style paleotopo update
      end do
   else
      call run_history()
   end if
   esl = res%esl

   ! converged horizontal: rebuild the converged surface load exactly as vilma_sle's
   ! commit_step does (grounded ice + ocean water) and read the spheroidal V field.
   ! resp%horizontal reuses the last begin_step's frozen drift, so v_lm is
   ! consistent with res%u_lm/res%N_lm.
   tmp = rho_ice*ice_now*(1.0_wp - C) + rho_water*(C*rsl)
   call sht_grid_analysis(sht, tmp, u_lm)                           ! u_lm reused as load_lm
   call response_horizontal(resp, sht, u_lm, v_lm)

   ! eustatic decomposition at the final state
   rho_ratio = rho_ice/rho_water
   bary  = -rho_ratio*sht_grid_surface_integral(sht, ice_now)/sht_grid_surface_integral(sht, C)
   allocate(ug(sht%nphi,sht%nlat), Ng(sht%nphi,sht%nlat))
   call sle_result_grids(res, sht, u=ug, N=Ng)
   shift = sht_grid_surface_integral(sht, C*(Ng - ug))/sht_grid_surface_integral(sht, C)
   write(*,'(a)') ''
   write(*,'(a,f10.4,a,f10.4,a,f10.4)') ' eustatic: dphi(esl)=', esl, '  barystatic=', bary, &
                                        '  ocean-mean(N-u)=', shift
   write(*,'(a,f8.4,a,f10.2)') ' ocean frac: migrated(C)=', res%ocean_frac, &
        '  C_int[sr]=', sht_grid_surface_integral(sht, C)

   ! spectral coefficients of the final converged scalar fields
   u_lm = res%u_lm;  N_lm = res%N_lm
   tmp = rsl;    call sht_grid_analysis(sht, tmp, rsl_lm)

   write(*,'(a)') ''
   write(*,'(3a)') ' Peak-normalized max error vs Martinec-2018 ', casename, ':'
   write(*,'(a)') '   fig     rows    u[%]  vth[%]  vph[%]    N[%]   ss[%]  sle[%]'
   call compare_all()

   write(*,'(a)') ''
   if (ok) then
      write(*,'(3a)') ' PASS: Martinec-2018 case ', casename, ' (SLE u + horizontal + geoid + SLE) matches'
   else
      write(*,'(3a)') ' FAIL: Martinec-2018 case ', casename, ' SLE validation did not all pass'
      call response_destroy(resp);  call sht_grid_destroy(sht);  call radial_fe_finalize();  error stop 1
   end if
   call response_destroy(resp);  call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   subroutine set_case()
      !! Parse the case argument and set the per-case configuration.
      integer :: nargs
      character(16) :: arg, arg2
      nargs = command_argument_count()
      if (nargs >= 1) then
         call get_command_argument(1, arg);  casename = trim(adjustl(arg))
      else
         casename = 'E2'
      end if
      ! The ocean load uses the subgrid sloping-coast term by default (the library
      ! default; giapy's reference is a subgrid coast, so the migrating cases E2/F1
      ! match it to ~1%). Optional 2nd arg "binary" forces the old binary coastline
      ! for comparison.
      binary_run = .false.
      if (nargs >= 2) then
         call get_command_argument(2, arg2)
         if (trim(adjustl(arg2)) == 'binary') binary_run = .true.
      end if
      ! Configurations from Martinec-2018 Table 4 (the SBK file letter = paper
      ! letter + 1): SBK C2 = paper B, D3 = C, E2 = D, F1 = E. Paper ice models:
      ! L1 (1500 m) / L2 (500 m), both at colat25/lon75; basins B1 (760/1200,
      ! colat100/lon320, shallow, far from ice) / B2 (3800/6000, colat35/lon25,
      ! deep, near ice). Time: T0 = Heaviside (full ice held 10 kyr); T1 = linear
      ! growth over 10 kyr then held 5 kyr. Ocean: SLE1 = fixed geometry, SLE2 =
      ! time-varying (migrating coastline).
      ! defaults common to D3/E2/F1 (E-series profiles); overridden per case
      pfix = [ICE_LON, 25.0_wp*DEG, 25.0_wp*DEG, 35.0_wp*DEG]
      spinup = .false.;  fixed_ocean_case = .false.;  heaviside_case = .false.
      nsteps = NSTEP                                       ! T1: 750 steps (grow+hold)
      select case (casename)
      case ('C2')                                ! paper B: L1 1500 m + B1 + T0 + SLE1 fixed
         H0 = 1500.0_wp
         BAS_COLAT = 100.0_wp*DEG;  BAS_LON = 320.0_wp*DEG
         BMAX = 760.0_wp;  B0 = 1200.0_wp
         pfix = [ICE_LON, 25.0_wp*DEG, -40.0_wp*DEG, 100.0_wp*DEG]   ! B1 profiles
         fixed_ocean_case = .true.;  heaviside_case = .true.
         nsteps = NGROW                                    ! T0: 500 steps (10 kyr), full ice
      case ('D3')                                ! paper C: L2 500 m + B2 deep + T1 + SLE1 fixed
         H0 = 500.0_wp
         BAS_COLAT = 35.0_wp*DEG;  BAS_LON = 25.0_wp*DEG
         BMAX = 3800.0_wp;  B0 = 6000.0_wp
         fixed_ocean_case = .true.
      case ('E2')                                ! paper D: L2 500 m + B2 deep + T1 + SLE2 migrating
         H0 = 500.0_wp
         BAS_COLAT = 35.0_wp*DEG;  BAS_LON = 25.0_wp*DEG
         BMAX = 3800.0_wp;  B0 = 6000.0_wp
      case ('F1')                                ! paper E: as E2 + paleotopo (B2*) spinup
         H0 = 500.0_wp
         BAS_COLAT = 35.0_wp*DEG;  BAS_LON = 25.0_wp*DEG
         BMAX = 3800.0_wp;  B0 = 6000.0_wp
         spinup = .true.
      case default
         write(*,'(3a)') ' FAIL: unknown case "', trim(casename), '" (expected C2/D3/E2/F1)'
         error stop 1
      end select
      caseprefix = casename//'_'
   end subroutine set_case

   subroutine run_history()
      !! Run the time history on the current topo0 from a relaxed, ice-free
      !! reference (zero Maxwell memory). T1: linear growth 0->full over NGROW
      !! steps (10 kyr), then held NHOLD steps (5 kyr). T0 (heaviside): full ice
      !! from the first step, held nsteps (= NGROW = 10 kyr). Leaves the converged
      !! final-step fields in rsl, C, ice_now and the diagnostics in res.
      integer :: istep
      real(wp) :: frac, alpha, h
      real(wp) :: hist(nsteps,5)
      call reset_memory()
      tot_inner = 0;  tot_outer = 0
      do istep = 1, nsteps
         if (heaviside_case) then
            frac = 1.0_wp                                   ! T0: full ice every step
         else
            frac = min(1.0_wp, real(istep,wp)/real(NGROW,wp))  ! T1 ramp 0->1, then held
         end if
         alpha = frac*ALPHA_MAX;  h = frac*H0
         call spherical_cap(sht, ICE_COLAT, ICE_LON, alpha, h, ice_now)
         ! ice-free reference: d_ice (load) = ice (flotation) = current absolute cap
         call sle_solve(sle, sht, resp, ice_now, ice_now, topo0, rsl, C, res)
         tot_inner = tot_inner + res%n_inner_last   ! inner iters in the last outer pass
         tot_outer = tot_outer + res%n_outer_done   ! coastline (outer) passes this step
         hist(istep,:) = [real(istep,wp), real(istep,wp)*dt/kyr, frac, &
                          res%esl, res%mass_resid]
         if (mod(istep,250) == 0 .or. istep == nsteps) &
            write(*,'(a,i4,a,f7.4,a,es10.2,a,f9.4,a,i3)') '   step ', istep, '  frac=', frac, &
               '  mass_resid=', res%mass_resid, '  esl=', res%esl, '  n_outer=', res%n_outer_done
      end do
      write(*,'(a,i0,a,i0,a,f5.2,a,f5.2,a)') '   SLE iteration tally over ', nsteps, &
         ' steps: outer=', tot_outer, ' (', real(tot_outer,wp)/real(nsteps,wp), &
         '/step), last-pass inner=', real(tot_inner,wp)/real(nsteps,wp), '/step'
      ! F1 calls this once per paleotopo spin-up pass; the file is replaced each
      ! time, so the one left on disk is the final, converged history.
      call dump_cols('sle_'//trim(casename)//'_history.txt', &
           'step t[kyr] ice_frac esl[m] mass_resid', hist)
   end subroutine run_history

   subroutine reset_memory()
      !! Return the response to the relaxed, ice-free reference (zero memory, t=0)
      !! so each history (and each F1 spinup pass) starts from the same state.
      call response_destroy(resp)
      call response_init_ve(resp, em, sht, dt)
   end subroutine reset_memory

   subroutine compare_all()
      integer :: k
      do k = 1, NFIG
         call compare_fig(k)
      end do
   end subroutine compare_all

   subroutine compare_fig(k)
      integer, intent(in) :: k
      integer,  parameter :: MAXR = 800
      real(wp) :: c1(MAXR), uref(MAXR), vthref(MAXR), vphref(MAXR)
      real(wp) :: nref(MAXR), ssref(MAXR), sleref(MAXR)
      real(wp) :: colat, lon, um, nm, ssm, slem, vthm, vphm
      real(wp) :: eu, evth, evph, en, ess, esle
      real(wp) :: pu, pvth, pvph, pn, pss, psle
      real(wp) :: prof(MAXR,13)
      integer  :: nrow, i

      call read_fig(trim(fig(k)), c1, uref, vthref, vphref, nref, ssref, sleref, nrow)
      eu = 0; evth = 0; evph = 0; en = 0; ess = 0; esle = 0
      pu   = maxval(abs(uref(1:nrow)));    pn   = maxval(abs(nref(1:nrow)))
      pss  = maxval(abs(ssref(1:nrow)));   psle = maxval(abs(sleref(1:nrow)))
      ! Horizontal is a VECTOR: normalize both components by the peak horizontal
      ! vector magnitude |(vth,vph)|, not each by its own peak. A profile through a
      ! symmetry axis nearly nulls one component, and per-component normalization
      ! would otherwise blow up a few-cm error on that tiny component.
      pvth = maxval(sqrt(vthref(1:nrow)**2 + vphref(1:nrow)**2))
      pvph = pvth
      do i = 1, nrow
         if (ptype(k) == 'Z') then        ! col1 = colatitude, fixed longitude
            colat = c1(i)*DEG;  lon = pfix(k)
         else                              ! col1 = 180+lon; sample lon = col1 (SHTns periodic)
            colat = pfix(k);    lon = c1(i)*DEG
         end if
         call sht_grid_eval_point(sht, u_lm,   colat, lon, um)
         call sht_grid_eval_point(sht, N_lm,   colat, lon, nm)
         call sht_grid_eval_point(sht, rsl_lm, colat, lon, slem)
         call sht_grid_eval_point_horiz(sht, v_lm, colat, lon, vthm, vphm)
         ssm = nm + esl
         prof(i,:) = [c1(i), um, uref(i), vthm, vthref(i), vphm, vphref(i), &
                      nm, nref(i), ssm, ssref(i), slem, sleref(i)]
         eu   = max(eu,   abs(um   - uref(i)))
         evth = max(evth, abs(vthm - vthref(i)))
         evph = max(evph, abs(vphm - vphref(i)))
         en   = max(en,   abs(nm   - nref(i)))
         ess  = max(ess,  abs(ssm  - ssref(i)))
         esle = max(esle, abs(slem - sleref(i)))
      end do
      eu = eu/pu;  evth = evth/pvth;  evph = evph/pvph
      en = en/pn;  ess = ess/pss;  esle = esle/psle
      write(*,'(3x,a,i7,6f8.2)') fig(k), nrow, 100*eu, 100*evth, 100*evph, &
                                 100*en, 100*ess, 100*esle
      ! The profile behind that error row: col1 is colatitude (ptype 'Z') or
      ! 180+longitude, matching the Martinec/SBK reference file convention.
      call dump_cols('sle_'//trim(casename)//'_'//trim(fig(k))//'_profile.txt', &
           'col1[deg] u_model u_ref vth_model vth_ref vph_model vph_ref '// &
           'N_model N_ref ss_model ss_ref sle_model sle_ref   [m]', prof(1:nrow,:))
      if (eu   > TOL_U)   ok = .false.
      if (evth > TOL_H)   ok = .false.
      if (evph > TOL_H)   ok = .false.
      if (en   > TOL_N)   ok = .false.
      if (ess  > TOL_SS)  ok = .false.
      if (esle > TOL_SLE) ok = .false.
   end subroutine compare_fig

   subroutine read_fig(name, c1, uref, vthref, vphref, nref, ssref, sleref, nrow)
      !! Read a 7-column *_SBK.dat profile: col1; col2=u; col3=th-, col4=ph-
      !! horizontal displacement; col5=-N*g; col6=sea-surface; col7=SLE. Reference
      !! geoid N = -col5/g_giapy.
      character(*), intent(in)  :: name
      real(wp),     intent(out) :: c1(:), uref(:), vthref(:), vphref(:)
      real(wp),     intent(out) :: nref(:), ssref(:), sleref(:)
      integer,      intent(out) :: nrow
      character(256) :: line
      real(wp) :: v(7)
      integer  :: u, ios
      open(newunit=u, file=DIR//caseprefix//name//'_SBK.dat', status='old', &
           action='read', iostat=ios)
      if (ios /= 0) then
         write(*,'(5a)') ' FAIL: cannot read ', DIR, caseprefix, name, '_SBK.dat'
         error stop 1
      end if
      nrow = 0
      do
         read(u,'(a)',iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '#') cycle
         read(line,*,iostat=ios) v
         if (ios /= 0) cycle
         nrow = nrow + 1
         c1(nrow)     = v(1)
         uref(nrow)   = v(2)
         vthref(nrow) = v(3)
         vphref(nrow) = v(4)
         nref(nrow)   = -v(5)/G_GIAPY
         ssref(nrow)  = v(6)
         sleref(nrow) = v(7)
      end do
      close(u)
   end subroutine read_fig


   subroutine dump_cols(name, header, a)
      !! Write a column table to $VILMA_BENCH_DUMP/<name> for the analysis scripts.
      !!
      !! No-op unless VILMA_BENCH_DUMP names a directory, so `make check` and any
      !! plain run behave exactly as before — the dump is opt-in and costs
      !! nothing when off. `header` names the columns and is written as a leading
      !! `#` comment line, so the file is self-describing and readable with any
      !! delimited-text reader.
      character(*), intent(in) :: name, header
      real(wp),     intent(in) :: a(:,:)          ! (nrow, ncol)
      character(512) :: dir, path
      integer :: u, i, st
      call get_environment_variable('VILMA_BENCH_DUMP', dir, status=st)
      if (st /= 0 .or. len_trim(dir) == 0) return
      path = trim(dir)//'/'//name
      open(newunit=u, file=trim(path), status='replace', action='write', iostat=st)
      if (st /= 0) then
         write(*,'(3a)') '   WARNING: cannot write ', trim(path), ' (dump skipped)'
         return
      end if
      write(u,'(2a)') '# ', header
      do i = 1, size(a,1)
         write(u,'(*(es24.16,1x))') a(i,:)
      end do
      close(u)
      write(*,'(3a,i0,a,i0,a)') '   dumped ', trim(path), ' (', size(a,1), ' x ', size(a,2), ')'
   end subroutine dump_cols

end program test_benchmark_sle
