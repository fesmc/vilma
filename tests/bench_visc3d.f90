program bench_visc3d
   !! Laterally-heterogeneous GIA benchmark (setup note by Klemann et al.,
   !! doc/refs/klemann-visc3d-benchmarks-setup.pdf): a √-profile ice cap on
   !! M3-L70-V01 with a smoothed-Heaviside low-viscosity column, incompressible
   !! variants only (A-i, B-i, C-i). Writes the protocol's cross-section and
   !! geocentre files to the working directory; there is no reference solution to
   !! pass/fail against.
   !!
   !!   ./bench_visc3d.x bench_visc3d.nml      (runme: -e visc3d -n examples/bench_visc3d.nml)
   !!
   !! The run config carries the &benchv3d group (test, forcing, time_end [years],
   !! restart_in, restart_out) and &vilma overrides of input/vilma_defaults.nml. Only lmax, l_toroidal and
   !! visc3d_tol are read from &vilma: the earth model, memory scheme, degree-1
   !! frame, grid and the absence of ocean and rotation are fixed by the benchmark.
   !!
   !! Setup, and the readings taken where the note is ambiguous:
   !!   - Structure S# centred at (θ_S, λ_S), half-width ϑ_S = 12° (text; the table
   !!     says 10), edge width d = 1°, depth 60–220 km. Δη = −1e3 is read as a
   !!     factor 10⁻³ in log space, so the centre reaches η_c = 1e18 Pa s:
   !!       log10 η = log10 η_ref + (log10 η_c − log10 η_ref)·H(ϑ),
   !!       H(ϑ) = 1/(1 + exp(2(ϑ−ϑ_S)/d)).
   !!     That is −3 dex in the mantle and −12 dex in the 60–70 km lithosphere,
   !!     which the note asks to reach the same centre value.
   !!   - Earth: M3-L70-V01 with interfaces added at 60 and 220 km depth so FE
   !!     nodes bound the structure; the 60–70 km lithosphere is Maxwell with
   !!     η = 1e30 (elastic elements cannot carry a lateral viscosity).
   !!   - Load L#: σ(ϑ) = ρ_ice h_L √((cosϑ − cosα)/(1 − cosα)), α = 10°,
   !!     h_L = 1500 m, ρ_ice = 931 kg/m³; "heav" = on at t = 0, "ramp" = thickness
   !!     linear over 0–10 kyr then held. No ocean, no rotation.
   !!   - Degree 1 in the CM frame. u_gcm = u_CF − u_CM = u_CF, the surface mean
   !!     of the displacement vector, (U₁ + 2V₁)/3 in Cartesian components.
   !!   - δφ = −F = g·N: the perturbation of the gravitational potential at the
   !!     surface, including the load's own potential.
   !!   - Memory scheme: explicit (fe), Δt from the Maxwell ceiling of η_c, the
   !!     same Δt for all tests so A/B/C differ only by the structure.
   !!   - Constants are the model's own (sec_per_year, G), which differ from the
   !!     note's by ≲1e-4 relative.
   !!
   !! Cross-section: from the structure centre to its antipode through the load
   !! centre, every 0.1°. For A/B (both at the pole) this is the λ = 0 meridian.
   !! A/B are axisymmetric and run with mmax = 0.
   !!
   !! Restart, for runs longer than one job (C-i at lmax 128 is ~10-19 h): a run
   !! goes from the step in restart_in ("" = t = 0) to time_end, and with
   !! restart_out it saves the memory state at time_end. Time is kept as a global
   !! step count, so a chain of segments takes exactly the steps (and the ramp
   !! loads) of one uninterrupted run. A resumed run appends to the protocol files
   !! and skips the epochs already written.
   use vilma_precision,       only: wp
   use vilma_constants,       only: pi, rho_ice, kyr, sec_per_year
   use vilma_earth_structure, only: earth_model, earth_layer, build_M3L70V01, RHEOL_MAXWELL
   use vilma_radial_fe,       only: radial_fe_finalize
   use vilma_response,        only: response, response_init_ve, response_apply, response_horizontal, &
                                    response_horizontal_toroidal, response_begin_step, response_commit_step, &
                                    response_enable_lateral_visc, response_set_dt, response_destroy
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_lmidx, &
                                    sht_grid_eval_point, sht_grid_eval_point_horiz
   use vilma_params,          only: vilma_param_class, vilma_par_load
   use vilma_control,         only: DEFAULTS_FILE
   use nml,                   only: nml_read, nml_set_verbose
   use ncio,                  only: nc_create, nc_write_dim, nc_write, nc_read, nc_size, &
                                    nc_write_attr, nc_read_attr
   implicit none

   character(*), parameter :: CODE   = 'VILMA2'
   character(*), parameter :: AUTHOR = 'Alexander Robinson'
   real(wp), parameter :: DEG = pi/180.0_wp, KM = 1.0e3_wp, YR = 0.001_wp*kyr
   real(wp), parameter :: WIDTH_S = 12.0_wp*DEG, EDGE_D = 1.0_wp*DEG
   real(wp), parameter :: DEPTH_TOP = 60.0_wp*KM, DEPTH_BOT = 220.0_wp*KM
   real(wp), parameter :: LOG_ETA_C = 18.0_wp       !! log10 η at the structure centre
   real(wp), parameter :: ETA_LITH  = 1.0e30_wp     !! 60–70 km "elastic" Maxwell layer
   real(wp), parameter :: ALPHA_L = 10.0_wp*DEG, H_L = 1500.0_wp
   real(wp), parameter :: T_RAMP = 10.0_wp*kyr
   real(wp), parameter :: CFL = 1.0_wp              !! Maxwell number μΔt/η at η_c
   real(wp), parameter :: DS_OUT = 0.1_wp           !! cross-section spacing [deg]
   integer,  parameter :: NOUT = 13
   real(wp), parameter :: T_OUT_KYR(NOUT) = [0.1_wp, 0.2_wp, 0.5_wp, 1.0_wp, 2.0_wp, 5.0_wp, &
                                             10.0_wp, 11.0_wp, 12.0_wp, 15.0_wp, 20.0_wp, 50.0_wp, 100.0_wp]

   character(len=512) :: cfg, fdisp, fgcm, restart_in = '', restart_out = ''
   character(len=16)  :: test = 'A', forcing = 'heav', tname
   integer  :: lmax, mmax, nsub, nstep, i, i0, iout, npt, ud, ug
   real(wp) :: t_end, dt, t, colat_s, lon_s, colat_l, lon_l, dlog_s
   real(wp) :: gcm(3), mu_max
   real(wp), allocatable :: pcol(:), plon(:), pdist(:)
   complex(wp), allocatable :: load_lm(:), slm(:), ulm(:), nlm(:), vlm(:), tlm(:)
   integer(8) :: c0, c1, crate
   type(sht_grid)    :: sht
   type(earth_model) :: em
   type(response)    :: ve
   type(vilma_param_class) :: par

   ! --- run config: &benchv3d + &vilma ----------------------------------------------
   if (command_argument_count() < 1) error stop 'usage: bench_visc3d.x <run-config.nml>'
   call get_command_argument(1, cfg)
   call vilma_par_load(par, trim(cfg), DEFAULTS_FILE)
   t_end = 100.0e3_wp                                ! [years] in the namelist
   call nml_set_verbose(.false.)
   call nml_read(trim(cfg), 'benchv3d', 'test',     test)
   call nml_read(trim(cfg), 'benchv3d', 'forcing',  forcing)
   call nml_read(trim(cfg), 'benchv3d', 'time_end', t_end)
   call nml_read(trim(cfg), 'benchv3d', 'restart_in',  restart_in)
   call nml_read(trim(cfg), 'benchv3d', 'restart_out', restart_out)
   t_end = t_end*sec_per_year
   lmax  = par%lmax
   if (forcing /= 'heav' .and. forcing /= 'ramp') error stop 'benchv3d: forcing must be heav or ramp'

   ! S#/L# placement (colatitude, longitude) and the structure's log10 η drop.
   select case (trim(test))
   case ('A');  colat_s =  0.0_wp;  lon_s =  0.0_wp;  colat_l =  0.0_wp;  lon_l =  0.0_wp;  dlog_s = 0.0_wp
   case ('B');  colat_s =  0.0_wp;  lon_s =  0.0_wp;  colat_l =  0.0_wp;  lon_l =  0.0_wp;  dlog_s = 1.0_wp
   case ('C');  colat_s = 35.0_wp;  lon_s = 25.0_wp;  colat_l = 30.0_wp;  lon_l = 25.0_wp;  dlog_s = 1.0_wp
   case default; error stop 'benchv3d: test must be A, B or C'
   end select
   colat_s = colat_s*DEG;  lon_s = lon_s*DEG;  colat_l = colat_l*DEG;  lon_l = lon_l*DEG
   tname = trim(test)//'-i'

   ! --- grid, earth, response ----------------------------------------------------
   mmax = merge(0, lmax, trim(test) /= 'C')          ! A, B axisymmetric about the pole
   call sht_grid_init(sht, lmax, nlat=2*lmax+2, nphi=merge(2, 4*lmax, mmax == 0), mmax=mmax)
   em = split_earth()
   ve%deg1_cm = .true.                               ! frame is part of the init gains
   call response_init_ve(ve, em, sht, YR)
   ve%toroidal   = par%l_toroidal                    ! both read before the structure is enabled
   ve%visc3d_tol = par%visc3d_tol
   call enable_structure(ve, dlog_s, mu_max)

   ! Explicit Δt from the Maxwell ceiling at η_c, rounded to divide 100 yr (the
   ! output epochs are multiples of 0.1 kyr).
   nsub = ceiling(100.0_wp*YR / (CFL*10.0_wp**LOG_ETA_C/mu_max))
   dt   = 100.0_wp*YR/real(nsub, wp)
   call response_set_dt(ve, dt)
   nstep = nint(t_end/dt)                            ! global index of the last step
   i0 = 0
   if (len_trim(restart_in) > 0) call read_restart(trim(restart_in), i0)
   if (i0 >= nstep) error stop 'benchv3d: time_end is not after the restart time'

   allocate(load_lm(sht%nlm), slm(sht%nlm), ulm(sht%nlm), nlm(sht%nlm), vlm(sht%nlm), tlm(sht%nlm))
   call cap_load(colat_l, lon_l, load_lm)
   call section_points(colat_s, lon_s, colat_l, lon_l, pcol, plon, pdist)
   npt = size(pcol)

   write(*,'(3a,i0,a,i0,a,f6.3,a,i0,a,i0,a)') ' 3D viscosity benchmark ', trim(tname)//' '//trim(forcing), &
        ': lmax=', lmax, ' mmax=', mmax, ' dt=', dt/YR, ' yr, steps ', i0, ' -> ', nstep, ''
   write(*,'(a,i0,a,i0,a,l1)') '   radial elements: ', ve%ne, ', genuinely 3-D: ', ve%ne3d, &
        ', toroidal: ', ve%toroidal
   ! Load mass from the degree-0 coefficient vs the closed form (2π/3)ρh(1−cosα)·2a².
   write(*,'(a,es14.6,a,es14.6,a)') '   load mass: ', &
        sqrt(4.0_wp*pi)*em%r_earth**2*real(load_lm(sht_grid_lmidx(sht,0,0))), ' kg (exact ', &
        4.0_wp*pi/3.0_wp*rho_ice*H_L*(1.0_wp - cos(ALPHA_L))*em%r_earth**2, ' kg)'

   ! --- output files ---------------------------------------------------------------
   fdisp = 'disp_'//CODE//'_'//trim(tname)//'_'//trim(forcing)//'.txt'
   fgcm  = 'gcm_'//CODE//'_'//trim(tname)//'_'//trim(forcing)//'.txt'
   if (i0 == 0) then
      open(newunit=ud, file=trim(fdisp), status='replace', action='write')
      open(newunit=ug, file=trim(fgcm),  status='replace', action='write')
      call write_header(ud, 'longitude, latitude, distance on cross section, time, '// &
                        'u_r, u_\vartheta, u_\varphi, \delta_\phi')
      call write_header(ug, 'time, u_x, u_y, u_z')
      write(ud,'(a)') '# units: deg, deg, deg, kyr, m, m, m, m^2 s^-2'
      write(ug,'(a)') '# units: kyr, m, m, m'
   else                                              ! resumed: the earlier segments wrote the rest
      open(newunit=ud, file=trim(fdisp), status='old', position='append', action='write')
      open(newunit=ug, file=trim(fgcm),  status='old', position='append', action='write')
   end if

   ! --- time loop: output at t_i uses σ(t_i) and memory τ(t_i); commit → τ(t_{i+1}) --
   call system_clock(c0, crate)
   iout = 1
   if (i0 > 0) then                                  ! a previous segment wrote every epoch up to
      do while (iout <= NOUT)                        ! and including its last step, i0
         if (T_OUT_KYR(iout)*kyr > (real(i0, wp) + 0.5_wp)*dt) exit
         iout = iout + 1
      end do
   end if
   do i = i0, nstep
      t = real(i, wp)*dt
      slm = load_factor(t)*load_lm
      call response_begin_step(ve, sht)
      call response_apply(ve, sht, slm, ulm, nlm)
      if (iout <= NOUT) then
         if (abs(t - T_OUT_KYR(iout)*kyr) < 0.5_wp*dt) then
            call response_horizontal(ve, sht, slm, vlm)
            call response_horizontal_toroidal(ve, sht, tlm)
            call write_section(ud, t)
            gcm = geocentre(ulm, vlm)
            write(ug,'(f8.2,3(a,es16.8))') t/kyr, char(9), gcm(1), char(9), gcm(2), char(9), gcm(3)
            call system_clock(c1)
            write(*,'(a,f7.2,a,es12.4,a,f9.1,a)') '   t=', t/kyr, ' kyr  u_r(centre)=', &
                 point_ur(1), ' m   wall=', real(c1-c0,wp)/real(crate,wp), ' s'
            iout = iout + 1
         end if
      end if
      if (i < nstep) call response_commit_step(ve, sht, slm)
   end do
   close(ud);  close(ug)
   write(*,'(4a)') '   wrote ', trim(fdisp), ', ', trim(fgcm)
   ! The memory is τ(t_end): the last step was evaluated but not committed, so a
   ! resumed run begins with exactly the step an uninterrupted one takes next.
   if (len_trim(restart_out) > 0) then
      call write_restart(trim(restart_out), nstep)
      write(*,'(3a,f8.2,a)') '   wrote restart ', trim(restart_out), ' at t=', real(nstep,wp)*dt/kyr, ' kyr'
   end if

   call response_destroy(ve);  call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   function split_earth() result(e)
      !! M3-L70-V01 with interfaces at the structure's top and bottom (60, 220 km
      !! depth), the 60–70 km lithosphere made Maxwell (η = ETA_LITH) so it can
      !! carry the lateral viscosity.
      type(earth_model) :: e
      e = build_M3L70V01()
      call split_layer(e, e%r_earth - DEPTH_TOP)
      call split_layer(e, e%r_earth - DEPTH_BOT)
      e%layers(2)%rheology = RHEOL_MAXWELL
      e%layers(2)%eta      = ETA_LITH
      e%name = 'M3-L70-V01-split60-220'
   end function split_earth

   subroutine split_layer(e, r)
      !! Insert an interface at radius r inside the layer containing it; both
      !! halves keep that layer's properties.
      type(earth_model), intent(inout) :: e
      real(wp),          intent(in)    :: r
      type(earth_layer), allocatable :: lay(:)
      integer :: k, n
      n = size(e%layers)
      do k = 1, n
         if (r > e%layers(k)%r_bot .and. r < e%layers(k)%r_top) exit
      end do
      if (k > n) error stop 'split_layer: radius not strictly inside a layer'
      allocate(lay(n+1))
      lay(1:k)   = e%layers(1:k)
      lay(k+1:)  = e%layers(k:)
      lay(k)%r_bot   = r                              ! upper half (surface-first)
      lay(k+1)%r_top = r                              ! lower half
      call move_alloc(lay, e%layers)
   end subroutine split_layer

   subroutine enable_structure(resp, frac, mu_max)
      !! Lateral log10 η perturbation per element: (LOG_ETA_C − log10 η_ref(e))·H(ϑ)
      !! inside 60–220 km depth, scaled by frac (0 for S0, 1 for S1/S2). Returns the
      !! largest μ over the structure's elements (for the Maxwell Δt ceiling).
      type(response), intent(inout) :: resp
      real(wp),       intent(in)    :: frac
      real(wp),       intent(out)   :: mu_max
      real(wp), allocatable :: pert(:,:,:), hs(:,:)
      real(wp) :: ps(3), p(3), cosd, logeta_ref
      integer  :: ie, j, k
      allocate(pert(sht%nphi, sht%nlat, resp%ne), hs(sht%nphi, sht%nlat))
      ps = unit_vec(colat_s, lon_s)
      do j = 1, sht%nlat
         do k = 1, sht%nphi
            p = unit_vec(sht%colat(j), sht%lon(k))
            cosd = max(-1.0_wp, min(1.0_wp, dot_product(ps, p)))
            hs(k,j) = 1.0_wp/(1.0_wp + exp(2.0_wp*(acos(cosd) - WIDTH_S)/EDGE_D))
         end do
      end do
      pert = 0.0_wp;  mu_max = 0.0_wp
      do ie = 1, resp%ne
         if (resp%r(ie)   < em%r_earth - DEPTH_BOT - 1.0_wp) cycle
         if (resp%r(ie+1) > em%r_earth - DEPTH_TOP + 1.0_wp) cycle
         if (resp%MkPerDt(ie) == 0.0_wp) error stop 'enable_structure: structure element is not Maxwell'
         logeta_ref = log10(resp%mu(ie)/resp%MkPerDt(ie))
         pert(:,:,ie) = frac*(LOG_ETA_C - logeta_ref)*hs
         mu_max = max(mu_max, resp%mu(ie))
      end do
      call response_enable_lateral_visc(resp, sht, pert)
   end subroutine enable_structure

   pure function unit_vec(colat, lon) result(v)
      real(wp), intent(in) :: colat, lon
      real(wp) :: v(3)
      v = [sin(colat)*cos(lon), sin(colat)*sin(lon), cos(colat)]
   end function unit_vec

   pure real(wp) function load_factor(t) result(f)
      !! Ice-thickness fraction of h_L at time t.
      real(wp), intent(in) :: t
      if (forcing == 'heav') then
         f = 1.0_wp
      else
         f = min(t/T_RAMP, 1.0_wp)
      end if
   end function load_factor

   subroutine cap_load(colat0, lon0, lm)
      !! SH coefficients [kg/m²] of the √-profile cap centred at (colat0, lon0).
      !! Zonal coefficients about its own axis by Gauss–Legendre quadrature with
      !! x = c + (1−c)s², which turns the √ edge into a polynomial in s — exact for
      !! nq ≥ lmax+2 — then rotated by the addition theorem,
      !!   f_lm = √(4π/(2l+1)) a_l · conj(Y_lm(colat0, lon0)).
      !! Y_lm at the centre comes from point evaluations of unit coefficients, so it
      !! carries SHTns' normalization (a real field stores m > 0 once: eval of 1 is
      !! 2 Re Y_lm, eval of i is −2 Im Y_lm).
      real(wp),    intent(in)  :: colat0, lon0
      complex(wp), intent(out) :: lm(:)
      real(wp), allocatable :: s(:), w(:), p(:,:)
      real(wp) :: c, a_l, x, e1, e2
      complex(wp) :: ylm
      complex(wp), allocatable :: unit(:)
      integer :: nq, q, l, m, idx
      c  = cos(ALPHA_L);  nq = lmax + 4
      call gauss_legendre01(nq, s, w)
      allocate(p(0:lmax, nq), unit(sht%nlm))
      do q = 1, nq
         x = c + (1.0_wp - c)*s(q)**2
         p(0,q) = 1.0_wp;  if (lmax >= 1) p(1,q) = x
         do l = 1, lmax - 1
            p(l+1,q) = (real(2*l+1,wp)*x*p(l,q) - real(l,wp)*p(l-1,q))/real(l+1,wp)
         end do
      end do
      lm = (0.0_wp, 0.0_wp);  unit = (0.0_wp, 0.0_wp)
      do l = 0, lmax
         ! a_l = 2π √((2l+1)/4π) ρ h ∫_c^1 √((x−c)/(1−c)) P_l dx,  dx = 2(1−c) s ds
         a_l = 0.0_wp
         do q = 1, nq
            a_l = a_l + w(q)*s(q)*p(l,q)*2.0_wp*(1.0_wp - c)*s(q)
         end do
         a_l = 2.0_wp*pi*sqrt(real(2*l+1,wp)/(4.0_wp*pi))*rho_ice*H_L*a_l
         do m = 0, min(l, sht%mmax)
            idx = sht_grid_lmidx(sht, l, m)
            unit(idx) = (1.0_wp, 0.0_wp);  call sht_grid_eval_point(sht, unit, colat0, lon0, e1)
            if (m == 0) then
               ylm = cmplx(e1, 0.0_wp, wp)
            else
               unit(idx) = (0.0_wp, 1.0_wp);  call sht_grid_eval_point(sht, unit, colat0, lon0, e2)
               ylm = cmplx(0.5_wp*e1, -0.5_wp*e2, wp)
            end if
            unit(idx) = (0.0_wp, 0.0_wp)
            lm(idx) = sqrt(4.0_wp*pi/real(2*l+1,wp))*a_l*conjg(ylm)
         end do
      end do
   end subroutine cap_load

   subroutine gauss_legendre01(n, x, w)
      !! n-point Gauss–Legendre nodes and weights on [0, 1] (Newton on P_n).
      integer,  intent(in) :: n
      real(wp), allocatable, intent(out) :: x(:), w(:)
      real(wp) :: z, z1, p1, p2, p3, pp
      integer  :: i, j, it
      allocate(x(n), w(n))
      do i = 1, (n+1)/2
         z = cos(pi*(real(i,wp) - 0.25_wp)/(real(n,wp) + 0.5_wp))
         do it = 1, 100
            p1 = 1.0_wp;  p2 = 0.0_wp
            do j = 1, n
               p3 = p2;  p2 = p1
               p1 = (real(2*j-1,wp)*z*p2 - real(j-1,wp)*p3)/real(j,wp)
            end do
            pp = real(n,wp)*(z*p1 - p2)/(z*z - 1.0_wp)
            z1 = z;  z = z1 - p1/pp
            if (abs(z - z1) < 1.0e-15_wp) exit
         end do
         x(i) = 0.5_wp*(1.0_wp - z);  x(n+1-i) = 0.5_wp*(1.0_wp + z)
         w(i) = 1.0_wp/((1.0_wp - z*z)*pp*pp);  w(n+1-i) = w(i)
      end do
   end subroutine gauss_legendre01

   subroutine section_points(cs, ls, cl, ll, pc, pl, pd)
      !! Great circle from the structure centre through the load centre to the
      !! antipode, every DS_OUT degrees: colatitude, longitude [rad], distance [deg].
      !! When the two centres coincide (A, B at the pole) the path heads along λ = 0.
      real(wp), intent(in) :: cs, ls, cl, ll
      real(wp), allocatable, intent(out) :: pc(:), pl(:), pd(:)
      real(wp) :: a(3), b(3), tng(3), p(3), s
      integer  :: n, k
      a = unit_vec(cs, ls);  b = unit_vec(cl, ll)
      tng = b - dot_product(a, b)*a
      if (norm2(tng) < 1.0e-12_wp) tng = unit_vec(cs + 0.5_wp*pi, ls) - &
                                         dot_product(a, unit_vec(cs + 0.5_wp*pi, ls))*a
      tng = tng/norm2(tng)
      n = nint(180.0_wp/DS_OUT) + 1
      allocate(pc(n), pl(n), pd(n))
      do k = 1, n
         pd(k) = real(k-1, wp)*DS_OUT
         s = pd(k)*DEG
         p = cos(s)*a + sin(s)*tng
         pc(k) = acos(max(-1.0_wp, min(1.0_wp, p(3))))
         pl(k) = modulo(atan2(p(2), p(1)), 2.0_wp*pi)
      end do
      ! At a pole the longitude (and so the (θ, φ) basis of u_θ, u_φ) is undefined;
      ! take the meridian the path arrives along (leaves along, at the start).
      do k = 1, n
         if (sin(pc(k)) < 1.0e-9_wp) pl(k) = pl(merge(k+1, k-1, k == 1))
      end do
   end subroutine section_points

   subroutine write_header(u, cols)
      integer,      intent(in) :: u
      character(*), intent(in) :: cols
      character(len=8) :: d
      call date_and_time(date=d)
      write(u,'(a)') '# '//CODE//' (VILMA2, incompressible)'
      write(u,'(a)') '# '//trim(tname)//' ('//trim(forcing)//'), lmax='//itoa(lmax)
      write(u,'(a)') '# '//AUTHOR
      write(u,'(a)') '# '//d(1:4)//'-'//d(5:6)//'-'//d(7:8)
      write(u,'(a)') '# '//cols
   end subroutine write_header

   function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=16) :: b
      write(b,'(i0)') i;  s = trim(b)
   end function itoa

   subroutine write_section(u, t)
      !! One epoch of the cross-section: u_r, u_θ, u_φ (spheroidal + toroidal), δφ.
      integer,  intent(in) :: u
      real(wp), intent(in) :: t
      real(wp) :: ur, nn, uth, uph
      integer  :: k
      do k = 1, npt
         call sht_grid_eval_point(sht, ulm, pcol(k), plon(k), ur)
         call sht_grid_eval_point(sht, nlm, pcol(k), plon(k), nn)
         call sht_grid_eval_point_horiz(sht, vlm, pcol(k), plon(k), uth, uph, t_lm=tlm)
         write(u,'(f9.4,a,f9.4,a,f6.1,a,f8.2,4(a,es16.8))') plon(k)/DEG, char(9), &
              90.0_wp - pcol(k)/DEG, char(9), pdist(k), char(9), t/kyr, char(9), ur, &
              char(9), uth, char(9), uph, char(9), ve%g*nn
      end do
   end subroutine write_section

   real(wp) function point_ur(k) result(ur)
      integer, intent(in) :: k
      call sht_grid_eval_point(sht, ulm, pcol(k), plon(k), ur)
   end function point_ur

   function geocentre(u_lm, v_lm) result(g)
      !! u_CF in the CM frame: the surface mean of u = U r̂ + ∇₁V. For degree 1,
      !! U = a·r̂ and V = b·r̂ give ⟨u⟩ = (a + 2b)/3, with the Cartesian vector of a
      !! real degree-1 field f = √(3/4π)(√2 Re f₁₁, −√2 Im f₁₁, f₁₀) (orthonormal,
      !! no Condon–Shortley phase; m > 0 stored once).
      complex(wp), intent(in) :: u_lm(:), v_lm(:)
      real(wp) :: g(3)
      g = (cart1(u_lm) + 2.0_wp*cart1(v_lm))/3.0_wp
   end function geocentre

   function cart1(f) result(v)
      complex(wp), intent(in) :: f(:)
      real(wp) :: v(3)
      complex(wp) :: f11
      f11 = (0.0_wp, 0.0_wp)
      if (sht%mmax >= 1) f11 = f(sht_grid_lmidx(sht, 1, 1))
      v = sqrt(3.0_wp/(4.0_wp*pi))*[sqrt(2.0_wp)*real(f11), -sqrt(2.0_wp)*aimag(f11), &
                                    real(f(sht_grid_lmidx(sht, 1, 0)))]
   end function cart1

   subroutine write_restart(fname, istep)
      !! Save the Maxwell memory (the only prognostic state of the explicit
      !! scheme: begin_step re-derives the drift from it) and the step index, with
      !! the case identity and Δt so a mismatched resume is refused.
      character(*), intent(in) :: fname
      integer,      intent(in) :: istep
      call nc_create(fname, overwrite=.true.)
      call nc_write_dim(fname, 'nlam', x=1, dx=1, nx=ve%nlam, units='1')
      call nc_write_dim(fname, 'ne',   x=1, dx=1, nx=ve%ne,   units='1')
      call nc_write_dim(fname, 'nk',   x=1, dx=1, nx=ve%nk,   units='1')
      call nc_write(fname, 'tau_a_re', ve%Are, dim1='nlam', dim2='ne', dim3='nk')
      call nc_write(fname, 'tau_a_im', ve%Aim, dim1='nlam', dim2='ne', dim3='nk')
      call nc_write(fname, 'tau_b_re', ve%Bre, dim1='nlam', dim2='ne', dim3='nk')
      call nc_write(fname, 'tau_b_im', ve%Bim, dim1='nlam', dim2='ne', dim3='nk')
      call nc_write(fname, 'tau_c_re', ve%Cre, dim1='nlam', dim2='ne', dim3='nk')
      call nc_write(fname, 'tau_c_im', ve%Cim, dim1='nlam', dim2='ne', dim3='nk')
      call nc_write_attr(fname, 'step', istep)
      call nc_write_attr(fname, 'dt_s', dt)
      call nc_write_attr(fname, 'lmax', lmax)
      call nc_write_attr(fname, 'case', trim(tname)//'_'//trim(forcing))
   end subroutine write_restart

   subroutine read_restart(fname, istep)
      !! Restore the memory saved by write_restart into the freshly built response
      !! (same test, forcing, lmax, Δt and channel count, or stop) and return the
      !! global step index to resume from.
      character(*), intent(in)  :: fname
      integer,      intent(out) :: istep
      character(len=64) :: cs
      integer  :: lm
      real(wp) :: dts
      call nc_read_attr(fname, 'case', cs)
      call nc_read_attr(fname, 'lmax', lm)
      call nc_read_attr(fname, 'dt_s', dts)
      call nc_read_attr(fname, 'step', istep)
      if (trim(cs) /= trim(tname)//'_'//trim(forcing)) error stop 'benchv3d: restart is from another case'
      if (lm /= lmax)  error stop 'benchv3d: restart lmax differs'
      if (dts /= dt)   error stop 'benchv3d: restart dt differs'
      if (nc_size(fname, 'nlam') /= ve%nlam .or. nc_size(fname, 'ne') /= ve%ne .or. &
          nc_size(fname, 'nk') /= ve%nk) error stop 'benchv3d: restart memory shape differs'
      call nc_read(fname, 'tau_a_re', ve%Are);  call nc_read(fname, 'tau_a_im', ve%Aim)
      call nc_read(fname, 'tau_b_re', ve%Bre);  call nc_read(fname, 'tau_b_im', ve%Bim)
      call nc_read(fname, 'tau_c_re', ve%Cre);  call nc_read(fname, 'tau_c_im', ve%Cim)
      write(*,'(3a,i0,a,f8.2,a)') '   restart ', fname, ': step ', istep, ' (t=', &
           real(istep,wp)*dt/kyr, ' kyr)'
   end subroutine read_restart

end program bench_visc3d
