program test_rotation_sle
   !! Rung-5c validation: rotational feedback coupled INTO the sea-level equation.
   !! The centrifugal potential of polar motion m perturbs the geoid and deforms the
   !! solid (Adhikari et al. 2016, eq. 8): N_rot = (1+k^T)Λ/g, u_rot = h^T Λ/g, with
   !! Λ = −Ω²a² sinθcosθ(m₁cosφ+m₂sinφ); it enters the SLE as s_rot = N_rot − u_rot.
   !! m in turn responds to the ice + ocean load, so {sea level, m} is a fixed point.
   !!
   !! No published SLE+rotation benchmark exists (Spada Test 3/2 gives only m(t)), so
   !! we validate by consistency + the analytic elastic relation:
   !!   (1) HOOK OFF: sle_solve with the rotation disabled is bit-for-bit the
   !!       no-rotation solve.
   !!   (2) FIELD vs Adhikari eq. 8: s_rot = (1+k^T_e − h^T_e)Λ/g (elastic), with Λ the
   !!       EXACT change of ½Ω²a² sin²θ' about the displaced pole ẑ + m₁x̂ + m₂ŷ, built
   !!       here from geometry, not from the formula the model uses (so it pins the sign).
   !!   (3) MASS: the rotation-coupled SLE still conserves ocean mass.
   !!   (4) FIXED POINT: sle_solve converges the rotation ↔ SLE iteration inside its
   !!       own loop, m is the polar motion of the solve's own converged load, and the
   !!       ocean feedback on m is a small correction to the ice-only polar motion.
   !!   (5) FINGERPRINT: s_rot is a degree-2 order-1, ~m-scale pattern.
   !!   (6) DIRECTION: the pole moves away from the cap (arg m = λ_c + π), and sea level
   !!       rises at the cap centre, which the pole's move takes farther from the axis.
   !!   (7) STEPPING: one open step sub-cycled n times (rotation_open_step / trial /
   !!       close_step) is bit-for-bit n steps of 1/n the length, a repeated trial is
   !!       pure, and rotation_save_state / restore_state round-trip exactly.
   !!   (8) INERTIA: the (2,1)-coefficient inertia equals a direct grid quadrature of
   !!       −a⁴∫σ sinθcosθ e^{iφ} dΩ, computed here.
   use vilma_precision,       only: wp
   use vilma_constants,       only: omega_earth
   use vilma_earth_structure, only: earth_model, build_M3L70V01
   use vilma_radial_fe,       only: radial_fe_finalize
   use vilma_response,        only: response_destroy, response, response_init_elastic, response_init_ve, response_init_null
   use vilma_sle,             only: sle_solve, sle_solver, sle_result
   use vilma_rotation,        only: rotation_destroy, rotation_s_rot, rotation_solve_m, rotation_begin_step, rotation_init, rotation_state, &
                                    rotation_open_step, rotation_trial, rotation_close_step, rotation_save_state, &
                                    rotation_restore_state, rotation_get_memory, rotation_ne, ROT_NCOMP, &
                                    rotation_inertia21
   use vilma_viscoelastic,    only: NLAM
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy, sht_grid_synthesis, &
                                    sht_grid_analysis, sht_grid_surface_integral
   implicit none

   real(wp), parameter :: deg = acos(-1.0_wp)/180.0_wp
   real(wp), parameter :: yr = 3.15576e7_wp
   real(wp), parameter :: rho_i = 931.0_wp
   real(wp), parameter :: hcap = 1500.0_wp, alpha = 10.0_wp*deg
   real(wp), parameter :: thetac = 25.0_wp*deg, lambdac = 75.0_wp*deg

   type(earth_model)      :: earth
   type(sht_grid)         :: sht
   type(response) :: resp
   type(sle_solver)       :: sle
   type(sle_result)       :: res0, res1
   type(rotation_state)   :: rot, rot_a, rot_b
   complex(wp), allocatable :: sig(:), ice_lm(:)
   real(wp), allocatable  :: w13(:,:), w23(:,:)
   complex(wp) :: Iq, Ic
   real(wp), allocatable  :: la(:,:,:), ta(:,:,:), lb(:,:,:), tb(:,:,:)
   complex(wp) :: m_a, m_b, m_chk
   real(wp)    :: h
   integer     :: k
   real(wp), allocatable  :: d_ice(:,:), ice(:,:), topo0(:,:), rsl(:,:), rsl0(:,:), C(:,:)
   real(wp), allocatable  :: srot(:,:), lam(:,:), load(:,:)
   real(wp) :: dt, expected, relfield, mcpl, mice, cp, dphase, sc
   integer  :: ilc, ipc
   complex(wp) :: m_ice
   integer  :: il, ip
   logical  :: ok

   ok = .true.
   earth = build_M3L70V01()
   call sht_grid_init(sht, 128, nlat=256, nphi=512)
   call response_init_elastic(resp, earth, 128)
   allocate(d_ice(sht%nphi,sht%nlat), ice(sht%nphi,sht%nlat), topo0(sht%nphi,sht%nlat))
   allocate(rsl(sht%nphi,sht%nlat), rsl0(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat))
   allocate(srot(sht%nphi,sht%nlat), lam(sht%nphi,sht%nlat))
   allocate(load(sht%nphi,sht%nlat))

   call build_cap(sht, d_ice)                 ! cap thickness [m]
   ice = d_ice
   where (d_ice > 0.0_wp)                      ! cap grounds on a continent; ocean elsewhere
      topo0 =  100.0_wp
   elsewhere
      topo0 = -2000.0_wp
   end where

   dt = 50.0_wp*yr
   call rotation_init(rot, earth, sht, dt)
   rot%enabled = .true.

   ! --- (1) hook off: a disabled rotation reproduces the plain SLE --------------
   rsl0 = 0.0_wp
   call sle_solve(sle, sht, resp, d_ice, ice, topo0, rsl0, C, res0)             ! no rotation
   rot%enabled = .false.;  rsl = 0.0_wp
   call sle_solve(sle, sht, resp, d_ice, ice, topo0, rsl,  C, res1, rot=rot, rot_dt=dt)
   rot%enabled = .true.
   write(*,'(a)') ' (1) hook off (rotation disabled) reproduces the plain SLE'
   write(*,'(a,es10.2,a,es10.2)') '      max|rsl diff| = ', maxval(abs(rsl - rsl0)), &
        '   |esl diff| = ', abs(res0%esl - res1%esl)
   if (maxval(abs(rsl - rsl0)) > 0.0_wp) then
      write(*,'(a)') '      FAIL: a disabled rotation changed the solution'; ok = .false.
   end if

   ! --- ice-only elastic polar motion (reference for the feedback size) --------
   call rotation_begin_step(rot, sht, dt)
   allocate(ice_lm(sht%nlm))
   load = rho_i*d_ice
   call sht_grid_analysis(sht, load, ice_lm)          ! overwrites load
   call rotation_solve_m(rot, sht, ice_lm)
   m_ice = rot%m;  mice = abs(m_ice)/deg

   ! --- (4) rotation <-> SLE fixed point inside sle_solve (elastic; at rest) ----
   ! Report-only: converge against the resting channels without committing, so rot
   ! is left holding the converged step. The m of the solve must be the polar motion
   ! of its own converged load (the closing trial).
   allocate(sig(sht%nlm))
   rsl = 0.0_wp
   call sle_solve(sle, sht, resp, d_ice, ice, topo0, rsl, C, res1, report_only=.true., &
                  sigma_lm=sig, rot=rot, rot_dt=dt)
   call rotation_s_rot(rot, sht, srot)
   mcpl = abs(rot%m)/deg
   m_chk = rot%m
   call rotation_open_step(rot, sht, 0.0_wp)
   call rotation_trial(rot, sht, sig)
   write(*,'(a)') ''
   write(*,'(a)') ' (4) rotation <-> SLE fixed point inside sle_solve'
   write(*,'(a,i0,a,es10.2,a,es10.2)') '      inner iterations ', res1%n_inner_last, &
        '   resid ', res1%resid, '   |m − m(load)|/|m| = ', abs(rot%m - m_chk)/abs(m_chk)
   if (abs(rot%m - m_chk) > 1.0e-10_wp*abs(m_chk)) then
      write(*,'(a)') '      FAIL: m is not the polar motion of the converged load'; ok = .false.
   end if
   rot%m = m_chk

   ! --- (3) mass conservation with rotation on --------------------------------
   write(*,'(a)') ''
   write(*,'(a,es10.2)') ' (3) ocean-mass residual (rotation on) = ', res1%mass_resid
   if (res1%mass_resid > 1.0e-6_wp) then
      write(*,'(a)') '      FAIL: rotation broke ocean-mass conservation'; ok = .false.
   end if

   ! --- (2) field vs Adhikari eq. 8: s_rot = (1+k^T_e − h^T_e)Λ/g, Λ exact -----
   ! Λ = ½Ω²a²(sin²θ' − sin²θ) about the displaced pole, cosθ' = (cosθ + sinθ(m₁cosφ +
   ! m₂sinφ))/√(1+|m|²). Its O(m²) remainder is ~|m| of the field, far below the tolerance.
   do il = 1, sht%nlat
      do ip = 1, sht%nphi
         cp = (cos(sht%colat(il)) + sin(sht%colat(il)) &
              * (real(rot%m,wp)*cos(sht%lon(ip)) + aimag(rot%m)*sin(sht%lon(ip)))) &
              / sqrt(1.0_wp + abs(rot%m)**2)
         lam(ip,il) = 0.5_wp*omega_earth**2*rot%a**2*(cos(sht%colat(il))**2 - cp**2)
      end do
   end do
   expected = (1.0_wp + rot%kTe - rot%hTe)/rot%g
   relfield = maxval(abs(srot - expected*lam))/maxval(abs(expected*lam))
   write(*,'(a)') ''
   write(*,'(a)') ' (2) rotational field vs Adhikari eq. 8 (elastic, exact Λ)'
   write(*,'(a,f8.4,a,f8.4)') '      k^T_e = ', rot%kTe, '   h^T_e = ', rot%hTe
   write(*,'(a,es10.2)') '      max|s_rot − (1+k−h)Λ/g| / max|(1+k−h)Λ/g| = ', relfield
   if (relfield > 1.0e-3_wp) then
      write(*,'(a)') '      FAIL: s_rot field off Adhikari eq. 8'; ok = .false.
   end if

   ! --- (4b) feedback size: ice-only vs coupled m -----------------------------
   write(*,'(a)') ''
   write(*,'(a,f10.6,a,f10.6,a,f6.2,a)') ' (4b) |m| ice-only = ', mice, ' deg ;  coupled = ', &
        mcpl, ' deg  (ocean feedback ', 100.0_wp*abs(mcpl-mice)/mice, ' %)'
   if (abs(mcpl - mice)/mice > 0.30_wp) then
      write(*,'(a)') '      FAIL: ocean feedback on m implausibly large (>30%)'; ok = .false.
   end if

   ! --- (5) fingerprint magnitude ---------------------------------------------
   write(*,'(a)') ''
   write(*,'(a,f8.3,a)') ' (5) max |s_rot| = ', maxval(abs(srot)), ' m'
   write(*,'(a,f8.4,a,f8.4)') '      k_s_fluid = ', rot%k_s_fluid, &
        '   k_s_flat (observed) = ', rot%k_s_flat
   if (maxval(abs(srot)) < 0.1_wp .or. maxval(abs(srot)) > 50.0_wp) then
      write(*,'(a)') '      FAIL: rotational fingerprint magnitude unphysical'; ok = .false.
   end if

   ! --- (6) direction: pole away from the cap, sea level up under it -----------
   dphase = modulo(atan2(aimag(m_ice), real(m_ice,wp)) - (lambdac + acos(-1.0_wp)) &
                   + acos(-1.0_wp), 2.0_wp*acos(-1.0_wp)) - acos(-1.0_wp)
   ilc = minloc(abs(sht%colat - thetac), 1);  ipc = minloc(abs(sht%lon - lambdac), 1)
   sc  = srot(ipc,ilc)
   write(*,'(a)') ''
   write(*,'(a,f8.3,a,f8.3,a)') ' (6) arg m − (λ_c + π) = ', dphase/deg, ' deg ;  s_rot at cap centre = ', sc, ' m'
   if (abs(dphase) > 1.0_wp*deg) then
      write(*,'(a)') '      FAIL: the pole does not move away from the load'; ok = .false.
   end if
   if (sc <= 0.0_wp) then
      write(*,'(a)') '      FAIL: rotational sea level does not rise where the pole moves away'; ok = .false.
   end if

   ! --- (7) stepping: sub-cycling, trial purity, save/restore ------------------
   call rotation_init(rot_a, earth, sht, dt);  rot_a%enabled = .true.
   call rotation_init(rot_b, earth, sht, dt);  rot_b%enabled = .true.
   h = 4.0_wp*rot_a%dt_fe_max
   call rotation_open_step(rot_a, sht, h)                ! one step, sub-cycled 4 times
   call rotation_trial(rot_a, sht, ice_lm)
   call rotation_trial(rot_a, sht, ice_lm)          ! a repeat trial must change nothing
   call rotation_close_step(rot_a, sht)
   do k = 1, 4                                      ! four steps of h/4
      call rotation_open_step(rot_b, sht, 0.25_wp*h)
      call rotation_trial(rot_b, sht, ice_lm)
      call rotation_close_step(rot_b, sht)
   end do
   allocate(la(NLAM, rotation_ne(rot_a), ROT_NCOMP), source=0.0_wp)
   allocate(ta, lb, tb, mold=la)
   call rotation_get_memory(rot_a, m_a, la, ta)
   call rotation_get_memory(rot_b, m_b, lb, tb)
   write(*,'(a)') ''
   write(*,'(a,i0,a,es10.2,a,es10.2,a,es10.2)') ' (7) sub-cycles ', rot_a%n_cycle, &
        ':  max|mem diff| = ', max(maxval(abs(la - lb)), maxval(abs(ta - tb))), '   |m diff| = ', abs(m_a - m_b), &
        '   |time diff| = ', abs(rot_a%time - rot_b%time)
   if (rot_a%n_cycle /= 4 .or. any(la /= lb) .or. any(ta /= tb) .or. m_a /= m_b &
       .or. rot_a%time /= rot_b%time) then
      write(*,'(a)') '      FAIL: a sub-cycled step differs from the equivalent short steps'; ok = .false.
   end if
   call rotation_save_state(rot_b)
   call rotation_open_step(rot_b, sht, h)
   call rotation_trial(rot_b, sht, 2.0_wp*ice_lm)
   call rotation_close_step(rot_b, sht)
   call rotation_restore_state(rot_b)
   call rotation_get_memory(rot_b, m_b, lb, tb)
   write(*,'(a,es10.2)') '      save/restore round trip: max|mem diff| = ', &
        max(maxval(abs(la - lb)), maxval(abs(ta - tb)))
   if (any(la /= lb) .or. any(ta /= tb) .or. m_a /= m_b .or. rot_a%time /= rot_b%time) then
      write(*,'(a)') '      FAIL: rotation_restore_state did not return the saved state'; ok = .false.
   end if
   call rotation_destroy(rot_a);  call rotation_destroy(rot_b)

   ! --- (8) inertia: (2,1) coefficient vs direct grid quadrature ---------------
   ! The converged solve load of (4) (ice + ocean + off-axis structure), in grid form.
   call sht_grid_synthesis(sht, sig, load)
   allocate(w13(sht%nphi,sht%nlat), w23(sht%nphi,sht%nlat))
   do il = 1, sht%nlat
      do ip = 1, sht%nphi
         w13(ip,il) = load(ip,il)*sin(sht%colat(il))*cos(sht%colat(il))*cos(sht%lon(ip))
         w23(ip,il) = load(ip,il)*sin(sht%colat(il))*cos(sht%colat(il))*sin(sht%lon(ip))
      end do
   end do
   Iq = cmplx(-rot%a**4*sht_grid_surface_integral(sht, w13), &
              -rot%a**4*sht_grid_surface_integral(sht, w23), wp)
   Ic = rotation_inertia21(sht, sig, rot%a)
   write(*,'(a)') ''
   write(*,'(a,es10.2)') ' (8) inertia: |I(coefficient) − I(quadrature)|/|I| = ', abs(Ic - Iq)/abs(Iq)
   if (abs(Ic - Iq) > 1.0e-10_wp*abs(Iq)) then
      write(*,'(a)') '      FAIL: the (2,1)-coefficient inertia differs from the quadrature'; ok = .false.
   end if

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: rotational feedback couples into the SLE (hook, field,'
      write(*,'(a)') '       mass, fixed point, fingerprint, direction, stepping, inertia)'
   else
      write(*,'(a)') ' FAIL: rotation-SLE coupling did not all pass'
      call sht_grid_destroy(sht);  call radial_fe_finalize();  error stop 1
   end if
   call rotation_destroy(rot);  call response_destroy(resp);  call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   subroutine build_cap(sht, thick)
      type(sht_grid), intent(in)  :: sht
      real(wp),       intent(out) :: thick(:,:)
      real(wp) :: ca, cg, frac
      integer  :: il, ip
      ca = cos(alpha);  thick = 0.0_wp
      do il = 1, sht%nlat
         do ip = 1, sht%nphi
            cg = cos(thetac)*cos(sht%colat(il)) &
               + sin(thetac)*sin(sht%colat(il))*cos(sht%lon(ip) - lambdac)
            if (cg >= ca) then
               frac = (cg - ca)/(1.0_wp - ca)
               thick(ip,il) = hcap*sqrt(max(frac, 0.0_wp))
            end if
         end do
      end do
   end subroutine build_cap

end program test_rotation_sle
