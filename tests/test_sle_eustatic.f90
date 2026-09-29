program test_sle_eustatic
   !! CONTROLLED self-check of the SLE eustatic / ocean-area handling, to settle
   !! whether our uniform offset dphi (=res%esl) is physically right or whether the
   !! ocean AREA is off. Geometry is chosen so the ocean area is KNOWN exactly:
   !!   * topo0 = +2000 m for colat<=90 (north land), -4000 m for colat>90 (south
   !!     ocean): ocean is the southern hemisphere, ocean fraction = 0.5 exactly.
   !!   * a grounded ice cap at the north pole (alpha=10 deg, h0=1000 m) -- on land,
   !!     far from the ocean, so no flotation/coastline complications.
   !! Then we decompose:
   !!   barystatic  = -rho_ratio * <ice> / <C>          (pure eustatic)
   !!   dphi(=esl)  = (ice_int - Cs_int)/C_int           (our mass-conserving offset)
   !!   def. term   = ocean_mean(N-u)  = barystatic - dphi
   !! Independent checks: ocean_frac == 0.5, mass_resid == 0. If ocean_frac is right
   !! and barystatic is sensible, the dphi-vs-barystatic gap is the (legitimate)
   !! ocean-mean deformation, NOT an ocean-area bug.
   use vilma_precision,       only: wp
   use vilma_constants,       only: pi, kyr, rho_ice, rho_water
   use vilma_earth_structure, only: earth_model, build_M3L70V01
   use vilma_radial_fe,       only: radial_fe_finalize
   use vilma_response,        only: response_destroy, response, response_init_elastic, response_init_ve, response_init_null
   use vilma_sht,             only: sht_grid_destroy, sht_grid_eval_point, sht_grid_analysis, sht_grid_surface_integral, sht_grid_init, sht_grid
   use vilma_sle,             only: sle_solve, sle_solver, sle_result, sle_result_grids
   use vilma_field,           only: spherical_cap
   implicit none

   integer, parameter :: LMAX = 48, NSTEP = 20
   real(wp), parameter :: DEG = pi/180.0_wp
   type(sht_grid)    :: sht
   type(earth_model) :: em
   type(response) :: resp
   type(sle_solver)  :: sle
   type(sle_result)  :: res
   real(wp), allocatable :: topo0(:,:), ice(:,:), rsl(:,:), C(:,:), tmp(:,:), ug(:,:), Ng(:,:)
   complex(wp), allocatable :: N_lm(:)
   real(wp) :: dt, rho_ratio, ice_int, C_int, bary, defo, Nsouth, ss_south
   integer  :: i, j, istep
   real(wp) :: th

   dt = 0.02_wp*kyr
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   em = build_M3L70V01()
   call response_init_ve(resp, em, sht, dt)
   rho_ratio = rho_ice/rho_water

   allocate(topo0(sht%nphi,sht%nlat), ice(sht%nphi,sht%nlat), &
            rsl(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat), tmp(sht%nphi,sht%nlat))
   allocate(N_lm(sht%nlm))

   do j = 1, sht%nlat
      th = sht%colat(j)
      do i = 1, sht%nphi
         topo0(i,j) = merge(2000.0_wp, -4000.0_wp, th <= 0.5_wp*pi)  ! N land / S ocean
      end do
   end do
   call spherical_cap(sht, 0.0_wp, 0.0_wp, 10.0_wp*DEG, 1000.0_wp, ice)  ! N-pole cap

   write(*,'(a)') ' SLE eustatic self-check: N-pole cap on land, S-hemisphere ocean'
   write(*,'(a)') '   step   ocean_frac   mass_resid      esl(dphi)[m]'
   do istep = 1, NSTEP
      call sle_solve(sle, sht, resp, ice, ice, topo0, rsl, C, res)
      if (mod(istep,5) == 0 .or. istep == 1) &
         write(*,'(i7,f12.5,es15.2,f16.5)') istep, res%ocean_frac, res%mass_resid, res%esl
   end do

   ! decomposition at the final state
   ice_int = -rho_ratio*sht_grid_surface_integral(sht, ice)
   C_int   = sht_grid_surface_integral(sht, C)
   bary    = ice_int/C_int                                  ! pure barystatic eustatic
   allocate(ug(sht%nphi,sht%nlat), Ng(sht%nphi,sht%nlat))
   call sle_result_grids(res, sht, u=ug, N=Ng)
   defo    = sht_grid_surface_integral(sht, C*(Ng - ug))/C_int      ! ocean-mean(N-u)
   N_lm    = res%N_lm
   call sht_grid_eval_point(sht, N_lm, pi, 0.0_wp, Nsouth)            ! geoid at S pole (deep ocean)
   ss_south = Nsouth + res%esl

   write(*,'(a)') ''
   write(*,'(a,f12.6)')  '   ocean fraction (expect 0.5000)      = ', res%ocean_frac
   write(*,'(a,f12.4)')  '   <ice> surface integral [m]          = ', sht_grid_surface_integral(sht, ice)
   write(*,'(a,f12.4)')  '   ocean area C_int [sr, expect 6.2832]= ', C_int
   write(*,'(a,f12.5)')  '   barystatic  -rho_r<ice>/<C> [m]     = ', bary
   write(*,'(a,f12.5)')  '   ocean-mean(N-u)  [m]                = ', defo
   write(*,'(a,f12.5)')  '   dphi (=res%esl) [m]                 = ', res%esl
   write(*,'(a,f12.5)')  '   bary - ocean_mean(N-u) [m]          = ', bary - defo
   write(*,'(a,es12.2)') '   |dphi - (bary - defo)| (consistency)= ', abs(res%esl - (bary - defo))
   write(*,'(a)') ''
   write(*,'(a,f12.5)')  '   geoid N at S pole [m]               = ', Nsouth
   write(*,'(a,f12.5)')  '   sea-surface (N+dphi) at S pole [m]  = ', ss_south
   write(*,'(a,f12.5)')  '   (barystatic for comparison) [m]     = ', bary

   call response_destroy(resp);  call sht_grid_destroy(sht);  call radial_fe_finalize()
end program test_sle_eustatic
