program test_flotation_load
   !! Locks in the grounded-ice LOAD masking in the SLE (vilma_sle): floating ice
   !! does not press its full weight on the bed -- it is borne by buoyancy and
   !! carried by the ocean term -- so the surface load uses rho_ice*d_ice*(1-C),
   !! NOT rho_ice*d_ice. Without the (1-C) mask, ice overhanging deep water
   !! over-loads and spuriously subsides the bed (this was the Martinec-2018 E2
   !! west-flank error).
   !!
   !! Test: two small 3 km ice caps over a deep ocean (topo=-5000). The NORTH cap
   !! (colat<15) sits on a shallow bed (topo=-200) and GROUNDS (rho_i*3000=2.79e6
   !! >= rho_w*200=2.0e5); the SOUTH cap (colat>165) sits on the deep bed and
   !! FLOATS (2.79e6 < rho_w*5000=5.0e6). Small caps keep the eustatic drop modest
   !! so the floating cap stays floating. The grounded cap must load and strongly
   !! subside the solid; the floating cap must not subside from its own weight.
   !! We assert |u_grounded| >> |u_floating|.
   use vilma_precision,       only: wp
   use vilma_constants,       only: pi, rho_ice, rho_water
   use vilma_earth_structure, only: earth_model, build_M3L70V01
   use vilma_radial_fe,       only: radial_fe_finalize
   use vilma_response,        only: response_destroy, response, response_init_elastic, response_init_ve, response_init_null
   use vilma_sht,             only: sht_grid, sht_grid_init, sht_grid_destroy
   use vilma_sle,             only: sle_solve, sle_solver, sle_result, sle_result_grids
   implicit none

   integer, parameter :: LMAX = 24
   type(sht_grid)         :: sht
   type(earth_model)      :: em
   type(response) :: resp
   type(sle_solver)       :: sle
   type(sle_result)       :: res
   real(wp), allocatable :: topo0(:,:), ice(:,:), rsl(:,:), C(:,:), ug(:,:)
   real(wp) :: u_ground, u_float, c_ground, c_float, thd
   integer  :: i, j, jg, jf
   logical  :: ok

   ok = .true.
   call sht_grid_init(sht, LMAX, nlat=2*LMAX, nphi=4*LMAX)
   em = build_M3L70V01()
   call response_init_elastic(resp, em, LMAX)
   allocate(topo0(sht%nphi,sht%nlat), ice(sht%nphi,sht%nlat), &
            rsl(sht%nphi,sht%nlat), C(sht%nphi,sht%nlat))

   do j = 1, sht%nlat
      thd = sht%colat(j)*180.0_wp/pi
      do i = 1, sht%nphi
         if (thd < 15.0_wp) then
            topo0(i,j) = -200.0_wp;   ice(i,j) = 3000.0_wp   ! N cap: grounds (shallow)
         else if (thd > 165.0_wp) then
            topo0(i,j) = -5000.0_wp;  ice(i,j) = 3000.0_wp   ! S cap: floats  (deep)
         else
            topo0(i,j) = -5000.0_wp;  ice(i,j) =    0.0_wp   ! deep ocean
         end if
      end do
   end do

   call sle_solve(sle, sht, resp, ice, ice, topo0, rsl, C, res)

   allocate(ug(sht%nphi,sht%nlat))
   call sle_result_grids(res, sht, u=ug)
   jg = row_near(7.0_wp);  jf = row_near(173.0_wp)
   u_ground = ug(1,jg);  c_ground = C(1,jg)
   u_float  = ug(1,jf);  c_float  = C(1,jf)

   write(*,'(a)')        ' Grounded-ice LOAD masking (floating ice does not load the bed)'
   write(*,'(a,f9.3,a,f4.1,a)') '      grounded N cap (colat  7): u =', u_ground, ' m  (C=', c_ground, ')'
   write(*,'(a,f9.3,a,f4.1,a)') '      floating S cap (colat173): u =', u_float,  ' m  (C=', c_float,  ')'

   if (c_ground > 0.5_wp) then
      write(*,'(a)') '      FAIL: zone G should be grounded (C=0)'; ok = .false.
   end if
   if (c_float < 0.5_wp) then
      write(*,'(a)') '      FAIL: zone F should be floating ocean (C=1)'; ok = .false.
   end if
   if (u_ground > -50.0_wp) then
      write(*,'(a)') '      FAIL: grounded 3 km ice should strongly subside the bed'; ok = .false.
   end if
   if (abs(u_float) > 0.2_wp*abs(u_ground)) then
      write(*,'(a)') '      FAIL: floating ice subsides the bed -> load not masked by (1-C)'; ok = .false.
   end if

   write(*,'(a)') ''
   if (ok) then
      write(*,'(a)') ' PASS: floating ice does not load the bed; grounded ice does'
   else
      write(*,'(a)') ' FAIL: grounded-ice load masking incorrect'
      call response_destroy(resp);  call sht_grid_destroy(sht);  call radial_fe_finalize();  error stop 1
   end if
   call response_destroy(resp);  call sht_grid_destroy(sht);  call radial_fe_finalize()

contains

   integer function row_near(colat_deg) result(jbest)
      real(wp), intent(in) :: colat_deg
      real(wp) :: d, dbest
      integer  :: j
      jbest = 1;  dbest = huge(1.0_wp)
      do j = 1, sht%nlat
         d = abs(sht%colat(j)*180.0_wp/pi - colat_deg)
         if (d < dbest) then;  dbest = d;  jbest = j;  end if
      end do
   end function row_near

end program test_flotation_load
