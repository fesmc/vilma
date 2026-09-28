module vilma_v1
   !! OPTIONAL VILMA-v1 backend for the vilma_coupling API (&vilma solver = "v1").
   !!
   !! VILMA-v1 (Martinec/Klemann; the CLIMBER-X i_geo=2 solid-earth backend) is swapped
   !! in behind the SAME driver, namelist, forcing, remap and output as the native
   !! VILMA solver, so an v2-vs-v1 comparison differs only in the solver. The
   !! contract is the one vilma_coupling already states: ice thickness in, relative sea
   !! level out. This module is the direct analogue of CLIMBER-X's src/geo/vilma.F90,
   !! ported to drive from the VILMA Gauss grid instead of the CLIMBER-X geo
   !! grid; read the two side by side.
   !!
   !! ===========================================================================
   !! VILMA-v1 IS NOT A DEPENDENCY OF VILMA.
   !! ===========================================================================
   !! It is a hand-installed, precompiled library (its own .mod files under
   !! include, plus lib/vega_pism.a) that is absent on most machines. Do not spell
   !! that include glob out here: the sources are compiled through the C
   !! preprocessor (-cpp for gfortran, -fpp for Intel), which reads the slash-star
   !! as the start of a C comment and swallows the rest of the file.
   !! Every reference to VILMA-v1 in this file sits
   !! inside `#ifdef VILMA_V1`, which only `make vilma vilma_v1=1 VILMA_V1_ROOT=<install>`
   !! defines. In the DEFAULT build this file compiles to a pure-Fortran stub that
   !! names no VILMA-v1 symbol and needs no VILMA-v1 include path; selecting
   !! solver="v1" then aborts at init, from vilma_v1_require, with an actionable
   !! message — never a link error and never a crash. See doc/vilma-v1-backend.md.
   !!
   !! --- Grids -----------------------------------------------------------------
   !! There are THREE grids in play, and the boundaries matter:
   !!   host lon-lat     — the driver's forcing grid. vilma_coupling/vilma_drive remap it
   !!                      to the model Gauss grid; this module never sees it.
   !!   model Gauss      — the VILMA Gauss-Legendre grid (sht, nphi x nlat, rows
   !!                      NORTH-first). Everything this module is handed and returns
   !!                      is on THIS grid, exactly as for the native solver, so all
   !!                      output files and diagnostics are directly comparable.
   !!   VILMA-v1 grid       — VILMA-v1's own Gauss-Legendre grid at vg%jmax, whose lon/lat
   !!                      axes are read from par%vilma_v1_grid_file (512 x 256 for the
   !!                      CLIMBER-X jmax=170 setup). VILMA-v1's `rsl` lives here, as
   !!                      (lat,lon) — note the transpose.
   !! This module owns the second remap leg and nothing else:
   !!   model Gauss --(conservative, coords "con")--> VILMA-v1   for ice thickness
   !!   VILMA-v1 --(bilinear, coords "bilinear")--> model Gauss   for relative sea level
   !! mirroring the conservative-in / bilinear-out pairing vilma_remap uses for the
   !! host leg and CLIMBER-X uses for its own VILMA-v1 coupling.
   !!
   !! --- Clock ------------------------------------------------------------------
   !! VILMA-v1's time axis is in kyr. CLIMBER-X feeds it calendar ka BP; this wrapper
   !! feeds it the MODEL time (vilma_coupling's se%time) divided by 1000, which for the
   !! standalone driver is the forcing file's own time axis (e.g. -26.0 -> 0.0 kyr
   !! for a deglaciation). VILMA-v1's viscoelastic response depends only on time
   !! DIFFERENCES, so the labelling affects only VILMA-v1's own diagnostic output files.
   !!
   !! --- Deferred setup ---------------------------------------------------------
   !! VILMA-v1's `setup` needs both the integration time step and the start time, and
   !! the vilma_coupling API supplies neither at init (the driver's coupling interval
   !! comes from the forcing axis, and the clock is set after init). So init does
   !! everything that does not need them — grids, maps, the reference NetCDFs — and
   !! `setup` is called on the FIRST advancing update, from that update's own dt and
   !! start time. A zero-length update (the driver's seed call) only records the
   !! entering ice, exactly as it does for the native solver.
   use vilma_precision, only: wp
   use vilma_constants, only: rad2deg, rho_ice, rho_water, sec_per_year
   use vilma_params,    only: vilma_param_class
   use vilma_sht,       only: sht_grid, sht_grid_surface_integral
   use coords,       only: grid_class, grid_init, map_class, map_init, map_field
   use ncio
   use iso_fortran_env, only: error_unit
#ifdef VILMA_V1
   ! --- the VILMA-v1 library (only under `make vilma_v1=1`) --------------------------
   ! explicit `only:` lists throughout, so no VILMA-v1 name can shadow one of ours
   use mod_struct_vg,  only: vg
   use mod_sle,        only: rsl
   use mod_io,         only: io_densi, io_visko, io_nc3in, io_tint, io_sliin, &
                             io_tmp, io_lis, io_wepoch, io_hist, io_surf, io_radii, &
                             io_nc3out, io_sliout, io_hsliout, io_rsl, io_dfgl, &
                             io_rsl_rs, io_dfgl_rs, io_deg1, io_oce, io_rpt, io_rslog, &
                             io_mos_indx, io_mos_amtrx, io_mos_acomp, io_mos_acompl, &
                             io_ve_struct, io_pjj, io_pefgh, io_nwl_struct, io_visc3drs, &
                             io_disp, io_stress, io_1dstress, io_rsrpt
   use mod_firstlevel, only: setup, time_evolution, close_evolution
#endif
   implicit none
   private

   public :: vilma_v1_backend
   public :: vilma_v1_init, vilma_v1_update, vilma_v1_finalize
   public :: vilma_v1_available, vilma_v1_require

   integer, parameter :: BIL_NEIGHBORS = 8     !! neighbour pool for the bilinear leg
   integer, parameter :: VILMA_V1_PATH_LEN = 120  !! observed limit: paths longer than
                                               !! this are silently truncated
   real(wp), parameter :: FOURPI = 12.566370614359172_wp

   type :: vilma_v1_backend
      !! Everything this wrapper owns between calls. The VILMA-v1 library itself keeps
      !! its state in its own module variables (it is a singleton), so only ONE
      !! vilma_v1_backend may be active in a process — the same restriction CLIMBER-X
      !! lives with.
      logical :: active  = .false.   !! init has run
      logical :: started = .false.   !! VILMA-v1 `setup` has been called
      integer :: nlon = 0, nlat = 0  !! VILMA-v1 grid dimensions (lon, lat)
      integer :: nphi = 0, ngauss = 0!! model Gauss grid dimensions
      real(wp) :: t_kyr = 0.0_wp     !! VILMA-v1 clock: end of the last completed interval [kyr]
      integer  :: nsub = 1           !! VILMA-v1 sub-steps per coupling interval (p%vilma_v1_nsub)

      character(len=512) :: out_dir = ""       !! VILMA-v1 scratch/output directory
      real(wp), allocatable :: lon(:), lat(:)  !! VILMA-v1 grid axes [degrees]

      type(grid_class) :: gauss                !! model Gauss grid, lat ASCENDING (coords order)
      type(grid_class) :: vgrid                !! VILMA-v1 grid
      type(map_class)  :: to_vilma_v1             !! model Gauss -> VILMA-v1, conservative
      type(map_class)  :: to_gauss             !! VILMA-v1 -> model Gauss, bilinear

      real(wp), allocatable :: h_ice_v(:,:)    !! current ice on the VILMA-v1 grid [m]
      real(wp), allocatable :: h_ref_v(:,:)    !! reference ice on the VILMA-v1 grid [m] (history slice 1)
      real(wp), allocatable :: work_v(:,:)     !! (nlon,nlat) scratch
   end type vilma_v1_backend

contains

   ! --- availability -----------------------------------------------------------

   pure logical function vilma_v1_available() result(ok)
      !! .true. only in a binary built with `make vilma vilma_v1=1`.
#ifdef VILMA_V1
      ok = .true.
#else
      ok = .false.
#endif
   end function vilma_v1_available

   subroutine vilma_v1_require()
      !! Abort with an actionable message when solver="v1" is selected in a
      !! binary that was not built with the backend. This is THE guard that keeps
      !! the default build honest: no link error, no silent no-op, no crash.
      if (vilma_v1_available()) return
      write(error_unit,'(a)') ''
      write(error_unit,'(a)') ' ======================================================================'
      write(error_unit,'(a)') '  VILMA: solver="v1" requested, but this binary has no VILMA-v1.'
      write(error_unit,'(a)') ' ======================================================================'
      write(error_unit,'(a)') '  The VILMA-v1 backend is OPTIONAL and is OFF by default, because VILMA-v1 is'
      write(error_unit,'(a)') '  a hand-installed precompiled library that is absent on most machines.'
      write(error_unit,'(a)') ''
      write(error_unit,'(a)') '  To use it, rebuild with the backend switched on:'
      write(error_unit,'(a)') ''
      write(error_unit,'(a)') '      make clean'
      write(error_unit,'(a)') '      make vilma vilma_v1=1 VILMA_V1_ROOT=/path/to/vilma1'
      write(error_unit,'(a)') ''
      write(error_unit,'(a)') '  where VILMA_V1_ROOT holds  include/*.mod  and  lib/vega_pism.a .'
      write(error_unit,'(a)') ''
      write(error_unit,'(a)') '  Otherwise set  solver = "v2"  in the &vilma namelist group to use'
      write(error_unit,'(a)') '  the native VILMA solver.  See doc/vilma-v1-backend.md.'
      write(error_unit,'(a)') ' ======================================================================'
      flush(error_unit)
      error stop 'solver="v1" requires building with `make vilma vilma_v1=1` and a VILMA-v1 install at VILMA_V1_ROOT'
   end subroutine vilma_v1_require

   ! --- lifecycle --------------------------------------------------------------

   subroutine vilma_v1_init(self, par, sht, z_bed_eq, h_ice_eq, h_ice)
      !! Prepare the backend. The three reference/state fields are on the MODEL
      !! GAUSS GRID (nphi,nlat, rows north-first), exactly as the native solver
      !! receives them; this routine remaps them onto VILMA-v1's grid and lays down the
      !! file-based environment VILMA-v1 reads. VILMA-v1's own `setup` is deferred to the
      !! first advancing update (see the module header).
      type(vilma_v1_backend),  intent(inout) :: self
      type(vilma_param_class), intent(in)    :: par
      type(sht_grid),       intent(in)    :: sht
      real(wp),             intent(in)    :: z_bed_eq(:,:)  !! relaxed bedrock [m]
      real(wp),             intent(in)    :: h_ice_eq(:,:)  !! reference grounded ice [m]
      real(wp),             intent(in)    :: h_ice(:,:)     !! start-slice ice [m]
#ifndef VILMA_V1
      call vilma_v1_require()     ! aborts; the arguments are unused in the stub build
      ! (referenced so the stub compiles warning-free without touching VILMA-v1)
      if (.false.) print *, self%active, par%vilma_v1_jmax, sht%lmax, &
                            size(z_bed_eq), size(h_ice_eq), size(h_ice)
#else
      real(wp), allocatable :: zeq_v(:,:)

      call vilma_v1_finalize(self)               ! clean slate

      self%out_dir = trim(par%vilma_v1_out_dir)
      call require_dir(self%out_dir)
      call require_dir(trim(self%out_dir)//'/restart')
      ! Check the LONGEST path the backend will build, which lives under
      ! out_dir/restart/ -- not out_dir itself. The guard used to test
      ! out_dir//'/mos_acompl.nc' (od+14) while the restart scratch names are
      ! od+22, so an out_dir in that 8-character window passed the check and then
      ! truncated into VILMA-v1's character(len=120) fields -- which surfaces as a
      ! bare "No such file or directory", the exact failure this guard prevents.
      call check_path_len(trim(self%out_dir)//'/restart/nwl_struct.nc')

      ! --- grids + the model-Gauss <-> VILMA-v1 map pair --------------------------
      call build_maps(self, par, sht)

      allocate(self%work_v(self%nlon, self%nlat))
      allocate(self%h_ice_v(self%nlon, self%nlat))
      allocate(self%h_ref_v(self%nlon, self%nlat))
      allocate(zeq_v(self%nlon, self%nlat))

      ! bed is geometry (no mass rescale); ice carries mass (rescaled, see to_vilma_v1)
      call gauss_to_vilma_v1(self, sht, z_bed_eq, zeq_v,        is_ice=.false.)
      call gauss_to_vilma_v1(self, sht, h_ice_eq, self%h_ref_v, is_ice=.true.)
      call gauss_to_vilma_v1(self, sht, h_ice,    self%h_ice_v, is_ice=.true.)
      call warn_if_not_relaxed(self)

      ! --- VILMA-v1 run configuration (vg) ---------------------------------------
      ! Mirrors CLIMBER-X vilma_init. vg%dt / vg%btime / vg%etime are filled on the
      ! first advancing update, which is also when `setup` runs.
      vg%jmin       = 0
      vg%jmax       = par%vilma_v1_jmax
      ! 1 = elastic structure from a polynomial PREM (densi.inp then supplies only
      ! layer boundaries and element sizes); 0 = use densi.inp's rho/mu columns.
      vg%l_prem     = par%vilma_v1_l_prem
      vg%l_mod      = merge(1, 0, par%l_visc_3d)   ! 0 = 1-D radial, 1 = read 3-D field
      vg%l_toro     = 0        ! toroidal loading is irrelevant for GIA
      ! Follow par%rotation rather than hardwiring it on: otherwise rotation=.false.
      ! gives a non-rotating VILMA against a rotating VILMA-v1, the highest-order
      ! physics term differing silently in an v2-v1 pair.
      vg%l_rot      = merge(31, 0, par%rotation)   ! rotational variations in the potential
      vg%l_grid     = 2        ! loading supplied as a spatial grid
      vg%l_envonly  = 0
      vg%ntime      = 10000000 ! no cap: the driver's window decides
      vg%l_wepoch   = 1        ! output epochs listed in wepochs.inp
      vg%l_load_hist = .false. ! only the current slice is kept in the ice-history file
      vg%restart    = .false.  ! VILMA restarts are not wired to VILMA-v1's (see doc)

      ! --- VILMA-v1's file environment -------------------------------------------
      ! Only what does not depend on the (not-yet-known) coupling interval: the
      ! reference NetCDFs, the ice-history file, and the index naming them. The
      ! epoch-dependent pieces (loadh.inp, wepochs.inp, io.tmp) are written in
      ! start_vilma_v1, once the first interval is known.
      call set_io_paths(self, par)
      call write_reference_nc(self, zeq_v)
      call write_load_hist_index(self)

      self%t_kyr   = 0.0_wp
      self%nsub    = max(1, par%vilma_v1_nsub)
      self%started = .false.
      self%active  = .true.

      write(*,'(a)')             ' ======================================================='
      write(*,'(a)')             '  solid-earth backend: VILMA-v1'
      write(*,'(a,i0)')          '    jmax          = ', vg%jmax
      write(*,'(a,i0,a,i0)')     '    VILMA-v1 grid    = ', self%nlon, ' x ', self%nlat
      write(*,'(a,i0,a,i0)')     '    model Gauss   = ', self%nphi, ' x ', self%ngauss
      write(*,'(a,a)')           '    work dir      = ', trim(self%out_dir)
      write(*,'(a,l1)')          '    3-D viscosity = ', par%l_visc_3d
      write(*,'(a)')             '    (VILMA-v1 setup deferred to the first advancing step)'
      write(*,'(a)')             ' ======================================================='
#endif
   end subroutine vilma_v1_init

   subroutine vilma_v1_update(self, sht, h_ice, dt_yr, rsl_gauss, t_remap, t_solve)
      !! Advance VILMA-v1 by dt_yr years under the ice thickness h_ice (MODEL GAUSS
      !! GRID) and return the relative sea level on the MODEL GAUSS GRID. t_remap
      !! and t_solve are accumulators (seconds) for the two phases, kept separate so
      !! the coupling cost and the solver cost can be reported independently.
      !!
      !! dt_yr <= 0 is a seed step: the entering ice is recorded and rsl is left
      !! untouched, matching what a zero-length step does for the native solver.
      type(vilma_v1_backend), intent(inout) :: self
      type(sht_grid),      intent(in)    :: sht
      real(wp),            intent(in)    :: h_ice(:,:)      !! (nphi,nlat) [m]
      real(wp),            intent(in)    :: dt_yr           !! interval [years]
      real(wp),            intent(inout) :: rsl_gauss(:,:)  !! (nphi,nlat) [m], filled here
      real(wp),            intent(inout) :: t_remap, t_solve
#ifndef VILMA_V1
      call vilma_v1_require()     ! unreachable: init already aborted
      if (.false.) print *, self%active, sht%lmax, size(h_ice), dt_yr, &
                            size(rsl_gauss), t_remap, t_solve
#else
      integer(kind=8) :: pc0, pc1, prate

      if (.not. self%active) error stop 'vilma_v1_update: backend not initialised'

      ! ice onto VILMA-v1's grid (conservative + global mass match)
      call system_clock(pc0, prate)
      call gauss_to_vilma_v1(self, sht, h_ice, self%h_ice_v, is_ice=.true.)
      call system_clock(pc1);  t_remap = t_remap + real(pc1-pc0,wp)/prate

      if (dt_yr <= 0.0_wp) then
         ! Seed step: the entering ice is now recorded in h_ice_v and nothing is
         ! advanced. VILMA-v1's `setup` needs a real interval, so it stays deferred.
         return
      end if

      call system_clock(pc0, prate)
      if (.not. self%started) call start_vilma_v1(self, dt_yr)

      ! A fresh interval: drop the per-interval restart scratch VILMA-v1 would otherwise
      ! pick up (exactly as CLIMBER-X's vilma_update does).
      call delete_if_present(io_rsl_rs%n)
      call delete_if_present(io_dfgl_rs%n)

      ! VILMA-v1 re-reads the load from its ice-history NetCDF every step (that is what
      ! vg%l_load_hist = .false. buys), so the current slice is published first.
      ! Once per COUPLING interval, not per sub-step: the load is held across the
      ! sub-steps, which is VILMA-v1's own convention and the documented difference
      ! from the native solver's linear ramp within an interval.
      call write_ice_slice(self)
      block
        real(wp) :: dt_sub
        integer  :: ksub
        dt_sub = dt_yr/real(self%nsub, wp)
        do ksub = 1, self%nsub
           vg%btime = vg%etime
           vg%etime = vg%btime + real(dt_sub, kind(vg%etime))*1.0d-3
           call time_evolution
        end do
      end block
      self%t_kyr = real(vg%etime, wp)
      call system_clock(pc1);  t_solve = t_solve + real(pc1-pc0,wp)/prate

      if (size(rsl,1) /= self%nlat .or. size(rsl,2) /= self%nlon) then
         write(error_unit,'(a,i0,a,i0,a,i0,a,i0,a)') &
              ' vilma_v1: VILMA-v1 returned rsl(', size(rsl,1), ',', size(rsl,2), &
              ') but vilma_v1_grid_file describes a (', self%nlat, ',', self%nlon, ') grid'
         write(error_unit,'(a)') '   vilma_v1_jmax and vilma_v1_grid_file must describe the same grid.'
         flush(error_unit)
         error stop 'vilma_v1_update: VILMA-v1 grid does not match vilma_v1_grid_file'
      end if

      ! relative sea level back onto the model Gauss grid. VILMA-v1 stores rsl as
      ! (lat,lon) in double precision; transpose to (lon,lat) before remapping.
      call system_clock(pc0)
      self%work_v = transpose(real(rsl, wp))
      call vilma_v1_to_gauss(self, self%work_v, rsl_gauss)
      call system_clock(pc1);  t_remap = t_remap + real(pc1-pc0,wp)/prate
#endif
   end subroutine vilma_v1_update

   subroutine vilma_v1_finalize(self)
      !! Release the wrapper's state and close VILMA-v1 down. Safe on a fresh object.
      type(vilma_v1_backend), intent(inout) :: self
#ifdef VILMA_V1
      if (self%active .and. self%started) call close_evolution
#endif
      if (allocated(self%lon))     deallocate(self%lon)
      if (allocated(self%lat))     deallocate(self%lat)
      if (allocated(self%h_ice_v)) deallocate(self%h_ice_v)
      if (allocated(self%h_ref_v)) deallocate(self%h_ref_v)
      if (allocated(self%work_v))  deallocate(self%work_v)
      self%active  = .false.
      self%started = .false.
      self%nlon = 0;  self%nlat = 0;  self%nphi = 0;  self%ngauss = 0
      self%t_kyr = 0.0_wp
   end subroutine vilma_v1_finalize

#ifdef VILMA_V1
   ! ===========================================================================
   ! Internals — compiled only in a VILMA-v1 build.
   ! ===========================================================================

   subroutine start_vilma_v1(self, dt_yr)
      !! Finish the configuration that needs the coupling interval, write VILMA-v1's
      !! stdin file, and call `setup`. Runs once, on the first advancing update.
      !!
      !! vg%dt is fixed here from the FIRST coupling interval and written into
      !! io.tmp, because that is what VILMA-v1's setup consumes. A forcing axis with a
      !! non-uniform cadence therefore hands VILMA-v1 sub-steps it was not configured
      !! for. That is NOT currently detected: there is no cadence check anywhere in
      !! this wrapper (see doc/vilma-v1-backend.md).
      type(vilma_v1_backend), intent(inout) :: self
      real(wp),            intent(in)    :: dt_yr

      ! VILMA-v1's dt is in SECONDS, and it is the SUB-step: VILMA-v1 checks this value
      ! against the shortest Maxwell time in the structure at setup and aborts if
      ! it is too large. It has no sub-stepping of its own, so the coupling
      ! interval is divided here instead.
      vg%dt = real(dt_yr/real(self%nsub, wp)*sec_per_year, kind(vg%dt))
      ! Start of the transient. `setup` takes etime as the initial epoch and the
      ! 9999 sentinel in btime means "no restart" (CLIMBER-X uses the same values).
      vg%btime = 9999.0d0
      vg%etime = real(self%t_kyr, kind(vg%etime))

      call write_stdin_file()
      call write_output_epochs(self)
      call create_load_history(self, self%t_kyr, self%t_kyr + dt_yr*1.0e-3_wp)

      write(*,'(a,es10.3,a,i0,a,es10.3,a,f10.4,a)') ' VILMA-v1 setup: coupling dt =', dt_yr, &
           ' yr / ', self%nsub, ' sub-step(s) =', dt_yr/real(self%nsub, wp), &
           ' yr, start epoch =', self%t_kyr, ' kyr'
      call setup
      self%started = .true.
      write(*,'(a)') ' VILMA-v1 setup complete.'
   end subroutine start_vilma_v1

   subroutine build_maps(self, par, sht)
      !! Read VILMA-v1's lon/lat axes and build the two remap legs between the model
      !! Gauss grid and VILMA-v1's grid. Weights are cached by coords under "maps" and
      !! keyed by the grid names, which carry their dimensions so a resolution change
      !! invalidates the cache — the same scheme vilma_remap uses.
      type(vilma_v1_backend),  intent(inout) :: self
      type(vilma_param_class), intent(in)    :: par
      type(sht_grid),       intent(in)    :: sht
      real(wp), allocatable :: lon_g(:), lat_g(:)
      character(len=64) :: gname, vname
      logical :: ok
      integer :: j

      inquire(file=trim(par%vilma_v1_grid_file), exist=ok)
      if (.not. ok) then
         write(error_unit,'(a)') ' vilma_v1: vilma_v1_grid_file not found: '//trim(par%vilma_v1_grid_file)
         write(error_unit,'(a)') '   it must carry VILMA-v1''s own lon/lat axes (see doc/vilma-v1-backend.md)'
         flush(error_unit)
         error stop 'vilma_v1_init: vilma_v1_grid_file not found'
      end if
      self%nlon = nc_size(trim(par%vilma_v1_grid_file), "lon")
      self%nlat = nc_size(trim(par%vilma_v1_grid_file), "lat")
      call check_grid_vs_jmax(self%nlon, self%nlat, par%vilma_v1_jmax, trim(par%vilma_v1_grid_file))
      allocate(self%lon(self%nlon), self%lat(self%nlat))
      call nc_read(trim(par%vilma_v1_grid_file), "lon", self%lon)
      call nc_read(trim(par%vilma_v1_grid_file), "lat", self%lat)

      self%nphi = sht%nphi;  self%ngauss = sht%nlat
      allocate(lon_g(self%nphi), lat_g(self%ngauss))
      do j = 1, self%nphi
         lon_g(j) = sht%lon(j)*rad2deg                  ! SHTns longitudes [0,360)
      end do
      do j = 1, self%ngauss
         ! SHTns rows run north -> south; coords wants an ASCENDING latitude axis.
         lat_g(j) = 90.0_wp - sht%colat(self%ngauss - j + 1)*rad2deg
      end do

      call execute_command_line("mkdir -p 'maps'")
      write(gname, '(a,i0,a,i0)') "gauss_", self%nphi, "x", self%ngauss
      write(vname, '(a,i0,a,i0)') "vilma_v1_", self%nlon, "x", self%nlat
      call grid_init(self%gauss, name=trim(gname), mtype="latlon", units="degrees", &
                     lon180=.false., x=lon_g,    y=lat_g)
      call grid_init(self%vgrid, name=trim(vname), mtype="latlon", units="degrees", &
                     lon180=.true.,  x=self%lon, y=self%lat)

      call map_init(self%to_vilma_v1, self%gauss, self%vgrid, method="con",      gen="coords", fldr="maps")
      call map_init(self%to_gauss, self%vgrid, self%gauss, method="bilinear", &
                    max_neighbors=BIL_NEIGHBORS, gen="coords", fldr="maps")
   end subroutine build_maps

   subroutine gauss_to_vilma_v1(self, sht, f_gauss, f_vilma_v1, is_ice)
      !! Model Gauss grid -> VILMA-v1 grid, conservatively.
      !!
      !! is_ice=.true. additionally (a) remaps a 0/1 ice mask and zeroes cells whose
      !! mask falls below 0.5, so ice is not smeared across the coastline — the trick
      !! CLIMBER-X's wrapper uses — and (b) rescales by a single global factor so the
      !! ice mass VILMA-v1 receives equals the mass on the model Gauss grid as SHTns
      !! integrates it. (b) is a DELIBERATE difference from CLIMBER-X, which does not
      !! rescale: it makes the two backends see identical ice mass, which is the
      !! whole point of driving them from one driver.
      type(vilma_v1_backend), intent(in)  :: self
      type(sht_grid),      intent(in)  :: sht
      real(wp),            intent(in)  :: f_gauss(:,:)   !! (nphi,nlat), rows north-first
      real(wp),            intent(out) :: f_vilma_v1(:,:)   !! (nlon,nlat_v)
      logical,             intent(in)  :: is_ice
      real(wp), allocatable :: asc(:,:), mask_g(:,:), mask_v(:,:)
      logical,  allocatable :: m2(:,:)
      real(wp) :: src_sr, dst_sr
      integer  :: j

      allocate(asc(self%nphi, self%ngauss), m2(self%nlon, self%nlat))
      do j = 1, self%ngauss                     ! north-first -> ascending latitude
         asc(:, j) = f_gauss(:, self%ngauss - j + 1)
      end do
      call map_field(self%to_vilma_v1, "f", asc, f_vilma_v1, stat="mean", mask2=m2)
      where (.not. m2) f_vilma_v1 = 0.0_wp         ! uncovered cells (global source: none)

      if (.not. is_ice) return

      ! coastline-preserving mask, as in CLIMBER-X's vilma wrapper
      allocate(mask_g(self%nphi, self%ngauss), mask_v(self%nlon, self%nlat))
      where (asc > 0.0_wp)
         mask_g = 1.0_wp
      elsewhere
         mask_g = 0.0_wp
      end where
      mask_v = 1.0_wp
      call map_field(self%to_vilma_v1, "m", mask_g, mask_v, stat="mean", mask2=m2)
      where (mask_v < 0.5_wp) f_vilma_v1 = 0.0_wp

      ! global mass match against the SHTns quadrature on the source grid
      src_sr = sht_grid_surface_integral(sht, f_gauss)
      dst_sr = sum(f_vilma_v1 * self%vgrid%area) / sum(self%vgrid%area) * FOURPI
      if (abs(dst_sr) > tiny(1.0_wp)) f_vilma_v1 = f_vilma_v1 * (src_sr/dst_sr)
   end subroutine gauss_to_vilma_v1

   subroutine vilma_v1_to_gauss(self, f_vilma_v1, f_gauss)
      !! VILMA-v1 grid -> model Gauss grid, bilinearly (a smooth field: relative sea
      !! level). No mass rescale: vilma_coupling reconstructs the bed from it as
      !! z_bed_eq - rsl, exactly as it does for the native solver.
      type(vilma_v1_backend), intent(in)  :: self
      real(wp),            intent(in)  :: f_vilma_v1(:,:)   !! (nlon,nlat_v)
      real(wp),            intent(out) :: f_gauss(:,:)   !! (nphi,nlat), rows north-first
      real(wp), allocatable :: asc(:,:)
      logical,  allocatable :: m2(:,:)
      integer :: j
      allocate(asc(self%nphi, self%ngauss), m2(self%nphi, self%ngauss))
      call map_field(self%to_gauss, "f", f_vilma_v1, asc, method="bilinear", mask2=m2)
      where (.not. m2) asc = 0.0_wp
      do j = 1, self%ngauss                     ! ascending latitude -> north-first
         f_gauss(:, j) = asc(:, self%ngauss - j + 1)
      end do
   end subroutine vilma_v1_to_gauss

   ! --- VILMA-v1's file environment ----------------------------------------------

   subroutine set_io_paths(self, par)
      !! Point every VILMA-v1 input/output unit at a path under vilma_v1_input_dir (static
      !! inputs) or vilma_v1_out_dir (everything the run produces). One-for-one with
      !! CLIMBER-X's vilma_init, except that paths are always joined with "/" and
      !! the restart-matrix scratch lives under <out_dir>/restart rather than a
      !! separate restart-input tree (VILMA does not wire VILMA-v1 restarts).
      type(vilma_v1_backend),  intent(in) :: self
      type(vilma_param_class), intent(in) :: par
      character(len=:), allocatable :: id, od, rd
      id = trim(par%vilma_v1_input_dir)
      od = trim(self%out_dir)
      rd = od//'/restart'

      ! static inputs
      io_densi%n  = id//'/densi.inp'                          ! radial discretisation + density
      io_visko%n  = id//'/'//trim(par%vilma_v1_visc_1d_file)     ! 1-D radial viscosity
      io_nc3in%n  = id//'/'//trim(par%vilma_v1_visc_3d_file)     ! 3-D viscosity (l_mod=1 only)
      io_tint%n   = id//'/tint.inp'                           ! read only when vg%dt == 0
      io_sliin%n  = id//'/SLI_data.inp'

      ! run-generated control + output
      io_tmp%n     = od//'/io.tmp'          ! VILMA-v1's stdin, written by write_stdin_file
      io_lis%n     = od//'/vega.lis'        ! traceable log
      io_wepoch%n  = od//'/wepochs.inp'     ! epochs at which VILMA-v1 writes its own output
      io_hist%n    = od//'/load_hist.inp'   ! index of the ice-load-history files
      io_surf%n    = od//'/vega1.nc'
      io_radii%n   = od//'/radii.dat'
      io_nc3out%n  = od//'/visc3d.nc'
      io_sliout%n  = od//'/SLI_data.out'
      io_hsliout%n = od//'/hSLI_data.out'
      io_rsl%n     = od//'/rsl.nc'
      io_dfgl%n    = od//'/dflag.nc'
      io_rsl_rs%n  = od//'/rsl_rs.nc'
      io_dfgl_rs%n = od//'/dflg_rs.nc'
      io_deg1%n    = od//'/vega_deg1.dat'
      io_oce%n     = od//'/vega_oce.dat'
      io_rpt%n     = od//'/vega_rpt.dat'
      io_rslog%n   = od//'/restart.log'

      ! Galerkin-system / structure scratch (VILMA-v1's restart arrays)
      io_mos_indx%n   = rd//'/mos_indx.nc'
      io_mos_amtrx%n  = rd//'/mos_amtrx.nc'
      io_mos_acomp%n  = rd//'/mos_acomp.nc'
      io_mos_acompl%n = rd//'/mos_acompl.nc'
      io_ve_struct%n  = rd//'/ve_struct.nc'
      io_pjj%n        = rd//'/pjj.nc'
      io_pefgh%n      = rd//'/pefgh.nc'
      io_nwl_struct%n = rd//'/nwl_struct.nc'
      io_visc3drs%n   = rd//'/visc3d.nc'
      io_disp%n       = rd//'/disp.nc'
      io_stress%n     = rd//'/stress.nc'
      io_1dstress%n   = rd//'/ctc_stress.nc'
      ! io_rsrpt is the rotation-potential state, touched by setup / w_restart /
      ! r_restart alike. CLIMBER-X never assigns it, so it silently falls back to
      ! VILMA-v1's default relative path "restart/rotpot.log" and depends on the
      ! process working directory having such a folder. Point it at ours instead.
      io_rsrpt%n      = rd//'/rotpot.log'

      ! a stale output file from a previous run in the same directory would be
      ! appended to rather than replaced
      call delete_if_present(io_surf%n)
      call delete_if_present(io_rsl%n)
      call delete_if_present(io_dfgl%n)
   end subroutine set_io_paths

   subroutine write_stdin_file()
      !! VILMA-v1's parameter file (io.tmp), which replaces its interactive stdin. The
      !! order of the ten records is fixed by the library; see CLIMBER-X's vilma.F90
      !! for the same block with the same comments.
      integer :: u
      open(newunit=u, file=trim(io_tmp%n), form='formatted', status='replace', action='write')
      write(u,fmt=*) vg%jmin, vg%jmax     ! 1. spectral resolution
      write(u,fmt=*) vg%dt                ! 2. time step [s] (0 => read from tint.inp)
      write(u,fmt=*) vg%l_prem            ! 3. polynomial PREM earth structure
      write(u,fmt=*) vg%l_mod             ! 4. lateral viscosity: 0 no / 1 read / 2 adjust
      write(u,fmt=*) vg%l_toro            ! 5. toroidal loading
      write(u,fmt=*) vg%l_rot             ! 6. rotational variations
      write(u,fmt=*) vg%l_grid            ! 7. load grid type (2 = spatial grid)
      write(u,fmt=*) vg%l_envonly         ! 8. stop after reading the environment
      write(u,fmt=*) vg%ntime             ! 9. number of integration steps
      write(u,fmt=*) vg%l_wepoch          ! 10. output epochs from wepochs.inp
      close(u)
   end subroutine write_stdin_file

   subroutine write_load_hist_index(self)
      !! load_hist.inp: three lines naming (1) the epoch index file, (2) the ice
      !! history NetCDF, (3) the reference topography and reference ice NetCDFs.
      type(vilma_v1_backend), intent(in) :: self
      integer :: u
      open(newunit=u, file=trim(io_hist%n), form='formatted', status='replace', action='write')
      write(u,fmt='(A)')             trim(self%out_dir)//'/loadh.inp'
      write(u,fmt='(A)')             trim(self%out_dir)//'/vilma_h_ice.nc'
      write(u,fmt='(A)',advance='no') trim(self%out_dir)//'/vilma_z_bed_eq.nc'//' '// &
                                      trim(self%out_dir)//'/vilma_h_ice_eq.nc'
      close(u)
   end subroutine write_load_hist_index

   subroutine create_load_history(self, t0_kyr, t1_kyr)
      !! Create the ice-load history file VILMA-v1 integrates against, and the ASCII
      !! epoch index that describes it. Called once, from start_vilma_v1, when the first
      !! coupling interval is known.
      !!
      !! WHY IT LOOKS LIKE THIS. With vg%l_load_hist = .false. the history holds
      !! exactly TWO slices, and VILMA-v1 re-reads them unconditionally on every step:
      !!
      !!   slice 1, epoch t0 : the REFERENCE ice load. VILMA-v1 requires the first
      !!                       referenced load to vanish against the reference ice
      !!                       file (vilma_h_ice_eq.nc) and aborts otherwise with
      !!                       "first referenced load file should vanish but shows
      !!                       range of ..." — so this slice is fixed at the reference
      !!                       and never rewritten.
      !!   slice 2, epoch t1 : the CURRENT ice load, rewritten before every advance
      !!                       by write_ice_slice.
      !!
      !! The epoch coordinate is written once and left alone, exactly as CLIMBER-X's
      !! wrapper does: VILMA-v1's internal time marches past t1 after the first interval
      !! and it then keeps using the last record, i.e. the load is held at the current
      !! slice across each interval. That is a real difference from the native solver,
      !! which ramps the load linearly from the previous slice to the new one within
      !! the interval — see doc/vilma-v1-backend.md, "known differences".
      !!
      !! loadh.inp's first row is `nlat nlon rho_ice rho_ocean` — LATITUDE FIRST, and
      !! nlat/nlon must equal the NetCDF lat/lon sizes or VILMA-v1's check_dim_ne aborts.
      !! The densities are VILMA's OWN rho_ice / rho_water, so both backends
      !! turn the same ice thickness into the same load; CLIMBER-X hard-codes
      !! 910/1020 there instead.
      type(vilma_v1_backend), intent(in) :: self
      real(wp),            intent(in) :: t0_kyr, t1_kyr
      character(len=512) :: fnm
      integer :: ncid, u

      fnm = trim(self%out_dir)//'/vilma_h_ice.nc'
      call nc_create(fnm)
      call nc_open(fnm, ncid, writable=.TRUE.)
      call nc_write_dim(fnm, "epoch", x=[t0_kyr, t1_kyr], units="ka BP", unlimited=.TRUE., ncid=ncid)
      call nc_write_dim(fnm, "lon", x=self%lon, axis="x", ncid=ncid)
      call nc_write_dim(fnm, "lat", x=self%lat, axis="y", ncid=ncid)
      call nc_write(fnm, "Ice", self%h_ref_v, dims=["lon  ","lat  ","epoch"], &
                    start=[1,1,1], count=[self%nlon, self%nlat, 1], &
                    long_name="Ice thickness", units="m", ncid=ncid)
      call nc_write(fnm, "Ice", self%h_ice_v, dims=["lon  ","lat  ","epoch"], &
                    start=[1,1,2], count=[self%nlon, self%nlat, 1], &
                    long_name="Ice thickness", units="m", ncid=ncid)
      call nc_close(ncid)

      open(newunit=u, file=trim(self%out_dir)//'/loadh.inp', form='formatted', &
           status='replace', action='write')
      write(u,fmt=*) self%nlat, self%nlon, nint(rho_ice), nint(rho_water)
      write(u,fmt=*) 2
      write(u,fmt=*) t0_kyr
      write(u,fmt=*) t1_kyr
      close(u)
   end subroutine create_load_history

   subroutine write_output_epochs(self)
      !! wepochs.inp: the epochs at which VILMA-v1 writes its own diagnostic output.
      !! Only the run bounds are listed — the fields this coupling needs come back
      !! in memory (mod_sle's rsl), not through VILMA-v1's files.
      type(vilma_v1_backend), intent(in) :: self
      integer :: u
      open(newunit=u, file=trim(io_wepoch%n), form='formatted', status='replace', action='write')
      write(u,fmt=*) self%t_kyr
      write(u,fmt=*) self%t_kyr
      close(u)
   end subroutine write_output_epochs

   subroutine write_reference_nc(self, zeq_v)
      !! The relaxed reference state on VILMA-v1's grid, named by line 3 of
      !! load_hist.inp: the topography VILMA-v1 measures sea level against, and the
      !! reference ice load the history's first record must reproduce. Written once,
      !! at init, and never touched again.
      type(vilma_v1_backend), intent(in) :: self
      real(wp),            intent(in) :: zeq_v(:,:)
      character(len=512) :: fnm
      integer :: ncid

      fnm = trim(self%out_dir)//'/vilma_z_bed_eq.nc'
      call nc_create(fnm)
      call nc_open(fnm, ncid, writable=.TRUE.)
      call nc_write_dim(fnm, "epoch", x=0.0_wp, units="ka BP", unlimited=.TRUE., ncid=ncid)
      call nc_write_dim(fnm, "lon", x=self%lon, axis="x", ncid=ncid)
      call nc_write_dim(fnm, "lat", x=self%lat, axis="y", ncid=ncid)
      call nc_write(fnm, "topo", zeq_v, dims=["lon  ","lat  ","epoch"], &
                    start=[1,1,1], count=[self%nlon, self%nlat, 1], &
                    long_name="equilibrium bedrock elevation", units="m", ncid=ncid)
      call nc_close(ncid)

      fnm = trim(self%out_dir)//'/vilma_h_ice_eq.nc'
      call nc_create(fnm)
      call nc_open(fnm, ncid, writable=.TRUE.)
      call nc_write_dim(fnm, "epoch", x=0.0_wp, units="ka BP", unlimited=.TRUE., ncid=ncid)
      call nc_write_dim(fnm, "lon", x=self%lon, axis="x", ncid=ncid)
      call nc_write_dim(fnm, "lat", x=self%lat, axis="y", ncid=ncid)
      call nc_write(fnm, "Ice", self%h_ref_v, dims=["lon  ","lat  ","epoch"], &
                    start=[1,1,1], count=[self%nlon, self%nlat, 1], &
                    long_name="reference ice thickness", units="m", ncid=ncid)
      call nc_close(ncid)
   end subroutine write_reference_nc

   subroutine write_ice_slice(self)
      !! Publish the current ice load as slice 2 of the ice-history file (slice 1
      !! stays the reference load). This is what VILMA-v1 picks up when it integrates
      !! the next interval; CLIMBER-X does the same with l_load_hist = .false.
      type(vilma_v1_backend), intent(in) :: self
      character(len=512) :: fnm
      integer :: ncid
      fnm = trim(self%out_dir)//'/vilma_h_ice.nc'
      call nc_open(fnm, ncid, writable=.TRUE.)
      call nc_write(fnm, "Ice", self%h_ice_v, dims=["lon  ","lat  ","epoch"], &
                    start=[1,1,2], count=[self%nlon, self%nlat, 1], &
                    long_name="Ice thickness", units="m", ncid=ncid)
      call nc_close(ncid)
   end subroutine write_ice_slice

   ! --- small helpers ----------------------------------------------------------

   subroutine delete_if_present(path)
      character(len=*), intent(in) :: path
      integer :: u, stat
      if (len_trim(path) == 0) return
      open(newunit=u, iostat=stat, file=trim(path), status='old')
      if (stat == 0) close(u, status='delete')
   end subroutine delete_if_present

   subroutine require_dir(path)
      character(len=*), intent(in) :: path
      call execute_command_line("mkdir -p '"//trim(path)//"'")
   end subroutine require_dir

   subroutine warn_if_not_relaxed(self)
      !! VILMA-v1's load history MUST begin from the reference state: its first record
      !! is checked against the reference ice file and the run aborts with
      !!     "first referenced load file should vanish but shows range of ..."
      !! if the two differ. create_load_history therefore pins history slice 1 to the
      !! reference ice, which means VILMA-v1 starts with ZERO load anomaly and zero
      !! viscous memory at t0 — it has no way to be handed a pre-existing memory
      !! state, and solid_earth_spinup has no VILMA-v1 analogue (vilma_coupling refuses it).
      !!
      !! So if the run's start-slice ice is NOT the reference ice (the usual case for
      !! i_eq=1, a present-day reference with an LGM start), VILMA-v1 will absorb the
      !! entire start-vs-reference anomaly as a jump in the FIRST interval, while the
      !! native solver measures a genuine departure from a relaxed reference. That is
      !! a physics difference, not a bug, and the run must not pretend otherwise —
      !! hence this warning rather than silence or an abort.
      type(vilma_v1_backend), intent(in) :: self
      real(wp) :: dmax
      dmax = maxval(abs(self%h_ice_v - self%h_ref_v))
      if (dmax <= 1.0_wp) return
      write(*,'(a)')            ' -----------------------------------------------------------------'
      write(*,'(a,f0.1,a)')     ' WARNING (solver="v1"): start-slice ice differs from the reference'
      write(*,'(a,f0.1,a)')     '   ice by up to ', dmax, ' m.'
      write(*,'(a)')            '   VILMA-v1 requires its load history to start FROM the reference state,'
      write(*,'(a)')            '   so it begins with zero load anomaly and zero viscous memory and will'
      write(*,'(a)')            '   take the whole start-vs-reference difference as a jump in the first'
      write(*,'(a)')            '   interval. The native solver does not. For a like-for-like comparison'
      write(*,'(a)')            '   use &ctl i_eq=0 (the start slice IS the reference) or start the'
      write(*,'(a)')            '   window where the two agree. See doc/vilma-v1-backend.md.'
      write(*,'(a)')            ' -----------------------------------------------------------------'
   end subroutine warn_if_not_relaxed

   subroutine check_grid_vs_jmax(nlon, nlat, jmax, fname)
      !! VILMA-v1 derives its Gauss-Legendre working grid from vg%jmax alone: nlon is the
      !! smallest power of two strictly greater than 3*jmax, and nlat = nlon/2 (so
      !! jmax=170 gives the N128 grid, 512 x 256). The spatial grid of the load and of
      !! `rsl`, on the other hand, comes from what WE declare in loadh.inp and supply
      !! in the NetCDFs. If the two disagree, VILMA-v1 either aborts deep inside
      !! check_dim_ne or, worse, quietly works on mismatched fields — so check here,
      !! where the message can name the two settings that must agree.
      integer,          intent(in) :: nlon, nlat, jmax
      character(len=*), intent(in) :: fname
      integer :: want_lon, i
      want_lon = 2
      do i = 1, 30
         if (want_lon > 3*jmax) exit
         want_lon = want_lon*2
      end do
      if (nlon == want_lon .and. nlat == want_lon/2) return
      write(error_unit,'(a)')      ' vilma_v1: VILMA-v1 grid mismatch.'
      write(error_unit,'(a,i0,a,i0,a,i0)') '   vilma_v1_jmax = ', jmax, &
           ' implies a Gauss grid of ', want_lon, ' x ', want_lon/2
      write(error_unit,'(a,i0,a,i0,a)')    '   but vilma_v1_grid_file describes ', nlon, ' x ', nlat, ':'
      write(error_unit,'(a)')      '     '//trim(fname)
      write(error_unit,'(a)')      '   Set &vilma vilma_v1_jmax and vilma_v1_grid_file consistently'
      write(error_unit,'(a)')      '   (nlon = smallest power of 2 > 3*jmax, nlat = nlon/2).'
      flush(error_unit)
      error stop 'vilma_v1_init: vilma_v1_jmax and vilma_v1_grid_file describe different grids'
   end subroutine check_grid_vs_jmax

   subroutine check_path_len(path)
      !! VILMA-v1 stores its file names in character(len=120); a longer path would be
      !! silently truncated into a file it cannot open. Fail loudly instead.
      character(len=*), intent(in) :: path
      if (len_trim(path) <= VILMA_V1_PATH_LEN) return
      write(error_unit,'(a,i0,a)') ' vilma_v1: vilma_v1_out_dir is too long — VILMA-v1 stores file names in ', &
           VILMA_V1_PATH_LEN, ' characters and this would be truncated:'
      write(error_unit,'(a)') '   '//trim(path)
      write(error_unit,'(a)') '   Set &vilma vilma_v1_out_dir to a shorter (e.g. relative) path.'
      flush(error_unit)
      error stop 'vilma_v1_init: vilma_v1_out_dir path too long for VILMA-v1'
   end subroutine check_path_len
#endif

end module vilma_v1
