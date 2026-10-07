module vilma_io
   !! netCDF input/output for the coupling state, following the yelmo convention:
   !! a general, table-driven write_step writes a chosen subset of variables to a
   !! file with an unlimited time axis, so several snapshots accumulate in one
   !! file (handy for debugging). A restart file is just the full variable set
   !! written this way; vilma_restart_read restores the prognostic state.
   !!
   !! The variable registry — names, dimensions, units, long_names — lives in a
   !! markdown table (input/vilma-variables.md) parsed by variable_io. The
   !! table is both the I/O metadata source and the human documentation.
   !!
   !! Prognostic (restored on restart): the Maxwell memory-stress fields tau_*
   !! (nlam,ne,nk), the model time, and the adaptive controller's next-step Δt
   !! suggestion dt_try (so a restarted run sub-steps the same path → bit-for-bit
   !! continuation). With rotation enabled, the polar motion m and both degree-2
   !! channels' memory stress (rot_*) are persisted too. Static (written once,
   !! checked on read): the reference state z_bed_eq, h_ice_eq. Diagnostic: h_ice,
   !! rsl, z_bed, C_ocean (lon,lat) — written for inspection and restored if present.
   use vilma_precision,    only: wp
   use vilma_constants,    only: rad2deg, sec_per_year
   use vilma_viscoelastic, only: NLAM, NLAM_TOR
   use vilma_response,     only: response_prime_sigma, response, response_init_elastic, response_init_ve, &
                              response_init_null, RESP_VE, &
                              response_horizontal, response_horizontal_toroidal
   use vilma_sht,          only: sht_grid_sph_synthesis, sht_grid_tor_synthesis
   use vilma_rotation,     only: rotation_ne, rotation_get_memory, rotation_set_memory, ROT_NCOMP
   use vilma_coupling,     only: solid_earth, solid_earth_sync_host
   use ncio
   use variable_io
   implicit none
   private

   public :: vilma_restart_write, vilma_restart_read, vilma_write_step, vilma_io_set_table
   public :: vilma_write_horizontal

   ! The full time-varying variable set a restart writes, split by response kind.
   ! COMMON_VARS apply to every kind; the prognostic memory differs: RESP_VE carries
   ! the Maxwell stress tensor + trapezoidal σ_n; RESP_ELASTIC/RESP_NULL are
   ! memoryless (COMMON_VARS only).
   ! Static reference fields are written once in the init branch.
   character(len=12), parameter :: COMMON_VARS(6) = [character(len=12) :: &
        "h_ice", "rsl", "z_bed", "C_ocean", "dt_try", "bsl"]
   character(len=12), parameter :: VE_MEM_VARS(9) = [character(len=12) :: &
        "tau_a_re", "tau_a_im", "tau_b_re", "tau_b_im", "tau_c_re", "tau_c_im", &
        "sigma_n_re", "sigma_n_im", "sigma_primed"]

   character(len=256), save              :: table_file = "input/vilma-variables.md"
   type(var_io_type), allocatable, save  :: vtable(:)

contains

   subroutine vilma_io_set_table(filename)
      !! Override the variable-io table path (default input/vilma-variables.md)
      !! and force a reload on next use.
      character(len=*), intent(in) :: filename
      table_file = filename
      if (allocated(vtable)) deallocate(vtable)
   end subroutine vilma_io_set_table

   subroutine ensure_table()
      if (.not. allocated(vtable)) call load_var_io_table(vtable, trim(table_file))
   end subroutine ensure_table

   ! --- writing ---------------------------------------------------------------

   subroutine vilma_restart_write(self, time, filename, folder, init)
      !! Write the full coupling state as a restart snapshot at `time` to
      !! trim(folder)//"/"//trim(filename) (folder is created if needed). Defaults:
      !! folder="." and filename="vilma_restart.nc". With init=.true. (default) the file
      !! is (re)created; with init=.false. the snapshot is appended at a new time index.
      type(solid_earth), intent(in) :: self
      real(wp),           intent(in) :: time
      character(len=*), optional, intent(in) :: filename, folder
      logical, optional,  intent(in) :: init
      character(len=:), allocatable :: fn, fd
      fn = "vilma_restart.nc";  if (present(filename)) fn = trim(filename)
      fd = ".";              if (present(folder))   fd = trim(folder)
      call execute_command_line("mkdir -p '"//fd//"'")
      call vilma_write_step(self, fd//"/"//fn, time, init=init)
   end subroutine vilma_restart_write

   subroutine vilma_write_step(self, filename, time, nms, init)
      !! General step writer. Writes the variables named in nms (default: all
      !! time-varying variables) at the time index for `time`, appending along
      !! the unlimited time axis. The model time should be passed as `time` for
      !! a restart to round-trip the clock.
      type(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename
      real(wp),           intent(in) :: time
      character(len=*), optional, intent(in) :: nms(:)
      logical,          optional, intent(in) :: init
      logical :: do_init
      integer :: ncid, n, q
      real(wp) :: tyr

      call ensure_table()
      do_init = .true.;  if (present(init)) do_init = init

      ! the coupling clock (se%time) and this output time axis are both in YEARS.
      tyr = time

      if (do_init) call write_init(self, filename, time)

      call nc_open(filename, ncid, writable=.true.)
      if (do_init) then
         n = 1                                    ! the file was created with slice 1
      else
         n = nc_time_index(filename, "time", tyr, ncid)
         call nc_write(filename, "time", tyr, dim1="time", start=[n], count=[1], ncid=ncid)
      end if

      if (present(nms)) then
         do q = 1, size(nms);  call write_one(self, filename, trim(nms(q)), n, ncid);  end do
      else
         block
            character(len=12), allocatable :: vars(:)
            select case (self%resp%kind)
            case (RESP_VE);    vars = [VE_MEM_VARS,    COMMON_VARS]
            case default;      vars = COMMON_VARS            ! elastic / null: memoryless
            end select
            do q = 1, size(vars);  call write_one(self, filename, trim(vars(q)), n, ncid);  end do
         end block
         ! rotation prognostic state (m + both channels' memory), if enabled
         if (self%rotation%enabled) call write_rotation(self, filename, n, ncid)
      end if

      call nc_close(ncid)
   end subroutine vilma_write_step

   subroutine vilma_write_horizontal(self, filename, time, init)
      !! Surface horizontal displacement [m] at `time`, into its own file (&ctl
      !! file_hor): the total, spheroidal + toroidal, and the toroidal part alone,
      !! as east/north components on the Gauss grid of file_out. Relative to the
      !! reference state, like rsl: the response evaluated at the converged load of
      !! the last interval (se%sigma_lm) with the relaxation drift of its last
      !! sub-step — the same state rsl and z_bed were solved in.
      !!
      !! The toroidal part is zero, exactly, for every response that carries no
      !! toroidal field (radially symmetric or laterally uniform viscosity,
      !! elastic); it is written anyway so the file has one layout.
      !!
      !! East/north from the colatitude-ordered (θ, φ) components: e_φ is east,
      !! e_θ points south, so north = −u_θ.
      type(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename
      real(wp),           intent(in) :: time      !! [years]
      logical,            intent(in) :: init      !! create the file (first slice)
      complex(wp), allocatable :: v_lm(:), t_lm(:)
      real(wp),    allocatable :: sth(:,:), sph(:,:), tth(:,:), tph(:,:)
      real(wp),    allocatable :: lon_deg(:), lat_deg(:)
      integer :: ncid, n, np, nl
      if (self%use_vilma_v1) error stop 'vilma_write_horizontal: the VILMA1 backend returns no horizontal field'
      call ensure_table()
      np = self%sht%nphi;  nl = self%sht%nlat
      allocate(v_lm(self%sht%nlm), t_lm(self%sht%nlm))
      allocate(sth(np,nl), sph(np,nl), tth(np,nl), tph(np,nl))
      call response_horizontal(self%resp, self%sht, self%sigma_lm, v_lm)
      call response_horizontal_toroidal(self%resp, self%sht, t_lm)
      call sht_grid_sph_synthesis(self%sht, v_lm, sth, sph)     ! ∇₁V: (u_θ, u_φ)
      call sht_grid_tor_synthesis(self%sht, t_lm, tth, tph)     ! e_r×∇₁W
      if (init) then
         lon_deg = self%sht%lon * rad2deg
         lat_deg = 90.0_wp - self%sht%colat * rad2deg
         call nc_create(filename, overwrite=.true.)
         call nc_write_dim(filename, "lon",  x=lon_deg, units="degrees_east")
         call nc_write_dim(filename, "lat",  x=lat_deg, units="degrees_north")
         call nc_write_dim(filename, "time", x=time, dx=1.0_wp, nx=1, units="years", unlimited=.true.)
      end if
      call nc_open(filename, ncid, writable=.true.)
      if (init) then
         n = 1
      else
         n = nc_time_index(filename, "time", time, ncid)
         call nc_write(filename, "time", time, dim1="time", start=[n], count=[1], ncid=ncid)
      end if
      call put_hor(filename, "u_east",      sph + tph,    n, ncid, np, nl)
      call put_hor(filename, "u_north",     -(sth + tth), n, ncid, np, nl)
      call put_hor(filename, "u_east_tor",  tph,          n, ncid, np, nl)
      call put_hor(filename, "u_north_tor", -tth,         n, ncid, np, nl)
      call nc_close(ncid)
   end subroutine vilma_write_horizontal

   subroutine put_hor(filename, name, dat, n, ncid, np, nl)
      character(len=*), intent(in) :: filename, name
      real(wp),         intent(in) :: dat(:,:)
      integer,          intent(in) :: n, ncid, np, nl
      type(var_io_type) :: v
      call find_var_io_in_table(v, name, vtable, with_error=.true.)
      call nc_write(filename, name, dat, ncid=ncid, dim1="lon", dim2="lat", dim3="time", &
           start=[1,1,n], count=[np, nl, 1], units=trim(v%units), long_name=trim(v%long_name))
   end subroutine put_hor

   subroutine write_rotation(self, filename, n, ncid)
      !! Write the rotation solver's prognostic state at time slice n: the polar
      !! motion m (two scalars) and both degree-2 channels' packed memory stress
      !! (NLAM, ne_rot, ROT_NCOMP). Serialization is owned by vilma_rotation; this just
      !! moves the packed arrays into the file.
      type(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename
      integer,            intent(in) :: n, ncid
      complex(wp) :: m
      real(wp), allocatable :: load_mem(:,:,:), tidal_mem(:,:,:)
      integer :: ne
      ne = rotation_ne(self%rotation)
      allocate(load_mem(NLAM, ne, ROT_NCOMP), tidal_mem(NLAM, ne, ROT_NCOMP))
      call rotation_get_memory(self%rotation, m, load_mem, tidal_mem)
      call put_scalar(filename, "rot_m_re", real(m, wp),  n, ncid)
      call put_scalar(filename, "rot_m_im", aimag(m),     n, ncid)
      call put_rotmem(filename, "rot_load_mem",  load_mem,  ne, n, ncid)
      call put_rotmem(filename, "rot_tidal_mem", tidal_mem, ne, n, ncid)
   end subroutine write_rotation

   subroutine put_rotmem(filename, name, dat, ne, n, ncid)
      !! Write a packed degree-2 channel memory field (nlam, ne_rot, nrc) at slice n.
      character(len=*), intent(in) :: filename, name
      real(wp),         intent(in) :: dat(:,:,:)
      integer,          intent(in) :: ne, n, ncid
      type(var_io_type) :: v
      call find_var_io_in_table(v, name, vtable, with_error=.true.)
      call nc_write(filename, name, dat, ncid=ncid, &
           dim1="nlam_rot", dim2="ne_rot", dim3="nrc", dim4="time", &
           start=[1,1,1,n], count=[NLAM, ne, ROT_NCOMP, 1], &
           units=trim(v%units), long_name=trim(v%long_name))
   end subroutine put_rotmem

   subroutine write_init(self, filename, time)
      !! Create the file and its dimensions (lon, lat, nlam, ne, nk, time) and
      !! write the static reference fields.
      type(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename
      real(wp),           intent(in) :: time
      real(wp), allocatable :: lon_deg(:), lat_deg(:)

      lon_deg = self%sht%lon * rad2deg
      lat_deg = 90.0_wp - self%sht%colat * rad2deg

      call nc_create(filename, overwrite=.true.)
      call nc_write_dim(filename, "lon",  x=lon_deg, units="degrees_east")
      call nc_write_dim(filename, "lat",  x=lat_deg, units="degrees_north")
      ! prognostic-memory dimensions depend on the response kind
      select case (self%resp%kind)
      case (RESP_VE)
         ! nlam = memory channels the response carries: NLAM, or NLAM+NLAM_TOR once
         ! a genuinely 3-D element has switched the toroidal field on
         call nc_write_dim(filename, "nlam", x=1, dx=1, nx=self%resp%nlam, units="1")
         call nc_write_dim(filename, "ne",   x=1, dx=1, nx=self%resp%ne,  units="1")
         ! nk = # deforming (l>=1) coefficients, in the ve_response degree-grouped order
         call nc_write_dim(filename, "nk",   x=1, dx=1, nx=self%resp%nk,  units="1")
         ! nlm = # spherical-harmonic (l,m) coefficients (carries the σ_n load vector)
         call nc_write_dim(filename, "nlm",  x=1, dx=1, nx=self%resp%nlm, units="1")
      end select
      ! rotation memory dimensions (independent of the response kind / lmax). The
      ! packed channel memory is (nlam_rot, ne_rot, nrc): degree-2 and radially
      ! symmetric, so always spheroidal, and its own dimension now that the
      ! response's nlam can be wider. Files written before the split carry it on
      ! "nlam"; the read goes by variable and count, so both forms restore.
      if (self%rotation%enabled) then
         call nc_write_dim(filename, "nlam_rot", x=1, dx=1, nx=NLAM, units="1")
         call nc_write_dim(filename, "ne_rot", x=1, dx=1, nx=rotation_ne(self%rotation), units="1")
         call nc_write_dim(filename, "nrc",    x=1, dx=1, nx=ROT_NCOMP,                  units="1")
      end if
      call nc_write_dim(filename, "time", x=time, dx=1.0_wp, nx=1, &
           units="years", unlimited=.true.)

      call put2d(self, filename, "z_bed_eq",  self%gg%z_bed_eq,  static=.true.)
      call put2d(self, filename, "h_ice_eq", self%gg%h_ice_eq, static=.true.)

      call write_visc_attrs(self, filename)
   end subroutine write_init

   subroutine write_visc_attrs(self, filename)
      !! Stamp the lateral-viscosity (3-D) configuration into the file as global
      !! attributes, so a restart/output file is self-certifying for reproducibility
      !! (independent of the namelist): whether 3-D was active, how many radial
      !! elements were genuinely laterally-varying, the 1-D/3-D split tolerance, and
      !! an order-sensitive checksum of the absolute log10(eta) field that was loaded.
      type(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename
      call nc_write_attr(filename, "l_visc_3d", merge(1, 0, self%resp%lat_visc))
      call nc_write_attr(filename, "ne3d",       self%resp%ne3d)
      call nc_write_attr(filename, "ne",         self%resp%ne)
      call nc_write_attr(filename, "visc3d_tol", real(self%resp%visc3d_tol, wp))
      if (allocated(self%earth%visc_3d)) &
         call nc_write_attr(filename, "visc_3d_checksum", visc_3d_checksum(self%earth%visc_3d))
   end subroutine write_visc_attrs

   pure function visc_3d_checksum(v) result(chk)
      !! Order-sensitive fingerprint of the (node, radius) absolute log10(eta) field:
      !! sum of v weighted by linear index, so two distinct fields (or a permuted one)
      !! give different values. A bit-identical field reproduces the value exactly.
      real(wp), intent(in) :: v(:,:)
      real(wp) :: chk
      integer  :: i, j, n1
      n1  = size(v, 1)
      chk = 0.0_wp
      do j = 1, size(v, 2)
         do i = 1, n1
            chk = chk + v(i,j) * real(i + n1*(j-1), wp)
         end do
      end do
   end function visc_3d_checksum

   subroutine write_one(self, filename, name, n, ncid)
      !! Dispatch a single variable name to its array + dimensions.
      type(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename, name
      integer,            intent(in) :: n, ncid
      select case (name)
      case ("tau_a_re"); call put3d(self, filename, name, self%resp%Are, n, ncid)
      case ("tau_a_im"); call put3d(self, filename, name, self%resp%Aim, n, ncid)
      case ("tau_b_re"); call put3d(self, filename, name, self%resp%Bre, n, ncid)
      case ("tau_b_im"); call put3d(self, filename, name, self%resp%Bim, n, ncid)
      case ("tau_c_re"); call put3d(self, filename, name, self%resp%Cre, n, ncid)
      case ("tau_c_im"); call put3d(self, filename, name, self%resp%Cim, n, ncid)
      case ("h_ice");    call put2d(self, filename, name, self%gg%h_ice, n=n, ncid=ncid)
      case ("rsl");      call put2d(self, filename, name, self%gg%rsl,   n=n, ncid=ncid)
      case ("z_bed");    call put2d(self, filename, name, self%gg%z_bed, n=n, ncid=ncid)
      case ("C_ocean");  call put2d(self, filename, name, self%gg%C,     n=n, ncid=ncid)
      case ("bsl");      call put_scalar(filename, name, self%bsl, n, ncid)
      ! Polar motion m = m1 + i*m2 [rad]. Already prognostic and already
      ! serialized into restarts by write_rotation; these two cases let the
      ! DIAGNOSTIC output file carry it as well, via an explicit nms list.
      ! Only meaningful when the rotation solver ran -- vilma_drive adds them to
      ! its list only if se%rotation%enabled, so a rotation=.false. run and the
      ! VILMA1 backend (which keeps its own internal rotation and never touches
      ! this state) omit the variables rather than writing a fabricated zero.
      case ("rot_m_re"); call put_scalar(filename, name, real(self%rotation%m, wp), n, ncid)
      case ("rot_m_im"); call put_scalar(filename, name, aimag(self%rotation%m),    n, ncid)
      case ("dt_try");   call put_scalar(filename, name, self%stepper%dt_try/sec_per_year, n, ncid)
      case ("sigma_n_re"); call put_sigma(self, filename, name, want_re=.true.,  n=n, ncid=ncid)
      case ("sigma_n_im"); call put_sigma(self, filename, name, want_re=.false., n=n, ncid=ncid)
      case ("sigma_primed")
         call put_scalar(filename, name, &
              merge(1.0_wp, 0.0_wp, allocated(self%resp%sigma_n) .and. self%resp%sigma_primed), &
              n, ncid)
      case default
         error stop 'vilma_io: unknown variable "'//trim(name)//'"'
      end select
   end subroutine write_one

   subroutine put3d(self, filename, name, dat, n, ncid)
      !! Write a (nlam,ne,nk) Maxwell-memory field at time slice n.
      type(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename, name
      real(wp),           intent(in) :: dat(:,:,:)
      integer,            intent(in) :: n, ncid
      type(var_io_type) :: v
      call find_var_io_in_table(v, name, vtable, with_error=.true.)
      call nc_write(filename, name, dat, ncid=ncid, &
           dim1="nlam", dim2="ne", dim3="nk", dim4="time", &
           start=[1,1,1,n], count=[self%resp%nlam, self%resp%ne, self%resp%nk, 1], &
           units=trim(v%units), long_name=trim(v%long_name))
   end subroutine put3d

   subroutine put_sigma(self, filename, name, want_re, n, ncid)
      !! Write the real or imaginary part of the trapezoidal start-of-step load σ_n
      !! (a spectral (nlm) vector) at time slice n; zeros if σ_n is not yet tracked.
      type(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename, name
      logical,            intent(in) :: want_re
      integer,            intent(in) :: n, ncid
      type(var_io_type)     :: v
      real(wp), allocatable :: dat(:)
      integer :: nlm
      nlm = self%resp%nlm
      allocate(dat(nlm))
      if (allocated(self%resp%sigma_n)) then
         if (want_re) then;  dat = real(self%resp%sigma_n, wp)
         else;               dat = aimag(self%resp%sigma_n)
         end if
      else
         dat = 0.0_wp
      end if
      call find_var_io_in_table(v, name, vtable, with_error=.true.)
      call nc_write(filename, name, dat, ncid=ncid, dim1="nlm", dim2="time", &
           start=[1,n], count=[nlm,1], units=trim(v%units), long_name=trim(v%long_name))
   end subroutine put_sigma

   subroutine put_scalar(filename, name, val, n, ncid)
      !! Write a single time-varying scalar (controller state) at time slice n.
      character(len=*), intent(in) :: filename, name
      real(wp),         intent(in) :: val
      integer,          intent(in) :: n, ncid
      type(var_io_type) :: v
      call find_var_io_in_table(v, name, vtable, with_error=.true.)
      call nc_write(filename, name, val, ncid=ncid, dim1="time", &
           start=[n], count=[1], units=trim(v%units), long_name=trim(v%long_name))
   end subroutine put_scalar

   subroutine put2d(self, filename, name, dat, n, ncid, static)
      !! Write a (lon,lat) field — static (no time) or at time slice n.
      type(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename, name
      real(wp),           intent(in) :: dat(:,:)
      integer, optional,  intent(in) :: n, ncid
      logical, optional,  intent(in) :: static
      type(var_io_type) :: v
      logical :: is_static
      is_static = .false.;  if (present(static)) is_static = static
      call find_var_io_in_table(v, name, vtable, with_error=.true.)
      if (is_static) then
         call nc_write(filename, name, dat, dim1="lon", dim2="lat", &
              units=trim(v%units), long_name=trim(v%long_name))
      else
         call nc_write(filename, name, dat, ncid=ncid, &
              dim1="lon", dim2="lat", dim3="time", &
              start=[1,1,n], count=[self%sht%nphi, self%sht%nlat, 1], &
              units=trim(v%units), long_name=trim(v%long_name))
      end if
   end subroutine put2d

   ! --- reading ---------------------------------------------------------------

   subroutine vilma_restart_read(self, filename, time)
      !! Restore the prognostic state into an already-initialised solid_earth
      !! (init() must have rebuilt the operators, spectrum and grid). Reads the
      !! response's memory state + model time at the time slice matching `time`
      !! (default: the last), dispatching on the response kind:
      !!
      !!   RESP_VE    — the Maxwell memory-stress tensor tau_* (+ trapezoidal σ_n).
      !!                Same-resolution restore is exact (bit-for-bit continuation);
      !!                a lower-resolution file upsamples via the degree-grouped block.
      !!   elastic/null — memoryless: only the diagnostics + clock are restored.
      !!
      !! The clock and the adaptive controller's dt_try seed are restored for every
      !! kind so the restarted run sub-steps the same path as the uninterrupted one.
      type(solid_earth), intent(inout) :: self
      character(len=*),   intent(in)    :: filename
      real(wp), optional, intent(in)    :: time
      real(wp), allocatable :: tvals(:)
      integer :: nt, n, np, nl

      np = self%sht%nphi;  nl = self%sht%nlat

      ! select the time slice. The file stores YEARS, matching the coupling clock
      ! self%time and the `time` selector argument.
      nt = nc_size(filename, "time")
      allocate(tvals(nt));  call nc_read(filename, "time", tvals)
      if (present(time)) then
         n = minloc(abs(tvals - time), dim=1)
         if (abs(tvals(n) - time) > 1.0e-6_wp*max(abs(time), 1.0_wp)) &
            error stop 'vilma_restart_read: requested time not present in file'
      else
         n = nt
      end if

      ! kind-specific prognostic memory + (where applicable) reference/diagnostics
      select case (self%resp%kind)
      case (RESP_VE);    call read_ve_state(self, filename, n, np, nl)
      case default;      call read_diagnostics(self, filename, n, np, nl)  ! memoryless
      end select

      ! restore the clock and the adaptive controller's step-size seed. self%time is
      ! in years (file units); the response's internal clock is in SI seconds.
      self%resp%time = tvals(n)*sec_per_year
      self%time      = tvals(n)
      call nc_read(filename, "dt_try", self%stepper%dt_try, start=[n], count=[1])
      self%stepper%dt_try = self%stepper%dt_try*sec_per_year   ! file stores years

      ! rotation prognostic state, if this run has rotation on AND the file carries
      ! it (an older, rotation-off restart lacks the rot_* vars → cold-start rotation)
      if (self%rotation%enabled .and. nc_exists_var(filename, "rot_m_re")) &
         call read_rotation(self, filename, n)

      ! the restart restores the Gauss-grid state (gg) + memory; re-derive the
      ! host-grid outputs (rsl, z_bed) from it so the host reads a consistent state.
      call solid_earth_sync_host(self)
   end subroutine vilma_restart_read

   subroutine read_rotation(self, filename, n)
      !! Restore the rotation solver's prognostic state at slice n: the polar motion
      !! m and both degree-2 channels' memory stress. The channels are rebuilt by
      !! init() (so ne_rot matches and is lmax-independent); only m + memory restore.
      type(solid_earth), intent(inout) :: self
      character(len=*),   intent(in)    :: filename
      integer,            intent(in)    :: n
      complex(wp) :: m
      real(wp), allocatable :: load_mem(:,:,:), tidal_mem(:,:,:)
      real(wp) :: mre, mim
      integer  :: ne, ne_f
      ne   = rotation_ne(self%rotation)
      ne_f = nc_size(filename, "ne_rot")
      if (ne_f /= ne) &
         error stop 'vilma_restart_read: rotation channel size (ne_rot) does not match the model'
      allocate(load_mem(NLAM, ne, ROT_NCOMP), tidal_mem(NLAM, ne, ROT_NCOMP))
      call nc_read(filename, "rot_load_mem",  load_mem,  start=[1,1,1,n], count=[NLAM, ne, ROT_NCOMP, 1])
      call nc_read(filename, "rot_tidal_mem", tidal_mem, start=[1,1,1,n], count=[NLAM, ne, ROT_NCOMP, 1])
      call nc_read(filename, "rot_m_re", mre, start=[n], count=[1])
      call nc_read(filename, "rot_m_im", mim, start=[n], count=[1])
      m = cmplx(mre, mim, wp)
      call rotation_set_memory(self%rotation, m, load_mem, tidal_mem)
      self%rotation%time = self%time
   end subroutine read_rotation

   subroutine read_ve_state(self, filename, n, np, nl)
      !! Restore RESP_VE memory at slice n. Same-resolution (file nk == model nk):
      !! exact restore of the tau_* tensor, σ_n and diagnostics, validated against
      !! the model reference, so the run continues bit-for-bit. Cross-resolution
      !! (file nk < model nk, e.g. an l64 spin-up into an l128 run): the memory is
      !! degree-grouped, so the shared low-degree block is contiguous — copy it and
      !! zero-pad the higher degrees (upsampling only); σ_n and the spatial diagnostics
      !! are NOT transferred (different grid), so the caller re-seeds with a zero-Δt solve.
      type(solid_earth), intent(inout) :: self
      character(len=*),   intent(in)    :: filename
      integer,            intent(in)    :: n, np, nl
      integer :: ne, nk, nk_f, nl_f, L
      logical :: cross_res, ok_block

      ne = self%resp%ne;  nk = self%resp%nk
      ! the radial mesh (ne) must always match; the horizontal resolution (nk) may
      ! differ -> cross-resolution restart.
      if (nc_size(filename, "ne") /= ne) &
         error stop 'vilma_restart_read: radial mesh (ne) does not match the model'
      ! Memory channels. A spheroidal file (NLAM) restores into a run that carries
      ! the toroidal channels too: they start at zero, which is exactly the
      ! toroidal state of every run that has only ever been spheroidal. The
      ! reverse would drop toroidal memory a laterally varying run built up, so
      ! it is refused rather than truncated.
      nl_f = nc_size(filename, "nlam")
      if (nl_f /= self%resp%nlam .and. .not. (nl_f == NLAM .and. self%resp%nlam == NLAM + NLAM_TOR)) then
         write(*,'(a,i0,a,i0,a)') ' vilma_restart_read: the file carries ', nl_f, &
              ' memory channels, the model ', self%resp%nlam, '.'
         if (nl_f > self%resp%nlam) write(*,'(a)') '   The file holds toroidal memory from a run with '// &
              'laterally varying viscosity; this run has none to put it in.'
         error stop 'vilma_restart_read: memory channel count does not match the model'
      end if
      nk_f      = nc_size(filename, "nk")
      cross_res = (nk_f /= nk)
      if (cross_res .and. nk_f > nk) &
         error stop 'vilma_restart_read: cross-resolution downsampling (file lmax > model) not supported'
      if (cross_res) then
         ! degree-grouped contiguity holds iff nk_f is the model's cumulative slot
         ! count through some degree L (i.e. same mmax/mres convention).
         ok_block = .false.
         do L = 1, self%resp%lmax
            if (self%resp%kbeg(L+1) - 1 == nk_f) then;  ok_block = .true.;  exit;  end if
         end do
         if (.not. ok_block) &
            error stop 'vilma_restart_read: cross-resolution restart incompatible (mmax/mres mismatch)'
      end if

      if (cross_res) then
         ! prognostic Maxwell memory: copy the shared low-degree block, zero-pad above.
         call get3d(filename, "tau_a_re", self%resp%Are, ne, nl_f, nk_f, n)
         call get3d(filename, "tau_a_im", self%resp%Aim, ne, nl_f, nk_f, n)
         call get3d(filename, "tau_b_re", self%resp%Bre, ne, nl_f, nk_f, n)
         call get3d(filename, "tau_b_im", self%resp%Bim, ne, nl_f, nk_f, n)
         call get3d(filename, "tau_c_re", self%resp%Cre, ne, nl_f, nk_f, n)
         call get3d(filename, "tau_c_im", self%resp%Cim, ne, nl_f, nk_f, n)
         ! σ_n is on the file grid and the trapezoidal load is cheap to re-derive; let
         ! the first step re-prime it (O(Δt); exact for the explicit fe scheme).
         self%resp%sigma_primed = .false.
      else
         call check_reference(self, filename, np, nl)
         ! prognostic Maxwell memory at slice n
         call get3d(filename, "tau_a_re", self%resp%Are, ne, nl_f, nk, n)
         call get3d(filename, "tau_a_im", self%resp%Aim, ne, nl_f, nk, n)
         call get3d(filename, "tau_b_re", self%resp%Bre, ne, nl_f, nk, n)
         call get3d(filename, "tau_b_im", self%resp%Bim, ne, nl_f, nk, n)
         call get3d(filename, "tau_c_re", self%resp%Cre, ne, nl_f, nk, n)
         call get3d(filename, "tau_c_im", self%resp%Cim, ne, nl_f, nk, n)
         call read_diagnostics(self, filename, n, np, nl)
         ! restore the trapezoidal start-of-step load σ_n (the other prognostic piece
         ! of the implicit scheme); without it the first step would re-derive σ_n to
         ! only the SLE solver tolerance, breaking bit-for-bit continuation.
         call restore_sigma(self, filename, n)
      end if
   end subroutine read_ve_state

   subroutine check_reference(self, filename, np, nl)
      !! Verify the file's static reference fields match the initialised model, so
      !! a restored memory state is not paired with a different reference.
      type(solid_earth), intent(in) :: self
      character(len=*),   intent(in) :: filename
      integer,            intent(in) :: np, nl
      real(wp), allocatable :: ref(:,:)
      real(wp), parameter :: tol = 1.0e-3_wp
      allocate(ref(np,nl))
      call nc_read(filename, "z_bed_eq", ref)
      if (maxval(abs(ref - self%gg%z_bed_eq)) > tol) &
         error stop 'vilma_restart_read: z_bed_eq does not match the initialised model'
      call nc_read(filename, "h_ice_eq", ref)
      if (maxval(abs(ref - self%gg%h_ice_eq)) > tol) &
         error stop 'vilma_restart_read: h_ice_eq does not match the initialised model'
   end subroutine check_reference

   subroutine read_diagnostics(self, filename, n, np, nl)
      !! Restore the (lon,lat) diagnostic state so the restarted object reports
      !! correctly: ice thickness, RSL, bedrock and the ocean function at slice n.
      type(solid_earth), intent(inout) :: self
      character(len=*),   intent(in)    :: filename
      integer,            intent(in)    :: n, np, nl
      call get2d(filename, "h_ice",   self%gg%h_ice, np, nl, n)
      call get2d(filename, "rsl",     self%gg%rsl,   np, nl, n)
      call get2d(filename, "z_bed",   self%gg%z_bed, np, nl, n)
      call get2d(filename, "C_ocean", self%gg%C,     np, nl, n)
   end subroutine read_diagnostics

   subroutine restore_sigma(self, filename, n)
      !! Restore σ_n at time slice n into the response, marking it primed, if the
      !! snapshot recorded a tracked σ_n (sigma_primed flag set).
      type(solid_earth), intent(inout) :: self
      character(len=*),   intent(in)    :: filename
      integer,            intent(in)    :: n
      real(wp), allocatable :: sre(:), sim(:)
      real(wp) :: primed
      integer  :: nlm
      call nc_read(filename, "sigma_primed", primed, start=[n], count=[1])
      if (primed <= 0.5_wp) return
      nlm = self%resp%nlm
      allocate(sre(nlm), sim(nlm))
      call nc_read(filename, "sigma_n_re", sre, start=[1,n], count=[nlm,1])
      call nc_read(filename, "sigma_n_im", sim, start=[1,n], count=[nlm,1])
      call response_prime_sigma(self%resp, cmplx(sre, sim, wp))
   end subroutine restore_sigma

   subroutine get3d(filename, name, dat, ne, nl_f, nk_f, n)
      !! Read a Maxwell-memory field at slice n into dat(nlam, ne, nk). The file
      !! holds nl_f ≤ nlam channels and nk_f ≤ nk degree-grouped slots; whatever it
      !! does not hold — the higher degrees of a cross-resolution restart, the
      !! toroidal channels of a spheroidal file — is zero. Read directly when the
      !! shapes match, so a same-run restart continues bit for bit.
      character(len=*), intent(in)    :: filename, name
      real(wp),         intent(inout) :: dat(:,:,:)
      integer,          intent(in)    :: ne, nl_f, nk_f, n
      real(wp), allocatable :: buf(:,:,:)
      if (nl_f == size(dat,1) .and. nk_f == size(dat,3)) then
         call nc_read(filename, name, dat, start=[1,1,1,n], count=[nl_f, ne, nk_f, 1])
         return
      end if
      allocate(buf(nl_f, ne, nk_f))
      call nc_read(filename, name, buf, start=[1,1,1,n], count=[nl_f, ne, nk_f, 1])
      dat = 0.0_wp
      dat(1:nl_f,:,1:nk_f) = buf
   end subroutine get3d

   subroutine get2d(filename, name, dat, np, nl, n)
      character(len=*), intent(in)  :: filename, name
      real(wp),         intent(out) :: dat(:,:)
      integer,          intent(in)  :: np, nl, n
      call nc_read(filename, name, dat, start=[1,1,n], count=[np, nl, 1])
   end subroutine get2d


end module vilma_io
