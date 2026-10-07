module vilma_control
   !! Standalone-program control record: the run-management settings the VILMA2
   !! *executables* need but the solid-Earth model itself does not. Loaded from one
   !! namelist group `&ctl`, separate from the physics/numerics record `&vilma`
   !! (vilma_params). These are the forcing/reference/output file paths, variable names,
   !! the time window, the online-remap toggle, the reference-equilibration selector
   !! (i_eq), and the restart-in path — everything a host model (CLIMBER-X) supplies
   !! through the API and its own time loop instead, and so never reads here.
   !!
   !! Split rationale: `vilma_param_class` is the model's configuration contract (what a
   !! host fills in memory); `vilma_ctl_class` is the standalone driver's I/O glue. The
   !! driver (vilma_drive) and the offline tools (vilma_remap, vilma_mkref) load
   !! BOTH groups — `&vilma` for the grid/physics, `&ctl` for what to read and write.
   use vilma_precision, only: wp
   use vilma_constants, only: sec_per_year
   use vilma_params,    only: expand_path
   use nml
   implicit none
   private

   public :: vilma_ctl_class, vilma_ctl_load, vilma_ctl_print
   public :: DEFAULTS_FILE

   !! Canonical physics defaults the executables load automatically. A run config
   !! given on the command line is overlaid on this (yelmo defaults_file convention).
   !! Relative to the run directory; input/ is linked into each rundir (see .runme).
   character(len=*), parameter :: DEFAULTS_FILE = "input/vilma_defaults.nml"

   type :: vilma_ctl_class
      ! --- ice-thickness forcing (lon,lat,time) ---------------------------------
      character(len=512) :: file_forcing = ""            !! ice-thickness forcing
      character(len=64)  :: name_ice     = "h_ice"       !! ice variable in file_forcing
      character(len=64)  :: name_time    = "time"        !! time axis [years] in file_forcing

      ! --- legacy Gauss-grid reference file (remap_input=.false.) ----------------
      character(len=512) :: file_ref     = ""            !! reference state (z_bed_eq, h_ice_eq)
      character(len=64)  :: name_zbed_eq = "z_bed_eq"    !! bed var (legacy 2D ref; or 3D bed in forcing)
      character(len=64)  :: name_hice_ref = "h_ice_ref"  !! ice var (legacy 2D ref)

      ! --- step output ----------------------------------------------------------
      character(len=512) :: file_out     = "vilma_out.nc"  !! step output
      ! Surface horizontal displacement, spheroidal + toroidal, written to its own
      ! file on the same grid and time axis as file_out; empty = not written. Its
      ! own file because each field is as large as rsl, and most runs need none.
      character(len=512) :: file_hor     = ""                  !! horizontal output (off when empty)

      ! --- forcing time window --------------------------------------------------
      ! Time fields are SI [s] in the record; the nml supplies them in YEARS and
      ! vilma_ctl_load converts on read (so the in-memory record is uniformly SI).
      real(wp) :: time_init = -huge(1.0_wp)  !! start time [years] (default: first forcing slice)
      real(wp) :: time_end  =  huge(1.0_wp)  !! end time   [years] (default: last  forcing slice)

      ! --- online lon-lat -> Gauss remap of the forcing -------------------------
      logical  :: remap_input = .true.   !! .true.: forcing is lon-lat, remap on the fly;
                                         !! .false.: forcing already on the Gauss grid (legacy)
      character(len=64) :: name_lon = "lon"  !! source longitude axis variable [deg]
      character(len=64) :: name_lat = "lat"  !! source latitude  axis variable [deg]

      ! --- reference / equilibration selector -----------------------------------
      ! The relaxed reference (z_bed_eq = SLE topo0, h_ice_eq) is set by i_eq,
      ! mirroring CLIMBER-X i_equilibrium (src/geo/geo.f90). RSL is measured against
      ! this reference, so i_eq=1 (present-day reference) yields rsl ~0 at the present
      ! day. Reference fields are lon-lat and remapped online.
      integer  :: i_eq = 1
         !! 0: start slice is the equilibrium (z_bed_eq=bed[k0], h_ice_eq=ice[k0]);
         !! 1: present-day reference from z_bed_ref_file / h_ice_ref_file (default);
         !! 2: equilibrium read from z_bed_eq_file / h_ice_eq_file;
         !! 3: z_bed_eq = present-day reference bed + rsl from rsl_restart_file.
      character(len=512) :: z_bed_ref_file = ""   !! present-day reference bed  (i_eq=1,3)
      character(len=512) :: h_ice_ref_file = ""   !! present-day reference ice  (i_eq=1,3)
      character(len=512) :: z_bed_eq_file  = ""   !! equilibrium bed (i_eq=2)
      character(len=512) :: h_ice_eq_file  = ""   !! equilibrium ice (i_eq=2)
      character(len=512) :: rsl_restart_file = "" !! rsl field added to ref bed (i_eq=3)
      character(len=64)  :: name_z_bed_ref = "bedrock_topography" !! bed var in the ref/eq files
      character(len=64)  :: name_h_ice_ref = "ice_thickness"      !! ice var in the ref/eq files
      character(len=64)  :: name_rsl       = "rsl"                !! rsl var in rsl_restart_file

      ! --- restart-in path ------------------------------------------------------
      ! The restart *capability* (vilma_restart_read/write) is the model's; the *path* is
      ! run management. A host (CLIMBER-X) derives its own restart path from its global
      ! restart_in_dir and never reads this field.
      character(len=512) :: restart_in_file = ""
         !! full-state restart (vilma_restart.nc) to resume from: solid_earth_init at the
         !! reference, then restore the saved memory/clock. A lower-resolution restart is
         !! interpolated up to the model grid.
   end type vilma_ctl_class

contains

   subroutine vilma_ctl_load(c, filename, group)
      !! Fill the control record from the `&ctl` group of `filename`. There is no
      !! separate &ctl defaults file (the canonical defaults — input/vilma_defaults.nml
      !! — carry only &vilma, the host API contract). Reads are therefore non-strict:
      !! a parameter present in `filename` overrides, one that is absent keeps the
      !! vilma_ctl_class in-code default. So a run config may set only the &ctl keys it
      !! needs. Override `group` to read a differently-named namelist. Time fields are
      !! given in YEARS and converted to SI here.
      type(vilma_ctl_class), intent(inout) :: c
      character(len=*),   intent(in)    :: filename
      character(len=*),   intent(in), optional :: group
      character(len=64)  :: g
      real(wp) :: time_init_yr, time_end_yr

      g = "ctl";  if (present(group)) g = group
      call nml_set_verbose(.false.)             ! vilma_ctl_print echoes a concise summary instead

      ! forcing
      call nml_read(filename, g, "file_forcing",  c%file_forcing)
      c%file_forcing = expand_path(c%file_forcing)
      call nml_read(filename, g, "name_ice",      c%name_ice)
      call nml_read(filename, g, "name_time",     c%name_time)

      ! legacy Gauss-grid reference file
      call nml_read(filename, g, "file_ref",      c%file_ref)
      c%file_ref = expand_path(c%file_ref)
      call nml_read(filename, g, "name_zbed_eq",  c%name_zbed_eq)
      call nml_read(filename, g, "name_hice_ref", c%name_hice_ref)

      ! output
      call nml_read(filename, g, "file_out",      c%file_out)
      c%file_out = expand_path(c%file_out)
      call nml_read(filename, g, "file_hor",      c%file_hor)
      if (len_trim(c%file_hor) > 0) c%file_hor = expand_path(c%file_hor)

      ! time window (YEARS in the nml -> SI seconds in the record)
      time_init_yr = c%time_init/sec_per_year
      time_end_yr  = c%time_end /sec_per_year
      call nml_read(filename, g, "time_init",     time_init_yr)
      call nml_read(filename, g, "time_end",      time_end_yr)
      c%time_init = time_init_yr*sec_per_year
      c%time_end  = time_end_yr *sec_per_year

      ! online remap
      call nml_read(filename, g, "remap_input",   c%remap_input)
      call nml_read(filename, g, "name_lon",      c%name_lon)
      call nml_read(filename, g, "name_lat",      c%name_lat)

      ! reference / equilibration selector
      call nml_read(filename, g, "i_eq",            c%i_eq)
      call nml_read(filename, g, "z_bed_ref_file",  c%z_bed_ref_file)
      c%z_bed_ref_file = expand_path(c%z_bed_ref_file)
      call nml_read(filename, g, "h_ice_ref_file",  c%h_ice_ref_file)
      c%h_ice_ref_file = expand_path(c%h_ice_ref_file)
      call nml_read(filename, g, "z_bed_eq_file",   c%z_bed_eq_file)
      c%z_bed_eq_file = expand_path(c%z_bed_eq_file)
      call nml_read(filename, g, "h_ice_eq_file",   c%h_ice_eq_file)
      c%h_ice_eq_file = expand_path(c%h_ice_eq_file)
      call nml_read(filename, g, "rsl_restart_file", c%rsl_restart_file)
      c%rsl_restart_file = expand_path(c%rsl_restart_file)
      call nml_read(filename, g, "name_z_bed_ref",  c%name_z_bed_ref)
      call nml_read(filename, g, "name_h_ice_ref",  c%name_h_ice_ref)
      call nml_read(filename, g, "name_rsl",        c%name_rsl)

      ! restart-in path
      call nml_read(filename, g, "restart_in_file", c%restart_in_file)
      c%restart_in_file = expand_path(c%restart_in_file)
   end subroutine vilma_ctl_load

   subroutine vilma_ctl_print(c, unit)
      !! Echo the active control configuration (to stdout, or `unit` if given).
      type(vilma_ctl_class), intent(in) :: c
      integer, intent(in), optional :: unit
      integer :: u
      u = 6;  if (present(unit)) u = unit

      write(u,'(a)')       ' [ctl] standalone-run configuration'
      write(u,'(a,a)')     '   forcing:  ', trim(c%file_forcing)
      write(u,'(a,l1,a,i0)') '   remap_input=', c%remap_input, '   i_eq=', c%i_eq
      write(u,'(a,es12.4,a,es12.4,a)') '   window:   t=[', c%time_init/sec_per_year, &
           ',', c%time_end/sec_per_year, '] yr'
      write(u,'(a,a)')     '   output:   ', trim(c%file_out)
      if (len_trim(c%file_hor) > 0) &
         write(u,'(a,a)')  '   horizontal: ', trim(c%file_hor)
      if (len_trim(c%restart_in_file) > 0) &
         write(u,'(a,a)')  '   restart:  ', trim(c%restart_in_file)
   end subroutine vilma_ctl_print

end module vilma_control
