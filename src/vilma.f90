module vilma
   !! VILMA2 umbrella module — single `use vilma` entry point that
   !! re-exports the public API of every component (including the parameter record
   !! vilma_param_class / vilma_par_load and the standalone driver vilma_run). Host
   !! models and tests should depend on this rather than the individual vilma_* modules.
   use vilma_precision,       only: wp, sp, dp
   use vilma_constants
   use vilma_params,          only: vilma_param_class, vilma_par_load, vilma_par_print
   use vilma_control,         only: vilma_ctl_class, vilma_ctl_load, vilma_ctl_print, DEFAULTS_FILE
   use vilma_sht,             only: sht_grid
   use vilma_earth_structure, only: earth_model, build_earth, build_M3L70V01
   use vilma_radial_fe,       only: radial_operator
   use vilma_viscoelastic,    only: ve_degree
   use vilma_sle,             only: sle_solver
   use vilma_rotation,        only: rotation_state
   use vilma_coupling,        only: solid_earth
   use vilma_io,              only: vilma_restart_write, vilma_restart_read, vilma_write_step, &
                                 vilma_io_set_table
   use vilma_drive,           only: vilma_run
   implicit none
   public

#ifndef VERSION
#define VERSION "unknown"
#endif
   character(len=*), parameter :: vilma_version = VERSION

end module vilma
