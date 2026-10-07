program vilma_main
   !! Standalone VILMA2 driver.
   !!
   !!   ./bin/vilma.x [run-config.nml]      (default: vilma.nml)
   !!
   !! The run config supplies the &ctl run-control group and any &vilma physics
   !! overrides; its &vilma is overlaid on the complete physics defaults in
   !! input/vilma_defaults.nml (loaded automatically — see DEFAULTS_FILE).
   use vilma_drive,   only: vilma_run
   use vilma_control, only: DEFAULTS_FILE
   implicit none
   character(len=512) :: cfg

   if (command_argument_count() >= 1) then
      call get_command_argument(1, cfg)
   else
      cfg = "vilma.nml"
   end if

   call vilma_run(trim(cfg), DEFAULTS_FILE)
end program vilma_main
