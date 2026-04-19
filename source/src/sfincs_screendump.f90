module sfincs_screendump
   !
   ! User-facing screen / log output for SFINCS.
   !
   ! Two large formatted blocks live here:
   !   - screendump_startup    : welcome banner + ASCII art + build info
   !   - screendump_processes  : yes/no "Processes" summary
   !
   ! The per-timestep progress reporter stays inline in sfincs_lib.f90.
   !
   use sfincs_log
   use sfincs_data
   !
   implicit none
   !
   private
   !
   public :: screendump_startup
   public :: screendump_processes
   !
contains
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine screendump_startup()
   !
   ! Welcome banner, ASCII logo and build-revision / build-date lines.
   ! Called once at the start of sfincs_initialize, after build_revision
   ! and build_date have been set in sfincs_data.
   !
   implicit none
   !
   call write_log('', 1)
   call write_log('------------ Welcome to SFINCS ------------', 1)
   call write_log('', 1)
   call write_log('  @@@@@  @@@@@@@ @@ @@  @@   @@@@   @@@@@ ', 1)
   call write_log(' @@@ @@@ @@@@@@@ @@ @@@ @@ @@@@@@@ @@@ @@@', 1)
   call write_log(' @@@     @@      @@ @@@ @@ @@   @@ @@@    ', 1)
   call write_log('  @@@@@  @@@@@@  @@ @@@@@@ @@       @@@@@ ', 1)
   call write_log('     @@@ @@      @@ @@ @@@ @@   @@     @@@', 1)
   call write_log(' @@@ @@@ @@      @@ @@  @@  @@@@@@ @@@ @@@', 1)
   call write_log('  @@@@@  @@      @@ @@   @   @@@@   @@@@@ ', 1)
   call write_log('', 1)
   call write_log('              ..............              ', 1)
   call write_log('          ......:@@@@@@@@:......          ', 1)
   call write_log('       ..::::..@@........@@.:::::..       ', 1)
   call write_log('     ..:::::..@@..::..::..@@.::::::..     ', 1)
   call write_log('    .::::::..@@............@@.:::::::.    ', 1)
   call write_log('   .::::::..@@..............@@.:::::::.   ', 1)
   call write_log('  .::::::::..@@............@@..::::::::.  ', 1)
   call write_log(' .:::::::::...@@.@..@@..@.@@..::::::::::. ', 1)
   call write_log(' .:::::::::...:@@@..@@..@@@:..:::::::::.. ', 1)
   call write_log(' ............@@.@@..@@..@@.@@............ ', 1)
   call write_log(' ^^^~~^^~~^^@@..............@@^^^~^^^~~^^ ', 1)
   call write_log(' .::::::::::@@..............@@.:::::::::. ', 1)
   call write_log('  .......:.@@.....@.....@....@@.:.......  ', 1)
   call write_log('   .::....@@......@.@@@.@....@@.....::.   ', 1)
   call write_log('    .:::~@@.:...:.@@...@@.:.:.@@~::::.    ', 1)
   call write_log('     .::~@@@@@@@@@@.....@@@@@@@@@~::.     ', 1)
   call write_log('       ..:~~~~~~~:.......:~~~~~~~:..      ', 1)
   call write_log('          ......................          ', 1)
   call write_log('              ..............              ', 1)
   call write_log('', 1)
   call write_log('------------------------------------------', 1)
   call write_log('', 1)
   call write_log('Build-Revision: '//trim(build_revision), 1)
   call write_log('Build-Date: '//trim(build_date), 1)
   call write_log('', 1)
   !
   end subroutine screendump_startup
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine screendump_processes()
   !
   ! "Processes" summary block listing which physical processes are
   ! enabled for this run. Reads the process flags from sfincs_data.
   !
   implicit none
   !
   call write_log('', 1)
   call write_log('------------------------------------------', 1)
   call write_log('Processes', 1)
   call write_log('------------------------------------------', 1)
   !
   if (subgrid) then
      call write_log('Subgrid topography   : yes', 1)
   else
      call write_log('Subgrid topography   : no', 1)
   endif
   !
   if (use_quadtree) then
      call write_log('Quadtree refinement  : yes', 1)
   else
      call write_log('Quadtree refinement  : no', 1)
   endif
   !
   if (advection) then
      call write_log('Advection            : yes', 1)
   else
      call write_log('Advection            : no', 1)
   endif
   !
   if (viscosity) then
      call write_log('Viscosity            : yes', 1)
   else
      call write_log('Viscosity            : no', 1)
   endif
   !
   if (coriolis) then
      call write_log('Coriolis             : yes', 1)
   else
      call write_log('Coriolis             : no', 1)
   endif
   !
   if (wind) then
      call write_log('Wind                 : yes', 1)
   else
      call write_log('Wind                 : no', 1)
   endif
   !
   if (patmos) then
      call write_log('Atmospheric pressure : yes', 1)
   else
      call write_log('Atmospheric pressure : no', 1)
   endif
   !
   if (precip) then
      call write_log('Precipitation        : yes', 1)
   else
      call write_log('Precipitation        : no', 1)
   endif
   !
   if (infiltration) then
      call write_log('Infiltration         : yes', 1)
   else
      call write_log('Infiltration         : no', 1)
   endif
   !
   if (snapwave) then
      call write_log('SnapWave             : yes', 1)
   else
      call write_log('SnapWave             : no', 1)
   endif
   !
   if (wavemaker) then
      call write_log('Wave paddles         : yes', 1)
   else
      call write_log('Wave paddles         : no', 1)
   endif
   !
   if (nonhydrostatic) then
      call write_log('Non-hydrostatic      : yes', 1)
   else
      ! call write_log('Non-hydrostatic         : no', 1)
   endif
   !
   if (bathtub) then
      call write_log('Bathtub              : yes', 1)
   else
      ! call write_log('Bathtub              : no', 1)
   endif
   !
   call write_log('------------------------------------------', 1)
   call write_log('', 1)
   !
   end subroutine screendump_processes
   !
   !-----------------------------------------------------------------------------------------------------!
   !
end module sfincs_screendump
