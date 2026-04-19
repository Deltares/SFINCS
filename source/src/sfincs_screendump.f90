module sfincs_screendump
   !
   ! User-facing screen / log output for SFINCS.
   !
   ! Formatted blocks that live here:
   !   - screendump_startup      : welcome banner + ASCII art + build info
   !   - screendump_processes    : yes/no "Processes" summary
   !   - screendump_progress     : per-timestep progress / ETA line
   !   - screendump_run_finished : end-of-run banner + timer summary +
   !                               average time step
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
   public :: screendump_progress
   public :: screendump_run_finished
   !
   ! Next percentage threshold at which the progress reporter prints a
   ! line. Incremented in steps of percdoneval (set from the
   ! 'percentage_done' input keyword). Zero-initialised so the first
   ! call prints at 0%.
   !
   real, save :: percdonenext = 0.0
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
   subroutine screendump_progress(t, t0, t1)
   !
   ! Per-timestep progress reporter. Prints a "NN% complete, TT.T s
   ! remaining ..." line each time the simulated-time percentage
   ! crosses the next percdoneval threshold. Remaining time is
   ! estimated from the wall-clock elapsed in the 'Simulation loop'
   ! timer.
   !
   use sfincs_timers, only: timer_elapsed
   !
   implicit none
   !
   real*8, intent(in) :: t
   real*4, intent(in) :: t0, t1
   !
   real                :: percdone, trun, trem
   character(len=256)  :: logstr
   !
   percdone = min(100 * (t - t0) / (t1 - t0), 100.0)
   !
   if (percdone >= percdonenext) then
      !
      ! percdoneval is increment of % to show to log, default=+5%
      !
      percdonenext = 1.0 * (int(percdone) + percdoneval)
      !
      trun = real(timer_elapsed('Simulation loop'), 4)
      trem = trun / max(0.01*percdone, 1.0e-6) - trun
      !
      if (int(percdone)>0) then
         !
         write(logstr,'(i4,a,f7.1,a)')int(percdone),'% complete, ',trem,' s remaining ...'
         call write_log(logstr, 1)
         !
      else
         !
         write(logstr,'(i4,a,f7.1,a)')int(percdone),'% complete,       - s remaining ...'
         call write_log(logstr, 1)
         !
      endif
      !
   endif
   !
   end subroutine screendump_progress
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine screendump_run_finished(dtavg)
   !
   ! End-of-run log block: "Simulation finished" banner, per-phase
   ! timer summary, and the average time step line. Called once from
   ! sfincs_finalize, after the simulation loop has stopped and
   ! dtavg has been averaged.
   !
   use sfincs_timers, only: timer_write_headers, timer_write_summary, timer_elapsed
   !
   implicit none
   !
   real, intent(in) :: dtavg
   !
   character(len=256) :: logstr
   !
   call write_log('', 1)
   call write_log('---------- Simulation finished -----------', 1)
   call write_log('', 1)
   !
   call timer_write_headers(1)
   !
   ! Per-phase timing summary. Percentages are relative to the total
   ! wall time of the simulation loop.
   !
   call timer_write_summary(1, timer_elapsed('Simulation loop'), 0.0005_8)
   !
   call write_log('', 1)
   !
   write(logstr,'(a,20f10.3)')           ' Average time step (s)  : ', dtavg
   call write_log(logstr, 1)
   !
   call write_log('', 1)
   !
   end subroutine screendump_run_finished
   !
   !-----------------------------------------------------------------------------------------------------!
   !
end module sfincs_screendump
