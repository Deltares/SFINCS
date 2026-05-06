module sfincs_log
   !
   ! User-facing log / screen output for SFINCS.
   !
   ! Owns sfincs.log (fid, open_log/close_log/write_log) and all the
   ! formatted blocks that the driver writes to it.
   !
   ! Rendering of named-timer data (headers, per-phase summary, the
   ! runtimes.txt payload) also lives here, so that sfincs_timers can
   ! remain a pure data module with no dependency on sfincs_log. This
   ! breaks what used to be a circular dependency between the two.
   !
   ! Subroutines:
   !
   !   open_log() / close_log() / write_log(str, to_screen)
   !     File handle management and the single-line writer. Called from
   !     every SFINCS module that emits user-facing output.
   !
   !   write_startup_log()
   !     Welcome banner + ASCII logo + build-revision / build-date lines.
   !     Called once from sfincs_initialize (sfincs_lib).
   !
   !   write_processes_log()
   !     "Processes" yes/no summary. Called once from sfincs_initialize.
   !
   !   write_progress_log(t, t0, t1)
   !     Per-timestep progress / ETA line. Called every time step from
   !     the main loop in sfincs_lib. Uses timer_elapsed('simulation').
   !
   !   write_finished_log(dtavg)
   !     End-of-run banner + per-phase timer summary + average time step.
   !     Called once from sfincs_finalize (sfincs_lib).
   !
   !   write_timer_headers_log(to_screen)
   !     Three-line "Total / Total simulation / Input" header block.
   !     Called from write_finished_log.
   !
   !   write_timer_summary_log(to_screen, total_wall, min_elapsed)
   !     Per-timer summary table (name, seconds, % of total, #calls).
   !     Walks the timer list via the iteration API on sfincs_timers.
   !     Called from write_finished_log.
   !
   !   write_runtimes_file(unit, filename)
   !     Writes the runtimes.txt payload (simulation-loop wall time,
   !     input wall time, and each phase timer) in the format the
   !     original inline code in sfincs_lib produced. Called from
   !     sfincs_finalize (sfincs_lib) when write_time_output is set.
   !
   !   fmt_real(val, decimals) result(s)
   !     Format a real value with the minimum necessary field width and
   !     a guaranteed leading zero for |val| < 1. Works around a quirk
   !     in ifx that drops the leading zero for the "f0.d" edit
   !     descriptor. Returns a 32-char string, left-justified; callers
   !     use trim(fmt_real(...)) when embedding it in a larger format.
   !
   use sfincs_timers
   !
   integer        :: fid
   character(256) :: logstr
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
   subroutine open_log()
      !
      ! Open sfincs.log on the module-local unit fid=777. Called once
      ! at the very start of sfincs_initialize (sfincs_lib).
      !
      implicit none
      !
      fid = 777
      open(unit = fid, file = 'sfincs.log')
      !
   end subroutine open_log
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine write_log(str, to_screen)
      !
      ! Write one line to sfincs.log, optionally echoed to stdout.
      ! Called from every SFINCS module that emits user-facing output.
      !
      implicit none
      !
      character(*), intent(in) :: str
      integer, intent(in)      :: to_screen
      !
      write(fid,'(a)') trim(str)
      !
      if (to_screen == 1) then
         write(*,'(a)') trim(str)
      endif
      !
   end subroutine write_log
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   function fmt_real(val, decimals) result(s)
      !
      ! Format a real with minimum width and a guaranteed leading zero
      ! for |val| < 1. ifx's "f0.d" descriptor drops the leading zero in
      ! that range, which is not standard-conforming; this helper rewrites
      ! the result so the log output always reads "0.6670" rather than
      ! ".6670".
      !
      ! Called from: write_src_structures_log_summary (sfincs_src_structures),
      ! urban_drainage log summary (sfincs_urban_drainage), and anywhere
      ! else a real needs to be embedded in a log line with the smallest
      ! reasonable field width.
      !
      implicit none
      !
      real,    intent(in) :: val
      integer, intent(in) :: decimals
      character(len=32)   :: s
      !
      character(len=16)   :: fmt
      !
      write(fmt,'(a,i0,a)') '(f0.', decimals, ')'
      write(s,fmt) val
      s = adjustl(s)
      !
      if (s(1:1) == '.') then
         !
         s = '0' // s(1:len_trim(s))
         !
      else if (s(1:2) == '-.') then
         !
         s = '-0' // trim(s(2:))
         !
      endif
      !
   end function fmt_real
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine close_log()
      !
      ! Close the sfincs.log file handle. Called once at the end of
      ! sfincs_finalize (sfincs_lib).
      !
      implicit none
      !
      close(fid)
      !
   end subroutine close_log
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine write_startup_log()
      !
      ! Welcome banner, ASCII logo and build-revision / build-date lines.
      ! Called once at the start of sfincs_initialize (sfincs_lib),
      ! after build_revision and build_date have been set in sfincs_data.
      !
      use sfincs_data
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
   end subroutine write_startup_log
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine write_processes_log()
      !
      ! "Processes" summary block listing which physical processes are
      ! enabled for this run. Reads the process flags from sfincs_data.
      ! Called once from sfincs_initialize (sfincs_lib).
      !
      use sfincs_data
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
      if (discharges) then
         call write_log('Discharges           : yes', 1)
      else
         call write_log('Discharges           : no', 1)
      endif
      !
      if (drainage_structures) then
         call write_log('Drainage structures  : yes', 1)
      else
         call write_log('Drainage structures  : no', 1)
      endif
      !
      if (urban_drainage) then
         call write_log('Urban drainage       : yes', 1)
      else
         call write_log('Urban drainage       : no', 1)
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
   end subroutine write_processes_log
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine write_progress_log(t, t0, t1)
      !
      ! Per-timestep progress reporter. Prints a "NN% complete, TT.T s
      ! remaining ..." line each time the simulated-time percentage
      ! crosses the next percdoneval threshold. Remaining time is
      ! estimated from the wall-clock elapsed in the 'simulation'
      ! timer.
      !
      ! Called every time step from the main loop in sfincs_lib.
      !
      use sfincs_data, only: percdoneval
      !
      implicit none
      !
      real*8, intent(in) :: t
      real*4, intent(in) :: t0, t1
      !
      real                :: percdone, trun, trem
      character(len=256)  :: logstr
      !
      percdone = min(100.0 * (real(t, 4) - t0) / (t1 - t0), 100.0)
      !
      if (percdone >= percdonenext) then
         !
         ! percdoneval is increment of % to show to log, default=+5%
         !
         percdonenext = 1.0 * (int(percdone) + percdoneval)
         !
         trun = real(timer_elapsed('simulation'), 4)
         trem = trun / max(0.01*percdone, 1.0e-6) - trun
         !
         if (int(percdone) > 0) then
            !
            write(logstr,'(i4,a,f7.1,a)') int(percdone),'% complete, ',trem,' s remaining ...'
            call write_log(logstr, 1)
            !
         else
            !
            write(logstr,'(i4,a,f7.1,a)') int(percdone),'% complete,       - s remaining ...'
            call write_log(logstr, 1)
            !
         endif
         !
      endif
      !
   end subroutine write_progress_log
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine write_finished_log(dtavg)
      !
      ! End-of-run log block: "Simulation finished" banner, per-phase
      ! timer summary, and the average time step line. Called once from
      ! sfincs_finalize (sfincs_lib), after the simulation loop has
      ! stopped and dtavg has been averaged.
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
      call write_timer_headers_log(1)
      !
      ! Per-phase timing summary. Percentages are relative to the total
      ! wall time of the simulation loop.
      !
      call write_timer_summary_log(1, timer_elapsed('simulation'), 0.0005_8)
      !
      call write_log('', 1)
      !
      write(logstr,'(a,20f10.3)') ' Average time step (s)  : ', dtavg
      call write_log(logstr, 1)
      !
      call write_log('', 1)
      !
   end subroutine write_finished_log
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine write_timer_headers_log(to_screen)
      !
      ! Write the three 'Total time / Total simulation time / Time in input' header
      ! lines to the log, using the 'input' and 'simulation' named timers.
      !
      ! Called from: write_finished_log.
      !
      integer, intent(in) :: to_screen
      !
      real(8) :: t_input
      real(8) :: t_loop
      !
      t_input = timer_elapsed('input')
      t_loop  = timer_elapsed('simulation')
      !
      write(logstr, '(a,f10.3)') ' Total time                   : ', t_input + t_loop
      call write_log(trim(logstr), to_screen)
      !
!      write(logstr, '(a,f10.3)') ' Total simulation time        : ', t_loop
!      call write_log(trim(logstr), to_screen)
      !
      write(logstr, '(a,f10.3)') ' Time in input                : ', t_input
      call write_log(trim(logstr), to_screen)
      !
   end subroutine write_timer_headers_log
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine write_timer_summary_log(to_screen, total_wall, min_elapsed)
      !
      ! Pretty-print a summary of all registered timers via write_log.
      ! Walks the timer list via the iteration API on sfincs_timers
      ! (timer_num_registered, timer_name_by_index, timer_elapsed_by_index,
      ! timer_count_by_index) so this module does not need to know about
      ! the internal storage of sfincs_timers.
      !
      ! to_screen   : passed to write_log (1 = also echo to stdout).
      ! total_wall  : reference total wall time used for the '%' column.
      !               If <= 0, the sum across all timers is used instead.
      ! min_elapsed : timers with accumulated time below this threshold (in s)
      !               are skipped. Pass a negative value to print every timer.
      !
      ! Called from: write_finished_log.
      !
      integer, intent(in) :: to_screen
      real(8), intent(in) :: total_wall
      real(8), intent(in) :: min_elapsed
      !
      real(8)        :: denom
      real(8)        :: t_el
      real(8)        :: pct
      integer        :: i
      integer        :: n
      integer        :: ncalls
      character(32)  :: tname
      character(256) :: line
      !
      if (total_wall > 0.0_8) then
         denom = total_wall
      else
         denom = max(timer_total_wall(), 1.0e-12_8)
      endif
      !
      n = timer_num_registered()
      !
      do i = 1, n
         !
         t_el = timer_elapsed_by_index(i)
         !
         if (t_el < min_elapsed) cycle
         !
         ! Skip input (was already added in header)
         !
         if (trim(timer_name_by_index(i)) == 'input') cycle
         !
         pct    = 100.0_8 * t_el / denom
         tname  = timer_name_by_index(i)
         !
         write(line, '(1x,a,1x,a,t31,a,f10.3,a,f5.1,a,a,a)') &
            'Time in', trim(tname), ': ', t_el, ' (', pct, '%)'
         !
         call write_log(trim(line), to_screen)
         !
      enddo
      !
   end subroutine write_timer_summary_log
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine write_runtimes_file(unit, filename)
      !
      ! Write the runtimes.txt payload: simulation-loop wall time, input wall time,
      ! and each phase timer, in the same order and with the same keys as the
      ! previous inline implementation in sfincs_lib.f90.
      !
      ! Called from: sfincs_finalize (sfincs_lib) when write_time_output
      ! is set.
      !
      integer,          intent(in) :: unit
      character(len=*), intent(in) :: filename
      !
      open(unit, file=filename)
      !
      write(unit, '(f10.3,a)') real(timer_elapsed('simulation'),     4), ' % total'
      write(unit, '(f10.3,a)') real(timer_elapsed('input'),               4), ' % input'
      write(unit, '(f10.3,a)') real(timer_elapsed('boundaries'),          4), ' % boundaries'
      write(unit, '(f10.3,a)') real(timer_elapsed('discharges'),          4), ' % discharges'
      write(unit, '(f10.3,a)') real(timer_elapsed('drainage structures'), 4), ' % drainage_structures'
      write(unit, '(f10.3,a)') real(timer_elapsed('urban drainage'),     4), ' % urban_drainage'
      write(unit, '(f10.3,a)') real(timer_elapsed('meteo fields'),        4), ' % meteo1'
      write(unit, '(f10.3,a)') real(timer_elapsed('meteo forcing'),       4), ' % meteo2'
      write(unit, '(f10.3,a)') real(timer_elapsed('infiltration'),        4), ' % infiltration'
      write(unit, '(f10.3,a)') real(timer_elapsed('momentum'),            4), ' % momentum'
      write(unit, '(f10.3,a)') real(timer_elapsed('structures'),          4), ' % structures'
      write(unit, '(f10.3,a)') real(timer_elapsed('continuity'),          4), ' % continuity'
      write(unit, '(f10.3,a)') real(timer_elapsed('output'),              4), ' % output'
      !
      close(unit)
      !
   end subroutine write_runtimes_file
   !
end module sfincs_log
