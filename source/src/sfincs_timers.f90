module sfincs_timers
   !
   ! Named wall-clock timers for SFINCS.
   !
   ! Lightweight replacement for the scattered tloop*/tstart_*/tend_*
   ! bookkeeping that used to live in each module. Timers are registered
   ! lazily: the first timer_start('name') with a new name creates it;
   ! subsequent calls find the existing record and accumulate.
   !
   ! All timing is done via system_clock with integer*8 counts, so 64-bit
   ! counters (typical on modern systems) do not wrap within any realistic
   ! SFINCS run.
   !
   ! Thread safety: timer_start / timer_stop are intended to be called
   ! from the serial driver, outside of !$omp parallel regions. They are
   ! NOT thread-safe.
   !
   use sfincs_log
   !
   implicit none
   !
   private
   !
   public :: timer_start
   public :: timer_stop
   public :: timer_reset
   public :: timer_elapsed
   public :: timer_count
   public :: timer_total_wall
   public :: timer_write_summary
   public :: timer_write_headers
   public :: timer_write_runtimes_file
   !
   integer, parameter :: name_len = 32
   integer, parameter :: max_timers = 64
   !
   type :: timer_record
      character(len=name_len) :: name         = ''
      real(8)                 :: accumulated  = 0.0_8
      integer(8)              :: last_start   = 0_8
      integer                 :: n_calls      = 0
      logical                 :: running      = .false.
      logical                 :: warned_start = .false.
      logical                 :: warned_stop  = .false.
   end type timer_record
   !
   type(timer_record), save :: timers(max_timers)
   integer,            save :: n_timers = 0
   logical,            save :: warned_full = .false.
   !
contains
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   integer function timer_find(name) result(idx)
   !
   ! Return the index of the timer with the given name, or 0 if not present.
   !
   character(len=*), intent(in) :: name
   integer                      :: i
   !
   idx = 0
   !
   do i = 1, n_timers
      !
      if (trim(timers(i)%name) == trim(name)) then
         idx = i
         return
      endif
      !
   enddo
   !
   end function timer_find
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   integer function timer_find_or_register(name) result(idx)
   !
   ! Return the index of the timer with the given name, creating a new
   ! record if it did not yet exist. Returns 0 if the table is full.
   !
   character(len=*), intent(in) :: name
   !
   idx = timer_find(name)
   !
   if (idx > 0) return
   !
   if (n_timers >= max_timers) then
      !
      if (.not. warned_full) then
         !
         call write_log(' Warning: sfincs_timers table full, timer ignored: '//trim(name), 1)
         warned_full = .true.
         !
      endif
      !
      idx = 0
      return
      !
   endif
   !
   n_timers = n_timers + 1
   idx      = n_timers
   !
   timers(idx)%name         = name
   timers(idx)%accumulated  = 0.0_8
   timers(idx)%last_start   = 0_8
   timers(idx)%n_calls      = 0
   timers(idx)%running      = .false.
   timers(idx)%warned_start = .false.
   timers(idx)%warned_stop  = .false.
   !
   end function timer_find_or_register
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine timer_start(name)
   !
   ! Start (or resume-and-accumulate-on-stop) the timer with the given name.
   ! Lazily registers a new timer on first call.
   !
   character(len=*), intent(in) :: name
   integer                      :: idx
   integer(8)                   :: c, rate
   !
   idx = timer_find_or_register(name)
   !
   if (idx == 0) return
   !
   if (timers(idx)%running) then
      !
      if (.not. timers(idx)%warned_start) then
         !
         call write_log(' Warning: timer_start on already-running timer: '//trim(name), 1)
         timers(idx)%warned_start = .true.
         !
      endif
      !
      return
      !
   endif
   !
   call system_clock(c, rate)
   !
   timers(idx)%last_start = c
   timers(idx)%running    = .true.
   !
   end subroutine timer_start
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine timer_stop(name)
   !
   ! Stop the timer and add the elapsed interval to its accumulated total.
   !
   character(len=*), intent(in) :: name
   integer                      :: idx
   integer(8)                   :: c, rate
   !
   idx = timer_find(name)
   !
   if (idx == 0) then
      !
      call write_log(' Warning: timer_stop on unknown timer: '//trim(name), 1)
      return
      !
   endif
   !
   if (.not. timers(idx)%running) then
      !
      if (.not. timers(idx)%warned_stop) then
         !
         call write_log(' Warning: timer_stop on non-running timer: '//trim(name), 1)
         timers(idx)%warned_stop = .true.
         !
      endif
      !
      return
      !
   endif
   !
   call system_clock(c, rate)
   !
   timers(idx)%accumulated = timers(idx)%accumulated + real(c - timers(idx)%last_start, 8) / real(rate, 8)
   timers(idx)%n_calls     = timers(idx)%n_calls + 1
   timers(idx)%running     = .false.
   !
   end subroutine timer_stop
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine timer_reset(name)
   !
   ! Reset a single timer's accumulated time and call count to zero.
   !
   character(len=*), intent(in) :: name
   integer                      :: idx
   !
   idx = timer_find(name)
   !
   if (idx == 0) return
   !
   timers(idx)%accumulated = 0.0_8
   timers(idx)%n_calls     = 0
   timers(idx)%running     = .false.
   !
   end subroutine timer_reset
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   real(8) function timer_elapsed(name) result(elapsed)
   !
   ! Accumulated wall time (in seconds) for the named timer.
   ! Returns 0 if the timer is unknown. If the timer is currently running,
   ! the interval since the most recent timer_start is included (without
   ! modifying the stored accumulated value).
   !
   character(len=*), intent(in) :: name
   integer                      :: idx
   integer(8)                   :: c, rate
   !
   idx = timer_find(name)
   !
   if (idx == 0) then
      elapsed = 0.0_8
      return
   endif
   !
   elapsed = timers(idx)%accumulated
   !
   if (timers(idx)%running) then
      !
      call system_clock(c, rate)
      !
      elapsed = elapsed + real(c - timers(idx)%last_start, 8) / real(rate, 8)
      !
   endif
   !
   end function timer_elapsed
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   integer function timer_count(name) result(n)
   !
   ! Number of completed start/stop cycles for the named timer.
   !
   character(len=*), intent(in) :: name
   integer                      :: idx
   !
   idx = timer_find(name)
   !
   if (idx == 0) then
      n = 0
      return
   endif
   !
   n = timers(idx)%n_calls
   !
   end function timer_count
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   real(8) function timer_total_wall() result(total)
   !
   ! Sum of accumulated wall time across all registered timers.
   !
   integer :: i
   !
   total = 0.0_8
   !
   do i = 1, n_timers
      total = total + timers(i)%accumulated
   enddo
   !
   end function timer_total_wall
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine timer_write_summary(to_screen, total_wall, min_elapsed)
   !
   ! Pretty-print a summary of all registered timers via write_log.
   !
   ! to_screen   : passed to write_log (1 = also echo to stdout).
   ! total_wall  : reference total wall time used for the '%' column.
   !               If <= 0, the sum across all timers is used instead.
   ! min_elapsed : timers with accumulated time below this threshold (in s)
   !               are skipped. Pass a negative value to print every timer.
   !
   integer, intent(in) :: to_screen
   real(8), intent(in) :: total_wall
   real(8), intent(in) :: min_elapsed
   !
   real(8)        :: denom
   real(8)        :: t_el
   real(8)        :: pct
   integer        :: i
   character(32)  :: call_label
   character(256) :: line
   !
   if (total_wall > 0.0_8) then
      denom = total_wall
   else
      denom = max(timer_total_wall(), 1.0e-12_8)
   endif
   !
   do i = 1, n_timers
      !
      t_el = timers(i)%accumulated
      !
      if (t_el < min_elapsed) cycle
      !
      pct = 100.0_8 * t_el / denom
      !
      if (timers(i)%n_calls == 1) then
         write(call_label, '(i0,a)') timers(i)%n_calls, ' call'
      else
         write(call_label, '(i0,a)') timers(i)%n_calls, ' calls'
      endif
      !
      write(line, '(1x,a,t25,a,f10.3,a,f5.1,a,a,a)') &
         trim(timers(i)%name), ': ', t_el, ' (', pct, '%, ', trim(call_label), ')'
      !
      call write_log(trim(line), to_screen)
      !
   enddo
   !
   end subroutine timer_write_summary
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine timer_write_headers(to_screen)
   !
   ! Write the three 'Total time / Total simulation time / Time in input' header
   ! lines to the log, using the 'Input' and 'Simulation loop' named timers.
   !
   integer, intent(in) :: to_screen
   !
   real(8) :: t_input
   real(8) :: t_loop
   !
   t_input = timer_elapsed('Input')
   t_loop  = timer_elapsed('Simulation loop')
   !
   write(logstr, '(a,f10.3)') ' Total time             : ', t_input + t_loop
   call write_log(trim(logstr), to_screen)
   !
   write(logstr, '(a,f10.3)') ' Total simulation time  : ', t_loop
   call write_log(trim(logstr), to_screen)
   !
   write(logstr, '(a,f10.3)') ' Time in input          : ', t_input
   call write_log(trim(logstr), to_screen)
   !
   end subroutine timer_write_headers
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine timer_write_runtimes_file(unit, filename)
   !
   ! Write the runtimes.txt payload: simulation-loop wall time, input wall time,
   ! and each phase timer, in the same order and with the same keys as the
   ! previous inline implementation in sfincs_lib.f90.
   !
   integer,          intent(in) :: unit
   character(len=*), intent(in) :: filename
   !
   open(unit, file=filename)
   !
   write(unit, '(f10.3,a)') real(timer_elapsed('Simulation loop'),     4), ' % total'
   write(unit, '(f10.3,a)') real(timer_elapsed('Input'),               4), ' % input'
   write(unit, '(f10.3,a)') real(timer_elapsed('Boundaries'),          4), ' % boundaries'
   write(unit, '(f10.3,a)') real(timer_elapsed('Discharges'),          4), ' % discharges'
   write(unit, '(f10.3,a)') real(timer_elapsed('Drainage structures'), 4), ' % drainage_structures'
   write(unit, '(f10.3,a)') real(timer_elapsed('Meteo fields'),        4), ' % meteo1'
   write(unit, '(f10.3,a)') real(timer_elapsed('Meteo forcing'),       4), ' % meteo2'
   write(unit, '(f10.3,a)') real(timer_elapsed('Infiltration'),        4), ' % infiltration'
   write(unit, '(f10.3,a)') real(timer_elapsed('Momentum'),            4), ' % momentum'
   write(unit, '(f10.3,a)') real(timer_elapsed('Structures'),          4), ' % structures'
   write(unit, '(f10.3,a)') real(timer_elapsed('Continuity'),          4), ' % continuity'
   write(unit, '(f10.3,a)') real(timer_elapsed('Output'),              4), ' % output'
   !
   close(unit)
   !
   end subroutine timer_write_runtimes_file
   !
end module sfincs_timers
