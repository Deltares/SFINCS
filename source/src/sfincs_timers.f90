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
   ! This module is a pure data module and deliberately has NO dependency
   ! on sfincs_log: rendering/pretty-printing of timer data lives in
   ! sfincs_log (write_timer_headers_log, write_timer_summary_log,
   ! write_runtimes_file), which walks the timer list via the iteration
   ! API below. Keeping the two modules separated avoids a circular
   ! dependency (sfincs_log already uses timer_elapsed internally).
   !
   ! Subroutines / functions:
   !
   !   timer_start(name) / timer_stop(name)
   !     Start / stop a named timer. Lazily registers on first start.
   !     Called from every phase in sfincs_lib (Input, Simulation loop,
   !     Boundaries, Momentum, Continuity, Output, ...) and from
   !     update_wave_field in sfincs_snapwave.
   !
   !   timer_reset(name)
   !     Zero a single timer. Currently unused by the main driver but
   !     kept as part of the public API.
   !
   !   timer_elapsed(name) / timer_count(name)
   !     Read accumulated wall time / call count for a named timer.
   !     Called from sfincs_log (write_progress_log,
   !     write_finished_log, write_timer_headers_log,
   !     write_runtimes_file).
   !
   !   timer_total_wall()
   !     Sum of accumulated wall time across all registered timers.
   !     Called from write_timer_summary_log in sfincs_log.
   !
   !   timer_num_registered()
   !     Number of timers currently registered. Called from
   !     write_timer_summary_log in sfincs_log.
   !
   !   timer_name_by_index(i) / timer_elapsed_by_index(i) /
   !   timer_count_by_index(i)
   !     Iteration API: read a timer's stored name / accumulated wall
   !     time / call count by index. Indices run 1 .. timer_num_registered().
   !     Called from write_timer_summary_log in sfincs_log.
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
   public :: timer_num_registered
   public :: timer_name_by_index
   public :: timer_elapsed_by_index
   public :: timer_count_by_index
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
      ! Called from: timer_find_or_register, timer_stop, timer_reset,
      ! timer_elapsed, timer_count (all within this module).
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
      ! Called from: timer_start (within this module).
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
            write(*, '(a)') ' Warning: sfincs_timers table full, timer ignored: '//trim(name)
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
      ! Called from: sfincs_lib (main driver, every phase) and
      ! update_wave_field in sfincs_snapwave.
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
            write(*, '(a)') ' Warning: timer_start on already-running timer: '//trim(name)
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
      ! Called from: sfincs_lib (main driver, every phase) and
      ! update_wave_field in sfincs_snapwave.
      !
      character(len=*), intent(in) :: name
      integer                      :: idx
      integer(8)                   :: c, rate
      !
      idx = timer_find(name)
      !
      if (idx == 0) then
         !
         write(*, '(a)') ' Warning: timer_stop on unknown timer: '//trim(name)
         return
         !
      endif
      !
      if (.not. timers(idx)%running) then
         !
         if (.not. timers(idx)%warned_stop) then
            !
            write(*, '(a)') ' Warning: timer_stop on non-running timer: '//trim(name)
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
      ! Called from: (currently no live callers; part of the public API.)
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
      ! Called from: sfincs_log (write_progress_log, write_finished_log,
      ! write_timer_headers_log, write_runtimes_file).
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
      ! Called from: (currently no live callers; part of the public API.)
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
      ! Called from: write_timer_summary_log in sfincs_log.
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
   integer function timer_num_registered() result(n)
      !
      ! Number of timers currently registered. Used by the rendering
      ! routines in sfincs_log to iterate over every timer without
      ! touching module-private state.
      !
      ! Called from: write_timer_summary_log in sfincs_log.
      !
      n = n_timers
      !
   end function timer_num_registered
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   function timer_name_by_index(i) result(name)
      !
      ! Return the stored name of the i-th registered timer, or an empty
      ! string for out-of-range i. Indices run 1 .. timer_num_registered().
      !
      ! Called from: write_timer_summary_log in sfincs_log.
      !
      integer, intent(in)     :: i
      character(len=name_len) :: name
      !
      if (i < 1 .or. i > n_timers) then
         name = ''
         return
      endif
      !
      name = timers(i)%name
      !
   end function timer_name_by_index
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   real(8) function timer_elapsed_by_index(i) result(elapsed)
      !
      ! Accumulated wall time of the i-th registered timer. Returns 0 for
      ! out-of-range i. Does NOT include a running-interval contribution
      ! (use timer_elapsed(name) if you need that — the rendering code
      ! runs after the simulation loop has stopped all timers).
      !
      ! Called from: write_timer_summary_log in sfincs_log.
      !
      integer, intent(in) :: i
      !
      if (i < 1 .or. i > n_timers) then
         elapsed = 0.0_8
         return
      endif
      !
      elapsed = timers(i)%accumulated
      !
   end function timer_elapsed_by_index
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   integer function timer_count_by_index(i) result(n)
      !
      ! Number of completed start/stop cycles of the i-th registered
      ! timer. Returns 0 for out-of-range i.
      !
      ! Called from: write_timer_summary_log in sfincs_log.
      !
      integer, intent(in) :: i
      !
      if (i < 1 .or. i > n_timers) then
         n = 0
         return
      endif
      !
      n = timers(i)%n_calls
      !
   end function timer_count_by_index
   !
end module sfincs_timers
