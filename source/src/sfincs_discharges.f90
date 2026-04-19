module sfincs_discharges
   !
   ! River point discharges: nr_discharge_points (x,y) locations from srcfile with matching
   ! time series qsrc_ts(:,:) from disfile, OR from a FEWS-style netCDF input
   ! via netsrcdisfile. Interpolates to the current model time every step,
   ! stores the interpolated value in qtsrc(nr_discharge_points) (for his output), and
   ! accumulates the per-cell discharge into the global qsrc(np) array used
   ! by sfincs_continuity.
   !
   ! Drainage structures (pumps, check valves, culverts, controlled gates)
   ! live in sfincs_src_structures. The two modules no longer share any
   ! arrays -- they cooperate only by both writing into qsrc(np).
   !
   ! Subroutines:
   !
   !   initialize_discharges()
   !     Read srcfile/disfile (ascii) or netsrcdisfile (netcdf), resolve
   !     each source to its quadtree cell, and allocate runtime state.
   !     Called from sfincs_lib at init time.
   !
   !   update_discharges(t, dt)
   !     Zero qsrc(np), interpolate the river discharge time series to the
   !     current time, and accumulate into qsrc at each source cell.
   !     Called from update_continuity (sfincs_continuity) once per
   !     time step, before update_src_structures.
   !
   !   count_tokens(line, ntok)
   !     Count whitespace-separated tokens in a string; used to decide
   !     between the 2-column (x y) and 3-column (x y name) src formats.
   !     Called from initialize_discharges (this module).
   !
   use sfincs_log
   use sfincs_error
   !
   implicit none
   !
   ! Module-level runtime state for river point discharges (moved from
   ! sfincs_data). The count, coordinate arrays, file-path strings, and
   ! qsrc_ts / tsrc / ntsrc stay in sfincs_data because they are also
   ! set by sfincs_input (keyword reader) or sfincs_ncinput (which this
   ! module uses, so a back-reference would be circular).
   !
   ! Public so downstream output modules (sfincs_output, sfincs_ncoutput)
   ! and the openacc bookkeeping module can reference them.
   !
   ! Name length (matches src_struc_name_len from sfincs_src_structures).
   !
   integer, parameter, public :: src_name_len = 128
   !
   ! Per-river-source names
   !
   character(len=src_name_len), dimension(:), allocatable, public :: src_name   ! (nr_discharge_points) user-supplied or auto-generated names for river point sources
   !
   ! Runtime state
   !
   integer,                                   public :: itsrclast   ! last-used bracket index into tsrc, for time-series interpolation
   real*4,    dimension(:), allocatable,      public :: qtsrc       ! (nr_discharge_points) interpolated discharge at current time, for his output
   integer*4, dimension(:), allocatable,      public :: nmindsrc    ! (nr_discharge_points) river source cell indices
   !
contains
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine initialize_discharges()
      !
      ! Read the river-point-discharge input and wire each source up to a grid
      ! cell. Two mutually-exclusive input paths:
      !   - srcfile (+ disfile): ascii, 2-column (x y) or 3-column (x y name)
      !     location file plus a separate time-series file.
      !   - netsrcdisfile: FEWS-style netcdf with locations and time series
      !     in one file (no per-point names; auto-generated).
      ! Allocates nmindsrc(nr_discharge_points) and qtsrc(nr_discharge_points),
      ! and populates shared tsrc/qsrc_ts arrays in sfincs_data.
      !
      ! Called from: sfincs_lib (once, at init time).
      !
      use sfincs_data
      use sfincs_ncinput
      use quadtree
      !
      implicit none
      !
      real*4    :: dummy
      integer   :: isrc, itsrc, nmq, n, stat, ntok
      logical   :: ok
      character(len=1024) :: line, line_trim
      character(len=src_name_len) :: name_tmp
      !
      nr_discharge_points = 0
      ntsrc               = 0
      itsrclast           = 1
      !
      if (srcfile(1:4) /= 'none') then
         !
         write(logstr,'(a)') 'Info    : reading discharges'
         call write_log(logstr, 0)
         !
         ok = check_file_exists(srcfile, 'River input locations src file', .true.)
         !
         open(500, file=trim(srcfile))
         !
         do while (.true.)
            !
            read(500, *, iostat=stat) dummy
            if (stat < 0) exit
            nr_discharge_points = nr_discharge_points + 1
            !
         enddo
         !
         rewind(500)
         !
      elseif (netsrcdisfile(1:4) /= 'none') then    ! FEWS-compatible NetCDF discharge time series
         !
         ok = check_file_exists(netsrcdisfile, 'Netcdf river input netsrcdis file', .true.)
         !
         call read_netcdf_discharge_data()   ! sets nr_discharge_points, ntsrc, xsrc, ysrc, qsrc_ts, tsrc
         !
         ! The netcdf discharge file does not carry per-point names; auto-generate
         ! the same way as the 2-column srcfile path.
         !
         if (nr_discharge_points > 0) then
            !
            allocate(src_name(nr_discharge_points))
            !
            src_name = ' '
            !
            do n = 1, nr_discharge_points
               !
               write(src_name(n), '(a,i4.4)') 'discharge_', n
               !
            enddo
            !
         endif
         !
         if ((tsrc(1) > (t0 + 1.0)) .or. (tsrc(ntsrc) < (t1 - 1.0))) then
            !
            write(logstr,'(a)') ' WARNING! Times in discharge file do not cover entire simulation period!'
            call write_log(logstr, 1)
            !
         endif
         !
      endif
      !
      if (nr_discharge_points <= 0) return
      !
      allocate(nmindsrc(nr_discharge_points))
      allocate(qtsrc(nr_discharge_points))
      !
      nmindsrc = 0
      qtsrc    = 0.0
      !
      ! Read src/dis contents for the srcfile case
      !
      if (srcfile(1:4) /= 'none') then
         !
         allocate(xsrc(nr_discharge_points))
         allocate(ysrc(nr_discharge_points))
         allocate(src_name(nr_discharge_points))
         !
         src_name = ' '
         !
         do n = 1, nr_discharge_points
            !
            read(500, '(a)') line
            line_trim = adjustl(line)
            !
            ! Count whitespace-separated tokens on the line.
            !
            call count_tokens(line_trim, ntok)
            !
            if (ntok == 2) then
               !
               read(line_trim, *) xsrc(n), ysrc(n)
               write(src_name(n), '(a,i4.4)') 'discharge_', n
               !
            elseif (ntok == 3) then
               !
               read(line_trim, *) xsrc(n), ysrc(n), name_tmp
               src_name(n) = adjustl(trim(name_tmp))
               !
            else
               !
               write(logstr,'(a,i0,a,i0,a)') ' Error ! src file line ', n, ' has ', ntok, &
                  ' tokens -- expected 2 (x y) or 3 (x y name) !'
               call write_log(logstr, 1)
               error = 1
               return
               !
            endif
            !
         enddo
         !
         close(500)
         !
         ! Read discharge time series
         !
         ok = check_file_exists(disfile, 'River discharge timeseries dis file', .true.)
         !
         open(502, file=trim(disfile))
         !
         do while (.true.)
            !
            read(502, *, iostat=stat) dummy
            if (stat < 0) exit
            ntsrc = ntsrc + 1
            !
         enddo
         !
         rewind(502)
         allocate(tsrc(ntsrc))
         allocate(qsrc_ts(nr_discharge_points, ntsrc))
         !
         do itsrc = 1, ntsrc
            !
            read(502, *) tsrc(itsrc), (qsrc_ts(isrc, itsrc), isrc = 1, nr_discharge_points)
            !
         enddo
         !
         close(502)
         !
         if ((tsrc(1) > (t0 + 1.0)) .or. (tsrc(ntsrc) < (t1 - 1.0))) then
            !
            write(logstr,'(a)') 'Warning! Times in discharge file do not cover entire simulation period !'
            call write_log(logstr, 1)
            !
            if (tsrc(1) > (t0 + 1.0)) then
               !
               write(logstr,'(a)') 'Warning! Adjusting first time in discharge time series !'
               call write_log(logstr, 1)
               tsrc(1) = t0 - 1.0
               !
            else
               !
               write(logstr,'(a)') 'Warning! Adjusting last time in discharge time series !'
               call write_log(logstr, 1)
               tsrc(ntsrc) = t1 + 1.0
               !
            endif
            !
         endif
         !
      endif
      !
      ! Map river sources to grid cells
      !
      do isrc = 1, nr_discharge_points
         !
         nmq = find_quadtree_cell(xsrc(isrc), ysrc(isrc))
         !
         if (nmq > 0) then
            !
            nmindsrc(isrc) = index_sfincs_in_quadtree(nmq)
            !
         endif
         !
      enddo
      !
      deallocate(xsrc)
      deallocate(ysrc)
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine update_discharges(t, dt)
      !
      ! Zero qsrc(np); interpolate the river discharge time series to t,
      ! store in qtsrc(1..nr_discharge_points), and accumulate into qsrc(nmindsrc(:)).
      !
      ! update_discharges is called BEFORE update_src_structures -- that is
      ! why it owns the qsrc zeroing. Both routines then additively write
      ! their contributions.
      !
      ! Called from: update_continuity (sfincs_continuity), once per time step.
      !
      use sfincs_data
      use sfincs_timers
      !
      implicit none
      !
      real*8  :: t
      real*4  :: dt
      !
      integer :: isrc, itsrc, nm, it_prev, it_next
      real*4  :: wt
      !
      call timer_start('Discharges')
      !
      ! Zero qsrc for this step. sfincs_src_structures will add to it next.
      !
      !$acc kernels present( qsrc )
      qsrc = 0.0
      !$acc end kernels
      !
      if (nr_discharge_points > 0) then
         !
         ! Locate the bracketing interval in tsrc and compute the interpolation
         ! weight once. Then run a single parallel loop that both interpolates
         ! qtsrc and accumulates it into qsrc.
         !
         it_prev = itsrclast
         it_next = itsrclast + 1
         !
         do itsrc = itsrclast, ntsrc
            !
            if (tsrc(itsrc) > t) then
               !
               it_prev = itsrc - 1
               it_next = itsrc
               itsrclast = it_prev
               exit
               !
            endif
            !
         enddo
         !
         ! Clamp to valid bracket. If t is outside [tsrc(1), tsrc(ntsrc)] (which
         ! can happen on the netcdf path, where the srcfile pre-padding is not
         ! applied), hold the endpoint value rather than read out of bounds.
         !
         it_prev = min(max(it_prev, 1), ntsrc - 1)
         it_next = it_prev + 1
         !
         wt = (t - tsrc(it_prev)) / (tsrc(it_next) - tsrc(it_prev))
         !
         ! Atomic accumulation because two river sources (or a river and a
         ! structure) can share a cell.
         !
         !$acc parallel loop present( qsrc, qtsrc, nmindsrc, qsrc_ts ) private( nm )
         !$omp parallel do private( nm ) schedule ( static )
         do isrc = 1, nr_discharge_points
            !
            qtsrc(isrc) = qsrc_ts(isrc, it_prev) + (qsrc_ts(isrc, it_next) - qsrc_ts(isrc, it_prev)) * wt
            nm = nmindsrc(isrc)
            !
            if (nm > 0) then
               !
               !$acc atomic update
               !$omp atomic
               qsrc(nm) = qsrc(nm) + qtsrc(isrc)
               !
            endif
            !
         enddo
         !$omp end parallel do
         !$acc end parallel loop
         !
      endif
      !
      call timer_stop('Discharges')
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine count_tokens(line, ntok)
      !
      ! Count the number of whitespace-separated tokens in a string.
      ! Whitespace = spaces and tabs. Empty string returns 0.
      !
      ! Called from: initialize_discharges (this module) to disambiguate the
      ! 2-column vs 3-column srcfile layout.
      !
      implicit none
      !
      character(len=*), intent(in)  :: line
      integer,          intent(out) :: ntok
      !
      integer :: i, n
      logical :: in_tok
      character(len=1) :: c
      !
      ntok   = 0
      in_tok = .false.
      n      = len_trim(line)
      !
      do i = 1, n
         !
         c = line(i:i)
         !
         if (c == ' ' .or. c == char(9)) then
            !
            in_tok = .false.
            !
         else
            !
            if (.not. in_tok) then
               !
               ntok   = ntok + 1
               in_tok = .true.
               !
            endif
            !
         endif
         !
      enddo
      !
   end subroutine
   !
end module
