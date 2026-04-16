module sfincs_discharges
   !
   ! River point discharges: nsrc (x,y) locations from srcfile with matching
   ! time series qsrc_ts(:,:) from disfile, OR from a FEWS-style netCDF input
   ! via netsrcdisfile. Interpolates to the current model time every step,
   ! stores the interpolated value in qtsrc(nsrc) (for his output), and
   ! accumulates the per-cell discharge into the global qsrc(np) array used
   ! by sfincs_continuity.
   !
   ! Drainage structures (pumps, check valves, culverts, controlled gates)
   ! live in sfincs_src_structures. The two modules no longer share any
   ! arrays -- they cooperate only by both writing into qsrc(np).
   !
   use sfincs_log
   use sfincs_error
   !
contains
   !
   subroutine initialize_discharges()
   !
   ! Read src/dis or netsrcdis. Allocate nmindsrc(nsrc), qtsrc(nsrc).
   !
   use sfincs_data
   use sfincs_ncinput
   use quadtree
   !
   implicit none
   !
   real*4    :: dummy
   integer   :: isrc, itsrc, nmq, n, stat
   logical   :: ok
   !
   nsrc      = 0
   ntsrc     = 0
   itsrclast = 1
   !
   if (srcfile(1:4) /= 'none') then
      !
      ok = check_file_exists(srcfile, 'Source points file', .true.)
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
         nsrc = nsrc + 1
         !
      enddo
      !
      rewind(500)
      !
   elseif (netsrcdisfile(1:4) /= 'none') then    ! FEWS-compatible NetCDF discharge time series
      !
      ok = check_file_exists(netsrcdisfile, 'Netcdf river input netsrcdis file', .true.)
      !
      call read_netcdf_discharge_data()   ! sets nsrc, ntsrc, xsrc, ysrc, qsrc_ts, tsrc
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
   if (nsrc <= 0) return
   !
   allocate(nmindsrc(nsrc))
   allocate(qtsrc(nsrc))
   !
   nmindsrc = 0
   qtsrc    = 0.0
   !
   ! --- Read src/dis contents for the srcfile case ---------------------
   !
   if (srcfile(1:4) /= 'none') then
      !
      allocate(xsrc(nsrc))
      allocate(ysrc(nsrc))
      !
      do n = 1, nsrc
         !
         read(500, *) xsrc(n), ysrc(n)
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
      allocate(qsrc_ts(nsrc, ntsrc))
      !
      do itsrc = 1, ntsrc
         !
         read(502, *) tsrc(itsrc), (qsrc_ts(isrc, itsrc), isrc = 1, nsrc)
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
   ! --- Map river sources to grid cells --------------------------------
   !
   do isrc = 1, nsrc
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
   !
   subroutine update_discharges(t, dt, tloop)
   !
   ! Zero qsrc(np); interpolate the river discharge time series to t,
   ! store in qtsrc(1..nsrc), and accumulate into qsrc(nmindsrc(:)).
   !
   ! update_discharges is called BEFORE update_src_structures -- that is
   ! why it owns the qsrc zeroing. Both routines then additively write
   ! their contributions.
   !
   use sfincs_data
   !
   implicit none
   !
   real*8  :: t
   real*4  :: dt
   real    :: tloop
   !
   integer :: count0, count1, count_rate, count_max
   integer :: isrc, itsrc, nm, it_prev, it_next
   real*4  :: wt
   !
   call system_clock(count0, count_rate, count_max)
   !
   ! Zero qsrc for this step. sfincs_src_structures will add to it next.
   !
   !$acc kernels present( qsrc )
   qsrc = 0.0
   !$acc end kernels
   !
   if (nsrc > 0) then
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
      wt = (t - tsrc(it_prev)) / (tsrc(it_next) - tsrc(it_prev))
      !
      ! Atomic accumulation because two river sources (or a river and a
      ! structure) can share a cell.
      !
      !$acc parallel loop present( qsrc, qtsrc, nmindsrc, qsrc_ts ) private( nm )
      !$omp parallel do private( nm ) schedule ( static )
      do isrc = 1, nsrc
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
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0 * (count1 - count0) / count_rate
   !
   end subroutine

end module
