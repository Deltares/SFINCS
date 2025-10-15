module sfincs_runup_gauges

contains

   subroutine read_rug_file()
   !
   ! Reads rug file
   !
   use sfincs_data
   use geometry
   use quadtree
   use sfincs_error
   use sfincs_log
   !
   implicit none
   !
   integer ip, nrows, ncols, irow, stat, irug, nmq, nrpmx
   !
   real      :: dummy
   real*4    :: xru0, yru0, xru1, yru1, x, y
   real*4    :: rlen, rdx, rdy, dxstep
   logical ok
   !
   character(len=256) :: cdummy
   !
   ! Read cross sections file
   !
   nr_runup_gauges = 0
   !
   if (rugfile(1:4) /= 'none') then
      ! 
      call write_log('Info    : reading run-up gauges', 0)
      !
      ok = check_file_exists(obsfile, 'Run-up gauges rug file', .true.)
      !      
      ! First count number of polylines
      !
      open(500, file=trim(rugfile))
      do while(.true.)
         read(500,*,iostat = stat)cdummy
         if (stat<0) exit         
         read(500,*,iostat = stat)nrows, ncols
         if (stat<0) exit
         if (nrows > 2) then
            call write_log('Warning : Invalid run-up gauge! May only contain begin and end point. Run-up gauge skipped.', 1)
         endif   
         nr_runup_gauges = nr_runup_gauges + 1
         do irow = 1, nrows
            read(500,*)dummy
         enddo
      enddo
      rewind(500)
      !
      if (nr_runup_gauges == 0) then
         !
         call write_log('Warning : no valid run-up gauges found', 1)
         return         
         !
      endif
      !
      ! Get nr points in longest runup gauge
      !
      nrpmx = 0
      !
      dxstep = 0.2 * dxyr(nref)
      !
      do irug = 1, nr_runup_gauges
         !
         read(500,*,iostat = stat)cdummy
         read(500,*,iostat = stat)nrows,ncols
         if (stat<0) exit
         read(500,*)xru0, yru0
         read(500,*)xru1, yru1
         rdx = xru1 - xru0
         rdy = yru1 - yru0
         rlen = sqrt(rdx**2 + rdy**2)
         nrpmx = max(nrpmx, int(rlen / dxstep) + 1)
         !
      enddo
      rewind(500)
      !
      allocate(runup_gauge_nm(nrpmx, nr_runup_gauges))
      allocate(runup_gauge_name(nr_runup_gauges))
      allocate(runup_gauge_nrp(nr_runup_gauges))
      !
      runup_gauge_nm = 0
      runup_gauge_nrp = 0
      !
      ! Loop through polylines
      !
      do irug = 1, nr_runup_gauges
         !
         read(500,*,iostat = stat)cdummy
         read(500,*,iostat = stat)nrows,ncols
         if (stat<0) exit
         !
         ! We only take the first two vertices
         !
         read(500,*)xru0, yru0
         read(500,*)xru1, yru1
         !
         runup_gauge_name(irug) = trim(cdummy)
         !
         rdx = xru1 - xru0
         rdy = yru1 - yru0
         rlen = sqrt(rdx**2 + rdy**2)
         runup_gauge_nrp(irug) = int(rlen / dxstep) + 1
         !
         do ip = 1, runup_gauge_nrp(irug)
            !
            x = xru0 + (ip - 1) * dxstep * rdx / rlen
            y = yru0 + (ip - 1) * dxstep * rdy / rlen
            !
            ! Index in quadtree
            !
            nmq = find_quadtree_cell(x, y)
            !
            if (nmq > 0) then
               !
               ! Index in SFINCS domain
               !
               runup_gauge_nm(ip, irug) = index_sfincs_in_quadtree(nmq)
               !
            endif   
            !
         enddo   
         !
      enddo 
      !
      close(500)
      !
   endif
   !
   end subroutine

   
   subroutine get_runup_levels(zru)
   !
   use sfincs_data
   !
   implicit none
   !
   real*4, dimension(:), allocatable :: zru
   !
   integer :: ip, irug, nm
   real*4  :: zbt
   !
   allocate(zru(nr_runup_gauges))
   !
   zru = -999.0
   !
   do irug = 1, nr_runup_gauges
      !
      do ip = 1, runup_gauge_nrp(irug)
         !
         nm = runup_gauge_nm(ip, irug)
         !
         if (subgrid) then
            !
            zbt = subgrid_z_zmin(nm)
            !
         else
            !
            zbt = zb(nm)
            !
         endif
         !
         if (zs(nm) > zbt + runup_gauge_depth) then
            !
            ! Point is okay
            !
            zru(irug) = zs(nm)
            !
         else
            !
            ! Shallower than runup_gauge_depth
            !
            !exit
            !
         endif
         !
      enddo
   enddo   
   !
   end subroutine
   
end module
