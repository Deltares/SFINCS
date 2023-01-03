module sfincs_crosssections

contains

   subroutine read_crs_file()
   !
   ! Reads crs file
   !
   use sfincs_data
   use geometry
   use quadtree
   !
   implicit none
   !
   integer ip, nm, nmu, ncrs, nrows, ncols, irow, stat, icrs, nr_points, iuv, indx, nr
   !
   real      :: dummy
   real*4    :: phiuv
   real*4    :: phic
   real*4    :: dphi
   character(len=256) :: cdummy
   !
   real*4, dimension(:),   allocatable :: xcrs
   real*4, dimension(:),   allocatable :: ycrs
   real*4,  dimension(2)                :: xp
   real*4,  dimension(2)                :: yp
   integer, dimension(:),   allocatable :: uv_indices   
   integer, dimension(:),   allocatable :: vertices
   !
   ! Read cross sections file
   !
   ncrs = 0
   nrcrosssections = 0
   !
   if (crsfile(1:4) /= 'none') then
      ! 
      write(*,*)'Reading cross sections ...'
      !
      ! First count number of polylines
      !
      open(500, file=trim(crsfile))
      do while(.true.)
         read(500,*,iostat = stat)cdummy
         read(500,*,iostat = stat)nrows,ncols
         if (stat<0) exit
         ncrs = ncrs + 1
         do irow = 1, nrows
            read(500,*)dummy
         enddo
      enddo
      rewind(500)
      !
      nrcrosssections = ncrs
      !
      allocate(crs_uv_index(1000, ncrs))
      allocate(crs_idir(1000,ncrs))
      allocate(crs_nr(ncrs))
      allocate(namecrs(ncrs))
      !
      ! Loop through polylines
      !
      do icrs = 1, ncrs
         !
         crs_nr(icrs) = 0
         nr = 0
         !
         read(500,*,iostat = stat)cdummy
         read(500,*,iostat = stat)nrows,ncols
         if (stat<0) exit
         allocate(xcrs(nrows))
         allocate(ycrs(nrows))
         do irow = 1, nrows
            read(500,*)xcrs(irow),ycrs(irow)
         enddo
         !
         namecrs(icrs) = trim(cdummy)
         !
         call find_uv_points_intersected_by_polyline(uv_indices, vertices, nr_points, xcrs, ycrs, nrows)
         !
         ! Now loop through uv points
         !
         do iuv = 1, nr_points
            !
            indx = uv_indices(iuv)
            irow = vertices(iuv)
            !
            nr                     = nr + 1
            crs_nr(icrs)           = nr
            crs_uv_index(nr, icrs) = indx
            !
            ! Check the angle at which these points cross
            !
            nm  = uv_index_z_nm(indx)
            nmu = uv_index_z_nmu(indx)
            !
            phiuv = atan2(z_yz(nmu) - z_yz(nm), z_xz(nmu) - z_xz(nm))
            phic  = atan2(ycrs(irow + 1) - ycrs(irow), xcrs(irow + 1) - xcrs(irow))
            dphi  = phiuv - phic
            if (dphi<0.0)  dphi = dphi + 2*pi
            if (dphi>2*pi) dphi = dphi - 2*pi
            !
            if (dphi<=pi) then
               crs_idir(nr, icrs) = 1
            else   
               crs_idir(nr, icrs) = -1
            endif   
            !
         enddo
         !
         deallocate(xcrs)
         deallocate(ycrs)
         !
      enddo
      !
      close(500)
      !
   endif
   !
   end subroutine

   
   subroutine get_discharges_through_crosssections(qq)
   !
   use sfincs_data
   !
   implicit none
   !
   real*4, dimension(:), allocatable :: qq
   !
   integer :: ip, icrs, iref, iuv, indx
   real*4  :: dxycrs
   !
   allocate(qq(nrcrosssections))
   !
   qq = 0.0
   !
   do icrs = 1, nrcrosssections
      !
      do ip = 1, crs_nr(icrs)
         !
         indx  = crs_uv_index(ip, icrs)
         iref  = uv_flags_iref(indx)
         iuv   = uv_flags_dir(indx) ! 0 is u, 1 is v
         !
         if (iuv==0) then
            !
            ! U
            !
            dxycrs = dyrm(iref) 
            !
         else   
            !
            ! V
            !
            if (crsgeo) then
               dxycrs = dxm(indx) 
            else
               dxycrs = dxrm(iref) 
            endif   
            !
         endif   
         !
         qq(icrs) = qq(icrs) + q(indx)*crs_idir(ip, icrs)*dxycrs
         !
      enddo
   enddo   
   !
   end subroutine
   
end module
