module sfincs_obspoints

contains
   ! 
   subroutine read_obs_points()
   !
   ! Reads obs files
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   real*4 xtmp, ytmp, dummy, di1, dj1
   !
   integer iobs, nm, m, n, nmq, stat, j1, j2, jdq, iref
   !
   character(len=256)        :: line
   character(len=256)        :: line2
   !
   real*4, dimension(:), allocatable :: value
   !
   ! Read observation points
   !
   nobs = 0
   !
   if (obsfile(1:4) /= 'none') then
      ! 
      write(*,*)'Reading observation points ...'
      !
      open(500, file=trim(obsfile))       
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         nobs = nobs + 1
      enddo
      rewind(500)
      allocate(xobs(nobs))
      allocate(yobs(nobs))
      allocate(zobs(nobs))
      allocate(hobs(nobs))
      allocate(nmindobs(nobs))
      allocate(mindobs(nobs))
      allocate(nindobs(nobs)) 
      allocate(idobs(nobs))  
      allocate(nameobs(nobs))     
      allocate(xgobs(nobs))
      allocate(ygobs(nobs))
      allocate(zbobs(nobs))      
      !
      allocate(nmwindobs(nobs))
      allocate(wobs(4, nobs))
      !
      allocate(value(2))
      !
      value(1) = 0.0
      value(2) = 0.0     
      !
      do n = 1, nobs
         !
         read(500,'(a)')line                           
         j1=index(line,"'")
         jdq=index(line,'"')
         if (j1 == 0 .and. jdq==0) then! no name supplied, give standard name
            j2 = 11
            nameobs(n) = ''
            write(nameobs(n)(1:j2), '(A8,I0.3)') 'station_', n       
         elseif (j1>0) then ! name supplied,         
            line2 = adjustl(trim(line(j1+1:256)))
            j2=index(line2,"'")      
            nameobs(n) = adjustl(trim(line2(1:j2-1)))
         else
            line2 = adjustl(trim(line(jdq+1:256)))
            j2=index(line2,'"')      
            nameobs(n) = adjustl(trim(line2(1:j2-1)))            
         endif 
         !
         read(line,*)(value(m), m = 1, 2)         
         xobs(n) = value(1)
         yobs(n) = value(2)
         ! 
      enddo
      !
      close(500)
      !
      ! Determine m and n indices of observation points
      !
      do iobs = 1, nobs 
         !
         nmindobs(iobs) = 0
         mindobs(iobs)  = 0
         nindobs(iobs)  = 0
         idobs(iobs)    = iobs                   
         xgobs(iobs)    = -999.0
         ygobs(iobs)    = -999.0
         zbobs(iobs)    = -999.0 
         !
         nmq = find_quadtree_cell(xobs(iobs), yobs(iobs))
         !
         if (nmq>0) then
            !
            nm             = index_sfincs_in_quadtree(nmq)
            nmindobs(iobs) = nm
            n              = z_index_z_n(nm)
            m              = z_index_z_m(nm)
            !
            xgobs(iobs)    = z_xz(nm)   
            ygobs(iobs)    = z_yz(nm)
            !
            if (subgrid) then
               zbobs(iobs) = subgrid_z_zmin(nm)
            else
               zbobs(iobs) = zb(nm)
            endif
            !
            iref = z_flags_iref(nm)
            !
!            write(*,'(a,i0,a,a,a,i0,a,i0,a,i0,a,i0,a,f0.3)')' Observation point ',iobs,' : "',trim(nameobs(iobs)),'" nm=',nm,' n=',n,' m=',m,' iref=',iref,' z=',zbobs(iobs)
            !
         else
            !
            write(*,'(a,a,a)')' Warning : observation point "', trim(nameobs(iobs)), '" falls outside model domain.'
            !
         endif   
         !
      enddo   
      !
   endif
   !
   end subroutine
   !
end module
