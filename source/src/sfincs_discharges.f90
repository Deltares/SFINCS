module sfincs_discharges

contains
   !
   subroutine read_discharges()
   !
   ! Reads discharge files
   !
   use sfincs_data
   use sfincs_ncinput
   !
   implicit none
   !
   real*4, dimension(:),     allocatable :: xsnk
   real*4, dimension(:),     allocatable :: ysnk
   !
   real*4 xtmp, ytmp, dummy
   !
   integer isrc, itsrc, idrn, nm, m, n, stat, j, iref
   !
   ! Read discharge points
   !
   nsrc  = 0
   ndrn  = 0
   ntsrc = 0
   itsrclast = 1
   !
   if (srcfile(1:4) /= 'none') then
      !
      write(*,*)'Reading discharges ...'
      open(500, file=trim(srcfile))
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         nsrc = nsrc + 1
      enddo
      rewind(500)
      !
   elseif (netsrcdisfile(1:4) /= 'none') then    ! FEWS compatible Netcdf discharge time-series input
      !
      call read_netcdf_discharge_data()  ! reads nsrc, ntsrc, xsrc, ysrc, qsrc, and tsrc
      !
   endif   
   !
   if (drnfile(1:4) /= 'none') then
      write(*,*)'Reading drainage file ...'
      open(501, file=trim(drnfile))
      do while(.true.)
         read(501,*,iostat = stat)dummy
         if (stat<0) exit
         ndrn = ndrn + 1
      enddo
      rewind(501)
   endif
   !
   nsrcdrn = nsrc + 2*ndrn
   !
   if (nsrcdrn>0) then
      allocate(nmindsrc(nsrcdrn))
      allocate(qtsrc(nsrcdrn))
      nmindsrc = 0
      qtsrc = 0.0
   endif
   !
   if (srcfile(1:4) /= 'none') then
      !
      ! Actually read src and dis files
      !
      allocate(xsrc(nsrc))
      allocate(ysrc(nsrc))
      !
      do n = 1, nsrc
         read(500,*)xsrc(n),ysrc(n)
      enddo
      close(500)
      !
      ! Read discharge time series
      !
      open(502, file=trim(disfile))
      do while(.true.)
         read(502,*,iostat = stat)dummy
         if (stat<0) exit
         ntsrc = ntsrc + 1
      enddo
      rewind(502)
      allocate(tsrc(ntsrc))
      allocate(qsrc(nsrc,ntsrc))
      do itsrc = 1, ntsrc
         read(502,*)tsrc(itsrc),(qsrc(isrc, itsrc), isrc = 1, nsrc)
      enddo
      close(502)
      !
   endif  
   !
   if (nsrc>0) then
      !
      ! Determine m and n indices of sources
      !
      do isrc = 1, nsrc
         !
         ! First determine temporary coordinates of source points in local coordinate system
         !
         xtmp =   cosrot*(xsrc(isrc) - x0) + sinrot*(ysrc(isrc) - y0)
         ytmp = - sinrot*(xsrc(isrc) - x0) + cosrot*(ysrc(isrc) - y0)
         !
         nmindsrc(isrc) = 0
         !
         ! Loop through refinement levels (from fine to coarse)
         !
         do iref = nref, 1, -1
            !
            m = int(xtmp/dxr(iref)) + 1
            n = int(ytmp/dyr(iref)) + 1
            !
            ! Find nm index of this src point
            !
            do nm = 1, np
               !
               if (z_flags_iref(nm) == iref) then
                  !
                  if (z_index_z_n(nm)==n .and. z_index_z_m(nm)==m) then
                     !
                     ! Point found
                     !
                     nmindsrc(isrc) = nm
                     !
                     exit
                     !
                  endif
               endif   
               !
               if (nmindsrc(isrc)>0) exit
               !
            enddo   
            !
         enddo   
         !
      enddo
      !
      ! Don't need coordinates anymore, and xsrc and ysrc may be used for drainage points as well
      !
      deallocate(xsrc)
      deallocate(ysrc)
      !
   endif   
   !
   ! And now the drainage points
   !
   if (ndrn>0) then
      !
      allocate(xsrc(ndrn))
      allocate(ysrc(ndrn))
      allocate(xsnk(ndrn))
      allocate(ysnk(ndrn))
      !
      allocate(drainage_type(ndrn))
      allocate(drainage_params(ndrn,5))
      !
      do idrn = 1, ndrn
         read(501,*)xsnk(idrn),ysnk(idrn),xsrc(idrn),ysrc(idrn),drainage_type(idrn),drainage_params(idrn,1),drainage_params(idrn,2),drainage_params(idrn,3),drainage_params(idrn,4),drainage_params(idrn,5)
      enddo
      close(501)
      !
      ! Determine m and n indices of source and sinks
      !
      do idrn = 1, ndrn
         !
         ! First determine temporary coordinates of sink points in local coordinate system
         !
         j = nsrc + idrn*2 - 1
         !
         xtmp =   cosrot*(xsnk(idrn) - x0) + sinrot*(ysnk(idrn) - y0)
         ytmp = - sinrot*(xsnk(idrn) - x0) + cosrot*(ysnk(idrn) - y0)
         !
         nmindsrc(j) = 0
         !
         ! Loop through refinement levels (from fine to coarse)
         !
         do iref = nref, 1, -1
            !
            n = int(ytmp/dyr(iref)) + 1
            m = int(xtmp/dxr(iref)) + 1
            !
            ! Find nm index of this src point
            !
            do nm = 1, np
               !
               if (z_flags_iref(nm) == iref) then
                  !
                  if (z_index_z_n(nm)==n .and. z_index_z_m(nm)==m) then
                     !
                     ! Point found
                     !
                     nmindsrc(j) = nm
                     !
                     exit
                     !
                  endif
               endif   
               !
               if (nmindsrc(j)>0) exit
               !
            enddo   
            !
         enddo   
         !
         ! Now the sources
         !
         j = nsrc + idrn*2
         !
         xtmp =   cosrot*(xsrc(idrn) - x0) + sinrot*(ysrc(idrn) - y0)
         ytmp = - sinrot*(xsrc(idrn) - x0) + cosrot*(ysrc(idrn) - y0)
         !
         nmindsrc(j) = 0         
         !
         ! Loop through refinement levels (from fine to coarse)
         !
         do iref = nref, 1, -1
            !
            m = int(xtmp/dxr(iref)) + 1
            n = int(ytmp/dyr(iref)) + 1
            !
            ! Find nm index of this src point
            !
            do nm = 1, np
               !
               if (z_flags_iref(nm) == iref) then
                  !
                  if (z_index_z_n(nm)==n .and. z_index_z_m(nm)==m) then
                     !
                     ! Point found
                     !
                     nmindsrc(j) = nm
                     !
                     exit
                     !
                  endif
               endif   
               !
               if (nmindsrc(j)>0) exit
               !
            enddo   
            !
         enddo   
         !
      enddo
      !
      deallocate(xsrc)
      deallocate(ysrc)
      deallocate(xsnk)
      deallocate(ysnk)
      !
   endif
   !
   end subroutine
   !
   !
   !
   subroutine update_discharges(t, dt, tloop)
   !
   ! Update discharges
   !
   use sfincs_data
   !
   implicit none
   !
   integer  :: count0
   integer  :: count1
   integer  :: count_rate
   integer  :: count_max
   real     :: tloop
   !
   real*8           :: t
   real*4           :: dt
   real*4           :: qq
   real*4           :: qq0
   !
   integer isrc, itsrc, idrn, jin, jout, nmin, nmout
   !
   call system_clock(count0, count_rate, count_max)
   !
   ! Compute instantaneous discharges from point sources
   !
   if (nsrc>0) then
      do itsrc = itsrclast, ntsrc
         ! Find first point in time series large than t
         if (tsrc(itsrc)>t) then
            do isrc = 1, nsrc
               qtsrc(isrc) = qsrc(isrc, itsrc - 1) + (qsrc(isrc, itsrc) - qsrc(isrc, itsrc - 1))*(t - tsrc(itsrc - 1))/(tsrc(itsrc) - tsrc(itsrc - 1))
            enddo
            itsrclast = itsrc - 1
            exit
         endif
      enddo
      !
      !$acc update device(qtsrc), async(1)
      !
   endif
   !
   if (ndrn>0) then
      !
      !$acc serial, present( z_volume, zs, zb, nmindsrc, qtsrc, drainage_type, drainage_params ), async(1)
      do idrn = 1, ndrn
         !
         jin  = nsrc + idrn*2 - 1
         jout = nsrc + idrn*2
         !
         nmin  = nmindsrc(jin)
         nmout = nmindsrc(jout)
         !
         qtsrc(jin)  = 0.0
         qtsrc(jout) = 0.0
         !
         if (nmin>0 .and. nmout>0) then
            !
            select case(drainage_type(idrn))
               !
               case(1)
                  !
                  ! Pump
                  !
                  qq = drainage_params(idrn,1)
                  !
                  if (qq>=0.0) then
                     !
                     ! Constant discharge, no need to change
                     !
               else
                  !
                  ! Time series
                  !
               endif   
               !
               if (subgrid) then
                  qq = max(min(qq, z_volume(nmin)/dt), 0.0)
               else
                  qq = max(min(qq, (zs(nmin) - zb(nmin))*area/dt), 0.0)
               endif
               !
               qtsrc(jin)  = -qq
               qtsrc(jout) = qq
               !
            case(2)
               !
               ! Culvert
               !
               qq0 = -qtsrc(jin) ! Previous time step, directed from intake to outfall
               !
               if (zs(nmin)>zs(nmout)) then
                  !
                  qq  = drainage_params(idrn,1)*sqrt(zs(nmin) - zs(nmout))
                  !
               else
                  !
                  qq  = -drainage_params(idrn,1)*sqrt(zs(nmout) - zs(nmin))
                  !
               endif
               !
               if (subgrid) then
                  if (qq>0.0) then
                     qq = min(qq, max(z_volume(nmin),0.0)/dt)
                  else
                     qq = max(qq, -max(z_volume(nmout),0.0)/dt)
                  endif
               else
                  if (qq>0.0) then
                     qq = min(qq, max((zs(nmin) - zb(nmin))*area,0.0)/dt)
                  else
                     qq = max(qq, -max((zs(nmout) - zb(nmout))*area,0.0)/dt)
                  endif
               endif
               !
               ! Add some relaxation
               !
               qq = 0.10*qq + 0.90*qq0
               !
               qtsrc(jin)  = -qq
               qtsrc(jout) =  qq
               !
            case(3)
               !
               ! Check valve
               !
               qq0 = -qtsrc(jin) ! Previous time step, directed from intake to outfall
               !
               if (zs(nmin)>zs(nmout)) then
                  !
                  qq  = drainage_params(idrn,1)*sqrt(zs(nmin) - zs(nmout))
                  !
               else
                  !
                  qq  = -drainage_params(idrn,1)*sqrt(zs(nmout) - zs(nmin))
                  !
               endif
               !
               if (subgrid) then
                  if (qq>0.0) then
                     qq = min(qq, max(z_volume(nmin),0.0)/dt)
                  else
                     qq = max(qq, -max(z_volume(nmout),0.0)/dt)
                  endif
               else
                  if (qq>0.0) then
                     qq = min(qq, max((zs(nmin) - zb(nmin))*area,0.0)/dt)
                  else
                     qq = max(qq, -max((zs(nmout) - zb(nmout))*area,0.0)/dt)
                  endif
               endif
               !
               ! Add some relaxation
               !
               qq = 0.10*qq + 0.90*qq0
               !
               ! Make sure it can only flow from intake to outfall point
               !
               qq = max(qq, 0.0)
               qtsrc(jin)  = -qq
               qtsrc(jout) =  qq
               !
            end select
         endif   
         !
      enddo
      !$acc end serial
      !
   endif
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine

end module
