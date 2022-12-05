module sfincs_meteo

contains

   subroutine read_meteo_data()  
   !
   ! Read different meteo input types
   !
   use sfincs_data
   use sfincs_spiderweb
   use sfincs_ncinput
   !
   implicit none
   !   
   integer it, irow, icol, stat, iw, icd
   !
   real*4 dummy, wnd, xx, yy
   !
   write(*,*)'Reading meteo data ...'
   !
   spw_precip = .false.
   !
   if (spwfile(1:4) /= 'none') then
      !
      write(*,*)'Reading spiderweb file ...'
      !  
      call read_spw_dimensions(spwfile,spw_nt,spw_nrows,spw_ncols,spw_radius,spw_nquant)
      !
      ! Allocate
      !
      allocate(spw_times(spw_nt))
      allocate(spw_xe(spw_nt))
      allocate(spw_ye(spw_nt))
      allocate(spw_vmag(spw_nt, spw_nrows, spw_ncols))
      allocate(spw_vdir(spw_nt, spw_nrows, spw_ncols))
      allocate(spw_pdrp(spw_nt, spw_nrows, spw_ncols))
      allocate(spw_prcp(spw_nt, spw_nrows, spw_ncols))
      allocate(spw_wu(spw_nt, spw_nrows, spw_ncols))
      allocate(spw_wv(spw_nt, spw_nrows, spw_ncols))
      allocate(spw_pdrp01(spw_nrows, spw_ncols))
      allocate(spw_prcp01(spw_nrows, spw_ncols))
      allocate(spw_wu01(spw_nrows, spw_ncols))
      allocate(spw_wv01(spw_nrows, spw_ncols))
      !
      call read_spw_file(spwfile,spw_nt,spw_nrows,spw_ncols,spw_radius,spw_times,spw_xe,spw_ye,spw_vmag,spw_vdir,spw_pdrp,spw_prcp,spw_nquant,trefstr)
      !
      dradspw = spw_radius/spw_nrows
      dphispw = 2*pi/spw_ncols
      !      
      ! Compute u and v components of wind
      !
      do it = 1, spw_nt
         do irow = 1, spw_nrows
            do icol = 1, spw_ncols
               spw_wu(it, irow, icol) = spw_vmag(it, irow, icol)*cos(pi*(270.0 - spw_vdir(it, irow, icol))/180)
               spw_wv(it, irow, icol) = spw_vmag(it, irow, icol)*sin(pi*(270.0 - spw_vdir(it, irow, icol))/180)
            enddo
         enddo
      enddo
      !
      if (spw_nquant==4) then
          spw_precip = .true.
      endif
      !
      if (utmzone/='nil') then
         !
         ! Convert spiderweb coordinates to utm zone
         ! 
         do it = 1, spw_nt
            call deg2utm(spw_ye(it),spw_xe(it),xx,yy,utmzone)
            spw_xe(it) = xx
            spw_ye(it) = yy
         enddo
         !
      endif
      !
   endif
   ! 
   if (amufile(1:4) /= 'none') then
      !
      write(*,*)'Reading amu and amv file ...'
      !  
      call read_amuv_dimensions(amufile,amuv_nt,amuv_nrows,amuv_ncols,amuv_x_llcorner,amuv_y_llcorner,amuv_dx,amuv_dy,amuv_nquant)
      !
      ! Allocate
      !
      allocate(amuv_times(amuv_nt))
      allocate(amuv_wu(amuv_nt, amuv_nrows, amuv_ncols))
      allocate(amuv_wv(amuv_nt, amuv_nrows, amuv_ncols))
      allocate(amuv_wu01(amuv_nrows, amuv_ncols))
      allocate(amuv_wv01(amuv_nrows, amuv_ncols))
      !
      call read_amuv_file(amufile,amuv_nt,amuv_nrows,amuv_ncols,amuv_times,amuv_wu,trefstr)
      call read_amuv_file(amvfile,amuv_nt,amuv_nrows,amuv_ncols,amuv_times,amuv_wv,trefstr)
      !
   elseif (netamuamvfile(1:4) /= 'none') then   ! FEWS compatible Netcdf amu&amv wind spatial input
      !
      call read_netcdf_amuv_data()
      !
      allocate(amuv_wu01(amuv_nrows, amuv_ncols)) !(has to be allocated somewhere)
      allocate(amuv_wv01(amuv_nrows, amuv_ncols))
      !
   endif
   ! 
   if (amprfile(1:4) /= 'none') then
      !
      write(*,*)'Reading ampr file ...'
      !  
      call read_amuv_dimensions(amprfile,ampr_nt,ampr_nrows,ampr_ncols,ampr_x_llcorner,ampr_y_llcorner,ampr_dx,ampr_dy,ampr_nquant)
      !
      ! Allocate
      !
      allocate(ampr_times(ampr_nt))
      allocate(ampr_pr(ampr_nt, ampr_nrows, ampr_ncols))
      allocate(ampr_pr01(ampr_nrows, ampr_ncols))
      !
      call read_amuv_file(amprfile,ampr_nt,ampr_nrows,ampr_ncols,ampr_times,ampr_pr,trefstr)
      !
   elseif (netamprfile(1:4) /= 'none') then   ! FEWS compatible Netcdf ampr precipitation spatial input
      !
      call read_netcdf_ampr_data()
      !
      allocate(ampr_pr01(ampr_nrows, ampr_ncols))!(has to be allocated somewhere)
      !
   endif
   !
   if (ampfile(1:4) /= 'none') then
      !
      write(*,*)'Reading amp file ...'
      !  
      call read_amuv_dimensions(ampfile,amp_nt,amp_nrows,amp_ncols,amp_x_llcorner,amp_y_llcorner,amp_dx,amp_dy,amp_nquant)
      !
      ! Allocate
      !
      allocate(amp_times(amp_nt))
      allocate(amp_patm(amp_nt, amp_nrows, amp_ncols))
      allocate(amp_patm01(amp_nrows, amp_ncols))
      !
      call read_amuv_file(ampfile,amp_nt,amp_nrows,amp_ncols,amp_times,amp_patm,trefstr)
      !
   elseif (netampfile(1:4) /= 'none') then   ! FEWS compatible Netcdf amp barometric pressure spatial input
      !
      call read_netcdf_amp_data()
      !
      allocate(amp_patm01(amp_nrows, amp_ncols)) ! (has to be allocated somewhere)
      !
   endif
   !
   if (wndfile(1:4) /= 'none') then
      !
      ! Wind in time series file 
      write(*,*)'Reading ', trim(wndfile)              
      !
      ntwnd = 0
      itwndlast = 1
      !
      open(500, file=trim(wndfile))       
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         ntwnd = ntwnd + 1
      enddo 
      rewind(500) 
      !
      allocate(twnd(ntwnd))
      allocate(wndmag(ntwnd))
      allocate(wnddir(ntwnd))
      !
      do it = 1, ntwnd         
         read(500,*)twnd(it),wndmag(it),wnddir(it)
         wnddir(it) = pi*(270.0 - wnddir(it))/180 - rotation ! Convert to cartesian, where the wind is blowing to and include grid rotation
      enddo
      !
      close(500)
      !
   endif   
   !
   if (prcpfile(1:4) /= 'none') then
      !
      ! Rainfall in time series file 
      write(*,*)'Reading ',trim(prcpfile)       
      !
      ntprcp = 0 
      itprcplast = 1
      !
      open(500, file=trim(prcpfile))       
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         ntprcp = ntprcp + 1
      enddo 
      rewind(500) 
      !
      allocate(tprcpt(ntprcp))
      allocate(tprcpv(ntprcp))
      !
      do it = 1, ntprcp  
         read(500,*)tprcpt(it),tprcpv(it)
      enddo
      !
      close(500)
      !
   endif   
   !
   ! Determine Cd coefficients
   !
   do iw = 1, 1000
      wnd = iw*0.1
      do icd = 1, cd_nr - 1
         cdval(iw) = cd_val(icd)
         if (wnd>=cd_wnd(icd) .and. wnd<cd_wnd(icd + 1)) then
            cdval(iw) = cd_val(icd) + (wnd - cd_wnd(icd))*(cd_val(icd + 1) - cd_val(icd))/(cd_wnd(icd + 1) - cd_wnd(icd))
            exit
         elseif (wnd>=cd_wnd(cd_nr)) then  
            cdval(iw) = cd_val(cd_nr)
            exit
         endif
      enddo   
   enddo
   !   
   end subroutine
   !
   !
   !
   subroutine update_spiderweb_data()   
   !
   use sfincs_data
   !
   implicit none
   !   
   real*4                           :: spw_xe01
   real*4                           :: spw_ye01
   real*4, dimension(4)             :: f
   !
   integer itspw, itw0, itw1, idstspw, iphispw, itw, irow, icol, nm, n, m, ip, ivm, iwa, idir, id
   integer, dimension(4)            :: ind1
   integer, dimension(4)            :: ind2
   !
   real*4 x, y, dstspw
   real*4 di1, dj1, phispw
   real*4 wup, wvp, prp, pcp, vmag, cd
   real*4 meteo_t
   real*4 twfac
   real*4 merge_frac
   real*4 wa
   real*4 z0marine
   real*4 z0l
   real*4 fd
   real*4 wdir
   real*4 facint
   !
   do itw = 1, 2
      !      
      ! Find time indices in spw file
      !
      if (itw==1) then
         meteo_t = meteo_t0
      else
         meteo_t = meteo_t1
      endif
      !
      do itspw = 1, spw_nt
         if (spw_times(itspw)<=meteo_t) then
            itw0 = itspw
         endif
      enddo
      itw1 = itw0 + 1
      itw1 = min(itw1, spw_nt)
      !
      meteo_t = min(meteo_t, spw_times(itw1))
      !
      twfac  = (meteo_t - spw_times(itw0))/max(spw_times(itw1) - spw_times(itw0), 1.0e-6)
      !
      ! Eye
      ! 
      spw_xe01 = spw_xe(itw0)*(1.0 - twfac) + spw_xe(itw1)*twfac
      spw_ye01 = spw_ye(itw0)*(1.0 - twfac) + spw_ye(itw1)*twfac
      !
      ! Wind, pressure, precipitation
      !
      do irow = 1, spw_nrows
         do icol = 1, spw_ncols
            !
            spw_wu01(irow, icol) = spw_wu(itw0, irow, icol)*(1.0 - twfac)   + spw_wu(itw1, irow, icol)*twfac
            spw_wv01(irow, icol) = spw_wv(itw0, irow, icol)*(1.0 - twfac)   + spw_wv(itw1, irow, icol)*twfac
            if (patmos) then
               spw_pdrp01(irow, icol) = spw_pdrp(itw0, irow, icol)*(1.0 - twfac) + spw_pdrp(itw1, irow, icol)*twfac
            endif
            if (spw_precip) then
               spw_prcp01(irow, icol) = spw_prcp(itw0, irow, icol)*(1.0 - twfac) + spw_prcp(itw1, irow, icol)*twfac
            endif
            !
         enddo
      enddo
      !
      ! Loop through grid points
      !
      ! TODO : add openmp and openacc code
      !
      do nm = 1, np
         !
         x = z_xz(nm)
         y = z_yz(nm)      
         !
         ! Compute tauwu0, tauwv0, pdrp0, prcp0 at spw_t0
         !
         if (crsgeo) then
            dstspw  = sqrt(((cos(spw_ye01*pi/180)*111111)*(x - spw_xe01))**2 + (111111*(y - spw_ye01))**2) ! Distance to eye
         else
            dstspw  = sqrt((x - spw_xe01)**2 + (y - spw_ye01)**2) ! Distance to eye
         endif
         !
         ! Initialize meteo data, but only if we don't also use background meteo
         !
         if (itw==1) then
            !
            if (amufile(1:4) == 'none' .and. netamuamvfile(1:4) == 'none') then            
               !
               tauwu0(nm) = 0.0
               tauwv0(nm) = 0.0
               !
               if (store_meteo) then
                  windu0(nm) = 0.0
                  windv0(nm) = 0.0
               endif   
               !
            endif
            !
            if (patmos) then
               if (ampfile(1:4) == 'none') then            
                  patm0(nm)  = gapres
               endif
            endif   
            !
            if (precip .and. spw_precip) then
              if (amprfile(1:4) == 'none' .and. netamprfile(1:4) == 'none') then            
                  prcp0(nm)  = 0.0 ! m/s
               endif
            endif
            !
         else
            !
            if (amufile(1:4) == 'none' .and. netamuamvfile(1:4) == 'none') then            
               !
               tauwu1(nm) = 0.0
               tauwv1(nm) = 0.0
               !
               if (store_meteo) then
                  windu1(nm) = 0.0
                  windv1(nm) = 0.0
               endif   
               !
            endif
            !
            if (patmos) then
               if (ampfile(1:4) == 'none') then            
                  patm1(nm)  = gapres
               endif
            endif   
            !
            if (precip .and. spw_precip) then
               if (amprfile(1:4) == 'none' .and. netamprfile(1:4) == 'none') then            
                  prcp1(nm)  = 0.0 ! m/s
               endif
            endif
            !
         endif
         !          
         if (dstspw>spw_radius) cycle ! Point outside spiderweb
         !
         ! Determine row indices
         !
         idstspw = int(dstspw/dradspw)
         idstspw = max(idstspw, 1)
         ind1(1) = idstspw
         ind1(2) = idstspw
         ind1(3) = idstspw + 1
         ind1(4) = idstspw + 1
         if (ind1(3)>spw_nrows) cycle             
         dj1     = (dstspw - dradspw*idstspw) / dradspw
         phispw  = 0.5*pi - atan2(y - spw_ye01, x - spw_xe01) ! Geographic
         phispw  = modulo(phispw, 2*pi)
         !
         ! Determine column indices
         !
         iphispw = min(int(phispw/dphispw) + 1, spw_ncols)
         ind2(1) = iphispw
         ind2(2) = iphispw + 1
         ind2(3) = iphispw + 1
         ind2(4) = iphispw
         if (ind2(2)>spw_ncols) ind2(2) = 1
         if (ind2(3)>spw_ncols) ind2(3) = 1
         di1     = (phispw - dphispw*(iphispw-1)) / dphispw
         !
         ! Weight factors
         !    
         f(1) = (1.0 - di1) * (1.0 - dj1)
         f(2) = (      di1) * (1.0 - dj1)
         f(3) = (      di1) * (      dj1) 
         f(4) = (1.0 - di1) * (      dj1)
         !
         ! Wind
         !
         wup = 0.0
         wvp = 0.0
         prp = 0.0
         pcp = 0.0
         !                 
         do ip = 1, 4
            !
            wup     = wup + f(ip)*spw_wu01(ind1(ip),ind2(ip))
            wvp     = wvp + f(ip)*spw_wv01(ind1(ip),ind2(ip))
            !
            if (patmos) then
               pcp  = pcp + f(ip)*spw_pdrp01(ind1(ip),ind2(ip))
            endif
            !
            if (precip .and. spw_precip) then
               prp  = prp + f(ip)*spw_prcp01(ind1(ip),ind2(ip))
            endif
            !
         enddo
         !
         ! Compute wind stress
         !
         vmag   = sqrt(wup**2 + wvp**2)
         !
         if (waveage>0) then
            !
            ! Determine Cd with wave age based on LGX method
            !
            wa = max(min(waveage, 1.99), 0.1001)
            !
            ivm = max(int(vmag/2.0), 1)
            iwa = int(wa/0.1)
            cd = 0.001*cdlgx(ivm, iwa)
            !
         else
            !
            cd = cdval(int(vmag*10) + 1)
            !
         endif
         !
         ! Merge frac for merging with background winds
         !
         if (dstspw>spw_merge_frac*spw_radius) then
            merge_frac = 2*(spw_radius - dstspw)/spw_radius
         else
            merge_frac = 1.0
         endif
         !
         if (wind_reduction) then
            !
            ! Reduction of wind speed over land (Westerink et al., 2008. A Basin- to Channel-Scale Unstructured Grid Hurricane Storm Surge Model Applied to Southern Louisiana)
            !
            z0marine = 0.001835*cd*vmag**2 ! 0.018*cd*vmag**2/9.81, should maybe use a look-up table for this
            !
            ! Determine z0land
            !
            wdir   = 270.0 - atan2(wvp, wup)*180/pi ! Nautical degrees
            if (wdir<0.0)    wdir = wdir + 360.0
            if (wdir>=360.0) wdir = wdir - 360.0
            idir   = int(wdir/30) + 1
            facint = (wdir - (idir - 1)*30.0) / 30
            if (idir<12) then
               z0l    = z0land(idir, nm) + (z0land(idir + 1, nm) - z0land(idir, nm))*facint
            else
               z0l    = z0land(  12, nm) + (z0land(       1, nm) - z0land(  12, nm))*facint
            endif   
            !
            if (z0l>z0marine) then
               !
               ! Should really use look-up table for this next very non-linear bit
               !
               id   = min(int(0.2*z0l/z0marine) + 1, 101)
               fd   = z0land_table(id) ! fd   = (z0marine/z0l)**0.0706
               vmag = fd*vmag
               wup  = fd*wup
               wvp  = fd*wvp
               !
            endif
            !
         endif
         !
         if (itw==1) then
            tauwu0(nm) = (1.0 - merge_frac)*tauwu0(nm) + merge_frac*vmag*( cosrot*wup + sinrot*wvp)*rhoa*cd/rhow
            tauwv0(nm) = (1.0 - merge_frac)*tauwv0(nm) + merge_frac*vmag*(-sinrot*wup + cosrot*wvp)*rhoa*cd/rhow
            if (patmos) then
               patm0(nm)  = (1.0 - merge_frac)*patm0(nm) + merge_frac*(gapres - pcp)
            endif
            if (precip .and. spw_precip) then
               prcp0(nm)  = (1.0 - merge_frac)*prcp0(nm) + merge_frac*prp/(1000*3600) ! m/s
            endif
            if (store_meteo) then
               windu0(nm) = (1.0 - merge_frac)*windu0(nm) + merge_frac*wup
               windv0(nm) = (1.0 - merge_frac)*windv0(nm) + merge_frac*wvp
            endif   
         else
            tauwu1(nm) = (1.0 - merge_frac)*tauwu1(nm) + merge_frac*vmag*( cosrot*wup + sinrot*wvp)*rhoa*cd/rhow
            tauwv1(nm) = (1.0 - merge_frac)*tauwv1(nm) + merge_frac*vmag*(-sinrot*wup + cosrot*wvp)*rhoa*cd/rhow
            if (patmos) then
               patm1(nm)  = (1.0 - merge_frac)*patm1(nm) + merge_frac*(gapres - pcp)
            endif
            if (precip .and. spw_precip) then
               prcp1(nm)  = (1.0 - merge_frac)*prcp1(nm) + merge_frac*prp/(1000*3600) ! m/s
            endif
            if (store_meteo) then
               windu1(nm) = (1.0 - merge_frac)*windu1(nm) + merge_frac*wup
               windv1(nm) = (1.0 - merge_frac)*windv1(nm) + merge_frac*wvp
            endif   
         endif
         !
      enddo                          
   enddo
   !
   end subroutine
   !
   !
   !
   subroutine update_amuv_data()   
   !
   use sfincs_data
   !
   implicit none
   !   
   real*4, dimension(4)             :: f
   !
   integer itspw, itw0, itw1, ix, iy, itw, irow, icol, nm, n, m, ip, ivm, iwa
   integer, dimension(4)            :: ind1
   integer, dimension(4)            :: ind2
   !
   real*4 x, y, yul, xll
   real*4 di1, dj1
   real*4 wup, wvp, vmag, cd
   real*4 meteo_t
   real*4 twfac
   real*4 wa
   !
   ! Find time indices in amu file
   !
   do itw = 1, 2
      !
      if (itw==1) then
         meteo_t = meteo_t0
      else
         meteo_t = meteo_t1
      endif
      !
      do itspw = 1, amuv_nt
         if (amuv_times(itspw)<=meteo_t) then
            itw0 = itspw
         endif
      enddo
      itw1 = itw0 + 1
      itw1 = min(itw1, amuv_nt)
      !
      meteo_t = min(meteo_t, amuv_times(itw1))
      !
      twfac  = (meteo_t - amuv_times(itw0))/max(amuv_times(itw1) - amuv_times(itw0), 1.0e-6)
      !
      do irow = 1, amuv_nrows
         do icol = 1, amuv_ncols
            !
            amuv_wu01(irow, icol) = amuv_wu(itw0, irow, icol)*(1.0 - twfac)   + amuv_wu(itw1, irow, icol)*twfac
            amuv_wv01(irow, icol) = amuv_wv(itw0, irow, icol)*(1.0 - twfac)   + amuv_wv(itw1, irow, icol)*twfac
            !
         enddo
      enddo
      !
      ! Loop through grid points
      !
      yul = amuv_y_llcorner + (amuv_nrows - 1)*amuv_dy + 0.5*amuv_dy
      xll = amuv_x_llcorner + 0.5*amuv_dx
      !
      do nm = 1, np
         !
         x = z_xz(nm)
         y = z_yz(nm)      
         !
         ! Determine row indices
         !
         iy = int((yul - y)/amuv_dy) + 1
         ind1(1) = iy
         ind1(2) = iy
         ind1(3) = iy + 1
         ind1(4) = iy + 1
         dj1     = (yul - (iy - 1)*amuv_dy - y) / amuv_dy
         !
         ! Determine column indices
         !
         ix = int((x - xll)/amuv_dx) + 1
         ind2(1) = ix
         ind2(2) = ix + 1
         ind2(3) = ix + 1
         ind2(4) = ix
         di1     = (x - (xll + (ix - 1)*amuv_dx)) / amuv_dx
         !
         if (iy<=1 .or. iy>=amuv_nrows .or. ix<=0 .or. ix>=amuv_ncols) then
            !
            ! Point outside rectangular amu/amv grid, this should really produce an error ...
            !
            if (itw==1) then
               tauwu0(nm) = 0.0
               tauwv0(nm) = 0.0
            else
               tauwu1(nm) = 0.0
               tauwv1(nm) = 0.0
            endif
            !
            cycle
            !
         endif
         !            
         ! Weight factors
         !    
         f(1) = (1.0 - di1) * (1.0 - dj1)
         f(2) = (      di1) * (1.0 - dj1)
         f(3) = (      di1) * (      dj1) 
         f(4) = (1.0 - di1) * (      dj1)
         !
         ! Wind
         !
         wup = 0.0
         wvp = 0.0
         !                 
         do ip = 1, 4
            !
            wup     = wup + f(ip)*amuv_wu01(ind1(ip),ind2(ip))
            wvp     = wvp + f(ip)*amuv_wv01(ind1(ip),ind2(ip))
            !
         enddo
         !
         ! Compute wind stress
         !
         vmag   = sqrt(wup**2 + wvp**2)
         !
         if (waveage>0) then
            !
            ! Determine Cd with wave age based on LGX method
            !
            wa = max(min(waveage, 1.99), 0.1001)
            !            
            ivm = max(int(vmag/2.0), 1)
            iwa = int(wa/0.1)
            cd = 0.001*cdlgx(ivm, iwa)                  
            !
         else
            !
            cd = cdval(int(vmag*10)+1)
            !
         endif
         !
         !
         if (radstr) then
            !
            ! Treat wind speeds as radiation stresses (temporary solution for including wave forcing)  
            !   
            if (itw==1) then
               tauwu0(nm) = ( cosrot*wup + sinrot*wvp)/rhow
               tauwv0(nm) = (-sinrot*wup + cosrot*wvp)/rhow
            else
               tauwu1(nm) = ( cosrot*wup + sinrot*wvp)/rhow
               tauwv1(nm) = (-sinrot*wup + cosrot*wvp)/rhow
            endif
            !
         else   
            !
            if (itw==1) then
               tauwu0(nm) = vmag*( cosrot*wup + sinrot*wvp)*rhoa*cd/rhow
               tauwv0(nm) = vmag*(-sinrot*wup + cosrot*wvp)*rhoa*cd/rhow
               if (store_meteo) then
                  windu0(nm) = wup
                  windv0(nm) = wvp
               endif   
            else
               tauwu1(nm) = vmag*( cosrot*wup + sinrot*wvp)*rhoa*cd/rhow
               tauwv1(nm) = vmag*(-sinrot*wup + cosrot*wvp)*rhoa*cd/rhow
               if (store_meteo) then
                  windu1(nm) = wup
                  windv1(nm) = wvp
               endif   
            endif
            !
         endif
         !
      enddo                          
   enddo
   !
   end subroutine


   subroutine update_amp_data()   
   !
   use sfincs_data
   !
   implicit none
   !   
   real*4, dimension(4)             :: f
   !
   integer itspw, itw0, itw1, ix, iy, itw, irow, icol, nm, n, m, ip
   integer, dimension(4)            :: ind1
   integer, dimension(4)            :: ind2
   !
   real*4 x, y, yul, xll
   real*4 di1, dj1
   real*4 pr
   real*4 meteo_t
   real*4 twfac
   !
   do itw = 1, 2
      !
      if (itw==1) then
         meteo_t = meteo_t0
      else
         meteo_t = meteo_t1
      endif
      !
      ! Find time indices in amp file
      !
      do itspw = 1, amp_nt
         if (amp_times(itspw)<=meteo_t) then
            itw0 = itspw
         endif
      enddo
      itw1 = itw0 + 1
      itw1 = min(itw1, amp_nt)
      !
      meteo_t = min(meteo_t, amp_times(itw1))
      !
      twfac  = (meteo_t - amp_times(itw0))/max(amp_times(itw1) - amp_times(itw0), 1.0e-6)
      !
      do irow = 1, amp_nrows
         do icol = 1, amp_ncols
            !
            amp_patm01(irow, icol) = amp_patm(itw0, irow, icol)*(1.0 - twfac)   + amp_patm(itw1, irow, icol)*twfac
            !
         enddo
      enddo
      !
      ! Loop through grid points
      !
      yul = amp_y_llcorner + (amp_nrows - 1)*amp_dy + 0.5*amp_dy
      xll = amp_x_llcorner + 0.5*amp_dx
      !
      do nm = 1, np
         !
         x = z_xz(nm)
         y = z_yz(nm)      
         !
         ! Determine row indices
         !
         iy = int((yul - y)/amp_dy) + 1
         ind1(1) = iy
         ind1(2) = iy
         ind1(3) = iy + 1
         ind1(4) = iy + 1
         dj1     = (yul - (iy - 1)*amp_dy - y) / amp_dy
         !
         ! Determine column indices
         !
         ix = int((x - xll)/amp_dx) + 1
         ind2(1) = ix
         ind2(2) = ix + 1
         ind2(3) = ix + 1
         ind2(4) = ix
         di1     = (x - (xll + (ix - 1)*amp_dx)) / amp_dx
         !
         if (iy<=1 .or. iy>=amp_nrows .or. ix<=0 .or. ix>=amp_ncols) then
            !
            ! Point outside rectangular amp grid
            !
            if (itw==1) then
               !
               patm0(nm) = gapres
               !
            else
               !
               patm1(nm) = gapres
               !
            endif
            !
            cycle
            !
         endif
         !            
         ! Weight factors
         !    
         f(1) = (1.0 - di1) * (1.0 - dj1)
         f(2) = (      di1) * (1.0 - dj1)
         f(3) = (      di1) * (      dj1) 
         f(4) = (1.0 - di1) * (      dj1)
         !
         ! Wind
         !
         pr = 0.0
         !                 
         do ip = 1, 4
            !
            pr = pr + f(ip)*amp_patm01(ind1(ip),ind2(ip))
            !
         enddo
         !
         if (itw==1) then
            patm0(nm)  = pr
         else
            patm1(nm)  = pr
         endif
         !
      enddo                          
      !
   enddo
   !
   end subroutine


   subroutine update_ampr_data()   
   !
   use sfincs_data
   !
   implicit none
   !   
   real*4, dimension(4)             :: f
   !
   integer itspw, itw0, itw1, ix, iy, itw, irow, icol, nm, n, m, ip
   integer, dimension(4)            :: ind1
   integer, dimension(4)            :: ind2
   !
   real*4 x, y, yul, xll
   real*4 di1, dj1
   real*4 prp
   real*4 meteo_t
   real*4 twfac
   !
   do itw = 1, 2
      !
      if (itw==1) then
         meteo_t = meteo_t0
      else
         meteo_t = meteo_t1
      endif
      !
      ! Find time indices in ampr file
      !
      do itspw = 1, ampr_nt
         if (ampr_times(itspw)<=meteo_t) then
            itw0 = itspw
         endif
      enddo
      itw1 = itw0 + 1
      itw1 = min(itw1, ampr_nt)
      !
      meteo_t = min(meteo_t, ampr_times(itw1))
      !
      twfac  = (meteo_t - ampr_times(itw0))/max(ampr_times(itw1) - ampr_times(itw0), 1.0e-6)
      !
      if (ampr_block) then
         !
         ! Use block function for rainfall data from ampr files (default)
         !
         twfac  = 0.0
         !
      endif   
      !
      do irow = 1, ampr_nrows
         do icol = 1, ampr_ncols
            !
            ampr_pr01(irow, icol) = ampr_pr(itw0, irow, icol)*(1.0 - twfac)   + ampr_pr(itw1, irow, icol)*twfac
            !
         enddo
      enddo
      !
      ! Loop through grid points
      !
      yul = ampr_y_llcorner + (ampr_nrows - 1)*ampr_dy + 0.5*ampr_dy
      xll = ampr_x_llcorner + 0.5*ampr_dx
      !
      do nm = 1, np
         !
         x = z_xz(nm)
         y = z_yz(nm)      
         !
         ! Determine row indices
         !
         iy = int((yul - y)/ampr_dy) + 1
         ind1(1) = iy
         ind1(2) = iy
         ind1(3) = iy + 1
         ind1(4) = iy + 1
         dj1     = (yul - (iy - 1)*ampr_dy - y) / ampr_dy
         !
         ! Determine column indices
         !
         ix = int((x - xll)/ampr_dx) + 1
         ind2(1) = ix
         ind2(2) = ix + 1
         ind2(3) = ix + 1
         ind2(4) = ix
         di1     = (x - (xll + (ix - 1)*ampr_dx)) / ampr_dx
         !
         if (iy<=1 .or. iy>=ampr_nrows .or. ix<=0 .or. ix>=ampr_ncols) then
            !
            ! Point outside rectangular ampr grid
            !
            if (itw==1) then
               !
               prcp0(nm) = 0.0
               !
            else
               !
               prcp0(nm) = 0.0
               !
            endif
            !
            cycle
            !
         endif
         !
         ! Weight factors
         !    
         f(1) = (1.0 - di1) * (1.0 - dj1)
         f(2) = (      di1) * (1.0 - dj1)
         f(3) = (      di1) * (      dj1) 
         f(4) = (1.0 - di1) * (      dj1)
         !
         ! Wind
         !
         prp = 0.0
         !                 
         do ip = 1, 4
            !
            prp     = prp + f(ip)*ampr_pr01(ind1(ip),ind2(ip))
            !
         enddo
         !
         if (itw==1) then
            prcp0(nm)  = prp/(1000*3600) ! m/s
         else
            prcp1(nm)  = prp/(1000*3600) ! m/s
         endif
         !
      enddo                          
   enddo
   !
   end subroutine


   subroutine update_meteo_forcing(t, dt, tloop)
   !
   ! Update wind stresses and precipitation (this happens every time step)
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
   real*8                           :: t
   real*4                           :: dt
   real*4                           :: twfact
   real*4                           :: onemintwfact
   real*4                           :: smfac
   real*4                           :: oneminsmfac
   integer                          :: nm, ib
   !
   call system_clock(count0, count_rate, count_max)
   !
   if (meteo3d) then
      !
      twfact  = (t - meteo_t0)/(meteo_t1 - meteo_t0)
      onemintwfact = 1.0 - twfact
      !
      !$omp parallel &
      !$omp private ( nm )
      !$omp do
      !$acc kernels, present(tauwu, tauwv,  tauwu0, tauwv0, tauwu1, tauwv1, &
      !$acc                  windu, windv, windu0, windv0, windu1, windv1, windmax, &
      !$acc                  patm, patm0, patm1, prcp, prcp0, prcp1 ), async(1)
      !$acc loop independent, private(nm)
      do nm = 1, np
         !
         if (wind) then
            !
            tauwu(nm) = tauwu0(nm)*onemintwfact + tauwu1(nm)*twfact
            tauwv(nm) = tauwv0(nm)*onemintwfact + tauwv1(nm)*twfact
            !
            if (store_meteo) then
               !
               windu(nm) = windu0(nm)*onemintwfact + windu1(nm)*twfact
               windv(nm) = windv0(nm)*onemintwfact + windv1(nm)*twfact
               !
               if (store_wind_max) then
                  windmax(nm) = max(windmax(nm), sqrt(windu(nm)**2 + windv(nm)**2))
               endif   
               !
            endif
            !
         endif   
         !
         if (patmos) then
            patm(nm)  = patm0(nm)*onemintwfact  + patm1(nm)*twfact  ! atmospheric pressure (Pa)
         endif   
         !
         if (precip) then
            !
            prcp(nm)    = prcp0(nm)*onemintwfact  + prcp1(nm)*twfact  ! rainfall in m/s !!!
            netprcp(nm) = prcp(nm)
            cumprcp(nm) = cumprcp(nm) + prcp(nm)*dt
            !
         endif   
         !
      enddo   
      !$omp end do
      !$omp end parallel
      !$acc end kernels
      !
      ! Apply spin-up factor
      !
      if (t<tspinup - 1.0e-3 .and. spinup_meteo) then
         !
         smfac = (t - t0)/(tspinup - t0)
         oneminsmfac = 1.0 - smfac
         !
         !$omp parallel &
         !$omp private ( nm )
         !$omp do
         !$acc kernels, present(tauwu, tauwv, patm, prcp ), async(1)
         !$acc loop independent, private(nm)
         do nm = 1, np
            !
            if (wind) then
               tauwu(nm) = tauwu(nm)*smfac
               tauwv(nm) = tauwv(nm)*smfac
            endif   
            !
            if (patmos) then
               patm(nm)  =patm(nm)*smfac + gapres*oneminsmfac
            endif   
            !
            if (precip) then
               netprcp(nm) = netprcp(nm)*smfac
            endif   
            !
         enddo   
         !$omp end do
         !$omp end parallel
         !$acc end kernels
         !
      endif         
      !   
      if (patmos .and. pavbnd>0.0) then
         !
         !$acc serial, present( patmb, nmindbnd, patm ), async(1) 
         do ib = 1, ngbnd
            patmb(ib) = patm(nmindbnd(ib))
         enddo
         !$acc end serial
         !
         ! patmb is used at boundary points in the CPU part of update_boundary_conditions (should try to make this faster)
         !
         !$acc update host(patmb), async(1)
         !
      endif
      !   
   endif
   !
   ! Wind from time series
   !
   if (wndfile(1:4) /= 'none') then
      !
      ! Wind from time series 
      !
      call update_wind_forcing_from_timeseries(t) 
      ! 
   endif
   !
   ! Rainfall from time series
   !
   if (prcpfile(1:4) /= 'none') then
      !
      call update_precipitation_from_timeseries(t, dt) 
      !
   endif
   !
   if (infiltration) then
      !
      ! Compute infiltration rates
      !
      call update_infiltration_map(dt)
      !
   endif   
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !         
   end subroutine


   subroutine update_wind_forcing_from_timeseries(t)
   !
   ! Update values at boundary points
   !
   use sfincs_data
   !
   implicit none
   !
   integer itw, nm
   !
   real*8                           :: t, twu, twv
   !
   real*4 cd, vmag, vdir, dr0, dr1, twfac
   !
   do itw = itwndlast, ntwnd ! Loop in time
      !
      if (twnd(itw)>t) then
         !
         twfac  = (t - twnd(itw - 1))/(twnd(itw) - twnd(itw - 1))
         !
         vmag = wndmag(itw - 1)*(1.0 - twfac) + wndmag(itw)*twfac
         dr0  = modulo(wnddir(itw - 1), 2*pi)
         dr1  = modulo(wnddir(itw ),    2*pi)
         if (dr1>dr0 + pi) then
             dr0 = dr0 + 2*pi
         elseif (dr0>dr1 + pi) then    
             dr1 = dr1 + 2*pi
         endif
         vdir = dr0*(1.0 - twfac) + dr1*twfac
         !
         cd = cdval(int(vmag*10)+1)
         !
         twu = vmag**2*cos(vdir)*rhoa*cd/rhow
         twv = vmag**2*sin(vdir)*rhoa*cd/rhow
         !
         !$omp parallel &
         !$omp private ( nm ) &
         !$omp shared ( tauwu,tauwv )
         !$omp do
         !$acc kernels, present( tauwu, tauwv ), async(1)
         !$acc loop independent, private(nm)
         do nm = 1, np
            tauwu(nm) = twu
            tauwv(nm) = twv
         enddo   
         !$acc end kernels
         !$omp end do
         !$omp end parallel
         !
         itwndlast = itw - 1
         !
         exit
         !
       endif
   enddo
   !
   end subroutine   

   subroutine update_precipitation_from_timeseries(t, dt)
   !
   ! Update values at boundary points
   !
   use sfincs_data
   !
   implicit none
   !
   integer itp, nm
   !
   real*8  :: t
   real*4  :: dt
   real*4  :: twfac
   real*4  :: ptmp
   !
   do itp = itprcplast, ntprcp ! Loop in time
      if (tprcpt(itp)>t) then
         exit
      endif
   enddo
   !
   twfac  = (t - tprcpt(itp - 1))/(tprcpt(itp) - tprcpt(itp - 1))
   !
   ptmp = (tprcpv(itp - 1)*(1.0 - twfac) + tprcpv(itp)*twfac)/(1000*3600) ! rain in m/s
   !
   !$omp parallel &
   !$omp private ( nm )
   !$omp do
   !$acc kernels present( prcp, cumprcp, netprcp )
   do nm = 1, np
      prcp(nm)    = ptmp
      netprcp(nm) = ptmp
      cumprcp(nm) = cumprcp(nm) + prcp(nm)*dt
   enddo   
   !$acc end kernels
   !$omp end do
   !$omp end parallel
   !
   itprcplast = itp - 1
   !
   end subroutine   


   subroutine update_infiltration_map(dt)
   !
   ! Update infiltration rates in each grid cell
   !
   use sfincs_data
   !
   implicit none
   !
   integer nm
   !
   real*4  :: Qq
   real*4  :: I
   real*4  :: dt   
   !
   if (inftype == 'con' .or. inftype == 'c2d') then
      !
      ! Infiltration rate map stays constant
      !
      !$omp parallel &
      !$omp private ( nm )
      !$omp do
      do nm = 1, np
         !
         qinfmap(nm) = qinffield(nm)
         !
         ! No infiltration if there is no water
         !  
         if (subgrid) then
            if (z_volume(nm)<=0.0) then
               qinfmap(nm) = 0.0
            endif
         else
            if (zs(nm)<=zb(nm)) then
               qinfmap(nm) = 0.0
            endif
         endif
         !
         cuminf(nm) = cuminf(nm) + qinfmap(nm)*dt
         netprcp(nm) = netprcp(nm) - qinfmap(nm)
         !
      enddo
      !$omp end do
      !$omp end parallel
      !
   elseif (inftype == 'cna') then
      !
      ! Determine infiltration rate with Curve Number (old method; no recovery)
      !
!      !$acc update host(cumprcp), async(1)
      !
      !$omp parallel &
      !$omp private ( Qq,I,nm )
      !$omp do
      do nm = 1, np
         !
         ! Check if Ia (0.2 x S) is larger than cumulative rainfall
         !
         if (cumprcp(nm) > sfacinf*qinffield(nm)) then ! qinffield is S
            ! 
            ! Compute runoff as function of rain
            !
            Qq  = (cumprcp(nm) - sfacinf*qinffield(nm))**2 / (cumprcp(nm) + (1.0 - sfacinf)*qinffield(nm))  ! cumulative runoff in m
            I   = cumprcp(nm) - Qq                      ! cumulative infiltration in m
            qinfmap(nm) = (I - cuminf(nm))/dt           ! infiltration in m/s
            !
         else
            !
            ! Everything still infiltrating
            !
            qinfmap(nm) = prcp(nm)
            !
         endif   
         ! 
         cuminf(nm) = cuminf(nm) + qinfmap(nm)*dt
         netprcp(nm) = netprcp(nm) - qinfmap(nm)
         !
      enddo
      !$omp end do
      !$omp end parallel
      !
      !$acc update device(qinfmap), async(1)
      !
      ! For now, curve number infiltration is done on the CPU 
!      !$acc update device(cuminf), async(1)
      !
      ! Provide update to user
      !write(*,'(a,f5.3,a,f10.3,a)') '  update from SCS-CN method: Average cumulative rainfall of ',sum(cumprcp)/size(cumprcp),' meter and infiltration rate of ',sum(qinfmap*3.6e3*1.0e3)/size(qinfmap) ,' mm/hour'
      !   
   elseif (inftype == 'cnb') then
      !
      ! Determine infiltration rate with Curve Number with recovery
!      !$acc update host(cumprcp), async(1)
      !
      !$omp parallel &
      !$omp private ( Qq,I,nm ) &
      !$omp shared ( scs_P1,scs_F1,scs_Se,scs_rain,cuminf,qinfmap,qinffield,qinffield2 )
      !$omp do
      do nm = 1, np
         !
         ! If there is precip in this grid cell for this time step  
         !
         if (prcp(nm) > 0.0) then
            !
            ! Is raining now
            !
            if (scs_rain(nm) == 1) then
               ! if it was raining before; do nothing
            else
               ! initalise these variables 
               scs_P1(nm)          = 0.0               ! cumulative rainfall for this 'event'
               scs_F1(nm)          = 0.0               ! cumulative infiltration for this 'event'
               scs_Se(nm)          = qinffield2(nm)    ! S for this 'event'
               scs_rain(nm)        = 1                 ! logic used to determine if there is an event ongoing
            endif
            ! 
            !  Compute cum rainfall
            ! 
            scs_P1(nm) = scs_P1(nm) + prcp(nm)*dt
            ! 
            ! Compute runoff
            ! 
            if (scs_P1(nm) > (sfacinf* scs_Se(nm)) ) then ! scs_Se is S
               Qq              = (scs_P1(nm) - sfacinf*scs_Se(nm))**2 / (scs_P1(nm) + (1.0 - sfacinf)*scs_Se(nm))  ! cumulative runoff in m
               I              = scs_P1(nm) - Qq                       ! cum infiltration this event
               qinfmap(nm)    = (I - scs_F1(nm))/dt                ! infiltration in m/s
               scs_F1(nm)     = I                                    ! cum infiltration this event
            else
               Qq              = 0.0                               ! no runoff
               scs_F1(nm)     = scs_P1(nm)                        ! all rainfall is infiltrated
               qinfmap(nm)    = prcp(nm)                          ! infiltration rate = rainfall rate
            endif
            ! 
            ! Compute "remaining S", but note that qinffield2 is not used in computation
            ! 
            qinffield2(nm)  = qinffield2(nm) - qinfmap(nm)*dt
            qinffield2(nm)  = max(qinffield2(nm), 0.0)
            qinfmap(nm)     = max(qinfmap(nm), 0.0)
            !
         else
            ! 
            ! It is not raining here
            ! 
            if (scs_rain(nm) == 1) then
               !
               ! if it was raining before
               ! change logic and set rate to 0
               !
               scs_rain(nm)   = 0
               qinfmap(nm)    = 0.0
               !
            endif
            !
            ! It is not raining and stil not raining
            ! compute recovery of S (not larger than Smax)
            ! note that qinffield2 is S and qinffield is Smax
            ! 
            qinffield2(nm) = qinffield2(nm) + (scs_kr(nm) * qinffield(nm) * dt / 3600)  ! scs_kr is recovery % in hours 
            qinffield2(nm) = min(qinffield2(nm), qinffield(nm))
            !
         endif
         ! 
         cuminf(nm) = cuminf(nm) + qinfmap(nm)*dt
         netprcp(nm) = netprcp(nm) - qinfmap(nm)
         !
      enddo
      !$omp end do
      !$omp end parallel
      !
      !$acc update device(qinfmap), async(1)
      !
   endif
   !
   end subroutine   

   
   subroutine update_meteo_fields(t, tloop)
   !
   ! Update values at boundary points
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
   integer  :: nm
   !
   real*8   :: t
   !
   call system_clock(count0, count_rate, count_max)
   !
   if (amufile(1:4) /= 'none' .or. netamuamvfile(1:4) /= 'none') then
      !
      call update_amuv_data()
      !
      !$acc update device(tauwu0,tauwu1,tauwv0,tauwv1), async(1)
      !
   endif
   !
   if ((ampfile(1:4) /= 'none' .or. netampfile(1:4) /= 'none') .and. patmos) then
      !
      call update_amp_data()
      !
      !$acc update device(patm0,patm1), async(1)
      !
   endif
   !
   if (amprfile(1:4) /= 'none' .or. netamprfile(1:4) /= 'none') then
      !
      call update_ampr_data()
      !
      !$acc update device(prcp0,prcp1), async(1)
      !
   endif
   !      
   if (spwfile(1:4) /= 'none') then
      !
      call update_spiderweb_data()
      !
      !$acc update device(tauwu0,tauwu1,tauwv0,tauwv1,patm0,patm1), async(1)
      !
      if (precip) then
         !
         !$acc update device(prcp0,prcp1), async(1)
         !
      endif
      !
   endif
   !
!   if (precip .and. store_cumulative_precipitation) then
!      !        
!      call update_cumprcp_map()
!      !
!   endif
!   !
!   if (infiltration2d) then
!      !
!      call update_infiltration_map()
!      !
!   endif
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !         
   end subroutine   

end module
