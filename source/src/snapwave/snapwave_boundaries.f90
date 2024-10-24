module snapwave_boundaries

contains   
   
   subroutine read_boundary_data()
   !
   ! Reads bnd, bhs, spw etc. files
   !Boundary conditions > lever tp, Hs en wave dir aand (ook dirspr?) en zet om naar spectrum op de randcellen (ntheta, nsigma) > JONSWAP met gamma als user waarde
   !Of tijdseries als bij FM > block per punt? > voorkeur van Maarten. tim-file
   
   !Voor b.c., sp2 files van XBeach code kopieren?    > een binarire file aanhouden > alvast een version header toevoegen oid
   !Maken matlab converter van sp2 naar binair
   !
   use snapwave_data
   use snapwave_ncinput   
   !
   implicit none
   !
   nwbnd  = 0
   ntwbnd = 0
   update_grid_boundary_points = .true.
   itwbndlast = 2
   !
   if (jonswapfile /= '') then
      !
      ! Read data from timeseries in single point
      !
      call read_boundary_data_singlepoint()
      !
   elseif (netsnapwavefile /= '') then
      ! 
      ! Read space- and time-varying data from netcdf file
      !       
      call read_netcdf_wave_boundary_data()      
      ! 
   else
      !
      ! Read space- and time-varying data
      !
      call read_boundary_data_timeseries()
      !
   endif   
   !
   ! Compute reference table between wave nwbd boundary support points and grid boundary points 
   !
   call find_boundary_indices()
   !
   ! Allocate later needed variables (interpolated values in time):
   !
   ! Allocatie vars for incident waves (same for all input variations):
   allocate(hst_bwv(nwbnd))
   allocate(tpt_bwv(nwbnd))
   allocate(wdt_bwv(nwbnd))
   allocate(dst_bwv(nwbnd))
   !allocate(zst_bwv(nwbnd)) 
   allocate(eet_bwv(ntheta, nwbnd))    
   !
   ! Determine IG wave height and period at boundary point(s) later using Herbers 1994 
   ! (depends on local -changing- water depth, so not calculatable a priori)
   !
   if (igwaves) then
      ! 
      ! Allocate vars for IG: 
      allocate(hst_bwv_ig(nwbnd))
      allocate(tpt_bwv_ig(nwbnd))          
      allocate(eet_bwv_ig(ntheta,nwbnd)) 
      ! 'nwbnd' was determined in read_boundary_data_singlepoint/read_boundary_data_timeseries/read_netcdf_wave_boundary_data
      ! 
      ! Also needed:
      allocate(deptht_bwv(nwbnd))      
      !
   endif    
   !
   ! Convert to cartesian, going-to, radians - independent of input type
   !
   wd_bwv = (270.0 - wd_bwv)*pi/180
   !   
   ! Convert directional spreading input to radians - independent of input type
   ds_bwv = ds_bwv * pi / 180      
   !
   ! Write number of input points sounds - independent of input type
   write(*,*)'Input wave boundary points found: ',nwbnd
   !
   ! Check length time-series - independent of input type
   !
   if ((t_bwv(1) > (t0 + 1.0)) .or. (t_bwv(ntwbnd) < (t1 - 1.0))) then
       ! 
       write(*,'(a)')' WARNING! Times in wave boundary conditions file do not cover entire simulation period!'
       !
   endif      
   !
   end subroutine

   
   subroutine read_boundary_enclosure ()
   !
   ! Reads enc file
   !
   use snapwave_data
   !
   implicit none
   !
   integer n, itb, ib, stat, ifreq
   !
   real*4 dummy,r
   !
   ! Read wave boundaries
   !
   write(*,*)'Reading wave boundary enclosure ', trim(encfile), ' ...'
   open(500, file=trim(encfile)) !as in bwvfile of SFINCS
   do while(.true.)
      read(500,*,iostat = stat)dummy
      if (stat<0) exit
      n_bndenc=n_bndenc+1
   enddo
   rewind(500)
   allocate(x_bndenc(n_bndenc))
   allocate(y_bndenc(n_bndenc))
   do n = 1, n_bndenc
      read(500,*)x_bndenc(n),y_bndenc(n)
   enddo
   close(500)
   !
   end subroutine

   subroutine read_boundary_data_singlepoint()
   !
   ! Reads jonswap and water level file for single point
   ! t(s) Hm0 Tp dir ms zst
   !
   use snapwave_data   
   !
   implicit none
   !
   integer :: irec
   integer :: nrec
   integer :: ier
   real*4  :: dum
   !   
   write(*,*)'Reading boundary file ', trim(jonswapfile), ' ...'
   !
   ! Read jonswap wave time series
   !
   open(11,file=jonswapfile)
   irec=0
   ier=0
   do while (ier==0)
      read(11,*,iostat=ier)dum
      irec = irec + 1
   enddo
   rewind(11)
   nrec = irec - 1
   !
   allocate(t_bwv(nrec))   
   allocate(hs_bwv(1, nrec))
   allocate(tp_bwv(1, nrec))
   allocate(wd_bwv(1, nrec))
   allocate(ds_bwv(1, nrec))
   allocate(zs_bwv(1, nrec))
   !   
   do irec=1,nrec
      read(11,*)t_bwv(irec),hs_bwv(1, irec),tp_bwv(1, irec),wd_bwv(1, irec),ds_bwv(1, irec),zs_bwv(1, irec)      
   enddo
   !
   close(11)
   !
   nwbnd = 1
   !
   ntwbnd = nrec
   !
   end subroutine 

   
   subroutine read_boundary_data_timeseries()
   !
   ! Reads bnd, bhs, spw etc. files
   ! Boundary conditions > lever tp, Hs en wave dir aand (ook dirspr?) en zet om naar spectrum op de randcellen (ntheta, nsigma) > JONSWAP met gamma als user waarde
   !
   use snapwave_data
   !
   implicit none
   !
   integer n, itb, ib, stat, ifreq
   !
   real*4 dummy,r
   !
   ! Read wave boundaries
   !
   if (bndfile(1:4) /= 'none') then    ! Normal ascii input files
      ! temporarily use this input of hs/tp/wavdir/dirspr in separate files, later change to DFM type tim-files
      !
      write(*,*)'Reading wave boundary locations ...'
      !
      open(500, file=trim(bndfile)) !as in bwvfile of SFINCS
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         nwbnd = nwbnd + 1
      enddo
      rewind(500)
      allocate(x_bwv(nwbnd))
      allocate(y_bwv(nwbnd))
      do n = 1, nwbnd
         read(500,*)x_bwv(n),y_bwv(n)
      enddo
      close(500)
      !
      ! Read wave boundaries
      !
      write(*,*)'Reading wave boundaries ...'
      !
      ! Wave time series
      !
      ! First find times in bhs file
      !
      open(500, file=trim(bhsfile))
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         ntwbnd = ntwbnd + 1
      enddo
      close(500)
      !
      allocate(t_bwv(ntwbnd))
      !
      ! Hs (significant wave height)
      ! Times in btp and bwd files must be the same as in bhs file!
      !
      open(500, file=trim(bhsfile))
      allocate(hs_bwv(nwbnd, ntwbnd))
      do itb = 1, ntwbnd
         read(500,*)t_bwv(itb),(hs_bwv(ib, itb), ib = 1, nwbnd)
      enddo
      close(500)
      !
      ! Tp (peak period)
      !
      open(500, file=trim(btpfile))
      allocate(tp_bwv(nwbnd, ntwbnd))
      do itb = 1, ntwbnd
         read(500,*)t_bwv(itb),(tp_bwv(ib, itb), ib = 1, nwbnd)
      enddo
      close(500)
      !
      ! Wd (wave direction)
      !
      open(500, file=trim(bwdfile))
      allocate(wd_bwv(nwbnd, ntwbnd))
      do itb = 1, ntwbnd
         read(500,*)t_bwv(itb),(wd_bwv(ib, itb), ib = 1, nwbnd)
      enddo
      close(500)
      !
      ! Ds (directional spreading)
      !
      open(500, file=trim(bdsfile))
      allocate(ds_bwv(nwbnd, ntwbnd))
      do itb = 1, ntwbnd
         read(500,*)t_bwv(itb),(ds_bwv(ib, itb), ib = 1, nwbnd)
      enddo
      close(500)
      !
      ! zs (water level) - TL: TODO CHECK STILL NEEDED?
      !if (trim(bzsfile) /= '') then 
      !   !
      !   open(500, file=trim(bzsfile))
      !   allocate(zs_bwv(nwbnd, ntwbnd))
      !   do itb = 1, ntwbnd
      !      read(500,*)t_bwv(itb),(zs_bwv(ib, itb), ib = 1, nwbnd)
      !   enddo
      !   close(500)
      !   !
      !else   
      !   !
      !   allocate(zs_bwv(nwbnd, ntwbnd))
      !   zs_bwv = 0.0         
      !   !
      !endif
      !
   endif
   !
   end subroutine
   !
   !
   !
   subroutine find_boundary_indices()
   !
   use snapwave_data
   !
   implicit none
   !
   ! For each grid boundary point (kcs=2) :
   !
   ! Determine indices and weights of boundary points
   ! For tide and surge, these are the indices and weights of the points in the bnd file
   ! For waves, these are the indices and weights of the points in the cst file
   !
   integer k, m, n, ib1, ib2, ib, ic
   !
   real xgb, ygb, dst1, dst2, dst
   !
   if (nwbnd>0 .and. nb>0) then  !nbnd
      !
      ! Count number of boundary points
      !
      ! Allocate boundary arrays
      !
      ! Water levels at boundary
      !
      ! Wave arrays
      !
      allocate(ind1_bwv_cst(nb))
      allocate(ind2_bwv_cst(nb))
      allocate(fac_bwv_cst(nb))
      !
      ! Find two closest boundary condition points for each boundary point
      ! And the two closest coastline points
      !
      nb = 0
      !
      ! Loop through all grid points
      !
      do k = 1, no_nodes
         !         
         ! Check if this point is a boundary point
         !
         if (msk(k) == 2) then
            !
            nb = nb + 1
            !
            ! set boundary cell indices
            !
            nmindbnd(nb) = k 
            !
            xgb = x(k)
            ygb = y(k)
            !
            ! Indices and weights for wave boundaries
            !
            if (nwbnd>1) then
               !
               dst1 = 1.0e10
               dst2 = 1.0e10
               ib1 = 0
               ib2 = 0
               !
               ! Loop through all water level boundary points
               !
               do ic = 1, nwbnd
                  !
                  ! Compute distance of this point to grid boundary point
                  !
                  dst = sqrt((x_bwv(ic) - xgb)**2 + (y_bwv(ic) - ygb)**2)
                  !
                  if (dst<dst1) then
                     !
                     ! Nearest point found
                     !
                     dst2 = dst1
                     ib2  = ib1
                     dst1 = dst
                     ib1  = ic
                     !
                  elseif (dst<dst2) then
                     !
                     ! Second nearest point found
                     !
                     dst2 = dst
                     ib2  = ic
                     !
                  endif
               enddo
               !
               ind1_bwv_cst(nb)  = ib1
               ind2_bwv_cst(nb)  = ib2
               fac_bwv_cst(nb) = dst2/(dst1 + dst2)
               !
            else
               !
               ind1_bwv_cst(nb)  = 1
               ind2_bwv_cst(nb)  = 1
               fac_bwv_cst(nb)   = 1.0
               !
            endif
            !        
         endif
      enddo
   endif
   !   
   end subroutine   
   
subroutine find_nearest_depth_for_boundary_points()
    !
    ! Find nearest grid index in (xgb,ygb) for every boundary input point (x_bwv(ib), y_bwv(ib))
    ! Output is: deptht_bwv    
    !
    use snapwave_data
    !
    implicit none
    !
    real*4  :: h1, h2, fac
    !
    real xgb, ygb, dst1, dst2, dst   
    integer k, ib1, ib2, ic
    !
    ! Loop through all water level boundary points
    !
    do ic = 1, nwbnd    
        ! Loop through all grid points
        !
	    dst1 = 1.0e10
	    dst2 = 1.0e10
	    ib1 = 0
	    ib2 = 0
        !    
        do k = 1, no_nodes
            !
	        xgb = x(k)
	        ygb = y(k)          
	        !
	        dst = sqrt((x_bwv(ic) - xgb)**2 + (y_bwv(ic) - ygb)**2)
	        !
	        if (dst<dst1) then
		        !
		        ! Nearest point found
		        !
		        dst2 = dst1
		        ib2  = ib1
		        dst1 = dst
		        ib1  = k
		        !
	        elseif (dst<dst2) then
		        !
		        ! Second nearest point found
		        !
		        dst2 = dst
		        ib2  = k
		        !
	        endif    
        enddo
        !
        if (ib2 == 0) then
            !
            write(*,*)'Warning: only 1 close grid point found for boundary input location (x,y): ',x_bwv(ic),y_bwv(ic)
            deptht_bwv(ic) = depth(ib1)
            !
        else
            !    
            h1  = depth(ib1)
            h2  = depth(ib2)
            fac = dst2/(dst1 + dst2)   
            deptht_bwv(ic) = h1*fac + h2*(1.0 - fac)
            !
        endif        
    enddo    
    !
end subroutine

subroutine update_boundary_conditions(t)
   !
   ! Update all wave boundary conditions
   !
   use snapwave_data
   !
   implicit none
   !
   real*8           :: t
   !
   ! Update boundary conditions at boundary points
   !
   call update_boundary_points(t)
   !
   ! Update wind forcing 
   !
   if (wind) then
      call update_wind_field()
   endif
   !
   ! Make directional grid around boundary mean wave/wind direction
   !
   thetamean = wdmean_bwv
   if (ntwbnd > 0) then
      call make_theta_grid(wdmean_bwv)
   else
      thetamean=u10dmean
      call make_theta_grid(u10dmean)
   endif 
   !
   ! Build spectra on the boundary support points
   !
   !call build_boundary_support_points_spectra() TODO - TL: later can clean up this code by also using this function
   !
   ! Update boundary conditions at grid points
   !
   call update_boundaries()
   !
end subroutine   
   
subroutine update_boundary_points(t)
   !
   ! Update vardens at boundary points
   !
   use snapwave_data
   use snapwave_infragravity
   !
   implicit none
   !
   integer i, k, ib, itb, itb0, itb1, itsp2now, itheta, i1, i2, ii, jj, ind
   !
   real*8  :: t, tb
   !
   real*4  :: tbfac
   real*4  :: hs, tps, wd, dsp, zst, thetamin, thetamax, E0, ms, modth, E0_ig
   real*4  :: hlocal
   logical :: always_update_bnd_spec
   !
   always_update_bnd_spec = .true.
   !
   ! Update directional spectra on boundary points from time series
   !
   ! First interpolate boundary conditions in timeseries to boundary points
   !
   update_grid_boundary_points = .true.
   !
   ! Interpolate boundary conditions in time
   !
   if (t_bwv(1) > (t - 1.0e-3)) then ! use first time in boundary conditions
      !
      itb0 = 1
      itb1 = 1
      tb   = t_bwv(itb0)
      !
   elseif (t_bwv(ntwbnd) < (t + 1.0e-3)) then  ! use last time in boundary conditions       
      !
      itb0 = ntwbnd
      itb1 = ntwbnd
      tb   = t_bwv(itb0)
      !
   else
      !
      do itb = itwbndlast, ntwbnd ! Loop in time
         if (t_bwv(itb) > (t + 1.0e-6)) then
            itb0 = itb - 1
            itb1 = itb
            tb   = t
            itwbndlast = itb - 1
            exit
         endif
      enddo 
      !
   endif            
   !
   tbfac  = (tb - t_bwv(itb0))/max(t_bwv(itb1) - t_bwv(itb0), 1.0e-6)
   !
   do ib = 1, nwbnd ! Loop along boundary points
      !
      hs  = hs_bwv(ib, itb0) + (hs_bwv(ib, itb1) - hs_bwv(ib, itb0))*tbfac
      tps = tp_bwv(ib, itb0) + (tp_bwv(ib, itb1) - tp_bwv(ib, itb0))*tbfac
      dsp = ds_bwv(ib, itb0) + (ds_bwv(ib, itb1) - ds_bwv(ib, itb0))*tbfac    !dirspr
      !zst = zs_bwv(ib, itb0) + (zs_bwv(ib, itb1) - zs_bwv(ib, itb0))*tbfac
      !
      ! Limit wave height (0 < hs < 25)
      !       
      if ((hs < 0.0) .or. (hs > 25.0)) then
      	  write(*,*)'DEBUG SnapWave - input wave height is outside acceptable range of 0-25 m: ',hs, ' and is therefore limited back to this range, please check whether your input is realistic!'
          hs = max(min(hs, 25.0), 0.0)
      endif  
      !
      ! Limit directional spreading (1 < ds < 60)
      !       
      if ((dsp < 3*pi/180) .or. (dsp > 60*pi/180)) then  
      	  write(*,*)'DEBUG SnapWave - input wave spreading is outside acceptable range of 3-60 degrees: ',dsp/pi*180, ' and is therefore limited back to this range, please check whether your input is realistic!'          
          dsp = max(min(dsp, 60*pi/180), 3*pi/180)
      endif      
      !
      ! Limit period (0.1 < tps < 25)
      !       
      if ((tps < 0.1) .or. (tps > 25.0)) then
      	  write(*,*)'DEBUG SnapWave - input wave period is outside acceptable range of 0.1-25 s: ',tps, ' and is therefore limited back to this range, please check whether your input is realistic!'
          tps = max(min(tps, 25.0), 0.1)
      endif      
      !
      call weighted_average(wd_bwv(ib, itb0), wd_bwv(ib, itb1), 1.0 - tbfac, 2, wd)  !wavdir
      !
      hst_bwv(ib) = hs
      tpt_bwv(ib) = tps                  
      wdt_bwv(ib) = wd
      dst_bwv(ib) = dsp
      !zst_bwv(ib) = zst
      !
   enddo
   !
!   do itb = itwbndlast, ntwbnd ! Loop in time
!      !
!      if (t_bwv(itb)>t) then
!         !
!         tbfac  = (t - t_bwv(itb - 1))/(t_bwv(itb) - t_bwv(itb - 1))
!         !
!         do ib = 1, nwbnd ! Loop along boundary points
!            !
!            hs    = hs_bwv(ib, itb - 1) + (hs_bwv(ib, itb) - hs_bwv(ib, itb - 1))*tbfac
!            tps   = tp_bwv(ib, itb - 1) + (tp_bwv(ib, itb) - tp_bwv(ib, itb - 1))*tbfac
!            dsp    = ds_bwv(ib, itb - 1) + (ds_bwv(ib, itb) - ds_bwv(ib, itb - 1))*tbfac    !dirspr
!            zst    = zs_bwv(ib, itb - 1) + (zs_bwv(ib, itb) - zs_bwv(ib, itb - 1))*tbfac  
!            !
!            call weighted_average(wd_bwv(ib, itb - 1), wd_bwv(ib, itb), 1.0 - tbfac, 2, wd)  !wavdir
!            !
!            hst_bwv(ib) = hs
!            tpt_bwv(ib) = tps                  
!            wdt_bwv(ib) = wd
!            dst_bwv(ib) = dsp
!            zst_bwv(ib) = zst
!            !
!         enddo
!         !
!         itwbndlast = itb
!         exit
!         !
!      endif
!   enddo
   !
   ! Now generate wave spectra at the boundary points
   !
   ! Average wave period and direction to determine theta grid
   !
   tpmean_bwv = sum(tpt_bwv)/size(tpt_bwv)
   !zsmean_bwv = sum(zst_bwv)/size(zst_bwv)
   !depth      = max(zsmean_bwv - zb,hmin) ! TL: For SFINCS we don't want this, because it overrides the real updated water depth we're inserting
   depth      = max(depth,hmin)
   wdmean_bwv = atan2(sum(sin(wdt_bwv)*hst_bwv)/sum(hst_bwv), sum(cos(wdt_bwv)*hst_bwv)/sum(hst_bwv))
   !
   ! Determine IG boundary conditions
   !
   if (igwaves) then
      ! 
      if (igherbers) then 
         !
         ! Get local water depth at boundary points (can change in time)        
         call find_nearest_depth_for_boundary_points() ! Output is: deptht_bwv    
         !
         do ib = 1, nwbnd ! Loop along boundary points
            !           
            ! Determine IG wave height and period at boundary
            !           
            call determine_ig_bc(x_bwv(ib), y_bwv(ib), hst_bwv(ib), tpt_bwv(ib), dst_bwv(ib), jonswapgam, deptht_bwv(ib), Tinc2ig, tpig_opt, hst_bwv_ig(ib), tpt_bwv_ig(ib))
            !           
            ! input, input, input, input, input, input, input, output, output
            !
         enddo   
         !
         tpmean_bwv_ig = sum(tpt_bwv_ig)/size(tpt_bwv_ig)       
         !
         write(*,*)'Herbers computed: hst_bwv_ig= ',hst_bwv_ig
         write(*,*)'Herbers computed: tpmean_bwv_ig= ',tpmean_bwv_ig      
         !     
      else
         !
         tpmean_bwv_ig = tpmean_bwv * Tinc2ig !TL: the old way using Tp inc to IG ratio
         !
      endif
   endif  
   !
   ! Determine theta grid and adjust w, prev and ds tables
   !
   ! Definition of directional grid
   !
   thetamean = wdmean_bwv
   !
   ind = nint(thetamean / dtheta) + 1
   !
   do itheta = 1, ntheta
!      i360(itheta) = mod2(itheta + ind - 10, 36)
      i360(itheta) = mod2(itheta + ind - (1 + ntheta / 2), ntheta * 2)
   enddo
   !
   do itheta = 1, ntheta
      !
      theta(itheta) = theta360(i360(itheta))
      !
      do k = 1, no_nodes
         w(1, itheta, k)    = w360(1, i360(itheta), k)
         w(2, itheta, k)    = w360(2, i360(itheta), k)
         prev(1, itheta, k) = prev360(1, i360(itheta), k)
         prev(2, itheta, k) = prev360(2, i360(itheta), k)
         ds(itheta, k)      = ds360(i360(itheta), k)
      enddo
      !
   enddo   
   !
   ! Build spectra on wave boundary support points
   !
   !write(*,*)' thetamean = ',thetamean*180./pi
   !write(*,'(a,18f7.1)')'theta = ',theta*180./pi
   do ib = 1, nwbnd ! Loop along boundary points
      !
      E0   = 0.0625 * rho * g * hst_bwv(ib)**2
      ms   = 1.0 / dst_bwv(ib)**2 - 1.0
      dist = (cos(theta - thetamean))**ms
      !    
      eet_bwv(:,ib) = dist/sum(dist)*E0/dtheta
      !
   enddo
   !
   ! Build IG spectra on wave boundary support points   
   if (igwaves) then   
      if (igherbers) then 
          do ib = 1, nwbnd ! Loop along boundary points    
             !          
             E0_ig   = 0.0625 * rho * g * hst_bwv_ig(ib)**2
             ms   = 1.0 / dst_bwv(ib)**2 - 1.0
             dist = (cos(theta - thetamean))**ms      
             !         
             eet_bwv_ig(:,ib) = dist / sum(dist) * E0_ig / dtheta          
             !      
          enddo
      endif
   endif
   !         
end subroutine update_boundary_points
!
subroutine update_wind_field()
   !
   ! Update wind field at all grid cells
   !
   use snapwave_data
   !
   implicit none
   !
   integer k
   real*4, dimension(:), allocatable :: windspread360k
   !
   allocate(windspread360k(ntheta360))
   windspread360k=0.0
   !
   ! Interpolate boundary conditions in timeseries to boundary points 
   ! TL: already done in SFINCS
   !
   ! average wind direction
   ! 
   u10dmean = atan2(sum(sin(u10dir)*u10),sum(cos(u10dir)*u10))  
   !
   ! Initialize the distribution array of wind input
   ! 
   do k = 1,no_nodes
      windspread360k = (cos(theta360-u10dir(k)))**2.0
      where(cos(theta360-u10dir(k))<0.0) windspread360k = 0.0
      windspread360k = (windspread360k/sum(windspread360k))/dtheta ! normalized and converted to input per rad
      windspread360(:,k) = windspread360k
   enddo
   !   
end subroutine update_wind_field
!
!    
subroutine make_theta_grid(central_theta)
   !
   ! make theta grid based on boundary mean wave direction
   !
   use snapwave_data
   !
   implicit none
   !
   real, intent(in) :: central_theta
   integer k, itheta, ind
   !
   ! Determine theta grid and adjust w, prev and ds tables
   !
   ! Definition of directional grid
   !   
   ind=nint(central_theta/dtheta)-ntheta/2;
   do itheta = 1, ntheta
      i360(itheta)=mod2(itheta+ind,ntheta360)
   enddo
   !
   do itheta = 1, ntheta
      !
      theta(itheta) = theta360(i360(itheta))
      !
      do k = 1, no_nodes
         w(1, itheta, k)    = w360(1, i360(itheta), k)
         w(2, itheta, k)    = w360(2, i360(itheta), k)
         prev(1, itheta, k) = prev360(1, i360(itheta), k)
         prev(2, itheta, k) = prev360(2, i360(itheta), k)
         ds(itheta, k)      = ds360(i360(itheta), k)
         !
         windspreadfac(itheta, k) = windspread360(i360(itheta), k)
      enddo
      !
   enddo  
   !
   if (wind==1) then
       !
       ! initialization of distribution array of wind input
       !
       do k = 1, no_nodes
          windspreadfac(:,k)=(cos(theta-u10dir(k)))**mwind
          where(cos(theta-u10dir(k))<0.0) windspreadfac(:,k)=0.0
          if (sum(windspreadfac(:,k))>0.0) then      
             windspreadfac(:,k) = (windspreadfac(:,k)/sum(windspreadfac(:,k)))/dtheta
          else
             windspreadfac(:,k)=0.0
          endif  
       enddo
       !
   endif
   !
end subroutine make_theta_grid
!
subroutine build_boundary_support_points_spectra()
   !
   ! Update directional spectra on boundary points from time series
   !
   use snapwave_data
   !
   implicit none
   !
   integer :: ib
   !
   real*4  :: E0, ms
   !  
   do ib = 1, nwbnd ! Loop along boundary points
      E0   = 0.0625*rho*g*hst_bwv(ib)**2
      ms = 1.0/dst_bwv(ib)**2-1.0
      dist = sign(1.0,cos(theta - thetamean))*abs(cos(theta - thetamean))**ms
      where (abs(mod(pi+theta - thetamean,2.0*pi)-pi)>0.999*pi/2.0) dist = 0.0
      eet_bwv(:,ib) = dist/sum(dist)*E0/dtheta
   enddo
   !
end subroutine build_boundary_support_points_spectra
!
subroutine update_boundaries()
   !
   ! Update values at boundary points
   !
   use snapwave_data
   !
   implicit none
   !
   integer ib, i, k
   !
   !      
   ! Set wave parameters in all boundary points on grid
   !
   ! Loop through grid boundary points
   ! Now for all grid boundary points do spatial interpolation of the wave spectra from boundary points at polygon
   !
   do ib = 1, nb
      !
      k = nmindbnd(ib)
      !   
      do i = 1, ntheta
         ! 
         ee(i,k) = eet_bwv(i,ind1_bwv_cst(ib))*fac_bwv_cst(ib)  + eet_bwv(i,ind2_bwv_cst(ib))*(1.0 - fac_bwv_cst(ib))
         !                 
      enddo
      !
   enddo
   !
   if (igwaves) then
      ! 
      if (igherbers) then 
         !
         do ib = 1, nb
            !
            k = nmindbnd(ib)       
            !
            do i = 1, ntheta
               !          
               ee_ig(i,k) = eet_bwv_ig(i,ind1_bwv_cst(ib))*fac_bwv_cst(ib)  + eet_bwv_ig(i,ind2_bwv_cst(ib))*(1.0 - fac_bwv_cst(ib)) 
               ! 
            enddo
            !
         enddo
      else !TL: the old way using eeinc2ig ratio times incident wave energy
         !
         do ib = 1, nb
            !
            k = nmindbnd(ib)       
            !          
            ee_ig(:, k)  = eeinc2ig*ee(:,k)
            !       
          enddo
      endif      
   endif          
   !
   end subroutine   
   
   
   subroutine weighted_average(val1,val2,fac,iopt,val3) ! as in sfincs_boundaries.f90
   !
   implicit none
   !
   integer, intent(in)    :: iopt
   real*4,  intent(in)    :: val1
   real*4,  intent(in)    :: val2
   real*4,  intent(in)    :: fac
   real*4,  intent(out)   :: val3
   !
   real*4                 :: u1
   real*4                 :: v1
   real*4                 :: u2
   real*4                 :: v2
   real*4                 :: u
   real*4                 :: v
   !
   if (iopt==1) then
      !
      ! Regular
      !
      val3 = val1*fac  + val2*(1.0 - fac)
      !
   else
      !
      ! Angles (input must be in radians!)
      !
      u1 = cos(val1)
      v1 = sin(val1)
      u2 = cos(val2)
      v2 = sin(val2)
      !
      u = u1*fac  + u2*(1.0 - fac)
      v = v1*fac  + v2*(1.0 - fac)
      !
      val3 = atan2(v, u)
      !
   endif
   !
   end subroutine
   
   function mod2 (a,b) result (c)
   integer a,b,c
   !
   c = mod(a,b)
   if (c==0) c = b
   if (c<0)  c = c + b
   
   end function
   
!   function mod2real (a,b) result (c)
!   real*4 :: a,b,c
!   !
!   c = amod(a,b)
!   if (c==0) c = b
!   if (c<0)  c = c + b
!   !   
!   end function  

end module