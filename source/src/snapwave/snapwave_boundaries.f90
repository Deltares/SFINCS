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
   ! Determine IG wave height and period at boundary point(s) later using Herbers 1994 
   ! (depends on local -changing- water depth, so not calculatable a priori)
   !
   if (igwaves) then
      ! 
      ! Allocate vars for IG: 
      allocate(hst_bwv_ig(nwbnd))
      allocate(tpt_bwv_ig(nwbnd))          
      allocate(eet_bwv_ig(ntheta,nwbnd)) 
      ! 'nwbnd' was determined in read_boundary_data_singlepoint/read_boundary_data_timeseries
      ! 
      ! Also needed:
      allocate(deptht_bwv(nwbnd))      
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
   write(*,*)'Reading wave boundary enclosure ...'
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
   wd_bwv = (270.0 - wd_bwv)*pi/180
   ds_bwv = ds_bwv*pi/180.

   !
   close(11)
   !
   nwbnd = 1
   !
   allocate(hst_bwv(nwbnd))
   allocate(tpt_bwv(nwbnd))
   allocate(wdt_bwv(nwbnd))
   allocate(dst_bwv(nwbnd))
   allocate(zst_bwv(nwbnd)) 
   allocate(eet_bwv(ntheta, nwbnd)) 
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
      allocate(hst_bwv(nwbnd))
      allocate(tpt_bwv(nwbnd))
      allocate(wdt_bwv(nwbnd))
      allocate(dst_bwv(nwbnd))
      allocate(zst_bwv(nwbnd)) 
      allocate(eet_bwv(ntheta,nwbnd)) 
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
      ! Convert to cartesian, going-to, radians
      wd_bwv = (270.0 - wd_bwv)*pi/180.
      !
      ! Ds (directional spreading)
      !
      open(500, file=trim(bdsfile))
      allocate(ds_bwv(nwbnd, ntwbnd))
      do itb = 1, ntwbnd
         read(500,*)t_bwv(itb),(ds_bwv(ib, itb), ib = 1, nwbnd)
      enddo
      close(500)
      ds_bwv = ds_bwv*pi/180.
      !
      ! zs (water level)
      if (trim(bzsfile) /= '') then 
         !
         open(500, file=trim(bzsfile))
         allocate(zs_bwv(nwbnd, ntwbnd))
         do itb = 1, ntwbnd
            read(500,*)t_bwv(itb),(zs_bwv(ib, itb), ib = 1, nwbnd)
         enddo
         close(500)
         !
      else   
         !
         allocate(zs_bwv(nwbnd, ntwbnd))
         zs_bwv = 0.0         
         !
      endif
      !
      write(*,*)'Input wave boundary points found: ',nwbnd
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
   ! Now do also the same for grid active point (kcs=1) in case of IG waves:
   ! Used for calculate e.g. relevant offshore wave steepness, representative waterdepth, incident wave breaking point, alphaig shoaling parameter depth cutoff
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
   ! Now do the same in case of IG waves, for normal active points (msk=1)
   if (igwaves) then   
       if (nwbnd>0) then  !nbnd
          !
          ! Count number of boundary points
          !
          ! Allocate active cell arrays
          !
          ! Water levels at boundary
          !
          ! Wave arrays
          !
          allocate(ind1_awv_cst(no_nodes))
          allocate(ind2_awv_cst(no_nodes))
          allocate(fac_awv_cst(no_nodes))
          !
          ! Find two closest boundary condition points for each active grid point
          ! And the two closest coastline points
          !
          ! Loop through all grid points
          !
          do k = 1, no_nodes
             !         
             ! Check if this point is a active grid point
             !
             if (msk(k) == 1) then
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
                   ind1_awv_cst(k)  = ib1
                   ind2_awv_cst(k)  = ib2
                   fac_awv_cst(k) = dst2/(dst1 + dst2)
                   !
                else
                   !
                   ind1_awv_cst(k)  = 1
                   ind2_awv_cst(k)  = 1
                   fac_awv_cst(k)   = 1.0
                   !
                endif
                !  
             else ! just to fill array, for now, set to 1 for not msk=1 points
                ind1_awv_cst(k)  = 1
                ind2_awv_cst(k)  = 1
                fac_awv_cst(k)   = 1.0                 
             endif
          enddo
       endif
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
!   write(*,*)'t=',t
   !
   call update_boundary_points(t)
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
   !
   implicit none
   !
   integer i, k, ib, itb, itsp2now, itheta, i1, i2, ii, jj, ind
   !
   real*8  :: t
   !
   real*4  :: tbfac
   real*4  :: hs, tps, wd, dsp, zst, thetamin, thetamax, E0, ms, modth, E0_ig
   real*4  :: jonswapgam, hlocal
   logical :: always_update_bnd_spec
   !
!   write(*,*)'dtheta',dtheta
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
   do itb = itwbndlast, ntwbnd ! Loop in time
!      write(*,*)itb,t_bwv(itb)
      !
      if (t_bwv(itb)>t) then
         !
         tbfac  = (t - t_bwv(itb - 1))/(t_bwv(itb) - t_bwv(itb - 1))
         !
         do ib = 1, nwbnd ! Loop along boundary points
            !
            hs    = hs_bwv(ib, itb - 1) + (hs_bwv(ib, itb) - hs_bwv(ib, itb - 1))*tbfac
            tps   = tp_bwv(ib, itb - 1) + (tp_bwv(ib, itb) - tp_bwv(ib, itb - 1))*tbfac
            dsp    = ds_bwv(ib, itb - 1) + (ds_bwv(ib, itb) - ds_bwv(ib, itb - 1))*tbfac    !dirspr
            zst    = zs_bwv(ib, itb - 1) + (zs_bwv(ib, itb) - zs_bwv(ib, itb - 1))*tbfac  
            !
            call weighted_average(wd_bwv(ib, itb - 1), wd_bwv(ib, itb), 1.0 - tbfac, 2, wd)  !wavdir
            !
            hst_bwv(ib) = hs
            tpt_bwv(ib) = tps                  
            wdt_bwv(ib) = wd
            dst_bwv(ib) = dsp
            zst_bwv(ib) = zst
            !
         enddo
         !
         itwbndlast = itb
         exit
         !
      endif
   enddo
   !
   ! Now generate wave spectra at the boundary points
   !
   ! Average wave period and direction to determine theta grid
   !
   tpmean_bwv = sum(tpt_bwv)/size(tpt_bwv)
   zsmean_bwv = sum(zst_bwv)/size(zst_bwv)
   !depth      = max(zsmean_bwv - zb,hmin) ! TL: For SFINCS we don't want this, because it overrides the real updated water depth we're inserting
   depth      = max(depth,hmin)
   wdmean_bwv = atan2(sum(sin(wdt_bwv)*hst_bwv)/sum(hst_bwv), sum(cos(wdt_bwv)*hst_bwv)/sum(hst_bwv))
   !
   ! Determine IG boundary conditions
   !
   if (igwaves) then
      ! 
      jonswapgam = 3.3 ! TODO: TL: later make spatially varying? > then as gam_bwv(ib) in 'determine_ig_bc'
      !jonswapgam = 20.0
      !
      ! Get local water depth at boundary points (can change in time)        
      call find_nearest_depth_for_boundary_points() ! Output is: deptht_bwv    
      !
      do ib = 1, nwbnd ! Loop along boundary points
         !           
         ! Determine IG wave height and period at boundary
         call determine_ig_bc(hst_bwv(ib), tpt_bwv(ib), dst_bwv(ib), jonswapgam, deptht_bwv(ib), Tinc2ig, hst_bwv_ig(ib), tpt_bwv_ig(ib)) 
         ! input, input, input, input, input, input, output, output
         !
      enddo   
      !
      tpmean_bwv_ig = sum(tpt_bwv_ig)/size(tpt_bwv_ig)       
      !
      write(*,*)'hst_bwv_ig= ',hst_bwv_ig
      write(*,*)'tpmean_bwv_ig= ',tpmean_bwv_ig      
      !      
   endif  
   !
   ! Determine theta grid and adjust w, prev and ds tables
   !
   ! Definition of directional grid
   !
   thetamean = wdmean_bwv
   !
   ind = nint(thetamean/dtheta) + 1
   do itheta = 1, ntheta
!      i360(itheta) = mod2(itheta + ind - 10, 36)
      i360(itheta) = mod2(itheta + ind - (1+ntheta/2), ntheta*2)
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
      E0   = 0.0625*rho*g*hst_bwv(ib)**2
      ms   = 1.0/dst_bwv(ib)**2-1
      dist = (cos(theta - thetamean))**ms
!      where (abs(mod(pi + theta - thetamean, 2*pi) - pi)>pi/2) dist = 0.0
!      if (t> 405000.0) then
!         write(*,'(a,i8,40e14.4)')'E',ib,E0,ms,dist
!         endif
      !do itheta = 1, ntheta
!         modth = mod2real(pi + theta(itheta) - thetamean, 2*pi)
!         if (abs(modth)>0.5*pi) then
!            dist(itheta) = 0.0
!         endif   
      !enddo   
      eet_bwv(:,ib) = dist/sum(dist)*E0/dtheta
!      write(*,'(a,18f7.0)')' ee          ',eet_bwv(:,ib)
      !
   enddo
   !
   ! Build IG spectra on wave boundary support points   
   if (igwaves) then   
      do ib = 1, nwbnd ! Loop along boundary points    
         !          
         E0_ig   = 0.0625*rho*g*hst_bwv_ig(ib)**2
         ms   = 1.0/dst_bwv(ib)**2-1
         dist = (cos(theta - thetamean))**ms      
         !         
         eet_bwv_ig(:,ib) = dist/sum(dist)*E0_ig/dtheta          
         !      
      enddo
      !
   endif
   !         
end subroutine
   
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
   endif          
   !
   ! Set representative incident wave height in all active points on grid > TL: to be removed
   ! 
   if (igwaves) then
      do k = 1, no_nodes
         if (inner(k)) then
            !
            H_rep(k) = max(0.1, (hs_bwv(1,ind1_awv_cst(k))*fac_awv_cst(k)  + hs_bwv(1,ind2_awv_cst(k))*(1.0 - fac_awv_cst(k))))
            ! Don't convert from Hm0 to Hrms yet!
            ! Note; make sure that H_rep(k) is never 0! (problem in calculating steepness_bc, reldepth and depthforcerelease later)
            ! 
         endif         
      enddo               
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

    subroutine determine_ig_bc(hsinc, tpinc, ds, jonswapgam, depth, Tinc2ig, hsig, tpig) 
    ! (input, input, input, input, input, input, output, output)
    !  
    ! Build boundwave offshore spectrum, and determine Hig0 and Tpig0, using Herbers 1994 as in XBeach implementation
    !
    ! Input hsinc is Hm0, not Hrms
    ! Output hsig is therefore also Hm0 (and expected like that in E0_ig   = 0.0625*rho*g*hst_bwv_ig(ib)**2 in subroutine update_boundary_points)
    !
    ! tpig is estimated as ... based on calculated spectrum
    !
    ! Input ds is already in rad
    !
    !use snapwave_data
    !   
    implicit none
    !
    real*4, intent(in)    :: hsinc, tpinc, ds, jonswapgam, depth, Tinc2ig
    real*4, intent(out)   :: hsig, tpig
    !
    real*4                :: pi, scoeff
    real*4                :: Tm01, Tm10, Tp, Tpsmooth
    integer               :: correctHm0, tpigopt
    !
    correctHm0 = 0 ! Choice between correcting Hm0 in build_jonswap if build 2D Vardens spectrum too low (1) or not (default, 0)
    !
    pi    = 4.*atan(1.)
    !
    ! Convert wave spreading in degrees (input) to S
    scoeff = (2/ds**2) - 1
    !  
    ! Call function that calculates Hig0 following Herbers, as also implemented in XBeach and secordspec2 in Matlab
    ! Loosely based on 3 step calculation in waveparams.F90 of XBeach (build_jonswap, build_etdir, build_boundw), here all in 1 subroutine calculate_herbers
    !
    call compute_herbers(hsig, Tm01, Tm10, Tp, Tpsmooth, hsinc, tpinc, scoeff, jonswapgam, depth, correctHm0) ![out,out,out,out,out, in,in,in,in,in,in]
    !   
    ! Catch NaN values (if depth=0 probably) or unrealistically large values above 2 meters
    if (hsig < 0.0) then
	    write(*,*)'DEBUG - computed hm0ig at boundary dropped below 0 m: ',hsig, ' and is therefore limited back to 0 m!'
	    hsig = max(hsig, 0.0)
    endif	
    if (hsig > 3.0) then
	    write(*,*)'DEBUG - computed hm0ig at boundary exceeds 2 meter: ',hsig, ' and is therefore limited back to 2 m!'
	    hsig = min(hsig, 3.0)
    endif	        
    !
    ! Choose what wave period option value for IG to choose:
    !
    tpigopt = 1 !TODO: TL: later make user definable
    !
    if (tpigopt == 1) then
        !
        tpig = Tm01
        !
    elseif (tpigopt == 2) then
        !
        tpig = tpinc * Tinc2ig
        !
    elseif (tpigopt == 3) then
        !        
        tpig = Tpsmooth
        !
    elseif (tpigopt == 4) then
        !        
        tpig = Tp            
        !
    elseif (tpigopt == 5) then
        !        
        tpig = Tm10 ! Tm-1,0
        !    
    elseif (tpigopt == 6) then ! Old default ratio of Tpig = 7 * Tpinc
        !        
        tpig = tpinc * 7.0
        !
    endif
    !
    ! Check on ratio tpig/tpinc whether it is deemed realistic
    if (tpig/tpinc < 2.0) then
	    write(*,*)'DEBUG - computed tpig/tpinc ratio at offshore boundary dropped below 2 and might be unrealistic! value: ',tpig/tpinc
    endif	    
    if (tpig/tpinc > 25.0) then
	    write(*,*)'DEBUG - computed tpig/tpinc ratio at offshore boundary increased above 25 and might be unrealistic! value: ',tpig/tpinc
    endif	         
    !   
    end subroutine
   
    subroutine compute_herbers(hsig, Tm01, Tm10, Tp, Tpsmooth, hsinc, tpinc, scoeff, jonswapgam, depth, correctHm0)
    !
    use interp        
    !
    implicit none
    !
    ! Incoming and outgoing variables
    real*4, intent(in)                      :: hsinc, tpinc, scoeff, jonswapgam, depth
    integer, intent(in)                     :: correctHm0
    real*4, intent(out)                     :: hsig, Tm01, Tm10, Tp, Tpsmooth    
    !
    ! Internal variables - for part 1: build_jonswap
    real*4                                  :: dfj, fp, fnyq, dang, iang, angtemp, hsinc_check1, df
    real*4, dimension(:), allocatable       :: temp, x, y, f, ang, Dd, Sf, findline, fgen    
    real*4, dimension(:,:), allocatable     :: S_array    
    integer                                 :: i=0, ii, nang, nfreq, peakf
    integer                                 :: firstp, lastp, M, K    
    real*4                                  :: pi, g, sprdthr
    !
    ! Internal variables - for part 2: build_etdir    
    real*4                                  :: kmax, pp
    real*4, dimension(400)                  :: P0, kk, phase, Sf0, Sf0org,S0org     !K=400 - OR allocate size later once K is defined
    real*4, dimension(202)                  :: P1, ang1   ! Are size 200+2        
    real*4, dimension(:), allocatable       :: theta, theta0, temp2, S0, Sn, dthetafin
    real*4, dimension(:,:), allocatable     :: Snew_array, Snew1_array, Ddnew        
    real*4                                  :: hm0now, s1, s2, modf, modang, hsinc_check2, hsinc_check2tmp    
    integer                                 :: F2, stepf, stepang, jj
    !
    ! Internal variables - for part 3: build_boundw    
    real*4                                  :: deltaf
    real*4, dimension(:), allocatable       :: w1, k1
    real*4, dimension(:), allocatable       :: Ebnd, fbnd
    real*4, dimension(:), allocatable       :: term1, term2, term2new, dif, chk1, chk2
    real*4, dimension(:,:), allocatable     :: Eforc, D, deltheta!, KKx, KKy, theta3
    real*4, dimension(:,:), allocatable     :: dphi3, k3, cg3!, Abnd   
    integer                                 :: valdensmax
    !
    ! Constants
    pi  = 4.*atan(1.)    
    g   = 9.81
    !   
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! TL: Here we go to the code of XBeach subroutine 'build_jonswap' from waveparams.F90 
    !
    dfj = 0.001
    fp = 1 / tpinc
    fnyq = fp * 2
    !
    ! Define number of frequency bins by defining an array of the necessary length using the Nyquist frequency and frequency step size
    allocate(temp(ceiling((fnyq-dfj)/dfj)))
    temp=(/(i,i=1,size(temp))/)   
    !
    ! Define array with actual equidistant frequency bins
    allocate(f(size(temp)))
    f=temp*dfj
    deallocate (temp)
    !
    ! Determine frequency bins relative to peak frequency
    allocate(x(size(f)))
    x=f/fp
    !
    ! Calculate unscaled and non-directional JONSWAP spectrum using peak-enhancement factor and pre-determined frequency bins
    allocate(y(size(f)))
    call jonswapgk(x,jonswapgam,y)
    !
    ! Determine scaled and non-directional JONSWAP spectrum using the JONSWAP characteristics
    y = (hsinc/(4.d0*sqrt(sum(y)*dfj)))**2*y
    deallocate (x)  
    ! 
    ! Define 200 directions relative to main angle running from -pi to pi
  
    allocate(temp(200))
    allocate(ang(200))    
    temp=(/(i,i=0,200-1)/) ! as in Matlab version have 200 values    
    ang=temp*((2*pi)/200.d0)-pi
    deallocate (temp)    
    !
    ! Determine directional step size: pi/200
    dang=ang(2)-ang(1)
    !
    ! Define 200 directional bins accordingly
    allocate (Dd(size(ang)))    
    !
    ! Calculate directional spreading based on cosine law    
    Dd = cos(ang/2)**(2*nint(scoeff) ) ! Robert: apparently nint is needed here, else MATH domain error TL: also the case for us
    !
    ! Scale directional spreading to have a surface of unity by dividing by it's own surface
    Dd = Dd / (sum(Dd)*dang)    
    !
    ! Define number of directional and frequency bins
    nang=size(ang)
    nfreq=size(y)        
    !
    ! Define two-dimensional variance density spectrum array and distribute variance density for each frequency over directional bins
    allocate(S_array(nfreq,nang))
    !
    do i=1,nang
        do ii=1,nfreq
            S_array(ii,i)=y(ii)*Dd(i)
        end do
    end do
    deallocate (y)    
    !
    ! Add check on Hm0 of created two-dimensional variance density spectrum :
    hsinc_check1 = 4*sqrt(sum(S_array)*dfj*dang)
    if (abs(hsinc - hsinc_check1) > 0.01) then
        write(*,*)'WARNING - computed Hm0,inc of 2D var dens spectrum differs from input! see subroutine build_jonswap in determine_ig_bc in module snapwave_boundaries.f90'
        write(*,*)'Input: ',hsinc,' , computed: ', hsinc_check1
    endif
    !
    ! Back integrate two-dimensional variance density spectrum over directions    
    allocate(Sf(size(f)))
    Sf = sum(S_array, DIM = 2)*dang    
    !
    ! Determine frequencies around peak frequency of one-dimensional non-directional variance density spectrum, based on factor sprdthr (0.08 by default)
    allocate(findline(size(Sf)))
    findline = 0.d0    
    sprdthr = 0.08
    call frange(Sf,firstp,lastp,findline,sprdthr)
    !findline = findloc(Sf, Sf>0.08*maxval(Sf)) ! TL: should probably be possible to do this more easily
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! TL: Here we go to the code of XBeach subroutine 'build_etdir' from waveparams.F90 
    !
    K = 400 ! as in Matlab script
    !
    ! Determine number of frequencies in discrete variance density spectrum to be included    
    M = int(sum(findline)) ! number of points in frequency range around peak
    !
    ! Define number of wave components to be used
    allocate (temp(K))
    temp=(/(i,i=0,K-1)/)
    !  
    ! Select equidistant wave components between the earlier selected range of frequencies around the peak frequency based on sprdthr
    allocate(fgen(K))
    fgen = temp * ((f(lastp) - f(firstp)) / (K-1)) + f(firstp)
    deallocate(temp)
    !
    ! Determine equidistant frequency step size
    df = (fgen(K) - fgen(1)) / (dble(K) - 1.d0)    
    !        
    ! 2D interpolate from S_array of (f,ang) grid to (fgen,ang)
    !
    allocate(Snew_array(K,size(ang)))    
    do ii=1,K
        do jj=1,size(ang)
            !
            call linear_interp_2d_real4(f,size(f), &    ! input frequency, size
                ang,size(ang), &                        ! input angles, size
                S_array, &                              ! input 2D variance density
                fgen(ii),ang(jj), &                     ! output frequency/angle
                Snew_array(ii,jj), &                    ! output variance density
                'interp',0.0)                           ! method and exception return value    
        enddo        
    enddo
    !
    Sn = sum(Snew_array, DIM = 2) * dang ! TL: similar to Sf    
    !
    allocate(Ddnew(K,size(ang)))    
    !
    do ii=1,K
        Ddnew(ii,:) = Snew_array(ii,:) / Sn(ii)
    enddo
    !
    ! Define random number between 0 and 1 for each wave component
    do i=1,K
        call RANDOM_NUMBER(P0(i)) 
        P0(i)=0.99*P0(i)+0.01/2 ! TL: as in Matlab: Define random number between ~0.025 and 0975 for each wave component
    enddo
    !
    ! Define direction for each wave component based on random number and linear interpolation of the probability density function
    allocate(theta(K))
    !
    if (scoeff >= 1000.0) then !longcrested waves
        theta = 0 
    else
        do ii=1,size(P0)
            !
            P1(1) = 0                   ! P1 is a function between 0 and 1, that looks like a cdf function
            ang1(1) = minval(ang)-dang  ! ang1 is a linear function between -pi and pi
            !
            do jj=1,size(ang)
                P1(jj+1) = sum(Ddnew(ii,1:jj)) * dang + 0.00001 * jj
                ang1(jj+1) = ang(jj)
            enddo            
            !
            P1(size(ang)+2) = P1(size(ang)+1) + 0.0001
            ang1(size(ang)+2) = ang(size(ang))+dang            
            ! TL: equivalent to Matlab implementation: P1 = [0 P P(length(ang))+0.0001] & ang1 = [min(ang)-dang ang max(ang)+dang];
            !                        
            ! Find the angle for each frequency component            
            call linear_interp_real4(P1,ang1,size(P1),P0(ii),theta(ii),F2)  ! theta(ii) is the output, finding the corresponding angle per random P0 value
            !
            ! TL: as in Matlab - for interpolation purposes:
            if (theta(ii) < minval(ang1)) then
                theta(ii) = theta(ii) + 2*pi
            endif
            if (theta(ii) > maxval(ang1)) then
                theta(ii) = theta(ii) - 2*pi
            endif            
            !
        enddo
    endif    
    !    
    ! TL: equivalent to Matlab - Snew = [S(:,noang) S S(:,1)]; % this is for interpolation purposes
    allocate(Snew1_array(K,size(ang)+2))    
    !
    do ii=1,K ! TL: =size(fgen)
        !
        Snew1_array(ii,1) = Snew_array(ii,size(ang))  
        Snew1_array(ii,size(ang)+2) = Snew_array(ii,1)
        !
        do jj=1,size(ang)
            !
            Snew1_array(ii,jj+1) = Snew_array(ii,jj)            
            !
        enddo        
    enddo
    !
    ! Determine variance density spectrum values for all relevant wave components
    ! around the peak frequency by interpolation of two-dimensional variance density spectrum array.    
    !
    ! 2D interpolate from S_array of (f,ang) grid to (fgen,ang)    
    allocate(S0(K))     
    !        
    do ii=1,K ! TL: =size(fgen)    
        !
        call linear_interp_2d_real4(fgen,size(fgen), &  ! input frequency, size
            ang1,size(ang1), &                          ! input angles, size
            Snew1_array, &                              ! input 2D variance density
            fgen(ii),theta(ii), &                       ! output frequency/angle
            S0(ii), &                                   ! output variance density
            'interp',0.0)                               ! method and exception return value    
        !
    enddo    
    !
    allocate(dthetafin(K))         
    dthetafin = Sn/S0					
    !
    ! Determine significant wave height using Hm0 = 4*sqrt(m0) using the one-dimensional non-directional variance density spectrum    
    hsinc_check2 = 4*sqrt(sum(S0 * dthetafin * df))
    hsinc_check2tmp = 4*sqrt(sum(Sn) * df)    
    !
    ! Add check on Hm0 of created two-dimensional variance density spectrum :
    !if (abs(hsinc_check2 - hsinc_check1) > 0.01) then
    if (abs(hsinc_check2 - hsinc_check2tmp) > 0.01) then        
        write(*,*)'WARNING - computed Hm0,inc of 2D var dens spectrum differs from input! see subroutine build_jonswap in determine_ig_bc in module snapwave_boundaries.f90'
        write(*,*)'Newly computed in part 2: ',hsinc_check2,' , while computed before: ', hsinc_check2tmp, ' and input was: ', hsinc
        !write(*,*)'Newly computed in part 2: ',hsinc_check2,' , computed in part 1: ', hsinc_check1        
    endif    
    !
    ! Correct spectra for wave height > TL: in Xbeach, not Matlab, option added for now
    if (correctHm0 == 1) then
        !
        S0 = (hsinc/hsinc_check2)**2 * S0
        Sn = (hsinc/hsinc_check2)**2 * Sn
        dthetafin = Sn/S0					    
        !
        ! To check again:
        hsinc_check2 = 4*sqrt(sum(S0 * dthetafin * df))
        hsinc_check2tmp = 4*sqrt(sum(Sn) * df)                 
        !
        write(*,*)'DEBUG - Hm0 of vardens corrected to: ',hsinc_check2,' and ', hsinc_check2tmp, ' so it is close to input: ', hsinc
        !        
    endif    
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! TL: Here we go to the code of XBeach subroutine 'build_boundw' from waveparams.F90 
    !
    ! Allocate two-dimensional variables for all combinations of interacting wave components to be filled triangular
    allocate(Eforc(K-1,K))      ! Herbers et al. (1994) eq. 1 - (rows = difference frequency, columns is interaction)
    allocate(D(K-1,K))          ! Herbers eq. A5.
    allocate(deltheta(K-1,K))   ! Difference angle between two primary wave components
    allocate(k3(K-1,K))         ! Wavenumber of difference wave
    allocate(cg3(K-1,K))
    !
    ! Allocate variables for angular velocity and wave numbers for wave components
    allocate(w1(size(fgen)))    ! Radial frequency of primary waves
    allocate(k1(size(fgen)))    ! Wave numbers of primary waves
    !
    ! Initialize variables as zero
    Eforc = 0
    D = 0
    deltheta = 0
    k3 = 0
    cg3 = 0
    w1=0
    k1=0    
    !
    ! Steps indepent of loop, that can be done once a priori
    ! Determine angular velocity of primary waves
    w1=2*pi*fgen
    !
    ! Determine wave numbers of primary waves
    call bc_disper(k1,w1,size(w1),depth,g)    
    !
    ! Determine for each wave component interactions with all other wave components
    ! as far as not processed yet (each loop step the number of interactions thus decrease with one)
    !
    do m=1,K-1
        !
        ! Determine difference frequency
        deltaf=m*df
        !
        ! Determine difference angles (pi already added)
        deltheta(m,1:K-m) = abs(theta(m+1:K)-theta(1:K-m))+pi
        !
        ! Determine difference wave numbers according to Van Dongeren et al. 2003 eq. 19
        k3(m,1:K-m) =sqrt(k1(1:K-m)**2+k1(m+1:K)**2+2*k1(1:K-m)*k1(m+1:K)*cos(deltheta(m,1:K-m))) 
        ! Note, dcos is for double precision, use cos for single precision (real*4)
        !
        ! Determine group velocity of difference waves
        cg3(m,1:K-m) = 2.d0*pi*deltaf/k3(m,1:K-m)
        !
        ! XBeach: Make sure that we don't blow up bound long wave when offshore boundary is too close to shore > not needed for us?
        !cg3(m,1:K-m) = min(cg3(m,1:K-m),par%nmax*sqrt(g/k3(m,1:K-m)*tanh(k3(m,1:K-m)*depth)))
        !
        ! Determine difference-interaction coefficient according to Herbers 1994 eq. A5
        allocate(term1(K-m),term2(K-m),term2new(K-m),dif(K-m),chk1(K-m),chk2(K-m))
        !
        term1 = (-w1(1:K-m))*w1(m+1:K)
        term2 = (-w1(1:K-m))+w1(m+1:K)
        term2new = cg3(m,1:K-m)*k3(m,1:K-m)
        dif = (abs(term2-term2new))
        !
        chk1  = cosh(k1(1:K-m)*depth)
        chk2  = cosh(k1(m+1:K)*depth)
        !
        D(m,1:K-m) = -g*k1(1:K-m)*k1(m+1:K)*cos(deltheta(m,1:K-m))/2.d0/term1+g*term2*(chk1*chk2)/ &
        ((g*k3(m,1:K-m)*tanh(k3(m,1:K-m)*depth)-(term2new)**2)*term1*cosh(k3(m,1:K-m)*depth))* &
        (term2*((term1)**2/g/g - k1(1:K-m)*k1(m+1:K)*cos(deltheta(m,1:K-m))) &
        - 0.50d0*((-w1(1:K-m))*k1(m+1:K)**2/(chk2**2)+w1(m+1:K)*k1(1:K-m)**2/(chk1**2)))
        !
        deallocate(term1,term2,term2new,dif,chk1,chk2)
        !
        ! Correct for surface elevation input and output instead of bottom pressure so it is consistent with Van Dongeren et al 2003 eq. 18
        D(m,1:K-m) = D(m,1:K-m)*cosh(k3(m,1:K-m)*depth)/(cosh(k1(1:K-m)*depth)*cosh(k1(m+1:K)*depth))
        !
        ! Exclude interactions with components smaller than or equal to current component according to lower limit Herbers 1994 eq. 1
        where(fgen<=m*df)
        D(m,:)=0.d0                                           ! Bas: redundant with initial determination of D ??
        endwhere
        !
        ! Determine energy of bound long wave according to Herbers 1994 eq. 1 based
        ! on difference-interaction coefficient and energy density spectra of
        ! primary waves
        Eforc(m,1:K-m) = 2*D(m,1:K-m)**2*S0(1:K-m)*S0(m+1:K)*dthetafin(1:K-m)*dthetafin(m+1:K)*df
        !
    end do
    !
    ! Allocate variables for energy of bound long wave
    allocate(Ebnd(K-1))
    !
    ! Sum over the components to get total forced wave at diff freq
    Ebnd = sum(Eforc,2)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Determine final parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    !
    ! What we actually want for SnapWave offshore IG bc: Hm0ig 
    hsig = 4*sqrt(sum(Ebnd)*df)   
    !
    ! Calculate representative value for IG wave period
    allocate (temp(K))
    temp=(/(i,i=0,K-1)/)
    !  
    ! Select equidistant wave components between the earlier selected range of frequencies around the peak frequency based on sprdthr
    allocate(fbnd(K-1))
    fbnd = temp * df
    deallocate(temp)
    !
    ! Calculate (mean) wave period based on one-dimensional non-directional variance density spectrum and factor trepfac
    call tpDcalc(Ebnd,fbnd,0.01,Tm01,Tm10,Tp,Tpsmooth)! [in,in,in,out,out,out,out]
    !XBeach default for trepfac = 0.01
    !
    end subroutine
            
    ! -----------------------------------------------------------
    ! --------- JONSWAP  unscaled JONSWAP spectrum --------------
    ! -----------------(used by compute_herbers)-----------------
    subroutine jonswapgk(x,gam,y)

    implicit none
    ! Required input: - x           : nondimensional frequency, divided by the peak frequency
    !                 - gam         : peak enhancement factor, optional parameter (DEFAULT 3.3)
    !                 - y is output : nondimensional relative spectral density, equal to one at the peak

    real*4, INTENT(IN)                  :: gam
    real*4,dimension(:), INTENT(IN)     :: x
    real*4,dimension(:), INTENT(INOUT)  :: y

    ! Internal variables
    real*4,dimension(size(x))           :: xa, sigma, fac1, fac2, fac3, temp

    xa=abs(x)

    where (xa==0)
        xa=1e-20
    end where

    sigma=xa

    where (sigma<1.)
        sigma=0.07
    end where

    where (sigma>=1.)
        sigma=0.09
    end where

    temp=0*xa+1

    fac1=xa**(-5)
    fac2=exp(-1.25*(xa**(-4)))
    fac3=(gam*temp)**(exp(-((xa-1)**2)/(2*(sigma**2))))

    y=fac1*fac2*fac3
    y=y/maxval(y)

    return

    end subroutine jonswapgk   
    
   ! -----------------------------------------------------------
   ! ---- Small subroutine to determine f-range round peak -----
   ! -----------------(used by compute_herbers)-----------------
   subroutine frange(Sf,firstp,lastp,findlineout,sprdthr)

      implicit none

      real*4, dimension(:), intent(in)        :: Sf
      integer, intent(out)                    :: firstp, lastp

      real*4, dimension(:), intent(out)       :: findlineout
      real*4, dimension(:),allocatable        :: temp, findline
      integer                                 :: i = 0
      real*4,intent(in)                       :: sprdthr
      
      allocate(findline(size(Sf)))
      findline=0*Sf                           ! find frequency range around peak

      where (Sf>sprdthr*maxval(Sf))
         findline=1
      end where

      firstp=maxval(maxloc(findline))         ! Picks the first "1" in temp

      allocate (temp(size(findline)))
      temp=(/(i,i=1,size(findline))/)
      lastp=maxval(maxloc(temp*findline))     ! Picks the last "1" in temp

      findlineout=findline
      deallocate(temp, findline)

   end subroutine frange    
    
   ! --------------------------------------------------------------
   ! --------------------- Dispersion relation --------------------
   ! -------------------(used by compute_herbers)------------------
   subroutine bc_disper(k1,w1,m,h,g)
      !          k  = wave number             (2 * pi / wave length)
      !          w  = wave angular frequency  (2 * pi / wave period)
      !          m  = size k and w vectors
      !          h  = water depth
      !          g  = gravitational acceleration constant, optional (DEFAULT 9.81)
      !
      !          absolute error in k*h < 5.0e-16 for all k*h
      !
      !
      !          original Matlab code by: G. Klopman, Delft Hydraulics, 6 Dec 1994

      integer, intent(in)                     :: m
      real*4,dimension(m),intent(in)          :: w1
      real*4,dimension(m),intent(out)         :: k1
      real*4, intent(in)                      :: h, g

      ! internal variables

      real*4,dimension(m)                     :: w2,q,thq,thq2,a,b,c,arg,sign
      integer                                 :: j
      real*4                                  :: hu

      w2 = w1**2*(h/g)
      q = w2/(1.0d0-exp(-(w2**(5.0d0/4.0d0))))**(2.0d0/5.0d0)

      do j=1,4
         thq  = tanh(q)
         thq2 = 1.0d0-thq**2
         a    = (1.0d0-q*thq)*thq2
         b    = thq + q*thq2
         c    = q*thq-w2
         where (abs(a*c)<(b**2*1.0e-8))
            arg = -c/b
         elsewhere
            arg  = (b**2)-4.0d0*a*c
            arg  = (-b + sqrt(arg))/(2.0d0*a)
         endwhere
         q    = q+arg
      end do

      where (w1>0.0d0)
         sign=1.0d0
      endwhere

      where (w1==0.0d0)
         sign=0.0d0
      endwhere

      where (w1<0.0d0)
         sign=-1.0d0
      endwhere

      k1 = sign*q/h

      where (k1==huge(hu))
         k1=0.0d0
      endwhere

      where (k1==-1.0d0*huge(hu))
         k1=0.0d0
      endwhere

      return

   end subroutine bc_disper   
   
   ! -----------------------------------------------------------
   ! ----------- Small subroutine to determine tpD -------------
   ! -----------------(used by compute_herbers)-----------------
   !
   ! TL: Slightly adapted from XBeach version to give more output,
   ! give back Tm01, Tm-1,0 and Tp
   subroutine tpDcalc(Sf,f,trepfac,Tm01,Tm10,Tp,Tpsmooth) ! [in,in,in,out,out,out,out]

      implicit none

      ! In:
      real*4, dimension(:), intent(in)        :: Sf, f
      real*4, intent(in)                      :: trepfac
      !
      ! Out:
      real*4, intent(out)                     :: Tm01
      real*4, intent(out)                     :: Tm10
      real*4, intent(out)                     :: Tp   
      real*4, intent(out)                     :: Tpsmooth            
      !
      ! Local:
      real*4, dimension(:),allocatable        :: temp, Sfsmooth
      integer                                 :: id  
      !
      allocate(temp(size(Sf)))
      temp=0.d0
      where (Sf>=trepfac*maxval(Sf))
         temp=1.d0
      end where
      !
      Tm01=sum(temp*Sf)/sum(temp*Sf*f)   ! Tm01
      !
      Tm10 = sum(temp*Sf/f)/sum(temp*Sf) ! Tm-1,0
      !
      id = int(MAXLOC(Sf, dim=1))
      Tp = 1.0 / f(id)                   ! Tp
      !    
      ! And Tp over smoothed spectrum:
      allocate(Sfsmooth(size(Sf)))
      !
      Sfsmooth = 0.0      
      do id=2,size(Sf)-1
         Sfsmooth(id) = (0.5*Sf(id-1)+Sf(id)+0.5*Sf(id+1))/2
      enddo               
      !
      id = int(MAXLOC(Sfsmooth, dim=1))
      Tpsmooth = 1.0 / f(id)                   ! Tp,smooth      
      !
   end subroutine tpDcalc
   
end module