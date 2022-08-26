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
      write(*,*)'Input boundary points found: ',nwbnd
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
   real*4  :: hs, tps, wd, dsp, zst, thetamin, thetamax, E0, ms, modth
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
   depth      = max(zsmean_bwv - zb,hmin)
   wdmean_bwv = atan2(sum(sin(wdt_bwv)*hst_bwv)/sum(hst_bwv), sum(cos(wdt_bwv)*hst_bwv)/sum(hst_bwv))
   !
   ! Determine theta grid and adjust w, prev and ds tables
   !
   ! Definition of directional grid
   !
   thetamean = wdmean_bwv
   
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
      E0   = 0.0625*rho*g*hst_bwv(ib)**2  !QUESTION TL: why 1/16 instead of 1/8?
      ms   = 1.0/dst_bwv(ib)**2-1
      dist = (cos(theta - thetamean))**ms
!      where (abs(mod(pi + theta - thetamean, 2*pi) - pi)>pi/2) dist = 0.0
!      if (t> 405000.0) then
!         write(*,'(a,i8,40e14.4)')'E',ib,E0,ms,dist
!         endif
      do itheta = 1, ntheta
!         modth = mod2real(pi + theta(itheta) - thetamean, 2*pi)
!         if (abs(modth)>0.5*pi) then
!            dist(itheta) = 0.0
!         endif   
      enddo   
      eet_bwv(:,ib) = dist/sum(dist)*E0/dtheta
!      write(*,'(a,18f7.0)')' ee          ',eet_bwv(:,ib)
   enddo
   
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
