module sfincs_boundaries
   
   use sfincs_log

contains

   subroutine read_boundary_data()
   !
   ! Reads bnd, bzs etc. files
   !
   use sfincs_ncinput
   use sfincs_data
   use sfincs_date
   use astro
   !
   implicit none
   !
   integer n, itb, ib, stat, ifreq, iok
   !
   real*4 dummy,r
   !
   ! Astro
   !
   character(8), dimension(:), allocatable :: tidal_component_names
   integer, dimension(6) :: i_date_time
   character(len=256)    :: line
   integer               :: n_sets, n_components, ios
   integer               :: current_set, current_component
   character(len=8)      :: cname
   real*4                :: a, p
   logical               :: has_a0
   !
   ! Read water level boundaries
   !
   nbnd      = 0
   ntbnd     = 0
   itbndlast = 1
   !
   if (netbndbzsbzifile(1:4) /= 'none') then    ! FEWS compatible Netcdf water level time-series input
      !
      call read_netcdf_boundary_data()
      !
      if ((t_bnd(1) > (t0 + 1.0)) .or. (t_bnd(ntbnd) < (t1 - 1.0))) then
         !
         write(logstr,'(a)')' WARNING! Times in boundary conditions file do not cover entire simulation period!'
         call write_log(logstr, 1)
         !
      endif   
      !
   elseif (bndfile(1:4) /= 'none') then    ! Normal ascii input files
      !
      call write_log('Info    : reading water level boundaries', 0)
      !
      open(500, file=trim(bndfile))
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         nbnd = nbnd + 1
      enddo
      rewind(500)
      allocate(x_bnd(nbnd))
      allocate(y_bnd(nbnd))
      do n = 1, nbnd
         read(500,*)x_bnd(n),y_bnd(n)
      enddo
      close(500)
      !
      ! Read water level boundary conditions file
      !
      if (bzsfile(1:4) /= 'none') then
         !
         open(500, file=trim(bzsfile))
         do while(.true.)
            read(500,*,iostat = stat)dummy
            if (stat<0) exit
            ntbnd = ntbnd + 1
         enddo
         rewind(500)
         !
         allocate(t_bnd(ntbnd))
         allocate(zs_bnd(nbnd,ntbnd))
         allocate(zst_bnd(nbnd))
         !
         do itb = 1, ntbnd
            read(500,*)t_bnd(itb),(zs_bnd(ib, itb), ib = 1, nbnd)
         enddo
         !
         close(500)
         !
      else
         !
         ntbnd = 2
         !
         allocate(t_bnd(ntbnd))
         allocate(zs_bnd(nbnd,ntbnd))
         allocate(zst_bnd(nbnd))
         !
         t_bnd(1) = t0
         t_bnd(2) = t1
         zs_bnd = 0.0
         !
      endif   
      !
      if (bzifile(1:4) /= 'none') then
         !
         ! Incoming infragravity waves
         !
         allocate(zsi_bnd(nbnd,ntbnd))
         allocate(zsit_bnd(nbnd))
         !
         open(500, file=trim(bzifile))
         do itb = 1, ntbnd
            read(500,*)dummy,(zsi_bnd(ib, itb), ib = 1, nbnd)
         enddo
         close(500)
         !
      endif
      !      
      if ((t_bnd(1) > (t0 + 1.0)) .or. (t_bnd(ntbnd) < (t1 - 1.0))) then
         ! 
         write(logstr,'(a)')'Warning! Times in boundary conditions file do not cover entire simulation period !'
         call write_log(logstr, 1)
         !
         if (t_bnd(1) > (t0 + 1.0)) then
            ! 
            write(logstr,'(a)')'Warning! Adjusting first time in boundary conditions time series !'
            call write_log(logstr, 1)
            !
            t_bnd(1) = t0 - 1.0
            !
         else
            ! 
            write(logstr,'(a)')'Warning! Adjusting last time in boundary conditions time series !'
            call write_log(logstr, 1)
            !
            t_bnd(ntbnd) = t1 + 1.0
            !
         endif
         !
      endif   
      !
   elseif (include_boundaries) then   
      !
      write(logstr,'(a)')'Warning! Boundary cells found in mask, without boundary conditions. Using water level of 0.0 m at these points.'
      call write_log(logstr, 1)
      !
      nbnd = 1
      allocate(x_bnd(nbnd))
      allocate(y_bnd(nbnd))
      x_bnd(1) = 0.0
      y_bnd(1) = 0.0
      !
      ntbnd = 2
      !
      allocate(t_bnd(ntbnd))
      allocate(zs_bnd(nbnd,ntbnd))
      allocate(zst_bnd(nbnd))
      !
      t_bnd(1) = t0
      t_bnd(2) = t1
      zs_bnd(1,1) = 0.0
      zs_bnd(2,1) = 0.0
      zs_bnd(1,2) = 0.0
      zs_bnd(2,2) = 0.0
      zst_bnd(1)  = 0.0      
      !
   endif
   !
   ! Check for 'weird' values
   !
   iok = 1
   ! 
   do ib = 1, nbnd
      do itb = 1, ntbnd
         !
         if (zs_bnd(ib, itb)<-99.0 .or. zs_bnd(ib, itb)>990.0) then
            !
            iok = 0
            !
            zs_bnd(ib, itb) = 0.0
            !
         endif
         !
      enddo 
   enddo    
   !
   if (iok == 0) then
      ! 
      write(logstr,'(a)')'Warning! Very low, very high or NaN values found in boundary conditions file ! These have now been replaced with zeros. Please check !'
      call write_log(logstr, 1)
      ! 
   endif
   !
   ! Read wave boundaries
   !
   nwbnd = 0
   ntwbnd = 0
   !
   if (bwvfile(1:4) /= 'none') then
      !
      write(logstr,'(a)')'Info    : reading wave boundaries'
      call write_log(logstr, 0)
      !
      ! Locations
      !
      open(500, file=trim(bwvfile))
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
      allocate(hst_bwv(nwbnd))
      allocate(l0t_bwv(nwbnd))
      allocate(tpt_bwv(nwbnd))
      allocate(wdt_bwv(nwbnd))
!      allocate(setupt_bwv(nwbnd))
      !
      ! Hs (significant wave height)
      ! Times in btp and bwd files must be the same as in bhs file!
      !
      open(500, file=trim(bhsfile))
      allocate(hs_bwv(nwbnd,ntwbnd))
      do itb = 1, ntwbnd
         read(500,*)t_bwv(itb),(hs_bwv(ib, itb), ib = 1, nwbnd)
      enddo
      close(500)
      !
      ! Tp (peak period)
      !
      open(500, file=trim(btpfile))
      allocate(tp_bwv(nwbnd,ntwbnd))
      do itb = 1, ntwbnd
         read(500,*)t_bwv(itb),(tp_bwv(ib, itb), ib = 1, nwbnd)
      enddo
      close(500)
      !
      if (bwdfile(1:4) /= 'none') then
         !
         ! Wd (wave direction)
         !
         open(500, file=trim(bwdfile))
         allocate(wd_bwv(nwbnd,ntwbnd))
         do itb = 1, ntwbnd
            read(500,*)t_bwv(itb),(wd_bwv(ib, itb), ib = 1, nwbnd)
         enddo
         close(500)
         wd_bwv = (270.0 - wd_bwv)*pi/180 ! Convert to cartesian going to
         !
      endif
      !
!      if (stufile(1:4) /= 'none') then
!         !
!         ! Set-up
!         !
!         open(500, file=trim(stufile))
!         allocate(setup_bwv(nwbnd, ntwbnd))
!         do itb = 1, ntwbnd
!            read(500,*)t_bwv(itb),(setup_bwv(ib, itb), ib = 1, nwbnd)
!         enddo
!         close(500)
!         !
!      endif
      !
   endif
!   !
!   ! Infragravity frequencies
!   !
!   allocate(freqig(nfreqsig))
!   allocate(costig(nfreqsig))
!   allocate(phiig(nfreqsig))
!   allocate(dphiig(nfreqsig))
!   dfreqig = (freqmaxig - freqminig)/nfreqsig
!   do ifreq = 1, nfreqsig
!      freqig(ifreq) = freqminig + ifreq*dfreqig - 0.5*dfreqig
!      call RANDOM_NUMBER(r)
!      phiig(ifreq) = r*2*3.1416
!!      call RANDOM_NUMBER(r)
!      dphiig(ifreq) = 1.0e-6*2*3.1416/freqig(ifreq)
!   enddo

   !
   ! If bca file is present, read it
   !
   nr_tidal_components = 0
   !
   if (bcafile(1:4) /= 'none' .and. nbnd >0) then
      !
      ! First pass: count number of sets and number of components in first set
      !
      ! NOTE 1 : The number of component sets in the bca file must match the number of points in the bnd file!
      ! NOTE 2 : The number of components in each set must be the same!
      !
      n_sets = 0
      n_components = 0
      has_a0 = .false.
      !
      open(500, file=trim(bcafile))
      !
      do
         !
         read(500, '(A)', iostat=ios) line
         !
         if (ios /= 0) exit
         !
         if (trim(line) == '[forcing]') then
            !
            n_sets = n_sets + 1
            !
         elseif (index(line, 'Quantity') > 0 .or. index(line, 'Unit') > 0 .or. index(line, 'Name') > 0 .or. index(line, 'Function') > 0) then
            !
            cycle  ! skip line
            !
         else
            !
            ! Try to read a data line
            !
            read(line, *, iostat=ios) cname, a, p
            !
            if (ios /= 0) cycle  ! not a data line
            !
            if (n_sets == 1) then  ! only count on first set
               !
               n_components = n_components + 1
               !
               ! Check if first component is A0
               !
               if (n_components == 1) then
                  !
                  if (cname(1:2) == 'a0' .or. cname(1:2) == 'A0') then
                     !
                     has_a0 = .true.
                     !
                  endif   
                  !
               endif   
               !
            endif   
            !
         endif
      enddo
      !
      ! Check whether n_sets is same as nbnd (otherwise give error)
      !
      if (n_sets /= nbnd) then
         !
         write(*,*)'ERROR! Number of astronomical tidal datasets in *.bca file does not match number of boundary points in *.bnd file !'
         !
      endif
      !
      rewind(500)
      !
      if (.not. has_a0) then
         !
         ! Add 1 for A0
         !
         n_components = n_components + 1
         !
      endif
      !
      ! Allocate memory for names, amplitude, phase and frequency
      !      
      nr_tidal_components = n_components
      !
      allocate(tidal_component_names(n_components))
      allocate(tidal_component_frequency(n_components))
      allocate(tidal_component_data(2, n_components, nbnd))
      !
      ! Initialize
      !
      tidal_component_data = 0.0
      tidal_component_names = ''
      !
      ! Second pass: read and parse data
      !
      current_set = 0
      current_component = 0
      !
      do
         read(500, '(A)', iostat=ios) line
         !
         if (ios /= 0) exit
         !
         if (trim(line) == '[forcing]') then
            !
            current_set = current_set + 1
            current_component = 0
            !
            if (.not. has_a0) then
               !
               current_component = 1
               tidal_component_names(1) = 'A0      '
               !
               ! Data has already be initialized at 0.0
               !
            endif   
            !
         elseif (index(line, 'Quantity') > 0 .or. index(line, 'Unit') > 0 .or. index(line, 'Name') > 0 .or. index(line, 'Function') > 0) then
            !
            cycle  ! skip metadata
            !
         else
            !
            ! Try to read a data line
            !
            read(line, *, iostat=ios) cname, a, p
            !
            if (ios /= 0) cycle  ! not a data line
            !
            current_component = current_component + 1
            tidal_component_names(current_component) = cname
            tidal_component_data(1, current_component, current_set) = a
            tidal_component_data(2, current_component, current_set) = p
            !
         endif
         !
      enddo
      !
      close(500)
      !      
      i_date_time = time_to_vector(0.0d0, trefstr)
      !
      !
      call update_nodal_factors(i_date_time, tidal_component_names, nr_tidal_components, nbnd, tidal_component_data, tidal_component_frequency)
      !
      tidal_component_frequency = tidal_component_frequency / 3600 ! Convert to rad/s      
      !
      !do ios = 1, nr_tidal_components
      !   write(*,'(a,20f16.3)')tidal_component_names(ios), (180.0 / pi) * tidal_component_frequency(ios) * 3600.0,tidal_component_data(1,ios,1), (180.0 / pi) * tidal_component_data(2,ios,1)
      !enddo   
      !
   endif      
   !
   end subroutine


   subroutine find_boundary_indices()
   !
   use sfincs_data
   !
   implicit none
   !
   ! For each grid boundary point (kcs=2) :
   !
   ! Determine indices and weights of boundary points
   ! For tide and surge, these are the indices and weights of the points in the bnd file
   ! For waves, these are the indices and weights of the points in the cst file
   !
   integer nm, m, n, nb, ib1, ib2, ib, ic
   !
   real x, y, dst1, dst2, dst
   !
   if (nbnd>0 .and. ngbnd>0) then
      !
      ! Count number of boundary points
      !
      ! Allocate boundary arrays
      !
      ! First water levels
      !
      ! Water level arrays
      !
      allocate(ind1_bnd_gbp(ngbnd))
      allocate(ind2_bnd_gbp(ngbnd))
      allocate(fac_bnd_gbp(ngbnd))
      !
      ! Water levels at boundary
      !
      if (waves) then
         !
         ! Wave arrays
         !
         allocate(ind1_cst_gbp(ngbnd))
         allocate(ind2_cst_gbp(ngbnd))
         allocate(fac_cst_gbp(ngbnd))
         !
      endif
      !
      ! Find two closest boundary condition points for each boundary point
      ! And the two closest coastline points
      !
      nb = 0
      !
      ! Loop through all grid points
      !
      do nm = 1, np
         !
         ! Check if this point is a boundary point
         !
         if (kcs(nm) > 1) then
            !
            nb = nb + 1
            !
            x = z_xz(nm)
            y = z_yz(nm)
            !
            ! Indices and weights for water level boundaries
            !
            if (nbnd>1) then
               !
               dst1 = 1.0e10
               dst2 = 1.0e10
               ib1 = 0
               ib2 = 0
               !
               ! Loop through all water level boundary points
               !
               do ib = 1, nbnd
                  !
                  ! Compute distance of this point to grid boundary point
                  !
                  dst = sqrt((x_bnd(ib) - x)**2 + (y_bnd(ib) - y)**2)
                  !
                  if (dst<dst1) then
                     !
                     ! Nearest point found
                     !
                     dst2 = dst1
                     ib2  = ib1
                     dst1 = dst
                     ib1  = ib
                     !
                  elseif (dst<dst2) then
                     !
                     ! Second nearest point found
                     !
                     dst2 = dst
                     ib2  = ib
                     !
                  endif
               enddo
               !
               ind1_bnd_gbp(nb)  = ib1
               ind2_bnd_gbp(nb)  = ib2
               fac_bnd_gbp(nb) = dst2/max(dst1 + dst2, 1.0e-9)
               !
            else
               !
               ind1_bnd_gbp(nb)  = 1
               ind2_bnd_gbp(nb)  = 1
               fac_bnd_gbp(nb)   = 1.0
               !
            endif
            !
            ! Indices and weights for wave boundaries
            !
            if (waves) then
               if (ncst>1) then
                  !
                  dst1 = 1.0e10
                  dst2 = 1.0e10
                  ib1 = 0
                  ib2 = 0
                  !
                  ! Loop through all water level boundary points
                  !
                  do ic = 1, ncst
                     !
                     ! Compute distance of this point to grid boundary point
                     !
                     dst = sqrt((x_cst(ic) - x)**2 + (y_cst(ic) - y)**2)
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
                  ind1_cst_gbp(nb)  = ib1
                  ind2_cst_gbp(nb)  = ib2
                  fac_cst_gbp(nb) = dst2/(dst1 + dst2)
                  !
               else
                  !
                  ind1_cst_gbp(nb)  = 1
                  ind2_cst_gbp(nb)  = 1
                  fac_cst_gbp(nb)   = 1.0
                  !
               endif
            endif
            !
         endif
      enddo
   endif
   !
   end subroutine


   subroutine update_boundary_points(t)
   !
   ! Update values at boundary points
   !
   use sfincs_data
   !
   implicit none
   !
   integer ib, itb, itb0, itb1, ic
   !
   real*8 t
   !
   real*4 zstb, tbfac, hs, tp, wd, tb
   !
   if (nbnd > 0) then
      !
      ! Start with updating values at boundary polylines
      !
      ! Water levels
      !
      ! Interpolate boundary conditions in time
      !
      if (t_bnd(1) > (t - 1.0e-3)) then ! use first time in boundary conditions
         !
         itb0 = 1
         itb1 = 1
         tb   = t_bnd(itb0)
         !
      elseif (t_bnd(ntbnd) < (t + 1.0e-3)) then  ! use last time in boundary conditions       
         !
         itb0 = ntbnd
         itb1 = ntbnd
         tb   = t_bnd(itb0)
         !
      else
         !
         do itb = itbndlast, ntbnd ! Loop in time
            if (t_bnd(itb) > (t + 1.0e-6)) then
               itb0 = itb - 1
               itb1 = itb
               tb   = t
               itbndlast = itb - 1
               exit
            endif
         enddo 
         !
      endif            
      !
      tbfac  = (tb - t_bnd(itb0))/max(t_bnd(itb1) - t_bnd(itb0), 1.0e-6)
      !
      do ib = 1, nbnd ! Loop along boundary points
         !
         ! Tide and surge
         !
         zstb = zs_bnd(ib, itb0) + (zs_bnd(ib, itb1) - zs_bnd(ib, itb0))*tbfac
         !
         ! Add astronomical tides
         !
         if (nr_tidal_components > 0) then
            !
            do ic = 1, nr_tidal_components
               !
               zstb = zstb + tidal_component_data(1, ic, ib) * cos(tidal_component_frequency(ic) * t - tidal_component_data(2, ic, ib))
               !
            enddo   
            !
         endif   
         !
         zst_bnd(ib) = zstb
         !
         if (bzifile(1:4) /= 'none') then
            !
            ! Incoming infragravity waves
            !
            zsit_bnd(ib) = zsi_bnd(ib, itb0) + (zsi_bnd(ib, itb1) - zsi_bnd(ib, itb0))*tbfac
            !
         endif
         !
      enddo
      !      
   endif
   !
   end subroutine

   
   
   subroutine update_boundary_conditions(t, dt)
   !
   ! Update values at boundary points
   !
   use sfincs_data
   !
   implicit none
   !
   integer ib, ifreq, ic
   !
   real*8                            :: t
   real                              :: dt
   !
   real*4 zst, hst, l0t, wdt, ksi, zsetup, zig, a, angdiff, angfac, hm0ig, tmmin01ig, fp
   real*4 zs0act
   real*4 smfac
   real*4 zs0smooth
   !
   ! Set water level in all boundary points on grid
   !
   do ib = 1, ngbnd
      !
      if (ibndtype(ib) == 1) then
         !
         ! Regular boundary point (otherwise (for kcs==3) zsb and zsb0 have already been initialized at zbmin)
         !
         ! Water levels (surge+tide)
         !
         if (nbnd>1) then
            !
            ! Interpolation of nearby points
            !
            zst   = zst_bnd(ind1_bnd_gbp(ib))*fac_bnd_gbp(ib)  + zst_bnd(ind2_bnd_gbp(ib))*(1.0 - fac_bnd_gbp(ib))
            !
         else
            !
            ! Just use the value of the one boundary point
            !
           zst   = zst_bnd(1)
            !
         endif
         !
         if (patmos .and. pavbnd>1.0) then
            !
            ! Barometric pressure correction
            !
            zst = zst + ( pavbnd - patmb(ib)) / (rhow*9.81)
            !
         endif
         !
         zsetup = 0.0
         zig    = 0.0
         !
         ! Incoming IG waves from file (this will overrule IG signal computed before)
         !
         if (bzifile(1:4) /= 'none') then
            !
            if (nbnd>1) then
               !
               ! Interpolation of nearby points
               !
               zig   = zsit_bnd(ind1_bnd_gbp(ib))*fac_bnd_gbp(ib)  + zsit_bnd(ind2_bnd_gbp(ib))*(1.0 - fac_bnd_gbp(ib))
               !
            else
               !
               ! Just use the value of the one boundary point
               !
               zig   = zsit_bnd(1)
               !
            endif
            !
         endif
         !
         if (t < (tspinup - 1.0e-3)) then
            !
            smfac = 1.0 - (t - t0)/(tspinup - t0)
            !
            zs0act = zst + zsetup
            call weighted_average(zini, zs0act, smfac, 1, zs0smooth)
            zsb0(ib) = zs0smooth           ! Still water level at kcs=2 point
            !
            zs0act = zst + zsetup + zig
            call weighted_average(zini, zs0act, smfac, 1, zs0smooth)
            zsb(ib) = zs0smooth            ! Total water level at kcs=2 point
            !
         else
            !
            zsb0(ib) = zst + zsetup        ! Still water level at kcs=2 point
            zsb(ib)  = zst + zsetup + zig  ! Total water level at kcs=2 point
            !
         endif
         !
         if (subgrid) then                  ! Check on waterlevels minimally equal to z_zmin
            zsb0(ib) = max(zsb0(ib), subgrid_z_zmin(nmindbnd(ib)))
            zsb(ib)  = max(zsb(ib),  subgrid_z_zmin(nmindbnd(ib)))
         else                               ! Check on waterlevels minimally equal to zb           
            zsb0(ib) = max(zsb0(ib), zb(nmindbnd(ib)))
            zsb(ib)  = max(zsb(ib),  zb(nmindbnd(ib)))
         endif
         !         
      endif
   enddo
   !
   end subroutine


   subroutine update_boundary_fluxes(dt)
   !
   ! Update fluxes qx and qy at boundary points
   !
   use sfincs_data
   !
   implicit none
   !
   integer ib, nm, nmi, nmb, iuv, indb, ip
   real*4  hnmb, dt, zsnmi, zsnmb, zs0nmb, facrel
   real*4  factime, one_minus_factime
   !
   real*4 ui, ub, dzuv, facint, zsuv, depthuv
   !
   !$acc update device( zsb0, zsb ), async(1)
   !
   factime = min(dt/btfilter, 1.0)
   one_minus_factime = 1.0 - factime
   facrel  = 1.0 - min(dt/btrelax, 1.0)
   !
   ! UV fluxes at boundaries
   !
   !$acc kernels present(index_kcuv2, nmikcuv2, nmbkcuv2, ibkcuv2, kcuv, zs, z_volume, q, uvmean, uv, zb, zbuv, zsb, zsb0, &
   !$acc                 subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_havg, subgrid_uv_havg_zmax, subgrid_z_zmin, ibuvdir, zsmax ), async(1)
   !$acc loop independent, private(ib)
   do ib = 1, nkcuv2
      !
      ip     = index_kcuv2(ib)  ! Index in uv array of kcuv=2 velocity point
      !
      nmi    = nmikcuv2(ib)     ! Index of kcs=1 point
      nmb    = nmbkcuv2(ib)     ! Index of kcs=2/3 boundary point
      indb   = ibkcuv2(ib)
      !
      zsnmi  = zs(nmi)          ! total water level inside model
      zsnmb  = zsb(indb)        ! total water level at boundary
      zs0nmb = zsb0(indb)       ! average water level inside model (without waves)
      !
      if (bndtype == 1) then
         !
         ! Weakly reflective boundary (default)
         !          
         if (subgrid) then
            !
            zsuv = max(zsnmb, zsnmi)
            !
            if (zsuv>=subgrid_uv_zmax(ip) - 1.0e-4) then
               !
               ! Entire cell is wet, no interpolation from table needed
               !
               depthuv  = subgrid_uv_havg_zmax(ip) + zsuv
               !
            elseif (zsuv>subgrid_uv_zmin(ip)) then
               !
               ! Interpolation required
               !
               dzuv    = (subgrid_uv_zmax(ip) - subgrid_uv_zmin(ip)) / (subgrid_nlevels - 1)
               iuv     = int((zsuv - subgrid_uv_zmin(ip))/dzuv) + 1
               facint  = (zsuv - (subgrid_uv_zmin(ip) + (iuv - 1)*dzuv) ) / dzuv
               depthuv = subgrid_uv_havg(iuv, ip) + (subgrid_uv_havg(iuv + 1, ip) - subgrid_uv_havg(iuv, ip))*facint
               !
            else
               !
               depthuv = 0.0
               !
            endif
            !
            hnmb   = max(depthuv, huthresh)
            zsnmb  = max(zsnmb,  subgrid_z_zmin(nmb))
            zs0nmb = max(zs0nmb, subgrid_z_zmin(nmb))
            !
         else
            !
            hnmb   = max(0.5 * (zsnmb + zsnmi) - zbuv(ip), huthresh)
            zsnmb  = max(zsnmb,  zb(nmb))
            zs0nmb = max(zs0nmb, zb(nmb))
            !
         endif      
         !
         if (hnmb<huthresh + 1.0e-6 .or. kcuv(ip)==3) then
            !
            ! Very shallow or also a structure point.
            !
            q(ip)      = 0.0
            uv(ip)     = 0.0
            uvmean(ib) = 0.0
            !
         else
            !
            ui = sqrt(g / hnmb) * (zsnmb - zs0nmb)
            ub = ibuvdir(ib) * (2 * ui - sqrt(g / hnmb) * (zsnmi - zs0nmb))
            !
            q(ip) = ub * hnmb + uvmean(ib)            
            !
            ! Riemann
            !
            ! R = ubnd + 2 * sqrt(g * hnmb) ! from bca or bzs/buv
            !
            ! hnmi = hnmb + zsnmi - zsnmb
            !
            ! ub = R - 2 * sqrt(g * hnmi)
            !
            ! q(ip) = ub * hnmb
            !
            if (subgrid) then
               !
               ! Sub-grid
               !
               if (z_volume(nmi)<=0.0) then
                  if (ibuvdir(ib)==1) then
                     q(ip) = max(q(ip), 0.0) ! Nothing can flow out
                  else
                     q(ip) = min(q(ip), 0.0) ! Nothing can flow out
                  endif
               endif
               !
               if (zsnmb - subgrid_z_zmin(nmb)<huthresh) then
                  if (ibuvdir(ib)==1) then
                     q(ip) = min(q(ip), 0.0) ! Nothing can flow in
                  else
                     q(ip) = max(q(ip), 0.0) ! Nothing can flow in
                  endif
               endif
               !
            else
               !
               ! Regular
               !
               if (zsnmi - zb(nmi)<=huthresh) then                  
                  if (ibuvdir(ib)==1) then
                     q(ip) = max(q(ip), 0.0) ! Nothing can flow out
                  else
                     q(ip) = min(q(ip), 0.0) ! Nothing can flow out
                  endif
               endif
               !
               if (zsnmb - zb(nmb)<huthresh) then
                  if (ibuvdir(ib)==1) then
                     q(ip) = min(q(ip), 0.0) ! Nothing can flow in
                  else
                     q(ip) = max(q(ip), 0.0) ! Nothing can flow in
                  endif
               endif
            endif
            !
            ! Limit velocities (this does not change fluxes, but may prevent advection term from exploding in the next time step)
            !
            uv(ip)  = max(min(q(ip)/hnmb, 4.0), -4.0)
            !
         endif
         !
         if (btfilter>=-1.0e-6) then
            !
            ! Added a little bit of relaxation in uvmean to avoid persistent jets shooting into the model
            ! Using: facrel = 1.0 - min(dt/btrelax, 1.0)
            ! Default: btrelax = 3600 s.
            !
            uvmean(ib) = factime * q(ip) + facrel * one_minus_factime * uvmean(ib)
            !
         else
            !
            uvmean(ib) = 0.0
            !
         endif
         !
         ! Set value on kcs=2 point to boundary condition
         !
         zs(nmb) = zsb(indb)
         !
         ! Store maximum water levels also on the boundary
         !
         if (store_maximum_waterlevel) then
            zsmax(nmb) = max(zsmax(nmb), zs(nmb))
         endif
         !
      endif
      !
   enddo
   !
   !$acc end kernels
   !
   end subroutine

   
   
   subroutine update_boundaries(t, dt, tloop)
   !
   ! Update all boundary conditions
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
   !
   call system_clock(count0, count_rate, count_max)
   !
   if (include_boundaries) then
      !
      if (nbnd>0) then
         !
         ! Update boundary conditions at boundary points
         !
         call update_boundary_points(t)
         !
         ! Update boundary conditions at grid points (water levels)
         !
         call update_boundary_conditions(t, dt)
         !
      endif
      !
      ! Update boundary fluxes()
      !
      call update_boundary_fluxes(dt)
      !
   endif
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine
   !
   !
   !
   subroutine weighted_average(val1,val2,fac,iopt,val3)
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

end module
