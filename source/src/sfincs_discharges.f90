module sfincs_discharges
   
   use sfincs_log
   use sfincs_error
   
   type :: NormalFlow
     real*4 :: n    ! velocity exponent
     real*4 :: stt  ! sediment total transport capacity
     real*4 :: normal_water_depth
     real*4 :: normal_flow_velocity
     real*4 :: friction_coeff
     real*4 :: water_depth
     real*4 :: slope_flow_velocity
     real*4 :: discharge
     real*4 :: t1
     real*4 :: t2
     real*4 :: breach_width_waterline
     real*4 :: breach_width_total
     real*4 :: gamma0
     real*4 :: breach_bottom
     real*4 :: breach_level
   end type NormalFlow


contains
   !
   subroutine read_discharges()
   !
   ! Reads discharge files
   !
   use sfincs_data
   use sfincs_ncinput
   use quadtree
   !
   implicit none
   !
   real*4, dimension(:),     allocatable :: xsnk
   real*4, dimension(:),     allocatable :: ysnk
   !
   real*4 dummy, xsnk_tmp, ysnk_tmp, xsrc_tmp, ysrc_tmp
   !
   integer isrc, itsrc, idrn, nm, m, n, stat, j, iref, nmq, npars
   !
   logical :: ok
   !
   character(len=256) :: drainage_line, message 
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
      ok = check_file_exists(srcfile, 'Source points file', .true.)
      !
      write(logstr,'(a)')'Info    : reading discharges'
      call write_log(logstr, 0)
      !
      ok = check_file_exists(srcfile, 'River input locations src file', .true.)      
      !
      open(500, file=trim(srcfile))
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat < 0) exit
         nsrc = nsrc + 1
      enddo
      rewind(500)
      !
   elseif (netsrcdisfile(1:4) /= 'none') then    ! FEWS compatible Netcdf discharge time-series input
      !
      ok = check_file_exists(netsrcdisfile, 'Netcdf river input netsrcdis file', .true.)       
      !
      call read_netcdf_discharge_data()  ! reads nsrc, ntsrc, xsrc, ysrc, qsrc, and tsrc
      !
      if ((tsrc(1) > (t0 + 1.0)) .or. (tsrc(ntsrc) < (t1 - 1.0))) then
         !
         write(logstr,'(a)')' WARNING! Times in discharge file do not cover entire simulation period!'
         call write_log(logstr, 1)
         !
      endif         
      !
   endif   
   !
   if (drnfile(1:4) /= 'none') then
      !
      write(logstr,'(a)')'Info    : reading drainage file'
      call write_log(logstr, 0)
      !
      ok = check_file_exists(drnfile, 'Drainage points drn file', .true.)
      !
      open(501, file=trim(drnfile))
      do while(.true.)
         read(501,*,iostat = stat)dummy
         if (stat < 0) exit
         ndrn = ndrn + 1
      enddo
      rewind(501)
   endif
   !
   nsrcdrn = nsrc + 2 * ndrn
   !
   if (nsrcdrn > 0) then
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
         read(500,*)xsrc(n), ysrc(n)
      enddo
      close(500)
      !
      ! Read discharge time series
      !
      ok = check_file_exists(disfile, 'River discharge timeseries dis file', .true.)      
      !      
      open(502, file=trim(disfile))
      do while(.true.)
         read(502,*,iostat = stat)dummy
         if (stat < 0) exit
         ntsrc = ntsrc + 1
      enddo
      rewind(502)
      allocate(tsrc(ntsrc))
      allocate(qsrc(nsrc,ntsrc))
      do itsrc = 1, ntsrc
         read(502,*)tsrc(itsrc), (qsrc(isrc, itsrc), isrc = 1, nsrc)
      enddo
      close(502)
      !
      if ((tsrc(1) > (t0 + 1.0)) .or. (tsrc(ntsrc) < (t1 - 1.0))) then
         ! 
         write(logstr,'(a)')'Warning! Times in discharge file do not cover entire simulation period !'
         call write_log(logstr, 1)
         !
         if (tsrc(1) > (t0 + 1.0)) then
            ! 
            write(logstr,'(a)')'Warning! Adjusting first time in discharge time series !'
            call write_log(logstr, 1)
            !
            tsrc(1) = t0 - 1.0
            !
         else
            ! 
            write(logstr,'(a)')'Warning! Adjusting last time in discharge time series !'
            call write_log(logstr, 1)
            !
            tsrc(ntsrc) = t1 + 1.0
            !
         endif
         !
      endif   
      !
   endif  
   !
   if (nsrc > 0) then
      !
      ! Determine m and n indices of sources
      !
      do isrc = 1, nsrc
         !
         ! Find cell in quadtree first
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
      write(logstr,'(a,a,a,i0,a)')' Reading ',trim(drnfile),' (', ndrn, ' drainage points found) ...'
      call write_log(logstr, 0)
      !         
      allocate(xsrc(ndrn))
      allocate(ysrc(ndrn))
      allocate(xsnk(ndrn))
      allocate(ysnk(ndrn))
      
      
      !
      allocate(drainage_type(ndrn))
      allocate(drainage_params(ndrn, 6))
      allocate(drainage_status(ndrn))
      allocate(drainage_distance(ndrn))
      allocate(drainage_fraction_open(ndrn))
      !
      drainage_params = 0.0
      drainage_distance = 0.0
      drainage_fraction_open = 1.0   ! initially fully open (should fix this based on zmin and zmax in params)
      drainage_status = 1            ! open (0=closed, 1=open, 2=closing, 3=opening)
      
      !
      do idrn = 1, ndrn
         !
         read(501, '(a)') drainage_line
         !
         ! First find out what type of drainage structure it is (integer 5th item in line)
         !
         read(drainage_line,*,iostat=stat)xsnk_tmp, ysnk_tmp, xsrc_tmp, ysrc_tmp, drainage_type(idrn)
         !
         npars = 0 ! Default (if npars stays 0, throw error)
         !
         if (drainage_type(idrn) == 1 .or. drainage_type(idrn) == 2 .or. drainage_type(idrn) == 3) then
            !
            ! Pump, culvert or check valve (1 parameter)
            !
            npars = 1
            !
         elseif (drainage_type(idrn) == 4 .or. drainage_type(idrn) == 5) then
            !
            ! Controlled gate (6 parameters : width, sill elevation, manning, zmin, zmax, closing time)
            !
            npars = 6
            !
         elseif (drainage_type(idrn)==6 .or. drainage_type(idrn)==7 .or. drainage_type(idrn)==8 .or. drainage_type(idrn)==9) then
            !
			! Dike breaching
            ! Drainage_type 6: Verheij with submerge or free flow formula (6 parameters: z_crest, tbreach, z_min, B0, t_0, dike_core)
            ! Drainage_type 7: Verheij with Bates or free flow formula (6 parameters: z_crest, tbreach, z_min, B0, t_0, dike_core)
            ! Drainage_type 8: Tadesse with submerge or free flow formula (6 parameters: time of breaching, breach duration, final_breach_width, z_min, crest level of the dike, initial breach width)
            ! Drainage_type 9: Visser with submerge or free flow formula (working progress)
            
            ! These are needed to remember the previous breach width and write the breach width and breach level as output
            allocate(breach_width(ndrn))
            allocate(breach_level_gather(ndrn))
             
            ! These are needed for Visser
            allocate(running_Visser_phase1(ndrn))
            allocate(running_Visser_phase2(ndrn))
            allocate(discharge_t1(ndrn))
            allocate(t1_Visser(ndrn))
            allocate(breach_bottom_Visser(ndrn))
            allocate(breach_width_waterline_Visser(ndrn))
            allocate(gamma0_Visser(ndrn))
            allocate(discharge_t2(ndrn))
            allocate(t2_Visser(ndrn))
            allocate(end_breaching(ndrn)) ! After phase 5 is done you can't go back to phase 4
            
            npars = 6
            !
         endif
         !
         if (npars == 0) then
            !
            write(logstr,'(a,i0,a)')'Drainage type ', drainage_type(idrn), ' not recognized !'
            call stop_sfincs(logstr, -1)
            !
         endif               
         !
         if (npars == 1) then
            !
            ! Pump, culvert or check valve
            !
            read(drainage_line,*,iostat=stat)xsnk(idrn), ysnk(idrn), xsrc(idrn), ysrc(idrn), drainage_type(idrn), drainage_params(idrn,1)
            !
         elseif (npars == 6) then
            !
            ! Controlled gate or dike breaching, needs 6 parameters
            !
            read(drainage_line,*,iostat=stat)xsnk(idrn), ysnk(idrn), xsrc(idrn), ysrc(idrn), drainage_type(idrn), drainage_params(idrn,1), drainage_params(idrn,2), drainage_params(idrn,3), drainage_params(idrn,4), drainage_params(idrn,5), drainage_params(idrn,6)
            !
         endif
         !
         if (stat /= 0) then
            !
            write(logstr,'(a,i0,a,i0,a)')'Drainage type ', drainage_type(idrn), ' requires ', npars, ' parameters !'
            call stop_sfincs(logstr, -1)
            !
         endif
         !
      enddo
      !
      close(501)
      !
      ! Determine nm indices of source and sinks
      !
      do idrn = 1, ndrn
         !
         ! Determine index of sink first
         !
         j = nsrc + idrn*2 - 1
         !
         nmq = find_quadtree_cell(xsnk(idrn), ysnk(idrn))
         !
         if (nmq > 0) then
            !
            nmindsrc(j) = index_sfincs_in_quadtree(nmq)
            !
         endif
         !
         ! And now the index of the source
         !
         j = nsrc + idrn * 2
         !
         nmq = find_quadtree_cell(xsrc(idrn), ysrc(idrn))
         !
         if (nmq > 0) then
            !
            nmindsrc(j) = index_sfincs_in_quadtree(nmq)
            !
         endif
         !
         ! Get coords of source and sink points, and compute distance between them
         ! This is needed for controlled gates (type 4)
         !
         xsnk_tmp = z_xz(nmindsrc(nsrc + idrn * 2 - 1))
         ysnk_tmp = z_yz(nmindsrc(nsrc + idrn * 2 - 1))
         xsrc_tmp = z_xz(nmindsrc(nsrc + idrn * 2))
         ysrc_tmp = z_yz(nmindsrc(nsrc + idrn * 2))
         !
         drainage_distance(idrn) = sqrt( (xsrc_tmp - xsnk_tmp)**2 + (ysrc_tmp - ysnk_tmp)**2 )
         !
      enddo
      !
      deallocate(xsrc)
      deallocate(ysrc)
      deallocate(xsnk)
      deallocate(ysnk)
      !
      ! Check if all sink/source points have found an index
      if (any(nmindsrc == 0)) then
         !
         write(logstr,'(a)')'Warning ! For some sink/source drainage points no matching active grid cell was found!'
         call write_log(logstr, 0)
         write(logstr,'(a)')'Warning ! These points will be skipped, please check your input!'
         call write_log(logstr, 0)
         !
      endif      
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
   use sfincs_Visser
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
   real*4           :: dzds, frac, wdt, zsill, zmin, zmax, mng, hgate, dfrac, tcls, topen, tclose
   real*4           :: m_afvoercoeff, dike_core, tbreach, z_crest, z_min, t_0, t_phase1, B0, Z, B, f1, f2, uc, dB, B_old, denom
   real*4           :: breach_duration, final_breach_width, t_end, dt_hr, tau_hr, H, dBdt
   real*4           :: stage, beta, beta1, beta0, breach_bottom, breach_level, test_value
   real*4           :: alpha, crest_level, W, polder_level, t0_Visser, crest_width
   real*4           :: d50, d90, Cf, Kappa, water_temperature, p, phi
   real*4           :: sediment_density, water_density, delta, dyn_viscosity, sediment_fall_velocity
   real*4           :: dstar, theta_crit, ni, k, breach_width_waterline, breach_width_total
   real*4           :: outside_water_level, polder_water_level, gamma0, gamma1, outside_level, h_breach, dike_width
   real*4           :: discharge_coeff, crit_water_depth, flow_velocity, water_depth, breach_width_avg_water_depth
   character*256 :: formula
   type(NormalFlow) :: results_t1,results_t2,results_t3,results_t4,results_t5
   
   
   !
   integer          :: idir
   !
   integer isrc, itsrc, idrn, jin, jout, nmin, nmout
   !
   call system_clock(count0, count_rate, count_max)
   !
   ! Compute instantaneous discharges from point sources
   !
   if (nsrc > 0) then
      do itsrc = itsrclast, ntsrc
         ! Find first point in time series large than t
         if (tsrc(itsrc) > t) then
            do isrc = 1, nsrc
               qtsrc(isrc) = qsrc(isrc, itsrc - 1) + (qsrc(isrc, itsrc) - qsrc(isrc, itsrc - 1)) * (t - tsrc(itsrc - 1)) / (tsrc(itsrc) - tsrc(itsrc - 1))
            enddo
            itsrclast = itsrc - 1
            exit
         endif
      enddo
      !
      !$acc update device(qtsrc)
      !
   endif
   !
   if (ndrn > 0) then
      !
      !$acc serial, present( z_volume, zs, zb, nmindsrc, qtsrc, drainage_type, drainage_params )
      do idrn = 1, ndrn
         !
         jin  = nsrc + idrn * 2 - 1
         jout = nsrc + idrn * 2
         !
         nmin  = nmindsrc(jin)
         nmout = nmindsrc(jout)
         !
         if (nmin > 0 .and. nmout > 0) then
            !
            select case(drainage_type(idrn))
               !
               case(1)
                  !
                  ! Pump
                  !
                  qq = drainage_params(idrn, 1)
                  !
               case(2)
                  !
                  ! Culvert
                  !
                  if (zs(nmin)>zs(nmout)) then
                     !
                     qq  = drainage_params(idrn, 1) * sqrt(zs(nmin) - zs(nmout))
                     !
                  else
                     !
                     qq  = -drainage_params(idrn, 1) * sqrt(zs(nmout) - zs(nmin))
                     !
                  endif
                  !
               case(3)
                  !
                  ! Check valve (same as culvert, but only works in one direction)
                  !
                  if (zs(nmin) > zs(nmout)) then
                     !
                     qq = drainage_params(idrn, 1) * sqrt(zs(nmin) - zs(nmout))
                     !
                  else
                     !
                     qq = -drainage_params(idrn, 1) * sqrt(zs(nmout) - zs(nmin))
                     !
                  endif
                  !
                  ! Make sure it can only flow from intake to outfall point
                  !
                  qq = max(qq, 0.0)
                  !
               case(4)
                  !
                  ! Controlled gate. Gate opens when water level at intake point is between zmin and zmax.
                  !
                  wdt   = drainage_params(idrn, 1)                        ! width
                  zsill = drainage_params(idrn, 2)                        ! sill elevation
                  mng   = drainage_params(idrn, 3)                        ! Manning's n
                  zmin  = drainage_params(idrn, 4)                        ! min water level for open
                  zmax  = drainage_params(idrn, 5)                        ! max water level for open
                  tcls  = drainage_params(idrn, 6)                        ! closing time (seconds)
                  !
                  dzds = (zs(nmout) - zs(nmin)) / drainage_distance(idrn) ! water level slope
                  frac = drainage_fraction_open(idrn)                     ! fraction open (from previous time step)
                  hgate = max(max(zs(nmin), zs(nmout)) - zsill, 0.0)      ! water depth
                  dfrac = dt / tcls                                       ! change in fraction open per time step
                  !
                  qq0 = -qtsrc(jin) / (wdt * max(frac, 0.001))            ! discharge (in m2/s) from previous time step, excluding fraction open
                  !
                  ! Get status of gate
                  !
                  if (drainage_status(idrn) == 0) then
                     !
                     ! Gate fully closed
                     !
                     if (zs(nmin) > zmin .and. zs(nmin) < zmax) then
                        !
                        ! Water level is in allowable range, so need to open the gate
                        !
                        drainage_status(idrn) = 3
                        !
                        write(logstr,'(a,i0,a,f0.1)')'INFO Gates - Opening structure ',idrn,' at t= ',t
                        call write_log(logstr, 0)                        
                        !
                     endif
                     !
                  elseif (drainage_status(idrn) == 1) then
                     !
                     ! Gate fully open
                     !
                     if (zs(nmin) <= zmin .or. zs(nmin) >= zmax) then
                        !
                        ! Water level is NOT in allowable range, so need to close the gate
                        !
                        drainage_status(idrn) = 2
                        !
                        write(logstr,'(a,i0,a,f0.1)')'INFO Gates - Closing structure ',idrn,' at t= ',t
                        call write_log(logstr, 0)                        
                        !                        
                     endif
                     !
                  endif                        
                  !
                  ! Update fraction open
                  !
                  if (drainage_status(idrn) == 2) then
                     !
                     ! Gate is closing
                     !
                     frac = frac - dfrac
                     !
                     if (frac < 0.0) then
                        !
                        ! Gate is now fully closed
                        !
                        frac = 0.0
                        drainage_status(idrn) = 0
                        !
                     endif
                     !
                  elseif (drainage_status(idrn) == 3) then
                     !
                     ! Gate is opening
                     !
                     frac = frac + dfrac
                     !
                     if (frac > 1.0) then
                        !
                        ! Gate is now fully open
                        !
                        frac = 1.0
                        drainage_status(idrn) = 1
                        !
                     endif
                     !
                  endif
                  !
                  drainage_fraction_open(idrn) = frac
                  !
                  ! Use Bates et al. (2010) formulation to include inertia effects
                  !
                  qq = (qq0 - g * hgate * dzds * dt) / (1.0 + g * mng**2 * dt * abs(qq0) / hgate**(7.0 / 3.0))
                  !
                  ! Multiply with width and fraction open to get discharge in m3/s
                  !
                  qq = qq * wdt * frac
                  !
               case(5)
                  !
                  ! Controlled gate. Gate opens and closes at set user input times (only, and once), still using closing time.
                  !
                  wdt   = drainage_params(idrn, 1)                        ! width
                  zsill = drainage_params(idrn, 2)                        ! sill elevation
                  mng   = drainage_params(idrn, 3)                        ! Manning's n
                  tclose = drainage_params(idrn, 4)                       ! time wrt tref for closing gate
                  topen  = drainage_params(idrn, 5)                       ! time wrt tref for opening gate
                  tcls  = drainage_params(idrn, 6)                        ! closing time (seconds)
                  !
                  dzds = (zs(nmout) - zs(nmin)) / drainage_distance(idrn) ! water level slope
                  frac = drainage_fraction_open(idrn)                     ! fraction open (from previous time step)
                  hgate = max(max(zs(nmin), zs(nmout)) - zsill, 0.0)      ! water depth
                  dfrac = dt / tcls                                       ! change in fraction open per time step
                  !
                  qq0 = -qtsrc(jin) / (wdt * max(frac, 0.001))            ! discharge (in m2/s) from previous time step, excluding fraction open
                  !
                  ! Get status of gate
                  !
                  if (drainage_status(idrn) == 0) then
                     !
                     ! Gate fully closed
                     !
                     if (t >= topen) then
                        !
                        ! Time has passed 'topen', so need to open the gate
                        !
                        drainage_status(idrn) = 3
                        !
                        write(logstr,'(a,i0,a,f0.1)')'INFO Gates - Opening structure ',idrn,' at t= ',t
                        call write_log(logstr, 0)                        
                        !
                     endif
                     !
                  elseif (drainage_status(idrn) == 1) then
                     !
                     ! Gate fully open
                     !
                     if (t >= tclose .and. t < topen) then
                        !
                        ! Time has passed 'tclose', so need to close the gate
                        !
                        drainage_status(idrn) = 2
                        !
                        write(logstr,'(a,i0,a,f0.1)')'INFO Gates - Closing structure ',idrn,' at t= ',t
                        call write_log(logstr, 0)                        
                        !                        
                     endif
                     !
                  endif                        
                  !
                  ! Update fraction open
                  !
                  if (drainage_status(idrn) == 2) then
                     !
                     ! Gate is closing
                     !
                     frac = frac - dfrac
                     !
                     if (frac < 0.0) then
                        !
                        ! Gate is now fully closed
                        !
                        frac = 0.0
                        drainage_status(idrn) = 0
                        !
                     endif
                     !
                  elseif (drainage_status(idrn) == 3) then
                     !
                     ! Gate is opening
                     !
                     frac = frac + dfrac
                     !
                     if (frac > 1.0) then
                        !
                        ! Gate is now fully open
                        !
                        frac = 1.0
                        drainage_status(idrn) = 1
                        !
                     endif
                     !
                  endif
                  !
                  drainage_fraction_open(idrn) = frac
                  !
                  ! Use Bates et al. (2010) formulation to include inertia effects
                  !
                  qq = (qq0 - g * hgate * dzds * dt) / (1.0 + g * mng**2 * dt * abs(qq0) / hgate**(7.0 / 3.0))
                  !
                  ! Multiply with width and fraction open to get discharge in m3/s
                  !
                  qq = qq * wdt * frac
                  ! 
               case(6)
                  !
                  ! Dike breaching based on Verheij (2003) and discharge through breach based on submerge and free flow equations
                  !   
                  z_crest  = drainage_params(idrn, 1)              ! initial crest level
                  tbreach   = drainage_params(idrn, 2)             ! time of breaching
                  z_min  = drainage_params(idrn, 3)                ! lowest elevation of breach
                  B0   = drainage_params(idrn, 4)                  ! initial breach width  
                  t_0  = drainage_params(idrn, 5)                  ! time to reach lowest breach elevation
                  dike_core = drainage_params(idrn, 6)             ! material of dike core (1 = sand and 2 = clay)
                  !
				  t_phase1 = tbreach + t_0
				  m_afvoercoeff = 1.0   ! afvoercoefficient
                  !
                  B_old = breach_width(idrn)
				  if (dike_core == 1.0) then
				    !
					! dike core made of sand
					!
					f1 = 1.3
					f2 = 0.04
					uc = 0.2
				  elseif (dike_core == 2.0) then
				    !
					! dike core made of clay
					!
					f1 = 1.3
					f2 = 0.04
					uc = 0.5        
                  endif
                  
                  !
                  ! Updating dike dimensions (crest height and breach width)
                  !
                  if (t >= tbreach) then
                    if (t < t_phase1) then
                        !
                        ! Start of phase 1: lowering of the crest
                        !
						breach_width(idrn) = B0 ! no widening of the breach yet
						Z = z_crest - (z_crest - z_min)*(t-tbreach)/t_0 ! lowering of the crest lineair with time
                        breach_level_gather(ndrn) = Z
                    elseif (t >= t_phase1) then
                        !
						! Start of phase 2: widening, Once phase 2 begins, crest is at minimum
						!
                        ! 
                        Z = z_min
                        breach_level_gather(ndrn) = Z

                        ! Choose downstream level per your earlier logic
                        if (zs(nmout) > z_min) then
                            H = zs(nmin) - zs(nmout)
                        else
                            H = zs(nmin) - z_min
                        endif

                        ! Prevent negative head (no widening if no driving head)
                        H = MAX(H, 0.0)
                        

                        ! Convert time since phase2 start to hours if your formulation expects hours
                        ! Your earlier (15) used /3600, so keep consistency here:
                        tau_hr = (t - t_phase1) / 3600.0
                        dt_hr  = dt / 3600.0

                        ! Denominator term: 1 + (f2*g/uc) * (t_i - t0)
                        denom = 1.0 + (f2 * g / uc) * tau_hr
                        denom = MAX(denom, 1.0e-12)   ! safety

                        ! dB/dt at time t_i  [units: m/hour if dt_hr used]
                        dBdt = (f1 * f2 / LOG(10.0)) * ( (g * H)**1.5 ) / (uc*uc) * (1.0 / denom)

                        ! No negative widening rate
                        dBdt = MAX(dBdt, 0.0)

                        ! update width
                        breach_width(idrn) = B_old + dBdt * dt_hr
                    endif
                    
                  else
                      ! Before breaching time, no changes
                      breach_level_gather(ndrn) = z_crest
                      breach_width(idrn) = 0.0
                  endif
                  
                  !
                  ! Now that the breaching geometry is updated, compute discharge through the breach
                  !
                  if (t >= tbreach) then
                      if (breach_level_gather(ndrn) > MAX(zs(nmin), zs(nmout))) then
                            !
                            ! Dike crest higher than out- and inside water level, so no flow
                            !
                            qq = 0.0
                      elseif (min(zs(nmin),zs(nmout)) > (2.0/3.0)*max(zs(nmin),zs(nmout))) then
                            !
                            ! Fully submerged flow
                            !
                            h_breach = max(max(zs(nmin),zs(nmout))- breach_level_gather(ndrn), 0.0)
				            qq = m_afvoercoeff * breach_width(idrn)  * h_breach * sqrt(2.0 * 9.81 * (max(MAX(zs(nmin), zs(nmout))-MIN(zs(nmin), zs(nmout)),0.0))) 
                            if (zs(nmout)>zs(nmin)) then
                                qq = -qq ! return flow
                            end if
                      else
                            !
                            ! Free flow
                            !
                            h_breach = max(max(zs(nmin),zs(nmout))- breach_level_gather(ndrn), 0.0)
                            qq = 1.71 * breach_width(idrn) * sqrt(9.81) * (h_breach)**1.5 
                            if (zs(nmout)>zs(nmin)) then
                              qq = -qq ! return flow
                            end if
                      endif
                  else
                      !
                      ! No discharge through dike if t<tbreach 
			          !
                      qq = 0.0 
                  endif

                   ! ---- write discharge to log ----

               case(7)
                  !
                   ! Dike breaching based on Verheij (2003) but with Bates discharge formulation
                  !   
                  z_crest  = drainage_params(idrn, 1)              ! initial crest level
                  tbreach   = drainage_params(idrn, 2)             ! time of breaching
                  z_min  = drainage_params(idrn, 3)                ! lowest elevation of breach
                  B0   = drainage_params(idrn, 4)                  ! initial breach width  
                  t_0  = drainage_params(idrn, 5)                  ! time to reach lowest breach elevation
                  dike_core = drainage_params(idrn, 6)             ! material of dike core (1 = sand and 2 = clay)
                  
                  mng = 0.03 ! Manning for the Bates formula, could be adjusted based on dike material or as input 
                  !
				  t_phase1 = tbreach + t_0
				  m_afvoercoeff = 1.0   ! afvoercoefficient, 
                  !
                  B_old = breach_width(idrn)
                  qq0 = -qtsrc(jin) / (breach_width(idrn))            ! discharge (in m2/s) from previous time step

                                    
				  if (dike_core == 1.0) then
				    !
					! dike core made of sand
					!
					f1 = 1.3
					f2 = 0.04
					uc = 0.2
				  elseif (dike_core == 2.0) then
				    !
					! dike core made of clay
					!
					f1 = 1.3
					f2 = 0.04
					uc = 0.5        
                  endif
                  
                  !
                  ! Updating dike dimensions (crest height and breach width)
                  !
                  if (t >= tbreach) then
                    if (t < t_phase1) then
                        !
                        ! Start of phase 1: lowering of the crest
                        !
						breach_width(idrn) = B0 ! no widening of the breach yet
						Z = z_crest - (z_crest - z_min)*(t-tbreach)/t_0 ! lowering of the crest lineair with time
                        breach_level_gather(ndrn) = Z
                    elseif (t >= t_phase1) then
                        ! 
                        ! Start of phase 2: widening, Once phase 2 begins, crest is at minimum
						!
                        Z = z_min
                        breach_level_gather(ndrn) = Z

                        ! Choose downstream level per your earlier logic
                        if (zs(nmout) > z_min) then
                            H = zs(nmin) - zs(nmout)
                        else
                            H = zs(nmin) - z_min
                        endif

                        ! Prevent negative head (no widening if no driving head)
                        H = MAX(H, 0.0)
                        

                        ! Convert time since phase2 start to hours if your formulation expects hours
                        ! Your earlier (15) used /3600, so keep consistency here:
                        tau_hr = (t - t_phase1) / 3600.0
                        dt_hr  = dt / 3600.0

                        ! Denominator term: 1 + (f2*g/uc) * (t_i - t0)
                        denom = 1.0 + (f2 * g / uc) * tau_hr
                        denom = MAX(denom, 1.0e-12)   ! safety

                        ! dB/dt at time t_i  [units: m/hour if dt_hr used]
                        dBdt = (f1 * f2 / LOG(10.0)) * ( (g * H)**1.5 ) / (uc*uc) * (1.0 / denom)

                        ! No negative widening rate
                        dBdt = MAX(dBdt, 0.0)

                        ! update width
                        breach_width(idrn) = B_old + dBdt * dt_hr

                    endif
                    
                    !
                    ! Now that the breaching geometry is updated, compute discharge through the breach
                    !
                    if (breach_level_gather(ndrn) > MAX(zs(nmin), zs(nmout))) then
                        !
                        ! Dike crest higher than out- and inside water level, so no flow
                        !
                        qq = 0.0
                    elseif (min(zs(nmin), zs(nmout)) > 2.0/3.0*max(zs(nmin), zs(nmout))) then
                        !
                        ! Fully submerged flow
                        !
                        dike_width = 10 ! TO DO: check this assumption, now assuming all dikes have a dike width of  10 m instead of using drainage_distance(idrn)
                        dzds = (zs(nmout) - zs(nmin)) / dike_width ! water level slope
                        h_breach = max(max(zs(nmin), zs(nmout)) - breach_level_gather(ndrn), 0.0)   ! water depth
                        !
                        ! Use Bates et al. (2010) formulation to include inertia effects
                        !
                        qq = (qq0 - g * h_breach * dzds * dt) / (1.0 + g * mng**2 * dt * abs(qq0) / h_breach**(7.0 / 3.0))
                        !
                        ! Multiply with width and fraction open to get discharge in m3/s
                        !
                        qq = qq * breach_width(idrn)

                    else
                        !
                        ! Free flow
                        !
                        h_breach = max(max(zs(nmin),zs(nmout))- breach_level_gather(ndrn), 0.0)
                        qq = 1.71 * breach_width(idrn) * sqrt(9.81) * (h_breach)**1.5
                        if (zs(nmout)>zs(nmin)) then
                            qq = -qq
                        end if
                    endif
                    
                  else
                    !
                    ! Before breaching time, no changes
                    !
                    breach_level_gather(ndrn) = z_crest
                    breach_width(idrn) = 0.0
                    qq = 0.0
                  endif
                   ! ---- write discharge to log ----
                  
               case(8)
                  !
                   ! Dike breaching based on Tadesse et al (2014), discharge through breach based on submerge and free flow equations
                  !
                  z_crest  = drainage_params(idrn, 1)              ! crest level of the dike
                  tbreach   = drainage_params(idrn, 2)             ! time of breaching, 
                  breach_duration   = drainage_params(idrn, 3)     ! breach_duration
                  z_min   = drainage_params(idrn, 4)               ! final_breach_level
                  final_breach_width = drainage_params(idrn, 5)    ! final_breach_width
                  B0 = drainage_params(idrn, 6)                    ! initial breach width
                  
                  
                  t_0  = breach_duration/10.0                      ! Time of phase 1 (lowering of crest), 1/10th of the total breach duration
                  t_phase1 = tbreach + t_0
                  t_end = tbreach + breach_duration
                  
                  !
                  ! Updating dike dimensions (crest height and breach width)
                  !
                  if (t < tbreach ) then
                      !
                      ! Before breaching time, no changes
                      !
                      Z = z_crest
                      breach_level_gather(ndrn) = Z
                      breach_width(idrn) = 0.0
                  elseif (t >= tbreach .and. t < t_phase1) then
                     !
                     ! Start of phase 1: lowering of the crest
                     !
                     Z = z_crest - (z_crest - z_min)*(t-tbreach)/t_0 ! lowering of the crest lineair with time at 1/10th of the breach duration
                     breach_level_gather(ndrn) = Z
                     B = B0 + (final_breach_width - B0)*(t-tbreach)/breach_duration ! widening of the breach
                     breach_width(ndrn) = B
                  elseif (t >= t_phase1 .and. t<=t_end) then
                     B = B0 + (final_breach_width - B0)*(t-tbreach)/breach_duration ! widening of the breach
                     breach_width(ndrn) = B
                     breach_level_gather(ndrn) = Z
                  else
                      ! After breach end time, final dimensions
                      Z = z_min
                      breach_level_gather(ndrn) = Z
                      breach_width(idrn) = final_breach_width
                  endif
                  
                  !
                  ! Now that the breaching geometry is updated, compute discharge through the breach
                  !
                  if (t >= tbreach) then
                      if (breach_level_gather(ndrn) > MAX(zs(nmin), zs(nmout))) then
                            !
                            ! Dike crest higher than out- and inside water level, so no flow
                            !
                            qq = 0.0
                      elseif (min(zs(nmin),zs(nmout)) > (2.0/3.0)*max(zs(nmin),zs(nmout))) then
                            !
                            ! Fully submerged flow
                            !
                            h_breach = max(max(zs(nmin),zs(nmout))- breach_level_gather(ndrn), 0.0)
				            qq = m_afvoercoeff * breach_width(idrn)  * h_breach * sqrt(2.0 * 9.81 * (max(MAX(zs(nmin), zs(nmout))-MIN(zs(nmin), zs(nmout)),0.0)))  
                            if (zs(nmout)>zs(nmin)) then
                              qq = -qq! return flow
                            end if
                      else
                            !
                            ! Free flow
                            !
                            h_breach = max(max(zs(nmin),zs(nmout))- breach_level_gather(ndrn), 0.0)
                            qq = 1.71 * breach_width(idrn) * sqrt(9.81) * (h_breach)**1.5
                            if (zs(nmout)>zs(nmin)) then
                              qq = -qq ! return flow
                            end if
                      endif
                  else
                      ! No discharge through dike if t<tbreach and water level is below crest level
			          qq = 0.0 
                  endif
                   
               case(9)                 
                  ! Breaching based on BRES model (Visser 1998)
                  outside_level= drainage_params(idrn, 1)       ! outside groundlevel, called Zw in Visser
                  polder_level = drainage_params(idrn, 2)       ! polder groundlevel, called Zp in Visser
                  crest_level = drainage_params(idrn, 3)        ! crest level of the dike
                  crest_width   = drainage_params(idrn, 4)      ! Crest width of the dike
                  breach_level = drainage_params(idrn, 5)       ! Initial breach level
                  t0_Visser = drainage_params(idrn, 6)          ! Time of initiation of breaching
                  
                  !
                  outside_water_level = zs(nmin)
                  polder_water_level = zs(nmout)
                  
                  !
                  ! Input parameters 
                  ! 
                  breach_bottom = 1.0 ! initial breach width
                  beta0 = 18.0  * 3.141592653589793 / 180.0 ! Inner slope angle of the dike (radians)
                  alpha =  32.0 * 3.141592653589793 / 180.0 ! Outside slope angle of the dike (radians)      
                  beta1 =  40 * 3.141592653589793 / 180.0 ! critical inside slope angle of the dike (radians)
                  d50 = 220.0/10**6 ! m
                  d90 = 350.0/10**6 ! m
                  Cf = 0.025 ! Friction coefficient
                  Kappa = 0.41 ! Von Karman constant
                  water_temperature = 17.0 ! degrees Celsius
                  water_density = 1025.0 ! kg/m3
                  p = 0.4 ! porosity
                  phi = 32.0 * 3.141592653589793 / 180.0! angle of repose (radian)
                  gamma0 = 32* 3.141592653589793 / 180.0 ! Breach slope (radian)
                  gamma1 = 60* 3.141592653589793 / 180.0 ! Critical breach slope (radian)
                  sediment_density = 2650
                  water_density = 1025
                  delta = (sediment_density - water_density) / water_density
                  dyn_viscosity = dynamic_viscosity(water_temperature, water_density)
                  sediment_fall_velocity = calc_sediment_fall_velocity(d50, delta, dyn_viscosity)
                  dstar = dimensionless_diameter(d50, delta, dyn_viscosity) 
                  formula = 'Shields'
                  theta_crit = critical_shields_parameter(dstar, formula)
                  ni = 0.48 ! Sheared porosity
                  k = 0.000375 ! permeability 
                  
                  !
                  ! Calculating breaching geometry for phase 1 and 2, only needs to be done once
                  !
                  if (t >= t0_Visser) then
                      if (running_Visser_phase1(ndrn) == 0) then
                            !
                            ! Phase 1, calculate only once
                            !
                            W = breach_crest_length(crest_width, crest_level, breach_level, alpha, beta0)
                            results_t1 = stage_1(t0_Visser, breach_bottom, polder_level, polder_water_level, breach_level, beta1, beta0, outside_water_level, gamma0, alpha, W, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity)
                            discharge_t1(ndrn) = results_t1%discharge
                            t1_Visser(ndrn) = results_t1%t1
                            breach_bottom_Visser(ndrn) = breach_bottom
                            breach_level_gather(ndrn) = breach_level
                            breach_width_waterline_Visser(ndrn) = results_t1%breach_width_waterline
                            breach_width(ndrn) = results_t1%breach_width_total
                            gamma0_Visser(ndrn) = gamma0
                            running_Visser_phase1(ndrn) = 1
                      end if
                      if (t >= t1_Visser(ndrn)) then
                          if (running_Visser_phase2(ndrn) == 0) then
                                !
                                ! Phase 2, calculate only once
                                !
                                W = breach_crest_length(crest_width, crest_level, breach_level, alpha, beta0) !! TO DO: change beta0 to beta1 for phase 2
                                results_t2 = stage_2(t1_Visser(ndrn),polder_level, breach_bottom, breach_level, polder_water_level,beta1, outside_water_level, gamma0, alpha, W, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity)
                                discharge_t2(ndrn) = results_t2%discharge
                                t2_Visser(ndrn) = results_t2%t2
                                breach_bottom_Visser(ndrn) = breach_bottom
                                breach_level_gather(ndrn) = breach_level
                                breach_width_waterline_Visser(ndrn) = results_t2%breach_width_waterline
                                breach_width(ndrn) = results_t2%breach_width_total
                                gamma0_Visser(ndrn) = gamma0
                                running_Visser_phase2(ndrn) = 1
                          end if
                      end if 
                  else
                      !
                      ! Before breaching time, no changes
                      !                      
                      breach_level_gather(ndrn) = crest_level
                      breach_width(idrn) = 0.0
                  end if 
                  
                  
                  !
                  ! Calculting breaching geometry for phase 3, 4 and 5
                  !
                  if (t > t0_Visser .AND. t> t2_Visser(ndrn) .AND. end_breaching(ndrn) == 0.0) then
                      breach_bottom = breach_bottom_Visser(ndrn)
                      breach_level = breach_level_gather(ndrn)
                      gamma0 = gamma0_Visser(ndrn)
                      breach_width_waterline = breach_width_waterline_Visser(ndrn)
                      breach_width_total = breach_width(ndrn)
                      crit_water_depth = calc_crit_water_depth(outside_water_level, breach_level, breach_bottom, gamma1)
                      if (breach_level - outside_level > 0) then ! if the breach is above the outside groundlevel
                            ! 
                            ! Phase 3
                            !
                            results_t3 = stage_3(dt,breach_width_total, breach_width_waterline, theta_crit, ni, dstar, k, sediment_density, water_density, outside_level, polder_level, breach_bottom, breach_level, polder_water_level,beta1, outside_water_level, gamma0,gamma1, alpha, W, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity)
                            breach_bottom_Visser(ndrn) = results_t3%breach_bottom
                            breach_level_gather(ndrn) = results_t3%breach_level
                            breach_width_waterline_Visser(ndrn) = results_t3%breach_width_waterline 
                            breach_width(ndrn) = results_t3%breach_width_total 
                            gamma0_Visser(ndrn) = results_t3%gamma0
                      
                      else if (polder_water_level - breach_level <= crit_water_depth) then
                            ! 
                            ! Phase 4
                            !                            
                            results_t4 = stage_4(dt, breach_width_total, breach_width_waterline, theta_crit, ni, dstar, k, sediment_density, water_density, outside_level, polder_level, breach_bottom, breach_level, polder_water_level, outside_water_level, gamma1, alpha, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity) 
                            breach_bottom_Visser(ndrn) = results_t4%breach_bottom
                            breach_level_gather(ndrn) = results_t4%breach_level
                            breach_width_waterline_Visser(ndrn) = results_t4%breach_width_waterline 
                            breach_width(ndrn) = results_t4%breach_width_total 
                            gamma0_Visser(ndrn) = gamma0
                      
                      else if (outside_water_level > polder_water_level .AND. outside_water_level>breach_level) then
                            ! 
                            ! Phase 5
                            !                            
                            results_t5 = stage_5(dt, breach_width_total, breach_width_waterline, theta_crit, beta1, ni, dstar, k, sediment_density, water_density, outside_level, polder_level, breach_bottom, breach_level, polder_water_level, outside_water_level, gamma1, alpha, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity) 
                            breach_bottom_Visser(ndrn) = results_t5%breach_bottom
                            breach_level_gather(ndrn) = results_t5%breach_level
                            breach_width_waterline_Visser(ndrn) = results_t5%breach_width_waterline 
                            breach_width(ndrn) = results_t5%breach_width_total 
                            gamma0_Visser(ndrn) = gamma0
                      else
                            end_breaching(ndrn) = 1.0
                      end if
                  end if
                  
                  
                  !
                  ! Now that the breaching geometry is updated, compute discharge through the breach
                  !
                  if (t >= t0_Visser) then
                      if (breach_level_gather(ndrn) > MAX(zs(nmin), zs(nmout))) then
                            !
                            ! Dike crest higher than out- and inside water level, so no flow
                            !
                            qq = 0.0
                      elseif (min(zs(nmin),zs(nmout)) > (2.0/3.0)*max(zs(nmin),zs(nmout))) then
                            !
                            ! Fully submerged flow
                            !
                            h_breach = max(max(zs(nmin),zs(nmout))- breach_level_gather(ndrn), 0.0)
				            qq = m_afvoercoeff * breach_width(idrn)  * h_breach * sqrt(2.0 * 9.81 * (max(MAX(zs(nmin), zs(nmout))-MIN(zs(nmin), zs(nmout)),0.0))) 
                            if (zs(nmout)>zs(nmin)) then
                                qq = -qq ! return flow
                            end if
                      else
                            !
                            ! Free flow
                            !
                            h_breach = max(max(zs(nmin),zs(nmout))- breach_level_gather(ndrn), 0.0)
                            qq = 1.71 * breach_width(idrn) * sqrt(9.81) * (h_breach)**1.5
                            if (zs(nmout)>zs(nmin)) then
                              qq = -qq ! return flow
                            end if
                      endif
                  else
                      !
                      ! No discharge through dike if t<tbreach
                      !
			          qq = 0.0 
                  endif

            end select
            !
            ! Add some relaxation
            ! structure_relax in seconds => gives ratio between new and old discharge (default 10s)
            !   
            qq = 1.0 / (structure_relax / dt) * qq + (1.0 - (1.0 / (structure_relax / dt))) * -qtsrc(jin)
            !   
            ! Limit discharge based on available volume in cell (regular or subgrid)
            !    
            if (subgrid) then
               !
               if (qq > 0.0) then
                  qq = min(qq, max(z_volume(nmin), 0.0) / dt)
               else
                  qq = max(qq, -max(z_volume(nmout), 0.0) / dt)
               endif
               !
            else
               !
               if (qq > 0.0) then
                  qq = min(qq, max((zs(nmin) - zb(nmin)) * cell_area(z_flags_iref(nmin)), 0.0) / dt)
               else
                  qq = max(qq, -max((zs(nmout) - zb(nmout)) * cell_area(z_flags_iref(nmout)), 0.0) / dt)
               endif
               !
            endif
            !      
            qtsrc(jin)  = -qq 
            qtsrc(jout) = qq
            !
         endif
         !
      enddo
      !$acc end serial
      !
   endif
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0 * (count1 - count0) / count_rate
   !
   end subroutine

end module
