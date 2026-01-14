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
         elseif (drainage_type(idrn)==6 .or. drainage_type(idrn)==7 .or. drainage_type(idrn)==8) then
            !
			! Dike breaching
            ! Drainage_type 6: Verheij (6 parameters: afvoercoefficient, material of dike core (1 = sand and 2 = clay), time of breaching, initial crest level, lowest elevation of breach, time to reach lowest breach elevation)
            ! Drainage_type 7: Tadesse (6 parameters: time of breaching, breach duration, final_breach_width, final_breach_level, crest level of the dike, initial breach width)

            !
            allocate(breach_width(ndrn))
            breach_width= 0.0
            allocate(parameters_Visser_phase1(ndrn))
            parameters_Visser_phase1 = 0
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
   real*4           :: breach_duration, final_breach_width, final_breach_level, B_0, t_end, dt_hr, tau_hr, H, dBdt
   real*4           :: stage, beta, beta1, beta0, breach_bottom, breach_level, test_value
   real*4           :: alpha, crest_level, W, polder_level, t0_Visser
   real*4           :: d50, d90, Cf, Kappa, water_temperature, p, phi
   real*4           :: sediment_density, water_density, delta, dyn_viscosity, sediment_fall_velocity
   real*4           :: discharge_t1, t1_Visser, outside_water_level, polder_water_level, gamma0,n_t1
   type(NormalFlow) :: results_t1
   
   
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
               case(7)
                  !
                   ! Dike breaching based on Verheij (2003)
                  !   
                  m_afvoercoeff   = drainage_params(idrn, 1)                  ! afvoercoefficient, 
                  dike_core = drainage_params(idrn, 2)                        ! material of dike core (1 = sand and 2 = clay), 
                  tbreach   = drainage_params(idrn, 3)                        ! time of breaching, 
                  z_crest  = drainage_params(idrn, 4)                         ! initial crest level, 
                  z_min  = drainage_params(idrn, 5)                           ! lowest elevation of breach, 
                  t_0  = drainage_params(idrn, 6)                             ! time to reach lowest breach elevation
                  !
				  t_phase1 = tbreach + t_0
				  B0 = 10.0   ! initial breach width
				  ! Z = z_crest ! Initial crest level of the breach
				  ! B=0.0       ! initially no breach width
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
                  ! Updating dike dimensions (crest height and breach width)
                  if (t >= tbreach .and. t < t_phase1) then
                        !
                        ! Start of phase 1: lowering of the crest
                        !
						breach_width(idrn) = B0 ! no widening of the breach yet
						Z = z_crest - (z_crest - z_min)*(t-tbreach)/t_0 ! lowering of the crest lineair with time
                  elseif (t >= t_phase1) then
						!
						! Start of phase 2: widening
						!
                        ! Once phase 2 begins, crest is at minimum
                        Z = z_min

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
                      
                      
                        ! if (zs(nmout)>z_min) then
                            ! Downstream water level above crest level
						!    dB = f1 * 9.81**0.5 * (zs(nmin) - zs(nmout))**1.5/(log(10.0)*uc) * LOG( 1.0+f2 * 9.81/uc * (t-t_phase1)/3600) ! widening of the breach
                        !    dB_dt = 
                        !else
                            ! Downstream water level not above crest level
                        !    dB = f1 * 9.81**0.5 * (zs(nmin) - z_min)**1.5/(log(10.0)*uc) * LOG( 1.0+f2 * 9.81/uc * (t-t_phase1)/3600) ! widening of the breach
                        !endif
                        !dB = MAX(dB, 0.0) ! breach width cannot decrease
                        !breach_width(idrn) = MAX(B_old, B0 + dB) 
                   endif
                  
				   if (t >= tbreach .and. zs(nmin)>Z) then
                    !   After breaching time and water level at intake point is above crest level
                    if (zs(nmout) > (2.0/3.0)*zs(nmin)) then
                        ! Fully submerged flow
				        qq = m_afvoercoeff * breach_width(idrn)  * (zs(nmin) - Z) * sqrt(2.0 * 9.81 * (zs(nmin) - Z)) 
                        write(logstr,'(A,F0.1,A,ES14.6,A,F0.2)')'Breach at time: ', t,  'submerged discharge flow: ',qq, ' Width of Breach: ', breach_width(idrn)
                        call write_log(logstr, 0) 
                    else
                        ! Free flow
                        qq = 1.71 * breach_width(idrn) * sqrt(9.81) * (zs(nmin) - Z)**1.5
                        write(logstr,'(A,F0.1,A,ES14.6,A,F0.2)')'Breach at time: ', t,  'free discharge flow: ',qq, ' Width of Breach: ', breach_width(idrn)
                         call write_log(logstr, 0)    
                    endif
                   else
                    ! No discharge through dike if t<tbreach and water level is below crest level
			        qq = 0.0 
                   endif
                   ! ---- write discharge to log ----

               case(8)
                  !
                  ! Dike breaching based on Tadesse et al (2014)
                  !
                  tbreach   = drainage_params(idrn, 1)                        ! time of breaching, 
                  breach_duration   = drainage_params(idrn, 2)                  ! breach_duration
                  final_breach_width = drainage_params(idrn, 3)                 ! final_breach_width
                  final_breach_level   = drainage_params(idrn, 4)               ! final_breach_level
                  z_crest  = drainage_params(idrn, 5)                           ! crest level of the dike
                  B_0 = drainage_params(idrn, 6)                                 ! initial breach width
                  
                  
                  t_0  = breach_duration/10.0                       ! Time of phase 1 (lowering of crest)

                  !
                  t_phase1 = tbreach + t_0
                  t_end = tbreach + breach_duration
                  B0 = 10.0   ! initial breach width
                  Z = z_crest ! Initial crest level of the breach
                  B=0.0       ! initially no breach width
                  !
                  ! Updating dike dimensions (crest height and breach width)
                  if (t >= tbreach .and. t < t_phase1) then
                     !
                     ! Start of phase 1: lowering of the crest
                     !
                     Z = z_crest - (z_crest - final_breach_level)*(t-tbreach)/t_0 ! lowering of the crest lineair with time at 1/10th of the breach duration
                     B = B_0 + (final_breach_width - B_0)*(t-tbreach)/breach_duration ! widening of the breach
                  elseif (t >= t_phase1 .and. t<=t_end) then
                     B = B_0 + (final_breach_width - B_0)*(t-tbreach)/breach_duration ! widening of the breach
                  endif

                  if (t >= tbreach .and. zs(nmin)>Z) then
                       !   After breaching time and water level at intake point is above crest level
                       if (zs(nmout) > 2/3*zs(nmin)) then
                          ! Fully submerged flow
                          qq = m_afvoercoeff * B * (zs(nmin) - Z) * sqrt(2.0 * 9.81 * (zs(nmin) - Z)) 
                       else
                          ! Free flow
                          qq = 1.71 * B * sqrt(9.81) * (zs(nmin) - Z)**1.5
                       endif
                  else
                      ! No discharge through dike if t<tbreach and water level is below crest level
                      qq = 0.0 
                  endif
                   
               case(6)

                   
                   
                  ! Breaching based on BRES model (Visser 1998)
                  beta0   = drainage_params(idrn, 1)      ! Inner slope angle of the dike (radians)
                  alpha   = drainage_params(idrn, 2)      ! Outside slope angle of the dike (radians)
                  crest_level = drainage_params(idrn, 3)                 
                  W   = drainage_params(idrn, 4)          ! Length of the weir to calculate the discharge coefficient
                  polder_level = drainage_params(idrn, 5)
                  t0_Visser = drainage_params(idrn, 6)          ! Time of initiation of breaching
                  
                  !! Assume starting breach dimensions
                  breach_bottom = 1 ! breach width
                  breach_level = 0.75 !!!!!! TO DO make it crest level dependent starting breach level
                  
                  !! Assume breach parameters
                  beta1 =  40 * 3.141592653589793 / 180.0 ! critical inside slope angle of the dike (radians)
                  d50 = 220 ! um
                  d90 = 350 ! um
                  Cf = 0.025 ! Friction coefficient
                  Kappa = 0.41 ! Von Karman constant
                  water_temperature = 17.0 ! degrees Celsius
                  water_density = 1025.0 ! kg/m3
                  p = 0.4 ! porosity
                  phi = 32.0 * 3.141592653589793 / 180.0! angle of repose (radian)
                  gamma0 = 32* 3.141592653589793 / 180.0 ! Breach slope (radian)
                  beta0 = beta0 * 3.141592653589793 / 180.0 ! radians
                  alpha =  alpha * 3.141592653589793 / 180.0 ! radians
                  
                  sediment_density = 2650
                  water_density = 1025
                  delta = (sediment_density - water_density) / water_density
                  dyn_viscosity = dynamic_viscosity(water_temperature, water_density)
                  sediment_fall_velocity = calc_sediment_fall_velocity(d50, delta, dyn_viscosity)
                  
                  outside_water_level = zs(nmin)
                  polder_water_level = zs(nmout)
                  
                  if (t >= t0_Visser) then
                      
                     if (parameters_Visser_phase1(ndrn) == 0) then
                        call write_log('------------ calculating phase 1 ------------', 1)
                        write(logstr,'(a,f6.3)')'Outside water level: ', zs(nmin)
                        call write_log(logstr,1)
                        results_t1 = stage_1(t0_Visser, breach_bottom, polder_level, polder_water_level, breach_level, beta1, beta0, outside_water_level, gamma0, alpha, W, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity)
                        discharge_t1 = results_t1%discharge
                        t1_Visser = results_t1%t1
                        
                        
                        write(logstr,'(a,f6.5,a,f6.5, a)')'Calculated t1 (', t1_Visser, ') with a discharge of (', discharge_t1, ' m3/s)'
                        call write_log(logstr,1)
                        parameters_Visser_phase1(ndrn) = 1
                     end if
                     !
                     ! Stage 0: No breaching
                     if (t < t1_Visser) then
                        !
                        qq = discharge_t1
                        !
                     !
                  !elseif (t >= t1 .and. t < t2) then
                     !
                     ! Stage 1: Initiation of breaching
                     !qq = testing_function(beta)
                  !elseif (t >= t2 .and. t < t3) then
                     !
                     ! Stage 2: Breach growth
                     !qq = testing_function(beta)
                     !
                  !elseif (t >= t3) then
                     !
                     ! Stage 3: Breach stabilization
                     !qq = testing_function(beta)
                     !
                      else
                         !
                         qq = 0.0
                      end if
                     !
                  end if
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
   
   function dynamic_viscosity(water_temperature, water_density) result (nu)
   !Computes the dynamic viscosity of water in Ns/m^2.
   implicit none
   REAL*4:: water_temperature, water_density
   REAL*4:: nu
   REAL*4:: temp_min, temp_max

   temp_min = 0
   temp_max = 40

   if (0 <= water_temperature <= 20) then
        nu = ((water_density + 1505) / (2500 * water_density) * 10**(13 / (10 - 0.081 * (20 - water_temperature)) - 4.3))

   elseif (20 < water_temperature <= 40) then
        nu = ((water_density + 1505) / (2500 * water_density) * 10**((1.33 * (20 - water_temperature) / (water_temperature + 104)) - 3.0))
   else
        error stop "Water temperature is out of range! The valid range is 0 to 40 degrees Celsius"
   end if 
   end function
   
   function calc_sediment_fall_velocity(d50, rel_density, dyn_viscosity) result (ws)
   !Computes the sediment fall velocity in m/s.
   implicit none
   REAL*4:: d50, rel_density, dyn_viscosity
   REAL*4:: ws
   REAL*4:: g, d50_min
   g = 9.81  ! m/s2, gravitational acceleration
   d50_min = 1e-6

   if (1e-6 <= d50 <= 1e-4) then
        ws = rel_density * g * d50**2 / (18 * dyn_viscosity)

   elseif (1e-4 < d50 < 1e-3) then
        ws = ((10 * dyn_viscosity) / d50* (((1 + 0.01 * rel_density * g * d50**3 * dyn_viscosity**(-2))**0.5) - 1))

   elseif (d50 >= 1e-3) then
        ws = 1.1 * (rel_density * g * d50)**0.5
   else
        error stop "The value for d50 is out of range!"
   end if
    
   end function
   
   function testing_function(random_test_value) result(test_value)
   implicit none
   real*4:: random_test_value
   real*4:: test_value
   
   test_value= 10.0*random_test_value
   
   end function 
   
   function calc_crit_water_depth(water_level, breach_level, breach_width_bottom, slope_angle_breach) result (crit_water_depth)
   implicit none
   real*4:: water_level, breach_level, breach_width_bottom, slope_angle_breach
   real*4:: crit_water_depth, crit_water_depth_iter, breach_width_avg_water_depth, breach_width_waterline
   real*4, parameter :: pi =  3.141592653589793  ! pi in single precision
       
   !! Not sure to add this below because there might be cases where the water level is below the breach level due to the tides
   !if (water_level <= breach_level) then
   !     error stop "Water level should be above breach level!"
   !end if
   
   ! Estimate critical water depth from frictionless flow (2/3=0.7 of water depth)
   crit_water_depth = 0.7 * (water_level - breach_level)

   ! Define an iteration variable to start while loop
   crit_water_depth_iter = 0.9 * crit_water_depth
   
   
   do while (abs((crit_water_depth - crit_water_depth_iter) / crit_water_depth) > 0.01)
    ! Set iteration variable to previous estimate
    crit_water_depth_iter = crit_water_depth
    
    if (breach_width_bottom <= 0) then
        error stop "Breach bottom width: {breach_width_bottom:.2f} is non-positive!"
    end if

    if (crit_water_depth < 0) then
        !error stop "Water depth ", crit_water_depth, " is negative! because outside water level is ", outside_water_level, " and breach level is ", breach_level
        write(logstr,'(a,f6.2,a,f6.2,a,f6.2, a)')"Warning: Critical water depth is negative (",crit_water_depth,") because outside water level (", water_level ,") is below breach level (",breach_level,"). Setting critical water depth to zero."
        error stop logstr
    end if

    if (slope_angle_breach > pi/2 .or. slope_angle_breach <= 0) then
        error stop "Breach slope angle is not within range, the range is: 0 < radians <= pi/2"
    end if

    ! Calculate widths
    breach_width_avg_water_depth = breach_width_bottom + (crit_water_depth / tan(slope_angle_breach))
    breach_width_waterline = breach_width_bottom + (2 * crit_water_depth / tan(slope_angle_breach))

    ! Update estimate of critical water depth
    crit_water_depth = (2.0 / (2.0 + breach_width_avg_water_depth / breach_width_waterline)) * &
                     (water_level - breach_level)
   end do
   end function 
   
   function calc_normal_flow_velocity(hydraulic_radius, slope_angle_dike, friction_coeff) result (normal_flow_velocity)
   ! Calculates the normal flow velocity on the dike inner slope
   implicit none
   real*4:: hydraulic_radius, slope_angle_dike, friction_coeff, normal_flow_velocity
   real*4 :: g
   g = 9.81
   normal_flow_velocity = (g * hydraulic_radius * sin(slope_angle_dike) / friction_coeff)**0.5
   end function
   
   function calc_crit_flow_velocity(crit_water_depth, breach_width_avg_water_depth, breach_width_waterline) result (crit_flow_velocity)
   !Calculates the critical flow velocity in the breach based on broad-crested weir flow
   implicit none
   real*4:: crit_water_depth, breach_width_avg_water_depth, breach_width_waterline, crit_flow_velocity
   real*4:: g
   g = 9.81
   crit_flow_velocity = (g * crit_water_depth * breach_width_avg_water_depth / breach_width_waterline)**0.5
   end function
   
   
   function calc_discharge(breach_width_avg_water_depth, water_depth, flow_velocity, discharge_coeff) result (discharge)
   !Calculates the discharge through the breach in stages I to V.
   !In stages I to IV, assumes critical flow over a weir.
   !In stage V, assumes free gravity flow.
   implicit none
   real*4:: breach_width_avg_water_depth, water_depth, flow_velocity, discharge_coeff, discharge
     
   discharge = discharge_coeff * breach_width_avg_water_depth * water_depth * flow_velocity
   end function
   
   function normal_flow_conditions(water_level, breach_level, breach_width_bottom, discharge, slope_angle_breach, slope_angle_dike, delta, d50, d90, friction_coeff, kappa) result (res)
   ! Calculates the normal flow conditions on the slope
   implicit none
   real*4:: water_level, breach_level, breach_width_bottom, discharge, slope_angle_breach, slope_angle_dike, delta, d50, d90, friction_coeff, kappa,g,normal_water_depth, normal_flow_velocity
   REAL*4:: breach_width_avg_water_depth, normal_water_depth_iter, hydraulic_radius, theta
   
    
   !type :: NormalFlow
   !  real*4 :: normal_water_depth
   !  real*4 :: normal_flow_velocity
   !  real*4 :: friction_coeff
   !end type NormalFlow

   ! Result
   type(NormalFlow) :: res


   g = 9.81
   ! Estimate the normal water depth
   normal_water_depth = 0.3 * (water_level - breach_level)

   breach_width_avg_water_depth = breach_width_bottom + (normal_water_depth / tan(slope_angle_breach))
   normal_water_depth = ((friction_coeff * (discharge / breach_width_avg_water_depth)**2) / (g * sin(slope_angle_dike)))**(1.0/3.0)

   ! Define an iteration variable to start the while loop
   normal_water_depth_iter = 0.9 * normal_water_depth

   do while (abs((normal_water_depth - normal_water_depth_iter) / normal_water_depth) > 0.01)
        normal_water_depth_iter = normal_water_depth
        breach_width_avg_water_depth = breach_width_bottom + (normal_water_depth / tan(slope_angle_breach)) 
        hydraulic_radius = calc_hydraulic_radius(normal_water_depth, breach_width_bottom, breach_width_avg_water_depth, slope_angle_breach)
        normal_flow_velocity = calc_normal_flow_velocity(hydraulic_radius, slope_angle_dike, friction_coeff)
        normal_water_depth = discharge / (normal_flow_velocity * breach_width_avg_water_depth)
        theta = calc_shields_parameter(normal_flow_velocity, delta, d50, friction_coeff)
        friction_coeff = calc_friction_coefficient(theta, hydraulic_radius, d90, kappa)
   end do
   
   res%normal_water_depth = normal_water_depth
   res%normal_flow_velocity = normal_flow_velocity
   res%friction_coeff = friction_coeff
   
   end function
   
   function calc_shields_parameter(flow_velocity, delta, d50, friction_coeff) result (theta)
   !Calculates the Shields mobility parameter.
   implicit none
   REAL*4:: flow_velocity, delta, d50, friction_coeff, g, theta
   g = 9.81
   if (friction_coeff <= 0) then
        error stop "Friction coefficient is non-positive!"
   end if

   if (flow_velocity <= 0)  then
        error stop "Flow velocity is non-positive!"
   end if 

   theta = (friction_coeff * flow_velocity**2) / (g * delta * d50)
   end function
   
   function calc_friction_coefficient(theta, hydraulic_radius, d90, kappa) result (friction_coeff)
   !Computes the friction coefficient of the sediment layer and corresponding Shields (mobility) parameter.
   implicit none
   REAL*4:: theta, hydraulic_radius, d90, kappa, friction_coeff, k
   
   if (theta < 0) then
        error stop "Shields parameter (theta) is negative!"
   end if

   if (hydraulic_radius <= 0) then
         error stop "Hydraulic radius is non-positive!"
   end if

   ! A limit for theta to avoid very high k
   theta = min(theta, 300.0)

   ! Compute layer roughness (k)
   if (theta < 1) then
        k = 3.0 * d90
   else
       k = 3.0 * d90 * theta
   end if

   ! New estimate of the friction coefficient
   friction_coeff = kappa**2 / ((log(12 * hydraulic_radius / k))**2.0)
   end function
   
   function calc_hydraulic_radius(water_depth, breach_width_bottom, breach_width_avg_water_depth, slope_angle_breach) result (hydraulic_radius)

   !Computes the hydraulic radius of the trapezoidal cross-section.
   implicit none
   REAL*4:: water_depth, breach_width_bottom, breach_width_avg_water_depth, slope_angle_breach, hydraulic_radius
   REAL*4:: area, wetted_perimeter
   real*4, parameter :: pi = 3.141592653589793  ! pi in single precision
   
   if (water_depth < 0.0) then
        error stop "Water depth is negative!"
   end if

   if (breach_width_bottom < 0) then
        error stop "Breach bottom width is negative!"
   end if

   if (slope_angle_breach > pi/2 .or. slope_angle_breach <= 0) then
        error stop "Breach slope angle is not within range, the range is: 0 < radians <= pi/2"
   end if 

   area = breach_width_avg_water_depth * water_depth
   wetted_perimeter = breach_width_bottom + (2 * water_depth / sin(slope_angle_breach))
   hydraulic_radius = area / wetted_perimeter
   end function
   
   function calc_breach_width_avg_water_depth(breach_width_bottom, water_depth, slope_angle_breach) result (breach_width_avg_water_depth)
   !Calculates the breach width at average water depth
   implicit none
   REAL*4:: breach_width_bottom, water_depth, slope_angle_breach, breach_width_avg_water_depth
   real*4, parameter :: pi =  3.141592653589793 ! pi in single precision
   if (breach_width_bottom <= 0) then
        error stop "Breach bottom width is non-positive!"
   end if

   if (water_depth < 0) then
        error stop "Water depth is negative!"
   end if

   if (slope_angle_breach > pi/2 .or. slope_angle_breach <= 0) then
        error stop "Breach slope angle is not within range, the range is: 0 < radians <= pi/2"
   end if

   breach_width_avg_water_depth = breach_width_bottom + (water_depth / tan(slope_angle_breach))
   end function
   
   function calc_breach_width_waterline(breach_width_bottom, water_depth, slope_angle_breach) result(breach_width_waterline)
   ! Calculates the breach width at the waterline
   implicit none
   REAL*4:: breach_width_bottom, water_depth, slope_angle_breach, breach_width_waterline
   real*4, parameter :: pi =  3.141592653589793  ! pi in single precision

    if (breach_width_bottom <= 0) then
         error stop "Breach bottom width is non-positive!"
   end if

   if (water_depth < 0) then
        error stop "Water depth is negative!"
   end if

   if (slope_angle_breach > pi/2 .or. slope_angle_breach <= 0) then
        error stop "Breach slope angle is not within range, the range is: 0 < radians <= pi/2"
   end if

   breach_width_waterline = breach_width_bottom + (2 * water_depth / tan(slope_angle_breach))
   end function
   
   function calc_discharge_coeff(formula, outer_slope_angle, inner_slope_angle, outside_water_level, breach_level, weir_length) result(discharge_coeff)
   
   implicit none
   character(len=* ) :: formula
   real*4 :: outer_slope_angle, inner_slope_angle, outside_water_level, breach_level, weir_length
   real*4 :: outer_slope, inner_slope, c1, c2, Cd, discharge_coeff, zeta, arg
   character(len=200) :: errmsg

   select case (trim(formula))
   case ('Chen2018')
        outer_slope = 1.0 / tan(outer_slope_angle)
        inner_slope = 1.0 / tan(inner_slope_angle)

        if (inner_slope >= 0.0 .and. inner_slope <= 0.8) then
            c1 = (-1.3 * outer_slope + 8.09) * ( 2.862 * inner_slope +  7.658) * 1.0e-3
            c2 = (-8.6 * outer_slope**2 + 7.9 * outer_slope + 493.5) * ( 5.33 * inner_slope + 95.74) * 1.0e-5
        else
            c1 = (-1.3 * outer_slope + 8.09) * (-1.797 * inner_slope + 11.355) * 1.0e-3
            c2 = (-8.6 * outer_slope**2 + 7.9 * outer_slope + 493.5) * (-4.24 * inner_slope + 103.28) * 1.0e-5
        end if

        arg = (outside_water_level - breach_level) / (breach_level + weir_length)
        if (arg <= 0.0) then
            errmsg = 'Chen2018: log argument <= 0'
            error stop trim(errmsg)
        end if

        Cd = c1 * log(arg) + c2
        discharge_coeff = Cd * sqrt(2.0)

   case ('Zerihun2020')
        zeta = (outside_water_level - breach_level) / weir_length
        Cd = (0.40- 0.215 * sin(outer_slope_angle)**(22.0/125.0)+ 0.13  * sin(inner_slope_angle)**( 3.0/ 20.0) + 0.134 * zeta / (1.0 + 0.596 * zeta) )
        discharge_coeff = Cd * sqrt(2.0)

   case default
        errmsg = "Invalid value 'formula'. Allowed: 'Chen2018', 'Zerihun2020'."
        error stop trim(errmsg)
   end select

   end function
    
   function calc_water_depth_width_avg(water_depth, breach_width_avg_water_depth, breach_width_waterline) result(water_depth_width_avg)
   !Calculates the width-averaged water depth in the breach.
   implicit none
   REAL*4:: water_depth, breach_width_avg_water_depth, breach_width_waterline, water_depth_width_avg
    

   if (water_depth < 0) then
        error stop "Water depth is negative!"
   end if

   if (breach_width_avg_water_depth > breach_width_waterline) then
        error stop "Water depth averaged breach width by definition cannot be larger than breach width at waterline"
   end if

   water_depth_width_avg = water_depth * (breach_width_avg_water_depth / breach_width_waterline)    
   end function
    
   function calc_froude_number(flow_velocity, water_depth_width_avg, slope_angle_dike) result (froude_number)
   !Computes the Froude number of the flow.
   implicit none
   REAL*4:: flow_velocity, water_depth_width_avg, slope_angle_dike, froude_number, g
   real*4, parameter :: pi =  3.141592653589793 ! pi in single precision
    
   g = 9.81

   if (flow_velocity < 0) then
        error stop "Flow velocity is negative!"
   end if

   if (water_depth_width_avg <= 0) then
        error stop "Water depth is non-positive!"
   end if

   if (slope_angle_dike > pi/2 .or. slope_angle_dike < 0) then
        error stop "Slope angle is not within range, the range is: 0 <= radians <= pi/2"
   end if

   froude_number = flow_velocity / (g * water_depth_width_avg * cos(slope_angle_dike))**0.5
    
   end function
    
   function calc_adaptation_length_flow(froude_number, water_depth, slope_angle_dike) result (adaptation_length)
   !Computes the length for the flow to achieve normal flow conditions.
   implicit none
   REAL*4:: froude_number, water_depth, slope_angle_dike, adaptation_length
   real*4, parameter :: pi = 3.141592653589793  ! pi in single precision
   if (froude_number < 1.0) then
       error stop "Froude number is less than 1.0!"
   end if

   if (water_depth <= 0.0) then
        error stop "Water depth is non-positive!"
   end if

   if (slope_angle_dike > pi/2 .or. slope_angle_dike < 0) then
        error stop "Slope angle is not within range, the range is: 0 <= radians <= pi/2"
   end if

   adaptation_length = 2.5 * (froude_number**2 - 1) * water_depth / tan(slope_angle_dike)
    
   end function
    
   function calc_adaptation_length_sediment(discharge, breach_width_avg_water_depth, sediment_fall_velocity, slope_angle_dike) result (adaptation_length)
   !Computes the length for the flow to achieve equilibrium sediment transport capacity.

   implicit none
   REAL*4:: discharge, breach_width_avg_water_depth, sediment_fall_velocity, slope_angle_dike, epsilon, adaptation_length
   real*4, parameter :: pi = 3.141592653589793  ! pi in single precision
    
   epsilon=1.0
    
   if (discharge < 0) then
        error stop "Discharge is negative!"
   end if

   if (slope_angle_dike > pi/2 .or. slope_angle_dike < 0) then
        error stop "Slope angle is not within range, the range is: 0 <= radians <= pi/2"
   end if

   adaptation_length = (epsilon * discharge / (breach_width_avg_water_depth * sediment_fall_velocity * cos(slope_angle_dike)))
   
   end function
    
   function flow_along_slope(discharge, crit_water_depth, normal_water_depth, breach_width_bottom, slope_angle_breach, adaptation_length_flow, loc_along_slope) result (res)
       
   !Computes the flow conditions at a location along the slope.
   implicit none
   REAL*4:: discharge, crit_water_depth, normal_water_depth, breach_width_bottom, slope_angle_breach, adaptation_length_flow, loc_along_slope
   REAL*4:: exponential, water_depth, breach_width_avg_water_depth, slope_flow_velocity
   !type :: NormalFlow
    !real*4 :: water_depth
    !real*4 :: slope_flow_velocity
   !end type NormalFlow

   ! Result
   type(NormalFlow) :: res
   
   ! Calculate the water depth on the slope at the location along the slope
   exponential = exp(-5 * loc_along_slope / adaptation_length_flow)
   water_depth = normal_water_depth + (crit_water_depth - normal_water_depth) * exponential

   ! Calculate the water depth averaged breach width
   breach_width_avg_water_depth = breach_width_bottom + (water_depth / tan(slope_angle_breach))

   ! Calculate the flow velocity on the slope at the location along the slope
   slope_flow_velocity = discharge / (breach_width_avg_water_depth * water_depth)
    
   res%water_depth = water_depth
   res%slope_flow_velocity = slope_flow_velocity

   end function
    
   function iter_friction_coefficient(friction_coeff, flow_velocity, hydraulic_radius, delta, d50, d90, kappa) result (friction_coeff_final)
   !Iteration function to computes the friction coefficient.
   implicit none
   REAL*4:: friction_coeff, flow_velocity, hydraulic_radius, delta, d50, d90, kappa, friction_coeff_iter, theta, g,friction_coeff_final
    
   g = 9.81
   ! Define an iteration variable to start the while loop
   friction_coeff_iter = 0.9 * friction_coeff

   do while (abs((friction_coeff - friction_coeff_iter) / friction_coeff) > 0.01)
        ! Set iteration variable to estimate to use for next iteration
        friction_coeff_iter = friction_coeff
        theta = calc_shields_parameter(flow_velocity, delta, d50, friction_coeff)
        friction_coeff = calc_friction_coefficient(theta, hydraulic_radius, d90, kappa)
   end do
   friction_coeff_final = friction_coeff
   end function
   
   function calc_breach_width_total(breach_width_bottom, crest_level, breach_level, slope_angle_breach) result (breach_width_total)
   !Calculates the breach width at the dike crest
   implicit none
   REAL*4:: breach_width_bottom, crest_level, breach_level, slope_angle_breach, breach_width_total
   real*4, parameter :: pi = 3.141592653589793  ! pi in single precision

   if (crest_level <= breach_level) then
        error stop "Crest level is below breach level!"
   end if

   if (breach_width_bottom < 0) then
        error stop "Breach bottom width is negative!"
   end if

   if (slope_angle_breach > pi/2 .or. slope_angle_breach <= 0) then
        error stop "Breach slope angle is not within range, the range is: 0 < radians <= pi/2"
   end if

   breach_width_total = breach_width_bottom + 2 * ((crest_level - breach_level) / tan(slope_angle_breach))
   end function
    
   !!!!! Sediment function
    
   function bagnold_visser(delta, p, d50, phi, ws, flow_velocity, slope_angle, friction_coeff) result (res)
   !Computes the sediment total transport capacity following the Bagnold-Visser (1989) formulation.
   !Formulation limits:
   !11 < theta < 106
   !2.8 < froude < 4.1
   !100 <= d50 <= 220 micron
   !0.36 < tan(slope_angle) < 0.62
   !1.2 < flow velocity < 3.5 m/s
   !0.007 < concentration < 0.28
    
   implicit none
   REAL*4:: delta, p, d50, phi, flow_velocity, slope_angle, friction_coeff, ws
   REAL*4:: sb_max, sb, ss, stt, g, n
   !type :: NormalFlow
   !  real*4 :: n
    ! real*4 :: stt
   !end type NormalFlow

   ! Result
   type(NormalFlow) :: res

   ! velocity exponent
   n = 4.0
   g = 9.81

   ! upper limit bed load transport
   sb_max = 2 * (1 - p) * d50 * flow_velocity

   if (slope_angle >= phi) then
        sb = sb_max
   else
        ! sediment bed load transport
        sb = (0.13 / ((tan(phi) - tan(slope_angle)) * cos(slope_angle)) * friction_coeff * flow_velocity**(n-1) / (delta * g))

        ! used sediment bed load transport
        sb = min(sb, sb_max)
   end if
    
   ! sediment suspended load transport
   ss = 0.01 * friction_coeff * flow_velocity**n / (delta * g * ws * cos(slope_angle)**2)

   ! sediment total transport
   stt = sb + ss
    
   res%n = n
   res%stt = stt
   end function   
   
   
   !!!! Stages
   
   function stage_1(t0, breach_bottom, polder_level, polder_water_level, breach_level, beta1, beta0, outside_water_level, gamma0, alpha, W, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity) result(results_t1)
   implicit none
   REAL*4:: t0, polder_level, polder_water_level, breach_bottom, breach_level, beta1, beta0
   REAL*4:: outside_water_level, gamma0
   REAL*4:: alpha, W, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity
   REAL*4:: beta, crit_water_depth
   REAL*4:: crit_breach_width_avg_water_depth, crit_breach_width_waterline
   REAL*4:: crit_flow_velocity, discharge_coeff, discharge
   REAL*4:: breach_width_avg_water_depth, breach_width_waterline, breach_width_total
   REAL*4:: normal_water_depth, normal_flow_velocity, friction_coeff
   REAL*4:: normal_breach_width_avg_water_depth, normal_breach_width_waterline
   REAL*4:: normal_water_depth_width_avg, froude_number
   REAL*4:: adap_length_flow, adap_length_sediment, length_slope
   REAL*4:: slope_loc
   REAL*4:: water_depth_slope, flow_velocity_slope
   REAL*4:: breach_width_avg_water_depth_slope, hydraulic_radius_slope, friction_coeff_slope
   REAL*4:: theta, stt, n
   REAL*4:: t1
   
   logical :: FLOWSLOPE,sediment_density,water_density
   
   type(NormalFlow) :: results_t1
   
   type(NormalFlow) :: nf
   !real(dp) :: normal_water_depth, normal_flow_velocity, friction_coeff_out
   
   type(NormalFlow) :: res
   !real(dp) :: water_depth, slope_flow_velocity
   
   type(NormalFlow) :: sediment_transport_capacity
   !real(dp) :: n, stt
   
   
   beta = (beta1 + beta0) / 2  ! average inner slope angle (radians)
   

   ! Compute critical flow depth
   crit_water_depth = calc_crit_water_depth(outside_water_level, breach_level, breach_bottom, gamma0)
   write(logstr,'(a,f6.4)') 'crit_water_depth:', crit_water_depth
   call write_log(logstr,1)
   
   ! Compute breach widths for critical flow depth
   crit_breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, crit_water_depth,gamma0)
   crit_breach_width_waterline = calc_breach_width_waterline(breach_bottom, crit_water_depth, gamma0)

   ! Compute critical flow velocity
   crit_flow_velocity = calc_crit_flow_velocity(crit_water_depth, crit_breach_width_avg_water_depth,  crit_breach_width_waterline)
   write(logstr,'(a,f6.4)') 'crit_flow_velocity:', crit_flow_velocity
   call write_log(logstr,1)
   
   discharge_coeff = calc_discharge_coeff('Zerihun2020', alpha, beta, outside_water_level, breach_level, W)
   write(logstr,'(a,f6.4)') 'discharge_coeff:', discharge_coeff
   call write_log(logstr,1)

   ! Compute discharge
   discharge = calc_discharge(crit_breach_width_avg_water_depth, crit_water_depth, crit_flow_velocity, discharge_coeff)
   write(logstr,'(a,f6.5)') 'discharge:', discharge
   call write_log(logstr,1)

   ! Compute breach widths
   breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, crit_water_depth, gamma0)
   breach_width_waterline = calc_breach_width_waterline(breach_bottom, crit_water_depth, gamma0)
   breach_width_total = calc_breach_width_total(breach_bottom, crest_level, breach_level, gamma0)

   nf = normal_flow_conditions(outside_water_level, breach_level, breach_bottom, discharge, gamma0, beta, delta, d50, d90, Cf, kappa)
    
   normal_water_depth   = nf%normal_water_depth
   
   write(logstr,'(a,f6.4)') 'normal_water_depth:', normal_water_depth
   call write_log(logstr,1)
   normal_flow_velocity = nf%normal_flow_velocity
   
   write(logstr,'(a,f6.4)') 'normal_flow_velocity:', normal_flow_velocity
   call write_log(logstr,1)
   friction_coeff   = nf%friction_coeff

   ! Compute breach width for normal flow depth
   normal_breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, normal_water_depth, gamma0)
   normal_breach_width_waterline = calc_breach_width_waterline(breach_bottom, normal_water_depth, gamma0)

   ! Compute the width-averaged normal water depth
   normal_water_depth_width_avg = calc_water_depth_width_avg(normal_water_depth, normal_breach_width_avg_water_depth, normal_breach_width_waterline)

   ! Compute Froude number
   froude_number = calc_froude_number(normal_flow_velocity, normal_water_depth_width_avg, beta)

   ! Compute adaptation length to reach normal flow conditions
   adap_length_flow = calc_adaptation_length_flow(froude_number, normal_water_depth, beta)

   ! Compute adaptation length to reach equilibrium sediment transport
   adap_length_sediment = calc_adaptation_length_sediment(discharge, normal_breach_width_avg_water_depth, sediment_fall_velocity, beta1)
   adap_length_sediment = adap_length_sediment * breach_width_waterline / breach_width_total

   ! Sediment capacity adaptation length cannot be smaller than normal flow
   ! adaptation length (as long as velocity increases, capacity increases)
   adap_length_sediment = max(adap_length_sediment, adap_length_flow)
   ! Average inner slope length from breach level (Zbr) to inner toe (here: Zp)
   length_slope = (breach_level - polder_level) / sin(beta)

   ! If normal flow conditions are not reached at the inner toe, set the
   ! location for flow conditions that determine sediment transport at inner toe.
   slope_loc = min(adap_length_flow, length_slope - (polder_water_level / sin(beta)))

   ! Determine the sediment transport. If FLOWSLOPE == True, determine sediment transport at slope_loc.
   if (FLOWSLOPE) then
        res = flow_along_slope(discharge, crit_water_depth, normal_water_depth, breach_bottom, gamma0, adap_length_flow, slope_loc)
        water_depth_slope = res%water_depth
        flow_velocity_slope = res%slope_flow_velocity
        
        breach_width_avg_water_depth_slope = calc_breach_width_avg_water_depth(breach_bottom, water_depth_slope,gamma0)
        hydraulic_radius_slope = calc_hydraulic_radius(water_depth_slope, breach_bottom, breach_width_avg_water_depth_slope, gamma0)
        friction_coeff_slope = iter_friction_coefficient(friction_coeff, flow_velocity_slope, hydraulic_radius_slope, delta, d50, d90, kappa)
        theta = calc_shields_parameter(flow_velocity_slope, delta, d50, friction_coeff_slope)
        sediment_transport_capacity = bagnold_visser(delta, p, d50, phi, sediment_fall_velocity, flow_velocity_slope, beta1, friction_coeff_slope)
        n= sediment_transport_capacity%n
        stt = sediment_transport_capacity%stt 
        ![stt, n] = sediment_transport_capacity(params, flow_velocity_slope, water_depth_slope, beta, adap_length_sediment, theta, friction_coeff_slope, params.formula_stc1)
        ! Compute time duration of stage I
        t1 = t0 + ((breach_width_total / breach_width_waterline) * (1 - p) * adap_length_sediment * (beta1 - beta0) * slope_loc / (stt * (1 + 5 * n * (slope_loc / adap_length_flow) *(normal_flow_velocity / flow_velocity_slope - 1))))
   else
        theta = calc_shields_parameter(normal_flow_velocity, delta, d50, friction_coeff) 
        sediment_transport_capacity = bagnold_visser(delta, p, d50, phi, sediment_fall_velocity, normal_flow_velocity, beta1, friction_coeff)
        n= sediment_transport_capacity%n
        stt = sediment_transport_capacity%stt
        ![stt, n] = sediment_transport_capacity(params, normal_flow_velocity, normal_water_depth, beta, adap_length_sediment, theta, friction_coeff, params.formula_stc1)
        ! Compute time duration of stage I
        t1 = t0 + ((breach_width_total / breach_width_waterline) * (1 - p) * adap_length_sediment * (beta1 - beta0) * slope_loc / stt)
   end if
   
   results_t1%discharge = discharge
   results_t1%t1  = t1
   end function

end module
