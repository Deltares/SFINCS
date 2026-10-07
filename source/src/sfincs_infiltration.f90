module sfincs_infiltration

   use sfincs_log
   use sfincs_error
    
contains

   subroutine initialize_infiltration()
   !
   use sfincs_data
   !
   implicit none
   !
   logical :: ok
   !
   character(len=3), parameter :: allowed_types(7) = &
        ['c2d', 'cna', 'cnb', 'gai', 'hor', 'bkt', 'gwt']

   logical :: inftype_exists   
   !
   ! INFILTRATION
   !
   ! Infiltration only works when rainfall is activated ! If you want infiltration without rainfall, use a precip file with 0.0s
   !
   ! Note, infiltration methods not designed to be stacked
   !
   infiltration   = .false.
   netcdf_infiltration   = .false.   
   !
   ! Seven infiltration flavors (inftype):
   !
   ! 1) 'con' - Spatially-uniform constant infiltration
   !    Requires: qinf (mm/hr in sfincs.inp)
   ! 2) 'c2d' - Spatially-varying constant infiltration
   !    Requires: qinffile or inffile
   ! 3) 'cna' - SCS Curve Number (old, no recovery)
   !    Requires: scsfile or inffile
   ! 4) 'cnb' - SCS Curve Number (new, with recovery)
   !    Requires: sefffile or inffile
   ! 5) 'gai' - Green-Ampt infiltration
   !    Requires: psifile or inffile
   ! 6) 'hor' - Modified Horton equation
   !    Requires: f0file or inffile
   ! 7) 'bkt' - Bucket model (linear reservoir, HBV/wflow style)
   !    Requires: inffile with bucket_smax, bucket_k and bucket_loss
   !
   ! cumprcp and cuminf are stored in the netcdf output if store_cumulative_precipitation == .true. (storecumprcp = 1)
   !
   ! We need to keep cumprcp and cuminf updated when:
   !   a) store_cumulative_precipitation == .true.
   ! or:
   !   b) inftype == 'cna' or inftype == 'cnb' (store_cumulative_precipitation is then forced to .true.)
   !   
   !!!!!!!!!!!!!!!!!!!!!
   ! Initializing steps:
   !!!!!!!!!!!!!!!!!!!!!
   !
   ! 1) First we determine infiltration type
   !
   if (precip) then
      !
      if (inftype == 'bkt' .and. inffile == 'none') then
         !
         call stop_sfincs('Error ! Bucket model requires inffile together with inftype = bkt !', 1)
         !
      endif
      !
      if (inffile  /= 'none') then
         !
         ! inftype is user defined, keyword: 'inftype' in sfincs.inp:
         !
         ! inftype is either: c2d, cna, cnb, gai, hor
         ! 'inftype = con' is not relevant for netcdf input
         !
         ! Check if specified type is correct
         !
         inftype_exists = any(inftype == allowed_types)
         !
         if (inftype_exists) then
            !
            infiltration = .true.
            netcdf_infiltration = .true.  
            !
            write(logstr,'(a,a)')'Info    : specified inftype is ', trim(inftype)
            call write_log(logstr, 0)
            !
            ! Curve Number methods need cumprcp and cuminf to be updated (same as binary input path)
            !
            if (inftype == 'cna' .or. inftype == 'cnb') then
               store_cumulative_precipitation = .true.
            endif
            !
         else
            !
            write(logstr,*)'Error    : infiltration input type ',trim(inftype),' is not part of supported types c2d cna cnb gai hor bkt !'
            call stop_sfincs(trim(logstr), 1)   
            !
         end if
         !
         !
      elseif (inftype == 'gwt') then
         !
         ! Groundwater table model with uniform parameters from sfincs.inp (no inffile)
         !
         infiltration = .true.
         !
      elseif (qinf > 0.0) then
         !
         ! Spatially-uniform constant infiltration (specified as +mm/hr)
         !
         inftype = 'con'
         infiltration = .true.
         ! 
      elseif (qinffile /= 'none') then
         !
         ! Spatially-varying constant infiltration
         !
         inftype = 'c2d'
         infiltration = .true.      
         !
      elseif (scsfile /= 'none') then
         !
         ! Spatially-varying infiltration with CN numbers (old)
         !
         inftype = 'cna'
         infiltration = .true.     
         store_cumulative_precipitation = .true.         
         !
      elseif (sefffile /= 'none') then  
         !
         ! Spatially-varying infiltration with CN numbers (new)
         !
         inftype = 'cnb'
         infiltration = .true.      
         store_cumulative_precipitation = .true.         
         !
      elseif (psifile /= 'none') then  
         !
         ! The Green-Ampt (GA) model for infiltration
         !
         inftype = 'gai'
         infiltration = .true.      
         !
      elseif (f0file /= 'none') then
         !
         ! The Horton Equation model for infiltration
         !
         inftype        = 'hor'
         infiltration   = .true.
         store_meteo    = .true.
         !
      endif
      !
      ! 2) We need cumprcp and cuminf
      !
      allocate(cumprcp(np))
      cumprcp = 0.0
      !
      allocate(cuminf(np))
      cuminf = 0.0
      !
      ! 3) Now allocate and read spatially-varying inputs 
      !
      if (infiltration) then
         !
         allocate(qinfmap(np))
         qinfmap = 0.0
         ! 
      endif
      !
      ! 4) Pre-check whether netcdf infiltration file exists - once
      !
      if (netcdf_infiltration) then
         !      
         write(logstr,'(a)')'Info    : turning on infiltration from netcdf input file'      
         call write_log(logstr, 0)
         !
         write(logstr,'(a,a)')'Info    : reading netcdf infiltration file ', trim(inffile)
         call write_log(logstr, 0)
         !
         ok = check_file_exists(inffile, 'Infiltration netcdf file', .true.)
         !
      endif
      !
      ! 5) Check whether infiltration input type (orignal vs netcdf) are correctly matched to grid type (regular vs quadtree)
      !
      if (infiltration .and. inftype /= 'con' .and. inftype /= 'gwt') then !constant uniform and gwt with uniform parameters work for both options
         !
         ! Netcdf infiltration works for both regular and quadtree grids
         ! (regular grids populate quadtree_nr_points and index_sfincs_in_quadtree
         !  via make_quadtree_from_indices)
         !
         if (.not. netcdf_infiltration) then
            !
            if (use_quadtree .eqv. .true.) then
               !
               call stop_sfincs('Error ! Infiltration input for quadtree mesh model can only be specified using the inffile Netcdf format! !', 1)
               !
            endif
            !
         endif
         !
      endif      
      !      
      ! 6) Read in data per type, either from ascii or general netcdf file
      !
      if (inftype == 'con') then
         !
         call initialize_infiltration_con()
         !
      elseif (inftype == 'c2d') then
         !
         call initialize_infiltration_c2d()
         !
      elseif (inftype == 'cna') then
         !
         call initialize_infiltration_cna()
         !
      elseif (inftype == 'cnb') then
         !
         call initialize_infiltration_cnb()
         !
      elseif (inftype == 'gai') then
         !
         call initialize_infiltration_gai()
         !
      elseif (inftype == 'hor') then
         !
         call initialize_infiltration_hor()
         !
      elseif (inftype == 'bkt') then
         !
         ! Bucket model (linear reservoir) - mimics hydrology models like wflow/HBV
         !
         call write_log('Info    : turning on process infiltration (via bucket model)', 0)
         !
         call initialize_bucket_model()
         !
      elseif (inftype == 'gwt') then
         !
         ! 0D groundwater table model (PRIMo-style)
         !
         call initialize_groundwater_table()
         !
      endif
      !
   else
      !
      ! Overrule input
      !
      store_cumulative_precipitation = .false.
      !
   endif
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
   real*4  :: dt
   !
   if (inftype == 'con' .or. inftype == 'c2d') then
      !
      call compute_infiltration_constant(dt)
      !
   elseif (inftype == 'cna') then
      !
      call compute_infiltration_cna(dt)
      !
   elseif (inftype == 'cnb') then
      !
      call compute_infiltration_cnb(dt)
      !
   elseif (inftype == 'gai') then
      !
      call compute_infiltration_gai(dt)
      !
   elseif (inftype == 'hor') then
      !
      call compute_infiltration_hor(dt)
      !
   elseif (inftype == 'bkt') then
      !
      call compute_bucket_drainage(dt)
      !
   elseif (inftype == 'gwt') then
      !
      call compute_infiltration_gwt(dt)
      !
   endif
   !
   ! Apply the resulting infiltration-rate field to the point-source field
   ! qsrc (m3/s). qinfmap is m/s, so multiply by cell area and subtract.
   ! qsrc already holds this step's prcp*area contribution (from
   ! update_meteo_forcing) plus any discharges / src-structures updates
   ! done earlier in update_continuity.
   !
   !$acc parallel loop present( qsrc, qinfmap, cell_area, cell_area_m2, z_flags_iref )
   !$omp parallel do default(shared) private(nm) schedule(static)
   do nm = 1, np
      !
      if (crsgeo) then
         qsrc(nm) = qsrc(nm) - qinfmap(nm) * cell_area_m2(nm)
      else
         qsrc(nm) = qsrc(nm) - qinfmap(nm) * cell_area(z_flags_iref(nm))
      endif
      !
   enddo
   !$omp end parallel do
   !
   end subroutine


   subroutine initialize_infiltration_con()
   !
   ! Spatially-uniform constant infiltration (specified as +mm/hr)
   !
   ! Note : Input directly in sfincs.inp, so no file needs to be read
   !
   use sfincs_data
   !
   implicit none
   !
   integer :: nm
   !
   call write_log('Info    : turning on spatially-uniform constant infiltration', 0)
   !
   allocate(qinffield(np))
   !
   ! Note : qinf has already been converted to m/s in sfincs_input.f90 !
   !
   do nm = 1, np
      if (subgrid) then
         if (subgrid_z_zmin(nm) > qinf_zmin) then
            qinffield(nm) = qinf
         else
            qinffield(nm) = 0.0
         endif
      else
         if (zb(nm) > qinf_zmin) then
            qinffield(nm) = qinf
         else
            qinffield(nm) = 0.0
         endif
      endif
   enddo
   !
   end subroutine


   subroutine initialize_infiltration_c2d()
   !
   ! Spatially-varying constant infiltration (specified as +mm/hr)
   !
   use sfincs_data
   !
   implicit none
   !
   call write_log('Info    : turning on spatially-varying constant infiltration', 0)
   !
   allocate(qinffield(np))
   qinffield = 0.0
   !
   call read_infiltration_field('qinf', qinffile, qinffield)
   !
   qinffield = qinffield / 3600 / 1000   ! convert to +m/s
   !
   end subroutine


   subroutine initialize_infiltration_cna()
   !
   ! Spatially-varying infiltration with CN numbers (old, no recovery)
   !
   use sfincs_data
   !
   implicit none
   !
   call write_log('Info    : turning on infiltration (via Curve Number method - A)', 0)
   !
   ! qinffield is S (input in inches)
   !
   allocate(qinffield(np))
   qinffield = 0.0
   !
   call read_infiltration_field('scs', scsfile, qinffield)
   !
   qinffield = qinffield * 0.0254   ! inches to m
   !
   end subroutine


   subroutine initialize_infiltration_cnb()
   !
   ! Spatially-varying infiltration with CN numbers (new, with recovery)
   !
   use sfincs_data
   !
   implicit none
   !
   call write_log('Info    : turning on infiltration (via Curve Number method - B)', 0)
   !
   ! Smax (stored in qinffield), Se and Ks
   !
   allocate(qinffield(np))
   qinffield = 0.0
   call read_infiltration_field('smax', smaxfile, qinffield)
   !
   allocate(scs_Se(np))
   scs_Se = 0.0
   call read_infiltration_field('seff', sefffile, scs_Se)
   !
   allocate(ksfield(np))
   ksfield = 0.0
   call read_infiltration_field('ks', ksfile, ksfield)
   !
   ! Compute recovery                     ! Equation 4-36
   !
   allocate(inf_kr(np))
   inf_kr = sqrt(ksfield/25.4) / 75       ! Note that we assume ksfield to be in mm/hr, convert it here to inch/hr (/25.4)
                                          ! /75 is conversion to recovery rate (in days)
   !
   ! Allocate support variables
   !
   allocate(scs_P1(np))
   scs_P1 = 0.0
   allocate(scs_F1(np))
   scs_F1 = 0.0
   allocate(rain_T1(np))
   rain_T1 = 0.0
   allocate(scs_S1(np))
   scs_S1 = 0.0
   allocate(scs_rain(np))
   scs_rain = 0
   !
   end subroutine


   subroutine initialize_infiltration_gai()
   !
   ! Spatially-varying infiltration with the Green-Ampt (GA) model
   !
   use sfincs_data
   !
   implicit none
   !
   call write_log('Info    : turning on process infiltration (via Green-Ampt)', 0)
   !
   ! Suction head at the wetting front (psi), maximum soil moisture deficit (sigma)
   ! and saturated hydraulic conductivity (ks)
   !
   allocate(GA_head(np))
   GA_head = 0.0
   call read_infiltration_field('psi', psifile, GA_head)
   !
   allocate(GA_sigma_max(np))
   GA_sigma_max = 0.0
   call read_infiltration_field('sigma', sigmafile, GA_sigma_max)
   !
   allocate(ksfield(np))
   ksfield = 0.0
   call read_infiltration_field('ks', ksfile, ksfield)
   !
   ! Compute recovery                         ! Equation 4-36
   !
   allocate(inf_kr(np))
   inf_kr     = sqrt(ksfield/25.4) / 75       ! Note that we assume ksfield to be in mm/hr, convert it here to inch/hr (/25.4)
                                              ! /75 is conversion to recovery rate (in days)
   !
   allocate(rain_T1(np))                      ! minimum amount of time that a soil must remain in recovery
   rain_T1    = 0.0
   !
   ! Allocate support variables
   !
   allocate(GA_sigma(np))                     ! variable for sigma_max_du
   GA_sigma   = GA_sigma_max
   allocate(GA_F(np))                         ! total infiltration
   GA_F       = 0.0
   allocate(GA_Lu(np))                        ! depth of upper soil recovery zone
   GA_Lu      = 4 * sqrt(25.4) * sqrt(ksfield) ! Equation 4-33
   !
   ! Input values for green-ampt are in mm and mm/hr, but computation is in m and m/s
   !
   GA_head    = GA_head / 1000                ! from mm to m
   GA_Lu      = GA_Lu / 1000                  ! from mm to m
   ksfield    = ksfield / 1000 / 3600         ! from mm/hr to m/s
   !
   ! First time step doesnt have an estimate yet
   !
   allocate(qinffield(np))
   qinffield = 0.0
   !
   end subroutine


   subroutine initialize_infiltration_hor()
   !
   ! Spatially-varying infiltration with the modified Horton Equation
   !
   use sfincs_data
   !
   implicit none
   !
   call write_log('Info    : turning on process infiltration (via modified Horton)', 0)
   !
   ! Final infiltration capacity (fc), initial infiltration capacity (f0) and
   ! empirical constant kd (1/hr) => note that kd is different than ks used in Curve Number and Green-Ampt
   !
   allocate(horton_fc(np))
   horton_fc = 0.0
   call read_infiltration_field('fc', fcfile, horton_fc)
   !
   allocate(horton_f0(np))
   horton_f0 = 0.0
   call read_infiltration_field('f0', f0file, horton_f0)
   !
   allocate(horton_kd(np))
   horton_kd = 0.0
   call read_infiltration_field('kd', kdfile, horton_kd)
   !
   write(logstr,'(a,a)')'Info    : Using constant recovery rate that is based on constant factor relative to ',trim(kdfile)
   call write_log(logstr, 0)
   !
   ! Prescribe the current estimate (for output only; initial capacity)
   !
   allocate(qinffield(np))
   qinffield = horton_f0/3600/1000
   !
   ! Estimate of time
   !
   allocate(rain_T1(np))
   rain_T1 = 0.0
   !
   end subroutine


   subroutine read_infiltration_field(varname, binfile, field)
   !
   ! Read one spatially-varying infiltration parameter, either from the netcdf
   ! inffile (variable varname) or from a legacy binary file (regular grids only)
   !
   use sfincs_data
   use sfincs_ncinput
   !
   implicit none
   !
   character(len=*), intent(in)         :: varname
   character(len=*), intent(in)         :: binfile
   real*4, dimension(np), intent(inout) :: field
   !
   character*256 :: ncvarname
   logical       :: ok
   !
   if (netcdf_infiltration) then
      !
      ncvarname = varname
      call read_netcdf_quadtree_to_sfincs(inffile, ncvarname, field)
      !
   else
      !
      write(logstr,'(a,a,a,a)')'Info    : reading ', trim(varname), ' file ', trim(binfile)
      call write_log(logstr, 0)
      !
      ok = check_file_exists(binfile, 'Infiltration '//trim(varname)//' file', .true.)
      !
      open(unit = 500, file = trim(binfile), form = 'unformatted', access = 'stream')
      read(500)field
      close(500)
      !
   endif
   !
   end subroutine


   subroutine compute_infiltration_constant(dt)
   !
   ! Constant infiltration (con and c2d): infiltration rate map stays constant
   !
   use sfincs_data
   !
   implicit none
   !
   real*4  :: dt
   !
   integer :: nm
   !
   !$omp parallel &
   !$omp private ( nm )
   !$omp do
   !$acc parallel present( qinfmap, qinffield, z_volume, zs, zb, cuminf )
   !$acc loop independent gang vector
   do nm = 1, np
      !
      qinfmap(nm) = qinffield(nm) ! Set spatially varying infiltration field
      !
      ! No infiltration if there is no water
      !  
      if (subgrid) then
         !
         if (z_volume(nm) <= 0.0) then
            qinfmap(nm) = 0.0
         endif
         !
      else
         !
         if (zs(nm) <= zb(nm)) then
            qinfmap(nm) = 0.0
         endif
         !
      endif
      !
      if (store_cumulative_precipitation) then
         !
         ! Compute cumulative infiltration
         !
         cuminf(nm) = cuminf(nm) + qinfmap(nm) * dt
         !
      endif
      !
   enddo
   !$omp end do
   !$omp end parallel
   !$acc end parallel
   !
   end subroutine


   subroutine compute_infiltration_cna(dt)
   !
   ! Infiltration rate with Curve Number (old method; no recovery)
   !
   use sfincs_data
   !
   implicit none
   !
   real*4  :: dt
   !
   integer :: nm
   real*4  :: Qq
   real*4  :: I
   !
   !$omp parallel &
   !$omp private ( Qq,I,nm )
   !$omp do
   !$acc parallel present( qinfmap, qinffield, prcp, cumprcp, cuminf )
   !$acc loop independent gang vector
   do nm = 1, np
      !
      ! Check if Ia (0.2 x S) is larger than cumulative rainfall
      !
      if (cumprcp(nm) > sfacinf * qinffield(nm)) then ! qinffield is S
         ! 
         ! Compute runoff as function of rain
         !
         Qq  = (cumprcp(nm) - sfacinf * qinffield(nm))**2 / (cumprcp(nm) + (1.0 - sfacinf) * qinffield(nm))  ! cumulative runoff in m
         I   = cumprcp(nm) - Qq                        ! cumulative infiltration in m
         qinfmap(nm) = (I - cuminf(nm)) / dt           ! infiltration in m/s
         !
      else
         !
         ! Everything still infiltrating
         !
         qinfmap(nm) = prcp(nm)
         !
      endif   
      !
      if (store_cumulative_precipitation) then
         !
         ! Compute cumulative infiltration
         !
         cuminf(nm) = cuminf(nm) + qinfmap(nm) * dt
         !
      endif
      !
   enddo
   !$omp end do
   !$omp end parallel
   !$acc end parallel
   !
   end subroutine


   subroutine compute_infiltration_cnb(dt)
   !
   ! Infiltration rate with Curve Number with recovery
   !
   use sfincs_data
   !
   implicit none
   !
   real*4  :: dt
   !
   integer :: nm
   real*4  :: Qq
   real*4  :: I
   !
   !$omp parallel &
   !$omp private ( Qq,I,nm )       
   !$omp do       
   !$acc parallel present( qinfmap, prcp, cuminf, scs_rain, scs_Se, scs_P1, scs_F1, scs_S1, rain_T1, qinffield, inf_kr )
   !$acc loop independent gang vector
   do nm = 1, np
      !
      ! If there is precip in this grid cell for this time step  
      !
      if (prcp(nm) > 0.0) then
         !
         ! Is raining now
         !
         if (scs_rain(nm) == 1) then
            !
            ! It was raining before; do nothing
            !
         else
            !
            ! Initalise these variables for new rainfall event
            !
            scs_P1(nm)          = 0.0               ! cumulative rainfall for this 'event'
            scs_F1(nm)          = 0.0               ! cumulative infiltration for this 'event'
            scs_S1(nm)          = scs_Se(nm)        ! S for this 'event'
            scs_rain(nm)        = 1                 ! logic used to determine if there is an event ongoing
            !
         endif
         ! 
         !  Compute cum rainfall
         ! 
         scs_P1(nm) = scs_P1(nm) + prcp(nm) * dt
         ! 
         ! Compute runoff
         ! 
         if (scs_P1(nm) > (sfacinf * scs_S1(nm)) ) then ! scs_S1 is S
            !
            Qq          = (scs_P1(nm) - (sfacinf * scs_S1(nm)))**2 / (scs_P1(nm) + (1.0 - sfacinf) * scs_S1(nm))  ! cumulative runoff in m
            I           = scs_P1(nm) - Qq                       ! cum infiltration this event
            qinfmap(nm) = (I - scs_F1(nm))/dt                   ! infiltration in m/s
            scs_F1(nm)  = I                                     ! cum infiltration this event
            !
         else
            !
            Qq          = 0.0                                   ! no runoff
            scs_F1(nm)  = scs_P1(nm)                            ! all rainfall is infiltrated
            qinfmap(nm) = prcp(nm)                              ! infiltration rate = rainfall rate
            !
         endif
         ! 
         ! Compute "remaining S", but note that scs_Se is not used in computation
         ! 
         scs_Se(nm)  = max(scs_Se(nm) - qinfmap(nm) * dt, 0.0)
         qinfmap(nm) = max(qinfmap(nm), 0.0)
         !
      else
         ! 
         ! It is not raining here
         !
         if (scs_rain(nm) == 1) then
            !
            ! if it was raining before; cange logic and set rate to 0
            !
            scs_rain(nm)   = 0
            qinfmap(nm)    = 0.0
            rain_T1(nm)    = 0.0
            !
         endif
         !
         ! Add to recovery time
         !
         rain_T1(nm) = rain_T1(nm) + dt / 3600
         !
         ! compute recovery of S if time is larger than this
         !
         if (rain_T1(nm) > (0.06 / inf_kr(nm)) ) then	! Equation 4-37 from SWMM
            !
            ! note that scs_Se is S and qinffield is Smax
            !
            scs_Se(nm)  = scs_Se(nm) + (inf_kr(nm) * qinffield(nm) * dt / 3600)  ! scs_kr is recovery in hours 
            scs_Se(nm)  = min(scs_Se(nm), qinffield(nm))
            !
         endif
         !
      endif
      !
      if (store_cumulative_precipitation) then
         !
         ! Compute cumulative infiltration
         !
         cuminf(nm) = cuminf(nm) + qinfmap(nm)*dt
         !
      endif
      !
   enddo
   !$omp end do
   !$omp end parallel
   !$acc end parallel
   !
   end subroutine


   subroutine compute_infiltration_gai(dt)
   !
   ! Infiltration rate with the Green-Ampt (GA) model
   !
   use sfincs_data
   !
   implicit none
   !
   real*4  :: dt
   !
   integer :: nm
   !
   !$omp parallel &
   !$omp private ( nm )
   !$omp do              
   !$acc parallel present( qinfmap, prcp, cuminf, rain_T1,  &
   !$acc                  ksfield, GA_head, GA_sigma, GA_sigma_max, GA_F, GA_Lu, inf_kr )
   !$acc loop independent gang vector
   do nm = 1, np
      !
      ! If there is precip in this grid cell for this time step?
      !
      if (prcp(nm) > 0.0) then
         !
         ! Is raining now
         !
         if (prcp(nm) < ksfield(nm)) then
            !
            ! Small amounts of rainfall - infiltration is same as soil
            !
            qinfmap(nm) = prcp(nm)                       ! infiltration is same as rainfall
            !
         else
            !
            ! Larger amounts of rainfall - Equation 4-27 from SWMM manual
            !
            if (GA_F(nm) < 1.0e-10) then
               !
               ! No cumulative infiltration yet (first timestep) - all rainfall infiltrates
               !
               qinfmap(nm) = prcp(nm)
               !
            else
               !
               qinfmap(nm) = (ksfield(nm) * (1.0 + (GA_head(nm) * GA_sigma(nm)) / GA_F(nm)))
               qinfmap(nm) = max(min(qinfmap(nm), prcp(nm)), 0.0)     ! never more than rainfall and never negative
               !
            endif
            !
         endif
         !
         ! Update sigma 
         !
         GA_sigma(nm) = max(GA_sigma(nm) - (qinfmap(nm) * dt / GA_Lu(nm)), 0.0)
         ! 
         ! Update others
         !
         GA_F(nm)    = GA_F(nm) + qinfmap(nm) * dt   ! internal cumulative rainfall from Green-Ampt
         rain_T1(nm) = 0.0                           ! recovery time not started
         !
      else
         ! 
         ! Not raining here
         !
         ! Add to recovery time
         !
         rain_T1(nm)     = rain_T1(nm) + dt / 3600      
         !
         ! Compute recovery of S if time is larger than this
         !
         if (rain_T1(nm) > (0.06 / inf_kr(nm)) ) then			! Equation 4-37 from SWMM
            ! 
            ! Update sigma 
            !
            GA_sigma(nm) = GA_sigma(nm) + (inf_kr(nm) * GA_sigma_max(nm) * dt / 3600)       ! Equation 4-35
            GA_sigma(nm) = min(GA_sigma(nm), GA_sigma_max(nm))                              ! never more than max
            !
            ! Update internal cumulative rainfall
            !
            GA_F(nm)    = max(GA_F(nm) - (inf_kr(nm) * GA_sigma_max(nm) * dt / 3600 * GA_Lu(nm)), 0.0)    ! Page 112 SWMM
            !
         endif
      endif
      ! 
      if (store_cumulative_precipitation) then
         !
         ! Compute cumulative infiltration
         !
         cuminf(nm)  = cuminf(nm) + qinfmap(nm) * dt
         !
      endif
      !
   enddo
   !$omp end do
   !$omp end parallel
   !$acc end parallel
   !
   end subroutine


   subroutine compute_infiltration_hor(dt)
   !
   ! Infiltration rate with the modified Horton model
   !
   use sfincs_data
   !
   implicit none
   !
   real*4  :: dt
   !
   integer :: nm
   real*4  :: Qq
   real*4  :: I
   real*4  :: hh_local
   !
   !$omp parallel &
   !$omp private  ( nm, Qq, I, hh_local )
   !$omp do              
   !$acc parallel present( qinfmap, prcp, cuminf, cell_area_m2, cell_area, z_flags_iref, z_volume, zs, zb, rain_T1,  &
   !$acc                  horton_kd, horton_fc, horton_f0 )
   !$acc loop independent gang vector
   do nm = 1, np
      !
      ! Get local water depth estimate
      !
      if (subgrid) then
         !
         if (crsgeo) then
            !
            hh_local = z_volume(nm) / cell_area_m2(nm)
            !
         else   
            !
            hh_local = z_volume(nm) / cell_area(z_flags_iref(nm))
            !
         endif
         !
      else
         !
         hh_local = zs(nm) - zb(nm)
         !
      endif
      !
      ! Check if there is water
      !
      if (hh_local> 0.0 .or. prcp(nm) > 0.0) then
         !
         ! Infiltrating here
         !
         ! Count how long this is already going.
         ! If rain_T1 was positive (recovery phase), reset it to 0 for this storm onset
         ! and do NOT apply the decrement yet — otherwise the first time step of a new
         ! storm would start with rain_T1 = -dt, underestimating infiltration capacity.
         !
         if (rain_T1(nm) > 0.0) then
            !
            rain_T1(nm) = 0.0
            !
         else
            !
            rain_T1(nm) = rain_T1(nm) - dt                                           ! negative amount of how long it is infiltrating
            !
         endif
         ! 
         ! Compute estimate of infiltration                                          ! Note that qinffield = horton_fc
         !
         I = exp(horton_kd(nm) * rain_T1(nm) / 3600)                                 ! note that horton_kd is factor in hours while dt is seconds
         !
         ! Stop keeping track of this when less than 1% left (same for time)
         !
         if (I < 0.01) then
            !
            I           = 0.0                           ! which reduces qinfmap to horton_fc
            rain_T1(nm) = rain_T1(nm) + dt              ! also make sure time doesnt further decrease
            qinfmap(nm) = horton_fc(nm) / 3600 / 1000   ! from mm/hr to m/s
            !
         else
            !
            qinfmap(nm) = (horton_fc(nm) + (horton_f0(nm) - horton_fc(nm)) * I) / 3600 / 1000 ! from mm/hr to m/s
            !
         endif
         !
         ! Check how much there can infiltrate
         !
         if (hh_local > 0.0) then
            !
            ! Qq = prcp(nm) * dt + (zs(nm) - zb(nm))  ! Qq is estimate in meter of how much water there is
            Qq = prcp(nm) * dt + hh_local             ! Qq is estimate in meter of how much water there is (MvO: using hh_local instead?)
            !
         else
            !
            Qq = prcp(nm) * dt                        ! if no water; only compare with rainfall
            !
         endif
         !
         ! Compare how much Horton wants to infiltrate
         !
         I = qinfmap(nm) * dt                                              ! I is estimate in meter of how much Horton allows
         !
         if (I > Qq) then
            !
            qinfmap(nm) = qinfmap(nm) * Qq / I                             ! scale Horton if capacity > available
            !
         endif
         !
      else
         !
         ! Not raining here NOR ponding
         !
         rain_T1(nm) = rain_T1(nm) + dt / horton_kr_kd                 ! positive amount of how long it is infiltrating
         qinfmap(nm) = 0.0
         !
      endif
      !
      if (store_cumulative_precipitation) then
         !
         ! Compute cumulative infiltration
         !
         cuminf(nm)  = cuminf(nm) + qinfmap(nm) * dt
         !
      endif
      !
   enddo
   !$omp end do
   !$omp end parallel
   !$acc end parallel
   !
   end subroutine


   subroutine initialize_bucket_model()
   !
   use sfincs_data
   use sfincs_ncinput
   !
   implicit none
   !
   character*256 :: varname
   !
   if (netcdf_infiltration) then
      !
      write(logstr,'(a)')'Info    : turning on bucket model (linear reservoir)'
      call write_log(logstr, 0)
      !
      allocate(bucket_capacity(np))
      allocate(bucket_k(np))
      allocate(bucket_volume(np))
      allocate(bucket_drain_rate(np))
      allocate(bucket_loss(np))
      allocate(bucket_runoff(np))
      !
      bucket_capacity   = 0.0
      bucket_k          = 0.0
      bucket_volume     = 0.0
      bucket_drain_rate = 0.0
      bucket_loss       = 0.0
      bucket_runoff     = 0.0
      !
      ! Read from inffile (netcdf) - works for both regular and quadtree grids
      ! (read_netcdf_quadtree_to_sfincs stops if a variable is missing)
      !
      varname = 'bucket_smax'
      call read_netcdf_quadtree_to_sfincs(inffile, varname, bucket_capacity)
      bucket_capacity = bucket_capacity / 1000.0   ! mm to m
      !
      varname = 'bucket_k'
      call read_netcdf_quadtree_to_sfincs(inffile, varname, bucket_k)
      bucket_k = bucket_k / 3600.0   ! 1/hr to 1/s
      !
      varname = 'bucket_loss'
      call read_netcdf_quadtree_to_sfincs(inffile, varname, bucket_loss)
      !
      write(logstr,'(a,f10.4,a)')'Info    : bucket max capacity = ', maxval(bucket_capacity) * 1000.0, ' mm'
      call write_log(logstr, 0)
      write(logstr,'(a,f10.4,a)')'Info    : bucket max k        = ', maxval(bucket_k) * 3600.0, ' 1/hr'
      call write_log(logstr, 0)
      write(logstr,'(a,f6.3)')'Info    : bucket loss fraction = ', maxval(bucket_loss)
      call write_log(logstr, 0)
      !
   else
      !
      ! Allocate minimal arrays for OpenACC compatibility
      !
      allocate(bucket_capacity(1))
      allocate(bucket_k(1))
      allocate(bucket_volume(1))
      allocate(bucket_drain_rate(1))
      allocate(bucket_loss(1))
      allocate(bucket_runoff(1))
      bucket_capacity   = 0.0
      bucket_k          = 0.0
      bucket_volume     = 0.0
      bucket_drain_rate = 0.0
      bucket_loss       = 0.0
      bucket_runoff     = 0.0
      !
   endif
   !
   end subroutine


   subroutine compute_bucket_drainage(dt)
   !
   ! Bucket model with loss: linear reservoir + loss fraction (HBV/wflow style)
   !
   ! Steps per cell:
   !   1. P_eff = P * (1 - loss)       -- fraction lost to ET/deep percolation
   !   2. Fill bucket with P_eff (up to Smax capacity)
   !   3. Drain bucket: S(t+dt) = S(t)*exp(-k*dt), drainage returned as runoff
   !   4. qinfmap = P - runoff         -- net removal from surface
   !
   ! In continuity: zs += prcp*dt - qinfmap*dt = bucket_runoff*dt
   ! => Only bucket drainage reaches the surface water level
   !
   ! Literature: Linear reservoir (Nash, 1957), HBV soil moisture bucket (Bergstrom, 1995)
   !
   use sfincs_data
   !
   implicit none
   !
   real*4           :: dt
   integer          :: nm
   real*4           :: exp_factor
   real*4           :: drain_vol
   real*4           :: P_eff
   real*4           :: available_cap
   real*4           :: actual_inflow
   real*4           :: precip_rate
   !
   !$omp parallel do private(nm, exp_factor, drain_vol, P_eff, available_cap, actual_inflow, precip_rate)
   !$acc parallel present( kcs, prcp, qinfmap, cuminf, bucket_volume, bucket_capacity, bucket_k, &
   !$acc                   bucket_drain_rate, bucket_loss, bucket_runoff )
   !$acc loop independent gang vector
   do nm = 1, np
      !
      if (kcs(nm) == 1 .and. bucket_k(nm) > 0.0) then
         !
         ! Step 1: Compute effective precipitation (after loss)
         !
         precip_rate = max(prcp(nm), 0.0)
         P_eff = precip_rate * (1.0 - bucket_loss(nm))             ! m/s after loss
         !
         ! Step 2: Fill bucket with effective precip (up to capacity)
         !
         if (bucket_capacity(nm) > 0.0) then
            available_cap = bucket_capacity(nm) - bucket_volume(nm)
            actual_inflow = min(P_eff * dt, available_cap)          ! m
         else
            ! No capacity limit (Smax = 0 means infinite)
            actual_inflow = P_eff * dt                              ! m
         endif
         bucket_volume(nm) = bucket_volume(nm) + actual_inflow
         !
         ! Step 3: Drain bucket (analytical linear reservoir)
         ! S(t+dt) = S(t) * exp(-k*dt), drainage = S(t) - S(t+dt)
         !
         exp_factor = exp(-bucket_k(nm) * dt)
         drain_vol = bucket_volume(nm) * (1.0 - exp_factor)        ! m drained this step
         bucket_volume(nm) = bucket_volume(nm) * exp_factor
         !
         ! Step 4: Bucket drainage becomes runoff returned to surface
         !
         bucket_runoff(nm) = drain_vol / dt                         ! m/s
         !
         ! Step 5: Set qinfmap = loss + what entered bucket - what drained back
         ! In continuity: zs += prcp*dt - qinfmap*dt
         ! Water balance: qinfmap = prcp*loss + actual_inflow/dt - bucket_runoff
         ! When bucket has room:  actual_inflow = P_eff*dt => qinfmap = prcp - bucket_runoff
         ! When bucket is full:   actual_inflow = 0       => qinfmap can be negative (drainage > inflow)
         !
         qinfmap(nm) = precip_rate * bucket_loss(nm) + actual_inflow / dt - bucket_runoff(nm)
         !
         bucket_drain_rate(nm) = bucket_runoff(nm)
         !
         if (store_cumulative_precipitation) then
            cuminf(nm) = cuminf(nm) + qinfmap(nm) * dt
         endif
         !
      else
         !
         qinfmap(nm) = 0.0
         bucket_drain_rate(nm) = 0.0
         bucket_runoff(nm) = 0.0
         !
      endif
      !
   enddo
   !$acc end parallel
   !$omp end parallel do
   !
   end subroutine


   subroutine initialize_groundwater_table()
   !
   ! 0D groundwater table model (Sanders et al. 2025, PRIMo, eqs. 5-8).
   !
   ! Per cell:  Sy dH/dt = phi * f - w,   w = kappa * max(H - H0, 0),   kappa = Keff / l0
   !
   ! Allocates the state and parameter arrays, reads optional fields from
   ! inffile, applies the uniform overrides from sfincs.inp, disables the
   ! aquifer in open-water and inactive cells, and seeds the water table
   ! (from gw_depth0, or from the restart file when present).
   !
   use sfincs_data
   use sfincs_ncinput
   use quadtree
   !
   implicit none
   !
   integer       :: nm
   integer       :: ipq
   integer       :: n_no_receiver
   character*256 :: varname
   logical       :: have_depth, have_fmax, have_phi, have_sy, have_keff, have_l0
   !
   real*4, dimension(:), allocatable :: gw_depth0
   real*4, dimension(:), allocatable :: gw_keff
   real*4, dimension(:), allocatable :: gw_l0
   integer*4, dimension(:), allocatable :: gw_receiver_qt
   !
   call write_log('Info    : turning on groundwater table model (0D, PRIMo-style)', 0)
   !
   allocate(gw_rise(np))
   allocate(gw_level(np))
   allocate(gw_level0(np))
   allocate(gw_zground(np))
   allocate(gw_fmax(np))
   allocate(gw_phi(np))
   allocate(gw_sy(np))
   allocate(gw_kappa(np))
   allocate(gw_seepage(np))
   allocate(gw_cumseep(np))
   allocate(gw_receiver(np))
   allocate(gw_depth0(np))
   allocate(gw_keff(np))
   allocate(gw_l0(np))
   !
   gw_rise = 0.0
   gw_level = 0.0
   gw_level0 = 0.0
   gw_zground = 0.0
   gw_fmax = 0.0
   gw_phi = 0.0
   gw_sy = 0.0
   gw_kappa = 0.0
   gw_seepage = 0.0
   gw_cumseep = 0.0
   gw_receiver = 0
   gw_depth0 = 0.0
   gw_keff = 0.0
   gw_l0 = 0.0
   !
   ! Ground level: zb on regular grids, lowest subgrid pixel on subgrid grids
   !
   do nm = 1, np
      !
      if (subgrid) then
         gw_zground(nm) = subgrid_z_zmin(nm)
      else
         gw_zground(nm) = zb(nm)
      endif
      !
   enddo
   !
   ! 1) Optional spatially-varying fields from inffile
   !
   have_depth = .false.
   have_fmax = .false.
   have_phi = .false.
   have_sy = .false.
   have_keff = .false.
   have_l0 = .false.
   !
   if (netcdf_infiltration) then
      !
      varname = 'gw_depth0'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_depth0)
         have_depth = .true.
      endif
      !
      varname = 'gw_fmax'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_fmax)
         gw_fmax = gw_fmax / 3.6e6   ! mm/hr to m/s
         have_fmax = .true.
      endif
      !
      varname = 'gw_phi'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_phi)
         have_phi = .true.
      endif
      !
      varname = 'gw_sy'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_sy)
         have_sy = .true.
      endif
      !
      varname = 'gw_keff'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_keff)
         gw_keff = gw_keff / 3.6e6   ! mm/hr to m/s
         have_keff = .true.
      endif
      !
      varname = 'gw_l0'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_l0)
         have_l0 = .true.
      endif
      !
   endif
   !
   ! 2) Uniform overrides from sfincs.inp (already converted to SI in sfincs_input)
   !
   if (gw_depth_ini >= 0.0) then
      gw_depth0 = gw_depth_ini
      have_depth = .true.
   endif
   !
   if (gw_fmax_uniform >= 0.0) then
      gw_fmax = gw_fmax_uniform
      have_fmax = .true.
   endif
   !
   if (gw_phi_uniform >= 0.0) then
      gw_phi = gw_phi_uniform
      have_phi = .true.
   endif
   !
   if (gw_sy_uniform >= 0.0) then
      gw_sy = gw_sy_uniform
      have_sy = .true.
   endif
   !
   if (gw_keff_uniform >= 0.0) then
      gw_keff = gw_keff_uniform
      have_keff = .true.
   endif
   !
   if (gw_l0_uniform >= 0.0) then
      gw_l0 = gw_l0_uniform
      have_l0 = .true.
   endif
   !
   ! 3) Defaults and checks
   !
   if (.not. have_depth) then
      gw_depth0 = 1.0
      call write_log('Warning : gw_depth0 not specified, using initial depth to groundwater of 1.0 m', 0)
   endif
   !
   if (.not. have_phi) then
      gw_phi = 1.0
   endif
   !
   if (.not. have_sy) then
      gw_sy = 0.3
   endif
   !
   if (.not. have_fmax) then
      call stop_sfincs('Error ! inftype = gwt requires gw_fmax in sfincs.inp or a gw_fmax field in inffile !', 1)
   endif
   !
   if (have_keff .and. have_l0) then
      !
      do nm = 1, np
         if (gw_l0(nm) > 0.0) then
            gw_kappa(nm) = gw_keff(nm) / gw_l0(nm)
         else
            gw_kappa(nm) = 0.0
         endif
      enddo
      !
   else
      !
      gw_kappa = 0.0
      call write_log('Warning : gw_keff and/or gw_l0 not specified, groundwater seepage is switched off', 0)
      !
   endif
   !
   ! 4) Seepage destination
   !
   select case (trim(gw_seepage_mode))
   case ('loss')
      gw_seepage_imode = 0
   case ('local')
      gw_seepage_imode = 1
   case ('receiver')
      gw_seepage_imode = 2
   case default
      call stop_sfincs('Error ! gw_seepage_mode should be receiver, loss or local !', 1)
   end select
   !
   if (gw_seepage_imode == 2) then
      !
      ! Receiver cells: gw_receiver in inffile holds the 1-based quadtree index of the
      ! surface-water cell that receives the seepage of each cell. Map to sfincs indices;
      ! cells without a valid receiver lose their seepage.
      !
      varname = 'gw_receiver'
      !
      if (.not. netcdf_infiltration) then
         call stop_sfincs('Error ! gw_seepage_mode = receiver requires inffile with a gw_receiver field !', 1)
      endif
      !
      if (.not. netcdf_quadtree_variable_exists(inffile, varname)) then
         call stop_sfincs('Error ! gw_seepage_mode = receiver requires a gw_receiver field in inffile !', 1)
      endif
      !
      allocate(gw_receiver_qt(np))
      gw_receiver_qt = 0
      !
      call read_netcdf_quadtree_to_sfincs_int(inffile, varname, gw_receiver_qt)
      !
      n_no_receiver = 0
      !
      do nm = 1, np
         !
         ipq = gw_receiver_qt(nm)
         !
         if (ipq >= 1 .and. ipq <= quadtree_nr_points) then
            gw_receiver(nm) = index_sfincs_in_quadtree(ipq)
         else
            gw_receiver(nm) = 0
         endif
         !
         if (gw_receiver(nm) <= 0 .and. gw_zground(nm) > qinf_zmin .and. kcs(nm) == 1) then
            n_no_receiver = n_no_receiver + 1
         endif
         !
      enddo
      !
      deallocate(gw_receiver_qt)
      !
      if (n_no_receiver > 0) then
         write(logstr,'(a,i0,a)')'Warning : ', n_no_receiver, ' groundwater cells have no valid receiver cell, their seepage is lost'
         call write_log(logstr, 0)
      endif
      !
   endif
   !
   ! 5) No aquifer in open-water (ground below qinf_zmin) or inactive cells.
   !    phi = 0 switches the cell off in compute_infiltration_gwt.
   !
   do nm = 1, np
      !
      if (gw_zground(nm) <= qinf_zmin .or. kcs(nm) /= 1 .or. gw_sy(nm) <= 0.0) then
         gw_phi(nm) = 0.0
         gw_kappa(nm) = 0.0
      endif
      !
      gw_sy(nm) = max(gw_sy(nm), 1.0e-3)
      !
   enddo
   !
   ! 6) Seed the water table. The seepage baseline always follows gw_depth0;
   !    a restart only replaces the current level.
   !
   do nm = 1, np
      gw_level0(nm) = gw_zground(nm) - max(gw_depth0(nm), 0.0)
      gw_level(nm) = gw_level0(nm)
   enddo
   !
   if (allocated(gw_rise_rst)) then
      !
      gw_rise = gw_rise_rst
      gw_level = gw_level0 + gw_rise
      deallocate(gw_rise_rst)
      !
   endif
   !
   write(logstr,'(a,f8.3,a,f8.3,a)')'Info    : gw initial depth    = ', minval(gw_depth0), ' - ', maxval(gw_depth0), ' m'
   call write_log(logstr, 0)
   write(logstr,'(a,f8.3,a,f8.3,a)')'Info    : gw fmax             = ', minval(gw_fmax) * 3.6e6, ' - ', maxval(gw_fmax) * 3.6e6, ' mm/hr'
   call write_log(logstr, 0)
   write(logstr,'(a,f8.3,a,f8.3)')'Info    : gw pervious fraction= ', minval(gw_phi), ' - ', maxval(gw_phi)
   call write_log(logstr, 0)
   write(logstr,'(a,f8.3,a,f8.3)')'Info    : gw specific yield   = ', minval(gw_sy), ' - ', maxval(gw_sy)
   call write_log(logstr, 0)
   write(logstr,'(a,e10.3,a,e10.3,a)')'Info    : gw kappa            = ', minval(gw_kappa), ' - ', maxval(gw_kappa), ' 1/s'
   call write_log(logstr, 0)
   write(logstr,'(a,a)')'Info    : gw seepage mode     = ', trim(gw_seepage_mode)
   call write_log(logstr, 0)
   !
   deallocate(gw_depth0)
   deallocate(gw_keff)
   deallocate(gw_l0)
   !
   end subroutine initialize_groundwater_table


   subroutine compute_infiltration_gwt(dt)
   !
   ! 0D groundwater table model, one time step.
   !
   !   depth to groundwater d = z_ground - H
   !   f = fmax             if d > 0 and water on the surface      (inundation case)
   !   f = min(prcp, fmax)  if d > 0 and no water on the surface   (rainfall case)
   !   f = 0                if d <= 0                              (saturation case)
   !   f is capped by the water available this step and by the aquifer space.
   !
   !   Surface sink  qinf_loc = phi * f   (mass-consistent with the aquifer gain)
   !   Aquifer       Sy dH/dt = qinf_loc - kappa * (H - H0), solved exactly over dt.
   !                 The state is the rise H - H0 (gw_rise) so that micrometre
   !                 changes are resolved in single precision; gw_level is derived.
   !                 with qinf_loc held constant, giving the seepage w from the balance.
   !
   !   Seepage goes to: nowhere (loss), the same cell (local), or a receiver
   !   cell via qsrc (receiver). qinfmap is picked up by update_infiltration_map.
   !
   use sfincs_data
   !
   implicit none
   !
   real*4  :: dt
   !
   integer :: nm
   integer :: nmr
   real*4  :: hh_local
   real*4  :: pr
   real*4  :: depth
   real*4  :: f
   real*4  :: qinf_loc
   real*4  :: rise
   real*4  :: rise_new
   real*4  :: efac
   real*4  :: w
   real*4  :: area_loc
   !
   !$omp parallel do private(nm, nmr, hh_local, pr, depth, f, qinf_loc, rise, rise_new, efac, w, area_loc) schedule(static)
   !$acc parallel present( kcs, prcp, zs, zb, z_volume, cell_area, cell_area_m2, z_flags_iref, subgrid_z_zmin, qinfmap, cuminf, qsrc, &
   !$acc                   gw_rise, gw_level, gw_level0, gw_zground, gw_fmax, gw_phi, gw_sy, gw_kappa, gw_seepage, gw_cumseep, gw_receiver )
   !$acc loop independent gang vector
   do nm = 1, np
      !
      qinf_loc = 0.0
      w = 0.0
      !
      if (gw_phi(nm) > 0.0) then
         !
         ! Water available on the surface (m)
         !
         if (subgrid) then
            !
            if (crsgeo) then
               hh_local = z_volume(nm) / cell_area_m2(nm)
            else
               hh_local = z_volume(nm) / cell_area(z_flags_iref(nm))
            endif
            !
         else
            !
            hh_local = zs(nm) - zb(nm)
            !
         endif
         !
         hh_local = max(hh_local, 0.0)
         pr = max(prcp(nm), 0.0)
         !
         depth = (gw_zground(nm) - gw_level0(nm)) - gw_rise(nm)
         !
         if (depth > 0.0) then
            !
            if (hh_local > 0.0) then
               f = gw_fmax(nm)
            else
               f = min(pr, gw_fmax(nm))
            endif
            !
            ! Never remove more than is available on the surface this step
            !
            f = min(f, (hh_local + pr * dt) / dt)
            !
            ! Never overfill the aquifer in one step
            !
            f = min(f, depth * gw_sy(nm) / (gw_phi(nm) * dt))
            f = max(f, 0.0)
            !
            qinf_loc = gw_phi(nm) * f
            !
         endif
         !
         ! Aquifer update: exact solution of Sy dH/dt = qinf_loc - kappa (H - H0) over dt
         !
         rise = gw_rise(nm)
         !
         if (gw_kappa(nm) > 0.0) then
            !
            efac = exp(-gw_kappa(nm) * dt / gw_sy(nm))
            rise_new = rise * efac + qinf_loc / gw_kappa(nm) * (1.0 - efac)
            w = max(qinf_loc - (rise_new - rise) * gw_sy(nm) / dt, 0.0)
            !
         else
            !
            rise_new = rise + qinf_loc * dt / gw_sy(nm)
            !
         endif
         !
         gw_rise(nm) = rise_new
         gw_level(nm) = gw_level0(nm) + rise_new
         !
      endif
      !
      gw_seepage(nm) = w
      gw_cumseep(nm) = gw_cumseep(nm) + w * dt
      !
      if (gw_seepage_imode == 1) then
         qinfmap(nm) = qinf_loc - w
      else
         qinfmap(nm) = qinf_loc
      endif
      !
      if (store_cumulative_precipitation) then
         cuminf(nm) = cuminf(nm) + qinfmap(nm) * dt
      endif
      !
      ! Route seepage to the receiving cell
      !
      if (gw_seepage_imode == 2 .and. w > 0.0) then
         !
         nmr = gw_receiver(nm)
         !
         if (nmr > 0) then
            !
            if (crsgeo) then
               area_loc = cell_area_m2(nm)
            else
               area_loc = cell_area(z_flags_iref(nm))
            endif
            !
            !$acc atomic update
            !$omp atomic
            qsrc(nmr) = qsrc(nmr) + w * area_loc
            !
         endif
         !
      endif
      !
   enddo
   !$acc end parallel
   !$omp end parallel do
   !
   end subroutine compute_infiltration_gwt

end module
