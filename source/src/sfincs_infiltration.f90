module sfincs_infiltration

   use sfincs_log
   use sfincs_error
    
contains

   subroutine initialize_infiltration()
   !
   use sfincs_data
   use sfincs_ncinput   
   !
   implicit none
   !
   integer :: nm
   !
   logical :: ok
   !
   character*256 :: varname
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
   ! Four options for infiltration:
   !
   ! 1) Spatially-uniform constant infiltration
   !    Requires: -
   ! 2) Spatially-varying constant infiltration
   !    Requires: qinfmap (does not require qinffield !)
   ! 3) Spatially-varying infiltration with CN numbers (old)
   !    Requires: cumprcp, cuminf, qinfmap, qinffield
   ! 4) Spatially-varying infiltration with CN numbers (new)
   !    Requires: qinfmap, qinffield, qinffield, ksfield, scs_P1, scs_F1, scs_Se and scs_rain (but not necessarily cuminf and cumprcp)
   ! 5) Spatially-varying infiltration with the Green-Ampt (GA) model
   !    Requires: qinfmap, qinffield, ksfield, GA_head, GA_sigma_max, GA_Lu
   ! 6) Spatially-varying infiltration with the modified Horton Equation 
   !    Requires: qinfmap, qinffield, horton_fc, horton_f0     
   !
   ! cumprcp and cuminf are stored in the netcdf output if store_cumulative_precipitation == .true. which is the default
   !
   ! We need to keep cumprcp and cuminf in memory when:
   !   a) store_cumulative_precipitation == .true.
   ! or:  
   !   b) inftype == 'cna' or inftype == 'cnb'
   !   
   !!!!!!!!!!!!!!!!!!!!!
   ! Initializing steps:
   !!!!!!!!!!!!!!!!!!!!!
   !
   ! 1) First we determine infiltration type
   !
   if (precip) then
      !
      if (infiltrationfile  /= 'none') then
         !
         ! inftype is user defined, keyword: 'inftype' in sfincs.inp:
         !
         ! inftype is either: c2d, cna, cnb, gai, hor
         ! 'inftype = con' is not relevant for netcdf input
         !
         infiltration = .true.
         netcdf_infiltration = .true.
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
         !
      elseif (sefffile /= 'none') then  
         !
         ! Spatially-varying infiltration with CN numbers (new)
         !
         inftype = 'cnb'
         infiltration = .true.      
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
         write(logstr,'(a,a)')'Info    : reading netcdf infiltration file ', trim(infiltrationfile)
         call write_log(logstr, 0)
         !
         ok = check_file_exists(infiltrationfile, 'Infiltration netcdf file', .true.)
         !
         write(logstr,'(a,a)')'Info    : specified inftype is ', trim(inftype)
         call write_log(logstr, 0)
         !
      endif
      !
      ! 5) Check whether infiltration input type (orignal vs netcdf) are correctly matched to grid type (regular vs quadtree)
      !
      if (infiltration .and. inftype /= 'con') then !constant uniform works for both options
         ! 
         if (netcdf_infiltration) then
            !            
            if (use_quadtree .eqv. .false.) then
               ! 
               call stop_sfincs('Error ! Netcdf infiltration input format can only be specified for quadtree mesh model !', 1)              
               !
            endif
            !
         else ! Original
            ! 
            if (use_quadtree .eqv. .true.) then
               ! 
               call stop_sfincs('Error ! Infiltration input for quadtree mesh model can only be specified using the infiltrationfile Netcdf format! !', 1)                       
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
         ! Spatially-uniform constant infiltration (specified as +mm/hr)
         !
         ! Note : Input directly in sfincs.inp, so no file needs to be read
         !
         write(logstr,'(a)')'Info    : turning on spatially-uniform constant infiltration'        
         call write_log(logstr, 0)
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
      elseif (inftype == 'c2d') then
         !
         ! Spatially-varying constant infiltration
         !
         write(logstr,'(a)')'Info    : turning on spatially-varying constant infiltration'      
         call write_log(logstr, 0)
         !
         allocate(qinffield(np))
         !
         qinffield = 0.0         
         !
         ! Read spatially-varying infiltration (specified in +mm/hr)
         !
         if (netcdf_infiltration) then
            !
            ! Call the generic quadtree nc file reader function
            varname = 'qinf'
            call read_netcdf_quadtree_to_sfincs(infiltrationfile, varname, qinffield) !ncfile, varname, varout)               
            !
         else ! from separate qinffile - only binary:
            ! 
            write(logstr,'(a,a)')'Info    : reading infiltration file ', trim(qinffile)
            call write_log(logstr, 0)
            !
            ok = check_file_exists(qinffile, 'Infiltration qinf file', .true.)
            !
            open(unit = 500, file = trim(qinffile), form = 'unformatted', access = 'stream')
            read(500)qinffield
            close(500)
            !
         endif
         !
         ! Generic needed conversion:
         !
         qinffield = qinffield / 3600 / 1000   ! convert to +m/s         
         !
      elseif (inftype == 'cna') then
         !
         ! Spatially-varying infiltration with CN numbers (old)
         !
         write(logstr,'(a)')'Info    : turning on infiltration (via Curve Number method - A)'               
         call write_log(logstr, 0)
         !
         allocate(qinffield(np))
         !
         qinffield = 0.0
         !
         if (netcdf_infiltration) then
            !
            ! Call the generic quadtree nc file reader function
            varname = 'scs'
            call read_netcdf_quadtree_to_sfincs(infiltrationfile, varname, qinffield)               
            !
         else ! from separate scsfile - only binary:      
            !
            write(logstr,'(a,a)')'Info    : reading scs file ',trim(scsfile)
            call write_log(logstr, 0)
            !
            ok = check_file_exists(scsfile, 'Infiltration scs file', .true.)
            !
            open(unit = 500, file = trim(scsfile), form = 'unformatted', access = 'stream')
            read(500)qinffield
            close(500)
            !
         endif
         !
         ! Generic needed conversion:         
         !
         qinffield = qinffield * 0.0254   ! to m
         ! already convert qinffield from inches to m here         
         !
      elseif (inftype == 'cnb') then  
         !
         ! Spatially-varying infiltration with CN numbers (new)
         !
         write(logstr,'(a)')'Info    : turning on infiltration (via Curve Number method - B)' 
         call write_log(logstr, 0)
         ! 
         ! Allocate Smax
         allocate(qinffield(np))
         qinffield = 0.0
         !
         if (netcdf_infiltration) then
            !
            ! Call the generic quadtree nc file reader function
            varname = 'smax'
            call read_netcdf_quadtree_to_sfincs(infiltrationfile, varname, qinffield)              
            !
         else ! from separate smaxfile - only binary:
            ! 
            write(logstr,'(a,a)')'Info    : reading smax file ',trim(smaxfile)
            call write_log(logstr, 0)
            !
            ok = check_file_exists(smaxfile, 'Infiltration smax file', .true.)
            !
            open(unit = 500, file = trim(smaxfile), form = 'unformatted', access = 'stream')
            read(500)qinffield
            close(500)
            ! 
         endif         
         !
         ! Allocate Se
         allocate(scs_Se(np))
         scs_Se = 0.0
         !
         if (netcdf_infiltration) then
            !
            ! Call the generic quadtree nc file reader function
            varname = 'seff'
            call read_netcdf_quadtree_to_sfincs(infiltrationfile, varname, scs_Se)              
            !
         else ! from separate sefffile - only binary:         
            !
            write(logstr,'(a,a)')'Info    : reading seff file ',trim(sefffile)
            call write_log(logstr, 0)
            !
            ok = check_file_exists(sefffile, 'Infiltration seff file', .true.)
            !
            open(unit = 501, file = trim(sefffile), form = 'unformatted', access = 'stream')
            read(501)scs_Se
            close(501)
            !
         endif
         !
         ! Allocate Ks
         !
         allocate(ksfield(np))
         ksfield = 0.0
         !
         if (netcdf_infiltration) then
            !
            ! Call the generic quadtree nc file reader function
            varname = 'ks'
            call read_netcdf_quadtree_to_sfincs(infiltrationfile, varname, ksfield)              
            !
         else ! from separate ksfile - only binary:         
            !
            write(logstr,'(a,a)')'Info    : reading ks file ',trim(ksfile)
            call write_log(logstr, 0)
            !
            ok = check_file_exists(ksfile, 'Infiltration ks file', .true.)
            !
            open(unit = 502, file = trim(ksfile), form = 'unformatted', access = 'stream')
            read(502)ksfield
            close(502)
            !
         endif
         !
         ! Generic needed conversion:         
         !         
         ! Compute recovery                     ! Equation 4-36
         !
         allocate(inf_kr(np))
         inf_kr = sqrt(ksfield/25.4) / 75       ! Note that we assume ksfield to be in mm/hr, convert it here to inch/hr (/25.4)
                                                ! /75 is conversion to recovery rate (in days)
         !
         ! Allocate support variables:
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
      elseif (inftype == 'gai') then  
         !
         ! Spatially-varying infiltration with the Green-Ampt (GA) model
         !
         call write_log('Info    : turning on process infiltration (via Green-Ampt)', 0)
         ! 
         ! Allocate suction head at the wetting front 
         !
         allocate(GA_head(np))
         GA_head = 0.0
         !
         if (netcdf_infiltration) then
            !
            ! Call the generic quadtree nc file reader function
            varname = 'psi'
            call read_netcdf_quadtree_to_sfincs(infiltrationfile, varname, GA_head)             
            !
         else ! from separate psifile - only binary:
             !
             write(logstr,'(a,a)')'Info    : reading psi file ',trim(psifile)
             call write_log(logstr, 0)
             !
             ok = check_file_exists(psifile, 'Infiltration psi file', .true.)
             !
             open(unit = 500, file = trim(psifile), form = 'unformatted', access = 'stream')
             read(500)GA_head
             close(500)
             !
         endif
         !
         ! Allocate maximum soil moisture deficit
         !
         allocate(GA_sigma_max(np))
         GA_sigma_max = 0.0
         !
         if (netcdf_infiltration) then
            !
            ! Call the generic quadtree nc file reader function
            varname = 'sigma'
            call read_netcdf_quadtree_to_sfincs(infiltrationfile, varname, GA_sigma_max)
            !
         else ! from separate sigmafile - only binary:
             !         
             write(logstr,'(a,a)')'Info    : reading sigma file ',trim(sigmafile)
             call write_log(logstr, 0)
             !
             ok = check_file_exists(sigmafile, 'Infiltration sigma file', .true.)
             !
             open(unit = 501, file = trim(sigmafile), form = 'unformatted', access = 'stream')
             read(501)GA_sigma_max
             close(501)
             !
         endif
         !
         ! Allocate saturated hydraulic conductivity
         !
         allocate(ksfield(np))
         ksfield = 0.0
         !
         if (netcdf_infiltration) then
            !
            ! Call the generic quadtree nc file reader function
            varname = 'ks'
            call read_netcdf_quadtree_to_sfincs(infiltrationfile, varname, ksfield)
            !
         else ! from separate ksfile - only binary:
             !             
             write(logstr,'(a,a)')'Info    : reading ks file ',trim(ksfile)
             call write_log(logstr, 0)
             !
             ok = check_file_exists(ksfile, 'Infiltration ks file', .true.)
             !
             open(unit = 502, file = trim(ksfile), form = 'unformatted', access = 'stream')
             read(502)ksfield
             close(502)
             !
         endif
         
         !
         ! Generic needed conversion:         
         !          
         ! Compute recovery                         ! Equation 4-36
         !
         allocate(inf_kr(np))
         inf_kr     = sqrt(ksfield/25.4) / 75       ! Note that we assume ksfield to be in mm/hr, convert it here to inch/hr (/25.4)
                                                    ! /75 is conversion to recovery rate (in days)
         
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
         ! Input values for green-ampt are in mm and mm/hr, but computation is in m a m/s
         !
         GA_head    = GA_head / 1000                  ! from mm to m
         GA_Lu      = GA_Lu / 1000                    ! from mm to m
         ksfield    = ksfield / 1000 / 3600           ! from mm/hr to m/s
         ! 
         ! First time step doesnt have an estimate yet
         !
         ! Allocate support variables:
         !
         allocate(qinffield(np))
         qinffield(nm) = 0.0
         !
      elseif (inftype == 'hor') then  
         !
         ! Spatially-varying infiltration with the modified Horton Equation 
         !
         call write_log('Info    : turning on process infiltration (via modified Horton)', 0)
         ! 
         ! Horton: final infiltration capacity (fc)
         ! Note that qinffield = horton_fc (/3600/1000, see below)
         !
         allocate(horton_fc(np))
         horton_fc = 0.0
         !
         if (netcdf_infiltration) then
            !
            ! Call the generic quadtree nc file reader function
            varname = 'fc'
            call read_netcdf_quadtree_to_sfincs(infiltrationfile, varname, horton_fc)
            !
         else ! from separate fcfile - only binary:
            !
            write(logstr,'(a,a)')'Info    : reading fc file ',trim(fcfile)
            call write_log(logstr, 0)
            !
            ok = check_file_exists(fcfile, 'Infiltration fc file', .true.)
            !
            open(unit = 500, file = trim(fcfile), form = 'unformatted', access = 'stream')
            read(500)horton_fc
            close(500)
            !
         endif
         !
         ! Horton: initial infiltration capacity (f0)
         allocate(horton_f0(np))
         horton_f0 = 0.0
         !
         if (netcdf_infiltration) then
            !
            ! Call the generic quadtree nc file reader function
            varname = 'f0'
            call read_netcdf_quadtree_to_sfincs(infiltrationfile, varname, horton_f0)
            !
         else ! from separate f0file - only binary:
            !         
            write(logstr,'(a,a)')'Info    : reading f0 file ',trim(f0file)
            call write_log(logstr, 0)
            !
            ok = check_file_exists(f0file, 'Infiltration f0 file', .true.)
            !
            open(unit = 501, file = trim(f0file), form = 'unformatted', access = 'stream')
            read(501)horton_f0
            close(501)
            !
         endif         
         !
         ! Empirical constant (1/hr) k => note that this is different than ks used in Curve Number and Green-Ampt
         allocate(horton_kd(np))
         horton_kd = 0.0
         !
         if (netcdf_infiltration) then
            !
            ! Call the generic quadtree nc file reader function
            varname = 'kd'
            call read_netcdf_quadtree_to_sfincs(infiltrationfile, varname, horton_kd)
            !
         else ! from separate kdfile - only binary:
            !             
            write(logstr,'(a,a)')'Info    : reading kd file ',trim(kdfile)
            call write_log(logstr, 0)
            !
            ok = check_file_exists(kdfile, 'Infiltration kd file', .true.)
            !
            open(unit = 502, file = trim(kdfile), form = 'unformatted', access = 'stream')
            read(502)horton_kd
            close(502)
            !
         endif
         !
         write(logstr,'(a,a)')'Info    : Using constant recovery rate that is based on constant factor relative to ',trim(kdfile)
         call write_log(logstr, 0)         
         !
         ! Generic needed conversion:         
         !         
         ! Prescribe the current estimate (for output only; initial capacity)
         qinffield = horton_f0/3600/1000
         !
         ! Allocate support variables:
         !
         ! Estimate of time
         allocate(rain_T1(np))
         rain_T1 = 0.0
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
   
   
   subroutine update_infiltration_map(dt, tloop)
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
   real*4  :: hh_local, a
   real*4  :: dt   
   !
   integer   :: count0
   integer   :: count1
   integer   :: count_rate
   integer   :: count_max
   real      :: tloop
   !
   call system_clock(count0, count_rate, count_max)
   !
   if (inftype == 'con' .or. inftype == 'c2d') then
      !
      ! Infiltration rate map stays constant
      !
      !$omp parallel &
      !$omp private ( nm )
      !$omp do
      !$acc parallel present( qinfmap, qinffield, z_volume, zs, zb, netprcp, cuminf )
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
         ! Compute nett precip
         !
         netprcp(nm) = netprcp(nm) - qinfmap(nm)
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
   elseif (inftype == 'cna') then
      !
      ! Determine infiltration rate with Curve Number (old method; no recovery)
      !
      !$omp parallel &
      !$omp private ( Qq,I,nm )
      !$omp do
      !$acc parallel present( qinfmap, qinffield, prcp, netprcp, cumprcp, cuminf )
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
         ! Compute nett precip
         !
         netprcp(nm) = netprcp(nm) - qinfmap(nm)
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
   elseif (inftype == 'cnb') then
      !
      ! Determine infiltration rate with Curve Number with recovery
      !
      !$omp parallel &
      !$omp private ( Qq,I,nm )       
      !$omp do       
      !$acc parallel present( qinfmap, prcp, netprcp, cuminf, scs_rain, scs_Se, scs_P1, scs_F1, scs_S1, rain_T1, qinffield, inf_kr )
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
         ! Compute nett precip
         !
         netprcp(nm) = netprcp(nm) - qinfmap(nm)
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
   elseif (inftype == 'gai') then
      !
      ! Determine infiltration rate with  with the Green-Ampt (GA) model
      !
      !$omp parallel &
      !$omp private ( nm )
      !$omp do              
      !$acc parallel present( qinfmap, prcp, netprcp, cuminf, rain_T1,  &
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
               qinfmap(nm) = (ksfield(nm) * (1.0 + (GA_head(np) * GA_sigma(np)) / GA_F(nm)))
               qinfmap(nm) = max(min(qinfmap(nm), prcp(nm)), 0.0)     ! never more than rainfall and and never negative
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
         ! Compute nett precip
         !
         !qinffield(nm)  = qinfmap(nm) ! Really ? Why ?
         netprcp(nm)    = netprcp(nm) - qinfmap(nm)
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
   elseif (inftype == 'hor') then
      !
      ! Determine infiltration rate with with the Horton model
      !
      !$omp parallel &
      !$omp private  ( nm, Qq, I, a, hh_local )
      !$omp do              
      !$acc parallel present( qinfmap, prcp, netprcp, cuminf, cell_area_m2, cell_area, z_flags_iref, z_volume, zs, zb, rain_T1,  &
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
            ! Count how long this is already going
            !
            if (rain_T1(nm) > 0.0) then
               !
               rain_T1(nm) = 0.0
               !
            endif
            !
            rain_T1(nm) = rain_T1(nm) - dt                                              ! negative amount of how long it is infiltrating
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
         ! Compute nett precip
         !
         netprcp(nm)    = netprcp(nm) - qinfmap(nm)
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
   endif
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0 * (count1 - count0) / count_rate
   !
   end subroutine   

end module
