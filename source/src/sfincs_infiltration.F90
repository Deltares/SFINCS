#define NF90(nf90call) call handle_err(nf90call,__FILE__,__LINE__)
module sfincs_infiltration
   !
   use netcdf    
   !
   type net_type_inf
       integer :: ncid
       integer :: np_dimid
       integer :: mask_varid, type_varid
       integer :: qinffield_varid, scs_varid, scs_se_varid, ksfield_varid
   end type      
   type(net_type_inf) :: net_file_inf      
   !
contains
    
    subroutine read_infiltration_file_original()
    !
    ! Six options for infiltration:
    !
    ! 1) Spatially-uniform constant infiltration
    !    Requires: -
    ! 2) Spatially-varying constant infiltration
    !    Requires: qinfmap (does not require qinffield !)
    ! 3) Spatially-varying infiltration with CN numbers (old)
    !    Requires: cumprcp, cuminf, qinfmap, qinffield
    ! 4) Spatially-varying infiltration with CN numbers (new)
    !    Requires: qinfmap, qinffield, ksfield, scs_P1, scs_F1, scs_Se and scs_rain (but not necessarily cuminf and cumprcp)
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
    ! First we determine infiltration type
    !
    use sfincs_data
    !
    implicit none
    !
    integer :: nm      
    !
    if (qinf > 0.0) then   
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
    if (precip) then
        !
        ! We need cumprcp and cuminf
        !
        allocate(cumprcp(np))
        cumprcp = 0.0
        !
        allocate(cuminf(np))
        cuminf = 0.0
        ! 
    endif
    !
    ! Now allocate and read spatially-varying inputs 
    !
    if (infiltration) then
        !
        allocate(qinfmap(np))
        qinfmap = 0.0
        ! 
    endif
    !
    if (inftype == 'con') then   
        !
        ! Spatially-uniform constant infiltration (specified as +mm/hr)
        !
        write(*,*)'Turning on process: Spatially-uniform constant infiltration'        
        !
        allocate(qinffield(np))
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
        write(*,*)'Turning on process: Spatially-varying constant infiltration'      
        !
        ! Read spatially-varying infiltration (only binary, specified in +mm/hr)
        !
        write(*,*)'Reading ', trim(qinffile), ' ...'
        allocate(qinffield(np))
        open(unit = 500, file = trim(qinffile), form = 'unformatted', access = 'stream')
        read(500)qinffield
        close(500)
        !
        do nm = 1, np
        qinffield(nm) = qinffield(nm)/3.6e3/1.0e3   ! convert to +m/s
        enddo
        !
    elseif (inftype == 'cna') then
        !
        ! Spatially-varying infiltration with CN numbers (old)
        !
        write(*,*)'Turning on process: Infiltration (via Curve Number method - A)'               
        !
        allocate(qinffield(np))
        qinffield = 0.0
        !
        write(*,*)'Reading ',trim(scsfile)
        open(unit = 500, file = trim(scsfile), form = 'unformatted', access = 'stream')
        read(500)qinffield
        close(500)
        !
        ! already convert qinffield from inches to m here
        do nm = 1, np
        qinffield(nm) = qinffield(nm)*0.0254   !to m
        enddo      
        !
    elseif (inftype == 'cnb') then  
        !
        ! Spatially-varying infiltration with CN numbers (new)
        !
        write(*,*)'Turning on process: Infiltration (via Curve Number method - B)'            
        ! 
        ! Allocate Smax
        allocate(qinffield(np))
        qinffield = 0.0
        write(*,*)'Reading ',trim(smaxfile)
        open(unit = 500, file = trim(smaxfile), form = 'unformatted', access = 'stream')
        read(500)qinffield
        close(500)
        !
        ! Allocate Se
        allocate(scs_Se(np))
        scs_Se = 0.0
        write(*,*)'Reading ',trim(sefffile)
        open(unit = 501, file = trim(sefffile), form = 'unformatted', access = 'stream')
        read(501)scs_Se
        close(501)
        !
        ! Compute recovery                     ! Equation 4-36        
        ! Allocate Ks
        allocate(ksfield(np))
        ksfield = 0.0
        write(*,*)'Reading ',trim(ksfile)
        open(unit = 502, file = trim(ksfile), form = 'unformatted', access = 'stream')
        read(502)ksfield
        close(502)
        !
        ! Compute recovery                     ! Equation 4-36
        allocate(inf_kr(np))
        inf_kr = sqrt(ksfield/25.4) / 75       ! Note that we assume ksfield to be in mm/hr, convert it here to inch/hr (/25.4)
                                            ! /75 is conversion to recovery rate (in days)
        !
        ! Allocate support variables
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
        write(*,*)'Turning on process: Infiltration (via Green-Ampt)'            
        ! 
        ! Allocate suction head at the wetting front 
        allocate(GA_head(np))
        GA_head = 0.0
        write(*,*)'Reading ',trim(psifile)
        open(unit = 500, file = trim(psifile), form = 'unformatted', access = 'stream')
        read(500)GA_head
        close(500)
        !
        ! Allocate maximum soil moisture deficit
        allocate(GA_sigma_max(np))
        GA_sigma_max = 0.0
        write(*,*)'Reading ',trim(sigmafile)
        open(unit = 501, file = trim(sigmafile), form = 'unformatted', access = 'stream')
        read(501)GA_sigma_max
        close(501)
        !
        ! Allocate saturated hydraulic conductivity
        allocate(ksfield(np))
        ksfield = 0.0
        write(*,*)'Reading ',trim(ksfile)
        open(unit = 502, file = trim(ksfile), form = 'unformatted', access = 'stream')
        read(502)ksfield
        close(502)
        !
        ! Compute recovery                         ! Equation 4-36
        allocate(inf_kr(np))
        inf_kr     = sqrt(ksfield/25.4) / 75       ! Note that we assume ksfield to be in mm/hr, convert it here to inch/hr (/25.4)
                                                ! /75 is conversion to recovery rate (in days)
         
        allocate(rain_T1(np))                      ! minimum amount of time that a soil must remain in recovery 
        rain_T1    = 0.0         
        ! Allocate support variables
        allocate(GA_sigma(np))                     ! variable for sigma_max_du
        GA_sigma   = GA_sigma_max
        allocate(GA_F(np))                         ! total infiltration
        GA_F       = 0.0
        allocate(GA_Lu(np))                        ! depth of upper soil recovery zone
        GA_Lu      = 4 *sqrt(25.4) * sqrt(ksfield) ! Equation 4-33
        !
        ! Input values for green-ampt are in mm and mm/hr, but computation is in m a m/s
        GA_head    = GA_head/1000                  ! from mm to m
        GA_Lu      = GA_Lu/1000                    ! from mm to m
        ksfield    = ksfield/1000/3600             ! from mm/hr to m/s
        ! 
        ! First time step doesnt have an estimate yet
        allocate(qinffield(np))
        qinffield(nm) = 0.00
        !
    elseif (inftype == 'hor') then  
        !
        ! Spatially-varying infiltration with the modified Horton Equation 
        !
        write(*,*)'Turning on process: Infiltration (via modified Horton)'            
        ! 
        ! Horton: final infiltration capacity (fc)
        ! Note that qinffield = horton_fc
        allocate(horton_fc(np))
        horton_fc = 0.0
        write(*,*)'Reading ',trim(fcfile)
        open(unit = 500, file = trim(fcfile), form = 'unformatted', access = 'stream')
        read(500)horton_fc
        close(500)
        !
        ! Horton: initial infiltration capacity (f0)
        allocate(horton_f0(np))
        horton_f0 = 0.0
        write(*,*)'Reading ',trim(f0file)
        open(unit = 501, file = trim(f0file), form = 'unformatted', access = 'stream')
        read(501)horton_f0
        close(501)
        !
        ! Prescribe the current estimate (for output only; initial capacity)
        qinffield = horton_f0/3600/1000
        !
        ! Empirical constant (1/hr) k => note that this is different than ks used in Curve Number and Green-Ampt
        allocate(horton_kd(np))
        horton_kd = 0.0
        write(*,*)'Reading ',trim(kdfile)
        open(unit = 502, file = trim(kdfile), form = 'unformatted', access = 'stream')
        read(502)horton_kd
        close(502)
        write(*,*)'Using constant recovery rate that is based on constant factor relative to ',trim(kdfile)
        !
        ! Estimate of time
        allocate(rain_T1(np))
        rain_T1 = 0.0
        !
    endif
    !            
    end subroutine read_infiltration_file_original
    
    
    subroutine read_infiltration_file_netcdf()
    !
    use sfincs_data
    !
    implicit none
    !
    !character*256, intent(in) :: netinfiltrationfile
    !    
    integer :: np_nc,ip, nm
    !
    integer*4                                       :: quadtree_nr_points
    !   
    integer*1,          dimension(:),   allocatable :: quadtree_mask
    integer*4, dimension(:),  allocatable :: z_index   
    real*4, dimension(:),     allocatable :: rtmpz   
    !
    write(*,*)'Reading Infiltration netCDF file on quadtree mesh...'
    !
    NF90(nf90_open(trim(netinfiltrationfile), NF90_CLOBBER, net_file_inf%ncid))
    !          
    ! Get dimensions id's: nr points  
    !
    NF90(nf90_inq_dimid(net_file_inf%ncid, "mesh2d_nFaces", net_file_inf%np_dimid))
    !
    ! Get dimensions sizes    
    !
    NF90(nf90_inquire_dimension(net_file_inf%ncid, net_file_inf%np_dimid, len = np_nc))   ! nr of cells
    !
    ! Get variable id's
    NF90(nf90_inq_varid(net_file_inf%ncid, 'mask',  net_file_inf%mask_varid))    
    !
    !NF90(nf90_inq_varid(net_file_inf%ncid, 'type',  net_file_inf%type_varid))                
    inftype = 'cna' !TL: later determine this based on the netcdf file
    !
    ! Allocate variables
    allocate(z_index(np_nc))   
    allocate(rtmpz(np_nc))    
    !allocate(kcs(np))
    allocate(quadtree_mask(np))        
    !
    ! Netcdf quadtree infiltrationfile contains data for the entire quadtree. So also for points with kcs==0 !
    ! This means that the data needs to be re-mapped to the active cell indices:
    if (use_quadtree) then
        !
        do ip = 1, np_nc !np_nc~quadtree_nr_points
            !
            nm = index_sfincs_in_quadtree(ip)
            !
            if (nm>0) then
                z_index(nm) = ip  
            endif 
            !    
        enddo
        !
    !else
    !    ! Regular grid > Tl: only keep regular if we want to support this in the end
    !    do nm = 1, np_nc
    !        z_index(nm) = nm
    !    enddo
    endif
    !
    ! Read Z points
    !    
    NF90(nf90_get_var(net_file_inf%ncid, net_file_inf%mask_varid,  rtmpz(:)))    
    !
    do ip = 1, np
        quadtree_mask(ip) = rtmpz(z_index(ip))
    enddo       
    !
    ! Check whether number of cells matches
    !if (np_nc /= np) then        
    !    write(*,'(a,i8,a,i8,a)')'Error! Number of cells in netinfiltrationfile ',np_nc,' does not match number of active cells in mesh ',np,' ! Abort reading infiltration input...',np_nc
    !    continue        
    !endif 
    !
    ! Check whether incoming mask is matches the one as retrieved from netcdf quadtree grid file (our general 'msk')
    if (all(quadtree_mask == kcs)) then
        print *, "The mask arrays are the same."
    else
        write(*,*)'Error! Mask of netinfiltrationfile does not match with actiave mask in mesh! Abort reading infiltration input...'
        continue
    endif    
    !
    ! Now read in infiltration type specific variables
    !
    if (inftype == 'cna') then
        !
        ! Spatially-varying infiltration with CN numbers (old)
        !
        write(*,*)'Turning on process: Infiltration (via Curve Number method - A)'               
        !
        ! Allocate qinffield
        allocate(qinffield(np))
        !
        ! Read variables for cna method
        !
        ! Get variable id's
        NF90(nf90_inq_varid(net_file_inf%ncid, 'scs',  net_file_inf%scs_varid))     
        ! Read Z points
        !    
        NF90(nf90_get_var(net_file_inf%ncid, net_file_inf%scs_varid,  rtmpz(:)))    
        do ip = 1, np
            qinffield(ip) = rtmpz(z_index(ip))
        enddo          
        !
        ! already convert qinffield from inches to m here
        do nm = 1, np
        qinffield(nm) = qinffield(nm)*0.0254   !to m
        enddo           
        !
    elseif (inftype == 'cnb') then
        !
        ! Spatially-varying infiltration with CN numbers (with recovery)
        !
        write(*,*)'Turning on process: Infiltration (via Curve Number method - B)'            
        ! 
        ! Allocate Smax, Se, Ks
        allocate(qinffield(np))
        allocate(scs_Se(np))
        allocate(ksfield(np))
        allocate(inf_kr(np))        
        !
        ! Read variables for cnb method
        !
        ! Get variable id's
        NF90(nf90_inq_varid(net_file_inf%ncid, 'qinffield',  net_file_inf%qinffield_varid))     
        NF90(nf90_inq_varid(net_file_inf%ncid, 'scs_se',  net_file_inf%scs_se_varid))     
        NF90(nf90_inq_varid(net_file_inf%ncid, 'ksfield',  net_file_inf%ksfield_varid))             
        ! Read Z points
        !    
        NF90(nf90_get_var(net_file_inf%ncid, net_file_inf%qinffield_varid,  rtmpz(:)))    
        do ip = 1, np
            qinffield(ip) = rtmpz(z_index(ip))
        enddo  
        !    
        NF90(nf90_get_var(net_file_inf%ncid, net_file_inf%scs_se_varid,  rtmpz(:)))    
        do ip = 1, np
            scs_Se(ip) = rtmpz(z_index(ip))
        enddo  
        !    
        NF90(nf90_get_var(net_file_inf%ncid, net_file_inf%ksfield_varid,  rtmpz(:)))    
        do ip = 1, np
            ksfield(ip) = rtmpz(z_index(ip))
        enddo          
        !
        ! Compute recovery                     ! Equation 4-36
        inf_kr = sqrt(ksfield/25.4) / 75       ! Note that we assume ksfield to be in mm/hr, convert it here to inch/hr (/25.4)
        !                                    ! /75 is conversion to recovery rate (in days)
        !
        ! Allocate support variables
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
    endif
    !
    infiltration = .true.   
    !    
    NF90(nf90_close(net_file_inf%ncid))       
    !
    ! Deallocate  vars
    deallocate(rtmpz)
    deallocate(z_index)
    deallocate(quadtree_mask)   
    !         
    end subroutine read_infiltration_file_netcdf
    
    
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
   real*4  :: hh_local, a
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
      !
      !$omp parallel &
      !$omp private ( Qq,I,nm )       
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
               scs_S1(nm)          = scs_Se(nm)        ! S for this 'event'
               scs_rain(nm)        = 1                 ! logic used to determine if there is an event ongoing
            endif
            ! 
            !  Compute cum rainfall
            ! 
            scs_P1(nm) = scs_P1(nm) + prcp(nm)*dt
            ! 
            ! Compute runoff
            ! 
            if (scs_P1(nm) > (sfacinf* scs_S1(nm)) ) then ! scs_S1 is S
                Qq              = (scs_P1(nm) - (sfacinf*scs_S1(nm)))**2 / (scs_P1(nm) + (1.0 - sfacinf)*scs_S1(nm))  ! cumulative runoff in m
                I              = scs_P1(nm) - Qq                       ! cum infiltration this event
                qinfmap(nm)    = (I - scs_F1(nm))/dt                   ! infiltration in m/s
                scs_F1(nm)     = I                                     ! cum infiltration this event
            else
                Qq              = 0.0                                  ! no runoff
                scs_F1(nm)     = scs_P1(nm)                            ! all rainfall is infiltrated
                qinfmap(nm)    = prcp(nm)                              ! infiltration rate = rainfall rate
            endif
            ! 
            ! Compute "remaining S", but note that scs_Se is not used in computation
            ! 
            scs_Se(nm)      = scs_Se(nm) - qinfmap(nm)*dt
            scs_Se(nm)      = max(scs_Se(nm), 0.0)
            qinfmap(nm)     = max(qinfmap(nm), 0.0)
            !
         else
            ! 
            ! It is not raining here
            if (scs_rain(nm) == 1) then
               !
               ! if it was raining before; cange logic and set rate to 0
               scs_rain(nm)   = 0
               qinfmap(nm)    = 0.0
               rain_T1(nm)    = 0.0
               !
            endif
            !
            ! Add to recovery time
            rain_T1(nm)    = rain_T1(nm) + dt / 3600
            !
            ! compute recovery of S if time is larger than this
            if (rain_T1(nm) >  (0.06 / inf_kr(nm)) ) then			! Equation 4-37 from SWMM
                ! note that scs_Se is S and qinffield is Smax
                scs_Se(nm)  = scs_Se(nm) + (inf_kr(nm) * qinffield(nm) * dt / 3600)  ! scs_kr is recovery in hours 
                scs_Se(nm)  = min(scs_Se(nm), qinffield(nm))
                !
            endif
            !
         endif
         ! 
         ! Compute cumulative values
         cuminf(nm)     = cuminf(nm) + qinfmap(nm)*dt
         netprcp(nm)    = netprcp(nm) - qinfmap(nm)
         !
      enddo
      !$omp end do
      !$omp end parallel 
      !
      !$acc update device(qinfmap), async(1)      
      !
   elseif (inftype == 'gai') then
      !
      ! Determine infiltration rate with  with the Green-Ampt (GA) model
      !$omp parallel &
      !$omp private ( nm )
      !$omp do              
      do nm = 1, np
         !
         ! If there is precip in this grid cell for this time step?
         if (prcp(nm) > 0.0) then
            !
            ! Is raining now
            if ( (prcp(nm)) < ksfield(nm)) then
                !
                ! Small amounts of rainfall - infiltration is same as soil
                qinfmap(nm)    = prcp(nm)                       ! infiltration is same as rainfall
                !
            else
                !
                ! Larger amounts of rainfall - Equation 4-27 from SWMM manual
                qinfmap(nm)    = (ksfield(nm) * (1 +  (GA_head(np)*GA_sigma(np)) / GA_F(nm)))
                qinfmap(nm)    = min(qinfmap(nm), prcp(nm))     ! never more than rainfall
                qinfmap(nm)    = max(qinfmap(nm), 0.0)          ! and never negative
                !
            endif
            !
            ! Update sigma 
            GA_sigma(nm) = GA_sigma(nm) - (qinfmap(nm)*dt/GA_Lu(nm))
            GA_sigma(nm) = max(GA_sigma(nm),0.00)
            ! 
            ! Update others
            GA_F(nm)    = GA_F(nm) + qinfmap(nm)*dt     ! internal cumulative rainfall from Green-Ampt
            rain_T1(nm) = 0.0                           ! recovery time not started
            !
         else
            ! 
            ! Not raining here
            !
            ! Add to recovery time
            rain_T1(nm)     = rain_T1(nm) + dt / 3600      
            !
            ! compute recovery of S if time is larger than this
            if (rain_T1(nm) >  (0.06 / inf_kr(nm)) ) then			! Equation 4-37 from SWMM
                ! 
                ! Update sigma 
                GA_sigma(nm) = GA_sigma(nm) + (inf_kr(nm) * GA_sigma_max(nm) * dt/3600)         ! Equation 4-35
                GA_sigma(nm) = min(GA_sigma(nm), GA_sigma_max(nm))                              ! never more than max
                !
                ! Update internal cumulative rainfall
                GA_F(nm)    = GA_F(nm) - (inf_kr(nm) * GA_sigma_max(nm) * dt/3600*GA_Lu(nm))    ! Page 112 SWMM
                GA_F(nm)    = max(GA_F(nm), 0.0)                                                ! never negative
            endif
         endif
         ! 
         ! Compute cumulative values
         qinffield(nm)  = qinfmap(nm)
         cuminf(nm)     = cuminf(nm) + qinfmap(nm)*dt
         netprcp(nm)    = netprcp(nm) - qinfmap(nm)
         !
      enddo
      !$omp end do
      !$omp end parallel       
      !
      !$acc update device(qinfmap), async(1)        
      !
   elseif (inftype == 'hor') then
      !
      ! Determine infiltration rate with with the Horton model
      !$omp parallel &
      !$omp private  ( Qq,I,nm,a,hh_local)       
      !$omp do              
      do nm = 1, np
         !
         ! Get local water depth estimate
         if (subgrid) then
             if (crsgeo) then
                a = cell_area_m2(nm)
             else   
                a = cell_area(z_flags_iref(nm))
             endif
             hh_local = z_volume(nm)/a
         else
             hh_local = zs(nm) - zb(nm)
         endif
         !
         ! Check if there is water
          if (hh_local> 0.0 .OR. prcp(nm) > 0.0) then
            !
            ! Infiltrating here
            !
            ! Count how long this is already going
            if (rain_T1(nm) > 0) then
                rain_T1 = 0.0
            endif
            rain_T1(nm)     = rain_T1(nm) - dt                                              ! negative amount of how long it is infiltrating
            ! 
            ! Compute estimate of infiltration                                              ! Note that qinffield = horton_fc
            I               = exp(horton_kd(nm)*rain_T1(nm)/3600)                             ! note that horton_kd is factor in hours while dt is seconds
            !
            ! Stop keeping track of this when less than 1% left (same for time)
            if (I < 0.01) then
                I              = 0.0                        ! which reduces qinfmap to horton_fc
                rain_T1(nm)     = rain_T1(nm) + dt          ! also make sure time doesnt further decrease
                qinfmap(nm)     = (horton_fc(nm))/3600/1000 ! from mm/hr to m/s
            else
                qinfmap(nm)     = (horton_fc(nm) + (horton_f0(nm) - horton_fc(nm))*I)/3600/1000 ! from mm/hr to m/s
            endif
            !
            ! Check how much there can infiltrate
            if (hh_local > 0.00) then
                Qq              = prcp(nm)*dt + (zs(nm) - zb(nm))                           ! Qq is estimate in meter of how much water there is
            else
                Qq              = prcp(nm)*dt                                               ! if no water; only compare with rainfall
            endif
            !
            ! Compare how much Horton want to infiltrate
            I               = qinfmap(nm) * dt                                              ! I is estimate in meter of how much Horton allows
            if (I > Qq) then
                qinfmap(nm) = qinfmap(nm) * Qq/I                                            ! scale Horton if capacity > available
            endif
            !
            else
                !
                ! Not raining here NOR ponding
                rain_T1(nm)     = rain_T1(nm) + dt/horton_kr_kd                                 ! positive amount of how long it is infiltrating
                qinfmap(nm)     = 0.0
                !
            endif
         ! 
         ! Compute cumulative values
         cuminf(nm)     = cuminf(nm) + qinfmap(nm)*dt
         netprcp(nm)    = netprcp(nm) - qinfmap(nm)
         !
      enddo
      !$omp end do
      !$omp end parallel 
      !
   endif
   !
   end subroutine   

   
   subroutine handle_err(status,file,line)
      !
      integer, intent ( in)    :: status
      character(*), intent(in) :: file
      integer, intent ( in)    :: line
      integer :: status2
      !   
      if(status /= nf90_noerr) then
         write(0,'("NETCDF ERROR: ",a,i6,":",a)') file,line,trim(nf90_strerror(status))
      end if
   end subroutine handle_err
   !   
end module
