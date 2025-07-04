module sfincs_infiltration

contains

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
   elseif (inftype == 'r2d') then
      !
      ! Run-off coefficient map
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
               qinfmap(nm) = 100.0
            endif
         else
            if (zs(nm)<=zb(nm)) then
               qinfmap(nm) = 100.0
            endif
         endif
         !
         cuminf(nm) = cuminf(nm) + (1.0-qinfmap(nm)/100.0)*netprcp(nm)*dt
         netprcp(nm) = netprcp(nm) - (1.0-qinfmap(nm)/100.0)*netprcp(nm)
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

end module
