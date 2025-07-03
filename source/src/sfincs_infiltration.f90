module sfincs_infiltration

contains

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
!      !$acc parallel present( qinfmap, qinffield, z_volume, zs, zb, netprcp,cuminf ), async(1)
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
!      !$acc parallel present( qinfmap, qinffield, prcp, netprcp, cumprcp, cuminf ), async(1)
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
      !$acc parallel present( qinfmap, prcp, netprcp, cuminf, scs_rain, scs_Se, scs_P1, scs_F1, scs_S1, rain_T1, qinffield, inf_kr ), async(1)
!      !$acc parallel present( qinfmap, prcp, netprcp, cuminf, scs_rain, scs_Se, scs_P1, scs_F1, scs_S1, rain_T1 ), async(1)
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
!      !$acc                  ksfield, GA_head, GA_sigma, GA_sigma_max, GA_F, GA_Lu, inf_kr ), async(1)
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
!      !$acc                  horton_kd, horton_fc, horton_f0 ), async(1)
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
