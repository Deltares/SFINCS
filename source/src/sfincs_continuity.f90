module sfincs_continuity

contains
   
   subroutine compute_water_levels(t, dt, tloop)
   !
   use sfincs_data
   !
   implicit none
   !
   real*4           :: dt
   real*8           :: t
   !
   integer          :: count0
   integer          :: count1
   integer          :: count_rate
   integer          :: count_max
   real             :: tloop
   !
   call system_clock(count0, count_rate, count_max)
   !
   if (subgrid) then
      !
      call compute_water_levels_subgrid(dt,t)
      !
   else
      !
      call compute_water_levels_regular(dt,t)
      !
   endif  
   !
   ! Put non-default store options in a separate subroutine (all but zsmax) to save computation time when not used (both regular and subgrid):
   !
   if ((store_maximum_velocity .eqv. .true.) .or. (store_maximum_flux .eqv. .true.) .or. (store_twet .eqv. .true.)) then    
      ! 
      call compute_store_variables(dt)       
      !    
   endif
   !   
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine

   
   subroutine compute_water_levels_regular(dt,t)
   !
   use sfincs_data
   !
   implicit none
   !
   real*4           :: dt
   real*8           :: t   
   !
   integer          :: nm
   integer          :: isrc
   !
   integer          :: iwm
   !
   integer          :: nmu
   integer          :: nmd
   integer          :: num
   integer          :: ndm
   !
   real*4           :: qnmu
   real*4           :: qnmd
   real*4           :: qnum
   real*4           :: qndm
   real*4           :: factime
   real*4           :: dvol    
   !
   if (snapwave) then ! need to compute filtered water levels for snapwave
      !
      factime = min(dt/wmtfilter, 1.0)
      !
   endif
   !
   !$acc parallel present( kcs, zs, zb, netprcp, prcp, q, qext, zsmax, zsm, maxzsm, &
   !$acc                   z_flags_iref, uv_flags_iref, &
   !$acc                   z_index_uv_md, z_index_uv_nd, z_index_uv_mu, z_index_uv_nu, &
   !$acc                   dxm, dxrm, dyrm, dxminv, dxrinv, dyrinv, cell_area_m2, cell_area,  &
   !$acc                   nmindsrc, qtsrc, &
   !$acc                   z_index_wavemaker, wavemaker_uvmean, wavemaker_nmd, wavemaker_nmu, wavemaker_ndm, wavemaker_num )
   !
   ! First discharges (don't do this parallel, as it's probably not worth it)
   !
   if (nsrcdrn > 0) then
      ! 
      !$acc loop
      do isrc = 1, nsrcdrn
         ! 
         nm = nmindsrc(isrc)
         ! 
         if (crsgeo) then
            ! 
            zs(nmindsrc(isrc)) = max(zs(nm) + qtsrc(isrc) * dt / cell_area_m2(nm), zb(nm))
            ! 
         else
            ! 
            zs(nmindsrc(isrc)) = max(zs(nm) + qtsrc(isrc) * dt / cell_area(z_flags_iref(nm)), zb(nm))
            ! 
         endif
         ! 
      enddo
      ! 
   endif
   !
   !$omp parallel &
   !$omp private ( nm,dvol,nmd,nmu,ndm,num,qnmd,qnmu,qndm,qnum,iwm)
   !$omp do schedule ( dynamic, 256 )
   !$acc loop gang vector
   do nm = 1, np
      ! 
      if (kcs(nm) == 1) then ! Regular point
         !
         if (precip) then
            !
            zs(nm) = zs(nm) + netprcp(nm) * dt
            !
         endif
         !
         if (use_qext) then
            !
            ! Add external source (e.g. from XMI coupling) 
            ! 
            zs(nm) = zs(nm) + qext(nm) * dt 
            !
         endif
         !
         nmd = z_index_uv_md(nm)
         nmu = z_index_uv_mu(nm)
         ndm = z_index_uv_nd(nm)
         num = z_index_uv_nu(nm)
         !
         if (crsgeo) then
            !
            ! Use cell width dxm (which varies with latitude)
            !
            zs(nm)   = zs(nm) + ( (q(nmd) - q(nmu)) / dxm(nm) + (q(ndm) - q(num)) * dyrinv(z_flags_iref(nm)) ) * dt
            !
            ! Should really be:
            !
            ! zs(nm)   = zs(nm) + ( (q(nmd) - q(nmu)) / dxm(nm) + (q(ndm) * f - q(num) / f) * dyrinv(z_flags_iref(nm)) ) * dt
            !
            ! Where f if a correction factor for the ratio between cell width at cell centre and cell bottom:
            !
            ! f = cos(y - 0.5*dy) / cos(y) (this holds both in northern and southern hemisphere)
            ! 
         else   
            !
            zs(nm)   = zs(nm) + ( (q(nmd) - q(nmu)) * dxrinv(z_flags_iref(nm)) + (q(ndm) - q(num)) * dyrinv(z_flags_iref(nm)) ) * dt
            !
         endif
         !
      endif
      !
      if (wavemaker) then
         !      
         if (kcs(nm) == 4) then
            !
            ! Wave maker point (seaward of wave maker)
            ! Here we use the mean flux at the location of the wave maker 
            !
            iwm = z_index_wavemaker(nm)
            !
            if (wavemaker_nmd(iwm) > 0) then
               !
               ! Wave paddle on the left
               !
               qnmd = wavemaker_uvmean(wavemaker_nmd(iwm))
               !
            else
               !
               qnmd = q(z_index_uv_md(nm))
               !
            endif   
            !
            if (wavemaker_nmu(iwm) > 0) then
               !
               ! Wave paddle on the right
               !
               qnmu = wavemaker_uvmean(wavemaker_nmu(iwm))
               !
            else
               !
               qnmu = q(z_index_uv_mu(nm))
               !
            endif   
            !
            if (wavemaker_ndm(iwm) > 0) then
               !
               ! Wave paddle below
               !
               qndm = wavemaker_uvmean(wavemaker_ndm(iwm))
               !
            else
               !
               qndm = q(z_index_uv_nd(nm))
               !
            endif   
            !
            if (wavemaker_num(iwm) > 0) then
               !
               ! Wave paddle above
               !
               qnum = wavemaker_uvmean(wavemaker_num(iwm))
               !
            else
               !
               qnum = q(z_index_uv_nu(nm))
               !
            endif   
            !
            zs(nm)   = zs(nm) + (((qnmd - qnmu) * dxrinv(z_flags_iref(nm)) + (qndm - qnum) * dyrinv(z_flags_iref(nm)))) * dt
            !
         endif   
         !
      endif
      !
      if (snapwave) then
         !
         ! Time-averaged water level used for SnapWave using exponential filter
         !
         ! Would double exponential filtering be better?
         !
         zsm(nm) = factime * zs(nm) + (1.0 - factime) * zsm(nm)
         !
         if (store_maximum_waterlevel) then
             !
             maxzsm(nm) = max(maxzsm(nm), zsm(nm))             
             !
         endif         
         !             
      endif 
      !
      ! No continuity update but keeping track of variables         
      ! zsmax used by default, therefore keep in standard continuity loop:
      !
      if (store_maximum_waterlevel) then
         !
         ! Store when the maximum water level changed
         !
         if (store_t_zsmax) then
             if (zs(nm) > zsmax(nm)) then
                 if ( (zs(nm) - zb(nm)) > huthresh) then
                    t_zsmax(nm) = t
                 endif
             endif
         endif
         !
         ! Store the maximum water level itself
         !
         zsmax(nm) = max(zsmax(nm), zs(nm))
         !
      endif
      !
   enddo
   !$omp end do
   !$omp end parallel
   !$acc end parallel
   !         
   end subroutine

   
   subroutine compute_water_levels_subgrid(dt,t)
   !
   use sfincs_data
   !
   implicit none
   !
   real*4           :: dt
   real*8           :: t
   !
   integer          :: nm
   integer          :: isrc
   !
   integer          :: iwm
   integer          :: ind
   !
   integer          :: nmu
   integer          :: nmd
   integer          :: num
   integer          :: ndm
   !
   real*4           :: factime
   real*4           :: dvol  
   real*4           :: dzsdt
   !
   real*4           :: qnmu
   real*4           :: qnmd
   real*4           :: qnum
   real*4           :: qndm
   !
   integer          :: iuv
   real*4           :: dzvol
   real*4           :: facint
   real*4           :: a
   real*4           :: dv
   real*4           :: zs00
   real*4           :: zs11
   !
   if (wavemaker) then
      !
      factime = min(dt / wmtfilter, 1.0)
      !
   endif   
   !
   ! First discharges (don't do this parallel, as it's probably not worth it)
   ! NVFORTAN turns this into a sequential loop (!$acc loop seq)
   !
   if (nsrcdrn > 0) then
      !
      !$acc serial present( z_volume, nmindsrc, qtsrc )
      do isrc = 1, nsrcdrn
         !
         nm = nmindsrc(isrc)
         !
         if ((z_volume(nm) >= 0) .or. ((qtsrc(isrc)<0.0) .and. (z_volume(nm) >= 0))) then
            z_volume(nm) = z_volume(nm) + qtsrc(isrc) * dt
         endif
         !
      enddo
      !$acc end serial
      !
   endif
   !
   !$omp parallel &
   !$omp private ( dvol,dzsdt,nmd,nmu,ndm,num,a,iuv,facint,dzvol,ind,iwm,qnmd,qnmu,qndm,qnum,dv,zs00,zs11 )
   !$omp do schedule ( dynamic, 256 )
   !$acc parallel present( kcs, zs, zs0, zb, z_volume, zsmax, zsm, maxzsm, zsderv, &
   !$acc                   subgrid_z_zmin,  subgrid_z_zmax, subgrid_z_dep, subgrid_z_volmax, &
   !$acc                   netprcp, prcp, q, qext, z_flags_iref, uv_flags_iref, &
   !$acc                   z_index_uv_md, z_index_uv_nd, z_index_uv_mu, z_index_uv_nu, &
   !$acc                   dxm, dxrm, dyrm, dxminv, dxrinv, dyrinv, cell_area_m2, cell_area, &   
   !$acc                   z_index_wavemaker, wavemaker_uvmean, wavemaker_nmd, wavemaker_nmu, wavemaker_ndm, wavemaker_num, storage_volume)
   !$acc loop gang vector
   do nm = 1, np
      !
      ! And now water level changes due to horizontal fluxes
      !
      dvol = 0.0
      !
      if (kcs(nm) == 1) then
         !
         nmd = z_index_uv_md(nm)
         nmu = z_index_uv_mu(nm)
         ndm = z_index_uv_nd(nm)
         num = z_index_uv_nu(nm)
         !
         if (crsgeo) then
            !
            ! dxm  = size of cell in x - direction (it varies for all cells)
            ! dyrm = size of cell in y - direction (it varies for all zoom levels)
            !
            dvol = dvol + ( (q(nmd) - q(nmu)) * dyrm(z_flags_iref(nm)) + (q(ndm) - q(num)) * dxm(nm) ) * dt
            !
            ! Should really be:
            !
            ! dvol = dvol + ( (q(nmd) - q(nmu)) * dyrm(z_flags_iref(nm)) + (q(ndm) * f - q(num) / f) * dxm(nm) ) * dt
            !
            ! where f if a correction factor for the ratio between cell width at cell centre and cell bottom:
            !
            ! f = cos(y - 0.5*dy) / cos(y) (this holds both in northern and southern hemisphere)
            !
            ! This assumes that we can use the same factor f for q(ndm) and q(num), i.e.:
            !
            ! cos(y - 0.5*dy) / cos(y) ~= cos(y + 0.5*dy) / cos(y) or: cos(y - 0.5*dy) ~= cos(y + 0.5*dy) which is pretty much true for dy < 1.0 degree
            !
         else
            !
            if (use_quadtree) then   
               !
               ! dxrm = size of cell in x - direction (it varies for all zoom levels)
               ! dyrm = size of cell in y - direction (it varies for all zoom levels)
               !               
               dvol = dvol + ( (q(nmd) - q(nmu)) * dyrm(z_flags_iref(nm)) + (q(ndm) - q(num)) * dxrm(z_flags_iref(nm)) ) * dt
               ! 
            else
               ! 
               dvol = dvol + ( (q(nmd) - q(nmu)) * dy + (q(ndm) - q(num)) * dx ) * dt
               ! 
            endif   
            !
         endif   
      endif ! kcs==1
      !
      if (wavemaker .and. kcs(nm) == 4) then
         !
         ! Wave maker point (seaward of wave maker)
         ! Here we use the mean flux at the location of the wave maker
         !
         iwm = z_index_wavemaker(nm)
         !
         if (wavemaker_nmd(iwm) > 0) then
            !
            ! Wave paddle on the left
            !
            qnmd = wavemaker_uvmean(wavemaker_nmd(iwm))
            !
         else
            !
            qnmd = q(z_index_uv_md(nm))
            !
         endif   
         !
         if (wavemaker_nmu(iwm) > 0) then
            !
            ! Wave paddle on the right
            !
            qnmu = wavemaker_uvmean(wavemaker_nmu(iwm))
            !
         else
            !
            qnmu = q(z_index_uv_mu(nm))
            !
         endif   
         !
         if (wavemaker_ndm(iwm) > 0) then
            !
            ! Wave paddle below
            !
            qndm = wavemaker_uvmean(wavemaker_ndm(iwm))
            !
         else
            !
            qndm = q(z_index_uv_nd(nm))
            !
         endif   
         !
         if (wavemaker_num(iwm) > 0) then
            !
            ! Wave paddle above
            !
            qnum = wavemaker_uvmean(wavemaker_num(iwm))
            !
         else
            !
            qnum = q(z_index_uv_nu(nm))
            !
         endif   
         !
         if (use_quadtree) then
            ! 
            dvol = dvol + ( (qnmd - qnmu) * dyrm(z_flags_iref(nm)) + (qndm - qnum) * dxrm(z_flags_iref(nm)) ) * dt
            ! 
         else
            ! 
            dvol = dvol + ( (qnmd - qnmu) * dy + (qndm - qnum) * dx ) * dt
            !
         endif   
         !
      endif
      !
      ! We got the volume change dvol in each active cell from fluxes
      ! Now first add precip and qext
      ! Then adjust for storage volume
      ! Then update the volume and compute new water level           
      !
      if (kcs(nm) == 1 .or. kcs(nm) == 4) then
         !
         ! Obtain cell area
         !
         if (crsgeo) then
            !
            a = cell_area_m2(nm)
            !
         else   
            !
            a = cell_area(z_flags_iref(nm))
            !
         endif
         !
         if (precip .or. use_qext) then
            !
            dzsdt = 0.0
            !   
            if (precip) then
               !
               ! Add nett rainfall 
               !
               dzsdt = dzsdt + netprcp(nm)
               !
            endif
            !
            if (use_qext) then
               !
               ! Add external source (e.g. from XMI coupling) 
               ! 
               dzsdt = dzsdt + qext(nm)
               !
            endif
            !
            ! dzsdt is still in m/s, so multiply with a * dt to get m^3
            !
            dvol = dvol + dzsdt * a * dt         
            !      
         endif
         !
         if (use_storage_volume) then
            !
            ! If water enters the cell through a point discharge, it will NOT end up in storage volume !  
            !   
            if (storage_volume(nm) > 1.0e-6 .and. dvol > 0.0) then
               !
               ! There is still some storage left, and water is entering the cell
               !
               ! Compute remaining storage volume
               !
               dv = storage_volume(nm) - dvol
               !
               ! Update storage volume (it cannot become negative)) 
               ! 
               storage_volume(nm) = max(dv, 0.0)
               !
               if (dv < 0.0) then
                  !
                  ! Overshoot, so add remaining volume to z_volume
                  !
                  dvol = - dv
                  ! 
               else                   
                  !
                  ! Everything went into storage  
                  ! 
                  dvol = 0.0 
                  !
               endif
               !
            endif
            !
         endif
         !
         ! Update volume 
         !
         z_volume(nm) = z_volume(nm) + dvol
         !
         if (wiggle_suppression) then
            !
            ! Store previous water level to determine gradient
            ! 
            zs00    = zs0(nm) ! previous time step
            zs11    = zs(nm)  ! current time step before updating
            zs0(nm) = zs11    ! next previous time step
            ! 
         endif
         ! 
         ! Obtain new water level from subgrid tables 
         !
         if (z_volume(nm) >= subgrid_z_volmax(nm) * 0.999) then
            !
            ! Entire cell is wet, no interpolation needed
            !
            zs(nm) = max(subgrid_z_zmax(nm), -20.0) + (z_volume(nm) - subgrid_z_volmax(nm)) / a
            !
         elseif (z_volume(nm) <= 1.0e-6) then
            !
            ! No water in this cell. Set zs to z_zmin.
            !
            zs(nm) = max(subgrid_z_zmin(nm), -20.0)
            !
         else
            !
            ! Interpolation from subgrid tables needed.
            !
            dzvol    = subgrid_z_volmax(nm) / (subgrid_nlevels - 1)
            iuv      = int(z_volume(nm) / dzvol) + 1
            facint   = (z_volume(nm) - (iuv - 1) * dzvol ) / dzvol
            zs(nm)   = subgrid_z_dep(iuv, nm) + (subgrid_z_dep(iuv + 1, nm) - subgrid_z_dep(iuv, nm)) * facint
            !
         endif
         !
         !
         if (wiggle_suppression) then 
            ! 
            zsderv(nm) = zs(nm) - 2 * zs11 + zs00
            ! 
         endif
         !
      endif
      !
      if (snapwave) then
         !
         ! Time-averaged water level used for SnapWave using exponential filter
         !
         ! Would double exponential filtering be better?
         !
         zsm(nm) = factime * zs(nm) + (1.0 - factime) * zsm(nm)
         !
         if (store_maximum_waterlevel) then
             !
             maxzsm(nm) = max(maxzsm(nm), zsm(nm))             
             !
         endif         
         !    
      endif 
      !
      ! No continuity update but keeping track of variables         
      ! zsmax used by default, therefore keep in standard continuity loop:         
      !
      if (store_maximum_waterlevel) then
         !
         ! Store when the maximum water level changed
         !
         if (store_t_zsmax) then
             if (zs(nm) > zsmax(nm)) then
                 if ( (zs(nm) - subgrid_z_zmin(nm)) > huthresh) then
                    t_zsmax(nm) = t
                 endif
             endif
         endif
         !
         ! Store the maximum water level itself
         !
         zsmax(nm) = max(zsmax(nm), zs(nm))
         !
      endif
      !
   enddo
   !$omp end do
   !$omp end parallel
   !         
   !$acc end parallel
   !         
   end subroutine
   
   subroutine compute_store_variables(dt)
   !
   use sfincs_data
   !
   implicit none
   !
   real*4           :: dt
   !   
   integer          :: nm
   !   
   integer          :: nmu
   integer          :: nmd
   integer          :: num
   integer          :: ndm
   !
   real*4           :: quz
   real*4           :: qvz
   real*4           :: qz
   real*4           :: uvz      
   !
   !$omp parallel &
   !$omp private ( nmd, nmu, ndm, num, quz, qvz, qz, uvz )
   !$omp do schedule ( dynamic, 256 )
   !$acc parallel present( kcs, zs, zb, subgrid_z_zmin, q, uv, vmax, qmax, twet, &
   !$acc                   z_index_uv_md, z_index_uv_nd, z_index_uv_mu, z_index_uv_nu )
   !$acc loop gang vector 
   do nm = 1, np
      !
      ! And now water level changes due to horizontal fluxes
      !
      qz = 0.0
      uvz = 0.0    
      !
      if (kcs(nm) == 1 .or. kcs(nm) == 4) then ! TL: kcs(nm)==4 also correct for regular?
         !      
         ! Regular point with four surrounding cells of the same size
         !
         nmd = z_index_uv_md(nm)
         nmu = z_index_uv_mu(nm)
         ndm = z_index_uv_nd(nm)
         num = z_index_uv_nu(nm)
         !
         if (store_maximum_velocity) then
            quz = (uv(nmd) + uv(nmu)) / 2   
            qvz = (uv(ndm) + uv(num)) / 2      
            uvz = sqrt(quz**2 + qvz**2)
         endif
         !
         if (store_maximum_flux) then
            quz = (q(nmd) + q(nmu)) / 2   
            qvz = (q(ndm) + q(num)) / 2      
            qz = sqrt(quz**2 + qvz**2)
         endif   
         !
         ! No continuity update but keeping track of variables
         ! 1. store vmax
         if (store_maximum_velocity) then
            !             
            vmax(nm) = max(vmax(nm), uvz)
            !
         endif      
         !
         ! 2. store qmax         
         if (store_maximum_flux) then
            !             
            qmax(nm) = max(qmax(nm), qz)
            !
         endif      
         !
         ! 3. store Twet
         if (store_twet) then
            if (subgrid) then
              !
              if ( (zs(nm) - subgrid_z_zmin(nm)) > twet_threshold) then
                 twet(nm) = twet(nm) + dt
              endif
              !
            else
              !
              if ( (zs(nm) - zb(nm)) > twet_threshold) then
                 !
                 twet(nm) = twet(nm) + dt
                 !
              endif
              !
            endif               
         endif         
         !
      endif   
   enddo   
   !$omp end do
   !$omp end parallel
   !$acc end parallel
   !       
   end subroutine
   
end module
