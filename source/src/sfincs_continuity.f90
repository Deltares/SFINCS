module sfincs_continuity

contains
   
   subroutine compute_water_levels(dt, tloop)
   !
   use sfincs_data
   !
   implicit none
   !
   real*4           :: dt
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
      call compute_water_levels_subgrid(dt)
      !
   else
      !
      call compute_water_levels_regular(dt)
      !
   endif  
   !
   ! Put non-default store options in a separate subroutine (all but zsmax) to save computation time when not used (both regular and subgrid):
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

   
   subroutine compute_water_levels_regular(dt)
   !
   use sfincs_data
   !
   implicit none
   !
   real*4           :: dt
   !
   integer          :: nm
   integer          :: isrc
   !
   integer          :: iwm
   integer          :: ileft
   integer          :: iright
   integer          :: itop
   integer          :: ibottom
   !
   integer          :: nmu
   integer          :: nmd
   integer          :: num
   integer          :: ndm
   integer          :: nmu1
   integer          :: nmd1
   integer          :: num1
   integer          :: ndm1
   integer          :: nmu2
   integer          :: nmd2
   integer          :: num2
   integer          :: ndm2
   !
   integer          :: ip
   !
   real*4           :: qnmu
   real*4           :: qnmd
   real*4           :: qnum
   real*4           :: qndm
   real*4           :: factime
   real*4           :: dvol    
   !
   real*4           :: zs0
   real*4           :: zsm0
   real*4           :: zsm1
   real*4           :: dzdt
   real*4           :: dzf 
   real*4           :: dzs
   real*4           :: zcor
   !
   if (snapwave) then ! need to compute filtered water levels for snapwave
      !
      factime = min(dt/wmtfilter, 1.0)
      !
   endif   
   !
   ! First discharges (don't do this parallel, as it's probably not worth it)
   ! Should try to do this in a smart way for openacc
   !
   if (nsrcdrn>0) then
      !$acc serial, present( zs,nmindsrc,qtsrc,zb,cell_area,z_flags_iref ), async(1)
      do isrc = 1, nsrcdrn
         nm = nmindsrc(isrc)
         zs(nmindsrc(isrc))   = max(zs(nm) + qtsrc(isrc)*dt/cell_area(z_flags_iref(nm)), zb(nm))
      enddo
      !$acc end serial
   endif   
   !
   !$omp parallel &
   !$omp private ( nm,dvol,nmd,nmu,ndm,num,qnmd,qnmu,qndm,qnum,iwm)
   !$omp do schedule ( dynamic, 256 )
   !$acc kernels present( kcs, zs, zb, netprcp, cumprcpt, prcp, q, zsmax, zsm, maxzsm, &
   !$acc                  z_flags_iref, uv_flags_iref, &
   !$acc                  z_index_uv_md, z_index_uv_nd, z_index_uv_mu, z_index_uv_nu, &
   !$acc                  dxm, dxrm, dyrm, dxminv, dxrinv, dyrinv, cell_area_m2, cell_area,  &
   !$acc                  z_index_wavemaker, wavemaker_uvmean, wavemaker_nmd, wavemaker_nmu, wavemaker_ndm, wavemaker_num), async(1)
   !$acc loop independent, private( nm )
   do nm = 1, np
      ! 
      if (kcs(nm)==1) then ! Regular point
         !
         if (precip) then
            !
            cumprcpt(nm) = cumprcpt(nm) + netprcp(nm)*dt
            !
            ! Add rain and/or infiltration only when cumulative effect over last interval exceeds 0.001 m
            ! Otherwise single precision may miss a lot of the rainfall/infiltration
            !
            if (cumprcpt(nm)>0.001 .or. cumprcpt(nm)<-0.001) then
               !
               zs(nm) = zs(nm) + cumprcpt(nm)
               cumprcpt(nm) = 0.0
!               zs(nm) = max(zs(nm), zb(nm)) ! don't allow negative water levels due to infiltration
               !
            endif
            !
         endif
         !
         nmd = z_index_uv_md(nm)
         nmu = z_index_uv_mu(nm)
         ndm = z_index_uv_nd(nm)
         num = z_index_uv_nu(nm)
         !
         if (crsgeo) then
            zs(nm)   = zs(nm) + (((q(nmd) - q(nmu))*dxminv(nmu) + (q(ndm) - q(num))*dyrinv(z_flags_iref(nm))))*dt
         else   
            zs(nm)   = zs(nm) + (((q(nmd) - q(nmu))*dxrinv(z_flags_iref(nm)) + (q(ndm) - q(num))*dyrinv(z_flags_iref(nm))))*dt
         endif   
         !
      endif
      !
      if (wavemaker) then
         !      
         if (kcs(nm)==4) then
            !
            ! Wave maker point (seaward of wave maker)
            ! Here we use the mean flux at the location of the wave maker 
            !
            iwm = z_index_wavemaker(nm)
            !
            if (wavemaker_nmd(iwm)>0) then
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
            if (wavemaker_nmu(iwm)>0) then
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
            if (wavemaker_ndm(iwm)>0) then
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
            if (wavemaker_num(iwm)>0) then
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
            zs(nm)   = zs(nm) + (((qnmd - qnmu)*dxrinv(z_flags_iref(nm)) + (qndm - qnum)*dyrinv(z_flags_iref(nm))))*dt
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
         zsm(nm) = factime*zs(nm) + (1.0 - factime)*zsm(nm)
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
         zsmax(nm) = max(zsmax(nm), zs(nm))
         !
      endif
      !
   enddo
   !$omp end do
   !$omp end parallel
   !$acc end kernels
   !$acc wait(1)
   !         
   end subroutine

   
   subroutine compute_water_levels_subgrid(dt)
   !
   use sfincs_data
   !
   implicit none
   !
   real*4           :: dt
   !
   integer          :: nm
   integer          :: isrc
   !
   integer          :: iwm
   integer          :: ileft
   integer          :: iright
   integer          :: itop
   integer          :: ibottom
   integer          :: ind
   !
   integer          :: nmu
   integer          :: nmd
   integer          :: num
   integer          :: ndm
   integer          :: nmu1
   integer          :: nmd1
   integer          :: num1
   integer          :: ndm1
   integer          :: nmu2
   integer          :: nmd2
   integer          :: num2
   integer          :: ndm2
   !
   integer          :: ip
   integer          :: n
   integer          :: m
   !
   real*4           :: qxnm
   real*4           :: qxnmd
   real*4           :: qynm
   real*4           :: qyndm
   real*4           :: factime
   real*4           :: dvol  
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
   real*4           :: uz
   real*4           :: vz
   real*4           :: dv
   !
   if (wavemaker) then
      !
      factime = min(dt/wmtfilter, 1.0)
      !
   endif   
   !
   ! First discharges (don't do this parallel, as it's probably not worth it)
   ! Should try to do this in a smart way for openacc
   !
   if (nsrcdrn>0) then
      !$acc serial, present( z_volume, nmindsrc, qtsrc ), async(1)
      do isrc = 1, nsrcdrn
         if (nmindsrc(isrc)>0) then ! should really let this happen
            z_volume(nmindsrc(isrc)) = max(z_volume(nmindsrc(isrc)) + qtsrc(isrc)*dt, 0.0)         
         endif   
      enddo
      !$acc end serial
   endif   
   !
   !$omp parallel &
   !$omp private ( dvol,nmd,nmu,ndm,num,a,iuv,facint,dzvol,ind,iwm,qnmd,qnmu,qndm,qnum,dv )
   !$omp do schedule ( dynamic, 256 )
   !$acc kernels present( kcs, zs, zb, z_volume, zsmax, zsm, maxzsm, &
   !$acc                  subgrid_z_zmin,  subgrid_z_zmax, subgrid_z_dep, subgrid_z_volmax, &
   !$acc                  netprcp, cumprcpt, prcp, q, z_flags_iref, uv_flags_iref, &
   !$acc                  z_index_uv_md, z_index_uv_nd, z_index_uv_mu, z_index_uv_nu, &
   !$acc                  dxm, dxrm, dyrm, dxminv, dxrinv, dyrinv, cell_area_m2, cell_area, &
   !$acc                  z_index_wavemaker, wavemaker_uvmean, wavemaker_nmd, wavemaker_nmu, wavemaker_ndm, wavemaker_num, storage_volume), async(1)
   !$acc loop independent, private( nm )
   do nm = 1, np
      !
      ! And now water level changes due to horizontal fluxes
      !
      dvol = 0.0
      !
      if (kcs(nm)==1) then
         !
         if (crsgeo) then
            a = cell_area_m2(nm)
         else   
            a = cell_area(z_flags_iref(nm))
         endif
         !
         if (precip) then
            !
            cumprcpt(nm) = cumprcpt(nm) + netprcp(nm)*dt
            !
            ! Add rain and/or infiltration only when cumulative effect over last interval exceeds 0.001 m
            ! Otherwise single precision may miss a lot of the rainfall/infiltration
            !
            if (cumprcpt(nm)>0.001 .or. cumprcpt(nm)<-0.001) then
               !
               dvol = dvol + cumprcpt(nm)*a
               cumprcpt(nm) = 0.0
               !
            endif
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
            dvol = dvol + ( (q(nmd) - q(nmu))*dyrm(uv_flags_iref(nm)) + (q(ndm) - q(num))*dxm(nm) ) * dt
            !
         else
            !
            if (use_quadtree) then   
               dvol = dvol + ( (q(nmd) - q(nmu))*dyrm(uv_flags_iref(nm)) + (q(ndm) - q(num))*dxrm(uv_flags_iref(nm)) ) * dt
            else
               dvol = dvol + ( (q(nmd) - q(nmu))*dy + (q(ndm) - q(num))*dx ) * dt
            endif   
            !
         endif   
      endif ! kcs==1    
      !
      if (wavemaker .and. kcs(nm)==4) then
         !
         ! Wave maker point (seaward of wave maker)
         ! Here we use the mean flux at the location of the wave maker 
         !
         iwm = z_index_wavemaker(nm)
         !
         if (wavemaker_nmd(iwm)>0) then
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
         if (wavemaker_nmu(iwm)>0) then
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
         if (wavemaker_ndm(iwm)>0) then
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
         if (wavemaker_num(iwm)>0) then
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
            dvol = dvol + ( (qnmd - qnmu)*dyrm(uv_flags_iref(nm)) + (qndm - qnum)*dxrm(uv_flags_iref(nm)) ) * dt
         else
            dvol = dvol + ( (qnmd - qnmu)*dy + (qndm - qnum)*dx ) * dt
         endif   
         !
      endif
      !
      ! We got the volume change dvol in each active cell
      ! Now update the volume and compute new water level           
      !
      if (kcs(nm) == 1 .or. kcs(nm) == 4) then 
         !      
         if (use_storage_volume) then
            !
            if (storage_volume(nm)>1.0e-6) then
               !
               ! Still some storage left
               !
               ! Compute remaining storage volume
               !
               dv = storage_volume(nm) - dvol
               !
               storage_volume(nm) =  max(dv, 0.0)
               !
               if (dv < 0.0) then
                  !
                  ! Overshoot, so add remaining volume to z_volume
                  !
                  z_volume(nm) = z_volume(nm) - dv
                  !
               endif
               !
            else
               ! 
               ! Storage is full 
               !
               z_volume(nm) = z_volume(nm) + dvol
               !
            endif
            !
         else
            !
            z_volume(nm) = z_volume(nm) + dvol
            !
         endif
         ! 
         ! Obtain new water level from subgrid tables 
         !
         if (z_volume(nm)>=subgrid_z_volmax(nm) - 1.0e-6) then
            !
            ! Entire cell is wet, no interpolation needed
            !
            zs(nm) = max(subgrid_z_zmax(nm), -20.0) + (z_volume(nm) - subgrid_z_volmax(nm))/a
            !
         elseif (z_volume(nm)<=1.0e-6) then
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
            iuv      = min(int(z_volume(nm)/dzvol) + 1, subgrid_nlevels - 1)
            facint   = (z_volume(nm) - (iuv - 1)*dzvol ) / dzvol
            zs(nm)   = subgrid_z_dep(iuv, nm) + (subgrid_z_dep(iuv + 1, nm) - subgrid_z_dep(iuv, nm))*facint                    
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
         zsm(nm) = factime*zs(nm) + (1.0 - factime)*zsm(nm)
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
      if (store_maximum_waterlevel) then
         !
         zsmax(nm) = max(zsmax(nm), zs(nm))
         !
      endif
      !
   enddo
   !$omp end do
   !$omp end parallel
   !$acc end kernels
   !$acc wait(1)
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
   !$acc kernels present( kcs, zs, zb, subgrid_z_zmin, q, uv, vmax, qmax, twet, &
   !$acc                  z_index_uv_md, z_index_uv_nd, z_index_uv_mu, z_index_uv_nu), async(2)
   !$acc loop independent, private( nm )   
   do nm = 1, np
      !
      ! And now water level changes due to horizontal fluxes
      !
      qz = 0.0
      uvz = 0.0    
      !
      if (kcs(nm)==1 .or. kcs(nm)==4) then ! TL: kcs(nm)==4 also correct for regular?
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
                 twet(nm) = twet(nm) + dt
              endif
              !
            endif               
         endif         
         !
      endif   
   enddo   
   !$omp end do
   !$omp end parallel
   !$acc end kernels
   !$acc wait(2)
   !       
   end subroutine
   
end module
