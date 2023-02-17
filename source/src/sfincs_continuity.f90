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
!   real*8 :: t
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
      !$acc serial, present( zs,nmindsrc,qtsrc,zb ), async(1)
      do isrc = 1, nsrcdrn
         nm = nmindsrc(isrc)
         zs(nmindsrc(isrc))   = max(zs(nm) + qtsrc(isrc)*dt/cell_area(z_flags_iref(nm)), zb(nm))
      enddo
      !$acc end serial
   endif   
   !
   !$omp parallel &
   !$omp private ( nm,dvol,nmd1,nmu1,ndm1,num1,nmd2,nmu2,ndm2,num2,nmd,nmu,ndm,num,qnmd,qnmu,qndm,qnum,iwm )
   !$omp do schedule ( dynamic, 256 )
   !$acc kernels present( kcs, zs, zb, netprcp, cumprcpt, prcp, q, zsmax, twet, zsm, &
   !$acc                  z_flags_type, z_flags_iref, uv_flags_iref, &
   !$acc                  z_index_uv_md1, z_index_uv_md2, z_index_uv_nd1, z_index_uv_nd2, z_index_uv_mu1, z_index_uv_mu2, z_index_uv_nu1, z_index_uv_nu2, &
   !$acc                  dxm, dxrm, dyrm, dxminv, dxrinv, dyrinv, cell_area_m2, cell_area,  &
   !$acc                  z_index_wavemaker, wavemaker_uvmean, wavemaker_nmd, wavemaker_nmu, wavemaker_ndm, wavemaker_num), async(1)
   !$acc loop independent, private( nm )
   do nm = 1, np
      !
      ! And now water level changes due to horizontal fluxes
      !
      if (kcs(nm)==1) then
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
         if (z_flags_type(nm) == 0) then
            !
            ! Regular point with four surrounding cells of the same size
            !
            nmd = z_index_uv_md1(nm)
            nmu = z_index_uv_mu1(nm)
            ndm = z_index_uv_nd1(nm)
            num = z_index_uv_nu1(nm)
            !
            if (crsgeo) then
               zs(nm)   = zs(nm) + (((q(nmd) - q(nmu))*dxminv(nmu) + (q(ndm) - q(num))*dyrinv(z_flags_iref(nm))))*dt
            else   
               zs(nm)   = zs(nm) + (((q(nmd) - q(nmu))*dxrinv(z_flags_iref(nm)) + (q(ndm) - q(num))*dyrinv(z_flags_iref(nm))))*dt
            endif   
            !
         else
            !
            if (use_quadtree) then   
               !
               ! More complicated
               !
               nmd1 = z_index_uv_md1(nm)
               nmu1 = z_index_uv_mu1(nm)
               ndm1 = z_index_uv_nd1(nm)
               num1 = z_index_uv_nu1(nm)
               nmd2 = z_index_uv_md2(nm)
               nmu2 = z_index_uv_mu2(nm)
               ndm2 = z_index_uv_nd2(nm)
               num2 = z_index_uv_nu2(nm)
               !
               dvol = 0.0
               !
               ! U
               !
               if (nmd1>0) then
                  dvol = dvol + q(nmd1)*dyrm(uv_flags_iref(nmd1))
               endif   
               !
               if (nmd2>0) then
                  dvol = dvol + q(nmd2)*dyrm(uv_flags_iref(nmd2))
               endif
               !
               if (nmu1>0) then
                  dvol = dvol - q(nmu1)*dyrm(uv_flags_iref(nmu1))
               endif   
               !
               if (nmu2>0) then
                  dvol = dvol - q(nmu2)*dyrm(uv_flags_iref(nmu2))
               endif
               !
               ! V
               !
               if (crsgeo) then
                  !
                  if (ndm1>0) then                     
                     dvol = dvol + q(ndm1)*dxm(ndm1)
                  endif   
                  !
                  if (ndm2>0) then
                     dvol = dvol + q(ndm2)*dxm(ndm2)
                  endif 
                  !
                  if (num1>0) then
                     dvol = dvol - q(num1)*dxm(num1)
                  endif   
                  !
                  if (num2>0) then
                     dvol = dvol - q(num2)*dxm(num2)
                  endif
                  !
               else
                  !
                  if (ndm1>0) then
                     dvol = dvol + q(ndm1)*dxrm(uv_flags_iref(ndm1))
                  endif   
                  !
                  if (ndm2>0) then
                    dvol = dvol + q(ndm2)*dxrm(uv_flags_iref(ndm2))
                  endif 
                  !
                  if (num1>0) then
                     dvol = dvol - q(num1)*dxrm(uv_flags_iref(num1))
                  endif   
                  !
                  if (num2>0) then
                     dvol = dvol - q(num2)*dxrm(uv_flags_iref(num2))
                  endif
                  !
               endif               
               !
            else
               !
               nmd = z_index_uv_md1(nm)
               nmu = z_index_uv_mu1(nm)
               ndm = z_index_uv_nd1(nm)
               num = z_index_uv_nu1(nm)
               !
               dvol = 0.0
               !
               ! U
               !
               if (nmd>0) then
                  dvol = dvol + q(nmd)*dyrm(1)
               endif   
               !
               if (nmu>0) then
                  dvol = dvol - q(nmu)*dyrm(1)
               endif   
               !
               ! V
               !
               if (crsgeo) then
                  !
                  if (ndm>0) then
                     dvol = dvol + q(ndm)*dxm(ndm)
                  endif   
                  !
                  if (num>0) then
                     dvol = dvol - q(num)*dxm(num)
                  endif   
                  !
               else
                  !
                  if (ndm>0) then
                     dvol = dvol + q(ndm)*dxrm(1)
                  endif   
                  !
                  if (num>0) then
                     dvol = dvol - q(num)*dxrm(1)
                  endif   
                  !
               endif   
               !
            endif   
            !
            if (crsgeo) then
               zs(nm) = zs(nm) + dvol*dt/cell_area_m2(nm)
            else                  
               zs(nm) = zs(nm) + dvol*dt/cell_area(z_flags_iref(nm))
            endif
            !
         endif
         !
         if (store_maximum_waterlevel) then
            zsmax(nm) = max(zsmax(nm), zs(nm))
         endif
         !
         ! No continuity update but keeping track of variables
         ! to do is vmax
         ! 1. store Twet
         if (store_twet) then
            if ( (zs(nm) - zb(nm)) > twet_threshold) then
                twet(nm) = twet(nm) + dt
            endif
          endif
         !
      endif
      !
      if (wavemaker) then
         !      
         if (kcs(nm)==4) then
            !
            ! Wave maker point
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
               qnmd = q(z_index_uv_md1(nm))
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
               qnmu = q(z_index_uv_mu1(nm))
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
               qndm = q(z_index_uv_nd1(nm))
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
               qnum = q(z_index_uv_nu1(nm))
               !
            endif   
            !
            zs(nm)   = zs(nm) + (((qnmd - qnmu)*dxrinv(z_flags_iref(nm)) + (qndm - qnum)*dyrinv(z_flags_iref(nm))))*dt
            !
         endif   
         !
         ! Would double exponential filtering be better?
         !
         zsm(nm) = factime*zs(nm) + (1.0 - factime)*zsm(nm)
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
   !
   real*4           :: qxnm
   real*4           :: qxnmd
   real*4           :: qynm
   real*4           :: qyndm
   real*4           :: factime
   real*4           :: dvol
   !
   integer          :: iuv
   real*4           :: dzvol
   real*4           :: facint
   real*4           :: a
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
   !$omp private ( dvol,nmd1,nmu1,ndm1,num1,nmd2,nmu2,ndm2,num2,nmd,nmu,ndm,num,a,iuv,facint,dzvol,ind )
   !$omp do schedule ( dynamic, 256 )
   !$acc kernels present( kcs, zs, zb, z_volume, zsmax, twet, zsm, &
   !$acc                  subgrid_z_zmin,  subgrid_z_zmax, subgrid_z_dep, subgrid_z_volmax, &
   !$acc                  netprcp, cumprcpt, prcp, q, z_flags_type, z_flags_iref, uv_flags_iref, &
   !$acc                  z_index_uv_md1, z_index_uv_md2, z_index_uv_nd1, z_index_uv_nd2, z_index_uv_mu1, z_index_uv_mu2, z_index_uv_nu1, z_index_uv_nu2, &
   !$acc                  dxm, dxrm, dyrm, dxminv, dxrinv, dyrinv, cell_area_m2, cell_area, &
   !$acc                  z_index_wavemaker, wavemaker_uvmean, wavemaker_nmd, wavemaker_nmu, wavemaker_ndm, wavemaker_num), async(1)
   !$acc loop independent, private( nm )
   do nm = 1, np
      !
      ! And now water level changes due to horizontal fluxes
      !
      if (kcs(nm)==1 .or. kcs(nm)==4) then
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
               z_volume(nm) = z_volume(nm) + cumprcpt(nm)*a
               cumprcpt(nm) = 0.0
               !
            endif
            !
         endif
         !
         if (z_flags_type(nm) == 0) then
            !
            ! Regular point with four surrounding cells of the same size
            !
            nmd = z_index_uv_md1(nm)
            nmu = z_index_uv_mu1(nm)
            ndm = z_index_uv_nd1(nm)
            num = z_index_uv_nu1(nm)
            !
            if (crsgeo) then
               !
               z_volume(nm) = z_volume(nm) + ( (q(nmd) - q(nmu))*dyrm(uv_flags_iref(nm)) + (q(ndm) - q(num))*dxm(nm) ) * dt
               !
            else
               !
               if (use_quadtree) then   
                  z_volume(nm) = z_volume(nm) + ( (q(nmd) - q(nmu))*dyrm(uv_flags_iref(nm)) + (q(ndm) - q(num))*dxrm(uv_flags_iref(nm)) ) * dt
               else
                  z_volume(nm) = z_volume(nm) + ( (q(nmd) - q(nmu))*dy + (q(ndm) - q(num))*dx ) * dt
               endif   
               !
            endif   
            !
         else
            !
            if (use_quadtree) then   
               !
               ! More complicated
               !
               nmd1 = z_index_uv_md1(nm)
               nmu1 = z_index_uv_mu1(nm)
               ndm1 = z_index_uv_nd1(nm)
               num1 = z_index_uv_nu1(nm)
               nmd2 = z_index_uv_md2(nm)
               nmu2 = z_index_uv_mu2(nm)
               ndm2 = z_index_uv_nd2(nm)
               num2 = z_index_uv_nu2(nm)
               !
               dvol = 0.0
               !
               ! U
               !
               if (nmd1>0) then
                  dvol = dvol + q(nmd1)*dyrm(uv_flags_iref(nmd1))
               endif   
               !
               if (nmd2>0) then
                  dvol = dvol + q(nmd2)*dyrm(uv_flags_iref(nmd2))
               endif
               !
               if (nmu1>0) then
                  dvol = dvol - q(nmu1)*dyrm(uv_flags_iref(nmu1))
               endif   
               !
               if (nmu2>0) then
                  dvol = dvol - q(nmu2)*dyrm(uv_flags_iref(nmu2))
               endif
               !
               ! V
               !
               if (crsgeo) then
                  !
                  if (ndm1>0) then                     
                     dvol = dvol + q(ndm1)*dxm(ndm1)
                  endif   
                  !
                  if (ndm2>0) then
                     dvol = dvol + q(ndm2)*dxm(ndm2)
                  endif 
                  !
                  if (num1>0) then
                     dvol = dvol - q(num1)*dxm(num1)
                  endif   
                  !
                  if (num2>0) then
                     dvol = dvol - q(num2)*dxm(num2)
                  endif
                  !
               else
                  !
                  if (ndm1>0) then
                     dvol = dvol + q(ndm1)*dxrm(uv_flags_iref(ndm1))
                  endif   
                  !
                  if (ndm2>0) then
                    dvol = dvol + q(ndm2)*dxrm(uv_flags_iref(ndm2))
                  endif 
                  !
                  if (num1>0) then
                     dvol = dvol - q(num1)*dxrm(uv_flags_iref(num1))
                  endif   
                  !
                  if (num2>0) then
                     dvol = dvol - q(num2)*dxrm(uv_flags_iref(num2))
                  endif
                  !
               endif               
               !
            else
               !
               nmd = z_index_uv_md1(nm)
               nmu = z_index_uv_mu1(nm)
               ndm = z_index_uv_nd1(nm)
               num = z_index_uv_nu1(nm)
               !
               dvol = 0.0
               !
               !
               ! U
               !
               if (nmd>0) then
                  dvol = dvol + q(nmd)*dyrm(1)
               endif   
               !
               if (nmu>0) then
                  dvol = dvol - q(nmu)*dyrm(1)
               endif   
               !
               ! V
               !
               if (crsgeo) then
                  !
                  if (ndm>0) then
                     dvol = dvol + q(ndm)*dxm(ndm)
                  endif   
                  !
                  if (num>0) then
                     dvol = dvol - q(num)*dxm(num)
                  endif   
                  !
               else
                  !
                  if (ndm>0) then
                     dvol = dvol + q(ndm)*dxrm(1)
                  endif   
                  !
                  if (num>0) then
                     dvol = dvol - q(num)*dxrm(1)
                  endif   
                  !
               endif   
               !
            endif   
            !
            z_volume(nm) = z_volume(nm) + dvol*dt
            !
         endif
         !
         if (z_volume(nm)>=subgrid_z_volmax(nm) - 1.0e-6) then
            !
            ! Entire cell is wet, no interpolation needed
            !
            zs(nm) = max(subgrid_z_zmax(nm), -20.0) + (z_volume(nm) - subgrid_z_volmax(nm))/a
            !
         elseif (z_volume(nm)<=1.0e-6) then
            !
            zs(nm) = max(subgrid_z_zmin(nm), -20.0)
            !
         else
            !
            dzvol    = subgrid_z_volmax(nm) / (subgrid_nbins - 1)
            iuv      = min(int(z_volume(nm)/dzvol) + 1, subgrid_nbins - 1)
            facint   = (z_volume(nm) - (iuv - 1)*dzvol ) / dzvol
            zs(nm)   = subgrid_z_dep(iuv, nm) + (subgrid_z_dep(iuv + 1, nm) - subgrid_z_dep(iuv, nm))*facint                    
            !
         endif
         !
         if (store_maximum_waterlevel) then
            !
            zsmax(nm) = max(zsmax(nm), zs(nm))
            !
         endif
         !
         !
      endif
      !       
      ! No continuity update but keeping track of variables
      ! to do is vmax
      ! 1. store Twet
      if (store_twet) then
        if ( (zs(nm) - subgrid_z_zmin(nm)) > twet_threshold) then
            twet(nm) = twet(nm) + dt
        endif
      endif
      !
      if (wavemaker) then
         !
         zsm(nm) = factime*zs(nm) + (1.0 - factime)*zsm(nm)
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
   
   
end module
