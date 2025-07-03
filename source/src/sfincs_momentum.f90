module sfincs_momentum
   !
   use sfincs_data
   !
   implicit none
   !     
contains
   !
   subroutine compute_fluxes(dt, tloop)
   !
   ! Computes fluxes over subgrid u and v points
   !
   integer   :: count0
   integer   :: count1
   integer   :: count_rate
   integer   :: count_max
   real      :: tloop
   !
   real*4    :: dt
   !
   integer   :: ip
   integer   :: nm
   integer   :: nmu
   integer   :: n
   integer   :: m

   integer   :: idir
   integer   :: iref
   integer   :: itype
   integer   :: iuv
   integer   :: ind
   integer   :: icuv
   integer   :: iuv1
   integer   :: iuv2
   !
   real*4    :: hu
   real*4    :: dxuvinv
   real*4    :: dxuv2inv
   real*4    :: dyuvinv
   real*4    :: dyuv2inv
   real*4    :: adv
   real*4    :: fcoriouv
   real*4    :: frc
   real*4    :: gammax
   real*4    :: facmax
   real*4    :: wsumax
   !
   real*4    :: qx_nm
   real*4    :: qx_nmd
   real*4    :: qx_nmu
   real*4    :: qy_nm
   real*4    :: qy_nmu
   real*4    :: qy_ndm
   real*4    :: qy_ndmu
   real*4    :: qfr
   real*4    :: qsm
   !   
   real*4    :: uu_nm
   real*4    :: uu_nmd
   real*4    :: uu_nmu
   real*4    :: uu_ndm
   real*4    :: uu_num
   real*4    :: vu
   real*4    :: qvt
   !
   real*4    :: zsu
   real*4    :: dzuv
   real*4    :: facint
   real*4    :: gnavg2
   real*4    :: fwmax
   real*4    :: zmax
   real*4    :: zmin
   real*4    :: one_minus_facint 
   !
   real*4    :: dqxudx
   real*4    :: dqyudy
   real*4    :: qu
   real*4    :: qd
   real*4    :: uu
   real*4    :: ud
   real*4    :: qy
   real*4    :: dzdx
   !
   real*4    :: hwet
   real*4    :: phi
   !
   real*4    :: mdrv
   real*4    :: hu73
   !
   real*4, parameter :: expo = 1.0 / 3.0
   !integer, parameter :: expo = 1
   !
   logical   :: iok
   !
   call system_clock(count0, count_rate, count_max)
   !
   min_dt = dtmax
   !
   !$acc update device(min_dt)
   !
   ! Copy flux and velocity from previous time step
   !
   ! For some reason, it is necessary to set num_gangs here! Without, the program launches only 1 gang, and everything becomes VERY slow!
   !
   !$acc parallel, present( kcuv, kfuv, zs, q, q0, uv, uv0, zsderv, &
   !$acc                    uv_flags_iref, uv_flags_type, uv_flags_dir, mask_adv, &
   !$acc                    subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_havg, subgrid_uv_nrep, subgrid_uv_pwet, &
   !$acc                    subgrid_uv_havg_zmax, subgrid_uv_nrep_zmax, subgrid_uv_fnfit, subgrid_uv_navg_w, &
   !$acc                    uv_index_z_nm, uv_index_z_nmu, uv_index_u_nmd, uv_index_u_nmu, uv_index_u_ndm, uv_index_u_num, &
   !$acc                    uv_index_v_ndm, uv_index_v_ndmu, uv_index_v_nm, uv_index_v_nmu, cuv_index_uv, cuv_index_uv1, cuv_index_uv2, &
   !$acc                    zb, zbuv, zbuvmx, tauwu, tauwv, patm, fwuv, gn2uv, dxminv, dxrinv, dyrinv, dxm2inv, dxr2inv, dyr2inv, &
!   !$acc                    dxrinvc, fcorio2d, nuvisc, min_dt, z_volume, gnapp2 )
   !$acc                    dxrinvc, fcorio2d, nuvisc, min_dt, z_volume, gnapp2 ) num_gangs( 1024 ) vector_length( 128 )
   !
   !$omp parallel &
   !$omp private ( ip )
   !$omp do
   !$acc loop independent gang vector
   do ip = 1, npuv + ncuv
      !
      q0(ip)  = q(ip)
      uv0(ip) = uv(ip)
      !
   enddo
   !$omp end do
   !$omp end parallel
   !
   ! Update fluxes
   !
   !$omp parallel &
   !$omp private ( ip,hu,qfr,qsm,qx_nm,nm,nmu,dzdx,frc,idir,itype,iref,dxuvinv,dxuv2inv,dyuvinv,dyuv2inv, &
   !$omp           qx_nmd,qx_nmu,qy_nm,qy_ndm,qy_nmu,qy_ndmu,uu_nm,uu_nmd,uu_nmu,uu_num,uu_ndm,vu, & 
   !$omp           fcoriouv,gnavg2,iok,zsu,dzuv,iuv,facint,fwmax,zmax,zmin,one_minus_facint,dqxudx,dqyudy,uu,ud,qu,qd,qy,hwet,phi,adv,mdrv,hu73 ) &
   !$omp reduction ( min : min_dt  )
   !$omp do schedule ( dynamic, 256 )
   !$acc loop reduction( min : min_dt ) independent gang vector
   do ip = 1, npuv
      !
      if (kcuv(ip) == 1 .or. kcuv(ip) == 6) then
         !
         ! Regular UV point (or a coastal lateral boundary point)
         !
         ! Indices of surrounding water level points
         !
         nm  = uv_index_z_nm(ip)
         nmu = uv_index_z_nmu(ip)
         !
         iok  = .false.
         !
         zsu = max(zs(nm), zs(nmu)) ! water level at u point
         !
         if (subgrid) then
            !
            zmin = subgrid_uv_zmin(ip)
            zmax = subgrid_uv_zmax(ip)
            !
            if (zsu>zmin + huthresh) then ! In the 'new' subgrid formulations, zmin is lowest pixel + huthresh. Huthresh was already applied when building the subgrid tables. In sfincs_domain, huthresh is set to 0.0 to acount for this.
               iok = .true.
            endif   
            !
         else
            !            
            if (zsu>zbuvmx(ip)) then ! zbuvmx = max(zb(nm), zb(nmu)) + huthresh
               iok = .true.
            endif   
            !            
         endif
         !
         if (iok) then
            !
            ! UV point is wet 
            !
            if (use_quadtree) then
               iref  = uv_flags_iref(ip) ! refinement level
               itype = uv_flags_type(ip) ! -1 is fine to coarse, 0 is normal, 1 is coarse to fine
            else
               iref  = 1
               itype = 0
            endif
            !
            idir  = uv_flags_dir(ip) ! 0 is u, 1 is v
            !
            ! Determine grid spacing (and coriolis factor fcoriouv)
            !
            if (crsgeo) then
               !
               ! Geographic coordinate system
               !
               if (itype==0) then
                  !
                  ! Regular
                  !
                  if (idir==0) then
                     !
                     ! U point
                     !
                     dxuvinv  = dxminv(ip)
                     dyuvinv  = dyrinv(iref)
                     dxuv2inv = dxm2inv(ip)
                     dyuv2inv = dyr2inv(iref)
                     !
                  else
                     !
                     ! V point
                     !
                     dxuvinv  = dyrinv(iref)
                     dyuvinv  = dxminv(ip)
                     dxuv2inv = dyr2inv(iref)
                     dyuv2inv = dxm2inv(ip)
                     !
                  endif
                  !
               else   
                  !
                  ! Fine to coarse or coarse to fine
                  !
                  if (idir==0) then
                     !
                     dxuvinv = 1.0 / (3*(1.0/dxminv(ip))/2)
                     dyuvinv = dyrinv(iref)
                     dxuv2inv = 0.0 ! no viscosity term
                     dyuv2inv = dyr2inv(iref)
                     !
                  else   
                     !
                     dxuvinv = 1.0 / (3*(1.0/dyrinv(iref))/2)
                     dyuvinv  = dxminv(ip)
                     dxuv2inv = 0.0 ! no viscosity term
                     dyuv2inv = dxm2inv(ip)
                     !
                  endif
                  !
               endif
               !
               fcoriouv = fcorio2d(nm)
               !
            else
               !
               ! Projected coordinate system
               !
               if (itype==0) then
                  !
                  ! Regular
                  !
                  if (idir==0) then
                     !
                     ! U point
                     !
                     dxuvinv  = dxrinv(iref)
                     dyuvinv  = dyrinv(iref)
                     dxuv2inv = dxr2inv(iref)
                     dyuv2inv = dyr2inv(iref)
                     !
                  else
                     !
                     ! V point
                     !
                     dxuvinv  = dyrinv(iref)
                     dyuvinv  = dxrinv(iref)
                     dxuv2inv = dyr2inv(iref)
                     dyuv2inv = dxr2inv(iref)
                     !
                  endif   
                  !
               else   
                  !
                  ! Fine to coarse or coarse to fine
                  !
                  if (idir==0) then
                     !
                     ! U point
                     !
                     dxuvinv  = dxrinvc(iref)
                     dyuvinv  = dyrinv(iref)
                     dxuv2inv = 0.0 ! no viscosity
                     dyuv2inv = dyr2inv(iref)
                     !
                  else
                     !
                     ! V point
                     !
                     dxuvinv  = dyrinv(iref)
                     dyuvinv  = dxrinv(iref)
                     dxuv2inv = 0.0
                     dyuv2inv = dxr2inv(iref)
                     !
                  endif   
                  !
               endif
               !
               fcoriouv = fcorio
               !
            endif
            !
            ! Get fluxes and velocities from previous time step
            !
            qx_nm = q0(ip)
            !
            if (advection .or. coriolis .or. viscosity .or. friction2d .or. thetasmoothing) then
               !
               ! Get the neighbors
               !
               qx_nmd  = q0(uv_index_u_nmd(ip))
               qx_nmu  = q0(uv_index_u_nmu(ip))
               qy_nm   = q0(uv_index_v_nm(ip))
               qy_nmu  = q0(uv_index_v_nmu(ip))
               qy_ndm  = q0(uv_index_v_ndm(ip))
               qy_ndmu = q0(uv_index_v_ndmu(ip))
               !
               uu_nm   = uv0(ip)
               uu_nmd  = uv0(uv_index_u_nmd(ip))
               uu_nmu  = uv0(uv_index_u_nmu(ip))
               uu_ndm  = uv0(uv_index_u_ndm(ip))
               uu_num  = uv0(uv_index_u_num(ip))
               vu      = (uv0(uv_index_v_ndm(ip)) + uv0(uv_index_v_ndmu(ip)) + uv0(uv_index_v_nm(ip)) + uv0(uv_index_v_nmu(ip))) / 4
               !
            endif
            !
            ! Wet fraction phi (for non-subgrid or original subgrid approach phi should be 1.0)
            !
            phi  = 1.0
            !
            ! Compute water depth at uv point
            !
            if (subgrid) then
               !
               if (zsu > zmax) then
                  !
                  ! Entire cell is wet, no interpolation from table needed for depth hu
                  !
                  hu = subgrid_uv_havg_zmax(ip) + zsu
                  !
                  if (wave_enhanced_roughness) then
                     !
                     ! Apparent roughness is computed in sfincs_wave_enhanced_roughness.f90. It is called by sfincs_bmi.f90.
                     ! Note: wave enhanced roughness is only done for uv points that are completely wet!
                     !
                     gnavg2 = gnapp2(ip)
                     !
                  else
                     ! 
                     ! Use fitting function for gnavg2 
                     !
                     gnavg2 = subgrid_uv_navg_w(ip) - (subgrid_uv_navg_w(ip) - subgrid_uv_nrep_zmax(ip)) / (subgrid_uv_fnfit(ip) * (zsu - zmax) + 1.0)
                     ! 
                  endif
                  !
               else
                  !
                  ! Interpolation required
                  !
                  dzuv   = (zmax - zmin) / (subgrid_nlevels - 1)                                                           ! level size (is storing this in memory faster?)
                  iuv    = min(int((zsu - zmin) / dzuv) + 1, subgrid_nlevels - 1)                                          ! index of level below zsu 
                  facint = (zsu - (zmin + (iuv - 1) * dzuv) ) / dzuv                                                       ! 1d interpolation coefficient
                  !
                  hu     = subgrid_uv_havg(iuv, ip) + (subgrid_uv_havg(iuv + 1, ip) - subgrid_uv_havg(iuv, ip)) * facint   ! grid-average depth
                  gnavg2 = subgrid_uv_nrep(iuv, ip) + (subgrid_uv_nrep(iuv + 1, ip) - subgrid_uv_nrep(iuv, ip)) * facint   ! representative g*n^2
                  phi    = subgrid_uv_pwet(iuv, ip) + (subgrid_uv_pwet(iuv + 1, ip) - subgrid_uv_pwet(iuv, ip)) * facint   ! wet fraction
                  !
               endif
               !
            else
               !
               hu     = max(zsu - zbuvmx(ip), huthresh)
               gnavg2 = gn2uv(ip)
               !
            endif
            !
            ! If coupled with a wave model (SnapWave or HurryWave), we may want to set an apparent roughness (which is computed elsewhere for each grid cell)
            !
            ! if (use_apparent_roughness) then
            !    !
            !    gnavg2 = 9.81 * (0.5 * (napp(nm) + napp(nmu))**2
            !    !
            ! endif
            !
            ! Compute wet average depth hwet (used in wind and wave forcing)
            !
            hwet = hu / phi
            !
            ! FORCING TERMS
            !
            ! Pressure term 
            !
            ! Apply slope limiter to dzdx (turned off by default)
            !
            if (slopelim < 9999.0) then
               !
               dzdx = min(max((zs(nmu) - zs(nm)) * dxuvinv, -slopelim), slopelim) 
               !
            else
               !
               dzdx = (zs(nmu) - zs(nm)) * dxuvinv
               !
            endif
            !
            frc = - g * hu * dzdx
            !
            ! Advection term
            !
            if (advection) then
               !
               ! Turn off advection next to open boundaries
               !
               if (mask_adv(ip) == 1) then
                  ! 
                  dqxudx = 0.0  
                  dqyudy = 0.0  
                  !
                  if (advection_scheme == 0) then
                     !
                     ! Original (not very good, but reasonably robust)
                     !
                     if (qx_nm > 1.0e-6) then
                        !
                        dqxudx = (qx_nm*uu_nm - qx_nmd*uu_nmd) * dxuvinv                  
                        ! 
                     elseif (qx_nm < -1.0e-6) then
                        !
                        dqxudx = (qx_nmu*uu_nmu - qx_nm*uu_nm) * dxuvinv
                        !
                     endif
                     !
                     qy = qy_ndm + qy_ndmu + qy_nm + qy_nmu
                     !
                     if ( qy > 1.0e-6 ) then
                        !
                        dqyudy = ( qy_ndm + qy_ndmu ) * (uu_nm - uu_ndm) * dyuvinv / 2
                        ! 
                     elseif ( qy < -1.0e-6 ) then
                        !
                        dqyudy = (qy_nm + qy_nmu) * (uu_num - uu_nm) * dyuvinv / 2
                        !
                     endif
                     !
                  elseif (advection_scheme == 1) then
                     !
                     ! 1st order upwind (upw1)
                     !
                     ! d qu u / dx = qu du / dx + u d qu / dx
                     ! d qv u / dy = qv du / dy + u d qv / dy
                     !
                     ! d qu u / dx
                     !
                     qd = (qx_nmd + qx_nm) / 2
                     qu = (qx_nm + qx_nmu) / 2
                     !
                     if (qd > 1.0e-6) then
                        ud = (uu_nmd + uu_nm) / 2
                        dqxudx = ( qd * (uu_nm - uu_nmd) + ud * (qx_nm - qx_nmd) ) * dxuvinv
                     endif
                     !
                     if (qu < -1.0e-6) then
                        uu = (uu_nm + uu_nmu) / 2
                        dqxudx = dqxudx + ( qu*(uu_nmu - uu_nm) + uu * (qx_nmu - qx_nm) ) * dxuvinv
                     endif
                     ! 
                     ! d qv u / dy
                     !
                     ! qv d u / d y
                     ! 
                     qu = (qy_nm + qy_nmu) / 2
                     qd = (qy_ndm + qy_ndmu) / 2
                     !
                     if (qd > 1.0e-6) then
                        dqyudy = qd * (uu_nm - uu_ndm) * dyuvinv
                     endif
                     !
                     if (qu < -1.0e-6) then
                        dqyudy = dqyudy + qu * (uu_num - uu_nm) * dyuvinv
                     endif
                     ! 
                     ! u d qv / dy
                     ! 
                     ud = (uu_nmd + uu_nm) / 2
                     uu = (uu_nm + uu_nmu) / 2
                     !
                     if (ud > 1.0e-6) then
                        dqyudy = dqyudy + ud * ( qy_nm - qy_ndm ) * dyuvinv
                     endif   
                     !
                     if (uu < -1.0e-6) then
                        dqyudy = dqyudy + uu * ( qy_nmu - qy_ndmu ) * dyuvinv
                     endif
                     !  
                  endif
                  !
                  adv = - phi * (dqxudx + dqyudy)
                  !
                  ! Limit advection term such that horizontal acceleration due to advection does not exceed advlim (default 1.0 m/s2)
                  ! Default advlim is 1.0 m/s2
                  !
                  adv = min(max(adv, - advlim * hu), advlim * hu) 
                  !
                  frc = frc + adv
                  !
               endif
               !   
            endif   
            !
            ! Viscosity term
            !
            if (viscosity) then
               !
               if (itype == 0) then
                  ! 
                  frc = frc + nuvisc(iref) * hu * ( (uu_nmu - 2*uu_nm + uu_nmd ) * dxuv2inv + (uu_num - 2*uu_nm + uu_ndm ) * dyuv2inv )
                  !
               else
                  !
                  ! Increase viscosity to prevent instabilities on refinement (related to advection?)
                  !
                  frc = frc + nuviscfac * nuvisc(iref) * hu * ( (uu_nmu - 2*uu_nm + uu_nmd ) * dxuv2inv + (uu_num - 2*uu_nm + uu_ndm ) * dyuv2inv )
                  !
               endif               
               !
            endif
            !
            ! Coriolis term
            !
            if (coriolis) then
               !
               if (idir==0) then
                  !
                  frc = frc + fcoriouv * hu * vu ! U
                  !
               else
                  !
                  frc = frc - fcoriouv * hu * vu ! V
                  !
               endif
               !
            endif            
            !
            ! Wind forcing
            !
            if (wind) then
               !
               if (hwet > 0.25) then
                  !
                  if (idir==0) then
                     !
                     frc = frc + phi * tauwu(nm)
                     !
                  else
                     !
                     frc = frc + phi * tauwv(nm)
                     !
                  endif   
                  !
               else
                  !
                  ! Reduce wind drag at water depths < 0.25 m
                  !
                  if (idir==0) then
                     !
                     ! frc = frc + phi * tauwu(nm) * hwet * 4
                     frc = frc + tauwu(nm) * hu * 4
                     !
                  else
                     !
                     ! frc = frc + phi * tauwv(nm) * hwet * 4
                     frc = frc + tauwv(nm) * hu * 4
                     !
                  endif   
                  !
               endif
               !
            endif            
            !
            ! Atmospheric pressure
            !
            if (patmos) then
               !
               frc = frc + hu * (patm(nm) - patm(nmu)) * dxuvinv / rhow
               !
            endif
            !
            ! Wave forcing
            !
            if (snapwave) then
               !
               ! Limited wave forces in shallow water 
               !
               ! facmax = 0.25*sqrt(g)*rhow*gammax**2
               ! fmax = facmax*hu*sqrt(hu)/tp/rhow (we already divided by rhow in sfincs_snapwave)
               !
               fwmax = 0.8 * hwet * sqrt(hwet) / 6
               !
               frc = frc + phi * sign(min(abs(fwuv(ip)), fwmax), fwuv(ip))
               !
            endif
            !
            ! Compute flux qfr used for friction term
            !
            if (kfuv(ip) == 0) then
               !
               ! This uv point just became wet, so estimate equilibrium flux 
               !
               qfr = sqrt(abs(dzdx) / (max(gnavg2, 1.0e-5) / 10)) * hu ** (5.0 / 3.0)
               !
            else   
               !
               if (friction2d) then
                  !
                  ! Computed friction term with both qx and qy
                  !
                  qfr = sqrt(qx_nm**2 + (hu * vu)**2)
                  !
               else
                  !
                  ! Computed friction term with only qx (original Bates et al. (2010))
                  !
                  qfr = abs(qx_nm) ! flux to be used in the friction term
                  !
               endif
               !
            endif
            !
            ! Apply some smoothing if theta < 1.0 (not recommended anymore!)
            !
            qsm = qx_nm
            !
            if (thetasmoothing) then
               ! 
               ! Apply theta smoothing 
               ! 
               if ( kcuv(uv_index_u_nmd(ip))==1 .and. kcuv(uv_index_u_nmu(ip))==1 ) then 
                  !
                  ! But only at regular points
                  ! 
                  if (abs(qx_nmd) > 1.0e-6 .and. abs(qx_nmu) > 1.0e-6) then
                     ! if (kfu(uv_index_u_nmd(ip))==1 .and. kfu(uv_index_u_nmu(ip))==1) then
                     !
                     ! And if both uv neighbors are active
                     ! 
                     qsm = theta*qx_nm + 0.5 * (1.0 - theta) * (qx_nmu + qx_nmd)
                     ! 
                  endif
                  ! 
               endif
               !
            endif            
            !
            ! Compute new flux for this uv point (Bates et al., 2010)
            ! 
            if (h73table) then
               !
               ! Get hu**(7/3) from look-up table
               !
               hu73 = power7over3(max(hu, 1.0e-6))
               !
            else
               !
               ! Compute hu**(7/3)
               !
               hu73 = hu**2 * hu**expo
               !
            endif
            ! 
            q(ip) = (qsm + frc * dt) / (1.0 + gnavg2 * dt * qfr / hu73)
            !
            if (subgrid .and. wiggle_suppression) then 
               !
               ! If the acceleration of water level in cell nm is large and positive and in nmu large and negative, or vice versa, apply limiter to the flux. Only for subgrid.
               !
               mdrv = abs(zsderv(nm) - zsderv(nmu)) - wiggle_threshold
               !
               if (mdrv > 0.0) then
                  !
                  q(ip) = q(ip) * wiggle_threshold / (wiggle_factor * mdrv + wiggle_threshold)
                  !
               endif
               !
            endif
            !
            ! Making sure that no water can flow out of a cell when its water depth is negative
            !            
            if (subgrid) then
               !
               if (z_volume(nm) < 0.0) then
                  q(ip) = min(q(ip), 0.0)
               endif
               !
               if (z_volume(nmu) < 0.0) then
                  q(ip) = max(q(ip), 0.0)
               endif
               !
            else
               !
               if (zs(nm) < zb(nm)) then
                  q(ip) = min(q(ip), 0.0)
               endif
               !
               if (zs(nmu) < zb(nmu)) then
                  q(ip) = max(q(ip), 0.0)
               endif
               !
            endif
            !
            ! Apply flux limiter (default 10 m/s)
            !
            q(ip) = min(max(q(ip), - hu * uvlim), hu * uvlim)
            !
            ! Compute velocity
            !
            uv(ip) = q(ip) / hu ! Store wet-averaged velocity
            !
            kfuv(ip) = 1
            !
            ! Determine minimum time step (alpha is added later on in sfincs_lib.f90) of all uv points
            ! Use maximum of sqrt(gh) and current velocity
            !
            min_dt = min(min_dt, 1.0 / ( max(sqrt(g * hu), abs(uv(ip)) ) * dxuvinv))
            !
         else
            !
            q(ip)  = 0.0
            uv(ip) = 0.0
            kfuv(ip) = 0
            !
         endif
         !
      endif
   enddo   
   !$omp end do
   !$omp end parallel
   !
   if (ncuv > 0) then
      !
      ! Loop through combined uv points and determine average uv and q
      ! The combined q and uv values are used in the continuity equation and in the netcdf output
      !
      !$omp parallel &
      !$omp private ( icuv )
      !$omp do
      !$acc loop independent gang vector
      do icuv = 1, ncuv
         !
         ! Average of the two uv points
         !
         q(cuv_index_uv(icuv))  = (q(cuv_index_uv1(icuv)) + q(cuv_index_uv2(icuv))) / 2
         uv(cuv_index_uv(icuv)) = (uv(cuv_index_uv1(icuv)) + uv(cuv_index_uv2(icuv))) / 2
         !
      enddo
      !$omp end do
      !$omp end parallel
      !
   endif
   !
   !$acc end parallel
   !
   !$acc update host(min_dt)
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine      
   !
   !
   function power7over3(hu) result(hu73)
   !
   ! Computes hu^(7/3) using a table for hu < 10^4 and hu^(7/3) for hu >= 10^4
   !
   real*4, intent(in) :: hu
   real*4             :: hu73
   !
   !!$acc routine(power7over3) seq
   !
   if (hu < 0.00001) then
      hu73 = x73(int(1e8 * hu), 1)
   elseif (hu < 0.0001) then
      hu73 = x73(int(1e7 * hu), 2)
   elseif (hu < 0.001) then
      hu73 = x73(int(1e6 * hu), 3)
   elseif (hu < 0.01) then
      hu73 = x73(int(1e5 * hu), 4)
   elseif (hu < 0.1) then
      hu73 = x73(int(1e4 * hu), 5)
   elseif (hu < 1.0) then
      hu73 = x73(int(1e3 * hu), 6)
   elseif (hu < 10.0) then
      hu73 = x73(int(100 * hu), 7)
   elseif (hu < 100.0) then
      hu73 = x73(int(10 * hu), 8)
   elseif (hu < 1000.0) then
      hu73 = x73(int(1 * hu), 9)
   elseif (hu < 10000.0) then
      hu73 = x73(int(0.1 * hu), 10)
   else
      hu73 = 10000.0 ** (7.0/3.0)
   endif    
   !
   end function
   !
end module
