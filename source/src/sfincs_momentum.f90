module sfincs_momentum
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! Momentum equation: computes the fluxes and velocities at the uv points every time step.
   !
   ! Subroutines / functions:
   !
   !   compute_fluxes(dt)
   !     Compute the fluxes q and velocities uv at all uv points (regular and subgrid).
   !     Called from: sfincs_update (sfincs_lib)
   !
   ! ------------ Copyright notice SFINCS ------------
   !
   ! Copyright (C)  Stichting Deltares, 2016-2026.
   !
   ! This file is part of SFINCS.
   !
   ! SFINCS is free software: you can redistribute it and/or modify
   ! it under the terms of the GNU General Public License as
   ! published by the Free Software Foundation version 3.
   !
   ! SFINCS is distributed in the hope that it will be useful,
   ! but WITHOUT ANY WARRANTY; without even the implied warranty of
   ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   ! GNU General Public License for more details.
   !
   ! You should have received a copy of the GNU General Public License
   ! along with SFINCS.  If not, see <http://www.gnu.org/licenses/>.
   !
   ! contact: software@deltares.nl
   ! Stichting Deltares
   ! P.O. Box 177
   ! 2600 MH Delft, The Netherlands
   !
   ! All indications, logos of, and references to, "SFINCS"
   ! remain the property of Stichting Deltares. All rights reserved.
   !
   ! ------------ Copyright notice SFINCS ------------
   !
   use sfincs_data
   use sfincs_timers
   !
   implicit none
   !
contains
   !
   subroutine compute_fluxes(dt)
   !
   ! Computes fluxes over subgrid u and v points
   !
   ! Input:
   !   dt : time step [s]
   !
   ! Called from: sfincs_update (sfincs_lib)
   !
   real*4    :: dt
   !
   integer   :: ip
   integer   :: nm
   integer   :: nmu
   integer   :: n
   integer   :: m
   !
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
   integer   :: iup                             ! upwind cell of the uv point
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
   real*4    :: q_abs                           ! |q0| at this uv point (floored), for the partner weights
   real*4    :: q_k                             ! |q0| at a partner face
   real*4    :: q_ref                           ! flux below which a partner gets a weight below 1
   real*4    :: wt_nmd, wt_nmu, wt_ndm, wt_num  ! weights of the four partner faces in the advection stencil
   !
   real*4    :: hwet
   real*4    :: phi
   !
   real*4    :: mdrv                            ! wiggle suppression driver
   real*4    :: phiz                            ! wet-area fraction of the least wet neighbouring cell
   real*4    :: fac                             ! under-relaxation factor for wiggle suppression
   real*4    :: hu73
   real*4    :: y_cbrt                          ! cube-root approximation hu^(1/3)
   integer*4 :: i_cbrt                          ! bit pattern of hu for the cube-root seed
   !
   real*4    :: min_dt_ip
   !
   real*4, dimension(:), allocatable :: swap                 ! work array to swap q and q0, uv and uv0
   !
   logical   :: iok
   !
   call timer_start('momentum')
   !
   min_dt = dtmax
   !
   if (timestep_analysis) then
       !
       ! Do in loop for updating on GPU
       !
       !$acc parallel, present( timestep_analysis_required_timestep )
       !$omp parallel &
       !$omp private ( ip )
       !$omp do
       !$acc loop gang vector       
       do ip = 1, npuv
          ! 
          timestep_analysis_required_timestep(ip) = dtmax ! Reset per-cell limits; dry cells will retain dtmax
          !
       enddo
       !$acc end parallel
       !$omp end do
       !$omp end parallel       
       !
   endif   
   !
   ! For some reason, it is necessary to set num_gangs here! Without, the program launches only 1 gang, and everything becomes VERY slow!
   !
   ! Flux and velocity from the previous time step. The arrays are swapped instead of copied:
   ! q0 and uv0 now hold the previous time step, q and uv hold older values and are overwritten
   ! below at every point. Points that are not computed here (open boundaries, structures, wave
   ! makers) get the previous-step value, as with a copy. The dummy point npuv + ncuv + 1 is 0.0
   ! in both sets of arrays.
   !
   call move_alloc(q0, swap)
   call move_alloc(q, q0)
   call move_alloc(swap, q)
   !
   call move_alloc(uv0, swap)
   call move_alloc(uv, uv0)
   call move_alloc(swap, uv)
   !
   !$acc parallel, present( kcuv, kfuv, zs, q, q0, uv, uv0, zsderv, z_wetfrac, &
   !$acc                    uv_flags_iref, uv_flags_type, uv_flags_dir, mask_adv, &
   !$acc                    subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_havg, subgrid_uv_nrep, subgrid_uv_pwet, &
   !$acc                    subgrid_uv_havg_zmax, subgrid_uv_nrep_zmax, subgrid_uv_fnfit, subgrid_uv_navg_w, &
   !$acc                    uv_index_z_nm, uv_index_z_nmu, uv_index_u_nmd, uv_index_u_nmu, uv_index_u_ndm, uv_index_u_num, &
   !$acc                    uv_index_v_ndm, uv_index_v_ndmu, uv_index_v_nm, uv_index_v_nmu, cuv_index_uv, cuv_index_uv1, cuv_index_uv2, &
   !$acc                    zb, zbuv, zbuvmx, tauwu, tauwv, patm, fwuv, gn2uv, dxminv, dxrinv, dyrinv, dxm2inv, dxr2inv, dyr2inv, &
   !$acc                    dxrinvc, dyrinvc, fcorio2d, nuvisc, z_volume, gnapp2, timestep_analysis_required_timestep ) num_gangs( 1024 ) vector_length( 128 )
   !$omp parallel &
   !$omp private ( ip,hu,qfr,qsm,qx_nm,nm,nmu,dzdx,frc,idir,itype,iref,dxuvinv,dxuv2inv,dyuvinv,dyuv2inv, &
   !$omp           qx_nmd,qx_nmu,qy_nm,qy_ndm,qy_nmu,qy_ndmu,uu_nm,uu_nmd,uu_nmu,uu_num,uu_ndm,vu, & 
   !$omp           fcoriouv,gnavg2,iok,zsu,dzuv,iuv,facint,fwmax,zmax,zmin,one_minus_facint,dqxudx,dqyudy,uu,ud,qu,qd,qy,hwet,phi,adv,mdrv,phiz,fac,hu73,y_cbrt,i_cbrt,min_dt_ip, &
   !$omp           q_abs,q_k,q_ref,wt_nmd,wt_nmu,wt_ndm,wt_num,iup ) &
   !$omp reduction ( min : min_dt  )
   !$omp do schedule ( dynamic, 256 )
   !$acc loop, reduction( min : min_dt ), gang, vector
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
         ! Water level at the uv point: the higher of the two water levels (default). With
         ! bates_upwind_level the level of the upwind cell is used instead, i.e. the
         ! cell the water comes from according to the sign of the previous-step velocity. The
         ! two differ where water flows against the surface slope (late uprush, run-down).
         !
         zsu = max(zs(nm), zs(nmu))
         !
         if (bates_upwind_level) then
            !
            if (uv0(ip) > 1.0e-6) then
               !
               iup = nm
               zsu = zs(iup)
               !
            elseif (uv0(ip) < -1.0e-6) then
               !
               iup = nmu
               zsu = zs(iup)
               !
            endif
            !
         endif
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
                     dxuvinv  = dyrinvc(iref)
                     dyuvinv  = dxrinv(iref)
                     dxuv2inv = 0.0 ! no viscosity
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
                  dzuv   = (zmax - zmin) / (subgrid_nlevels - 1)                                                          ! level size (is storing this in memory faster?)
                  iuv    = min(int((zsu - zmin) / dzuv) + 1, subgrid_nlevels - 1)                                         ! index of level below zsu 
                  facint = (zsu - (zmin + (iuv - 1) * dzuv) ) / dzuv                                                        ! 1d interpolation coefficient
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
                     ! Weights of the four partner faces (advection_partner_weight, off by default in
                     ! this scheme). A dry partner face has u = 0, so the velocity differences below treat
                     ! every staircase step of a bank as a real 90 degree corner and slow the face down.
                     ! Each difference is therefore scaled with the larger of the depth ratio and the flux
                     ! ratio between partner and face (at most 1):
                     ! wt = min(1, max(hu_k / hu, |q_k| / |q|)) with the partner depth hu_k = q_k / u_k,
                     ! written as |q_k| / q_ref with q_ref = min(|u_k| * hu, |q|). Only the terms with du/dx
                     ! and du/dy are weighted. The terms with dqu/dx and dqv/dy add up to u times the flux
                     ! divergence of the upwind cell, which is zero in steady flow.
                     !
                     wt_nmd = 1.0
                     wt_nmu = 1.0
                     wt_ndm = 1.0
                     wt_num = 1.0
                     !
                     if (advection_partner_weight) then
                        !
                        q_abs = max(abs(qx_nm), 1.0e-9)
                        !
                        q_k   = abs(qx_nmd)
                        q_ref = min(abs(uu_nmd) * hu, q_abs)
                        if (q_k < q_ref)  wt_nmd = q_k / q_ref
                        if (q_k <= 0.0)   wt_nmd = 0.0
                        !
                        q_k   = abs(qx_nmu)
                        q_ref = min(abs(uu_nmu) * hu, q_abs)
                        if (q_k < q_ref)  wt_nmu = q_k / q_ref
                        if (q_k <= 0.0)   wt_nmu = 0.0
                        !
                        q_k   = abs(q0(uv_index_u_ndm(ip)))
                        q_ref = min(abs(uu_ndm) * hu, q_abs)
                        if (q_k < q_ref)  wt_ndm = q_k / q_ref
                        if (q_k <= 0.0)   wt_ndm = 0.0
                        !
                        q_k   = abs(q0(uv_index_u_num(ip)))
                        q_ref = min(abs(uu_num) * hu, q_abs)
                        if (q_k < q_ref)  wt_num = q_k / q_ref
                        if (q_k <= 0.0)   wt_num = 0.0
                        !
                     endif
                     !
                     ! d qu u / dx
                     !
                     qd = (qx_nmd + qx_nm) / 2
                     qu = (qx_nm + qx_nmu) / 2
                     !
                     if (qd > 1.0e-6) then
                        ud = (uu_nmd + uu_nm) / 2
                        dqxudx = ( qd * wt_nmd * (uu_nm - uu_nmd) + ud * (qx_nm - qx_nmd) ) * dxuvinv
                     endif
                     !
                     if (qu < -1.0e-6) then
                        uu = (uu_nm + uu_nmu) / 2
                        dqxudx = dqxudx + ( qu * wt_nmu * (uu_nmu - uu_nm) + uu * (qx_nmu - qx_nm) ) * dxuvinv
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
                        dqyudy = qd * wt_ndm * (uu_nm - uu_ndm) * dyuvinv
                     endif
                     !
                     if (qu < -1.0e-6) then
                        dqyudy = dqyudy + qu * wt_num * (uu_num - uu_nm) * dyuvinv
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
                  ! Limit advection term such that horizontal acceleration due to advection does not exceed advlim
                  ! (default 9999.9 m/s2, i.e. no limit)
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
               ! old: fwmax = 0.8 * hwet * sqrt(hwet) / 15
               ! fix for lab cases: fwmax = 999
               !
               fwmax = fwmaxfac * hwet * sqrt(hwet)     
               ! Note, fwmaxfac is determined in sfincs_snapwave every 'update_wave_field' call
               !
               frc = frc + phi * sign(min(abs(fwuv(ip)), fwmax), fwuv(ip))               
               !
            endif
            !
            ! hu**(1/3) for the Manning friction, via a fast cube root: an integer bit-hack seed
            ! (magic constant 709921077) refined by one Newton iteration
            ! (y <- y - (y^3 - hu) / (3 y^2)). About 0.1 % accurate over the depth range, no pow,
            ! no table. y_cbrt gives hu**(7/3) below and hu**(5/3) in the newly-wet estimate.
            !
            i_cbrt = transfer(hu, i_cbrt)
            i_cbrt = i_cbrt / 3 + 709921077
            y_cbrt = transfer(i_cbrt, y_cbrt)
            y_cbrt = y_cbrt - (y_cbrt * y_cbrt * y_cbrt - hu) / (3.0 * y_cbrt * y_cbrt)
            !
            ! Compute flux qfr used for friction term
            !
            if (kfuv(ip) == 0) then
               !
               ! This uv point just became wet, so estimate equilibrium flux (hu**(5/3) = hu * y_cbrt**2)
               !
               qfr = sqrt(abs(dzdx) / (max(gnavg2, 1.0e-5) / 10)) * hu * y_cbrt * y_cbrt
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
            hu73 = hu * hu * y_cbrt   ! hu**(7/3)
            !
            q(ip) = (qsm + frc * dt) / (1.0 + gnavg2 * dt * qfr / hu73)
            !
            if (subgrid .and. wiggle_suppression) then
               !
               ! Wiggle suppression for cells in which only a few subgrid pixels are wet (large dzs/dV).
               ! There the local gravity-wave speed is sqrt(g * hu / phiz) with phiz the wet-area
               ! fraction of the cell, which exceeds the speed sqrt(g * hu) the time step was set for.
               ! Under-relax the flux increment with fac = phiz / alfa**2 (floored at wiggle_facmin).
               ! Steady flow (q = q0) is not affected, 2*dt modes are damped. Fully wet cells give fac = 1.
               ! With wiggle_detect, only relax when the level accelerations of the two cells diverge
               ! by more than wiggle_threshold (anti-phase sloshing across this uv point).
               !
               phiz = min(z_wetfrac(nm), z_wetfrac(nmu))
               !
               fac  = min(max(phiz / alfa**2, wiggle_facmin), 1.0)
               !
               if (fac < 1.0) then
                  !
                  if (wiggle_detect) then
                     !
                     mdrv = abs(zsderv(nm) - zsderv(nmu)) - wiggle_threshold
                     !
                     if (mdrv > 0.0) then
                        q(ip) = qsm + fac * (q(ip) - qsm)
                     endif
                     !
                  else
                     !
                     q(ip) = qsm + fac * (q(ip) - qsm)
                     !
                  endif
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
            ! Compute wet-averaged velocity. hu can become extremely small in subgrid mode,
            ! so use huvmin to prevent unrealistically high velocities.
            !
            uv(ip) = q(ip) / max(hu, huvmin)
            !
            kfuv(ip) = 1
            !
            ! Determine minimum time step (alpha is added later on in sfincs_lib.f90) of all uv points
            ! Use maximum of sqrt(gh) and current velocity
            !
            min_dt_ip = 1.0 / ( max(sqrt(g * hu), abs(uv(ip)) ) * dxuvinv)
            !
            min_dt = min(min_dt, min_dt_ip)
            !
            ! Compute timestep per grid cell
            !
            if (timestep_analysis) then
                !
                timestep_analysis_required_timestep(ip) = min_dt_ip
                !
                
            endif            
            !
         else
            !
            q(ip)  = 0.0
            uv(ip) = 0.0
            kfuv(ip) = 0
            !
         endif
         !
      else
         !
         ! Not computed here (open boundary, structure, wave maker or inactive point): keep the
         ! value of the previous time step
         !
         q(ip)  = q0(ip)
         uv(ip) = uv0(ip)
         !
      endif
   enddo
   !$omp end do
   !$omp end parallel
   !$acc end parallel
   !
   if (ncuv > 0) then
      !
      ! Loop through combined uv points and determine average uv and q
      ! The combined q and uv values are used in the continuity equation and in the netcdf output
      !
      !$omp parallel &
      !$omp private ( icuv )
      !$omp do
      !$acc parallel, present( q, uv, cuv_index_uv, cuv_index_uv1, cuv_index_uv2  )
      !$acc loop gang vector
      do icuv = 1, ncuv
         !
         ! Average of the two uv points
         !
         q(cuv_index_uv(icuv))  = (q(cuv_index_uv1(icuv)) + q(cuv_index_uv2(icuv))) / 2
         uv(cuv_index_uv(icuv)) = (uv(cuv_index_uv1(icuv)) + uv(cuv_index_uv2(icuv))) / 2
         !
      enddo
      !$acc end parallel
      !$omp end do
      !$omp end parallel
      !
   endif
   !
   call timer_stop('momentum')
   !
   end subroutine
   !
end module
