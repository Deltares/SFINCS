module sfincs_momentum_velocity
   !
   use sfincs_data
   !
   implicit none
   !
contains
   !
   subroutine compute_fluxes_velocity(dt, tloop)
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
   !
   integer   :: idir
   integer   :: iref
   integer   :: itype
   integer   :: iuv
   integer   :: icuv
   !
   real*4    :: hu
   real*4    :: dxuvinv
   real*4    :: dxuv2inv
   real*4    :: dyuvinv
   real*4    :: dyuv2inv
   real*4    :: adv
   real*4    :: fcoriouv
   real*4    :: frc
   !
   real*4    :: qfr
   !
   real*4    :: uu_nm
   real*4    :: uu_nmd
   real*4    :: uu_nmu
   real*4    :: uu_ndm
   real*4    :: uu_num
   real*4    :: vu
   !
   real*4    :: zsu
   real*4    :: dzuv
   real*4    :: facint
   real*4    :: gnavg2
   real*4    :: fwmax
   real*4    :: zmax
   real*4    :: zmin
   !
   real*4    :: dqxudx
   real*4    :: dqyudy
   real*4    :: qu
   real*4    :: qd
   real*4    :: uu
   real*4    :: ud
   real*4    :: dzdx
   !
   real*4    :: hwet
   real*4    :: phi
   !
   real*4    :: hu73
   !
   real*4    :: zs2w, zs1e, dnm, dnmu, zrec     ! NEOWAVE advection work vars
   real*4    :: zbavg                           ! NEOWAVE average still bed at u-point = 0.5*(zb(nm)+zb(nmu))
   integer   :: ipw, ipe
   !
   real*4    :: min_dt_ip
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
   ! Copy velocity from previous time step (velocity form carries uv, not q)
   !
   !$acc parallel, present( uv, uv0 )
   !$omp parallel &
   !$omp private ( ip )
   !$omp do
   !$acc loop gang vector
   do ip = 1, npuv + ncuv
      !
      uv0(ip) = uv(ip)
      !
   enddo
   !$acc end parallel
   !$omp end do
   !$omp end parallel
   !
   !$acc parallel, present( kcuv, kfuv, zs, q, uv, uv0, zsderv, &
   !$acc                    uv_flags_iref, uv_flags_type, uv_flags_dir, mask_adv, &
   !$acc                    subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_havg, subgrid_uv_nrep, subgrid_uv_pwet, &
   !$acc                    subgrid_uv_havg_zmax, subgrid_uv_nrep_zmax, subgrid_uv_fnfit, subgrid_uv_navg_w, &
   !$acc                    uv_index_z_nm, uv_index_z_nmu, uv_index_u_nmd, uv_index_u_nmu, uv_index_u_ndm, uv_index_u_num, &
   !$acc                    uv_index_v_ndm, uv_index_v_ndmu, uv_index_v_nm, uv_index_v_nmu, cuv_index_uv, cuv_index_uv1, cuv_index_uv2, &
   !$acc                    zb, zbuv, zbuvmx, tauwu, tauwv, patm, fwuv, gn2uv, dxminv, dxrinv, dyrinv, dxm2inv, dxr2inv, dyr2inv, &
   !$acc                    dxrinvc, dyrinvc, fcorio2d, nuvisc, z_volume, cell_area, cell_area_m2, z_flags_iref, gnapp2, x73, timestep_analysis_required_timestep ) num_gangs( 1024 ) vector_length( 128 )
   !$omp parallel &
   !$omp private ( ip,hu,qfr,nm,nmu,dzdx,frc,idir,itype,iref,dxuvinv,dxuv2inv,dyuvinv,dyuv2inv, &
   !$omp           uu_nm,uu_nmd,uu_nmu,uu_num,uu_ndm,vu, &
   !$omp           fcoriouv,gnavg2,iok,zsu,dzuv,iuv,facint,fwmax,zmax,zmin,dqxudx,dqyudy,uu,ud,qu,qd,hwet,phi,adv,hu73,min_dt_ip,zs2w,zs1e,dnm,dnmu,zrec,zbavg,ipw,ipe ) &
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
         ! NEOWAVE upwind surface at the u-point (Yamazaki et al. 2009): take the surface from
         ! the cell the flow comes from (sign of the previous-step velocity); tie-break to the
         ! higher surface. This is the upwind zeta of the Mader flux (eq 13).
         !
         if (uv0(ip) > 1.0e-6) then
            zsu = zs(nm)
         elseif (uv0(ip) < -1.0e-6) then
            zsu = zs(nmu)
         else
            zsu = max(zs(nm), zs(nmu))
         endif
         !
         if (subgrid) then
            !
            zmin = subgrid_uv_zmin(ip)
            zmax = subgrid_uv_zmax(ip)
            !
            if (zsu > zmin) then ! In the 'new' subgrid formulations, zmin is lowest pixel + huthresh. Huthresh was already applied when building the subgrid tables. In sfincs_domain, huthresh is set to 0.0 to acount for this.
               iok = .true.
            endif
            !
         else
            !
            ! NEOWAVE flow depth at the u-point: D = upwind surface + average still depth (h = -zb),
            ! i.e. zsu - 0.5*(zb(nm)+zb(nmu)). Uses the AVERAGE bed (not the max as in the Bates
            ! conveyance), so the waterline can advance up a slope -> stronger run-up.
            !
            zbavg = 0.5 * (zb(nm) + zb(nmu))
            !
            if (zsu - zbavg > huthresh) then
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
            ! Get velocities from the previous time step
            !
            if (advection .or. coriolis .or. viscosity .or. friction2d) then
               !
               ! Get the neighbors
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
               hu     = max(zsu - zbavg, huthresh)   ! NEOWAVE flow depth D = upwind zeta + avg still depth
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
            adv = 0.0
            !
            if (advection) then
               !
               if (mask_adv(ip) == 1) then
                  !
                  ! NEOWAVE momentum-conserved advection (Yamazaki, Kowalik & Cheung 2009,
                  ! eqs 18/20/22), VELOCITY form. Streamwise advective speeds from the Mader
                  ! upwind-zeta flux  FLU = mean(U)*(upwind surface zeta + still-depth h), h=-zb;
                  ! zeta reconstructed 2nd-order (upwind face, +/-2 stencil) when co-directional,
                  ! else 1st-order. Cross term advects u by the transverse velocity vu (upwind).
                  ! 'adv' is a VELOCITY tendency [m/s^2] (added to frc/hu below).
                  !
                  ! Cell flow depths D_nm (west), D_nmu (east). NON-subgrid: D = max(zeta - zb, 0)
                  ! with the Mader upwind-reconstructed surface zeta (2nd-order when co-directional).
                  ! SUBGRID: zb is not a valid conveyance bed -> use the subgrid cell-averaged water
                  ! depth from the cell volume, D = z_volume / cell_area (the upwind-zeta surface
                  ! reconstruction is dropped; the subgrid tables already supply dissipation).
                  !
                  if (subgrid) then
                     !
                     if (crsgeo) then
                        dnm  = z_volume(nm)  / cell_area_m2(nm)
                        dnmu = z_volume(nmu) / cell_area_m2(nmu)
                     else
                        dnm  = z_volume(nm)  / cell_area(z_flags_iref(nm))
                        dnmu = z_volume(nmu) / cell_area(z_flags_iref(nmu))
                     endif
                     dnm  = max(dnm,  0.0)
                     dnmu = max(dnmu, 0.0)
                     qd   = 0.5 * (uu_nmd + uu_nm) * dnm                       ! FLU_p (west cell)
                     qu   = 0.5 * (uu_nm + uu_nmu) * dnmu                      ! FLU_n (east cell)
                     !
                  else
                     !
                     ipw = uv_index_u_nmd(ip)
                     ipe = uv_index_u_nmu(ip)
                     zs2w = zs(nm)  ; if (ipw > 0) zs2w = zs(uv_index_z_nm(ipw))
                     zs1e = zs(nmu) ; if (ipe > 0) zs1e = zs(uv_index_z_nmu(ipe))
                     if (uu_nmd >= 0.0) then
                        zrec = 0.5 * (zs2w + zs(nm))
                     else
                        zrec = zs(nm)
                     endif
                     qd = 0.5 * (uu_nmd + uu_nm) * max(zrec - zb(nm), 0.0)     ! FLU_p (west cell)
                     if (uu_nmu <= 0.0) then
                        zrec = 0.5 * (zs(nmu) + zs1e)
                     else
                        zrec = zs(nmu)
                     endif
                     qu = 0.5 * (uu_nm + uu_nmu) * max(zrec - zb(nmu), 0.0)    ! FLU_n (east cell)
                     dnm  = max(zs(nm)  - zb(nm),  0.0)
                     dnmu = max(zs(nmu) - zb(nmu), 0.0)
                     !
                  endif
                  !
                  dnm = max(dnm + dnmu, huthresh)       ! D_nm + D_nmu
                  ud  = max(2.0 * qd / dnm, 0.0)        ! U_p  (advective speed, from west)
                  uu  = min(2.0 * qu / dnm, 0.0)        ! U_n  (advective speed, from east)
                  dqxudx = ( ud * (uu_nm - uu_nmd) + uu * (uu_nmu - uu_nm) ) * dxuvinv
                  !
                  if (vu >= 0.0) then                   ! cross term: u advected by vu (upwind)
                     dqyudy = vu * (uu_nm - uu_ndm) * dyuvinv
                  else
                     dqyudy = vu * (uu_num - uu_nm) * dyuvinv
                  endif
                  !
                  adv = - phi * (dqxudx + dqyudy)        ! velocity tendency [m/s^2]
                  adv = min(max(adv, -advlim), advlim)   ! limit advective acceleration
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
               fwmax = 0.8 * hwet * sqrt(hwet) / 15
               !
               frc = frc + phi * sign(min(abs(fwuv(ip)), fwmax), fwuv(ip))
               !
            endif
            !
            ! Convert the accumulated flux forces (pressure, viscosity, coriolis, wind, atm,
            ! wave) to a VELOCITY tendency (/hu) and add the velocity-form advection -> the whole
            ! momentum equation is now velocity-form (NEOWAVE); uv is updated directly below.
            !
            frc = frc / max(hu, huvmin) + adv
            !
            ! Friction flux proxy qfr (velocity form: qfr = |u|*hu, so the implicit factor
            ! gnavg2*qfr/hu^(7/3) = gnavg2*|u|/hu^(4/3) is the velocity-form Manning term; no q0).
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
                  ! Both velocity components: qfr = hu*sqrt(u^2 + v^2)
                  !
                  qfr = hu * sqrt(uv0(ip)**2 + vu**2)
                  !
               else
                  !
                  ! Streamwise velocity only: qfr = |u|*hu
                  !
                  qfr = abs(uv0(ip)) * hu
                  !
               endif
               !
            endif
            !
            ! Update the velocity for this uv point (Yamazaki velocity form; semi-implicit Manning)
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
            ! VELOCITY-form update: solve uv directly. frc is now a velocity tendency;
            ! the implicit friction factor matches the flux form (|q|/hu^(7/3) = |u|/hu^(4/3)).
            !
            uv(ip) = (uv0(ip) + frc * dt) / (1.0 + gnavg2 * dt * qfr / hu73)
            !
            ! velocity limiter (default 10 m/s)
            uv(ip) = min(max(uv(ip), - uvlim), uvlim)
            !
            ! no flow out of a cell that is (going) dry
            if (zs(nm)  < zb(nm))  uv(ip) = min(uv(ip), 0.0)
            if (zs(nmu) < zb(nmu)) uv(ip) = max(uv(ip), 0.0)
            !
            ! Continuity flux from the updated velocity and the conveyance depth.
            !
            q(ip) = uv(ip) * hu
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
