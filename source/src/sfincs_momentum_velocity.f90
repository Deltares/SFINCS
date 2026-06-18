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
   real*4    :: ufr
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
   real*4    :: vp, vn                          ! flux-based cross-advective speeds (Yamazaki eq. 22)
   real*4    :: qu
   real*4    :: qd
   real*4    :: un                              ! U_n advective speed (from east)
   real*4    :: up                              ! U_p advective speed (from west)
   real*4    :: dnminv                          ! 1/(D_nm + D_nmu) reused for both advective speeds
   real*4    :: dzdx
   !
   real*4    :: hwet
   real*4    :: phi
   !
   real*4    :: hu43
   real*4    :: y_cbrt                          ! cube-root approximation hu^(1/3)
   integer*4 :: i_cbrt                          ! bit pattern of hu/y_cbrt for the cube-root seed
   !
   real*4    :: zs2w, zs1e, dnm, dnmu, zrec     ! advection work vars
   real*4    :: zbup                            ! upwind still bed at u-point (bed of the cell the flow comes from)
   integer   :: ipw, ipe
   real*4    :: zbnm, zbnmu                     ! bed at the west/east cell for subgrid advection (not necessarily zb)
   real*4    :: mdrv                            ! subgrid wiggle suppression driver
   !
   real*4    :: min_dt_ip
   !
   logical   :: iwet
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
   ! Copy velocity and flux from the previous time step (the velocity form
   ! advects uv0; q0 provides the consistent previous-step fluxes for the
   ! momentum-conserving cross-advection speeds)
   !
   !$acc parallel, present( uv, uv0, q, q0 )
   !$omp parallel &
   !$omp private ( ip )
   !$omp do
   !$acc loop gang vector
   do ip = 1, npuv + ncuv
      !
      uv0(ip) = uv(ip)
      q0(ip)  = q(ip)
      !
   enddo
   !$acc end parallel
   !$omp end do
   !$omp end parallel
   !
   !$omp parallel &
   !$omp private ( ip,hu,ufr,nm,nmu,dzdx,frc,idir,itype,iref,dxuvinv,dxuv2inv,dyuvinv,dyuv2inv, &
   !$omp           uu_nm,uu_nmd,uu_nmu,uu_num,uu_ndm,vu, &
   !$omp           fcoriouv,gnavg2,iwet,zsu,dzuv,iuv,facint,fwmax,zmax,zmin,dqxudx,dqyudy,un,up,vp,vn, &
   !$omp           dnminv,qu,qd,hwet,phi,adv,hu43,y_cbrt,i_cbrt,min_dt_ip,zs2w,zs1e,dnm,dnmu,zrec,zbup,ipw,ipe,zbnm,zbnmu ) &
   !$omp reduction ( min : min_dt  )
   !$omp do schedule ( dynamic, 256 )
   !$acc parallel, present( kcuv, kfuv, zs, q, q0, uv, uv0, &
   !$acc                    uv_flags_iref, uv_flags_type, uv_flags_dir, mask_adv, &
   !$acc                    subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_havg, subgrid_uv_nrep, subgrid_uv_pwet, &
   !$acc                    subgrid_uv_havg_zmax, subgrid_uv_nrep_zmax, subgrid_uv_fnfit, subgrid_uv_navg_w, &
   !$acc                    uv_index_z_nm, uv_index_z_nmu, uv_index_u_nmd, uv_index_u_nmu, uv_index_u_ndm, uv_index_u_num, &
   !$acc                    uv_index_v_ndm, uv_index_v_ndmu, uv_index_v_nm, uv_index_v_nmu, cuv_index_uv, cuv_index_uv1, cuv_index_uv2, &
   !$acc                    zb, tauwu, tauwv, patm, fwuv, gn2uv, dxminv, dxrinv, dyrinv, dxm2inv, dxr2inv, dyr2inv, &
   !$acc                    dxrinvc, dyrinvc, fcorio2d, nuvisc, z_volume, cell_area, cell_area_m2, z_flags_iref, gnapp2, timestep_analysis_required_timestep ) num_gangs( 1024 ) vector_length( 128 )
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
         iwet  = .false.
         !
         ! Upwind surface at the u-point: take the surface from the cell the flow comes from
         ! (sign of the previous-step velocity). The upwind bed (zbup) is only needed for the
         ! regular-grid conveyance and is set in the non-subgrid branch below.
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
            if (zsu > zmin) then ! In the subgrid formulations, zmin is lowest pixel + huthresh. Huthresh was already applied when building the subgrid tables. In sfincs_domain, huthresh is set to 0.0 to acount for this.
               iwet = .true.
            endif
            !
         else
            !
            ! Flow depth at the u-point: D = upwind surface - upwind bed = zsu - zbup, i.e. the water
            ! depth in the cell the flow comes from. The upwind bed follows the same flow direction
            ! that selected zsu. The face is wet when that upwind depth exceeds huthresh. This avoids
            ! the average-bed depth overshoot on steep downslopes while still allowing run-up.
            !
            if (uv0(ip) > 1.0e-6) then
               zbup = zb(nm)
            elseif (uv0(ip) < -1.0e-6) then
               zbup = zb(nmu)
            else
               if (zs(nm) >= zs(nmu)) then
                  zbup = zb(nm)
               else
                  zbup = zb(nmu)
               endif
            endif
            !
            if (zsu - zbup > huthresh) then
               iwet = .true.
            endif
            !
         endif
         !
         if (iwet) then
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
               hu     = zsu - zbup    ! Flow depth D = upwind zeta - upwind bed
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
            ! Velocity form: build frc directly as an acceleration [m/s^2]. Forces that scale
            ! with depth (pressure, viscosity, Coriolis, atm) are written WITHOUT hu -- their hu
            ! would only be divided out again. Only the surface stresses (wind, waves) keep a /hu.
            !
            frc = - g * dzdx
            !
            if (advection) then
               !
               if (mask_adv(ip) == 1) then
                  !
                  ! Momentum-conserved advection (Yamazaki, Kowalik & Cheung 2009,
                  ! eqs 18/20/22), VELOCITY form. Streamwise advective speeds from the Mader
                  ! upwind-zeta flux  FLU = mean(U)*(upwind surface zeta + still-depth h), h=-zb;
                  ! zeta reconstructed 2nd-order (upwind face, +/-2 stencil) when co-directional,
                  ! else 1st-order. Cross term: flux-based two-sided upwind (eq. 22, see below).
                  ! 'adv' is a VELOCITY tendency [m/s^2] (added straight into frc).
                  !
                  ! One path for regular and subgrid bathymetry: on subgrid models zb
                  ! is the EFFECTIVE bed zs - z_volume/area maintained every step in
                  ! continuity (zb_effective; the file zb is not a valid conveyance
                  ! bed), so D = zs - zb is the subgrid cell-mean depth there. Then
                  ! zrec - zbnm = D_cell + 0.5*(upwind dzs) -- centered depth when
                  ! flat, with the upwind-surface boost at a front (subgrid-consistent
                  ! Mader reconstruction).
                  !
                  dnm   = max(zs(nm)  - zb(nm),  0.0)
                  dnmu  = max(zs(nmu) - zb(nmu), 0.0)
                  !
                  zbnm  = zb(nm)
                  zbnmu = zb(nmu)
                  !
                  ipw = uv_index_u_nmd(ip)
                  ipe = uv_index_u_nmu(ip)
                  !
                  ! Only use the +/-2 stencil when the neighbor is a REGULAR uv point
                  ! (<= npuv). At a quadtree refinement boundary uv_index_u_nm* points to a
                  ! COMBINED uv point (index > npuv), for which uv_index_z_* is out of bounds
                  ! (those arrays are sized npuv); fall back to 1st order there. Note ipw/ipe
                  ! are never 0 (sfincs_domain sets a missing neighbor to ip), so the old
                  ! "> 0" test never triggered the fallback.
                  !
                  if (ipw > 0 .and. ipw <= npuv) then
                     zs2w = zs(uv_index_z_nm(ipw))
                  else
                     zs2w = zs(nm)
                  endif
                  !
                  if (ipe > 0 .and. ipe <= npuv) then
                     zs1e = zs(uv_index_z_nmu(ipe))
                  else
                     zs1e = zs(nmu)
                  endif
                  !
                  if (uu_nmd >= 0.0) then
                     zrec = 0.5 * (zs2w + zs(nm))
                  else
                     zrec = zs(nm)
                  endif
                  !
                  qd = 0.5 * (uu_nmd + uu_nm) * max(zrec - zbnm, 0.0)       ! FLU_p (west cell)
                  !
                  if (uu_nmu <= 0.0) then
                     zrec = 0.5 * (zs(nmu) + zs1e)
                  else
                     zrec = zs(nmu)
                  endif
                  !
                  qu = 0.5 * (uu_nm + uu_nmu) * max(zrec - zbnmu, 0.0)      ! FLU_n (east cell)
                  !
                  dnm    = max(dnm + dnmu, huthresh)    ! D_nm + D_nmu
                  dnminv = 2.0 / dnm
                  !
                  up  = max(qd * dnminv, 0.0)           ! U_p  (advective speed, from west)
                  un  = min(qu * dnminv, 0.0)           ! U_n  (advective speed, from east)                  
                  !
                  dqxudx = ( up * (uu_nm - uu_nmd) + un * (uu_nmu - uu_nm) ) * dxuvinv
                  !
                  ! Cross-advection of u by v -- momentum-conserving two-sided upwind
                  ! (Yamazaki et al. 2009, eq. 22): advective speeds from the y-FLUXES
                  ! through the faces below/above the u-point (mean of the two flanking
                  ! v-point fluxes, previous time step q0), normalized by the same total
                  ! depth as the streamwise term (dnminv). Transport INTO the point from
                  ! either side contributes; at a collision line (v converging from both
                  ! sides) both terms stay active, where the velocity-average form below
                  ! gave ~zero cross-advection.
                  !
                  vp = max( 0.5 * (q0(uv_index_v_ndm(ip)) + q0(uv_index_v_ndmu(ip))) * dnminv, 0.0 )
                  vn = min( 0.5 * (q0(uv_index_v_nm(ip))  + q0(uv_index_v_nmu(ip)))  * dnminv, 0.0 )
                  !
                  dqyudy = ( vp * (uu_nm - uu_ndm) + vn * (uu_num - uu_nm) ) * dyuvinv
                  !
                  adv = - phi * (dqxudx + dqyudy)        ! velocity tendency [m/s^2]
                  !
                  frc = frc + min(max(adv, -advlim), advlim)   ! add limited advective acceleration
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
                  frc = frc + nuvisc(iref) * ( (uu_nmu - 2*uu_nm + uu_nmd ) * dxuv2inv + (uu_num - 2*uu_nm + uu_ndm ) * dyuv2inv )
                  !
               else
                  !
                  ! Increase viscosity to prevent instabilities on refinement (related to advection?)
                  !
                  frc = frc + nuviscfac * nuvisc(iref) * ( (uu_nmu - 2*uu_nm + uu_nmd ) * dxuv2inv + (uu_num - 2*uu_nm + uu_ndm ) * dyuv2inv )
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
                  frc = frc + fcoriouv * vu ! U
                  !
               else
                  !
                  frc = frc - fcoriouv * vu ! V
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
                  ! Wind stress is a surface force -> acceleration = stress / depth
                  !
                  if (idir==0) then
                     !
                     frc = frc + phi * tauwu(nm) / max(hu, huvmin)
                     !
                  else
                     !
                     frc = frc + phi * tauwv(nm) / max(hu, huvmin)
                     !
                  endif
                  !
               else
                  !
                  ! Reduce wind drag at water depths < 0.25 m (tauw*hu*4 / hu = tauw*4)
                  !
                  if (idir==0) then
                     !
                     frc = frc + tauwu(nm) * 4
                     !
                  else
                     !
                     frc = frc + tauwv(nm) * 4
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
               frc = frc + (patm(nm) - patm(nmu)) * dxuvinv / rhow
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
               ! Wave force is a surface force -> acceleration = force / depth
               !
               frc = frc + phi * sign(min(abs(fwuv(ip)), fwmax), fwuv(ip)) / max(hu, huvmin)
               !
            endif
            !
            ! hu**(1/3) and hu**(4/3) for the velocity-form Manning friction, via a fast cube
            ! root: an integer bit-hack seed (magic constant 709921077) refined by one Newton
            ! iteration (y <- y - (y^3 - hu)/(3 y^2)). ~0.1% accurate over the depth range,
            ! no pow, no table. y_cbrt = hu^(1/3) is reused for the newly-wet estimate below.
            !
            i_cbrt = transfer(hu, i_cbrt)
            i_cbrt = i_cbrt / 3 + 709921077
            y_cbrt = transfer(i_cbrt, y_cbrt)
            y_cbrt = y_cbrt - (y_cbrt * y_cbrt * y_cbrt - hu) / (3.0 * y_cbrt * y_cbrt)
            hu43   = hu * y_cbrt
            !
            ! Friction velocity proxy ufr (velocity form: the implicit Manning factor is
            ! gnavg2*ufr/hu^(4/3) with ufr the friction-driving velocity magnitude).
            !
            if (kfuv(ip) == 0) then
               !
               ! This uv point just became wet, so estimate the equilibrium velocity
               ! (hu^(2/3) = (hu^(1/3))^2 = y_cbrt^2, reusing the cube root above)
               !
               ufr = sqrt(abs(dzdx) / (max(gnavg2, 1.0e-5) / 10)) * y_cbrt * y_cbrt
               !
            else
               !
               if (friction2d) then
                  !
                  ! Both velocity components: ufr = sqrt(u^2 + v^2)
                  !
                  ufr = sqrt(uv0(ip)**2 + vu**2)
                  !
               else
                  !
                  ! Streamwise velocity only: ufr = |u|
                  !
                  ufr = abs(uv0(ip))
                  !
               endif
               !
            endif
            !
            ! Velocity update. frc is a velocity tendency and the
            ! implicit friction factor is the velocity-form Manning term gnavg2*|u|/hu^(4/3).
            !
            uv(ip) = (uv0(ip) + frc * dt) / (1.0 + gnavg2 * dt * ufr / hu43)
            !
            if (subgrid .and. wiggle_suppression) then 
               !
               ! If the acceleration of water level in cell nm is large and positive and in nmu large and negative, or vice versa, apply limiter to the flux. Only for subgrid.
               !
               mdrv = abs(zsderv(nm) - zsderv(nmu)) - wiggle_threshold
               !
               if (mdrv > 0.0) then
                  !
                  uv(ip) = uv(ip) * wiggle_threshold / (wiggle_factor * mdrv + wiggle_threshold)
                  !
               endif
               !
            endif
            !
            ! Velocity limiter (default 10 m/s)            
            !
            uv(ip) = min(max(uv(ip), - uvlim), uvlim)
            !
            ! No flow out of a cell that is (going) dry
            !
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
end module
