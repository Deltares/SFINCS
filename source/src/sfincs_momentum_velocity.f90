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
   real*4    :: umax                            ! local characteristic speed bound sqrt(g*hu) + |u|
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
   real*4    :: phiz                            ! wet-area fraction of the least wet neighbouring cell
   real*4    :: fac                             ! under-relaxation factor for wiggle suppression
   !
   real*4    :: min_dt_ip
   !
   integer   :: iup                             ! upwind cell of the uv point
   integer   :: idn                             ! downstream cell of the uv point
   real*4    :: w                               ! face regime weight (0 level-driven, 1 slope-driven)
   real*4    :: s_s                             ! bed-plane slope of the upwind cell along the face normal
   real*4    :: hwet_up                         ! level-based wet depth of the upwind cell (floored)
   real*4    :: d_up                            ! cell-mean depth of the upwind cell
   real*4    :: zb_face                         ! bed at the face (subgrid uv zmin)
   real*4    :: tol_edge                        ! tolerance of the face bed above the lowest pixel of the upwind cell
   real*4    :: r_slope                         ! bed drop over the face stencil / wet depth
   real*4    :: hu_flux                         ! (blended) depth for the flux and time step
   real*4    :: hu_fric                         ! (blended) depth for the friction
   real*4    :: dzdx_eff                        ! driving surface gradient
   !
   integer   :: iside                          ! face counter for the z_wface diagnostic
   integer   :: iface
   integer   :: jface
   real*4    :: wq_sum
   real*4    :: q_sum
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
   !$omp           fcoriouv,gnavg2,iwet,zsu,dzuv,iuv,facint,fwmax,zmax,zmin,dqxudx,dqyudy,un,up,vp,vn,umax, &
   !$omp           dnminv,qu,qd,hwet,phi,adv,mdrv,phiz,fac,hu43,y_cbrt,i_cbrt,min_dt_ip,zs2w,zs1e,dnm,dnmu,zrec,zbup,ipw,ipe,zbnm,zbnmu, &
   !$omp           iup,idn,w,s_s,hwet_up,d_up,zb_face,tol_edge,r_slope,hu_flux,hu_fric,dzdx_eff ) &
   !$omp reduction ( min : min_dt  )
   !$omp do schedule ( dynamic, 256 )
   !$acc parallel, present( kcuv, kfuv, zs, q, q0, uv, uv0, zsderv, z_wetfrac, &
   !$acc                    uv_flags_iref, uv_flags_type, uv_flags_dir, mask_adv, &
   !$acc                    subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_havg, subgrid_uv_nrep, subgrid_uv_pwet, &
   !$acc                    subgrid_uv_havg_zmax, subgrid_uv_nrep_zmax, subgrid_uv_fnfit, subgrid_uv_navg_w, &
   !$acc                    uv_index_z_nm, uv_index_z_nmu, uv_index_u_nmd, uv_index_u_nmu, uv_index_u_ndm, uv_index_u_num, &
   !$acc                    uv_index_v_ndm, uv_index_v_ndmu, uv_index_v_nm, uv_index_v_nmu, cuv_index_uv, cuv_index_uv1, cuv_index_uv2, &
   !$acc                    zb, tauwu, tauwv, patm, fwuv, gn2uv, dxminv, dxrinv, dyrinv, dxm2inv, dxr2inv, dyr2inv, &
   !$acc                    z_hwet, subgrid_z_dzbdm, subgrid_z_dzbdn, subgrid_z_zmin, w_uv, iup_uv, &
   !$acc                    dxrinvc, dyrinvc, fcorio2d, nuvisc, z_volume, cell_area, cell_area_m2, z_flags_iref, gnapp2, timestep_analysis_required_timestep ) num_gangs( 1024 ) vector_length( 128 )
   !$acc loop, reduction( min : min_dt ), gang, vector
   do ip = 1, npuv
      !
      if (store_slope_regime) then
         w_uv(ip)   = 0.0
         iup_uv(ip) = 0
      endif
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
         ! Upwind cell of the uv point, from the sign of the previous-step velocity, else the
         ! higher water level; zsu and, for regular grids, zbup follow it
         !
         if (uv0(ip) > 1.0e-6) then
            iup = nm
            idn = nmu
         elseif (uv0(ip) < -1.0e-6) then
            iup = nmu
            idn = nm
         else
            if (zs(nm) >= zs(nmu)) then
               iup = nm
               idn = nmu
            else
               iup = nmu
               idn = nm
            endif
         endif
         !
         zsu = zs(iup)
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
            zbup = zb(iup)
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
            ! Level-driven vs slope-driven flow (slope_driven_flow, subgrid only). A face is
            ! slope-driven when a sheet leaves the upwind cell over its low edge into a downstream
            ! water surface below the face bed (free overfall or dry neighbour). The face weight w
            ! follows from the ratio of the bed drop over the face stencil along the face normal
            ! and the level-based wet depth of the upwind cell. With w > 0 the flux and friction
            ! depths and g*n^2 are blended with those of a sheet over the upwind cell (cell-mean
            ! depth, deep-water g*n^2), and the driving gradient is at least the bed slope.
            ! With w = 0 everything below is as without the option.
            !
            w        = 0.0
            hu_flux  = hu
            hu_fric  = hu
            dzdx_eff = dzdx
            !
            if (slope_driven_flow) then
               !
               ! Bed-plane slope of the upwind cell along the face normal (positive = bed rising
               ! in +m / +n, same sign convention as dzdx)
               !
               if (uv_flags_dir(ip) == 0) then
                  s_s = subgrid_z_dzbdm(iup)
               else
                  s_s = subgrid_z_dzbdn(iup)
               endif
               !
               hwet_up  = max(z_hwet(iup), slope_driven_hmin)
               d_up     = max(real(zs(iup) - zb(iup), 4), 0.0)
               zb_face  = subgrid_uv_zmin(ip)
               tol_edge = max(5.0 * slope_driven_hmin, 0.25 * hwet_up)
               !
               ! Face at the low edge of the upwind cell, downstream surface below the face bed,
               ! and a bed slope along the face normal (nested, no short-circuit in Fortran)
               !
               if (zb_face - subgrid_z_zmin(iup) <= tol_edge) then
                  !
                  if (zs(idn) < zb_face) then
                     !
                     if (abs(s_s) > 1.0e-6) then
                        !
                        r_slope = abs(s_s) / (dxuvinv * hwet_up)
                        w = r_slope**2 / (r_slope**2 + slope_driven_ratio0**2)
                        !
                     endif
                     !
                  endif
                  !
               endif
               !
               if (w > 1.0e-6) then
                  !
                  hu_flux = (1.0 - w) * hu + w * d_up
                  hu_fric = (1.0 - w) * hu + w * max(d_up, slope_driven_hmin)
                  gnavg2  = (1.0 - w) * gnavg2 + w * subgrid_uv_navg_w(ip)
                  !
                  ! The downstream surface cannot exert hydrostatic back-pressure: the driving
                  ! gradient is at least the bed slope in the downslope direction
                  !
                  if (dzdx * s_s > 0.0) then
                     dzdx_eff = sign(max(abs(dzdx), abs(s_s)), dzdx)
                  endif
                  !
               else
                  !
                  w = 0.0
                  !
               endif
               !
               ! Store the face weight and the upwind cell for the diagnostic z_wface
               !
               if (store_slope_regime) then
                  w_uv(ip)   = w
                  iup_uv(ip) = iup
               endif
               !
            endif
            !
            ! Compute wet average depth hwet (used in wind and wave forcing)
            !
            hwet = hu_flux / phi
            !
            ! Velocity form: build frc directly as an acceleration [m/s^2]. Forces that scale
            ! with depth (pressure, viscosity, Coriolis, atm) are written WITHOUT hu -- their hu
            ! would only be divided out again. Only the surface stresses (wind, waves) keep a /hu.
            !
            frc = - g * dzdx_eff
            !
            if (advection) then
               !
!               if (mask_adv(ip) == 1) then
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
                  ! Clamp the advective speeds to the local characteristic speed
                  ! sqrt(g*hu) + |u|, the bound the time-step limiter guarantees to
                  ! resolve. On subgrid grids the cell-mean depth D used in dnminv can
                  ! be far smaller than the conveyance depth hu, so q/D can otherwise
                  ! exceed the advective CFL and blow up (e.g. steep valleys).
                  !
                  umax = sqrt(g * hu) + abs(uv0(ip))
                  !
                  up  = min(max(qd * dnminv, 0.0),  umax)   ! U_p  (advective speed, from west)
                  un  = max(min(qu * dnminv, 0.0), -umax)   ! U_n  (advective speed, from east)
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
                  vp = min(max( 0.5 * (q0(uv_index_v_ndm(ip)) + q0(uv_index_v_ndmu(ip))) * dnminv, 0.0 ),  umax)
                  vn = max(min( 0.5 * (q0(uv_index_v_nm(ip))  + q0(uv_index_v_nmu(ip)))  * dnminv, 0.0 ), -umax)
                  !
                  dqyudy = ( vp * (uu_nm - uu_ndm) + vn * (uu_num - uu_nm) ) * dyuvinv
                  !
                  adv = - phi * (dqxudx + dqyudy)        ! velocity tendency [m/s^2]
                  !
                  if (subgrid .and. advection_fade_power > 0) then
                     !
                     ! Fade advection on partly wet subgrid faces (phi = wet fraction of the face)
                     !
                     adv = adv * max(phi, 0.0)**advection_fade_power
                     !
                  endif
                  !
                  if (w > 0.0) then
                     !
                     ! No inertia in the kinematic (slope-driven) limit
                     !
                     adv = (1.0 - w) * adv
                     !
                  endif
                  !
!                  frc = frc + min(max(adv, -advlim), advlim)   ! add limited advective acceleration
                  frc = frc + adv   ! advective speeds are already clamped to umax above
                  !
!               endif
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
            i_cbrt = transfer(hu_fric, i_cbrt)
            i_cbrt = i_cbrt / 3 + 709921077
            y_cbrt = transfer(i_cbrt, y_cbrt)
            y_cbrt = y_cbrt - (y_cbrt * y_cbrt * y_cbrt - hu_fric) / (3.0 * y_cbrt * y_cbrt)
            hu43   = hu_fric * y_cbrt
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
               ! Wiggle suppression for cells in which only a few subgrid pixels are wet (large dzs/dV).
               ! There the local gravity-wave speed is sqrt(g * hu / phiz) with phiz the wet-area
               ! fraction of the cell, which exceeds the speed sqrt(g * hu) the time step was set for.
               ! Under-relax the velocity increment with fac = phiz / alfa**2 (floored at wiggle_facmin).
               ! Steady flow (uv = uv0) is not affected, 2*dt modes are damped. Fully wet cells give fac = 1.
               ! With wiggle_detect, only relax when the level accelerations of the two cells diverge
               ! by more than wiggle_threshold (anti-phase sloshing across this uv point).
               !
               phiz = min(z_wetfrac(nm), z_wetfrac(nmu))
               !
               fac  = min(max(phiz / alfa**2, wiggle_facmin), 1.0)
               !
               if (w > 0.0) then
                  !
                  ! Relaxation fades out in slope-driven faces
                  !
                  fac = 1.0 - (1.0 - w) * (1.0 - fac)
                  !
               endif
               !
               if (fac < 1.0) then
                  !
                  if (wiggle_detect) then
                     !
                     mdrv = abs(zsderv(nm) - zsderv(nmu)) - wiggle_threshold
                     !
                     if (mdrv > 0.0) then
                        uv(ip) = uv0(ip) + fac * (uv(ip) - uv0(ip))
                     endif
                     !
                  else
                     !
                     uv(ip) = uv0(ip) + fac * (uv(ip) - uv0(ip))
                     !
                  endif
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
            q(ip) = uv(ip) * hu_flux
            !
            kfuv(ip) = 1
            !
            ! Determine minimum time step (alpha is added later on in sfincs_lib.f90) of all uv points
            ! Use maximum of sqrt(gh) and current velocity
            !
            min_dt_ip = 1.0 / ( max(sqrt(g * hu_flux), abs(uv(ip)) ) * dxuvinv)
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
   if (store_slope_regime) then
      !
      ! Diagnostic z_wface (map output only, store_slope_regime): face regime weight used above,
      ! flux-weighted over the outflow faces of each cell, i.e. the faces for which the cell is the
      ! upwind cell (combined quadtree uv points via their two sub-faces). A cell without outflow
      ! faces gets 0. Cell-parallel, no atomics.
      !
      !$omp parallel &
      !$omp private ( nm,iside,iface,jface,icuv,wq_sum,q_sum )
      !$omp do schedule ( dynamic, 256 )
      !$acc parallel present( z_wface, w_uv, iup_uv, q, z_index_uv_md, z_index_uv_mu, z_index_uv_nd, z_index_uv_nu, &
      !$acc                   cuv_index_uv1, cuv_index_uv2 )
      !$acc loop gang vector
      do nm = 1, np
         !
         wq_sum = 0.0
         q_sum  = 0.0
         !
         do iside = 1, 4
            !
            if (iside == 1) then
               iface = z_index_uv_md(nm)
            elseif (iside == 2) then
               iface = z_index_uv_mu(nm)
            elseif (iside == 3) then
               iface = z_index_uv_nd(nm)
            else
               iface = z_index_uv_nu(nm)
            endif
            !
            if (iface <= npuv) then
               !
               if (iup_uv(iface) == nm) then
                  wq_sum = wq_sum + w_uv(iface) * abs(q(iface))
                  q_sum  = q_sum + abs(q(iface))
               endif
               !
            elseif (iface <= npuv + ncuv) then
               !
               icuv  = iface - npuv
               jface = cuv_index_uv1(icuv)
               !
               if (jface <= npuv .and. iup_uv(min(jface, npuv)) == nm) then
                  wq_sum = wq_sum + w_uv(jface) * abs(q(jface))
                  q_sum  = q_sum + abs(q(jface))
               endif
               !
               jface = cuv_index_uv2(icuv)
               !
               if (jface <= npuv .and. iup_uv(min(jface, npuv)) == nm) then
                  wq_sum = wq_sum + w_uv(jface) * abs(q(jface))
                  q_sum  = q_sum + abs(q(jface))
               endif
               !
            endif
            !
         enddo
         !
         z_wface(nm) = wq_sum / max(q_sum, 1.0e-12)
         !
      enddo
      !$omp end do
      !$omp end parallel
      !$acc end parallel
      !
   endif
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine
   !
end module
