   module sfincs_momentum
   
   contains

   subroutine compute_fluxes(dt, min_dt, tloop)
   !
   ! Computes fluxes over subgrid u and v points
   !
   use sfincs_data
   !
   implicit none
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

   integer   :: idir
   integer   :: iref
   integer   :: itype
   integer   :: iuv
   integer   :: ind
   !
   real*4    :: huv
   real*4    :: dxuvinv
   real*4    :: dxuv2inv
   real*4    :: dyuvinv
   real*4    :: dyuv2inv
   real*4    :: adv
   real*4    :: fcoriouv
   real*4    :: frc
   real*4    :: min_dt
   real*4    :: gammax
   real*4    :: facmax
   real*4    :: wsumax
   real*4    :: tp
   !
   real*4    :: qx0_nm
   real*4    :: qx0_nmd
   real*4    :: qx0_nmu
   real*4    :: qx0_ndm
   real*4    :: qx0_num
   real*4    :: qsm
   !   
   real*4    :: qy0_nm
   real*4    :: qy0_nmu
   real*4    :: qy0_ndm
   real*4    :: qy0_ndmu
   !
   real*4    :: u0_nm
   real*4    :: u0_nmd
   real*4    :: u0_nmu
   real*4    :: u0_ndm
   real*4    :: u0_num
   real*4    :: zsuv
   real*4    :: dzuv
   real*4    :: facint
   real*4    :: gnavg2
   real*4    :: v0_nm
   real*4    :: v0_nmu
   real*4    :: v0_ndm
   real*4    :: v0_ndmu
   real*4    :: fwmax
   real*4    :: zmax
   real*4    :: zmin
   real*4    :: one_minus_facint 
   !
   real*4, parameter :: expo = 1.0/3.0
!   integer, parameter :: expo = 1
   !
   logical   :: iadv, ivis, icorio, iok
   !
   call system_clock(count0, count_rate, count_max)
   !
   min_dt = dtmax
   !
   !$acc update device(min_dt), async(1)
   !
!   q0(:)  = q(:)
!   uv0(:) = uv(:)
   !
   !$omp parallel &
   !$omp private ( ip ) &
   !$omp shared ( q,q0,uv,uv0 )
   !$omp do
   !$acc kernels, present(q, q0, uv, uv0), async(1)
   !$acc loop independent, private(nm)
   do ip = 1, npuv
      q0(ip)  = q(ip)
      uv0(ip) = uv(ip)
   enddo
   !$acc end kernels
   !$omp end do
   !$omp end parallel
   !
   !$omp parallel &
   !$omp private ( ip,huv,qsm,qx0_nm,nm,nmu,frc,adv,idir,itype,iadv,ivis,icorio,iref,dxuvinv,dxuv2inv,dyuvinv,dyuv2inv,qx0_nmd,qx0_nmu,qx0_ndm,qx0_num,u0_nm,u0_nmd,u0_nmu,u0_num,u0_ndm,v0_ndm,v0_ndmu,v0_nm,v0_nmu,fcoriouv,gnavg2,iok,zsuv,dzuv,iuv,facint,fwmax,zmax,zmin,one_minus_facint ) &
   !$omp shared ( q,uv ) &
   !$omp reduction ( min : min_dt )
   !$omp do schedule ( dynamic, 256 )
   !$acc kernels, present( kcuv, zs, q, q0, uv, uv0, min_dt, &
   !$acc                   uv_flags_iref, uv_flags_type, uv_flags_vis, uv_flags_adv, uv_flags_dir, &
   !$acc                   subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_hrep, subgrid_uv_navg, subgrid_uv_hrep_zmax, subgrid_uv_navg_zmax, &
   !$acc                   uv_index_z_nm, uv_index_z_nmu, uv_index_u_nmd, uv_index_u_nmu, uv_index_u_ndm, uv_index_u_num, &
   !$acc                   uv_index_v_ndm, uv_index_v_ndmu, uv_index_v_nm, uv_index_v_nmu, &
   !$acc                   zb, zbuv, zbuvmx, tauwu, tauwv, patm, fwuv, gn2uv, dxminv, dxrinv, dyrinv, dxm2inv, dxr2inv, dyr2inv, dxrinvc, fcorio2d ), async(1)
   !$acc loop independent, reduction( min : min_dt )
   do ip = 1, npuv
      !
      if (kcuv(ip)==1) then
         !
         ! Indices of surrounding water level points
         !
         nm  = uv_index_z_nm(ip)
         nmu = uv_index_z_nmu(ip)
         !
         iok  = .false.
         !
         zsuv = max(zs(nm), zs(nmu)) ! water level at u point
         !
         if (subgrid) then
            !
            zmin = subgrid_uv_zmin(ip)
            zmax = subgrid_uv_zmax(ip)
            !            
            if (zsuv>zmin + huthresh) then
               iok = .true.
            endif   
            !
         else
            !            
            if (zsuv>zbuvmx(ip)) then
               iok = .true.
            endif   
            !            
         endif   
         !
         if (iok) then
            !
            if (use_quadtree) then
               iref  = uv_flags_iref(ip) ! refinement level
               itype = uv_flags_type(ip) ! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            else
               iref  = 1
               itype = 0
            endif
            !
            idir  = uv_flags_dir(ip) ! 0 is u, 1 is v
            !
            ! Check if viscosity term is computed for this point
            !
            if (viscosity) then
               if (uv_flags_vis(ip)==1) then
                  ivis = .true.
               else
                  ivis = .false.
               endif   
            else
               ivis = .false.
            endif                  
            !
            ! Check if advection term is computed for this point
            !
            if (advection) then
               if (uv_flags_adv(ip)==1) then
                  iadv = .true.
               else
                  iadv = .false.
               endif   
            else
               iadv = .false.
            endif
            !
            ! Check if Coriolis term is computed for this point
            !
            if (coriolis) then
               !
               if (uv_flags_adv(ip)==1) then
                  icorio = .true.
               else
                  icorio = .false.
               endif   
               !
            else
               !
               icorio = .false.
               !
            endif
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
                     dyuvinv  = dxminv(ip)
                     dxuvinv  = dyrinv(iref)
                     dyuv2inv = dxm2inv(ip)
                     dxuv2inv = dyr2inv(iref)
                     !
                  endif   
               else   
                  !
                  ! Fine to coarse or coarse to fine
                  !
                  if (idir==0) then
                     dxuvinv = 1.0 / (3*(1.0/dxminv(ip))/2)
                  else   
                     dxuvinv = 1.0 / (3*(1.0/dyrinv(iref))/2)
                  endif   
                  !
                  ! Should not need the other values ...
                  ! dxuv2inv = 0.0 ! Should not need this one as viscosity term is not computed for this uv point
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
                  dxuvinv  = dxrinvc(iref)
                  !
                  ! Should not need the other values ...
                  ! dxuv2inv = 0.0 ! Should not need this one as viscosity term is not computed for this uv point
                  !
               endif
               !
               fcoriouv = fcorio
               !
            endif
            !
            ! Get fluxes and velocities from previous time step
            !
            qx0_nm = q0(ip)
            !
            if (iadv) then
               !
               ! First term
               !
               qx0_nmd  = q0(uv_index_u_nmd(ip))
               qx0_nmu  = q0(uv_index_u_nmu(ip))
               u0_nm    = uv0(ip)
               u0_nmd   = uv0(uv_index_u_nmd(ip))
               u0_nmu   = uv0(uv_index_u_nmu(ip))
               !
               ! Second term
               !
               qx0_ndm  = q0(uv_index_u_ndm(ip))
               qx0_num  = q0(uv_index_u_num(ip))
               v0_ndm   = uv0(uv_index_v_ndm(ip))
               v0_ndmu  = uv0(uv_index_v_ndmu(ip))
               v0_nm    = uv0(uv_index_v_nm(ip))
               v0_nmu   = uv0(uv_index_v_nmu(ip))
               !
            endif
            !
            if (ivis) then
               !
               if (.not. iadv) then
                  u0_nm  = uv0(ip)
                  u0_nmd = uv0(uv_index_u_nmd(ip))
                  u0_nmu = uv0(uv_index_u_nmu(ip))
               endif
               !
               u0_ndm  = uv0(uv_index_u_ndm(ip))
               u0_num  = uv0(uv_index_u_num(ip))
               !
            endif
            !
            if (icorio .and. .not.iadv) then
               !
               v0_nm    = uv0(uv_index_v_nm(ip))
               v0_ndm   = uv0(uv_index_v_ndm(ip))
               v0_ndmu  = uv0(uv_index_v_ndmu(ip))
               v0_nmu   = uv0(uv_index_v_nmu(ip))
               !
            endif
            !
            if (theta<0.9999) then ! for backward compatibility
               qx0_nmd  = q0(uv_index_u_nmd(ip))
               qx0_nmu  = q0(uv_index_u_nmu(ip))
            endif
            !
            ! Compute water depth at uv point
            !
            if (subgrid) then
               !
               if (zsuv>zmax - 1.0e-4) then
                  !
                  ! Entire cell is wet, no interpolation from table needed
                  !
                  huv    = subgrid_uv_hrep_zmax(ip) + zsuv
                  gnavg2 = subgrid_uv_navg_zmax(ip)
                  !
               else
                  !
                  ! Interpolation required
                  !
                  dzuv   = (subgrid_uv_zmax(ip) - subgrid_uv_zmin(ip)) / (subgrid_nbins - 1)
                  iuv    = int((zsuv - subgrid_uv_zmin(ip))/dzuv) + 1
                  facint = (zsuv - (subgrid_uv_zmin(ip) + (iuv - 1)*dzuv) ) / dzuv
                  huv    = subgrid_uv_hrep(iuv, ip) + (subgrid_uv_hrep(iuv + 1, ip) - subgrid_uv_hrep(iuv, ip))*facint
                  gnavg2 = subgrid_uv_navg(iuv, ip) + (subgrid_uv_navg(iuv + 1, ip) - subgrid_uv_navg(iuv, ip))*facint
                  !
               endif
               !
               huv    = max(huv, huthresh)
               !
            else
               !
               huv    = max(zsuv - zbuv(ip), huthresh)
               gnavg2 = gn2uv(ip)
               !
            endif
            !
            ! Determine minimum time step (alpha is added later on in sfincs_lib.f90) of all uv points
            !
            min_dt = min(min_dt, 1.0/(sqrt(g*huv)*dxuvinv))
            !
            ! FORCING TERMS
            !
            ! Pressure term
            !
            frc = -g*huv*(zs(nmu) - zs(nm))*dxuvinv
            !
            ! Advection term
            !
            if (iadv) then
               !
               ! 1D upwind advection slightly more robust than central scheme
               !
               if (qx0_nm>1.0e-6) then
                  !
                  adv = ( qx0_nm*u0_nm - qx0_nmd*u0_nmd ) * dxuvinv                  
                  !
               elseif (qx0_nm<-1.0e-6) then
                  !
                  adv = ( qx0_nmu*u0_nmu - qx0_nm*u0_nm ) * dxuvinv
                  !
               else
                  !
                  adv = 0.0
                  !
               endif
               !
               if ( (v0_nmu + v0_nm > 0.0) .and. (v0_ndm + v0_ndmu > 0.0) ) then
                  !
                  adv = adv + ( (v0_ndm + v0_ndmu) * (qx0_nm - qx0_ndm) ) * dyuvinv / 2
                  !
               elseif ( (v0_nmu + v0_nm < 0.0) .and. (v0_ndm + v0_ndmu < 0.0) ) then
                  !
                  adv = adv + ( (v0_nm + v0_nmu) * (qx0_num - qx0_nm) ) * dyuvinv / 2
                  !
               endif   
               !   
               ! Let's try without the advection limiter  
               ! adv = min(max(adv, -advlim), advlim)
               !
               frc = frc - adv
               !   
            endif   
            !
            ! Viscosity term
            !
            if (ivis) then
               !
               frc = frc + nuvisc * huv * ( (u0_nmu - 2*u0_nm + u0_nmd )*dxuv2inv + (u0_num - 2*u0_nm + u0_ndm )*dyuv2inv )
               !
            endif
            !
            ! Coriolis term
            !
            if (icorio) then
               !
               if (idir==0) then
                  !
                  frc = frc + 0.25*fcoriouv*huv*(v0_nm + v0_ndm + v0_nmu + v0_ndmu)! U
                  !
               else
                  !
                  frc = frc - 0.25*fcoriouv*huv*(v0_nm + v0_ndm + v0_nmu + v0_ndmu)! U
                  !
               endif
               !
            endif            
            !
            ! Wind forcing
            !
            if (wind) then
               !
               if (huv>0.25) then
                  !
                  if (idir==0) then
                     !
                     frc = frc + tauwu(nm)
                     !
                  else
                     !
                     frc = frc + tauwv(nm)
                     !
                  endif   
                  !
               else
                  !
                  ! Reduce wind drag at small water depths
                  !
                  if (idir==0) then
                     !
                     frc = frc + tauwu(nm)*huv*4
                     !
                  else
                     !
                     frc = frc + tauwv(nm)*huv*4
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
               frc = frc - huv*(patm(nmu) - patm(nm))*dxuvinv/rhow
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
               ! fmax = facmax*huv*sqrt(huv)/tp/rhow (we already divided by rhow in sfincs_snapwave)
               !
               fwmax = 0.8*huv*sqrt(huv)/15
               !
               frc = frc + sign(min(abs(fwuv(ip)), fwmax), fwuv(ip))
               !
            endif   
            !
            ! Apply some smoothing if theta < 1.0 (not recommended anymore!)
            !
            if (theta<0.9999) then
               qsm = theta*qx0_nm + 0.5*(1.0 - theta)*(qx0_nmu + qx0_nmd)             
            else
               qsm = qx0_nm
            endif            
            !
            ! Compute new flux for this uv point (Bates et al., 2010)
            ! 
            qx0_nm = (qsm + frc*dt) / (1.0 + gnavg2*dt*abs(qx0_nm)/(huv**2*huv**expo))
            !
!            if (zs(nm)<=zb(nm)) then
!               !
!               qx0_nm = min(qx0_nm, 0.0)
!               !
!            endif
!            !
!            if (zs(nmu)<=zb(nmu)) then
!               !
!               qx0_nm = max(qx0_nm, 0.0)
!               !
!            endif
            !
            q(ip) = qx0_nm
            !
            ! Limit velocities (this does not change fluxes, but may prevent advection term to explode in the next time step)
            !
            uv(ip) = max(min(qx0_nm/huv, 4.0), -4.0)
            !
         else
            !
            q(ip)  = 0.0
            uv(ip) = 0.0
            !
         endif
         !
      endif
   enddo   
   !$omp end do
   !$omp end parallel
   !$acc end kernels
   !
   !$acc update host(min_dt), async(1)
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !         
   end subroutine

   
   subroutine compute_fluxes_simple(dt, min_dt, tloop)
   !
   ! Computes fluxes over subgrid u and v points
   !
   use sfincs_data
   !
   implicit none
   !
   integer  :: count0
   integer  :: count1
   integer  :: count_rate
   integer  :: count_max
   real     :: tloop
   !
   real*4    :: dt
   !
   integer   :: ip
   integer   :: nm
   integer   :: nmu
   integer   :: idir
   integer   :: iref
   integer   :: itype
   integer   :: iuv
   integer   :: ind
   !
   real*4    :: huv
   real*4    :: dxuvinv
   real*4    :: dxuv2inv
   real*4    :: dyuvinv
   real*4    :: dyuv2inv
   real*4    :: adv
   real*4    :: fcoriouv
   real*4    :: frc
   real*4    :: min_dt
   real*4    :: gammax
   real*4    :: facmax
   real*4    :: wsumax
   real*4    :: tp
   !
   real*4    :: qx0_nm
   real*4    :: qx0_nmd
   real*4    :: qx0_nmu
   real*4    :: qy0_nm
   real*4    :: qy0_nmu
   real*4    :: qy0_ndm
   real*4    :: qy0_ndmu
   real*4    :: u0_nm
   real*4    :: u0_nmd
   real*4    :: u0_nmu
   real*4    :: u0_ndm
   real*4    :: u0_num
   real*4    :: zsuv
   real*4    :: dzuv
   real*4    :: facint
   real*4    :: gnavg2
   real*4    :: v0_nm
   real*4    :: v0_nmu
   real*4    :: v0_ndm
   real*4    :: v0_ndmu
   real*4, dimension(1:subgrid_nbins) :: ahrep
   real*4, dimension(1:subgrid_nbins) :: anavg
   real*4    :: zmin
   real*4    :: zmax
   real*4    :: iuvr
   real*4    :: one_minus_facint
   !
   real*4, parameter :: expo=1.0/3.0
   !
   logical   :: iadv, ivis, icorio, iok
   !
   call system_clock(count0, count_rate, count_max)
   !
   min_dt = 1.0e6
   !
   !$omp parallel &
   !$omp private ( ip,huv,qx0_nm,nm,nmu,frc,dxuvinv,gnavg2,iok,zsuv,dzuv,iuv,facint,ahrep,anavg,zmin,zmax,one_minus_facint,iuvr,ind ) &
   !$omp shared ( kcuv,q,zs,subgrid_uv_zmin,subgrid_uv_zmax,subgrid_uv_hrep,subgrid_uv_navg,uv_index_z_nm,uv_index_z_nmu,subgrid_uv_hrep_zmax,subgrid_uv_navg_zmax ) &
   !$omp reduction ( min : min_dt )
   !$omp do
   do ip = 1, npuv
      !
      if (kcuv(ip)==1) then
         !
         ! Indices of surrounding water level points
         !
         nm  = uv_index_z_nm(ip)
         nmu = uv_index_z_nmu(ip)
         !
         iok  = .false.
         !
         zsuv = max(zs(nm), zs(nmu)) ! water level at u point
         !
         if (subgrid) then
            !            
            if (zsuv>subgrid_uv_zmin(ip) + huthresh) then
               iok = .true.
            endif   
            !
         else
            !            
            if (zsuv>zbuvmx(ip)) then
               iok = .true.
            endif   
            !            
         endif   
         !
         if (iok) then
            !
            dxuvinv = 1.0/dx
            !
            ! Compute water depth at uv point
            !
            if (subgrid) then
               !
               zmin = subgrid_uv_zmin(ip)
               zmax = subgrid_uv_zmax(ip)
               !            
               if (zsuv>zmax - 0.001) then
                  !
                  ! Entire cell is wet, no interpolation from table needed
                  !
                  huv    = subgrid_uv_hrep_zmax(ip) + zsuv
                  gnavg2 = subgrid_uv_navg_zmax(ip)
                  !
               else
                  !
                  ! Interpolation required
                  !
                  dzuv   = (subgrid_uv_zmax(ip) - subgrid_uv_zmin(ip)) / (subgrid_nbins - 1)
                  iuv    = int((zsuv - subgrid_uv_zmin(ip))/dzuv) + 1
                  facint = (zsuv - (subgrid_uv_zmin(ip) + (iuv - 1)*dzuv) ) / dzuv
                  huv    = subgrid_uv_hrep(iuv, ip) + (subgrid_uv_hrep(iuv + 1, ip) - subgrid_uv_hrep(iuv, ip))*facint
                  gnavg2 = subgrid_uv_navg(iuv, ip) + (subgrid_uv_navg(iuv + 1, ip) - subgrid_uv_navg(iuv, ip))*facint
                  !
               endif
               !
            else
               !
               huv  = max(zsuv - zbuv(ip), huthresh)
               gnavg2 = gn2uv(ip)
               !
            endif
            !
            ! Determine minimum time step (alpha is added later on in sfincs_lib.f90) of all uv points
            !   
            qx0_nm = q(ip)
            !
            ! FORCING TERMS
            !
            ! Pressure term
            !
            frc = -g*huv*(zs(nmu) - zs(nm))*dxuvinv
            !
            ! Compute new flux for this uv point (Bates et al., 2010)
            !
            qx0_nm  = (qx0_nm + frc*dt ) / (1.0 + gnavg2*dt*abs(qx0_nm)/(huv**2*huv**expo))
            !
            q(ip) = qx0_nm
            !
            min_dt = min(min_dt, 1.0/(sqrt(g*huv)*dxuvinv))
            !
         else
            !
            q(ip) = 0.0
            !
         endif
         !
      endif
   enddo   
   !$omp end do
   !$omp end parallel
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !         
   end subroutine
   
end module
