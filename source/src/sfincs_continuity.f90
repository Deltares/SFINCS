module sfincs_continuity
   !
   ! Water-level / volume update stage of the SFINCS time step. Runs after
   ! sfincs_momentum has produced the face fluxes q on each cell edge and
   ! is responsible for closing the volume balance on every active cell.
   !
   ! See the breakdown at the top of update_continuity for exactly
   ! which terms are accumulated into qsrc, which operate on qinfmap, and
   ! which come from the hydrodynamic fluxes q already computed upstream.
   !
   ! Data flow per step:
   !   input  : q(nuv), qext(np) (optional BMI), river/structure state,
   !            qinfmap, storage_volume (subgrid), zs/z_volume at t
   !   output : zs (all paths) and z_volume (subgrid path) advanced to
   !            t+dt; optional zsmax, vmax, qmax, twet accumulators
   !
   ! Subroutines:
   !
   !   update_continuity(t, dt)
   !     Main per-timestep entry. Orchestrates river discharges, drainage
   !     structures, optional BMI qext, infiltration, and dispatches the
   !     water-level update. Called from sfincs_lib (main time-stepping
   !     loop).
   !
   !   compute_water_levels_regular(dt, t)
   !     Non-subgrid (bathtub / simple bathy) water-level update. Called
   !     from update_continuity.
   !
   !   compute_water_levels_subgrid(dt, t)
   !     Subgrid-tables water-level update with storage-volume
   !     bookkeeping. Called from update_continuity.
   !
   !   compute_store_variables(dt)
   !     Update optional per-cell vmax / qmax / twet diagnostics. Called
   !     from update_continuity only when any of store_maximum_velocity,
   !     store_maximum_flux or store_twet is enabled.
   !
contains
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine update_continuity(t, dt)
      !
      ! Unified continuity update: orchestrates all water balance terms
      ! for one time step. Modifies qsrc in place (zeroed and
      ! re-accumulated), advances zs (and z_volume on the subgrid path),
      ! and optionally updates the store_* running maxima.
      !
      ! Called from: sfincs_lib (main time-stepping loop).
      !
      ! Sources and sinks (all accumulated into qsrc, in m3/s):
      !    1. River discharges (+/-)               => update_discharges (zeros and accumulates qsrc)
      !    2. Drainage structures (+/-)            => update_src_structures (adds to qsrc)
      !    3. Precipitation (+)                    => update_meteo_forcing (precip * cell area)
      !    4. Infiltration rate field qinfmap (-)  => update_infiltration_map (infiltration * cell area)
      !                                              (flavors: con, c2d, cna, cnb, gai, hor, bkt)
      !    5. Urban drainage                       => update_urban_drainage (adds to qsrc)
      !    6. External source/sink qext (+/-)      => added to qsrc here (BMI coupling)      
      !
      ! Hydrodynamic fluxes q                      => computed in sfincs_momentum
      !
      ! compute_water_levels_{regular,subgrid} then updates zs/z_volume using:
      !    - qsrc * dt                             => point source/sink contribution
      !    - div(q) * dt                           => horizontal flux divergence
      !    - storage volume                        => absorbs excess volume (subgrid only)
      !
      use sfincs_data
      use sfincs_timers
      use sfincs_infiltration
      use sfincs_discharges
      use sfincs_src_structures
      use sfincs_urban_drainage
      !
      implicit none
      !
      real*8           :: t
      real*4           :: dt
      !
      integer          :: nm
      !
      ! 1. River discharges => update_discharges (adds to qsrc)
      !
      call update_discharges(t, dt)
      !
      ! 2. Drainage structures (pumps/gates/culverts/...) => update_src_structures (adds to qsrc)
      !
      call update_src_structures(t, dt)
      !
      ! 3. Precipitation => update_meteo_forcing (adds to qsrc)
      !
      ! 4. Compute infiltration rates => qinfmap(adds to qsrc)
      !
      if (infiltration) then
         !
         call update_infiltration_map(dt)
         !
      endif
      !
      ! 5. Urban drainage => update_urban_drainage (adds to qsrc)
      !
      if (urban_drainage) then
         !
         call update_urban_drainage(t, dt)
         !
      endif
      !
      ! 6. External source/sink (+/-) => add qext to qsrc (set via BMI coupling)
      !
      if (use_qext) then
         !
         !$omp parallel &
         !$omp private ( nm )
         !$omp do
         !$acc loop gang vector
         do nm = 1, np
            !
            qsrc(nm) = qsrc(nm) + qext(nm)
            !
         enddo
         !$acc end loop
         !$omp end parallel
         !
      endif
      !
      ! Update water levels: applies qsrc * dt and flux divergence to zs/z_volume
      !
      call timer_start('Continuity')
      !
      if (subgrid) then
         !
         call compute_water_levels_subgrid(dt, t)
         !
      else
         !
         call compute_water_levels_regular(dt, t)
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
      call timer_stop('Continuity')
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine compute_water_levels_regular(dt, t)
      !
      ! Advance zs(np) by dt on the non-subgrid (bathtub / simple bathy)
      ! path. Applies cell-wise qsrc contributions and the horizontal
      ! flux divergence, handles the optional wavemaker cells (kcs == 4),
      ! updates the snapwave-filtered water level zsm, and accumulates
      ! zsmax / t_zsmax when requested.
      !
      ! Called from: update_continuity (when subgrid is false).
      !
      use sfincs_data
      !
      implicit none
      !
      real*4           :: dt
      real*8           :: t
      !
      integer          :: nm
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
      !
      if (snapwave) then ! need to compute filtered water levels for snapwave
         !
         factime = min(dt / wavemaker_filter_time, 1.0)
         !
      endif
      !
      !$acc parallel present( kcs, zs, zb, prcp, q, qext, qinfmap, qdrain_rate, zsmax, zsm, maxzsm, &
      !$acc                   z_flags_iref, uv_flags_iref, &
      !$acc                   z_index_uv_md, z_index_uv_nd, z_index_uv_mu, z_index_uv_nu, &
      !$acc                   dxm, dxrm, dyrm, dxminv, dxrinv, dyrinv, cell_area_m2, cell_area,  &
      !$acc                   qsrc, &
      !$acc                   z_index_wavemaker, wavemaker_uvmean, wavemaker_nmd, wavemaker_nmu, wavemaker_ndm, wavemaker_num )
      !
      !$omp parallel &
      !$omp private ( nm, nmd, nmu, ndm, num, qnmd, qnmu, qndm, qnum, iwm )
      !$omp do schedule ( dynamic, 256 )
      !$acc loop gang vector
      do nm = 1, np
         !
         if (kcs(nm) == 1) then ! Regular point
            !
            ! Apply cell-wise discharges qsrc (rivers, drainage structures, qext)
            !
            if (qsrc(nm) /= 0.0) then
               !
               if (crsgeo) then
                  zs(nm) = max(zs(nm) + qsrc(nm) * dt / cell_area_m2(nm), zb(nm))
               else
                  zs(nm) = max(zs(nm) + qsrc(nm) * dt / cell_area(z_flags_iref(nm)), zb(nm))
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
         ! Reset qsrc to zero for the next time step
         !
         qsrc(nm) = 0.0
         !
      enddo
      !$omp end do
      !$omp end parallel
      !$acc end parallel
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine compute_water_levels_subgrid(dt,t)
      !
      ! Advance z_volume(np) and zs(np) by dt on the subgrid-tables path.
      ! Accumulates the cell volume change dvol from qsrc and the
      ! horizontal flux divergence, routes excess volume through
      ! storage_volume (when use_storage_volume is set), updates
      ! z_volume, and recovers the new water level via subgrid table
      ! interpolation. Also handles wavemaker cells (kcs == 4), the
      ! snapwave-filtered zsm, and zsmax / t_zsmax accumulation.
      !
      ! Called from: update_continuity (when subgrid is true).
      !
      use sfincs_data
      !
      implicit none
      !
      real*4           :: dt
      real*8           :: t
      !
      integer          :: nm
      !
      integer          :: iwm
      !
      integer          :: nmu
      integer          :: nmd
      integer          :: num
      integer          :: ndm
      !
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
      real*4           :: dv
      real*4           :: zs00
      real*4           :: zs11
      !
      if (wavemaker) then
         !
         factime = min(dt / wavemaker_filter_time, 1.0)
         !
      endif
      !
      !$omp parallel &
      !$omp private ( dvol, nmd, nmu, ndm, num, a, iuv, facint, dzvol, iwm, &
      !$omp           qnmd, qnmu, qndm, qnum, dv, zs00, zs11 )
      !$omp do schedule ( dynamic, 256 )
      !$acc parallel present( kcs, zs, zs0, zb, z_volume, zsmax, zsm, maxzsm, zsderv, &
      !$acc                   subgrid_z_zmin,  subgrid_z_zmax, subgrid_z_dep, subgrid_z_volmax, &
      !$acc                   prcp, q, qext, qinfmap, qdrain_rate, z_flags_iref, uv_flags_iref, &
      !$acc                   z_index_uv_md, z_index_uv_nd, z_index_uv_mu, z_index_uv_nu, &
      !$acc                   dxm, dxrm, dyrm, dxminv, dxrinv, dyrinv, cell_area_m2, cell_area, &
      !$acc                   z_index_wavemaker, wavemaker_uvmean, &
      !$acc                   wavemaker_nmd, wavemaker_nmu, wavemaker_ndm, wavemaker_num, storage_volume)
      !$acc loop gang vector
      do nm = 1, np
         !
         ! And now water level changes due to horizontal fluxes
         !
         dvol = 0.0
         !
         if (kcs(nm) == 1) then
            !
            ! Apply cell-wise discharges qsrc (rivers, drainage structures, qext)
            !
            if (qsrc(nm) /= 0.0) then
               !
               dvol = dvol + qsrc(nm) * dt
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
            ! C5. Storage volume
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
         ! Reset qsrc to zero for the next time step
         !
         qsrc(nm) = 0.0
         !
      enddo
      !$omp end do
      !$omp end parallel
      !
      !$acc end parallel
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine compute_store_variables(dt)
      !
      ! Update the optional per-cell running diagnostics vmax, qmax and
      ! twet from the edge-centred uv / q fields at the current step.
      ! Cell-centred vmax / qmax are reconstructed as the 2D magnitude
      ! of the mean of the four surrounding edge values. twet is the
      ! cumulative time a cell has been wet above twet_threshold. Kept
      ! in a separate routine to avoid the overhead when unused.
      !
      ! Called from: update_continuity (only when any of
      ! store_maximum_velocity, store_maximum_flux or store_twet is
      ! enabled).
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
   !
end module
