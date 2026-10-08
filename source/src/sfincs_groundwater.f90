module sfincs_groundwater
   !
   ! Groundwater table model (after Sanders et al. 2025, PRIMo). Each active cell
   ! carries a water table that rises by infiltration and, when gw_conductivity and
   ! gw_aquifer_thickness are given, drains by lateral flow to its neighbours and to open water. Switched on with
   ! groundwater = 1 underneath any infiltration method (inftype).
   !
   ! The aquifer has its own clock. Between groundwater steps the table does not
   ! change; the surface side only keeps books:
   !
   !    every hydrodynamic step (update_groundwater, called from update_continuity):
   !       gw_recharge          += qinfmap * dt           infiltrated depth since the last groundwater step
   !       gw_infiltration_cap   = (gw_space - gw_recharge) / dt     largest rate the method may use next step
   !       qsrc                 += gw_surface_exchange * A   exchange rate set at the last groundwater step
!
   !    every gw_dt seconds (default 60 s, groundwater_step):
   !       H += gw_recharge / Sy                           vertical, in one lump
   !       lateral substeps (lateral flow on)             Sy A dH = sum C_f (H_nb - H) dt
   !       H > ground: excess exfiltrates                 open-water cells receive / supply
   !       gw_surface_exchange = delivered depth / gw_dt  rate for the coming interval
   !       gw_space = (ground - H) * Sy, gw_recharge = 0
   !
   ! The infiltration methods limit their rate with gw_infiltration_cap before updating
   ! their own state, so saturation is exact within an interval and their bookkeeping
   ! agrees with the water table. That array is the only coupling with the infiltration
   ! module: the groundwater model also runs without any infiltration method, and then
   ! only exchanges water with the surface through lateral flow. State and parameters
   ! live in sfincs_data (gw_*).
!
   !   initialize_groundwater()
   !     Allocates the gw_* arrays, reads optional fields from inffile, applies the
   !     uniform gw_* overrides, classifies cells, seeds the water table (also from a
   !     restart), face conductances and the stable step. Called from sfincs_lib.
!
   !   update_groundwater(dt)
   !     Per hydrodynamic step bookkeeping (above); calls groundwater_step when the
   !     accumulated time reaches gw_dt. Called from update_continuity after the
   !     infiltration update.
!
   !   groundwater_step(dtgw)
   !     The physics: recharge, lateral flow, exfiltration, exchange rate, new space.
   !
   use sfincs_log
   use sfincs_error
   !
contains

   subroutine initialize_groundwater()
   !
   use sfincs_data
   use sfincs_ncinput
   !
   implicit none
   !
   integer       :: nm
   character*256 :: varname
   logical       :: have_depth, have_sy
   logical       :: have_level, wet_is_open_water
   integer       :: n_wet_open
   integer       :: n_floored
   real*4, dimension(:), allocatable :: gw_initial_level
   real*4, dimension(:), allocatable :: gw_level_ini
   real*4, dimension(:), allocatable :: gw_initial_depth
   real*4, dimension(:), allocatable :: gw_conductivity
   real*4, dimension(:), allocatable :: gw_aquifer_thickness
   real*4, dimension(:), allocatable :: csum
   logical       :: have_k, have_b
   integer       :: ip, nmu, iref, itype, idir
   real*4        :: tr1, tr2, trf, dist, width
   !
   if (inftype == 'bkt') then
      call stop_sfincs('Error ! The bucket model (inftype = bkt) is itself a storage model and cannot be combined with groundwater = 1 !', 1)
   endif
   !
   call write_log('Info    : turning on groundwater table model', 0)
   !
   allocate(gw_level(np))
   allocate(gw_level_ini(np))
   allocate(gw_ground_level(np))
   allocate(gw_infiltration_cap(np))
   allocate(gw_specific_yield(np))
   allocate(gw_space(np))
   allocate(gw_recharge(np))
   allocate(gw_cumseep(np))
   allocate(gw_initial_depth(np))
   allocate(gw_conductivity(np))
   allocate(gw_aquifer_thickness(np))
   allocate(gw_cell_type(np))
   allocate(gw_surface_exchange(np))
   allocate(gw_lateral_volume(np))
   !
   gw_level = 0.0
   gw_level_ini = 0.0
   gw_ground_level = 0.0
   gw_infiltration_cap = 0.0
   gw_specific_yield = 0.0
   gw_space = 0.0
   gw_recharge = 0.0
   gw_cumseep = 0.0
   gw_initial_depth = 0.0
   gw_conductivity = 0.0
   gw_aquifer_thickness = 0.0
   gw_cell_type = 0
   gw_surface_exchange = 0.0
   gw_lateral_volume = 0.0
   gw_time_acc = 0.0
   gw_dt_stable = 0.0
   !
   ! Ground level: zb on regular grids, lowest subgrid pixel on subgrid grids
   !
   do nm = 1, np
      !
      if (subgrid) then
         gw_ground_level(nm) = subgrid_z_zmin(nm)
      else
         gw_ground_level(nm) = zb(nm)
      endif
      !
   enddo
   !
   ! 1) Optional spatially-varying fields from inffile
   !
   have_depth = .false.
   have_level = .false.
   have_sy = .false.
   have_k = .false.
   have_b = .false.
   !
   if (netcdf_infiltration) then
      !
      varname = 'gw_initial_depth'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_initial_depth)
         have_depth = .true.
      endif
      !
      varname = 'gw_initial_level'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         allocate(gw_initial_level(np))
         gw_initial_level = 0.0
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_initial_level)
         have_level = .true.
      endif
      !
      varname = 'gw_specific_yield'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_specific_yield)
         have_sy = .true.
      endif
      !
      varname = 'gw_conductivity'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_conductivity)
         gw_conductivity = gw_conductivity / 86400.0   ! m/day to m/s
         have_k = .true.
      endif
      !
      varname = 'gw_aquifer_thickness'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_aquifer_thickness)
         have_b = .true.
      endif
      !
   endif
   !
   ! 2) Uniform overrides from sfincs.inp (already converted to SI in sfincs_input)
   !
   if (gw_initial_depth_uniform >= 0.0) then
      gw_initial_depth = gw_initial_depth_uniform
      have_depth = .true.
   endif
   !
   if (gw_initial_level_uniform > -998.0) then
      if (.not. allocated(gw_initial_level)) allocate(gw_initial_level(np))
      gw_initial_level = gw_initial_level_uniform
      have_level = .true.
   endif
   !
   if (gw_specific_yield_uniform >= 0.0) then
      gw_specific_yield = gw_specific_yield_uniform
      have_sy = .true.
   endif
   !
   if (gw_conductivity_uniform >= 0.0) then
      gw_conductivity = gw_conductivity_uniform
      have_k = .true.
   endif
   !
   if (gw_aquifer_thickness_uniform >= 0.0) then
      gw_aquifer_thickness = gw_aquifer_thickness_uniform
      have_b = .true.
   endif
   !
   ! 3) Defaults and checks
   !
   if (have_level) then
      !
      ! Absolute water-table elevation given: depth follows from the ground level
      !
      do nm = 1, np
         gw_initial_depth(nm) = gw_ground_level(nm) - gw_initial_level(nm)
      enddo
      !
      if (have_depth) then
         call write_log('Info    : gw_initial_level given, gw_initial_depth is ignored', 0)
      endif
      !
      deallocate(gw_initial_level)
      !
   elseif (.not. have_depth) then
      !
      gw_initial_depth = 1.0
      call write_log('Warning : gw_initial_depth not specified, using initial depth to groundwater of 1.0 m', 0)
      !
   endif
   !
   if (.not. have_sy) then
      gw_specific_yield = 0.3
   endif
   !
   ! Lateral flow is on when both the conductivity and the aquifer thickness are given
   !
   gw_lateral = have_k .and. have_b
   !
   if (gw_lateral) then
      call write_log('Info    : 2D lateral groundwater flow on (gw_conductivity and gw_aquifer_thickness given)', 0)
   elseif (have_k .or. have_b) then
      call stop_sfincs('Error ! lateral groundwater flow needs both gw_conductivity and gw_aquifer_thickness (keyword or inffile field) !', 1)
   else
      call write_log('Info    : no lateral groundwater flow, the aquifer acts as storage only', 0)
   endif
   !
   ! 4) Cell types: 0 inactive, 1 aquifer, 2 open water (ground below qinf_zmin, or
   !    wet at the start). Initially wet cells are open water unless switched off, or
   !    when starting from a restart file, where the wet cells may be a flooded floodplain.
   !
   wet_is_open_water = (gw_initial_wet_open_water == 1)
   !
   if (wet_is_open_water .and. rstfile(1:4) /= 'none') then
      wet_is_open_water = .false.
      call write_log('Info    : restart run, initially wet cells are not treated as open water for groundwater', 0)
   endif
   !
   n_wet_open = 0
   n_floored = 0
   !
   do nm = 1, np
      !
      if (kcs(nm) /= 1) then
         gw_cell_type(nm) = 0
      elseif (gw_ground_level(nm) <= qinf_zmin) then
         gw_cell_type(nm) = 2
      elseif (wet_is_open_water .and. real(zs(nm), 4) > gw_ground_level(nm) + 1.0e-3) then
         gw_cell_type(nm) = 2
         n_wet_open = n_wet_open + 1
      elseif (gw_specific_yield(nm) <= 0.0) then
         gw_cell_type(nm) = 0
      else
         gw_cell_type(nm) = 1
      endif
      !
      gw_specific_yield(nm) = max(gw_specific_yield(nm), 1.0e-3)
      !
   enddo
   !
   ! 5) Seed the water table from gw_initial_depth (a restart replaces the level)
   !
   do nm = 1, np
      !
      if (gw_cell_type(nm) == 2) then
         gw_level_ini(nm) = max(real(zs(nm), 4), gw_ground_level(nm))
      else
         !
         gw_level_ini(nm) = gw_ground_level(nm) - max(gw_initial_depth(nm), 0.0)
         !
         ! The table cannot sit below the regional water level zsini (sea, canals),
         ! unless an explicit level was given or zsini represents standing flood water
         ! (gw_initial_wet_open_water = 0). Capped at the ground level.
         !
         if (wet_is_open_water .and. .not. have_level .and. gw_level_ini(nm) < zini) then
            gw_level_ini(nm) = min(zini, gw_ground_level(nm))
            n_floored = n_floored + 1
         endif
         !
      endif
      !
      gw_level(nm) = gw_level_ini(nm)
      !
   enddo
   !
   if (n_floored > 0) then
      write(logstr,'(a,i0,a,f8.3,a)')'Info    : ', n_floored, ' cells had their initial water table raised to zsini = ', zini, ' m'
      call write_log(logstr, 0)
   endif
   !
   if (n_wet_open > 0) then
      write(logstr,'(a,i0,a)')'Info    : ', n_wet_open, ' cells wet at the start are treated as open water for groundwater'
      call write_log(logstr, 0)
   endif
   !
   if (allocated(gw_level_rst)) then
      !
      gw_level = gw_level_rst
      deallocate(gw_level_rst)
      !
      if (allocated(gw_recharge_rst)) then
         gw_recharge = gw_recharge_rst
         deallocate(gw_recharge_rst)
      endif
      !
      if (gw_time_acc_rst >= 0.0) then
         gw_time_acc = gw_time_acc_rst    ! keep the phase of the groundwater clock
      endif
      !
   endif
   !
   ! 6) Lateral flow: face conductance T * width / distance per uv point, and the
   !    stable explicit substep 0.5 * Sy * A / sum(conductances) over aquifer cells.
   !
   if (gw_lateral) then
      !
      allocate(gw_face_conductance(npuv))
      allocate(csum(np))
      gw_face_conductance = 0.0
      csum = 0.0
      !
      do ip = 1, npuv
         !
         if (kcuv(ip) /= 1) cycle
         !
         nm = uv_index_z_nm(ip)
         nmu = uv_index_z_nmu(ip)
         !
         if (nm < 1 .or. nmu < 1) cycle
         !
         tr1 = 0.0
         tr2 = 0.0
         if (gw_cell_type(nm) == 1) tr1 = gw_conductivity(nm) * gw_aquifer_thickness(nm)
         if (gw_cell_type(nmu) == 1) tr2 = gw_conductivity(nmu) * gw_aquifer_thickness(nmu)
         !
         if (gw_cell_type(nm) == 1 .and. gw_cell_type(nmu) == 1) then
            if (tr1 + tr2 > 0.0) then
               trf = 2.0 * tr1 * tr2 / (tr1 + tr2)
            else
               trf = 0.0
            endif
         elseif (gw_cell_type(nm) == 1 .and. gw_cell_type(nmu) == 2) then
            trf = tr1
         elseif (gw_cell_type(nm) == 2 .and. gw_cell_type(nmu) == 1) then
            trf = tr2
         else
            trf = 0.0
         endif
         !
         if (trf <= 0.0) cycle
         !
         if (use_quadtree) then
            iref = uv_flags_iref(ip)
            itype = uv_flags_type(ip)
         else
            iref = 1
            itype = 0
         endif
         !
         idir = uv_flags_dir(ip)
         !
         if (crsgeo) then
            if (idir == 0) then
               dist = 1.0 / dxminv(ip)
               width = dyrm(iref)
            else
               dist = 1.0 / dyrinv(iref)
               width = 1.0 / dxminv(ip)
            endif
            if (itype /= 0) dist = 1.5 * dist
         else
            if (idir == 0) then
               width = dyrm(iref)
               if (itype == 0) then
                  dist = 1.0 / dxrinv(iref)
               else
                  dist = 1.0 / dxrinvc(iref)
               endif
            else
               width = dxrm(iref)
               if (itype == 0) then
                  dist = 1.0 / dyrinv(iref)
               else
                  dist = 1.0 / dyrinvc(iref)
               endif
            endif
         endif
         !
         gw_face_conductance(ip) = trf * width / dist
         !
         csum(nm) = csum(nm) + gw_face_conductance(ip)
         csum(nmu) = csum(nmu) + gw_face_conductance(ip)
         !
      enddo
      !
      gw_dt_stable = 1.0e9
      !
      do nm = 1, np
         !
         if (gw_cell_type(nm) /= 1 .or. csum(nm) <= 0.0) cycle
         !
         if (crsgeo) then
            gw_dt_stable = min(gw_dt_stable, 0.5 * gw_specific_yield(nm) * cell_area_m2(nm) / csum(nm))
         else
            gw_dt_stable = min(gw_dt_stable, 0.5 * gw_specific_yield(nm) * cell_area(z_flags_iref(nm)) / csum(nm))
         endif
         !
      enddo
      !
      deallocate(csum)
      !
      write(logstr,'(a,e10.3,a,e10.3,a)')'Info    : gw transmissivity   = ', minval(gw_conductivity * gw_aquifer_thickness, mask = gw_cell_type == 1), ' - ', maxval(gw_conductivity * gw_aquifer_thickness, mask = gw_cell_type == 1), ' m2/s'
      call write_log(logstr, 0)
      !
   else
      !
      allocate(gw_face_conductance(1))
      gw_face_conductance = 0.0
      gw_dt_stable = 1.0e9
      !
   endif
   !
   ! 7) Groundwater time step: gw_dt is an upper bound (default 60 s) so that the
   !    exchange with the surface stays responsive; the stable explicit step of the
   !    lateral flow is used when that is smaller.
   !
   if (gw_dt <= 0.0) then
      gw_dt = 60.0
   endif
   !
   gw_dt = min(gw_dt, gw_dt_stable)
   !
   write(logstr,'(a,f10.3,a)')'Info    : gw time step        = ', gw_dt, ' s'
   call write_log(logstr, 0)
   !
   if (gw_lateral) then
      write(logstr,'(a,f10.3,a,i0,a)')'Info    : gw stable substep   = ', gw_dt_stable, ' s (', max(1, ceiling(gw_dt / max(gw_dt_stable, 1.0e-6))), ' substep(s) per groundwater step)'
      call write_log(logstr, 0)
   endif
   !
   ! 8) Space above the table and the first infiltration cap (conservative: dtmax)
   !
   do nm = 1, np
      !
      if (gw_cell_type(nm) == 1) then
         gw_space(nm) = real(max(real(gw_ground_level(nm), 8) - gw_level(nm), 0.0d0), 4) * gw_specific_yield(nm)
         gw_infiltration_cap(nm) = max(gw_space(nm) - gw_recharge(nm), 0.0) / max(dtmax, 1.0e-3)
      elseif (gw_cell_type(nm) == 2) then
         gw_infiltration_cap(nm) = 0.0
      else
         gw_infiltration_cap(nm) = 1.0e30
      endif
      !
   enddo
   !
   write(logstr,'(a,f8.3,a,f8.3,a)')'Info    : gw initial depth    = ', minval(gw_initial_depth), ' - ', maxval(gw_initial_depth), ' m'
   call write_log(logstr, 0)
   write(logstr,'(a,f8.3,a,f8.3)')'Info    : gw specific yield   = ', minval(gw_specific_yield), ' - ', maxval(gw_specific_yield)
   call write_log(logstr, 0)
   !
   deallocate(gw_initial_depth)
   deallocate(gw_level_ini)
   deallocate(gw_conductivity)
   deallocate(gw_aquifer_thickness)
   !
   end subroutine initialize_groundwater


   subroutine update_groundwater(dt)
   !
   ! Bookkeeping on the surface side, every hydrodynamic step, after the infiltration
   ! update (qinfmap already limited by gw_infiltration_cap and applied to qsrc):
   !
   !    aquifer cells   : accumulate the infiltrated depth, add the exchange set at the
   !                      last groundwater step to qsrc, update the cap for the next step
   !    open-water cells: add the exchange to qsrc (seepage in, or recharge out)
   !    inactive cells  : nothing
!
   ! Then advance the aquifer when the accumulated time reaches gw_dt.
   !
   use sfincs_data
   !
   implicit none
   !
   real*4  :: dt
   !
   integer :: nm
   real*4  :: qin
   real*4  :: area_loc
   !
   !$omp parallel do private(nm, qin, area_loc) schedule(static)
   !$acc parallel present( gw_cell_type, gw_recharge, gw_space, gw_infiltration_cap, gw_surface_exchange, gw_cumseep, qinfmap, qsrc, &
   !$acc                   cell_area, cell_area_m2, z_flags_iref )
   !$acc loop independent gang vector
   do nm = 1, np
      !
      if (gw_cell_type(nm) == 0) cycle
      !
      if (gw_cell_type(nm) == 1 .and. infiltration) then
         !
         qin = max(qinfmap(nm), 0.0)
         gw_recharge(nm) = gw_recharge(nm) + qin * dt
         gw_infiltration_cap(nm) = max(gw_space(nm) - gw_recharge(nm), 0.0) / dt
         !
      endif
      !
      ! Exchange with the surface (positive = water delivered to the surface)
      !
      if (gw_surface_exchange(nm) /= 0.0) then
         !
         if (crsgeo) then
            area_loc = cell_area_m2(nm)
         else
            area_loc = cell_area(z_flags_iref(nm))
         endif
         !
         qsrc(nm) = qsrc(nm) + gw_surface_exchange(nm) * area_loc
         !
         if (store_cumulative_precipitation) then
            gw_cumseep(nm) = gw_cumseep(nm) + gw_surface_exchange(nm) * dt
         endif
         !
      endif
      !
   enddo
!$acc end parallel
   !$omp end parallel do
   !
   gw_time_acc = gw_time_acc + dt
   !
   if (gw_time_acc >= gw_dt - 1.0e-3) then
      !
      call groundwater_step(gw_time_acc, dt)
      gw_time_acc = 0.0
      !
   endif
   !
   end subroutine update_groundwater


   subroutine groundwater_step(dtgw, dt)
   !
   ! Advance the aquifer over dtgw (the time since the last groundwater step):
   !
   !    1) recharge: H += gw_recharge / Sy
   !    2) lateral flow (when on): explicit substeps within the stability limit,
   !       Sy A dH = sum over faces of C_f (H_nb - H) dts; open-water cells hold H = zs
   !    3) exfiltration: H above the ground returns the excess to the surface
   !    4) the delivered depth becomes the exchange rate gw_surface_exchange for the
   !       coming interval; the space above the table and the cap for the next
   !       hydrodynamic step (length dt) are renewed; gw_recharge is reset
   !
   use sfincs_data
   !
   implicit none
   !
   real*4  :: dtgw
   real*4  :: dt
   !
   integer :: nsub
   integer :: isub
   integer :: ip
   integer :: nm
   integer :: nmu
   real*4  :: dts
   real*4  :: hnm
   real*4  :: hnmu
   real*4  :: qf
   real*4  :: area_loc
   real*4  :: excess
   !
   ! 1) Recharge, and reset the delivered depth (kept in gw_surface_exchange as a depth
   !    until step 4 turns it into a rate)
   !
   !$omp parallel do private(nm) schedule(static)
   !$acc parallel present( gw_cell_type, gw_level, gw_recharge, gw_specific_yield, gw_surface_exchange, zs )
   !$acc loop independent gang vector
   do nm = 1, np
      !
      gw_surface_exchange(nm) = 0.0
      !
      if (gw_cell_type(nm) == 1) then
         gw_level(nm) = gw_level(nm) + gw_recharge(nm) / gw_specific_yield(nm)
         gw_recharge(nm) = 0.0
      elseif (gw_cell_type(nm) == 2) then
         gw_level(nm) = zs(nm)
      endif
      !
   enddo
   !$acc end parallel
   !$omp end parallel do
   !
   ! 2) Lateral flow
   !
   if (gw_lateral) then
      !
      nsub = max(1, ceiling(dtgw / max(gw_dt_stable, 1.0e-6)))
      dts = dtgw / nsub
      !
      do isub = 1, nsub
         !
         !$acc parallel loop present( gw_lateral_volume )
         !$omp parallel do private(nm)
         do nm = 1, np
            gw_lateral_volume(nm) = 0.0
         enddo
         !$omp end parallel do
         !
         ! Face fluxes (m3 over the substep, positive from nm to nmu)
         !
         !$omp parallel do private(ip, nm, nmu, hnm, hnmu, qf)
         !$acc parallel present( gw_face_conductance, gw_cell_type, gw_level, zs, gw_lateral_volume, uv_index_z_nm, uv_index_z_nmu )
         !$acc loop independent gang vector
         do ip = 1, npuv
            !
            if (gw_face_conductance(ip) <= 0.0) cycle
            !
            nm = uv_index_z_nm(ip)
            nmu = uv_index_z_nmu(ip)
            !
            hnm = real(gw_level(nm), 4)
            hnmu = real(gw_level(nmu), 4)
            !
            qf = gw_face_conductance(ip) * (hnm - hnmu) * dts
            !
            !$acc atomic update
            !$omp atomic
            gw_lateral_volume(nm) = gw_lateral_volume(nm) - qf
            !$acc atomic update
            !$omp atomic
            gw_lateral_volume(nmu) = gw_lateral_volume(nmu) + qf
            !
         enddo
         !$acc end parallel
         !$omp end parallel do
         !
         ! Apply to the water table; open-water cells collect what they receive
         !
         !$omp parallel do private(nm, area_loc)
         !$acc parallel present( gw_cell_type, gw_lateral_volume, gw_level, gw_specific_yield, gw_surface_exchange, cell_area, cell_area_m2, z_flags_iref )
         !$acc loop independent gang vector
         do nm = 1, np
            !
            if (gw_cell_type(nm) == 0) cycle
            !
            if (crsgeo) then
               area_loc = cell_area_m2(nm)
            else
               area_loc = cell_area(z_flags_iref(nm))
            endif
            !
            if (gw_cell_type(nm) == 1) then
               gw_level(nm) = gw_level(nm) + gw_lateral_volume(nm) / (gw_specific_yield(nm) * area_loc)
            else
               gw_surface_exchange(nm) = gw_surface_exchange(nm) + gw_lateral_volume(nm) / area_loc
            endif
            !
         enddo
         !$acc end parallel
         !$omp end parallel do
         !
      enddo
      !
   endif
   !
   ! 3) Exfiltration where the table is above the ground, 4) exchange rate, space and cap
   !
   !$omp parallel do private(nm, excess) schedule(static)
   !$acc parallel present( gw_cell_type, gw_level, gw_ground_level, gw_specific_yield, gw_surface_exchange, gw_space, gw_infiltration_cap )
   !$acc loop independent gang vector
   do nm = 1, np
      !
      if (gw_cell_type(nm) == 1) then
         !
         if (gw_level(nm) > gw_ground_level(nm)) then
            excess = real(gw_level(nm) - gw_ground_level(nm), 4) * gw_specific_yield(nm)
            gw_surface_exchange(nm) = gw_surface_exchange(nm) + excess
            gw_level(nm) = gw_ground_level(nm)
         endif
         !
         gw_space(nm) = real(max(real(gw_ground_level(nm), 8) - gw_level(nm), 0.0d0), 4) * gw_specific_yield(nm)
         gw_infiltration_cap(nm) = gw_space(nm) / dt
         !
      endif
      !
      if (gw_cell_type(nm) /= 0) then
         gw_surface_exchange(nm) = gw_surface_exchange(nm) / dtgw
      endif
      !
   enddo
   !$acc end parallel
   !$omp end parallel do
   !
   end subroutine groundwater_step

end module sfincs_groundwater
