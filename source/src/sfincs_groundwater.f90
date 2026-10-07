module sfincs_groundwater
   !
   ! Groundwater table model (Sanders et al. 2025, PRIMo). Each active cell carries a
   ! water table that rises by infiltration and drains by seepage. Selected with
   ! inftype = gwt; the infiltration module dispatches here and applies qinfmap.
   !
   ! State and parameters live in sfincs_data (gw_* arrays and keywords).
   !
   !   initialize_groundwater()
   !     Allocates the gw_* arrays, reads optional fields from inffile, applies
   !     the uniform gw_* overrides, maps receiver cells, seeds the water table
   !     (also from a restart). Called from initialize_infiltration.
   !
   !   update_groundwater(dt)
   !     One step: vertical exchange with the surface (sets qinfmap), seepage
   !     sink to the receiver cell / lost / local. Called from
   !     update_infiltration_map once per time step.
   !
   use sfincs_log
   use sfincs_error
   !
contains

   subroutine initialize_groundwater()
   !
   ! 0D groundwater table model (Sanders et al. 2025, PRIMo, eqs. 5-8).
   !
   ! Per cell:  Sy dH/dt = phi * f - w,   w = kappa * max(H - H0, 0),   kappa = Keff / l0
   !
   ! Allocates the state and parameter arrays, reads optional fields from
   ! inffile, applies the uniform overrides from sfincs.inp, disables the
   ! aquifer in open-water and inactive cells, and seeds the water table
   ! (from gw_depth0, or from the restart file when present).
   !
   use sfincs_data
   use sfincs_ncinput
   use quadtree
   !
   implicit none
   !
   integer       :: nm
   integer       :: ipq
   integer       :: n_no_receiver
   character*256 :: varname
   logical       :: have_depth, have_fmax, have_phi, have_sy, have_keff, have_l0
   !
   real*4, dimension(:), allocatable :: gw_depth0
   real*4, dimension(:), allocatable :: gw_keff
   real*4, dimension(:), allocatable :: gw_l0
   integer*4, dimension(:), allocatable :: gw_receiver_qt
   !
   call write_log('Info    : turning on groundwater table model (0D, PRIMo-style)', 0)
   !
   allocate(gw_rise(np))
   allocate(gw_level(np))
   allocate(gw_level0(np))
   allocate(gw_zground(np))
   allocate(gw_fmax(np))
   allocate(gw_phi(np))
   allocate(gw_sy(np))
   allocate(gw_kappa(np))
   allocate(gw_seepage(np))
   allocate(gw_cumseep(np))
   allocate(gw_receiver(np))
   allocate(gw_depth0(np))
   allocate(gw_keff(np))
   allocate(gw_l0(np))
   !
   gw_rise = 0.0
   gw_level = 0.0
   gw_level0 = 0.0
   gw_zground = 0.0
   gw_fmax = 0.0
   gw_phi = 0.0
   gw_sy = 0.0
   gw_kappa = 0.0
   gw_seepage = 0.0
   gw_cumseep = 0.0
   gw_receiver = 0
   gw_depth0 = 0.0
   gw_keff = 0.0
   gw_l0 = 0.0
   !
   ! Ground level: zb on regular grids, lowest subgrid pixel on subgrid grids
   !
   do nm = 1, np
      !
      if (subgrid) then
         gw_zground(nm) = subgrid_z_zmin(nm)
      else
         gw_zground(nm) = zb(nm)
      endif
      !
   enddo
   !
   ! 1) Optional spatially-varying fields from inffile
   !
   have_depth = .false.
   have_fmax = .false.
   have_phi = .false.
   have_sy = .false.
   have_keff = .false.
   have_l0 = .false.
   !
   if (netcdf_infiltration) then
      !
      varname = 'gw_depth0'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_depth0)
         have_depth = .true.
      endif
      !
      varname = 'gw_fmax'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_fmax)
         gw_fmax = gw_fmax / 3.6e6   ! mm/hr to m/s
         have_fmax = .true.
      endif
      !
      varname = 'gw_phi'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_phi)
         have_phi = .true.
      endif
      !
      varname = 'gw_sy'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_sy)
         have_sy = .true.
      endif
      !
      varname = 'gw_keff'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_keff)
         gw_keff = gw_keff / 3.6e6   ! mm/hr to m/s
         have_keff = .true.
      endif
      !
      varname = 'gw_l0'
      if (netcdf_quadtree_variable_exists(inffile, varname)) then
         call read_netcdf_quadtree_to_sfincs(inffile, varname, gw_l0)
         have_l0 = .true.
      endif
      !
   endif
   !
   ! 2) Uniform overrides from sfincs.inp (already converted to SI in sfincs_input)
   !
   if (gw_depth_ini >= 0.0) then
      gw_depth0 = gw_depth_ini
      have_depth = .true.
   endif
   !
   if (gw_fmax_uniform >= 0.0) then
      gw_fmax = gw_fmax_uniform
      have_fmax = .true.
   endif
   !
   if (gw_phi_uniform >= 0.0) then
      gw_phi = gw_phi_uniform
      have_phi = .true.
   endif
   !
   if (gw_sy_uniform >= 0.0) then
      gw_sy = gw_sy_uniform
      have_sy = .true.
   endif
   !
   if (gw_keff_uniform >= 0.0) then
      gw_keff = gw_keff_uniform
      have_keff = .true.
   endif
   !
   if (gw_l0_uniform >= 0.0) then
      gw_l0 = gw_l0_uniform
      have_l0 = .true.
   endif
   !
   ! 3) Defaults and checks
   !
   if (.not. have_depth) then
      gw_depth0 = 1.0
      call write_log('Warning : gw_depth0 not specified, using initial depth to groundwater of 1.0 m', 0)
   endif
   !
   if (.not. have_phi) then
      gw_phi = 1.0
   endif
   !
   if (.not. have_sy) then
      gw_sy = 0.3
   endif
   !
   if (.not. have_fmax) then
      call stop_sfincs('Error ! inftype = gwt requires gw_fmax in sfincs.inp or a gw_fmax field in inffile !', 1)
   endif
   !
   if (have_keff .and. have_l0) then
      !
      do nm = 1, np
         if (gw_l0(nm) > 0.0) then
            gw_kappa(nm) = gw_keff(nm) / gw_l0(nm)
         else
            gw_kappa(nm) = 0.0
         endif
      enddo
      !
   else
      !
      gw_kappa = 0.0
      call write_log('Warning : gw_keff and/or gw_l0 not specified, groundwater seepage is switched off', 0)
      !
   endif
   !
   ! 4) Seepage destination
   !
   select case (trim(gw_seepage_mode))
   case ('loss')
      gw_seepage_imode = 0
   case ('local')
      gw_seepage_imode = 1
   case ('receiver')
      gw_seepage_imode = 2
   case default
      call stop_sfincs('Error ! gw_seepage_mode should be receiver, loss or local !', 1)
   end select
   !
   if (gw_seepage_imode == 2) then
      !
      ! Receiver cells: gw_receiver in inffile holds the 1-based quadtree index of the
      ! surface-water cell that receives the seepage of each cell. Map to sfincs indices;
      ! cells without a valid receiver lose their seepage.
      !
      varname = 'gw_receiver'
      !
      if (.not. netcdf_infiltration) then
         call stop_sfincs('Error ! gw_seepage_mode = receiver requires inffile with a gw_receiver field !', 1)
      endif
      !
      if (.not. netcdf_quadtree_variable_exists(inffile, varname)) then
         call stop_sfincs('Error ! gw_seepage_mode = receiver requires a gw_receiver field in inffile !', 1)
      endif
      !
      allocate(gw_receiver_qt(np))
      gw_receiver_qt = 0
      !
      call read_netcdf_quadtree_to_sfincs_int(inffile, varname, gw_receiver_qt)
      !
      n_no_receiver = 0
      !
      do nm = 1, np
         !
         ipq = gw_receiver_qt(nm)
         !
         if (ipq >= 1 .and. ipq <= quadtree_nr_points) then
            gw_receiver(nm) = index_sfincs_in_quadtree(ipq)
         else
            gw_receiver(nm) = 0
         endif
         !
         if (gw_receiver(nm) <= 0 .and. gw_zground(nm) > qinf_zmin .and. kcs(nm) == 1) then
            n_no_receiver = n_no_receiver + 1
         endif
         !
      enddo
      !
      deallocate(gw_receiver_qt)
      !
      if (n_no_receiver > 0) then
         write(logstr,'(a,i0,a)')'Warning : ', n_no_receiver, ' groundwater cells have no valid receiver cell, their seepage is lost'
         call write_log(logstr, 0)
      endif
      !
   endif
   !
   ! 5) No aquifer in open-water (ground below qinf_zmin) or inactive cells.
   !    phi = 0 switches the cell off in update_groundwater.
   !
   do nm = 1, np
      !
      if (gw_zground(nm) <= qinf_zmin .or. kcs(nm) /= 1 .or. gw_sy(nm) <= 0.0) then
         gw_phi(nm) = 0.0
         gw_kappa(nm) = 0.0
      endif
      !
      gw_sy(nm) = max(gw_sy(nm), 1.0e-3)
      !
   enddo
   !
   ! 6) Seed the water table. The seepage baseline always follows gw_depth0;
   !    a restart only replaces the current level.
   !
   do nm = 1, np
      gw_level0(nm) = gw_zground(nm) - max(gw_depth0(nm), 0.0)
      gw_level(nm) = gw_level0(nm)
   enddo
   !
   if (allocated(gw_rise_rst)) then
      !
      gw_rise = gw_rise_rst
      gw_level = gw_level0 + gw_rise
      deallocate(gw_rise_rst)
      !
   endif
   !
   write(logstr,'(a,f8.3,a,f8.3,a)')'Info    : gw initial depth    = ', minval(gw_depth0), ' - ', maxval(gw_depth0), ' m'
   call write_log(logstr, 0)
   write(logstr,'(a,f8.3,a,f8.3,a)')'Info    : gw fmax             = ', minval(gw_fmax) * 3.6e6, ' - ', maxval(gw_fmax) * 3.6e6, ' mm/hr'
   call write_log(logstr, 0)
   write(logstr,'(a,f8.3,a,f8.3)')'Info    : gw pervious fraction= ', minval(gw_phi), ' - ', maxval(gw_phi)
   call write_log(logstr, 0)
   write(logstr,'(a,f8.3,a,f8.3)')'Info    : gw specific yield   = ', minval(gw_sy), ' - ', maxval(gw_sy)
   call write_log(logstr, 0)
   write(logstr,'(a,e10.3,a,e10.3,a)')'Info    : gw kappa            = ', minval(gw_kappa), ' - ', maxval(gw_kappa), ' 1/s'
   call write_log(logstr, 0)
   write(logstr,'(a,a)')'Info    : gw seepage mode     = ', trim(gw_seepage_mode)
   call write_log(logstr, 0)
   !
   deallocate(gw_depth0)
   deallocate(gw_keff)
   deallocate(gw_l0)
   !
   end subroutine initialize_groundwater


   subroutine update_groundwater(dt)
   !
   ! Groundwater table model, one time step (0D vertical exchange + seepage).
   !
   !   depth to groundwater d = z_ground - H
   !   f = fmax             if d > 0 and water on the surface      (inundation case)
   !   f = min(prcp, fmax)  if d > 0 and no water on the surface   (rainfall case)
   !   f = 0                if d <= 0                              (saturation case)
   !   f is capped by the water available this step and by the aquifer space.
   !
   !   Surface sink  qinf = phi * f   (mass-consistent with the aquifer gain)
   !   Aquifer       Sy dH/dt = qinf - kappa * (H - H0), solved exactly over dt with
   !                 qinf held constant; the seepage w follows from the balance.
   !                 The state is the rise H - H0 (gw_rise) so that micrometre
   !                 changes are resolved in single precision; gw_level is derived.
   !
   !   Seepage goes to: nowhere (loss), the same cell (local), or a receiver
   !   cell via qsrc (receiver). qinfmap is picked up by update_infiltration_map.
   !
   use sfincs_data
   !
   implicit none
   !
   real*4  :: dt
   !
   integer :: nm
   integer :: nmr
   real*4  :: hh_local
   real*4  :: pr
   real*4  :: depth
   real*4  :: f
   real*4  :: qinf_loc
   real*4  :: rise
   real*4  :: rise_new
   real*4  :: afac
   real*4  :: rfac
   real*4  :: w
   real*4  :: area_loc
   !
   !$omp parallel do private(nm, nmr, hh_local, pr, depth, f, qinf_loc, rise, rise_new, afac, rfac, w, area_loc) schedule(static)
   !$acc parallel present( kcs, prcp, zs, zb, z_volume, cell_area, cell_area_m2, z_flags_iref, subgrid_z_zmin, qinfmap, cuminf, qsrc, &
   !$acc                   gw_rise, gw_level, gw_level0, gw_zground, gw_fmax, gw_phi, gw_sy, gw_kappa, gw_seepage, gw_cumseep, gw_receiver )
   !$acc loop independent gang vector
   do nm = 1, np
      !
      qinf_loc = 0.0
      w = 0.0
      !
      if (gw_phi(nm) > 0.0) then
         !
         ! Water available on the surface (m)
         !
         if (subgrid) then
            !
            if (crsgeo) then
               hh_local = z_volume(nm) / cell_area_m2(nm)
            else
               hh_local = z_volume(nm) / cell_area(z_flags_iref(nm))
            endif
            !
         else
            !
            hh_local = zs(nm) - zb(nm)
            !
         endif
         !
         hh_local = max(hh_local, 0.0)
         if (precip) then
            pr = max(prcp(nm), 0.0)
         else
            pr = 0.0
         endif
         !
         depth = (gw_zground(nm) - gw_level0(nm)) - gw_rise(nm)
         !
         if (depth > 0.0) then
            !
            if (hh_local > 0.0) then
               f = gw_fmax(nm)
            else
               f = min(pr, gw_fmax(nm))
            endif
            !
            ! Never remove more than is available on the surface this step
            !
            f = min(f, (hh_local + pr * dt) / dt)
            !
            ! Never overfill the aquifer in one step
            !
            f = min(f, depth * gw_sy(nm) / (gw_phi(nm) * dt))
            f = max(f, 0.0)
            !
            qinf_loc = gw_phi(nm) * f
            !
         endif
         !
         ! Aquifer update: exact solution of Sy dH/dt = qinf_loc - kappa (H - H0) over dt
         !
         rise = gw_rise(nm)
         !
         if (gw_kappa(nm) > 0.0) then
            !
            ! a = kappa dt / Sy, r = (1 - exp(-a)) / a. The seepage over the step is
            ! w = q (1 - r) + kappa rise r, which is exact for both a -> 0 (w -> kappa rise)
            ! and a -> inf (steady state q / kappa) and avoids the 1 - exp(-a)
            ! cancellation in single precision for very small time steps.
            !
            afac = gw_kappa(nm) * dt / gw_sy(nm)
            !
            if (afac > 1.0e-4) then
               rfac = (1.0 - exp(-afac)) / afac
            else
               rfac = 1.0 - 0.5 * afac + afac * afac / 6.0
            endif
            !
            w = max(qinf_loc * (1.0 - rfac) + gw_kappa(nm) * rise * rfac, 0.0)
            rise_new = rise + (qinf_loc - w) * dt / gw_sy(nm)
            !
         else
            !
            rise_new = rise + qinf_loc * dt / gw_sy(nm)
            !
         endif
         !
         gw_rise(nm) = rise_new
         gw_level(nm) = gw_level0(nm) + rise_new
         !
      endif
      !
      gw_seepage(nm) = w
      gw_cumseep(nm) = gw_cumseep(nm) + w * dt
      !
      if (gw_seepage_imode == 1) then
         qinfmap(nm) = qinf_loc - w
      else
         qinfmap(nm) = qinf_loc
      endif
      !
      if (store_cumulative_precipitation) then
         cuminf(nm) = cuminf(nm) + qinfmap(nm) * dt
      endif
      !
      ! Route seepage to the receiving cell
      !
      if (gw_seepage_imode == 2 .and. w > 0.0) then
         !
         nmr = gw_receiver(nm)
         !
         if (nmr > 0) then
            !
            if (crsgeo) then
               area_loc = cell_area_m2(nm)
            else
               area_loc = cell_area(z_flags_iref(nm))
            endif
            !
            !$acc atomic update
            !$omp atomic
            qsrc(nmr) = qsrc(nmr) + w * area_loc
            !
         endif
         !
      endif
      !
   enddo
   !$acc end parallel
   !$omp end parallel do
   !
   end subroutine update_groundwater

end module sfincs_groundwater
