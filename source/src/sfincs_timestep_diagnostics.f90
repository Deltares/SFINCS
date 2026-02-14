module sfincs_timestep_diagnostics
   !
   use sfincs_data
   use sfincs_log
   !
   implicit none
   !
   integer, parameter :: UNIT_TIMESTEP_DIAGNOSTICS = 971
   !
   logical, save :: timestep_diag_initialized = .false.
   real*8,  save :: timestep_diag_next_write_time = 0.0d0
   !
   integer*4, save :: timestep_diag_total_steps = 0
   integer*4, save :: timestep_diag_dtmax_steps = 0
   !
   integer*4, dimension(:), allocatable, save :: timestep_diag_times_limiting_uv
   integer*4, dimension(:), allocatable, save :: timestep_diag_times_limiting_wave
   integer*4, dimension(:), allocatable, save :: timestep_diag_times_limiting_vel
   real*4,    dimension(:), allocatable, save :: timestep_diag_min_limited_dt_uv
   !
contains
   !
   subroutine timestep_diagnostics_initialize()
   !
   implicit none
   !
   if (.not. timestep_diagnostics) then
      return
   endif
   !
   if (timestep_diag_initialized) then
      return
   endif
   !
   timestep_diag_total_steps = 0
   timestep_diag_dtmax_steps = 0
   !
   allocate(timestep_diag_times_limiting_uv(npuv))
   allocate(timestep_diag_times_limiting_wave(npuv))
   allocate(timestep_diag_times_limiting_vel(npuv))
   allocate(timestep_diag_min_limited_dt_uv(npuv))
   !
   timestep_diag_times_limiting_uv   = 0
   timestep_diag_times_limiting_wave = 0
   timestep_diag_times_limiting_vel  = 0
   timestep_diag_min_limited_dt_uv   = 1.0e30
   !
   timestep_diag_next_write_time = real(t0, kind(timestep_diag_next_write_time))
   !
   open(unit = UNIT_TIMESTEP_DIAGNOSTICS, status = 'replace', file = 'timestep_diagnostics.csv')
   write(UNIT_TIMESTEP_DIAGNOSTICS,'(a)') &
      'nt,t,dt_used,dt_next,min_dt,dt_uv_min,reason,ip,nm,nmu,x,y,hu,abs_uv,c_wave,dx,iref,itype,idir'
   !
   timestep_diag_initialized = .true.
   !
   call write_log('Info    : timestep diagnostics enabled (timestep_diagnostics.csv, timestep_diagnostics_domain.csv)', 0)
   !
   end subroutine
   !
   subroutine timestep_diagnostics_update(nt, t, dt_used)
   !
   integer, intent(in) :: nt
   real*8, intent(in) :: t
   real*4, intent(in) :: dt_used
   !
   integer*4 :: ip_lim, nm_lim, nmu_lim
   integer*4 :: iref_lim, itype_lim, idir_lim
   integer*4 :: reason
   real*4    :: dt_uv_min
   real*4    :: hu_lim, abs_uv_lim, c_wave_lim, dx_lim
   real*4    :: x_lim, y_lim
   real*4    :: dt_next
   !
   if (.not. timestep_diagnostics .and. .not. timestep_analysis) then
      return
   endif
   !
   if (timestep_diagnostics) then
      if (.not. timestep_diag_initialized) then
         call timestep_diagnostics_initialize()
      endif
   endif
   !
   !$acc update host(zs, uv)
   !
   call find_timestep_limiter_uv(dt_uv_min, reason, ip_lim, nm_lim, nmu_lim, x_lim, y_lim, hu_lim, abs_uv_lim, &
                                 c_wave_lim, dx_lim, iref_lim, itype_lim, idir_lim, timestep_analysis, timestep_diagnostics)
   !
   if (timestep_diagnostics) then
      timestep_diag_total_steps = timestep_diag_total_steps + 1
      !
      if (reason == 2) then
         timestep_diag_dtmax_steps = timestep_diag_dtmax_steps + 1
      else
         timestep_diag_times_limiting_uv(ip_lim) = timestep_diag_times_limiting_uv(ip_lim) + 1
         if (reason == 0) then
            timestep_diag_times_limiting_wave(ip_lim) = timestep_diag_times_limiting_wave(ip_lim) + 1
         elseif (reason == 1) then
            timestep_diag_times_limiting_vel(ip_lim) = timestep_diag_times_limiting_vel(ip_lim) + 1
         endif
         timestep_diag_min_limited_dt_uv(ip_lim) = min(timestep_diag_min_limited_dt_uv(ip_lim), min_dt)
      endif
   endif
   !
   dt_next = alfa * min_dt
   !
   if (timestep_diagnostics) then
      if (dt_timestep_diagnostics <= 0.0) then
         call write_timestep_diagnostics_line(nt, t, dt_used, dt_next, dt_uv_min, reason, ip_lim, nm_lim, nmu_lim, &
                                              x_lim, y_lim, hu_lim, abs_uv_lim, c_wave_lim, dx_lim, &
                                              iref_lim, itype_lim, idir_lim)
      else
         if (t + 1.0d-9 >= timestep_diag_next_write_time) then
            call write_timestep_diagnostics_line(nt, t, dt_used, dt_next, dt_uv_min, reason, ip_lim, nm_lim, nmu_lim, &
                                                 x_lim, y_lim, hu_lim, abs_uv_lim, c_wave_lim, dx_lim, &
                                                 iref_lim, itype_lim, idir_lim)
            do while (t + 1.0d-9 >= timestep_diag_next_write_time)
               timestep_diag_next_write_time = timestep_diag_next_write_time + &
                                               real(dt_timestep_diagnostics, kind(timestep_diag_next_write_time))
            end do
         endif
      endif
   endif
   !
   if (timestep_analysis) then
      call update_timestep_analysis(dt_next)
   endif
   !
   end subroutine
   !
   subroutine timestep_diagnostics_finalize()
   !
   implicit none
   !
   integer*4 :: ip, nm, nmu
   real*4    :: xuv, yuv
   real*4    :: dt_min_alpha
   !
   if (.not. timestep_diagnostics) then
      return
   endif
   !
   if (timestep_diag_initialized) then
      close(UNIT_TIMESTEP_DIAGNOSTICS)
   endif
   !
   open(unit = UNIT_TIMESTEP_DIAGNOSTICS + 1, status = 'replace', file = 'timestep_diagnostics_domain.csv')
   write(UNIT_TIMESTEP_DIAGNOSTICS + 1,'(a)') &
      'ip,nm,nmu,x,y,times_limiting,times_limiting_wave,times_limiting_vel,dt_min,min_dt_min'
   !
   if (timestep_diag_dtmax_steps > 0) then
      dt_min_alpha = alfa * dtmax
      write(UNIT_TIMESTEP_DIAGNOSTICS + 1,'(i0,a,i0,a,i0,a,es16.8,a,es16.8,a,i0,a,i0,a,i0,a,es16.8,a,es16.8)') &
         0, ',', 0, ',', 0, ',', 0.0, ',', 0.0, ',', timestep_diag_dtmax_steps, ',', 0, ',', 0, ',', dt_min_alpha, ',', dtmax
   endif
   !
   do ip = 1, npuv
      if (timestep_diag_times_limiting_uv(ip) <= 0) then
         cycle
      endif
      nm  = uv_index_z_nm(ip)
      nmu = uv_index_z_nmu(ip)
      xuv = 0.5 * (z_xz(nm) + z_xz(nmu))
      yuv = 0.5 * (z_yz(nm) + z_yz(nmu))
      dt_min_alpha = alfa * timestep_diag_min_limited_dt_uv(ip)
      write(UNIT_TIMESTEP_DIAGNOSTICS + 1,'(i0,a,i0,a,i0,a,es16.8,a,es16.8,a,i0,a,i0,a,i0,a,es16.8,a,es16.8)') &
         ip, ',', nm, ',', nmu, ',', xuv, ',', yuv, ',', timestep_diag_times_limiting_uv(ip), ',', &
         timestep_diag_times_limiting_wave(ip), ',', timestep_diag_times_limiting_vel(ip), ',', dt_min_alpha, ',', &
         timestep_diag_min_limited_dt_uv(ip)
   enddo
   !
   close(UNIT_TIMESTEP_DIAGNOSTICS + 1)
   !
   if (allocated(timestep_diag_times_limiting_uv)) deallocate(timestep_diag_times_limiting_uv)
   if (allocated(timestep_diag_times_limiting_wave)) deallocate(timestep_diag_times_limiting_wave)
   if (allocated(timestep_diag_times_limiting_vel)) deallocate(timestep_diag_times_limiting_vel)
   if (allocated(timestep_diag_min_limited_dt_uv)) deallocate(timestep_diag_min_limited_dt_uv)
   !
   timestep_diag_initialized = .false.
   !
   end subroutine
   !
   subroutine write_timestep_diagnostics_line(nt, t, dt_used, dt_next, dt_uv_min, reason, ip_lim, nm_lim, nmu_lim, x_lim, y_lim, &
                                              hu_lim, abs_uv_lim, c_wave_lim, dx_lim, iref_lim, itype_lim, idir_lim)
   !
   integer, intent(in) :: nt
   real*8,  intent(in) :: t
   real*4,  intent(in) :: dt_used
   real*4,  intent(in) :: dt_next
   real*4,  intent(in) :: dt_uv_min
   integer*4, intent(in) :: reason, ip_lim, nm_lim, nmu_lim, iref_lim, itype_lim, idir_lim
   real*4,  intent(in) :: x_lim, y_lim, hu_lim, abs_uv_lim, c_wave_lim, dx_lim
   !
   write(UNIT_TIMESTEP_DIAGNOSTICS,'(i0,a,es16.8,a,es16.8,a,es16.8,a,es16.8,a,es16.8,a,i0,a,i0,a,i0,a,i0,a,es16.8,a,&
&es16.8,a,es16.8,a,es16.8,a,es16.8,a,es16.8,a,i0,a,i0,a,i0)') &
      nt, ',', t, ',', dt_used, ',', dt_next, ',', min_dt, ',', dt_uv_min, ',', reason, ',', ip_lim, ',', nm_lim, ',', &
      nmu_lim, ',', x_lim, ',', y_lim, ',', hu_lim, ',', abs_uv_lim, ',', c_wave_lim, ',', dx_lim, ',', iref_lim, ',', &
      itype_lim, ',', idir_lim
   !
   end subroutine
   !
   subroutine find_timestep_limiter_uv(dt_uv_min, reason, ip_lim, nm_lim, nmu_lim, x_lim, y_lim, hu_lim, abs_uv_lim, &
                                       c_wave_lim, dx_lim, iref_lim, itype_lim, idir_lim, compute_cell_mins, &
                                       compute_limiter_details)
   !
   real*4,    intent(out) :: dt_uv_min
   integer*4, intent(out) :: reason, ip_lim, nm_lim, nmu_lim
   real*4,    intent(out) :: x_lim, y_lim, hu_lim, abs_uv_lim, c_wave_lim, dx_lim
   integer*4, intent(out) :: iref_lim, itype_lim, idir_lim
   logical,   intent(in)  :: compute_cell_mins
   logical,   intent(in)  :: compute_limiter_details
   !
   integer*4 :: ip, nm, nmu
   integer*4 :: iref, itype, idir
   integer*4 :: iuv
   real*4    :: zsu
   real*4    :: zmin, zmax
   real*4    :: dzuv, facint
   real*4    :: hu
   real*4    :: dxuvinv
   real*4    :: absu, cwave, speed, dtip
   logical   :: iok
   real*4    :: best_dt_uv
   !
   best_dt_uv = 1.0e30
   ip_lim     = 0
   nm_lim     = 0
   nmu_lim    = 0
   reason     = 2
   dt_uv_min  = -1.0
   x_lim      = 0.0
   y_lim      = 0.0
   hu_lim     = 0.0
   abs_uv_lim = 0.0
   c_wave_lim = 0.0
   dx_lim     = 0.0
   iref_lim   = 0
   itype_lim  = 0
   idir_lim   = 0
   !
   if (.not. compute_cell_mins .and. .not. compute_limiter_details) then
      return
   endif
   !
   if (compute_cell_mins) then
      if (allocated(min_timestep)) then
         min_timestep = 1.0e30
      endif
   endif
   !
   do ip = 1, npuv
      !
      if (.not. (kcuv(ip) == 1 .or. kcuv(ip) == 6)) then
         cycle
      endif
      !
      nm  = uv_index_z_nm(ip)
      nmu = uv_index_z_nmu(ip)
      !
      zsu = max(real(zs(nm), kind(zsu)), real(zs(nmu), kind(zsu)))
      !
      iok = .false.
      !
      if (subgrid) then
         zmin = subgrid_uv_zmin(ip)
         if (zsu > zmin + huthresh) then
            iok = .true.
         endif
      else
         if (zsu > zbuvmx(ip)) then
            iok = .true.
         endif
      endif
      !
      if (.not. iok) then
         cycle
      endif
      !
      if (use_quadtree) then
         iref  = uv_flags_iref(ip)
         itype = uv_flags_type(ip)
      else
         iref  = 1
         itype = 0
      endif
      !
      idir = uv_flags_dir(ip)
      !
      if (crsgeo) then
         !
         if (itype == 0) then
            if (idir == 0) then
               dxuvinv = dxminv(ip)
            else
               dxuvinv = dyrinv(iref)
            endif
         else
            if (idir == 0) then
               dxuvinv = 1.0 / (3*(1.0/dxminv(ip))/2)
            else
               dxuvinv = 1.0 / (3*(1.0/dyrinv(iref))/2)
            endif
         endif
         !
      else
         !
         if (itype == 0) then
            if (idir == 0) then
               dxuvinv = dxrinv(iref)
            else
               dxuvinv = dyrinv(iref)
            endif
         else
            if (idir == 0) then
               dxuvinv = dxrinvc(iref)
            else
               dxuvinv = dyrinvc(iref)
            endif
         endif
         !
      endif
      !
      if (subgrid) then
         !
         zmax = subgrid_uv_zmax(ip)
         !
         if (zsu > zmax) then
            hu = subgrid_uv_havg_zmax(ip) + zsu
         else
            dzuv = (zmax - zmin) / (subgrid_nlevels - 1)
            iuv = min(int((zsu - zmin) / dzuv) + 1, subgrid_nlevels - 1)
            facint = (zsu - (zmin + (iuv - 1) * dzuv) ) / dzuv
            hu = subgrid_uv_havg(iuv, ip) + (subgrid_uv_havg(iuv + 1, ip) - subgrid_uv_havg(iuv, ip)) * facint
         endif
         !
      else
         !
         hu = max(zsu - zbuvmx(ip), huthresh)
         !
      endif
      !
      absu  = abs(uv(ip))
      cwave = sqrt(g * hu)
      speed = max(cwave, absu)
      !
      dtip = 1.0 / (speed * dxuvinv)
      !
      if (compute_cell_mins) then
         if (allocated(min_timestep)) then
            min_timestep(nm)  = min(min_timestep(nm),  dtip)
            min_timestep(nmu) = min(min_timestep(nmu), dtip)
         endif
      endif
      !
      if (.not. compute_limiter_details) then
         cycle
      endif
      !
      if (dtip < best_dt_uv) then
         best_dt_uv = dtip
         ip_lim     = ip
         nm_lim     = nm
         nmu_lim    = nmu
         x_lim      = 0.5 * (z_xz(nm) + z_xz(nmu))
         y_lim      = 0.5 * (z_yz(nm) + z_yz(nmu))
         hu_lim     = hu
         abs_uv_lim = absu
         c_wave_lim = cwave
         dx_lim     = 1.0 / dxuvinv
         iref_lim   = iref
         itype_lim  = itype
         idir_lim   = idir
         if (cwave >= absu) then
            reason = 0
         else
            reason = 1
         endif
      endif
      !
   enddo
   !
   if (.not. compute_limiter_details) then
      return
   endif
   !
   if (ip_lim > 0) then
      dt_uv_min = best_dt_uv
   endif
   !
   if (ip_lim <= 0) then
      reason = 2
      return
   endif
   !
   if (best_dt_uv > dtmax + max(1.0e-12, 1.0e-6 * dtmax)) then
      reason  = 2
      ip_lim  = 0
      nm_lim  = 0
      nmu_lim = 0
      x_lim   = 0.0
      y_lim   = 0.0
      hu_lim  = 0.0
      abs_uv_lim = 0.0
      c_wave_lim = 0.0
      dx_lim     = 0.0
      iref_lim   = 0
      itype_lim  = 0
      idir_lim   = 0
   endif
   !
   end subroutine
   !
   subroutine update_timestep_analysis(dt_next)
   !
   real*4, intent(in) :: dt_next
   !
   integer*4 :: nm
   integer*4 :: count_wet
   real*4    :: dt_cell
   real*4    :: dt_tol
   !
   if (.not. timestep_analysis) then
      return
   endif
   !
   if (.not. allocated(min_timestep)) then
      return
   endif
   if (.not. allocated(average_timestep)) then
      return
   endif
   if (.not. allocated(times_wet)) then
      return
   endif
   if (.not. allocated(times_limiting)) then
      return
   endif
   !
   dt_tol = max(1.0e-12, 1.0e-6 * max(dt_next, 1.0e-12))
   !
   do nm = 1, np
      !
      if (min_timestep(nm) > 1.0e29) then
         min_timestep(nm) = 0.0
         cycle
      endif
      !
      dt_cell = alfa * min_timestep(nm)
      !
      count_wet = times_wet(nm) + 1
      times_wet(nm) = count_wet
      !
      average_timestep(nm) = average_timestep(nm) + (dt_cell - average_timestep(nm)) / real(count_wet, kind(average_timestep))
      !
      if (dt_cell <= dt_next + dt_tol) then
         times_limiting(nm) = times_limiting(nm) + 1
      endif
      !
   enddo
   !
   end subroutine
   !
end module sfincs_timestep_diagnostics
