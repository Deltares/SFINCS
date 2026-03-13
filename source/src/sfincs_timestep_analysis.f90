module sfincs_timestep_analysis
   !
   use sfincs_data
   use sfincs_log
   !
   implicit none
   !   
   contains
   !
   subroutine initialize_timestep_analysis()
      ! 
      allocate(timestep_analysis_average_required_timestep(npuv))
      allocate(timestep_analysis_required_timestep(npuv))
      allocate(timestep_analysis_times_limiting(npuv))
      allocate(timestep_analysis_times_wet(npuv))      
      !
      timestep_analysis_average_required_timestep = 0.0
      timestep_analysis_required_timestep = dtmax  
      timestep_analysis_times_wet = 0
      timestep_analysis_times_limiting = 0
      !
   endsubroutine
       
   
   
   subroutine timestep_analysis_update(min_dt)    
      !
      ! Update running statistics immediately after compute_fluxes.
      ! timestep_analysis_required_timestep was initialised to dtmax in compute_fluxes; wet cells
      ! have a computed value < dtmax, dry cells retain dtmax and are skipped.
      !
      real*4, intent(in) :: min_dt
      integer            :: ip
      !
      !$omp parallel &
      !$omp private ( ip )
      !$omp do
      !$acc parallel, present( timestep_analysis_average_required_timestep, timestep_analysis_required_timestep, timestep_analysis_times_wet, timestep_analysis_times_limiting, kcuv)
      !$acc loop gang vector
      !
      do ip = 1, npuv
         !
         if (kcuv(ip) == 1 .or. kcuv(ip) == 6) then
            !
            ! Update running average for wet cells
            !
            timestep_analysis_average_required_timestep(ip) = timestep_analysis_average_required_timestep(ip) + timestep_analysis_required_timestep(ip) 
            !
            ! Note - final conversion to average is done in timestep_analysis_finalize()
            !
            timestep_analysis_times_wet(ip) = timestep_analysis_times_wet(ip) + 1
            !
            ! Check if this cell was limiting the global timestep
            !
            if (timestep_analysis_required_timestep(ip) <= min_dt + 1.0e-6) then            
               ! 
               timestep_analysis_times_limiting(ip) = timestep_analysis_times_limiting(ip) + 1
               !
            endif
            !
         endif         
         !
      enddo
      !$omp end do
      !$omp end parallel
      !$acc end parallel         
      !
   end subroutine timestep_analysis_update

   
   subroutine timestep_analysis_finalize(nt)
   !
   ! Compute average time step per cell and times limiting once at end of simulation
   !
   use sfincs_data
   !
   implicit none
   !
   integer, intent(in) :: nt
   real*4  :: dtm
   integer :: nm, nmd1, nmu1, ndm1, num1, nmd2, nmu2, ndm2, num2
   integer :: tmsl, maximum_times_limiting
   !
   ! Create arrays that will be written to map file
   !
   allocate(timestep_analysis_average_required_timestep_per_cell(np))
   allocate(timestep_analysis_percentage_limiting_per_cell(np))
   !
   timestep_analysis_average_required_timestep_per_cell = 0.0
   timestep_analysis_percentage_limiting_per_cell = 0.0
   !
   ! Copy arrays from GPU to CPU
   !
   !$acc update host(timestep_analysis_average_required_timestep)
   !$acc update host(times_wet)
   !$acc update host(times_limiting)       
   !
   do nm = 1, np
      !
      if (kcs(nm) > 0) then
         !
         ! Compute average required time step for the 4 (or potentially 8!) neighboring U/V points
         !
         nmd1 = z_index_uv_md1(nm)
         nmd2 = z_index_uv_md2(nm)
         nmu1 = z_index_uv_mu1(nm)
         nmu2 = z_index_uv_mu2(nm)
         ndm1 = z_index_uv_nd1(nm)
         ndm2 = z_index_uv_nd2(nm)
         num1 = z_index_uv_nu1(nm)
         num2 = z_index_uv_nu2(nm)
         !
         dtm = 1.0e6
         tmsl = 0
         !
         if (nmd1 > 0) then
            if (kcuv(nmd1) == 1 .or. kcuv(nmd1) == 6) then
               dtm = min(dtm, timestep_analysis_average_required_timestep(nmd1) / max(timestep_analysis_times_wet(nmd1), 1))
               tmsl = max(tmsl, timestep_analysis_times_limiting(nmd1))
            endif
         endif
         !
         if (nmd2 > 0) then
            if (kcuv(nmd2) == 1 .or. kcuv(nmd2) == 6) then
               dtm = min(dtm, timestep_analysis_average_required_timestep(nmd2) / max(timestep_analysis_times_wet(nmd2), 1))
               tmsl = max(tmsl, timestep_analysis_times_limiting(nmd2))
            endif
         endif
         !
         if (nmu1 > 0) then
            if (kcuv(nmu1) == 1 .or. kcuv(nmu1) == 6) then
               dtm = min(dtm, timestep_analysis_average_required_timestep(nmu1) / max(timestep_analysis_times_wet(nmu1), 1))
               tmsl = max(tmsl, timestep_analysis_times_limiting(nmu1))
            endif
         endif
         !
         if (nmu2 > 0) then
            if (kcuv(nmu2) == 1 .or. kcuv(nmu2) == 6) then
               dtm = min(dtm, timestep_analysis_average_required_timestep(nmu2) / max(timestep_analysis_times_wet(nmu2), 1))
               tmsl = max(tmsl, timestep_analysis_times_limiting(nmu2))
            endif
         endif
         !
         if (ndm1 > 0) then
            if (kcuv(ndm1) == 1 .or. kcuv(ndm1) == 6) then
               dtm = min(dtm, timestep_analysis_average_required_timestep(ndm1) / max(timestep_analysis_times_wet(ndm1), 1))
               tmsl = max(tmsl, timestep_analysis_times_limiting(ndm1))
            endif
         endif
         !
         if (ndm2 > 0) then
            if (kcuv(ndm2) == 1 .or. kcuv(ndm2) == 6) then
               dtm = min(dtm, timestep_analysis_average_required_timestep(ndm2) / max(timestep_analysis_times_wet(ndm2), 1))
               tmsl = max(tmsl, timestep_analysis_times_limiting(ndm2))
            endif
         endif
         !
         if (num1 > 0) then
            if (kcuv(num1) == 1 .or. kcuv(num1) == 6) then
               dtm = min(dtm, timestep_analysis_average_required_timestep(num1) / max(timestep_analysis_times_wet(num1), 1))
               tmsl = max(tmsl, timestep_analysis_times_limiting(num1))
            endif
         endif
         !
         if (num2 > 0) then
            if (kcuv(num2) == 1 .or. kcuv(num2) == 6) then
               dtm = min(dtm, timestep_analysis_average_required_timestep(num2) / max(timestep_analysis_times_wet(num2), 1))
               tmsl = max(tmsl, timestep_analysis_times_limiting(num2))
            endif
         endif
         !
         if (dtm < 1.0e5) then
            !
            ! Regular point with at least one normal wet neighboring U/V point
            !
            timestep_analysis_average_required_timestep_per_cell(nm) = min(dtm * alfa, dtmax)
            !
         else
            !
            ! No normal wet neighboring U/V points, likely a dry cell or near boundary, set to -1 to indicate not wet or not valid
            !
            timestep_analysis_average_required_timestep_per_cell(nm) = -1.0
            !
         endif
         !
         timestep_analysis_percentage_limiting_per_cell(nm) = 100.0 * tmsl / nt
         !
      endif
      ! 
   enddo
   !
   end subroutine timestep_analysis_finalize
   
   
   subroutine timestep_analysis_write_log()    
      !      
      integer :: ip, nm, nmu
      real*4 :: xuv, yuv, percentage_limiting, max_times_limiting, max_times_wet
      !      
      if (maxval(timestep_analysis_times_limiting) > 0.0) then
          !
          ip = maxloc(timestep_analysis_times_limiting, dim=1)
          !
          max_times_limiting = maxval(timestep_analysis_times_limiting)
          !
          max_times_wet = maxval(timestep_analysis_times_wet)           
          !
          percentage_limiting = max_times_limiting / max_times_wet * 100.0 ! percentage limiting
          !
      else
          !
          ! No limiting grid cells, all limited by dtmax
          !
          ip = 0    
          !
      endif
      !       
      nm  = uv_index_z_nm(ip)
      nmu = uv_index_z_nmu(ip)
      !
      ! Coordinates of flux link (check if flux link above...)
      !
      xuv = 0.5 * (z_xz(nm) + z_xz(nmu))
      yuv = 0.5 * (z_yz(nm) + z_yz(nmu))
      !
      call write_log('----------- Time step analysis -----------', 1)                  
      call write_log('', 1)
      write(logstr,'(a,20i12)')             ' Index of most time step limiting uv point : ', ip
      call write_log(logstr, 1) 
      write(logstr,'(a,f10.1,a)')           ' Percentage of time point was limiting     : ', percentage_limiting, ' %'
      call write_log(logstr, 1)           
      write(logstr,'(a,20f12.3)')           ' x coordinate                              : ', xuv
      call write_log(logstr, 1) 
      write(logstr,'(a,20f12.3)')           ' y coordinate                              : ', yuv
      call write_log(logstr, 1)
      write(logstr,'(a)') ''
      call write_log(logstr, 1)
      !        
   end subroutine timestep_analysis_write_log
   !
end module sfincs_timestep_analysis
