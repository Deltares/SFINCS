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
             if (timestep_analysis_required_timestep(ip) < dtmax) then
                ! 
                ! Update running average for wet cells
                !
                timestep_analysis_average_required_timestep(ip) = timestep_analysis_average_required_timestep(ip) + timestep_analysis_required_timestep(ip) 
                ! Note - final conversion to average is done in ncoutput_write_timestep_analysis()
                !
                timestep_analysis_times_wet(ip) = timestep_analysis_times_wet(ip) + 1
                !
                ! Check if this cell was limiting the global timestep
                !
                if (timestep_analysis_required_timestep(ip) <= (min_dt+1.0e-6)) then            
                   ! 
                   timestep_analysis_times_limiting(ip) = timestep_analysis_times_limiting(ip) + 1
                   !
                endif
                !
             endif
             !
         endif         
      enddo
      !$omp end do
      !$omp end parallel
      !$acc end parallel         
      !
   end subroutine timestep_analysis_update
   
   
   
   subroutine finalize_timestep_analysis()    
      !      
      integer :: ip, nm, nmu
      real*4 :: xuv, yuv, percentage_limiting, max_times_limiting, max_times_wet
      !      
      if (maxval(timestep_analysis_times_limiting) > 0.0) then
          !
          ip      = maxloc(timestep_analysis_times_limiting, dim=1)
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
      write(logstr,'(a,f10.3,a)')           ' Max timestep limiting  : ', percentage_limiting, ' % '
      call write_log(logstr, 1)           
      !write(logstr,'(a,i,a,f10.3,a,f10.3)')           ' uv : ', ip, ' x : ', xuv, ' y : ', yuv
      !call write_log(logstr, 1)
      write(logstr,'(a,20i10)')               '    uv  : ', ip
      call write_log(logstr, 1) 
      write(logstr,'(a,20f10.3)')           '    x   : ', xuv
      call write_log(logstr, 1) 
      write(logstr,'(a,20f10.3)')           '    y   : ', yuv
      call write_log(logstr, 1)          
      !        
   end subroutine finalize_timestep_analysis
   !
end module sfincs_timestep_analysis
