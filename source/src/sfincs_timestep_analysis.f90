module sfincs_timestep_analysis
   !
   ! Modules
   use sfincs_data
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
      ! min_timestep was initialised to dtmax in compute_fluxes; wet cells
      ! have a computed value < dtmax, dry cells retain dtmax and are skipped.
      !
      real*4, intent(in) :: min_dt
      integer            :: nm
      !
      !$omp parallel &
      !$omp private ( ip )
      !$omp do
      !$acc parallel, present( timestep_analysis_required_timestep, timestep_analysis_required_timestep, timestep_analysis_times_wet, timestep_analysis_times_limiting, kcuv)
      !$acc loop gang vector
      !
      do ip = 1, npuv
         !
         if (kcuv(ip) == 1 .or. kcuv(ip) == 6) then
             !
             if (min_timestep(ip) < dtmax) then
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
   !
end module sfincs_timestep_analysis
