module sfincs_timestep_diag
   !
   ! Modules
   use sfincs_data
   !
   contains
   !
   subroutine compute_cell_min_timestep(nm, hu, uv_ip, dxuvinv)
      !$acc routine seq
      !
      ! Store per-cell limiting timestep (called from inside momentum parallel loop)
      !
      integer, intent(in) :: nm
      real*4,  intent(in) :: hu, uv_ip, dxuvinv
      !
      min_timestep(nm) = 1.0 / (max(sqrt(g * hu), abs(uv_ip)) * dxuvinv)
      !
   end subroutine compute_cell_min_timestep

   
   
   subroutine timestep_diagnostics_update(dt)
      !
      ! Update running statistics immediately after compute_fluxes.
      ! min_timestep was initialised to dtmax in compute_fluxes; wet cells
      ! have a computed value < dtmax, dry cells retain dtmax and are skipped.
      !
      real*4, intent(in) :: dt
      integer            :: nm
      !
      do nm = 1, np
         !
         if (min_timestep(nm) < dtmax) then
            ! 
            ! Update running average for wet cells
            !
            average_timestep(nm) = (average_timestep(nm) * real(times_wet(nm), kind(min_timestep)) + alfa * min_timestep(nm)) / real(times_wet(nm) + 1, kind(min_timestep))
            !
            times_wet(nm) = times_wet(nm) + 1
            !
            ! Check if this cell was limiting the global timestep
            !
            if (alfa * min_timestep(nm) <= dt) then
               ! 
               times_limiting(nm) = times_limiting(nm) + 1
               !
            endif
            !
         endif
         !
      enddo
      !
   end subroutine timestep_diagnostics_update
   !
end module sfincs_timestep_diag
