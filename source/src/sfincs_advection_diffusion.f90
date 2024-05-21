module sfincs_advection_diffusion
   
contains

   subroutine compute_tracer_fluxes(dt, min_dt, tloop)
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
   integer   :: itracer
   integer   :: icuv
   integer   :: nm
   integer   :: nmu
   integer   :: n
   integer   :: m
   !
   call system_clock(count0, count_rate, count_max)
   !
   ! Update fluxes
   !
   !$omp parallel &
   !$omp private ( ip, itracer, icuv, nm, nmu) &
   !$omp do schedule ( dynamic, 256 )
   !$acc kernels, present( kcuv, q, trconc, trflux, uv_index_z_nm, uv_index_z_nmu )
   !$acc loop independent
   do ip = 1, npuv
      !
      if (kcuv(ip) == 1) then
         !
         ! Indices of surrounding water level points
         !
         nm  = uv_index_z_nm(ip)
         nmu = uv_index_z_nmu(ip)
         !
		 do itracer = 1, ntracer
		    ! 
            if (q(ip) > 0.0) then
   	           ! 
               trflux(itracer, ip) = q(ip) * trconc(itracer, nm) + (trconc(itracer, nm) - trconc(itracer, nmu)) * dico
		       ! 
            else
   	           ! 
               trflux(itracer, ip) = q(ip) * trconc(itracer, nmu) + (trconc(itracer, nm) - trconc(itracer, nmu)) * dico
   	           ! 
            endif
		    ! 
		 enddo 
         !
	  endif
   enddo	  
   !$omp end do
   !$omp end parallel
   !$acc end kernels
   !
   if (ncuv > 0) then
      !
      ! Loop through combined uv points and determine average uv and q
      ! The combined q and uv values are used in the continuity equation and in the netcdf output
      !
      !$omp parallel &
      !$omp private ( icuv, itracer )
      !$omp do
      !$acc kernels, present( trflux, cuv_index_uv, cuv_index_uv1, cuv_index_uv2 )
      !$acc loop independent
      do icuv = 1, ncuv
         do itracer = 1, ntracer
            !
            ! Average of the two uv points
            !
            trflux(itracer, cuv_index_uv(icuv))  = (itracer, trflux(cuv_index_uv1(icuv)) + trflux(itracer, cuv_index_uv2(icuv))) / 2
            !
         enddo
      enddo	  
      !$acc end kernels
      !$omp end do
      !$omp end parallel
      !
   endif
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine      
   !
end module
