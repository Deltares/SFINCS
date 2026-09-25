module sfincs_snapwave
   !
   use sfincs_log
   use sfincs_error
   use snapwave_input, only: read_snapwave_input
   !    
   implicit none
   !     
   integer                                   :: snapwave_no_nodes
   integer                                   :: snapwave_no_cells
   real*8,    dimension(:),   allocatable    :: snapwave_x
   real*8,    dimension(:),   allocatable    :: snapwave_y
   real*4,    dimension(:),   allocatable    :: snapwave_z
   real*4,    dimension(:),   allocatable    :: snapwave_mask   
   real*4,    dimension(:),   allocatable    :: snapwave_depth
   real*4,    dimension(:),   allocatable    :: snapwave_H
   real*4,    dimension(:),   allocatable    :: snapwave_H_ig
   real*4,    dimension(:),   allocatable    :: snapwave_Tp
   real*4,    dimension(:),   allocatable    :: snapwave_Tp_ig   
   real*4,    dimension(:),   allocatable    :: snapwave_mean_direction
   real*4,    dimension(:),   allocatable    :: snapwave_u10
   real*4,    dimension(:),   allocatable    :: snapwave_u10dir   
   real*4,    dimension(:),   allocatable    :: snapwave_Fx
   real*4,    dimension(:),   allocatable    :: snapwave_Fy
   real*4,    dimension(:),   allocatable    :: snapwave_Dw
   real*4,    dimension(:),   allocatable    :: snapwave_Df 
   real*4,    dimension(:),   allocatable    :: snapwave_Dwig
   real*4,    dimension(:),   allocatable    :: snapwave_Dfig
   real*4,    dimension(:),   allocatable    :: snapwave_cg
   real*4,    dimension(:),   allocatable    :: snapwave_beta
   real*4,    dimension(:),   allocatable    :: snapwave_srcig
   real*4,    dimension(:),   allocatable    :: snapwave_alphaig   
   integer,   dimension(:,:), allocatable    :: snapwave_connected_nodes
   integer*4, dimension(:),   allocatable    :: index_snapwave_in_sfincs
   integer*4, dimension(:),   allocatable    :: index_sfincs_in_snapwave
   integer*4, dimension(:),   allocatable    :: index_sw_in_qt ! used in sfincs_ncoutput (copy of index_snapwave_in_quadtree from snapwave_data)
   real*4                                    :: snapwave_hsmean
   real*4                                    :: snapwave_tpmean
   real*4                                    :: snapwave_tpigmean   
   real*4                                    :: snapwave_fwmaxfac   
   !
contains
   !
   subroutine couple_snapwave(crsgeo)
   !
   use snapwave_data
   use snapwave_domain
   use snapwave_boundaries
   !
   implicit none
   !
   logical       :: crsgeo
   !
   build_revision = '$Rev: git SFINCS_SnapWave:main' 
   build_date     = '$Date: 2026-06-10'
   !
   call write_log('', 1)
   call write_log('----------- Welcome to SnapWave ---------', 1)
   call write_log('', 1)
   call write_log('   @@@@@   @@  @@  @@@@@@  @@@@@@   @@@  ', 1)
   call write_log('  @@@ @@@  @@@ @@  @@@@@@  @@@@@@   @@@  ', 1)
   call write_log('  @@@      @@@ @@  @@  @@  @@  @@   @@@  ', 1)
   call write_log('   @@@@@   @@@@@@  @@@@@@  @@@@@@   @@@  ', 1)
   call write_log('      @@@  @@ @@@  @@  @@  @@            ', 1)
   call write_log('  @@@ @@@  @@  @@  @@  @@  @@       @@@  ', 1)
   call write_log('   @@@@@   @@   @  @@  @@  @@       @@@  ', 1)
   call write_log('', 1)
   call write_log('             .......:.......             ', 1)
   call write_log('         ...:::::::::::::::::...         ', 1)
   call write_log('      ..:::::::............::::::..      ', 1)
   call write_log('    ..::::::.....:@@@@@@@@....:::::..    ', 1)
   call write_log('   .::::::...~@@@@@@@@@@@@@@~..::::::.   ', 1)
   call write_log('  .::::::..:@@@@@@@@@@@@@@@@@@:.::::::.  ', 1)
   call write_log(' .:::::..:@@@@@@@@@@@@@@@@@@@@@:.::::::. ', 1)
   call write_log('.::::..:@@@@@@@@@@@@@@^......:@@.:::::::.', 1)
   call write_log('.::...:@@@@@@@@@@@@@@@.:::::..^^.:::::::.', 1)
   call write_log('::.:@@@@@@@@@@@@@@@@@@..::::::..:::::::::', 1)
   call write_log('..:@@@@@@@@@@@@@@@@@@@@^..............::.', 1)
   call write_log('..:@@@@@@@@@@@@@@@@@@@@@@@^:..:~^~^~:..:.', 1)
   call write_log(' .:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@:. ', 1)
   call write_log('  .@@~^~@@@@~^~@@@@~^~@@@@~^~@@@@~^~@@.  ', 1)
   call write_log('   ...................................   ', 1)
   call write_log('    ..:::::::::::::::::::::::::::::..    ', 1)
   call write_log('      ..:::::::::::::::::::::::::..      ', 1)
   call write_log('         ...:::::::::::::::::...         ', 1)
   call write_log('             .......:.......             ', 1)
   call write_log('', 1)
   call write_log('-----------------------------------------', 1)   
   call write_log('', 1)
   call write_log('Build-Revision: '//trim(build_revision), 1)
   call write_log('Build-Date: '//trim(build_date), 1)   
   call write_log('', 1)
   call write_log('------ Preparing model simulation --------', 1)
   call write_log('', 1)   
   !   
   ! Check whether SFINCS grid is spherical (T) or cartesian (F), and prescribe to SnapWave as variable 'sferic' -  spherical (1) or cartesian (0) grid
   if (crsgeo) then
      sferic  = 1 
      write(logstr,*)'SnapWave: Input grid interpreted as spherical coordinates, sferic= ',sferic     
      call write_log(logstr, 0)   
   endif   
   !
   call read_snapwave_input()            ! Reads snapwave.inp
   !
   call initialize_snapwave_domain()     ! Read mesh, finds upwind neighbors, etc.
   !
   call read_boundary_data()
   !
   !call find_boundary_indices() ! > is already called in read_boundary_data()
   !
   call write_log('', 1)
   !
   snapwave_no_nodes = no_nodes
   !
   allocate(snapwave_z(no_nodes))
   allocate(snapwave_depth(no_nodes))
   allocate(snapwave_mask(no_nodes))  
   allocate(snapwave_u10(no_nodes))   
   allocate(snapwave_u10dir(no_nodes))      
   !
   snapwave_z     = zb
   snapwave_depth = 0.0
   !
   if (wind) then
      snapwave_u10 = 0.0
      snapwave_u10dir = 0.0
   endif
   !
   snapwave_hsmean = 0.0
   snapwave_tpmean = 0.0
   snapwave_tpigmean = 0.0   
   !   
   call find_matching_cells(index_quadtree_in_snapwave, index_snapwave_in_quadtree)
   !
   ! Copy final snapwave mask from snapwave_domain for output in sfincs_ncoutput
   !
   snapwave_mask = msk   
   !
   call write_log('------------------------------------------', 1)
   call write_log('SnapWave Processes', 1)
   call write_log('------------------------------------------', 1)
   if (igwaves) then
      call write_log('SnapWave IG waves                  : yes', 1)
   else   
      call write_log('SnapWave IG waves                  : no', 1)
   endif
   if (iterative_srcig) then
      call write_log('SnapWave implicit IG source term   : yes', 1)
   else   
      call write_log('SnapWave implicit IG source term   : no', 1)
   endif
   if (igherbers) then 
      call write_log('SnapWave IG bc using Herbers       : yes', 1)
   else   
      call write_log('SnapWave IG bc using Herbers       : no', 1)
   endif
   if (wind) then
      call write_log('SnapWave wind growth               : yes', 1)
   else   
      call write_log('SnapWave wind growth               : no', 1)
   endif
   if (vegetation) then
      call write_log('SnapWave vegetation                : yes', 1)
   else   
      call write_log('SnapWave vegetation                : no', 1)
   endif
   !
   call write_log('------------------------------------------', 1)
   !
   end subroutine
   

   subroutine find_matching_cells(index_quadtree_in_snapwave, index_snapwave_in_quadtree)
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   integer, dimension(snapwave_no_nodes),  intent(in) :: index_quadtree_in_snapwave
   integer, dimension(quadtree_nr_points), intent(in) :: index_snapwave_in_quadtree
   !
   integer :: ipsw, ipsf, iq, ip, counter
   !
   real*4  :: xsw, ysw, dstmin, dst, min_distance 
   !
   real*4 :: distances(np)
   integer :: closest_index(1)
   !
   logical :: nearest_warning
   !   
   allocate(index_sfincs_in_snapwave(snapwave_no_nodes))
   allocate(index_snapwave_in_sfincs(np))
   allocate(index_sw_in_qt(quadtree_nr_points))
   !
   nearest_warning = .false.
   !
   index_sfincs_in_snapwave = 0
   index_snapwave_in_sfincs = 0
   index_sw_in_qt = 0
   counter = 0
   distances = 0.0
   min_distance = 0.0
   !
   ! Loop through SnapWave points
   !
   do ipsw = 1, snapwave_no_nodes
      !
      iq   = index_quadtree_in_snapwave(ipsw)
      ipsf = index_sfincs_in_quadtree(iq)
      !
      if (ipsf == 0 ) then
         !
         ! SFINCS not active at this SnapWave node, so find the nearest SFINCS point
         !
         counter = counter + 1
         !
         nearest_warning = .true. ! to print warning to screen that 'extrapolation' is performed
         !
         if (snapwave_use_nearest) then 
             !
             xsw = quadtree_xz(iq)
             ysw = quadtree_yz(iq)
             !
             dstmin = 1.0e6
             !
             ! Calculate the distance for each coordinate
             !$omp parallel &
             !$omp private ( ip, dst )
             !$omp do
             do ip = 1, np
                 !
                 dst = sqrt((z_xz(ip) - xsw)**2 + (z_yz(ip) - ysw)**2)
                 !
                 distances(ip) = dst
                 !
             enddo
             !$omp end do         
             !$omp end parallel             
             !
             ! Find the minimum distance
             min_distance = minval(distances)         
             !
             if (min_distance < dstmin) then
                 !
                 ! Find the index of the minimum distance
                 closest_index = minloc(distances)
                 !
                 ! To conform shapes
                 ipsf = closest_index(1)
                 !
             endif       
             !        
         endif         
      endif
      !         
      index_sfincs_in_snapwave(ipsw) = ipsf
      !
      index_sw_in_qt(iq) = ipsw      
      !
   enddo             
   !
   ! Loop through SFINCS points
   !
   do ipsf = 1, np
      !
      iq   = index_quadtree_in_sfincs(ipsf)
      ipsw = index_snapwave_in_quadtree(iq)
      index_snapwave_in_sfincs(ipsf) = ipsw
      !
   enddo   
   !
   ! Print warning message
   !
   if (nearest_warning) then
      if (snapwave_use_nearest) then
          write(logstr,'(a,i0,a)')'SnapWave: Info   : ',counter,' SnapWave node(s) do not have a matching SFINCS point, so water depth and wind conditions from the nearest SFINCS point within 1000 km are used for SnapWave calculation '
      else
          write(logstr,'(a,i0,a)')'SnapWave: Info   : ',counter,' SnapWave node(s) do not have a matching SFINCS point, water level at these points is set to 0.0 '          
      endif      
      ! 
      call write_log(logstr, 0)
      !
   endif   
   !
   end subroutine

   
   subroutine update_wave_field(t)
   !
   use sfincs_data
   use sfincs_timers
   use omp_lib
   !
   implicit none
   !
   real*4   :: u10, u10dir
   !
   real*4,    dimension(:), allocatable       :: fwx0
   real*4,    dimension(:), allocatable       :: fwy0
   integer   :: ip, nm, nmu, idir
   real*8    :: t
   real(8)   :: t3, t4
   !
   t3 = omp_get_wtime()
   !
   call timer_start('SnapWave')
   !
   allocate(fwx0(np))
   allocate(fwy0(np))
   !
   fwx0 = 0.0
   fwy0 = 0.0
   !
   ! Determine SnapWave water depth
   !
   do nm = 1, snapwave_no_nodes
      !
      ip = index_sfincs_in_snapwave(nm) ! matching index in SFINCS mesh
      !
      if (ip > 0) then
         !
         ! A matching SFINCS point is found
         !

         if (wavemaker) then
            !
            snapwave_depth(nm) = max(zsm(ip) - snapwave_z(nm), 0.00001)
            !
         else   
            !
            snapwave_depth(nm) = max(zs(ip) - snapwave_z(nm), 0.00001)      
            !
         endif   
         !
      else
         !
         ! Use 0.0 water level
         !
         snapwave_depth(nm) = max(0.0 - snapwave_z(nm), 0.00001)      
         !
      endif   
      !
   enddo   
   !
   ! Determine SnapWave wind
   !
   if (wind) then ! =We have wind inputs given to SFINCS
      !
      if (snapwavewind) then ! =We have windgrowth in SnapWave turned on 
          !
          do nm = 1, snapwave_no_nodes
             !
             ip = index_sfincs_in_snapwave(nm) ! matching index in SFINCS mesh
             !
             if (ip>0) then
                !
                ! A matching SFINCS point is found
                !
                ! Convert to umag & dir, as in ncoutput_update_his: 
                !
                u10 = sqrt(windu(ip)**2 + windv(ip)**2)
                !
                u10dir = atan2(windv(ip), windu(ip))*180/pi
                !
	            if (u10dir<0.0) u10dir = u10dir + 360.0
                if (u10dir>360.0) u10dir = u10dir - 360.0    
                !
                snapwave_u10(nm) = max(u10, 0.0)     
                snapwave_u10dir(nm) = u10dir / 180.0 * pi ! from nautical coming from in degrees to cartesian going to in radians
                !
             else
                !
                ! Use 0.0 wind speed and direction
                !
                snapwave_u10(nm) = 0.0
                snapwave_u10dir(nm) = 0.0            
                !
             endif   
             !
          enddo   
          !
      endif
      !
   endif   
   !
   call compute_snapwave(t)
   !
   do nm = 1, np
      !
      ip = index_snapwave_in_sfincs(nm) ! matching index in SFINCS mesh
      !
      if (ip>0) then
         !
         hm0(nm)    = snapwave_H(ip)
         hm0_ig(nm) = snapwave_H_ig(ip)
         fwx0(nm)   = snapwave_Fx(ip)
         fwy0(nm)   = snapwave_Fy(ip)
         !
      else
         !
         ! SnapWave point outside active SFINCS domain
         !
         hm0(nm)    = 0.0
         hm0_ig(nm) = 0.0
         fwx0(nm)   = 0.0
         fwy0(nm)   = 0.0
         !
      endif
      !
   enddo
   !
   hm0 = hm0 * sqrt(2.0)
   hm0_ig = hm0_ig * sqrt(2.0)
   !
   do ip = 1, npuv
      !
      nm   = uv_index_z_nm(ip)
      nmu  = uv_index_z_nmu(ip)
      idir = uv_flags_dir(ip) ! 0 is u, 1 is v
      !
      ! Should do better averaging for uv points that go from fine to coarse
      !
      if (idir == 0) then
         !
         ! U point
         !         
         fwuv(ip) = waveforces_ratio * (0.5 * (cosrot * fwx0(nm) + sinrot * fwy0(nm)) + 0.5 * ( cosrot * fwx0(nmu) + sinrot * fwy0(nmu))) / rhow
         ! waveforces_ratio = 1.0 by default, but can be set to 0 to avoid double counting incident setup if wavemaker_hinc true
      else
         !
         ! V point
         !         
         fwuv(ip) = waveforces_ratio * (0.5 * (-sinrot * fwx0(nm) + cosrot * fwy0(nm)) + 0.5 * (-sinrot * fwx0(nmu) + cosrot * fwy0(nmu))) / rhow
         !         
      endif   
      !
   enddo
   !
   !$acc update device(fwuv)
   !
   ! Set wave forces fwmaxfac factor
   fwmaxfac = snapwave_fwmaxfac
   !
   call timer_stop('SnapWave')
   !
   t4 = omp_get_wtime()
   !
   write(logstr,'(a,f10.1,a,f6.2,a)')'Computing SnapWave at t = ', t, ' s took ', t4 - t3, ' seconds'
   call write_log(logstr, 0)
   !
   end subroutine


   subroutine compute_snapwave(t)
   !
   use snapwave_data
   use snapwave_solver
   use snapwave_boundaries
   !
   real*8    :: t
   integer   :: k
   ! 
   depth = snapwave_depth
   !
   zb = snapwave_z   
   !
   u10 = snapwave_u10
   u10dir = snapwave_u10dir   
   !   
   ! TL: we use depth now in boundary conditions for Herbers bc determination of Hm0ig, in this order we use updated values of depth through SFINCS
   !
   call update_boundary_conditions(t) ! SnapWave boundary conditions
   !
   call compute_wave_field()
   !
   snapwave_H                     = H
   snapwave_H_ig                  = H_ig
   snapwave_Tp                    = Tp
   snapwave_Tp_ig                 = Tp_ig   
   snapwave_mean_direction        = modulo(270.0 - thetam * 180 / pi + 360.0, 360.0)
   snapwave_Dw                    = Dw
   snapwave_Df                    = Df
   snapwave_Dwig                  = Dw_ig
   snapwave_Dfig                  = Df_ig
   snapwave_cg                    = cg
   snapwave_beta                  = beta
   snapwave_srcig                 = srcig
   snapwave_alphaig               = alphaig   
   !   
   ! Convert wave force to correct unit [Dw/C] as expected by SFINCS, assumed to be piecewise (seems to work)
   snapwave_Fx                    = Fx * rho * depth
   snapwave_Fy                    = Fy * rho * depth
   !
   ! Pre-alculate wave forces limiter factor
   snapwave_fwmaxfac = 0.25 * sqrt(g) * rho * gammax**2 / tpmean_bwv    
   !
   ! FIXME - should we limit snapwave_fwmaxfac to a certain range?
   !
   ! Loop over points and set Tp, cg, direction, spreading to 0 where H and/or H_ig are zero
   ! TL: needed because e.g. Tp is set to Tpini initially, so shows values even if cell remains dry with H=0
   do k = 1, no_nodes
       if (snapwave_H(k) <= 0.0) then
           snapwave_Tp(k) = 0.0
           snapwave_mean_direction(k) = 0.0
           snapwave_cg(k) = 0.0
       endif
       !
       if (snapwave_H_ig(k) <= 0.0) then
           snapwave_Tp_ig(k) = 0.0        
       endif       
   enddo   
   !
   ! Wave periods from SnapWave, used in e.g. wavemakers - TL: moved behind call update_boundary_conditions & compute_wave_field so values at first timestep are not 0
   !
   snapwave_hsmean = hsmean_bwv
   snapwave_tpmean = tpmean_bwv
   !
   ! Do quick check whether incoming Tpig value seems realistic, before using it:
   if (igwaves) then
       !
       snapwave_tpigmean = tpmean_bwv_ig      
       ! 
   endif
   !
   end subroutine

   
end module
