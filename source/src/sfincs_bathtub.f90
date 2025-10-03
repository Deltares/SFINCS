module sfincs_bathtub

integer*4, dimension(:),   allocatable :: bathtub_i1
integer*4, dimension(:),   allocatable :: bathtub_i2
real*4,    dimension(:),   allocatable :: bathtub_w1
integer*4                              :: bathtub_snapwave_nwbnd
real*4, allocatable  :: bathtub_snapwave_x_bwv(:)
real*4, allocatable  :: bathtub_snapwave_y_bwv(:)
real*4, allocatable  :: bathtub_snapwave_hs_bwv(:)
integer, allocatable :: bathtub_snapwave_i1(:)
integer, allocatable :: bathtub_snapwave_i2(:)
real*4, allocatable  :: bathtub_snapwave_w1(:)

contains
   !
   subroutine initialize_bathtub()
   !
   ! Creates bathtub boundary conditions and determines weighting factors
   !
   use sfincs_data
   use geometry
   use quadtree
   use sfincs_error
   use sfincs_log
   !
   implicit none
   !
   integer   :: ip, ib, itb
   integer   :: i1, i2
   real*4    :: w1, w2
   real*4    :: h1, h2
   real*8    :: t8
   !
   if (bathtub_snapwave) then
      !
      ! First read wave boundary conditions. Same as in snapwave.
      ! We cannot use snapwave_data and snapwave_boundaries unfortunately, as these share common variable names.
      !
      call read_snapwave_boundary_data() ! snapwave_boundaries
      !
      ! Determine weights and indices for each wave boundary point
      !
      allocate(bathtub_snapwave_i1(nbnd))
      allocate(bathtub_snapwave_i2(nbnd))
      allocate(bathtub_snapwave_w1(nbnd))
      !
      ! Loop along boundary points
      !
      do ib = 1, nbnd
         !
         call interp_segment(bathtub_snapwave_x_bwv, bathtub_snapwave_y_bwv, bathtub_snapwave_nwbnd, x_bnd(ib), y_bnd(ib), i1, i2, w1, w2)
         !
         bathtub_snapwave_i1(ib) = i1
         bathtub_snapwave_i2(ib) = i2
         bathtub_snapwave_w1(ib) = w1      
         !
      enddo
      !
      ! Loop through times in bzs file and add f*Hm0 (typically f=0.2)
      !
      do itb = 1, ntbnd
         !
         ! Need to convert t_bnd to real*8
         !
         t8 = t_bnd(itb)
         !
         call update_snapwave_boundary_data(t8)
         !
         ! Loop through boundary points
         !
         do ib = 1, nbnd            
            !
            h1 = bathtub_snapwave_hs_bwv(bathtub_snapwave_i1(ib))
            h2 = bathtub_snapwave_hs_bwv(bathtub_snapwave_i2(ib))
            w1 = bathtub_snapwave_w1(ib)
            w2 = 1.0 - w1
            !
            zs_bnd(ib, itb) = zs_bnd(ib, itb) + bathtub_fac_hs * (w1 * h1 + w2 * h2)
            !
         enddo   
         !
      enddo   
      !
   endif   
   !
   ! Determine indices and weight for each grid point
   !
   allocate(bathtub_i1(np))
   allocate(bathtub_i2(np))
   allocate(bathtub_w1(np))
   !
   do ip = 1, np
      !
      call interp_segment(x_bnd, y_bnd, nbnd, z_xz(ip), z_yz(ip), i1, i2, w1, w2)
      !
      bathtub_i1(ip) = i1
      bathtub_i2(ip) = i2
      bathtub_w1(ip) = w1
      !
   enddo
   !
   meteo3d = .false.
   wind = .false.
   store_meteo = .false.
   store_wind = .false.
   store_wind_max = .false.
   precip = .false.
   patmos = .false.
   snapwave = .false.
   infiltration = .false.
   store_velocity = .false.
   store_maximum_velocity = .false.
   !
   end subroutine

   
   subroutine bathtub_compute_water_levels(tloop)
   !
   use sfincs_data
   use geometry
   !
   implicit none
   !
   integer  :: count0
   integer  :: count1
   integer  :: count_rate
   integer  :: count_max
   real     :: tloop
   !
   integer :: nm, i1, i2
   real*4  :: zbt, w1, w2
   !
   call system_clock(count0, count_rate, count_max)
   !
   !$omp parallel &
   !$omp private ( nm, i1, i2, w1, w2 )
   !$omp do schedule ( dynamic, 256 )
   do nm = 1, np
      !
      i1 = bathtub_i1(nm)
      i2 = bathtub_i2(nm)
      w1 = bathtub_w1(nm)
      w2 = 1.0 - w1
      !
      zbt = w1 * zst_bnd(i1) + w2 * zst_bnd(i2)
      !
      if (subgrid) then
         !
         zs(nm) = max(subgrid_z_zmin(nm), zbt)
         !
      else
         !
         zs(nm) = max(zb(nm), zbt)
         !
      endif
      !
      if (store_maximum_waterlevel) then
         !
         ! Store the maximum water level
         !
         zsmax(nm) = max(zsmax(nm), zs(nm))
         !
      endif      
      !
   enddo   
   !$omp end do
   !$omp end parallel
   !
   !$acc update device( zs, zsmax )
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0 * (count1 - count0) / count_rate
   !
   end subroutine

   subroutine read_snapwave_boundary_data()
   !
   use snapwave_data
   use snapwave_boundaries
   use sfincs_snapwave
   !
   ! Read SnapWave input file from sfincs.inp (this will store boundary file names in snapwave_data)
   !
   call read_snapwave_input()
   !
   ! Use read_boundary data from snapwave_boundaries
   !
   call read_boundary_data()   
   !
   allocate(bathtub_snapwave_x_bwv(nwbnd))
   allocate(bathtub_snapwave_y_bwv(nwbnd))
   allocate(bathtub_snapwave_hs_bwv(nwbnd))
   !
   bathtub_snapwave_nwbnd = nwbnd
   bathtub_snapwave_x_bwv = x_bwv
   bathtub_snapwave_y_bwv = y_bwv
   !
   end subroutine

   
   subroutine update_snapwave_boundary_data(tb)
   !
   use snapwave_data
   use snapwave_boundaries   
   !
   implicit none
   !
   real*8 :: tb
   !
   call update_boundary_points(tb, .true.)
   !
   bathtub_snapwave_hs_bwv = hst_bwv
   !
   end subroutine
   
end module
