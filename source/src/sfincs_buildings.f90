!! This module implements building polygon routines that:
!! 1. Collect rainfall on building footprint
!! 2. Block flow with thin dams on edges
!! 3. Redistribute collected water to perimeter cells OR gutter outlets
!! 4. Optional detention volume per building
!! 5. Optional gutter system with configurable spacing
!! IMDC, 2025

module sfincs_buildings

   use sfincs_data  ! Import building data structures and grid parameters

   implicit none

   ! Module contains only subroutines, no data
   ! All data is in sfincs_data module

contains

   ! ========================================================================
   ! Initialize building module
   ! Sets up grid parameters and allocates masks
   ! ========================================================================
   subroutine initialize_buildings()

      implicit none
      integer :: ncells

      ! Grid parameters are already in sfincs_data
      ! Just allocate the grid masks
      ncells = mmax*nmax

      if (.not. allocated(is_building_cell)) then
         allocate (is_building_cell(ncells))
         allocate (is_perimeter_cell(ncells))
         allocate (is_gutter_cell(ncells))
         allocate (building_id(ncells))
      end if

      ! Initialize
      is_building_cell = .false.
      is_perimeter_cell = .false.
      is_gutter_cell = .false.
      building_id = 0

      write (*, '(A)') 'Building module initialized'

   end subroutine initialize_buildings

   ! ========================================================================
   !  Read building polygons from file
   ! ========================================================================
   subroutine read_building_file(filename, ios_out)

      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ios_out

      integer :: iunit, ios_bld, i, j, ibld
      character(len=1024) :: line
      character(len=256) :: name
      integer :: npoints, ncols
      real*8 :: xtmp, ytmp

      ios_out = 0

      ! Open file
      open (newunit=iunit, file=trim(filename), status='old', action='read', iostat=ios_bld)
      if (ios_bld /= 0) then
         write (*, '(A,A)') 'ERROR: Could not open building file: ', trim(filename)
         ios_out = ios_bld
         return
      end if

      write (*, '(A,A)') 'Reading building file: ', trim(filename)

      ! Count buildings
      nbuildings = 0
      do
         read (iunit, '(A)', iostat=ios_bld) line
         if (ios_bld /= 0) exit
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '!') cycle

         read (line, *, iostat=ios_bld) name, npoints, ncols
         if (ios_bld == 0 .and. npoints > 0 .and. ncols >= 2) then
            nbuildings = nbuildings + 1
            do i = 1, npoints
               read (iunit, *, iostat=ios_bld)
               if (ios_bld /= 0) exit
            end do
         end if
      end do

      if (nbuildings == 0) then
         write (*, '(A)') 'WARNING: No buildings found in file'
         close (iunit)
         return
      end if

      write (*, '(A,I0,A)') 'Found ', nbuildings, ' buildings'

      ! Allocate buildings array (in sfincs_data)
	  
      allocate (buildings(nbuildings))

      ! Read building data
      rewind (iunit)

      ibld = 0
      do
         read (iunit, '(A)', iostat=ios_bld) line
         if (ios_bld /= 0) exit
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '!') cycle

         read (line, *, iostat=ios_bld) name, npoints, ncols
         if (ios_bld /= 0 .or. npoints <= 0 .or. ncols < 2) cycle

         ibld = ibld + 1
         if (ibld > nbuildings) exit

         buildings(ibld)%name = trim(name)
         buildings(ibld)%npoints = npoints

         allocate (buildings(ibld)%x(npoints))
         allocate (buildings(ibld)%y(npoints))

         do i = 1, npoints
            if (ncols == 2) then
               read (iunit, *, iostat=ios_bld) xtmp, ytmp
            else
               read (iunit, *, iostat=ios_bld) xtmp, ytmp
            end if

            if (ios_bld /= 0) then
               write (*, '(A,A)') 'ERROR: Reading coordinates for building ', trim(name)
               ios_out = ios_bld
               close (iunit)
               return
            end if

            buildings(ibld)%x(i) = xtmp
            buildings(ibld)%y(i) = ytmp
         end do

         if (buildings(ibld)%x(1) /= buildings(ibld)%x(npoints) .or. &
             buildings(ibld)%y(1) /= buildings(ibld)%y(npoints)) then
            write (*, '(A,A,A)') 'WARNING: Building ', trim(name), ' polygon not closed'
         end if

         ! Initialize all building properties
         buildings(ibld)%accumulated_rainfall = 0.0d0
         buildings(ibld)%accumulated_runoff = 0.0d0
         buildings(ibld)%ncells_inside = 0
         buildings(ibld)%ncells_perimeter = 0
         buildings(ibld)%ncells_gutter = 0
         buildings(ibld)%total_area = 0.0d0
         buildings(ibld)%perimeter_length = 0.0d0
         
         ! Initialize detention and gutter properties (defaults: no detention, no gutters)
         buildings(ibld)%detention_volume = 0.0d0
         buildings(ibld)%current_detention = 0.0d0
         buildings(ibld)%gutter_spacing = 0.0d0  ! 0 = uniform distribution (no gutters)

      end do

      close (iunit)

      write (*, '(A,I0,A)') 'Successfully read ', nbuildings, ' buildings'

   end subroutine read_building_file

   ! ========================================================================
   ! Read building properties file (.bpr)
   ! ========================================================================
   
   subroutine read_building_properties_file(filename, ios_out)

      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(out) :: ios_out

      integer :: iunit, ios_bpr, ibld, nmatched
      character(len=1024) :: line
      character(len=256) :: name
      real*8 :: det_vol, gut_spacing

      ios_out = 0
      nmatched = 0

      ! Open file
      open (newunit=iunit, file=trim(filename), status='old', action='read', iostat=ios_bpr)
      if (ios_bpr /= 0) then
         write (*, '(A,A)') 'ERROR: Could not open building properties file: ', trim(filename)
         ios_out = ios_bpr
         return
      end if

      write (*, '(A,A)') 'Reading building properties file: ', trim(filename)

      ! Read properties line by line
      do
         read (iunit, '(A)', iostat=ios_bpr) line
         if (ios_bpr /= 0) exit
         if (len_trim(line) == 0) cycle
         if (line(1:1) == '!' .or. line(1:1) == '#') cycle  ! Skip comments

         ! Parse: name, detention_volume, gutter_spacing
         read (line, *, iostat=ios_bpr) name, det_vol, gut_spacing
         if (ios_bpr /= 0) then
            write (*, '(A,A)') 'WARNING: Could not parse line: ', trim(line)
            cycle
         end if

         ! Find matching building by name
         do ibld = 1, nbuildings
            if (trim(buildings(ibld)%name) == trim(name)) then
               buildings(ibld)%detention_volume = det_vol
               buildings(ibld)%gutter_spacing = gut_spacing
               nmatched = nmatched + 1
                              
               ! Set global flags if any building uses these features
               if (det_vol > 0.0d0) use_building_detention = .true.
               if (gut_spacing > 0.0d0) use_building_gutters = .true.
               
               exit
            end if
         end do

      end do

      close (iunit)

      write (*, '(A,I0,A,I0,A)') 'Matched properties for ', nmatched, ' of ', nbuildings, ' buildings'
      
      if (use_building_detention) then
         write (*, '(A)') 'INFO: Building detention volumes enabled'
      end if
      if (use_building_gutters) then
         write (*, '(A)') 'INFO: Building gutter system enabled'
      end if

   end subroutine read_building_properties_file

   ! ========================================================================
   ! Identify which grid cells are inside buildings
   ! ========================================================================
   subroutine identify_building_cells()

      implicit none
      integer :: i, j, m, n, ibld, ip, ncells
      real*8 :: xc, yc
      logical :: inside
      integer, allocatable :: temp_cells(:)

      write (*, '(A)') 'Identifying building cells...'

      ncells = mmax*nmax

      do ibld = 1, nbuildings

         allocate (temp_cells(ncells))
         buildings(ibld)%ncells_inside = 0

         do n = 1, nmax
            do m = 1, mmax

               ! Cell center coordinates (accounting for grid rotation)
               xc = x0 + cosrot*(dble(m) - 0.5d0)*dble(dx) - sinrot*(dble(n) - 0.5d0)*dble(dy)
               yc = y0 + sinrot*(dble(m) - 0.5d0)*dble(dx) + cosrot*(dble(n) - 0.5d0)*dble(dy)

               call point_in_polygon(xc, yc, &
                                     buildings(ibld)%x, &
                                     buildings(ibld)%y, &
                                     buildings(ibld)%npoints, &
                                     inside)

               if (inside) then
                  ! Column-major indexing
                  ip = (m - 1)*nmax + n

                  ! Validate array index
                  if (ip < 1 .or. ip > np) then
                     cycle
                  end if

                  is_building_cell(ip) = .true.
                  building_id(ip) = ibld

                  ! Mark building cell as inactive (no flow computation)
                  kcs(ip) = 0

                  buildings(ibld)%ncells_inside = buildings(ibld)%ncells_inside + 1
                  temp_cells(buildings(ibld)%ncells_inside) = ip
               end if

            end do
         end do

         if (buildings(ibld)%ncells_inside > 0) then
            allocate (buildings(ibld)%cells_inside(buildings(ibld)%ncells_inside))
            buildings(ibld)%cells_inside = temp_cells(1:buildings(ibld)%ncells_inside)
         end if

         deallocate (temp_cells)

         buildings(ibld)%total_area = dble(buildings(ibld)%ncells_inside)*dx*dy
         
         ! Calculate perimeter length from polygon vertices
         call calculate_perimeter_length(ibld)

      end do

      write (*, '(A)') 'Building cell identification complete'

   end subroutine identify_building_cells

   ! ========================================================================
   ! Calculate perimeter length of a building polygon
   ! ========================================================================
   subroutine calculate_perimeter_length(ibld)

      implicit none
      integer, intent(in) :: ibld
      integer :: i
      real*8 :: length, dx_seg, dy_seg

      length = 0.0d0
      
      do i = 1, buildings(ibld)%npoints - 1
         dx_seg = buildings(ibld)%x(i+1) - buildings(ibld)%x(i)
         dy_seg = buildings(ibld)%y(i+1) - buildings(ibld)%y(i)
         length = length + sqrt(dx_seg*dx_seg + dy_seg*dy_seg)
      end do
      
      buildings(ibld)%perimeter_length = length

   end subroutine calculate_perimeter_length

   ! ========================================================================
   ! Ray-casting algorithm for point-in-polygon test
   ! ========================================================================
   subroutine point_in_polygon(xp, yp, xv, yv, nv, inside)

      implicit none
      real*8, intent(in) :: xp, yp
      real*8, intent(in) :: xv(:), yv(:)
      integer, intent(in) :: nv
      logical, intent(out) :: inside

      integer :: i, j
      logical :: intersect
      real*8 :: xi, yi, xj, yj

      inside = .false.
      j = nv

      do i = 1, nv
         xi = xv(i)
         yi = yv(i)
         xj = xv(j)
         yj = yv(j)

         intersect = ((yi > yp) .neqv. (yj > yp)) .and. &
                     (xp < (xj - xi)*(yp - yi)/(yj - yi) + xi)

         if (intersect) inside = .not. inside

         j = i
      end do

   end subroutine point_in_polygon

   ! ========================================================================
   ! Identify perimeter cells around buildings
   ! ========================================================================
   subroutine identify_perimeter_cells()

      implicit none
      integer :: ibld, ic, ip, m, n, mp, np, ipp
      integer :: dm, dn
      integer, allocatable :: temp_perimeter(:)
      integer :: ncells

      write (*, '(A)') 'Identifying perimeter cells...'

      ncells = mmax*nmax

      do ibld = 1, nbuildings

         allocate (temp_perimeter(ncells))
         buildings(ibld)%ncells_perimeter = 0

         do ic = 1, buildings(ibld)%ncells_inside
            ip = buildings(ibld)%cells_inside(ic)

            ! Reverse mapping using column-major convention
            m = (ip - 1)/nmax + 1
            n = ip - (m - 1)*nmax

            do dn = -1, 1
               do dm = -1, 1
                  if (abs(dm) + abs(dn) /= 1) cycle

                  mp = m + dm
                  np = n + dn

                  if (mp < 1 .or. mp > mmax) cycle
                  if (np < 1 .or. np > nmax) cycle

                  ! Forward mapping using column-major convention
                  ipp = (mp - 1)*nmax + np

                  if (.not. is_building_cell(ipp)) then
                     if (.not. is_perimeter_cell(ipp)) then
                        buildings(ibld)%ncells_perimeter = &
                           buildings(ibld)%ncells_perimeter + 1
                        temp_perimeter(buildings(ibld)%ncells_perimeter) = ipp
                        is_perimeter_cell(ipp) = .true.
                     end if
                  end if

               end do
            end do
         end do

         if (buildings(ibld)%ncells_perimeter > 0) then
            allocate (buildings(ibld)%cells_perimeter(buildings(ibld)%ncells_perimeter))
            buildings(ibld)%cells_perimeter = &
               temp_perimeter(1:buildings(ibld)%ncells_perimeter)
         end if

         deallocate (temp_perimeter)

      end do

      write (*, '(A)') 'Perimeter cell identification complete'

   end subroutine identify_perimeter_cells

! ========================================================================
! Identify gutter outlet cells based on gutter spacing
! ========================================================================
subroutine identify_gutter_cells()

   implicit none
   integer :: ibld, ic, ip, m, n, j, min_idx
   integer, allocatable :: temp_gutter(:), sorted_perimeter(:)
   real*8 :: spacing, cum_dist, next_gutter_dist
   real*8 :: xc, yc, xp_prev, yp_prev, seg_length
   real*8 :: x_centroid, y_centroid
   real*8, allocatable :: angles(:)
   real*8 :: temp_angle
   integer :: temp_ip
   logical :: first_cell

   write (*, '(A)') 'Identifying gutter outlet cells...'

   do ibld = 1, nbuildings

      ! Skip if no gutter spacing defined
      if (buildings(ibld)%gutter_spacing <= 0.0d0) then
         buildings(ibld)%ncells_gutter = 0
         cycle
      end if
      
      ! Skip if no perimeter cells
      if (buildings(ibld)%ncells_perimeter == 0) then
         buildings(ibld)%ncells_gutter = 0
         cycle
      end if

      spacing = buildings(ibld)%gutter_spacing

      ! ================================================================
      ! Calculate building centroid from polygon vertices
      ! (exclude last point if polygon is closed, i.e., last = first)
      ! ================================================================
      if (buildings(ibld)%npoints > 1) then
         x_centroid = sum(buildings(ibld)%x(1:buildings(ibld)%npoints-1)) / &
                      dble(buildings(ibld)%npoints - 1)
         y_centroid = sum(buildings(ibld)%y(1:buildings(ibld)%npoints-1)) / &
                      dble(buildings(ibld)%npoints - 1)
      else
         ! Fallback for single point (shouldn't happen)
         x_centroid = buildings(ibld)%x(1)
         y_centroid = buildings(ibld)%y(1)
      end if

      ! ================================================================
      ! Create sorted copy of perimeter cells by angle from centroid
      ! This ensures we walk along the perimeter in order
      ! ================================================================
      allocate(sorted_perimeter(buildings(ibld)%ncells_perimeter))
      allocate(angles(buildings(ibld)%ncells_perimeter))
      
      do ic = 1, buildings(ibld)%ncells_perimeter
         ip = buildings(ibld)%cells_perimeter(ic)
         sorted_perimeter(ic) = ip
         
         ! Get cell center coordinates
         m = (ip - 1)/nmax + 1
         n = ip - (m - 1)*nmax
         xc = x0 + cosrot*(dble(m) - 0.5d0)*dble(dx) - sinrot*(dble(n) - 0.5d0)*dble(dy)
         yc = y0 + sinrot*(dble(m) - 0.5d0)*dble(dx) + cosrot*(dble(n) - 0.5d0)*dble(dy)
         
         ! Calculate angle from centroid (radians, -pi to +pi)
         angles(ic) = atan2(yc - y_centroid, xc - x_centroid)
      end do
      
      ! ================================================================
      ! Sort perimeter cells by angle (selection sort - efficient for small arrays)
      ! ================================================================
      do ic = 1, buildings(ibld)%ncells_perimeter - 1
         min_idx = ic
         do j = ic + 1, buildings(ibld)%ncells_perimeter
            if (angles(j) < angles(min_idx)) min_idx = j
         end do
         if (min_idx /= ic) then
            ! Swap angles
            temp_angle = angles(ic)
            angles(ic) = angles(min_idx)
            angles(min_idx) = temp_angle
            
            ! Swap cell indices
            temp_ip = sorted_perimeter(ic)
            sorted_perimeter(ic) = sorted_perimeter(min_idx)
            sorted_perimeter(min_idx) = temp_ip
         end if
      end do

      ! ================================================================
      ! Walk along sorted perimeter cells and place gutters at intervals
      ! ================================================================
      allocate(temp_gutter(buildings(ibld)%ncells_perimeter))
      buildings(ibld)%ncells_gutter = 0

      cum_dist = 0.0d0
      next_gutter_dist = spacing / 2.0d0  ! First gutter at half spacing
      first_cell = .true.
      xp_prev = 0.0d0
      yp_prev = 0.0d0

      do ic = 1, buildings(ibld)%ncells_perimeter
         ip = sorted_perimeter(ic)
         
         ! Get cell center coordinates
         m = (ip - 1)/nmax + 1
         n = ip - (m - 1)*nmax
         xc = x0 + cosrot*(dble(m) - 0.5d0)*dble(dx) - sinrot*(dble(n) - 0.5d0)*dble(dy)
         yc = y0 + sinrot*(dble(m) - 0.5d0)*dble(dx) + cosrot*(dble(n) - 0.5d0)*dble(dy)

         if (.not. first_cell) then
            seg_length = sqrt((xc - xp_prev)**2 + (yc - yp_prev)**2)
            cum_dist = cum_dist + seg_length
         end if

         ! Place gutter when cumulative distance reaches threshold
         if (cum_dist >= next_gutter_dist) then
            buildings(ibld)%ncells_gutter = buildings(ibld)%ncells_gutter + 1
            temp_gutter(buildings(ibld)%ncells_gutter) = ip
            is_gutter_cell(ip) = .true.
            next_gutter_dist = next_gutter_dist + spacing
         end if

         xp_prev = xc
         yp_prev = yc
         first_cell = .false.
      end do

      ! ================================================================
      ! Ensure at least one gutter per building
      ! ================================================================
      if (buildings(ibld)%ncells_gutter == 0) then
         buildings(ibld)%ncells_gutter = 1
         temp_gutter(1) = sorted_perimeter(1)
         is_gutter_cell(temp_gutter(1)) = .true.
      end if

      ! ================================================================
      ! Store result in building structure
      ! ================================================================
      if (buildings(ibld)%ncells_gutter > 0) then
         allocate(buildings(ibld)%cells_gutter(buildings(ibld)%ncells_gutter))
         buildings(ibld)%cells_gutter = temp_gutter(1:buildings(ibld)%ncells_gutter)
      end if

      ! Clean up temporary arrays
      deallocate(temp_gutter)
      deallocate(sorted_perimeter)
      deallocate(angles)

   end do

   write (*, '(A)') 'Gutter cell identification complete'

end subroutine identify_gutter_cells


   ! ========================================================================
   ! Accumulate rainfall on building cells
   ! ========================================================================
   subroutine accumulate_building_rainfall(precip, dt)

      implicit none
      real*4, intent(inout) :: precip(:)
      real*8, intent(in) :: dt

      integer :: ibld, ic, ip
      real*8 :: rain_rate, rain_volume

      do ibld = 1, nbuildings

         do ic = 1, buildings(ibld)%ncells_inside
            ip = buildings(ibld)%cells_inside(ic)

            rain_rate = precip(ip)

            if (rain_rate > 0.0d0) then
               rain_volume = rain_rate*dt*dx*dy

               buildings(ibld)%accumulated_rainfall = &
                  buildings(ibld)%accumulated_rainfall + rain_volume

               total_rain_on_buildings = total_rain_on_buildings + rain_volume
            end if

            ! CRITICAL: Set to zero for building cells
            precip(ip) = 0.0d0

            ! Also zero netprcp if allocated (used in continuity equation)
            if (allocated(netprcp)) then
               netprcp(ip) = 0.0
            end if

         end do

      end do

   end subroutine accumulate_building_rainfall

   ! ========================================================================
   ! Redistribute accumulated water to perimeter cells or gutters
   ! ========================================================================
   subroutine redistribute_building_water(zs, dt)

      implicit none
      real*8, intent(inout) :: zs(:)
      real*8, intent(in) :: dt

      integer :: ibld, ic, ip
      real*8 :: volume_to_distribute, volume_per_cell, dh
      real*8 :: cell_area, available_detention, volume_to_detention

      cell_area = dx*dy

      do ibld = 1, nbuildings

         ! Skip if no accumulated rainfall
         if (buildings(ibld)%accumulated_rainfall <= 1.0d-100) cycle

         volume_to_distribute = buildings(ibld)%accumulated_rainfall

         ! ================================================================
         ! Handle detention volume (if enabled for this building)
         ! ================================================================
         if (use_building_detention .and. buildings(ibld)%detention_volume > 0.0d0) then
            
            ! Calculate available detention capacity
            available_detention = buildings(ibld)%detention_volume - buildings(ibld)%current_detention
            
            if (available_detention > 0.0d0) then
               ! Store water in detention (up to available capacity)
               volume_to_detention = min(volume_to_distribute, available_detention)
               buildings(ibld)%current_detention = buildings(ibld)%current_detention + volume_to_detention
               volume_to_distribute = volume_to_distribute - volume_to_detention
               
               ! Update global detention tracking
               total_detention_stored = total_detention_stored + volume_to_detention
            end if
            
         end if

         ! Skip redistribution if all water went to detention
         if (volume_to_distribute <= 1.0d-12) then
            buildings(ibld)%accumulated_rainfall = 0.0d0
            cycle
         end if

         ! ================================================================
         ! Redistribute remaining water
         ! ================================================================
         if (use_building_gutters .and. buildings(ibld)%ncells_gutter > 0) then
            ! Use gutter outlets
            if (buildings(ibld)%ncells_gutter == 0) cycle

            volume_per_cell = volume_to_distribute / dble(buildings(ibld)%ncells_gutter)
            dh = volume_per_cell / cell_area

            do ic = 1, buildings(ibld)%ncells_gutter
               ip = buildings(ibld)%cells_gutter(ic)
               zs(ip) = zs(ip) + real(dh, 4)
            end do

         else
            ! Use uniform distribution to all perimeter cells 
            if (buildings(ibld)%ncells_perimeter == 0) cycle

            volume_per_cell = volume_to_distribute / dble(buildings(ibld)%ncells_perimeter)
            dh = volume_per_cell / cell_area

            do ic = 1, buildings(ibld)%ncells_perimeter
               ip = buildings(ibld)%cells_perimeter(ic)
               zs(ip) = zs(ip) + real(dh, 4)
            end do

         end if

         ! Update statistics
         total_runoff_from_buildings = total_runoff_from_buildings + volume_to_distribute
         buildings(ibld)%accumulated_runoff = buildings(ibld)%accumulated_runoff + volume_to_distribute

         ! Reset accumulated rainfall
         buildings(ibld)%accumulated_rainfall = 0.0d0

      end do

   end subroutine redistribute_building_water

   ! ========================================================================
   ! Get building edge coordinates for thin dam creation
   ! ========================================================================
   subroutine get_building_edges(ibld, npoints, x, y)

      implicit none
      integer, intent(in) :: ibld
      integer, intent(out) :: npoints
      real*8, allocatable, intent(out) :: x(:), y(:)

      if (ibld < 1 .or. ibld > nbuildings) then
         write (*, *) 'ERROR: Invalid building index in get_building_edges'
         npoints = 0
         return
      end if

      npoints = buildings(ibld)%npoints
      allocate (x(npoints))
      allocate (y(npoints))

      x = buildings(ibld)%x
      y = buildings(ibld)%y

   end subroutine get_building_edges

end module sfincs_buildings
