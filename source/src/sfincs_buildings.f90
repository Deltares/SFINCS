module sfincs_buildings
   !
   ! Building drainage module for SFINCS
   ! 1. Collect rainfall on building footprint
   ! 2. Block flow with thin dams on edges
   ! 3. Redistribute collected water to perimeter cells OR gutter outlets
   ! 4. Optional detention volume per building
   ! 5. Optional gutter system with configurable spacing
   ! IMDC, 2025 (restructured by Deltares)
   !
   use sfincs_log
   !
   implicit none
   !
   type building_type
      character(len=256) :: name
      integer :: npoints
      real*8, allocatable :: x(:)
      real*8, allocatable :: y(:)
      integer, allocatable :: cells_inside(:)
      integer :: ncells_inside
      integer, allocatable :: cells_perimeter(:)
      integer :: ncells_perimeter
      integer, allocatable :: cells_gutter(:)
      integer :: ncells_gutter
      real*8 :: total_area
      real*8 :: accumulated_rainfall
      real*8 :: accumulated_runoff
      real*8 :: detention_volume
      real*8 :: current_detention
      real*8 :: gutter_spacing
      real*8 :: perimeter_length
   end type building_type
   !
   integer :: nbuildings = 0
   type(building_type), allocatable :: buildings(:)
   !
   logical, allocatable :: is_building_cell(:)
   logical, allocatable :: is_perimeter_cell(:)
   logical, allocatable :: is_gutter_cell(:)
   integer, allocatable :: building_id(:)
   !
   real*8 :: total_rain_on_buildings = 0.0d0
   real*8 :: total_runoff_from_buildings = 0.0d0
   real*8 :: total_detention_stored = 0.0d0
   !
   logical :: use_building_detention = .false.
   logical :: use_building_gutters = .false.
   !
contains
   !
   subroutine setup_buildings()
   !
   ! Read files, identify cells, create thin dams.
   ! Called once from sfincs_domain after mesh initialization.
   !
   use sfincs_data, only: has_buildings, bldfile, bprfile
   !
   implicit none
   !
   integer :: ios
   !
   call write_log('Info    : Setting up buildings ...', 1)
   !
   call initialize_building_masks()
   !
   call read_building_file(bldfile, ios)
   !
   if (ios /= 0) then
      call write_log('ERROR: Failed to read building file', 0)
      has_buildings = .false.
      return
   end if
   !
   if (len_trim(bprfile) > 0 .and. trim(bprfile) /= 'none') then
      call read_building_properties_file(bprfile, ios)
      if (ios /= 0) then
         call write_log('WARNING: Failed to read building properties file', 0)
         call write_log('         Using default values (no detention, uniform distribution)', 0)
      end if
   end if
   !
   call identify_building_cells()
   call identify_perimeter_cells()
   !
   if (use_building_gutters) then
      call identify_gutter_cells()
   end if
   !
   call create_building_thin_dams()
   !
   call write_log('Info    : Building setup complete', 0)
   !
   end subroutine setup_buildings
   !
   !
   subroutine update_buildings(precip, zs, dt)
   !
   ! Accumulate rainfall and redistribute. Called each timestep.
   !
   implicit none
   !
   real*4, intent(inout) :: precip(:)
   real*8, intent(inout) :: zs(:)
   real*8, intent(in) :: dt
   !
   if (nbuildings == 0) return
   !
   call accumulate_building_rainfall(precip, dt)
   call redistribute_building_water(zs, dt)
   !
   end subroutine update_buildings
   !
   !
   subroutine initialize_building_masks()
   !
   use sfincs_data, only: mmax, nmax
   !
   implicit none
   !
   integer :: ncells
   !
   ncells = mmax * nmax
   !
   if (.not. allocated(is_building_cell)) then
      allocate(is_building_cell(ncells))
      allocate(is_perimeter_cell(ncells))
      allocate(is_gutter_cell(ncells))
      allocate(building_id(ncells))
   end if
   !
   is_building_cell = .false.
   is_perimeter_cell = .false.
   is_gutter_cell = .false.
   building_id = 0
   !
   end subroutine initialize_building_masks
   !
   !
   subroutine read_building_file(filename, ios_out)
   !
   implicit none
   !
   character(len=*), intent(in) :: filename
   integer, intent(out) :: ios_out
   !
   integer :: iunit, ios_bld, i, ibld
   character(len=1024) :: line
   character(len=256) :: name
   integer :: npoints, ncols
   real*8 :: xtmp, ytmp
   !
   ios_out = 0
   !
   open(newunit=iunit, file=trim(filename), status='old', action='read', iostat=ios_bld)
   if (ios_bld /= 0) then
      write(logstr,'(A,A)') 'ERROR: Could not open building file: ', trim(filename)
      call write_log(logstr, 0)
      ios_out = ios_bld
      return
   end if
   !
   write(logstr,'(A,A)') 'Info    : Reading building file: ', trim(filename)
   call write_log(logstr, 0)
   !
   ! Count buildings
   !
   nbuildings = 0
   do
      read(iunit, '(A)', iostat=ios_bld) line
      if (ios_bld /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '!') cycle
      read(line, *, iostat=ios_bld) name, npoints, ncols
      if (ios_bld == 0 .and. npoints > 0 .and. ncols >= 2) then
         nbuildings = nbuildings + 1
         do i = 1, npoints
            read(iunit, *, iostat=ios_bld)
            if (ios_bld /= 0) exit
         end do
      end if
   end do
   !
   if (nbuildings == 0) then
      call write_log('WARNING: No buildings found in file', 0)
      close(iunit)
      return
   end if
   !
   write(logstr,'(A,I0,A)') 'Info    : Found ', nbuildings, ' buildings'
   call write_log(logstr, 0)
   !
   allocate(buildings(nbuildings))
   !
   rewind(iunit)
   !
   ibld = 0
   do
      read(iunit, '(A)', iostat=ios_bld) line
      if (ios_bld /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '!') cycle
      read(line, *, iostat=ios_bld) name, npoints, ncols
      if (ios_bld /= 0 .or. npoints <= 0 .or. ncols < 2) cycle
      !
      ibld = ibld + 1
      if (ibld > nbuildings) exit
      !
      buildings(ibld)%name = trim(name)
      buildings(ibld)%npoints = npoints
      allocate(buildings(ibld)%x(npoints))
      allocate(buildings(ibld)%y(npoints))
      !
      do i = 1, npoints
         read(iunit, *, iostat=ios_bld) xtmp, ytmp
         if (ios_bld /= 0) then
            write(logstr,'(A,A)') 'ERROR: Reading coordinates for building ', trim(name)
            call write_log(logstr, 0)
            ios_out = ios_bld
            close(iunit)
            return
         end if
         buildings(ibld)%x(i) = xtmp
         buildings(ibld)%y(i) = ytmp
      end do
      !
      if (buildings(ibld)%x(1) /= buildings(ibld)%x(npoints) .or. &
          buildings(ibld)%y(1) /= buildings(ibld)%y(npoints)) then
         write(logstr,'(A,A,A)') 'WARNING: Building ', trim(name), ' polygon not closed'
         call write_log(logstr, 0)
      end if
      !
      buildings(ibld)%accumulated_rainfall = 0.0d0
      buildings(ibld)%accumulated_runoff = 0.0d0
      buildings(ibld)%ncells_inside = 0
      buildings(ibld)%ncells_perimeter = 0
      buildings(ibld)%ncells_gutter = 0
      buildings(ibld)%total_area = 0.0d0
      buildings(ibld)%perimeter_length = 0.0d0
      buildings(ibld)%detention_volume = 0.0d0
      buildings(ibld)%current_detention = 0.0d0
      buildings(ibld)%gutter_spacing = 0.0d0
      !
   end do
   !
   close(iunit)
   !
   write(logstr,'(A,I0,A)') 'Info    : Successfully read ', nbuildings, ' buildings'
   call write_log(logstr, 0)
   !
   end subroutine read_building_file
   !
   !
   subroutine read_building_properties_file(filename, ios_out)
   !
   implicit none
   !
   character(len=*), intent(in) :: filename
   integer, intent(out) :: ios_out
   !
   integer :: iunit, ios_bpr, ibld, nmatched
   character(len=1024) :: line
   character(len=256) :: name
   real*8 :: det_vol, gut_spacing
   !
   ios_out = 0
   nmatched = 0
   !
   open(newunit=iunit, file=trim(filename), status='old', action='read', iostat=ios_bpr)
   if (ios_bpr /= 0) then
      write(logstr,'(A,A)') 'ERROR: Could not open building properties file: ', trim(filename)
      call write_log(logstr, 0)
      ios_out = ios_bpr
      return
   end if
   !
   write(logstr,'(A,A)') 'Info    : Reading building properties file: ', trim(filename)
   call write_log(logstr, 0)
   !
   do
      read(iunit, '(A)', iostat=ios_bpr) line
      if (ios_bpr /= 0) exit
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '!' .or. line(1:1) == '#') cycle
      !
      read(line, *, iostat=ios_bpr) name, det_vol, gut_spacing
      if (ios_bpr /= 0) then
         write(logstr,'(A,A)') 'WARNING: Could not parse line: ', trim(line)
         call write_log(logstr, 0)
         cycle
      end if
      !
      do ibld = 1, nbuildings
         if (trim(buildings(ibld)%name) == trim(name)) then
            buildings(ibld)%detention_volume = det_vol
            buildings(ibld)%gutter_spacing = gut_spacing
            nmatched = nmatched + 1
            if (det_vol > 0.0d0) use_building_detention = .true.
            if (gut_spacing > 0.0d0) use_building_gutters = .true.
            exit
         end if
      end do
      !
   end do
   !
   close(iunit)
   !
   write(logstr,'(A,I0,A,I0,A)') 'Info    : Matched properties for ', nmatched, ' of ', nbuildings, ' buildings'
   call write_log(logstr, 0)
   if (use_building_detention) call write_log('Info    : Building detention volumes enabled', 0)
   if (use_building_gutters) call write_log('Info    : Building gutter system enabled', 0)
   !
   end subroutine read_building_properties_file
   !
   !
   subroutine identify_building_cells()
   !
   use sfincs_data, only: mmax, nmax, np, x0, y0, dx, dy, cosrot, sinrot, kcs
   !
   implicit none
   !
   integer :: m, n, ibld, ip, ncells
   real*8 :: xc, yc
   logical :: inside
   integer, allocatable :: temp_cells(:)
   !
   call write_log('Info    : Identifying building cells ...', 0)
   !
   ncells = mmax * nmax
   !
   do ibld = 1, nbuildings
      !
      allocate(temp_cells(ncells))
      buildings(ibld)%ncells_inside = 0
      !
      do n = 1, nmax
         do m = 1, mmax
            !
            xc = x0 + cosrot*(dble(m) - 0.5d0)*dble(dx) - sinrot*(dble(n) - 0.5d0)*dble(dy)
            yc = y0 + sinrot*(dble(m) - 0.5d0)*dble(dx) + cosrot*(dble(n) - 0.5d0)*dble(dy)
            !
            call point_in_polygon(xc, yc, buildings(ibld)%x, buildings(ibld)%y, &
                                  buildings(ibld)%npoints, inside)
            !
            if (inside) then
               ip = (m - 1)*nmax + n
               if (ip < 1 .or. ip > np) cycle
               !
               is_building_cell(ip) = .true.
               building_id(ip) = ibld
               !
               ! Set kcs=0 to make cell inactive: suppresses continuity update
               ! AND prevents any boundary or inflow assignment on building cells.
               !
               kcs(ip) = 0
               !
               buildings(ibld)%ncells_inside = buildings(ibld)%ncells_inside + 1
               temp_cells(buildings(ibld)%ncells_inside) = ip
            end if
            !
         end do
      end do
      !
      if (buildings(ibld)%ncells_inside > 0) then
         allocate(buildings(ibld)%cells_inside(buildings(ibld)%ncells_inside))
         buildings(ibld)%cells_inside = temp_cells(1:buildings(ibld)%ncells_inside)
      end if
      !
      deallocate(temp_cells)
      !
      buildings(ibld)%total_area = dble(buildings(ibld)%ncells_inside) * dx * dy
      call calculate_perimeter_length(ibld)
      !
   end do
   !
   call write_log('Info    : Building cell identification complete', 0)
   !
   end subroutine identify_building_cells
   !
   !
   subroutine calculate_perimeter_length(ibld)
   !
   implicit none
   !
   integer, intent(in) :: ibld
   integer :: i
   real*8 :: length, dx_seg, dy_seg
   !
   length = 0.0d0
   do i = 1, buildings(ibld)%npoints - 1
      dx_seg = buildings(ibld)%x(i+1) - buildings(ibld)%x(i)
      dy_seg = buildings(ibld)%y(i+1) - buildings(ibld)%y(i)
      length = length + sqrt(dx_seg*dx_seg + dy_seg*dy_seg)
   end do
   buildings(ibld)%perimeter_length = length
   !
   end subroutine calculate_perimeter_length
   !
   !
   subroutine point_in_polygon(xp, yp, xv, yv, nv, inside)
   !
   implicit none
   !
   real*8, intent(in) :: xp, yp
   real*8, intent(in) :: xv(:), yv(:)
   integer, intent(in) :: nv
   logical, intent(out) :: inside
   !
   integer :: i, j
   logical :: intersect
   !
   inside = .false.
   j = nv
   do i = 1, nv
      ! Guard against horizontal edges (yv(j)==yv(i)) before dividing:
      !
      if ((yv(i) > yp) .neqv. (yv(j) > yp)) then
         if (xp < (xv(j) - xv(i))*(yp - yv(i))/(yv(j) - yv(i)) + xv(i)) &
            inside = .not. inside
      end if
      j = i
   end do
   !
   end subroutine point_in_polygon
   !
   !
   subroutine identify_perimeter_cells()
   !
   use sfincs_data, only: mmax, nmax
   !
   implicit none
   !
   integer :: ibld, ic, ip, m, n, mp, np_local, ipp
   integer :: dm, dn
   integer, allocatable :: temp_perimeter(:)
   integer :: ncells
   !
   call write_log('Info    : Identifying perimeter cells ...', 0)
   !
   ncells = mmax * nmax
   !
   do ibld = 1, nbuildings
      !
      allocate(temp_perimeter(ncells))
      buildings(ibld)%ncells_perimeter = 0
      !
      do ic = 1, buildings(ibld)%ncells_inside
         ip = buildings(ibld)%cells_inside(ic)
         m = (ip - 1)/nmax + 1
         n = ip - (m - 1)*nmax
         !
         do dn = -1, 1
            do dm = -1, 1
               if (abs(dm) + abs(dn) /= 1) cycle
               mp = m + dm
               np_local = n + dn
               if (mp < 1 .or. mp > mmax) cycle
               if (np_local < 1 .or. np_local > nmax) cycle
               !
               ipp = (mp - 1)*nmax + np_local
               !
               if (.not. is_building_cell(ipp) .and. .not. is_perimeter_cell(ipp)) then
                  buildings(ibld)%ncells_perimeter = buildings(ibld)%ncells_perimeter + 1
                  temp_perimeter(buildings(ibld)%ncells_perimeter) = ipp
                  is_perimeter_cell(ipp) = .true.
               end if
            end do
         end do
      end do
      !
      if (buildings(ibld)%ncells_perimeter > 0) then
         allocate(buildings(ibld)%cells_perimeter(buildings(ibld)%ncells_perimeter))
         buildings(ibld)%cells_perimeter = temp_perimeter(1:buildings(ibld)%ncells_perimeter)
      end if
      !
      deallocate(temp_perimeter)
      !
   end do
   !
   call write_log('Info    : Perimeter cell identification complete', 0)
   !
   end subroutine identify_perimeter_cells
   !
   !
   subroutine identify_gutter_cells()
   !
   use sfincs_data, only: mmax, nmax, x0, y0, dx, dy, cosrot, sinrot
   !
   implicit none
   !
   integer :: ibld, ic, ip, m, n, j, min_idx
   integer, allocatable :: temp_gutter(:), sorted_perimeter(:)
   real*8 :: spacing, cum_dist, next_gutter_dist
   real*8 :: xc, yc, xp_prev, yp_prev
   real*8 :: x_centroid, y_centroid
   real*8, allocatable :: angles(:)
   real*8 :: temp_angle
   integer :: temp_ip
   logical :: first_cell
   !
   call write_log('Info    : Identifying gutter outlet cells ...', 0)
   !
   do ibld = 1, nbuildings
      !
      if (buildings(ibld)%gutter_spacing <= 0.0d0) then
         buildings(ibld)%ncells_gutter = 0
         cycle
      end if
      !
      if (buildings(ibld)%ncells_perimeter == 0) then
         buildings(ibld)%ncells_gutter = 0
         cycle
      end if
      !
      spacing = buildings(ibld)%gutter_spacing
      !
      ! Building centroid from polygon vertices (exclude closing point)
      !
      if (buildings(ibld)%npoints > 1) then
         x_centroid = sum(buildings(ibld)%x(1:buildings(ibld)%npoints-1)) / &
                      dble(buildings(ibld)%npoints - 1)
         y_centroid = sum(buildings(ibld)%y(1:buildings(ibld)%npoints-1)) / &
                      dble(buildings(ibld)%npoints - 1)
      else
         x_centroid = buildings(ibld)%x(1)
         y_centroid = buildings(ibld)%y(1)
      end if
      !
      ! Sort perimeter cells by angle from centroid
      !
      allocate(sorted_perimeter(buildings(ibld)%ncells_perimeter))
      allocate(angles(buildings(ibld)%ncells_perimeter))
      !
      do ic = 1, buildings(ibld)%ncells_perimeter
         ip = buildings(ibld)%cells_perimeter(ic)
         sorted_perimeter(ic) = ip
         m = (ip - 1)/nmax + 1
         n = ip - (m - 1)*nmax
         xc = x0 + cosrot*(dble(m) - 0.5d0)*dble(dx) - sinrot*(dble(n) - 0.5d0)*dble(dy)
         yc = y0 + sinrot*(dble(m) - 0.5d0)*dble(dx) + cosrot*(dble(n) - 0.5d0)*dble(dy)
         angles(ic) = atan2(yc - y_centroid, xc - x_centroid)
      end do
      !
      ! Selection sort by angle
      !
      do ic = 1, buildings(ibld)%ncells_perimeter - 1
         min_idx = ic
         do j = ic + 1, buildings(ibld)%ncells_perimeter
            if (angles(j) < angles(min_idx)) min_idx = j
         end do
         if (min_idx /= ic) then
            temp_angle = angles(ic)
            angles(ic) = angles(min_idx)
            angles(min_idx) = temp_angle
            temp_ip = sorted_perimeter(ic)
            sorted_perimeter(ic) = sorted_perimeter(min_idx)
            sorted_perimeter(min_idx) = temp_ip
         end if
      end do
      !
      ! Walk sorted perimeter and place gutters at intervals
      !
      allocate(temp_gutter(buildings(ibld)%ncells_perimeter))
      buildings(ibld)%ncells_gutter = 0
      cum_dist = 0.0d0
      next_gutter_dist = spacing / 2.0d0
      first_cell = .true.
      !
      do ic = 1, buildings(ibld)%ncells_perimeter
         ip = sorted_perimeter(ic)
         m = (ip - 1)/nmax + 1
         n = ip - (m - 1)*nmax
         xc = x0 + cosrot*(dble(m) - 0.5d0)*dble(dx) - sinrot*(dble(n) - 0.5d0)*dble(dy)
         yc = y0 + sinrot*(dble(m) - 0.5d0)*dble(dx) + cosrot*(dble(n) - 0.5d0)*dble(dy)
         !
         if (.not. first_cell) then
            cum_dist = cum_dist + sqrt((xc - xp_prev)**2 + (yc - yp_prev)**2)
         end if
         !
         if (cum_dist >= next_gutter_dist) then
            buildings(ibld)%ncells_gutter = buildings(ibld)%ncells_gutter + 1
            temp_gutter(buildings(ibld)%ncells_gutter) = ip
            is_gutter_cell(ip) = .true.
            next_gutter_dist = next_gutter_dist + spacing
         end if
         !
         xp_prev = xc
         yp_prev = yc
         first_cell = .false.
      end do
      !
      ! Ensure at least one gutter per building
      !
      if (buildings(ibld)%ncells_gutter == 0) then
         buildings(ibld)%ncells_gutter = 1
         temp_gutter(1) = sorted_perimeter(1)
         is_gutter_cell(temp_gutter(1)) = .true.
      end if
      !
      if (buildings(ibld)%ncells_gutter > 0) then
         allocate(buildings(ibld)%cells_gutter(buildings(ibld)%ncells_gutter))
         buildings(ibld)%cells_gutter = temp_gutter(1:buildings(ibld)%ncells_gutter)
      end if
      !
      deallocate(temp_gutter, sorted_perimeter, angles)
      !
   end do
   !
   call write_log('Info    : Gutter cell identification complete', 0)
   !
   end subroutine identify_gutter_cells
   !
   !
   subroutine create_building_thin_dams()
   !
   use sfincs_data, only: kcuv
   use quadtree, only: find_uv_points_intersected_by_polyline
   !
   implicit none
   !
   integer :: ibld, iuv, nrows, nr_points, indx
   real*4, allocatable :: xthd(:), ythd(:)
   integer, allocatable :: uv_indices(:), vertices(:)
   !
   write(logstr,'(A,I0,A)') 'Info    : Processing ', nbuildings, ' buildings as thin dams'
   call write_log(logstr, 0)
   !
   do ibld = 1, nbuildings
      !
      nrows = buildings(ibld)%npoints
      allocate(xthd(nrows), ythd(nrows))
      !
      xthd = real(buildings(ibld)%x, 4)
      ythd = real(buildings(ibld)%y, 4)
      !
      call find_uv_points_intersected_by_polyline(uv_indices, vertices, nr_points, xthd, ythd, nrows)
      !
      do iuv = 1, nr_points
         indx = uv_indices(iuv)
         kcuv(indx) = 0
      end do
      !
      deallocate(xthd, ythd)
      if (allocated(uv_indices)) deallocate(uv_indices)
      if (allocated(vertices))   deallocate(vertices)
      !
   end do
   !
   write(logstr,'(A,I0,A)') 'Info    : Processed ', nbuildings, ' buildings as thin dams'
   call write_log(logstr, 0)
   !
   end subroutine create_building_thin_dams
   !
   !
   subroutine accumulate_building_rainfall(precip, dt)
   !
   ! precip is rainfall rate in m/s; dt in seconds; dx/dy in m.
   ! rain_volume = rate [m/s] * dt [s] * area [m^2] => volume in m^3.
   !
   use sfincs_data, only: dx, dy, netprcp
   !
   implicit none
   !
   real*4, intent(inout) :: precip(:)
   real*8, intent(in) :: dt
   !
   integer :: ibld, ic, ip
   real*8 :: rain_volume
   !
   do ibld = 1, nbuildings
      do ic = 1, buildings(ibld)%ncells_inside
         ip = buildings(ibld)%cells_inside(ic)
         !
         if (precip(ip) > 0.0) then
            rain_volume = dble(precip(ip)) * dt * dx * dy
            buildings(ibld)%accumulated_rainfall = buildings(ibld)%accumulated_rainfall + rain_volume
            total_rain_on_buildings = total_rain_on_buildings + rain_volume
         end if
         !
         precip(ip) = 0.0
         if (allocated(netprcp)) netprcp(ip) = 0.0
         !
      end do
   end do
   !
   end subroutine accumulate_building_rainfall
   !
   !
   subroutine redistribute_building_water(zs, dt)
   !
   use sfincs_data, only: dx, dy
   !
   implicit none
   !
   real*8, intent(inout) :: zs(:)
   real*8, intent(in) :: dt
   !
   integer :: ibld, ic, ip
   real*8 :: volume_to_distribute, volume_per_cell, dh
   real*8 :: cell_area, available_detention, volume_to_detention
   !
   cell_area = dx * dy
   !
   do ibld = 1, nbuildings
      !
      if (buildings(ibld)%accumulated_rainfall <= 1.0d-100) cycle
      !
      volume_to_distribute = buildings(ibld)%accumulated_rainfall
      !
      ! Detention storage
      !
      if (use_building_detention .and. buildings(ibld)%detention_volume > 0.0d0) then
         available_detention = buildings(ibld)%detention_volume - buildings(ibld)%current_detention
         if (available_detention > 0.0d0) then
            volume_to_detention = min(volume_to_distribute, available_detention)
            buildings(ibld)%current_detention = buildings(ibld)%current_detention + volume_to_detention
            volume_to_distribute = volume_to_distribute - volume_to_detention
            total_detention_stored = total_detention_stored + volume_to_detention
         end if
      end if
      !
      if (volume_to_distribute <= 1.0d-12) then
         buildings(ibld)%accumulated_rainfall = 0.0d0
         cycle
      end if
      !
      ! Redistribute to gutters or perimeter
      !
      if (use_building_gutters .and. buildings(ibld)%ncells_gutter > 0) then
         volume_per_cell = volume_to_distribute / dble(buildings(ibld)%ncells_gutter)
         dh = volume_per_cell / cell_area
         do ic = 1, buildings(ibld)%ncells_gutter
            ip = buildings(ibld)%cells_gutter(ic)
            zs(ip) = zs(ip) + dh
         end do
      else
         if (buildings(ibld)%ncells_perimeter == 0) cycle
         volume_per_cell = volume_to_distribute / dble(buildings(ibld)%ncells_perimeter)
         dh = volume_per_cell / cell_area
         do ic = 1, buildings(ibld)%ncells_perimeter
            ip = buildings(ibld)%cells_perimeter(ic)
            zs(ip) = zs(ip) + dh
         end do
      end if
      !
      total_runoff_from_buildings = total_runoff_from_buildings + volume_to_distribute
      buildings(ibld)%accumulated_runoff = buildings(ibld)%accumulated_runoff + volume_to_distribute
      buildings(ibld)%accumulated_rainfall = 0.0d0
      !
   end do
   !
   end subroutine redistribute_building_water
   !
   !
   subroutine finalize_buildings()
   !
   ! Deallocate all module-level and per-building allocatables.
   ! Called from sfincs_finalize to avoid leaks in DLL/BMI use.
   !
   implicit none
   !
   integer :: ibld
   !
   if (allocated(buildings)) then
      do ibld = 1, nbuildings
         if (allocated(buildings(ibld)%x))              deallocate(buildings(ibld)%x)
         if (allocated(buildings(ibld)%y))              deallocate(buildings(ibld)%y)
         if (allocated(buildings(ibld)%cells_inside))   deallocate(buildings(ibld)%cells_inside)
         if (allocated(buildings(ibld)%cells_perimeter))deallocate(buildings(ibld)%cells_perimeter)
         if (allocated(buildings(ibld)%cells_gutter))   deallocate(buildings(ibld)%cells_gutter)
      end do
      deallocate(buildings)
   end if
   !
   if (allocated(is_building_cell))  deallocate(is_building_cell)
   if (allocated(is_perimeter_cell)) deallocate(is_perimeter_cell)
   if (allocated(is_gutter_cell))    deallocate(is_gutter_cell)
   if (allocated(building_id))       deallocate(building_id)
   !
   nbuildings = 0
   use_building_detention = .false.
   use_building_gutters   = .false.
   total_rain_on_buildings    = 0.0d0
   total_runoff_from_buildings = 0.0d0
   total_detention_stored     = 0.0d0
   !
   end subroutine finalize_buildings
   !
end module sfincs_buildings
