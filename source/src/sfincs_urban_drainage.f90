module sfincs_urban_drainage
   !
   ! Simple urban-drainage sink/source model for SFINCS.
   !
   ! Each zone is a polygon in horizontal plane. Cells inside the polygon
   ! drain at a design rate capped by available water; flow is bidirectional
   ! so the outfall cell can push water back into the cells (tide / surge)
   ! unless a check valve is specified. All flow for a zone is collected at a
   ! single outfall cell, so the per-zone net flux is added as a point
   ! source / sink there.
   !
   ! Per-cell discharge (drain from cell to outfall, positive sign):
   !
   !    dzs = zs(nm) - zs(outfall)
   !    if dzs > 0:
   !       q = min( qmax(nm), max(zs(nm)-zb(nm),0) * cell_area(nm) / dt )
   !       gated further by h_threshold on cell water depth
   !    else:
   !       q = -backflow_coef(nm) * sqrt(-dzs), capped at -qmax(nm)
   !       suppressed if the zone has a check valve
   !
   ! Per-cell design-head:
   !
   !    dh_design(nm) = max( zb(nm) - zb(outfall), dh_design_min )
   !    backflow_coef(nm) = qmax(nm) / sqrt(dh_design(nm))
   !
   ! qmax from the design precipitation rate:
   !
   !    qmax(nm) = design_precip_mm_hr * 1e-3 / 3600 * cell_area(nm)  [m3/s]
   !
   ! Alternatively the user may supply max_outfall_rate [m3/s] for the
   ! whole zone (exclusive with design_precip per zone); design_precip is
   ! then derived as
   !
   !    design_precip_mm_hr = max_outfall_rate / zone_area * 1000 * 3600
   !
   ! which distributes the capacity proportionally to cell area.
   !
   ! Subroutines:
   !
   !   initialize_urban_drainage()
   !     Top-level driver. Calls read_urban_drainage, loads polygons, marks
   !     cells per zone (last zone wins on overlap), snaps outfall coords to
   !     the nearest active cell, precomputes per-cell qmax and
   !     backflow_coef. Called from sfincs_lib (once at init time).
   !
   !   read_urban_drainage(filename, ierr)
   !     Parses the *.urb TOML file into the per-zone arrays. Called from
   !     initialize_urban_drainage (this module).
   !
   !   update_urban_drainage(t, dt)
   !     Per-time-step entry: accumulates signed discharges into qsrc, and
   !     adds the outfall contribution at each zone's outfall cell. Called
   !     from update_continuity (sfincs_continuity).
   !
   !   write_urban_drainage_log_summary()
   !     Prints a one-block-per-zone summary (name, polygon file, cell
   !     count, total area, design precip, total qmax, thresholds, outfall)
   !     to the log. Called from initialize_urban_drainage (this module).
   !
   use sfincs_log
   use sfincs_error
   use sfincs_polygons
   !
   implicit none
   !
   private
   !
   public :: initialize_urban_drainage
   public :: update_urban_drainage
   !
   ! Per-zone runtime state. Sized nr_urban_drainage_zones.
   !
   integer, public :: nr_urban_drainage_zones = 0
   !
   character(len=64), dimension(:), allocatable, public :: urb_zone_name
   character(len=64), dimension(:), allocatable, public :: urb_zone_type
   character(len=256), dimension(:), allocatable        :: urb_zone_polygon_file
   !
   real*4,  dimension(:), allocatable, public :: urb_zone_outfall_x
   real*4,  dimension(:), allocatable, public :: urb_zone_outfall_y
   real*4,  dimension(:), allocatable, public :: urb_zone_design_precip      ! mm/hr (either given directly or derived from max_outfall_rate)
   real*4,  dimension(:), allocatable, public :: urb_zone_max_outfall_rate   ! m3/s; 0.0 if input was design_precip instead
   real*4,  dimension(:), allocatable, public :: urb_zone_h_threshold        ! m ponding threshold
   real*4,  dimension(:), allocatable, public :: urb_zone_dh_design_min      ! m floor on design head
   logical, dimension(:), allocatable, public :: urb_zone_include_outfall
   logical, dimension(:), allocatable, public :: urb_zone_check_valve
   !
   integer, dimension(:), allocatable, public :: urban_drainage_outfall_index  ! cell index, 0 if none
   real*4,  dimension(:), allocatable, public :: urban_drainage_q_outfall     ! m3/s per zone, per step
   real*4,  dimension(:), allocatable, public :: urb_zone_area                ! m2, sum of cell areas in zone
   integer, dimension(:), allocatable, public :: urb_zone_n_cells             ! number of cells in zone
   real*4,  dimension(:), allocatable, public :: urb_zone_qmax_total          ! m3/s, sum of per-cell qmax
   !
   ! Per-cell runtime state. Sized np.
   !
   integer, dimension(:), allocatable, public :: urban_drainage_zone_indices  ! 0 if not in any zone
   real*4,  dimension(:), allocatable, public :: urban_drainage_qmax          ! m3/s cap per cell
   real*4,  dimension(:), allocatable, public :: urban_drainage_backflow_coef ! qmax / sqrt(dh_design)
   real*4,  dimension(:), allocatable, public :: urban_drainage_cumulative_volume       ! m3 accumulated (optional)
   !
contains
   !
   subroutine initialize_urban_drainage()
      !
      ! Top-level initializer for urban drainage. Parses *.urb TOML file,
      ! loads polygons, stamps cells per zone (last-wins on overlap), snaps
      ! outfall coords to the nearest active cell, and precomputes per-cell
      ! qmax and backflow coefficients.
      !
      ! Sets sfincs_data::urban_drainage = .true. when at least one zone is
      ! loaded and has at least one participating cell. Otherwise leaves it
      ! .false. and returns early.
      !
      ! Called from: sfincs_lib (once, at init time, after
      ! initialize_src_structures).
      !
      use sfincs_data
      use quadtree
      !
      implicit none
      !
      integer                       :: ierr, ipoly, iz, nm, io
      integer                       :: n_cells_in_zones, n_outfalls, nmq
      real*4                        :: area_nm, dzb, dh_min
      type(t_polygon), allocatable  :: polygons(:)
      logical, allocatable          :: inside(:)
      character(len=256)            :: last_file
      integer                       :: ip
      !
      urban_drainage = .false.
      !
      if (urbfile(1:4) == 'none') return
      !
      call write_log('Info    : reading urban drainage file ...', 0)
      !
      call read_urban_drainage(trim(urbfile), ierr)
      !
      if (ierr /= 0) then
         call stop_sfincs('Error ! Failed to read urban drainage TOML file.', -1)
         return
      endif
      !
      if (nr_urban_drainage_zones <= 0) then
         call write_log('Info    : urban drainage file contains no zones; feature disabled', 0)
         return
      endif
      !
      ! Allocate per-zone snapped outfall index and per-step accumulator.
      !
      allocate(urban_drainage_outfall_index(nr_urban_drainage_zones))
      allocate(urban_drainage_q_outfall(nr_urban_drainage_zones))
      allocate(urb_zone_area(nr_urban_drainage_zones))
      allocate(urb_zone_n_cells(nr_urban_drainage_zones))
      allocate(urb_zone_qmax_total(nr_urban_drainage_zones))
      urban_drainage_outfall_index = 0
      urban_drainage_q_outfall     = 0.0
      urb_zone_area                = 0.0
      urb_zone_n_cells             = 0
      urb_zone_qmax_total          = 0.0
      !
      ! Allocate per-cell state.
      !
      allocate(urban_drainage_zone_indices(np))
      allocate(urban_drainage_qmax(np))
      allocate(urban_drainage_backflow_coef(np))
      allocate(urban_drainage_cumulative_volume(np))
      urban_drainage_zone_indices      = 0
      urban_drainage_qmax              = 0.0
      urban_drainage_backflow_coef     = 0.0
      urban_drainage_cumulative_volume = 0.0
      !
      ! Stamp cells per zone. Polygons are cached per unique file so that
      ! multiple zones sharing a polygon file only trigger one file read.
      ! Within a file each polygon name is matched against urb_zone_name.
      !
      allocate(inside(np))
      last_file = ''
      !
      do iz = 1, nr_urban_drainage_zones
         !
         if (trim(urb_zone_polygon_file(iz)) == '') then
            write(logstr,'(a,a,a)')' Error ! Urban drainage zone "', trim(urb_zone_name(iz)), &
                 '" has no polygon_file'
            call stop_sfincs(trim(logstr), -1)
         endif
         !
         if (trim(urb_zone_polygon_file(iz)) /= trim(last_file)) then
            if (allocated(polygons)) then
               do ip = 1, size(polygons)
                  if (allocated(polygons(ip)%x)) deallocate(polygons(ip)%x)
                  if (allocated(polygons(ip)%y)) deallocate(polygons(ip)%y)
               enddo
               deallocate(polygons)
            endif
            call read_tek_polygons(trim(urb_zone_polygon_file(iz)), polygons, ierr)
            if (ierr /= 0) then
               write(logstr,'(a,a)')' Error ! Failed to read urban drainage polygon file ', &
                    trim(urb_zone_polygon_file(iz))
               call stop_sfincs(trim(logstr), -1)
            endif
            last_file = urb_zone_polygon_file(iz)
         endif
         !
         ! Find a polygon in this file whose name matches this zone's name.
         !
         ipoly = 0
         do ip = 1, size(polygons)
            if (trim(polygons(ip)%name) == trim(urb_zone_name(iz))) then
               ipoly = ip
               exit
            endif
         enddo
         !
         if (ipoly == 0) then
            write(logstr,'(a,a,a,a)')' Error ! No polygon named "', trim(urb_zone_name(iz)), &
                 '" found in file ', trim(urb_zone_polygon_file(iz))
            call stop_sfincs(trim(logstr), -1)
         endif
         !
         ! Test all active cell centers against this polygon.
         !
         inside = .false.
         call points_in_polygon_omp(z_xz, z_yz, np, polygons(ipoly), inside)
         !
         ! Last-zone-wins: overwrite zone_indices wherever inside is true.
         !
         do nm = 1, np
            if (inside(nm)) urban_drainage_zone_indices(nm) = iz
         enddo
         !
      enddo
      !
      if (allocated(polygons)) then
         do ip = 1, size(polygons)
            if (allocated(polygons(ip)%x)) deallocate(polygons(ip)%x)
            if (allocated(polygons(ip)%y)) deallocate(polygons(ip)%y)
         enddo
         deallocate(polygons)
      endif
      deallocate(inside)
      !
      ! Snap outfall coordinates to nearest active cell.
      !
      n_outfalls = 0
      !
      do iz = 1, nr_urban_drainage_zones
         !
         if (.not. urb_zone_include_outfall(iz)) cycle
         !
         nmq = find_quadtree_cell(urb_zone_outfall_x(iz), urb_zone_outfall_y(iz))
         if (nmq > 0) urban_drainage_outfall_index(iz) = index_sfincs_in_quadtree(nmq)
         !
         if (urban_drainage_outfall_index(iz) <= 0) then
            write(logstr,'(a,a,a)')' Warning : outfall for zone "', trim(urb_zone_name(iz)), &
                 '" could not be snapped to an active cell; zone contributions will be discarded'
            call write_log(logstr, 0)
         else
            n_outfalls = n_outfalls + 1
         endif
         !
      enddo
      !
      ! Precompute per-cell qmax and backflow coef. Done in two passes so
      ! that zones specified via max_outfall_rate can derive their
      ! design_precip from the now-known total zone area. This keeps the
      ! per-cell qmax formula uniform across both input styles.
      !
      ! Pass 1: accumulate area and cell count per zone.
      !
      n_cells_in_zones = 0
      !
      do nm = 1, np
         !
         iz = urban_drainage_zone_indices(nm)
         if (iz == 0) cycle
         !
         if (crsgeo) then
            area_nm = cell_area_m2(nm)
         else
            area_nm = cell_area(z_flags_iref(nm))
         endif
         !
         urb_zone_area(iz)    = urb_zone_area(iz) + area_nm
         urb_zone_n_cells(iz) = urb_zone_n_cells(iz) + 1
         n_cells_in_zones     = n_cells_in_zones + 1
         !
      enddo
      !
      ! Derive design_precip for zones that were given max_outfall_rate.
      ! design_precip [mm/hr] = max_outfall_rate / area * 1000 * 3600.
      !
      do iz = 1, nr_urban_drainage_zones
         !
         if (urb_zone_max_outfall_rate(iz) > 0.0) then
            !
            if (urb_zone_area(iz) <= 0.0) then
               write(logstr,'(a,a,a)')' Error ! urban_drainage_zone "', trim(urb_zone_name(iz)), &
                    '" has max_outfall_rate set but zero participating cells; cannot derive design_precip'
               call stop_sfincs(trim(logstr), -1)
            endif
            !
            urb_zone_design_precip(iz) = urb_zone_max_outfall_rate(iz) / urb_zone_area(iz) * 1000.0 * 3600.0
            !
         endif
         !
      enddo
      !
      ! Pass 2: compute per-cell qmax and backflow coefficient.
      !
      do nm = 1, np
         !
         iz = urban_drainage_zone_indices(nm)
         if (iz == 0) cycle
         !
         if (crsgeo) then
            area_nm = cell_area_m2(nm)
         else
            area_nm = cell_area(z_flags_iref(nm))
         endif
         !
         ! mm/hr -> m/s then m3/s
         !
         urban_drainage_qmax(nm) = urb_zone_design_precip(iz) * 1.0e-3 / 3600.0 * area_nm
         !
         io = urban_drainage_outfall_index(iz)
         !
         if (io > 0) then
            dh_min = urb_zone_dh_design_min(iz)
            dzb    = max(zb(nm) - zb(io), dh_min)
            urban_drainage_backflow_coef(nm) = urban_drainage_qmax(nm) / sqrt(dzb)
         endif
         !
         urb_zone_qmax_total(iz) = urb_zone_qmax_total(iz) + urban_drainage_qmax(nm)
         !
      enddo
      !
      write(logstr,'(a,i0,a,i0,a,i0,a)')' Info    : urban drainage: ', nr_urban_drainage_zones, &
           ' zone(s), ', n_cells_in_zones, ' cell(s) assigned, ', n_outfalls, ' outfall(s)'
      call write_log(logstr, 0)
      !
      if (n_cells_in_zones > 0) urban_drainage = .true.
      !
      call write_urban_drainage_log_summary()
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine update_urban_drainage(t, dt)
      !
      ! Per-time-step entry: accumulate signed discharges into qsrc for
      ! cells inside drainage zones, and deposit the summed per-zone flux at
      ! each zone's outfall cell.
      !
      ! Sign convention: qd > 0 means water leaves the cell (drains to the
      ! outfall). qsrc(nm) -= qd subtracts that flux from the cell and the
      ! same amount is added back at the outfall.
      !
      ! Called from: update_continuity (sfincs_continuity), once per time
      ! step, after update_src_structures.
      !
      use sfincs_data
      !
      implicit none
      !
      real*8, intent(in) :: t
      real*4, intent(in) :: dt
      !
      integer :: nm, iz, io
      real*4  :: dzs, qd, area_nm, h_cell
      !
      if (nr_urban_drainage_zones <= 0) return
      !
      !$acc kernels present(urban_drainage_q_outfall)
      urban_drainage_q_outfall = 0.0
      !$acc end kernels
      !
      !$acc parallel loop present( qsrc, zs, zb, cell_area, cell_area_m2, z_flags_iref, &
      !$acc                        urban_drainage_zone_indices, urban_drainage_outfall_index, &
      !$acc                        urban_drainage_qmax, urban_drainage_backflow_coef, &
      !$acc                        urban_drainage_q_outfall, urban_drainage_cumulative_volume, &
      !$acc                        urb_zone_h_threshold, urb_zone_check_valve )
      !$omp parallel do default(shared) &
      !$omp private(nm, iz, io, dzs, qd, area_nm, h_cell) schedule(static)
      do nm = 1, np
         !
         iz = urban_drainage_zone_indices(nm)
         if (iz == 0) cycle
         !
         io = urban_drainage_outfall_index(iz)
         if (io <= 0) cycle
         !
         dzs = zs(nm) - zs(io)
         !
         if (dzs > 0.0) then
            !
            ! Drain from cell. Gate on cell ponding depth above grate.
            !
            h_cell = zs(nm) - zb(nm)
            if (h_cell <= urb_zone_h_threshold(iz)) cycle
            !
            if (crsgeo) then
               area_nm = cell_area_m2(nm)
            else
               area_nm = cell_area(z_flags_iref(nm))
            endif
            !
            qd = min(urban_drainage_qmax(nm), h_cell * area_nm / dt)
            !
         else
            !
            ! Backflow from outfall. Blocked by a check valve.
            !
            if (urb_zone_check_valve(iz)) cycle
            !
            qd = -urban_drainage_backflow_coef(nm) * sqrt(-dzs)
            if (qd < -urban_drainage_qmax(nm)) qd = -urban_drainage_qmax(nm)
            !
         endif
         !
         ! qsrc(nm) is unique per iteration (loop is over nm), no race.
         ! The race is on urban_drainage_q_outfall(iz): multiple threads
         ! (or gangs on device) may process cells belonging to the same
         ! zone, so guard the zone-accumulator with atomic.
         !
         qsrc(nm) = qsrc(nm) - qd
         !
         !$acc atomic update
         !$omp atomic
         urban_drainage_q_outfall(iz) = urban_drainage_q_outfall(iz) + qd
         !
         urban_drainage_cumulative_volume(nm) = urban_drainage_cumulative_volume(nm) + qd * dt
         !
      enddo
      !$omp end parallel do
      !
      ! Second pass: add each zone's net flux back at the outfall cell.
      ! Atomic guards against multiple zones snapping to the same outfall
      ! cell (rare but possible).
      !
      !$acc parallel loop present( qsrc, urban_drainage_outfall_index, urban_drainage_q_outfall )
      do iz = 1, nr_urban_drainage_zones
         !
         io = urban_drainage_outfall_index(iz)
         if (io <= 0) cycle
         !
         !$acc atomic update
         qsrc(io) = qsrc(io) + urban_drainage_q_outfall(iz)
         !
      enddo
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine write_urban_drainage_log_summary()
      !
      ! Emit a one-block-per-zone description of every parsed urban drainage
      ! zone to the log file. Intended for operator review at init time.
      !
      ! Called from: initialize_urban_drainage (this module), once after
      ! cells have been stamped and per-zone totals have been accumulated.
      !
      implicit none
      !
      integer :: iz
      !
      if (nr_urban_drainage_zones <= 0) return
      !
      call write_log('------------------------------------------', 0)
      call write_log('Urban drainage zones', 0)
      call write_log('------------------------------------------', 0)
      !
      write(logstr,'(a,i0,a)')'Added ', nr_urban_drainage_zones, ' urban drainage zone(s)'
      call write_log(logstr, 0)
      call write_log('', 0)
      !
      do iz = 1, nr_urban_drainage_zones
         !
         write(logstr,'(a,i0,a)')'Zone ', iz, ':'
         call write_log(logstr, 0)
         !
         write(logstr,'(a,a)')         '  name:            ', trim(urb_zone_name(iz))
         call write_log(logstr, 0)
         !
         if (len_trim(urb_zone_type(iz)) > 0) then
            write(logstr,'(a,a)')      '  type:            ', trim(urb_zone_type(iz))
            call write_log(logstr, 0)
         endif
         !
         write(logstr,'(a,a)')         '  polygon_file:    ', trim(urb_zone_polygon_file(iz))
         call write_log(logstr, 0)
         !
         write(logstr,'(a,i0)')        '  cells_assigned:  ', urb_zone_n_cells(iz)
         call write_log(logstr, 0)
         !
         write(logstr,'(a,f0.1,a)')    '  area:            ', urb_zone_area(iz), ' (m2)'
         call write_log(logstr, 0)
         !
         if (urb_zone_max_outfall_rate(iz) > 0.0) then
            write(logstr,'(a,f0.4,a)') '  max_outfall_rate:', urb_zone_max_outfall_rate(iz), ' (m3/s)'
            call write_log(logstr, 0)
            write(logstr,'(a,f0.2,a)') '  design_precip:   ', urb_zone_design_precip(iz), &
                 ' (mm/hr, derived from max_outfall_rate)'
            call write_log(logstr, 0)
         else
            write(logstr,'(a,f0.2,a)') '  design_precip:   ', urb_zone_design_precip(iz), ' (mm/hr)'
            call write_log(logstr, 0)
         endif
         !
         write(logstr,'(a,f0.4,a)')    '  qmax_total:      ', urb_zone_qmax_total(iz), ' (m3/s)'
         call write_log(logstr, 0)
         !
         write(logstr,'(a,f0.3,a)')    '  h_threshold:     ', urb_zone_h_threshold(iz), ' (m)'
         call write_log(logstr, 0)
         !
         write(logstr,'(a,f0.3,a)')    '  dh_design_min:   ', urb_zone_dh_design_min(iz), ' (m)'
         call write_log(logstr, 0)
         !
         if (urb_zone_include_outfall(iz)) then
            write(logstr,'(a,f0.3,a,f0.3,a)')'  outfall:         (', urb_zone_outfall_x(iz), ', ', &
                 urb_zone_outfall_y(iz), ')'
            call write_log(logstr, 0)
            !
            if (urban_drainage_outfall_index(iz) > 0) then
               write(logstr,'(a,i0)')  '  outfall_index:   ', urban_drainage_outfall_index(iz)
               call write_log(logstr, 0)
            else
               call write_log('  outfall_index:   (no active cell snapped)', 0)
            endif
         else
            call write_log('  outfall:         (disabled)', 0)
         endif
         !
         if (urb_zone_check_valve(iz)) then
            call write_log('  check_valve:     true', 0)
         else
            call write_log('  check_valve:     false', 0)
         endif
         !
         call write_log('', 0)
         !
      enddo
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine read_urban_drainage(filename, ierr)
      !
      ! Parse the *.urb TOML file into the per-zone arrays.
      !
      ! Schema:
      !
      !    [[urban_drainage_zone]]
      !    name            = "area 1"            ! required, string (matches polygon name)
      !    type            = "drainage"          ! optional, free-form tag (reserved)
      !    polygon_file    = "zones.tek"         ! required
      !    outfall_x       = 950.0               ! required if include_outfall = true
      !    outfall_y       = 150.0               ! required if include_outfall = true
      !    design_precip   = 20.0                ! required if max_outfall_rate absent, mm/hr
      !    max_outfall_rate = 6.0                ! alternative to design_precip, m3/s total zone capacity
      !                                          ! exactly one of {design_precip, max_outfall_rate} must be given
      !    h_threshold     = 0.0                 ! optional, m (default 0.0)
      !    dh_design_min   = 0.1                 ! optional, m (default 0.1)
      !    include_outfall = true                ! optional (default true)
      !    check_valve     = true                ! optional (default false)
      !
      ! Called from: initialize_urban_drainage (this module).
      !
      use tomlf
      !
      implicit none
      !
      character(len=*), intent(in)  :: filename
      integer,          intent(out) :: ierr
      !
      type(toml_table), allocatable :: top
      type(toml_error), allocatable :: err
      type(toml_array), pointer     :: arr_zones
      type(toml_table), pointer     :: tbl_zone
      character(len=:), allocatable :: name_str, type_str, poly_str
      integer                       :: nz, i, stat
      real*4                        :: r4_tmp
      real(kind=8)                  :: r8_tmp
      logical                       :: l_tmp, found
      !
      ierr = 0
      !
      call toml_load(top, filename, error=err)
      if (allocated(err)) then
         write(logstr,'(a,a,a,a)')' Error ! Failed to parse TOML file ', trim(filename), ': ', &
              trim(err%message)
         call write_log(logstr, 1)
         ierr = 1
         return
      endif
      !
      if (.not. allocated(top)) then
         write(logstr,'(a,a)')' Error ! Could not load TOML file ', trim(filename)
         call write_log(logstr, 1)
         ierr = 1
         return
      endif
      !
      nullify(arr_zones)
      call get_value(top, 'urban_drainage_zone', arr_zones, requested=.false., stat=stat)
      !
      if (.not. associated(arr_zones)) then
         nr_urban_drainage_zones = 0
         return
      endif
      !
      if (.not. is_array_of_tables(arr_zones)) then
         write(logstr,'(a,a)')' Error ! Key "urban_drainage_zone" must be an array of tables in ', &
              trim(filename)
         call write_log(logstr, 1)
         ierr = 1
         return
      endif
      !
      nz = len(arr_zones)
      nr_urban_drainage_zones = nz
      !
      if (nz == 0) return
      !
      allocate(urb_zone_name(nz))
      allocate(urb_zone_type(nz))
      allocate(urb_zone_polygon_file(nz))
      allocate(urb_zone_outfall_x(nz))
      allocate(urb_zone_outfall_y(nz))
      allocate(urb_zone_design_precip(nz))
      allocate(urb_zone_max_outfall_rate(nz))
      allocate(urb_zone_h_threshold(nz))
      allocate(urb_zone_dh_design_min(nz))
      allocate(urb_zone_include_outfall(nz))
      allocate(urb_zone_check_valve(nz))
      !
      urb_zone_name            = ''
      urb_zone_type            = ''
      urb_zone_polygon_file    = ''
      urb_zone_outfall_x       = 0.0
      urb_zone_outfall_y       = 0.0
      urb_zone_design_precip    = 0.0
      urb_zone_max_outfall_rate = 0.0
      urb_zone_h_threshold      = 0.0
      urb_zone_dh_design_min   = 0.1
      urb_zone_include_outfall = .true.
      urb_zone_check_valve     = .false.
      !
      do i = 1, nz
         !
         nullify(tbl_zone)
         call get_value(arr_zones, i, tbl_zone, stat=stat)
         if (.not. associated(tbl_zone)) then
            write(logstr,'(a,i0,a)')' Error ! urban_drainage_zone entry ', i, ' is not a table'
            call write_log(logstr, 1)
            ierr = 1
            return
         endif
         !
         if (allocated(name_str)) deallocate(name_str)
         call get_value(tbl_zone, 'name', name_str, stat=stat)
         if (.not. allocated(name_str)) then
            write(logstr,'(a,i0)')' Error ! Missing required "name" in urban_drainage_zone entry ', i
            call write_log(logstr, 1)
            ierr = 1
            return
         endif
         urb_zone_name(i) = name_str
         !
         if (allocated(type_str)) deallocate(type_str)
         call get_value(tbl_zone, 'type', type_str, stat=stat)
         if (allocated(type_str)) urb_zone_type(i) = type_str
         !
         if (allocated(poly_str)) deallocate(poly_str)
         call get_value(tbl_zone, 'polygon_file', poly_str, stat=stat)
         if (.not. allocated(poly_str)) then
            write(logstr,'(a,a,a)')' Error ! Missing required "polygon_file" in urban_drainage_zone "', &
                 trim(urb_zone_name(i)), '"'
            call write_log(logstr, 1)
            ierr = 1
            return
         endif
         urb_zone_polygon_file(i) = poly_str
         !
         call get_value(tbl_zone, 'outfall_x', r8_tmp, stat=stat)
         if (stat == 0) urb_zone_outfall_x(i) = real(r8_tmp, 4)
         !
         call get_value(tbl_zone, 'outfall_y', r8_tmp, stat=stat)
         if (stat == 0) urb_zone_outfall_y(i) = real(r8_tmp, 4)
         !
         ! Exactly one of design_precip / max_outfall_rate must be given.
         ! has_key distinguishes "absent" from "present but 0.0", so a user
         ! who really wants a zero-capacity zone can still write it.
         !
         block
            logical :: has_precip, has_rate
            has_precip = tbl_zone%has_key('design_precip')
            has_rate   = tbl_zone%has_key('max_outfall_rate')
            if (has_precip .and. has_rate) then
               write(logstr,'(a,a,a)')' Error ! urban_drainage_zone "', trim(urb_zone_name(i)), &
                    '" has both "design_precip" and "max_outfall_rate"; specify only one'
               call write_log(logstr, 1)
               ierr = 1
               return
            endif
            if (.not. has_precip .and. .not. has_rate) then
               write(logstr,'(a,a,a)')' Error ! urban_drainage_zone "', trim(urb_zone_name(i)), &
                    '" needs "design_precip" (mm/hr) or "max_outfall_rate" (m3/s)'
               call write_log(logstr, 1)
               ierr = 1
               return
            endif
            if (has_precip) then
               call get_value(tbl_zone, 'design_precip', r8_tmp, stat=stat)
               urb_zone_design_precip(i) = real(r8_tmp, 4)
            else
               call get_value(tbl_zone, 'max_outfall_rate', r8_tmp, stat=stat)
               urb_zone_max_outfall_rate(i) = real(r8_tmp, 4)
            endif
         end block
         !
         call get_value(tbl_zone, 'h_threshold', r8_tmp, stat=stat)
         if (stat == 0) urb_zone_h_threshold(i) = real(r8_tmp, 4)
         !
         call get_value(tbl_zone, 'dh_design_min', r8_tmp, stat=stat)
         if (stat == 0) urb_zone_dh_design_min(i) = real(r8_tmp, 4)
         !
         call get_value(tbl_zone, 'include_outfall', l_tmp, stat=stat)
         if (stat == 0) urb_zone_include_outfall(i) = l_tmp
         !
         call get_value(tbl_zone, 'check_valve', l_tmp, stat=stat)
         if (stat == 0) urb_zone_check_valve(i) = l_tmp
         !
         ! Minimal sanity check on outfall: if include_outfall is true, outfall
         ! coords must be specified (non-zero default is a weak check; keep it
         ! simple by warning rather than failing — snap will catch bad values).
         !
         if (urb_zone_include_outfall(i)) then
            found = (urb_zone_outfall_x(i) /= 0.0 .or. urb_zone_outfall_y(i) /= 0.0)
            if (.not. found) then
               write(logstr,'(a,a,a)')' Warning : urban_drainage_zone "', trim(urb_zone_name(i)), &
                    '" has include_outfall = true but outfall_x, outfall_y both 0.0'
               call write_log(logstr, 0)
            endif
         endif
         !
         if (urb_zone_dh_design_min(i) <= 0.0) urb_zone_dh_design_min(i) = 0.1
         !
      enddo
      !
      ! Keep the compiler from warning about unused variables in case get_value
      ! signatures drift; r4_tmp is reserved for future per-zone scalars.
      !
      r4_tmp = 0.0
      if (r4_tmp < 0.0) continue
      !
   end subroutine
   !
end module sfincs_urban_drainage
