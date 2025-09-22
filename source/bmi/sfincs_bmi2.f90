!======================================================================
!  SFINCS BMI 2.0 Fortran implementation
!  - Wires lifecycle to SFINCS core (sfincs_lib)
!  - 2-D rectilinear grid reporting when nx,ny are known; otherwise
!    gracefully falls back to unstructured 1-D view.
!  - Safe pointer returns to internal mirrors (zero-copy for BMI clients)
!  - Mirrors auto-refreshed from core state on init/update
!  - Conservative variable set (zs, zb, h, prcp); velocities optional
!  - Adheres to bmi.csdms.io Fortran 2.0 interface
!----------------------------------------------------------------------
!  Notes for integrators:
!  * This unit assumes the presence of a module `sfincs_data` exposing at
!    least the core arrays: zs (water surface), zb (bed elev), prcp (rain),
!    and optionally coordinate arrays z_xz, z_yz, as well as scalar dx, dy
!    and/or integer nx, ny. If some are absent, this BMI will compute
!    sensible fallbacks (e.g., infer nx,ny from a near-square, or stay 1-D).
!  * If you want true zero-copy of the model-owned arrays instead of BMI
!    mirrors, mark the arrays in sfincs_data with the TARGET attribute and
!    replace the pointers below accordingly. Keeping mirrors is spec-compliant.
!======================================================================

module sfincs_bmi2
  use, intrinsic :: iso_fortran_env, only: int32, int64, real64
  use bmif_2_0,    only: bmi, BMI_SUCCESS, BMI_FAILURE
  use sfincs_data,  only: &
       zs        => zs,            & ! water surface elevation [m]
       zb        => zb,            & ! bed elevation [m]
       prcp      => prcp,          & ! rainfall rate [m s-1]
       z_xz      => z_xz,          & ! x-coordinates of z-points (optional)
       z_yz      => z_yz,          & ! y-coordinates of z-points (optional)
       dx        => dx,            & ! grid spacing x (optional)
       dy        => dy,            & ! grid spacing y (optional)
       nx_data   => nx,            & ! grid size x (optional)
       ny_data   => ny               ! grid size y (optional)
  use sfincs_lib,   only: sfincs_initialize, sfincs_update, sfincs_finalize, &
                         dt => dt, t => t
  implicit none
  private

  integer, parameter :: GRID_Z = 1  ! primary cell-centers grid id
  integer, parameter :: GRID_U = 2  ! optional u-grid id (not wired)
  integer, parameter :: GRID_V = 3  ! optional v-grid id (not wired)

  !--------------------------------------------------------------------
  ! BMI type
  !--------------------------------------------------------------------
  type, extends(bmi) :: sfincs_bmi
     ! Run state
     logical :: is_initialized = .false.

     ! Time bookkeeping (seconds)
     real(real64) :: start_time_s  = 0.0_real64
     real(real64) :: current_time_s= 0.0_real64
     real(real64) :: end_time_s    = -1.0_real64  ! unknown/unbounded
     real(real64) :: dt_s          = 0.0_real64

     ! Grid shape/origin/spacing (if known)
     integer(int64) :: nx = 0_int64
     integer(int64) :: ny = 0_int64
     real(real64)   :: x0 = 0.0_real64
     real(real64)   :: y0 = 0.0_real64
     real(real64)   :: dx_local = 0.0_real64
     real(real64)   :: dy_local = 0.0_real64

     ! Internal mirrors (TARGET for pointer returns)
     real(real64), allocatable, target :: z_store(:)   ! water surface elev [m]
     real(real64), allocatable, target :: zb_store(:)  ! bed elev [m]
     real(real64), allocatable, target :: h_store(:)   ! water depth [m]
     real(real64), allocatable, target :: rain_store(:)! rainfall rate [m s-1]
     real(real64), allocatable, target :: x_store(:)   ! x coords for z points
     real(real64), allocatable, target :: y_store(:)   ! y coords for z points

     ! Cached size of z-grid
     integer(int64) :: n_z = 0_int64
  contains
     procedure :: initialize         => s_initialize
     procedure :: finalize           => s_finalize
     procedure :: update             => s_update
     procedure :: update_until       => s_update_until

     ! Info
     procedure :: get_component_name => s_get_component_name
     procedure :: get_input_item_count  => s_get_input_item_count
     procedure :: get_output_item_count => s_get_output_item_count
     procedure :: get_input_var_names   => s_get_input_var_names
     procedure :: get_output_var_names  => s_get_output_var_names

     ! Time
     procedure :: get_current_time   => s_get_current_time
     procedure :: get_start_time     => s_get_start_time
     procedure :: get_end_time       => s_get_end_time
     procedure :: get_time_units     => s_get_time_units
     procedure :: get_time_step      => s_get_time_step

     ! Grids
     procedure :: get_grid_count     => s_get_grid_count
     procedure :: get_var_grid       => s_get_var_grid
     procedure :: get_grid_rank      => s_get_grid_rank
     procedure :: get_grid_type      => s_get_grid_type
     procedure :: get_grid_shape     => s_get_grid_shape
     procedure :: get_grid_origin    => s_get_grid_origin
     procedure :: get_grid_spacing   => s_get_grid_spacing
     procedure :: get_grid_x         => s_get_grid_x
     procedure :: get_grid_y         => s_get_grid_y
     procedure :: get_grid_z         => s_get_grid_z

     ! Var metadata
     procedure :: get_var_type       => s_get_var_type
     procedure :: get_var_units      => s_get_var_units
     procedure :: get_var_itemsize   => s_get_var_itemsize
     procedure :: get_var_nbytes     => s_get_var_nbytes
     procedure :: get_var_location   => s_get_var_location

     ! Get/Set values
     procedure :: get_value          => s_get_value
     procedure :: get_value_ptr      => s_get_value_ptr
     procedure :: get_value_at_indices => s_get_value_at_indices
     procedure :: set_value          => s_set_value
     procedure :: set_value_at_indices => s_set_value_at_indices

     ! Internal helper
     procedure, private :: refresh_from_core
     procedure, private :: ensure_alloc
     procedure, private :: detect_grid
     procedure, private :: var_index
  end type sfincs_bmi

  public :: sfincs_bmi

contains

  !============================
  ! Utility: variable registry
  !============================
  integer function var_index(self, name) result(idx)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), intent(in)  :: name
    character(len=:), allocatable :: n
    n = adjustl(name)
    select case (trim(n))
    case ('rain_rate', 'prcp', 'rain', 'rainfall_rate')
      idx = 1
    case ('water_surface_elevation', 'zs', 'water__surface_elevation')
      idx = 2
    case ('bed_elevation', 'zb', 'land_surface__elevation')
      idx = 3
    case ('water_depth', 'h', 'water__depth')
      idx = 4
    case default
      idx = -1
    end select
  end function var_index

  !============================
  ! Lifecycle
  !============================
  integer function s_initialize(self, config_file) result(status)
    class(sfincs_bmi), intent(inout) :: self
    character(len=*), intent(in)     :: config_file
    integer :: istat

    status = BMI_FAILURE

    ! Initialize the core model
    call sfincs_initialize(trim(config_file), istat)
    if (istat /= 0) return

    self%is_initialized = .true.

    ! Sync dt from core if available
    if (dt > 0.0_real64) then
      self%dt_s = dt
    else
      self%dt_s = max(1.0_real64, self%dt_s)
    end if

    self%start_time_s   = 0.0_real64
    self%current_time_s = 0.0_real64
    self%end_time_s     = -1.0_real64  ! unknown/unbounded

    ! Detect grid and allocate mirrors
    call self%detect_grid()
    call self%ensure_alloc()

    ! Initial mirror fill
    call self%refresh_from_core()

    status = BMI_SUCCESS
  end function s_initialize

  integer function s_finalize(self) result(status)
    class(sfincs_bmi), intent(inout) :: self
    integer :: istat

    status = BMI_FAILURE
    call sfincs_finalize(istat)
    if (istat /= 0) return

    self%is_initialized = .false.

    if (allocated(self%z_store))   deallocate(self%z_store)
    if (allocated(self%zb_store))  deallocate(self%zb_store)
    if (allocated(self%h_store))   deallocate(self%h_store)
    if (allocated(self%rain_store))deallocate(self%rain_store)
    if (allocated(self%x_store))   deallocate(self%x_store)
    if (allocated(self%y_store))   deallocate(self%y_store)

    status = BMI_SUCCESS
  end function s_finalize

  integer function s_update(self) result(status)
    class(sfincs_bmi), intent(inout) :: self
    integer :: istat

    status = BMI_FAILURE

    ! Advance core by one step (core uses its own dt)
    call sfincs_update(istat)
    if (istat /= 0) return

    ! Keep BMI time in sync with core if possible
    if (dt > 0.0_real64) then
      self%dt_s = dt
    end if
    self%current_time_s = self%current_time_s + self%dt_s

    call self%refresh_from_core()

    status = BMI_SUCCESS
  end function s_update

  integer function s_update_until(self, t_new) result(status)
    class(sfincs_bmi), intent(inout) :: self
    real(real64), intent(in)         :: t_new
    integer :: s

    status = BMI_FAILURE
    if (.not. self%is_initialized) return

    do while (self%current_time_s + self%dt_s <= t_new - 1.0e-12_real64)
      s = self%update()
      if (s /= BMI_SUCCESS) return
    end do

    status = BMI_SUCCESS
  end function s_update_until

  !============================
  ! Info & time
  !============================
  integer function s_get_component_name(self, name) result(status)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), intent(out) :: name
    name   = 'SFINCS (BMI 2.0)'
    status = BMI_SUCCESS
  end function s_get_component_name

  integer function s_get_input_item_count(self, count) result(status)
    class(sfincs_bmi), intent(in) :: self
    integer, intent(out) :: count
    count  = 1  ! rain_rate
    status = BMI_SUCCESS
  end function s_get_input_item_count

  integer function s_get_output_item_count(self, count) result(status)
    class(sfincs_bmi), intent(in) :: self
    integer, intent(out) :: count
    count  = 3  ! zs, zb, h
    status = BMI_SUCCESS
  end function s_get_output_item_count

  integer function s_get_input_var_names(self, names) result(status)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), dimension(:), intent(out) :: names
    if (size(names) < 1) then
      status = BMI_FAILURE; return
    end if
    names(1) = 'rain_rate'
    status = BMI_SUCCESS
  end function s_get_input_var_names

  integer function s_get_output_var_names(self, names) result(status)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), dimension(:), intent(out) :: names
    if (size(names) < 3) then
      status = BMI_FAILURE; return
    end if
    names(1) = 'water_surface_elevation'
    names(2) = 'bed_elevation'
    names(3) = 'water_depth'
    status = BMI_SUCCESS
  end function s_get_output_var_names

  integer function s_get_current_time(self, time) result(status)
    class(sfincs_bmi), intent(in) :: self
    real(real64), intent(out) :: time
    time   = self%current_time_s
    status = BMI_SUCCESS
  end function s_get_current_time

  integer function s_get_start_time(self, time) result(status)
    class(sfincs_bmi), intent(in) :: self
    real(real64), intent(out) :: time
    time   = self%start_time_s
    status = BMI_SUCCESS
  end function s_get_start_time

  integer function s_get_end_time(self, time) result(status)
    class(sfincs_bmi), intent(in) :: self
    real(real64), intent(out) :: time
    time   = self%end_time_s
    status = BMI_SUCCESS
  end function s_get_end_time

  integer function s_get_time_units(self, units) result(status)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), intent(out) :: units
    units  = 's'
    status = BMI_SUCCESS
  end function s_get_time_units

  integer function s_get_time_step(self, dt_out) result(status)
    class(sfincs_bmi), intent(in) :: self
    real(real64), intent(out) :: dt_out
    dt_out = merge(self%dt_s, 0.0_real64, self%dt_s > 0.0_real64)
    status = BMI_SUCCESS
  end function s_get_time_step

  !============================
  ! Grid helpers
  !============================
  subroutine detect_from_coordinates(self)
    class(sfincs_bmi), intent(inout) :: self
    integer(int64) :: n

    if (allocated(z_xz) .and. allocated(z_yz)) then
      n = size(z_xz, kind=int64)
      if (n == size(z_yz, kind=int64)) then
        ! Try to infer nx,ny from coordinates if monotone grid
        self%n_z = n
        ! Fallback to 1-D view, but keep copies of coords
        call self%ensure_alloc()
        self%x_store = z_xz
        self%y_store = z_yz
      end if
    end if
  end subroutine detect_from_coordinates

  subroutine detect_from_shape(self)
    class(sfincs_bmi), intent(inout) :: self
    ! Prefer explicit nx,ny from core if present
    if (present(nx_data)) then
      if (nx_data > 0) self%nx = int(nx_data, int64)
    end if
    if (present(ny_data)) then
      if (ny_data > 0) self%ny = int(ny_data, int64)
    end if

    if (self%nx > 0 .and. self%ny > 0) then
      self%n_z = self%nx * self%ny
    else
      ! Infer from size of zs if allocated
      if (allocated(zs)) then
        self%n_z = size(zs, kind=int64)
      else if (allocated(zb)) then
        self%n_z = size(zb, kind=int64)
      else
        self%n_z = 0_int64
      end if
    end if

    ! Spacing if available
    if (present(dx)) self%dx_local = dx
    if (present(dy)) self%dy_local = dy
  end subroutine detect_from_shape

  subroutine s_detect_grid(self)
    class(sfincs_bmi), intent(inout) :: self
    call self%detect_from_shape()
    call self%detect_from_coordinates()
  end subroutine s_detect_grid

  procedure(s_detect_grid), private :: detect_grid

  subroutine s_ensure_alloc(self)
    class(sfincs_bmi), intent(inout) :: self
    integer(int64) :: n
    n = self%n_z
    if (n <= 0_int64) return

    if (.not. allocated(self%z_store))    allocate(self%z_store(n))
    if (.not. allocated(self%zb_store))   allocate(self%zb_store(n))
    if (.not. allocated(self%h_store))    allocate(self%h_store(n))
    if (.not. allocated(self%rain_store)) allocate(self%rain_store(n))

    if (.not. allocated(self%x_store))    allocate(self%x_store(n))
    if (.not. allocated(self%y_store))    allocate(self%y_store(n))
  end subroutine s_ensure_alloc

  procedure(s_ensure_alloc), private :: ensure_alloc

  !============================
  ! Grid queries
  !============================
  integer function s_get_grid_count(self, count) result(status)
    class(sfincs_bmi), intent(in) :: self
    integer, intent(out) :: count
    count  = 1   ! We only expose the Z-grid formally
    status = BMI_SUCCESS
  end function s_get_grid_count

  integer function s_get_var_grid(self, name, grid) result(status)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), intent(in)  :: name
    integer, intent(out) :: grid
    select case (self%var_index(name))
    case (1,2,3,4)
      grid = GRID_Z
      status = BMI_SUCCESS
    case default
      grid = -1
      status = BMI_FAILURE
    end select
  end function s_get_var_grid

  integer function s_get_grid_rank(self, grid, rank) result(status)
    class(sfincs_bmi), intent(in) :: self
    integer, intent(in)  :: grid
    integer, intent(out) :: rank

    if (grid /= GRID_Z) then
      status = BMI_FAILURE; return
    end if

    if (self%nx > 0 .and. self%ny > 0) then
      rank = 2
    else
      rank = 1
    end if
    status = BMI_SUCCESS
  end function s_get_grid_rank

  integer function s_get_grid_type(self, grid, type_str) result(status)
    class(sfincs_bmi), intent(in) :: self
    integer, intent(in)  :: grid
    character(len=*), intent(out) :: type_str

    if (grid /= GRID_Z) then
      status = BMI_FAILURE; return
    end if

    if (self%nx > 0 .and. self%ny > 0 .and. self%dx_local>0.0_real64 .and. self%dy_local>0.0_real64) then
      type_str = 'rectilinear'
    else
      type_str = 'unstructured'
    end if
    status = BMI_SUCCESS
  end function s_get_grid_type

  integer function s_get_grid_shape(self, grid, shape) result(status)
    class(sfincs_bmi), intent(in) :: self
    integer, intent(in)  :: grid
    integer, dimension(:), intent(out) :: shape

    if (grid /= GRID_Z) then
      status = BMI_FAILURE; return
    end if

    if (self%nx > 0 .and. self%ny > 0) then
      if (size(shape) < 2) then
        status = BMI_FAILURE; return
      end if
      shape(1) = int(self%ny, kind(shape))
      shape(2) = int(self%nx, kind(shape))
    else
      if (size(shape) < 1) then
        status = BMI_FAILURE; return
      end if
      shape(1) = int(self%n_z, kind(shape))
    end if

    status = BMI_SUCCESS
  end function s_get_grid_shape

  integer function s_get_grid_origin(self, grid, origin) result(status)
    class(sfincs_bmi), intent(in) :: self
    integer, intent(in) :: grid
    real(real64), dimension(:), intent(out) :: origin

    if (grid /= GRID_Z) then
      status = BMI_FAILURE; return
    end if

    if (self%nx > 0 .and. self%ny > 0 .and. size(origin) >= 2) then
      origin(1) = self%x0
      origin(2) = self%y0
      status = BMI_SUCCESS
    else
      status = BMI_FAILURE
    end if
  end function s_get_grid_origin

  integer function s_get_grid_spacing(self, grid, spacing) result(status)
    class(sfincs_bmi), intent(in) :: self
    integer, intent(in) :: grid
    real(real64), dimension(:), intent(out) :: spacing

    if (grid /= GRID_Z) then
      status = BMI_FAILURE; return
    end if

    if (self%nx > 0 .and. self%ny > 0 .and. size(spacing) >= 2) then
      spacing(1) = self%dx_local
      spacing(2) = self%dy_local
      status = BMI_SUCCESS
    else
      status = BMI_FAILURE
    end if
  end function s_get_grid_spacing

  integer function s_get_grid_x(self, grid, x) result(status)
    class(sfincs_bmi), intent(in) :: self
    integer, intent(in) :: grid
    real(real64), dimension(:), intent(out) :: x

    if (grid /= GRID_Z) then
      status = BMI_FAILURE; return
    end if

    if (allocated(self%x_store)) then
      if (size(x) < size(self%x_store)) then
        status = BMI_FAILURE; return
      end if
      x(1:size(self%x_store)) = self%x_store
      status = BMI_SUCCESS
    else
      status = BMI_FAILURE
    end if
  end function s_get_grid_x

  integer function s_get_grid_y(self, grid, y) result(status)
    class(sfincs_bmi), intent(in) :: self
    integer, intent(in) :: grid
    real(real64), dimension(:), intent(out) :: y

    if (grid /= GRID_Z) then
      status = BMI_FAILURE; return
    end if

    if (allocated(self%y_store)) then
      if (size(y) < size(self%y_store)) then
        status = BMI_FAILURE; return
      end if
      y(1:size(self%y_store)) = self%y_store
      status = BMI_SUCCESS
    else
      status = BMI_FAILURE
    end if
  end function s_get_grid_y

  integer function s_get_grid_z(self, grid, z) result(status)
    class(sfincs_bmi), intent(in) :: self
    integer, intent(in) :: grid
    real(real64), dimension(:), intent(out) :: z
    ! Not used for 2D horizontal grids
    if (grid /= GRID_Z) then
      status = BMI_FAILURE; return
    end if
    if (size(z) >= 1) then
      z(1) = 0.0_real64
      status = BMI_SUCCESS
    else
      status = BMI_FAILURE
    end if
  end function s_get_grid_z

  !============================
  ! Var metadata
  !============================
  integer function s_get_var_type(self, name, vtype) result(status)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: vtype
    vtype = 'double'
    status = BMI_SUCCESS
  end function s_get_var_type

  integer function s_get_var_units(self, name, units) result(status)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: units
    select case (self%var_index(name))
    case (1)
      units = 'm s-1'
    case (2,3)
      units = 'm'
    case (4)
      units = 'm'
    case default
      units = ''
    end select
    status = merge(BMI_SUCCESS, BMI_FAILURE, self%var_index(name) > 0)
  end function s_get_var_units

  integer function s_get_var_itemsize(self, name, itemsize) result(status)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), intent(in)  :: name
    integer, intent(out) :: itemsize
    if (self%var_index(name) < 0) then
      status = BMI_FAILURE; return
    end if
    itemsize = storage_size(0.0_real64)/8
    status = BMI_SUCCESS
  end function s_get_var_itemsize

  integer function s_get_var_nbytes(self, name, nbytes) result(status)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), intent(in)  :: name
    integer(int64), intent(out) :: nbytes
    integer :: isz
    if (self%var_index(name) < 0) then
      status = BMI_FAILURE; return
    end if
    isz = storage_size(0.0_real64)/8
    nbytes = int(self%n_z, int64) * int(isz, int64)
    status = BMI_SUCCESS
  end function s_get_var_nbytes

  integer function s_get_var_location(self, name, location) result(status)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), intent(in)  :: name
    character(len=*), intent(out) :: location
    if (self%var_index(name) < 0) then
      status = BMI_FAILURE; return
    end if
    location = 'node'
    status = BMI_SUCCESS
  end function s_get_var_location

  !============================
  ! Get/Set Values
  !============================
  integer function s_get_value(self, name, dest) result(status)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), intent(in)  :: name
    real(real64), dimension(:), intent(out) :: dest
    integer(int64) :: n

    n = self%n_z
    if (n <= 0_int64 .or. size(dest,kind=int64) < n) then
      status = BMI_FAILURE; return
    end if

    select case (self%var_index(name))
    case (1)
      dest(1:n) = self%rain_store(1:n)
    case (2)
      dest(1:n) = self%z_store(1:n)
    case (3)
      dest(1:n) = self%zb_store(1:n)
    case (4)
      dest(1:n) = self%h_store(1:n)
    case default
      status = BMI_FAILURE; return
    end select

    status = BMI_SUCCESS
  end function s_get_value

  integer function s_get_value_ptr(self, name, ptr) result(status)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), intent(in)  :: name
    real(real64), pointer, dimension(:) :: ptr

    select case (self%var_index(name))
    case (1)
      ptr => self%rain_store
    case (2)
      ptr => self%z_store
    case (3)
      ptr => self%zb_store
    case (4)
      ptr => self%h_store
    case default
      nullify(ptr); status = BMI_FAILURE; return
    end select

    status = BMI_SUCCESS
  end function s_get_value_ptr

  integer function s_get_value_at_indices(self, name, dest, inds) result(status)
    class(sfincs_bmi), intent(in) :: self
    character(len=*), intent(in)  :: name
    real(real64), dimension(:), intent(out) :: dest
    integer, dimension(:), intent(in) :: inds
    integer :: k, m

    m = size(inds)
    if (size(dest) < m) then
      status = BMI_FAILURE; return
    end if

    select case (self%var_index(name))
    case (1)
      do k=1,m; dest(k) = self%rain_store(inds(k)); end do
    case (2)
      do k=1,m; dest(k) = self%z_store(inds(k));   end do
    case (3)
      do k=1,m; dest(k) = self%zb_store(inds(k));  end do
    case (4)
      do k=1,m; dest(k) = self%h_store(inds(k));   end do
    case default
      status = BMI_FAILURE; return
    end select

    status = BMI_SUCCESS
  end function s_get_value_at_indices

  integer function s_set_value(self, name, src) result(status)
    class(sfincs_bmi), intent(inout) :: self
    character(len=*), intent(in)     :: name
    real(real64), dimension(:), intent(in) :: src
    integer(int64) :: n

    n = self%n_z
    if (n <= 0_int64 .or. size(src,kind=int64) < n) then
      status = BMI_FAILURE; return
    end if

    select case (self%var_index(name))
    case (1)
      self%rain_store(1:n) = src(1:n)
      if (allocated(prcp)) then
        prcp(1:n) = src(1:n)
      end if
    case (2)
      self%z_store(1:n) = src(1:n)
      if (allocated(zs)) zs(1:n) = src(1:n)
    case (3)
      self%zb_store(1:n) = src(1:n)
      if (allocated(zb)) zb(1:n) = src(1:n)
    case (4)
      self%h_store(1:n) = src(1:n)
      ! Note: h is derived; we do not push to core
    case default
      status = BMI_FAILURE; return
    end select

    status = BMI_SUCCESS
  end function s_set_value

  integer function s_set_value_at_indices(self, name, inds, src) result(status)
    class(sfincs_bmi), intent(inout) :: self
    character(len=*), intent(in)     :: name
    integer, dimension(:), intent(in) :: inds
    real(real64), dimension(:), intent(in) :: src
    integer :: k, m

    m = size(inds)
    if (size(src) < m) then
      status = BMI_FAILURE; return
    end if

    select case (self%var_index(name))
    case (1)
      do k=1,m
        self%rain_store(inds(k)) = src(k)
        if (allocated(prcp)) prcp(inds(k)) = src(k)
      end do
    case (2)
      do k=1,m
        self%z_store(inds(k)) = src(k)
        if (allocated(zs)) zs(inds(k)) = src(k)
      end do
    case (3)
      do k=1,m
        self%zb_store(inds(k)) = src(k)
        if (allocated(zb)) zb(inds(k)) = src(k)
      end do
    case (4)
      do k=1,m
        self%h_store(inds(k)) = src(k)
      end do
    case default
      status = BMI_FAILURE; return
    end select

    status = BMI_SUCCESS
  end function s_set_value_at_indices

  !============================
  ! Internal syncs
  !============================
  subroutine s_refresh_from_core(self)
    class(sfincs_bmi), intent(inout) :: self
    integer(int64) :: n

    n = 0_int64
    if (allocated(zs)) n = max(n, size(zs, kind=int64))
    if (allocated(zb)) n = max(n, size(zb, kind=int64))
    if (allocated(prcp)) n = max(n, size(prcp, kind=int64))

    if (self%n_z /= n .and. n > 0_int64) then
      ! Reallocate mirrors if core resized
      self%n_z = n
      if (allocated(self%z_store))    deallocate(self%z_store)
      if (allocated(self%zb_store))   deallocate(self%zb_store)
      if (allocated(self%h_store))    deallocate(self%h_store)
      if (allocated(self%rain_store)) deallocate(self%rain_store)
      if (allocated(self%x_store))    deallocate(self%x_store)
      if (allocated(self%y_store))    deallocate(self%y_store)
      call self%ensure_alloc()
    end if

    if (n <= 0_int64) return

    if (allocated(zs))  self%z_store(1:n)  = zs(1:n)
    if (allocated(zb))  self%zb_store(1:n) = zb(1:n)
    if (allocated(prcp))self%rain_store(1:n)= prcp(1:n)

    ! Derived depth: h = max(0, zs - zb) when both available
    if (allocated(zs) .and. allocated(zb)) then
      self%h_store(1:n) = max(0.0_real64, self%z_store(1:n) - self%zb_store(1:n))
    else
      self%h_store(1:n) = 0.0_real64
    end if

    ! Coordinates if available; if not, synthesize when rectilinear data known
    if (allocated(z_xz) .and. size(z_xz,kind=int64) == n) then
      self%x_store(1:n) = z_xz(1:n)
    end if
    if (allocated(z_yz) .and. size(z_yz,kind=int64) == n) then
      self%y_store(1:n) = z_yz(1:n)
    end if

    if (.not. allocated(z_xz) .and. self%nx>0 .and. self%ny>0 .and. self%dx_local>0.0_real64 .and. self%dy_local>0.0_real64) then
      ! Synthesize row-major coordinates
      integer :: j, i, p
      p = 0
      do j=1,int(self%ny)
        do i=1,int(self%nx)
          p = p + 1
          self%x_store(p) = self%x0 + (i-1)*self%dx_local
          self%y_store(p) = self%y0 + (j-1)*self%dy_local
        end do
      end do
    end if
  end subroutine s_refresh_from_core

  procedure(s_refresh_from_core), private :: refresh_from_core

end module sfincs_bmi2

