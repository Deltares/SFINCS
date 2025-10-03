module sfincs_bmi2
  use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double
  use bmif_2_0,   only: bmi, BMI_SUCCESS, BMI_FAILURE, BMI_MAX_COMPONENT_NAME
  use sfincs_data, only: &
       zs, zb, prcp, &                ! primary state (water level, bed, rainfall)
       z_xz, z_yz                     ! z-grid coordinates
  implicit none

  private
  public :: sfincs_bmi

  !========================
  ! Module-scope mirror buffers (TARGET so BMI can pointer-alias them)
  !========================
  integer,        save :: nz = 0, nu = 0, nv = 0
  real(c_double), save, target, allocatable :: zb_store(:)   ! bed elevation [m]
  real(c_double), save, target, allocatable :: z_store(:)    ! water surface elevation [m]
  real(c_double), save, target, allocatable :: h_store(:)    ! water depth [m] (derived: max(0, zs - zb))
  real(c_double), save, target, allocatable :: un_store(:)   ! u-velocity [m/s] (no source -> zero, size 0)
  real(c_double), save, target, allocatable :: vn_store(:)   ! v-velocity [m/s] (no source -> zero, size 0)
  real(c_double), save, target, allocatable :: rain_store(:) ! rainfall [m s-1] (mirror of prcp)

  real(c_double), save, target, allocatable :: xz_store(:), yz_store(:)
  ! u/v grid coords not available in sfincs_data → leave empty mirrors
  real(c_double), save, target, allocatable :: xu_store(:), yu_store(:)
  real(c_double), save, target, allocatable :: xv_store(:), yv_store(:)

  ! Time bookkeeping (no dt exported from sfincs_data; keep internal clock)
  real(c_double), save :: start_time_s   = 0.0d0
  real(c_double), save :: end_time_s     = 0.0d0
  real(c_double), save :: current_time_s = 0.0d0
  real(c_double), save :: dt_s           = 60.0d0   ! default; adjust as you wish

  type, extends(bmi) :: sfincs_bmi
    character(len=:), pointer :: component_name => null()
    character(len=:), pointer :: time_units     => null()
    character(len=:), pointer :: input_names(:) => null()
    character(len=:), pointer :: output_names(:)=> null()

    ! BMI-visible state (aliases to the TARGET mirrors above)
    real(c_double), pointer :: bed(:) => null()
    real(c_double), pointer :: z(:)  => null()
    real(c_double), pointer :: h(:)  => null()
    real(c_double), pointer :: un(:) => null()
    real(c_double), pointer :: vn(:) => null()

    real(c_double), pointer :: xz(:) => null()
    real(c_double), pointer :: yz(:) => null()
    real(c_double), pointer :: xu(:) => null()
    real(c_double), pointer :: yu(:) => null()
    real(c_double), pointer :: xv(:) => null()
    real(c_double), pointer :: yv(:) => null()

    real(c_double), pointer :: rain(:) => null()
  contains
    ! Lifecycle
    procedure :: initialize                 => sfincs_initialize
    procedure :: update                     => sfincs_update
    procedure :: update_until               => sfincs_update_until
    procedure :: finalize                   => sfincs_finalize

    ! Names
    procedure :: get_component_name         => sfincs_get_component_name
    procedure :: get_input_item_count       => sfincs_get_input_item_count
    procedure :: get_output_item_count      => sfincs_get_output_item_count
    procedure :: get_input_var_names        => sfincs_get_input_var_names
    procedure :: get_output_var_names       => sfincs_get_output_var_names

    ! Time
    procedure :: get_start_time             => sfincs_get_start_time
    procedure :: get_end_time               => sfincs_get_end_time
    procedure :: get_current_time           => sfincs_get_current_time
    procedure :: get_time_step              => sfincs_get_time_step
    procedure :: get_time_units             => sfincs_get_time_units

    ! Vars ↔ grids + var metadata
    procedure :: get_var_grid               => sfincs_get_var_grid
    procedure :: get_var_type               => sfincs_get_var_type
    procedure :: get_var_units              => sfincs_get_var_units
    procedure :: get_var_itemsize           => sfincs_get_var_itemsize
    procedure :: get_var_nbytes             => sfincs_get_var_nbytes
    procedure :: get_var_location           => sfincs_get_var_location

    ! Grid info
    procedure :: get_grid_type              => sfincs_get_grid_type
    procedure :: get_grid_rank              => sfincs_get_grid_rank
    procedure :: get_grid_size              => sfincs_get_grid_size
    procedure :: get_grid_shape             => sfincs_get_grid_shape
    procedure :: get_grid_spacing           => sfincs_get_grid_spacing
    procedure :: get_grid_origin            => sfincs_get_grid_origin
    procedure :: get_grid_x                 => sfincs_get_grid_x
    procedure :: get_grid_y                 => sfincs_get_grid_y
    procedure :: get_grid_z                 => sfincs_get_grid_z

    ! Optional connectivity
    procedure :: get_grid_node_count        => sfincs_get_grid_node_count
    procedure :: get_grid_edge_count        => sfincs_get_grid_edge_count
    procedure :: get_grid_face_count        => sfincs_get_grid_face_count
    procedure :: get_grid_nodes_per_face    => sfincs_get_grid_nodes_per_face
    procedure :: get_grid_edge_nodes        => sfincs_get_grid_edge_nodes
    procedure :: get_grid_face_edges        => sfincs_get_grid_face_edges
    procedure :: get_grid_face_nodes        => sfincs_get_grid_face_nodes

    ! Value getters/setters
    procedure :: get_value_int              => sfincs_get_value_int
    procedure :: get_value_float            => sfincs_get_value_float
    procedure :: get_value_double           => sfincs_get_value_double

    procedure :: get_value_ptr_int          => sfincs_get_value_ptr_int
    procedure :: get_value_ptr_float        => sfincs_get_value_ptr_float
    procedure :: get_value_ptr_double       => sfincs_get_value_ptr_double

    procedure :: get_value_at_indices_int   => sfincs_get_value_at_indices_int
    procedure :: get_value_at_indices_float => sfincs_get_value_at_indices_float
    procedure :: get_value_at_indices_double=> sfincs_get_value_at_indices_double

    procedure :: set_value_int              => sfincs_set_value_int
    procedure :: set_value_float            => sfincs_set_value_float
    procedure :: set_value_double           => sfincs_set_value_double

    procedure :: set_value_at_indices_int   => sfincs_set_value_at_indices_int
    procedure :: set_value_at_indices_float => sfincs_set_value_at_indices_float
    procedure :: set_value_at_indices_double=> sfincs_set_value_at_indices_double
  end type sfincs_bmi

contains

!----------------------------------------------------------------------
! Pull the latest arrays from sfincs_data into the BMI-facing mirrors.
! This makes the BMI succeed with size=0 when SFINCS hasn't allocated yet,
! and use real SFINCS arrays as soon as they exist.
!----------------------------------------------------------------------
subroutine refresh_from_sfincs()
  use, intrinsic :: iso_c_binding, only: c_double
  implicit none
  integer :: n

  ! -- water surface (zs) --
  if (allocated(zs)) then
    n = size(zs)
    if (.not. allocated(z_store) .or. size(z_store) /= n) then
      if (allocated(z_store)) deallocate(z_store)
      allocate(z_store(n))
    end if
    z_store = real(zs, c_double)
    nz = n
  else
    if (allocated(z_store)) deallocate(z_store)
    nz = 0
  end if

  ! -- bed elevation (zb) --
  if (allocated(zb)) then
    n = size(zb)
    if (.not. allocated(zb_store) .or. size(zb_store) /= n) then
      if (allocated(zb_store)) deallocate(zb_store)
      allocate(zb_store(n))
    end if
    zb_store = real(zb, c_double)
  else
    if (allocated(zb_store)) deallocate(zb_store)
  end if

  ! -- depth from zs & zb (clip at 0) --
  if (allocated(zs) .and. allocated(zb)) then
    n = size(zs)
    if (.not. allocated(h_store) .or. size(h_store) /= n) then
      if (allocated(h_store)) deallocate(h_store)
      allocate(h_store(n))
    end if
    h_store = max(0.0d0, real(zs, c_double) - real(zb, c_double))
  else
    if (allocated(h_store)) h_store = 0.0d0
  end if

  ! -- rainfall (prcp) --
  if (allocated(prcp)) then
    n = size(prcp)
    if (.not. allocated(rain_store) .or. size(rain_store) /= n) then
      if (allocated(rain_store)) deallocate(rain_store)
      allocate(rain_store(n))
    end if
    rain_store = real(prcp, c_double)
  else
    if (allocated(rain_store)) deallocate(rain_store)
  end if

  ! -- z-grid coordinates --
  if (allocated(z_xz)) then
    n = size(z_xz)
    if (.not. allocated(xz_store) .or. size(xz_store) /= n) then
      if (allocated(xz_store)) deallocate(xz_store)
      allocate(xz_store(n))
    end if
    xz_store = real(z_xz, c_double)
    nz = max(nz, n)
  else
    if (allocated(xz_store)) deallocate(xz_store)
  end if

  if (allocated(z_yz)) then
    n = size(z_yz)
    if (.not. allocated(yz_store) .or. size(yz_store) /= n) then
      if (allocated(yz_store)) deallocate(yz_store)
      allocate(yz_store(n))
    end if
    yz_store = real(z_yz, c_double)
  else
    if (allocated(yz_store)) deallocate(yz_store)
  end if

  ! -- u/v not exposed from sfincs_data; keep mirrors empty
  nu = 0; nv = 0
  if (allocated(xu_store)) deallocate(xu_store)
  if (allocated(yu_store)) deallocate(yu_store)
  if (allocated(xv_store)) deallocate(xv_store)
  if (allocated(yv_store)) deallocate(yv_store)
  if (allocated(un_store)) deallocate(un_store)
  if (allocated(vn_store)) deallocate(vn_store)
end subroutine refresh_from_sfincs

!======================================================================
!                       LIFECYCLE
!======================================================================
function sfincs_initialize(this, config_file) result(status)
  use, intrinsic :: iso_c_binding, only: c_double
  implicit none
  class(sfincs_bmi), intent(out) :: this
  character(len=*),  intent(in)  :: config_file
  integer                        :: status
  integer :: nmini, i

  ! -- names and I/O lists
  if (.not. associated(this%component_name)) then
    allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%component_name)
    this%component_name = 'SFINCS-BMI'
  end if
  if (.not. associated(this%time_units)) then
    allocate(character(len=16) :: this%time_units)
    this%time_units = 's'
  end if
  if (.not. associated(this%input_names)) then
    allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%input_names(1))
    this%input_names(1) = 'rain_rate'
  end if
  if (.not. associated(this%output_names)) then
    allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%output_names(3))
    this%output_names = [ character(len=BMI_MAX_COMPONENT_NAME) :: &
         'water_surface_elevation', 'water_depth', 'bed_elevation' ]
  end if

  ! -- pull arrays from SFINCS if already allocated
  call refresh_from_sfincs()

  ! -- provide a tiny fallback if SFINCS hasn't filled anything yet
  if (nz <= 0) then
    nmini = 4
    nz = nmini

    if (allocated(z_store))  deallocate(z_store)
    if (allocated(h_store))  deallocate(h_store)
    if (allocated(zb_store)) deallocate(zb_store)
    if (allocated(xz_store)) deallocate(xz_store)
    if (allocated(yz_store)) deallocate(yz_store)
    if (allocated(rain_store)) deallocate(rain_store)

    allocate(z_store(nz), h_store(nz), zb_store(nz), xz_store(nz), yz_store(nz), rain_store(nz))

    do i = 1, nz
      xz_store(i)   = dble(i-1)
      yz_store(i)   = 0.0d0
      zb_store(i)   = 0.0d0
      z_store(i)    = 0.0d0
      h_store(i)    = max(0.0d0, z_store(i) - zb_store(i))
      rain_store(i) = 0.0d0
    end do

    nu = 0; nv = 0
  end if

  ! -- (re)associate BMI pointers
  if (allocated(z_store))   this%z   => z_store
  if (allocated(h_store))   this%h   => h_store
  if (allocated(xz_store))  this%xz  => xz_store
  if (allocated(yz_store))  this%yz  => yz_store
  if (allocated(rain_store))this%rain=> rain_store
  ! u/v remain disassociated (not provided by sfincs_data)

  ! -- time bookkeeping (simple defaults)
  start_time_s   = 0.0d0
  current_time_s = start_time_s
  end_time_s     = 0.0d0
  dt_s           = 60.0d0

  status = BMI_SUCCESS
end function sfincs_initialize

function sfincs_update(this) result(status)
  class(sfincs_bmi), intent(inout) :: this
  integer                          :: status
  call refresh_from_sfincs()
  current_time_s = current_time_s + dt_s
  status = BMI_SUCCESS
end function sfincs_update

function sfincs_update_until(this, time) result(status)
  class(sfincs_bmi), intent(inout) :: this
  double precision,   intent(in)   :: time
  integer                          :: status
  do while (current_time_s + dt_s <= time)
    status = sfincs_update(this)
    if (status /= BMI_SUCCESS) return
  end do
  status = BMI_SUCCESS
end function sfincs_update_until

function sfincs_finalize(this) result(status)
  class(sfincs_bmi), intent(inout) :: this
  integer                          :: status

  if (allocated(z_store))  deallocate(z_store)
  if (allocated(h_store))  deallocate(h_store)
  if (allocated(un_store)) deallocate(un_store)
  if (allocated(vn_store)) deallocate(vn_store)

  if (allocated(xz_store)) deallocate(xz_store)
  if (allocated(yz_store)) deallocate(yz_store)
  if (allocated(xu_store)) deallocate(xu_store)
  if (allocated(yu_store)) deallocate(yu_store)
  if (allocated(xv_store)) deallocate(xv_store)
  if (allocated(yv_store)) deallocate(yv_store)

  if (associated(this%component_name)) deallocate(this%component_name)
  if (associated(this%time_units))     deallocate(this%time_units)
  if (associated(this%input_names))    deallocate(this%input_names)
  if (associated(this%output_names))   deallocate(this%output_names)

  if (allocated(rain_store)) deallocate(rain_store)

  status = BMI_SUCCESS
end function sfincs_finalize

!======================================================================
!                 COMPONENT / I/O VARIABLE NAMES
!======================================================================
function sfincs_get_component_name(this, name) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*), pointer, intent(out) :: name
  integer :: status
  if (associated(this%component_name)) then
    name => this%component_name
    status = BMI_SUCCESS
  else
    nullify(name)
    status = BMI_FAILURE
  end if
end function sfincs_get_component_name

function sfincs_get_input_var_names(this, names) result(status)
  class(sfincs_bmi), intent(in)          :: this
  character(len=*), pointer, intent(out) :: names(:)
  integer                                 :: status
  if (associated(this%input_names)) then
    names => this%input_names
    status = BMI_SUCCESS
  else
    nullify(names)
    status = BMI_FAILURE
  end if
end function sfincs_get_input_var_names

function sfincs_get_output_var_names(this, names) result(status)
  class(sfincs_bmi), intent(in)          :: this
  character(len=*), pointer, intent(out) :: names(:)
  integer                                 :: status
  if (associated(this%output_names)) then
    names => this%output_names
    status = BMI_SUCCESS
  else
    nullify(names)
    status = BMI_FAILURE
  end if
end function sfincs_get_output_var_names

function sfincs_get_input_item_count(this, count) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(out):: count
  integer                       :: status
  if (associated(this%input_names)) then
    count = size(this%input_names)
  else
    count = 0
  end if
  status = BMI_SUCCESS
end function sfincs_get_input_item_count

function sfincs_get_output_item_count(this, count) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(out):: count
  integer                       :: status
  if (associated(this%output_names)) then
    count = size(this%output_names)
  else
    count = 0
  end if
  status = BMI_SUCCESS
end function sfincs_get_output_item_count

!======================================================================
!                               TIME
!======================================================================
function sfincs_get_start_time(this, time) result(status)
  class(sfincs_bmi), intent(in) :: this
  double precision, intent(out) :: time
  integer                       :: status
  time = start_time_s
  status = BMI_SUCCESS
end function sfincs_get_start_time

function sfincs_get_end_time(this, time) result(status)
  class(sfincs_bmi), intent(in) :: this
  double precision, intent(out) :: time
  integer                       :: status
  time = end_time_s
  status = BMI_SUCCESS
end function sfincs_get_end_time

function sfincs_get_current_time(this, time) result(status)
  class(sfincs_bmi), intent(in) :: this
  double precision, intent(out) :: time
  integer                       :: status
  time = current_time_s
  status = BMI_SUCCESS
end function sfincs_get_current_time

function sfincs_get_time_step(this, time_step) result(status)
  class(sfincs_bmi), intent(in) :: this
  double precision, intent(out) :: time_step
  integer                       :: status
  time_step = dt_s
  status = BMI_SUCCESS
end function sfincs_get_time_step

function sfincs_get_time_units(this, units) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(out):: units
  integer                       :: status
  units = ' '
  units(1:1) = 's'
  status = BMI_SUCCESS
end function sfincs_get_time_units

!======================================================================
!                   VAR ↔ GRID + VAR METADATA
!======================================================================
function sfincs_get_var_grid(this, name, grid) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  integer,           intent(out):: grid
  integer                       :: status
  select case (trim(name))
  case('water_surface_elevation','water_depth','bed_elevation','rain_rate'); grid = 1
  case('velocity_x');                                                       grid = 2
  case('velocity_y');                                                       grid = 3
  case default; grid = -1; status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_var_grid

function sfincs_get_var_type(this, name, type) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  character(len=*),  intent(out):: type
  integer                       :: status
  type = ' '
  type(1:6) = 'double'
  status = BMI_SUCCESS
end function sfincs_get_var_type

function sfincs_get_var_units(this, name, units) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  character(len=*),  intent(out):: units
  integer                       :: status
  select case (trim(name))
  case('water_surface_elevation','water_depth','bed_elevation')
    call assign_trim(units, 'm')
  case('rain_rate')
    call assign_trim(units, 'm s-1')
  case('velocity_x','velocity_y')
    call assign_trim(units, 'm s-1')
  case default
    call assign_trim(units, '1'); status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_var_units

function sfincs_get_var_itemsize(this, name, size) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  integer,           intent(out):: size
  integer                       :: status
  size = 8
  status = BMI_SUCCESS
end function sfincs_get_var_itemsize

function sfincs_get_var_nbytes(this, name, nbytes) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  integer,           intent(out):: nbytes
  integer :: status, n
  n = var_size(this, name)
  if (n < 0) then
    nbytes = 0
    status = BMI_SUCCESS   ! <- was FAILURE; make it SUCCESS with zero bytes
    return
  end if
  nbytes = 8*n
  status = BMI_SUCCESS
end function


function sfincs_get_var_location(this, name, location) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  character(len=*),  intent(out):: location
  integer                       :: status
  select case (trim(name))
  case('water_surface_elevation','water_depth','bed_elevation','rain_rate'); call assign_trim(location,'node')
  case('velocity_x','velocity_y');                                           call assign_trim(location,'edge')
  case default; call assign_trim(location,'unknown'); status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_var_location

!======================================================================
!                            GRIDS
!======================================================================
function sfincs_get_grid_type(this, grid, type) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  character(len=*),  intent(out):: type
  integer                       :: status
  call assign_trim(type, 'unstructured')
  status = BMI_SUCCESS
end function sfincs_get_grid_type

function sfincs_get_grid_rank(this, grid, rank) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: rank
  integer                       :: status
  rank = 1
  status = BMI_SUCCESS
end function sfincs_get_grid_rank

function sfincs_get_grid_size(this, grid, size) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: size
  integer                       :: status
  select case (grid)
  case (1); size = nz
  case (2); size = nu
  case (3); size = nv
  case default;   size = 0; status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_grid_size

function sfincs_get_grid_shape(this, grid, shape) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: shape(:)
  integer                       :: status
  if (size(shape) >= 1) then
    select case (grid)
    case (1); shape(1) = nz
    case (2); shape(1) = nu
    case (3); shape(1) = nv
    case default;   status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  else
    status = BMI_FAILURE
  end if
end function sfincs_get_grid_shape

function sfincs_get_grid_spacing(this, grid, spacing) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  double precision,  intent(out):: spacing(:)
  integer                       :: status
  if (size(spacing) >= 1) spacing(1) = 0.0d0
  status = BMI_SUCCESS
end function sfincs_get_grid_spacing

function sfincs_get_grid_origin(this, grid, origin) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  double precision,  intent(out):: origin(:)
  integer                       :: status
  if (size(origin) >= 1) origin(1) = 0.0d0
  status = BMI_SUCCESS
end function sfincs_get_grid_origin

function sfincs_get_grid_x(this, grid, x) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  double precision,  intent(out):: x(:)
  integer                       :: status
  integer :: n
  select case (grid)
  case (1); n = min(size(x), nz); if (n>0 .and. associated(this%xz)) x(1:n) = this%xz(1:n)
  case (2); n = 0   ! no u-grid coords available
  case (3); n = 0   ! no v-grid coords available
  case default; status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_grid_x

function sfincs_get_grid_y(this, grid, y) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  double precision,  intent(out):: y(:)
  integer                       :: status
  integer :: n
  select case (grid)
  case (1); n = min(size(y), nz); if (n>0 .and. associated(this%yz)) y(1:n) = this%yz(1:n)
  case (2); n = 0
  case (3); n = 0
  case default; status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_grid_y

function sfincs_get_grid_z(this, grid, z) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  double precision,  intent(out):: z(:)
  integer                       :: status
  if (size(z) >= 1) z(1) = 0.0d0
  status = BMI_SUCCESS
end function sfincs_get_grid_z

function sfincs_get_grid_node_count(this, grid, count) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: count
  integer                       :: status
  select case(grid)
  case(1); count = nz
  case(2); count = nu
  case(3); count = nv
  case default;  count = 0; status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_grid_node_count

function sfincs_get_grid_edge_count(this, grid, count) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: count
  integer                       :: status
  count = 0
  status = BMI_SUCCESS
end function sfincs_get_grid_edge_count

function sfincs_get_grid_face_count(this, grid, count) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: count
  integer                       :: status
  count = 0
  status = BMI_SUCCESS
end function sfincs_get_grid_face_count

function sfincs_get_grid_nodes_per_face(this, grid, nodes_per_face) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: nodes_per_face(:)
  integer                       :: status
  if (size(nodes_per_face) > 0) nodes_per_face = 0
  status = BMI_SUCCESS
end function sfincs_get_grid_nodes_per_face

function sfincs_get_grid_edge_nodes(this, grid, edge_nodes) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: edge_nodes(:)
  integer                       :: status
  if (size(edge_nodes) > 0) edge_nodes = 0
  status = BMI_SUCCESS
end function sfincs_get_grid_edge_nodes

function sfincs_get_grid_face_edges(this, grid, face_edges) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: face_edges(:)
  integer                       :: status
  if (size(face_edges) > 0) face_edges = 0
  status = BMI_SUCCESS
end function sfincs_get_grid_face_edges

function sfincs_get_grid_face_nodes(this, grid, face_nodes) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: face_nodes(:)
  integer                       :: status
  if (size(face_nodes) > 0) face_nodes = 0
  status = BMI_SUCCESS
end function sfincs_get_grid_face_nodes

!======================================================================
!                        VALUE GETTERS
!======================================================================
function sfincs_get_value_int(this, name, dest) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  integer,           intent(inout) :: dest(:)
  integer                       :: status
  status = BMI_FAILURE
end function sfincs_get_value_int

function sfincs_get_value_float(this, name, dest) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  real,              intent(inout) :: dest(:)
  integer                       :: status
  integer :: n
  select case (trim(name))
  case('water_surface_elevation')
    n = min(size(dest), size(this%z));  if (n>0 .and. associated(this%z))  dest(1:n) = real(this%z(1:n),  c_float)
  case('water_depth')
    n = min(size(dest), size(this%h));  if (n>0 .and. associated(this%h))  dest(1:n) = real(this%h(1:n),  c_float)
  case('bed_elevation')
    if (allocated(zb)) then
      n = min(size(dest), size(zb)); if (n>0) dest(1:n) = real(real(zb(1:n),c_double), c_float)
    else
      status = BMI_FAILURE; return
    end if
  case('velocity_x')
    if (associated(this%un)) then
      n = min(size(dest), size(this%un)); if (n>0) dest(1:n) = real(this%un(1:n), c_float)
    else
      n = 0
    end if
  case('velocity_y')
    if (associated(this%vn)) then
      n = min(size(dest), size(this%vn)); if (n>0) dest(1:n) = real(this%vn(1:n), c_float)
    else
      n = 0
    end if
  case('rain_rate')
    if (associated(this%rain)) then
      n = min(size(dest), size(this%rain)); if (n>0) dest(1:n) = real(this%rain(1:n), c_float)
    else
      status = BMI_FAILURE; return
    end if
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_value_float

function sfincs_get_value_double(this, name, dest) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  double precision,  intent(inout)  :: dest(:)
  integer                       :: status
  integer :: n
  select case (trim(name))
  case('water_surface_elevation')
    n = min(size(dest), size(this%z));  if (n>0 .and. associated(this%z))  dest(1:n) = this%z(1:n)
  case('water_depth')
    n = min(size(dest), size(this%h));  if (n>0 .and. associated(this%h))  dest(1:n) = this%h(1:n)
  !case('bed_elevation')
  !  if (allocated(zb)) then
  !    n = min(size(dest), size(zb)); if (n>0) dest(1:n) = real(zb(1:n),c_double)
  !  else
  !    status = BMI_FAILURE; return
  !  end if
  case('velocity_x')
    if (associated(this%un)) then
      n = min(size(dest), size(this%un)); if (n>0) dest(1:n) = this%un(1:n)
    end if
  case('velocity_y')
    if (associated(this%vn)) then
      n = min(size(dest), size(this%vn)); if (n>0) dest(1:n) = this%vn(1:n)
    end if
  case('rain_rate')
    if (associated(this%rain)) then
      n = min(size(dest), size(this%rain)); if (n>0) dest(1:n) = this%rain(1:n)
    else
      status = BMI_FAILURE; return
    end if
  case('bed_elevation')
    if (allocated(zb_store)) then
      n = min(size(dest), size(zb_store))
      if (n>0) dest(1:n) = zb_store(1:n)
      status = BMI_SUCCESS
    else
      status = BMI_FAILURE
  end if
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_value_double

!----------------------------------------------------------------------
!                       POINTER GETTERS (BMI 2.0)
!----------------------------------------------------------------------
function sfincs_get_value_ptr_int(this, name, dest_ptr) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  integer, pointer,  intent(inout) :: dest_ptr(:)
  integer                       :: status
  nullify(dest_ptr)
  status = BMI_FAILURE
end function sfincs_get_value_ptr_int

function sfincs_get_value_ptr_float(this, name, dest_ptr) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  real, pointer,     intent(inout) :: dest_ptr(:)
  integer                       :: status
  nullify(dest_ptr)
  status = BMI_FAILURE
end function sfincs_get_value_ptr_float

function sfincs_get_value_ptr_double(this, name, dest_ptr) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  double precision, pointer, intent(inout) :: dest_ptr(:)
  integer                       :: status
  select case (trim(name))
  case('water_surface_elevation'); if (associated(this%z))   dest_ptr => this%z
  case('water_depth');              if (associated(this%h))   dest_ptr => this%h
  case('velocity_x');               if (associated(this%un))  dest_ptr => this%un
  case('velocity_y');               if (associated(this%vn))  dest_ptr => this%vn
  case('rain_rate');                if (associated(this%rain))dest_ptr => this%rain
  case('bed_elevation');            nullify(dest_ptr); status=BMI_FAILURE; return  ! cannot point to zb (not TARGET)
  !case('bed_elevation');  n = merge(size(zb_store), -1, allocated(zb_store))
  case default
    nullify(dest_ptr); status = BMI_FAILURE; return
  end select
  if (.not. associated(dest_ptr)) then
    status = BMI_FAILURE
  else
    status = BMI_SUCCESS
  end if
end function sfincs_get_value_ptr_double

!----------------------------------------------------------------------
!                 GET VALUE AT INDICES (1-based inds)
!----------------------------------------------------------------------
function sfincs_get_value_at_indices_int(this, name, dest, inds) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  integer,           intent(inout) :: dest(:)
  integer,           intent(in) :: inds(:)
  integer                       :: status
  status = BMI_FAILURE
end function sfincs_get_value_at_indices_int

function sfincs_get_value_at_indices_float(this, name, dest, inds) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  real,              intent(inout) :: dest(:)
  integer,           intent(in) :: inds(:)
  integer                       :: status
  integer :: k, n, i0
  n = min(size(dest), size(inds))
  do k = 1, n
    i0 = inds(k)
    select case (trim(name))
    case('water_surface_elevation'); if (i0>=1 .and. i0<=size(this%z))  dest(k) = real(this%z(i0),  c_float)
    case('water_depth');              if (i0>=1 .and. i0<=size(this%h))  dest(k) = real(this%h(i0),  c_float)
    case('velocity_x');               if (i0>=1 .and. associated(this%un) .and. i0<=size(this%un)) dest(k) = real(this%un(i0), c_float)
    case('velocity_y');               if (i0>=1 .and. associated(this%vn) .and. i0<=size(this%vn)) dest(k) = real(this%vn(i0), c_float)
    case('rain_rate');                if (i0>=1 .and. i0<=size(this%rain)) dest(k) = real(this%rain(i0), c_float)
    case('bed_elevation');            if (allocated(zb) .and. i0>=1 .and. i0<=size(zb)) dest(k) = real(real(zb(i0),c_double), c_float)
    case default
      status = BMI_FAILURE; return
    end select
  end do
  status = BMI_SUCCESS
end function sfincs_get_value_at_indices_float

function sfincs_get_value_at_indices_double(this, name, dest, inds) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  double precision,   intent(inout) :: dest(:)
  integer,            intent(in) :: inds(:)
  integer                       :: status
  integer :: k, n, i0
  n = min(size(dest), size(inds))
  do k = 1, n
    i0 = inds(k)
    select case (trim(name))
    case('water_surface_elevation'); if (i0>=1 .and. i0<=size(this%z))  dest(k) = this%z(i0)
    case('water_depth');              if (i0>=1 .and. i0<=size(this%h))  dest(k) = this%h(i0)
    case('velocity_x');               if (i0>=1 .and. associated(this%un) .and. i0<=size(this%un)) dest(k) = this%un(i0)
    case('velocity_y');               if (i0>=1 .and. associated(this%vn) .and. i0<=size(this%vn)) dest(k) = this%vn(i0)
    case('rain_rate');                if (i0>=1 .and. i0<=size(this%rain)) dest(k) = this%rain(i0)
    case('bed_elevation');            if (allocated(zb) .and. i0>=1 .and. i0<=size(zb)) dest(k) = real(zb(i0),c_double)
    case default
      status = BMI_FAILURE; return
    end select
  end do
  status = BMI_SUCCESS
end function sfincs_get_value_at_indices_double

!======================================================================
!                        VALUE SETTERS
!======================================================================
function sfincs_set_value_int(this, name, src) result(status)
  class(sfincs_bmi), intent(inout) :: this
  character(len=*),  intent(in)    :: name
  integer,           intent(in)    :: src(:)
  integer                       :: status
  status = BMI_FAILURE
end function sfincs_set_value_int

function sfincs_set_value_float(this, name, src) result(status)
  class(sfincs_bmi), intent(inout) :: this
  character(len=*),  intent(in)    :: name
  real,              intent(in)    :: src(:)
  integer                       :: status
  integer :: n
  select case (trim(name))
  case('rain_rate')
    if (allocated(prcp)) then
      n = min(size(src), size(prcp))
      if (n > 0) prcp(1:n) = real(src(1:n))
      call refresh_from_sfincs()
      status = BMI_SUCCESS
      return
    end if
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_FAILURE
end function sfincs_set_value_float

function sfincs_set_value_double(this, name, src) result(status)
  class(sfincs_bmi), intent(inout) :: this
  character(len=*),  intent(in)    :: name
  double precision,  intent(in)    :: src(:)
  integer                       :: status
  integer :: n
  select case (trim(name))
  case('water_surface_elevation')
    if (allocated(zs)) then
      n = min(size(src), size(zs))
      if (n > 0) zs(1:n) = real(src(1:n))
      call refresh_from_sfincs()
      status = BMI_SUCCESS; return
    end if
  case('bed_elevation')
    if (allocated(zb)) then
      n = min(size(src), size(zb))
      if (n > 0) zb(1:n) = real(src(1:n))
      call refresh_from_sfincs()
      status = BMI_SUCCESS; return
    end if
  case('rain_rate')
    if (allocated(prcp)) then
      n = min(size(src), size(prcp))
      if (n > 0) prcp(1:n) = real(src(1:n))
      call refresh_from_sfincs()
      status = BMI_SUCCESS; return
    end if
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_FAILURE
end function sfincs_set_value_double

function sfincs_set_value_at_indices_int(this, name, inds, src) result(status)
  class(sfincs_bmi), intent(inout) :: this
  character(len=*),  intent(in)    :: name
  integer,           intent(in)    :: inds(:)
  integer,           intent(in)    :: src(:)
  integer                       :: status
  status = BMI_FAILURE
end function sfincs_set_value_at_indices_int

function sfincs_set_value_at_indices_float(this, name, inds, src) result(status)
  class(sfincs_bmi), intent(inout) :: this
  character(len=*),  intent(in)    :: name
  integer,           intent(in)    :: inds(:)
  real,              intent(in)    :: src(:)
  integer                       :: status
  integer :: k, n, i0
  select case (trim(name))
  case('rain_rate')
    if (allocated(prcp)) then
      n = min(size(src), size(inds))
      do k = 1, n
        i0 = inds(k)
        if (i0>=1 .and. i0<=size(prcp)) prcp(i0) = src(k)
      end do
      call refresh_from_sfincs()
      status = BMI_SUCCESS; return
    end if
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_FAILURE
end function sfincs_set_value_at_indices_float

function sfincs_set_value_at_indices_double(this, name, inds, src) result(status)
  class(sfincs_bmi), intent(inout) :: this
  character(len=*),  intent(in)    :: name
  integer,           intent(in)    :: inds(:)
  double precision,  intent(in)    :: src(:)
  integer                       :: status
  integer :: k, n, i0
  select case (trim(name))
  case('water_surface_elevation')
    if (allocated(zs)) then
      n = min(size(src), size(inds))
      do k = 1, n
        i0 = inds(k)
        if (i0>=1 .and. i0<=size(zs)) zs(i0) = real(src(k))
      end do
      call refresh_from_sfincs()
      status = BMI_SUCCESS; return
    end if
  case('bed_elevation')
    if (allocated(zb)) then
      n = min(size(src), size(inds))
      do k = 1, n
        i0 = inds(k)
        if (i0>=1 .and. i0<=size(zb)) zb(i0) = real(src(k))
      end do
      call refresh_from_sfincs()
      status = BMI_SUCCESS; return
    end if
  case('rain_rate')
    if (allocated(prcp)) then
      n = min(size(src), size(inds))
      do k = 1, n
        i0 = inds(k)
        if (i0>=1 .and. i0<=size(prcp)) prcp(i0) = real(src(k))
      end do
      call refresh_from_sfincs()
      status = BMI_SUCCESS; return
    end if
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_FAILURE
end function sfincs_set_value_at_indices_double

!======================================================================
!                         SMALL HELPERS
!======================================================================
subroutine assign_trim(lhs, rhs)
  character(len=*), intent(inout) :: lhs
  character(len=*), intent(in)    :: rhs
  lhs = ' '
  lhs(1:min(len(lhs), len_trim(rhs))) = rhs(1:min(len(lhs), len_trim(rhs)))
end subroutine assign_trim

integer function var_size(this, name) result(n)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  select case (trim(name))
  case('water_surface_elevation');  if (associated(this%z))  then; n = size(this%z);  else; n = 0; end if
  case('water_depth');              if (associated(this%h))  then; n = size(this%h);  else; n = 0; end if
  case('bed_elevation');            if (allocated(zb))       then; n = size(zb);      else; n = 0; end if
  case('velocity_x');               if (associated(this%un)) then; n = size(this%un); else; n = 0; end if
  case('velocity_y');               if (associated(this%vn)) then; n = size(this%vn); else; n = 0; end if
  case('rain_rate');                if (associated(this%rain)) then; n = size(this%rain); else; n = 0; end if
  case default; n = 0
  end select
end function

end module sfincs_bmi2

