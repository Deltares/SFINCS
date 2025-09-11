module sfincs_bmi2
  use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit
  use bmif_2_0,   only: bmi, BMI_SUCCESS, BMI_FAILURE, BMI_MAX_COMPONENT_NAME
  use sfincs_data, only: &
       ! grid sizes (we expose z-grid = np, uv-grid = npuv)
       np, npuv,                           &
       ! coordinates for z-grid
       z_xz, z_yz,                         &
       ! bathy and state
       zb, zs,                              &
       ! meteo inputs
       prcp

  implicit none
  private

  ! =========
  ! clocks
  ! =========
  real(c_double), save :: start_time_s   = 0.0d0
  real(c_double), save :: end_time_s     = 0.0d0
  real(c_double), save :: current_time_s = 0.0d0
  real(c_double), save :: dt_s           = 60.0d0

  ! =========
  ! BMI-visible velocity arrays (placeholder; 0 until solver provides)
  ! =========
  real(c_float),   allocatable, target, save :: velx_bmi(:)  ! on uv grid
  real(c_float),   allocatable, target, save :: vely_bmi(:)  ! on uv grid

  ! =========
  ! BMI wrapper type
  ! =========
  type, public, extends(bmi) :: sfincs_bmi
    ! time caches (not required, but convenient)
    real(c_double) :: current_time = 0.0d0
    real(c_double) :: start_time   = 0.0d0
    real(c_double) :: end_time     = 0.0d0
    real(c_double) :: time_step    = 60.0d0

    ! pointer-len strings for component metadata & I/O names
    character(len=:), pointer :: component_name => null()
    character(len=:), pointer :: time_units     => null()
    character(len=:), pointer :: input_names(:) => null()
    character(len=:), pointer :: output_names(:)=> null()

    ! Pointers to exposed arrays (we will point to targets)
    ! z-grid (np)
    real(c_float), pointer :: zs_ptr(:)   => null()
    real(c_float), pointer :: h_ptr(:)    => null()   ! computed depth if we allocate a shadow
    ! rainfall input (np)
    real(c_float), pointer :: prcp_ptr(:) => null()
    ! u/v-velocity (uv-grid)
    real(c_float), pointer :: u_ptr(:)    => null()
    real(c_float), pointer :: v_ptr(:)    => null()

    ! z-grid coordinates (np)
    real(c_float), pointer :: xz_ptr(:)   => null()
    real(c_float), pointer :: yz_ptr(:)   => null()
  contains
    ! lifecycle
    procedure :: initialize           => sfincs_initialize
    procedure :: update               => sfincs_update
    procedure :: update_until         => sfincs_update_until
    procedure :: finalize             => sfincs_finalize

    ! names & counts
    procedure :: get_component_name   => sfincs_get_component_name
    procedure :: get_input_item_count => sfincs_get_input_item_count
    procedure :: get_output_item_count=> sfincs_get_output_item_count
    procedure :: get_input_var_names  => sfincs_get_input_var_names
    procedure :: get_output_var_names => sfincs_get_output_var_names

    ! time
    procedure :: get_start_time       => sfincs_get_start_time
    procedure :: get_end_time         => sfincs_get_end_time
    procedure :: get_current_time     => sfincs_get_current_time
    procedure :: get_time_step        => sfincs_get_time_step
    procedure :: get_time_units       => sfincs_get_time_units

    ! var metadata
    procedure :: get_var_grid         => sfincs_get_var_grid
    procedure :: get_var_type         => sfincs_get_var_type
    procedure :: get_var_units        => sfincs_get_var_units
    procedure :: get_var_itemsize     => sfincs_get_var_itemsize
    procedure :: get_var_nbytes       => sfincs_get_var_nbytes
    procedure :: get_var_location     => sfincs_get_var_location

    ! grid info (we use 1D unstructured point sets)
    procedure :: get_grid_type        => sfincs_get_grid_type
    procedure :: get_grid_rank        => sfincs_get_grid_rank
    procedure :: get_grid_size        => sfincs_get_grid_size
    procedure :: get_grid_shape       => sfincs_get_grid_shape
    procedure :: get_grid_spacing     => sfincs_get_grid_spacing
    procedure :: get_grid_origin      => sfincs_get_grid_origin
    procedure :: get_grid_x           => sfincs_get_grid_x
    procedure :: get_grid_y           => sfincs_get_grid_y
    procedure :: get_grid_z           => sfincs_get_grid_z

    ! connectivity helpers (return empty)
    procedure :: get_grid_node_count  => sfincs_get_grid_node_count
    procedure :: get_grid_edge_count  => sfincs_get_grid_edge_count
    procedure :: get_grid_face_count  => sfincs_get_grid_face_count
    procedure :: get_grid_nodes_per_face => sfincs_get_grid_nodes_per_face
    procedure :: get_grid_edge_nodes  => sfincs_get_grid_edge_nodes
    procedure :: get_grid_face_edges  => sfincs_get_grid_face_edges
    procedure :: get_grid_face_nodes  => sfincs_get_grid_face_nodes

    ! value getters/setters
    procedure :: get_value_int        => sfincs_get_value_int
    procedure :: get_value_float      => sfincs_get_value_float
    procedure :: get_value_double     => sfincs_get_value_double

    procedure :: get_value_ptr_int    => sfincs_get_value_ptr_int
    procedure :: get_value_ptr_float  => sfincs_get_value_ptr_float
    procedure :: get_value_ptr_double => sfincs_get_value_ptr_double

    procedure :: get_value_at_indices_int    => sfincs_get_value_at_indices_int
    procedure :: get_value_at_indices_float  => sfincs_get_value_at_indices_float
    procedure :: get_value_at_indices_double => sfincs_get_value_at_indices_double

    procedure :: set_value_int        => sfincs_set_value_int
    procedure :: set_value_float      => sfincs_set_value_float
    procedure :: set_value_double     => sfincs_set_value_double

    procedure :: set_value_at_indices_int    => sfincs_set_value_at_indices_int
    procedure :: set_value_at_indices_float  => sfincs_set_value_at_indices_float
    procedure :: set_value_at_indices_double => sfincs_set_value_at_indices_double
  end type sfincs_bmi

  public :: sfincs_bmi
  public :: register_bmi   ! C binding name expected by your adapter

contains

!===========================
! C-visible factory symbol
!===========================
subroutine register_bmi(comp) bind(C, name="register_bmi")
  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
  type(c_ptr), intent(out) :: comp
  type(sfincs_bmi), pointer :: self
  allocate(self)
  comp = c_loc(self)
end subroutine register_bmi

!===========================
! initialize
!===========================
function sfincs_initialize(this, config_file) result(status)
  class(sfincs_bmi), intent(out) :: this
  character(len=*),  intent(in)  :: config_file
  integer                        :: status
  integer                        :: rc, fu
  logical                        :: exists
  character(len=1000)            :: line

  ! namelist
  double precision :: model_start_time, model_end_time, time_step_size
  integer          :: num_time_steps
  namelist /sfincs/ model_start_time, model_end_time, time_step_size, num_time_steps

  ! defaults
  model_start_time = 0.0d0
  model_end_time   = -1.0d0
  time_step_size   = 60.0d0
  num_time_steps   = -1

  if (len_trim(config_file) > 0) then
    inquire(file=trim(config_file), exist=exists)
    if (.not. exists) then
      write(stderr,'(3a)') 'Error: input file "', trim(config_file), '" does not exist.'
      status = BMI_FAILURE
      return
    end if
    open (action='read', file=trim(config_file), iostat=rc, newunit=fu)
    if (rc /= 0) then
      write(stderr,'(3a,i0)') 'Error: cannot open "', trim(config_file), '", iostat=', rc
      status = BMI_FAILURE
      return
    end if
    read (unit=fu, nml=sfincs, iostat=rc)
    if (rc /= 0) then
      backspace(fu)
      read(fu,'(A)') line
      write(stderr,'(A)') 'Invalid line in namelist: '//trim(line)
      write(stderr,'(3a)') 'Error: invalid namelist in "', trim(config_file), '".'
      close(fu)
      status = BMI_FAILURE
      return
    end if
    close(fu)
  end if

  ! derive end time if needed
  if (model_end_time < 0.0d0 .and. num_time_steps >= 0) then
    model_end_time = model_start_time + dble(num_time_steps) * time_step_size
  end if
  if (time_step_size <= 0.0d0) then
    write(stderr,'(a)') 'Error: time_step_size must be > 0.'
    status = BMI_FAILURE
    return
  end if
  if (model_end_time < model_start_time) then
    write(stderr,'(a)') 'Error: model_end_time must be >= model_start_time.'
    status = BMI_FAILURE
    return
  end if
  if (model_end_time == model_start_time) then
    ! default to 24 steps
    model_end_time = model_start_time + 24.0d0 * time_step_size
  end if

  ! apply clocks
  start_time_s   = model_start_time
  end_time_s     = model_end_time
  dt_s           = time_step_size
  current_time_s = start_time_s

  this%start_time   = start_time_s
  this%end_time     = end_time_s
  this%time_step    = dt_s
  this%current_time = current_time_s

  ! ensure velocity placeholder arrays are allocated (uv grid)
  if (npuv > 0) then
    if (.not. allocated(velx_bmi)) allocate(velx_bmi(npuv)); velx_bmi = 0.0_c_float
    if (.not. allocated(vely_bmi)) allocate(vely_bmi(npuv)); vely_bmi = 0.0_c_float
  end if

  ! pointer-associate BMI pointers to SFINCS data
  if (np > 0) then
    ! state (z-grid)
    this%zs_ptr  => zs              ! zs is real*4, OK for c_float
    ! water depth we compute on the fly in getters; do not allocate here

    ! inputs
    this%prcp_ptr => prcp

    ! coords
    this%xz_ptr => z_xz
    this%yz_ptr => z_yz
  end if

  if (npuv > 0) then
    this%u_ptr => velx_bmi
    this%v_ptr => vely_bmi
  end if

  ! strings and var lists
  if (.not. associated(this%component_name)) then
    allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%component_name)
  end if
  this%component_name = 'SFINCS-BMI'

  if (.not. associated(this%time_units)) then
    allocate(character(len=16) :: this%time_units)
  end if
  this%time_units = 's'

  ! inputs
  if (.not. associated(this%input_names)) then
    allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%input_names(1))
  end if
  this%input_names(1) = 'rain_rate'

  ! outputs required by your NGen side:
  if (.not. associated(this%output_names)) then
    allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%output_names(4))
  end if
  this%output_names = [ character(len=BMI_MAX_COMPONENT_NAME) :: &
       'water_surface_elevation', 'water_depth', 'velocity_x', 'velocity_y' ]

  status = BMI_SUCCESS
end function sfincs_initialize

!===========================
! update
!===========================
function sfincs_update(this) result(status)
  class(sfincs_bmi), intent(inout) :: this
  integer                          :: status
  ! Here you would call the real SFINCS stepper. For now, just advance time.
  current_time_s = current_time_s + dt_s
  this%current_time = current_time_s
  status = BMI_SUCCESS
end function sfincs_update

!===========================
! update_until
!===========================
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

!===========================
! finalize
!===========================
function sfincs_finalize(this) result(status)
  class(sfincs_bmi), intent(inout) :: this
  integer                          :: status

  if (associated(this%component_name)) deallocate(this%component_name)
  if (associated(this%time_units))     deallocate(this%time_units)
  if (associated(this%input_names))    deallocate(this%input_names)
  if (associated(this%output_names))   deallocate(this%output_names)

  if (allocated(velx_bmi)) deallocate(velx_bmi)
  if (allocated(vely_bmi)) deallocate(vely_bmi)

  nullify(this%zs_ptr, this%h_ptr, this%prcp_ptr, this%u_ptr, this%v_ptr, this%xz_ptr, this%yz_ptr)

  status = BMI_SUCCESS
end function sfincs_finalize

!===========================
! names & counts
!===========================
function sfincs_get_component_name(this, name) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*), intent(out) :: name
  integer                       :: status
  call assign_trim(name, 'SFINCS-BMI')
  status = BMI_SUCCESS
end function sfincs_get_component_name

function sfincs_get_input_item_count(this, count) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(out):: count
  integer                       :: status
  count = 1
  status = BMI_SUCCESS
end function sfincs_get_input_item_count

function sfincs_get_output_item_count(this, count) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(out):: count
  integer                       :: status
  count = 4
  status = BMI_SUCCESS
end function sfincs_get_output_item_count

function sfincs_get_input_var_names(this, names) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*), intent(out) :: names(:)
  integer                       :: status
  integer :: n
  n = min(size(names), 1)
  if (n >= 1) call assign_trim(names(1), 'rain_rate')
  status = BMI_SUCCESS
end function sfincs_get_input_var_names

function sfincs_get_output_var_names(this, names) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*), intent(out) :: names(:)
  integer                       :: status
  integer :: n
  n = min(size(names), 4)
  if (n >= 1) call assign_trim(names(1), 'water_surface_elevation')
  if (n >= 2) call assign_trim(names(2), 'water_depth')
  if (n >= 3) call assign_trim(names(3), 'velocity_x')
  if (n >= 4) call assign_trim(names(4), 'velocity_y')
  status = BMI_SUCCESS
end function sfincs_get_output_var_names

!===========================
! time
!===========================
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
  call assign_trim(units, 's')
  status = BMI_SUCCESS
end function sfincs_get_time_units

!===========================
! var metadata
!===========================
function sfincs_get_var_grid(this, name, grid) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  integer,           intent(out):: grid
  integer                       :: status
  select case (trim(name))
  case('water_surface_elevation','water_depth','rain_rate')
    grid = 1      ! z-grid
  case('velocity_x','velocity_y')
    grid = 2      ! uv-grid
  case default
    grid = -1; status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_var_grid

function sfincs_get_var_type(this, name, type) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  character(len=*),  intent(out):: type
  integer                       :: status
  ! SFINCS internal uses real*4; expose as 'float'
  call assign_trim(type, 'float')
  status = BMI_SUCCESS
end function sfincs_get_var_type

function sfincs_get_var_units(this, name, units) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  character(len=*),  intent(out):: units
  integer                       :: status
  select case (trim(name))
  case('water_surface_elevation','water_depth'); call assign_trim(units, 'm')
  case('rain_rate');                             call assign_trim(units, 'm s-1')
  case('velocity_x','velocity_y');               call assign_trim(units, 'm s-1')
  case default;                                  call assign_trim(units, '1'); status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_var_units

function sfincs_get_var_itemsize(this, name, size) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  integer,           intent(out):: size
  integer                       :: status
  size = 4   ! float
  status = BMI_SUCCESS
end function sfincs_get_var_itemsize

function sfincs_get_var_nbytes(this, name, nbytes) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  integer,           intent(out):: nbytes
  integer                       :: status
  integer :: n
  n = var_size(this, name)
  if (n < 0) then
    nbytes = 0; status = BMI_FAILURE; return
  end if
  nbytes = 4*n
  status = BMI_SUCCESS
end function sfincs_get_var_nbytes

function sfincs_get_var_location(this, name, location) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  character(len=*),  intent(out):: location
  integer                       :: status
  select case (trim(name))
  case('water_surface_elevation','water_depth','rain_rate')
    call assign_trim(location, 'node') ! z-grid nodes
  case('velocity_x','velocity_y')
    call assign_trim(location, 'edge') ! uv-grid edges
  case default
    call assign_trim(location, 'unknown'); status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_var_location

!===========================
! grid info (1D unstructured)
!===========================
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
  select case(grid)
  case(1); size = np
  case(2); size = npuv
  case default; size = 0; status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_grid_size

function sfincs_get_grid_shape(this, grid, shape) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: shape(:)
  integer                       :: status
  if (size(shape) < 1) then
    status = BMI_FAILURE; return
  end if
  select case(grid)
  case(1); shape(1) = np
  case(2); shape(1) = npuv
  case default; status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
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
  select case(grid)
  case(1)
    n = min(size(x), np)
    if (n > 0) x(1:n) = real(this%xz_ptr(1:n), kind=c_double)
  case(2)
    ! UV-grid coordinates not available; return zeros
    if (size(x) > 0) x = 0.0d0
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_grid_x

function sfincs_get_grid_y(this, grid, y) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  double precision,  intent(out):: y(:)
  integer                       :: status
  integer :: n
  select case(grid)
  case(1)
    n = min(size(y), np)
    if (n > 0) y(1:n) = real(this%yz_ptr(1:n), kind=c_double)
  case(2)
    if (size(y) > 0) y = 0.0d0
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_grid_y

function sfincs_get_grid_z(this, grid, z) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  double precision,  intent(out):: z(:)
  integer                       :: status
  if (size(z) >= 1) z = 0.0d0
  status = BMI_SUCCESS
end function sfincs_get_grid_z

function sfincs_get_grid_node_count(this, grid, count) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: count
  integer                       :: status
  select case(grid)
  case(1); count = np
  case(2); count = npuv
  case default; count = 0; status = BMI_FAILURE; return
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

!===========================
! getters
!===========================
function sfincs_get_value_int(this, name, dest) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  integer,           intent(inout) :: dest(:)
  integer                       :: status
  ! no int fields exposed
  dest = 0
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
    n = min(size(dest), np); if (n>0) dest(1:n) = real(this%zs_ptr(1:n), kind=kind(dest))
  case('water_depth')
    n = min(size(dest), np)
    if (n > 0) dest(1:n) = max(0.0, real(this%zs_ptr(1:n), kind=kind(dest)) - real(zb(1:n), kind=kind(dest)))
  case('velocity_x')
    n = min(size(dest), npuv); if (n>0) dest(1:n) = real(this%u_ptr(1:n), kind=kind(dest))
  case('velocity_y')
    n = min(size(dest), npuv); if (n>0) dest(1:n) = real(this%v_ptr(1:n), kind=kind(dest))
  case('rain_rate')
    n = min(size(dest), np); if (n>0) dest(1:n) = real(this%prcp_ptr(1:n), kind=kind(dest))
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_value_float

function sfincs_get_value_double(this, name, dest) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  double precision,  intent(inout) :: dest(:)
  integer                       :: status
  integer :: n
  select case (trim(name))
  case('water_surface_elevation')
    n = min(size(dest), np); if (n>0) dest(1:n) = real(this%zs_ptr(1:n), kind=c_double)
  case('water_depth')
    n = min(size(dest), np)
    if (n > 0) dest(1:n) = max(0.0d0, real(this%zs_ptr(1:n), kind=c_double) - real(zb(1:n), kind=c_double))
  case('velocity_x')
    n = min(size(dest), npuv); if (n>0) dest(1:n) = real(this%u_ptr(1:n), kind=c_double)
  case('velocity_y')
    n = min(size(dest), npuv); if (n>0) dest(1:n) = real(this%v_ptr(1:n), kind=c_double)
  case('rain_rate')
    n = min(size(dest), np); if (n>0) dest(1:n) = real(this%prcp_ptr(1:n), kind=c_double)
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_value_double

!===========================
! pointer getters
!===========================
function sfincs_get_value_ptr_int(this, name, dest_ptr) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  integer, pointer,  intent(out):: dest_ptr(:)
  integer                       :: status
  nullify(dest_ptr)
  status = BMI_FAILURE
end function sfincs_get_value_ptr_int

function sfincs_get_value_ptr_float(this, name, dest_ptr) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  real, pointer,     intent(out):: dest_ptr(:)
  integer                       :: status
  ! Default real kind may not match c_float; be conservative and fail.
  nullify(dest_ptr)
  status = BMI_FAILURE
end function sfincs_get_value_ptr_float

function sfincs_get_value_ptr_double(this, name, dest_ptr) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  double precision, pointer, intent(out) :: dest_ptr(:)
  integer                       :: status
  ! Our exposed arrays are float; do not provide double pointers.
  nullify(dest_ptr)
  status = BMI_FAILURE
end function sfincs_get_value_ptr_double

!===========================
! get_value_at_indices
!===========================
function sfincs_get_value_at_indices_int(this, name, dest, inds) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  integer,           intent(inout) :: dest(:)
  integer,           intent(in) :: inds(:)
  integer                       :: status
  dest = 0
  status = BMI_FAILURE
end function sfincs_get_value_at_indices_int

function sfincs_get_value_at_indices_float(this, name, dest, inds) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  real,              intent(inout) :: dest(:)
  integer,           intent(in) :: inds(:)
  integer                       :: status
  integer :: n, k, i0
  n = min(size(dest), size(inds))
  do k = 1, n
    i0 = inds(k)
    select case (trim(name))
    case('water_surface_elevation'); if (i0>=1 .and. i0<=np)   dest(k) = real(this%zs_ptr(i0), kind=kind(dest))
    case('water_depth');              if (i0>=1 .and. i0<=np)   dest(k) = max(0.0, real(this%zs_ptr(i0), kind=kind(dest)) - real(zb(i0), kind=kind(dest)))
    case('velocity_x');               if (i0>=1 .and. i0<=npuv) dest(k) = real(this%u_ptr(i0), kind=kind(dest))
    case('velocity_y');               if (i0>=1 .and. i0<=npuv) dest(k) = real(this%v_ptr(i0), kind=kind(dest))
    case('rain_rate');                if (i0>=1 .and. i0<=np)   dest(k) = real(this%prcp_ptr(i0), kind=kind(dest))
    case default
      status = BMI_FAILURE; return
    end select
  end do
  status = BMI_SUCCESS
end function sfincs_get_value_at_indices_float

function sfincs_get_value_at_indices_double(this, name, dest, inds) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  double precision,  intent(inout) :: dest(:)
  integer,           intent(in) :: inds(:)
  integer                       :: status
  integer :: n, k, i0
  n = min(size(dest), size(inds))
  do k = 1, n
    i0 = inds(k)
    select case (trim(name))
    case('water_surface_elevation'); if (i0>=1 .and. i0<=np)   dest(k) = real(this%zs_ptr(i0), kind=c_double)
    case('water_depth');              if (i0>=1 .and. i0<=np)   dest(k) = max(0.0d0, real(this%zs_ptr(i0), kind=c_double) - real(zb(i0), kind=c_double))
    case('velocity_x');               if (i0>=1 .and. i0<=npuv) dest(k) = real(this%u_ptr(i0), kind=c_double)
    case('velocity_y');               if (i0>=1 .and. i0<=npuv) dest(k) = real(this%v_ptr(i0), kind=c_double)
    case('rain_rate');                if (i0>=1 .and. i0<=np)   dest(k) = real(this%prcp_ptr(i0), kind=c_double)
    case default
      status = BMI_FAILURE; return
    end select
  end do
  status = BMI_SUCCESS
end function sfincs_get_value_at_indices_double

!===========================
! setters
!===========================
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
  case('water_surface_elevation')
    n = min(size(src), np);   if (n>0) this%zs_ptr(1:n)   = real(src(1:n), kind=c_float)
  case('water_depth')
    ! typically read-only; ignore or compute back to zs by adding zb
    n = min(size(src), np);   if (n>0) this%zs_ptr(1:n)   = max(0.0_c_float, real(src(1:n),kind=c_float) + real(zb(1:n),kind=c_float))
  case('velocity_x')
    n = min(size(src), npuv); if (n>0) this%u_ptr(1:n)    = real(src(1:n), kind=c_float)
  case('velocity_y')
    n = min(size(src), npuv); if (n>0) this%v_ptr(1:n)    = real(src(1:n), kind=c_float)
  case('rain_rate')
    n = min(size(src), np);   if (n>0) this%prcp_ptr(1:n) = real(src(1:n), kind=c_float)
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_set_value_float

function sfincs_set_value_double(this, name, src) result(status)
  class(sfincs_bmi), intent(inout) :: this
  character(len=*),  intent(in)    :: name
  double precision,  intent(in)    :: src(:)
  integer                       :: status
  integer :: n
  select case (trim(name))
  case('water_surface_elevation')
    n = min(size(src), np);   if (n>0) this%zs_ptr(1:n)   = real(src(1:n), kind=c_float)
  case('water_depth')
    n = min(size(src), np);   if (n>0) this%zs_ptr(1:n)   = max(0.0_c_float, real(src(1:n),kind=c_float) + real(zb(1:n),kind=c_float))
  case('velocity_x')
    n = min(size(src), npuv); if (n>0) this%u_ptr(1:n)    = real(src(1:n), kind=c_float)
  case('velocity_y')
    n = min(size(src), npuv); if (n>0) this%v_ptr(1:n)    = real(src(1:n), kind=c_float)
  case('rain_rate')
    n = min(size(src), np);   if (n>0) this%prcp_ptr(1:n) = real(src(1:n), kind=c_float)
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
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
  n = min(size(src), size(inds))
  select case (trim(name))
  case('water_surface_elevation')
    do k = 1, n
      i0 = inds(k); if (i0>=1 .and. i0<=np) this%zs_ptr(i0) = real(src(k), kind=c_float)
    end do
  case('water_depth')
    do k = 1, n
      i0 = inds(k); if (i0>=1 .and. i0<=np) this%zs_ptr(i0) = max(0.0_c_float, real(src(k),kind=c_float) + real(zb(i0),kind=c_float))
    end do
  case('velocity_x')
    do k = 1, n
      i0 = inds(k); if (i0>=1 .and. i0<=npuv) this%u_ptr(i0) = real(src(k), kind=c_float)
    end do
  case('velocity_y')
    do k = 1, n
      i0 = inds(k); if (i0>=1 .and. i0<=npuv) this%v_ptr(i0) = real(src(k), kind=c_float)
    end do
  case('rain_rate')
    do k = 1, n
      i0 = inds(k); if (i0>=1 .and. i0<=np) this%prcp_ptr(i0) = real(src(k), kind=c_float)
    end do
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_set_value_at_indices_float

function sfincs_set_value_at_indices_double(this, name, inds, src) result(status)
  class(sfincs_bmi), intent(inout) :: this
  character(len=*),  intent(in)    :: name
  integer,           intent(in)    :: inds(:)
  double precision,  intent(in)    :: src(:)
  integer                       :: status
  integer :: k, n, i0
  n = min(size(src), size(inds))
  select case (trim(name))
  case('water_surface_elevation')
    do k = 1, n
      i0 = inds(k); if (i0>=1 .and. i0<=np) this%zs_ptr(i0) = real(src(k), kind=c_float)
    end do
  case('water_depth')
    do k = 1, n
      i0 = inds(k); if (i0>=1 .and. i0<=np) this%zs_ptr(i0) = max(0.0_c_float, real(src(k),kind=c_float) + real(zb(i0),kind=c_float))
    end do
  case('velocity_x')
    do k = 1, n
      i0 = inds(k); if (i0>=1 .and. i0<=npuv) this%u_ptr(i0) = real(src(k), kind=c_float)
    end do
  case('velocity_y')
    do k = 1, n
      i0 = inds(k); if (i0>=1 .and. i0<=npuv) this%v_ptr(i0) = real(src(k), kind=c_float)
    end do
  case('rain_rate')
    do k = 1, n
      i0 = inds(k); if (i0>=1 .and. i0<=np) this%prcp_ptr(i0) = real(src(k), kind=c_float)
    end do
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_set_value_at_indices_double

!===========================
! helpers
!===========================
subroutine assign_trim(lhs, rhs)
  character(len=*), intent(out) :: lhs
  character(len=*), intent(in)  :: rhs
  integer :: n
  lhs = ' '
  n = min(len(lhs), len_trim(rhs))
  if (n > 0) lhs(1:n) = rhs(1:n)
end subroutine assign_trim

integer function var_size(this, name) result(n)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  select case (trim(name))
  case('water_surface_elevation'); n = np
  case('water_depth');              n = np
  case('velocity_x');               n = npuv
  case('velocity_y');               n = npuv
  case('rain_rate');                n = np
  case default;                     n = -1
  end select
end function var_size

end module sfincs_bmi2

