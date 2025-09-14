module sfincs_bmi2
  use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double, c_sizeof, c_ptr, c_null_ptr
  use bmif_2_0,   only: bmi, BMI_SUCCESS, BMI_FAILURE, BMI_MAX_COMPONENT_NAME
  use sfincs_data
  implicit none

  private
  public :: sfincs_bmi

  !========================
  ! Module-scope storage
  !========================
  ! NOTE: These are simple stand-ins that the BMI pointer components will associate to
  !       (when TARGET is available). For now we keep internal TARGET arrays and copy.
  !
  ! z-grid (nodes/cell centers)
  integer,        save :: nz = 0
  real(c_double), save, target, allocatable :: z_store(:)   ! water surface elevation [m]
  real(c_double), save, target, allocatable :: h_store(:)   ! water depth [m]
  real(c_double), save, target, allocatable :: xz_store(:), yz_store(:)

  ! u-grid (staggered x-velocity points)
  integer,        save :: nu = 0
  real(c_double), save, target, allocatable :: un_store(:)  ! u-velocity [m/s]
  real(c_double), save, target, allocatable :: xu_store(:), yu_store(:)

  ! v-grid (staggered y-velocity points)
  integer,        save :: nv = 0
  real(c_double), save, target, allocatable :: vn_store(:)  ! v-velocity [m/s]
  real(c_double), save, target, allocatable :: xv_store(:), yv_store(:)

  ! rainfall input over z-grid (nodes/cell centers), units m s-1
  real(c_double), save, target, allocatable :: rain_store(:)

  !------------------------
  ! Time bookkeeping
  !------------------------
  real(c_double), save :: start_time_s   = 0.0d0
  real(c_double), save :: end_time_s     = 0.0d0
  real(c_double), save :: current_time_s = 0.0d0
  real(c_double), save :: dt_s           = 60.0d0  ! default 60 s; override in initialize

  !========================
  ! BMI wrapper type
  !========================
  type, extends(bmi) :: sfincs_bmi

    real(c_double)                    :: current_time = 0.0d0
    real(c_double)                    :: start_time   = 0.0d0
    real(c_double)                    :: end_time     = 0.0d0
    real(c_double)                    :: time_step    = 1.0d0

    ! pointer state buffers (never TARGET on components)
    integer(c_int),    pointer        :: state_i(:)  => null()
    real(c_float),     pointer        :: state_r(:)  => null()
    real(c_double),    pointer        :: state_d(:)  => null()

    real(c_double),    pointer        :: rain(:)     => null()

    ! Use POINTER (not ALLOCATABLE) so we can test with ASSOCIATED()
    character(len=:), pointer     :: component_name => null()
    character(len=:), pointer     :: time_units     => null()
    character(len=:), pointer     :: input_names(:) => null()
    character(len=:), pointer     :: output_names(:)=> null()

    ! state pointers (point to the module-scope TARGET arrays above)
    real(c_double), pointer :: z(:)  => null()
    real(c_double), pointer :: h(:)  => null()
    real(c_double), pointer :: un(:) => null()
    real(c_double), pointer :: vn(:) => null()

    ! coord pointers
    real(c_double), pointer :: xz(:) => null()
    real(c_double), pointer :: yz(:) => null()
    real(c_double), pointer :: xu(:) => null()
    real(c_double), pointer :: yu(:) => null()
    real(c_double), pointer :: xv(:) => null()
    real(c_double), pointer :: yv(:) => null()

  contains
    ! Lifecycle
    procedure :: initialize                 => sfincs_initialize
    procedure :: update                     => sfincs_update
    procedure :: update_until               => sfincs_update_until
    procedure :: finalize                   => sfincs_finalize

    ! Component / I/O names (BMI 2.0 here expects POINTER dummies in the base)
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

    ! Grid info (unstructured point sets)
    procedure :: get_grid_type              => sfincs_get_grid_type
    procedure :: get_grid_rank              => sfincs_get_grid_rank
    procedure :: get_grid_size              => sfincs_get_grid_size
    procedure :: get_grid_shape             => sfincs_get_grid_shape
    procedure :: get_grid_spacing           => sfincs_get_grid_spacing
    procedure :: get_grid_origin            => sfincs_get_grid_origin
    procedure :: get_grid_x                 => sfincs_get_grid_x
    procedure :: get_grid_y                 => sfincs_get_grid_y
    procedure :: get_grid_z                 => sfincs_get_grid_z

    ! Optional connectivity helpers
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

    ! BMI 2.0 in your build expects FORTRAN POINTER arrays here (not C_PTR)
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

!======================================================================
!                       LIFECYCLE
!======================================================================
function sfincs_initialize(this, config_file) result(status)
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit
  class(sfincs_bmi), intent(out)   :: this
  character(len=*),  intent(in)    :: config_file
  integer                         :: status
  integer                         :: i
  ! Local clock helpers
  double precision :: dt_local, start_local, end_local
  integer          :: n_steps

  ! 1) Read config if provided; otherwise set defaults
  if (len_trim(config_file) > 0) then
    call read_init_config(this, config_file, status)
    if (status /= BMI_SUCCESS) return
  else
    start_time_s   = 0.0d0
    dt_s           = 60.0d0
    end_time_s     = 24.0d0 * 3600.0d0
    current_time_s = start_time_s
  end if

  ! 2) Sanity/derivations similar to SCHISM initialize
  dt_local    = dt_s
  start_local = start_time_s
  end_local   = end_time_s

  if (dt_local <= 0.0d0) then
    write(stderr,'(a)') 'Error: dt_seconds (time_step_size) must be > 0.'
    status = BMI_FAILURE
    return
  end if

  ! If end time wasn’t set, default to 24 steps
  if (end_local <= start_local) then
    n_steps   = 24
    end_local = start_local + dble(n_steps) * dt_local
  end if

  if (end_local < start_local) then
    write(stderr,'(a)') 'Error: end_time < start_time after derivation.'
    status = BMI_FAILURE
    return
  end if

  ! Apply back to module clock
  start_time_s   = start_local
  end_time_s     = end_local
  dt_s           = dt_local
  current_time_s = start_time_s

  ! 3) Allocate simple demo domain (replace with real SFINCS wiring later)
  if (nz <= 0) then
    nz = 10_c_int
    allocate(z_store(nz), h_store(nz), xz_store(nz), yz_store(nz))
    do i = 1, nz
      xz_store(i) = real(i-1, c_double)
      yz_store(i) = 0.0d0
    end do
    z_store = 0.0d0
    h_store = 0.0d0
  end if

  if (nu <= 0) then
    nu = 10_c_int
    allocate(un_store(nu), xu_store(nu), yu_store(nu))
    do i = 1, nu
      xu_store(i) = real(i-1, c_double) + 0.5d0
      yu_store(i) = 0.0d0
    end do
    un_store = 0.0d0
  end if

  if (nv <= 0) then
    nv = 10_c_int
    allocate(vn_store(nv), xv_store(nv), yv_store(nv))
    do i = 1, nv
      xv_store(i) = real(i-1, c_double)
      yv_store(i) = 0.5d0
    end do
    vn_store = 0.0d0
  end if

  ! Rainfall field on z-grid
  if (.not. allocated(rain_store)) then
    allocate(rain_store(nz))
    rain_store = 0.0d0
  end if
  this%rain => rain_store

  ! 4) Pointer-associate state to module TARGET arrays
  this%z  => z_store
  this%h  => h_store
  this%un => un_store
  this%vn => vn_store
  this%xz => xz_store; this%yz => yz_store
  this%xu => xu_store; this%yu => yu_store
  this%xv => xv_store; this%yv => yv_store

  ! 5) Allocate/assign scalar/pointer-len strings
  if (.not. associated(this%component_name)) then
    allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%component_name)
    this%component_name = 'SFINCS-BMI'
  end if

  if (.not. associated(this%time_units)) then
    allocate(character(len=16) :: this%time_units)
    this%time_units = 's'
  end if

  ! 6) Inputs/outputs advertised to NGen
  if (.not. associated(this%input_names)) then
    allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%input_names(1))
    this%input_names(1) = 'rain_rate'
  end if

  if (.not. associated(this%output_names)) then
    allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%output_names(4))
    this%output_names = [ character(len=BMI_MAX_COMPONENT_NAME) :: &
         'water_surface_elevation', 'water_depth', 'velocity_x', 'velocity_y' ]
  end if

  status = BMI_SUCCESS
end function sfincs_initialize

function sfincs_update(this) result(status)
  class(sfincs_bmi), intent(inout) :: this
  integer                          :: status
  integer                          :: i

  ! Placeholder: shallow update (advect nothing, just a toy trend)
  do i = 1, size(this%z)
    this%z(i) = this%z(i) + 0.0d0
  end do

  do i = 1, size(this%h)
    this%h(i) = max(0.0d0, this%h(i) + rain_store(i) * dt_s)
  end do

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
  ! BMI core will provide/allocate these; just alias our internal list
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
  if (associated(this%time_units)) then
    call assign_trim(units, this%time_units)
    status = BMI_SUCCESS
  else
    status = BMI_FAILURE
  end if
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
  case('water_surface_elevation','water_depth','rain_rate'); grid = 1   ! z-grid
  case('velocity_x');                                        grid = 2   ! u-grid
  case('velocity_y');                                        grid = 3   ! v-grid
  case default;                                              grid = -1; status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
end function sfincs_get_var_grid

function sfincs_get_var_type(this, name, type) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  character(len=*),  intent(out):: type
  integer                       :: status
  call assign_trim(type, 'double')   ! expose all as double for simplicity
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
  size = 8    ! bytes per double
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
  nbytes = 8*n
  status = BMI_SUCCESS
end function sfincs_get_var_nbytes

function sfincs_get_var_location(this, name, location) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  character(len=*),  intent(out):: location
  integer                       :: status
  select case (trim(name))
  case('water_surface_elevation','water_depth','rain_rate'); call assign_trim(location,'node')
  case('velocity_x','velocity_y');                           call assign_trim(location,'edge')
  case default;                                              call assign_trim(location,'unknown'); status = BMI_FAILURE; return
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
  ! Unstructured: spacing undefined; return zeros
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
  case (1); n = min(size(x), nz); if (n>0) x(1:n) = this%xz(1:n)
  case (2); n = min(size(x), nu); if (n>0) x(1:n) = this%xu(1:n)
  case (3); n = min(size(x), nv); if (n>0) x(1:n) = this%xv(1:n)
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
  case (1); n = min(size(y), nz); if (n>0) y(1:n) = this%yz(1:n)
  case (2); n = min(size(y), nu); if (n>0) y(1:n) = this%yu(1:n)
  case (3); n = min(size(y), nv); if (n>0) y(1:n) = this%yv(1:n)
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

! Connectivity: provide flattened 1-D outputs as requested (names + rank 1)
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
  integer,           intent(out):: edge_nodes(:)   ! flattened pairs
  integer                       :: status
  if (size(edge_nodes) > 0) edge_nodes = 0
  status = BMI_SUCCESS
end function sfincs_get_grid_edge_nodes

function sfincs_get_grid_face_edges(this, grid, face_edges) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: face_edges(:)   ! flattened pairs
  integer                       :: status
  if (size(face_edges) > 0) face_edges = 0
  status = BMI_SUCCESS
end function sfincs_get_grid_face_edges

function sfincs_get_grid_face_nodes(this, grid, face_nodes) result(status)
  class(sfincs_bmi), intent(in) :: this
  integer,           intent(in) :: grid
  integer,           intent(out):: face_nodes(:)   ! flattened list
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
  integer :: n
  select case (trim(name))
  case default
    status = BMI_FAILURE; return
  end select
  n = 0
  status = BMI_SUCCESS
end function sfincs_get_value_int

function sfincs_get_value_float(this, name, dest) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  real,              intent(inout) :: dest(:)
  integer                       :: status
  integer :: n
  select case (trim(name))
  case('water_surface_elevation'); n = min(size(dest), size(this%z));  if (n>0) dest(1:n) = real(this%z(1:n),  c_float)
  case('water_depth');              n = min(size(dest), size(this%h));  if (n>0) dest(1:n) = real(this%h(1:n),  c_float)
  case('velocity_x');               n = min(size(dest), size(this%un)); if (n>0) dest(1:n) = real(this%un(1:n), c_float)
  case('velocity_y');               n = min(size(dest), size(this%vn)); if (n>0) dest(1:n) = real(this%vn(1:n), c_float)
  case('rain_rate');                n = min(size(dest), size(this%rain)); if (n>0) dest(1:n) = real(this%rain(1:n), c_float)
  case default; status = BMI_FAILURE; return
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
  case('water_surface_elevation'); n = min(size(dest), size(this%z));  if (n>0) dest(1:n) = this%z(1:n)
  case('water_depth');              n = min(size(dest), size(this%h));  if (n>0) dest(1:n) = this%h(1:n)
  case('velocity_x');               n = min(size(dest), size(this%un)); if (n>0) dest(1:n) = this%un(1:n)
  case('velocity_y');               n = min(size(dest), size(this%vn)); if (n>0) dest(1:n) = this%vn(1:n)
  case('rain_rate');                n = min(size(dest), size(this%rain)); if (n>0) dest(1:n) = this%rain(1:n)
  case default; status = BMI_FAILURE; return
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
  ! No default-kind integer state; return NULL pointer and failure
  nullify(dest_ptr)
  status = BMI_FAILURE
end function sfincs_get_value_ptr_int

function sfincs_get_value_ptr_float(this, name, dest_ptr) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  real, pointer,     intent(inout) :: dest_ptr(:)
  integer                       :: status
  ! No default-kind real state; return NULL pointer and failure
  nullify(dest_ptr)
  status = BMI_FAILURE
end function sfincs_get_value_ptr_float

function sfincs_get_value_ptr_double(this, name, dest_ptr) result(status)
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in) :: name
  double precision, pointer, intent(inout) :: dest_ptr(:)
  integer                       :: status
  select case (trim(name))
  case('water_surface_elevation'); dest_ptr => this%z
  case('water_depth');              dest_ptr => this%h
  case('velocity_x');               dest_ptr => this%un
  case('velocity_y');               dest_ptr => this%vn
  case('rain_rate');                dest_ptr => this%rain
  case default
    nullify(dest_ptr); status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
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
  integer :: k, n
  n = min(size(dest), size(inds))
  do k = 1, n
    dest(k) = 0
  end do
  status = BMI_SUCCESS
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
    case('velocity_x');               if (i0>=1 .and. i0<=size(this%un)) dest(k) = real(this%un(i0), c_float)
    case('velocity_y');               if (i0>=1 .and. i0<=size(this%vn)) dest(k) = real(this%vn(i0), c_float)
    case('rain_rate');                if (i0>=1 .and. i0<=size(this%rain)) dest(k) = real(this%rain(i0), c_float)
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
    case('velocity_x');               if (i0>=1 .and. i0<=size(this%un)) dest(k) = this%un(i0)
    case('velocity_y');               if (i0>=1 .and. i0<=size(this%vn)) dest(k) = this%vn(i0)
    case('rain_rate');                if (i0>=1 .and. i0<=size(this%rain)) dest(k) = this%rain(i0)
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
  case('velocity_x'); n = min(size(src), size(this%un));    if (n>0) this%un(1:n)         = real(src(1:n), c_double)
  case('velocity_y'); n = min(size(src), size(this%vn));    if (n>0) this%vn(1:n)         = real(src(1:n), c_double)
  case('rain_rate');  n = min(size(src), size(rain_store)); if (n>0) rain_store(1:n)      = real(src(1:n), c_double)
  case default; status = BMI_FAILURE; return
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
  case('water_surface_elevation');  n = min(size(src), size(this%z));    if (n>0) this%z(1:n)     = src(1:n)
  case('water_depth');              n = min(size(src), size(this%h));    if (n>0) this%h(1:n)     = src(1:n)
  case('velocity_x');               n = min(size(src), size(this%un));   if (n>0) this%un(1:n)    = src(1:n)
  case('velocity_y');               n = min(size(src), size(this%vn));   if (n>0) this%vn(1:n)    = src(1:n)
  case('rain_rate');                n = min(size(src), size(rain_store)); if (n>0) rain_store(1:n) = src(1:n)
  case default; status = BMI_FAILURE; return
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
  select case (trim(name))
  case('velocity_x')
    n = min(size(src), size(inds))
    do k = 1, n
      i0 = inds(k)
      if (i0>=1 .and. i0<=size(this%un)) this%un(i0) = real(src(k), c_double)
    end do
  case('velocity_y')
    n = min(size(src), size(inds))
    do k = 1, n
      i0 = inds(k)
      if (i0>=1 .and. i0<=size(this%vn)) this%vn(i0) = real(src(k), c_double)
    end do
  case('rain_rate')
    n = min(size(src), size(inds))
    do k = 1, n
      i0 = inds(k)
      if (i0>=1 .and. i0<=size(rain_store)) rain_store(i0) = real(src(k), c_double)
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
  select case (trim(name))
  case('water_surface_elevation')
    n = min(size(src), size(inds))
    do k = 1, n
      i0 = inds(k)
      if (i0>=1 .and. i0<=size(this%z)) this%z(i0) = src(k)
    end do
  case('water_depth')
    n = min(size(src), size(inds))
    do k = 1, n
      i0 = inds(k)
      if (i0>=1 .and. i0<=size(this%h)) this%h(i0) = src(k)
    end do
  case('velocity_x')
    n = min(size(src), size(inds))
    do k = 1, n
      i0 = inds(k)
      if (i0>=1 .and. i0<=size(this%un)) this%un(i0) = src(k)
    end do
  case('velocity_y')
    n = min(size(src), size(inds))
    do k = 1, n
      i0 = inds(k)
      if (i0>=1 .and. i0<=size(this%vn)) this%vn(i0) = src(k)
    end do
  case('rain_rate')
    n = min(size(src), size(inds))
    do k = 1, n
      i0 = inds(k)
      if (i0>=1 .and. i0<=size(rain_store)) rain_store(i0) = src(k)
    end do
  case default
    status = BMI_FAILURE; return
  end select
  status = BMI_SUCCESS
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
  case('water_surface_elevation');  n = size(this%z)
  case('water_depth');              n = size(this%h)
  case('velocity_x');               n = size(this%un)
  case('velocity_y');               n = size(this%vn)
  case('rain_rate');                n = size(rain_store)
  case default;                     n = -1
  end select
end function var_size

!----------------------------------------------------------------------
! Config reader (NAMELIST)
!----------------------------------------------------------------------
subroutine read_init_config(this, config_file, status)
  use, intrinsic :: iso_fortran_env, only: stderr => error_unit
  implicit none
  class(sfincs_bmi), intent(inout) :: this
  character(len=*),  intent(in)    :: config_file
  integer,           intent(out)   :: status

  ! locals
  integer :: rc, fu
  logical :: exists
  character(len=1000) :: line

  ! namelist items (provide sensible defaults)
  double precision :: model_start_time, model_end_time, time_step_size
  integer          :: num_time_steps

  namelist /sfincs/ model_start_time, model_end_time, time_step_size, num_time_steps

  ! defaults if not present in file
  model_start_time = 0.0d0
  model_end_time   = -1.0d0
  time_step_size   = 60.0d0
  num_time_steps   = -1

  ! ensure file exists
  inquire(file=trim(config_file), exist=exists)
  if (.not. exists) then
    write (stderr, '(3a)') 'Error: input file "', trim(config_file), '" does not exist.'
    status = BMI_FAILURE
    return
  end if

  ! open and read NAMELIST
  open (action='read', file=trim(config_file), iostat=rc, newunit=fu)
  if (rc /= 0) then
    write (stderr,'(3a,i0)') 'Error: cannot open "', trim(config_file), '", iostat=', rc
    status = BMI_FAILURE
    return
  end if

  read (nml=sfincs, iostat=rc, unit=fu)
  if (rc /= 0) then
    backspace(fu)
    read(fu,fmt='(A)') line
    write(stderr,'(A)') 'Invalid line in namelist: '//trim(line)
    write(stderr,'(3a)') 'Error: invalid namelist in "', trim(config_file), '".'
    close(fu)
    status = BMI_FAILURE
    return
  end if
  close(fu)

  ! derive end time from num_time_steps if needed
  if (model_end_time < 0.0d0 .and. num_time_steps >= 0) then
    model_end_time = model_start_time + dble(num_time_steps) * time_step_size
  end if

  ! basic validation
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

  ! apply to module-scope clock (used by the BMI methods)
  start_time_s   = model_start_time
  end_time_s     = model_end_time
  dt_s           = time_step_size
  current_time_s = start_time_s

  status = BMI_SUCCESS
end subroutine read_init_config

subroutine assert_ok(condition, msg)
  logical, intent(in) :: condition
  character(len=*), intent(in), optional :: msg
  if (.not. condition) then
    if (present(msg)) then
      print *, "Assertion failed:", trim(msg)
    else
      print *, "Assertion failed."
    end if
    ! For production, prefer returning BMI_FAILURE from caller instead of stopping.
  end if
end subroutine assert_ok

end module sfincs_bmi2

