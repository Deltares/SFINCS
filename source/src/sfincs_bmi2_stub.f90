module sfincs_bmi2
  use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double
  use bmif_2_0,  only: bmi, BMI_SUCCESS, BMI_FAILURE, BMI_MAX_COMPONENT_NAME
  implicit none
  private

  ! Concrete model that implements the bmif_2_0 abstract type
  type, public, extends(bmi) :: sfincs_bmi
    ! minimal internal state for a compilable stub implementation
    real(c_double)                    :: current_time = 0.0d0
    real(c_double)                    :: start_time   = 0.0d0
    real(c_double)                    :: end_time     = 0.0d0
    real(c_double)                    :: time_step    = 1.0d0

    ! pointer state buffers (never TARGET on components)
    integer(c_int),    pointer        :: state_i(:)  => null()
    real(c_float),     pointer        :: state_r(:)  => null()
    real(c_double),    pointer        :: state_d(:)  => null()

    ! component strings as deferred-length pointers
    character(len=:),  pointer        :: component_name => null()
    character(len=:),  pointer        :: time_units     => null()
    character(len=:),  pointer        :: input_names(:)  => null()
    character(len=:),  pointer        :: output_names(:) => null()
  contains
    ! lifecycle
    procedure :: initialize                => sfincs_initialize
    procedure :: update                    => sfincs_update
    procedure :: update_until              => sfincs_update_until
    procedure :: finalize                  => sfincs_finalize

    ! time
    procedure :: get_start_time            => sfincs_get_start_time
    procedure :: get_end_time              => sfincs_get_end_time
    procedure :: get_current_time          => sfincs_get_current_time
    procedure :: get_time_step             => sfincs_get_time_step
    procedure :: get_time_units            => sfincs_get_time_units

    ! component & var name helpers
    procedure :: get_component_name        => sfincs_get_component_name
    procedure :: get_input_item_count      => sfincs_get_input_item_count
    procedure :: get_output_item_count     => sfincs_get_output_item_count
    procedure :: get_input_var_names       => sfincs_get_input_var_names
    procedure :: get_output_var_names      => sfincs_get_output_var_names

    ! variable metadata
    procedure :: get_var_grid              => sfincs_get_var_grid
    procedure :: get_var_type              => sfincs_get_var_type
    procedure :: get_var_units             => sfincs_get_var_units
    procedure :: get_var_itemsize          => sfincs_get_var_itemsize
    procedure :: get_var_nbytes            => sfincs_get_var_nbytes
    procedure :: get_var_location          => sfincs_get_var_location

    ! grid metadata (unstructured 1D)
    procedure :: get_grid_type             => sfincs_get_grid_type
    procedure :: get_grid_rank             => sfincs_get_grid_rank
    procedure :: get_grid_size             => sfincs_get_grid_size
    procedure :: get_grid_shape            => sfincs_get_grid_shape
    procedure :: get_grid_spacing          => sfincs_get_grid_spacing
    procedure :: get_grid_origin           => sfincs_get_grid_origin
    procedure :: get_grid_x                => sfincs_get_grid_x
    procedure :: get_grid_y                => sfincs_get_grid_y
    procedure :: get_grid_z                => sfincs_get_grid_z
    procedure :: get_grid_node_count       => sfincs_get_grid_node_count
    procedure :: get_grid_edge_count       => sfincs_get_grid_edge_count
    procedure :: get_grid_face_count       => sfincs_get_grid_face_count
    procedure :: get_grid_nodes_per_face   => sfincs_get_grid_nodes_per_face
    procedure :: get_grid_edge_nodes       => sfincs_get_grid_edge_nodes
    procedure :: get_grid_face_edges       => sfincs_get_grid_face_edges
    procedure :: get_grid_face_nodes       => sfincs_get_grid_face_nodes

    ! value getters (dest is intent(inout) per bmif_2_0)
    procedure :: get_value_int             => sfincs_get_value_int
    procedure :: get_value_float           => sfincs_get_value_float
    procedure :: get_value_double          => sfincs_get_value_double
    procedure :: get_value_ptr_int         => sfincs_get_value_ptr_int
    procedure :: get_value_ptr_float       => sfincs_get_value_ptr_float
    procedure :: get_value_ptr_double      => sfincs_get_value_ptr_double
    procedure :: get_value_at_indices_int  => sfincs_get_value_at_indices_int
    procedure :: get_value_at_indices_float=> sfincs_get_value_at_indices_float
    procedure :: get_value_at_indices_double=> sfincs_get_value_at_indices_double

    ! value setters
    procedure :: set_value_int             => sfincs_set_value_int
    procedure :: set_value_float           => sfincs_set_value_float
    procedure :: set_value_double          => sfincs_set_value_double
    procedure :: set_value_at_indices_int  => sfincs_set_value_at_indices_int
    procedure :: set_value_at_indices_float=> sfincs_set_value_at_indices_float
    procedure :: set_value_at_indices_double=> sfincs_set_value_at_indices_double
  end type sfincs_bmi

contains
  !========================
  ! lifecycle
  !========================
  function sfincs_initialize(this, config_file) result(status)
    class(sfincs_bmi), intent(out) :: this
    character(len=*),  intent(in)    :: config_file
    integer(c_int)                   :: status

    integer :: nvars

    ! Initialize strings
    if (.not. associated(this%component_name)) then
      allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%component_name)
    end if
    this%component_name = 'SFINCS-BMI'

    if (.not. associated(this%time_units)) then
      allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%time_units)
    end if
    this%time_units = 's'

    ! Minimal names: empty arrays to be conformant
    if (.not. associated(this%input_names)) then
      allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%input_names(0))
    end if
    if (.not. associated(this%output_names)) then
      allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%output_names(0))
    end if

    ! Minimal state buffers
    if (.not. associated(this%state_i)) allocate(this%state_i(1)); this%state_i = 0
    if (.not. associated(this%state_r)) allocate(this%state_r(1)); this%state_r = 0.0_c_float
    if (.not. associated(this%state_d)) allocate(this%state_d(1)); this%state_d = 0.0d0

    this%start_time   = 0.0d0
    this%current_time = 0.0d0
    this%end_time     = 0.0d0
    this%time_step    = 1.0d0

    status = BMI_SUCCESS
  end function sfincs_initialize

  function sfincs_update(this) result(status)
    class(sfincs_bmi), intent(inout) :: this
    integer(c_int)                   :: status

    this%current_time = this%current_time + this%time_step
    status = BMI_SUCCESS
  end function sfincs_update

  function sfincs_update_until(this, time) result(status)
    class(sfincs_bmi), intent(inout) :: this
    real(c_double),   intent(in)     :: time
    integer(c_int)                   :: status

    do while (this%current_time + this%time_step <= time)
      this%current_time = this%current_time + this%time_step
    end do
    status = BMI_SUCCESS
  end function sfincs_update_until

  function sfincs_finalize(this) result(status)
    class(sfincs_bmi), intent(inout) :: this
    integer(c_int)                   :: status

    if (associated(this%state_i))      deallocate(this%state_i)
    if (associated(this%state_r))      deallocate(this%state_r)
    if (associated(this%state_d))      deallocate(this%state_d)
    if (associated(this%component_name)) deallocate(this%component_name)
    if (associated(this%time_units))     deallocate(this%time_units)
    if (associated(this%input_names))    deallocate(this%input_names)
    if (associated(this%output_names))   deallocate(this%output_names)

    nullify(this%state_i, this%state_r, this%state_d, &
            this%component_name, this%time_units, this%input_names, this%output_names)

    status = BMI_SUCCESS
  end function sfincs_finalize

  !========================
  ! time
  !========================
  function sfincs_get_start_time(this, time) result(status)
    class(sfincs_bmi), intent(in)  :: this
    real(c_double),    intent(out) :: time
    integer(c_int)                 :: status
    time   = this%start_time
    status = BMI_SUCCESS
  end function

  function sfincs_get_end_time(this, time) result(status)
    class(sfincs_bmi), intent(in)  :: this
    real(c_double),    intent(out) :: time
    integer(c_int)                 :: status
    time   = this%end_time
    status = BMI_SUCCESS
  end function

  function sfincs_get_current_time(this, time) result(status)
    class(sfincs_bmi), intent(in)  :: this
    real(c_double),    intent(out) :: time
    integer(c_int)                 :: status
    time   = this%current_time
    status = BMI_SUCCESS
  end function

  function sfincs_get_time_step(this, time_step) result(status)
    class(sfincs_bmi), intent(in)  :: this
    real(c_double),    intent(out) :: time_step
    integer(c_int)                 :: status
    time_step = this%time_step
    status    = BMI_SUCCESS
  end function

  function sfincs_get_time_units(this, units) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(out):: units     ! non-pointer dummy per bmif_2_0
    integer(c_int)                :: status
    if (associated(this%time_units)) then
      units = this%time_units
    else
      units = 's'
    end if
    status = BMI_SUCCESS
  end function

  !========================
  ! component & var names
  !========================
  function sfincs_get_component_name(this, name) result(status)
    class(sfincs_bmi), intent(in)           :: this
    character(len=*), pointer, intent(out)  :: name
    integer(c_int)                          :: status
    if (associated(this%component_name)) then
      name => this%component_name
      status = BMI_SUCCESS
    else
      nullify(name)
      status = BMI_FAILURE
    end if
  end function

  function sfincs_get_input_item_count(this, count) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(out):: count
    integer(c_int)                :: status
    if (associated(this%input_names)) then
      count = size(this%input_names)
    else
      count = 0_c_int
    end if
    status = BMI_SUCCESS
  end function

  function sfincs_get_output_item_count(this, count) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(out):: count
    integer(c_int)                :: status
    if (associated(this%output_names)) then
      count = size(this%output_names)
    else
      count = 0_c_int
    end if
    status = BMI_SUCCESS
  end function

  function sfincs_get_input_var_names(this, names) result(status)
    class(sfincs_bmi), intent(in)             :: this
    character(len=*), pointer, intent(out)    :: names(:)
    integer(c_int)                            :: status
    if (associated(this%input_names)) then
      names => this%input_names
      status = BMI_SUCCESS
    else
      nullify(names)
      status = BMI_FAILURE
    end if
  end function

  function sfincs_get_output_var_names(this, names) result(status)
    class(sfincs_bmi), intent(in)             :: this
    character(len=*), pointer, intent(out)    :: names(:)
    integer(c_int)                            :: status
    if (associated(this%output_names)) then
      names => this%output_names
      status = BMI_SUCCESS
    else
      nullify(names)
      status = BMI_FAILURE
    end if
  end function

  !========================
  ! variable metadata (stubs)
  !========================
  function sfincs_get_var_grid(this, name, grid) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    integer(c_int),    intent(out):: grid
    integer(c_int)                :: status
    grid   = 0_c_int
    status = BMI_SUCCESS
  end function

  function sfincs_get_var_type(this, name, type) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    character(len=*),  intent(out):: type
    integer(c_int)                :: status
    type   = 'double'
    status = BMI_SUCCESS
  end function

  function sfincs_get_var_units(this, name, units) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    character(len=*),  intent(out):: units
    integer(c_int)                :: status
    units  = '1'
    status = BMI_SUCCESS
  end function

  function sfincs_get_var_itemsize(this, name, size) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    integer(c_int),    intent(out):: size
    integer(c_int)                :: status
    size   = 8_c_int
    status = BMI_SUCCESS
  end function

  function sfincs_get_var_nbytes(this, name, nbytes) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    integer(c_int),    intent(out):: nbytes
    integer(c_int)                :: status
    nbytes = 8_c_int
    status = BMI_SUCCESS
  end function

  function sfincs_get_var_location(this, name, location) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    character(len=*),  intent(out):: location
    integer(c_int)                :: status
    location = 'node'
    status   = BMI_SUCCESS
  end function

  !========================
  ! grid metadata (stubs for unstructured 1D)
  !========================
  function sfincs_get_grid_type(this, grid, type) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    character(len=*),  intent(out):: type
    integer(c_int)                :: status
    type   = 'unstructured'
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_rank(this, grid, rank) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: rank
    integer(c_int)                :: status
    rank   = 1_c_int
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_size(this, grid, size) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: size
    integer(c_int)                :: status
    size   = 1_c_int
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_shape(this, grid, shape) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: shape(:)
    integer(c_int)                :: status
    if (size(shape) >= 1) shape(1) = 1_c_int
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_spacing(this, grid, spacing) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    real(c_double),    intent(out):: spacing(:)
    integer(c_int)                :: status
    if (size(spacing) >= 1) spacing(1) = 1.0d0
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_origin(this, grid, origin) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    real(c_double),    intent(out):: origin(:)
    integer(c_int)                :: status
    if (size(origin) >= 1) origin(1) = 0.0d0
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_x(this, grid, x) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    real(c_double),    intent(out):: x(:)
    integer(c_int)                :: status
    if (size(x) >= 1) x(1) = 0.0d0
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_y(this, grid, y) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    real(c_double),    intent(out):: y(:)
    integer(c_int)                :: status
    if (size(y) >= 1) y(1) = 0.0d0
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_z(this, grid, z) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    real(c_double),    intent(out):: z(:)
    integer(c_int)                :: status
    if (size(z) >= 1) z(1) = 0.0d0
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_node_count(this, grid, count) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: count
    integer(c_int)                :: status
    count  = 1_c_int
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_edge_count(this, grid, count) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: count
    integer(c_int)                :: status
    count  = 0_c_int
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_face_count(this, grid, count) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: count
    integer(c_int)                :: status
    count  = 0_c_int
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_nodes_per_face(this, grid, nodes_per_face) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: nodes_per_face(:)
    integer(c_int)                :: status
    if (size(nodes_per_face) > 0) nodes_per_face = 0_c_int
    status = BMI_SUCCESS
  end function

  ! NOTE: bmif_2_0 requires these connectivity arrays to be 1-D
  function sfincs_get_grid_edge_nodes(this, grid, edge_nodes) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: edge_nodes(:)
    integer(c_int)                :: status
    if (size(edge_nodes) > 0) edge_nodes = 0_c_int
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_face_edges(this, grid, face_edges) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: face_edges(:)
    integer(c_int)                :: status
    if (size(face_edges) > 0) face_edges = 0_c_int
    status = BMI_SUCCESS
  end function

  function sfincs_get_grid_face_nodes(this, grid, face_nodes) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: face_nodes(:)
    integer(c_int)                :: status
    if (size(face_nodes) > 0) face_nodes = 0_c_int
    status = BMI_SUCCESS
  end function

  !========================
  ! value getters (dest is INOUT per bmif_2_0)
  !========================
  function sfincs_get_value_int(this, name, dest) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    integer(c_int),    intent(inout) :: dest(:)
    integer(c_int)                   :: status
    if (associated(this%state_i)) then
      dest = 0_c_int
      dest(1:min(size(dest), size(this%state_i))) = this%state_i(1:min(size(dest), size(this%state_i)))
      status = BMI_SUCCESS
    else
      status = BMI_FAILURE
    end if
  end function

  function sfincs_get_value_float(this, name, dest) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    real(c_float),     intent(inout) :: dest(:)
    integer(c_int)                   :: status
    if (associated(this%state_r)) then
      dest = 0.0_c_float
      dest(1:min(size(dest), size(this%state_r))) = this%state_r(1:min(size(dest), size(this%state_r)))
      status = BMI_SUCCESS
    else
      status = BMI_FAILURE
    end if
  end function

  function sfincs_get_value_double(this, name, dest) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    real(c_double),    intent(inout) :: dest(:)
    integer(c_int)                   :: status
    if (associated(this%state_d)) then
      dest = 0.0d0
      dest(1:min(size(dest), size(this%state_d))) = this%state_d(1:min(size(dest), size(this%state_d)))
      status = BMI_SUCCESS
    else
      status = BMI_FAILURE
    end if
  end function

  function sfincs_get_value_ptr_int(this, name, dest_ptr) result(status)
    class(sfincs_bmi), intent(in)            :: this
    character(len=*),  intent(in)            :: name
    integer(c_int),    pointer, intent(inout)  :: dest_ptr(:)
    integer(c_int)                            :: status
    if (associated(this%state_i)) then
      dest_ptr => this%state_i
      status = BMI_SUCCESS
    else
      nullify(dest_ptr)
      status = BMI_FAILURE
    end if
  end function

  function sfincs_get_value_ptr_float(this, name, dest_ptr) result(status)
    class(sfincs_bmi), intent(in)            :: this
    character(len=*),  intent(in)            :: name
    real(c_float),     pointer, intent(inout)  :: dest_ptr(:)
    integer(c_int)                            :: status
    if (associated(this%state_r)) then
      dest_ptr => this%state_r
      status = BMI_SUCCESS
    else
      nullify(dest_ptr)
      status = BMI_FAILURE
    end if
  end function

  function sfincs_get_value_ptr_double(this, name, dest_ptr) result(status)
    class(sfincs_bmi), intent(in)             :: this
    character(len=*),  intent(in)             :: name
    real(c_double),    pointer, intent(inout)   :: dest_ptr(:)
    integer(c_int)                             :: status
    if (associated(this%state_d)) then
      dest_ptr => this%state_d
      status = BMI_SUCCESS
    else
      nullify(dest_ptr)
      status = BMI_FAILURE
    end if
  end function


  !========================
  ! value setters
  !========================
  function sfincs_set_value_int(this, name, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    integer(c_int),    intent(in)    :: src(:)
    integer(c_int)                   :: status
    integer                          :: need
    need = max(1, size(src))
    if (.not. associated(this%state_i)) then
      allocate(this%state_i(need))
    else if (size(this%state_i) /= need) then
      deallocate(this%state_i); allocate(this%state_i(need))
    end if
    this%state_i = src
    status = BMI_SUCCESS
  end function

  function sfincs_set_value_float(this, name, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    real(c_float),     intent(in)    :: src(:)
    integer(c_int)                   :: status
    integer                          :: need
    need = max(1, size(src))
    if (.not. associated(this%state_r)) then
      allocate(this%state_r(need))
    else if (size(this%state_r) /= need) then
      deallocate(this%state_r); allocate(this%state_r(need))
    end if
    this%state_r = src
    status = BMI_SUCCESS
  end function

  function sfincs_set_value_double(this, name, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    real(c_double),    intent(in)    :: src(:)
    integer(c_int)                   :: status
    integer                          :: need
    need = max(1, size(src))
    if (.not. associated(this%state_d)) then
      allocate(this%state_d(need))
    else if (size(this%state_d) /= need) then
      deallocate(this%state_d); allocate(this%state_d(need))
    end if
    this%state_d = src
    status = BMI_SUCCESS
  end function

 function sfincs_get_value_at_indices_int(this, name, dest, inds) result(status)
  use iso_c_binding, only: c_int
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in)    :: name
  integer(c_int),    intent(inout) :: dest(:)
  integer(c_int),    intent(in)    :: inds(:)
  integer(c_int)                   :: status
  integer                          :: k, n

  n = min(size(dest), size(inds))
  do k = 1, n
    if (associated(this%state_i) .and. inds(k) >= 1 .and. inds(k) <= size(this%state_i)) then
      dest(k) = this%state_i(inds(k))
    else
      dest(k) = 0_c_int
    end if
  end do

  status = BMI_SUCCESS
end function sfincs_get_value_at_indices_int


function sfincs_get_value_at_indices_float(this, name, dest, inds) result(status)
  use iso_c_binding, only: c_int, c_float
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in)    :: name
  real(c_float),     intent(inout) :: dest(:)
  integer(c_int),    intent(in)    :: inds(:)
  integer(c_int)                   :: status
  integer                          :: k, n

  n = min(size(dest), size(inds))
  do k = 1, n
    if (associated(this%state_r) .and. inds(k) >= 1 .and. inds(k) <= size(this%state_r)) then
      dest(k) = this%state_r(inds(k))
    else
      dest(k) = 0.0_c_float
    end if
  end do

  status = BMI_SUCCESS
end function sfincs_get_value_at_indices_float


function sfincs_get_value_at_indices_double(this, name, dest, inds) result(status)
  use iso_c_binding, only: c_int, c_double
  class(sfincs_bmi), intent(in) :: this
  character(len=*),  intent(in)    :: name
  real(c_double),    intent(inout) :: dest(:)
  integer(c_int),    intent(in)    :: inds(:)
  integer(c_int)                   :: status
  integer                          :: k, n

  n = min(size(dest), size(inds))
  do k = 1, n
    if (associated(this%state_d) .and. inds(k) >= 1 .and. inds(k) <= size(this%state_d)) then
      dest(k) = this%state_d(inds(k))
    else
      dest(k) = 0.0d0
    end if
  end do

  status = BMI_SUCCESS
end function sfincs_get_value_at_indices_double


function sfincs_set_value_at_indices_int(this, name, inds, src) result(status)
  use iso_c_binding, only: c_int
  class(sfincs_bmi), intent(inout) :: this
  character(len=*),  intent(in)    :: name
  integer(c_int),    intent(in)    :: inds(:)
  integer(c_int),    intent(in)    :: src(:)
  integer(c_int)                   :: status
  integer                          :: need, k, n, oldsz
  integer(c_int),    pointer       :: tmp(:)

  need = max(1, merge(maxval(inds), 1, size(inds) > 0))
  if (.not. associated(this%state_i)) then
    allocate(this%state_i(need)); this%state_i = 0_c_int
  else if (size(this%state_i) < need) then
    oldsz = size(this%state_i)
    allocate(tmp(need)); tmp = 0_c_int
    tmp(1:oldsz) = this%state_i
    deallocate(this%state_i)
    this%state_i => tmp
  end if

  n = min(size(src), size(inds))
  do k = 1, n
    if (inds(k) >= 1 .and. inds(k) <= size(this%state_i)) this%state_i(inds(k)) = src(k)
  end do
  status = BMI_SUCCESS
end function sfincs_set_value_at_indices_int


function sfincs_set_value_at_indices_float(this, name, inds, src) result(status)
  use iso_c_binding, only: c_int, c_float
  class(sfincs_bmi), intent(inout) :: this
  character(len=*),  intent(in)    :: name
  integer(c_int),    intent(in)    :: inds(:)
  real(c_float),     intent(in)    :: src(:)
  integer(c_int)                   :: status
  integer                          :: need, k, n, oldsz
  real(c_float),     pointer       :: tmp(:)

  need = max(1, merge(maxval(inds), 1, size(inds) > 0))
  if (.not. associated(this%state_r)) then
    allocate(this%state_r(need)); this%state_r = 0.0_c_float
  else if (size(this%state_r) < need) then
    oldsz = size(this%state_r)
    allocate(tmp(need)); tmp = 0.0_c_float
    tmp(1:oldsz) = this%state_r
    deallocate(this%state_r)
    this%state_r => tmp
  end if

  n = min(size(src), size(inds))
  do k = 1, n
    if (inds(k) >= 1 .and. inds(k) <= size(this%state_r)) this%state_r(inds(k)) = src(k)
  end do
  status = BMI_SUCCESS
end function sfincs_set_value_at_indices_float

function sfincs_set_value_at_indices_double(this, name, inds, src) result(status)
  use iso_c_binding, only: c_int, c_double
  class(sfincs_bmi), intent(inout) :: this
  character(len=*),  intent(in)    :: name
  integer(c_int),    intent(in)    :: inds(:)
  real(c_double),    intent(in)    :: src(:)
  integer(c_int)                   :: status
  integer                          :: need, k, n, oldsz
  real(c_double),    pointer       :: tmp(:)

  need = max(1, merge(maxval(inds), 1, size(inds) > 0))
  if (.not. associated(this%state_d)) then
    allocate(this%state_d(need)); this%state_d = 0.0d0
  else if (size(this%state_d) < need) then
    oldsz = size(this%state_d)
    allocate(tmp(need)); tmp = 0.0d0
    tmp(1:oldsz) = this%state_d
    deallocate(this%state_d)
    this%state_d => tmp
  end if

  n = min(size(src), size(inds))
  do k = 1, n
    if (inds(k) >= 1 .and. inds(k) <= size(this%state_d)) this%state_d(inds(k)) = src(k)
  end do
  status = BMI_SUCCESS
end function sfincs_set_value_at_indices_double


end module sfincs_bmi2

