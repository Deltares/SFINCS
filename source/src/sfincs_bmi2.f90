module sfincs_bmi2
  use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double, c_sizeof
  use bmif_2_0,   only: bmi, BMI_SUCCESS, BMI_FAILURE, BMI_MAX_COMPONENT_NAME
  implicit none
  private

  !========================
  ! Module-scope storage
  !========================
  ! NOTE: These are simple stand-ins that the BMI pointer components will associate to.
  !       When wiring to real SFINCS, keep type components as POINTERs and just
  !       pointer-assign them to SFINCS arrays having TARGET in their own modules.

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
  type, public, extends(bmi) :: sfincs_bmi

    real(c_double)                    :: current_time = 0.0d0
    real(c_double)                    :: start_time   = 0.0d0
    real(c_double)                    :: end_time     = 0.0d0
    real(c_double)                    :: time_step    = 1.0d0

    ! pointer state buffers (never TARGET on components)
    integer(c_int),    pointer        :: state_i(:)  => null()
    real(c_float),     pointer        :: state_r(:)  => null()
    real(c_double),    pointer        :: state_d(:)  => null()
    
    character(len=:), allocatable :: component_name
    character(len=:), allocatable :: time_units
    character(len=:), allocatable :: input_names(:)
    character(len=:), allocatable :: output_names(:)

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
    procedure :: initialize                 => sfincs_initialize
    procedure :: update                     => sfincs_update
    procedure :: update_until               => sfincs_update_until
    procedure :: finalize                   => sfincs_finalize

    ! Component / I/O names
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

    ! Optional connectivity helpers exposed as flattened lists
    procedure :: get_grid_node_count        => sfincs_get_grid_node_count
    procedure :: get_grid_edge_count        => sfincs_get_grid_edge_count
    procedure :: get_grid_face_count        => sfincs_get_grid_face_count
    procedure :: get_grid_nodes_per_face    => sfincs_get_grid_nodes_per_face
    procedure :: get_grid_edge_nodes        => sfincs_get_grid_edge_nodes
    procedure :: get_grid_face_edges        => sfincs_get_grid_face_edges
    procedure :: get_grid_face_nodes        => sfincs_get_grid_face_nodes

    ! Value getters/setters (int/float/double)
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
!======================================================================
!                       LIFECYCLE
!======================================================================
  function sfincs_initialize(this, config_file) result(status)
    class(sfincs_bmi), intent(out) :: this
    character(len=*),  intent(in)    :: config_file
    integer(c_int)                   :: status
    integer                          :: i

    ! Allocate simple demo domain if not already; replace with real SFINCS init later
    if (nz <= 0) then
      nz = 10_c_int
      allocate(z_store(nz), h_store(nz), xz_store(nz), yz_store(nz))
      do i = 1, nz
        xz_store(i) = real(i-1, c_double)
        yz_store(i) = real(0, c_double)
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

    ! Pointer-associate type components
    this%z  => z_store
    this%h  => h_store
    this%un => un_store
    this%vn => vn_store
    this%xz => xz_store; this%yz => yz_store
    this%xu => xu_store; this%yu => yu_store
    this%xv => xv_store; this%yv => yv_store

    if (.not. allocated(this%component_name)) then
      allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%component_name)
      this%component_name = 'SFINCS-BMI'
    end if

    if (.not. allocated(this%time_units)) then
      allocate(character(len=16) :: this%time_units)
      this%time_units = 's'
    end if

    if (.not. allocated(this%input_names)) then
      ! zero-length list is okay; give it a dummy len so ALLOCATE is valid
      allocate(character(len=1) :: this%input_names(0))
    end if

    if (.not. allocated(this%output_names)) then
      allocate(character(len=BMI_MAX_COMPONENT_NAME) :: this%output_names(4))
      this%output_names = [ character(len=BMI_MAX_COMPONENT_NAME) :: &
           'water_surface_elevation', 'water_depth', 'velocity_x', 'velocity_y' ]
    end if

    ! Timebook
    current_time_s = 0.0d0
    start_time_s   = 0.0d0
    end_time_s     = 86400.0d0  ! 1 day default
    if (len_trim(config_file) > 0) then
      ! TODO: parse config_file to set dt, times, sizes; call real SFINCS init
    end if

    status = BMI_SUCCESS
  contains
    function make_name(s) result(p)
      character(len=*), intent(in) :: s
      character(len=:), pointer    :: p
      allocate(character(len=len_trim(s)) :: p)
      p = trim(s)
    end function make_name
  end function sfincs_initialize

  function sfincs_update(this) result(status)
    class(sfincs_bmi), intent(inout) :: this
    integer(c_int)                   :: status
    integer                          :: i

    ! Placeholder: shallow update (advect nothing, just a toy trend)
    do i = 1, size(this%z)
      this%z(i) = this%z(i) + 0.0d0
      this%h(i) = max(0.0d0, this%h(i) + 0.0d0)
    end do
    current_time_s = current_time_s + dt_s
    status = BMI_SUCCESS
  end function sfincs_update

  function sfincs_update_until(this, time) result(status)
    class(sfincs_bmi), intent(inout) :: this
    real(c_double),   intent(in)     :: time
    integer(c_int)                   :: status
    do while (current_time_s + dt_s <= time)
      status = sfincs_update(this)
      if (status /= BMI_SUCCESS) return
    end do
    status = BMI_SUCCESS
  end function sfincs_update_until

  function sfincs_finalize(this) result(status)
    class(sfincs_bmi), intent(inout) :: this
    integer(c_int)                   :: status

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

    if (allocated(this%component_name)) deallocate(this%component_name)
    if (allocated(this%time_units))     deallocate(this%time_units)

    if (allocated(this%input_names))  call deallocate(this%input_names)
    if (allocated(this%output_names)) call deallocate(this%output_names)

    !nullify(this%z, this%h, this%un, this%vn, this%xz, this%yz, this%xu, this%yu, this%xv, this%yv)
    !nullify(this%component_name, this%time_units, this%input_names, this%output_names)

    status = BMI_SUCCESS
  contains
    subroutine dealloc_name_list(list)
      character(len=:), pointer :: list(:)
      integer :: k
      do k = 1, size(list)
      !if (associated(list(k))) deallocate(list(k))
      end do
      !deallocate(list)
    end subroutine dealloc_name_list
  end function sfincs_finalize

!======================================================================
!                 COMPONENT / I/O VARIABLE NAMES
!======================================================================

  function sfincs_get_component_name(this, name) result(status)
    class(sfincs_bmi), intent(in)           :: this
    character(len=*), pointer, intent(out)  :: name
    integer(c_int)                          :: status
    if (allocated(this%component_name)) then
      name = trim(this%component_name)
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
    if (allocated(this%input_names)) then
      count = int(size(this%input_names), c_int)
    else
      count = 0_c_int
    end if
    status = BMI_SUCCESS
  end function sfincs_get_input_item_count

  function sfincs_get_output_item_count(this, count) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(out):: count
    integer(c_int)                :: status
    if (allocated(this%output_names)) then
      count = int(size(this%output_names), c_int)
    else
      count = 0_c_int
    end if
    status = BMI_SUCCESS
  end function sfincs_get_output_item_count

  function sfincs_get_input_var_names(this, names) result(status)
    class(sfincs_bmi), intent(in)        :: this
    character(len=*), pointer, intent(out)    :: names(:)
    integer(c_int)                        :: status
    integer                               :: k, n
    if (.not. allocated(this%input_names)) then
      status = BMI_FAILURE; return
    end if
    n = min(size(names), size(this%input_names))
    do k = 1, n
      call assign_trim(names(k), this%input_names(k))
    end do
    status = BMI_SUCCESS
  end function sfincs_get_input_var_names

  function sfincs_get_output_var_names(this, names) result(status)
    class(sfincs_bmi), intent(in)        :: this
    character(len=*), pointer, intent(out)    :: names(:)
    integer(c_int)                        :: status
    integer                               :: k, n
    if (.not. allocated(this%output_names)) then
      status = BMI_FAILURE; return
    end if
    n = min(size(names), size(this%output_names))
    do k = 1, n
      call assign_trim(names(k), this%output_names(k))
    end do
    status = BMI_SUCCESS
  end function sfincs_get_output_var_names

!======================================================================
!                               TIME
!======================================================================
  function sfincs_get_start_time(this, time) result(status)
    class(sfincs_bmi), intent(in) :: this
    real(c_double),    intent(out):: time
    integer(c_int)                :: status
    time = start_time_s
    status = BMI_SUCCESS
  end function sfincs_get_start_time

  function sfincs_get_end_time(this, time) result(status)
    class(sfincs_bmi), intent(in) :: this
    real(c_double),    intent(out):: time
    integer(c_int)                :: status
    time = end_time_s
    status = BMI_SUCCESS
  end function sfincs_get_end_time

  function sfincs_get_current_time(this, time) result(status)
    class(sfincs_bmi), intent(in) :: this
    real(c_double),    intent(out):: time
    integer(c_int)                :: status
    time = current_time_s
    status = BMI_SUCCESS
  end function sfincs_get_current_time

  function sfincs_get_time_step(this, time_step) result(status)
    class(sfincs_bmi), intent(in) :: this
    real(c_double),    intent(out):: time_step
    integer(c_int)                :: status
    time_step = dt_s
    status = BMI_SUCCESS
  end function sfincs_get_time_step

  function sfincs_get_time_units(this, units) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(out) :: units
    integer(c_int)                   :: status
    if (allocated(this%time_units)) then
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
    integer(c_int),    intent(out):: grid
    integer(c_int)                :: status
    select case (trim(name))
    case('water_surface_elevation','water_depth'); grid = 1_c_int   ! z-grid
    case('velocity_x');                              grid = 2_c_int   ! u-grid
    case('velocity_y');                              grid = 3_c_int   ! v-grid
    case default;                                    grid = -1_c_int; status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_get_var_grid

  function sfincs_get_var_type(this, name, type) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    character(len=*),  intent(out) :: type
    integer(c_int)                   :: status
    call assign_trim(type, 'double')   ! expose all as double for simplicity
    status = BMI_SUCCESS
  end function sfincs_get_var_type

  function sfincs_get_var_units(this, name, units) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    character(len=*),  intent(out) :: units
    integer(c_int)                   :: status
    select case (trim(name))
    case('water_surface_elevation','water_depth'); call assign_trim(units, 'm')
    case('velocity_x','velocity_y');               call assign_trim(units, 'm s-1')
    case default;                                  call assign_trim(units, '1'); status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_get_var_units

  function sfincs_get_var_itemsize(this, name, size) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    integer(c_int),    intent(out):: size
    integer(c_int)                :: status
    size = int(c_sizeof(0.0d0), c_int)  ! 8 bytes
    status = BMI_SUCCESS
  end function sfincs_get_var_itemsize

  function sfincs_get_var_nbytes(this, name, nbytes) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    integer(c_int),    intent(out):: nbytes
    integer(c_int)                :: status
    integer :: n
    n = var_size(this, name)
    if (n < 0) then
      nbytes = 0_c_int; status = BMI_FAILURE; return
    end if
    nbytes = int(8*n, c_int)
    status = BMI_SUCCESS
  end function sfincs_get_var_nbytes

  function sfincs_get_var_location(this, name, location) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    character(len=*),  intent(out) :: location
    integer(c_int)                   :: status
    select case (trim(name))
    case('water_surface_elevation','water_depth'); call assign_trim(location,'node')
    case('velocity_x','velocity_y');               call assign_trim(location,'edge')
    case default;                                  call assign_trim(location,'unknown'); status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_get_var_location

!======================================================================
!                            GRIDS
!======================================================================
  function sfincs_get_grid_type(this, grid, type) result(status)
    class(sfincs_bmi), intent(in)    :: this
    integer(c_int),    intent(in)    :: grid
    character(len=*),  intent(out) :: type
    integer(c_int)                   :: status
    call assign_trim(type, 'unstructured')
    status = BMI_SUCCESS
  end function sfincs_get_grid_type

  function sfincs_get_grid_rank(this, grid, rank) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: rank
    integer(c_int)                :: status
    rank = 1_c_int
    status = BMI_SUCCESS
  end function sfincs_get_grid_rank

  function sfincs_get_grid_size(this, grid, size) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: size
    integer(c_int)                :: status
    select case (grid)
    case (1_c_int); size = int(nz, c_int)
    case (2_c_int); size = int(nu, c_int)
    case (3_c_int); size = int(nv, c_int)
    case default;   size = 0_c_int; status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_get_grid_size

  function sfincs_get_grid_shape(this, grid, shape) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out) :: shape(:)
    integer(c_int)                :: status
    if (size(shape) >= 1) then
      select case (grid)
      case (1_c_int); shape(1) = int(nz, c_int)
      case (2_c_int); shape(1) = int(nu, c_int)
      case (3_c_int); shape(1) = int(nv, c_int)
      case default;   status = BMI_FAILURE; return
      end select
      status = BMI_SUCCESS
    else
      status = BMI_FAILURE
    end if
  end function sfincs_get_grid_shape

  function sfincs_get_grid_spacing(this, grid, spacing) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    real(c_double),    intent(out) :: spacing(:)
    integer(c_int)                :: status
    ! Unstructured: spacing undefined; return zeros
    if (size(spacing) >= 1) spacing(1) = 0.0d0
    status = BMI_SUCCESS
  end function sfincs_get_grid_spacing

  function sfincs_get_grid_origin(this, grid, origin) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    real(c_double),    intent(out) :: origin(:)
    integer(c_int)                :: status
    if (size(origin) >= 1) origin(1) = 0.0d0
    status = BMI_SUCCESS
  end function sfincs_get_grid_origin

  function sfincs_get_grid_x(this, grid, x) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    real(c_double),    intent(out) :: x(:)
    integer(c_int)                :: status
    integer :: n
    select case (grid)
    case (1_c_int); n = min(size(x), nz); if (n>0) x(1:n) = this%xz(1:n)
    case (2_c_int); n = min(size(x), nu); if (n>0) x(1:n) = this%xu(1:n)
    case (3_c_int); n = min(size(x), nv); if (n>0) x(1:n) = this%xv(1:n)
    case default; status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_get_grid_x

  function sfincs_get_grid_y(this, grid, y) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    real(c_double),    intent(out) :: y(:)
    integer(c_int)                :: status
    integer :: n
    select case (grid)
    case (1_c_int); n = min(size(y), nz); if (n>0) y(1:n) = this%yz(1:n)
    case (2_c_int); n = min(size(y), nu); if (n>0) y(1:n) = this%yu(1:n)
    case (3_c_int); n = min(size(y), nv); if (n>0) y(1:n) = this%yv(1:n)
    case default; status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_get_grid_y

  function sfincs_get_grid_z(this, grid, z) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    real(c_double),    intent(out) :: z(:)
    integer(c_int)                :: status
    if (size(z) >= 1) z(1) = 0.0d0
    status = BMI_SUCCESS
  end function sfincs_get_grid_z

  ! Connectivity: provide flattened 1-D outputs as requested (names + rank 1)
  function sfincs_get_grid_node_count(this, grid, count) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: count
    integer(c_int)                :: status
    select case(grid)
    case(1_c_int); count = int(nz, c_int)
    case(2_c_int); count = int(nu, c_int)
    case(3_c_int); count = int(nv, c_int)
    case default;  count = 0_c_int; status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_get_grid_node_count

  function sfincs_get_grid_edge_count(this, grid, count) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: count
    integer(c_int)                :: status
    count = 0_c_int
    status = BMI_SUCCESS
  end function sfincs_get_grid_edge_count

  function sfincs_get_grid_face_count(this, grid, count) result(status)
    class(sfincs_bmi), intent(in) :: this
    integer(c_int),    intent(in) :: grid
    integer(c_int),    intent(out):: count
    integer(c_int)                :: status
    count = 0_c_int
    status = BMI_SUCCESS
  end function sfincs_get_grid_face_count

  function sfincs_get_grid_nodes_per_face(this, grid, nodes_per_face) result(status)
    class(sfincs_bmi), intent(in)    :: this
    integer(c_int),    intent(in)    :: grid
    integer(c_int),    intent(out) :: nodes_per_face(:)
    integer(c_int)                   :: status
    if (size(nodes_per_face) > 0) nodes_per_face = 0_c_int
    status = BMI_SUCCESS
  end function sfincs_get_grid_nodes_per_face

  function sfincs_get_grid_edge_nodes(this, grid, edge_nodes) result(status)
    class(sfincs_bmi), intent(in)    :: this
    integer(c_int),    intent(in)    :: grid
    integer(c_int),    intent(out) :: edge_nodes(:)   ! flattened pairs
    integer(c_int)                   :: status
    if (size(edge_nodes) > 0) edge_nodes = 0_c_int
    status = BMI_SUCCESS
  end function sfincs_get_grid_edge_nodes

  function sfincs_get_grid_face_edges(this, grid, face_edges) result(status)
    class(sfincs_bmi), intent(in)    :: this
    integer(c_int),    intent(in)    :: grid
    integer(c_int),    intent(out) :: face_edges(:)   ! flattened pairs
    integer(c_int)                   :: status
    if (size(face_edges) > 0) face_edges = 0_c_int
    status = BMI_SUCCESS
  end function sfincs_get_grid_face_edges

  function sfincs_get_grid_face_nodes(this, grid, face_nodes) result(status)
    class(sfincs_bmi), intent(in)    :: this
    integer(c_int),    intent(in)    :: grid
    integer(c_int),    intent(out) :: face_nodes(:)   ! flattened list
    integer(c_int)                   :: status
    if (size(face_nodes) > 0) face_nodes = 0_c_int
    status = BMI_SUCCESS
  end function sfincs_get_grid_face_nodes

!======================================================================
!                        VALUE GETTERS
!======================================================================
  function sfincs_get_value_int(this, name, dest) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    integer(c_int),    intent(inout) :: dest(:)
    integer(c_int)                   :: status
    integer :: n
    select case (trim(name))
    case default
      status = BMI_FAILURE; return
    end select
    n = 0
    status = BMI_SUCCESS
  end function sfincs_get_value_int

  function sfincs_get_value_float(this, name, dest) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    real(c_float),     intent(inout) :: dest(:)
    integer(c_int)                   :: status
    integer :: n
    select case (trim(name))
    case('water_surface_elevation'); n = min(size(dest), size(this%z));  if (n>0) dest(1:n) = real(this%z(1:n), c_float)
    case('water_depth');              n = min(size(dest), size(this%h));  if (n>0) dest(1:n) = real(this%h(1:n), c_float)
    case('velocity_x');               n = min(size(dest), size(this%un)); if (n>0) dest(1:n) = real(this%un(1:n), c_float)
    case('velocity_y');               n = min(size(dest), size(this%vn)); if (n>0) dest(1:n) = real(this%vn(1:n), c_float)
    case default; status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_get_value_float

  function sfincs_get_value_double(this, name, dest) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    real(c_double),    intent(inout) :: dest(:)
    integer(c_int)                   :: status
    integer :: n
    select case (trim(name))
    case('water_surface_elevation'); n = min(size(dest), size(this%z));  if (n>0) dest(1:n) = this%z(1:n)
    case('water_depth');              n = min(size(dest), size(this%h));  if (n>0) dest(1:n) = this%h(1:n)
    case('velocity_x');               n = min(size(dest), size(this%un)); if (n>0) dest(1:n) = this%un(1:n)
    case('velocity_y');               n = min(size(dest), size(this%vn)); if (n>0) dest(1:n) = this%vn(1:n)
    case default; status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_get_value_double

!----------------------------------------------------------------------
!                       POINTER GETTERS
!----------------------------------------------------------------------

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
    select case (trim(name))
    case('water_surface_elevation'); dest_ptr => this%z
    case('water_depth');              dest_ptr => this%h
    case('velocity_x');               dest_ptr => this%un
    case('velocity_y');               dest_ptr => this%vn
    case default; nullify(dest_ptr); status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_get_value_ptr_double

!----------------------------------------------------------------------
!                 GET VALUE AT INDICES (1-based inds)
!----------------------------------------------------------------------
  function sfincs_get_value_at_indices_int(this, name, dest, inds) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    integer(c_int),    intent(inout) :: dest(:)
    integer(c_int),    intent(in)    :: inds(:)
    integer(c_int)                   :: status
    integer :: k, n
    n = min(size(dest), size(inds))
    do k = 1, n
      dest(k) = 0_c_int
    end do
    status = BMI_SUCCESS
  end function sfincs_get_value_at_indices_int

  function sfincs_get_value_at_indices_float(this, name, dest, inds) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    real(c_float),     intent(inout) :: dest(:)
    integer(c_int),    intent(in)    :: inds(:)
    integer(c_int)                   :: status
    integer :: k, n, i0
    n = min(size(dest), size(inds))
    do k = 1, n
      i0 = inds(k)
      select case (trim(name))
      case('water_surface_elevation'); if (i0>=1 .and. i0<=size(this%z))  dest(k) = real(this%z(i0),  c_float)
      case('water_depth');              if (i0>=1 .and. i0<=size(this%h))  dest(k) = real(this%h(i0),  c_float)
      case('velocity_x');               if (i0>=1 .and. i0<=size(this%un)) dest(k) = real(this%un(i0), c_float)
      case('velocity_y');               if (i0>=1 .and. i0<=size(this%vn)) dest(k) = real(this%vn(i0), c_float)
      case default
        status = BMI_FAILURE; return
      end select
    end do
    status = BMI_SUCCESS
  end function sfincs_get_value_at_indices_float

  function sfincs_get_value_at_indices_double(this, name, dest, inds) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    real(c_double),    intent(inout) :: dest(:)
    integer(c_int),    intent(in)    :: inds(:)
    integer(c_int)                   :: status
    integer :: k, n, i0
    n = min(size(dest), size(inds))
    do k = 1, n
      i0 = inds(k)
      select case (trim(name))
      case('water_surface_elevation'); if (i0>=1 .and. i0<=size(this%z))  dest(k) = this%z(i0)
      case('water_depth');              if (i0>=1 .and. i0<=size(this%h))  dest(k) = this%h(i0)
      case('velocity_x');               if (i0>=1 .and. i0<=size(this%un)) dest(k) = this%un(i0)
      case('velocity_y');               if (i0>=1 .and. i0<=size(this%vn)) dest(k) = this%vn(i0)
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
    integer(c_int),    intent(in)    :: src(:)
    integer(c_int)                   :: status
    status = BMI_FAILURE
  end function sfincs_set_value_int

  function sfincs_set_value_float(this, name, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    real(c_float),     intent(in)    :: src(:)
    integer(c_int)                   :: status
    integer :: n
    select case (trim(name))
    case('velocity_x'); n = min(size(src), size(this%un)); if (n>0) this%un(1:n) = real(src(1:n), c_double)
    case('velocity_y'); n = min(size(src), size(this%vn)); if (n>0) this%vn(1:n) = real(src(1:n), c_double)
    case default; status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_set_value_float

  function sfincs_set_value_double(this, name, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    real(c_double),    intent(in)    :: src(:)
    integer(c_int)                   :: status
    integer :: n
    select case (trim(name))
    case('water_surface_elevation'); n = min(size(src), size(this%z));  if (n>0) this%z(1:n)  = src(1:n)
    case('water_depth');              n = min(size(src), size(this%h));  if (n>0) this%h(1:n)  = src(1:n)
    case('velocity_x');               n = min(size(src), size(this%un)); if (n>0) this%un(1:n) = src(1:n)
    case('velocity_y');               n = min(size(src), size(this%vn)); if (n>0) this%vn(1:n) = src(1:n)
    case default; status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_set_value_double

  function sfincs_set_value_at_indices_int(this, name, inds, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    integer(c_int),    intent(in)    :: inds(:)
    integer(c_int),    intent(in)    :: src(:)
    integer(c_int)                   :: status
    status = BMI_FAILURE
  end function sfincs_set_value_at_indices_int

  function sfincs_set_value_at_indices_float(this, name, inds, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    integer(c_int),    intent(in)    :: inds(:)
    real(c_float),     intent(in)    :: src(:)
    integer(c_int)                   :: status
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
    case default
      status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_set_value_at_indices_float

  function sfincs_set_value_at_indices_double(this, name, inds, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    integer(c_int),    intent(in)    :: inds(:)
    real(c_double),    intent(in)    :: src(:)
    integer(c_int)                   :: status
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
    case('water_surface_elevation'); n = size(this%z)
    case('water_depth');              n = size(this%h)
    case('velocity_x');               n = size(this%un)
    case('velocity_y');               n = size(this%vn)
    case default;                     n = -1
    end select
  end function var_size

end module sfincs_bmi2

