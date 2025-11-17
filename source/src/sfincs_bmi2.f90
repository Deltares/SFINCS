module sfincs_bmi2
  !! BMI 2.0 wrapper for SFINCS (rectilinear grid, real implementation)
  !! Exposes outputs: zs, zb, depth (= max(zs - zb, 0))
  use, intrinsic :: iso_fortran_env, only: real32, real64
  use bmif_2_0,  only: bmi, BMI_SUCCESS, BMI_FAILURE
  use sfincs_data, only: nmax, mmax, np, zs, zb, t0, t1   ! adjust if needed
  use sfincs_lib,  only: sfincs_initialize, sfincs_update, sfincs_finalize, t, dt
  implicit none
  private

  integer, parameter :: GRID_ID = 1

  character(len=*), parameter :: VAR_ZS    = 'zs'
  character(len=*), parameter :: VAR_ZB    = 'zb'
  character(len=*), parameter :: VAR_DEPTH = 'depth'

  ! Module-level TARGET buffers (for BMI pointer getters)
  real(real32), target, allocatable, save :: g_zs_f(:), g_zb_f(:), g_depth_f(:)
  real(real64), target, allocatable, save :: g_zs_d(:), g_zb_d(:), g_depth_d(:)

  ! Module-level strings / arrays used for pointer-string returns
  character(len=:), allocatable, target, save :: g_component_name
  character(len=:), allocatable, target, save :: g_time_units
  character(len=:), allocatable, target, save :: g_input_names(:), g_output_names(:)
  character(len=:), allocatable, target, save :: g_loc_node, g_grid_rect, g_units_m

  type, extends(bmi) :: sfincs_bmi
    ! --- Grid & time ---
    integer      :: nx      = 0
    integer      :: ny      = 0
    real(real64) :: dx      = 0.0d0
    real(real64) :: dy      = 0.0d0
    real(real64) :: x0      = 0.0d0
    real(real64) :: y0      = 0.0d0
    real(real64) :: rotation= 0.0d0
    real(real64) :: t       = 0.0d0
    real(real64) :: dt      = 0.0d0
    real(real64) :: t_start = 0.0d0
    real(real64) :: t_end   = 0.0d0

    ! --- Model state (allocatable; not TARGET, but mirrored from sfincs_data)
    real(real32), allocatable :: zs    (:)
    real(real32), allocatable :: zb    (:)
    real(real32), allocatable :: depth (:)

    ! --- Double mirrors
    real(real64), allocatable :: zs_d    (:)
    real(real64), allocatable :: zb_d    (:)
    real(real64), allocatable :: depth_d (:)

    ! --- BMI pointer-return strings/arrays
    character(len=:), pointer :: component_name => null()
    character(len=:), pointer :: input_names(:)  => null()
    character(len=:), pointer :: output_names(:) => null()

    logical :: is_initialized   = .false.
    logical :: depth_is_derived = .true.

  contains
    ! lifecycle
    procedure :: initialize                 => sfincs_bmi_initialize
    procedure :: update                     => sfincs_bmi_update
    procedure :: update_until               => sfincs_bmi_update_until
    procedure :: finalize                   => sfincs_bmi_finalize

    ! model info
    procedure :: get_component_name         => sfincs_bmi_get_component_name
    procedure :: get_input_item_count       => sfincs_bmi_get_input_item_count
    procedure :: get_output_item_count      => sfincs_bmi_get_output_item_count
    procedure :: get_input_var_names        => sfincs_bmi_get_input_var_names
    procedure :: get_output_var_names       => sfincs_bmi_get_output_var_names

    ! time
    procedure :: get_start_time             => sfincs_bmi_get_start_time
    procedure :: get_end_time               => sfincs_bmi_get_end_time
    procedure :: get_current_time           => sfincs_bmi_get_current_time
    procedure :: get_time_step              => sfincs_bmi_get_time_step
    procedure :: get_time_units             => sfincs_bmi_get_time_units

    ! var metadata
    procedure :: get_var_grid               => sfincs_bmi_get_var_grid
    procedure :: get_var_type               => sfincs_bmi_get_var_type
    procedure :: get_var_units              => sfincs_bmi_get_var_units
    procedure :: get_var_itemsize           => sfincs_bmi_get_var_itemsize
    procedure :: get_var_nbytes             => sfincs_bmi_get_var_nbytes
    procedure :: get_var_location           => sfincs_bmi_get_var_location

    ! grid metadata (rectilinear + coords)
    procedure :: get_grid_rank              => sfincs_bmi_get_grid_rank
    procedure :: get_grid_size              => sfincs_bmi_get_grid_size
    procedure :: get_grid_type              => sfincs_bmi_get_grid_type
    procedure :: get_grid_shape             => sfincs_bmi_get_grid_shape
    procedure :: get_grid_spacing           => sfincs_bmi_get_grid_spacing
    procedure :: get_grid_origin            => sfincs_bmi_get_grid_origin
    procedure :: get_grid_x                 => sfincs_bmi_get_grid_x
    procedure :: get_grid_y                 => sfincs_bmi_get_grid_y
    procedure :: get_grid_z                 => sfincs_bmi_get_grid_z

    ! Unstructured hooks present in bmif_2_0 (stubbed here)
    procedure :: get_grid_edge_count        => sfincs_bmi_get_grid_edge_count
    procedure :: get_grid_face_count        => sfincs_bmi_get_grid_face_count
    procedure :: get_grid_node_count        => sfincs_bmi_get_grid_node_count
    procedure :: get_grid_edge_nodes        => sfincs_bmi_get_grid_edge_nodes
    procedure :: get_grid_face_nodes        => sfincs_bmi_get_grid_face_nodes
    procedure :: get_grid_nodes_per_face    => sfincs_bmi_get_grid_nodes_per_face
    procedure :: get_grid_face_edges        => sfincs_bmi_get_grid_face_edges

    ! values (float/double)
    procedure :: get_value_float            => sfincs_bmi_get_value_float
    procedure :: set_value_float            => sfincs_bmi_set_value_float
    procedure :: get_value_ptr_float        => sfincs_bmi_get_value_ptr_float
    procedure :: get_value_at_indices_float => sfincs_bmi_get_value_at_indices_float
    procedure :: set_value_at_indices_float => sfincs_bmi_set_value_at_indices_float

    procedure :: get_value_double           => sfincs_bmi_get_value_double
    procedure :: set_value_double           => sfincs_bmi_set_value_double
    procedure :: get_value_ptr_double       => sfincs_bmi_get_value_ptr_double
    procedure :: get_value_at_indices_double=> sfincs_bmi_get_value_at_indices_double
    procedure :: set_value_at_indices_double=> sfincs_bmi_set_value_at_indices_double

    ! integer suite (stubs; bmif_2_0 defers get_value_at_indices_int)
    procedure :: get_value_int              => sfincs_bmi_get_value_int
    procedure :: set_value_int              => sfincs_bmi_set_value_int
    procedure :: get_value_ptr_int          => sfincs_bmi_get_value_ptr_int
    procedure :: get_value_at_indices_int   => sfincs_bmi_get_value_at_indices_int
    procedure :: set_value_at_indices_int   => sfincs_bmi_set_value_at_indices_int
  end type sfincs_bmi

  public :: sfincs_bmi

contains

  !============================
  !== Lifecycle
  !============================
  function sfincs_bmi_initialize(this, config_file) result(status)
    !! BMI initialize: call real sfincs_initialize, then mirror state into BMI arrays.
    class(sfincs_bmi), intent(out) :: this
    character(len=*),  intent(in)  :: config_file
    integer :: status
    integer :: ierr
    integer :: ncell
    integer :: L
    character(len=:), allocatable :: cfg

    ! --------------------------
    ! 1) Call real SFINCS init
    ! --------------------------
    if (len_trim(config_file) == 0) then
      ! Default to local sfincs.inp if BMI caller passes empty string
      cfg = 'sfincs.inp'
    else
      cfg = trim(config_file)
    end if

    ierr = sfincs_initialize(cfg)
    if (ierr /= 0) then
      status = BMI_FAILURE
      return
    end if

    ! --------------------------
    ! 2) Mirror time info
    ! --------------------------
    this%t_start = t0   ! from sfincs_data
    this%t_end   = t1   ! from sfincs_data
    this%t       = t    ! from sfincs_lib (current time)

    ! Choose a BMI "logical" time step: small fraction of full range.
    if (this%t_end > this%t_start) then
      this%dt = max(1.0d0, (this%t_end - this%t_start) / 100.d0)
    else
      this%dt = 1.0d0
    end if

    ! --------------------------
    ! 3) Grid info
    ! --------------------------
    ! Use np as the number of cells. Represent as a 2D grid with ny=1, nx=np.
    ! This keeps things simple but still gives rank=2 as required by BMI.
    ncell    = np
    this%ny  = 1
    this%nx  = ncell
    this%dx  = 1.0d0
    this%dy  = 1.0d0
    this%x0  = 0.0d0
    this%y0  = 0.0d0
    this%rotation = 0.0d0

    ! --------------------------
    ! 4) Allocate BMI state and pull from core
    ! --------------------------
    allocate(this%zs(ncell), this%zb(ncell), this%depth(ncell))
    allocate(this%zs_d(ncell), this%zb_d(ncell), this%depth_d(ncell))

    call sync_from_sfincs_core(this)

    this%depth_is_derived = .true.

    ! --------------------------
    ! 5) Pointer-return strings/arrays
    ! --------------------------
    if (.not. allocated(g_component_name)) allocate(character(len=13) :: g_component_name)
    g_component_name = 'SFINCS BMI 2.0'

    if (.not. allocated(g_time_units))     allocate(character(len=1) :: g_time_units)
    g_time_units = 's'

    if (.not. allocated(g_units_m))        allocate(character(len=1) :: g_units_m)
    g_units_m = 'm'

    if (.not. allocated(g_loc_node))       allocate(character(len=4) :: g_loc_node)
    g_loc_node   = 'node'

    if (.not. allocated(g_grid_rect))      allocate(character(len=19) :: g_grid_rect)
    g_grid_rect  = 'uniform_rectilinear'

    ! No BMI inputs for now
    if (.not. allocated(g_input_names))    allocate(character(len=1) :: g_input_names(0))

    ! outputs: zs, zb, depth
    if (.not. allocated(g_output_names)) then
      L = max(len(VAR_ZS), max(len(VAR_ZB), len(VAR_DEPTH)))
      allocate(character(len=L) :: g_output_names(3))
    end if
    g_output_names(1) = VAR_ZS
    g_output_names(2) = VAR_ZB
    g_output_names(3) = VAR_DEPTH

    ! Bind pointer-return fields
    this%component_name => g_component_name
    this%input_names    => g_input_names
    this%output_names   => g_output_names

    ! Module-level TARGET buffers for pointer getters
    call ensure_module_ptr_buffers(this)

    this%is_initialized = .true.
    status = BMI_SUCCESS
  end function sfincs_bmi_initialize

  function sfincs_bmi_update(this) result(status)
    !! BMI update: advance model by ~this%dt using real sfincs_update.
    class(sfincs_bmi), intent(inout) :: this
    integer :: status
    integer :: ierr
    double precision :: dtrange

    if (.not. this%is_initialized) then
      status = BMI_FAILURE
      return
    end if

    write(*,*) 'BMI-DEBUG: entering sfincs_bmi_update'
    write(*,*) 'BMI-DEBUG: this%t_before = ', this%t
    write(*,*) 'BMI-DEBUG: this%dt       = ', this%dt

    ! SFINCS semantics: sfincs_update(dtrange) advances from t to t + dtrange.
    ! We use BMI's logical dt for dtrange.
    dtrange = this%dt
    write(*,*) 'BMI-DEBUG: calling sfincs_update with dtrange = ', dble(this%dt)
    ierr = sfincs_update(dtrange)
    write(*,*) 'BMI-DEBUG: sfincs_update returned ierr=', ierr
    write(*,*) 'BMI-DEBUG: core t, dt after update: t =', t, ' dt =', dt
    if (ierr /= 0) then
      status = BMI_FAILURE
      return
    end if

    ! Mirror time from SFINCS core back into BMI
    this%t = t          ! t from sfincs_lib
    ! t_start/t_end remain as originally set from t0,t1

    ! Pull updated zs/zb/depth from core into BMI buffers
    call sync_from_sfincs_core(this)

    status = BMI_SUCCESS
  end function sfincs_bmi_update

  function sfincs_bmi_update_until(this, time) result(status)
    !! BMI update_until: advance SFINCS until current_time >= time.
    class(sfincs_bmi), intent(inout) :: this
    double precision,  intent(in)    :: time
    integer :: status
    integer :: ierr
    double precision :: dtrange

    if (.not. this%is_initialized) then
      status = BMI_FAILURE
      return
    end if

    ! No-op if target time is not in the future.
    if (time <= this%t + 1.0d-12) then
      status = BMI_SUCCESS
      return
    end if

    dtrange = time - this%t

    ierr = sfincs_update(dtrange)
    if (ierr /= 0) then
      status = BMI_FAILURE
      return
    end if

    ! Mirror time and state from core
    this%t = t
    call sync_from_sfincs_core(this)

    status = BMI_SUCCESS
  end function sfincs_bmi_update_until

  function sfincs_bmi_finalize(this) result(status)
    !! BMI finalize: call real sfincs_finalize, then clean BMI buffers.
    class(sfincs_bmi), intent(inout) :: this
    integer :: status
    integer :: ierr

    if (.not. this%is_initialized) then
      ! Even if not initialized, treat as success for BMI semantics
      status = BMI_SUCCESS
      return
    end if

    ! Call the real SFINCS finalizer
    ierr = sfincs_finalize()

    ! BMI-level cleanup
    if (allocated(this%zs))      deallocate(this%zs)
    if (allocated(this%zb))      deallocate(this%zb)
    if (allocated(this%depth))   deallocate(this%depth)
    if (allocated(this%zs_d))    deallocate(this%zs_d)
    if (allocated(this%zb_d))    deallocate(this%zb_d)
    if (allocated(this%depth_d)) deallocate(this%depth_d)

    nullify(this%component_name, this%input_names, this%output_names)
    this%is_initialized   = .false.
    this%depth_is_derived = .true.

    if (ierr /= 0) then
      status = BMI_FAILURE
    else
      status = BMI_SUCCESS
    end if
  end function sfincs_bmi_finalize

  !============================
  !== Model info
  !============================

  function sfincs_bmi_get_component_name(this, name) result(status)
    class(sfincs_bmi),         intent(in)  :: this
    character(len=:), pointer, intent(out) :: name
    integer :: status
    if (associated(this%component_name)) then
      name => this%component_name; status = BMI_SUCCESS
    else
      nullify(name); status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_component_name

  function sfincs_bmi_get_input_var_names(this, names) result(status)
    class(sfincs_bmi),         intent(in)  :: this
    character(len=:), pointer, intent(out) :: names(:)
    integer :: status
    if (associated(this%input_names)) then
      names => this%input_names; status = BMI_SUCCESS
    else
      nullify(names); status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_input_var_names

  function sfincs_bmi_get_output_var_names(this, names) result(status)
    class(sfincs_bmi),         intent(in)  :: this
    character(len=:), pointer, intent(out) :: names(:)
    integer :: status

    ! Always return a fresh array of 3 names: ['zs','zb','depth']
    allocate(character(len=len(VAR_DEPTH)) :: names(3))
    names(1) = VAR_ZS
    names(2) = VAR_ZB
    names(3) = VAR_DEPTH
    status = BMI_SUCCESS
  end function sfincs_bmi_get_output_var_names

  function sfincs_bmi_get_input_item_count(this, count) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,          intent(out) :: count
    integer :: status
    count = 0; status = BMI_SUCCESS
  end function sfincs_bmi_get_input_item_count

  function sfincs_bmi_get_output_item_count(this, count) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,          intent(out) :: count
    integer :: status
    count = 3; status = BMI_SUCCESS
  end function sfincs_bmi_get_output_item_count

  !============================
  !== Time
  !============================
  function sfincs_bmi_get_start_time(this, time) result(status)
    class(sfincs_bmi), intent(in)  :: this
    double precision,  intent(out) :: time
    integer :: status
    time = this%t_start; status = BMI_SUCCESS
  end function sfincs_bmi_get_start_time

  function sfincs_bmi_get_end_time(this, time) result(status)
    class(sfincs_bmi), intent(in)  :: this
    double precision,  intent(out) :: time
    integer :: status
    time = this%t_end; status = BMI_SUCCESS
  end function sfincs_bmi_get_end_time

  function sfincs_bmi_get_current_time(this, time) result(status)
    class(sfincs_bmi), intent(in)  :: this
    double precision,  intent(out) :: time
    integer :: status
    time = this%t; status = BMI_SUCCESS
  end function sfincs_bmi_get_current_time

  function sfincs_bmi_get_time_step(this, time_step) result(status)
    class(sfincs_bmi), intent(in)  :: this
    double precision,  intent(out) :: time_step
    integer :: status
    time_step = this%dt; status = BMI_SUCCESS
  end function sfincs_bmi_get_time_step

  function sfincs_bmi_get_time_units(this, units) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(out) :: units
    integer :: status
    if (.not. allocated(g_time_units)) then
      allocate(character(len=1) :: g_time_units); g_time_units = 's'
    end if
    units = g_time_units
    status = BMI_SUCCESS
  end function sfincs_bmi_get_time_units

  !============================
  !== Var metadata
  !============================
  function sfincs_bmi_get_var_grid(this, name, grid) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    integer,           intent(out) :: grid
    integer :: status
    if (is_known_var(name)) then
      grid = GRID_ID; status = BMI_SUCCESS
    else
      grid = -1; status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_var_grid

  function sfincs_bmi_get_var_type(this, name, type) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    character(len=*),  intent(out) :: type
    integer :: status
    if (is_known_var(name)) then
      type = 'real'; status = BMI_SUCCESS
    else
      type = ''; status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_var_type

  function sfincs_bmi_get_var_units(this, name, units) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    character(len=*),  intent(out) :: units
    integer :: status
    if (.not. allocated(g_units_m)) then
      allocate(character(len=1) :: g_units_m); g_units_m = 'm'
    end if
    select case (trim(name))
    case (VAR_ZS, VAR_ZB, VAR_DEPTH)
      units = g_units_m; status = BMI_SUCCESS
    case default
      units = ''; status = BMI_FAILURE
    end select
  end function sfincs_bmi_get_var_units

  function sfincs_bmi_get_var_itemsize(this, name, size) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    integer,           intent(out) :: size
    integer :: status
    if (is_known_var(name)) then
      size = 4; status = BMI_SUCCESS
    else
      size = 0; status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_var_itemsize

  function sfincs_bmi_get_var_nbytes(this, name, nbytes) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    integer,           intent(out) :: nbytes
    integer :: status, sz
    status = this%get_var_itemsize(name, sz)
    if (status /= BMI_SUCCESS) then
      nbytes = 0; return
    end if
    nbytes = sz * this%nx * this%ny
  end function sfincs_bmi_get_var_nbytes

  function sfincs_bmi_get_var_location(this, name, location) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    character(len=*),  intent(out) :: location
    integer :: status
    if (.not. is_known_var(name)) then
      location = ''; status = BMI_FAILURE; return
    end if
    if (.not. allocated(g_loc_node)) then
      allocate(character(len=4) :: g_loc_node); g_loc_node = 'node'
    end if
    location = g_loc_node
    status = BMI_SUCCESS
  end function sfincs_bmi_get_var_location

  !============================
  !== Grid metadata (rectilinear)
  !============================
  function sfincs_bmi_get_grid_rank(this, grid, rank) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: rank
    integer :: status
    if (grid == GRID_ID) then
      rank = 2; status = BMI_SUCCESS
    else
      rank = 0; status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_grid_rank

  function sfincs_bmi_get_grid_size(this, grid, size) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: size
    integer :: status
    if (grid == GRID_ID) then
      size = this%nx * this%ny; status = BMI_SUCCESS
    else
      size = 0; status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_grid_size

  function sfincs_bmi_get_grid_type(this, grid, type) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    character(len=*),  intent(out) :: type
    integer :: status
    if (grid == GRID_ID) then
      type = 'uniform_rectilinear'; status = BMI_SUCCESS
    else
      type = ''; status = BMI_FAILURE
    end if
  end function sfincs_bmi_get_grid_type

  function sfincs_bmi_get_grid_shape(this, grid, shape) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: shape(:)
    integer :: status
    if (grid /= GRID_ID .or. size(shape) < 2) then
      status = BMI_FAILURE
    else
      shape(1) = this%ny
      shape(2) = this%nx
      status = BMI_SUCCESS
    end if
  end function sfincs_bmi_get_grid_shape

  function sfincs_bmi_get_grid_spacing(this, grid, spacing) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    double precision,  intent(out) :: spacing(:)
    integer :: status
    if (grid /= GRID_ID .or. size(spacing) < 2) then
      status = BMI_FAILURE
    else
      spacing(1) = this%dy
      spacing(2) = this%dx
      status = BMI_SUCCESS
    end if
  end function sfincs_bmi_get_grid_spacing

  function sfincs_bmi_get_grid_origin(this, grid, origin) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    double precision,  intent(out) :: origin(:)
    integer :: status
    if (grid /= GRID_ID .or. size(origin) < 2) then
      status = BMI_FAILURE
    else
      origin(1) = this%y0
      origin(2) = this%x0
      status = BMI_SUCCESS
    end if
  end function sfincs_bmi_get_grid_origin

  function sfincs_bmi_get_grid_x(this, grid, x) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    double precision,  intent(out) :: x(:)
    integer :: status, ncell, i, j, idx
    if (grid /= GRID_ID) then
      status = BMI_FAILURE; return
    end if
    ncell = this%nx * this%ny
    if (size(x) < ncell) then
      status = BMI_FAILURE; return
    end if
    idx = 0
    do j = 1, this%ny
      do i = 1, this%nx
        idx = idx + 1
        x(idx) = this%x0 + (i-0.5d0)*this%dx
      end do
    end do
    status = BMI_SUCCESS
  end function sfincs_bmi_get_grid_x

  function sfincs_bmi_get_grid_y(this, grid, y) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    double precision,  intent(out) :: y(:)
    integer :: status, ncell, i, j, idx
    if (grid /= GRID_ID) then
      status = BMI_FAILURE; return
    end if
    ncell = this%nx * this%ny
    if (size(y) < ncell) then
      status = BMI_FAILURE; return
    end if
    idx = 0
    do j = 1, this%ny
      do i = 1, this%nx
        idx = idx + 1
        y(idx) = this%y0 + (j-0.5d0)*this%dy
      end do
    end do
    status = BMI_SUCCESS
  end function sfincs_bmi_get_grid_y

  function sfincs_bmi_get_grid_z(this, grid, z) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    double precision,  intent(out) :: z(:)
    integer :: status, ncell
    if (grid /= GRID_ID) then
      status = BMI_FAILURE; return
    end if
    ncell = this%nx * this%ny
    if (size(z) < ncell) then
      status = BMI_FAILURE; return
    end if
    z(1:ncell) = 0.0d0   ! 2-D horizontal grid (flat z)
    status = BMI_SUCCESS
  end function sfincs_bmi_get_grid_z

  !============================
  !== Unstructured-only (stubs matching bmif_2_0)
  !============================
  function sfincs_bmi_get_grid_edge_count(this, grid, count) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: count
    integer :: status
    count = 0; status = BMI_FAILURE
  end function sfincs_bmi_get_grid_edge_count

  function sfincs_bmi_get_grid_face_count(this, grid, count) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: count
    integer :: status
    count = 0; status = BMI_FAILURE
  end function sfincs_bmi_get_grid_face_count

  function sfincs_bmi_get_grid_node_count(this, grid, count) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: count
    integer :: status
    count = 0; status = BMI_FAILURE
  end function sfincs_bmi_get_grid_node_count

  function sfincs_bmi_get_grid_edge_nodes(this, grid, edge_nodes) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: edge_nodes(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_get_grid_edge_nodes

  function sfincs_bmi_get_grid_face_nodes(this, grid, face_nodes) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: face_nodes(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_get_grid_face_nodes

  function sfincs_bmi_get_grid_nodes_per_face(this, grid, nodes_per_face) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: nodes_per_face(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_get_grid_nodes_per_face

  function sfincs_bmi_get_grid_face_edges(this, grid, face_edges) result(status)
    class(sfincs_bmi), intent(in)  :: this
    integer,           intent(in)  :: grid
    integer,           intent(out) :: face_edges(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_get_grid_face_edges

  !============================
  !== Values (FLOAT)
  !============================
  function sfincs_bmi_get_value_float(this, name, dest) result(status)
    class(sfincs_bmi), intent(in)    :: this
    character(len=*),  intent(in)    :: name
    real(real32),      intent(inout) :: dest(:)
    integer :: status

    select case (trim(name))
    case (VAR_ZS)
      if (size(dest) < size(this%zs)) then
        status = BMI_FAILURE
        return
      end if
      dest = this%zs

    case (VAR_ZB)
      if (size(dest) < size(this%zb)) then
        status = BMI_FAILURE
        return
      end if
      dest = this%zb

    case (VAR_DEPTH)
      if (size(dest) < size(this%zs)) then
        ! use zs size as canonical grid size
        status = BMI_FAILURE
        return
      end if

      if (this%depth_is_derived) then
        ! Depth is derived from zs and zb: compute on the fly
        dest = max(this%zs - this%zb, 0.0_real32)
      else
        ! Depth has been explicitly set by the user; use stored array
        if (.not. allocated(this%depth)) then
          status = BMI_FAILURE
          return
        end if
        dest = this%depth
      end if

    case default
      status = BMI_FAILURE
      return
    end select

    status = BMI_SUCCESS
  end function sfincs_bmi_get_value_float

  function sfincs_bmi_set_value_float(this, name, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    real(real32),      intent(in)    :: src(:)
    integer :: status

    select case (trim(name))
    case (VAR_ZS)
      if (size(src) < size(this%zs)) then
        status = BMI_FAILURE; return
      end if
      this%zs = src
      call sync_double_buffers(this)
      call sync_module_ptr_buffers(this)

    case (VAR_ZB)
      if (size(src) < size(this%zb)) then
        status = BMI_FAILURE; return
      end if
      this%zb = src
      call sync_double_buffers(this)
      call sync_module_ptr_buffers(this)

    case (VAR_DEPTH)
      if (size(src) < size(this%depth)) then
        status = BMI_FAILURE; return
      end if
      this%depth = src
      this%depth_d = real(this%depth, kind=real64)
      this%depth_is_derived = .false.
      call sync_module_ptr_buffers(this)

    case default
      status = BMI_FAILURE
      return
    end select

    status = BMI_SUCCESS
  end function sfincs_bmi_set_value_float

  function sfincs_bmi_get_value_ptr_float(this, name, dest_ptr) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    real(real32),      pointer, intent(inout) :: dest_ptr(:)
    integer :: status
    nullify(dest_ptr)
    select case (trim(name))
    case (VAR_ZS);    dest_ptr => g_zs_f
    case (VAR_ZB);    dest_ptr => g_zb_f
    case (VAR_DEPTH); dest_ptr => g_depth_f
    case default; status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_bmi_get_value_ptr_float

  function sfincs_bmi_get_value_at_indices_float(this, name, dest, inds) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    real(real32),      intent(inout) :: dest(:)
    integer,           intent(in)  :: inds(:)
    integer :: status, k, n
    real(real32), pointer :: p(:) => null()
    status = this%get_value_ptr_float(name, p); if (status /= BMI_SUCCESS) return
    n = size(inds); if (size(dest)<n) then; status=BMI_FAILURE; return; endif
    do k = 1, n; dest(k) = p(inds(k)); end do
    status = BMI_SUCCESS
  end function sfincs_bmi_get_value_at_indices_float

  function sfincs_bmi_set_value_at_indices_float(this, name, inds, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    integer,           intent(in)    :: inds(:)
    real(real32),      intent(in)    :: src(:)
    integer :: status, k, n
    real(real32), pointer :: p(:) => null()

    status = this%get_value_ptr_float(name, p)
    if (status /= BMI_SUCCESS) return

    n = size(inds)
    if (size(src) < n) then
      status = BMI_FAILURE
      return
    end if

    do k = 1, n
      p(inds(k)) = src(k)
    end do

    select case (trim(name))
    case (VAR_ZS)
      this%zs = g_zs_f
      call sync_double_buffers(this)

    case (VAR_ZB)
      this%zb = g_zb_f
      call sync_double_buffers(this)

    case (VAR_DEPTH)
      this%depth = g_depth_f
      this%depth_d = real(this%depth, kind=real64)
      this%depth_is_derived = .false.

    end select

    call sync_module_ptr_buffers(this)
    status = BMI_SUCCESS
  end function sfincs_bmi_set_value_at_indices_float

  !============================
  !== Values (DOUBLE)
  !============================
  function sfincs_bmi_get_value_double(this, name, dest) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    real(real64),      intent(inout) :: dest(:)
    integer :: status

    select case (trim(name))
    case (VAR_ZS)
      if (size(dest) < size(this%zs_d)) then
        status = BMI_FAILURE; return
      end if
      dest = this%zs_d

    case (VAR_ZB)
      if (size(dest) < size(this%zb_d)) then
        status = BMI_FAILURE; return
      end if
      dest = this%zb_d

    case (VAR_DEPTH)
      if (size(dest) < size(this%depth_d)) then
        status = BMI_FAILURE; return
      end if
      ! Always use the maintained double buffer for depth
      dest = this%depth_d

    case default
      status = BMI_FAILURE
      return
    end select

    status = BMI_SUCCESS
  end function sfincs_bmi_get_value_double

  function sfincs_bmi_set_value_double(this, name, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    real(real64),      intent(in)    :: src(:)
    integer :: status

    select case (trim(name))
    case (VAR_ZS)
      if (.not. allocated(this%zs_d) .or. size(src) < size(this%zs_d)) then
        status = BMI_FAILURE; return
      end if
      this%zs_d = src
      ! zs/zb define depth
      call sync_float_buffers(this)
      call sync_module_ptr_buffers(this)

    case (VAR_ZB)
      if (.not. allocated(this%zb_d) .or. size(src) < size(this%zb_d)) then
        status = BMI_FAILURE; return
      end if
      this%zb_d = src
      ! zs/zb define depth
      call sync_float_buffers(this)
      call sync_module_ptr_buffers(this)

    case (VAR_DEPTH)
      if (.not. allocated(this%depth_d) .or. size(src) < size(this%depth_d)) then
        status = BMI_FAILURE; return
      end if
      ! Here depth is explicitly set by the user: do NOT recompute from zs/zb
      this%depth_d = src
      this%depth   = real(src, kind=real32)
      this%depth_is_derived = .false.
      call sync_module_ptr_buffers(this)

    case default
      status = BMI_FAILURE
      return
    end select

    status = BMI_SUCCESS
  end function sfincs_bmi_set_value_double

  function sfincs_bmi_get_value_ptr_double(this, name, dest_ptr) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    real(real64),      pointer, intent(inout) :: dest_ptr(:)
    integer :: status
    nullify(dest_ptr)
    select case (trim(name))
    case (VAR_ZS);    dest_ptr => g_zs_d
    case (VAR_ZB);    dest_ptr => g_zb_d
    case (VAR_DEPTH); dest_ptr => g_depth_d
    case default; status = BMI_FAILURE; return
    end select
    status = BMI_SUCCESS
  end function sfincs_bmi_get_value_ptr_double

  function sfincs_bmi_get_value_at_indices_double(this, name, dest, inds) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    real(real64),      intent(inout) :: dest(:)
    integer,           intent(in)  :: inds(:)
    integer :: status, k, n
    real(real64), pointer :: p(:) => null()
    status = this%get_value_ptr_double(name, p); if (status /= BMI_SUCCESS) return
    n = size(inds); if (size(dest)<n) then; status=BMI_FAILURE; return; endif
    do k = 1, n; dest(k) = p(inds(k)); end do
    status = BMI_SUCCESS
  end function sfincs_bmi_get_value_at_indices_double

  function sfincs_bmi_set_value_at_indices_double(this, name, inds, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    integer,           intent(in)    :: inds(:)
    real(real64),      intent(in)    :: src(:)
    integer :: status, k, n
    real(real64), pointer :: p(:) => null()

    status = this%get_value_ptr_double(name, p)
    if (status /= BMI_SUCCESS) return

    n = size(inds)
    if (size(src) < n) then
      status = BMI_FAILURE
      return
    end if

    do k = 1, n
      p(inds(k)) = src(k)
    end do

    select case (trim(name))
    case (VAR_ZS)
      ! p => g_zs_d; keep zs_d in sync then recompute floats + depth
      this%zs_d = g_zs_d
      call sync_float_buffers(this)

    case (VAR_ZB)
      ! p => g_zb_d; keep zb_d in sync then recompute floats + depth
      this%zb_d = g_zb_d
      call sync_float_buffers(this)

    case (VAR_DEPTH)
      ! p => g_depth_d; user explicitly sets depth only
      this%depth_d = g_depth_d
      this%depth   = real(this%depth_d, kind=real32)
      this%depth_is_derived = .false.
      ! IMPORTANT: no recompute from zs/zb here

    case default
      status = BMI_FAILURE
      return
    end select

    call sync_module_ptr_buffers(this)
    status = BMI_SUCCESS
  end function sfincs_bmi_set_value_at_indices_double
  
  !============================
  !== Integer suite (stubs)
  !============================
  function sfincs_bmi_get_value_int(this, name, dest) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    integer,           intent(inout) :: dest(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_get_value_int

  function sfincs_bmi_set_value_int(this, name, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    integer,           intent(in)    :: src(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_set_value_int

  function sfincs_bmi_get_value_ptr_int(this, name, dest_ptr) result(status)
    class(sfincs_bmi), intent(in) :: this
    character(len=*),  intent(in) :: name
    integer,           pointer, intent(inout) :: dest_ptr(:)
    integer :: status
    nullify(dest_ptr)
    status = BMI_FAILURE
  end function sfincs_bmi_get_value_ptr_int

  function sfincs_bmi_get_value_at_indices_int(this, name, dest, inds) result(status)
    class(sfincs_bmi), intent(in)  :: this
    character(len=*),  intent(in)  :: name
    integer,           intent(inout) :: dest(:)
    integer,           intent(in)  :: inds(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_get_value_at_indices_int

  function sfincs_bmi_set_value_at_indices_int(this, name, inds, src) result(status)
    class(sfincs_bmi), intent(inout) :: this
    character(len=*),  intent(in)    :: name
    integer,           intent(in)    :: inds(:)
    integer,           intent(in)    :: src(:)
    integer :: status
    status = BMI_FAILURE
  end function sfincs_bmi_set_value_at_indices_int

  !============================
  !== Helpers
  !============================
  logical function is_known_var(name)
    character(len=*), intent(in) :: name
    select case (trim(name))
    case (VAR_ZS, VAR_ZB, VAR_DEPTH); is_known_var = .true.
    case default; is_known_var = .false.
    end select
  end function is_known_var

  ! Mirror core SFINCS arrays into BMI arrays
  subroutine sync_from_sfincs_core(this)
    !! Copy zs/zb from sfincs_data into BMI state and refresh mirrors/pointers.
    class(sfincs_bmi), intent(inout) :: this
    integer :: ncell

    ! We assume sfincs_data exports:
    !   integer :: np
    !   real(real32) :: zs(:), zb(:)
    !
    ! and that np is the length of zs/zb.
    ncell = np

    if (.not. allocated(this%zs)) then
      allocate(this%zs(ncell), this%zb(ncell), this%depth(ncell))
    end if
    if (.not. allocated(this%zs_d)) then
      allocate(this%zs_d(ncell), this%zb_d(ncell), this%depth_d(ncell))
    end if

    ! Copy real model state into BMI arrays
    this%zs    = zs(1:ncell)
    this%zb    = zb(1:ncell)
    this%depth = max(this%zs - this%zb, 0.0_real32)

    ! Maintain double-precision mirrors and module-level pointer buffers
    call sync_double_buffers(this)
    call ensure_module_ptr_buffers(this)
  end subroutine sync_from_sfincs_core

  subroutine ensure_module_ptr_buffers(this)
    class(sfincs_bmi), intent(in) :: this
    integer :: n
    n = this%nx * this%ny
    if (.not. allocated(g_zs_f))    allocate(g_zs_f(n))
    if (.not. allocated(g_zb_f))    allocate(g_zb_f(n))
    if (.not. allocated(g_depth_f)) allocate(g_depth_f(n))
    if (.not. allocated(g_zs_d))    allocate(g_zs_d(n))
    if (.not. allocated(g_zb_d))    allocate(g_zb_d(n))
    if (.not. allocated(g_depth_d)) allocate(g_depth_d(n))
    call sync_module_ptr_buffers(this)
  end subroutine ensure_module_ptr_buffers

  subroutine sync_module_ptr_buffers(this)
    class(sfincs_bmi), intent(in) :: this
    if (allocated(g_zs_f))    g_zs_f    = this%zs
    if (allocated(g_zb_f))    g_zb_f    = this%zb
    if (allocated(g_depth_f)) g_depth_f = this%depth
    if (allocated(g_zs_d))    g_zs_d    = this%zs_d
    if (allocated(g_zb_d))    g_zb_d    = this%zb_d
    if (allocated(g_depth_d)) g_depth_d = this%depth_d
  end subroutine sync_module_ptr_buffers

  subroutine sync_double_buffers(this)
    class(sfincs_bmi), intent(inout) :: this
    if (.not. allocated(this%zs_d)) return
    this%zs_d = real(this%zs, kind=real64)
    this%zb_d = real(this%zb, kind=real64)
    this%depth_d = max(this%zs_d - this%zb_d, 0.0_real64)
    this%depth_is_derived = .true.
  end subroutine sync_double_buffers

  subroutine sync_float_buffers(this)
    class(sfincs_bmi), intent(inout) :: this
    if (.not. allocated(this%zs)) return
    this%zs = real(this%zs_d, kind=real32)
    this%zb = real(this%zb_d, kind=real32)
    this%depth_d = max(this%zs_d - this%zb_d, 0.0_real64)
    this%depth   = real(this%depth_d, kind=real32)
    this%depth_is_derived = .true.
  end subroutine sync_float_buffers

end module sfincs_bmi2

