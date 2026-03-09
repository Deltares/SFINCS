! ============================================================================
! Robust BMI 2.0 test harness for SFINCS (Fortran)
!
! - Requires your module:   use sfincs_bmi2
! - Requires BMI interface: use bmif_2_0
! - Usage:
!      ./sfincs_bmi2_driver_test <config_file>     (runs full test)
!      ./sfincs_bmi2_driver_test -v                (prints component name)
! - Notes:
!   * Pointers returned by BMI (names, component name, etc.) are NEVER
!     deallocated here; we just check association and print values.
!   * “dest” arrays in getters are allocated here and passed with INTENT(INOUT).
!   * Edge/face arrays are rank-2 (2, :) per BMI.
!   * Every call checked: prints “OK” or “FAIL(status=…)” and continues.
! ============================================================================

program sfincs_bmi2_driver_test
  use bmif_2_0
  use sfincs_bmi2
  implicit none

  ! -------------------------------
  ! Model & status
  ! -------------------------------
  type(sfincs_bmi) :: M
  integer(c_int)   :: s
  integer          :: argc
  character(len=512) :: arg, config

  ! -------------------------------
  ! Names/vars (POINTER dummies where BMI expects pointers)
  ! -------------------------------
  character(len=BMI_MAX_COMPONENT_NAME), pointer :: comp_name  => null()
  character(len=BMI_MAX_VAR_NAME),      pointer :: in_names(:) => null()
  character(len=BMI_MAX_VAR_NAME),      pointer :: out_names(:)=> null()
  integer :: n_in, n_out, n_all, i

  ! -------------------------------
  ! Time
  ! -------------------------------
  double precision :: t0, t, t_end, dt
  character(len=32) :: t_units

  ! -------------------------------
  ! Per-var descriptors
  ! -------------------------------
  character(len=BMI_MAX_VAR_NAME) :: vname
  character(len=64)  :: vtype
  character(len=32)  :: vunits
  integer(c_int)     :: itemsize, nbytes
  integer(c_int)     :: grid, grid_rank, grid_size
  character(len=64)  :: grid_type
  integer(c_int), allocatable :: shape(:)
  double precision,  allocatable :: spacing(:), origin(:), gx(:), gy(:), gz(:)

  ! -------------------------------
  ! Value buffers (for get/set)
  ! We allocate per-variable using grid_size where possible.
  ! -------------------------------
  integer(c_int),    allocatable, target :: buf_i(:)
  real(c_float),     allocatable, target :: buf_r(:)
  real(c_double),    allocatable, target :: buf_d(:)

  integer(c_int),    pointer :: p_i(:) => null()
  real(c_float),     pointer :: p_r(:) => null()
  real(c_double),    pointer :: p_d(:) => null()

  ! At-indices helpers
  integer(c_int), allocatable :: inds(:)
  integer                    :: n_pick, k

  ! Edge/face helpers (rank-2)
  integer(c_int) :: n_edges, n_faces, npf
  integer(c_int), allocatable :: edge_nodes(:,:)        ! (2, n_edges)
  integer(c_int), allocatable :: face_edges(:,:)        ! (2, ?)
  integer(c_int), allocatable :: face_nodes(:,:)        ! (2, ?)
  integer(c_int), allocatable :: nodes_per_face(:)      ! helper

  ! --------------------------------------------------------------------------
  ! CLI: config path or -v
  ! --------------------------------------------------------------------------
  argc = command_argument_count()
  if (argc < 1) then
    print *, 'Usage: sfincs_bmi2_driver_test [-v] <config_file>'
    stop 1
  end if
  call get_command_argument(1, arg)

  if (trim(arg) == '-v') then
    s = M%initialize('')             ! if your initialize requires a file, skip gracefully
    if (s == BMI_SUCCESS) then
      s = M%get_component_name(comp_name)
      if (s == BMI_SUCCESS .and. associated(comp_name)) then
        print *, trim(comp_name)
      else
        print *, 'SFINCS (BMI)'
      end if
      call quiet_finalize(M)
    else
      print *, 'SFINCS (BMI)'
    end if
    stop
  else
    config = trim(arg)
  end if

  ! --------------------------------------------------------------------------
  ! Initialize
  ! --------------------------------------------------------------------------
  call say('initialize("'//trim(config)//'")')
  s = M%initialize(config);            call chk('initialize', s)

  ! --------------------------------------------------------------------------
  ! Component name (POINTER dummy)
  ! --------------------------------------------------------------------------
  call say('get_component_name')
  s = M%get_component_name(comp_name); call chk('get_component_name', s)
  if (associated(comp_name)) then
    print *, '  component = ', trim(comp_name)
  else
    print *, '  component pointer not associated (OK if unimplemented)'
  end if

  ! --------------------------------------------------------------------------
  ! Input/output counts and names (POINTER dummies)
  ! --------------------------------------------------------------------------
  call say('get_input_item_count')
  s = M%get_input_item_count(n_in);    call chk('get_input_item_count', s)
  print *, '  n_input  = ', n_in

  call say('get_output_item_count')
  s = M%get_output_item_count(n_out);  call chk('get_output_item_count', s)
  print *, '  n_output = ', n_out

  call say('get_input_var_names (POINTER)')
  s = M%get_input_var_names(in_names); call chk('get_input_var_names', s)
  if (associated(in_names)) then
    do i=1, size(in_names)
      print *, '   input(', i, ')= ', trim(in_names(i))
    end do
  else
    print *, '   (no input names or not associated)'
  end if

  call say('get_output_var_names (POINTER)')
  s = M%get_output_var_names(out_names); call chk('get_output_var_names', s)
  if (associated(out_names)) then
    do i=1, size(out_names)
      print *, '  output(', i, ')= ', trim(out_names(i))
    end do
  else
    print *, '   (no output names or not associated)'
  end if

  n_all = size_with_null(in_names) + size_with_null(out_names)

  ! --------------------------------------------------------------------------
  ! Time info
  ! --------------------------------------------------------------------------
  call say('time metadata')
  s = M%get_start_time(t0);            call chk('get_start_time',   s)
  s = M%get_current_time(t);           call chk('get_current_time', s)
  s = M%get_end_time(t_end);           call chk('get_end_time',     s)
  s = M%get_time_step(dt);             call chk('get_time_step',    s)
  s = M%get_time_units(t_units);       call chk('get_time_units',   s)

  print '(A,1PE13.5,2A)', '  start = ', t0,   ' [', trim(t_units), ']'
  print '(A,1PE13.5,2A)', '  now   = ', t,    ' [', trim(t_units), ']'
  print '(A,1PE13.5,2A)', '  end   = ', t_end,' [', trim(t_units), ']'
  print '(A,1PE13.5,2A)', '  dt    = ', dt,   ' [', trim(t_units), ']'

  ! --------------------------------------------------------------------------
  ! Per-variable metadata & grid info
  ! --------------------------------------------------------------------------
  if (n_all > 0) then
    print *, '--- Variable metadata & grid checks ---'
  end if

  do i = 1, n_all
    vname = pick_name(i, in_names, out_names)
    if (len_trim(vname) == 0) cycle

    call say('var "'//trim(vname)//'"')

    s = M%get_var_type  (vname, vtype);   call chk('get_var_type',   s)
    s = M%get_var_units (vname, vunits);  call chk('get_var_units',  s)
    s = M%get_var_itemsize(vname, itemsize); call chk('get_var_itemsize', s)
    s = M%get_var_nbytes  (vname, nbytes);   call chk('get_var_nbytes',   s)

    print *, '  type=', trim(vtype), ' units=', trim(vunits), ' itemsize=', itemsize, ' nbytes=', nbytes

    ! Grid descriptor
    grid      = -1_c_int
    grid_rank = -1_c_int
    grid_size = -1_c_int
    grid_type = ''

    s = M%get_var_grid(vname, grid);      call opt('get_var_grid', s)
    if (s /= BMI_SUCCESS) cycle

    s = M%get_grid_type(grid, grid_type); call opt('get_grid_type', s)
    s = M%get_grid_rank(grid, grid_rank); call opt('get_grid_rank', s)
    s = M%get_grid_size(grid, grid_size); call opt('get_grid_size', s)

    print *, '   grid id=', grid, ' type=', trim(grid_type), ' rank=', grid_rank, ' size=', grid_size

    if (grid_rank > 0) then
      allocate(shape(grid_rank), spacing(grid_rank), origin(grid_rank), stat=s)
      if (s == 0) then
        s = M%get_grid_shape  (grid, shape);   call opt('get_grid_shape',   s)
        s = M%get_grid_spacing(grid, spacing); call opt('get_grid_spacing', s)
        s = M%get_grid_origin (grid, origin);  call opt('get_grid_origin',  s)
        call print_intvec ('     shape', shape)
        call print_realvec('   spacing', spacing)
        call print_realvec('    origin', origin)
      end if
    end if

    if (grid_size > 0) then
      allocate(gx(grid_size), gy(grid_size), gz(grid_size), stat=s)
      if (s == 0) then
        s = M%get_grid_x(grid, gx); call opt('get_grid_x', s)
        s = M%get_grid_y(grid, gy); call opt('get_grid_y', s)
        s = M%get_grid_z(grid, gz); call opt('get_grid_z', s)
        if (s == BMI_SUCCESS) then
          print *, '   x(1)=', gx(min(1, size(gx))), ' y(1)=', gy(min(1,size(gy))), ' z(1)=', gz(min(1,size(gz)))
        end if
      end if
    end if

    ! Edge/face topology (rank-2 dummies: (2, :))
    call topo_checks(M, grid)

    ! Value I/O checks (typed by vtype string)
    call value_roundtrip(M, vname, trim(lower(vtype)), grid_size)
  end do

  ! --------------------------------------------------------------------------
  ! Run forward in time
  ! --------------------------------------------------------------------------
  call say('time integration')
  s = M%get_current_time(t);  call chk('get_current_time', s)
  s = M%get_end_time(t_end);  call chk('get_end_time',     s)
  s = M%get_time_step(dt);    call chk('get_time_step',    s)

  do while (t < t_end - 0.5d0*dt)
    s = M%update();           call chk('update', s)
    s = M%get_current_time(t);call chk('get_current_time', s)
  end do
  if (t < t_end) then
    s = M%update_until(t_end);call chk('update_until', s)
  end if
  s = M%get_current_time(t);  call chk('get_current_time', s)
  print *, '  final time = ', t, ' [', trim(t_units), ']'

  ! --------------------------------------------------------------------------
  ! Finalize
  ! --------------------------------------------------------------------------
  call say('finalize')
  s = M%finalize();           call chk('finalize', s)
  print *, 'OK: all done.'

contains

  ! ========== helpers =========================================================

  subroutine say(msg)
    character(len=*), intent(in) :: msg
    print *, '--- ', trim(msg)
  end subroutine say

  subroutine chk(where, status)
    character(len=*), intent(in) :: where
    integer(c_int),   intent(in) :: status
    if (status == BMI_SUCCESS) then
      print *, '    ', trim(where), ' -> OK'
    else
      print *, '    ', trim(where), ' -> FAIL (status=', status, ')'
      stop 2
    end if
  end subroutine chk

  subroutine opt(where, status)
    character(len=*), intent(in) :: where
    integer(c_int),   intent(in) :: status
    if (status == BMI_SUCCESS) then
      print *, '    ', trim(where), ' -> OK'
    else
      print *, '    ', trim(where), ' -> (optional) FAIL (status=', status, ')'
    end if
  end subroutine opt

  integer function size_with_null(ptr_names) result(n)
    character(len=*), pointer, intent(in) :: ptr_names(:)
    if (associated(ptr_names)) then
      n = size(ptr_names)
    else
      n = 0
    end if
  end function size_with_null

  character(len=BMI_MAX_VAR_NAME) function pick_name(idx, ins, outs) result(name)
    integer, intent(in) :: idx
    character(len=*), pointer, intent(in) :: ins(:)
    character(len=*), pointer, intent(in) :: outs(:)
    integer :: ni, no
    ni = size_with_null(ins)
    no = size_with_null(outs)
    name = ''
    if (idx <= ni .and. ni > 0) then
      name = ins(idx)
    else if (idx > ni .and. (idx - ni) <= no .and. no > 0) then
      name = outs(idx - ni)
    end if
  end function pick_name

  subroutine print_intvec(tag, a)
    character(len=*), intent(in) :: tag
    integer(c_int),   intent(in) :: a(:)
    integer :: j
    write(*,'(A)',advance='no') '  '//trim(tag)//' = ['
    do j=1,size(a)
      write(*,'(I0)',advance='no') a(j)
      if (j < size(a)) write(*,'(A)',advance='no') ','
    end do
    write(*,'(A)') ']'
  end subroutine print_intvec

  subroutine print_realvec(tag, a)
    character(len=*), intent(in) :: tag
    double precision, intent(in) :: a(:)
    integer :: j
    write(*,'(A)',advance='no') '  '//trim(tag)//' = ['
    do j=1,size(a)
      write(*,'(1PE12.4)',advance='no') a(j)
      if (j < size(a)) write(*,'(A)',advance='no') ','
    end do
    write(*,'(A)') ']'
  end subroutine print_realvec

  pure character(len=len(s)) function lower(s)
    character(len=*), intent(in) :: s
    integer :: i, ich
    lower = s
    do i=1,len(s)
      ich = iachar(s(i:i))
      if (ich >= iachar('A') .and. ich <= iachar('Z')) then
        lower(i:i) = achar(ich + 32)
      end if
    end do
  end function lower

  subroutine value_roundtrip(M, name, vtype_lc, gsize)
    class(sfincs_bmi), intent(inout) :: M
    character(len=*),  intent(in)    :: name
    character(len=*),  intent(in)    :: vtype_lc
    integer(c_int),    intent(in)    :: gsize

    integer(c_int) :: s, n
    integer :: j

    n = max(1, gsize)

    select case (vtype_lc)
    case ('int','integer','int32')
      call say('get/set "'//trim(name)//'" (integer)')
      call alloc_or_resize_int(buf_i, n)

      ! get_value
      s = M%get_value(name, buf_i);                    call opt('get_value(int)', s)
      if (s == BMI_SUCCESS) print *, '    first=', buf_i(1)

      ! set_value (write back what we read)
      if (s == BMI_SUCCESS) then
        do j=1,n
          buf_i(j) = buf_i(j) + 1_c_int
        end do
        s = M%set_value(name, buf_i);                  call opt('set_value(int)', s)
        s = M%get_value(name, buf_i);                  call opt('get_value(int) after set', s)
        if (s == BMI_SUCCESS) print *, '    new first=', buf_i(1)
      end if

      ! get_value_ptr
      p_i => null()
      s = M%get_value_ptr(name, p_i);                  call opt('get_value_ptr(int)', s)
      if (s == BMI_SUCCESS .and. associated(p_i)) then
        print *, '    ptr first=', p_i(1)
      end if

      ! at-indices (if n >= 3)
      if (n >= 3) then
        allocate(inds(3)); inds = [1_c_int, 2_c_int, 3_c_int]
        s = M%get_value_at_indices(name, buf_i, inds); call opt('get_value_at_indices(int)', s)
        if (s == BMI_SUCCESS) print *, '    inds OK (n=3)'
        if (allocated(inds)) deallocate(inds)
      end if

    case ('real','float','real32')
      call say('get/set "'//trim(name)//'" (real32)')
      call alloc_or_resize_real(buf_r, n)

      s = M%get_value(name, buf_r);                    call opt('get_value(float)', s)
      if (s == BMI_SUCCESS) print *, '    first=', buf_r(1)

      if (s == BMI_SUCCESS) then
        do j=1,n
          buf_r(j) = buf_r(j) + 1.0_c_float
        end do
        s = M%set_value(name, buf_r);                  call opt('set_value(float)', s)
        s = M%get_value(name, buf_r);                  call opt('get_value(float) after set', s)
        if (s == BMI_SUCCESS) print *, '    new first=', buf_r(1)
      end if

      p_r => null()
      s = M%get_value_ptr(name, p_r);                  call opt('get_value_ptr(float)', s)
      if (s == BMI_SUCCESS .and. associated(p_r)) then
        print *, '    ptr first=', p_r(1)
      end if

      if (n >= 3) then
        allocate(inds(3)); inds = [1_c_int, 2_c_int, 3_c_int]
        s = M%get_value_at_indices(name, buf_r, inds); call opt('get_value_at_indices(float)', s)
        if (s == BMI_SUCCESS) print *, '    inds OK (n=3)'
        if (allocated(inds)) deallocate(inds)
      end if

    case ('double','real64','float64','double precision')
      call say('get/set "'//trim(name)//'" (real64)')
      call alloc_or_resize_double(buf_d, n)

      s = M%get_value(name, buf_d);                    call opt('get_value(double)', s)
      if (s == BMI_SUCCESS) print *, '    first=', buf_d(1)

      if (s == BMI_SUCCESS) then
        do j=1,n
          buf_d(j) = buf_d(j) + 1.0d0
        end do
        s = M%set_value(name, buf_d);                  call opt('set_value(double)', s)
        s = M%get_value(name, buf_d);                  call opt('get_value(double) after set', s)
        if (s == BMI_SUCCESS) print *, '    new first=', buf_d(1)
      end if

      p_d => null()
      s = M%get_value_ptr(name, p_d);                  call opt('get_value_ptr(double)', s)
      if (s == BMI_SUCCESS .and. associated(p_d)) then
        print *, '    ptr first=', p_d(1)
      end if

      if (n >= 3) then
        allocate(inds(3)); inds = [1_c_int, 2_c_int, 3_c_int]
        s = M%get_value_at_indices(name, buf_d, inds); call opt('get_value_at_indices(double)', s)
        if (s == BMI_SUCCESS) print *, '    inds OK (n=3)'
        if (allocated(inds)) deallocate(inds)
      end if

    case default
      call say('skip values for "'//trim(name)//'" (unknown type="'//trim(vtype_lc)//'")')
    end select
  end subroutine value_roundtrip

  subroutine alloc_or_resize_int(a, n)
    integer(c_int), allocatable, intent(inout), target :: a(:)
    integer, intent(in) :: n
    if (allocated(a)) then
      if (size(a) /= n) then
        deallocate(a)
        allocate(a(n))
      end if
    else
      allocate(a(n))
    end if
    a = 0_c_int
  end subroutine alloc_or_resize_int

  subroutine alloc_or_resize_real(a, n)
    real(c_float), allocatable, intent(inout), target :: a(:)
    integer, intent(in) :: n
    if (allocated(a)) then
      if (size(a) /= n) then
        deallocate(a)
        allocate(a(n))
      end if
    else
      allocate(a(n))
    end if
    a = 0.0_c_float
  end subroutine alloc_or_resize_real

  subroutine alloc_or_resize_double(a, n)
    real(c_double), allocatable, intent(inout), target :: a(:)
    integer, intent(in) :: n
    if (allocated(a)) then
      if (size(a) /= n) then
        deallocate(a)
        allocate(a(n))
      end if
    else
      allocate(a(n))
    end if
    a = 0.0d0
  end subroutine alloc_or_resize_double

  subroutine topo_checks(M, grid)
    class(sfincs_bmi), intent(inout) :: M
    integer(c_int),    intent(in)    :: grid
    integer(c_int) :: s, ne, nf, max_npf, j

    if (grid < 0_c_int) return

    ! Edge count
    ne = 0_c_int
    s = M%get_grid_edge_count(grid, ne);         call opt('get_grid_edge_count', s)
    if (s == BMI_SUCCESS .and. ne > 0) then
      allocate(edge_nodes(2, ne))
      edge_nodes = 0_c_int
      s = M%get_grid_edge_nodes(grid, edge_nodes); call opt('get_grid_edge_nodes(2,ne)', s)
      if (s == BMI_SUCCESS) then
        print *, '    edges: first edge nodes =', edge_nodes(1,1), edge_nodes(2,1)
      end if
    end if

    ! Face count
    nf = 0_c_int
    s = M%get_grid_face_count(grid, nf);         call opt('get_grid_face_count', s)
    if (s == BMI_SUCCESS .and. nf > 0) then
      ! nodes per face (max)
      max_npf = 0_c_int
      allocate(nodes_per_face(nf))
      nodes_per_face = 0_c_int
      s = M%get_grid_nodes_per_face(grid, nodes_per_face); call opt('get_grid_nodes_per_face', s)
      if (s == BMI_SUCCESS) then
        max_npf = 0_c_int
        do j=1,nf
          if (nodes_per_face(j) > max_npf) max_npf = nodes_per_face(j)
        end do
        if (max_npf <= 0_c_int) max_npf = 2_c_int
      else
        max_npf = 2_c_int
      end if

      allocate(face_edges(2, nf))
      allocate(face_nodes(2, nf))  ! allocate minimal 2×nf; model may only fill first K rows
      face_edges = 0_c_int
      face_nodes = 0_c_int

      s = M%get_grid_face_edges(grid, face_edges); call opt('get_grid_face_edges(2,nf)', s)
      s = M%get_grid_face_nodes(grid, face_nodes); call opt('get_grid_face_nodes(2,nf)', s)

      if (s == BMI_SUCCESS) then
        print *, '    faces: first face edges =', face_edges(1,1), face_edges(2,1)
        print *, '    faces: first face nodes =', face_nodes(1,1), face_nodes(2,1)
      end if
    end if
  end subroutine topo_checks

  subroutine quiet_finalize(M)
    class(sfincs_bmi), intent(inout) :: M
    integer(c_int) :: s
    s = M%finalize()
  end subroutine quiet_finalize

end program sfincs_bmi2_driver_test

