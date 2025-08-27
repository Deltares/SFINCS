! --------------------------------------------------------------------
! SFINCS BMI 2.0 driver test
! - Mirrors Snow17 test coverage
! - Matches your sfincs_bmi2.f90 interface (status-returning functions)
! - Uses pointer actuals for name/time strings (your dummies are POINTER)
! - Connectivity arrays are rank-1 (length 2*n)
! --------------------------------------------------------------------
program sfincs_bmi2_driver_test
  use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double
  use bmif_2_0
  use sfincs_bmi2
  implicit none

  ! -------- model ------------
  type(sfincs_bmi) :: m
  integer(c_int)   :: s

  ! -------- config -----------
  character(len=256) :: arg
  character(len=*), parameter :: default_cfg = 'sfincs_config.txt'

  ! -------- model info -------
  integer(c_int) :: nin, nout, nparams
  integer(c_int) :: i, j

  ! pointer actuals for pointer-dummy getters
  integer, parameter :: NAME_LEN = 2048
  character(len=NAME_LEN), target  :: cname_buf, tunits_buf
  character(len=NAME_LEN), pointer :: cname => null(), tunits => null()

  character(len=NAME_LEN), target  :: inputs_buf(256), outputs_buf(256), params_buf(256)
  character(len=NAME_LEN), pointer :: names_inputs(:)  => null()
  character(len=NAME_LEN), pointer :: names_outputs(:) => null()
  character(len=NAME_LEN), pointer :: names_params(:)  => null()
  character(len=NAME_LEN)          :: name

  ! var metadata
  character(len=NAME_LEN) :: vtype, vunits, vloc
  integer(c_int)          :: itemsize, nbytes, grid

  ! time
  real(c_double) :: t0, t, tend, dt, time_until, cur

  ! grid info (we keep 1-D to match your stub)
  integer(c_int) :: rank, gsize, nnode, nedge, nface
  integer(c_int) :: shape(1), nodes_per_face(1)
  real(c_double) :: spacing(1), origin(1)
  real(c_double), allocatable :: gx(:), gy(:), gz(:)

  ! topology (rank-1, length 2*n)
  integer(c_int), allocatable :: edge_nodes(:), face_edges(:), face_nodes(:)

  ! value buffers (typed)
  integer(c_int),     allocatable :: vi(:), inds(:)
  real(c_float),      allocatable :: vf(:)
  real(c_double),     allocatable :: vd(:)
  integer(c_int),     pointer     :: pi(:) => null()
  real(c_float),      pointer     :: pf(:) => null()
  real(c_double),     pointer     :: pd(:) => null()

  ! ---- parameters list (custom; adjust if you expose parameters in BMI later)
  nparams = 0
  nullify(names_params)

  ! ----------------- init -----------------
  call get_command_argument(1, arg)
  if (len_trim(arg) == 0) arg = default_cfg

  s = m%initialize(trim(arg)); call chk('initialize', s)

  ! ----------------- info -----------------
  cname  => cname_buf ; cname  = ''
  s = m%get_component_name(cname);    call chk('get_component_name', s)
  write(*,'(A,1X,A)') 'Component:', trim(cname)

  s = m%get_input_item_count(nin) ;   call chk('get_input_item_count', s)
  s = m%get_output_item_count(nout);  call chk('get_output_item_count', s)
  write(*,'(A,I0)') 'Input vars :', nin
  write(*,'(A,I0)') 'Output vars:', nout

  names_inputs  => inputs_buf(1:max(1,nin))
  names_outputs => outputs_buf(1:max(1,nout))

  s = m%get_input_var_names(names_inputs) ;   call chk('get_input_var_names', s)
  s = m%get_output_var_names(names_outputs);  call chk('get_output_var_names', s)

  if (nin  > 0) then
    write(*,*) 'Inputs:'
    do i=1,nin;  write(*,'(I3,2X,A)') i, trim(names_inputs(i));  end do
  end if
  if (nout > 0) then
    write(*,*) 'Outputs:'
    do i=1,nout; write(*,'(I3,2X,A)') i, trim(names_outputs(i)); end do
  end if

  tunits => tunits_buf ; tunits = ''
  s = m%get_time_units(tunits); call chk('get_time_units', s)

  s = m%get_start_time(t0);    call chk('get_start_time', s)
  s = m%get_current_time(t);   call chk('get_current_time', s)
  s = m%get_end_time(tend);    call chk('get_end_time', s)
  s = m%get_time_step(dt);     call chk('get_time_step', s)

  write(*,'(A,1X,A)') 'Time units:', trim(tunits)
  write(*,'("t0=",F10.3,"  t=",F10.3,"  tend=",F10.3,"  dt=",F10.3)') t0, t, tend, dt

  ! ----------------- var metadata -----------------
  call describe_vars('Inputs',  m, names_inputs,  nin)
  call describe_vars('Outputs', m, names_outputs, nout)

  ! ----------------- advance time -----------------
  time_until = min(t + 3.0d0*dt, tend)
  s = m%update_until(time_until); call chk('update_until', s)
  s = m%get_current_time(t);     call chk('get_current_time', s)
  write(*,'(A,F10.3)') 't after update_until:', t

  write(*,*) 'Stepping with update() until end_time...'
  do
    s = m%get_current_time(cur); call chk('get_current_time', s)
    if (cur >= tend) exit
    s = m%update();              call chk('update', s)
  end do
  s = m%get_current_time(t);     call chk('get_current_time', s)
  write(*,'(A,F10.3)') 't at end:', t

  ! ----------------- get/set tests (typed) -----------------
  call test_values(m, names_inputs,  nin)
  call test_values(m, names_outputs, nout)

  ! ----------------- grid tests -----------------
  if (nout > 0) then
    name = trim(names_outputs( max(1, min(1,nout)) ))
    s = m%get_var_grid(name, grid); call chk('get_var_grid', s)

    s = m%get_grid_rank(grid, rank);       call chk('get_grid_rank', s)
    s = m%get_grid_size(grid, gsize);      call chk('get_grid_size', s)
    s = m%get_grid_shape(grid, shape);     call chk('get_grid_shape', s)
    s = m%get_grid_spacing(grid, spacing); call chk('get_grid_spacing', s)
    s = m%get_grid_origin(grid, origin);   call chk('get_grid_origin', s)

    write(*,'(A,I0)')   'grid rank =', rank
    write(*,'(A,I0)')   'grid size =', gsize
    write(*,'(A,I0)')   'shape(1)  =', shape(1)
    write(*,'(A,F8.3)') 'spacing(1)=', spacing(1)
    write(*,'(A,F8.3)') 'origin(1) =', origin(1)

    allocate(gx(max(1,shape(1))), gy(max(1,shape(1))), gz(max(1,shape(1))))
    s = m%get_grid_x(grid, gx); call maybe(s,'get_grid_x')
    s = m%get_grid_y(grid, gy); call maybe(s,'get_grid_y')
    s = m%get_grid_z(grid, gz); call maybe(s,'get_grid_z')
    deallocate(gx, gy, gz)

    call maybe(m%get_grid_node_count(grid, nnode), 'get_grid_node_count')
    call maybe(m%get_grid_edge_count(grid, nedge), 'get_grid_edge_count')
    call maybe(m%get_grid_face_count(grid, nface), 'get_grid_face_count')
    write(*,'(A,3I0)') 'counts (node,edge,face)=', nnode, nedge, nface

    if (nedge > 0) then
      allocate(edge_nodes(2*nedge))
      s = m%get_grid_edge_nodes(grid, edge_nodes); call maybe(s,'get_grid_edge_nodes')
      deallocate(edge_nodes)
    end if
    if (nface > 0) then
      allocate(face_edges(2*nface), face_nodes(2*nface))
      s = m%get_grid_face_edges(grid, face_edges); call maybe(s,'get_grid_face_edges')
      s = m%get_grid_face_nodes(grid, face_nodes); call maybe(s,'get_grid_face_nodes')
      deallocate(face_edges, face_nodes)
      s = m%get_grid_nodes_per_face(grid, nodes_per_face); call maybe(s,'get_grid_nodes_per_face')
    end if
  end if

  ! ----------------- finalize -----------------
  s = m%finalize(); call chk('finalize', s)
  write(*,*) 'All done testing SFINCS BMI.'

contains

  subroutine chk(what, st)
    character(len=*), intent(in) :: what
    integer(c_int),   intent(in) :: st
    if (st /= BMI_SUCCESS) then
      write(*,'(A,": status=",I0)') trim(what), st
      stop 2
    end if
  end subroutine chk

  subroutine maybe(st, what)
    integer(c_int),   intent(in) :: st
    character(len=*), intent(in) :: what
    if (st /= BMI_SUCCESS) then
      write(*,'(A,": returned ",I0," (continuing)")') trim(what), st
    end if
  end subroutine maybe

  subroutine describe_vars(label, M, names, n)
    character(len=*), intent(in)       :: label
    class(sfincs_bmi), intent(inout)   :: M
    character(len=*),  intent(in)      :: names(:)
    integer(c_int),    intent(in)      :: n
    integer(c_int) :: i, s, grid, itemsize, nbytes
    character(len=NAME_LEN) :: vtype, vunits, vloc
    if (n <= 0) return
    write(*,'(A)') trim(label)//': metadata'
    do i = 1, n
      if (len_trim(names(i)) == 0) cycle
      s = M%get_var_type   (names(i), vtype);    call chk('get_var_type', s)
      s = M%get_var_units  (names(i), vunits);   call chk('get_var_units', s)
      s = M%get_var_itemsize(names(i), itemsize);call chk('get_var_itemsize', s)
      s = M%get_var_nbytes (names(i), nbytes);   call chk('get_var_nbytes', s)
      s = M%get_var_location(names(i), vloc);    call chk('get_var_location', s)
      s = M%get_var_grid   (names(i), grid);     call chk('get_var_grid', s)
      write(*,'(I3,2X,A)') i, trim(names(i))
      write(*,'(A,1X,A)')   '  type   =', trim(vtype)
      write(*,'(A,1X,A)')   '  units  =', trim(vunits)
      write(*,'(A,1X,A)')   '  loc    =', trim(vloc)
      write(*,'(A,I0)')     '  itemsz =', itemsize
      write(*,'(A,I0)')     '  nbytes =', nbytes
      write(*,'(A,I0)')     '  grid   =', grid
    end do
  end subroutine describe_vars

  subroutine test_values(M, names, n)
    class(sfincs_bmi), intent(inout) :: M
    character(len=*),  intent(in)    :: names(:)
    integer(c_int),    intent(in)    :: n

    integer(c_int) :: i, s, grid, gsize, k
    character(len=NAME_LEN) :: vtype

    do i = 1, n
      if (len_trim(names(i)) == 0) cycle

      s = M%get_var_type(names(i), vtype); call chk('get_var_type', s)
      s = M%get_var_grid(names(i), grid);  call chk('get_var_grid', s)
      s = M%get_grid_size(grid, gsize);    call chk('get_grid_size', s)

      select case (trim(adjustl(vtype)))
      case ('int','integer','int32','int64')
        allocate(vi(max(1,min(8,gsize))))
        vi = 0_c_int
        call maybe(M%get_value_int(names(i), vi), 'get_value_int')
        nullify(pi); call maybe(M%get_value_ptr_int(names(i), pi), 'get_value_ptr_int')
        do k=1,size(vi); vi(k)=k; end do
        call maybe(M%set_value_int(names(i), vi), 'set_value_int')

        allocate(inds(min(3, size(vi)))); do k=1,size(inds); inds(k)=k; end do
        call maybe(M%set_value_at_indices_int(names(i), inds, vi(1:size(inds))), 'set_value_at_indices_int')
        if (associated(pi)) nullify(pi)
        deallocate(vi, inds)

      case ('float','real','single')
        allocate(vf(max(1,min(8,gsize))))
        vf = 0.0_c_float
        call maybe(M%get_value_float(names(i), vf), 'get_value_float')
        nullify(pf); call maybe(M%get_value_ptr_float(names(i), pf), 'get_value_ptr_float')
        do k=1,size(vf); vf(k)=real(k,c_float); end do
        call maybe(M%set_value_float(names(i), vf), 'set_value_float')

        allocate(inds(min(3, size(vf)))); do k=1,size(inds); inds(k)=k; end do
        call maybe(M%set_value_at_indices_float(names(i), inds, vf(1:size(inds))), 'set_value_at_indices_float')
        if (associated(pf)) nullify(pf)
        deallocate(vf, inds)

      case ('double','float64','real*8','real(8)')
        allocate(vd(max(1,min(8,gsize))))
        vd = 0.0d0
        call maybe(M%get_value_double(names(i), vd), 'get_value_double')
        nullify(pd); call maybe(M%get_value_ptr_double(names(i), pd), 'get_value_ptr_double')
        do k=1,size(vd); vd(k)=real(k,c_double); end do
        call maybe(M%set_value_double(names(i), vd), 'set_value_double')

        allocate(inds(min(3, size(vd)))); do k=1,size(inds); inds(k)=k; end do
        call maybe(M%set_value_at_indices_double(names(i), inds, vd(1:size(inds))), 'set_value_at_indices_double')
        if (associated(pd)) nullify(pd)
        deallocate(vd, inds)

      case default
        write(*,'(A,1X,A)') 'Skipping unknown var type for', trim(names(i))
      end select
    end do
  end subroutine test_values

end program sfincs_bmi2_driver_test

