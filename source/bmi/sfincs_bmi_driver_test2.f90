! Robust BMI driver matching your current sfincs_bmi2.f90 interfaces
program test_sfincs_bmi_full
  use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double
  use bmif_2_0
  use sfincs_bmi2
  implicit none

  type(sfincs_bmi) :: M
  integer(c_int)   :: s, nin, nout, grid, rank, gsize, nnode, nedge, nface
  integer(c_int)   :: i, j, n, sz
  real(c_double)   :: t, dt, t0, t1, tend

  ! POINTER actuals for pointer-dummy getters
  character(len=2048), target  :: cname_buf, tunits_buf
  character(len=2048), pointer :: cname => null(), tunits => null()

  ! name lists (your get_*_var_names dummies are POINTER)
  character(len=2048), target  :: in_names_buf(64), out_names_buf(64)
  character(len=2048), pointer :: in_names(:)  => null()
  character(len=2048), pointer :: out_names(:) => null()

  ! small work arrays for value get/set tests
  integer(c_int),     allocatable :: wi(:)
  real(c_float),      allocatable :: wf(:)
  real(c_double),     allocatable :: wd(:)
  integer(c_int),     allocatable :: inds(:)
  integer(c_int),     pointer     :: pi(:)  => null()
  real(c_float),      pointer     :: pf(:)  => null()
  real(c_double),     pointer     :: pd(:)  => null()

  ! grid work buffers
  integer(c_int) :: shape(1), nodes_per_face(1)
  real(c_double) :: spacing(1), origin(1)
  real(c_double), allocatable :: gx(:), gy(:), gz(:)

  ! connectivity as RANK-1 arrays (2*n)
  integer(c_int), allocatable :: edge_nodes(:), face_edges(:), face_nodes(:)

  character(len=2048) :: vtype, vunits, vloc

  ! -- Initialize
  s = M%initialize('sfincs_config.txt'); call chk('initialize', s)

  ! -- Component & time info
  cname  => cname_buf   ; cname  = ''
  s = M%get_component_name(cname); call chk('get_component_name', s)
  write(*,'(A,1X,A)') 'component:', trim(cname)

  tunits => tunits_buf  ; tunits = ''
  s = M%get_time_units(tunits);     call chk('get_time_units', s)
  write(*,'(A,1X,A)') 'time units:', trim(tunits)

  s = M%get_start_time(t0);         call chk('get_start_time', s)
  s = M%get_current_time(t);        call chk('get_current_time', s)
  s = M%get_end_time(tend);         call chk('get_end_time', s)
  s = M%get_time_step(dt);          call chk('get_time_step', s)
  write(*,'("t0=",F10.3," t=",F10.3," tend=",F10.3," dt=",F10.3)') t0, t, tend, dt

  ! -- Var names & counts
  s = M%get_input_item_count(nin);   call chk('get_input_item_count', s)
  s = M%get_output_item_count(nout); call chk('get_output_item_count', s)

  nullify(in_names)   ; in_names  => in_names_buf(1:max(1,nin))
  nullify(out_names)  ; out_names => out_names_buf(1:max(1,nout))

  s = M%get_input_var_names(in_names);   call chk('get_input_var_names', s)
  s = M%get_output_var_names(out_names); call chk('get_output_var_names', s)

  if (nin  > 0) then
    write(*,'(A)') 'Inputs:';  do i=1,nin;  write(*,'(I3,2X,A)') i, trim(in_names(i));  end do
  end if
  if (nout > 0) then
    write(*,'(A)') 'Outputs:'; do i=1,nout; write(*,'(I3,2X,A)') i, trim(out_names(i)); end do
  end if

  call test_var_list(M, in_names,  nin)
  call test_var_list(M, out_names, nout)

  ! -- Advance model
  s = M%update();                call chk('update', s)
  s = M%get_current_time(t);     call chk('get_current_time', s)
  write(*,'(A,F10.3)') 't after one update:', t

  t1 = min(t + 3.0d0*dt, tend)
  s  = M%update_until(t1);       call chk('update_until', s)
  s  = M%get_current_time(t);    call chk('get_current_time', s)
  write(*,'(A,F10.3)') 't after update_until:', t

  ! -- Finalize
  s = M%finalize(); call chk('finalize', s)

contains

  subroutine chk(what, st)
    character(len=*), intent(in) :: what
    integer(c_int),   intent(in) :: st
    if (st /= BMI_SUCCESS) then
      write(*,'(A,": status=",I0)') trim(what), st
      stop 2
    end if
  end subroutine chk

  subroutine maybe(status, what)
    integer(c_int), intent(in) :: status
    character(len=*), intent(in) :: what
    if (status /= BMI_SUCCESS) then
      write(*,'(A,": returned ",I0," (continuing)")') trim(what), status
    end if
  end subroutine maybe

  subroutine test_var_list(M, names, nitems)
    class(sfincs_bmi), intent(inout)        :: M
    character(len=*),  intent(in)           :: names(:)
    integer(c_int),    intent(in)           :: nitems

    integer(c_int) :: i, s, grid, rank, gsize, nnode, nedge, nface
    integer(c_int) :: shape(1), nodes_per_face(1)
    real(c_double) :: spacing(1), origin(1)
    character(len=2048) :: vtype, vunits, vloc

    real(c_double), allocatable :: gx(:), gy(:), gz(:)
    integer(c_int), allocatable :: edge_nodes(:), face_edges(:), face_nodes(:)

    integer(c_int),     allocatable :: wi(:)
    real(c_float),      allocatable :: wf(:)
    real(c_double),     allocatable :: wd(:)
    integer(c_int),     allocatable :: inds(:)
    integer(c_int),     pointer     :: pi(:)  => null()
    real(c_float),      pointer     :: pf(:)  => null()
    real(c_double),     pointer     :: pd(:)  => null()

    do i = 1, nitems
      if (len_trim(names(i)) == 0) cycle

      ! var → grid & metadata
      s = M%get_var_grid(names(i), grid);        call chk('get_var_grid', s)
      s = M%get_var_type(names(i), vtype);       call chk('get_var_type', s)
      s = M%get_var_units(names(i), vunits);     call chk('get_var_units', s)
      s = M%get_var_itemsize(names(i), gsize);   call chk('get_var_itemsize', s)
      s = M%get_var_nbytes(names(i), gsize);     call chk('get_var_nbytes', s)
      s = M%get_var_location(names(i), vloc);    call chk('get_var_location', s)

      write(*,'(A,1X,A)')   'var:', trim(names(i))
      write(*,'(A,I0)')     '  grid   =', grid
      write(*,'(A,1X,A)')   '  type   =', trim(vtype)
      write(*,'(A,1X,A)')   '  units  =', trim(vunits)
      write(*,'(A,1X,A)')   '  loc    =', trim(vloc)
      write(*,'(A,I0)')     '  nbytes =', gsize

      ! grid basics
      s = M%get_grid_rank(grid, rank);        call chk('get_grid_rank', s)
      s = M%get_grid_size(grid, gsize);       call chk('get_grid_size', s)
      s = M%get_grid_shape(grid, shape);      call chk('get_grid_shape', s)
      s = M%get_grid_spacing(grid, spacing);  call chk('get_grid_spacing', s)
      s = M%get_grid_origin(grid, origin);    call chk('get_grid_origin', s)
      write(*,'(A,I0)')     '  grid rank =', rank
      write(*,'(A,I0)')     '  grid size =', gsize
      write(*,'(A,I0)')     '  shape(1)  =', shape(1)
      write(*,'(A,F8.3)')   '  spacing(1)=', spacing(1)
      write(*,'(A,F8.3)')   '  origin(1) =', origin(1)

      ! coords
      allocate(gx(max(1,shape(1))))
      allocate(gy(max(1,shape(1))))
      allocate(gz(max(1,shape(1))))
      s = M%get_grid_x(grid, gx); call chk('get_grid_x', s)
      s = M%get_grid_y(grid, gy); call chk('get_grid_y', s)
      s = M%get_grid_z(grid, gz); call chk('get_grid_z', s)
      deallocate(gx, gy, gz)

      ! topology counts
      call maybe(M%get_grid_node_count(grid, nnode), 'get_grid_node_count')
      call maybe(M%get_grid_edge_count(grid, nedge), 'get_grid_edge_count')
      call maybe(M%get_grid_face_count(grid, nface), 'get_grid_face_count')
      write(*,'(A,3I0)') '  counts (node,edge,face)=', nnode, nedge, nface

      ! connectivity — YOUR INTERFACES EXPECT RANK-1 ARRAYS of length 2*n
      if (nedge > 0) then
        allocate(edge_nodes(2*nedge))
        s = M%get_grid_edge_nodes(grid, edge_nodes); call chk('get_grid_edge_nodes', s)
        deallocate(edge_nodes)
      end if
      if (nface > 0) then
        allocate(face_edges(2*nface))
        s = M%get_grid_face_edges(grid, face_edges); call chk('get_grid_face_edges', s)
        deallocate(face_edges)

        allocate(face_nodes(2*nface))
        s = M%get_grid_face_nodes(grid, face_nodes); call chk('get_grid_face_nodes', s)
        deallocate(face_nodes)

        s = M%get_grid_nodes_per_face(grid, nodes_per_face); call chk('get_grid_nodes_per_face', s)
      end if

      ! Value GETs
      allocate(wi(8), wf(8), wd(8))
      wi = 0; wf = 0.0_c_float; wd = 0.0d0
      call maybe(M%get_value_int   (names(i), wi), 'get_value_int')
      call maybe(M%get_value_float (names(i), wf), 'get_value_float')
      call maybe(M%get_value_double(names(i), wd), 'get_value_double')

      ! Pointer GETs
      nullify(pi, pf, pd)
      call maybe(M%get_value_ptr_int   (names(i), pi), 'get_value_ptr_int')
      call maybe(M%get_value_ptr_float (names(i), pf), 'get_value_ptr_float')
      call maybe(M%get_value_ptr_double(names(i), pd), 'get_value_ptr_double')
      if (associated(pi)) nullify(pi)
      if (associated(pf)) nullify(pf)
      if (associated(pd)) nullify(pd)

      ! SET full arrays
      do j=1,size(wi); wi(j)=j; end do
      do j=1,size(wf); wf(j)=real(j,c_float); end do
      do j=1,size(wd); wd(j)=real(j,c_double); end do
      call maybe(M%set_value_int   (names(i), wi), 'set_value_int')
      call maybe(M%set_value_float (names(i), wf), 'set_value_float')
      call maybe(M%set_value_double(names(i), wd), 'set_value_double')

      ! SET at indices
      allocate(inds(min(3, max(1, size(wd)))))
      do j = 1, size(inds); inds(j) = j; end do
      call maybe(M%set_value_at_indices_int   (names(i), inds, wi(1:size(inds))), 'set_value_at_indices_int')
      call maybe(M%set_value_at_indices_float (names(i), inds, wf(1:size(inds))), 'set_value_at_indices_float')
      call maybe(M%set_value_at_indices_double(names(i), inds, wd(1:size(inds))), 'set_value_at_indices_double')
      deallocate(wi, wf, wd, inds)

    end do
  end subroutine test_var_list

end program test_sfincs_bmi_full

