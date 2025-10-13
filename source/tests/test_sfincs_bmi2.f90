program test_sfincs_bmi2
  use, intrinsic :: iso_fortran_env, only: real32, real64
  use bmif_2_0, only: BMI_SUCCESS, BMI_FAILURE
  use sfincs_bmi2, only: sfincs_bmi
  implicit none

  integer :: n_pass = 0, n_fail = 0

  ! Run then summarize (must be before CONTAINS)
  call run_all_tests()
  call summary_and_stop()

contains
  ! ---------- tiny test utils ----------
  subroutine check(cond, msg)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: msg
    if (cond) then
      n_pass = n_pass + 1
    else
      n_fail = n_fail + 1
      write(*,'(A)') 'FAIL: '//trim(msg)
    end if
  end subroutine check

  subroutine check_status_is(status, expect, msg)
    integer, intent(in) :: status, expect
    character(len=*), intent(in) :: msg
    call check(status == expect, trim(msg)//' [status='//trim(int_to_str(status))//']')
  end subroutine check_status_is

  function int_to_str(i) result(s)
    integer, intent(in) :: i
    character(len=32) :: s
    write(s,'(I0)') i
  end function int_to_str

  function rel_eq(a, b, tol) result(ok)
    real(real64), intent(in) :: a, b, tol
    logical :: ok
    ok = (abs(a-b) <= tol*(1.0d0 + max(abs(a),abs(b))))
  end function rel_eq

  subroutine summary_and_stop()
    write(*,'(A,I0,A,I0)') 'Tests: pass=', n_pass, '  fail=', n_fail
    if (n_fail == 0) then
      stop 0
    else
      stop 1
    end if
  end subroutine summary_and_stop

  ! =========================================================
  ! ================   MAIN TEST SEQUENCE   =================
  ! =========================================================
  subroutine run_all_tests()
    type(sfincs_bmi) :: m
    integer :: s, nx, ny, n, i
    integer :: count, grid, rank, sizei, shp(2)
    double precision :: spacing(2), origin(2)
    double precision, allocatable :: x(:), y(:), z(:)
    character(len=:), pointer :: comp_name
    character(len=:), pointer :: in_names(:), out_names(:)
    character(len=32) :: time_units, var_units, var_type, var_loc, grid_type
    double precision :: t0, t1, tcur, dt
    integer :: nbytes, itemsize
    real(real32), allocatable :: fbuf(:), fbuf2(:)
    real(real64), allocatable :: dbuf(:), dbuf2(:)
    real(real32), pointer     :: fptr(:) => null()
    real(real64), pointer     :: dptr(:) => null()
    integer,      allocatable :: inds(:)
    integer,      pointer     :: iptr(:) => null()   ! for get_value_ptr_int (pointer dummy)
    integer :: tmpi

    ! ---------- initialize ----------
    s = m%initialize(''); call check_status_is(s, BMI_SUCCESS, 'initialize')

    nx = 10; ny = 10; n = nx*ny
    call check(m%nx==nx .and. m%ny==ny, 'grid dims nx,ny')
    call check(rel_eq(m%dx,100.d0,1d-12) .and. rel_eq(m%dy,100.d0,1d-12), 'spacing dx,dy')
    call check(rel_eq(m%x0,0.d0,1d-12) .and. rel_eq(m%y0,0.d0,1d-12), 'origin x0,y0')

    ! ---------- component & names ----------
    s = m%get_component_name(comp_name); call check_status_is(s,BMI_SUCCESS,'get_component_name')
    call check(associated(comp_name), 'component_name associated')
    call check(len(comp_name) >= 1, 'component_name length > 0')

    s = m%get_input_item_count(count); call check_status_is(s,BMI_SUCCESS,'get_input_item_count')
    call check(count == 0, 'input_item_count == 0')

    s = m%get_output_item_count(count); call check_status_is(s,BMI_SUCCESS,'get_output_item_count')
    call check(count == 3, 'output_item_count == 3')

    s = m%get_input_var_names(in_names);  call check_status_is(s,BMI_SUCCESS,'get_input_var_names')
    call check(associated(in_names), 'input_names associated')
    call check(size(in_names) == 0, 'input_names size == 0')

    s = m%get_output_var_names(out_names); call check_status_is(s,BMI_SUCCESS,'get_output_var_names')
    call check(associated(out_names), 'output_names associated')
    call check(size(out_names) == 3, 'output_names size == 3')
    call check(trim(out_names(1))=='zs' .and. trim(out_names(2))=='zb' .and. trim(out_names(3))=='depth', 'output_names are zs,zb,depth')

    ! ---------- time ----------
    s = m%get_start_time(t0);       call check_status_is(s,BMI_SUCCESS,'get_start_time')
    s = m%get_end_time(t1);         call check_status_is(s,BMI_SUCCESS,'get_end_time')
    s = m%get_current_time(tcur);   call check_status_is(s,BMI_SUCCESS,'get_current_time')
    s = m%get_time_step(dt);        call check_status_is(s,BMI_SUCCESS,'get_time_step')
    s = m%get_time_units(time_units); call check_status_is(s,BMI_SUCCESS,'get_time_units')
    call check(rel_eq(t0,0.d0,1d-12), 't_start = 0')
    call check(rel_eq(t1,3600.d0,1d-12), 't_end = 3600')
    call check(rel_eq(tcur,t0,1d-12), 't_cur = t_start')
    call check(rel_eq(dt,1.d0,1d-12), 'dt = 1')
    call check(trim(time_units) == 's', 'time_units = s')

    ! ---------- var metadata ----------
    s = m%get_var_grid('zs', grid);          call check_status_is(s,BMI_SUCCESS,'get_var_grid zs'); call check(grid==1,'var_grid zs == 1')
    s = m%get_var_grid('zb', grid);          call check_status_is(s,BMI_SUCCESS,'get_var_grid zb')
    s = m%get_var_grid('depth', grid);       call check_status_is(s,BMI_SUCCESS,'get_var_grid depth')
    s = m%get_var_grid('unknown', grid);     call check_status_is(s,BMI_FAILURE,'get_var_grid unknown')

    s = m%get_var_type('zs', var_type);      call check_status_is(s,BMI_SUCCESS,'get_var_type zs');   call check(trim(var_type)=='real','var_type=real')
    s = m%get_var_units('zs', var_units);    call check_status_is(s,BMI_SUCCESS,'get_var_units zs');  call check(trim(var_units)=='m','var_units=m')
    s = m%get_var_itemsize('zs', itemsize);  call check_status_is(s,BMI_SUCCESS,'get_var_itemsize zs'); call check(itemsize==4, 'itemsize=4')
    s = m%get_var_nbytes('zs', nbytes);      call check_status_is(s,BMI_SUCCESS,'get_var_nbytes zs'); call check(nbytes==4*n,'nbytes=4*n')
    s = m%get_var_location('zs', var_loc);   call check_status_is(s,BMI_SUCCESS,'get_var_location zs'); call check(trim(var_loc)=='node','var_loc=node')

    s = m%get_var_type('unknown', var_type);   call check_status_is(s,BMI_FAILURE,'get_var_type unknown')
    s = m%get_var_units('unknown', var_units); call check_status_is(s,BMI_FAILURE,'get_var_units unknown')
    s = m%get_var_itemsize('unknown', itemsize); call check_status_is(s,BMI_FAILURE,'get_var_itemsize unknown')
    s = m%get_var_nbytes('unknown', nbytes);  call check(nbytes==0,'get_var_nbytes unknown -> 0')
    s = m%get_var_location('unknown', var_loc); call check_status_is(s,BMI_FAILURE,'get_var_location unknown')

    ! ---------- grid metadata ----------
    s = m%get_grid_rank(1, rank);     call check_status_is(s,BMI_SUCCESS,'get_grid_rank'); call check(rank==2,'rank=2')
    s = m%get_grid_size(1, sizei);    call check_status_is(s,BMI_SUCCESS,'get_grid_size'); call check(sizei==n,'size=n')
    s = m%get_grid_type(1, grid_type);call check_status_is(s,BMI_SUCCESS,'get_grid_type'); call check(trim(grid_type)=='uniform_rectilinear','grid_type=uniform_rectilinear')

    shp = -1
    s = m%get_grid_shape(1, shp);     call check_status_is(s,BMI_SUCCESS,'get_grid_shape'); call check(all(shp==[ny,nx]),'shape=(ny,nx)')

    spacing = -1.d0
    s = m%get_grid_spacing(1, spacing); call check_status_is(s,BMI_SUCCESS,'get_grid_spacing')
    call check(rel_eq(spacing(1),100.d0,1d-12) .and. rel_eq(spacing(2),100.d0,1d-12), 'spacing=(dy,dx)')

    origin = -1.d0
    s = m%get_grid_origin(1, origin);   call check_status_is(s,BMI_SUCCESS,'get_grid_origin')
    call check(rel_eq(origin(1),0.d0,1d-12) .and. rel_eq(origin(2),0.d0,1d-12), 'origin=(y0,x0)')

    allocate(x(n), y(n), z(n))
    x = -999; y = -999; z = -999
    s = m%get_grid_x(1, x);             call check_status_is(s,BMI_SUCCESS,'get_grid_x')
    s = m%get_grid_y(1, y);             call check_status_is(s,BMI_SUCCESS,'get_grid_y')
    s = m%get_grid_z(1, z);             call check_status_is(s,BMI_SUCCESS,'get_grid_z')
    call check(rel_eq(x(1),  50.d0,1d-12),'x(1)=50')
    call check(rel_eq(y(1),  50.d0,1d-12),'y(1)=50')
    call check(all(z == 0.d0), 'z all zeros')

    s = m%get_grid_rank(2, rank);       call check_status_is(s,BMI_FAILURE,'get_grid_rank wrong grid')

    ! ---------- unstructured stubs ----------
    s = m%get_grid_edge_count(1, tmpi);       call check_status_is(s,BMI_FAILURE,'get_grid_edge_count (stub)')
    s = m%get_grid_face_count(1, tmpi);       call check_status_is(s,BMI_FAILURE,'get_grid_face_count (stub)')
    s = m%get_grid_node_count(1, tmpi);       call check_status_is(s,BMI_FAILURE,'get_grid_node_count (stub)')
    allocate(inds(4)); inds = 0
    s = m%get_grid_edge_nodes(1, inds);       call check_status_is(s,BMI_FAILURE,'get_grid_edge_nodes (stub)')
    s = m%get_grid_face_nodes(1, inds);       call check_status_is(s,BMI_FAILURE,'get_grid_face_nodes (stub)')
    s = m%get_grid_nodes_per_face(1, inds);   call check_status_is(s,BMI_FAILURE,'get_grid_nodes_per_face (stub)')
    s = m%get_grid_face_edges(1, inds);       call check_status_is(s,BMI_FAILURE,'get_grid_face_edges (stub)')
    deallocate(inds)

    ! ---------- values: FLOAT ----------
    allocate(fbuf(n), fbuf2(n))
    fbuf = -123.0_real32
    s = m%get_value_float('zs', fbuf);  call check_status_is(s,BMI_SUCCESS,'get_value_float zs (zeros)')
    call check(all(fbuf == 0.0_real32),'zs initially zero')

    s = m%get_value_ptr_float('zs', fptr); call check_status_is(s,BMI_SUCCESS,'get_value_ptr_float zs')
    call check(associated(fptr),'fptr associated')
    call check(size(fptr)==n,'fptr size')

    fbuf = 2.0_real32
    s = m%set_value_float('zs', fbuf);  call check_status_is(s,BMI_SUCCESS,'set_value_float zs')
    fbuf2 = -7.0_real32
    s = m%get_value_float('zs', fbuf2); call check_status_is(s,BMI_SUCCESS,'get_value_float zs after set')
    call check(all(fbuf2 == 2.0_real32),'zs now 2.0')

    allocate(inds(n))
    do i=1,n; inds(i)=i; end do
    fbuf = 1.0_real32
    s = m%set_value_at_indices_float('zb', inds, fbuf); call check_status_is(s,BMI_SUCCESS,'set_value_at_indices_float zb -> 1.0')
    fbuf2 = -9.0_real32
    s = m%get_value_float('zb', fbuf2); call check_status_is(s,BMI_SUCCESS,'get_value_float zb after set_at_indices')
    call check(all(fbuf2 == 1.0_real32),'zb now 1.0')

    fbuf2 = -5.0_real32
    inds = [( (i), i=1,10 )]
    s = m%get_value_at_indices_float('zs', fbuf2(1:10), inds(1:10))
    call check_status_is(s,BMI_SUCCESS,'get_value_at_indices_float zs subset')
    call check(all(fbuf2(1:10) == 2.0_real32), 'subset zs values are 2.0')

    s = m%get_value_float('unknown', fbuf);      call check_status_is(s,BMI_FAILURE,'get_value_float unknown')
    s = m%get_value_ptr_float('unknown', fptr);  call check_status_is(s,BMI_FAILURE,'get_value_ptr_float unknown')
    deallocate(fbuf, fbuf2, inds)

    ! ---------- values: DOUBLE ----------
    allocate(dbuf(n), dbuf2(n))
    dbuf = -123.0_real64
    s = m%get_value_double('zs', dbuf);  call check_status_is(s,BMI_SUCCESS,'get_value_double zs')
    call check(all(dbuf == 2.0_real64),'zs double mirror 2.0')

    dbuf = -9.0_real64
    s = m%get_value_double('depth', dbuf); call check_status_is(s,BMI_SUCCESS,'get_value_double depth')
    call check(all(dbuf == 1.0_real64), 'depth double == 1.0')

    s = m%get_value_ptr_double('zb', dptr); call check_status_is(s,BMI_SUCCESS,'get_value_ptr_double zb')
    call check(associated(dptr),'dptr associated'); call check(size(dptr)==n,'dptr size')

    dbuf2 = 3.5_real64
    s = m%set_value_double('zs', dbuf2); call check_status_is(s,BMI_SUCCESS,'set_value_double zs -> 3.5')

    allocate(fbuf(n))
    s = m%get_value_float('zs', fbuf); call check_status_is(s,BMI_SUCCESS,'get_value_float zs after set_double')
    call check(all(abs(real(fbuf,real64)-3.5_real64) < 1e-6_real64),'zs float ~ 3.5 after double set')
    s = m%get_value_double('depth', dbuf); call check_status_is(s,BMI_SUCCESS,'get_value_double depth after set_double')
    call check(all(abs(dbuf-2.5_real64) < 1e-12_real64),'depth double == 2.5')

    do i=1,n; dbuf(i)=7.0_real64; end do
    allocate(inds(n)); do i=1,n; inds(i)=i; end do
    s = m%set_value_at_indices_double('depth', inds, dbuf); call check_status_is(s,BMI_SUCCESS,'set_value_at_indices_double depth -> 7.0')
    deallocate(inds)

    s = m%get_value_double('depth', dbuf2); call check_status_is(s,BMI_SUCCESS,'get_value_double depth after set_at_indices_double')
    call check(all(abs(dbuf2-7.0_real64) < 1e-12_real64),'depth double == 7.0')

    deallocate(fbuf, dbuf, dbuf2)

    ! ---------- integer suite (stubs -> failure) ----------
    allocate(inds(5)); inds = [1,2,3,4,5]
    s = m%get_value_int('zs', inds);                  call check_status_is(s,BMI_FAILURE,'get_value_int (stub)')
    s = m%set_value_int('zs', inds);                  call check_status_is(s,BMI_FAILURE,'set_value_int (stub)')
    nullify(iptr)
    s = m%get_value_ptr_int('zs', iptr);              call check_status_is(s,BMI_FAILURE,'get_value_ptr_int (stub)')
    call check(.not. associated(iptr), 'iptr not associated (stub failure path)')
    s = m%get_value_at_indices_int('zs', inds, inds); call check_status_is(s,BMI_FAILURE,'get_value_at_indices_int (stub)')
    s = m%set_value_at_indices_int('zs', inds, inds); call check_status_is(s,BMI_FAILURE,'set_value_at_indices_int (stub)')
    deallocate(inds)

    ! ---------- update / update_until ----------
    s = m%get_current_time(tcur); call check_status_is(s,BMI_SUCCESS,'current time before update')
    s = m%update();               call check_status_is(s,BMI_SUCCESS,'update')
    s = m%get_current_time(dt);   call check_status_is(s,BMI_SUCCESS,'current time after update')
    call check(rel_eq(dt, tcur+1.d0, 1d-12),'time advanced by dt=1')

    s = m%update_until(5.d0);     call check_status_is(s,BMI_SUCCESS,'update_until(5)')
    s = m%get_current_time(tcur); call check_status_is(s,BMI_SUCCESS,'current time after update_until')
    call check(rel_eq(tcur,5.d0,1d-12), 'time == 5')

    s = m%update_until(2.d0);     call check_status_is(s,BMI_SUCCESS,'update_until earlier time is no-op')
    s = m%get_current_time(dt);   call check_status_is(s,BMI_SUCCESS,'time after no-op update_until')
    call check(rel_eq(dt,tcur,1d-12),'time unchanged after no-op')

    ! ---------- finalize ----------
    s = m%finalize(); call check_status_is(s,BMI_SUCCESS,'finalize')
    s = m%update();   call check_status_is(s,BMI_FAILURE,'update after finalize should fail')
  end subroutine run_all_tests

end program test_sfincs_bmi2

