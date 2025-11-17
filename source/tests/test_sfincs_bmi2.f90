program test_sfincs_bmi2
  use, intrinsic :: iso_fortran_env, only: real32, real64
  use bmif_2_0, only: BMI_SUCCESS, BMI_FAILURE
  use sfincs_bmi2, only: sfincs_bmi
  implicit none

  logical, parameter :: VERBOSE = .true.
  integer :: n_pass = 0, n_fail = 0

  ! Run then summarize
  call run_all_tests()
  call summary_and_stop()

contains

  ! ---------- tiny test utils ----------
  subroutine check(cond, msg)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: msg
    if (cond) then
      n_pass = n_pass + 1
      if (VERBOSE) then
        write(*,'("PASS: ",A)') trim(msg)
      end if
    else
      n_fail = n_fail + 1
      write(*,'("FAIL: ",A)') trim(msg)
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
    s = adjustl(s)
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
    type(sfincs_bmi) :: m, m_uninit
    integer :: s, nx, ny, n, i
    integer :: count, grid, rank, sizei, shp(2)
    double precision :: spacing(2), origin(2)
    double precision, allocatable :: x(:), y(:), z(:)
    character(len=:), pointer :: comp_name
    character(len=:), pointer :: in_names(:), out_names(:)
    character(len=32) :: time_units, var_units, var_type, var_loc, grid_type
    double precision :: t0, t1, tcur, dt, t_before, t_after
    integer :: nbytes, itemsize
    real(real32), allocatable :: fbuf(:), fbuf2(:), fsmall(:)
    real(real64), allocatable :: dbuf(:), dbuf2(:)
    real(real32), pointer     :: fptr(:) => null()
    real(real64), pointer     :: dptr(:) => null()
    integer,      allocatable :: inds(:)
    integer,      pointer     :: iptr(:) => null()   ! for get_value_ptr_int (pointer dummy)
    integer :: tmpi

    ! ---------- update BEFORE initialize should fail ----------
    s = m_uninit%update()
    call check_status_is(s, BMI_FAILURE, 'update before initialize should fail')

    ! ---------- initialize ----------
    s = m%initialize('')   ! config path handled inside BMI; empty is OK for dummy or real wiring
    call check_status_is(s, BMI_SUCCESS, 'initialize')

    nx = m%nx
    ny = m%ny
    n  = nx*ny
    call check(nx > 0 .and. ny > 0, 'grid dims nx,ny are positive')
    call check(m%dx > 0.d0 .and. m%dy > 0.d0, 'spacing dx,dy > 0')

    ! ---------- component & names ----------
    s = m%get_component_name(comp_name)
    call check_status_is(s, BMI_SUCCESS, 'get_component_name')
    call check(associated(comp_name), 'component_name associated')
    call check(len(comp_name) >= 1, 'component_name length > 0')

    s = m%get_input_item_count(count)
    call check_status_is(s, BMI_SUCCESS, 'get_input_item_count')
    call check(count == 0, 'input_item_count == 0')

    s = m%get_output_item_count(count)
    call check_status_is(s, BMI_SUCCESS, 'get_output_item_count')
    call check(count == 3, 'output_item_count == 3')

    s = m%get_input_var_names(in_names)
    call check_status_is(s, BMI_SUCCESS, 'get_input_var_names')
    call check(associated(in_names), 'input_names associated')
    call check(size(in_names) == 0, 'input_names size == 0')

    s = m%get_output_var_names(out_names)
    call check_status_is(s, BMI_SUCCESS, 'get_output_var_names')
    call check(associated(out_names), 'output_names associated')
    call check(size(out_names) == 3, 'output_names size == 3')
    call check(trim(out_names(1))=='zs' .and. trim(out_names(2))=='zb' .and. trim(out_names(3))=='depth', &
               'output_names are zs,zb,depth')

    ! ---------- time ----------
    s = m%get_start_time(t0);       call check_status_is(s, BMI_SUCCESS, 'get_start_time')
    s = m%get_end_time(t1);         call check_status_is(s, BMI_SUCCESS, 'get_end_time')
    s = m%get_current_time(tcur);   call check_status_is(s, BMI_SUCCESS, 'get_current_time')
    s = m%get_time_step(dt);        call check_status_is(s, BMI_SUCCESS, 'get_time_step')
    s = m%get_time_units(time_units); call check_status_is(s, BMI_SUCCESS, 'get_time_units')

    call check(dt > 0.d0,      'dt > 0')
    call check(t1 >= t0,       't_end >= t_start')
    call check(tcur >= t0 .and. tcur <= t1 + 1.d-8, 't_cur in [t_start, t_end + eps]')
    call check(trim(time_units) == 's', 'time_units = s')

    ! ---------- var metadata ----------
    s = m%get_var_grid('zs', grid)
    call check_status_is(s, BMI_SUCCESS, 'get_var_grid zs')
    call check(grid == 1, 'var_grid zs == 1')

    s = m%get_var_grid('zb', grid);    call check_status_is(s, BMI_SUCCESS, 'get_var_grid zb')
    s = m%get_var_grid('depth', grid); call check_status_is(s, BMI_SUCCESS, 'get_var_grid depth')
    s = m%get_var_grid('unknown', grid); call check_status_is(s, BMI_FAILURE, 'get_var_grid unknown')

    s = m%get_var_type('zs', var_type)
    call check_status_is(s, BMI_SUCCESS, 'get_var_type zs')
    call check(trim(var_type) == 'real', 'var_type=real')

    s = m%get_var_units('zs', var_units)
    call check_status_is(s, BMI_SUCCESS, 'get_var_units zs')
    call check(trim(var_units) == 'm', 'var_units=m')

    s = m%get_var_itemsize('zs', itemsize)
    call check_status_is(s, BMI_SUCCESS, 'get_var_itemsize zs')
    call check(itemsize == 4, 'itemsize=4')

    s = m%get_var_nbytes('zs', nbytes)
    call check_status_is(s, BMI_SUCCESS, 'get_var_nbytes zs')
    call check(nbytes == 4*n, 'nbytes=4*n')

    s = m%get_var_location('zs', var_loc)
    call check_status_is(s, BMI_SUCCESS, 'get_var_location zs')
    call check(trim(var_loc) == 'node', 'var_loc=node')

    ! Unknown variable failure paths
    s = m%get_var_type('unknown', var_type)
    call check_status_is(s, BMI_FAILURE, 'get_var_type unknown')
    s = m%get_var_units('unknown', var_units)
    call check_status_is(s, BMI_FAILURE, 'get_var_units unknown')
    s = m%get_var_itemsize('unknown', itemsize)
    call check_status_is(s, BMI_FAILURE, 'get_var_itemsize unknown')
    s = m%get_var_nbytes('unknown', nbytes)
    call check(nbytes == 0, 'get_var_nbytes unknown -> 0')
    s = m%get_var_location('unknown', var_loc)
    call check_status_is(s, BMI_FAILURE, 'get_var_location unknown')

    ! ---------- grid metadata ----------
    s = m%get_grid_rank(1, rank)
    call check_status_is(s, BMI_SUCCESS, 'get_grid_rank')
    call check(rank == 2, 'rank=2')

    s = m%get_grid_size(1, sizei)
    call check_status_is(s, BMI_SUCCESS, 'get_grid_size')
    call check(sizei == n, 'size=n')

    s = m%get_grid_type(1, grid_type)
    call check_status_is(s, BMI_SUCCESS, 'get_grid_type')
    call check(trim(grid_type) == 'uniform_rectilinear', 'grid_type=uniform_rectilinear')

    shp = -1
    s = m%get_grid_shape(1, shp)
    call check_status_is(s, BMI_SUCCESS, 'get_grid_shape')
    call check(all(shp == [ny, nx]), 'shape=(ny,nx)')

    spacing = -1.d0
    s = m%get_grid_spacing(1, spacing)
    call check_status_is(s, BMI_SUCCESS, 'get_grid_spacing')
    call check(rel_eq(spacing(1), m%dy, 1d-12) .and. rel_eq(spacing(2), m%dx, 1d-12), 'spacing=(dy,dx)')

    origin = -1.d0
    s = m%get_grid_origin(1, origin)
    call check_status_is(s, BMI_SUCCESS, 'get_grid_origin')
    call check(rel_eq(origin(1), m%y0, 1d-12) .and. rel_eq(origin(2), m%x0, 1d-12), 'origin=(y0,x0)')

    allocate(x(n), y(n), z(n))
    x = -999.d0; y = -999.d0; z = -999.d0
    s = m%get_grid_x(1, x); call check_status_is(s, BMI_SUCCESS, 'get_grid_x')
    s = m%get_grid_y(1, y); call check_status_is(s, BMI_SUCCESS, 'get_grid_y')
    s = m%get_grid_z(1, z); call check_status_is(s, BMI_SUCCESS, 'get_grid_z')

    call check(size(x) == n .and. size(y) == n .and. size(z) == n, 'grid coordinate sizes = n')
    call check(all(z == 0.d0), 'z all zeros (2D grid)')

    ! Basic sanity on x,y monotonicity (if grid is at least 2x2)
    if (nx > 1) then
      call check(rel_eq(x(2) - x(1), m%dx, 1d-6), 'x spacing ~ dx for first row')
    end if
    if (ny > 1) then
      call check(rel_eq(y(nx+1) - y(1), m%dy, 1d-6), 'y spacing ~ dy across rows')
    end if

    s = m%get_grid_rank(2, rank)
    call check_status_is(s, BMI_FAILURE, 'get_grid_rank wrong grid')

    ! ---------- unstructured stubs ----------
    s = m%get_grid_edge_count(1, tmpi);       call check_status_is(s, BMI_FAILURE, 'get_grid_edge_count (stub)')
    s = m%get_grid_face_count(1, tmpi);       call check_status_is(s, BMI_FAILURE, 'get_grid_face_count (stub)')
    s = m%get_grid_node_count(1, tmpi);       call check_status_is(s, BMI_FAILURE, 'get_grid_node_count (stub)')
    allocate(inds(4)); inds = 0
    s = m%get_grid_edge_nodes(1, inds);       call check_status_is(s, BMI_FAILURE, 'get_grid_edge_nodes (stub)')
    s = m%get_grid_face_nodes(1, inds);       call check_status_is(s, BMI_FAILURE, 'get_grid_face_nodes (stub)')
    s = m%get_grid_nodes_per_face(1, inds);   call check_status_is(s, BMI_FAILURE, 'get_grid_nodes_per_face (stub)')
    s = m%get_grid_face_edges(1, inds);       call check_status_is(s, BMI_FAILURE, 'get_grid_face_edges (stub)')
    deallocate(inds)

    ! ---------- values: FLOAT ----------
    allocate(fbuf(n), fbuf2(n))
    fbuf = -123.0_real32
    s = m%get_value_float('zs', fbuf)
    call check_status_is(s, BMI_SUCCESS, 'get_value_float zs (initial)')
    ! no assumption about value, just test round-trip consistency

    s = m%get_value_ptr_float('zs', fptr)
    call check_status_is(s, BMI_SUCCESS, 'get_value_ptr_float zs')
    call check(associated(fptr), 'fptr associated')
    call check(size(fptr) == n, 'fptr size')

    ! Too-small buffer for get_value_float should fail
    allocate(fsmall(max(1, n-1)))
    s = m%get_value_float('zs', fsmall)
    call check_status_is(s, BMI_FAILURE, 'get_value_float with too-small dest should fail')
    deallocate(fsmall)

    ! set_value_float then get_value_float
    fbuf = 2.0_real32
    s = m%set_value_float('zs', fbuf)
    call check_status_is(s, BMI_SUCCESS, 'set_value_float zs -> 2.0')

    fbuf2 = -7.0_real32
    s = m%get_value_float('zs', fbuf2)
    call check_status_is(s, BMI_SUCCESS, 'get_value_float zs after set')
    call check(all(fbuf2 == 2.0_real32), 'zs now 2.0')

    ! Set zb via set_value_at_indices_float -> 1.0
    allocate(inds(n))
    do i = 1, n
      inds(i) = i
    end do
    fbuf = 1.0_real32
    s = m%set_value_at_indices_float('zb', inds, fbuf)
    call check_status_is(s, BMI_SUCCESS, 'set_value_at_indices_float zb -> 1.0')

    fbuf2 = -9.0_real32
    s = m%get_value_float('zb', fbuf2)
    call check_status_is(s, BMI_SUCCESS, 'get_value_float zb after set_at_indices')
    call check(all(fbuf2 == 1.0_real32), 'zb now 1.0')

    ! depth should now be zs - zb = 1.0 (or clipped at >=0)
    fbuf2 = -5.0_real32
    s = m%get_value_float('depth', fbuf2)
    call check_status_is(s, BMI_SUCCESS, 'get_value_float depth after zs/zb set')
    write(*,'(A,1X,F20.16,1X,A,1X,F20.16)') 'TEST (float) depth ~ 1.0: got=', real(fbuf2(1),real64), &
                                            'expected=', 1.0_real64
    call check(all(abs(real(fbuf2,real64) - 1.0_real64) < 1e-6_real64), 'depth ~ 1.0 (float)')

    ! Subset get_value_at_indices_float
    if (n >= 10) then
      fbuf2 = -5.0_real32
      inds(1:10) = [( (i), i=1,10 )]
      s = m%get_value_at_indices_float('zs', fbuf2(1:10), inds(1:10))
      call check_status_is(s, BMI_SUCCESS, 'get_value_at_indices_float zs subset')
      call check(all(fbuf2(1:10) == 2.0_real32), 'subset zs values are 2.0')
    end if

    s = m%get_value_float('unknown', fbuf)
    call check_status_is(s, BMI_FAILURE, 'get_value_float unknown')
    s = m%get_value_ptr_float('unknown', fptr)
    call check_status_is(s, BMI_FAILURE, 'get_value_ptr_float unknown')

    deallocate(fbuf, fbuf2, inds)

    ! ---------- values: DOUBLE ----------
    allocate(dbuf(n), dbuf2(n))
    dbuf = -123.0_real64
    s = m%get_value_double('zs', dbuf)
    call check_status_is(s, BMI_SUCCESS, 'get_value_double zs')
    ! Should mirror zs = 2.0
    call check(all(abs(dbuf - 2.0_real64) < 1e-12_real64), 'zs double mirror 2.0')

    ! --- depth double after zs/zb set: expect 1.0 ---
    dbuf = -9.0_real64
    s = m%get_value_double('depth', dbuf)
    call check_status_is(s, BMI_SUCCESS, 'get_value_double depth')
    write(*,'(A,1X,F20.16,1X,A,1X,F20.16)') 'TEST depth double == 1.0: got=', dbuf(1), &
                                            'expected=', 1.0_real64
    write(*,'(A,1X,F20.16,1X,F20.16)') '  min/max depth:', minval(dbuf), maxval(dbuf)
    call check(all(abs(dbuf - 1.0_real64) < 1e-10_real64), 'depth double == 1.0')

    s = m%get_value_ptr_double('zb', dptr)
    call check_status_is(s, BMI_SUCCESS, 'get_value_ptr_double zb')
    call check(associated(dptr), 'dptr associated')
    call check(size(dptr) == n, 'dptr size')

    ! set_value_double then get_value_float / depth consistency
    dbuf2 = 3.5_real64
    s = m%set_value_double('zs', dbuf2)
    call check_status_is(s, BMI_SUCCESS, 'set_value_double zs -> 3.5')

    allocate(fbuf(n))
    s = m%get_value_float('zs', fbuf)
    call check_status_is(s, BMI_SUCCESS, 'get_value_float zs after set_double')
    call check(all(abs(real(fbuf,real64) - 3.5_real64) < 1e-6_real64), 'zs float ~ 3.5 after double set')

    ! --- depth double after zs double-set: expect 2.5 ---
    s = m%get_value_double('depth', dbuf)
    call check_status_is(s, BMI_SUCCESS, 'get_value_double depth after set_double')
    write(*,'(A,1X,F20.16,1X,A,1X,F20.16)') 'TEST depth double == 2.5: got=', dbuf(1), &
                                            'expected=', 2.5_real64
    write(*,'(A,1X,F20.16,1X,F20.16)') '  min/max depth:', minval(dbuf), maxval(dbuf)
    ! zb was 1.0, so depth ~ 2.5
    call check(all(abs(dbuf - 2.5_real64) < 1e-12_real64), 'depth double == 2.5')

    ! set depth via set_value_at_indices_double
    do i = 1, n
      dbuf(i) = 7.0_real64
    end do
    allocate(inds(n))
    do i = 1, n
      inds(i) = i
    end do
    s = m%set_value_at_indices_double('depth', inds, dbuf)
    call check_status_is(s, BMI_SUCCESS, 'set_value_at_indices_double depth -> 7.0')
    deallocate(inds)

    ! --- depth double after explicit depth set: expect 7.0 ---
    s = m%get_value_double('depth', dbuf2)
    call check_status_is(s, BMI_SUCCESS, 'get_value_double depth after set_at_indices_double')
    write(*,'(A,1X,F20.16,1X,A,1X,F20.16)') 'TEST depth double == 7.0: got=', dbuf2(1), &
                                            'expected=', 7.0_real64
    write(*,'(A,1X,F20.16,1X,F20.16)') '  min/max depth:', minval(dbuf2), maxval(dbuf2)
    call check(all(abs(dbuf2 - 7.0_real64) < 1e-12_real64), 'depth double == 7.0')

    deallocate(fbuf, dbuf, dbuf2)

    ! ---------- integer suite (stubs -> failure) ----------
    allocate(inds(5)); inds = [1,2,3,4,5]
    s = m%get_value_int('zs', inds)
    call check_status_is(s, BMI_FAILURE, 'get_value_int (stub)')
    s = m%set_value_int('zs', inds)
    call check_status_is(s, BMI_FAILURE, 'set_value_int (stub)')
    nullify(iptr)
    s = m%get_value_ptr_int('zs', iptr)
    call check_status_is(s, BMI_FAILURE, 'get_value_ptr_int (stub)')
    call check(.not. associated(iptr), 'iptr not associated (stub failure path)')
    s = m%get_value_at_indices_int('zs', inds, inds)
    call check_status_is(s, BMI_FAILURE, 'get_value_at_indices_int (stub)')
    s = m%set_value_at_indices_int('zs', inds, inds)
    call check_status_is(s, BMI_FAILURE, 'set_value_at_indices_int (stub)')
    deallocate(inds)

    ! ---------- update / update_until ----------
    s = m%get_current_time(t_before)
    call check_status_is(s, BMI_SUCCESS, 'current time before update')

    s = m%update()
    call check_status_is(s, BMI_SUCCESS, 'update')

    s = m%get_current_time(t_after)
    call check_status_is(s, BMI_SUCCESS, 'current time after update')

    ! SFINCS uses adaptive time-stepping, so we only require that time advances forward,
    ! not that it equals t_before + dt exactly.
    call check(t_after > t_before, 'time advanced forward (variable-step model)')

    ! New test: time should still be within [t_start, t_end + eps]
    call check(t_after >= t0 .and. t_after <= t1 + 1d-8, &
               'time after update within [t_start, t_end + eps]')

    ! update_until forward in time
    s = m%update_until(t_after + 5.d0)
    call check_status_is(s, BMI_SUCCESS, 'update_until(t+5)')
    s = m%get_current_time(tcur)
    call check_status_is(s, BMI_SUCCESS, 'current time after update_until forward')
    call check(tcur >= t_after + 5.d0 - 1d-8, 'time >= target after update_until forward')

    ! update_until earlier time should be a no-op
    s = m%update_until(tcur - 2.d0)
    call check_status_is(s, BMI_SUCCESS, 'update_until earlier time (no-op)')
    s = m%get_current_time(t_after)
    call check_status_is(s, BMI_SUCCESS, 'time after no-op update_until')
    call check(rel_eq(t_after, tcur, 1d-12), 'time unchanged after no-op')

    ! ---------- finalize ----------
    s = m%finalize()
    call check_status_is(s, BMI_SUCCESS, 'finalize')

    s = m%update()
    call check_status_is(s, BMI_FAILURE, 'update after finalize should fail')
  end subroutine run_all_tests

end program test_sfincs_bmi2

