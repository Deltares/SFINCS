program test_sfincs_bmi2_real
  use, intrinsic :: iso_fortran_env, only: real32, real64
  use bmif_2_0, only: BMI_SUCCESS, BMI_FAILURE
  use sfincs_bmi2, only: sfincs_bmi
  implicit none

  logical, parameter :: VERBOSE = .true.
  integer :: n_pass = 0, n_fail = 0

  call run_all_real_tests()
  call summary_and_stop()

contains

  subroutine check(cond, msg)
    logical, intent(in) :: cond
    character(len=*), intent(in) :: msg
    if (cond) then
      n_pass = n_pass + 1
      if (VERBOSE) write(*,'("PASS: ",A)') trim(msg)
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
    ok = (abs(a-b) <= tol*(1.0_real64 + max(abs(a),abs(b))))
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
  ! ============ REAL SFINCS RUN THROUGH BMI ================
  ! =========================================================
  subroutine run_all_real_tests()
    type(sfincs_bmi) :: m
    integer :: s, nx, ny, n
    double precision :: t0, t1, tcur, dt
    double precision, allocatable :: zs0(:), zs_end(:), zb(:), depth(:)
    double precision, allocatable :: x(:), y(:)
    integer :: grid, sizei, rank
    character(len=32) :: grid_type
    integer, allocatable :: inds(:)
    double precision :: t_before, t_after, target_t
    double precision, allocatable :: zs_snapshot(:)
    ! *** IMPORTANT: point this to your real input ***
    character(len=*), parameter :: cfg = '/home/mohammed.karim/SFINCS/XMI_EXAMPLE/sfincs.inp'

    ! --- Initialize with real config ---
    s = m%initialize(cfg)
    call check_status_is(s, BMI_SUCCESS, 'real: initialize with sfincs.inp')

    nx = m%nx
    ny = m%ny
    n  = nx*ny
    call check(nx > 0 .and. ny > 0, 'real: nx,ny > 0')
    call check(m%dx > 0.d0 .and. m%dy > 0.d0, 'real: dx,dy > 0')

    ! --- Time info ---
    s = m%get_start_time(t0);     call check_status_is(s, BMI_SUCCESS, 'real: get_start_time')
    s = m%get_end_time(t1);       call check_status_is(s, BMI_SUCCESS, 'real: get_end_time')
    s = m%get_current_time(tcur); call check_status_is(s, BMI_SUCCESS, 'real: get_current_time')
    s = m%get_time_step(dt);      call check_status_is(s, BMI_SUCCESS, 'real: get_time_step')
    call check(t1 > t0, 'real: t_end > t_start')
    call check(dt > 0.d0, 'real: dt > 0')

    ! --- Grid metadata sanity (real run) ---
    s = m%get_var_grid('zs', grid)
    call check_status_is(s, BMI_SUCCESS, 'real: get_var_grid zs')
    s = m%get_grid_size(grid, sizei)
    call check_status_is(s, BMI_SUCCESS, 'real: get_grid_size')
    call check(sizei == n, 'real: get_grid_size == nx*ny')

    s = m%get_grid_rank(grid, rank)
    call check_status_is(s, BMI_SUCCESS, 'real: get_grid_rank')
    call check(rank == 2, 'real: rank=2')

    s = m%get_grid_type(grid, grid_type)
    call check_status_is(s, BMI_SUCCESS, 'real: get_grid_type')

    allocate(x(n), y(n))
    s = m%get_grid_x(grid, x); call check_status_is(s, BMI_SUCCESS, 'real: get_grid_x')
    s = m%get_grid_y(grid, y); call check_status_is(s, BMI_SUCCESS, 'real: get_grid_y')

    call check(any(x /= 0.d0), 'real: x grid not all zero')
    call check(any(y /= 0.d0), 'real: y grid not all zero')

    ! --- Snapshot initial zs/zb/depth ---
    allocate(zs0(n), zb(n), depth(n))
    s = m%get_value_double('zs', zs0);   call check_status_is(s, BMI_SUCCESS, 'real: get_value_double zs initial')
    s = m%get_value_double('zb', zb);    call check_status_is(s, BMI_SUCCESS, 'real: get_value_double zb initial')
    s = m%get_value_double('depth', depth); call check_status_is(s, BMI_SUCCESS, 'real: get_value_double depth initial')

    call check(any(abs(zs0) < 1.0e6_real64), 'real: zs initial not NaN/inf everywhere')
    call check(any(abs(zb) < 1.0e6_real64),  'real: zb initial not NaN/inf everywhere')

    ! --- Run full real simulation with repeated BMI updates ---
    t_before = tcur
    do
      s = m%update()
      call check_status_is(s, BMI_SUCCESS, 'real: update loop step')
      s = m%get_current_time(tcur)
      call check_status_is(s, BMI_SUCCESS, 'real: get_current_time in loop')
      if (tcur >= t1 - 1.0e-8_real64) exit
    end do
    t_after = tcur

    call check(t_after > t_before, 'real: t advanced over entire run')

    allocate(zs_end(n))
    s = m%get_value_double('zs', zs_end)
    call check_status_is(s, BMI_SUCCESS, 'real: get_value_double zs end')

    ! *** FIXED: no d-exponent + kind combo ***
    call check(any(abs(zs_end - zs0) > 1.0e-6_real64), &
               'real: zs changed over simulation (non-trivial run)')

    ! --- Check update_until semantics on real run ---
    target_t = min(t_before + 10.0d0, t1)
    s = m%finalize()
    call check_status_is(s, BMI_SUCCESS, 'real: finalize after full run')

    ! Reinitialize to test update_until from t0
    s = m%initialize(cfg)
    call check_status_is(s, BMI_SUCCESS, 'real: re-initialize for update_until test')
    s = m%get_current_time(tcur)
    call check_status_is(s, BMI_SUCCESS, 'real: get_current_time after re-init')

    s = m%update_until(target_t)
    call check_status_is(s, BMI_SUCCESS, 'real: update_until(target_t)')
    s = m%get_current_time(t_after)
    call check_status_is(s, BMI_SUCCESS, 'real: get_current_time after update_until')
    call check(t_after >= target_t - 1.0e-8_real64, 'real: time >= target_t after update_until')

    ! --- Depth consistency check: depth ~ zs - zb at a snapshot ---
    allocate(zs_snapshot(n))
    s = m%get_value_double('zs', zs_snapshot)
    call check_status_is(s, BMI_SUCCESS, 'real: get_value_double zs snapshot')
    s = m%get_value_double('zb', zb)
    call check_status_is(s, BMI_SUCCESS, 'real: get_value_double zb snapshot')
    s = m%get_value_double('depth', depth)
    call check_status_is(s, BMI_SUCCESS, 'real: get_value_double depth snapshot')

    call check(any(abs(depth - (zs_snapshot - zb)) < 1.0e-2_real64), &
               'real: depth ~ zs - zb in at least some cells')

    ! --- Final finalize ---
    s = m%finalize()
    call check_status_is(s, BMI_SUCCESS, 'real: final finalize')
  end subroutine run_all_real_tests

end program test_sfincs_bmi2_real

