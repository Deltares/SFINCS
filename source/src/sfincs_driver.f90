program sfincs_driver
  use bmif_2_0
  use sfincs_bmi2
  implicit none

  !--- Locals
  type(sfincs_bmi)                               :: M
  integer                                        :: s, argc
  character(len=512)                             :: arg
  character(len=512)                             :: config_file

  double precision                               :: t0, t, t_end, dt
  character(len=32)                              :: units
  character(len=BMI_MAX_COMPONENT_NAME), pointer :: cname => null()

  !--- Parse CLI
  argc = command_argument_count()
  if (argc < 1) then
    print *, 'Usage: sfincs_driver [-v] <config_file>'
    stop 1
  end if

  call get_command_argument(1, arg)
  if (trim(arg) == '-v') then
    ! Print component name (version banner) and exit, SCHISM-style
    s = M%initialize('')             ! if your initialize requires a file, you can skip this
    if (s == BMI_SUCCESS) then
      s = M%get_component_name(cname)
      if (s == BMI_SUCCESS .and. associated(cname)) then
        print *, trim(cname)
      else
        print *, 'SFINCS (BMI)';  ! fallback
      end if
      call M%finalize()
    else
      print *, 'SFINCS (BMI)'
    end if
    stop 0
  else
    config_file = trim(arg)
  end if

  !--- Initialize model
  print *, 'Initializing SFINCS (BMI) with config: ', trim(config_file)
  s = M%initialize(config_file);      call chk('initialize', s)

  !--- Time metadata
  s = M%get_start_time(t0);           call chk('get_start_time', s)
  s = M%get_current_time(t);          call chk('get_current_time', s)
  s = M%get_end_time(t_end);          call chk('get_end_time', s)
  s = M%get_time_step(dt);            call chk('get_time_step', s)
  s = M%get_time_units(units);        call chk('get_time_units', s)

  print '(A,1PE13.5,2A)', ' start_time = ', t0, ' [', trim(units), ']'
  print '(A,1PE13.5,2A)', '   end_time = ', t_end, ' [', trim(units), ']'
  print '(A,1PE13.5,2A)', ' time_step  = ', dt, ' [', trim(units), ']'
  print '(A,1PE13.5,2A)', ' current    = ', t, ' [', trim(units), ']'

  !--- Integrate to end
  print *, 'Running...'
  do while (t < t_end - 0.5d0*dt)
    s = M%update();                   call chk('update', s)
    s = M%get_current_time(t);        call chk('get_current_time', s)
  end do

  ! one last step if needed (guard against floating rounding)
  if (t < t_end) then
    s = M%update_until(t_end);        call chk('update_until', s)
    s = M%get_current_time(t);        call chk('get_current_time', s)
  end if

  print '(A,1PE13.5,2A)', ' done. final time = ', t, ' [', trim(units), ']'

  !--- Finalize
  print *, 'Finalizing...'
  s = M%finalize();                   call chk('finalize', s)
  print *, 'DONE.'

contains

  subroutine chk(where, status)
    character(len=*), intent(in) :: where
    integer(c_int),   intent(in) :: status
    if (status /= BMI_SUCCESS) then
      print *, 'FATAL: ', trim(where), ' returned status=', status
      stop 2
    end if
  end subroutine chk

end program sfincs_driver

