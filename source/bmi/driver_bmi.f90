! ================================================================================
! Driver for the SFINCS model via BMI 2.0
!
! Modeled after the Snow17 driver you shared: initialize -> run updates -> finalize
! All BMI calls here are functions returning integer status codes (BMI_SUCCESS/BMI_FAILURE).
! ================================================================================

program multi_driver
  use bmif_2_0
  use sfincs_bmi2
  implicit none

  ! ---- model ----
  type(sfincs_bmi) :: m

  ! ---- run control ----
  character(len=400) :: config_file
  integer            :: arg_status
  integer            :: status

  ! ---- time ----
  double precision :: current_time
  double precision :: end_time
  double precision :: dt
  integer          :: nsteps

  ! -------------------- initialize --------------------
  print *, 'Initializing SFINCS via BMI...'
  call get_command_argument(1, config_file, status=arg_status)
  if (arg_status /= 0 .or. len_trim(config_file) == 0) then
    config_file = 'sfincs_config.txt'
    print *, 'No config filename supplied; using default: ', trim(config_file)
  end if

  status = m%initialize(trim(config_file))
  if (status /= BMI_SUCCESS) then
    print *, 'ERROR: initialize() failed with status = ', status
    stop 1
  end if

  ! -------------------- time info ---------------------
  status = m%get_current_time(current_time)
  if (status /= BMI_SUCCESS) then
    print *, 'ERROR: get_current_time() failed: ', status; stop 1
  end if

  status = m%get_end_time(end_time)
  if (status /= BMI_SUCCESS) then
    print *, 'ERROR: get_end_time() failed: ', status; stop 1
  end if

  status = m%get_time_step(dt)
  if (status /= BMI_SUCCESS) then
    print *, 'ERROR: get_time_step() failed: ', status; stop 1
  end if

  if (dt > 0.d0) then
    nsteps = max(0, int( (end_time - current_time)/dt + 0.5d0 ))
  else
    nsteps = 0
  end if

  print *, 'Run window (s): start =', current_time, ' end =', end_time, ' dt =', dt
  print *, 'Estimated timesteps:', nsteps

  ! -------------------- run loop ----------------------
  print *, 'Running...'
  do while (current_time .le. end_time)
    status = m%update()
    if (status /= BMI_SUCCESS) then
      print *, 'ERROR: update() failed with status = ', status
      exit
    end if
    status = m%get_current_time(current_time)
    if (status /= BMI_SUCCESS) then
      print *, 'ERROR: get_current_time() failed with status = ', status
      exit
    end if
  end do

  ! -------------------- finalize ----------------------
  print *, 'Finalizing...'
  status = m%finalize()
  if (status /= BMI_SUCCESS) then
    print *, 'WARNING: finalize() returned status = ', status
  end if
  print *, 'DONE.'

end program multi_driver

