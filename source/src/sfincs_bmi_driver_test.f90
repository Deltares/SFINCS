program sfincs_bmi_driver_test

  use sfincs_bmi_module
  use bmif_2_0
  implicit none

  type(Sfincs_BMI) :: model
  real(real64) :: time, dt, endtime
  character(len=128) :: config_file
  integer :: i

  config_file = "sfincs_config.txt"  ! Replace with your real input

  call model%initialize(trim(config_file))
  print *, "Initialized model"

  time = model%get_current_time()
  dt = model%get_time_step()
  endtime = model%get_end_time()

  print *, "Start time: ", time
  print *, "Time step: ", dt
  print *, "End time: ", endtime

  do i = 1, 10
    call model%update()
    time = model%get_current_time()
    print *, "Step ", i, " Time: ", time
  end do

  call model%finalize()
  print *, "Finalized model"

end program sfincs_bmi_driver_test

