use bmif_2_0_mod, only: Bmi_F90_Interface

type, extends(Bmi_F90_Interface) :: Sfincs_BMI
  ! model-specific state variables go here, e.g.:
  real(real64), allocatable :: some_state(:)
  logical :: is_initialized = .false.
end type Sfincs_BMI

contains

subroutine initialize(self, config_file) bind(c)
  class(Sfincs_BMI), intent(inout) :: self
  character(len=*), intent(in) :: config_file
  ! parse config_file, allocate state, set up grid and initial time
  self%is_initialized = .true.
end subroutine initialize

subroutine update(self) bind(c)
  class(Sfincs_BMI), intent(inout) :: self
  call self%advance_time(self%time_step)
end subroutine update

subroutine update_until(self, target_time) bind(c)
  class(Sfincs_BMI), intent(inout) :: self
  real(real64), intent(in) :: target_time
  do while (self%current_time < target_time)
    call self%update()
  end do
end subroutine update_until

subroutine finalize(self) bind(c)
  class(Sfincs_BMI), intent(inout) :: self
  ! deallocate state and clean up
  self%is_initialized = .false.
end subroutine finalize


function get_current_time(self) result(time) bind(c)
  class(Sfincs_BMI), intent(in) :: self
  real(real64) :: time
  time = self%current_time
end function

function get_start_time(self) result(time) bind(c)
  class(Sfincs_BMI), intent(in) :: self
  real(real64) :: time
  time = self%start_time
end function

function get_end_time(self) result(time) bind(c)
  class(Sfincs_BMI), intent(in) :: self
  real(real64) :: time
  time = self%end_time
end function

function get_time_step(self) result(dt) bind(c)
  class(Sfincs_BMI), intent(in) :: self
  real(real64) :: dt
  dt = self%time_step
end function

function get_time_units(self) result(units) bind(c)
  class(Sfincs_BMI), intent(in) :: self
  character(len=64) :: units
  units = 'seconds' ! or as appropriate
end function


function get_var_nbytes(self, name) result(nbytes) bind(c)
  class(Sfincs_BMI), intent(in) :: self
  character(len=*), intent(in) :: name
  integer(c_size_t) :: nbytes
  if (trim(name) == 'temperature') then
    nbytes = size(self%temperature) * storage_size(self%temperature)/8
  else
    nbytes = 0
  end if
end function get_var_nbytes

subroutine get_value(self, name, dst) bind(c)
  class(Sfincs_BMI), intent(inout) :: self
  character(len=*), intent(in) :: name
  type(c_ptr), intent(out) :: dst
  if (trim(name) == 'temperature') then
    call c_f_pointer(dst, self%temperature)
  else
    dst = c_null_ptr
  end if
end subroutine get_value

