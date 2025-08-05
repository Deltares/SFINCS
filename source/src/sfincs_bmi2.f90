! BMI 2.0-compliant wrapper for SFINCS
module sfincs_bmi_module
  use iso_c_binding
  use iso_fortran_env, only: real64
  use bmif_2_0_mod, only: Bmi_F90_Interface
  use sfincs_core, only: Sfincs_State, init_sfincs, update_sfincs, finalize_sfincs
  implicit none

  type, extends(Bmi_F90_Interface) :: Sfincs_BMI
    type(Sfincs_State) :: state
  end type Sfincs_BMI

contains

  subroutine initialize(self, config_file) bind(c)
    class(Sfincs_BMI), intent(inout) :: self
    character(len=*), intent(in) :: config_file
    call init_sfincs(self%state, trim(config_file))
  end subroutine

  subroutine update(self) bind(c)
    class(Sfincs_BMI), intent(inout) :: self
    call update_sfincs(self%state)
  end subroutine

  subroutine update_until(self, time) bind(c)
    class(Sfincs_BMI), intent(inout) :: self
    real(real64), intent(in) :: time
    do while (self%state%current_time < time)
      call update_sfincs(self%state)
    end do
  end subroutine

  subroutine finalize(self) bind(c)
    class(Sfincs_BMI), intent(inout) :: self
    call finalize_sfincs(self%state)
  end subroutine

  function get_current_time(self) result(time) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    real(real64) :: time
    time = self%state%current_time
  end function

  function get_start_time(self) result(time) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    real(real64) :: time
    time = 0.0
  end function

  function get_end_time(self) result(time) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    real(real64) :: time
    time = 999999.0d0
  end function

  function get_time_step(self) result(dt) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    real(real64) :: dt
    dt = self%state%dt
  end function

  function get_time_units(self) result(units) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    character(len=64) :: units
    units = 'seconds'
  end function

  subroutine get_value(self, name, dest) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    character(len=*), intent(in) :: name
    type(c_ptr), intent(in) :: dest
    real(real64), pointer :: p(:)

    select case (trim(name))
    case ("zs")
      call c_f_pointer(dest, p, [size(self%state%zs)])
      p = reshape(self%state%zs, [size(self%state%zs)])
    case ("qx")
      call c_f_pointer(dest, p, [size(self%state%qx)])
      p = reshape(self%state%qx, [size(self%state%qx)])
    case ("qy")
      call c_f_pointer(dest, p, [size(self%state%qy)])
      p = reshape(self%state%qy, [size(self%state%qy)])
    end select
  end subroutine

  subroutine get_value_ptr(self, name, ptr) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    character(len=*), intent(in) :: name
    type(c_ptr), intent(out) :: ptr

    select case (trim(name))
    case ("zs")
      ptr = c_loc(self%state%zs)
    case ("qx")
      ptr = c_loc(self%state%qx)
    case ("qy")
      ptr = c_loc(self%state%qy)
    case default
      ptr = c_null_ptr
    end select
  end subroutine

  function get_var_units(self, name) result(units) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    character(len=*), intent(in) :: name
    character(len=64) :: units

    select case (trim(name))
    case ("zs")
      units = "m"
    case ("qx", "qy")
      units = "m2/s"
    case default
      units = ""
    end select
  end function

  function get_var_type(self, name) result(vartype) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    character(len=*), intent(in) :: name
    character(len=64) :: vartype
    vartype = "double"
  end function

  function get_input_item_count(self) result(count) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    integer :: count
    count = 0
  end function

  function get_output_item_count(self) result(count) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    integer :: count
    count = 3
  end function

  subroutine get_output_var_name(self, index, name) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    integer, intent(in) :: index
    character(len=64), intent(out) :: name

    select case(index)
    case(0)
      name = "zs"
    case(1)
      name = "qx"
    case(2)
      name = "qy"
    case default
      name = ""
    end select
  end subroutine

  function get_var_grid(self, name) result(grid_id) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    character(len=*), intent(in) :: name
    integer :: grid_id
    grid_id = 0
  end function

  function get_grid_shape(self, grid, shape) result(status) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    integer, intent(in) :: grid
    integer, dimension(2), intent(out) :: shape
    integer :: status
    shape(1) = self%state%m
    shape(2) = self%state%n
    status = 0
  end function

  function get_grid_rank(self, grid) result(rank) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    integer, intent(in) :: grid
    integer :: rank
    rank = 2
  end function

  function get_grid_type(self, grid) result(grid_type) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    integer, intent(in) :: grid
    character(len=64) :: grid_type
    grid_type = "uniform_rectilinear"
  end function

  subroutine set_value(self, name, src) bind(c)
    class(Sfincs_BMI), intent(inout) :: self
    character(len=*), intent(in) :: name
    type(c_ptr), intent(in) :: src
    real(real64), pointer :: p(:)

    select case (trim(name))
    case ("zs")
      call c_f_pointer(src, p, [size(self%state%zs)])
      self%state%zs = reshape(p, shape(self%state%zs))
    case ("qx")
      call c_f_pointer(src, p, [size(self%state%qx)])
      self%state%qx = reshape(p, shape(self%state%qx))
    case ("qy")
      call c_f_pointer(src, p, [size(self%state%qy)])
      self%state%qy = reshape(p, shape(self%state%qy))
    end select
  end subroutine

    function get_var_names(self) result(names) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    character(len=128), dimension(5) :: names

    names(1) = "zs"
    names(2) = "qx"
    names(3) = "qy"
    names(4) = "elevation"
    names(5) = "mask"
  end function

  function get_var_type(self, name) result(type_str) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    character(len=128), intent(in) :: name
    character(len=128) :: type_str

    select case (trim(name))
    case ("zs", "qx", "qy", "elevation", "mask")
      type_str = "double"
    case default
      type_str = "unknown"
    end select
  end function

  subroutine get_value(self, name, dest) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    character(len=*), intent(in) :: name
    real(real64), intent(out) :: dest(:)

    select case (trim(name))
    case ("zs")
      dest = reshape(self%state%zs, [size(dest)])
    case ("qx")
      dest = reshape(self%state%qx, [size(dest)])
    case ("qy")
      dest = reshape(self%state%qy, [size(dest)])
    case ("elevation")
      dest = reshape(self%state%elevation, [size(dest)])
    case ("mask")
      dest = reshape(self%state%mask, [size(dest)])
    end select
  end subroutine

  subroutine get_value_ptr(self, name, ptr) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    character(len=*), intent(in) :: name
    real(real64), pointer :: ptr(:,:)

    select case (trim(name))
    case ("zs")
      ptr => self%state%zs
    case ("qx")
      ptr => self%state%qx
    case ("qy")
      ptr => self%state%qy
    case ("elevation")
      ptr => self%state%elevation
    case ("mask")
      ptr => self%state%mask
    end select
  end subroutine

  function get_grid_shape(self, grid_id, shape) result(ierr) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    integer, intent(in) :: grid_id
    integer, intent(out) :: shape(2)
    integer :: ierr

    shape(1) = self%state%m
    shape(2) = self%state%n
    ierr = 0
  end function

  function get_grid_type(self, grid_id) result(type_str) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    integer, intent(in) :: grid_id
    character(len=128) :: type_str
    type_str = "uniform_rectilinear"
  end function

  function get_var_grid(self, name) result(grid_id) bind(c)
    class(Sfincs_BMI), intent(in) :: self
    character(len=*), intent(in) :: name
    integer :: grid_id
    grid_id = 0
  end function


end module sfincs_bmi_module

