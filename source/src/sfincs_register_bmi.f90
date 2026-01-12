module sfincs_register_bmi_mod
  use, intrinsic :: iso_c_binding, only: c_ptr, c_null_ptr, c_f_pointer, c_loc
  use iso_c_bmif_2_0,              only: box
  use sfincs_bmi2,                 only: sfincs_bmi

  implicit none
  private
  public :: register_bmi

  ! The actual SFINCS BMI implementation object:
  type(sfincs_bmi), target, save :: g_bmi

  ! The ABI "box" object NGen expects the handle to point at:
  type(box), target, save :: g_box

contains

  function register_bmi(model) bind(C, name="register_bmi") result(res)
    !! NGen expects: void *register_bmi(void *model)
    !!
    !! 'model' points to NGen's "handle" storage (i.e., address of a C pointer).
    !! We must store into that handle a pointer to g_box, and make g_box%ptr => g_bmi.
    type(c_ptr), value :: model
    type(c_ptr)        :: res

    type(c_ptr), pointer :: handle

    ! Treat incoming void* as a pointer to a C pointer (the handle)
    call c_f_pointer(model, handle)

    ! Link box -> actual BMI object
    g_box%ptr => g_bmi

    ! Write the handle = &g_box
    handle = c_loc(g_box)

    res = c_null_ptr
  end function register_bmi

end module sfincs_register_bmi_mod

