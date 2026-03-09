module sfincs_dll
  !!
  !! XMI-style shim around the sfincs_bmi type, exposing a simple C API:
  !!   int initialize(const char* cfg);
  !!   int update();
  !!   int finalize();
  !!   int get_value_ptr(const char* name, void** ptr);
  !!   int get_var_size(const char* name);
  !!
  !! This links to the existing sfincs_bmi2.f90 BMI wrapper.
  !!
  use, intrinsic :: iso_c_binding, only: &
       c_int, c_ptr, c_char, c_null_ptr, c_null_char, &
       c_f_pointer, c_loc, c_associated
  use, intrinsic :: iso_fortran_env, only: real32

  use bmif_2_0,    only: BMI_SUCCESS, BMI_FAILURE
  use sfincs_bmi2, only: sfincs_bmi

  implicit none

  ! Single global BMI instance
  type(sfincs_bmi), save :: M

contains

  !===================================================================
  !  C-callable initialize(config_file)
  !===================================================================
  function initialize(cfg_c) bind(C, name="initialize") result(ierr)
    type(c_ptr), value :: cfg_c        !! C const char* (NUL-terminated)
    integer(c_int)     :: ierr

    character(len=:), allocatable :: cfg
    integer :: status

    call cstring_to_fortran(cfg_c, cfg)

    if (.not. allocated(cfg)) then
      ! NULL -> let BMI handle default config
      status = M%initialize('')
    else
      status = M%initialize(trim(cfg))
    end if

    ierr = int(status, c_int)
  end function initialize

  !===================================================================
  !  C-callable update()
  !===================================================================
  function update() bind(C, name="update") result(ierr)
    integer(c_int) :: ierr
    integer :: status

    status = M%update()
    ierr   = int(status, c_int)
  end function update

  !===================================================================
  !  C-callable finalize()
  !===================================================================
  function finalize() bind(C, name="finalize") result(ierr)
    integer(c_int) :: ierr
    integer :: status

    status = M%finalize()
    ierr   = int(status, c_int)
  end function finalize

  !===================================================================
  !  C-callable get_value_ptr(name, ptr)
  !
  !  NOTE:
  !    - Supports variables "zs", "zb", "depth"
  !    - Returns pointer to REAL*4 (float) array
  !===================================================================
  function get_value_ptr(name_c, ptr_c) bind(C, name="get_value_ptr") result(ierr)
    type(c_ptr), value :: name_c   !! C const char*
    type(c_ptr)        :: ptr_c    !! C void* (Fortran: C_PTR passed by reference)
    integer(c_int)     :: ierr

    character(len=:), allocatable :: name
    real(real32), pointer :: fptr(:) => null()
    integer :: stat

    call cstring_to_fortran(name_c, name)

    if (.not. allocated(name)) then
      ptr_c = c_null_ptr
      ierr  = int(BMI_FAILURE, c_int)
      return
    end if

    select case (trim(name))
    case ('zs')
       stat = M%get_value_ptr_float('zs', fptr)
       if (stat == BMI_SUCCESS .and. associated(fptr)) then
         ptr_c = c_loc(fptr(1))
         ierr  = int(BMI_SUCCESS, c_int)
       else
         ptr_c = c_null_ptr
         ierr  = int(BMI_FAILURE, c_int)
       end if

    case ('zb')
       stat = M%get_value_ptr_float('zb', fptr)
       if (stat == BMI_SUCCESS .and. associated(fptr)) then
         ptr_c = c_loc(fptr(1))
         ierr  = int(BMI_SUCCESS, c_int)
       else
         ptr_c = c_null_ptr
         ierr  = int(BMI_FAILURE, c_int)
       end if

    case ('depth')
       stat = M%get_value_ptr_float('depth', fptr)
       if (stat == BMI_SUCCESS .and. associated(fptr)) then
         ptr_c = c_loc(fptr(1))
         ierr  = int(BMI_SUCCESS, c_int)
       else
         ptr_c = c_null_ptr
         ierr  = int(BMI_FAILURE, c_int)
       end if

    case default
       ptr_c = c_null_ptr
       ierr  = int(BMI_FAILURE, c_int)
    end select

  end function get_value_ptr

  !===================================================================
  !  C-callable get_var_size(name)
  !===================================================================
  function get_var_size(name_c) bind(C, name="get_var_size") result(out)
    type(c_ptr), value :: name_c   !! C const char*
    integer(c_int)     :: out

    character(len=:), allocatable :: name

    call cstring_to_fortran(name_c, name)

    if (.not. allocated(name)) then
      out = 0_c_int
      return
    end if

    if (trim(name) == 'zs' .or. trim(name) == 'zb' .or. trim(name) == 'depth') then
      out = int(M%nx * M%ny, c_int)
    else
      out = 0_c_int
    end if
  end function get_var_size

  !===================================================================
  !  Helper: Convert C string (char*) -> allocatable Fortran string
  !===================================================================
  subroutine cstring_to_fortran(cstr, fstr)
    type(c_ptr), value :: cstr
    character(len=:), allocatable :: fstr

    character(kind=c_char), pointer :: p_chars(:)
    integer :: maxlen, n, i
    character(len=1, kind=c_char) :: ch_c
    character(len=1)              :: ch_f

    ! If NULL, leave fstr unallocated so callers can treat as "no string"
    if (.not. c_associated(cstr)) then
      return
    end if

    ! Arbitrary upper bound; just to scan to first NUL
    maxlen = 10000
    call c_f_pointer(cstr, p_chars, [maxlen])

    n = 0
    do i = 1, maxlen
       if (p_chars(i) == c_null_char) exit
       n = n + 1
    end do

    if (n <= 0) then
      allocate(character(len=0) :: fstr)
      return
    end if

    allocate(character(len=n) :: fstr)
    do i = 1, n
      ch_c = p_chars(i)
      ! transfer from c_char (1 byte) to default character(1)
      ch_f = transfer(ch_c, ch_f)
      fstr(i:i) = ch_f
    end do
  end subroutine cstring_to_fortran

end module sfincs_dll

