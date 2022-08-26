   module sfincs_bmi
   use iso_c_binding
   use sfincs_lib
   use sfincs_ncoutput
   use sfincs_data
   
   implicit none

   public :: initialize
   public :: update
   public :: finalize
   public :: get_var_shape
   public :: get_var_type
   public :: get_var_rank
   public :: set_var
   public :: get_start_time
   public :: get_end_time
   public :: get_time_step
   public :: get_current_time

   private

   integer(c_int), bind(C, name="maxstrlen") :: maxstrlen = 1024
   integer(c_int), bind(C, name="maxdims") :: maxdims = 6
   
   contains
   
!-----------------------------------------------------------------------------------------------------!
   
   function initialize(c_config_file) result(ierr) bind(C, name="initialize")
   !DEC$ ATTRIBUTES DLLEXPORT :: initialize
   
   character(kind=c_char),intent(in)    :: c_config_file(maxstrlen)
   character(len=strlen(c_config_file)) :: config_file

   integer :: ierr

   ierr = sfincs_initialize(config_file)

   end function initialize
   
!-----------------------------------------------------------------------------------------------------!
   
   function update(dt) result(ierr) bind(C, name="update")
   !DEC$ ATTRIBUTES DLLEXPORT :: update
   
   real(kind=c_double), value, intent(in)  :: dt
   integer(kind=c_int)              :: ierr
   
   ierr = sfincs_update(dt)
   
   end function update
   
!-----------------------------------------------------------------------------------------------------!
   
   function finalize() result(ierr) bind(C, name="finalize")
   !DEC$ ATTRIBUTES DLLEXPORT :: finalize
   
   integer(kind=c_int)              :: ierr
   
   ierr = sfincs_finalize()
   
   end function finalize
   
!-----------------------------------------------------------------------------------------------------!
   
   subroutine get_var_shape(c_var_name, var_shape) bind(C, name="get_var_shape")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var_shape
   
   character(kind=c_char), intent(in) :: c_var_name(*)
   integer(c_int), intent(inout)      :: var_shape(maxdims)
   character(len=strlen(c_var_name))  :: var_name

   var_name = char_array_to_string(c_var_name, strlen(c_var_name))
   var_shape = (/0, 0, 0, 0, 0, 0/)
   
   select case(var_name)
   case("xg", "yg","zs","zb","u","v")     
      ! inverted shapes (fortran to c)
      var_shape(2) = nmax
      var_shape(1) = mmax
   end select
   
   end subroutine get_var_shape
   
!-----------------------------------------------------------------------------------------------------!
   
   subroutine get_var_type(c_var_name, c_type) bind(C, name="get_var_type")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var_type
   
   character(kind=c_char), intent(in)  :: c_var_name(*)
   character(kind=c_char), intent(out) ::  c_type(maxstrlen)
   character(len=maxstrlen)            :: type_name, var_name

   var_name = char_array_to_string(c_var_name, strlen(c_var_name))

   select case(var_name)
   case("xg", "yg", "zs", "zb", "u", "v")      
      type_name = "float"  
   end select
   
   c_type = string_to_char_array(trim(type_name), len(trim(type_name)))
   
   end subroutine get_var_type
   
!-----------------------------------------------------------------------------------------------------!
   
   subroutine get_var(c_var_name, x) bind(C, name="get_var")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var
   
   character(kind=c_char), intent(in)       :: c_var_name(*)
   type(c_ptr), intent(inout)               :: x
   real(c_float), target, allocatable, save :: xf(:,:)
   integer                                  :: nm

   ! The fortran name of the attribute name
   character(len=strlen(c_var_name)) :: var_name
   
   ! Store the name
   var_name = char_array_to_string(c_var_name,strlen(c_var_name))
   
   if(allocated(xf)) then
      deallocate(xf)   
   endif
   
   select case(var_name)
   case("xg") 
!      x = c_loc(xg)
   case("yg") 
!      x = c_loc(yg)
   case("zs")    
      allocate(xf(nmax,mmax))
      xf = FILL_VALUE       ! set to fill value
      do nm = 1, np
!         xf(index_v_n(nm),index_v_m(nm)) = zs(nm)
      enddo
      x = c_loc(xf)
   case("zb") 
      allocate(xf(nmax,mmax))
      xf = FILL_VALUE       ! set to fill value
      do nm = 1, np
!         xf(index_v_n(nm),index_v_m(nm)) = zb(nm)
      enddo
      x = c_loc(xf)
   case("u") 
      allocate(xf(nmax,mmax))
      xf = FILL_VALUE       ! set to fill value
      do nm = 1, np
!         xf(index_v_n(nm),index_v_m(nm)) = u(nm)
      enddo
      x = c_loc(xf)
   case("v") 
      allocate(xf(nmax,mmax))
      xf = FILL_VALUE       ! set to fill value
      do nm = 1, np
!         xf(index_v_n(nm),index_v_m(nm)) = v(nm)
      enddo
      x = c_loc(xf)
   end select

   end subroutine get_var
   
!-----------------------------------------------------------------------------------------------------!
   
   subroutine get_var_rank(c_var_name, rank) bind(C, name="get_var_rank")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var_rank
   
   character(kind=c_char), intent(in) :: c_var_name(*)
   integer(c_int), intent(out) :: rank

   ! The fortran name of the attribute name
   character(len=strlen(c_var_name)) :: var_name
   

   var_name = char_array_to_string(c_var_name, strlen(c_var_name))
   
   select case(var_name)
   case("xg","yg","zs","zb","u","v") 
      rank = 2
   end select

   end subroutine get_var_rank
   
!-----------------------------------------------------------------------------------------------------!   

   subroutine set_var(c_var_name, xptr) bind(C, name="set_var")
   !DEC$ ATTRIBUTES DLLEXPORT :: set_var

   character(kind=c_char), intent(in) :: c_var_name(*)
   type(c_ptr), value, intent(in) :: xptr

   real(c_float), pointer  :: x_1d_float_ptr(:)
   real(c_float), pointer  :: x_2d_float_ptr(:,:)

   ! The fortran name of the attribute name
   character(len=strlen(c_var_name)) :: var_name
   integer :: nm

   var_name = char_array_to_string(c_var_name, strlen(c_var_name))
   
   select case(var_name)
   
   case("xg")
   call c_f_pointer(xptr, x_2d_float_ptr, (/ nmax, mmax /))
!   xg = x_2d_float_ptr

   case("yg")
   call c_f_pointer(xptr, x_2d_float_ptr, (/ nmax, mmax /))
!   yg = x_2d_float_ptr
   
   case("zs")
   call c_f_pointer(xptr, x_2d_float_ptr, (/ nmax, mmax /))
   do nm = 1, np
!      zs(nm) =x_2d_float_ptr( index_v_n(nm), index_v_m(nm))
   enddo
   
   end select
         
   end subroutine set_var

!-----------------------------------------------------------------------------------------------------!  
      
   subroutine get_start_time(tstart) bind(C, name="get_start_time")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_start_time
   real(c_double), intent(out) :: tstart

   tstart = t0
   end subroutine get_start_time

!-----------------------------------------------------------------------------------------------------!  
      
   subroutine get_end_time(tend) bind(C, name="get_end_time")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_end_time
   
   real(c_double), intent(out) :: tend

   tend = t1
   end subroutine get_end_time
   
!-----------------------------------------------------------------------------------------------------!  
      
   subroutine get_time_step(deltat) bind(C, name="get_time_step")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_time_step

   real(c_double), intent(out) :: deltat

   deltat = dt
   end subroutine get_time_step
   
!-----------------------------------------------------------------------------------------------------!  
      
   subroutine get_current_time(tcurrent) bind(C, name="get_current_time")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_current_time
   
   real(c_double), intent(out) :: tcurrent

   tcurrent = t
   end subroutine get_current_time
   
!-----------------------------------------------------------------------------------------------------!     
   ! private functions
   integer(c_int) pure function strlen(char_array)
   character(c_char), intent(in) :: char_array(maxstrlen)
   integer :: inull, i
   strlen = 0
   do i = 1, size(char_array)
      if (char_array(i) .eq. C_NULL_CHAR) then
         strlen = i-1
         exit
      end if
   end do
   end function strlen
!!!!!!!!!!!!!!!!!!!!!!!   
   pure function char_array_to_string(char_array, length)
   integer(c_int), intent(in) :: length
   character(c_char),intent(in) :: char_array(length)
   character(len=length) :: char_array_to_string
   integer :: i
   do i = 1, length
      char_array_to_string(i:i) = char_array(i)
   enddo
   end function char_array_to_string
!!!!!!!!!!!!!!!!!!!!!!!      
   pure function string_to_char_array(string, length)
   integer(c_int),intent(in) :: length
   character(len=length), intent(in) :: string
   character(kind=c_char,len=1) :: string_to_char_array(length+1)
   integer :: i
   do i = 1, length
      string_to_char_array(i) = string(i:i)
   enddo
   string_to_char_array(length+1) = C_NULL_CHAR
   end function string_to_char_array

   end module sfincs_bmi