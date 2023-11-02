   module sfincs_bmi
   use iso_c_binding
   use sfincs_lib
   use sfincs_ncoutput
   use sfincs_data
   
   implicit none

   public :: initialize
   public :: finalize

   public :: update_until
   public :: update
   
   public :: get_start_time
   public :: get_end_time
   public :: get_current_time
   public :: get_time_step
   public :: get_time_units
   
   public :: get_value
   public :: get_value_ptr
   public :: get_var_shape
   public :: get_var_type
   public :: get_var_rank
   
   public :: set_value
   
   public :: get_grid_type
   public :: get_grid_rank
   !public :: get_grid_shape
   public :: get_grid_size
   public :: get_grid_x
   public :: get_grid_y
   
   private

   integer(c_int), bind(C, name="maxstrlen") :: maxstrlen = 1024
   
   type return_code_type
     integer(c_int) :: success = 0
     integer(c_int) :: failure = 1
   end type return_code_type
   
   type (return_code_type) :: ret_code
   
   contains
   
!-----------------------------------------------------------------------------------------------------!
   
   function initialize(c_config_file) result(ierr) bind(C, name="initialize")
   !DEC$ ATTRIBUTES DLLEXPORT :: initialize
   
   character(kind=c_char), intent(in)   :: c_config_file(*)
   character(len=strlen(c_config_file)) :: config_file
   integer(kind=c_int) :: ierr

   config_file = char_array_to_string(c_config_file, strlen(c_config_file))
   ierr = sfincs_initialize(config_file)

   end function initialize
   
!-----------------------------------------------------------------------------------------------------!
   
   function finalize() result(ierr) bind(C, name="finalize")
   !DEC$ ATTRIBUTES DLLEXPORT :: finalize
   
   integer(kind=c_int) :: ierr
   
   ierr = sfincs_finalize()
   
   end function finalize

!-----------------------------------------------------------------------------------------------------!
   
   function update_until(t) result(ierr) bind(C, name="update_until")
   !DEC$ ATTRIBUTES DLLEXPORT :: update_until
   
   real(kind=c_double), value, intent(in)  :: t
   integer(kind=c_int)                     :: ierr
   
   ierr = sfincs_update(t)
   
   end function update_until

!-----------------------------------------------------------------------------------------------------!
   
   function update() result(ierr) bind(C, name="update")
   !DEC$ ATTRIBUTES DLLEXPORT :: update
   
   integer(kind=c_int) :: ierr
   
   ierr = sfincs_update(real(dt, 8)) ! this is bad
   
   end function update
   
!-----------------------------------------------------------------------------------------------------!  
      
   function get_start_time(start_time) result(ierr) bind(C, name="get_start_time")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_start_time
   
   real(c_double), intent(out) :: start_time
   integer(kind=c_int) :: ierr

   start_time = t0
   ierr = ret_code%success
   
   end function get_start_time

!-----------------------------------------------------------------------------------------------------!  
      
   function get_end_time(end_time) result(ierr) bind(C, name="get_end_time")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_end_time
   
   real(c_double), intent(out) :: end_time
   integer(kind=c_int) :: ierr

   end_time = t1
   ierr = ret_code%success
   
   end function get_end_time
   
!-----------------------------------------------------------------------------------------------------!  
      
   function get_current_time(current_time) result(ierr) bind(C, name="get_current_time")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_current_time
   
   real(c_double), intent(out) :: current_time
   integer(kind=c_int) :: ierr

   current_time = t
   ierr = ret_code%success
   
   end function get_current_time
   
!-----------------------------------------------------------------------------------------------------!  
      
   function get_time_step(time_step) result(ierr) bind(C, name="get_time_step")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_time_step

   real(c_double), intent(out) :: time_step
   integer(kind=c_int) :: ierr

   time_step = dt
   ierr = ret_code%success
   
   end function get_time_step
   
!-----------------------------------------------------------------------------------------------------!  
      
   function get_time_units(time_units) result(ierr) bind(C, name="get_time_units")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_time_units

   character(kind=c_char), intent(out) :: time_units(maxstrlen)
   integer(kind=c_int) :: ierr
   
   time_units = string_to_char_array("s", 1)
   ierr = ret_code%success
   
   end function get_time_units
   
!-----------------------------------------------------------------------------------------------------!
   
   function get_value(var_name, ptr) result(ierr) bind(C, name="get_value")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_value
   
   character(kind=c_char), intent(inout) :: var_name(*)
   type(c_ptr), intent(inout) :: ptr
   integer(kind=c_int) :: ierr
 
   ! The fortran name of the attribute name
   character(len=strlen(var_name)) :: f_var_name
   real(c_float), pointer:: f_ptr(:)
   integer :: i
   
   ierr = ret_code%success
   
   ! get the name
   f_var_name = char_array_to_string(var_name, strlen(var_name))
      
   select case(f_var_name)
   case("z_xz") ! x grid cell centre
	 call c_f_pointer(ptr, f_ptr, [np])
     do i = 1, np
       f_ptr(i) = z_xz(i)
	 end do
   case("z_yz") ! y grid cell centre
	 call c_f_pointer(ptr, f_ptr, [np])
     do i = 1, np
       f_ptr(i) = z_yz(i)
	 end do
   case("zs") ! water level
	 call c_f_pointer(ptr, f_ptr, [np])
     do i = 1, np
       f_ptr(i) = zs(i)
	 end do
   case("zb") ! bed level
   	 call c_f_pointer(ptr, f_ptr, [np]) 
	 if(subgrid) then
	   do i = 1, np
         f_ptr(i) = subgrid_z_zmin(i)
	   end do
	 else
	   do i = 1, np
         f_ptr(i) = subgrid_z_zmin(i)
	   end do
     end if
   case("qtsrc")
	 call c_f_pointer(ptr, f_ptr, [np])
     do i = 1, np
       f_ptr(i) = qtsrc(i)
	 end do
   case("zst_bnd")
	 call c_f_pointer(ptr, f_ptr, [np])
     do i = 1, np
       f_ptr(i) = zst_bnd(i)
	 end do
   case default
     write(*,*) 'get_value error'
     ierr = ret_code%failure
     ptr = c_null_ptr
   end select
 
   end function get_value
   
!-----------------------------------------------------------------------------------------------------!
   
   function get_value_ptr(var_name, ptr) result(ierr) bind(C, name="get_value_ptr")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_value_ptr
   
   character(kind=c_char), intent(inout) :: var_name(*)
   type(c_ptr), intent(inout) :: ptr
   integer(kind=c_int) :: ierr

   ! The fortran name of the attribute name
   character(len=strlen(var_name)) :: f_var_name
   
   ierr = ret_code%success
   
   ! Store the name
   f_var_name = char_array_to_string(var_name, strlen(var_name))
   
   select case(f_var_name)
   case("z_xz") ! x grid cell centre
       ptr = c_loc(z_xz)
   case("z_yz") ! y grid cell centre
       ptr = c_loc(z_yz)
   case("zs") ! water level
      ptr = c_loc(zs)
   case("zb") ! bed level
      if(subgrid) then
        ptr = c_loc(subgrid_z_zmin)
      else
        ptr = c_loc(zb)
      end if
   case("qtsrc")
     ptr = c_loc(qtsrc)
   case("zst_bnd")
     ptr = c_loc(zst_bnd)
   case default
	 write(*,*) 'get_value_ptr error'
     ierr = ret_code%failure
     ptr = c_null_ptr
   end select

   end function get_value_ptr
   
   !-----------------------------------------------------------------------------------------------------!
   
   function get_var_shape(var_name, var_shape) result(ierr) bind(C, name="get_var_shape")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var_shape
   
   character(kind=c_char), intent(in) :: var_name(*)
   integer(c_int), intent(inout) :: var_shape(1)
   integer(kind=c_int) :: ierr
   
   character(len=strlen(var_name))  :: f_var_name

   f_var_name = char_array_to_string(var_name, strlen(var_name))
   var_shape = (/0/)
   
   select case(f_var_name)
   case("z_xz", "z_yz","zs","zb","qtsrc","zst_bnd")     
      ! inverted shapes (fortran to c)
      var_shape(1) = np
      ierr = ret_code%success
   case default
     write(*,*) 'get_var_shape error'
     ierr = ret_code%failure
   end select
   
   end function get_var_shape
   
!-----------------------------------------------------------------------------------------------------!
   
   function get_var_type(var_name, var_type) result(ierr) bind(C, name="get_var_type")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var_type
   
   character(kind=c_char), intent(in)  :: var_name(*)
   character(kind=c_char), intent(out) :: var_type(maxstrlen)
   integer(kind=c_int) :: ierr
   
   character(len=maxstrlen) :: f_var_name, f_var_type

   f_var_name = char_array_to_string(var_name, strlen(var_name))

   select case(f_var_name)
   case("z_xz", "z_yz", "zs", "zb", "qtsrc", "zst_bnd")      
      f_var_type = "float"
      var_type = string_to_char_array(trim(f_var_type), len(trim(f_var_type)))
      ierr = ret_code%success
   case default
     write(*,*) 'get_var_type error'
     ierr = ret_code%failure
   end select
   
   end function get_var_type
   
!-----------------------------------------------------------------------------------------------------!
   
   function get_var_rank(var_name, rank) result(ierr) bind(C, name="get_var_rank")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var_rank
   
   character(kind=c_char), intent(in) :: var_name(*)
   integer(c_int), intent(out) :: rank
   integer(kind=c_int) :: ierr

   character(len=strlen(var_name)) :: f_var_name

   f_var_name = char_array_to_string(var_name, strlen(var_name))
   
   select case(f_var_name)
   case("z_xz", "z_yz", "zs", "zb" ,"qtsrc", "zst_bnd") 
      rank = 1
      ierr = ret_code%success
   case default
     write(*,*) 'get_var_rank error'
     ierr = ret_code%failure
   end select

   end function get_var_rank
   
!-----------------------------------------------------------------------------------------------------!   

   function set_value(var_name, ptr) result(ierr) bind(C, name="set_value")
   !DEC$ ATTRIBUTES DLLEXPORT :: set_value

   character(kind=c_char), intent(in) :: var_name(*)
   type(c_ptr), value, intent(in) :: ptr
   integer(kind=c_int) :: ierr

   real(c_float), pointer  :: f_ptr(:)

   ! The fortran name of the attribute name
   character(len=strlen(var_name)) :: f_var_name
   integer :: i

   ierr = ret_code%success
   
   f_var_name = char_array_to_string(var_name, strlen(var_name))
   
   call c_f_pointer(ptr, f_ptr, [np])
   
   select case(f_var_name)
   case("zs")
     do i = 1, np
       zs(i) = f_ptr(i)
     end do
   case("zb")
     do i = 1, np
       zb(i) = f_ptr(i)
     end do
   case("qtsrc")
     do i = 1, np
       qtsrc(i) = f_ptr(i)
     end do
   case("zst_bnd")
     do i = 1, np
       zst_bnd(i) = f_ptr(i)
     end do
   case default
     write(*,*) 'set_value error'
     ierr = ret_code%failure
   end select
         
   end function set_value
   
!-----------------------------------------------------------------------------------------------------!

   function get_grid_type(grid_type) result(ierr) bind(C, name="get_grid_type")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_grid_rtype
   
   character(kind=c_char), intent(out) :: grid_type(maxstrlen)
   integer(kind=c_int) :: ierr
   
   character(len=maxstrlen) :: string
   
   if(use_quadtree) then
     string = "unstructered"
   else
     string = "rectilinear"
   end if
   grid_type = string_to_char_array(trim(string), len(trim(string)))
   ierr = ret_code%success
   
   end function get_grid_type
   
   !-----------------------------------------------------------------------------------------------------!
   
   function get_grid_rank(grid_rank) result(ierr) bind(C, name="get_grid_rank")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_grid_rank
   
   integer(c_int), intent(out) :: grid_rank
   integer(kind=c_int) :: ierr
   
   grid_rank = 2
   ierr = ret_code%success

   end function get_grid_rank
   
   !-----------------------------------------------------------------------------------------------------!
   
   function get_grid_size(grid_size) result(ierr) bind(C, name="get_grid_size")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_grid_size

   integer(c_int), intent(out) :: grid_size
   integer(kind=c_int) :: ierr
   
   grid_size = np
   ierr = ret_code%success
    
   end function get_grid_size
   
   !-----------------------------------------------------------------------------------------------------!
   
   function get_grid_x(ptr) result(ierr) bind(C, name="get_grid_x")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_grid_x
   
   type(c_ptr), intent(inout) :: ptr
   integer(kind=c_int) :: ierr
   
   real(c_float), pointer:: f_ptr(:)

   integer i, sz
      
   if(allocated(z_xz)) then
     sz = size(z_xz)
     call c_f_pointer(ptr, f_ptr, [sz])
     do i = 1, sz
       f_ptr(i) = z_xz(i)
     end do
     ierr = ret_code%success
   else
     write(*,*) 'get_grid_x error'
     ptr = c_null_ptr
     ierr = ret_code%failure
   end if
   
   end function get_grid_x
   
   !-----------------------------------------------------------------------------------------------------!
   
   function get_grid_y(ptr) result(ierr) bind(C, name="get_grid_y")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_grid_y
   
   type(c_ptr), intent(inout) :: ptr
   integer(kind=c_int) :: ierr
   
   real(c_float), pointer:: f_ptr(:)

   integer i, sz
      
   if(allocated(z_yz)) then
     sz = size(z_yz)
     call c_f_pointer(ptr, f_ptr, [sz])
     do i = 1, sz
       f_ptr(i) = z_yz(i)
     end do
     ierr = ret_code%success
   else
     write(*,*) 'get_grid_y error'
     ptr = c_null_ptr
     ierr = ret_code%failure
   end if
   
   end function get_grid_y
   
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