   module sfincs_bmi
   use iso_c_binding
   use sfincs_lib
   use sfincs_ncoutput
   use sfincs_data
   
   implicit none

   public :: initialize
   public :: finalize
   public :: update
   
   public :: get_start_time
   public :: get_end_time
   public :: get_current_time
   public :: get_time_step
   public :: get_time_units
   
   public :: get_var
   public :: get_var_shape
   public :: get_var_type
   public :: get_var_rank
   
   public :: set_var
   
   public :: get_grid_type
   public :: get_grid_rank
   !public :: get_grid_shape
   public :: get_grid_size
   public :: get_grid_x
   public :: get_grid_y
   
   private

   integer(c_int), bind(C, name="maxstrlen") :: maxstrlen = 1024
   
   contains
   
!-----------------------------------------------------------------------------------------------------!
   
   function initialize(c_config_file) result(ierr) bind(C, name="initialize")
   !DEC$ ATTRIBUTES DLLEXPORT :: initialize
   
   character(kind=c_char), intent(in)   :: c_config_file(*)
   character(len=strlen(c_config_file)) :: config_file
   integer(kind=c_int) :: ierr
   
   write(*,*) "BMI init start"

   config_file = char_array_to_string(c_config_file, strlen(c_config_file))
   write(*,*) "config file: ", config_file
   ierr = sfincs_initialize(config_file)

   write(*,*) "BMI init end"
   
   end function initialize
   
!-----------------------------------------------------------------------------------------------------!
   
   function finalize() result(ierr) bind(C, name="finalize")
   !DEC$ ATTRIBUTES DLLEXPORT :: finalize
   
   integer(kind=c_int) :: ierr
   
   ierr = sfincs_finalize()
   
   end function finalize
   
!-----------------------------------------------------------------------------------------------------!
   
   function update(dt) result(ierr) bind(C, name="update")
   !DEC$ ATTRIBUTES DLLEXPORT :: update
   
   real(kind=c_double), value, intent(in)  :: dt
   integer(kind=c_int)                     :: ierr
   
   ierr = sfincs_update(dt)
   
   end function update
   
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
      
   subroutine get_current_time(tcurrent) bind(C, name="get_current_time")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_current_time
   
   real(c_double), intent(out) :: tcurrent

   tcurrent = t
   
   end subroutine get_current_time
   
!-----------------------------------------------------------------------------------------------------!  
      
   subroutine get_time_step(deltat) bind(C, name="get_time_step")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_time_step

   real(c_double), intent(out) :: deltat

   deltat = dt
   
   end subroutine get_time_step
   
!-----------------------------------------------------------------------------------------------------!  
      
   subroutine get_time_units(c_time_units) bind(C, name="get_time_units")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_time_units

   character(kind=c_char), intent(out) :: c_time_units(maxstrlen)
   
   c_time_units = string_to_char_array("s", 1)
   
   end subroutine get_time_units
   
!-----------------------------------------------------------------------------------------------------!
   
   subroutine get_var(c_var_name, x) bind(C, name="get_var")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var
   
   character(kind=c_char), intent(in)       :: c_var_name(*)
   type(c_ptr), intent(inout)               :: x

   ! The fortran name of the attribute name
   character(len=strlen(c_var_name)) :: var_name
   
   ! Store the name
   var_name = char_array_to_string(c_var_name,strlen(c_var_name))
   
   select case(var_name)
   case("z_xz") ! x grid cell centre z_xz
       x = c_loc(z_xz)
	! x grid cell centre z_yz
   case("z_yz")
       x = c_loc(z_yz)
   ! water level zs
   case("zs")
      x = c_loc(zs)
   case("zb") ! bed level
      if(subgrid) then
        x = c_loc(subgrid_z_zmin)
      else
        x = c_loc(zb)
      end if
   case("qtsrc")
     x = c_loc(qtsrc)
   case("zst_bnd")
     x = c_loc(zst_bnd)
   case default
	 write(*,*) 'get_var error'
     ! nullptr
   end select

   end subroutine get_var
   
   !-----------------------------------------------------------------------------------------------------!
   
   subroutine get_var_shape(c_var_name, var_shape) bind(C, name="get_var_shape")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var_shape
   
   character(kind=c_char), intent(in) :: c_var_name(*)
   integer(c_int), intent(inout)      :: var_shape(1)
   character(len=strlen(c_var_name))  :: var_name

   var_name = char_array_to_string(c_var_name, strlen(c_var_name))
   var_shape = (/0/)
   
   select case(var_name)
   case("z_xz", "z_yz","zs","zb","qtsrc","zst_bnd")     
      ! inverted shapes (fortran to c)
      var_shape(1) = np
   case default
     write(*,*) 'get_var_shape error'
     ! nullptr
   end select
   
   end subroutine get_var_shape
   
!-----------------------------------------------------------------------------------------------------!
   
   subroutine get_var_type(c_var_name, c_var_type) bind(C, name="get_var_type")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var_type
   
   character(kind=c_char), intent(in)  :: c_var_name(*)
   character(kind=c_char), intent(out) :: c_var_type(maxstrlen)
   character(len=maxstrlen)            :: var_name, var_type

   var_name = char_array_to_string(c_var_name, strlen(c_var_name))

   select case(var_name)
   case("z_xz", "z_yz","zs","zb","qtsrc","zst_bnd")      
      var_type = "float"
   case default
     write(*,*) 'get_var_type error'
     ! nullptr
   end select
   
   c_var_type = string_to_char_array(trim(var_type), len(trim(var_type)))
   
   end subroutine get_var_type
   
!-----------------------------------------------------------------------------------------------------!
   
   subroutine get_var_rank(c_var_name, rank) bind(C, name="get_var_rank")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var_rank
   
   character(kind=c_char), intent(in) :: c_var_name(*)
   integer(c_int), intent(out) :: rank

   ! The fortran name of the attribute name
   character(len=strlen(c_var_name)) :: var_name
   

   var_name = char_array_to_string(c_var_name, strlen(c_var_name))
   
   select case(var_name)
   case("z_xz", "z_yz","zs","zb","qtsrc","zst_bnd") 
      rank = 1
   case default
     write(*,*) 'get_var_rank error'
     !rank  = some_invalid_val
   end select

   end subroutine get_var_rank
   
!-----------------------------------------------------------------------------------------------------!   

   subroutine set_var(c_var_name, c_var_ptr) bind(C, name="set_var")
   !DEC$ ATTRIBUTES DLLEXPORT :: set_var

   character(kind=c_char), intent(in) :: c_var_name(*)
   type(c_ptr), value, intent(in) :: c_var_ptr

   real(c_float), pointer  :: f_var_ptr(:)

   ! The fortran name of the attribute name
   character(len=strlen(c_var_name)) :: var_name
   integer :: i

   var_name = char_array_to_string(c_var_name, strlen(c_var_name))
   
   call c_f_pointer(c_var_ptr, f_var_ptr, [np])
   
   select case(var_name)
   case("zs")
     do i = 1, np
       f_var_ptr(i) = zs(i)
     end do
   case("zb")
     do i = 1, np
       f_var_ptr(i) = zb(i)
     end do
   case("qtsrc")
     do i = 1, np
       f_var_ptr(i) = qtsrc(i)
     end do
   case("zst_bnd")
     do i = 1, np
       f_var_ptr(i) = zst_bnd(i)
     end do
   case default
     write(*,*) 'set_var error'
     !nullptr
   end select
         
   end subroutine set_var
   
!-----------------------------------------------------------------------------------------------------!
   
   !subroutine get_grid_x(grid_x_c_ptr) bind(C, name="get_grid_x")
   !!DEC$ ATTRIBUTES DLLEXPORT :: get_grid_x
   !
   !type(c_ptr), intent(inout) :: grid_x_c_ptr
   !real(c_float), pointer:: grid_x_f_ptr(:)
   !integer i, z_xz_size
   !   
   !if(allocated(z_xz)) then
   !  z_xz_size = size(z_xz)
   !  call c_f_pointer(grid_x_c_ptr, grid_x_f_ptr, [z_xz_size])
   !  do i = 1, z_xz_size
   !    grid_x_f_ptr(i) = z_xz(i)
   !  end do
   !end if
   !
   !end subroutine get_grid_x

   subroutine get_grid_type(grid_type) bind(C, name="get_grid_type")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_grid_rtype
   character(kind=c_char), intent(out) :: grid_type(maxstrlen)
   character(len=maxstrlen) :: string
   if(use_quadtree) then
     string = "unstructered"
   else
     string = "rectilinear"
   end if
   grid_type = string_to_char_array(trim(string), len(trim(string)))
   end subroutine get_grid_type
   
   subroutine get_grid_rank(grid_rank) bind(C, name="get_grid_rank")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_grid_rank
   integer(c_int), intent(out) :: grid_rank
   grid_rank = 2
   end subroutine get_grid_rank
   
   subroutine get_grid_size(grid_size) bind(C, name="get_grid_size")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_grid_size
   integer(c_int), intent(out) :: grid_size
   grid_size = np
   end subroutine get_grid_size
   
   subroutine get_grid_x(grid_x_c_ptr) bind(C, name="get_grid_x")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_grid_x
   type(c_ptr), intent(inout) :: grid_x_c_ptr
   if(allocated(z_xz)) then
     grid_x_c_ptr = c_loc(z_xz(1))
   end if
   end subroutine get_grid_x
   
   
   subroutine get_grid_y(grid_y_c_ptr) bind(C, name="get_grid_y")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_grid_y
   type(c_ptr), intent(inout) :: grid_y_c_ptr
   if(allocated(z_yz)) then
     grid_y_c_ptr = c_loc(z_yz(1))
   end if
   end subroutine get_grid_y
   
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