module sfincs_bmi
   use iso_c_binding
   use sfincs_lib
   use sfincs_ncoutput
   use sfincs_data

   implicit none

   public :: initialize
   public :: finalize
   public :: get_component_name

   public :: update_until
   public :: update

   public :: get_start_time
   public :: get_end_time
   public :: get_current_time
   public :: get_time_step
   public :: get_time_units

   public :: get_value
   public :: get_value_at_indices
   public :: get_value_ptr
   public :: get_var_shape
   public :: get_var_type
   public :: get_var_rank

   public :: set_value
   public :: set_value_at_indices

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

   function get_component_name(component_name) result(ierr) bind(C, name="get_component_name")
      !DEC$ ATTRIBUTES DLLEXPORT :: get_component_name

      character(kind=c_char), intent(out) :: component_name(maxstrlen)
      integer(kind=c_int) :: ierr

      character(len=maxstrlen) :: string

      string = "Sfincs hydrodynamic model (C)"
      component_name = string_to_char_array(trim(string), len(trim(string)))
      ierr = ret_code%success

   end function get_component_name

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

   function get_value(var_name, dest, n) result(ierr) bind(C, name="get_value")
      !DEC$ ATTRIBUTES DLLEXPORT :: get_value

      integer(kind=c_int) :: ierr, n
      character(kind=c_char), intent(in) :: var_name(*)
      real(c_float), dimension(n) :: dest

      character(len=strlen(var_name)) :: f_var_name
      f_var_name = char_array_to_string(var_name, strlen(var_name))

      ierr = ret_code%success

      select case(f_var_name)
       case("z_xz") ! x grid cell centre
         dest = z_xz
       case("z_yz") ! y grid cell centre
         dest = z_yz
       case("zs") ! water level
         dest = zs
       case("zb") ! bed level
         if(subgrid) then
            dest = subgrid_z_zmin
         else
            dest = subgrid_z_zmin
         end if
       case("qsrc_1")
         dest = qsrc(:, 1)
       case("qsrc_2")
         dest = qsrc(:,2)
       case("xsrc")
         dest = xsrc ! TODO this gives segfault
       case("ysrc")
         dest = ysrc ! TODO this gives segfault
       case("tsrc")
         dest = tsrc
       case("zst_bnd")
         dest = zst_bnd
       case default
         write(*,*) 'get_value error'
         ierr = ret_code%failure
      end select

   end function get_value

!-----------------------------------------------------------------------------------------------------!

   function get_value_at_indices(var_name, dest, inds, count) result(ierr) bind(C, name="get_value_at_indices")
      !DEC$ ATTRIBUTES DLLEXPORT :: get_value

      character(kind=c_char), intent(in) :: var_name(*)
      type(c_ptr), intent(inout) :: dest
      type(c_ptr), value, intent(in) :: inds
      integer(kind=c_int), value, intent(in) :: count
      integer(kind=c_int) :: ierr

      character(len=strlen(var_name)) :: f_var_name
      real(kind=c_float), pointer:: f_dest(:)
      integer(kind=c_int), pointer  :: f_inds(:)
      integer :: i

      ierr = ret_code%success

      ! get the name
      f_var_name = char_array_to_string(var_name, strlen(var_name))

      call c_f_pointer(dest, f_dest, [count])
      call c_f_pointer(inds, f_inds, [count])

      select case(f_var_name)
       case("z_xz") ! x grid cell centre
         do i = 1, count
            f_dest(i) = z_xz(f_inds(i))
         end do
       case("z_yz") ! y grid cell centre
         do i = 1, count
            f_dest(i) = z_yz(f_inds(i))
         end do
       case("zs") ! water level
         do i = 1, count
            f_dest(i) = zs(f_inds(i))
         end do
       case("zb") ! bed level
         if(subgrid) then
            do i = 1, count
               f_dest(i) = subgrid_z_zmin(f_inds(i))
            end do
         else
            do i = 1, count
               f_dest(i) = subgrid_z_zmin(f_inds(i))
            end do
         end if
       case("qsrc_1")
         do i = 1, count
            f_dest(i) = qsrc(f_inds(i), 1)
         end do
       case("qsrc_2")
         do i = 1, count
            f_dest(i) = qsrc(f_inds(i), 2)
         end do
       case("xsrc") ! water level
         do i = 1, count
            f_dest(i) = xsrc(f_inds(i))
         end do
       case("ysrc") ! water level
         do i = 1, count
            f_dest(i) = ysrc(f_inds(i))
         end do
       case("tsrc") ! water level
         do i = 1, count
            f_dest(i) = tsrc(f_inds(i))
         end do
       case("zst_bnd")
         do i = 1, count
            f_dest(i) = zst_bnd(f_inds(i))
         end do
       case default
         write(*,*) 'get_value_at_indices error'
         ierr = ret_code%failure
         dest = c_null_ptr
      end select

   end function get_value_at_indices

!-----------------------------------------------------------------------------------------------------!

   function get_value_ptr(var_name, ptr) result(ierr) bind(C, name="get_value_ptr")
      !DEC$ ATTRIBUTES DLLEXPORT :: get_value_ptr

      character(kind=c_char), intent(in) :: var_name(*)
      type(c_ptr), intent(inout) :: ptr
      integer(kind=c_int) :: ierr

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
       case("qsrc_1")
         ptr = c_loc(qsrc(:, 1))
       case("qsrc_2")
         ptr = c_loc(qsrc(:, 2))
       case("xsrc")
         ptr = c_loc(xsrc)
       case("ysrc")
         ptr = c_loc(ysrc)
       case("tsrc")
         ptr = c_loc(tsrc)
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

      ierr = ret_code%success

      f_var_name = char_array_to_string(var_name, strlen(var_name))
      var_shape = (/0/)

      select case(f_var_name)
       case("z_xz", "z_yz","zs","zb","zst_bnd")
         var_shape(1) = np
       case("qsrc_1", "qsrc_2", "xsrc", "ysrc")
         var_shape(1) = nsrc
       case("tsrc")
         var_shape(1) = ntsrc
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

      ierr = ret_code%success
      select case(f_var_name)
       case("z_xz", "z_yz", "zs", "zb" ,"qsrc_1", "qsrc_2", "xsrc", "ysrc", "tsrc", "zst_bnd")
         f_var_type = "float"
         var_type = string_to_char_array(trim(f_var_type), len(trim(f_var_type)))
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
       case("z_xz", "z_yz", "zs", "zb" ,"qsrc_1", "qsrc_2", "xsrc", "ysrc", "tsrc", "zst_bnd")
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

      select case(f_var_name)
       case("zs")
         call c_f_pointer(ptr, f_ptr, [np])
         do i = 1, np
            zs(i) = f_ptr(i)
         end do
       case("zb")
         call c_f_pointer(ptr, f_ptr, [np])
         do i = 1, np
            zb(i) = f_ptr(i)
         end do
       case("qsrc_1")
         call c_f_pointer(ptr, f_ptr, [nsrc])
         do i = 1, nsrc
            qsrc(i, 1) = f_ptr(i)
         end do
       case("qsrc_2")
         call c_f_pointer(ptr, f_ptr, [nsrc])
         do i = 1, nsrc
            qsrc(i, 2) = f_ptr(i)
         end do
       case("tsrc")
         call c_f_pointer(ptr, f_ptr, [ntsrc])
         do i = 1, ntsrc
            tsrc(i) = f_ptr(i)
         end do
       case("zst_bnd")
         call c_f_pointer(ptr, f_ptr, [np])
         do i = 1, np
            zst_bnd(i) = f_ptr(i)
         end do
       case default
         write(*,*) 'set_value error'
         ierr = ret_code%failure
      end select

   end function set_value

!-----------------------------------------------------------------------------------------------------!

   function set_value_at_indices(var_name, inds, count, src) result(ierr) bind(C, name="set_value_at_indices")
      !DEC$ ATTRIBUTES DLLEXPORT :: set_value_at_indices

      character(kind=c_char), intent(in) :: var_name(*)
      type(c_ptr), value, intent(in) :: inds
      type(c_ptr), value, intent(in) :: src
      integer(kind=c_int), value, intent(in) :: count
      integer(kind=c_int) :: ierr


      character(len=strlen(var_name)) :: f_var_name
      integer(kind=c_int), pointer  :: f_inds(:)
      real(kind=c_float), pointer  :: f_src(:)
      integer :: i

      ierr = ret_code%success

      f_var_name = char_array_to_string(var_name, strlen(var_name))

      call c_f_pointer(src, f_src, [count])
      call c_f_pointer(inds, f_inds, [count])

      select case(f_var_name)
       case("zs")
         do i = 1, count
            zs(f_inds(i)) = f_src(i)
         end do
       case("zb")
         do i = 1, count
            zb(f_inds(i)) = f_src(i)
         end do
       case("qsrc_1")
         do i = 1, count
            qsrc(f_inds(i), 1) = f_src(i)
         end do
       case("qsrc_2")
         do i = 1, count
            qsrc(f_inds(i), 2) = f_src(i)
         end do
       case("tsrc")
         do i = 1, count
            tsrc(f_inds(i)) = f_src(i)
         end do
       case("zst_bnd")
         do i = 1, count
            zst_bnd(f_inds(i)) = f_src(i)
         end do
       case default
         write(*,*) 'set_value_at_indices error'
         ierr = ret_code%failure
      end select

   end function set_value_at_indices

!-----------------------------------------------------------------------------------------------------!

   function get_grid_type(grid_type) result(ierr) bind(C, name="get_grid_type")
      !DEC$ ATTRIBUTES DLLEXPORT :: get_grid_rtype

      character(kind=c_char), intent(out) :: grid_type(maxstrlen)
      integer(kind=c_int) :: ierr

      character(len=maxstrlen) :: string

      if(use_quadtree) then
         string = "unstructered"
      else
        ! Note: even if input is regular mesh, internal data structure is quadtree(?),
         string = "unstructured"
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
