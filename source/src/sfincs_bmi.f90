module sfincs_bmi
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! Basic Model Interface (BMI) of SFINCS. C-bound functions that let an external program
   ! initialize, step and finalize the model, and get pointers to (and metadata of) its
   ! variables.
   !
   ! Subroutines / functions:
   !
   !   initialize()
   !     Initialize the model in BMI mode.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   update()
   !     Advance the model by a single time step.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   update_until(t_target)
   !     Advance the model over the time interval t_target.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   finalize()
   !     Finalize the model.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   get_value_ptr(c_var_name, c_data)
   !     Return a C pointer to the model variable with the given name.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   get_var_shape(c_var_name, var_shape)
   !     Return the shape of the model variable with the given name.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   get_var_type(c_var_name, c_type)
   !     Return the data type of the model variable with the given name.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   get_var_rank(c_var_name, rank)
   !     Return the rank of the model variable with the given name.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   set_logical(c_flag_name, ival)
   !     Set a logical model flag by name (e.g. switch the external source qext on or off).
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   get_start_time(tstart)
   !     Return the simulation start time.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   get_end_time(tend)
   !     Return the simulation end time.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   get_time_step(deltat)
   !     Return the current time step.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   get_current_time(tcurrent)
   !     Return the current simulation time.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   bmi_update_bed_level()
   !     Apply the external bed level change dzbext (C name update_bed_level).
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   update_apparent_roughness()
   !     Update the apparent (wave-enhanced) roughness at the uv points.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   get_sfincs_cell_index(x, y, indx)
   !     Return the (1-based) index of the SFINCS cell that contains point (x, y).
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   get_sfincs_cell_indices(x, y, indx, n)
   !     Return the (1-based) SFINCS cell indices for n points (x, y).
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   get_sfincs_cell_area(indx, area)
   !     Return the area of the cell with the given index.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   get_last_bmi_error(c_error)
   !     Return the last BMI error message.
   !     Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !   strlen(char_array, max_len)
   !     Length of a null-terminated C character array.
   !     Called from: get_value_ptr (sfincs_bmi), get_var_rank (sfincs_bmi),
   !       get_var_shape (sfincs_bmi), get_var_type (sfincs_bmi), set_logical (sfincs_bmi)
   !
   !   char_array_to_string(char_array, length)
   !     Convert a C character array to a Fortran string.
   !     Called from: get_value_ptr (sfincs_bmi), get_var_rank (sfincs_bmi),
   !       get_var_shape (sfincs_bmi), get_var_type (sfincs_bmi), set_logical (sfincs_bmi)
   !
   !   string_to_char_array(string, length)
   !     Convert a Fortran string to a null-terminated C character array.
   !     Called from: get_last_bmi_error (sfincs_bmi), get_var_type (sfincs_bmi)
   !
   ! ------------ Copyright notice SFINCS ------------
   !
   ! Copyright (C)  Stichting Deltares, 2016-2026.
   !
   ! This file is part of SFINCS.
   !
   ! SFINCS is free software: you can redistribute it and/or modify
   ! it under the terms of the GNU General Public License as
   ! published by the Free Software Foundation version 3.
   !
   ! SFINCS is distributed in the hope that it will be useful,
   ! but WITHOUT ANY WARRANTY; without even the implied warranty of
   ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   ! GNU General Public License for more details.
   !
   ! You should have received a copy of the GNU General Public License
   ! along with SFINCS.  If not, see <http://www.gnu.org/licenses/>.
   !
   ! contact: software@deltares.nl
   ! Stichting Deltares
   ! P.O. Box 177
   ! 2600 MH Delft, The Netherlands
   !
   ! All indications, logos of, and references to, "SFINCS"
   ! remain the property of Stichting Deltares. All rights reserved.
   !
   ! ------------ Copyright notice SFINCS ------------
   !
   use iso_c_binding
   use sfincs_lib
   use sfincs_data
   use sfincs_domain
   use sfincs_wave_enhanced_roughness
   !   
   implicit none
   private
   !
   ! routines
   !
   public :: initialize
   public :: update
   public :: update_until
   public :: finalize
   public :: get_value_ptr
   public :: get_var_shape
   public :: get_var_type
   public :: get_var_rank
   public :: get_start_time
   public :: get_end_time
   public :: get_time_step
   public :: get_current_time
   public :: bmi_update_bed_level
   public :: update_apparent_roughness
   public :: get_sfincs_cell_index
   public :: get_sfincs_cell_indices
   public :: get_sfincs_cell_area
   !
   ! constants
   !
   public :: BMI_LENVARADDRESS
   public :: BMI_LENVARTYPE
   public :: BMI_LENGRIDTYPE
   public :: BMI_LENCOMPONENTNAME
   public :: BMI_LENVERSION
   public :: BMI_LENERRMESSAGE

   integer(c_int), bind(C, name="BMI_LENVARADDRESS") :: BMI_LENVARADDRESS = 64 !< max. length for the variable's address C-string, also the variable name
   !DIR$ ATTRIBUTES DLLEXPORT :: BMI_LENVARADDRESS
   integer(c_int), bind(C, name="BMI_LENVARTYPE") :: BMI_LENVARTYPE = 64 !< max. length for variable type C-strings
   !DIR$ ATTRIBUTES DLLEXPORT :: BMI_LENVARTYPE
   integer(c_int), bind(C, name="BMI_LENGRIDTYPE") :: BMI_LENGRIDTYPE = 64 !< max. length for grid type C-strings
   !DIR$ ATTRIBUTES DLLEXPORT :: BMI_LENGRIDTYPE
   integer(c_int), bind(C, name="BMI_LENCOMPONENTNAME") :: BMI_LENCOMPONENTNAME = 64 !< component name length, i.e. 'SFINCS'
   !DIR$ ATTRIBUTES DLLEXPORT :: BMI_LENCOMPONENTNAME
   integer(c_int), bind(C, name="BMI_LENVERSION") :: BMI_LENVERSION = 256 !< length of version string, e.g. '3.1' or '16.4 release candidate'
   !DIR$ ATTRIBUTES DLLEXPORT :: BMI_LENVERSION
   integer(c_int), bind(C, name="BMI_LENERRMESSAGE") :: BMI_LENERRMESSAGE = 1024 !< max. length of error message
   !DIR$ ATTRIBUTES DLLEXPORT :: BMI_LENERRMESSAGE

contains
   
   function initialize() result(ierr) bind(C, name="initialize")
   !DEC$ ATTRIBUTES DLLEXPORT :: initialize
      !
      ! Initialize the model in BMI mode.
      !
      ! Returns:
      !   ierr : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !
      integer(kind=c_int) :: ierr
      bmi  = .true.
      use_qext = .true.
      ierr = sfincs_initialize()

   end function initialize
   
   function update() result(ierr) bind(C, name="update")
   !DEC$ ATTRIBUTES DLLEXPORT :: update
      !
      ! Advance the model by a single time step.
      !
      ! Returns:
      !   ierr : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !
      integer(kind=c_int) :: ierr
      real*8 :: unused_t
      ! Run single time step      
      unused_t = -1.0
      ierr = sfincs_update(unused_t)
   
   end function update
   
   function update_until(t_target) result(ierr) bind(C, name="update_until")
   !DEC$ ATTRIBUTES DLLEXPORT :: update_until  
      !
      ! Advance the model over the time interval t_target.
      !
      ! Input:
      !   t_target : length of the time interval to run [s]
      !
      ! Returns:
      !   ierr     : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !
      real(kind=c_double), value, intent(in) :: t_target
      integer(kind=c_int) :: ierr

      real(kind=c_double) :: delta_t
      
      delta_t = t_target - t
      ierr = sfincs_update(delta_t)
      
   end function update_until
   
   function finalize() result(ierr) bind(C, name="finalize")
   !DEC$ ATTRIBUTES DLLEXPORT :: finalize
   !
   ! Finalize the model.
   !
   ! Returns:
   !   ierr : error code (0 = success)
   !
   ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   integer(kind=c_int) :: ierr
   
   ierr = sfincs_finalize()
   
   end function finalize

   function get_value_ptr(c_var_name, c_data) result(ierr) &
      bind(C, name="get_value_ptr")
   !DIR$ ATTRIBUTES DLLEXPORT :: get_value_ptr
      !
      ! Return a C pointer to the model variable with the given name.
      !
      ! Input:
      !   c_var_name : name of the variable (C string)
      !
      ! Input/output:
      !   c_data     : C pointer to the variable
      !
      ! Returns:
      !   ierr       : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !
      character(kind=c_char), intent(in) :: c_var_name(*)
      type(c_ptr), intent(inout) :: c_data
      integer(kind=c_int) :: ierr
      character(len=strlen(c_var_name, BMI_LENVARADDRESS)) :: var_name
      integer :: c_strlen

      ierr = 0

      c_strlen = strlen(c_var_name, BMI_LENVARADDRESS)
      var_name = char_array_to_string(c_var_name, c_strlen)

      select case(var_name)
      case("z_xz")
         c_data = c_loc(z_xz)
      case("z_yz")
         c_data = c_loc(z_yz)
      case("zs")
         c_data = c_loc(zs)
      case("zb")
         c_data = c_loc(zb)
      case("subgrid_z_zmin")
         c_data = c_loc(subgrid_z_zmin)
      case("qext")
         !
         ! qext(np): external source/sink per cell in m3/s (positive = water
         ! added to the cell). It is added to qsrc in update_continuity.
         ! NOTE: before the qsrc restructuring qext was a rate in m/s; couplers
         ! that still provide m/s must now multiply by the cell area.
         !
         c_data = c_loc(qext)
      case("dzbext")
         if (.not. allocated(dzbext)) then
            allocate(dzbext(np))
            dzbext = 0.0
         endif
         c_data = c_loc(dzbext)
      case("uorb")
         c_data = c_loc(uorb)
      case default
         c_data = c_null_ptr
         ierr = -1
      end select

   end function get_value_ptr
   
   function get_var_shape(c_var_name, var_shape) result(ierr) &
      bind(C, name="get_var_shape")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var_shape
      !
      ! Return the shape of the model variable with the given name.
      !
      ! Input:
      !   c_var_name : name of the variable (C string)
      !
      ! Input/output:
      !   var_shape  : shape of the variable
      !
      ! Returns:
      !   ierr       : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !
      character(kind=c_char), intent(in) :: c_var_name(*)
      integer(c_int), intent(inout) :: var_shape(*)
      character(len=strlen(c_var_name, BMI_LENVARADDRESS)) :: var_name
      integer(kind=c_int) :: ierr
      !
      ! The caller (xmipy) allocates var_shape with exactly get_var_rank entries,
      ! which is 1 for every variable exposed here. Only ever write var_shape(1):
      ! writing a fixed 6-element shape overruns the 4-byte buffer and corrupts
      ! the heap (intermittent 0xc0000374 in update/finalize).
      !
      ierr = 0

      var_name = char_array_to_string(c_var_name, strlen(c_var_name, BMI_LENVARADDRESS))
      var_shape(1) = 0

      select case(var_name)
      case("z_xz", "z_yz", "zs", "zb", "uorb")
         var_shape(1) = size(zs)
      case("subgrid_z_zmin")
         var_shape(1) = size(subgrid_z_zmin)
      case("z_index_z_n", "z_index_z_m")
         var_shape(1) = size(z_index_z_n)
      case("qext")
         var_shape(1) = size(qext)
      case("dzbext")
         if (.not. allocated(dzbext)) then
            allocate(dzbext(np))
            dzbext = 0.0
         endif
         var_shape(1) = size(dzbext)
      case default
         ierr = -1
      end select
   
   end function get_var_shape
   
   function get_var_type(c_var_name, c_type) result(ierr) &
      bind(C, name="get_var_type")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var_type   
      !
      ! Return the data type of the model variable with the given name.
      !
      ! Input:
      !   c_var_name : name of the variable (C string)
      !
      ! Output:
      !   c_type     : data type of the variable (C string)
      !
      ! Returns:
      !   ierr       : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !
      character(kind=c_char), intent(in) :: c_var_name(*)
      character(kind=c_char), intent(out) ::  c_type(BMI_LENVARTYPE)
      character(len=BMI_LENVARADDRESS-1) :: var_name
      character(len=BMI_LENVARTYPE-1) :: type_name

      integer(kind=c_int)              :: ierr

      ierr = 0
      var_name = char_array_to_string(c_var_name, strlen(c_var_name, BMI_LENVARADDRESS))

      select case(var_name)
      case("z_xz", "z_yz", "zb", "subgrid_z_zmin", "qext", "dzbext", "uorb")
         type_name = "float"
      case("zs")
         type_name = "double"
      case("z_index_z_n", "z_index_z_m")
         type_name = "integer"
      case default
         ierr = -1
      end select
      
      c_type = string_to_char_array(trim(type_name), len(trim(type_name)))
   
   end function get_var_type
   
   function get_var_rank(c_var_name, rank) result(ierr) &
      bind(C, name="get_var_rank")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_var_rank
      !
      ! Return the rank of the model variable with the given name.
      !
      ! Input:
      !   c_var_name : name of the variable (C string)
      !
      ! Output:
      !   rank       : rank of the variable
      !
      ! Returns:
      !   ierr       : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !
      character(kind=c_char), intent(in) :: c_var_name(*)
      integer(c_int), intent(out) :: rank
      integer(kind=c_int) :: ierr

      character(len=BMI_LENVARADDRESS-1) :: var_name
      
      ierr = 0

      var_name = char_array_to_string(c_var_name, strlen(c_var_name, BMI_LENVARADDRESS))
      
      select case(var_name)
      case("z_xz", "z_yz", "zs", "zb", "subgrid_z_zmin", "qext", "dzbext", "uorb")
         rank = 1
      case default
         ierr = -1
      end select

   end function get_var_rank

   function set_logical(c_flag_name, ival) result(ierr) bind(C, name="set_logical")
   !DEC$ ATTRIBUTES DLLEXPORT :: set_logical
      !
      ! Set a logical model flag by name (e.g. switch the external source qext on or off).
      !
      ! Input:
      !   c_flag_name : name of the flag (C string)
      !   ival        : value to set: 0 is false, anything else is true
      !
      ! Returns:
      !   ierr        : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !
      character(kind=c_char), intent(in) :: c_flag_name(*)
      integer(kind=c_int) :: ival(*)
      integer(kind=c_int) :: ierr
      character(len=BMI_LENVARADDRESS-1) :: flag_name
      logical :: bval
      !
      ierr = 0
      !
      flag_name = char_array_to_string(c_flag_name, strlen(c_flag_name, BMI_LENVARADDRESS))      
      !
      if (ival(1) == 0) then
         bval = .false.
      else
         bval = .true.
      endif
      !
      select case(flag_name)
      case("qext")
         !
         ! Switch the external source/sink qext (m3/s per cell) on or off
         !
         use_qext = bval
         !write(*,*)'use_qext = ', use_qext
      case("dzbext")
         !
         ! Lazily allocate the external delta-bed-level array on first enable.
         ! Once allocated we keep it around; toggling the flag off later just
         ! disables the update path without freeing memory (same pattern as
         ! qext is handled elsewhere).
         !
         if (bval .and. .not. allocated(dzbext)) then
            !
            allocate(dzbext(np))
            dzbext = 0.0
            !
         endif
         !
         use_dzbext = bval
         !
      case default
         ierr = -1
      end select

   end function set_logical
      
   function get_start_time(tstart) result(ierr) bind(C, name="get_start_time")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_start_time
      !
      ! Return the simulation start time.
      !
      ! Output:
      !   tstart : simulation start time [s]
      !
      ! Returns:
      !   ierr   : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !
      real(c_double), intent(out) :: tstart
      integer(kind=c_int) :: ierr

      tstart = t0
      ierr = 0

   end function get_start_time
      
   function get_end_time(tend) result(ierr) bind(C, name="get_end_time")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_end_time   
      !
      ! Return the simulation end time.
      !
      ! Output:
      !   tend : simulation end time [s]
      !
      ! Returns:
      !   ierr : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !
      real(c_double), intent(out) :: tend
      integer(kind=c_int) :: ierr

      tend = t1
      ierr = 0

   end function get_end_time
      
   function get_time_step(deltat) result(ierr) bind(C, name="get_time_step")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_time_step
      !
      ! Return the current time step.
      !
      ! Output:
      !   deltat : current time step [s]
      !
      ! Returns:
      !   ierr   : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !
      real(c_double), intent(out) :: deltat
      integer(kind=c_int) :: ierr

      deltat = dt
      ierr = 0

   end function get_time_step
      
   function get_current_time(tcurrent) result(ierr) bind(C, name="get_current_time")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_current_time
      !
      ! Return the current simulation time.
      !
      ! Output:
      !   tcurrent : current simulation time [s]
      !
      ! Returns:
      !   ierr     : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !
      real(c_double), intent(out) :: tcurrent
      integer(kind=c_int) :: ierr

      tcurrent = t
      ierr = 0

   end function get_current_time

   function bmi_update_bed_level() result(ierr) bind(C, name="update_bed_level")
   ! Apply dzbext to the bed-level arrays (zb or subgrid_z_zmin/zmax and
   ! subgrid_uv_zmin/zmax), rebuild zbuvmx for non-subgrid runs, lift the water
   ! surface over the moving cells and refresh the bed level at observation points.
   !
   ! Returns:
   !   ierr : error code (0 = success)
   !
   ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !DEC$ ATTRIBUTES DLLEXPORT :: bmi_update_bed_level
      integer(kind=c_int) :: ierr
      call update_bed_level()
      ierr = 0

   end function bmi_update_bed_level
   
   function update_apparent_roughness() result(ierr) bind(C, name="update_apparent_roughness")
   ! Update apparent roughness at uv points
   !
   ! Returns:
   !   ierr : error code (0 = success)
   !
   ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !DEC$ ATTRIBUTES DLLEXPORT :: update_apparent_roughness
      integer(kind=c_int) :: ierr
      !
      ! This should be called after uorb has been updated
      !      
      call update_wave_enhanced_roughness() ! sits in sfincs_wave_enhanced_roughness.f90
      !
      ierr = 0
      !    
   end function update_apparent_roughness
   
   function get_sfincs_cell_index(x, y, indx) result(ierr) bind(C, name="get_sfincs_cell_index")
   ! Return (1-based) cell index of SFINCS domain
   !
   ! Input:
   !   x    : x coordinate(s) of the point(s)
   !   y    : y coordinate(s) of the point(s)
   !
   ! Output:
   !   indx : index of the cell that contains the point (1-based)
   !
   ! Returns:
   !   ierr : error code (0 = success)
   !
   ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !DEC$ ATTRIBUTES DLLEXPORT :: get_sfincs_cell_index
      use quadtree
      !         
      real(c_double), intent(in)  :: x
      real(c_double), intent(in)  :: y
      integer(c_int), intent(out) :: indx
      integer(kind=c_int)         :: ierr
      integer                     :: nmq
      !
      !
      nmq = find_quadtree_cell(real(x, kind=4), real(y, kind=4))
      indx = index_sfincs_in_quadtree(nmq)
      !
      ierr = 0
      !
   end function get_sfincs_cell_index

   function get_sfincs_cell_indices(x, y, indx, n) result(ierr) bind(C, name="get_sfincs_cell_indices")
   ! Return (1-based) cell index of SFINCS domain
   !
   ! Input:
   !   x    : x coordinate(s) of the point(s)
   !   y    : y coordinate(s) of the point(s)
   !   n    : number of points
   !
   ! Output:
   !   indx : index of the cell that contains each point (1-based)
   !
   ! Returns:
   !   ierr : error code (0 = success)
   !
   ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
   !
   !DEC$ ATTRIBUTES DLLEXPORT :: get_sfincs_cell_indices
      use quadtree
      !         
      real(c_double), intent(in)  :: x(n)
      real(c_double), intent(in)  :: y(n)
      integer(c_int), intent(out) :: indx(n)
      integer(c_int), value       :: n
      integer(kind=c_int)         :: ierr
      integer                     :: nmq
      integer                     :: iq
      !
      do iq = 1, n
         !
         nmq = find_quadtree_cell(real(x(iq), kind=4), real(y(iq), kind=4))
         !
         if (nmq > 0.0) then
            !
            indx(iq) = index_sfincs_in_quadtree(nmq)
            !
         else
            !
            indx(iq) = 0
            !
         endif   
         !
      enddo
      !
      ierr = 0
      !
   end function get_sfincs_cell_indices
   
   function get_sfincs_cell_area(indx, area) result(ierr) bind(C, name="get_sfincs_cell_area")
   !DEC$ ATTRIBUTES DLLEXPORT :: get_sfincs_cell_area
      !
      ! Return the area of the cell with the given index.
      !
      ! Input:
      !   indx : cell index (1-based)
      !
      ! Output:
      !   area : area of the cell [m2]
      !
      ! Returns:
      !   ierr : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !         
      integer(kind=c_int), intent(in) :: indx
      real(c_double), intent(out)     :: area
      integer(kind=c_int)             :: ierr
      real*4                          :: a
      !
      if (crsgeo) then
         !
         a = cell_area_m2(indx)
         !
      else   
         !
         a = cell_area(z_flags_iref(indx))
         !
      endif
      !
      area = real(a, kind=8)
      !
      ierr = 0
      !
   end function get_sfincs_cell_area
   
   !> @brief Get the last error in the BMI as a character array
   !! with size BMI_LENERRMESSAGE
   !<
   function get_last_bmi_error(c_error) result(ierr) &
      bind(C, name="get_last_bmi_error")
   !DIR$ ATTRIBUTES DLLEXPORT :: get_last_bmi_error
      !
      ! Return the last BMI error message.
      !
      ! Output:
      !   c_error : last error message (C string)
      !
      ! Returns:
      !   ierr    : error code (0 = success)
      !
      ! Called from: external BMI caller, through the C interface (e.g. a model coupler)
      !
      character(kind=c_char, len=1), intent(out) :: c_error(BMI_LENERRMESSAGE)
      integer(kind=c_int) :: ierr

      c_error = string_to_char_array("error handling not implemented", 30)
      ierr = 0

   end function get_last_bmi_error
   
   !> @brief Returns the string length without the trailing null character
   !<
   pure function strlen(char_array, max_len) result(string_length)
      !
      ! Length of a null-terminated C character array.
      !
      ! Input:
      !   char_array    : C character array
      !   max_len       : maximum length to search
      !
      ! Returns:
      !   string_length : length up to the terminating null character
      !
      ! Called from: get_value_ptr (sfincs_bmi), get_var_rank (sfincs_bmi),
      !   get_var_shape (sfincs_bmi), get_var_type (sfincs_bmi), set_logical (sfincs_bmi)
      !
      integer, intent(in) :: max_len
      character(c_char), intent(in) :: char_array(max_len) !< C-style character string
      integer :: string_length !< Fortran string length
      integer :: i

      string_length = 0
      do i = 1, size(char_array)
         if (char_array(i) .eq. C_NULL_CHAR) then
            string_length = i - 1
            exit
         end if
      end do

   end function strlen

   !> @brief Convert C-style string to Fortran character string
   !<
   pure function char_array_to_string(char_array, length) result(f_string)
      !
      ! Convert a C character array to a Fortran string.
      !
      ! Input:
      !   char_array : C character array to convert
      !   length     : string length without the terminating null character
      !
      ! Returns:
      !   f_string   : Fortran string
      !
      ! Called from: get_value_ptr (sfincs_bmi), get_var_rank (sfincs_bmi),
      !   get_var_shape (sfincs_bmi), get_var_type (sfincs_bmi), set_logical (sfincs_bmi)
      !
      integer(c_int), intent(in) :: length !< string length without terminating null character
      character(c_char), intent(in) :: char_array(length) !< string to convert
      character(len=length) :: f_string !< Fortran fixed length character string
      integer :: i

      do i = 1, length
         f_string(i:i) = char_array(i)
      end do

   end function char_array_to_string

   !> @brief Convert Fortran string to C-style character string
   !<
   pure function string_to_char_array(string, length) result(c_array)
      !
      ! Convert a Fortran string to a null-terminated C character array.
      !
      ! Input:
      !   string  : Fortran string to convert
      !   length  : length of the Fortran string
      !
      ! Returns:
      !   c_array : null-terminated C character array
      !
      ! Called from: get_last_bmi_error (sfincs_bmi), get_var_type (sfincs_bmi)
      !
      integer(c_int), intent(in) :: length !< Fortran string length
      character(len=length), intent(in) :: string !< string to convert
      character(kind=c_char, len=1) :: c_array(length + 1) !< C-style character string
      integer :: i

      do i = 1, length
         c_array(i) = string(i:i)
      end do
      c_array(length + 1) = C_NULL_CHAR

   end function string_to_char_array

end module sfincs_bmi