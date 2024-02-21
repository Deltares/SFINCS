module sfincs_error

contains
   ! 
   function check_file_exists(file_name, file_type) result(exists)
   !
   ! Checks if file exists
   ! Returns true or false 
   !
   use sfincs_data
   !
   implicit none
   !
   character(len=*) :: file_name
   character(len=*) :: file_type
   logical          :: exists
   !
   exists = .true. 
   ! 
   inquire( file=trim(file_name), exist=exists )
   ! 
   if (.not. exists) then
      !
      error = 2
	  ! 
      write(error_message,'(5a)')'Error! ', trim(file_type), ' "', trim(file_name), '" not found!' 
      !
   endif    
   !
   end function
   !
end module
