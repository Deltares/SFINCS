module sfincs_error
   !
   use sfincs_log
   !
   contains
   !
   subroutine stop_sfincs(message, error_code)
   !
   implicit none
   !
   integer, intent(in)          :: error_code
   character(len=*), intent(in) :: message
   !
   logstr = 'Error   : ' // trim(message)
   call write_log(logstr, 1)
   stop error_code
   !   
   end
   ! 
   function check_file_exists(file_name, file_type, stop_on_error) result(exists)
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
   character(256)   :: message
   logical          :: exists
   logical          :: stop_on_error
   !
   exists = .true. 
   ! 
   inquire( file=trim(file_name), exist=exists)
   ! 
   if (.not. exists) then
      !
      if (stop_on_error) then
         !
         write(message,'(4a)')trim(file_type), ' "', trim(file_name), '" not found! SFINCS has stopped!' 
         call stop_sfincs(trim(message), 2)
         !
      endif   
      !
   endif    
   !
   end function
   !
end module
