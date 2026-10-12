module sfincs_error
   !
   ! Error handling: stopping SFINCS with a message, and checking that input files exist.
   !
   ! Subroutines / functions:
   !
   !   stop_sfincs(message, error_code)
   !     Write an error message to the log and stop SFINCS with the given error code.
   !     Called from: check_file_exists (sfincs_error) and the input and initialization
   !       routines of many modules
   !
   !   check_file_exists(file_name, file_type, stop_on_error)
   !     Check whether a file exists; optionally stop SFINCS when it does not.
   !     Called from: the input-reading routines of nearly all modules
   !
   use sfincs_log
   !
   contains
   !
   subroutine stop_sfincs(message, error_code)
   !
   ! Write an error message to the log and stop SFINCS with the given error code.
   !
   ! Input:
   !   message    : error message to write to the log
   !   error_code : exit code
   !
   ! Called from: check_file_exists (sfincs_error) and the input and initialization routines
   !   of many modules
   !
   implicit none
   !
   integer, intent(in)          :: error_code
   character(len=*), intent(in) :: message
   !
   write(logstr,'(a, a)')'Error   : ', message
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
   ! Input:
   !   file_name     : name of the file
   !   file_type     : description of the file, used in the error message
   !   stop_on_error : stop SFINCS when the file does not exist
   !
   ! Returns:
   !   exists        : true when the file exists
   !
   ! Called from: the input-reading routines of nearly all modules
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
