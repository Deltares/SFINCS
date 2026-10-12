module sfincs_error
   !
   ! Headers auto-generated, ask main developers in case of doubts
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
   ! All indications and logos of, and references to, "Delft3D"
   ! and "Deltares" are registered trademarks of Stichting Deltares,
   ! and remain the property of Stichting Deltares. All rights reserved.
   !
   ! ------------ Copyright notice SFINCS ------------
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
