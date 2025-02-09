module sfincs_log
   !
   integer        :: fid
   character(256) :: logstr
   !
contains   

   subroutine open_log()
   !
   implicit none
   !
   fid = 777
   open(unit = fid, file = 'sfincs.log')
   !
   end subroutine

   subroutine write_log(str, to_screen)
   !
   implicit none
   !
   character(*), intent(in) :: str
   integer, intent(in)      :: to_screen
   !
   write(fid,'(a)')trim(str)
   !
   if (to_screen==1) then
      write(*,'(a)')trim(str)
   endif
   !
   end subroutine

   subroutine close_log()
   !
   implicit none
   !
   close(fid)
   !
   end subroutine

end module
