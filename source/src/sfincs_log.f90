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

   function fmt_real(val, decimals) result(s)
      !
      ! Format a real with minimum width and a guaranteed leading zero
      ! for |val| < 1. ifx's "f0.d" descriptor drops the leading zero in
      ! that range, which is not standard-conforming; this helper rewrites
      ! the result so the log output always reads "0.6670" rather than
      ! ".6670".
      !
      ! Called from: write_src_structures_log_summary (sfincs_src_structures),
      ! urban_drainage log summary (sfincs_urban_drainage), and anywhere
      ! else a real needs to be embedded in a log line with the smallest
      ! reasonable field width.
      !
      implicit none
      !
      real,    intent(in) :: val
      integer, intent(in) :: decimals
      character(len=32)   :: s
      !
      character(len=16)   :: fmt
      !
      write(fmt,'(a,i0,a)') '(f0.', decimals, ')'
      write(s,fmt) val
      s = adjustl(s)
      !
      if (s(1:1) == '.') then
         !
         s = '0' // s(1:len_trim(s))
         !
      else if (s(1:2) == '-.') then
         !
         s = '-0' // trim(s(2:))
         !
      endif
      !
   end function fmt_real

end module
