module sfincs_read

contains
   
   subroutine read_real_input(fileid,keyword,value,default)
   !
   character(*), intent(in) :: keyword
   character(len=256)       :: keystr
   character(len=256)       :: valstr
   character(len=256)       :: line
   integer, intent(in)      :: fileid
   real*4, intent(out)      :: value
   real*4, intent(in)       :: default
   integer j,stat,ilen
   !
   value = default
   !
   rewind(fileid)   
   !
   do while(.true.)
      !
      read(fileid,'(a)',iostat = stat)line
      !
      if (stat==-1) exit
      !
      call read_line(line, keystr, valstr)
      !
      if (trim(keystr)==trim(keyword)) then
         !
         read(valstr,*)value         
         !
         exit
         !
      endif
      !
   enddo 
   !
   end  subroutine  

   subroutine read_real_array_input(fileid,keyword,value,default,nr)
   !
   character(*), intent(in) :: keyword
   character(len=256)       :: keystr
   character(len=256)       :: valstr
   character(len=256)       :: line
   integer, intent(in)      :: fileid
   integer, intent(in)      :: nr
   real*4, dimension(:), intent(out), allocatable :: value
   real*4, intent(in)       :: default
   integer j,stat, m,ilen
   !
   allocate(value(nr))
   !
   value = default
   !
   rewind(fileid)   
   !
   do while(.true.)
      !
      read(fileid,'(a)',iostat = stat)line
      !
      if (stat==-1) exit
      !
      call read_line(line, keystr, valstr)
      !
      if (trim(keystr)==trim(keyword)) then
         !
         read(valstr,*)(value(m), m = 1, nr)
         !
         exit
         !
      endif
      !
   enddo 
   !
   end  subroutine  

   
   subroutine read_int_input(fileid,keyword,value,default)
   !
   character(*), intent(in) :: keyword
   character(len=256)       :: keystr
   character(len=256)       :: valstr
   character(len=256)       :: line
   integer, intent(in)      :: fileid
   integer, intent(out)     :: value
   integer, intent(in)      :: default
   integer j,stat,ilen
   !
   value = default
   !
   rewind(fileid)   
   !
   do while(.true.)
      !
      read(fileid,'(a)',iostat = stat)line
      !
      if (stat==-1) exit
      !
      call read_line(line, keystr, valstr)
      !
      if (trim(keystr)==trim(keyword)) then
         !
         read(valstr,*)value         
         !
         exit
         !
      endif
      !
   enddo 
   !
   end subroutine

   
   subroutine read_char_input(fileid,keyword,value,default)
   !
   character(*), intent(in)  :: keyword
   character(len=256)        :: keystr0
   character(len=256)        :: keystr
   character(len=256)        :: valstr
   character(len=256)        :: line
   integer, intent(in)       :: fileid
   character(*), intent(in)  :: default
   character(*), intent(out) :: value
   integer j,stat,ilen,jn
   !
   value = default
   !
   rewind(fileid)   
   !
   do while(.true.)
      !
      read(fileid,'(a)',iostat = stat)line
      !
      if (stat==-1) exit
      !
      call read_line(line, keystr, valstr)
      !
      if (trim(keystr)==trim(keyword)) then
         !
         value = valstr
         !
         exit
         !
      endif
      !
   enddo 
   !
   end subroutine 


   subroutine read_logical_input(fileid,keyword,value,default)
   !
   character(*), intent(in)  :: keyword
   character(len=256)        :: keystr0
   character(len=256)        :: keystr
   character(len=256)        :: valstr
   character(len=256)        :: line
   integer, intent(in)       :: fileid
   logical, intent(in)       :: default
   logical, intent(out)      :: value
   integer j,stat,ilen
   !
   value = default
   !
   rewind(fileid)   
   !
   do while(.true.)
      !
      read(fileid,'(a)',iostat = stat)line
      !
      if (stat==-1) exit
      !
      call read_line(line, keystr, valstr)
      !
      if (trim(keystr)==trim(keyword)) then
         !
         if (valstr(1:1) == '1' .or. valstr(1:1) == 'y' .or. valstr(1:1) == 'Y' .or. valstr(1:1) == 't' .or. valstr(1:1) == 'T') then
            value = .true.
         else
            value = .false.
         endif 
         !
         exit
         !
      endif
      !
   enddo 
   !
   end subroutine 

   subroutine read_line(line0, keystr, valstr)
   !
   ! Reads line from input file, returns keyword and value strings
   !
   character(*), intent(in)  :: line0
   character(len=256)        :: line
   character(*), intent(out) :: keystr
   character(*), intent(out) :: valstr
   integer j, ilen, jn
   !
   keystr = ''
   valstr = '' 
   !
   ! Change tabs into spaces.
   !
   call notabs(line0, line, ilen)
   !
   ! Look for line ending character. Remove it if it exists.
   !
   jn = index(line, '\r')      
   !
   if (jn > 0) then
      !
      ! New line character detected (probably sfincs.inp with windows line endings, running in linux)
      !
      line = line(1 : jn - 1)            
      ! 
   endif
   !
   ! Remove leading and trailing spaces.
   !
   line = trim(line)
   !
   if (line(1:1) == '#' .or. line(1:1) == '!' .or. line(1:1) == '@') return
   !
   ! Find "="
   !
   j  = index(line, '=')
   !
   if (j == 0) return
   !
   keystr = trim(line(1:j-1))
   !
   valstr = trim(line(j+1:))
   !
   ! Remove comments
   !
   jn = index(valstr, '#')
   !
   if (jn > 0) then
      !
      valstr = trim(valstr(1 : jn - 1)) 
      ! 
   endif
   !
   valstr = adjustl(trim(valstr))
   !
   end subroutine 


   subroutine notabs(INSTR,OUTSTR,ILEN)
   ! @(#) convert tabs in input to spaces in output while maintaining columns, assuming a tab is set every 8 characters
   !
   ! USES:
   !       It is often useful to expand tabs in input files to simplify further processing such as tokenizing an input line.
   !       Some FORTRAN compilers hate tabs in input files; some printers; some editors will have problems with tabs
   ! AUTHOR:
   !       John S. Urban
   !
   ! SEE ALSO: 
   !       GNU/Unix commands expand(1) and unexpand(1) 
   !
   use ISO_FORTRAN_ENV, only : ERROR_UNIT     ! get unit for standard error. if not supported yet,  define ERROR_UNIT for your system (typically 0)
   character(len=*),intent(in)   :: INSTR     ! input line to scan for tab characters
   character(len=*),intent(out)  :: OUTSTR    ! tab-expanded version of INSTR produced
   integer,intent(out)           :: ILEN      ! column position of last character put into output string

   integer,parameter             :: TABSIZE=8 ! assume a tab stop is set every 8th column
   character(len=1)              :: c         ! character read from stdin
   integer                       :: ipos      ! position in OUTSTR to put next character of INSTR
   integer                       :: lenin     ! length of input string trimmed of trailing spaces
   integer                       :: lenout    ! number of characters output string can hold
   integer                       :: i10       ! counter that advances thru input string INSTR one character at a time
   !   
   IPOS=1                                  ! where to put next character in output string OUTSTR
   lenin=len(INSTR)                        ! length of character variable INSTR
   lenin=len_trim(INSTR(1:lenin))          ! length of INSTR trimmed of trailing spaces
   lenout=len(OUTSTR)                      ! number of characters output string OUTSTR can hold
   OUTSTR=" "                              ! this SHOULD blank-fill string, a buggy machine required a loop to set all characters
   !
   do i10=1,lenin                          ! look through input string one character at a time
      c=INSTR(i10:i10)
      if(ichar(c) == 9)then                ! test if character is a tab (ADE (ASCII Decimal Equivalent) of tab character is 9)
         IPOS = IPOS + (TABSIZE - (mod(IPOS-1,TABSIZE)))
      else                                 ! c is anything else other than a tab insert it in output string
         if(IPOS > lenout)then
            write(ERROR_UNIT,*)"*notabs* output string overflow"
            exit
         else
            OUTSTR(IPOS:IPOS)=c
            IPOS=IPOS+1
         endif
      endif
   enddo
   !
   ILEN=len_trim(OUTSTR(:IPOS))  ! trim trailing spaces
   return
   !
   end subroutine notabs

   
end module
