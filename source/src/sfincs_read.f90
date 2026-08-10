module sfincs_read
   !
   ! Keyword readers for SFINCS input files. The legacy read_*_input
   ! helpers are kept for callers outside sfincs_input (snapwave,
   ! spiderweb); new code should use the generic get_keyword(...)
   ! interface, which supports deprecated-alias lists.
   !
   use sfincs_log, only: write_log
   !
   interface get_keyword
      module procedure get_keyword_real
      module procedure get_keyword_int
      module procedure get_keyword_char
      module procedure get_keyword_logical
   end interface
   !

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


   
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine get_keyword_real(fileid, keyword, value, default, legacy)
      !
      ! Read one real*4 keyword. Tries `keyword` first; if absent, walks
      ! the optional `legacy` list of deprecated aliases and emits a
      ! one-line deprecation warning per matched alias. Falls back to
      ! `default` when nothing matches.
      !
      ! Called from: read_sfincs_input (sfincs_input).
      !
      implicit none
      !
      integer,                           intent(in)            :: fileid
      character(*),                      intent(in)            :: keyword
      real*4,                            intent(out)           :: value
      real*4,                            intent(in)            :: default
      character(len=*), dimension(:),    intent(in), optional  :: legacy
      !
      character(len=256) :: valstr
      logical            :: found
      integer            :: i
      !
      call find_value(fileid, keyword, valstr, found)
      if (found) then
         read(valstr, *) value
         return
      endif
      !
      if (present(legacy)) then
         !
         do i = 1, size(legacy)
            !
            call find_value(fileid, trim(legacy(i)), valstr, found)
            !
            if (found) then
               read(valstr, *) value
               call warn_legacy(trim(legacy(i)), keyword)
               return
            endif
            !
         enddo
         !
      endif
      !
      value = default
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine get_keyword_int(fileid, keyword, value, default, legacy)
      !
      ! Read one integer keyword. See get_keyword_real for the semantics.
      !
      ! Called from: read_sfincs_input (sfincs_input).
      !
      implicit none
      !
      integer,                           intent(in)            :: fileid
      character(*),                      intent(in)            :: keyword
      integer,                           intent(out)           :: value
      integer,                           intent(in)            :: default
      character(len=*), dimension(:),    intent(in), optional  :: legacy
      !
      character(len=256) :: valstr
      logical            :: found
      integer            :: i
      !
      call find_value(fileid, keyword, valstr, found)
      if (found) then
         read(valstr, *) value
         return
      endif
      !
      if (present(legacy)) then
         !
         do i = 1, size(legacy)
            !
            call find_value(fileid, trim(legacy(i)), valstr, found)
            !
            if (found) then
               read(valstr, *) value
               call warn_legacy(trim(legacy(i)), keyword)
               return
            endif
            !
         enddo
         !
      endif
      !
      value = default
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine get_keyword_char(fileid, keyword, value, default, legacy)
      !
      ! Read one character-string keyword. See get_keyword_real for the
      ! semantics. The entire right-hand side (after trailing comments
      ! are stripped) becomes `value`.
      !
      ! Called from: read_sfincs_input (sfincs_input).
      !
      implicit none
      !
      integer,                           intent(in)            :: fileid
      character(*),                      intent(in)            :: keyword
      character(*),                      intent(out)           :: value
      character(*),                      intent(in)            :: default
      character(len=*), dimension(:),    intent(in), optional  :: legacy
      !
      character(len=256) :: valstr
      logical            :: found
      integer            :: i
      !
      call find_value(fileid, keyword, valstr, found)
      if (found) then
         value = valstr
         return
      endif
      !
      if (present(legacy)) then
         !
         do i = 1, size(legacy)
            !
            call find_value(fileid, trim(legacy(i)), valstr, found)
            !
            if (found) then
               value = valstr
               call warn_legacy(trim(legacy(i)), keyword)
               return
            endif
            !
         enddo
         !
      endif
      !
      value = default
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine get_keyword_logical(fileid, keyword, value, default, legacy)
      !
      ! Read one logical keyword. Accepts `1`, `y`, `Y`, `t`, `T` as true;
      ! anything else (including absence → `default`, and `0`, `n`, `N`,
      ! `f`, `F`) as false.
      !
      ! Called from: read_sfincs_input (sfincs_input).
      !
      implicit none
      !
      integer,                           intent(in)            :: fileid
      character(*),                      intent(in)            :: keyword
      logical,                           intent(out)           :: value
      logical,                           intent(in)            :: default
      character(len=*), dimension(:),    intent(in), optional  :: legacy
      !
      character(len=256) :: valstr
      logical            :: found
      integer            :: i
      !
      call find_value(fileid, keyword, valstr, found)
      if (found) then
         value = parse_logical(valstr)
         return
      endif
      !
      if (present(legacy)) then
         !
         do i = 1, size(legacy)
            !
            call find_value(fileid, trim(legacy(i)), valstr, found)
            !
            if (found) then
               value = parse_logical(valstr)
               call warn_legacy(trim(legacy(i)), keyword)
               return
            endif
            !
         enddo
         !
      endif
      !
      value = default
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine get_keyword_real_array(fileid, keyword, value, default, nr, legacy)
      !
      ! Read one whitespace-separated real*4 array keyword. Allocates
      ! `value(nr)` on the way in and fills it from the matching line.
      ! Same fallback semantics as get_keyword_real.
      !
      ! Called from: read_sfincs_input (sfincs_input).
      !
      implicit none
      !
      integer,                              intent(in)            :: fileid
      character(*),                         intent(in)            :: keyword
      integer,                              intent(in)            :: nr
      real*4,                               intent(in)            :: default
      real*4, dimension(:),                 intent(out), allocatable :: value
      character(len=*), dimension(:),       intent(in),  optional :: legacy
      !
      character(len=256) :: valstr
      logical            :: found
      integer            :: i, m
      !
      allocate(value(nr))
      !
      call find_value(fileid, keyword, valstr, found)
      if (found) then
         read(valstr, *) (value(m), m = 1, nr)
         return
      endif
      !
      if (present(legacy)) then
         !
         do i = 1, size(legacy)
            !
            call find_value(fileid, trim(legacy(i)), valstr, found)
            !
            if (found) then
               read(valstr, *) (value(m), m = 1, nr)
               call warn_legacy(trim(legacy(i)), keyword)
               return
            endif
            !
         enddo
         !
      endif
      !
      value = default
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine find_value(fileid, keyword, valstr, found)
      !
      ! Scan an already-open sfincs.inp once for the given `keyword`.
      ! Returns the raw right-hand-side value string and whether the key
      ! was matched.
      !
      ! Called from: get_keyword_real / get_keyword_int / get_keyword_char /
      ! get_keyword_logical / get_keyword_real_array.
      !
      implicit none
      !
      integer,      intent(in)  :: fileid
      character(*), intent(in)  :: keyword
      character(*), intent(out) :: valstr
      logical,      intent(out) :: found
      !
      character(len=256) :: keystr
      character(len=256) :: line
      integer            :: stat
      !
      found  = .false.
      valstr = ''
      !
      rewind(fileid)
      !
      do while (.true.)
         !
         read(fileid, '(a)', iostat=stat) line
         if (stat == -1) exit
         !
         call read_line(line, keystr, valstr)
         !
         if (trim(keystr) == trim(keyword)) then
            found = .true.
            return
         endif
         !
      enddo
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine warn_legacy(legacy_key, new_key)
      !
      ! Emit a one-line deprecation warning to the log. Called whenever
      ! a legacy keyword alias was matched; the user can silence this
      ! by migrating the keyword in their sfincs.inp.
      !
      ! Called from: get_keyword_real / get_keyword_int / get_keyword_char /
      ! get_keyword_logical / get_keyword_real_array.
      !
      implicit none
      !
      character(*), intent(in) :: legacy_key
      character(*), intent(in) :: new_key
      !
      character(len=512) :: msg
      !
      write(msg, '(a,a,a,a,a)') ' Warning : sfincs.inp keyword "', trim(legacy_key), &
         '" is deprecated, use "', trim(new_key), '" instead'
      call write_log(trim(msg), 1)
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   function parse_logical(valstr) result(value)
      !
      ! Map an sfincs.inp value string to a logical. `1`, `y`, `Y`, `t`,
      ! `T` at position 1 are true; everything else is false.
      !
      ! Called from: get_keyword_logical (this module).
      !
      implicit none
      !
      character(*), intent(in) :: valstr
      logical                  :: value
      !
      value = (valstr(1:1) == '1' .or. valstr(1:1) == 'y' .or. valstr(1:1) == 'Y' .or. &
               valstr(1:1) == 't' .or. valstr(1:1) == 'T')
      !
   end function
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine read_line(line0, keystr, valstr)
      !
      ! Split one `key = value` line into key and value substrings.
      ! Strips leading/trailing whitespace, any tab characters (replaced
      ! by spaces via notabs), and a trailing `#`-delimited inline
      ! comment. Blank lines and lines starting with `#`, `!`, or `@`
      ! return empty strings.
      !
      ! Called from: find_value, and the legacy read_*_input helpers above.
      !
      implicit none
      !
      character(*), intent(in)  :: line0
      character(*), intent(out) :: keystr
      character(*), intent(out) :: valstr
      !
      character(len=256) :: line
      integer            :: j, ilen, jn
      !
      keystr = ''
      valstr = ''
      !
      ! Expand tabs to spaces in-place.
      !
      call notabs(line0, line, ilen)
      !
      ! Remove Windows-style CR (carriage return) line ending if present.
      !
      jn = index(line, achar(13))
      if (jn > 0) line = line(1:jn - 1)
      line = trim(line)
      !
      if (line(1:1) == '#' .or. line(1:1) == '!' .or. line(1:1) == '@') return
      !
      j = index(line, '=')
      if (j == 0) return
      !
      keystr = trim(line(1:j - 1))
      valstr = trim(line(j + 1:))
      !
      ! Strip inline comment after `#`.
      !
      jn = index(valstr, '#')
      if (jn > 0) valstr = trim(valstr(1:jn - 1))
      !
      valstr = adjustl(trim(valstr))
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine notabs(instr, outstr, ilen)
      !
      ! Expand embedded tab characters into spaces while keeping columns
      ! aligned (tab stops every 8 characters). Lets downstream tokenizers
      ! treat `key<tab>=<tab>value` and `key = value` identically.
      !
      ! Author: John S. Urban. See also GNU/Unix commands expand(1) /
      ! unexpand(1).
      !
      ! Called from: read_line (this module).
      !
      use iso_fortran_env, only : error_unit
      !
      implicit none
      !
      character(len=*), intent(in)  :: instr     ! input line (may contain tab characters)
      character(len=*), intent(out) :: outstr    ! tab-expanded output
      integer,          intent(out) :: ilen      ! column position of last character written
      !
      integer, parameter :: tabsize = 8          ! tab stops every 8th column
      character(len=1)   :: c
      integer            :: ipos                 ! position in outstr for next character
      integer            :: lenin                ! length of instr (trailing blanks trimmed)
      integer            :: lenout               ! capacity of outstr
      integer            :: i10                  ! cursor through instr
      !
      ipos   = 1
      lenin  = len(instr)
      lenin  = len_trim(instr(1:lenin))
      lenout = len(outstr)
      outstr = ' '
      !
      do i10 = 1, lenin
         !
         c = instr(i10:i10)
         !
         if (ichar(c) == 9) then
            !
            ! Tab character: advance ipos to the next tab stop.
            !
            ipos = ipos + (tabsize - (mod(ipos - 1, tabsize)))
            !
         else
            !
            if (ipos > lenout) then
               write(error_unit, *) '*notabs* output string overflow'
               exit
            else
               outstr(ipos:ipos) = c
               ipos              = ipos + 1
            endif
            !
         endif
         !
      enddo
      !
      ilen = len_trim(outstr(:ipos))
      !
   end subroutine
   !
end module
