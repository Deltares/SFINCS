module snapwave_input

contains

   subroutine read_snapwave_input()
   !
   ! Reads snapwave.inp
   !
   use snapwave_data, only: &
      nmax, mmax, dx, dy, x0, y0, rotation, &
      trefstr, tstartstr, tstopstr, timestep, &
      gamma, alpha, hmin, gridfile, sferic, &
      fw0, fw0_ig, dt, tol, dtheta, &
      upwfile, mskfile, indfile, depfile, obsfile, &
      outputformat, map_filename, his_filename, &
      tstart, tstop, restart, &
      snapwave_jonswapfile, snapwave_bndfile, snapwave_encfile, &
      snapwave_bhsfile, snapwave_btpfile, snapwave_bwdfile
   use snapwave_date
   !
   implicit none
   !
   integer :: dtsec
   integer :: irestart
   character(len=256) :: legacy_bdsfile
   character(len=256) :: legacy_bzsfile
   !
   write(*,*) 'Reading input file ...'
   !
   open(500, file='snapwave.inp')
   !
   ! Input section
   !
   call read_int_input(500,'nmax',nmax,0)
   call read_int_input(500,'mmax',mmax,0)
   call read_real_input(500,'dx',dx,0.0)
   call read_real_input(500,'dy',dy,0.0)
   call read_real_input(500,'x0',x0,0.0)
   call read_real_input(500,'y0',y0,0.0)
   call read_real_input(500,'rotation',rotation,0.0)
   call read_char_input(500,'tref',trefstr,'20000101 000000')
   call read_char_input(500,'tstart',tstartstr,'20000101 000000')
   call read_char_input(500,'tstop',tstopstr,'20000101 000000')
   call read_real_input(500,'timestep',timestep,3600.0)
   call read_real_input(500,'gamma',gamma,0.7)
   call read_real_input(500,'alpha',alpha,1.0)
   call read_real_input(500,'hmin',hmin,0.1)
   call read_char_input(500,'gridfile',gridfile,'.txt')
   call read_int_input(500,'sferic',sferic,0)
   call read_real_input(500,'fw',fw0,0.01)
   call read_real_input(500,'fwig',fw0_ig,0.015)
   call read_real_input(500,'dt',dt,36000.0)
   call read_real_input(500,'tol',tol,10.0)
   call read_real_input(500,'dtheta',dtheta,10.0)
!   call read_int_input(500,'ntheta',ntheta,36)
!   call read_int_input(500,'nHrel',nHrel,1)
!   call read_char_input(500,'hhtabname',hhtabname,'')
!   call read_char_input(500,'Htabname',Htabname,'')
!   call read_char_input(500,'Dwtabname',Dwtabname,'')
!   call read_char_input(500,'Ftabname',Ftabname,'')
!   call read_char_input(500,'Cgtabname',Cgtabname,'')
!   call read_char_input(500,'cthetafactabname',cthetafactabname,'')
!   call read_char_input(500,'waterlevelfile',waterlevelfile,'')

   ! New Deltares-prefixed SnapWave filenames
   call read_char_input(500,'jonswapfile',snapwave_jonswapfile,'')
   call read_char_input(500,'bndfile',snapwave_bndfile,'')
   call read_char_input(500,'encfile',snapwave_encfile,'')
   call read_char_input(500,'bhsfile',snapwave_bhsfile,'')
   call read_char_input(500,'btpfile',snapwave_btpfile,'')
   call read_char_input(500,'bwdfile',snapwave_bwdfile,'')

   ! Legacy keys: keep parsing them so older snapwave.inp files do not fail,
   ! even if the merged Deltares branch no longer stores them in snapwave_data.
   call read_char_input(500,'bdsfile',legacy_bdsfile,'')
   call read_char_input(500,'bzsfile',legacy_bzsfile,'')

   call read_char_input(500,'upwfile',upwfile,'')
   call read_char_input(500,'mskfile',mskfile,'')
   call read_char_input(500,'indfile',indfile,'')
   call read_char_input(500,'depfile',depfile,'')
   call read_char_input(500,'obsfile',obsfile,'')
   call read_char_input(500,'outputformat',outputformat,'bin')
   call read_char_input(500,'map_file',map_filename,'')
   call read_char_input(500,'his_file',his_filename,'')
   call read_int_input(500,'restart',irestart,1)
   !
   close(500)
   !
   call time_difference(trefstr, tstartstr, dtsec)  ! time difference in seconds between tstart and tref
   tstart = dtsec*1.0
   call time_difference(trefstr, tstopstr, dtsec)
   tstop = dtsec*1.0
   !
   mmax = mmax + 2  ! Original mmax and nmax are for number of cells in bathy grid. Add two dummy rows.
   nmax = nmax + 2
   !
   restart = .true.
   if (irestart == 0) restart = .false.
   !
   end subroutine read_snapwave_input



   subroutine read_real_input(fileid,keyword,value,default)
   !
   character(*), intent(in) :: keyword
   character(len=256)       :: keystr
   character(len=256)       :: valstr
   character(len=256)       :: line
   integer, intent(in)      :: fileid
   real*4, intent(out)      :: value
   real*4, intent(in)       :: default
   integer :: j, stat
   !
   value = default
   rewind(fileid)
   do while(.true.)
      read(fileid,'(a)',iostat = stat) line
      if (stat < 0) exit
      j = index(line,'=')
      if (j <= 0) cycle
      keystr = trim(line(1:j-1))
      if (trim(keystr) == trim(keyword)) then
         valstr = trim(line(j+1:256))
         read(valstr,*) value
         exit
      endif
   enddo
   !
   end subroutine read_real_input


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
   integer :: j, stat, m
   !
   allocate(value(nr))
   !
   value = default
   rewind(fileid)
   do while(.true.)
      read(fileid,'(a)',iostat = stat) line
      if (stat < 0) exit
      j = index(line,'=')
      if (j <= 0) cycle
      keystr = trim(line(1:j-1))
      if (trim(keystr) == trim(keyword)) then
         valstr = trim(line(j+1:256))
         read(valstr,*) (value(m), m = 1, nr)
         exit
      endif
   enddo
   !
   end subroutine read_real_array_input


   subroutine read_int_input(fileid,keyword,value,default)
   !
   character(*), intent(in) :: keyword
   character(len=256)       :: keystr
   character(len=256)       :: valstr
   character(len=256)       :: line
   integer, intent(in)      :: fileid
   integer, intent(out)     :: value
   integer, intent(in)      :: default
   integer :: j, stat
   !
   value = default
   rewind(fileid)
   do while(.true.)
      read(fileid,'(a)',iostat = stat) line
      if (stat < 0) exit
      j = index(line,'=')
      if (j <= 0) cycle
      keystr = trim(line(1:j-1))
      if (trim(keystr) == trim(keyword)) then
         valstr = trim(line(j+1:256))
         read(valstr,*) value
         exit
      endif
   enddo
   !
   end subroutine read_int_input


   subroutine read_char_input(fileid,keyword,value,default)
   !
   character(*), intent(in)  :: keyword
   character(len=256)        :: keystr
   character(len=256)        :: valstr
   character(len=256)        :: line
   integer, intent(in)       :: fileid
   character(*), intent(in)  :: default
   character(*), intent(out) :: value
   integer :: j, stat
   !
   value = default
   rewind(fileid)
   do while(.true.)
      read(fileid,'(a)',iostat = stat) line
      if (stat < 0) exit
      j = index(line,'=')
      if (j <= 0) cycle
      keystr = trim(line(1:j-1))
      if (trim(keystr) == trim(keyword)) then
         valstr = adjustl(trim(line(j+1:256)))
         value = valstr
         exit
      endif
   enddo
   !
   end subroutine read_char_input

end module snapwave_input
