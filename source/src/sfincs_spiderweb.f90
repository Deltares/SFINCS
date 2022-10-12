module sfincs_spiderweb

contains

   subroutine read_spw_file(filename,nt,nrows,ncols,spwrad,time,xe,ye,vmag,vdir,pdrp,prcp,nquant,trefstr)
   !
   implicit none
   !
   integer j, j2, it, ip, n, m, nquant, stat, nheader, id
   integer dtsec
   integer iyspw,imspw,idspw
   character*4 cyspw 
   character*2 cmspw 
   character*2 cdspw 
   !
   character*256, intent(in) :: filename
   character*256 :: cdummy
   character*256 :: line
   character*256 :: keystr
   character*256 :: valstr
   character*15  :: datespw
   character*15  :: trefstr
   character*256 :: timstr
   !
   integer, intent(in) :: nt
   integer, intent(in) :: nrows
   integer, intent(in) :: ncols
   !
   real*4, intent(in)                 :: spwrad
   real*4, dimension(nt), intent(out) :: time 
   real*4, dimension(nt), intent(out) :: xe 
   real*4, dimension(nt), intent(out) :: ye 
   real*4, dimension(nt,nrows,ncols), intent(out) :: vmag 
   real*4, dimension(nt,nrows,ncols), intent(out) :: vdir 
   real*4, dimension(nt,nrows,ncols), intent(out) :: pdrp 
   real*4, dimension(nt,nrows,ncols), intent(out) :: prcp 
   !
   vmag = 0.0
   vdir = 0.0
   pdrp = 0.0
   prcp = 0.0
   !
   ! Open SPW file
   !
   open(888,file=filename)   
   !
   ! Read header
   nheader = 0
   do ip = 1, 30     !read first 30 lines to find first TIME block
      read(888,'(a)')line
      !
      id=index(line,'TIME')
      !
      if (id == 1) then
         nheader = ip - 1
         exit
      endif
   enddo      
   !
   if (nheader ==0) then !in case not possible to find header
       !
       if (nquant==4) then
           nheader = 18
       else
           nheader = 16
       endif
       !
       write(*,*)'Debug: not possible to automatically find number of headers in .spw file, switch to hardcoded value of nheader = ',nheader
       !       
   endif
   !      
   rewind(888)   ! needed otherwise doesn't start reading in begin of file     
   !
   do ip = 1, nheader
      read(888,*)cdummy
   enddo
   !
   do it = 1, nt
      !
      ! Read time
      !
      read(888,'(a)')line
      !
      call compute_time_in_seconds(line,trefstr,dtsec)
!      j=index(line,'=')      
!      j2=index(line,'m')      
!      keystr = trim(line(1:j-1))
!      valstr = trim(line(j+1:j2-1))
!      timstr = trim(line(j2+14:j2+24))
!      read(valstr,*)time(it)
!      !
!      ! Compute difference in seconds between spiderweb reference time and simulation start time
!      !
!      read(timstr,'(A,1X,A,1X,A)')cyspw,cmspw,cdspw
!      datespw = cyspw // cmspw // cdspw // ' 000000'
!      call time_difference(datespw,trefstr,dtsec)
!      !
!      time(it) = time(it)*60 - dtsec*1.0 ! Convert to seconds w.r.t. reference time
      time(it) = dtsec*1.0 ! Convert to seconds w.r.t. reference time
      !
      ! Xe
      !
      read(888,'(a)')line
      j=index(line,'=')      
      keystr = trim(line(1:j-1))
      valstr = trim(line(j+1:j+12))
      read(valstr,*)xe(it)
      !
      ! Ye
      !
      read(888,'(a)')line
      j=index(line,'=')      
      keystr = trim(line(1:j-1))
      valstr = trim(line(j+1:j+12))
      read(valstr,*)ye(it)
      !
      ! Peye (dummy)
      !
      read(888,'(a)')line
      !
      do n = 1, nrows 
         read(888,*)(vmag(it,n,m), m = 1, ncols)
      enddo             
      do n = 1, nrows 
         read(888,*)(vdir(it,n,m), m = 1, ncols)
      enddo             
      do n = 1, nrows 
         read(888,*)(pdrp(it,n,m), m = 1, ncols)
      enddo             
      !
      if (nquant==4) then 
         do n = 1, nrows 
            read(888,*)(prcp(it,n,m), m = 1, ncols)
         enddo             
      endif
      !
   enddo
   !
   close(888)
   !
   end subroutine

   subroutine read_amuv_file(filename,nt,nrows,ncols,time,uv,trefstr)
   !
   implicit none
   !
   integer j, j2, it, ip, n, m, nquant, stat, nheader, id
   integer dtsec
   integer iyspw,imspw,idspw
   character*4 cyspw 
   character*2 cmspw 
   character*2 cdspw 
   !
   character*256, intent(in) :: filename
   character*256 :: cdummy
   character*256 :: line
   character*256 :: keystr
   character*256 :: valstr
   character*15  :: datespw
   character*15  :: trefstr
   character*256 :: timstr
   !
   integer, intent(in) :: nt
   integer, intent(in) :: nrows
   integer, intent(in) :: ncols
   !
   real*4, dimension(nt), intent(out) :: time 
   real*4, dimension(nt,nrows,ncols), intent(out) :: uv
   !
   uv = 0.0
   !
   ! Open AMU/AMV file
   !
   open(888,file=filename)   
   !
   ! Read header
   nheader = 0
   do ip = 1, 30     !read first 30 lines to find first TIME block
      read(888,'(a)')line
      !
      id=index(line,'TIME')
      !
      if (id == 1) then
         nheader = ip - 1 !-2 as in Hurrywave?
         exit
      endif
   enddo      
   !
   if (nheader ==0) then !in case not possible to find header
       !
       nheader = 13
       !
       write(*,*)'Debug: not possible to automatically find number of headers in amu/amv/amp file, switch to hardcoded value of nheader = ',nheader
       !       
   endif
   !      
   rewind(888)   ! needed otherwise doesn't start reading in begin of file     
   !
   do ip = 1, nheader
      read(888,*)cdummy
   enddo
   !
   do it = 1, nt
      !
      ! Read time
      !
      read(888,'(a)')line
      !
      call compute_time_in_seconds(line,trefstr,dtsec)
      !
!      write(*,*)line
!      j=index(line,'=')      
!      j2=index(line,'m')      
!      keystr = trim(line(1:j-1))
!      valstr = trim(line(j+1:j2-1))
!      timstr = trim(line(j2+14:j2+24))
!      read(valstr,*)time(it)
!      !
!      ! Compute difference in seconds between spiderweb reference time and simulation start time
!      !
!      read(timstr,'(A,1X,A,1X,A)')cyspw,cmspw,cdspw
!      datespw = cyspw // cmspw // cdspw // ' 000000'
!      call time_difference(datespw,trefstr,dtsec)
!      !
!      time(it) = time(it)*60 - dtsec*1.0 ! Convert to seconds w.r.t. reference time
      time(it) = dtsec*1.0
      !
      do n = 1, nrows 
         read(888,*)(uv(it,n,m), m = 1, ncols)
      enddo             
      !
   enddo
   !
   close(888)
   !
   end subroutine
   
   subroutine read_spw_dimensions(filename,nt,nrows,ncols,spwrad,nquant)
      !
      implicit none
      !
      integer j, it, ip, n, m, nquant, stat, nheader, id    
      !
      character*256, intent(in) :: filename
      character*256 :: cdummy
      character*256 :: line
      character*256 :: keystr
      character*256 :: valstr
      !
      integer, intent(out) :: nt
      integer, intent(out) :: nrows
      integer, intent(out) :: ncols
      real*4,  intent(out) :: spwrad
      !
      ! Read input file
      !
      open(888,file=filename)   
      !
      call read_int_input(888,'n_cols',ncols,0)
      call read_int_input(888,'n_rows',nrows,0)
      call read_int_input(888,'n_quantity',nquant,0)
      call read_real_input(888,'spw_radius',spwrad,0.0)
      !   
      close(888)
      !
      ! Open SPW file
      !
      open(888,file=filename)   
      !
      ! Read header
      nheader = 0
      do ip = 1, 30     !read first 30 lines to find first TIME block
         read(888,'(a)')line
         !
         id=index(line,'TIME')
         !
         if (id == 1) then
            nheader = ip - 1 !-2 as in Hurrywave?
            exit
         endif
      enddo      
      !
      if (nheader ==0) then !in case not possible to find header
          !
          if (nquant==4) then
              nheader = 18
          else
              nheader = 16
          endif
          !
          write(*,*)'Debug: not possible to automatically find number of headers in .spw file, switch to hardcoded value of nheader = ',nheader
          !       
      endif
      !      
      rewind(888)   ! needed otherwise doesn't start reading in begin of file     
      !
      do ip = 1, nheader
         read(888,*)cdummy
      enddo
      !           
      ! Read times
      !
      nt = 0
      !
      do while(.true.)  
         read(888,'(a)',iostat = stat)line         
         if (stat<0) exit
         nt = nt + 1      
         j=index(line,'=')      
         keystr = trim(line(1:j-1))
         valstr = trim(line(j+1:j+12))
         read(888,*)cdummy
         read(888,*)cdummy
         read(888,*)cdummy
         do ip = 1, nrows*nquant
            read(888,*)cdummy
         enddo      
      enddo
      !
      close(888)   
     !      
   end subroutine
   !
   !
   !
   subroutine read_amuv_dimensions(filename,nt,nrows,ncols,x_llcorner,y_llcorner,dx,dy,nquant)
      !
      implicit none
      !
      integer j, it, ip, n, m, nquant, stat, nheader, id    
      !
      character*256, intent(in) :: filename
      character*256 :: cdummy
      character*256 :: line
      character*256 :: keystr
      character*256 :: valstr
      !
      integer, intent(out) :: nt
      integer, intent(out) :: nrows
      integer, intent(out) :: ncols
      real*4,  intent(out) :: x_llcorner
      real*4,  intent(out) :: y_llcorner
      real*4,  intent(out) :: dx
      real*4,  intent(out) :: dy
      !
      ! Read input file
      !
      open(888, file=filename)   
      !
      call read_int_input(888,'n_cols',ncols,0)
      call read_int_input(888,'n_rows',nrows,0)
      call read_int_input(888,'n_quantity',nquant,0)
      call read_real_input(888,'x_llcorner',x_llcorner,0.0)
      call read_real_input(888,'y_llcorner',y_llcorner,0.0)
      call read_real_input(888,'dx',dx,0.0)
      call read_real_input(888,'dy',dy,0.0)
      !   
      close(888)
      !
      ! Open AMU file
      !
      open(888, file=filename)   
      !
      ! Read header
      nheader = 0
      do ip = 1, 30     !read first 30 lines to find first TIME block
         !
         read(888,'(a)')line
         !
         id=index(line,'TIME')
         !
         if (id == 1) then
            nheader = ip - 1
            exit
         endif
         !
      enddo
      !
      if (nheader ==0) then !in case not possible to find header
          !
          nheader = 13
          !
          write(*,*)'Debug: not possible to automatically find number of headers in amu/amv/amp file, switch to hardcoded value of nheader = ',nheader
          !       
      endif
      !      
      rewind(888)   ! needed otherwise doesn't start reading in begin of file  
      !
      do ip = 1, nheader
         read(888,*)cdummy
      enddo
      !
      ! Read times
      !
      nt = 0
      !
      do while(.true.)  
         read(888,'(a)',iostat = stat)line
         if (stat<0) exit
         nt = nt + 1      
         j=index(line,'=')      
         keystr = trim(line(1:j-1))
         valstr = trim(line(j+1:j+12))
         do ip = 1, nrows*nquant
            read(888,*)cdummy
         enddo      
      enddo
      !
      close(888)   
     !
   end subroutine

   
   subroutine read_real_input(fileid,keyword,value,default)
   !
   character(*), intent(in) :: keyword
   character(len=256)       :: keystr
   character(len=256)       :: valstr
   character(len=256)       :: line
   integer, intent(in)      :: fileid
   real*4, intent(out)      :: value
   real*4, intent(in)       :: default
   integer j,stat
   !
   value = default
   rewind(fileid)   
   do while(.true.)
      read(fileid,'(a)',iostat = stat)line
      if (stat<0) exit
      j=index(line,'=')      
      keystr = trim(line(1:j-1))
      if (trim(keystr)==trim(keyword)) then
         valstr = trim(line(j+1:256))
         read(valstr,*)value
         exit
      endif
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
   integer j,stat
   !
   value = default
   rewind(fileid)   
   do while(.true.)
      read(fileid,'(a)',iostat = stat)line
      if (stat<0) exit
      j=index(line,'=')      
      keystr = trim(line(1:j-1))
      if (trim(keystr)==trim(keyword)) then
         valstr = trim(line(j+1:256))
         read(valstr,*)value         
         exit
      endif
   enddo 
   !
   end subroutine

   
   subroutine compute_time_in_seconds(line,trefstr,dtsec)
   !
   use sfincs_date
   !   
   character(*), intent(in) :: line
   character*15, intent(in) :: trefstr
   integer, intent(out)     :: dtsec
   !
   integer j, j2, iopt
   character*4 cyspw 
   character*2 cmspw 
   character*2 cdspw 
   !
   character*256 :: cdummy
   character*256 :: keystr
   character*256 :: valstr
   character*15  :: datespw
   character*256 :: timstr
   real          :: tim   
   !
   j=index(line,'=')      
   j2=index(line,'minutes')    
   iopt = 1
   if (j2==0) then
      iopt = 2
      j2=index(line,'hours')      
   endif    
   keystr = trim(line(1:j-1))
   valstr = trim(line(j+1:j2-1))
   if (iopt==1) then
      !
      ! minutes
      !
      timstr = trim(line(j2+14:j2+24))
      !
   else
      !
      ! hours
      !
      timstr = trim(line(j2+12:j2+22))
      !
   endif
   !
   read(valstr,*)tim
   !
   ! Compute difference in seconds between spiderweb reference time and simulation start time
   !
   read(timstr,'(A,1X,A,1X,A)')cyspw,cmspw,cdspw
   datespw = cyspw // cmspw // cdspw // ' 000000'
   call time_difference(datespw,trefstr,dtsec)
   !
   if (iopt==1) then
      dtsec = tim*60 - dtsec ! Convert to seconds w.r.t. reference time
   else
      dtsec = tim*3600 - dtsec ! Convert to seconds w.r.t. reference time
   endif
   !
   end subroutine
   
end module        