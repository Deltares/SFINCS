module sfincs_input

contains

   subroutine read_sfincs_input()
   !
   ! Reads sfincs.inp
   !
   use sfincs_data
   use sfincs_date
   !
   implicit none
   !
   integer dtsec
   !
   character*256 wmsigstr 
   !   
   write(*,*)'Reading input file ...'
   !
   open(500, file='sfincs.inp')   
   !
   call read_int_input(500,'mmax',mmax,0)
   call read_int_input(500,'nmax',nmax,0)
   call read_real_input(500,'dx',dx,0.0)
   call read_real_input(500,'dy',dy,0.0)
   call read_real_input(500,'x0',x0,0.0)
   call read_real_input(500,'y0',y0,0.0)
   call read_real_input(500,'rotation',rotation,0.0)
   call read_char_input(500,'tref',trefstr,'20000101 000000')
   call read_char_input(500,'tstart',tstartstr,'20000101 000000')
   call read_char_input(500,'tstop',tstopstr,'20000101 000000')
   call read_real_input(500,'tspinup',tspinup,0.0)
   call read_real_input(500,'t0out',t0out,-999.0)
   call read_real_input(500,'t1out',t1out,-999.0)
   call read_real_input(500,'dtout',dtmapout,0.0)
   call read_real_input(500,'dtmaxout',dtmaxout,999999.0)
   call read_real_input(500,'dtrstout',dtrstout,0.0)
   call read_real_input(500,'trstout',trst,-999.0)
   call read_real_input(500,'dthisout',dthisout,600.0)
   call read_real_input(500,'dtwave',dtwave,3600.0)
   call read_real_input(500,'dtwnd',dtwindupd,1800.0)
   call read_real_input(500,'alpha',alfa,0.50)
   call read_real_input(500,'theta',theta,1.0)
   call read_real_input(500,'manning',manning,0.04)
   call read_real_input(500,'manning_land',manning_land,-999.0)
   call read_real_input(500,'manning_sea',manning_sea,-999.0)
   call read_real_input(500,'rgh_lev_land',rghlevland,0.0)
   call read_real_input(500,'zsini',zini,0.0)
   call read_real_input(500,'qinf',qinf,0.0)
   call read_real_input(500,'igperiod',tig,120.0)
   call read_real_input(500,'dtmax',dtmax,60.0)
   call read_real_input(500,'huthresh',huthresh,0.05)
   call read_real_input(500,'rhoa',rhoa,1.25)
   call read_real_input(500,'rhow',rhow,1024.0)
   call read_real_input(500,'maxlev',max_elev,99999.0)
   call read_char_input(500,'inputformat',inputtype,'asc')
   call read_char_input(500,'outputformat',outputtype,'asc')
   call read_char_input(500,'outputtype_map',outputtype_map,'nil')
   call read_char_input(500,'outputtype_his',outputtype_his,'nil')
   call read_int_input(500,'bndtype',bndtype,1)
   call read_int_input(500,'advection',iadvection,0)
   call read_int_input(500,'nfreqsig',nfreqsig,100)
   call read_real_input(500,'freqminig',freqminig,0.0)
   call read_real_input(500,'freqmaxig',freqmaxig,0.1)
   call read_real_input(500,'latitude',latitude,0.0)
   call read_real_input(500,'pavbnd',pavbnd,0.0)
   call read_real_input(500,'gapres',gapres,101200.0)
   call read_int_input(500,'baro',baro,1)
   call read_char_input(500,'utmzone',utmzone,'nil')
   call read_int_input(500,'epsg',epsg,0)
   call read_char_input(500,'epsg',epsg_code,'nil')      
   call read_real_input(500,'stopdepth',stopdepth,100.0)
   call read_real_input(500,'advlim',advlim,9999.9)
   call read_real_input(500,'qinf_zmin',qinf_zmin,0.0)
   call read_real_input(500,'horton_decay',horton_decay,0.00005)
   call read_real_input(500,'horton_time',horton_time,0.0)
   call read_real_input(500,'btfilter',btfilter,60.0)
   call read_real_input(500,'sfacinf',sfacinf,0.2)
   call read_int_input(500,'radstr',iradstr,0)
   call read_int_input(500,'crsgeo',igeo,0)
   call read_int_input(500,'coriolis',icoriolis,1)
   call read_int_input(500,'amprblock',iamprblock,1)
   call read_real_input(500,'spwmergefrac',spw_merge_frac,0.5)
   call read_int_input(500,'global',iglobal,0)
   call read_real_input(500,'nuvisc',nuviscinp,-999.0)
   call read_real_input(500,'nuviscdim',nuviscdim,1.0)   
   call read_int_input(500,'viscosity',iviscosity,0)
   call read_int_input(500,'spinup_meteo', ispinupmeteo, 0)
   call read_real_input(500,'waveage',waveage,-999.0)
   call read_int_input(500,'snapwave', isnapwave, 0)
   call read_int_input(500,'dtoutfixed', ioutfixed, 1)
   call read_real_input(500,'wmtfilter',wmtfilter,600.0)
   call read_real_input(500,'wmfred',wavemaker_freduv,0.99)
   call read_char_input(500,'wmsignal',wmsigstr,'spectrum')   
   !
   ! Domain
   !
   call read_char_input(500,'qtrfile',qtrfile,'none')
   call read_char_input(500,'depfile',depfile,'none')
   call read_char_input(500,'inifile',zsinifile,'none')
   call read_char_input(500,'rstfile',rstfile,'none')
   call read_char_input(500,'mskfile',mskfile,'none')
   call read_char_input(500,'indexfile',indexfile,'none')
   call read_char_input(500,'cstfile',cstfile,'none')
   call read_char_input(500,'sbgfile',sbgfile,'none')
   call read_char_input(500,'thdfile',thdfile,'none')
   call read_char_input(500,'weirfile',weirfile,'none')
   call read_char_input(500,'manningfile',manningfile,'none')   
   call read_char_input(500,'drnfile',drnfile,'none')
   !
   ! Forcing
   !
   call read_char_input(500,'bndfile',bndfile,'none')
   call read_char_input(500,'bzsfile',bzsfile,'none')
   call read_char_input(500,'bzifile',bzifile,'none')
   call read_char_input(500,'bwvfile',bwvfile,'none')
   call read_char_input(500,'bhsfile',bhsfile,'none')
   call read_char_input(500,'btpfile',btpfile,'none')
   call read_char_input(500,'bwdfile',bwdfile,'none')
   call read_char_input(500,'bdsfile',bdsfile,'none')
   call read_char_input(500,'wfpfile',wfpfile,'none')
   call read_char_input(500,'whifile',whifile,'none')
   call read_char_input(500,'wtifile',wtifile,'none')
   call read_char_input(500,'wstfile',wstfile,'none')
   call read_char_input(500,'srcfile',srcfile,'none')
   call read_char_input(500,'disfile',disfile,'none')
   call read_char_input(500,'spwfile',spwfile,'none')
   call read_char_input(500,'wndfile',wndfile,'none')
   call read_char_input(500,'precipfile',prcpfile,'none')
   call read_char_input(500,'amufile',amufile,'none')
   call read_char_input(500,'amvfile',amvfile,'none')
   call read_char_input(500,'ampfile',ampfile,'none')
   call read_char_input(500,'amprfile',amprfile,'none')
   call read_char_input(500,'qinffile',qinffile,'none')
   call read_char_input(500,'scsfile',scsfile,'none')
   call read_char_input(500,'scsfile_Smax',scsfile_Smax,'none')
   call read_char_input(500,'scsfile_Se',scsfile_Se,'none')
   call read_char_input(500,'scsfile_kr',scsfile_kr,'none')
   call read_char_input(500,'z0lfile',z0lfile,'none')
   call read_char_input(500,'wvmfile',wvmfile,'none')
   !
   ! Netcdf input
   call read_char_input(500,'netbndbzsbzifile',netbndbzsbzifile,'none')  
   call read_char_input(500,'netsrcdisfile',netsrcdisfile,'none')  
   call read_char_input(500,'netamuamvfile',netamuamvfile,'none')                  
   call read_char_input(500,'netamprfile',netamprfile,'none')      
   call read_char_input(500,'netampfile',netampfile,'none')      
   !
   ! Output
   call read_char_input(500,'obsfile',obsfile,'none')
   call read_char_input(500,'crsfile',crsfile,'none')
   call read_int_input(500,'storevelmax',storevelmax,0)
   call read_int_input(500,'storevel',storevel,0)
   call read_int_input(500,'storecumprcp',storecumprcp,0)
   call read_int_input(500,'storetwet',storetwet,0)
   call read_int_input(500,'storehsubgrid',storehsubgrid,0)
   call read_real_input(500,'twet_threshold',twet_threshold,0.01)
   call read_int_input(500,'store_tsunami_arrival_time',itsunamitime,0)
   call read_real_input(500,'tsunami_arrival_threshold',tsunami_arrival_threshold,0.01)
   call read_int_input(500,'storeqdrain',storeqdrain,1)
   call read_int_input(500,'storezvolume',storezvolume,0)
   call read_int_input(500,'writeruntime',wrttimeoutput,0)
   call read_int_input(500,'debug',idebug,0)
   call read_int_input(500,'storemeteo',storemeteo,0)
   call read_int_input(500,'storemaxwind',iwindmax,0)
   call read_int_input(500,'storefw', istorefw, 0)
   call read_int_input(500,'storewavdir', istorewavdir, 0)
   !
   ! Wind drag
   !
   call read_int_input(500,'cdnrb',cd_nr,0)
   !
   if (cd_nr==0) then
      !
      ! Use defaults
      !
      cd_nr = 3
      !
      allocate(cd_wnd(cd_nr))
      allocate(cd_val(cd_nr))
      !
      cd_wnd(1) =   0.0
      cd_wnd(2) =  28.0
      cd_wnd(3) =  50.0
      cd_val(1) = 0.0010
      cd_val(2) = 0.0025
      cd_val(3) = 0.0015
      !
   else
      !
      ! Use defaults
      !
      call read_real_array_input(500,'cdwnd',cd_wnd,0.0,cd_nr)
      call read_real_array_input(500,'cdval',cd_val,0.0,cd_nr)
      !
   endif   
   !
   ! Try new keywords for sfincs.inp file (ensure backward compatibility)
   !
   if (dtmapout==0.0) then
      call read_real_input(500,'dtmapout',dtmapout,0.0)
   endif   
   !
   close(500)
   !
   ! Check whether epsg code has been specified:
   if (epsg == 0) then
       write(*,*)'WARNING: no epsg code defined' 
   endif   
   !
   ! Compute simulation time
   !
   call time_difference(trefstr,tstartstr,dtsec)  ! time difference in seconds between tstart and tref
   t0 = dtsec*1.0 ! time difference in seconds between tstop and tstart
   call time_difference(trefstr,tstopstr,dtsec)
   t1 = dtsec*1.0 ! time difference in seconds between tstop and tstart
   tspinup = t0 + tspinup
   !
   g         = 9.81
   pi        = 3.14159
   gn2       = 9.81*0.02*0.02 ! Only to be used in subgrid
   !
   qinf = qinf/(3600*1000)
   !
   rotation = rotation*pi/180
   cosrot   = cos(rotation)
   sinrot   = sin(rotation)
   !
   area  = dx*dy
   !
   dxy   = min(dx, dy)
   dxinv = 1.0/dx
   dyinv = 1.0/dy
   !
   manning2d = .false.
   imanning2d = 0
   if (manningfile/='none') then
      manning2d = .true. 
      imanning2d = 1
   endif   
   !
   ! Coriolis parameter
   !
   fcorio = 2*7.2921e-05*sin(latitude*pi/180)
   if (latitude>0.01 .or. latitude<-0.01) then
      coriolis = .true.
   else
      coriolis = .false.
   endif
   !
   use_coriolis = .true.
   crsgeo = .false.
   if (igeo==1) then
      coriolis = .true.
      crsgeo   = .true.
      if (icoriolis==1) then
         use_coriolis = .true.
         write(*,*)'Turning on process: Coriolis'               
      else
         use_coriolis = .false.
      endif   
      write(*,*)'Input grid interpreted as spherical coordinates, crsgeo= ',crsgeo
      !
   endif
   !
   ! Output
   !
   if (t0out<-900.0) then
      t0out = t0
   endif    
   t0out = max(t0out, t0)
   if (t1out<-900.0) then
      t1out = t1
   endif    
   !
   store_maximum_waterlevel = .false.
   if (dtmaxout>0.0) then
      store_maximum_waterlevel = .true.
   endif
   !
   store_maximum_velocity = .false.
   if (storevelmax==1 .and. dtmaxout>0.0) then
      store_maximum_velocity = .true.
   endif
   !
   store_velocity = .false.
   if (storevel==1) then
      store_velocity = .true.
   endif
   !
   store_hsubgrid = .false.
   if (storehsubgrid==1) then
      store_hsubgrid = .true.
   endif   
   !
   store_meteo = .false.
   store_wind_max = .false.
   if (storemeteo==1) then
      store_meteo = .true.
      if (iwindmax==1) then
         store_wind_max = .true.
      endif
   endif
   !
   store_twet = .false.
   if (storetwet==1) then
      store_twet = .true.
      twet      = 0.0
   endif
   !
   store_cumulative_precipitation = .false.
   if (storecumprcp==1) then
      store_cumulative_precipitation = .true.
   endif
   !   
   if (storeqdrain==0) then
      store_qdrain = .false.
   else
      store_qdrain = .true.
   endif
   !   
   write_time_output = .false.
   if (wrttimeoutput==1) then
      write_time_output = .true.
   endif
   !
   debug = .false.
   if (idebug==1) then
      debug = .true.
   endif   
   !
   radstr = .false.
   if (iradstr==1) then
      radstr = .true.
   endif   
   !
   if ((outputtype_map == 'nil') .OR. (outputtype_his == 'nil')) then
        outputtype_map = outputtype
        outputtype_his = outputtype
   endif
   !
   ampr_block = .true. ! Default use data in ampr file as block rather than linear interpolation
   if (iamprblock==0) then
      ampr_block = .false.
   endif 
   !
   global = .false. ! Default use data in ampr file as block rather than linear interpolation
   if (iglobal==1) then
      global = .true.
   endif   
   !
   if (sbgfile(1:4) /= 'none') then
      !
      subgrid = .true.
      isubgrid = 1
      write(*,*)'Info : Running SFINCS in subgrid mode ...'         
      !
   else
      !
      subgrid = .false.
      isubgrid = 0
      write(*,*)'Info : Running SFINCS in regular mode ...'               
      !
   endif
   !
   store_zvolume = .false.
   if (subgrid) then
       if (storezvolume==1) then
          store_zvolume = .true.
       endif
   endif
   !
   store_tsunami_arrival_time = .false. ! Default use data in ampr file as block rather than linear interpolation
   if (itsunamitime==1) then
      store_tsunami_arrival_time = .true.
   endif      
   !
   viscosity = .false.   
   if (nuviscinp>0.0 .or. iviscosity==1)then
      !
      if (nuviscinp>0.0) then ! if nuvisc given by user, this overrules the dimensionless default use of nuvisc, for backwards compatability
         !
         nuvisc = nuviscinp
         viscosity = .true.         
         !
      else ! user defined viscosity=1:
         !
         if (qtrfile(1:4) == 'none') then ! this works only for a regular grid model, for quadtree use 'nuvisc' option
            viscosity = .true.         
            !
            if (crsgeo == .true.) then ! simplified conversion for spherical grids
                nuvisc = max(nuviscdim * min(dx,dy) / 0.001, 0.0) ! take min of dx and dy, don't allow to be negative  
                ! dx = 1 degree ~ 100km
                ! dx = 0.001 degree~ 100m > nuvisc = 1.0              
            else                 
                nuvisc = max(nuviscdim * min(dx,dy) / 100.0, 0.0) ! take min of dx and dy, don't allow to be negative    
                ! dx = 50 > nuvisc = 0.5
                ! dx = 100 > nuvisc = 1.0
                ! dx = 500 > nuvisc = 5.0
                ! nuviscdim = 1.0
                ! nuvisc = nuviscdim * dx / 100 
            endif
         endif          
      endif    
      !
      if (viscosity == .true.) then
         write(*,*)'Turning on process: Viscosity, with nuvisc= ',nuvisc
      endif      
      !
   endif   
   !
   spinup_meteo = .true. ! Default use data in ampr file as block rather than linear interpolation
   if (ispinupmeteo==0) then
      spinup_meteo = .false.
   endif      
   !
   snapwave = .false.
   if (isnapwave==1) then
      snapwave = .true.
   endif      
   !
   fixed_output_intervals = .true.
   if (ioutfixed==0) then
      fixed_output_intervals = .false.
   endif      
   !
   advection = .false.
   if (iadvection>0) then
      advection = .true.
   endif      
   !
   store_wave_forces = .false.
   if (istorefw==1) then
      store_wave_forces = .true.
   endif      
   !
   wavemaker = .false.
   wavemaker_spectrum = .true.
   if (wvmfile(1:4) /= 'none') then
      wavemaker = .true.
      iwavemaker = 1
      if (wmsigstr(1:3) == 'mon') then
         ! Monochromatic
         wavemaker_spectrum = .false.
      endif   
   endif
   !
   store_wave_direction = .false.
   if (istorewavdir==1) then
      store_wave_direction = .true.
   endif      
   !
   end subroutine

   
   
   subroutine read_real_input(fileid,keyword,value,default)
   !
   character(*), intent(in) :: keyword
   character(len=256)       :: keystr
   character(len=256)       :: keystr0
   character(len=256)       :: valstr
   character(len=256)       :: line
   integer, intent(in)      :: fileid
   real*4, intent(out)      :: value
   real*4, intent(in)       :: default
   integer j,stat,ilen
   !
   value = default
   rewind(fileid)   
   do while(.true.)
      read(fileid,'(a)',iostat = stat)line
      if (stat==-1) exit
      j=index(line,'=')      
      keystr0 = trim(line(1:j-1))
      call notabs(keystr0,keystr,ilen)
      if (trim(keystr)==trim(keyword)) then
         valstr = trim(line(j+1:256))
         read(valstr,*)value
         exit
      endif
   enddo 
   !
   end  subroutine  

   subroutine read_real_array_input(fileid,keyword,value,default,nr)
   !
   character(*), intent(in) :: keyword
   character(len=256)       :: keystr0
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
   rewind(fileid)   
   do while(.true.)
      read(fileid,'(a)',iostat = stat)line
      if (stat==-1) exit
      j=index(line,'=')      
      keystr0 = trim(line(1:j-1))
      call notabs(keystr0,keystr,ilen)
      if (trim(keystr)==trim(keyword)) then
         valstr = trim(line(j+1:256))
         read(valstr,*)(value(m), m = 1, nr)
         exit
      endif
   enddo 
   !
   end  subroutine  

   
   subroutine read_int_input(fileid,keyword,value,default)
   !
   character(*), intent(in) :: keyword
   character(len=256)       :: keystr
   character(len=256)       :: keystr0
   character(len=256)       :: valstr
   character(len=256)       :: line
   integer, intent(in)      :: fileid
   integer, intent(out)     :: value
   integer, intent(in)      :: default
   integer j,stat,ilen
   !
   value = default
   rewind(fileid)   
   do while(.true.)
      read(fileid,'(a)',iostat = stat)line
      if (stat==-1) exit
      j=index(line,'=')      
      keystr0 = trim(line(1:j-1))
      call notabs(keystr0,keystr,ilen)
      if (trim(keystr)==trim(keyword)) then
         valstr = trim(line(j+1:256))
         read(valstr,*)value         
         exit
      endif
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
   integer j,stat,ilen
   !
   value = default
   rewind(fileid)   
   do while(.true.)
      read(fileid,'(a)',iostat = stat)line
      if (stat==-1) exit
      j=index(line,'=')      
      keystr0 = trim(line(1:j-1))
      call notabs(keystr0,keystr,ilen)
      if (trim(keystr)==trim(keyword)) then
         valstr = adjustl(trim(line(j+1:256)))
         value = valstr
         exit
      endif
   enddo 
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
