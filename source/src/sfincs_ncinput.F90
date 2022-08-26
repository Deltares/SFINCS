#define NF90(nf90call) call handle_err(nf90call,__FILE__,__LINE__)   
module sfincs_ncinput
   !
   use netcdf       
   !
   implicit none
   !
   type net_type_bndbzsbzi
       integer :: ncid
       integer :: points_dimid
       integer :: time_dimid
       integer :: x_varid, y_varid
       integer :: time_varid, zs_varid, zi_varid
   end type
   type net_type_srcdis
       integer :: ncid
       integer :: points_dimid
       integer :: time_dimid
       integer :: x_varid, y_varid
       integer :: time_varid, q_varid
   end type
   type net_type_amuv
       integer :: ncid
       integer :: ncols_dimid, nrows_dimid
       integer :: time_dimid
       integer :: wx_varid, wy_varid
       integer :: time_varid
       integer :: wu_varid, wv_varid      
   end type   
   type net_type_amp
       integer :: ncid
       integer :: ncols_dimid, nrows_dimid
       integer :: time_dimid
       integer :: px_varid, py_varid
       integer :: time_varid
       integer :: patm_varid     
   end type      
   type net_type_ampr
       integer :: ncid
       integer :: ncols_dimid, nrows_dimid
       integer :: time_dimid
       integer :: px_varid, py_varid
       integer :: time_varid
       integer :: prcp_varid     
   end type      
   !
   type(net_type_bndbzsbzi) :: net_file_bndbzsbzi        
   type(net_type_srcdis)    :: net_file_srcdis        
   type(net_type_amuv)      :: net_file_amuv     
   type(net_type_amp)       :: net_file_amp 
   type(net_type_ampr)      :: net_file_ampr              

   contains   
   
   subroutine read_netcdf_boundary_data()
   !
   use sfincs_date   
   use netcdf
   use sfincs_data   
   !
   implicit none   
   !
   ! Wanted variable names for Fews compatible netcdf input
   !
   character (len=256)            :: x_varname
   character (len=256)            :: y_varname
   character (len=256), parameter :: time_varname = 'time'
   character (len=256), parameter :: zs_varname   = 'zs'
   character (len=256), parameter :: zi_varname   = 'zi'   
   character (len=256), parameter :: units        = 'units'  
   !
!   if (crsgeo) then
!      x_varname    = 'lon'
!      y_varname    = 'lat'  
!   else
      x_varname    = 'x'
      y_varname    = 'y'  
!   endif   
   !
   write(*,*)'Reading FEWS compatible netcdf type water level boundaries (bnd, bzs, bzi)...'
   !
   ! Actual reading of data
   !
   NF90(nf90_open(trim(netbndbzsbzifile), NF90_CLOBBER, net_file_bndbzsbzi%ncid))
   !          
   ! Get dimensions id's: time, stations      
   NF90(nf90_inq_dimid(net_file_bndbzsbzi%ncid, "time", net_file_bndbzsbzi%time_dimid))
   NF90(nf90_inq_dimid(net_file_bndbzsbzi%ncid, "stations", net_file_bndbzsbzi%points_dimid))
   !
   ! Get dimensions sizes: time, stations      
   NF90(nf90_inquire_dimension(net_file_bndbzsbzi%ncid, net_file_bndbzsbzi%time_dimid, len = ntbnd))   !nr of timesteps in file
   NF90(nf90_inquire_dimension(net_file_bndbzsbzi%ncid, net_file_bndbzsbzi%points_dimid, len = nbnd))  !nr of boundary points     
   !
   ! Get variable id's
   NF90(nf90_inq_varid(net_file_bndbzsbzi%ncid, x_varname, net_file_bndbzsbzi%x_varid) )  ! Has to be in the same UTM zone as SFINCS grid
   NF90(nf90_inq_varid(net_file_bndbzsbzi%ncid, y_varname, net_file_bndbzsbzi%y_varid) )
   NF90(nf90_inq_varid(net_file_bndbzsbzi%ncid, time_varname, net_file_bndbzsbzi%time_varid) )
   NF90(nf90_inq_varid(net_file_bndbzsbzi%ncid, zs_varname, net_file_bndbzsbzi%zs_varid) )    
   !try: block
   NF90(nf90_inq_varid(net_file_bndbzsbzi%ncid, zi_varname, net_file_bndbzsbzi%zi_varid) )   ! add try statement ?       
   !     exit try ! Exit the try block in case of normal execution
   !     !continue ! Jump here in case of an error
   !     write(*,*)'No bziwaves specified in netbndbzsbzifile'
   !end block try
   !
   ! Allocate variables   
   allocate(x_bnd(nbnd))
   allocate(y_bnd(nbnd)) 
   allocate(t_bnd(ntbnd))
   allocate(zs_bnd(nbnd,ntbnd))
   allocate(zst_bnd(nbnd))  
   !
   ! Read values
   NF90(nf90_get_var(net_file_bndbzsbzi%ncid, net_file_bndbzsbzi%x_varid, x_bnd(:)) )
   NF90(nf90_get_var(net_file_bndbzsbzi%ncid, net_file_bndbzsbzi%y_varid, y_bnd(:)) )
   NF90(nf90_get_var(net_file_bndbzsbzi%ncid, net_file_bndbzsbzi%time_varid, t_bnd(:)) ) 
   NF90(nf90_get_var(net_file_bndbzsbzi%ncid, net_file_bndbzsbzi%zs_varid, zs_bnd(:,:)) )   
   !   
   if (net_file_bndbzsbzi%zi_varid /= 0) then     ! Allocate and read if bzi data is present
      bziwaves = .true.   ! turn on flag
      write(*,*)'   Incident waves are forced...'      
      !      
      allocate(zsi_bnd(nbnd,ntbnd))
      allocate(zsit_bnd(nbnd))
      !      
      NF90(nf90_get_var(net_file_bndbzsbzi%ncid, net_file_bndbzsbzi%zi_varid, zsi_bnd(:,:)) )              
   endif      
   !
   ! Read time attibute
   NF90(nf90_get_att(net_file_bndbzsbzi%ncid, net_file_bndbzsbzi%time_varid, UNITS, treftimefews))
   !      
   ! Convert input time to sfincs time wrt treftime, time in same timezone as sfincs (UTC)
   t_bnd = convert_fewsdate(t_bnd, ntbnd, treftimefews, trefstr)
   !
   NF90(nf90_close(net_file_bndbzsbzi%ncid))       
   ! 
   end subroutine

   
   
   subroutine read_netcdf_discharge_data()
   !
   use sfincs_date   
   use netcdf
   use sfincs_data   
   !
   implicit none   
   !
   ! Wanted variable names for Fews compatible netcdf input
   !
   character (len=256)            :: x_varname
   character (len=256)            :: y_varname
   character (len=256), parameter :: time_varname = 'time'
   character (len=256), parameter :: q_varname    = 'discharge'
   character (len=256), parameter :: units        = 'units'  
   !
!   if (crsgeo) then
!      x_varname    = 'lon'
!      y_varname    = 'lat'  
!   else
      x_varname    = 'x'
      y_varname    = 'y'  
!   endif   
   !
   write(*,*)'Reading FEWS compatible netcdf type discharges ...'
   !
   ! Actual reading of data
   !
   NF90(nf90_open(trim(netsrcdisfile), NF90_CLOBBER, net_file_srcdis%ncid))
   !          
   ! Get dimensions id's: time, stations      
   NF90(nf90_inq_dimid(net_file_srcdis%ncid, "time", net_file_srcdis%time_dimid))
   NF90(nf90_inq_dimid(net_file_srcdis%ncid, "stations", net_file_srcdis%points_dimid))
   !
   ! Get dimensions sizes: time, stations      
   NF90(nf90_inquire_dimension(net_file_srcdis%ncid, net_file_srcdis%time_dimid,   len = ntsrc))   !nr of timesteps in file
   NF90(nf90_inquire_dimension(net_file_srcdis%ncid, net_file_srcdis%points_dimid, len = nsrc))  !nr of discharge points     
   !
   ! Get variable id's
   NF90(nf90_inq_varid(net_file_srcdis%ncid, x_varname,    net_file_srcdis%x_varid) )  ! Has to be in the same UTM zone as SFINCS grid
   NF90(nf90_inq_varid(net_file_srcdis%ncid, y_varname,    net_file_srcdis%y_varid) )
   NF90(nf90_inq_varid(net_file_srcdis%ncid, time_varname, net_file_srcdis%time_varid) )
   NF90(nf90_inq_varid(net_file_srcdis%ncid, q_varname,    net_file_srcdis%q_varid) )    
   !
   ! Allocate variables   
   allocate(xsrc(nsrc))
   allocate(ysrc(nsrc)) 
   allocate(tsrc(ntsrc))
   allocate(qsrc(nsrc,ntsrc))
   !
   ! Read values
   NF90(nf90_get_var(net_file_srcdis%ncid, net_file_srcdis%x_varid,    xsrc(:)) )
   NF90(nf90_get_var(net_file_srcdis%ncid, net_file_srcdis%y_varid,    ysrc(:)) )
   NF90(nf90_get_var(net_file_srcdis%ncid, net_file_srcdis%time_varid, tsrc(:)) ) 
   NF90(nf90_get_var(net_file_srcdis%ncid, net_file_srcdis%q_varid,    qsrc(:,:)) )   
   !   
   ! Read time attibute
   !
   NF90(nf90_get_att(net_file_srcdis%ncid, net_file_srcdis%time_varid, UNITS, treftimefews))
   !      
   ! Convert input time to sfincs time wrt treftime, time in same timezone as sfincs (UTC)
   !
   tsrc = convert_fewsdate(tsrc, ntsrc, treftimefews, trefstr)
   !
   NF90(nf90_close(net_file_srcdis%ncid))       
   ! 
   end subroutine

   
   
   subroutine read_netcdf_amuv_data()
   ! Output is made exactly the same as original read_amuv_dimensions & read_amuv_file subroutines but then with data given by netcdf file
   !
   use sfincs_date   
   use netcdf
   use sfincs_data   
   !
   implicit none   
   !
   integer nt   
   !
   real*4, dimension(:,:,:),   allocatable :: wu 
   real*4, dimension(:,:,:),   allocatable :: wv
   real*4, dimension(:,:,:),   allocatable :: wutmp 
   real*4, dimension(:,:,:),   allocatable :: wvtmp   
   real*4, dimension(:,:,:),   allocatable :: amuv_wutmp 
   real*4, dimension(:,:,:),   allocatable :: amuv_wvtmp     
   !
   real*4, dimension(:),     allocatable :: wx
   real*4, dimension(:),     allocatable :: wy
   !
   ! Wanted variable names for Fews compatible netcdf input for amu/amv
   !
   character (len=256)            :: x_varname
   character (len=256)            :: y_varname
   character (len=256), parameter :: time_varname   = 'time'
   character (len=256), parameter :: wv_varname     = 'northward_wind'    
   character (len=256), parameter :: wu_varname     = 'eastward_wind'     
   character (len=256), parameter :: units          = 'units'     
   !
!   if (crsgeo) then
!      x_varname    = 'lon'
!      y_varname    = 'lat'  
!   else
      x_varname    = 'x'
      y_varname    = 'y'  
!   endif   
   !
   write(*,*)'Reading FEWS compatible NetCDF type wind input ...'
   !
   ! Actual reading of data
   !
   NF90(nf90_open(trim(netamuamvfile), NF90_CLOBBER, net_file_amuv%ncid))     
   !
   ! Get dimensions id's: time, stations      
   NF90(nf90_inq_dimid(net_file_amuv%ncid, "time",    net_file_amuv%time_dimid))
   NF90(nf90_inq_dimid(net_file_amuv%ncid, x_varname, net_file_amuv%ncols_dimid))
   NF90(nf90_inq_dimid(net_file_amuv%ncid, y_varname, net_file_amuv%nrows_dimid))
   !
   ! Get dimensions sizes: time, cols, rows      
   NF90(nf90_inquire_dimension(net_file_amuv%ncid, net_file_amuv%time_dimid,  len = amuv_nt))     ! nr of timesteps in file
   NF90(nf90_inquire_dimension(net_file_amuv%ncid, net_file_amuv%ncols_dimid, len = amuv_ncols))  ! nr of cols    
   NF90(nf90_inquire_dimension(net_file_amuv%ncid, net_file_amuv%nrows_dimid, len = amuv_nrows))  ! nr of rows    
   !
   ! Get variable id's
   NF90(nf90_inq_varid(net_file_amuv%ncid, x_varname,    net_file_amuv%wx_varid) )  ! Has to be in the same UTM zone as SFINCS grid
   NF90(nf90_inq_varid(net_file_amuv%ncid, y_varname,    net_file_amuv%wy_varid) )  ! Also, has to be a rectilinear raster, not curvilinear! 
   NF90(nf90_inq_varid(net_file_amuv%ncid, time_varname, net_file_amuv%time_varid) )
   NF90(nf90_inq_varid(net_file_amuv%ncid, wu_varname,   net_file_amuv%wu_varid) )      
   NF90(nf90_inq_varid(net_file_amuv%ncid, wv_varname,   net_file_amuv%wv_varid) )         
   !
   ! Allocate
   !
   allocate(amuv_times(amuv_nt))
   allocate(wu(amuv_ncols, amuv_nrows, amuv_nt))
   allocate(wv(amuv_ncols, amuv_nrows, amuv_nt))    
   allocate(wutmp(amuv_ncols, amuv_nrows, 1))
   allocate(wvtmp(amuv_ncols, amuv_nrows, 1))       
   allocate(amuv_wu(amuv_nt, amuv_nrows, amuv_ncols))
   allocate(amuv_wv(amuv_nt, amuv_nrows, amuv_ncols))   
   allocate(amuv_wutmp(1, amuv_nrows, amuv_ncols))
   allocate(amuv_wvtmp(1, amuv_nrows, amuv_ncols))      
   !
   allocate(wx(amuv_ncols))
   allocate(wy(amuv_nrows))     
   !
   ! Read values 
   NF90(nf90_get_var(net_file_amuv%ncid, net_file_amuv%time_varid, amuv_times(:)))
   NF90(nf90_get_var(net_file_amuv%ncid, net_file_amuv%wx_varid, wx(:)))
   NF90(nf90_get_var(net_file_amuv%ncid, net_file_amuv%wy_varid, wy(:)))   
   !
   do nt = 1, amuv_nt 
      !      
      ! Read wind data per time frame to reduce chance of memory errors       
      !
      NF90(nf90_get_var(net_file_amuv%ncid, net_file_amuv%wu_varid, wutmp, start = (/ 1, 1, nt /), count = (/ amuv_ncols, amuv_nrows, 1 /))) ! be aware of start indices
      NF90(nf90_get_var(net_file_amuv%ncid, net_file_amuv%wv_varid, wvtmp, start = (/ 1, 1, nt /), count = (/ amuv_ncols, amuv_nrows, 1 /))) ! be aware of start indices      
      !                 
      ! First reshape wind input in the dimensions nt, nrows, ncols        
      !
      amuv_wutmp = reshape( wutmp, (/ 1, amuv_nrows, amuv_ncols /), ORDER = (/ 3, 2, 1 /))     
      amuv_wvtmp = reshape( wvtmp, (/ 1, amuv_nrows, amuv_ncols /), ORDER = (/ 3, 2, 1 /))                  
      !
      ! Then flip order of rows so it is 'nrows' to 1 as in reading of original ampr file instead of 1 to 'nrows'
      !
      amuv_wutmp = amuv_wutmp(:,amuv_nrows:1:-1,:)
      amuv_wvtmp = amuv_wvtmp(:,amuv_nrows:1:-1,:)      
      !
      amuv_wu(nt,:,:) = amuv_wutmp(1,:,:)  
      amuv_wv(nt,:,:) = amuv_wvtmp(1,:,:)            
      !
   enddo                   
   !
   ! Determine other needed characteristics from grid input
   !
   amuv_x_llcorner = wx(1)
   amuv_y_llcorner = wy(1)
   !
   amuv_dx = wx(2) - wx(1)   ! Input grid has to be a rectilinear, not curvilinear! So amuv_dx/dy is constant
   amuv_dy = wy(2) - wy(1)   
   !
   ! Read time attibute
   !
   NF90(nf90_get_att(net_file_amuv%ncid, net_file_amuv%time_varid, UNITS, treftimefews))
   !      
   ! Convert input time to sfincs time wrt treftime, time in same timezone as sfincs (UTC)
   !
   amuv_times = convert_fewsdate(amuv_times, amuv_nt, treftimefews, trefstr)  
   !
   NF90(nf90_close(net_file_amuv%ncid))   
   !
   end subroutine

   
   
   subroutine read_netcdf_amp_data()
   ! 
   ! Output is made exactly the same as original read_amp_dimensions & read_amp_file subroutines but then with data given by netcdf file
   ! 
   use sfincs_date   
   use netcdf
   use sfincs_data   
   !
   implicit none   
   !
   integer nt
   !
   real*4, dimension(:,:,:),   allocatable :: patmtmp 
   real*4, dimension(:,:,:),   allocatable :: amp_patmtmp 
   !
   real*4, dimension(:),     allocatable :: px
   real*4, dimension(:),     allocatable :: py
   !
   ! Wanted variable names for Fews compatible netcdf input for prcp
   !
   character (len=256)            :: x_varname
   character (len=256)            :: y_varname
   character (len=256), parameter :: time_varname   = 'time'
   character (len=256), parameter :: patm_varname   = 'barometric_pressure'
   character (len=256), parameter :: units          = 'units'     
   !
!   if (crsgeo) then
!      x_varname    = 'lon'
!      y_varname    = 'lat'  
!   else
      x_varname    = 'x'
      y_varname    = 'y'  
!   endif   
   !
   write(*,*)'Reading FEWS compatible NetCDF type barometric pressure input ...'
   !
   ! Actual reading of data
   !
   NF90(nf90_open(trim(netampfile), NF90_CLOBBER, net_file_amp%ncid))     
   !
   ! Get dimensions id's: time, stations      
   !         
   NF90(nf90_inq_dimid(net_file_amp%ncid, time_varname,  net_file_amp%time_dimid))
   NF90(nf90_inq_dimid(net_file_amp%ncid, x_varname,     net_file_amp%ncols_dimid))
   NF90(nf90_inq_dimid(net_file_amp%ncid, y_varname,     net_file_amp%nrows_dimid))
   !
   ! Get dimensions sizes: time, cols, rows      
   !         
   NF90(nf90_inquire_dimension(net_file_amp%ncid, net_file_amp%time_dimid,  len = amp_nt))      ! nr of timesteps in file
   NF90(nf90_inquire_dimension(net_file_amp%ncid, net_file_amp%ncols_dimid, len = amp_ncols))   ! nr of cols    
   NF90(nf90_inquire_dimension(net_file_amp%ncid, net_file_amp%nrows_dimid, len = amp_nrows))   ! nr of rows    
   !
   ! Get variable id's
   !
   NF90(nf90_inq_varid(net_file_amp%ncid, x_varname,    net_file_amp%px_varid) )  ! Has to be the same as SFINCS grid
   NF90(nf90_inq_varid(net_file_amp%ncid, y_varname,    net_file_amp%py_varid) )  ! Also, has to be a rectilinear raster, not curvilinear! 
   NF90(nf90_inq_varid(net_file_amp%ncid, time_varname, net_file_amp%time_varid) )
   NF90(nf90_inq_varid(net_file_amp%ncid, patm_varname, net_file_amp%patm_varid) )      
   !   
   ! Allocate
   !
   allocate(amp_times(amp_nt))
   allocate(patmtmp(amp_ncols, amp_nrows, 1))
   allocate(amp_patm(amp_nt, amp_nrows, amp_ncols))   
   allocate(amp_patmtmp(1, amp_nrows, amp_ncols))   
   !
   allocate(px(amp_ncols))
   allocate(py(amp_nrows))   
   !
   ! Read values 
   !
   NF90(nf90_get_var(net_file_amp%ncid, net_file_amp%time_varid, amp_times(:)))
   NF90(nf90_get_var(net_file_amp%ncid, net_file_amp%px_varid, px(:)))
   NF90(nf90_get_var(net_file_amp%ncid, net_file_amp%py_varid, py(:)))
   !
   do nt = 1, amp_nt 
      !      
      ! Read precipitation data per time frame to reduce chance of memory errors       
      !         
      NF90(nf90_get_var(net_file_amp%ncid, net_file_amp%patm_varid, patmtmp, start = (/ 1, 1, nt /), count = (/ amp_ncols, amp_nrows, 1 /))) ! be aware of start indices
      !                 
      ! First reshape prcp input in the dimensions nt, nrows, ncols        
      amp_patmtmp = reshape( patmtmp, (/ 1, amp_nrows, amp_ncols /), ORDER = (/ 3, 2, 1 /))            
      !
      ! Then flip order of rows so it is 'nrows' to 1 as in reading of original amp file instead of 1 to 'nrows'
      !
      amp_patmtmp = amp_patmtmp(:,amp_nrows:1:-1,:)
      !
      amp_patm(nt,:,:) = amp_patmtmp(1,:,:)      
      !
   enddo
   !
   ! Determine other needed characteristics from grid input
   !
   amp_x_llcorner = px(1)
   amp_y_llcorner = py(1)
   !
   amp_dx = px(2) - px(1)   ! Input grid has to be a rectilinear, not curvilinear! So amp_dx/dy is constant
   amp_dy = py(2) - py(1)   
   !   
   ! Read time attibute
   !
   NF90(nf90_get_att(net_file_amp%ncid, net_file_amp%time_varid, UNITS, treftimefews))
   !         
   ! Convert input time to sfincs time wrt treftime, time in same timezone as sfincs (UTC)
   !         
   amp_times = convert_fewsdate(amp_times, amp_nt, treftimefews, trefstr)
   !       
   NF90(nf90_close(net_file_amp%ncid))     
   !
   end subroutine

   
   
   subroutine read_netcdf_ampr_data()
   ! 
   ! Output is made exactly the same as original read_amuv_dimensions & read_amuv_file subroutines but then with data given by netcdf file
   ! 
   use sfincs_date   
   use netcdf
   use sfincs_data   
   !
   implicit none   
   !
   integer nt
   !
   real*4, dimension(:,:,:),   allocatable :: prtmp 
   real*4, dimension(:,:,:),   allocatable :: ampr_prtmp 
   !
   real*4, dimension(:),     allocatable :: px
   real*4, dimension(:),     allocatable :: py
   !
   ! Wanted variable names for Fews compatible netcdf input for prcp
   !
   character (len=256)            :: x_varname
   character (len=256)            :: y_varname
   character (len=256), parameter :: time_varname   = 'time'
   character (len=256), parameter :: prcp_varname   = 'Precipitation'
   character (len=256), parameter :: units          = 'units'     
   !
!   if (crsgeo) then
!      x_varname    = 'lon'
!      y_varname    = 'lat'  
!   else
      x_varname    = 'x'
      y_varname    = 'y'  
!   endif   
   !
   write(*,*)'Reading FEWS compatible NetCDF type precipitation input ...'
   !
   ! Actual reading of data
   !
   NF90(nf90_open(trim(netamprfile), NF90_CLOBBER, net_file_ampr%ncid))     
   !
   ! Get dimensions id's: time, stations
   !
   NF90(nf90_inq_dimid(net_file_ampr%ncid, "time",    net_file_ampr%time_dimid))
   NF90(nf90_inq_dimid(net_file_ampr%ncid, x_varname, net_file_ampr%ncols_dimid))
   NF90(nf90_inq_dimid(net_file_ampr%ncid, y_varname, net_file_ampr%nrows_dimid))
   !
   ! Get dimensions sizes: time, cols, rows      
   !
   NF90(nf90_inquire_dimension(net_file_ampr%ncid, net_file_ampr%time_dimid,  len = ampr_nt))      !nr of timesteps in file
   NF90(nf90_inquire_dimension(net_file_ampr%ncid, net_file_ampr%ncols_dimid, len = ampr_ncols))   !nr of cols    
   NF90(nf90_inquire_dimension(net_file_ampr%ncid, net_file_ampr%nrows_dimid, len = ampr_nrows))   !nr of rows    
   !
   ! Get variable id's
   !
   NF90(nf90_inq_varid(net_file_ampr%ncid, x_varname,    net_file_ampr%px_varid) )  ! Has to be same as SFINCS grid
   NF90(nf90_inq_varid(net_file_ampr%ncid, y_varname,    net_file_ampr%py_varid) )  ! Also, has to be a rectilinear raster, not curvilinear! 
   NF90(nf90_inq_varid(net_file_ampr%ncid, time_varname, net_file_ampr%time_varid) )
   NF90(nf90_inq_varid(net_file_ampr%ncid, prcp_varname, net_file_ampr%prcp_varid) )      
   !   
   ! Allocate
   !
   allocate(ampr_times(ampr_nt))
   allocate(prtmp(ampr_ncols, ampr_nrows,1))
   allocate(ampr_pr(ampr_nt, ampr_nrows, ampr_ncols))   
   allocate(ampr_prtmp(1, ampr_nrows, ampr_ncols))   
   !
   allocate(px(ampr_ncols))
   allocate(py(ampr_nrows))   
   !
   ! Read values 
   !
   NF90(nf90_get_var(net_file_ampr%ncid, net_file_ampr%time_varid, ampr_times(:)))
   NF90(nf90_get_var(net_file_ampr%ncid, net_file_ampr%px_varid, px(:)))
   NF90(nf90_get_var(net_file_ampr%ncid, net_file_ampr%py_varid, py(:)))
   !
   do nt = 1, ampr_nt 
      !      
      ! Read precipitation data per time frame to reduce chance of memory errors       
      NF90(nf90_get_var(net_file_ampr%ncid, net_file_ampr%prcp_varid, prtmp, start = (/ 1, 1, nt /), count = (/ ampr_ncols, ampr_nrows, 1 /))) ! be aware of start indices
      !                 
      ! First reshape prcp input in the dimensions nt, nrows, ncols        
      ampr_prtmp = reshape( prtmp, (/ 1, ampr_nrows, ampr_ncols /), ORDER = (/ 3, 2, 1 /))            
      !
      ! Then flip order of rows so it is 'nrows' to 1 as in reading of original ampr file instead of 1 to 'nrows'
      ampr_prtmp = ampr_prtmp(:,ampr_nrows:1:-1,:)
      !
      ampr_pr(nt,:,:) = ampr_prtmp(1,:,:)      
      !
   enddo                   
   ! Determine other needed characteristics from grid input
   ampr_x_llcorner = px(1)
   ampr_y_llcorner = py(1)
   !
   ampr_dx = px(2) - px(1)   ! Input grid has to be a rectilinear, not curvilinear! So ampr_dx/dy is constant
   ampr_dy = py(2) - py(1)   
   !   
   ! Read time attibute
   NF90(nf90_get_att(net_file_ampr%ncid, net_file_ampr%time_varid, UNITS, treftimefews))
   !         
   ! Convert input time to sfincs time wrt treftime, time in same timezone as sfincs (UTC)
   ampr_times = convert_fewsdate(ampr_times, ampr_nt, treftimefews, trefstr)
   !       
   NF90(nf90_close(net_file_ampr%ncid))     
   !
   end subroutine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine handle_err(status,file,line)
      !
      integer, intent ( in)    :: status
      character(*), intent(in) :: file
      integer, intent ( in)    :: line
      integer :: status2
   
      if(status /= nf90_noerr) then
      !   !UNIT=6 for stdout and UNIT=0 for stderr.
         write(0,'("NETCDF ERROR: ",a,i6,":",a)') file,line,trim(nf90_strerror(status))
      !   write(0,*) 'closing file'
      !   !status2 = nf90_close(net_file_bndbzs%ncid) ! even uit gezet omdat hij net_file_bndbzs en ncid niet kent omdat het locale variabelen zijn
      !   !if (status2 /= nf90_noerr) then
      !   !   write(0,*) 'NETCDF ERROR: ', __FILE__,__LINE__,trim(nf90_strerror(status2))
      !   !end if
      end if
   end subroutine handle_err
   !
   end module