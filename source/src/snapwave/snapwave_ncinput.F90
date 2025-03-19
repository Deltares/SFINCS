#define NF90(nf90call) call handle_err(nf90call,__FILE__,__LINE__)   
module snapwave_ncinput
   !
   use sfincs_log
   use netcdf       
   !
   implicit none
   !
   type net_type_snapwave
       integer :: ncid
       integer :: points_dimid
       integer :: time_dimid
       integer :: x_varid, y_varid
       integer :: time_varid, hs_varid, tp_varid, wd_varid, ds_varid
   end type   
   !          
   type(net_type_snapwave)  :: net_file_snapwave              
   !
   contains   
    
   
    
   subroutine read_netcdf_wave_boundary_data()
   !
   use sfincs_date   
   use netcdf
   use snapwave_data      
   !
   implicit none   
   !
   ! Wanted variable names for snapwave netcdf input (FEWS style)
   !
   character (len=256)            :: x_varname
   character (len=256)            :: y_varname
   character (len=256), parameter :: time_varname = 'time'
   character (len=256), parameter :: hs_varname   = 'hs'
   character (len=256), parameter :: tp_varname   = 'tp'   
   character (len=256), parameter :: wd_varname   = 'wd'   
   character (len=256), parameter :: ds_varname   = 'ds'      
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
   write(logstr,*)'Reading netcdf type SnapWave input boundaries (bnd, bhs, btp, bwd, bds)...'
   call write_log(logstr, 1)   
   !
   ! Actual reading of data
   !
   NF90(nf90_open(trim(netsnapwavefile), NF90_CLOBBER, net_file_snapwave%ncid))
   !          
   ! Get dimensions id's: time, stations      
   NF90(nf90_inq_dimid(net_file_snapwave%ncid, "time", net_file_snapwave%time_dimid))
   NF90(nf90_inq_dimid(net_file_snapwave%ncid, "stations", net_file_snapwave%points_dimid))
   !
   ! Get dimensions sizes: time, stations      
   NF90(nf90_inquire_dimension(net_file_snapwave%ncid, net_file_snapwave%time_dimid, len = ntwbnd))   !nr of timesteps in file
   NF90(nf90_inquire_dimension(net_file_snapwave%ncid, net_file_snapwave%points_dimid, len = nwbnd))  !nr of boundary points     
   !
   ! Get variable id's
   NF90(nf90_inq_varid(net_file_snapwave%ncid, x_varname, net_file_snapwave%x_varid) )  ! Has to be in the same UTM zone as SnapWave grid
   NF90(nf90_inq_varid(net_file_snapwave%ncid, y_varname, net_file_snapwave%y_varid) )
   NF90(nf90_inq_varid(net_file_snapwave%ncid, time_varname, net_file_snapwave%time_varid) )
   NF90(nf90_inq_varid(net_file_snapwave%ncid, hs_varname, net_file_snapwave%hs_varid) )    
   NF90(nf90_inq_varid(net_file_snapwave%ncid, tp_varname, net_file_snapwave%tp_varid) )    
   NF90(nf90_inq_varid(net_file_snapwave%ncid, wd_varname, net_file_snapwave%wd_varid) )    
   NF90(nf90_inq_varid(net_file_snapwave%ncid, ds_varname, net_file_snapwave%ds_varid) )  
   !
   ! Allocate variables   
   !
   allocate(x_bwv(nwbnd))
   allocate(y_bwv(nwbnd)) 
   allocate(t_bwv(ntwbnd))
   allocate(hs_bwv(nwbnd, ntwbnd))
   allocate(tp_bwv(nwbnd, ntwbnd))
   allocate(wd_bwv(nwbnd, ntwbnd))
   allocate(ds_bwv(nwbnd, ntwbnd))
   !
   ! Read values
   NF90(nf90_get_var(net_file_snapwave%ncid, net_file_snapwave%x_varid, x_bwv(:)) )
   NF90(nf90_get_var(net_file_snapwave%ncid, net_file_snapwave%y_varid, y_bwv(:)) )
   NF90(nf90_get_var(net_file_snapwave%ncid, net_file_snapwave%time_varid, t_bwv(:)) ) 
   NF90(nf90_get_var(net_file_snapwave%ncid, net_file_snapwave%hs_varid, hs_bwv(:,:)) )   
   NF90(nf90_get_var(net_file_snapwave%ncid, net_file_snapwave%tp_varid, tp_bwv(:,:)) )   
   NF90(nf90_get_var(net_file_snapwave%ncid, net_file_snapwave%wd_varid, wd_bwv(:,:)) )   
   NF90(nf90_get_var(net_file_snapwave%ncid, net_file_snapwave%ds_varid, ds_bwv(:,:)) )      
   !
   ! Read time attibute
   NF90(nf90_get_att(net_file_snapwave%ncid, net_file_snapwave%time_varid, UNITS, treftimefews))
   !      
   ! Convert input time to sfincs time wrt treftime, time in same timezone as sfincs (UTC)
   t_bwv = convert_fewsdate(t_bwv, ntwbnd, treftimefews, trefstr)
   !
   NF90(nf90_close(net_file_snapwave%ncid))       
   ! 
   end subroutine
   
   
   
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
      end if
   end subroutine handle_err
   !
   end module    