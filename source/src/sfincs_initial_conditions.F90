module sfincs_initial_conditions
   !
   use sfincs_data
   use sfincs_log
   use sfincs_error
   use netcdf       
   !         
   real*8, dimension(:),   allocatable :: inizs
   real*4, dimension(:),   allocatable :: inizs4   
   real*4, dimension(:),   allocatable :: iniq
   !
contains
   !
   subroutine set_initial_conditions()
      !
      use sfincs_ncinput      
      !
      ! Initialize SFINCS variables (qx, qy, zs etc.)
      !
      implicit none
      !
      integer    :: nm, nmu, ip, iuv
      real*4     :: facint
      real*4     :: dzuv
      real*4     :: zmax
      real*4     :: zmin
      real*4     :: huv
      real*4     :: zsuv   
      ! 
      logical   :: iok
      character*256 :: varname      
      !
      allocate(inizs(np))
      allocate(inizs4(np))      
      allocate(iniq(npuv+ncuv+1)) ! Note: q is actually larger than npuv! It has size npuv + ncuv + 1
      !
      inizs = zini
      inizs4 = zini      
      iniq  = 0.0
      !
      ! Check the type of initial conditions
      !
      if (rstfile(1:4) /= 'none') then
         !
         ! Binary restart file
         ! Note - older type real*4 for zs 
         !
         write(logstr,'(a,a)')'Info    : reading restart file ', trim(rstfile)
         call write_log(logstr, 0)
         !
         iok = check_file_exists(rstfile, 'Restart file', .true.)
         !
         call read_binary_restart_file() ! Note - older type real*4 for zs
         !
      elseif (zsinifile(1:4) /= 'none') then 
         !
         ! Read binary (!) initial water level file
         ! Note - older type real*4 for zs 
         !
         write(logstr,'(a,a)')'Info    : reading initial conditions file ', trim(zsinifile)
         call write_log(logstr, 0)
         !
         iok = check_file_exists(zsinifile, 'Binary initial conditions ini file', .true.)
         !
         call read_zsini_file()
         !
      elseif (ncinifile(1:4) /= 'none') then 
         !
         ! Read netcdf (!) initial water level file
         ! Note - newer type real*8 for zs 
         !
         ! NetCDF file
         !
         write(logstr,'(a,a)')'Info    : reading NetCDF initial conditions file ', trim(zsinifile)
         call write_log(logstr, 0)
         !
         iok = check_file_exists(ncinifile, 'NetCDF initial conditions file', .true.)
         !
         !call read_nc_ini_file()
         ! Call the generic quadtree nc file reader function - in real*8
         varname = 'zs'
         call read_netcdf_quadtree_to_sfincs_real8(ncinifile, varname, inizs) !ncfile, varname, varout)          
         !
      else
         !
         ! No initial conditions provided
         !
      endif        
      !
      ! Water levels
      !
      do nm = 1, np
         !
         if (subgrid) then
            zs(nm) = max(subgrid_z_zmin(nm), inizs(nm)) ! Water level at zini or bed level (whichever is higher)
         else
            zs(nm) = max(zb(nm), inizs(nm)) ! Water level at zini or bed level (whichever is higher)
         endif
         !
      enddo      
      !
      ! Flux q
      !
      do ip = 1, npuv
         !
         q(ip) = iniq(ip)
         !
         ! Also need to compute initial uv
         ! Use same method as used in sfincs_momentum.f90
         !
         nm  = uv_index_z_nm(ip)
         nmu = uv_index_z_nmu(ip)
         !
         zsuv = max(zs(nm), zs(nmu)) ! water level at uv point 
         !
         iok = .false.
         !
         if (subgrid) then
            !
            zmin = subgrid_uv_zmin(ip)
            zmax = subgrid_uv_zmax(ip)
            !            
            if (zsuv>zmin + huthresh) then
               iok = .true.
            endif   
            !
         else
            !            
            if (zsuv>zbuvmx(ip)) then
               iok = .true.
            endif   
            !            
         endif   
         !
         if (iok) then
            !
            if (subgrid) then
               !
               if (zsuv>zmax - 1.0e-4) then
                  !
                  ! Entire cell is wet, no interpolation from table needed
                  !
                  huv    = subgrid_uv_havg_zmax(ip) + zsuv
                  !
               else
                  !
                  ! Interpolation required
                  !
                  dzuv   = (subgrid_uv_zmax(ip) - subgrid_uv_zmin(ip)) / (subgrid_nlevels - 1)
                  iuv    = int((zsuv - subgrid_uv_zmin(ip))/dzuv) + 1
                  facint = (zsuv - (subgrid_uv_zmin(ip) + (iuv - 1)*dzuv) ) / dzuv
                  huv    = subgrid_uv_havg(iuv, ip) + (subgrid_uv_havg(iuv + 1, ip) - subgrid_uv_havg(iuv, ip))*facint
                  !                      
               endif
               !
               huv    = max(huv, huthresh)
               !
            else
               !
               huv    = max(zsuv - zbuv(ip), huthresh)
               !
            endif
            !
            uv(ip)   = max(min(q(ip)/huv, 4.0), -4.0)   
            !
         endif   
         !
      enddo
      !
      deallocate(inizs)
      deallocate(iniq)
      !
   end subroutine
   !
   !
   ! 
   subroutine read_binary_restart_file()
      !
      ! Binary restart file
      !
      implicit none
      !
      integer    :: rsttype
      real*4     :: rdummy
      !
      open(unit = 500, file = trim(rstfile), form = 'unformatted', access = 'stream')
      !
      ! Restartfile flavours:
      ! 1: zs, q, uvmean  
      ! 2: zs, q 
      ! 3: zs  - 
      ! 4: zs, q, uvmean and cnb infiltration (writing scs_Se)
      ! 5: zs, q, uvmean and gai infiltration (writing GA_sigma & GA_F)
      ! 6: zs, q, uvmean and hor infiltration (writing rain_T1)         
      !
      read(500)rdummy
      read(500)rsttype
      read(500)rdummy
      !
      write(logstr,'(a,i0)')'Info    : found rsttype = ', rsttype
      call write_log(logstr, 0)
      !
      if (rsttype < 1 .or. rsttype > 6) then
         !
         ! Give warning, rstfile input rsttype not recognized
         !
         write(logstr,'(a,i0)')'Warning! rstfile not recognized, skipping restart file input! rsttype should be 1-6, but found rsttype = ', rsttype 
         call write_log(logstr, 1)
         !          
         close(500)      
         !          
      else
         ! Always read in inizs
         !
         read(500)rdummy
         read(500)inizs4
         read(500)rdummy
         !      
         ! Read fluxes q
         !
         if (rsttype==1 .or. rsttype==2 .or. rsttype==4 .or. rsttype==5 .or. rsttype==6) then     
            read(500)rdummy
            read(500)iniq
            read(500)rdummy
            !
            read(500)rdummy
            read(500)uvmean
            read(500)rdummy
         endif
         !
         if (rsttype==4) then ! Infiltration method cnb 
            !
            read(500)rdummy                 
            read(500)scs_Se
            ! write(*,*)'Reading scs_Se from rstfile, overwrites input values of: ',trim(sefffile)
            !
         elseif (rsttype==5) then ! Infiltration method gai    
            !
            read(500)rdummy
            read(500)GA_sigma
            read(500)rdummy
            read(500)GA_F
            write(logstr,'(a,a)')'Info    : reading GA_sigma from rstfile, overwrites input values of ', trim(sigmafile)
            call write_log(logstr, 0)
            !
         elseif (rsttype==6) then ! Infiltration method horton
            !
            read(500)rdummy                               
            read(500)rain_T1
            write(logstr,'(a,a)')'Info    : reading rain_T1 from rstfile, complements input values of ', trim(fcfile) 
            call write_log(logstr, 0)
            !              
         endif          
         !
         close(500)      
         !
         ! remap zs from real*4 to real*8
         inizs = inizs4
         !
      endif
      !
   end subroutine
   !
   subroutine read_zsini_file()
      !
      ! Binary zs ini file
      !
      implicit none
      !
      write(logstr,'(a,a)')'Info    : reading zsini file ', trim(zsinifile)
      call write_log(logstr, 0)
      !
      call write_log('Warning : binary ini files from SFINCS v2.1.1 and older are not compatible with SFINCS v2.1.2+, remake your inifile containing zs as real*8 double precision', 0)
      !
      open(unit = 500, file = trim(zsinifile), form = 'unformatted', access = 'stream')
      read(500)inizs4
      close(500)       
      !
      ! remap from real*4 to real*8
      inizs = inizs4
      !
   end subroutine
   !
end module
