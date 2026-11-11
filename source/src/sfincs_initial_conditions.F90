module sfincs_initial_conditions
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! Initial conditions: sets the initial water levels and fluxes, from sfincs.inp, an
   ! initial water level file or a restart file.
   !
   ! Subroutines / functions:
   !
   !   set_initial_conditions()
   !     Set the initial water levels, fluxes and velocities for the chosen type of initial
   !     condition.
   !     Called from: initialize_hydro (sfincs_domain)
   !
   !   read_binary_restart_file()
   !     Read the binary restart file.
   !     Called from: set_initial_conditions (sfincs_initial_conditions)
   !
   !   read_zsini_file()
   !     Read the binary initial water level file.
   !     Called from: set_initial_conditions (sfincs_initial_conditions)
   !
   !   set_restart_infiltration()
   !     Copy the infiltration state from the restart file into the infiltration arrays.
   !     Called from: sfincs_initialize (sfincs_lib)
   !
   ! ------------ Copyright notice SFINCS ------------
   !
   ! Copyright (C)  Stichting Deltares, 2016-2026.
   !
   ! This file is part of SFINCS.
   !
   ! SFINCS is free software: you can redistribute it and/or modify
   ! it under the terms of the GNU General Public License as
   ! published by the Free Software Foundation version 3.
   !
   ! SFINCS is distributed in the hope that it will be useful,
   ! but WITHOUT ANY WARRANTY; without even the implied warranty of
   ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   ! GNU General Public License for more details.
   !
   ! You should have received a copy of the GNU General Public License
   ! along with SFINCS.  If not, see <http://www.gnu.org/licenses/>.
   !
   ! contact: software@deltares.nl
   ! Stichting Deltares
   ! P.O. Box 177
   ! 2600 MH Delft, The Netherlands
   !
   ! All indications, logos of, and references to, "SFINCS"
   ! remain the property of Stichting Deltares. All rights reserved.
   !
   ! ------------ Copyright notice SFINCS ------------
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
   ! Infiltration state from the restart file. The restart file is read before the infiltration arrays are
   ! allocated in initialize_infiltration, so the state is kept here until set_restart_infiltration is called.
   !
   integer                             :: rst_inftype = 0 ! rsttype with infiltration state (4, 5 or 6), 0 if none
   real*4, dimension(:),   allocatable :: rst_inf1        ! scs_Se (4), GA_sigma (5) or rain_T1 (6)
   real*4, dimension(:),   allocatable :: rst_inf2        ! GA_F (5)
   !
contains
   !
   subroutine set_initial_conditions()
      !
      ! Set the initial water levels, fluxes and velocities for the chosen type of initial
      ! condition.
      !
      ! Called from: initialize_hydro (sfincs_domain)
      !
      use sfincs_ncinput      
      use sfincs_read, only: is_netcdf_file
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
         write(logstr,'(a)')'Info    : initializing initial water level input'      
         call write_log(logstr, 0)
         !
         iok = check_file_exists(zsinifile, 'Initial conditions ini file', .true.)
         !
         if (is_netcdf_file(zsinifile)) then
            !
            ! Read netcdf (!) initial water level file
            ! Note - newer type real*8 for zs 
            !
            ! NetCDF file
            !
            write(logstr,'(a,a)')'Info    : reading NetCDF initial conditions file ', trim(zsinifile)
            call write_log(logstr, 0)
            !
            !call read_nc_ini_file()
            ! Call the generic quadtree nc file reader function - in real*8
            varname = 'zs'
            call read_netcdf_quadtree_to_sfincs_real8(zsinifile, varname, inizs) !ncfile, varname, varout)                
            !
         else                          
            !
            ! Read binary (!) initial water level file
            ! Note - older type real*4 for zs 
            !
            write(logstr,'(a,a)')'Info    : reading binary initial conditions file ', trim(zsinifile)
            call write_log(logstr, 0)             
            !
            call read_zsini_file()
            !
         endif
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
      ! Called from: set_initial_conditions (sfincs_initial_conditions)
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
         ! Infiltration state is stored in rst_inf1 / rst_inf2 and copied in set_restart_infiltration.
         ! Each record is surrounded by two record markers (read into rdummy).
         !
         if (rsttype==4) then ! Infiltration method cnb
            !
            allocate(rst_inf1(np))
            !
            read(500)rdummy
            read(500)rst_inf1
            !
            rst_inftype = 4
            !
         elseif (rsttype==5) then ! Infiltration method gai
            !
            allocate(rst_inf1(np))
            allocate(rst_inf2(np))
            !
            read(500)rdummy
            read(500)rst_inf1
            read(500)rdummy
            read(500)rdummy
            read(500)rst_inf2
            !
            rst_inftype = 5
            !
         elseif (rsttype==6) then ! Infiltration method horton
            !
            allocate(rst_inf1(np))
            !
            read(500)rdummy
            read(500)rst_inf1
            !
            rst_inftype = 6
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
      ! Called from: set_initial_conditions (sfincs_initial_conditions)
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
   !
   subroutine set_restart_infiltration()
      !
      ! Copy the infiltration state read from the restart file into the infiltration arrays. This must be
      ! called after initialize_infiltration, which allocates and initialises these arrays.
      !
      ! Called from: sfincs_initialize (sfincs_lib)
      !
      implicit none
      !
      if (rst_inftype == 4) then
         !
         if (allocated(scs_Se)) then
            !
            scs_Se = rst_inf1
            call write_log('Info    : scs_Se taken from rstfile', 0)
            !
         else
            !
            call write_log('Warning : rstfile contains cnb infiltration state, but cnb infiltration is not used, state is ignored', 1)
            !
         endif
         !
      elseif (rst_inftype == 5) then
         !
         if (allocated(GA_sigma) .and. allocated(GA_F)) then
            !
            GA_sigma = rst_inf1
            GA_F     = rst_inf2
            write(logstr,'(a,a)')'Info    : GA_sigma and GA_F taken from rstfile, overwrites input values of ', trim(sigmafile)
            call write_log(logstr, 0)
            !
         else
            !
            call write_log('Warning : rstfile contains gai infiltration state, but gai infiltration is not used, state is ignored', 1)
            !
         endif
         !
      elseif (rst_inftype == 6) then
         !
         if (allocated(rain_T1)) then
            !
            rain_T1 = rst_inf1
            write(logstr,'(a,a)')'Info    : rain_T1 taken from rstfile, complements input values of ', trim(fcfile)
            call write_log(logstr, 0)
            !
         else
            !
            call write_log('Warning : rstfile contains hor infiltration state, but hor infiltration is not used, state is ignored', 1)
            !
         endif
         !
      endif
      !
      if (allocated(rst_inf1)) deallocate(rst_inf1)
      if (allocated(rst_inf2)) deallocate(rst_inf2)
      !
      rst_inftype = 0
      !
   end subroutine set_restart_infiltration
   !
end module
