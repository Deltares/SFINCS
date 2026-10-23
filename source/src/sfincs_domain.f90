module sfincs_domain
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! Model domain set-up: builds the computational mesh (regular or quadtree), indices,
   ! masks and neighbours, and initializes bathymetry, boundaries, roughness, storage volume
   ! and the hydrodynamic state.
   !
   ! Subroutines / functions:
   !
   !   initialize_domain()
   !     Initialize the complete model domain by calling the routines below in order.
   !     Called from: sfincs_initialize (sfincs_lib)
   !
   !   initialize_processes()
   !     Set the flags of the physical processes that are active in this run.
   !     Called from: initialize_domain (sfincs_domain)
   !
   !   initialize_mesh()
   !     Build the SFINCS mesh: cell and uv indices, flags and neighbours.
   !     Called from: initialize_domain (sfincs_domain)
   !
   !   initialize_bathymetry()
   !     Read the bed levels (dep file, quadtree file or subgrid tables) and set them at
   !     cells and uv points.
   !     Called from: initialize_domain (sfincs_domain)
   !
   !   update_bed_level()
   !     Apply the external bed level change dzbext to the bed level arrays, lift the water
   !     surface over the moving cells and refresh the bed level at observation points.
   !     Called from: bmi_update_bed_level (sfincs_bmi)
   !
   !   initialize_boundaries()
   !     Count the boundary cells per type and set up the boundary index arrays.
   !     Called from: initialize_domain (sfincs_domain)
   !
   !   initialize_roughness()
   !     Set the Manning roughness at the uv points (uniform, land/sea or from file).
   !     Called from: initialize_domain (sfincs_domain)
   !
   !   initialize_storage_volume()
   !     Read the spatially-varying storage volume per cell (e.g. for green infrastructure).
   !     Called from: initialize_domain (sfincs_domain)
   !
   !   set_advection_mask()
   !     Set the mask that switches off advection at uv points next to boundaries and wave
   !     makers.
   !     Called from: initialize_domain (sfincs_domain)
   !
   !   initialize_hydro()
   !     Allocate and initialize the hydrodynamic variables (zs, q, uv, etc.).
   !     Called from: initialize_domain (sfincs_domain)
   !
   !   deallocate_quadtree()
   !     Deallocate the quadtree arrays that are no longer needed after initialization.
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

   use sfincs_log
   use sfincs_error

   !
   ! Moving bed (update_bed_level): cells with less water than this are not stamped and do
   ! not receive filtered volume; the threshold is max(huthresh, movingbed_hwet_min), as the
   ! non-hydrostatic solver's huthresh_nh.
   !
   real*4, parameter :: movingbed_hwet_min = 0.01

contains

   subroutine initialize_domain()
   !
   ! Initialize the complete model domain by calling the routines below in order.
   !
   ! Called from: sfincs_initialize (sfincs_lib)
   !
   use sfincs_data
   use quadtree
   use sfincs_infiltration   
   use sfincs_vegetation
   use sfincs_timestep_analysis   
   !
   implicit none
   !
   call initialize_processes() 
   !
   call initialize_mesh()
   !
   call initialize_bathymetry()
   !
   call initialize_boundaries()
   !
   call initialize_roughness()
   !
   call initialize_storage_volume()
   !
   call initialize_vegetation()
   !   
   call initialize_hydro()
   !
   if (timestep_analysis) then
      !
      call initialize_timestep_analysis()
      ! 
   endif   
   !
   if (quadtree_nr_levels == 1 .and. .not. use_quadtree_output) then
      !
      ! Only one refinement level found in quadtree file. 
      ! Reverting to use_quadtree is false (default), with use_quadtree_output = false (default).
      ! This means netcdf output will be written to a regular grid.
      !
      use_quadtree = .false.
      !
   endif
   !
   call set_advection_mask()
   !
   call fill_h73_tables()
   !
   end subroutine


   subroutine initialize_processes()
   !
   ! Initialize physical processes
   !
   ! Called from: initialize_domain (sfincs_domain)
   !
   use sfincs_data
   !
   implicit none
   !
   ! Set some flags
   !
   use_quadtree = .false.
   !
   if (qtrfile(1:4) /= 'none') then
      !
      use_quadtree = .true.
      call write_log('Info    : using quadtree mesh', 0)
      !
   else
      !
      call write_log('Info    : using regular mesh', 0)
      !
   endif
   !
   ! Always turn off waves at boundaries. Using wave makers for that sort of thing now.
   !   
   waves = .false.
   !
   wind    = .false.
   precip  = .false.
   patmos  = .false.
   meteo3d = .false.
   !
   boundaries_in_mask = .false.
   water_level_boundaries_in_mask = .false.
   outflow_boundaries_in_mask = .false.
   downstream_river_boundaries_in_mask = .false.
   neumann_boundaries_in_mask = .false.
   !
   if (spwfile(1:4) /= 'none' .or. wndfile(1:4) /= 'none' .or. amufile(1:4) /= 'none' .or. netamuamvfile(1:4) /= 'none' .or. netspwfile(1:4) /= 'none') then
      !
      wind = .true.
      call write_log('Info    : turning on wind', 0)
      !
      if (spw_precip) then
         precip = .true.
         call write_log('Info    : turning on precipitation from spwfile', 0)
      endif
      !
   endif
   !
   if ((ampfile(1:4) /= 'none' .or. netampfile(1:4) /= 'none' .or. spwfile(1:4) /= 'none' .or. netspwfile(1:4) /= 'none') .and. baro==1) then
      patmos = .true.
      call write_log('Info    : turning on atmospheric pressure', 0)
   endif
   !
   if (prcpfile(1:4) /= 'none' .or. amprfile(1:4) /= 'none' .or. netamprfile(1:4) /= 'none' .or. spw_precip) then
       precip = .true.
      call write_log('Info    : turning on precipitation', 0)
   endif
   !
   ! Check for time-spatially varying meteo
   !
   if (spwfile(1:4) /= 'none' .or. amufile(1:4) /= 'none' .or. amprfile(1:4) /= 'none' .or. ampfile(1:4) /= 'none' .or. netampfile(1:4) /= 'none' .or. netamprfile(1:4) /= 'none' .or. netamuamvfile(1:4) /= 'none' .or. netspwfile(1:4) /= 'none') then
      meteo3d = .true.
      call write_log('Info    : using gridded meteo data', 0)
   endif
   !
   wind_reduction = .false.
   if (spwfile(1:4) /= 'none' .or. netspwfile(1:4) /= 'none') then
      if (z0lfile(1:4) /= 'none') then
         wind_reduction = .true.
         call write_log('Info    : turning on wind reduction over land', 0)
      endif
   endif            
   !
   ! If not meteo3d, then also no storage of meteo data in netcdf output
   !
   if (.not. meteo3d) then
      store_meteo = .false.
   endif
   !
   end subroutine


   subroutine initialize_mesh()
   !
   ! Initialize SFINCS domain (indices, flags and neighbors)
   !
   ! Called from: initialize_domain (sfincs_domain)
   !
   use sfincs_data
   use sfincs_nonhydrostatic, only: mask_nonh
   use quadtree
   !
   implicit none
   !
   integer ip, n, m, nm, nmu, num, iref, npu, npv, ndm, nmd, iq
   integer ndmu, nmdu, numu
   integer irefuv, npuvtotal, icuv
   logical ok
   !
   ! Temporary arrays
   !
   integer*1, dimension(:),   allocatable :: z_flags_md
   integer*4, dimension(:),   allocatable :: z_index_z_md1
   integer*4, dimension(:),   allocatable :: z_index_z_md2
   integer*1, dimension(:),   allocatable :: z_flags_mu
   integer*4, dimension(:),   allocatable :: z_index_z_mu1
   integer*4, dimension(:),   allocatable :: z_index_z_mu2
   integer*1, dimension(:),   allocatable :: z_flags_nd
   integer*4, dimension(:),   allocatable :: z_index_z_nd1
   integer*4, dimension(:),   allocatable :: z_index_z_nd2
   integer*1, dimension(:),   allocatable :: z_flags_nu
   integer*4, dimension(:),   allocatable :: z_index_z_nu1
   integer*4, dimension(:),   allocatable :: z_index_z_nu2
   integer*4, dimension(:),   allocatable :: z_index_cuv_md
   integer*4, dimension(:),   allocatable :: z_index_cuv_mu
   integer*4, dimension(:),   allocatable :: z_index_cuv_nd
   integer*4, dimension(:),   allocatable :: z_index_cuv_nu
   !
   integer*1, dimension(:),   allocatable :: msk
   !   
   ! For 'old' input
   !
   integer*4,          dimension(:),   allocatable :: indices
   !
   real*4 :: dxymin
   !
   logical :: okay
   !
   ! READ MESH
   !
   ! 2 options:
   !
   ! 1) quadtree grid (netcdf file)
   ! 2) 'old' inputs with index file
   !
   if (use_quadtree) then
      ! 
      ! Read quadtree file
      !
      call quadtree_read_file(qtrfile, snapwave, nonhydrostatic)
      !
   else
      !
      if (inputtype=='asc') then
         !
         ! Read ASCII input
         !
         call write_log('Info : reading ASCII input ...', 0)
         !
         allocate(kcsg(nmax,  mmax))       ! KCS
         !
         kcsg = 0      
         !
         ! Read mask file
         !
         ok = check_file_exists(mskfile, 'Mask file', .true.)
         !
         open(500, file=trim(mskfile))
         do n = 1, nmax
            read(500,*)(kcsg(n, m), m = 1, mmax)
         enddo
         close(500)
         !
         ! Count active points
         !
         np = 0
         do m = 1, mmax
            do n = 1, nmax
               if (kcsg(n,m)>0) then
                  np = np + 1
               endif
            enddo
         enddo
         !
         allocate(indices(np))
         !
         ip = 0
         do m = 1, mmax
            do n = 1, nmax
               if (kcsg(n,m)>0) then
                  ip = ip + 1
                  indices(ip) = (m - 1)*nmax + n
               endif
            enddo
         enddo
         !
      else
         !
         ! Read 'old' binary index file and make quadtree with those
         !
         ! Read index file (first time, only read number of active points)
         !
         write(logstr,'(a,a)')'Info    : reading index file ', trim(indexfile)
         call write_log(logstr, 0)
         !
         ok = check_file_exists(indexfile, 'Index file', .true.)
         !
         open(unit = 500, file = trim(indexfile), form = 'unformatted', access = 'stream')
         read(500)np
         allocate(indices(np))
         read(500)indices
         close(500)
         !
      endif
      !
      call make_quadtree_from_indices(np, indices, nmax, mmax, x0, y0, dx, dy, rotation)
      !
   endif   
   !   
   rotation = quadtree_rotation
   x0       = quadtree_x0
   y0       = quadtree_y0
   nmax     = quadtree_nmax
   mmax     = quadtree_mmax
   dx       = quadtree_dx
   dy       = quadtree_dy
   cosrot   = cos(rotation)
   sinrot   = sin(rotation)
   !
   nref     = quadtree_nr_levels
   !
   ! Allocate some temporary arrays
   !
   allocate(msk(quadtree_nr_points))
   !
   ! Allocate some permanent arrays
   !
   allocate(index_sfincs_in_quadtree(quadtree_nr_points))
   index_sfincs_in_quadtree = 0     
   !
   if (inputtype=='asc') then
      !
      ip = 0
      do m = 1, mmax
         do n = 1, nmax
            if (kcsg(n, m) > 0) then
               ip = ip + 1
               msk(ip) = kcsg(n, m)
            endif
         enddo
      enddo
      !
   else
      !
      if (use_quadtree .and. quadtree_netcdf) then
         !
         ! Mask is stored in netcdf file 
         !
         msk = quadtree_mask
         !
      else
         !
         ! Mask in separate file
         !
         ! Read mask file (must have same number of cells as quadtree file !)
         !
         write(logstr,'(a,a)')'Info    : reading mask file : ', trim(mskfile)
         call write_log(logstr, 0)
         !
         ok = check_file_exists(mskfile, 'Mask file', .true.)
         !
         open(unit = 500, file = trim(mskfile), form = 'unformatted', access = 'stream')
         read(500)msk
         close(500)
         !
      endif
      !
   endif
   !
   ! Now first get rid of all the points that have msk==0
   !
   ! First count them
   !
   np = 0
   !
   do ip = 1, quadtree_nr_points
      if (msk(ip) > 0) then
         np = np + 1
      endif
   enddo
   !
   ! Found the number of active points 
   !
   ! Now allocate a whole bunch of arrays
   !
   ! Permanent
   !
   allocate(kcs(np))
   !
   allocate(index_quadtree_in_sfincs(np))
   index_sfincs_in_quadtree = 0
   !
   allocate(z_flags_iref(np))
   !
   z_flags_iref = 0
   !      
   ! Temporary index arrays
   !
   allocate(z_flags_md(np))
   allocate(z_flags_mu(np))
   allocate(z_flags_nd(np))
   allocate(z_flags_nu(np))
   allocate(z_index_z_md1(np))
   allocate(z_index_z_md2(np))
   allocate(z_index_z_mu1(np))
   allocate(z_index_z_mu2(np))
   allocate(z_index_z_nd1(np))
   allocate(z_index_z_nd2(np))
   allocate(z_index_z_nu1(np))
   allocate(z_index_z_nu2(np))
   !
   allocate(z_index_z_n(np))
   allocate(z_index_z_m(np))
   !
   z_flags_md    = 0
   z_flags_mu    = 0
   z_flags_nd    = 0
   z_flags_nu    = 0
   z_index_z_md1 = 0
   z_index_z_md2 = 0
   z_index_z_mu1 = 0
   z_index_z_mu2 = 0
   z_index_z_nd1 = 0
   z_index_z_nd2 = 0
   z_index_z_nu1 = 0
   z_index_z_nu2 = 0
   z_index_z_n   = 0
   z_index_z_m   = 0
   !
   ! Permanent index arrays
   !      
   allocate(z_index_uv_md1(np))
   allocate(z_index_uv_md2(np))
   allocate(z_index_uv_mu1(np))
   allocate(z_index_uv_mu2(np))
   allocate(z_index_uv_nd1(np))
   allocate(z_index_uv_nd2(np))
   allocate(z_index_uv_nu1(np))
   allocate(z_index_uv_nu2(np))
   !
   allocate(z_index_uv_md(np))
   allocate(z_index_uv_mu(np))
   allocate(z_index_uv_nd(np))
   allocate(z_index_uv_nu(np))
   !
   allocate(z_index_cuv_md(np))
   allocate(z_index_cuv_mu(np))
   allocate(z_index_cuv_nd(np))
   allocate(z_index_cuv_nu(np))
   !
   z_index_uv_md1 = 0
   z_index_uv_md2 = 0
   z_index_uv_mu1 = 0
   z_index_uv_mu2 = 0
   z_index_uv_nd1 = 0
   z_index_uv_nd2 = 0
   z_index_uv_nu1 = 0
   z_index_uv_nu2 = 0
   z_index_uv_md  = 0
   z_index_uv_mu  = 0
   z_index_uv_nd  = 0
   z_index_uv_nu  = 0
   z_index_z_n    = 0
   z_index_z_m    = 0
   !
   ! Re-mapping indices
   !
   np = 0
   !
   do ip = 1, quadtree_nr_points
      if (msk(ip)>0) then
         !
         np = np + 1
         index_sfincs_in_quadtree(ip) = np
         index_quadtree_in_sfincs(np) = ip
         !
      endif
   enddo
   !
   ! Copy and re-map indices in temporary arrays
   !
   do ip = 1, np
      !
      iq = index_quadtree_in_sfincs(ip)
      !
      kcs(ip)           = msk(iq)
      !
      z_flags_iref(ip)  = quadtree_level(iq)
      !
      z_index_z_n(ip)   = quadtree_n(iq)
      z_index_z_m(ip)   = quadtree_m(iq)
      !
      z_flags_nu(ip)    = quadtree_nu(iq)
      z_flags_nd(ip)    = quadtree_nd(iq)
      z_flags_mu(ip)    = quadtree_mu(iq)
      z_flags_md(ip)    = quadtree_md(iq)
      !
      ! Neighbors
      !
      if (quadtree_nu1(iq)>0) then
         z_index_z_nu1(ip)   = index_sfincs_in_quadtree(quadtree_nu1(iq))
      endif
      if (quadtree_nu2(iq)>0) then
         z_index_z_nu2(ip)   = index_sfincs_in_quadtree(quadtree_nu2(iq))
      endif
      if (quadtree_nd1(iq)>0) then
         z_index_z_nd1(ip)   = index_sfincs_in_quadtree(quadtree_nd1(iq))
      endif
      if (quadtree_nd2(iq)>0) then
         z_index_z_nd2(ip)   = index_sfincs_in_quadtree(quadtree_nd2(iq))
      endif
      if (quadtree_mu1(iq)>0) then
         z_index_z_mu1(ip)   = index_sfincs_in_quadtree(quadtree_mu1(iq))
      endif
      if (quadtree_mu2(iq)>0) then
         z_index_z_mu2(ip)   = index_sfincs_in_quadtree(quadtree_mu2(iq))
      endif
      if (quadtree_md1(iq)>0) then
         z_index_z_md1(ip)   = index_sfincs_in_quadtree(quadtree_md1(iq))
      endif
      if (quadtree_md2(iq)>0) then
         z_index_z_md2(ip)   = index_sfincs_in_quadtree(quadtree_md2(iq))
      endif
      !
   enddo   
   ! 
   ! Loop through quadtree to count the number of u and v points
   !
   npu = 0
   npv = 0
   !
   do nm = 1, np
      !
      ! Right
      !
      if (z_flags_mu(nm)<1) then
         !
         if (z_index_z_mu1(nm)>0) then
            npu = npu + 1
         endif   
         !
      else
         !
         if (z_index_z_mu1(nm)>0) then
            npu = npu + 1
         endif   
         !
         if (z_index_z_mu2(nm)>0) then
            npu = npu + 1
         endif   
         !            
      endif   
      !
      ! Above
      !
      if (z_flags_nu(nm)<1) then
         !
         if (z_index_z_nu1(nm)>0) then
            npv = npv + 1
         endif   
         !
      else
         !
         if (z_index_z_nu1(nm)>0) then
            npv = npv + 1
         endif   
         !
         if (z_index_z_nu2(nm)>0) then
            npv = npv + 1
         endif   
         !            
      endif   
      !
   enddo 
   !
   npuv = npu + npv
   !
   ncuv = 0 
   !
   if (quadtree_nr_levels>1) then 
      !
      ! Loop through all the points again and find 'combined' uv points at refinement boundaries
      !
      do nm = 1, np
         !
         ! Left
         !
         if (z_flags_md(nm) == 1) then
            !
            ! Finer on the left
            !
            if (z_index_z_md1(nm)>0 .or. z_index_z_md2(nm)>0) then
               ncuv = ncuv + 1
            endif 
            !            
         endif   
         !
         ! Right
         !
         if (z_flags_mu(nm) == 1) then
            !
            ! Finer on the right
            !
            if (z_index_z_mu1(nm)>0 .or. z_index_z_mu2(nm)>0) then
               ncuv = ncuv + 1
            endif
            !            
         endif   
         !
         ! Below
         !
         if (z_flags_nd(nm) == 1) then
            !
            ! Finer below
            !
            if (z_index_z_nd1(nm)>0 .or. z_index_z_nd2(nm)>0) then
               ncuv = ncuv + 1
            endif
            !            
         endif   
         !
         ! Above
         !
         if (z_flags_nu(nm) == 1) then
            !
            ! Finer above
            !
            if (z_index_z_nu1(nm)>0 .or. z_index_z_nu2(nm)>0) then
               ncuv = ncuv + 1
            endif
            !            
         endif   
         !
      enddo 
      ! 
      write(logstr,'(a,i0)')'Info    : number of quadtree refinement transitions : ', ncuv
      call write_log(logstr, 0)
      !
   endif
   !
   npuvtotal = npuv + ncuv ! total number of points in uv arrays (including combined uv points)
   !
   ! Set neighbors of cells to dummy uv points 
   !
   z_index_uv_md  = npuv + ncuv + 1
   z_index_uv_mu  = npuv + ncuv + 1
   z_index_uv_nd  = npuv + ncuv + 1
   z_index_uv_nu  = npuv + ncuv + 1
   ! 
   ! Now allocate and fill arrays with combined uv points
   !
   allocate(cuv_index_uv(ncuv))  ! index in uv array of combined point
   allocate(cuv_index_uv1(ncuv)) ! index in uv array of first uv point
   allocate(cuv_index_uv2(ncuv)) ! index in uv array of second uv point
   !
   ! Initialize cuv_index_uv1 and cuv_index_uv2 at dummy uv index (npuv + ncuv + 1)
   !
   cuv_index_uv1 = npuv + ncuv + 1
   cuv_index_uv2 = npuv + ncuv + 1
   !
   irefuv = 0 
   !
   do nm = 1, np
      !
      ! Left
      !
      if (z_flags_md(nm) == 1) then
         !
         ! Finer on the left
         !
         if (z_index_z_md1(nm)>0 .or. z_index_z_md2(nm)>0) then
            irefuv = irefuv + 1
            cuv_index_uv(irefuv)  = npuv + irefuv ! used to update uv and q in momentum
            z_index_cuv_md(nm)    = irefuv ! temporary
         endif
         !            
      endif   
      !
      ! Right
      !
      if (z_flags_mu(nm) == 1) then
         !
         ! Finer on the right
         !
         if (z_index_z_mu1(nm)>0 .or. z_index_z_mu2(nm)>0) then
            irefuv = irefuv + 1
            cuv_index_uv(irefuv)  = npuv + irefuv
            z_index_cuv_mu(nm)    = irefuv ! temporary
         endif
         !            
      endif   
      !
      ! Below
      !
      if (z_flags_nd(nm) == 1) then
         !
         ! Finer below
         !
         if (z_index_z_nd1(nm)>0 .or. z_index_z_nd2(nm)>0) then
            irefuv = irefuv + 1
            cuv_index_uv(irefuv)  = npuv + irefuv ! needed in momentum
            z_index_cuv_nd(nm)    = irefuv ! temporary         
         endif
         !            
      endif   
      !
      ! Above
      !
      if (z_flags_nu(nm) == 1) then
         !
         ! Finer above
         !
         if (z_index_z_nu1(nm)>0 .or. z_index_z_nu2(nm)>0) then
            irefuv = irefuv + 1  
            cuv_index_uv(irefuv)  = npuv + irefuv
            z_index_cuv_nu(nm)    = irefuv
         endif
         !            
      endif   
      !
   enddo 
   !
   ! UV-points
   !
   allocate(kcuv(npuv + ncuv + 1))
   allocate(uv_index_z_nm(npuv))
   allocate(uv_index_z_nmu(npuv))
   allocate(uv_flags_iref(npuv))
   allocate(uv_flags_dir(npuv))
   allocate(uv_flags_type(npuv))
   !
   uv_flags_iref = 0
   uv_flags_dir  = 0
   uv_flags_type = 0
   !
   kcuv = 0
   uv_index_z_nm  = 0
   uv_index_z_nmu = 0
   !
   ! If advection, viscosity, coriolis or thetasmoothing, allocate extra arrays for 8 uv neighbor indices
   !
   if (advection .or. coriolis .or. viscosity .or. thetasmoothing .or. friction2d) then
      !
      allocate(uv_index_u_nmd(npuv))
      allocate(uv_index_u_nmu(npuv))
      allocate(uv_index_u_num(npuv))
      allocate(uv_index_u_ndm(npuv))
      allocate(uv_index_v_ndm(npuv))
      allocate(uv_index_v_nm(npuv))
      allocate(uv_index_v_ndmu(npuv))
      allocate(uv_index_v_nmu(npuv))
      !
      uv_index_u_nmd  = 0
      uv_index_u_nmu  = 0
      uv_index_u_num  = 0
      uv_index_u_ndm  = 0
      uv_index_v_ndm  = 0
      uv_index_v_nm   = 0
      uv_index_v_ndmu = 0
      uv_index_v_nmu  = 0
      !
   endif
   !
   ip = 0 ! Number of uv points
   !
   ! Loop through cells to set indices neighboring cells for each uv point
   ! Also fill temporary arrays with indices of all 8 uv points for each cell 
   ! Also set some flags for each uv point
   !
   do nm = 1, np
      !
      ! Right
      !
      if (z_flags_mu(nm)==0) then
         !
         ! Same level
         !
         if (z_index_z_mu1(nm)>0) then ! and there is a neighbor
            !
            ip = ip + 1
            !
            nmu = z_index_z_mu1(nm)
            !
            uv_index_z_nm(ip)     = nm
            uv_index_z_nmu(ip)    = nmu
            z_index_uv_mu(nm)     = ip 
            z_index_uv_md(nmu)    = ip
            z_index_uv_mu1(nm)     = ip 
            z_index_uv_md1(nmu)    = ip
            !
            ! Flags to describe uv point 
            !
            uv_flags_iref(ip)                  = z_flags_iref(nm) 
            uv_flags_dir(ip)                   = 0 ! u point
            uv_flags_type(ip)                  = 0 ! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            !
         endif   
         !
      elseif (z_flags_mu(nm)==-1) then
         !
         ! Coarser to the right
         !
         if (z_index_z_mu1(nm)>0) then ! and there is a neighbor
            !
            ip  = ip + 1
            nmu = z_index_z_mu1(nm)
            ! 
            uv_index_z_nm(ip)     = nm
            uv_index_z_nmu(ip)    = nmu
            z_index_uv_mu(nm)     = ip
            z_index_uv_mu1(nm)    = ip
            !
            ! Set the combined uv point of the neighbor to the right (icuv is index of combined uv point in cuv array)
            !
            icuv = z_index_cuv_md(nmu)
            !
            ! set indices
            !
            if (z_index_z_md1(z_index_z_mu1(nm)) == nm) then
               !
               ! Lower cell
               !
               cuv_index_uv1(icuv) = ip
               z_index_uv_md1(nmu) = ip
               !
            else
               !
               ! Upper cell
               !
               cuv_index_uv2(icuv) = ip              
               z_index_uv_md2(nmu) = ip
               !
            endif
            !
            ! Combined uv point for cell on the right
            !  
            z_index_uv_md(nmu) = cuv_index_uv(icuv)
            !
            ! Flags to describe uv point 
            !
            uv_flags_iref(ip)                  = z_flags_iref(nm) 
            uv_flags_dir(ip)                   = 0  ! u point
            uv_flags_type(ip)                  = -1 ! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            !
         endif
         !
      else
         !
         ! Finer to the right
         !
         if (z_index_z_mu1(nm)>0) then
            !
            ! Lower cell
            !
            ip                   = ip + 1
            nmu                  = z_index_z_mu1(nm) ! z index of lower fine cell on the right
            uv_index_z_nm(ip)    = nm
            uv_index_z_nmu(ip)   = nmu
            z_index_uv_mu1(nm)   = ip
            z_index_uv_md(nmu)   = ip
            z_index_uv_md1(nmu)  = ip
            !
            icuv = z_index_cuv_mu(nm) ! (icuv is index of combined uv point in cuv array)
            !  
            ! Set the combined uv point of this cell
            !
            z_index_uv_mu(nm) = cuv_index_uv(icuv)
            !
            cuv_index_uv1(icuv) = ip              
            !
            ! Flags to describe uv point 
            !
            uv_flags_iref(ip)                  = z_flags_iref(nm) + 1
            uv_flags_dir(ip)                   = 0 ! u point
            uv_flags_type(ip)                  = 1! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            !
         endif   
         !
         if (z_index_z_mu2(nm)>0) then
            !
            ! Upper cell
            !
            ip                   = ip + 1
            nmu                  = z_index_z_mu2(nm) ! z index of upper fine cell on the right
            uv_index_z_nm(ip)    = nm
            uv_index_z_nmu(ip)   = nmu
            z_index_uv_mu2(nm)   = ip
            z_index_uv_md(nmu)   = ip
            z_index_uv_md1(nmu)  = ip
            !
            icuv = z_index_cuv_mu(nm) ! (icuv is index of combined uv point in cuv array)
            !  
            ! Set the combined uv point of this cell (this may already have been done in the lower cell)
            !
            z_index_uv_mu(nm) = cuv_index_uv(icuv)
            !
            cuv_index_uv2(icuv) = ip              
            !
            ! Flags to describe uv point 
            !
            uv_flags_iref(ip)                  = z_flags_iref(nm) + 1
            uv_flags_dir(ip)                   = 0 ! u point
            uv_flags_type(ip)                  = 1! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            !
         endif   
         !            
      endif   
      !
      ! Above
      !
      if (z_flags_nu(nm)==0) then
         !
         ! Same level
         !
         if (z_index_z_nu1(nm)>0) then ! and there is a neighbor
            !
            ip = ip + 1
            !
            num = z_index_z_nu1(nm)
            !
            uv_index_z_nm(ip)     = nm
            uv_index_z_nmu(ip)    = num
            z_index_uv_nu(nm)     = ip 
            z_index_uv_nd(num)    = ip
            z_index_uv_nu1(nm)     = ip 
            z_index_uv_nd1(num)    = ip
            !
            ! Flags to describe uv point 
            !
            uv_flags_iref(ip)                  = z_flags_iref(nm) 
            uv_flags_dir(ip)                   = 1 ! u point
            uv_flags_type(ip)                  = 0 ! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            !
         endif   
         !
      elseif (z_flags_nu(nm)==-1) then
         !
         ! Coarser above
         !
         if (z_index_z_nu1(nm)>0) then ! and there is a neighbor
            !
            ip  = ip + 1
            num = z_index_z_nu1(nm)
            ! 
            uv_index_z_nm(ip)     = nm
            uv_index_z_nmu(ip)    = num
            z_index_uv_nu(nm)     = ip
            z_index_uv_nu1(nm)    = ip
            !
            ! Set the combined uv point of the neighbor above (icuv is index of combined uv point in cuv array)
            !
            icuv = z_index_cuv_nd(num)
            !
            ! set indices
            !
            if (z_index_z_nd1(z_index_z_nu1(nm)) == nm) then
               !
               ! Lower cell
               !
               cuv_index_uv1(icuv) = ip
               z_index_uv_nd1(num) = ip
               !
            else
               !
               ! Upper cell
               !
               cuv_index_uv2(icuv) = ip              
               z_index_uv_nd2(num) = ip
               !
            endif
            !
            ! Combined uv point for cell above
            !  
            z_index_uv_nd(num) = cuv_index_uv(icuv)
            !
            ! Flags to describe uv point 
            !
            uv_flags_iref(ip)                  = z_flags_iref(nm) 
            uv_flags_dir(ip)                   = 1 ! v point
            uv_flags_type(ip)                  = -1! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            !
         endif  
         !
      else
         !
         ! Finer above
         !
         if (z_index_z_nu1(nm)>0) then
            !
            ! Left cell
            !
            ip                   = ip + 1
            num                  = z_index_z_nu1(nm) ! z index of left fine cell above
            uv_index_z_nm(ip)    = nm
            uv_index_z_nmu(ip)   = num
            z_index_uv_nd(num)   = ip
            z_index_uv_nu1(nm)   = ip
            z_index_uv_nd1(num)  = ip
            !
            icuv = z_index_cuv_nu(nm) ! (icuv is index of combined uv point in cuv array)
            !  
            ! Set the combined uv point of this cell
            !
            z_index_uv_nu(nm) = cuv_index_uv(icuv)
            !
            cuv_index_uv1(icuv) = ip              
            !
            ! Flags to describe uv point 
            !
            uv_flags_iref(ip)                  = z_flags_iref(nm) + 1
            uv_flags_dir(ip)                   = 1 ! v point
            uv_flags_type(ip)                  = 1! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            !
         endif   
         !
         if (z_index_z_nu2(nm)>0) then
            !
            ! Right cell
            !
            ip                   = ip + 1
            num                  = z_index_z_nu2(nm) ! z index of left fine cell above
            uv_index_z_nm(ip)    = nm
            uv_index_z_nmu(ip)   = num
            z_index_uv_nd(num)   = ip
            z_index_uv_nu2(nm)   = ip
            z_index_uv_nd1(num)  = ip
            !
            icuv = z_index_cuv_nu(nm) ! (icuv is index of combined uv point in cuv array)
            !  
            ! Set the combined uv point of this cell
            !
            z_index_uv_nu(nm) = cuv_index_uv(icuv)
            !
            cuv_index_uv2(icuv) = ip              
            !
            ! Flags to describe uv point 
            !
            uv_flags_iref(ip)                  = z_flags_iref(nm) + 1
            uv_flags_dir(ip)                   = 1 ! v point
            uv_flags_type(ip)                  = 1 ! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            !
         endif   
         !            
      endif   
      !
   enddo
   !
   ! And now set the indices of the 8 uv neighbors of the uv points (u_nmu, u_nmd, u_num,  u_ndm,
   !                                                                 v_ndm, v_nm,  v_ndmu, v_nmu)
   if (advection .or. coriolis .or. viscosity .or. friction2d .or. thetasmoothing) then
      !
      do ip = 1, npuv
         !
         if (uv_flags_dir(ip) == 0) then
            !
            ! U point
            !
            nm  = uv_index_z_nm(ip)  ! index of cell to the left of this uv point
            nmu = uv_index_z_nmu(ip) ! index of cell to the right of this uv point
            !
            ! Left (u_nmd)
            ! 
            uv_index_u_nmd(ip) = z_index_uv_md(nm)
            !
            ! Right (u_nmu)
            ! 
            uv_index_u_nmu(ip) = z_index_uv_mu(nmu)
            !
            ! Below (u_ndm)
            !
            if (z_flags_nd(nm)>-1) then ! The cell below is coarser or same level. Easy.
               !    
               ndm = z_index_z_nd1(nm) ! z-index of cell below
               !
               if (ndm>0) then ! And it exists
                  !
                  uv_index_u_ndm(ip) = z_index_uv_mu(ndm)
                  !
               endif 
               !
            else
               ! 
               ! Cell below is finer. Difficult.
               ! 
               ! Check if cell ndmu is same level. If so, use z_index_uv_md of that cell. Otherwise, use z_index_uv_mu of cell below. 
               !    
               ndm = z_index_z_nd2(nm) ! z-index of cell below
               !
               ! Check if cell below exists
               ! 
               if (ndm>0) then
                  !
                  if (z_flags_mu(ndm)==0) then ! cell ndmu is also finer level (could also be same level, but not coarser)
                     !
                     uv_index_u_ndm(ip) = z_index_uv_mu(ndm)
                     !
                  else ! cell ndmu is same level
                     !
                     ndmu = z_index_z_mu1(ndm) ! z-index of cell below-right
                     !
                     if (ndmu>0) then
                        !
                        uv_index_u_ndm(ip) = z_index_uv_md(ndmu)
                        !
                     endif  
                     !
                  endif  
               endif    
            endif   
            !
            ! Above (u_num)
            !
            if (z_flags_nu(nm)>-1) then ! The cell above is coarser or same level. Easy.
               !    
               num = z_index_z_nu1(nm) ! z-index of cell above
               !
               if (num>0) then ! And it exists
                  !
                  uv_index_u_num(ip) = z_index_uv_mu(num)
                  !
               endif 
               !
            else
               ! 
               ! Cell above is finer. Difficult.
               ! 
               ! Check if cell numu is same level. If so, use z_index_uv_md of that cell. Otherwise, use z_index_uv_mu of cell above. 
               !    
               num = z_index_z_nu2(nm) ! z-index of cell above
               !
               ! Check if cell above exists
               ! 
               if (num>0) then
                  !
                  if (z_flags_mu(num)==0) then ! cell numu is also finer level (could also same level, but not coarser)
                     !
                     uv_index_u_num(ip) = z_index_uv_mu(num)
                     !
                  else ! cell numu is same level
                     !
                     numu = z_index_z_mu1(num) ! z-index of cell below-right
                     !
                     if (numu>0) then
                        !
                        uv_index_u_num(ip) = z_index_uv_md(numu)
                        !
                     endif  
                     !
                  endif  
               endif    
            endif   
            !
            ! Left (v_nm)
            !
            uv_index_v_nm(ip) = z_index_uv_nu(nm)
            !
            ! Left-below (v_ndm)
            !
            uv_index_v_ndm(ip) = z_index_uv_nd(nm)
            !
            ! Right (v_nmu)
            !
            uv_index_v_nmu(ip) = z_index_uv_nu(nmu)
            !
            ! Right-below (v_ndmu)
            !
            uv_index_v_ndmu(ip) = z_index_uv_nd(nmu)
            !
         else
            !
            ! V point (mirror-image of U point)
            !
            nm  = uv_index_z_nm(ip)  ! index of cell below this uv point
            num = uv_index_z_nmu(ip) ! index of cell above this uv point
            !
            ! Left (u_nmd)
            ! 
            uv_index_u_nmd(ip) = z_index_uv_nd(nm)
            !
            ! Right (u_nmu)
            ! 
            uv_index_u_nmu(ip) = z_index_uv_nu(num)
            !
            ! Below (u_ndm)
            !
            if (z_flags_md(nm)>-1) then ! The cell left is coarser or same level. Easy.
               !    
               nmd = z_index_z_md1(nm) ! z-index of cell left
               !
               if (nmd>0) then ! And it exists
                  !
                  uv_index_u_ndm(ip) = z_index_uv_nu(nmd)
                  !
               endif 
               !
            else
               ! 
               ! Cell below is finer. Difficult.
               ! 
               ! Check if cell ndmu is same level. If so, use z_index_uv_md of that cell. Otherwise, use z_index_uv_mu of cell below. 
               !    
               nmd = z_index_z_md2(nm) ! z-index of cell left
               !
               ! Check if cell below exists
               ! 
               if (nmd>0) then
                  !
                  if (z_flags_nu(nmd)==0) then ! cell ndmu is also finer level (could also be same level, but not coarser)
                     !
                     uv_index_u_ndm(ip) = z_index_uv_nu(nmd)
                     !
                  else ! cell nmdu is same level
                     !
                     nmdu = z_index_z_nu1(nmd) ! z-index of cell below-right
                     !
                     if (nmdu>0) then
                        !
                        uv_index_u_ndm(ip) = z_index_uv_nd(nmdu)
                        !
                     endif  
                     !
                  endif  
               endif    
            endif   
            !
            ! Above (u_num)
            !
            if (z_flags_mu(nm)>-1) then ! The cell right coarser or same level. Easy.
               !    
               nmu = z_index_z_mu1(nm) ! z-index of cell right
               !
               if (nmu>0) then ! And it exists
                  !
                  uv_index_u_num(ip) = z_index_uv_nu(nmu)
                  !
               endif 
               !
            else
               ! 
               ! Cell above is finer. Difficult.
               ! 
               ! Check if cell numu is same level. If so, use z_index_uv_md of that cell. Otherwise, use z_index_uv_mu of cell above. 
               !    
               nmu = z_index_z_mu2(nm) ! z-index of cell right
               !
               ! Check if cell above exists
               ! 
               if (nmu>0) then
                  !
                  if (z_flags_nu(nmu)==0) then ! cell numu is also finer level (could also same level, but not coarser)
                     !
                     uv_index_u_num(ip) = z_index_uv_nu(nmu)
                     !
                  else ! cell numu is same level
                     !
                     numu = z_index_z_nu1(nmu) ! z-index of cell below-right
                     !
                     if (numu>0) then
                        !
                        uv_index_u_num(ip) = z_index_uv_nd(numu)
                        !
                     endif  
                     !
                  endif  
               endif    
            endif   
            !
            ! Left (v_nm)
            !
            uv_index_v_nm(ip) = z_index_uv_mu(nm)
            !
            ! Left-below (v_ndm)
            !
            uv_index_v_ndm(ip) = z_index_uv_md(nm)
            !
            ! Right (v_nmu)
            !
            uv_index_v_nmu(ip) = z_index_uv_mu(num)
            !
            ! Right-below (v_ndmu)
            !
            uv_index_v_ndmu(ip) = z_index_uv_md(num)
            !
         endif
         !
         ! Missing uv neighbors
         !
         ! For u dir points, set neighbor to this uv point itself
         !    
         if (uv_index_u_nmd(ip) == 0) uv_index_u_nmd(ip) = ip
         if (uv_index_u_nmu(ip) == 0) uv_index_u_nmu(ip) = ip
         if (uv_index_u_ndm(ip) == 0) uv_index_u_ndm(ip) = ip
         if (uv_index_u_num(ip) == 0) uv_index_u_num(ip) = ip
         !
         ! For v dir points, set neighbor to dummy index
         !    
         if (uv_index_v_ndm(ip) == 0) uv_index_v_ndm(ip) = npuvtotal + 1
         if (uv_index_v_ndmu(ip) == 0) uv_index_v_ndmu(ip) = npuvtotal + 1
         if (uv_index_v_nm(ip) == 0) uv_index_v_nm(ip) = npuvtotal + 1
         if (uv_index_v_nmu(ip) == 0) uv_index_v_nmu(ip) = npuvtotal + 1
         ! 
      enddo ! uv point
      ! 
   endif ! advection, viscosity, coriolis or thetasmoothing
   !
   if (nonhydrostatic) then
      !
      ! Intialize mask for cells with non-hydrostatic corrections
      !
      allocate(mask_nonh(np))
      !
      mask_nonh = 0
      !
      do nm = 1, np
         !
         okay = .true.
         !
         if (kcs(nm) /= 1) cycle ! not a regular point
         !
         ! Check that the point has regular (kcs = 1) neighbors. In the general
         ! 2D case all four (left/right/bottom/top) are required. For a true 1-D
         ! model (nmax == 1) cells legitimately have no bottom/top neighbor, so
         ! only left/right are required there; the non-hydrostatic solver then
         ! degrades to a 1-D (x-only) pressure Poisson problem, which all the
         ! downstream assembly already supports through its per-direction guards.
         !
         ! Left
         !
         if (z_flags_md(nm) /= 0) cycle ! not a regular neighbor
         !
         nmd = z_index_z_md1(nm)
         !
         if (nmd == 0) cycle ! no neighbor
         !
         if (kcs(nmd) < 1 .or. kcs(nmd) > 3) cycle ! allow regular (1) or open-boundary (2/3) neighbour
         !
         ! Right
         !
         if (z_flags_mu(nm) /= 0) cycle ! not a regular neighbor
         !
         nmu = z_index_z_mu1(nm)
         !
         if (nmu == 0) cycle ! no neighbor
         !
         if (kcs(nmu) < 1 .or. kcs(nmu) > 3) cycle ! allow regular (1) or open-boundary (2/3) neighbour
         !
         if (nmax > 1) then
            !
            ! Bottom
            !
            if (z_flags_nd(nm) /= 0) cycle ! not a regular neighbor
            !
            ndm = z_index_z_nd1(nm)
            !
            if (ndm == 0) cycle ! no neighbor
            !
            if (kcs(ndm) < 1 .or. kcs(ndm) > 3) cycle ! allow regular (1) or open-boundary (2/3) neighbour
            !
            ! Top
            !
            if (z_flags_nu(nm) /= 0) cycle ! not a regular neighbor
            !
            num = z_index_z_nu1(nm)
            !
            if (num == 0) cycle ! no neighbor
            !
            if (kcs(num) < 1 .or. kcs(num) > 3) cycle ! allow regular (1) or open-boundary (2/3) neighbour
            !
         endif
         !
         ! This cell has the required regular neighbors, so copy from quadtree mask
         !
         mask_nonh(nm) = quadtree_nonh_mask(index_quadtree_in_sfincs(nm))
         !
      enddo
      !
   endif   
   !
   ! Okay, got all the quadtree cells including indices, neighbors flags 
   !
   call write_log('------------------------------------------', 1)
   call write_log('Computational mesh', 1)
   call write_log('------------------------------------------', 1)
   write(logstr,'(a,i9)')'Number of active z points      : ', np
   call write_log(logstr, 1)   
   write(logstr,'(a,i9)')'Number of active u/v points    : ', npuv
   call write_log(logstr, 1)
   !
   if (use_quadtree) then
      do iref = 1, nref
         write(logstr,'(a,i3,a,i9)')'Number of cells in level ', iref, '   : ', quadtree_last_point_per_level(iref) - quadtree_first_point_per_level(iref) + 1
         call write_log(logstr, 1)
      enddo
   endif
   !
   ! Time to de-allocate some arrays
   ! 
   if (allocated(z_flags_md)) deallocate(z_flags_md)
   if (allocated(z_index_z_md1)) deallocate(z_index_z_md1)
   if (allocated(z_index_z_md2)) deallocate(z_index_z_md2)
   if (allocated(z_flags_mu)) deallocate(z_flags_mu)
   if (allocated(z_index_z_mu1)) deallocate(z_index_z_mu1)
   if (allocated(z_index_z_mu2)) deallocate(z_index_z_mu2)
   if (allocated(z_flags_nd)) deallocate(z_flags_nd)
   if (allocated(z_index_z_nd1)) deallocate(z_index_z_nd1)
   if (allocated(z_index_z_nd2)) deallocate(z_index_z_nd2)
   if (allocated(z_flags_nu)) deallocate(z_flags_nu)
   if (allocated(z_index_z_nu1)) deallocate(z_index_z_nu1)
   if (allocated(z_index_z_nu2)) deallocate(z_index_z_nu2)
   if (allocated(z_index_cuv_md)) deallocate(z_index_cuv_md)
   if (allocated(z_index_cuv_mu)) deallocate(z_index_cuv_mu)
   if (allocated(z_index_cuv_nd)) deallocate(z_index_cuv_nd)
   if (allocated(z_index_cuv_nu)) deallocate(z_index_cuv_nu)
   if (allocated(msk)) deallocate(msk)
   if (allocated(indices)) deallocate(indices)
   !
   ! MESH GEOMETRY
   !
   ! Refinement levels
   !
   allocate(dxr(nref))
   allocate(dyr(nref))
   allocate(dxyr(nref))
   !
   ! Determine grid spacing for different refinement levels
   !
   do iref = 1, nref
      !
      dxr(iref) = dx/(2**(iref - 1))
      dyr(iref) = dy/(2**(iref - 1))
      dxyr(iref) = min(dxr(iref), dyr(iref))      
      !
   enddo 
   !
   ! Calculate x,y coordinates of cell centres
   !
   allocate(z_xz(np))
   allocate(z_yz(np))
   !
   do nm = 1, np
      !
      n    = z_index_z_n(nm)
      m    = z_index_z_m(nm)
      iref = z_flags_iref(nm)
      !
      z_xz(nm) = x0 + cosrot * (1.0 * (m - 0.5)) * dxr(iref) - sinrot * (1.0 * (n - 0.5)) * dyr(iref)
      z_yz(nm) = y0 + sinrot * (1.0 * (m - 0.5)) * dxr(iref) + cosrot * (1.0 * (n - 0.5)) * dyr(iref)
      !
   enddo
   !
   ! And now the grid spacing in metres
   !
   dxymin = 1.0e6
   !
   if (crsgeo) then
      !
      ! dxr and dyr are now in degrees. Must convert to metres.
      ! For dxr, we need a value for every uv point. For dyr, just multiply by 111111.1
      !
      allocate(dxminv(npuv))   
      allocate(dxm2inv(npuv))
      !
      allocate(dyrm(nref))
      allocate(dyrinv(nref))
      allocate(dyr2inv(nref))
      allocate(dyrinvc(nref))
      !      
      allocate(cell_area_m2(np))
      allocate(dxm(np)) 
      allocate(fcorio2d(np))
      !
      dyrinvc = 0.0
      !
      fcorio2d = 0.0
      !
      do ip = 1, npuv
         !
         iref = uv_flags_iref(ip)
         nm   = uv_index_z_nm(ip)
         !
         dxminv(ip)  = 1.0 / ( dxr(iref) * 111111.1 * cos(z_yz(nm) * pi / 180) )
         dxm2inv(ip) = dxminv(ip)**2
         !
         ! Minimum grid spacing to determine dtmax 
         !  
         dxymin = min(dxymin, 1.0 / dxminv(ip))
         !
      enddo   
      !
      do iref = 1, nref
         !
         dyrm(iref)    = 111111.1*dyr(iref)
         dyrinv(iref)  = 1.0/dyrm(iref)
         dyr2inv(iref) = dyrinv(iref)**2
         !
         if (iref>1) then
            !         
            dyrinvc(iref) = 1.0 / (3 * dyrm(iref) / 2)
            !
         endif   
         !
         ! Minimum grid spacing to determine dtmax 
         !  
         dxymin = min(dxymin, dyrm(iref))
         !
      enddo   
      !
      ! Determine cell areas
      !
      do nm = 1, np
         !
         iref = z_flags_iref(nm)
         !
         dxm(nm)          = dxr(iref) * 111111.1 * cos(z_yz(nm) * pi / 180)
         cell_area_m2(nm) = dxm(nm) * dyrm(iref) 
         !         
      enddo   
      !
      ! Spatially varying Coriolis
      !
      if (coriolis) then
         do nm = 1, np            
            fcorio2d(nm) = 2 * 7.2921e-05 * sin(z_yz(nm) * pi / 180)
         enddo
      endif   
      !
   else 
      !
      ! Projected
      !
      allocate(dxrm(nref))
      allocate(dxrinv(nref))
      allocate(dxr2inv(nref))
      allocate(dxrinvc(nref))
      !
      allocate(dyrm(nref))
      allocate(dyrinv(nref))
      allocate(dyr2inv(nref))
      allocate(dyrinvc(nref))
      !
      allocate(cell_area(nref))
      !
      dxrinvc = 0.0
      dyrinvc = 0.0
      !
      do iref = 1, nref
         !
         dxrm(iref)      = dxr(iref)
         dxrinv(iref)    = 1.0 / dxrm(iref)
         dxr2inv(iref)   = dxrinv(iref)**2
         !
         dyrm(iref)      = dyr(iref)
         dyrinv(iref)    = 1.0 / dyrm(iref)
         dyr2inv(iref)   = dyrinv(iref)**2
         !
         cell_area(iref) = dxrm(iref)*dyrm(iref)
         !
         if (iref>1) then
            !         
            dxrinvc(iref) = 1.0 / (3 * dxrm(iref) / 2)
            dyrinvc(iref) = 1.0 / (3 * dyrm(iref) / 2)
            !
         endif   
         !
         ! Minimum grid spacing to determine dtmax 
         !  
         dxymin = min(dxymin, dxr(iref))
         dxymin = min(dxymin, dyr(iref))
         !
      enddo   
      !
      ! Determine minimum and maximum time step
      !
      ! dtmax will ensure that time step does not become too large because all cells are dry
      !
      dtmax = min(dtmax, alfa * dxymin / (sqrt(9.81 * hmin_cfl)))
      !
      ! Compute dtmin (if dt falls below this value, the simulation will stop).
      !
      ! No longer use long wave celerity, OLD implementation with 'stopdepth'
      !
      ! dtmin = alfa * dxymin / (sqrt(9.81 * stopdepth))
      !
      ! Current velocity
      !
      ! dtmin = min(dtmin, alfa * dxymin / (1.25 * abs(uvmax)))
      !
      ! Would be better to stop the simulation if water level change exceeds a certain value
      ! But use uvmax for now (which is typically set to a very high value)
      !       
      dtmin = alfa * dxymin / (1.25 * abs(uvmax))
      !
   endif   
   !
   ! Determine viscosity per refinement level
   !
   allocate(nuvisc(nref))  
   !
   if (viscosity) then
        !
        if (crsgeo) then
            !
            do iref = 1, nref
                !
                nuvisc(iref) = max(nuviscdim * dxyr(iref) / 0.00001, 0.0) ! take min of dx and dy, don't allow to be negative 
                !
                ! with default nuviscdim = 0.01 (dimensionless per meter grid cell length):                
                ! dx = 1 degree ~ 100km > nuvisc = 1000.0
                ! dx = 0.001 degree~ 100m > nuvisc = 1.0
                !
            enddo           
            !
        else
            !
            do iref = 1, nref
                !       
                nuvisc(iref) = max(nuviscdim * dxyr(iref), 0.0) ! take min of dx and dy, don't allow to be negative    
                !
                ! with default nuviscdim = 0.01 (dimensionless per meter grid cell length): 
                ! dx = 50 > nuvisc = 0.5
                ! dx = 100 > nuvisc = 1.0
                ! dx = 500 > nuvisc = 5.0
                !
            enddo                    
            ! 
        endif
        !
        if (use_quadtree) then
            write(logstr,'(a,e10.3,a,e10.3,a)')'Info    : viscosity - nuvisc  = [', minval(nuvisc),'- ',maxval(nuvisc),']'           
            call write_log(logstr, 0)
        else
            write(logstr,'(a,e10.3)')'Info    : viscosity - nuvisc  = ', maxval(nuvisc)           
            call write_log(logstr, 0)
        endif            
        !
   endif
   !
   end subroutine

   subroutine initialize_bathymetry()
   !
   ! Read the bed levels (dep file, quadtree file or subgrid tables) and set them at cells
   ! and uv points.
   !
   ! Called from: initialize_domain (sfincs_domain)
   !
   use sfincs_data
   use sfincs_subgrid
   use quadtree
   !
   implicit none
   !
   real*4, dimension(:,:),   allocatable :: zbg
   !
   integer :: ip
   integer :: nm
   integer :: nmu
   integer :: n
   integer :: m
   !
   logical :: ok
   !
   ! DEPTHS
   !
   if (.not. subgrid) then
      !
      ! Depths
      !   
      allocate(zb(np))
      allocate(zbuv(npuv))
      allocate(zbuvmx(npuv))
      !
      ! Read binary depth file
      !
      if (use_quadtree) then
         !
         ! If dep file not provided, use z from quadtree
         !
         if (depfile(1:4) /= 'none') then
            !
            write(logstr,'(a,a)')'Info    : reading dep file : ',trim(depfile)
            call write_log(logstr, 0)
            !
            ok = check_file_exists(depfile, 'Dep file', .true.)
            !
            open(unit = 500, file = trim(depfile), form = 'unformatted', access = 'stream')
            read(500)zb
            close(500)
            !
         else
            !
            do ip = 1, np
               zb(ip) = quadtree_zz(index_quadtree_in_sfincs(ip))
            enddo   
            !
         endif
         !
      else
         !
         write(logstr,'(a,a)')'Reading dep file : ',trim(depfile)
         call write_log(logstr, 0)
         !
         if (inputtype=='asc') then
            !
            allocate(zbg(nmax, mmax))
            !
            ok = check_file_exists(depfile, 'Dep file', .true.)
            !
            open(500, file=trim(depfile))
            do n = 1, nmax
               read(500,*)(zbg(n, m), m = 1, mmax)
            enddo
            close(500)
            !
            ip = 0
            do m = 1, mmax
               do n = 1, nmax
                  if (kcsg(n, m) > 0) then
                     ip = ip + 1
                     zb(ip)      = zbg(n, m)
                  endif                  
               enddo
            enddo
            !
            deallocate(zbg)
            !
         else
            !
            ok = check_file_exists(depfile, 'Dep file', .true.)
            !
            open(unit = 500, file = trim(depfile), form = 'unformatted', access = 'stream')
            read(500)zb
            close(500)
            !
         endif 
      endif   
      !
      ! Determine depths at u and v points
      !
      do ip = 1, npuv
         !
         nm  = uv_index_z_nm(ip)
         nmu = uv_index_z_nmu(ip)
         !
         if (uv_flags_type(ip)==0) then
            !
            ! Same refinement level
            !
            zbuv(ip) = 0.5*(zb(nm) + zb(nmu))
            !
         elseif (uv_flags_type(ip)==-1) then
            !
            ! Fine to coarse
            !
            zbuv(ip) = 0.3333*zb(nm) + 0.6667*zb(nmu)
            !
         else
            !
            ! Coarse to fine
            !
            zbuv(ip) = 0.6667*zb(nm) + 0.3333*zb(nmu)
            !
         endif
         !
         zbuvmx(ip) = max(zb(nm), zb(nmu)) + huthresh
         !
      enddo 
      !
   else
      !
      ! Read subgrid data
      !
      call read_subgrid_file()
      !
      ! In case of nonh or the velocity-form momentum scheme, we also need a per-cell
      ! bed level. Subgrid models have no reliable file zb, so use the EFFECTIVE bed
      ! implied by the subgrid tables at full wetness, zb = z_zmax - z_volmax/area:
      ! the cell-mean bed. From the first time step onward continuity maintains the
      ! water-level-dependent effective bed zb = zs - z_volume/area (zb_effective).
      ! The non-hydrostatic initialization freezes its w_b bed slopes from the
      ! cell-mean bed computed here.
      !
      if (nonhydrostatic .or. momentum_scheme == 1) then
         !
         allocate(zb(np))
         !
         do ip = 1, np
            !
            if (crsgeo) then
               zb(ip) = subgrid_z_zmax(ip) - subgrid_z_volmax(ip) / cell_area_m2(ip)
            else
               zb(ip) = subgrid_z_zmax(ip) - subgrid_z_volmax(ip) / cell_area(z_flags_iref(ip))
            endif
            !
         enddo
         !
      endif
      !
   endif
   !
   if (allocated(kcsg)) deallocate(kcsg)
   !
   end subroutine

   subroutine update_bed_level()
   !
   ! Apply the externally-supplied delta bed level (dzbext) to the kernel's
   ! bed-level arrays and refresh derived quantities at uv points.
   !
   ! Non-subgrid mode:
   !    zb           = zb           + dzbext           (cell centres)
   !    zbuvmx(ip)   = max(zb(nm), zb(nmu)) + huthresh (uv points, rebuilt)
   !
   ! Subgrid mode:
   !    subgrid_z_zmin/zmax shift rigidly by dzbext at the cell centre.
   !    subgrid_uv_zmin/zmax shift rigidly by the average dzbext of the two
   !    neighbouring cells, with subgrid_uv_zmin clamped from below by the
   !    larger of the two updated cell-centre subgrid_z_zmin values so the uv
   !    minimum can never sit below either neighbour's minimum.
   !
   ! Surface lift (always): the increment is depth-filtered into dzbext_filtered
   ! (filter_moving_bed_increment: wet cells, round-off increments dropped,
   ! Gaussian of movingbed_filter x depth of the moving cells, 0 = raw, volume of
   ! the rising and sinking parts kept) and the free surface of the wet cells is
   ! lifted with it (stamp_moving_bed_surface):
   !    non-subgrid: zs = zs + dzbext_filtered
   !    subgrid:     the shifted tables already imply zs + dzbext at the unchanged
   !                 z_volume; the volume is corrected by V(zs + dzbext_filtered) -
   !                 V(zs + dzbext) in the shifted table and zs set to the target
   ! Dry cells and cells whose increment is round-off follow the bed in both
   ! modes (their depth is unchanged). A moving bed therefore needs no mass
   ! bookkeeping by the caller on any grid type (a caller must not lift zs as
   ! well). dzbext_filtered is also the non-hydrostatic moving-bed source
   ! (d(zb)/dt = dzbext_filtered/dt in sfincs_nonhydrostatic, step 2b), so the
   ! surface lift and the pressure source share the one filter width.
   !
   ! The bed level at observation points (zbobs) is refreshed as well.
   !
   ! The caller (Python via BMI) owns the lifecycle of dzbext: this routine
   ! does not zero it out after applying.  When use_dzbext is .false. the
   ! routine still rebuilds zbuvmx in non-subgrid mode.
   !
   ! Called from: bmi_update_bed_level (sfincs_bmi)
   !
   use sfincs_data
   !
   integer :: ip
   integer :: nm
   integer :: nmu
   integer :: ilevel
   integer :: iobs
   real*4  :: avg_dzb
   !
   if (.not. subgrid) then
      !
      ! Non-subgrid path: shift zb, then rebuild zbuvmx for every uv point.
      !
      if (use_dzbext) then
         !
         zb(:) = zb(:) + dzbext(:)
         !
      endif
      !
      do ip = 1, npuv
         !
         nm  = uv_index_z_nm(ip)
         nmu = uv_index_z_nmu(ip)
         !
         zbuvmx(ip) = max(zb(nm), zb(nmu)) + huthresh
         !
      enddo
      !
   else
      !
      ! Subgrid path: only do anything when an external delta has been set.
      !
      if (use_dzbext) then
         !
         ! Cell-centre arrays shift rigidly with dzbext.
         !
         subgrid_z_zmin(:) = subgrid_z_zmin(:) + dzbext(:)
         subgrid_z_zmax(:) = subgrid_z_zmax(:) + dzbext(:)
         !
         ! The per-level volume->water-level table holds ABSOLUTE levels, so it
         ! must translate with the bed as well -- otherwise partially-wet cells
         ! (0 < z_volume < volmax, the interpolation branch in
         ! compute_water_levels_subgrid) would not follow the moving bed.  A rigid
         ! vertical shift leaves the volume bins (volmax / dzvol) unchanged, so
         ! only the levels move.  Keep the -20 m floor used in the level lookups.
         !
         do ilevel = 1, subgrid_nlevels
            subgrid_z_dep(ilevel, :) = max(subgrid_z_dep(ilevel, :) + dzbext(:), -20.0)
         enddo
         !
         ! UV-point arrays shift by the average delta of the two neighbours,
         ! then clamp uv_zmin from below by the higher of the two updated
         ! cell-centre minima.
         !
         do ip = 1, npuv
            !
            nm  = uv_index_z_nm(ip)
            nmu = uv_index_z_nmu(ip)
            !
            avg_dzb = 0.5 * (dzbext(nm) + dzbext(nmu))
            !
            subgrid_uv_zmin(ip) = subgrid_uv_zmin(ip) + avg_dzb
            subgrid_uv_zmax(ip) = subgrid_uv_zmax(ip) + avg_dzb
            !
            subgrid_uv_zmin(ip) = max(subgrid_uv_zmin(ip), &
                                      max(subgrid_z_zmin(nm), subgrid_z_zmin(nmu)))
            !
            ! The fully-wet face depth in momentum is hu = subgrid_uv_havg_zmax + zsu, with
            ! subgrid_uv_havg_zmax = havg(zmax) - zmax precomputed when the tables are read.
            ! Refresh it so it follows the shifted zmax; otherwise every fully-wet face keeps
            ! an hu that is off by the bed change, which goes negative on flat faces where the
            ! bed subsided (uv = q / max(hu, huvmin) then divides by zero and dt collapses).
            !
            subgrid_uv_havg_zmax(ip) = subgrid_uv_havg(subgrid_nlevels, ip) - subgrid_uv_zmax(ip)
            !
         enddo
         !
      endif
      !
   endif
   !
   ! Depth-filtered increment of the wet cells (the non-hydrostatic moving-bed
   ! source reads it as well), then the surface lift.
   !
   if (use_dzbext) then
      !
      call filter_moving_bed_increment()
      !
      call stamp_moving_bed_surface()
      !
   endif
   !
   ! Refresh the bed level at observation points, so that the water depth
   ! written to the his file (hobs = zs - zbobs) follows the moving bed.
   ! Points outside the domain (nmindobs = 0) keep their -999.0 value.
   !
   if (use_dzbext .and. nobs > 0 .and. allocated(zbobs)) then
      !
      do iobs = 1, nobs
         !
         nm = nmindobs(iobs)
         !
         if (nm > 0) then
            !
            if (subgrid) then
               zbobs(iobs) = subgrid_z_zmin(nm)
            else
               zbobs(iobs) = zb(nm)
            endif
            !
         endif
         !
      enddo
      !
   endif
   !
   end subroutine update_bed_level

   !
   !
   subroutine filter_moving_bed_increment()
   !
   ! Builds dzbext_filtered from the caller's bed increment dzbext: the increment of
   ! the wet cells that move (moving_bed_cell_moves: depth > max(huthresh,
   ! movingbed_hwet_min) and increment above nh_dzb_ulps float32 ulps of the bed
   ! level; smaller increments are round-off of the caller's float64 target bed
   ! minus the float32 kernel bed, and as a bed velocity on the 1e-6 s start-up
   ! step such noise is an impulse, see sfincs_nonhydrostatic step 2b),
   ! depth-filtered with depth_filter_cell_field when movingbed_filter > 0. Other
   ! cells get 0.
   !
   ! Called from: update_bed_level
   !
   use sfincs_data
   !
   integer :: nm
   real*4, dimension(:), allocatable :: fin
   !
   if (.not. allocated(dzbext_filtered)) then
      !
      allocate(dzbext_filtered(np))
      dzbext_filtered = 0.0
      !
   endif
   !
   allocate(fin(np))
   !
   !$omp parallel do schedule ( static ) private ( nm )
   do nm = 1, np
      !
      fin(nm) = 0.0
      !
      if (moving_bed_cell_moves(nm)) then
         !
         fin(nm) = dzbext(nm)
         !
      endif
      !
   enddo
   !$omp end parallel do
   !
   if (movingbed_filter > 0.0) then
      !
      call depth_filter_cell_field(fin, dzbext_filtered, movingbed_filter)
      !
   else
      !
      dzbext_filtered = fin
      !
   endif
   !
   deallocate(fin)
   !
   end subroutine filter_moving_bed_increment
   !
   !
   subroutine stamp_moving_bed_surface()
   !
   ! Lifts the free surface of the wet cells with the filtered bed increment (which
   ! extends beyond the moving footprint) plus, where the cell's own increment is
   ! round-off (below the filter's motion threshold, moving_bed_cell_moves, but
   ! possibly still a change of the float32 bed), that raw increment, so the
   ! depth there is unchanged. Dry cells follow the bed with the raw increment
   ! (their depth, zero or below the threshold, is unchanged: a dry cell whose bed
   ! sinks must not fill itself).
   !
   ! Subgrid: update_bed_level has shifted the level tables by dzbext, so the
   ! unchanged z_volume now maps to zs + dzbext (every branch of the table lookup
   ! in compute_water_levels_subgrid is a rigid shift). The volume is corrected by
   ! V(zs + dzbext_filtered) - V(zs + dzbext) of the shifted table
   ! (subgrid_volume_at_level: a x dzs above z_zmax, inverse table lookup below
   ! it as compute_initial_subgrid_volumes), zs is set to the target and the
   ! effective bed refreshed as continuity does (zb_effective).
   !
   ! Host arrays only, as the rest of update_bed_level (on a GPU build the device
   ! copies of zs, z_volume and zb are not refreshed here).
   !
   ! Called from: update_bed_level
   !
   use sfincs_data
   !
   integer :: nm
   real*4  :: a
   real*4  :: hwet
   real*4  :: dzs
   real*8  :: zs_raw
   real*8  :: zs_new
   !
   hwet = max(huthresh, movingbed_hwet_min)
   !
   if (.not. subgrid) then
      !
      !$omp parallel do schedule ( static ) private ( nm, dzs )
      do nm = 1, np
         !
         if (moving_bed_cell_depth(nm) > hwet) then
            !
            dzs = stamp_increment(nm)
            !
            zs(nm) = zs(nm) + dzs
            !
         else
            !
            zs(nm) = zs(nm) + dzbext(nm)
            !
         endif
         !
      enddo
      !$omp end parallel do
      !
   else
      !
      !$omp parallel do schedule ( static ) private ( nm, a, dzs, zs_raw, zs_new )
      do nm = 1, np
         !
         a = moving_bed_cell_area(nm)
         !
         if (moving_bed_cell_depth(nm) > hwet) then
            !
            dzs          = stamp_increment(nm)
            zs_raw       = zs(nm) + dzbext(nm)
            zs_new       = zs(nm) + dzs
            z_volume(nm) = z_volume(nm) + (subgrid_volume_at_level(nm, zs_new, a) - subgrid_volume_at_level(nm, zs_raw, a))
            zs(nm)       = zs_new
            !
         else
            !
            zs(nm) = zs(nm) + dzbext(nm)
            !
         endif
         !
         if (zb_effective) then
            !
            zb(nm) = zs(nm) - max(z_volume(nm) / a, 0.0d0)
            !
         endif
         !
      enddo
      !$omp end parallel do
      !
   endif
   !
   contains
      !
      function stamp_increment(i) result(dz)
         !
         ! Surface increment of a wet cell: the filtered field, plus the cell's
         ! own raw increment when that is round-off (not in the filtered field,
         ! but possibly a change of the float32 bed)
         !
         integer, intent(in) :: i
         real*4 :: dz
         !
         dz = dzbext_filtered(i)
         !
         if (.not. moving_bed_cell_moves(i)) dz = dz + dzbext(i)
         !
      end function stamp_increment
      !
   end subroutine stamp_moving_bed_surface
   !
   !
   subroutine depth_filter_cell_field(fin, fout, width)
   !
   ! Depth filter (Kajiura-type attenuation) of a cell field, e.g. a bed increment:
   ! Gaussian of standard deviation width x hc, where hc is the |fin|-weighted mean
   ! depth of the non-zero cells, with a floor of one cell. Separable x/y passes
   ! along the same-level neighbour chains (the single faces z_index_uv_md1/mu1/
   ! nd1/nu1 of type 0), scattered from the non-zero cells only, so the cost scales
   ! with the footprint and the kernel radius (serial host routine: the scatter
   ! writes to neighbour cells). The rising (> 0) and sinking (< 0) parts are
   ! filtered separately, masked to the wet cells and rescaled to their own
   ! area-weighted totals, so each keeps its volume (a translating slide has a net
   ! volume of ~0 and would lose its shape in a single rescaling).
   !
   ! Quadtree grids: the cells of one refinement level form a regular lattice with
   ! holes. The chains stay within the level (they stop at refinement transitions,
   ! at the domain edge and at holes), the width in cells follows the level's
   ! spacing (per cell for crsgeo), so every level is filtered as a regular grid
   ! and the rescaling restores what the truncated chains lose. Mass is not
   ! exchanged across a transition: a feature that straddles one is smoothed on
   ! each side separately (a step of the filtered field at the transition, of the
   ! order of the smoothing itself), and the part lost at a transition or edge is
   ! returned to the whole part by the rescaling rather than locally. A non-zero
   ! cell without a neighbour in a direction keeps its value in that pass (nothing
   ! to spread to); cells without any neighbour are reported once in the log.
   !
   ! Called from: filter_moving_bed_increment (bed increment of a moving bed)
   !
   use sfincs_data
   use sfincs_log
   !
   real*4, dimension(:), intent(in)  :: fin
   real*4, dimension(:), intent(out) :: fout
   real*4,               intent(in)  :: width
   !
   integer, parameter :: rmax = 500
   !
   integer :: nm
   integer :: j
   integer :: k
   integer :: r
   integer :: ipart
   integer :: nochain
   real*4  :: hc
   real*4  :: sig
   real*4  :: sigprev
   real*4  :: hwet
   real*4  :: fac
   real*8  :: sw8
   real*8  :: swh8
   real*8  :: tot08
   real*8  :: tot8
   real*4, dimension(0:rmax) :: w
   real*4, dimension(:), allocatable :: part
   real*4, dimension(:), allocatable :: tmp
   real*4, dimension(:), allocatable :: outp
   logical, save :: warned = .false.
   !
   fout = 0.0
   hwet = max(huthresh, movingbed_hwet_min)
   !
   ! |fin|-weighted mean depth of the moving cells: the width scale
   !
   sw8     = 0.0d0
   swh8    = 0.0d0
   nochain = 0
   !
   do nm = 1, np
      !
      if (fin(nm) /= 0.0) then
         !
         sw8  = sw8 + abs(fin(nm))
         swh8 = swh8 + abs(fin(nm)) * max(moving_bed_cell_depth(nm), 0.0)
         !
         if (nb_mu(nm) == 0 .and. nb_md(nm) == 0 .and. nb_nu(nm) == 0 .and. nb_nd(nm) == 0) then
            !
            nochain = nochain + 1
            !
         endif
         !
      endif
      !
   enddo
   !
   if (sw8 <= 0.0d0) return
   !
   hc = real(swh8 / sw8, 4)
   !
   if (nochain > 0 .and. .not. warned) then
      !
      write(logstr,'(a,i0,a)')'Warning : moving-bed depth filter: ', nochain, ' moving cells have no same-level neighbour; their increment is not filtered'
      call write_log(logstr, 0)
      warned = .true.
      !
   endif
   !
   allocate(part(np))
   allocate(tmp(np))
   allocate(outp(np))
   !
   do ipart = 1, 2
      !
      ! Rising (ipart = 1) or sinking (ipart = 2) part and its volume
      !
      tot08 = 0.0d0
      !
      do nm = 1, np
         !
         if (ipart == 1) then
            !
            part(nm) = max(fin(nm), 0.0)
            !
         else
            !
            part(nm) = min(fin(nm), 0.0)
            !
         endif
         !
         tot08 = tot08 + part(nm) * moving_bed_cell_area(nm)
         !
      enddo
      !
      if (tot08 == 0.0d0) cycle
      !
      ! x pass: scatter every non-zero cell along its west/east chain. The width in
      ! cells is the cell's own (level or, for crsgeo, latitude dependent) spacing,
      ! floor one cell; the weights are rebuilt when it changes.
      !
      tmp     = 0.0
      sigprev = -1.0
      r       = 0
      !
      do nm = 1, np
         !
         if (part(nm) == 0.0) cycle
         !
         if (nb_mu(nm) == 0 .and. nb_md(nm) == 0) then
            !
            tmp(nm) = tmp(nm) + part(nm)
            cycle
            !
         endif
         !
         sig = width * max(hc * moving_bed_dxinv(nm), 1.0)
         !
         if (sig /= sigprev) then
            !
            r = min(ceiling(4.0 * sig), rmax)
            call gaussian_weights(sig, r, w)
            sigprev = sig
            !
         endif
         !
         tmp(nm) = tmp(nm) + w(0) * part(nm)
         j       = nm
         !
         do k = 1, r
            !
            j = nb_mu(j)
            if (j == 0) exit
            tmp(j) = tmp(j) + w(k) * part(nm)
            !
         enddo
         !
         j = nm
         !
         do k = 1, r
            !
            j = nb_md(j)
            if (j == 0) exit
            tmp(j) = tmp(j) + w(k) * part(nm)
            !
         enddo
         !
      enddo
      !
      ! y pass (2-D grids): scatter along the south/north chains
      !
      if (nmax > 1) then
         !
         outp    = 0.0
         sigprev = -1.0
         r       = 0
         !
         do nm = 1, np
            !
            if (tmp(nm) == 0.0) cycle
            !
            if (nb_nu(nm) == 0 .and. nb_nd(nm) == 0) then
               !
               outp(nm) = outp(nm) + tmp(nm)
               cycle
               !
            endif
            !
            sig = width * max(hc * dyrinv(z_flags_iref(nm)), 1.0)
            !
            if (sig /= sigprev) then
               !
               r = min(ceiling(4.0 * sig), rmax)
               call gaussian_weights(sig, r, w)
               sigprev = sig
               !
            endif
            !
            outp(nm) = outp(nm) + w(0) * tmp(nm)
            j        = nm
            !
            do k = 1, r
               !
               j = nb_nu(j)
               if (j == 0) exit
               outp(j) = outp(j) + w(k) * tmp(nm)
               !
            enddo
            !
            j = nm
            !
            do k = 1, r
               !
               j = nb_nd(j)
               if (j == 0) exit
               outp(j) = outp(j) + w(k) * tmp(nm)
               !
            enddo
            !
         enddo
         !
      else
         !
         outp = tmp
         !
      endif
      !
      ! Wet cells only, rescaled to the volume of this part
      !
      tot8 = 0.0d0
      !
      do nm = 1, np
         !
         if (moving_bed_cell_depth(nm) <= hwet) outp(nm) = 0.0
         tot8 = tot8 + outp(nm) * moving_bed_cell_area(nm)
         !
      enddo
      !
      if (tot8 == 0.0d0) cycle
      !
      fac = real(tot08 / tot8, 4)
      !
      do nm = 1, np
         !
         fout(nm) = fout(nm) + outp(nm) * fac
         !
      enddo
      !
   enddo
   !
   deallocate(part, tmp, outp)
   !
   contains
      !
      function nb_mu(i) result(inb)
         !
         ! Same-level neighbour to the right (0 = none: domain edge, hole or refinement transition)
         !
         integer, intent(in) :: i
         integer :: inb
         integer :: ip
         !
         inb = 0
         ip  = z_index_uv_mu1(i)
         !
         if (ip > 0) then
            !
            if (uv_flags_type(ip) == 0) inb = uv_index_z_nmu(ip)
            !
         endif
         !
      end function nb_mu
      !
      function nb_md(i) result(inb)
         !
         integer, intent(in) :: i
         integer :: inb
         integer :: ip
         !
         inb = 0
         ip  = z_index_uv_md1(i)
         !
         if (ip > 0) then
            !
            if (uv_flags_type(ip) == 0) inb = uv_index_z_nm(ip)
            !
         endif
         !
      end function nb_md
      !
      function nb_nu(i) result(inb)
         !
         integer, intent(in) :: i
         integer :: inb
         integer :: ip
         !
         inb = 0
         ip  = z_index_uv_nu1(i)
         !
         if (ip > 0) then
            !
            if (uv_flags_type(ip) == 0) inb = uv_index_z_nmu(ip)
            !
         endif
         !
      end function nb_nu
      !
      function nb_nd(i) result(inb)
         !
         integer, intent(in) :: i
         integer :: inb
         integer :: ip
         !
         inb = 0
         ip  = z_index_uv_nd1(i)
         !
         if (ip > 0) then
            !
            if (uv_flags_type(ip) == 0) inb = uv_index_z_nm(ip)
            !
         endif
         !
      end function nb_nd
      !
      subroutine gaussian_weights(sg, rr, ww)
         !
         ! One-sided normalized Gaussian weights ww(0:rr), ww(0) + 2 sum ww(1:rr) = 1.
         !
         real*4,  intent(in)  :: sg
         integer, intent(in)  :: rr
         real*4,  intent(out) :: ww(0:)
         !
         integer :: kk
         real*4  :: wsum
         !
         if (sg < 0.1 .or. rr < 1) then
            !
            ww(0:rr) = 0.0
            ww(0)    = 1.0
            return
            !
         endif
         !
         do kk = 0, rr
            !
            ww(kk) = exp(-0.5 * (real(kk) / sg)**2)
            !
         enddo
         !
         wsum     = ww(0) + 2.0 * sum(ww(1:rr))
         ww(0:rr) = ww(0:rr) / wsum
         !
      end subroutine gaussian_weights
      !
   end subroutine depth_filter_cell_field
   !
   !
   function moving_bed_cell_moves(nm) result(moves)
   !
   ! Whether the bed increment dzbext of a cell counts as motion of a wet cell:
   ! depth > max(huthresh, movingbed_hwet_min) and |dzbext| above nh_dzb_ulps
   ! float32 ulps of the bed level (zb, or subgrid_z_zmin on subgrid grids).
   !
   use sfincs_data
   use sfincs_nonhydrostatic, only: nh_dzb_ulps
   !
   integer, intent(in) :: nm
   logical :: moves
   !
   real*4 :: zbed
   !
   if (subgrid) then
      !
      zbed = subgrid_z_zmin(nm)
      !
   else
      !
      zbed = zb(nm)
      !
   endif
   !
   moves = abs(dzbext(nm)) > nh_dzb_ulps * spacing(max(abs(zbed), 1.0)) .and. &
           moving_bed_cell_depth(nm) > max(huthresh, movingbed_hwet_min)
   !
   end function moving_bed_cell_moves
   !
   !
   function moving_bed_cell_depth(nm) result(h)
   !
   ! Water depth of a cell for the moving-bed stamp and filter: zs - zb on regular
   ! grids, the volume-mean depth z_volume/area on subgrid grids (independent of zb,
   ! which is the effective bed only with zb_effective).
   !
   use sfincs_data
   !
   integer, intent(in) :: nm
   real*4 :: h
   !
   if (subgrid) then
      !
      h = real(z_volume(nm), 4) / moving_bed_cell_area(nm)
      !
   else
      !
      h = real(zs(nm) - zb(nm), 4)
      !
   endif
   !
   end function moving_bed_cell_depth
   !
   !
   function moving_bed_cell_area(nm) result(a)
   !
   use sfincs_data
   !
   integer, intent(in) :: nm
   real*4 :: a
   !
   if (crsgeo) then
      !
      a = cell_area_m2(nm)
      !
   else
      !
      a = cell_area(z_flags_iref(nm))
      !
   endif
   !
   end function moving_bed_cell_area
   !
   !
   function moving_bed_dxinv(nm) result(rdx)
   !
   ! 1/dx in metres of a cell: per refinement level when projected, per cell
   ! (latitude dependent) when crsgeo
   !
   use sfincs_data
   !
   integer, intent(in) :: nm
   real*4 :: rdx
   !
   if (crsgeo) then
      !
      rdx = 1.0 / dxm(nm)
      !
   else
      !
      rdx = dxrinv(z_flags_iref(nm))
      !
   endif
   !
   end function moving_bed_dxinv
   !
   !
   function subgrid_volume_at_level(nm, zsv, a) result(vol)
   !
   ! Cell volume of the subgrid table at water level zsv (inverse of the lookup in
   ! compute_water_levels_subgrid, as compute_initial_subgrid_volumes)
   !
   use sfincs_data
   !
   integer, intent(in) :: nm
   real*8,  intent(in) :: zsv
   real*4,  intent(in) :: a
   real*8 :: vol
   !
   integer :: ivol
   integer :: ilevel
   real*4  :: dzvol
   real*4  :: facint
   !
   if (zsv >= subgrid_z_zmax(nm)) then
      !
      vol = subgrid_z_volmax(nm) + a * (zsv - max(subgrid_z_zmax(nm), -20.0))
      !
   else
      !
      ivol = 1
      !
      do ilevel = 2, subgrid_nlevels
         !
         if (subgrid_z_dep(ilevel, nm) > zsv) then
            !
            ivol = ilevel - 1
            exit
            !
         endif
         !
      enddo
      !
      dzvol  = subgrid_z_volmax(nm) / (subgrid_nlevels - 1)
      facint = real(zsv - subgrid_z_dep(ivol, nm), 4) / max(subgrid_z_dep(ivol + 1, nm) - subgrid_z_dep(ivol, nm), 0.001)
      facint = min(max(facint, 0.0), 1.0)
      vol    = (ivol - 1) * dzvol + facint * dzvol
      !
   endif
   !
   end function subgrid_volume_at_level

   subroutine initialize_boundaries()
   !
   ! Count the boundary cells per type and set up the boundary index arrays.
   !
   ! Called from: initialize_domain (sfincs_domain)
   !
   use sfincs_data
   !
   implicit none
   !
   integer :: ip
   integer :: ib
   integer :: nm
   integer :: nmu
   integer :: ikcuv2
   !
   ! BOUNDARIES
   !
   ! First count number of boundary cells
   !
   ! kcs = 2 : water level
   ! kcs = 3 : outflow
   ! kcs = 5 : river
   ! kcs = 6 : lateral (coastal, dzs/dx = 0.0), in which case we set kcuv = 6 and we do not add this point to the boundary u/v points ! We do add this point to the grid boundary points (later on).
   !
   ngbnd = 0
   !
   do nm = 1, np
      if (kcs(nm) == 2 .or. kcs(nm) == 3 .or. kcs(nm) == 5 .or. kcs(nm) == 6) then
         !
         boundaries_in_mask = .true.
         ngbnd = ngbnd + 1
         !
         if (kcs(nm) == 2) water_level_boundaries_in_mask = .true.
         if (kcs(nm) == 3) outflow_boundaries_in_mask = .true.
         if (kcs(nm) == 5) downstream_river_boundaries_in_mask = .true.
         if (kcs(nm) == 6) neumann_boundaries_in_mask = .true.
         ! Neumann - Need nm indices of internal points, determined in sfincs_boundaries.f90
         !
      endif
   enddo
   !
   ! Now allocate them
   !
   allocate(nmindbnd(ngbnd))
   allocate(zsb(ngbnd))
   allocate(zsb0(ngbnd))
   !
   zsb  = 0.0  ! Total water level at boundary grid point
   zsb0 = 0.0  ! Filtered water level at boundary grid point
   !
   ! And now set the nm index for each boundary point
   !
   ngbnd = 0
   !
   do nm = 1, np
      !
      if (kcs(nm) == 2 .or. kcs(nm) == 3 .or. kcs(nm) == 5 .or. kcs(nm) == 6) then
         !
         ! This is a boundary point
         !
         ngbnd = ngbnd + 1
         nmindbnd(ngbnd) = nm
         !
      endif
   enddo
   !
   ! UV boundary points
   !
   nkcuv2 = 0 ! Number of points with kcuv=2 or kcuv=6
   kcuv   = 0
   !
   do ip = 1, npuv
      !
      nm  = uv_index_z_nm(ip)
      nmu = uv_index_z_nmu(ip)
      !
      if (kcs(nm) > 1 .and. kcs(nmu) > 1)  then
         !
         ! Two surrounding boundary points
         !
         kcuv(ip) = 0 ! Is this really necessary ?
         !
      elseif (kcs(nm) * kcs(nmu) == 0) then
         !
         ! One or two surrounding inactive points, so make this velocity point inactive
         !
         kcuv(ip) = 0
         !
      elseif ((kcs(nm) == 1 .and. kcs(nmu) == 2) .or. (kcs(nm) == 2 .and. kcs(nmu) == 1))  then
         !
         ! One surrounding water level boundary point, so make this velocity point a boundary point
         !
         kcuv(ip) = 2
         nkcuv2 = nkcuv2 + 1
         !
      elseif ((kcs(nm) == 1 .and. kcs(nmu) == 3) .or. (kcs(nm) == 3 .and. kcs(nmu) == 1))  then
         !
         ! One surrounding outflow boundary point, so make this velocity point a boundary point
         !
         kcuv(ip) = 2
         nkcuv2 = nkcuv2 + 1
         !
      elseif ((kcs(nm) == 1 .and. kcs(nmu) == 5) .or. (kcs(nm) == 5 .and. kcs(nmu) == 1))  then
         !
         ! One surrounding downstream river boundary point, so make this velocity point a boundary point
         !
         kcuv(ip) = 2
         nkcuv2 = nkcuv2 + 1
         !
      elseif ((kcs(nm) == 1 .and. kcs(nmu) == 6) .or. (kcs(nm) == 6 .and. kcs(nmu) == 1))  then
         !
         ! Lateral coastal boundary
         !
         kcuv(ip) = 6
         nkcuv2 = nkcuv2 + 1
         !         
      else
         !
         kcuv(ip) = 1
         !
      endif
      !
   enddo
   !
   ! Loop through boundary velocity points (kcuv=2)
   !
   allocate(index_kcuv2(nkcuv2))
   allocate(nmbkcuv2(nkcuv2))
   allocate(nmikcuv2(nkcuv2))
   allocate(ibuvdir(nkcuv2))
   allocate(uvmean(nkcuv2))
   allocate(ibkcuv2(nkcuv2))
   !
   uvmean = 0.0
   !
   ikcuv2 = 0
   !
   do ip = 1, npuv
      !
      if (kcuv(ip) == 2 .or. kcuv(ip) == 6) then
         !
         ikcuv2 = ikcuv2 + 1       ! Counter for kcuv==2 points
         index_kcuv2(ikcuv2) = ip
         !
         nm  = uv_index_z_nm(ip)
         nmu = uv_index_z_nmu(ip)
         !
         if (kcs(nm) == 1) then
            !
            ! Boundary point on the right
            !
            nmbkcuv2(ikcuv2) = nmu
            nmikcuv2(ikcuv2) = nm
            ibuvdir(ikcuv2)  = -1
            !
         else
            !
            ! Boundary point on the left
            !
            nmbkcuv2(ikcuv2) = nm
            nmikcuv2(ikcuv2) = nmu
            ibuvdir(ikcuv2)  = 1
            !
         endif         
         !
         ! Find index of boundary point in boundary array
         !         
         do ib = 1, ngbnd
            !
            if (nmindbnd(ib) == nmbkcuv2(ikcuv2)) then
               !
               ibkcuv2(ikcuv2) = ib ! index in grid boundary array of this uv boundary point
               !
               ! Check if this is an outflow boundary (kcs=3)
               !
               if (kcs(nmindbnd(ib)) == 3) then
                  !
                  ! Set water levels for this point (this only happens here, so these do not change throughout the simulation)
                  !
                  if (subgrid) then
                     zsb0(ib) = subgrid_z_zmin(nmindbnd(ib))
                     zsb(ib)  = zsb0(ib)
                  else    
                     zsb0(ib) = zb(nmindbnd(ib))
                     zsb(ib)  = zsb0(ib)
                  endif
               endif   
               !
               exit
               !
            endif
         enddo
         !
      endif
   enddo
   !
   end subroutine

   subroutine initialize_roughness()
   !
   ! Set the Manning roughness at the uv points (uniform, land/sea or from file).
   !
   ! Called from: initialize_domain (sfincs_domain)
   !
   use sfincs_data
   use sfincs_ncinput   
   use sfincs_read, only: is_netcdf_file
   !
   implicit none
   !
   !
   integer :: ip
   integer :: nm
   integer :: nmu
   logical :: ok
   !
   character*256 :: varname
   !
   ! FRICTION COEFFICIENTS (only for regular bathymetry, as for subgrid the Manning's n values are stored in the tables)
   !
   if (.not. subgrid) then
      !
      allocate(gn2uv(npuv))
      !
      gn2uv = 9.81*0.02*0.02
      !
      if (manningfile(1:4) /= 'none') then 
         !
         ! Read spatially-varying friction
         ! File is either binary or netcdf
         !
         allocate(rghfield(np))
         !
         write(logstr,'(a,a)')'Info    : reading roughness file ',trim(manningfile)
         call write_log(logstr, 0)
         !
         ok = check_file_exists(manningfile, 'Roughness file', .true.)
         !
         if (is_netcdf_file(manningfile)) then
            !
            ! Call the generic quadtree nc file reader function
            varname = 'manning'
            call read_netcdf_quadtree_to_sfincs(manningfile, varname, rghfield) !ncfile, varname, varout)      
            !
         else
            !
            ! Read from binary file
            !               
            open(unit = 500, file = trim(manningfile), form = 'unformatted', access = 'stream')
            read(500)rghfield
            close(500)
            !
         endif
         !
         do ip = 1, npuv
            nm  = uv_index_z_nm(ip)
            nmu = uv_index_z_nmu(ip)
            if (kcuv(ip)>0) then
               gn2uv(ip) = 9.81*(0.5*(rghfield(nm) + rghfield(nmu)))**2
            endif
         enddo   
         !
      else
         !
         ! Manning's n values for land and sea
         !
         if (manning_land<0.0) then
            gn2land = 9.81*manning**2
         else
            gn2land = 9.81*manning_land**2
         endif
         !
         if (manning_sea<0.0) then
            gn2sea = 9.81*manning**2
         else
            gn2sea = 9.81*manning_sea**2
         endif
         !
         do ip = 1, npuv
            if (kcuv(ip)>0) then
               if (zbuv(ip)>rghlevland) then
                  gn2uv(ip) = gn2land
               else
                  gn2uv(ip) = gn2sea
               endif
            endif
         enddo
         !
      endif
   else
      !
      ! Give warning if manningfile is supplied, but also a subgrid file
      !
      if (manningfile(1:4) /= 'none') then 
         !
         call write_log('Warning   : manningfile input will be ignored because SFINCS will use the friction information from sbgfile!', 1)
         !
      endif
       !
   endif
   !
   end subroutine


   subroutine initialize_storage_volume()
   !
   ! Read the spatially-varying storage volume per cell (e.g. for green infrastructure).
   !
   ! Called from: initialize_domain (sfincs_domain)
   !
   use sfincs_data
   use sfincs_ncinput
   use sfincs_read, only: is_netcdf_file
   !
   implicit none
   !
   logical :: ok
   character*256 :: varname
   !
   if (use_storage_volume) then 
      !
      ! Spatially-varying storage volume for green infra
      !
      allocate(storage_volume(np))
      !
      write(logstr,'(a)')'Info    : turning on Storage Green Infrastructure'      
      call write_log(logstr, 0)
      !
      ! Read spatially-varying storage per cell (m3)
      !
      write(logstr,'(a,a)')'Info    : reading vol file ',trim(volfile)
      call write_log(logstr, 0)
      !
      ok = check_file_exists(volfile, 'Storage volume vol file', .true.)
      !
      if (is_netcdf_file(volfile)) then
         !
         ! Call the generic quadtree nc file reader function
         varname = 'vol'
         call read_netcdf_quadtree_to_sfincs(volfile, varname, storage_volume) !ncfile, varname, varout)            
         !
      else
         !
         ! Read from binary file
         !      
         open(unit = 500, file = trim(volfile), form = 'unformatted', access = 'stream')
         read(500)storage_volume
         close(500)
         !
      endif   
      !
   endif
   !
   end subroutine


   subroutine set_advection_mask()
   !
   ! Set the mask that switches off advection at uv points next to boundaries and wave
   ! makers.
   !
   ! Called from: initialize_domain (sfincs_domain)
   !
   use sfincs_data
   !
   implicit none
   !
   integer :: icuv, ip
   !
   ! First set kcuv for combined points
   !
   if (ncuv > 0) then
      !
      ! Loop through combined uv points to set kcuv
      !
      do icuv = 1, ncuv
         !
         if (kcuv(cuv_index_uv1(icuv)) == 0 .or. kcuv(cuv_index_uv2(icuv)) == 0) then
            kcuv(cuv_index_uv(icuv)) = 0
         elseif (kcuv(cuv_index_uv1(icuv)) == 2 .or. kcuv(cuv_index_uv2(icuv)) == 2) then
            kcuv(cuv_index_uv(icuv)) = 2
         elseif (kcuv(cuv_index_uv1(icuv)) == 4 .or. kcuv(cuv_index_uv2(icuv)) == 4) then
            kcuv(cuv_index_uv(icuv)) = 4
         else
            kcuv(cuv_index_uv(icuv)) = 1
         endif
         !
      enddo
      !
   endif
   !
   ! And now set the mask
   ! The advection term at uv points that have neighbors that are 1) boundary points or 2) wave maker points will be set to 0.0
   !
   if (advection) then
      !
      ! Should only do this for quadtree or when nmax>1 (i.e. not for 1D profile models)
      !
      allocate(mask_adv(npuv))
      !
      mask_adv = 1 
      !
      if (advection_mask) then 
         if (use_quadtree .or. nmax>1) then
            !
            do ip = 1, npuv
               !
               if ( kcuv(uv_index_u_nmu(ip)) /= 1 .or. &
                    kcuv(uv_index_u_nmd(ip)) /= 1 .or. &
                    kcuv(uv_index_v_nm(ip)) /= 1  .or. &
                    kcuv(uv_index_v_nmu(ip)) /= 1 .or. &
                    kcuv(uv_index_v_ndm(ip)) /= 1 .or. &
                    kcuv(uv_index_v_ndmu(ip)) /= 1 .or. &
                    kcuv(uv_index_u_ndm(ip)) /= 1 .or. &
                    kcuv(uv_index_u_num(ip)) /= 1 ) then
                  !
                  mask_adv(ip) = 0
                  !
               endif
               !  
            enddo 
            ! 
         endif
      endif 
      ! 
   endif
   !
   end subroutine


   subroutine initialize_hydro()
   !
   ! Initialize SFINCS variables (qx, qy, zs etc.)
   !
   ! Called from: initialize_domain (sfincs_domain)
   !
   use sfincs_data
   use sfincs_subgrid
   use sfincs_initial_conditions
   !
   implicit none
   !
   integer    :: icuv
   !
   allocate(zs(np))
   !
   ! q, q0, uv, and uv0 also contain the combined values at quadtree transition cells. Plus 1 for dummy values set to 0.0.
   !
   allocate(q(npuv + ncuv + 1))
   allocate(q0(npuv + ncuv + 1))
   allocate(uv(npuv + ncuv + 1))
   allocate(uv0(npuv + ncuv + 1))
   !
   allocate(kfuv(npuv))
   !
   ! Cell-wise discharge accumulator (point sources + drainage structures),
   ! read by sfincs_continuity.
   !
   allocate(qsrc(np))
   !
   zs   = 0.0
   q    = 0.0
   q0   = 0.0
   uv   = 0.0
   uv0  = 0.0
   qsrc = 0.0
   !
   kfuv = 0 
   !
   if (wiggle_suppression) then 
      ! 
      allocate(zs0(np))
      allocate(zsderv(np))
      zs0 = 0.0
      zsderv = 0.0
      ! 
   endif
   !
   if (compute_wetfrac) then
      !
      allocate(z_wetfrac(np))
      z_wetfrac = 1.0
      !
   endif
   !
   if (slope_driven_flow) then
      !
      allocate(z_hwet(np))
      z_hwet = 0.0
      !
   endif
   !
   if (store_slope_regime) then
      !
      allocate(z_wface(np))
      allocate(w_uv(npuv))
      allocate(iup_uv(npuv))
      z_wface = 0.0
      w_uv = 0.0
      iup_uv = 0
      !
   endif
   !
   if (snapwave) then
      !
      allocate(hm0(np))
      allocate(hm0_ig(np))
      allocate(fwuv(npuv))
      !
      hm0    = 0.0
      hm0_ig = 0.0
      fwuv   = 0.0
      !
   endif
   !
   if (wavemaker .or. snapwave) then !TL: zsm also used in sfincs_continuity if 'snapwave=true' > todo: check if needed, or only for wavemaker 
      allocate(zsm(np))
      !
      if (store_maximum_waterlevel) then
          allocate(maxzsm(np))      
      endif      
   endif
   !
   if (store_maximum_waterlevel) then
      allocate(zsmax(np))
   endif
   !
   if (store_zvolume_max) then
      allocate(zvolmax(np))
   endif
   !
   if (store_maximum_velocity) then
      allocate(vmax(np))
   endif
   !
   if (store_maximum_flux) then
      allocate(qmax(np))
   endif
   !
   if (store_twet) then
       allocate(twet(np))
   endif
   !
   if (store_t_zsmax) then
       allocate(t_zsmax(np))
   endif
   !
   if (store_tsunami_arrival_time) then
      allocate(tsunami_arrival_time(np))
   endif
   !
   ! Set initial conditions (found in module sfincs_initial_conditions)
   !
   call set_initial_conditions()
   !
   if (wiggle_suppression) then 
      zs0 = zs
   endif
   !
   ! Loop through combined uv points and determine average uv and q (this also happens at the end of sfincs_momentum.f90).
   !
   do icuv = 1, ncuv
      !
      ! Average of the two uv points
      !
      q(cuv_index_uv(icuv))  = 0.5*(q(cuv_index_uv1(icuv)) + q(cuv_index_uv2(icuv)))
      uv(cuv_index_uv(icuv)) = 0.5*(uv(cuv_index_uv1(icuv)) + uv(cuv_index_uv2(icuv)))
      !
   enddo
   !
   if (wavemaker .or. snapwave) then !TL: zsm also used in sfincs_continuity if 'snapwave=true' > todo: check if needed, or only for wavemaker
      !      
      zsm = zs
      !
      if (store_maximum_waterlevel) then
         maxzsm = -999.0
      endif
      !
   endif   
   !
   if (store_maximum_waterlevel) then
      zsmax = -999.0
   endif
   !
   if (store_zvolume_max) then
      zvolmax = 0.0
   endif
   !
   if (store_maximum_velocity) then
      vmax = -999.0
   endif
   !
   if (store_twet) then
      twet = 0.0
   endif
   !
   if (store_t_zsmax) then
      t_zsmax = -999.0
   endif
   !
   ! NOTE: do NOT reset uv here. set_initial_conditions() (called above) computes the initial
   ! velocity uv = q/hu from the restart/initial flux. The flux-form (Bates) scheme carries q and
   ! recomputes uv every step, so this stray reset was harmless there, but the velocity-form
   ! (momentum_scheme=velocity) scheme carries uv0 and needs the initial velocity preserved.
   !
   if (wind) then
      !
      allocate(tauwu(np))
      allocate(tauwv(np))
      allocate(tauwu0(np))
      allocate(tauwv0(np))
      allocate(tauwu1(np))
      allocate(tauwv1(np))
      tauwu  = 0.0
      tauwv  = 0.0
      tauwu0 = 0.0
      tauwv0 = 0.0
      tauwu1 = 0.0
      tauwv1 = 0.0
      !
      if (store_wind) then
         !
         allocate(windu(np))
         allocate(windv(np))
         allocate(windu0(np))
         allocate(windv0(np))
         allocate(windu1(np))
         allocate(windv1(np))
         !
         windu  = 0.0
         windv  = 0.0
         windu0 = 0.0
         windv0 = 0.0
         windu1 = 0.0
         windv1 = 0.0
         !
         if (store_wind_max) then
            allocate(windmax(np))
            windmax = 0.0
         endif
         !
      endif   
      !
   endif
   !
   if (patmos) then
      !
      allocate(patm(np))
      allocate(patm0(np))
      allocate(patm1(np))
      patm    = 0.0
      patm0   = 0.0
      patm1   = 0.0
      if (pavbnd>0) then
         allocate(patmb(ngbnd))
         patmb = 101200.0
      endif
      !
   endif
   !
   if (precip) then
      !
      allocate(prcp(np))
      allocate(prcp0(np))
      allocate(prcp1(np))
      prcp    = 0.0
      prcp0   = 0.0
      prcp1   = 0.0
      !
      allocate(netprcp(np))      
      netprcp  = 0.0
      !
   endif
   !
   if (use_qext) then
      !
      ! Allocating and initializing qext (used e.g. for BMI coupling with ground water model)
      !
      allocate(qext(np))
      qext = 0.0
      !
   endif    
   !
   if (subgrid) then
      !
      ! Compute initial volumes
      !
      call compute_initial_subgrid_volumes()
      !
   endif
   !
   if (wave_enhanced_roughness) then
      !
      allocate(uorb(np))
      allocate(gnapp2(npuv))
      !
      uorb = 0.0      
      gnapp2 = 0.0
      !
   endif   
   !
   end subroutine



   subroutine deallocate_quadtree()
   !
   ! De-allocate arrays that are no longer necessary
   !
   ! Called from: sfincs_initialize (sfincs_lib)
   !
   use sfincs_data
   use quadtree
   !
   ! In sfincs_data
   !
   ! The following arrays are still needed for the timestep analysis (see sfincs_ncoutput.F90)
   !
   ! The single-face arrays z_index_uv_md1/mu1/nd1/nu1 are kept: with uv_flags_type they give
   ! the same-level neighbour of every cell, which the moving-bed depth filter
   ! (depth_filter_cell_field) walks at run time.
   !
   if (.not. timestep_analysis) then
      !
      if (allocated(z_index_uv_md2)) deallocate(z_index_uv_md2)
      if (allocated(z_index_uv_mu2)) deallocate(z_index_uv_mu2)
      if (allocated(z_index_uv_nd2)) deallocate(z_index_uv_nd2)
      if (allocated(z_index_uv_nu2)) deallocate(z_index_uv_nu2)
      !
   endif
   !
   ! In quadtree
   !
   if (allocated(quadtree_md)) deallocate(quadtree_md)
   if (allocated(quadtree_md1)) deallocate(quadtree_md1)
   if (allocated(quadtree_md2)) deallocate(quadtree_md2)
   if (allocated(quadtree_mu)) deallocate(quadtree_mu)
   if (allocated(quadtree_mu1)) deallocate(quadtree_mu1)
   if (allocated(quadtree_mu2)) deallocate(quadtree_mu2)
   if (allocated(quadtree_nd)) deallocate(quadtree_nd)
   if (allocated(quadtree_nd1)) deallocate(quadtree_nd1)
   if (allocated(quadtree_nd2)) deallocate(quadtree_nd2)
   if (allocated(quadtree_nu)) deallocate(quadtree_nu)
   if (allocated(quadtree_nu1)) deallocate(quadtree_nu1)
   if (allocated(quadtree_nu2)) deallocate(quadtree_nu2)
   if (allocated(quadtree_n)) deallocate(quadtree_n)
   if (allocated(quadtree_m)) deallocate(quadtree_m)
   if (allocated(quadtree_n_oddeven)) deallocate(quadtree_n_oddeven)
   if (allocated(quadtree_m_oddeven)) deallocate(quadtree_m_oddeven)
   if (allocated(quadtree_xz)) deallocate(quadtree_xz)
   if (allocated(quadtree_yz)) deallocate(quadtree_yz)
   if (allocated(quadtree_zz)) deallocate(quadtree_zz)
   if (allocated(quadtree_mask)) deallocate(quadtree_mask)
   if (allocated(quadtree_snapwave_mask)) deallocate(quadtree_snapwave_mask)
   !
   if (.not. bmi) then 
      !
      ! These may be used when finding grid cell when using BMI
      !
      if (allocated(quadtree_level)) deallocate(quadtree_level)
      if (allocated(quadtree_nm_indices)) deallocate(quadtree_nm_indices)
      if (allocated(quadtree_first_point_per_level)) deallocate(quadtree_first_point_per_level)  
      if (allocated(quadtree_last_point_per_level)) deallocate(quadtree_last_point_per_level)      
      if (allocated(quadtree_dxr)) deallocate(quadtree_dxr)
      if (allocated(quadtree_dyr)) deallocate(quadtree_dyr)
      !
   endif
   !
   end subroutine
   

end module
