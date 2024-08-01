#define NF90(nf90call) call handle_err(nf90call,__FILE__,__LINE__)
module sfincs_subgrid
   !
   use netcdf       
   !
   type net_type_subgrid
       integer :: ncid
       integer :: np_dimid
       integer :: npuv_dimid
       integer :: nlevels_dimid
       integer :: z_zmin_varid, z_zmax_varid, z_volmax_varid, z_dep_varid
       integer :: uv_zmin_varid, uv_zmax_varid, uv_fnfit_varid, uv_navg_w_varid, uv_nrep_varid, uv_havg_varid, uv_pwet_varid
   end type      
   type(net_type_subgrid) :: net_file_sbg
   !
contains

   subroutine read_subgrid_file()
   !
   use sfincs_data
   !
   implicit none
   !
   ! Check what sort of file we're dealing with
   !
   NF90(nf90_open(trim(sbgfile), NF90_CLOBBER, net_file_sbg%ncid))
   !
   if (net_file_sbg%ncid > 0) then
      !
      NF90(nf90_close(net_file_sbg%ncid))
      !
      ! Netcdf format with havg and nrep
      ! 
      huthresh = 0.0
      ! 
      call read_subgrid_file_netcdf()
      !
   else
      !
      ! Binary format with hrep and navg
      !
      write(*,*)'Warning : subgrid file has the "old" binary format, the simulation will continue, but we do recommended switching to the new Netcdf subgrid input format!'            ! 
      !
      call read_subgrid_file_original()
      !
   endif
   ! 
   end subroutine


   subroutine read_subgrid_file_netcdf()
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   real*4, dimension(:,:),   allocatable :: zbg
   real*4, dimension(:),     allocatable :: rtmp
   real*4, dimension(:),     allocatable :: rtmpz
   real*4, dimension(:),     allocatable :: rtmpuv
   real*4, dimension(:,:),   allocatable :: rtmpz2
   real*4, dimension(:,:),   allocatable :: rtmpuv2
   integer*4, dimension(:),  allocatable :: uv_index
   integer*4, dimension(:),  allocatable :: z_index
   !
   integer :: idummy
   integer :: ip
   integer :: ipuv
   integer :: ilevel
   integer :: nm
   integer :: nmu
   integer :: n
   integer :: m
   integer :: npzq
   integer :: npuvq
   integer :: npuvs
   integer :: np_nc
   integer :: npuv_nc
   !
   ! Read subgrid data
   !
   write(*,*)'Reading sub-grid netCDF file ...'
   !
   NF90(nf90_open(trim(sbgfile), NF90_CLOBBER, net_file_sbg%ncid))
   !          
   ! Using Quadtree file.
   ! This means that the subgrid file contains data for the entire quadtree. So also for points with kcs==0 !
   ! This also means that the data needs to be re-mapped to the active cell indices.
   !
   ! Get dimensions sizes    
   !
   NF90(nf90_inq_dimid(net_file_sbg%ncid, 'levels', net_file_sbg%nlevels_dimid))
   NF90(nf90_inq_dimid(net_file_sbg%ncid, 'np',   net_file_sbg%np_dimid))
   NF90(nf90_inq_dimid(net_file_sbg%ncid, 'npuv', net_file_sbg%npuv_dimid))
   !
   NF90(nf90_inquire_dimension(net_file_sbg%ncid, net_file_sbg%np_dimid,    len = np_nc)) ! nr of z points
   NF90(nf90_inquire_dimension(net_file_sbg%ncid, net_file_sbg%npuv_dimid,  len = npuv_nc)) ! nr of uv points
   NF90(nf90_inquire_dimension(net_file_sbg%ncid, net_file_sbg%nlevels_dimid, len = subgrid_nlevels)) ! nr of levels
   !
   ! Get variable id's
   !
   NF90(nf90_inq_varid(net_file_sbg%ncid, 'z_zmin',     net_file_sbg%z_zmin_varid))
   NF90(nf90_inq_varid(net_file_sbg%ncid, 'z_zmax',     net_file_sbg%z_zmax_varid))
   NF90(nf90_inq_varid(net_file_sbg%ncid, 'z_volmax',   net_file_sbg%z_volmax_varid))
   NF90(nf90_inq_varid(net_file_sbg%ncid, 'z_level',    net_file_sbg%z_dep_varid))
   NF90(nf90_inq_varid(net_file_sbg%ncid, 'uv_zmin',    net_file_sbg%uv_zmin_varid))
   NF90(nf90_inq_varid(net_file_sbg%ncid, 'uv_zmax',    net_file_sbg%uv_zmax_varid))
   NF90(nf90_inq_varid(net_file_sbg%ncid, 'uv_navg',    net_file_sbg%uv_navg_w_varid))
   NF90(nf90_inq_varid(net_file_sbg%ncid, 'uv_ffit',    net_file_sbg%uv_fnfit_varid))
   NF90(nf90_inq_varid(net_file_sbg%ncid, 'uv_havg',    net_file_sbg%uv_havg_varid))
   NF90(nf90_inq_varid(net_file_sbg%ncid, 'uv_nrep',    net_file_sbg%uv_nrep_varid))
   NF90(nf90_inq_varid(net_file_sbg%ncid, 'uv_pwet',    net_file_sbg%uv_pwet_varid))
   !
   ! ! Check if dimensions match
   ! !
   ! if (np /= np_nc) then
   !    write(*,'(a,i8,a,i8,a)')'Error! Number of cells in subgrid file ',np_nc,' does not match number of active cells in mesh ',np,' !'
   ! endif
   ! if (npuv /= npuv_nc) then
   !    write(*,'(a,i8,a,i8,a)')'Error! Number of velocity points in subgrid file ',npuv_nc,' does not match number of velocity points in mesh ',npuv,' !'
   ! endif
   !
   npzq  = np_nc
   npuvq = npuv_nc
   !
   allocate(subgrid_z_zmin(np))
   allocate(subgrid_z_zmax(np))
   allocate(subgrid_z_volmax(np))
   allocate(subgrid_z_dep(subgrid_nlevels, np))
   allocate(subgrid_uv_zmin(npuv))
   allocate(subgrid_uv_zmax(npuv))
   allocate(subgrid_uv_navg_w(npuv))
   allocate(subgrid_uv_fnfit(npuv))
   allocate(subgrid_uv_havg(subgrid_nlevels, npuv))
   allocate(subgrid_uv_nrep(subgrid_nlevels, npuv))
   allocate(subgrid_uv_pwet(subgrid_nlevels, npuv))
   allocate(subgrid_uv_havg_zmax(npuv))
   allocate(subgrid_uv_nrep_zmax(npuv))
   !
   allocate(rtmpz(npzq))
   allocate(rtmpuv(npuvq))
   allocate(rtmpz2(subgrid_nlevels, npzq))
   allocate(rtmpuv2(subgrid_nlevels, npuvq))
   allocate(uv_index(npuv))
   allocate(z_index(np))
   !
   write(*,*)'Number of subgrid levels : ',subgrid_nlevels
   !
   ! Need to make a new temporary re-mapping index array for the u and v points
   ! This is needed for reading in the subgrid file which has values for the entire quadtree grid
   !
   npuvq = 0
   npuvs = 0
   !
   if (use_quadtree) then
      !
      do ip = 1, quadtree_nr_points
         !
         nm = index_sfincs_in_quadtree(ip)
         !
         if (nm>0) then
            z_index(nm) = ip  
         endif 
         !
         ! Right
         !
         if (quadtree_mu(ip)<1) then
            if (quadtree_mu1(ip)>0) then
               npuvq = npuvq + 1
               if (nm>0) then
                  if (z_index_uv_mu1(nm)>0) then
                     npuvs = npuvs + 1
                     uv_index(npuvs) = npuvq
                  endif   
               endif
            endif   
         else
            if (quadtree_mu1(ip)>0) then
               npuvq = npuvq + 1
               if (nm>0) then
                  if (z_index_uv_mu1(nm)>0) then
                     npuvs = npuvs + 1
                     uv_index(npuvs) = npuvq
                  endif   
               endif   
            endif   
            if (quadtree_mu2(ip)>0) then
               npuvq = npuvq + 1
               if (nm>0) then
                  if (z_index_uv_mu2(nm)>0) then
                     npuvs = npuvs + 1
                     uv_index(npuvs) = npuvq
                  endif   
               endif   
            endif   
         endif   
         !
         ! Above
         !
         if (quadtree_nu(ip)<1) then
            if (quadtree_nu1(ip)>0) then
               npuvq = npuvq + 1
               if (nm>0) then
                  if (z_index_uv_nu1(nm)>0) then
                     npuvs = npuvs + 1
                     uv_index(npuvs) = npuvq
                  endif   
               endif   
            endif   
         else 
            if (quadtree_nu1(ip)>0) then
               npuvq = npuvq + 1
               if (nm>0) then
                  if (z_index_uv_nu1(nm)>0) then
                     npuvs = npuvs + 1
                     uv_index(npuvs) = npuvq
                  endif   
               endif   
            endif   
            if (quadtree_nu2(ip)>0) then
               npuvq = npuvq + 1
               if (nm>0) then
                  if (z_index_uv_nu2(nm)>0) then
                     npuvs = npuvs + 1
                     uv_index(npuvs) = npuvq
                  endif   
               endif   
            endif   
         endif   
         !
      enddo
      !
   else
      !
      ! Regular grid
      !
      ipuv = 0
      !  
      do nm = 1, np
         !
         z_index(nm) = nm
         !
         if (z_index_uv_mu1(nm)>0) then
            !
            ipuv = ipuv + 1
            ! 
            uv_index(ipuv) = ipuv
            !
         endif
         !
         if (z_index_uv_nu1(nm)>0) then
            !
            ipuv = ipuv + 1
            ! 
            uv_index(ipuv) = ipuv
            !
         endif
         !
      enddo   
      !
   endif
   !
   ! Read Z points
   !
   NF90(nf90_get_var(net_file_sbg%ncid, net_file_sbg%z_zmin_varid, rtmpz(:)))
   !
   do ip = 1, np
      subgrid_z_zmin(ip) = rtmpz(z_index(ip))
   enddo   
   !
   NF90(nf90_get_var(net_file_sbg%ncid, net_file_sbg%z_zmax_varid, rtmpz(:)))
   !
   do ip = 1, np
      subgrid_z_zmax(ip) = rtmpz(z_index(ip))
   enddo   
   !
   NF90(nf90_get_var(net_file_sbg%ncid, net_file_sbg%z_volmax_varid, rtmpz(:)))
   !
   do ip = 1, np
      subgrid_z_volmax(ip) = rtmpz(z_index(ip))
   enddo   
   !
   NF90(nf90_get_var(net_file_sbg%ncid, net_file_sbg%z_dep_varid, rtmpz2(:,:)))
   !
   do ip = 1, np
      do ilevel = 1, subgrid_nlevels
         subgrid_z_dep(ilevel, ip) = rtmpz2(ilevel, z_index(ip))
      enddo   
   enddo
   !
   ! Read UV points
   !
   NF90(nf90_get_var(net_file_sbg%ncid, net_file_sbg%uv_zmin_varid, rtmpuv(:) ))
   !
   do ip = 1, npuv
      subgrid_uv_zmin(ip) = rtmpuv(uv_index(ip))
   enddo   
   !
   NF90(nf90_get_var(net_file_sbg%ncid, net_file_sbg%uv_zmax_varid, rtmpuv(:) ))
   !
   do ip = 1, npuv
      subgrid_uv_zmax(ip) = rtmpuv(uv_index(ip))
   enddo   
   !
   do ip = 1, npuv
      subgrid_uv_zmax(ip) = max(subgrid_uv_zmax(ip), subgrid_uv_zmin(ip) + 0.01)
   enddo   
   !
   NF90(nf90_get_var(net_file_sbg%ncid, net_file_sbg%uv_fnfit_varid, rtmpuv(:) ))
   !
   do ip = 1, npuv
      subgrid_uv_fnfit(ip) = rtmpuv(uv_index(ip))
   enddo   
   !
   NF90(nf90_get_var(net_file_sbg%ncid, net_file_sbg%uv_navg_w_varid, rtmpuv(:) ))
   !
   do ip = 1, npuv
      subgrid_uv_navg_w(ip) = g * max(rtmpuv(uv_index(ip)), 0.0001) ** 2
   enddo   
   !
   NF90(nf90_get_var(net_file_sbg%ncid, net_file_sbg%uv_havg_varid, rtmpuv2(:,:) ))
   !
   do ip = 1, npuv
      do ilevel = 1, subgrid_nlevels
         subgrid_uv_havg(ilevel, ip) = rtmpuv2(ilevel, uv_index(ip))
      enddo   
   enddo
   !
   NF90(nf90_get_var(net_file_sbg%ncid, net_file_sbg%uv_nrep_varid, rtmpuv2(:,:) ))
   !
   do ip = 1, npuv
      do ilevel = 1, subgrid_nlevels
         !
         ! Already convert here to gn^2
         !
         subgrid_uv_nrep(ilevel, ip) = g*max(rtmpuv2(ilevel, uv_index(ip)), 0.0001)**2
         !
      enddo   
   enddo
   !
   NF90(nf90_get_var(net_file_sbg%ncid, net_file_sbg%uv_pwet_varid, rtmpuv2(:,:) ))
   !
   do ip = 1, npuv
      do ilevel = 1, subgrid_nlevels
         subgrid_uv_pwet(ilevel, ip) = rtmpuv2(ilevel, uv_index(ip))
      enddo   
   enddo
   !
   NF90(nf90_close(net_file_sbg%ncid))
   !
   ! Make sure that zmin at uv point is always higher or equal to zmin of neighboring cells (this is always the case if the pre-processing was done right)
   !
   do ip = 1, npuv
      nm = uv_index_z_nm(ip)
      nmu = uv_index_z_nmu(ip)
      subgrid_uv_zmin(ip) = max(subgrid_uv_zmin(ip), subgrid_z_zmin(nm))
      subgrid_uv_zmin(ip) = max(subgrid_uv_zmin(ip), subgrid_z_zmin(nmu))
   enddo   
   !
   ! Make sure zmax is always bigger than zmin
   !
   do nm = 1, np
      if (subgrid_z_zmax(nm) - subgrid_z_zmin(nm) < 0.001) subgrid_z_zmax(nm) = subgrid_z_zmax(nm) + 0.001
   enddo
   !
   do nm = 1, npuv
      if (subgrid_uv_zmax(nm) - subgrid_uv_zmin(nm) < 0.001) subgrid_uv_zmax(nm) = subgrid_uv_zmax(nm) + 0.001
   enddo
   !
   ! Make arrays for subgrid_uv_havg_zmax and subgrid_uv_nrep_zmax for faster searching
   !
   do ip = 1, npuv
      subgrid_uv_havg_zmax(ip) = subgrid_uv_havg(subgrid_nlevels, ip) - subgrid_uv_zmax(ip)
      subgrid_uv_nrep_zmax(ip) = subgrid_uv_nrep(subgrid_nlevels, ip)
   enddo    
   ! 
   deallocate(rtmpz)
   deallocate(rtmpz2)
   deallocate(rtmpuv)
   deallocate(rtmpuv2)
   deallocate(uv_index)
   deallocate(z_index)
   !
   end subroutine


   subroutine read_subgrid_file_original()
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   real*4, dimension(:,:),   allocatable :: zbg
   real*4, dimension(:),     allocatable :: rtmp
   real*4, dimension(:),     allocatable :: rtmpz
   real*4, dimension(:),     allocatable :: rtmpuv
   integer*4, dimension(:),  allocatable :: uv_index_qt_in_sf
   !
   integer :: idummy
   integer :: ioption
   integer :: ip
   integer :: ilevel
   integer :: nm
   integer :: nmu
   integer :: n
   integer :: m
   integer :: npzq
   integer :: npuvq
   integer :: npuvs
   !
   ! Read subgrid data
   !
   if (use_quadtree) then
      !
      ! Using Quadtree file.
      ! This means that the subgrid file contains data for the entire quadtree. So also for points with kcs==0 !
      ! This also means that the data needs to be re-mapped to the active cell indices.
      !
      write(*,*)'Reading ',trim(sbgfile), ' ...'
      open(unit = 500, file = trim(sbgfile), form = 'unformatted', access = 'stream')
      read(500)idummy ! version
      read(500)npzq ! nr cells
      read(500)npuvq ! nr uv points
      read(500)subgrid_nlevels
      ! In the 'old' binary format, the number of bins was written to the sbg file. Add 1 to get to nr of levels.
      subgrid_nlevels = subgrid_nlevels + 1
      allocate(subgrid_z_zmin(np))
      allocate(subgrid_z_zmax(np))
      allocate(subgrid_z_volmax(np))
      allocate(subgrid_z_dep(subgrid_nlevels, np))
      allocate(subgrid_uv_zmin(npuv))
      allocate(subgrid_uv_zmax(npuv))
      allocate(subgrid_uv_havg(subgrid_nlevels, npuv))
      allocate(subgrid_uv_nrep(subgrid_nlevels, npuv))
      allocate(subgrid_uv_havg_zmax(npuv))
      allocate(subgrid_uv_nrep_zmax(npuv))
      allocate(subgrid_uv_navg_w(npuv))
      allocate(subgrid_uv_fnfit(npuv))
      allocate(subgrid_uv_pwet(subgrid_nlevels, npuv))
      subgrid_uv_fnfit = 0.0
      subgrid_uv_pwet  = 1.0
      !
      allocate(rtmpz(npzq))
      allocate(rtmpuv(npuvq))
      allocate(uv_index_qt_in_sf(npuv))
      !
      write(*,*)'Number of subgrid levels : ',subgrid_nlevels
      !
      ! Need to make a new temporary re-mapping index array for the u and v points
      ! This is needed for reading in the subgrid file which has values for the entire quadtree grid
      !
      npuvq = 0
      npuvs = 0
      !
      do ip = 1, quadtree_nr_points
         !
         nm = index_sfincs_in_quadtree(ip)
         !
         ! Right
         !
         if (quadtree_mu(ip)<1) then
            if (quadtree_mu1(ip)>0) then
               npuvq = npuvq + 1
               if (nm>0) then
                  if (z_index_uv_mu1(nm)>0) then
                     npuvs = npuvs + 1
                     uv_index_qt_in_sf(npuvs) = npuvq
                  endif   
               endif
            endif   
         else
            if (quadtree_mu1(ip)>0) then
               npuvq = npuvq + 1
               if (nm>0) then
                  if (z_index_uv_mu1(nm)>0) then
                     npuvs = npuvs + 1
                     uv_index_qt_in_sf(npuvs) = npuvq
                  endif   
               endif   
            endif   
            if (quadtree_mu2(ip)>0) then
               npuvq = npuvq + 1
               if (nm>0) then
                  if (z_index_uv_mu2(nm)>0) then
                     npuvs = npuvs + 1
                     uv_index_qt_in_sf(npuvs) = npuvq
                  endif   
               endif   
            endif   
         endif   
         !
         ! Above
         !
         if (quadtree_nu(ip)<1) then
            if (quadtree_nu1(ip)>0) then
               npuvq = npuvq + 1
               if (nm>0) then
                  if (z_index_uv_nu1(nm)>0) then
                     npuvs = npuvs + 1
                     uv_index_qt_in_sf(npuvs) = npuvq
                  endif   
               endif   
            endif   
         else
            if (quadtree_nu1(ip)>0) then
               npuvq = npuvq + 1
               if (nm>0) then
                  if (z_index_uv_nu1(nm)>0) then
                     npuvs = npuvs + 1
                     uv_index_qt_in_sf(npuvs) = npuvq
                  endif   
               endif   
            endif   
            if (quadtree_nu2(ip)>0) then
               npuvq = npuvq + 1
               if (nm>0) then
                  if (z_index_uv_nu2(nm)>0) then
                     npuvs = npuvs + 1
                     uv_index_qt_in_sf(npuvs) = npuvq
                  endif   
               endif   
            endif   
         endif   
         !
      enddo
      !
      ! Read Z points
      !
      read(500)rtmpz
      do nm = 1, np
         subgrid_z_zmin(nm) = rtmpz(index_quadtree_in_sfincs(nm))
      enddo   
      !
      read(500)rtmpz
      do nm = 1, np
         subgrid_z_zmax(nm) = rtmpz(index_quadtree_in_sfincs(nm))
      enddo   
      !
      read(500)rtmpz
      ! do nm = 1, np
      !    subgrid_z_zmean(nm) = rtmpz(index_quadtree_in_sfincs(nm))
      ! enddo   
      !
      read(500)rtmpz
      do nm = 1, np
         subgrid_z_volmax(nm) = rtmpz(index_quadtree_in_sfincs(nm))
      enddo   
      !
      subgrid_z_dep(1, :) = max(subgrid_z_zmin(:), -20.0)
      do ilevel = 1, subgrid_nlevels - 1
         read(500)rtmpz
         do nm = 1, np
            subgrid_z_dep(ilevel + 1, nm) = rtmpz(index_quadtree_in_sfincs(nm))
         enddo   
      enddo
      !
      ! Read UV points
      !
      read(500)rtmpuv
      do nm = 1, npuv
         subgrid_uv_zmin(nm) = rtmpuv(uv_index_qt_in_sf(nm))
      enddo   
      !
      read(500)rtmpuv
      do nm = 1, npuv
         subgrid_uv_zmax(nm) = rtmpuv(uv_index_qt_in_sf(nm))
      enddo   
      !
      ! Initialize subgrid_uv_havg with huthresh
      !
      subgrid_uv_havg = huthresh
      do ilevel = 1, subgrid_nlevels - 1
         read(500)rtmpuv
         do nm = 1, npuv
            subgrid_uv_havg(ilevel + 1, nm) = max(rtmpuv(uv_index_qt_in_sf(nm)), huthresh)
         enddo   
      enddo
      !
      do ilevel = 1, subgrid_nlevels - 1
         read(500)rtmpuv
         do nm = 1, npuv
            ! Already convert here to gn^2
            subgrid_uv_nrep(ilevel + 1, nm) = g*max(rtmpuv(uv_index_qt_in_sf(nm)), 0.005)**2
         enddo   
      enddo
      ! Set bottom navg equal to one level above
      subgrid_uv_nrep(1, :) = subgrid_uv_nrep(2, :)
      !
      close(500)
      !
      deallocate(rtmpz)
      deallocate(rtmpuv)
      deallocate(uv_index_qt_in_sf)
      !
   else
      !
      ! Subgrid file on regular grid
      !
      write(*,*)'Reading ', trim(sbgfile), ' ...'
      open(unit = 500, file = trim(sbgfile), form = 'unformatted', access = 'stream')
      read(500)idummy ! nr cells
      read(500)ioption ! option
      read(500)subgrid_nlevels
      subgrid_nlevels = subgrid_nlevels + 1
      write(*,*)'Number of subgrid levels : ',subgrid_nlevels - 1
      allocate(subgrid_z_zmin(np))
      allocate(subgrid_z_zmax(np))
      allocate(subgrid_z_volmax(np))
      allocate(subgrid_z_dep(subgrid_nlevels, np))
      allocate(subgrid_uv_zmin(npuv))
      allocate(subgrid_uv_zmax(npuv))
      allocate(subgrid_uv_havg(subgrid_nlevels, npuv))
      allocate(subgrid_uv_nrep(subgrid_nlevels, npuv))
      allocate(subgrid_uv_havg_zmax(npuv))
      allocate(subgrid_uv_nrep_zmax(npuv))
      allocate(subgrid_uv_navg_w(npuv))
      allocate(subgrid_uv_fnfit(npuv))
      allocate(subgrid_uv_pwet(subgrid_nlevels, npuv))
      subgrid_uv_fnfit = 0.0
      subgrid_uv_pwet  = 1.0
      !
      allocate(rtmpz(np))
      !
      read(500)subgrid_z_zmin
      read(500)subgrid_z_zmax
      read(500)subgrid_z_volmax
      !
      ! ! Make sure zmax is always bigger than zmin (do this later on)
      ! !
      ! do nm = 1, np
      !    if (subgrid_z_zmax(nm) - subgrid_z_zmin(nm) < 0.01) subgrid_z_zmax(nm) = subgrid_z_zmax(nm) + 0.001
      ! enddo
      !
      subgrid_z_dep(1,:) = max(subgrid_z_zmin(:), -20.0)
      do ilevel = 1, subgrid_nlevels - 1
         read(500)rtmpz
         do nm = 1, np
            subgrid_z_dep(ilevel + 1, nm) = rtmpz(nm)
         enddo
      enddo   
      !
      ! U zmin
      !
      read(500) rtmpz
      !
      do nm = 1, np
         !
         ip = z_index_uv_mu1(nm) ! index of uv point
         !
         if (ip>0) then
            !
            subgrid_uv_zmin(ip) = rtmpz(nm)
            !
         endif
         !
      enddo   
      !         
      ! U zmax
      !
      read(500)rtmpz
      !
      do nm = 1, np
         !
         ip = z_index_uv_mu1(nm) ! index of uv point
         !
         if (ip>0) then
            !
            subgrid_uv_zmax(ip) = rtmpz(nm)
            !
         endif
         !
      enddo   
      !         
      ! U dhdz (not used anymore in this SFINCS version)
      !
      read(500)rtmpz
      !
      ! U hrep
      !
      subgrid_uv_havg(1,:) = huthresh
      do ilevel = 1, subgrid_nlevels - 1
         read(500)rtmpz
         do nm = 1, np
            ip = z_index_uv_mu1(nm) ! index of uv point
            if (ip>0) then
               subgrid_uv_havg(ilevel + 1, ip) = rtmpz(nm)
            endif
         enddo
      enddo         
      !
      ! U navg
      !
      do ilevel = 1, subgrid_nlevels - 1
         !
         read(500)rtmpz
         !            
         do nm = 1, np
            !
            ip = z_index_uv_mu1(nm) ! index of uv point
            !
            if (ip>0) then
               !
               subgrid_uv_nrep(ilevel + 1, ip) = g*max(rtmpz(nm), 0.0001)**2
               !
            endif
            !
         enddo
         !
      enddo
      !
      do nm = 1, np
         subgrid_uv_nrep(1, nm) = subgrid_uv_nrep(2, nm)
      enddo 
      !
      ! V zmin
      !
      read(500)rtmpz
      !
      do nm = 1, np
         !
         ip = z_index_uv_nu1(nm) ! index of uv point
         !
         if (ip>0) then
            !
            subgrid_uv_zmin(ip) = rtmpz(nm)
            !
         endif
         !
      enddo   
      !         
      ! V zmax
      !
      read(500)rtmpz
      !
      do nm = 1, np
         !
         ip = z_index_uv_nu1(nm) ! index of uv point
         !
         if (ip>0) then
            !
            subgrid_uv_zmax(ip) = rtmpz(nm)
            !
         endif
         !
      enddo   
      !         
      ! V dhdz (not used anymore in this SFINCS version)
      !
      read(500)rtmpz
      !
      ! V hrep
      !
      subgrid_uv_havg(1,:) = huthresh
      do ilevel = 1, subgrid_nlevels - 1
         read(500)rtmpz
         do nm = 1, np
            ip = z_index_uv_nu1(nm) ! index of uv point
            if (ip>0) then
               subgrid_uv_havg(ilevel + 1, ip) = rtmpz(nm)
            endif
         enddo
      enddo         
      !
      ! V navg
      !
      do ilevel = 1, subgrid_nlevels - 1
         !
         read(500)rtmpz
         !            
         do nm = 1, np
            !
            ip = z_index_uv_nu1(nm) ! index of uv point
            !
            if (ip>0) then
               !
               subgrid_uv_nrep(ilevel + 1, ip) = g*max(rtmpz(nm), 0.0001)**2
               !
            endif
            !
         enddo
         !
      enddo
      !
      do nm = 1, np
         subgrid_uv_nrep(1, nm) = subgrid_uv_nrep(2, nm)
      enddo 
      !
      do nm = 1, np
         subgrid_uv_nrep(1, nm) = subgrid_uv_nrep(2, nm)
      enddo
      !
      close(500)
      !
      deallocate(rtmpz)         
      !
   endif   
   !
   ! Make sure that zmin at uv point is always higher or equal to zmin of neighboring cells
   !
   do ip = 1, npuv
      nm = uv_index_z_nm(ip)
      nmu = uv_index_z_nmu(ip)
      subgrid_uv_zmin(ip) = max(subgrid_uv_zmin(ip), subgrid_z_zmin(nm))
      subgrid_uv_zmin(ip) = max(subgrid_uv_zmin(ip), subgrid_z_zmin(nmu))
   enddo   
   !
   ! Make sure zmax is always bigger than zmin
   !
   do nm = 1, np
      if (subgrid_z_zmax(nm) - subgrid_z_zmin(nm) < 0.001) subgrid_z_zmax(nm) = subgrid_z_zmax(nm) + 0.001
   enddo
   !
   do nm = 1, npuv
      if (subgrid_uv_zmax(nm) - subgrid_uv_zmin(nm) < 0.001) subgrid_uv_zmax(nm) = subgrid_uv_zmax(nm) + 0.001
   enddo
   !
   ! Make arrays for subgrid_uv_havg_zmax and subgrid_uv_nrep_zmax for faster searching
   !
   do ip = 1, npuv
      subgrid_uv_havg_zmax(ip) = subgrid_uv_havg(subgrid_nlevels, ip) - subgrid_uv_zmax(ip)
      subgrid_uv_nrep_zmax(ip) = subgrid_uv_nrep(subgrid_nlevels, ip)
      subgrid_uv_navg_w(ip)    = subgrid_uv_nrep_zmax(ip)
   enddo    
   !
   end subroutine



   subroutine compute_initial_subgrid_volumes()
   !
   use sfincs_data
   !
   implicit none
   !
   integer    :: nm, m, n, ivol, ilevel, ip
   real*4     :: dzvol
   real*4     :: facint
   real*4     :: one_minus_facint 
   !
   ! Compute initial volumes
   !
   allocate(z_volume(np))
   !
   do nm = 1, np
      !
      if (zs(nm)>=subgrid_z_zmax(nm)) then
         !
         ! Entire cell is wet, no interpolation from table needed
         !
         if (crsgeo) then
            z_volume(nm) = subgrid_z_volmax(nm) + cell_area_m2(nm)*(zs(nm) - max(subgrid_z_zmax(nm), -20.0))
         else   
            z_volume(nm) = subgrid_z_volmax(nm) + cell_area(z_flags_iref(nm))*(zs(nm) - max(subgrid_z_zmax(nm), -20.0))
         endif
         !
      else   
         !
         ! Interpolation required
         !
         ivol = 1
         do ilevel = 2, subgrid_nlevels
            if (subgrid_z_dep(ilevel, nm)>zs(nm)) then
               ivol = ilevel - 1
               exit
            endif
         enddo
         !
         dzvol  = subgrid_z_volmax(nm) / (subgrid_nlevels - 1)
         facint = (zs(nm) - subgrid_z_dep(ivol, nm)) / max(subgrid_z_dep(ivol + 1, nm) - subgrid_z_dep(ivol, nm), 0.001)
         z_volume(nm) = (ivol - 1)*dzvol + facint*dzvol
         !
      endif
      !
      if (use_storage_volume) then
         !
         ! Make sure wet cells do not have initial storage
         !
         storage_volume(nm) = max(storage_volume(nm) - z_volume(nm), 0.0)
         !
      endif
      !
   enddo
   !
   end subroutine

   subroutine handle_err(status,file,line)
      !
      integer, intent ( in)    :: status
      character(*), intent(in) :: file
      integer, intent ( in)    :: line
      integer :: status2
      !   
      if (status /= nf90_noerr) then
         !
         ! Do not give error message if this is not a valid netcdf file (it's possible that we tried to open a file with the 'old' binary format)
         !  
         if (status /= -51) then
            write(0,'("NETCDF ERROR: ",a,i6,":",a)') file, line, trim(nf90_strerror(status) )
         endif 
         !  
      end if
      !  
   end subroutine handle_err

end module
