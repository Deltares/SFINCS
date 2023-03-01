module sfincs_domain

contains

   subroutine initialize_domain()
   !
   ! Initialize SFINCS domain (indices, flags and depths)
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   integer ip, n, m, nm, nmu, num, ib, iref, npu, npv, ipuv, ipu, ipv, iud, iuu, ndm, nmd, j, iq, nmx, ibin, ind
   integer ikcuv2
   integer npac, idummy
   integer npuvq, npuvs, npzq
   !
   ! Temporary arrays
   !
   integer*1, dimension(:,:), allocatable :: kcsg
   !
   real*4, dimension(:,:),   allocatable :: zbg
   real*4, dimension(:),     allocatable :: rtmp
   real*4, dimension(:),     allocatable :: rtmpz
   real*4, dimension(:),     allocatable :: rtmpuv
   real*4, dimension(:),     allocatable :: rghfield
   !
   integer*1, dimension(:),   allocatable :: z_index_md
   integer*4, dimension(:),   allocatable :: z_index_md1
   integer*4, dimension(:),   allocatable :: z_index_md2
   integer*1, dimension(:),   allocatable :: z_index_mu
   integer*4, dimension(:),   allocatable :: z_index_mu1
   integer*4, dimension(:),   allocatable :: z_index_mu2
   integer*1, dimension(:),   allocatable :: z_index_nd
   integer*4, dimension(:),   allocatable :: z_index_nd1
   integer*4, dimension(:),   allocatable :: z_index_nd2
   integer*1, dimension(:),   allocatable :: z_index_nu
   integer*4, dimension(:),   allocatable :: z_index_nu1
   integer*4, dimension(:),   allocatable :: z_index_nu2
   !
   integer*1,          dimension(:),   allocatable :: msk
   !   
   integer*1,          dimension(:),   allocatable :: u_iref
   integer*1,          dimension(:),   allocatable :: u_type
   integer*4,          dimension(:),   allocatable :: u_index_nm
   integer*4,          dimension(:),   allocatable :: u_index_nmu
   integer*1,          dimension(:),   allocatable :: v_iref
   integer*1,          dimension(:),   allocatable :: v_type
   integer*4,          dimension(:),   allocatable :: v_index_nm
   integer*4,          dimension(:),   allocatable :: v_index_num
   !
   integer*4, dimension(:),   allocatable :: uv_index_qt_in_sf
   !
   ! For 'old' input
   !
   integer*4,          dimension(:),   allocatable :: indices
   !
   real*4 :: ylat
   real*4 :: dxymin
   !
   ! Set some flags
   !
   use_quadtree = .false.
   if (qtrfile(1:4) /= 'none') then
      use_quadtree = .true.
      write(*,*)'Info : Preparing SFINCS grid on quadtree mesh ...'   
   else
       write(*,*)'Info : Preparing SFINCS grid on regular mesh ...'      
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
   include_boundaries = .false.   
   !
   if (spwfile(1:4) /= 'none' .or. wndfile(1:4) /= 'none' .or. amufile(1:4) /= 'none' .or. netamuamvfile(1:4) /= 'none') then
      !
      wind = .true.
      write(*,*)'Turning on process: Wind'
      !
      if (spw_precip) then
          precip = .true.
          write(*,*)'Turning on process: Precipitation'
      endif
   endif
   !
   if ((ampfile(1:4) /= 'none' .or. netampfile(1:4) /= 'none' .or. spwfile(1:4) /= 'none') .and. baro==1) then
       patmos = .true.
       write(*,*)'Turning on process: Atmospheric pressure'       
   endif
   !
   if (prcpfile(1:4) /= 'none' .or. amprfile(1:4) /= 'none' .or. netamprfile(1:4) /= 'none' .or. spw_precip) then
       precip = .true.
       write(*,*)'Turning on process: Precipitation'       
   endif
   !
   ! Check for time-spatially varying meteo
   !
   if (spwfile(1:4) /= 'none' .or. amufile(1:4) /= 'none' .or. amprfile(1:4) /= 'none' .or. ampfile(1:4) /= 'none' .or. netampfile(1:4) /= 'none' .or. netamprfile(1:4) /= 'none' .or. netamuamvfile(1:4) /= 'none') then
      meteo3d = .true.
      write(*,*)'Turning on flag: meteo3d'
   endif
   !
   wind_reduction = .false.
   if (spwfile(1:4) /= 'none') then
      if (z0lfile(1:4) /= 'none') then
         wind_reduction = .true.
         write(*,*)'Turning on process: wind reduction over land'
      endif
   endif            
   !
   ! If not meteo3d, then also no storage of meteo data in netcdf output
   !
   if (.not. meteo3d) then
      store_meteo = .false.
   endif
   !
   ! READ MESH
   !
   ! 2 options:
   !
   ! 1) quadtree grid
   ! 2) 'old' inputs with index file
   !
   if (use_quadtree) then
      ! 
      ! Read quadtree file
      !
      write(*,*)'Reading quadtree file ', trim(qtrfile), ' ...'
      !
      call quadtree_read_file(qtrfile)
      !
   else
      !
      if (inputtype=='asc') then
         !
         ! Read ASCII input
         !
         write(*,*)'Reading ASCII input ...'
         !
         allocate(kcsg(nmax,  mmax))       ! KCS
         !
         kcsg = 0      
         !
         ! Read mask file
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
         write(*,*)'Reading ', trim(indexfile), ' ...'
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
   ! Read mask file (must have same number of cells as quadtree file !)
   !
   write(*,*)'Reading ', trim(mskfile), ' ...'
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
      open(unit = 500, file = trim(mskfile), form = 'unformatted', access = 'stream')
      read(500)msk
      close(500)
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
   allocate(z_flags_type(np))
   !
   z_flags_type = 0
   z_flags_iref = 0
   !      
   ! Temporary index arrays
   !
   allocate(z_index_md(np))
   allocate(z_index_mu(np))
   allocate(z_index_nd(np))
   allocate(z_index_nu(np))
   allocate(z_index_md1(np))
   allocate(z_index_md2(np))
   allocate(z_index_mu1(np))
   allocate(z_index_mu2(np))
   allocate(z_index_nd1(np))
   allocate(z_index_nd2(np))
   allocate(z_index_nu1(np))
   allocate(z_index_nu2(np))
   !
   allocate(z_index_z_n(np))
   allocate(z_index_z_m(np))
   !
   z_index_md    = 0
   z_index_mu    = 0
   z_index_nd    = 0
   z_index_nu    = 0
   z_index_md1   = 0
   z_index_md2   = 0
   z_index_mu1   = 0
   z_index_mu2   = 0
   z_index_nd1   = 0
   z_index_nd2   = 0
   z_index_nu1   = 0
   z_index_nu2   = 0
   z_index_z_n   = 0
   z_index_z_m   = 0
!   z_index_z_nm  = 0
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
   z_index_uv_md1   = 0
   z_index_uv_md2   = 0
   z_index_uv_mu1   = 0
   z_index_uv_mu2   = 0
   z_index_uv_nd1   = 0
   z_index_uv_nd2   = 0
   z_index_uv_nu1   = 0
   z_index_uv_nu2   = 0
   z_index_z_n      = 0
   z_index_z_m      = 0
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
      z_index_nu(ip)    = quadtree_nu(iq)
      z_index_nd(ip)    = quadtree_nd(iq)
      z_index_mu(ip)    = quadtree_mu(iq)
      z_index_md(ip)    = quadtree_md(iq)
      !
      ! Neighbors
      !
      if (quadtree_nu1(iq)>0) then
         z_index_nu1(ip)   = index_sfincs_in_quadtree(quadtree_nu1(iq))
      endif
      if (quadtree_nu2(iq)>0) then
         z_index_nu2(ip)   = index_sfincs_in_quadtree(quadtree_nu2(iq))
      endif
      if (quadtree_nd1(iq)>0) then
         z_index_nd1(ip)   = index_sfincs_in_quadtree(quadtree_nd1(iq))
      endif
      if (quadtree_nd2(iq)>0) then
         z_index_nd2(ip)   = index_sfincs_in_quadtree(quadtree_nd2(iq))
      endif
      if (quadtree_mu1(iq)>0) then
         z_index_mu1(ip)   = index_sfincs_in_quadtree(quadtree_mu1(iq))
      endif
      if (quadtree_mu2(iq)>0) then
         z_index_mu2(ip)   = index_sfincs_in_quadtree(quadtree_mu2(iq))
      endif
      if (quadtree_md1(iq)>0) then
         z_index_md1(ip)   = index_sfincs_in_quadtree(quadtree_md1(iq))
      endif
      if (quadtree_md2(iq)>0) then
         z_index_md2(ip)   = index_sfincs_in_quadtree(quadtree_md2(iq))
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
      if (z_index_mu(nm)<1) then
         !
         if (z_index_mu1(nm)>0) then
            npu = npu + 1
         endif   
         !
      else
         !
         if (z_index_mu1(nm)>0) then
            npu = npu + 1
         endif   
         !
         if (z_index_mu2(nm)>0) then
            npu = npu + 1
         endif   
         !            
      endif   
      !
      ! Above
      !
      if (z_index_nu(nm)<1) then
         !
         if (z_index_nu1(nm)>0) then
            npv = npv + 1
         endif   
         !
      else
         !
         if (z_index_nu1(nm)>0) then
            npv = npv + 1
         endif   
         !
         if (z_index_nu2(nm)>0) then
            npv = npv + 1
         endif   
         !            
      endif   
      !
   enddo 
   !
   npuv = npu + npv
   !
   ! UV-points
   !
   allocate(kcuv(npuv))
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
!   if (advection .or. coriolis) then
      allocate(uv_flags_adv(npuv))
      uv_flags_adv = 0
!   endif   
   !
!   if (viscosity) then
      allocate(uv_flags_vis(npuv))
      uv_flags_vis = 0
!   endif   
   !      
!   if (advection) then
      allocate(uv_index_u_nmd(npuv))
      allocate(uv_index_u_nmu(npuv))
      allocate(uv_index_u_num(npuv))
      allocate(uv_index_u_ndm(npuv))
      allocate(uv_index_v_ndm(npuv))
      allocate(uv_index_v_nm(npuv))
      allocate(uv_index_v_ndmu(npuv))
      allocate(uv_index_v_nmu(npuv))
      uv_index_u_nmd  = 0
      uv_index_u_nmu  = 0
      uv_index_u_num  = 0
      uv_index_u_ndm  = 0
      uv_index_v_ndm  = 0
      uv_index_v_nm   = 0
      uv_index_v_ndmu = 0
      uv_index_v_nmu  = 0
!   endif
   !
!   if (viscosity .and. .not.advection) then
!      allocate(uv_index_u_nmd(npuv))
!      allocate(uv_index_u_nmu(npuv))
!      allocate(uv_index_u_num(npuv))
!      allocate(uv_index_u_ndm(npuv))
!      uv_index_u_nmd  = 0
!      uv_index_u_nmu  = 0
!      uv_index_u_num  = 0
!      uv_index_u_ndm  = 0
!   endif
   !
!   if (coriolis .and. .not.advection) then
!      allocate(uv_index_v_nm(npuv))
!      uv_index_v_nm   = 0
!   endif
   !
   ip = 0
   !
   do nm = 1, np
      !
      ! Right
      !
      if (z_index_mu(nm)==0) then
         !
         ! Same level
         !
         if (z_index_mu1(nm)>0) then
            !
            ip                                 = ip + 1
            uv_index_z_nm(ip)                  = nm
            uv_index_z_nmu(ip)                 = z_index_mu1(nm)
            z_index_uv_mu1(nm)                 = ip 
            z_index_uv_md1(z_index_mu1(nm))    = ip
            uv_flags_iref(ip)                  = z_flags_iref(nm) 
            uv_flags_dir(ip)                   = 0 ! u point
            uv_flags_type(ip)                  = 0 ! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            uv_flags_vis(ip)                   = 1 ! enable viscosity
            if (advection .or. coriolis) uv_flags_adv(ip) = 1 ! enable advection
            !
         endif   
         !
      elseif (z_index_mu(nm)==-1) then
         !
         ! Coarser to the right
         !
         if (z_index_mu1(nm)>0) then
            !
            ip                                 = ip + 1
            uv_index_z_nm(ip)                  = nm
            uv_index_z_nmu(ip)                 = z_index_mu1(nm)
            z_index_uv_mu1(nm)                 = ip 
            if (z_index_md1(z_index_mu1(nm)) == nm) then
                z_index_uv_md1(z_index_mu1(nm)) = ip
            else
                z_index_uv_md2(z_index_mu1(nm)) = ip
            endif   
            uv_flags_iref(ip)                  = z_flags_iref(nm) 
            uv_flags_dir(ip)                   = 0 ! u point
            uv_flags_type(ip)                  = -1! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            uv_flags_vis(ip)                   = 1 ! enable viscosity
            if (advection .or. coriolis) uv_flags_adv(ip) = 1 ! enable advection
            !
         endif   
         !
      else
         !
         ! Finer to the right
         !
         if (z_index_mu1(nm)>0) then
            !
            ip                                 = ip + 1
            uv_index_z_nm(ip)                  = nm
            uv_index_z_nmu(ip)                 = z_index_mu1(nm)
            z_index_uv_mu1(nm)                 = ip 
            z_index_uv_md1(z_index_mu1(nm))    = ip
            uv_flags_iref(ip)                  = z_flags_iref(nm) + 1
            uv_flags_dir(ip)                   = 0 ! u point
            uv_flags_type(ip)                  = 1! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            uv_flags_vis(ip)                   = 1 ! enable viscosity
            if (advection .or. coriolis) uv_flags_adv(ip) = 1 ! enable advection
            !
         endif   
         !
         if (z_index_mu2(nm)>0) then
            !
            ip                                 = ip + 1
            uv_index_z_nm(ip)                  = nm
            uv_index_z_nmu(ip)                 = z_index_mu2(nm)
            z_index_uv_mu2(nm)                 = ip 
            z_index_uv_md1(z_index_mu2(nm))    = ip
            uv_flags_iref(ip)                  = z_flags_iref(nm) + 1
            uv_flags_dir(ip)                   = 0 ! u point
            uv_flags_type(ip)                  = 1 ! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            uv_flags_vis(ip)                   = 1 ! enable viscosity
            if (advection .or. coriolis) uv_flags_adv(ip) = 1 ! enable advection
            !
         endif   
         !            
      endif   
      !
      ! Above
      !
      if (z_index_nu(nm)==0) then
         !
         ! Same level
         !
         if (z_index_nu1(nm)>0) then
            !
            ip                                 = ip + 1
            uv_index_z_nm(ip)                  = nm
            uv_index_z_nmu(ip)                 = z_index_nu1(nm)
            z_index_uv_nu1(nm)                 = ip 
            z_index_uv_nd1(z_index_nu1(nm))    = ip
            uv_flags_iref(ip)                  = z_flags_iref(nm) 
            uv_flags_dir(ip)                   = 1 ! v point
            uv_flags_type(ip)                  = 0 ! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            uv_flags_vis(ip)                   = 1 ! enable viscosity
            if (advection .or. coriolis) uv_flags_adv(ip) = 1 ! enable advection
            !
         endif   
         !
      elseif (z_index_nu(nm)==-1) then
         !
         ! Coarser above
         !
         if (z_index_nu1(nm)>0) then
            !
            ip                                 = ip + 1
            uv_index_z_nm(ip)                  = nm
            uv_index_z_nmu(ip)                 = z_index_nu1(nm)
            z_index_uv_nu1(nm)                 = ip 
            if (z_index_nd1(z_index_nu1(nm)) == nm) then
                z_index_uv_nd1(z_index_nu1(nm)) = ip
            else
                z_index_uv_nd2(z_index_nu1(nm)) = ip
            endif   
            uv_flags_iref(ip)                  = z_flags_iref(nm) 
            uv_flags_dir(ip)                   = 1 ! u point
            uv_flags_type(ip)                  = -1! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            uv_flags_vis(ip)                   = 1 ! enable viscosity
            if (advection .or. coriolis) uv_flags_adv(ip) = 1 ! enable advection
            !
         endif   
         !
      else
         !
         ! Finer above
         !
         if (z_index_nu1(nm)>0) then
            !
            ip                                 = ip + 1
            uv_index_z_nm(ip)                  = nm
            uv_index_z_nmu(ip)                 = z_index_nu1(nm)
            z_index_uv_nu1(nm)                 = ip 
            z_index_uv_nd1(z_index_nu1(nm))    = ip
            uv_flags_iref(ip)                  = z_flags_iref(nm) + 1
            uv_flags_dir(ip)                   = 1 ! v point
            uv_flags_type(ip)                  = 1! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            uv_flags_vis(ip)                   = 1 ! enable viscosity
            if (advection .or. coriolis) uv_flags_adv(ip) = 1 ! enable advection
            
         endif   
         !
         if (z_index_nu2(nm)>0) then
            !
            ip                                 = ip + 1
            uv_index_z_nm(ip)                  = nm
            uv_index_z_nmu(ip)                 = z_index_nu2(nm)
            z_index_uv_nu2(nm)                 = ip 
            z_index_uv_nd1(z_index_nu2(nm))    = ip
            uv_flags_iref(ip)                  = z_flags_iref(nm) + 1
            uv_flags_dir(ip)                   = 1 ! v point
            uv_flags_type(ip)                  = 1 ! -1 is fine too coarse, 0 is normal, 1 is coarse to fine
            uv_flags_vis(ip)                   = 1 ! enable viscosity
            if (advection .or. coriolis) uv_flags_adv(ip) = 1 ! enable advection
            !
         endif   
         !            
      endif   
      !
   enddo 
   !
   ! And now for the uv neighbors of the uv points
   !
   do ip = 1, npuv
      !
      if (uv_flags_dir(ip)==0) then
         !
         ! U point
         !
         nm = uv_index_z_nm(ip)
         !
         ! left
         !
         if (z_index_md(nm)<1) then
            !
            ! Coarser or same level to the left
            !
            uv_index_u_nmd(ip) = z_index_uv_md1(nm)
            !
         endif
         !
         ! right
         !
         if (z_index_mu(nm)==0) then
            !
            ! Same level to the right
            !
            nmu = z_index_mu1(nm)
            !
            if (nmu>0) then
               if (z_index_mu(nmu)<1) then
                  !
                  ! Not finer to the right of nmu
                  !
                  uv_index_u_nmu(ip) = z_index_uv_mu1(nmu)
                  !
               endif
            endif
            !
         endif
         !
         ! below
         !
         if (z_index_nd(nm)==0) then
            !
            ! Same level below
            !
            ndm = z_index_nd1(nm)
            !
            if (ndm>0) then
               if (z_index_mu(ndm)<1) then
                  !
                  ! Not finer to the right of ndm
                  !
                  uv_index_u_ndm(ip) = z_index_uv_mu1(ndm)
                  !
               endif
            endif
            !
         elseif (z_index_nd(nm)==-1) then
            !
            ! Coarser below
            !
            ndm = z_index_nd1(nm)
            !
            if (z_index_nu2(ndm)==nm) then
               if (z_index_mu(ndm)==1) then
                  !
                  ! Finer to the right of ndm
                  !
                  uv_index_u_ndm(ip) = z_index_uv_mu2(ndm)
                  !
               endif
            endif
            !
         endif
         !
         ! 6
         !
         if (z_index_nu(nm)==0) then
            !
            ! Same level above
            !
            num = z_index_nu1(nm)
            !
            if (num>0) then
               if (z_index_mu(num)<1) then
                  !
                  ! Not finer to the right of num
                  !
                  uv_index_u_num(ip) = z_index_uv_mu1(num)
                  !
               endif
            endif
            !
         elseif (z_index_nu(nm)==-1) then
            !
            ! Coarser above
            !
            num = z_index_nu1(nm)
            !
            if (z_index_nd2(num)==nm) then
               if (z_index_mu(num)==1) then
                  !
                  ! Finer to the right of num
                  !
                  uv_index_u_num(ip) = z_index_uv_mu1(num)
                  !
               endif
            endif   
            !
         endif
         !
         if (advection .or. coriolis) then
            !
            ! 7
            !
            if (z_index_nd(nm)<1) then
               !
               ! Same or coarser below
               !
               uv_index_v_ndm(ip) = z_index_uv_nd1(nm)
               !
            endif   
            !
            ! 8
            !
            if (z_index_mu(nm)==0) then
               !
               ! Same on the right
               !
               nmu = z_index_mu1(nm)
               !
               if (nmu>0) then
                  !
                  if (z_index_nd(nmu)<1) then
                     !
                     ! Same or coarser below
                     !
                     uv_index_v_ndmu(ip) = z_index_uv_nd1(nmu)
                     !
                  endif   
                  !
               endif
               !
            endif
            !               
            ! 9 
            !
            if (z_index_nu(nm)<1) then
               !
               ! Same or coarser above
               !
               uv_index_v_nm(ip) = z_index_uv_nu1(nm)
               !
            endif   
            !
            ! 10
            !
            if (z_index_mu(nm)==0) then
               !
               ! Same to the right
               !
               nmu = z_index_mu1(nm)
               !
               if (nmu>0) then
                  !
                  if (z_index_nu(nmu)<1) then
                     !
                     ! Same or coarser above
                     !
                     uv_index_v_nmu(ip) = z_index_uv_nu1(nmu)
                     !
                  endif   
                  !
               endif
               !
            endif
            !                  
         endif   
         !
      else
         !
         ! V point
         !
         nm = uv_index_z_nm(ip)
         !
         ! left
         !
         if (z_index_nd(nm)<1) then
            !
            ! Coarser or same level to the left
            !
            uv_index_u_nmd(ip) = z_index_uv_nd1(nm)
            !
         endif
         !
         ! right
         !
         if (z_index_nu(nm)==0) then
            !
            ! Same level to the right
            !
            num = z_index_nu1(nm)
            !
            if (num>0) then
               if (z_index_nu(num)<1) then
                  !
                  ! Not finer to the right of nmu
                  !
                  uv_index_u_nmu(ip) = z_index_uv_nu1(num)
                  !
               endif
            endif
            !
         endif
         !
         ! below
         !
         if (z_index_md(nm)==0) then
            !
            ! Same level below
            !
            nmd = z_index_md1(nm)
            !
            if (nmd>0) then
               if (z_index_nu(nmd)<1) then
                  !
                  ! Not finer to the right of ndm
                  !
                  uv_index_u_ndm(ip) = z_index_uv_nu1(nmd)
                  !
               endif
            endif
            !
         elseif (z_index_md(nm)==-1) then
            !
            ! Coarser below
            !
            nmd = z_index_md1(nm)
            !
            if (z_index_mu2(nmd)==nm) then
               if (z_index_nu(nmd)==1) then
                  !
                  ! Finer to the right of ndm
                  !
                  uv_index_u_ndm(ip) = z_index_uv_nu2(nmd)
                  !
               endif
            endif
            !
         endif
         !
         ! 6
         !
         if (z_index_mu(nm)==0) then
            !
            ! Same level above
            !
            nmu = z_index_mu1(nm)
            !
            if (nmu>0) then
               if (z_index_nu(nmu)<1) then
                  !
                  ! Not finer to the right of num
                  !
                  uv_index_u_num(ip) = z_index_uv_nu1(nmu)
                  !
               endif
            endif
            !
         elseif (z_index_mu(nm)==-1) then
            !
            ! Coarser above
            !
            nmu = z_index_mu1(nm)
            !
            if (z_index_md2(nmu)==nm) then
               if (z_index_nu(nmu)==1) then
                  !
                  ! Finer to the right of num
                  !
                  uv_index_u_num(ip) = z_index_uv_nu1(nmu)
                  !
               endif
            endif   
            !
         endif
         !
         if (advection .or. coriolis) then
            !
            ! 7
            !
            if (z_index_md(nm)<1) then
               !
               ! Same or coarser below
               !
               uv_index_v_ndm(ip) = z_index_uv_md1(nm)
               !
            endif   
            !
            ! 8
            !
            if (z_index_nu(nm)==0) then
               !
               ! Same on the right
               !
               num = z_index_nu1(nm)
               !
               if (num>0) then
                  !
                  if (z_index_nd(num)<1) then
                     !
                     ! Same or coarser below
                     !
                     uv_index_v_ndmu(ip) = z_index_uv_md1(num)
                     !
                  endif   
                  !
               endif
               !
            endif
            !               
            ! 9 
            !
            if (z_index_mu(nm)<1) then
               !
               ! Same or coarser above
               !
               uv_index_v_nm(ip) = z_index_uv_mu1(nm)
               !
            endif   
            !
            ! 10
            !
            if (z_index_nu(nm)==0) then
               !
               ! Same to the right
               !
               num = z_index_nu1(nm)
               !
               if (num>0) then
                  !
                  if (z_index_mu(num)<1) then
                     !
                     ! Same or coarser above
                     !
                     uv_index_v_nmu(ip) = z_index_uv_mu1(num)
                     !
                  endif   
                  !
               endif
               !
            endif
            !                  
         endif   
         !
      endif   
      !
   enddo
   !
   ! Set flags for cells that do not have 4 normal neighbors (used in continuity)
   !
   do nm = 1, np
      !
      z_flags_type(nm) = 0
      !
      if (use_quadtree) then   
         !
         if (z_index_uv_md1(nm)==0 .or. z_index_uv_md2(nm)>0) then
            z_flags_type(nm) = 1
         endif
         !
         if (z_index_uv_mu1(nm)==0 .or. z_index_uv_mu2(nm)>0) then
            z_flags_type(nm) = 1
         endif
         !
         if (z_index_uv_nd1(nm)==0 .or. z_index_uv_nd2(nm)>0) then
            z_flags_type(nm) = 1
         endif
         !
         if (z_index_uv_nu1(nm)==0 .or. z_index_uv_nu2(nm)>0) then
            z_flags_type(nm) = 1
         endif
         !
      else
         !
         if (z_index_uv_md1(nm)==0) then
            z_flags_type(nm) = 1
         endif
         !
         if (z_index_uv_mu1(nm)==0) then
            z_flags_type(nm) = 1
         endif
         !
         if (z_index_uv_nd1(nm)==0) then
            z_flags_type(nm) = 1
         endif
         !
         if (z_index_uv_nu1(nm)==0) then
            z_flags_type(nm) = 1
         endif
         !
      endif   
      !
   enddo      
!   !
!   ! Set flags to turn off viscosity in some cells
!   !
!   do ip = 1, npuv
!      do j = 3, 6
!         if (uv_index(j, ip)==0) then
!            uv_flags(4, ip) = 0
!         endif   
!      enddo   
!   enddo      
!   !
!   ! Set flags to turn off advection in some cells
!   !
!   if (advection .or. coriolis) then
!      do ip = 1, npuv
!         do j = 3, 10
!            if (uv_index(j, ip)==0) then
!               uv_flags(5, ip) = 0
!            endif   
!         enddo   
!      enddo      
!   endif
   !
   ! Okay, got all the quadtree cells including indices and flags 
   !
   write(*,*)'Number of active z points    : ', np
   write(*,*)'Number of active u/v points  : ', npuv
   !
   if (use_quadtree) then
      do iref = 1, nref
         write(*,'(a,i3,a,i8)')' Number of cells in level ', iref, ' : ', quadtree_last_point_per_level(iref) - quadtree_first_point_per_level(iref) + 1
      enddo
   endif
   !
   ! MESH GEOMETRY
   !
   ! Refinement levels
   !
   allocate(dxr(nref))
   allocate(dyr(nref))
   !
   ! Determine grid spacing for different refinement levels
   !
   do iref = 1, nref
      !
      dxr(iref) = dx/(2**(iref - 1))
      dyr(iref) = dy/(2**(iref - 1))
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
      z_xz(nm) = x0 + cosrot*(1.0*(m - 0.5))*dxr(iref) - sinrot*(1.0*(n - 0.5))*dyr(iref)
      z_yz(nm) = y0 + sinrot*(1.0*(m - 0.5))*dxr(iref) + cosrot*(1.0*(n - 0.5))*dyr(iref)
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
      allocate(dxm(npuv)) 
      allocate(dxminv(npuv))   
      allocate(dxm2inv(npuv))
      !
      allocate(dyrm(nref))
      allocate(dyrinv(nref))
      allocate(dyr2inv(nref))
      allocate(dyrinvc(nref))
      !      
      allocate(cell_area_m2(np))
      !
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
         dxm(ip)     = dxr(iref)*111111.1*cos(z_yz(nm)*pi/180)
         dxminv(ip)  = 1.0/dxm(ip)
         dxm2inv(ip) = dxminv(ip)**2
         !
         ! Minimum grid spacing to determine dtmax 
         !  
         dxymin = min(dxymin, dxm(iref))
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
            dyrinvc(iref) = 1.0*(3*dyrm(iref)/2)
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
         cell_area_m2(nm) = dxr(iref)*111111.1*cos(z_yz(nm)*pi/180)*dyrm(iref) 
         !         
      enddo   
      !
      ! Spatially varying Coriolis
      !
      if (use_coriolis) then
         do nm = 1, np            
            fcorio2d(nm) = 2*7.2921e-05*sin(z_yz(nm)*pi/180)
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
         dxrinv(iref)    = 1.0/dxrm(iref)
         dxr2inv(iref)   = dxrinv(iref)**2
         !
         dyrm(iref)      = dyr(iref)
         dyrinv(iref)    = 1.0/dyrm(iref)
         dyr2inv(iref)   = dyrinv(iref)**2
         !
         cell_area(iref) = dxrm(iref)*dyrm(iref)
         !
         if (iref>1) then
            !         
            dxrinvc(iref) = 1.0/(3*dxrm(iref)/2)
            dyrinvc(iref) = 1.0*(3*dyrm(iref)/2)
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
   endif   
   !
   ! Determine minimum and maximum time step
   !
   dtmax = min(dtmax, dxymin/(sqrt(9.81*0.1)))
   dtmin = alfa*dxymin/(sqrt(9.81*stopdepth)) ! If dt falls below this value, the simulation will stop
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
            write(*,*)'Reading ',trim(depfile)
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
         write(*,*)'Reading ',trim(depfile)
         !
         if (inputtype=='asc') then
            !
            allocate(zbg(nmax, mmax))
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
            deallocate(kcsg)
            !
         else
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
         read(500)subgrid_nbins
         subgrid_nbins = subgrid_nbins + 1
         allocate(subgrid_z_zmin(np))
         allocate(subgrid_z_zmax(np))
         allocate(subgrid_z_zmean(np))
         allocate(subgrid_z_volmax(np))
         allocate(subgrid_z_dep(subgrid_nbins, np))
         allocate(subgrid_uv_zmin(npuv))
         allocate(subgrid_uv_zmax(npuv))
         allocate(subgrid_uv_hrep(subgrid_nbins, npuv))
         allocate(subgrid_uv_navg(subgrid_nbins, npuv))
         allocate(subgrid_uv_hrep_zmax(npuv))
         allocate(subgrid_uv_navg_zmax(npuv))
         !
         allocate(rtmpz(npzq))
         allocate(rtmpuv(npuvq))
         allocate(uv_index_qt_in_sf(npuv))
         !
         write(*,*)'Number of subgrid bins : ',subgrid_nbins - 1
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
         do nm = 1, np
            subgrid_z_zmean(nm) = rtmpz(index_quadtree_in_sfincs(nm))
         enddo   
         !
         read(500)rtmpz
         do nm = 1, np
            subgrid_z_volmax(nm) = rtmpz(index_quadtree_in_sfincs(nm))
         enddo   
         !
         subgrid_z_dep(1, :) = max(subgrid_z_zmin(:), -20.0)
         do ibin = 1, subgrid_nbins - 1
            read(500)rtmpz
            do nm = 1, np
               subgrid_z_dep(ibin + 1, nm) = rtmpz(index_quadtree_in_sfincs(nm))
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
         ! Initialize subgrid_uv_hrep with huthresh
         !
         subgrid_uv_hrep = huthresh
         do ibin = 1, subgrid_nbins - 1
            read(500)rtmpuv
            do nm = 1, npuv
               subgrid_uv_hrep(ibin + 1, nm) = max(rtmpuv(uv_index_qt_in_sf(nm)), huthresh)
            enddo   
         enddo
         !
         do ibin = 1, subgrid_nbins - 1
            read(500)rtmpuv
            do nm = 1, npuv
               ! Already convert here to gn^2
               subgrid_uv_navg(ibin + 1, nm) = g*max(rtmpuv(uv_index_qt_in_sf(nm)), 0.005)**2
            enddo   
         enddo
         ! Set bottom navg equal to one bin above
         subgrid_uv_navg(1, :) = subgrid_uv_navg(2, :)
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
         read(500)idummy ! option
         read(500)subgrid_nbins
         subgrid_nbins = subgrid_nbins + 1
         allocate(subgrid_z_zmin(np))
         allocate(subgrid_z_zmax(np))
         allocate(subgrid_z_volmax(np))
         allocate(subgrid_z_dep(subgrid_nbins, np))
         allocate(subgrid_uv_zmin(npuv))
         allocate(subgrid_uv_zmax(npuv))
         allocate(subgrid_uv_hrep(subgrid_nbins, npuv))
         allocate(subgrid_uv_navg(subgrid_nbins, npuv))
         allocate(subgrid_uv_hrep_zmax(npuv))
         allocate(subgrid_uv_navg_zmax(npuv))
         !
         allocate(rtmpz(np))
         !
         read(500)subgrid_z_zmin
         read(500)subgrid_z_zmax
         read(500)subgrid_z_volmax
         !
         ! Make sure zmax is always bigger than zmin
         !
         do nm = 1, np
            if (subgrid_z_zmax(nm) - subgrid_z_zmin(nm) < 0.01) subgrid_z_zmax(nm) = subgrid_z_zmax(nm) + 0.01
         enddo
         !
         subgrid_z_dep(1,:) = max(subgrid_z_zmin(:), -20.0)
         do ibin = 1, subgrid_nbins - 1
            read(500)rtmpz
            do nm = 1, np
               subgrid_z_dep(ibin + 1, nm) = rtmpz(nm)
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
         subgrid_uv_hrep(1,:) = huthresh
         do ibin = 1, subgrid_nbins - 1
            read(500)rtmpz
            do nm = 1, np
               ip = z_index_uv_mu1(nm) ! index of uv point
               if (ip>0) then
                  subgrid_uv_hrep(ibin + 1, ip) = rtmpz(nm)
               endif
            enddo
         enddo         
         !
         ! U navg
         !
         do ibin = 1, subgrid_nbins - 1
            !
            read(500)rtmpz
            !            
            do nm = 1, np
               !
               ip = z_index_uv_mu1(nm) ! index of uv point
               !
               if (ip>0) then
                  !
                  subgrid_uv_navg(ibin + 1, ip) = g*max(rtmpz(nm), 0.005)**2
                  !
               endif
               !
            enddo
            !
         enddo
         !
         do nm = 1, np
            subgrid_uv_navg(1, nm) = subgrid_uv_navg(2, nm)
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
         !
         subgrid_uv_hrep(1,:) = huthresh
         !
         do ibin = 1, subgrid_nbins - 1
            !
            read(500)rtmpz
            !            
            do nm = 1, np
               !
               ip = z_index_uv_nu1(nm) ! index of uv point
               !
               if (ip>0) then
                  !
                  subgrid_uv_hrep(ibin + 1, ip ) = max(rtmpz(nm), huthresh)
                  !
               endif
               !
            enddo
            !
         enddo
         !
         ! V navg
         !
         !
         do ibin = 1, subgrid_nbins - 1
            !
            read(500)rtmpz
            !            
            do nm = 1, np
               !
               ip = z_index_uv_nu1(nm) ! index of uv point
               !
               if (ip>0) then
                 !
                  subgrid_uv_navg(ibin + 1, ip) = g*max(rtmpz(nm), 0.005)**2
                  !
               endif
               !
            enddo
            !
         enddo
         !
         do nm = 1, np
            subgrid_uv_navg(1, nm) = subgrid_uv_navg(2, nm)
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
         if (subgrid_z_zmax(nm) - subgrid_z_zmin(nm) < 0.01) subgrid_z_zmax(nm) = subgrid_z_zmax(nm) + 0.01
      enddo
      !
      do nm = 1, npuv
         if (subgrid_uv_zmax(nm) - subgrid_uv_zmin(nm) < 0.01) subgrid_uv_zmax(nm) = subgrid_uv_zmax(nm) + 0.01
      enddo
      !
      ! Make arrays for subgrid_uv_hrep_zmax and subgrid_uv_navg_zmax for faster searching
      !
      do ip = 1, npuv
         subgrid_uv_hrep_zmax(ip) = subgrid_uv_hrep(subgrid_nbins, ip) - subgrid_uv_zmax(ip)
         subgrid_uv_navg_zmax(ip) = subgrid_uv_navg(subgrid_nbins, ip)
      enddo    
      !
   endif
   !
   ! BOUNDARIES
   !
   ! First count number of boundary cells
   !
   ngbnd = 0
   do nm = 1, np
      if (kcs(nm)==2 .or. kcs(nm)==3) then
         !
         include_boundaries = .true.
         ngbnd = ngbnd + 1
         !
      endif
   enddo
   !
   ! Now allocate them
   !
   allocate(nmindbnd(ngbnd))
   allocate(zsb(ngbnd))
   allocate(zsb0(ngbnd))
   allocate(ibndtype(ngbnd)) ! 0 for outflow boundary, 1 for regular boundary
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
      if (kcs(nm)==2 .or. kcs(nm)==3) then
         !
         ! This is a boundary point
         !
         ngbnd = ngbnd + 1
         nmindbnd(ngbnd) = nm
         !
         ! Determine whether this is a regular or outflow grid point
         !
         if (kcs(nm)==2) then
            !
            ! Regular boundary point
            !
            ibndtype(ngbnd) = 1
            !
         else
            !
            ! Outflow boundary point (set kcs back to 2)
            !
            ibndtype(ngbnd) = 0
!            kcs(nm) = 2
            !
         endif   
      endif
   enddo
   !
   ! UV boundary points
   !
   nkcuv2 = 0 ! Number of points with kcuv=2
   kcuv   = 0
   !
   do ip = 1, npuv
      !
      nm  = uv_index_z_nm(ip)
      nmu = uv_index_z_nmu(ip)
      !
      if (kcs(nm)>1 .and. kcs(nmu)>1)  then
         !
         ! Two surrounding boundary points
         !
         kcuv(ip) = 0 ! Is this really necessary
         !
      elseif (kcs(nm)*kcs(nmu) == 0) then
            !
            ! One or two surrounding inactive points, so make this velocity point inactive
            !
            kcuv(ip) = 0
            !
      elseif ((kcs(nm)==1 .and. kcs(nmu)>1) .or. (kcs(nm)>1 .and. kcs(nmu)==1))  then
         !
         ! One surrounding boundary points, so make this velocity point a boundary point
         !
         kcuv(ip) = 2
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
      if (kcuv(ip)==2) then
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
            if (nmindbnd(ib)==nmbkcuv2(ikcuv2)) then
               !
               ibkcuv2(ikcuv2) = ib ! index in grid boundary array of this uv boundary point
               !
               ! Set water levels for this point (this only happens here)
               !
               if (ibndtype(ib)==0) then
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
   ! Set flags to turn off viscosity in some cells
   !
   if (nmax/=1) then
      do ip = 1, npuv
         ! 
         if (uv_index_u_nmd(ip) == 0)  uv_flags_vis(ip) = 0  
         if (uv_index_u_nmu(ip) == 0)  uv_flags_vis(ip) = 0
         if (uv_index_u_num(ip) == 0)  uv_flags_vis(ip) = 0
         if (uv_index_u_ndm(ip) == 0)  uv_flags_vis(ip) = 0
         !      
      enddo      
   endif 
   !
   ! Turn off advection in some cells (next to inactive points and boundaries)
   !
   if ((advection .or. coriolis) .and. nmax/=1) then
      !
      do ip = 1, npuv
         !
         ! Turn off advection near inactive points
         !
         if (uv_index_u_nmd(ip) == 0)  uv_flags_adv(ip) = 0
         if (uv_index_u_nmu(ip) == 0)  uv_flags_adv(ip) = 0
         if (uv_index_u_num(ip) == 0)  uv_flags_adv(ip) = 0
         if (uv_index_u_ndm(ip) == 0)  uv_flags_adv(ip) = 0
         if (uv_index_v_ndm(ip) == 0)  uv_flags_adv(ip) = 0
         if (uv_index_v_nm(ip) == 0)   uv_flags_adv(ip) = 0
         if (uv_index_v_nmu(ip) == 0)  uv_flags_adv(ip) = 0
         if (uv_index_v_ndmu(ip) == 0) uv_flags_adv(ip) = 0
         !
         ! Turn off advection near the boundaries
         !
         if (uv_index_u_nmd(ip)>0) then
            if (kcuv(uv_index_u_nmd(ip)) == 2)  uv_flags_adv(ip) = 0
         endif   
         if (uv_index_u_nmu(ip)>0) then
            if (kcuv(uv_index_u_nmu(ip)) == 2)  uv_flags_adv(ip) = 0
         endif   
         if (uv_index_u_num(ip)>0) then
            if (kcuv(uv_index_u_num(ip)) == 2)  uv_flags_adv(ip) = 0
         endif   
         if (uv_index_u_ndm(ip)>0) then
            if (kcuv(uv_index_u_ndm(ip)) == 2)  uv_flags_adv(ip) = 0
         endif   
         if (uv_index_v_ndm(ip)>0) then
            if (kcuv(uv_index_v_ndm(ip)) == 2)  uv_flags_adv(ip) = 0
         endif   
         if (uv_index_v_nm(ip)>0) then
            if (kcuv(uv_index_v_nm(ip)) == 2)   uv_flags_adv(ip) = 0
         endif   
         if (uv_index_v_nmu(ip)>0) then
            if (kcuv(uv_index_v_nmu(ip)) == 2)  uv_flags_adv(ip) = 0
         endif   
         if (uv_index_v_ndmu(ip)>0) then
            if (kcuv(uv_index_v_ndmu(ip)) == 2) uv_flags_adv(ip) = 0
         endif   
         !
      enddo
      !
   endif
   !
   ! FRICTION COEFFICIENTS
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
         !
         allocate(rghfield(np))
         write(*,*)'Reading ',trim(manningfile)
         open(unit = 500, file = trim(manningfile), form = 'unformatted', access = 'stream')
         read(500)rghfield
         close(500)
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
   endif
   !
   ! INFILTRATION
   !
   ! Infiltration only works when rainfall is activated ! If you want infiltration without rainfall, use a precip file with 0.0s
   !
   ! Note, infiltration methods not designed to be stacked
   !
   infiltration   = .false.
   !
   ! Four options for infiltration:
   !
   ! 1) Spatially-uniform constant infiltration
   !    Requires: -
   ! 2) Spatially-varying constant infiltration
   !    Requires: qinfmap (does not require qinffield !)
   ! 3) Spatially-varying infiltration with CN numbers (old)
   !    Requires: cumprcp, cuminf, qinfmap, qinffield
   ! 4) Spatially-varying infiltration with CN numbers (new)
   !    Requires: qinfmap, qinffield, qinffield, scs_kr, scs_P1, scs_F1, scs_Se and scs_rain (but not necessarily cuminf and cumprcp)
   !
   ! cumprcp and cuminf are stored in the netcdf output if store_cumulative_precipitation == .true. which is the default
   !
   ! We need to keep cumprcp and cuminf in memory when:
   !   a) store_cumulative_precipitation == .true.
   ! or:  
   !   b) inftype == 'cna' or inftype == 'cnb'
   !
   ! First we determine precipitation type
   !
   if (precip) then
      !
      if (qinf > 0.0) then   
         !
         ! Spatially-uniform constant infiltration (specified as +mm/hr)
         !
         inftype = 'con'
         infiltration = .true.
         ! 
      elseif (qinffile /= 'none') then
         !
         ! Spatially-varying constant infiltration
         !
         inftype = 'c2d'
         infiltration = .true.      
         !
      elseif (scsfile /= 'none') then
         !
         ! Spatially-varying infiltration with CN numbers (old)
         !
         inftype = 'cna'
         infiltration = .true.      
         !
      elseif (scsfile_Se /= 'none') then  
         !
         ! Spatially-varying infiltration with CN numbers (new)
         !
         inftype = 'cnb'
         infiltration = .true.      
         !
      endif
      !
      if (precip) then
         !
         ! We need cumprcp and cuminf
         !
         allocate(cumprcp(np))
         cumprcp = 0.0
         !
         allocate(cuminf(np))
         cuminf = 0.0
         ! 
      endif
      !
      ! Now allocate and read spatially-varying inputs 
      !
      if (infiltration) then
         !
         allocate(qinfmap(np))
         qinfmap = 0.0
         ! 
      endif
      !
      if (inftype == 'con') then   
         !
         ! Spatially-uniform constant infiltration (specified as +mm/hr)
         !
         write(*,*)'Turning on process: Spatially-uniform constant infiltration'        
         !
         allocate(qinffield(np))
         do nm = 1, np
             if (subgrid) then
                 if (subgrid_z_zmin(nm) > qinf_zmin) then
                    qinffield(nm) = qinf
                 else
                    qinffield(nm) = 0.0
                 endif
             else
                 if (zb(nm) > qinf_zmin) then
                    qinffield(nm) = qinf
                 else
                    qinffield(nm) = 0.0
                 endif
             endif
         enddo
         !
      elseif (inftype == 'c2d') then
         !
         ! Spatially-varying constant infiltration
         !
         write(*,*)'Turning on process: Spatially-varying constant infiltration'      
         !
         ! Read spatially-varying infiltration (only binary, specified in +mm/hr)
         !
         write(*,*)'Reading ', trim(qinffile), ' ...'
         allocate(qinffield(np))
         open(unit = 500, file = trim(qinffile), form = 'unformatted', access = 'stream')
         read(500)qinffield
         close(500)
         !
         do nm = 1, np
            qinffield(nm) = qinffield(nm)/3.6e3/1.0e3   ! convert to +m/s
         enddo
         !
      elseif (inftype == 'cna') then
         !
         ! Spatially-varying infiltration with CN numbers (old)
         !
         write(*,*)'Turning on process: Infiltration (via CN method - A)'               
         !
         allocate(qinffield(np))
         qinffield = 0.0
         !
         write(*,*)'Reading ',trim(scsfile)
         open(unit = 500, file = trim(scsfile), form = 'unformatted', access = 'stream')
         read(500)qinffield
         close(500)
         !
         ! already convert qinffield from inches to m here
         do nm = 1, np
            qinffield(nm) = qinffield(nm)*0.0254   !to m
         enddo      
         !
      elseif (inftype == 'cnb') then  
         !
         ! Spatially-varying infiltration with CN numbers (new)
         !
         write(*,*)'Turning on process: Infiltration (via CN method - B)'            
         ! 
         ! Allocate Smax
         ! 
         allocate(qinffield(np))
         qinffield = 0.0
         write(*,*)'Reading ',trim(scsfile_Smax)
         open(unit = 500, file = trim(scsfile_Smax), form = 'unformatted', access = 'stream')
         read(500)qinffield
         close(500)
         !
         ! Allocate Se
         !
         allocate(qinffield2(np))
         qinffield2 = 0.0
         write(*,*)'Reading ',trim(scsfile_Se)
         open(unit = 501, file = trim(scsfile_Se), form = 'unformatted', access = 'stream')
         read(501)qinffield2
         close(501)
         !
         ! Allocate kr
         !
         allocate(scs_kr(np))
         scs_kr = 0.0
         write(*,*)'Reading ',trim(scsfile_kr)
         open(unit = 502, file = trim(scsfile_kr), form = 'unformatted', access = 'stream')
         read(502)scs_kr
         close(502)
         !
         ! Allocate support variables
         !
         allocate(scs_P1(np))
         scs_P1 = 0.0
         allocate(scs_F1(np))
         scs_F1 = 0.0
         allocate(scs_Se(np))
         scs_Se = 0.0
         allocate(scs_rain(np))
         scs_rain = 0
         !
      endif
      !
   else
      !
      ! Overrule input
      !
      store_cumulative_precipitation = .false.
      !
   endif
   !
   end subroutine


   subroutine initialize_hydro()
   !
   ! Initialize SFINCS variables (qx, qy, zs etc.)
   !
   use sfincs_data
   !
   implicit none
   !
   real*4, dimension(:,:), allocatable :: inilev
   real*4, dimension(:),   allocatable :: inizs
   real*4, dimension(:),   allocatable :: iniq
   !
   integer    :: nm, m, n, ivol, num, nmu, ibin, rsttype, ip, ind, iuv
   real*4     :: dzvol
   real*4     :: facint
   real*4     :: rdummy
   real*4     :: dzuv
   real*4     :: zmax
   real*4     :: zmin
   real*4     :: one_minus_facint 
   real*4     :: huv
   real*4     :: zsuv   
   !
   logical   :: iok
   !
   allocate(zs(np))
   allocate(q(npuv))
   allocate(q0(npuv))
   allocate(uv(npuv))
   allocate(uv0(npuv))
   allocate(kfuv(npuv)) ! Not needed anymore?
   !
   zs   = 0.0
   q    = 0.0
   q0   = 0.0
   uv   = 0.0
   uv0  = 0.0
   kfuv = 0
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
      if (store_wave_forces) then
         allocate(fwx(np))
         allocate(fwy(np))
         fwx = 0.0
         fwy = 0.0
      endif
      !
      if (store_wave_direction) then
         allocate(mean_wave_direction(np))
         allocate(wave_directional_spreading(np))
         mean_wave_direction        = 0.0
         wave_directional_spreading = 0.0
      endif   
      !
   endif   
   !
   if (wavemaker) then 
      allocate(zsm(np))
   endif
   !
   if (store_maximum_waterlevel) then
      allocate(zsmax(np))
   endif
   !
   if (store_maximum_velocity) then
      allocate(vmax(np))
   endif
   if (store_twet) then
       allocate(twet(np))
   endif
   !
   if (store_tsunami_arrival_time) then
      allocate(tsunami_arrival_time(np))
   endif
   !
   ! Make u and v points always active when they are boundary points
   !
   do ip = 1, npuv
      if (kcuv(ip)==2) kfuv(ip) = 1 ! Get rid of this
   enddo
   !
   ! Initialize water level points
   !
   if (rstfile(1:4) /= 'none') then
      !
      ! Binary restart file
      !
      allocate(inizs(np))
      allocate(iniq(npuv))
      !
      write(*,*)'Reading restart file ', trim(rstfile), ' ...'
      !
      open(unit = 500, file = trim(rstfile), form = 'unformatted', access = 'stream')
      !
      ! Type of restart - 1: zs, qx, qy, umean and vmean  - 2: zs, qx, qy - 3: zs
      !
      read(500)rdummy
      read(500)rsttype
      read(500)rdummy
      !
      ! Always read in inizs
      !
      read(500)rdummy
      read(500)inizs
      read(500)rdummy
      !      
      if (rsttype==1 .or. rsttype==2) then     
         read(500)rdummy
         read(500)iniq
         read(500)rdummy
      endif
      !
      if (rsttype==1) then     
         read(500)rdummy
         read(500)uvmean
         read(500)rdummy
      endif
      !   
      close(500)      
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
      if (rsttype==1 .or. rsttype==2) then
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
                     huv    = subgrid_uv_hrep_zmax(ip) + zsuv
                     !
                  else
                     !
                     ! Interpolation required
                     !
                     dzuv   = (subgrid_uv_zmax(ip) - subgrid_uv_zmin(ip)) / (subgrid_nbins - 1)
                     iuv    = int((zsuv - subgrid_uv_zmin(ip))/dzuv) + 1
                     facint = (zsuv - (subgrid_uv_zmin(ip) + (iuv - 1)*dzuv) ) / dzuv
                     huv    = subgrid_uv_hrep(iuv, ip) + (subgrid_uv_hrep(iuv + 1, ip) - subgrid_uv_hrep(iuv, ip))*facint
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
               uv(ip)   = max(min(q(ip)/huv, 5.0), -5.0)   
               !
            endif   
            !
         enddo
         !         
      endif   
      !
      deallocate(inizs)
      deallocate(iniq)
      !
   elseif (zsinifile(1:4) /= 'none') then ! Read binary (!) initial water level file
      !
      allocate(inizs(np))
      !       
      write(*,*)'Reading ',trim(zsinifile)
      open(unit = 500, file = trim(zsinifile), form = 'unformatted', access = 'stream')
      read(500)inizs
      close(500)       
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
      deallocate(inizs)    
      !
   else
      !
      ! No initial conditions file
      !
      if (subgrid) then
         do nm = 1, np
            zs(nm)   = max(subgrid_z_zmin(nm), zini) ! Water level at zini or bed level (whichever is higher)
         enddo
      else
         do nm = 1, np
            zs(nm)   = max(zb(nm), zini) ! Water level at zini or bed level (whichever is higher)
         enddo
      endif
      !
   endif
   !
   if (wavemaker) then
      !      
      zsm = zs
      !
   endif   
   !
   if (store_maximum_waterlevel) then
      zsmax = -999.0
   endif
   !
   if (store_maximum_velocity) then
      vmax = 0.0
   endif
   !
   if (store_twet) then
      twet = 0.0
   endif
   !
   uv = 0.0
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
      if (store_meteo) then
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
      allocate(cumprcpt(np))      
      cumprcpt = 0.0
      allocate(netprcp(np))      
      netprcp  = 0.0
      !
   endif
   !
   if (subgrid) then
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
            do ibin = 2, subgrid_nbins
               if (subgrid_z_dep(ibin, nm)>zs(nm)) then
                  ivol = ibin - 1
                  exit
               endif
            enddo
            !
            dzvol  = subgrid_z_volmax(nm) / (subgrid_nbins - 1)
            facint = (zs(nm) - subgrid_z_dep(ivol, nm)) / max(subgrid_z_dep(ivol + 1, nm) - subgrid_z_dep(ivol, nm), 0.001)
            z_volume(nm) = (ivol - 1)*dzvol + facint*dzvol
            !
         endif
         !
      enddo
   endif
   !
   end subroutine
   !
end module
