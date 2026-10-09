   module sfincs_structures

   use sfincs_error
   use sfincs_log
   
   contains

   subroutine read_structures()
   !
   ! Reads all structures file
   !
   use sfincs_data
   !
   implicit none
   !
   ! First thin dams (kcs will be set to 0)
   !
   call read_thin_dams()
   !
   ! And now structures (kcs will be set to 3)
   !
   nrstructures = 0
   !
   if (weirfile(1:4) /= 'none') then
      !
      ! Weir (type = 1), 2 parameters
      !
      ! Parameter 1 : Weir level (w.r.t. ref datum)
      ! Parameter 2 : Cd
      !
      call read_structure_file(weirfile, 1, 2)
      !
   endif
   !
   ! Flow-blocking features: no structures, the subgrid uv tables are raised to the crest level
   !
   nr_flow_blocking_points = 0
   !
   if (flow_blocking_file(1:4) /= 'none') then
      !
      call read_flow_blocking_file()
      !
   endif
   !
   end subroutine


   subroutine read_structure_file(filename,itype,npars)
   !
   ! Reads structures file
   !
   use sfincs_data
   use quadtree
   use geometry
   !
   implicit none
   !
   integer ip, nm, nmu, nthd, nrows, ncols, irow, stat, ithd, n, m, ii, itype
   integer nstruc, npars, ipar, iuv, nuv, irefnm, irefnmu
   !
   real      :: dummy
   real      :: xuv1
   real      :: yuv1
   real      :: xuv2
   real      :: yuv2
   real      :: xuv
   real      :: yuv
   real      :: xstr1
   real      :: ystr1
   real      :: xstr2
   real      :: ystr2
   real      :: dst1
   real      :: dst2
   real      :: wfac
   real      :: height
   real      :: x1
   real      :: y1
   real      :: x2
   real      :: y2
   real      :: x0a
   real      :: y0a
   real      :: x0b
   real      :: y0b
   real      :: d   
   integer   :: nr_points
   integer   :: indx
   integer   :: iref
   integer   :: idir
   real*4    :: xpa,xpb,ypa,ypb
   real*4    :: dxs
   logical   :: okay
   !
   character :: cdummy
   character*256 :: filename
   !
   real*4, dimension(:),   allocatable   :: xthd
   real*4, dimension(:),   allocatable   :: ythd
   real*4, dimension(:,:), allocatable   :: pars
   real*4, dimension(:,:), allocatable   :: strucpars
   real*4, dimension(:),   allocatable   :: lngth
   real*4, dimension(2)                  :: xp
   real*4, dimension(2)                  :: yp
   !
   integer,   dimension(:), allocatable  :: ipol
   integer,   dimension(:), allocatable  :: uv_indices   
   integer,   dimension(:), allocatable  :: vertices 
   integer*1, dimension(:), allocatable  :: istruc
   !
   write(logstr,'(a)')'Info    : reading weir file'
   call write_log(logstr, 0)
   !
   okay = check_file_exists(filename, 'Weir file', .true.)
   !
   ! Read structures file
   !
   nthd   = 0
   nstruc = nrstructures
   !
   allocate(istruc(npuv))
   istruc = 0
   !
   allocate(strucpars(npars,npuv))
   strucpars = 0.0
   strucpars(1,:) = -9999.0 ! height
   !
   allocate(lngth(npuv))
   lngth = 0.0
   !
   ! First count number of polylines
   !
   open(500, file=trim(filename))
   do while(.true.)
      read(500,*,iostat = stat)cdummy
      if (stat<0) exit     
      read(500,*,iostat = stat)nrows,ncols
      if (stat<0) exit
      nthd = nthd + 1
      do irow = 1, nrows
         read(500,*)dummy
      enddo
   enddo
   rewind(500)
   !
   ! Loop through polylines
   !
   do ithd = 1, nthd
      !
      read(500,*,iostat = stat)cdummy
      read(500,*,iostat = stat)nrows,ncols
      if (stat<0) exit
      allocate(xthd(nrows))
      allocate(ythd(nrows))
      allocate(pars(npars,nrows))
      select case(npars)
      case(1)
         do irow = 1, nrows
            read(500,*)xthd(irow),ythd(irow),pars(1,irow)
         enddo
      case(2)
         do irow = 1, nrows
            read(500,*)xthd(irow),ythd(irow),pars(1,irow),pars(2,irow)
         enddo
      case(3)
         do irow = 1, nrows
            read(500,*)xthd(irow),ythd(irow),pars(1,irow),pars(2,irow),pars(3,irow)
         enddo
      end select
      !
      call find_uv_points_intersected_by_polyline(uv_indices, vertices, nr_points, xthd, ythd, nrows)
      !
      ! Now loop through uv points to get the length and determine parameters
      !
      do iuv = 1, nr_points
         !
         indx = uv_indices(iuv)
         irow = vertices(iuv)
         istruc(indx) = 1
         !
         ! Projecting uv corner points onto structure segment to determine length of this uv point
         !
         nm   = uv_index_z_nm(indx)
         nmu  = uv_index_z_nmu(indx)
         !
         ! Need to find the two corner points that make up the velocity point
         !
         iref    = uv_flags_iref(indx)
         idir    = uv_flags_dir(indx)
         irefnm  = z_flags_iref(nm)
         irefnmu = z_flags_iref(nmu)
         !
         if (crsgeo) then
            ! Geographic coordinate system
            if (idir==0) then
               ! U point
               dxs = dyrm(iref)
            else
               ! V point
               dxs = 1.0 / dxminv(indx)
            endif   
         else   
            ! Projected coordinate system
            if (idir==0) then
               ! U point
               dxs  = dyrm(iref)
            else
               ! V point
               dxs  = dxrm(iref)
            endif
         endif
         !
         if (irefnm>=irefnmu) then
            !
            ! Normal point or fine to coarse
            !
            if (idir==0) then
               xuv  = z_xz(nm) + 0.5*dxs*cos(rotation) ! coordinates of uv point
               yuv  = z_yz(nm) + 0.5*dxs*sin(rotation)
            else   
               xuv  = z_xz(nm) + 0.5*dxs*cos(rotation + 0.5*pi) ! coordinates of uv point
               yuv  = z_yz(nm) + 0.5*dxs*sin(rotation + 0.5*pi)
            endif
            !
         else
            !
            if (idir==0) then
               xuv  = z_xz(nmu) + 0.5*dxs*cos(rotation - pi) ! coordinates of uv point
               yuv  = z_yz(nmu) + 0.5*dxs*sin(rotation - pi)
            else   
               xuv  = z_xz(nmu) + 0.5*dxs*cos(rotation - 0.5*pi) ! coordinates of uv point
               yuv  = z_yz(nmu) + 0.5*dxs*sin(rotation - 0.5*pi)
            endif
            !
         endif   
         !
         if (idir==0) then
            !
            ! U point
            !
            xuv1 = xuv + 0.5*dxs*cos(rotation - 0.5*pi)
            yuv1 = yuv + 0.5*dxs*sin(rotation - 0.5*pi)
            xuv2 = xuv + 0.5*dxs*cos(rotation + 0.5*pi)
            yuv2 = yuv + 0.5*dxs*sin(rotation + 0.5*pi)
            !
         else
            !
            ! V point
            !
            xuv1 = xuv + 0.5*dxs*cos(rotation - pi)
            yuv1 = yuv + 0.5*dxs*sin(rotation - pi)
            xuv2 = xuv + 0.5*dxs*cos(rotation)
            yuv2 = yuv + 0.5*dxs*sin(rotation)
            !
         endif   
         !
         ! Length
         !
         d = distance_between_points_projected_on_line_segment(xuv1, yuv1, xuv2, yuv2, xthd(irow), ythd(irow), xthd(irow + 1), ythd(irow + 1), 999999.0)
         !
         lngth(indx) = lngth(indx) + d
         !
         dst1   = sqrt((xuv - xthd(irow))**2 + (yuv - ythd(irow))**2)
         dst2   = sqrt((xuv - xthd(irow + 1))**2 + (yuv - ythd(irow + 1))**2)
         wfac   = dst1 / (dst1 + dst2)
         !
         height = pars(1, irow)*(1.0 - wfac) + pars(1, irow + 1)*wfac
         !
         strucpars(1, indx) = max(strucpars(1, indx), height)
         ! 
         do ipar = 2, npars
            strucpars(ipar, indx) = pars(ipar, irow)*(1.0 - wfac) + pars(ipar, irow + 1)*wfac
         enddo
         !
      enddo
      !
      deallocate(xthd)
      deallocate(ythd)
      deallocate(pars)
      !
   enddo
   !
   close(500)
   !
   ! Count number of structures
   !
   nstruc = 0
   !
   do ip = 1, npuv
      !
      okay = .false.
      !
      if (istruc(ip)==1) then
         !
         height = strucpars(1, ip)
         !
         if (subgrid) then
            if (height>subgrid_uv_zmin(ip)) then
               okay = .true.
            endif   
         else
            if (height>zbuvmx(ip)) then
               okay = .true.
            endif   
         endif   
      endif      
      !
      if (okay) then
         !
         nstruc = nstruc + 1
         !
      else
         !
         istruc(ip) = 0
         !
      endif   
      !
   enddo   
   !
   allocate(structure_uv_index(nstruc))
   allocate(structure_type(nstruc))
   allocate(structure_length(nstruc))
   allocate(structure_parameters(npars, nstruc))
   !
   nstruc = 0
   !
   do ip = 1, npuv
      !
      if (istruc(ip)==1) then
         !
         nstruc = nstruc + 1
         structure_type(nstruc)     = itype
         structure_uv_index(nstruc) = ip
         !
         iref  = uv_flags_iref(ip) ! refinement level
         idir  = uv_flags_dir(ip)
         !
         if (crsgeo) then
            !
            ! Geographic coordinate system
            !
            if (idir==0) then
               !
               ! U point
               !
               dxs = dyrm(iref)
               !
            else
               !
               ! V point
               !
               dxs = 1.0 / dxminv(ip)
               !
            endif   
            !
         else   
            !
            ! Projected coordinate system
            !
            if (idir==0) then
               !
               ! U point
               !
               dxs  = dyrm(iref)
               !
            else
               !
               ! V point
               !
               dxs  = dxrm(iref)
               !
            endif
            !
         endif
         !         
         structure_length(nstruc)   = lngth(ip)/dxs ! Turn into structure length relative to grid spacing
         !
         do ipar = 1, npars
            structure_parameters(ipar, nstruc) = strucpars(ipar, ip)
         enddo
         !
         kcuv(ip) = 3
         !
      endif   
      !
   enddo   
   !
   nrstructures = nstruc
   !
   write(logstr,'(a,i0,a)')'Info    : ', nrstructures, ' structure u/v points found'
   call write_log(logstr, 0)
   !
   end subroutine

   
   
   subroutine read_thin_dams()
   !
   ! Reads thd file
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   integer ip, nm, nthd, nstruc, nrows, ncols, irow, stat, ithd, nr_points, total_nr_points, iuv, indx
   !
   real      :: dst, dstmin, xxx, yyy, dstx, dsty
   real      :: dummy
   character :: cdummy
   logical   :: ok
   !
   real*4,  dimension(:),   allocatable :: xthd
   real*4,  dimension(:),   allocatable :: ythd
   real*4,  dimension(2)                :: xp
   real*4,  dimension(2)                :: yp
   integer, dimension(:),   allocatable :: uv_indices   
   integer, dimension(:),   allocatable :: vertices
   integer*1, dimension(:), allocatable  :: istruc  
   !
   ! Read thin dams file
   !
   nthd = 0
   total_nr_points = 0
   !
   allocate(istruc(npuv))
   istruc = 0   
   !
   if (thdfile(1:4) /= 'none') then
      !
      ok = check_file_exists(thdfile, 'Thin dams thd file', .true.)
      !
      ! First count number of polylines
      !
      open(500, file=trim(thdfile))
      do while(.true.)
         read(500,*,iostat = stat)cdummy
         if (stat<0) exit
         read(500,*,iostat = stat)nrows,ncols
         if (stat<0) exit
         nthd = nthd + 1
         do irow = 1, nrows
            read(500,*)dummy
         enddo
      enddo
      rewind(500)
      !
      ! Loop through polylines
      !
      write(logstr,'(a,a,a,i0,a)')'Info    : reading ', trim(thdfile), ' (', nthd, ' thin dams found)'
      call write_log(logstr, 0)
      !
      do ithd = 1, nthd
         !
         read(500,*,iostat = stat)cdummy
         if (stat<0) exit
         read(500,*,iostat = stat)nrows,ncols
         if (stat<0) exit
         allocate(xthd(nrows))
         allocate(ythd(nrows))
         do irow = 1, nrows
            read(500,*)xthd(irow),ythd(irow)
         enddo
         !
         call find_uv_points_intersected_by_polyline(uv_indices, vertices, nr_points, xthd, ythd, nrows)
         !
         ! Now loop through uv points to get the length and determine parameters
         !
         do iuv = 1, nr_points
            !
            indx = uv_indices(iuv)
            kcuv(indx) = 0
            istruc(indx) = 1            
            !
         enddo
         !
         deallocate(xthd)
         deallocate(ythd)
         !
         total_nr_points = total_nr_points + nr_points
         !    
      enddo
      !
      close(500)
      !
      nrthindams = total_nr_points
      !      
      allocate(thindam_uv_index(nrthindams))
      !
      nstruc = 0
      !
      do ip = 1, npuv
         !
         if (istruc(ip)==1) then
            !
            nstruc = nstruc + 1
            thindam_uv_index(nstruc) = ip
            !
         endif
      enddo      
      !
      write(logstr,'(a,i0,a)')'Info    : ', total_nr_points,' structure u/v points found'
      call write_log(logstr, 0)
      !
   endif   
   !
   end subroutine

   
   subroutine give_structure_information(struc_info)
   !
   ! Subroutine to provide structure information to output
   use sfincs_data
   !
   implicit none
   !
   integer                      :: n, m, istruc, nm, nmu, ip
   real*4, dimension(:,:), allocatable   :: struc_info
   real*4, dimension(:,:), allocatable :: xg
   real*4, dimension(:,:), allocatable :: yg
   
   ! Make empty struc_info
   allocate(struc_info(nrstructures,3))
   struc_info = 0.0
   ! 
   ! Define the grid 
   !
   allocate(xg(mmax + 1, nmax + 1))
   allocate(yg(mmax + 1, nmax + 1))
   !
   do n = 1, nmax + 1
       do m = 1, mmax + 1
           xg(m, n) = x0 + cosrot*(1.0*(m - 1))*dx - sinrot*(1.0*(n - 1))*dy
           yg(m, n) = y0 + sinrot*(1.0*(m - 1))*dx + cosrot*(1.0*(n - 1))*dy
       enddo
   enddo
   !
   ! Get coordinates and height
   !
   do istruc = 1, nrstructures
       !
       ! Get index
       ip       = structure_uv_index(istruc)
       nmu      = uv_index_z_nmu(ip)
       nm       = uv_index_z_nm(ip)
       ! 
       ! Get coordinates of face
       struc_info(istruc,1) = z_xz(nmu)*0.5 + z_xz(nm)*0.5
       struc_info(istruc,2) = z_yz(nmu)*0.5 + z_yz(nm)*0.5
       !
       ! Save height
       struc_info(istruc,3) = structure_parameters(1, istruc)
   enddo
   !
   end subroutine
   
   
   subroutine give_thindam_information(struc_info)
   !
   ! Subroutine to provide structure information to output
   use sfincs_data
   !
   implicit none
   !
   integer                      :: n, m, istruc, nm, nmu, ip
   real*4, dimension(:,:), allocatable   :: struc_info
   real*4, dimension(:,:), allocatable :: xg
   real*4, dimension(:,:), allocatable :: yg
   
   ! Make empty struc_info
   allocate(struc_info(nrthindams,2))
   struc_info = 0.0
   ! 
   ! Define the grid 
   !
   allocate(xg(mmax + 1, nmax + 1))
   allocate(yg(mmax + 1, nmax + 1))
   !
   do n = 1, nmax + 1
       do m = 1, mmax + 1
           xg(m, n) = x0 + cosrot*(1.0*(m - 1))*dx - sinrot*(1.0*(n - 1))*dy
           yg(m, n) = y0 + sinrot*(1.0*(m - 1))*dx + cosrot*(1.0*(n - 1))*dy
       enddo
   enddo
   !
   ! Get coordinates
   !
   do istruc = 1, nrthindams
       !
       ! Get index
       ip       = thindam_uv_index(istruc)
       nmu      = uv_index_z_nmu(ip)
       nm       = uv_index_z_nm(ip)
       ! 
       ! Get coordinates of face
       struc_info(istruc,1) = z_xz(nmu)*0.5 + z_xz(nm)*0.5
       struc_info(istruc,2) = z_yz(nmu)*0.5 + z_yz(nm)*0.5
       !
   enddo
   !
   end subroutine

   
   subroutine compute_fluxes_over_structures()
   !
   ! Computes fluxes over structures (THIS HAS TO BE SERIOUSLY IMPROVED!!!)
   !
   use sfincs_data
   use sfincs_timers
!   use quadtree
   !
   implicit none
   !
   integer                       :: ip
   integer*4                     :: nm
   integer*4                     :: nmu
   integer                       :: istruc
   integer                       :: ikf
   integer                       :: idir
   !
   real*4                       :: zsnm
   real*4                       :: zsnmu
   real*4                       :: cweir
   real*4                       :: Cd
   real*4                       :: m
   real*4                       :: h1
   real*4                       :: h2
   real*4                       :: qstruc
   !
   call timer_start('structures')
   !
   !$acc parallel, present(zs, q, uv, structure_uv_index, uv_index_z_nm, uv_index_z_nmu, structure_parameters, structure_type, structure_length)
   !$acc loop independent gang vector
   do istruc = 1, nrstructures
      !
      ip     = structure_uv_index(istruc)
      !
      q(ip)  = 0.0
      uv(ip) = 0.0
      !
      nm  = uv_index_z_nm(ip)
      nmu = uv_index_z_nmu(ip)
      !
      zsnm  = zs(nm)
      zsnmu = zs(nmu)      
      !
      if (zsnm<structure_parameters(1, istruc) .and. zsnmu<structure_parameters(1, istruc)) then         
         !
         cycle ! No flow over structure
         !
      endif   
      !
      select case(structure_type(istruc))
         case(1)
            !
            ! Broad-crested weir
            !
            cweir = 1.7049
            m = 0.0 ! Modular limit ...
            Cd = structure_parameters(2, istruc)
            !
            ikf   = 0
            !
            ikf   = 1  ! use flux now
            !            
            if  (zsnm>zsnmu) then
               idir = 1
               h1 = zsnm  - structure_parameters(1, istruc)
               h2 = zsnmu - structure_parameters(1, istruc)
            else
               idir = -1
               h1 = zsnmu - structure_parameters(1, istruc)
               h2 = zsnm  - structure_parameters(1, istruc)
            endif
            !
            if (h2 > 2.0 / 3.0 * h1) then
               !
               ! fully submerged
               !
               qstruc = Cd*h2*sqrt(2.0 * 9.81 * (h1 - h2)/(1.0 - m))
               !
            else
               !
               ! free flow
               !
               qstruc = Cd*cweir*h1**1.5
               !
            endif 
            !
         case(2)
         case(3)
      end select
      !         
      qstruc = qstruc * structure_length(istruc)
      !
      q(ip)  = qstruc*idir ! Add relaxation here !!!
      !
   enddo
   !$acc end parallel
   !
   call timer_stop('structures')
   !
   end subroutine
   !
   !
   subroutine read_flow_blocking_file()
      !
      ! Reads the flow-blocking features file (elevated roads, railways and levees
      ! that the grid does not resolve), snaps the polylines to the uv points and
      ! raises the subgrid uv tables (or zbuvmx without subgrid) to the crest level.
      ! The file has the weir layout: name, "nrows ncols", rows of x y crest [manning].
      !
      ! Called from: sfincs_structures -> read_structures
      !
      use sfincs_data
      use sfincs_subgrid
      use quadtree
      !
      implicit none
      !
      integer       :: ip
      integer       :: nm
      integer       :: nmu
      integer       :: irow
      integer       :: ipol
      integer       :: iuv
      integer       :: npol
      integer       :: nrows
      integer       :: ncols
      integer       :: nr_points
      integer       :: nraised
      integer       :: stat
      real*4        :: dummy
      real*4        :: xuv
      real*4        :: yuv
      real*4        :: dst1
      real*4        :: dst2
      real*4        :: wfac
      real*4        :: zcrest
      integer       :: nskip_zmin
      logical       :: okay
      character*256 :: cdummy
      !
      real*4,  dimension(:), allocatable :: xpol
      real*4,  dimension(:), allocatable :: ypol
      real*4,  dimension(:), allocatable :: zpol
      real*4,  dimension(:), allocatable :: manning_pol
      real*4,  dimension(:), allocatable :: uv_zcrest
      real*4,  dimension(:), allocatable :: uv_manning
      integer, dimension(:), allocatable :: uv_indices
      integer, dimension(:), allocatable :: vertices
      !
      write(logstr,'(a)')'Info    : reading flow-blocking features file'
      call write_log(logstr, 0)
      !
      okay = check_file_exists(flow_blocking_file, 'Flow-blocking features file', .true.)
      !
      allocate(uv_zcrest(npuv))
      allocate(uv_manning(npuv))
      !
      uv_zcrest  = -1.0e10
      uv_manning = flow_blocking_manning
      !
      ! Count polylines
      !
      npol = 0
      !
      open(500, file=trim(flow_blocking_file))
      !
      do while (.true.)
         !
         read(500,*,iostat=stat)cdummy
         if (stat<0) exit
         read(500,*,iostat=stat)nrows,ncols
         if (stat<0) exit
         npol = npol + 1
         !
         do irow = 1, nrows
            !
            read(500,*)dummy
            !
         enddo
         !
      enddo
      !
      rewind(500)
      !
      ! Loop through polylines and find the highest crest at each uv point
      !
      nskip_zmin = 0
      !
      do ipol = 1, npol
         !
         read(500,*,iostat=stat)cdummy
         read(500,*,iostat=stat)nrows,ncols
         if (stat<0) exit
         !
         allocate(xpol(nrows))
         allocate(ypol(nrows))
         allocate(zpol(nrows))
         allocate(manning_pol(nrows))
         !
         manning_pol = flow_blocking_manning
         !
         do irow = 1, nrows
            !
            if (ncols>=4) then
               !
               read(500,*)xpol(irow), ypol(irow), zpol(irow), manning_pol(irow)
               !
            else
               !
               read(500,*)xpol(irow), ypol(irow), zpol(irow)
               !
            endif
            !
         enddo
         !
         call find_uv_points_intersected_by_polyline(uv_indices, vertices, nr_points, xpol, ypol, nrows)
         !
         do iuv = 1, nr_points
            !
            ip   = uv_indices(iuv)
            irow = vertices(iuv)
            nm   = uv_index_z_nm(ip)
            nmu  = uv_index_z_nmu(ip)
            !
            ! Crest level interpolated along the segment to the face centre
            !
            xuv  = 0.5*(z_xz(nm) + z_xz(nmu))
            yuv  = 0.5*(z_yz(nm) + z_yz(nmu))
            dst1 = sqrt((xuv - xpol(irow))**2 + (yuv - ypol(irow))**2)
            dst2 = sqrt((xuv - xpol(irow + 1))**2 + (yuv - ypol(irow + 1))**2)
            wfac = dst1 / max(dst1 + dst2, 1.0e-6)
            !
            zcrest = zpol(irow)*(1.0 - wfac) + zpol(irow + 1)*wfac
            !
            if (zcrest>uv_zcrest(ip)) then
               !
               uv_zcrest(ip)  = zcrest
               uv_manning(ip) = manning_pol(irow)*(1.0 - wfac) + manning_pol(irow + 1)*wfac
               !
            endif
            !
         enddo
         !
         deallocate(xpol)
         deallocate(ypol)
         deallocate(zpol)
         deallocate(manning_pol)
         !
      enddo
      !
      close(500)
      !
      ! Raise the faces whose tables sit below the crest
      !
      nraised = 0
      !
      do ip = 1, npuv
         !
         if (uv_zcrest(ip)>-1.0e9) then
            !
            if (subgrid) then
               !
               if (uv_zcrest(ip)>subgrid_uv_zmin(ip)) then
                  !
                  ! A face with a table bed level below the given minimum is left
                  ! alone (channels under bridge decks left in the DTM)
                  !
                  if (subgrid_uv_zmin(ip)<flow_blocking_min_zmin) then
                     !
                     nskip_zmin = nskip_zmin + 1
                     !
                  else
                     !
                     call set_subgrid_uv_blocking(ip, uv_zcrest(ip), uv_manning(ip))
                     nraised = nraised + 1
                     !
                  endif
                  !
               endif
               !
            else
               !
               if (uv_zcrest(ip) + huthresh>zbuvmx(ip)) then
                  !
                  if (zbuvmx(ip)<flow_blocking_min_zmin) then
                     !
                     nskip_zmin = nskip_zmin + 1
                     !
                  else
                     !
                     zbuvmx(ip) = uv_zcrest(ip) + huthresh
                     nraised    = nraised + 1
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
      nr_flow_blocking_points = nraised
      !
      write(logstr,'(a,i0,a,i0,a)')'Info    : ', npol, ' flow-blocking polylines read, ', nraised, ' u/v points raised to crest level'
      call write_log(logstr, 0)
      write(logstr,'(a,i0,a)')'Info    : flow-blocking u/v points skipped because the bed level is below flow_blocking_min_zmin : ', nskip_zmin
      call write_log(logstr, 0)
      !
      deallocate(uv_zcrest)
      deallocate(uv_manning)
      !
   end subroutine read_flow_blocking_file
   !
end module
