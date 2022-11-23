module quadtree

   integer*4                                       :: quadtree_nr_points
   integer*1                                       :: quadtree_nr_levels
   real*4                                          :: quadtree_x0
   real*4                                          :: quadtree_y0
   real*4                                          :: quadtree_dx
   real*4                                          :: quadtree_dy
   integer*4                                       :: quadtree_nmax
   integer*4                                       :: quadtree_mmax
   real*4                                          :: quadtree_rotation
   real*4                                          :: quadtree_cosrot
   real*4                                          :: quadtree_sinrot
   integer*1,          dimension(:),   allocatable :: quadtree_level
   integer*1,          dimension(:),   allocatable :: quadtree_md
   integer*4,          dimension(:),   allocatable :: quadtree_md1
   integer*4,          dimension(:),   allocatable :: quadtree_md2
   integer*1,          dimension(:),   allocatable :: quadtree_mu
   integer*4,          dimension(:),   allocatable :: quadtree_mu1
   integer*4,          dimension(:),   allocatable :: quadtree_mu2
   integer*1,          dimension(:),   allocatable :: quadtree_nd
   integer*4,          dimension(:),   allocatable :: quadtree_nd1
   integer*4,          dimension(:),   allocatable :: quadtree_nd2
   integer*1,          dimension(:),   allocatable :: quadtree_nu
   integer*4,          dimension(:),   allocatable :: quadtree_nu1
   integer*4,          dimension(:),   allocatable :: quadtree_nu2
   integer*4,          dimension(:),   allocatable :: quadtree_n
   integer*4,          dimension(:),   allocatable :: quadtree_m
   integer*1,          dimension(:),   allocatable :: quadtree_n_oddeven
   integer*1,          dimension(:),   allocatable :: quadtree_m_oddeven
   real*4,             dimension(:),   allocatable :: quadtree_xz
   real*4,             dimension(:),   allocatable :: quadtree_yz
   real*4,             dimension(:),   allocatable :: quadtree_zz
   integer*4,          dimension(:),   allocatable :: quadtree_nm_indices
   integer*4,          dimension(:),   allocatable :: quadtree_first_point_per_level      
   integer*4,          dimension(:),   allocatable :: quadtree_last_point_per_level      
   real*4,             dimension(:),   allocatable :: quadtree_dxr
   real*4,             dimension(:),   allocatable :: quadtree_dyr
   
contains

   subroutine quadtree_read_file(qtrfile)
   !
   ! Reads quadtree file
   !
   implicit none
   !
   character*256, intent(in)                       :: qtrfile
   !
   real*4,             dimension(:),   allocatable :: dxr
   real*4,             dimension(:),   allocatable :: dyr
   !
   integer*1 :: iversion
   integer :: np, nm, n, m, iref, nmx, ireflast, ip, iepsg
   real*4  :: pi, cosrot, sinrot
   !
   pi = 3.141592653589793
   !
   ! Read quadtree file (first time, only read number of active points)
   !
   open(unit = 500, file = trim(qtrfile), form = 'unformatted', access = 'stream')
   read(500)iversion
   read(500)iepsg
   read(500)np
   read(500)quadtree_nr_levels
   close(500)
   !
   allocate(quadtree_level(np))
   allocate(quadtree_md(np))
   allocate(quadtree_md1(np))
   allocate(quadtree_md2(np))
   allocate(quadtree_mu(np))
   allocate(quadtree_mu1(np))
   allocate(quadtree_mu2(np))
   allocate(quadtree_nd(np))
   allocate(quadtree_nd1(np))
   allocate(quadtree_nd2(np))
   allocate(quadtree_nu(np))
   allocate(quadtree_nu1(np))
   allocate(quadtree_nu2(np))
   allocate(quadtree_n(np))
   allocate(quadtree_m(np))
   allocate(quadtree_n_oddeven(np))
   allocate(quadtree_m_oddeven(np))
   allocate(quadtree_xz(np))
   allocate(quadtree_yz(np))
   allocate(quadtree_zz(np))
   !
   ! Read whole quadtree file
   !
   open(unit = 500, file = trim(qtrfile), form = 'unformatted', access = 'stream')
   !
   read(500)iversion
   read(500)iepsg
   read(500)quadtree_nr_points
   read(500)quadtree_nr_levels
   quadtree_nr_levels = quadtree_nr_levels
   !
   read(500)quadtree_x0
   read(500)quadtree_y0
   read(500)quadtree_dx
   read(500)quadtree_dy
   read(500)quadtree_rotation
   !
   read(500)quadtree_level
   read(500)quadtree_n
   read(500)quadtree_m
   !
   read(500)quadtree_nu
   read(500)quadtree_nu1
   read(500)quadtree_nu2
   read(500)quadtree_mu
   read(500)quadtree_mu1
   read(500)quadtree_mu2
   read(500)quadtree_nd
   read(500)quadtree_nd1
   read(500)quadtree_nd2
   read(500)quadtree_md
   read(500)quadtree_md1
   read(500)quadtree_md2
   !
   read(500)quadtree_zz
   !
   close(500)
   !
   quadtree_rotation = quadtree_rotation*pi/180
   cosrot            = cos(quadtree_rotation)
   sinrot            = sin(quadtree_rotation)
   quadtree_cosrot   = cosrot
   quadtree_sinrot   = sinrot
   !
   allocate(dxr(quadtree_nr_levels))
   allocate(dyr(quadtree_nr_levels))
   allocate(quadtree_dxr(quadtree_nr_levels))
   allocate(quadtree_dyr(quadtree_nr_levels))
   !
   do iref = 1, quadtree_nr_levels
      dxr(iref) = quadtree_dx / 2**(iref - 1)
      dyr(iref) = quadtree_dy / 2**(iref - 1)
   enddo   
   !
   quadtree_dxr = dxr
   quadtree_dyr = dyr
   !
   quadtree_nmax = 0
   quadtree_mmax = 0
   !
   do nm = 1, quadtree_nr_points
      !
      n    = quadtree_n(nm)
      m    = quadtree_m(nm)
      iref = quadtree_level(nm)
      !
      quadtree_nmax = max(quadtree_nmax, int( (1.0*(n - 1) + 0.01) / (2**(iref - 1))) + 2)
      quadtree_mmax = max(quadtree_mmax, int( (1.0*(m - 1) + 0.01) / (2**(iref - 1))) + 2)
      !
      quadtree_xz(nm) = quadtree_x0 + cosrot*(1.0*(m - 0.5))*dxr(iref) - sinrot*(1.0*(n - 0.5))*dyr(iref)
      quadtree_yz(nm) = quadtree_y0 + sinrot*(1.0*(m - 0.5))*dxr(iref) + cosrot*(1.0*(n - 0.5))*dyr(iref)
      !
      ! Check for odd or even
      !
      if (mod(n, 2)>0) then
         quadtree_n_oddeven(nm) = 1         
      else
         quadtree_n_oddeven(nm) = 2         
      endif   
      !
      if (mod(m, 2)>0) then
         quadtree_m_oddeven(nm) = 1         
      else
         quadtree_m_oddeven(nm) = 2         
      endif   
      !
   enddo
   !
   ! Make some arrays for easy searching
   !
   allocate(quadtree_nm_indices(quadtree_nr_points))
   allocate(quadtree_first_point_per_level(quadtree_nr_levels))
   allocate(quadtree_last_point_per_level(quadtree_nr_levels))
   quadtree_first_point_per_level = 0
   quadtree_last_point_per_level = 0
   quadtree_nm_indices = 0
   !
   ! First count
   !
   ireflast = 0
   !
   do ip = 1, quadtree_nr_points
      !
      iref = quadtree_level(ip)
      n    = quadtree_n(ip)
      m    = quadtree_m(ip)
      nmx  = quadtree_nmax*2**(iref - 1)
      nm   = (m - 1)*nmx + n
      !
      quadtree_nm_indices(ip) = nm
      !
      if (iref>ireflast) then
         !
         ! Found new level
         !
         quadtree_first_point_per_level(iref) = ip
         ireflast = iref
         !
      endif   
      !
      quadtree_last_point_per_level(iref) = ip
      !
   enddo   
   !
   end subroutine

   
   subroutine make_quadtree_from_indices(np, indices, nmax, mmax, x0, y0, dx, dy, rotation)
   !
   implicit none
   !
   integer*4 :: np
   integer*4 :: nmax
   integer*4 :: mmax
   real*4    :: x0 
   real*4    :: y0 
   real*4    :: dx 
   real*4    :: dy
   real*4    :: rotation
   integer*4, dimension(1:np), intent(in) :: indices   
   !
   integer*4,          dimension(:),   allocatable :: index_v_m
   integer*4,          dimension(:),   allocatable :: index_v_n
   integer*4,          dimension(:),   allocatable :: index_v_nm
   integer*4,          dimension(:),   allocatable :: index_v_nmu
   integer*4,          dimension(:),   allocatable :: index_v_num
   integer*4,          dimension(:),   allocatable :: index_v_nmd
   integer*4,          dimension(:),   allocatable :: index_v_ndm
   integer*4,          dimension(:,:), allocatable :: index_g_nm
   !
   integer*4  :: nm, n, m, ip, nmx, iref, ireflast
   real*4 :: cosrot
   real*4 :: sinrot
   !
   logical :: global
   !
   global = .false.
   !
   allocate(quadtree_level(np))
   allocate(quadtree_md(np))
   allocate(quadtree_md1(np))
   allocate(quadtree_md2(np))
   allocate(quadtree_mu(np))
   allocate(quadtree_mu1(np))
   allocate(quadtree_mu2(np))
   allocate(quadtree_nd(np))
   allocate(quadtree_nd1(np))
   allocate(quadtree_nd2(np))
   allocate(quadtree_nu(np))
   allocate(quadtree_nu1(np))
   allocate(quadtree_nu2(np))
   allocate(quadtree_n(np))
   allocate(quadtree_m(np))
   allocate(quadtree_n_oddeven(np))
   allocate(quadtree_m_oddeven(np))
   allocate(quadtree_xz(np))
   allocate(quadtree_yz(np))
   !
   allocate(index_v_m(np))
   allocate(index_v_n(np))
   allocate(index_v_nm(np))
   allocate(index_v_nmu(np))
   allocate(index_v_nmd(np))
   allocate(index_v_num(np))
   allocate(index_v_ndm(np))
   allocate(index_g_nm(nmax, mmax))
   !
   index_v_m   = 0
   index_v_n   = 0
   index_v_nm  = 0
   index_v_nmu = 0
   index_v_nmd = 0
   index_v_num = 0
   index_v_ndm = 0
   index_g_nm  = 0
   !
   do nm = 1, np
      !
      m = int((indices(nm) - 1)/nmax) + 1
      n = indices(nm) - (m - 1)*nmax
      !
      index_v_m(nm)    = m
      index_v_n(nm)    = n
      index_v_nm(nm)   = nm
      index_g_nm(n, m) = nm
      !
   enddo
   !
   ! Find neighboring cell indices
   !
   do nm = 1, np
      !
      m = index_v_m(nm)
      n = index_v_n(nm)
      !
      if (global) then
         !
         ! Let's try a global SFINCS model and see how that works out!
         !
         if (m==mmax - 1) then
            !
            ! Boundary point on the right
            !
            index_v_nmu(nm) = index_g_nm(n, 2)
            !
         else   
            !
            ! Normal point 
            !            
            index_v_nmu(nm) = index_g_nm(n, m + 1)
            !
         endif   
         !
         if (m==2) then
            !
            ! Boundary point on the left
            !
            index_v_nmd(nm) = index_g_nm(n, mmax - 1)
            !
         else
            !
            ! Normal point 
            !
            index_v_nmd(nm) = index_g_nm(n, m - 1)
            !
         endif   
         !            
         index_v_num(nm) = index_g_nm(n + 1, m)
         index_v_ndm(nm) = index_g_nm(n - 1, m)
         !
      else
         !
         if (m<mmax) then
            index_v_nmu(nm) = index_g_nm(n, m + 1)
         endif
         !
         if (m>1) then            
            index_v_nmd(nm) = index_g_nm(n, m - 1)
         endif
         !
         if (n<nmax) then
            index_v_num(nm) = index_g_nm(n + 1, m)
         endif
         !
         if (n>1) then
            index_v_ndm(nm) = index_g_nm(n - 1, m)
         endif
         !
      endif
      !
   enddo
   !
   quadtree_nr_points = np
   quadtree_nr_levels = 1
   quadtree_nmax      = nmax
   quadtree_mmax      = mmax
   quadtree_x0        = x0
   quadtree_y0        = y0
   quadtree_dx        = dx
   quadtree_dy        = dy
   quadtree_rotation  = rotation
   !
   cosrot             = cos(quadtree_rotation)
   sinrot             = sin(quadtree_rotation)
   quadtree_cosrot    = cosrot
   quadtree_sinrot    = sinrot
   !
   allocate(quadtree_dxr(1))
   allocate(quadtree_dyr(1))
   !
   quadtree_dxr(1) = dx
   quadtree_dyr(1) = dy
   !   
   do ip = 1, np
      quadtree_level(ip) = 1
      quadtree_md(ip)    = 0
      quadtree_md1(ip)   = index_v_nmd(ip)
      quadtree_md2(ip)   = 0
      quadtree_mu(ip)    = 0
      quadtree_mu1(ip)   = index_v_nmu(ip)
      quadtree_mu2(ip)   = 0
      quadtree_nd(ip)    = 0
      quadtree_nd1(ip)   = index_v_ndm(ip)
      quadtree_nd2(ip)   = 0
      quadtree_nu(ip)    = 0
      quadtree_nu1(ip)   = index_v_num(ip)
      quadtree_nu2(ip)   = 0
      quadtree_n(ip)     = index_v_n(ip)
      quadtree_m(ip)     = index_v_m(ip)
      quadtree_n_oddeven(ip) = 1
      quadtree_m_oddeven(ip) = 1      
      quadtree_xz(ip) = quadtree_x0 + cosrot*(1.0*index_v_m(ip) - 0.5)*dx - sinrot*(1.0*index_v_n(ip) - 0.5)*dy
      quadtree_yz(ip) = quadtree_y0 + sinrot*(1.0*index_v_m(ip) - 0.5)*dx + cosrot*(1.0*index_v_n(ip) - 0.5)*dy
      !
   enddo   
   !
   ! Make some arrays for easy searching
   !
   allocate(quadtree_nm_indices(quadtree_nr_points))
   allocate(quadtree_first_point_per_level(quadtree_nr_levels))
   allocate(quadtree_last_point_per_level(quadtree_nr_levels))
   quadtree_first_point_per_level = 0
   quadtree_last_point_per_level = 0
   quadtree_nm_indices = 0
   !
   ireflast = 0
   !
   do ip = 1, quadtree_nr_points
      !
      iref = quadtree_level(ip)
      n    = quadtree_n(ip)
      m    = quadtree_m(ip)
      nmx  = quadtree_nmax*2**(iref - 1)
      nm   = (m - 1)*nmx + n
      !
      quadtree_nm_indices(ip) = nm
      !
      if (iref>ireflast) then
         !
         ! Found new level
         !
         quadtree_first_point_per_level(iref) = ip
         ireflast = iref
         !
      endif   
      !
      quadtree_last_point_per_level(iref) = ip
      !
   enddo   
   !         
   end subroutine 

   subroutine find_cells_intersected_by_line(cell_indices,nr_cells,xpa, ypa, xpb, ypb)
   !
   ! Builds array cell_indices of all cells intersected by line segment ([xpa ypa],[xpb ypb])
   !
   use geometry
   !
   implicit none
   !
   real*4, intent(in)                                :: xpa, ypa, xpb, ypb
   integer*4, dimension(:), allocatable, intent(out) :: cell_indices
   integer :: nr_cells
   !   
   integer*4, dimension(:), allocatable :: n0, m0, n1, m1
   integer*4, dimension(:), allocatable :: ind
   !
   real*4    :: xtmpa, xtmpb, ytmpa, ytmpb   
   real*4    :: xx0, yy0, xx1, yy1   
   integer*4 :: na, nb, ma, mb, nmx, ip, iref, n, m
   logical   :: incell
   !
   allocate(n0(quadtree_nr_levels))
   allocate(n1(quadtree_nr_levels))
   allocate(m0(quadtree_nr_levels))
   allocate(m1(quadtree_nr_levels))
   !
   allocate(ind(quadtree_nr_points))
   ind = 0
   !
   ! Determine n0,m0,n1, and m1 for each level
   !
   xtmpa =   quadtree_cosrot*(xpa - quadtree_x0) + quadtree_sinrot*(ypa - quadtree_y0)
   ytmpa = - quadtree_sinrot*(xpa - quadtree_x0) + quadtree_cosrot*(ypa - quadtree_y0)
   xtmpb =   quadtree_cosrot*(xpb - quadtree_x0) + quadtree_sinrot*(ypb - quadtree_y0)
   ytmpb = - quadtree_sinrot*(xpb - quadtree_x0) + quadtree_cosrot*(ypb - quadtree_y0)
   !   
   do iref = quadtree_nr_levels, 1, -1
      !
      nmx = quadtree_nmax*2**(iref - 1)
      !
      na  = int(ytmpa/quadtree_dyr(iref)) + 1
      ma  = int(xtmpa/quadtree_dxr(iref)) + 1
      nb  = int(ytmpb/quadtree_dyr(iref)) + 1
      mb  = int(xtmpb/quadtree_dxr(iref)) + 1
      !
      n0(iref) = min(na, nb) - 2
      n1(iref) = max(na, nb) + 2
      m0(iref) = min(ma, mb) - 2
      m1(iref) = max(ma, mb) + 2
      !
   enddo   
   !
   nr_cells = 0
   !
   do ip = 1, quadtree_nr_points
      !
      iref = quadtree_level(ip)
      !
      n = quadtree_n(ip)
      m = quadtree_m(ip)
      !
      ! Check if cell can be intersected at all
      !
      if (n>=n0(iref) .and. n<=n1(iref) .and. m>=m0(iref) .and. m<=m1(iref)) then
         !
         xx0 = (m - 1)*quadtree_dxr(iref)
         yy0 = (n - 1)*quadtree_dyr(iref)
         xx1 = (m    )*quadtree_dxr(iref)
         yy1 = (n    )*quadtree_dyr(iref)
         !
         incell = check_line_through_square(xx0, yy0, xx1, yy1, xtmpa, ytmpa, xtmpb, ytmpb)
         !
         if (incell) then
            !
            nr_cells      = nr_cells + 1
            ind(nr_cells) = ip
            !
         endif   
         !
      endif   
      !
   enddo   
   !
   allocate(cell_indices(nr_cells))
   cell_indices = 0
   !
   do ip = 1, nr_cells
      cell_indices(ip) = ind(ip)
   enddo   
   !
   end subroutine

   
   
   function find_quadtree_cell(x, y) result (indx)
   !
   implicit none
   !
   real*4, intent(in)   :: x
   real*4, intent(in)   :: y
   integer              :: indx
   !
   real*4               :: xtmp
   real*4               :: ytmp
   !
   integer              :: iref
   integer              :: nmx
   integer              :: n, m, nm, i1, i2, j
   !
   indx = 0
   !
   xtmp =   quadtree_cosrot*(x - quadtree_x0) + quadtree_sinrot*(y - quadtree_y0)
   ytmp = - quadtree_sinrot*(x - quadtree_x0) + quadtree_cosrot*(y - quadtree_y0)
   !   
   do iref = quadtree_nr_levels, 1, -1
      !
      nmx = quadtree_nmax*2**(iref - 1)
      !
      n  = int(ytmp/quadtree_dyr(iref)) + 1
      m  = int(xtmp/quadtree_dxr(iref)) + 1
      nm = (m - 1)*nmx + n
      !
      ! Find nm index of this point
      !
      i1 = quadtree_first_point_per_level(iref)
      i2 = quadtree_last_point_per_level(iref)
      !
      if (nm>0) then      
         !
         indx = binary_search(quadtree_nm_indices(i1:i2), i2 - i1 + 1, nm)
         !
      endif
      !
      if (indx>0) then
         !
         ! Now check if this is the point we're actually interested in (could also be another point with another refinement level with the same nm index)
         !
         if (quadtree_level(i1 + indx - 1) == iref) then
            !
            ! This must be our point
            !
            indx = indx + i1 - 1
            exit
            !
         endif
         !
      endif
      !
   enddo
   !
   end function

   function find_quadtree_cell_by_index(n, m, iref) result (indx)
   !
   implicit none
   !
   integer              :: indx
   !
   integer              :: iref
   integer              :: nmx
   integer              :: n, m, nm, i1, i2, j
   !
   indx = 0
   !
   nmx = quadtree_nmax*2**(iref - 1)
   !
   nm = (m - 1)*nmx + n
   !
   ! Find nm index of this point
   !
   i1 = quadtree_first_point_per_level(iref)
   i2 = quadtree_last_point_per_level(iref)
   !
   if (nm>0) then      
      !
      indx = binary_search(quadtree_nm_indices(i1:i2), i2 - i1 + 1, nm)
      indx = indx + i1 - 1
      !
   endif
   !
   end function
   

   subroutine find_uv_points_intersected_by_polyline(uv_indices,vertices,nint,xpol,ypol,npol)
   !
   ! Builds array of uv-points that are intersected by polyline
   !
   use geometry
   use sfincs_data
   !
   implicit none
   !
   integer     :: npol
   real*4, dimension(npol),              intent(in)  :: xpol,ypol
   integer*4, dimension(:), allocatable, intent(out) :: uv_indices
   integer*4, dimension(:), allocatable, intent(out) :: vertices
   integer :: nint
   !   
   integer*4, dimension(:), allocatable :: ind
   !
   real*4    :: xtmpa, xtmpb, ytmpa, ytmpb   
   real*4    :: xpa, ypa, xpb, ypb
   real*4    :: xx0, yy0, xx1, yy1   
   real*4    :: xuv1, yuv1, xuv2, yuv2
   integer*4 :: na, nb, ma, mb, nmx, ip, iref, n, m, nm, nmu, ipol, m0, n0, m1, n1
   logical   :: iok
   !
   allocate(ind(npuv))
   !
   ind  = 0
   nint = 0
   !
   ! Loop through segments along polyline
   !
   do ipol = 1, npol - 1
      !
      xpa = xpol(ipol)
      ypa = ypol(ipol)
      xpb = xpol(ipol + 1)
      ypb = ypol(ipol + 1)
      !
      ! First rotate line segement to quadtree coordinate system in order to find nearby cells
      !
      xtmpa =   quadtree_cosrot*(xpa - x0) + quadtree_sinrot*(ypa - y0)
      ytmpa = - quadtree_sinrot*(xpa - x0) + quadtree_cosrot*(ypa - y0)
      xtmpb =   quadtree_cosrot*(xpb - x0) + quadtree_sinrot*(ypb - y0)
      ytmpb = - quadtree_sinrot*(xpb - x0) + quadtree_cosrot*(ypb - y0)
      !
      do iref = nref, 1, -1
         !
         ! Determine which cell may be crossed by line segment
         !
         nmx = quadtree_nmax*2**(iref - 1)
         !
         na  = int(ytmpa/quadtree_dyr(iref)) + 1
         ma  = int(xtmpa/quadtree_dxr(iref)) + 1
         nb  = int(ytmpb/quadtree_dyr(iref)) + 1
         mb  = int(xtmpb/quadtree_dxr(iref)) + 1
         !
         n0 = min(na, nb) - 1
         n1 = max(na, nb) + 1
         m0 = min(ma, mb) - 1
         m1 = max(ma, mb) + 1
         !
         do m = m0, m1
            do n = n0, n1            
               !
               nm = find_quadtree_cell_by_index(n, m, iref)
               !
               if (nm==0) then
                  cycle
               endif   
               !
               ! Right (same level or coarser)
               !
               if (z_index_uv_mu1(nm)>0) then
                  !
                  ip  = z_index_uv_mu1(nm)
                  nmu = uv_index_z_nmu(ip)
                  iok = cross(xpa,ypa,xpb,ypb,z_xz(nm),z_yz(nm),z_xz(nmu),z_yz(nmu))
                  !
                  if (iok) then
                     ind(ip) = ipol
                     nint = nint + 1
                  endif
                  !
               endif   
               !
               ! Right (finer)
               !
               if (z_index_uv_mu2(nm)>0) then
                  !
                  ip  = z_index_uv_mu2(nm)
                  nmu = uv_index_z_nmu(ip)
                  iok = cross(xpa,ypa,xpb,ypb,z_xz(nm),z_yz(nm),z_xz(nmu),z_yz(nmu))
                  !
                  if (iok) then
                     ind(ip) = ipol
                     nint = nint + 1
                  endif
                  !
               endif   
               !
               ! Top (same level or coarser)
               !
               if (z_index_uv_nu1(nm)>0) then
                  !
                  ip  = z_index_uv_nu1(nm)
                  nmu = uv_index_z_nmu(ip)
                  iok = cross(xpa,ypa,xpb,ypb,z_xz(nm),z_yz(nm),z_xz(nmu),z_yz(nmu))
                  !
                  if (iok) then
                     ind(ip) = ipol
                     nint = nint + 1
                  endif
                  !
               endif   
               !
               ! Top (finer)
               !
               if (z_index_uv_nu2(nm)>0) then
                  !
                  ip  = z_index_uv_nu2(nm)
                  nmu = uv_index_z_nmu(ip)
                  iok = cross(xpa,ypa,xpb,ypb,z_xz(nm),z_yz(nm),z_xz(nmu),z_yz(nmu))
                  !
                  if (iok) then
                     ind(ip) = ipol
                     nint = nint + 1
                  endif
                  !
               endif   
               !
            enddo   
         enddo   
      enddo   
   enddo   
   !
   allocate(uv_indices(nint))
   allocate(vertices(nint))
   !
   nint = 0
   !
   do ip = 1, npuv
      !
      if (ind(ip)>0) then
         !
         nint = nint + 1
         uv_indices(nint) = ip
         vertices(nint)   = ind(ip)
         !
      endif
   enddo
   !
   end subroutine   
      
   function binary_search(x,n,val) result (indx)
   !
   implicit none
   !
   integer, intent(in)  :: n
   integer, intent(in)  :: x(n)
   integer, intent(in)  :: val
   integer              :: indx
   !
   integer              :: low
   integer              :: high
   integer              :: middle
   !
   logical              :: found
   !
   found = .false.
   low  = 0
   high = n
   indx = 0
   !
   do while(low <= high .and. .not.found)
      !
      middle = (low + high)/2
      ! 
      if (middle==0) then
         !
         low = middle + 1
         !
      elseif (val == x(middle) ) then
         !
         found = .true.
         indx  = middle
         !
      elseif (val < x(middle)) then
         !
         high = middle - 1
         !
      else
         !
         low = middle + 1
         !
      endif
      !
   enddo
   !
   end function
   
   
   subroutine quadtree_deallocate()
   !
   deallocate(quadtree_level)
   deallocate(quadtree_md)
   deallocate(quadtree_md1)
   deallocate(quadtree_md2)
   deallocate(quadtree_mu)
   deallocate(quadtree_mu1)
   deallocate(quadtree_mu2)
   deallocate(quadtree_nd)
   deallocate(quadtree_nd1)
   deallocate(quadtree_nd2)
   deallocate(quadtree_nu)
   deallocate(quadtree_nu1)
   deallocate(quadtree_nu2)
   deallocate(quadtree_n)
   deallocate(quadtree_m)
   deallocate(quadtree_n_oddeven)
   deallocate(quadtree_m_oddeven)
   deallocate(quadtree_xz)
   deallocate(quadtree_yz)
   deallocate(quadtree_nm_indices)
   deallocate(quadtree_first_point_per_level)
   deallocate(quadtree_last_point_per_level)
   deallocate(quadtree_dxr)
   deallocate(quadtree_dyr)
   !
   end subroutine   

end module
