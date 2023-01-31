   module sfincs_wavemaker

   contains

   subroutine read_wavemaker_polylines()
   !
   ! Reads polylines files
   ! Computes cross sections
   ! Determines interpolation weights
   ! Interpolates bathymetry
   !
   use sfincs_data
   use quadtree
!   use sfincs_interpolation
   !
   implicit none
   !
   integer*4 :: npol
   integer*4 :: nrows
   integer*4 :: ncols
   integer*4 :: stat
   integer*4 :: ipol
   integer*4 :: irow   
!   integer*4 :: ncross
!   integer*4 :: ncs
!   integer*4 :: ic
!   integer*4 :: icross
!   integer*4 :: ipc
   integer*4 :: i
   integer*4 :: j
   integer*4 :: m
   integer*4 :: n
   integer*4 :: nm
   integer*4 :: nmu
!   integer*4 :: nrpts
   integer*4 :: ind
!   integer*4 :: indwm
   integer*4 :: nrw
   integer*4 :: ip
   integer*4 :: iwm
   integer*4 :: iz
   integer*4 :: nr_cells
   integer*4 :: nok
   integer*4 :: ifreq
   integer*4 :: itb
   !
   real*4    :: dummy
   real*4    :: phip
   !
   real*4    :: r
!   real*4    :: dxcross
!   real*4    :: dst
!   real*4    :: xxx
!   real*4    :: yyy
!   real*4    :: lsec
!   real*4    :: phi
   !
   character :: cdummy
   character*256 :: filename
   !
   real*4,    dimension(:),     allocatable :: xpol
   real*4,    dimension(:),     allocatable :: ypol
   real*4,    dimension(:),     allocatable :: phi
!   integer*4, dimension(:),     allocatable :: npts
   !
!   integer*4, dimension(:,:), allocatable   :: icode
!   integer*4, dimension(:,:), allocatable   :: idummy
!   real*8,    dimension(:,:), allocatable   :: xzd
!   real*8,    dimension(:,:), allocatable   :: yzd
   !
   integer*4, dimension(:),     allocatable :: cell_indices
   integer*4, dimension(:),     allocatable :: indwm
   !
   logical :: iok
   !
   integer ib1, ib2, ib, ic, nmb
   !
   real x, y, dst1, dst2, dst
   !
   allocate(indwm(np))
   allocate(phi(np))
   !
   indwm = 0
   phi   = 0.0
   !
   nrw = 0
   !
   write(*,*)'Reading wavemaker polyline file ...'
   !
   ! Loop through all polylines
   !
   open(500, file=trim(wvmfile))
   !
   do while(.true.)
      !
      read(500,*,iostat = stat)cdummy
      read(500,*,iostat = stat)nrows,ncols
      if (stat<0) exit
      allocate(xpol(nrows))
      allocate(ypol(nrows))
      do irow = 1, nrows
         read(500,*)xpol(irow),ypol(irow)
      enddo
      !
      do irow = 1, nrows - 1
         !
         ! Determine angle with respect to grid orientation
         !
         phip = atan2(ypol(irow + 1) - ypol(irow), xpol(irow + 1) - xpol(irow)) + 0.5*pi
         phip = phip - rotation
         if (phip>=2*pi) phip = phip - 2*pi
         if (phip<0.0)   phip = phip + 2*pi
         !
         call find_cells_intersected_by_line(cell_indices, nr_cells, xpol(irow), ypol(irow), xpol(irow + 1), ypol(irow + 1))
         !
         do j = 1, nr_cells
            !
            ip = index_sfincs_in_quadtree(cell_indices(j))
            !
            if (indwm(ip) == 0) then
               !
               indwm(ip) = 1 ! set temporary flag to 1
               phi(ip)   = phip
               nrw       = nrw + 1
               !
            endif   
            !
         enddo   
         !
      enddo   
      ! 
      deallocate(xpol)
      deallocate(ypol)
      !
   enddo
   !
   close(500)
   !
   ! Now get rid of cells that have neighbor closer to shore that is also a wavemaker point
   !
   nok = 0
   !
   do ip = 1, np
      !
      ! Check if these cells have neighbor closer to shore that is also a wavemaker point
      !
      if (indwm(ip)==1) then
         !
         iok = .false.
         !
         if (phi(ip)>=0.0 .and. phi(ip)<0.5*pi) then
            !
            ! Check right and above
            !
            nmu = z_index_uv_mu1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (indwm(iz)==0) then
                  !
                  ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                  !
                  iok = .true.
                  !
               endif   
            endif
            !
            if (.not. iok) then
               !
               nmu = z_index_uv_mu2(ip)    ! index of 2nd uv neighbor to the right
               !
               if (nmu>0) then
                  !
                  iz = uv_index_z_nmu(nmu)
                  !
                  if (indwm(iz)==0) then
                     !
                     ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                     !
                     iok = .true.
                     !
                  endif   
               endif
            endif
            !
            if (.not. iok) then
               !
               nmu = z_index_uv_nu1(ip)    ! index of 1st uv neighbor above
               !
               if (nmu>0) then
                  !
                  iz = uv_index_z_nmu(nmu)
                  !
                  if (indwm(iz)==0) then
                     !
                     ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                     !
                     iok = .true.
                     !
                  endif   
               endif
            endif
            !
            if (.not. iok) then
               !
               nmu = z_index_uv_nu2(ip)    ! index of 2nd uv neighbor above
               !
               if (nmu>0) then
                  !
                  iz = uv_index_z_nmu(nmu)
                  !
                  if (indwm(iz)==0) then
                     !
                     ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                     !
                     iok = .true.
                     !
                  endif   
               endif
            endif
            !
         elseif (phi(ip)>=0.5*pi .and. phi(ip)<pi) then
            !
            ! Check left and above
            !
            nmu = z_index_uv_md1(ip)    ! index of 1st uv neighbor to the left
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (indwm(iz)==0) then
                  !
                  ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                  !
                  iok = .true.
                  !
               endif   
            endif
            !
            if (.not. iok) then
               !
               nmu = z_index_uv_md2(ip)    ! index of 2nd uv neighbor to the left
               !
               if (nmu>0) then
                  !
                  iz = uv_index_z_nm(nmu)
                  !
                  if (indwm(iz)==0) then
                     !
                     ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                     !
                     iok = .true.
                     !
                  endif   
               endif
            endif
            !
            if (.not. iok) then
               !
               nmu = z_index_uv_nu1(ip)    ! index of 1st uv neighbor above
               !
               if (nmu>0) then
                  !
                  iz = uv_index_z_nmu(nmu)
                  !
                  if (indwm(iz)==0) then
                     !
                     ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                     !
                     iok = .true.
                     !
                  endif   
               endif
            endif
            !
            if (.not. iok) then
               !
               nmu = z_index_uv_nu2(ip)    ! index of 2nd uv neighbor above
               !
               if (nmu>0) then
                  !
                  iz = uv_index_z_nmu(nmu)
                  !
                  if (indwm(iz)==0) then
                     !
                     ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                     !
                     iok = .true.
                     !
                  endif   
               endif
            endif
            !
         elseif (phi(ip)>=1.0*pi .and. phi(ip)<1.5*pi) then
            !
            ! Check left and below
            !
            nmu = z_index_uv_md1(ip)    ! index of 1st uv neighbor to the left
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(ip)
               !
               if (indwm(iz)==0) then
                  !
                  ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                  !
                  iok = .true.
                  !
               endif   
            endif
            !
            if (.not. iok) then
               !
               nmu = z_index_uv_md2(ip)    ! index of 2nd uv neighbor to the left
               !
               if (nmu>0) then
                  !
                  iz = uv_index_z_nm(nmu)
                  !
                  if (indwm(iz)==0) then
                     !
                     ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                     !
                     iok = .true.
                     !
                  endif   
               endif
            endif
            !
            if (.not. iok) then
               !
               nmu = z_index_uv_nd1(ip)    ! index of 1st uv neighbor below
               !
               if (nmu>0) then
                  !
                  iz = uv_index_z_nm(nmu)
                  !
                  if (indwm(iz)==0) then
                     !
                     ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                     !
                     iok = .true.
                     !
                  endif   
               endif
            endif
            !
            if (.not. iok) then
               !
               nmu = z_index_uv_nd2(ip)    ! index of 2nd uv neighbor below
               !
               if (nmu>0) then
                  !
                  iz = uv_index_z_nm(nmu)
                  !
                  if (indwm(iz)==0) then
                     !
                     ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                     !
                     iok = .true.
                     !
                  endif   
               endif
            endif
            !
         else
            !
            ! Check right and below
            !
            nmu = z_index_uv_mu1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (indwm(iz)==0) then
                  !
                  ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                  !
                  iok = .true.
                  !
               endif   
            endif
            !
            if (.not. iok) then
               !
               nmu = z_index_uv_mu2(ip)    ! index of 2nd uv neighbor to the right
               !
               if (nmu>0) then
                  !
                  iz = uv_index_z_nmu(nmu)
                  !
                  if (indwm(iz)==0) then
                     !
                     ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                     !
                     iok = .true.
                     !
                  endif   
               endif
            endif
            !
            if (.not. iok) then
               !
               nmu = z_index_uv_nd1(ip)    ! index of 1st uv neighbor below
               !
               if (nmu>0) then
                  !
                  iz = uv_index_z_nm(nmu)
                  !
                  if (indwm(iz)==0) then
                     !
                     ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                     !
                     iok = .true.
                     !
                  endif   
               endif
            endif
            !
            if (.not. iok) then
               !
               nmu = z_index_uv_nd2(ip)    ! index of 2nd uv neighbor below
               !
               if (nmu>0) then
                  !
                  iz = uv_index_z_nm(nmu)
                  !
                  if (indwm(iz)==0) then
                     !
                     ! This neighbor is not a wavemaker, so this point is a valid wave maker point
                     !
                     iok = .true.
                     !
                  endif   
               endif
            endif
         endif   
         !
         if (iok) then
            !
            ! This is a valid wave maker point
            !
            kcs(ip) = 4
            !
            nok = nok + 1
            !
         endif   
         !
      endif      
   enddo   
   !
   allocate(z_index_wavemaker(np))
   z_index_wavemaker = 0
   !
   allocate(wavemaker_nmd(nok))
   allocate(wavemaker_nmu(nok))
   allocate(wavemaker_ndm(nok))
   allocate(wavemaker_num(nok))
   wavemaker_nmd = 0
   wavemaker_nmu = 0
   wavemaker_ndm = 0
   wavemaker_num = 0
   !
   ! Now set the uv wave maker points (first count them)
   !
   iwm = 0
   nok = 0
   !
   do ip = 1, np
      !
      if (kcs(ip)==4) then
         !
         nok = nok + 1
         !
         z_index_wavemaker(ip) = nok
         !
         if (phi(ip)>=0.0 .and. phi(ip)<0.5*pi) then
            !
            ! Check right
            !
            nmu = z_index_uv_mu1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_mu2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
            ! Check above
            !
            nmu = z_index_uv_nu1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_nu2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
         elseif (phi(ip)>=0.5*pi .and. phi(ip)<1.0*pi) then
            !
            ! Check left
            !
            nmu = z_index_uv_md1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_md2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
            ! Check above
            !
            nmu = z_index_uv_nu1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_nu2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
         elseif (phi(ip)>=1.0*pi .and. phi(ip)<1.5*pi) then
            !
            ! Check left
            !
            nmu = z_index_uv_md1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_md2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
            ! Check below
            !
            nmu = z_index_uv_nd1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_nd2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
         else
            !
            ! Check right
            !
            nmu = z_index_uv_mu1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_mu2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
            ! Check below
            !
            nmu = z_index_uv_nd1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_nd2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !
               endif
               !
            endif   
            !
         endif
      endif
   enddo   
   !
   ! Allocate arrays
   !
   write(*,*)'Number of wavemaker u/v points : ', iwm
   !
   wavemaker_nr_uv_points = iwm   
   !
   allocate(wavemaker_index_uv(iwm))
   allocate(wavemaker_index_nmi(iwm))
   allocate(wavemaker_index_nmb(iwm))
   allocate(wavemaker_idir(iwm))
   allocate(wavemaker_angfac(iwm))
   allocate(wavemaker_uvmean(iwm))
   allocate(wavemaker_uvtrend(iwm))
   !
   wavemaker_uvmean  = 0.0
   wavemaker_uvtrend = 0.0
   !
   ! Okay, we counted the number of uv wavemaker points
   ! Now let's set them (same procedure)
   !
   iwm = 0
   nok = 0
   !
   do ip = 1, np
      !
      if (kcs(ip)==4) then
         !
         nok = nok + 1
         !
         if (phi(ip)>=0.0 .and. phi(ip)<0.5*pi) then
            !
            ! Check right
            !
            nmu = z_index_uv_mu1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = 1
                  wavemaker_angfac(iwm)    = max(cos(phi(ip) - 0.0), 0.0)
                  !
                  wavemaker_nmu(nok) = iwm
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_mu2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = 1
                  wavemaker_angfac(iwm)    = max(cos(phi(ip) - 0.0), 0.0)
                  !
               endif
               !
            endif   
            !
            ! Check above
            !
            nmu = z_index_uv_nu1(ip)    ! index of 1st uv neighbor above
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = 1
                  wavemaker_angfac(iwm)    = max(sin(phi(ip) - 0.0), 0.0)
                  !
                  wavemaker_num(nok) = iwm
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_nu2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = 1
                  wavemaker_angfac(iwm)    = max(sin(phi(ip) - 0.0), 0.0)
                  !
               endif
               !
            endif   
            !
         elseif (phi(ip)>=0.5*pi .and. phi(ip)<1.0*pi) then
            !
            ! Check left
            !
            nmu = z_index_uv_md1(ip)    ! index of 1st uv neighbor to the left
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = -1
                  wavemaker_angfac(iwm)    = max(cos(pi - phi(ip)), 0.0)
                  !
                  wavemaker_nmd(nok) = iwm
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_md2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = -1
                  wavemaker_angfac(iwm)    = max(cos(pi - phi(ip)), 0.0)
                  !
               endif
               !
            endif   
            !
            ! Check above
            !
            nmu = z_index_uv_nu1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = 1
                  wavemaker_angfac(iwm)    = max(sin(phi(ip) - 0.0), 0.0)
                  !
                  wavemaker_num(nok) = iwm
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_nu2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = 1
                  wavemaker_angfac(iwm)    = max(sin(phi(ip) - 0.0), 0.0)
                  !
               endif
               !
            endif   
            !
         elseif (phi(ip)>=1.0*pi .and. phi(ip)<1.5*pi) then
            !
            ! Check left
            !
            nmu = z_index_uv_md1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = -1
                  wavemaker_angfac(iwm)    = max(cos(pi - phi(ip)), 0.0)
                  !
                  wavemaker_nmd(nok) = iwm
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_md2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = -1
                  wavemaker_angfac(iwm)    = max(cos(pi - phi(ip)), 0.0)
                  !
               endif
               !
            endif   
            !
            ! Check below
            !
            nmu = z_index_uv_nd1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = -1
                  wavemaker_angfac(iwm)    = max(sin(pi - phi(ip)), 0.0)
                  !
                  wavemaker_ndm(nok) = iwm
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_nd2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = -1
                  wavemaker_angfac(iwm)    = max(sin(pi - phi(ip)), 0.0)
                  !
               endif
               !
            endif   
         else
            !
            ! Check right
            !
            nmu = z_index_uv_mu1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = 1
                  wavemaker_angfac(iwm)    = max(cos(phi(ip) - 0.0), 0.0)
                  !
                  wavemaker_nmu(nok) = iwm
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_mu2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = 1
                  wavemaker_angfac(iwm)    = max(cos(phi(ip) - 0.0), 0.0)
                  !
               endif
               !
            endif   
            !
            ! Check below
            !
            nmu = z_index_uv_nd1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = -1
                  wavemaker_angfac(iwm)    = max(sin(pi - phi(ip)), 0.0)
                  !
                  wavemaker_ndm(nok) = iwm
                  !
               endif
               !
            endif   
            !
            nmu = z_index_uv_nd2(ip)    ! index of 2nd uv neighbor to the right
            !
            if (nmu>0) then
               !
               iz = uv_index_z_nm(nmu)
               !
               if (kcs(iz) == 1) then
                  !
                  iwm = iwm + 1
                  !               
                  wavemaker_index_uv(iwm)  = nmu
                  wavemaker_index_nmi(iwm) = iz
                  wavemaker_index_nmb(iwm) = ip
                  wavemaker_idir(iwm)      = -1
                  wavemaker_angfac(iwm)    = max(sin(pi - phi(ip)), 0.0)
                  !
               endif
               !
            endif   
            !
         endif
      endif
   enddo
   !
   ! Set flags for kcuv points
   !
   do iwm = 1, wavemaker_nr_uv_points
      !
      ip = wavemaker_index_uv(iwm)
      kcuv(ip) = 4
      !
   enddo   
   !
   ! Turn off advection near wave maker points
   !
   if (advection .or. coriolis) then
      !
      do ip = 1, npuv
         !
         if (kcuv(uv_index_u_nmd(ip)) == 4)  uv_flags_adv(ip) = 0
         if (kcuv(uv_index_u_nmu(ip)) == 4)  uv_flags_adv(ip) = 0
         if (kcuv(uv_index_u_num(ip)) == 4)  uv_flags_adv(ip) = 0
         if (kcuv(uv_index_u_ndm(ip)) == 4)  uv_flags_adv(ip) = 0
         if (kcuv(uv_index_v_ndm(ip)) == 4)  uv_flags_adv(ip) = 0
         if (kcuv(uv_index_v_nm(ip)) == 4)   uv_flags_adv(ip) = 0
         if (kcuv(uv_index_v_nmu(ip)) == 4)  uv_flags_adv(ip) = 0
         if (kcuv(uv_index_v_ndmu(ip)) == 4) uv_flags_adv(ip) = 0
         !
      enddo      
   endif
   !
   ! In case of forcing by boundary condition files, determine indices in bwv file
   !
   ! Read wave maker forcing points
   !
   nwmfp  = 0     ! Total number of wave maker forcing points
   ntwmfp = 0     ! Number of time steps in wave maker forcing time series
   itwmfplast = 1 ! Last time point read in time series file 
   !
   wavemaker_timeseries = .false.
   !
   if (wfpfile(1:4) /= 'none') then
      !
      wavemaker_timeseries = .true.
      !
      write(*,*)'Reading wave conditions at wave makers ...'
      !
      ! Locations
      !
      open(500, file=trim(wfpfile))
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         nwmfp = nwmfp + 1
      enddo
      rewind(500)
      allocate(x_wmfp(nwmfp))
      allocate(y_wmfp(nwmfp))
      do n = 1, nwbnd
         read(500,*)x_wmfp(n),y_wmfp(n)
      enddo
      close(500)
      !
      ! Wave time series
      !
      ! First find times in whi file
      !
      open(500, file=trim(whifile))
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         ntwmfp = ntwmfp + 1
      enddo
      close(500)
      !
      allocate(wmf_time(ntwmfp))
      !
      ! Instantaneous values at wave maker forcing points
      allocate(wmf_hm0_ig_t(nwmfp))
      allocate(wmf_tp_ig_t(nwmfp))
      allocate(wmf_setup_t(nwmfp))
      !
      ! Hm0 IG (significant wave height)
      ! Times in wti and wst files must be the same as in whi file!
      !
      open(500, file=trim(whifile))
      allocate(wmf_hm0_ig(nwmfp, ntwmfp))
      do itb = 1, ntwmfp
         read(500,*)wmf_time(itb),(wmf_hm0_ig(ib, itb), ib = 1, nwmfp)
      enddo
      close(500)
      !
      ! Tp IG (peak period)
      !
      open(500, file=trim(wtifile))
      allocate(wmf_tp_ig(nwmfp, ntwmfp))
      do itb = 1, ntwmfp
         read(500,*)wmf_time(itb),(wmf_tp_ig(ib, itb), ib = 1, nwmfp)
      enddo
      close(500)
      !
      ! Set-up
      !
      allocate(wmf_setup(nwmfp, ntwmfp))
      wmf_setup = 0.0
      if (wstfile(1:4) /= 'none') then
         open(500, file=trim(wstfile))
         do itb = 1, ntwmfp
            read(500,*)wmf_time(itb),(wmf_setup(ib, itb), ib = 1, nwmfp)
         enddo
         close(500)
      endif
      !
      ! Now determine weights and indices of wave maker forcing points for each uv point  
      !
      allocate(wavemaker_index_wmfp1(wavemaker_nr_uv_points))
      allocate(wavemaker_index_wmfp2(wavemaker_nr_uv_points))
      allocate(wavemaker_fac_wmfp(wavemaker_nr_uv_points))
      !   
      do iwm = 1, wavemaker_nr_uv_points
         !
         nmb    = wavemaker_index_nmb(iwm)
         !
         x = z_xz(nmb) ! x-coordinate of cell centre behind wave maker u/v point
         y = z_yz(nmb) ! x-coordinate of cell centre behind wave maker u/v point
         !
         if (nwmfp>1) then ! More than one wave maker forcing point
            !
            dst1 = 1.0e10
            dst2 = 1.0e10
            ib1 = 0
            ib2 = 0
            !
            ! Loop through all water level boundary points
            !
            do ic = 1, nwmfp
               !
               ! Compute distance of this point to grid boundary point
               !
               dst = sqrt((x_wmfp(ic) - x)**2 + (y_wmfp(ic) - y)**2)
               !
               if (dst<dst1) then
                  !
                  ! Nearest point found
                  !
                  dst2 = dst1
                  ib2  = ib1
                  dst1 = dst
                  ib1  = ic
                  !
               elseif (dst<dst2) then
                  !
                  ! Second nearest point found
                  !
                  dst2 = dst
                  ib2  = ic
                  !
               endif
            enddo
            !
            wavemaker_index_wmfp1(iwm) = ib1
            wavemaker_index_wmfp2(iwm) = ib2
            wavemaker_fac_wmfp(iwm)    = dst2/(dst1 + dst2)
            !
         else
            !
            wavemaker_index_wmfp1(iwm) = 1
            wavemaker_index_wmfp2(iwm) = 1
            wavemaker_fac_wmfp(iwm)    = 1.0
            !
         endif
         !
      enddo
      !
   endif   
   !
   ! Infragravity frequencies
   !
   allocate(freqig(nfreqsig))
   allocate(costig(nfreqsig))
   allocate(phiig(nfreqsig))
   allocate(dphiig(nfreqsig))
   dfreqig = (freqmaxig - freqminig)/nfreqsig
   do ifreq = 1, nfreqsig
      freqig(ifreq) = freqminig + ifreq*dfreqig - 0.5*dfreqig
      call RANDOM_NUMBER(r)
      phiig(ifreq) = r*2*3.1416
!      call RANDOM_NUMBER(r)
      dphiig(ifreq) = 1.0e-6*2*3.1416/freqig(ifreq)
   enddo
   !
   end subroutine


   subroutine update_wavemaker_fluxes(t, dt, tloop)
   !
   ! Update fluxes qx and qy at wave maker points
   !
   use sfincs_data
   use sfincs_snapwave
   !
   implicit none
   !
   integer ib, nm, nmi, nmb, iuv, indb, idir, ip, ifreq, itb
   real*4  hnmb, dt, zsnmi, zsnmb, zs0nmb, zwav
   real*4  alpha, beta
   real*8  t
   real*4  sinth
   real*4  a, fp
   real*4 tbfac, hs, tp, tpsum, setup
   !
   real*4 ui, ub, dzuv, facint, zsuv, depthuv, uvm0
   !
   integer  :: count0
   integer  :: count1
   integer  :: count_rate
   integer  :: count_max
   real     :: tloop
   !
   call system_clock(count0, count_rate, count_max)
   !
   ! First update wave heights at boundary points (only need if conditions come from time series file)
   !
   if (wavemaker_timeseries) then
      !
      ! Interpolate boundary conditions in time
      !
      do itb = itwmfplast, ntwmfp ! Loop in time
         !
         if (wmf_time(itb)>t) then
            !
            tbfac  = (t - wmf_time(itb - 1))/(wmf_time(itb) - wmf_time(itb - 1))
            !
            tpsum = 0.0
            !
            do ib = 1, nwmfp ! Loop along forcing points
               !
               hs    = wmf_hm0_ig(ib, itb - 1) + (wmf_hm0_ig(ib, itb) - wmf_hm0_ig(ib, itb - 1))*tbfac
               tp    = wmf_tp_ig(ib, itb - 1)  + (wmf_tp_ig(ib, itb)  - wmf_tp_ig(ib, itb - 1))*tbfac
               setup = wmf_setup(ib, itb - 1)  + (wmf_setup(ib, itb)  - wmf_setup(ib, itb - 1))*tbfac
               !
               wmf_hm0_ig_t(ib) = hs
               wmf_setup_t(ib)  = setup
               !
               tpsum = tpsum + tp
               !
            enddo
            !
            tp = tpsum / nwmfp ! Take average Tp from boundary points
            !
            itwmfplast = itb - 1
            exit
            !
         endif
      enddo
   else
      !
      ! Use mean peak period from SnapWave boundary conditions
      !
      tp = snapwave_tpmean*6
      !
   endif      
   !
   alpha   = min(dt/wmtfilter, 1.0)
   beta    = min(dt/(0.2*wmtfilter), 1.0)
   !
   if (wavemaker_spectrum) then
      !
      ! First update phases
      !
      do ifreq = 1, nfreqsig
         phiig(ifreq) = phiig(ifreq) + dphiig(ifreq)*dt
         costig(ifreq) = cos(2*pi*t*freqig(ifreq) + phiig(ifreq))
      enddo
      !
      fp = 1.0/tp ! Wave period
      !
      ! Now spectrum and wave excitation
      !
      zwav = 0.0
      !
      do ifreq = 1, nfreqsig
         a = 0.125*(1.0**2)*(fp**(-2))*freqig(ifreq)*exp(-freqig(ifreq)/fp)
         zwav = zwav + costig(ifreq)*sqrt(a*dfreqig)
      enddo
      !
   else
      !
      ! Monochromatic signal
      !
      zwav = 0.5*sin(2*pi*t/tp)
      !
   endif   
   !
   if (t<tspinup) then
      !
      zwav = zwav * (t - t0)/(tspinup - t0)
      !
   endif      
   !
   !$acc kernels present( wavemaker_index_uv, wavemaker_index_nmi, wavemaker_index_nmb, &
   !$acc                  zs, q, hm0_ig, zb, zbuv, &
   !$acc                  wmf_hm0_ig_t, wmf_setup_t, wavemaker_index_wmfp1, wavemaker_index_wmfp2, wavemaker_fac_wmfp, &
   !$acc                  wavemaker_uvmean, wavemaker_idir, wavemaker_angfac, wavemaker_freduv, wavemaker_uvtrend, & 
   !$acc                  subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_hrep, subgrid_uv_navg, subgrid_z_zmin, subgrid_uv_hrep_zmax), async(1)
   ! 
   ! UV fluxes at boundaries
   !
   !$acc loop independent, private(ib)
   do ib = 1, wavemaker_nr_uv_points
      !
      ip     = wavemaker_index_uv(ib)
      nmi    = wavemaker_index_nmi(ib)
      nmb    = wavemaker_index_nmb(ib)
      !
      zsnmi  = zs(nmi)                    ! total water level on wave side inside model
      !
      if (wavemaker_timeseries) then
         !
         ! Take wave height from boundary conditions file (weighted average of two nearby forcing points)
         !
         hs    = wmf_hm0_ig_t(wavemaker_index_wmfp1(ib))*wavemaker_fac_wmfp(ib) + wmf_hm0_ig_t(wavemaker_index_wmfp2(ib))*(1.0 - wavemaker_fac_wmfp(ib))
         setup = wmf_setup_t(wavemaker_index_wmfp1(ib))*wavemaker_fac_wmfp(ib)  + wmf_setup_t(wavemaker_index_wmfp2(ib))*(1.0 - wavemaker_fac_wmfp(ib))
         !
         zs0nmb = zs(nmb) + setup            ! average water level inside model without waves (this should be zs)
         zsnmb  = zs0nmb + zwav*hs           ! total water level in wave maker (i.e. mean water level plus wave)         
         !
      else
         !
         ! Take wave height from SnapWave
         !
         zs0nmb = zs(nmb)                    ! average water level inside model without waves (this should be zs)
         zsnmb  = zs0nmb + zwav*hm0_ig(nmb)  ! total water level in wave maker (i.e. mean water level plus wave)
         !
      endif   
      !
      if (subgrid) then
         !
         zsuv = max(zsnmb, zsnmi)
         !
         if (zsuv>=subgrid_uv_zmax(nm) - 1.0e-3) then
            !
            ! Entire cell is wet, no interpolation from table needed
            !
            depthuv  = subgrid_uv_hrep_zmax(ip) + zsuv
            !
         elseif (zsuv>subgrid_uv_zmin(nm)) then
            !
            ! Interpolation required
            !            
            dzuv    = (subgrid_uv_zmax(ip) - subgrid_uv_zmin(ip)) / (subgrid_nbins - 1)
            iuv     = int((zsuv - subgrid_uv_zmin(ip))/dzuv) + 1
            facint  = (zsuv - (subgrid_uv_zmin(ip) + (iuv - 1)*dzuv) ) / dzuv
            depthuv = subgrid_uv_hrep(iuv, ip) + (subgrid_uv_hrep(iuv + 1, ip) - subgrid_uv_hrep(iuv, ip))*facint
            !
         else
            !
            depthuv = 0.0
            !
         endif
         !
         hnmb   = max(depthuv, huthresh)
         zsnmb  = max(zsnmb,  subgrid_z_zmin(nmb))
         zs0nmb = max(zs0nmb, subgrid_z_zmin(nmb))
         !
      else
         !
         hnmb   = max(0.5*(zsnmb + zsnmi) - zbuv(ip), huthresh)
         zsnmb  = max(zsnmb,  zb(nmb))
         zs0nmb = max(zs0nmb, zb(nmb))
         !
      endif
      !
      ! Weakly reflective boundary
      !
      if (hnmb<huthresh) then
         !
         ! Very shallow
         !
         q(ip) = 0.0
         wavemaker_uvmean(ib) = 0.0
         !
      else
         !
         ui = sqrt(g/hnmb)*(zsnmb - zs0nmb)
         ub = wavemaker_idir(ib) * (2*ui - sqrt(g/hnmb)*(zsnmi - zs0nmb)) * wavemaker_angfac(ib)
         !
         q(ip) = ub*hnmb + wavemaker_uvmean(ib)
         !
      endif
      !
      if (wmtfilter>=0.0) then
         !
         ! Use double exponential time filter
         !
         uvm0 = wavemaker_uvmean(ib) ! Previous time step
         wavemaker_uvmean(ib)  = alpha*q(ip)   + wavemaker_freduv*(1.0 - alpha)*(wavemaker_uvmean(ib) + wavemaker_uvtrend(ib))
         wavemaker_uvtrend(ib) = beta*(wavemaker_uvmean(ib) - uvm0) + (1.0 - beta)*wavemaker_uvtrend(ib)
         !
      else
         !
         wavemaker_uvmean(ib) = 0.0
         !
      endif
      !
   enddo
   !
   !$acc end kernels
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine
      
   end module
