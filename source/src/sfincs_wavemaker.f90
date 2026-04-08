   module sfincs_wavemaker

   use sfincs_log
   use sfincs_error

   contains

   subroutine initialize_wavemakers()
   !
   ! Reads polylines files
   ! Computes cross sections
   ! Determines interpolation weights
   ! Interpolates bathymetry
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   integer*4 :: nrows
   integer*4 :: ncols
   integer*4 :: stat
   integer*4 :: ipol
   integer*4 :: irow   
   integer*4 :: j
   integer*4 :: n
   integer*4 :: nmu
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
   !
   character :: cdummy
   !
   real*4,    dimension(:),     allocatable :: xpol
   real*4,    dimension(:),     allocatable :: ypol
   real*4,    dimension(:),     allocatable :: phi
   !
   integer*4, dimension(:),     allocatable :: cell_indices
   integer*4, dimension(:),     allocatable :: indwm
   !
   real*4, dimension(:),     allocatable :: wavemaker_xfp
   real*4, dimension(:),     allocatable :: wavemaker_yfp   
   !
   logical :: iok, ok
   !
   integer ib1, ib2, ib, ic, nmb, nrwvm
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
   nrwvm = 0
   !
   write(logstr,*)'Reading wavemaker polyline file ...'
   call write_log(logstr, 0)
   !
   ok = check_file_exists(wavemaker_wvmfile, 'Wave maker wvm file', .true.)
   !
   open(500, file=trim(wavemaker_wvmfile))
   do while(.true.)
      read(500,*,iostat = stat)cdummy
      if (stat<0) exit      
      read(500,*,iostat = stat)nrows,ncols
      if (stat<0) exit
      nrwvm = nrwvm + 1
      do irow = 1, nrows
         read(500,*)dummy
      enddo
   enddo
   rewind(500)
   !
   ! Loop through polylines
   !
   write(logstr,*)'Number of wavemaker polylines found : ', nrwvm
   call write_log(logstr, 0)   
   !
   do ipol = 1, nrwvm
      !
      read(500,*,iostat = stat)cdummy
      if (stat<0) exit
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
         phip = atan2(ypol(irow + 1) - ypol(irow), xpol(irow + 1) - xpol(irow)) + 0.5 * pi
         phip = phip - rotation
         if (phip >= 2 * pi) phip = phip - 2 * pi
         if (phip < 0.0)   phip = phip + 2 * pi
         !
         call find_cells_intersected_by_line(cell_indices, nr_cells, xpol(irow), ypol(irow), xpol(irow + 1), ypol(irow + 1))
         !
         do j = 1, nr_cells
            !
            ip = index_sfincs_in_quadtree(cell_indices(j))
            !
            if (ip > 0) then
               !
               if (indwm(ip) == 0) then
                  !
                  indwm(ip) = 1 ! set temporary flag to 1
                  phi(ip)   = phip
                  nrw       = nrw + 1
                  !
               endif   
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
      if (indwm(ip) == 1) then
         !
         iok = .false.
         !
         if (phi(ip) >= 0.0 .and. phi(ip) < 0.5 * pi) then
            !
            ! Check right and above
            !
            nmu = z_index_uv_mu1(ip)    ! index of 1st uv neighbor to the right
            !
            if (nmu > 0) then
               !
               iz = uv_index_z_nmu(nmu)
               !
               if (indwm(iz) == 0) then
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
               if (nmu > 0) then
                  !
                  iz = uv_index_z_nmu(nmu)
                  !
                  if (indwm(iz) == 0) then
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
         if (iok .and. kcs(ip) == 1) then
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
         if (phi(ip) >= 0.0 .and. phi(ip) < 0.5*pi) then
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
   write(logstr,*)'Number of wavemaker u/v points : ', iwm
   call write_log(logstr, 0)
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
   write(logstr,*)'Setting wave makers ...'
   call write_log(logstr, 0)   
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
   call write_log(logstr, 0)   
   !
   do iwm = 1, wavemaker_nr_uv_points
      !
      ip = wavemaker_index_uv(iwm)
      kcuv(ip) = 4
      !
   enddo   
   !
   ! In case of forcing by boundary condition files, determine indices in bwv file
   !
   ! Read wave maker forcing points
   !
   wavemaker_nr_forcing_points  = 0     ! Total number of wave maker forcing points
   wavemaker_nr_forcing_timesteps = 0     ! Number of time steps in wave maker forcing time series
   wavemaker_itlast = 1 ! Last time point read in time series file 
   !
   wavemaker_timeseries = .false.
   !
   if (wavemaker_wfpfile(1:4) /= 'none') then
      !
      wavemaker_timeseries = .true.
      !
      write(logstr,*)'Reading wave conditions at wave makers ...'
      call write_log(logstr, 0)      
      !
      ! Locations
      !
      ok = check_file_exists(wavemaker_wfpfile, 'Wave maker wfp file', .true.)
      !
      open(500, file=trim(wavemaker_wfpfile))
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         wavemaker_nr_forcing_points = wavemaker_nr_forcing_points + 1
      enddo
      rewind(500)
      allocate(wavemaker_xfp(wavemaker_nr_forcing_points))
      allocate(wavemaker_yfp(wavemaker_nr_forcing_points))
      do n = 1, wavemaker_nr_forcing_points
         read(500,*)wavemaker_xfp(n),wavemaker_yfp(n)
      enddo
      close(500)
      !
      ! Wave time series
      !
      ! First find times in whi file
      !
      ok = check_file_exists(wavemaker_wfpfile, 'Wave maker wfp file', .true.)
      !      
      open(500, file=trim(wavemaker_whifile))
      do while(.true.)
         read(500,*,iostat = stat)dummy
         if (stat<0) exit
         wavemaker_nr_forcing_timesteps = wavemaker_nr_forcing_timesteps + 1
      enddo
      close(500)
      !
      allocate(wavemaker_forcing_time(wavemaker_nr_forcing_timesteps))
      !
      ! Hm0 IG (significant wave height)
      ! Times in wti and wst files must be the same as in whi file!
      !
      open(500, file=trim(wavemaker_whifile))
      allocate(wavemaker_forcing_hm0_ig(wavemaker_nr_forcing_points, wavemaker_nr_forcing_timesteps))
      do itb = 1, wavemaker_nr_forcing_timesteps
         read(500,*)wavemaker_forcing_time(itb),(wavemaker_forcing_hm0_ig(ib, itb), ib = 1, wavemaker_nr_forcing_points)
      enddo
      close(500)
      !
      ! Tp IG (peak period)
      !
      ok = check_file_exists(wavemaker_wtifile, 'Wave maker wti file', .true.)
      !      
      open(500, file=trim(wavemaker_wtifile))
      allocate(wavemaker_forcing_tp_ig(wavemaker_nr_forcing_points, wavemaker_nr_forcing_timesteps))
      do itb = 1, wavemaker_nr_forcing_timesteps
         read(500,*)wavemaker_forcing_time(itb),(wavemaker_forcing_tp_ig(ib, itb), ib = 1, wavemaker_nr_forcing_points)
      enddo
      close(500)
      !
      ! Set-up
      !
      allocate(wavemaker_forcing_setup(wavemaker_nr_forcing_points, wavemaker_nr_forcing_timesteps))
      wavemaker_forcing_setup = 0.0
      if (wavemaker_wstfile(1:4) /= 'none') then
         !
         ok = check_file_exists(wavemaker_wstfile, 'Wave maker wst file', .true.)
         ! 
         open(500, file=trim(wavemaker_wstfile))
         do itb = 1, wavemaker_nr_forcing_timesteps
            read(500,*)wavemaker_forcing_time(itb),(wavemaker_forcing_setup(ib, itb), ib = 1, wavemaker_nr_forcing_points)
         enddo
         close(500)
      endif
      !
      if ((wavemaker_forcing_time(1) > (t0 + 1.0)) .or. (wavemaker_forcing_time(wavemaker_nr_forcing_timesteps) < (t1 - 1.0))) then
         ! 
         write(logstr,'(a)')' WARNING! Times in wave maker time series file do not cover entire simulation period !'
         call write_log(logstr, 0)         
         ! 
         if (wavemaker_forcing_time(1) > (t0 + 1.0)) then
            ! 
            write(logstr,'(a)')' WARNING! Adjusting first time in wave maker time series !'
            call write_log(logstr, 0)                                 
            !
            wavemaker_forcing_time(1) = t0 - 1.0
            !
         else
            ! 
            write(logstr,'(a)')' WARNING! Adjusting last time in wave maker time series !'
            call write_log(logstr, 0)                     
            !
            wavemaker_forcing_time(wavemaker_nr_forcing_timesteps) = t1 + 1.0
            !
         endif
         !
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
         if (wavemaker_nr_forcing_points>1) then ! More than one wave maker forcing point
            !
            dst1 = 1.0e10
            dst2 = 1.0e10
            ib1 = 0
            ib2 = 0
            !
            ! Loop through all water level boundary points
            !
            do ic = 1, wavemaker_nr_forcing_points
               !
               ! Compute distance of this point to grid boundary point
               !
               dst = sqrt((wavemaker_xfp(ic) - x)**2 + ( wavemaker_yfp(ic) - y)**2)
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
   allocate(wavemaker_freq_ig(wavemaker_nfreqs_ig))
   allocate(wavemaker_cost_ig(wavemaker_nfreqs_ig))
   allocate(wavemaker_phi_ig(wavemaker_nfreqs_ig))
   allocate(wavemaker_dphi_ig(wavemaker_nfreqs_ig))
   wavemaker_dfreq_ig = (wavemaker_freqmax_ig - wavemaker_freqmin_ig) / wavemaker_nfreqs_ig
   do ifreq = 1, wavemaker_nfreqs_ig
      wavemaker_freq_ig(ifreq) = wavemaker_freqmin_ig + ifreq * wavemaker_dfreq_ig - 0.5 * wavemaker_dfreq_ig
      call random_number(r)
      wavemaker_phi_ig(ifreq) = r * 2 * 3.1416
      wavemaker_dphi_ig(ifreq) = 1.0e-6 * 2 * 3.1416 / wavemaker_freq_ig(ifreq)
   enddo
   !
   if (wavemaker_hinc) then
      !
      allocate(wavemaker_freq_inc(wavemaker_nfreqs_inc))
      allocate(wavemaker_cost_inc(wavemaker_nfreqs_inc))
      allocate(wavemaker_phi_inc(wavemaker_nfreqs_inc))
      allocate(wavemaker_dphi_inc(wavemaker_nfreqs_inc))
      wavemaker_dfreq_inc = (wavemaker_freqmax_inc - wavemaker_freqmin_inc) / wavemaker_nfreqs_inc
      do ifreq = 1, wavemaker_nfreqs_inc
         wavemaker_freq_inc(ifreq) = wavemaker_freqmin_inc + ifreq * wavemaker_dfreq_inc - 0.5 * wavemaker_dfreq_inc
         call random_number(r)
         wavemaker_phi_inc(ifreq) = r * 2 * 3.1416
         wavemaker_dphi_inc(ifreq) = 1.0e-6 * 2 * 3.1416 / wavemaker_freq_inc(ifreq)
      enddo
      !
   endif   
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
   integer :: ib, nmi, nmb, iuv, ip, ifreq, itb, itb0, itb1, kst
   real*4  :: hnmb, dt, zsnmi, zsnmb, zs0nmb, zwav_ig, zwav_inc
   real*4  :: alpha, beta
   real*8  :: t, tb
   real*4  :: tbfac, hs, tp_ig, tp_inc, tpsum, setup, fm_ig, a, fm_inc
   real*4  :: wave_steepness, betas, zinc, zig, dwvm, ztot, hm0_inc
   real*4  :: ui, ub, dzuv, facint, zsuv, depthuv, uvm0
   !
   integer  :: count0
   integer  :: count1
   integer  :: count_rate
   integer  :: count_max
   real     :: tloop
   !
   real*4, dimension(:),     allocatable :: wavemaker_forcing_hm0_ig_t
   real*4, dimension(:),     allocatable :: wavemaker_forcing_tp_ig_t
   real*4, dimension(:),     allocatable :: wavemaker_forcing_setup_t   
   !
   call system_clock(count0, count_rate, count_max)
   !
   ! Factors for double-exponential filtering
   !
   alpha = min(dt / wavemaker_filter_time, 1.0)
   beta  = min(dt / (0.2 * wavemaker_filter_time), 1.0)
   !
   ! For time series forcing, we now update values at the forcing points and determine Tp_ig
   ! For forcing with SnapWave, we only need to determine Tp_ig
   !
   if (wavemaker_timeseries) then
      !
      ! Only IG wave forcing supported at the moment !
      !
      allocate(wavemaker_forcing_hm0_ig_t(wavemaker_nr_forcing_points))
      allocate(wavemaker_forcing_tp_ig_t(wavemaker_nr_forcing_points))
      allocate(wavemaker_forcing_setup_t(wavemaker_nr_forcing_points))
      !
      ! Interpolate boundary conditions in time
      !
      if (wavemaker_forcing_time(1) > t - 1.0e-3) then ! use first time in boundary conditions
         !
         itb0 = 1
         itb1 = 1
         tb   = wavemaker_forcing_time(itb0)
         !
      elseif (wavemaker_forcing_time(wavemaker_nr_forcing_timesteps) < t + 1.0e-3) then  ! use last time in boundary conditions       
         !
         itb0 = wavemaker_nr_forcing_timesteps
         itb1 = wavemaker_nr_forcing_timesteps
         tb   = wavemaker_forcing_time(itb0)
         !
      else
         !
         do itb = wavemaker_itlast, wavemaker_nr_forcing_timesteps ! Loop in time
            if (wavemaker_forcing_time(itb) > t + 1.0e-6) then
               itb0 = itb - 1
               itb1 = itb
               tb   = t
               wavemaker_itlast = itb - 1
               exit
            endif
         enddo 
         !
      endif            
      !
      tbfac  = (tb - wavemaker_forcing_time(itb0)) / max(wavemaker_forcing_time(itb1) - wavemaker_forcing_time(itb0), 1.0e-6)
      !
      tpsum = 0.0
      !
      do ib = 1, wavemaker_nr_forcing_points ! Loop along forcing points
         !
         hs    = wavemaker_forcing_hm0_ig(ib, itb0) + (wavemaker_forcing_hm0_ig(ib, itb1) - wavemaker_forcing_hm0_ig(ib, itb0)) * tbfac
         tp_ig = wavemaker_forcing_tp_ig(ib, itb0)  + (wavemaker_forcing_tp_ig(ib, itb1)  - wavemaker_forcing_tp_ig(ib, itb0)) * tbfac
         setup = wavemaker_forcing_setup(ib, itb0)  + (wavemaker_forcing_setup(ib, itb1)  - wavemaker_forcing_setup(ib, itb0)) * tbfac
         !
         wavemaker_forcing_hm0_ig_t(ib) = hs
         wavemaker_forcing_setup_t(ib)  = setup
         !
         tpsum = tpsum + tp_ig
         !
      enddo
      !
      tp_ig = tpsum / wavemaker_nr_forcing_points ! Take average Tp from boundary points
      tp_inc = 10.0 ! update this!
      !
   else
      !
      ! Use mean peak period from SnapWave boundary conditions
      !
      tp_ig = snapwave_tpigmean ! TL: Now calculated in SnapWave, different options for using a period based on Herbers spectrum (snapwave_tpig_opt, if snapwave_use_herbers=1, or user defined snapwave_Tinc2ig ratio (if snapwave_use_herbers = 0)
      !
      ! We may want to use Herbers for computation of IG waves in SnapWave, but we want to have control over peak IG period at wave makers.
      !
      if (wavemaker_Tinc2ig > 0.0) then
         !
         ! Use factor on mean Tp_inc at boundaries
         !
         tp_ig = snapwave_tpmean * wavemaker_Tinc2ig
         !
!      elseif (wavemaker_surfslope > 0.0) then ! Dean a
!         !
!         ! Turn this option off now, because snapwave_hsmean is not available in current branch
!         ! Will need to be updated if we want to use this option, but it is not a priority at the moment
!         !
!         ! Estimate surfzone slope from Dean's a, using gambr = 1.0
!         !
!         betas = snapwave_hsmean / (snapwave_hsmean / (1.0 * wavemaker_surfslope))**(3.0 / 2.0)
!         !
!         wave_steepness = snapwave_hsmean / (1.56 * snapwave_tpmean**2)
!         !
!         ! From empirical run-up equation (van Ormondt et al., 2021), but slightly adjusted
!         !
!         tp_ig = snapwave_tpmean * max(1.86 * betas**-0.43 * wave_steepness**0.07, 5.0)
!         !
      endif
      !
      tp_inc = max(snapwave_tpmean, wavemaker_tpmin)
      ! 
   endif      
   !
   ! Now determine zwav_ig and zwav_inc based on spectrum or monochromatic signal.
   ! Time series of zwav_ig and zwav_inc will be used to modulate water level at wave maker points.
   ! They both give at Hm0 of 1.0 m, and therefore need to be scaled with the data at the wave maker points (either from time series or SnapWave boundary conditions)
   !
   zwav_ig = 0.0
   zwav_inc = 0.0
   !
   if (wavemaker_spectrum) then
      !
      ! Infragravity waves
      !
      if (wavemaker_hig) then
         !
         fm_ig = 1.0 / tp_ig ! Wave period
         !
         ! Now spectrum and wave excitation
         !
         do ifreq = 1, wavemaker_nfreqs_ig
            !
            ! Update phase
            !
            wavemaker_phi_ig(ifreq) = modulo(wavemaker_phi_ig(ifreq) + wavemaker_dphi_ig(ifreq) * dt, 2 * pi)
            wavemaker_cost_ig(ifreq) = cos(2 * pi * t * wavemaker_freq_ig(ifreq) + wavemaker_phi_ig(ifreq))         
            !
            ! Use this spectral shape instead
            !
            a = 0.125 * (fm_ig**-2) * wavemaker_freq_ig(ifreq) * (exp(-wavemaker_freq_ig(ifreq) / fm_ig))
            !
            zwav_ig = zwav_ig + wavemaker_cost_ig(ifreq) * sqrt(a * wavemaker_dfreq_ig)
            !
         enddo
         !
      endif
      !
      if (wavemaker_hinc) then
         !
         fm_inc = 1.0 / tp_inc ! Wave period
         !
         do ifreq = 1, wavemaker_nfreqs_ig
            !
            wavemaker_phi_inc(ifreq) = modulo(wavemaker_phi_inc(ifreq) + wavemaker_dphi_inc(ifreq) * dt, 2 * pi)
            wavemaker_cost_inc(ifreq) = cos(2 * pi * t * wavemaker_freq_inc(ifreq) + wavemaker_phi_inc(ifreq))
            !
            ! The ISSC spectrum (also known as Bretschneider or modified Pierson-Moskowitz)
            !
            a = 0.625 * (fm_inc**4) * (wavemaker_freq_inc(ifreq)**-5) * (exp(-1.25 * (wavemaker_freq_inc(ifreq) / fm_inc)**-4))
            !
            zwav_inc = zwav_inc + wavemaker_cost_inc(ifreq) * sqrt(a * wavemaker_dfreq_inc)  
            !
         enddo
         !
         !zwav_inc = 0.5 * sin(2 * pi * t / tp_inc)
         !
         ! Saw tooth
         !
         !zwav_inc = - mod(t, tp_inc) / tp_inc + 0.5
         !
         ! Let zwav_inc be modulated by zwav_ig (i.e. higher incident waves at the peaks of the IG wave)
         !
         !zwav_inc = zwav_inc * sqrt(max(zwav + 1.0, 0.0)) ! this assumes zwav is somewhere between -0.5 and +0.5
         !
      endif   
      !
   else
      !
      ! Monochromatic signal
      !
      if (wavemaker_hinc) then
         !
         zwav_ig = 0.5 * sin(2 * pi * t / tp_ig)
         !
      endif
      !
      if (wavemaker_hinc) then
         !
         zwav_inc = 0.5 * sin(2 * pi * t / tp_inc)
         !
      endif   
      !
   endif   
   !
   if (t < tspinup) then
      !
      zwav_ig = zwav_ig * (t - t0) / (tspinup - t0)
      zwav_inc = zwav_inc * (t - t0) / (tspinup - t0)
      !
   endif
   !
   ! UV fluxes at wave makers
   !
   ! No OMP acceleration here?
   !
   !$acc parallel present( wavemaker_index_uv, wavemaker_index_nmi, wavemaker_index_nmb, &
   !$acc                  zs, q, hm0, hm0_ig, zb, zbuv, subgrid_z_zmax, &
   !$acc                  wavemaker_forcing_hm0_ig_t, wavemaker_forcing_setup_t, wavemaker_index_wmfp1, wavemaker_index_wmfp2, wavemaker_fac_wmfp, &
   !$acc                  wavemaker_uvmean, wavemaker_idir, wavemaker_angfac, wavemaker_uvtrend, & 
   !$acc                  subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_havg, subgrid_uv_nrep, subgrid_z_zmin, subgrid_uv_havg_zmax)
   !$acc loop independent gang vector
   do ib = 1, wavemaker_nr_uv_points
      !
      ip     = wavemaker_index_uv(ib)
      nmi    = wavemaker_index_nmi(ib)
      nmb    = wavemaker_index_nmb(ib)
      !
      zsnmi  = zs(nmi)                    ! total water level on wave side inside model
      !
      ! Now determine total water levels (zs0nmb and zsnmb) on boundary (i.e. wave maker) side,
      ! which is based on mean water level plus wave height
      !
      if (wavemaker_timeseries) then
         !
         ! Take wave height from boundary conditions file (weighted average of two nearby forcing points)
         !
         hs    = wavemaker_forcing_hm0_ig_t(wavemaker_index_wmfp1(ib)) * wavemaker_fac_wmfp(ib) + wavemaker_forcing_hm0_ig_t(wavemaker_index_wmfp2(ib)) * (1.0 - wavemaker_fac_wmfp(ib))
         setup = wavemaker_forcing_setup_t(wavemaker_index_wmfp1(ib)) * wavemaker_fac_wmfp(ib)  + wavemaker_forcing_setup_t(wavemaker_index_wmfp2(ib)) * (1.0 - wavemaker_fac_wmfp(ib))
         !
         zs0nmb = zs(nmb) + setup            ! average water level inside model without waves (this should be zs)
         zsnmb  = zs0nmb + zwav_ig * hs      ! total water level in wave maker (i.e. mean water level plus wave)         
         !
      else
         !
         ! Take wave height from SnapWave
         !
         zs0nmb = zs(nmb) ! average water level inside model without waves
         !
         zig    = wavemaker_hm0_ig_factor * zwav_ig * hm0_ig(nmb)
         zinc   = wavemaker_hm0_inc_factor * zwav_inc * hm0(nmb)
         !
         ! Compute water depth including IG wave
         !
         if (subgrid) then
            dwvm   = max(zs0nmb + zig - subgrid_z_zmax(nmb), 0.0) ! depth at wave maker
         else
            dwvm   = max(zs0nmb + zig - zb(nmb), 0.0) ! depth at wave maker
         endif
         !
         ! Limit incident wave height (not IG?)
         !
         !zinc = min(zinc,  wavemaker_gammax * dwvm)
         !
         !zsnmb  = zs0nmb + zig + zinc ! total water level in wave maker (i.e. mean water level plus wave)         
         !
         zsnmb  = zs0nmb + min(zinc + zig,  wavemaker_gammax * dwvm) ! total water level in wave maker (i.e. mean water level plus wave)
         !
         !
      endif   
      !
      if (subgrid) then
         !
         zsuv = max(zsnmb, zsnmi)
         !
         if (zsuv >= subgrid_uv_zmax(ip) - 1.0e-3) then
            !
            ! Entire cell is wet, no interpolation from table needed
            !
            depthuv  = subgrid_uv_havg_zmax(ip) + zsuv
            !
         elseif (zsuv > subgrid_uv_zmin(ip)) then
            !
            ! Interpolation required
            !            
            dzuv    = (subgrid_uv_zmax(ip) - subgrid_uv_zmin(ip)) / (subgrid_nlevels - 1)
            iuv     = int((zsuv - subgrid_uv_zmin(ip)) / dzuv) + 1
            facint  = (zsuv - (subgrid_uv_zmin(ip) + (iuv - 1) * dzuv) ) / dzuv
            depthuv = subgrid_uv_havg(iuv, ip) + (subgrid_uv_havg(iuv + 1, ip) - subgrid_uv_havg(iuv, ip)) * facint
            !
         else
            !
            depthuv = 0.0
            !
         endif
         !
         hnmb   = depthuv
         zsnmb  = max(zsnmb,  subgrid_z_zmin(nmb))
         zs0nmb = max(zs0nmb, subgrid_z_zmin(nmb))
         !
      else
         !
         hnmb   = 0.5 * (zsnmb + zsnmi) - zbuv(ip)
         zsnmb  = max(zsnmb, zb(nmb))
         zs0nmb = max(zs0nmb, zb(nmb))
         !
      endif
      !
      ! Use weakly reflective boundary condition 
      !
      if (hnmb < wavemaker_hmin) then
         !
         ! Very shallow
         !
         q(ip) = 0.0
         wavemaker_uvmean(ib) = 0.0
         !
      else
         !
         ui = sqrt(g / hnmb) * (zsnmb - zs0nmb)
         ub = wavemaker_idir(ib) * (2 * ui - sqrt(g / hnmb) * (zsnmi - zs0nmb)) * wavemaker_angfac(ib)
         !
         q(ip) = ub * hnmb + wavemaker_uvmean(ib)
         !
      endif
      !
      if (wavemaker_filter_time >= 0.0) then
         !
         ! Use double exponential time filter
         !
         uvm0 = wavemaker_uvmean(ib) ! Previous time step
         wavemaker_uvmean(ib)  = alpha * q(ip) + wavemaker_filter_fred * (1.0 - alpha) * (wavemaker_uvmean(ib) + wavemaker_uvtrend(ib))
         wavemaker_uvtrend(ib) = beta * (wavemaker_uvmean(ib) - uvm0) + (1.0 - beta) * wavemaker_uvtrend(ib)
         !
      else
         !
         wavemaker_uvmean(ib) = 0.0
         !
      endif
      !
   enddo
   !$acc end parallel
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine
      
   end module
