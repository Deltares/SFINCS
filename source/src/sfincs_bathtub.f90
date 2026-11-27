module sfincs_bathtub
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! Bathtub flooding mode: water levels in all cells are set directly from the interpolated
   ! boundary water levels, without solving the momentum and continuity equations.
   ! With SnapWave, the boundary water levels are raised by a factor times the SnapWave
   ! boundary wave height; the full SnapWave solver can run alongside (bathtub_snapwave).
   !
   ! Subroutines / functions:
   !
   !   initialize_bathtub()
   !     Set up the interpolation weights per cell and switch off unused processes.
   !     Called from: sfincs_initialize (sfincs_lib)
   !
   !   initialize_bathtub_wave_setup()
   !     Add bathtub_fac_hs * Hs from the SnapWave boundary conditions to the boundary water levels.
   !     Called from: sfincs_initialize (sfincs_lib)
   !
   !   bathtub_compute_water_levels()
   !     Set the water level in every cell from the weighted boundary water levels.
   !     Called from: sfincs_update (sfincs_lib)
   !
   !   read_snapwave_boundary_data(read_files)
   !     Read (if needed) and copy the SnapWave wave boundary points for use in bathtub mode.
   !     Called from: initialize_bathtub_wave_setup (sfincs_bathtub)
   !
   !   update_snapwave_boundary_data(tb)
   !     Interpolate the SnapWave boundary wave heights to time tb.
   !     Called from: initialize_bathtub_wave_setup (sfincs_bathtub)
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

integer*4, dimension(:),   allocatable :: bathtub_i1
integer*4, dimension(:),   allocatable :: bathtub_i2
real*4,    dimension(:),   allocatable :: bathtub_w1
integer*4                              :: bathtub_snapwave_nwbnd
real*4, allocatable  :: bathtub_snapwave_x_bwv(:)
real*4, allocatable  :: bathtub_snapwave_y_bwv(:)
real*4, allocatable  :: bathtub_snapwave_hs_bwv(:)
integer, allocatable :: bathtub_snapwave_i1(:)
integer, allocatable :: bathtub_snapwave_i2(:)
real*4, allocatable  :: bathtub_snapwave_w1(:)

contains
   !
   !
   subroutine initialize_bathtub()
      !
      ! Determine the interpolation weights from the boundary points to every cell, and
      ! switch off the processes that bathtub mode does not use.
      !
      ! Called from: sfincs_initialize (sfincs_lib)
      !
      use sfincs_data
      use geometry
      !
      implicit none
      !
      integer :: ip
      integer :: i1, i2
      real*4  :: w1, w2
      !
      ! Determine indices and weight for each grid point
      !
      allocate(bathtub_i1(np))
      allocate(bathtub_i2(np))
      allocate(bathtub_w1(np))
      !
      do ip = 1, np
         !
         call interp_segment(x_bnd, y_bnd, nbnd, z_xz(ip), z_yz(ip), i1, i2, w1, w2)
         !
         bathtub_i1(ip) = i1
         bathtub_i2(ip) = i2
         bathtub_w1(ip) = w1
         !
      enddo
      !
      ! initialize_domain sets the meteo flags again from the input files, so switch them
      ! off here. Wind is only kept as input for SnapWave wind growth.
      !
      if (.not. snapwavewind) then
         !
         meteo3d    = .false.
         wind       = .false.
         store_wind = .false.
         !
      endif
      !
      store_meteo            = .false.
      store_wind_max         = .false.
      precip                 = .false.
      patmos                 = .false.
      infiltration           = .false.
      store_velocity         = .false.
      store_maximum_velocity = .false.
      !
   end subroutine initialize_bathtub
   !
   !
   subroutine initialize_bathtub_wave_setup()
      !
      ! Raise the boundary water levels by bathtub_fac_hs * Hs, with Hs interpolated from the
      ! SnapWave wave boundary conditions to every boundary point and boundary time.
      ! When the full SnapWave solver runs, its boundary conditions have already been read
      ! by couple_snapwave and are reused.
      !
      ! Called from: sfincs_initialize (sfincs_lib)
      !
      use sfincs_data
      use geometry
      use sfincs_log
      !
      implicit none
      !
      integer :: ib, itb
      integer :: i1, i2
      real*4  :: w1, w2
      real*4  :: h1, h2
      real*8  :: t8
      !
      ! Read the wave boundary conditions, unless couple_snapwave has already done so
      !
      call read_snapwave_boundary_data(.not. snapwave)
      !
      if (bathtub_snapwave_nwbnd == 0) then
         !
         ! No wave boundary conditions given: no wave contribution to the bathtub water levels
         !
         call write_log('Warning : bathtub mode with SnapWave but without wave boundary conditions, wave contribution to bathtub water levels is turned off', 1)
         !
         return
         !
      endif
      !
      ! Determine weights and indices for each wave boundary point
      !
      allocate(bathtub_snapwave_i1(nbnd))
      allocate(bathtub_snapwave_i2(nbnd))
      allocate(bathtub_snapwave_w1(nbnd))
      !
      do ib = 1, nbnd
         !
         call interp_segment(bathtub_snapwave_x_bwv, bathtub_snapwave_y_bwv, bathtub_snapwave_nwbnd, x_bnd(ib), y_bnd(ib), i1, i2, w1, w2)
         !
         bathtub_snapwave_i1(ib) = i1
         bathtub_snapwave_i2(ib) = i2
         bathtub_snapwave_w1(ib) = w1
         !
      enddo
      !
      ! Loop through times in bzs file and add f*Hm0 (typically f=0.2)
      !
      do itb = 1, ntbnd
         !
         t8 = t_bnd(itb)
         !
         call update_snapwave_boundary_data(t8)
         !
         do ib = 1, nbnd
            !
            h1 = bathtub_snapwave_hs_bwv(bathtub_snapwave_i1(ib))
            h2 = bathtub_snapwave_hs_bwv(bathtub_snapwave_i2(ib))
            w1 = bathtub_snapwave_w1(ib)
            w2 = 1.0 - w1
            !
            zs_bnd(ib, itb) = zs_bnd(ib, itb) + bathtub_fac_hs * (w1 * h1 + w2 * h2)
            !
         enddo
         !
      enddo
      !
   end subroutine initialize_bathtub_wave_setup
   !
   !
   subroutine bathtub_compute_water_levels()
   !
   ! Set the water level in every cell from the weighted boundary water levels.
   !
   ! Called from: sfincs_update (sfincs_lib)
   !
   use sfincs_data
   use sfincs_timers
   use geometry
   !
   implicit none
   !
   integer :: nm, i1, i2
   real*4  :: zbt, w1, w2
   !
   call timer_start('continuity')
   !
   !$omp parallel &
   !$omp private ( nm, i1, i2, w1, w2, zbt )
   !$omp do schedule ( dynamic, 256 )
   do nm = 1, np
      !
      i1 = bathtub_i1(nm)
      i2 = bathtub_i2(nm)
      w1 = bathtub_w1(nm)
      w2 = 1.0 - w1
      !
      zbt = w1 * zst_bnd(i1) + w2 * zst_bnd(i2)
      !
      if (subgrid) then
         !
         zs(nm) = max(subgrid_z_zmin(nm), zbt)
         !
      else
         !
         zs(nm) = max(zb(nm), zbt)
         !
      endif
      !
      if (store_maximum_waterlevel) then
         !
         ! Store the maximum water level
         !
         zsmax(nm) = max(zsmax(nm), zs(nm))
         !
      endif      
      !
   enddo   
   !$omp end do
   !$omp end parallel
   !
   !$acc update device( zs, zsmax )
   !
   call timer_stop('continuity')
   !
   end subroutine
   !
   !
   subroutine read_snapwave_boundary_data(read_files)
      !
      ! Copy the SnapWave wave boundary points for use in bathtub mode, after reading the
      ! SnapWave input and wave boundary conditions when they have not been read yet.
      !
      ! Input:
      !   read_files : read the SnapWave input and wave boundary conditions (false when
      !                couple_snapwave has already read them)
      !
      ! Called from: initialize_bathtub_wave_setup (sfincs_bathtub)
      !
      use snapwave_data
      use snapwave_boundaries
      use snapwave_input, only: read_snapwave_input
      !
      implicit none
      !
      logical, intent(in) :: read_files
      !
      bathtub_snapwave_nwbnd = 0
      !
      if (read_files) then
         !
         ! Read SnapWave input file from sfincs.inp (this will store boundary file names in snapwave_data)
         !
         call read_snapwave_input()
         !
         ! Without wave boundary conditions there is nothing to read (wave part is turned off in initialize_bathtub_wave_setup)
         !
         if (snapwave_jonswapfile == 'none' .and. netsnapwavefile == 'none' .and. snapwave_bndfile(1:4) == 'none') return
         !
         call read_boundary_data()
         !
      endif
      !
      if (nwbnd == 0) return
      !
      allocate(bathtub_snapwave_x_bwv(nwbnd))
      allocate(bathtub_snapwave_y_bwv(nwbnd))
      allocate(bathtub_snapwave_hs_bwv(nwbnd))
      !
      bathtub_snapwave_nwbnd = nwbnd
      bathtub_snapwave_x_bwv = x_bwv
      bathtub_snapwave_y_bwv = y_bwv
      !
   end subroutine read_snapwave_boundary_data
   !
   !
   subroutine update_snapwave_boundary_data(tb)
      !
      ! Interpolate the SnapWave boundary wave heights to time tb. The time index of the
      ! SnapWave boundary conditions is restored afterwards, so that a full SnapWave run that
      ! follows starts its search from the beginning of the series.
      !
      ! Input:
      !   tb : time to interpolate the wave boundary conditions to [s since tref]
      !
      ! Called from: initialize_bathtub_wave_setup (sfincs_bathtub)
      !
      use snapwave_data
      use snapwave_boundaries
      !
      implicit none
      !
      real*8  :: tb
      !
      integer :: itwbndlast_saved
      !
      itwbndlast_saved = itwbndlast
      !
      call update_boundary_points(tb, .true.)
      !
      itwbndlast = itwbndlast_saved
      !
      bathtub_snapwave_hs_bwv = hst_bwv
      !
   end subroutine update_snapwave_boundary_data
   
end module
