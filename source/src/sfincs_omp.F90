module sfincs_omp
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! OpenMP runtime settings.
   !
   ! Subroutines / functions:
   !
   !   set_openmp_wait_policy(message)
   !     Keep the worker threads spinning between parallel regions (Intel runtime).
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
   implicit none
   !
contains
   !
   subroutine set_openmp_wait_policy(message)
      !
      ! SFINCS enters thousands of short parallel regions per second. The Intel
      ! OpenMP runtime puts its worker threads to sleep KMP_BLOCKTIME (200 ms by
      ! default) after a region and wakes them through the operating system for
      ! the next one; on small grids that costs more than the work itself.
      ! Unless the user set KMP_BLOCKTIME or OMP_WAIT_POLICY, the threads are
      ! told to keep spinning, the equivalent of KMP_BLOCKTIME=infinite. This is
      ! a property of the runtime library the compiler links (libiomp5 for ifort
      ! and ifx on every platform), not of the processor or the operating system.
      ! gfortran's libgomp has no such call; set OMP_WAIT_POLICY=active there.
      ! Spinning threads occupy their cores, so run one thread per physical core.
      !
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
      use omp_lib, only: kmp_set_blocktime
#endif
      !
      implicit none
      !
      character(len=*), intent(out) :: message
      !
      integer :: status_blocktime
      integer :: status_wait_policy
      !
      ! 2 000 000 ms (33 minutes) is the largest value the Intel runtime accepts without an overflow
      ! notice (it stores microseconds in a 32-bit integer); no gap between two parallel regions of a
      ! SFINCS run comes anywhere near it, so the threads never sleep
      !
      integer, parameter :: blocktime_ms = 2000000
      !
      call get_environment_variable('KMP_BLOCKTIME',   status=status_blocktime)
      call get_environment_variable('OMP_WAIT_POLICY', status=status_wait_policy)
      !
      if (status_blocktime == 0 .or. status_wait_policy == 0) then
         !
         message = '---- OpenMP wait policy taken from the environment (KMP_BLOCKTIME / OMP_WAIT_POLICY) ----'
         !
      else
         !
#if defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
         call kmp_set_blocktime(blocktime_ms)
         message = '---- OpenMP threads keep spinning between parallel regions (KMP_BLOCKTIME 2000000 ms) ----'
#else
         message = '---- OpenMP wait policy left to the runtime (set OMP_WAIT_POLICY=active for small grids) ----'
#endif
         !
      endif
      !
   end subroutine set_openmp_wait_policy
   !
end module sfincs_omp
