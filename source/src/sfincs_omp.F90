module sfincs_omp
   !
   ! OpenMP runtime settings.
   !
   ! Subroutines / functions:
   !
   !   set_openmp_wait_policy(message)
   !     Keep the worker threads spinning between parallel regions (Intel runtime).
   !     Called from: sfincs_initialize (sfincs_lib)
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
      ! the next one; on small grids that costs more than the work itself
      ! (Bellinge, 18 000 cells and 1020 1D nodes: 70 to 110 s instead of 12 s).
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
