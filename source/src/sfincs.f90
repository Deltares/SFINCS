program sfincs
   !
   ! Main program of the SFINCS executable. Runs a complete simulation by calling
   ! sfincs_initialize, sfincs_update and sfincs_finalize from sfincs_lib, and exits with
   ! the error code they return.
   !
   use sfincs_data 
   use sfincs_lib
   !
   implicit none
   !
   integer :: ierr
   double precision    :: deltat
   !
   deltat = -999.9
   ierr = 0
   !
   if (ierr == 0) then
      !
      ! Set BMI flags to false 
      ! 
      bmi = .false. 
      use_qext = .false.
      ! 
      ierr = sfincs_initialize()
      ! 
   endif
   !
   if (ierr == 0) then
      ! 
      ! if deltat < 999.0 update until end
      !
      ierr = sfincs_update(deltat)
      ! 
   endif
   !
   ! Always finalize, especially in case of error
   !
   ierr = sfincs_finalize()
   !
   call exit(ierr)
   !
end program
