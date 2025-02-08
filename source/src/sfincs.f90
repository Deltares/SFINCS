program sfincs
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
