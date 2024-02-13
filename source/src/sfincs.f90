program sfincs
   !
   use sfincs_lib
   !
   implicit none
   !
   integer :: ierr
   character(len=1024) :: config_file
   double precision    :: deltat
   !
   deltat = -1.0
   ierr = 0
   !
   if (ierr==0) then
      ierr = sfincs_initialize(config_file)
   endif
   !
   if (ierr==0) then
      ierr = sfincs_update(deltat)
   endif
   ! Output:
   !ierr = 0 ! simulation ran succesfully
   !ierr = 1 ! simulation was stopped because of instabilities   
   !
   ! Always finalize, especially in case of error:
   !if (ierr==0) then
   ierr = sfincs_finalize()
   !endif
   !
end program
