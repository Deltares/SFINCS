program sfincs
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! Main program of the SFINCS executable. Runs a complete simulation by calling
   ! sfincs_initialize, sfincs_update and sfincs_finalize from sfincs_lib, and exits with
   ! the error code they return.
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
      use_dzbext = .false.
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
