module sfincs_source_revision
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! Source revision the code was built from (Source-Revision in sfincs.log and in the
   ! netCDF output). The autotools build generates sfincs_source_revision.inc with the
   ! output of git describe (or the SFINCS_GIT_REVISION override) and defines
   ! HAVE_SOURCE_REVISION_INC, see src/Makefile.am. Other builds report 'unknown'.
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
#ifdef HAVE_SOURCE_REVISION_INC
#include "sfincs_source_revision.inc"
#else
   character(len=*), parameter :: build_source_revision = 'unknown'
#endif
   !
end module sfincs_source_revision
