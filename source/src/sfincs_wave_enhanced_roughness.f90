module sfincs_wave_enhanced_roughness
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! Wave-enhanced bed roughness: apparent roughness at the uv points due to wave-current
   ! interaction.
   !
   ! Subroutines / functions:
   !
   !   update_wave_enhanced_roughness()
   !     Compute the apparent roughness gnapp2 at the uv points from the wave conditions.
   !     Called from: update_apparent_roughness (sfincs_bmi)
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

contains

   subroutine update_wave_enhanced_roughness()
   !
   ! Compute the apparent roughness gnapp2 at the uv points from the wave conditions.
   !
   ! Called from: update_apparent_roughness (sfincs_bmi)
   !
   use sfincs_data
   !
   implicit none
   !
   integer :: ip, nm, nmu
   real*4  :: Uw, Uc, zsu, hu, n_app, n_base, cd, cdeff, uu, vu
   !
   if (.not. wave_enhanced_roughness) return
   !
   !$omp parallel &
   !$omp private ( ip, nm, nmu, Uw,uu, vu, Uc, zsu, hu, n_base, cd, cdeff, n_app )
   !$omp do schedule ( dynamic, 256 )
   do ip = 1, npuv
      !
      if (kcuv(ip) == 1 .or. kcuv(ip) == 6) then
         !
         ! Regular UV point (or a coastal lateral boundary point)
         !
         ! Indices of surrounding water level points
         !
         nm  = uv_index_z_nm(ip)
         nmu = uv_index_z_nmu(ip)
         !
         gnapp2(ip) = subgrid_uv_navg_w(ip) ! subgrid_uv_navg_w(ip) is already converted to g * n**2 !
         !
         Uw = 0.5 * (uorb(nm) + uorb(nmu))
         ! 
         uu = uv(ip)  
         vu = (uv(uv_index_v_ndm(ip)) + uv(uv_index_v_ndmu(ip)) + uv(uv_index_v_nm(ip)) + uv(uv_index_v_nmu(ip))) / 4
         !
         Uc = max(sqrt(uu**2 + vu**2) , 0.25)
         !
         if (Uw < 0.1) cycle
         !
         zsu = max(zs(nm), zs(nmu))         
         hu = subgrid_uv_havg_zmax(ip) + zsu
         !
         if (hu < 0.1) cycle
         !
         ! Base friction factor from Manning's n
         !
         n_base = sqrt(subgrid_uv_navg_w(ip) / g)
         !
         cd = g * n_base**2 / hu**(1.0 / 3.0)
         !
         ! Use Ruessink (2001) formulation
         !
         cdeff = cd * sqrt((1.16 * Uw)**2 + Uc**2) / Uc
         !
         ! Use the Grant & Madsen (1979) or Soulsby (1997) formulation instead
         !
         ! cdeff = cd * (1.4 * Uw + Uc) / Uc
         !
         n_app = sqrt(cdeff * hu**(1.0/3.0) / g)
         !                  
         gnapp2(ip) = g * n_app**2
         !
      endif
      !
   enddo
   !$omp end do
   !$omp end parallel
   !
   ! gnapp2 is only computed on the CPU, so must update GPU memory now
   !
   !$acc update device(gnapp2)
   !
   end subroutine
   
end module
