module sfincs_wave_enhanced_roughness

contains

   subroutine update_wave_enhanced_roughness()
   !
   use sfincs_data
   !
   implicit none
   !
   integer :: ip, nm, nmu
   real*4  :: Uw, Uc, zsu, hu, n_app, n_base, cd, cdeff, fenh, uu, vu
   !
   if (.not. wave_enhanced_roughness) return
   !
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
         Uc = max(sqrt(uu*2 + vu**2) , 0.25)
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
         cdeff = cd * sqrt((1.16*Uw)**2 + Uc**2) / Uc
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
   !
   end subroutine
   
end module
