module sfincs_output

   use sfincs_ncoutput
   
   contains

   subroutine initialize_output(tmapout,tmaxout,thisout, trstout)
   !
   use sfincs_data
   !
   implicit none
   !
   real*8   :: tmapout
   real*8   :: tmaxout
   real*8   :: trstout
   real*8   :: thisout
   !
   if (dtmapout>1.0e-6) then
      tmapout     = t0out
      if (outputtype_map == 'net') then
         if (use_quadtree) then
            call ncoutput_quadtree_map_init()
         else
            call ncoutput_regular_map_init()
         endif         
      else    
         call open_map_output()       
      endif      
   else
       tmapout     = 1.0e9
   endif
   !
   if (dtmaxout>1.0e-6) then
      tmaxout     = t0out + dtmaxout
      if (outputtype_map /= 'net') then
         call open_max_output()   ! For netcdf output this is written to mapfile
      endif
   else
      tmaxout     = 1.0e9
   endif
   !
   if (dtrstout>1.0e-6) then
      !
      ! Interval given
      !
      trstout     = t0out + dtrstout
      !
   elseif (trst>0.0) then
      !
      ! Restart time given
      !
      trstout = trst
      !
   else   
      !
      ! No restart file required
      !
      trstout     = 1.0e9
      !
   endif
   !
   if (dthisout>1.0e-6 .and. (nobs>0 .or. nrcrosssections>0)) then
      thisout     = t0
      if (outputtype_his == 'net') then    
         call ncoutput_his_init()
      else      
         call open_his_output()
      endif
   else
     thisout     = 1.0e9
   endif
   !   
   end subroutine

   
   subroutine write_output(t,write_map,write_his,write_max,write_rst,ntmapout,ntmaxout,nthisout,tloop)
   !
   use sfincs_data
   !
   implicit none
   !
   integer  :: count0
   integer  :: count1
   integer  :: count_rate
   integer  :: count_max
   real     :: tloop
   !
   logical  :: write_map   
   logical  :: write_max
   logical  :: write_his   
   logical  :: write_rst   
   !
   integer  :: ntmapout 
   integer  :: ntmaxout
   integer  :: nthisout      
   integer  :: nm
   !
   real*8   :: t
   !
   call system_clock(count0, count_rate, count_max)
   !
   ! Time-varying water level output maps
   !
   ! Update CPU memory
   !
   if (write_map .or. write_his .or. write_rst) then
      !
      !$acc update host(zs)
      !
      if (store_cumulative_precipitation) then
         !      
         !$acc update host(prcp)
         !
      endif   
      !
      if (store_velocity) then
         !      
         !$acc update host(uv)
         !
      endif   
      !
      if (write_rst) then
         !$acc update host(q)
      endif   
      !      
   endif
   !
   if (write_map) then
      !
      if (outputtype_map == 'net') then
         !
         if (use_quadtree) then
            call ncoutput_update_quadtree_map(t, ntmapout)
         else
            call ncoutput_update_regular_map(t, ntmapout)
         endif
         !
      else
         !
         call write_map_output()
         !
      endif
      !
   endif
   !      
   ! Maximum water level output maps
   !
   if (write_max) then
      !
      !$acc update host(zsmax)
      !
      if (store_maximum_velocity) then
         !$acc update host(vmax)
      endif
      !
      if (store_twet) then
         !$acc update host(twet)
      endif
      !
      if (store_cumulative_precipitation) then
         !$acc update host(cumprcp)
      endif
      !
      if (precip .and. scsfile(1:4) /= 'none') then !!! also include other infiltration 
         !$acc update host(cuminf)
      endif
      !
      if (outputtype_map == 'net') then
         !
         if (use_quadtree) then
            call ncoutput_update_quadtree_max(t, ntmaxout)
         else   
            call ncoutput_update_max(t, ntmaxout)
         endif   
         !      
      else
         !      
         call write_max_output()
         !      
      endif
      !
      if (store_maximum_waterlevel) then
         zsmax = -999.0 ! Set zsmax back to a small value
         !$acc update device(zsmax)
      endif
      !
      if (store_maximum_velocity) then
         vmax = 0.0 ! Set vmax back to 0.0
         !$acc update device(vmax)
      endif
      !
!      if (precip .and. store_cumulative_precipitation) then
!         cumprcp = 0.0 ! Set cumprcp back to a 0.0
!         !$acc update device(cumprcp)
!      endif            
      !
      if (store_twet) then
         twet = 0.0 ! Set twet back to 0.0
         !$acc update device(twet)
      endif
      !      
   endif
   !
   !      
   ! Binary restart file
   !
   if (write_rst) then
      !
      call write_rst_file(t)
      !
   endif
   !      
   ! Water level time series
   !
   if (write_his .and. (nobs>0 .or. nrcrosssections>0)) then
      !      
      if (outputtype_his == 'net') then
         !      
         call ncoutput_update_his(t,nthisout)
         !      
      else
         !      
         call write_his_output(t)
         !      
      endif
      !
   endif
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !   
   end subroutine
   
   subroutine finalize_output(t, ntmaxout, tloopoutput)
   !
   use sfincs_data
   !
   implicit none
   !
   integer  :: ntmaxout
   real*8   :: t   
   real     :: tloopoutput 
   !   
   if (dtmaxout>1.e-6 .and. ntmaxout == 0) then
      write(*,'(a)')''       
      write(*,*)'Info : Write maximum values of final timestep since t=dtmaxout was not reached yet...'
      ntmaxout = 1
      call write_output(t,.false.,.false.,.true.,.false.,0,ntmaxout,0,tloopoutput)
   endif
   !
   if (outputtype_map == 'net') then
      !
      call ncoutput_map_finalize()
      !
   else
      !
      call close_map_output()
      !
      call close_max_output()
      !
   endif
   !
   if (outputtype_his == 'net') then 
      !
      call ncoutput_his_finalize()
      !
   endif
   !
   end subroutine
   
   
   subroutine open_map_output()
   !
   use sfincs_data
   !
   implicit none
   !
   if (trim(outputtype_map)=='asc') then
      open(unit = 900, file = 'zs.txt')
      open(unit = 901, file = 'u.txt')
      open(unit = 902, file = 'v.txt')
      if (store_cumulative_precipitation) then
          open(unit = 903, file = 'cumprcp.txt')
      endif      
   else
      open(unit = 900, status = 'replace', file = 'zs.dat', form = 'unformatted')
      !
      if (store_zvolume) then
        open(unit = 908, status = 'replace', file = 'z_volume.dat', form = 'unformatted')
      endif
   endif
   !
   end subroutine

   
   
   subroutine open_max_output()
   !
   use sfincs_data
   !
   implicit none
   !
   if (store_maximum_waterlevel) then
    open(unit = 850, status = 'replace', file = 'zsmax.dat', form = 'unformatted')
   endif
   !
   if (store_maximum_velocity) then
      open(unit = 851, status = 'replace', file = 'velmax.dat', form = 'unformatted')
   endif
   ! 
   if (store_twet) then
      open(unit = 852, status = 'replace', file = 'tmax.dat', form = 'unformatted')
   endif
   ! 
   if (store_cumulative_precipitation) then
      open(unit = 853, status = 'replace', file = 'cumprcp.dat', form = 'unformatted')
   endif
   !
   if (infiltration) then
      open(unit = 854, status = 'replace', file = 'cuminf.dat', form = 'unformatted')
   endif
   !
   end subroutine
   

   subroutine write_map_output()
   !
   use sfincs_data
   !
   implicit none
   !
   integer                      :: nm, n, m
   !
   real*4, dimension(:,:), allocatable :: zsg
   !
!   if (trim(outputtype_map)=='asc') then
!      !
!      allocate(zsg(nmax,mmax))
!      zsg = 0.0
!      !
!      do nm = 1, np
!         n = index_v_n(nm)
!         m = index_v_m(nm)
!         zsg(n, m) = zs(nm)
!      enddo
!      do n = 2, nmax - 1
!         write(900,'(10000f10.4)')(zsg(n,m), m = 2, mmax - 1)
!      enddo             
!      !
!      if (store_cumulative_precipitation) then
!          zsg = 0.0
!          !
!          do nm = 1, np
!             n = index_v_n(nm)
!             m = index_v_m(nm)
!             zsg(n, m) = cumprcp(nm)
!          enddo
!          do n = 2, nmax - 1
!             write(903,'(10000f10.4)')(zsg(n,m), m = 2, mmax - 1)
!          enddo             
!      !          
!      endif    
      !      
!      zsg = 0.0
!      do nm = 1, np
!         zsg(index_v_n(nm),index_v_m(nm)) = u(nm)
!      enddo
!      do n = 2, nmax - 1
!         write(901,'(10000f9.3)')(zsg(n,m), m = 2, mmax - 1)
!      enddo             
!      !      
!      zsg = 0.0
!      do nm = 1, np
!         zsg(index_v_n(nm),index_v_m(nm)) = v(nm)
!      enddo
!      do n = 2, nmax - 1
!         write(902,'(10000f9.3)')(zsg(n,m), m = 2, mmax - 1)
!      enddo
      !
!      if (wind) then
!         do nm = 1, np
!            zsg(index_v_n(nm),index_v_m(nm)) = wndu(nm)
!         enddo
!         do n = 2, nmax - 1
!            write(903,'(10000f9.3)')(zsg(n,m), m = 2, mmax - 1)
!         enddo
!         do nm = 1, np
!            zsg(index_v_n(nm),index_v_m(nm)) = wndv(nm)
!         enddo
!         do n = 2, nmax - 1
!            write(904,'(10000f9.3)')(zsg(n,m), m = 2, mmax - 1)
!         enddo
!      endif
!      !
!   else
      !
      write(900)zs
      !
      if (store_zvolume) then
          write(908)z_volume
      endif
      !
!   endif
   !
   end subroutine


   subroutine write_max_output()
   !
   use sfincs_data
   !
   implicit none
   !
   integer                      :: nm
   !
   if (store_maximum_waterlevel) then
      !
      ! Set zsmax in dry points to -999.0
      !
      if (subgrid) then
         !
         do nm = 1, np
            !
            if (zsmax(nm) < subgrid_z_zmin(nm) + huthresh) then
               zsmax(nm) = -999.0
            endif   
            !
         enddo
         !
      else
         !
         do nm = 1, np
            !
            if (zsmax(nm) < zb(nm) + huthresh) then
               zsmax(nm) = -999.0
            endif   
            !
         enddo
         !
      endif      
      !
      write(850)zsmax
      !
   endif
   !
   if (store_maximum_velocity) then
      write(851)vmax
   endif
   !
   if (store_twet) then
      write(852)twet
   endif
   !
   if (store_cumulative_precipitation) then
      write(853)cumprcp
   endif
   ! 
   if (infiltration) then
      write(854)cuminf
   endif
   ! 
   end subroutine
   
   
   subroutine close_map_output()
   !
   use sfincs_data
   !
   implicit none
   !
   if (dtmapout>1.e-6) then
      close(900)
      close(901)
      close(902)
      !
      if (store_zvolume) then
         close(908)
      endif
   endif   
   !
   end subroutine

   subroutine close_max_output()
   !
   use sfincs_data
   !
   implicit none
   !
   if (store_maximum_waterlevel) then
      close(850)
   endif   
   if (store_maximum_velocity) then
      close(851)
   endif
   !
   if (store_twet) then
      close(852)
   endif
   !
   if (store_cumulative_precipitation) then
      close(853)
   endif
   !
   if (infiltration) then
      close(854)
   endif
   !
   end subroutine

   subroutine open_his_output()
   !
   use sfincs_data
   !
   implicit none
   !
   open(unit = 950, file = 'zst.txt')
   close(unit = 950 ,status='delete')
   if (precip) then
      open(unit = 965, file = 'prcp.txt')
      close(unit = 965 ,status='delete')
   endif
   if (nrcrosssections>0) then
      open(unit = 966, file = trim('qt.txt'))
      close(unit = 966 ,status='delete')
   endif
   if (nsrcdrn>0) then
      open(unit = 970, file = trim('qdrain.txt'))
      close(unit = 970 ,status='delete')
   endif
   !
   ! Delete existing files
   !
   end subroutine


   subroutine write_his_output(t)
   !
   use sfincs_data
   use sfincs_crosssections
   !
   implicit none
   !
   integer                      :: nm
   integer                      :: iobs   
   integer                      :: icrs   
   integer                      :: ip   
   !
   real*8                             :: t
   real*4, dimension(nobs)            :: tprcp
   real*4, dimension(:), allocatable  :: qq
   !
   do iobs = 1, nobs
      nm = nmindobs(iobs)
      if (nm>0) then
         zobs(iobs)  = zs(nm)
         hobs(iobs)  = zs(nm) - zbobs(iobs)                 
         if (precip) then
            tprcp(iobs) = prcp(nm)*1000*3600
         endif   
      else
         zobs(iobs) = -999.0
         hobs(iobs) = -999.0
         if (precip) then
            tprcp(iobs) = -999.0
         endif   
      endif
   enddo
   !
   open(unit = 950, file = trim('zst.txt'), access='append')
   write(950,'(f12.1,10000f9.3)')t,(zobs(iobs), iobs = 1, nobs)
   close(950)
   if (precip) then
      open(unit = 965, file = trim('prcp.txt'), access='append')
      write(965,'(f12.1,10000f9.3)')t,(tprcp(iobs), iobs = 1, nobs)
      close(965)
   endif
   !   
   if (nrcrosssections>0) then
      !
      !$acc update host(q)
      !      
      ! Determine fluxes through cross sections
      !
      call get_discharges_through_crosssections(qq)
      !
      open(unit = 966, file = trim('qt.txt'), access='append')
      write(966,'(f12.1,10000f12.3)')t,(q(icrs), icrs = 1, nrcrosssections)
      close(966)
      !
   endif
   !
   if (ndrn>0 .and. store_qdrain) then
      !$acc update host(qtsrc)
      open(unit = 970, file = trim('qdrain.txt'), access='append')
      write(970,'(f12.1,10000f9.3)')t,(qtsrc(iobs), iobs = nsrc + 1, nsrcdrn, 2)
      close(970)
   endif
   !
   end subroutine
   
   subroutine write_rst_file(t)
   !
   use sfincs_data
   use sfincs_date
   !
   implicit none
   !
   character*256 :: file_name
   character*15  :: tstring
   real*8        :: t
   real*8        :: tt
   !
   tt = 1.0*int(t)
   tstring = time_to_string(tt, trefstr)
   write(file_name,'(A,A,A)')'sfincs.',tstring,'.rst'
   !
   open(unit = 911, status = 'replace', file = trim(file_name), form = 'unformatted')
   write(911)1    ! 1: zs, qx, qy, umean and vmean  - 2: zs, qx, qy - 3: zs 
   write(911)zs
   write(911)q
   write(911)uvmean
   close(911)
   !
   end subroutine

   
   subroutine write_tsunami_arrival_file()
   !
   use sfincs_data
   !
   implicit none
   !
   open(unit = 911, status = 'replace', file = 'tsunami_arrival_time.dat', form = 'unformatted')
   write(911)tsunami_arrival_time
   close(911)   
   !
   end subroutine
   
end module
