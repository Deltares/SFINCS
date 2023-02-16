module sfincs_snapwave
   !
   implicit none
   !     
   integer                                   :: snapwave_no_nodes
   integer                                   :: snapwave_no_cells
   real*8,    dimension(:),   allocatable    :: snapwave_x
   real*8,    dimension(:),   allocatable    :: snapwave_y
   real*4,    dimension(:),   allocatable    :: snapwave_z
   real*4,    dimension(:),   allocatable    :: snapwave_depth
   real*4,    dimension(:),   allocatable    :: snapwave_H
   real*4,    dimension(:),   allocatable    :: snapwave_H_ig
   real*4,    dimension(:),   allocatable    :: snapwave_mean_direction
   real*4,    dimension(:),   allocatable    :: snapwave_directional_spreading
   real*4,    dimension(:),   allocatable    :: snapwave_Fx
   real*4,    dimension(:),   allocatable    :: snapwave_Fy
   integer,   dimension(:,:), allocatable    :: snapwave_connected_nodes
   integer*4, dimension(:),   allocatable    :: index_snapwave_in_sfincs
   integer*4, dimension(:),   allocatable    :: index_sfincs_in_snapwave
   integer*4, dimension(:),   allocatable    :: index_sw_in_qt ! used in sfincs_ncoutput (copy of index_snapwave_in_quadtree from snapwave_data)
   real*4                                    :: snapwave_tpmean
   !
contains
   !
   subroutine couple_snapwave()
   !
   use snapwave_data
   use snapwave_domain
   use snapwave_boundaries
   !
   implicit none
   !
   integer :: ipsw, iq
   !
   call read_snapwave_input()            ! Reads snapwave.inp
   !
   call initialize_snapwave_domain()     ! Read mesh, finds upwind neighbors, etc.
   !
   call read_boundary_data()
   !
   snapwave_no_nodes = no_nodes
   !
   allocate(snapwave_z(no_nodes))
   allocate(snapwave_depth(no_nodes))
   !
   snapwave_z     = zb
   snapwave_depth = 0.0
   !
   snapwave_tpmean = 0.0
   !   
   call find_matching_cells(index_quadtree_in_snapwave, index_snapwave_in_quadtree)
   !
   ! No longer need any of these
   !
   end subroutine
   

   subroutine find_matching_cells(index_quadtree_in_snapwave, index_snapwave_in_quadtree)
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   integer, dimension(snapwave_no_nodes),  intent(in) :: index_quadtree_in_snapwave
   integer, dimension(quadtree_nr_points), intent(in) :: index_snapwave_in_quadtree
   !
   integer :: nm, n, m, ipsw, ipsf, iq
   !
   allocate(index_sfincs_in_snapwave(snapwave_no_nodes))
   allocate(index_snapwave_in_sfincs(np))
   allocate(index_sw_in_qt(quadtree_nr_points))
   !
   index_sfincs_in_snapwave = 0
   index_snapwave_in_sfincs = 0
   index_sw_in_qt = 0
   !
   ! Loop through SnapWave points
   !
   do ipsw = 1, snapwave_no_nodes
      iq   = index_quadtree_in_snapwave(ipsw)
      ipsf = index_sfincs_in_quadtree(iq)
      index_sfincs_in_snapwave(ipsw) = ipsf
      index_sw_in_qt(iq) = ipsw
   enddo   
   !
   ! Loop through SFINCS points
   !
   do ipsf = 1, np
      iq   = index_quadtree_in_sfincs(ipsf)
      ipsw = index_snapwave_in_quadtree(iq)
      index_snapwave_in_sfincs(ipsf) = ipsw
   enddo   
   !
   end subroutine

   
   subroutine update_wave_field(t, tloop)
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
   real*4,    dimension(snapwave_no_nodes)    :: zssw
   real*4,    dimension(snapwave_no_cells, 4) :: hc
   real*4,    dimension(:), allocatable       :: fwx0
   real*4,    dimension(:), allocatable       :: fwy0
   integer   :: ip, ii, m, n, nm, nmu, idir
   real*4    :: f
   real*8    :: t
   !
   call system_clock(count0, count_rate, count_max)
   !
   allocate(fwx0(np))
   allocate(fwy0(np))
   !
   fwx0 = 0.0
   fwy0 = 0.0
   !
   ! Determine SnapWave water depth
   !
   do nm = 1, snapwave_no_nodes
      !
      ip = index_sfincs_in_snapwave(nm) ! matching index in SFINCS mesh
      !
      if (ip>0) then
         !
         ! A matching SFINCS point is found
         !
         if (wavemaker) then
            !
            snapwave_depth(nm) = max(zsm(ip) - snapwave_z(nm), 0.01)      
            !
         else   
            !
            snapwave_depth(nm) = max(zs(ip) - snapwave_z(nm), 0.01)      
            !
         endif   
         !
      else
         !
         ! Use 0.0 water level
         !
         snapwave_depth(nm) = max(0.0 - snapwave_z(nm), 0.1)      
         !
      endif   
      !
   enddo   
   !
   call compute_snapwave(t)
   !
   do nm = 1, np
      !
      ip = index_snapwave_in_sfincs(nm) ! matching index in SFINCS mesh
      !
      if (ip>0) then
         !
         hm0(nm)    = snapwave_H(ip)   
         hm0_ig(nm) = snapwave_H_ig(ip)   
         fwx0(nm)   = snapwave_Fx(ip)   
         fwy0(nm)   = snapwave_Fy(ip)   
         if (store_wave_direction) then
            mean_wave_direction(nm)        = 270.0 - snapwave_mean_direction(ip)*180/pi   
            wave_directional_spreading(nm) = snapwave_directional_spreading(ip)*180/pi   
         endif
         !
      else
         !
         ! SnapWave point outside active SFINCS domain
         !
         hm0(nm)    = 0.0
         hm0_ig(nm) = 0.0
         fwx0(nm)   = 0.0
         fwy0(nm)   = 0.0   
         if (store_wave_direction) then
            mean_wave_direction(nm)        = 0.0
            wave_directional_spreading(nm) = 0.0  
         endif
         !
      endif   
      !
      if (store_wave_forces) then
         !
         fwx(nm) = fwx0(nm)
         fwy(nm) = fwy0(nm)
         !
      endif   
      !
   enddo   
   !   
   hm0 = hm0*sqrt(2.0)
   hm0_ig = hm0_ig*sqrt(2.0)
   !
   do ip = 1, npuv
      !
      nm   = uv_index_z_nm(ip)
      nmu  = uv_index_z_nmu(ip)
      idir = uv_flags_dir(ip) ! 0 is u, 1 is v
      !
      ! Should do better averaging for uv points that go from fine to coarse
      !
      if (idir == 0) then
         !
         ! U point
         !         
         fwuv(ip) = (0.5*(cosrot*fwx0(nm) + sinrot*fwy0(nm)) + 0.5*( cosrot*fwx0(nmu) + sinrot*fwy0(nmu)))/rhow
         !         
      else
         !
         ! V point
         !         
         fwuv(ip) = (0.5*(-sinrot*fwx0(nm) + cosrot*fwy0(nm)) + 0.5*(-sinrot*fwx0(nmu) + cosrot*fwy0(nmu)))/rhow
         !         
      endif   
      !
!      fwuv = 0.0
      !
   enddo
   !
   !$acc update device(fwuv), async(1)
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine


   subroutine compute_snapwave(t)
   !
   use snapwave_data
   use snapwave_solver
   use snapwave_boundaries
   !
   real*8    :: t
   !
   call update_boundary_conditions(t) ! SnapWave boundary conditions
   !
   ! Mean wave height for writing to netcdf his file
   !
   snapwave_tpmean = tpmean_bwv
   !
   depth = snapwave_depth
   !
   call compute_wave_field()
   !
   snapwave_H                     = H
   snapwave_H_ig                  = H_ig
   snapwave_mean_direction        = thetam
   snapwave_directional_spreading = thetam
   snapwave_Fx                    = Fx
   snapwave_Fy                    = Fy   
   !
   end subroutine

   
   
   subroutine read_snapwave_input()
   !
   ! Reads snapwave data from sfincs.inp
   !
   use snapwave_data   
   !
   implicit none
   !
   integer :: irestart, iig
   !
   open(500, file='sfincs.inp')   
   !
   ! Input section
   !
   call read_real_input(500,'snapwave_gamma',gamma,0.7)
   call read_real_input(500,'snapwave_alpha',snapwave_alpha,1.0)
   call read_real_input(500,'snapwave_hmin',hmin,0.1)
   call read_real_input(500,'snapwave_fw',fw0,0.01)
   call read_real_input(500,'snapwave_fwig',fw0_ig,0.015)
   call read_real_input(500,'snapwave_dt',dt,36000.0)
   call read_real_input(500,'snapwave_tol',tol,10.0)
   call read_real_input(500,'snapwave_dtheta',dtheta,10.0)
   call read_real_input(500,'snapwave_crit',crit,0.01)
   call read_int_input(500,'snapwave_igwaves',iig,1)
   call read_int_input(500,'snapwave_nrsweeps',nr_sweeps,1)
!   call read_int_input(500,'ntheta',ntheta,36)
!   call read_int_input(500,'nHrel',nHrel,1)
!   call read_char_input(500,'hhtabname',hhtabname,'')
!   call read_char_input(500,'Htabname',Htabname,'')
!   call read_char_input(500,'Dwtabname',Dwtabname,'')
!   call read_char_input(500,'Ftabname',Ftabname,'')
!   call read_char_input(500,'Cgtabname',Cgtabname,'')
!   call read_char_input(500,'cthetafactabname',cthetafactabname,'')
!   call read_char_input(500,'waterlevelfile',waterlevelfile,'')
   call read_char_input(500,'snapwave_jonswapfile',jonswapfile,'')
   call read_char_input(500,'snapwave_bndfile',bndfile,'')
   call read_char_input(500,'snapwave_encfile',encfile,'')
   call read_char_input(500,'snapwave_bhsfile',bhsfile,'')
   call read_char_input(500,'snapwave_btpfile',btpfile,'')
   call read_char_input(500,'snapwave_bwdfile',bwdfile,'')
   call read_char_input(500,'snapwave_bdsfile',bdsfile,'') 
   call read_char_input(500,'snapwave_upwfile',upwfile,'')
   call read_char_input(500,'snapwave_mskfile',mskfile,'')
   call read_char_input(500,'snapwave_depfile',depfile,'none')   
   call read_char_input(500,'snapwave_ncfile', gridfile,'snapwave_net.nc')   
   !
   close(500)
   !
   igwaves             = .true.
   if (iig==0) then
      igwaves = .false.
      write(*,*)'SnapWave: IG waves turned OFF!'
   else
      write(*,*)'SnapWave: IG waves turned ON!'
   endif
   !
   restart           = .true.
   coupled_to_sfincs = .true.
   bzsfile           = ''
   !
   end subroutine 

   
    
   subroutine read_real_input(fileid,keyword,value,default)
   !
   character(*), intent(in) :: keyword
   character(len=256)       :: keystr
   character(len=256)       :: valstr
   character(len=256)       :: line
   integer, intent(in)      :: fileid
   real*4, intent(out)      :: value
   real*4, intent(in)       :: default
   integer j,stat
   !
   value = default
   rewind(fileid)   
   do while(.true.)
      read(fileid,'(a)',iostat = stat)line
      if (stat<0) exit
      j=index(line,'=')      
      keystr = trim(line(1:j-1))
      if (trim(keystr)==trim(keyword)) then
         valstr = trim(line(j+1:256))
         read(valstr,*)value
         exit
      endif
   enddo 
   !
   end  subroutine  

   subroutine read_real_array_input(fileid,keyword,value,default,nr)
   !
   character(*), intent(in) :: keyword
   character(len=256)       :: keystr
   character(len=256)       :: valstr
   character(len=256)       :: line
   integer, intent(in)      :: fileid
   integer, intent(in)      :: nr
   real*4, dimension(:), intent(out), allocatable :: value
   real*4, intent(in)       :: default
   integer j,stat, m
   !
   allocate(value(nr))
   !
   value = default
   rewind(fileid)   
   do while(.true.)
      read(fileid,'(a)',iostat = stat)line
      if (stat<0) exit
      j=index(line,'=')      
      keystr = trim(line(1:j-1))
      if (trim(keystr)==trim(keyword)) then
         valstr = trim(line(j+1:256))
         read(valstr,*)(value(m), m = 1, nr)
         exit
      endif
   enddo 
   !
   end  subroutine  

   
   subroutine read_int_input(fileid,keyword,value,default)
   !
   character(*), intent(in) :: keyword
   character(len=256)       :: keystr
   character(len=256)       :: valstr
   character(len=256)       :: line
   integer, intent(in)      :: fileid
   integer, intent(out)     :: value
   integer, intent(in)      :: default
   integer j,stat
   !
   value = default
   rewind(fileid)   
   do while(.true.)
      read(fileid,'(a)',iostat = stat)line
      if (stat<0) exit
      j=index(line,'=')      
      keystr = trim(line(1:j-1))
      if (trim(keystr)==trim(keyword)) then
         valstr = trim(line(j+1:256))
         read(valstr,*)value         
         exit
      endif
   enddo 
   !
   end subroutine

   
   subroutine read_char_input(fileid,keyword,value,default)
   !
   character(*), intent(in)  :: keyword
   character(len=256)        :: keystr
   character(len=256)        :: valstr
   character(len=256)        :: line
   integer, intent(in)       :: fileid
   character(*), intent(in)  :: default
   character(*), intent(out) :: value
   integer j,stat
   !
   value = default
   rewind(fileid)   
   do while(.true.)
      read(fileid,'(a)',iostat = stat)line
      if (stat<0) exit
      j=index(line,'=')      
      keystr = trim(line(1:j-1))
      if (trim(keystr)==trim(keyword)) then
         valstr = adjustl(trim(line(j+1:256)))
         value = valstr
         exit
      endif
   enddo 
   !
   end subroutine 
   
end module
