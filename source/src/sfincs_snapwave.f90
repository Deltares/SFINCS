module sfincs_snapwave
   !
   implicit none
   !     
   integer                                   :: snapwave_no_nodes
   integer                                   :: snapwave_no_cells
   real*8,    dimension(:),   allocatable    :: snapwave_x
   real*8,    dimension(:),   allocatable    :: snapwave_y
   real*4,    dimension(:),   allocatable    :: snapwave_z
   real*4,    dimension(:),   allocatable    :: snapwave_mask   
   real*4,    dimension(:),   allocatable    :: snapwave_depth
   real*4,    dimension(:),   allocatable    :: snapwave_H
   real*4,    dimension(:),   allocatable    :: snapwave_H_ig
   real*4,    dimension(:),   allocatable    :: snapwave_Tp
   real*4,    dimension(:),   allocatable    :: snapwave_Tp_ig   
   real*4,    dimension(:),   allocatable    :: snapwave_mean_direction
   real*4,    dimension(:),   allocatable    :: snapwave_directional_spreading
   real*4,    dimension(:),   allocatable    :: snapwave_u10
   real*4,    dimension(:),   allocatable    :: snapwave_u10dir   
   real*4,    dimension(:),   allocatable    :: snapwave_Fx
   real*4,    dimension(:),   allocatable    :: snapwave_Fy
   real*4,    dimension(:),   allocatable    :: snapwave_Dw
   real*4,    dimension(:),   allocatable    :: snapwave_Df 
   real*4,    dimension(:),   allocatable    :: snapwave_Dwig
   real*4,    dimension(:),   allocatable    :: snapwave_Dfig
   real*4,    dimension(:),   allocatable    :: snapwave_cg
   real*4,    dimension(:),   allocatable    :: snapwave_Qb
   real*4,    dimension(:),   allocatable    :: snapwave_beta
   real*4,    dimension(:),   allocatable    :: snapwave_srcig
   real*4,    dimension(:),   allocatable    :: snapwave_alphaig   
   integer,   dimension(:,:), allocatable    :: snapwave_connected_nodes
   integer*4, dimension(:),   allocatable    :: index_snapwave_in_sfincs
   integer*4, dimension(:),   allocatable    :: index_sfincs_in_snapwave
   integer*4, dimension(:),   allocatable    :: index_sw_in_qt ! used in sfincs_ncoutput (copy of index_snapwave_in_quadtree from snapwave_data)
   real*4                                    :: snapwave_tpmean
   real*4                                    :: snapwave_tpigmean   
   !
contains
   !
   subroutine couple_snapwave(crsgeo)
   !
   use snapwave_data
   use snapwave_domain
   use snapwave_boundaries
   !
   implicit none
   !
   logical       :: crsgeo
   !
   build_revision = '$Rev: svn 108-branch:SnapWave_IG'
   build_date     = '$Date: 2024-10-16'
   !
   write(*,'(a)')''
   write(*,*)'----------- Welcome to SnapWave ---------'   
   write(*,'(a)')''
   write(*,*)'   @@@@@   @@  @@  @@@@@@  @@@@@@   @@@  '
   write(*,*)'  @@@ @@@  @@@ @@  @@@@@@  @@@@@@   @@@  '
   write(*,*)'  @@@      @@@ @@  @@  @@  @@  @@   @@@  '
   write(*,*)'   @@@@@   @@@@@@  @@@@@@  @@@@@@   @@@  '
   write(*,*)'      @@@  @@ @@@  @@  @@  @@            '
   write(*,*)'  @@@ @@@  @@  @@  @@  @@  @@       @@@  '
   write(*,*)'   @@@@@   @@   @  @@  @@  @@       @@@  '
   write(*,'(a)')''   
   write(*,*)'             .......:.......             '
   write(*,*)'         ...:::::::::::::::::...         '
   write(*,*)'      ..:::::::............::::::..      '
   write(*,*)'    ..::::::.....:@@@@@@@@....:::::..    '
   write(*,*)'   .::::::...~@@@@@@@@@@@@@@~..::::::.   '
   write(*,*)'  .::::::..:@@@@@@@@@@@@@@@@@@:.::::::.  '
   write(*,*)' .:::::..:@@@@@@@@@@@@@@@@@@@@@:.::::::. '
   write(*,*)'.::::..:@@@@@@@@@@@@@@^......:@@.:::::::.'
   write(*,*)'.::...:@@@@@@@@@@@@@@@.:::::..^^.:::::::.'
   write(*,*)'::.:@@@@@@@@@@@@@@@@@@..::::::..:::::::::'
   write(*,*)'..:@@@@@@@@@@@@@@@@@@@@^..............::.'
   write(*,*)'..:@@@@@@@@@@@@@@@@@@@@@@@^:..:~^~^~:..:.'
   write(*,*)' .:@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@:. '
   write(*,*)'  .@@~^~@@@@~^~@@@@~^~@@@@~^~@@@@~^~@@.  '
   write(*,*)'   ...................................   '
   write(*,*)'    ..:::::::::::::::::::::::::::::..    '
   write(*,*)'      ..:::::::::::::::::::::::::..      '
   write(*,*)'         ...:::::::::::::::::...         '
   write(*,*)'             .......:.......             '
   write(*,'(a)')''
   write(*,*)'-----------------------------------------'   
   write(*,'(a)')''
   write(*,*)'Build-Revision: ',trim(build_revision)
   write(*,*)'Build-Date:     ',trim(build_date)
   write(*,'(a)')''   
   !   
   ! Check whether SFINCS grid is spherical (T) or cartesian (F), and prescribe to SnapWave as variable 'sferic' -  spherical (1) or cartesian (0) grid
   if (crsgeo) then
      sferic  = 1 
      write(*,*)'SnapWave: Input grid interpreted as spherical coordinates, sferic= ',sferic      
   endif   
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
   allocate(snapwave_mask(no_nodes))  
   allocate(snapwave_u10(no_nodes))   
   allocate(snapwave_u10dir(no_nodes))      
   !
   snapwave_z     = zb
   snapwave_depth = 0.0
   !
   if (wind) then
      snapwave_u10 = 0.0
      snapwave_u10dir = 0.0
   endif
   !
   snapwave_tpmean = 0.0
   snapwave_tpigmean = 0.0   

   !   
   call find_matching_cells(index_quadtree_in_snapwave, index_snapwave_in_quadtree)
   !
   ! Copy final snapwave mask from snapwave_domain for output in sfincs_ncoutput
   snapwave_mask = msk   
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
   integer :: ipsw, ipsf, iq
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
   real*4   :: u10, u10dir
   !   
   real*4,    dimension(:), allocatable       :: fwx0
   real*4,    dimension(:), allocatable       :: fwy0
   real*4,    dimension(:), allocatable       :: dw0
   real*4,    dimension(:), allocatable       :: df0   
   real*4,    dimension(:), allocatable       :: dwig0
   real*4,    dimension(:), allocatable       :: dfig0   
   real*4,    dimension(:), allocatable       :: cg0   
   real*4,    dimension(:), allocatable       :: qb0   
   real*4,    dimension(:), allocatable       :: beta0 
   real*4,    dimension(:), allocatable       :: srcig0      
   real*4,    dimension(:), allocatable       :: alphaig0   
   integer   :: ip, nm, nmu, idir
   real*8    :: t
   !
   call system_clock(count0, count_rate, count_max)
   !
   allocate(fwx0(np))
   allocate(fwy0(np))
   allocate(dw0(np))
   allocate(df0(np))   
   allocate(dwig0(np))
   allocate(dfig0(np))  
   allocate(cg0(np))  
   allocate(qb0(np))   
   allocate(beta0(np))   
   allocate(srcig0(np))      
   allocate(alphaig0(np))      
   !
   fwx0 = 0.0
   fwy0 = 0.0
   dw0 = 0.0
   df0 = 0.0   
   dwig0 = 0.0
   dfig0 = 0.0
   cg0 = 0.0
   qb0 = 0.0
   beta0 = 0.0
   srcig0 = 0.0
   alphaig0 = 0.0   
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
            snapwave_depth(nm) = max(zsm(ip) - snapwave_z(nm), 0.00001)      
            !
         else   
            !
            snapwave_depth(nm) = max(zs(ip) - snapwave_z(nm), 0.00001)      
            !
         endif   
         !
      else
         !
         ! Use 0.0 water level
         !
         snapwave_depth(nm) = max(0.0 - snapwave_z(nm), 0.00001)      
         !
      endif   
      !
   enddo   
   !
   ! Determine SnapWave wind
   !
   if (store_wind) then
      !
      do nm = 1, snapwave_no_nodes
         !
         ip = index_sfincs_in_snapwave(nm) ! matching index in SFINCS mesh
         !
         if (ip>0) then
            !
            ! A matching SFINCS point is found
            !
            ! Convert to umag & dir, as in ncoutput_update_his: 
            !
            u10 = sqrt(windu(ip)**2 + windv(ip)**2)
            !
            u10dir = atan2(windv(ip), windu(ip))*180/pi
            !
	        if (u10dir<0.0) u10dir = u10dir + 360.0
            if (u10dir>360.0) u10dir = u10dir - 360.0    
            !
            snapwave_u10(nm) = max(u10, 0.0)     
            snapwave_u10dir(nm) = u10dir / 180.0 * pi ! from nautical coming from in degrees to cartesian going to in radians
            !
         else
            !
            ! Use 0.0 wind speed and direction
            !
            snapwave_u10(nm) = 0.0
            snapwave_u10dir(nm) = 0.0            
            !
         endif   
         !
      enddo   
      !
   endif
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
         sw_tp(nm)    = snapwave_Tp(ip)         
         sw_tp_ig(nm) = snapwave_Tp_ig(ip)         
         fwx0(nm)   = snapwave_Fx(ip)   
         fwy0(nm)   = snapwave_Fy(ip) 
         dw0(nm)    = snapwave_Dw(ip)   
         df0(nm)    = snapwave_Df(ip)     
         dwig0(nm)  = snapwave_Dwig(ip)   
         dfig0(nm)  = snapwave_Dfig(ip)
         cg0(nm)    = snapwave_cg(ip)
         qb0(nm)    = snapwave_Qb(ip)
         beta0(nm)  = snapwave_beta(ip)
         srcig0(nm) = snapwave_srcig(ip)
         alphaig0(nm) = snapwave_alphaig(ip)
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
         sw_tp(nm)  = 0.0
         sw_tp_ig(nm) = 0.0         
         fwx0(nm)   = 0.0
         fwy0(nm)   = 0.0   
         dw0(nm)    = 0.0
         df0(nm)    = 0.0         
         dwig0(nm)  = 0.0
         dfig0(nm)  = 0.0
         cg0(nm)    = 0.0
         qb0(nm)    = 0.0
         beta0(nm)  = 0.0
         srcig0(nm) = 0.0
         alphaig0(nm) = 0.0         
         if (store_wave_direction) then
            mean_wave_direction(nm)        = 0.0
            wave_directional_spreading(nm) = 0.0  
         endif
         !
      endif   
      !
      if (store_wave_forces) then
         !
         fwx(nm)        = fwx0(nm)
         fwy(nm)        = fwy0(nm)
         dw(nm)         = dw0(nm)
         df(nm)         = df0(nm)         
         dwig(nm)       = dwig0(nm)
         dfig(nm)       = dfig0(nm)
         cg(nm)         = cg0(nm)   
         qb(nm)         = qb0(nm)         
         betamean(nm)   = beta0(nm)         
         srcig(nm)      = srcig0(nm)         
         alphaig(nm)    = alphaig0(nm)                  
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
   depth = snapwave_depth
   !
   zb = snapwave_z   
   !
   u10 = snapwave_u10
   u10dir = snapwave_u10dir   
   !   
   ! TL: we use depth now in boundary conditions for Herbers bc determination of Hm0ig, in this order we use updated values of depth through SFINCS
   !
   call update_boundary_conditions(t) ! SnapWave boundary conditions
   !
   call compute_wave_field()
   !
   snapwave_H                     = H
   snapwave_H_ig                  = H_ig
   snapwave_Tp                    = Tp
   snapwave_Tp_ig                 = Tp_ig   
   snapwave_mean_direction        = thetam
   snapwave_directional_spreading = thetam  ! TL: CORRECT? > is not spreading but mean direction?
   snapwave_Dw                    = Dw
   snapwave_Df                    = Df
   snapwave_Dwig                  = Dw_ig
   snapwave_Dfig                  = Df_ig
   snapwave_cg                    = cg
   snapwave_Qb                    = Qb
   snapwave_beta                  = beta
   snapwave_srcig                 = srcig
   snapwave_alphaig               = alphaig   
   !
   ! Convert wave force to correct unit [Dw/C] as expected by SFINCS, assumed to be piecewise (seems to work)
   snapwave_Fx                    = Fx * rho * depth
   snapwave_Fy                    = Fy * rho * depth
   !
   ! Wave periods from SnapWave, used in e.g. wavemakers - TL: moved behind call update_boundary_conditions & compute_wave_field so values at first timestep are not 0
   snapwave_tpmean = tpmean_bwv
   !
   ! Do quick check whether incoming Tpig value seems realistic, before using it:
   if (igwaves) then
       !
       snapwave_tpigmean = tpmean_bwv_ig      
       !   
       if (snapwave_tpigmean < 10.0) then
           ! These warnings should not occur here
!	       write(*,*)'DEBUG SFINCS_SnapWave - incoming tp for IG wave at wavemaker might be unrealistically small! value: ',snapwave_tpigmean
       elseif (snapwave_tpigmean > 250.0) then
!	       write(*,*)'DEBUG SFINCS_SnapWave - incoming tp for IG wave at wavemaker might be unrealistically large! value: ',snapwave_tpigmean
       endif	 
   endif
   ! TL: NOTE - in first timestep run of SnapWave tp = 0, therefore excluded that case from the check     
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
   open(500, file='sfincs.inp')   
   !
   ! Input section
   !
   call read_real_input(500,'snapwave_gamma',gamma,0.7)
   call read_real_input(500,'snapwave_alpha',alpha,1.0)
   call read_real_input(500,'snapwave_hmin',hmin,0.1)
   call read_real_input(500,'snapwave_fw',fw0,0.01)
   call read_real_input(500,'snapwave_fwig',fw0_ig,0.015)
   call read_real_input(500,'snapwave_dt',dt,36000.0)
   call read_real_input(500,'snapwave_tol',tol,1000.0)
   call read_real_input(500,'snapwave_dtheta',dtheta,10.0)
   call read_real_input(500,'snapwave_crit',crit,0.00001) !TL: Old default was 0.01
   call read_int_input(500,'snapwave_nrsweeps',nr_sweeps,4)
   call read_int_input(500,'snapwave_niter',niter, 10) !TL: Old default was 40  
   call read_int_input(500,'snapwave_baldock_opt',baldock_opt,1)     
   call read_real_input(500,'snapwave_baldock_ratio',baldock_ratio,0.2)
   call read_real_input(500,'rgh_lev_land',rghlevland,0.0)
   call read_real_input(500,'snapwave_fw_ratio',fwratio,1.0)
   call read_real_input(500,'snapwave_fwig_ratio',fwigratio,1.0)
   call read_real_input(500,'snapwave_Tpini',Tpini,1.0)
   call read_int_input (500,'snapwave_mwind',mwind,2)      
   call read_real_input(500,'sigmin',sigmin,8.0*atan(1.0)/25.0)
   call read_real_input(500,'sigmax',sigmax,8.0*atan(1.0)/1.0)   
   call read_int_input (500,'jadcgdx',jadcgdx,1)
   call read_real_input(500,'c_dispT',c_dispT,1.0)   
   call read_real_input(500,'sector',sector,180.0)
   !
   ! Settings related to IG waves:   
   call read_int_input(500,'snapwave_igwaves',igwaves_opt,1)   
   call read_real_input(500,'snapwave_alpha_ig',alpha_ig,1.0) !TODO choose whether snapwave_alphaig or snapwave_gamma_ig  
   call read_real_input(500,'snapwave_gammaig',gamma_ig,0.2)   
   call read_real_input(500,'snapwave_shinc2ig',shinc2ig,1.0)                   ! Ratio of how much of the calculated IG wave source term, is subtracted from the incident wave energy (0-1, 1=default=all energy as sink)
   call read_real_input(500,'snapwave_alphaigfac',alphaigfac,1.0)               ! Multiplication factor for IG shoaling source/sink term         
   call read_real_input(500,'snapwave_baldock_ratio_ig',baldock_ratio_ig,0.2)       
   call read_int_input(500,'snapwave_ig_opt',ig_opt,1)     
   call read_int_input(500,'snapwave_iterative_srcig',iterative_srcig,1)        ! Option whether to calculate IG source/sink term in iterative lower (better, but potentially slower, 1=default), or effectively based on previous timestep (faster, potential mismatch, =0)
   !
   ! IG boundary conditions options:
   call read_int_input(500,'snapwave_use_herbers',herbers_opt,1)    ! Choice whether you want IG Hm0&Tp be calculated by herbers (=1, default), or want to specify user defined values (0> then snapwave_eeinc2ig & snapwave_Tinc2ig are used) 
   call read_int_input(500,'snapwave_tpig_opt',tpig_opt,1) ! IG wave period option based on Herbers calculated spectrum, only used if snapwave_use_herbers = 1. Options are: 1=Tm01 (default), 2=Tpsmooth, 3=Tp, 4=Tm-1,0   
   call read_real_input(500,'snapwave_jonswapgamma',jonswapgam,3.3)  ! JONSWAP gamma value for determination offshore spectrum and IG wave conditions using Herbers, default=3.3, only used if snapwave_use_herbers = 1   
   call read_real_input(500,'snapwave_eeinc2ig',eeinc2ig,0.01)  ! Only used if snapwave_use_herbers = 0       
   call read_real_input(500,'snapwave_Tinc2ig',Tinc2ig,7.0)  ! Only used if snapwave_use_herbers = 0
   !
   ! Wind
   !
   call read_int_input(500,'snapwave_wind',wind_opt,0)   ! Flag whether to include windgrowth in SnapWave (1) or not (0, default)
   !
   ! Vegetation input
   !
   call read_int_input(500, 'vegetation', vegetation_opt, 0)
   !
   ! Input files
   call read_char_input(500,'snapwave_jonswapfile',jonswapfile,'')
   call read_char_input(500,'snapwave_bndfile',bndfile,'')
   call read_char_input(500,'snapwave_encfile',encfile,'')
   call read_char_input(500,'snapwave_bhsfile',bhsfile,'')
   call read_char_input(500,'snapwave_btpfile',btpfile,'')
   call read_char_input(500,'snapwave_bwdfile',bwdfile,'')
   call read_char_input(500,'snapwave_bdsfile',bdsfile,'') 
   call read_char_input(500,'snapwave_upwfile',upwfile,'snapwave.upw')
   call read_char_input(500,'snapwave_mskfile',mskfile,'')
   call read_char_input(500,'snapwave_depfile',depfile,'none')   
   call read_char_input(500,'snapwave_ncfile', gridfile,'snapwave_net.nc')   
   call read_char_input(500,'netsnapwavefile',netsnapwavefile,'')
   call read_char_input(500,'tref',trefstr,'20000101 000000')   ! Read again > needed in sfincs_ncinput.F90   
   !
   close(500)
   !
   igwaves          = .true.
   igherbers        = .false.    
   if (igwaves_opt==0) then
      igwaves       = .false.
      write(*,*)'SnapWave: IG waves turned OFF!'
   else
      write(*,*)'SnapWave: IG waves turned ON!'
      !
      if (herbers_opt==0) then
         write(*,*)'SnapWave: IG bc determination using Herbers turned OFF! --> Use eeinc2ig= ',eeinc2ig,' and snapwave_Tinc2ig= ',Tinc2ig
      else
         igherbers     = .true.          
         write(*,*)'SnapWave: IG bc determination using Herbers turned ON!'
      endif      
      !      
   endif
   !
   wind          = .true.
   if (wind_opt==0) then
      wind       = .false.
      write(*,*)'SnapWave: wind growth turned OFF!'
   else
      write(*,*)'SnapWave: wind growth turned ON!'
   endif   
   !
   vegetation          = .true.
   if (vegetation_opt==0) then
      vegetation       = .false.
      write(*,*)'SnapWave: vegetation turned OFF!'
   else
      write(*,*)'SnapWave: vegetation turned ON!'
   endif   
   !
   if (nr_sweeps /= 1 .and. nr_sweeps /= 4) then
      nr_sweeps = 4
      write(*,*)'SnapWave: Warning! nr_sweeps must be 1 or 4! Now set to 4.'
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
