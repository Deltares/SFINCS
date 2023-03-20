module sfincs_lib
   !
   use sfincs_spiderweb
   use sfincs_input
   use sfincs_domain
   use sfincs_structures
   use sfincs_boundaries
   use sfincs_obspoints
   use sfincs_crosssections
   use sfincs_discharges
   use sfincs_meteo
   use sfincs_data
   use sfincs_date
   use sfincs_output
   use sfincs_ncinput
   use sfincs_ncoutput
   use sfincs_momentum
   use sfincs_continuity
   use sfincs_snapwave
   use sfincs_wavemaker
   !
   implicit none
   !
   public :: sfincs_initialize
   public :: sfincs_update
   public :: sfincs_finalize
   !
   ! exposed to get start time, current time and dt
   public :: t
   public :: dt
   !
   private
   !
   ! TODO: deze ook in sfincs_data.f90 initialiseren
   integer*8  :: count0
   integer*8  :: count00
   integer*8  :: countdt0
   integer*8  :: countdt1
   integer*8  :: count1
   integer*8  :: count_rate
   integer*8  :: count_max
   integer    :: nt
   !
   integer  :: ntmapout
   integer  :: ntmaxout
   integer  :: nthisout
   !
   real*8   :: t
   real*8   :: tout
   real*4   :: dt
   real*4   :: min_dt
   real*8   :: tmapout
   real*8   :: tmaxout
   real*8   :: trstout
   real*8   :: thisout
   !
   real*4   :: maxdepth
   real*4   :: maxmaxdepth
   real*4   :: maxvmax
   real*4   :: twaveupd
   !
   logical  :: write_map   
   logical  :: write_max
   logical  :: write_his   
   logical  :: write_rst   
   !
   logical  :: update_meteo
   logical  :: update_waves
   !
   real :: tstart, tfinish, tloopflux, tloopcont, tloopstruc, tloopbnd, tloopsrc, tloopwnd1, tloopwnd2, tloopoutput, tloopsnapwave, tloopwavemaker
   real :: time_per_timestep
   real :: tinput
   real :: percdone,percdonenext,trun,trem
   !
   contains
   !
   function sfincs_initialize(config_file) result(ierr)
   !
   character(len=*) :: config_file
   integer :: ierr
   !
   ierr = 0
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   build_revision = '$Rev: v2.0.1-beta$'
   build_date     = '$Date: 2023-02-24$'
   !
   write(*,'(a)')''   
   write(*,*)'----------- Welcome to SFINCS -----------'   
   write(*,'(a)')''   
   write(*,*)' @@@@@  @@@@@@@ @@ @@  @@   @@@@   @@@@@ '
   write(*,*)'@@@ @@@ @@@@@@@ @@ @@@ @@ @@@@@@@ @@@ @@@'
   write(*,*)'@@@     @@      @@ @@@ @@ @@   @@ @@@    '
   write(*,*)' @@@@@  @@@@@@  @@ @@@@@@ @@       @@@@@ '
   write(*,*)'    @@@ @@      @@ @@ @@@ @@   @@     @@@'
   write(*,*)'@@@ @@@ @@      @@ @@  @@  @@@@@@ @@@ @@@'
   write(*,*)' @@@@@  @@      @@ @@   @   @@@@   @@@@@ '   
   write(*,'(a)')''      
   write(*,*)'             ..............              ' 
   write(*,*)'         ......:@@@@@@@@:......          '
   write(*,*)'      ..::::..@@........@@.:::::..       '
   write(*,*)'    ..:::::..@@..::..::..@@.::::::..     '
   write(*,*)'   .::::::..@@............@@.:::::::.    '
   write(*,*)'  .::::::..@@..............@@.:::::::.   '
   write(*,*)' .::::::::..@@............@@..::::::::.  '
   write(*,*)'.:::::::::...@@.@..@@..@.@@..::::::::::. '
   write(*,*)'.:::::::::...:@@@..@@..@@@:..:::::::::.. '
   write(*,*)'............@@.@@..@@..@@.@@............ '
   write(*,*)'^^^~~^^~~^^@@..............@@^^^~^^^~~^^ '
   write(*,*)'.::::::::::@@..............@@.:::::::::. '
   write(*,*)' .......:.@@.....@.....@....@@.:.......  '
   write(*,*)'  .::....@@......@.@@@.@....@@.....::.   '
   write(*,*)'   .:::~@@.:...:.@@...@@.:.:.@@~::::.    '
   write(*,*)'    .::~@@@@@@@@@@.....@@@@@@@@@~::.     '
   write(*,*)'      ..:~~~~~~~:.......:~~~~~~~:..      '
   write(*,*)'         ......................          '
   write(*,*)'             ..............              '
   write(*,'(a)')''
   write(*,*)'-----------------------------------------'
   write(*,'(a)')''
   write(*,*)'Build-Revision: ',trim(build_revision)
   write(*,*)'Build-Date:     ',trim(build_date)
   write(*,'(a)')''
   !
   call system_clock(count0, count_rate, count_max)
   !
   call read_sfincs_input()     ! Reads sfincs.inp
   !
   call read_meteo_data()       ! Reads meteo data (amu, amv, spw file etc.)
   !
   call initialize_domain()     ! Reads dep, msk, index files. Creates index, flag and depth arrays
   !
   call initialize_hydro()      ! Initializes water levels, fluxes, flags
   !
   call read_structures()       ! Reads thd files. Sets kcuv to zero where necessary
   !
   call read_boundary_data()    ! Reads bnd, bzs, etc files
   !
   call read_coastline()        ! Reads cst file
   !
   call find_boundary_indices()
   !
   call read_obs_points()       ! Reads obs file
   !
   call read_crs_file()         ! Reads cross sections
   !
   call read_discharges()       ! Reads dis and src file
   !
   if (snapwave) then
      !
      write(*,*)'Coupling SnapWave ...'
      !
      call couple_snapwave()
      !
   endif   
   !
   if (wavemaker) then
      !
      call read_wavemaker_polylines()
      !
   endif   
   !
   call system_clock(count1, count_rate, count_max)
   !
   tinput  = 1.0*(count1 - count0)/count_rate
   !
   ! Initialize some parameters
   !
   t           = t0     ! start time
   tout        = t0
   dt          = 1.0e-6 ! First time step very small
   dtavg       = 0.0    ! average time step
   maxdepth    = 999.0  ! maximum depth over time step
   maxmaxdepth = 0.0    ! maximum depth over entire simulation
   min_dt       = 0.0    ! minimum time step from compute_fluxes
   nt          = 0      ! number of time steps
   ntmapout    = 0      ! number of map time steps
   ntmaxout    = 0      ! number of max time steps
   nthisout    = 0      ! number of his time steps
   twindupd    = t0       ! time to update meteo
   twaveupd    = t0       ! time to update waves
   !
   write_map    = .false. ! write map output
   write_max    = .false. ! write max output
   write_his    = .false. ! write his output
   write_rst    = .false. ! write restart file
   !
   update_meteo = .false. ! update meteo fields
   update_waves = .false. ! update wave fields
   !
   tloopflux      = 0.0      ! flux
   tloopcont      = 0.0
   tloopstruc     = 0.0
   tloopbnd       = 0.0
   tloopsrc       = 0.0
   tloopwnd1      = 0.0
   tloopwnd2      = 0.0
   tloopsnapwave  = 0.0
   tloopwavemaker = 0.0
   !
   write(*,*)'Initializing output ...'
   !
   call initialize_output(tmapout, tmaxout, thisout, trstout)
   !
   !call acc_init( acc_device_nvidia )
   !
   end function sfincs_initialize
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   function sfincs_update(dtrange) result(ierr)
   !
   double precision, intent(in)  :: dtrange
   integer                       :: ierr
   integer                       :: nmx
   integer                       :: mmx
   integer                       :: nm
   real*4                        :: hmx
   !
   ierr = -1
   !
   !$acc data, copyin( kcs, kcuv, zs, q, q0, uv, uv0, zb, zbuv, zbuvmx, zsmax, twet, zsm, z_volume, &
   !$acc               z_flags_iref, z_flags_type, uv_flags_iref, uv_flags_type, uv_flags_vis, uv_flags_adv, uv_flags_dir, &
   !$acc               index_kcuv2, nmikcuv2, nmbkcuv2, ibkcuv2, zsb, zsb0, ibuvdir, uvmean, &
   !$acc               subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_hrep, subgrid_uv_navg, subgrid_uv_hrep_zmax, subgrid_uv_navg_zmax, &
   !$acc               subgrid_z_zmin,  subgrid_z_zmax, subgrid_z_dep, subgrid_z_volmax, &
   !$acc               z_index_uv_md1, z_index_uv_md2, z_index_uv_nd1, z_index_uv_nd2, z_index_uv_mu1, z_index_uv_mu2, z_index_uv_nu1, z_index_uv_nu2, &
   !$acc               uv_index_z_nm, uv_index_z_nmu, uv_index_u_nmd, uv_index_u_nmu, uv_index_u_ndm, uv_index_u_num, &
   !$acc               uv_index_v_ndm, uv_index_v_ndmu, uv_index_v_nm, uv_index_v_nmu, &
   !$acc               nmindsrc, qtsrc, drainage_type, drainage_params, &
   !$acc               z_index_wavemaker, wavemaker_uvmean, wavemaker_nmd, wavemaker_nmu, wavemaker_ndm, wavemaker_num, &
   !$acc               fwuv, &
   !$acc               tauwu, tauwv, tauwu0, tauwv0, tauwu1, tauwv1, &
   !$acc               windu, windv, windu0, windv0, windu1, windv1, windmax, & 
   !$acc               patm, patm0, patm1, patmb, nmindbnd, &
   !$acc               prcp, prcp0, prcp1, cumprcp, cumprcpt, netprcp, prcp, q, qinfmap, cuminf, & 
   !$acc               dxminv, dxrinv, dyrinv, dxm2inv, dxr2inv, dyr2inv, dxrinvc, dxm, dxrm, dyrm, cell_area_m2, cell_area, &
   !$acc               gn2uv, fcorio2d, min_dt ) 
   !
   ! Set target time: if dt range is negative, do not modify t1
   !
   if ( dtrange > 0.0 ) then
       t1 = t + dtrange
   endif
   !
   ! Start computational loop
   !
   write(*,'(a)')''   
   write(*,*)'---------- Starting simulation ----------'
   write(*,'(a)')''   
   !
   call system_clock(count00, count_rate, count_max)
   !
   do while (t<t1)
      !
      call system_clock(countdt0, count_rate, count_max)
      !
      write_map = .false.
      write_his = .false.
      write_max = .false.
      write_rst = .false.
      !
      ! New time step
      !
      nt = nt + 1
      dt = alfa*min_dt ! min_dt was computed in sfincs_momentum.f90 without alfa
      !
      ! A bit unclear why this happens, but large jumps in the time step lead to weird oscillations.
      ! In the 'original' sfincs v11 version, this behavior was supressed by the use of theta.
      ! Avoid this, by not not changing time step dt (used in momentum equation), but only changing dtt,
      ! which is used in the time updating and continuity equation.
      !
      ! Update time
      !
      t = t + dt
      dtavg = dtavg + dt
      !
      ! Check whether map output is required at this time step (if so, change dt)
      !
      if (t >= tmapout) then
         !
         write_map = .true.
         ntmapout  = ntmapout + 1
         tout      = tmapout 
         tmapout   = tmapout + dtmapout
         !
      endif
      !
      ! Check whether max map output is required at this time step
      !
      if (t >= tmaxout) then
         !
         write_max = .true.
         ntmaxout  = ntmaxout + 1    ! now also keep track of nr of max output
         tout      = tmaxout 
         tmaxout   = tmaxout + dtmaxout
         !
      endif
      !
      ! Check whether restart output is required at this time step
      !
      if (t >= trstout) then
         !
         write_rst = .true.
         tout      = trstout 
         !
         if (dtrstout>1.0e-6) then
            !
            ! Restart time interval in input file
            !
            trstout  = trstout + dtrstout
            !
         else
            !
            ! Restart time provided in input file
            !
            trstout = 1.0e9
            !
         endif   
         !
      endif
      !
      ! Check whether history output is required at this time step
      !
      if (t>=thisout) then
         !
         write_his = .true.
         nthisout  = nthisout + 1
         tout      = thisout 
         thisout   = thisout + dthisout
         !
      endif
      !
      if (debug .and. t>=t0out) then
         !
         ! Write every time step to map and his file when in debug mode
         !
         tout = t
         write_map = .true.
         ntmapout  = ntmapout + 1
         write_his = .true.
         nthisout  = nthisout + 1
         !
      endif
      !
      ! Check whether spatially varying meteo data needs to be updated in this time step
      !
      update_meteo = .false.
      !
      if (meteo3d) then
         if (t >= twindupd) then
            !
            meteo_t0 = twindupd
            meteo_t1 = twindupd + dtwindupd
            twindupd   = twindupd + dtwindupd ! Next time at which meteo needs to be updated
            update_meteo = .true.
            !
         endif
      endif
      !
      update_waves = .false.
      !
      if (snapwave) then
         if (t >= twaveupd) then
            twaveupd = twaveupd + dtwave ! Next time at which meteo needs to be updated
            update_waves = .true.
         endif
      endif      
      !
      ! Update meteo fields
      !
      if (wind .or. patmos .or. precip) then
         !
         if (meteo3d .and. update_meteo)  then
            !             
            ! Update spatially-varying meteo (this does not happen every time step)
            ! Read and interpolate to grid
            !
            call update_meteo_fields(t, tloopwnd1)
            !
         endif
         !
         ! Update forcing used in momentum and continuity equations (this does happen every time step)
         ! Interpolate in time
         !
         call update_meteo_forcing(t, dt, tloopwnd2)
         !
      endif
      !
      ! Update boundary conditions
      !
      call update_boundaries(t, dt, tloopbnd)
      !
      ! Update discharges
      !
      call update_discharges(t, dt, tloopsrc)
      !
      if (snapwave .and. update_waves) then
         !
         write(*,'(a,f10.1,a)')'Computing SnapWave at t = ', t, ' s'
         !
         call update_wave_field(t, tloopsnapwave)
         !
!         if (wavemaker) then
!            !
!            call update_wavemaker_points(tloopwavemaker)   
!            !
!         endif   
         !
      endif   
      !
      ! And now for the real computations !
      !
      ! First compute fluxes
      !
      call compute_fluxes(dt, min_dt, tloopflux)
      !
      if (wavemaker) then
         !
         call update_wavemaker_fluxes(t, dt, tloopwavemaker)
         !         
      endif   
      !
      if (nrstructures>0) then
         !
         call compute_fluxes_over_structures(tloopstruc)
         !
      endif
      !
      ! Update water levels
      !
      call compute_water_levels(dt, tloopcont)
      !
      ! OUTPUT
      !
      if (write_map .or. write_his .or. write_max .or. write_rst) then
         !
         if (.not. fixed_output_intervals) tout = t
         !
         call write_output(tout, write_map, write_his, write_max, write_rst, ntmapout, ntmaxout, nthisout, tloopoutput)
         !
      endif
      !
      ! Stop loop in case of instabilities (make sure water depth does not exceed stopdepth)
      !
      if (dt<dtmin .and. nt>1) then
         !
         write(*,'(a,f0.1,a)')'Maximum depth of ', stopdepth, ' m reached!!! Simulation stopped.'
         !
         ! Write map output at last time step 
         !
         call write_output(t, .true., .true., .true., .false., ntmapout + 1, ntmaxout + 1, nthisout + 1, tloopoutput)
         !
         t = t1 + 1.0
         !
      endif
      !
      percdone = min(100*(t - t0)/(t1 - t0), 100.0)
      !
      if (percdone>=percdonenext) then
         !
         percdonenext = 1.0*(int(percdone) + 5)
         call system_clock(count1, count_rate, count_max)
         trun  = 1.0*(count1 - count00)/count_rate
         trem = trun / max(0.01*percdone, 1.0e-6) - trun
         if (int(percdone)>0) then         
            write(*,'(i4,a,f7.1,a)')int(percdone),'% complete, ',trem,' s remaining ...'
         else
            write(*,'(i4,a,f7.1,a)')int(percdone),'% complete,       - s remaining ...'
         endif   
         !
      endif
      !
   enddo
   !
   !
   !$acc end data
   !
   ierr = 0
   !
   end function sfincs_update
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   function sfincs_finalize() result(ierr)
   !
   integer :: ierr
   !
   ierr = -1
   !
   call system_clock(count1, count_rate, count_max)
   !
   tstart_all = 0.0
   tfinish_all = 1.0*(count1- count00)/count_rate
   !
   call finalize_output(t,ntmaxout,tloopoutput)
   !
   if (store_tsunami_arrival_time) then
      !
      call write_tsunami_arrival_file()   
      !
   endif
   !
   dtavg = dtavg/nt
   !
   write(*,'(a)')''
   write(*,*)'---------- Simulation finished ----------'   
   write(*,'(a)')''
   write(*,'(a,f10.3)')          ' Total time             : ',tinput + tfinish_all - tstart_all
   write(*,'(a,f10.3)')          ' Total simulation time  : ',tfinish_all - tstart_all
   write(*,'(a,f10.3)')          ' Time in input          : ',tinput
   if (include_boundaries) then
      write(*,'(a,f10.3,a,f5.1,a)') ' Time in boundaries     : ',tloopbnd,' (',100*tloopbnd/(tfinish_all - tstart_all),'%)'
   endif   
   if (nsrc>0 .or. ndrn>0) then
      write(*,'(a,f10.3,a,f5.1,a)') ' Time in discharges     : ',tloopsrc,' (',100*tloopsrc/(tfinish_all - tstart_all),'%)'
   endif   
   if (meteo3d)  then
      write(*,'(a,f10.3,a,f5.1,a)') ' Time in meteo fields   : ',tloopwnd1,' (',100*tloopwnd1/(tfinish_all - tstart_all),'%)'
   endif   
   if (wind .or. patmos .or. precip) then
      write(*,'(a,f10.3,a,f5.1,a)') ' Time in meteo forcing  : ',tloopwnd2,' (',100*tloopwnd2/(tfinish_all - tstart_all),'%)'
   endif   
   write(*,'(a,f10.3,a,f5.1,a)') ' Time in momentum       : ',tloopflux,' (',100*tloopflux/(tfinish_all - tstart_all),'%)'
   if (nrstructures>0) then
      write(*,'(a,f10.3,a,f5.1,a)') ' Time in structures     : ',tloopstruc,' (',100*tloopstruc/(tfinish_all - tstart_all),'%)'
   endif   
   write(*,'(a,f10.3,a,f5.1,a)') ' Time in continuity     : ',tloopcont,' (',100*tloopcont/(tfinish_all - tstart_all),'%)'
   write(*,'(a,f10.3,a,f5.1,a)') ' Time in output         : ',tloopoutput,' (',100*tloopoutput/(tfinish_all - tstart_all),'%)'
   if (snapwave) then
      write(*,'(a,f10.3,a,f5.1,a)') ' Time in SnapWave       : ',tloopsnapwave,' (',100*tloopsnapwave/(tfinish_all - tstart_all),'%)'
   endif
   if (wavemaker) then
      write(*,'(a,f10.3,a,f5.1,a)') ' Time in wave maker     : ',tloopwavemaker,' (',100*tloopwavemaker/(tfinish_all - tstart_all),'%)'
   endif
   write(*,'(a)')''
   write(*,'(a,20f10.3)')        ' Average time step (s)  : ',dtavg
   write(*,'(a)')''
   !
   if (write_time_output) then
      open(123,file='runtimes.txt')
      write(123,'(f10.3,a)')tfinish_all - tstart_all,' % total'
      write(123,'(f10.3,a)')tinput,' % input'
      write(123,'(f10.3,a)')tloopbnd,' % boundaries'
      write(123,'(f10.3,a)')tloopsrc,' % discharges'
      write(123,'(f10.3,a)')tloopwnd1,' % meteo1'
      write(123,'(f10.3,a)')tloopwnd2,' % meteo2'
      write(123,'(f10.3,a)')tloopflux,' % momentum'
      write(123,'(f10.3,a)')tloopstruc,' % structures'
      write(123,'(f10.3,a)')tloopcont,' % continuity'
      write(123,'(f10.3,a)')tloopoutput,' % output'
      close(123)
   endif
   !
   write(*,*)'---------- Closing off SFINCS -----------'   
   !
!   call finalize_parameters()
   !
   ierr = 0
   !
   end function sfincs_finalize
   !
end module sfincs_lib
