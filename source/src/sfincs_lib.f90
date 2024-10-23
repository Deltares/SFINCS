module sfincs_lib
   !
   use omp_lib
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
   use sfincs_infiltration
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
   error = 0 ! Error code. This is now only set to 1 in case of instabilities. Could also use other error codes, e.g. for missing files.
   !
   ierr = 0 ! Always 0 or 1  ! Always 0 or 1 
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   build_revision = '$Rev: v2.1.1-Dollerup - double precision coordinates'
   build_date     = '$Date: 2024-10-22'
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
   call initialize_domain()     ! Reads dep, msk, index files, creates index, flag and depth arrays, initializes hydro quantities
   !
   call read_structures()       ! Reads thd files and sets kcuv to zero where necessary
   !
   call read_boundary_data()    ! Reads bnd, bzs, etc files
   !
   ! call read_coastline()        ! Reads cst file. Do we still do this ?
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
      write(*,*)'Coupling with SnapWave ...'
      !
      call couple_snapwave(crsgeo)
      !
   endif   
   !
   if (wavemaker) then
      !
      call read_wavemaker_polylines()
      !
   endif
   !
   call set_advection_mask()
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
   min_dt      = 0.0    ! minimum time step from compute_fluxes
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
   tloopflux      = 0.0
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
   ! Quadtree no longer needed, so deallocate (this is done in sfincs_domain.f90)
   ! 
   call deallocate_quadtree()
   !
   !call acc_init( acc_device_nvidia )
   !
   ierr = error
   !
   end function sfincs_initialize
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   function sfincs_update(dtrange) result(ierr)
   !
   double precision, intent(in)  :: dtrange
   integer                       :: ierr
   !
   ierr = 0
   !
   ! Copy arrays to GPU memory
   ! 
   !$acc data, copyin( kcs, kfuv, kcuv, zs, zs0, zsderv, q, q0, uv, uv0, zb, zbuv, zbuvmx, zsmax, maxzsm, qmax, vmax, twet, zsm, z_volume, &
   !$acc               z_flags_iref, uv_flags_iref, uv_flags_type, uv_flags_dir, mask_adv, &
   !$acc               index_kcuv2, nmikcuv2, nmbkcuv2, ibkcuv2, zsb, zsb0, ibuvdir, uvmean, &
   !$acc               subgrid_uv_zmin, subgrid_uv_zmax, subgrid_uv_havg, subgrid_uv_nrep, subgrid_uv_pwet, &
   !$acc               subgrid_uv_havg_zmax, subgrid_uv_nrep_zmax, subgrid_uv_fnfit, subgrid_uv_navg_w, &
   !$acc               subgrid_z_zmin,  subgrid_z_zmax, subgrid_z_dep, subgrid_z_volmax, &
   !$acc               z_index_uv_md, z_index_uv_nd, z_index_uv_mu, z_index_uv_nu, &
   !$acc               uv_index_z_nm, uv_index_z_nmu, uv_index_u_nmd, uv_index_u_nmu, uv_index_u_ndm, uv_index_u_num, &
   !$acc               uv_index_v_ndm, uv_index_v_ndmu, uv_index_v_nm, uv_index_v_nmu, &
   !$acc               nmindsrc, qtsrc, drainage_type, drainage_params, &
   !$acc               z_index_wavemaker, wavemaker_uvmean, wavemaker_nmd, wavemaker_nmu, wavemaker_ndm, wavemaker_num, &
   !$acc               fwuv, &
   !$acc               tauwu, tauwv, tauwu0, tauwv0, tauwu1, tauwv1, &
   !$acc               windu, windv, windu0, windv0, windu1, windv1, windmax, & 
   !$acc               patm, patm0, patm1, patmb, nmindbnd, &
   !$acc               prcp, prcp0, prcp1, cumprcp, cumprcpt, netprcp, prcp, qinfmap, cuminf, & 
   !$acc               dxminv, dxrinv, dyrinv, dxm2inv, dxr2inv, dyr2inv, dxrinvc, dxm, dxrm, dyrm, cell_area_m2, cell_area, &
   !$acc               gn2uv, fcorio2d, min_dt, storage_volume, nuvisc, &
   !$acc               cuv_index_uv, cuv_index_uv1, cuv_index_uv2 )
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
   write(*,'(a,i0,a,i0,a)')' ---- Using ', omp_get_max_threads(), ' of ', omp_get_num_procs(), ' available threads -----'   
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
      dt = alfa * min_dt ! min_dt was computed in sfincs_momentum.f90 without alfa
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
         tout      = max(tmapout, t - dt) 
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
         tout      = max(tmaxout, t - dt) 
         !
         if (t < t1) then 
            tmaxout   = tmaxout + dtmaxout       
            ! in case the last 'dt' made us exactly past tstop time 't1', 
            ! then we don't want to flag later another dtmax output timestep in 'finalize_output' check,
            ! so if t > t1 don't add 'dtmaxout' again            
         endif         
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
         tout      = max(thisout, t - dt) 
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
         !
         call update_meteo_forcing(t, dt, tloopwnd2)
         !
         ! Update infiltration
         !
         if (infiltration) then
             !
             ! Compute infiltration rates
             !
             call update_infiltration_map(dt)
             !
         endif
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
         call timer(t3)          
         !
         call update_wave_field(t, tloopsnapwave)
         !
         call timer(t4)                   
         write(*,'(a,f10.1,a,f6.2,a)')'Computing SnapWave at t = ', t, ' s took ', t4 - t3, ' seconds'         
         !
         ! Maybe we'll add moving wave makers back at some point
         !
         ! if (wavemaker) then
         !    !
         !    call update_wavemaker_points(tloopwavemaker)   
         !    !
         ! endif   
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
         ! if (.not. fixed_output_intervals) tout = t
         !
         call write_output(tout, write_map, write_his, write_max, write_rst, ntmapout, ntmaxout, nthisout, tloopoutput)
         !
      endif
      !      
      ! Stop loop in case of instabilities (make sure water depth does not exceed stopdepth)
      !
      if (dt<dtmin .and. nt>1) then
         !
         error = 1
         write(error_message,'(a,f0.1,a)')'Error! Maximum depth of ', stopdepth, ' m reached!!! Simulation stopped.'
         !
         ! Write map output at last time step 
         !
         ntmaxout = ntmaxout + 1 ! Max sure that max output is not called again through 'finalize_output' 
         !
         call write_output(t, .true., .true., .true., .false., ntmapout + 1, ntmaxout, nthisout + 1, tloopoutput)
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
   end function sfincs_update
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   function sfincs_finalize() result(ierr)
   !
   integer :: ierr
   !
   call system_clock(count1, count_rate, count_max)
   !
   tstart_all = 0.0
   tfinish_all = 1.0*(count1- count00)/count_rate
   !
   call finalize_output(t,ntmaxout,tloopoutput,tmaxout)
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
   ! call finalize_parameters()
   !
   ierr = 0
   !
   if (error > 0) then
      !
      ! Stop depth was exceeded, or e.g. input file missing
      !
      write(*,*)''
      write(*,*)trim(error_message)
      ierr = error
      !
   endif
   !
   end function sfincs_finalize
   !
end module sfincs_lib
