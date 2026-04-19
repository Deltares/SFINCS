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
   use sfincs_runup_gauges
   use sfincs_discharges
   use sfincs_src_structures
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
   use sfincs_nonhydrostatic
   use sfincs_bathtub
   use sfincs_openacc
   use sfincs_log
   use sfincs_timers
   use sfincs_timestep_analysis
   use sfincs_screendump
   !
   implicit none
   !
   public :: sfincs_initialize
   public :: sfincs_update
   public :: sfincs_finalize
   !
   ! exposed to get start time, current time and dt
   !
   public :: t
   public :: dt
   !
   private
   !
   integer    :: nt
   !
   integer  :: ntmapout
   integer  :: ntmaxout
   integer  :: nthisout
   !
   real*8   :: t
   real*8   :: tout
   real*4   :: dt
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
   real :: time_per_timestep
   !
   contains
   !
   function sfincs_initialize() result(ierr)
   !
   integer :: ierr
   !
   call open_log()   
   !
   error = 0 ! Error code. This is now only set to 1 in case of instabilities. Could also use other error codes, e.g. for missing files.
   !
   ierr = 0 ! Always 0 or 1  ! Always 0 or 1 
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   build_revision = "$Rev: v2.3.2 Mount Faber + branch-redo-infiltration + urban_drainage + discharges + timers + screendump"
   build_date     = "$Date: 2026-04-19"
   !
   call screendump_startup()
   !
   call timer_start('Input')
   !
   call write_log('------ Preparing model simulation --------', 1)
   call write_log('', 1)
   !
   call write_log('Reading input file ...', 0) 
   call read_sfincs_input()     ! Reads sfincs.inp
   !
   if (.not. bathtub) then
      !
      call write_log('Reading meteo data ...', 0)
      !
      call read_meteo_data()       ! Reads meteo data (amu, amv, spw file etc.)
      !
   endif   
   !
   call write_log('Preparing domain ...', 0)
   !
   call initialize_domain()     ! Reads dep, msk, index files, creates index, flag and depth arrays, initializes hydro quantities
   !
   call read_structures()       ! Reads thd files and sets kcuv to zero where necessary
   !
   call read_boundary_data()    ! Reads bnd, bzs, etc files
   !
   call find_boundary_indices()
   !
   call read_obs_points()       ! Reads obs file
   !
   call read_crs_file()         ! Reads cross sections
   !
   call read_rug_file()         ! Read runup gauge file
   !
   call initialize_discharges()       ! Reads dis and src file (river point discharges)
   !
   call initialize_src_structures()   ! Reads drn file (pumps / culverts / check valves / gates)
   !
   if (nonhydrostatic) then
      !
      ! Initialize non-hydrostatic solver
      !
      call write_log('Initialize non-hydrostatic solver ...', 0) 
      !
      call initialize_nonhydrostatic()
      !
   endif   
   !
   if (wavemaker) then
      !
      call initialize_wavemakers()
      !
   endif
   !
   if (bathtub) then
      !
      call write_log('Initialize bathtub mode ...', 0) 
      !
      call initialize_bathtub()
      !
   endif   
   !
   call screendump_processes()
   !
   if (snapwave) then
      !
      call write_log('Coupling with SnapWave ...', 1)
      call couple_snapwave(crsgeo)
      !
   endif   
   !
   call timer_stop('Input')
   !
   ! Initialize some parameters
   !
   t           = t0     ! start time
   tout        = t0
   dt          = 1.0e-6 ! First time step very small
   min_dt      = 1.0e-6 ! First time step very small
   dtavg       = 0.0    ! average time step
   maxdepth    = 999.0  ! maximum depth over time step
   maxmaxdepth = 0.0    ! maximum depth over entire simulation
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
   call write_log('Initializing output ...', 0)
   !
   call initialize_output(tmapout, tmaxout, thisout, trstout)
   !
   ! Quadtree no longer needed, so deallocate (this is done in sfincs_domain.f90)
   ! 
   call deallocate_quadtree()
   !
   call initialize_openacc() ! Enter data region
   !
   ierr = error
   !
   call write_log('', 1)
   call write_log('---------- Starting simulation -----------', 1)
   write(logstr,'(a,i0,a,i0,a)')'---- Using ', omp_get_max_threads(), ' of ', omp_get_num_procs(), ' available threads ----'   
   call write_log(logstr, 1)
   call write_log('', 1)
   !
   call timer_start('Simulation loop')
   !
   end function sfincs_initialize
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   function sfincs_update(dtrange) result(ierr)
   !
   double precision, intent(in)  :: dtrange
   integer                       :: ierr
   real*8                        :: tend !< end of update interval
   real*4                        :: dtchk !< dt to check for instability
   logical                       :: single_time_step
   !
   ierr = 0
   !
   ! Set target time: if dt range is negative, do not modify t1
   !
   single_time_step = .false.
   !   
   if (dtrange < -999.0) then
      !
      ! Regular (used in executable, called by sfincs.f90)
      !
      tend = t1
      !      
   elseif (dtrange <= 0.0) then
      ! 
      ! XMI run for one time step (called from sfincs_bmi.f90 -> update)
      !
      ! Set tend to t1
      ! 
      tend = t1
      !
      ! But make sure only one iteration is done
      !
      single_time_step = .true.
      !
   else !  ( dtrange > 0.0 ) 
      ! 
      ! XMI run for time interval dtrange (called from sfincs_bmi.f90 -> update_until)
      !
      tend = t + dtrange
      !
   endif
   !
   ! Start computational loop
   !
   do while (t < tend)
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
      dtchk = alfa * min_dt
      !
      if (bathtub) then
         !
         ! In bathtub mode, use fixed time step
         !
         if (nt > 1) then
            !
            dt = bathtub_dt
            dtchk = bathtub_dt
            !
         endif   
         !
      endif
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
            !
            tmaxout   = tmaxout + dtmaxout
            !
            ! in case the last 'dt' made us exactly past tstop time 't1', 
            ! then we don't want to flag later another dtmax output timestep in 'finalize_output' check,
            ! so if t > t1 don't add 'dtmaxout' again
            !
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
      if (t >= thisout) then
         !
         write_his = .true.
         nthisout  = nthisout + 1
         tout      = max(thisout, t - dt) 
         thisout   = thisout + dthisout
         !
      endif
      !
      if (debug .and. t >= t0out) then
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
            call update_meteo_fields(t)
            !
         endif
         !
         ! Update forcing used in momentum and continuity equations (this does happen every time step)
         !
         call update_meteo_forcing(t, dt)
         !
      endif
      !
      ! Update boundary conditions
      !
      call update_boundaries(t, dt)
      !
      ! Update discharges (river sources) and src-point structures (pumps/gates/...)
      !
      call update_discharges(t, dt)
      call update_src_structures(t, dt)
      !
      if (snapwave .and. update_waves) then
         !
         ! Update wave fields from SnapWave coupling (this happens at intervals of dtwave)
         !
         call update_wave_field(t)
         !
      endif
      !
      if (bathtub) then
         !
         ! In bathtub mode, only update water levels based on boundary conditions
         !
         call bathtub_compute_water_levels()
         !
      else
         !
         ! And now for the real computations ! Unless we're in bathtub mode
         !
         ! First compute fluxes
         !
         call compute_fluxes(dt)
         !
         if (timestep_analysis) then
             !
             call timestep_analysis_update(min_dt)
             !
         endif           
         !
         if (wavemaker) then
            !
            call update_wavemaker_fluxes(t, dt)
            !         
         endif   
         !
         if (nrstructures>0) then
            !
            call compute_fluxes_over_structures()
            !
         endif
         !      
         if (nonhydrostatic) then
            !
            if (t < nh_tstop) then ! Check if non-hydrostatic corrections still need to be made
               !
               ! Apply non-hydrostatic pressure corrections to q and uv
               !
               call compute_nonhydrostatic(dt)
               !
            endif   
            !
         endif
         !      
         ! Update water levels
         !
         ! Update continuity (discharges, infiltration, drainage, water levels)
         !
         call update_continuity(t, dt)
         !
      endif   
      !
      ! OUTPUT
      !      
      if (write_map .or. write_his .or. write_max .or. write_rst) then
         !
         ! if (.not. fixed_output_intervals) tout = t
         !
         call write_output(tout, write_map, write_his, write_max, write_rst, ntmapout, ntmaxout, nthisout)
         !
      endif
      !      
      ! Stop loop in case of instabilities (make sure time step 'dtmin' does not get too small compared to 'uvmax' flow velocity)
      !
      if (dtchk < dtmin .and. nt > 1) then
         !
         error = 1
         !
         write(error_message,'(a,f0.4,a,f0.1,a)')'Error! Minimum time step of ', dtmin, ' s reached ! Current velocity exceeded uvmax ', uvmax, ' m/s. Simulation stopped.'
         !
         ! Write map output at last time step 
         !
         ntmaxout = ntmaxout + 1 ! Max sure that max output is not called again through 'finalize_output' 
         !
         call write_output(t, .true., .true., .true., .false., ntmapout + 1, ntmaxout, nthisout + 1)
         !
         t = t1 + 1.0
         !
      endif
      !
      call screendump_progress(t, t0, t1)
      !
      if (single_time_step) then
         !
         ! Update was called with XMI update so only run one time step.
         ! Do this by setting tend to t0, so next iteration does not run.
         !
         tend = t0
         !
      endif      
      !
   enddo
   !
   end function sfincs_update
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   function sfincs_finalize() result(ierr)
   !
   integer :: ierr
   !
   call timer_stop('Simulation loop')
   !
   if (timestep_analysis) then
      !
      call timestep_analysis_finalize(nt)
      !
   endif
   !
   call finalize_output(t, ntmaxout, tmaxout)
   !
   call finalize_openacc() ! Exit data region
   !
   dtavg = dtavg / (nt - 1)
   !
   call screendump_run_finished(dtavg)
   !
   if (timestep_analysis) then
      !
      call timestep_analysis_write_log()
      !
   endif
   !
   if (write_time_output) then
      !
      call timer_write_runtimes_file(123, 'runtimes.txt')
      !
   endif
   !
   call write_log('----------- Closing off SFINCS -----------', 1)
   !
   ! call finalize_parameters()
   !
   ierr = 0
   !
   if (error > 0) then
      !
      ! Stop depth was exceeded, or e.g. input file missing
      !
      call write_log('', 1)
      call write_log(trim(error_message), 1)
      ierr = error
      !
   endif
   !
   call close_log()
   !
   end function sfincs_finalize
   !
end module sfincs_lib
