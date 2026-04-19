module sfincs_input
   !
   ! Parser for the SFINCS main input file `sfincs.inp` plus the small set
   ! of primitive helpers that read one keyword at a time from that file.
   !
   ! `sfincs.inp` is a flat keyword / value text file (one `key = value`
   ! pair per line, comment lines start with `#`, `!`, or `@`). A read is
   ! performed by rewinding the file, scanning until the matching key is
   ! found, and extracting the value string. When the key is absent, the
   ! caller's supplied default is returned.
   !
   ! The module does not own the variables it fills — it writes directly
   ! into module-level state declared in sfincs_data, sfincs_src_structures,
   ! sfincs_discharges, etc.
   !
   ! Subroutines:
   !
   !   read_sfincs_input()
   !     Main driver. Opens sfincs.inp, calls the per-type helpers below
   !     once per keyword, then derives secondary flags (e.g. crsgeo vs.
   !     Coriolis, subgrid vs. regular, bathtub overrides). Called once
   !     from sfincs_initialize (sfincs_lib).
   !
   !   read_real_input(fileid, keyword, value, default)
   !     Read one real*4 keyword. Called from read_sfincs_input.
   !
   !   read_real_array_input(fileid, keyword, value, default, nr)
   !     Read one space-separated real*4 array keyword. Called from
   !     read_sfincs_input.
   !
   !   read_int_input(fileid, keyword, value, default)
   !     Read one integer keyword. Called from read_sfincs_input.
   !
   !   read_char_input(fileid, keyword, value, default)
   !     Read one character-string keyword. Called from read_sfincs_input.
   !
   !   read_logical_input(fileid, keyword, value, default)
   !     Read one logical keyword. Accepts `1/0`, `y/n`, `t/f` (upper or
   !     lower) as true/false. Called from read_sfincs_input.
   !
   !   read_line(line0, keystr, valstr)
   !     Strip tab/line-ending noise, split `key = value` on the first `=`,
   !     strip any trailing `# ...` inline comment. Called from each of
   !     the read_*_input helpers.
   !
   !   notabs(instr, outstr, ilen)
   !     Expand embedded tab characters into spaces preserving 8-column
   !     tab stops. Called from read_line.
   !

contains
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine read_sfincs_input()
      !
      ! Top-level driver: open sfincs.inp, pull every keyword the solver
      ! knows about into the appropriate module-level variable, then
      ! compute derived flags (CRS / Coriolis, subgrid vs. regular,
      ! wavemaker modes, bathtub overrides, advection scheme).
      !
      ! Called from: sfincs_initialize (sfincs_lib).
      !
      use sfincs_data
      use sfincs_date
      use sfincs_log
      use sfincs_error
      use sfincs_src_structures, only: drnfile
      use sfincs_discharges,     only: srcfile, disfile, netsrcdisfile
      !
      implicit none
      !
      integer*8                                 :: dtsec
      logical                                   :: ok
      character(len=256)                        :: wmsigstr
      character(len=256)                        :: advstr
      character(len=256)                        :: removed_input
      !
      ok = check_file_exists('sfincs.inp', 'SFINCS input file', .true.)
      !
      open(500, file='sfincs.inp')
      !
      ! Grid geometry and time window
      !
      call read_int_input(500,     'mmax',                            mmax,                            0)                 ! number of grid cells in m-direction
      call read_int_input(500,     'nmax',                            nmax,                            0)                 ! number of grid cells in n-direction
      call read_real_input(500,    'dx',                              dx,                              0.0)               ! cell size in m-direction (m)
      call read_real_input(500,    'dy',                              dy,                              0.0)               ! cell size in n-direction (m)
      call read_real_input(500,    'x0',                              x0,                              0.0)               ! grid origin x (m or deg)
      call read_real_input(500,    'y0',                              y0,                              0.0)               ! grid origin y (m or deg)
      call read_real_input(500,    'rotation',                        rotation,                        0.0)               ! grid rotation (deg, counter-clockwise from east)
      call read_char_input(500,    'tref',                            trefstr,                         'none')            ! reference time (yyyymmdd HHMMSS); defaults to tstart
      call read_char_input(500,    'tstart',                          tstartstr,                       '20000101 000000') ! simulation start time
      call read_char_input(500,    'tstop',                           tstopstr,                        '20000101 000000') ! simulation stop time
      call read_real_input(500,    'tspinup',                         tspinup,                         0.0)               ! spin-up interval after t0 (s)
      call read_real_input(500,    't0out',                           t0out,                           -999.0)            ! output start time (s rel. tref); -999 = t0
      call read_real_input(500,    't1out',                           t1out,                           -999.0)            ! output stop time  (s rel. tref); -999 = t1
      call read_real_input(500,    'dtout',                           dtmapout,                        0.0)               ! map output interval (s); 0 = no map output
      call read_real_input(500,    'dtmaxout',                        dtmaxout,                        9999999.0)         ! zsmax etc. interval (s); 0 = end-of-run only
      call read_real_input(500,    'dtrstout',                        dtrstout,                        0.0)               ! restart interval (s); 0 = no periodic restart
      call read_real_input(500,    'trstout',                         trst,                            -999.0)            ! single restart time (s rel. tref); -999 = unused
      call read_real_input(500,    'dthisout',                        dthisout,                        600.0)             ! his output interval (s)
      call read_real_input(500,    'dtwave',                          dtwave,                          3600.0)            ! SnapWave update interval (s)
      call read_real_input(500,    'dtwnd',                           dtwindupd,                       1800.0)            ! 2D meteo update interval (s)
      !
      ! Solver and physical constants
      !
      call read_real_input(500,    'alpha',                           alfa,                            0.50)              ! CFL Courant factor
      call read_real_input(500,    'theta',                           theta,                           1.0)               ! semi-implicit theta; <1 adds smoothing
      call read_real_input(500,    'hmin_cfl',                        hmin_cfl,                        0.1)               ! minimum depth used in CFL check (m)
      call read_real_input(500,    'manning',                         manning,                         0.04)              ! uniform Manning n (s/m^(1/3))
      call read_real_input(500,    'manning_land',                    manning_land,                    -999.0)            ! Manning n above rghlevland (s/m^(1/3))
      call read_real_input(500,    'manning_sea',                     manning_sea,                     -999.0)            ! Manning n below rghlevland (s/m^(1/3))
      call read_real_input(500,    'rgh_lev_land',                    rghlevland,                      0.0)               ! bed level separating land/sea friction (m)
      call read_real_input(500,    'zsini',                           zini,                            0.0)               ! initial water level (m)
      call read_real_input(500,    'qinf',                            qinf,                            0.0)               ! uniform infiltration rate (mm/hr); converted below
      call read_real_input(500,    'dtmax',                           dtmax,                           60.0)              ! upper bound on computational dt (s)
      call read_real_input(500,    'huthresh',                        huthresh,                        0.05)              ! wet/dry depth threshold (m)
      call read_real_input(500,    'huvmin',                          huvmin,                          0.0)               ! minimum depth for uv = q / max(hu, huvmin) (output + advection)
      call read_real_input(500,    'rhoa',                            rhoa,                            1.25)              ! air density (kg/m3)
      call read_real_input(500,    'rhow',                            rhow,                            1024.0)            ! water density (kg/m3)
      call read_char_input(500,    'inputformat',                     inputtype,                       'bin')             ! legacy bin/asc toggle for binary inputs
      call read_char_input(500,    'outputformat',                    outputtype,                      'net')             ! global output format (bin/asc/net)
      call read_char_input(500,    'outputtype_map',                  outputtype_map,                  'nil')             ! map-file output format (nil = follow outputformat)
      call read_char_input(500,    'outputtype_his',                  outputtype_his,                  'nil')             ! his-file output format (nil = follow outputformat)
      call read_int_input(500,     'nc_deflate_level',                nc_deflate_level,                2)                 ! netCDF deflate level (0-9)
      call read_int_input(500,     'bndtype',                         bndtype,                         1)                 ! boundary condition type
      call read_logical_input(500, 'advection',                       advection,                       .true.)            ! enable momentum advection terms
      call read_real_input(500,    'latitude',                        latitude,                        0.0)               ! reference latitude for projected Coriolis (deg)
      call read_real_input(500,    'pavbnd',                          pavbnd,                          0.0)               ! atmospheric pressure applied at boundary (Pa)
      call read_real_input(500,    'gapres',                          gapres,                          101200.0)          ! atmospheric reference pressure (Pa)
      call read_int_input(500,     'baro',                            baro,                            1)                 ! include atmospheric-pressure gradient (1=on, 0=off)
      call read_char_input(500,    'utmzone',                         utmzone,                         'nil')             ! UTM zone string (e.g. '17N')
      call read_int_input(500,     'epsg',                            epsg,                            0)                 ! EPSG integer code for the grid
      call read_char_input(500,    'epsg',                            epsg_code,                       'nil')             ! EPSG as string (fallback)
      call read_real_input(500,    'advlim',                          advlim,                          1.0)               ! cap on advection term
      call read_real_input(500,    'slopelim',                        slopelim,                        9999.9)            ! cap on bed-slope water-level gradient
      call read_real_input(500,    'qinf_zmin',                       qinf_zmin,                       0.0)               ! minimum bed level for infiltration to apply (m)
      call read_real_input(500,    'btfilter',                        btfilter,                        60.0)              ! bathtub filter time scale (s)
      call read_real_input(500,    'sfacinf',                         sfacinf,                         0.2)               ! SCS initial-abstraction fraction (0.2S)
      call read_logical_input(500, 'radstr',                          radstr,                          .false.)           ! radiation-stress forcing from SnapWave
      call read_logical_input(500, 'crsgeo',                          crsgeo,                          .false.)           ! interpret grid coords as geographic (WGS84)
      call read_logical_input(500, 'coriolis',                        coriolis,                        .true.)            ! include Coriolis force
      call read_logical_input(500, 'amprblock',                       ampr_block,                      .true.)            ! treat 2D rainfall as block (true) or linearly interpolated (false)
      call read_real_input(500,    'spwmergefrac',                    spw_merge_frac,                  0.5)               ! merge factor for spiderweb wind composite
      call read_logical_input(500, 'usespwprecip',                    use_spw_precip,                  .true.)            ! use precipitation field from spiderweb file
      call read_logical_input(500, 'global',                          global,                          .false.)           ! treat grid as global (wrap in x)
      call read_real_input(500,    'nuvisc',                          nuviscdim,                       0.01)              ! viscosity coefficient (m2/s)
      call read_logical_input(500, 'viscosity',                       viscosity,                       .false.)           ! enable horizontal viscosity term
      call read_logical_input(500, 'spinup_meteo',                    spinup_meteo,                    .false.)           ! ramp wind/pressure from zero during tspinup
      call read_real_input(500,    'waveage',                         waveage,                         -999.0)            ! wave age (for SnapWave wind growth)
      call read_logical_input(500, 'snapwave',                        snapwave,                        .false.)           ! enable coupled SnapWave wave solver
      call read_logical_input(500, 'dtoutfixed',                      fixed_output_intervals,          .true.)            ! snap map/his to exact intervals (true) or let them drift with dt (false)
      !
      ! Wave maker parameters. Old 3-letter keywords (wvmfile, wfpfile,
      ! whifile, wtifile, wstfile, wmtfilter, wmfred, wmsignal, wmhmin,
      ! nfreqsinc/ig, freq*min/max*inc/ig) are retained for backward
      ! compatibility; the wavemaker_<name> keywords below override them
      ! when supplied.
      !
      call read_char_input(500,    'wavemaker_wvmfile',               wavemaker_wvmfile,               'none')            ! wavemaker polyline file
      if (wavemaker_wvmfile(1:4) == 'none') &
         call read_char_input(500, 'wvmfile',                         wavemaker_wvmfile,               'none')            ! legacy keyword
      call read_char_input(500,    'wavemaker_wfpfile',               wavemaker_wfpfile,               'none')            ! wavemaker forcing points file
      if (wavemaker_wfpfile(1:4) == 'none') &
         call read_char_input(500, 'wfpfile',                         wavemaker_wfpfile,               'none')            ! legacy keyword
      call read_char_input(500,    'wavemaker_whifile',               wavemaker_whifile,               'none')            ! wavemaker wave-height time series file
      if (wavemaker_whifile(1:4) == 'none') &
         call read_char_input(500, 'whifile',                         wavemaker_whifile,               'none')            ! legacy keyword
      call read_char_input(500,    'wavemaker_wtifile',               wavemaker_wtifile,               'none')            ! wavemaker wave-period time series file
      if (wavemaker_wtifile(1:4) == 'none') &
         call read_char_input(500, 'wtifile',                         wavemaker_wtifile,               'none')            ! legacy keyword
      call read_char_input(500,    'wavemaker_wstfile',               wavemaker_wstfile,               'none')            ! wavemaker wave set-up time series file
      if (wavemaker_wstfile(1:4) == 'none') &
         call read_char_input(500, 'wstfile',                         wavemaker_wstfile,               'none')            ! legacy keyword
      !
      call read_real_input(500,    'wmtfilter',                       wavemaker_filter_time,           600.0)             ! wavemaker filter time scale (s, legacy keyword)
      call read_real_input(500,    'wavemaker_filter_time',           wavemaker_filter_time,           wavemaker_filter_time) ! override with new keyword if present
      call read_real_input(500,    'wmfred',                          wavemaker_filter_fred,           0.99)              ! wavemaker filter fred (legacy keyword)
      call read_real_input(500,    'wavemaker_filter_fred',           wavemaker_filter_fred,           wavemaker_filter_fred) ! override with new keyword if present
      call read_char_input(500,    'wmsignal',                        wmsigstr,                        'spectrum')        ! wavemaker signal type (legacy keyword)
      call read_char_input(500,    'wavemaker_signal',                wmsigstr,                        trim(wmsigstr))    ! override with new keyword if present
      call read_real_input(500,    'wmhmin',                          wavemaker_hmin,                  0.1)               ! wavemaker minimum depth for wave generation (legacy keyword)
      call read_real_input(500,    'wavemaker_hmin',                  wavemaker_hmin,                  wavemaker_hmin)    ! override with new keyword if present
      call read_int_input(500,     'nfreqsinc',                       wavemaker_nfreqs_inc,            100)               ! wavemaker number of incident-wave frequencies (legacy)
      call read_int_input(500,     'wavemaker_nfreqs_inc',            wavemaker_nfreqs_inc,            wavemaker_nfreqs_inc) ! override
      call read_real_input(500,    'freqmininc',                      wavemaker_freqmin_inc,           0.04)              ! wavemaker incident-wave min frequency (Hz, legacy)
      call read_real_input(500,    'wavemaker_freqmin_inc',           wavemaker_freqmin_inc,           wavemaker_freqmin_inc) ! override
      call read_real_input(500,    'freqmaxinc',                      wavemaker_freqmax_inc,           1.0)               ! wavemaker incident-wave max frequency (Hz, legacy)
      call read_real_input(500,    'wavemaker_freqmax_inc',           wavemaker_freqmax_inc,           wavemaker_freqmax_inc) ! override
      call read_int_input(500,     'nfreqsig',                        wavemaker_nfreqs_ig,             100)               ! wavemaker number of IG-wave frequencies (legacy)
      call read_int_input(500,     'wavemaker_nfreqs_ig',             wavemaker_nfreqs_ig,             wavemaker_nfreqs_ig) ! override
      call read_real_input(500,    'freqminig',                       wavemaker_freqmin_ig,            0.0)               ! wavemaker IG-wave min frequency (Hz, legacy)
      call read_real_input(500,    'wavemaker_freqmin_ig',            wavemaker_freqmin_ig,            wavemaker_freqmin_ig) ! override
      call read_real_input(500,    'freqmaxig',                       wavemaker_freqmax_ig,            0.1)               ! wavemaker IG-wave max frequency (Hz, legacy)
      call read_real_input(500,    'wavemaker_freqmax_ig',            wavemaker_freqmax_ig,            wavemaker_freqmax_ig) ! override
      !
      call read_real_input(500,    'wavemaker_tinc2ig',               wavemaker_tinc2ig,               -1.0)              ! wavemaker ig/inc period ratio (<=0 uses Herbers)
      call read_real_input(500,    'wavemaker_surfslope',             wavemaker_surfslope,             -1.0)              ! wavemaker surf-zone slope for empirical Tp_ig (van Ormondt et al., 2021)
      call read_real_input(500,    'wavemaker_hm0_ig_factor',         wavemaker_hm0_ig_factor,         1.0)               ! wavemaker IG Hm0 scaling factor
      call read_real_input(500,    'wavemaker_hm0_inc_factor',        wavemaker_hm0_inc_factor,        1.0)               ! wavemaker incident Hm0 scaling factor
      call read_real_input(500,    'wavemaker_gammax',                wavemaker_gammax,                1.0)               ! wavemaker maximum Hrms/h
      call read_real_input(500,    'wavemaker_tpmin',                 wavemaker_tpmin,                 1.0)               ! wavemaker minimum Tp (s)
      call read_logical_input(500, 'wavemaker_hig',                   wavemaker_hig,                   .true.)            ! wavemaker include IG waves
      call read_logical_input(500, 'wavemaker_hinc',                  wavemaker_hinc,                  .false.)           ! wavemaker include incident waves
      !
      ! Numerical parameters
      !
      call read_char_input(500,    'advection_scheme',                advstr,                          'upw1')            ! advection scheme label ('upw1' = 1st-order upwind, 'original' = legacy)
      call read_real_input(500,    'btrelax',                         btrelax,                         3600.0)            ! bathtub relaxation time (s)
      call read_logical_input(500, 'wiggle_suppression',              wiggle_suppression,              .true.)            ! suppress spurious free-surface oscillations
      call read_real_input(500,    'structure_relax',                 structure_relax,                 4.0)               ! drainage-structure state-machine smoothing steps
      call read_real_input(500,    'wiggle_factor',                   wiggle_factor,                   0.1)               ! wiggle-suppression amplitude factor
      call read_real_input(500,    'wiggle_threshold',                wiggle_threshold,                0.1)               ! wiggle-suppression trigger threshold
      call read_real_input(500,    'uvlim',                           uvlim,                           10.0)              ! clipping velocity for momentum (m/s)
      call read_real_input(500,    'uvmax',                           uvmax,                           1000.0)            ! error-trigger velocity for momentum (m/s)
      call read_logical_input(500, 'friction2d',                      friction2d,                      .true.)            ! apply friction at every UV point (true) or cell-wise (false)
      call read_logical_input(500, 'advection_mask',                  advection_mask,                  .true.)            ! mask advection near dry cells
      call read_real_input(500,    'nuviscfac',                       nuviscfac,                       100.0)             ! multiplier on nuvisc near "difficult" points
      call read_logical_input(500, 'nonh',                            nonhydrostatic,                  .false.)           ! enable non-hydrostatic pressure corrector
      call read_real_input(500,    'nh_fnudge',                       nh_fnudge,                       0.9)               ! non-hydrostatic nudging factor
      call read_real_input(500,    'nh_tstop',                        nh_tstop,                        -999.0)            ! non-hydrostatic stop time (s rel. tref); -999 = t1+999
      call read_real_input(500,    'nh_tol',                          nh_tol,                          0.001)             ! non-hydrostatic solver tolerance
      call read_int_input(500,     'nh_itermax',                      nh_itermax,                      100)               ! non-hydrostatic solver max iterations
      call read_logical_input(500, 'h73table',                        h73table,                        .false.)           ! tabulate h^(7/3) for friction
      call read_real_input(500,    'rugdepth',                        runup_gauge_depth,               0.05)              ! runup gauge trigger depth (m)
      call read_logical_input(500, 'wave_enhanced_roughness',         wave_enhanced_roughness,         .false.)           ! augment bed roughness with wave orbital velocity
      call read_logical_input(500, 'use_bcafile',                     use_bcafile,                     .true.)            ! use tidal components from bca file
      call read_real_input(500,    'factor_wind',                     factor_wind,                     1.0)               ! scaling factor on wind forcing
      call read_real_input(500,    'factor_pres',                     factor_pres,                     1.0)               ! scaling factor on atmospheric pressure
      call read_real_input(500,    'factor_prcp',                     factor_prcp,                     1.0)               ! scaling factor on precipitation
      call read_real_input(500,    'factor_spw_size',                 factor_spw_size,                 1.0)               ! scaling factor on spiderweb radius
      call read_logical_input(500, 'bathtub',                         bathtub,                         .false.)           ! run in bathtub (no momentum) mode
      call read_real_input(500,    'bathtub_fachs',                   bathtub_fac_hs,                  0.2)               ! bathtub Hs multiplier
      call read_real_input(500,    'bathtub_dt',                      bathtub_dt,                      -999.0)            ! bathtub time step (s); -999 = use dtmapout
      !
      ! Domain files
      !
      call read_char_input(500,    'qtrfile',                         qtrfile,                         'none')            ! quadtree netCDF file
      call read_char_input(500,    'depfile',                         depfile,                         'none')            ! bed-level (depth) file
      call read_char_input(500,    'inifile',                         zsinifile,                       'none')            ! initial water-level file
      call read_char_input(500,    'rstfile',                         rstfile,                         'none')            ! restart input file
      call read_char_input(500,    'mskfile',                         mskfile,                         'none')            ! active-cell mask file
      call read_char_input(500,    'indexfile',                       indexfile,                       'none')            ! index-to-active-cell mapping file
      call read_char_input(500,    'cstfile',                         cstfile,                         'none')            ! coastline polyline file
      call read_char_input(500,    'sbgfile',                         sbgfile,                         'none')            ! subgrid tables netCDF file
      call read_char_input(500,    'thdfile',                         thdfile,                         'none')            ! thin dams polyline file
      call read_char_input(500,    'weirfile',                        weirfile,                        'none')            ! weirs polyline file
      call read_char_input(500,    'manningfile',                     manningfile,                     'none')            ! spatially-varying Manning n file
      call read_char_input(500,    'drnfile',                         drnfile,                         'none')            ! drainage-structures (pumps/gates/culverts) TOML file
      call read_char_input(500,    'urbfile',                         urbfile,                         'none')            ! urban drainage zones TOML file
      call read_char_input(500,    'volfile',                         volfile,                         'none')            ! depression-storage volume file
      !
      ! Forcing files (ascii / binary)
      !
      call read_char_input(500,    'bndfile',                         bndfile,                         'none')            ! water-level boundary points
      call read_char_input(500,    'bzsfile',                         bzsfile,                         'none')            ! water-level boundary time series
      call read_char_input(500,    'bcafile',                         bcafile,                         'none')            ! tidal components per boundary point
      call read_char_input(500,    'bzifile',                         bzifile,                         'none')            ! IG wave boundary time series
      call read_char_input(500,    'bdrfile',                         bdrfile,                         'none')            ! downstream river boundary file
      call read_char_input(500,    'srcfile',                         srcfile,                         'none')            ! river-point source locations
      call read_char_input(500,    'disfile',                         disfile,                         'none')            ! river-point discharge time series
      call read_char_input(500,    'spwfile',                         spwfile,                         'none')            ! spiderweb tropical-cyclone file
      call read_char_input(500,    'wndfile',                         wndfile,                         'none')            ! uniform wind time series
      call read_char_input(500,    'prcfile',                         prcpfile,                        'none')            ! uniform precipitation time series
      if (prcpfile(1:4) == 'none') then
         !
         call read_char_input(500, 'precipfile',                      prcpfile,                        'none')            ! legacy keyword for prcfile
         !
      endif
      call read_char_input(500,    'amufile',                         amufile,                         'none')            ! 2D wind u-component file
      call read_char_input(500,    'amvfile',                         amvfile,                         'none')            ! 2D wind v-component file
      call read_char_input(500,    'ampfile',                         ampfile,                         'none')            ! 2D atmospheric pressure file
      call read_char_input(500,    'amprfile',                        amprfile,                        'none')            ! 2D precipitation rate file
      call read_char_input(500,    'z0lfile',                         z0lfile,                         'none')            ! 2D land roughness (z0) file
      !
      ! NetCDF-format forcing files (FEWS-style)
      !
      call read_char_input(500,    'netbndbzsbzifile',                netbndbzsbzifile,                'none')            ! combined bnd/bzs/bzi netCDF file
      call read_char_input(500,    'netsrcdisfile',                   netsrcdisfile,                   'none')            ! combined src/dis netCDF file
      call read_char_input(500,    'netamuamvfile',                   netamuamvfile,                   'none')            ! combined amu/amv netCDF file
      call read_char_input(500,    'netamprfile',                     netamprfile,                     'none')            ! 2D precipitation netCDF file
      call read_char_input(500,    'netampfile',                      netampfile,                      'none')            ! 2D pressure netCDF file
      call read_char_input(500,    'netspwfile',                      netspwfile,                      'none')            ! netCDF spiderweb file
      !
      ! Infiltration and losses
      !
      call read_char_input(500,    'infiltrationfile',                infiltrationfile,                'none')            ! infiltration parameters TOML file
      call read_char_input(500,    'infiltrationtype',                inftype,                         'none')            ! infiltration flavor (con, c2d, cna, cnb, gai, hor, bkt)
      call read_char_input(500,    'bucketfile',                      removed_input,                   '__removed_keyword_not_present__')
      if (trim(removed_input) /= '__removed_keyword_not_present__') then
         !
         write(logstr,'(a)') 'Error    : keyword bucketfile has been removed. Use infiltrationfile together with infiltrationtype = bkt.'
         call stop_sfincs(trim(logstr), 1)
         !
      endif
      call read_char_input(500,    'bucket_loss_frac',                removed_input,                   '__removed_keyword_not_present__')
      if (trim(removed_input) /= '__removed_keyword_not_present__') then
         !
         write(logstr,'(a)') 'Error    : keyword bucket_loss_frac has been removed. Add bucket_loss to infiltrationfile instead.'
         call stop_sfincs(trim(logstr), 1)
         !
      endif
      !
      ! Legacy binary infiltration inputs (kept for backward compatibility).
      !
      call read_char_input(500,    'qinffile',                        qinffile,                        'none')            ! binary spatially-varying infiltration field
      call read_char_input(500,    'scsfile',                         scsfile,                         'none')            ! SCS curve-number S field (legacy binary)
      call read_char_input(500,    'smaxfile',                        smaxfile,                        'none')            ! SCS max storage S field (legacy binary)
      call read_char_input(500,    'sefffile',                        sefffile,                        'none')            ! SCS effective storage S_e field (legacy binary)
      call read_char_input(500,    'psifile',                         psifile,                         'none')            ! Green-Ampt suction head (legacy binary, mm)
      call read_char_input(500,    'sigmafile',                       sigmafile,                       'none')            ! Green-Ampt maximum moisture deficit (legacy binary)
      call read_char_input(500,    'ksfile',                          ksfile,                          'none')            ! Green-Ampt saturated hydraulic conductivity (legacy binary, mm/hr)
      call read_char_input(500,    'f0file',                          f0file,                          'none')            ! Horton initial infiltration capacity F0 (legacy binary)
      call read_char_input(500,    'fcfile',                          fcfile,                          'none')            ! Horton asymptotic infiltration rate Fc (legacy binary)
      call read_char_input(500,    'kdfile',                          kdfile,                          'none')            ! Horton decay constant k (legacy binary, 1/hr)
      call read_real_input(500,    'horton_kr_kd',                    horton_kr_kd,                    10.0)              ! Horton recovery/decay ratio (recovery is kr_kd times slower than decay)
      !
      ! Output files
      !
      call read_char_input(500,    'obsfile',                         obsfile,                         'none')            ! observation-point locations file
      call read_char_input(500,    'crsfile',                         crsfile,                         'none')            ! cross-section polyline file
      call read_char_input(500,    'rugfile',                         rugfile,                         'none')            ! runup-gauge locations file
      call read_logical_input(500, 'storevelmax',                     store_maximum_velocity,          .false.)           ! store maximum flow velocity on dtmaxout interval (only if dtmaxout > 0)
      call read_logical_input(500, 'storefluxmax',                    store_maximum_flux,              .false.)           ! store maximum flux on dtmaxout interval (only if dtmaxout > 0)
      call read_logical_input(500, 'storevel',                        store_velocity,                  .false.)           ! store velocity on dtout interval
      call read_logical_input(500, 'storecumprcp',                    store_cumulative_precipitation,  .false.)           ! store cumulative precipitation + infiltration on dtmaxout interval
      call read_logical_input(500, 'storetwet',                       store_twet,                      .false.)           ! store per-cell wet duration
      call read_logical_input(500, 'storetzsmax',                     store_t_zsmax,                   .false.)           ! store time stamp of zsmax occurrence
      call read_logical_input(500, 'storehsubgrid',                   store_hsubgrid,                  .false.)           ! store hmax in subgrid mode (zsmax - subgrid_z_zmin)
      call read_logical_input(500, 'storehmean',                      store_hmean,                     .false.)           ! store hmax as subgrid-mean depth instead of max (requires storehsubgrid)
      call read_real_input(500,    'twet_threshold',                  twet_threshold,                  0.01)              ! water-depth threshold counting a cell as wet (storetwet)
      call read_logical_input(500, 'store_tsunami_arrival_time',      store_tsunami_arrival_time,      .false.)           ! store tsunami arrival time per cell
      call read_real_input(500,    'tsunami_arrival_threshold',       tsunami_arrival_threshold,       0.01)              ! water-depth threshold for tsunami arrival
      call read_logical_input(500, 'timestep_analysis',               timestep_analysis,               .false.)           ! write per-cell timestep limiter diagnostics
      call read_logical_input(500, 'storeqdrain',                     store_qdrain,                    .true.)            ! store per-drainage-structure discharge in his file
      call read_logical_input(500, 'store_river_discharge',           store_river_discharge,           .false.)           ! store per-river-point discharge in his file
      call read_logical_input(500, 'store_urban_drainage_discharge',  store_urban_drainage_discharge,  .false.)           ! store per-urban-zone outfall discharge in his file
      call read_logical_input(500, 'store_cumulative_urban_drainage', store_cumulative_urban_drainage, .false.)           ! store cumulative urban drainage depth per cell in map file
      call read_logical_input(500, 'storezvolume',                    store_zvolume,                   .false.)           ! store subgrid cell volume (requires subgrid)
      call read_logical_input(500, 'storestoragevolume',              store_storagevolume,             .false.)           ! store remaining storage volume (requires subgrid + volfile)
      call read_logical_input(500, 'writeruntime',                    write_time_output,               .false.)           ! write runtimes.txt at end of simulation
      call read_logical_input(500, 'debug',                           debug,                           .false.)           ! debug output at every time step
      call read_logical_input(500, 'storemeteo',                      store_meteo,                     .false.)           ! store 2D meteo forcing fields in map file
      call read_logical_input(500, 'storemaxwind',                    store_wind_max,                  .false.)           ! store maximum wind speed (requires storemeteo)
      call read_logical_input(500, 'storefw',                         store_wave_forces,               .false.)           ! store wave-radiation forces
      call read_logical_input(500, 'storewavdir',                     store_wave_direction,            .false.)           ! store wave direction
      call read_logical_input(500, 'regular_output_on_mesh',          use_quadtree_output,             .false.)           ! write quadtree output on regular m/n mesh
      call read_logical_input(500, 'store_dynamic_bed_level',         store_dynamic_bed_level,         .false.)           ! store time-varying bed level (subgrid)
      call read_logical_input(500, 'snapwave_use_nearest',            snapwave_use_nearest,            .true.)            ! use nearest-neighbour lookup for SnapWave boundary points
      call read_int_input(500,     'percentage_done',                 percdoneval,                     5)                 ! progress-reporter interval (% complete)
      !
      ! Limit progress reporter to (0, 100]%
      !
      percdoneval = max(min(percdoneval, 100), 0)
      !
      ! Coupled SnapWave solver parameters
      !
      call read_logical_input(500, 'snapwave_wind',                   snapwavewind,                    .false.)           ! feed wind into SnapWave (implies storing wind speed/direction)
      call read_real_input(500,    'snapwave_waveforces_factor',      waveforces_factor,               1.0)               ! multiplier on SnapWave wave forces
      !
      ! Wind drag coefficient table
      !
      call read_int_input(500,     'cdnrb',                           cd_nr,                           0)                 ! number of wind-drag breakpoints (0 = use defaults)
      !
      if (cd_nr == 0) then
         !
         ! Standard Smith & Banke-style Cd curve: constant at low wind,
         ! linear rise, then plateau at high wind.
         !
         cd_nr = 3
         !
         allocate(cd_wnd(cd_nr))
         allocate(cd_val(cd_nr))
         !
         cd_wnd(1) =   0.0
         cd_wnd(2) =  28.0
         cd_wnd(3) =  50.0
         cd_val(1) = 0.0010
         cd_val(2) = 0.0025
         cd_val(3) = 0.0025
         !
      else
         !
         call read_real_array_input(500, 'cdwnd', cd_wnd, 0.0, cd_nr)
         call read_real_array_input(500, 'cdval', cd_val, 0.0, cd_nr)
         !
      endif
      !
      ! Late retry of dtmapout for older sfincs.inp files that put it after
      ! keywords which could have shifted its position.
      !
      if (dtmapout == 0.0) then
         !
         call read_real_input(500, 'dtmapout', dtmapout, 0.0)
         !
      endif
      !
      close(500)
      !
      if (epsg == 0) then
         !
         call write_log('Warning : no EPSG code defined', 0)
         !
      endif
      !
      ! If tref not provided, assume tref = tstart.
      !
      if (trefstr(1:4) == 'none') then
         !
         trefstr = tstartstr
         !
         write(logstr, *) 'Warning : no tref provided, set to tstart: ', trefstr
         call write_log(logstr, 1)
         !
      endif
      !
      ! Compute simulation time span in seconds, relative to tref.
      !
      call time_difference(trefstr, tstartstr, dtsec)
      t0 = dtsec * 1.0
      call time_difference(trefstr, tstopstr, dtsec)
      t1 = dtsec * 1.0
      tspinup = t0 + tspinup
      !
      g   = 9.81
      pi  = 3.14159
      gn2 = 9.81 * 0.02 * 0.02                     ! only used in subgrid mode
      !
      qinf = qinf / (3600 * 1000)                  ! mm/hr -> m/s
      !
      rotation = rotation * pi / 180
      cosrot   = cos(rotation)
      sinrot   = sin(rotation)
      !
      area  = dx * dy
      dxy   = min(dx, dy)
      dxinv = 1.0 / dx
      dyinv = 1.0 / dy
      !
      manning2d = .false.
      if (manningfile /= 'none') manning2d = .true.
      !
      ! CRS and Coriolis parameter
      !
      fcorio = 0.0
      !
      if (.not. crsgeo) then
         !
         ! Projected: compute fcorio from latitude; zero it out at the
         ! equator or if the user did not set a latitude at all.
         !
         fcorio = 2 * 7.2921e-05 * sin(latitude * pi / 180)
         !
         if (latitude < 0.01 .and. latitude > -0.01) coriolis = .false.
         !
      else
         !
         ! Geographic: fcorio2d is filled in later by sfincs_domain.
         !
      endif
      !
      if (crsgeo) then
         call write_log('Info    : input grid interpreted as geographic coordinates', 0)
      else
         call write_log('Info    : input grid interpreted as projected coordinates', 0)
      endif
      !
      if (coriolis) then
         call write_log('Info    : turning on Coriolis', 0)
      else
         call write_log('Info    : turning off Coriolis', 0)
      endif
      !
      if (.not. crsgeo .and. .not. coriolis) then
         call write_log('Info    : no Coriolis, as latitude is not specified in sfincs.inp', 0)
      endif
      !
      ! Map/his output window: default to tstart..tstop.
      !
      if (t0out < -900.0) t0out = t0
      t0out = max(t0out, t0)
      if (t1out < -900.0) t1out = t1
      !
      store_maximum_waterlevel = .false.
      if (dtmaxout > 0.0) store_maximum_waterlevel = .true.
      !
      ! Apply gates to the flags now that the full set of inputs has been read.
      !
      if (dtmaxout <= 0.0) then
         store_maximum_velocity = .false.
         store_maximum_flux     = .false.
      endif
      !
      ! storemeteo implies store_wind (SFINCS needs 2D wind to feed the
      ! meteo map output); storemaxwind is only meaningful if we are
      ! storing the wind in the first place.
      !
      if (store_meteo)      store_wind     = .true.
      if (.not. store_wind) store_wind_max = .false.
      !
      ! SnapWave with wind-growth needs wind stored and a dedicated
      ! snapwavewind flag; snapwavewind is ignored when SnapWave is off.
      !
      if (.not. snapwave) snapwavewind = .false.
      if (snapwavewind)   store_wind   = .true.
      !
      if (viscosity) call write_log('Info    : turning on process: Viscosity', 0)
      !
      ! Map/his format fallback: inherit the global outputformat when either
      ! per-file format was left at 'nil'.
      !
      if ((outputtype_map == 'nil') .or. (outputtype_his == 'nil')) then
         outputtype_map = outputtype
         outputtype_his = outputtype
      endif
      !
      if (sbgfile(1:4) /= 'none') then
         subgrid = .true.
         call write_log('Info    : running SFINCS with subgrid bathymetry', 0)
      else
         subgrid = .false.
         call write_log('Info    : running SFINCS with regular bathymetry', 0)
      endif
      !
      if (subgrid .eqv. .true. .and. store_hsubgrid .eqv. .true. .and. store_hmean .eqv. .false.) then
         !
         call write_log('Info    : storing maximum depth in subgrid cell for hmax output', 0)
         !
      elseif (subgrid .eqv. .true. .and. store_hsubgrid .eqv. .true. .and. store_hmean .eqv. .true.) then
         !
         call write_log('Info    : storing mean depth in subgrid cell for hmax output', 0)
         !
      endif
      !
      ! store_zvolume / store_storagevolume are subgrid-only.
      !
      if (.not. subgrid) store_zvolume = .false.
      !
      thetasmoothing = .false.
      if (theta < 0.9999) thetasmoothing = .true.      ! use 0.9999 instead of 1.0 for numerical robustness
      !
      wavemaker          = .false.
      wavemaker_spectrum = .true.
      !
      if (wavemaker_wvmfile(1:4) /= 'none') then
         !
         wavemaker = .true.
         !
         call write_log('Info    : turning on process: Dynamic waves', 0)
         !
         if (wmsigstr(1:3) == 'mon') then
            !
            wavemaker_spectrum = .false.
            !
            call write_log('Info    : use monochromatic wave spectrum', 0)
            !
         endif
         !
      endif
      !
      use_storage_volume = .false.
      !
      if (volfile(1:4) /= 'none') then
         if (subgrid) then
            use_storage_volume = .true.
         else
            call write_log('Warning : storage volume only supported for subgrid topographies!', 1)
            store_storagevolume = .false.
         endif
      else
         store_storagevolume = .false.
      endif
      !
      if (advection) then
         !
         ! Default scheme is 1st-order upwind; 'original' keeps the legacy form.
         !
         advection_scheme = 1
         !
         call write_log('Info    : turning on advection', 0)
         !
         if (trim(advstr) == 'original') then
            advection_scheme = 0
            call write_log('Info    : advection scheme : Original', 0)
         elseif (trim(advstr) == 'upw1') then
            advection_scheme = 1
            call write_log('Info    : advection scheme : first-order upwind', 0)
         else
            write(logstr, *) 'Warning : advection scheme ', trim(advstr), ' not recognized! Using default upw1 instead!'
            call write_log(logstr, 1)
         endif
         !
      endif
      !
      if (nonhydrostatic) then
         !
         if (nh_tstop > 0.0) then
            nh_tstop = t0 + nh_tstop
         else
            nh_tstop = t1 + 999.0
         endif
         !
      endif
      !
      if (bathtub) then
         !
         call write_log('Info    : turning on process: Bathtub flooding', 0)
         !
         ! Time step defaults to dtmapout when the user does not set it.
         !
         if (bathtub_dt < 0.0) bathtub_dt = dtmapout
         !
         dthisout = bathtub_dt
         !
         ! Turn off processes not needed for bathtub flooding. Forcing the
         ! input file paths to 'none' makes each initialize_* routine take
         ! its standard early-return path; that way the counters
         ! (nr_discharge_points, nr_src_structures, nr_urban_drainage_zones)
         ! and derived logicals (discharges, drainage_structures,
         ! urban_drainage) stay consistent with the "no input" state.
         !
         srcfile       = 'none'
         disfile       = 'none'
         netsrcdisfile = 'none'
         drnfile       = 'none'
         urbfile       = 'none'
         !
         meteo3d        = .false.
         wind           = .false.
         store_meteo    = .false.
         store_wind     = .false.
         store_wind_max = .false.
         precip         = .false.
         patmos         = .false.
         if (snapwave) bathtub_snapwave = .true.
         snapwave               = .false.
         infiltration           = .false.
         store_velocity         = .false.
         store_maximum_velocity = .false.
         !
      endif
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine read_real_input(fileid, keyword, value, default)
      !
      ! Read a single real*4 keyword from an already-open sfincs.inp. The
      ! file is rewound on each call; scanning is linear. If the keyword
      ! is not found, `value` is set to `default`.
      !
      ! Called from: read_sfincs_input (this module).
      !
      implicit none
      !
      character(*), intent(in)  :: keyword
      integer,      intent(in)  :: fileid
      real*4,       intent(out) :: value
      real*4,       intent(in)  :: default
      !
      character(len=256) :: keystr
      character(len=256) :: valstr
      character(len=256) :: line
      integer            :: stat
      !
      value = default
      !
      rewind(fileid)
      !
      do while (.true.)
         !
         read(fileid, '(a)', iostat=stat) line
         if (stat == -1) exit
         !
         call read_line(line, keystr, valstr)
         !
         if (trim(keystr) == trim(keyword)) then
            !
            read(valstr, *) value
            exit
            !
         endif
         !
      enddo
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine read_real_array_input(fileid, keyword, value, default, nr)
      !
      ! Read one whitespace-separated real*4 array keyword. Allocates
      ! `value(nr)` on the way in and fills it from the matching line; if
      ! the keyword is absent, every slot is left at `default`.
      !
      ! Called from: read_sfincs_input (this module).
      !
      implicit none
      !
      character(*), intent(in)                        :: keyword
      integer,      intent(in)                        :: fileid
      integer,      intent(in)                        :: nr
      real*4,       intent(in)                        :: default
      real*4, dimension(:), intent(out), allocatable  :: value
      !
      character(len=256) :: keystr
      character(len=256) :: valstr
      character(len=256) :: line
      integer            :: m, stat
      !
      allocate(value(nr))
      !
      value = default
      !
      rewind(fileid)
      !
      do while (.true.)
         !
         read(fileid, '(a)', iostat=stat) line
         if (stat == -1) exit
         !
         call read_line(line, keystr, valstr)
         !
         if (trim(keystr) == trim(keyword)) then
            !
            read(valstr, *) (value(m), m = 1, nr)
            exit
            !
         endif
         !
      enddo
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine read_int_input(fileid, keyword, value, default)
      !
      ! Read a single integer keyword. Same scanning contract as
      ! read_real_input.
      !
      ! Called from: read_sfincs_input (this module).
      !
      implicit none
      !
      character(*), intent(in)  :: keyword
      integer,      intent(in)  :: fileid
      integer,      intent(out) :: value
      integer,      intent(in)  :: default
      !
      character(len=256) :: keystr
      character(len=256) :: valstr
      character(len=256) :: line
      integer            :: stat
      !
      value = default
      !
      rewind(fileid)
      !
      do while (.true.)
         !
         read(fileid, '(a)', iostat=stat) line
         if (stat == -1) exit
         !
         call read_line(line, keystr, valstr)
         !
         if (trim(keystr) == trim(keyword)) then
            !
            read(valstr, *) value
            exit
            !
         endif
         !
      enddo
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine read_char_input(fileid, keyword, value, default)
      !
      ! Read a single character-string keyword. The entire right-hand side
      ! (after stripping any trailing `# ...` comment) becomes `value`.
      !
      ! Called from: read_sfincs_input (this module).
      !
      implicit none
      !
      character(*), intent(in)  :: keyword
      integer,      intent(in)  :: fileid
      character(*), intent(in)  :: default
      character(*), intent(out) :: value
      !
      character(len=256) :: keystr
      character(len=256) :: valstr
      character(len=256) :: line
      integer            :: stat
      !
      value = default
      !
      rewind(fileid)
      !
      do while (.true.)
         !
         read(fileid, '(a)', iostat=stat) line
         if (stat == -1) exit
         !
         call read_line(line, keystr, valstr)
         !
         if (trim(keystr) == trim(keyword)) then
            !
            value = valstr
            exit
            !
         endif
         !
      enddo
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine read_logical_input(fileid, keyword, value, default)
      !
      ! Read a single logical keyword. Accepts `1`, `y`, `Y`, `t`, `T` as
      ! true; anything else (including absence plus fallback to `default`,
      ! `0`, `n`, `N`, `f`, `F`) as false.
      !
      ! Called from: read_sfincs_input (this module).
      !
      implicit none
      !
      character(*), intent(in)  :: keyword
      integer,      intent(in)  :: fileid
      logical,      intent(in)  :: default
      logical,      intent(out) :: value
      !
      character(len=256) :: keystr
      character(len=256) :: valstr
      character(len=256) :: line
      integer            :: stat
      !
      value = default
      !
      rewind(fileid)
      !
      do while (.true.)
         !
         read(fileid, '(a)', iostat=stat) line
         if (stat == -1) exit
         !
         call read_line(line, keystr, valstr)
         !
         if (trim(keystr) == trim(keyword)) then
            !
            if (valstr(1:1) == '1' .or. valstr(1:1) == 'y' .or. valstr(1:1) == 'Y' .or. &
                valstr(1:1) == 't' .or. valstr(1:1) == 'T') then
               value = .true.
            else
               value = .false.
            endif
            !
            exit
            !
         endif
         !
      enddo
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine read_line(line0, keystr, valstr)
      !
      ! Split one `key = value` line into key and value substrings.
      ! Strips leading/trailing whitespace, any tab characters (replaced
      ! by spaces via notabs), and a trailing `#`-delimited inline
      ! comment. Blank lines and lines starting with `#`, `!`, or `@`
      ! return empty strings.
      !
      ! Called from: read_real_input / read_real_array_input /
      ! read_int_input / read_char_input / read_logical_input.
      !
      implicit none
      !
      character(*), intent(in)  :: line0
      character(*), intent(out) :: keystr
      character(*), intent(out) :: valstr
      !
      character(len=256) :: line
      integer            :: j, ilen, jn
      !
      keystr = ''
      valstr = ''
      !
      ! Expand tabs to spaces in-place.
      !
      call notabs(line0, line, ilen)
      !
      ! Remove Windows-style `\r` line ending if present.
      !
      jn = index(line, '\r')
      if (jn > 0) line = line(1:jn - 1)
      !
      line = trim(line)
      !
      if (line(1:1) == '#' .or. line(1:1) == '!' .or. line(1:1) == '@') return
      !
      j = index(line, '=')
      if (j == 0) return
      !
      keystr = trim(line(1:j - 1))
      valstr = trim(line(j + 1:))
      !
      ! Strip inline comment after `#`.
      !
      jn = index(valstr, '#')
      if (jn > 0) valstr = trim(valstr(1:jn - 1))
      !
      valstr = adjustl(trim(valstr))
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine notabs(instr, outstr, ilen)
      !
      ! Expand embedded tab characters into spaces while keeping columns
      ! aligned (tab stops every 8 characters). Lets downstream tokenizers
      ! treat `key<tab>=<tab>value` and `key = value` identically.
      !
      ! Author: John S. Urban. See also GNU/Unix commands expand(1) /
      ! unexpand(1).
      !
      ! Called from: read_line (this module).
      !
      use iso_fortran_env, only : error_unit
      !
      implicit none
      !
      character(len=*), intent(in)  :: instr     ! input line (may contain tab characters)
      character(len=*), intent(out) :: outstr    ! tab-expanded output
      integer,          intent(out) :: ilen      ! column position of last character written
      !
      integer, parameter :: tabsize = 8          ! tab stops every 8th column
      character(len=1)   :: c
      integer            :: ipos                 ! position in outstr for next character
      integer            :: lenin                ! length of instr (trailing blanks trimmed)
      integer            :: lenout               ! capacity of outstr
      integer            :: i10                  ! cursor through instr
      !
      ipos   = 1
      lenin  = len(instr)
      lenin  = len_trim(instr(1:lenin))
      lenout = len(outstr)
      outstr = ' '
      !
      do i10 = 1, lenin
         !
         c = instr(i10:i10)
         !
         if (ichar(c) == 9) then
            !
            ! Tab character: advance ipos to the next tab stop.
            !
            ipos = ipos + (tabsize - (mod(ipos - 1, tabsize)))
            !
         else
            !
            if (ipos > lenout) then
               write(error_unit, *) '*notabs* output string overflow'
               exit
            else
               outstr(ipos:ipos) = c
               ipos              = ipos + 1
            endif
            !
         endif
         !
      enddo
      !
      ilen = len_trim(outstr(:ipos))
      !
   end subroutine
   !
end module
