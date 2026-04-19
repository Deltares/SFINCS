module sfincs_input
   !
   ! Parser for the SFINCS main input file `sfincs.inp` plus a small set
   ! of primitive helpers that read one keyword at a time from that file.
   !
   ! `sfincs.inp` is a flat keyword / value text file (one `key = value`
   ! pair per line, comment lines start with `#`, `!`, or `@`). Reads go
   ! through the generic `get_keyword(...)` interface below. Every call
   ! accepts an optional array of deprecated keyword aliases — the new
   ! key is tried first, then each legacy alias in order, then the
   ! caller's supplied default. A one-line warning is written to the log
   ! the first time an alias is matched, so users see that their .inp is
   ! using a deprecated name.
   !
   ! The module does not own the variables it fills — it writes directly
   ! into module-level state declared in sfincs_data, sfincs_src_structures,
   ! sfincs_discharges, etc.
   !
   ! Subroutines:
   !
   !   read_sfincs_input()
   !     Main driver. Opens sfincs.inp, pulls every keyword it knows
   !     about via get_keyword(), then derives secondary flags (CRS /
   !     Coriolis, subgrid vs. regular, bathtub overrides). Called once
   !     from sfincs_initialize (sfincs_lib).
   !
   !   get_keyword(fileid, keyword, value, default [, legacy])
   !     Generic interface; resolves by the type of `value`. Type-
   !     specific module procedures do the actual scan:
   !        get_keyword_real    (real*4 scalar)
   !        get_keyword_int     (integer scalar)
   !        get_keyword_char    (character*(*) scalar)
   !        get_keyword_logical (logical scalar)
   !     Called from read_sfincs_input.
   !
   !   get_keyword_real_array(fileid, keyword, value, default, nr [, legacy])
   !     Variant for whitespace-separated real arrays (sized nr).
   !     Called from read_sfincs_input.
   !
   !   find_value(fileid, keyword, valstr, found)
   !     Scan the file once for a single keyword; return the raw value
   !     string and whether it was found. Called from each
   !     get_keyword_* procedure.
   !
   !   warn_legacy(legacy_key, new_key)
   !     Emit a deprecation warning to the log. Called from each
   !     get_keyword_* procedure when a legacy alias has been resolved.
   !
   !   read_line(line0, keystr, valstr)
   !     Strip tab/line-ending noise, split `key = value` on the first
   !     `=`, strip any trailing `# ...` inline comment. Called from
   !     find_value.
   !
   !   notabs(instr, outstr, ilen)
   !     Expand embedded tab characters into spaces preserving 8-column
   !     tab stops. Called from read_line.
   !
   use sfincs_log, only: write_log, logstr
   !
   implicit none
   !
   private
   public :: read_sfincs_input
   public :: get_keyword, get_keyword_real_array
   !
   interface get_keyword
      module procedure get_keyword_real
      module procedure get_keyword_int
      module procedure get_keyword_char
      module procedure get_keyword_logical
   end interface
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
      use sfincs_error
      use sfincs_src_structures, only: drnfile
      use sfincs_discharges,     only: srcfile, disfile, netsrcdisfile
      !
      implicit none
      !
      integer*8          :: dtsec
      logical            :: ok
      character(len=256) :: wmsigstr
      character(len=256) :: advstr
      character(len=256) :: removed_input
      !
      ok = check_file_exists('sfincs.inp', 'SFINCS input file', .true.)
      !
      open(500, file='sfincs.inp')
      !
      ! Grid geometry and time window
      !
      call get_keyword(500, 'mmax',                            mmax,                            0)                 ! number of grid cells in m-direction
      call get_keyword(500, 'nmax',                            nmax,                            0)                 ! number of grid cells in n-direction
      call get_keyword(500, 'dx',                              dx,                              0.0)               ! cell size in m-direction (m)
      call get_keyword(500, 'dy',                              dy,                              0.0)               ! cell size in n-direction (m)
      call get_keyword(500, 'x0',                              x0,                              0.0)               ! grid origin x (m or deg)
      call get_keyword(500, 'y0',                              y0,                              0.0)               ! grid origin y (m or deg)
      call get_keyword(500, 'rotation',                        rotation,                        0.0)               ! grid rotation (deg, counter-clockwise from east)
      call get_keyword(500, 'tref',                            trefstr,                         'none')            ! reference time (yyyymmdd HHMMSS); defaults to tstart
      call get_keyword(500, 'tstart',                          tstartstr,                       '20000101 000000') ! simulation start time
      call get_keyword(500, 'tstop',                           tstopstr,                        '20000101 000000') ! simulation stop time
      call get_keyword(500, 'tspinup',                         tspinup,                         0.0)               ! spin-up interval after t0 (s)
      call get_keyword(500, 't0out',                           t0out,                           -999.0)            ! output start time (s rel. tref); -999 = t0
      call get_keyword(500, 't1out',                           t1out,                           -999.0)            ! output stop time  (s rel. tref); -999 = t1
      call get_keyword(500, 'dtout',                           dtmapout,                        0.0)               ! map output interval (s); 0 = no map output
      call get_keyword(500, 'dtmaxout',                        dtmaxout,                        9999999.0)         ! zsmax etc. interval (s); 0 = end-of-run only
      call get_keyword(500, 'dtrstout',                        dtrstout,                        0.0)               ! restart interval (s); 0 = no periodic restart
      call get_keyword(500, 'trstout',                         trst,                            -999.0)            ! single restart time (s rel. tref); -999 = unused
      call get_keyword(500, 'dthisout',                        dthisout,                        600.0)             ! his output interval (s)
      call get_keyword(500, 'dtwave',                          dtwave,                          3600.0)            ! SnapWave update interval (s)
      call get_keyword(500, 'dtwnd',                           dtwindupd,                       1800.0)            ! 2D meteo update interval (s)
      !
      ! Solver and physical constants
      !
      call get_keyword(500, 'alpha',                           alfa,                            0.50)              ! CFL Courant factor
      call get_keyword(500, 'theta',                           theta,                           1.0)               ! semi-implicit theta; <1 adds smoothing
      call get_keyword(500, 'hmin_cfl',                        hmin_cfl,                        0.1)               ! minimum depth used in CFL check (m)
      call get_keyword(500, 'manning',                         manning,                         0.04)              ! uniform Manning n (s/m^(1/3))
      call get_keyword(500, 'manning_land',                    manning_land,                    -999.0)            ! Manning n above rghlevland (s/m^(1/3))
      call get_keyword(500, 'manning_sea',                     manning_sea,                     -999.0)            ! Manning n below rghlevland (s/m^(1/3))
      call get_keyword(500, 'rgh_lev_land',                    rghlevland,                      0.0)               ! bed level separating land/sea friction (m)
      call get_keyword(500, 'zsini',                           zini,                            0.0)               ! initial water level (m)
      call get_keyword(500, 'qinf',                            qinf,                            0.0)               ! uniform infiltration rate (mm/hr); converted below
      call get_keyword(500, 'dtmax',                           dtmax,                           60.0)              ! upper bound on computational dt (s)
      call get_keyword(500, 'huthresh',                        huthresh,                        0.05)              ! wet/dry depth threshold (m)
      call get_keyword(500, 'huvmin',                          huvmin,                          0.0)               ! minimum depth for uv = q / max(hu, huvmin)
      call get_keyword(500, 'rhoa',                            rhoa,                            1.25)              ! air density (kg/m3)
      call get_keyword(500, 'rhow',                            rhow,                            1024.0)            ! water density (kg/m3)
      call get_keyword(500, 'inputformat',                     inputtype,                       'bin')             ! legacy bin/asc toggle for binary inputs
      call get_keyword(500, 'outputformat',                    outputtype,                      'net')             ! global output format (bin/asc/net)
      call get_keyword(500, 'outputtype_map',                  outputtype_map,                  'nil')             ! map-file output format (nil = follow outputformat)
      call get_keyword(500, 'outputtype_his',                  outputtype_his,                  'nil')             ! his-file output format (nil = follow outputformat)
      call get_keyword(500, 'nc_deflate_level',                nc_deflate_level,                2)                 ! netCDF deflate level (0-9)
      call get_keyword(500, 'bndtype',                         bndtype,                         1)                 ! boundary condition type
      call get_keyword(500, 'advection',                       advection,                       .true.)            ! enable momentum advection terms
      call get_keyword(500, 'latitude',                        latitude,                        0.0)               ! reference latitude for projected Coriolis (deg)
      call get_keyword(500, 'pavbnd',                          pavbnd,                          0.0)               ! atmospheric pressure applied at boundary (Pa)
      call get_keyword(500, 'gapres',                          gapres,                          101200.0)          ! atmospheric reference pressure (Pa)
      call get_keyword(500, 'baro',                            baro,                            1)                 ! include atmospheric-pressure gradient (1=on, 0=off)
      call get_keyword(500, 'utmzone',                         utmzone,                         'nil')             ! UTM zone string (e.g. '17N')
      call get_keyword(500, 'epsg',                            epsg,                            0)                 ! EPSG integer code for the grid
      call get_keyword(500, 'epsg',                            epsg_code,                       'nil')             ! EPSG as string (fallback)
      call get_keyword(500, 'advlim',                          advlim,                          1.0)               ! cap on advection term
      call get_keyword(500, 'slopelim',                        slopelim,                        9999.9)            ! cap on bed-slope water-level gradient
      call get_keyword(500, 'qinf_zmin',                       qinf_zmin,                       0.0)               ! minimum bed level for infiltration to apply (m)
      call get_keyword(500, 'btfilter',                        btfilter,                        60.0)              ! bathtub filter time scale (s)
      call get_keyword(500, 'sfacinf',                         sfacinf,                         0.2)               ! SCS initial-abstraction fraction (0.2S)
      call get_keyword(500, 'radstr',                          radstr,                          .false.)           ! radiation-stress forcing from SnapWave
      call get_keyword(500, 'crsgeo',                          crsgeo,                          .false.)           ! interpret grid coords as geographic (WGS84)
      call get_keyword(500, 'coriolis',                        coriolis,                        .true.)            ! include Coriolis force
      call get_keyword(500, 'amprblock',                       ampr_block,                      .true.)            ! treat 2D rainfall as block (true) or linearly interpolated (false)
      call get_keyword(500, 'spwmergefrac',                    spw_merge_frac,                  0.5)               ! merge factor for spiderweb wind composite
      call get_keyword(500, 'usespwprecip',                    use_spw_precip,                  .true.)            ! use precipitation field from spiderweb file
      call get_keyword(500, 'global',                          global,                          .false.)           ! treat grid as global (wrap in x)
      call get_keyword(500, 'nuvisc',                          nuviscdim,                       0.01)              ! viscosity coefficient (m2/s)
      call get_keyword(500, 'viscosity',                       viscosity,                       .false.)           ! enable horizontal viscosity term
      call get_keyword(500, 'spinup_meteo',                    spinup_meteo,                    .false.)           ! ramp wind/pressure from zero during tspinup
      call get_keyword(500, 'waveage',                         waveage,                         -999.0)            ! wave age (for SnapWave wind growth)
      call get_keyword(500, 'snapwave',                        snapwave,                        .false.)           ! enable coupled SnapWave wave solver
      call get_keyword(500, 'dtoutfixed',                      fixed_output_intervals,          .true.)            ! snap map/his to exact intervals (true) or let them drift with dt (false)
      !
      ! Wave maker parameters. Each call tries the modern `wavemaker_*`
      ! keyword first and falls back to the legacy shortname in the
      ! trailing array; a deprecation warning is emitted per legacy
      ! match.
      !
      call get_keyword(500, 'wavemaker_wvmfile',               wavemaker_wvmfile,               'none',    ['wvmfile'])    ! wavemaker polyline file
      call get_keyword(500, 'wavemaker_wfpfile',               wavemaker_wfpfile,               'none',    ['wfpfile'])    ! wavemaker forcing points file
      call get_keyword(500, 'wavemaker_whifile',               wavemaker_whifile,               'none',    ['whifile'])    ! wavemaker wave-height time series file
      call get_keyword(500, 'wavemaker_wtifile',               wavemaker_wtifile,               'none',    ['wtifile'])    ! wavemaker wave-period time series file
      call get_keyword(500, 'wavemaker_wstfile',               wavemaker_wstfile,               'none',    ['wstfile'])    ! wavemaker wave set-up time series file
      call get_keyword(500, 'wavemaker_filter_time',           wavemaker_filter_time,           600.0,     ['wmtfilter'])  ! wavemaker filter time scale (s)
      call get_keyword(500, 'wavemaker_filter_fred',           wavemaker_filter_fred,           0.99,      ['wmfred'])     ! wavemaker filter fred
      call get_keyword(500, 'wavemaker_signal',                wmsigstr,                        'spectrum',['wmsignal'])   ! wavemaker signal type (spectrum or monochromatic)
      call get_keyword(500, 'wavemaker_hmin',                  wavemaker_hmin,                  0.1,       ['wmhmin'])     ! wavemaker minimum depth for wave generation (m)
      call get_keyword(500, 'wavemaker_nfreqs_inc',            wavemaker_nfreqs_inc,            100,       ['nfreqsinc'])  ! wavemaker number of incident-wave frequencies
      call get_keyword(500, 'wavemaker_freqmin_inc',           wavemaker_freqmin_inc,           0.04,      ['freqmininc']) ! wavemaker incident-wave min frequency (Hz)
      call get_keyword(500, 'wavemaker_freqmax_inc',           wavemaker_freqmax_inc,           1.0,       ['freqmaxinc']) ! wavemaker incident-wave max frequency (Hz)
      call get_keyword(500, 'wavemaker_nfreqs_ig',             wavemaker_nfreqs_ig,             100,       ['nfreqsig'])   ! wavemaker number of IG-wave frequencies
      call get_keyword(500, 'wavemaker_freqmin_ig',            wavemaker_freqmin_ig,            0.0,       ['freqminig'])  ! wavemaker IG-wave min frequency (Hz)
      call get_keyword(500, 'wavemaker_freqmax_ig',            wavemaker_freqmax_ig,            0.1,       ['freqmaxig'])  ! wavemaker IG-wave max frequency (Hz)
      call get_keyword(500, 'wavemaker_tinc2ig',               wavemaker_tinc2ig,               -1.0)                      ! wavemaker ig/inc period ratio (<=0 uses Herbers)
      call get_keyword(500, 'wavemaker_surfslope',             wavemaker_surfslope,             -1.0)                      ! wavemaker surf-zone slope for empirical Tp_ig (van Ormondt et al., 2021)
      call get_keyword(500, 'wavemaker_hm0_ig_factor',         wavemaker_hm0_ig_factor,         1.0)                       ! wavemaker IG Hm0 scaling factor
      call get_keyword(500, 'wavemaker_hm0_inc_factor',        wavemaker_hm0_inc_factor,        1.0)                       ! wavemaker incident Hm0 scaling factor
      call get_keyword(500, 'wavemaker_gammax',                wavemaker_gammax,                1.0)                       ! wavemaker maximum Hrms/h
      call get_keyword(500, 'wavemaker_tpmin',                 wavemaker_tpmin,                 1.0)                       ! wavemaker minimum Tp (s)
      call get_keyword(500, 'wavemaker_hig',                   wavemaker_hig,                   .true.)                    ! wavemaker include IG waves
      call get_keyword(500, 'wavemaker_hinc',                  wavemaker_hinc,                  .false.)                   ! wavemaker include incident waves
      !
      ! Numerical parameters
      !
      call get_keyword(500, 'advection_scheme',                advstr,                          'upw1')            ! advection scheme label ('upw1' = 1st-order upwind, 'original' = legacy)
      call get_keyword(500, 'btrelax',                         btrelax,                         3600.0)            ! bathtub relaxation time (s)
      call get_keyword(500, 'wiggle_suppression',              wiggle_suppression,              .true.)            ! suppress spurious free-surface oscillations
      call get_keyword(500, 'structure_relax',                 structure_relax,                 4.0)               ! drainage-structure state-machine smoothing steps
      call get_keyword(500, 'wiggle_factor',                   wiggle_factor,                   0.1)               ! wiggle-suppression amplitude factor
      call get_keyword(500, 'wiggle_threshold',                wiggle_threshold,                0.1)               ! wiggle-suppression trigger threshold
      call get_keyword(500, 'uvlim',                           uvlim,                           10.0)              ! clipping velocity for momentum (m/s)
      call get_keyword(500, 'uvmax',                           uvmax,                           1000.0)            ! error-trigger velocity for momentum (m/s)
      call get_keyword(500, 'friction2d',                      friction2d,                      .true.)            ! apply friction at every UV point (true) or cell-wise (false)
      call get_keyword(500, 'advection_mask',                  advection_mask,                  .true.)            ! mask advection near dry cells
      call get_keyword(500, 'nuviscfac',                       nuviscfac,                       100.0)             ! multiplier on nuvisc near "difficult" points
      call get_keyword(500, 'nonh',                            nonhydrostatic,                  .false.)           ! enable non-hydrostatic pressure corrector
      call get_keyword(500, 'nh_fnudge',                       nh_fnudge,                       0.9)               ! non-hydrostatic nudging factor
      call get_keyword(500, 'nh_tstop',                        nh_tstop,                        -999.0)            ! non-hydrostatic stop time (s rel. tref); -999 = t1+999
      call get_keyword(500, 'nh_tol',                          nh_tol,                          0.001)             ! non-hydrostatic solver tolerance
      call get_keyword(500, 'nh_itermax',                      nh_itermax,                      100)               ! non-hydrostatic solver max iterations
      call get_keyword(500, 'h73table',                        h73table,                        .false.)           ! tabulate h^(7/3) for friction
      call get_keyword(500, 'rugdepth',                        runup_gauge_depth,               0.05)              ! runup gauge trigger depth (m)
      call get_keyword(500, 'wave_enhanced_roughness',         wave_enhanced_roughness,         .false.)           ! augment bed roughness with wave orbital velocity
      call get_keyword(500, 'use_bcafile',                     use_bcafile,                     .true.)            ! use tidal components from bca file
      call get_keyword(500, 'factor_wind',                     factor_wind,                     1.0)               ! scaling factor on wind forcing
      call get_keyword(500, 'factor_pres',                     factor_pres,                     1.0)               ! scaling factor on atmospheric pressure
      call get_keyword(500, 'factor_prcp',                     factor_prcp,                     1.0)               ! scaling factor on precipitation
      call get_keyword(500, 'factor_spw_size',                 factor_spw_size,                 1.0)               ! scaling factor on spiderweb radius
      call get_keyword(500, 'bathtub',                         bathtub,                         .false.)           ! run in bathtub (no momentum) mode
      call get_keyword(500, 'bathtub_fachs',                   bathtub_fac_hs,                  0.2)               ! bathtub Hs multiplier
      call get_keyword(500, 'bathtub_dt',                      bathtub_dt,                      -999.0)            ! bathtub time step (s); -999 = use dtmapout
      !
      ! Domain files
      !
      call get_keyword(500, 'qtrfile',                         qtrfile,                         'none')            ! quadtree netCDF file
      call get_keyword(500, 'depfile',                         depfile,                         'none')            ! bed-level (depth) file
      call get_keyword(500, 'inifile',                         zsinifile,                       'none')            ! initial water-level file
      call get_keyword(500, 'rstfile',                         rstfile,                         'none')            ! restart input file
      call get_keyword(500, 'mskfile',                         mskfile,                         'none')            ! active-cell mask file
      call get_keyword(500, 'indexfile',                       indexfile,                       'none')            ! index-to-active-cell mapping file
      call get_keyword(500, 'cstfile',                         cstfile,                         'none')            ! coastline polyline file
      call get_keyword(500, 'sbgfile',                         sbgfile,                         'none')            ! subgrid tables netCDF file
      call get_keyword(500, 'thdfile',                         thdfile,                         'none')            ! thin dams polyline file
      call get_keyword(500, 'weirfile',                        weirfile,                        'none')            ! weirs polyline file
      call get_keyword(500, 'manningfile',                     manningfile,                     'none')            ! spatially-varying Manning n file
      call get_keyword(500, 'drnfile',                         drnfile,                         'none')            ! drainage-structures (pumps/gates/culverts) TOML file
      call get_keyword(500, 'urbfile',                         urbfile,                         'none')            ! urban drainage zones TOML file
      call get_keyword(500, 'volfile',                         volfile,                         'none')            ! depression-storage volume file
      !
      ! Forcing files (ascii / binary)
      !
      call get_keyword(500, 'bndfile',                         bndfile,                         'none')            ! water-level boundary points
      call get_keyword(500, 'bzsfile',                         bzsfile,                         'none')            ! water-level boundary time series
      call get_keyword(500, 'bcafile',                         bcafile,                         'none')            ! tidal components per boundary point
      call get_keyword(500, 'bzifile',                         bzifile,                         'none')            ! IG wave boundary time series
      call get_keyword(500, 'bdrfile',                         bdrfile,                         'none')            ! downstream river boundary file
      call get_keyword(500, 'srcfile',                         srcfile,                         'none')            ! river-point source locations
      call get_keyword(500, 'disfile',                         disfile,                         'none')            ! river-point discharge time series
      call get_keyword(500, 'spwfile',                         spwfile,                         'none')            ! spiderweb tropical-cyclone file
      call get_keyword(500, 'wndfile',                         wndfile,                         'none')            ! uniform wind time series
      call get_keyword(500, 'prcfile',                         prcpfile,                        'none',    ['precipfile']) ! uniform precipitation time series
      call get_keyword(500, 'amufile',                         amufile,                         'none')            ! 2D wind u-component file
      call get_keyword(500, 'amvfile',                         amvfile,                         'none')            ! 2D wind v-component file
      call get_keyword(500, 'ampfile',                         ampfile,                         'none')            ! 2D atmospheric pressure file
      call get_keyword(500, 'amrfile',                         amprfile,                        'none',    ['amprfile'])            ! 2D precipitation rate file
      call get_keyword(500, 'z0lfile',                         z0lfile,                         'none')            ! 2D land roughness (z0) file
      !
      ! NetCDF-format forcing files (FEWS-style)
      !
      call get_keyword(500, 'netbndbzsbzifile',                netbndbzsbzifile,                'none')            ! combined bnd/bzs/bzi netCDF file
      call get_keyword(500, 'netsrcdisfile',                   netsrcdisfile,                   'none')            ! combined src/dis netCDF file
      call get_keyword(500, 'netamuamvfile',                   netamuamvfile,                   'none')            ! combined amu/amv netCDF file
      call get_keyword(500, 'netamprfile',                     netamprfile,                     'none')            ! 2D precipitation netCDF file
      call get_keyword(500, 'netampfile',                      netampfile,                      'none')            ! 2D pressure netCDF file
      call get_keyword(500, 'netspwfile',                      netspwfile,                      'none')            ! netCDF spiderweb file
      !
      ! Infiltration and losses
      !
      call get_keyword(500, 'infiltrationfile',                infiltrationfile,                'none')            ! infiltration parameters TOML file
      call get_keyword(500, 'infiltrationtype',                inftype,                         'none')            ! infiltration flavor (con, c2d, cna, cnb, gai, hor, bkt)
      !
      ! Legacy binary infiltration inputs (kept for backward compatibility).
      !
      call get_keyword(500, 'qinffile',                        qinffile,                        'none')            ! binary spatially-varying infiltration field
      call get_keyword(500, 'scsfile',                         scsfile,                         'none')            ! SCS curve-number S field (legacy binary)
      call get_keyword(500, 'smaxfile',                        smaxfile,                        'none')            ! SCS max storage S field (legacy binary)
      call get_keyword(500, 'sefffile',                        sefffile,                        'none')            ! SCS effective storage S_e field (legacy binary)
      call get_keyword(500, 'psifile',                         psifile,                         'none')            ! Green-Ampt suction head (legacy binary, mm)
      call get_keyword(500, 'sigmafile',                       sigmafile,                       'none')            ! Green-Ampt maximum moisture deficit (legacy binary)
      call get_keyword(500, 'ksfile',                          ksfile,                          'none')            ! Green-Ampt saturated hydraulic conductivity (legacy binary, mm/hr)
      call get_keyword(500, 'f0file',                          f0file,                          'none')            ! Horton initial infiltration capacity F0 (legacy binary)
      call get_keyword(500, 'fcfile',                          fcfile,                          'none')            ! Horton asymptotic infiltration rate Fc (legacy binary)
      call get_keyword(500, 'kdfile',                          kdfile,                          'none')            ! Horton decay constant k (legacy binary, 1/hr)
      call get_keyword(500, 'horton_kr_kd',                    horton_kr_kd,                    10.0)              ! Horton recovery/decay ratio
      !
      ! Output
      !
      call get_keyword(500, 'obsfile',                         obsfile,                         'none')            ! observation-point locations file
      call get_keyword(500, 'crsfile',                         crsfile,                         'none')            ! cross-section polyline file
      call get_keyword(500, 'rugfile',                         rugfile,                         'none')            ! runup-gauge locations file
      call get_keyword(500, 'store_maximum_waterlevel',        store_maximum_waterlevel,        .true.)            ! store maximum water level on dtmaxout interval (only if dtmaxout > 0)
      call get_keyword(500, 'storevelmax',                     store_maximum_velocity,          .false.)           ! store maximum flow velocity on dtmaxout interval (only if dtmaxout > 0)
      call get_keyword(500, 'storefluxmax',                    store_maximum_flux,              .false.)           ! store maximum flux on dtmaxout interval (only if dtmaxout > 0)
      call get_keyword(500, 'storevel',                        store_velocity,                  .false.)           ! store velocity on dtout interval
      call get_keyword(500, 'storecumprcp',                    store_cumulative_precipitation,  .false.)           ! store cumulative precipitation + infiltration on dtmaxout interval
      call get_keyword(500, 'storetwet',                       store_twet,                      .false.)           ! store per-cell wet duration
      call get_keyword(500, 'storetzsmax',                     store_t_zsmax,                   .false.)           ! store time stamp of zsmax occurrence
      call get_keyword(500, 'storehsubgrid',                   store_hsubgrid,                  .false.)           ! store hmax in subgrid mode (zsmax - subgrid_z_zmin)
      call get_keyword(500, 'storehmean',                      store_hmean,                     .false.)           ! store hmax as subgrid-mean depth instead of max (requires storehsubgrid)
      call get_keyword(500, 'twet_threshold',                  twet_threshold,                  0.01)              ! water-depth threshold counting a cell as wet (storetwet)
      call get_keyword(500, 'store_tsunami_arrival_time',      store_tsunami_arrival_time,      .false.)           ! store tsunami arrival time per cell
      call get_keyword(500, 'tsunami_arrival_threshold',       tsunami_arrival_threshold,       0.01)              ! water-depth threshold for tsunami arrival
      call get_keyword(500, 'timestep_analysis',               timestep_analysis,               .false.)           ! write per-cell timestep limiter diagnostics
      call get_keyword(500, 'storeqdrain',                     store_qdrain,                    .true.)            ! store per-drainage-structure discharge in his file
      call get_keyword(500, 'store_river_discharge',           store_river_discharge,           .false.)           ! store per-river-point discharge in his file
      call get_keyword(500, 'store_urban_drainage_discharge',  store_urban_drainage_discharge,  .false.)           ! store per-urban-zone outfall discharge in his file
      call get_keyword(500, 'store_cumulative_urban_drainage', store_cumulative_urban_drainage, .false.)           ! store cumulative urban drainage depth per cell in map file
      call get_keyword(500, 'storezvolume',                    store_zvolume,                   .false.)           ! store subgrid cell volume (requires subgrid)
      call get_keyword(500, 'storestoragevolume',              store_storagevolume,             .false.)           ! store remaining storage volume (requires subgrid + volfile)
      call get_keyword(500, 'writeruntime',                    write_time_output,               .false.)           ! write runtimes.txt at end of simulation
      call get_keyword(500, 'debug',                           debug,                           .false.)           ! debug output at every time step
      call get_keyword(500, 'storemeteo',                      store_meteo,                     .false.)           ! store 2D meteo forcing fields in map file
      call get_keyword(500, 'storemaxwind',                    store_wind_max,                  .false.)           ! store maximum wind speed (requires storemeteo)
      call get_keyword(500, 'storefw',                         store_wave_forces,               .false.)           ! store wave-radiation forces
      call get_keyword(500, 'storewavdir',                     store_wave_direction,            .false.)           ! store wave direction
      call get_keyword(500, 'regular_output_on_mesh',          use_quadtree_output,             .false.)           ! write quadtree output on regular m/n mesh
      call get_keyword(500, 'store_dynamic_bed_level',         store_dynamic_bed_level,         .false.)           ! store time-varying bed level (subgrid)
      call get_keyword(500, 'snapwave_use_nearest',            snapwave_use_nearest,            .true.)            ! use nearest-neighbour lookup for SnapWave boundary points
      call get_keyword(500, 'percentage_done',                 percdoneval,                     5)                 ! progress-reporter interval (% complete)
      !
      ! Coupled SnapWave solver parameters
      !
      call get_keyword(500, 'snapwave_wind',                   snapwavewind,                    .false.)           ! feed wind into SnapWave (implies storing wind speed/direction)
      call get_keyword(500, 'snapwave_waveforces_factor',      waveforces_factor,               1.0)               ! multiplier on SnapWave wave forces
      !
      ! Wind drag coefficient table
      !
      call get_keyword(500, 'cdnrb',                           cd_nr,                           0)                 ! number of wind-drag breakpoints (0 = use defaults)
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
         call get_keyword_real_array(500, 'cdwnd', cd_wnd, 0.0, cd_nr)
         call get_keyword_real_array(500, 'cdval', cd_val, 0.0, cd_nr)
         !
      endif
      !
      close(500)
      !
      ! Done with reading input
      !
      ! Now do some post-processing and consistency checks on the inputs, and emit
      !
      ! Limit progress reporter to (0, 100]%
      !
      percdoneval = max(min(percdoneval, 100), 0)
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
      if (dtmaxout > 0.0) store_maximum_waterlevel = .true.
      !
      ! Apply gates to the flags now that the full set of inputs has been read.
      !
      if (dtmaxout <= 0.0) then
         !
         ! Are there more to be added here?
         !
         store_maximum_waterlevel = .false.
         store_maximum_velocity   = .false.
         store_maximum_flux       = .false.
         !
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
      ! Map/his format fallback: inherit the global outputformat when either
      ! per-file format was left at 'nil'.
      !
      if ((outputtype_map == 'nil') .or. (outputtype_his == 'nil')) then
         outputtype_map = outputtype
         outputtype_his = outputtype
      endif
      !
      if (sbgfile(1:4) /= 'none') then
         !
         subgrid = .true.
         !
      else
         !
         subgrid = .false.
         !
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
         if (wmsigstr(1:3) == 'mon') then
            !
            wavemaker_spectrum = .false.
            !
            call write_log('Info    : use monochromatic wave spectrum for wave makers', 0)
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
   subroutine get_keyword_real(fileid, keyword, value, default, legacy)
      !
      ! Read one real*4 keyword. Tries `keyword` first; if absent, walks
      ! the optional `legacy` list of deprecated aliases and emits a
      ! one-line deprecation warning per matched alias. Falls back to
      ! `default` when nothing matches.
      !
      ! Called from: read_sfincs_input (this module).
      !
      implicit none
      !
      integer,                           intent(in)            :: fileid
      character(*),                      intent(in)            :: keyword
      real*4,                            intent(out)           :: value
      real*4,                            intent(in)            :: default
      character(len=*), dimension(:),    intent(in), optional  :: legacy
      !
      character(len=256) :: valstr
      logical            :: found
      integer            :: i
      !
      call find_value(fileid, keyword, valstr, found)
      if (found) then
         read(valstr, *) value
         return
      endif
      !
      if (present(legacy)) then
         !
         do i = 1, size(legacy)
            !
            call find_value(fileid, trim(legacy(i)), valstr, found)
            !
            if (found) then
               read(valstr, *) value
               call warn_legacy(trim(legacy(i)), keyword)
               return
            endif
            !
         enddo
         !
      endif
      !
      value = default
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine get_keyword_int(fileid, keyword, value, default, legacy)
      !
      ! Read one integer keyword. See get_keyword_real for the semantics.
      !
      ! Called from: read_sfincs_input (this module).
      !
      implicit none
      !
      integer,                           intent(in)            :: fileid
      character(*),                      intent(in)            :: keyword
      integer,                           intent(out)           :: value
      integer,                           intent(in)            :: default
      character(len=*), dimension(:),    intent(in), optional  :: legacy
      !
      character(len=256) :: valstr
      logical            :: found
      integer            :: i
      !
      call find_value(fileid, keyword, valstr, found)
      if (found) then
         read(valstr, *) value
         return
      endif
      !
      if (present(legacy)) then
         !
         do i = 1, size(legacy)
            !
            call find_value(fileid, trim(legacy(i)), valstr, found)
            !
            if (found) then
               read(valstr, *) value
               call warn_legacy(trim(legacy(i)), keyword)
               return
            endif
            !
         enddo
         !
      endif
      !
      value = default
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine get_keyword_char(fileid, keyword, value, default, legacy)
      !
      ! Read one character-string keyword. See get_keyword_real for the
      ! semantics. The entire right-hand side (after trailing comments
      ! are stripped) becomes `value`.
      !
      ! Called from: read_sfincs_input (this module).
      !
      implicit none
      !
      integer,                           intent(in)            :: fileid
      character(*),                      intent(in)            :: keyword
      character(*),                      intent(out)           :: value
      character(*),                      intent(in)            :: default
      character(len=*), dimension(:),    intent(in), optional  :: legacy
      !
      character(len=256) :: valstr
      logical            :: found
      integer            :: i
      !
      call find_value(fileid, keyword, valstr, found)
      if (found) then
         value = valstr
         return
      endif
      !
      if (present(legacy)) then
         !
         do i = 1, size(legacy)
            !
            call find_value(fileid, trim(legacy(i)), valstr, found)
            !
            if (found) then
               value = valstr
               call warn_legacy(trim(legacy(i)), keyword)
               return
            endif
            !
         enddo
         !
      endif
      !
      value = default
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine get_keyword_logical(fileid, keyword, value, default, legacy)
      !
      ! Read one logical keyword. Accepts `1`, `y`, `Y`, `t`, `T` as true;
      ! anything else (including absence → `default`, and `0`, `n`, `N`,
      ! `f`, `F`) as false.
      !
      ! Called from: read_sfincs_input (this module).
      !
      implicit none
      !
      integer,                           intent(in)            :: fileid
      character(*),                      intent(in)            :: keyword
      logical,                           intent(out)           :: value
      logical,                           intent(in)            :: default
      character(len=*), dimension(:),    intent(in), optional  :: legacy
      !
      character(len=256) :: valstr
      logical            :: found
      integer            :: i
      !
      call find_value(fileid, keyword, valstr, found)
      if (found) then
         value = parse_logical(valstr)
         return
      endif
      !
      if (present(legacy)) then
         !
         do i = 1, size(legacy)
            !
            call find_value(fileid, trim(legacy(i)), valstr, found)
            !
            if (found) then
               value = parse_logical(valstr)
               call warn_legacy(trim(legacy(i)), keyword)
               return
            endif
            !
         enddo
         !
      endif
      !
      value = default
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine get_keyword_real_array(fileid, keyword, value, default, nr, legacy)
      !
      ! Read one whitespace-separated real*4 array keyword. Allocates
      ! `value(nr)` on the way in and fills it from the matching line.
      ! Same fallback semantics as get_keyword_real.
      !
      ! Called from: read_sfincs_input (this module).
      !
      implicit none
      !
      integer,                              intent(in)            :: fileid
      character(*),                         intent(in)            :: keyword
      integer,                              intent(in)            :: nr
      real*4,                               intent(in)            :: default
      real*4, dimension(:),                 intent(out), allocatable :: value
      character(len=*), dimension(:),       intent(in),  optional :: legacy
      !
      character(len=256) :: valstr
      logical            :: found
      integer            :: i, m
      !
      allocate(value(nr))
      !
      call find_value(fileid, keyword, valstr, found)
      if (found) then
         read(valstr, *) (value(m), m = 1, nr)
         return
      endif
      !
      if (present(legacy)) then
         !
         do i = 1, size(legacy)
            !
            call find_value(fileid, trim(legacy(i)), valstr, found)
            !
            if (found) then
               read(valstr, *) (value(m), m = 1, nr)
               call warn_legacy(trim(legacy(i)), keyword)
               return
            endif
            !
         enddo
         !
      endif
      !
      value = default
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine find_value(fileid, keyword, valstr, found)
      !
      ! Scan an already-open sfincs.inp once for the given `keyword`.
      ! Returns the raw right-hand-side value string and whether the key
      ! was matched.
      !
      ! Called from: get_keyword_real / get_keyword_int / get_keyword_char /
      ! get_keyword_logical / get_keyword_real_array.
      !
      implicit none
      !
      integer,      intent(in)  :: fileid
      character(*), intent(in)  :: keyword
      character(*), intent(out) :: valstr
      logical,      intent(out) :: found
      !
      character(len=256) :: keystr
      character(len=256) :: line
      integer            :: stat
      !
      found  = .false.
      valstr = ''
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
            found = .true.
            return
         endif
         !
      enddo
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   subroutine warn_legacy(legacy_key, new_key)
      !
      ! Emit a one-line deprecation warning to the log. Called whenever
      ! a legacy keyword alias was matched; the user can silence this
      ! by migrating the keyword in their sfincs.inp.
      !
      ! Called from: get_keyword_real / get_keyword_int / get_keyword_char /
      ! get_keyword_logical / get_keyword_real_array.
      !
      implicit none
      !
      character(*), intent(in) :: legacy_key
      character(*), intent(in) :: new_key
      !
      character(len=512) :: msg
      !
      write(msg, '(a,a,a,a,a)') ' Warning : sfincs.inp keyword "', trim(legacy_key), &
         '" is deprecated, use "', trim(new_key), '" instead'
      call write_log(trim(msg), 1)
      !
   end subroutine
   !
   !-----------------------------------------------------------------------------------------------------!
   !
   function parse_logical(valstr) result(value)
      !
      ! Map an sfincs.inp value string to a logical. `1`, `y`, `Y`, `t`,
      ! `T` at position 1 are true; everything else is false.
      !
      ! Called from: get_keyword_logical (this module).
      !
      implicit none
      !
      character(*), intent(in) :: valstr
      logical                  :: value
      !
      value = (valstr(1:1) == '1' .or. valstr(1:1) == 'y' .or. valstr(1:1) == 'Y' .or. &
               valstr(1:1) == 't' .or. valstr(1:1) == 'T')
      !
   end function
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
      ! Called from: find_value.
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
