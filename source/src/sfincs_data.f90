module sfincs_data

      !!! Revision data
      character*256 :: build_revision, build_date
      !!!!
      !!! Time variables
      real    :: tstart_all, tfinish_all
      real*4  :: dtavg
      !!!
      !!! Constants
      !!!
      real*4 g                                   ! gravitational constant g
      real*4 pi                                  ! pi
      !!!
      !!! Input data from sfincs.inp
      !!!
      real*4 alfa                                ! courant factor
      real*4 manning                             ! constant manning
      real*4 manning_land                        ! manning on land
      real*4 manning_sea                         ! manning on water
      real*4 rghlevland
      real*4 gn2land
      real*4 gn2sea
      real*4 gn2
      real*4 t0
      real*4 t1
      real*4 dx
      real*4 dy
      real*4 dxinv
      real*4 dyinv
      real*4 area
      real*4 dxy
      real*4 t0out
      real*4 t1out
      real*4 dtmapout
      real*4 dtmaxout
      real*4 trst
      real*4 dtrstout
      real*4 dthisout
      real*4 dtwindupd
      real*4 theta
      real*4 dtmax
      real*4 dtmin
      real*4 zini
      real*4 x0
      real*4 y0
      real*4 rotation
      real*4 huthresh
      real*4 max_elev
      real*4 qinf
      real*4 tig
      real*4 tspinup
      real*4 cosrot
      real*4 sinrot
      real*4 rhoa
      real*4 rhow
      real*4 latitude
      real*4 fcorio
      real*4 pavbnd
      real*4 gapres
      real*4 stopdepth
      real*4 advlim
      real*4 twet_threshold
      real*4, dimension(:), allocatable :: cd_wnd
      real*4, dimension(:), allocatable :: cd_val
      real*4 qinf_zmin
      real*4 btfilter
      real*4 sfacinf
      real*4 dym
      real*4 dym2
      real*4 nuvisc
      real*4 nuviscdim
      real*4 nuviscinp      
      real*4 spw_merge_frac
      real*4 tsunami_arrival_threshold
      real*4 dtwave
      real*4 wmtfilter
      real*4 horton_kr_kd
      !
      real*4 freqminig
      real*4 freqmaxig
      integer nfreqsig
      !
      integer mmax
      integer nmax
      integer bndtype
      integer cd_nr
      integer baro
      !
      character*256 :: depfile
      character*256 :: mskfile
      character*256 :: bndfile
      character*256 :: bzsfile
      character*256 :: bzifile
      character*256 :: bwvfile
      character*256 :: bhsfile
      character*256 :: btpfile
      character*256 :: bwdfile
      character*256 :: bdsfile
      character*256 :: wfpfile
      character*256 :: whifile
      character*256 :: wtifile
      character*256 :: wstfile
      character*256 :: obsfile
      character*256 :: crsfile
      character*256 :: srcfile
      character*256 :: disfile
      character*256 :: drnfile
      character*256 :: zsinifile
      character*256 :: rstfile
      character*256 :: indexfile
      character*256 :: bindepfile
      character*256 :: binmskfile
      character*256 :: manningfile
      character*256 :: cstfile
      character*256 :: spwfile
      character*256 :: amufile
      character*256 :: amvfile
      character*256 :: ampfile
      character*256 :: amprfile
      character*256 :: wndfile
      character*256 :: prcpfile
      character*256 :: sbgfile
      character*256 :: thdfile
      character*256 :: weirfile
      character*256 :: qinffile
      character*256 :: netbndbzsbzifile
      character*256 :: netsrcdisfile
      character*256 :: netamuamvfile
      character*256 :: netampfile
      character*256 :: netamprfile
      character*256 :: scsfile
      character*256 :: smaxfile
      character*256 :: sefffile
      character*256 :: psifile
      character*256 :: sigmafile
      character*256 :: ksfile
      character*256 :: f0file
      character*256 :: fcfile
      character*256 :: kdfile
      character*256 :: z0lfile
      character*256 :: wvmfile
      character*256 :: qtrfile
      character*256 :: volfile
      !
      character*256 :: trefstr_iso8601
      character*41  :: treftimefews
      character*15  :: trefstr
      character*15  :: tstartstr
      character*15  :: tstopstr
      !
      character*3   :: inputtype
      character*3   :: outputtype
      character*3   :: outputtype_map
      character*3   :: outputtype_his
      character*3   :: utmzone
      character*3   :: inftype
      integer       :: epsg
      character*15  :: epsg_code
      !
      logical       :: waves
      logical       :: wind
      logical       :: patmos
      logical       :: precip
      logical       :: spw_precip
      logical       :: meteo3d
      logical       :: subgrid
      logical       :: manning2d ! spatially-varying roughness
      logical       :: coriolis
      logical       :: store_cumulative_precipitation
      logical       :: store_maximum_waterlevel
      logical       :: store_maximum_waterdepth
      logical       :: store_maximum_velocity
      logical       :: store_maximum_flux      
      logical       :: store_velocity
      logical       :: store_twet
      logical       :: store_hsubgrid
      logical       :: store_qdrain
      logical       :: store_zvolume
      logical       :: store_meteo
      logical       :: store_wind_max
      logical       :: store_wave_forces
      logical       :: store_wave_direction
      logical       :: useqxy0
      logical       :: usehuv
      logical       :: usehuv0
      logical       :: write_time_output
      logical       :: bziwaves
      logical       :: infiltration
      logical       :: debug
      logical       :: radstr
      logical       :: crsgeo
      logical       :: use_coriolis
      logical       :: ampr_block
      logical       :: global
      logical       :: store_tsunami_arrival_time
      logical       :: viscosity
      logical       :: spinup_meteo
      logical       :: snapwave
      logical       :: wind_reduction
      logical       :: include_boundaries
      logical       :: use_uv
      logical       :: wavemaker
      logical       :: wavemaker_mobile
      logical       :: use_quadtree
      logical       :: interpolate_zst
      logical       :: advection
      logical       :: fixed_output_intervals
      logical       :: use_storage_volume
      !!!
      !!! sfincs_input.f90 switches
      integer storevelmax
      integer storefluxmax
      integer storevel
      integer storecumprcp
      integer storetwet
      integer storeqdrain
      integer storezvolume
      integer storemeteo
      integer storehsubgrid
      integer wrttimeoutput
      integer idebug
      integer iradstr
      integer igeo
      integer icoriolis
      integer iamprblock
      integer iglobal
      integer itsunamitime
      integer ispinupmeteo
      integer isnapwave
      integer iwindmax
      integer ioutfixed
      integer iadvection
      integer istorefw
      integer istorewavdir   
      integer imanning2d
      integer iviscosity   
      integer isubgrid  
      integer iwavemaker      
      integer iwavemaker_spectrum  
      integer ispwprecip
      !!!
      !!! Static data
      !!!
      integer*4 :: np
      integer*4 :: npuv
      integer*4 :: nkcuv2
      !
      ! Internal wave maker
      !
      integer*4 :: nkcs4
      integer*4 :: nkcuv4
      !
      ! Indices
      !
      integer*4,          dimension(:),   allocatable :: nmindbnd
!      integer*4,          dimension(:,:), allocatable :: z_index
      integer*4,          dimension(:), allocatable   :: z_index_z_n
      integer*4,          dimension(:), allocatable   :: z_index_z_m
!      integer*4,          dimension(:), allocatable   :: z_index_z_nm
      integer*4,          dimension(:),   allocatable :: z_index_uv_md1
      integer*4,          dimension(:),   allocatable :: z_index_uv_md2
      integer*4,          dimension(:),   allocatable :: z_index_uv_mu1
      integer*4,          dimension(:),   allocatable :: z_index_uv_mu2
      integer*4,          dimension(:),   allocatable :: z_index_uv_nd1
      integer*4,          dimension(:),   allocatable :: z_index_uv_nd2
      integer*4,          dimension(:),   allocatable :: z_index_uv_nu1
      integer*4,          dimension(:),   allocatable :: z_index_uv_nu2
      !
      integer*4,          dimension(:),   allocatable :: uv_index_z_nm
      integer*4,          dimension(:),   allocatable :: uv_index_z_nmu
      integer*4,          dimension(:),   allocatable :: uv_index_u_nmu
      integer*4,          dimension(:),   allocatable :: uv_index_u_nmd
      integer*4,          dimension(:),   allocatable :: uv_index_u_num
      integer*4,          dimension(:),   allocatable :: uv_index_u_ndm
      integer*4,          dimension(:),   allocatable :: uv_index_v_ndm
      integer*4,          dimension(:),   allocatable :: uv_index_v_nm
      integer*4,          dimension(:),   allocatable :: uv_index_v_ndmu
      integer*4,          dimension(:),   allocatable :: uv_index_v_nmu
      !
      ! Flags
      !
      integer*1,          dimension(:),   allocatable :: z_flags_iref
      integer*1,          dimension(:),   allocatable :: z_flags_type
      !
      integer*1,          dimension(:),   allocatable :: uv_flags_iref
      integer*1,          dimension(:),   allocatable :: uv_flags_type
      integer*1,          dimension(:),   allocatable :: uv_flags_dir
      integer*1,          dimension(:),   allocatable :: uv_flags_adv
      integer*1,          dimension(:),   allocatable :: uv_flags_vis      
      !
      integer*1,          dimension(:),   allocatable :: kcs
      integer*1,          dimension(:),   allocatable :: kcuv
      integer*1,          dimension(:),   allocatable :: kfuv
      integer*1,          dimension(:),   allocatable :: scs_rain   ! logic if previous time step was raining
      !
      ! Quadtree
      !
      integer*1 :: nref
      integer*4, dimension(:),   allocatable :: index_sfincs_in_quadtree
      integer*4, dimension(:),   allocatable :: index_quadtree_in_sfincs
      !
      ! Refinement levels
      !
      real*4, dimension(:),   allocatable :: dxr
      real*4, dimension(:),   allocatable :: dyr
      real*4, dimension(:),   allocatable :: dxrm
      real*4, dimension(:),   allocatable :: dyrm
      real*4, dimension(:),   allocatable :: dxrinv
      real*4, dimension(:),   allocatable :: dyrinv
      real*4, dimension(:),   allocatable :: dxr2inv
      real*4, dimension(:),   allocatable :: dyr2inv
      real*4, dimension(:),   allocatable :: dxrinvc
      real*4, dimension(:),   allocatable :: dyrinvc
      real*4, dimension(:),   allocatable :: cell_area
      !
!      integer*4, dimension(:),   allocatable :: nr_points_per_level      
!      integer*4, dimension(:,:), allocatable :: nm_indices_per_level
!      integer*4, dimension(:,:), allocatable :: cell_indices_per_level
      !
      ! Cell sizes
      !
      real*4, dimension(:),   allocatable :: dxm
      real*4, dimension(:),   allocatable :: dxminv
      real*4, dimension(:),   allocatable :: dxm2inv
      !
      ! Z-points
      !
      real*4,             dimension(:),   allocatable, target :: z_xz
      real*4,             dimension(:),   allocatable, target :: z_yz
      real*4,             dimension(:),   allocatable :: cell_area_m2
      !
      ! UV-points
      !
!      integer*4,          dimension(:,:), allocatable :: uv_index
!      integer*4,          dimension(:), allocatable   :: uv_index_nm
!      integer*4,          dimension(:), allocatable   :: uv_index_nmu
!      integer*1,          dimension(:,:), allocatable :: uv_flags
      !
      real*4, dimension(:),   allocatable, target :: zb
      real*4, dimension(:),   allocatable :: zbuv
      real*4, dimension(:),   allocatable :: zbuvmx
      real*4, dimension(:),   allocatable :: gn2uv
      !
      ! Rainfall and infiltration
      !      
      real*4, dimension(:),   allocatable :: qinffield       ! infiltration map prescribed by user and typically not updated
      real*4, dimension(:),   allocatable :: ksfield         ! saturated hydraulic conductivity
      real*4, dimension(:),   allocatable :: qinfmap         ! infiltration used in the computation
      real*4, dimension(:),   allocatable :: rain_T1         ! time that it is not raining - used both in cnb and gai
      real*4, dimension(:),   allocatable :: inf_kr          ! recovery concept - used both in cnb and gai
      real*4, dimension(:),   allocatable :: scs_Se          ! effective S (Se) used in SCS method
      real*4, dimension(:),   allocatable :: scs_P1          ! cumulative rainfall of this 'event'
      real*4, dimension(:),   allocatable :: scs_F1          ! infiltration of this 'event'
      real*4, dimension(:),   allocatable :: scs_S1          ! S for this 'event'
      real*4, dimension(:),   allocatable :: GA_head         ! the soil suction head in mm
      real*4, dimension(:),   allocatable :: GA_sigma_max    ! the maximum soil capacity (porosity) in [-]
      real*4, dimension(:),   allocatable :: GA_sigma        ! the current soil capacity (porosity) in [-]
      real*4, dimension(:),   allocatable :: GA_F            ! cumulative infiltration for green-ampt
      real*4, dimension(:),   allocatable :: GA_Lu           ! depth of upper soil recovery zone (computed from ksfield)
      real*4, dimension(:),   allocatable :: horton_fc       ! final infiltration capacity - Horton
      real*4, dimension(:),   allocatable :: horton_f0       ! initial infiltration capacity 
      real*4, dimension(:),   allocatable :: horton_kd       ! coefficient of the exponential term
      !
      ! Storage volume
      !
      real*4, dimension(:),   allocatable :: storage_volume  ! Storage volume green infra
      !
      ! Wind reduction for spiderweb winds
      !
      real*4, dimension(:,:), allocatable :: z0land          ! z0 values over land for spiderweb wind speed reduction   
      real*4, dimension(:),   allocatable :: z0land_table    ! z0land_over_z0marine look-up table
      !
      ! Coriolis
      !
      real*4, dimension(:),   allocatable :: fcorio2d
      !
      ! Boundary velocity points
      !
      integer*4, dimension(:),     allocatable :: index_kcuv2
      integer*4, dimension(:),     allocatable :: nmbkcuv2
      integer*4, dimension(:),     allocatable :: nmikcuv2
      integer*1, dimension(:),     allocatable :: ibuvdir
      !
      ! Mean velocity at boundaries
      !
      real*4, dimension(:),   allocatable :: uvmean
      !
      ! Wave maker points
      !
      integer*4                            :: wavemaker_nr_uv_points
      integer*4, dimension(:), allocatable :: wavemaker_index_uv
      integer*4, dimension(:), allocatable :: wavemaker_index_nmi
      integer*4, dimension(:), allocatable :: wavemaker_index_nmb
      real*4,    dimension(:), allocatable :: wavemaker_fac_wmfp
      integer*1, dimension(:), allocatable :: wavemaker_idir
      real*4,    dimension(:), allocatable :: wavemaker_angfac
      real*4,    dimension(:), allocatable :: wavemaker_uvmean
      real*4,    dimension(:), allocatable :: wavemaker_uvtrend
      integer*4, dimension(:), allocatable :: wavemaker_nmu
      integer*4, dimension(:), allocatable :: wavemaker_nmd
      integer*4, dimension(:), allocatable :: wavemaker_num
      integer*4, dimension(:), allocatable :: wavemaker_ndm
      integer*4, dimension(:), allocatable :: z_index_wavemaker
      real*4                               :: wavemaker_freduv      
      logical                              :: wavemaker_spectrum
      !
      ! Wave maker time series
      !
      integer*4                             :: nwmfp
      integer*4                             :: ntwmfp
      integer*4                             :: itwmfplast
      logical                               :: wavemaker_timeseries
      integer*4, dimension(:), allocatable  :: wavemaker_index_wmfp1
      integer*4, dimension(:), allocatable  :: wavemaker_index_wmfp2
      real*4, dimension(:),     allocatable :: x_wmfp
      real*4, dimension(:),     allocatable :: y_wmfp
      real*4, dimension(:),     allocatable :: wmf_time
      real*4, dimension(:,:),   allocatable :: wmf_hm0_ig
      real*4, dimension(:,:),   allocatable :: wmf_tp_ig
      real*4, dimension(:,:),   allocatable :: wmf_setup
      real*4, dimension(:),     allocatable :: wmf_hm0_ig_t
      real*4, dimension(:),     allocatable :: wmf_tp_ig_t
      real*4, dimension(:),     allocatable :: wmf_setup_t
      !
      
!      integer*4                              :: wavemaker_nr_cross
!      integer*4                              :: wavemaker_nr_along
!      real*4                                 :: wavemaker_dx_cross
!      real*4                                 :: wavemaker_dx_along
 !     real*8,    dimension(:),   allocatable :: wavemaker_cross_x
 !     real*8,    dimension(:),   allocatable :: wavemaker_cross_y
 !     real*8,    dimension(:,:), allocatable :: wavemaker_cross_w
 !     integer*4, dimension(:,:), allocatable :: wavemaker_cross_i
 !     real*4,    dimension(:),   allocatable :: wavemaker_cross_phi
      !
      ! NM
      !
!      integer*4, dimension(:),   allocatable :: wavemaker_index_z_in_nm     ! points to index of wm maker water level points (size=np) 
      !
      ! WMZ
      !
!      integer*4, dimension(:),   allocatable :: wavemaker_index_u0_in_z    ! points to index of wm maker u points to the left (size = nr wave maker water level points) 
!      integer*4, dimension(:),   allocatable :: wavemaker_index_u1_in_z    ! points to index of wm maker u points to the right (size = nr wave maker water level points) 
 !     integer*4, dimension(:),   allocatable :: wavemaker_index_v0_in_z    ! points to index of wm maker v points below (size = nr wave maker water level points) 
!      integer*4, dimension(:),   allocatable :: wavemaker_index_v1_in_z    ! points to index of wm maker v points above (size = nr wave maker water level points) 
!      integer*4, dimension(:),   allocatable :: wavemaker_index_nm_in_z    ! points to nm index (size=nr wave maker water level points) 
!      real*4,    dimension(:),   allocatable :: wavemaker_zwav      
      !
      ! WMU
      !
!      integer*4, dimension(:),   allocatable :: wavemaker_index_nm_in_u
!      integer*4, dimension(:),   allocatable :: wavemaker_index_nmi_in_u
!      integer*4, dimension(:),   allocatable :: wavemaker_index_nmb_in_u
!      integer*4, dimension(:),   allocatable :: wavemaker_index_z_in_u
!      integer*4, dimension(:),   allocatable :: wavemaker_idir_u
!      real*4,    dimension(:),   allocatable :: wavemaker_qxm      
!      real*4,    dimension(:),   allocatable :: wavemaker_angfac_u
      !
      ! WMV
      !
!      integer*4, dimension(:),   allocatable :: wavemaker_index_nm_in_v
!      integer*4, dimension(:),   allocatable :: wavemaker_index_nmi_in_v
!      integer*4, dimension(:),   allocatable :: wavemaker_index_nmb_in_v
!      integer*4, dimension(:),   allocatable :: wavemaker_index_z_in_v
!      integer*4, dimension(:),   allocatable :: wavemaker_idir_v
!      real*4,    dimension(:),   allocatable :: wavemaker_qym
!      real*4,    dimension(:),   allocatable :: wavemaker_angfac_v
      !
      ! General grid
      !
!      real*4, dimension(:,:), allocatable, target :: xg, yg, xz, yz
      !
      ! Sub-grid
      !
      integer                             :: subgrid_nbins
      !
      real*4, dimension(:),   allocatable, target :: subgrid_z_zmin
      real*4, dimension(:),   allocatable :: subgrid_z_zmax
      real*4, dimension(:),   allocatable :: subgrid_z_zmean
      real*4, dimension(:),   allocatable :: subgrid_z_volmax
      real*4, dimension(:,:), allocatable :: subgrid_z_dep
      !
      real*4, dimension(:),   allocatable :: subgrid_uv_zmin
      real*4, dimension(:),   allocatable :: subgrid_uv_zmax
      real*4, dimension(:,:), allocatable :: subgrid_uv_hrep
      real*4, dimension(:,:), allocatable :: subgrid_uv_navg
      real*4, dimension(:),   allocatable :: subgrid_uv_hrep_zmax
      real*4, dimension(:),   allocatable :: subgrid_uv_navg_zmax
      !
      ! Dynamic data on the grid
      !
      real*4, dimension(:),   allocatable :: zsmax
      real*4, dimension(:),   allocatable :: vmax
      real*4, dimension(:),   allocatable :: qmax
      real*4, dimension(:),   allocatable, target :: zs
      real*4, dimension(:),   allocatable :: zsm
      real*4, dimension(:),   allocatable :: q
      real*4, dimension(:),   allocatable :: q0
      real*4, dimension(:),   allocatable :: uv
      real*4, dimension(:),   allocatable :: uv0
      real*4, dimension(:),   allocatable :: z_volume
      real*4, dimension(:),   allocatable :: twet
      real*4, dimension(:),   allocatable :: tsunami_arrival_time
      !
      real*4, dimension(:),   allocatable :: tauwu
      real*4, dimension(:),   allocatable :: tauwv
      real*4, dimension(:),   allocatable :: patm
      real*4, dimension(:),   allocatable :: prcp
      real*4, dimension(:),   allocatable :: cumprcp
      real*4, dimension(:),   allocatable :: netprcp
      real*4, dimension(:),   allocatable :: cumprcpt
      real*4, dimension(:),   allocatable :: cuminf
      real*4, dimension(:),   allocatable :: tauwu0
      real*4, dimension(:),   allocatable :: tauwu1
      real*4, dimension(:),   allocatable :: tauwv0
      real*4, dimension(:),   allocatable :: tauwv1
      real*4, dimension(:),   allocatable :: patm0
      real*4, dimension(:),   allocatable :: patm1
      real*4, dimension(:),   allocatable :: prcp0
      real*4, dimension(:),   allocatable :: prcp1 
      !
      real*4, dimension(:),   allocatable :: windu
      real*4, dimension(:),   allocatable :: windv
      real*4, dimension(:),   allocatable :: windu0
      real*4, dimension(:),   allocatable :: windv0
      real*4, dimension(:),   allocatable :: windu1
      real*4, dimension(:),   allocatable :: windv1
      real*4, dimension(:),   allocatable :: windmax
      !
      real*4, dimension(:),   allocatable :: fwx
      real*4, dimension(:),   allocatable :: fwy
      real*4, dimension(:),   allocatable :: hm0
      real*4, dimension(:),   allocatable :: hm0_ig
      real*4, dimension(:),   allocatable :: fwuv
      real*4, dimension(:),   allocatable :: mean_wave_direction
      real*4, dimension(:),   allocatable :: wave_directional_spreading
!      real*4, dimension(:),   allocatable :: tauwavv
      !
      !!!
      !!! Boundary data
      !!!
      integer ntb
      integer nbnd
      integer ngbnd, itbndlast, ntbnd
      integer nwbnd, itwbndlast, ntwbnd
      !
      ! Grid boundary points
      !
      integer*4,          dimension(:),     allocatable :: ind1_bnd_gbp
      integer*4,          dimension(:),     allocatable :: ind2_bnd_gbp
      integer*4,          dimension(:),     allocatable :: ind1_cst_gbp
      integer*4,          dimension(:),     allocatable :: ind2_cst_gbp
      real*4,             dimension(:),     allocatable :: fac_bnd_gbp
      real*4,             dimension(:),     allocatable :: fac_cst_gbp
      real*4,             dimension(:),     allocatable :: zsb
      real*4,             dimension(:),     allocatable :: zsb0
      real*4,             dimension(:),     allocatable :: patmb
      integer*4,          dimension(:),     allocatable :: ibkcuv2
      integer*1,          dimension(:),     allocatable :: ibndtype
      !
      ! Water level boundary points
      !
      real*4, dimension(:),     allocatable :: x_bnd
      real*4, dimension(:),     allocatable :: y_bnd
      real*4, dimension(:),     allocatable :: t_bnd
      real*4, dimension(:,:),   allocatable :: zs_bnd
      real*4, dimension(:,:),   allocatable :: zsi_bnd
      real*4, dimension(:),     allocatable, target :: zst_bnd
      real*4, dimension(:),     allocatable :: zsit_bnd
      !
      ! Wave boundary points
      !
      real*4, dimension(:),     allocatable :: x_bwv
      real*4, dimension(:),     allocatable :: y_bwv
      real*4, dimension(:),     allocatable :: t_bwv
      real*4, dimension(:,:),   allocatable :: hs_bwv
      real*4, dimension(:,:),   allocatable :: tp_bwv
      real*4, dimension(:,:),   allocatable :: wd_bwv
      real*4, dimension(:),     allocatable :: hst_bwv
      real*4, dimension(:),     allocatable :: tpt_bwv
      real*4, dimension(:),     allocatable :: l0t_bwv
      real*4, dimension(:),     allocatable :: wdt_bwv
      !
      ! cst points
      !
      integer                               :: ncst
      real*4, dimension(:),     allocatable :: x_cst
      real*4, dimension(:),     allocatable :: y_cst
      real*4, dimension(:),     allocatable :: slope_cst
      real*4, dimension(:),     allocatable :: angle_cst
      real*4, dimension(:),     allocatable :: zsetup_cst
      real*4, dimension(:),     allocatable :: zig_cst
      real*4, dimension(:),     allocatable :: fac_bwv_cst
      integer*4, dimension(:),  allocatable :: ind1_bwv_cst
      integer*4, dimension(:),  allocatable :: ind2_bwv_cst
      !
      ! IG frequencies
      !
      real*4, dimension(:),     allocatable :: freqig
      real*4, dimension(:),     allocatable :: costig
      real*4, dimension(:),     allocatable :: phiig
      real*4, dimension(:),     allocatable :: dphiig
      real*4                                :: dfreqig
      !!!
      !!! Meteo data
      !!!
      integer                                 :: ntwnd
      integer                                 :: ntprcp
      integer                                 :: itwndlast
      integer                                 :: itprcplast
      integer                                 :: spw_nt
      integer                                 :: spw_nrows
      integer                                 :: spw_ncols
      integer                                 :: spw_nquant
      real*4                                  :: spw_radius
      real*4, dimension(:),     allocatable   :: spw_times
      real*4, dimension(:),     allocatable   :: spw_xe
      real*4, dimension(:),     allocatable   :: spw_ye
      real*4, dimension(:,:,:),   allocatable :: spw_vmag
      real*4, dimension(:,:,:),   allocatable :: spw_vdir
      real*4, dimension(:,:,:),   allocatable :: spw_wu
      real*4, dimension(:,:,:),   allocatable :: spw_wv
      real*4, dimension(:,:,:),   allocatable :: spw_pdrp
      real*4, dimension(:,:,:),   allocatable :: spw_prcp
      real*4, dimension(:,:),     allocatable :: spw_wu01
      real*4, dimension(:,:),     allocatable :: spw_wv01
      real*4, dimension(:,:),     allocatable :: spw_pdrp01
      real*4, dimension(:,:),     allocatable :: spw_prcp01
      !
      integer                                 :: amuv_nt
      integer                                 :: amuv_nrows
      integer                                 :: amuv_ncols
      real*4, dimension(:),     allocatable   :: amuv_times
      real*4, dimension(:,:,:),   allocatable :: amuv_wu
      real*4, dimension(:,:,:),   allocatable :: amuv_wv
      real*4, dimension(:,:),     allocatable :: amuv_wu01
      real*4, dimension(:,:),     allocatable :: amuv_wv01
      integer                                 :: amuv_nquant
      real*4                                  :: amuv_x_llcorner
      real*4                                  :: amuv_y_llcorner
      real*4                                  :: amuv_dx
      real*4                                  :: amuv_dy
      !
      integer                                 :: amp_nt
      integer                                 :: amp_nrows
      integer                                 :: amp_ncols
      real*4, dimension(:),     allocatable   :: amp_times
      real*4, dimension(:,:,:),   allocatable :: amp_patm
      real*4, dimension(:,:),   allocatable   :: amp_patm01
      integer                                 :: amp_nquant
      real*4                                  :: amp_x_llcorner
      real*4                                  :: amp_y_llcorner
      real*4                                  :: amp_dx
      real*4                                  :: amp_dy
      !
      integer                                 :: ampr_nt
      integer                                 :: ampr_nrows
      integer                                 :: ampr_ncols
      real*4, dimension(:),     allocatable   :: ampr_times
      real*4, dimension(:,:,:),   allocatable :: ampr_pr
      real*4, dimension(:,:),   allocatable   :: ampr_pr01
      integer                                 :: ampr_nquant
      real*4                                  :: ampr_x_llcorner
      real*4                                  :: ampr_y_llcorner
      real*4                                  :: ampr_dx
      real*4                                  :: ampr_dy
      !
      real*4, dimension(:),     allocatable   :: twnd
      real*4, dimension(:),     allocatable   :: wndmag
      real*4, dimension(:),     allocatable   :: wnddir
      real*4, dimension(:),     allocatable   :: tprcpt
      real*4, dimension(:),     allocatable   :: tprcpv
      real*4, dimension(1000)                 :: cdval
      !
      real*4 meteo_t0
      real*4 meteo_t1
      !
      real*4 twindupd
      real*4 dradspw
      real*4 dphispw
      !!!
      !!! Observation points
      !!!
      integer                               :: nobs
      integer*4, dimension(:),  allocatable :: nmindobs, mindobs, nindobs, idobs, nmwindobs
      real*4, dimension(:),     allocatable :: xobs
      real*4, dimension(:),     allocatable :: yobs
      real*4, dimension(:),     allocatable :: zobs
      real*4, dimension(:),     allocatable :: hobs
      real*4, dimension(:),     allocatable :: xgobs
      real*4, dimension(:),     allocatable :: ygobs
      real*4, dimension(:),     allocatable :: zbobs
      real*4, dimension(:,:),   allocatable :: wobs
      character*256, dimension(:), allocatable :: nameobs
      !!!
      !!! Discharges and drainage
      !!!
      integer                               :: nsrc
      integer                               :: ndrn
      integer                               :: nsrcdrn
      integer                               :: ntsrc
      integer                               :: itsrclast
      real*4, dimension(:),     allocatable, target :: tsrc
      real*4, dimension(:,:),   allocatable, target :: qsrc
      real*4, dimension(:),     allocatable :: qtsrc
      integer*4, dimension(:),  allocatable :: nmindsrc
      integer*1, dimension(:),  allocatable :: drainage_type
      real*4, dimension(:,:),   allocatable :: drainage_params
      real*4, dimension(:),     allocatable, target :: xsrc
      real*4, dimension(:),     allocatable, target :: ysrc
      !!!
      !!! Structures
      !!!
      integer                                  :: nrstructures
      integer,   dimension(:),     allocatable :: structure_uv_index
      integer*1, dimension(:),     allocatable :: structure_type
      real*4,    dimension(:,:),   allocatable :: structure_parameters
      real*4,    dimension(:),     allocatable :: structure_length
!      real*4,    dimension(:),     allocatable :: struc_x
!      real*4,    dimension(:),     allocatable :: struc_y
!      real*4,    dimension(:),     allocatable :: struc_height
      !!!
      !!! Cross-sections
      !!!
      integer                                  :: nrcrosssections
      integer, dimension(:),     allocatable   :: crs_nr
      integer, dimension(:,:),   allocatable   :: crs_uv_index
      integer, dimension(:,:),   allocatable   :: crs_idir
      character*256, dimension(:), allocatable :: namecrs
      !
      real*4 :: waveage
      !
      ! LGX Cd table
      !
      real*4, dimension(50, 20) :: cdlgx = reshape((/ 0.957, 1.241, 1.473, 1.682, 1.875, 2.057, 2.235, 2.414, 2.588, 2.767, 2.933, 3.115, 3.286, 3.421, 3.441, 3.462, 3.462, 3.472, 3.462, 3.451, 3.441, 3.421, 3.402, 3.372, 3.343, 3.314, 3.286, 3.258, 3.221, 3.185, 3.159, 3.124, 3.081, 3.047, 3.014, 2.973, 2.941, 2.902, 2.871, 2.833, 2.796, 2.767, 2.731, 2.696, 2.661, 2.627, 2.594, 2.562, 2.530, 2.498, &
                                              1.000, 1.320, 1.586, 1.832, 2.067, 2.294, 2.523, 2.745, 2.973, 3.194, 3.431, 3.662, 3.869, 3.857, 3.834, 3.798, 3.764, 3.718, 3.662, 3.608, 3.555, 3.502, 3.441, 3.382, 3.324, 3.267, 3.212, 3.150, 3.098, 3.039, 2.981, 2.933, 2.879, 2.826, 2.774, 2.724, 2.675, 2.627, 2.581, 2.536, 2.486, 2.443, 2.402, 2.361, 2.322, 2.283, 2.246, 2.204, 2.168, 2.134, &
                                              0.992, 1.325, 1.609, 1.871, 2.124, 2.379, 2.627, 2.879, 3.132, 3.392, 3.662, 3.942, 4.068, 4.004, 3.942, 3.881, 3.798, 3.729, 3.651, 3.576, 3.492, 3.411, 3.333, 3.258, 3.185, 3.115, 3.039, 2.973, 2.902, 2.841, 2.774, 2.710, 2.648, 2.594, 2.536, 2.480, 2.425, 2.373, 2.322, 2.278, 2.230, 2.184, 2.139, 2.095, 2.053, 2.016, 1.976, 1.938, 1.900, 1.863, &
                                              0.960, 1.295, 1.583, 1.851, 2.114, 2.379, 2.641, 2.910, 3.176, 3.462, 3.752, 4.055, 4.120, 4.042, 3.942, 3.845, 3.752, 3.651, 3.555, 3.462, 3.362, 3.277, 3.185, 3.098, 3.014, 2.933, 2.856, 2.774, 2.703, 2.634, 2.562, 2.498, 2.437, 2.373, 2.316, 2.262, 2.204, 2.154, 2.100, 2.053, 2.007, 1.959, 1.917, 1.875, 1.832, 1.793, 1.756, 1.717, 1.682, 1.648, &
                                              0.915, 1.243, 1.528, 1.797, 2.062, 2.322, 2.594, 2.863, 3.150, 3.441, 3.741, 4.055, 4.120, 4.004, 3.881, 3.775, 3.651, 3.544, 3.431, 3.324, 3.212, 3.115, 3.022, 2.925, 2.833, 2.752, 2.668, 2.588, 2.511, 2.437, 2.373, 2.305, 2.241, 2.179, 2.119, 2.062, 2.007, 1.955, 1.908, 1.859, 1.812, 1.767, 1.724, 1.685, 1.645, 1.606, 1.568, 1.534, 1.499, 1.465, &
                                              0.866, 1.180, 1.454, 1.717, 1.976, 2.235, 2.498, 2.774, 3.055, 3.353, 3.662, 3.979, 4.055, 3.930, 3.798, 3.662, 3.534, 3.411, 3.286, 3.168, 3.055, 2.957, 2.856, 2.752, 2.661, 2.575, 2.492, 2.408, 2.333, 2.262, 2.189, 2.124, 2.062, 1.998, 1.942, 1.888, 1.836, 1.782, 1.734, 1.689, 1.641, 1.599, 1.559, 1.520, 1.482, 1.446, 1.408, 1.375, 1.342, 1.311, &
                                              0.807, 1.105, 1.367, 1.622, 1.871, 2.124, 2.384, 2.648, 2.925, 3.212, 3.523, 3.845, 3.979, 3.834, 3.685, 3.534, 3.402, 3.267, 3.141, 3.022, 2.902, 2.796, 2.689, 2.588, 2.498, 2.408, 2.322, 2.246, 2.168, 2.095, 2.025, 1.963, 1.900, 1.840, 1.782, 1.731, 1.678, 1.628, 1.583, 1.537, 1.493, 1.454, 1.413, 1.375, 1.340, 1.304, 1.272, 1.239, 1.209, 1.178, &
                                              0.747, 1.029, 1.274, 1.514, 1.749, 1.994, 2.241, 2.498, 2.767, 3.047, 3.343, 3.662, 3.881, 3.718, 3.555, 3.411, 3.258, 3.124, 2.989, 2.871, 2.752, 2.641, 2.536, 2.437, 2.344, 2.251, 2.168, 2.090, 2.016, 1.946, 1.879, 1.812, 1.753, 1.696, 1.641, 1.590, 1.540, 1.490, 1.446, 1.403, 1.362, 1.323, 1.285, 1.249, 1.215, 1.182, 1.150, 1.120, 1.092, 1.064, &
                                              0.684, 0.946, 1.178, 1.398, 1.622, 1.847, 2.081, 2.322, 2.581, 2.848, 3.132, 3.431, 3.752, 3.597, 3.431, 3.277, 3.124, 2.981, 2.848, 2.724, 2.601, 2.492, 2.384, 2.289, 2.194, 2.109, 2.025, 1.950, 1.875, 1.805, 1.742, 1.678, 1.622, 1.565, 1.511, 1.462, 1.413, 1.370, 1.325, 1.285, 1.245, 1.209, 1.172, 1.139, 1.107, 1.075, 1.045, 1.017, 0.989, 0.963, &
                                              0.682, 0.862, 1.075, 1.281, 1.488, 1.699, 1.917, 2.144, 2.379, 2.634, 2.894, 3.176, 3.482, 3.482, 3.305, 3.141, 2.989, 2.841, 2.710, 2.581, 2.461, 2.350, 2.246, 2.154, 2.062, 1.976, 1.896, 1.820, 1.749, 1.678, 1.615, 1.556, 1.499, 1.446, 1.395, 1.347, 1.302, 1.258, 1.217, 1.180, 1.143, 1.107, 1.073, 1.040, 1.009, 0.980, 0.951, 0.926, 0.900, 0.875, &
                                              0.682, 0.852, 0.984, 1.162, 1.350, 1.543, 1.742, 1.950, 2.168, 2.402, 2.648, 2.910, 3.194, 3.353, 3.176, 3.014, 2.856, 2.710, 2.575, 2.449, 2.333, 2.220, 2.119, 2.025, 1.933, 1.851, 1.771, 1.699, 1.628, 1.565, 1.502, 1.446, 1.390, 1.340, 1.290, 1.245, 1.202, 1.160, 1.122, 1.083, 1.049, 1.016, 0.989, 0.968, 0.947, 0.929, 0.909, 0.892, 0.875, 0.858, &
                                              0.682, 0.852, 0.984, 1.100, 1.213, 1.387, 1.565, 1.756, 1.955, 2.163, 2.390, 2.634, 2.886, 3.168, 3.055, 2.886, 2.724, 2.581, 2.449, 2.322, 2.204, 2.100, 1.998, 1.904, 1.816, 1.734, 1.658, 1.586, 1.520, 1.457, 1.400, 1.345, 1.292, 1.243, 1.196, 1.164, 1.135, 1.107, 1.082, 1.057, 1.032, 1.009, 0.989, 0.968, 0.947, 0.929, 0.909, 0.892, 0.875, 0.858, &
                                              0.682, 0.852, 0.984, 1.100, 1.205, 1.302, 1.395, 1.562, 1.742, 1.929, 2.134, 2.344, 2.581, 2.833, 2.933, 2.759, 2.601, 2.461, 2.327, 2.204, 2.086, 1.981, 1.884, 1.793, 1.706, 1.628, 1.556, 1.485, 1.421, 1.372, 1.332, 1.295, 1.258, 1.226, 1.194, 1.164, 1.135, 1.107, 1.082, 1.057, 1.032, 1.009, 0.989, 0.968, 0.947, 0.929, 0.909, 0.892, 0.875, 0.858, &
                                              0.682, 0.852, 0.984, 1.100, 1.205, 1.302, 1.395, 1.485, 1.574, 1.696, 1.875, 2.062, 2.267, 2.492, 2.731, 2.641, 2.486, 2.344, 2.209, 2.090, 1.976, 1.875, 1.778, 1.689, 1.606, 1.556, 1.505, 1.459, 1.413, 1.372, 1.332, 1.295, 1.258, 1.226, 1.194, 1.164, 1.135, 1.107, 1.082, 1.057, 1.032, 1.009, 0.989, 0.968, 0.947, 0.929, 0.909, 0.892, 0.875, 0.858, &
                                              0.682, 0.852, 0.984, 1.100, 1.205, 1.302, 1.395, 1.485, 1.574, 1.658, 1.742, 1.828, 1.968, 2.158, 2.367, 2.530, 2.373, 2.230, 2.105, 1.985, 1.871, 1.797, 1.731, 1.668, 1.612, 1.556, 1.505, 1.459, 1.413, 1.372, 1.332, 1.295, 1.258, 1.226, 1.194, 1.164, 1.135, 1.107, 1.082, 1.057, 1.032, 1.009, 0.989, 0.968, 0.947, 0.929, 0.909, 0.892, 0.875, 0.858, &
                                              0.682, 0.852, 0.984, 1.100, 1.205, 1.302, 1.395, 1.485, 1.574, 1.658, 1.742, 1.828, 1.908, 1.989, 2.071, 2.204, 2.267, 2.129, 2.034, 1.946, 1.871, 1.797, 1.731, 1.668, 1.612, 1.556, 1.505, 1.459, 1.413, 1.372, 1.332, 1.295, 1.258, 1.226, 1.194, 1.164, 1.135, 1.107, 1.082, 1.057, 1.032, 1.009, 0.989, 0.968, 0.947, 0.929, 0.909, 0.892, 0.875, 0.858, &
                                              0.682, 0.852, 0.984, 1.100, 1.205, 1.302, 1.395, 1.485, 1.574, 1.658, 1.742, 1.828, 1.908, 1.989, 2.071, 2.154, 2.225, 2.124, 2.034, 1.946, 1.871, 1.797, 1.731, 1.668, 1.612, 1.556, 1.505, 1.459, 1.413, 1.372, 1.332, 1.295, 1.258, 1.226, 1.194, 1.164, 1.135, 1.107, 1.082, 1.057, 1.032, 1.009, 0.989, 0.968, 0.947, 0.929, 0.909, 0.892, 0.875, 0.858, &
                                              0.682, 0.852, 0.984, 1.100, 1.205, 1.302, 1.395, 1.485, 1.574, 1.658, 1.742, 1.828, 1.908, 1.989, 2.071, 2.154, 2.225, 2.124, 2.034, 1.946, 1.871, 1.797, 1.731, 1.668, 1.612, 1.556, 1.505, 1.459, 1.413, 1.372, 1.332, 1.295, 1.258, 1.226, 1.194, 1.164, 1.135, 1.107, 1.082, 1.057, 1.032, 1.009, 0.989, 0.968, 0.947, 0.929, 0.909, 0.892, 0.875, 0.858, &
                                              0.682, 0.852, 0.984, 1.100, 1.205, 1.302, 1.395, 1.485, 1.574, 1.658, 1.742, 1.828, 1.908, 1.989, 2.071, 2.154, 2.225, 2.124, 2.034, 1.946, 1.871, 1.797, 1.731, 1.668, 1.612, 1.556, 1.505, 1.459, 1.413, 1.372, 1.332, 1.295, 1.258, 1.226, 1.194, 1.164, 1.135, 1.107, 1.082, 1.057, 1.032, 1.009, 0.989, 0.968, 0.947, 0.929, 0.909, 0.892, 0.875, 0.858, &
                                              0.682, 0.852, 0.984, 1.100, 1.205, 1.302, 1.395, 1.485, 1.574, 1.658, 1.742, 1.828, 1.908, 1.989, 2.071, 2.154, 2.225, 2.124, 2.034, 1.946, 1.871, 1.797, 1.731, 1.668, 1.612, 1.556, 1.505, 1.459, 1.413, 1.372, 1.332, 1.295, 1.258, 1.226, 1.194, 1.164, 1.135, 1.107, 1.082, 1.057, 1.032, 1.009, 0.989, 0.968, 0.947, 0.929, 0.909, 0.892, 0.875, 0.858 /), shape(cdlgx))

   contains


   subroutine initialize_parameters()
      !
      implicit none
      !
      ! Parameters for sfincs.f90
      integer                       :: nt
      integer                       :: itmapout
      integer                       :: itmaxout
      integer                       :: itrstout
      integer                       :: ithisout
      integer                       :: iwndupd
      !
      real*8                       :: t
      real*4                       :: dt
      real*4                       :: tmapout
      real*4                       :: tmaxout
      real*4                       :: trstout
      real*4                       :: thisout
      real*4                       :: maxdepth
      real*4                       :: maxmaxdepth
      real*4                       :: twindupd
      !
      ! I don't think this is used anywhere ...
      !
      real :: tstart, tfinish, tloop2, tloop3, tloopstruc, tloopbnd, tloopsrc, tloopwnd1
      real :: tinput
      !
      t           = t0     ! start time
      dt          = 1.0e-6 ! First time step very small
      dtavg       = 0.0
      maxdepth    = 999.0
      maxmaxdepth = 0.0
      nt          = 0
      itmapout    = 0
      itmaxout    = 0
      itrstout    = 0
      ithisout    = 0
      twindupd    = t0
      !
      tloop2      = 0.0
      tloop3      = 0.0
      tloopstruc  = 0.0
      tloopbnd    = 0.0
      tloopsrc    = 0.0
      tloopwnd1   = 0.0
      !
      !
   end subroutine


   subroutine finalize_parameters()
   !
   implicit none
   !
    ! memory cleanup
!    if(allocated(indices)) deallocate(indices)
    if(allocated(nmindbnd)) deallocate(nmindbnd)
!    if(allocated(index_v_m)) deallocate(index_v_m)
!    if(allocated(index_v_n)) deallocate(index_v_n)
!    if(allocated(index_v_nm)) deallocate(index_v_nm)
!    if(allocated(index_v_nmu)) deallocate(index_v_nmu)
!    if(allocated(index_v_num)) deallocate(index_v_num)
!    if(allocated(index_v_nmd)) deallocate(index_v_nmd)
!    if(allocated(index_v_ndm)) deallocate(index_v_ndm)
!    if(allocated(index_g_nm)) deallocate(index_g_nm)
    !
    ! Flags
    !
!    if(allocated(kcu)) deallocate(kcu)
!    if(allocated(kcv)) deallocate(kcv)
!    if(allocated(kcs)) deallocate(kcs)
!    if(allocated(kcs_incl3)) deallocate(kcs_incl3)
    !
    if(allocated(zb)) deallocate(zb)
!    if(allocated(zbu)) deallocate(zbu)
!    if(allocated(zbv)) deallocate(zbv)
!    if(allocated(zbumx)) deallocate(zbumx)
!    if(allocated(zbvmx)) deallocate(zbvmx)
!    if(allocated(gn2u)) deallocate(gn2u)
!    if(allocated(gn2v)) deallocate(gn2v)
    if(allocated(qinfmap)) deallocate(qinfmap)
    if(allocated(qinffield)) deallocate(qinffield)
    if(allocated(ksfield)) deallocate(ksfield)
    if(allocated(scs_Se)) deallocate(scs_Se)

    !
    ! Boundary velocity points
    !
!    if(allocated(nmkcu2)) deallocate(nmkcu2)
!    if(allocated(nmbkcu2)) deallocate(nmbkcu2)
!    if(allocated(nmikcu2)) deallocate(nmikcu2)
!    if(allocated(ibudir)) deallocate(ibudir)
!    if(allocated(nmkcv2)) deallocate(nmkcv2)
!    if(allocated(nmbkcv2)) deallocate(nmbkcv2)
!    if(allocated(nmikcv2)) deallocate(nmikcv2)
!    if(allocated(ibvdir)) deallocate(ibvdir)
!    if(allocated(ibkcu2)) deallocate(ibkcu2)
!    if(allocated(ibkcv2)) deallocate(ibkcv2)
!    if(allocated(dirfac_kcu2)) deallocate(dirfac_kcu2)
!    if(allocated(dirfac_kcv2)) deallocate(dirfac_kcv2)
    !
    ! Mean velocity at boundaries
    !
    if(allocated(uvmean)) deallocate(uvmean)
    !
    ! General grid
    !
!    if(allocated(xg)) deallocate(xg)
!    if(allocated(yg)) deallocate(yg)
!    if(allocated(xz)) deallocate(xz)
!    if(allocated(yz)) deallocate(yz)
    !
    ! Sub-grid
    !
    if(allocated(subgrid_z_zmin)) deallocate(subgrid_z_zmin)
    if(allocated(subgrid_z_zmax)) deallocate(subgrid_z_zmax)
    if(allocated(subgrid_z_dep)) deallocate(subgrid_z_dep)
    if(allocated(subgrid_z_volmax)) deallocate(subgrid_z_volmax)
!    if(allocated(subgrid_u_zmin)) deallocate(subgrid_u_zmin)
!    if(allocated(subgrid_u_zmax)) deallocate(subgrid_u_zmax)
!    if(allocated(subgrid_u_hrep)) deallocate(subgrid_u_hrep)
!    if(allocated(subgrid_u_navg)) deallocate(subgrid_u_navg)
!    if(allocated(subgrid_u_dhdz)) deallocate(subgrid_u_dhdz)
!    if(allocated(subgrid_v_zmin)) deallocate(subgrid_v_zmin)
!    if(allocated(subgrid_v_zmax)) deallocate(subgrid_v_zmax)
!    if(allocated(subgrid_v_hrep)) deallocate(subgrid_v_hrep)
!    if(allocated(subgrid_v_dhdz)) deallocate(subgrid_v_dhdz)
!    if(allocated(subgrid_v_navg)) deallocate(subgrid_v_navg)
    !!!
    !!! Dynamic data on the grid
    !!!
    if(allocated(zsmax)) deallocate(zsmax)
    if(allocated(vmax)) deallocate(vmax)
    if(allocated(qmax)) deallocate(qmax)
    if(allocated(zs)) deallocate(zs)
    !if(allocated(z_volume)) deallocate(z_volume) > this one seems to cause an error, not sure why
    if(allocated(q)) deallocate(q)
    if(allocated(q0)) deallocate(q0)
    if(allocated(uv)) deallocate(uv)
    if(allocated(uv0)) deallocate(uv0)
    if(allocated(twet)) deallocate(twet)
    !
!    if(allocated(huu)) deallocate(huu)
!    if(allocated(hvv)) deallocate(hvv)
!    if(allocated(huu0)) deallocate(huu0)
!    if(allocated(hvv0)) deallocate(hvv0)
    !
    if(allocated(tauwu)) deallocate(tauwu)
    if(allocated(tauwv)) deallocate(tauwv)
    if(allocated(patm)) deallocate(patm)
    if(allocated(prcp)) deallocate(prcp)
    if(allocated(cumprcp)) deallocate(cumprcp)
    if(allocated(cumprcpt)) deallocate(cumprcpt)
    if(allocated(tauwu0)) deallocate(tauwu0)
    if(allocated(tauwu1)) deallocate(tauwu1)
    if(allocated(tauwv0)) deallocate(tauwv0)
    if(allocated(tauwv1)) deallocate(tauwv1)
    if(allocated(patm0)) deallocate(patm0)
    if(allocated(patm1)) deallocate(patm1)
    if(allocated(prcp0)) deallocate(prcp0)
    if(allocated(prcp1)) deallocate(prcp1)
    !
    if(allocated(kfuv)) deallocate(kfuv)
    !
    ! Grid boundary points
    !
    if(allocated(ind1_bnd_gbp)) deallocate(ind1_bnd_gbp)
    if(allocated(ind2_bnd_gbp)) deallocate(ind2_bnd_gbp)
    if(allocated(ind1_cst_gbp)) deallocate(ind1_cst_gbp)
    if(allocated(ind2_cst_gbp)) deallocate(ind2_cst_gbp)
    if(allocated(fac_bnd_gbp))  deallocate(fac_bnd_gbp)
    if(allocated(fac_cst_gbp))  deallocate(fac_cst_gbp)
    if(allocated(zsb))          deallocate(zsb)
    if(allocated(zsb0))         deallocate(zsb0)
    if(allocated(ibndtype))     deallocate(ibndtype)
    !
    ! Water level boundary points
    !
    if(allocated(x_bnd)) deallocate(x_bnd)
    if(allocated(y_bnd)) deallocate(y_bnd)
    if(allocated(t_bnd)) deallocate(t_bnd)
    if(allocated(zs_bnd)) deallocate(zs_bnd)
    if(allocated(zsi_bnd)) deallocate(zsi_bnd)
    if(allocated(zst_bnd)) deallocate(zst_bnd)
    if(allocated(zsit_bnd)) deallocate(zsit_bnd)
    !
    ! Wave boundary points
    !
    if(allocated(x_bwv)) deallocate(x_bwv)
    if(allocated(y_bwv)) deallocate(y_bwv)
    if(allocated(t_bwv)) deallocate(t_bwv)
    if(allocated(hs_bwv)) deallocate(hs_bwv)
    if(allocated(tp_bwv)) deallocate(tp_bwv)
    if(allocated(wd_bwv)) deallocate(wd_bwv)
    if(allocated(hst_bwv)) deallocate(hst_bwv)
    if(allocated(tpt_bwv)) deallocate(tpt_bwv)
    if(allocated(l0t_bwv)) deallocate(l0t_bwv)
    if(allocated(wdt_bwv)) deallocate(wdt_bwv)
    !
    ! cst points
    !
    if(allocated(x_cst)) deallocate(x_cst)
    if(allocated(y_cst)) deallocate(y_cst)
    if(allocated(slope_cst)) deallocate(slope_cst)
    if(allocated(angle_cst)) deallocate(angle_cst)
    if(allocated(zsetup_cst)) deallocate(zsetup_cst)
    if(allocated(zig_cst)) deallocate(zig_cst)
    if(allocated(fac_bwv_cst)) deallocate(fac_bwv_cst)
    if(allocated(ind1_bwv_cst)) deallocate(ind1_bwv_cst)
    if(allocated(ind2_bwv_cst)) deallocate(ind2_bwv_cst)
    !
    ! IG frequencies
    !
    if(allocated(freqig)) deallocate(freqig)
    if(allocated(costig)) deallocate(costig)
    if(allocated(phiig)) deallocate(phiig)
    if(allocated(dphiig)) deallocate(dphiig)
    !
    if(allocated(spw_times)) deallocate(spw_times)
    if(allocated(spw_xe)) deallocate(spw_xe)
    if(allocated(spw_ye)) deallocate(spw_ye)
    if(allocated(spw_vmag)) deallocate(spw_vmag)
    if(allocated(spw_vdir)) deallocate(spw_vdir)
    if(allocated(spw_wu)) deallocate(spw_wu)
    if(allocated(spw_wv)) deallocate(spw_wv)
    if(allocated(spw_pdrp)) deallocate(spw_pdrp)
    if(allocated(spw_prcp)) deallocate(spw_prcp)
    if(allocated(spw_wu01)) deallocate(spw_wu01)
    if(allocated(spw_wv01)) deallocate(spw_wv01)
    if(allocated(spw_pdrp01)) deallocate(spw_pdrp01)
    if(allocated(spw_prcp01)) deallocate(spw_prcp01)
    !
    if(allocated(amuv_times)) deallocate(amuv_times)
    if(allocated(amuv_wu)) deallocate(amuv_wu)
    if(allocated(amuv_wv)) deallocate(amuv_wv)
    if(allocated(amuv_wu01)) deallocate(amuv_wu01)
    if(allocated(amuv_wv01)) deallocate(amuv_wv01)
    !
    if(allocated(amp_times)) deallocate(amp_times)
    if(allocated(amp_patm)) deallocate(amp_patm)
    if(allocated(amp_patm01)) deallocate(amp_patm01)
    !
    if(allocated(ampr_times)) deallocate(ampr_times)
    if(allocated(ampr_pr)) deallocate(ampr_pr)
    if(allocated(ampr_pr01)) deallocate(ampr_pr01)
    !
    if(allocated(twnd)) deallocate(twnd)
    if(allocated(wndmag)) deallocate(wndmag)
    if(allocated(wnddir)) deallocate(wnddir)
    if(allocated(tprcpt)) deallocate(tprcpt)
    if(allocated(tprcpv)) deallocate(tprcpv)
    !!!
    !!! Obs Points
    !!!
    if(allocated(nmindobs)) deallocate(nmindobs)
    if(allocated(mindobs)) deallocate(mindobs)
    if(allocated(nindobs)) deallocate(nindobs)
    if(allocated(idobs)) deallocate(idobs)
    if(allocated(xobs)) deallocate(xobs)
    if(allocated(yobs)) deallocate(yobs)
    if(allocated(zobs)) deallocate(zobs)
    if(allocated(hobs)) deallocate(hobs)
    if(allocated(ygobs)) deallocate(xgobs)
    if(allocated(ygobs)) deallocate(ygobs)
    if(allocated(zbobs)) deallocate(zbobs)
    !!!
    !!! Discharges
    !!!
    if(allocated(tsrc)) deallocate(tsrc)
    if(allocated(qsrc)) deallocate(qsrc)
    if(allocated(qtsrc)) deallocate(qtsrc)
    if(allocated(nmindsrc)) deallocate(nmindsrc)
    !!!
    !!! Structures
    !!!
!    if(allocated(struc_uv)) deallocate(struc_uv)
!    if(allocated(struc_nm)) deallocate(struc_nm)
!    if(allocated(struc_type)) deallocate(struc_type)
!    if(allocated(struc_pars)) deallocate(struc_pars)
    !
    if(allocated(cd_wnd)) deallocate(cd_wnd)
    if(allocated(cd_val)) deallocate(cd_val)
    !
    end subroutine

end module
