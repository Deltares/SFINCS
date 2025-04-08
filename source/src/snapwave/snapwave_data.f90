module snapwave_data
   !
   !!! Revision data
   character*256 :: build_revision, build_date
   !
   real*8,  dimension(:),       allocatable    :: x, y                    ! x,y coordinates of grid (double precision)
   real*4,  dimension(:),       allocatable    :: xs, ys                  ! x,y coordinates of grid (single precision)
   real*8                                      :: xmn, ymn                
   integer*1, dimension(:),     allocatable    :: msk                     ! mask array (0=land, 1=inner point, 2=boundary)
   real*4,  dimension(:),       allocatable    :: zb                      ! bed level
   real*4,  dimension(:),       allocatable    :: depth                   ! water depth
   real*4,  dimension(:),       allocatable    :: dhdx, dhdy              ! water depth gradients
   real*4,  dimension(:),       allocatable    :: kwav, nwav              ! wave number, ratio Cg/C
   real*4,  dimension(:),       allocatable    :: kwav_ig, nwav_ig        ! wave number, ratio Cg/C
   real*4,  dimension(:),       allocatable    :: C, Cg                   ! wave celerity, group velocity0
   real*4,  dimension(:),       allocatable    :: C_ig, Cg_ig             ! wave celerity, group velocity0
   real*4,  dimension(:),       allocatable    :: fw                      ! friction coefficient
   real*4,  dimension(:),       allocatable    :: fw_ig                   ! friction coefficient
   real*4,  dimension(:),       allocatable    :: H, H_ig                 ! rms wave height
   real*4,  dimension(:),       allocatable    :: Tp                      ! peak wave period
   real*4,  dimension(:),       allocatable    :: Tp_ig                   ! infragravity peak wave period   
   real*4,  dimension(:),       allocatable    :: Dw,Df                   ! dissipation due to breaking, bed friction
   real*4,  dimension(:),       allocatable    :: Dw_ig,Df_ig             ! dissipation due to breaking, bed friction for IG   
   real*4,  dimension(:),       allocatable    :: F                       ! wave force Dw/C/rho/depth
   real*4,  dimension(:),       allocatable    :: Fx, Fy                  ! wave force Dw/C/rho/depth
   real*4,  dimension(:),       allocatable    :: u10, u10dir             ! wind speed and wind direction   
   real*4                                      :: u10dmean=-999.          ! average wind direction   
   real*4,  dimension(:),       allocatable    :: thetam                  ! mean wave direction
   real*4                                      :: thetamean=-999.         ! mean wave direction
   real*4,  dimension(:),       allocatable    :: buf                     ! buffer for writing output to netcdf
   integer, dimension(:,:),     allocatable    :: kp                      ! surrounding points for each grid point (unstructured)
   logical, dimension(:),       allocatable    :: inner                   ! mask for inner points (not on any boundary)
   integer, dimension(:),       allocatable    :: neumannconnected        ! neumann boundary indices connected to inner points
   integer                                     :: noneumannpts            ! number of neumann boundary points
   real*4,  dimension(:),       allocatable    :: theta                   ! wave angles,sine and cosine of wave angles
   real*4,  dimension(:),       allocatable    :: dist                    ! relative distribution of energy over theta bins
   real*4,  dimension(:),       allocatable    :: theta360                ! wave angles,sine and cosine of wave angles
   integer, dimension(:),       allocatable    :: i360                    ! reference between partial thea grid and full 360 deg theta grid
   real*4,  dimension(:,:),     allocatable    :: windspread360           ! wind input distribution array full 360 deg theta grid
   real*4,  dimension(:,:),     allocatable    :: windspreadfac           ! wind input distribution array   
   integer, dimension(:,:,:),   allocatable    :: prev                    ! two upwind grid points per grid point and wave direction
   integer, dimension(:,:,:),   allocatable    :: prev360                 ! two upwind grid points per grid point and wave direction
   real*4,  dimension(:,:,:),   allocatable    :: w                       ! weights of upwind grid points, 2 per grid point and per wave direction
   real*4,  dimension(:,:,:),   allocatable    :: w360                    ! weights of upwind grid points, 2 per grid point and per wave direction
   real*4,  dimension(:,:),     allocatable    :: ds                      ! distance to interpolated upwind point, per grid point and direction
   real*4,  dimension(:,:),     allocatable    :: ds360                   ! distance to interpolated upwind point, per grid point and direction
   real*4,  dimension(:,:),     allocatable    :: ctheta                  ! refraction speed, per grid point and direction
   real*4,  dimension(:,:),     allocatable    :: ctheta_ig               ! refraction speed, per grid point and direction
   real*4,  dimension(:,:),     allocatable    :: ctheta360               ! refraction speed, per grid point and direction
!   real*4,  dimension(:),       allocatable    :: xn,yn,zn               ! coordinates of nodes of unstructured grid
   real*4,  dimension(:),       allocatable    :: dzdx,dzdy               ! bed slopes at nodes of unstructured grid
   integer, dimension(:,:),     allocatable    :: face_nodes         ! node numbers connected to each cell
   integer, dimension(:,:),     allocatable    :: edge_nodes             ! node numbers connected to each edge
   real*4,  dimension(:),       allocatable    :: bndindx
   real*4,  dimension(:),       allocatable    :: tau
!   real*4,  dimension(:,:),     allocatable    :: Fluxtab
!   real*4,  dimension(:),       allocatable    :: F0tab
   real*4                                      :: Hmax
   real*4,  dimension(:),       allocatable    :: sinhkh
   real*4,  dimension(:),       allocatable    :: Hmx
   real*4,  dimension(:),       allocatable    :: sinhkh_ig
   real*4,  dimension(:),       allocatable    :: Hmx_ig
   real*4,  dimension(:,:),     allocatable    :: ee                      ! directional energy density
   real*4,  dimension(:,:),     allocatable    :: ee_ig                   ! directional infragravity energy density
   !
   real*4,  dimension(:,:),     allocatable    :: aa                      ! directional action density
   real*4,  dimension(:),       allocatable    :: sig                     ! mean frequency
   real*4,  dimension(:,:),     allocatable    :: WsorE                   ! wind input energy
   real*4,  dimension(:,:),     allocatable    :: WsorA                   ! wind input action
   real*4,  dimension(:),       allocatable    :: SwE                     ! directionally integrated wind input energy
   real*4,  dimension(:),       allocatable    :: SwA                     ! directionally integrated wind input wave action   
   !
   real*4,  dimension(:),       allocatable    :: Qb
   real*4,  dimension(:),       allocatable    :: beta
   real*4,  dimension(:),       allocatable    :: srcig
   real*4,  dimension(:),       allocatable    :: alphaig   
   !
   integer*4,  dimension(:),     allocatable    :: index_snapwave_in_quadtree
   integer*4,  dimension(:),     allocatable    :: index_quadtree_in_snapwave
   
   character*256 :: trefstr_iso8601
   character*41  :: treftimefews
   character*15  :: trefstr
   character*15  :: tstartstr
   character*15  :: tstopstr
   !
   ! Boundary conditions (single point, time-varying)
   !
   character*232                               :: snapwave_jonswapfile   ! filename of time-varying wave and water level data
   !
   ! Boundary conditions (space- and time-varying)
   !
   character*256                               :: snapwave_bndfile
   character*256                               :: snapwave_encfile
   character*256                               :: snapwave_bhsfile
   character*256                               :: snapwave_btpfile
   character*256                               :: snapwave_bwdfile
   character*256                               :: snapwave_bdsfile
   character*256                               :: netsnapwavefile   
   !
   integer                                     :: nwbnd                   ! number of support points wave boundary 
   integer                                     :: ntwbnd                  ! number of time points wave boundary 
   real*4                                      :: tpmean_bwv              ! mean tp over boundary points for given time
   real*4                                      :: wdmean_bwv              ! mean wave direction for given time, used to make theta grid
   real*4                                      :: zsmean_bwv              ! mean water level for given time, used to make theta grid
   real*8,  dimension(:),     allocatable      :: x_bwv                   ! x coordinates of boundary points
   real*8,  dimension(:),     allocatable      :: y_bwv                   ! y coordinates of boundary points
   real*4,  dimension(:),     allocatable      :: t_bwv                   ! times (s) of wave boundary conditions
   integer                                     :: n_bndenc                ! number of boundary enclosure points
   real*8,  dimension(:),     allocatable      :: x_bndenc                ! x coordinates of boundary enclosure points
   real*8,  dimension(:),     allocatable      :: y_bndenc                ! y coordinates of boundary enclosure points
   real*4,  dimension(:),     allocatable      :: hst_bwv                 ! wave height at boundary points for given time
   real*4,  dimension(:),     allocatable      :: tpt_bwv                 ! wave period at boundary points for given time
   real*4,  dimension(:),     allocatable      :: wdt_bwv                 ! wave direction at boundary points for given time
   real*4,  dimension(:),     allocatable      :: dst_bwv                 ! directional spreading at boundary points for given time
   real*4,  dimension(:),     allocatable      :: zst_bwv                 ! water level at boundary points for given time
   real*4,  dimension(:,:),   allocatable      :: eet_bwv                 ! directional spectra at boundary points for given time
   real*4,  dimension(:),     allocatable      :: deptht_bwv                ! water depth at boundary points for given time   
   !
   real*4,  dimension(:,:),     allocatable    :: hs_bwv                  ! wave height for all boundary locations and time points
   real*4,  dimension(:,:),     allocatable    :: tp_bwv                  ! wave period for all boundary locations and time points
   real*4,  dimension(:,:),     allocatable    :: wd_bwv                  ! wave direction (nautical deg) for all boundary locations and time points
   real*4,  dimension(:,:),     allocatable    :: ds_bwv                  ! directional spreading (deg) for all boundary locations and time points
   real*4,  dimension(:,:),     allocatable    :: zs_bwv                  ! water level for all boundary locations and time points
   !
   ! IG:
   real*4                                      :: tpmean_bwv_ig           ! mean IG tp over boundary points for given time   
   real*4,  dimension(:),       allocatable    :: hst_bwv_ig              ! IG wave height at boundary points for given time   
   real*4,  dimension(:),       allocatable    :: tpt_bwv_ig              ! IG wave period at boundary points for given time   
   real*4,  dimension(:,:),     allocatable    :: hs_bwv_ig               ! IG wave height for all boundary locations and time points  
   real*4,  dimension(:,:),     allocatable    :: tp_bwv_ig               ! IG wave period for all boundary locations and time points
   real*4,  dimension(:,:),     allocatable    :: eet_bwv_ig              ! directional IG spectra at boundary points for given time   
   real*4                                      :: jonswapgam              ! JONSWAP gamma value for determination offshore spectrum and IG wave conditions using Herbers, default=3.3
   !
   logical                                     :: update_grid_boundary_points
   !
   integer*4                                       :: itwbndlast
   integer*4,          dimension(:),   allocatable :: nmindbnd            ! index of grid point at grid boundary
   integer*4,          dimension(:),   allocatable :: ind1_bwv_cst        ! index to closest wave boundary point
   integer*4,          dimension(:),   allocatable :: ind2_bwv_cst        ! index to second closest wave boundary point
   real*4,             dimension(:),   allocatable :: fac_bwv_cst         ! weight of closest wave boundary point
   integer*4,          dimension(:),   allocatable :: nmindact            ! index of grid point at active grid    
   !   
   integer*4 mmax
   integer*4 nmax
   real*4 dx
   real*4 dy
   real*4 x0
   real*4 y0
   real*4 rotation
   real*4 cosrot
   real*4 sinrot
   !
   ! Timing
   !
   real*8 :: tstart
   real*8 :: tstop
   real*4 :: timestep
   !
   character*232                               :: map_filename
   character*232                               :: his_filename
   ! 
   ! Local input variables
   !
!   integer                                   :: nx,ny           ! number of grid cells in two directions
!   integer                                   :: m,n             ! number of grid points in two directions
   real*4                                    :: dtheta          ! theta grid resolution
   real*4                                    :: sector          ! theta grid sector   
   real*4                                    :: fw0             ! uniform wave friction factor
   real*4                                    :: fw0_ig          ! uniform wave friction factor (ig waves)
   real*4                                    :: Tpini           ! initial condition for the wave period
   real*4                                    :: fwcutoff        ! depth below which to apply space-varying fw
   real*4                                    :: alpha,gamma     ! coefficients in Baldock wave breaking dissipation model
   real*4                                    :: gammax          ! max wave height/water depth ratio
   integer                                   :: baldock_opt     ! option of Baldock wave breaking dissipation model (opt=1 is without gamma&depth, else is including)
   real*4                                    :: baldock_ratio   ! option controlling from what depth wave breaking should take place: (Hk>baldock_ratio*Hmx(k)), default baldock_ratio=0.2 
   ! TODO - TL: bring back baldock_ratio?
 
   real*4                                    :: hmin            ! minimum water depth
   character*256                             :: gridfile        ! name of gridfile (Delft3D .grd format)
   integer                                   :: sferic          ! sferical (1) or cartesian (0) grid
   integer                                   :: niter           ! maximum number of iterations   
   character*232                             :: depfile         ! name of bathymetry file (Delft3D .dep format)
   character*232                             :: fwfile          ! name of bed friction factor file (Delft3D .dep format)
   character*232                             :: upwfile         ! name of upwind neighbors file
   character*232                             :: mskfile         ! name of mask file
   character*232                             :: indfile         ! name of index file
   character*232                             :: obsfile         ! name of observation points file
   integer                                   :: nobs
   real*8,    dimension(:),    allocatable   :: xobs
   real*8,    dimension(:),    allocatable   :: yobs
   integer ,  dimension(:,:),  allocatable   :: irefobs
   integer ,  dimension(:),    allocatable   :: nrefobs
   real*8  ,  dimension(:,:),  allocatable   :: wobs            
   character*32, dimension(:), allocatable   :: nameobs         ! names of observation points 
   real*4,    dimension(:),    allocatable   :: hm0obs
   real*4,    dimension(:),    allocatable   :: hm0xobs
   real*4,    dimension(:),    allocatable   :: hm0yobs
   real*4,    dimension(:),    allocatable   :: tpobs
   real*4,    dimension(:),    allocatable   :: wdobs
   real*4                                    :: dt              ! time step (no limitation)
   real*4                                    :: tol             ! tolerance(m) for boundary points
   integer                                   :: no_nodes        ! number of unstructured grid nodes
   integer                                   :: no_faces        ! number of unstructured grid cells
   integer                                   :: no_edges        ! number of unstructured grid edges
!   integer                                   :: k,j,k1,k2
!   integer                                   :: irec,nrec
   integer                                   :: ntab = 10
   integer                                   :: nHrel
!   real*4                                    :: xmin,xmax       ! minimum and maximum x coordinates
!   real*4                                    :: ymin,ymax       ! minimum and maximum y coordinates
!   real*4                                    :: ddir,alfa
   character*232                             :: Htabname,Dwtabname,Ftabname,Cgtabname,cthetafactabname,hhtabname
   character*232                             :: neumannfile
   real*4, dimension(:), allocatable         :: xb,yb,xneu,yneu
   real*4                                    :: rghlevland       ! Elevation separation as in SFINCS for simple elevation varying roughness
   real*4                                    :: fwratio          ! Above 'rghlevland' elevation of zb, the friction for incident waves is multiplied with value 'fwratio'
   real*4                                    :: fwigratio        ! Above 'rghlevland' elevation of zb, the friction for IG waves is multiplied with value 'fwratio'      
   !
   character*3                               :: outputformat
   integer                                   :: ja_save_each_iter       ! logical to save output after each iteration or not   
   !   
   ! Local constants
   !   
   real*4                                    :: rho             ! water density
   real*4                                    :: pi              ! water density
   real*4                                    :: g               ! acceleration of gravity
   real*4                                    :: t0,t1,t2,t3,t4  ! timers
   integer                                   :: nb
   integer                                   :: np
   integer                                   :: ntheta
   integer                                   :: ntheta360
   !
   ! Wind input constants
   !
   integer                                   :: jadcgdx
   real*4                                    :: c_dispT             
   integer                                   :: mwind
   real*4                                    :: Tini
   real*4                                    :: sigmin
   real*4                                    :: sigmax   
   integer                                   :: wind_opt     ! option of wind growth on (1) or off (0)         
   !
   ! Vegetation parameters
   !
   integer                                      :: vegetation_opt
   character*232                                :: vegmapfile   ! name of vegetation map file (Delft3D .dep format)
   integer                                      :: nveg         ! Number of vegetation species used [-]
   integer                                      :: no_secveg    ! Number of sections used in vertical schematization of vegetation [-]
   integer                                      :: no_secvegmax
   real*4,  dimension(:,:), allocatable         :: veg_ah       ! Height of vertical sections used in vegetation schematization [m wrt zb_ini (zb0)]
   real*4,  dimension(:,:), allocatable         :: veg_Cd       ! Bulk drag coefficient [-]
   real*4,  dimension(:,:), allocatable         :: veg_bstems   ! Width/diameter of individual vegetation stems [m]
   real*4,  dimension(:,:), allocatable         :: veg_Nstems   ! Number of vegetation stems per unit horizontal area [m-2]
   !
   real*4,  dimension(:),   allocatable         :: Dveg
   !
   ! Infragravity parameters
   !
   integer                                   :: ig_opt          ! option of IG wave settings (1 = default = conservative shoaling based dSxx and Baldock breaking)   
   real*4                                    :: alpha_ig,gamma_ig ! coefficients in Baldock wave breaking dissipation model for IG waves
   real*4                                    :: shinc2ig        ! Ratio of how much of the calculated IG wave source term, is subtracted from the incident wave energy (0-1, 0=default)
   real*4                                    :: alphaigfac      ! Multiplication factor for IG shoaling source/sink term, default = 1.0
   real*4                                    :: eeinc2ig        ! ratio of incident wave energy as first estimate of IG wave energy at boundary
   real*4                                    :: Tinc2ig         ! ratio compared to period Tinc to estimate Tig 
   !
   real*4                                    :: baldock_ratio_ig   ! option controlling from what depth wave breaking should take place for IG waves: (Hk>baldock_ratio*Hmx(k)), default baldock_ratio_ig=0.2    
   integer                                   :: igwaves_opt     ! option of IG waves on (1) or off (0)      
   integer                                   :: iterative_srcig_opt ! option whether IG source/sink term should be calculated in the iterative loop again (iterative_srcig = 1, is a bit slower), 
                                                                ! ... or just a priori based on effectively incident wave energy from previous timestep only   
   integer                                   :: herbers_opt     ! Choice whether you want IG Hm0&Tp be calculated by herbers (=1, default), or want to specify user defined values (0> then snapwave_eeinc2ig & snapwave_Tinc2ig are used) 
   integer                                   :: tpig_opt        ! IG wave period option based on Herbers calculated spectrum, only used if herbers_opt = 1. Options are: 1=Tm01 (default), 2=Tpsmooth, 3=Tp, 4=Tm-1,0    
   !
   ! Switches
   logical                                   :: igwaves             ! switch whether include IG or not
   logical                                   :: wind                ! switch whether include wind or not
   logical                                   :: vegetation          ! switch whether include veggie or not      
   logical                                   :: igherbers   
   logical                                   :: iterative_srcig   
   real*4                                    :: crit
   integer                                   :: nr_sweeps
   !
   logical                                   :: restart
   logical                                   :: coupled_to_sfincs
   !
   integer                                   :: nr_quadtree_points
   !
end module snapwave_data
