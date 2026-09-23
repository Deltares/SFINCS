module snapwave_input
   !
   ! Reads the SnapWave settings from sfincs.inp into snapwave_data.
   !
   ! Subroutines:
   !
   !   read_snapwave_input()
   !     Reads all snapwave_* keywords via get_keyword() and derives the
   !     SnapWave process flags (IG waves, Herbers, wind, vegetation).
   !     Called from couple_snapwave (sfincs_snapwave) and
   !     read_snapwave_boundary_data (sfincs_bathtub).
   !
   use sfincs_log, only: write_log, logstr
   !
   implicit none
   !
   private
   public :: read_snapwave_input
   !
contains
   !
   subroutine read_snapwave_input()
   !
   ! Reads snapwave data from sfincs.inp
   !
   use snapwave_data   
   use sfincs_read   
   !
   implicit none
   !
   open(500, file='sfincs.inp')   
   !
   ! Input section
   !
   call get_keyword(500, 'snapwave_gamma',                  gamma,                0.7)
   call get_keyword(500, 'snapwave_gammax',                 gammax,               999.0)
   call get_keyword(500, 'snapwave_alpha',                  alpha,                1.0)
   call get_keyword(500, 'snapwave_hmin',                   hmin,                 0.1)
   call get_keyword(500, 'snapwave_fw',                     fw0,                  0.01)
   call get_keyword(500, 'snapwave_fwig',                   fw0_ig,               0.015)
   call get_keyword(500, 'snapwave_dt',                     dt,                   36000.0)
   call get_keyword(500, 'snapwave_tol',                    tol,                  1000.0)
   call get_keyword(500, 'snapwave_dtheta',                 dtheta,               10.0)
   call get_keyword(500, 'snapwave_crit',                   crit,                 0.001)
   call get_keyword(500, 'snapwave_nrsweeps',               nr_sweeps,            4)
   call get_keyword(500, 'snapwave_niter',                  niter,                10)                          !TL: Old default was 40
   !call get_keyword(500, 'snapwave_baldock_opt',           baldock_opt,          1)
   call get_keyword(500, 'snapwave_baldock_ratio',          baldock_ratio,        0.2)
   call get_keyword(500, 'snapwave_baldock_exponent',       baldock_exponent,     2)                           ! Exponent for multiplying the Baldock dissipation with a factor 'f = (Hloc / Hmax)**iexp' to enhance breaking when H > Hmax, with iexp = 0 (means unused), 1 or 2 (default). Generally, only active for steep coastlines, where Baldock dissipation can be too low in the surf zone.
   call get_keyword(500, 'rgh_lev_land',                    rghlevland,           0.0)
   call get_keyword(500, 'snapwave_fw_ratio',               fwratio,              1.0)
   call get_keyword(500, 'snapwave_fwig_ratio',             fwigratio,            1.0)
   call get_keyword(500, 'snapwave_Tpini',                  Tpini,                1.0)
   call get_keyword(500, 'snapwave_mwind',                  mwind,                2)
   call get_keyword(500, 'snapwave_sigmin',                 sigmin,               8.0 * atan(1.0) / 25.0)
   call get_keyword(500, 'snapwave_sigmax',                 sigmax,               8.0 * atan(1.0) / 1.0)
   call get_keyword(500, 'snapwave_jadcgdx',                jadcgdx,              1)
   call get_keyword(500, 'snapwave_c_dispT',                c_dispT,              1.0)
   call get_keyword(500, 'snapwave_sector',                 sector,               180.0)
   call get_keyword(500, 'snapwave_relax_factor_DoverA',    relax_factor_DoverA,  0.25)                        ! underrelaxation factor for DoverA (set to 1.0 to disable)
   call get_keyword(500, 'snapwave_relax_factor_DoverE',    relax_factor_DoverE,  0.25)                        ! underrelaxation factor for DoverE (set to 1.0 to disable)
   !
   ! Settings related to IG waves:   
   call get_keyword(500, 'snapwave_igwaves',                igwaves_opt,          1)
   call get_keyword(500, 'snapwave_alpha_ig',               alpha_ig,             1.0)                         !TODO choose whether snapwave_alphaig or snapwave_gamma_ig
   call get_keyword(500, 'snapwave_gammaig',                gamma_ig,             0.7)                         ! Wave breaking parameter for IG waves, default=0.7
   call get_keyword(500, 'snapwave_gamma_fac_br',           gamma_fac_br,         0.45)                        ! factor times gamma that is used to determine the maximum incident wave breaking point in the surf zone using local incident wave height over water depth ratio, among others used to set the IG source term to 0 shallower than this point
   call get_keyword(500, 'snapwave_shinc2ig',               shinc2ig,             1.0)                         ! Ratio of how much of the calculated IG wave source term, is subtracted from the incident wave energy (0-1, 1=default=all energy as sink)
   call get_keyword(500, 'snapwave_alphaigfac',             alphaigfac,           1.0)                         ! Multiplication factor for IG shoaling source/sink term
   call get_keyword(500, 'snapwave_baldock_ratio_ig',       baldock_ratio_ig,     0.2)
   call get_keyword(500, 'snapwave_ig_opt',                 ig_opt,               1)
   call get_keyword(500, 'snapwave_iterative_srcig',        iterative_srcig_opt,  0)                           ! Option whether to calculate IG source/sink term in iterative lower (better, but potentially slower, 1=default), or effectively based on previous timestep (faster, potential mismatch, =0)
   !
   ! IG boundary conditions options:
   call get_keyword(500, 'snapwave_use_herbers',            herbers_opt,          1)                           ! Choice whether you want IG Hm0&Tp be calculated by herbers (=1, default), or want to specify user defined values (0> then snapwave_eeinc2ig & snapwave_Tinc2ig are used)
   call get_keyword(500, 'snapwave_tpig_opt',               tpig_opt,             1)                           ! IG wave period option based on Herbers calculated spectrum, only used if snapwave_use_herbers = 1. Options are: 1=Tm01 (default), 2=Tpsmooth, 3=Tp, 4=Tm-1,0
   call get_keyword(500, 'snapwave_jonswapgamma',           jonswapgam,           3.3)                         ! JONSWAP gamma value for determination offshore spectrum and IG wave conditions using Herbers, default=3.3, only used if snapwave_use_herbers = 1
   call get_keyword(500, 'snapwave_eeinc2ig',               eeinc2ig,             0.01)                        ! Only used if snapwave_use_herbers = 0
   call get_keyword(500, 'snapwave_Tinc2ig',                Tinc2ig,              7.0)                         ! Only used if snapwave_use_herbers = 0
   !
   ! Wind
   !
   call get_keyword(500, 'snapwave_wind',                   wind_opt,             0)                           ! Flag whether to include windgrowth in SnapWave (1) or not (0, default)
   !
   ! Vegetation input
   !
   call get_keyword(500, 'snapwave_vegetation',             vegetation_opt,       0)
   !
   ! Input files
   !
   call get_keyword(500, 'snapwave_jonswapfile',            snapwave_jonswapfile, 'none')
   call get_keyword(500, 'snapwave_bndfile',                snapwave_bndfile,     'none')
   call get_keyword(500, 'snapwave_encfile',                snapwave_encfile,     'none')
   call get_keyword(500, 'snapwave_bhsfile',                snapwave_bhsfile,     'none')
   call get_keyword(500, 'snapwave_btpfile',                snapwave_btpfile,     'none')
   call get_keyword(500, 'snapwave_bwdfile',                snapwave_bwdfile,     'none')
   call get_keyword(500, 'snapwave_bdsfile',                snapwave_bdsfile,     'none')
   call get_keyword(500, 'snapwave_upwfile',                upwfile,              'snapwave.upw')
   call get_keyword(500, 'snapwave_mskfile',                mskfile,              'none')
   call get_keyword(500, 'snapwave_depfile',                depfile,              'none')
   call get_keyword(500, 'snapwave_ncfile',                 gridfile,             'snapwave_net.nc')
   call get_keyword(500, 'netsnapwavefile',                 netsnapwavefile,      'none')
   call get_keyword(500, 'storesnapwavegrid',               storesnapwavegrid,    .false.)
   call get_keyword(500, 'tref',                            trefstr,              '20000101 000000')           ! Read again > needed in sfincs_ncinput.F90
   !
   close(500)
   !
   igwaves          = .true.
   igherbers        = .false.
   iterative_srcig  = .false.   
   !
   if (igwaves_opt==0) then
      !
      igwaves       = .false.
      !
   else
      ! 
      if (iterative_srcig_opt==1) then
         iterative_srcig = .true.
      endif      
      !
      if (herbers_opt==0) then
         !
         write(logstr,*)'SnapWave: IG bc using use eeinc2ig= ',eeinc2ig,' and snapwave_Tinc2ig= ',Tinc2ig
         call write_log(logstr, 0)         
         !
      else
         !
         igherbers     = .true.          
         !
      endif
      !
   endif
   !
   wind = .true.
   if (wind_opt==0) then
      wind = .false.
   endif   
   !
   vegetation = .true.
   !
   if (vegetation_opt == 0) then
      vegetation = .false.
   endif   
   !
   if (nr_sweeps /= 1 .and. nr_sweeps /= 4) then
      !
      nr_sweeps = 4
      call write_log('SnapWave: Warning! nr_sweeps must be 1 or 4! Now set to 4.', 1)
      !
   endif
   ! 
   restart           = .true.
   coupled_to_sfincs = .true.
   !
   end subroutine read_snapwave_input

end module snapwave_input
