module snapwave_solver
   
      use sfincs_log
    
      implicit none
   contains
   
   subroutine compute_wave_field()
      !
      use snapwave_data
      !
      implicit none
      !
      !real*8, intent(in)   :: time > TL: not used in this implementation
      !
      real*4               :: tpb
      !
      real*4, parameter    :: waveps = 1e-5
      !real*4, dimension(:), allocatable   :: sig
      real*4, dimension(:), allocatable   :: sigm_ig
      real*4, dimension(:), allocatable   :: expon
      !
      integer   :: itheta
      integer   :: k
      !
      !allocate(sig(no_nodes))
      allocate(sigm_ig(no_nodes))
      !
      g     = 9.81
      pi    = 4.*atan(1.)
      !
      call timer(t0)
      !
      if (.not.restart) then
         !
         ! Set energies to 0.0; note that boundary values have been set in update_boundaries
         !
         do k = 1, no_nodes
            if (inner(k)) then
               ee(:,k) = waveps               
            endif
         enddo
         !
         ee_ig = waveps
         !
         restart=1 !TODO TL: CHECK > we need this turned on right now for IG...
         !
      endif
      !
      ! Initialize wave period
      !
      do k = 1, no_nodes
         if (inner(k)) then
            Tp(k) = Tpini
         endif
         if (neumannconnected(k)>0) then
            Tp(neumannconnected(k))=Tpini
         endif
      enddo
      !
      ! Compute celerities and refraction speed
      !
      Tp     = max(tpmean_bwv,Tpini)   ! to check voor windgroei
      sig    = 2.0*pi/Tp
      Tp_ig   = tpmean_bwv_ig! TL: now determined in snapwave_boundaries.f90 instead of Tinc2ig*Tp
      sigm_ig = 2.0*pi/Tp_ig !TODO - TL: Question do we want Tp_ig now as contant, or also spatially varying like Tp ?
      !
      expon   = -(sig*sqrt(depth/g))**(2.5)
      kwav    = sig**2/g*(1.0-exp(expon))**(-0.4)
      C       = sig/kwav
      nwav    = 0.5+kwav*depth/sinh(min(2*kwav*depth,50.0))
      Cg      = nwav*C
      !
      if (igwaves) then
         cg_ig   = Cg
         expon   = -(sigm_ig*sqrt(depth/g))**(2.5)
         kwav_ig  = sig**2/g*(1.0-exp(expon))**(-0.4)
      else
         cg_ig   = 0.0
         kwav_ig = 0.0 
      endif
      !
      do k = 1, no_nodes
         sinhkh(k)    = sinh(min(kwav(k)*depth(k), 50.0))
         Hmx(k)       = gamma*depth(k)
      enddo
      if (igwaves) then          
         do k = 1, no_nodes             
            Hmx_ig(k)    = 0.88/kwav_ig(k)*tanh(gamma_ig*kwav_ig(k)*depth(k)/0.88) ! Note - uses gamma_ig
         enddo
      else
         Hmx_ig    = 0.0
      endif
      !   
      do itheta = 1, ntheta
         ctheta(itheta,:)    = sig/sinh(min(2.0*kwav*depth, 50.0))*(dhdx*sin(theta(itheta)) - dhdy*cos(theta(itheta)))
      enddo
      !
      if (igwaves) then
         do itheta = 1, ntheta
            ctheta_ig(itheta,:) = sigm_ig/sinh(min(2.0*kwav_ig*depth, 50.0))*(dhdx*sin(theta(itheta)) - dhdy*cos(theta(itheta)))
         enddo
      else
         ctheta_ig = 0.0
      endif
      !
      ! Limit unrealistic refraction speed to 1/2 pi per wave period
      !
      do k = 1, no_nodes
         ctheta(:,k)    = sign(1.0, ctheta(:,k))*min(abs(ctheta(:, k)), sig(k)/4)
      enddo
      !
      if (igwaves) then
           do k=1, no_nodes
              ctheta_ig(:,k) = sign(1.0, ctheta_ig(:,k))*min(abs(ctheta_ig(:, k)), sigm_ig(k)/4.0)
           enddo
      endif
      !
      ! Solve the directional wave energy balance on an unstructured grid
      !
      call timer(t2)
      !
      call solve_energy_balance2Dstat (x,y,dhdx, dhdy, no_nodes,inner, &
                                         w, ds, prev,   &
                                         neumannconnected,       &
                                         theta,ntheta,thetamean, &
                                         depth,kwav,cg,ctheta,fw,     &
                                         Tp,Tp_ig,dt,rho,alpha,gamma, gammax, &
                                         wind,   &
                                         H,Dw,F,Df,thetam,sinhkh,&
                                         Hmx, ee, windspreadfac, u10, niter, crit, &
                                         hmin, baldock_ratio, baldock_ratio_ig, &
                                         aa, sig, jadcgdx, sigmin, sigmax,&
                                         c_dispT, WsorE, WsorA, SwE, SwA, Tpini, &
                                         igwaves,kwav_ig, cg_ig,H_ig,ctheta_ig,Hmx_ig, ee_ig,fw_ig, &
                                         beta, srcig, alphaig, Dw_ig, Df_ig, &
                                         vegetation, no_secveg, veg_ah, veg_bstems, veg_Nstems, veg_Cd, Dveg, &
                                         zb, nwav, ig_opt, alpha_ig, gamma_ig, eeinc2ig, Tinc2ig, alphaigfac, shinc2ig, iterative_srcig)
      !
      call timer(t3)
      !
      Fx = F*cos(thetam)
      Fy = F*sin(thetam)
      !   
   end subroutine
   
   
   subroutine solve_energy_balance2Dstat(x,y,dhdx, dhdy, no_nodes,inner, &
                                         w, ds, prev,   &
                                         neumannconnected,       &
                                         theta,ntheta,thetamean, &
                                         depth,kwav,cg,ctheta,fw,     &
                                         Tp,T_ig,dt,rho,alfa,gamma, gammax, &
                                         wind,   &
                                         H,Dw,F,Df,thetam,sinhkh,&
                                         Hmx, ee, windspreadfac, u10, niter, crit, &
                                         hmin, baldock_ratio, baldock_ratio_ig, &       
                                         aa, sig, jadcgdx, sigmin, sigmax,&
                                         c_dispT, WsorE, WsorA, SwE, SwA, Tpini, &
                                         igwaves,kwav_ig, cg_ig,H_ig,ctheta_ig,Hmx_ig, ee_ig,fw_ig, &
                                         betamean, srcig, alphaig, Dw_ig, Df_ig, &       
                                         vegetation, no_secveg, veg_ah, veg_bstems, veg_Nstems, veg_Cd, Dveg, &
                                         zb, nwav, ig_opt, alfa_ig, gamma_ig, eeinc2ig, Tinc2ig, alphaigfac, shinc2ig, iterative_srcig)
   !
   use snapwave_windsource
   !use snapwave_ncoutput ! TL: removed, we don't use this in SF+SW
   !
   implicit none
   !
   ! In/output variables and arrays
   !
   real*8, dimension(no_nodes),intent(in)           :: x,y                    ! x,y coordinates of grid
   real*4, dimension(no_nodes),intent(in)           :: dhdx, dhdy             ! bed level gradients in x and y direction
   integer, intent(in)                              :: no_nodes,ntheta        ! number of grid points, number of directions
   real*4,  dimension(2,ntheta,no_nodes),intent(in) :: w                      ! weights of upwind grid points, 2 per grid point and per wave direction
   integer, dimension(2,ntheta,no_nodes),intent(in) :: prev                   ! two upwind grid points per grid point and wave direction
   real*4, dimension(ntheta,no_nodes), intent(in)   :: ds                     ! distance to interpolated upwind point, per grid point and direction
   real*4, intent(in)                               :: thetamean              ! mean offshore wave direction (rad)
   real*4, dimension(ntheta), intent(in)            :: theta                  ! distribution of wave angles and offshore wave energy density
   logical, dimension(no_nodes), intent(inout)      :: inner                  ! mask of inner grid points (not on boundary)
   logical, intent(in)                              :: wind                   ! logical wind on/off
   logical, intent(in)                              :: igwaves                ! logical IG waves on/off
   logical, intent(in)                              :: iterative_srcig        ! option whether IG source/sink term should be calculated implicitly (iterative_srcig = 1), or just a priori based on effectively incident wave energy from previous timestep only (iterative_srcig=0, explicitly, bit faster)
   integer, dimension(no_nodes),intent(in)          :: neumannconnected       ! number of neumann boundary point if connected to inner point
   real*4, dimension(no_nodes), intent(in)          :: depth                  ! water depth
   real*4, dimension(no_nodes), intent(in)          :: zb                     ! actual bed level
   real*4, dimension(no_nodes), intent(inout)       :: kwav                   ! wave number
   real*4, dimension(no_nodes), intent(in)          :: kwav_ig                ! wave number
   real*4, dimension(no_nodes), intent(inout)       :: cg                     ! group velocity
   real*4, dimension(ntheta,no_nodes), intent(inout):: ctheta                 ! refractioon speed
   real*4, dimension(no_nodes), intent(in)          :: cg_ig                  ! group velocity
   real*4, dimension(no_nodes), intent(in)          :: nwav                   ! wave number n      
   real*4, dimension(ntheta,no_nodes), intent(inout):: ee                     ! 
   real*4, dimension(ntheta,no_nodes), intent(inout):: ee_ig                  ! 
   real*4, dimension(ntheta,no_nodes), intent(in)   :: ctheta_ig              ! refractioon speed
   real*4, dimension(no_nodes), intent(in)          :: fw                     ! wave friction factor
   real*4, dimension(no_nodes), intent(in)          :: fw_ig                  ! wave friction factor
   real*4, dimension(no_nodes), intent(out)         :: betamean               ! Mean local bed slope parameter  
   real*4, dimension(no_nodes), intent(out)         :: srcig                  ! Directionally averaged incident wave sink/infragravity source term 
   real*4, dimension(no_nodes), intent(out)         :: alphaig                ! Mean IG shoaling parameter alpha    
   real*4, intent(in)                               :: dt                     ! time step (s)
   real*4, intent(in)                               :: rho                    ! water density
   real*4, intent(in)                               :: alfa,gamma, gammax     ! coefficients in Baldock wave breaking dissipation
   real*4, intent(in)                               :: baldock_ratio          ! option controlling from what depth wave breaking should take place: (Hk>baldock_ratio*Hmx(k)), default baldock_ratio=0.2
   real*4, intent(in)                               :: baldock_ratio_ig       ! option controlling from what depth wave breaking should take place for IG waves: (Hk_ig>baldock_ratio_ig*Hmx_ig(k)), default baldock_ratio_ig=0.2     
   real*4, dimension(no_nodes), intent(inout)       :: H                      ! wave height - TODO - TL - CHECK > inout needed to have updated 'H' for determining srcig
   real*4, dimension(no_nodes), intent(out)         :: H_ig                   ! wave height
   real*4, dimension(no_nodes), intent(out)         :: Dw                     ! wave breaking dissipation
   real*4, dimension(no_nodes), intent(out)         :: Dw_ig                  ! wave breaking dissipation IG   
   real*4, dimension(no_nodes), intent(out)         :: F                      ! wave force Dw/C/rho/h   
   real*4, dimension(no_nodes), intent(out)         :: Df                     ! wave friction dissipation
   real*4, dimension(no_nodes), intent(out)         :: Df_ig                  ! wave friction dissipation IG      
   real*4, dimension(no_nodes), intent(out)         :: thetam                 ! mean wave direction
   real*4, dimension(no_nodes), intent(inout)       :: sinhkh                 ! sinh(k*depth)
   real*4, dimension(no_nodes), intent(inout)       :: Hmx                    ! Hmax
   real*4, dimension(no_nodes), intent(inout)       :: Tp                     ! Peak wave period
   real*4, dimension(no_nodes), intent(in)          :: T_ig                   ! IG wave period      
   real*4, dimension(no_nodes), intent(in)          :: Hmx_ig                 ! Hmax
   real*4, dimension(no_nodes), intent(in)          :: u10                    ! wind speed and direction
   integer,                     intent(in)          :: niter                  ! max number of iterations
   real*4,                      intent(in)          :: crit                   ! relative accuracy for stopping criterion
   integer                                          :: ig_opt                 ! option of IG wave settings (1 = default = conservative shoaling based dSxx and Baldock breaking)
   !
   ! wind source vars
   !
   integer, intent(in)                                 :: jadcgdx               ! logical yes/no               ! logical yes/no
   real*4, intent(in)                                  :: sigmin, sigmax, c_dispT
   real*4, dimension(ntheta, no_nodes), intent(in)     :: windspreadfac        !< [-] distribution array for wind input
   real*4, dimension(ntheta,no_nodes),  intent(inout)  :: aa    
   real*4, dimension(ntheta,no_nodes),  intent(out)    :: WsorE, WsorA
   real*4, dimension(no_nodes),         intent(out)    :: SwE, SwA
   real*4, dimension(no_nodes),         intent(inout)  :: sig
   real*4, intent(in)                                  :: Tpini
   !
   ! vegetation vars
   !
   logical, intent(in)                                  :: vegetation               ! logical yes/no
   real*4, dimension(no_nodes), intent(out)             :: Dveg                     ! dissipation by vegetation: N.B. spatial field!
   integer, intent(in)                                  :: no_secveg                ! number of sections in the vertical 
   real*4, dimension(no_nodes,no_secveg), intent(in)    :: veg_ah                   ! Height of vertical sections used in vegetation schematization [m wrt zb_ini (zb0)]
   real*4, dimension(no_nodes,no_secveg), intent(in)    :: veg_bstems               ! Width/diameter of individual vegetation stems [m]
   real*4, dimension(no_nodes,no_secveg), intent(in)    :: veg_Nstems               ! Number of vegetation stems per unit horizontal area [m-2]
   real*4, dimension(no_nodes,no_secveg), intent(in)    :: veg_Cd                   ! Bulk drag coefficient [-]     
   real*4                                               :: Dvegk                    ! dissipation by vegetation: N.B. scalar value!
   real*4, dimension(no_nodes)                          :: Fvw                      ! vegetation wave drag force   
   real*4, dimension(no_nodes,50)                       :: unl                      ! non-linear wave orbital velocity time series, in 50 points per wave length
   !
   !
   ! Local variables and arrays
   !
   integer, dimension(:), allocatable         :: ok                     ! mask for fully iterated points
   real*4                                     :: eemax,dtheta           ! maximum wave energy density, directional resolution
   real*4                                     :: uorbi
   integer                                    :: sweep,iter            ! sweep number, number of iterations
   integer                                    :: k,k1,k2,count,kn,itheta ! counters (k is grid index)
   integer, dimension(:,:), allocatable       :: indx                   ! index for grid sorted per sweep direction
   real*4, dimension(:,:), allocatable        :: eeold               ! wave energy density, energy density previous iteration
   real*4, dimension(:), allocatable          :: Eold                   ! mean wave energy, previous iteration
   real*4, dimension(:,:), allocatable        :: srcig_local            ! Energy source/sink term because of IG wave energy transfer from incident waves
   real*4, dimension(:,:), allocatable        :: beta_local             ! Local bed slope based on bed level per direction      
   real*4, dimension(:,:), allocatable        :: alphaig_local          ! Local infragravity wave shoaling parameter alpha
   real*4, dimension(:,:), allocatable        :: depthprev              ! water depth at upwind intersection point per direction     
   real*4, dimension(:), allocatable          :: dee                    ! difference with energy previous iteration
   real*4, dimension(:), allocatable          :: eeprev, cgprev         ! energy density and group velocity at upwind intersection point
   real*4, dimension(:), allocatable          :: eeprev_ig, cgprev_ig   ! energy density and group velocity at upwind intersection point
   real*4, dimension(:), allocatable          :: A,B,C,R                ! coefficients in the tridiagonal matrix solved per point
   real*4, dimension(:), allocatable          :: B_aa,R_aa,aaprev       ! coefficients in the tridiagonal matrix solved per point
   real*4, dimension(:), allocatable          :: A_ig,B_ig,C_ig,R_ig    ! coefficients in the tridiagonal matrix solved per point
   real*4, dimension(:), allocatable          :: DoverE                 ! ratio of mean wave dissipation over mean wave energy
   real*4, dimension(:), allocatable          :: DoverA                 ! ratio of mean wave dissipation over mean wave energy
   real*4, dimension(:), allocatable          :: DoverE_ig              ! ratio of mean wave dissipation over mean wave energy
   real*4, dimension(:), allocatable          :: E                      ! mean wave energy
   real*4, dimension(:), allocatable          :: E_ig                   ! mean wave energy
   real*4, dimension(:), allocatable          :: diff                   ! maximum difference of wave energy relative to previous iteration
   real*4, dimension(:), allocatable          :: ra                     ! coordinate in sweep direction
   !real*4, dimension(:), allocatable          :: sig
   real*4, dimension(:), allocatable          :: sigm_ig
   integer, dimension(4)                      :: shift
   real*4                                     :: pi = 4.*atan(1.0)
   real*4                                     :: g=9.81
   real*4                                     :: hmin                   ! minimum water depth! TL: make user changeable also here according to 'snapwave_hmin' in sfincs.inp   
   real*4                                     :: fac=1.0             ! underrelaxation factor for DoverA
   real*4                                     :: oneoverdt
   real*4                                     :: oneover2dtheta
   real*4                                     :: rhog8
   real*4                                     :: Dfk
   real*4                                     :: Dwk
   real*4                                     :: Ek
   real*4                                     :: Hk
   real*4                                     :: percok
   real*4                                     :: error
   real*4                                     :: Dfk_ig
   real*4                                     :: Dwk_ig
   real*4                                     :: Ek_ig
   real*4                                     :: Hk_ig
   real*4                                     :: alfa_ig,gamma_ig  ! coefficients in Baldock wave breaking dissipation model for IG waves
   real*4                                     :: eeinc2ig          ! ratio of incident wave energy as first estimate of IG wave energy at boundary
   real*4                                     :: Tinc2ig           ! ratio compared to period Tinc to estimate Tig
   real*4                                     :: alphaigfac        ! Multiplication factor for IG shoaling source/sink term, default = 1.0
   real*4                                     :: shinc2ig          ! Ratio of how much of the calculated IG wave source term, is subtracted from the incident wave energy (0-1, 0=default)   
   integer, save                              :: callno=1
   !
   real*4, dimension(ntheta)                  :: sinth,costh            ! distribution of wave angles and offshore wave energy density   
   !
   !local wind source vars
   !
   real*4                                     :: Ak
   real*4                                     :: DwT
   real*4                                     :: DwAk
   real*4                                     :: ndissip      !       
   real*4                                     :: depthlimfac=1.0
   real*4                                     :: waveps=0.0001
   !
   ! Allocate local arrays
   !
   waveps         = 0.0001
   allocate(ok(no_nodes)); ok=0
   allocate(indx(no_nodes,4)); indx=0
   allocate(eeold(ntheta,no_nodes)); eeold=0.0
   allocate(dee(ntheta)); dee=0.0
   allocate(eeprev(ntheta)); eeprev=0.0
   allocate(cgprev(ntheta)); cgprev=0.0
   allocate(A(ntheta)); A=0.0
   allocate(B(ntheta)); B=0.0
   allocate(C(ntheta)); C=0.0
   allocate(R(ntheta)); R=0.0
   allocate(DoverE(no_nodes)); DoverE=0.0
   allocate(E(no_nodes)); E=waveps
   allocate(Eold(no_nodes)); Eold=0.0   
   !
   if (igwaves) then
      allocate(A_ig(ntheta)); A_ig=0.0
      allocate(B_ig(ntheta)); B_ig=0.0
      allocate(C_ig(ntheta)); C_ig=0.0
      allocate(R_ig(ntheta)); R_ig=0.0
      allocate(eeprev_ig(ntheta)); eeprev_ig=0.0
      allocate(cgprev_ig(ntheta)); cgprev_ig=0.0
      allocate(DoverE_ig(no_nodes)); DoverE_ig=0.0
      allocate(E_ig(no_nodes)); E_ig=waveps
      !allocate(T_ig(no_nodes)); T_ig=0.0
      allocate(sigm_ig(no_nodes)); sigm_ig=0.0
      allocate(depthprev(ntheta,no_nodes)); depthprev=0.0                     
      allocate(beta_local(ntheta,no_nodes)); beta_local=0.0        
      allocate(alphaig_local(ntheta,no_nodes)); alphaig_local=0.0
   endif
   !
   if (wind) then
      allocate(B_aa(ntheta)); B_aa=0.0
      allocate(R_aa(ntheta)); R_aa=0.0
      allocate(DoverA(no_nodes)); DoverA=0.0
      allocate(aaprev(ntheta)); aaprev=0.0
   endif
   !
   allocate(diff(no_nodes)); diff=0.0
   allocate(ra(no_nodes)); ra=0.0
   allocate(srcig_local(ntheta,no_nodes)); srcig_local = 0.0
   !
   do itheta = 1, ntheta
      sinth(itheta) = sin(theta(itheta))
      costh(itheta) = cos(theta(itheta))
   enddo
   !   
   df   = 0.0
   dw   = 0.0
   F    = 0.0
   !
   ok             = 0
   indx           = 0
   eemax          = maxval(ee)
   dtheta         = theta(2) - theta(1)
   if (dtheta<0.)   dtheta = dtheta + 2.*pi
   if (wind) then
      sig         = 2*pi/Tpini
   else
      sig         = 2*pi/Tp
   endif
   oneoverdt      = 1.0/dt
   oneover2dtheta = 1.0/2.0/dtheta
   rhog8          = 0.125*rho*g
   thetam         = 0.0
   !H              = 0.0 ! TODO - TL: CHeck > needed for restart for IG > set to 0 now in snapwave_domain.f90
   Dveg           = 0.0
   Fvw            = 0.0
   unl            = 0.0
   !
   if (igwaves) then
      !T_ig = Tinc2ig*Tp
      sigm_ig        = 2*pi/T_ig
      DoverE_ig      = 0.0
   endif
   !
   if (wind) then
      DoverA = 0.0      
      ndissip = 3.0
      WsorE = 0.0
      WsorA = 0.0
      Ak = waveps/sigmax
   endif
   !   
   ! Sort coordinates in sweep directions
   !
   shift = [0,1,-1,2]
   do sweep = 1, 4
      !
      ra = x*cos(thetamean + 0.5*pi*shift(sweep)) + y*sin(thetamean + 0.5*pi*shift(sweep))
      call hpsort_eps_epw(no_nodes, ra , indx(:, sweep), 1.0e-6)
      !
   enddo
   !      
   ! Set inner to false for all points at grid edge or adjacent to dry point
   !
   do k=1,no_nodes
      ! 
      do itheta = 1, ntheta
         !
         k1 = prev(1, itheta, k)
         k2 = prev(2, itheta, k)
         !
         if (k1 * k2 == 0) then
            !
            ! No upwind point (is this check not already done somewhere before?)
            !
            inner(k) = .false.
            !
         elseif (depth(k1) < hmin .or. depth(k2) < hmin .or. (k1 == 1 .and. k2 == 1)) then
            !
            ! Do not change inner here! It should be static! In a next update of the wave fields, these points may be wet.
            !
            !inner(k) = .false.
            !
            !exit
            !
         endif
      enddo
   enddo
   !
   !
   ! 0-a) Set boundary and initial conditions      
   !
   do k = 1, no_nodes
      !
      ! Boundary condition at sea side (uniform)
      !
      if (.not.inner(k)) then
         !
         ee(:,k)=max(ee(:,k),waveps)
         E(k)      = sum(ee(:, k))*dtheta
         H(k)      = sqrt(8*E(k)/rho/g)
         thetam(k) = atan2(sum(ee(:, k)*sin(theta)), sum(ee(:, k)*cos(theta)))
         !
         !ee_ig(:, k)  = eeinc2ig*ee(:,k) !TODO TL: determined in snapwave_boundaries.f90        
         !
         if (igwaves) then             
            E_ig(k)      = sum(ee_ig(:, k))*dtheta
            H_ig(k)      = sqrt(8*E_ig(k)/rho/g)
         endif
         ! 
         if (wind) then
            sig(k)    = 2*pi/Tp(k)
           !aa(:,k)   = max(aa(:,k),waveps/sig(k))
            aa(:,k)   = max(ee(:,k),waveps)/sig(k)
            Ak        = E(k)/sig(k)
         endif
         !
      endif
   enddo
   !
   ! 0-b) Determine IG source/sink term
   !
   if (igwaves) then
      !        
      ! As defined in Leijnse, van Ormondt, van Dongeren, Aerts & Muis et al. 2024 
      !      
      ! Actual determining of source term: 
      !          
      call determine_infragravity_source_sink_term(inner, no_nodes, ntheta, w, ds, prev, cg_ig, nwav, depth, zb, H, ee, ee_ig, eeprev, eeprev_ig, cgprev, ig_opt, alphaigfac, alphaig_local, beta_local, srcig_local) 
      !         
      ! inout: alphaig_local, srcig_local - eeprev, eeprev_ig, cgprev, beta_local
      ! in: the rest
      !
      ! NOTE - This is based on the energy in the precious SnapWave timestep 'ee' and 'ee_ig', and waveheight 'H', which should therefore be made available. 
      !
   endif             
   !
   ! 0-c) Set initial condition at inner cells  
   !
   do k = 1, no_nodes   
      ! 
      if (inner(k)) then
         !
         if (wind) then
             ee(:,k) = waveps
             sig(k)  = 2*pi/Tpini
             aa(:,k) = ee(:,k)/sig(k)
         else
             ee(:,k) = waveps
         endif
         !
         ! Make sure DoverE is filled based on previous ee
         Ek       = sum(ee(:, k))*dtheta      
         Hk       = min(sqrt(Ek/rhog8), gamma*depth(k))
         Ek       = rhog8*Hk**2
         if (.not. wind) then
            uorbi    = 0.5*sig(k)*Hk/sinhkh(k)
            Dfk      = 0.28*rho*fw(k)*uorbi**3
            call baldock(rho, g, alfa, gamma, depth(k), Hk, Tp(k), 1, Dwk, Hmx(k))
            DoverE(k) = (Dwk + Dfk)/max(Ek, 1.0e-6)
         endif
         !
         if (wind) then
            call compute_celerities(depth(k), sig(k), sinth, costh, ntheta, gamma, dhdx(k), dhdy(k), sinhkh(k), Hmx(k), kwav(k), cg(k), ctheta(:,k))  
            uorbi    = 0.5*sig(k)*Hk/sinhkh(k)  
            Dfk      = 0.28*rho*fw(k)*uorbi**3
            call baldock(rho, g, alfa, gamma, depth(k), Hk, 2.0*pi/sig(k), 1, Dwk, Hmx(k))
            DoverE(k) = (Dwk + Dfk)/max(Ek, 1.0e-6)      
            !
            ! initial conditions are not equal to bc conditions  
            DwT = - c_dispT/(1.0 -ndissip)*(2.*pi)/sig(k)**2*cg(k)*kwav(k) * DoverE(k)
            DwAk = 0.5/pi * (E(k)*DwT+2.0*pi*Ak*DoverE(k) )
            DoverA(k) = DwAk/max(Ak,1e-6) 
         endif
         !
      endif
      !
   enddo        
   !
   ! Start iteration
   !
   do iter=1,niter
      !
      sweep = mod(iter, 4) !TODO - TL: problem that we don't have option for sweep = 1 anymore?
      !
      if (sweep==0) then
         sweep = 4
      endif
      !
      if (sweep==1) then
         eeold = ee
         do k = 1, no_nodes
            Eold(k) = sum(eeold(:, k))
         enddo
         !
         if (igwaves) then
            !        
            if (iterative_srcig) then 
                ! Update H(k) based on updated ee(:,k), as used in IG source term to determine alphaig
                ! 
                do k = 1, no_nodes   
                   ! 
                   if (inner(k)) then             
                       !
                       H(k)      = sqrt(8*sum(ee(:, k))*dtheta/rho/g)                   
                       !
                   endif               
                enddo
                !
                ! Actual determining of source term - every first sweep of iteration
                !          
                call determine_infragravity_source_sink_term(inner, no_nodes, ntheta, w, ds, prev, cg_ig, nwav, depth, zb, H, ee, ee_ig, eeprev, eeprev_ig, cgprev, ig_opt, alphaigfac, alphaig_local, beta_local, srcig_local) 
                !    
            endif
            !
         endif   
         !
      endif
      !
      !  Loop over all points depending on sweep direction
      !
      do count = 1, no_nodes
         !
         k=indx(count, sweep)
         !
         if (inner(k)) then
            if (depth(k)>1.1*hmin) then
               !
               if (ok(k) == 0)  then
                  !
                  ! Only perform computations on wet inner points that are not yet converged (ok)
                  !
                  do itheta = 1, ntheta
                     !
                     k1 = prev(1, itheta, k)
                     k2 = prev(2, itheta, k)
                     !
                     eeprev(itheta) = w(1, itheta, k)*ee(itheta, k1) + w(2, itheta, k)*ee(itheta, k2)
                     cgprev(itheta) = w(1, itheta, k)*cg(k1) + w(2, itheta, k)*cg(k2)
                     !
                     if (igwaves) then
                        eeprev_ig(itheta) = w(1, itheta, k)*ee_ig(itheta, k1) + w(2, itheta, k)*ee_ig(itheta, k2)
                        cgprev_ig(itheta) = w(1, itheta, k)*cg_ig(k1) + w(2, itheta, k)*cg_ig(k2)
                     endif
                     !
                     if (wind) then
                        aaprev(itheta) = w(1, itheta, k)*aa(itheta, k1) + w(2, itheta, k)*aa(itheta, k2)
                     endif
                     !
                  enddo
                  !
                  Ek = sum(eeprev)*dtheta     ! to check                
                  !
                  depthlimfac = max(1.0, (sqrt(Ek/rhog8)/(gammax*depth(k)))**2.0)
                  Hk = min(sqrt(Ek/rhog8), gamma*depth(k))
                  Ek = Ek/depthlimfac
                  !
                  if (wind) then
                    !
                    Ak = sum(aaprev)*dtheta
                    !
                    Ak = Ak/depthlimfac                   
                    ee(:,k) = ee(:,k) / depthlimfac
                    aa(:,k) = aa(:,k) / depthlimfac                   
                    sig(k)  = Ek/Ak
                    sig(k)  = max(sig(k),sigmin)
                    sig(k)  = min(sig(k),sigmax)
                    Ak      = Ek/sig(k) ! to avoid small T in windinput
                    if (wind) then
                       aaprev=min(aaprev,eeprev/sigmin)
                       aaprev=max(aaprev,eeprev/sigmax)
                    endif
                    !                    
                    call compute_celerities(depth(k), sig(k), sinth, costh, ntheta, gamma, dhdx(k), dhdy(k), sinhkh(k), Hmx(k), kwav(k), cg(k), ctheta(:,k))    
                  endif
                  ! 
                  ! Fill DoverE 
                  uorbi    = 0.5*sig(k)*Hk/sinhkh(k)
                  Dfk      = 0.28*rho*fw(k)*uorbi**3
                  !if (Hk>0.) then !
                  if (Hk>baldock_ratio*Hmx(k)) then
                     call baldock(rho, g, alfa, gamma, depth(k), Hk, 2*pi/sig(k) , 1, Dwk, Hmx(k))
                  else
                     Dwk   = 0.
                  endif
                  !                  
                  if (vegetation) then
                      call vegatt(sig(k), no_nodes, kwav(k), no_secveg, veg_ah(k,:), veg_bstems(k,:), veg_Nstems(k,:), veg_Cd(k,:), depth(k), rho, g, Hk, Dvegk)
                      !call vegatt(sig(k), no_nodes, kwav(k), no_secveg, (/6.5/), (/0.3/), (/0.7/), (/1.0/), depth(k), rho, g, Hk, Dvegk)                    
                  else
                      Dvegk = 0.
                  endif
                  !
                  DoverE(k) = (Dwk + Dfk + Dvegk)/max(Ek, 1.0e-6)
                  !if (Dvegk > 0.0) then
                  !   write(logstr,*)'k ',k,'depth(k)',depth(k),'Hk ',Hk,'Ek ',Ek,'Dwk ', Dwk,'Dfk ', Dfk,'Dvegk ', Dvegk,'DoverE veggie ', DoverE(k),'DoverE org ', (Dwk + Dfk)/max(Ek, 1.0e-6)
                  !   call write_log(logstr, 0)                                    
                  !endif
                  !
                  if (wind) then
                     !
                     if (iter==1) then
                        call windinput(u10(k), rho, g, depth(k), ntheta, windspreadfac(:,k), Ek, Ak, cg(k), eeprev, aaprev, ds(:,k), WsorE(:,k), WsorA(:,k), jadcgdx)
                     else
                        call windinput(u10(k), rho, g, depth(k), ntheta, windspreadfac(:,k), Ek, Ak, cg(k), ee(:,k), aa(:,k), ds(:,k), WsorE(:,k), WsorA(:,k), jadcgdx)   
                     endif
                     !
                     DwT = - c_dispT/(1.0 -ndissip)*(2.0*pi)/sig(k)**2*cg(k)*kwav(k) * DoverE(k) 
                     DwAk = 1/2.0/pi * (E(k)*DwT+2.0*pi*Ak*DoverE(k) )
                     !
                     if (iter==1) then
                        DoverA(k) = DwAk/max(Ak,1e-6) 
                     else
                        DoverA(k) = (1.0-fac)*DoverA(k)+fac*DwAk/max(Ak,1.e-6) 
                     endif
                     !                     
                     call compute_celerities(depth(k), sig(k), sinth, costh, ntheta, gamma, dhdx(k), dhdy(k), sinhkh(k), Hmx(k), kwav(k), cg(k), ctheta(:,k))  
                     !
                  endif
                  !
                   do itheta = 1, ntheta
                     !
                     R(itheta) = oneoverdt*ee(itheta, k) + cgprev(itheta)*eeprev(itheta)/ds(itheta, k) - srcig_local(itheta, k) * shinc2ig
                     !
                  enddo                  
                  !
                  do itheta = 2, ntheta - 1
                     !
                     A(itheta) = -ctheta(itheta - 1, k)*oneover2dtheta
                     B(itheta) = oneoverdt + cg(k)/ds(itheta,k) + DoverE(k)
                     C(itheta) = ctheta(itheta + 1, k)*oneover2dtheta
                     !
                  enddo
                  !
                  A(1) = -ctheta(ntheta, k)*oneover2dtheta
                  B(1) = oneoverdt + cg(k)/ds(1,k) + DoverE(k)
                  C(1) = ctheta(2, k)*oneover2dtheta
                  !
                  A(ntheta) = -ctheta(ntheta - 1, k)*oneover2dtheta
                  B(ntheta) = oneoverdt + cg(k)/ds(ntheta,k) + DoverE(k)
                  C(ntheta) = ctheta(1, k)*oneover2dtheta
                  !
                  ! Solve tridiagonal system per point
                  !
                  if (wind) then
                     do itheta = 2, ntheta - 1
                        B_aa(itheta) = oneoverdt + cg(k)/ds(itheta,k) + DoverA(k)       
                        R_aa(itheta) = (oneoverdt)*aa(itheta, k) + cgprev(itheta)*aaprev(itheta)/ds(itheta, k)
                     enddo
                     !
                     if (ctheta(1,k)<0) then
                        B_aa(1) = oneoverdt - ctheta(1, k)/dtheta + cg(k)/ds(1, k) + DoverA(k)
                        R_aa(1) = (oneoverdt)*aa(1, k) + cgprev(1)*aaprev(1)/ds(1, k) 
                     else
                        B_aa(1) = oneoverdt + cg(k)/ds(1, k) + DoverA(k)
                        R_aa(1) = (oneoverdt)*aa(1, k) + cgprev(1)*aaprev(1)/ds(1, k)
                     endif
                     !
                     if (ctheta(ntheta, k)>0) then
                        B_aa(ntheta) = oneoverdt + ctheta(ntheta, k)/dtheta + cg(k)/ds(ntheta, k) + DoverA(k)
                        R_aa(ntheta) = (oneoverdt )*aa(ntheta,k) + cgprev(ntheta)*aaprev(ntheta)/ds(ntheta, k)
                     else
                        B_aa(ntheta) = oneoverdt + cg(k)/ds(ntheta, k) + DoverA(k)
                        R_aa(ntheta) = (oneoverdt)*aa(ntheta,k) + cgprev(ntheta)*aaprev(ntheta)/ds(ntheta, k)
                     endif
                     R(:)    = R(:)    + WsorE(:,k)
                     R_aa(:) = R_aa(:) + WsorA(:,k)
                     !
                     call solve_tridiag(A, B, C, R, ee(:,k), ntheta)
                     call solve_tridiag(A,B_aa,C,R_aa,aa(:,k),ntheta)
                     ee(:, k) = max(ee(:, k), waveps)
                     aa(:,k) = max(aa(:,k),waveps/sigmax)
                     aa(:,k) = max(aa(:,k),waveps/sig(k))
                     !
                     Ek       = sum(ee(:, k))*dtheta     
                     Ak       = sum(aa(:,k))*dtheta
                     !
                     depthlimfac = max(1.0, (sqrt(Ek/rhog8)/(gammax*depth(k)))**2.0)
                     Hk = sqrt(Ek/rhog8/depthlimfac)
                     Ek = Ek/depthlimfac
                     Ak = Ak/depthlimfac
                     ee(:,k) = ee(:,k)/depthlimfac
                     aa(:,k) = aa(:,k)/depthlimfac
                     !
                     sig(k)  = Ek/Ak
                     sig(k)  = max(sig(k),sigmin)
                     sig(k)  = min(sig(k),sigmax)
                     call compute_celerities(depth(k), sig(k), sinth, costh, ntheta, gamma, dhdx(k), dhdy(k), sinhkh(k), Hmx(k), kwav(k), cg(k), ctheta(:,k))   
                     if (sig(k)<0.1) then
                         a=1
                     endif
                  else
                     !
                     ! Solve tridiagonal system per point
                     !
                     call solve_tridiag(A, B, C, R, ee(:,k), ntheta)
                     ee(:, k) = max(ee(:, k),waveps)
                     !
                  endif !wind 
                  !
                  ! IG
                  !
                  if (igwaves) then
                     Ek_ig       = sum(eeprev_ig)*dtheta                  
                     !Hk_ig       = sqrt(Ek_ig/rhog8) !org trunk
                     Hk_ig       = min(sqrt(Ek_ig/rhog8), gamma_ig*depth(k))  !TL: Question - why not this one?                     
                     Ek_ig       = rhog8*Hk_ig**2
                     ! 
                     ! Bottom friction Henderson and Bowen (2002) - D = 0.015*rhow*(9.81/depth(k))**1.5*(Hk/sqrt(8.0))*Hk_ig**2/8
                     !
                     Dfk_ig      = fw_ig(k)*0.0361*(9.81/depth(k))**1.5*Hk*Ek_ig
                     !
                     ! Dissipation of infragravity waves
                     !
                     if (Hk_ig>baldock_ratio_ig*Hmx_ig(k)) then
                        call baldock(rho, g, alfa_ig, gamma_ig, depth(k), Hk_ig, T_ig(k), 1, Dwk_ig, Hmx_ig(k))
                     else
                        Dwk_ig   = 0.
                     endif
                     !
                     DoverE_ig(k) = (Dwk_ig + Dfk_ig)/max(Ek_ig, 1.0e-6) ! org trunk
                     !DoverE_ig(k) = (1.0 - fac)*DoverE_ig(k) + fac*(Dwk_ig + Dfk_ig)/max(Ek_ig, 1.0e-6) ! TODO - TL CHECK - why not with relaxation anymore?                 
                     !
                     do itheta = 1, ntheta
                        !
                        R_ig(itheta) = oneoverdt*ee_ig(itheta, k) + cgprev_ig(itheta)*eeprev_ig(itheta)/ds(itheta, k) + srcig_local(itheta, k) !TL: new version
                        !
                     enddo
                     !
                     do itheta = 2, ntheta - 1
                        !
                        A_ig(itheta) = -ctheta_ig(itheta - 1, k)*oneover2dtheta
                        B_ig(itheta) = oneoverdt + cg_ig(k)/ds(itheta,k) + DoverE_ig(k)
                        C_ig(itheta) = ctheta_ig(itheta + 1, k)*oneover2dtheta
                        !
                     enddo
                     !
                     if (ctheta_ig(1,k)<0) then
                        A_ig(1) = 0.0
                        B_ig(1) = oneoverdt - ctheta_ig(1, k)/dtheta + cg_ig(k)/ds(1, k) + DoverE_ig(k)
                        C_ig(1) = ctheta_ig(2, k)/dtheta
                     else
                        A_ig(1)=0.0
                        B_ig(1)=1.0/dt + cg_ig(k)/ds(1, k) + DoverE_ig(k)
                        C_ig(1)=0.0
                     endif
                     !
                     if (ctheta_ig(ntheta, k)>0) then
                        A_ig(ntheta) = -ctheta_ig(ntheta - 1, k)/dtheta
                        B_ig(ntheta) = oneoverdt + ctheta_ig(ntheta, k)/dtheta + cg_ig(k)/ds(ntheta, k) + DoverE_ig(k)
                        C_ig(ntheta) = 0.0
                     else
                        A_ig(ntheta) = 0.0
                        B_ig(ntheta) = oneoverdt + cg_ig(k)/ds(ntheta, k) + DoverE_ig(k)
                        C_ig(ntheta) = 0.0
                     endif
                     !
                     ! Solve tridiagonal system per point
                     !
                     call solve_tridiag(A_ig, B_ig, C_ig, R_ig, ee_ig(:,k), ntheta)
                     ee_ig(:, k) = max(ee_ig(:, k), 0.0)
                     !  
                  else
                     !
                     ee_ig(:, k) = 0.0
                     !
                  endif
                  !
               endif
               !
            else
               !
               ee(:, k) = 0.0
               if (wind) then
                  aa(:,k)  = 0.0
               endif
               ee_ig(:, k) = 0.0
               !               
            endif 
            !
         endif
         !
         if (neumannconnected(k)/=0) then
            kn = neumannconnected(k)
            sinhkh(kn) = sinhkh(k)
            kwav(kn) = kwav(k)
            Hmx(kn) = Hmx(k)
            ee(:, kn) = ee(:, k)
            ee_ig(:, kn) = ee_ig(:, k) ! TL: Added Neumann option for IG            
            ctheta(:, kn) = ctheta(:, k)
            cg(kn) = cg(k)
            if (wind) then
               sig(kn) = sig(k)
               Tp(kn) = 2.0*pi/sig(k)
               WsorE(:,kn) = WsorE(:,k)
               WsorA(:,kn) = WsorA(:,k)
               aa(:,kn) = aa(:,k)
            endif
            Df(kn) = Df(k)
            Dw(kn) = Dw(k)
         endif
         !
      enddo
      !
      if (sweep==4) then
         !
         ! Check convergence after all 4 sweeps
         !
         do k = 1, no_nodes
            !
            dee     = ee(:, k) - eeold(:, k)
            diff(k) = maxval(abs(dee))
            !
            if (diff(k)/eemax<crit ) then
               ok(k) = 1
            endif   
            !
         enddo
         !
         ! Percentage of converged points
         !
         percok = sum(ok)/dble(no_nodes)*100.0
         eemax  = maxval(ee)
         !
         ! Relative maximum error
         !
         error = maxval(diff)/eemax
         !
         write(logstr,'(a,i6,a,f10.5,a,f7.2)')'   iteration ',iter/4 ,' error = ',error,'   %ok = ',percok
         call write_log(logstr, 0)         
         !
         if (error<crit) then
            write(logstr,'(a,i6,a,f10.5,a,f7.2)')'   converged at iteration ',iter/4 ,' error = ',error,'   %ok = ',percok
            call write_log(logstr, 0)            
            exit
         else
             if (iter == niter) then !made it to the end without reaching 'error<crit', still want output
                write(logstr,'(a,i6,a,f10.5,a,f7.2)')'    ended at iteration ',iter/4 ,' error = ',error,'   %ok = ',percok
                call write_log(logstr, 0)                
             endif
         endif
         !
      endif
      !
      ! write output every iteration > TL: removed in this implementation
      !
   enddo
   !
   do k=1,no_nodes
      !
      ! Compute directionally integrated parameters for output
      !
         if (depth(k)>1.1*hmin) then
            !
            E(k)      = sum(ee(:,k))*dtheta
            H(k)      = sqrt(8*E(k)/rho/g)
            thetam(k) = atan2(sum(ee(:, k)*sin(theta)),sum(ee(:, k)*cos(theta)))
            if (wind) then
               uorbi     = 0.5*sig(k)*Hk/sinhkh(k)
            else
               uorbi     = pi*H(k)/Tp(k)/sinhkh(k) !! to check if we want array
            endif
            Df(k)     = 0.28*rho*fw(k)*uorbi**3
            !                     
            if (wind) then
               call baldock(rho, g, alfa, gamma, depth(k), H(k), 2.0*pi/sig(k), 1, Dw(k), Hmx(k))
            else
               call baldock(rho, g, alfa, gamma, depth(k), H(k), Tp(k), 1, Dw(k), Hmx(k))                                   
            endif
            !
            if (vegetation) then
                !
                ! Compute wave dissipation due to vegetation
                call vegatt(sig(k), no_nodes, kwav(k), no_secveg, veg_ah(k,:), veg_bstems(k,:), veg_Nstems(k,:), veg_Cd(k,:), depth(k), rho, g, H(k), Dveg(k)) 
                !
                ! Compute the non-linear wave velocity time series (unl) using a wave shape model
                call swvegnonlin(no_nodes, kwav(k), depth(k), H(k), g, Tp(k), unl(k,:))
                ! NOTE - TODO: double check whether we want to call this in or outside the 'do k=1,no_nodes' loop!
                !
                ! Now also call 'momeqveg' to compute wave drag force due to vegetation
                call momeqveg(no_nodes, no_secveg, veg_ah(k,:), veg_bstems(k,:), veg_Nstems(k,:), veg_Cd(k,:), depth(k), rho, H(k), Tp(k), unl(k,:), Fvw(k))
                ! NOTE - TL: for now replaced 'Trep' by 'Tp(k)' 
                !
            else
                Dveg(k) = 0.
                Fvw(k) = 0.                
            endif
	        !
            !F(k) = Dw(k)*kwav(k)/sig(k)/rho/depth(k)
            !F(k) = (Dw(k) + Df(k))*kwav(k)/sig(k)/rho/depth(k) 	     
            F(k) = (Dw(k) + Dveg(k))*kwav(k)/sig(k)/rho/depth(k)
            F(k) = F(k) + Fvw(k)
	    !F(k) = (Dw(k) + Df(k))*kwav(k)/sigm ! TODO TL: before was this, now multiplied with rho*depth(k) in sfincs_snapwave.f90  
            !
            if (vegetation) then                
                if (Dveg(k) > 0.0) then
                    write(logstr,*)'k ',k,'depth(k)',depth(k),'H(k) ',H(k),'Dw(k) ', Dw(k),'Dveg(k) ', Dveg(k),'Hmx(k) ', Hmx(k),'kwav(k) ', kwav(k),'sig(k) ', sig(k), 'thetam(k)',thetam(k),'F(k) ',F(k)
                    call write_log(logstr, 0)  
                endif
            else
                if (H(k) > 0.0) then                
                    write(logstr,*)'k ',k,'depth(k)',depth(k),'H(k) ',H(k),'Dw(k) ', Dw(k),'Hmx(k) ', Hmx(k),'kwav(k) ', kwav(k),'sig(k) ', sig(k), 'thetam(k)',thetam(k),'F(k) ',F(k)
                    call write_log(logstr, 0)                      
                endif                
            endif             
            !
            if (igwaves) then
               !
               ! IG wave height
               !
               E_ig(k)      = sum(ee_ig(:,k))*dtheta
               H_ig(k)      = sqrt(8*E_ig(k)/rho/g)
               !
               Df_ig(k)     = fw_ig(k)*0.0361*(9.81/depth(k))**1.5*H(k)*E_ig(k) !org               
               !
               call baldock(rho, g, alfa_ig, gamma_ig, depth(k), H_ig(k), T_ig(k), 1, Dw_ig(k), Hmx_ig(k))               
               !
               ! average beta, alphaig, srcig over directions
               betamean(k) = sum(beta_local(:,k))/ntheta ! real mean
               !betamean(k)   = maxval(beta_local(:,k))
               !               
               alphaig(k) = sum(alphaig_local(:,k))/ntheta ! real mean 
               !
               srcig(k)   = sum(srcig_local(:,k)) /ntheta ! real mean                    
               !
            endif
            !
            if (wind) then
                Tp(k)  = 2.0*pi/sig(k)
                SwE(k) = sum(WsorE(:,k))*dtheta    
                SwA(k) = sum(WsorA(:,k))*dtheta    !(2*pi*sum(WsorA(:,k))*dtheta - 2*pi/sig(k)*SwE(k)) / E(k)
            endif
         endif
      !
   enddo
   callno=callno+1
   !
   end subroutine solve_energy_balance2Dstat


   subroutine solve_tridiag(a, b, c, d, x, n)
   !
   implicit none
   !	 a - sub-diagonal (means it is the diagonal below the main diagonal)
   !	 b - the main diagonal
   !	 c - sup-diagonal (means it is the diagonal above the main diagonal)
   !	 d - right part
   !	 x - the answer
   !	 n - number of equations

!   integer,parameter :: r8 = kind(1.)
   !
   integer,intent(in) :: n
   real*4,dimension(n),intent(in) :: a, b, c, d
   real*4,dimension(n),intent(out) :: x
   real*4,dimension(n) :: cp,dp
   real*4 :: m
   integer i
   !
   ! initialize c-prime and d-prime
   !
   cp(1) = c(1)/b(1)
   dp(1) = d(1)/b(1)
   ! solve for vectors c-prime and d-prime
   do i = 2, n
      m = b(i) - cp(i - 1)*a(i)
      cp(i) = c(i)/m
      dp(i) = (d(i) - dp(i - 1)*a(i))/m
   enddo
   ! initialize x
   x(n) = dp(n)
   ! solve for x from the vectors c-prime and d-prime
   do i = n-1, 1, -1
      x(i) = dp(i) - cp(i)*x(i + 1)
   end do
   !
   end subroutine solve_tridiag

   subroutine baldock (rho,g,alfa,gamma,depth,H,T,opt,Dw,Hmax)
   !
   real*4, intent(in)                :: rho
   real*4, intent(in)                :: g
   real*4, intent(in)                :: alfa
   real*4, intent(in)                :: gamma
   real*4, intent(in)                :: depth
   real*4, intent(in)                :: H
   real*4, intent(in)                :: T
   integer, intent(in)               :: opt
   real*4, intent(out)               :: Dw
   real*4, intent(in)                :: Hmax
   real*4                            :: Hloc
   !
   ! Compute dissipation according to Baldock
   !
   Hloc=max(H,1.e-6)
   if (opt==1) then
      Dw=0.28*alfa*rho*g/T*exp(-(Hmax/Hloc)**2)*(Hmax**2+Hloc**2)
      !Dw=0.25*alfa*rho*g/T*exp(-(Hmax/Hloc)**2)*(Hmax**2+Hloc**2) !Old
   else
      Dw=0.28*alfa*rho*g/T*exp(-(Hmax/Hloc)**2)*(Hmax**3+Hloc**3)/gamma/depth
      !Dw=0.25*alfa*rho*g/T*exp(-(Hmax/Hloc)**2)*(Hmax**3+Hloc**3)/gamma/depth !Old     
   endif
   !
   end subroutine baldock
   
   subroutine determine_infragravity_source_sink_term(inner, no_nodes, ntheta, w, ds, prev, cg_ig, nwav, depth, zb, H, ee, ee_ig, eeprev, eeprev_ig, cgprev, ig_opt, alphaigfac, alphaig_local, beta_local, srcig_local)
    !   
    implicit none
    !  
    ! Incoming variables
    logical, dimension(no_nodes), intent(in)         :: inner           ! mask of inner grid points (not on boundary)    
    integer, intent(in)                              :: no_nodes,ntheta ! number of grid points, number of directions  
    real*4,  dimension(2,ntheta,no_nodes),intent(in) :: w               ! weights of upwind grid points, 2 per grid point and per wave direction
    real*4, dimension(ntheta,no_nodes), intent(in)   :: ds              ! distance to interpolated upwind point, per grid point and direction   
    integer, dimension(2,ntheta,no_nodes),intent(in) :: prev            ! two upwind grid points per grid point and wave direction    
    real*4, dimension(no_nodes), intent(in)          :: cg_ig           ! group velocity
    real*4, dimension(no_nodes), intent(in)          :: nwav            ! wave number n  
    real*4, dimension(no_nodes), intent(in)          :: depth           ! water depth
    real*4, dimension(no_nodes), intent(in)          :: zb              ! actual bed level       
    real*4, dimension(no_nodes), intent(in)          :: H               ! wave height        
    real*4, dimension(ntheta,no_nodes), intent(in)   :: ee              ! energy density
    real*4, dimension(ntheta,no_nodes), intent(in)   :: ee_ig           ! energy density infragravity waves
    integer, intent(in)                              :: ig_opt          ! option of IG wave settings (1 = default = conservative shoaling based dSxx and Baldock breaking)    
    real*4, intent(in)                               :: alphaigfac      ! Multiplication factor for IG shoaling source/sink term, default = 1.0
    !
    ! Inout variables
    real*4, dimension(:,:), intent(inout)            :: alphaig_local   ! Local infragravity wave shoaling parameter alpha
    real*4, dimension(:,:), intent(inout)            :: srcig_local     ! Energy source/sink term because of IG wave shoaling
    real*4, dimension(:), intent(inout)              :: eeprev, cgprev  ! energy density and group velocity at upwind intersection point
    real*4, dimension(:), intent(inout)              :: eeprev_ig       ! energy density  at upwind intersection point 
    real*4, dimension(ntheta,no_nodes), intent(inout):: beta_local      ! Local bed slope based on bed level per direction              
    !
    ! Internal variables
    integer                                          :: itheta          ! directional counter
    integer                                          :: k               ! counters (k is grid index)    
    integer                                          :: k1,k2           ! upwind counters (k is grid index)
    real*4                                           :: gam             ! local gamma (Hinc / depth ratio)   
    real*4, dimension(ntheta,no_nodes)               :: depthprev       ! water depth at upwind intersection point         
    real*4, dimension(ntheta,no_nodes)               :: Sxx             ! Radiation Stress
    real*4, dimension(:), allocatable                :: Sxxprev         ! radiation stress at upwind intersection point  
    real*4, dimension(:), allocatable                :: Hprev           ! Incident wave height at upwind intersection point  
    real*4                                           :: dSxx            ! difference in Radiation stress
    real*4                                           :: Sxx_cons        ! conservative estimate of radiation stress using conservative shoaling     
    !   
    ! Allocate internal variables
    allocate(Sxxprev(ntheta))       
    allocate(Hprev(ntheta))  
    !
    Sxx = 0.0
    !
    do k = 1, no_nodes
        !
        if (inner(k)) then    
            !
            ! Compute exchange source term inc to ig waves - per direction      
            !
            do itheta = 1, ntheta
                !
                k1 = prev(1, itheta, k)
                k2 = prev(2, itheta, k)
                !
                if (k1>0 .and. k2>0) then ! IMPORTANT - for some reason (k1*k2)>0 is not reliable always, resulting in directions being uncorrectly skipped!!!    
                    !
                    ! First calculate upwind direction dependent variables
                    depthprev(itheta,k)     = w(1, itheta, k)*depth(k1) + w(2, itheta, k)*depth(k2)           
                    !            
                    beta_local(itheta,k)  = max((w(1, itheta, k)*(zb(k) - zb(k1)) + w(2, itheta, k)*(zb(k) - zb(k2)))/ds(itheta, k), 0.0)
                    !
                    ! Notes:
                    ! - use actual bed level now for slope, because depth changes because of wave setup/tide/surge
                    ! - in zb, depth is negative > therefore zb(k) minus zb(k1)
                    ! - beta=0 means a horizontal or decreasing slope > need alphaig=0 then in IG src/sink term
                    !
                    !betan_local(itheta,k) = (beta/sigm_ig)*sqrt(9.81/max(depth(k), hmin)) ! TL: in case in the future we would need the normalised bed slope again   
                    !
                    ! TL - Note: cg_ig = cg
                    cgprev(itheta)      = w(1, itheta, k)*cg_ig(k1) + w(2, itheta, k)*cg_ig(k2)
                    !              
                    Sxx(itheta,k1)      = ((2.0 * max(0.0,min(1.0,nwav(k1)))) - 0.5) * ee(itheta, k1) ! limit so value of nwav is between 0 and 1
                    Sxx(itheta,k2)      = ((2.0 * max(0.0,min(1.0,nwav(k2)))) - 0.5) * ee(itheta, k2) ! limit so value of nwav is between 0 and 1
                    !
                    Sxxprev(itheta)     = w(1, itheta, k)*Sxx(itheta,k1) + w(2, itheta, k)*Sxx(itheta,k2)
                    !
                    eeprev(itheta)      = w(1, itheta, k)*ee(itheta, k1) + w(2, itheta, k)*ee(itheta, k2)  
                    eeprev_ig(itheta)   = w(1, itheta, k)*ee_ig(itheta, k1) + w(2, itheta, k)*ee_ig(itheta, k2)      
                    !
                    Hprev(itheta)       = w(1, itheta, k)*H(k1) + w(2, itheta, k)*H(k2)     
                    !     
                    ! Determine relative waterdepth 'gam'
                    !
                    gam = max(0.5*(Hprev(itheta)/depthprev(itheta,k) + H(k)/depth(k)), 0.0) ! mean gamma over current and upwind point
                    !
                    ! Determine dSxx and IG source/sink term 'srcig'
                    !
                    if (ig_opt == 1 .or. ig_opt == 2) then 
                        !
                        ! Calculate shoaling parameter alpha_ig following Leijnse et al. (2024)
                        !  
                        call estimate_shoaling_parameter_alphaig(beta_local(itheta,k), gam, alphaig_local(itheta,k)) ! [input, input, output]
                        !                
                        ! Now calculate source term component
                        !         
                        ! Newest dSxx/dx based method, using estimate of Sxx(k) using conservative shoaling
                        if (Sxxprev(itheta)<=0.0) then 
                            !
                            srcig_local(itheta, k) = 0.0 !Avoid big jumps in dSxx that can happen if a upwind point is a boundary point with Hinc=0
                            !
                        else
                        !              
                        if (ig_opt == 1) then ! Option using conservative shoaling for dSxx/dx
                            !
                            ! Calculate Sxx based on conservative shoaling of upwind point's energy: 
                            ! Sxx_cons = E(i-1) * Cg(i-1) / Cg * (2 * n(i) - 0.5)
                            Sxx_cons = eeprev(itheta) * cgprev(itheta) / cg_ig(k) * ((2.0 * max(0.0,min(1.0,nwav(k)))) - 0.5)
                            ! Note - limit so value of nwav is between 0 and 1, and Sxx therefore doesn't become NaN for nwav=Infinite  
                            !
                            dSxx = Sxx_cons - Sxxprev(itheta)
                            !
                        elseif (ig_opt == 2) then ! Option taking actual difference for dSxx/dx
                            !
                            dSxx = Sxx(itheta,k) - Sxxprev(itheta)                        
                            !
                        endif
                        !
                        dSxx = max(dSxx, 0.0)
                        !
                        srcig_local(itheta, k) = alphaigfac * alphaig_local(itheta,k) * sqrt(eeprev_ig(itheta)) * cgprev(itheta) / depthprev(itheta,k) * dSxx / ds(itheta, k)
                        !                         
                        endif                      
                        !                                        
                    else  ! TL: option to add future parameterisations here for e.g. coral reef type coasts
                        !
                        srcig_local(itheta, k) = 0.0
                        !
                    endif
                    !                  
                    srcig_local(itheta, k)  = max(srcig_local(itheta, k), 0.0)
                    !   
                endif
                !
            enddo  
            !
        endif
        !
    enddo    
    !    
   end subroutine determine_infragravity_source_sink_term   
   
   subroutine estimate_shoaling_parameter_alphaig(beta, gam, alphaig)
   real*4, intent(in)                :: beta
   real*4, intent(in)                :: gam 
   real*4, intent(out)               :: alphaig
   !
   real*4                            :: beta1, beta2, beta3, beta4, beta5, beta6, beta7        
   !
   ! Estimate shoaling parameter alphaig - as in Leijnse et al. (2024)
   !
   ! Determine constants 
   beta1 = 0.016993
   beta2 = 0.5
   beta3 = 17.7104
   beta4 = 1
   beta5 = 0.7
   beta6 = 0.11841
   beta7 = 0.34037
   !
   if (beta <= 0.0) then
       !
       alphaig = 0.0
       !
   else ! positively increasing local bed slope beta
       !
       if (gam > 0.0 .and. gam < beta7) then !deep water
          !
          alphaig = exp(-beta3 * beta ** beta4) * ((beta5 - gam) * beta6 + (beta7 - gam) * (beta1 / beta ** beta2)) 
          !
       elseif (gam >= beta7) then ! shallow water - for gam>0.7 the fit automatically goes to 0
          !
          alphaig = exp(-beta3 * beta ** beta4) * (max(beta5 - gam, 0.0)) * beta6           
          !
       else ! for safety, but negative gamma should not occur
          ! 
          alphaig = 0.0
          ! 
       endif
       !
   endif
   !
   ! Limit alphaig between [0, 1] to prevent large overshoots in case of low gamma and very small beta
   !
   alphaig = max(alphaig, 0.0)
   alphaig = min(alphaig, 1.0)   
   !             
   end subroutine estimate_shoaling_parameter_alphaig 
   
   subroutine hpsort_eps_epw (n, ra, ind, eps)
   !---------------------------------------------------------------------
   ! sort an array ra(1:n) into ascending order using heapsort algorithm,
   ! and considering two elements being equal if their values differ
   ! for less than "eps".
   ! n is input, ra is replaced on output by its sorted rearrangement.
   ! create an index table (ind) by making an exchange in the index array
   ! whenever an exchange is made on the sorted data array (ra).
   ! in case of equal values in the data array (ra) the values in the
   ! index array (ind) are used to order the entries.
   ! if on input ind(1)  = 0 then indices are initialized in the routine,
   ! if on input ind(1) != 0 then indices are assumed to have been
   !                initialized before entering the routine and these
   !                indices are carried around during the sorting process
   !
   ! no work space needed !
   ! free us from machine-dependent sorting-routines !
   !
   ! adapted from Numerical Recipes pg. 329 (new edition)
   !
   ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
   ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
   !
   ! This file is distributed under the terms of the GNU General Public
   ! License. See the file `LICENSE' in the root directory of the
   ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
   !
   ! Adapted from flib/hpsort_eps
   !---------------------------------------------------------------------
   implicit none
   !-input/output variables
   integer, intent(in)   :: n
   real*4, intent(in)  :: eps
   integer :: ind (n)
   real*4 :: ra (n)
   !-local variables
   integer :: i, ir, j, l, iind
   real*4 :: rra
   !
   ! initialize index array
   IF (ind (1) .eq.0) then
      DO i = 1, n
         ind (i) = i
      ENDDO
   ENDIF
   ! nothing to order
   IF (n.lt.2) return
   ! initialize indices for hiring and retirement-promotion phase
   l = n / 2 + 1

   ir = n

   sorting: do

      ! still in hiring phase
      IF ( l .gt. 1 ) then
         l    = l - 1
         rra  = ra (l)
         iind = ind (l)
         ! in retirement-promotion phase.
      ELSE
         ! clear a space at the end of the array
         rra  = ra (ir)
         !
         iind = ind (ir)
         ! retire the top of the heap into it
         ra (ir) = ra (1)
         !
         ind (ir) = ind (1)
         ! decrease the size of the corporation
         ir = ir - 1
         ! done with the last promotion
         IF ( ir .eq. 1 ) then
            ! the least competent worker at all !
            ra (1)  = rra
            !
            ind (1) = iind
            exit sorting
         ENDIF
      ENDIF
      ! wheter in hiring or promotion phase, we
      i = l
      ! set up to place rra in its proper level
      j = l + l
      !
      DO while ( j .le. ir )
         IF ( j .lt. ir ) then
            ! compare to better underling
            IF ( hslt( ra (j),  ra (j + 1) ) ) then
               j = j + 1
               !else if ( .not. hslt( ra (j+1),  ra (j) ) ) then
               ! this means ra(j) == ra(j+1) within tolerance
               !  if (ind (j) .lt.ind (j + 1) ) j = j + 1
            ENDIF
         ENDIF
         ! demote rra
         IF ( hslt( rra, ra (j) ) ) then
            ra (i) = ra (j)
            ind (i) = ind (j)
            i = j
            j = j + j
            !else if ( .not. hslt ( ra(j) , rra ) ) then
            !this means rra == ra(j) within tolerance
            ! demote rra
            ! if (iind.lt.ind (j) ) then
            !    ra (i) = ra (j)
            !    ind (i) = ind (j)
            !    i = j
            !    j = j + j
            ! else
            ! set j to terminate do-while loop
            !    j = ir + 1
            ! endif
            ! this is the right place for rra
         ELSE
            ! set j to terminate do-while loop
            j = ir + 1
         ENDIF
      ENDDO
      ra (i) = rra
      ind (i) = iind

   END DO sorting

   contains

   !  internal function
   !  compare two real number and return the result

   logical function hslt( a, b )
      real*4 :: a, b
      IF( abs(a-b) <  eps ) then
         hslt = .false.
      ELSE
         hslt = ( a < b )
      end if
   end function hslt

   !
   end subroutine hpsort_eps_epw

   subroutine timer(t)
   real*4,intent(out)               :: t
   integer*4                        :: count,count_rate,count_max
   call system_clock (count,count_rate,count_max)
   t = real(count)/count_rate
   end subroutine timer

   subroutine vegatt(sigm, no_nodes, kwav, no_secveg, veg_ah, veg_bstems, veg_Nstems, veg_Cd, depth, rho, g, H, Dveg) 
        use snapwave_domain
		
        implicit none
        
		! declare variables
		real*4, intent(in)                               :: sigm            ! wave frequency (per cell)
        integer, intent(in)                              :: no_nodes        ! number of unstructured grid nodes
        integer, intent(in)                              :: no_secveg       ! number of sections in the vertical
        real*4, dimension(no_secveg), intent(in)         :: veg_ah          ! Height of vertical sections used in vegetation schematization [m wrt zb_ini (zb0)]  (per cell)
        real*4, dimension(no_secveg), intent(in)         :: veg_bstems      ! Width/diameter of individual vegetation stems [m] (per cell)
        real*4, dimension(no_secveg), intent(in)         :: veg_Nstems      ! Number of vegetation stems per unit horizontal area [m-2] (per cell)
        real*4, dimension(no_secveg), intent(in)         :: veg_Cd          ! Bulk drag coefficient [-] (per cell)
        real*4, intent(in)                               :: depth           ! bed level, water depth (per cell)
        real*4, intent(in)                               :: rho
        real*4, intent(in)                               :: g
        real*4, intent(in)                               :: H               ! wave height (per cell)
		real*4                                           :: Cdterm
        integer                                          :: m
        real*4, intent(in)                               :: kwav           ! wave number (per cell)
        real*4, intent(out)                              :: Dveg           ! dissipation by vegetation (per cell)
		
        ! Set dissipation in vegetation to zero everywhere for a start
		Dveg = 0.d0
        
        ! XB contains a check for porous in-canopy model -> not relevant for Snapwave?
		! ...
	
		! First compute drag coefficient (if not user-defined)

		if (no_secveg > 0) then ! only in case vegetation is present
			do m=1,no_secveg ! for each vertical vegetation section
				if (veg_Cd(m) < 0.d0) then ! If Cd is not user specified: call subroutine of M. Bendoni (see below)
                    !
					!call bulkdragcoeff(veg_ah(m),m,Cdterm,no_nodes,no_secveg,depth,H,kwav,veg_bstems(m),sigm) ! bulkdragcoeff(ahveg(k,m)+zb0(k)-zb(k),m,k,Cdterm) <- no bed level change implemented in Snapwave
					!write(logstr,*)'Cd is not user specified: using m. bendoni bulkdragcoefficient to compute cd: ',cdterm
                    !veg_Cd(m) = Cdterm
                    !                     
                    write(logstr,*)'SnapWave ERROR - Cd is not specified for layer: ',m    
                    call write_log(logstr, 0)
                    !
                    !
				endif
			enddo
		endif
		!
		! Attenuation by vegetation is computed in wave action balance (swvegatt) and the momentum balance (momeqveg);
		! 1) Short wave dissipation by vegetation
        call swvegatt(sigm, no_nodes, kwav, no_secveg, veg_ah, veg_bstems, veg_Nstems, veg_Cd, depth, rho, g, H, Dveg)
        !
		! 2) Mom.Eq.: Long wave dissipation, mean flow dissipation, nonlinear short wave effects, effect of emerged vegetation -> not implemented
		! call momeqveg(zb0) !< not relevant in Snapwave because no long waves and currents computed
		!
    end subroutine vegatt
    
    subroutine swvegatt(sigm, no_nodes, kwav, no_secveg, veg_ah, veg_bstems, veg_Nstems, veg_Cd, depth, rho, g, H, Dveg)! Short wave dissipation by vegetation
        !use snapwave_data
        !use snapwave_domain
        !
        implicit none
	    !
		! declare variables
		integer, intent(in)                             :: no_nodes        ! number of unstructured grid nodes
        integer, intent(in)                             :: no_secveg
		real*4,intent(in)                      	        :: sigm  ! 
        real*4, dimension(no_secveg), intent(in)        :: veg_ah          ! Height of vertical sections used in vegetation schematization [m wrt zb_ini (zb0)]
        real*4, dimension(no_secveg), intent(in)        :: veg_bstems      ! Width/diameter of individual vegetation stems [m]
        real*4, dimension(no_secveg), intent(in)        :: veg_Nstems      ! Number of vegetation stems per unit horizontal area [m-2]
        real*4, dimension(no_secveg), intent(in)        :: veg_Cd          ! Bulk drag coefficient [-]
        real*4, intent(in)                              :: depth           ! bed level, water depth
        real*4, intent(in)                              :: rho
        real*4, intent(in)                              :: g
        real*4, intent(in)                              :: H               ! wave height
		! 
		! local variables
		real*4                                      :: pi              ! 3.14159
        integer                                     :: k,m  ! indices of actual x,y point
        !
		real*4                                      :: aht,hterm,htermold,Dvgt,ahtold
		real*4            		                    :: Dvg,kmr!,kwav
        real*4, intent(in)                          :: kwav!,k
        !
        real*4, intent(out)                         :: Dveg
		!
		pi = 4.d0*atan(1.d0)
		kmr = min(max(kwav, 0.01d0), 100.d0)
		!
		! Set dissipation in vegetation to zero everywhere for a start
		Dvg = 0.d0
        Dvgt = 0.d0
        htermold = 0.d0
        ahtold = 0.d0
        if (no_secveg>0) then ! only if vegetation is present
            do m=1,no_secveg
	            !
                ! Determine height of vegetation section (restricted to current bed level)
                !aht = veg(ind)%ah(m)+ahtold !+s%zb0(k,j)-s%zb(k,j)!(max(veg(ind)%zv(m)+s%zb0(k,j),s%zb(k,j)))
                aht = veg_ah(m)+ahtold
	            ! 
                ! restrict vegetation height to local water depth
                aht = min(aht, depth)
	            !
                ! compute hterm based on ah
                hterm = (sinh(kmr*aht)**3+3*sinh(kmr*aht))/(3.d0*kmr*cosh(kmr* depth)**3) !
	            !
                ! compute dissipation based on aht and correct for lower elevated dissipation layers (following Suzuki et al. 2012)
                Dvgt = 0.5d0/sqrt(pi)*rho*veg_Cd(m)*veg_bstems(m)*veg_Nstems(m)*(0.5d0*kmr*g/sigm)**3*(hterm-htermold)*H**3
	            !
                ! save hterm to htermold to correct possibly in next vegetation section
                htermold = hterm
                ahtold   = aht
	            !
                ! add dissipation current vegetation section
                Dvg = Dvg + Dvgt
            enddo
        endif
		Dveg = Dvg
    end subroutine swvegatt
    
	subroutine bulkdragcoeff(ahh, m, Cdterm, no_nodes, no_secveg, depth, H, kwav, veg_bstems, sigm)!(ahh,m,i,Cdterm)
		!    Michele Bendoni: subroutine to calculate bulk drag coefficient for short wave
		!    energy dissipation based on the Keulegan-Carpenter number (adapted from XBeach)
		!    Ozeren et al. (2013) or Mendez and Losada (2004)	 
        !
		implicit none
	    !
		real*4,  intent(out)                            :: Cdterm
		real*4,  intent(in)                             :: ahh              ! [m] plant (total) height
		integer, intent(in)                             :: m
        integer, intent(in)                             :: no_nodes         ! number of unstructured grid nodes
        integer, intent(in)                             :: no_secveg
        real*4, intent(in)                              :: depth            ! bed level, water depth
        real*4, intent(in)                              :: H                ! wave height
		real*4, intent(in)                              :: kwav             ! wave number
        real*4, intent(in)                              :: veg_bstems       ! Width/diameter of individual vegetation stems [m]
        real*4, intent(in)                              :: sigm             ! [rad/s] mean frequency
	    !	
		! Local variables
		real*4                                          :: pi                        ! 3.14159
        real*4                                          :: alfav            ! [-] ratio between plant height and water depth
		real*4                                          :: um               ! [m/s] typical velocity acting on the plant
		real*4                                          :: Tp               ! [s] reference wave period
		real*4                                          :: KC               ! [-] Keulegan-Carpenter number
		real*4                                          :: Q                ! [-] modified Keulegan-Carpenter number
		integer                                         :: myflag           ! 1 => Ozeren et al. (2013); 2 => Mendez and Losada (2004)
		!
		myflag = 2
        pi = 4.d0*atan(1.d0)
		!
		! Representative wave period
		Tp = 2*pi/sigm 
		!
		! Coefficient alfa
		if (ahh>=depth) then
			alfav = 1.d0
		else
			alfav = ahh/depth
		endif
		!
		! Representative orbital velocity
		! (Could we also use urms here?)
		um = 0.5d0*H*sigm*cosh(kwav*alfav*depth)/sinh(kwav*depth) 
		!
		! Keulegan-Carpenter number
		KC = um*Tp/veg_bstems
		!
		! Bulk drag coefficient
		if (myflag == 1) then
			!
			! Approach from Ozeren et al. (2013), eq?
			!
			if (KC>=10.d0) then
				Cdterm = 0.036d0+50.d0/(KC**0.926d0)
			else
				Cdterm = 0.036d0+50.d0/(10.d0**0.926d0)
			endif
		elseif (myflag == 2) then
			!
			! Approach from Mendez and Losada (2004), eq. 40
			! Only applicable for Laminaria Hyperborea (kelp)???
			!
			Q = KC/(alfav**0.76d0)
			if (Q>=7) then
				Cdterm = exp(-0.0138*Q)/(Q**0.3d0)
			else
				Cdterm = exp(-0.0138*7)/(7**0.3d0)
			endif
		endif
		!
    end subroutine bulkdragcoeff
    
subroutine momeqveg(no_nodes, no_secveg, veg_ah, veg_bstems, veg_Nstems, veg_Cd, depth, rho, H, Trep, unl, Fvw)
    !
    implicit none
    !
    ! Inputs
    integer, intent(in) :: no_nodes, no_secveg
    real*4, intent(in) :: depth ,rho, H, Trep
    real*4, dimension(no_secveg), intent(in) :: veg_ah, veg_bstems, veg_Nstems, veg_Cd
    real*4, dimension(50), intent(in) :: unl
    !
    ! Output
    real*4, intent(out) :: Fvw
    !
    ! Local variables
    integer :: m, t
    real*4 :: dt, hvegeff, Fvgnlt, integral
    real*4 :: Cd, b, N
    !
    ! Initialize output force
    !
    Fvw = 0.0
    !
    ! Time step within wave period
    !
    dt = Trep / 50.0
    !
    ! Loop over vertical vegetation sections
    do m = 1 , no_secveg
        ! Effective submerged height of vegetation section
        hvegeff = min(veg_ah(m), depth)
        ! Read vegetation parameters
        Cd = veg_Cd(m)
        b = veg_bstems(m)
        N = veg_Nstems(m)
        ! Integrate vegetation drag over wave period using unl
        integral = 0.0
        do t = 1, 50 !50=PPWL
            integral = integral + (0.5 * Cd * b * N * hvegeff * unl(t) * abs(unl(t) ) ) * dt
        enddo
        ! Convert to force per unit mass and sum
        Fvgnlt = integral / depth / rho
        Fvw = Fvw + Fvgnlt
    enddo
end subroutine momeqveg    
   
subroutine swvegnonlin(no_nodes, kwav, depth, H, g, Trep, unl)
    ! input= no_nodes, kwav(k), H(k), depth(k), g, Tp(k), unl(k,:)
    !
    ! Based on Deltares' XBeach SurfBeat' subroutine: swvegnonlin
    !
    implicit none
    !
    integer :: no_nodes, k
    integer :: irf, ih0, it0, jrf, ih1, it1
    integer , save :: nh , nt !TL: NOTE - NOT familiar with THIS_IMAGE 'save' statement, for now keep
    real*4 :: p ,q , f0 , f1 , f2 , f3
    real*4, save :: dh , dt
    real*4, dimension(no_nodes) :: kmr , Urs , phi , w1 , w2
    real*4, dimension(8) , save :: urf0
    real*4, dimension(50) , save :: urf2 , urf
    real*4, dimension(50 ,8), save :: cs , sn , urf1
    real*4, dimension(:), save , allocatable :: h0, t0
    real*4, dimension(no_nodes, 50),intent(out) :: unl ! NOTE - TL: we don't use 'etaw0' in the end?
    !
    real*4  :: pi = 4.*atan(1.0)   
    real*4, intent(in) :: kwav, depth, g, H, Trep ! depth = the 'hh' of XBeach  
    real*4, dimension(:,:,:), allocatable :: RFveg
    !
    ! Compute net drag force due to wave skewness based on Rienecker & Fenton (1981)
    !
    ! load Ad's RF-table (update for depth averaged velocities?)
    !include 'RFveg.inc' 
    !include 'RFtable.inp'
    !
    ! Prepare interpolation of RF table
    if (.not. allocated(h0)) then
        allocate(h0(no_nodes))
        allocate(t0(no_nodes))
        allocate(RFveg(no_nodes))        
        dh = 0.0
        dt = 1.25
        nh = floor(0.54/ dh)
        nt = floor(25 / dt )
        do irf =1 ,8
            do jrf =1 ,50
                cs ( jrf , irf ) = cos (( jrf * 2 * pi / 50) * irf )
                sn ( jrf , irf ) = sin (( jrf * 2 * pi / 50) * irf )
            enddo
        enddo
    endif    
    !
    h0 = min(nh * dh, max(dh, min(H, depth) / depth) )
    t0 = min(nt * dt, max(dt, Trep * sqrt (g / depth) ) )
    !
    ! Initialize
    urf0 = 0
    urf1 = 0
    urf2 = 0
    urf = 0
    w1 = 0
    w2 = 0
    phi = 0
    Urs = 0
    kmr = 0
    !
    ! Compute phase and weights for Ruessink wave shape
    kmr = min(max(kwav, 0.01), 100.0)
    Urs = H / (kmr * kmr * (depth **3) )
    phi = pi /2 * (1 - tanh (0.815/(Urs **0.672) ) )
    w1 = 1 - phi /( pi /2)
    w2 = 1 - w1    
    !
    ! Interpolate RF table and compute velocity profiles
    do k =1, no_nodes
        !
        ih0 = floor( h0(k) / dh)
        it0 = floor( t0(k) / dt)
        ih1 = min(ih0 + 1, nh)
        it1 = min(it0 + 1, nt)
        p = ( h0(k) - ih0 * dh) / dh
        q = ( t0(k) - it0 * dt) / dt
        f0 = (1 - p) * (1 - q)
        f1 = p * (1 - q)
        f2 = q * (1 - p)
        f3 = p * q
        !
        do irf = 1, 8
            urf0(irf) = f0 * RFveg(irf + 3, ih0, it0) + f1 * RFveg(irf + 3, ih1, it0) + f2 * RFveg(irf+3, ih0, it1) + f3 * RFveg(irf + 3, ih1, it1)
        enddo    
        !
        do irf = 1, 8
            urf1(:, irf) = urf0(irf)
        enddo
        !
        urf1 = urf1 * (w1(k) * cs + w2(k) * sn )
        urf2 = sum(urf1, 2)
        unl(k,:) = urf2 * sqrt(g * depth )
        !etaw0(k,:) = unl0 (i ,j ,:) * sqrt (max( depth(k ) ,0 ) / g ) #TL: not used
    enddo   
    !
end subroutine swvegnonlin

end module snapwave_solver 
