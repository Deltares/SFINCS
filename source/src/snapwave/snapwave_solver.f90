module snapwave_solver
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! SnapWave solver: stationary 2D wave energy balance for incident and infragravity waves,
   ! including breaking, bottom friction, wind growth and vegetation. Infragravity source
   ! terms are in snapwave_infragravity, vegetation routines in snapwave_vegetation.
   !
   ! Subroutines / functions:
   !
   !   compute_wave_field()
   !     Compute the wave field for the current conditions: celerities, energy balance and
   !     derived output.
   !     Called from: compute_snapwave (sfincs_snapwave)
   !
   !   solve_energy_balance2Dstat(...)
   !     Solve the stationary 2D wave energy (and wave action) balance by sweeping over the
   !     nodes.
   !     Called from: compute_wave_field (snapwave_solver)
   !
   !   solve_tridiag(a, b, c, d, x, n)
   !     Solve a tridiagonal system of equations.
   !     Called from: solve_energy_balance2Dstat (snapwave_solver)
   !
   !   baldock(rho, g, alfa, gamma, depth, H, T, iexp, Dw, Hmax)
   !     Wave breaking dissipation according to Baldock et al. (1998).
   !     Called from: solve_energy_balance2Dstat (snapwave_solver)
   !
   !   timer(t)
   !     Return the wall-clock time in seconds from system_clock.
   !     Called from: compute_wave_field (snapwave_solver)
   !
   ! ------------ Copyright notice SnapWave ------------
   !
   ! Copyright (C)  Dano Roelvink, IHE Delft, 2020-2026.
   !
   ! This file is part of SnapWave.
   !
   ! SnapWave is free software: you can redistribute it and/or modify
   ! it under the terms of the GNU Lesser General Public License as
   ! published by the Free Software Foundation version 2.1.
   !
   ! SnapWave is distributed in the hope that it will be useful,
   ! but WITHOUT ANY WARRANTY; without even the implied warranty of
   ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   ! GNU Lesser General Public License for more details.
   !
   ! You should have received a copy of the GNU Lesser General Public
   ! License along with SnapWave.  If not, see <http://www.gnu.org/licenses/>.
   !
   ! This file has been modified by Stichting Deltares for coupling
   ! with SFINCS. The modifications are distributed under the same
   ! licence (GNU Lesser General Public License version 2.1).
   !
   ! Modifications Copyright (C)  Stichting Deltares, 2021-2026.
   !
   ! contact: software@deltares.nl
   ! Stichting Deltares
   ! P.O. Box 177
   ! 2600 MH Delft, The Netherlands
   !
   ! ------------ Copyright notice SnapWave ------------
   !
   
   use sfincs_log
    
   implicit none

contains
   
   subroutine compute_wave_field()
      !
      ! Compute the wave field for the current conditions: celerities, energy balance and
      ! derived output.
      !
      ! Called from: compute_snapwave (sfincs_snapwave)
      !
      use snapwave_data
      !
      implicit none
      !
      !real*8, intent(in)   :: time > TL: not used in this implementation
      !
      real*4, parameter    :: waveps = 1e-5
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
      pi    = 4 * atan(1.0)
      !
      call timer(t0)
      !
      if (.not.restart) then
         !
         ! Set energies to 0.0; note that boundary values have been set in update_boundaries
         !
         !$omp parallel do schedule(static)
         do k = 1, no_nodes
            if (inner(k)) then
               ee(:, k) = waveps
            endif
         enddo
         !$omp end parallel do
         !
         ee_ig = waveps
         !
      endif
      !
      ! Initialize wave period
      !
      !$omp parallel do schedule(static)
      do k = 1, no_nodes
         !
         if (inner(k)) then
            !
            Tp(k) = Tpini
            !
         endif
         !
         if (neumannconnected(k) > 0) then
            !
            Tp(neumannconnected(k)) = Tpini
            !
         endif
         !
      enddo
      !$omp end parallel do
      !
      ! Compute celerities and refraction speed
      !
      Tp      = max(tpmean_bwv, Tpini)   ! to check voor windgroei
      sig     = 2.0 * pi / Tp
      Tp_ig   = tpmean_bwv_ig! TL: now determined in snapwave_boundaries.f90 instead of Tinc2ig*Tp
      sigm_ig = 2.0 * pi / Tp_ig !TODO - TL: Question do we want Tp_ig now as contant, or also spatially varying like Tp ?
      expon   = - (sig * sqrt(depth / g))**2.5
      kwav    = sig**2 / g * (1.0 - exp(expon))**-0.4
      C       = sig / kwav
      nwav    = 0.5 + kwav * depth / sinh(min(2 * kwav * depth, 50.0))
      Cg      = nwav * C
      !
      if (igwaves) then
         !
         cg_ig   = Cg
         expon   = -(sigm_ig * sqrt(depth / g))**2.5
         kwav_ig = sigm_ig**2 / g * (1.0 - exp(expon))**-0.4
         !
      else
         !
         cg_ig   = 0.0
         kwav_ig = 0.0 
         !
      endif
      !
      ! Set Hmx and sinh(kh) for regular waves, and Hmx_ig for IG waves. Note that we use the same gamma for regular and IG waves, but this can be easily changed if needed.
      !
      !$omp parallel do schedule(static)
      do k = 1, no_nodes
         !
         sinhkh(k) = sinh(min(kwav(k) * depth(k), 50.0))
         !Hmx(k)    = 0.88 / kwav(k) * tanh(gamma * kwav(k) * depth(k) / 0.88)
         Hmx(k)    = gamma * depth(k)
         !
         if (igwaves) then
            !
            ! Why is this different from Hmx for regular waves where we use gamma * h?
            !
            Hmx_ig(k) = 0.88 / kwav_ig(k) * tanh(gamma_ig * kwav_ig(k) * depth(k) / 0.88) ! Note - uses gamma_ig
            !Hmx_ig(k)    = gamma_ig * depth(k)
            !
         endif
         !
      enddo
      !$omp end parallel do
      !   
      do itheta = 1, ntheta
         !
         ctheta(itheta,:)    = sig / sinh(min(2 * kwav * depth, 50.0)) * (dhdx * sin(theta(itheta)) - dhdy * cos(theta(itheta)))
         !
      enddo
      !
      if (igwaves) then
         !
         do itheta = 1, ntheta
            ctheta_ig(itheta,:) = sigm_ig / sinh(min(2 * kwav_ig * depth, 50.0)) * (dhdx * sin(theta(itheta)) - dhdy * cos(theta(itheta)))
         enddo
         !
      else
         !
         ctheta_ig = 0.0
         !
      endif
      !
      ! Limit unrealistic refraction speed to 1/2 pi per wave period
      !
      !$omp parallel do schedule(static)
      do k = 1, no_nodes
         ! 
         ctheta(:,k)    = sign(1.0, ctheta(:,k)) * min(abs(ctheta(:, k)), sig(k) / 4)
         !
      enddo
      !$omp end parallel do
      !
      if (igwaves) then
         !$omp parallel do schedule(static)
         do k=1, no_nodes
            ! 
            ctheta_ig(:,k) = sign(1.0, ctheta_ig(:,k)) * min(abs(ctheta_ig(:, k)), sigm_ig(k) / 4)
            !
         enddo
         !$omp end parallel do
      endif
      !
      ! Solve the directional wave energy balance on an unstructured grid
      !
      call timer(t2)
      !
      call solve_energy_balance2Dstat (x,y,dhdx, dhdy, no_nodes,inner, &
                                         w, ds, prev, &
                                         neumannconnected, &
                                         theta,ntheta,thetamean, &
                                         depth,kwav,cg,ctheta,fw, &
                                         Tp,Tp_ig,dt,rho,alpha,gamma, gammax, &
                                         wind, &
                                         H,Dw,F,Df,thetam,sinhkh, &
                                         Hmx, ee, windspreadfac, u10, niter, crit, &
                                         hmin, baldock_ratio, baldock_ratio_ig, baldock_exponent, &
                                         aa, sig, jadcgdx, sigmin, sigmax, &
                                         DoverE, relax_factor_DoverE, relax_factor_DoverA, &             
                                         c_dispT, WsorE, WsorA, SwE, SwA, Tpini, &
                                         igwaves, kwav_ig, cg_ig,H_ig,ctheta_ig,Hmx_ig, ee_ig,fw_ig, &
                                         beta, srcig, alphaig, Dw_ig, Df_ig, &
                                         vegetation, no_secveg, veg_ah, veg_bstems, veg_Nstems, veg_Cd, Dveg, &
                                         zb, nwav, ig_opt, alpha_ig, gamma_ig, gamma_fac_br, eeinc2ig, Tinc2ig, alphaigfac, shinc2ig, iterative_srcig)
      !
      call timer(t3)
      !
      Fx = F * cos(thetam)
      Fy = F * sin(thetam)
      !
   end subroutine
   
   
   
   subroutine solve_energy_balance2Dstat(x,y,dhdx, dhdy, no_nodes,inner, &
                                         w, ds, prev, &
                                         neumannconnected, &
                                         theta,ntheta,thetamean, &
                                         depth,kwav,cg,ctheta,fw, &
                                         Tp,T_ig,dt,rho,alfa,gamma, gammax, &
                                         wind, &
                                         H,Dw,F,Df,thetam,sinhkh, &
                                         Hmx, ee, windspreadfac, u10, niter, crit, &
                                         hmin, baldock_ratio, baldock_ratio_ig, baldock_exponent, &       
                                         aa, sig, jadcgdx, sigmin, sigmax, &
                                         DoverE, relax_factor_DoverE, relax_factor_DoverA, &   
                                         c_dispT, WsorE, WsorA, SwE, SwA, Tpini, &
                                         igwaves,kwav_ig, cg_ig,H_ig,ctheta_ig,Hmx_ig, ee_ig,fw_ig, &
                                         betamean, srcig, alphaig, Dw_ig, Df_ig, &       
                                         vegetation, no_secveg, veg_ah, veg_bstems, veg_Nstems, veg_Cd, Dveg, &
                                         zb, nwav, ig_opt, alfa_ig, gamma_ig, gamma_fac_br, eeinc2ig, Tinc2ig, alphaigfac, shinc2ig, iterative_srcig)
   !
   ! Solve the stationary 2D wave energy (and wave action) balance by sweeping over the
   ! nodes.
   !
   ! Input:
   !   x                   : x coordinates of the nodes
   !   y                   : y coordinates of the nodes
   !   dhdx                : bed slope in x direction
   !   dhdy                : bed slope in y direction
   !   no_nodes            : number of grid nodes
   !   w                   : weights of the two upwind nodes, per direction and node
   !   ds                  : distance to the upwind point, per direction and node
   !   prev                : indices of the two upwind nodes, per direction and node
   !   neumannconnected    : Neumann boundary node connected to each inner node (0 if none)
   !   theta               : wave directions [rad]
   !   ntheta              : number of wave directions
   !   thetamean           : mean offshore wave direction [rad]
   !   depth               : water depth [m]
   !   fw                  : wave friction factor [-]
   !   T_ig                : infragravity wave period [s]
   !   dt                  : time step [s]
   !   rho                 : water density [kg/m3]
   !   alfa                : coefficient alpha of the Baldock breaking model [-]
   !   gamma               : breaker index gamma [-]
   !   gammax              : maximum ratio of wave height to water depth [-]
   !   wind                : wind growth on or off
   !   windspreadfac       : directional distribution of the wind input [-]
   !   u10                 : wind speed at 10 m [m/s]
   !   niter               : maximum number of iterations
   !   crit                : relative accuracy of the stopping criterion
   !   hmin                : minimum water depth [m]
   !   baldock_ratio       : ratio H/Hmax above which incident waves break
   !   baldock_ratio_ig    : ratio H/Hmax above which infragravity waves break
   !   baldock_exponent    : exponent of the extra factor on the Baldock dissipation
   !   jadcgdx             : 1 to include the group velocity gradient in the wind input, 0
   !                         not to
   !   sigmin              : minimum wave frequency [rad/s]
   !   sigmax              : maximum wave frequency [rad/s]
   !   relax_factor_DoverE : under-relaxation factor for D/E (1.0 disables)
   !   relax_factor_DoverA : under-relaxation factor for D/A (1.0 disables)
   !   c_dispT             : coefficient of the wave period dissipation
   !   Tpini               : initial peak period [s]
   !   igwaves             : infragravity waves on or off
   !   kwav_ig             : infragravity wave number [rad/m]
   !   cg_ig               : infragravity group velocity [m/s]
   !   ctheta_ig           : infragravity refraction speed [rad/s]
   !   Hmx_ig              : maximum infragravity wave height [m]
   !   fw_ig               : infragravity friction factor [-]
   !   vegetation          : vegetation on or off
   !   no_secveg           : number of vertical vegetation sections
   !   veg_ah              : height of the vertical vegetation sections [m]
   !   veg_bstems          : stem diameter per section [m]
   !   veg_Nstems          : number of stems per unit area, per section [1/m2]
   !   veg_Cd              : bulk drag coefficient per section [-]
   !   zb                  : bed level [m]
   !   nwav                : ratio of group to phase velocity n [-]
   !   ig_opt              : infragravity wave option
   !   alfa_ig             : coefficient alpha of the Baldock model for infragravity waves
   !                         [-]
   !   gamma_ig            : breaker index for infragravity waves [-]
   !   gamma_fac_br        : factor on gamma for the maximum incident wave height in the
   !                         shoaling term
   !   eeinc2ig            : ratio of infragravity to incident energy at the boundary
   !   Tinc2ig             : ratio of infragravity to incident period
   !   alphaigfac          : factor on the infragravity shoaling source term
   !   shinc2ig            : fraction of the infragravity source term that is taken from the
   !                         incident waves
   !   iterative_srcig     : compute the infragravity source term implicitly within the
   !                         iterations
   !
   ! Input/output:
   !   inner               : mask of the inner nodes (not on a boundary)
   !   kwav                : wave number [rad/m]
   !   cg                  : group velocity [m/s]
   !   ctheta              : refraction speed, per direction and node [rad/s]
   !   Tp                  : peak wave period [s]
   !   sinhkh              : sinh(k h)
   !   Hmx                 : maximum wave height [m]
   !   ee                  : wave energy density, per direction and node
   !   aa                  : wave action density, per direction and node
   !   sig                 : mean wave frequency [rad/s]
   !   DoverE              : ratio of dissipation to energy D/E
   !   ee_ig               : infragravity energy density, per direction and node
   !
   ! Output:
   !   H                   : wave height [m]
   !   Dw                  : wave breaking dissipation [W/m2]
   !   F                   : wave force Dw/C/rho/h
   !   Df                  : bottom friction dissipation [W/m2]
   !   thetam              : mean wave direction [rad]
   !   WsorE               : wind input of wave energy, per direction and node
   !   WsorA               : wind input of wave action, per direction and node
   !   SwE                 : directionally integrated wind input of energy
   !   SwA                 : directionally integrated wind input of action
   !   H_ig                : infragravity wave height [m]
   !   betamean            : mean local bed slope parameter
   !   srcig               : directionally averaged infragravity source term
   !   alphaig             : mean infragravity shoaling parameter
   !   Dw_ig               : infragravity breaking dissipation [W/m2]
   !   Df_ig               : infragravity friction dissipation [W/m2]
   !   Dveg                : wave dissipation by vegetation [W/m2]
   !
   ! Called from: compute_wave_field (snapwave_solver)
   !
   use snapwave_windsource
   use snapwave_infragravity, only: determine_infragravity_source_sink_term
   use snapwave_vegetation, only: vegatt, momeqveg, swvegnonlin
   use hpsort_eps, only: hpsort_eps_epw
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
   real*4, dimension(no_nodes), intent(out)         :: H                      ! wave height
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
   integer,                     intent(in)          :: ig_opt                 ! option of IG wave settings (1 = default = conservative shoaling based dSxx as in Leijnse et al. 2024)  
   real*4,                      intent(in)          :: relax_factor_DoverA    ! underrelaxation factor for DoverA (set to 1.0 to disable)
   real*4,                      intent(in)          :: relax_factor_DoverE    ! underrelaxation factor for DoverE (set to 1.0 to disable)   
   !
   ! wind source vars
   !
   integer, intent(in)                                 :: jadcgdx               ! logical yes/no               ! logical yes/no
   real*4, intent(in)                                  :: sigmin, sigmax, c_dispT
   real*4, dimension(ntheta, no_nodes), intent(in)     :: windspreadfac        !< [-] distribution array for wind input
   real*4, dimension(ntheta,no_nodes),  intent(inout)  :: aa    
   real*4, dimension(no_nodes),         intent(inout)  :: DoverE
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
   real*4, dimension(no_nodes,50)                       :: etaw0                    ! non-linear sea surface time series, in 50 points per wave length
   !
   !
   ! Local variables and arrays
   !
   integer, dimension(:), allocatable         :: ok                     ! mask for fully iterated points
   integer, dimension(:), allocatable         :: ok_ig                  ! mask for fully iterated IG points   
   real*4                                     :: dtheta                 ! directional resolution
   real*4                                     :: eemax                  ! maximum wave energy density   
   real*4                                     :: eemax_ig               ! maximum IG wave energy density
   real*4                                     :: uorbi
   integer                                    :: sweep,iter            ! sweep number, number of iterations
   integer                                    :: k,k1,k2,count,kn,itheta ! counters (k is grid index)
   integer, dimension(:,:), allocatable       :: indx                   ! index for grid sorted per sweep direction
   real*4, dimension(:,:), allocatable        :: eeold               ! wave energy density, energy density previous iteration
   real*4, dimension(:,:), allocatable        :: eeold_ig            ! IG wave energy density, energy density previous iteration
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
   real*4, dimension(:), allocatable          :: DoverA                 ! ratio of mean wave dissipation over mean wave energy
   real*4, dimension(:), allocatable          :: DoverE_ig              ! ratio of mean wave dissipation over mean wave energy
   real*4, dimension(:), allocatable          :: E                      ! mean wave energy
   real*4, dimension(:), allocatable          :: E_ig                   ! mean wave energy
   real*4, dimension(:), allocatable          :: diff                   ! maximum difference of wave energy relative to previous iteration
   real*4, dimension(:), allocatable          :: diff_ig                ! maximum difference of IG wave energy relative to previous iteration
   real*4, dimension(:), allocatable          :: ra                     ! coordinate in sweep direction
   real*4, dimension(:), allocatable          :: sigm_ig
   integer, dimension(4)                      :: shift
   real*4                                     :: pi = 4.0 * atan(1.0)
   real*4                                     :: g = 9.81
   real*4                                     :: hmin                   ! minimum water depth! TL: make user changeable also here according to 'snapwave_hmin' in sfincs.inp   
   real*4                                     :: oneoverdt
   real*4                                     :: oneover2dtheta
   real*4                                     :: rhog8
   real*4                                     :: Dfk
   real*4                                     :: Dwk
   real*4                                     :: Ek
   real*4                                     :: Hk
   real*4                                     :: percok
   real*4                                     :: percok_ig              ! percentage of converged IG points   
   real*4                                     :: error
   real*4                                     :: error_ig               ! relative maximum IG wave error
   real*4                                     :: Dfk_ig
   real*4                                     :: Dwk_ig
   real*4                                     :: Ek_ig
   real*4                                     :: Hk_ig
   real*4                                     :: alfa_ig,gamma_ig  ! coefficients in Baldock wave breaking dissipation model for IG waves
   real*4                                     :: gamma_fac_br    ! factor times gamma that is used to determine the maximum incident wave breaking point in the surf zone using local incident wave height over water depth ratio, among others used to set the IG source term to 0 shallower than this point   
   real*4                                     :: eeinc2ig          ! ratio of incident wave energy as first estimate of IG wave energy at boundary
   real*4                                     :: Tinc2ig           ! ratio compared to period Tinc to estimate Tig
   real*4                                     :: alphaigfac        ! Multiplication factor for IG shoaling source/sink term, default = 1.0
   real*4                                     :: shinc2ig          ! Ratio of how much of the calculated IG wave source term, is subtracted from the incident wave energy (0-1, 0=default)   
   integer, save                              :: callno=1
   integer, intent(in)                        :: baldock_exponent  ! Exponent for multiplying the Baldock dissipation with a factor 'f = (Hloc / Hmax)**iexp' to enhance breaking when H > Hmax, with iexp = 0 (default, means unused), 1 or 2 
   !
   real*4, dimension(ntheta)                  :: sinth, costh            ! distribution of wave angles and offshore wave energy density   
   !
   ! Local wind source vars
   !
   real*4                                     :: Ak
   real*4                                     :: DwT
   real*4                                     :: DwAk
   real*4                                     :: ndissip    
   real*4                                     :: depthlimfac
   real*4                                     :: waveps=0.0001
   !
   ! Allocate local arrays
   !
   waveps = 0.0001
   !
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
      allocate(sigm_ig(no_nodes)); sigm_ig=0.0
      allocate(depthprev(ntheta,no_nodes)); depthprev=0.0
      allocate(beta_local(ntheta,no_nodes)); beta_local=0.0
      allocate(alphaig_local(ntheta,no_nodes)); alphaig_local=0.0
      allocate(eeold_ig(ntheta,no_nodes)); eeold_ig=0.0
      allocate(diff_ig(no_nodes)); diff_ig=0.0
      allocate(ok_ig(no_nodes)); ok_ig=0
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
   Df   = 0.0
   Dw   = 0.0
   Dveg = 0.0
   F    = 0.0
   !
   ok     = 0
   indx   = 0
   eemax  = maxval(ee)
   if (igwaves) then
      ok_ig   = 0
      eemax_ig = maxval(ee_ig)
   endif
   dtheta = theta(2) - theta(1)
   !
   if (dtheta < 0.0) dtheta = dtheta + 2*pi
   !
   if (wind) then
      !
      sig = 2 * pi / Tpini
      !
   else
      !
      sig = 2 * pi / Tp
      !
   endif
   !
   oneoverdt      = 1.0 / dt
   oneover2dtheta = 1.0 / 2.0 / dtheta
   rhog8          = 0.125 * rho * g
   thetam         = 0.0
   !H              = 0.0 ! TODO - TL: CHeck > needed for restart for IG > set to 0 now in snapwave_domain.f90
   Dveg           = 0.0
   Fvw            = 0.0
   unl            = 0.0
   etaw0          = 0.0
   !
   if (igwaves) then
      !
      sigm_ig        = 2 * pi / T_ig
      DoverE_ig      = 0.0
      !
   endif
   !
   if (wind) then
      !
      DoverA = 0.0
      ndissip = 3.0
      WsorE = 0.0
      WsorA = 0.0
      Ak = waveps / sigmax
      !
      ! Re-initialise Tp at inner / neumann-connected cells only;
      ! boundary cells must keep their prescribed Tp (set by
      ! update_boundaries) so the wind iteration starts from the
      ! correct boundary forcing.
      !
      do k = 1, no_nodes
         !
         if (inner(k)) Tp(k) = Tpini
         if (neumannconnected(k) > 0) Tp(neumannconnected(k)) = Tpini
         !
      enddo
      !
      do k = 1, no_nodes
         !
         sig(k) = 2 * pi / Tp(k)
         sig(k) = min(max(sig(k), sigmin), sigmax)
         !
         if (.not. inner(k)) then
            aa(:, k) = max(ee(:, k), waveps) / sig(k)
         endif
         !
         call compute_celerities(depth(k), sig(k), sinth, costh, ntheta, gamma, dhdx(k), dhdy(k), sinhkh(k), Hmx(k), kwav(k), cg(k), ctheta(:,k))
         !
      enddo
      !
   endif
   !   
   ! Sort coordinates in sweep directions (can we not do this already in snapwave_domain?)
   !
   shift = [0, 1, -1, 2]
   !
   do sweep = 1, 4
      !
      ra = x * cos(thetamean + 0.5 * pi * shift(sweep)) + y * sin(thetamean + 0.5 * pi * shift(sweep))
      call hpsort_eps_epw(no_nodes, ra , indx(:, sweep), 1.0e-6)
      !
   enddo
   !
   ! Set inner to false for all points at grid edge
   !
   do k = 1, no_nodes
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
         elseif (k1==1 .and. k2==1) then ! TL: for now still needed for a working IG solver
            !
			   inner(k) = .false.
            exit                       
            !
         endif
         !
      enddo
   enddo  
   !
   ! Start iteration
   !
   do iter = 1, niter * 4
      !
      sweep = mod(iter, 4) !TODO - TL: problem that we don't have option for sweep = 1 anymore?
      !
      if (sweep == 0) then
         sweep = 4
      endif
      !
      !write(*,*)'iter:', iter, 'sweep:', sweep
      !
      ! At start of each sweep, compute E, H, E_ig and H_ig, and aa
      !
      !$omp parallel do private(kn) schedule(static)
      do k = 1, no_nodes
         !
         if (ok(k) == 1) cycle
         !
         ee(:, k) = max(ee(:, k), waveps)
         !
         E(k)     = sum(ee(:, k)) * dtheta
         H(k)      = sqrt(8 * E(k) / rho / g)
         !
         if (igwaves) then
            !
            E_ig(k) = sum(ee_ig(:, k)) * dtheta
            H_ig(k) = sqrt(8 * E_ig(k) / rho / g)
            !
         endif
         !
         ! Set Neumann boundaries
         !
         if (neumannconnected(k) /= 0) then
            !
            ! Do we really need all of these? Hmx?
            !
            kn = neumannconnected(k) ! Index of internal point
            !
            sinhkh(k) = sinhkh(kn)
            kwav(k) = kwav(kn)
            Hmx(k) = Hmx(kn)
            ee(:, k) = ee(:, kn)
            ee_ig(:, k) = ee_ig(:, kn)
            ctheta(:, k) = ctheta(:, kn)
            cg(k) = cg(kn)
            !
            if (wind) then
               !
               sig(k) = sig(kn)
               Tp(k) = 2 * pi / sig(kn)
               WsorE(:, k) = WsorE(:, kn)
               WsorA(:, k) = WsorA(:, kn)
               aa(:, k) = aa(:, kn)
               !
            endif
            !
            Df(k) = Df(kn)
            Dw(k) = Dw(kn)
            !
         endif
         !
      enddo
      !$omp end parallel do
      !
      if (sweep == 1) then
         !
         eeold = ee
         !
         if (igwaves) then
             !
             eeold_ig = ee_ig
             !
         endif
         !
         !$omp parallel do schedule(static)
         do k = 1, no_nodes
            !
            Eold(k) = sum(eeold(:, k))
            !
         enddo
         !$omp end parallel do
         !
      endif
      !
      ! Update IG source and sink terms
      !
      if (igwaves) then
         !
         ! Do this in each first sweep, or (in case of iterative_srcig) in every sweep 
         !
         if (sweep == 1 .or. iterative_srcig) then
            !
            ! Determining of IG source term as defined in Leijnse et al. 2024 
            !      
            call determine_infragravity_source_sink_term(inner, no_nodes, ntheta, w, ds, prev, dtheta, cg_ig, nwav, depth, zb, H, &
                                                         ee, ee_ig, ig_opt, alphaigfac, &
                                                         alphaig_local, beta_local, srcig_local, gamma, gamma_fac_br)
            !
         endif
         !
      endif
      !
      !  Loop over all points depending on sweep direction
      !
      do count = 1, no_nodes
         !
         k = indx(count, sweep)
         !
         ! Skip non-inner (boundary) nodes
         !
         if (.not. inner(k)) cycle
         !
         ! For nodes below minimum depth: zero out energy and skip
         !
         if (depth(k) <= hmin) then
            !
            ee(:, k)    = 0.0
            if (wind) aa(:, k) = 0.0
            ee_ig(:, k) = 0.0
            cycle
            !
         endif
         !
         ! Skip nodes that have already converged
         !
         if (ok(k) == 1) cycle
         !
         ! Retrieve integrated quantities computed in the start-of-sweep pre-loop
         !
         Ek = E(k)
         Hk = H(k)
         !
         if (igwaves) then
            !
            Ek_ig = E_ig(k)
            Hk_ig = H_ig(k)
            !
         endif
         !
         ! --- Step 1: Upwind energy, group velocity (and action for wind) ------------
         !
         do itheta = 1, ntheta
            !
            k1 = prev(1, itheta, k)
            k2 = prev(2, itheta, k)
            !
            eeprev(itheta) = w(1, itheta, k) * ee(itheta, k1) + w(2, itheta, k) * ee(itheta, k2)
            cgprev(itheta) = w(1, itheta, k) * cg(k1)         + w(2, itheta, k) * cg(k2)
            !
            if (igwaves) then
               eeprev_ig(itheta) = w(1, itheta, k) * ee_ig(itheta, k1) + w(2, itheta, k) * ee_ig(itheta, k2)
               cgprev_ig(itheta) = w(1, itheta, k) * cg_ig(k1)         + w(2, itheta, k) * cg_ig(k2)
            endif
            !
            if (wind) then
               aaprev(itheta) = w(1, itheta, k) * aa(itheta, k1) + w(2, itheta, k) * aa(itheta, k2)
               aaprev(itheta) = min(aaprev(itheta), eeprev(itheta) / sigmin)
               aaprev(itheta) = max(aaprev(itheta), eeprev(itheta) / sigmax)
            endif
            !
         enddo
         !
         ! --- Step 2: Pre-solve sig and celerities from upwind Ek/Ak (wind only) ----
         !
         ! The upwind Ek/Ak provides the best pre-solve estimate of local sig, so that
         ! both the ee and aa matrices are assembled with a consistent cg(k).
         ! Post-solve, sig is updated from the solved Ek/Ak in Step 6.
         !
         if (wind) then
            !
            Ek     = sum(eeprev) * dtheta
            Ak     = sum(aaprev) * dtheta
            sig(k) = max(min(Ek / Ak, sigmax), sigmin)
            Ak     = Ek / sig(k)
            aaprev = min(aaprev, eeprev / sigmin)
            aaprev = max(aaprev, eeprev / sigmax)
            call compute_celerities(depth(k), sig(k), sinth, costh, ntheta, gamma, &
                                    dhdx(k), dhdy(k), sinhkh(k), Hmx(k), kwav(k), cg(k), ctheta(:, k))
            !
         endif
         !
         ! --- Step 3: Source and sink terms ------------------------------------------
         !
         ! Bottom friction
         !
         uorbi = 0.5 * sig(k) * Hk / sinhkh(k)
         Dfk   = 0.28 * rho * fw(k) * uorbi**3
         !
         ! Wave breaking (Baldock)
         ! First check if wave breaking could occur based on Baldock criterion (Hk > baldock_ratio * Hmx(k))
         !
         if (Hk > baldock_ratio * Hmx(k)) then
            !
            call baldock(rho, g, alfa, gamma, depth(k), Hk, 2 * pi / sig(k), baldock_exponent, Dwk, Hmx(k))
            !
         else
            !
            ! No wave breaking according to Baldock criterion
            !
            Dwk = 0.0
            !
         endif
         !
         ! Vegetation
         !
         if (vegetation) then
            !
            call vegatt(sig(k), no_nodes, kwav(k), no_secveg, veg_ah(k, :), veg_bstems(k, :), &
                        veg_Nstems(k, :), veg_Cd(k, :), depth(k), rho, g, Hk, Dvegk)
            !
         else
            !
            Dvegk = 0.0
            !
         endif
         !
         ! Energy dissipation ratio (with under-relaxation)
         !
         DoverE(k) = (1.0 - relax_factor_DoverE) * DoverE(k) &
                   + relax_factor_DoverE * (Dwk + Dfk + Dvegk) / max(Ek, 1.0e-6)
         !
         Df(k)   = Dfk
         Dw(k)   = Dwk
         Dveg(k) = Dvegk
         !
         ! Wind source terms and action dissipation ratio DoverA (wind only)
         !
         if (wind) then
            !
            if (iter == 1) then
               call windinput(u10(k), rho, g, depth(k), ntheta, windspreadfac(:, k), Ek, Ak, cg(k), &
                              eeprev, aaprev, ds(:, k), WsorE(:, k), WsorA(:, k), jadcgdx)
            else
               call windinput(u10(k), rho, g, depth(k), ntheta, windspreadfac(:, k), Ek, Ak, cg(k), &
                              ee(:, k), aa(:, k), ds(:, k), WsorE(:, k), WsorA(:, k), jadcgdx)
            endif
            !
            DwT  = -c_dispT / (1.0 - ndissip) * (2.0 * pi) / sig(k)**2 * cg(k) * kwav(k) * DoverE(k)
            DwAk = 1.0 / (2.0 * pi) * (E(k) * DwT + 2.0 * pi * Ak * DoverE(k))
            !
            if (iter == 1) then
               DoverA(k) = DwAk / max(Ak, 1.0e-6)
            else
               DoverA(k) = (1.0 - relax_factor_DoverA) * DoverA(k) &
                         + relax_factor_DoverA * DwAk / max(Ak, 1.0e-6)
            endif
            !
         endif
         !
         ! --- Step 4: Assemble and solve energy balance (ee) -------------------------
         !
         do itheta = 1, ntheta
            !
            R(itheta) = oneoverdt * ee(itheta, k) + cgprev(itheta) * eeprev(itheta) / ds(itheta, k) &
                      - srcig_local(itheta, k) * shinc2ig
            !
         enddo
         !
         do itheta = 2, ntheta - 1
            !
            A(itheta) = -ctheta(itheta - 1, k) * oneover2dtheta
            B(itheta) =  oneoverdt + cg(k) / ds(itheta, k) + DoverE(k)
            C(itheta) =  ctheta(itheta + 1, k) * oneover2dtheta
            !
         enddo
         !
         ! Directional edge BCs: first-order upwind, open at the grid ends, 
         ! same upwind scheme as the IG (ee_ig) and action (aa) balances
         !
         if (ctheta(1, k) < 0.0) then
            A(1) = 0.0
            B(1) = oneoverdt - ctheta(1, k) / dtheta + cg(k) / ds(1, k) + DoverE(k)
            C(1) = ctheta(2, k) / dtheta
         else
            A(1) = 0.0
            B(1) = oneoverdt + cg(k) / ds(1, k) + DoverE(k)
            C(1) = 0.0
         endif
         !
         if (ctheta(ntheta, k) > 0.0) then
            A(ntheta) = -ctheta(ntheta - 1, k) / dtheta
            B(ntheta) =  oneoverdt + ctheta(ntheta, k) / dtheta + cg(k) / ds(ntheta, k) + DoverE(k)
            C(ntheta) =  0.0
         else
            A(ntheta) = 0.0
            B(ntheta) = oneoverdt + cg(k) / ds(ntheta, k) + DoverE(k)
            C(ntheta) = 0.0
         endif
         !
         if (wind) R(:) = R(:) + WsorE(:, k)
         !
         call solve_tridiag(A, B, C, R, ee(:, k), ntheta)
         ee(:, k) = max(ee(:, k), waveps)
         !
         ! --- Step 5: Assemble and solve action balance (aa, wind only) ---------------
         !
         ! A and C are the same as for ee (refraction terms don't change).
         ! Only B_aa differs: DoverA instead of DoverE, plus upwind BC for ctheta endpoints.
         !
         if (wind) then
            !
            do itheta = 2, ntheta - 1
               B_aa(itheta) = oneoverdt + cg(k) / ds(itheta, k) + DoverA(k)
               R_aa(itheta) = oneoverdt * aa(itheta, k) + cgprev(itheta) * aaprev(itheta) / ds(itheta, k)
            enddo
            !
            if (ctheta(1, k) < 0.0) then
               B_aa(1) = oneoverdt - ctheta(1, k) / dtheta + cg(k) / ds(1, k) + DoverA(k)
            else
               B_aa(1) = oneoverdt + cg(k) / ds(1, k) + DoverA(k)
            endif
            !
            R_aa(1) = oneoverdt * aa(1, k) + cgprev(1) * aaprev(1) / ds(1, k)
            !
            if (ctheta(ntheta, k) > 0.0) then
               B_aa(ntheta) = oneoverdt + ctheta(ntheta, k) / dtheta + cg(k) / ds(ntheta, k) + DoverA(k)
            else
               B_aa(ntheta) = oneoverdt + cg(k) / ds(ntheta, k) + DoverA(k)
            endif
            !
            R_aa(ntheta) = oneoverdt * aa(ntheta, k) + cgprev(ntheta) * aaprev(ntheta) / ds(ntheta, k)
            R_aa(:) = R_aa(:) + WsorA(:, k)
            !
            call solve_tridiag(A, B_aa, C, R_aa, aa(:, k), ntheta)
            !
            aa(:, k) = max(aa(:, k), waveps / sigmax)
            aa(:, k) = max(aa(:, k), waveps / sig(k))
            !
         endif
         !
         ! --- Step 6: Depth-limit energy (and action), update sig and celerities -----
         !
         Ek          = sum(ee(:, k)) * dtheta
         depthlimfac = max(1.0, (sqrt(Ek / rhog8) / (gammax * depth(k)))**2)
         ee(:, k)    = ee(:, k) / depthlimfac
         !
!         if (wind) then
!            !
!            Ek       = Ek / depthlimfac
!            Ak       = sum(aa(:, k)) * dtheta
!            Ak       = Ak / depthlimfac
!            aa(:, k) = aa(:, k) / depthlimfac
!            sig(k)   = max(min(Ek / Ak, sigmax), sigmin)
!            !
!            call compute_celerities(depth(k), sig(k), sinth, costh, ntheta, gamma, &
!                                    dhdx(k), dhdy(k), sinhkh(k), Hmx(k), kwav(k), cg(k), ctheta(:, k))
!            !
!         endif
         !
         ! --- Step 7: IG wave balance (optional) -------------------------------------
         !
         if (igwaves) then
            !
            ! Update incident Hk from post-solve ee (needed for IG bottom friction)
            !
            Hk     = sqrt(8.0 * sum(ee(:, k)) * dtheta / rho / g)
            Dfk_ig = fw_ig(k) * 0.0361 * (9.81 / depth(k))**1.5 * Hk * Ek_ig
            !
            ! IG wave breaking (Baldock)
            !
            if (Hk_ig > baldock_ratio_ig * Hmx_ig(k)) then
               call baldock(rho, g, alfa_ig, gamma_ig, depth(k), Hk_ig, T_ig(k), baldock_exponent, Dwk_ig, Hmx_ig(k))
            else
               Dwk_ig = 0.0
            endif
            !
            Df_ig(k)     = Dfk_ig
            Dw_ig(k)     = Dwk_ig
            !
            ! Not using underrelaxation for IG dissipation for now, but we could add this if needed (relax_factor_DoverE_ig)
            DoverE_ig(k) = (Dwk_ig + Dfk_ig) / max(Ek_ig, 1.0e-6)
            ! 
            !
            ! IG RHS
            !
            do itheta = 1, ntheta
               R_ig(itheta) = oneoverdt * ee_ig(itheta, k) &
                            + cgprev_ig(itheta) * eeprev_ig(itheta) / ds(itheta, k) &
                            + srcig_local(itheta, k)
            enddo
            !
            ! IG matrix with directional boundary conditions
            !
            do itheta = 2, ntheta - 1
               A_ig(itheta) = -ctheta_ig(itheta - 1, k) * oneover2dtheta
               B_ig(itheta) =  oneoverdt + cg_ig(k) / ds(itheta, k) + DoverE_ig(k)
               C_ig(itheta) =  ctheta_ig(itheta + 1, k) * oneover2dtheta
            enddo
            !
            if (ctheta_ig(1, k) < 0.0) then
               A_ig(1) = 0.0
               B_ig(1) = oneoverdt - ctheta_ig(1, k) / dtheta + cg_ig(k) / ds(1, k) + DoverE_ig(k)
               C_ig(1) = ctheta_ig(2, k) / dtheta
            else
               A_ig(1) = 0.0
               B_ig(1) = oneoverdt + cg_ig(k) / ds(1, k) + DoverE_ig(k)
               C_ig(1) = 0.0
            endif
            !
            if (ctheta_ig(ntheta, k) > 0.0) then
               A_ig(ntheta) = -ctheta_ig(ntheta - 1, k) / dtheta
               B_ig(ntheta) =  oneoverdt + ctheta_ig(ntheta, k) / dtheta + cg_ig(k) / ds(ntheta, k) + DoverE_ig(k)
               C_ig(ntheta) =  0.0
            else
               A_ig(ntheta) = 0.0
               B_ig(ntheta) = oneoverdt + cg_ig(k) / ds(ntheta, k) + DoverE_ig(k)
               C_ig(ntheta) = 0.0
            endif
            !
            call solve_tridiag(A_ig, B_ig, C_ig, R_ig, ee_ig(:, k), ntheta)
            ee_ig(:, k) = max(ee_ig(:, k), 0.0)
            !
            ! Depth-limit IG energy
            !
            depthlimfac = max(1.0, (sqrt(sum(ee_ig(:, k)) * dtheta / rhog8) / (gammax * depth(k)))**2)
            ee_ig(:, k) = ee_ig(:, k) / depthlimfac
            !
         else
            !
            ee_ig(:, k) = 0.0
            !
         endif
         !
      enddo
      !
      if (sweep==4) then
         !
         ! Check convergence after all 4 sweeps
         !
         !$omp parallel do private(dee) schedule(static)
         do k = 1, no_nodes
            !
            dee     = ee(:, k) - eeold(:, k)
            diff(k) = maxval(abs(dee))
            !
            if (diff(k) / eemax < crit) then
               ok(k) = 1
            endif
            !
         enddo
         !$omp end parallel do
         !
         ! Percentage of converged points
         !
         percok = sum(ok) / dble(no_nodes) * 100.0
         eemax  = maxval(ee)
         !
         ! Relative maximum error
         !
         error = maxval(diff) / eemax
         !
         ! Check convergence of IG waves
         !
         if (igwaves) then
            !
            !$omp parallel do private(dee) schedule(static)
            do k = 1, no_nodes
               !
               dee        = ee_ig(:, k) - eeold_ig(:, k)
               diff_ig(k) = maxval(abs(dee))
               !
               if (diff_ig(k) / eemax_ig < crit) then
                  ok_ig(k) = 1
               endif
               !
            enddo
            !$omp end parallel do
            !
            percok_ig = sum(ok_ig) / dble(no_nodes) * 100.0
            eemax_ig  = maxval(ee_ig)
            error_ig  = maxval(diff_ig) / eemax_ig
            !
            write(logstr,'(a,i6,a,f10.5,a,f7.2,a,f10.5,a,f7.2)')'   iteration ', iter / 4 , &
               ' error = ', error,'   %ok = ', percok,'   error_ig = ', error_ig,'   %ok_ig = ', percok_ig
            call write_log(logstr, 0)
            !
            if ((error < crit .or. percok > 99.0)) then!.and. (error_ig < crit .or. percok_ig > 99.0)) then
               !
               write(logstr,'(a,i6,a,f10.5,a,f7.2,a,f10.5,a,f7.2)')'   converged at iteration ', iter / 4 , &
                  ' error = ', error,'   %ok = ', percok,'   error_ig = ', error_ig,'   %ok_ig = ', percok_ig
               call write_log(logstr, 0)
               exit
               !
            elseif (iter == niter * 4) then ! Made it to the end without reaching 'error<crit', still want output
               !
               write(logstr,'(a,i6,a,f10.5,a,f7.2,a,f10.5,a,f7.2)')'    ended at iteration ', iter / 4 , &
                  ' error = ', error,'   %ok = ', percok,'   error_ig = ', error_ig,'   %ok_ig = ', percok_ig
               call write_log(logstr, 0)
               !
            endif
            !
         else
            !
            write(logstr,'(a,i6,a,f10.5,a,f7.2)')'   iteration ', iter / 4 ,' error = ', error,'   %ok = ', percok
            call write_log(logstr, 0)
            !
            if (error < crit .or. percok > 99.0) then
               !
               write(logstr,'(a,i6,a,f10.5,a,f7.2)')'   converged at iteration ', iter / 4 ,' error = ',error,'   %ok = ', percok
               call write_log(logstr, 0)
               exit
               !
            elseif (iter == niter * 4) then ! Made it to the end without reaching 'error<crit', still want output
               !
               write(logstr,'(a,i6,a,f10.5,a,f7.2)')'    ended at iteration ', iter / 4 ,' error = ', error,'   %ok = ', percok
               call write_log(logstr, 0)
               !
            endif
            !
         endif
         !
      endif
      !
   enddo
   !
   !$omp parallel do schedule(static)
   do k = 1, no_nodes
      !
      ! Compute some directionally integrated parameters for output
      !
      if (depth(k) > hmin) then
         !
         E(k) = sum(ee(:, k)) * dtheta
         H(k) = sqrt(8 * E(k) / rho / g)
         !
         thetam(k) = atan2(sum(ee(:, k) * sin(theta)), sum(ee(:, k) * cos(theta)))
         !
         if (igwaves) then
            !
            ! IG wave height
            !
            E_ig(k) = sum(ee_ig(:, k)) * dtheta
            H_ig(k) = sqrt(8 * E_ig(k) / rho / g)
            !
            ! average beta, alphaig, srcig over directions
            !
            betamean(k) = sum(beta_local(:, k)) / ntheta ! real mean
            !
            alphaig(k) = sum(alphaig_local(:, k)) / ntheta ! real mean
            !
            srcig(k)   = sum(srcig_local(:, k)) / ntheta ! real mean
            !
         endif
         !
         if (wind) then
            !
            Tp(k)  = 2 * pi / sig(k)
            SwE(k) = sum(WsorE(:, k)) * dtheta
            SwA(k) = sum(WsorA(:, k)) * dtheta    !(2*pi*sum(WsorA(:,k))*dtheta - 2*pi/sig(k)*SwE(k)) / E(k)
            !
         endif
      endif
      !
   enddo
   !$omp end parallel do
   !
   ! For vegetation, need to compute unl using updated H(k)
   !
   if (vegetation) then
      ! 
      ! Compute the non-linear wave velocity time series (unl) using a wave shape model > only needs to be called once per calling SnapWave
      !          
      call swvegnonlin(no_nodes, kwav, depth, H, g, Tp, unl, etaw0)       
      !       
   endif   
   !
   ! Update wave forces - depends on vegetation
   !$omp parallel do schedule(static)
   do k = 1, no_nodes
      !
      ! Compute some directionally integrated parameters for output
      !
      if (depth(k) > hmin) then   
         !
         ! Update wave breaking dissipation > FIXME - should be called again using updated H(k) or not? 
         !if (wind) then
         !   call baldock(rho, g, alfa, gamma, depth(k), H(k), 2.0*pi/sig(k), 1, Dw(k), Hmx(k))
         !else
         !   call baldock(rho, g, alfa, gamma, depth(k), H(k), Tp(k), 1, Dw(k), Hmx(k))
         !endif
         !
         ! Update wave forces
         F(k) = Dw(k) * kwav(k) / sig(k) / rho / depth(k)
         !
         ! In case of vegetation, we also need to update the vegetation dissipation and compute wave drag force due to vegetation
         !
         if (vegetation) then
            !
            ! Compute wave dissipation due to vegetation > FIXME - correct that it should be called again to update Dveg based on updated H(k)?
            !
            call vegatt(sig(k), no_nodes, kwav(k), no_secveg, veg_ah(k,:), veg_bstems(k,:), veg_Nstems(k,:), veg_Cd(k,:), depth(k), rho, g, H(k), Dveg(k)) 
            !
            F(k) = F(k) + ( Dveg(k) * kwav(k) / sig(k) / rho / depth(k) )
            !
            ! Now also call 'momeqveg' to compute wave drag force due to vegetation Fvw
            !
            call momeqveg(no_nodes, no_secveg, veg_ah(k,:), veg_bstems(k,:), veg_Nstems(k,:), veg_Cd(k,:), depth(k), rho, H(k), Tp(k), unl(k,:), Fvw(k))
            ! NOTE - TL: for now replaced 'Trep' by 'Tp(k)'       
            !
            F(k) = F(k) + Fvw(k) 
            !
         endif
         !
      endif
      !
   enddo
   !$omp end parallel do
   !   
   callno = callno + 1
   !
   end subroutine solve_energy_balance2Dstat


   subroutine solve_tridiag(a, b, c, d, x, n)
   !
   ! Solve a tridiagonal system of equations.
   !
   ! Input:
   !   a : sub-diagonal
   !   b : main diagonal
   !   c : super-diagonal
   !   d : right-hand side
   !   n : number of equations
   !
   ! Output:
   !   x : solution
   !
   ! Called from: solve_energy_balance2Dstat (snapwave_solver)
   !
   implicit none
   !
   !	 a - sub-diagonal (means it is the diagonal below the main diagonal)
   !	 b - the main diagonal
   !	 c - sup-diagonal (means it is the diagonal above the main diagonal)
   !	 d - right part
   !	 x - the answer
   !	 n - number of equations
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
   cp(1) = c(1) / b(1)
   dp(1) = d(1) / b(1)
   !
   ! solve for vectors c-prime and d-prime
   !
   do i = 2, n
      m = b(i) - cp(i - 1) * a(i)
      cp(i) = c(i) / m
      dp(i) = (d(i) - dp(i - 1) * a(i)) / m
   enddo
   !
   ! initialize x
   !
   x(n) = dp(n)
   !
   ! solve for x from the vectors c-prime and d-prime
   !
   do i = n - 1, 1, - 1
      x(i) = dp(i) - cp(i) * x(i + 1)
   end do
   !
   end subroutine solve_tridiag

   subroutine baldock (rho, g, alfa, gamma, depth, H, T, iexp, Dw, Hmax)
   !
   ! Wave breaking dissipation according to Baldock et al. (1998).
   !
   ! Input:
   !   rho   : water density [kg/m3]
   !   g     : gravitational acceleration [m/s2]
   !   alfa  : coefficient alpha of the Baldock breaking model [-]
   !   gamma : breaker index gamma [-]
   !   depth : water depth [m]
   !   H     : wave height [m]
   !   T     : wave period [s]
   !   iexp  : exponent of the extra factor on the dissipation
   !   Hmax  : maximum wave height [m]
   !
   ! Output:
   !   Dw    : wave breaking dissipation [W/m2]
   !
   ! Called from: solve_energy_balance2Dstat (snapwave_solver)
   !
   real*4, intent(in)                :: rho
   real*4, intent(in)                :: g
   real*4, intent(in)                :: alfa
   real*4, intent(in)                :: gamma
   real*4, intent(in)                :: depth
   real*4, intent(in)                :: H
   real*4, intent(in)                :: T
   integer, intent(in)               :: iexp
   real*4, intent(out)               :: Dw
   real*4, intent(in)                :: Hmax
   real*4                            :: Hloc
   real*4                            :: f
   !
   ! Compute dissipation according to Baldock
   !
   Hloc = max(H, 1.e-6)
   !
   if (iexp > 0 .and. Hloc > Hmax) then
      !
      ! Add extra dissipation when Hloc exceeds Hmax.
      ! This is needed at very steep coast lines, where Baldock dissipation cannot always keep up with
      ! the wave height increase due to shoaling. The extra dissipation is added by multiplying
      ! the Baldock dissipation with a factor f, which is larger than 1 when Hloc > Hmax.
      !
      f = (Hloc / Hmax)**iexp
      !
   else  
      !
      f = 1.0
      !
   endif
   !
   Dw = 0.28 * alfa * rho * g / T * exp( - (Hmax / Hloc)**2) * (Hmax**2 + Hloc**2) * f
   !
   ! Other options for wave breaking dissipation (not used, but left here for reference)
   !
   ! Dw = 0.28 * alfa * rho * g / T * exp( - (Hmax / Hloc)**2) * (Hmax**3 + Hloc**3) / gamma / depth
   !
   end subroutine baldock
   
   
   subroutine timer(t)
   !
   ! Return the wall-clock time in seconds from system_clock.
   !
   ! Output:
   !   t : wall-clock time [s]
   !
   ! Called from: compute_wave_field (snapwave_solver)
   !
   real*4,intent(out)               :: t
   integer*8                        :: count,count_rate,count_max
   call system_clock (count,count_rate, count_max)
   t = real(dble(count) / dble(count_rate), 4)
   end subroutine timer

end module snapwave_solver 
