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
   !
   ! Strip schedule of one ordered sweep: nodes(ptr(s):ptr(s+1)-1) are the inner nodes of strip s,
   ! in sweep order. Strips lie side by side across the sweep direction and are swept in parallel.
   !
   type t_sweep_strips
      integer, dimension(:), allocatable :: ptr
      integer, dimension(:), allocatable :: nodes
   end type t_sweep_strips
   !
   ! Maximum number of directional bins for the per-node work buffers. These are fixed-size so
   ! they stay on the thread stack; the build uses heap arrays, which would otherwise malloc
   ! every automatic array on each node update and serialise the threads on the heap lock.
   !
   integer, parameter :: ntheta_max = 1024

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
      sigm_ig = 2.0 * pi / Tp_ig !TODO - TL: Question do we want Tp_ig now as constant, or also spatially varying like Tp ?
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
         !
         ! Neumann nodes take cg_ig of their inner node, as cg in the sweep, so that the IG
         ! energy copied to them is not amplified by cg_ig(Neumann) / cg_ig(inner)
         !
         do k = 1, no_nodes
            !
            if (neumannconnected(k) > 0) cg_ig(k) = cg_ig(neumannconnected(k))
            !
         enddo
         !
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
   use sfincs_error,  only: stop_sfincs
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
   real*4, dimension(:), allocatable                    :: Fvw                      ! vegetation wave drag force (vegetation only)
   real*4, dimension(:,:), allocatable                  :: unl                      ! non-linear wave orbital velocity time series, in 50 points per wave length (vegetation only)
   real*4, dimension(:,:), allocatable                  :: etaw0                    ! non-linear sea surface time series, in 50 points per wave length (vegetation only)
   !
   !
   ! Local variables and arrays
   !
   integer, dimension(:), allocatable         :: ok                     ! mask for fully iterated points
   integer, dimension(:), allocatable         :: ok_ig                  ! mask for fully iterated IG points   
   real*4                                     :: dtheta                 ! directional resolution
   real*4                                     :: eemax                  ! maximum wave energy density   
   real*4                                     :: eemax_ig               ! maximum IG wave energy density
   integer                                    :: sweep,iter            ! sweep number, number of iterations
   integer                                    :: k,k1,k2,count,kn,itheta ! counters (k is grid index)
   integer, dimension(:,:), allocatable       :: indx                   ! index for grid sorted per sweep direction
   real*4, dimension(:,:), allocatable        :: eeold               ! wave energy density, energy density previous iteration
   real*4, dimension(:,:), allocatable        :: eeold_ig            ! IG wave energy density, energy density previous iteration
   real*4, dimension(:), allocatable          :: Eold                   ! mean wave energy, previous iteration
   real*4, dimension(:,:), allocatable        :: srcig_local            ! Energy source/sink term because of IG wave energy transfer from incident waves
   real*4, dimension(:,:), allocatable        :: beta_local             ! Local bed slope based on bed level per direction      
   real*4, dimension(:,:), allocatable        :: alphaig_local          ! Local infragravity wave shoaling parameter alpha
   real*4, dimension(:), allocatable          :: DoverA                 ! ratio of mean wave dissipation over mean wave energy
   real*4, dimension(:), allocatable          :: DoverE_ig              ! ratio of mean wave dissipation over mean wave energy
   real*4, dimension(:), allocatable          :: E                      ! mean wave energy
   real*4, dimension(:), allocatable          :: E_ig                   ! mean wave energy
   real*4, dimension(:), allocatable          :: diff                   ! maximum difference of wave energy relative to previous iteration
   real*4, dimension(:), allocatable          :: diff_ig                ! maximum difference of IG wave energy relative to previous iteration
   real*4, dimension(:,:), allocatable        :: ra                     ! coordinate in sweep direction, per sweep
   real*4, dimension(:), allocatable          :: sigm_ig
   integer, dimension(4)                      :: shift
   real*4                                     :: pi = 4.0 * atan(1.0)
   real*4                                     :: g = 9.81
   real*4                                     :: hmin                   ! minimum water depth! TL: make user changeable also here according to 'snapwave_hmin' in sfincs.inp   
   real*4                                     :: oneoverdt
   real*4                                     :: oneover2dtheta
   real*4                                     :: rhog8
   real*4                                     :: dtot                   ! sum of the dissipation terms at a node before scaling
   real*4                                     :: percok
   real*4                                     :: percok_ig              ! percentage of converged IG points   
   real*4                                     :: error
   real*4                                     :: error_ig               ! relative maximum IG wave error
   real*4                                     :: eemax_new              ! maximum wave energy density in the current convergence check
   real*4                                     :: diffmax                ! maximum difference with the previous iteration
   integer                                    :: nok                    ! number of converged points
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
   real*4                                     :: ndissip    
   real*4                                     :: waveps=0.0001
   !
   ! Strip-parallel sweeps
   !
   type(t_sweep_strips), dimension(4)         :: strips                 ! strip schedule per sweep
   integer, dimension(:,:), allocatable       :: strip_of               ! strip of each node per axis (0 = not inner)
   integer, dimension(:,:), allocatable       :: halo_nodes             ! inner nodes read by another strip, per axis
   integer, dimension(2)                      :: halo_n                 ! number of halo nodes per axis
   integer                                    :: nstrip                 ! number of strips used in this call
   integer                                    :: istrip, ii, axis
   real*4, dimension(:,:), allocatable        :: ee_s, ee_ig_s, aa_s    ! start-of-sweep copies of halo nodes
   real*4, dimension(:), allocatable          :: cg_s                   ! start-of-sweep copy of cg of halo nodes
   !
   ! Allocate local arrays
   !
   waveps = 0.0001
   !
   allocate(ok(no_nodes)); ok=0
   allocate(indx(no_nodes,4)); indx=0
   allocate(eeold(ntheta,no_nodes)) ! set in the first sweep
   allocate(E(no_nodes)); E=waveps
   allocate(Eold(no_nodes)); Eold=0.0   
   !
   if (igwaves) then
      allocate(DoverE_ig(no_nodes)); DoverE_ig=0.0
      allocate(E_ig(no_nodes)); E_ig=waveps
      allocate(sigm_ig(no_nodes)); sigm_ig=0.0
      allocate(beta_local(ntheta,no_nodes)); beta_local=0.0
      allocate(alphaig_local(ntheta,no_nodes)); alphaig_local=0.0
      allocate(eeold_ig(ntheta,no_nodes)) ! set in the first sweep
      allocate(diff_ig(no_nodes)); diff_ig=0.0
      allocate(ok_ig(no_nodes)); ok_ig=0
   endif
   !
   if (wind) then
      allocate(DoverA(no_nodes)); DoverA=0.0
   endif
   !
   allocate(diff(no_nodes)); diff=0.0
   allocate(ra(no_nodes, 4)); ra=0.0
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
   !
   ! Vegetation work arrays (no_nodes x 50) are only needed with vegetation
   !
   if (vegetation) then
      allocate(Fvw(no_nodes)); Fvw = 0.0
      allocate(unl(no_nodes, 50)); unl = 0.0
      allocate(etaw0(no_nodes, 50)); etaw0 = 0.0
   endif
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
   ! The four sorts are independent, so they run concurrently (each sort itself is serial,
   ! so the order is the same as when they run one after the other)
   !
   !$omp parallel do private(sweep) schedule(static, 1)
   do sweep = 1, 4
      !
      ra(:, sweep) = x * cos(thetamean + 0.5 * pi * shift(sweep)) + y * sin(thetamean + 0.5 * pi * shift(sweep))
      call hpsort_eps_epw(no_nodes, ra(:, sweep), indx(:, sweep), 1.0e-6)
      !
   enddo
   !$omp end parallel do
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
   ! DIAGNOSTIC: strip-decomposition statistics of the ordered sweeps (first call only)
   !
   if (ntheta > ntheta_max) then
      write(logstr, '(a,i0,a,i0)') 'Error ! SnapWave: number of directional bins ', ntheta, ' exceeds maximum ', ntheta_max
      call stop_sfincs(logstr, 1)
   endif
   !
   if (callno == 1) call compute_sweep_strip_stats(x, y, no_nodes, ntheta, prev, w, inner, indx, depth, hmin, thetamean)
   !
   ! Strip schedule of the four sweeps (indx and inner are fixed during this call)
   !
   allocate(strip_of(no_nodes, 2), halo_nodes(no_nodes, 2))
   call build_sweep_strips(no_nodes, ntheta, prev, inner, depth, hmin, indx, &
                           nstrip, strips, strip_of, halo_n, halo_nodes)
   !
   if (callno == 1 .and. nstrip > 1) then
      write(logstr, '(a,i0,a,i0,a,i0,a)') 'SnapWave: parallel sweep with ', nstrip, ' strips (', &
         halo_n(1), ' + ', halo_n(2), ' halo nodes)'
      call write_log(logstr, 0)
   endif
   !
   if (nstrip > 1) then
      allocate(ee_s(ntheta, no_nodes)); ee_s = 0.0
      allocate(cg_s(no_nodes)); cg_s = 0.0
      if (igwaves) then
         allocate(ee_ig_s(ntheta, no_nodes)); ee_ig_s = 0.0
      endif
      if (wind) then
         allocate(aa_s(ntheta, no_nodes)); aa_s = 0.0
      endif
   endif
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
            ! Neumann nodes are not inner nodes and are skipped in the sweep, so DoverE and
            ! Dveg are copied as well. The wave force scales Dw with DoverE after the solve.
            !
            Df(k)     = Df(kn)
            Dw(k)     = Dw(kn)
            Dveg(k)   = Dveg(kn)
            DoverE(k) = DoverE(kn)
            !
         endif
         !
      enddo
      !$omp end parallel do
      !
      if (sweep == 1) then
         !
         !$omp parallel do schedule(static)
         do k = 1, no_nodes
            !
            eeold(:, k) = ee(:, k)
            if (igwaves) eeold_ig(:, k) = ee_ig(:, k)
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
      ! Strip-parallel Gauss-Seidel sweep. The domain is split into strips lying side by side
      ! across the sweep direction; each strip is swept serially in sweep order by one thread.
      ! Upwind values from the same strip are read live (as in the serial sweep); upwind values
      ! from another strip are read from a start-of-sweep copy (lagged by one sweep), so strips
      ! are independent. With nstrip = 1 this is exactly the serial sweep.
      !
      axis = 1
      if (sweep == 2 .or. sweep == 3) axis = 2
      !
      if (nstrip > 1) then
         !
         !$omp parallel do private(ii, k) schedule(static)
         do ii = 1, halo_n(axis)
            k = halo_nodes(ii, axis)
            ee_s(:, k) = ee(:, k)
            cg_s(k)    = cg(k)
            if (igwaves) ee_ig_s(:, k) = ee_ig(:, k)
            if (wind)    aa_s(:, k)    = aa(:, k)
         enddo
         !$omp end parallel do
         !
      endif
      !
      !$omp parallel do if (nstrip > 1) schedule(dynamic, 1) private(istrip, ii)
      do istrip = 1, nstrip
         do ii = strips(sweep)%ptr(istrip), strips(sweep)%ptr(istrip + 1) - 1
            call update_node(strips(sweep)%nodes(ii), strip_of(strips(sweep)%nodes(ii), axis))
         enddo
      enddo
      !$omp end parallel do
      !
      if (sweep==4) then
         !
         ! Check convergence after all 4 sweeps
         !
         ! The maxima and the count are reductions over all nodes (exact in any order). The
         ! test uses eemax of the previous check, the new eemax is used from the next check on.
         !
         eemax_new = -huge(1.0)
         diffmax   = -huge(1.0)
         nok       = 0
         !
         !$omp parallel do private(k, itheta) reduction(max:eemax_new, diffmax) reduction(+:nok) schedule(static)
         do k = 1, no_nodes
            !
            diff(k) = 0.0
            do itheta = 1, ntheta
               diff(k)   = max(diff(k), abs(ee(itheta, k) - eeold(itheta, k)))
               eemax_new = max(eemax_new, ee(itheta, k))
            enddo
            diffmax = max(diffmax, diff(k))
            !
            if (diff(k) / eemax < crit) then
               ok(k) = 1
            endif
            !
            nok = nok + ok(k)
            !
         enddo
         !$omp end parallel do
         !
         ! Percentage of converged points
         !
         percok = nok / dble(no_nodes) * 100.0
         eemax  = eemax_new
         !
         ! Relative maximum error
         !
         error = diffmax / eemax
         !
         ! Check convergence of IG waves
         !
         if (igwaves) then
            !
            eemax_new = -huge(1.0)
            diffmax   = -huge(1.0)
            nok       = 0
            !
            !$omp parallel do private(k, itheta) reduction(max:eemax_new, diffmax) reduction(+:nok) schedule(static)
            do k = 1, no_nodes
               !
               diff_ig(k) = 0.0
               do itheta = 1, ntheta
                  diff_ig(k) = max(diff_ig(k), abs(ee_ig(itheta, k) - eeold_ig(itheta, k)))
                  eemax_new  = max(eemax_new, ee_ig(itheta, k))
               enddo
               diffmax = max(diffmax, diff_ig(k))
               !
               if (diff_ig(k) / eemax_ig < crit) then
                  ok_ig(k) = 1
               endif
               !
               nok = nok + ok_ig(k)
               !
            enddo
            !$omp end parallel do
            !
            percok_ig = nok / dble(no_nodes) * 100.0
            eemax_ig  = eemax_new
            error_ig  = diffmax / eemax_ig
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
   !$omp parallel do private(dtot) schedule(static)
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
         ! Make the breaking dissipation consistent with the solved energy balance.
         !
         ! Dw(k) was stored in the sweep from the wave height before the solve, whereas the
         ! solve removes DoverE(k) * E(k) from the balance (the dissipation is implicit in the
         ! energy). Where the solution has not converged these differ. At nodes with a depth
         ! just above hmin the pre-solve wave height can be many times Hmx, while the solved
         ! wave height is a few cm. With the factor (H / Hmx)**baldock_exponent the pre-solve
         ! breaking dissipation is then orders of magnitude larger than the energy flux that
         ! arrives at the node, and the wave force Dw * k / sig below followed that value.
         !
         ! Dw(k) is therefore capped at its share of what the solve removed, which cannot
         ! exceed the incoming energy flux. The share is taken from the pre-solve terms. The
         ! cap only reduces Dw(k): where the solve removed more than the pre-solve terms
         ! (e.g. because DoverE is under-relaxed), Dw(k) is left as it is. At converged
         ! nodes nothing changes.
         !
         ! Only Dw(k) is scaled:
         ! - Df(k) is not used after the solve (output only) and is left as it is, so at a
         !   node that has not converged it can still belong to the pre-solve wave height.
         ! - Dveg(k) is recomputed below from the solved wave height when vegetation is on.
         !
         ! Possible improvements:
         ! - Store the solved dissipation per term in the sweep (DoverE split into its
         !   breaking, friction and vegetation parts, times the solved energy), so that all
         !   three are consistent with the output wave height and nothing needs scaling here.
         ! - Compute the wave force from the gradients of the radiation stresses instead of
         !   from the dissipation (F = -dSxx/dx - dSxy/dy, ...). The force then follows
         !   from the solved energy field alone, includes the part that is not related to
         !   dissipation (shoaling, refraction), and cannot exceed what the wave field holds.
         !   That is where this should go in the end.
         !
         dtot = Dw(k) + Df(k) + Dveg(k)
         !
         if (dtot > 0.0) then
            !
            Dw(k) = Dw(k) * min(DoverE(k) * E(k) / dtot, 1.0)
            !
         endif
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
   contains
   !
   subroutine update_node(k, sk)
   !
   ! Node update with work buffers on the calling thread's stack (see ntheta_max).
   !
   ! Input:
   !   k  : node index
   !   sk : strip of node k for the current sweep axis
   !
   ! Called from: solve_energy_balance2Dstat (snapwave_solver), strip loop
   !
   integer, intent(in)            :: k, sk
   real*4, dimension(ntheta_max, 15) :: wrk
   !
   call update_node_body(k, sk, wrk(:, 1), wrk(:, 2), wrk(:, 3), wrk(:, 4), wrk(:, 5), &
                         wrk(:, 6), wrk(:, 7), wrk(:, 8), wrk(:, 9), wrk(:, 10), &
                         wrk(:, 11), wrk(:, 12), wrk(:, 13), wrk(:, 14), wrk(:, 15))
   !
   end subroutine update_node
   !
   subroutine update_node_body(k, sk, eeprev, cgprev, eeprev_ig, cgprev_ig, aaprev, &
                               A, B, C, R, A_ig, B_ig, C_ig, R_ig, B_aa, R_aa)
   !
   ! Update wave energy (ee), and optionally IG energy (ee_ig) and wave action (aa), of node k in
   ! the current sweep. Upwind values from strip sk (or from non-inner nodes) are read live,
   ! upwind values from other strips from the start-of-sweep copies (*_s).
   !
   ! All per-node scalars and work arrays are local, so they live on the calling thread's stack.
   ! Everything else is host-associated from solve_energy_balance2Dstat; only entries of node k
   ! are written.
   !
   ! Input:
   !   k  : node index
   !   sk : strip of node k for the current sweep axis
   !
   ! Work arrays (explicit shape, storage provided by update_node):
   !   eeprev .. R_aa : per-direction scratch of length ntheta
   !
   ! Called from: update_node (snapwave_solver)
   !
   integer, intent(in)       :: k, sk
   !
   integer                   :: k1, k2, l1, l2, itheta
   real*4                    :: ee1, ee2, cg1, cg2, eig1, eig2, aa1, aa2
   real*4                    :: Ek, Hk, Ek_ig, Hk_ig, Ak, uorbi, Dfk, Dwk, Dvegk, DwT, DwAk, depthlimfac, Dfk_ig, Dwk_ig
   real*4, dimension(ntheta) :: eeprev, cgprev, eeprev_ig, cgprev_ig, aaprev ! explicit-shape dummies
   real*4, dimension(ntheta) :: A, B, C, R, A_ig, B_ig, C_ig, R_ig, B_aa, R_aa
   !
      ! Skip non-inner (boundary) nodes
      !
      if (.not. inner(k)) return
      !
      ! For nodes below minimum depth: zero out energy and skip
      !
      if (depth(k) <= hmin) then
         !
         ee(:, k)    = 0.0
         if (wind) aa(:, k) = 0.0
         ee_ig(:, k) = 0.0
         return
         !
      endif
      !
      ! Skip nodes that have already converged
      !
      if (ok(k) == 1) return
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
         ! Upwind node in the same strip (or not inner, i.e. never updated in the sweep): live value.
         ! Upwind node in another strip: start-of-sweep copy.
         !
         l1 = strip_of(k1, axis)
         l2 = strip_of(k2, axis)
         !
         if (l1 == sk .or. l1 == 0) then
            ee1 = ee(itheta, k1)
            cg1 = cg(k1)
         else
            ee1 = ee_s(itheta, k1)
            cg1 = cg_s(k1)
         endif
         !
         if (l2 == sk .or. l2 == 0) then
            ee2 = ee(itheta, k2)
            cg2 = cg(k2)
         else
            ee2 = ee_s(itheta, k2)
            cg2 = cg_s(k2)
         endif
         !
         eeprev(itheta) = w(1, itheta, k) * ee1 + w(2, itheta, k) * ee2
         cgprev(itheta) = w(1, itheta, k) * cg1 + w(2, itheta, k) * cg2
         !
         if (igwaves) then
            !
            if (l1 == sk .or. l1 == 0) then
               eig1 = ee_ig(itheta, k1)
            else
               eig1 = ee_ig_s(itheta, k1)
            endif
            !
            if (l2 == sk .or. l2 == 0) then
               eig2 = ee_ig(itheta, k2)
            else
               eig2 = ee_ig_s(itheta, k2)
            endif
            !
            eeprev_ig(itheta) = w(1, itheta, k) * eig1 + w(2, itheta, k) * eig2
            cgprev_ig(itheta) = w(1, itheta, k) * cg_ig(k1) + w(2, itheta, k) * cg_ig(k2)
            !
         endif
         !
         if (wind) then
            !
            if (l1 == sk .or. l1 == 0) then
               aa1 = aa(itheta, k1)
            else
               aa1 = aa_s(itheta, k1)
            endif
            !
            if (l2 == sk .or. l2 == 0) then
               aa2 = aa(itheta, k2)
            else
               aa2 = aa_s(itheta, k2)
            endif
            !
            aaprev(itheta) = w(1, itheta, k) * aa1 + w(2, itheta, k) * aa2
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
   end subroutine update_node_body
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
   real*4,dimension(ntheta_max) :: cp,dp ! fixed size: stays on the stack (n <= ntheta_max)
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
   !
   !
   subroutine compute_sweep_strip_stats(x, y, no_nodes, ntheta, prev, w, inner, indx, depth, hmin, thetamean)
   !
   ! DIAGNOSTIC (snapwave_stripstats = 1): cost of a strip decomposition of the ordered sweeps.
   !
   ! For a parallel sweep, each sweep's domain would be split into nstrip strips lying side by
   ! side across the sweep direction (equal numbers of active nodes per strip). Each strip is
   ! swept serially in the existing sweep order; upwind values from another strip are read from
   ! a start-of-sweep copy (lagged by one sweep).
   !
   ! For nstrip = 8, 24, 48 and each sweep, this logs:
   !   - the share of "forward" upwind reads (upwind node earlier in the sweep order, i.e. read
   !     as an updated value in the serial sweep) that would cross a strip boundary and become
   !     lagged, both as a count and weighted with the interpolation weights w;
   !   - the share of active nodes with at least one such lagged read;
   !   - for reference, the share of "backward" reads, which are lagged in the serial sweep already.
   !
   ! Active nodes: inner and depth > hmin. Reads from non-inner nodes are never lagged (those
   ! nodes are not updated during the sweep).
   !
   ! Called from: solve_energy_balance2Dstat (snapwave_solver), first call only
   !
   use snapwave_data,  only: stripstats
   use hpsort_eps,     only: hpsort_eps_epw
   !
   implicit none
   !
   integer,                                 intent(in) :: no_nodes, ntheta
   real*8,  dimension(no_nodes),            intent(in) :: x, y
   integer, dimension(2, ntheta, no_nodes), intent(in) :: prev
   real*4,  dimension(2, ntheta, no_nodes), intent(in) :: w
   logical, dimension(no_nodes),            intent(in) :: inner
   integer, dimension(no_nodes, 4),         intent(in) :: indx
   real*4,  dimension(no_nodes),            intent(in) :: depth
   real*4,                                  intent(in) :: hmin
   real*4,                                  intent(in) :: thetamean
   !
   integer, dimension(3), parameter   :: nstrips = [8, 24, 48]
   integer, dimension(4), parameter   :: shift   = [0, 1, -1, 2]
   real*4,                parameter   :: pi      = 3.14159265
   !
   integer, dimension(:), allocatable :: pos, strip, ind
   real*4,  dimension(:), allocatable :: rp
   logical, dimension(:), allocatable :: active
   integer                            :: is, nstrip, sweep, count, k, j, i, itheta, nact, iact, nlagnode
   integer*8                          :: nfwd, nlag, nback
   real*8                             :: wfwd, wlag
   logical                            :: lagged
   real*4                             :: phi
   !
   if (stripstats /= 1) return
   !
   allocate(pos(no_nodes), strip(no_nodes), ind(no_nodes), rp(no_nodes), active(no_nodes))
   !
   nact = 0
   do k = 1, no_nodes
      active(k) = inner(k) .and. depth(k) > hmin
      if (active(k)) nact = nact + 1
   enddo
   !
   call write_log('', 0)
   write(logstr, '(a,i0,a,i0,a,i0)') 'SnapWave strip stats: nodes = ', no_nodes, ', active = ', nact, ', ntheta = ', ntheta
   call write_log(logstr, 0)
   call write_log(' nstrip  sweep  nodes/strip  %fwd_reads_lagged  %fwd_weight_lagged  %nodes_with_lag  %backward_reads(serial)', 0)
   !
   do is = 1, size(nstrips)
      !
      nstrip = nstrips(is)
      !
      do sweep = 1, 4
         !
         do count = 1, no_nodes
            pos(indx(count, sweep)) = count
         enddo
         !
         ! Strips: sort on the coordinate perpendicular to the sweep direction, cut into equal shares of active nodes
         !
         phi = thetamean + 0.5 * pi * shift(sweep)
         do k = 1, no_nodes
            rp(k) = real(-x(k) * sin(phi) + y(k) * cos(phi), 4)
         enddo
         ind(1) = 0
         call hpsort_eps_epw(no_nodes, rp, ind, 1.0e-6)
         !
         strip = 0
         iact  = 0
         do count = 1, no_nodes
            k = ind(count)
            if (.not. inner(k)) cycle
            strip(k) = 1 + min(nstrip - 1, int(dble(iact) * nstrip / dble(max(nact, 1))))
            if (active(k)) iact = iact + 1
         enddo
         !
         nfwd     = 0
         nlag     = 0
         nback    = 0
         wfwd     = 0.0d0
         wlag     = 0.0d0
         nlagnode = 0
         !
         do k = 1, no_nodes
            !
            if (.not. active(k)) cycle
            !
            lagged = .false.
            !
            do itheta = 1, ntheta
               do i = 1, 2
                  !
                  j = prev(i, itheta, k)
                  if (j <= 0 .or. j == k) cycle
                  if (.not. inner(j)) cycle
                  if (w(i, itheta, k) <= 0.0) cycle
                  !
                  if (pos(j) < pos(k)) then
                     nfwd = nfwd + 1
                     wfwd = wfwd + w(i, itheta, k)
                     if (strip(j) /= strip(k)) then
                        nlag   = nlag + 1
                        wlag   = wlag + w(i, itheta, k)
                        lagged = .true.
                     endif
                  else
                     nback = nback + 1
                  endif
                  !
               enddo
            enddo
            !
            if (lagged) nlagnode = nlagnode + 1
            !
         enddo
         !
         write(logstr, '(i7,i7,f13.1,f19.2,f20.2,f17.2,f25.2)') nstrip, sweep, dble(nact) / nstrip, &
            100.0d0 * dble(nlag) / dble(max(nfwd, 1_8)), 100.0d0 * wlag / max(wfwd, 1.0d-12), &
            100.0d0 * dble(nlagnode) / dble(max(nact, 1)), 100.0d0 * dble(nback) / dble(max(nfwd + nback, 1_8))
         call write_log(logstr, 0)
         !
      enddo
      !
   enddo
   !
   call write_log('', 0)
   !
   deallocate(pos, strip, ind, rp, active)
   !
   end subroutine compute_sweep_strip_stats
   !
   !
   subroutine build_sweep_strips(no_nodes, ntheta, prev, inner, depth, hmin, indx, &
                                 nstrip, strips, strip_of, halo_n, halo_nodes)
   !
   ! Build the strip schedule of the four ordered sweeps.
   !
   ! Number of strips (snapwave_nstrips): 1 = serial sweep, >1 = fixed, 0 = automatic:
   ! min(number of OpenMP threads, active nodes * ntheta / snapwave_strip_minwork), at least 1.
   ! The work per node scales with ntheta, so models with many directions get more strips.
   !
   ! For each of the two sweep axes (axis 1: sweeps 1 and 4 along thetamean; axis 2: sweeps 2
   ! and 3, perpendicular to it), inner nodes are ordered on the coordinate perpendicular to the
   ! sweep direction and cut into nstrip strips with equal numbers of active nodes. That order
   ! is the sweep order of the perpendicular sweep: sweep 2 (thetamean + pi/2) for axis 1, and
   ! sweep 4 (thetamean + pi) for axis 2, so no extra sort is needed. Per sweep,
   ! the inner nodes of each strip are listed in sweep order (indx). Halo nodes are inner nodes
   ! that are an upwind node of an active node in another strip; only these need a
   ! start-of-sweep copy.
   !
   ! Input:
   !   no_nodes, ntheta : number of nodes and directions
   !   prev             : two upwind nodes per direction and node
   !   inner            : inner-node mask
   !   depth, hmin      : water depth and minimum depth
   !   indx             : node order per sweep
   !
   ! Output:
   !   nstrip           : number of strips used
   !   strips           : strip schedule per sweep
   !   strip_of         : strip of each node per axis (0 = not inner)
   !   halo_n           : number of halo nodes per axis
   !   halo_nodes       : halo nodes per axis
   !
   ! Called from: solve_energy_balance2Dstat (snapwave_solver)
   !
   use snapwave_data, only: nstrips, strip_minwork
   !$ use omp_lib,    only: omp_get_max_threads
   !
   implicit none
   !
   integer,                                 intent(in)    :: no_nodes, ntheta
   integer, dimension(2, ntheta, no_nodes), intent(in)    :: prev
   logical, dimension(no_nodes),            intent(in)    :: inner
   real*4,  dimension(no_nodes),            intent(in)    :: depth
   real*4,                                  intent(in)    :: hmin
   integer, dimension(no_nodes, 4),         intent(in)    :: indx
   integer,                                 intent(out)   :: nstrip
   type(t_sweep_strips), dimension(4),      intent(inout) :: strips
   integer, dimension(no_nodes, 2),         intent(out)   :: strip_of
   integer, dimension(2),                   intent(out)   :: halo_n
   integer, dimension(no_nodes, 2),         intent(out)   :: halo_nodes
   !
   integer, dimension(2), parameter   :: axis_order = [2, 4]          ! sweep whose order runs across axis 1 and axis 2
   integer, dimension(4), parameter   :: sweep_axis = [1, 2, 2, 1]
   !
   integer, dimension(:), allocatable :: next, mark
   logical, dimension(:), allocatable :: active
   integer                            :: nthreads, nact, nin, iact, a, sweep, count, k, j, i, itheta, s
   !
   allocate(active(no_nodes))
   !
   nact = 0
   nin  = 0
   do k = 1, no_nodes
      active(k) = inner(k) .and. depth(k) > hmin
      if (active(k)) nact = nact + 1
      if (inner(k))  nin  = nin + 1
   enddo
   !
   ! Number of strips
   !
   nthreads = 1
   !$ nthreads = omp_get_max_threads()
   !
   if (nstrips == 0) then
      nstrip = max(1, min(nthreads, int(dble(nact) * ntheta / dble(max(strip_minwork, 1)))))
   else
      nstrip = max(1, min(nstrips, max(nact, 1)))
   endif
   !
   ! Strip of each node per axis
   !
   strip_of = 0
   halo_n   = 0
   !
   if (nstrip == 1) then
      !
      do k = 1, no_nodes
         if (inner(k)) strip_of(k, :) = 1
      enddo
      !
   else
      !
      allocate(mark(no_nodes))
      !
      do a = 1, 2
         !
         iact = 0
         do count = 1, no_nodes
            k = indx(count, axis_order(a))
            if (.not. inner(k)) cycle
            strip_of(k, a) = 1 + min(nstrip - 1, int(dble(iact) * nstrip / dble(max(nact, 1))))
            if (active(k)) iact = iact + 1
         enddo
         !
         ! Halo nodes: inner upwind nodes of an active node in another strip
         !
         mark = 0
         !
         !$omp parallel do private(k, itheta, i, j) schedule(static)
         do k = 1, no_nodes
            if (.not. active(k)) cycle
            do itheta = 1, ntheta
               do i = 1, 2
                  j = prev(i, itheta, k)
                  if (j <= 0) cycle
                  if (strip_of(j, a) /= 0 .and. strip_of(j, a) /= strip_of(k, a)) then
                     !$omp atomic write
                     mark(j) = 1
                  endif
               enddo
            enddo
         enddo
         !$omp end parallel do
         !
         do k = 1, no_nodes
            if (mark(k) == 1) then
               halo_n(a) = halo_n(a) + 1
               halo_nodes(halo_n(a), a) = k
            endif
         enddo
         !
      enddo
      !
      deallocate(mark)
      !
   endif
   !
   ! Per sweep: inner nodes per strip, in sweep order (counting sort on strip)
   !
   allocate(next(nstrip + 1))
   !
   do sweep = 1, 4
      !
      a = sweep_axis(sweep)
      !
      if (allocated(strips(sweep)%ptr))   deallocate(strips(sweep)%ptr)
      if (allocated(strips(sweep)%nodes)) deallocate(strips(sweep)%nodes)
      allocate(strips(sweep)%ptr(nstrip + 1))
      allocate(strips(sweep)%nodes(max(nin, 1)))
      !
      strips(sweep)%ptr = 0
      do k = 1, no_nodes
         s = strip_of(k, a)
         if (s > 0) strips(sweep)%ptr(s + 1) = strips(sweep)%ptr(s + 1) + 1
      enddo
      strips(sweep)%ptr(1) = 1
      do s = 1, nstrip
         strips(sweep)%ptr(s + 1) = strips(sweep)%ptr(s + 1) + strips(sweep)%ptr(s)
      enddo
      !
      next = strips(sweep)%ptr
      do count = 1, no_nodes
         k = indx(count, sweep)
         s = strip_of(k, a)
         if (s > 0) then
            strips(sweep)%nodes(next(s)) = k
            next(s) = next(s) + 1
         endif
      enddo
      !
   enddo
   !
   deallocate(next, active)
   !
   end subroutine build_sweep_strips

end module snapwave_solver 
