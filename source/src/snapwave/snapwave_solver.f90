module snapwave_solver

   contains
   
   subroutine compute_wave_field()
      !
      use snapwave_data
      use snapwave_results
      !
      implicit none
      !
      real*4    :: hrmsb
      real*4    :: tpb
      real*4    :: wdirb
      real*4    :: msb
      !
      real*4    :: thetamin
      real*4    :: thetamax
      real*4    :: E0
      real*4    :: sigm
      real*4    :: sigm_ig
      real*4    :: Tpb_ig
      !
      integer   :: itheta
      integer   :: k
      integer   :: i1, i2, ii, jj
      !
      Tpb     = tpmean_bwv
      sigm    = 2*pi/Tpb
      Tpb_ig  = 7*Tpb
      sigm_ig = 2*pi/Tpb_ig
      !
!      call timer(t0)
      !
      ! Compute celerities and refraction speed
      !
!      write(*,*)'max depth=',maxval(depth),' min depth = ',minval(depth)
      call disper_approx(depth, Tpb,    kwav,    nwav,    C,    Cg,    no_nodes)
      call disper_approx(depth, Tpb_ig, kwav_ig, nwav_ig, C_ig, Cg_ig, no_nodes)
      cg_ig = cg   !TL: ??? is Cg = cg? And why put cg_ig equal to normal cg?
      !
      do k = 1, no_nodes
         sinhkh(k)    = sinh(min(kwav(k)*depth(k), 100.0))
         Hmx(k)       = 0.88/kwav(k)*tanh(gamma*kwav(k)*depth(k)/0.88)
         sinhkh_ig(k) = sinh(min(kwav_ig(k)*depth(k), 50.0))
         Hmx_ig(k)    = 0.88/kwav_ig(k)*tanh(gamma*kwav_ig(k)*depth(k)/0.88)
      enddo
      !   
      do itheta = 1, ntheta
         !
         ctheta(itheta,:)    = sigm/sinh(min(2.0*kwav*depth, 50.0))*(dhdx*sin(theta(itheta)) - dhdy*cos(theta(itheta)))
         !
         ctheta_ig(itheta,:) = sigm_ig/sinh(min(2.0*kwav_ig*depth, 100.0))*(dhdx*sin(theta(itheta)) - dhdy*cos(theta(itheta)))
         !
      enddo
      !
      ! Limit unrealistic refraction speed to 1/2 pi per wave period
      !
      ctheta    = sign(1.0, ctheta)*min(abs(ctheta), sigm/4)
      ctheta_ig = sign(1.0, ctheta_ig)*min(abs(ctheta_ig), sigm_ig/4)
      !
      ! Solve the directional wave energy balance on an unstructured grid
      !
      if (.not.restart) then
         !
         ! Set energies to 0.0; note that boundary values have been set in update_boundaries
         !
         do k = 1, no_nodes
            if (inner(k)) then
               ee(:,k) = 0.0
            endif
         enddo
         !
         ee_ig = 0.0
         !
         restart = .true.
         !
      endif
      !
!      call timer(t2)
      !
      call solve_energy_balance2Dstat (x,y,no_nodes,w,ds,inner,prev,neumannconnected,       &
                                       theta,ntheta,thetamean,                                    &
                                       depth,kwav,kwav_ig,cg,cg_ig,ctheta,ctheta_ig,fw,fw_ig,Tpb,50000.,rho,snapwave_alpha,gamma,                 &
                                       H,H_ig,Dw,F,Df,thetam,sinhkh,sinhkh_ig,Hmx,Hmx_ig, ee, ee_ig, igwaves, nr_sweeps, crit, hmin)
      !
!      call timer(t3)
      !
      Fx = F*cos(thetam)
      Fy = F*sin(thetam)
      !
!      write(*,*)'computation:               ', t3 - t2, ' seconds'
      !
   end subroutine
   
   
   
   subroutine solve_energy_balance2Dstat(x,y,no_nodes,w,ds,inner,prev,neumannconnected,       &
                                         theta,ntheta,thetamean,                                    &
                                         depth,kwav,kwav_ig,cg,cg_ig,ctheta,ctheta_ig,fw,fw_ig,T,dt,rho,alfa,gamma,                 &
                                         H,H_ig,Dw,F,Df,thetam,sinhkh,sinhkh_ig,Hmx,Hmx_ig, ee, ee_ig, igwaves, nr_sweeps, crit, hmin)
   !
   implicit none
   !
   ! In/output variables and arrays
   !
   real*8, dimension(no_nodes),intent(in)           :: x,y                    ! x,y coordinates of grid
   integer, intent(in)                        :: no_nodes,ntheta              ! number of grid points, number of directions
   real*4,  dimension(2,ntheta,no_nodes),intent(in) :: w                      ! weights of upwind grid points, 2 per grid point and per wave direction
   real*4, dimension(ntheta,no_nodes), intent(in)   :: ds                     ! distance to interpolated upwind point, per grid point and direction
   real*4, intent(in)                         :: thetamean              ! mean offshore wave direction (rad)
   real*4, dimension(ntheta), intent(in)      :: theta                  ! distribution of wave angles and offshore wave energy density
!   real*4, dimension(ntheta), intent(in)      :: ee0_ig                 ! distribution of wave angles and offshore wave energy density
   logical, dimension(no_nodes), intent(in)         :: inner                  ! mask of inner grid points (not on boundary)
   integer, dimension(2,ntheta,no_nodes),intent(in) :: prev                   ! two upwind grid points per grid point and wave direction
   integer, dimension(no_nodes),intent(in)          :: neumannconnected       ! number of neumann boundary point if connected to inner point
   real*4, dimension(no_nodes), intent(in)          :: depth                     ! water depth
   real*4, dimension(no_nodes), intent(in)          :: kwav                   ! wave number
   real*4, dimension(no_nodes), intent(in)          :: kwav_ig                ! wave number
   real*4, dimension(no_nodes), intent(in)          :: cg                     ! group velocity
   real*4, dimension(ntheta,no_nodes), intent(in)   :: ctheta                 ! refractioon speed
   real*4, dimension(no_nodes), intent(in)          :: cg_ig                  ! group velocity
   real*4, dimension(ntheta,no_nodes), intent(inout)   :: ee                  ! 
   real*4, dimension(ntheta,no_nodes), intent(inout)   :: ee_ig                 ! 
   real*4, dimension(ntheta,no_nodes), intent(in)   :: ctheta_ig              ! refractioon speed
   real*4, dimension(no_nodes), intent(in)          :: fw                     ! wave friction factor
   real*4, dimension(no_nodes), intent(in)          :: fw_ig                  ! wave friction factor
   real*4, intent(in)                         :: T                      ! wave period
   real*4                                     :: T_ig                   ! wave period
   real*4, intent(in)                         :: dt                     ! time step (s)
   real*4, intent(in)                         :: rho                    ! water density
   real*4, intent(in)                         :: alfa,gamma             ! coefficients in Baldock wave breaking dissipation
   real*4, dimension(no_nodes), intent(out)         :: H                      ! wave height
   real*4, dimension(no_nodes), intent(out)         :: H_ig                      ! wave height
   real*4, dimension(no_nodes), intent(out)         :: Dw                     ! wave breaking dissipation
   real*4, dimension(no_nodes), intent(out)         :: F                      ! wave force Dw/C/rho/h
   real*4, dimension(no_nodes), intent(out)         :: Df                     ! wave friction dissipation
   real*4, dimension(no_nodes), intent(out)         :: thetam                 ! mean wave direction
   real*4, dimension(no_nodes), intent(in)          :: sinhkh                 ! sinh(k*depth)
   real*4, dimension(no_nodes), intent(in)          :: Hmx                    ! Hmax
   real*4, dimension(no_nodes), intent(in)          :: sinhkh_ig              ! sinh(k*depth)
   real*4, dimension(no_nodes), intent(in)          :: Hmx_ig                 ! Hmax
   !
   ! Local variables and arrays
   !
   integer, dimension(:), allocatable         :: ok                     ! mask for fully iterated points
   real*4                                     :: eemax,dtheta           ! maximum wave energy density, directional resolution
   real*4                                     :: uorbi
   integer                                    :: sweep,niter,iter            ! sweep number, number of iterations
   integer                                    :: k,k1,k2,i,ind(1),count,kn,itheta ! counters (k is grid index)
   integer                                    :: indint
   integer, dimension(:,:), allocatable       :: indx                   ! index for grid sorted per sweep direction
!   real*4, dimension(:,:), allocatable        :: ee,eeold               ! wave energy density, energy density previous iteration
!   real*4, dimension(:,:), allocatable        :: ee_ig                  ! wave energy density
   real*4, dimension(:,:), allocatable        :: eeold               ! wave energy density, energy density previous iteration
   real*4, dimension(:,:), allocatable        :: srcsh                  ! 
   real*4, dimension(:), allocatable          :: dee                    ! difference with energy previous iteration
   real*4, dimension(:), allocatable          :: eeprev, cgprev         ! energy density and group velocity at upwind intersection point
   real*4, dimension(:), allocatable          :: eeprev_ig, cgprev_ig         ! energy density and group velocity at upwind intersection point
   real*4, dimension(:), allocatable          :: A,B,C,R                ! coefficients in the tridiagonal matrix solved per point
   real*4, dimension(:), allocatable          :: A_ig,B_ig,C_ig,R_ig    ! coefficients in the tridiagonal matrix solved per point
   real*4, dimension(:), allocatable          :: DoverE                 ! ratio of mean wave dissipation over mean wave energy
   real*4, dimension(:), allocatable          :: DoverE_ig              ! ratio of mean wave dissipation over mean wave energy
   real*4, dimension(:), allocatable          :: E                      ! mean wave energy
   real*4, dimension(:), allocatable          :: E_ig                   ! mean wave energy
   real*4, dimension(:), allocatable          :: diff                   ! maximum difference of wave energy relative to previous iteration
   real*4, dimension(:), allocatable          :: ra                     ! coordinate in sweep direction
   integer, dimension(:), allocatable         :: tempneu                ! work array
   integer, dimension(4)                      :: shift
   real*4                                     :: pi = 4.*atan(1.0)
   real*4                                     :: g=9.81
   real*4                                     :: sigm
   real*4                                     :: hmin                   ! minimum water depth! TL: make user changeable also here according to 'snapwave_hmin' in sfincs.inp
   real*4                                     :: fac=0.25             ! underrelaxation factor for DoverE
   real*4                                     :: oneoverdt
   real*4                                     :: oneover2dtheta
   real*4                                     :: rhog8
   real*4                                     :: Dfk
   real*4                                     :: Dwk
   real*4                                     :: Ek
   real*4                                     :: Hk
   real*4                                     :: Hk0
   real*4                                     :: Fk
   real*4                                     :: percok
   real*4                                     :: error
   real*4                                     :: Dfk_ig
   real*4                                     :: Dwk_ig
   real*4                                     :: Ek_ig
   real*4                                     :: Hk_ig
   real*4                                     :: Fk_ig
   real*4                                     :: finc2ig
   real*4                                     :: shinc2ig
   real*4                                     :: fbr
   real*4                                     :: fsh
   real*4                                     :: sigm_ig
   real*4                                     :: beta
   real*4                                     :: betan
   character*20                               :: fname
   integer, save                              :: callno=1
   !
   real*4, dimension(ntheta)                  :: sinth,costh            ! distribution of wave angles and offshore wave energy density   
   !
   integer                                    :: nr_sweeps
   logical                                    :: igwaves
   real*4                                     :: crit           ! relative accuracy
   real*4                                     :: gamma_ig
   !
   integer  :: count0
   integer  :: count1
   integer  :: count_rate
   integer  :: count_max
   real     :: tloop
   !
!   igwaves   = .false.
!   nr_sweeps = 1
!   write(*,*)theta
!   write(*,*)theta*180/pi
   !
   call system_clock(count0, count_rate, count_max)
   !
   ! Allocate local arrays
   !
   allocate(ok(no_nodes))
   allocate(indx(no_nodes,nr_sweeps))
!   allocate(ee(ntheta,no_nodes))
!   allocate(ee_ig(ntheta,no_nodes))
   allocate(eeold(ntheta,no_nodes))
   allocate(dee(ntheta))
   allocate(eeprev(ntheta))
   allocate(cgprev(ntheta))
   allocate(A(ntheta))
   allocate(B(ntheta))
   allocate(C(ntheta))
   allocate(R(ntheta))
   allocate(DoverE(no_nodes))
   allocate(E(no_nodes))
   allocate(diff(no_nodes))
   allocate(ra(no_nodes))
   allocate(srcsh(ntheta,no_nodes))
   !
   if (igwaves) then
      !
      allocate(eeprev_ig(ntheta))
      allocate(cgprev_ig(ntheta))
      allocate(A_ig(ntheta))
      allocate(B_ig(ntheta))
      allocate(C_ig(ntheta))
      allocate(R_ig(ntheta))
      allocate(DoverE_ig(no_nodes))
      allocate(E_ig(no_nodes))
      !
   endif
   !   
   do itheta = 1, ntheta
      sinth(itheta) = sin(theta(itheta))
      costh(itheta) = cos(theta(itheta))
   enddo   
   !
   finc2ig = 0.20
   !
   T_ig = 7.0*T
   !   
   df   = 0.0
   dw   = 0.0
   F    = 0.0
   srcsh = 0.0
   !
   ok             = 0
   indx           = 0
   !eemax          = 1.0
   eemax          = maxval(ee)
!   dtheta         = theta(2) - theta(1)
   dtheta         = pi/ntheta
   niter          = 40
   sigm           = 2*pi/T
   oneoverdt      = 1.0/dt
   oneover2dtheta = 1.0/2.0/dtheta
   rhog8          = 0.125*rho*g
   shinc2ig       = 0.8
   thetam         = 0.0
   H              = 0.0
   H_ig           = 0.0
   !
   gamma_ig = gamma
   !   
   if (igwaves) then
      !
      DoverE_ig      = 0.0
      sigm_ig        = 2*pi/T_ig
      !
   endif
   !
   ! Sort coordinates in sweep directions
   !
   shift = [0,1,-1,2]
   !
   do sweep = 1, nr_sweeps
      !
      ra = x*cos(thetamean + 0.5*pi*shift(sweep)*pi) + y*sin(thetamean + 0.5*pi*shift(sweep)*pi)
      call hpsort_eps_epw(no_nodes, ra , indx(:, sweep), 1.0e-6)
      !
   enddo
   !
   ! Boundary condition at sea side (uniform)
   !
   do k = 1, no_nodes
      !
      if (.not.inner(k)) then         
         !
         E(k)      = sum(ee(:, k))*dtheta
         H(k)      = sqrt(8*E(k)/rho/g)
         thetam(k) = atan2(sum(ee(:, k)*sin(theta)), sum(ee(:, k)*cos(theta)))
         !
         if (igwaves) then
            !
            ee_ig(:, k)  = 0.01*ee(:,k)
            E_ig(k)      = sum(ee_ig(:, k))*dtheta
            H_ig(k)      = sqrt(8*E_ig(k)/rho/g)
            !
         endif
         !
      endif
      !
      if (inner(k)) then
         !
         if (igwaves) then
            !
            ! Compute exchange source term inc to ig waves
            !
            do itheta = 1, ntheta
               !
               k1 = prev(1, itheta, k)
               k2 = prev(2, itheta, k)
               !
               if (k1*k2>0) then
                  !
                  beta  = max((w(1, itheta, k)*(depth(k1) - depth(k)) + w(2, itheta, k)*(depth(k2) - depth(k)))/ds(itheta, k), 0.0)
                  betan = (beta/sigm_ig)*sqrt(9.81/max(depth(k), hmin))
                  fbr   = 1.0
                  fsh   = fbr*exp(-4.0*sqrt(betan))
                  !
                  cgprev(itheta) = w(1, itheta, k)*cg_ig(k1) + w(2, itheta, k)*cg_ig(k2)
                  !         
                  srcsh(itheta, k)  = - shinc2ig*fsh*((cg(k) - cgprev(itheta))/ds(itheta, k))
                  srcsh(itheta, k)  = max(srcsh(itheta, k), 0.0)
                  !
               endif
               !
            enddo
            !
         endif   
         !
         ! Make sure DoverE is filled based on previous ee
         !
         Ek       = sum(ee(:, k))*dtheta                  
         Hk       = min(sqrt(Ek/rhog8), gamma*depth(k))
         Ek       = rhog8*Hk**2
         uorbi    = pi*Hk/T/sinhkh(k)
         Dfk      = 0.28*rho*fw(k)*uorbi**3
         call baldock(g, rho, alfa, gamma, kwav(k), depth(k), Hk, T, 1, Dwk, Hmx(k))
         DoverE(k) = (Dwk + Dfk)/max(Ek, 1.0e-6)
         !
         if (igwaves) then
            !
            Ek_ig    = sum(ee_ig(:, k))*dtheta                  
            Hk_ig    = min(sqrt(Ek_ig/rhog8), gamma*depth(k))
            Ek_ig    = rhog8*Hk_ig**2
            Dfk_ig      = fw_ig(k)*0.0361*(9.81/depth(k))**(3.0/2.0)*Hk*Ek_ig
            call baldock(g, rho, alfa, gamma, kwav_ig(k), depth(k), Hk_ig, T_ig, 1, Dwk_ig, Hmx_ig(k))
            DoverE_ig(k) = (Dwk_ig + Dfk_ig)/max(Ek_ig, 1.0e-6)
            !
         endif
         !                
      endif
      !
   enddo
   !
   call system_clock(count1, count_rate, count_max)
   tloop = 1.0*(count1 - count0)/count_rate
   !
!   write(*,'(a,e14.4)')'Initialization (s) : ',tloop
   !
   ! Start iteration
   !
   call system_clock(count0, count_rate, count_max)
   !
   do iter=1,niter
      !
      sweep = mod(iter, nr_sweeps)
      !
      if (sweep==0) then
         sweep = nr_sweeps
      endif
      !
      if (sweep==1) then  !==1) then
         eeold=ee
      endif
      !
      !  Loop over all points depending on sweep direction
      !
      do count = 1, no_nodes
         !
         k = indx(count, sweep)
         !
         if (inner(k)) then
            !
            if (depth(k)>1.1*hmin) then
               !
               if (ok(k) == 0) then
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
!                     write(*,*)k,itheta,eeprev(itheta),cgprev(itheta)
                     !
                     if (igwaves) then
                        eeprev_ig(itheta) = w(1, itheta, k)*ee_ig(itheta, k1) + w(2, itheta, k)*ee_ig(itheta, k2)
                        cgprev_ig(itheta) = w(1, itheta, k)*cg_ig(k1) + w(2, itheta, k)*cg_ig(k2)
                     endif
                     !
                  enddo
                  !
                  do itheta = 2, ntheta - 1
                     !
                     A(itheta) = -ctheta(itheta - 1, k)*oneover2dtheta
                     B(itheta) = oneoverdt + cg(k)/ds(itheta,k) + DoverE(k)
                     C(itheta) = ctheta(itheta + 1, k)*oneover2dtheta
                     R(itheta) = oneoverdt*ee(itheta, k) + cgprev(itheta)*eeprev(itheta)/ds(itheta, k) - srcsh(itheta, k)*ee(itheta,k)
                     !
                  enddo
                  !
                  if (ctheta(1,k)<0) then
                     A(1) = 0.0
                     B(1) = oneoverdt - ctheta(1, k)/dtheta + cg(k)/ds(1, k) + DoverE(k)
                     C(1) = ctheta(2, k)/dtheta
                     R(1) = oneoverdt*ee(1, k) + cgprev(1)*eeprev(1)/ds(1, k) - srcsh(1, k)*ee(1, k)
                  else
                     A(1) = 0.0
                     B(1) = 1.0/dt
                     C(1) = 0.0
                     R(1) = 0.0
                  endif
                  !
                  if (ctheta(ntheta, k)>0) then
                     A(ntheta) = -ctheta(ntheta - 1, k)/dtheta
                     B(ntheta) = oneoverdt + ctheta(ntheta, k)/dtheta + cg(k)/ds(ntheta, k) + DoverE(k)
                     C(ntheta) = 0.0
                     R(ntheta) = oneoverdt*ee(ntheta,k) + cgprev(ntheta)*eeprev(ntheta)/ds(ntheta, k) - srcsh(ntheta, k)*ee(ntheta, k)
                  else
                     A(ntheta) = 0.0
                     B(ntheta) = oneoverdt
                     C(ntheta) = 0.0
                     R(ntheta) = 0.0
                  endif
                  !
                  ! Solve tridiagonal system per point
                  !
                  call solve_tridiag(A, B, C, R, ee(:,k), ntheta)
                  !
                  ee(:, k) = max(ee(:, k), 0.0)                  
                  Ek       = sum(ee(:, k))*dtheta                  
!                  write(*,*)'Ek',k,Ek,DoverE(k),dtheta
                  Hk0      = sqrt(Ek/rhog8)
                  Hk       = min(sqrt(Ek/rhog8), gamma*depth(k))
                  Ek       = rhog8*Hk**2
                  uorbi    = pi*Hk/T/sinhkh(k)
                  Dfk      = 0.28*rho*fw(k)*uorbi**3
!      if (k==23064) then
!         write(*,'(a,20e14.4)')'a',Hk,depth(k),Hmx(k)
!      endif   
                  !
                  ! Dissipation of incident waves
                  !
                  if (Hk>0.2*Hmx(k)) then
                     !
!                     call baldock(g,rho,alfa,gamma,kwav(k),depth(k),Hk,T,1,Dwk,Hmx(k))
                     call baldock(g,rho,alfa,gamma,kwav(k),depth(k),Hk0,T,1,Dwk,Hmx(k))
                     !
                     DoverE(k) = (1.0 - fac)*DoverE(k) + fac*(Dwk + Dfk)/max(Ek, 1.0e-6)
                     !
                     do itheta = 1, ntheta
                        B(itheta) = oneoverdt + cg(k)/ds(itheta,k) + DoverE(k)
                     enddo
                     !                     
                     ! Solve tridiagonal system per point
                     !                     
                     call solve_tridiag(A, B, C, R, ee(:,k), ntheta)
                     !                     
                     ee(:,k) = max(ee(:, k), 0.0)
                     Ek      = sum(ee(:, k))*dtheta
                     Hk      = min(sqrt(Ek/rhog8), gamma*depth(k))
                     Ek      = rhog8*Hk**2
                     uorbi   = pi*Hk/T/sinhkh(k)
                     Dfk     = 0.28*rho*fw(k)*uorbi**3
!      if (k==23064) then
!         write(*,'(a,20e14.4)')'b',Hk,depth(k),Hmx(k),Dwk
!      endif   
                     !                     
                     call baldock(g, rho, alfa, gamma, kwav(k), depth(k), Hk, T, 1, Dwk, Hmx(k))
                     Dw(k) = Dwk
                     !                     
                     DoverE(k) = (1.0 - fac)*DoverE(k) + fac*(Dwk + Dfk)/max(Ek, 1.0e-6)
                     !                     
                  else
                     !
                     Dwk       = 0.0
!                     Dw(k)    = 0.0
                     DoverE(k) = Dfk/max(Ek, 1.0e-6)
                     !
                  endif
                  !
                  if (igwaves) then
                     !
                     ! IG
                     !
                     do itheta = 2, ntheta - 1
                        !
                        A_ig(itheta) = -ctheta_ig(itheta - 1, k)*oneover2dtheta
                        B_ig(itheta) = oneoverdt + cg_ig(k)/ds(itheta,k) + DoverE_ig(k)
                        C_ig(itheta) = ctheta_ig(itheta + 1, k)*oneover2dtheta
                        R_ig(itheta) = oneoverdt*ee_ig(itheta, k) + cgprev_ig(itheta)*eeprev_ig(itheta)/ds(itheta, k) + srcsh(itheta, k)*ee(itheta,k)
                        !
                     enddo
                     !
                     if (ctheta_ig(1,k)<0) then
                        A_ig(1) = 0.0
                        B_ig(1) = oneoverdt - ctheta_ig(1, k)/dtheta + cg_ig(k)/ds(1, k) + DoverE_ig(k)
                        C_ig(1) = ctheta_ig(2, k)/dtheta
                        R_ig(1) = oneoverdt*ee_ig(1, k) + cgprev_ig(1)*eeprev_ig(1)/ds(1, k) + srcsh(1, k)*ee(1,k)
                     else
                        A_ig(1)=0.0
                        B_ig(1)=1.0/dt
                        C_ig(1)=0.0
                        R_ig(1)=0.0
                     endif
                     !
                     if (ctheta_ig(ntheta, k)>0) then
                        A_ig(ntheta) = -ctheta_ig(ntheta - 1, k)/dtheta
                        B_ig(ntheta) = oneoverdt + ctheta_ig(ntheta, k)/dtheta + cg_ig(k)/ds(ntheta, k) + DoverE_ig(k)
                        C_ig(ntheta) = 0.0
                        R_ig(ntheta) = oneoverdt*ee_ig(ntheta,k) + cgprev_ig(ntheta)*eeprev_ig(ntheta)/ds(ntheta,k) + srcsh(ntheta, k)*ee(ntheta,k)
                     else
                        A_ig(ntheta) = 0.0
                        B_ig(ntheta) = oneoverdt
                        C_ig(ntheta) = 0.0
                        R_ig(ntheta) = 0.0
                     endif
                     !
                     ! Solve tridiagonal system per point
                     !
                     call solve_tridiag(A_ig, B_ig, C_ig, R_ig, ee_ig(:,k), ntheta)
                     !
                     ee_ig(:, k) = max(ee_ig(:, k), 0.0)
                     Ek_ig       = sum(ee_ig(:, k))*dtheta                  
!                     Hk_ig       = min(sqrt(Ek_ig/rhog8), gamma*depth(k))
                     Hk_ig       = sqrt(Ek_ig/rhog8)
                     Ek_ig       = rhog8*Hk_ig**2
                     ! 
                     ! Bottom friction Henderson and Bowen (2002) - D = 0.015*rhow*(9.81/depth(k))**1.5*(Hk/sqrt(8.0))*Hk_ig**2/8
                     !
                     Dfk_ig      = fw_ig(k)*0.0361*(9.81/depth(k))**(3.0/2.0)*Hk*Ek_ig
                     !
                     ! Breaking of infragravity waves (should probably find another formulation for this)
                     !
                     if (Hk_ig>0.2*Hmx_ig(k)) then
                        !
                        call baldock(g, rho, alfa, gamma_ig, kwav_ig(k), depth(k), Hk_ig, T_ig, 1, Dwk_ig, Hmx_ig(k))
                        !
                        DoverE_ig(k) = (1.0 - fac)*DoverE_ig(k) + fac*(Dwk_ig + Dfk_ig)/max(Ek_ig, 1.0e-6)
                        !
                        do itheta = 1, ntheta
                           B_ig(itheta) = oneoverdt + cg_ig(k)/ds(itheta,k) + DoverE_ig(k)
                        enddo
                        !                     
                        ! Solve tridiagonal system per point
                        !                     
                        call solve_tridiag(A_ig, B_ig, C_ig, R_ig, ee_ig(:,k), ntheta)
                        !                     
                        ee_ig(:,k) = max(ee_ig(:, k), 0.0)
                        Ek_ig      = sum(ee_ig(:, k))*dtheta
!                        Hk_ig      = min(sqrt(Ek_ig/rhog8), gamma*depth(k))
                        Hk_ig      = sqrt(Ek_ig/rhog8)
                        Ek_ig      = rhog8*Hk_ig**2
                        Dfk_ig      = fw_ig(k)*0.0361*(9.81/depth(k))**(3.0/2.0)*Hk*Ek_ig
                        !                     
                        call baldock(g, rho, alfa, gamma_ig, kwav_ig(k), depth(k), Hk_ig, T_ig, 1, Dwk_ig, Hmx_ig(k))
                        !                     
                        DoverE_ig(k) = (1.0 - fac)*DoverE_ig(k) + fac*(Dwk_ig + Dfk_ig)/max(Ek_ig, 1.0)
                        !                     
                     else
                        !
                        Dwk_ig       = 0.0
                        DoverE_ig(k) = (Dfk_ig)/max(Ek_ig, 1.0)
                        !
                     endif
                  endif
                  !                     
               endif
               !
            else
               !
               ee(:, k) = 0.0
               if (igwaves) then
                  ee_ig(:, k) = 0.0
               endif
               ok(k) = 1
               !               
            endif 
            !
         endif
         !
      enddo
      !
      if (sweep==nr_sweeps) then
         !
         ! Check convergence after all sweeps
         !
         do k = 1, no_nodes
            !
            dee     = ee(:, k) - eeold(:, k)
            diff(k) = maxval(abs(dee))
            !
            if (diff(k)/eemax<crit) then
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
         call system_clock(count1, count_rate, count_max)
         tloop = 1.0*(count1 - count0)/count_rate
         !
!         write(*,'(a,i6,a,f10.5,a,f7.2,a,e14.4)')'iteration ',iter/4 ,' error = ',error,'   %ok = ',percok,' time = ',tloop
         !
         if (error<crit) then
            exit
         endif
         !
      endif
      !
   enddo
   !
   do k=1,no_nodes
      !
      ! Compute directionally integrated parameters for output
      !
      if (inner(k)) then
         if (depth(k)>1.1*hmin) then
            !
            E(k)      = sum(ee(:,k))*dtheta
!            write(*,*)k,E(k)
            H(k)      = sqrt(8*E(k)/rho/g)
!      if (k==23064) then
!         write(*,'(a,20e14.4)')'final',H(k),depth(k)
!      endif   
            thetam(k) = atan2(sum(ee(:, k)*sin(theta)),sum(ee(:, k)*cos(theta)))
            uorbi     = pi*H(k)/T/sinhkh(k)
            Df(k)     = 0.28*rho*fw(k)*uorbi**3
            !                     
            call baldock(g, rho, alfa, gamma, kwav(k), depth(k), H(k), T, 1, Dw(k), Hmx(k))
            F(k) = (Dw(k) + Df(k))*kwav(k)/sigm
            !
            if (igwaves) then
               !
               ! IG wave height
               !
               E_ig(k)      = sum(ee_ig(:,k))*dtheta
               H_ig(k)      = sqrt(8*E_ig(k)/rho/g)
               !
            endif
            !
         endif
      endif
      !
   enddo
   !
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

   subroutine baldock (rho,g,alfa,gamma,k,depth,H,T,opt,Dw,Hmax)
   !
   real*4, intent(in)                :: rho
   real*4, intent(in)                :: g
   real*4, intent(in)                :: alfa
   real*4, intent(in)                :: gamma
   real*4, intent(in)                :: k
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
!   Hmax=0.88/k*tanh(gamma*k*depth/0.88)
   Hloc=max(H,1.e-6)
   if (opt==1) then
      Dw=0.25*alfa*rho*g/T*exp(-(Hmax/Hloc)**2)*(Hmax**2+Hloc**2)
   else
      Dw=0.25*alfa*rho*g/T*exp(-(Hmax/Hloc)**2)*(Hmax**3+Hloc**3)/gamma/depth
   endif
   !
   end subroutine baldock

   
   
   subroutine disper_approx(h,T,k,n,C,Cg,no_nodes)
   integer, intent(in)                :: no_nodes
   real*4, dimension(no_nodes), intent(in)  :: h
   real*4, intent(in)                 :: T
   real*4, dimension(no_nodes), intent(out) :: k,n,C,Cg
   real*4                             :: sigma,g,pi
   !
   g     = 9.81
   pi    = 4.*atan(1.)
   sigma = 2.*pi/T
   k     = sigma**2/g*(1-exp(-(sigma*sqrt(h/g))**(2.5)))**(-0.4)
   C     = sigma/k
   n     = 0.5+k*h/sinh(min(2*k*h,50.d0))
   Cg    = n*C
   !
   end subroutine disper_approx

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
   t = dble(count)/count_rate
   end subroutine timer

       
   subroutine writebinary(filename,var,m,n)
   character(*), intent(in)           :: filename
   integer, intent(in)                :: m,n
   real*4, dimension(m,n), intent(in) :: var

   integer                            :: i,j

   open(11,file=filename,status='replace',form='unformatted')
   do j=1,n
      write(11)(var(i,j),i=1,m)
   enddo
   close(11)
   
    end subroutine writebinary
   
end module snapwave_solver 
