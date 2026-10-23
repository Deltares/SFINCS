module snapwave_vegetation
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! Effect of vegetation on incident waves: dissipation, wave force on the vegetation and
   ! nonlinear wave shape. Adapted from XBeach (vegetation.F90).
   !
   ! Subroutines / functions:
   !
   !   vegatt(...)
   !     Wave energy dissipation due to vegetation.
   !     Called from: solve_energy_balance2Dstat (snapwave_solver)
   !
   !   swvegatt(...)
   !     Short-wave dissipation due to vegetation, per vertical vegetation section.
   !     Called from: vegatt (snapwave_vegetation)
   !
   !   bulkdragcoeff(...)
   !     Bulk drag coefficient for vegetation from the Keulegan-Carpenter number.
   !     Called from: not called (the call in vegatt is commented out)
   !
   !   momeqveg(...)
   !     Wave-induced force on vegetation from the nonlinear orbital velocities.
   !     Called from: solve_energy_balance2Dstat (snapwave_solver)
   !
   !   swvegnonlin(no_nodes, kwav, depth, H, g, Trep, unl, etaw0)
   !     Nonlinear wave shape (orbital velocity and surface elevation) from the
   !     Rienecker-Fenton table.
   !     Called from: solve_energy_balance2Dstat (snapwave_solver)
   !
   ! ------------ Copyright notice SnapWave ------------
   !
   ! Copyright (C)  Stichting Deltares, 2024-2026.
   !
   ! This file is part of SnapWave as coupled with SFINCS.
   !
   ! The routines in this file are adapted from XBeach (vegetation.F90),
   ! Copyright (C)  UNESCO-IHE, Deltares and Delft University, 2011-2026,
   ! distributed under the GNU Lesser General Public License version 2.1
   ! or (at your option) any later version.
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
   ! contact: software@deltares.nl
   ! Stichting Deltares
   ! P.O. Box 177
   ! 2600 MH Delft, The Netherlands
   !
   ! ------------ Copyright notice SnapWave ------------
   !
   use sfincs_log
   !
   implicit none
   !
contains
   !
   !
   subroutine vegatt(sigm, no_nodes, kwav, no_secveg, veg_ah, veg_bstems, veg_Nstems, veg_Cd, depth, rho, g, H, Dveg) 
        !
        ! Wave energy dissipation due to vegetation.
        !
        ! Input:
        !   sigm       : mean wave frequency [rad/s]
        !   no_nodes   : number of grid nodes
        !   kwav       : wave number [rad/m]
        !   no_secveg  : number of vertical vegetation sections
        !   veg_ah     : height of the vertical vegetation sections [m]
        !   veg_bstems : stem diameter per section [m]
        !   veg_Nstems : number of stems per unit area, per section [1/m2]
        !   veg_Cd     : bulk drag coefficient per section [-]
        !   depth      : water depth [m]
        !   rho        : water density [kg/m3]
        !   g          : gravitational acceleration [m/s2]
        !   H          : wave height [m]
        !
        ! Output:
        !   Dveg       : wave dissipation by vegetation [W/m2]
        !
        ! Called from: solve_energy_balance2Dstat (snapwave_solver)
        !
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
   !
   !
    subroutine swvegatt(sigm, no_nodes, kwav, no_secveg, veg_ah, veg_bstems, veg_Nstems, veg_Cd, depth, rho, g, H, Dveg)! Short wave dissipation by vegetation
        !
        ! Short-wave dissipation due to vegetation, per vertical vegetation section.
        !
        ! Input:
        !   sigm       : mean wave frequency [rad/s]
        !   no_nodes   : number of grid nodes
        !   kwav       : wave number [rad/m]
        !   no_secveg  : number of vertical vegetation sections
        !   veg_ah     : height of the vertical vegetation sections [m]
        !   veg_bstems : stem diameter per section [m]
        !   veg_Nstems : number of stems per unit area, per section [1/m2]
        !   veg_Cd     : bulk drag coefficient per section [-]
        !   depth      : water depth [m]
        !   rho        : water density [kg/m3]
        !   g          : gravitational acceleration [m/s2]
        !   H          : wave height [m]
        !
        ! Output:
        !   Dveg       : wave dissipation by vegetation [W/m2]
        !
        ! Called from: vegatt (snapwave_vegetation)
        !
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
   !
   !
	subroutine bulkdragcoeff(ahh, m, Cdterm, no_nodes, no_secveg, depth, H, kwav, veg_bstems, sigm) !(ahh,m,i,Cdterm)
      !
		!    Michele Bendoni: subroutine to calculate bulk drag coefficient for short wave
		!    energy dissipation based on the Keulegan-Carpenter number (adapted from XBeach)
		!    Ozeren et al. (2013) or Mendez and Losada (2004)	 
      !
      ! Input:
      !   ahh        : total plant height [m]
      !   m          : index of the vegetation section
      !   no_nodes   : number of grid nodes
      !   no_secveg  : number of vertical vegetation sections
      !   depth      : water depth [m]
      !   H          : wave height [m]
      !   kwav       : wave number [rad/m]
      !   veg_bstems : stem diameter per section [m]
      !   sigm       : mean wave frequency [rad/s]
      !
      ! Output:
      !   Cdterm     : bulk drag coefficient [-]
      !
      ! Called from: not called (the call in vegatt is commented out)
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
      !	
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
		!
		Tp = 2*pi/sigm 
		!
		! Coefficient alfa
		!
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
		!
		if (myflag == 1) then
			!
			! Approach from Ozeren et al. (2013), eq?
			!
			if (KC>=10.d0) then
				Cdterm = 0.036d0+50.d0/(KC**0.926d0)
			else
				Cdterm = 0.036d0+50.d0/(10.d0**0.926d0)
			endif
         !
		elseif (myflag == 2) then
			!
			! Approach from Mendez and Losada (2004), eq. 40
			! Only applicable for Laminaria Hyperborea (kelp)???
			!
			Q = KC/(alfav**0.76d0)
			!
			if (Q>=7) then
				Cdterm = exp(-0.0138*Q)/(Q**0.3d0)
			else
				Cdterm = exp(-0.0138*7)/(7**0.3d0)
			endif
			!
		endif
		!
    end subroutine bulkdragcoeff
   !
   !
subroutine momeqveg(no_nodes, no_secveg, veg_ah, veg_bstems, veg_Nstems, veg_Cd, depth, rho, H, Trep, unl, Fvw)
    !
    ! Wave-induced force on vegetation from the nonlinear orbital velocities.
    !
    ! Input:
    !   no_nodes   : number of grid nodes
    !   no_secveg  : number of vertical vegetation sections
    !   veg_ah     : height of the vertical vegetation sections [m]
    !   veg_bstems : stem diameter per section [m]
    !   veg_Nstems : number of stems per unit area, per section [1/m2]
    !   veg_Cd     : bulk drag coefficient per section [-]
    !   depth      : water depth [m]
    !   rho        : water density [kg/m3]
    !   H          : wave height [m]
    !   Trep       : representative wave period [s]
    !   unl        : nonlinear orbital velocity over one wave period [m/s]
    !
    ! Output:
    !   Fvw        : wave-induced force on the vegetation
    !
    ! Called from: solve_energy_balance2Dstat (snapwave_solver)
    !
    ! INput: no_nodes, no_secveg, veg_ah(k,:), veg_bstems(k,:), veg_Nstems(k,:), veg_Cd(k,:), depth(k), rho, H(k), Tp(k), unl(k,:), Fvw(k)
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
        Fvgnlt = -integral / depth / Trep !> units match with F(k) m/s2
        
        Fvw = Fvw + Fvgnlt
    enddo
    !    
end subroutine momeqveg    
   !
   !
subroutine swvegnonlin(no_nodes, kwav, depth, H, g, Trep, unl, etaw0)
    !
    ! Nonlinear wave shape (orbital velocity and surface elevation) from the
    ! Rienecker-Fenton table.
    !
    ! Input:
    !   no_nodes : number of grid nodes
    !   kwav     : wave number [rad/m]
    !   depth    : water depth [m]
    !   H        : wave height [m]
    !   g        : gravitational acceleration [m/s2]
    !   Trep     : representative wave period [s]
    !
    ! Output:
    !   unl      : nonlinear orbital velocity over one wave period, per node [m/s]
    !   etaw0    : surface elevation over one wave period, per node [m]
    !
    ! Called from: solve_energy_balance2Dstat (snapwave_solver)
    !
    use snapwave_RFtable
    !
    implicit none
    !
    integer, intent(in) :: no_nodes
    real*4, dimension(no_nodes), intent(in) :: kwav 
    real*4, dimension(no_nodes), intent(in) :: depth
    real*4, dimension(no_nodes), intent(in) :: H
    real*4, intent(in) :: g
    real*4, dimension(no_nodes), intent(in) :: Trep  
    real*4, dimension(no_nodes, 50),intent(out) :: unl
    real*4, dimension(no_nodes, 50),intent(out) :: etaw0 
    
    real*4, dimension(:), save , allocatable :: h0, t0
    integer, save :: nh , nt ! save as it only needs to be done at first call
    real*4, save :: dh , dt
    real*4, dimension(50 ,8), save :: cs , sn ! MvdL: what is this fixed dimension 8 and 50?

    real*4, dimension(8) :: urf0
    real*4, dimension(50) :: urf2 , urf
    real*4, dimension(50, 8) :: urf1

    real*4, dimension(no_nodes) :: kmr , Urs , phi , w1 , w2
    real*4 :: p ,q , f0 , f1 , f2 , f3

    integer :: k, irf, ih0, it0, jrf, ih1, it1
    !
    real*4  :: pi = 4.*atan(1.0)   
    
    real*4, dimension(:,:,:), allocatable :: RFveg
    !
    allocate(RFveg(11,18,20))
    !
    ! Based on Deltares' XBeach SurfBeat' subroutine: swvegnonlin
    ! Subroutine to compute a net drag force due to wave skewness. Based on (matlab based) roller model with veggies by Ad.
    !
    ! Background:
    ! The drag force (Fveg) is a function of u*abs(u), which is zero for linear waves. For non-linear, skewed waves the
    ! depth-averaged velocity integrated over the wave period is zero. However, due to the sharp peaks and flat troughs
    ! the integral of u*abs(u) is non-zero, and can significantly reduce wave setup, or even lead to set-down (e.g. Dean & Bender,2006).
    !
    ! Here we use a method based on Rienecker & Fenton (1981), similar to the method used for onshore sediment transport due to wave asymmetry/
    ! skewness (see also morphevolution.F90 + Van Thiel de Vries Phd thesis par 6.2.3).
    !
    ! load Ad's RF-table (update for depth averaged velocities?)
    call load_RFtable(RFveg)
    !
    ! Initialize/Prepare for interpolation of RF-value from RFveg-table
    if (.not. allocated(h0)) then
        allocate(h0(no_nodes))
        allocate(t0(no_nodes))
        dh = 0.03
        dt = 1.25
        nh = floor(0.54/ dh)
        nt = floor(25 / dt )
        !construct velocity profile based on cosine/sine functions / Fourier components
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
    ! Now compute weight factors (w1,w2) for relative contribution of cosine and sine functions (for w1 = 1: only cosines ->
    ! fully skewed Stokes wave, for w2 = 1: only sines -> fully asymmetric wave) based on Ruessink.
    kmr = min(max(kwav, 0.01), 100.0)
    Urs = H / (kmr * kmr * (depth **3) )
    
    ! Compute phase and weight factors
    phi = pi /2 * (1 - tanh (0.815/(Urs **0.672) ) )
    w1 = 1 - phi /( pi /2)
    w2 = 1 - w1    
    !   
    ! Interpolate RieneckerFenton velocity from RFveg table from Ad
    ! in ftab-dimension, only read 4:11 and sum later
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
        ! Compute velocity amplitude per component
        do irf = 1, 8
            urf0(irf) = f0 * RFveg(irf + 3, ih0, it0) + f1 * RFveg(irf + 3, ih1, it0) + f2 * RFveg(irf+3, ih0, it1) + f3 * RFveg(irf + 3, ih1, it1)
        enddo    
        ! fill velocity amplitude matrix urf1([50 time points, 8 components])
        do irf = 1, 8
            urf1(:, irf) = urf0(irf)
        enddo
        !
        ! Compute velocity profile matrix per component
        urf1 = urf1 * (w1(k) * cs + w2(k) * sn )
        !
        ! Add velocity components
        urf2 = sum(urf1, 2)
        !
        ! Scale the results to get velocity profile over wave period
        unl(k,:) = urf2 * sqrt(g * depth(k) )
        etaw0(k,:) = unl(k,:)*sqrt(max(depth(k),0.d0)/g)
    enddo   
    !   
end subroutine swvegnonlin
   !
   !
end module snapwave_vegetation
