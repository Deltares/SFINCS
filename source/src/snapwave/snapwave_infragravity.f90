module snapwave_infragravity

    use sfincs_log
    
    !TL: Module to calculate IG boundary conditions at offshore boundary input points, following Herbers et al. (1994) 
    ! and based on implementations in Matlab (secordspec) and XBeach (waveparams.f90)
    contains   
   
    !-----------------------------Main subroutine----------------------------!
    
    subroutine determine_ig_bc(x_bwv, y_bwv, hsinc, tpinc, ds, jonswapgam, depth, Tinc2ig, tpig_opt, hsig, tpig) 
    ! (input, input, input, input, input, input, input, output, output)
    !  
    ! Build boundwave offshore spectrum, and determine Hig0 and Tpig0, using Herbers 1994 as in XBeach implementation
    !
    ! Input hsinc is Hm0, not Hrms
    ! Output hsig is therefore also Hm0 (and expected like that in E0_ig   = 0.0625*rho*g*hst_bwv_ig(ib)**2 in subroutine update_boundary_points)
    !
    ! Output tpig is estimated in different ways based on calculated spectrum > choice which one to use
    !
    ! Input ds is already in rad
    !
    !use snapwave_data
    !   
    implicit none
    !
    real*4, intent(in)    :: hsinc, tpinc, ds, jonswapgam, Tinc2ig
    real*8, intent(in)    :: x_bwv, y_bwv
    integer, intent(in)   :: tpig_opt
    real*4, intent(inout) :: depth
    real*4, intent(out)   :: hsig, tpig
    !
    real*4                :: pi, scoeff
    real*4                :: Tm01, Tm10, Tp, Tpsmooth
    integer               :: correctHm0
    !
    correctHm0 = 0 ! Choice between correcting Hm0 in build_jonswap if build 2D Vardens spectrum too low (1) or not (default, 0)
    ! TL: Question - should we make this user definable?
    !
    pi    = 4.*atan(1.)
    !
    ! Convert wave spreading in degrees (input) to S
    scoeff = (2/ds**2) - 1
    !  
    ! Call function that calculates Hig0 following Herbers, as also implemented in XBeach and secordspec2 in Matlab
    ! Loosely based on 3 step calculation in waveparams.F90 of XBeach (build_jonswap, build_etdir, build_boundw), here all in 1 subroutine calculate_herbers
    !
    if (depth < 5.0) then
        !
	    write(logstr,*)'ERROR SnapWave - depth at boundary input point ',x_bwv, y_bwv,' dropped below 5 m: ',depth
        call write_log(logstr, 1)     
        !
        write(logstr,*)'This might lead to large values of Hm0ig as bc, especially when directional spreading is low! Please specify input in deeper water. '
        call write_log(logstr, 1)   
        !
        write(logstr,*)'Depth set back to 5 meters for stability, simulation will continue.'
        call write_log(logstr, 1)   
        !        
        depth = 5.0
        !
    endif	
    !
    call compute_herbers(hsig, Tm01, Tm10, Tp, Tpsmooth, hsinc, tpinc, scoeff, jonswapgam, depth, correctHm0) ![out,out,out,out,out, in,in,in,in,in,in]
    !   
    ! Catch NaN values (if depth=0 probably) or unrealistically large values above 2 meters
    if (hsig < 0.0) then
        !
	    write(logstr,*)'DEBUG SnapWave - computed hm0ig at boundary dropped below 0 m: ',hsig, ' and is therefore limited back to 0 m!'
        call write_log(logstr, 1)        
	    hsig = max(hsig, 0.0)
        !
    endif
    !
    if (hsig > 3.0) then
        !
	    write(logstr,*)'DEBUG SnapWave - computed hm0ig at boundary exceeds 3 meter: ',hsig, ' - please check whether this might be realistic!'
        call write_log(logstr, 1)        	    
        !write(*,*)'DEBUG - computed hm0ig at boundary exceeds 3 meter: ',hsig, ' and is therefore limited back to 3 m!'
	    !hsig = min(hsig, 3.0)
        !
    endif
    !
    ! Choose what wave period option value for IG to choose:
    ! Options: 1=Tm01, 2=Tpsmooth, 3=Tp, 4=Tm-1,0
    !
    if (tpig_opt == 1) then ! Default
        !
        tpig = Tm01
        !
    elseif (tpig_opt == 2) then
        !        
        tpig = Tpsmooth
        !
    elseif (tpig_opt == 3) then
        !        
        tpig = Tp            
        !
    elseif (tpig_opt == 4) then
        !        
        tpig = Tm10 ! Tm-1,0
        !    
    endif
    !
    ! Check on ratio tpig/tpinc whether it is deemed realistic
    if (tpig/tpinc < 2.0) then
	    write(logstr,*)'DEBUG SnapWave - computed tpig/tpinc ratio at offshore boundary dropped below 2 and might be unrealistic! value: ',tpig/tpinc
        call write_log(logstr, 1)        
    elseif (tpig/tpinc > 20.0) then
	    write(logstr,*)'DEBUG SnapWave - computed tpig/tpinc ratio at offshore boundary increased above 20 and might be unrealistic! value: ',tpig/tpinc
        call write_log(logstr, 1)        
    endif	         
    !
    end subroutine
   
    !-------------------------Supporting subroutines-------------------------!
    
    subroutine compute_herbers(hsig, Tm01, Tm10, Tp, Tpsmooth, hsinc, tpinc, scoeff, jonswapgam, depth, correctHm0)
    !
    use interp        
    !
    implicit none
    !
    ! Incoming and outgoing variables
    real*4, intent(in)                      :: hsinc, tpinc, scoeff, jonswapgam, depth
    integer, intent(in)                     :: correctHm0
    real*4, intent(out)                     :: hsig, Tm01, Tm10, Tp, Tpsmooth    
    !
    ! Internal variables - for part 1: build_jonswap
    real*4                                  :: dfj, fp, fnyq, dang, iang, angtemp, hsinc_check1, df
    real*4, dimension(:), allocatable       :: temp, x, y, f, ang, Dd, Sf, findline, fgen    
    real*4, dimension(:,:), allocatable     :: S_array    
    integer                                 :: i=0, ii, nang, nfreq, peakf
    integer                                 :: firstp, lastp, M, K    
    real*4                                  :: pi, g, sprdthr
    !
    ! Internal variables - for part 2: build_etdir    
    real*4                                  :: kmax, pp
    real*4, dimension(400)                  :: P0, kk, phase, Sf0, Sf0org,S0org     !K=400 - OR allocate size later once K is defined
    real*4, dimension(202)                  :: P1, ang1   ! Are size 200+2        
    real*4, dimension(:), allocatable       :: theta, theta0, temp2, S0, Sn, dthetafin
    real*4, dimension(:,:), allocatable     :: Snew_array, Snew1_array, Ddnew        
    real*4                                  :: hm0now, s1, s2, modf, modang, hsinc_check2, hsinc_check2tmp    
    integer                                 :: F2, stepf, stepang, jj
    !
    ! Internal variables - for part 3: build_boundw    
    real*4                                  :: deltaf
    real*4, dimension(:), allocatable       :: w1, k1
    real*4, dimension(:), allocatable       :: Ebnd, fbnd
    real*4, dimension(:), allocatable       :: term1, term2, term2new, dif, chk1, chk2
    real*4, dimension(:,:), allocatable     :: Eforc, D, deltheta!, KKx, KKy, theta3
    real*4, dimension(:,:), allocatable     :: dphi3, k3, cg3!, Abnd   
    integer                                 :: valdensmax
    !
    ! Constants
    pi  = 4.*atan(1.)    
    g   = 9.81
    !   
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! TL: Here we go to the code of XBeach subroutine 'build_jonswap' from waveparams.F90 
    !
    dfj = 0.001
    fp = 1 / tpinc
    fnyq = fp * 2
    !
    ! Define number of frequency bins by defining an array of the necessary length using the Nyquist frequency and frequency step size
    allocate(temp(ceiling((fnyq-dfj)/dfj)))
    temp=(/(i,i=1,size(temp))/)   
    !
    ! Define array with actual equidistant frequency bins
    allocate(f(size(temp)))
    f=temp*dfj
    deallocate (temp)
    !
    ! Determine frequency bins relative to peak frequency
    allocate(x(size(f)))
    x=f/fp
    !
    ! Calculate unscaled and non-directional JONSWAP spectrum using peak-enhancement factor and pre-determined frequency bins
    allocate(y(size(f)))
    call jonswapgk(x,jonswapgam,y)
    !
    ! Determine scaled and non-directional JONSWAP spectrum using the JONSWAP characteristics
    y = (hsinc/(4.d0*sqrt(sum(y)*dfj)))**2*y
    deallocate (x)  
    ! 
    ! Define 200 directions relative to main angle running from -pi to pi
  
    allocate(temp(200))
    allocate(ang(200))    
    temp=(/(i,i=0,200-1)/) ! as in Matlab version have 200 values    
    ang=temp*((2*pi)/200.d0)-pi
    deallocate (temp)    
    !
    ! Determine directional step size: pi/200
    dang=ang(2)-ang(1)
    !
    ! Define 200 directional bins accordingly
    allocate (Dd(size(ang)))    
    !
    ! Calculate directional spreading based on cosine law    
    Dd = cos(ang/2)**(2*nint(scoeff) ) ! Robert: apparently nint is needed here, else MATH domain error TL: also the case for us
    !
    ! Scale directional spreading to have a surface of unity by dividing by it's own surface
    Dd = Dd / (sum(Dd)*dang)    
    !
    ! Define number of directional and frequency bins
    nang=size(ang)
    nfreq=size(y)        
    !
    ! Define two-dimensional variance density spectrum array and distribute variance density for each frequency over directional bins
    allocate(S_array(nfreq,nang))
    !
    do i=1,nang
        do ii=1,nfreq
            S_array(ii,i)=y(ii)*Dd(i)
        end do
    end do
    deallocate (y)    
    !
    ! Add check on Hm0 of created two-dimensional variance density spectrum :
    hsinc_check1 = 4*sqrt(sum(S_array)*dfj*dang)
    if (abs(hsinc - hsinc_check1) > 0.01) then
        write(logstr,*)'WARNING - computed Hm0,inc of 2D var dens spectrum differs from input! see subroutine build_jonswap in determine_ig_bc in module snapwave_boundaries.f90'
        call write_log(logstr, 1)           
        write(logstr,*)'Input: ',hsinc,' , computed: ', hsinc_check1
        call write_log(logstr, 1)           
    endif
    !
    ! Back integrate two-dimensional variance density spectrum over directions    
    allocate(Sf(size(f)))
    Sf = sum(S_array, DIM = 2)*dang    
    !
    ! Determine frequencies around peak frequency of one-dimensional non-directional variance density spectrum, based on factor sprdthr (0.08 by default)
    allocate(findline(size(Sf)))
    findline = 0.d0    
    sprdthr = 0.08
    call frange(Sf,firstp,lastp,findline,sprdthr)
    !findline = findloc(Sf, Sf>0.08*maxval(Sf)) ! TL: should probably be possible to do this more easily
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! TL: Here we go to the code of XBeach subroutine 'build_etdir' from waveparams.F90 
    !
    K = 400 ! as in Matlab script
    !
    ! Determine number of frequencies in discrete variance density spectrum to be included    
    M = int(sum(findline)) ! number of points in frequency range around peak
    !
    ! Define number of wave components to be used
    allocate (temp(K))
    temp=(/(i,i=0,K-1)/)
    !  
    ! Select equidistant wave components between the earlier selected range of frequencies around the peak frequency based on sprdthr
    allocate(fgen(K))
    fgen = temp * ((f(lastp) - f(firstp)) / (K-1)) + f(firstp)
    deallocate(temp)
    !
    ! Determine equidistant frequency step size
    df = (fgen(K) - fgen(1)) / (dble(K) - 1.d0)    
    !        
    ! 2D interpolate from S_array of (f,ang) grid to (fgen,ang)
    !
    allocate(Snew_array(K,size(ang)))    
    do ii=1,K
        do jj=1,size(ang)
            !
            call linear_interp_2d_real4(f,size(f), &    ! input frequency, size
                ang,size(ang), &                        ! input angles, size
                S_array, &                              ! input 2D variance density
                fgen(ii),ang(jj), &                     ! output frequency/angle
                Snew_array(ii,jj), &                    ! output variance density
                'interp',0.0)                           ! method and exception return value    
        enddo        
    enddo
    !
    Sn = sum(Snew_array, DIM = 2) * dang ! TL: similar to Sf    
    !
    allocate(Ddnew(K,size(ang)))    
    !
    do ii=1,K
        Ddnew(ii,:) = Snew_array(ii,:) / Sn(ii)
    enddo
    !
    ! Define random number between 0 and 1 for each wave component
    do i=1,K
        call RANDOM_NUMBER(P0(i)) 
        P0(i)=0.99*P0(i)+0.01/2 ! TL: as in Matlab: Define random number between ~0.025 and 0975 for each wave component
    enddo
    !
    ! Define direction for each wave component based on random number and linear interpolation of the probability density function
    allocate(theta(K))
    !
    if (scoeff >= 1000.0) then !longcrested waves
        theta = 0 
    else
        do ii=1,size(P0)
            !
            P1(1) = 0                   ! P1 is a function between 0 and 1, that looks like a cdf function
            ang1(1) = minval(ang)-dang  ! ang1 is a linear function between -pi and pi
            !
            do jj=1,size(ang)
                P1(jj+1) = sum(Ddnew(ii,1:jj)) * dang + 0.00001 * jj
                ang1(jj+1) = ang(jj)
            enddo            
            !
            P1(size(ang)+2) = P1(size(ang)+1) + 0.0001
            ang1(size(ang)+2) = ang(size(ang))+dang            
            ! TL: equivalent to Matlab implementation: P1 = [0 P P(length(ang))+0.0001] & ang1 = [min(ang)-dang ang max(ang)+dang];
            !                        
            ! Find the angle for each frequency component            
            call linear_interp_real4(P1,ang1,size(P1),P0(ii),theta(ii),F2)  ! theta(ii) is the output, finding the corresponding angle per random P0 value
            !
            ! TL: as in Matlab - for interpolation purposes:
            if (theta(ii) < minval(ang1)) then
                theta(ii) = theta(ii) + 2*pi
            endif
            if (theta(ii) > maxval(ang1)) then
                theta(ii) = theta(ii) - 2*pi
            endif            
            !
        enddo
    endif    
    !    
    ! TL: equivalent to Matlab - Snew = [S(:,noang) S S(:,1)]; % this is for interpolation purposes
    allocate(Snew1_array(K,size(ang)+2))    
    !
    do ii=1,K ! TL: =size(fgen)
        !
        Snew1_array(ii,1) = Snew_array(ii,size(ang))  
        Snew1_array(ii,size(ang)+2) = Snew_array(ii,1)
        !
        do jj=1,size(ang)
            !
            Snew1_array(ii,jj+1) = Snew_array(ii,jj)            
            !
        enddo        
    enddo
    !
    ! Determine variance density spectrum values for all relevant wave components
    ! around the peak frequency by interpolation of two-dimensional variance density spectrum array.    
    !
    ! 2D interpolate from S_array of (f,ang) grid to (fgen,ang)    
    allocate(S0(K))     
    !        
    do ii=1,K ! TL: =size(fgen)    
        !
        call linear_interp_2d_real4(fgen,size(fgen), &  ! input frequency, size
            ang1,size(ang1), &                          ! input angles, size
            Snew1_array, &                              ! input 2D variance density
            fgen(ii),theta(ii), &                       ! output frequency/angle
            S0(ii), &                                   ! output variance density
            'interp',0.0)                               ! method and exception return value    
        !
    enddo    
    !
    allocate(dthetafin(K))         
    dthetafin = Sn/S0					
    !
    ! Determine significant wave height using Hm0 = 4*sqrt(m0) using the one-dimensional non-directional variance density spectrum    
    hsinc_check2 = 4*sqrt(sum(S0 * dthetafin * df))
    hsinc_check2tmp = 4*sqrt(sum(Sn) * df)    
    !
    ! Add check on Hm0 of created two-dimensional variance density spectrum :
    !if (abs(hsinc_check2 - hsinc_check1) > 0.01) then
    if (abs(hsinc_check2 - hsinc_check2tmp) > 0.01) then        
        write(logstr,*)'WARNING SnapWave - computed Hm0,inc of 2D var dens spectrum differs from input! see subroutine build_jonswap in determine_ig_bc in module snapwave_boundaries.f90'
        call write_log(logstr, 1)        
        write(logstr,*)'Newly computed in part 2: ',hsinc_check2,' , while computed before: ', hsinc_check2tmp, ' and input was: ', hsinc
        call write_log(logstr, 1)           
        !write(*,*)'Newly computed in part 2: ',hsinc_check2,' , computed in part 1: ', hsinc_check1        
    endif    
    !
    ! Correct spectra for wave height > TL: in Xbeach, not Matlab, option added for now
    if (correctHm0 == 1) then
        !
        S0 = (hsinc/hsinc_check2)**2 * S0
        Sn = (hsinc/hsinc_check2)**2 * Sn
        dthetafin = Sn/S0					    
        !
        ! To check again:
        hsinc_check2 = 4*sqrt(sum(S0 * dthetafin * df))
        hsinc_check2tmp = 4*sqrt(sum(Sn) * df)                 
        !
        write(logstr,*)'DEBUG SnapWave - Hm0 of vardens corrected to: ',hsinc_check2,' and ', hsinc_check2tmp, ' so it is close to input: ', hsinc
        call write_log(logstr, 1)           
        !        
    endif    
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! TL: Here we go to the code of XBeach subroutine 'build_boundw' from waveparams.F90 
    !
    ! Allocate two-dimensional variables for all combinations of interacting wave components to be filled triangular
    allocate(Eforc(K-1,K))      ! Herbers et al. (1994) eq. 1 - (rows = difference frequency, columns is interaction)
    allocate(D(K-1,K))          ! Herbers eq. A5.
    allocate(deltheta(K-1,K))   ! Difference angle between two primary wave components
    allocate(k3(K-1,K))         ! Wavenumber of difference wave
    allocate(cg3(K-1,K))
    !
    ! Allocate variables for angular velocity and wave numbers for wave components
    allocate(w1(size(fgen)))    ! Radial frequency of primary waves
    allocate(k1(size(fgen)))    ! Wave numbers of primary waves
    !
    ! Initialize variables as zero
    Eforc = 0
    D = 0
    deltheta = 0
    k3 = 0
    cg3 = 0
    w1=0
    k1=0    
    !
    ! Steps indepent of loop, that can be done once a priori
    ! Determine angular velocity of primary waves
    w1=2*pi*fgen
    !
    ! Determine wave numbers of primary waves
    call bc_disper(k1,w1,size(w1),depth,g)    
    !
    ! Determine for each wave component interactions with all other wave components
    ! as far as not processed yet (each loop step the number of interactions thus decrease with one)
    !
    do m=1,K-1
        !
        ! Determine difference frequency
        deltaf=m*df
        !
        ! Determine difference angles (pi already added)
        deltheta(m,1:K-m) = abs(theta(m+1:K)-theta(1:K-m))+pi
        !
        ! Determine difference wave numbers according to Van Dongeren et al. 2003 eq. 19
        k3(m,1:K-m) =sqrt(k1(1:K-m)**2+k1(m+1:K)**2+2*k1(1:K-m)*k1(m+1:K)*cos(deltheta(m,1:K-m))) 
        ! Note, dcos is for double precision, use cos for single precision (real*4)
        !
        ! Determine group velocity of difference waves
        cg3(m,1:K-m) = 2.d0*pi*deltaf/k3(m,1:K-m)
        !
        ! XBeach: Make sure that we don't blow up bound long wave when offshore boundary is too close to shore > not needed for us?
        !cg3(m,1:K-m) = min(cg3(m,1:K-m),par%nmax*sqrt(g/k3(m,1:K-m)*tanh(k3(m,1:K-m)*depth)))
        !
        ! Determine difference-interaction coefficient according to Herbers 1994 eq. A5
        allocate(term1(K-m),term2(K-m),term2new(K-m),dif(K-m),chk1(K-m),chk2(K-m))
        !
        term1 = (-w1(1:K-m))*w1(m+1:K)
        term2 = (-w1(1:K-m))+w1(m+1:K)
        term2new = cg3(m,1:K-m)*k3(m,1:K-m)
        dif = (abs(term2-term2new))
        !
        chk1  = cosh(k1(1:K-m)*depth)
        chk2  = cosh(k1(m+1:K)*depth)
        !
        D(m,1:K-m) = -g*k1(1:K-m)*k1(m+1:K)*cos(deltheta(m,1:K-m))/2.d0/term1+g*term2*(chk1*chk2)/ &
        ((g*k3(m,1:K-m)*tanh(k3(m,1:K-m)*depth)-(term2new)**2)*term1*cosh(k3(m,1:K-m)*depth))* &
        (term2*((term1)**2/g/g - k1(1:K-m)*k1(m+1:K)*cos(deltheta(m,1:K-m))) &
        - 0.50d0*((-w1(1:K-m))*k1(m+1:K)**2/(chk2**2)+w1(m+1:K)*k1(1:K-m)**2/(chk1**2)))
        !
        deallocate(term1,term2,term2new,dif,chk1,chk2)
        !
        ! Correct for surface elevation input and output instead of bottom pressure so it is consistent with Van Dongeren et al 2003 eq. 18
        D(m,1:K-m) = D(m,1:K-m)*cosh(k3(m,1:K-m)*depth)/(cosh(k1(1:K-m)*depth)*cosh(k1(m+1:K)*depth))
        !
        ! Exclude interactions with components smaller than or equal to current component according to lower limit Herbers 1994 eq. 1
        where(fgen<=m*df)
        D(m,:)=0.d0                                           ! Bas: redundant with initial determination of D ??
        endwhere
        !
        ! Determine energy of bound long wave according to Herbers 1994 eq. 1 based
        ! on difference-interaction coefficient and energy density spectra of
        ! primary waves
        Eforc(m,1:K-m) = 2*D(m,1:K-m)**2*S0(1:K-m)*S0(m+1:K)*dthetafin(1:K-m)*dthetafin(m+1:K)*df
        !
    end do
    !
    ! Allocate variables for energy of bound long wave
    allocate(Ebnd(K-1))
    !
    ! Sum over the components to get total forced wave at diff freq
    Ebnd = sum(Eforc,2)
    !
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Determine final parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    !
    ! What we actually want for SnapWave offshore IG bc: Hm0ig 
    hsig = 4*sqrt(sum(Ebnd)*df)   
    !
    ! Calculate representative value for IG wave period
    allocate (temp(K))
    temp=(/(i,i=0,K-1)/)
    !  
    ! Select equidistant wave components between the earlier selected range of frequencies around the peak frequency based on sprdthr
    allocate(fbnd(K-1))
    fbnd = temp * df
    deallocate(temp)
    !
    ! Calculate (mean) wave period based on one-dimensional non-directional variance density spectrum and factor trepfac
    call tpDcalc(Ebnd,fbnd,0.01,Tm01,Tm10,Tp,Tpsmooth)! [in,in,in,out,out,out,out]
    !XBeach default for trepfac = 0.01
    !
    end subroutine
            
    ! -----------------------------------------------------------
    ! --------- JONSWAP  unscaled JONSWAP spectrum --------------
    ! -----------------(used by compute_herbers)-----------------
    subroutine jonswapgk(x,gam,y)

    implicit none
    ! Required input: - x           : nondimensional frequency, divided by the peak frequency
    !                 - gam         : peak enhancement factor, optional parameter (DEFAULT 3.3)
    !                 - y is output : nondimensional relative spectral density, equal to one at the peak

    real*4, INTENT(IN)                  :: gam
    real*4,dimension(:), INTENT(IN)     :: x
    real*4,dimension(:), INTENT(INOUT)  :: y

    ! Internal variables
    real*4,dimension(size(x))           :: xa, sigma, fac1, fac2, fac3, temp

    xa=abs(x)

    where (xa==0)
        xa=1e-20
    end where

    sigma=xa

    where (sigma<1.)
        sigma=0.07
    end where

    where (sigma>=1.)
        sigma=0.09
    end where

    temp=0*xa+1

    fac1=xa**(-5)
    fac2=exp(-1.25*(xa**(-4)))
    fac3=(gam*temp)**(exp(-((xa-1)**2)/(2*(sigma**2))))

    y=fac1*fac2*fac3
    y=y/maxval(y)

    return

    end subroutine jonswapgk   
    
   ! -----------------------------------------------------------
   ! ---- Small subroutine to determine f-range round peak -----
   ! -----------------(used by compute_herbers)-----------------
   subroutine frange(Sf,firstp,lastp,findlineout,sprdthr)

      implicit none

      real*4, dimension(:), intent(in)        :: Sf
      integer, intent(out)                    :: firstp, lastp

      real*4, dimension(:), intent(out)       :: findlineout
      real*4, dimension(:),allocatable        :: temp, findline
      integer                                 :: i = 0
      real*4,intent(in)                       :: sprdthr
      
      allocate(findline(size(Sf)))
      findline=0*Sf                           ! find frequency range around peak

      where (Sf>sprdthr*maxval(Sf))
         findline=1
      end where

      firstp=maxval(maxloc(findline))         ! Picks the first "1" in temp

      allocate (temp(size(findline)))
      temp=(/(i,i=1,size(findline))/)
      lastp=maxval(maxloc(temp*findline))     ! Picks the last "1" in temp

      findlineout=findline
      deallocate(temp, findline)

   end subroutine frange    
    
   ! --------------------------------------------------------------
   ! --------------------- Dispersion relation --------------------
   ! -------------------(used by compute_herbers)------------------
   subroutine bc_disper(k1,w1,m,h,g)
      !          k  = wave number             (2 * pi / wave length)
      !          w  = wave angular frequency  (2 * pi / wave period)
      !          m  = size k and w vectors
      !          h  = water depth
      !          g  = gravitational acceleration constant, optional (DEFAULT 9.81)
      !
      !          absolute error in k*h < 5.0e-16 for all k*h
      !
      !
      !          original Matlab code by: G. Klopman, Delft Hydraulics, 6 Dec 1994

      integer, intent(in)                     :: m
      real*4,dimension(m),intent(in)          :: w1
      real*4,dimension(m),intent(out)         :: k1
      real*4, intent(in)                      :: h, g

      ! internal variables

      real*4,dimension(m)                     :: w2,q,thq,thq2,a,b,c,arg,sign
      integer                                 :: j
      real*4                                  :: hu

      w2 = w1**2*(h/g)
      q = w2/(1.0d0-exp(-(w2**(5.0d0/4.0d0))))**(2.0d0/5.0d0)

      do j=1,4
         thq  = tanh(q)
         thq2 = 1.0d0-thq**2
         a    = (1.0d0-q*thq)*thq2
         b    = thq + q*thq2
         c    = q*thq-w2
         where (abs(a*c)<(b**2*1.0e-8))
            arg = -c/b
         elsewhere
            arg  = (b**2)-4.0d0*a*c
            arg  = (-b + sqrt(arg))/(2.0d0*a)
         endwhere
         q    = q+arg
      end do

      where (w1>0.0d0)
         sign=1.0d0
      endwhere

      where (w1==0.0d0)
         sign=0.0d0
      endwhere

      where (w1<0.0d0)
         sign=-1.0d0
      endwhere

      k1 = sign*q/h

      where (k1==huge(hu))
         k1=0.0d0
      endwhere

      where (k1==-1.0d0*huge(hu))
         k1=0.0d0
      endwhere

      return

   end subroutine bc_disper   
   
   ! -----------------------------------------------------------
   ! ----------- Small subroutine to determine tpD -------------
   ! -----------------(used by compute_herbers)-----------------
   !
   ! TL: Slightly adapted from XBeach version to give more output,
   ! give back Tm01, Tm-1,0 and Tp
   subroutine tpDcalc(Sf,f,trepfac,Tm01,Tm10,Tp,Tpsmooth) ! [in,in,in,out,out,out,out]

      implicit none

      ! In:
      real*4, dimension(:), intent(in)        :: Sf, f
      real*4, intent(in)                      :: trepfac
      !
      ! Out:
      real*4, intent(out)                     :: Tm01
      real*4, intent(out)                     :: Tm10
      real*4, intent(out)                     :: Tp   
      real*4, intent(out)                     :: Tpsmooth            
      !
      ! Local:
      real*4, dimension(:),allocatable        :: temp, Sfsmooth
      integer                                 :: id  
      !
      allocate(temp(size(Sf)))
      temp=0.d0
      where (Sf>=trepfac*maxval(Sf))
         temp=1.d0
      end where
      !
      Tm01=sum(temp*Sf)/sum(temp*Sf*f)   ! Tm01
      !
      Tm10 = sum(temp*Sf/f)/sum(temp*Sf) ! Tm-1,0
      !
      id = int(MAXLOC(Sf, dim=1))
      Tp = 1.0 / f(id)                   ! Tp
      !    
      ! And Tp over smoothed spectrum:
      allocate(Sfsmooth(size(Sf)))
      !
      Sfsmooth = 0.0      
      do id=2,size(Sf)-1
         Sfsmooth(id) = (0.5*Sf(id-1)+Sf(id)+0.5*Sf(id+1))/2
      enddo               
      !
      id = int(MAXLOC(Sfsmooth, dim=1))
      Tpsmooth = 1.0 / f(id)                   ! Tp,smooth      
      !
   end subroutine tpDcalc
               
end module    