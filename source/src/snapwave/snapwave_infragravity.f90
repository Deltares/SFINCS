module snapwave_infragravity
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! Infragravity wave boundary conditions at the offshore boundary points, following
   ! Herbers et al. (1994), and the infragravity source/sink term in the energy balance
   ! (Leijnse et al., 2024).
   !
   ! Subroutines / functions:
   !
   !   determine_ig_bc(...)
   !     Determine the infragravity wave height and period at a boundary point from the
   !     incident wave conditions.
   !     Called from: update_boundary_points (snapwave_boundaries)
   !
   !   compute_herbers(...)
   !     Compute the bound infragravity wave spectrum from a JONSWAP spectrum with the
   !     Herbers et al. (1994) method.
   !     Called from: determine_ig_bc (snapwave_infragravity)
   !
   !   jonswapgk(x, gam, y)
   !     Evaluate the JONSWAP spectral shape for peak enhancement factor gam.
   !     Called from: compute_herbers (snapwave_infragravity)
   !
   !   frange(Sf, firstp, lastp, findlineout, sprdthr)
   !     Find the range of frequencies where the spectrum exceeds the threshold sprdthr.
   !     Called from: compute_herbers (snapwave_infragravity)
   !
   !   bc_disper(k1, w1, m, h, g)
   !     Solve the linear dispersion relation for the wave numbers k1 at frequencies w1.
   !     Called from: compute_herbers (snapwave_infragravity)
   !
   !   tpDcalc(Sf, f, trepfac, Tm01, Tm10, Tp, Tpsmooth)
   !     Compute the wave periods Tm01, Tm-1,0, Tp and smoothed Tp from a frequency
   !     spectrum.
   !     Called from: compute_herbers (snapwave_infragravity)
   !
   !   determine_infragravity_source_sink_term(...)
   !     Determine the infragravity source/sink term from the incident wave energy (Leijnse
   !     et al., 2024).
   !     Called from: solve_energy_balance2Dstat (snapwave_solver)
   !
   !   estimate_shoaling_parameter_alphaig(beta, gam, alphaig)
   !     Estimate the infragravity shoaling parameter alphaig from the bed slope and
   !     relative depth.
   !     Called from: determine_infragravity_source_sink_term (snapwave_infragravity)
   !
   ! ------------ Copyright notice SnapWave ------------
   !
   ! Copyright (C)  Stichting Deltares, 2023-2026.
   !
   ! This file is part of SnapWave as coupled with SFINCS.
   !
   ! Parts of this file are adapted from XBeach (waveparams.F90),
   ! Copyright (C)  UNESCO-IHE, Deltares and Delft University, 2007-2026,
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
    
    !TL: Module to calculate IG boundary conditions at offshore boundary input points, following Herbers et al. (1994) 
    ! and based on implementations in Matlab (secordspec) and XBeach (waveparams.f90)
    contains   
   
    !-----------------------------Main subroutine----------------------------!
    
    subroutine determine_ig_bc(x_bwv, y_bwv, hsinc, tpinc, ds, jonswapgam, depth, Tinc2ig, tpig_opt, hsig, tpig) 
    !
    ! Determine the infragravity wave height and period at a boundary point from the
    ! incident wave conditions.
    !
    ! Input:
    !   x_bwv      : x coordinate of the boundary point
    !   y_bwv      : y coordinate of the boundary point
    !   hsinc      : incident wave height Hm0 [m]
    !   tpinc      : incident peak period [s]
    !   ds         : directional spreading [rad]
    !   jonswapgam : JONSWAP peak enhancement factor [-]
    !   Tinc2ig    : ratio of infragravity to incident period, used when positive
    !   tpig_opt   : which period of the infragravity spectrum to return as tpig
    !
    ! Input/output:
    !   depth      : water depth at the boundary point [m]
    !
    ! Output:
    !   hsig       : infragravity wave height [m]
    !   tpig       : infragravity wave period [s]
    !
    ! Called from: update_boundary_points (snapwave_boundaries)
    !
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
    ! Loosely based on 3 step calculation in waveparams.F90 of XBeach (build_jonswap, build_etdir, build_boundw), here all in 1 subroutine compute_herbers
    !
    if (hsinc / depth > 0.5) then
       !
	   write(logstr, *)'ERROR SnapWave - Hs over depth at boundary input point ', x_bwv, ',', y_bwv,' is larger then 0.5: ', hsinc / depth
       call write_log(logstr, 0)   
	   write(logstr, *)'This may lead to large values of Hm0ig as bc, especially when directional spreading is low! Please specify input in deeper water.'
       call write_log(logstr, 0)   
       write(logstr,*)'Depth set back to 2.0 * hsinc meters for stability, simulation will continue.'
       call write_log(logstr, 0)   
       !        
       depth = 2.0 * hsinc
       !
    elseif (depth > 200.0) then
       !
       ! Limit depth to 200 m. Larger depth can result in NaNs. @Tim, why?
       !
       depth = 200.0
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
    elseif (hsig > 3.0) then
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
    ! Compute the bound infragravity wave spectrum from a JONSWAP spectrum with the Herbers
    ! et al. (1994) method.
    !
    ! Input:
    !   hsinc      : incident wave height Hm0 [m]
    !   tpinc      : incident peak period [s]
    !   scoeff     : directional spreading coefficient s [-]
    !   jonswapgam : JONSWAP peak enhancement factor [-]
    !   depth      : water depth [m]
    !   correctHm0 : 1 to correct Hm0 of the constructed spectrum, 0 not to
    !
    ! Output:
    !   hsig       : infragravity wave height [m]
    !   Tm01       : mean period Tm01 [s]
    !   Tm10       : mean period Tm-1,0 [s]
    !   Tp         : peak period [s]
    !   Tpsmooth   : smoothed peak period [s]
    !
    ! Called from: determine_ig_bc (snapwave_infragravity)
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
    !
    ! Evaluate the JONSWAP spectral shape for peak enhancement factor gam.
    !
    ! Input:
    !   x   : nondimensional frequency f/fp
    !   gam : peak enhancement factor [-]
    !
    ! Input/output:
    !   y   : JONSWAP spectral shape at x
    !
    ! Called from: compute_herbers (snapwave_infragravity)
    !

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
      !
      ! Find the range of frequencies where the spectrum exceeds the threshold sprdthr.
      !
      ! Input:
      !   Sf          : frequency spectrum
      !   sprdthr     : threshold as a fraction of the spectral peak
      !
      ! Output:
      !   firstp      : index of the first frequency above the threshold
      !   lastp       : index of the last frequency above the threshold
      !   findlineout : 1 where the spectrum exceeds the threshold, else 0
      !
      ! Called from: compute_herbers (snapwave_infragravity)
      !

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
      !
      ! Called from: compute_herbers (snapwave_infragravity)
      !

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
      !
      ! Compute the wave periods Tm01, Tm-1,0, Tp and smoothed Tp from a frequency spectrum.
      !
      ! Input:
      !   Sf       : frequency spectrum
      !   f        : frequencies [Hz]
      !   trepfac  : threshold as a fraction of the spectral peak, used for the smoothed
      !              peak period
      !
      ! Output:
      !   Tm01     : mean period Tm01 [s]
      !   Tm10     : mean period Tm-1,0 [s]
      !   Tp       : peak period [s]
      !   Tpsmooth : smoothed peak period [s]
      !
      ! Called from: compute_herbers (snapwave_infragravity)
      !

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
   !
   !
   subroutine determine_infragravity_source_sink_term(inner, no_nodes, ntheta, w, ds, prev, dtheta, cg_ig, nwav, depth, zb, H, ee, ee_ig, ig_opt, alphaigfac, alphaig_local, beta_local, srcig_local, gamma, gamma_fac_br)
    !
    ! Determining of IG source term as defined in Leijnse et al. 2024
    !
    ! inout: alphaig_local, srcig_local, beta_local
    ! in: the rest
    !
    ! NOTE - This is based on the energy in the previous SnapWave timestep 'ee' and 'ee_ig', and waveheight 'H', which should therefore be made available.
    !
    ! Input:
    !   inner         : mask of the inner nodes (not on a boundary)
    !   no_nodes      : number of grid nodes
    !   ntheta        : number of wave directions
    !   w             : weights of the two upwind nodes, per direction and node
    !   ds            : distance to the upwind point, per direction and node
    !   prev          : indices of the two upwind nodes, per direction and node
    !   dtheta        : directional resolution [rad]
    !   cg_ig         : infragravity group velocity [m/s]
    !   nwav          : ratio of group to phase velocity n [-]
    !   depth         : water depth [m]
    !   zb            : bed level [m]
    !   H             : wave height [m]
    !   ee            : wave energy density, per direction and node
    !   ee_ig         : infragravity energy density, per direction and node
    !   ig_opt        : infragravity wave option
    !   alphaigfac    : factor on the infragravity shoaling source term
    !   gamma         : breaker index gamma [-]
    !   gamma_fac_br  : factor on gamma for the maximum incident wave height
    !
    ! Input/output:
    !   alphaig_local : local infragravity shoaling parameter, per direction and node
    !   beta_local    : local bed slope, per direction and node
    !   srcig_local   : infragravity source/sink term, per direction and node
    !
    ! Called from: solve_energy_balance2Dstat (snapwave_solver)
    !   
    implicit none
    !  
    ! Incoming variables
    !   
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
    real*4, intent(in)                               :: dtheta          ! directional resolution  
    real*4, intent(in)                               :: gamma           ! coefficients in Baldock wave breaking dissipation    
    real*4, intent(in)                               :: gamma_fac_br    ! factor times gamma that is used to determine the maximum incident wave breaking point in the surf zone using local incident wave height over water depth ratio, among others used to set the IG source term to 0 shallower than this point    
    !
    ! Inout variables
    !
    real*4, dimension(:,:), intent(inout)            :: alphaig_local   ! Local infragravity wave shoaling parameter alpha
    real*4, dimension(:,:), intent(inout)            :: srcig_local     ! Energy source/sink term because of IG wave shoaling
    real*4, dimension(ntheta,no_nodes), intent(inout):: beta_local      ! Local bed slope based on bed level per direction
    !
    ! Internal variables
    !
    integer                                          :: itheta          ! directional counter
    integer                                          :: k               ! counters (k is grid index)
    integer                                          :: k1,k2           ! upwind counters (k is grid index)
    real*4                                           :: gam             ! local gamma (Hinc / depth ratio)
    real*4, dimension(ntheta,no_nodes)               :: depthprev       ! water depth at upwind intersection point
    real*4, dimension(no_nodes)                      :: Sxx             ! radiation Stress
    real*4, dimension(ntheta)                        :: Sxxprev         ! radiation stress at upwind point
    real*4, dimension(ntheta)                        :: Hprev           ! wave height at upwind point
    real*4, dimension(ntheta)                        :: cgprev          ! group velocity at upwind point
    real*4, dimension(ntheta)                        :: Eprev           ! Mean incident wave energy at upwind intersection point      
    real*4, dimension(ntheta)                        :: Eprev_ig        ! Mean infragravity wave energy at upwind intersection point      
    real*4, dimension(no_nodes)                      :: E_local         ! mean wave energy waves - just local               
    real*4, dimension(no_nodes)                      :: E_ig_local      ! mean wave energy infragravity waves - just local    
    real*4                                           :: dSxx            ! difference in Radiation stress
    real*4                                           :: Sxx_cons        ! conservative estimate of radiation stress using conservative shoaling
    real*4                                           :: transition_factor ! Transition factor for letting srcig go to zero smoothly, around gamma*gamma_fac_br
    real*4                                           :: transition_factor_width_1 ! Width factor of generalized (FermiDirac style) transfer function with adjustable midpoint and width
    real*4                                           :: transition_factor_width_2 ! Width factor of generalized (FermiDirac style) transfer function with adjustable midpoint and width    
    real*4                                           :: gamma_fac_br_transition ! Transitioned version of gamma_fac_br, so that for steep slopes it remains 1.0
    real*4                                           :: beta_limit_1    ! Cut-off beta_local for end of validity alphaig formulation of Leijnse et al. 2024
    real*4                                           :: beta_limit_2    ! Beta_local limit for transition function    
    !
    ! Set internal variables
    !
    Sxx = 0.0
    Hprev = 0.0
    Eprev = 0.0
    Eprev_ig = 0.0
    !
    E_local    = 0.0
    E_ig_local = 0.0
    !
    ! Used is generalized (FermiDirac style) transfer function with adjustable midpoint and width
    !
    transition_factor_width_1 = 0.005 
    transition_factor_width_2 = 0.002
    beta_limit_1 = 0.07
    !beta_limit_2 = beta_limit_1 - 0.01  
    beta_limit_2 = beta_limit_1 - 0.02    
    !
    ! Pre-compute Sxx for all nodes
    !
    !$omp parallel do schedule(static)
    do k = 1, no_nodes
        !
        if (inner(k)) then !TODO: check whether should be on only 'inner' or not
            !
            ! Update E (not saved from previous timestep)
            !
            E_local(k)         = sum(ee(:,k)) * dtheta                
            !
            ! Update E_ig (not saved from previous timestep)
            !
            E_ig_local(k)      = sum(ee_ig(:, k)) * dtheta   
            !
        endif
        !
        Sxx(k) = ((2.0 * max(0.0, min(1.0, nwav(k)))) - 0.5) * E_local(k)
        !
    enddo
    !$omp end parallel do
    !
    ! Main loop: compute IG source/sink term per node.
    ! All writes target the column (itheta, k), so the loop is data-independent across k.
    ! Per-k scratch arrays (cgprev, Eprev, Eprev_ig, Sxxprev, Hprev) are
    ! listed as private so each thread gets its own copy on the stack.
    !
    !$omp parallel do &
    !$omp&   private(itheta, k1, k2, gam, dSxx, Sxx_cons, &
    !$omp&           gamma_fac_br_transition, transition_factor, &
    !$omp&           cgprev, Eprev, Eprev_ig, Sxxprev, Hprev) &
    !$omp&   schedule(static)
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
                if (k1 > 0 .and. k2 > 0) then ! IMPORTANT - for some reason (k1*k2)>0 is not reliable always, resulting in directions being uncorrectly skipped!!!
                    !
                    ! First calculate upwind direction dependent variables
                    !
                    depthprev(itheta,k)   = w(1, itheta, k) * depth(k1) + w(2, itheta, k) * depth(k2)
                    !
                    beta_local(itheta,k)  = max((w(1, itheta, k) * (zb(k) - zb(k1)) + w(2, itheta, k) * (zb(k) - zb(k2))) / ds(itheta, k), 0.0)
                    !
                    ! Notes:
                    ! - use actual bed level now for slope, because depth changes because of wave setup/tide/surge
                    ! - in zb, depth is negative > therefore zb(k) minus zb(k1)
                    ! - beta=0 means a horizontal or decreasing slope > need alphaig=0 then in IG src/sink term
                    !
                    !betan_local(itheta,k) = (beta/sigm_ig)*sqrt(9.81/max(depth(k), hmin)) ! TL: in case in the future we would need the normalised bed slope again
                    !
                    cgprev(itheta)    = w(1, itheta, k) * cg_ig(k1) + w(2, itheta, k) * cg_ig(k2)
                    !
                    ! Sxx is pre-computed for all nodes before this parallel loop
                    !
                    Sxxprev(itheta)       = w(1, itheta, k) * Sxx(k1) + w(2, itheta, k) * Sxx(k2)
                    !
                    Eprev(itheta)       = w(1, itheta, k) * E_local(k1) + w(2, itheta, k) * E_local(k2)     
                    Eprev_ig(itheta)    = w(1, itheta, k) * E_ig_local(k1) + w(2, itheta, k) * E_ig_local(k2)  
                    !
                    Hprev(itheta)         = w(1, itheta, k) * H(k1) + w(2, itheta, k) * H(k2)
                    !
                    ! Determine relative waterdepth 'gam'
                    !
                    gam = max(0.5 * (Hprev(itheta) / depthprev(itheta,k) + H(k) / depth(k)), 0.0) ! mean gamma over current and upwind point
                    !
                    ! Determine dSxx and IG source/sink term 'srcig'
                    !
                    if (ig_opt == 1 .or. ig_opt == 2) then 
                        !
                        ! Calculate shoaling parameter alpha_ig following Leijnse et al. (2024)
                        !  
                        call estimate_shoaling_parameter_alphaig(beta_local(itheta, k), gam, alphaig_local(itheta, k)) ! [input, input, output]
                        !                
                        ! Now calculate source term component
                        !         
                        ! Newest dSxx/dx based method, using estimate of Sxx(k) using conservative shoaling
                        !
                        if (Sxxprev(itheta) <= 0.0) then 
                           !
                           srcig_local(itheta, k) = 0.0 !Avoid big jumps in dSxx that can happen if a upwind point is a boundary point with Hinc=0
                           !
                        else
                           !              
                           if (ig_opt == 1) then ! Option using conservative shoaling for dSxx/dx
                               !
                               ! Calculate Sxx based on conservative shoaling of upwind point's energy: 
                               ! Sxx_cons = E(i-1) * Cg(i-1) / Cg * (2 * n(i) - 0.5)
                               !
                               Sxx_cons = Eprev(itheta) * cgprev(itheta) / cg_ig(k) * ((2.0 * max(0.0, min(1.0, nwav(k)))) - 0.5)
                               !
                               ! Note - limit so value of nwav is between 0 and 1, and Sxx therefore doesn't become NaN for nwav=Infinite
                               !
                               dSxx = Sxx_cons - Sxxprev(itheta)
                               !
                           elseif (ig_opt == 2) then ! Option taking actual difference for dSxx/dx
                               !
                               dSxx = Sxx(k) - Sxxprev(itheta)
                               !
                           endif
                           !
                           dSxx = max(dSxx, 0.0)
                           !
                           ! Base on E_prev_ig instead of eeprev_ig(itheta) > no bins but total energy
                           ! NOTE - already here multiplied with ee(itheta,k), for direct inclusion in 'R'-term
                           srcig_local(itheta, k) = alphaigfac * alphaig_local(itheta,k) * sqrt(Eprev_ig(itheta)) * cgprev(itheta) / depthprev(itheta,k) * dSxx / ds(itheta, k) /max(E_local(k), 1.0e-6) * ee(itheta,k)                           
                           !
                           ! Limit srcig to 0 after waves start (significantly) breaking, as defined here as gam=Hrms,inc / h > (gamma_fac_br * gamma)
                           ! Ergo, it is assumed that after this point IG waves are free, and no bound wave forcing is happening anymore, so srcig should be 0 from here on
                           !                           
                           ! Let srcig transition to 0 more smoothly using fac_transition that reduced from 1 to 0 around gamma_fac_br * snapwave_gamma
                           ! But, only for beta_local < 0.07, so adjust based on beta_local so that transition_factor = 1.0 for Beta_local = 0.07
                           !
                           gamma_fac_br_transition = gamma_fac_br + ((1-gamma_fac_br) / (1 + exp(- (beta_local(itheta,k) - beta_limit_2) / transition_factor_width_2)))
                           !
                           transition_factor = 1.0 - (1.0 / (1.0 + exp(- (gam - (gamma_fac_br_transition * gamma)) / transition_factor_width_1)))
                           !
                           srcig_local(itheta, k) = transition_factor * srcig_local(itheta, k)                           
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
    !$omp end parallel do
    !
   end subroutine determine_infragravity_source_sink_term   
   !
   !
   subroutine estimate_shoaling_parameter_alphaig(beta, gam, alphaig)
   !
   ! Estimate the infragravity shoaling parameter alphaig from the bed slope and relative
   ! depth - as in Leijnse et al. (2024).
   !
   ! Input:
   !   beta    : local bed slope [-]
   !   gam     : ratio of wave height to water depth [-]
   !
   ! Output:
   !   alphaig : infragravity shoaling parameter
   !
   ! Called from: determine_infragravity_source_sink_term (snapwave_infragravity)
   !
   implicit none
   !
   real*4, intent(in)                :: beta
   real*4, intent(in)                :: gam 
   real*4, intent(out)               :: alphaig
   !
   real*4                            :: beta1, beta2, beta3, beta4, beta5, beta6, beta7        
   !
   ! Estimate shoaling parameter alphaig
   !
   ! Determine constants 
   !
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
   !

               
end module    