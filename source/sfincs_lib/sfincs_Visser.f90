   module sfincs_Visser
        ! not added to github, not sure why
       use sfincs_error
       use sfincs_log
       type :: NormalFlow 
            real*4 :: n    ! velocity exponent
            real*4 :: stt  ! sediment total transport capacity
            real*4 :: normal_water_depth
            real*4 :: normal_flow_velocity
            real*4 :: friction_coeff
            real*4 :: water_depth
            real*4 :: slope_flow_velocity
            real*4 :: discharge
            real*4 :: t1
            real*4 :: t2
            real*4 :: breach_width_waterline
            real*4 :: breach_width_total
            real*4 :: gamma0
            real*4 :: breach_bottom
            real*4 :: breach_level
       end type NormalFlow
       
       contains

             
       function critical_shields_parameter(dimless_diameter, formula) result(theta_crit)
   
       implicit none
       real*4 :: dimless_diameter
       character*256 :: formula
       real*4 :: theta_crit


       select case (formula)
       case ('Shields')
            theta_crit = 0.24 / dimless_diameter + 0.055 * (1.0 - exp(-0.020 * dimless_diameter))
            !write(logstr,'(a,G12.6,a,G12.6)') 'theta_crit', theta_crit, ' dimless_diameter:',dimless_diameter
            !call write_log(logstr,1)
       case ('Soulsby')
            theta_crit = 0.30 / (1.0 + 1.2 * dimless_diameter) + 0.055 * (1.0 - exp(-0.020 * dimless_diameter))

       case default
            ! Mirror Python's ValueError behavior
            write(logstr, '(a)') "Invalid value for formula: '"//trim(formula)// "'! Allowed options are 'Shields' and 'Soulsby'."
            error stop logstr
       end select
  
       end function
       
       function dimensionless_diameter (d50, rel_density, dyn_viscosity) result (dstar)
       implicit none
       real*4:: d50, rel_density, dyn_viscosity, g
       real*4:: d50_min, dstar
   
       g = 9.81  ! m/s2, gravitational acceleration
   
       d50_min = 1e-6

       if (d50 > d50_min) then
            dstar = d50 * (rel_density * g / dyn_viscosity ** 2.0) ** (1.0 / 3.0)
            !write(logstr,'(a,G12.6,a,G12.6,a,G12.6,a,G12.6)') 'dstar: ', dstar, ' rel_density:',rel_density, ' d50:',d50, ' dyn_viscosity:',dyn_viscosity
            !call write_log(logstr,1)
       else
            write(logstr,'(a,f12.6,a,f12.6, a)') "The value for d50: ",d50," is out of range! d50 should be larger than ", d50_min*1e6, " micrometre."
            error stop logstr
       end if
       end function

   
       function breach_crest_length(crest_width, crest_level, breach_level, slope_angle_outer, slope_angle_inner) result (W)
       implicit none
       REAL*4:: crest_width, crest_level, breach_level, slope_angle_outer, slope_angle_inner
       REAL*4:: W, pi
   
       pi = 3.141592653589793
   
       if (crest_width < 0) then
            write(logstr,'(a,f6.4,a)') "Crest width: ",crest_width," is negative!"
            error stop logstr
       end if
   
       if (crest_level <= breach_level) then
           write(logstr,'(a,f6.4,a,f6.4)') "Crest level: ", crest_level," is below breach level: ",breach_level
           error stop logstr
       end if
   
       if (slope_angle_outer > pi/2 .or. slope_angle_outer <= 0) then
            write(logstr,'(a,f6.4,a)') "Dike outer slope angle: ",slope_angle_outer," is not within range! The range is: 0 < radians <= pi/2"
            error stop logstr
       end if
   
       if (slope_angle_inner > pi/2 .or. slope_angle_inner <= 0) then
            write(logstr,'(a,f6.4,a)') "Dike inner slope angle: ",slope_angle_inner," is not within range! The range is: 0 < radians <= pi/2"
            error stop logstr
       end if
    
       W = (crest_width + (crest_level - breach_level) * (1 / tan(slope_angle_outer) + 1 / tan(slope_angle_inner)))
       end function
   
   
       function dynamic_viscosity(water_temperature, water_density) result (nu)
       !Computes the dynamic viscosity of water in Ns/m^2.
       implicit none
       REAL*4:: water_temperature, water_density
       REAL*4:: nu
       REAL*4:: temp_min, temp_max

       temp_min = 0
       temp_max = 40

       if ((0 <= water_temperature) .AND. (water_temperature <= 20)) then
            nu = ((water_density + 1505) / (2500 * water_density) * 10**(13 / (10 - 0.081 * (20 - water_temperature)) - 4.3))

       elseif ((20 < water_temperature) .AND. (water_temperature<= 40)) then
            nu = ((water_density + 1505) / (2500 * water_density) * 10**((1.33 * (20 - water_temperature) / (water_temperature + 104)) - 3.0))
       else
            error stop "Water temperature is out of range! The valid range is 0 to 40 degrees Celsius"
       end if 
       end function
   
       function calc_sediment_fall_velocity(d50, rel_density, dyn_viscosity) result (ws)
       !Computes the sediment fall velocity in m/s.
       implicit none
       REAL*4:: d50, rel_density, dyn_viscosity
       REAL*4:: ws
       REAL*4:: g, d50_min
       g = 9.81  ! m/s2, gravitational acceleration
       d50_min = 1e-6
       !write(logstr,'(a,f6.3,a,f6.3)')'d50_min: ', d50_min, ' rel_density: ', rel_density
       !call write_log(logstr,1)
       !write(logstr,'(a,G12.3,a,f6.3)')'dyn_viscosity: ', dyn_viscosity, ' rel_density: ', rel_density
       !call write_log(logstr,1)

       if ((1e-6 <= d50) .AND. (d50<= 1e-4)) then
            ws = rel_density * g * d50**2 / (18 * dyn_viscosity)

       elseif ((1e-4 < d50) .AND. (d50 < 1e-3)) then
            ws = ((10 * dyn_viscosity) / d50* (((1 + 0.01 * rel_density * g * d50**3 * dyn_viscosity**(-2))**0.5) - 1))

       elseif (d50 >= 1e-3) then
            ws = 1.1 * (rel_density * g * d50)**0.5
       else
            error stop "The value for d50 is out of range!"
       end if
       !write(logstr,'(a,G12.3)')'ws: ', ws
       !call write_log(logstr,1)
    
       end function
       
       function calc_crit_water_depth(water_level, breach_level, breach_width_bottom, slope_angle_breach) result (crit_water_depth)
       implicit none
       real*4:: water_level, breach_level, breach_width_bottom, slope_angle_breach
       real*4:: crit_water_depth, crit_water_depth_iter, breach_width_avg_water_depth, breach_width_waterline
       real*4, parameter :: pi =  3.141592653589793  ! pi in single precision
       
       !! Not sure to add this below because there might be cases where the water level is below the breach level due to the tides
       !if (water_level <= breach_level) then
       !     error stop "Water level should be above breach level!"
       !end if
   
       ! Estimate critical water depth from frictionless flow (2/3=0.7 of water depth)
       crit_water_depth = 0.7 * (water_level - breach_level)
   
       !write(logstr,'(a,f6.4,a,f6.4,a,f6.4)') 'START crit_water_depth: ', crit_water_depth, ' water_level: ', water_level, ' breach_level:', breach_level
       !call write_log(logstr,1)

       ! Define an iteration variable to start while loop
       crit_water_depth_iter = 0.9 * crit_water_depth
       !write(logstr,'(a,f6.4,a,f6.3,a,f6.3)') 'START crit_water_depth_iter: ', crit_water_depth_iter, ' crit_water_depth: ', crit_water_depth
       !call write_log(logstr,1)
   
   
       do while (abs((crit_water_depth - crit_water_depth_iter) / crit_water_depth) > 0.01)
            ! Set iteration variable to previous estimate
            crit_water_depth_iter = crit_water_depth
    
            if (breach_width_bottom <= 0) then
                write(logstr, '(a,f12.4,a,f12.4,a,f12.4)') "Breach bottom width: ", breach_width_bottom, " is non-positive! Water_level: ", water_level," breach_level: ", breach_level
                error stop logstr
                error stop ""
            end if

            if (crit_water_depth < 0) then
                !error stop "Water depth ", crit_water_depth, " is negative! because outside water level is ", outside_water_level, " and breach level is ", breach_level
                write(logstr,'(a,f6.4,a,f6.4,a,f6.4, a)')"Warning: Critical water depth is negative (",crit_water_depth,") because outside water level (", water_level ,") is below breach level (",breach_level,"). Setting critical water depth to zero."
                error stop logstr
            end if

            if (slope_angle_breach > pi/2 .or. slope_angle_breach <= 0) then
                write(logstr, '(a,f6.4,a)') "In function calc_crit_water_depth, Breach slope angle ", slope_angle_breach," is not within range, the range is: 0 < radians <= pi/2"
                error stop logstr
            end if

            ! Calculate widths
            breach_width_avg_water_depth = breach_width_bottom + (crit_water_depth / tan(slope_angle_breach))
            breach_width_waterline = breach_width_bottom + (2 * crit_water_depth / tan(slope_angle_breach))

            ! Update estimate of critical water depth
            crit_water_depth = (2.0 / (2.0 + breach_width_avg_water_depth / breach_width_waterline)) * &
                             (water_level - breach_level)
            !write(logstr,'(a,f6.4,a,f6.4,a,f6.4,a,f6.4)') 'crit_water_depth: ', crit_water_depth, ' breach_level: ', breach_level, ' breach_width_avg_water_depth:', breach_width_avg_water_depth, ' breach_width_waterline:', breach_width_waterline
            !call write_log(logstr,1)
       end do
       end function 
   
       function calc_normal_flow_velocity(hydraulic_radius, slope_angle_dike, friction_coeff) result (normal_flow_velocity)
       ! Calculates the normal flow velocity on the dike inner slope
       implicit none
       real*4:: hydraulic_radius, slope_angle_dike, friction_coeff, normal_flow_velocity
       real*4 :: g
       g = 9.81
       normal_flow_velocity = (g * hydraulic_radius * sin(slope_angle_dike) / friction_coeff)**0.5
       end function
       
   
       function calc_crit_flow_velocity(crit_water_depth, breach_width_avg_water_depth, breach_width_waterline) result (crit_flow_velocity)
       !Calculates the critical flow velocity in the breach based on broad-crested weir flow
       implicit none
       real*4:: crit_water_depth, breach_width_avg_water_depth, breach_width_waterline, crit_flow_velocity
       real*4:: g
       g = 9.81
       crit_flow_velocity = (g * crit_water_depth * breach_width_avg_water_depth / breach_width_waterline)**0.5
       end function
   
   
       function calc_discharge(breach_width_avg_water_depth, water_depth, flow_velocity, discharge_coeff) result (discharge)
       !Calculates the discharge through the breach in stages I to V.
       !In stages I to IV, assumes critical flow over a weir.
       !In stage V, assumes free gravity flow.
       implicit none
       real*4:: breach_width_avg_water_depth, water_depth, flow_velocity, discharge_coeff, discharge
     
       discharge = discharge_coeff * breach_width_avg_water_depth * water_depth * flow_velocity
       end function
       
       
       function normal_flow_conditions(water_level, breach_level, breach_width_bottom, discharge, slope_angle_breach, slope_angle_dike, delta, d50, d90, friction_coeff, kappa) result (res)
       ! Calculates the normal flow conditions on the slope
       implicit none
       real*4:: water_level, breach_level, breach_width_bottom, discharge, slope_angle_breach, slope_angle_dike, delta, d50, d90, friction_coeff, kappa,g,normal_water_depth, normal_flow_velocity
       REAL*4:: breach_width_avg_water_depth, normal_water_depth_iter, hydraulic_radius, theta
   
    
       !type :: NormalFlow
       !  real*4 :: normal_water_depth
       !  real*4 :: normal_flow_velocity
       !  real*4 :: friction_coeff
       !end type NormalFlow

       ! Result
       type(NormalFlow) :: res


       g = 9.81
       ! Estimate the normal water depth
       normal_water_depth = 0.3 * (water_level - breach_level)
       !write(logstr,'(a,G12.6,a,G12.6)') 'Start** water_level:', water_level, ' breach_level:', breach_level
       !call write_log(logstr,1)

       breach_width_avg_water_depth = breach_width_bottom + (normal_water_depth / tan(slope_angle_breach))
       normal_water_depth = ((friction_coeff * (discharge / breach_width_avg_water_depth)**2) / (g * sin(slope_angle_dike)))**(1.0/3.0)
       !write(logstr,'(a,G12.6,a,G12.6)') 'Start** breach_width_avg_water_depth:', breach_width_avg_water_depth, ' normal_water_depth:', normal_water_depth
       !call write_log(logstr,1)


       ! Define an iteration variable to start the while loop
       normal_water_depth_iter = 0.9 * normal_water_depth

       do while (abs((normal_water_depth - normal_water_depth_iter) / normal_water_depth) > 0.01)
            normal_water_depth_iter = normal_water_depth
            breach_width_avg_water_depth = breach_width_bottom + (normal_water_depth / tan(slope_angle_breach)) 
            hydraulic_radius = calc_hydraulic_radius(normal_water_depth, breach_width_bottom, breach_width_avg_water_depth, slope_angle_breach)
            normal_flow_velocity = calc_normal_flow_velocity(hydraulic_radius, slope_angle_dike, friction_coeff)
            normal_water_depth = discharge / (normal_flow_velocity * breach_width_avg_water_depth)
            theta = calc_shields_parameter(normal_flow_velocity, delta, d50, friction_coeff)
            friction_coeff = calc_friction_coefficient(theta, hydraulic_radius, d90, kappa)
       end do

   
       res%normal_water_depth = normal_water_depth
       res%normal_flow_velocity = normal_flow_velocity
       res%friction_coeff = friction_coeff
   
       end function
   
       function calc_shields_parameter(flow_velocity, delta, d50, friction_coeff) result (theta)
       !Calculates the Shields mobility parameter.
       implicit none
       REAL*4:: flow_velocity, delta, d50, friction_coeff, g, theta
       g = 9.81
   
       !write(logstr,'(a,G12.6,a, G12.6)') 'friction_coeff:', friction_coeff, ' flow_velocity:', flow_velocity
       !call write_log(logstr,1)
       !write(logstr,'(a,G12.6,a,f6.4)') ' d50:', d50, ' delta:', delta
       !call write_log(logstr,1)
   
       if (friction_coeff <= 0) then
            error stop "Friction coefficient is non-positive!"
       end if

       if (flow_velocity <= 0)  then
            error stop "Flow velocity is non-positive!"
       end if 

       theta = (friction_coeff * flow_velocity**2) / (g * delta * d50)
       end function
   
       function calc_friction_coefficient(theta, hydraulic_radius, d90, kappa) result (friction_coeff)
       !Computes the friction coefficient of the sediment layer and corresponding Shields (mobility) parameter.
       implicit none
       REAL*4:: theta, hydraulic_radius, d90, kappa, friction_coeff, k
   
       if (theta < 0) then
            error stop "Shields parameter (theta) is negative!"
       end if

       if (hydraulic_radius <= 0) then
             error stop "Hydraulic radius is non-positive!"
       end if

       ! A limit for theta to avoid very high k
       theta = min(theta, 300.0)

       ! Compute layer roughness (k)
       if (theta < 1) then
            k = 3.0 * d90
       else
           k = 3.0 * d90 * theta
       end if

       ! New estimate of the friction coefficient
       friction_coeff = kappa**2 / ((log(12 * hydraulic_radius / k))**2.0)
   
       !write(logstr,'(a,G12.6,a,f6.4,a, G12.6)') 'friction_coeff:', friction_coeff, ' theta:', theta, ' hydraulic_radius:', hydraulic_radius
       !call write_log(logstr,1)
       !write(logstr,'(a,G12.6,a,f6.4)') ' d90:', d90, ' kappa:', kappa
       !call write_log(logstr,1)
   
       end function
   
       function calc_hydraulic_radius(water_depth, breach_width_bottom, breach_width_avg_water_depth, slope_angle_breach) result (hydraulic_radius)

       !Computes the hydraulic radius of the trapezoidal cross-section.
       implicit none
       REAL*4:: water_depth, breach_width_bottom, breach_width_avg_water_depth, slope_angle_breach, hydraulic_radius
       REAL*4:: area, wetted_perimeter
       real*4, parameter :: pi = 3.141592653589793  ! pi in single precision
   
       if (water_depth < 0.0) then
            write(logstr, '(a,f6.4)') "Water depth is negative in calc_hydraulic_radius: ", water_depth
            error stop logstr
       end if

       if (breach_width_bottom < 0) then
            error stop "Breach bottom width is negative!"
       end if

       if (slope_angle_breach > pi/2 .or. slope_angle_breach <= 0) then
            write(logstr, '(a,f6.4,a)') "In function calc_hydraulic_radius, Breach slope angle ", slope_angle_breach," is not within range, the range is: 0 < radians <= pi/2"
            error stop logstr
       end if 

       area = breach_width_avg_water_depth * water_depth
       wetted_perimeter = breach_width_bottom + (2 * water_depth / sin(slope_angle_breach))
       hydraulic_radius = area / wetted_perimeter
       end function
   
       function calc_breach_width_avg_water_depth(breach_width_bottom, water_depth, slope_angle_breach) result (breach_width_avg_water_depth)
       !Calculates the breach width at average water depth
       implicit none
       REAL*4:: breach_width_bottom, water_depth, slope_angle_breach, breach_width_avg_water_depth
       real*4, parameter :: pi =  3.141592653589793 ! pi in single precision
       if (breach_width_bottom <= 0) then
            error stop "Breach bottom width is non-positive!"
       end if

       if (water_depth < 0) then
            write(logstr, '(a,G12.4)') "Water depth is negative in calc_breach_width_avg_water_depth: ", water_depth
            error stop logstr
       end if

       if (slope_angle_breach > pi/2 .or. slope_angle_breach <= 0) then
           write(logstr, '(a,f6.4,a)') "In function calc_breach_width_avg_water_depth, Breach slope angle ", slope_angle_breach," is not within range, the range is: 0 < radians <= pi/2"
           error stop logstr
       end if

       breach_width_avg_water_depth = breach_width_bottom + (water_depth / tan(slope_angle_breach))
       end function
   
       function calc_breach_width_waterline(breach_width_bottom, water_depth, slope_angle_breach) result(breach_width_waterline)
       ! Calculates the breach width at the waterline
       implicit none
       REAL*4:: breach_width_bottom, water_depth, slope_angle_breach, breach_width_waterline
       real*4, parameter :: pi =  3.141592653589793  ! pi in single precision

        if (breach_width_bottom <= 0) then
             error stop "Breach bottom width is non-positive!"
       end if

       if (water_depth < 0) then
            write(logstr, '(a,f6.4)') "Water depth is negative in calc_breach_width_waterline: ", water_depth
            error stop logstr
       end if

       if (slope_angle_breach > pi/2 .or. slope_angle_breach <= 0) then
           write(logstr, '(a,f6.4,a)') "In function calc_breach_width_waterline, Breach slope angle ", slope_angle_breach," is not within range, the range is: 0 < radians <= pi/2"
           error stop logstr
       end if

       breach_width_waterline = breach_width_bottom + (2 * water_depth / tan(slope_angle_breach))
       end function
   
       function calc_discharge_coeff(formula, outer_slope_angle, inner_slope_angle, outside_water_level, breach_level, weir_length) result(discharge_coeff)
   
       implicit none
       character(len=* ) :: formula
       real*4 :: outer_slope_angle, inner_slope_angle, outside_water_level, breach_level, weir_length
       real*4 :: outer_slope, inner_slope, c1, c2, Cd, discharge_coeff, zeta, arg
       character(len=200) :: errmsg

       select case (trim(formula))
       case ('Chen2018')
            outer_slope = 1.0 / tan(outer_slope_angle)
            inner_slope = 1.0 / tan(inner_slope_angle)

            if (inner_slope >= 0.0 .and. inner_slope <= 0.8) then
                c1 = (-1.3 * outer_slope + 8.09) * ( 2.862 * inner_slope +  7.658) * 1.0e-3
                c2 = (-8.6 * outer_slope**2 + 7.9 * outer_slope + 493.5) * ( 5.33 * inner_slope + 95.74) * 1.0e-5
            else
                c1 = (-1.3 * outer_slope + 8.09) * (-1.797 * inner_slope + 11.355) * 1.0e-3
                c2 = (-8.6 * outer_slope**2 + 7.9 * outer_slope + 493.5) * (-4.24 * inner_slope + 103.28) * 1.0e-5
            end if

            arg = (outside_water_level - breach_level) / (breach_level + weir_length)
            if (arg <= 0.0) then
                errmsg = 'Chen2018: log argument <= 0'
                error stop trim(errmsg)
            end if

            Cd = c1 * log(arg) + c2
            discharge_coeff = Cd * sqrt(2.0)

       case ('Zerihun2020')
            zeta = (outside_water_level - breach_level) / weir_length
            Cd = (0.40- 0.215 * sin(outer_slope_angle)**(22.0/125.0)+ 0.13  * sin(inner_slope_angle)**( 3.0/ 20.0) + 0.134 * zeta / (1.0 + 0.596 * zeta) )
            discharge_coeff = Cd * sqrt(2.0)

       case default
            errmsg = "Invalid value 'formula'. Allowed: 'Chen2018', 'Zerihun2020'."
            error stop trim(errmsg)
       end select

       end function
    
       function calc_water_depth_width_avg(water_depth, breach_width_avg_water_depth, breach_width_waterline) result(water_depth_width_avg)
       !Calculates the width-averaged water depth in the breach.
       implicit none
       REAL*4:: water_depth, breach_width_avg_water_depth, breach_width_waterline, water_depth_width_avg
    

       if (water_depth < 0) then
            write(logstr, '(a,f6.4)') "Water depth is negative in calc_water_depth_width_avg: ", water_depth
            error stop logstr
       end if

       if (breach_width_avg_water_depth > breach_width_waterline) then
            error stop "Water depth averaged breach width by definition cannot be larger than breach width at waterline"
       end if

       water_depth_width_avg = water_depth * (breach_width_avg_water_depth / breach_width_waterline)    
       end function
    
       function calc_froude_number(flow_velocity, water_depth_width_avg, slope_angle_dike) result (froude_number)
       !Computes the Froude number of the flow.
       implicit none
       REAL*4:: flow_velocity, water_depth_width_avg, slope_angle_dike, froude_number, g
       real*4, parameter :: pi =  3.141592653589793 ! pi in single precision
    
       g = 9.81

       if (flow_velocity < 0) then
            error stop "Flow velocity is negative!"
       end if

       if (water_depth_width_avg <= 0) then
            error stop "Water depth is non-positive!"
       end if

       if (slope_angle_dike > pi/2 .or. slope_angle_dike < 0) then
            error stop "Slope angle is not within range, the range is: 0 <= radians <= pi/2"
       end if

       froude_number = flow_velocity / (g * water_depth_width_avg * cos(slope_angle_dike))**0.5
    
       end function
    
       function calc_adaptation_length_flow(froude_number, water_depth, slope_angle_dike) result (adaptation_length)
       !Computes the length for the flow to achieve normal flow conditions.
       implicit none
       REAL*4:: froude_number, water_depth, slope_angle_dike, adaptation_length
       real*4, parameter :: pi = 3.141592653589793  ! pi in single precision
       if (froude_number < 1.0) then
           error stop "Froude number is less than 1.0!"
       end if

       if (water_depth <= 0.0) then
            error stop "Water depth is non-positive!"
       end if

       if (slope_angle_dike > pi/2 .or. slope_angle_dike < 0) then
            error stop "Slope angle is not within range, the range is: 0 <= radians <= pi/2"
       end if

       adaptation_length = 2.5 * (froude_number**2 - 1) * water_depth / tan(slope_angle_dike)
    
       end function
    
       function calc_adaptation_length_sediment(discharge, breach_width_avg_water_depth, sediment_fall_velocity, slope_angle_dike) result (adaptation_length)
       !Computes the length for the flow to achieve equilibrium sediment transport capacity.

       implicit none
       REAL*4:: discharge, breach_width_avg_water_depth, sediment_fall_velocity, slope_angle_dike, epsilon, adaptation_length
       real*4, parameter :: pi = 3.141592653589793  ! pi in single precision
    
       epsilon=1.0
    
       if (discharge < 0) then
            error stop "Discharge is negative!"
       end if

       if (slope_angle_dike > pi/2 .or. slope_angle_dike < 0) then
            error stop "Slope angle is not within range, the range is: 0 <= radians <= pi/2"
       end if

       adaptation_length = (epsilon * discharge / (breach_width_avg_water_depth * sediment_fall_velocity * cos(slope_angle_dike)))
       !write(logstr,'(a,G12.6,a,G12.6,a,G12.6)') 'Calc adaptation length:', adaptation_length, ' sediment_fall_velocity:',sediment_fall_velocity, 'slope_angle_dike', slope_angle_dike
       !call write_log(logstr,1)
       end function
    
       function flow_along_slope(discharge, crit_water_depth, normal_water_depth, breach_width_bottom, slope_angle_breach, adaptation_length_flow, loc_along_slope) result (res)
       
       !Computes the flow conditions at a location along the slope.
       implicit none
       REAL*4:: discharge, crit_water_depth, normal_water_depth, breach_width_bottom, slope_angle_breach, adaptation_length_flow, loc_along_slope
       REAL*4:: exponential, water_depth, breach_width_avg_water_depth, slope_flow_velocity
       !type :: NormalFlow
        !real*4 :: water_depth
        !real*4 :: slope_flow_velocity
       !end type NormalFlow

       ! Result
       type(NormalFlow) :: res
   
       ! Calculate the water depth on the slope at the location along the slope
       exponential = exp(-5 * loc_along_slope / adaptation_length_flow)
       water_depth = normal_water_depth + (crit_water_depth - normal_water_depth) * exponential

       ! Calculate the water depth averaged breach width
       breach_width_avg_water_depth = breach_width_bottom + (water_depth / tan(slope_angle_breach))

       ! Calculate the flow velocity on the slope at the location along the slope
       slope_flow_velocity = discharge / (breach_width_avg_water_depth * water_depth)
    
       res%water_depth = water_depth
       res%slope_flow_velocity = slope_flow_velocity

       end function
    
       function iter_friction_coefficient(friction_coeff, flow_velocity, hydraulic_radius, delta, d50, d90, kappa) result (friction_coeff_final)
       !Iteration function to computes the friction coefficient.
       implicit none
       REAL*4:: friction_coeff, flow_velocity, hydraulic_radius, delta, d50, d90, kappa, friction_coeff_iter, theta, g,friction_coeff_final
    
       g = 9.81
       ! Define an iteration variable to start the while loop
       friction_coeff_iter = 0.9 * friction_coeff

       do while (abs((friction_coeff - friction_coeff_iter) / friction_coeff) > 0.01)
            ! Set iteration variable to estimate to use for next iteration
            friction_coeff_iter = friction_coeff
            theta = calc_shields_parameter(flow_velocity, delta, d50, friction_coeff)
            friction_coeff = calc_friction_coefficient(theta, hydraulic_radius, d90, kappa)
       end do
       friction_coeff_final = friction_coeff
       end function
   
       function calc_breach_width_total(breach_width_bottom, crest_level, breach_level, slope_angle_breach) result (breach_width_total)
       !Calculates the breach width at the dike crest
       implicit none
       REAL*4:: breach_width_bottom, crest_level, breach_level, slope_angle_breach, breach_width_total
       real*4, parameter :: pi = 3.141592653589793  ! pi in single precision

       if (crest_level <= breach_level) then
            error stop "Crest level is below breach level!"
       end if

       if (breach_width_bottom < 0) then
            error stop "Breach bottom width is negative!"
       end if

       if (slope_angle_breach > pi/2 .or. slope_angle_breach <= 0) then
           write(logstr, '(a,f6.4,a)') "In function calc_breach_width_total, Breach slope angle ", slope_angle_breach," is not within range, the range is: 0 < radians <= pi/2"
           error stop logstr
       end if

       breach_width_total = breach_width_bottom + 2 * ((crest_level - breach_level) / tan(slope_angle_breach))
       end function
   
       function submerged_weir_flow (outside_water_level, breach_level, polder_water_level, discharge_coeff) result(discharge_coeff_reduced)
       implicit none
       REAL*4:: outside_water_level, breach_level, polder_water_level, discharge_coeff, discharge_coeff_reduced
       REAL*4:: submergence_ratio, discharge_reduction
   
       submergence_ratio = (polder_water_level - breach_level) / (outside_water_level - breach_level)
        ! The equation for discharge_reduction is originally valid for submergence_ratio >= 0.84. The result should go from
        ! 1.0 to 0.0. Here we apply it to all submergence ratios from 0 to 1. Up to submergence_ratio 0.34 the result grows
        ! from 0 to 1, hence this falls outside the applicable range. We only use the equation above submergence_ratio 0.34,
        ! and the result maximum must be 1.0.
       if (submergence_ratio <= 0.35) then
            discharge_coeff_reduced = discharge_coeff
       else
           discharge_reduction = ((2.41 * (outside_water_level - breach_level)**0.03 * (1.0 - submergence_ratio)**1.53) / ((-log(submergence_ratio))**1.20 * 9.81**0.5))
           discharge_reduction = min(1.0, discharge_reduction)
           discharge_coeff_reduced = discharge_coeff * discharge_reduction
       end if
       end function
    
       !!!!! Sediment function
    
       function bagnold_visser(delta, p, d50, phi, ws, flow_velocity, slope_angle, friction_coeff) result (res)
       !Computes the sediment total transport capacity following the Bagnold-Visser (1989) formulation.
       !Formulation limits:
       !11 < theta < 106
       !2.8 < froude < 4.1
       !100 <= d50 <= 220 micron
       !0.36 < tan(slope_angle) < 0.62
       !1.2 < flow velocity < 3.5 m/s
       !0.007 < concentration < 0.28
    
       implicit none
       REAL*4:: delta, p, d50, phi, flow_velocity, slope_angle, friction_coeff, ws
       REAL*4:: sb_max, sb, ss, stt, g, n
       !type :: NormalFlow
       !  real*4 :: n
        ! real*4 :: stt
       !end type NormalFlow

       ! Result
       type(NormalFlow) :: res

       ! velocity exponent
       n = 4.0
       g = 9.81

       ! upper limit bed load transport
       sb_max = 2.0 * (1.0 - p) * d50 * flow_velocity

       if (slope_angle >= phi) then
            sb = sb_max
       else
            ! sediment bed load transport
            sb = (0.13 / ((tan(phi) - tan(slope_angle)) * cos(slope_angle)) * friction_coeff * flow_velocity**(n-1.0) / (delta * g))

            ! used sediment bed load transport
            sb = min(sb, sb_max)
       end if
    
       ! sediment suspended load transport
       ss = 0.01 * friction_coeff * flow_velocity**n / (delta * g * ws * cos(slope_angle)**2.0)
       !write(logstr,'(a,G12.6,a,G12.6,a,G12.6)') 'phi:', phi, ' slope_angle:',slope_angle, 'ws', ws
       !call write_log(logstr,1)
   
       !write(logstr,'(a,G12.6,a,G12.6)') 'ss:', ss, ' sb:',sb
       !call write_log(logstr,1)

       ! sediment total transport
       stt = sb + ss
    
       res%n = n
       res%stt = stt
       end function   
   
       function vanrhee_simp (theta, theta_crit, adap_length_sediment, water_depth, ni, p, delta, d50, dstar, k, rhos, rhow) result (res)
       implicit none
       real *4:: theta, theta_crit, adap_length_sediment, water_depth, ni, p, delta, d50, dstar, k, rhos, rhow
       real *4:: rel_theta, dilatancy_factor, alpha, erosion_velocity, concentration, stt, n, g
       type(NormalFlow) :: res
   
       g = 9.81
       ! Compute relative Shields parameter
       rel_theta = (theta - theta_crit) / theta_crit
       !write(logstr,'(a,G12.6,a,G12.6,a,G12.6)') 'rel_theta:', rel_theta, ' theta:',theta, ' theta_crit:',theta_crit
       !call write_log(logstr,1)
       ! Negative relative theta means no erosion, thus no transport
       if (rel_theta <= 0) then
            stt = 0
            n = 1.0
       else
            ! Compute dilatancy factor
            dilatancy_factor = (ni - p) / (1.0 - ni) * 1.0 / (delta * (1.0 - p))
            !write(logstr,'(a,G12.6,a,G12.6)') 'dilatancy_factor:', dilatancy_factor, ' ni:',ni
            !call write_log(logstr,1)

            ! Compute alpha factor, an adaptation from the Van Rijn formulation
            alpha = 0.00033 * (delta * g * d50)**0.5 / (1 - p)
            !write(logstr,'(a,G12.6,a,G12.6)') 'alpha:', alpha, ' delta:',delta
            !call write_log(logstr,1)
            ! Compute bed erosion velocity
            erosion_velocity = (alpha**2.0 * dstar**0.6 * rel_theta**3.0 * (k / dilatancy_factor)**3.0)**(1.0/5.0)
            !write(logstr,'(a,G18.6,a,G12.6)') 'k:', k, ' alpha:',alpha
            !call write_log(logstr,1)
            !write(logstr,'(a,G12.6,a,G12.6)') 'erosion_velocity:', erosion_velocity, ' dstar:',dstar
            !call write_log(logstr,1)
            ! Limit the erosion velocity so concentration in water column does not exceed sheared porosity
            concentration = erosion_velocity * rhos * (1 - p) / (water_depth * rhow)
            !write(logstr,'(a,G12.6,a,G12.6)') 'concentration:', concentration, ' water_depth:',water_depth
            !call write_log(logstr,1)
            if (concentration > (1 - ni)) then
                erosion_velocity = erosion_velocity*( (1 - ni) / concentration)
            end if
        
            ! Compute sediment total transport capacity
            stt = erosion_velocity * adap_length_sediment
            !write(logstr,'(a,G12.6,a,G12.6)') 'stt:', stt, ' adap_length_sediment:',adap_length_sediment
            !call write_log(logstr,1)
            write(logstr,'(a,G12.6,a,G18.6,a,G12.6,a,G12.6)') 'stt:', stt,' adap_length_sediment:', adap_length_sediment,'concentration:', concentration, ' erosion_velocity:',erosion_velocity
            call write_log(logstr,0)
            n=1.0
       end if
       res%n = n
       res%stt = stt
       end function
   
       !!!! Stages
   
       function stage_1(t0, breach_bottom, polder_level, polder_water_level, breach_level, beta1, beta0, outside_water_level, gamma0, alpha, W, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity) result(results_t1)
       implicit none
       REAL*4:: t0, polder_level, polder_water_level, breach_bottom, breach_level, beta1, beta0
       REAL*4:: outside_water_level, gamma0
       REAL*4:: alpha, W, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity
       REAL*4:: beta, crit_water_depth
       REAL*4:: crit_breach_width_avg_water_depth, crit_breach_width_waterline
       REAL*4:: crit_flow_velocity, discharge_coeff, discharge
       REAL*4:: breach_width_avg_water_depth, breach_width_waterline, breach_width_total
       REAL*4:: normal_water_depth, normal_flow_velocity, friction_coeff
       REAL*4:: normal_breach_width_avg_water_depth, normal_breach_width_waterline
       REAL*4:: normal_water_depth_width_avg, froude_number
       REAL*4:: adap_length_flow, adap_length_sediment, length_slope
       REAL*4:: slope_loc
       REAL*4:: water_depth_slope, flow_velocity_slope
       REAL*4:: breach_width_avg_water_depth_slope, hydraulic_radius_slope, friction_coeff_slope
       REAL*4:: theta, stt, n
       REAL*4:: t1
   
       logical :: FLOWSLOPE,sediment_density,water_density
   
       type(NormalFlow) :: results_t1
   
       type(NormalFlow) :: nf
       !real(dp) :: normal_water_depth, normal_flow_velocity, friction_coeff_out
   
       type(NormalFlow) :: res
       !real(dp) :: water_depth, slope_flow_velocity
   
       type(NormalFlow) :: sediment_transport_capacity
       !real(dp) :: n, stt
   
   
       beta = (beta1 + beta0) / 2  ! average inner slope angle (radians)
   

       ! Compute critical flow depth
       crit_water_depth = calc_crit_water_depth(outside_water_level, breach_level, breach_bottom, gamma0)
       !write(logstr,'(a,f6.4)') 'crit_water_depth:', crit_water_depth
       !call write_log(logstr,1)
   
       ! Compute breach widths for critical flow depth
       crit_breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, crit_water_depth,gamma0)
       crit_breach_width_waterline = calc_breach_width_waterline(breach_bottom, crit_water_depth, gamma0)

       ! Compute critical flow velocity
       crit_flow_velocity = calc_crit_flow_velocity(crit_water_depth, crit_breach_width_avg_water_depth,  crit_breach_width_waterline)
       write(logstr,'(a,f6.4)') 'crit_flow_velocity:', crit_flow_velocity
       call write_log(logstr,0)
   
       discharge_coeff = calc_discharge_coeff('Zerihun2020', alpha, beta, outside_water_level, breach_level, W)
       write(logstr,'(a,f6.4)') 'discharge_coeff:', discharge_coeff
       call write_log(logstr,0)

       ! Compute discharge
       discharge = calc_discharge(crit_breach_width_avg_water_depth, crit_water_depth, crit_flow_velocity, discharge_coeff)
       write(logstr,'(a,f12.5)') 'discharge:', discharge
       call write_log(logstr,0)

       ! Compute breach widths
       breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, crit_water_depth, gamma0)
       breach_width_waterline = calc_breach_width_waterline(breach_bottom, crit_water_depth, gamma0)
       breach_width_total = calc_breach_width_total(breach_bottom, crest_level, breach_level, gamma0)
   
       write(logstr,'(a,f12.4,a,f12.4,a,f12.4)') 'breach_width_avg_water_depth:', breach_width_avg_water_depth, ' breach_width_waterline:', breach_width_waterline, ' breach_width_total:', breach_width_total
       call write_log(logstr,0)

       nf = normal_flow_conditions(outside_water_level, breach_level, breach_bottom, discharge, gamma0, beta, delta, d50, d90, Cf, kappa)
    
       normal_water_depth   = nf%normal_water_depth
   
       write(logstr,'(a,f12.4)') 'normal_water_depth:', normal_water_depth
       call write_log(logstr,0)
       normal_flow_velocity = nf%normal_flow_velocity
   
       write(logstr,'(a,f12.4)') 'normal_flow_velocity:', normal_flow_velocity
       call write_log(logstr,0)
       friction_coeff   = nf%friction_coeff

       ! Compute breach width for normal flow depth
       normal_breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, normal_water_depth, gamma0)
       normal_breach_width_waterline = calc_breach_width_waterline(breach_bottom, normal_water_depth, gamma0)
   
       !write(logstr,'(a,f6.4,a,f6.4,a,f6.4)') 'friction_coeff:', friction_coeff, ' normal_breach_width_avg_water_depth:', normal_breach_width_avg_water_depth, ' normal_breach_width_waterline:', normal_breach_width_waterline
       !call write_log(logstr,1)

       ! Compute the width-averaged normal water depth
       normal_water_depth_width_avg = calc_water_depth_width_avg(normal_water_depth, normal_breach_width_avg_water_depth, normal_breach_width_waterline)

       ! Compute Froude number
       froude_number = calc_froude_number(normal_flow_velocity, normal_water_depth_width_avg, beta)

       ! Compute adaptation length to reach normal flow conditions
       adap_length_flow = calc_adaptation_length_flow(froude_number, normal_water_depth, beta)

       !write(logstr,'(a,f6.4,a,f6.4,a,f6.4)') 'normal_water_depth_width_avg:', normal_water_depth_width_avg, ' froude_number:', froude_number, ' adap_length_flow:', adap_length_flow
       !call write_log(logstr,1)
       ! Compute adaptation length to reach equilibrium sediment transport
       adap_length_sediment = calc_adaptation_length_sediment(discharge, normal_breach_width_avg_water_depth, sediment_fall_velocity, beta1)
       adap_length_sediment = adap_length_sediment * breach_width_waterline / breach_width_total

       ! Sediment capacity adaptation length cannot be smaller than normal flow
       ! adaptation length (as long as velocity increases, capacity increases)
       adap_length_sediment = max(adap_length_sediment, adap_length_flow)
       ! Average inner slope length from breach level (Zbr) to inner toe (here: Zp)
       length_slope = (breach_level - polder_level) / sin(beta)

       ! If normal flow conditions are not reached at the inner toe, set the
       ! location for flow conditions that determine sediment transport at inner toe.
       slope_loc = min(adap_length_flow, length_slope - (polder_water_level / sin(beta)))
   
       !write(logstr,'(a,f6.4,a,f6.4,a,f6.4)') 'adap_length_sediment:', adap_length_sediment, ' length_slope:', length_slope, ' slope_loc:', slope_loc
       !call write_log(logstr,1)

       ! Determine the sediment transport. If FLOWSLOPE == True, determine sediment transport at slope_loc.
       if (FLOWSLOPE) then
            res = flow_along_slope(discharge, crit_water_depth, normal_water_depth, breach_bottom, gamma0, adap_length_flow, slope_loc)
            water_depth_slope = res%water_depth
            flow_velocity_slope = res%slope_flow_velocity
        
            breach_width_avg_water_depth_slope = calc_breach_width_avg_water_depth(breach_bottom, water_depth_slope,gamma0)
            hydraulic_radius_slope = calc_hydraulic_radius(water_depth_slope, breach_bottom, breach_width_avg_water_depth_slope, gamma0)
            friction_coeff_slope = iter_friction_coefficient(friction_coeff, flow_velocity_slope, hydraulic_radius_slope, delta, d50, d90, kappa)
            theta = calc_shields_parameter(flow_velocity_slope, delta, d50, friction_coeff_slope)
            sediment_transport_capacity = bagnold_visser(delta, p, d50, phi, sediment_fall_velocity, flow_velocity_slope, beta, friction_coeff_slope)
            n= sediment_transport_capacity%n
            stt = sediment_transport_capacity%stt 
            ![stt, n] = sediment_transport_capacity(params, flow_velocity_slope, water_depth_slope, beta, adap_length_sediment, theta, friction_coeff_slope, params.formula_stc1)
            ! Compute time duration of stage I
            t1 = t0 + ((breach_width_total / breach_width_waterline) * (1 - p) * adap_length_sediment * (beta1 - beta0) * slope_loc / (stt * (1 + 5 * n * (slope_loc / adap_length_flow) *(normal_flow_velocity / flow_velocity_slope - 1))))
       else
            theta = calc_shields_parameter(normal_flow_velocity, delta, d50, friction_coeff) 
            sediment_transport_capacity = bagnold_visser(delta, p, d50, phi, sediment_fall_velocity, normal_flow_velocity, beta, friction_coeff)
            n= sediment_transport_capacity%n
            stt = sediment_transport_capacity%stt
            !write(logstr,'(a,f6.5,a,f6.5)') 'theta:', theta, ' sediment_transport_capacity:', sediment_transport_capacity
            !call write_log(logstr,1)
        
            write(logstr,'(a,f6.4,a,G12.6)') ' n: ', n, ' stt:', stt
            call write_log(logstr,0)
            ![stt, n] = sediment_transport_capacity(params, normal_flow_velocity, normal_water_depth, beta, adap_length_sediment, theta, friction_coeff, params.formula_stc1)
            ! Compute time duration of stage I
            t1 = t0 + ((breach_width_total / breach_width_waterline) * (1 - p) * adap_length_sediment * (beta1 - beta0) * slope_loc / stt)
            !write(logstr,'(a,f10.4,a,f10.4)') 't1:', t1, ' t0:', t0
            !call write_log(logstr,1)
        
            !write(logstr,'(a,f6.4,a,f6.4,a,f6.4)') 'p:', p, ' beta1:', beta1, ' beta0:', beta0
            !call write_log(logstr,1)
       end if
   
       results_t1%discharge = discharge
       results_t1%t1  = t1
       results_t1%breach_width_waterline = breach_width_waterline
       results_t1%breach_width_total  = breach_width_total
       end function
   
       function stage_2(t1, polder_level, breach_bottom, breach_level, polder_water_level,beta1, outside_water_level, gamma0, alpha, W, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity) result(results_t2) !  
       implicit none
       REAL*4:: t1, polder_level, polder_water_level, breach_bottom, breach_level, beta1
       REAL*4:: outside_water_level, gamma0
       REAL*4:: alpha, W, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity
       REAL*4:: beta, crit_water_depth
       REAL*4:: crit_breach_width_avg_water_depth, crit_breach_width_waterline
       REAL*4:: crit_flow_velocity, discharge_coeff, discharge
       REAL*4:: breach_width_avg_water_depth, breach_width_waterline, breach_width_total
       REAL*4:: normal_water_depth, normal_flow_velocity, friction_coeff
       REAL*4:: normal_breach_width_avg_water_depth, normal_breach_width_waterline
       REAL*4:: normal_water_depth_width_avg, froude_number
       REAL*4:: adap_length_flow, adap_length_sediment, length_slope
       REAL*4:: slope_loc
       REAL*4:: water_depth_slope, flow_velocity_slope
       REAL*4:: breach_width_avg_water_depth_slope, hydraulic_radius_slope, friction_coeff_slope
       REAL*4:: theta, stt, n
       REAL*4:: t2
   
       logical :: FLOWSLOPE,sediment_density,water_density
   
       type(NormalFlow) :: results_t2
   
       type(NormalFlow) :: nf
       !real(dp) :: normal_water_depth, normal_flow_velocity, friction_coeff_out
   
       type(NormalFlow) :: res
       !real(dp) :: water_depth, slope_flow_velocity
   
       type(NormalFlow) :: sediment_transport_capacity
       !real(dp) :: n, stt
   

       ! Compute critical flow depth
       crit_water_depth = calc_crit_water_depth(outside_water_level, breach_level, breach_bottom, gamma0)
       !write(logstr,'(a,f6.4)') 'crit_water_depth:', crit_water_depth
       !call write_log(logstr,1)
   
       ! Compute breach widths for critical flow depth
       crit_breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, crit_water_depth,gamma0)
       crit_breach_width_waterline = calc_breach_width_waterline(breach_bottom, crit_water_depth, gamma0)

       ! Compute critical flow velocity
       crit_flow_velocity = calc_crit_flow_velocity(crit_water_depth, crit_breach_width_avg_water_depth,  crit_breach_width_waterline)
       !write(logstr,'(a,f6.4)') 'crit_flow_velocity:', crit_flow_velocity
       !call write_log(logstr,1)
   
       discharge_coeff = calc_discharge_coeff('Zerihun2020', alpha, beta1, outside_water_level, breach_level, W/2)
       !write(logstr,'(a,f6.4)') 'discharge_coeff:', discharge_coeff
       !call write_log(logstr,1)

       ! Compute discharge
       discharge = calc_discharge(crit_breach_width_avg_water_depth, crit_water_depth, crit_flow_velocity, discharge_coeff)
       !write(logstr,'(a,f6.5)') 'discharge:', discharge
       ! call write_log(logstr,1)

       ! Compute breach widths
       breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, crit_water_depth, gamma0)
       breach_width_waterline = calc_breach_width_waterline(breach_bottom, crit_water_depth, gamma0)
       breach_width_total = calc_breach_width_total(breach_bottom, crest_level, breach_level, gamma0)
   
       !write(logstr,'(a,f6.4,a,f6.4,a,f6.4)') 'breach_width_avg_water_depth:', breach_width_avg_water_depth, ' breach_width_waterline:', breach_width_waterline, ' breach_width_total:', breach_width_total
       !call write_log(logstr,1)

       nf = normal_flow_conditions(outside_water_level, breach_level, breach_bottom, discharge, gamma0, beta1, delta, d50, d90, Cf, kappa)
    
       normal_water_depth   = nf%normal_water_depth
   
       !write(logstr,'(a,f6.4)') 'normal_water_depth:', normal_water_depth
       !call write_log(logstr,1)
       normal_flow_velocity = nf%normal_flow_velocity
   
       !write(logstr,'(a,f6.4)') 'normal_flow_velocity:', normal_flow_velocity
       !call write_log(logstr,1)
       friction_coeff   = nf%friction_coeff

       ! Compute breach width for normal flow depth
       normal_breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, normal_water_depth, gamma0)
       normal_breach_width_waterline = calc_breach_width_waterline(breach_bottom, normal_water_depth, gamma0)
   
       !write(logstr,'(a,f6.4,a,f6.4,a,f6.4)') 'friction_coeff:', friction_coeff, ' normal_breach_width_avg_water_depth:', normal_breach_width_avg_water_depth, ' normal_breach_width_waterline:', normal_breach_width_waterline
       !call write_log(logstr,1)

       ! Compute the width-averaged normal water depth
       normal_water_depth_width_avg = calc_water_depth_width_avg(normal_water_depth, normal_breach_width_avg_water_depth, normal_breach_width_waterline)

       ! Compute Froude number
       froude_number = calc_froude_number(normal_flow_velocity, normal_water_depth_width_avg, beta1)

       ! Compute adaptation length to reach normal flow conditions
       adap_length_flow = calc_adaptation_length_flow(froude_number, normal_water_depth, beta1)

       !write(logstr,'(a,f6.4,a,f6.4,a,f6.4)') 'normal_water_depth_width_avg:', normal_water_depth_width_avg, ' froude_number:', froude_number, ' adap_length_flow:', adap_length_flow
       !call write_log(logstr,1)
       ! Compute adaptation length to reach equilibrium sediment transport
       adap_length_sediment = calc_adaptation_length_sediment(discharge, normal_breach_width_avg_water_depth, sediment_fall_velocity, beta1)
       adap_length_sediment = adap_length_sediment * breach_width_waterline / breach_width_total

       ! Sediment capacity adaptation length cannot be smaller than normal flow
       ! adaptation length (as long as velocity increases, capacity increases)
       adap_length_sediment = max(adap_length_sediment, adap_length_flow)
       ! Average inner slope length from breach level (Zbr) to inner toe (here: Zp)
       length_slope = (breach_level - polder_level) / sin(beta1)
       !write(logstr,'(a,f12.4,a,f6.4)') 'polder_level:', polder_level, ' beta1:', beta1
       !call write_log(logstr,1)
       ! If normal flow conditions are not reached at the inner toe, set the
       ! location for flow conditions that determine sediment transport at inner toe.
       slope_loc = min(adap_length_flow, max(length_slope - (polder_water_level / sin(beta1)),0.0))
   
       !write(logstr,'(a,f6.4,a,f6.4,a,f6.4)') 'adap_length_sediment:', adap_length_sediment, ' length_slope:', length_slope, ' slope_loc:', slope_loc
       !call write_log(logstr,1)

       ! Determine the sediment transport. If FLOWSLOPE == True, determine sediment transport at slope_loc.
       if (FLOWSLOPE) then
            res = flow_along_slope(discharge, crit_water_depth, normal_water_depth, breach_bottom, gamma0, adap_length_flow, slope_loc)
            water_depth_slope = res%water_depth
            flow_velocity_slope = res%slope_flow_velocity
        
            breach_width_avg_water_depth_slope = calc_breach_width_avg_water_depth(breach_bottom, water_depth_slope,gamma0)
            hydraulic_radius_slope = calc_hydraulic_radius(water_depth_slope, breach_bottom, breach_width_avg_water_depth_slope, gamma0)
            friction_coeff_slope = iter_friction_coefficient(friction_coeff, flow_velocity_slope, hydraulic_radius_slope, delta, d50, d90, kappa)
            theta = calc_shields_parameter(flow_velocity_slope, delta, d50, friction_coeff_slope)
            sediment_transport_capacity = bagnold_visser(delta, p, d50, phi, sediment_fall_velocity, flow_velocity_slope, beta1, friction_coeff_slope)
            n= sediment_transport_capacity%n
            stt = sediment_transport_capacity%stt 
       else
            theta = calc_shields_parameter(normal_flow_velocity, delta, d50, friction_coeff) 
            sediment_transport_capacity = bagnold_visser(delta, p, d50, phi, sediment_fall_velocity, normal_flow_velocity, beta1, friction_coeff)
            n= sediment_transport_capacity%n
            stt = sediment_transport_capacity%stt
            !write(logstr,'(a,f6.5,a,f6.5)') 'theta:', theta, ' sediment_transport_capacity:', sediment_transport_capacity
            !call write_log(logstr,1)
        
            write(logstr,'(a,f6.4,a,G12.6)') ' n: ', n, ' stt:', stt
            call write_log(logstr,0)

       t2 = t1 + ((breach_width_total / breach_width_waterline) * W * (1 - p) * adap_length_sediment* sin(beta1) / stt)
       write(logstr,'(a,G12.6,a,G12.6)') ' t1: ', t1, ' W: ', W
       call write_log(logstr,0)
       end if
   
       results_t2%discharge = discharge
       results_t2%t2  = t2
       results_t2%breach_width_waterline = breach_width_waterline
       results_t2%breach_width_total  = breach_width_total
       end function
   
       function stage_3(dt, breach_width_total, breach_width_waterline, theta_crit, ni, dstar, k, rhos, rhow, outside_level, polder_level, breach_bottom, breach_level, polder_water_level,beta1, outside_water_level, gamma0,gamma1, alpha, W, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity) result(results_t3) !  
       implicit none
       REAL*4:: dt, theta_crit, ni, dstar, k, rhos, rhow, outside_level, polder_level, polder_water_level, breach_bottom, breach_level, beta1
       REAL*4:: outside_water_level, gamma0, gamma1,gamma
       REAL*4:: alpha, W, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity
       REAL*4:: beta, crit_water_depth
       REAL*4:: crit_breach_width_avg_water_depth, crit_breach_width_waterline
       REAL*4:: crit_flow_velocity, discharge_coeff, discharge
       REAL*4:: breach_width_avg_water_depth, breach_width_waterline, breach_width_total
       REAL*4:: normal_water_depth, normal_flow_velocity, friction_coeff
       REAL*4:: normal_breach_width_avg_water_depth, normal_breach_width_waterline
       REAL*4:: normal_water_depth_width_avg, froude_number
       REAL*4:: adap_length_flow, adap_length_sediment, length_slope
       REAL*4:: slope_loc
       REAL*4:: water_depth_slope, flow_velocity_slope
       REAL*4:: breach_width_avg_water_depth_slope, hydraulic_radius_slope, friction_coeff_slope
       REAL*4:: theta, stt, n, stt2, n2
   
       logical :: FLOWSLOPE,sediment_density,water_density
   
       type(NormalFlow) :: results_t3
   
       type(NormalFlow) :: nf
       !real(dp) :: normal_water_depth, normal_flow_velocity, friction_coeff_out
   
       type(NormalFlow) :: res
       !real(dp) :: water_depth, slope_flow_velocity
   
       type(NormalFlow) :: sediment_transport_capacity, sediment_transport_capacity2
       !real(dp) :: n, stt
   
   
       crit_water_depth = calc_crit_water_depth(outside_water_level, breach_level, breach_bottom, gamma0)
       !write(logstr,'(a,f6.4)') 'crit_water_depth:', crit_water_depth
       !call write_log(logstr,1)
        
       if ((polder_water_level - breach_level) >= crit_water_depth) then
           write(logstr,'(a,f12.4,a,f12.4,a,f12.4)') 'Stopping simulation, polder water level: ', polder_water_level, ' breach level: ', breach_level, ' crit_water_depth: ', crit_water_depth
            call write_log(logstr,1)
            error stop "Subcritical flow in stage 3. Stopping simulation."
       end if
   
       if (crit_water_depth < 0) then
            write(logstr, '(a,f6.4)') "Water depth is negative in calc_breach_width_avg_water_depth stage 3: ", crit_water_depth
            error stop logstr
       end if
       crit_breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, crit_water_depth,gamma0)
       crit_breach_width_waterline = calc_breach_width_waterline(breach_bottom, crit_water_depth, gamma0)

       ! Compute critical flow velocity
       crit_flow_velocity = calc_crit_flow_velocity(crit_water_depth, crit_breach_width_avg_water_depth,  crit_breach_width_waterline)
       !write(logstr,'(a,f6.4)') 'crit_flow_velocity:', crit_flow_velocity
       !call write_log(logstr,1)
   
       discharge_coeff = calc_discharge_coeff('Zerihun2020', alpha, beta1, outside_water_level, breach_level, crest_level-breach_level)
       !write(logstr,'(a,f6.4)') 'discharge_coeff:', discharge_coeff
       !call write_log(logstr,1)
        
       discharge_coeff = submerged_weir_flow (outside_water_level, breach_level, polder_water_level, discharge_coeff)
        
       ! Compute discharge
       discharge = calc_discharge(crit_breach_width_avg_water_depth, crit_water_depth, crit_flow_velocity, discharge_coeff)
       !write(logstr,'(a,f6.5)') 'discharge:', discharge
       !call write_log(logstr,1)
        
       nf = normal_flow_conditions(outside_water_level, breach_level, breach_bottom, discharge, gamma0, beta1, delta, d50, d90, Cf, kappa)
    
       normal_water_depth   = nf%normal_water_depth
   
       !write(logstr,'(a,f6.4)') 'normal_water_depth:', normal_water_depth
       !call write_log(logstr,1)
       normal_flow_velocity = nf%normal_flow_velocity
   
       !write(logstr,'(a,f6.4)') 'normal_flow_velocity:', normal_flow_velocity
       !call write_log(logstr,1)
       friction_coeff   = nf%friction_coeff

       ! Compute breach width for normal flow depth
       normal_breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, normal_water_depth, gamma0)
       normal_breach_width_waterline = calc_breach_width_waterline(breach_bottom, normal_water_depth, gamma0)
   
       !write(logstr,'(a,f6.4,a,f6.4,a,f6.4)') 'friction_coeff:', friction_coeff, ' normal_breach_width_avg_water_depth:', normal_breach_width_avg_water_depth, ' normal_breach_width_waterline:', normal_breach_width_waterline
       !call write_log(logstr,1)

       ! Compute the width-averaged normal water depth
       normal_water_depth_width_avg = calc_water_depth_width_avg(normal_water_depth, normal_breach_width_avg_water_depth, normal_breach_width_waterline)

       ! Compute Froude number
       froude_number = calc_froude_number(normal_flow_velocity, normal_water_depth_width_avg, beta1)

       ! Compute adaptation length to reach normal flow conditions
       adap_length_flow = calc_adaptation_length_flow(froude_number, normal_water_depth, beta1)

       !write(logstr,'(a,f6.4,a,f6.4,a,f6.4)') 'normal_water_depth_width_avg:', normal_water_depth_width_avg, ' froude_number:', froude_number, ' adap_length_flow:', adap_length_flow
       !call write_log(logstr,1)
        
       ! Compute adaptation length to reach equilibrium sediment transport
       adap_length_sediment = calc_adaptation_length_sediment(discharge, normal_breach_width_avg_water_depth, sediment_fall_velocity, beta1)
       adap_length_sediment = adap_length_sediment * breach_width_waterline / breach_width_total
       !write(logstr,'(a,G12.4,a,G12.4,a,G12.4)') 'adap_length_sediment:', adap_length_sediment, ' breach_width_waterline:', breach_width_waterline, ' breach_width_total:', breach_width_total
       !call write_log(logstr,1)
       ! Sediment capacity adaptation length cannot be smaller than normal flow
       ! adaptation length (as long as velocity increases, capacity increases)
       adap_length_sediment = max(adap_length_sediment, adap_length_flow)
       ! Average inner slope length from breach level (Zbr) to inner toe (here: Zp)
       length_slope = (breach_level - polder_level) / sin(beta1)
       !write(logstr,'(a,f12.4,a,f6.4)') 'polder_level:', polder_level, ' beta1:', beta1
       !call write_log(logstr,1)
       ! If normal flow conditions are not reached at the inner toe, set the
       ! location for flow conditions that determine sediment transport at inner toe.
       slope_loc = min(adap_length_flow, max(length_slope - (polder_water_level / sin(beta1)),0.0))
   
       write(logstr,'(a,G12.4,a,G12.4,a,G12.4)') 'adap_length_sediment: ', adap_length_sediment, ' breach_width_waterline:', breach_width_waterline, ' breach_width_total:', breach_width_total
       call write_log(logstr,0)
        
       ! Determine the sediment transport. If FLOWSLOPE == True, determine sediment transport at slope_loc.
       if (FLOWSLOPE) then
            res = flow_along_slope(discharge, crit_water_depth, normal_water_depth, breach_bottom, gamma0, adap_length_flow, slope_loc)
            water_depth_slope = res%water_depth
            flow_velocity_slope = res%slope_flow_velocity
        
            breach_width_avg_water_depth_slope = calc_breach_width_avg_water_depth(breach_bottom, water_depth_slope,gamma0)
            hydraulic_radius_slope = calc_hydraulic_radius(water_depth_slope, breach_bottom, breach_width_avg_water_depth_slope, gamma0)
            friction_coeff_slope = iter_friction_coefficient(friction_coeff, flow_velocity_slope, hydraulic_radius_slope, delta, d50, d90, kappa)
            theta = calc_shields_parameter(flow_velocity_slope, delta, d50, friction_coeff_slope)
            
            
            sediment_transport_capacity = vanrhee_simp(theta, theta_crit, adap_length_sediment, water_depth_slope, ni, p, delta, d50, dstar, k, rhos, rhow) !params, theta, adap_length_sediment, water_depth
            n= sediment_transport_capacity%n
            stt = sediment_transport_capacity%stt 
       else
            theta = calc_shields_parameter(normal_flow_velocity, delta, d50, friction_coeff) 
            
            sediment_transport_capacity = vanrhee_simp(theta, theta_crit, adap_length_sediment, normal_water_depth, ni, p, delta, d50, dstar, k, rhos, rhow)
            n= sediment_transport_capacity%n
            stt = sediment_transport_capacity%stt
            !write(logstr,'(a,f6.5,a,f6.5)') 'theta:', theta, ' sediment_transport_capacity:', sediment_transport_capacity
            !call write_log(logstr,1)
        
            !write(logstr,'(a,f6.4,a,G12.6)') ' n:', n, ' stt:', stt
            !call write_log(logstr,1)
       end if
   
       ! Compute new breach level
       breach_level = breach_level - ((breach_width_waterline / breach_width_total)* sin(alpha) / sin(alpha + beta1) * stt * dt / ((1 - p) * adap_length_sediment))

       ! Compute new breach side slope angle (radians)
       if (gamma0 < gamma1) then
            gamma = atan((crest_level - breach_level) / (0.5 * (breach_width_total - breach_bottom)))
            gamma0 = min(gamma, gamma1)
       end if
        
       ! Increase breach width due to erosion of the side slopes at breach entry (upstream)
       theta = calc_shields_parameter(crit_flow_velocity, delta, d50, friction_coeff)
       sediment_transport_capacity2 = vanrhee_simp(theta, theta_crit, adap_length_sediment, crit_water_depth, ni, p, delta, d50, dstar, k, rhos, rhow)
       n2= sediment_transport_capacity2%n
       stt2 = sediment_transport_capacity2%stt
            

       breach_bottom = (breach_bottom + 2 * (crit_water_depth / (crest_level - breach_level)) * stt2 * dt / ((1 - p) * adap_length_sediment * tan(gamma0)))

       ! Compute breach widths
       if (crit_water_depth < 0) then
            write(logstr, '(a,f6.4)') "Water depth is negative in calc_breach_width_avg_water_depth stage 3: ", crit_water_depth
            error stop logstr
       end if
       breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, crit_water_depth, gamma0)
       breach_width_waterline = calc_breach_width_waterline(breach_bottom, crit_water_depth, gamma0)
       breach_width_total = calc_breach_width_total(breach_bottom, crest_level, breach_level, gamma0)
        
   
       results_t3%discharge = discharge
       results_t3%breach_width_waterline = breach_width_waterline
       results_t3%breach_width_total  = breach_width_total
       results_t3%breach_level  = breach_level
       results_t3%breach_bottom = breach_bottom
       results_t3%gamma0 = gamma0
       end function
   
       function stage_4(dt, breach_width_total, breach_width_waterline, theta_crit, ni, dstar, k, rhos, rhow, outside_level, polder_level, breach_bottom, breach_level, polder_water_level, outside_water_level, gamma1, alpha, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity) result(results_t4) !  
       implicit none
       REAL*4:: dt, theta_crit, ni, dstar, k, rhos, rhow, outside_level, polder_level, polder_water_level, breach_bottom, breach_level
       REAL*4:: outside_water_level, gamma1
       REAL*4:: alpha, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity
       REAL*4:: crit_water_depth
       REAL*4:: crit_breach_width_avg_water_depth, crit_breach_width_waterline
       REAL*4:: crit_flow_velocity, discharge_coeff, discharge
       REAL*4:: breach_width_avg_water_depth, breach_width_waterline, breach_width_total
       REAL*4:: adap_length_sediment
       REAL*4:: flow_velocity
       REAL*4:: hydraulic_radius, friction_coeff
       REAL*4:: theta, stt, n
   
       type(NormalFlow) :: results_t4
   
       type(NormalFlow) :: sediment_transport_capacity   
   
       crit_water_depth = calc_crit_water_depth(outside_water_level, breach_level, breach_bottom, gamma1)
       write(logstr,'(a,f6.4)') 'crit_water_depth:', crit_water_depth
       call write_log(logstr,0)
        
       if ((polder_water_level - breach_level) >= crit_water_depth) then
           write(logstr,'(a,f12.4,a,f12.4,a,f12.4)') 'Stopping simulation, polder water level: ', polder_water_level, ' breach level: ', breach_level, ' crit_water_depth: ', crit_water_depth
            call write_log(logstr,1)
            error stop "Subcritical flow in stage 3. Stopping simulation."
       end if
        
   
       if (crit_water_depth < 0) then
            write(logstr, '(a,f6.4)') "Water depth is negative in calc_breach_width_avg_water_depth stage 4: ", crit_water_depth
            error stop logstr
       end if
       crit_breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, crit_water_depth,gamma1)
       crit_breach_width_waterline = calc_breach_width_waterline(breach_bottom, crit_water_depth, gamma1)

       ! Compute critical flow velocity
       crit_flow_velocity = calc_crit_flow_velocity(crit_water_depth, crit_breach_width_avg_water_depth,  crit_breach_width_waterline)
       !write(logstr,'(a,f6.4)') 'crit_flow_velocity:', crit_flow_velocity
       !call write_log(logstr,1)
   
       discharge_coeff = 1.3
        
       ! Compute discharge
       discharge = calc_discharge(crit_breach_width_avg_water_depth, crit_water_depth, crit_flow_velocity, discharge_coeff)
       !write(logstr,'(a,f6.5)') 'discharge:', discharge
       !call write_log(logstr,1)
   
       ! Compute flow velocity in breach
       flow_velocity = discharge / (crit_breach_width_avg_water_depth * crit_water_depth)
    
        
       hydraulic_radius = calc_hydraulic_radius(crit_water_depth, breach_bottom, crit_breach_width_avg_water_depth, gamma1)
       friction_coeff = iter_friction_coefficient(Cf, flow_velocity, hydraulic_radius, delta, d50, d90, kappa)
       theta = calc_shields_parameter(flow_velocity, delta, d50, friction_coeff)
   
       adap_length_sediment = (0.4 * crit_water_depth / (crest_level - breach_level) * flow_velocity *crit_water_depth / sediment_fall_velocity)
            
            
       sediment_transport_capacity = vanrhee_simp(theta, theta_crit, adap_length_sediment, crit_water_depth, ni, p, delta, d50, dstar, k, rhos, rhow) !params, theta, adap_length_sediment, water_depth
       n= sediment_transport_capacity%n
       stt = sediment_transport_capacity%stt 
        
       breach_bottom = (breach_bottom + 2.0 * (crit_water_depth / (crest_level - breach_level)) * stt * dt / ((1 - p) * adap_length_sediment * tan(gamma1)))

       ! Compute breach widths
       if (crit_water_depth < 0) then
            write(logstr, '(a,f6.4)') "#2 Water depth is negative in calc_breach_width_avg_water_depth stage 4: ", crit_water_depth
            error stop logstr
       end if
       breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, crit_water_depth, gamma1)
       breach_width_waterline = calc_breach_width_waterline(breach_bottom, crit_water_depth, gamma1)
       breach_width_total = calc_breach_width_total(breach_bottom, crest_level, breach_level, gamma1)
    
       results_t4%discharge = discharge
       results_t4%breach_width_waterline = breach_width_waterline
       results_t4%breach_width_total  = breach_width_total
       results_t4%breach_level  = breach_level
       results_t4%breach_bottom = breach_bottom
       end function
   
       function stage_5(dt, breach_width_total, breach_width_waterline, theta_crit, beta1, ni, dstar, k, rhos, rhow, outside_level, polder_level, breach_bottom, breach_level, polder_water_level, outside_water_level, gamma1, alpha, crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity) result(results_t5) !  
       implicit none
       REAL*4:: dt, theta_crit, ni, dstar, k, rhos, rhow, outside_level, polder_level, polder_water_level, breach_bottom, breach_level
       REAL*4:: outside_water_level, gamma1, alpha, beta1
       REAL*4:: crest_level, d50, d90, Cf, kappa, delta, p, phi, sediment_fall_velocity
       REAL*4:: water_depth
       REAL*4:: discharge_coeff, discharge
       REAL*4:: breach_width_avg_water_depth, breach_width_waterline, breach_width_total
       REAL*4:: adap_length_sediment
       REAL*4:: flow_velocity
       REAL*4:: hydraulic_radius, friction_coeff
       REAL*4:: theta, stt, n, g
   
       type(NormalFlow) :: results_t5
   
       type(NormalFlow) :: sediment_transport_capacity   
   
       ! Compute flow velocity in breach
       g =9.81
       flow_velocity = (2.0 * g * (outside_water_level - polder_water_level)) ** 0.5
   
       ! Calculate the water depth in the breach
       water_depth = polder_water_level - breach_level
       write(logstr,'(a,f6.4)') 'water_depth:', water_depth
       call write_log(logstr,0)
   
       if (water_depth < 0) then
            write(logstr, '(a,f6.4)') "Water depth is negative in calc_breach_width_avg_water_depth stage 5: ", water_depth
            error stop logstr
       end if
       
       breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, water_depth, gamma1)
   
       if (breach_level > polder_level) then
           discharge_coeff = calc_discharge_coeff('Zerihun2020', alpha, beta1, outside_water_level, breach_level, crest_level-breach_level)
           discharge_coeff = submerged_weir_flow (outside_water_level, breach_level, polder_water_level, discharge_coeff)
       else
           discharge_coeff = 1.3
       end if
   
       ! Calculate the discharge
       discharge = discharge_coeff * breach_width_avg_water_depth * water_depth * flow_velocity
    
       theta = calc_shields_parameter(flow_velocity, delta, d50, Cf)
   
       if (theta > theta_crit) then
            ! Compute flow velocity in breach
            flow_velocity = discharge / (breach_width_avg_water_depth * water_depth)
            hydraulic_radius = calc_hydraulic_radius(water_depth, breach_bottom, breach_width_avg_water_depth, gamma1)
            friction_coeff = iter_friction_coefficient(Cf, flow_velocity, hydraulic_radius, delta, d50, d90, kappa)
            theta = calc_shields_parameter(flow_velocity, delta, d50, friction_coeff)
            adap_length_sediment = (0.4 * water_depth / (crest_level - breach_level) * flow_velocity * water_depth / sediment_fall_velocity)
       
            sediment_transport_capacity = vanrhee_simp(theta, theta_crit, adap_length_sediment, water_depth, ni, p, delta, d50, dstar, k, rhos, rhow) !params, theta, adap_length_sediment, water_depth
            n= sediment_transport_capacity%n
            stt = sediment_transport_capacity%stt 
        
            breach_bottom = (breach_bottom + 2.0 * (water_depth / (crest_level - breach_level)) * stt * dt / ((1 - p) * adap_length_sediment * tan(gamma1)))
       end if     
        
   

       ! Compute breach widths
       if (water_depth < 0) then
            write(logstr, '(a,f6.4)') " #2 Water depth is negative in calc_breach_width_avg_water_depth stage 5: ", water_depth
            error stop logstr
       end if
       breach_width_avg_water_depth = calc_breach_width_avg_water_depth(breach_bottom, water_depth, gamma1)
       breach_width_waterline = calc_breach_width_waterline(breach_bottom, water_depth, gamma1)
       breach_width_total = calc_breach_width_total(breach_bottom, crest_level, breach_level, gamma1)
    
       results_t5%discharge = discharge
       results_t5%breach_width_waterline = breach_width_waterline
       results_t5%breach_width_total  = breach_width_total
       results_t5%breach_level  = breach_level
       results_t5%breach_bottom = breach_bottom
       end function
   end module