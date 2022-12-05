Input parameters
=====

Different parameters for model input and output of SFINCS can be changed in **sfincs.inp**, see below. 

Traditionally SFINCS neglects the advection term in the SFINCS-LIE version ('advection = 0'). 
For super-critical flow conditions or modelling waves, the SFINCS-SSWE version can be used ('advection = 1' for 1D modelling and 'advection = 2' for 2D modelling) for better performance.        
     
Parameters for model input
-----

	mmax
	  :description:		Number of grid cells in x-direction
	  :units:		-
	  :default:		0
	  :min:			1
	  :max:			Inf (recommended is to limit the total number of active cells to max 3 million)
	nmax
	  :description:		Number of grid cells in y-direction
	  :units:		-
	  :default:		0
	  :min:			1
	  :max:			Inf (recommended is to limit the total number of active cells to max 3 million)	  
	dx
	  :description:		Grid size in x-direction
	  :units:		m
	  :default:		0
	  :min:			1.0e-3
	  :max:			Inf (recommended is a maximum grid size of 1000 meters)
	dy
	  :description:		Grid size in y-direction
	  :units:		m
	  :default:		0
	  :min:			1.0e-3
	  :max:			Inf (recommended is a maximum grid size of 1000 meters)	  
	x0
	  :description:		X-coordinate of first grid cell corner (1,1), thus not center of grid cell, in projected UTM zone.
	  :units:		m in projected UTM zone
	  :default:		0
	  :min:			0
	  :max:			Inf 
	y0
	  :description:		Y-coordinate of first grid cell corner (1,1), thus not center of grid cell, in projected UTM zone.
	  :units:		m in projected UTM zone
	  :default:		0
	  :min:			0
	  :max:			Inf 	  
	rotation
	  :description:		Rotation of the grid in degrees from the x-axis (east) in anti-clockwise direction
	  :units:		degrees
	  :default:		0
	  :min:			0
	  :max:			359.999 	  
	advection
	  :description:		setting for advection. 0 for no advection scheme (SFINCS-LIE), 1 for 1D advection scheme for modelling in 1D OR 2 for 2D advection scheme for modelling in 2D (SFINCS-SSWE).
	  :units:		-
	  :default:		0
	  :min:			0
	  :max:			2
	alpha	
	  :description:		CFL-condition reduction. Decrease for additional numerical stability, minimum value is 0.1 and maximum is 0.75.
	  :units:		-	
	  :default:		0.75		
	  :min:			0.1 (recommended)	
	  :max:			0.75 (recommended)		  
	huthresh	
	  :description:		Minimum flow depth limiter.
	  :units:		m
	  :default:		0.05
	  :min:			0.001 (recommended)
	  :max:			0.1 (recommended)
	theta
	  :description:		Smoothing factor in momentum equation. Advised not too change and to use 1.0 by default.
	  :units:		-
	  :default:		1.0
	  :min:			0.8
	  :max:			1.0
	zsini
	  :description:		Initial water level.
	  :units:		m above reference level
	  :default:		0
	  :min:			-Inf
	  :max:			Inf
	qinf
	  :description:		Infiltration rate, specify in +mm/hr.
	  :units:		mm/hr
	  :default:		0
	  :min:			0
	  :max:			100  
	manning
	  :description:		Uniform manning roughness, specify in s/m^(1/3).
	  :units:		s/m^(1/3)
	  :default:		0.04
	  :min:			0
	  :max:			0.1 (advised)  	
	rgh_level_land
	  :description:		Elevation level to distinguish land and sea roughness (when using 'manning_land' and 'manning_sea').
	  :units:		m above reference level
	  :default:		0
	  :min:			-Inf
	  :max:			Inf  		  
	manning_land
	  :description:		Varying manning roughness based on elevation (above 'rgh_level_land', overules uniform 'manning', specify in s/m^(1/3).
	  :units:		s/m^(1/3)
	  :default:		-999 (=not used)
	  :min:			0
	  :max:			0.1 (advised)    		  
	manning_sea
	  :description:		Varying manning roughness based on elevation (below 'rgh_level_land', overules uniform 'manning', specify in s/m^(1/3).
	  :units:		s/m^(1/3)
	  :default:		-999 (=not used)
	  :min:			0
	  :max:			0.1 (advised)    	  
	  
More parameters for model input (only for advanced users)
-----

	bndtype        
	  :description:		Boundary type for interpretation of 'sfincs.bzs' time-series. bndtype=1 is for water levels, bndtype=2 (advanced) is for horizontal velocities (in m/s) and bndtype=3 (advanced) for horizontal discharges (in m2/s).
	  :units:		-
	  :default:		1
	  :min:			1
	  :max:			3
	rhoa
	  :description:		Density of the air
	  :units:		kg/m3
	  :default:		1.25
	  :min:			-
	  :max:			-
	rhow
	  :description:		Density of the water
	  :units:		kg/m3
	  :default:		1024
	  :min:			-
	  :max:			-
	stopdepth
	  :description:		Water depth anywhere in the domain after which the simulation is classified as unstable and stopped
	  :units:		m
	  :default:		100
	  :min:			0
	  :max:			Inf	  
	advlim
	  :description:		Advection limiter when advection>0 to limit the magnitude of the advection term when calculating fluxes between cells.
	  :units:		-
	  :default:		9999
	  :min:			0
	  :max:			9999
	dtmax
	  :description:		Maximum internal time step to be used
	  :units:		s
	  :default:		60
	  :min:			1.0e-3
	  :max:			Inf
	dtmin
	  :description:		Minimum internal time step to be used
	  :units:		s
	  :default:		1.0e-3
	  :min:			1.0e-3
	  :max:			Inf	  
	tspinup
	  :description:		Duration of internal spinup period before tstart
	  :units:		s
	  :default:		60
	  :min:			0
	  :max:			Inf
	  
	**Drag coefficients:**
	
	cdnrb
	  :description:		Number of specified break points
	  :units:		-
	  :default:		3
	  :min:			2
	  :max:			-	
	cdwnd	  
	  :description:		Wind speed break points (including 0)
	  :units:		-
	  :default:		0  28  50
	  :min:			2 values
	  :max:			-
	  :description:		Drag coefficient brak points
	  :units:		-
	  :default:		0.001 0.0025 0.0015
	  :min:			2 values
	  :max:			-	  
	
Different parameters influencing the given output by SFINCS can be changed, see below. 

Parameters for model output
-----

	tref
	  :description:		Reference date in 'yyyymmdd HHMMSS'
	  :units:		-
	  :default:		20000101 000000
	tstart
	  :description:		Start date in 'yyyymmdd HHMMSS'
	  :units:		-	
	  :default:		20000101 000000				  
	tstop
	  :description:		Stop date in 'yyyymmdd HHMMSS'
	  :units:		m
	  :default:		20000101 000000
	trstout
	  :description:		Specific time in seconds since 'tref' for binary restart file output being written away, turned of by default.
	  :units:		s
	  :default:		-999.0	  	  
	dtout
	  :description:		Time-step global map output.
	  :units:		s
	  :default:		0
	dthisout
	  :description:		Time-step observation points output.
	  :units:		s
	  :default:		600
	dtmaxout
	  :description:		Time-step interval of global map output of maximum water level. If not specified, the maximum over the entire simulation is calculated. If no output is wanted, specify 'dtmaxout = 0'.
	  :units:		s
	  :default:		999999
	  :min:			0
	  :max:			'tstop - start in seconds'  
	dtrstout
	  :description:		Time-step for binary restart file output being written away, turned of by default.
	  :units:		s
	  :default:		0	  	  
	dtwnd
	  :description:		Time-interval wind update (only for spiderweb)
	  :units:		s
	  :default:		1800
	outputformat
	  :description:		Choice whether the SFINCS model output is given in binary 'bin', ascii 'asc' or netcdf files 'net' (default). In case of netcdf output, global output is given in 'sfincs_map.nc', point output in 'sfincs_his.nc' in case observation points are specified.
	  :units:		-
	  :default:		net
	twet_threshold
	  :description:		Threshold value of water depth to count cell as flooded for keeping track of wet cells with storetwet = 1
	  :units:		m
	  :default:		0.01	  
	storetwet
	  :description:		Flag to turn on writing away duration that a cell was wet during simulation (storetwet = 1)
	  :units:		-
	  :default:		0	  
	storevel
	  :description:		Flag to turn on writing away velocities on 'dtout' interval during simulation (storevel = 1)
	  :units:		-
	  :default:		0	
	storevelmax
	  :description:		Flag to turn on writing away maximum velocities on 'dtmaxout' interval during simulation (storevelmax = 1)
	  :units:		-
	  :default:		0	
	storecumprcp
	  :description:		Flag to turn on writing away cumulative precipitation on 'dtmaxout' interval during simulation (storecumprcp = 1)
	  :units:		-
	  :default:		0		
	storehsubgrid
	  :description:		Flag to turn on writing away unaccurate water depth estimate for subgrid mode on 'dtmaxout' interval during simulation (storehsubgrid = 1)
	  :units:		-
	  :default:		0		    	  
	storeqdrain
	  :description:		Flag to turn on writing away drainage discharge during simulation (storeqdrain = 1)
	  :units:		-
	  :default:		0	
	storezvolume
	  :description:		Flag to turn on writing away water volumes for the subgrid mode during simulation (storezvolume = 1)
	  :units:		-
	  :default:		0		  
	storemeteo
	  :description:		Flag to turn on writing away meteo input data during simulation (storemeteo = 1)
	  :units:		-
	  :default:		0	
	storemaxwind
	  :description:		Flag to turn on writing away maximum wind speed during simulation (storemaxwind = 1)
	  :units:		-
	  :default:		0	
	debug
	  :description:		Flag to turn on writing away every timestep to output as debug mode (debug = 1)
	  :units:		-
	  :default:		0	

Input files
=====	 

SFINCS consists of many different input files, this overview gives a description, whether they are required or not, unit and format (bin = binary, asc = ascii and net = netcdf).

.. figure:: ./figures/SFINCS_documentation_figure1.png
   :width: 800px
   :align: center

   Overview of input file of SFINCS with indication whther they are required or not	
	

Domain
-----

	sfincs.inp
	  :description:		General input file of SFINCS describing all model settings, the domain, forcing and structures.
	  :required:		yes
	  :format:		asc	 
	depfile = sfincs.dep
	  :description:		Elevation (bathymetry and topography) at grid cell centres above a reference level. 
	  :units:		m above reference level
	  :required:		yes in case of regular mode, no in case of subgrid mode
	  :format:		bin or asc
	mskfile = sfincs.msk
	  :description:		This mask indicates for every cell whether it is an inactive cell (msk=0), active cell (msk=1), boundary cell (msk=2) or outflow boundary cell msk=3).
	  :units:		-
	  :required:		yes	  
	  :format:		bin or asc
	indexfile = sfincs.ind
	  :description:		File describing the indices of active grid cells within the overall grid. Not used by SFINCS with ascii input.
	  :units:		-
	  :required:		Only if 'inputformat = bin'
	  :format:		bin	  
	mskfile = sfincs.msk
	  :description:		This mask indicates for every cell whether it is an inactive cell (msk=0), active cell (msk=1), water level boundary cell (msk=2) or outflow boundary cell msk=3).
	  :units:		-
	  :required:		yes	  
	  :format:		bin or asc	  
	manningfile = sfincs.man
	  :description:		For spatially varying friction values per cell use the manningfile option, with the same grid based input as the depfile using a binary file. Not used by SFINCS in subgrid mode.
	  :units:		s/m^(1/3)
	  :required:		no in case of regular mode, ignored in case of subgrid mode	  
	  :format:		bin	 
	qinffile = sfincs.qinf
	  :description:		For spatially varying constant in time infiltration values per cell use the qinffile option, with the same grid based input as the depfile using a binary file.
	  :units:		mm/hr
	  :required:		no	  
	  :format:		bin	  
	scsfile = sfincs.scs
	  :description:		For spatially varying infiltration values per cell using the Curve Number method use the scsfile option, with the same grid based input as the depfile using a binary file.
	  :units:		-
	  :required:		no	  
	  :format:		bin	  	  
	sbgfile = sfincs.sbg
	  :description:		File containing subgrid tables, only needed by SFINCS if you want to run your model in the subgrid mode.
	  :units:		-
	  :required:		Only for running SFINCS in subgrid mode	  
	  :format:		bin		  
	obsfile = sfincs.obs
	  :description:		To get output time-series at individual point locations, observations points have to be specified.
	  :units:		m in projected UTM zone
	  :required:		no (only if point output is wanted)
	  :format:		asc		  
	crsfile = sfincs.crs
	  :description:		To get output time-series of discharge through a cross-section, cross_sections have to be specified.
	  :units:		m in projected UTM zone
	  :required:		no (only if cross-section output is wanted)
	  :format:		tekal		  	  
	inifile = sfincs.ini
	  :description:		For spatially varying initial water level per cell, with the same grid based input as the depfile using a **binary file**. In older version this was an ascii file still, not from official release v0.0.1 onwards!
	  :units:		m above reference level
	  :required:		no
	  :format:		**bin**		
	rstfile = sfincs.rst
	  :description:		More advanced restartfile that can also contain fluxes and velocities. As produced by SFINCS if dtrstout > 0 OR trstout > 0. Type of restart - 1: zs, qx, qy, umean and vmean  - 2: zs, qx, qy - 3: zs
	  :units:		-
	  :required:		no 
	  :format:		bin	

Forcing - Water levels and waves
-----

	bndfile = sfincs.bnd
	  :description:		To specify water-level time-series to the boundary cells (msk=2), first the input locations have to be specified in 'sfincs.bnd'.
	  :units:		m in projected UTM zone	  
	  :required:		Only when specifying water levels and waves.
	  :format:		asc	 
	bzsfile = sfincs.bzs
	  :description:		In the file 'sfincs.bzs' the (slowly varying) water level time-series are specified per input location. 
	  :units:		m above reference level
	  :required:		Only when specifying water levels.
	  :format:		asc	 	
	bzifile = sfincs.bzi
	  :description:		Tn the file 'sfincs.bzi' the quickly varying water level time-series due to incoming waves are specified per input location. Do note that the input timestep should be the same in both the bzs and bzi files!
	  :units:		m around mean water level of bzsfile
	  :required:		Only when specifying waves.
	  :format:		asc		
	netbndbzsbzifile = sfincs_netbndbzsbzifile.nc
	  :description:		To specify all bnd, bzs (and bzi) input in 1 FEWS compatible netcdf input file. Specify either the netcdf version or ascii, not both.
	  :units:		m in projected UTM zone, m above reference level & m around mean water level of bzsfile
	  :required:		Only when specifying water levels and waves using netcdf input file.
	  :format:		net	 
	  
Forcing - Discharges
-----

	srcfile = sfincs.src
	  :description:		To specify discharge points, first the input locations have to be specified in 'sfincs.src'.
	  :units:		m in projected UTM zone
	  :required:		Only when specifying discharges.
	  :format:		asc	 
	disfile = sfincs.dis
	  :description:		In the file 'sfincs.dis' the discharge time-series are specified per input location. 
	  :units:		m^3/s
	  :required:		Only when specifying discharges.
	  :format:		asc	 	
	netsrcdisfile = sfincs_netsrcdisfile.nc
	  :description:		To specify all src & dis input in 1 FEWS compatible netcdf input file. Specify either the netcdf version or ascii, not both.
	  :units:		m in projected UTM zone, discharge in m3/s
	  :required:		Only when specifying discharges.
	  :format:		net	 

Forcing - Meteo
-----

	spwfile = sfincs.spw
	  :description:		Spiderweb file including wind speed, direction, pressure (and possibly rainfall).
	  :units:		coordinates: m in projected UTM zone, data: m/s, wind_from_direction in degrees, p_drop in Pa (and precipitation in mm/hr).
	  :required:		no
	  :format:		asc	 
	amufile = sfincs.amu
	  :description:		Delft3D-meteo ascii type input of wind speed in x-direction.
	  :units:		coordinates: m in projected UTM zone, data: m/s
	  :required:		no
	  :format:		asc	 	
	amvfile = sfincs.amv
	  :description:		Delft3D-meteo ascii type input of wind speed in y-direction.
	  :units:		coordinates: m in projected UTM zone, data: m/s
	  :required:		no
	  :format:		asc	 	  
	ampfile = sfincs.amp
	  :description:		Delft3D-meteo ascii type input of atmospheric pressure.
	  :units:		coordinates: m in projected UTM zone, data: Pa
	  :required:		no
	  :format:		asc
	amprfile = sfincs.ampr
	  :description:		Delft3D-meteo ascii type input of precipitation intensity.
	  :units:		coordinates: m in projected UTM zone, data: mm/hr
	  :required:		no
	  :format:		asc	 
	wndfile = sfincs.wnd
	  :description:		Spatially uniform wind 
	  :units:		wind speed in m/s, wind direction in nautical from where the wind is coming
	  :required:		no
	  :format:		asc	 	 
	precipfile = sfincs.prcp
	  :description:		Spatially uniform precipitation
	  :units:		mm/hr
	  :required:		no
	  :format:		asc	
	netamuamvfile = sfincs_netamuamvfile.nc
	  :description:		FEWS type netcdf meteo input with wind speed in both x-&y-direction in m/s.
	  :units:		coordinates: m in projected UTM zone, data: m/s
	  :required:		no
	  :format:		net	 	
	netampfile = sfincs_netampfile.nc
	  :description:		FEWS type netcdf meteo input with atmospheric pressure in Pa.
	  :units:		coordinates: m in projected UTM zone, data: Pa
	  :required:		no
	  :format:		net	 		
	netamprfile = sfincs_netamprfile.nc
	  :description:		FEWS type netcdf meteo input with precipitation in mm/hr.
	  :units:		coordinates: m in projected UTM zone, data: mm/hr
	  :required:		no
	  :format:		net	 			    
	  
Structures
-----

	thdfile = sfincs.thd
	  :description:		With a thin dam flow through certain grid cells is completely blocked (i.e. an infinitely high wall).
	  :units:		coordinates: m in projected UTM zone.
	  :required:		no
	  :format:		asc	 
	weirfile = sfincs.weir
	  :description:		Weirs are in principle the same as a thin dam, but then with a certain height (levee).
	  :units:		coordinates: m in projected UTM zone, elevation in m above reference level, weir formula coefficient in [-]
	  :required:		no
	  :format:		asc	 
	drnfile = sfincs.drn
	  :description:		Drainage pumps and culverts are both specified using the same format file, put with a different indication of the type (type=1 is drainage pump, type=2 is culvert).
	  :units:		coordinates: m in projected UTM zone, discharges in m^3/s.
	  :required:		no
	  :format:		asc	 
	  