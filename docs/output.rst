Output messages
=====

Interpreting the information on the screen
-----

When running SFINCS, some information will be written to the screen, and the same information is written to the log file 'sfincs.log' in the model folder.
This consists of what version was run, the input files that were read, the number of active points, and an overview of which processes are turned on.
After '---------- Starting simulation -----------' appears the initialisation phase is finished.

Once '0% complete,       - s remaining ...' appears, the model has actually started computing.
The following times the percentual progress % is shown (every 'percentage_done' %, default 5%), with a rough estimate of the time remaining until completion of the model run.

Once '---------- Simulation finished -----------' appears, your model has run successfully and the model output files are written.
Additionally some information is written to the screen regarding the total runtime, the time spent per section (input, boundaries, momentum, continuity, output, etc.) and the average time step.
A very small average time step compared to what you expect can indicate that the model has encountered instabilities somewhere.
With 'writeruntime = 1' the runtime per section is also written to a file 'runtimes.txt' (from SFINCS v2026.02 Hautacam release onwards).
Hereafter SFINCS is closed off, ready to start a new simulation.

Normally SFINCS runs using OpenMP, utilising all computing power of available cores (e.g. all 4).
This does not make it efficient to runs multiple SFINCS runs in parallel, it is advised to just run multiple simulations in series.

**Example**

.. code-block:: text


	 ----------- Welcome to SFINCS -----------
	
	  @@@@@  @@@@@@@ @@ @@  @@   @@@@   @@@@@ 
	 @@@ @@@ @@@@@@@ @@ @@@ @@ @@@@@@@ @@@ @@@
	 @@@     @@      @@ @@@ @@ @@   @@ @@@    
	  @@@@@  @@@@@@  @@ @@@@@@ @@       @@@@@ 
	     @@@ @@      @@ @@ @@@ @@   @@     @@@
	 @@@ @@@ @@      @@ @@  @@  @@@@@@ @@@ @@@
	  @@@@@  @@      @@ @@   @   @@@@   @@@@@ 

	              ..............              
	          ......:@@@@@@@@:......          
	       ..::::..@@........@@.:::::..       
	     ..:::::..@@..::..::..@@.::::::..     
    	.::::::..@@............@@.:::::::.    
	   .::::::..@@..............@@.:::::::.   
	  .::::::::..@@............@@..::::::::.  
	 .:::::::::...@@.@..@@..@.@@..::::::::::. 
	 .:::::::::...:@@@..@@..@@@:..:::::::::.. 
	 ............@@.@@..@@..@@.@@............ 
	 ^^^~~^^~~^^@@..............@@^^^~^^^~~^^ 
	 .::::::::::@@..............@@.:::::::::. 
	  .......:.@@.....@.....@....@@.:.......  
	   .::....@@......@.@@@.@....@@.....::.   
    	.:::~@@.:...:.@@...@@.:.:.@@~::::.    
	     .::~@@@@@@@@@@.....@@@@@@@@@~::.     
	       ..:~~~~~~~:.......:~~~~~~~:..      
	          ......................          
	              ..............              
	 -----------------------------------------

	By running a SFINCS pre-compiled executable
	provided by Deltares, you accept the
	conditions of the Deltares SFINCS FreeWare
	Software License Agreement, see:
	https://download.deltares.nl/sfincs

 	-----------------------------------------

 	Build-Revision: $Rev: v2.5.0 Hautacam
 	Build-Date:     $Date: 2026-10-08

	------ Preparing model simulation --------

	------------------------------------------
	Computational mesh
	------------------------------------------
	Number of active z points      :      4055
	Number of active u/v points    :      7972

	------------------------------------------
	Processes
	------------------------------------------
	Subgrid topography   : yes
	Quadtree refinement  : no
	Advection            : yes
	Viscosity            : yes
	Coriolis             : yes
	Wind                 : no
	Atmospheric pressure : no
	Precipitation        : yes
	Infiltration         : no
	SnapWave             : no
	Wave paddles         : no
	Non-hydrostatic      : no
	------------------------------------------

	---------- Starting simulation -----------

	   0% complete,       - s remaining ...
	   5% complete,     0.9 s remaining ...
	  10% complete,     0.9 s remaining ...
	  15% complete,     0.9 s remaining ...
	  20% complete,     0.9 s remaining ...
	  25% complete,     0.8 s remaining ...
	  30% complete,     0.8 s remaining ...
	  35% complete,     0.7 s remaining ...
	  40% complete,     0.7 s remaining ...
	  45% complete,     0.6 s remaining ...
	  50% complete,     0.6 s remaining ...
	  55% complete,     0.5 s remaining ...
	  60% complete,     0.5 s remaining ...
	  65% complete,     0.4 s remaining ...
	  70% complete,     0.4 s remaining ...
	  75% complete,     0.3 s remaining ...
	  80% complete,     0.2 s remaining ...
	  85% complete,     0.2 s remaining ...
	  90% complete,     0.1 s remaining ...
	  95% complete,     0.1 s remaining ...
	 100% complete,     0.0 s remaining ...

	Info : Write maximum values at final timestep since t=dtmaxout was not reached yet...

	---------- Simulation finished -----------

	 Total time                   :      1.211
	 Time in input                :      0.013
	 Time in boundaries           :      0.042 (  3.5%)
	 Time in momentum             :      0.881 ( 73.5%)
	 Time in continuity           :      0.207 ( 17.3%)
	 Time in output               :      0.055 (  4.6%)

	 Average time step (s)  :     22.031

	----------- Closing off SFINCS -----------


Possible error messages and possible solutions
-----

In case the following message is written to the screen, it means that something in the simulation has gone wrong.

.. code-block:: text

	Error! Minimum time step of 0.0050 s reached ! Current velocity exceeded uvmax 1000.0 m/s. Simulation stopped.

This means that the time step has become smaller than the minimum time step 'dtmin', which SFINCS computes from the maximum allowed flow velocity 'uvmax' (default 1000 m/s): dtmin = alpha * dx_min / (1.25 * uvmax).
In other words, a too large flow velocity has occurred somewhere in the domain, indicating that the model has become unstable and some input is probably not optimal.
Before stopping, SFINCS writes the map, his and max output of this last time step, so you can look up where the instability occurred (e.g. by plotting 'zs', 'h' or 'zsmax').

**NOTE - Before SFINCS v2.1.1 Dollerup, the simulation was stopped on a maximum water depth ('stopdepth', message 'Maximum depth of 100.0 m reached!!! Simulation stopped.'). This keyword is no longer used.**

Possible problems can be:

* The provided elevation file has very rapid changes in elevation, that locally lead to large water level gradients and fluxes. Possible solution: locally smooth the elevation data and provide this as a new depfile.

* In general the internal timesteps of SFINCS might be too large. Possible solution: reduce timesteps by supplying a lower value of alpha (e.g. 0.4, default is 0.5), set a higher value for 'hmin_cfl' or set a low enough value of 'dtmax'.

* To analyze which locations are limiting the global timestep, enable `timestep_analysis = 1` in `sfincs.inp` to write variables `average_required_timestep` and `percentage_limiting_timestep` to the sfincs_map.nc file. See the Timestep analysis section in the Model settings for more details.

* When only forcing discharges in a for the rest entirely dry domain, the initial time steps can be too coarse to account for the needed timesteps when the discharge starts to flow. Possible solution: Make sure that part of the river/domain initially has water (limiting the time step) by specifying either 'zsini' or an 'inifile'.

* When forcing waves, the bzifile time-series might contain too rapid changes in water level, the internal timesteps of SFINCS are too large. Possible solution: reduce timesteps by supplying a lower value of alpha (e.g. 0.4, default is 0.5).

* **Tip to check your model**: specify netcdf output and load in the sfincs_map.nc file (e.g. QGIS, Quickplot, Panoply, Matlab, Python) and have a look at the variables 'zb' and 'msk'. Then you can see how SFINCS has interpreted the provided elevation and mask input. Do map plots of these variables look weird? Probably something in your input file is not entirely correct!

* When more stability is needed still, have a look at the input parameter options of 'advlim' or 'hmin_cfl'.

Besides model instabilities, other recurring problems might be:

* A specified (forcing) file/parameters is not read in > check whether you specified the name (e.g. netamuamvfile   = netamuamv.nc ) with **ONLY SPACES** in between the keyword and argument. SFINCS does not interpret a mixture of spaces and tabs well. This may cause a file or parameter to be read in as 'none', whereafter this is not used in the model simulation as wanted.

* Also, check whether a certain expected forcing is coming through. SFINCS writes an overview of the 'Processes' to the screen and to 'sfincs.log' (see the example above), so if you force rainfall but 'Precipitation' says 'no', something probably went wrong with the input file. The same holds for 'Subgrid topography', 'Quadtree refinement', 'Advection', 'Viscosity', 'Coriolis', 'Wind', 'Atmospheric pressure', 'Infiltration', 'Wave paddles' and 'Non-hydrostatic'.

* SFINCS also gives input about certain files after reading in data and how these are snapped/interpreted. For instance for subgrid; "Number of subgrid levels : XXX", weirfile; "XXX structure u/v points found", wavemaker; "Number of wavemaker polylines found : XXX", observation points; "Warning : observation point XXX falls outside model domain.' - compare whether this is as expected.

Tips and tricks
-----

* **Tip to view netcdf output while SFINCS is still running**: by default the HDF5 library (used for writing sfincs_map.nc and sfincs_his.nc) locks the output files while SFINCS writes to them, so other programs (e.g. Quickplot, Panoply, Matlab, Python/xarray) cannot open them until the simulation has finished. Disabling HDF5 file locking, by setting the environment variable HDF5_USE_FILE_LOCKING to FALSE before starting SFINCS, allows you to inspect the output during the simulation:

.. code-block:: text

	Windows (command prompt or run.bat):
	set HDF5_USE_FILE_LOCKING=FALSE
	sfincs.exe

	Linux:
	export HDF5_USE_FILE_LOCKING=FALSE
	sfincs

	Docker:
	docker run -e HDF5_USE_FILE_LOCKING=FALSE -v <model_folder>:/data deltares/sfincs-cpu

Note that the file being read may not yet contain the latest output timestep, and that the reading program might also need HDF5_USE_FILE_LOCKING=FALSE in its own environment.


Output description
=====

Fill values: dry cells and undefined quantities
-----

All variables in sfincs_map.nc and sfincs_his.nc have a fill value of -99999.0, set as the '_FillValue' attribute, so that most tools (e.g. xarray, QGIS, Matlab ncread with fill value handling) read these as NaN.
From **SFINCS v2.5.0 (Hautacam release) onwards**, the map, maximum and his output use one rule: the fill value is written where a cell or observation point is dry, where there are no waves, or where a quantity is undefined.
Cells that are inactive (msk = 0) or fall outside the quadtree/regular output grid are fill value as well.

A cell or observation point is considered dry when its water depth is at or below 'huthresh' (default 0.05 m), with the water depth computed as zs - zb (regular) or zs - subgrid_z_zmin (subgrid).

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Variables
     - Fill value where
   * - zs, h, point_zs, point_h
     - the cell or point is dry
   * - zsmax, hmax
     - the cell has not been wet during the simulation (maximum water depth at or below 'huthresh')
   * - u, v, point_u, point_v, point_uvmag
     - the cell or point is dry
   * - point_uvdir
     - the point is dry, or there is no flow (direction undefined)
   * - vmax, qmax
     - the cell has not been wet during the simulation (only applied when 'zsmax' is stored as well)
   * - t_zsmax
     - no maximum water level time was recorded in the cell (t_zsmax not larger than 0)
   * - point_wind_dir
     - there is no wind (direction undefined)
   * - hm0, hm0ig (map and his)
     - there are no waves (H = 0, e.g. on dry SnapWave nodes)
   * - tp, wavdir, cg, wave forces, beta, dw, df, srcig, alphaig
     - there are no incident waves (H = 0), as these are undefined there
   * - tpig, dwig, dfig
     - there are no infragravity waves (H_ig = 0)
   * - snapwavedepth
     - the SnapWave depth is at or below the SnapWave 'hmin'

**Note for post-processing**: dry cells therefore do not contain a water depth or velocity of 0, but NaN. Use NaN-aware functions (e.g. numpy.nanmax, xarray's default skipna) when computing statistics, and fill the NaNs with 0 yourself where a zero depth or velocity is needed (e.g. ds['h'].fillna(0.0)).

Parameters netcdf file global (sfincs_map.nc)
-----

In case of netcdf output, the given parameters mean the following.
For a quadtree model (see :ref:`quadtree`), the cell centre coordinates are 'mesh2d_face_x' and 'mesh2d_face_y' instead of 'x' and 'y', and all cell variables have the dimension 'mesh2d_nFaces' (UGRID conventions).

	x
	  :description:		x coordinate of cell centers in projected reference system
	  :standard_name:	projection_x_coordinate
	  :units:		m in projected reference system	  
	y
	  :description:		y coordinate of cell centers in projected reference system
	  :standard_name:	projection_y_coordinate	  
	  :units:		m in projected reference system  
	zb
	  :description:		Bed level elevation (in case of subgrid version of SFINCS, this elevation is not used in the model but the sbgfile with subgrid tables is used instead).
	  :standard_name:	altitude	  
	  :units:		m above reference level
	msk
	  :description:		Mask of active cells as interpreted by SFINCS: inactive (0), active (1), water level boundary (2) or outflow boundary (3) cell.
	  :standard_name:	land_binary_mask
	  :units:		-
	time
	  :description:		Time of global map output.
	  :standard_name:	time	  
	  :units:		seconds since 'tref'	  
	zs
	  :description:		Instantaneous water level per 'dtmapout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_surface_height_above_mean_sea_level	  
	  :units:		m above reference level
	h
	  :description:		Instantaneous water depth per 'dtmapout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	depth	  
	  :units:		m
	u
	  :description:		Instantaneous flow velocity in u-direction per 'dtmapout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_water_x_velocity	  
	  :units:		m/s
	v
	  :description:		Instantaneous flow velocity in v-direction per 'dtmapout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_water_y_velocity	  
	  :units:		m/s	
	subgrid_volume
	  :description:		Instantaneous subgrid volume per 'dtmapout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	subgrid_volume_in_cell	  
	  :units:		m^3
	storage_volume
	  :description:		Instantaneous storage volume per 'dtmapout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	storage_volume_in_cell	  
	  :units:		m^3	  
	timemax
	  :description:		Time of global map output per 'dtmaxout' timestep.
	  :standard_name:	time	  
	  :units:		seconds since 'tref'	  
	zsmax
	  :description:		Maximum water level per 'dtmaxout' timestep, only given if dtmaxout>0, corresponding with netcdf variable 'timemax'.
	  :standard_name:	maximum of sea_surface_height_above_mean_sea_level
	  :units:		m above reference level
	hmax
	  :description:		Maximum water depth per 'dtmaxout' timestep, only given if dtmaxout>0, corresponding with netcdf variable 'timemax'. In subgrid mode only given with 'storehsubgrid = 1', see 'Input parameters'.
	  :standard_name:	sea_floor_depth_below_sea_surface
	  :units:		m
	zvolmax
	  :description:		Maximum subgrid volume in cell per 'dtmaxout' timestep, only given if dtmaxout>0, subgrid is used and storezvolmax = 1 (default), corresponding with netcdf variable 'timemax'.
	  :standard_name:	maximum_subgrid_volume_in_cell
	  :units:		m^3
	t_zsmax
	  :description:		Time of max water level per cell and per 'dtmaxout' timestep, only given if dtmaxout>0 and 'storetzsmax = 1', corresponding with netcdf variable 'timemax'.
	  :standard_name:	t_zsmax
	  :units:		seconds since 'tref'
	tmax
	  :description:		Duration that a cell was wet (water depth above 'twet_threshold') per 'dtmaxout' timestep, only given if 'storetwet = 1', corresponding with netcdf variable 'timemax'.
	  :standard_name:	duration_wet
	  :units:		s
	vmax
	  :description:		Maximum flow velocity proxy per 'dtmaxout' timestep, only given if dtmaxout>0 and 'storevelmax = 1', corresponding with netcdf variable 'timemax'.
	  :standard_name:	maximum_flow_velocity
	  :units:		m/s
	qmax
	  :description:		Maximum flow flux proxy per 'dtmaxout' timestep, only given if dtmaxout>0 and 'storefluxmax = 1', corresponding with netcdf variable 'timemax'.
	  :standard_name:	maximum_flux	  
	  :units:		m^2/s	  	  	  
	cuminf
	  :description:		Cumulative infiltration depth since the start of the simulation, only given if 'storecumprcp = 1' and infiltration is used, corresponding with netcdf variable 'timemax'.
	  :units:		m
	cumprcp
	  :description:		Cumulative precipitation depth since the start of the simulation, only given if 'storecumprcp = 1', corresponding with netcdf variable 'timemax'.
	  :units:		m
	urban_drainage_cumulative_depth
	  :description:		Cumulative urban drainage depth (drained volume / cell area) since the start of the simulation, only given if 'store_cumulative_urban_drainage = 1' and 'urbfile' is specified, corresponding with netcdf variable 'timemax'.
	  :units:		m
	wind_u, wind_v
	  :description:		Wind speed components per 'dtmapout' timestep, only given if 'storemeteo = 1' and wind input is given.
	  :units:		m/s
	windmax
	  :description:		Maximum wind speed per 'dtmaxout' timestep, only given if 'storemeteo = 1', 'storemaxwind = 1' and spatially varying wind input is given.
	  :units:		m/s
	surface_air_pressure
	  :description:		Atmospheric pressure per 'dtmapout' timestep, only given if 'storemeteo = 1' and pressure input is given.
	  :units:		Pa
	precipitation_rate
	  :description:		Precipitation rate per 'dtmapout' timestep, only given if 'storemeteo = 1' and precipitation input is given.
	  :units:		mm/hr
	bucket_volume
	  :description:		Current bucket storage per 'dtmapout' timestep, only given for inftype = bkt.
	  :units:		m
	inp
	  :description:		Copy of all the supplied input to SFINCS from 'sfincs.inp'.
	  :units:		-
	total_runtime
	  :description:		Total model runtime in seconds, as displayed by SFINCS to the screen.
	  :units:		s	
	average_dt
	  :description:		Model average timestep in seconds, as displayed by SFINCS to the screen.
	  :units:		s

Parameters netcdf file observation points (sfincs_his.nc)
-----	

This file is only created if at least one of the following is supplied: observation points ('obsfile'), cross-sections, weirs, drainage structures, run-up gauges, river discharge points combined with 'store_river_discharge = 1', or urban drainage zones combined with 'store_urban_drainage_discharge = 1'.

	point_x
	  :description:		x coordinate of interpreted observation points in projected reference system
	  :standard_name:	projection_x_coordinate
	  :units:		m in projected reference system	  
	point_y
	  :description:		y coordinate of interpreted observation points in projected reference system
	  :standard_name:	projection_y_coordinate	  
	  :units:		m in projected reference system
	station_x
	  :description:		x coordinate of specified observation points in projected reference system
	  :standard_name:	projection_x_coordinate
	  :units:		m in projected reference system	  
	station_y
	  :description:		y coordinate of specified observation points in projected reference system
	  :standard_name:	projection_y_coordinate	  
	  :units:		m in projected reference system	  
	structure_x
	  :description:		x coordinate of snapped location on SFINCS grid of weirs in projected reference system
	  :standard_name:	projection_x_coordinate	  
	  :units:		m in projected reference system	 	  
	structure_y
	  :description:		y coordinate of snapped location on SFINCS grid of weirs in projected reference system
	  :standard_name:	projection_y_coordinate	  
	  :units:		m in projected reference system	 
	structure_height
	  :description:		weir height on snapped location on SFINCS grid of weirs
	  :standard_name:	structure_height
	  :units:		m above reference level 	  	  
	thindam_x
	  :description:		x coordinate of snapped location on SFINCS grid of thin dams in projected reference system
	  :standard_name:	projection_x_coordinate	  
	  :units:		m in projected reference system	 	  
	thindam_y
	  :description:		y coordinate of snapped location on SFINCS grid of thin dams in projected reference system
	  :standard_name:	projection_y_coordinate	  
	  :units:		m in projected reference system	 	  
	point_zb
	  :description:		Bed level elevation of observation points.
	  :standard_name:	altitude	  
	  :units:		m above reference level
	time
	  :description:		Time of his output.
	  :standard_name:	time	  
	  :units:		seconds since 'tref'	  
	point_zs
	  :description:		Instantaneous water level per 'dthisout' timestep of observation points, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_surface_height_above_mean_sea_level	  
	  :units:		m above reference level
	point_h
	  :description:		Instantaneous water depth per 'dthisout' timestep of observation points, corresponding with netcdf variable 'time'.
	  :standard_name:	point_h	  
	  :units:		m
	point_u
	  :description:		Instantaneous flow velocity in u-direction per 'dthisout' timestep of observation points, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_water_x_velocity
	  :units:		m/s	
	point_v
	  :description:		Instantaneous flow velocity in v-direction per 'dthisout' timestep of observation points, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_water_y_velocity
	  :units:		m/s	  	
	point_uvmag
	  :description:		Instantaneous absolute flow velocity per 'dthisout' timestep of observation points, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_water_velocity
	  :units:		m/s	  	
	point_uvdir
	  :description:		Instantaneous flow direction (bearing) per 'dthisout' timestep of observation points, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_water_velocity_direction
	  :units:		degrees wrt north	 	        
	point_prcp
	  :description:		Instantaneous precipitation rate per 'dthisout' timestep, corresponding with netcdf variable 'time'.
	  :units:		mm/hr
	point_cumprcp
	  :description:		Cumulative precipitation depth per 'dthisout' timestep, only given if 'storecumprcp = 1', corresponding with netcdf variable 'time'.
	  :units:		m
	point_qinf
	  :description:		Instantaneous infiltration rate per 'dthisout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	point_qinf
	  :units:		mm/hr
	point_S
	  :description:		Current infiltration state per 'dthisout' timestep: Se (m) for inftype = cnb, soil moisture deficit for gai, infiltration capacity (mm/hr) for hor, bucket storage (m) for bkt.
	  :units:		depends on inftype
	crosssection_discharge
	  :description:		Discharge through cross-section per 'dthisout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	discharge	  
	  :units:		m3/s
	drainage_discharge
	  :description:		Discharge through drainage structure per 'dthisout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	discharge	  
	  :units:		m3/s
	river_name
	  :description:		Name of each river discharge point, taken from the optional 3rd column of the 'srcfile', or auto-generated as 'discharge_0001', 'discharge_0002', etc. when no name is given (always the case for 'netsrcdisfile'). Only written when 'store_river_discharge = 1'.
	  :units:		-
	river_discharge
	  :description:		Discharge per river discharge point per 'dthisout' timestep, corresponding with netcdf variable 'time'. This is the prescribed discharge from the dis/netsrcdis input, linearly interpolated to the model time (instantaneous value, not time-averaged). It is not corrected afterwards, so a point that falls outside the active model domain is still reported, even though its discharge does not enter the model. Only written when 'store_river_discharge = 1'.
	  :standard_name:	discharge
	  :units:		m3/s
