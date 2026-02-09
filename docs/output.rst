Output messages
=====

Interpreting the information on the screen
-----

When running SFINCS, some information will be written to the screen.
This consists of what version was run (number/date not always completely up to date), some files that were read in and the number of active points.
After 'Starting computation ...' appears the initialisation phase is finished.

Once '0% complete,     Inf s remaining ...' appears, the model has actually started computing.
The following times the percentual progress % is shown, a rough estimate of the time remaining untill completion of the model run is given.

Once '---Simulation is finished---', your model has run succesfully and is writing away the model output files.
Additionaly some information is written to the screen regarding total runtime, time consumption per section, the average time step, and the maximum occured water depth in the entire computation.
If you know the initial water depth, this can give an indication whether the model has encountered instabilities or not.
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

 	-----------------------------------------

 	Build-Revision: $Rev: v0.0.1-alpha$
 	Build-Date:     $Date: 2022-10-25$

	Reading input file ...
 	Info : Running SFINCS in subgrid mode ...
 	Reading meteo data ...
 	Info : Preparing SFINCS grid on regular mesh ...
 	Reading sfincs.ind ...
 	Reading sfincs.msk ...
 	Number of active z points    :         4055
 	Number of active u/v points  :         7972
 	Reading sfincs.sbg ...
 	Reading water level boundaries ...
 	Reading observation points ...
 	Initializing output ...

	 ---------- Starting simulation ----------

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

	 Info : Write maximum values of final timestep since t=dtmaxout was not reached 
	 yet...

	 ---------- Simulation finished ----------
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

	 Total time             :      1.211
	 Total simulation time  :      1.198
	 Time in input          :      0.013
	 Time in boundaries     :      0.042 (  3.5%)
	 Time in momentum       :      0.881 ( 73.5%)
	 Time in continuity     :      0.207 ( 17.3%)
	 Time in output         :      0.055 (  4.6%)

	 Average time step (s)  :     22.031

	 ---------- Closing off SFINCS -----------


Possible error messages and possible solutions
-----

In case the following message is written to the screen, it means that something in the simulation has gone wrong.

.. code-block:: text

	Maximum depth of 100.0 m reached!!! Simulation stopped.

This means that a too large water depth has occured somewhere in the domain, indicating that some input is probably not optimal.

Possible problems can be:

* The provided elevation file has very rapid changes in elevation, that locally lead to large water level gradients and fluxes. Possible solution: locally smooth the elevation data and provide this as a new depfile.

* In general the internal timesteps of SFINCS might be too large. Possible solution: reduce timesteps by supplying a lower value of alpha (e.g. 0.5), set a higher value for 'hmin_cfl' or set a low enough value of 'dtmax'.

* When only forcing discharges in a for the rest entirely dry domain, the initial time steps can be too coarse to account for the needed timesteps when the discharge starts to flow. Possible solution: Make sure that part of the river/domain initially has water (limiting the time step) by specifying either 'zsini' or an 'inifile'.

* When forcing waves, the bzifile time-series might contain too rapid changes in water level, the internal timesteps of SFINCS are too large. Possible solution: reduce timesteps by supplying a lower value of alpha (e.g. 0.5).

* **Tip to check your model**: specify netcdf output and load in the sfincs_map.nc file (e.g. Quickplot, Panoply, Matlab, Python) and have a look at the variables 'zb' and 'msk'. Then you can see how SFINCS has interpreted the prodivided depfile and mskfile. Does map plots of these variables look weird? Probably something in your input file is not entirely correct!

* When more stability is needed still, have a look at the input parameter options of 'advlim' or 'hmin_cfl'.

Besides model instabilities, other recurring problems might be:

* A specified (forcing) file/parameters is not read in > check whether you specified the name (e.g. netamuamvfile   = netamuamv.nc ) with **ONLY SPACES** in between the keyword and argument. SFINCS does not interpret a mixture of spaces and tabs well. This may cause a file or parameter to be read in as 'none', whereafter this is not used in the model simulation as wanted.

* Also, check whether a certain expected forcing is coming through. SFINCS displays messages like "Turning on process: Precipitation", so if you force rainfall is this message is not visible in your log-file, something probably went wrong with the input file. Also for "Advection scheme", "Wind", "Atmospheric pressure", "Coriolis", "Viscosity", "Dynamic waves", "Infiltration XXX-type", "Precipitation from spwfile", "Storage Green Infrastructure".

* SFINCS also gives input about certain files after reading in data and how these are snapped/interpreted. For instance for subgrid; "Number of subgrid levels : XXX", weirfile; "XXX structure u/v points found", wavemaker; "Number of wavemaker polylines found : XXX", observation points; "Warning : observation point XXX falls outside model domain.' - compare whether this is as expected.


Output description
=====

Time step diagnostics (CSV)
-----

If you set 'timestep_diagnostics = 1' in **sfincs.inp**, SFINCS writes two optional diagnostic CSV files to the run directory:

* **timestep_diagnostics.csv**: time series of the UV-point that limits the global timestep. The write interval can be controlled with 'dt_timestep_diagnostics' (seconds); set to 0 to write every timestep.
* **timestep_diagnostics_domain.csv**: end-of-run domain summary for cells that were adjacent to a limiting UV-point at least once (includes how often they limited and the minimum limiting timestep).
* In **timestep_diagnostics.csv**, the `reason` column indicates limiter type: `0` = wave-speed limited (`sqrt(g*h)`), `1` = velocity limited (`|u|`), `2` = limited by `dtmax` (no UV-point limiter for that step).

Parameters netcdf file global (sfincs_map.nc)
-----

In case of netcdf output, the given parameters mean the following:

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
	  :description:		Time-step global map output.
	  :standard_name:	land_binary_mask	  
	  :units:		-
	time
	  :description:		Time of global map output.
	  :standard_name:	time	  
	  :units:		seconds since 'tref'	  
	zs
	  :description:		Instantaneous water level per 'dtout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_surface_height_above_mean_sea_level	  
	  :units:		m above reference level
	h
	  :description:		Instantaneous water depth per 'dtout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	depth	  
	  :units:		m
	u
	  :description:		Instantaneous flow velocity in u-direction per 'dtout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_water_x_velocity	  
	  :units:		m/s
	v
	  :description:		Instantaneous flow velocity in v-direction per 'dtout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_water_y_velocity	  
	  :units:		m/s	
	subgrid_volume
	  :description:		Instantaneous subgrid volume per 'dtout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	subgrid_volume_in_cell	  
	  :units:		m^3
	storage_volume
	  :description:		Instantaneous storage volume per 'dtout' timestep, corresponding with netcdf variable 'time'.
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
	t_zsmax
	  :description:		Time of max water level per cell and per 'dtmaxout' timestep, only given if dtmaxout>0, corresponding with netcdf variable 'timemax'.
	  :standard_name:	maximum of sea_surface_height_above_mean_sea_level	  
	  :units:		m above reference level	  
	vmax
	  :description:		Maximum flow velocity proxy per 'dtmaxout' timestep, only given if dtmaxout>0, corresponding with netcdf variable 'timemax'.
	  :standard_name:	maximum_flow_velocity	  
	  :units:		m/
	qmax
	  :description:		Maximum flow flux proxy per 'dtmaxout' timestep, only given if dtmaxout>0, corresponding with netcdf variable 'timemax'.
	  :standard_name:	maximum_flux	  
	  :units:		m^2/s	  	  	  
	cuminf
	  :description:		Cumulative infiltration depth over whole simulation.
	  :units:		m	  
	cumprcp
	  :description:		Cumulative precipitation depth over whole simulation.
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

This file is only created if observation points are supplied in the 'obsfile', or if weirs/cross-sections are supplied.

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
	  :description:		weir height on snapped location on SFINCS grid of weirs in projected reference system
	  :standard_name:	projection_x_coordinate	  
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
	  :description:		Instantaneous absolute flow velocity per 'dthisout' timestep of observation points, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_water_velocity_direction
	  :units:		degrees wrt north	 	        
	point_prcp
	  :description:		Instantaneous precipitation rate 'dthisout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_surface_height_above_mean_sea_level	  
	  :units:		m above reference level
	point_qinf
	  :description:		Instantaneous infiltration rate per 'dthisout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	point_qinf	  
	  :units:		m
	point_hm0
	  :description:		Significant short-wave height per 'dthisout' timestep (only when SnapWave is enabled).
	  :standard_name:	hm0_wave_height
	  :units:		m
	point_hm0ig
	  :description:		Significant infragravity-wave height per 'dthisout' timestep (only when SnapWave is enabled).
	  :standard_name:	hm0_ig_wave_height
	  :units:		m
	point_tp
	  :description:		Peak wave period per 'dthisout' timestep (only when SnapWave is enabled).
	  :standard_name:	peak_wave_period
	  :units:		s
	point_tpig
	  :description:		Peak infragravity-wave period per 'dthisout' timestep (only when SnapWave is enabled).
	  :standard_name:	ig_peak_wave_period
	  :units:		s
	point_wavdir
	  :description:		Mean wave direction per 'dthisout' timestep (only when SnapWave is enabled and wave direction output is stored).
	  :standard_name:	mean_wave_direction
	  :units:		degrees
	point_dirspr
	  :description:		Wave directional spreading per 'dthisout' timestep (only when SnapWave is enabled and wave direction output is stored).
	  :standard_name:	wave_directional_spreading
	  :units:		degrees
	point_zsm
	  :description:		Filtered water level per 'dthisout' timestep (only when SnapWave is enabled and wavemaker is used).
	  :standard_name:	filtered_water_level
	  :units:		m
	point_dw, point_df, point_dwig, point_dfig, point_cg, point_beta, point_srcig, point_alphaig
	  :description:		Optional SnapWave force-related diagnostics at observation points (only when SnapWave is enabled and wave-force output is stored).
	  :units:		see variable attributes in sfincs_his.nc
	snapwave_his_legacy_names
	  :description:		If enabled in sfincs.inp, SFINCS also writes legacy SnapWave his variable names (hm0, hm0ig, tp, tpig, wavdir, dirspr, zsm, dw, df, dwig, dfig, cg, beta, srcig, alphaig) in addition to the point_ names for backward compatibility.
	crosssection_discharge
	  :description:		Discharge through cross-section per 'dthisout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	discharge	  
	  :units:		m3/s
	drainage_discharge
	  :description:		Discharge through drainage structure per 'dthisout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	discharge	  
	  :units:		m3/s
