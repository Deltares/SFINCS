Output messages
======================

Interpreting the information on the screen
----------------------

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

	 ---Welcome to SFINCS---

	 Build-Revision: $Rev: 220 $
	 Build-Date: $Date: 2021-01-15$

	 Reading input file ...
	 Reading sfincs.ind
	 Reading sfincs.dep
	 Reading sfincs.msk
	      330418  active points
	 Determining m and n indices ...
	 Finding neighboring cell indices ...
	 Reading sfincs.man
	 Reading sfincs.scs
	 Reading water level boundaries ...
	 Reading observation points ...
	 Starting computation ...
	   0% complete,     Inf s remaining ...
	   5% complete,    42.9 s remaining ...
	  10% complete,    40.1 s remaining ...
	  15% complete,    36.1 s remaining ...
	  20% complete,    33.3 s remaining ...
	  25% complete,    31.1 s remaining ...
	  30% complete,    29.4 s remaining ...
	  35% complete,    27.3 s remaining ...
	  40% complete,    25.4 s remaining ...
	  45% complete,    23.3 s remaining ...
	  50% complete,    21.2 s remaining ...
	  55% complete,    19.1 s remaining ...
	  60% complete,    16.8 s remaining ...
	  65% complete,    14.7 s remaining ...
	  70% complete,    12.5 s remaining ...
	  75% complete,    10.4 s remaining ...
	  80% complete,     8.2 s remaining ...
	  85% complete,     6.1 s remaining ...
	  90% complete,     4.1 s remaining ...
	  95% complete,     2.0 s remaining ...
	 100% complete,     0.0 s remaining ...
	 ---Simulation is finished---

	 Total run time          :     40.624
	 Time in input           :      0.118
	 Time in boundary loop   :      1.445 (  3.6%)
	 Time in discharges      :      0.004 (  0.0%)
	 Time in meteo (1)       :      0.016 (  0.0%)
	 Time in meteo (2)       :      5.293 ( 13.0%)
	 Time in flux loop       :     22.789 ( 56.1%)
	 Time in structures loop :      0.000 (  0.0%)
	 Time in continuity loop :      6.351 ( 15.6%)
	 Time in output          :      4.682 ( 11.5%)

	 Average time step (s)   :      3.329
	 Maximum depth (m)       :     17.391

	 ---Closing off SFINCS---


Possible error messages and possible solutions
----------------------

In case the following message is written to the screen, it means that something in the simulation has gone wrong.

.. code-block:: text

	Maximum depth of 100.0 m reached!!! Simulation stopped.
	Maximum depth occurs at (n,m)=(  705,  690), (x,y)=(  511087.5,14252488.0).

This means that a too large water depth has occured somewhere in the domain, indicating that some input is probably not optimal.
As bonus, the grid cell indices and x&y location is given for faster debugging.

Possible problems can be:

- The provided elevation file has very rapid changes in elevation, that locally lead to large water level gradients and fluxes. Possible solution: locally smooth the elevation data and provide this as a new depfile

- In general the internal timesteps of SFINCS might be too large. Possible solution: reduce timesteps by supplying a lower value of alpha (e.g. 0.5) or set a low enough value of 'dtmax'.

- Sometimes a simulation might contain too large water depths are start in too deep water. This can potentially create problems as SFINCS is intented as a shallow water model.

- When only forcing discharges in a for the rest entirely dry domain, the initial time steps can be too coarse to account for the needed timesteps when the discharge starts to flow. Possible solution: Make sure that part of the river/domain initially has water (limiting the time step) by specifying either 'zsini' or an 'inifile'.

- When forcing waves, the bzifile time-series might contain too rapid changes in water level, the internal timesteps of SFINCS are too large. Possible solution: reduce timesteps by supplying a lower value of alpha (e.g. 0.5).

- **Tip to check your model**: specify netcdf output and load in the sfincs_map.nc file (e.g. Quickplot, Panoply, Matlab, Python) and have a look at the variables 'zb' and 'msk'. Then you can see how SFINCS has interpreted the prodivided depfile and mskfile. Does map plots of these variables look weird? Probably something in your input file is not entirely correct!


Besides model instabilities, other recurring problems might be:

- A specified (forcing) file/parameters is not read in > check whether you specified the name (e.g. netamuamvfile   = netamuamv.nc ) with **ONLY SPACES** in between the keyword and argument. SFINCS does not interpret a mixture of spaces and tabs well. This may cause a file or parameter to be read in as 'none', whereafter this is not used in the model simulation as wanted.


Output description
======================

Parameters netcdf file global (sfincs_map.nc)
----------------------

In case of netcdf output, the given parameters mean the following:

	x
	  :description:		x coordinate of cell centers in projected reference system
	  :standard_name:	projection_x_coordinate
	  :units:		m in projected reference system	  
	y
	  :description:		y coordinate of cell centers in projected reference system
	  :standard_name:	projection_y_coordinate	  
	  :units:		m in projected reference system
	edge_x
	  :description:		x coordinate of cell corners in projected reference system
	  :standard_name:	projection_x_coordinate
	  :units:		m in projected reference system	  
	edge_y
	  :description:		y coordinate of cell corners in projected reference system
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
	timemax
	  :description:		Time of global map output per 'dtmaxout' timestep.
	  :standard_name:	time	  
	  :units:		seconds since 'tref'	  
	zsmax
	  :description:		Maximum water level per 'dtmaxout' timestep, only given if dtmaxout>0, corresponding with netcdf variable 'timemax'.
	  :standard_name:	maximum of sea_surface_height_above_mean_sea_level	  
	  :units:		m above reference level
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
----------------------		

This file is only created if observation points are supplied in the 'obsfile'.

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
	  :standard_name:	depth	  
	  :units:		m
	point_prcp
	  :description:		Instantaneous precipitation rate 'dthisout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	sea_surface_height_above_mean_sea_level	  
	  :units:		m above reference level
	point_qinf
	  :description:		Instantaneous infiltration rate per 'dthisout' timestep, corresponding with netcdf variable 'time'.
	  :standard_name:	depth	  
	  :units:		m
