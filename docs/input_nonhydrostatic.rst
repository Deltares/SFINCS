User manual - non-hydrostatic and moving bed
============================================

Overview
--------

For short, dispersive waves (tsunamis from small sources, landslide-generated waves, solitary waves on beaches) the hydrostatic shallow-water equations are not accurate enough.
SFINCS has an optional non-hydrostatic pressure correction for these cases, and the bed level can be changed during a run through the BMI/XMI interface to model a moving seafloor or landslide.
All keywords mentioned below are listed in the tab 'Input parameters'.
These options are meant for advanced users; they are normally combined with 'momentum_scheme = velocity'.

Non-hydrostatic pressure
------------------------

Method
^^^^^^

The non-hydrostatic option ('nonhydrostatic = 1') adds a pressure that varies linearly over a single layer (one-layer Keller-box approach).
After each hydrostatic step a pressure equation is solved and the velocities and fluxes are corrected with its gradient.
The pressure operator is symmetric and positive definite.
It is assembled on a 5-point stencil and solved with a Jacobi-preconditioned conjugate-gradient method (tolerance 'nonh_tol', maximum 'nonh_itermax' iterations), warm-started from the previous time step.
The solver runs in parallel with OpenMP (CPU) and OpenACC (GPU), on regular and quadtree grids, for projected and geographic coordinates.
On quadtree grids, cells at refinement transitions stay hydrostatic, and an optional variable 'nonh_mask' in the quadtree netcdf file restricts the correction to part of the domain.
The non-hydrostatic pressure is written to the map file as 'pnonh'.

Bottom kinematic condition
^^^^^^^^^^^^^^^^^^^^^^^^^^

The bottom kinematic condition is carried implicitly by the layer term of the pressure operator (Escalante, Morales de Luna & Castro, 2018), so no separate bottom vertical velocity is needed.
A moving bed enters as the known source d(zb)/dt (see 'nonh_movingbed' below).

Dispersion correction
^^^^^^^^^^^^^^^^^^^^^

The strict one-layer Keller box has a phase speed close to the linear (Airy) value up to kd = 2.5, but it has a frequency cutoff: shorter waves have no group velocity and can stay trapped on a slope.
SFINCS therefore replaces the vertical momentum equation of the layer with an implicit Madsen-Sorensen-type form.
This gives the dispersion relation

.. math::

	\omega^2 = g h k^2 \frac{1 + \beta (kh)^2}{1 + (\beta + 1/3)(kh)^2}

without a cutoff.
With the default 'nonh_pade_beta = 0.04' the linear phase speed is within about 3 % of Airy up to kd = 6 (beta = 1/15 is the [2/2] Pade approximant, beta = 0 the [0/2] approximant).
The correction keeps the operator symmetric and positive definite and adds no time step restriction.

Breaking and filtering
^^^^^^^^^^^^^^^^^^^^^^

Breaking waves are treated by switching off the non-hydrostatic pressure locally:

* 'nonh_brfr' (default 0.5): Froude criterion as in NEOWAVE; off where abs(U)/sqrt(g D) > nonh_brfr, on again below 0.3 * nonh_brfr. 0 = no breaking.
* 'nonh_treform' (default 1 s) lets released cells recover the pressure gradually.

Further options: 'nonh_filter' (default 0.5) is a 2dx filter on the pressure after each solve, 'nonh_fnudge' (default 1) sets the fraction of the correction that is applied and 'nonh_tstop' switches the correction off after a given time.

.. code-block:: text

	momentum_scheme  = velocity
	nonhydrostatic   = 1
	nonh_pade_beta   = 0.04

Moving bed via BMI/XMI
----------------------

Applying a bed change
^^^^^^^^^^^^^^^^^^^^^

When SFINCS is run as a library through the BMI/XMI interface (e.g. with 'xmipy' or the 'SfincsXmi' class of 'cht_sfincs'), the bed level can be changed between time steps.
The caller writes a bed level increment per cell into the array 'dzbext' and calls 'update_bed_level'.
The kernel then:

* shifts the bed level (zb on regular grids; the subgrid tables, minimum/maximum bed levels and uv-point tables on subgrid grids),
* lifts the free surface of the wet cells with a depth-filtered copy of the increment, on regular, subgrid (by redistributing the cell volume) and quadtree grids. Dry cells follow the bed, so their depth does not change. The caller must not lift the water level as well,
* updates the bed level at observation points, so that the water depth in sfincs_his.nc follows the moving bed.

The filter is a Gaussian with a standard deviation of 'movingbed_filter' (default 0.8) times the mean depth of the moving cells, which mimics the Kajiura attenuation of short bed features by the water column (a single non-hydrostatic layer lacks this attenuation).
The rising and sinking parts keep their own volume.
'movingbed_filter = 0' applies the raw increment.
The earlier keywords 'nonh_movingbed_filter', 'movingbed_stamp' and 'movingbed_stamp_filter' have been removed (the surface lift is now always applied).
Increments smaller than a few float32 round-off units (ulps) of the local bed level are ignored, because they result from the difference between the caller's double-precision target bed and the single-precision kernel bed.
The kernel does not reset 'dzbext' after applying it; the caller owns the array.
'update_bed_level' updates the host arrays only, so it is meant for CPU builds.

With 'nonh_movingbed = 1' the filtered rate d(zb)/dt is also added as a kinematic source in the non-hydrostatic pressure solve.
For landslides 'nonh_movingbed = no' (default) is recommended: the surface lift alone gives the generated waves.
Use 'store_dynamic_bed_level = 1' to write the time-varying bed level 'zb' to the map file.

Example
^^^^^^^

Recommended sfincs.inp settings for a landslide (rigid slide pushed through 'dzbext'):

.. code-block:: text

	momentum_scheme          = velocity
	advection                = 1
	nonhydrostatic           = 1
	nonh_pade_beta           = 0.04
	nonh_movingbed           = no
	movingbed_filter         = 0.8
	store_dynamic_bed_level  = 1

A minimal Python loop with 'cht_sfincs', where bed_fn(t, x, y, zb0) returns the absolute bed level at time t (seconds since tstart).
The slide at its initial position should already be part of the model bed (depfile or subgrid tables), so that only its displacement is pushed:

.. code-block:: python

	import numpy as np
	from cht_sfincs.xmi import SfincsXmi

	sfx = SfincsXmi("sfincs_dll.dll", "path/to/model")
	sfx.initialize()
	sfx.get_domain()   # x, y, zs and zb (z_zmin in subgrid mode) as live arrays
	sfx.get_dzbext()   # switches use_dzbext on and exposes sfx.dzbext

	x = np.asarray(sfx.xz[:], float)
	y = np.asarray(sfx.yz[:], float)
	zb0 = np.asarray(sfx.zb[:], float).copy()

	t_start = sfx.get_start_time()
	t_end = sfx.get_end_time()
	t = sfx.get_current_time()
	while t < t_end:
	    zb_new = bed_fn(t - t_start, x, y, zb0)
	    sfx.dzbext[:] = (zb_new - np.asarray(sfx.zb[:], float)).astype(np.float32)
	    sfx.update_bed_level()   # shift bed and lift water surface
	    t = sfx.run_timestep()
	sfx.finalize()

Once the slide has stopped, set 'dzbext' to zero and finish the run with a single 'update_until' call.

References
^^^^^^^^^^

* Escalante, C., Morales de Luna, T. and Castro, M.J. (2018). Non-hydrostatic pressure shallow flows: GPU implementation using finite volume and finite difference scheme. Applied Mathematics and Computation, 338, 631-659.
* Madsen, P.A. and Sorensen, O.R. (1992). A new form of the Boussinesq equations with improved linear dispersion characteristics. Part 2. A slowly-varying bathymetry. Coastal Engineering, 18, 183-204.
