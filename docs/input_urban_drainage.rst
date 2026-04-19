Urban Drainage
==============

Overview
--------

Urban drainage is a simple bulk sink/source model for two kinds of lumped drainage infrastructure: buried pipe networks that discharge to a receiving water body (**piped drainage**) and pumps that remove water from the model and store it underground (**injection wells**). Each **drainage zone** is a polygon in the horizontal plane and has exactly one ``type``:

- ``piped_drainage`` — cells inside the polygon drain to a single outfall cell through a conceptual pipe network. Flow is bidirectional: during high water at the outfall (tide or surge), water can push back into the zone cells unless a check valve is configured. The per-zone net flux is deposited as a point source/sink at the outfall cell.
- ``injection_well`` — cells inside the polygon are pumped down at a fixed total rate, split evenly across the cells, and the extracted water is *removed from the model* (it does not reappear at an outfall). Pumping stops when the cumulative injected volume reaches the well's maximum capacity.

The approach is deliberately coarse: there is no pipe network, no hydraulic routing, no pressure head other than the difference between cell water level and outfall water level (piped_drainage) and no subsurface storage model (injection_well beyond the single capacity cap). It is intended for compound-flood applications where detailed geometry is unknown but municipal-scale design parameters (rainfall intensity, outfall pipe capacity, or pump rate + well capacity) are available.

**IMPORTANT** — urban drainage does not represent physical pipes or wells. It is a mass-balance abstraction: water disappears from urban cells, and for ``piped_drainage`` reappears summed at the outfall cell. It does not block or route flow between cells.

Inputs
------

The feature is activated by the ``urbfile`` keyword in ``sfincs.inp``:

.. code-block:: text

	urbfile                         = sfincs.urb
	store_urban_drainage_discharge  = 1
	store_cumulative_urban_drainage = 1

``store_urban_drainage_discharge`` writes per-zone time series to ``sfincs_his.nc``: the zone total discharge and (for injection wells) the cumulative injection volume. ``store_cumulative_urban_drainage`` writes the cumulative drained depth (m) per cell to ``sfincs_map.nc``.

The ``.urb`` file is a TOML document with one or more ``[[urban_drainage_zone]]`` entries.

Zone definition
---------------

Every zone has three required keys regardless of type: ``name``, ``type``, and ``polygon_file``. The rest depends on the type.

Piped drainage example
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: toml

	[[urban_drainage_zone]]
	name             = "downtown"
	type             = "piped_drainage"
	polygon_file     = "zones.tek"
	outfall_x        = 950.0
	outfall_y        = 150.0
	design_precip    = 20.0
	check_valve      = true

	[[urban_drainage_zone]]
	name             = "harbor_district"
	type             = "piped_drainage"
	polygon_file     = "zones.tek"
	outfall_x        = 1020.0
	outfall_y        = 180.0
	max_outfall_rate = 6.0

Injection well example
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: toml

	[[urban_drainage_zone]]
	name             = "north_well_field"
	type             = "injection_well"
	polygon_file     = "zones.tek"
	injection_rate   = 0.5
	maximum_capacity = 5000.0

Common keys (both types)
^^^^^^^^^^^^^^^^^^^^^^^^

``name`` (required, string)
	Zone name. Must match a polygon name in ``polygon_file``. Used as the station identifier in ``sfincs_his.nc`` when discharge output is enabled.

``type`` (required, string)
	One of ``"piped_drainage"`` or ``"injection_well"``. Selects the per-zone physics and the set of remaining required keys.

``polygon_file`` (required, string)
	Path to a Delft3D-style ``.tek`` polygon file. Multiple zones can share the same file — each zone's ``name`` is matched against polygon names inside the file. See "Polygon file format" below.

``h_threshold`` (optional, m, default ``0.0``)
	Depth over which the drainage rate ramps linearly from zero to ``q_max``. At cell ponding depth ``h_cell = 0`` the drainage is zero; at ``h_cell >= h_threshold`` it is at full ``q_max``; in between it is ``(h_cell / h_threshold) * q_max``. Smooths the discharge time series compared to a hard on/off gate. Typical values: 0.02–0.05 m. Set to ``0.0`` to reproduce the hard-cap behaviour (full ``q_max`` for any ``h_cell > 0``).

Piped drainage keys
^^^^^^^^^^^^^^^^^^^

``design_precip`` (conditional, mm/hr)
	Design rainfall intensity the zone's drainage is sized for. Per-cell capacity is ``qmax = design_precip * cell_area / 3.6e6`` [m³/s]. Typical municipal values: 10–20 mm/hr for suburban residential, 20–40 mm/hr for dense city centre. **Exactly one of** ``design_precip`` **or** ``max_outfall_rate`` **must be provided for piped_drainage zones.**

``max_outfall_rate`` (conditional, m³/s)
	Total zone outfall capacity. Useful when you know what the outfall pipe can deliver but not the design storm it was sized for. SFINCS derives ``design_precip = max_outfall_rate / zone_area * 3.6e6`` from the zone's total polygon-covered area, so per-cell capacity is distributed proportionally to cell area. **Exactly one of** ``design_precip`` **or** ``max_outfall_rate`` **must be provided for piped_drainage zones.**

``outfall_x``, ``outfall_y`` (required when ``include_outfall = true``)
	Coordinates of the single point where all zone discharge is summed and deposited. Snapped to the nearest active cell. If no active cell can be found, zone contributions are silently discarded and a warning is logged.

``include_outfall`` (optional, bool, default ``true``)
	Set to ``false`` to disable the outfall deposit step. Flow still leaves (or enters) cells, but does not reappear anywhere — treats the zone as an unconnected sink. Mostly useful for sensitivity tests.

``check_valve`` (optional, bool, default ``false``)
	When ``true``, the zone only drains outward. Backflow from the outfall into the cells (bay flooding through the pipe) is suppressed. Represents a flap valve / tide gate at the outfall.

``dh_design_min`` (optional, m, default ``0.1``)
	Floor on the per-cell design head used to compute the backflow coefficient. Per-cell backflow discharge is

	.. math::
		Q_{back}(nm) = \frac{q_{max}(nm)}{\sqrt{\max(z_b(nm) - z_b(outfall),\,\Delta h_{design,min})}} \cdot \sqrt{z_s(outfall) - z_s(nm)}

	so that a cell at the outfall bed elevation, or below it, doesn't produce an infinite backflow coefficient.

Injection well keys
^^^^^^^^^^^^^^^^^^^

``injection_rate`` (required, m³/s)
	Total pumping rate across the zone. Distributed over the zone cells by cell area so the sum across zone cells is exactly ``injection_rate`` and quadtree refinement inside the polygon does not shift the per-cell flux relative to cell area:

	.. math::
		q_{max}(nm) = \text{injection\_rate} \cdot \frac{A(nm)}{A_{zone}}

``maximum_capacity`` (required, m³)
	Total volume the injection well can accept over the simulation. SFINCS tracks ``cumulative_injection_volume`` over time; once it reaches ``maximum_capacity`` pumping is skipped for that zone (flow drops to zero). There is a potential one-time-step overshoot at the transition step (per-cell flux is not scaled to hit the cap exactly).

Polygon file format
-------------------

Zones are defined in a Delft3D-style ``.tek`` file — one or more named polygon blocks:

.. code-block:: text

	downtown
	6 2
	900.0  100.0
	900.0  200.0
	1000.0 200.0
	1000.0 100.0
	950.0   80.0
	900.0  100.0
	harbor_district
	5 2
	1000.0 150.0
	1000.0 250.0
	1100.0 250.0
	1100.0 150.0
	1000.0 150.0

Each block has a name line, a ``nrows ncols`` line, and ``nrows`` vertex lines. If the last vertex does not equal the first, SFINCS closes the ring automatically at read time. A cell center falling inside multiple polygons is assigned to the **last** zone encountered — overlap warnings are not emitted, so order your zones with intent.

Flow formulas
-------------

For each active cell ``nm`` inside zone ``iz``, with ``h_cell`` the cell ponding depth (``zs - subgrid_z_zmin`` in subgrid mode, ``zs - zb`` otherwise) and the rate ramp

.. math::
	r = \min(h_{cell} / h_{threshold},\; 1) \text{ if } h_{threshold} > 0, \text{ else } 1

**Piped drainage** (``outfall cell io``):

.. math::
	\Delta z_s = z_s(nm) - z_s(io)

Outflow (``Δz_s > 0`` and ``h_cell > 0``):

.. math::
	q = \min\left(r \cdot q_{max}(nm),\; \frac{h_{cell} \cdot A(nm)}{\Delta t}\right)

Backflow (``Δz_s < 0`` and check valve off):

.. math::
	q = -\frac{q_{max}(nm)}{\sqrt{\max(z_b(nm) - z_b(io),\,\Delta h_{design,min})}} \cdot \sqrt{-\Delta z_s}

capped at ``-q_{max}(nm)``. With ``check_valve = true`` backflow is skipped entirely.

The zone's per-step net flux is deposited at the outfall cell, so mass is conserved (up to the outfall-snap warning above).

**Injection well** (no outfall):

.. math::
	q = \min\left(r \cdot q_{max}(nm),\; \frac{h_{cell} \cdot A(nm)}{\Delta t}\right)

where :math:`q_{max}(nm) = \text{injection\_rate} \cdot A(nm) / A_{zone}` (area-weighted split; sums to ``injection_rate`` across the zone). Flow is positive (water leaves cells) only; there is no backflow. Pumping is skipped entirely once :math:`\text{cumulative\_injection\_volume}(iz) \geq \text{maximum\_capacity}(iz)`.

Outputs
-------

**``sfincs_his.nc``** — when ``store_urban_drainage_discharge = 1``:

``urban_drainage_discharge(urban_drainage_zones, time)``
	Per-zone total discharge in m³/s. Positive means net outflow from the cells (to outfall or to injection well); negative means net inflow (backflow from the outfall, piped_drainage with check valve off). Named ``urban drainage zone total discharge`` in the long_name attribute.

``cumulative_injection_volume(urban_drainage_zones, time)``
	Per-zone cumulative injection volume in m³. Tracked for all zones, but only physically meaningful for ``injection_well`` zones; ``piped_drainage`` zones keep this at 0.0 (there is no underground storage).

``urban_drainage_zone_name(urban_drainage_zones)``
	Zone names, in the order they appear in the ``.urb`` file.

**``sfincs_map.nc``** — when ``store_cumulative_urban_drainage = 1``:

``urban_drainage_cumulative_depth(m, n, timemax)`` (regular) or ``(nmesh2d_face, timemax)`` (quadtree)
	Cumulative drained volume divided by cell area (m), written at the ``dtmaxout`` interval. Positive means net outflow from the cell over the simulation; negative means net inflow.

At init time a per-zone summary block is written to the SFINCS log listing zone name, type, polygon file, number of cells assigned, total area, design precipitation / max outfall rate / qmax (piped_drainage) or injection rate / maximum capacity (injection_well), thresholds, outfall coords + snapped cell index (piped_drainage), and check-valve state.
