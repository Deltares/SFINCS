# SFINCS developer guide

How the SFINCS Fortran code is organised, and where new code goes. Read this before adding a feature. The rules a finished change must meet are in [CODING_REVIEW.md](CODING_REVIEW.md); formatting is in [CODING_STYLE.md](CODING_STYLE.md).

The module headers in `source/src` describe each routine and its callers. This guide gives the overview that the headers do not. Keep it up to date: a change that adds a module or adds, moves or removes a call in the initialisation or the time loop also updates this file.

AI agents: in case of doubt (where a call belongs, which input format to use, whether to support regular grids), ask the user explicitly instead of guessing.

## Principles

- **Keep it simple.** Choose the simplest solution that does the job. No extra options, abstraction layers, generic frameworks or "for later" code that the feature does not need now. A simpler solution is easier to review, to run on the GPU and to maintain, and is almost always better.
- **Performance always matters.** SFINCS is used for large domains, ensembles and real-time forecasting, on CPU and GPU. Code in the time loop runs millions of times: avoid work that can be done once at initialisation, avoid allocation and file I/O in the loop, skip inactive processes as a whole, and keep loops GPU-friendly. Measure the run time before and after a change to the time loop (CPU, and GPU where possible). Extra complexity for speed is only worth it when the gain is measured.
- **Reuse what exists.** Before writing a new routine, look for one that already does it: the keyword readers in `sfincs_input.f90` and `sfincs_read.f90`, the netCDF and quadtree readers in `sfincs_ncinput.F90`, the polygon helpers in `utils/`, the date utilities in `sfincs_date.f90`, `write_log` and `stop_sfincs`, the timers. Extend an existing routine where that is natural rather than writing a near-copy.

## 1. Code map

Driver and data:

| File | Role |
|---|---|
| `sfincs.f90` | Main program: calls `sfincs_initialize`, `sfincs_update` and `sfincs_finalize`. |
| `sfincs_lib.f90` | Driver: initialisation sequence, time loop and finalisation. Used by the executable and the BMI. |
| `sfincs_bmi.f90` | Basic Model Interface for coupling with other programs. Public interface: changes are breaking. |
| `sfincs_data.f90` | All global variables: parameters, flags, grid and index arrays, model state. |
| `sfincs_input.f90` | Reads `sfincs.inp` (`get_keyword`), with defaults. |
| `sfincs_domain.f90` | Sets up the domain: process flags, mesh, bathymetry, boundaries, roughness, storage, vegetation, initial state. |
| `sfincs_openacc.f90` | Copies arrays to and from the GPU at the start and end of the run. |

Solver:

| File | Role |
|---|---|
| `sfincs_momentum.f90` | `compute_fluxes`: fluxes `q` and velocities `uv` at the uv points. Terms that act on the flow between cells go here. |
| `sfincs_continuity.f90` | `update_continuity`: collects sources and sinks in `qsrc` and updates the water levels (regular or subgrid). |
| `sfincs_nonhydrostatic.f90` | Non-hydrostatic pressure correction (regular grids only). |
| `sfincs_bathtub.f90` | Bathtub mode: water levels without solving the flow. |

Grid and input:

| File | Role |
|---|---|
| `sfincs_quadtree.F90` | Reads the quadtree netCDF file (`qtrfile`) and builds the mesh. |
| `sfincs_subgrid.F90` | Reads the subgrid tables. |
| `sfincs_initial_conditions.F90` | Initial water levels, including restart files. |
| `sfincs_ncinput.F90` | netCDF readers, including the generic quadtree readers for spatial input. |
| `sfincs_read.f90`, `sfincs_spiderweb.f90` | Generic keyword readers, and legacy ASCII spiderweb and gridded meteo readers. |

Forcing and processes:

| File | Role |
|---|---|
| `sfincs_boundaries.f90` | Open boundary conditions (water levels, waves on the boundary). |
| `sfincs_meteo.f90` | Wind, pressure and precipitation. |
| `sfincs_discharges.f90` | River point discharges. |
| `sfincs_src_structures.f90` | Drainage structures (pumps, culverts, gates, dike breaches), TOML input (`drnfile`, `dkbfile`). |
| `sfincs_rule_expression.f90` | Rule language for opening and closing drainage structures. |
| `sfincs_structures.f90` | Weirs and thin dams on uv points. |
| `sfincs_infiltration.f90` | Infiltration methods (constant, Curve Number, Green-Ampt, Horton, bucket). |
| `sfincs_urban_drainage.f90` | Urban drainage per zone. |
| `sfincs_vegetation.f90` | Vegetation per cell, TOML lookup table (`vegtype`). |
| `sfincs_wavemaker.f90` | Wave makers along polylines. |
| `sfincs_wave_enhanced_roughness.f90` | Wave-enhanced bed roughness. |
| `sfincs_snapwave.f90` | Coupling with the SnapWave wave solver (`snapwave/`). |
| `sfincs_advection_diffusion.f90` | Tracer advection and diffusion (experimental). |

Output and utilities:

| File | Role |
|---|---|
| `sfincs_output.f90` | Decides when map, max, his and restart output is due. |
| `sfincs_ncoutput.F90`, `sfincs_ncoutput_helpers.F90` | Writes the netCDF output and restart files. |
| `sfincs_obspoints.f90`, `sfincs_crosssections.f90`, `sfincs_runup_gauges.f90` | Observation points, cross-sections and run-up gauges. |
| `sfincs_log.f90`, `sfincs_error.f90`, `sfincs_timers.f90` | Log output (`write_log`), stopping with an error (`stop_sfincs`), timers. |
| `sfincs_date.f90`, `sfincs_timestep_analysis.f90`, `utils/` | Date utilities, time-step analysis, polygon helpers. |

`snapwave/` holds the SnapWave wave solver (LGPL). `source/third_party_open/` holds external code that is not modified.

## 2. Initialisation order

`sfincs_initialize` in `sfincs_lib.f90` runs these steps in order:

1. `read_sfincs_input`: all keywords from `sfincs.inp`.
2. `read_meteo_data` (in bathtub mode only for SnapWave wind growth).
3. `initialize_domain` (`sfincs_domain.f90`):
   1. `initialize_processes`: sets process flags (wind, precipitation, quadtree, ...) from the input.
   2. `initialize_mesh`: reads the quadtree file or builds the regular grid; indices and neighbours.
   3. `initialize_bathymetry`: bed levels, subgrid tables.
   4. `initialize_boundaries`, `initialize_roughness`, `initialize_storage_volume`, `initialize_vegetation`.
   5. `initialize_hydro`: initial water levels and fluxes (including restart).
4. Structures, boundary data, observation points, cross-sections, run-up gauges.
5. Processes: `initialize_infiltration`, `initialize_discharges`, `initialize_src_structures`, `initialize_urban_drainage`, then non-hydrostatic, wave makers and bathtub when switched on.
6. SnapWave coupling, when switched on. Then, in bathtub mode with SnapWave, `initialize_bathtub_wave_setup` adds `bathtub_fachs` times the boundary wave height to the boundary water levels.
7. `initialize_output`.
8. `deallocate_quadtree`: the quadtree neighbour arrays (`quadtree_md`, `quadtree_mu`, ...) are freed here. Code that needs them must run before this point.
9. `initialize_openacc`: all device arrays are copied to the GPU. Every array used on the GPU must be allocated by now.

A new process that has its own input is initialised in step 5. A new spatial field that belongs to the domain (like roughness or storage volume) goes in step 3.

## 3. Time loop

`sfincs_update` in `sfincs_lib.f90` runs, per time step:

1. Output due at this time (map, max, restart, his) is flagged.
2. Meteo fields and SnapWave are updated when their update interval is reached.
3. `update_meteo_forcing`: wind stress, pressure, precipitation.
4. `update_boundaries`: boundary water levels and velocities.
5. In bathtub mode, `bathtub_compute_water_levels`, then `update_wave_field` (SnapWave, with `bathtub_snapwave`). Otherwise:
   1. `update_wave_field` (SnapWave).
   2. `compute_fluxes` (`sfincs_momentum.f90`): fluxes at the uv points.
   3. Wave-maker fluxes, fluxes over structures (weirs, thin dams), non-hydrostatic correction.
   4. `update_continuity` (`sfincs_continuity.f90`): river discharges, drainage structures, infiltration and urban drainage add to `qsrc`; then the water levels are updated.
6. `write_output` when output is due.

Where a new term goes:

- **Source or sink of water** (m3/s per cell): add it to `qsrc` from an `update_<feature>` routine called in `update_continuity`, next to `update_discharges` and `update_infiltration_map`. `qsrc` is reset after each water-level update.
- **Force on the flow** (friction, stress, drag at uv points): in `compute_fluxes`, or as a quantity at the uv points that `compute_fluxes` uses.
- **Flux over a structure or line**: next to `compute_fluxes_over_structures`, between momentum and continuity.
- **Forcing that changes in time** (read from file): updated before momentum, like `update_meteo_forcing` and `update_boundaries`, and copied to the GPU with `!$acc update device`.

## 4. Adding a new process

1. Create `source/src/sfincs_<feature>.f90` with its own module and `initialize_<feature>` / `update_<feature>` routines. `sfincs_infiltration.f90` and `sfincs_urban_drainage.f90` show the pattern. ([CODING_REVIEW.md](CODING_REVIEW.md) section 7)
2. Add global variables and the process switch to `sfincs_data.f90`. (section 1)
3. Add keywords to `sfincs_input.f90` with defaults and a trailing comment. Choose the input format by the policy: NetCDF for spatial or time-varying input, quadtree models only, TOML for structures and vegetation types. (section 4)
4. Call `initialize_<feature>` in the initialisation order (section 2 above) and `update_<feature>` in the time loop (section 3 above), skipped as a whole when the process is off. (sections 7 and 8)
5. Add arrays used on the GPU to `sfincs_openacc.f90`, and put `!$acc` and `!$omp` directives on the new loops. (section 2)
6. Add output variables, and state variables to the restart file. (sections 1 and 5)
7. Add the file to `Makefile.am` (after every module it uses) and to `sfincs_lib.vfproj`. (section 3)
8. Document keywords and output in `docs/`, and add a conceptual and a real-world test to the testbed. (sections 4, 5 and 9)
9. Add the new module to the code map in this file, and the new calls to sections 2 and 3. (section 7)

Then check the whole change against [CODING_REVIEW.md](CODING_REVIEW.md) before committing or opening a pull request.

## 5. Fixing a bug

- Reproduce it first, preferably with a small testbed case.
- Read the module and its callers before changing anything; the `Called from:` lines in the headers help.
- Fix the root cause, not the symptom, and look for the same pattern elsewhere, in particular in the regular-grid and quadtree versions and the regular and subgrid versions of a routine.
- Check CPU and GPU behaviour: a fix in a loop may need a matching change in its `!$acc` directive or in `sfincs_openacc.f90`.

## 6. Opening a pull request

Fill in [.github/pull_request_template.md](../.github/pull_request_template.md) completely: the description (what changes and why, which keywords and defaults, the effect on results and run time) and every checklist item. The `pr-checklist` workflow fails as long as an item is unticked.

Tick an item only when it is true:

- **Done:** tick it, only when it was actually done. Compilation boxes only for builds that were really run on that platform (Windows, Linux CPU, Linux GPU). Testbed and test boxes only when the testbed was run or the tests were added.
- **Not relevant:** tick it and add "n/a" behind it, e.g. the documentation items when no documentation or user-facing behaviour changed, the test items when no model functionality was added or changed, the HydroMT-SFINCS item when no new type of input is needed, the figure item for small changes.
- **Relevant but not done:** leave it unticked and say in the description what is missing (e.g. "not compiled on Linux GPU"). Never tick it to make the check pass.

AI agents: never tick **"I have reviewed my own changes."** That box is for the human developer, after reading the diff. Leave it unticked, so the pull request cannot pass the checklist check before a person has reviewed it.
