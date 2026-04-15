# SFINCS Fortran Model — Claude Agent Configuration

## Project Overview

SFINCS (Super-Fast INundation of CoastS) is a reduced-complexity hydrodynamic model written in Fortran 90/95.
It solves the depth-averaged shallow water equations on regular grids, quadtree grids, and subgrid-based meshes
for coastal/fluvial/pluvial flood simulations.

## Repository Layout

```
source/
  src/
    sfincs.f90                  # Main program entry point
    sfincs_bmi.f90              # Basic Model Interface (BMI) — public API
    sfincs_domain.f90           # High-level orchestration (initialize/update/finalize)
    sfincs_data.f90             # Global variables and arrays
    sfincs_input.f90            # Input parameter parsing (sfincs.inp)
    sfincs_read.f90             # File readers (grid, bathymetry, forcings)
    sfincs_momentum.f90         # Momentum equations (u/v velocity)
    sfincs_continuity.f90       # Continuity equation (water level h/zs)
    sfincs_boundaries.f90       # Water level and wave boundary conditions
    sfincs_discharges.f90       # Discharge source/sink points
    sfincs_meteo.f90            # Wind, pressure, precipitation forcing
    sfincs_output.f90           # NetCDF/binary output
    sfincs_structures.f90       # Weirs, thin dams, drainage structures
    sfincs_infiltration.f90     # Green-Ampt / curve number infiltration
    sfincs_vegetation.f90       # Vegetation drag
    sfincs_snapwave.f90         # SnapWave coupling interface
    sfincs_advection_diffusion.f90
    sfincs_nonhydrostatic.f90
    sfincs_crosssections.f90
    sfincs_obspoints.f90
    sfincs_runup_gauges.f90
    sfincs_wavemaker.f90
    sfincs_wave_enhanced_roughness.f90
    sfincs_spiderweb.f90
    sfincs_openacc.f90          # GPU acceleration
    sfincs_timestep_analysis.f90
    sfincs_bathtub.f90
    sfincs_date.f90
    sfincs_log.f90
    sfincs_error.f90
    snapwave/                   # SnapWave wave propagation submodule (Fortran)
  third_party_open/             # External libraries (netcdf-fortran, Delft3D utils, bicgstab)
docs/                           # Sphinx documentation (.rst files)
```

## Key Conventions

### Code Style
- Fortran 90/95; modules with `use` statements (never `include` for data sharing)
- `implicit none` in every subroutine and module procedure
- All global state lives in `sfincs_data.f90` — no new global variables in other modules
- Variable naming: grid coordinates `xcor`/`ycor`, bathymetry `zb`/`dep`, water level `zs`/`h`, velocities `u`/`v`
- Mask arrays: `msk` (active cells), `kfu`/`kfv` (velocity point masks)
- New input parameters must be added to `sfincs_input.f90` with a default value and read from `sfincs.inp`

### Module Integration Pattern
When adding a new physics process:
1. Create `sfincs_<process>.f90` with `initialize_<process>()`, `update_<process>()`, `finalize_<process>()` subroutines
2. Add `use sfincs_<process>` inside the relevant subroutine in `sfincs_domain.f90` (NOT at module level)
3. Call `initialize_<process>()` from `initialize_domain()` in `sfincs_domain.f90`
4. Call `update_<process>()` at the correct point in the time loop
5. Wire input keywords through `sfincs_input.f90`
6. **Never break the BMI interface** in `sfincs_bmi.f90` — this is the public API used by coupling frameworks

### BMI Interface
`sfincs_bmi.f90` is the stable public API. Functions: `initialize`, `update`, `update_until`, `finalize`, `get_var`, `set_var`, `get_time_step`, `get_current_time`, etc. Do not rename or remove existing BMI functions.

### Build System
CMake-based. Build targets: `sfincs` (executable), `sfincs_dll` (shared library for BMI coupling), `sfincs_lib` (static library).

### NetCDF Output
All gridded output uses NetCDF via `sfincs_output.f90`. Follow existing patterns for defining dimensions, variables, and attributes. Coordinate reference system metadata must be included.

## Roles and Slash Commands

### `/review-fortran` — Code Reviewer
Review Fortran changes for:
- Correctness of numerical schemes (finite differences, time stepping, stability)
- Array bounds safety and allocation/deallocation hygiene
- Consistency with `sfincs_data.f90` variable naming and global state
- BMI interface integrity (no breaking changes to `sfincs_bmi.f90`)
- Module `use` placement (use inside subroutines, not at module level — except `sfincs_data`)
- Memory leaks (allocated arrays must be deallocated in `finalize_*` routines)
- OpenACC GPU directives correctness (in `sfincs_openacc.f90`)
- Regression risk to existing functionality

### `/develop-fortran` — Code Developer
When implementing new features:
- Follow the module integration pattern above
- Add input keywords with defaults to `sfincs_input.f90`
- Write the new module as `sfincs_<feature>.f90`
- Integrate via `sfincs_domain.f90`, not directly in `sfincs.f90`
- Do not modify third_party_open/ libraries
- Keep numerical schemes consistent with the existing explicit time-stepping approach
