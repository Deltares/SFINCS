# SFINCS code review checklist

Use this when reviewing a pull request with Fortran changes under `source/`. [CODING_STYLE.md](CODING_STYLE.md) covers how the code should look; this checklist covers whether a change is correct, complete and wired in everywhere it needs to be.

AI agents using this checklist: in case of doubt (whether a rule applies, which input format or grid type a new option should support, whether a difference is a regression), ask the user explicitly instead of guessing or picking a default.

Not every item applies to every PR. Go through the sections that the diff touches, and report findings as **blocking** (wrong results, crashes, broken builds, missing GPU data) or **non-blocking** (clean-up, documentation, style).

## 1. Declarations and data

- Every new variable is declared, with the right kind (`real*4` vs `real*8`, `integer` vs `integer*8` for large counts). No reliance on implicit typing: `implicit none` is present in every new subroutine and function.
- Mixed-precision arithmetic is deliberate (e.g. `real*4` state updated with a `real*8` time step).
- Global (model-state) variables are declared in `sfincs_data.f90` only. Locals stay local.
- New `use` statements go inside the subroutine that needs them, not at module level, except `use sfincs_data`.
- New arrays are allocated with the right size (`np`, `npuv`, `nrsrc`, ...) and initialised before first use. Check the path where the feature is switched off: arrays referenced in a kernel still need to be allocated (often with size 1) when the process is inactive.
- Before re-allocating, arrays are guarded with `if (allocated(x)) deallocate(x)` when the routine can run more than once (e.g. a model re-initialised through BMI).
- New allocatable arrays are added to `finalize_parameters` in `sfincs_data.f90`, so that the deallocation list stays complete.
- Counters, flags and accumulated values (maxima, cumulative volumes, timers) are reset at initialisation, not only at declaration.
- New state variables that evolve in time (e.g. storage volumes, infiltration state, structure status) are written to the restart file in `sfincs_ncoutput.F90` and read back in `sfincs_initial_conditions.F90`. Otherwise a restarted run silently starts that process from zero.

## 2. GPU (OpenACC) and OpenMP

- Every new array used inside an `!$acc` region is on the device: added to the `!$acc enter data` in `initialize_openacc` and to the matching `!$acc exit data delete` in `finalize_openacc` (both in `sfincs_openacc.f90`).
- Arrays changed on the host during the run (forcing read from file, boundary updates, BMI `set_value`) are copied to the device with `!$acc update device(...)` before the next kernel uses them.
- Arrays computed on the device and written to output are copied back with `!$acc update host(...)` before writing (see `sfincs_ncoutput.F90`). Check both map and his output, and BMI `get_value`.
- New loops in the time loop have an `!$acc` directive (and `!$omp` where the CPU loop is parallel), consistent with neighbouring loops. A loop without a directive inside a device region silently runs on the host with stale data.
- Loop-local scalars are `private`; sums, maxima and minima use a `reduction` clause. No two iterations write the same element (race conditions, e.g. when looping over velocity points and updating cell values).
- No subroutine calls, `write` statements or allocation inside device kernels, unless the routine is declared `!$acc routine`.
- The same behaviour on CPU and GPU: if the testbed is run on only one of them, ask which.
- Reductions (sums, maxima over cells) can differ slightly between CPU and GPU, and between GPU runs, because the summation order is not fixed. Differences at that level are not regressions; larger differences are.

## 3. Build files

- New source files are added to `source/src/Makefile.am`, in dependency order: a module must be listed after every module it `use`s (the Linux build compiles in list order).
- The file name in `Makefile.am` matches the file on disk exactly, including the extension case (`.F90` is preprocessed, `.f90` is not, and Linux file names are case-sensitive).
- New source files are added to `source/sfincs_lib/sfincs_lib.vfproj` for the Windows (Intel ifx) build.
- Code that needs the preprocessor (`#ifdef`, `#include`) is in a `.F90` file.
- Compiler-specific constructs compile with all three compilers: Intel ifx (Windows), gfortran (Linux CPU Docker) and nvfortran (Linux GPU).

## 4. Input

- New keywords are read in `sfincs_input.f90` with `get_keyword` and a sensible default, so that existing models run unchanged when the keyword is absent.
- New SnapWave keywords are read in `sfincs_snapwave.f90`, with defaults.
- Keyword names are lowercase and consistent with existing names (e.g. `*file` for file names).
- New input files (forcing, infiltration and other spatial or time-varying input) are NetCDF only, using the readers in `sfincs_ncinput.F90` where they fit (e.g. the quadtree readers for spatial input). Existing ASCII and binary formats stay supported for backward compatibility, but no new ASCII or binary formats or readers are added.
- The same applies to the grid: new functionality only has to be available for models with quadtree netCDF grid input (`qtrfile`, typically `sfincs.nc`). Regular-grid input (`depfile`, `mskfile`, `indexfile`) stays supported for backward compatibility, but new options do not need to support it and no new regular-grid input files are added.
- Exception: parameter and lookup-table input that is not a grid or time series uses the TOML format. New options in these areas are added to the existing TOML file, not to a new file or a NetCDF file:
  - Drainage structures (`drnfile`: pumps, culverts, gates, dike breaches and their control rules) and dike-breach parameters (`dkbfile`), read in `sfincs_src_structures.f90`. Legacy fixed-column `.drn` files are converted to TOML on reading, so new options do not need to be added to the legacy format. Document new keys in `docs/input_structures.rst`.
  - Vegetation types (`vegtype`): the lookup table with the characteristics per vegetation type, read in `sfincs_vegetation.f90`. The vegetation type per cell itself is spatial input (`vegfile`).
- New input files are checked for existence, and missing or malformed files produce a clear error through `sfincs_error`/`sfincs_log`, not a crash.
- Renamed or removed keywords keep the old name working as an alias, with a deprecation warning in the log.
- Changed defaults are deliberate and mentioned in the PR description, since existing models will give different results.
- New keywords are documented in `docs/parameters.rst` (and the relevant `docs/input*.rst` page).
- When a new type of input is needed, an issue is opened in HydroMT-SFINCS.

## 5. Output

- New NetCDF variables have `units`, `long_name` and, where one exists, a `standard_name`, and the correct dimensions for regular and quadtree grids.
- New output is written to both map and his files where relevant, with the right fill value for dry or inactive points.
- Output-only arrays are copied from the device before writing (see section 2).
- New output variables are documented in `docs/output.rst` and, when switched by a keyword, in `docs/parameters.rst`.

## 6. Numerics

- Divisions are protected against zero (water depth, area, distance, wave period). Check dry cells and cells at the wet/dry threshold `huthresh`.
- Array indices are correct (1-based) and stay in bounds at domain edges, boundary cells and quadtree refinement-level transitions.
- The change works on quadtree grids, with and without subgrid tables, including at refinement-level transitions. Support for regular-grid input is not required (see section 4), but a new option that does not support it stops with a clear error (via `stop_sfincs`) when used with a regular-grid model, rather than giving wrong results.
- Changes to existing functionality keep working for regular-grid models.
- Explicit time stepping and the CFL condition are respected; new terms do not make the scheme unstable for large time steps.
- Mass is conserved: new sources and sinks are included in the continuity update and in volume bookkeeping.

## 7. Integration

- New processes are initialised and wired through `sfincs_domain.f90` and the time loop in `sfincs_lib.f90`, not directly in `sfincs.f90`.
- Changes to the public interface of `sfincs_bmi.f90` are flagged as breaking changes.
- Effects on SnapWave coupling (`sfincs_snapwave.f90`) are considered.
- New timers are added to `sfincs_timers.f90` and reported in the log.
- Errors stop the model through `stop_sfincs` (`sfincs_error.f90`), never with a bare `stop`: when SFINCS runs as a library through BMI, `stop` kills the host process.
- Messages to the user go through `write_log` (`sfincs_log.f90`), not `print` or `write(*,*)`, so they also end up in the log file.
- Headers list the correct callers (`Called from:`), including for routines whose callers changed.

## 8. Performance and housekeeping

- No allocation, file I/O or work that could be done once at initialisation inside the time loop (`sfincs_update`).
- Inactive processes cost nothing: the call is skipped as a whole with an `if` on the process switch, not checked per cell.
- No debug output, hard-coded paths or large commented-out blocks left behind.
- Third-party code copied into the repository has a licence compatible with GPL-3 (LGPL for `snapwave/`), and its origin is stated in the file header.

## 9. Testing and regressions

- List the modules affected by the change and the testbed cases that cover them.
- The testbed has been run; differences in reference results are expected and explained.
- New functionality comes with a conceptual test and a real-world test case in the testbed.
- Bug fixes come with a case that reproduces the bug, where practical.
