# SFINCS Fortran coding style

Complete formatting conventions for all SFINCS and HurryWave Fortran edits — apply to every edit, not just reformatting tasks.

## Indentation

- Module-level declarations and `contains` keyword: 3-space indent
- Subroutine/function body: 6-space indent (3 spaces relative to the `subroutine` keyword)
- Each nested control block (`if`/`do`/etc.) adds 3 more spaces

## Assignments

Align `=` signs within any logical block of assignments. The one exception is when a single outlier name is dramatically longer than the rest and would create excessive whitespace — in that case, break it into a separate sub-group or leave it unaligned:

```fortran
      ! wrong
      x1 = 313
      xx2 = 1324
      a = 123

      ! correct
      x1  = 313
      xx2 = 1324
      a   = 123

      ! correct — logical block; align even when name lengths vary
      nmindobs(iobs) = nm
      n              = quadtree_n(nm)
      m              = quadtree_m(nm)
      iref           = quadtree_level(nm)
      zobs(iobs)     = quadtree_zb(nm)
```

## Variable naming

All-lowercase snake_case. No camelCase, PascalCase, or SHOUTY_CASE:

```fortran
      drn_type, nm, isrc, it_prev, qsrc_ts   ! correct
      drnType, NM, iSrc, itPrev              ! wrong
```

Exceptions only for external library API names or well-established mathematical constants.

## Module structure

A blank `!` line goes directly after `module <name>` and directly after `contains`:

```fortran
module sfincs_continuity
   !
   use ...
   !
contains
   !
   subroutine ...
```

## Subroutine separators

Between subroutines/functions inside a module, use 2 blank `!` lines. Same banner before the first subroutine after `contains`:

```fortran
   end subroutine foo
   !
   !
   subroutine bar()
```

No other banner/divider lines anywhere (`! ---`, `! ===`, `! ***`, etc.). Use plain blank `!` lines for visual separation instead.

## Blank `!` spacers inside routines

**if-blocks**: `!` after every `then`, `!` before every `endif`/`else`/`elseif`, and `!` between logical groups inside the block:

```fortran
      if (condition) then
         !
         x = 1
         y = 2
         !
         if (other) then
            !
            z = 3
            !
         endif
         !
         call foo()
         !
      endif
```

**do-loops**: `!` after every `do`, `!` before every `enddo`. For regular (non-directive) loops also add `!` before `do` and after `enddo`. Skip the surrounding `!` when an `!$omp` / `!$acc` directive sits immediately adjacent — the directive already provides the visual break:

```fortran
      do nm = 1, np
         !
         ... body ...
         !
      enddo
```

**Exception — tight nested loop-open lines**: do NOT add `!` between consecutive bare `do` lines:

```fortran
      do nm = 1, np        ! correct: no ! between these
         do i = 1, ntheta
            do j = 1, nsigma
```

**Descriptive comments**: bracket any descriptive `! text` comment with blank `!` lines on both sides — never place code directly above or below:

```fortran
      !
      ! Compute friction term using Manning's n.
      !
      qq = ...
```

## Documentation headers

Note: most existing headers were generated automatically (October 2026) and have not all been reviewed. The code is the reference; correct a header when you find it wrong.

**Module header**: start with the line `! Headers auto-generated, ask main developers in case of doubts` (also in every new file), then list each public subroutine/function with a one-line description and its caller(s):

```fortran
module sfincs_continuity
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! Continuity solver: updates water levels each time step.
   !
   ! Subroutines:
   !
   !   update_continuity(t, dt)
   !     Main per-timestep entry
   !     Called from: sfincs_lib
   !
   !   compute_water_levels_regular(dt)
   !     Regular-grid water-level update
   !     Called from: update_continuity
   !
```

**Licence block**: the module header of every file in `source/src` ends with the block below, unchanged. Its first and last lines are the only banner lines allowed:

```fortran
   ! ------------ Copyright notice SFINCS ------------
   !
   ! Copyright (C)  Stichting Deltares, 2016-2026.
   !
   ! This file is part of SFINCS.
   !
   ! SFINCS is free software: you can redistribute it and/or modify
   ! it under the terms of the GNU General Public License as
   ! published by the Free Software Foundation version 3.
   !
   ! SFINCS is distributed in the hope that it will be useful,
   ! but WITHOUT ANY WARRANTY; without even the implied warranty of
   ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   ! GNU General Public License for more details.
   !
   ! You should have received a copy of the GNU General Public License
   ! along with SFINCS.  If not, see <http://www.gnu.org/licenses/>.
   !
   ! contact: software@deltares.nl
   ! Stichting Deltares
   ! P.O. Box 177
   ! 2600 MH Delft, The Netherlands
   !
   ! All indications and logos of, and references to, "Delft3D"
   ! and "Deltares" are registered trademarks of Stichting Deltares,
   ! and remain the property of Stichting Deltares. All rights reserved.
   !
   ! ------------ Copyright notice SFINCS ------------
   !
```

**Per-subroutine header**: describe what it does and list (module -> subroutine) who calls it:

```fortran
   subroutine update_continuity(t, dt)
      !
      ! Advance water levels over one time step.
      !
      ! Called from: sfincs_lib
      !
```

**Arguments**: when a subroutine/function has arguments, describe them between the description and `Called from:`, one line per argument with its meaning and unit. Group them under `Input:`, `Input/output:`, `Output:` and (for a function result) `Returns:`, and align the colons. Routines without arguments get no such block:

```fortran
   subroutine read_spw_dimensions(filename, nt, nrows, ncols, spwrad, nquant)
      !
      ! Read the dimensions (times, rows, columns, radius) of an ASCII spiderweb file.
      !
      ! Input:
      !   filename : name of the meteo file
      !
      ! Output:
      !   nt       : number of time blocks in the file
      !   nrows    : number of rows
      !   ncols    : number of columns
      !   spwrad   : radius of the spiderweb grid [m]
      !   nquant   : number of quantities per time block
      !
      ! Called from: read_meteo_data (sfincs_meteo)
      !
```

## Naming of end statements

Always name `end subroutine` and `end module`:

```fortran
   end subroutine update_continuity

end module sfincs_continuity
```

## implicit none

Every subroutine and function must have `implicit none`.

## Application scope

Apply all rules to every edit — bug fixes, new features, and reformatting alike. When editing existing code, fix style locally in the blocks you touch; do not sweep-format untouched files unless the user explicitly asks.

**Why:** User-defined rules developed across SFINCS and HurryWave sessions; consolidated here so both projects share one authoritative source.
