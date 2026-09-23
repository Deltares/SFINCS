# Groundwater solver (`gwflow`)

## What it is

`gwflow = 1` adds a 2D depth-averaged (Dupuit) unconfined aquifer to SFINCS:

    Sy dh/dt = div( K (h - zbase) grad h ) + recharge + exchange

With `semi_implicit = 1` the aquifer is solved as extra rows in the **same** semi-implicit matrix
as the surface (surface head in rows 1..N, aquifer head in rows N+1..2N), with a symmetric
exchange term coupling each surface/aquifer pair. Transmissivity `K*(h - zbase)` depends on the
solution, so the system is mildly nonlinear; that nonlinearity is lagged one Picard (outer)
iterate to keep the matrix symmetric and let CG solve it. With `semi_implicit = 0` the aquifer is
stepped explicitly, sub-cycled within one surface timestep to respect a diffusion-number limit
(`gw_numax`).

Source: `source/src/sfincs_groundwater.f90` (the solver), `source/src/sfincs_semi_implicit.f90`
(aquifer rows in the coupled matrix, `gw_` prefix), `source/src/sfincs_infiltration.f90`
(infiltration-to-recharge hookup).

## Keywords

All defaults and units below are read directly from `source/src/sfincs_input.f90`.

| keyword | type | default | unit | meaning |
|---|---|---|---|---|
| `gwflow` | int (0/1) | 0 | - | turns the aquifer on |
| `gw_kh` | real | 1.0e-2 | m/s | uniform horizontal hydraulic conductivity |
| `gw_sy` | real | 0.05 | - | uniform specific yield |
| `gw_zbase` | real | 0.0 | m | uniform aquifer base elevation |
| `gw_recharge` | real | 0.0 | m/s | uniform background recharge rate (negative = a sink, e.g. drainage) |
| `gw_leakance` | real | 1.0e-5 | 1/s | surface/aquifer exchange conductance per unit area |
| `gw_theta` | real | 0.75 | - | time weighting of the aquifer diffusion term (matches `theta_si`'s default) |
| `gw_numax` | real | 0.25 | - | diffusion-number cap, `Nu = K*b*dt/(Sy*dx^2)`, that limits the explicit path's sub-step |
| `gw_zsini` | real | -999.0 | m | initial aquifer head; if left at default the initial head is the initial surface level `zs`, clamped to `gw_zbase` |
| `gw_headfile` | char | 'none' | - | spatial initial head field (see Files) |
| `gw_rechargefile` | char | 'none' | - | spatial recharge field (see Files) |
| `gw_khfile` | char | 'none' | - | spatial conductivity field (see Files) |
| `gw_syfile` | char | 'none' | - | spatial specific-yield field (see Files) |
| `gw_zbasefile` | char | 'none' | - | spatial aquifer-base field (see Files) |
| `gw_bnd_from_zs` | int (0/1) | 0 | - | at open-boundary cells (`kcs == 2`) the aquifer head is reset to the current surface level `zs` every outer iterate, instead of staying at whatever it was given |
| `gw_from_infiltration` | int (0/1) | 0 | - | routes infiltrated water into `gw_recharge` instead of discarding it; requires `gwflow = 1` and cannot be combined with `gw_rechargefile` (both would write `gw_recharge`) |
| `gw_seepage_fac` | real | 1.0 | - | seepage-face strength as a multiple of the volume stored above the ceiling; 1.0 removes exactly that volume in one timestep; 0.0 disables the seepage face (kept only to reproduce the old defect where water above the ceiling was simply lost) |
| `gw_storage_mode` | int (0/1/2) | 1 | - | storage convention above the ground under standing water: `0` = pre-2026-09-23 status quo (the two branches disagree, see Coupling); `1` = confined storativity at `gw_ss` above the ground in both branches; `2` = non-subgrid adopts the subgrid rule (capped at the ground, nothing stored above it) |
| `gw_ss` | real | 1.0e-4 | 1/m | confined/elastic storativity applied above the ground when `gw_storage_mode = 1`; unused otherwise |
| `gw_tolouter` | real | 1.0e-5 | m | outer-loop tolerance for rows carrying a lagged coupling term (an active seepage face, or a one-sided surface/aquifer exchange); the bulk surface tolerances (`si_tolouter`, `si_outer_frac`) cannot see these rows because the aquifer moves ~1e-4 m/step |
| `gw_zdrain` | real | -999.0 | m | uniform drain level; the drain is off unless this or `gw_zdrainfile` is set **and** `gw_cdrain > 0`. A uniform value drains *every* cell whose head exceeds it, sea bed and levees included |
| `gw_zdrainfile` | char | 'none' | - | spatial drain-level field (see Files); use `-999` outside the drained area (e.g. outside a polder's ditches) so nothing else drains |
| `gw_cdrain` | real | 0.0 | 1/s | drain conductance per unit area, same units as `gw_leakance` |

Two related keywords, not `gw_`-prefixed, that groundwater runs need to know about:

| keyword | type | default | unit | meaning |
|---|---|---|---|---|
| `semi_implicit` | logical | `.false.` | - | selects the coupled semi-implicit path (aquifer rows in the surface matrix) over the explicit sub-cycled path |
| `huthresh` | real | 0.05 | m | wet/dry film-depth threshold on the surface; a cell must hold this much water before it can pass any on, so it quantises groundwater-driven flooding (see "Two things to know", below) |

## Files

`gw_headfile`, `gw_rechargefile`, `gw_khfile`, `gw_syfile`, `gw_zbasefile`, `gw_zdrainfile` are
all a **flat `real*4` binary stream over active cells in internal order** — the same convention
`manningfile` uses (`sfincs_domain.f90:2005`). Make them by writing a `real*4` array of length
`np` (the number of active cells) in that cell order, unformatted stream access, no header.

`gw_khfile`, `gw_syfile` and `gw_zbasefile` are read **before** the initial head, because the
head is clamped to `gw_zbase` at initialization; reading the base field afterwards would clamp
against the uniform value instead and silently leave the head below the base wherever the field
sits higher. `gw_headfile` cells with `kcs == 2` (open boundary) are not solved for, so whatever
head they are given there is what they keep — that is how a fixed-head aquifer boundary is set.
A negative value in `gw_rechargefile` is a sink (e.g. a polder's ditch network, represented as
distributed drainage rather than a point sink). `gw_zdrainfile` follows the same convention;
`-999` marks a cell with no drain.

A conductivity `<= 0` or a specific yield `< 0` anywhere active stops the run rather than
silently decoupling a cell from its neighbours or breaking the symmetric-positive-definite
property the CG solver depends on.

## Coupling

- **Leakance exchange with the surface.** `Q = C * (max(zs, zref) - max(h, zref))`, positive
  from surface into aquifer, where `zref` is the bed (or, under subgrid, the cell's lowest
  subgrid point). This is the MODFLOW river/drain form: a dry surface above a deep water table
  exchanges nothing; a wet surface over a deep table infiltrates at a rate set by the ponded
  depth; a dry surface over a shallow table seeps water out. The part of this that is symmetric
  in `(zs, h)` is assembled implicitly in the matrix; the rest is lagged onto the right-hand
  side, so the pair reproduces the exact `Q` at outer convergence.
- **Seepage face.** The seepage ceiling is `max(ground, zs)`: the ground beneath standing water
  is saturated, so the table is allowed to sit as high as the free surface there. When the head
  exceeds that ceiling, `Q = cseep * (h - zceil)` carries the row's net inflow out of the aquifer
  as surface water rather than storing it; at `gw_seepage_fac = 1.0` that removes in one timestep
  exactly the volume that would otherwise have piled up above the ceiling.
- **Storage ceiling, and a known inconsistency between the two branches.** The two branches do not
  agree on what a cell stores under standing water. Without subgrid, storage is capped at
  `max(zb, zs)`, so a cell under a 1 m pond holds `Sy * 1 m * A` of extra water above its own bed,
  in a place where there is no soil. With subgrid, storage accrues only over the dry area
  `acell - awet(head)`, which is zero above the crest, so the same cell holds nothing. The
  non-subgrid convention also ejects `Sy * A * dzs` in a single step on a falling tide and refills
  it through leakance on the rise, which leaves the tidal-mean head under a submerged bed about
  0.09 m below mean sea level on the polder case. Removing the storage above the ground was
  attempted on 2026-09-02 and rejected: with the storage derivative at zero above the ground and
  the seepage face still switched off below `max(ground, zs)`, a cell in that band has neither
  storage nor seepage and its row loses its diagonal, which cost 3866 m3 of a 4320 m3 recharge on
  the sloping seepage-face case. Lowering the seepage ceiling to the ground closes that band but
  pins the head at the bed under a pond, which is wrong for a submerged aquifer. The convention is
  still open; see `gw_cases/RESULTS.md`.
- **`gw_storage_mode`.** Two candidates against that status quo are available behind this keyword,
  measured but not yet chosen between. `gw_storage_mode = 0` (default) is the status quo above,
  bit for bit. `gw_storage_mode = 1` gives the cell a small, real storage capacity above the
  ground instead of either extreme: below the ground it is `Sy` per metre as always, and above it
  the pond plays no part — the cell instead stores `gw_ss` (a confined/elastic storativity) per
  metre of confined aquifer thickness, in both branches alike, so they agree with each other for
  the first time. `gw_storage_mode = 2` makes the non-subgrid branch adopt the subgrid rule as-is:
  capped at the ground alone, nothing stored above it, no pond term. Both leave the seepage
  ceiling (`max(ground, zs)`) untouched — only what happens to storage between the ground and that
  ceiling changes.
  Measured 2026-09-22 over the whole case matrix (`plans/2026-09-22-groundwater-solver-defects-RESULTS.md`
  in the project folder): mode 2 fails the way the rejected 2026-09-02 attempt did (explicit and
  semi-implicit disagree by up to 1.8 m in the band with neither storage nor seepage); mode 1
  makes the two branches agree wherever storage was their difference (ceiling pond 0.232 m in
  both, seepslope 1374 m3 in both, the exchange case meets at 0.60 m in both) and puts the
  tidal-mean head under the sea on the polder case at mean sea level to a millimetre (mode 0:
  0.084 m below), because a falling tide no longer ejects `Sy * A * dzs`. Its cost is on the
  semi-implicit path, where the seepage-switch chatter that was confined to subgrid rows reaches
  the non-subgrid rows too (thousands of stalled outer steps per polder run, closure still
  within 0.01 %). Mode 1 is the default since 2026-09-23.
- **How the surface receives its share.** Without subgrid the level the semi-implicit solve
  returns is the state, and the surface row already carried rain, `qext`, the exchange and the
  seepage. With subgrid the state is the cell volume: the continuity re-integrates `z_volume`
  from the back-substituted fluxes and inverts the table, so the solve records the source
  volume it applied to each surface row (`si_qsrc`, m3) and the continuity adds it. From
  2026-09-02 to 2026-09-21 that hand-off was missing: a semi-implicit subgrid run received no
  rain, no `qext`, no exchange and no seepage on the surface (see Known limits). On the explicit
  path the aquifer sub-steps hand their net volume to the surface once per step (`gw_qsurf`).
- **Explicit path under a pond (subgrid).** The subgrid storage curve is flat between the highest
  pixel and the pond level, so a saturated submerged cell's head is set by the exchange alone.
  The explicit sub-step takes the symmetric part of the exchange implicitly in the cell's own
  head (Newton on `V(h) + C dt h = vol + C dt zs`, iterated to tolerance) and, for a cell that was
  saturated and submerged at the start of the sub-step, holds the head at the pond level and
  refills from the pond whatever the lateral flow took. Sub-stepping the flat band explicitly is
  otherwise a Jacobi sweep of an elliptic problem whose off-diagonal (`K b w / dx`) exceeds its
  diagonal (`leakance * awet`): the compound test case ran a crest/pond checkerboard under its
  pond with 65,000 m3 cycling through the exchange where the coupled solve moved 2,800. What the
  rule gives up is the rate-limited drawdown the coupled solve shows at a pond edge draining
  into the dry slope (6 cm on the compound case, 1e-4 m elsewhere).
- **Infiltration under a thin sheet depends on cell size with subgrid tables.** The exchange
  acts over the table's wet area, `C * awet(zs)`, and a rain sheet a few millimetres deep wets a
  small fraction of a cell whose pixels span a slope: a 20 m cell on the compound slope is about
  half as wet as a 10 m cell at the same sheet depth, so it infiltrates about half as much. That
  is the convention doing what it says, not a defect, but it makes the recharge under a rain
  sheet a function of the level of refinement. Measured: on the compound case the quadtree
  upslope head is 0.10 m below the regular grid's after 24 h; on the tidal island 29.3 % of the
  rain reaches the aquifer with subgrid against 39.8 % without, on either grid. Without subgrid
  the exchange acts over the whole cell as soon as it is wet, so it does not see the sheet depth
  at all. The one knob a user has is resolution over the infiltrating area: a finer level there
  brings the wet fraction, and the recharge, up towards the non-subgrid value.
- **Drain boundary.** `Q = gw_cdrain * A * max(h - gw_zdrain, 0)`, out of the aquifer and out of
  the model entirely (representing pumped ditch water), implicit on the aquifer diagonal,
  switched at the outer iterate the same way the seepage face is.

## Outputs and the water balance

Aquifer head is written as `gw_head` in `sfincs_map.nc` (`source/src/sfincs_ncoutput.f90`),
alongside `zs`, on both the regular and quadtree paths, whenever `gwflow` is on. Unlike `zs` it
is written wherever the cell is active, not masked below `huthresh` — a water table under dry
ground is the normal state of an aquifer and is usually the thing being examined.

With `gwflow` on, SFINCS prints a water balance to `sfincs.log` at the end of the run
(`gw_budget_report`, `sfincs_groundwater.f90`). Every term is signed as water **entering** the
aquifer:

| log line | meaning |
|---|---|
| `recharge` | volume from `gw_recharge` (negative where it is being used as a sink) |
| `exchange w/ surface` | net volume through the leakance exchange with the surface (negative = net exfiltration) |
| `lateral boundary` | flux across faces into cells that are not solved for (fixed-head boundary cells) |
| `ceiling seepage` | volume forced out through the seepage face |
| `drainage` | volume removed by the drain boundary (negative, or zero if the drain is off) |
| `Throughput` | sum of the *magnitudes* of every elementary contribution (face by face, step by step) — the scale the closure error is judged against, not the net inflow |
| `Closure error` | `(final storage - initial storage) - (recharge + exchange + lateral boundary + ceiling seepage + drainage)`, in m3 and as a percentage of throughput |
| `Storage change` | `final storage - initial storage`, which the five signed terms above must sum to |

Throughput, not the signed net, is the right denominator: a tidal aquifer takes water in on the
flood and gives it back on the ebb, so a signed boundary total near zero can sit next to an
enormous throughput, and closing against the small signed number manufactures a closure error
that is really just cancellation.

## Two things a user must know

1. **A uniform `gw_zdrain` drains every cell** whose head exceeds it — sea bed and levees
   included, not just the intended ditch network. A run with a real, spatially limited drain
   network needs `gw_zdrainfile`, set to the drain level under the drains and `-999` (no drain)
   everywhere else.
2. **`huthresh` sets the film, not the water.** A cell must hold `huthresh` of surface water
   before it can pass any on to a neighbour, so it puts a floor on any depth read off
   groundwater-driven flooding (seepage, drain outflow). Read depths only at `huthresh <= 0.01 m`;
   volumes and heads are reliable at any `huthresh`.

## Conceptual test cases

The suite's entry point is `D:/ClaudeProjects/sfincs-dev/gw_cases/README.md` (one exe folder,
one runner, one `reference.json` gate); the long-form records are the RESULTS files it points to.
Since 2026-09-21 the cases are three tiers, and every family in the first two is built in all
eight combinations of grid (regular / quadtree with one refinement), solver (explicit /
semi-implicit) and storage (without / with subgrid tables), `{uni|qt}_{exp|si}_{nosbg|sbg}`:

| tier | families | runs |
|---|---|---|
| `01_analytical` | Dupuit, Edelman (1 m and 0.1 m steps), Ferris, two-zone Dupuit | 40 |
| `02_coupling` | basin (infiltration to recharge), ceiling, exchange, sloping seepage face, compound (+ a `gwflow = 0` control) | 41 |
| `03_application` | polder (levee seepage, drains, pump failure), tidal island (2D, with and without subgrid and rain) | 16 |

Each family isolates something the others cannot see (the Dupuit parabola can be exact in
discharge while the transient runs 25% slow, because at steady state `theta` cancels and only
Edelman would catch it). Per family the reference is exact where one exists:

| family | what it constrains | accuracy reached (all 8 combinations unless noted) |
|---|---|---|
| Dupuit steady seepage | the lateral operator, both grids, the refinement transition | discharge +0.008% (regular) / +0.115% (quadtree); flux jump at the transition 0.000% |
| Edelman step response | the transient, and its linearisation (the error must scale tenfold from a 1 m to a 0.1 m step) | normalised max error 0.019 at 1 m, 0.0020 at 0.1 m |
| Ferris tidal wave | decay and phase together | fitted decay length 99.16 m vs 99.42 m exact (-0.27%), phase length +0.27% |
| Two-zone Dupuit | face transmissivity under a tenfold K contrast | flux +1.39% (regular) / +1.78% (quadtree, zone 2 coarse), interface step 0.000% |
| Closed basin | infiltration routed into the aquifer | recharge within 0.02% of rate x time x area on every mesh |
| Topographic ceiling | recharge a full aquifer cannot store becomes surface water | pond 0.1933 m (non-subgrid, `V/(A(1+Sy))`) and 0.2360 m (subgrid, `V/A` + 4 mm of table relief); explicit = semi-implicit to 0.2 mm within each |
| Exchange relaxation | the surface/aquifer coupling is symmetric | tau within 0.09% (non-subgrid); with subgrid the head meets the pond at once, by convention |
| Sloping seepage face | seeped water reaches the right surface cell | both paths close to 2e-4%; heads agree to 0.6 mm, ponded volumes 1145 / 2290 m3 (regular / quadtree, non-subgrid) |
| Compound | rain, runoff, infiltration, lateral flow and seepage in one step | 120 mm of rain against 119.98 mm stored (non-subgrid), 120.6 mm (subgrid) |

Across the matrix (`check_matrix.py`): explicit = semi-implicit to 1e-6 m on the analytical
families and 1e-4 m on the coupling families; subgrid = non-subgrid bit-identical on a flat bed;
quadtree = regular to 5 mm (Dupuit), 13 mm (two-zone), 8 mm (seepage face). The pairs where the
two storage conventions differ by design (ceiling, seepage face, compound, exchange with subgrid)
are reported and not gated. The quadtree strip is one level-1 row 2dx wide, so its volumes and
discharges are double the regular ones by construction.

The 2026-09-21 matrix exposed three defects, all fixed the same day and recorded in
`plans/2026-09-21-gw-cases-matrix-tier1-2-RESULTS.md`: the semi-implicit subgrid hand-off
(above), the semi-implicit budget booking seepage against the wrong ceiling with subgrid (17 m3,
1.3% on the ceiling case, while the state was right), and the explicit path's crest/pond
checkerboard under a pond (above). A fourth was in the test tables, not the model: a face table
whose sill lies below the neighbouring cell's lowest pixel reads zero depth on a face SFINCS
flags wet, and the predictor's friction term divided 0 by 0; the builder now applies the
hydromt convention, and the predictor floors the depth as its lookup-table branch already did.

Internal water-balance closure reaches machine precision (1e-5 to 1e-7 %) on every case whose
ceiling does not move during the run. Cases with an active seepage face or a moving ceiling
close to a few thousandths of a percent, because the seepage/ceiling switch is evaluated at the
previous outer iterate, a scheme choice.

## Known limits (unresolved as of this writing)

- (Resolved 2026-09-21.) Semi-implicit + subgrid received no surface sources at all between
  `df5b749` (2026-09-02) and this fix: rain, `qext`, the aquifer exchange and the seepage face
  were on the surface row of the pressure solve but never reached `z_volume`, which is the
  state on the subgrid path. `df5b749` had removed the continuity's rain term as a supposed
  double count and put nothing in its place; the crash it was chasing was a NaN from a face
  table with zero depth on a wet face (see the test-case section), not a doubled volume. Any
  semi-implicit subgrid run with rainfall, `qext` or `gwflow` from 2026-09-02 to 2026-09-21
  under-delivered its sources (a basin under 20 mm/h of rain infiltrated nothing); those runs,
  the Harvey semi-implicit runs with subgrid included, have to be redone.
- The two branches store different things under standing water (Sy per metre of pond above the
  bed without subgrid, nothing above the highest pixel with subgrid). The convention is still
  open; the explicit and semi-implicit paths agree with each other within each. The subgrid
  island budgets, which read -72 to -83 % of throughput on 8cadae6, close to 0.005 % since
  2026-09-21: that error was the missing surface hand-off, not the convention.
- The seepage face is asymmetric under standing water: head above `zs` is ejected in one step on
  a falling tide, but refill on a rising tide goes through leakance, so the tidal-mean head under
  the sea sits measurably below MSL.
- Upstream (not harmonic-mean) face transmissivity biases flux about 1.4% high at a tenfold K
  contrast; accepted, because a harmonic mean would freeze any dried cell permanently.
- Infiltration-to-recharge is verified for the constant-rate, Green-Ampt, Curve Number (`cna`)
  and Horton schemes (`gw_cases/basin/`, all within 0.02 %); Curve Number with recovery (`cnb`,
  `sefffile`) carries the same hookup untested.
- No confined aquifers, wells, or exponential conductivity with depth (all present in wflow).
- No separate infiltration/exfiltration leakance.
