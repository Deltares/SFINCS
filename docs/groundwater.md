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

From `D:/ClaudeProjects/sfincs-dev/gw_cases/RESULTS.md`. Each case isolates something the
others cannot see (e.g. the Dupuit parabola can be exact in discharge while the transient runs
25% slow, because at steady state `theta` cancels and only Edelman would catch it).

| case | what it constrains | accuracy reached |
|---|---|---|
| Dupuit steady seepage | the lateral operator, at steady state | head within 0.06%, discharge within 0.003%, uniform flux |
| Edelman step response | the transient — the time weighting | 0.26% of amplitude at delta/h0 = 1% |
| Ferris tidal wave | decay and phase together | amplitude within 0.4%, phase within 0.4%, fitted decay length 99.50 m vs 99.42 m exact |
| Exchange relaxation | that the surface/aquifer coupling is symmetric | tau within 0.09%, equilibrium level within 1.8e-3 m |
| Two-zone Dupuit | face transmissivity under heterogeneous K | heads within 0.25%, flux within 1.9%, interface step <1% |
| Topographic ceiling | what happens to recharge a full aquifer cannot store | both paths close to 0.05%, and agree on head and pond to 0.3 mm |
| Sloping seepage face | that seeped water reaches the right surface cell, within the step | both close; heads agree to 0.014 m, ponded volume to 2.2% |
| Closed basin | infiltration routed into the aquifer | recharge within 0.02% of the closed form |

The isolated ceiling case (`gw_cases/ceiling/`) ponds at **0.1933 m** rather than the
464 m3 / 2000 m2 = 0.232 m the raw excess implies, because the non-subgrid storage cap lets the
saturated cell hold part of that excess as groundwater above its own bed. The factor is
`1 / (1 + Sy)`. See the storage-ceiling note above.

Internal water-balance closure, once the balance was implemented to measure the same volumes the
solver already computed rather than reconstructing them in Python, reaches machine precision
(1e-5 to 1e-7 %) on every case whose ceiling does not move during the run. The three cases with an
active seepage face or a moving ceiling (`ceiling`, `seepslope`, `polder`) close to a few
thousandths of a percent instead, because the seepage/ceiling switch is evaluated at the previous
outer iterate — a scheme choice, not an error.

## Known limits (unresolved as of this writing)

- SI + subgrid + precipitation crashes at the first timestep, with or without `gwflow` — a
  surface-solver defect, not root-caused in source.
- The subgrid aquifer budget does not close on a moving tidal shore (island cases, -72% to -81%
  of throughput); the non-subgrid budget on the same kind of case closes to 1e-5%. Likely a
  budget-measurement gap rather than a physical leak, not yet instrumented further.
- The seepage face is asymmetric under standing water: head above `zs` is ejected in one step on
  a falling tide, but refill on a rising tide goes through leakance, so the tidal-mean head under
  the sea sits measurably below MSL.
- Upstream (not harmonic-mean) face transmissivity biases flux about 1.4% high at a tenfold K
  contrast; accepted, because a harmonic mean would freeze any dried cell permanently.
- Infiltration-to-recharge is verified only for the constant-rate and Green-Ampt schemes; Curve
  Number and Horton carry the same hookup untested.
- No confined aquifers, wells, or exponential conductivity with depth (all present in wflow).
- No separate infiltration/exfiltration leakance.
