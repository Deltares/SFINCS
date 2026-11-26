Numerical implementation
=======

Introduction 
----------------------

Based on Leijnse (2018)


Flow-related processes
-----

Shallow water equations
%%%%%


Numerical grid
%%%%%


Momentum equations
%%%%%



Advection
%%%%%



Continuity equation
%%%%%



Stability conditions
%%%%%


Wave-related processes
-----


Swash zone modelling approach
%%%%%




Weakly reflective generating-absorbing boundary condition
%%%%%



Wave generation
%%%%%
(Work in progress)



Wave-induced setup
%%%%%
(Work in progress)


SnapWave convergence criterion
%%%%%

SnapWave solves the stationary wave energy balance iteratively: each iteration consists of four sweeps through the grid, and each wave energy node is updated from its upwind nodes.
The dissipation ratio D/E of a node is evaluated from the energy of that node in the previous sweep and is under-relaxed between sweeps (snapwave_relax_factor_DoverE).
After every iteration a node is marked as converged, and skipped in the remaining iterations of that SnapWave call, when its energy change is smaller than snapwave_crit times the largest energy in the domain.
The call stops when the largest change satisfies this criterion, when 99% of the nodes have converged, or after snapwave_niter iterations.

Because the change of a node is compared with the largest energy in the domain, a node with little energy always looks converged.
This occurred at nodes without wave energy at the start of a SnapWave call: at the first call of a simulation, and at nodes that became wet since the previous call.
Such a node first receives waves without dissipation and overshoots, the large dissipation ratio that follows (kept high by the under-relaxation) then reduces its energy to almost zero, and in that state the node was marked as converged and kept for the rest of the call.
This resulted in too low wave heights behind shoals, reefs and dikes, and in nearshore areas that become wet during a storm.

From SFINCS v2.5.1 onwards, a node is therefore not marked as converged as long as it is 'cold' (its energy is below 1% of the energy that flows into it in the current sweep), or as long as its under-relaxed D/E differs more than 10% from its unrelaxed value.
A SnapWave Neumann boundary node (a lateral boundary node with a zero gradient in the wave field) copies the energy of its inner node, and is now marked as converged only when that inner node is.
At convergence the solution is unchanged; the change only prevents nodes from being skipped too early.
In the tested cases this costs about 20% more SnapWave iterations.


Other processes
-----



Infiltration
%%%%%



Precipitation
%%%%%


Discharge points
%%%%%


Wind forcing
%%%%%


Model limitations
-----



Computational efficiency
-----

