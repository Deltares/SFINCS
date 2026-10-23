! Non-hydrostatic pressure projection (single-layer Keller-box), GPU-ready.
!
! The depth-averaged momentum/continuity set is closed with a linearly-varying
! non-hydrostatic pressure pnh (bed value; depth mean pnh/2). Eliminating the
! provisional velocities gives a symmetric-positive-definite pressure operator
!
!     A = (dt/rho) G^T diag(hu) G  +  diag( kbfac dt / (rho H) ) ,
!
! built on the momentum stencil (G and the divergence are exact negative
! transposes), solved with a Chronopoulos-Gear (single-reduction) conjugate
! gradient, warm-started from the previous step. Works on regular and quadtree
! grids (refinement-transition cells are excluded from the nonh mask and stay
! hydrostatic, so every nonh face joins two same-level cells) and on projected
! and geographic (crsgeo) grids via per-face/per-row metric spacing.
!
! The solver machinery is designed for memory bandwidth and GPU offload:
!
!  1) ASSEMBLED 5-POINT STENCIL. The symmetric operator is assembled ONCE PER
!     SOLVE into four off-diagonal arrays (nh_aw/ae/as/an) with STATIC neighbour
!     indices (nh_nbw/nbe/nbs/nbn). The matvec is then a single loop:
!         w(i) = r(i) + aw*r(nbw) + ae*r(nbe) + as*r(nbs) + an*r(nbn)
!     instead of a per-iteration two-pass face form.
!
!  2) SYMMETRIC JACOBI SCALING BAKED IN. The system is solved as
!         (S A S) y = S b,   pnh = S y,   S = diag(1/sqrt(diag(A))),
!     so the scaled operator has UNIT diagonal: no preconditioner work in the
!     iteration, and <r,z> = <r,r> removes one reduction term.
!
!  3) BREAKING DIRICHLET BAKED IN. The row AND column of a breaking cell are
!     zeroed at assembly time (nh_sclo = nh_scl * activity, symmetric), the
!     scaled diagonal stays 1, and b=0 there -> pnh=0 is held exactly with NO
!     per-iteration branch (faces stay active for the flux correction).
!
!  4) GHOST ELEMENT 0. CG vectors run (0:nrows) with element 0 pinned at zero;
!     a missing neighbour has index 0 and a zero coefficient (nh_sclo(0)=0), so
!     the matvec has NO conditionals at all -- fully vectorizable on CPU,
!     coalesced on GPU.
!
!  5) TWO FUSED KERNELS PER ITERATION:
!         kernel B:  w = A r  fused with the dots  gam=<r,r>, del=<r,w>
!         kernel A:  p=r+beta*p ; s=w+beta*s ; y+=alpha*p ; r-=alpha*s
!     One host synchronisation point per iteration (the reduction). Dot products
!     accumulate in real*8 (the real*4 vectors are unchanged).
!
!  6) DUAL DIRECTIVES. Every kernel carries both !$omp (CPU, schedule(static))
!     and !$acc (GPU) directives, following the sfincs_momentum.f90 convention.
!     All module arrays are made device-resident at initialization (enter data);
!     per-step work then runs entirely on the device. The host copy of pnh is
!     refreshed only at map output times ('update host(pnh) if_present' in
!     sfincs_output, alongside the other output transfers).
!
! NOTE: the convergence test uses the SCALED residual norm ||S r|| / ||S b||
! (the natural Jacobi-CG norm); nonh_tol applies to that.
!
! On subgrid models zb is the EFFECTIVE bed zs - z_volume/area maintained every
! step in continuity (zb_effective), so all depth checks here see the volume-
! consistent depth.
!
! BOTTOM KINEMATIC CONDITION. The layer term abf in the face coefficients
! makes the transpose G^T (hu u) = u . grad(zb) - 0.5 h div(u), so the
! operator carries the bed-slope part of w_b implicitly and symmetrically
! (one-layer Keller-box constraint with the bottom condition eliminated,
! Escalante, Morales de Luna & Castro 2018, eq. 6). There is therefore no
! separate w_b: ws is the depth-mean vertical velocity (w_s + w_b)/2, updated
! by the pressure only, and a moving bed (nonh_movingbed) enters as the known
! source d(zb)/dt on the right-hand side. The source is the depth-filtered bed
! increment dzbext_filtered that update_bed_level (sfincs_domain) produces
! from dzbext (movingbed_filter, Gaussian std = factor x depth of the moving
! cells, because a single layer lacks the Kajiura attenuation of short bed
! features); the same filtered increment lifts the free surface there, so the
! mass and the pressure source of a moving bed are consistent on every grid
! type.
!
! BOUSSINESQ DISPERSION CORRECTION. Linearised on a flat bed the plain
! Keller-box scheme is, with the vertical factor kappa,
!
!     u_t = -g eta_x - (1/(2 rho)) p_x ,   h wbar_t = kappa p / rho ,
!     wbar = -(h/2) u_x  (constraint) ,     eta_t = -h u_x ,
!
! so p = -(rho h^2 / (2 kappa)) u_xt and
!
!     omega^2 = g h k^2 / (1 + (kh)^2 / (4 kappa)) ,
!
! a [0/2]-type approximant of Airy with a frequency cutoff omega_max =
! 2 sqrt(kappa g / h) above which the group velocity vanishes: short waves on
! a slope are trapped between the beach and their turning depth and never
! disperse away (LS2 landslide benchmark). The vertical momentum equation of
! the layer is therefore replaced by the Madsen & Sorensen (1992) type form
!
!     h wbar_t = (3/4) ( p - beta h^2 laplacian(p) ) / rho ,
!
! which with the constraint and the unchanged momentum equation gives, with
! X = (kh)^2, p = -rho h^2 k omega u / (2 (3/4) (1 + beta X)) and
!
!     omega^2 = g h k^2 (1 + beta X) / (1 + (beta + 1/3) X) .
!
! beta = 0 is the [0/2] Pade approximant 1/(1 + X/3) of tanh(kh)/kh (the
! strict Keller box has X/4 instead), beta = 1/15 the [2/2] approximant
! (1 + X/15)/(1 + 2X/5): c(k) within 1 % of Airy to kh ~ 3 and 5 % to kh ~ 6,
! no cutoff, group velocity -> sqrt(g h / 6) for kh -> inf. The default is
! beta = 0.04 (nonh_pade_beta; NTHMP test bed sfincs03: BP04-BP06 solitary
! waves and the LS1/LS2 landslides).
! In the discrete operator this is the diagonal with the Keller-box factor
! 3/4 and the face products of G^T diag(hu) G weighted by 1 + 3 beta: the
! extra part, (3/4) beta (dt/rho) G^T diag(hu) G, is a weighted graph
! Laplacian (symmetric, positive semi-definite, zero row sum), so the operator
! stays SPD and the CG solver, the momentum predictor and the flux correction
! are untouched. ws is taken from the constraint,
! ws = G^T (hu u^{n+1}) + d(zb)/dt, plus the residual A p - b of the solve:
! the sum is the vertical equation with the solved pressure, so the residual
! enters the next right-hand side and is cancelled by the next solve, while
! the pressure filter (nonh_filter) acts as applied. (The constraint value
! alone commits the residual of every solve to the flow: at the 1 % tolerance
! one impulsive solve left a 2 cm spurious wave on the LS1 landslide case.)
!
! Stability: the predictor stays hydrostatic and the whole correction is
! implicit, so one linear step is the forward-backward pair
!
!     u^{n+1}   = u^n - i dt g k eta^n (1 + beta X) / (1 + (beta + 1/3) X) ,
!     eta^{n+1} = eta^n - i dt h k u^{n+1} ,
!
! stable for dt omega_eff <= 2 with omega_eff <= k sqrt(g h): the hydrostatic
! CFL limit, no extra restriction. The same relation follows from an explicit
! term beta g h^2 d/dx laplacian(eta) in the momentum predictor, but that
! predictor increment is (1 + beta X) times the hydrostatic one, ~10^3 at the
! grid scale for h/dx = 75, which the projection must cancel to within the CG
! tolerance; that form blew up on the LS2 static-pulse test (max |zs| 2.8 m
! from a 2 cm bed pulse, also with nonh_tol = 1e-5) and was dropped.
!
module sfincs_nonhydrostatic
   !
   ! Headers auto-generated, ask main developers in case of doubts
   !
   ! Non-hydrostatic pressure correction (regular and quadtree grids). See the notes above
   ! the module statement.
   !
   ! Subroutines / functions:
   !
   !   initialize_nonhydrostatic()
   !     Set up the indices and the sparse matrix structure of the non-hydrostatic solver.
   !     Called from: sfincs_initialize (sfincs_lib)
   !
   !   compute_nonhydrostatic(dt)
   !     Solve for the non-hydrostatic pressure and correct the fluxes and velocities.
   !     Called from: sfincs_update (sfincs_lib)
   !
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
   ! All indications, logos of, and references to, "SFINCS"
   ! remain the property of Stichting Deltares. All rights reserved.
   !
   ! ------------ Copyright notice SFINCS ------------
   !   
   implicit none
   !
   ! Input parameters (sfincs.inp, read in sfincs_input)
   !
   real*4    :: nonh_fnudge      ! fraction of the pressure correction applied to uv/q
   real*4    :: nonh_tstop       ! time (s, absolute after input processing) to stop applying nonh corrections
   integer   :: nonh_itermax     ! max CG iterations per solve
   real*4    :: nonh_filter      ! global 2dx pnh filter weight; 0 = off
   real*4    :: nonh_brfr        ! Froude breaking onset (NEOWAVE); 0 = off (no breaking)
   real*4    :: nonh_treform     ! breaking reformation time scale (s): released cells recover gradually, brfac += dt/treform; 0 = instant
   real*4    :: nonh_tol         ! CG relative tolerance (scaled residual norm)
   logical   :: nonh_movingbed   ! add d(zb)/dt (= dzbext_filtered/dt) as the moving-bed source of the layer constraint, so a moving seafloor radiates a depth-filtered (Kajiura-like) surface response
   real*4    :: nonh_pade_beta   ! coefficient beta of the Boussinesq dispersion correction (default 0.04; 1/15 = Pade [2/2] of Airy, 0 = [0/2])
   !
   ! Keller-box factor (3/4) and weight of the face products of the operator
   ! (1 + 3 beta) of the Boussinesq correction (set in initialize_nonhydrostatic)
   !
   real*4    :: nh_kbfac
   real*4    :: nh_pdf
   !
   ! Non-hydrostatic cell mask (filled from the quadtree file in sfincs_domain)
   !
   integer*1, dimension(:), allocatable :: mask_nonh
   !
   integer, dimension(:), allocatable   :: nm_index_of_row
   integer, dimension(:), allocatable   :: row_index_of_nm
   integer, dimension(:), allocatable   :: uv_index_of_nhuv
   !
   real*4, dimension(:),   allocatable  :: pnh        ! non-hydrostatic bed pressure   (nrows)
   real*4, dimension(:),   allocatable  :: ws         ! depth-mean vertical velocity   (nrows)
   !
   real*4,  dimension(:), allocatable   :: Dnm        ! layer depth per row            (nrows)
   !
   real*4,  dimension(:), allocatable   :: nh_dvert   ! vertical-accel diagonal term   (nrows)
   real*4,  dimension(:), allocatable   :: nh_brfac   ! breaking STATE 1=full nonh .. 0=hydrostatic (nrows; hysteresis memory)
   real*4,  dimension(:), allocatable   :: nh_bract   ! breaking ACTIVITY used by the pressure operator (nrows; = nh_brfac)
   !
   integer, dimension(:), allocatable   :: nh_faceuv  ! full uv index of each nh face  (nhuv)
   integer, dimension(:), allocatable   :: nh_faceL   ! row index of left/bottom cell  (nhuv, 0 = boundary)
   integer, dimension(:), allocatable   :: nh_faceR   ! row index of right/top cell    (nhuv, 0 = boundary)
   real*4,  dimension(:), allocatable   :: nh_cf      ! 0.5/dx static gradient weight  (nhuv)
   real*4,  dimension(:), allocatable   :: nh_cR      ! gradient coeff of p(R) per face(nhuv)
   real*4,  dimension(:), allocatable   :: nh_cL      ! gradient coeff of p(L) per face(nhuv)
   real*4,  dimension(:), allocatable   :: nh_hu      ! water depth at face this step  (nhuv)
   !
   ! Per-row incident-face map. For each row, the (up to 4) non-hydrostatic faces
   ! touching it, and implicitly which side the row is on:
   !   slot 1 = left  face -> row is the R-cell -> use nh_cR ; "other" cell = nh_faceL
   !   slot 2 = right face -> row is the L-cell -> use nh_cL ; "other" cell = nh_faceR
   !   slot 3 = below face -> row is the R-cell -> use nh_cR ; "other" cell = nh_faceL
   !   slot 4 = above face -> row is the L-cell -> use nh_cL ; "other" cell = nh_faceR
   ! 0 = no face on that side.
   !
   integer, dimension(:,:), allocatable :: nh_cellface  ! (4, nrows)
   !
   ! Static 5-point neighbour map: row index of the West/East/South/North
   ! neighbour of each row (0 = no nonh neighbour on that side -> ghost row 0).
   !
   integer, dimension(:), allocatable   :: nh_nbw, nh_nbe, nh_nbs, nh_nbn      ! (nrows)
   !
   ! Assembled, symmetrically-scaled off-diagonals of the operator (unit diagonal),
   ! rebuilt every solve (wet/dry, depths and breaking change each step).
   !
   real*4, dimension(:), allocatable    :: nh_aw, nh_ae, nh_as, nh_an          ! (nrows)
   real*4, dimension(:), allocatable    :: nh_scl     ! 1/sqrt(diag(A))          (nrows)
   real*4, dimension(:), allocatable    :: nh_sclo    ! scl * breaking-activity  (0:nrows), ghost 0
   !
   ! CG vectors, ghost element 0 pinned at zero (never written by any kernel).
   !
   real*4, dimension(:), allocatable    :: nh_y       ! scaled pressure iterate  (0:nrows)
   real*4, dimension(:), allocatable    :: nh_r       ! residual                 (0:nrows)
   real*4, dimension(:), allocatable    :: nh_p       ! search direction         (0:nrows)
   real*4, dimension(:), allocatable    :: nh_s       ! s = A p (recurrence)     (0:nrows)
   real*4, dimension(:), allocatable    :: nh_w       ! w = A r                  (0:nrows)
   real*4, dimension(:), allocatable    :: nh_b       ! scaled right-hand side   (0:nrows)
   !
   ! Moving-bed rate d(zb)/dt = dzbext_filtered/dt of this step per row: the known
   ! source of the layer-consistent constraint (nonh_movingbed).
   ! Bed increments below nh_dzb_ulps float32 ulps of the bed level (of 1 m for
   ! |zb| < 1 m) are round-off of the caller's bed difference, not motion; they are
   ! dropped by update_bed_level (sfincs_domain) when it builds dzbext_filtered.
   !
   real*4, parameter :: nh_dzb_ulps = 4.0
   !
   real*4, dimension(:), allocatable    :: nh_dzbdt   ! d(zb)/dt this step       (nrows)
   !
   real*4    :: huthresh_nh
   !
   integer   :: nrows
   integer   :: nhuv
   !
   ! CG iteration diagnostics (cumulative over the run; reported in the timing summary)
   !
   integer*8 :: nh_iter_total = 0     ! sum of CG iterations over all solves
   integer*8 :: nh_solve_count = 0    ! number of CG solves (active time steps)
   integer   :: nh_iter_max   = 0     ! most CG iterations in any single solve
   !
contains
   !
   subroutine initialize_nonhydrostatic()
   !
   ! Builds the row<->cell maps, the nonh face maps and the static neighbour
   ! map, allocates the solver work arrays and pushes everything
   ! the per-step kernels touch onto the GPU.
   !
   ! Called from: sfincs_initialize (sfincs_lib)
   !
   use sfincs_data
   use sfincs_log
   !
   implicit none
   !
   integer :: nm, irow, ip, inhuv, rl, rr, i, f
   !
   allocate(row_index_of_nm(np))
   !
   row_index_of_nm = 0
   !
   ! Count number of rows (= nr cols) in matrix
   !
   nrows = 0
   !
   do nm = 1, np
      !
      if (mask_nonh(nm) == 1) then
         !
         nrows = nrows + 1
         !
      endif
      !
   enddo
   !
   allocate(pnh(nrows))
   allocate(ws(nrows))
   allocate(Dnm(nrows))
   allocate(nm_index_of_row(nrows))
   allocate(nh_dvert(nrows))
   allocate(nh_brfac(nrows))
   allocate(nh_bract(nrows))
   !
   pnh      = 0.0
   ws       = 0.0
   Dnm      = 0.0
   nh_dvert = 0.0
   nh_brfac = 1.0
   nh_bract = 1.0
   nm_index_of_row = 0
   !
   huthresh_nh = max(huthresh, 0.01)
   !
   ! Map row indices to nm indices and vice versa
   !
   irow = 0
   !
   do nm = 1, np
      if (mask_nonh(nm) == 1) then
         irow = irow + 1
         row_index_of_nm(nm) = irow
         nm_index_of_row(irow) = nm
      endif
   enddo
   !
   ! Find velocity points needed for nh computations: any uv point that touches
   ! a nonh cell. First count, then store.
   !
   nhuv = 0
   !
   do ip = 1, npuv
      !
      if (mask_nonh(uv_index_z_nm(ip)) == 1 .or. mask_nonh(uv_index_z_nmu(ip)) == 1) then
         !
         nhuv = nhuv + 1
         !
      endif
      !
   enddo
   !
   allocate(uv_index_of_nhuv(nhuv))
   uv_index_of_nhuv = 0
   !
   inhuv = 0
   !
   do ip = 1, npuv
      !
      if (mask_nonh(uv_index_z_nm(ip)) == 1 .or. mask_nonh(uv_index_z_nmu(ip)) == 1) then
         !
         inhuv = inhuv + 1
         uv_index_of_nhuv(inhuv) = ip
         !
      endif
      !
   enddo
   !
   ! Face arrays
   !
   allocate(nh_faceuv(nhuv))
   allocate(nh_faceL(nhuv))
   allocate(nh_faceR(nhuv))
   allocate(nh_cf(nhuv))
   allocate(nh_cR(nhuv))
   allocate(nh_cL(nhuv))
   allocate(nh_hu(nhuv))
   nh_faceuv = 0
   nh_faceL  = 0
   nh_faceR  = 0
   nh_cf     = 0.0
   nh_cR     = 0.0
   nh_cL     = 0.0
   nh_hu     = 0.0
   !
   ! Static per-face data. nh_cf = 0.5 / dx is the (depth-averaged) gradient
   ! weight; the factor 0.5 reflects that the linearly-varying non-hydrostatic
   ! pressure has a depth mean of p_bed/2. The spacing is taken PER FACE so it is
   ! correct across quadtree refinement levels and for geographic (crsgeo) grids:
   ! x-faces use the per-uv-point metre spacing dxminv (crsgeo) or the per-level
   ! dxrinv (projected), y-faces the per-level dyrinv (uniform in both modes).
   ! nh_cf is a single value per face, used identically by the gradient and its
   ! transpose divergence, so the operator A = G^T diag(hu) G stays symmetric
   ! regardless of spacing variation.
   !
   do ip = 1, nhuv
      !
      inhuv = uv_index_of_nhuv(ip)
      nh_faceuv(ip) = inhuv
      nh_faceL(ip)  = row_index_of_nm(uv_index_z_nm(inhuv))
      nh_faceR(ip)  = row_index_of_nm(uv_index_z_nmu(inhuv))
      !
      if (uv_flags_dir(inhuv) == 0) then
         if (crsgeo) then
            nh_cf(ip) = 0.5 * dxminv(inhuv)
         else
            nh_cf(ip) = 0.5 * dxrinv(uv_flags_iref(inhuv))
         endif
      else
         nh_cf(ip) = 0.5 * dyrinv(uv_flags_iref(inhuv))
      endif
      !
   enddo
   !
   ! Per-row incident-face map (see header). Walk the faces once: each face is
   ! the RIGHT face of its L-cell and the LEFT face of its R-cell (x), or the
   ! ABOVE face of its L-cell and the BELOW face of its R-cell (y).
   !
   allocate(nh_cellface(4, nrows))
   nh_cellface = 0
   do ip = 1, nhuv
      inhuv = nh_faceuv(ip)
      rl    = nh_faceL(ip)
      rr    = nh_faceR(ip)
      if (uv_flags_dir(inhuv) == 0) then
         if (rl > 0) nh_cellface(2, rl) = ip   ! face is the L-cell's RIGHT face
         if (rr > 0) nh_cellface(1, rr) = ip   ! face is the R-cell's LEFT face
      else
         if (rl > 0) nh_cellface(4, rl) = ip   ! face is the L-cell's ABOVE face
         if (rr > 0) nh_cellface(3, rr) = ip   ! face is the R-cell's BELOW face
      endif
   enddo
   !
   ! Static 5-point neighbour map: neighbour row of each row through its (up to 4)
   ! incident nonh faces. nh_faceL/R are 0 on a domain/open-boundary side ->
   ! neighbour 0 -> ghost element.
   !
   allocate(nh_nbw(nrows))
   allocate(nh_nbe(nrows))
   allocate(nh_nbs(nrows))
   allocate(nh_nbn(nrows))
   !
   do i = 1, nrows
      nh_nbw(i) = 0 ; f = nh_cellface(1, i) ; if (f > 0) nh_nbw(i) = nh_faceL(f)
      nh_nbe(i) = 0 ; f = nh_cellface(2, i) ; if (f > 0) nh_nbe(i) = nh_faceR(f)
      nh_nbs(i) = 0 ; f = nh_cellface(3, i) ; if (f > 0) nh_nbs(i) = nh_faceL(f)
      nh_nbn(i) = 0 ; f = nh_cellface(4, i) ; if (f > 0) nh_nbn(i) = nh_faceR(f)
   enddo
   !
   ! Solver work arrays: assembled stencil, scaling and CG vectors (ghost
   ! element 0). First-touch initialization (NUMA page placement matches the
   ! static loop partition the per-step kernels use); ghost stays 0 forever.
   !
   allocate(nh_aw(nrows))
   allocate(nh_ae(nrows))
   allocate(nh_as(nrows))
   allocate(nh_an(nrows))
   allocate(nh_scl(nrows))
   allocate(nh_sclo(0:nrows))
   allocate(nh_y(0:nrows))
   allocate(nh_r(0:nrows))
   allocate(nh_p(0:nrows))
   allocate(nh_s(0:nrows))
   allocate(nh_w(0:nrows))
   allocate(nh_b(0:nrows))
   allocate(nh_dzbdt(nrows))
   !
   !$omp parallel do schedule ( static ) private ( i )
   do i = 1, nrows
      nh_dzbdt(i) = 0.0
      nh_aw(i)   = 0.0
      nh_ae(i)   = 0.0
      nh_as(i)   = 0.0
      nh_an(i)   = 0.0
      nh_scl(i)  = 1.0
      nh_sclo(i) = 1.0
      nh_y(i)    = 0.0
      nh_r(i)    = 0.0
      nh_p(i)    = 0.0
      nh_s(i)    = 0.0
      nh_w(i)    = 0.0
      nh_b(i)    = 0.0
   enddo
   !$omp end parallel do
   nh_sclo(0) = 0.0
   nh_y(0)    = 0.0
   nh_r(0)    = 0.0
   nh_p(0)    = 0.0
   nh_s(0)    = 0.0
   nh_w(0)    = 0.0
   nh_b(0)    = 0.0
   !
   ! Device residency for everything the per-step kernels read or write. The
   ! global sfincs_data arrays (zs, zb, uv, q, kfuv, kcs, z_index_uv_*,
   ! uv_index_z_*) are already on the device via initialize_openacc.
   !
   !$acc enter data copyin( nh_nbw, nh_nbe, nh_nbs, nh_nbn, &
   !$acc                    nh_aw, nh_ae, nh_as, nh_an, nh_scl, nh_sclo, &
   !$acc                    nh_y, nh_r, nh_p, nh_s, nh_w, nh_b, nh_dzbdt, &
   !$acc                    pnh, ws, Dnm, nm_index_of_row, &
   !$acc                    nh_dvert, nh_brfac, nh_bract, &
   !$acc                    nh_faceuv, nh_faceL, nh_faceR, nh_cellface, &
   !$acc                    nh_cf, nh_cR, nh_cL, nh_hu )
   !
   write(logstr,'(a,i0,a,i0,a)')'Non-hydrostatic solver (assembled-stencil CG): ', nrows, ' cells, ', nhuv, ' faces'
   call write_log(logstr, 0)
   !
   ! Keller-box factor and face weight of the Boussinesq correction (see header)
   !
   nh_kbfac = 0.75
   nh_pdf   = 1.0 + 3.0 * nonh_pade_beta
   !
   write(logstr,'(a,f8.5,a,f6.3,a,f6.3)')'Non-hydrostatic Boussinesq correction: beta = ', nonh_pade_beta, ', Keller-box factor ', nh_kbfac, ', face weight ', nh_pdf
   call write_log(logstr, 0)
   !
   end subroutine


   subroutine compute_nonhydrostatic(dt)
   !
   ! Per-step non-hydrostatic pressure projection (see module header). Every
   ! numbered kernel below is one !$omp / !$acc loop.
   !
   ! Input:
   !   dt : time step [s]
   !
   ! Called from: sfincs_update (sfincs_lib)
   !
   use sfincs_data
   use sfincs_timers
   !
   implicit none
   !
   real*4    :: dt
   !
   integer   :: i, ip, ipuv, nm, nmu, iuvl, iuvr, j, f, iter, cnt
   real*4    :: dtrho, kbfac, pdf, ta, hu, Dnm1, Dnmu, abf, tf, dval, sq, braw, accv, sumn, brrec, act
   real*4    :: qr, ql, dzdt, gf, pL, pR, unh
   real*4    :: alpha, beta, fmb
   real*8    :: gam8, del8, bn8, gamold8, alpha8, beta8, pap8, bnorm8
   logical   :: bndok
   !
   if (nrows == 0) return
   !
   call timer_start('non-hydrostatic')
   !
   ! Keller-box vertical factor 3/4 and face-product weight pdf = 1 + 3 beta of
   ! the Boussinesq correction (see header).
   !
   kbfac = nh_kbfac
   pdf   = nh_pdf
   dtrho = dt / rhow
   !
   ! breaking recovery per step: instant unless a reformation time is set
   !
   brrec = 1.0
   if (nonh_treform > 0.0) brrec = dt / nonh_treform
   !
   ! 1) Layer depth per row and the vertical-acceleration diagonal.
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i, nm )
   do i = 1, nrows
      nm = nm_index_of_row(i)
      Dnm(i)      = max(zs(nm) - zb(nm), huthresh_nh)
      nh_dvert(i) = kbfac * dt / (rhow * Dnm(i))
   enddo
   !$omp end parallel do
   !
   ! 1b) Breaking criterion -> nh_brfac (1 = full nonh, <1 = hydrostatic pnh=0).
   !     Froude (NEOWAVE; Yamazaki, Kowalik & Cheung 2009, eqs 33/34) when
   !     nonh_brfr > 0; no breaking otherwise.
   !
   if (nonh_brfr > 0.0) then
      !
      ! Froude criterion: a cell is set hydrostatic (pnh = 0) when the local flow
      ! Froude number Fr = |U| / sqrt(g D) exceeds the onset nonh_brfr (~0.5), and
      ! reactivated when Fr drops below 0.3*nonh_brfr -- NEOWAVE's hysteresis ratio.
      ! Purely LOCAL (no smoothing / neighbour-spread / hole-filling needed) and
      ! tracks a bore through its life.
      !
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i, nm, iuvl, iuvr, ql, qr, dzdt )
      do i = 1, nrows
         nm = nm_index_of_row(i)
         ql = 0.0 ; qr = 0.0
         iuvl = z_index_uv_md(nm) ; iuvr = z_index_uv_mu(nm)
         if (iuvl > 0) ql = uv(iuvl)
         if (iuvr > 0) qr = uv(iuvr)
         dzdt = (0.5 * (ql + qr))**2                       ! u_cell^2
         if (nmax > 1) then
            ql = 0.0 ; qr = 0.0
            iuvl = z_index_uv_nd(nm) ; iuvr = z_index_uv_nu(nm)
            if (iuvl > 0) ql = uv(iuvl)
            if (iuvr > 0) qr = uv(iuvr)
            dzdt = dzdt + (0.5 * (ql + qr))**2             ! + v_cell^2
         endif
         dzdt = sqrt(dzdt) / sqrt(g * Dnm(i))              ! Fr = |U| / sqrt(g D)
         if (nh_brfac(i) >= 1.0) then                      ! was not breaking
            if (dzdt > nonh_brfr) nh_brfac(i) = 0.0          ! onset
         else                                              ! was breaking -> release only when slow
            if (dzdt < 0.3 * nonh_brfr) then
               nh_brfac(i) = min(1.0, nh_brfac(i) + brrec)  ! recover (gradually if nonh_treform > 0)
            else
               nh_brfac(i) = 0.0
            endif
         endif
      enddo
      !$omp end parallel do
      !
   else
      !
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i )
      do i = 1, nrows
         nh_brfac(i) = 1.0
      enddo
      !$omp end parallel do
      !
   endif
   !
   ! 1c) Breaking ACTIVITY for the pressure operator: nh_bract = nh_brfac. The
   !     hysteresis STATE stays in nh_brfac.
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i )
   do i = 1, nrows
      nh_bract(i) = nh_brfac(i)
   enddo
   !$omp end parallel do
   !
   ! 2) Per-face gradient coefficients (wet/dry, layer term abf, conveyance depth).
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( ip, ipuv, nm, nmu, hu, Dnm1, Dnmu, abf, bndok )
   do ip = 1, nhuv
      ipuv = nh_faceuv(ip)
      nm   = uv_index_z_nm(ipuv)
      nmu  = uv_index_z_nmu(ipuv)
      nh_cR(ip) = 0.0
      nh_cL(ip) = 0.0
      nh_hu(ip) = 0.0
      ! Open (kcs 2/3) boundary faces are skipped (handled hydrostatically);
      ! closed walls / Neumann keep kfuv = 0 -> correct solid-wall BC.
      bndok = (kcs(nm) == 1 .and. kcs(nmu) == 1)
      if (kfuv(ipuv) == 1 .and. bndok) then
         if (zs(nm) > zb(nm) + huthresh_nh .and. zs(nmu) > zb(nmu) + huthresh_nh) then
            hu = max(zs(nm), zs(nmu)) - 0.5 * (zb(nm) + zb(nmu))
            if (hu > huthresh_nh) then
               Dnm1 = max(zs(nm)  - zb(nm),  huthresh_nh)
               Dnmu = max(zs(nmu) - zb(nmu), huthresh_nh)
               abf  = ( (zs(nmu) + zb(nmu)) - (zs(nm) + zb(nm)) ) / (Dnm1 + Dnmu)
               nh_cR(ip) = nh_cf(ip) * (1.0 + abf)
               nh_cL(ip) = nh_cf(ip) * (abf - 1.0)
               ! Conveyance depth q/uv keeps the corrected q and uv consistent
               if (uv(ipuv) /= 0.0) then
                  nh_hu(ip) = q(ipuv) / uv(ipuv)
               else
                  nh_hu(ip) = hu
               endif
            endif
         endif
      endif
   enddo
   !$omp end parallel do
   !
   ! 2b) Moving-bed source of the layer-consistent constraint: d(zb)/dt =
   !     dzbext_filtered/dt of this step per wet row. dzbext_filtered
   !     is built by update_bed_level (sfincs_domain) from the caller's dzbext:
   !     round-off increments dropped (a few float32 ulps of the bed level; divided
   !     by the 1e-6 s start-up step such noise is a 1e-2 m/s bed velocity whose
   !     impulse pressure the CG does not clean up, a static slide then radiated a
   !     6 cm wave train on LS1), wet cells only, depth-filtered (movingbed_filter).
   !     Host loop (dzbext_filtered is a host array), pushed to the device. fmb
   !     switches the source on in the assembly; nh_dzbdt is only read when fmb = 1.
   !     NOTE (GPU builds): this loop reads the HOST zs and zb, which are not
   !     refreshed from the device every step; on a GPU build they can be stale.
   !
   fmb = 0.0
   !
   if (nonh_movingbed .and. use_dzbext .and. dt > 0.0) then
      !
      if (allocated(dzbext_filtered)) then
         !
         fmb = 1.0
         !
         !$omp parallel do schedule ( static ) private ( i, nm )
         do i = 1, nrows
            !
            nm          = nm_index_of_row(i)
            nh_dzbdt(i) = 0.0
            !
            if (zs(nm) - zb(nm) > huthresh_nh) then
               !
               nh_dzbdt(i) = dzbext_filtered(nm) / dt
               !
            endif
            !
         enddo
         !$omp end parallel do
         !
         !$acc update device(nh_dzbdt)
         !
      endif
      !
   endif
   !
   ! 3) ASSEMBLY (once per solve). Per row: raw off-diagonals a_ij = dtrho*hu*cR*cL
   !    (one product per incident face, same both ways -> symmetric), the diagonal
   !    dval = dvert + sum dtrho*hu*coef^2, the scaling scl = 1/sqrt(dval), the
   !    breaking-masked scaling sclo, the SCALED right-hand side
   !    b = scl * [ G^T (hu u*) - ws + d(zb)/dt ], and the warm start
   !    y0 = pnh * sqrt(dval) (scaled previous pressure; 0 at breaking cells).
   !    The face product ta = dtrho*hu*cR*cL is weighted by pdf = 1 + 3 beta and
   !    the diagonal gets the matching -(pdf - 1) ta, i.e. the added part is a
   !    weighted graph Laplacian (zero row sum, SPD kept).
   !
   !    Layer-consistent bottom kinematic condition: with the layer term abf the
   !    transpose G^T (hu u) already equals u . grad(zb) - 0.5 h div(u), i.e. it
   !    contains the bed-slope part of w_b implicitly and symmetrically
   !    (Keller-box constraint with the bottom condition eliminated, Escalante
   !    et al. 2018 eq. 6). ws holds the depth-mean vertical velocity
   !    (w_s + w_b)/2 and the constraint is G^T (hu u^{n+1}) = ws^{n+1} -
   !    d(zb)/dt, so b = scl * [ G^T (hu u*) - ws + d(zb)/dt ] with no separate
   !    w_b.
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i, dval, braw, f, tf, ta, sq, act )
   do i = 1, nrows
      dval = nh_dvert(i)
      !
      braw = - ws(i) + fmb * nh_dzbdt(i)
      !
      f = nh_cellface(1, i)
      if (f > 0) then
         tf = dtrho * nh_hu(f)
         ta = tf * nh_cR(f) * nh_cL(f)
         nh_aw(i) = ta * pdf
         dval = dval + tf * nh_cR(f)**2 - (pdf - 1.0) * ta
         braw = braw + nh_cR(f) * nh_hu(f) * uv(nh_faceuv(f))
      else
         nh_aw(i) = 0.0
      endif
      f = nh_cellface(2, i)
      if (f > 0) then
         tf = dtrho * nh_hu(f)
         ta = tf * nh_cL(f) * nh_cR(f)
         nh_ae(i) = ta * pdf
         dval = dval + tf * nh_cL(f)**2 - (pdf - 1.0) * ta
         braw = braw + nh_cL(f) * nh_hu(f) * uv(nh_faceuv(f))
      else
         nh_ae(i) = 0.0
      endif
      f = nh_cellface(3, i)
      if (f > 0) then
         tf = dtrho * nh_hu(f)
         ta = tf * nh_cR(f) * nh_cL(f)
         nh_as(i) = ta * pdf
         dval = dval + tf * nh_cR(f)**2 - (pdf - 1.0) * ta
         braw = braw + nh_cR(f) * nh_hu(f) * uv(nh_faceuv(f))
      else
         nh_as(i) = 0.0
      endif
      f = nh_cellface(4, i)
      if (f > 0) then
         tf = dtrho * nh_hu(f)
         ta = tf * nh_cL(f) * nh_cR(f)
         nh_an(i) = ta * pdf
         dval = dval + tf * nh_cL(f)**2 - (pdf - 1.0) * ta
         braw = braw + nh_cL(f) * nh_hu(f) * uv(nh_faceuv(f))
      else
         nh_an(i) = 0.0
      endif
      sq = sqrt(dval)              ! dval >= dvert > 0 always (Dnm capped)
      nh_scl(i) = 1.0 / sq
      !
      ! Breaking activity act in [0,1]: act = 1 full nonh; act = 0 holds pnh = 0
      ! (row+column zeroed via sclo, unit scaled diagonal). Fractional act (from
      ! nonh_treform) ramps the pressure back in smoothly; the operator stays
      ! SPD (T A T + diag(1 - act^2), act <= 1).
      !
      act = nh_bract(i)
      nh_sclo(i) = nh_scl(i) * act
      nh_b(i)    = braw * nh_scl(i) * act
      nh_y(i)    = pnh(i) * sq * act
   enddo
   !$omp end parallel do
   !
   ! 3b) Scale the off-diagonals: a_ij <- a_ij * sclo(i) * sclo(j). Needs the
   !     completed sclo of the NEIGHBOUR -> separate kernel. sclo(0) = 0 zeroes
   !     coefficients to missing neighbours (the scaled diagonal is exactly 1).
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i )
   do i = 1, nrows
      nh_aw(i) = nh_aw(i) * nh_sclo(i) * nh_sclo(nh_nbw(i))
      nh_ae(i) = nh_ae(i) * nh_sclo(i) * nh_sclo(nh_nbe(i))
      nh_as(i) = nh_as(i) * nh_sclo(i) * nh_sclo(nh_nbs(i))
      nh_an(i) = nh_an(i) * nh_sclo(i) * nh_sclo(nh_nbn(i))
   enddo
   !$omp end parallel do
   !
   ! 4) CG solve of the scaled system (unit diagonal -> no preconditioner work).
   !    Chronopoulos-Gear single-reduction CG: per iteration ONE fused
   !    matvec+dots kernel and ONE fused vector-update kernel.
   !
   ! 4a) initial residual: w = A y0 (branch-free ghost matvec), r = b - w
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i )
   do i = 1, nrows
      nh_w(i) = nh_y(i) + nh_aw(i) * nh_y(nh_nbw(i)) &
                        + nh_ae(i) * nh_y(nh_nbe(i)) &
                        + nh_as(i) * nh_y(nh_nbs(i)) &
                        + nh_an(i) * nh_y(nh_nbn(i))
   enddo
   !$omp end parallel do
   !
   bn8 = 0.0d0
   !$acc parallel loop default(present) reduction(+:bn8)
   !$omp parallel do schedule ( static ) private ( i ) reduction(+:bn8)
   do i = 1, nrows
      nh_r(i) = nh_b(i) - nh_w(i)
      nh_p(i) = 0.0
      nh_s(i) = 0.0
      bn8 = bn8 + nh_b(i) * nh_b(i)
   enddo
   !$omp end parallel do
   !
   bnorm8 = sqrt(bn8)
   if (bnorm8 <= 0.0d0) bnorm8 = 1.0d0
   !
   ! 4b) iterate. Convergence on the scaled residual: ||S r|| / ||S b|| < nonh_tol.
   !
   iter    = 0
   gamold8 = 1.0d0
   alpha8  = 1.0d0
   !
   do
      !
      ! kernel B: w = A r fused with gam = <r,r>, del = <r,w> (real*8 accumulators)
      !
      gam8 = 0.0d0
      del8 = 0.0d0
      !$acc parallel loop default(present) reduction(+:gam8,del8)
      !$omp parallel do schedule ( static ) private ( i, accv ) reduction(+:gam8,del8)
      do i = 1, nrows
         accv = nh_r(i) + nh_aw(i) * nh_r(nh_nbw(i)) &
                        + nh_ae(i) * nh_r(nh_nbe(i)) &
                        + nh_as(i) * nh_r(nh_nbs(i)) &
                        + nh_an(i) * nh_r(nh_nbn(i))
         nh_w(i) = accv
         gam8 = gam8 + nh_r(i) * nh_r(i)
         del8 = del8 + nh_r(i) * accv
      enddo
      !$omp end parallel do
      !
      if (sqrt(gam8) < nonh_tol * bnorm8) exit
      if (iter >= nonh_itermax) exit
      !
      ! scalar recurrences (Chronopoulos-Gear): beta, pAp, alpha
      !
      if (iter == 0) then
         beta8 = 0.0d0
         pap8  = del8
      else
         beta8 = gam8 / gamold8
         pap8  = del8 - beta8 * gam8 / alpha8     ! = <p_new, A p_new>
      endif
      if (pap8 <= 0.0d0) exit                     ! breakdown guard (SPD -> pAp > 0)
      alpha8  = gam8 / pap8
      gamold8 = gam8
      alpha   = real(alpha8, 4)
      beta    = real(beta8, 4)
      !
      ! kernel A: p = r + beta p ; s = w + beta s ; y += alpha p ; r -= alpha s
      !
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i )
      do i = 1, nrows
         nh_p(i) = nh_r(i) + beta * nh_p(i)
         nh_s(i) = nh_w(i) + beta * nh_s(i)
         nh_y(i) = nh_y(i) + alpha * nh_p(i)
         nh_r(i) = nh_r(i) - alpha * nh_s(i)
      enddo
      !$omp end parallel do
      !
      iter = iter + 1
      !
   enddo
   !
   ! diagnostics (reported in the end-of-run timing summary)
   !
   nh_iter_total  = nh_iter_total + iter
   nh_solve_count = nh_solve_count + 1
   if (iter > nh_iter_max) nh_iter_max = iter
   !
   ! 4c) unscale: pnh = S y (breaking rows have y = 0 -> pnh = 0 held exactly)
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i )
   do i = 1, nrows
      pnh(i) = nh_scl(i) * nh_y(i)
   enddo
   !$omp end parallel do
   !
   ! 5b) Optional global 2dx filter on pnh (nonh_filter > 0): one Jacobi
   !     neighbour-mean pass; snapshot into nh_w first (race-free). Damps the
   !     persistent grid-scale mode without touching the resolved smooth pressure.
   !
   if (nonh_filter > 0.0) then
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i )
      do i = 1, nrows
         nh_w(i) = pnh(i)
      enddo
      !$omp end parallel do
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i, sumn, cnt, j )
      do i = 1, nrows
         sumn = 0.0 ; cnt = 0
         j = nh_nbw(i) ; if (j > 0) then ; sumn = sumn + nh_w(j) ; cnt = cnt + 1 ; endif
         j = nh_nbe(i) ; if (j > 0) then ; sumn = sumn + nh_w(j) ; cnt = cnt + 1 ; endif
         j = nh_nbs(i) ; if (j > 0) then ; sumn = sumn + nh_w(j) ; cnt = cnt + 1 ; endif
         j = nh_nbn(i) ; if (j > 0) then ; sumn = sumn + nh_w(j) ; cnt = cnt + 1 ; endif
         if (cnt > 0) pnh(i) = (1.0 - nonh_filter) * nh_w(i) + nonh_filter * sumn / real(cnt)
      enddo
      !$omp end parallel do
   endif
   !
   ! 6) Correct fluxes / velocities with the pressure gradient -(dt/rho) G pnh.
   !    Each face writes only its own uv/q point -> race-free.
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( ip, ipuv, pL, pR, gf, unh )
   do ip = 1, nhuv
      if (nh_cR(ip) == 0.0 .and. nh_cL(ip) == 0.0) cycle
      ipuv = nh_faceuv(ip)
      pR = 0.0
      pL = 0.0
      if (nh_faceR(ip) > 0) pR = pnh(nh_faceR(ip))
      if (nh_faceL(ip) > 0) pL = pnh(nh_faceL(ip))
      gf  = nh_cR(ip) * pR + nh_cL(ip) * pL
      unh = - dtrho * gf
      uv(ipuv) = (1.0 - nonh_fnudge) * uv(ipuv) + nonh_fnudge * (uv(ipuv) + unh)
      q(ipuv)  = (1.0 - nonh_fnudge) * q(ipuv) + nonh_fnudge * (q(ipuv) + nh_hu(ip) * unh)
   enddo
   !$omp end parallel do
   !
   ! 7) Update the depth-mean vertical velocity ws for the next step's forcing
   !    (the wet checks use the current, possibly effective, zb -- they are
   !    depth checks).
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i, nm, f, accv )
   do i = 1, nrows
      nm = nm_index_of_row(i)
      !
      ! Non-hydrostatically dry cell: carry no nonh state (clean re-wetting).
      !
      if (zs(nm) - zb(nm) <= huthresh_nh) then
         ws(i)  = 0.0
         pnh(i) = 0.0
         nh_brfac(i) = 1.0
         cycle
      endif
      !
      ! Layer-consistent form: ws is the depth-mean vertical velocity
      ! (w_s + w_b)/2 of the Keller box; the bed kinematics are carried by
      ! the operator (layer term abf), so there is no separate w_b.
      !
      ! Boussinesq form: ws from the constraint G^T (hu u^{n+1}) + d(zb)/dt
      ! with the corrected velocities, PLUS the residual A p - b of this
      ! solve (scaled CG residual nh_r, unscaled with 1/scl; zero at
      ! breaking rows). The sum equals ws^n + A p - L p_f, i.e. the
      ! vertical equation with the solved pressure p and the horizontal
      ! correction with the applied (filtered) pressure p_f, so the next
      ! right-hand side is -(A p - b) plus the new forcing: the residual
      ! of every solve is cancelled by the next one, while the effect of
      ! the pnh filter is accepted as applied.
      ! Without the residual term the residual of every solve stays in
      ! the flow for good: at the 1 % tolerance one impulsive solve (LS1,
      ! 5 mm bed jump in a step) was a 2 cm spurious wave at the gauges.
      !
      accv = fmb * nh_dzbdt(i)
      f = nh_cellface(1, i)
      if (f > 0) accv = accv + nh_cR(f) * nh_hu(f) * uv(nh_faceuv(f))
      f = nh_cellface(2, i)
      if (f > 0) accv = accv + nh_cL(f) * nh_hu(f) * uv(nh_faceuv(f))
      f = nh_cellface(3, i)
      if (f > 0) accv = accv + nh_cR(f) * nh_hu(f) * uv(nh_faceuv(f))
      f = nh_cellface(4, i)
      if (f > 0) accv = accv + nh_cL(f) * nh_hu(f) * uv(nh_faceuv(f))
      ws(i) = accv - nh_r(i) / nh_scl(i)
      !
   enddo
   !$omp end parallel do
   !
   call timer_stop('non-hydrostatic')
   !
   end subroutine

end module
