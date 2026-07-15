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
! consistent depth. The bed SLOPES for the bottom kinematic condition w_b are
! FROZEN at initialization (nh_slbed): file bed on regular grids, subgrid
! cell-mean bed z_zmax - z_volmax/area on subgrid models.
!
module sfincs_nonhydrostatic
   !
   implicit none
   !
   ! Input parameters (sfincs.inp, read in sfincs_input)
   !
   real*4    :: nonh_fnudge      ! fraction of the pressure correction applied to uv/q
   real*4    :: nonh_tstop       ! time (s, absolute after input processing) to stop applying nonh corrections
   integer   :: nonh_itermax     ! max CG iterations per solve
   integer   :: nonh_fadein      ! open-boundary nonh fade-in width (cells); 0 = off
   real*4    :: nonh_filter      ! global 2dx pnh filter weight; 0 = off
   real*4    :: nonh_dzbmax      ! cap on |d(zb)/dx| in the bottom kinematic w_b; 0 = no cap
   real*4    :: nonh_brsteep     ! HFA steepness breaking onset; 0 = off
   real*4    :: nonh_brfr        ! Froude breaking onset (NEOWAVE); 0 = off -> use nonh_brsteep
   integer   :: nonh_brsmooth    ! breaking-flag smoothing passes: ramps pnh out over ~brsmooth+1 cells at the breaking-zone edges; 0 = sharp (one face)
   integer   :: nonh_slsmooth    ! frozen-bed-slope smoothing passes: bounds bed curvature d2zb/dx2 at slope breaks (e.g. island toe); 0 = off
   real*4    :: nonh_treform     ! breaking reformation time scale (s): released cells recover gradually, brfac += dt/treform; 0 = instant
   real*4    :: nonh_smoothbnd   ! localized 2dx pnh smoothing strength (fade-in / shallow zones)
   real*4    :: nonh_smoothdep   ! depth below which the localized smoothing also acts; 0 = off
   real*4    :: nonh_disp        ! Keller-box vertical factor (dispersion tuning; 1.0 ~ Airy, 2.0 strict)
   real*4    :: nonh_pmax        ! depth limiter: |pnh| <= nonh_pmax * rho g H; 0 = off
   real*4    :: nonh_tol         ! CG relative tolerance (scaled residual norm)
   logical   :: nonh_movingbed   ! add d(zb)/dt (= dzbext/dt) to the bottom kinematic w_b, so a moving seafloor radiates a depth-filtered (Kajiura-like) surface response
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
   real*4, dimension(:),   allocatable  :: ws         ! surface vertical velocity      (nrows)
   real*4, dimension(:),   allocatable  :: wb         ! bottom vertical velocity       (nrows)
   real*4, dimension(:),   allocatable  :: wb0        ! previous bottom v. velocity    (nrows)
   !
   real*4,  dimension(:), allocatable   :: Dnm        ! layer depth per row            (nrows)
   !
   ! Per-row inverse grid spacing in metres (handles quadtree refinement levels and
   ! geographic crsgeo lat-dependence). nh_dxr/nh_dyr are the per-CELL 1/dx, 1/dy used
   ! by the breaking dzdt divergence and the frozen bed slopes; nh_cf (below) is the
   ! per-FACE 0.5/dx used in the pressure operator.
   !
   real*4,  dimension(:), allocatable   :: nh_dxr     ! per-row 1/dx in metres         (nrows)
   real*4,  dimension(:), allocatable   :: nh_dyr     ! per-row 1/dy in metres         (nrows)
   !
   ! Bed slopes for the bottom kinematic condition w_b (step 7), FROZEN at
   ! initialization (file bed on regular grids, subgrid cell-mean bed
   ! z_zmax - z_volmax/area on subgrid models) and capped at +/- nonh_dzbmax.
   ! They must NOT follow the per-step effective bed zs - z_volume/area that
   ! subgrid models maintain (zb_effective): that bed moves with the wet fraction
   ! every step and its slope jitter would feed straight into w_b/ws. The real
   ! bed slope is constant in time.
   !
   real*4,  dimension(:,:), allocatable :: nh_slbed   ! (4, nrows) slots: 1 left(md), 2 right(mu), 3 below(nd), 4 above(nu)
   !
   real*4,  dimension(:), allocatable   :: nh_dvert   ! vertical-accel diagonal term   (nrows)
   real*4,  dimension(:), allocatable   :: nh_fade    ! open-boundary nonh fade-in 0..1 (nrows)
   real*4,  dimension(:), allocatable   :: nh_brfac   ! breaking STATE 1=full nonh .. 0=hydrostatic (nrows; hysteresis memory, unsmoothed)
   real*4,  dimension(:), allocatable   :: nh_bract   ! breaking ACTIVITY used by the pressure operator (nrows; = nh_brfac, optionally smoothed)
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
   ! Builds the row<->cell maps, the nonh face maps, the static neighbour map,
   ! the frozen bed slopes and the open-boundary fade-in, allocates the solver
   ! work arrays and pushes everything the per-step kernels touch onto the GPU.
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
   allocate(wb(nrows))
   allocate(wb0(nrows))
   allocate(Dnm(nrows))
   allocate(nm_index_of_row(nrows))
   allocate(nh_dxr(nrows))
   allocate(nh_dyr(nrows))
   allocate(nh_dvert(nrows))
   allocate(nh_fade(nrows))
   allocate(nh_brfac(nrows))
   allocate(nh_bract(nrows))
   !
   pnh      = 0.0
   ws       = 0.0
   wb       = 0.0
   wb0      = 0.0
   Dnm      = 0.0
   nh_dvert = 0.0
   nh_fade  = 1.0
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
         !
         ! Per-cell inverse grid spacing in metres. y-spacing is uniform per refinement
         ! level (dyrinv) in both projected and geographic mode. x-spacing is uniform per
         ! level (dxrinv) when projected, but lat-dependent per cell (1/dxm) when crsgeo.
         !
         nh_dyr(irow) = dyrinv(z_flags_iref(nm))
         if (crsgeo) then
            nh_dxr(irow) = 1.0 / dxm(nm)
         else
            nh_dxr(irow) = dxrinv(z_flags_iref(nm))
         endif
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
   ! Frozen, capped per-face bed slopes for the bottom kinematic condition
   ! w_b = u . d(zb)/dx (step 7). Computed here from the initialization-time zb
   ! (file bed on regular grids; subgrid cell-mean bed from sfincs_domain); on
   ! subgrid models zb is overwritten with the water-level-dependent effective
   ! bed every step (zb_effective), so the slopes must be taken now. The cap at +/- nonh_dzbmax
   ! clips near-vertical wall steps that would otherwise produce an enormous
   ! spurious w_b when a thin wall cell wets (2dt explicit instability); real
   ! bed slopes are well below it.
   !
   allocate(nh_slbed(4, nrows))
   nh_slbed = 0.0
   !
   block
      integer :: iuv, nmn
      real*4  :: slc
      !
      slc = nonh_dzbmax
      if (slc <= 0.0) slc = huge(1.0)   ! 0 = no cap
      !
      do irow = 1, nrows
         !
         nm = nm_index_of_row(irow)
         !
         iuv = z_index_uv_md(nm)
         if (iuv > 0) then
            nmn = uv_index_z_nm(iuv)
            nh_slbed(1, irow) = max(-slc, min(slc, (zb(nmn) - zb(nm)) * nh_dxr(irow)))
         endif
         !
         iuv = z_index_uv_mu(nm)
         if (iuv > 0) then
            nmn = uv_index_z_nmu(iuv)
            nh_slbed(2, irow) = max(-slc, min(slc, (zb(nm) - zb(nmn)) * nh_dxr(irow)))
         endif
         !
         if (nmax > 1) then
            !
            iuv = z_index_uv_nd(nm)
            if (iuv > 0) then
               nmn = uv_index_z_nm(iuv)
               nh_slbed(3, irow) = max(-slc, min(slc, (zb(nmn) - zb(nm)) * nh_dyr(irow)))
            endif
            !
            iuv = z_index_uv_nu(nm)
            if (iuv > 0) then
               nmn = uv_index_z_nmu(iuv)
               nh_slbed(4, irow) = max(-slc, min(slc, (zb(nm) - zb(nmn)) * nh_dyr(irow)))
            endif
            !
         endif
         !
      enddo
      !
   end block
   !
   ! Optional smoothing of the frozen bed slopes (nonh_slsmooth passes of a
   ! [1 2 1]/4 stencil over the nonh neighbours, per direction). The single-layer
   ! bottom kinematic condition w_b = u.d(zb)/dx makes the non-hydrostatic source
   ! sensitive to the bed CURVATURE d2(zb)/dx2: at a slope break (e.g. the toe of
   ! the conical island, where d(zb)/dx jumps over one cell) the per-cell w_b
   ! jumps by ~ u.dx.d2(zb)/dx2, an impulsive grid-scale forcing that radiates a
   ! trailing train of short dispersive waves. Smoothing ramps the slope over a
   ! few cells, bounding the curvature (~ slope/(N.dx) instead of slope/dx) and
   ! suppressing the wiggles at the source. Slopes are still frozen; this is the
   ! curvature analogue of the nonh_dzbmax slope cap. Default 0 = no smoothing.
   !
   if (nonh_slsmooth > 0) then
      block
         real*4, dimension(:,:), allocatable :: sltmp
         integer :: ps, ii, s, jw, je, js, jn
         real*4  :: accs
         allocate(sltmp(4, nrows))
         do ps = 1, nonh_slsmooth
            sltmp = nh_slbed
            do ii = 1, nrows
               jw = nh_nbw(ii) ; je = nh_nbe(ii)
               do s = 1, 2                       ! x-slopes smoothed along x-neighbours
                  accs = sltmp(s, ii)
                  if (jw > 0) accs = accs + 0.25 * (sltmp(s, jw) - sltmp(s, ii))
                  if (je > 0) accs = accs + 0.25 * (sltmp(s, je) - sltmp(s, ii))
                  nh_slbed(s, ii) = accs
               enddo
               js = nh_nbs(ii) ; jn = nh_nbn(ii)
               do s = 3, 4                       ! y-slopes smoothed along y-neighbours
                  accs = sltmp(s, ii)
                  if (js > 0) accs = accs + 0.25 * (sltmp(s, js) - sltmp(s, ii))
                  if (jn > 0) accs = accs + 0.25 * (sltmp(s, jn) - sltmp(s, ii))
                  nh_slbed(s, ii) = accs
               enddo
            enddo
         enddo
         deallocate(sltmp)
      end block
   endif
   !
   ! Open-boundary non-hydrostatic FADE-IN. Over nonh_fadein cells
   ! inside each open (water level / outflow) boundary, ramp the non-hydrostatic
   ! coupling from 0 (at the boundary) to full (nonh_fadein cells in). The incident
   ! wave then enters hydrostatically at the boundary and the nonh turns on
   ! gradually, avoiding the sharp pnh=0 -> full-nonh transition that reflects an
   ! incoming dispersive wave. nh_fade is 1 everywhere by default (no fade).
   !
   if (nonh_fadein > 0) then
      !
      block
         integer, dimension(:), allocatable :: idist
         integer :: pass, nmL, nmR, rL, rR
         !
         allocate(idist(nrows))
         idist = nonh_fadein + 1
         do ip = 1, nhuv
            inhuv = uv_index_of_nhuv(ip)
            nmL = uv_index_z_nm(inhuv)
            nmR = uv_index_z_nmu(inhuv)
            rL  = nh_faceL(ip)
            rR  = nh_faceR(ip)
            if (rL > 0 .and. rR == 0) then
               if (kcs(nmR) == 2 .or. kcs(nmR) == 3) idist(rL) = 0
            endif
            if (rR > 0 .and. rL == 0) then
               if (kcs(nmL) == 2 .or. kcs(nmL) == 3) idist(rR) = 0
            endif
         enddo
         do pass = 1, nonh_fadein
            do ip = 1, nhuv
               rL = nh_faceL(ip)
               rR = nh_faceR(ip)
               if (rL > 0 .and. rR > 0) then
                  if (idist(rL) + 1 < idist(rR)) idist(rR) = idist(rL) + 1
                  if (idist(rR) + 1 < idist(rL)) idist(rL) = idist(rR) + 1
               endif
            enddo
         enddo
         do irow = 1, nrows
            nh_fade(irow) = min(1.0, real(idist(irow)) / real(nonh_fadein))
         enddo
         deallocate(idist)
      end block
      !
   endif
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
   !
   !$omp parallel do schedule ( static ) private ( i )
   do i = 1, nrows
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
   !$acc                    nh_y, nh_r, nh_p, nh_s, nh_w, nh_b, &
   !$acc                    pnh, ws, wb, wb0, Dnm, nm_index_of_row, &
   !$acc                    nh_dvert, nh_fade, nh_brfac, nh_bract, nh_dxr, nh_dyr, nh_slbed, &
   !$acc                    nh_faceuv, nh_faceL, nh_faceR, nh_cellface, &
   !$acc                    nh_cf, nh_cR, nh_cL, nh_hu )
   !
   write(logstr,'(a,i0,a,i0,a)')'Non-hydrostatic solver (assembled-stencil CG): ', nrows, ' cells, ', nhuv, ' faces'
   call write_log(logstr, 0)
   !
   end subroutine


   subroutine compute_nonhydrostatic(dt, tloop)
   !
   ! Per-step non-hydrostatic pressure projection (see module header). Every
   ! numbered kernel below is one !$omp / !$acc loop.
   !
   use sfincs_data
   !
   implicit none
   !
   integer   :: count0, count1, count_rate, count_max
   real      :: tloop
   real*4    :: dt
   !
   integer   :: i, ip, ipuv, nm, nmn, nmu, iuv, iuvl, iuvr, j, jl, jr, f, iter, ipass, cnt, ineighbour
   real*4    :: dtrho, kbfac, hu, Dnm1, Dnmu, abf, fdep, tf, dval, sq, braw, accv, sumn, brrec, act
   real*4    :: qr, ql, dzdt, wmax, breform, pcap, gf, pL, pR, unh
   real*4    :: alpha, beta
   real*8    :: gam8, del8, bn8, gamold8, alpha8, beta8, pap8, bnorm8
   logical   :: bndok
   !
   call system_clock(count0, count_rate, count_max)
   !
   if (nrows == 0) return
   !
   ! Keller-box vertical factor (nonh_disp): 2.0 = strict linear-pressure single layer;
   ! lowering it (<2) strengthens the short-wave pressure feedback, flattening
   ! c(k) toward Airy (~1.0 optimal). Default 1.0.
   !
   kbfac = nonh_disp
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
   !     nonh_brfr > 0, else steepness/HFA (Smit, Zijlema & Stelling 2013) when
   !     nonh_brsteep > 0.
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
   elseif (nonh_brsteep > 0.0) then
      !
      ! Steepness/HFA criterion (gradual XBeach hydrostatic-front approximation).
      ! Detect steepening from the surface rate of rise built from the flux
      ! divergence: d(zs)/dx = dzdt / sqrt(g h), threshold nonh_brsteep. Only the
      ! rising front (dzdt > 0) is reduced; a breaking cell stays breaking until
      ! the surface falls.
      !
      breform = 0.25 * nonh_brsteep        ! XBeach reformsteep default (neighbour-spread threshold)
      !
      ! (a) raw surface rate-of-rise dzdt = -div(q) into nh_w
      !
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i, nm, iuvl, iuvr, ql, qr, dzdt )
      do i = 1, nrows
         nm   = nm_index_of_row(i)
         iuvr = z_index_uv_mu(nm) ; iuvl = z_index_uv_md(nm)
         qr = 0.0 ; ql = 0.0
         if (iuvr > 0) qr = q(iuvr)
         if (iuvl > 0) ql = q(iuvl)
         dzdt = - (qr - ql) * nh_dxr(i)
         if (nmax > 1) then
            iuvr = z_index_uv_nu(nm) ; iuvl = z_index_uv_nd(nm)
            qr = 0.0 ; ql = 0.0
            if (iuvr > 0) qr = q(iuvr)
            if (iuvl > 0) ql = q(iuvl)
            dzdt = dzdt - (qr - ql) * nh_dyr(i)
         endif
         nh_w(i) = dzdt
      enddo
      !$omp end parallel do
      !
      ! (b) two [1 2 1]/4 smoothing passes (2dx noise -> connected front)
      !
      do ipass = 1, 2
         !$acc parallel loop default(present)
         !$omp parallel do schedule ( static ) private ( i, accv, j )
         do i = 1, nrows
            accv = nh_w(i)
            j = nh_nbw(i) ; if (j > 0) accv = accv + 0.25 * (nh_w(j) - nh_w(i))
            j = nh_nbe(i) ; if (j > 0) accv = accv + 0.25 * (nh_w(j) - nh_w(i))
            j = nh_nbs(i) ; if (j > 0) accv = accv + 0.25 * (nh_w(j) - nh_w(i))
            j = nh_nbn(i) ; if (j > 0) accv = accv + 0.25 * (nh_w(j) - nh_w(i))
            nh_p(i) = accv
         enddo
         !$omp end parallel do
         !$acc parallel loop default(present)
         !$omp parallel do schedule ( static ) private ( i )
         do i = 1, nrows
            nh_w(i) = nh_p(i)
         enddo
         !$omp end parallel do
      enddo
      !
      ! (c) snapshot last step's breaking state (no within-sweep contamination)
      !
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i )
      do i = 1, nrows
         nh_r(i) = nh_brfac(i)
      enddo
      !$omp end parallel do
      !
      ! (d) hysteretic state machine on the smoothed dzdt (onset / neighbour-spread /
      !     release only when the surface falls)
      !
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i, dzdt, ineighbour, j, wmax )
      do i = 1, nrows
         dzdt = nh_w(i)
         ineighbour = 0
         j = nh_nbw(i) ; if (j > 0) then ; if (nh_r(j) < 1.0) ineighbour = 1 ; endif
         j = nh_nbe(i) ; if (j > 0) then ; if (nh_r(j) < 1.0) ineighbour = 1 ; endif
         j = nh_nbs(i) ; if (j > 0) then ; if (nh_r(j) < 1.0) ineighbour = 1 ; endif
         j = nh_nbn(i) ; if (j > 0) then ; if (nh_r(j) < 1.0) ineighbour = 1 ; endif
         wmax = sqrt(g * Dnm(i))
         if (nh_r(i) >= 1.0) then                                   ! was not breaking
            if (dzdt > nonh_brsteep * wmax) then
               nh_brfac(i) = 0.0
            elseif (dzdt > breform * wmax .and. ineighbour == 1) then
               nh_brfac(i) = 0.0
            else
               nh_brfac(i) = 1.0
            endif
         else                                                       ! was breaking
            if (dzdt < 0.0) then
               nh_brfac(i) = min(1.0, nh_r(i) + brrec)              ! release only when surface falls (gradually if nonh_treform > 0)
            else
               nh_brfac(i) = 0.0
            endif
         endif
      enddo
      !$omp end parallel do
      !
      ! (e) close interior holes (a non-breaking cell flanked on both x-sides, or
      !     both y-sides, by breaking cells is filled in; 3 passes)
      !
      do ipass = 1, 3
         !$acc parallel loop default(present)
         !$omp parallel do schedule ( static ) private ( i )
         do i = 1, nrows
            nh_r(i) = nh_brfac(i)
         enddo
         !$omp end parallel do
         !$acc parallel loop default(present)
         !$omp parallel do schedule ( static ) private ( i, ineighbour, jl, jr )
         do i = 1, nrows
            if (nh_r(i) >= 1.0) then
               ineighbour = 0
               jl = nh_nbw(i) ; if (jl > 0) then ; if (nh_r(jl) < 1.0) ineighbour = ineighbour + 1 ; endif
               jr = nh_nbe(i) ; if (jr > 0) then ; if (nh_r(jr) < 1.0) ineighbour = ineighbour + 2 ; endif
               if (ineighbour == 3) then
                  nh_brfac(i) = 0.0
               endif
               if (nmax > 1 .and. nh_brfac(i) >= 1.0) then
                  ineighbour = 0
                  jl = nh_nbs(i) ; if (jl > 0) then ; if (nh_r(jl) < 1.0) ineighbour = ineighbour + 1 ; endif
                  jr = nh_nbn(i) ; if (jr > 0) then ; if (nh_r(jr) < 1.0) ineighbour = ineighbour + 2 ; endif
                  if (ineighbour == 3) nh_brfac(i) = 0.0
               endif
            endif
         enddo
         !$omp end parallel do
      enddo
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
   ! 1c) Breaking ACTIVITY for the pressure operator: nh_bract = nh_brfac,
   !     optionally smoothed over the nonh neighbours (nonh_brsmooth passes of a
   !     [1 2 1]/4 stencil). The binary flag switches the pressure off across a
   !     single face, which prints a kink on the surface at the rear edge of the
   !     breaking region; the smoothed activity ramps the pressure out over a few
   !     cells instead. The hysteresis STATE stays in nh_brfac (unsmoothed).
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i )
   do i = 1, nrows
      nh_bract(i) = nh_brfac(i)
   enddo
   !$omp end parallel do
   !
   do ipass = 1, nonh_brsmooth
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i )
      do i = 1, nrows
         nh_w(i) = nh_bract(i)
      enddo
      !$omp end parallel do
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i, accv, j )
      do i = 1, nrows
         accv = nh_w(i)
         j = nh_nbw(i) ; if (j > 0) accv = accv + 0.25 * (nh_w(j) - nh_w(i))
         j = nh_nbe(i) ; if (j > 0) accv = accv + 0.25 * (nh_w(j) - nh_w(i))
         j = nh_nbs(i) ; if (j > 0) accv = accv + 0.25 * (nh_w(j) - nh_w(i))
         j = nh_nbn(i) ; if (j > 0) accv = accv + 0.25 * (nh_w(j) - nh_w(i))
         nh_bract(i) = accv
      enddo
      !$omp end parallel do
   enddo
   !
   ! 2) Per-face gradient coefficients (wet/dry, layer term abf, conveyance depth,
   !    open-boundary fade-in).
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( ip, ipuv, nm, nmu, hu, Dnm1, Dnmu, abf, fdep, bndok )
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
               if (nonh_fadein > 0) then
                  fdep = 1.0
                  if (nh_faceL(ip) > 0) fdep = min(fdep, nh_fade(nh_faceL(ip)))
                  if (nh_faceR(ip) > 0) fdep = min(fdep, nh_fade(nh_faceR(ip)))
                  nh_cR(ip) = nh_cR(ip) * fdep
                  nh_cL(ip) = nh_cL(ip) * fdep
               endif
            endif
         endif
      endif
   enddo
   !$omp end parallel do
   !
   ! 3) ASSEMBLY (once per solve). Per row: raw off-diagonals a_ij = dtrho*hu*cR*cL
   !    (one product per incident face, same both ways -> symmetric), the diagonal
   !    dval = dvert + sum dtrho*hu*coef^2, the scaling scl = 1/sqrt(dval), the
   !    breaking-masked scaling sclo, the SCALED right-hand side
   !    b = scl * [ G^T (hu u*) - (ws + wb0 - 2 wb) ], and the warm start
   !    y0 = pnh * sqrt(dval) (scaled previous pressure; 0 at breaking cells).
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i, dval, braw, f, tf, sq, act )
   do i = 1, nrows
      dval = nh_dvert(i)
      braw = - (ws(i) + wb0(i) - 2.0 * wb(i))
      f = nh_cellface(1, i)
      if (f > 0) then
         tf = dtrho * nh_hu(f)
         nh_aw(i) = tf * nh_cR(f) * nh_cL(f)
         dval = dval + tf * nh_cR(f)**2
         braw = braw + nh_cR(f) * nh_hu(f) * uv(nh_faceuv(f))
      else
         nh_aw(i) = 0.0
      endif
      f = nh_cellface(2, i)
      if (f > 0) then
         tf = dtrho * nh_hu(f)
         nh_ae(i) = tf * nh_cL(f) * nh_cR(f)
         dval = dval + tf * nh_cL(f)**2
         braw = braw + nh_cL(f) * nh_hu(f) * uv(nh_faceuv(f))
      else
         nh_ae(i) = 0.0
      endif
      f = nh_cellface(3, i)
      if (f > 0) then
         tf = dtrho * nh_hu(f)
         nh_as(i) = tf * nh_cR(f) * nh_cL(f)
         dval = dval + tf * nh_cR(f)**2
         braw = braw + nh_cR(f) * nh_hu(f) * uv(nh_faceuv(f))
      else
         nh_as(i) = 0.0
      endif
      f = nh_cellface(4, i)
      if (f > 0) then
         tf = dtrho * nh_hu(f)
         nh_an(i) = tf * nh_cL(f) * nh_cR(f)
         dval = dval + tf * nh_cL(f)**2
         braw = braw + nh_cL(f) * nh_hu(f) * uv(nh_faceuv(f))
      else
         nh_an(i) = 0.0
      endif
      sq = sqrt(dval)              ! dval >= dvert > 0 always (Dnm capped)
      nh_scl(i) = 1.0 / sq
      !
      ! Breaking activity act in [0,1]: act = 1 full nonh; act = 0 holds pnh = 0
      ! (row+column zeroed via sclo, unit scaled diagonal). Fractional act (from
      ! nonh_brsmooth / nonh_treform) ramps the pressure out smoothly; the
      ! operator stays SPD (T A T + diag(1 - act^2), act <= 1).
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
   ! 5c) Localized 2dx smoothing in the marginal-nonh zones (fade-in zone and/or
   !     shallow run-up water, D < nonh_smoothdep); zero in the resolved interior.
   !
   if (nonh_smoothbnd > 0.0 .and. (nonh_fadein > 0 .or. nonh_smoothdep > 0.0)) then
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i )
      do i = 1, nrows
         nh_w(i) = pnh(i)
      enddo
      !$omp end parallel do
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i, sumn, cnt, j, fdep, gf )
      do i = 1, nrows
         sumn = 0.0 ; cnt = 0
         j = nh_nbw(i) ; if (j > 0) then ; sumn = sumn + nh_w(j) ; cnt = cnt + 1 ; endif
         j = nh_nbe(i) ; if (j > 0) then ; sumn = sumn + nh_w(j) ; cnt = cnt + 1 ; endif
         j = nh_nbs(i) ; if (j > 0) then ; sumn = sumn + nh_w(j) ; cnt = cnt + 1 ; endif
         j = nh_nbn(i) ; if (j > 0) then ; sumn = sumn + nh_w(j) ; cnt = cnt + 1 ; endif
         if (cnt > 0) then
            fdep = 0.0
            if (nonh_fadein > 0)      fdep = 1.0 - nh_fade(i)
            if (nonh_smoothdep > 0.0) fdep = max(fdep, (nonh_smoothdep - Dnm(i)) / nonh_smoothdep)
            gf = nonh_smoothbnd * max(0.0, min(1.0, fdep))
            pnh(i) = (1.0 - gf) * nh_w(i) + gf * sumn / real(cnt)
         endif
      enddo
      !$omp end parallel do
   endif
   !
   ! 5e) Depth limiter: |pnh| <= nonh_pmax * rho g H. The total bed pressure
   !     rho*g*H + pnh then stays >= 0 for nonh_pmax <= 1 (no "suction"); resolved
   !     waves have |pnh| << rho*g*H so are untouched -- this only clips spikes
   !     at steep wall run-up, breaking fronts and grid-scale modes.
   !
   if (nonh_pmax > 0.0) then
      !$acc parallel loop default(present)
      !$omp parallel do schedule ( static ) private ( i, pcap )
      do i = 1, nrows
         pcap = nonh_pmax * rhow * g * Dnm(i)
         pnh(i) = max(-pcap, min(pcap, pnh(i)))
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
   ! 7) Update surface/bottom vertical velocities for the next step's forcing
   !    (bottom kinematic condition w_b = u . d(zb)/dx with the slopes FROZEN at
   !    initialization in nh_slbed; the wet checks use the current, possibly
   !    effective, zb -- they are depth checks).
   !
   !$acc parallel loop default(present)
   !$omp parallel do schedule ( static ) private ( i, nm, iuv, nmn )
   do i = 1, nrows
      nm = nm_index_of_row(i)
      !
      ! Non-hydrostatically dry cell: carry no nonh state (clean re-wetting).
      !
      if (zs(nm) - zb(nm) <= huthresh_nh) then
         wb0(i) = 0.0
         wb(i)  = 0.0
         ws(i)  = 0.0
         pnh(i) = 0.0
         nh_brfac(i) = 1.0
         cycle
      endif
      !
      wb0(i) = wb(i)
      wb(i)  = 0.0
      !
      ! w_b over faces to non-hydrostatically wet neighbours only.
      !
      iuv = z_index_uv_md(nm)
      if (kfuv(iuv) > 0) then
         nmn = uv_index_z_nm(iuv)
         if (zs(nmn) - zb(nmn) > huthresh_nh) then
            wb(i) = wb(i) - 0.5 * uv(iuv) * nh_slbed(1, i)
         endif
      endif
      iuv = z_index_uv_mu(nm)
      if (kfuv(iuv) > 0) then
         nmn = uv_index_z_nmu(iuv)
         if (zs(nmn) - zb(nmn) > huthresh_nh) then
            wb(i) = wb(i) - 0.5 * uv(iuv) * nh_slbed(2, i)
         endif
      endif
      if (nmax > 1) then
         iuv = z_index_uv_nd(nm)
         if (kfuv(iuv) > 0) then
            nmn = uv_index_z_nm(iuv)
            if (zs(nmn) - zb(nmn) > huthresh_nh) then
               wb(i) = wb(i) - 0.5 * uv(iuv) * nh_slbed(3, i)
            endif
         endif
         iuv = z_index_uv_nu(nm)
         if (kfuv(iuv) > 0) then
            nmn = uv_index_z_nmu(iuv)
            if (zs(nmn) - zb(nmn) > huthresh_nh) then
               wb(i) = wb(i) - 0.5 * uv(iuv) * nh_slbed(4, i)
            endif
         endif
      endif
      !
      ws(i) = ws(i) - (wb(i) - wb0(i)) + (kbfac * dt / (rhow * Dnm(i))) * pnh(i)
   enddo
   !$omp end parallel do
   !
   ! Optional moving-bed source: add d(zb)/dt = dzbext/dt to the bottom kinematic
   ! w_b (full condition w_b = d(zb)/dt + u.d(zb)/dx), so a moving seafloor
   ! radiates a depth-filtered (Kajiura-like) surface response instead of being
   ! stamped onto zs.  Kept in a SEPARATE flag-guarded loop so the step-7 loop
   ! above is byte-for-byte the original when the feature is off.  Adding the term
   ! here (wb += d(zb)/dt, ws -= d(zb)/dt) is algebraically identical to carrying
   ! it inside step 7 (it cancels through the ws = ws - (wb - wb0) update).
   if (nonh_movingbed .and. use_dzbext .and. dt > 0.0) then
      !$omp parallel do schedule ( static ) private ( i, nm )
      do i = 1, nrows
         nm = nm_index_of_row(i)
         if (zs(nm) - zb(nm) <= huthresh_nh) cycle
         wb(i) = wb(i) + dzbext(nm) / dt
         ws(i) = ws(i) - dzbext(nm) / dt
      enddo
      !$omp end parallel do
   endif
   !
   call system_clock(count1, count_rate, count_max)
   tloop = tloop + 1.0*(count1 - count0)/count_rate
   !
   end subroutine

end module
