module sfincs_semi_implicit
   !
   ! Semi-implicit pressure treatment for SFINCS (Casulli-style theta-method)
   !
   ! Supports regular and quadtree grids. Subgrid is not supported yet: it makes the
   ! system mildly nonlinear, because the storage term V(eta) lands on the diagonal and
   ! its derivative is the wet area, which depends on the solution.
   !
   ! Solves a Helmholtz equation for free surface elevation:
   !   eta^{n+1} - theta * g * dt^2 * div(H * grad(eta^{n+1})) / A = RHS
   !
   ! NOTE the single power of theta. The assembled coefficient carries theta once, from
   ! si_coeff in sfincs_momentum.f90; the RHS applies no (1-theta) weighting, so continuity
   ! is fully implicit while momentum is theta-weighted. Earlier comments here claimed
   ! theta^2, which the code has never done.
   !
   ! Uses Conjugate Gradient solver with SSOR preconditioning.
   ! The system is SPD so CG is the natural choice.
   !
   use sfincs_data
   use sfincs_groundwater
   !
   implicit none
   !
   private
   public :: initialize_semi_implicit, assemble_and_solve_pressure, backsubstitute_fluxes_si
   public :: get_tloop_si, get_si_iter_avg, get_si_iter_max
   public :: get_si_outer_avg, get_si_outer_max, get_si_outer_capped
   public :: get_si_outer_profile
   !
   ! CG solver work arrays (allocated once in initialize, reused each timestep)
   !
   real*4, dimension(:), allocatable :: cg_r      ! residual
   real*4, dimension(:), allocatable :: cg_z      ! preconditioned residual
   real*4, dimension(:), allocatable :: cg_p      ! search direction
   real*4, dimension(:), allocatable :: cg_Ap     ! matrix-vector product
   real*4, dimension(:), allocatable :: cg_diag   ! diagonal (for Jacobi preconditioner)
   !
   ! Per-row face list. The assembly walks the faces belonging to each row instead of
   ! four hard-coded directions, so it does not care how many neighbours a row has.
   ! On a regular grid every row has up to 4. On a quadtree a row next to a refinement
   ! transition has up to 8, because a coarse-to-fine face is stored as TWO UV points
   ! (z_index_uv_mu1/mu2 and friends, sfincs_data.f90:302-309).
   ! Keeping the outer loop over rows rather than over faces matters: a bare face loop
   ! would race on the diagonal under OpenMP, and fixing that with atomics would make
   ! the summation order non-deterministic between runs.
   !
   integer, dimension(:), allocatable :: si_row_face_ptr    ! nrows_si+1, CSR-style
   integer, dimension(:), allocatable :: si_row_face_ip     ! UV point index
   integer, dimension(:), allocatable :: si_row_face_slot   ! slot in si_AA, 0 = Dirichlet
   integer, dimension(:), allocatable :: si_row_face_isy    ! 0 = x-direction, 1 = y
   integer, dimension(:), allocatable :: si_row_face_bnd    ! cell to read zs from when Dirichlet
   integer, dimension(:), allocatable :: si_row_face_gwslot ! slot in si_AA, aquifer block
   real*4,  dimension(:), allocatable :: si_row_face_dinv   ! 1 / centre-to-centre distance
   !
   ! Position of the diagonal within each CSR row. This is
   ! what lets the preconditioner split lower from upper without a fixed stencil.
   !
   integer, dimension(:), allocatable :: si_diag_ptr     ! nrows_tot
   !
   ! Groundwater block bookkeeping: where each aquifer face and each exchange entry
   ! lands in si_AA, and the total row count including the aquifer half.
   !
   integer, dimension(:,:), allocatable :: gw_face_slot  ! (8, nrows_si)
   integer, dimension(:), allocatable :: si_exch_ptr     ! nrows_tot
   integer :: nrows_tot
   !
   ! What the LAST assembly pass actually put into the matrix, kept so the water budget can
   ! measure the operator that was applied rather than re-deriving it.
   !
   ! Re-deriving is not good enough. Transmissivity and the exchange conductance are both lagged
   ! one Picard iterate behind the solution, so evaluating them from the converged head gives a
   ! slightly different operator -- on the Edelman step response that alone put the reported
   ! closure at 0.03% when the scheme itself was conserving to round-off. A budget that measures
   ! a different operator from the one the solver used cannot tell a leak from its own lag.
   !
   real*4, dimension(:), allocatable :: gw_coeff_applied  ! nfaces_si, transmissivity * width / d * dt
   real*4, dimension(:), allocatable :: gw_cexch_applied  ! nrows_si, exchange conductance * dt
   real*4, dimension(:), allocatable :: gw_qexpl_applied  ! nrows_si, lagged exchange remainder
   real*4, dimension(:), allocatable :: gw_cseep_applied  ! nrows_si, seepage conductance * dt
   !
   ! Lower bound on the wet-area derivative, as a fraction of cell area. Sets the worst
   ! spread the coupled diagonal can take, and so the conditioning CG has to cope with.
   !
   real*4, parameter :: awet_floor = 0.01
   !
   ! Current outer-iterate water level, used only when subgrid is on. The subgrid
   ! storage relation V(eta) sits on the diagonal and its derivative is the wetted
   ! area, which depends on the level being solved for, so the system is mildly
   ! nonlinear and needs an outer iteration.
   !
   real*8, dimension(:), allocatable :: si_eta_k       ! nrows_si
   !
   ! Symmetric (Jacobi) scaling factors, 1/sqrt(diagonal). The subgrid diagonal spans
   ! several orders of magnitude between isolated dry cells and well-connected wet
   ! ones, which stalls CG. Scaling preserves symmetry and definiteness.
   !
   real*4, dimension(:), allocatable :: si_scale       ! nrows_si
   !
   ! Per-timestep parts of the surface rows that do NOT depend on the outer iterate. Built once
   ! before the outer loop in assemble_and_solve_pressure; the loop adds only the storage term.
   !
   real*4, dimension(:), allocatable :: si_cface       ! nrows_si, sum of face conductances on the row
   real*4, dimension(:), allocatable :: si_cbnd        ! nrows_si, the Dirichlet share of si_cface
   real*8, dimension(:), allocatable :: si_rhs_const   ! nrows_si, flux divergence + sources + known boundary levels
   real*8, dimension(:), allocatable :: si_vol_n       ! nrows_si, subgrid volume at the old-time level
   !
   ! Diagonally scaled copy of the matrix handed to CG. si_AA itself stays unscaled so that the
   ! hoisted off-diagonals survive from one outer iteration to the next.
   !
   real*4, dimension(:), allocatable :: si_AA_scaled   ! nnz_si
   !
   ! Timing and diagnostics
   !
   real    :: tloop_si
   integer :: si_iter_total       ! total CG iterations across all timesteps
   integer :: si_solve_count      ! number of solver calls
   integer :: si_iter_max_seen    ! max iterations in any single solve
   integer :: si_outer_total      ! total nonlinear outer iterations (subgrid)
   integer :: si_outer_max_seen   ! worst outer count in any timestep
   integer :: si_outer_capped     ! timesteps that hit si_maxouter
   integer :: si_outer_stagnant   ! outer loops stopped on stagnation
   !
   ! Per outer-iteration index: how many timesteps reached it, and the summed nbad and
   ! dmax_outer there. Says how fast the nonlinear iteration actually converges.
   integer,   dimension(:), allocatable :: si_prof_count
   integer*8, dimension(:), allocatable :: si_prof_nbad
   real*8,    dimension(:), allocatable :: si_prof_dmax
   !
   ! Outer-loop convergence: the field counts as converged when at most this fraction of rows
   ! still changes by more than si_tolouter. A handful of rows sitting on a wet/dry threshold
   ! flip between two states every iteration and would keep a max-norm test from ever passing.
   real*4, parameter :: si_outer_frac = 1.0e-3
   integer :: si_solve_count_outer ! timesteps with an outer loop
   !
contains
   !
   subroutine initialize_semi_implicit()
   !
   ! Build the sparse matrix structure (CSR format) for the pressure Helmholtz system.
   ! The sparsity pattern is static; only the values change each timestep.
   ! Rows are variable length: 5 entries on a regular grid, up to 9 next to a quadtree
   ! refinement transition.
   !
   implicit none
   !
   integer :: nm, ip, irow, icol, k
   integer :: j, nb, nfaces_si, maxrow, ndiag, kn, nperm
   real*4  :: dref
   !
   integer, dimension(:), allocatable :: row_count
   integer, dimension(:), allocatable :: perm       ! old CSR slot -> new CSR slot
   integer, dimension(:), allocatable :: col_tmp
   integer, dimension(:,:), allocatable :: face_ip    ! (8, nrows_si) UV point per face slot
   integer, dimension(:,:), allocatable :: face_nb    ! (8, nrows_si) neighbour row, 0 = Dirichlet
   integer, dimension(:,:), allocatable :: face_bnd   ! (8, nrows_si) neighbour cell index
   integer, dimension(:,:), allocatable :: face_isy   ! (8, nrows_si) 0 = x, 1 = y
   integer, dimension(:,:), allocatable :: face_slot  ! (8, nrows_si) slot in si_AA
   !
   !
   ! Allocate row mapping arrays
   !
   allocate(si_row_of_nm(np))
   si_row_of_nm = 0
   !
   ! Count number of active interior cells (kcs==1) for the linear system
   ! Boundary cells (kcs==2) are Dirichlet and not unknowns
   !
   nrows_si = 0
   !
   do nm = 1, np
      if (kcs(nm) == 1) then
         nrows_si = nrows_si + 1
      endif
   enddo
   !
   ! Total unknowns. With groundwater the aquifer head is a second block of the same size, so
   ! every solver work array below is sized for both.
   !
   if (gwflow) then
      nrows_tot = 2 * nrows_si
   else
      nrows_tot = nrows_si
   endif
   !
   ! Allocate solver arrays
   !
   allocate(si_nm_of_row(nrows_si))
   allocate(si_row_ptr(nrows_tot + 1))
   allocate(face_ip(8, nrows_si))
   allocate(face_nb(8, nrows_si))
   allocate(face_bnd(8, nrows_si))
   allocate(face_isy(8, nrows_si))
   allocate(face_slot(8, nrows_si))
   allocate(si_rhs(nrows_tot))
   allocate(si_x(nrows_tot))
   allocate(si_dx(nrows_tot))
   allocate(si_b(nrows_tot))
   ! Sized exactly like q and uv (sfincs_domain.f90:2196), NOT npuv.
   !
   ! A quadtree creates ncuv combined uv points that live past npuv, and div_qstar below
   ! reads z_index_uv_md/mu/nd/nu, which point at those combined points next to a refinement
   ! transition. The +1 is the sentinel slot sfincs_domain.f90:1251 assigns to unset indices.
   ! Allocating only npuv reads past the end -- silently in Release, and it did.
   !
   allocate(si_q_star(npuv + ncuv + 1))
   allocate(si_coeff(npuv + ncuv + 1))
   !
   ! CG work arrays (allocated once, reused every timestep)
   !
   allocate(cg_r(nrows_tot))
   allocate(cg_z(nrows_tot))
   allocate(cg_p(nrows_tot))
   allocate(cg_Ap(nrows_tot))
   allocate(cg_diag(nrows_tot))
   allocate(si_eta_k(nrows_tot))
   allocate(si_scale(nrows_tot))
   allocate(si_cface(nrows_si))
   allocate(si_cbnd(nrows_si))
   allocate(si_rhs_const(nrows_si))
   allocate(si_vol_n(nrows_si))
   !
   si_nm_of_row = 0
   si_row_ptr = 0
   face_ip = 0
   face_nb = 0
   face_bnd = 0
   face_isy = 0
   face_slot = 0
   si_rhs = 0.0
   si_x = 0.0
   si_q_star = 0.0
   si_coeff = 0.0
   cg_r = 0.0
   cg_z = 0.0
   cg_p = 0.0
   cg_Ap = 0.0
   cg_diag = 0.0
   si_eta_k = 0.0
   si_scale = 1.0
   !
   tloop_si = 0.0
   si_iter_total = 0
   si_solve_count = 0
   si_iter_max_seen = 0
   si_outer_total = 0
   si_outer_max_seen = 0
   si_outer_capped = 0
   si_outer_stagnant = 0
   si_solve_count_outer = 0
   !
   allocate(si_prof_count(si_maxouter))
   allocate(si_prof_nbad(si_maxouter))
   allocate(si_prof_dmax(si_maxouter))
   si_prof_count = 0
   si_prof_nbad  = 0
   si_prof_dmax  = 0.0d0
   !
   ! Build row <-> nm mapping
   !
   irow = 0
   do nm = 1, np
      if (kcs(nm) == 1) then
         irow = irow + 1
         si_row_of_nm(nm) = irow
         si_nm_of_row(irow) = nm
      endif
   enddo
   !
   ! Enumerate each row's faces from the per-cell UV index arrays.
   !
   ! A cell has up to EIGHT faces, not four: sfincs_data.f90:302-309 holds md1/md2, mu1/mu2,
   ! nd1/nd2 and nu1/nu2, and on a quadtree a coarse cell facing refinement uses both slots
   ! on that side (sfincs_domain.f90:806-860 emits two UV points there, each pairing the
   ! coarse cell with one fine cell). On a regular grid only the *1 slots are set, so the
   ! same code produces the old four-face stencil.
   !
   ! The slot order below is deliberate: 1-4 come before the diagonal and 5-8 after it, so a
   ! regular grid reproduces the historical CSR order left, bottom, centre, top, right.
   !
   do irow = 1, nrows_si
      !
      nm = si_nm_of_row(irow)
      !
      face_ip(1, irow) = z_index_uv_md1(nm)   ! left
      face_ip(2, irow) = z_index_uv_md2(nm)   ! left, second fine neighbour
      face_ip(3, irow) = z_index_uv_nd1(nm)   ! bottom
      face_ip(4, irow) = z_index_uv_nd2(nm)   ! bottom, second
      face_ip(5, irow) = z_index_uv_nu1(nm)   ! top
      face_ip(6, irow) = z_index_uv_nu2(nm)   ! top, second
      face_ip(7, irow) = z_index_uv_mu1(nm)   ! right
      face_ip(8, irow) = z_index_uv_mu2(nm)   ! right, second
      !
      do j = 1, 8
         !
         ip = face_ip(j, irow)
         if (ip <= 0) cycle
         !
         ! The neighbour is whichever end of the face is not this cell.
         !
         if (uv_index_z_nm(ip) == nm) then
            nb = uv_index_z_nmu(ip)
         else
            nb = uv_index_z_nm(ip)
         endif
         !
         face_nb(j, irow)  = si_row_of_nm(nb)   ! 0 when the neighbour is a Dirichlet cell
         face_bnd(j, irow) = nb                 ! cell to read zs from in that case
         !
         if (j <= 2 .or. j >= 7) then
            face_isy(j, irow) = 0               ! x-direction
         else
            face_isy(j, irow) = 1               ! y-direction
         endif
         !
      enddo
      !
   enddo
   !
   ! Build the CSR sparsity pattern. Rows are variable length: 5 entries on a regular grid,
   ! up to 9 next to a refinement transition.
   !
   ! With groundwater the system carries two unknowns per cell: surface head in rows
   ! 1..nrows_si, aquifer head in rows nrows_si+1..2*nrows_si. A groundwater row has the same
   ! lateral connectivity as its surface counterpart, and both gain one exchange entry linking
   ! them. Nothing below assumes a stencil size, so this is a row-count change and not a
   ! restructuring.
   !
   allocate(row_count(nrows_tot))
   !
   do irow = 1, nrows_si
      row_count(irow) = 1                       ! the diagonal is always present
      do j = 1, 8
         if (face_nb(j, irow) > 0) row_count(irow) = row_count(irow) + 1
      enddo
      if (gwflow) then
         row_count(irow) = row_count(irow) + 1              ! exchange with the aquifer
         row_count(nrows_si + irow) = row_count(irow)       ! same lateral pattern + exchange
      endif
   enddo
   !
   si_row_ptr(1) = 1
   do irow = 1, nrows_tot
      si_row_ptr(irow + 1) = si_row_ptr(irow) + row_count(irow)
   enddo
   !
   nnz_si = si_row_ptr(nrows_tot + 1) - 1
   !
   deallocate(row_count)
   !
   allocate(si_col_idx(nnz_si))
   allocate(si_diag_ptr(nrows_tot))
   if (gwflow) then
      allocate(gw_face_slot(8, nrows_si))
      allocate(si_exch_ptr(nrows_tot))
      gw_face_slot = 0
      si_exch_ptr = 0
   endif
   !
   si_col_idx = 0
   si_diag_ptr = 0
   !
   ! Columns are emitted in face order, not sorted ascending. After the row builder, every row
   ! is reordered once into [lower | diagonal | upper] (see below); the SSOR sweeps rely on it.
   !
   do irow = 1, nrows_si
      !
      k = si_row_ptr(irow) - 1
      !
      do j = 1, 4                               ! faces emitted before the diagonal
         if (face_nb(j, irow) > 0) then
            k = k + 1
            si_col_idx(k) = face_nb(j, irow)
            face_slot(j, irow) = k
         endif
      enddo
      !
      k = k + 1
      si_col_idx(k) = irow
      si_diag_ptr(irow) = k
      !
      do j = 5, 8                               ! faces emitted after the diagonal
         if (face_nb(j, irow) > 0) then
            k = k + 1
            si_col_idx(k) = face_nb(j, irow)
            face_slot(j, irow) = k
         endif
      enddo
      !
      if (gwflow) then
         k = k + 1
         si_col_idx(k) = nrows_si + irow         ! exchange, upper block
         si_exch_ptr(irow) = k
      endif
      !
   enddo
   !
   if (gwflow) then
      !
      do irow = 1, nrows_si
         !
         k = si_row_ptr(nrows_si + irow) - 1
         !
         k = k + 1
         si_col_idx(k) = irow                    ! exchange, lower block
         si_exch_ptr(nrows_si + irow) = k
         !
         do j = 1, 4
            if (face_nb(j, irow) > 0) then
               k = k + 1
               si_col_idx(k) = nrows_si + face_nb(j, irow)
               gw_face_slot(j, irow) = k
            endif
         enddo
         !
         k = k + 1
         si_col_idx(k) = nrows_si + irow
         si_diag_ptr(nrows_si + irow) = k
         !
         do j = 5, 8
            if (face_nb(j, irow) > 0) then
               k = k + 1
               si_col_idx(k) = nrows_si + face_nb(j, irow)
               gw_face_slot(j, irow) = k
            endif
         enddo
         !
      enddo
      !
   endif
   !
   allocate(si_AA(nnz_si))
   si_AA = 0.0
   allocate(si_AA_scaled(nnz_si))
   si_AA_scaled = 0.0
   !
   ! Order every row as [lower | diagonal | upper] so the SSOR sweeps run over bounded index
   ! ranges instead of testing each column against the row. The row builder emits faces 1-4
   ! before the diagonal and 5-8 after it; on a regular grid that already is the split, on a
   ! quadtree a neighbour across a refinement transition can land on the wrong side. A stable
   ! partition keeps the relative order within each side, so wherever nothing moves the sums
   ! in the sweeps round exactly as before. Every table that points into si_AA is remapped.
   !
   allocate(perm(nnz_si))
   allocate(col_tmp(nnz_si))
   col_tmp = si_col_idx
   nperm = 0
   !
   do irow = 1, nrows_tot
      kn = si_row_ptr(irow) - 1
      do k = si_row_ptr(irow), si_row_ptr(irow + 1) - 1
         if (col_tmp(k) < irow) then
            kn = kn + 1
            perm(k) = kn
         endif
      enddo
      do k = si_row_ptr(irow), si_row_ptr(irow + 1) - 1
         if (col_tmp(k) == irow) then
            kn = kn + 1
            perm(k) = kn
         endif
      enddo
      do k = si_row_ptr(irow), si_row_ptr(irow + 1) - 1
         if (col_tmp(k) > irow) then
            kn = kn + 1
            perm(k) = kn
         endif
      enddo
      do k = si_row_ptr(irow), si_row_ptr(irow + 1) - 1
         if (perm(k) /= k) then
            nperm = nperm + 1
            exit
         endif
      enddo
   enddo
   !
   do k = 1, nnz_si
      si_col_idx(perm(k)) = col_tmp(k)
   enddo
   do irow = 1, nrows_tot
      si_diag_ptr(irow) = perm(si_diag_ptr(irow))
   enddo
   do irow = 1, nrows_si
      do j = 1, 8
         if (face_slot(j, irow) > 0) face_slot(j, irow) = perm(face_slot(j, irow))
      enddo
   enddo
   if (gwflow) then
      do irow = 1, nrows_tot
         si_exch_ptr(irow) = perm(si_exch_ptr(irow))
      enddo
      do irow = 1, nrows_si
         do j = 1, 8
            if (gw_face_slot(j, irow) > 0) gw_face_slot(j, irow) = perm(gw_face_slot(j, irow))
         enddo
      enddo
   endif
   !
   write(*,'(a,i0,a,i0,a)') ' Semi-implicit CSR: ', nperm, ' of ', nrows_tot, ' rows reordered to lower/diag/upper'
   !
   deallocate(perm)
   deallocate(col_tmp)
   !
   ! Flatten into the per-row face list the assembly walks.
   !
   allocate(si_row_face_ptr(nrows_si + 1))
   si_row_face_ptr = 0
   !
   k = 0
   do irow = 1, nrows_si
      do j = 1, 8
         if (face_ip(j, irow) > 0) k = k + 1
      enddo
   enddo
   nfaces_si = k
   !
   allocate(si_row_face_ip(nfaces_si))
   allocate(si_row_face_slot(nfaces_si))
   allocate(si_row_face_isy(nfaces_si))
   allocate(si_row_face_bnd(nfaces_si))
   allocate(si_row_face_gwslot(nfaces_si))
   allocate(si_row_face_dinv(nfaces_si))
   si_row_face_gwslot = 0
   si_row_face_dinv = 0.0
   !
   if (gwflow) then
      allocate(gw_coeff_applied(nfaces_si))
      allocate(gw_cexch_applied(nrows_si))
      allocate(gw_qexpl_applied(nrows_si))
      gw_coeff_applied = 0.0
      gw_cexch_applied = 0.0
      gw_qexpl_applied = 0.0
      allocate(gw_cseep_applied(nrows_si))
      gw_cseep_applied = 0.0
   endif
   !
   k = 0
   do irow = 1, nrows_si
      si_row_face_ptr(irow) = k + 1
      do j = 1, 8
         if (face_ip(j, irow) > 0) then
            k = k + 1
            si_row_face_ip(k)   = face_ip(j, irow)
            si_row_face_slot(k) = face_slot(j, irow)   ! 0 means Dirichlet neighbour
            si_row_face_isy(k)  = face_isy(j, irow)
            si_row_face_bnd(k)  = face_bnd(j, irow)
            if (gwflow) si_row_face_gwslot(k) = gw_face_slot(j, irow)
            !
            ! Centre-to-centre distance for this face, matching what momentum uses.
            ! Same refinement level: the cell size at that level. Across a transition:
            ! 1.5 x the FINE cell size, which is coarse-half plus fine-half
            ! (sfincs_momentum.f90:239 computes exactly this).
            !
            ip = face_ip(j, irow)
            if (face_isy(j, irow) == 0) then
               dref = 1.0 / dxrinv(uv_flags_iref(ip))
            else
               dref = 1.0 / dyrinv(uv_flags_iref(ip))
            endif
            if (uv_flags_type(ip) /= 0) dref = 1.5 * dref
            si_row_face_dinv(k) = 1.0 / dref
            !
         endif
      enddo
   enddo
   si_row_face_ptr(nrows_si + 1) = k + 1
   !
   ! Check the CSR invariants the preconditioner relies on. It splits lower from upper by
   ! comparing the column index against the row index, so every row must contain exactly
   ! one entry on the diagonal and si_diag_ptr must point at it. If the row builder ever
   ! emits a duplicate column or misses the diagonal, SSOR would silently stop being a
   ! valid preconditioner rather than failing visibly.
   !
   maxrow = 0
   do irow = 1, nrows_tot
      maxrow = max(maxrow, si_row_ptr(irow + 1) - si_row_ptr(irow))
      ndiag = 0
      do k = si_row_ptr(irow), si_row_ptr(irow + 1) - 1
         if (si_col_idx(k) == irow) ndiag = ndiag + 1
      enddo
      if (ndiag /= 1) then
         write(*,*) 'Error: semi-implicit row ', irow, ' has ', ndiag, ' diagonal entries'
         stop
      endif
      if (si_col_idx(si_diag_ptr(irow)) /= irow) then
         write(*,*) 'Error: semi-implicit si_diag_ptr does not point at the diagonal, row ', irow
         stop
      endif
      do k = si_row_ptr(irow), si_diag_ptr(irow) - 1
         if (si_col_idx(k) >= irow) then
            write(*,*) 'Error: semi-implicit row ', irow, ' has an upper entry before the diagonal'
            stop
         endif
      enddo
      do k = si_diag_ptr(irow) + 1, si_row_ptr(irow + 1) - 1
         if (si_col_idx(k) <= irow) then
            write(*,*) 'Error: semi-implicit row ', irow, ' has a lower entry after the diagonal'
            stop
         endif
      enddo
   enddo
   !
   ! The maximum row length is the useful number: 5 on a regular grid, above 5 once
   ! quadtree connectivity is picked up (a coarse cell next to refinement reaches 9).
   ! If it stays at 5 on a quadtree model, the refinement transitions are not being seen.
   !
   write(*,'(a,i10,a,i10,a,i3,a,i10)') ' Semi-implicit: rows ', nrows_tot, &
      '  nonzeros ', nnz_si, '  max row ', maxrow, '  faces ', nfaces_si
   !
   deallocate(face_ip)
   deallocate(face_nb)
   deallocate(face_bnd)
   deallocate(face_isy)
   deallocate(face_slot)
   !
   end subroutine initialize_semi_implicit
   !
   !
   subroutine assemble_and_solve_pressure(dt)
   !
   ! Assemble the Helmholtz pressure matrix and RHS, then solve.
   !
   ! The system is:
   !   eta(nm) - theta^2 * g * dt^2 * sum_neighbors[ h_uv / (dx^2) * (eta_nb - eta_nm) ] = RHS
   !
   ! where RHS = zs(nm) - dt * div(q_star) + dt * sources
   !
   ! The matrix is SPD (symmetric positive definite).
   !
   implicit none
   !
   real*4, intent(in) :: dt
   !
   integer :: irow, nm, ip, iter
   integer :: nmd, nmu, ndm, num
   integer :: count0, count1, count_rate, count_max
   real*4  :: relres
   real*4  :: div_qstar
   integer :: kface, islot, iouter, nbad, nbad_prev, nbad_tol
   real*4  :: coeff_face
   real*4  :: acell, awet_n, awet_k, diag_store, dmax_outer, dchg, scale_i
   real*8  :: vol_n, vol_k
   integer :: jrow
   real*4  :: cexch, tface, gdvol, qexpl, cseep
   real*8  :: gvol_n, gvol_k, hk, resid, zceil, xi
   integer :: nmb
   real*4  :: diag
   real*4  :: dxr_val, dyr_val
   real*8  :: bv_rech, bv_exch, bv_bnd, bv_ceil, bv_gross, bv_term
   real*4  :: csum_face, cbnd_face
   real*8  :: rhs_c
   !
   call system_clock(count0, count_rate, count_max)
   !
   ! Zero the matrix and RHS
   !
   si_AA = 0.0
   si_rhs = 0.0
   !
   ! Initialize solution with current water levels as initial guess
   !
   !$omp parallel do private(irow) schedule(static)
   do irow = 1, nrows_si
      si_x(irow) = zs(si_nm_of_row(irow))
   enddo
   !$omp end parallel do
   !
   nbad_prev = huge(nbad)
   nbad_tol  = ceiling(si_outer_frac * nrows_tot)
   !
   ! Seed the outer iterate from the current water level. On the first pass the subgrid
   ! branch then reduces to exactly the linear form, so a subgrid model starts from the same
   ! place a non-subgrid one would.
   !
   ! Unconditionally, for EVERY model. This used to be guarded by (subgrid .or. gwflow),
   ! because si_eta_k was only a linearisation point and a linear model has nothing to
   ! linearise. It is now also the point the residual is expanded about -- si_rhs holds
   ! b - rowsum*si_eta_k -- so a plain model that left it at zero got a residual that was not
   ! a residual at all, and blew up to 1e16 m on the first Bates test.
   !
   do irow = 1, nrows_si
      si_eta_k(irow) = zs(si_nm_of_row(irow))
   enddo
   !
   if (gwflow) then
      !
      ! Snapshot the head at time level n BEFORE the outer loop starts.
      !
      ! gw_head itself is overwritten every outer iteration, because the transmissivity lag and
      ! the storage linearisation both read the latest iterate. The old-time storage term must
      ! NOT: if it reads gw_head it sees the current iterate, the storage residual cancels to
      ! zero, and the aquifer jumps to the steady state of the diffusion operator on the first
      ! timestep and then sits there with dh/dt identically zero.
      !
      ! Where the aquifer meets a prescribed-water-level boundary, its head is that water
      ! level. This lets an ordinary bzs time series drive a tidal aquifer boundary, instead of
      ! groundwater needing forcing machinery of its own. Off by default, because a boundary
      ! cell with no bzs forcing carries zs = zsini, which would overwrite a prescribed head.
      !
      if (gw_bnd_from_zs) then
         do nm = 1, np
            if (kcs(nm) == 2) gw_head(nm) = real(zs(nm))
         enddo
      endif
      !
      do nm = 1, np
         gw_head_n(nm) = gw_head(nm)
      enddo
      !
      do irow = 1, nrows_si
         si_eta_k(nrows_si + irow) = gw_head(si_nm_of_row(irow))
         si_x(nrows_si + irow)     = gw_head(si_nm_of_row(irow))
      enddo
      !
   endif
   !
   !
   ! Per-timestep assembly of everything that does NOT depend on the outer iterate.
   !
   ! si_coeff and si_q_star come from the momentum predictor and are fixed for the step, and so is
   ! the geometry. So every off-diagonal entry of a surface row, the face sum on its diagonal, the
   ! flux-divergence and source part of its right-hand side, and the known boundary levels are the
   ! same in every outer iteration. Measured on Harvey the outer loop ran 2.7 times per step to feed
   ! 14 CG iterations, and rebuilt all of this each time. It is built once here; the loop below adds
   ! only the storage term, which is the one thing the iterate changes.
   !
   ! The Dirichlet term coeff * (zs_bnd - eta_k) is split: coeff * zs_bnd goes into the constant
   ! part, and the row's total boundary conductance si_cbnd multiplies -eta_k inside the loop.
   !
   ! si_AA stays UNSCALED from here on. The diagonal scaling before each solve writes into
   ! si_AA_scaled instead of overwriting in place, which is what lets these off-diagonals survive
   ! from one outer iteration to the next.
   !
   !$omp parallel do private(irow, nm, nmd, nmu, ndm, num, div_qstar, acell, kface, ip, islot, &
   !$omp                      coeff_face, csum_face, cbnd_face, rhs_c, vol_n, awet_n) &
   !$omp schedule(static)
   do irow = 1, nrows_si
      !
      nm = si_nm_of_row(irow)
      !
      nmd = z_index_uv_md(nm)
      nmu = z_index_uv_mu(nm)
      ndm = z_index_uv_nd(nm)
      num = z_index_uv_nu(nm)
      !
      ! Flux divergence of q_star. Same formula as compute_water_levels_regular in
      ! sfincs_continuity.f90; inflow positive.
      !
      if (crsgeo) then
         div_qstar = (si_q_star(nmd) - si_q_star(nmu)) / dxm(nm) &
                   + (si_q_star(ndm) - si_q_star(num)) * dyrinv(z_flags_iref(nm))
         acell = cell_area_m2(nm)
      else
         div_qstar = (si_q_star(nmd) - si_q_star(nmu)) * dxrinv(z_flags_iref(nm)) &
                   + (si_q_star(ndm) - si_q_star(num)) * dyrinv(z_flags_iref(nm))
         acell = cell_area(z_flags_iref(nm))
      endif
      !
      ! Volumetric, like every other term in the row.
      !
      rhs_c = dble(acell) * dble(dt) * dble(div_qstar)
      if (precip)   rhs_c = rhs_c + acell * dt * netprcp(nm)
      if (use_qext) rhs_c = rhs_c + acell * dt * qext(nm)
      !
      ! Walk this row's faces. Grid-agnostic: nothing here assumes there are four.
      ! si_coeff is per unit face width, so multiply by the width of THIS face, taken from the
      ! face's own refinement level. At a quadtree transition the coarse cell's two fine faces get
      ! half the coarse width each, which keeps A(coarse,fine) equal to A(fine,coarse).
      !
      csum_face = 0.0
      cbnd_face = 0.0
      !
      do kface = si_row_face_ptr(irow), si_row_face_ptr(irow + 1) - 1
         !
         ip    = si_row_face_ip(kface)
         islot = si_row_face_slot(kface)
         !
         if (si_row_face_isy(kface) == 0) then
            coeff_face = si_coeff(ip) * dt * dyrm(uv_flags_iref(ip))
         else
            coeff_face = si_coeff(ip) * dt * dxrm(uv_flags_iref(ip))
         endif
         !
         csum_face = csum_face + coeff_face
         !
         if (islot > 0) then
            ! Interior neighbour: the off-diagonal entry, constant for the step
            si_AA(islot) = -coeff_face
         else
            ! Boundary neighbour (kcs==2): known level to the constant right-hand side, and the
            ! conductance remembered so the loop can subtract coeff * eta_k against it
            cbnd_face = cbnd_face + coeff_face
            rhs_c     = rhs_c + dble(coeff_face) * dble(zs(si_row_face_bnd(kface)))
         endif
         !
      enddo
      !
      si_cface(irow)     = csum_face
      si_cbnd(irow)      = cbnd_face
      si_rhs_const(irow) = rhs_c
      !
      ! Subgrid volume at the OLD level, needed by the storage residual every iteration
      !
      if (subgrid) then
         call subgrid_storage(nm, zs(nm), vol_n, awet_n)
         si_vol_n(irow) = vol_n
      endif
      !
   enddo
   !$omp end parallel do
   !
   do iouter = 1, si_maxouter
   !
   ! Assemble the iterate-dependent part of the matrix and RHS row by row
   !
   !$omp parallel do private(irow, nm, nmd, nmu, ndm, num, dxr_val, dyr_val, &
   !$omp                      div_qstar, diag, ip, kface, islot, coeff_face, &
   !$omp                      acell, vol_n, vol_k, awet_n, awet_k, diag_store, cexch, qexpl, &
   !$omp                      cseep, zceil) &
   !$omp schedule(static)
   do irow = 1, nrows_si
      !
      nm = si_nm_of_row(irow)
      !
      ! Everything that does not depend on the outer iterate -- off-diagonals, the face sum
      ! si_cface, the boundary conductance si_cbnd and the constant right-hand side si_rhs_const --
      ! was assembled once for this timestep before the loop. What is left is the storage term.
      !
      if (crsgeo) then
         acell = cell_area_m2(nm)
      else
         acell = cell_area(z_flags_iref(nm))
      endif
      !
      ! Storage term.
      !
      ! Without subgrid, V = A*eta so the storage contributes A on the diagonal and A*eta^n on
      ! the RHS, which is what the two branches below reduce to.
      !
      ! With subgrid the relation is nonlinear. Expand it about the current outer iterate:
      !
      !   V(eta^{n+1}) ~= V(eta^k) + A_wet(eta^k) * (eta^{n+1} - eta^k)
      !
      ! This is Newton written in residual/increment form. It is NOT the same as lagging A_wet as
      ! a multiplicative coefficient, which is a period-2 limit cycle at every timestep. The
      ! nonlinearity is diagonal-only, so the matrix stays symmetric and CG is retained.
      !
      ! si_rhs holds the DIAGONAL-REDUCED RESIDUAL, not the right-hand side -- see the header of
      ! the residual loop below. The storage term contributes awet_k * eta^k to the right-hand
      ! side and awet_k to the diagonal, so the two cancel exactly and what is left is
      ! vol_n - vol_k: the volume the cell actually has to shed this step. The Dirichlet faces
      ! reduce the same way: their known levels sit in si_rhs_const and their conductance times
      ! this row's iterate is subtracted here.
      !
      ! Volumetric: every term is a volume per timestep, not a level. Diagonal scaling before the
      ! solve removes the resulting O(1e5) magnitudes, so CG is unaffected.
      !
      if (subgrid) then
         call subgrid_storage(nm, si_eta_k(irow), vol_k, awet_k)
         diag_store   = awet_k
         si_rhs(irow) = (si_vol_n(irow) - vol_k) + si_rhs_const(irow) &
                      - dble(si_cbnd(irow)) * si_eta_k(irow)
      else
         diag_store   = acell
         si_rhs(irow) = dble(acell) * (zs(nm) - si_eta_k(irow)) + si_rhs_const(irow) &
                      - dble(si_cbnd(irow)) * si_eta_k(irow)
      endif
      !
      ! Set diagonal: storage derivative plus the per-step face sum
      !
      si_AA(si_diag_ptr(irow)) = diag_store + si_cface(irow)
      !
      ! Exchange with the aquifer, written symmetrically.
      !
      ! The SAME conductance goes into A(i, N+i) and A(N+i, i) and onto both diagonals. A spike
      ! measured the coupled matrix at exactly 0.000e+00 symmetry residual with this
      ! construction; anything one-sided leaves CG converging confidently to a wrong answer
      ! rather than failing.
      !
      if (gwflow) then
         call gw_exchange_terms(nm, si_eta_k(irow), si_eta_k(nrows_si + irow), cexch, qexpl)
         cexch = cexch * dt
         si_AA(si_exch_ptr(irow)) = -cexch
         si_AA(si_diag_ptr(irow)) = si_AA(si_diag_ptr(irow)) + cexch
         si_rhs(irow) = si_rhs(irow) - dble(dt) * dble(qexpl)
         !
         ! Seepage face, arriving as surface water. LAGGED, not implicit.
         !
         ! Q_seep depends on the aquifer head but not on the surface level, so writing it
         ! implicitly on both rows would give A(i, N+i) = -(cexch + cseep)*dt against
         ! A(N+i, i) = -cexch*dt. That asymmetry costs us CG, and a spike measured the coupled
         ! matrix at exactly 0.000e+00 symmetry residual with the current construction. So the
         ! aquifer keeps the implicit half, where stability needs it, and the surface takes the
         ! lagged half, where it only needs to be conservative. At outer convergence the two are
         ! the same number; if the loop exits on stagnation they differ, and the water balance
         ! reports the difference rather than hiding it.
         !
         ! si_rhs is volumetric here, like every other term in this row, so a discharge times dt
         ! is the right thing to add.
         !
         call gw_seepage_terms(nm, si_eta_k(irow), si_eta_k(nrows_si + irow), dt, cseep, zceil)
         si_rhs(irow) = si_rhs(irow) &
                      + dble(cseep) * dble(dt) * (si_eta_k(nrows_si + irow) - zceil)
      endif
      !
   enddo
   !$omp end parallel do
   !
   ! Aquifer rows.
   !
   !    Sy dh/dt = div( K b grad h ) + R + exchange
   !
   ! written volumetrically and expanded about the current outer iterate, so the storage term
   ! contributes its derivative to the diagonal and the residual to the right-hand side. The
   ! transmissivity K*b is LAGGED from the previous iterate (Picard), which keeps the aquifer
   ! block symmetric and lets CG stay. That lag is safe here because at SFINCS timesteps the head
   ! moves ~1e-4 m per step; the spike measured roughly 150x margin before it fails.
   !
   if (gwflow) then
      !
      !$omp parallel do private(irow, nm, jrow, kface, ip, islot, acell, tface, coeff_face, &
      !$omp                     diag, cexch, gvol_n, gvol_k, gdvol, hk, nmb, qexpl, &
      !$omp                     cseep, zceil) schedule(static)
      do irow = 1, nrows_si
         !
         nm   = si_nm_of_row(irow)
         jrow = nrows_si + irow
         hk   = si_eta_k(jrow)
         !
         if (crsgeo) then
            acell = cell_area_m2(nm)
         else
            acell = cell_area(z_flags_iref(nm))
         endif
         !
         ! gvol_n is the storage at time n and takes the time-n ceiling, which is what zs(nm)
         ! still holds. gvol_k is the new iterate and takes the new iterate's surface level: the
         ! ceiling moves with the pond, and a cell whose table is pinned to the ground stores
         ! Sy*A more per metre the pond deepens. Freezing the ceiling at zs^n instead ejects that
         ! water as seepage and then finds it back in the storage a step later, unaccounted.
         !
         call gw_cell_storage(nm, gw_head_n(nm), gvol_n, gdvol)
         call gw_cell_storage(nm, hk, gvol_k, gdvol, si_eta_k(irow))
         !
         diag = gdvol
         ! gdvol * hk against gdvol on the diagonal: cancels, leaving gvol_n - gvol_k.
         !
         si_rhs(jrow) = (gvol_n - gvol_k) &
                      + dble(acell) * dble(dt) * dble(gw_recharge(nm))
         !
         do kface = si_row_face_ptr(irow), si_row_face_ptr(irow + 1) - 1
            !
            ip    = si_row_face_ip(kface)
            islot = si_row_face_gwslot(kface)
            !
            ! The cell on the other side of this face, whether or not it is an unknown.
            !
            nmb = uv_index_z_nm(ip)
            if (nmb == nm) nmb = uv_index_z_nmu(ip)
            !
            call gw_face_transmissivity(ip, tface)
            !
            if (si_row_face_isy(kface) == 0) then
               coeff_face = tface * dyrm(uv_flags_iref(ip)) * si_row_face_dinv(kface) * dt
            else
               coeff_face = tface * dxrm(uv_flags_iref(ip)) * si_row_face_dinv(kface) * dt
            endif
            !
            gw_coeff_applied(kface) = coeff_face
            !
            ! theta-weighted in time: theta on the new level, (1 - theta) on the old.
            !
            ! The explicit half is NOT optional. Leaving it out integrates
            ! Sy dh/dt = theta * L(h^{n+1}), so the aquifer runs at an effective diffusivity of
            ! theta*D. A steady case cannot see this -- theta multiplies both sides and cancels,
            ! which is why the Dupuit parabola came out right to 0.003% while the transient was
            ! 25% slow. The Edelman step response fits D_eff/D = 0.75 = theta exactly.
            !
            diag = diag + gw_theta * coeff_face
            !
            si_rhs(jrow) = si_rhs(jrow) &
                         + dble((1.0 - gw_theta) * coeff_face) * (gw_head_n(nmb) - gw_head_n(nm))
            !
            if (islot > 0) then
               si_AA(islot) = -gw_theta * coeff_face
            else
               si_rhs(jrow) = si_rhs(jrow) &
                            + dble(gw_theta * coeff_face) * (gw_head(nmb) - hk)
            endif
            !
         enddo
         !
         call gw_exchange_terms(nm, si_eta_k(irow), hk, cexch, qexpl)
         cexch = cexch * dt
         si_AA(si_exch_ptr(jrow)) = -cexch
         diag = diag + cexch
         si_rhs(jrow) = si_rhs(jrow) + dble(dt) * dble(qexpl)
         !
         gw_cexch_applied(irow) = cexch
         gw_qexpl_applied(irow) = qexpl
         !
         ! Seepage face. Implicit here, because this is where stability needs it: cseep*dt sits on
         ! the diagonal and grows the outflow as the head rises, which is what stops the runaway.
         ! The matching term on the surface row is lagged -- see the surface loop for why.
         !
         call gw_seepage_terms(nm, si_eta_k(irow), hk, dt, cseep, zceil)
         diag = diag + cseep * dt
         ! cseep*dt*zceil against cseep*dt on the diagonal: cancels to the head's excess
         ! above the ceiling, which is what the seepage face actually responds to.
         !
         si_rhs(jrow) = si_rhs(jrow) - dble(cseep) * dble(dt) * (hk - zceil)
         !
         gw_cseep_applied(irow) = cseep * dt
         !
         si_AA(si_diag_ptr(jrow)) = diag
         !
      enddo
      !$omp end parallel do
      !
   endif
   !
   ! Symmetric diagonal scaling: solve (D A D) y = D b with D = diag(1/sqrt(diag A)),
   ! then recover x = D y. Every scaled diagonal becomes exactly 1, so CG sees a system
   ! whose conditioning no longer depends on the spread between dry and wet cells.
   !
   !$omp parallel do private(irow) schedule(static)
   do irow = 1, nrows_tot
      si_scale(irow) = 1.0 / sqrt(max(si_AA(si_diag_ptr(irow)), 1.0e-20))
   enddo
   !$omp end parallel do
   !
   ! Solve for the INCREMENT, and never form the datum.
   !
   ! CG's stopping test is ||r|| / ||b||, and if the unknown is the absolute level then ||b|| is
   ! dominated by rows that carry no information: a dry cell contributes acell * zs, so a domain
   ! sitting 50 m above datum inflates ||b|| by three orders of magnitude and the tolerance stops
   ! constraining anything. With groundwater that is fatal -- the aquifer's per-step head change
   ! is ~6e-5 m while the effective tolerance would permit ~2e-3 m, so CG returns "converged"
   ! after zero iterations and the water table freezes mid-transient.
   !
   ! Solving A dx = b - A x0 fixes the tolerance, but computing that difference NUMERICALLY does
   ! not get the datum out. Written out,
   !
   !    b - A x0  =  [ b - rowsum(i) * x0(i) ]  -  sum_j A_ij * ( x0(j) - x0(i) )
   !
   ! the second sum is harmless: every coefficient multiplies a level DIFFERENCE between
   ! neighbours, millimetres at most, so a real*4 coefficient's ~6e-8 relative error is scaled by
   ! the gradient. The bracket is where the datum lives, and it cannot be rescued by widening
   ! anything: rowsum is an accumulation of real*4 coefficients, x0(i) carries the datum, and one
   ! rounding of the diagonal against a 10 m datum is already 6e-7 m against a residual whose
   ! true size is ~5e-6 m. Widening si_AA to real*8 was tried and did NOT fix it -- it improved
   ! the datum-free cases and left the datum ones worse, which is the signature of paying for
   ! precision instead of removing a cancellation.
   !
   ! So the bracket is never computed. It cancels analytically, term by term: the storage term
   ! contributes gdvol*h^k to b and gdvol to the diagonal; each Dirichlet coefficient puts its
   ! known level in b and itself on the diagonal; the seepage conductance does the same. Every
   ! one of those pairs reduces to a level difference, and the assembly above writes si_rhs as
   ! the already-reduced quantity. si_rhs is therefore the RESIDUAL, not the right-hand side.
   !
   ! What is left here is the neighbour sum, and the exchange coefficient, which cancels between
   ! the diagonal and the coupling entry and so needs no special treatment.
   !
   ! Done BEFORE the diagonal scaling, because the scaling breaks the difference form: the scaled
   ! unknown is x(j)/s(j) with s varying row to row, so x_scaled(j) - x_scaled(i) is no longer a
   ! level difference. Form the residual unscaled, then scale the one number that comes out.
   !
   ! The scaled matrix row is written in the same pass: it needs exactly the entries the
   ! residual just touched (si_AA, si_col_idx, si_scale of row and column), so the row's data
   ! is in cache and the matrix is walked once instead of twice. Into a separate array:
   ! si_AA must stay unscaled, because its surface off-diagonals are assembled once per
   ! timestep and reused by every outer iteration. The loop only reads si_x, so writing
   ! si_AA_scaled here is race-free.
   !
   !$omp parallel do private(irow, kface, resid, xi, scale_i) schedule(static)
   do irow = 1, nrows_tot
      xi      = si_x(irow)
      scale_i = si_scale(irow)
      resid   = si_rhs(irow)
      do kface = si_row_ptr(irow), si_row_ptr(irow + 1) - 1
         resid = resid - dble(si_AA(kface)) * (si_x(si_col_idx(kface)) - xi)
         si_AA_scaled(kface) = si_AA(kface) * scale_i * si_scale(si_col_idx(kface))
      enddo
      !
      si_b(irow)  = real(resid * dble(scale_i))
      si_dx(irow) = 0.0
   enddo
   !$omp end parallel do
   !
   ! Solve using CG with SSOR preconditioning
   !
   !
   call cg_solve(nrows_tot, nnz_si, si_AA_scaled, si_col_idx, si_row_ptr, &
                  si_b, si_dx, si_tol, si_maxiter, iter, relres)
   !
   ! Add the increment, unscaled. CG solved for dx in the scaled unknown x/s, so the level
   ! change is dx*s.
   !
   ! si_x used to be round-tripped through si_scale (x -> x/s -> (x/s + dx)*s). When si_x was
   ! real*4 the two forms differed measurably (Dupuit slightly better, Edelman clearly worse
   ! without the round trip) and both were dominated by real*4 noise. si_x is real*8 now, so
   ! the difference is 2e-16 relative, and the round trip cost a serial pass over si_x that
   ! the residual loop could not be fused with (it reads neighbours' si_x).
   !
   !$omp parallel do private(irow) schedule(static)
   do irow = 1, nrows_tot
      si_x(irow) = si_x(irow) + dble(si_dx(irow)) * dble(si_scale(irow))
   enddo
   !$omp end parallel do
   !
   si_iter_total = si_iter_total + iter
   si_solve_count = si_solve_count + 1
   si_iter_max_seen = max(si_iter_max_seen, iter)
   !
   ! Advance the expansion point. Unconditional for the same reason it is seeded
   ! unconditionally: the residual is expanded about si_eta_k whether or not anything here is
   ! nonlinear. For a linear model the loop exits after this one pass, so it costs a copy.
   !
   ! The change is measured BEFORE the copy -- it is the distance between the new solution and
   ! the point the residual was expanded about, which is exactly what the convergence test below
   ! wants, and it is identically zero if the copy happens first.
   !
   dmax_outer = 0.0
   nbad = 0
   !$omp parallel do private(irow, dchg) reduction(+:nbad) reduction(max:dmax_outer) schedule(static)
   do irow = 1, nrows_tot
      dchg = real(abs(si_x(irow) - si_eta_k(irow)))
      dmax_outer = max(dmax_outer, dchg)
      if (dchg > si_tolouter) nbad = nbad + 1
   enddo
   !$omp end parallel do
   !
   si_prof_count(iouter) = si_prof_count(iouter) + 1
   si_prof_nbad(iouter)  = si_prof_nbad(iouter) + nbad
   si_prof_dmax(iouter)  = si_prof_dmax(iouter) + dmax_outer
   !
   si_eta_k(1:nrows_tot) = si_x(1:nrows_tot)
   !
   ! Nonlinear outer iteration, subgrid only. Without subgrid the system is linear and one
   ! pass is exact, so this costs nothing there.
   !
   if (subgrid .or. gwflow) then
      !
      ! Feed the aquifer head back so the next iterate relags transmissivity and storage.
      !
      if (gwflow) then
         do irow = 1, nrows_si
            gw_head(si_nm_of_row(irow)) = si_x(nrows_si + irow)
         enddo
         !
      endif
      !
      si_outer_total = si_outer_total + 1
      !
      ! Stop when the bulk of the field has converged, or when that stops improving.
      !
      ! The max-norm is not a usable test here: a handful of cells sitting on a wet/dry
      ! threshold keep flipping between two states, so the maximum change never falls below
      ! tolerance even though the field has converged. Measured on Harvey, a pure max-norm
      ! test ran to the 50-iteration cap on 13479 of 13480 timesteps while iterations 4-50
      ! changed the gauge RMSE by less than 0.01 m -- a factor 10.5 in runtime for nothing.
      ! A 10 %-improvement test on that same max-norm was tried next and turned out to be
      ! decided by rounding: the maximum is set by whichever threshold cell flipped hardest,
      ! an O(1 m) quantity that is random from one iteration to the next, so the outer count
      ! (and the answer) changed with the thread count and with any reordering of the sums.
      !
      ! So judge the bulk instead: nbad is the number of rows still moving by more than
      ! si_tolouter. Converged when nbad is down to a fraction si_outer_frac of the rows.
      ! Stop as well when an iteration fails to at least halve nbad: what is left then is
      ! the population that does not contract at all. Measured on Harvey (profile in the
      ! log): 38 % of rows moving after iteration 1, 6.6 % after 2, 5.0 % after 3, then a
      ! plateau of ~4.3 % that flips by ~0.25 m every iteration for as long as the loop
      ! runs -- cells at a kink of the subgrid storage curve, on which the Newton update
      ! oscillates. Iterating on them changes the gauges by nothing (RMSE identical to four
      ! decimals between 2 and 11 iterations) and costs a full solve each. An integer count
      ! compared against half of itself is only perturbed by rounding when the count sits
      ! within one row of the threshold, not whenever the worst cell wobbles.
      !
      if (iouter > 1 .and. 2 * nbad > nbad_prev) then
         si_outer_stagnant = si_outer_stagnant + 1
      else
         nbad_prev = nbad
         if (nbad > nbad_tol .and. iouter < si_maxouter) cycle
      endif
      !
      si_outer_max_seen = max(si_outer_max_seen, iouter)
      si_solve_count_outer = si_solve_count_outer + 1
      if (iouter >= si_maxouter) si_outer_capped = si_outer_capped + 1
      !
   endif
   !
   exit
   !
   enddo
   !
   ! Aquifer water budget for this timestep.
   !
   ! The implicit path's fluxes live inside the matrix, so they have to be recovered from the
   ! converged heads by one extra pass over the same face list the assembly walks. That pass is
   ! cheap, and it is the only way to be sure the budget measures what the solver did rather
   ! than what it was supposed to do.
   !
   ! Three things this has to match exactly or the closure number means nothing:
   !   - the surface level is si_eta_k(irow), NOT zs(nm). zs is not updated until
   !     backsubstitute_fluxes_si runs, so reading it here would use last step's level.
   !   - the lateral flux is theta-weighted between gw_head_n and gw_head, the same split the
   !     assembly applies. Using only the new-time heads leaves a (1 - theta) share of every
   !     lateral flux unaccounted, which on Edelman is the whole signal.
   !   - a face is a boundary face when its aquifer slot is zero, which is precisely the branch
   !     where the assembly moved the neighbour's head to the right-hand side.
   !
   ! There is no ceiling term: the implicit path enforces the ceiling through the storage
   ! relation rather than by moving a discrete excess volume to the surface.
   !
   if (gwflow) then
      !
      bv_rech  = 0.0d0
      bv_exch  = 0.0d0
      bv_bnd   = 0.0d0
      bv_ceil  = 0.0d0
      bv_gross = 0.0d0
      !
      do irow = 1, nrows_si
         !
         nm = si_nm_of_row(irow)
         !
         if (crsgeo) then
            acell = cell_area_m2(nm)
         else
            acell = cell_area(z_flags_iref(nm))
         endif
         !
         bv_term  = dble(acell) * dble(gw_recharge(nm)) * dble(dt)
         bv_rech  = bv_rech + bv_term
         bv_gross = bv_gross + abs(bv_term)
         !
         ! Exchange, with the conductance and the lagged remainder the matrix carried, applied to
         ! the converged levels. si_eta_k(irow) is the new surface level: zs is not updated until
         ! backsubstitute_fluxes_si runs, so reading zs here would use last step's level.
         !
         bv_term  = dble(gw_cexch_applied(irow)) &
                  * (dble(si_eta_k(irow)) - dble(gw_head(nm))) &
                  + dble(gw_qexpl_applied(irow)) * dble(dt)
         bv_exch  = bv_exch + bv_term
         bv_gross = bv_gross + abs(bv_term)
         !
         ! Seepage out of the aquifer at the ceiling. Signed as every other term is: negative
         ! because it LEAVES. Applied CONDUCTANCE against CONVERGED levels, exactly as the
         ! exchange term above -- the on/off switch the matrix made is honoured through
         ! gw_cseep_applied, and the ceiling is re-evaluated at the level the surface actually
         ! reached.
         !
         ! The ceiling has to be the converged one, not the lagged one the assembly used. A
         ! saturated cell under a deepening pond stores Sy*A more per metre of pond, so the row's
         ! storage term and its seepage term BOTH shift when the ceiling moves, by Sy*A*dzs each,
         ! and the two shifts cancel in the head -- which is why the ceiling case gets the right
         ! water level either way. They do not cancel in the budget: gw_total_storage measures
         ! against the new ceiling, so the seepage has to as well. Measured against the lagged
         ! ceiling it over-reports by Sy*A*dzs per step, which on the ceiling case is 77.3 m3 over
         ! the run and reads as a 5.8 % closure error against a state that is in fact correct.
         !
         call gw_seepage_terms(nm, si_eta_k(irow), gw_head(nm), dt, cseep, zceil)
         bv_term  = -dble(gw_cseep_applied(irow)) * (dble(gw_head(nm)) - dble(zceil))
         bv_ceil  = bv_ceil + bv_term
         bv_gross = bv_gross + abs(bv_term)
         !
         ! Lateral flux across faces whose neighbour is not an unknown -- exactly the branch where
         ! the assembly moved the neighbour's head to the right-hand side. Theta-weighted the same
         ! way the assembly weighted it: using only the new-time heads would leave a (1 - theta)
         ! share of every boundary flux unaccounted, which on Edelman is the whole signal.
         !
         do kface = si_row_face_ptr(irow), si_row_face_ptr(irow + 1) - 1
            !
            if (si_row_face_gwslot(kface) /= 0) cycle
            !
            ip  = si_row_face_ip(kface)
            nmb = uv_index_z_nm(ip)
            if (nmb == nm) nmb = uv_index_z_nmu(ip)
            !
            bv_term = dble(gw_coeff_applied(kface)) &
                    * (dble(gw_theta) * (dble(gw_head(nmb)) - dble(gw_head(nm))) &
                     + (1.0d0 - dble(gw_theta)) &
                     * (dble(gw_head_n(nmb)) - dble(gw_head_n(nm))))
            bv_bnd   = bv_bnd + bv_term
            bv_gross = bv_gross + abs(bv_term)
            !
         enddo
         !
      enddo
      !
      call gw_budget_add(bv_rech, bv_exch, bv_bnd, bv_ceil, bv_gross)
      !
   endif
   !
   call system_clock(count1, count_rate, count_max)
   tloop_si = tloop_si + 1.0 * (count1 - count0) / count_rate
   !
   end subroutine assemble_and_solve_pressure
   !
   !
   subroutine subgrid_storage(nm, eta, vol, awet)
   !
   ! Volume and wetted area of a subgrid cell at water level eta.
   !
   ! Mirrors the table inversion in sfincs_ncoutput.F90:4160-4181. The subgrid table stores
   ! level as a function of volume, binned uniformly in volume, so awet = dV/deta is
   ! piecewise CONSTANT and jumps at bin edges. That discontinuous derivative is why the
   ! outer iteration below has to be written in residual form: lagging awet as a
   ! multiplicative coefficient does not converge at any timestep.
   !
   implicit none
   !
   integer, intent(in)  :: nm
   real*8,  intent(in)  :: eta
   real*8,  intent(out) :: vol
   real*4,  intent(out) :: awet
   !
   integer :: ilevel, ivol
   real*4  :: acell, dzvol, dz, facint, zmn, zmx
   !
   if (crsgeo) then
      acell = cell_area_m2(nm)
   else
      acell = cell_area(z_flags_iref(nm))
   endif
   !
   zmn = max(subgrid_z_zmin(nm), -20.0)
   zmx = max(subgrid_z_zmax(nm), -20.0)
   !
   if (eta >= zmx) then
      !
      ! Cell fully wet: storage grows with the full cell area
      !
      vol  = subgrid_z_volmax(nm) + acell * (eta - zmx)
      awet = acell
      !
   elseif (eta <= zmn) then
      !
      ! Cell dry. Use the wet area of the FIRST table bin: that is the area the first water
      ! entering the cell would occupy, and it is continuous with the interior branch below.
      !
      ! An earlier version used a token 1e-6*acell here to keep the row non-singular. That
      ! put a millionfold spread on the diagonal, CG hit its iteration ceiling every step,
      ! and the solution never moved off its initial guess -- gauges stayed frozen at their
      ! starting level through the whole event.
      !
      dzvol = subgrid_z_volmax(nm) / (subgrid_nlevels - 1)
      dz    = max(subgrid_z_dep(2, nm) - subgrid_z_dep(1, nm), 0.001)
      vol   = 0.0
      awet  = dzvol / dz
      !
   else
      !
      ivol = 1
      do ilevel = 2, subgrid_nlevels
         if (subgrid_z_dep(ilevel, nm) > eta) then
            ivol = ilevel - 1
            exit
         endif
      enddo
      !
      dzvol  = subgrid_z_volmax(nm) / (subgrid_nlevels - 1)
      dz     = max(subgrid_z_dep(ivol + 1, nm) - subgrid_z_dep(ivol, nm), 0.001)
      facint = (eta - subgrid_z_dep(ivol, nm)) / dz
      vol    = (ivol - 1) * dzvol + facint * dzvol
      awet   = dzvol / dz
      !
   endif
   !
   ! Keep the derivative inside a sane range. A flat bench in the table makes dz tiny and
   ! awet enormous; a steep one makes it vanish. Either wrecks the conditioning of the
   ! coupled matrix. Clamping is safe because awet is only the linearisation derivative --
   ! the residual form carries V(eta) exactly on the RHS, so this changes the convergence
   ! rate of the outer loop, not the solution it converges to.
   !
   awet = min(max(awet, awet_floor * acell), acell)
   !
   end subroutine subgrid_storage
   !
   !
   subroutine cg_solve(n, nnz, val, col_ind, row_ptr, b, x, tol, maxiter, iter, relres)
   !
   ! Conjugate Gradient solver with SSOR preconditioning.
   !
   ! For SPD systems. Uses the known 5-point stencil structure
   ! via si_diag_ptr for efficient SSOR forward/backward sweeps.
   !
   ! All scalar accumulators use double precision for robustness
   ! with large cell counts (>1M cells).
   !
   implicit none
   !
   integer, intent(in)    :: n, nnz, col_ind(nnz), row_ptr(n+1), maxiter
   real*4,  intent(in)    :: val(nnz), b(n), tol
   real*4,  intent(inout) :: x(n)
   integer, intent(out)   :: iter
   real*4,  intent(out)   :: relres
   !
   integer  :: i, k, idx
   real*8   :: rz, rz_new, pAp, bnorm2, rnorm2
   real*4   :: alpha, beta, tmp
   real*4   :: omega
   !
   omega = 1.5  ! SSOR relaxation parameter (1.0 = SGS, 1.5 = typical SSOR)
   !
   ! Extract diagonal
   !
   !$omp parallel do private(i) schedule(static)
   do i = 1, n
      cg_diag(i) = val(si_diag_ptr(i))
   enddo
   !$omp end parallel do
   !
   ! Initial residual: r = b - A*x (using stencil structure, not CSR traversal)
   !
   bnorm2 = 0.0d0
   rnorm2 = 0.0d0
   !$omp parallel do private(i, k, idx, tmp) reduction(+:bnorm2, rnorm2) schedule(static)
   do i = 1, n
      tmp = b(i)
      do k = row_ptr(i), row_ptr(i + 1) - 1
         tmp = tmp - val(k) * x(col_ind(k))
      enddo
      cg_r(i) = tmp
      bnorm2 = bnorm2 + dble(b(i)) * dble(b(i))
      rnorm2 = rnorm2 + dble(tmp) * dble(tmp)
   enddo
   !$omp end parallel do
   !
   bnorm2 = max(bnorm2, 1.0d-60)
   relres = real(sqrt(rnorm2 / bnorm2))
   !
   if (relres <= tol) then
      iter = 0
      return
   endif
   !
   ! Apply SSOR preconditioner: z = M^{-1} r
   !
   call apply_ssor_precond(n, val, cg_r, cg_z, omega)
   !
   ! p = z, rz = r.z
   !
   rz = 0.0d0
   !$omp parallel do private(i) reduction(+:rz) schedule(static)
   do i = 1, n
      cg_p(i) = cg_z(i)
      rz = rz + dble(cg_r(i)) * dble(cg_z(i))
   enddo
   !$omp end parallel do
   !
   ! CG iterations
   !
   do iter = 1, maxiter
      !
      ! Ap = A * p and pAp = p . Ap (fused)
      !
      pAp = 0.0d0
      !$omp parallel do private(i, k) reduction(+:pAp) schedule(static)
      do i = 1, n
         cg_Ap(i) = 0.0
         do k = row_ptr(i), row_ptr(i + 1) - 1
            cg_Ap(i) = cg_Ap(i) + val(k) * cg_p(col_ind(k))
         enddo
         pAp = pAp + dble(cg_p(i)) * dble(cg_Ap(i))
      enddo
      !$omp end parallel do
      !
      if (pAp <= 0.0d0) exit
      alpha = real(rz / pAp)
      !
      ! x += alpha*p, r -= alpha*Ap, compute rnorm2
      !
      rnorm2 = 0.0d0
      !$omp parallel do private(i) reduction(+:rnorm2) schedule(static)
      do i = 1, n
         x(i) = x(i) + alpha * cg_p(i)
         cg_r(i) = cg_r(i) - alpha * cg_Ap(i)
         rnorm2 = rnorm2 + dble(cg_r(i)) * dble(cg_r(i))
      enddo
      !$omp end parallel do
      !
      relres = real(sqrt(rnorm2 / bnorm2))
      if (relres <= tol) exit
      !
      ! Apply SSOR preconditioner: z = M^{-1} r
      !
      call apply_ssor_precond(n, val, cg_r, cg_z, omega)
      !
      ! rz_new = r . z
      !
      rz_new = 0.0d0
      !$omp parallel do private(i) reduction(+:rz_new) schedule(static)
      do i = 1, n
         rz_new = rz_new + dble(cg_r(i)) * dble(cg_z(i))
      enddo
      !$omp end parallel do
      !
      if (rz <= 0.0d0) exit
      beta = real(rz_new / rz)
      rz = rz_new
      !
      !$omp parallel do private(i) schedule(static)
      do i = 1, n
         cg_p(i) = cg_z(i) + beta * cg_p(i)
      enddo
      !$omp end parallel do
      !
   enddo
   !
   end subroutine cg_solve
   !
   !
   subroutine apply_ssor_precond(n, val, r, z, omega)
   !
   ! Symmetric Successive Over-Relaxation (SSOR) preconditioner.
   !
   ! Works on the CSR structure directly. Every row is stored as [lower | diagonal | upper]
   ! (enforced once at setup and checked there), so the lower part of row i is the slots
   ! si_row_ptr(i) .. si_diag_ptr(i)-1 and the upper part si_diag_ptr(i)+1 .. si_row_ptr(i+1)-1.
   ! No assumption about how many entries a row has, which is what a quadtree row needs.
   !
   ! SSOR = (D/omega + L) * D^{-1} * (D/omega + U)
   !
   ! Forward sweep:  (D/omega + L) * y = r
   ! Backward sweep: (D/omega + U) * z = D/omega * y
   !
   ! Both sweeps are inherently sequential. For the Helmholtz system SSOR typically
   ! reduces CG iterations by 3-5x compared to Jacobi.
   !
   implicit none
   !
   integer, intent(in)  :: n
   real*4,  intent(in)  :: val(*), r(n), omega
   real*4,  intent(out) :: z(n)
   !
   integer :: i, k
   real*4  :: diag_i, tmp
   !
   ! Forward sweep: (D/omega + L) * z = r
   ! Ascending k over the lower slots only.
   !
   do i = 1, n
      !
      diag_i = val(si_diag_ptr(i))
      tmp = r(i)
      !
      do k = si_row_ptr(i), si_diag_ptr(i) - 1
         tmp = tmp - val(k) * z(si_col_idx(k))
      enddo
      !
      z(i) = omega * tmp / diag_i
      !
   enddo
   !
   ! Scale: z = D/omega * z (prepare for backward sweep)
   !
   !$omp parallel do private(i) schedule(static)
   do i = 1, n
      z(i) = val(si_diag_ptr(i)) / omega * z(i)
   enddo
   !$omp end parallel do
   !
   ! Backward sweep: (D/omega + U) * z_new = z_old
   ! Descending k over the upper slots only. Descending rather than ascending on purpose:
   ! with the historical emission order (left, bottom, centre, top, right) it visits right
   ! before top, which is the order the old slot-indexed version used, so the sum rounds
   ! identically.
   !
   do i = n, 1, -1
      !
      diag_i = val(si_diag_ptr(i))
      tmp = z(i)
      !
      do k = si_row_ptr(i + 1) - 1, si_diag_ptr(i) + 1, -1
         tmp = tmp - val(k) * z(si_col_idx(k))
      enddo
      !
      z(i) = omega * tmp / diag_i
      !
   enddo
   !
   end subroutine apply_ssor_precond
   !
   !
   subroutine backsubstitute_fluxes_si(dt)
   !
   ! After solving the pressure system, update fluxes and water levels.
   !
   ! q(ip) = q_star(ip) - si_coeff(ip) * (eta_new(nmu) - eta_new(nm))
   ! zs(nm) = si_x(row)  for interior cells
   !
   implicit none
   !
   !
   real*4, intent(in) :: dt
   !
   integer :: ip, nm, nmu_z, irow
   real*4  :: eta_nm, eta_nmu, hu
   !
   ! Back-substitute fluxes
   !
   !$omp parallel &
   !$omp private ( ip, nm, nmu_z, eta_nm, eta_nmu, hu )
   !$omp do schedule ( dynamic, 256 )
   do ip = 1, npuv
      !
      ! Strictly the wet flag. The old test also accepted si_coeff(ip) > 0.0, which let DRY
      ! faces through: sfincs_momentum.f90:793-797 sets q = 0 and kfuv = 0 for a dry face but
      ! leaves si_q_star and si_coeff at whatever they held when the face was last wet, so
      ! backsubstitution overwrote that zero with a spurious flux built from stale
      ! coefficients. Every face that ever wetted then kept injecting water forever.
      !
      if (kfuv(ip) == 1) then
         !
         nm    = uv_index_z_nm(ip)
         nmu_z = uv_index_z_nmu(ip)
         !
         ! Get new water levels: from solver for interior, from zs for boundary
         !
         if (kcs(nm) == 1 .and. si_row_of_nm(nm) > 0) then
            eta_nm = si_x(si_row_of_nm(nm))
         else
            eta_nm = real(zs(nm))
         endif
         !
         if (kcs(nmu_z) == 1 .and. si_row_of_nm(nmu_z) > 0) then
            eta_nmu = si_x(si_row_of_nm(nmu_z))
         else
            eta_nmu = real(zs(nmu_z))
         endif
         !
         ! Update flux: q = q_star - si_coeff * (eta_new(nmu) - eta_new(nm))
         !
         q(ip) = si_q_star(ip) - si_coeff(ip) * (eta_nmu - eta_nm)
         !
         ! Stop water flowing out of a cell that has no water in it.
         !
         ! sfincs_momentum.f90:713 skips this when semi_implicit, on the understanding that
         ! backsubstitution applies it instead, so the criterion has to match what momentum
         ! would have used. With subgrid that is the negative-VOLUME test, not a level
         ! comparison: the subgrid table is inverted from z_volume, so letting a volume go
         ! negative sends the interpolation index to garbage and the run dies inside
         ! compute_water_levels_subgrid.
         !
         if (subgrid) then
            !
            if (z_volume(nm) < 0.0) then
               q(ip) = min(q(ip), 0.0)
            endif
            !
            if (z_volume(nmu_z) < 0.0) then
               q(ip) = max(q(ip), 0.0)
            endif
            !
         else
            !
            if (eta_nm < zb(nm)) then
               q(ip) = min(q(ip), 0.0)
            endif
            !
            if (eta_nmu < zb(nmu_z)) then
               q(ip) = max(q(ip), 0.0)
            endif
            !
         endif
         !
         ! Apply flux limiter
         !
         ! With subgrid there is no zbuvmx; the face reference level is subgrid_uv_zmin,
         ! which already has huthresh folded in when the tables are built (see the comment
         ! at sfincs_momentum.f90:177).
         !
         if (subgrid) then
            hu = max(real(max(zs(nm), zs(nmu_z))) - subgrid_uv_zmin(ip), huthresh)
         else
            hu = max(real(max(zs(nm), zs(nmu_z))) - zbuvmx(ip), huthresh)
         endif
         q(ip) = min(max(q(ip), -hu * uvlim), hu * uvlim)
         !
         ! Update velocity
         !
         uv(ip) = q(ip) / max(hu, huvmin)
         !
      endif
      !
   enddo
   !$omp end do
   !$omp end parallel
   !
   ! Update water levels from solver solution
   !
   !$omp parallel do private(irow, nm) schedule(static)
   do irow = 1, nrows_si
      !
      nm = si_nm_of_row(irow)
      zs(nm) = dble(si_x(irow))
      !
   enddo
   !$omp end parallel do
   !
   end subroutine backsubstitute_fluxes_si
   !
   !
   function get_tloop_si() result(t)
      real :: t
      t = tloop_si
   end function get_tloop_si
   !
   function get_si_iter_avg() result(avg)
      real :: avg
      if (si_solve_count > 0) then
         avg = real(si_iter_total) / real(si_solve_count)
      else
         avg = 0.0
      endif
   end function get_si_iter_avg
   !
   function get_si_iter_max() result(m)
      integer :: m
      m = si_iter_max_seen
   end function get_si_iter_max
   !
   function get_si_outer_avg() result(avg)
      real :: avg
      if (si_outer_total == 0) then
         avg = 0.0
      else
         avg = 1.0 * si_outer_total / max(si_solve_count_outer, 1)
      endif
   end function get_si_outer_avg
   !
   function get_si_outer_max() result(m)
      integer :: m
      m = si_outer_max_seen
   end function get_si_outer_max
   !
   function get_si_outer_capped() result(m)
      integer :: m
      m = si_outer_capped
   end function get_si_outer_capped
   !
   subroutine get_si_outer_profile(i, nsteps, mean_nbad, mean_dmax, nrows)
      ! Profile of outer iteration i: timesteps that reached it, mean rows still moving by
      ! more than si_tolouter, mean max change (m), and the row count for the fraction.
      integer, intent(in)  :: i
      integer, intent(out) :: nsteps, nrows
      real*4,  intent(out) :: mean_nbad, mean_dmax
      nsteps = si_prof_count(i)
      nrows  = nrows_tot
      if (nsteps > 0) then
         mean_nbad = real(si_prof_nbad(i)) / nsteps
         mean_dmax = real(si_prof_dmax(i) / nsteps)
      else
         mean_nbad = 0.0
         mean_dmax = 0.0
      endif
   end subroutine get_si_outer_profile
   !
end module sfincs_semi_implicit
