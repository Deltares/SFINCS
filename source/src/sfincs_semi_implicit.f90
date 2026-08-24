module sfincs_semi_implicit
   !
   ! Semi-implicit pressure treatment for SFINCS (Casulli-style theta-method)
   !
   ! Phase 1: Regular grid only (no subgrid, no quadtree)
   ! Solves Helmholtz equation for free surface elevation:
   !   eta^{n+1} - theta^2 * g * dt^2 * div(H * grad(eta^{n+1})) / A = RHS
   !
   ! Uses Conjugate Gradient solver with SSOR preconditioning.
   ! The system is SPD so CG is the natural choice.
   !
   use sfincs_data
   !
   implicit none
   !
   private
   public :: initialize_semi_implicit, assemble_and_solve_pressure, backsubstitute_fluxes_si
   public :: get_tloop_si, get_si_iter_avg, get_si_iter_max
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
   integer, dimension(:), allocatable :: si_row_face_side   ! 1=left 2=right 3=bottom 4=top
   !
   ! Position of the diagonal within each CSR row. This is
   ! what lets the preconditioner split lower from upper without a fixed stencil.
   !
   integer, dimension(:), allocatable :: si_diag_ptr     ! nrows_si
   !
   ! Timing and diagnostics
   !
   real    :: tloop_si
   integer :: si_iter_total       ! total CG iterations across all timesteps
   integer :: si_solve_count      ! number of solver calls
   integer :: si_iter_max_seen    ! max iterations in any single solve
   !
contains
   !
   subroutine initialize_semi_implicit()
   !
   ! Build the sparse matrix structure (CSR format) for the pressure Helmholtz system.
   ! The sparsity pattern is static (5-point stencil on regular grid).
   ! Only the values change each timestep.
   !
   ! Phase 1: regular grid only (subgrid=.false., use_quadtree=.false.)
   !
   implicit none
   !
   integer :: nm, ip, irow, icol, k
   integer :: nmu_z, nmd_z, num_z, ndm_z
   integer :: inb
   integer :: idir, nfaces_si, maxrow, ndiag
   !
   integer, dimension(:), allocatable :: row_count
   integer, dimension(:,:), allocatable :: si_nm_index  ! neighbor row indices (4, nrows_si)
   integer, dimension(:,:), allocatable :: slot_of_dir ! direction -> slot in si_AA, build only
   !
   ! Phase 1 guard: only regular grid without subgrid
   !
   if (subgrid) then
      write(*,*) 'Error: semi_implicit is not yet supported with subgrid=.true. (Phase 2)'
      stop
   endif
   !
   if (use_quadtree) then
      write(*,*) 'Error: semi_implicit is not yet supported with use_quadtree=.true. (Phase 3)'
      stop
   endif
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
   ! Allocate solver arrays
   !
   allocate(si_nm_of_row(nrows_si))
   allocate(si_row_ptr(nrows_si + 1))
   allocate(si_uv_index(4, nrows_si))
   allocate(si_nm_index(4, nrows_si))
   allocate(slot_of_dir(4, nrows_si))
   allocate(si_rhs(nrows_si))
   allocate(si_x(nrows_si))
   allocate(si_q_star(npuv))
   allocate(si_coeff(npuv))
   !
   ! CG work arrays (allocated once, reused every timestep)
   !
   allocate(cg_r(nrows_si))
   allocate(cg_z(nrows_si))
   allocate(cg_p(nrows_si))
   allocate(cg_Ap(nrows_si))
   allocate(cg_diag(nrows_si))
   !
   si_nm_of_row = 0
   si_row_ptr = 0
   si_uv_index = 0
   si_nm_index = 0
   slot_of_dir = 0
   si_rhs = 0.0
   si_x = 0.0
   si_q_star = 0.0
   si_coeff = 0.0
   cg_r = 0.0
   cg_z = 0.0
   cg_p = 0.0
   cg_Ap = 0.0
   cg_diag = 0.0
   !
   tloop_si = 0.0
   si_iter_total = 0
   si_solve_count = 0
   si_iter_max_seen = 0
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
   ! Find neighboring UV points and neighbor row indices for each row
   ! Uses the same stencil convention as sfincs_nonhydrostatic.f90:
   !
   !              4 (top/nu)
   !       +------|------+
   !       |             |
   !       |       irow  |
   !    1  -      + 5    -  2
   !  (md) |       (c)   | (mu)
   !       |             |
   !       +------|------+
   !              3 (bottom/nd)
   !
   do ip = 1, npuv
      !
      nm  = uv_index_z_nm(ip)
      nmu_z = uv_index_z_nmu(ip)
      !
      ! Only process UV points that touch at least one interior cell
      !
      if (kcs(nm) == 1 .or. kcs(nmu_z) == 1) then
         !
         if (uv_flags_dir(ip) == 0) then
            !
            ! x-direction UV point
            !
            if (kcs(nm) == 1) then
               irow = si_row_of_nm(nm)
               inb  = si_row_of_nm(nmu_z)  ! 0 if nmu is boundary
               si_uv_index(2, irow) = ip   ! right UV
               si_nm_index(2, irow) = inb  ! right neighbor row (0 if boundary)
            endif
            !
            if (kcs(nmu_z) == 1) then
               irow = si_row_of_nm(nmu_z)
               inb  = si_row_of_nm(nm)     ! 0 if nm is boundary
               si_uv_index(1, irow) = ip   ! left UV
               si_nm_index(1, irow) = inb  ! left neighbor row (0 if boundary)
            endif
            !
         else
            !
            ! y-direction UV point
            !
            if (kcs(nm) == 1) then
               irow = si_row_of_nm(nm)
               inb  = si_row_of_nm(nmu_z)
               si_uv_index(4, irow) = ip   ! top UV
               si_nm_index(4, irow) = inb
            endif
            !
            if (kcs(nmu_z) == 1) then
               irow = si_row_of_nm(nmu_z)
               inb  = si_row_of_nm(nm)
               si_uv_index(3, irow) = ip   ! bottom UV
               si_nm_index(3, irow) = inb
            endif
            !
         endif
         !
      endif
      !
   enddo
   !
   ! Build CSR sparsity pattern
   ! Order: left(1), bottom(3), center(5), top(4), right(2)
   ! Same ordering as sfincs_nonhydrostatic.f90
   !
   ! Pass 1: count the entries in each row, then prefix-sum into si_row_ptr.
   ! Rows are variable length. Nothing below assumes five entries, which is what makes
   ! room for a quadtree row (up to eight neighbours plus the diagonal) later.
   !
   allocate(row_count(nrows_si))
   !
   do irow = 1, nrows_si
      row_count(irow) = 1                     ! the diagonal is always present
      do idir = 1, 4
         if (si_nm_index(idir, irow) > 0) row_count(irow) = row_count(irow) + 1
      enddo
   enddo
   !
   si_row_ptr(1) = 1
   do irow = 1, nrows_si
      si_row_ptr(irow + 1) = si_row_ptr(irow) + row_count(irow)
   enddo
   !
   nnz_si = si_row_ptr(nrows_si + 1) - 1
   !
   deallocate(row_count)
   !
   allocate(si_col_idx(nnz_si))
   allocate(si_diag_ptr(nrows_si))
   !
   si_col_idx = 0
   si_diag_ptr = 0
   !
   ! Pass 2: fill each row.
   ! Entries are emitted in the historical order left(1), bottom(3), center(5), top(4),
   ! right(2). Columns are therefore NOT sorted ascending, and they do not need to be:
   ! the preconditioner splits lower from upper by comparing si_col_idx against the row
   ! index, which works for any ordering. Keeping the original order keeps the CSR layout
   ! byte-for-byte what it was, so the matrix-vector product sums in the same sequence and
   ! the refactor stays bit-identical.
   !
   do irow = 1, nrows_si
      !
      k = si_row_ptr(irow) - 1
      !
      icol = si_nm_index(1, irow)             ! left
      if (icol > 0) then
         k = k + 1
         si_col_idx(k) = icol
         slot_of_dir(1, irow) = k
      endif
      !
      icol = si_nm_index(3, irow)             ! bottom
      if (icol > 0) then
         k = k + 1
         si_col_idx(k) = icol
         slot_of_dir(3, irow) = k
      endif
      !
      k = k + 1                               ! centre
      si_col_idx(k) = irow
      si_diag_ptr(irow) = k
      !
      icol = si_nm_index(4, irow)             ! top
      if (icol > 0) then
         k = k + 1
         si_col_idx(k) = icol
         slot_of_dir(4, irow) = k
      endif
      !
      icol = si_nm_index(2, irow)             ! right
      if (icol > 0) then
         k = k + 1
         si_col_idx(k) = icol
         slot_of_dir(2, irow) = k
      endif
      !
   enddo
   !
   allocate(si_AA(nnz_si))
   si_AA = 0.0
   !
   ! Build the per-row face list.
   ! Faces are listed in the original direction order (left, right, bottom, top) so the
   ! diagonal accumulates in exactly the same sequence as the previous hard-coded version
   ! and the assembly stays bit-for-bit identical on a regular grid.
   ! When quadtree support lands this loop is the only place that needs to change: it
   ! walks whatever faces a row has, and nothing downstream assumes there are four.
   !
   allocate(si_row_face_ptr(nrows_si + 1))
   si_row_face_ptr = 0
   !
   k = 0
   do irow = 1, nrows_si
      do idir = 1, 4
         if (si_uv_index(idir, irow) > 0) k = k + 1
      enddo
   enddo
   nfaces_si = k
   !
   allocate(si_row_face_ip(nfaces_si))
   allocate(si_row_face_slot(nfaces_si))
   allocate(si_row_face_side(nfaces_si))
   !
   k = 0
   do irow = 1, nrows_si
      si_row_face_ptr(irow) = k + 1
      do idir = 1, 4
         ip = si_uv_index(idir, irow)
         if (ip > 0) then
            k = k + 1
            si_row_face_ip(k)   = ip
            si_row_face_slot(k) = slot_of_dir(idir, irow)   ! 0 means Dirichlet neighbour
            si_row_face_side(k) = idir
         endif
      enddo
   enddo
   si_row_face_ptr(nrows_si + 1) = k + 1
   !
   ! Report the sparse structure. The maximum row length is the useful number: 5 means a
   ! plain regular grid, and anything above 5 means quadtree connectivity is being picked
   ! up (a coarse cell next to refinement reaches up to 9). If that stays at 5 on a
   ! quadtree model, the face list is not seeing the refinement transitions.
   !
   ! Check the CSR invariants the preconditioner relies on. It splits lower from upper by
   ! comparing the column index against the row index, so every row must contain exactly
   ! one entry on the diagonal, si_diag_ptr must point at it, and every other entry must
   ! be strictly off-diagonal. If the quadtree row builder ever emits a duplicate column
   ! or misses the diagonal, SSOR would silently stop being a valid preconditioner.
   !
   maxrow = 0
   do irow = 1, nrows_si
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
   enddo
   !
   write(*,'(a,i10,a,i10,a,i3,a,i10)') ' Semi-implicit: rows ', nrows_si, &
      '  nonzeros ', nnz_si, '  max row ', maxrow, '  faces ', nfaces_si
   !
   deallocate(si_nm_index)
   deallocate(slot_of_dir)
   !
   end subroutine initialize_semi_implicit
   !
   !
   subroutine assemble_and_solve_pressure(dt)
   !
   ! Assemble the Helmholtz pressure matrix and RHS, then solve.
   !
   ! For regular grid (Phase 1), the system is:
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
   integer :: kface, iside, islot
   real*4  :: coeff_face
   real*4  :: diag
   real*4  :: dxr_val, dyr_val
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
      si_x(irow) = real(zs(si_nm_of_row(irow)))
   enddo
   !$omp end parallel do
   !
   ! Assemble matrix and RHS row by row
   !
   !$omp parallel do private(irow, nm, nmd, nmu, ndm, num, dxr_val, dyr_val, &
   !$omp                      div_qstar, diag, ip, kface, iside, islot, coeff_face) &
   !$omp schedule(static)
   do irow = 1, nrows_si
      !
      nm = si_nm_of_row(irow)
      !
      ! Get UV indices for flux divergence
      !
      nmd = z_index_uv_md(nm)  ! left UV
      nmu = z_index_uv_mu(nm)  ! right UV
      ndm = z_index_uv_nd(nm)  ! bottom UV
      num = z_index_uv_nu(nm)  ! top UV
      !
      ! Grid spacing (for regular grid, use reference level)
      !
      if (crsgeo) then
         dxr_val = dxm(nm)
         dyr_val = 1.0 / dyrinv(z_flags_iref(nm))
      else
         dxr_val = 1.0 / dxrinv(z_flags_iref(nm))
         dyr_val = 1.0 / dyrinv(z_flags_iref(nm))
      endif
      !
      ! Compute flux divergence of q_star for RHS
      ! Same formula as sfincs_continuity.f90 compute_water_levels_regular
      !
      if (crsgeo) then
         div_qstar = (si_q_star(nmd) - si_q_star(nmu)) / dxm(nm) &
                   + (si_q_star(ndm) - si_q_star(num)) * dyrinv(z_flags_iref(nm))
      else
         div_qstar = (si_q_star(nmd) - si_q_star(nmu)) * dxrinv(z_flags_iref(nm)) &
                   + (si_q_star(ndm) - si_q_star(num)) * dyrinv(z_flags_iref(nm))
      endif
      !
      ! RHS = current zs + dt * div(q_star) + dt * sources
      ! (div_qstar already has the sign convention: inflow positive)
      !
      si_rhs(irow) = real(zs(nm)) + dt * div_qstar
      !
      ! Include precipitation in the pressure system so that the solver
      ! accounts for the added volume when computing fluxes
      !
      if (precip) then
         si_rhs(irow) = si_rhs(irow) + dt * netprcp(nm)
      endif
      !
      ! Include external sources (e.g. from BMI/XMI coupling)
      !
      if (use_qext) then
         si_rhs(irow) = si_rhs(irow) + dt * qext(nm)
      endif
      !
      ! Now assemble matrix coefficients
      ! The coefficient comes from substituting the momentum into continuity.
      ! For regular grid: cell_area = dxr * dyr
      !
      diag = 1.0
      !
      ! Walk this row's faces. Grid-agnostic: nothing here assumes there are four.
      ! Sides 1 and 2 are x-direction, 3 and 4 are y-direction.
      ! For a Dirichlet neighbour (kcs==2) the coefficient stays on the diagonal and the
      ! known water level moves to the RHS. Which of the two cells on the face is the
      ! boundary depends on the side: for left/bottom it is uv_index_z_nm, for right/top
      ! it is uv_index_z_nmu.
      !
      do kface = si_row_face_ptr(irow), si_row_face_ptr(irow + 1) - 1
         !
         ip    = si_row_face_ip(kface)
         iside = si_row_face_side(kface)
         islot = si_row_face_slot(kface)
         !
         if (iside <= 2) then
            coeff_face = si_coeff(ip) * dt / dxr_val
         else
            coeff_face = si_coeff(ip) * dt / dyr_val
         endif
         !
         diag = diag + coeff_face
         !
         if (islot > 0) then
            ! Interior neighbour
            si_AA(islot) = -coeff_face
         else
            ! Boundary neighbour (kcs==2): known eta, move to RHS
            if (iside == 1 .or. iside == 3) then
               si_rhs(irow) = si_rhs(irow) + coeff_face * real(zs(uv_index_z_nm(ip)))
            else
               si_rhs(irow) = si_rhs(irow) + coeff_face * real(zs(uv_index_z_nmu(ip)))
            endif
         endif
         !
      enddo
      !
      ! Set diagonal
      !
      si_AA(si_diag_ptr(irow)) = diag
      !
   enddo
   !$omp end parallel do
   !
   ! Solve using CG with Jacobi preconditioning
   !
   call cg_solve(nrows_si, nnz_si, si_AA, si_col_idx, si_row_ptr, &
                  si_rhs, si_x, si_tol, si_maxiter, iter, relres)
   !
   si_iter_total = si_iter_total + iter
   si_solve_count = si_solve_count + 1
   si_iter_max_seen = max(si_iter_max_seen, iter)
   !
   call system_clock(count1, count_rate, count_max)
   tloop_si = tloop_si + 1.0 * (count1 - count0) / count_rate
   !
   end subroutine assemble_and_solve_pressure
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
   ! Works on the CSR structure directly. Lower and upper are separated by comparing the
   ! column index against the row index, so this makes no assumption about how many
   ! entries a row has or what order they are stored in — which is what a quadtree row
   ! needs. The diagonal is located via si_diag_ptr rather than a fixed stencil slot.
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
   integer :: i, k, icol
   real*4  :: diag_i, tmp
   !
   ! Forward sweep: (D/omega + L) * z = r
   ! Ascending k, taking entries whose column lies below the diagonal.
   !
   do i = 1, n
      !
      diag_i = val(si_diag_ptr(i))
      tmp = r(i)
      !
      do k = si_row_ptr(i), si_row_ptr(i + 1) - 1
         icol = si_col_idx(k)
         if (icol < i) tmp = tmp - val(k) * z(icol)
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
   ! Descending k, taking entries whose column lies above the diagonal.
   ! Descending rather than ascending on purpose: with the historical emission order
   ! (left, bottom, centre, top, right) it visits right before top, which is the order
   ! the previous slot-indexed version used, so the sum rounds identically.
   !
   do i = n, 1, -1
      !
      diag_i = val(si_diag_ptr(i))
      tmp = z(i)
      !
      do k = si_row_ptr(i + 1) - 1, si_row_ptr(i), -1
         icol = si_col_idx(k)
         if (icol > i) tmp = tmp - val(k) * z(icol)
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
      if (kfuv(ip) == 1 .or. si_coeff(ip) > 0.0) then
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
         ! Apply dry-cell limiter
         !
         if (eta_nm < zb(nm)) then
            q(ip) = min(q(ip), 0.0)
         endif
         !
         if (eta_nmu < zb(nmu_z)) then
            q(ip) = max(q(ip), 0.0)
         endif
         !
         ! Apply flux limiter
         !
         hu = max(real(max(zs(nm), zs(nmu_z))) - zbuvmx(ip), huthresh)
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
end module sfincs_semi_implicit
