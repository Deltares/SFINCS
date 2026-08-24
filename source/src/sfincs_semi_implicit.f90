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
   integer, dimension(:), allocatable :: si_row_face_isy    ! 0 = x-direction, 1 = y
   integer, dimension(:), allocatable :: si_row_face_bnd    ! cell to read zs from when Dirichlet
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
   ! The sparsity pattern is static; only the values change each timestep.
   ! Rows are variable length: 5 entries on a regular grid, up to 9 next to a quadtree
   ! refinement transition.
   !
   implicit none
   !
   integer :: nm, ip, irow, icol, k
   integer :: j, nb, nfaces_si, maxrow, ndiag
   !
   integer, dimension(:), allocatable :: row_count
   integer, dimension(:,:), allocatable :: face_ip    ! (8, nrows_si) UV point per face slot
   integer, dimension(:,:), allocatable :: face_nb    ! (8, nrows_si) neighbour row, 0 = Dirichlet
   integer, dimension(:,:), allocatable :: face_bnd   ! (8, nrows_si) neighbour cell index
   integer, dimension(:,:), allocatable :: face_isy   ! (8, nrows_si) 0 = x, 1 = y
   integer, dimension(:,:), allocatable :: face_slot  ! (8, nrows_si) slot in si_AA
   !
   ! Subgrid needs a nonlinear outer iteration that does not exist yet.
   !
   if (subgrid) then
      write(*,*) 'Error: semi_implicit is not yet supported with subgrid=.true. (Phase 2)'
      stop
   endif
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
   ! Allocate solver arrays
   !
   allocate(si_nm_of_row(nrows_si))
   allocate(si_row_ptr(nrows_si + 1))
   allocate(face_ip(8, nrows_si))
   allocate(face_nb(8, nrows_si))
   allocate(face_bnd(8, nrows_si))
   allocate(face_isy(8, nrows_si))
   allocate(face_slot(8, nrows_si))
   allocate(si_rhs(nrows_si))
   allocate(si_x(nrows_si))
   ! Sized exactly like q and uv (sfincs_domain.f90:2196), NOT npuv.
   !
   ! A quadtree creates ncuv combined uv points that live past npuv, and div_qstar below
   ! reads z_index_uv_md/mu/nd/nu, which point at those combined points next to a
   ! refinement transition. The +1 is the sentinel slot sfincs_domain.f90:1251 assigns to
   ! unset indices. Allocating only npuv read past the end -- silently, in Release.
   !
   allocate(si_q_star(npuv + ncuv + 1))
   allocate(si_coeff(npuv + ncuv + 1))
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
   allocate(row_count(nrows_si))
   !
   do irow = 1, nrows_si
      row_count(irow) = 1                       ! the diagonal is always present
      do j = 1, 8
         if (face_nb(j, irow) > 0) row_count(irow) = row_count(irow) + 1
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
   ! Columns are NOT sorted ascending. They do not need to be: the preconditioner separates
   ! lower from upper by comparing the column index against the row index.
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
   enddo
   !
   allocate(si_AA(nnz_si))
   si_AA = 0.0
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
   ! The maximum row length is the useful number: 5 on a regular grid, above 5 once
   ! quadtree connectivity is picked up (a coarse cell next to refinement reaches 9).
   ! If it stays at 5 on a quadtree model, the refinement transitions are not being seen.
   !
   write(*,'(a,i10,a,i10,a,i3,a,i10)') ' Semi-implicit: rows ', nrows_si, &
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
   integer :: kface, islot
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
   !$omp                      div_qstar, diag, ip, kface, islot, coeff_face) &
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
         islot = si_row_face_slot(kface)
         !
         if (si_row_face_isy(kface) == 0) then
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
            ! Boundary neighbour (kcs==2): known eta, move to RHS.
            ! si_row_face_bnd already holds whichever end of the face is not this cell,
            ! so there is no left/right special case to get wrong.
            si_rhs(irow) = si_rhs(irow) + coeff_face * real(zs(si_row_face_bnd(kface)))
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
