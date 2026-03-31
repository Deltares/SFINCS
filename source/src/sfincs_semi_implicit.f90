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
   !
   integer, dimension(:), allocatable :: col_idx0
   integer, dimension(:,:), allocatable :: si_nm_index  ! neighbor row indices (4, nrows_si)
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
   allocate(si_index_sparse(5, nrows_si))
   allocate(si_row_ptr(nrows_si + 1))
   allocate(col_idx0(5 * nrows_si))
   allocate(si_uv_index(4, nrows_si))
   allocate(si_nm_index(4, nrows_si))
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
   si_index_sparse = 0
   si_row_ptr = 0
   col_idx0 = 0
   si_uv_index = 0
   si_nm_index = 0
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
   k = 0
   !
   do irow = 1, nrows_si
      !
      ! Left neighbor
      !
      icol = si_nm_index(1, irow)
      if (icol > 0) then
         k = k + 1
         col_idx0(k) = icol
         si_index_sparse(1, irow) = k
         if (si_row_ptr(irow) == 0) si_row_ptr(irow) = k
      endif
      !
      ! Bottom neighbor
      !
      icol = si_nm_index(3, irow)
      if (icol > 0) then
         k = k + 1
         col_idx0(k) = icol
         si_index_sparse(3, irow) = k
         if (si_row_ptr(irow) == 0) si_row_ptr(irow) = k
      endif
      !
      ! Center (always present)
      !
      k = k + 1
      col_idx0(k) = irow
      si_index_sparse(5, irow) = k
      if (si_row_ptr(irow) == 0) si_row_ptr(irow) = k
      !
      ! Top neighbor
      !
      icol = si_nm_index(4, irow)
      if (icol > 0) then
         k = k + 1
         col_idx0(k) = icol
         si_index_sparse(4, irow) = k
         if (si_row_ptr(irow) == 0) si_row_ptr(irow) = k
      endif
      !
      ! Right neighbor
      !
      icol = si_nm_index(2, irow)
      if (icol > 0) then
         k = k + 1
         col_idx0(k) = icol
         si_index_sparse(2, irow) = k
         if (si_row_ptr(irow) == 0) si_row_ptr(irow) = k
      endif
      !
   enddo
   !
   nnz_si = k
   si_row_ptr(nrows_si + 1) = k + 1
   !
   allocate(si_col_idx(nnz_si))
   si_col_idx(1:nnz_si) = col_idx0(1:nnz_si)
   !
   allocate(si_AA(nnz_si))
   si_AA = 0.0
   !
   deallocate(col_idx0)
   deallocate(si_nm_index)
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
   real*4  :: coeff_left, coeff_right, coeff_bottom, coeff_top
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
   !$omp                      div_qstar, diag, ip, coeff_left, coeff_right, coeff_bottom, coeff_top) &
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
      ! Left (direction 1): UV point at left face, x-direction
      !
      ip = si_uv_index(1, irow)
      if (ip > 0) then
         coeff_left = si_coeff(ip) * dt / dxr_val
         diag = diag + coeff_left
         !
         if (si_index_sparse(1, irow) > 0) then
            ! Interior neighbor
            si_AA(si_index_sparse(1, irow)) = -coeff_left
         else
            ! Boundary neighbor (kcs==2): known eta, move to RHS
            si_rhs(irow) = si_rhs(irow) + coeff_left * real(zs(uv_index_z_nm(ip)))
         endif
      endif
      !
      ! Right (direction 2): UV point at right face, x-direction
      !
      ip = si_uv_index(2, irow)
      if (ip > 0) then
         coeff_right = si_coeff(ip) * dt / dxr_val
         diag = diag + coeff_right
         !
         if (si_index_sparse(2, irow) > 0) then
            si_AA(si_index_sparse(2, irow)) = -coeff_right
         else
            ! Boundary: right UV, boundary cell is uv_index_z_nmu(ip)
            si_rhs(irow) = si_rhs(irow) + coeff_right * real(zs(uv_index_z_nmu(ip)))
         endif
      endif
      !
      ! Bottom (direction 3): UV point at bottom face, y-direction
      !
      ip = si_uv_index(3, irow)
      if (ip > 0) then
         coeff_bottom = si_coeff(ip) * dt / dyr_val
         diag = diag + coeff_bottom
         !
         if (si_index_sparse(3, irow) > 0) then
            si_AA(si_index_sparse(3, irow)) = -coeff_bottom
         else
            si_rhs(irow) = si_rhs(irow) + coeff_bottom * real(zs(uv_index_z_nm(ip)))
         endif
      endif
      !
      ! Top (direction 4): UV point at top face, y-direction
      !
      ip = si_uv_index(4, irow)
      if (ip > 0) then
         coeff_top = si_coeff(ip) * dt / dyr_val
         diag = diag + coeff_top
         !
         if (si_index_sparse(4, irow) > 0) then
            si_AA(si_index_sparse(4, irow)) = -coeff_top
         else
            si_rhs(irow) = si_rhs(irow) + coeff_top * real(zs(uv_index_z_nmu(ip)))
         endif
      endif
      !
      ! Set diagonal
      !
      si_AA(si_index_sparse(5, irow)) = diag
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
   ! (si_index_sparse) for efficient SSOR forward/backward sweeps.
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
      cg_diag(i) = val(si_index_sparse(5, i))
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
   ! Uses the known 5-point stencil structure stored in si_index_sparse:
   !   1 = left, 2 = right, 3 = bottom, 4 = top, 5 = center (diagonal)
   !
   ! SSOR = (D/omega + L) * D^{-1} * (D/omega + U)
   !
   ! Forward sweep:  (D/omega + L) * y = r
   ! Backward sweep: (D/omega + U) * z = D/omega * y
   !
   ! This is inherently sequential, but the stencil structure avoids
   ! CSR index searches. For the Helmholtz system, SSOR typically
   ! reduces CG iterations by 3-5x compared to Jacobi.
   !
   implicit none
   !
   integer, intent(in)  :: n
   real*4,  intent(in)  :: val(*), r(n), omega
   real*4,  intent(out) :: z(n)
   !
   integer :: i, idx
   real*4  :: diag_i, tmp, dom
   !
   ! Forward sweep: (D/omega + L) * z = r
   ! L = lower triangular part = left (1) and bottom (3) neighbors
   ! These neighbors have LOWER row indices (irow < i)
   !
   do i = 1, n
      !
      diag_i = val(si_index_sparse(5, i))
      dom = diag_i / omega
      tmp = r(i)
      !
      ! Left neighbor (index 1) — lower triangle
      !
      idx = si_index_sparse(1, i)
      if (idx > 0) then
         tmp = tmp - val(idx) * z(si_col_idx(idx))
      endif
      !
      ! Bottom neighbor (index 3) — lower triangle
      !
      idx = si_index_sparse(3, i)
      if (idx > 0) then
         tmp = tmp - val(idx) * z(si_col_idx(idx))
      endif
      !
      z(i) = omega * tmp / diag_i
      !
   enddo
   !
   ! Scale: z = D/omega * z (prepare for backward sweep)
   !
   !$omp parallel do private(i) schedule(static)
   do i = 1, n
      z(i) = val(si_index_sparse(5, i)) / omega * z(i)
   enddo
   !$omp end parallel do
   !
   ! Backward sweep: (D/omega + U) * z_new = z_old
   ! U = upper triangular part = right (2) and top (4) neighbors
   ! These neighbors have HIGHER row indices (irow > i)
   !
   do i = n, 1, -1
      !
      diag_i = val(si_index_sparse(5, i))
      dom = diag_i / omega
      tmp = z(i)
      !
      ! Right neighbor (index 2) — upper triangle
      !
      idx = si_index_sparse(2, i)
      if (idx > 0) then
         tmp = tmp - val(idx) * z(si_col_idx(idx))
      endif
      !
      ! Top neighbor (index 4) — upper triangle
      !
      idx = si_index_sparse(4, i)
      if (idx > 0) then
         tmp = tmp - val(idx) * z(si_col_idx(idx))
      endif
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
