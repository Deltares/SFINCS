module bicgstab_solver
   use omp_lib
   implicit none
   contains

   ! Computes y = A * x for a sparse matrix in CSR format
   subroutine spmv(n, row_ptr, col_idx, values, vec, result)
      integer, intent(in) :: n, row_ptr(n+1), col_idx(*)
      real*4, intent(in) :: values(*), vec(n)
      real*4, intent(out) :: result(n)
      integer :: i, j

      !$omp parallel do private(j)
      do i = 1, n
         result(i) = 0.0
         do j = row_ptr(i), row_ptr(i+1) - 1
            result(i) = result(i) + values(j) * vec(col_idx(j))
         end do
      end do
      !$omp end parallel do
   end subroutine spmv

   ! BiCGStab solver for sparse matrices in CSR format
   ! Inputs:
   !   n       - Number of unknowns (size of the system)
   !   row_ptr - Row pointers for the CSR matrix (size n+1)
   !   col_idx - Column indices for the CSR matrix
   !   values  - Non-zero values of the CSR matrix
   !   x       - Initial guess (updated with solution)
   !   b       - Right-hand side vector
   !   tol     - Convergence tolerance
   !   max_iter- Maximum number of iterations
   subroutine bicgstab(n, row_ptr, col_idx, values, x, b, tol, max_iter)
      implicit none
      integer, intent(in) :: n, max_iter
      real*4, intent(in) :: tol
      integer, intent(in) :: row_ptr(n+1), col_idx(*)
      real*4, intent(in) :: values(*)
      real*4, intent(inout) :: x(n)
      real*4, intent(in) :: b(n)
      real*4 :: r(n), r0(n), p(n), v(n), s(n), t(n), x_new(n)
      real*4 :: rho, rho_old, alpha, omega, beta, norm_r, b_norm
      integer :: i, k

      ! Initial residual
      call spmv(n, row_ptr, col_idx, values, x, v)
      r = b - v
      r0 = r
      p = 0.0
      v = 0.0
      omega = 1.0
      rho_old = 1.0
      alpha = 1.0

      b_norm = sqrt(dot_product(b, b))
      if (b_norm < 1.0e-7) return

      do k = 1, max_iter
         rho = dot_product(r0, r)
         if (abs(rho) < 1.0e-7) exit
         beta = (rho / rho_old) * (alpha / omega)

         !$omp parallel do
         do i = 1, n
            p(i) = r(i) + beta * (p(i) - omega * v(i))
         end do
         !$omp end parallel do

         ! Matrix-vector product v = A*p
         call spmv(n, row_ptr, col_idx, values, p, v)
         if (abs(dot_product(r0, v)) < 1.0e-7) exit
         alpha = rho / dot_product(r0, v)

         ! Compute s
         !$omp parallel do
         do i = 1, n
            s(i) = r(i) - alpha * v(i)
         end do
         !$omp end parallel do

         norm_r = sqrt(dot_product(s, s))
         if (norm_r / b_norm < tol) then
            !$omp parallel do
            do i = 1, n
               x(i) = x(i) + alpha * p(i)
            end do
            !$omp end parallel do
            exit
         end if

         ! Matrix-vector product t = A*s
         call spmv(n, row_ptr, col_idx, values, s, t)
         if (abs(dot_product(t, t)) < 1.0e-7) exit
         omega = dot_product(t, s) / (dot_product(t, t) + 1.0e-7)

         ! Update solution
         !$omp parallel do
         do i = 1, n
            x_new(i) = x(i) + alpha * p(i) + omega * s(i)
         end do
         !$omp end parallel do

         ! Check convergence
         norm_r = sqrt(dot_product(x_new - x, x_new - x)) / b_norm
         x = x_new
         if (norm_r < tol) exit

         !$omp parallel do
         do i = 1, n
            r(i) = s(i) - omega * t(i)
         end do
         !$omp end parallel do

         rho_old = rho
      end do

   end subroutine bicgstab

end module bicgstab_solver