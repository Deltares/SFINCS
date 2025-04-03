module bicgstab_solver
   
   use omp_lib
   implicit none

contains

   ! Computes y = A * x for a sparse matrix in CSR format
   subroutine spmv(n, row_ptr, col_idx, values, vec, result)
      integer, intent(in) :: n, row_ptr(n+1), col_idx(*)
      real*4, intent(in) :: values(*)
      real*4, intent(in) :: vec(n)
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
   !   n        - Number of unknowns (size of the system)
   !   row_ptr  - Row pointers for the CSR matrix (size n+1)
   !   col_idx  - Column indices for the CSR matrix
   !   values   - Non-zero values of the CSR matrix
   !   x        - Initial guess (updated with solution)
   !   b        - Right-hand side vector
   !   tol      - Convergence tolerance
   !   max_iter - Maximum number of iterations
   !   nnz      - number of values in matrix
   !   use_preconditioner - (0/1)
   subroutine bicgstab_mine(n, row_ptr, col_idx, values, x, b, tol, max_iter, nnz, use_preconditioner)
      implicit none
      integer, intent(in) :: n, max_iter, nnz, use_preconditioner
      real*4, intent(in) :: tol
      integer, intent(in) :: row_ptr(n+1), col_idx(*)
      real*4, intent(in) :: values(*)
      real*4, intent(inout) :: x(n)
      real*4, intent(in) :: b(n)
      real*4 :: r(n), r0(n), p(n), v(n), s(n), t(n)
      real*4 :: rho, rho_old, alpha, omega, beta, norm_r, b_norm
      integer :: i, k

      ! Initial residual
      call spmv(n, row_ptr, col_idx, values, x, v)
      r = b - v
      r0 = r

      ! Generate a random vector r0 (uniform random numbers in the range [0, 1])
      !call random_seed()  ! Initialize the random number generator
      !call random_number(r0)
      !r0 = r0 * sqrt(dot_product(r, r))
      
      rho = dot_product(r0, r)
      p = r0
      rho_old = 1.0      

      b_norm = sqrt(dot_product(b, b))
      if (b_norm < 1.0e-7) return

      do k = 1, max_iter

         !print *, 'Iteration:', k, 'rho:', rho

         ! 1)
         call spmv(n, row_ptr, col_idx, values, p, v)
         
         
         if (abs(dot_product(r0, v)) < 1.0e-7) exit
         
         ! 2)
         alpha = rho / (dot_product(r0, v))

         ! Skipping 3) 
         ! h = x + alpha * p
         
         ! 4)
         !$omp parallel do
         do i = 1, n
            s(i) = r(i) - alpha * v(i)
         end do
         !$omp end parallel do
         
         ! 5) If h is small enough, i.e. if s is small enough, then set x = h and quit
         !write(*,*)'dps',sqrt(dot_product(s, s))
         
         ! 6)
         call spmv(n, row_ptr, col_idx, values, s, t) 

         ! 7)
         omega = dot_product(t, s) / (dot_product(t, t) + 1.0e-7)

         ! 8) do not use h from step 3 but  x + alpha * p directly
         !$omp parallel do
         do i = 1, n
            x(i) = x(i) + alpha * p(i) + omega * s(i)
         end do
         !$omp end parallel do

         ! 9)
         !$omp parallel do
         do i = 1, n
            r(i) = s(i) - omega * t(i)
         end do
         !$omp end parallel do
         
         ! 10) If x is accurate enough, i.e. if r is small enough, then quit
         ! Check convergence
         norm_r = sqrt(dot_product(r, r))
         if (norm_r < tol) exit

         ! 11)
         rho = dot_product(r0, r)
         
         ! 12)      
         beta = (rho / rho_old) * (alpha / omega)

         ! 13)
         !$omp parallel do
         do i = 1, n
            p(i) = r(i) + beta * (p(i) - omega * v(i))
         end do
         !$omp end parallel do

         rho_old = rho

      end do

   end subroutine bicgstab_mine
      
end module bicgstab_solver
