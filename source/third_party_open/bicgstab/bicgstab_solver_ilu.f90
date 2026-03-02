module bicgstab_solver_ilu
   !
   use omp_lib
   implicit none
   private
   public :: bicgstab_solve
   !
   integer, parameter :: sp = kind(1.0)
   !
contains
   !
   subroutine bicgstab_solve(n, val, col_ind, row_ptr, b, x, tol, maxiter, iter, relres, use_precond)
      !
      integer, intent(in) :: n, col_ind(:), row_ptr(:), maxiter
      real(sp), intent(in) :: val(:), b(:), tol
      real(sp), intent(inout) :: x(:)
      integer, intent(out) :: iter
      real(sp), intent(out) :: relres
      logical, intent(in) :: use_precond
      !
      real(sp), allocatable :: r(:), r0(:), p(:), v(:), s(:), t(:)
      real(sp), allocatable :: phat(:), shat(:)
      real(sp), allocatable :: LU(:)
      integer, allocatable :: diag_ptr(:)
      real(sp) :: rho, alpha, omega, beta, rs_old
      integer :: i
      real(sp) :: bnorm, rnorm
      !
      allocate(r(n), r0(n), p(n), v(n), s(n), t(n))
      allocate(phat(n), shat(n))
      !
      if (use_precond) then
         !
         allocate(LU(size(val)), diag_ptr(n))
         call ilu0_factor(n, val, col_ind, row_ptr, LU, diag_ptr)
         !
      endif
      !
      call spmv(n, val, col_ind, row_ptr, x, r)
      r = b - r
      r0 = r
      p = 0.0_sp
      v = 0.0_sp
      rho = 1.0_sp
      alpha = 1.0_sp
      omega = 1.0_sp
      !
      bnorm = max(norm2(b), 1.0e-30_sp)
      rnorm = norm2(r)
      relres = rnorm / bnorm
      !
      iter = 0
      !
      do while (relres > tol .and. iter < maxiter)
         !
         iter = iter + 1
         rs_old = dot_product(r0, r)
         !
         if (rs_old == 0.0_sp) exit
         !
         if (iter == 1) then
            p = r
         else
            beta = (rs_old / rho) * (alpha / omega)
            p = r + beta * (p - omega * v)
         endif
         !
         if (use_precond) then
            call apply_precond(n, p, phat, LU, col_ind, row_ptr, diag_ptr)
         else
            phat = p
         endif
         !
         call spmv(n, val, col_ind, row_ptr, phat, v)
         alpha = rs_old / dot_product(r0, v)
         s = r - alpha * v
         !
         if (use_precond) then
            call apply_precond(n, s, shat, LU, col_ind, row_ptr, diag_ptr)
         else
            shat = s
         endif
         !
         call spmv(n, val, col_ind, row_ptr, shat, t)
         omega = dot_product(t, s) / dot_product(t, t)
         !
         x = x + alpha * phat + omega * shat
         r = s - omega * t
         !
         rnorm = norm2(r)
         relres = rnorm / bnorm
         !
         if (relres < tol .or. omega == 0.0_sp) exit
         !
         rho = rs_old
         !
      enddo
      !
      deallocate(r, r0, p, v, s, t, phat, shat)
      !
      if (use_precond) deallocate(LU, diag_ptr)
      !
   end subroutine bicgstab_solve

   subroutine ilu0_factor(n, val, col_ind, row_ptr, LU, diag_ptr)
      !
      integer, intent(in) :: n, col_ind(:), row_ptr(:)
      real(sp), intent(in) :: val(:)
      real(sp), intent(out) :: LU(:)
      integer, intent(out) :: diag_ptr(:)
      integer :: i, j, k, m, row_start, row_end, col_k
      real(sp) :: val_k
      !
      LU = val
      diag_ptr = 0
      !
      do i = 1, n
         do k = row_ptr(i), row_ptr(i + 1) - 1
            !
            if (col_ind(k) == i) then
               diag_ptr(i) = k
               exit
            endif
            !
         enddo
      enddo
      !
      do i = 2, n
         !
         row_start = row_ptr(i)
         !
         row_end = row_ptr(i + 1) - 1
         !
         do k = row_start, row_end
            !
            col_k = col_ind(k)
            !
            if (col_k >= i) cycle
            !
            val_k = LU(k) / LU(diag_ptr(col_k))
            LU(k) = val_k
            !
            do j = row_ptr(col_k), row_ptr(col_k + 1) - 1
               !
               if (col_ind(j) <= col_k) cycle
               !
               do m = row_start, row_end
                  !
                  if (col_ind(m) == col_ind(j)) then
                     LU(m) = LU(m) - val_k * LU(j)
                     exit
                  endif
               enddo
            enddo
         enddo
      enddo
      !
   end subroutine ilu0_factor

   subroutine apply_precond(n, x, y, LU, col_ind, row_ptr, diag_ptr)
      integer, intent(in) :: n, col_ind(:), row_ptr(:), diag_ptr(:)
      real(sp), intent(in) :: x(:), LU(:)
      real(sp), intent(out) :: y(:)
      real(sp), allocatable :: z(:)
      integer :: i, k
      real(sp) :: sum
      !
      allocate(z(n))
      !
      ! Forward solve: L y = x
      !
      do i = 1, n
         !
         sum = x(i)
         !
         do k = row_ptr(i), diag_ptr(i) - 1
            sum = sum - LU(k) * y(col_ind(k))
         enddo
         !
         y(i) = sum
         !
      enddo
      !
      ! Backward solve: U z = y
      !
      do i = n, 1, -1
         !
         sum = y(i)
         !
         do k = diag_ptr(i) + 1, row_ptr(i + 1) - 1
            !
            sum = sum - LU(k) * z(col_ind(k))
            !
         enddo
         !
         z(i) = sum / LU(diag_ptr(i))
         !
      enddo
      !
      y = z
      deallocate(z)
      !
   endsubroutine apply_precond

   subroutine spmv(n, val, col_ind, row_ptr, x, y)
      !
      integer, intent(in) :: n, col_ind(:), row_ptr(:)
      real(sp), intent(in) :: val(:), x(:)
      real(sp), intent(out) :: y(:)
      integer :: i, k
      !
      y = 0.0_sp
      !
      !$omp parallel &
      !$omp private ( i, k )
      !$omp do schedule ( dynamic, 256 )
      do i = 1, n
         do k = row_ptr(i), row_ptr(i + 1) - 1
            y(i) = y(i) + val(k) * x(col_ind(k))
         enddo
      enddo
      !$omp end do
      !$omp end parallel
      !
   end subroutine spmv

   function norm2(v) result (snorm2)
      ! 
      real(sp), intent(in) :: v(:)
      integer :: i
      real(sp) :: local_sum
      real(sp) :: snorm2
      !
      local_sum = 0.0_sp
      !
      !$omp parallel &
      !$omp private ( i ) &
      !$omp reduction(+:local_sum)
      !$omp do
      do i = 1, size(v)
         local_sum = local_sum + v(i) * v(i)
      enddo
      !$omp end do
      !$omp end parallel
      !
      snorm2 = sqrt(local_sum)
      !
  end function norm2

end module bicgstab_solver_ilu
