module params
  implicit none
  integer,parameter :: dp = kind(1.0d0)
  integer,parameter :: size = 1024

  integer,parameter :: iter_max = 10000
  real(dp),parameter :: tol = 1.0d-10
end module params

module subs
  use params
  implicit none
contains
  subroutine init(size, a, x, b)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(out) :: a
    real(dp),dimension(size),intent(out) :: x, b
    integer :: i, j

    !$omp parallel
    !$omp workshare
    a = 0.0d0
    x = 0.0d0
    b = 1.0d0
    !$omp end workshare
    !$omp do
    do i = 1, size
       a(i,i) = 1.0d0*i
    end do
    !$omp end do
    !$omp end parallel

    a(size, 1)    = 2.0d0
    a(1,    size) = -2.0d0
    a(7,    2)    = 3.0d0
    a(2,    7)    = -3.0d0
    a(3,    6)    = 4.0d0
    a(6,    3)    = -4.0d0
    a(4,    5)    = 5.0d0
    a(5,    4)    = -5.0d0

    write(6, *) "size:", size
#ifdef _DEBUG
    write(6, *) "matrix A:"
    do i = 1, size
       write(6,'(8(1pe14.5))') (a(i, j), j = 1, size)
    end do

    write(6, *)
    write(6, *) "right hand side vector b:"
    do i = 1, size
       write(6, '(1pe14.5)') b(i)
    end do
#endif
  end subroutine init

  ! ax = A*x
  subroutine a_dot_x(size, a, x, ax)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(in) :: a
    real(dp),dimension(size),intent(in) :: x
    real(dp),dimension(size),intent(out) :: ax
    integer :: i, j

    !$omp parallel
    !$omp workshare
    ax = 0.0d0
    !$omp end workshare
    !$omp do
    do i = 1, size
       do j = 1, size
          ax(i) = ax(i) + a(i, j)*x(j) 
       end do
    end do
    !$omp end do
    !$omp end parallel
  end subroutine a_dot_x

  ! inner product xy = x*y
  subroutine x_dot_y(size, x, y, xy)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size),intent(in) :: x, y
    real(dp),intent(out) :: xy
    integer :: i
    
    xy = 0.0d0
    !$omp parallel do reduction(+:xy)
    do i = 1, size
       xy = xy + x(i)*y(i)
    end do
    
  end subroutine x_dot_y

  ! calc r = b - A*x
  subroutine b_minus_ax(size, a, x, b, r)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(in) :: a
    real(dp),dimension(size),intent(in) :: x, b
    real(dp),dimension(size),intent(out) :: r
    real(dp),dimension(size) :: ax
    integer :: i

    call a_dot_x(size, a, x, ax)
    !$omp parallel do
    do i = 1, size
       r(i) = b(i) - ax(i)
    end do
  end subroutine b_minus_ax

  ! need to solve asymmetric matrix
  ! http://www.jicfus.jp/wiki/index.php?Bi-CGSTAB%20%E6%B3%95
  ! https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
  subroutine bicgstab(size, a, x, b)
    implicit none
    integer,intent(in) :: size
    real(dp),dimension(size, size),intent(in) :: a
    real(dp),dimension(size),intent(out) :: x
    real(dp),dimension(size),intent(in) :: b
    integer :: i, j, iter
    real(dp),dimension(size) :: r, r_new, v, v_new, p, p_new
    real(dp),dimension(size) :: r0, h, s, t
    real(dp),dimension(size) :: rr
    real(dp),dimension(size) :: xx, xx_new
    real(dp) :: rho, rho_new, omega, omega_new
    real(dp) :: alpha, beta
    real(dp) :: r0r, t2, r0v, ts, res

    !$omp parallel
    !$omp workshare
    xx = b ! initial guess
    !$omp end workshare
    !$omp end parallel
    call b_minus_ax(size, a, xx, b, r)
    !$omp parallel
    !$omp workshare
    r0 = r
    !$omp end workshare
    !$omp end parallel
    call x_dot_y(size, r0, r, r0r)

    if (r0r == 0.0d0) then
       write(6, *) "r*r0 is zero."
       stop
    end if
    
    rho   = 1.0d0
    alpha = 1.0d0
    omega = 1.0d0

    !$omp parallel
    !$omp workshare
    v = 0.0d0
    p = 0.0d0
    !$omp end workshare
    !$omp end parallel

    do iter = 1, iter_max
       call x_dot_y(size, r0, r, rho_new)
       beta = (rho_new/rho)*(alpha/omega)
       !$omp parallel do
       do i = 1, size
          p_new(i) = r(i) + beta*(p(i)-omega*v(i))
       end do
       call a_dot_x(size, a, p_new, v_new)
       call x_dot_y(size, r0, v_new, r0v)
       alpha = rho_new/r0v
       !$omp parallel do
       do i = 1, size
          h(i) = xx(i) + alpha*p_new(i)
       end do
       call b_minus_ax(size, a, h, b, rr)
       call x_dot_y(size, rr, rr, res)
       res = sqrt(res)
       if (res <= tol) then
          !$omp parallel
          !$omp workshare
          xx_new = h
          !$omp end workshare
          !$omp end parallel
          exit
       end if
       !$omp parallel do
       do i = 1, size
          s(i) = r(i) - alpha*v_new(i)
       end do
       call a_dot_x(size, a, s, t)
       call x_dot_y(size, t, s, ts)
       call x_dot_y(size, t, t, t2)
       omega_new = ts/t2
       !$omp parallel do
       do i = 1, size
          xx_new(i) = h(i) + omega_new*s(i)
       end do
       call b_minus_ax(size, a, xx_new, b, rr)
       call x_dot_y(size, rr, rr, res)
       res = sqrt(res)
       if (res <= tol) exit
       !$omp parallel do
       do i = 1, size
          r_new(i) = s(i) - omega_new*t(i)
       end do

       !$omp parallel
       !$omp workshare
       xx    = xx_new
       p     = p_new
       r     = r_new
       v     = v_new
       !$omp end workshare
       !$omp end parallel
       rho   = rho_new
       omega = omega_new
    end do ! iter

    if (iter>=iter_max .and. res>tol) then
       write(6, *) "did not converge."
       write(6, *) "iter_max, res:", iter_max, res
       stop
    end if

    ! converged
    write(6, *) "BiCGSTAB method converged."
    write(6, *) "iter, res:", iter, res
    !$omp parallel
    !$omp workshare
    x = xx_new
    !$omp end workshare
    !$omp end parallel
  end subroutine bicgstab
end module subs

program main
  use subs
  implicit none
  real(dp),dimension(size, size) :: a
  real(dp),dimension(size) :: x, b
  real(dp),dimension(size) ::r
  real(dp) :: time
  integer :: c(2), c_rate, c_max
  real(dp) :: res
  integer :: i, j

  call init(size, a, x, b)
  call system_clock(c(1), c_rate, c_max)
  call bicgstab(size, a, x, b)
  call system_clock(c(2))
  time = 1.0d0*(c(2)-c(1))/c_rate
  write(6, *) "time[s]:", time
#ifdef _DEBUG
  write(6, *) "solution vector:"
  do i = 1, size
     write(6, '(1pe14.5)') x(i)
  end do
#endif
  write(6, *) "check the result: calc res = b - A*x"
  call b_minus_ax(size, a, x, b, r)
  call x_dot_y(size, r, r, res)
  res = sqrt(res)
  write(6, *) "residual:", res
  
  stop
end program main
