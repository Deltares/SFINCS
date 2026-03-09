! Core model logic for SFINCS used by BMI wrapper
module sfincs_core
  use iso_fortran_env, only: real64
  implicit none

  type :: Sfincs_State
    integer :: m, n
    real(real64) :: dx, dy
    real(real64) :: dt, current_time
    real(real64), allocatable :: zs(:,:), qx(:,:), qy(:,:)
    real(real64), allocatable :: elevation(:,:), mask(:,:)
  end type Sfincs_State

contains

  subroutine init_sfincs(state, config_file)
    type(Sfincs_State), intent(inout) :: state
    character(len=*), intent(in) :: config_file

    ! For demonstration, fixed values are used
    state%m = 100
    state%n = 100
    state%dx = 1.0
    state%dy = 1.0
    state%dt = 60.0
    state%current_time = 0.0

    allocate(state%zs(state%m, state%n))
    allocate(state%qx(state%m, state%n))
    allocate(state%qy(state%m, state%n))
    allocate(state%elevation(state%m, state%n))
    allocate(state%mask(state%m, state%n))

    state%zs = 0.0
    state%qx = 0.0
    state%qy = 0.0
    state%elevation = 1.0
    state%mask = 1.0
  end subroutine init_sfincs

  subroutine update_sfincs(state)
    type(Sfincs_State), intent(inout) :: state
    integer :: i, j

    ! Dummy update: add 0.001 m/s flow to qx and qy, advance zs
    do i = 1, state%m
      do j = 1, state%n
        state%qx(i,j) = state%qx(i,j) + 0.001
        state%qy(i,j) = state%qy(i,j) + 0.001
        state%zs(i,j) = state%zs(i,j) + state%dt * (state%qx(i,j) + state%qy(i,j)) * 1.0e-6
      end do
    end do

    state%current_time = state%current_time + state%dt
  end subroutine update_sfincs

  subroutine finalize_sfincs(state)
    type(Sfincs_State), intent(inout) :: state

    if (allocated(state%zs)) deallocate(state%zs)
    if (allocated(state%qx)) deallocate(state%qx)
    if (allocated(state%qy)) deallocate(state%qy)
    if (allocated(state%elevation)) deallocate(state%elevation)
    if (allocated(state%mask)) deallocate(state%mask)
  end subroutine finalize_sfincs

end module sfincs_core

