program test_bmi
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit
  use bmif_2_0,    only: bmi, BMI_SUCCESS, BMI_FAILURE, BMI_MAX_COMPONENT_NAME
  use sfincs_bmi2, only: sfincs_bmi
  implicit none

  type(sfincs_bmi) :: m
  integer :: s, i

  ! Must match BMI dummy: character(len=*), pointer
  character(len=BMI_MAX_COMPONENT_NAME), pointer :: comp_name
  character(len=BMI_MAX_COMPONENT_NAME), pointer :: in_names(:)
  character(len=BMI_MAX_COMPONENT_NAME), pointer :: out_names(:)
  integer :: nin, nout

  ! time
  double precision :: t0, t1, dt
  character(len=16) :: tunits

  call banner('=== SFINCS BMI TEST START ===')

  ! init with defaults
  s = m%initialize(''); call ok('initialize', s)

  nullify(comp_name, in_names, out_names)

  s = m%get_component_name(comp_name); call ok('get_component_name', s)
  if (s==BMI_SUCCESS .and. associated(comp_name)) write(output_unit,'(a)') 'component = '//trim(comp_name)

  s = m%get_input_item_count(nin);   call ok('get_input_item_count', s)
  s = m%get_output_item_count(nout); call ok('get_output_item_count', s)

  s = m%get_input_var_names(in_names);   call ok('get_input_var_names', s)
  s = m%get_output_var_names(out_names); call ok('get_output_var_names', s)

  write(output_unit,'(a,i0)') 'num inputs  = ', nin
  do i=1, nin
    write(output_unit,'(a,i0,2a)') '  in[',i,'] = ', trim(in_names(i))
  end do

  write(output_unit,'(a,i0)') 'num outputs = ', nout
  do i=1, nout
    write(output_unit,'(a,i0,2a)') '  out[',i,'] = ', trim(out_names(i))
  end do

  ! time API
  s = m%get_start_time(t0);     call ok('get_start_time', s)
  s = m%get_end_time(t1);       call ok('get_end_time',   s)
  s = m%get_time_step(dt);      call ok('get_time_step',  s)
  s = m%get_time_units(tunits); call ok('get_time_units', s)
  s = m%get_current_time(t0);   call ok('get_current_time', s)
  write(output_unit,'(a,f12.3,1x,a)') 'time step (', dt, ') units: '//trim(tunits)
  write(output_unit,'(a,f12.3)') 'start time : ', t0
  write(output_unit,'(a,f12.3)') 'end   time : ', t1

  ! variables to exercise
  call test_var('rain_rate',               m)
  call test_var('water_surface_elevation', m)
  call test_var('water_depth',             m)
  call test_var('velocity_x',              m)
  call test_var('velocity_y',              m)
  ! If you exposed 'bed_elevation' as an output, also:
  call test_var('bed_elevation',           m)

  ! drive model with rain, then check depth increased
  call drive_with_rain(m)

  s = m%finalize(); call ok('finalize', s)
  call banner('=== SFINCS BMI TEST DONE ===')

contains

  subroutine test_var(vname, m)
    character(len=*), intent(in) :: vname
    class(bmi),      intent(inout) :: m
    integer :: s, grid, nbytes, rank, gsize, isz
    character(len=32) :: vtype, vunits, vloc, gtype
    double precision, pointer :: dptr(:)
    double precision, allocatable :: d(:), gx(:), gy(:), gz(:)
    integer, allocatable :: shape(:)

    ! var metadata
    s = m%get_var_type(vname, vtype);     call ok('get_var_type('//trim(vname)//')', s)
    s = m%get_var_units(vname, vunits);   call ok('get_var_units('//trim(vname)//')', s)
    s = m%get_var_location(vname, vloc);  call ok('get_var_location('//trim(vname)//')', s)

    ! nbytes (guarded print)
    nbytes = -1
    s = m%get_var_nbytes(vname, nbytes)
    call ok('get_var_nbytes('//trim(vname)//')', s)
    if (s == BMI_SUCCESS) then
      write(output_unit,'(*(g0,1x))') '  var "',trim(vname),'" type=',trim(vtype),' units=',trim(vunits),' nbytes=',nbytes
    else
      write(output_unit,'(*(g0,1x))') '  var "',trim(vname),'" type=',trim(vtype),' units=',trim(vunits),' nbytes=(n/a)'
    end if

    ! grid meta
    s = m%get_var_grid(vname, grid);      call ok('get_var_grid('//trim(vname)//')', s)
    s = m%get_grid_type(grid, gtype);     call ok('get_grid_type', s)
    s = m%get_grid_rank(grid, rank);      call ok('get_grid_rank', s)
    s = m%get_grid_size(grid, gsize);     call ok('get_grid_size', s)
    write(output_unit,'(*(g0,1x))') '  grid id=',grid,' type=',trim(gtype),' rank=',rank,' size=',gsize

    allocate(shape(max(1,rank)))
    s = m%get_grid_shape(grid, shape);    call ok('get_grid_shape', s)
    write(output_unit,'(*(g0,1x))') '  shape:', (shape(i), i=1,size(shape))

    ! coords
    allocate(gx(gsize), gy(gsize), gz(1))
    s = m%get_grid_x(grid, gx);           call ok('get_grid_x', s)
    s = m%get_grid_y(grid, gy);           call ok('get_grid_y', s)
    s = m%get_grid_z(grid, gz);           call ok('get_grid_z', s)
    if (gsize > 0) write(output_unit,'(*(g0,1x))') '  x[1],x[end]=', gx(1), gx(gsize), '  y[1],y[end]=', gy(1), gy(gsize)

    ! itemsize cross-check (optional)
    isz = -1
    s = m%get_var_itemsize(vname, isz)
    if (s==BMI_SUCCESS .and. nbytes>=0 .and. gsize>=0) then
      if (isz>=0 .and. nbytes /= isz*gsize) then
        write(error_unit,'(*(g0,1x))') 'WARN: nbytes mismatch for ',trim(vname), ' got=', nbytes, ' expected=', isz*gsize
      end if
    end if

    ! pointer path (double)
    nullify(dptr)
    s = m%get_value_ptr_double(vname, dptr)
    if (s == BMI_SUCCESS .and. associated(dptr)) then
      write(output_unit,'(*(g0,1x))') '  ptr OK len=', size(dptr)
    else
      write(output_unit,'(a)') '  ptr not available'
    end if

    ! safe copy
    allocate(d(max(1,gsize)))
    d = -999.0d0
    s = m%get_value_double(vname, d(1:gsize)); call ok('get_value_double('//trim(vname)//')', s)
    if (gsize>0) write(output_unit,'(*(g0,1x))') '  sample "', trim(vname), '": ', d(1)

    deallocate(d, gx, gy, gz, shape)
  end subroutine test_var


  subroutine drive_with_rain(m)
    class(bmi), intent(inout) :: m
    integer :: s, grid, k, gsize
    double precision, allocatable :: rain(:), depth(:)
    double precision :: dt
    integer :: inds(3)
    double precision :: d3(3)

    s = m%get_var_grid('rain_rate', grid); call ok('get_var_grid(rain_rate)', s)
    s = m%get_grid_size(grid, gsize);      call ok('get_grid_size(rain)', s)

    allocate(rain(gsize), depth(gsize))

    do k=1, gsize
      rain(k) = 1.0d-4 * dble(k)/dble(max(1,gsize))
    end do
    s = m%set_value_double('rain_rate', rain); call ok('set_value_double(rain_rate)', s)

    s = m%get_time_step(dt); call ok('get_time_step', s)

    s = m%get_value_double('water_depth', depth); call ok('get_value_double(water_depth)', s)
    if (gsize>0) write(output_unit,'(a,f10.4)') '  depth[1] before = ', depth(1)

    do k=1,3
      s = m%update(); call ok('update', s)
    end do

    s = m%get_value_double('water_depth', depth); call ok('get_value_double(water_depth)', s)
    if (gsize>0) write(output_unit,'(a,f10.4)') '  depth[1] after  = ', depth(1)

    if (gsize >= 3) then
      inds = [1, (gsize+1)/2, gsize]
      d3 = -999.0d0
      s = m%get_value_at_indices_double('water_depth', d3, inds)
      call ok('get_value_at_indices_double(water_depth)', s)
      write(output_unit,'(*(g0,1x))') '  depth idx ', inds(1), ' = ', d3(1)
      write(output_unit,'(*(g0,1x))') '  depth idx ', inds(2), ' = ', d3(2)
      write(output_unit,'(*(g0,1x))') '  depth idx ', inds(3), ' = ', d3(3)
    end if

    deallocate(rain, depth)
  end subroutine drive_with_rain


  subroutine ok(where, status)
    character(len=*), intent(in) :: where
    integer,          intent(in) :: status
    if (status == BMI_SUCCESS) then
      write(output_unit,'(a)') 'OK:   '//trim(where)
    else
      write(error_unit,'(a)')  'FAIL: '//trim(where)//'  status='//trim(itoa(status))
    end if
  end subroutine ok

  subroutine banner(msg)
    character(len=*), intent(in) :: msg
    write(output_unit,'(a)') '-- '//trim(msg)//' --'
  end subroutine banner

  pure function itoa(i) result(s)
    integer, intent(in) :: i
    character(len=16) :: s
    write(s,'(i0)') i
    s = adjustl(s)
  end function itoa

end program test_bmi
