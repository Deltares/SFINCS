module sfincs_vegetation

    use sfincs_log
    use sfincs_error

contains

    subroutine initialize_vegetation()
    !
    use sfincs_data
    use sfincs_ncinput
    !
    implicit none
    !
    integer :: nm, flag_idx, type_idx
    integer :: max_layers, nflags
    !
    logical :: ok
    !
    character*256 :: varname
    !
    integer,           allocatable :: flag_values(:)
    character(len=64), allocatable :: flag_meanings(:)
    integer,           allocatable :: type_to_lookup(:)   ! maps flag integer -> lookup row index
    real*4,            allocatable :: veg_cd_lookup(:,:)
    real*4,            allocatable :: veg_ah_lookup(:,:)
    real*4,            allocatable :: veg_bstems_lookup(:,:)
    real*4,            allocatable :: veg_nstems_lookup(:,:)
    !
    if (use_quadtree .eqv. .false.) then
        call stop_sfincs('Error ! Netcdf vegetation input format can only be specified for quadtree mesh model !', 1)
    endif
    !
    if (store_vegetation) then
        !
        write(logstr,'(a,a)')'Info    : reading vegetation file ',trim(veggiefile)
        call write_log(logstr, 0)
        !
        ok = check_file_exists(veggiefile, 'Vegetation file', .true.)
        !
        if (trim(veggiefile_toml) == 'none') then
            !
            ! Old path: all four parameter arrays stored directly in the NetCDF per grid cell
            !
            varname = 'nsec'
            call read_netcdf_quadtree_get_dimension(veggiefile, varname, vegetation_vertical_segments)
            !
            if (vegetation_vertical_segments > 64) then
                call stop_sfincs('Error ! vegetation_vertical_segments exceeds 64, check vegetationfile !', 1)
            endif
            !
            allocate(vegetation_stems_cd(np, vegetation_vertical_segments))
            allocate(vegetation_stems_height(np, vegetation_vertical_segments))
            allocate(vegetation_stems_diameter(np, vegetation_vertical_segments))
            allocate(vegetation_stems_density(np, vegetation_vertical_segments))
            !
            vegetation_stems_cd          = 0.0
            vegetation_stems_height = 0.0
            vegetation_stems_diameter  = 0.0
            vegetation_stems_density = 0.0
            !
            varname = 'snapwave_veg_Cd'
            call read_netcdf_quadtree_to_sfincs(veggiefile, varname, vegetation_stems_cd)
            !
            varname = 'snapwave_veg_ah'
            call read_netcdf_quadtree_to_sfincs(veggiefile, varname, vegetation_stems_height)
            !
            varname = 'snapwave_veg_bstems'
            call read_netcdf_quadtree_to_sfincs(veggiefile, varname, vegetation_stems_diameter)
            !
            varname = 'snapwave_veg_Nstems'
            call read_netcdf_quadtree_to_sfincs(veggiefile, varname, vegetation_stems_density)
            !
        else
            !
            ! New path: veggiefile holds only integer vegetation_type per cell (CF flag conventions)
            !           veggiefile_toml holds the parameter lookup table keyed by type name
            !
            ok = check_file_exists(veggiefile_toml, 'Vegetation TOML file', .true.)
            !
            write(logstr,'(a,a)')'Info    : reading vegetation TOML file ',trim(veggiefile_toml)
            call write_log(logstr, 0)
            !
            ! Read CF flag attributes: flag_values=[0,1,2,...] and flag_meanings=["none","seagrass",...]
            !
            varname = 'vegetation_type'
            call read_netcdf_flag_meanings(veggiefile, varname, &
                flag_values, flag_meanings, nflags)
            !
            write(logstr,'(a,i0,a)')'Info    : found ',nflags,' vegetation type(s) in NetCDF flag attributes'
            call write_log(logstr, 0)
            !
            ! Read integer vegetation_type per active cell
            !
            allocate(vegetation_type_index(np))
            vegetation_type_index = 0
            call read_netcdf_quadtree_integer(veggiefile, varname, vegetation_type_index)
            !
            ! Parse TOML lookup table: fills veg_*_lookup(nflags, max_layers), zero-padded
            !
            call read_vegetation_toml(flag_meanings, nflags, &
                veg_cd_lookup, veg_ah_lookup, veg_bstems_lookup, veg_nstems_lookup, max_layers)
            !
            vegetation_vertical_segments = max_layers
            !
            allocate(vegetation_stems_cd(np, vegetation_vertical_segments))
            allocate(vegetation_stems_height(np, vegetation_vertical_segments))
            allocate(vegetation_stems_diameter(np, vegetation_vertical_segments))
            allocate(vegetation_stems_density(np, vegetation_vertical_segments))
            !
            vegetation_stems_cd          = 0.0
            vegetation_stems_height = 0.0
            vegetation_stems_diameter  = 0.0
            vegetation_stems_density = 0.0
            !
            ! Build a direct map from flag integer value to lookup row to avoid per-cell linear search
            !
            allocate(type_to_lookup(0:maxval(flag_values)))
            type_to_lookup = -1
            do flag_idx = 1, nflags
                type_to_lookup(flag_values(flag_idx)) = flag_idx
            enddo
            !
            ! Expand lookup table to per-cell arrays
            !
            do nm = 1, np
                type_idx = vegetation_type_index(nm)
                if (type_idx == 0) cycle   ! type 0 = no vegetation
                !
                if (type_idx < 0 .or. type_idx > maxval(flag_values)) then
                    write(logstr,'(a,i0,a)')'Error    : vegetation_type ', type_idx, &
                        ' in NetCDF is outside the range of flag_values !'
                    call stop_sfincs(trim(logstr), 1)
                endif
                !
                flag_idx = type_to_lookup(type_idx)
                !
                if (flag_idx < 1) then
                    write(logstr,'(a,i0,a)')'Error    : vegetation_type ', type_idx, &
                        ' found in NetCDF but has no entry in the vegetation TOML file !'
                    call stop_sfincs(trim(logstr), 1)
                endif
                !
                vegetation_stems_cd(nm,:)           = veg_cd_lookup(flag_idx,:)
                vegetation_stems_height(nm,:) = veg_ah_lookup(flag_idx,:)
                vegetation_stems_diameter(nm,:)  = veg_bstems_lookup(flag_idx,:)
                vegetation_stems_density(nm,:) = veg_nstems_lookup(flag_idx,:)
                !
            enddo
            !
            deallocate(vegetation_type_index)
            deallocate(flag_values, flag_meanings, type_to_lookup)
            deallocate(veg_cd_lookup, veg_ah_lookup, veg_bstems_lookup, veg_nstems_lookup)
            !
        endif
        !
    endif
    !
    end subroutine


    subroutine read_vegetation_toml(flag_meanings, nflags, &
        veg_cd_lookup, veg_ah_lookup, veg_bstems_lookup, veg_nstems_lookup, max_layers)
    ! Parse the TOML vegetation lookup table and return parameter arrays (nflags x max_layers).
    ! Shorter types are zero-padded to max_layers so the output is rectangular.
    !
    use tomlf
    use sfincs_data
    !
    implicit none
    !
    integer,            intent(in)  :: nflags
    character(len=64),  intent(in)  :: flag_meanings(nflags)
    real*4, allocatable, intent(out) :: veg_cd_lookup(:,:)
    real*4, allocatable, intent(out) :: veg_ah_lookup(:,:)
    real*4, allocatable, intent(out) :: veg_bstems_lookup(:,:)
    real*4, allocatable, intent(out) :: veg_nstems_lookup(:,:)
    integer,            intent(out) :: max_layers
    !
    integer, parameter :: max_layers_limit = 64
    !
    type(toml_table), allocatable :: table
    type(toml_error), allocatable :: parse_error
    type(toml_table), pointer     :: veg_table, type_table
    type(toml_array), pointer     :: cd_arr, ah_arr, bstems_arr, nstems_arr
    !
    real*4, allocatable :: tmp_cd(:,:), tmp_ah(:,:), tmp_bstems(:,:), tmp_nstems(:,:)
    integer, allocatable :: nlayers_per_type(:)
    !
    integer    :: it, il, nlayers, stat
    real(kind=8) :: rval
    character(len=64) :: typename
    !
    allocate(tmp_cd(nflags, max_layers_limit))
    allocate(tmp_ah(nflags, max_layers_limit))
    allocate(tmp_bstems(nflags, max_layers_limit))
    allocate(tmp_nstems(nflags, max_layers_limit))
    allocate(nlayers_per_type(nflags))
    tmp_cd = 0.0;  tmp_ah = 0.0;  tmp_bstems = 0.0;  tmp_nstems = 0.0
    nlayers_per_type = 0
    !
    call toml_load(table, trim(veggiefile_toml), error=parse_error)
    !
    if (allocated(parse_error)) then
        write(logstr,'(a,a,a,a)')'Error    : failed to parse vegetation TOML file ', &
            trim(veggiefile_toml), ': ', trim(parse_error%message)
        call stop_sfincs(trim(logstr), 1)
    endif
    !
    if (.not. allocated(table)) then
        call stop_sfincs('Error    : vegetation TOML file is empty or could not be read !', 1)
    endif
    !
    ! Expect a [vegetation_type] table at the root level
    !
    nullify(veg_table)
    call get_value(table, 'vegetation_type', veg_table, stat=stat)
    !
    if (.not. associated(veg_table)) then
        call stop_sfincs('Error    : vegetation TOML file missing [vegetation_type] section !', 1)
    endif
    !
    ! Loop over all names from CF flag_meanings and read their parameter arrays
    !
    do it = 1, nflags
        typename = trim(flag_meanings(it))
        !
        ! Skip the no-vegetation entry; its cells are handled by the zero initialisation above
        !
        if (trim(typename) == 'none' .or. trim(typename) == '') cycle
        !
        nullify(type_table)
        call get_value(veg_table, trim(typename), type_table, stat=stat)
        !
        if (.not. associated(type_table)) then
            write(logstr,'(a,a,a)')'Warning : vegetation type "',trim(typename), &
                '" listed in NetCDF flag_meanings but absent from TOML - treated as no vegetation'
            call write_log(logstr, 0)
            cycle
        endif
        !
        ! --- snapwave_veg_Cd ---
        !
        nullify(cd_arr)
        call get_value(type_table, 'vegetation_stems_cd', cd_arr, stat=stat)
        !
        if (.not. associated(cd_arr)) then
            write(logstr,'(a,a,a)')'Error    : vegetation type "',trim(typename), &
                '" in TOML is missing vegetation_stems_cd array !'
            call stop_sfincs(trim(logstr), 1)
        endif
        !
        nlayers = len(cd_arr)
        !
        if (nlayers > max_layers_limit) then
            write(logstr,'(a,a,a,i0,a,i0,a)')'Error    : vegetation type "',trim(typename), &
                '" has ',nlayers,' layers which exceeds the limit of ',max_layers_limit,' !'
            call stop_sfincs(trim(logstr), 1)
        endif
        !
        nlayers_per_type(it) = nlayers
        !
        do il = 1, nlayers
            call get_value(cd_arr, il, rval, stat=stat)
            tmp_cd(it, il) = real(rval, kind=4)
        enddo
        !
        ! --- snapwave_veg_ah ---
        !
        nullify(ah_arr)
        call get_value(type_table, 'vegetation_stems_height', ah_arr, stat=stat)
        !
        if (.not. associated(ah_arr)) then
            write(logstr,'(a,a,a)')'Error    : vegetation type "',trim(typename), &
                '" in TOML is missing vegetation_stems_height array !'
            call stop_sfincs(trim(logstr), 1)
        endif
        !
        if (len(ah_arr) /= nlayers) then
            write(logstr,'(a,a,a)')'Error    : vegetation type "',trim(typename), &
                '" in TOML has inconsistent array lengths across parameters !'
            call stop_sfincs(trim(logstr), 1)
        endif
        !
        do il = 1, nlayers
            call get_value(ah_arr, il, rval, stat=stat)
            tmp_ah(it, il) = real(rval, kind=4)
        enddo
        !
        ! --- snapwave_veg_bstems ---
        !
        nullify(bstems_arr)
        call get_value(type_table, 'vegetation_stems_diameter', bstems_arr, stat=stat)
        !
        if (.not. associated(bstems_arr)) then
            write(logstr,'(a,a,a)')'Error    : vegetation type "',trim(typename), &
                '" in TOML is missing vegetation_stems_diameter array !'
            call stop_sfincs(trim(logstr), 1)
        endif
        !
        if (len(bstems_arr) /= nlayers) then
            write(logstr,'(a,a,a)')'Error    : vegetation type "',trim(typename), &
                '" in TOML has inconsistent array lengths across parameters !'
            call stop_sfincs(trim(logstr), 1)
        endif
        !
        do il = 1, nlayers
            call get_value(bstems_arr, il, rval, stat=stat)
            tmp_bstems(it, il) = real(rval, kind=4)
        enddo
        !
        ! --- snapwave_veg_Nstems ---
        !
        nullify(nstems_arr)
        call get_value(type_table, 'vegetation_stems_density', nstems_arr, stat=stat)
        !
        if (.not. associated(nstems_arr)) then
            write(logstr,'(a,a,a)')'Error    : vegetation type "',trim(typename), &
                '" in TOML is missing vegetation_stems_density array !'
            call stop_sfincs(trim(logstr), 1)
        endif
        !
        if (len(nstems_arr) /= nlayers) then
            write(logstr,'(a,a,a)')'Error    : vegetation type "',trim(typename), &
                '" in TOML has inconsistent array lengths across parameters !'
            call stop_sfincs(trim(logstr), 1)
        endif
        !
        do il = 1, nlayers
            call get_value(nstems_arr, il, rval, stat=stat)
            tmp_nstems(it, il) = real(rval, kind=4)
        enddo
        !
        write(logstr,'(a,a,a,i0,a)')'Info    : vegetation type "',trim(typename), &
            '" loaded with ',nlayers,' vertical layer(s)'
        call write_log(logstr, 0)
        !
    enddo
    !
    max_layers = maxval(nlayers_per_type)
    !
    if (max_layers == 0) then
        call stop_sfincs('Error    : no valid vegetation types found in TOML file !', 1)
    endif
    !
    allocate(veg_cd_lookup(nflags, max_layers))
    allocate(veg_ah_lookup(nflags, max_layers))
    allocate(veg_bstems_lookup(nflags, max_layers))
    allocate(veg_nstems_lookup(nflags, max_layers))
    !
    veg_cd_lookup     = tmp_cd(:, 1:max_layers)
    veg_ah_lookup     = tmp_ah(:, 1:max_layers)
    veg_bstems_lookup = tmp_bstems(:, 1:max_layers)
    veg_nstems_lookup = tmp_nstems(:, 1:max_layers)
    !
    deallocate(tmp_cd, tmp_ah, tmp_bstems, tmp_nstems, nlayers_per_type)
    !
    end subroutine


end module
