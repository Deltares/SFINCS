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
    integer :: nm
    !
    logical :: ok
    !
    character*256 :: varname
    !
    if (use_quadtree .eqv. .false.) then
        ! 
        call stop_sfincs('Error ! Netcdf vegetation input format can only be specified for quadtree mesh model !', 1)              
        !
    endif
    !
    if (store_vegetation) then !either SFINCS and/or SnapWave needs veggie input
        !              
        write(logstr,'(a,a)')'Info    : reading vegetation file ',trim(veggiefile)
        call write_log(logstr, 0)
        !
        ok = check_file_exists(veggiefile, 'Vegetation file', .true.)        
        !        
        ! Get dimension of vertical sections 
        !
        ! Call the generic quadtree nc file reader function
        varname = 'nsec'
        !varname = 'vegetation_vertical_segments' ! TODO: change input into this
        call read_netcdf_quadtree_get_dimension(veggiefile, varname, vegetation_vertical_segments) !ncfile, varname, varout)
        !        
        ! allocate variables
        allocate(vegetation_cd(np, vegetation_vertical_segments))
        allocate(vegetation_stems_height(np, vegetation_vertical_segments)) !=vegetation_ah
        allocate(vegetation_stems_width(np, vegetation_vertical_segments)) !=vegetation_bstems
        allocate(vegetation_stems_density(np, vegetation_vertical_segments)) !=vegetation_Nstems
        !
        vegetation_cd = 0.0
        vegetation_stems_height = 0.0
        vegetation_stems_width = 0.0
        vegetation_stems_density = 0.0        
        !
        ! Call the generic quadtree nc file reader function
        varname = 'snapwave_veg_Cd'
        !varname = 'vegegation_cd' ! TODO: change naming netcdf file into this
        call read_netcdf_quadtree_to_sfincs(veggiefile, varname, vegetation_cd) !ncfile, varname, varout)
        !
        ! Call the generic quadtree nc file reader function
        varname = 'snapwave_veg_ah'
        !varname = 'vegetation_stems_height' ! TODO: change naming netcdf file into this
        call read_netcdf_quadtree_to_sfincs(veggiefile, varname, vegetation_stems_height) !ncfile, varname, varout)
        !
        ! Call the generic quadtree nc file reader function
        varname = 'snapwave_veg_bstems'
        !varname = 'vegetation_stems_width' ! TODO: change naming netcdf file into this
        call read_netcdf_quadtree_to_sfincs(veggiefile, varname, vegetation_stems_width) !ncfile, varname, varout)
        !
        ! Call the generic quadtree nc file reader function
        varname = 'snapwave_veg_Nstems'
        !varname = 'vegetation_stems_density' ! TODO: change naming netcdf file into this
        call read_netcdf_quadtree_to_sfincs(veggiefile, varname, vegetation_stems_density) !ncfile, varname, varout)

    endif
    !
    end subroutine    
        
end module    