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
    integer :: nm, iveg
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
        !varname = 'vegetation_vertical_segments' ! TODO: change naming netcdf file into this
        call read_netcdf_quadtree_get_dimension(veggiefile, varname, vegetation_vertical_segments) !ncfile, varname, varout)
        !        
        if (vegetation_vertical_segments > 4) then
            !
            call stop_sfincs('Error ! Maximum allowed vertical sections in vegetationfile is 4 !', 1)
        elseif(vegetation_vertical_segments == 0) then
            !
            call stop_sfincs('Error ! Prescribed vertical sections in vegetationfile is 0 !', 1)                        
            !
        endif        
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
        !
    endif
    !
    ! For SFINCS precalculate cd * bstems * Nstems for each vertical section - FIXME, already convert from np to uv?
    !
    if (vegetation) then
        !
        allocate(vegetation_stems_cd_width_density(np, vegetation_vertical_segments))
        !
        !allocate(vegetation_fvm(np, vegetation_vertical_segments)) 
        allocate(vegetation_fvm(npuv, vegetation_vertical_segments))            
        !
        vegetation_stems_cd_width_density = 0.0
        vegetation_fvm = 0.0
        fvm = 0.0              
        !
        do nm = 1, np
            !
            do iveg = 1, vegetation_vertical_segments
                !
                vegetation_stems_cd_width_density(nm,iveg) = vegetation_cd(nm,iveg) * vegetation_stems_width(nm,iveg) * vegetation_stems_density(nm,iveg)
                !
            enddo            
            !    
        enddo                 
        !            
    endif        
    !
    end subroutine    
        
end module    