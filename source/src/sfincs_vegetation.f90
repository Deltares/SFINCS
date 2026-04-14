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
    integer :: nm, nmu, ip, iveg, k
    real*4  :: dh_veg, h_k, section_bottom, section_top
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
    ! For SFINCS precalculate cd * bstems * Nstems for each vertical section, as well as the vegetation height on uv points
    !
    ! TODO - do this now as pre-processing table    
    !
    if (vegetation) then
        !
        allocate(vegetation_stems_cd_width_density_uv(npuv, vegetation_vertical_segments)) ! product of cd, width and density on uv points
        !
        allocate(vegetation_stems_height_uv(npuv, vegetation_vertical_segments)) ! vegetation height on uv points
        !
        vegetation_stems_cd_width_density_uv = 0.0
        vegetation_stems_height_uv = 0.0
        !
        ! Interpolate vegetation properties from z-points to uv-points
        !
        do ip = 1, npuv
            !
            nm  = uv_index_z_nm(ip)
            nmu = uv_index_z_nmu(ip)
            !
            do iveg = 1, vegetation_vertical_segments
                !
                vegetation_stems_height_uv(ip,iveg) = 0.5*(vegetation_stems_height(nm,iveg)+vegetation_stems_height(nmu,iveg))
                !
                vegetation_stems_cd_width_density_uv(ip,iveg) = 0.5 * (0.5 * (vegetation_cd(nm,iveg) + vegetation_cd(nmu,iveg))) * (0.5 * (vegetation_stems_width(nm,iveg) + vegetation_stems_width(nmu,iveg))) * (0.5 * (vegetation_stems_density(nm,iveg) + vegetation_stems_density(nmu,iveg)))  / rhow
                !
                ! vegetation_stems_cd_width_density = 0.5 * cd * stems_width * stems_density / rhow, so everything that is precalculatable
                !
            enddo
            !
        enddo
        !
        ! Pre-compute lookup table: cumulative sum of cd*width*density at vegetation_nlookup equidistant depth levels
        ! Sections are stacked from the bed (consistent with swvegatt in snapwave_solver.f90):
        !   vegetation_cd_sum_table(ip, k) = sum_iveg( cd_wd(ip,iveg) * max(0, min(section_top_iveg, h_k) - section_bottom_iveg) )
        ! In compute_fluxes: fvm = table_lookup(ip, hu) * uv0 * |uv0|  (no inner do-loop needed)
        !
        allocate(vegetation_lookup_hmin_uv(npuv))
        allocate(vegetation_lookup_hmax_uv(npuv))
        allocate(vegetation_lookup_dh_uv(npuv))
        allocate(vegetation_cd_sum_table(npuv, 0:vegetation_nlookup))
        allocate(vegetation_cd_slope_table(npuv, 0:vegetation_nlookup-1))
        !
        vegetation_lookup_hmin_uv  = 0.0
        vegetation_lookup_hmax_uv  = 0.0
        vegetation_lookup_dh_uv    = 0.0
        vegetation_cd_sum_table    = 0.0
        vegetation_cd_slope_table  = 0.0
        !
        !$omp parallel do private( ip, k, iveg, dh_veg, h_k, section_bottom, section_top ) schedule( static )
        do ip = 1, npuv
            !
            ! Sections are stacked from the bed upward (consistent with swvegatt in snapwave_solver):
            !   section_bottom(iveg) = sum of all previous section heights
            !   section_top(iveg)    = section_bottom + vegetation_stems_height_uv(ip, iveg)
            ! hmax = total vegetation height = sum of all section thicknesses
            ! hmin = height of the bottom of the lowest section = 0 (all sections start at the bed)
            !
            vegetation_lookup_hmin_uv(ip) = 0.0
            vegetation_lookup_hmax_uv(ip) = sum(vegetation_stems_height_uv(ip,:))
            !
            if (vegetation_lookup_hmax_uv(ip) > 0.0) then
                dh_veg = vegetation_lookup_hmax_uv(ip) / real(vegetation_nlookup)
                vegetation_lookup_dh_uv(ip) = dh_veg
                do k = 1, vegetation_nlookup
                    h_k = k * dh_veg
                    section_bottom = 0.0
                    do iveg = 1, vegetation_vertical_segments
                        section_top = section_bottom + vegetation_stems_height_uv(ip, iveg)
                        vegetation_cd_sum_table(ip, k) = vegetation_cd_sum_table(ip, k) + vegetation_stems_cd_width_density_uv(ip, iveg) * max(0.0, min(section_top, h_k) - section_bottom)
                        section_bottom = section_top
                    enddo
                enddo
            endif
            !
        enddo
        !$omp end parallel do
        !
        ! Pre-compute slope between consecutive table entries to avoid the subtraction in compute_fluxes
        !
        !$omp parallel do private( ip, k ) schedule( static )
        do ip = 1, npuv
            do k = 0, vegetation_nlookup - 1
                vegetation_cd_slope_table(ip, k) = vegetation_cd_sum_table(ip, k+1) - vegetation_cd_sum_table(ip, k)
            enddo
        enddo
        !$omp end parallel do
        !
    endif
    !
    end subroutine    
        
end module    