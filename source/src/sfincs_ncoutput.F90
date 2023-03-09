#define NF90(nf90call) call handle_err(nf90call,__FILE__,__LINE__)     
module sfincs_ncoutput
   !
   use netcdf 
   !
   implicit none
   !
   type map_type
      !
      integer :: ncid   
      integer :: n_dimid, m_dimid, corner_n_dimid, corner_m_dimid
      integer :: time_dimid 
      integer :: timemax_dimid 
      integer :: runtime_dimid
      integer :: corner_x_varid, corner_y_varid, face_x_varid, face_y_varid, crs_varid, grid_varid 
      integer :: zb_varid, msk_varid, qinf_varid
      integer :: time_varid, timemax_varid
      integer :: zs_varid, zsmax_varid, h_varid, u_varid, v_varid, tmax_varid, Seff_varid 
      integer :: hmax_varid, vmax_varid, cumprcp_varid, cuminf_varid, windmax_varid
      integer :: patm_varid, wind_u_varid, wind_v_varid, precip_varid        
      integer :: hm0_varid, hm0ig_varid
      integer :: fwx_varid, fwy_varid
      integer :: zsm_varid
      integer :: inp_varid, total_runtime_varid, average_dt_varid
      !
      integer :: mesh2d_varid
      integer :: mesh2d_node_x_varid, mesh2d_node_y_varid
      integer :: mesh2d_face_nodes_varid
      integer :: nmesh2d_node_dimid
      integer :: nmesh2d_face_dimid
      integer :: max_nmesh2d_face_nodes_dimid
      !
   end type
   !
   type his_type
      !
      integer :: ncid   
      integer :: time_dimid 
      integer :: points_dimid, pointnamelength_dimid
      integer :: crosssections_dimid, structures_dimid
      integer :: runtime_dimid
      integer :: point_x_varid, point_y_varid, station_x_varid, station_y_varid, crs_varid, qinf_varid, S_varid  
      integer :: station_id_varid, station_name_varid
      integer :: crosssection_id_varid, crosssection_name_varid
      integer :: structure_height_varid, structure_x_varid, structure_y_varid
      integer :: zb_varid
      integer :: time_varid
      integer :: zs_varid, h_varid, u_varid, v_varid, prcp_varid, discharge_varid, uvmag_varid, uvdir_varid
      integer :: patm_varid, wind_speed_varid, wind_dir_varid
      integer :: inp_varid, total_runtime_varid, average_dt_varid  
      integer :: hm0_varid, hm0ig_varid, zsm_varid, tp_varid, wavdir_varid, dirspr_varid
      !
   end type
   !
   type(map_type) :: map_file
   type(his_type) :: his_file
   !
   real*4, parameter :: FILL_VALUE = -99999.0

contains

   subroutine ncoutput_regular_map_init()
   !
   ! 1. Initialise dimensions/variables/attributes
   ! 2. write grid/msk/zb to file
   !
   use sfincs_date
   use sfincs_data   
   !
   implicit none   
   !   
   integer                      :: nm, n, m, ntmx
   !
   real*4, dimension(:,:), allocatable :: zsg
   real*4, dimension(:,:), allocatable :: xz
   real*4, dimension(:,:), allocatable :: yz
   real*4, dimension(:,:), allocatable :: xg
   real*4, dimension(:,:), allocatable :: yg
   !
   !write(*,*) trim(nf90_inq_libvers())
   !
   NF90(nf90_create('sfincs_map.nc', NF90_CLOBBER, map_file%ncid))
   !
   ! Create dimensions
   ! grid, time, points
   ! do mmax/nmax-2 to not write away dummy cells
   NF90(nf90_def_dim(map_file%ncid, 'n', nmax, map_file%n_dimid)) ! rows 
   NF90(nf90_def_dim(map_file%ncid, 'm', mmax, map_file%m_dimid)) ! columns
   NF90(nf90_def_dim(map_file%ncid, 'corner_n', nmax + 1, map_file%corner_n_dimid)) ! rows of corners
   NF90(nf90_def_dim(map_file%ncid, 'corner_m', mmax + 1, map_file%corner_m_dimid)) ! columns of corners   
   NF90(nf90_def_dim(map_file%ncid, 'time', NF90_UNLIMITED, map_file%time_dimid)) ! time
   ntmx = max(int((t1out - t0out)/dtmaxout), 1)
   NF90(nf90_def_dim(map_file%ncid, 'timemax', ntmx, map_file%timemax_dimid)) ! time
   NF90(nf90_def_dim(map_file%ncid, 'runtime', 1, map_file%runtime_dimid)) ! total_runtime, average_dt       
   !
   ! Some metadata attributes 
   NF90(nf90_put_att(map_file%ncid,nf90_global, "Conventions", "Conventions = 'CF-1.6, SGRID-0.3")) 
   NF90(nf90_put_att(map_file%ncid,nf90_global, "Build-Revision-Date-Netcdf-library", trim(nf90_inq_libvers()))) ! version of netcdf library
   NF90(nf90_put_att(map_file%ncid,nf90_global, "Producer", "SFINCS model: Super-Fast INundation of CoastS"))
   NF90(nf90_put_att(map_file%ncid,nf90_global, "Build-Revision", trim(build_revision))) 
   NF90(nf90_put_att(map_file%ncid,nf90_global, "Build-Date", trim(build_date)))
   NF90(nf90_put_att(map_file%ncid,nf90_global, "title", "SFINCS map netcdf output"))   
   !
   ! add input params for reproducability
   !
   call ncoutput_add_params(map_file%ncid,map_file%inp_varid)   
   !
   !! Create variables
   !
   ! Domain
   !
   NF90(nf90_def_var(map_file%ncid, 'x', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid/), map_file%face_x_varid)) ! location of zb, zs etc. in cell centre
   NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, '_FillValue', FILL_VALUE))         
   if (crsgeo) then
      NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'units', 'degrees'))
      NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'standard_name', 'longitude'))
   else
      NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'standard_name', 'projection_x_coordinate'))
   endif
   NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'long_name', 'face_x'))   
   NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'grid_mapping', 'crs'))   
   NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'grid', 'sfincsgrid'))   
   !
   NF90(nf90_def_var(map_file%ncid, 'y', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid/), map_file%face_y_varid)) ! location of zb, zs etc. in cell centre  
   NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, '_FillValue', FILL_VALUE))            
   if (crsgeo) then
      NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'units', 'degrees'))
      NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'standard_name', 'latitude'))
   else
      NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'standard_name', 'projection_y_coordinate'))
   endif
   NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'long_name', 'face_y'))
   NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'grid_mapping', 'crs'))   
   NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'grid', 'sfincsgrid'))      
   !
   NF90(nf90_def_var(map_file%ncid, 'corner_x', NF90_FLOAT, (/map_file%corner_m_dimid, map_file%corner_n_dimid/), map_file%corner_x_varid)) ! location of u points in cell corner
   NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, '_FillValue', FILL_VALUE))         
   if (crsgeo) then
      NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'units', 'degrees'))
      NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'standard_name', 'longitude'))
   else
      NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'standard_name', 'projection_x_coordinate'))
   endif
   NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'long_name', 'corner_x'))
   NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'grid_mapping', 'crs'))   
   NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'grid', 'sfincsgrid'))     
   !
   NF90(nf90_def_var(map_file%ncid, 'corner_y', NF90_FLOAT, (/map_file%corner_m_dimid, map_file%corner_n_dimid/), map_file%corner_y_varid)) ! location of v points in cell corner
   NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, '_FillValue', FILL_VALUE))         
   if (crsgeo) then
      NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'units', 'degrees'))
      NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'standard_name', 'latitude'))
   else
      NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'standard_name', 'projection_y_coordinate'))
   endif
   NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'long_name', 'corner_y'))
   NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'grid_mapping', 'crs'))   
   NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'grid', 'sfincsgrid'))   
   !
   NF90(nf90_def_var(map_file%ncid, 'crs', NF90_INT, map_file%crs_varid)) ! For EPSG code
   NF90(nf90_put_att(map_file%ncid, map_file%crs_varid, 'EPSG', '-'))      
   NF90(nf90_put_att(map_file%ncid, map_file%crs_varid, 'epsg_code', 'EPSG:' // trim(epsg_code) ))   !--> add epsg_code like FEWS wants   
   !
   NF90(nf90_def_var(map_file%ncid, 'sfincsgrid', NF90_INT, map_file%grid_varid)) ! For neat grid clarification
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'cf_role', 'grid_topology'))
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'topology_dimension', 2))
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'node_dimensions', 'n m')) !or n m?   
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'face_dimensions', 'n: m:')) 
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'corner_dimensions', 'corner_n: corner_m:'))   
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'face_coordinates', 'x y'))
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'corner_coordinates', 'corner_x corner_y'))   
   !
   NF90(nf90_def_var(map_file%ncid, 'msk', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid/), map_file%msk_varid)) ! input msk value in cell centre
   NF90(nf90_put_att(map_file%ncid, map_file%msk_varid, '_FillValue', FILL_VALUE))      
   NF90(nf90_put_att(map_file%ncid, map_file%msk_varid, 'units', '-'))
   NF90(nf90_put_att(map_file%ncid, map_file%msk_varid, 'standard_name', 'land_binary_mask')) ! land_binary_mask but with added boundary=2
   NF90(nf90_put_att(map_file%ncid, map_file%msk_varid, 'long_name', 'msk_active_cells')) 
   NF90(nf90_put_att(map_file%ncid, map_file%msk_varid, 'description', 'inactive=0, active=1, normal_boundary=2, outflow_boundary=3'))    
   NF90(nf90_put_att(map_file%ncid, map_file%msk_varid, 'coordinates', 'x y'))
   !
   ! Infiltration map
   !
   if (infiltration) then
       NF90(nf90_def_var(map_file%ncid, 'qinf', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid/), map_file%qinf_varid)) 
       NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, '_FillValue', FILL_VALUE))
       NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'coordinates', 'x y'))
       if (inftype == 'cna') then
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'standard_name', 'S')) 
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'long_name', 'moisture storage (S) capacity as computed from the curve number')) 
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'units', 'm'))
       elseif (inftype == 'cnb') then
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'standard_name', 'Smax')) 
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'long_name', 'maximum moisture storage (Smax) capacity as computed from the curve number')) 
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'units', 'm'))
       else
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'standard_name', 'qinf')) 
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'long_name', 'infiltration rate - constant in time')) 
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'units', 'mm h-1'))
       endif
   endif
   !
   NF90(nf90_def_var(map_file%ncid, 'zb', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid/), map_file%zb_varid)) ! bed level in cell centre
   NF90(nf90_put_att(map_file%ncid, map_file%zb_varid, '_FillValue', FILL_VALUE))   
   NF90(nf90_put_att(map_file%ncid, map_file%zb_varid, 'units', 'm'))
   NF90(nf90_put_att(map_file%ncid, map_file%zb_varid, 'standard_name', 'altitude'))
   NF90(nf90_put_att(map_file%ncid, map_file%zb_varid, 'long_name', 'bed_level_above_reference_level'))   
   NF90(nf90_put_att(map_file%ncid, map_file%zb_varid, 'coordinates', 'x y'))   
   !
   ! Time variables   
   !
   trefstr_iso8601 = date_to_iso8601(trefstr)
   NF90(nf90_def_var(map_file%ncid, 'time', NF90_FLOAT, (/map_file%time_dimid/), map_file%time_varid)) ! time
   NF90(nf90_put_att(map_file%ncid, map_file%time_varid, 'units', 'seconds since ' // trim(trefstr_iso8601) ))  ! time stamp following ISO 8601
   NF90(nf90_put_att(map_file%ncid, map_file%time_varid, 'standard_name', 'time'))     
   NF90(nf90_put_att(map_file%ncid, map_file%time_varid, 'long_name', 'time_in_seconds_since_' // trim(trefstr_iso8601) ))  
   !
   ! Time varying map output
   !
   NF90(nf90_def_var(map_file%ncid, 'zs', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%zs_varid)) ! time-varying water level map
   NF90(nf90_put_att(map_file%ncid, map_file%zs_varid, '_FillValue', FILL_VALUE))
   NF90(nf90_put_att(map_file%ncid, map_file%zs_varid, 'units', 'm'))
   NF90(nf90_put_att(map_file%ncid, map_file%zs_varid, 'standard_name', 'sea_surface_height_above_mean_sea_level')) 
   NF90(nf90_put_att(map_file%ncid, map_file%zs_varid, 'long_name', 'water_level'))  
   NF90(nf90_put_att(map_file%ncid, map_file%zs_varid, 'coordinates', 'x y'))
   !
   if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then   
      NF90(nf90_def_var(map_file%ncid, 'h', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%h_varid)) ! time-varying water depth map
      NF90(nf90_put_att(map_file%ncid, map_file%h_varid, '_FillValue', FILL_VALUE))
      NF90(nf90_put_att(map_file%ncid, map_file%h_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%h_varid, 'standard_name', 'depth')) 
      NF90(nf90_put_att(map_file%ncid, map_file%h_varid, 'long_name', 'water_depth'))     
      NF90(nf90_put_att(map_file%ncid, map_file%h_varid, 'coordinates', 'x y'))
   endif
   !
   ! Velocity is optional
   !
   if (store_velocity) then
      !
      NF90(nf90_def_var(map_file%ncid, 'u', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%u_varid)) ! time-varying u map 
      NF90(nf90_put_att(map_file%ncid, map_file%u_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(map_file%ncid, map_file%u_varid, 'units', 'm s-1'))
      NF90(nf90_put_att(map_file%ncid, map_file%u_varid, 'standard_name', 'eastward_sea_water_velocity')) ! not truly eastward when rotated, eastward_sea_water_velocity
      NF90(nf90_put_att(map_file%ncid, map_file%u_varid, 'long_name', 'flow_velocity_x_direction'))     
      NF90(nf90_put_att(map_file%ncid, map_file%u_varid, 'coordinates', 'x y'))
      !
      NF90(nf90_def_var(map_file%ncid, 'v', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%v_varid)) ! time-varying u map 
      NF90(nf90_put_att(map_file%ncid, map_file%v_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(map_file%ncid, map_file%v_varid, 'units', 'm s-1'))
      NF90(nf90_put_att(map_file%ncid, map_file%v_varid, 'standard_name', 'northward_sea_water_velocity')) ! not truly eastward when rotated, eastward_sea_water_velocity
      NF90(nf90_put_att(map_file%ncid, map_file%v_varid, 'long_name', 'flow_velocity_y_direction'))     
      NF90(nf90_put_att(map_file%ncid, map_file%v_varid, 'coordinates', 'x y'))
   endif
   !
   ! Store S_effective (only for CN method with recovery)
   !
   if (inftype == 'cnb') then
      NF90(nf90_def_var(map_file%ncid, 'Seff', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%Seff_varid)) ! time-varying S
      NF90(nf90_put_att(map_file%ncid, map_file%Seff_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(map_file%ncid, map_file%Seff_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%Seff_varid, 'standard_name', 'Se')) 
      NF90(nf90_put_att(map_file%ncid, map_file%Seff_varid, 'long_name', 'current moisture storage (Se) capacity'))     
      NF90(nf90_put_att(map_file%ncid, map_file%Seff_varid, 'coordinates', 'corner_x corner_y'))
   endif
   !
   ! Time varying spatial output max
   !
   if (store_maximum_waterlevel) then
      NF90(nf90_def_var(map_file%ncid, 'timemax', NF90_FLOAT, (/map_file%timemax_dimid/), map_file%timemax_varid)) ! time
      NF90(nf90_put_att(map_file%ncid, map_file%timemax_varid, 'units', 'seconds since ' // trim(trefstr_iso8601) ))  ! time stamp following ISO 8601
      NF90(nf90_put_att(map_file%ncid, map_file%timemax_varid, 'standard_name', 'time'))     
      NF90(nf90_put_att(map_file%ncid, map_file%timemax_varid, 'long_name', 'time_in_seconds_since_' // trim(trefstr_iso8601) ))  
   endif
   !
   if (store_maximum_waterlevel) then
      NF90(nf90_def_var(map_file%ncid, 'zsmax', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%timemax_dimid/), map_file%zsmax_varid)) ! time-varying maximum water level map
      NF90(nf90_put_att(map_file%ncid, map_file%zsmax_varid, '_FillValue', FILL_VALUE))
      NF90(nf90_put_att(map_file%ncid, map_file%zsmax_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%zsmax_varid, 'standard_name', 'maximum of sea_surface_height_above_mean_sea_level')) 
      NF90(nf90_put_att(map_file%ncid, map_file%zsmax_varid, 'long_name', 'maximum_water_level'))
      NF90(nf90_put_att(map_file%ncid, map_file%zsmax_varid, 'coordinates', 'x y'))
   endif
   !
   if (store_cumulative_precipitation) then
      NF90(nf90_def_var(map_file%ncid, 'cumprcp', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%timemax_dimid/), map_file%cumprcp_varid)) ! time-varying maximum water level map
      NF90(nf90_put_att(map_file%ncid, map_file%cumprcp_varid, '_FillValue', FILL_VALUE))          
      NF90(nf90_put_att(map_file%ncid, map_file%cumprcp_varid, 'units', 'mm'))
      NF90(nf90_put_att(map_file%ncid, map_file%cumprcp_varid, 'long_name', 'cumulative_precipitation_depth')) 
      NF90(nf90_put_att(map_file%ncid, map_file%cumprcp_varid, 'standard_name', 'cumulative_precipitation_depth')) 
      NF90(nf90_put_att(map_file%ncid, map_file%cumprcp_varid, 'cell_methods', 'time: sum'))       
      NF90(nf90_put_att(map_file%ncid, map_file%cumprcp_varid, 'coordinates', 'x y'))
   endif
   !
   if (store_twet) then
      NF90(nf90_def_var(map_file%ncid, 'tmax', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%timemax_dimid/), map_file%tmax_varid)) ! time-varying duration wet cell
      NF90(nf90_put_att(map_file%ncid, map_file%tmax_varid, '_FillValue', FILL_VALUE))
      NF90(nf90_put_att(map_file%ncid, map_file%tmax_varid, 'units', 'seconds'))
      NF90(nf90_put_att(map_file%ncid, map_file%tmax_varid, 'standard_name', 'duration cell is considered wet')) 
      NF90(nf90_put_att(map_file%ncid, map_file%tmax_varid, 'long_name', 'duration_wet_cell'))  
      NF90(nf90_put_att(map_file%ncid, map_file%tmax_varid, 'cell_methods', 'time: sum'))    
      NF90(nf90_put_att(map_file%ncid, map_file%tmax_varid, 'coordinates', 'x y'))
   endif
   !
   if (store_maximum_waterlevel) then
      if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then   
         NF90(nf90_def_var(map_file%ncid, 'hmax', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%timemax_dimid/), map_file%hmax_varid)) ! time-varying maximum water depth map
         NF90(nf90_put_att(map_file%ncid, map_file%hmax_varid, '_FillValue', FILL_VALUE))   
         NF90(nf90_put_att(map_file%ncid, map_file%hmax_varid, 'units', 'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%hmax_varid, 'standard_name', 'sea_floor_depth_below_sea_surface')) 
         NF90(nf90_put_att(map_file%ncid, map_file%hmax_varid, 'long_name', 'maximum_water_depth')) 
         NF90(nf90_put_att(map_file%ncid, map_file%hmax_varid, 'cell_methods', 'time: maximum'))    
         NF90(nf90_put_att(map_file%ncid, map_file%hmax_varid, 'coordinates', 'x y'))
      endif
   endif
   !
   if (store_maximum_velocity) then
      NF90(nf90_def_var(map_file%ncid, 'vmax', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%timemax_dimid/), map_file%vmax_varid)) ! maximum flow velocity map
      NF90(nf90_put_att(map_file%ncid, map_file%vmax_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(map_file%ncid, map_file%vmax_varid, 'units', 'm s-1'))
      NF90(nf90_put_att(map_file%ncid, map_file%vmax_varid, 'standard_name', 'maximum_flow_velocity')) ! no standard name available
      NF90(nf90_put_att(map_file%ncid, map_file%vmax_varid, 'long_name', 'maximum_flow_velocity')) 
      NF90(nf90_put_att(map_file%ncid, map_file%vmax_varid, 'cell_methods', 'time: maximum'))
      NF90(nf90_put_att(map_file%ncid, map_file%vmax_varid, 'coordinates', 'x y'))
   endif
   !
   if (store_cumulative_precipitation) then
       NF90(nf90_def_var(map_file%ncid, 'cuminf', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%timemax_dimid/), map_file%cuminf_varid)) ! time-varying cumulative infiltration map
       NF90(nf90_put_att(map_file%ncid, map_file%cuminf_varid, '_FillValue', FILL_VALUE))          
       NF90(nf90_put_att(map_file%ncid, map_file%cuminf_varid, 'units', 'm'))
       NF90(nf90_put_att(map_file%ncid, map_file%cuminf_varid, 'long_name', 'cumulative_infiltration_depth')) 
       NF90(nf90_put_att(map_file%ncid, map_file%cuminf_varid, 'cell_methods', 'time: sum'))     
       NF90(nf90_put_att(map_file%ncid, map_file%cuminf_varid, 'coordinates', 'x y'))
   endif
   !
   if (store_meteo) then  
      !
      if (wind) then
         !
         NF90(nf90_def_var(map_file%ncid, 'wind_u', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%wind_u_varid)) ! cumulative precipitation map
         NF90(nf90_put_att(map_file%ncid, map_file%wind_u_varid, '_FillValue', FILL_VALUE))          
         NF90(nf90_put_att(map_file%ncid, map_file%wind_u_varid, 'units', 'm s-1'))
         NF90(nf90_put_att(map_file%ncid, map_file%wind_u_varid, 'standard_name', 'eastward_wind'))
         NF90(nf90_put_att(map_file%ncid, map_file%wind_u_varid, 'long_name', 'wind_speed_u')) 
         NF90(nf90_put_att(map_file%ncid, map_file%wind_u_varid, 'coordinates', 'x y'))
         !
         NF90(nf90_def_var(map_file%ncid, 'wind_v', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%wind_v_varid)) ! cumulative precipitation map
         NF90(nf90_put_att(map_file%ncid, map_file%wind_v_varid, '_FillValue', FILL_VALUE))          
         NF90(nf90_put_att(map_file%ncid, map_file%wind_v_varid, 'units', 'm s-1'))
         NF90(nf90_put_att(map_file%ncid, map_file%wind_v_varid, 'standard_name', 'northward_wind'))
         NF90(nf90_put_att(map_file%ncid, map_file%wind_v_varid, 'long_name', 'wind_speed_v')) 
         NF90(nf90_put_att(map_file%ncid, map_file%wind_v_varid, 'coordinates', 'x y'))
         !
         if (meteo3d .and. store_wind_max) then  
            NF90(nf90_def_var(map_file%ncid, 'windmax', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid/), map_file%windmax_varid)) ! maximum wind speed 
            NF90(nf90_put_att(map_file%ncid, map_file%windmax_varid, '_FillValue', FILL_VALUE))          
            NF90(nf90_put_att(map_file%ncid, map_file%windmax_varid, 'units', 'm s-1'))
            NF90(nf90_put_att(map_file%ncid, map_file%windmax_varid, 'long_name', 'maximum_wind_speed')) 
            NF90(nf90_put_att(map_file%ncid, map_file%windmax_varid, 'cell_methods', 'time: sum'))       
            NF90(nf90_put_att(map_file%ncid, map_file%windmax_varid, 'coordinates', 'x y'))
        endif
         !
      endif
      !
      if (patmos) then
         !
         NF90(nf90_def_var(map_file%ncid, 'surface_air_pressure', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%patm_varid)) ! cumulative precipitation map
         NF90(nf90_put_att(map_file%ncid, map_file%patm_varid, '_FillValue', FILL_VALUE))          
         NF90(nf90_put_att(map_file%ncid, map_file%patm_varid, 'units', 'N m-2'))
         NF90(nf90_put_att(map_file%ncid, map_file%patm_varid, 'standard_name', 'surface_air_pressure'))
         NF90(nf90_put_att(map_file%ncid, map_file%patm_varid, 'long_name', 'surface_air_pressure')) 
         NF90(nf90_put_att(map_file%ncid, map_file%patm_varid, 'coordinates', 'x y'))
         !
      endif   
      !
      if (precip) then
         !
         NF90(nf90_def_var(map_file%ncid, 'precipitation_rate', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%precip_varid)) ! cumulative precipitation map
         NF90(nf90_put_att(map_file%ncid, map_file%precip_varid, '_FillValue', FILL_VALUE))          
         NF90(nf90_put_att(map_file%ncid, map_file%precip_varid, 'units', 'mm h-1'))
         NF90(nf90_put_att(map_file%ncid, map_file%precip_varid, 'standard_name', 'precipitation_rate'))
         NF90(nf90_put_att(map_file%ncid, map_file%precip_varid, 'long_name', 'precipitation_rate')) 
         NF90(nf90_put_att(map_file%ncid, map_file%precip_varid, 'coordinates', 'x y'))
         !
      endif   
      !
   endif
   !
   if (snapwave) then  
      !
      NF90(nf90_def_var(map_file%ncid, 'hm0', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%hm0_varid))
      NF90(nf90_put_att(map_file%ncid, map_file%hm0_varid, '_FillValue', FILL_VALUE))          
      NF90(nf90_put_att(map_file%ncid, map_file%hm0_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%hm0_varid, 'standard_name', 'hm0_wave_height'))
      NF90(nf90_put_att(map_file%ncid, map_file%hm0_varid, 'long_name', 'Hm0 wave height')) 
      NF90(nf90_put_att(map_file%ncid, map_file%hm0_varid, 'coordinates', 'x y'))
      !
      NF90(nf90_def_var(map_file%ncid, 'hm0ig', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%hm0ig_varid))
      NF90(nf90_put_att(map_file%ncid, map_file%hm0ig_varid, '_FillValue', FILL_VALUE))          
      NF90(nf90_put_att(map_file%ncid, map_file%hm0ig_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%hm0ig_varid, 'standard_name', 'hm0_ig_wave_height'))
      NF90(nf90_put_att(map_file%ncid, map_file%hm0ig_varid, 'long_name', 'Hm0 infragravity wave height')) 
      NF90(nf90_put_att(map_file%ncid, map_file%hm0ig_varid, 'coordinates', 'x y'))
      !
      if (store_wave_forces) then
         !
         NF90(nf90_def_var(map_file%ncid, 'fwx', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%fwx_varid))
         NF90(nf90_put_att(map_file%ncid, map_file%fwx_varid, '_FillValue', FILL_VALUE))          
         NF90(nf90_put_att(map_file%ncid, map_file%fwx_varid, 'units', 'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%fwx_varid, 'standard_name', 'wave_force_x'))
         NF90(nf90_put_att(map_file%ncid, map_file%fwx_varid, 'long_name', 'Wave force in x-direction')) 
         NF90(nf90_put_att(map_file%ncid, map_file%fwx_varid, 'coordinates', 'x y'))
         !
         NF90(nf90_def_var(map_file%ncid, 'fwy', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%fwy_varid))
         NF90(nf90_put_att(map_file%ncid, map_file%fwy_varid, '_FillValue', FILL_VALUE))          
         NF90(nf90_put_att(map_file%ncid, map_file%fwy_varid, 'units', 'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%fwy_varid, 'standard_name', 'wave_force_y'))
         NF90(nf90_put_att(map_file%ncid, map_file%fwy_varid, 'long_name', 'Wave force in y-direction')) 
         NF90(nf90_put_att(map_file%ncid, map_file%fwy_varid, 'coordinates', 'x y'))
         !
      endif
      !
      if (wavemaker) then
         !
         NF90(nf90_def_var(map_file%ncid, 'zsm', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid, map_file%time_dimid/), map_file%zsm_varid))
         NF90(nf90_put_att(map_file%ncid, map_file%zsm_varid, '_FillValue', FILL_VALUE))          
         NF90(nf90_put_att(map_file%ncid, map_file%zsm_varid, 'units', 'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%zsm_varid, 'standard_name', 'mean_water_level'))
         NF90(nf90_put_att(map_file%ncid, map_file%zsm_varid, 'long_name', 'Filtered water level')) 
         NF90(nf90_put_att(map_file%ncid, map_file%zsm_varid, 'coordinates', 'x y'))
         !
      endif
      !
   endif
   !
   ! Add for final output:
   !
   NF90(nf90_def_var(map_file%ncid, 'total_runtime', NF90_FLOAT, (/map_file%runtime_dimid/),map_file%total_runtime_varid))
   NF90(nf90_put_att(map_file%ncid, map_file%total_runtime_varid, 'units', 's'))   
   NF90(nf90_put_att(map_file%ncid, map_file%total_runtime_varid, 'long_name', 'total_model_runtime_in_seconds'))
   !
   NF90(nf90_def_var(map_file%ncid, 'average_dt', NF90_FLOAT, (/map_file%runtime_dimid/), map_file%average_dt_varid))
   NF90(nf90_put_att(map_file%ncid, map_file%total_runtime_varid, 'units', 's'))   
   NF90(nf90_put_att(map_file%ncid, map_file%total_runtime_varid, 'long_name', 'model_average_timestep_in_seconds'))   
   ! 
   ! Finish definitions
   NF90(nf90_enddef(map_file%ncid))
   ! 
   ! Write grid to file
   !
   allocate(xz(mmax, nmax))
   allocate(yz(mmax, nmax))
   allocate(xg(mmax + 1, nmax + 1))
   allocate(yg(mmax + 1, nmax + 1))
   !
   do n = 1, nmax
      do m = 1, mmax
         xz(m, n) = x0 + cosrot*(1.0*(m - 0.5))*dx - sinrot*(1.0*(n - 0.5))*dy
         yz(m, n) = y0 + sinrot*(1.0*(m - 0.5))*dx + cosrot*(1.0*(n - 0.5))*dy
      enddo
   enddo   
   !
   do n = 1, nmax + 1
      do m = 1, mmax + 1
         xg(m, n) = x0 + cosrot*(1.0*(m - 1))*dx - sinrot*(1.0*(n - 1))*dy
         yg(m, n) = y0 + sinrot*(1.0*(m - 1))*dx + cosrot*(1.0*(n - 1))*dy
      enddo
   enddo   
   !
   NF90(nf90_put_var(map_file%ncid, map_file%face_x_varid, xz, (/1, 1/))) ! write xz of faces 
   !
   NF90(nf90_put_var(map_file%ncid, map_file%face_y_varid, yz, (/1, 1/))) ! write yz of faces   
   ! 
   ! now for cell corners
   NF90(nf90_put_var(map_file%ncid, map_file%corner_x_varid, xg, (/1, 1/))) ! write xg of corners
   !
   NF90(nf90_put_var(map_file%ncid, map_file%corner_y_varid, yg, (/1, 1/))) ! write yg of corners   
   !
   ! Write epsg, msk & bed level already to file
   !
   NF90(nf90_put_var(map_file%ncid, map_file%crs_varid, epsg))
   !   
   allocate(zsg(mmax, nmax))
   !
   zsg = FILL_VALUE       
   !
   do nm = 1, np
      !
      n    = z_index_z_n(nm)
      m    = z_index_z_m(nm)
      !      
      if (subgrid) then
         zsg(m, n) = subgrid_z_zmin(nm)
      else   
         zsg(m, n) = zb(nm)
      endif   
      !
   enddo
   !
   NF90(nf90_put_var(map_file%ncid, map_file%zb_varid, zsg, (/1, 1/))) ! write zb
   !
   zsg = 0 ! initialise as inactive points       
   !
   do nm = 1, np
      !
      n    = z_index_z_n(nm)
      m    = z_index_z_m(nm)
      !      
      zsg(m, n) = kcs(nm)
      !
   enddo
   !
   NF90(nf90_put_var(map_file%ncid, map_file%msk_varid, zsg, (/1, 1/))) ! write msk     
   !   
   ! Write infiltration map
   !
   if (infiltration) then
      !
      zsg = FILL_VALUE       
      !
      do nm = 1, np
         !
         n    = z_index_z_n(nm)
         m    = z_index_z_m(nm)
         !      
         if (inftype == 'con' .or. inftype == 'c2d') then
            zsg(m, n) = qinffield(nm)*3600*1000
         else
            zsg(m, n) = qinffield(nm)
         endif
         ! 
      enddo
      !
      NF90(nf90_put_var(map_file%ncid, map_file%qinf_varid, zsg, (/1, 1/))) ! write infiltration map
      !
   endif
   !
   ! write away intermediate data
   !
   NF90(nf90_sync(map_file%ncid)) !write away intermediate data
   !   
   end subroutine

   
   subroutine ncoutput_quadtree_map_init()
   !
   ! 1. Initialise dimensions/variables/attributes
   ! 2. write grid/msk/zb to file
   !
   use sfincs_date
   use sfincs_data   
   use quadtree   
   !
   implicit none   
   !   
   integer    :: nm, nmq, n, m, nn, ntmx, n_nodes, n_faces, ip, iref
   integer*8  :: two
   real*4     :: dxdy
   !
   real*4,    dimension(:),   allocatable :: nodes_x
   real*4,    dimension(:),   allocatable :: nodes_y
   integer*4, dimension(:,:), allocatable :: face_nodes
   real*4,    dimension(:),   allocatable :: vtmp
   integer*4, dimension(:),   allocatable :: vtmpi
   !
   ! Very lazy for now
   !
   n_faces = quadtree_nr_points
   n_nodes = quadtree_nr_points*4
   !
   nodes_x = 0.0
   nodes_y = 0.0
   face_nodes = 0
   !
   two = 2
   !
   allocate(nodes_x(n_nodes))
   allocate(nodes_y(n_nodes))
   allocate(face_nodes(4, n_faces))
   allocate(vtmp(n_faces))
   allocate(vtmpi(n_faces))
   !
   nn = 0
   !
   do nmq = 1, quadtree_nr_points
      !
      n = quadtree_n(nmq)
      m = quadtree_m(nmq)
         !
         iref = quadtree_level(nmq)
         dxdy  = quadtree_dxr(iref)
         !         
         nn = nn + 1
         !
         nodes_x(nn) = x0 + cosrot*(m - 1)*dxdy - sinrot*(n - 1)*dxdy
         nodes_y(nn) = y0 + sinrot*(m - 1)*dxdy + cosrot*(n - 1)*dxdy
         face_nodes(1, nmq) = nn
         !         
         nn = nn + 1
         !
         nodes_x(nn) = x0 + cosrot*(m    )*dxdy - sinrot*(n - 1)*dxdy
         nodes_y(nn) = y0 + sinrot*(m    )*dxdy + cosrot*(n - 1)*dxdy
         face_nodes(2, nmq) = nn
         !         
         nn = nn + 1
         !
         nodes_x(nn) = x0 + cosrot*(m    )*dxdy - sinrot*(n    )*dxdy
         nodes_y(nn) = y0 + sinrot*(m    )*dxdy + cosrot*(n    )*dxdy
         face_nodes(3, nmq) = nn
         !         
         nn = nn + 1
         !
         nodes_x(nn) = x0 + cosrot*(m - 1)*dxdy - sinrot*(n    )*dxdy
         nodes_y(nn) = y0 + sinrot*(m - 1)*dxdy + cosrot*(n    )*dxdy
         face_nodes(4, nmq) = nn
         !
!      endif
   enddo   
   !  
   NF90(nf90_create('sfincs_map.nc', ior(NF90_CLOBBER, NF90_64BIT_OFFSET), map_file%ncid))
   !
   ! Create dimensions
   ! grid, time, points
   ! do mmax/nmax-2 to not write away dummy cells
   !
   NF90(nf90_def_dim(map_file%ncid, 'nmesh2d_node', n_nodes, map_file%nmesh2d_node_dimid))
   NF90(nf90_def_dim(map_file%ncid, 'nmesh2d_face', n_faces, map_file%nmesh2d_face_dimid))
   NF90(nf90_def_dim(map_file%ncid, 'max_nmesh2d_face_nodes', 4, map_file%max_nmesh2d_face_nodes_dimid))
   !
   ! Time
   !
   NF90(nf90_def_dim(map_file%ncid, 'time', NF90_UNLIMITED, map_file%time_dimid)) ! time
   ntmx = max(int((t1out - t0out)/dtmaxout), 1)
   NF90(nf90_def_dim(map_file%ncid, 'timemax', ntmx, map_file%timemax_dimid)) ! time
   NF90(nf90_def_dim(map_file%ncid, 'runtime', 1, map_file%runtime_dimid)) ! total_runtime, average_dt       
   !
   ! Some metadata attributes 
   !
   NF90(nf90_put_att(map_file%ncid,nf90_global, "Conventions", "Conventions = 'CF-1.8 UGRID-1.0 Deltares-0.10'")) 
   NF90(nf90_put_att(map_file%ncid,nf90_global, "Build-Revision-Date-Netcdf-library", trim(nf90_inq_libvers()))) ! version of netcdf library
   NF90(nf90_put_att(map_file%ncid,nf90_global, "Producer", "SFINCS model: Super-Fast INundation of CoastS"))
   NF90(nf90_put_att(map_file%ncid,nf90_global, "Build-Revision", trim(build_revision))) 
   NF90(nf90_put_att(map_file%ncid,nf90_global, "Build-Date", trim(build_date)))
   NF90(nf90_put_att(map_file%ncid,nf90_global, "title", "SFINCS map netcdf output"))   
   !
   ! Add input params for reproducability
   !
   call ncoutput_add_params(map_file%ncid,map_file%inp_varid)   
   !
   ! Create variables
   !
   ! Domain
   !
   NF90(nf90_def_var(map_file%ncid, 'mesh2d', NF90_INT, (/map_file%runtime_dimid/), map_file%mesh2d_varid)) ! location of zb, zs etc. in cell centre
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'cf_role', 'mesh_topology'))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'long_name', 'Topology data of 2D network'))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'topology_dimension', 2))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'node_coordinates', 'mesh2d_node_x mesh2d_node_y'))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'node_dimension', 'nmesh2d_node'))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'max_face_nodes_dimension', 'max_nmesh2d_face_nodes'))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'face_node_connectivity', 'mesh2d_face_nodes'))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'face_dimension', 'nmesh2d_face'))
   !
   if (crsgeo) then
      !
      NF90(nf90_def_var(map_file%ncid, 'mesh2d_node_x', NF90_FLOAT, (/map_file%nmesh2d_node_dimid/), map_file%mesh2d_node_x_varid)) ! location of zb, zs etc. in cell centre
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'units', 'degrees'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'standard_name', 'longitude'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'long_name', 'longitude'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'mesh', 'mesh2d'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'location', 'node'))
      !
      NF90(nf90_def_var(map_file%ncid, 'mesh2d_node_y', NF90_FLOAT, (/map_file%nmesh2d_node_dimid/), map_file%mesh2d_node_y_varid)) ! location of zb, zs etc. in cell centre
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'units', 'degrees'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'standard_name', 'latitude'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'long_name', 'latitude'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'mesh', 'mesh2d'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'location', 'node'))
      !
   else
      !
      NF90(nf90_def_var(map_file%ncid, 'mesh2d_node_x', NF90_FLOAT, (/map_file%nmesh2d_node_dimid/), map_file%mesh2d_node_x_varid)) ! location of zb, zs etc. in cell centre
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'standard_name', 'projection_x_coordinate'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'long_name', 'x-coordinate of mesh nodes'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'mesh', 'mesh2d'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'location', 'node'))
      !
      NF90(nf90_def_var(map_file%ncid, 'mesh2d_node_y', NF90_FLOAT, (/map_file%nmesh2d_node_dimid/), map_file%mesh2d_node_y_varid)) ! location of zb, zs etc. in cell centre
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'standard_name', 'projection_y_coordinate'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'long_name', 'y-coordinate of mesh nodes'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'mesh', 'mesh2d'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'location', 'node'))
      !
   endif
   !
   NF90(nf90_def_var(map_file%ncid, 'mesh2d_face_nodes', NF90_INT, (/map_file%max_nmesh2d_face_nodes_dimid, map_file%nmesh2d_face_dimid/), map_file%mesh2d_face_nodes_varid)) ! location of zb, zs etc. in cell centre
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'cf_role', 'face_node_connectivity'))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'mesh', 'mesh2d'))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'location', 'face'))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'long_name', 'Mapping from every face to its corner nodes (counterclockwise)'))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'start_index', 1))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, '_FillValue', -999))
   !
   NF90(nf90_def_var(map_file%ncid, 'crs', NF90_INT, map_file%crs_varid)) ! For EPSG code
   NF90(nf90_put_att(map_file%ncid, map_file%crs_varid, 'EPSG', '-'))
   !
   NF90(nf90_def_var(map_file%ncid, 'zb', NF90_FLOAT, (/map_file%nmesh2d_face_dimid/), map_file%zb_varid)) ! bed level in cell centre
   NF90(nf90_put_att(map_file%ncid, map_file%zb_varid, '_FillValue', FILL_VALUE))   
   NF90(nf90_put_att(map_file%ncid, map_file%zb_varid, 'units', 'm'))
   NF90(nf90_put_att(map_file%ncid, map_file%zb_varid, 'standard_name', 'altitude'))
   NF90(nf90_put_att(map_file%ncid, map_file%zb_varid, 'long_name', 'bed_level_above_reference_level'))   
   !
   NF90(nf90_def_var(map_file%ncid, 'msk', NF90_INT, (/map_file%nmesh2d_face_dimid/), map_file%msk_varid)) ! input msk value in cell centre
   NF90(nf90_put_att(map_file%ncid, map_file%msk_varid, '_FillValue', -999))      
   NF90(nf90_put_att(map_file%ncid, map_file%msk_varid, 'units', '-'))
   NF90(nf90_put_att(map_file%ncid, map_file%msk_varid, 'standard_name', 'mask'))
   NF90(nf90_put_att(map_file%ncid, map_file%msk_varid, 'long_name', 'msk_active_cells')) 
   NF90(nf90_put_att(map_file%ncid, map_file%msk_varid, 'description', 'inactive=0, active=1, normal_boundary=2, outflow_boundary=3, wavemaker=4'))    
   !
   ! Time variables   
   !
   trefstr_iso8601 = date_to_iso8601(trefstr)
   !
   NF90(nf90_def_var(map_file%ncid, 'time', NF90_FLOAT, (/map_file%time_dimid/), map_file%time_varid)) ! time
   NF90(nf90_put_att(map_file%ncid, map_file%time_varid, 'units', 'seconds since ' // trim(trefstr_iso8601) ))  ! time stamp following ISO 8601
   NF90(nf90_put_att(map_file%ncid, map_file%time_varid, 'standard_name', 'time'))     
   NF90(nf90_put_att(map_file%ncid, map_file%time_varid, 'long_name', 'time_in_seconds_since_' // trim(trefstr_iso8601) ))  
   !
   ! Time varying map output
   NF90(nf90_def_var(map_file%ncid, 'zs', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%time_dimid/), map_file%zs_varid)) ! time-varying water level map
   NF90(nf90_put_att(map_file%ncid, map_file%zs_varid, '_FillValue', FILL_VALUE))
   NF90(nf90_put_att(map_file%ncid, map_file%zs_varid, 'units', 'm'))
   NF90(nf90_put_att(map_file%ncid, map_file%zs_varid, 'standard_name', 'sea_surface_height_above_mean_sea_level')) 
   NF90(nf90_put_att(map_file%ncid, map_file%zs_varid, 'long_name', 'water_level'))  
   !
   if (store_velocity) then
      !
      NF90(nf90_def_var(map_file%ncid, 'u', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%time_dimid/), map_file%u_varid)) ! time-varying u map 
      NF90(nf90_put_att(map_file%ncid, map_file%u_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(map_file%ncid, map_file%u_varid, 'units', 'm s-1'))
      NF90(nf90_put_att(map_file%ncid, map_file%u_varid, 'standard_name', 'sea_water_x_velocity')) ! not truly eastward when rotated, eastward_sea_water_velocity
      NF90(nf90_put_att(map_file%ncid, map_file%u_varid, 'long_name', 'flow_velocity_x_direction_in_cell_edge'))     
      !
      NF90(nf90_def_var(map_file%ncid, 'v', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%time_dimid/), map_file%v_varid)) ! time-varying u map 
      NF90(nf90_put_att(map_file%ncid, map_file%v_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(map_file%ncid, map_file%v_varid, 'units', 'm s-1'))
      NF90(nf90_put_att(map_file%ncid, map_file%v_varid, 'standard_name', 'sea_water_y_velocity')) ! not truly eastward when rotated, eastward_sea_water_velocity
      NF90(nf90_put_att(map_file%ncid, map_file%v_varid, 'long_name', 'flow_velocity_y_direction_in_cell_edge'))     
   endif   
   !
   ! Time varying spatial output
   if (store_maximum_waterlevel) then
      NF90(nf90_def_var(map_file%ncid, 'timemax', NF90_FLOAT, (/map_file%timemax_dimid/), map_file%timemax_varid)) ! time
      NF90(nf90_put_att(map_file%ncid, map_file%timemax_varid, 'units', 'seconds since ' // trim(trefstr_iso8601) ))  ! time stamp following ISO 8601
      NF90(nf90_put_att(map_file%ncid, map_file%timemax_varid, 'standard_name', 'time'))     
      NF90(nf90_put_att(map_file%ncid, map_file%timemax_varid, 'long_name', 'time_in_seconds_since_' // trim(trefstr_iso8601) ))  
   endif
   !
   if (store_maximum_waterlevel) then
      NF90(nf90_def_var(map_file%ncid, 'zsmax', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%timemax_dimid/), map_file%zsmax_varid)) ! time-varying maximum water level map
      NF90(nf90_put_att(map_file%ncid, map_file%zsmax_varid, '_FillValue', FILL_VALUE))
      NF90(nf90_put_att(map_file%ncid, map_file%zsmax_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%zsmax_varid, 'standard_name', 'maximum of sea_surface_height_above_mean_sea_level')) 
      NF90(nf90_put_att(map_file%ncid, map_file%zsmax_varid, 'long_name', 'maximum_water_level'))        
   endif
   !
   if (store_twet) then
      NF90(nf90_def_var(map_file%ncid, 'tmax', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%timemax_dimid/), map_file%tmax_varid)) ! time-varying duration wet cell
      NF90(nf90_put_att(map_file%ncid, map_file%tmax_varid, '_FillValue', FILL_VALUE))
      NF90(nf90_put_att(map_file%ncid, map_file%tmax_varid, 'units', 'seconds'))
      NF90(nf90_put_att(map_file%ncid, map_file%tmax_varid, 'standard_name', 'duration cell is considered wet')) 
      NF90(nf90_put_att(map_file%ncid, map_file%tmax_varid, 'long_name', 'duration_wet_cell'))   
      NF90(nf90_put_att(map_file%ncid, map_file%tmax_varid, 'cell_methods', 'time: sum'))
   endif
   !
   if (store_maximum_velocity) then
      NF90(nf90_def_var(map_file%ncid, 'vmax', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%timemax_dimid/), map_file%vmax_varid)) ! maximum flow velocity map
      NF90(nf90_put_att(map_file%ncid, map_file%vmax_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(map_file%ncid, map_file%vmax_varid, 'units', 'm s-1'))
      NF90(nf90_put_att(map_file%ncid, map_file%vmax_varid, 'standard_name', 'maximum_flow_velocity')) ! no standard name available
      NF90(nf90_put_att(map_file%ncid, map_file%vmax_varid, 'long_name', 'maximum_flow_velocity')) 
      NF90(nf90_put_att(map_file%ncid, map_file%vmax_varid, 'cell_methods', 'time: maximum'))
   endif
   !
   ! Store cumulative rainfall 
   !
   if (store_cumulative_precipitation) then  
      !   
      ! Cumulative precipitation
      !
      NF90(nf90_def_var(map_file%ncid, 'cumprcp', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%timemax_dimid/), map_file%cumprcp_varid)) ! cumulative precipitation map
      NF90(nf90_put_att(map_file%ncid, map_file%cumprcp_varid, '_FillValue', FILL_VALUE))          
      NF90(nf90_put_att(map_file%ncid, map_file%cumprcp_varid, 'units', 'mm'))
      NF90(nf90_put_att(map_file%ncid, map_file%cumprcp_varid, 'long_name', 'cumulative_precipitation_depth')) 
      NF90(nf90_put_att(map_file%ncid, map_file%cumprcp_varid, 'cell_methods', 'time: sum'))       
      !   
      ! Cumulative infiltration
      !
      NF90(nf90_def_var(map_file%ncid, 'cuminf', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%timemax_dimid/), map_file%cuminf_varid)) ! cumulative infiltration map
      NF90(nf90_put_att(map_file%ncid, map_file%cuminf_varid, '_FillValue', FILL_VALUE))          
      NF90(nf90_put_att(map_file%ncid, map_file%cuminf_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%cuminf_varid, 'long_name', 'cumulative_infiltration_depth')) 
      NF90(nf90_put_att(map_file%ncid, map_file%cuminf_varid, 'cell_methods', 'time: sum'))     
      !
   endif
   !
   ! Store maximum wind speed
   !
   if (wind .and. store_wind_max .and. meteo3d) then 
      NF90(nf90_def_var(map_file%ncid, 'windmax', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%timemax_dimid/), map_file%windmax_varid)) ! maximum wind speed m/s
      NF90(nf90_put_att(map_file%ncid, map_file%windmax_varid, '_FillValue', FILL_VALUE))
      NF90(nf90_put_att(map_file%ncid, map_file%windmax_varid, 'units', 'm s-1'))
      NF90(nf90_put_att(map_file%ncid, map_file%windmax_varid, 'long_name', 'maximum_wind_speed')) 
      NF90(nf90_put_att(map_file%ncid, map_file%windmax_varid, 'cell_methods', 'time: maximum'))         
   endif
   !
   if (snapwave) then  
      !
      NF90(nf90_def_var(map_file%ncid, 'hm0', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%time_dimid/), map_file%hm0_varid))
      NF90(nf90_put_att(map_file%ncid, map_file%hm0_varid, '_FillValue', FILL_VALUE))          
      NF90(nf90_put_att(map_file%ncid, map_file%hm0_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%hm0_varid, 'standard_name', 'hm0_wave_height'))
      NF90(nf90_put_att(map_file%ncid, map_file%hm0_varid, 'long_name', 'Hm0 wave height')) 
      !
      NF90(nf90_def_var(map_file%ncid, 'hm0ig', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%time_dimid/), map_file%hm0ig_varid))
      NF90(nf90_put_att(map_file%ncid, map_file%hm0ig_varid, '_FillValue', FILL_VALUE))          
      NF90(nf90_put_att(map_file%ncid, map_file%hm0ig_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%hm0ig_varid, 'standard_name', 'hm0_ig_wave_height'))
      NF90(nf90_put_att(map_file%ncid, map_file%hm0ig_varid, 'long_name', 'Hm0 infragravity wave height')) 
      !
      if (store_wave_forces) then
         !
         NF90(nf90_def_var(map_file%ncid, 'fwx', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%time_dimid/), map_file%fwx_varid))
         NF90(nf90_put_att(map_file%ncid, map_file%fwx_varid, '_FillValue', FILL_VALUE))          
         NF90(nf90_put_att(map_file%ncid, map_file%fwx_varid, 'units', 'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%fwx_varid, 'standard_name', 'wave_force_x'))
         NF90(nf90_put_att(map_file%ncid, map_file%fwx_varid, 'long_name', 'Wave force in x-direction')) 
         !
         NF90(nf90_def_var(map_file%ncid, 'fwy', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%time_dimid/), map_file%fwy_varid))
         NF90(nf90_put_att(map_file%ncid, map_file%fwy_varid, '_FillValue', FILL_VALUE))          
         NF90(nf90_put_att(map_file%ncid, map_file%fwy_varid, 'units', 'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%fwy_varid, 'standard_name', 'wave_force_y'))
         NF90(nf90_put_att(map_file%ncid, map_file%fwy_varid, 'long_name', 'Wave force in y-direction')) 
         !  
      endif
      !
      if (wavemaker) then
         !
         NF90(nf90_def_var(map_file%ncid, 'zsm', NF90_FLOAT, (/map_file%nmesh2d_face_dimid, map_file%time_dimid/), map_file%zsm_varid))
         NF90(nf90_put_att(map_file%ncid, map_file%zsm_varid, '_FillValue', FILL_VALUE))          
         NF90(nf90_put_att(map_file%ncid, map_file%zsm_varid, 'units', 'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%zsm_varid, 'standard_name', 'filtered_water_level'))
         NF90(nf90_put_att(map_file%ncid, map_file%zsm_varid, 'long_name', 'Filtered water level')) 
         !
      endif   
      !
   endif
   !
   if (infiltration) then
       NF90(nf90_def_var(map_file%ncid, 'qinf', NF90_FLOAT, (/map_file%nmesh2d_face_dimid/), map_file%qinf_varid)) 
       NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, '_FillValue', FILL_VALUE))     
       if (inftype == 'cna') then
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'standard_name', 'S')) 
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'long_name', 'moisture storage (S) capacity as computed from the curve number')) 
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'units', 'm'))
       elseif (inftype == 'cnb') then
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'standard_name', 'Smax')) 
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'long_name', 'maximum moisture storage (Smax) capacity as computed from the curve number')) 
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'units', 'm'))
       else
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'standard_name', 'qinf')) 
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'long_name', 'infiltration rate - constant in time')) 
           NF90(nf90_put_att(map_file%ncid, map_file%qinf_varid, 'units', 'mm h-1'))
       endif
   endif
   !
   ! Add for final output:
   NF90(nf90_def_var(map_file%ncid, 'total_runtime', NF90_FLOAT, (/map_file%runtime_dimid/),map_file%total_runtime_varid))
   NF90(nf90_put_att(map_file%ncid, map_file%total_runtime_varid, 'units', 's'))   
   NF90(nf90_put_att(map_file%ncid, map_file%total_runtime_varid, 'long_name', 'total_model_runtime_in_seconds'))
   !
   NF90(nf90_def_var(map_file%ncid, 'average_dt', NF90_FLOAT, (/map_file%runtime_dimid/), map_file%average_dt_varid))
   NF90(nf90_put_att(map_file%ncid, map_file%total_runtime_varid, 'units', 's'))   
   NF90(nf90_put_att(map_file%ncid, map_file%total_runtime_varid, 'long_name', 'model_average_timestep_in_seconds'))   
   ! 
   ! Finish definitions
   NF90(nf90_enddef(map_file%ncid))
   !
   NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_node_x_varid, nodes_x)) ! write node x 
   !
   NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_node_y_varid, nodes_y)) ! write node y
   ! 
   NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_face_nodes_varid, face_nodes))
   !
!   ! now for cell edges
!   NF90(nf90_put_var(map_file%ncid, map_file%face_x_varid, xz(2:nmax+1-1, 2:mmax+1-1), (/1, 1/))) ! write xz of edges
!   !
!   NF90(nf90_put_var(map_file%ncid, map_file%face_y_varid, yz(2:nmax+1-1, 2:mmax+1-1), (/1, 1/))) ! write yz of edges   
   !
   ! Write epsg, msk & bed level already to file
   !
   NF90(nf90_put_var(map_file%ncid, map_file%crs_varid, epsg))
   !
   vtmp = FILL_VALUE
   !
   if (subgrid) then
      do nmq = 1, quadtree_nr_points
         nm = index_sfincs_in_quadtree(nmq)
         if (nm>0) then
            vtmp(nmq) = subgrid_z_zmin(nm)
         endif
      enddo 
      NF90(nf90_put_var(map_file%ncid, map_file%zb_varid, vtmp))
   else
      do nmq = 1, quadtree_nr_points
         nm = index_sfincs_in_quadtree(nmq)
         if (nm>0) then
            vtmp(nmq) = zb(nm)
         endif
      enddo 
      NF90(nf90_put_var(map_file%ncid, map_file%zb_varid, vtmp))
   endif
   !
   vtmpi = 0
   !
   do nmq = 1, quadtree_nr_points
      nm = index_sfincs_in_quadtree(nmq)
      if (nm>0) then
         vtmpi(nmq) = kcs(nm)
      endif
   enddo 
   !
   NF90(nf90_put_var(map_file%ncid, map_file%msk_varid, vtmpi)) ! write msk 
   !   
   ! Write infiltration map
   !   
   vtmp = FILL_VALUE
   !
   if (infiltration) then
      !
      if (inftype == 'con' .or. inftype == 'c2d') then
         do nmq = 1, quadtree_nr_points
            nm = index_sfincs_in_quadtree(nmq)
            if (nm>0) then
               vtmp(nmq) = qinfmap(nm)*3600*1000
            endif
         enddo 
      else
         do nmq = 1, quadtree_nr_points
            nm = index_sfincs_in_quadtree(nmq)
            if (nm>0) then
               vtmp(nmq) = qinffield(nm)
            endif
         enddo 
      endif
      !
      NF90(nf90_put_var(map_file%ncid, map_file%qinf_varid, vtmp)) ! write infiltration map
      !
   endif
   !
   ! write away intermediate data
   !
   NF90(nf90_sync(map_file%ncid)) !write away intermediate data
   !
   end subroutine
   
   
   subroutine ncoutput_his_init()
   ! 1. Initialise dimensions/variables/attributes
   ! 2. write grid/msk/zb to file
   !
   use sfincs_date
   use sfincs_data   
   use sfincs_structures
   !
   implicit none   
   !
   integer                      :: nm, n, m, istruc, struc_nm, npars
   
   integer*4,          dimension(:),   allocatable :: index_u_m
   integer*4,          dimension(:),   allocatable :: index_u_n
   integer*4,          dimension(:),   allocatable :: index_v_m
   integer*4,          dimension(:),   allocatable :: index_v_n
     
   real*4, dimension(:,:), allocatable :: struc_info
   real*4, dimension(:), allocatable :: struc_x
   real*4, dimension(:), allocatable :: struc_y
   real*4, dimension(:), allocatable :: struc_height
   
   !
   if (nobs==0 .and. nrcrosssections==0) then ! If no observation points his file is not created        
        return
   endif
   !
   NF90(nf90_create('sfincs_his.nc', NF90_CLOBBER, his_file%ncid))
   !
   ! Create dimensions
   ! time, stations
   NF90(nf90_def_dim(his_file%ncid, 'time', NF90_UNLIMITED, his_file%time_dimid)) ! time   
   !
   if (nobs>0) then   
      NF90(nf90_def_dim(his_file%ncid, 'stations', nobs, his_file%points_dimid)) ! nr of observation points
   else
      NF90(nf90_def_dim(his_file%ncid, 'stations', 1, his_file%points_dimid)) ! easiest to initiate dimension as 1    
   endif
   !
   if (nrcrosssections>0) then   
      NF90(nf90_def_dim(his_file%ncid, 'crosssections', nrcrosssections, his_file%crosssections_dimid)) ! nr of crosssections
   endif
   !
   if (nrstructures>0) then   
      NF90(nf90_def_dim(his_file%ncid, 'structures', nrstructures, his_file%structures_dimid)) ! nr of structures
   endif   
   !
   NF90(nf90_def_dim(his_file%ncid, 'pointnamelength', 256, his_file%pointnamelength_dimid)) ! length of station_name per obs point  
   NF90(nf90_def_dim(his_file%ncid, 'runtime', 1, his_file%runtime_dimid)) ! total_runtime, average_dt    
   !
   !     
   ! Some metadata attributes 
   NF90(nf90_put_att(his_file%ncid,nf90_global, "Conventions", "Conventions = 'CF-1.6, SGRID-0.3")) 
   NF90(nf90_put_att(his_file%ncid,nf90_global, "Build-Revision-Date-Netcdf-library", trim(nf90_inq_libvers()))) ! version of netcdf library
   NF90(nf90_put_att(his_file%ncid,nf90_global, "Producer", "SFINCS model: Super-Fast INundation of CoastS"))
   NF90(nf90_put_att(his_file%ncid,nf90_global, "Build-Revision", trim(build_revision))) 
   NF90(nf90_put_att(his_file%ncid,nf90_global, "Build-Date", trim(build_date)))
   NF90(nf90_put_att(his_file%ncid,nf90_global, "title", "SFINCS his point netcdf output"))     
   !
   ! add input params for reproducability   
   call ncoutput_add_params(his_file%ncid,his_file%inp_varid)   
   !
   ! Create variables
   ! Point identifier
   NF90(nf90_def_var(his_file%ncid, 'station_id', NF90_FLOAT, (/his_file%points_dimid/), his_file%station_id_varid)) 
   !NF90(nf90_put_att(his_file%ncid, his_file%station_id_varid, 'units', '-')) !not wanted in fews
   !
   NF90(nf90_def_var(his_file%ncid, 'station_name', NF90_CHAR, (/his_file%pointnamelength_dimid, his_file%points_dimid/), his_file%station_name_varid))
   !NF90(nf90_put_att(his_file%ncid, his_file%station_name_varid, 'units', '-')) !not wanted in fews   
   !
   if (nrcrosssections>0) then
      NF90(nf90_def_var(his_file%ncid, 'crosssection_name', NF90_CHAR, (/his_file%pointnamelength_dimid, his_file%crosssections_dimid/), his_file%crosssection_name_varid))
   endif   
   !
   !NF90(nf90_put_att(his_file%ncid, his_file%station_name_varid, 'units', '-')) !not wanted in fews   
   !
   ! Domain
   NF90(nf90_def_var(his_file%ncid, 'station_x', NF90_FLOAT, (/his_file%points_dimid/), his_file%station_x_varid))   ! non snapped input coordinate 
   NF90(nf90_put_att(his_file%ncid, his_file%station_x_varid, 'units', 'm'))
   NF90(nf90_put_att(his_file%ncid, his_file%station_x_varid, 'standard_name', 'projection_x_coordinate'))
   NF90(nf90_put_att(his_file%ncid, his_file%station_x_varid, 'long_name', 'original_x_coordinate_of_station'))
   NF90(nf90_put_att(his_file%ncid, his_file%station_x_varid, 'grid_mapping', 'crs'))   
   !NF90(nf90_put_att(his_file%ncid, his_file%station_x_varid, 'grid', 'sfincsgrid'))   !keep this?
   !
   NF90(nf90_def_var(his_file%ncid, 'station_y', NF90_FLOAT, (/his_file%points_dimid/), his_file%station_y_varid)) 
   NF90(nf90_put_att(his_file%ncid, his_file%station_y_varid, 'units', 'm'))
   NF90(nf90_put_att(his_file%ncid, his_file%station_y_varid, 'standard_name', 'projection_y_coordinate'))
   NF90(nf90_put_att(his_file%ncid, his_file%station_y_varid, 'long_name', 'original_y_coordinate_of_station'))
   NF90(nf90_put_att(his_file%ncid, his_file%station_y_varid, 'grid_mapping', 'crs'))   
   !
   NF90(nf90_def_var(his_file%ncid, 'point_x', NF90_FLOAT, (/his_file%points_dimid/), his_file%point_x_varid))  ! snapped coordinate as used in sfincs
   NF90(nf90_put_att(his_file%ncid, his_file%point_x_varid, 'units', 'm'))
   NF90(nf90_put_att(his_file%ncid, his_file%point_x_varid, 'standard_name', 'projection_x_coordinate'))
   NF90(nf90_put_att(his_file%ncid, his_file%point_x_varid, 'long_name', 'point_x'))    
   NF90(nf90_put_att(his_file%ncid, his_file%point_x_varid, 'grid_mapping', 'crs'))   
   !NF90(nf90_put_att(his_file%ncid, his_file%point_x_varid, 'grid', 'sfincsgrid'))   !keep this?
   !
   NF90(nf90_def_var(his_file%ncid, 'point_y', NF90_FLOAT, (/his_file%points_dimid/), his_file%point_y_varid)) 
   NF90(nf90_put_att(his_file%ncid, his_file%point_y_varid, 'units', 'm'))
   NF90(nf90_put_att(his_file%ncid, his_file%point_y_varid, 'standard_name', 'projection_y_coordinate'))
   NF90(nf90_put_att(his_file%ncid, his_file%point_y_varid, 'long_name', 'point_y'))
   NF90(nf90_put_att(his_file%ncid, his_file%point_y_varid, 'grid_mapping', 'crs'))   
   !NF90(nf90_put_att(his_file%ncid, his_file%point_y_varid, 'grid', 'sfincsgrid'))   !keep this?   
   !
   NF90(nf90_def_var(his_file%ncid, 'crs', NF90_INT, his_file%crs_varid)) ! For EPSG code
   NF90(nf90_put_att(his_file%ncid, his_file%crs_varid, 'EPSG', '-'))     
   NF90(nf90_put_att(his_file%ncid, his_file%crs_varid, 'epsg_code', 'EPSG:' // trim(epsg_code) ))   !--> add epsg_code like FEWS wants   
   !
   NF90(nf90_def_var(his_file%ncid, 'point_zb', NF90_FLOAT, (/his_file%points_dimid/), his_file%zb_varid)) ! bed level in cell centre, for points
   NF90(nf90_put_att(his_file%ncid, his_file%zb_varid, '_FillValue', FILL_VALUE))   
   NF90(nf90_put_att(his_file%ncid, his_file%zb_varid, 'units', 'm'))
   NF90(nf90_put_att(his_file%ncid, his_file%zb_varid, 'standard_name', 'altitude'))
   NF90(nf90_put_att(his_file%ncid, his_file%zb_varid, 'long_name', 'bed_level_above_reference_level'))
   NF90(nf90_put_att(his_file%ncid, his_file%zb_varid, 'coordinates', 'station_id station_name point_x point_y'))
   !
   if (nrstructures>0) then
      !
      NF90(nf90_def_var(his_file%ncid, 'structure_x', NF90_FLOAT, (/his_file%structures_dimid/), his_file%structure_x_varid))  ! snapped coordinate as used in sfincs
      NF90(nf90_put_att(his_file%ncid, his_file%structure_x_varid, 'units', 'm'))
      NF90(nf90_put_att(his_file%ncid, his_file%structure_x_varid, 'standard_name', 'projection_x_coordinate'))
      NF90(nf90_put_att(his_file%ncid, his_file%structure_x_varid, 'long_name', 'structure_x'))    
      NF90(nf90_put_att(his_file%ncid, his_file%structure_x_varid, 'grid_mapping', 'crs'))       
      !
      NF90(nf90_def_var(his_file%ncid, 'structure_y', NF90_FLOAT, (/his_file%structures_dimid/), his_file%structure_y_varid))  ! snapped coordinate as used in sfincs
      NF90(nf90_put_att(his_file%ncid, his_file%structure_y_varid, 'units', 'm'))
      NF90(nf90_put_att(his_file%ncid, his_file%structure_y_varid, 'standard_name', 'projection_y_coordinate'))
      NF90(nf90_put_att(his_file%ncid, his_file%structure_y_varid, 'long_name', 'structure_y'))    
      NF90(nf90_put_att(his_file%ncid, his_file%structure_y_varid, 'grid_mapping', 'crs'))   
      !
      NF90(nf90_def_var(his_file%ncid, 'structure_height', NF90_FLOAT, (/his_file%structures_dimid/), his_file%structure_height_varid)) ! structure height 
      NF90(nf90_put_att(his_file%ncid, his_file%structure_height_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(his_file%ncid, his_file%structure_height_varid, 'units', 'm'))
      NF90(nf90_put_att(his_file%ncid, his_file%structure_height_varid, 'standard_name', 'altitude'))
      NF90(nf90_put_att(his_file%ncid, his_file%structure_height_varid, 'long_name', 'interpolated_structure_height_above_reference_level'))      
      !
   endif         
   !
   ! Time variables 
   trefstr_iso8601 = date_to_iso8601(trefstr)
   !  
   NF90(nf90_def_var(his_file%ncid, 'time', NF90_FLOAT, (/his_file%time_dimid/), his_file%time_varid)) ! time
   NF90(nf90_put_att(his_file%ncid, his_file%time_varid, 'units', 'seconds since ' // trim(trefstr_iso8601) ))  ! time stamp following ISO 8601
   NF90(nf90_put_att(his_file%ncid, his_file%time_varid, 'standard_name', 'time'))     
   NF90(nf90_put_att(his_file%ncid, his_file%time_varid, 'long_name', 'time_in_seconds_since_' // trim(trefstr_iso8601) ))   !--> add  trefstr from sfincs_input
   !
   ! Time varying map output
   !
   NF90(nf90_def_var(his_file%ncid, 'point_zs', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%zs_varid)) ! time-varying water level point
   NF90(nf90_put_att(his_file%ncid, his_file%zs_varid, '_FillValue', FILL_VALUE))
   NF90(nf90_put_att(his_file%ncid, his_file%zs_varid, 'units', 'm'))
   NF90(nf90_put_att(his_file%ncid, his_file%zs_varid, 'standard_name', 'sea_surface_height_above_mean_sea_level'))
   NF90(nf90_put_att(his_file%ncid, his_file%zs_varid, 'long_name', 'water_level'))
   NF90(nf90_put_att(his_file%ncid, his_file%zs_varid, 'coordinates', 'station_id station_name point_x point_y'))
   !
   if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then
      NF90(nf90_def_var(his_file%ncid, 'point_h', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%h_varid)) ! time-varying water depth map
      NF90(nf90_put_att(his_file%ncid, his_file%h_varid, '_FillValue', FILL_VALUE))
      NF90(nf90_put_att(his_file%ncid, his_file%h_varid, 'units', 'm'))
      NF90(nf90_put_att(his_file%ncid, his_file%h_varid, 'standard_name', 'depth')) 
      NF90(nf90_put_att(his_file%ncid, his_file%h_varid, 'long_name', 'water_depth'))     
      NF90(nf90_put_att(his_file%ncid, his_file%h_varid, 'coordinates', 'station_id station_name point_x point_y'))
   endif
   !
   if (store_velocity) then
      !
      NF90(nf90_def_var(his_file%ncid, 'point_u', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%u_varid)) ! time-varying u point 
      NF90(nf90_put_att(his_file%ncid, his_file%u_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(his_file%ncid, his_file%u_varid, 'units', 'm s-1'))
      NF90(nf90_put_att(his_file%ncid, his_file%u_varid, 'standard_name', 'sea_water_x_velocity')) ! not truly eastward when rotated, eastward_sea_water_velocity
      NF90(nf90_put_att(his_file%ncid, his_file%u_varid, 'long_name', 'flow_velocity_x_direction'))     
      NF90(nf90_put_att(his_file%ncid, his_file%u_varid, 'coordinates', 'station_id station_name point_x point_y'))
      !
      NF90(nf90_def_var(his_file%ncid, 'point_v', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%v_varid)) ! time-varying u point 
      NF90(nf90_put_att(his_file%ncid, his_file%v_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(his_file%ncid, his_file%v_varid, 'units', 'm s-1'))
      NF90(nf90_put_att(his_file%ncid, his_file%v_varid, 'standard_name', 'sea_water_x_velocity')) ! not truly northward when rotated, northward_sea_water_velocity
      NF90(nf90_put_att(his_file%ncid, his_file%v_varid, 'long_name', 'flow_velocity_y_direction'))     
      NF90(nf90_put_att(his_file%ncid, his_file%v_varid, 'coordinates', 'station_id station_name point_x point_y'))
      !
      NF90(nf90_def_var(his_file%ncid, 'point_uvmag', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%uvmag_varid)) ! time-varying flow velocity magnitude 
      NF90(nf90_put_att(his_file%ncid, his_file%uvmag_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(his_file%ncid, his_file%uvmag_varid, 'units', 'm s-1'))
      NF90(nf90_put_att(his_file%ncid, his_file%uvmag_varid, 'standard_name', 'sea_water_velocity'))
      NF90(nf90_put_att(his_file%ncid, his_file%uvmag_varid, 'long_name', 'flow_velocity_magnitude'))     
      NF90(nf90_put_att(his_file%ncid, his_file%uvmag_varid, 'coordinates', 'station_id station_name point_x point_y'))
      !
      NF90(nf90_def_var(his_file%ncid, 'point_uvdir', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%uvdir_varid)) ! time-varying flow velocity direction 
      NF90(nf90_put_att(his_file%ncid, his_file%uvdir_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(his_file%ncid, his_file%uvdir_varid, 'units', 'degrees'))
      NF90(nf90_put_att(his_file%ncid, his_file%uvdir_varid, 'standard_name', 'sea_water_velocity_direction'))
      NF90(nf90_put_att(his_file%ncid, his_file%uvdir_varid, 'long_name', 'flow_velocity_direction'))     
      NF90(nf90_put_att(his_file%ncid, his_file%uvdir_varid, 'coordinates', 'station_id station_name point_x point_y'))
      !
   endif
   !
   ! Add infiltration
   !
   if (infiltration) then
      !
      NF90(nf90_def_var(his_file%ncid, 'point_qinf', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%qinf_varid)) ! time-varying infiltration point 
      NF90(nf90_put_att(his_file%ncid, his_file%qinf_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(his_file%ncid, his_file%qinf_varid, 'units', 'mm hr-1'))
      NF90(nf90_put_att(his_file%ncid, his_file%qinf_varid, 'long_name', 'infiltration_rate'))        
      NF90(nf90_put_att(his_file%ncid, his_file%qinf_varid, 'coordinates', 'station_id station_name point_x point_y'))
      !
   endif
   ! 
   ! More output for CN method with recovery
   !
   if (inftype == 'cnb') then
      !
      NF90(nf90_def_var(his_file%ncid, 'point_S', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%S_varid)) ! time-varying S
      NF90(nf90_put_att(his_file%ncid, his_file%S_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(his_file%ncid, his_file%S_varid, 'units', 'm'))
      NF90(nf90_put_att(his_file%ncid, his_file%S_varid, 'long_name', 'current moisture storage (Se) capacity')) 
      NF90(nf90_put_att(his_file%ncid, his_file%S_varid, 'coordinates', 'station_id station_name point_x point_y'))
      !
   endif
   !
   if (snapwave) then  
      !
      NF90(nf90_def_var(his_file%ncid, 'hm0', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%hm0_varid)) ! time-varying water level point
      NF90(nf90_put_att(his_file%ncid, his_file%hm0_varid, '_FillValue', FILL_VALUE))
      NF90(nf90_put_att(his_file%ncid, his_file%hm0_varid, 'units', 'm'))
      NF90(nf90_put_att(his_file%ncid, his_file%hm0_varid, 'standard_name', 'hm0_wave_height')) 
      NF90(nf90_put_att(his_file%ncid, his_file%hm0_varid, 'long_name', 'Hm0 wave height'))  
      NF90(nf90_put_att(his_file%ncid, his_file%hm0_varid, 'coordinates', 'station_id station_name point_x point_y'))
      !
      NF90(nf90_def_var(his_file%ncid, 'hm0ig', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%hm0ig_varid)) ! time-varying water level point
      NF90(nf90_put_att(his_file%ncid, his_file%hm0ig_varid, '_FillValue', FILL_VALUE))
      NF90(nf90_put_att(his_file%ncid, his_file%hm0ig_varid, 'units', 'm'))
      NF90(nf90_put_att(his_file%ncid, his_file%hm0ig_varid, 'standard_name', 'hm0_ig_wave_height')) 
      NF90(nf90_put_att(his_file%ncid, his_file%hm0ig_varid, 'long_name', 'Hm0 infragravity wave height'))  
      NF90(nf90_put_att(his_file%ncid, his_file%hm0ig_varid, 'coordinates', 'station_id station_name point_x point_y'))
      !
      NF90(nf90_def_var(his_file%ncid, 'tp', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%tp_varid)) ! time-varying water level point
      NF90(nf90_put_att(his_file%ncid, his_file%tp_varid, '_FillValue', FILL_VALUE))
      NF90(nf90_put_att(his_file%ncid, his_file%tp_varid, 'units', 's'))
      NF90(nf90_put_att(his_file%ncid, his_file%tp_varid, 'standard_name', 'peak_wave_period')) 
      NF90(nf90_put_att(his_file%ncid, his_file%tp_varid, 'long_name', 'Peak wave period'))  
      NF90(nf90_put_att(his_file%ncid, his_file%tp_varid, 'coordinates', 'station_id station_name point_x point_y'))
      !
      if (store_wave_direction) then
         !
         NF90(nf90_def_var(his_file%ncid, 'wavdir', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%wavdir_varid)) ! time-varying water level point
         NF90(nf90_put_att(his_file%ncid, his_file%wavdir_varid, '_FillValue', FILL_VALUE))
         NF90(nf90_put_att(his_file%ncid, his_file%wavdir_varid, 'units', 's'))
         NF90(nf90_put_att(his_file%ncid, his_file%wavdir_varid, 'standard_name', 'mean_wave_direction')) 
         NF90(nf90_put_att(his_file%ncid, his_file%wavdir_varid, 'long_name', 'Mean wave direction'))  
         NF90(nf90_put_att(his_file%ncid, his_file%wavdir_varid, 'coordinates', 'station_id station_name point_x point_y'))
         !
         NF90(nf90_def_var(his_file%ncid, 'dirspr', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%dirspr_varid)) ! time-varying water level point
         NF90(nf90_put_att(his_file%ncid, his_file%dirspr_varid, '_FillValue', FILL_VALUE))
         NF90(nf90_put_att(his_file%ncid, his_file%dirspr_varid, 'units', 's'))
         NF90(nf90_put_att(his_file%ncid, his_file%dirspr_varid, 'standard_name', 'wave_directional_spreading')) 
         NF90(nf90_put_att(his_file%ncid, his_file%dirspr_varid, 'long_name', 'Wave directional spreading'))  
         NF90(nf90_put_att(his_file%ncid, his_file%dirspr_varid, 'coordinates', 'station_id station_name point_x point_y'))
         !
      endif   
      !
      if (wavemaker) then
         !
         NF90(nf90_def_var(his_file%ncid, 'zsm', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%zsm_varid)) ! time-varying water level point
         NF90(nf90_put_att(his_file%ncid, his_file%zsm_varid, '_FillValue', FILL_VALUE))
         NF90(nf90_put_att(his_file%ncid, his_file%zsm_varid, 'units', 'm'))
         NF90(nf90_put_att(his_file%ncid, his_file%zsm_varid, 'standard_name', 'filtered_water_level')) 
         NF90(nf90_put_att(his_file%ncid, his_file%zsm_varid, 'long_name', 'Filtered water level'))  
         NF90(nf90_put_att(his_file%ncid, his_file%zsm_varid, 'coordinates', 'station_id station_name point_x point_y'))
         !
      endif   
   endif
   !
   if (store_meteo) then      
      !
      if (wind) then
         !
         NF90(nf90_def_var(his_file%ncid, 'point_wind_speed', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%wind_speed_varid)) ! time-varying patm point 
         NF90(nf90_put_att(his_file%ncid, his_file%wind_speed_varid, '_FillValue', FILL_VALUE))   
         NF90(nf90_put_att(his_file%ncid, his_file%wind_speed_varid, 'units', 'm s-1'))
         NF90(nf90_put_att(his_file%ncid, his_file%wind_speed_varid, 'long_name', 'wind_speed'))        
         NF90(nf90_put_att(his_file%ncid, his_file%wind_speed_varid, 'coordinates', 'station_id station_name point_x point_y'))
         !
         NF90(nf90_def_var(his_file%ncid, 'point_wind_direction', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%wind_dir_varid)) ! time-varying patm point 
         NF90(nf90_put_att(his_file%ncid, his_file%wind_dir_varid, '_FillValue', FILL_VALUE))   
         NF90(nf90_put_att(his_file%ncid, his_file%wind_dir_varid, 'units', 'degrees'))
         NF90(nf90_put_att(his_file%ncid, his_file%wind_dir_varid, 'long_name', 'wind_direction'))        
         NF90(nf90_put_att(his_file%ncid, his_file%wind_dir_varid, 'coordinates', 'station_id station_name point_x point_y'))
         !
      endif
      !
      if (patmos) then      
         !
         NF90(nf90_def_var(his_file%ncid, 'point_patm', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%patm_varid)) ! time-varying patm point 
         NF90(nf90_put_att(his_file%ncid, his_file%patm_varid, '_FillValue', FILL_VALUE))   
         NF90(nf90_put_att(his_file%ncid, his_file%patm_varid, 'units', 'Pa'))
         NF90(nf90_put_att(his_file%ncid, his_file%patm_varid, 'long_name', 'surface_air_pressure'))        
         NF90(nf90_put_att(his_file%ncid, his_file%patm_varid, 'coordinates', 'station_id station_name point_x point_y'))
         !
      endif
      !
      if (precip) then      
         NF90(nf90_def_var(his_file%ncid, 'point_prcp', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%prcp_varid)) ! time-varying prcp point 
         NF90(nf90_put_att(his_file%ncid, his_file%prcp_varid, '_FillValue', FILL_VALUE))   
         NF90(nf90_put_att(his_file%ncid, his_file%prcp_varid, 'units', 'mm hr-1'))
         NF90(nf90_put_att(his_file%ncid, his_file%prcp_varid, 'long_name', 'precipitation_rate'))        
         NF90(nf90_put_att(his_file%ncid, his_file%prcp_varid, 'coordinates', 'station_id station_name point_x point_y'))
      endif
      ! 
   endif   
   !   
   if (nrcrosssections>0) then
      !
      NF90(nf90_def_var(his_file%ncid, 'crosssection_discharge', NF90_FLOAT, (/his_file%crosssections_dimid, his_file%time_dimid/), his_file%discharge_varid)) ! time-varying crossection discharge 
      NF90(nf90_put_att(his_file%ncid, his_file%discharge_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(his_file%ncid, his_file%discharge_varid, 'units', 'm3 s-1'))
      NF90(nf90_put_att(his_file%ncid, his_file%discharge_varid, 'long_name', 'discharge'))
      NF90(nf90_put_att(his_file%ncid, his_file%discharge_varid, 'coordinates', 'crosssection_name'))
      !
   endif
   !
   ! Add for final output:
   NF90(nf90_def_var(his_file%ncid, 'total_runtime', NF90_FLOAT, (/his_file%runtime_dimid/), his_file%total_runtime_varid))
   NF90(nf90_put_att(his_file%ncid, his_file%total_runtime_varid, 'units', 's'))   
   NF90(nf90_put_att(his_file%ncid, his_file%total_runtime_varid, 'long_name', 'total_model_runtime_in_seconds'))
   !
   NF90(nf90_def_var(his_file%ncid, 'average_dt', NF90_FLOAT, (/his_file%runtime_dimid/), his_file%average_dt_varid))
   NF90(nf90_put_att(his_file%ncid, his_file%average_dt_varid, 'units', 's'))   
   NF90(nf90_put_att(his_file%ncid, his_file%average_dt_varid, 'long_name', 'model_average_timestep_in_seconds'))   
   !    
   ! Finish definitions
   NF90(nf90_enddef(his_file%ncid))
   !
   NF90(nf90_put_var(his_file%ncid, his_file%crs_varid, epsg))  ! write epsg
   !
   NF90(nf90_put_var(his_file%ncid, his_file%station_id_varid, idobs))  ! write station_id   
   !   
   NF90(nf90_put_var(his_file%ncid, his_file%station_name_varid, nameobs))  ! write station_name      ! , (/1, nobs/)
   !
   if (nrcrosssections>0) then
      NF90(nf90_put_var(his_file%ncid, his_file%crosssection_name_varid, namecrs))  ! write station_name      ! , (/1, nobs/)
   endif   
   !
   if (nrstructures>0) then
      !
      ! Allocate structure
      call give_structure_information(struc_info)
      allocate(struc_x(nrstructures))
      allocate(struc_y(nrstructures))
      allocate(struc_height(nrstructures))
      !
      do istruc = 1, nrstructures
         !
         struc_x(istruc)        = struc_info(istruc,1)
         struc_y(istruc)        = struc_info(istruc,2)
         struc_height(istruc)   = struc_info(istruc,3)
         !
      enddo
      !
      NF90(nf90_put_var(his_file%ncid,  his_file%structure_x_varid, struc_x)) ! write structure_x, input struc_x
      !
      NF90(nf90_put_var(his_file%ncid,  his_file%structure_y_varid, struc_y)) ! write structure_y, input struc_y
      !
      NF90(nf90_put_var(his_file%ncid,  his_file%structure_height_varid, struc_height)) ! write structure_height, input struc_height
      !      
   endif
   !   
   NF90(nf90_put_var(his_file%ncid, his_file%station_x_varid, xobs)) ! write station_x, input xobs
   !    
   NF90(nf90_put_var(his_file%ncid, his_file%station_y_varid, yobs)) ! write station_y, input yobs
   !   
   NF90(nf90_put_var(his_file%ncid, his_file%point_x_varid, xgobs)) ! write point_x, now actual value on grid is written rather than input xobs
   !    
   NF90(nf90_put_var(his_file%ncid, his_file%point_y_varid, ygobs)) ! write point_y, now actual value on grid is written rather than input yobs
   !  
   NF90(nf90_put_var(his_file%ncid, his_file%zb_varid, zbobs)) ! write point_zb
   !     
   NF90(nf90_sync(his_file%ncid)) !write away intermediate data
   !
   end subroutine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  
   subroutine ncoutput_update_regular_map(t,ntmapout)
   !
   ! Write time, zs, u, v  
   !
   use sfincs_data   
   !
   implicit none   
   !
   real*8                       :: t  
   real*4                       :: uz, vz
   !
   integer  :: ntmapout       
   !
   integer                      :: nm, n, m, nmd1, nmu1, ndm1, num1
   real*4, dimension(:,:), allocatable :: zsg
   real*4, dimension(:,:), allocatable :: zsgu
   real*4, dimension(:,:), allocatable :: zsgv
   !
   NF90(nf90_put_var(map_file%ncid, map_file%time_varid, t, (/ntmapout/))) ! write time
   !  
   allocate(zsg(mmax, nmax))
   !
   zsg = FILL_VALUE       ! set to fill value
   !
   do nm = 1, np
      !
      n    = z_index_z_n(nm)
      m    = z_index_z_m(nm)
      !      
      if (subgrid) then
         if ((zs(nm) - subgrid_z_zmin(nm)) > huthresh) then
            zsg(m, n) = zs(nm)
         endif
      else
         if ((zs(nm) - zb(nm)) > huthresh) then
            zsg(m, n) = zs(nm)
         endif
      endif
   enddo
   !
   NF90(nf90_put_var(map_file%ncid, map_file%zs_varid, zsg, (/1, 1, ntmapout/))) ! write zs
   ! 
   zsg = FILL_VALUE       
   !
   do nm = 1, np
      !
      n    = z_index_z_n(nm)
      m    = z_index_z_m(nm)
      !      
      if (subgrid) then
         zsg(m, n) = zs(nm) - subgrid_z_zmin(nm)
      else
         zsg(m, n) = zs(nm) - zb(nm)
      endif
   enddo
   !
   if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then   
      ! 
      NF90(nf90_put_var(map_file%ncid, map_file%h_varid, zsg, (/1, 1, ntmapout/))) ! write h
      !
   endif
   !            
   if (store_velocity) then
      !
      allocate(zsgu(mmax, nmax))
      allocate(zsgv(mmax, nmax))
      zsgu = FILL_VALUE
      zsgv = FILL_VALUE
      !
      do nm = 1, np
         !
         n    = z_index_z_n(nm)
         m    = z_index_z_m(nm)
         !
         nmd1 = z_index_uv_md1(nm)
         nmu1 = z_index_uv_mu1(nm)
         uz = 0.0
         if (nmd1>0) then
            uz = uz + 0.5*uv(nmd1)
         endif   
         if (nmu1>0) then
            uz = uz + 0.5*uv(nmu1)
         endif   
         !
         ndm1 = z_index_uv_nd1(nm)
         num1 = z_index_uv_nu1(nm)
         vz = 0.0
         if (ndm1>0) then
            vz = vz + 0.5*uv(ndm1)
         endif   
         if (num1>0) then
            vz = vz + 0.5*uv(num1)
         endif   
         !         
         zsgu(m, n) = cosrot*uz - sinrot*vz            
         zsgv(m, n) = sinrot*uz + cosrot*vz
         !
      enddo   
      !
      NF90(nf90_put_var(map_file%ncid, map_file%u_varid, zsgu, (/1, 1, ntmapout/)))
      !
      NF90(nf90_put_var(map_file%ncid, map_file%v_varid, zsgv, (/1, 1, ntmapout/)))
      !
      deallocate(zsgu)
      deallocate(zsgv)
      !
   endif
   !
   ! Store S_effective (only for CN method with recovery)
   !
   if (inftype == 'cnb') then
      !
      zsg = FILL_VALUE
      !
      do nm = 1, np
         !
         n    = z_index_z_n(nm)
         m    = z_index_z_m(nm)
         !      
         zsg(m, n) = qinffield2(nm)
         !
      enddo
      !
      NF90(nf90_put_var(map_file%ncid, map_file%Seff_varid, zsg, (/1, 1, ntmapout/)))
      !
   endif
   !           
   if (store_meteo) then
      !
      if (wind) then
         !
         zsg = FILL_VALUE
         !
         do nm = 1, np
            !
            n    = z_index_z_n(nm)
            m    = z_index_z_m(nm)
            !      
            zsg(m, n) = windu(nm)
            !
         enddo
         !
         NF90(nf90_put_var(map_file%ncid, map_file%wind_u_varid, zsg, (/1, 1, ntmapout/)))
         !
         zsg = FILL_VALUE
         !
         do nm = 1, np
            !
            n    = z_index_z_n(nm)
            m    = z_index_z_m(nm)
            !      
            zsg(m, n) = windv(nm)
            !
         enddo
         !
         NF90(nf90_put_var(map_file%ncid, map_file%wind_v_varid, zsg, (/1, 1, ntmapout/)))
         !
      endif
      !
      if (patmos) then
         !
         zsg = FILL_VALUE
         !
         do nm = 1, np
            !
            n    = z_index_z_n(nm)
            m    = z_index_z_m(nm)
            !      
            zsg(m, n) = patm(nm)
            !
         enddo
         !
         NF90(nf90_put_var(map_file%ncid, map_file%patm_varid, zsg, (/1, 1, ntmapout/)))
         !
      endif   
      !
      if (precip) then
         !
         zsg = FILL_VALUE
         !
         do nm = 1, np
            !
            n    = z_index_z_n(nm)
            m    = z_index_z_m(nm)
            !      
            zsg(m, n) = prcp(nm)*3600000
            !
         enddo
         !
         NF90(nf90_put_var(map_file%ncid, map_file%precip_varid, zsg, (/1, 1, ntmapout/)))
         !
      endif   
      !
   endif
   !
   if (snapwave) then
      ! 
      zsg = FILL_VALUE       
      !
      do nm = 1, np
         !
         n    = z_index_z_n(nm)
         m    = z_index_z_m(nm)
         !      
         zsg(m, n) = hm0(nm)
         !
      enddo
      !
      NF90(nf90_put_var(map_file%ncid, map_file%hm0_varid, zsg, (/1, 1, ntmapout/))) ! write h
      ! 
      zsg = FILL_VALUE       
      !
      do nm = 1, np
         !
         n    = z_index_z_n(nm)
         m    = z_index_z_m(nm)
         !      
         zsg(m, n) = hm0_ig(nm)
         !
      enddo
      !
      NF90(nf90_put_var(map_file%ncid, map_file%hm0ig_varid, zsg, (/1, 1, ntmapout/))) ! write h
      !            
      if (store_wave_forces) then
         !      
         zsg = FILL_VALUE       
         !
         do nm = 1, np
            !
            n    = z_index_z_n(nm)
            m    = z_index_z_m(nm)
            !      
            zsg(m, n) = fwx(nm)
            !
         enddo
         ! 
         NF90(nf90_put_var(map_file%ncid, map_file%fwx_varid, zsg, (/1, 1, ntmapout/))) ! write h
         !            
         zsg = FILL_VALUE       
         ! 
         do nm = 1, np
            !
            n    = z_index_z_n(nm)
            m    = z_index_z_m(nm)
            !      
            zsg(m, n) = fwy(nm)
            !
         enddo
         !
         NF90(nf90_put_var(map_file%ncid, map_file%fwy_varid, zsg, (/1, 1, ntmapout/))) ! write h
         !
      endif
      !            
      if (wavemaker) then
         !
         zsg = FILL_VALUE       
         !
         do nm = 1, np
            !
            n    = z_index_z_n(nm)
            m    = z_index_z_m(nm)
            !       
            zsg(m, n) = zsm(nm)
            !
         enddo
         !
         NF90(nf90_put_var(map_file%ncid, map_file%zsm_varid, zsg, (/1, 1, ntmapout/))) ! write h
         !
      endif
      !
   endif   
   !           
   NF90(nf90_sync(map_file%ncid)) !write away intermediate data ! TL: in first test it seems to be faster to let the file update than keep in memory
   !
   end subroutine

   
   subroutine ncoutput_update_quadtree_map(t,ntmapout)
      !
      ! Write time, zs, u, v  
      !
      use sfincs_data   
      use sfincs_snapwave
      use quadtree
!      use snapwave_data
      !
      implicit none   
      !
      real*8                       :: t  
      !
      integer  :: ntmapout       
      !
      integer   :: nm, nmq, n, m, nmu1, nmd1, num1, ndm1, nmu2, nmd2, num2, ndm2
      real*4    :: uz, vz, u1, u2, v1, v2, sq2
      !
      real*4, dimension(:), allocatable :: utmp, vtmp 
      !
      NF90(nf90_put_var(map_file%ncid, map_file%time_varid, t, (/ntmapout/))) ! write time
      !  
      allocate(utmp(quadtree_nr_points))
      allocate(vtmp(quadtree_nr_points))
      !
      sq2 = sqrt(2.0)
      !
      vtmp = FILL_VALUE
      !
      do nmq = 1, quadtree_nr_points
        !
        nm = index_sfincs_in_quadtree(nmq)
        !
        if (kcs(nm)>0) then
           if (subgrid) then
              if ( (zs(nm) - subgrid_z_zmin(nm)) > huthresh) then
                 vtmp(nmq) = zs(nm)
              endif
           else
              if ( (zs(nm) - zb(nm)) > huthresh) then
                 vtmp(nmq) = zs(nm)
              endif
           endif
        endif
      enddo
      !
      NF90(nf90_put_var(map_file%ncid, map_file%zs_varid, vtmp, (/1, ntmapout/))) ! write zs
      ! 
!      zsg = FILL_VALUE       
!      !
!      do nm = 1, np
!            if (subgrid) then
!                zsg(z_index_nm(1,nm),z_index_nm(2,nm)) = zs(nm) - subgrid_z_zmin(nm)
!            else
!                zsg(z_index_nm(1,nm),z_index_nm(2,nm)) = zs(nm) - zb(nm)
!            endif
!      enddo
      !
!      if (subgrid) then
!         NF90(nf90_put_var(map_file%ncid, map_file%h_varid, zsg(2:nmax-1, 2:mmax-1), (/1, 1, ntmapout/))) ! write h
      !            
      !!      
      if (store_velocity) then
         !
         utmp = FILL_VALUE
         vtmp = FILL_VALUE
         !
         do nmq = 1, quadtree_nr_points
            !
            nm = index_sfincs_in_quadtree(nmq)
            !
            if (nm>0) then
               !
               if (z_flags_type(nm) == 0) then
                  !
                  ! Regular point with four surrounding cells of the same size
                  !
                  nmd1 = z_index_uv_md1(nm)
                  nmu1 = z_index_uv_mu1(nm)
                  ndm1 = z_index_uv_nd1(nm)
                  num1 = z_index_uv_mu1(nm)
                  uz  = 0.5*(uv(nmd1) + uv(nmu1))
                  vz  = 0.5*(uv(ndm1) + uv(num1))
                  !
               else
                  !
                  ! Bit more complicated
                  !
                  nmd1 = z_index_uv_md1(nm)
                  nmu1 = z_index_uv_mu1(nm)
                  ndm1 = z_index_uv_nd1(nm)
                  num1 = z_index_uv_mu1(nm)
                  nmd2 = z_index_uv_md2(nm)
                  nmu2 = z_index_uv_mu2(nm)
                  ndm2 = z_index_uv_nd2(nm)
                  num2 = z_index_uv_nu2(nm)
                  !
                  ! Left
                  !
                  u1 = 0.0
                  if (nmd1>0 .and. nmd2==0) then   
                     u1 = uv(nmd1)
                  elseif (nmd1>0 .and. nmd2>0) then
                     u1 = 0.5*uv(nmd1) + 0.5*uv(nmd2)
                  elseif (nmd1==0 .and. nmd2>0) then   
                     u1 = uv(nmd2)
                  endif   
                  !
                  ! Right
                  !
                  u2 = 0.0
                  if (nmu1>0  .and. nmu2==0) then   
                     u2 = uv(nmu1)
                  elseif (nmu1>0 .and. nmu2>0) then
                     u2 = 0.5*uv(nmu1) + 0.5*uv(nmu2)
                  elseif (nmu2==0  .and. nmu2>0) then   
                     u2 = uv(nmu2)
                  endif   
                  !
                  ! Bottom
                  !
                  v1 = 0.0
                  if (ndm1>0 .and. ndm2==0) then   
                     v1 = uv(ndm1)
                  elseif (ndm1>0 .and. ndm2>0) then
                     v1 = 0.5*uv(ndm1) + 0.5*uv(ndm2)
                  elseif (ndm1==0 .and. ndm2>0) then   
                     v1 = uv(ndm2)
                  endif   
                  !
                  ! Right
                  !
                  v2 = 0.0
                  if (num1>0  .and. num2==0) then   
                     v2 = uv(num1)
                  elseif (num1>0 .and. num2>0) then
                     v2 = 0.5*uv(num1) + 0.5*uv(num2)
                  elseif (num2==0  .and. num2>0) then   
                     v2 = uv(num2)
                  endif   
                  !
                  uz = 0.5*(u1 + u2)
                  vz = 0.5*(v1 + v2)
                  !
               endif   
               !
               utmp(nmq) = cosrot*uz - sinrot*vz            
               vtmp(nmq) = sinrot*uz + cosrot*vz  
               !
            endif
            !
         enddo
         !
         NF90(nf90_put_var(map_file%ncid, map_file%u_varid, utmp, (/1, ntmapout/)))
         NF90(nf90_put_var(map_file%ncid, map_file%v_varid, vtmp, (/1, ntmapout/)))
         !
      endif
      !
      if (snapwave) then
         !
         vtmp = FILL_VALUE
         !
         do nmq = 1, quadtree_nr_points
            !
            nm = index_sw_in_qt(nmq)
            !
            if (nm>0) then
               ! 
               vtmp(nmq) = snapwave_H(nm)*sq2
               !
            endif   
         enddo
         !
         NF90(nf90_put_var(map_file%ncid, map_file%hm0_varid, vtmp, (/1, ntmapout/)))
         !
         vtmp = FILL_VALUE
         !
         do nmq = 1, quadtree_nr_points
            !
            nm = index_sw_in_qt(nmq)
            !
            if (nm>0) then
               ! 
               vtmp(nmq) = snapwave_H_ig(nm)*sq2
               ! 
            endif   
         enddo
         !
         NF90(nf90_put_var(map_file%ncid, map_file%hm0ig_varid, vtmp, (/1, ntmapout/)))
         !            
         if (store_wave_forces) then
            !      
            utmp = FILL_VALUE
            vtmp = FILL_VALUE
            !
            do nmq = 1, quadtree_nr_points
               !
               nm = index_sfincs_in_quadtree(nmq)
               !
               if (nm>0) then
                  utmp(nmq) = fwx(nm)
                  vtmp(nmq) = fwy(nm)
               endif
               !
            enddo
            !
            NF90(nf90_put_var(map_file%ncid, map_file%fwx_varid, utmp, (/1, ntmapout/)))
            NF90(nf90_put_var(map_file%ncid, map_file%fwy_varid, vtmp, (/1, ntmapout/)))
            !
         endif
         !
         if (wavemaker) then
            !
            vtmp = FILL_VALUE
            !
            do nmq = 1, quadtree_nr_points
               !
               nm = index_sfincs_in_quadtree(nmq)
               !
               if (nm>0) then
                  vtmp(nmq) = zsm(nm)
               endif
               !
            enddo
            !
            NF90(nf90_put_var(map_file%ncid, map_file%zsm_varid, vtmp, (/1, ntmapout/)))
            !
         endif            
         !            
      endif         
      !           
      NF90(nf90_sync(map_file%ncid)) !write away intermediate data ! TL: in first test it seems to be faster to let the file update than keep in memory
      !      
   end subroutine
   
   
   
   subroutine ncoutput_update_his(t,nthisout)
   ! Write time, zs, u, v, prcp of points 
   !
   use sfincs_data   
   use sfincs_crosssections
   use sfincs_snapwave
   !
   implicit none   
   !
   integer :: iobs, nm, icrs, ip
   !
   integer :: nthisout      
   integer :: nmd1, nmu1, ndm1, num1, nmd2, nmu2, ndm2, num2
   !
   real*4                  :: u1, u2, v1, v2, uz, vz
   real*8                  :: t
!   real*4, dimension(nobs) :: zobs, hobs
   real*4, dimension(nobs) :: uobs
   real*4, dimension(nobs) :: vobs   
   real*4, dimension(nobs) :: uvmag   
   real*4, dimension(nobs) :: uvdir   
   real*4, dimension(nobs) :: tprcp
   real*4, dimension(nobs) :: tqinf
   real*4, dimension(nobs) :: tS_effective
   real*4, dimension(nobs) :: tpatm
   real*4, dimension(nobs) :: twndmag
   real*4, dimension(nobs) :: twnddir
   real*4, dimension(nobs) :: hm0obs
   real*4, dimension(nobs) :: hm0igobs
   real*4, dimension(nobs) :: zsmobs
   real*4, dimension(nobs) :: tpobs
   real*4, dimension(nobs) :: wavdirobs
   real*4, dimension(nobs) :: dirsprobs
   real*4, dimension(:), allocatable :: qq
   !
   zobs         = FILL_VALUE
   zsmobs       = FILL_VALUE
   hobs         = FILL_VALUE         
   hm0obs       = FILL_VALUE         
   hm0igobs     = FILL_VALUE         
   tpobs        = FILL_VALUE         
   wavdirobs    = FILL_VALUE         
   dirsprobs    = FILL_VALUE         
   tprcp        = FILL_VALUE
   tqinf        = FILL_VALUE
   tS_effective = FILL_VALUE
   tpatm        = FILL_VALUE
   twndmag      = FILL_VALUE
   twnddir      = FILL_VALUE
   !
   do iobs = 1, nobs ! determine zs and prcp of obervation points at required timestep
      !
      nm = nmindobs(iobs)
      !
      if (nm>0) then
         !
         zobs(iobs)  = zs(nm)
         !
         if (subgrid) then
            hobs(iobs)  = zs(nm) - subgrid_z_zmin(nm)
         else
            hobs(iobs)  = zs(nm) - zb(nm)
         endif
         !
         if (store_velocity) then
            !
            if (z_flags_type(nm) == 0) then
               !
               ! Regular point with four surrounding cells of the same size
               !
               nmd1 = z_index_uv_md1(nm)
               nmu1 = z_index_uv_mu1(nm)
               ndm1 = z_index_uv_nd1(nm)
               num1 = z_index_uv_mu1(nm)
               uz  = 0.5*(uv(nmd1) + uv(nmu1))
               vz  = 0.5*(uv(ndm1) + uv(num1))
               !
            else
               !
               ! Bit more complicated
               !
               nmd1 = z_index_uv_md1(nm)
               nmu1 = z_index_uv_mu1(nm)
               ndm1 = z_index_uv_nd1(nm)
               num1 = z_index_uv_mu1(nm)
               nmd2 = z_index_uv_md2(nm)
               nmu2 = z_index_uv_mu2(nm)
               ndm2 = z_index_uv_nd2(nm)
               num2 = z_index_uv_nu2(nm)
               !
               ! Left
               !
               u1 = 0.0
               if (nmd1>0 .and. nmd2==0) then   
                  u1 = uv(nmd1)
               elseif (nmd1>0 .and. nmd2>0) then
                  u1 = 0.5*uv(nmd1) + 0.5*uv(nmd2)
               elseif (nmd1==0 .and. nmd2>0) then   
                  u1 = uv(nmd2)
               endif   
               !
               ! Right
               !
               u2 = 0.0
               if (nmu1>0  .and. nmu2==0) then   
                  u2 = uv(nmu1)
               elseif (nmu1>0 .and. nmu2>0) then
                  u2 = 0.5*uv(nmu1) + 0.5*uv(nmu2)
               elseif (nmu2==0  .and. nmu2>0) then   
                  u2 = uv(nmu2)
               endif   
               !
               ! Bottom
               !
               v1 = 0.0
               if (ndm1>0 .and. ndm2==0) then   
                  v1 = uv(ndm1)
               elseif (ndm1>0 .and. ndm2>0) then
                  v1 = 0.5*uv(ndm1) + 0.5*uv(ndm2)
               elseif (ndm1==0 .and. ndm2>0) then   
                  v1 = uv(ndm2)
               endif   
               !
               ! Right
               !
               v2 = 0.0
               if (num1>0  .and. num2==0) then   
                  v2 = uv(num1)
               elseif (num1>0 .and. num2>0) then
                  v2 = 0.5*uv(num1) + 0.5*uv(num2)
               elseif (num2==0  .and. num2>0) then   
                  v2 = uv(num2)
               endif   
               !
               uz = 0.5*(u1 + u2)
               vz = 0.5*(v1 + v2)
               !
            endif   
            !
            uobs(iobs)  = cosrot*uz - sinrot*vz                         
            vobs(iobs)  = sinrot*uz + cosrot*vz
            uvmag(iobs) = sqrt(uobs(iobs)**2 + vobs(iobs)**2)
            uvdir(iobs) = atan2(vobs(iobs), uobs(iobs))*180/pi
            !
         endif
         !
         if (infiltration) then
            !
            tqinf(iobs) = qinfmap(nm)*3.6e3*1.0e3 ! show as mm/hr
            ! 
            ! Output for CN method
            !
            if (inftype == 'cnb') then
               tS_effective(iobs) = qinffield2(nm)
            endif
            !
         endif  
         !
         if (store_meteo) then
            !
            if (wind) then
               !
               twndmag(iobs) = sqrt(windu(nm)**2 + windv(nm)**2)
               twnddir(iobs) = 270.0 - atan2(windv(nm), windu(nm))*180/pi
               if (twnddir(iobs)<0.0) twnddir(iobs) = twnddir(iobs) + 360.0
               if (twnddir(iobs)>360.0) twnddir(iobs) = twnddir(iobs) - 360.0
               !
            endif   
            !
            if (patmos) then
               !
               tpatm(iobs) = patm(nm)
               !
            endif   
            !
            if (precip) then
               !
               tprcp(iobs) = prcp(nm)*3600000 ! show as mm/hr
               !
            endif
            !         
         endif
         !
         if (snapwave) then
            !
            hm0obs(iobs)   = hm0(nm)
            hm0igobs(iobs) = hm0_ig(nm)
            tpobs(iobs)    = snapwave_tpmean
            !
            if (store_wave_direction) then
               !
               wavdirobs(iobs)   = mean_wave_direction(nm)
               dirsprobs(iobs)   = wave_directional_spreading(nm)
               !
            endif
            !            
            if (wavemaker) then
               !
               zsmobs(iobs)   = zsm(nm)
               !
            endif   
            !
         endif   
         !
      endif
   enddo   
   !   
   NF90(nf90_put_var(his_file%ncid, his_file%time_varid, t, (/nthisout/))) ! write time
   !   
   NF90(nf90_put_var(his_file%ncid, his_file%zs_varid, zobs, (/1, nthisout/))) ! write point_zs
   !
   if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then   
      !
      NF90(nf90_put_var(his_file%ncid, his_file%h_varid, hobs, (/1, nthisout/))) ! write point_h   
      !
   endif
   !
   if (infiltration) then
      !
      NF90(nf90_put_var(his_file%ncid, his_file%qinf_varid, tqinf, (/1, nthisout/))) ! write qinf
      !
      if (inftype == 'cnb') then
         NF90(nf90_put_var(his_file%ncid, his_file%S_varid, tS_effective, (/1, nthisout/))) ! write S
      endif
      !
   endif
   !
   if (snapwave) then  
      !
      NF90(nf90_put_var(his_file%ncid, his_file%hm0_varid, hm0obs, (/1, nthisout/)))
      NF90(nf90_put_var(his_file%ncid, his_file%hm0ig_varid, hm0igobs, (/1, nthisout/)))
      NF90(nf90_put_var(his_file%ncid, his_file%tp_varid, tpobs, (/1, nthisout/)))
      !
      if (store_wave_direction) then
         !
         NF90(nf90_put_var(his_file%ncid, his_file%wavdir_varid, wavdirobs, (/1, nthisout/)))
         NF90(nf90_put_var(his_file%ncid, his_file%dirspr_varid, dirsprobs, (/1, nthisout/)))
         !
      endif
      !            
      if (wavemaker) then
         !
         NF90(nf90_put_var(his_file%ncid, his_file%zsm_varid, zsmobs, (/1, nthisout/)))
         !
      endif
      !
   endif
   !
   if (store_meteo) then
      !
      if (wind) then
         !
         NF90(nf90_put_var(his_file%ncid, his_file%wind_speed_varid, twndmag, (/1, nthisout/)))
         NF90(nf90_put_var(his_file%ncid, his_file%wind_dir_varid,   twnddir, (/1, nthisout/)))
        !
      endif      
      !
      if (patmos) then
         !
         NF90(nf90_put_var(his_file%ncid, his_file%patm_varid, tpatm, (/1, nthisout/))) ! write prcp
         !
      endif   
      !
      if (precip) then
         !
         NF90(nf90_put_var(his_file%ncid, his_file%prcp_varid, tprcp, (/1, nthisout/))) ! write prcp
         !
      endif
      !   
   endif
   !
   if (nrcrosssections>0) then
      !
      !$acc update host(q)
      !      
      ! Get fluxes through cross sections
      !
      call get_discharges_through_crosssections(qq)
      !
      NF90(nf90_put_var(his_file%ncid, his_file%discharge_varid, qq, (/1, nthisout/))) ! write prcp
      !
   endif
   !
   if (store_velocity) then
      !
      NF90(nf90_put_var(his_file%ncid, his_file%u_varid, uobs, (/1, nthisout/)))
      NF90(nf90_put_var(his_file%ncid, his_file%v_varid, vobs, (/1, nthisout/)))   
      NF90(nf90_put_var(his_file%ncid, his_file%uvmag_varid, uvmag, (/1, nthisout/)))   
      NF90(nf90_put_var(his_file%ncid, his_file%uvdir_varid, uvdir, (/1, nthisout/)))   
      !
   endif
   !
   NF90(nf90_sync(his_file%ncid)) !write away intermediate data ! TL: in first test it seems to be faster to let the file update than keep in memory
   !   
   end subroutine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   !
   subroutine ncoutput_update_max(t,ntmaxout)
   !
   ! write zsmax per dtmaxout
   !
   use sfincs_data   
   !
   implicit none   
   !
   integer :: nm, n, m
   !
   real*8                       :: t  
   !
   integer  :: ntmaxout   
   !
   real*4, dimension(:,:), allocatable :: zstmp
   !
   allocate(zstmp(mmax, nmax))
   !
   zstmp = FILL_VALUE
   !
   if (subgrid) then   
      !
      do nm = 1, np
         !
         n    = z_index_z_n(nm)
         m    = z_index_z_m(nm)
         !      
         if ( (zsmax(nm) - subgrid_z_zmin(nm)) > huthresh) then
            zstmp(m, n) = zsmax(nm) 
         endif
      enddo
   else
      do nm = 1, np       
         !
         n    = z_index_z_n(nm)
         m    = z_index_z_m(nm)
         !      
         if ( (zsmax(nm) - zb(nm)) > huthresh) then
            zstmp(m, n) = zsmax(nm) 
         endif      
      enddo
   endif
   !
   NF90(nf90_put_var(map_file%ncid, map_file%timemax_varid, t, (/ntmaxout/))) ! write time_max
   NF90(nf90_put_var(map_file%ncid, map_file%zsmax_varid, zstmp, (/1, 1, ntmaxout/))) ! write zsmax      
   !
   ! Write maximum water depth (optional)
   if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then
      !
      zstmp = FILL_VALUE
      !
      if (subgrid) then   
         do nm = 1, np
            !
            n    = z_index_z_n(nm)
            m    = z_index_z_m(nm)             
            !
            if ( (zsmax(nm) - subgrid_z_zmin(nm)) > huthresh) then
               zstmp(m, n) = zsmax(nm) - subgrid_z_zmin(nm)
            endif
            !
         enddo
      else
         do nm = 1, np       
            !
            n    = z_index_z_n(nm)
            m    = z_index_z_m(nm)
            !      
            if ( (zsmax(nm) - zb(nm)) > huthresh) then
               zstmp(m, n) = zsmax(nm) - zb(nm)
            endif      
         enddo
      endif 
      NF90(nf90_put_var(map_file%ncid, map_file%hmax_varid, zstmp, (/1, 1, ntmaxout/))) ! write hmax   
   endif
   !
   ! Write cumulative rainfall
   ! 
   if (store_cumulative_precipitation) then  
      !
      ! Precipitation
      ! 
      zstmp = FILL_VALUE
      !
      do nm = 1, np       
         !
         n    = z_index_z_n(nm)
         m    = z_index_z_m(nm)
         !      
         zstmp(m, n) = cumprcp(nm)*1000
         !
      enddo
      !
      NF90(nf90_put_var(map_file%ncid, map_file%cumprcp_varid, zstmp, (/1, 1, ntmaxout/))) ! write zsmax   
      !
      ! Infiltration
      !
      zstmp = FILL_VALUE
      !
      do nm = 1, np
         n              = z_index_z_n(nm)
         m              = z_index_z_m(nm)
         zstmp(m, n)    = cuminf(nm) 
      enddo
      !
      NF90(nf90_put_var(map_file%ncid, map_file%cuminf_varid, zstmp, (/1, 1, ntmaxout/))) ! write cuminf   
      !
   endif
   !
   ! Maximum flow velocity
   !
   if (store_maximum_velocity) then
      zstmp = FILL_VALUE
      do nm = 1, np   
         n              = z_index_z_n(nm)
         m              = z_index_z_m(nm)
         zstmp(m, n)    = vmax(nm)
      enddo
      NF90(nf90_put_var(map_file%ncid, map_file%vmax_varid, zstmp, (/1, 1, ntmaxout/))) ! write vmax
   endif
   !
   ! Duration wet cell
   !
   if (store_twet) then
      zstmp = FILL_VALUE
      do nm = 1, np
         n              = z_index_z_n(nm)
         m              = z_index_z_m(nm)
         zstmp(m, n)     = twet(nm) 
      enddo
      NF90(nf90_put_var(map_file%ncid, map_file%tmax_varid, zstmp, (/1, 1, ntmaxout/))) ! write tmax   
   endif
   !
   ! Maximum wind speed
   if (wind .and. store_wind_max .and. meteo3d) then 
      zstmp = FILL_VALUE
      do nm = 1, np
         n              = z_index_z_n(nm)
         m              = z_index_z_m(nm)
         zstmp(m, n)    = windmax(nm) 
      enddo
      NF90(nf90_put_var(map_file%ncid, map_file%windmax_varid, zstmp, (/1, 1, ntmaxout/))) ! write windmax   
   endif
   !   
   end subroutine


   
   subroutine ncoutput_update_quadtree_max(t,ntmaxout)
   !
   ! write zsmax per dtmaxout
   !
   use sfincs_data  
   use sfincs_snapwave
   use quadtree
   !
   implicit none   
   !
   integer                              :: nmq, nm, n, m, ntmaxout
   real*8                               :: t  
   !
   real*4, dimension(:), allocatable    :: zstmp
   allocate(zstmp(quadtree_nr_points))
   !
   zstmp = FILL_VALUE
   !
   ! Write maximum water level
   do nmq = 1, quadtree_nr_points
       !
       nm = index_sfincs_in_quadtree(nmq)
       !
       if (kcs(nm)>0) then
           if (subgrid) then
               if ( (zs(nm) - subgrid_z_zmin(nm)) > huthresh) then
                   zstmp(nmq) = zs(nm)
               endif
           else
              if ( (zs(nm) - zb(nm)) > huthresh) then
                  zstmp(nmq) = zs(nm)
              endif
           endif
       endif
   enddo
   !
   NF90(nf90_put_var(map_file%ncid, map_file%timemax_varid, t, (/ntmaxout/)))       ! write time_max
   NF90(nf90_put_var(map_file%ncid, map_file%zsmax_varid, zstmp, (/1, ntmaxout/)))  ! write zsmax   
   !
   ! Write cumulative precipitation
   if (store_cumulative_precipitation) then  
       zstmp = FILL_VALUE       
       do nmq = 1, quadtree_nr_points
           nm = index_sfincs_in_quadtree(nmq)
           if (kcs(nm)>0) then
               zstmp(nmq) = cumprcp(nm)*1000
           endif
       enddo
       NF90(nf90_put_var(map_file%ncid, map_file%cumprcp_varid, zstmp, (/1, ntmaxout/))) ! write cumprcp
   endif   
   !
   ! Duration wet cell
   if (store_twet) then
       zstmp = FILL_VALUE       
       do nmq = 1, quadtree_nr_points
           nm = index_sfincs_in_quadtree(nmq)
           if (kcs(nm)>0) then
               zstmp(nmq) = twet(nm)
           endif
       enddo
      NF90(nf90_put_var(map_file%ncid, map_file%tmax_varid, zstmp, (/1, ntmaxout/))) ! write tmax   
   endif
   !
   ! Maximum wind speed
   if (wind .and. store_wind_max .and. meteo3d) then 
       zstmp = FILL_VALUE       
       do nmq = 1, quadtree_nr_points
           nm = index_sfincs_in_quadtree(nmq)
           if (kcs(nm)>0) then
               zstmp(nmq) = windmax(nm)
           endif
       enddo
      NF90(nf90_put_var(map_file%ncid, map_file%windmax_varid, zstmp, (/1, ntmaxout/))) ! write windmax   
   endif
   !   
   ! Cumulative infiltration
   if (infiltration) then
      zstmp = FILL_VALUE
      do nmq = 1, quadtree_nr_points
          nm = index_sfincs_in_quadtree(nmq)
          if (kcs(nm)>0) then
              zstmp(nmq) = cuminf(nm)
          endif
      enddo
      NF90(nf90_put_var(map_file%ncid, map_file%cuminf_varid, zstmp, (/1, ntmaxout/))) ! write cuminf   
   endif
   !
   end subroutine   
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   !
   subroutine ncoutput_map_finalize() 
   ! Add total runtime, dtavg to file and close
   !
   use sfincs_data
   !   
   implicit none   
   !   
   NF90(nf90_put_var(map_file%ncid, map_file%total_runtime_varid, tfinish_all - tstart_all)) 
   NF90(nf90_put_var(map_file%ncid, map_file%average_dt_varid,  dtavg)) 
   !
   NF90(nf90_close(map_file%ncid))
   !
   end subroutine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   subroutine ncoutput_his_finalize()
   ! Add total runtime, dtavg to file and close
   !
   use sfincs_data
   !   
   implicit none   
   !   
   if (nobs==0 .and. nrcrosssections==0) then ! If no observation points nor cross-sections his file is not created        
        return
   endif   
   !
   NF90(nf90_put_var(his_file%ncid, his_file%total_runtime_varid, tfinish_all - tstart_all)) 
   NF90(nf90_put_var(his_file%ncid, his_file%average_dt_varid,  dtavg)) 
   !   
   NF90(nf90_close(his_file%ncid))
   !
   end subroutine
   !  
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   subroutine ncoutput_add_params(ncid, varid)
   ! Add user params to netcdf file (both map & his)
   use sfincs_data   
   !
   implicit none   
   !
   integer :: ncid, varid
   !
   ! add user input params for reproducability
        ! User input
        NF90(nf90_def_var(ncid, 'inp', NF90_INT, varid)) ! use order of sfincs_input.f90
        NF90(nf90_put_att(ncid, varid, 'mmax', mmax)) !without dummy cells
        NF90(nf90_put_att(ncid, varid, 'nmax', nmax))
        NF90(nf90_put_att(ncid, varid, 'dx', dx)) 
        NF90(nf90_put_att(ncid, varid, 'dy', dy))  
        NF90(nf90_put_att(ncid, varid, 'x0', x0))        
        NF90(nf90_put_att(ncid, varid, 'y0', y0))        
        NF90(nf90_put_att(ncid, varid, 'rotation', rotation))   
        NF90(nf90_put_att(ncid, varid, 'tref', trefstr))   
        NF90(nf90_put_att(ncid, varid, 'tstart', tstartstr))   
        NF90(nf90_put_att(ncid, varid, 'tstop', tstopstr))   
        NF90(nf90_put_att(ncid, varid, 'tspinup', tspinup))   
        NF90(nf90_put_att(ncid, varid, 't0out', t0out))   
        NF90(nf90_put_att(ncid, varid, 't1out', t1out))   
        NF90(nf90_put_att(ncid, varid, 'dtout', dtmapout))    
        NF90(nf90_put_att(ncid, varid, 'dtmaxout', dtmaxout))    
        NF90(nf90_put_att(ncid, varid, 'dtrstout', dtrstout))
        NF90(nf90_put_att(ncid, varid, 'trstout', trst))        
        NF90(nf90_put_att(ncid, varid, 'dthisout', dthisout))
        NF90(nf90_put_att(ncid, varid, 'dtwave',dtwave))        
        NF90(nf90_put_att(ncid, varid, 'dtwnd',dtwindupd))        
        NF90(nf90_put_att(ncid, varid, 'alpha',alfa))        
        NF90(nf90_put_att(ncid, varid, 'theta',theta))        
        NF90(nf90_put_att(ncid, varid, 'manning',manning))        
        NF90(nf90_put_att(ncid, varid, 'manning_land',manning_land))
        NF90(nf90_put_att(ncid, varid, 'manning_sea',manning_sea))        
        NF90(nf90_put_att(ncid, varid, 'rgh_lev_land',rghlevland))        
        NF90(nf90_put_att(ncid, varid, 'zsini',zini))        
        NF90(nf90_put_att(ncid, varid, 'qinf',qinf))        
        NF90(nf90_put_att(ncid, varid, 'igperiod',tig))        
        NF90(nf90_put_att(ncid, varid, 'dtmax',dtmax))        
        NF90(nf90_put_att(ncid, varid, 'dtmin',dtmin))        
        NF90(nf90_put_att(ncid, varid, 'huthresh',huthresh))        
        NF90(nf90_put_att(ncid, varid, 'rhoa',rhoa))        
        NF90(nf90_put_att(ncid, varid, 'rhow',rhow))        
        NF90(nf90_put_att(ncid, varid, 'maxlev',max_elev))        
        NF90(nf90_put_att(ncid, varid, 'inputformat',inputtype))        
        NF90(nf90_put_att(ncid, varid, 'outputformat',outputtype))   
        NF90(nf90_put_att(ncid, varid, 'outputtype_map',outputtype_map))
        NF90(nf90_put_att(ncid, varid, 'outputtype_his',outputtype_his))
        NF90(nf90_put_att(ncid, varid, 'bndtype',bndtype))
        NF90(nf90_put_att(ncid, varid, 'advection',iadvection))  
        NF90(nf90_put_att(ncid, varid, 'nfreqsig',nfreqsig))  
        NF90(nf90_put_att(ncid, varid, 'freqminig',freqminig))  
        NF90(nf90_put_att(ncid, varid, 'freqmaxig',freqmaxig))  
        NF90(nf90_put_att(ncid, varid, 'latitude',latitude))  
        NF90(nf90_put_att(ncid, varid, 'pavbnd',pavbnd))  
        NF90(nf90_put_att(ncid, varid, 'gapres',gapres))  
        NF90(nf90_put_att(ncid, varid, 'baro',baro))  
        NF90(nf90_put_att(ncid, varid, 'utmzone',utmzone))  
        NF90(nf90_put_att(ncid, varid, 'epsg',epsg))  
        NF90(nf90_put_att(ncid, varid, 'epsg_code',epsg_code))  
        NF90(nf90_put_att(ncid, varid, 'stopdepth',stopdepth))  
        NF90(nf90_put_att(ncid, varid, 'advlim',advlim))  
        NF90(nf90_put_att(ncid, varid, 'qinf_zmin',qinf_zmin))    
        NF90(nf90_put_att(ncid, varid, 'horton_decay', horton_decay))  
        NF90(nf90_put_att(ncid, varid, 'horton_time', horton_time))     
        NF90(nf90_put_att(ncid, varid, 'btfilter', btfilter))                     
        NF90(nf90_put_att(ncid, varid, 'sfacinf', sfacinf))  
        NF90(nf90_put_att(ncid, varid, 'radstr', iradstr)) 
        NF90(nf90_put_att(ncid, varid, 'crsgeo',igeo)) 
        NF90(nf90_put_att(ncid, varid, 'coriolis',icoriolis)) 
        NF90(nf90_put_att(ncid, varid, 'amprblock',iamprblock)) 
        NF90(nf90_put_att(ncid, varid, 'spwmergefrac',spw_merge_frac)) 
        NF90(nf90_put_att(ncid, varid, 'global',iglobal)) 
        NF90(nf90_put_att(ncid, varid, 'nuvisc',nuvisc)) 
        NF90(nf90_put_att(ncid, varid, 'spinup_meteo', ispinupmeteo)) 
        NF90(nf90_put_att(ncid, varid, 'waveage',waveage)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave', isnapwave)) 
        NF90(nf90_put_att(ncid, varid, 'wmtfilter', wmtfilter))         
        NF90(nf90_put_att(ncid, varid, 'wmfred',wavemaker_freduv))         
        !
        ! Domain
        !
        NF90(nf90_put_att(ncid, varid, 'qtrfile',qtrfile))        
        NF90(nf90_put_att(ncid, varid, 'depfile',depfile))        
        NF90(nf90_put_att(ncid, varid, 'inifile',zsinifile))        
        NF90(nf90_put_att(ncid, varid, 'rstfile',rstfile))        
        NF90(nf90_put_att(ncid, varid, 'mskfile',mskfile))        
        NF90(nf90_put_att(ncid, varid, 'indexfile',indexfile))        
        NF90(nf90_put_att(ncid, varid, 'cstfile',cstfile))             
        NF90(nf90_put_att(ncid, varid, 'sbgfile',sbgfile))        
        NF90(nf90_put_att(ncid, varid, 'thdfile',thdfile))        
        NF90(nf90_put_att(ncid, varid, 'weirfile',weirfile))        
        NF90(nf90_put_att(ncid, varid, 'manningfile',manningfile))    
        NF90(nf90_put_att(ncid, varid, 'drnfile',drnfile))    
        !
        ! Forcing
        !
        NF90(nf90_put_att(ncid, varid, 'bndfile',bndfile))        
        NF90(nf90_put_att(ncid, varid, 'bzsfile',bzsfile))        
        NF90(nf90_put_att(ncid, varid, 'bzifile',bzifile))        
        NF90(nf90_put_att(ncid, varid, 'bwvfile',bwvfile))        
        NF90(nf90_put_att(ncid, varid, 'bhsfile',bhsfile))        
        NF90(nf90_put_att(ncid, varid, 'btpfile',btpfile))        
        NF90(nf90_put_att(ncid, varid, 'bwdfile',bwdfile))  
        NF90(nf90_put_att(ncid, varid, 'bdsfile',bdsfile))   
        NF90(nf90_put_att(ncid, varid, 'wfpfile',wfpfile))   
        NF90(nf90_put_att(ncid, varid, 'whifile',whifile))   
        NF90(nf90_put_att(ncid, varid, 'wtifile',wtifile))   
        NF90(nf90_put_att(ncid, varid, 'wstfile',wstfile))   
        NF90(nf90_put_att(ncid, varid, 'srcfile',srcfile))        
        NF90(nf90_put_att(ncid, varid, 'disfile',disfile))        
        NF90(nf90_put_att(ncid, varid, 'spwfile',spwfile))        
        NF90(nf90_put_att(ncid, varid, 'wndfile',wndfile))        
        NF90(nf90_put_att(ncid, varid, 'precipfile',prcpfile))        
        NF90(nf90_put_att(ncid, varid, 'amufile',amufile))        
        NF90(nf90_put_att(ncid, varid, 'amvfile',amvfile))  
        NF90(nf90_put_att(ncid, varid, 'ampfile',amvfile))              
        NF90(nf90_put_att(ncid, varid, 'amprfile',amprfile))  
        NF90(nf90_put_att(ncid, varid, 'qinffile',qinffile))   
        NF90(nf90_put_att(ncid, varid, 'scsfile',scsfile)) 
        NF90(nf90_put_att(ncid, varid, 'scsfile_Smax',scsfile_Smax)) 
        NF90(nf90_put_att(ncid, varid, 'scsfile_Se',scsfile_Se)) 
        NF90(nf90_put_att(ncid, varid, 'scsfile_kr',scsfile_kr)) 
        NF90(nf90_put_att(ncid, varid, 'z0lfile',z0lfile)) 
        NF90(nf90_put_att(ncid, varid, 'wvmfile',wvmfile)) 
        !
        ! Netcdf input
        NF90(nf90_put_att(ncid, varid, 'netbndbzsbzifile',netbndbzsbzifile))
        NF90(nf90_put_att(ncid, varid, 'netsrcdisfile',netsrcdisfile))     
        NF90(nf90_put_att(ncid, varid, 'netamuamvfile',netamuamvfile))     
        NF90(nf90_put_att(ncid, varid, 'netamprfile',netamprfile))  
        NF90(nf90_put_att(ncid, varid, 'netampfile',netampfile))        
        !
        ! Output
        !
        NF90(nf90_put_att(ncid, varid, 'obsfile',obsfile))   
        NF90(nf90_put_att(ncid, varid, 'nobs',nobs))                    
        NF90(nf90_put_att(ncid, varid, 'crsfile',crsfile))   
        !
        NF90(nf90_put_att(ncid, varid, 'storevelmax',storevelmax)) 
        NF90(nf90_put_att(ncid, varid, 'storevel',storevel)) 
        NF90(nf90_put_att(ncid, varid, 'storecumprcp',storecumprcp)) 
        NF90(nf90_put_att(ncid, varid, 'storetwet',storetwet)) 
        NF90(nf90_put_att(ncid, varid, 'storehsubgrid',storehsubgrid)) 
        NF90(nf90_put_att(ncid, varid, 'twet_threshold',twet_threshold)) 
        NF90(nf90_put_att(ncid, varid, 'store_tsunami_arrival_time',itsunamitime)) 
        NF90(nf90_put_att(ncid, varid, 'tsunami_arrival_threshold',tsunami_arrival_threshold)) 
        NF90(nf90_put_att(ncid, varid, 'storeqdrain',storeqdrain)) 
        NF90(nf90_put_att(ncid, varid, 'storezvolume',storezvolume)) 
        NF90(nf90_put_att(ncid, varid, 'writeruntime',wrttimeoutput)) 
        NF90(nf90_put_att(ncid, varid, 'debug',idebug)) 
        NF90(nf90_put_att(ncid, varid, 'storemeteo',storemeteo)) 
        NF90(nf90_put_att(ncid, varid, 'storemaxwind',iwindmax)) 
        NF90(nf90_put_att(ncid, varid, 'storefw',istorefw))         
        NF90(nf90_put_att(ncid, varid, 'storewavdir', istorewavdir)) 
        !
        NF90(nf90_put_att(ncid, varid, 'cdnrb', cd_nr))   
        NF90(nf90_put_att(ncid, varid, 'cdwnd', cd_wnd))        
        NF90(nf90_put_att(ncid, varid, 'cdval', cd_val))  
        !
        ! Internal code switches - note, you can't store logicals in netcdf, only integers for these type of switches
        NF90(nf90_put_att(ncid, varid, 'manning2d', imanning2d))   
        NF90(nf90_put_att(ncid, varid, 'subgrid', isubgrid))   
        NF90(nf90_put_att(ncid, varid, 'viscosity', iviscosity))   
        NF90(nf90_put_att(ncid, varid, 'wavemaker', iwavemaker))         
        NF90(nf90_put_att(ncid, varid, 'wavemaker_spectrum', iwavemaker_spectrum))         
   end subroutine
   !   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !   
   subroutine handle_err(status,file,line)
      !
      integer, intent ( in)    :: status
      character(*), intent(in) :: file
      integer, intent ( in)    :: line
      integer :: status2

      if(status /= nf90_noerr) then
         !UNIT=6 for stdout and UNIT=0 for stderr.
         write(0,'("NETCDF ERROR: ",a,i6,":",a)') file,line,trim(nf90_strerror(status))
         write(0,*) 'closing file'
         status2 = nf90_close(map_file%ncid)
         if (status2 /= nf90_noerr) then
            write(0,*) 'NETCDF ERROR: ', __FILE__,__LINE__,trim(nf90_strerror(status2))
         end if
      end if
   end subroutine handle_err
   !
   end module
