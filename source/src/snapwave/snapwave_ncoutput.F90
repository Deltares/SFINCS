! Unstructured mesh
#define NF90(nf90call) call handle_err(nf90call,__FILE__,__LINE__)     
module snapwave_ncoutput
   !
   use netcdf 
   use snapwave_date
   !
   implicit none
   !
   type map_type
       integer :: ncid   
       integer :: mesh2d_nNodes_dimid, mesh2d_nEdges_dimid, Two_dimid, mesh2d_nFaces_dimid, mesh2d_nMax_face_nodes_dimid, mesh2d_nMax_edge_nodes_dimid, time_dimid
       integer :: mesh2d_node_x_varid
       integer :: mesh2d_node_y_varid
       integer :: mesh2d_node_z_varid
       integer :: mesh2d_face_nodes_varid
       integer :: mesh2d_edge_nodes_varid
       integer :: crs_varid, grid_varid 
       integer :: time_varid
       integer :: hm0_varid
       integer :: hm0_ig_varid
       integer :: tp_varid
       integer :: wd_varid
   end type
   type his_type
       integer :: ncid   
       integer :: time_dimid 
       integer :: points_dimid, pointnamelength_dimid
       integer :: runtime_dimid
       integer :: point_x_varid, point_y_varid, station_x_varid, station_y_varid, crs_varid  
       integer :: station_id_varid, station_name_varid
       !integer :: zb_varid
       integer :: time_varid
       integer :: hm0_varid, tp_varid, wavdir_varid, dirspr_varid
       integer :: vmag_varid, vdir_varid
       integer :: inp_varid, total_runtime_varid, average_dt_varid   
   end type
   type(map_type) :: map_file
   type(his_type) :: his_file
   !
   real*4, parameter :: FILL_VALUE = -999999  
   !
contains
   !
   subroutine ncoutput_init()
   !
   use snapwave_data
   !
   if (map_filename/='') call ncoutput_map_init()
   if (his_filename/='') call ncoutput_his_init()
   !
   end subroutine
   !
   subroutine ncoutput_update(t,it)
   !
   use snapwave_data
   !
   real*8       :: t
   integer      :: it
   !
   if (map_filename/='') call ncoutput_update_map(t,it)
   if (his_filename/='') call ncoutput_update_his(t,it)
   !
   end subroutine
   !
   subroutine ncoutput_finalize()
   !
   use snapwave_data
   !
   if (map_filename/='') call ncoutput_map_finalize()
   if (his_filename/='') call ncoutput_his_finalize()
   !
   end subroutine
   !
   subroutine ncoutput_map_init()
   !
   ! 1. Initialise dimensions/variables/attributes
   ! 2. write grid/msk/zb to file
   !
   use snapwave_data   
   !
   implicit none   
   !   
   integer                      :: k
   !
   real*4, dimension(:), allocatable :: zsg
   !
   allocate(buf(no_nodes))
   !
   NF90(nf90_create(map_filename, NF90_CLOBBER, map_file%ncid))
   !
   ! Create dimensions
   ! grid, time, points
   NF90(nf90_def_dim(map_file%ncid, 'mesh2d_nNodes', no_nodes, map_file%mesh2d_nNodes_dimid))
   NF90(nf90_def_dim(map_file%ncid, 'mesh2d_nFaces', no_faces, map_file%mesh2d_nFaces_dimid))
   NF90(nf90_def_dim(map_file%ncid, 'mesh2d_nMax_face_nodes', 4, map_file%mesh2d_nMax_face_nodes_dimid))
   NF90(nf90_def_dim(map_file%ncid, 'mesh2d_nEdges', no_edges, map_file%mesh2d_nEdges_dimid))
   NF90(nf90_def_dim(map_file%ncid, 'mesh2d_nMax_edge_nodes', 2, map_file%mesh2d_nMax_edge_nodes_dimid))
   NF90(nf90_def_dim(map_file%ncid, 'time', NF90_UNLIMITED, map_file%time_dimid)) ! time
!   NF90(nf90_def_dim(map_file%ncid, 'runtime', 1, map_file%runtime_dimid)) ! total_runtime, average_dt       
   !
   ! Some metadata attributes 
   NF90(nf90_put_att(map_file%ncid,nf90_global, "Conventions", "Conventions = 'CF-1.6, SGRID-0.3")) 
   NF90(nf90_put_att(map_file%ncid,nf90_global, "Build-Revision-Date-Netcdf-library", trim(nf90_inq_libvers()))) ! version of netcdf library
   NF90(nf90_put_att(map_file%ncid,nf90_global, "Producer", "SnapWave"))
!   NF90(nf90_put_att(map_file%ncid,nf90_global, "Build-Revision", trim(build_revision))) 
!   NF90(nf90_put_att(map_file%ncid,nf90_global, "Build-Date", trim(build_date)))
   NF90(nf90_put_att(map_file%ncid,nf90_global, "title", "SnapWave map netcdf output"))   
   !
   ! add input params for reproducability
!   call ncoutput_add_params(map_file%ncid,map_file%inp_varid)   
   !
   ! Create variables
   ! Domain
   NF90(nf90_def_var(map_file%ncid, 'mesh2d', NF90_INT, map_file%grid_varid)) ! For neat grid clarification
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'cf_role', 'mesh_topology'))
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'long_name', 'Topology data of 2D network'))
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'topology_dimension', 2))
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'node_coordinates', 'mesh2d_node_x mesh2d_node_y'))
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'node_dimension', 'mesh2d_nNodes'))   
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'face_node_connectivity', 'mesh2d_face_nodes')) 
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'face_dimension', 'mesh2d_nFaces')) 
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'max_face_nodes_dimension', 'mesh2d_nMax_face_nodes')) 
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'edge_node_connectivity', 'mesh2d_edge_nodes')) 
   NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'edge_dimension', 'mesh2d_nEdges')) 
  !NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'edge_face_connectivity', 'mesh2d_edge_faces')) 
  !NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'face_coordinates', 'mesh2d_face_x mesh2d_face_y')) 
   !
   NF90(nf90_def_var(map_file%ncid, 'mesh2d_face_nodes', NF90_INT, (/map_file%mesh2d_nMax_face_nodes_dimid,map_file%mesh2d_nFaces_dimid/), map_file%mesh2d_face_nodes_varid)) ! location of zb, zs etc. in cell centre
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'cf_role', 'face_node_connectivity'))         
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'mesh', 'mesh2d'))   
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'location', 'face'))   
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'long_name', 'Mapping from every face to its corner nodes (counterclockwise)'))   
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'start_index', 1))         
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, '_FillValue', -999))         
   !   
   if (sferic==0) then
      NF90(nf90_def_var(map_file%ncid, 'mesh2d_node_x', NF90_FLOAT, (/map_file%mesh2d_nNodes_dimid/), map_file%mesh2d_node_x_varid)) ! location of zb, zs etc. in cell centre
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, '_FillValue', FILL_VALUE))         
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'standard_name', 'projection_x_coordinate'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'long_name', 'x-coordinate of mesh nodes'))   
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'location', 'node'))   
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'mesh', 'mesh2d'))   
      !
      NF90(nf90_def_var(map_file%ncid, 'mesh2d_node_y', NF90_FLOAT, (/map_file%mesh2d_nNodes_dimid/), map_file%mesh2d_node_y_varid)) ! location of zb, zs etc. in cell centre
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, '_FillValue', FILL_VALUE))         
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'units', 'm'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'standard_name', 'projection_y_coordinate'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'long_name', 'y-coordinate of mesh nodes'))   
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'location', 'node'))   
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'mesh', 'mesh2d'))   
      !
   else
      NF90(nf90_def_var(map_file%ncid, 'mesh2d_node_x', NF90_FLOAT, (/map_file%mesh2d_nNodes_dimid/), map_file%mesh2d_node_x_varid)) ! location of zb, zs etc. in cell centre
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, '_FillValue', FILL_VALUE))         
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'units', 'degrees_east'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'standard_name', 'longitude'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'long_name', 'x-coordinate of mesh nodes'))   
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'location', 'node'))   
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'mesh', 'mesh2d'))   
      !
      NF90(nf90_def_var(map_file%ncid, 'mesh2d_node_y', NF90_FLOAT, (/map_file%mesh2d_nNodes_dimid/), map_file%mesh2d_node_y_varid)) ! location of zb, zs etc. in cell centre
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, '_FillValue', FILL_VALUE))         
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'units', 'degrees_north'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'standard_name', 'latitude'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'long_name', 'y-coordinate of mesh nodes'))   
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'location', 'node'))   
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'mesh', 'mesh2d'))   
      !
   endif
   NF90(nf90_def_var(map_file%ncid, 'mesh2d_node_z', NF90_FLOAT, (/map_file%mesh2d_nNodes_dimid/), map_file%mesh2d_node_z_varid)) ! location of zb, zs etc. in cell centre
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_z_varid, '_FillValue', FILL_VALUE))         
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_z_varid, 'units', 'm'))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_z_varid, 'standard_name', 'projection_z_coordinate'))
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_z_varid, 'long_name', 'z-coordinate of mesh nodes'))   
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_z_varid, 'location', 'node'))   
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_z_varid, 'mesh', 'mesh2d'))   
   !
   NF90(nf90_def_var(map_file%ncid, 'mesh2d_edge_nodes', NF90_INT, (/map_file%mesh2d_nEdges_dimid, map_file%mesh2d_nMax_edge_nodes_dimid/), map_file%mesh2d_edge_nodes_varid)) ! location of zb, zs etc. in cell centre
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_edge_nodes_varid, 'cf_role', 'edge_node_connectivity'))         
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_edge_nodes_varid, 'mesh', 'mesh2d'))   
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_edge_nodes_varid, 'location', 'edge'))   
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_edge_nodes_varid, 'long_name', 'Mapping from every edge to its adjacent cells '))   
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_edge_nodes_varid, 'start_index', 1))         
   NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_edge_nodes_varid, '_FillValue', -999))         
   !
   NF90(nf90_def_var(map_file%ncid, 'crs', NF90_INT, map_file%crs_varid)) ! For EPSG code
   NF90(nf90_put_att(map_file%ncid, map_file%crs_varid, 'EPSG', '-'))      
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
   NF90(nf90_def_var(map_file%ncid, 'hm0', NF90_FLOAT, (/map_file%mesh2d_nNodes_dimid, map_file%time_dimid/), map_file%hm0_varid)) ! time-varying wave height map
   NF90(nf90_put_att(map_file%ncid, map_file%hm0_varid, '_FillValue', FILL_VALUE))
   NF90(nf90_put_att(map_file%ncid, map_file%hm0_varid, 'units', 'm'))
   NF90(nf90_put_att(map_file%ncid, map_file%hm0_varid, 'standard_name', 'significant_wave_height')) 
   NF90(nf90_put_att(map_file%ncid, map_file%hm0_varid, 'long_name', 'Wave height Hm0'))  
   !
   NF90(nf90_def_var(map_file%ncid, 'hm0_ig', NF90_FLOAT, (/map_file%mesh2d_nNodes_dimid, map_file%time_dimid/), map_file%hm0_ig_varid)) ! time-varying wave height map
   NF90(nf90_put_att(map_file%ncid, map_file%hm0_ig_varid, '_FillValue', FILL_VALUE))
   NF90(nf90_put_att(map_file%ncid, map_file%hm0_ig_varid, 'units', 'm'))
   NF90(nf90_put_att(map_file%ncid, map_file%hm0_ig_varid, 'standard_name', 'significant_infragravity_wave_height')) 
   NF90(nf90_put_att(map_file%ncid, map_file%hm0_ig_varid, 'long_name', 'Infragravity wave height Hm0ig'))  
   !
   NF90(nf90_def_var(map_file%ncid, 'tp', NF90_FLOAT, (/map_file%mesh2d_nNodes_dimid, map_file%time_dimid/), map_file%tp_varid)) ! time-varying wave period map
   NF90(nf90_put_att(map_file%ncid, map_file%tp_varid, '_FillValue', FILL_VALUE))
   NF90(nf90_put_att(map_file%ncid, map_file%tp_varid, 'units', 's'))
   NF90(nf90_put_att(map_file%ncid, map_file%tp_varid, 'standard_name', 'peak period')) 
   NF90(nf90_put_att(map_file%ncid, map_file%tp_varid, 'long_name', 'Peak period Tp'))  
   !
   NF90(nf90_def_var(map_file%ncid, 'wd', NF90_FLOAT, (/map_file%mesh2d_nNodes_dimid, map_file%time_dimid/), map_file%wd_varid)) ! time-varying wave direction map
   NF90(nf90_put_att(map_file%ncid, map_file%wd_varid, '_FillValue', FILL_VALUE))
   NF90(nf90_put_att(map_file%ncid, map_file%wd_varid, 'units', 'deg N'))
   NF90(nf90_put_att(map_file%ncid, map_file%wd_varid, 'standard_name', 'mean wave direction')) 
   NF90(nf90_put_att(map_file%ncid, map_file%wd_varid, 'long_name', 'Mean wave direction'))  
   !
   ! Add for final output:
!   NF90(nf90_def_var(map_file%ncid, 'total_runtime', NF90_FLOAT, (/map_file%runtime_dimid/),map_file%total_runtime_varid))
!   NF90(nf90_put_att(map_file%ncid, map_file%total_runtime_varid, 'units', 's'))   
!   NF90(nf90_put_att(map_file%ncid, map_file%total_runtime_varid, 'long_name', 'total_model_runtime_in_seconds'))
   !
   ! Finish definitions
   !
   NF90(nf90_enddef(map_file%ncid))
   ! 
   ! Write grid to file
   !
   NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_face_nodes_varid, face_nodes, (/1, 1/))) 
   !  
   NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_node_x_varid, x, (/1/))) 
   NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_node_y_varid, y, (/1/))) 
   NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_node_z_varid, zb, (/1/))) 
   ! 
   ! now for cell corners
!   NF90(nf90_put_var(map_file%ncid, map_file%corner_x_varid, xg(1:nmax - 1, 1:mmax - 1), (/1, 1/))) ! write xz of corners
   !
!   NF90(nf90_put_var(map_file%ncid, map_file%corner_y_varid, yg(1:nmax - 1, 1:mmax - 1), (/1, 1/))) ! write yz of corners   
   !
   ! Write epsg, msk & bed level already to file
   !
!   NF90(nf90_put_var(map_file%ncid, map_file%crs_varid, epsg))
   !
   zsg = 0 ! initialise as inactive points       
   !   
   ! write away intermediate data
   !
   NF90(nf90_sync(map_file%ncid)) !write away intermediate data
   !   
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
    subroutine ncoutput_his_init()
   ! 1. Initialise dimensions/variables/attributes
   ! 2. write grid/msk/zb to file
   !
!   use hurrywave_domain
!   use hurrywave_obspoints
   use snapwave_date
   use snapwave_data   
   !
   implicit none   
   !
   if (nobs==0) then ! If no observation points his file is not created        
        return
   endif
   !
   NF90(nf90_create(his_filename, NF90_CLOBBER, his_file%ncid))
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
   NF90(nf90_def_dim(his_file%ncid, 'pointnamelength', 256, his_file%pointnamelength_dimid)) ! length of station_name per obs point  
   NF90(nf90_def_dim(his_file%ncid, 'runtime', 1, his_file%runtime_dimid)) ! total_runtime, average_dt    
   !     
   ! Some metadata attributes 
   NF90(nf90_put_att(his_file%ncid,nf90_global, "Conventions", "Conventions = 'CF-1.6, SGRID-0.3")) 
   NF90(nf90_put_att(his_file%ncid,nf90_global, "Build-Revision-Date-Netcdf-library", trim(nf90_inq_libvers()))) ! version of netcdf library
   NF90(nf90_put_att(his_file%ncid,nf90_global, "Producer", "SnapWave"))
   !NF90(nf90_put_att(his_file%ncid,nf90_global, "Build-Revision", trim(build_revision))) 
   !NF90(nf90_put_att(his_file%ncid,nf90_global, "Build-Date", trim(build_date)))
   NF90(nf90_put_att(his_file%ncid,nf90_global, "title", "HurryWave his point netcdf output"))     
   !
   ! add input params for reproducability
   !call ncoutput_add_params(his_file%ncid,his_file%inp_varid)   
   !
   ! Create variables
   ! Point identifier
   NF90(nf90_def_var(his_file%ncid, 'station_id', NF90_FLOAT, (/his_file%points_dimid/), his_file%station_id_varid)) 
   !NF90(nf90_put_att(his_file%ncid, his_file%station_id_varid, 'units', '-')) !not wanted in fews
   !
   NF90(nf90_def_var(his_file%ncid, 'station_name', NF90_CHAR, (/his_file%pointnamelength_dimid, his_file%points_dimid/), his_file%station_name_varid))
   !NF90(nf90_put_att(his_file%ncid, his_file%station_name_varid, 'units', '-')) !not wanted in fews   
   !
   ! Domain
   NF90(nf90_def_var(his_file%ncid, 'station_x', NF90_FLOAT, (/his_file%points_dimid/), his_file%station_x_varid))   ! non snapped input coordinate 
   NF90(nf90_put_att(his_file%ncid, his_file%station_x_varid, 'units', 'm'))
   NF90(nf90_put_att(his_file%ncid, his_file%station_x_varid, 'standard_name', 'projection_x_coordinate'))
   NF90(nf90_put_att(his_file%ncid, his_file%station_x_varid, 'long_name', 'original_x_coordinate_of_station'))
   NF90(nf90_put_att(his_file%ncid, his_file%station_x_varid, 'grid_mapping', 'crs'))   
   !
   NF90(nf90_def_var(his_file%ncid, 'station_y', NF90_FLOAT, (/his_file%points_dimid/), his_file%station_y_varid)) 
   NF90(nf90_put_att(his_file%ncid, his_file%station_y_varid, 'units', 'm'))
   NF90(nf90_put_att(his_file%ncid, his_file%station_y_varid, 'standard_name', 'projection_y_coordinate'))
   NF90(nf90_put_att(his_file%ncid, his_file%station_y_varid, 'long_name', 'original_y_coordinate_of_station'))
   NF90(nf90_put_att(his_file%ncid, his_file%station_y_varid, 'grid_mapping', 'crs'))   
   !
   NF90(nf90_def_var(his_file%ncid, 'crs', NF90_INT, his_file%crs_varid)) ! For EPSG code
   NF90(nf90_put_att(his_file%ncid, his_file%crs_varid, 'EPSG', '-'))   
   !
   !NF90(nf90_def_var(his_file%ncid, 'point_zb', NF90_FLOAT, (/his_file%points_dimid/), his_file%zb_varid)) ! bed level in cell centre, for points
   !NF90(nf90_put_att(his_file%ncid, his_file%zb_varid, '_FillValue', FILL_VALUE))   
   !NF90(nf90_put_att(his_file%ncid, his_file%zb_varid, 'units', 'm'))
   !NF90(nf90_put_att(his_file%ncid, his_file%zb_varid, 'standard_name', 'altitude'))
   !NF90(nf90_put_att(his_file%ncid, his_file%zb_varid, 'long_name', 'bed_level_above_reference_level'))  
   !
   ! Time variables 
   trefstr_iso8601 = date_to_iso8601(trefstr)
   !  
   NF90(nf90_def_var(his_file%ncid, 'time', NF90_FLOAT, (/his_file%time_dimid/), his_file%time_varid)) ! time
   NF90(nf90_put_att(his_file%ncid, his_file%time_varid, 'units', 'seconds since ' // trim(trefstr_iso8601) ))  ! time stamp following ISO 8601
   NF90(nf90_put_att(his_file%ncid, his_file%time_varid, 'standard_name', 'time'))     
   NF90(nf90_put_att(his_file%ncid, his_file%time_varid, 'long_name', 'time_in_seconds_since_' // trim(trefstr_iso8601) ))   !--> add  trefstr from hurrywave_input
   !
   ! Time varying his output
   NF90(nf90_def_var(his_file%ncid, 'point_hm0', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%hm0_varid)) ! time-varying hm0
   NF90(nf90_put_att(his_file%ncid, his_file%hm0_varid, '_FillValue', FILL_VALUE))
   NF90(nf90_put_att(his_file%ncid, his_file%hm0_varid, 'units', 'm'))
   NF90(nf90_put_att(his_file%ncid, his_file%hm0_varid, 'standard_name', 'sea_surface_wave_significant_height')) 
   NF90(nf90_put_att(his_file%ncid, his_file%hm0_varid, 'long_name', 'Significant wave height Hm0 (m)'))  
   !   
   NF90(nf90_def_var(his_file%ncid, 'point_tp', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%tp_varid)) ! time-varying tp
   NF90(nf90_put_att(his_file%ncid, his_file%tp_varid, '_FillValue', FILL_VALUE))
   NF90(nf90_put_att(his_file%ncid, his_file%tp_varid, 'units', 's'))
   NF90(nf90_put_att(his_file%ncid, his_file%tp_varid, 'standard_name', 'sea_surface_wave_period_at_variance_spectral_density_maximum')) 
   NF90(nf90_put_att(his_file%ncid, his_file%tp_varid, 'long_name', 'Peak wave period Tp (s)'))     
   !   
   NF90(nf90_def_var(his_file%ncid, 'point_wavdir', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%wavdir_varid)) ! time-varying wavdir
   NF90(nf90_put_att(his_file%ncid, his_file%wavdir_varid, '_FillValue', FILL_VALUE))
   NF90(nf90_put_att(his_file%ncid, his_file%wavdir_varid, 'units', 'degree'))
   NF90(nf90_put_att(his_file%ncid, his_file%wavdir_varid, 'standard_name', 'sea_surface_wave_from_direction_at_variance_spectral_density_maximum')) 
   NF90(nf90_put_att(his_file%ncid, his_file%wavdir_varid, 'long_name', 'Peak wave direction (deg, nautical)'))    ! indeed peak wave dir?
   !   
   !NF90(nf90_def_var(his_file%ncid, 'point_dirspr', NF90_FLOAT, (/his_file%points_dimid, his_file%time_dimid/), his_file%dirspr_varid)) ! time-varying wavdir
   !NF90(nf90_put_att(his_file%ncid, his_file%dirspr_varid, '_FillValue', FILL_VALUE))
   !NF90(nf90_put_att(his_file%ncid, his_file%dirspr_varid, 'units', 'degree'))
   !NF90(nf90_put_att(his_file%ncid, his_file%dirspr_varid, 'standard_name', 'sea_surface_wave_directional_spread')) 
   !NF90(nf90_put_att(his_file%ncid, his_file%dirspr_varid, 'long_name', 'Wave directional spread (deg)'))  
  
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
   !NF90(nf90_put_var(his_file%ncid, his_file%crs_varid, crs_epsg))  ! write epsg
   !
   !NF90(nf90_put_var(his_file%ncid, his_file%station_id_varid, idobs))  ! write station_id   
   !   
   NF90(nf90_put_var(his_file%ncid, his_file%station_name_varid, nameobs))  ! write station_name      ! , (/1, nobs/)
   !   
   NF90(nf90_put_var(his_file%ncid, his_file%station_x_varid, xobs)) ! write station_x, input xobs
   !    
   NF90(nf90_put_var(his_file%ncid, his_file%station_y_varid, yobs)) ! write station_y, input yobs
   !   
   !NF90(nf90_put_var(his_file%ncid, his_file%point_x_varid, xgobs)) ! write point_x, now actual value on grid is written rather than input xobs
   !    
   NF90(nf90_sync(his_file%ncid)) !write away intermediate data
   !
   end subroutine
   !
   subroutine ncoutput_update_map(t,ntmapout)
      !
      ! Write time, zs, u, v  
      !
      use snapwave_data   
      use snapwave_results
      !
      implicit none   
      !
      real*8                       :: t  
      !
      integer  :: ntmapout       
      !
      integer  :: k
      !
      !
      NF90(nf90_put_var(map_file%ncid, map_file%time_varid, t, (/ntmapout/))) ! write time
      !
      buf = H*sqrt(2.0)
      !where (depth<0.1) buf=-999.
      NF90(nf90_put_var(map_file%ncid, map_file%hm0_varid, buf, (/1, ntmapout/))) ! write Hm0
      !
      buf = H_ig*sqrt(2.0)
      where (depth<0.1) buf = -999.0
      NF90(nf90_put_var(map_file%ncid, map_file%hm0_ig_varid, buf, (/1, ntmapout/))) ! write Hm0_ig
      !           
      buf = tpmean_bwv*H**0
      where (depth<0.1) buf=-999.0
      NF90(nf90_put_var(map_file%ncid, map_file%tp_varid, buf, (/1, ntmapout/))) ! write Tp
      ! 
!      buf = amod(270.0 - thetam*180/pi + 360.0, 360.0)
      buf = 0.0
      where (depth<0.1) buf = -999.0
      NF90(nf90_put_var(map_file%ncid, map_file%wd_varid, buf, (/1, ntmapout/))) ! write wave direction
      !           
      NF90(nf90_sync(map_file%ncid)) !write away intermediate data ! TL: in first test it seems to be faster to let the file update than keep in memory
      !      
   end subroutine
   !
   !
   subroutine ncoutput_update_his(t, nthisout)
   ! Write time, hm0, tp, wavdir, dirspr of point output
   !
!   use hurrywave_input
!   use hurrywave_domain
!   use hurrywave_obspoints
   use snapwave_data   
   !
   implicit none   
   real*8        :: t
   integer       :: nthisout
   !
   NF90(nf90_put_var(his_file%ncid, his_file%time_varid, t, (/nthisout/))) ! write time
   !   
   NF90(nf90_put_var(his_file%ncid, his_file%hm0_varid, hm0obs, (/1, nthisout/))) ! write point_hm0
   !
   NF90(nf90_put_var(his_file%ncid, his_file%tp_varid, tpmean_bwv*hm0obs**0, (/1, nthisout/))) ! write point_tp  
   !   
   NF90(nf90_put_var(his_file%ncid, his_file%wavdir_varid, wdobs, (/1, nthisout/))) ! write point_wavdir
   !
   !NF90(nf90_put_var(his_file%ncid, his_file%dirspr_varid, dirsprobs, (/1, nthisout/))) ! write point_tp      
   !
   NF90(nf90_sync(his_file%ncid)) !write away intermediate data ! TL: in first test it seems to be faster to let the file update than keep in memory
   !   
   end subroutine
   !
   subroutine ncoutput_map_finalize() 
   ! Add total runtime, dtavg to file and close
   !
   use snapwave_data
   !   
   implicit none   
   !   
!   NF90(nf90_put_var(map_file%ncid, map_file%total_runtime_varid, tfinish_all - tstart_all)) 
!   NF90(nf90_put_var(map_file%ncid, map_file%average_dt_varid,  dtavg)) 
   !
   NF90(nf90_close(map_file%ncid))
   !
   end subroutine

   
!
   subroutine ncoutput_his_finalize()
   ! Add total runtime, dtavg to file and close
   !
   use snapwave_data
   !   
   implicit none   
   !   
   if (nobs==0) then ! If no observation points his file is not created        
        return
   endif   
   !
!   NF90(nf90_put_var(his_file%ncid, his_file%total_runtime_varid, tfinish_all - tstart_all)) 
!   NF90(nf90_put_var(his_file%ncid, his_file%average_dt_varid,  dtavg)) 
   !   
   NF90(nf90_close(his_file%ncid))
   !
   end subroutine
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
   subroutine nc_read_net()
   use snapwave_data
   !
   integer                    :: ierror,j,k
   integer                    :: idfile
   integer                    :: iddim_no_nodes, iddim_no_faces, iddim_nmax_face_nodes, iddim_no_edges, iddim_nmax_edge_nodes
   integer                    :: idvar_node_x, idvar_node_y, idvar_node_z, idvar_face_nodes
   character(NF90_MAX_NAME)   :: string

   !
   ! open grid file
   !
   ierror = nf90_open(gridfile, NF90_NOWRITE, idfile);         call nc_check_err(ierror, "opening file", gridfile)
   !
   !
   ierror = nf90_inq_dimid(idfile, 'mesh2d_nNodes'          , iddim_no_nodes        ); call nc_check_err(ierror, "inq_dimid mesh2d_nNodes"         , gridfile)
   ierror = nf90_inq_dimid(idfile, 'mesh2d_nFaces'          , iddim_no_faces        ); call nc_check_err(ierror, "inq_dimid mesh2d_nFaces"         , gridfile)
   ierror = nf90_inq_dimid(idfile, 'mesh2d_nEdges'          , iddim_no_edges        ); call nc_check_err(ierror, "inq_dimid mesh2d_nEdges"         , gridfile)
   !
   ierror = nf90_inquire_dimension(idfile, iddim_no_nodes        , string, no_nodes        ); call nc_check_err(ierror, "inq_dim nNodes"         , gridfile)
   ierror = nf90_inquire_dimension(idfile, iddim_no_faces        , string, no_faces        ); call nc_check_err(ierror, "inq_dim nFaces"         , gridfile)
   ierror = nf90_inquire_dimension(idfile, iddim_no_edges        , string, no_edges        ); call nc_check_err(ierror, "inq_dim nEdgess"        , gridfile)

   ierror = nf90_inq_varid(idfile, 'mesh2d_node_x'    , idvar_node_x      )                 ; call nc_check_err(ierror, "inq_varid node_x", gridfile)
   ierror = nf90_inq_varid(idfile, 'mesh2d_node_y'    , idvar_node_y      )                 ; call nc_check_err(ierror, "inq_varid node_y", gridfile)
   ierror = nf90_inq_varid(idfile, 'mesh2d_node_z'    , idvar_node_z      )                 ; call nc_check_err(ierror, "inq_varid node_z", gridfile)
   ierror = nf90_inq_varid(idfile, 'mesh2d_face_nodes', idvar_face_nodes  )                 ; call nc_check_err(ierror, "inq_varid face_nodes", gridfile)

    ierror = nf90_get_att(idfile, idvar_node_x,  'standard_name', string); call nc_check_err(ierror, "get_att node_x standard_name", gridfile)
    if (string == 'longitude') then
       sferic = 1
    else
       sferic = 0
    endif
    !
    ! Allocate arrays
    !
    allocate(x(no_nodes))
    allocate(y(no_nodes))
    allocate(xs(no_nodes))
    allocate(ys(no_nodes))
    allocate(zb(no_nodes))
    allocate(msk(no_nodes))
    allocate(face_nodes(4,no_faces))
    allocate(edge_nodes(2,no_edges))
    !
    ierror = nf90_get_var(idfile, idvar_node_x     , x          , start=(/ 1 /)   , count=(/ no_nodes /) ); call nc_check_err(ierror, "get_var x", gridfile)
    ierror = nf90_get_var(idfile, idvar_node_y     , y          , start=(/ 1 /)   , count=(/ no_nodes /) ); call nc_check_err(ierror, "get_var x", gridfile)
    ierror = nf90_get_var(idfile, idvar_node_z     , zb         , start=(/ 1 /)   , count=(/ no_nodes /) ); call nc_check_err(ierror, "get_var z", gridfile)
    ierror = nf90_get_var(idfile, idvar_face_nodes , face_nodes , start=(/ 1 /)   , count=(/ 4, no_faces /) ); call nc_check_err(ierror, "get_var face_nodes", gridfile)
    !
    if (abs(y(1))>90.) sferic = 0  ! Fix to deal with case that network wasconverted to UTM but admin not adjusted
    !
    msk=1
    !
    NF90(nf90_close(idfile))
    
    open(111,file='test.txt')
    write(111,*)no_nodes,no_faces, no_edges
    do k=1,no_nodes
        write(111,'(2e15.8,f11.3,i2)')x(k),y(k),zb(k),msk(k)
    enddo
    do k=1,no_faces
        write(111,'(4i7)')(face_nodes(j,k),j=1,4)
    enddo
    close(111)
        
   
   end subroutine
   !
   subroutine nc_check_err(ierror, description, filename)
!----- GPL ---------------------------------------------------------------------
!                                                                               
!  Copyright (C)  Stichting Deltares, 2011-2019.                                
!                                                                               
!  This program is free software: you can redistribute it and/or modify         
!  it under the terms of the GNU General Public License as published by         
!  the Free Software Foundation version 3.                                      
!                                                                               
!  This program is distributed in the hope that it will be useful,              
!  but WITHOUT ANY WARRANTY; without even the implied warranty of               
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
!  GNU General Public License for more details.                                 
!                                                                               
!  You should have received a copy of the GNU General Public License            
!  along with this program.  If not, see <http://www.gnu.org/licenses/>.        
!                                                                               
!  contact: delft3d.support@deltares.nl                                         
!  Stichting Deltares                                                           
!  P.O. Box 177                                                                 
!  2600 MH Delft, The Netherlands                                               
!                                                                               
!  All indications and logos of, and references to, "Delft3D" and "Deltares"    
!  are registered trademarks of Stichting Deltares, and remain the property of  
!  Stichting Deltares. All rights reserved.                                     
!                                                                               
!-------------------------------------------------------------------------------
!  $Id: nc_check_err.f90 62962 2019-01-16 11:38:35Z mourits $
!  $HeadURL: https://svn.oss.deltares.nl/repos/delft3d/trunk/src/engines_gpl/flow2d3d/packages/data/src/general/nc_check_err.f90 $
!!--description-----------------------------------------------------------------
! NONE
!!--pseudo code and references--------------------------------------------------
! NONE
!!--declarations----------------------------------------------------------------
    use netcdf
    !
    implicit none
!
! Global variables
!
    integer     , intent(in) :: ierror
    character(*), intent(in) :: description
    character(*), intent(in) :: filename
!
! Local variables
!
!
!! executable statements -------------------------------------------------------
!
    if (ierror /= nf90_noerr) then
        write (*,'(6a)') 'ERROR ', trim(description), '. NetCDF file : "', trim(filename), '". Error message:', nf90_strerror(ierror)
    endif
end subroutine nc_check_err
   
end module
