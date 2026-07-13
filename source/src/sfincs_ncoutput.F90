#define NF90(nf90call) call handle_err(nf90call,__FILE__,__LINE__)
!
! ============================================================================
! sfincs_ncoutput — NetCDF output for the SFINCS map file (sfincs_map.nc)
!                   and the his point file (sfincs_his.nc).
!
! Handles regular and quadtree grids through a single set of helpers.
! Caller code does not need separate `if (use_quadtree) … else …` branches
! for standard cell-centered or point-station outputs.
!
! ----------------------------------------------------------------------------
! HOW TO ADD A NEW OUTPUT VARIABLE
! ----------------------------------------------------------------------------
!
! In the recipes below, 'waterlevel' is just a placeholder — replace it with
! the real name of the variable you are adding.
!
! Map-file (cell-centered) output, time-varying — e.g. a new field 'waterlevel':
!   1. Add `integer :: waterlevel_varid` to `map_type` below.
!   2. In ncoutput_map_init, define the var:
!         call def_time_cell_float('waterlevel', map_file%waterlevel_varid, &
!              'units', 'long_name', standard_name='cf_standard_name')
!   3. In ncoutput_update_map, write it each timestep:
!         call write_cell_var(map_file%ncid, map_file%waterlevel_varid, &
!              waterlevel, ntmapout)
!
! Map-file output, max-aggregated — e.g. 'waterlevel_max' (zsmax / hmax / vmax style):
!   1. Add `integer :: waterlevel_max_varid` to `map_type`.
!   2. In ncoutput_map_init:
!         call def_maxtime_cell_float('waterlevel_max', map_file%waterlevel_max_varid, &
!              'units', 'long_name', cell_methods='time: maximum')
!   3. In ncoutput_update_max:
!         call write_cell_var(map_file%ncid, map_file%waterlevel_max_varid, &
!              waterlevel_max, ntmaxout, check_kcs=.true.)
!
! Map-file output, static (single value per cell, written once) — e.g. 'soil':
!   1. Add `integer :: soil_varid` to `map_type`.
!   2. In ncoutput_map_init:
!         call def_static_cell_float('soil', map_file%soil_varid, 'units', &
!              'long_name', standard_name='...')
!   3. In ncoutput_map_init's static-write block (after nf90_enddef):
!         call put_static_cell_float(map_file%ncid, map_file%soil_varid, soil, FILL_VALUE)
!
! Map-file output, static integer mask (stored as integer on quadtree,
! float on regular grid) — e.g. a 'valid_cell' flag:
!   1. Add `integer :: valid_cell_varid` to `map_type`.
!   2. In ncoutput_map_init:
!         call def_static_cell_int('valid_cell', map_file%valid_cell_varid, &
!              'long_name', units='1', standard_name='valid_cell_mask', &
!              description='inactive=0, active=1')
!   3. In ncoutput_map_init's static-write block:
!         call put_static_cell_mask(map_file%ncid, map_file%valid_cell_varid, &
!              real(valid_cell, 4))    ! cast int*1/int*4 source to real*4
!
! His-file (point/station) output, time-varying — e.g. 'point_waterlevel':
!   1. Add `integer :: waterlevel_varid` to `his_type`.
!   2. In ncoutput_his_init:
!         call def_time_point_float('point_waterlevel', his_file%waterlevel_varid, &
!              'units', 'long_name', standard_name='...')
!   3. In ncoutput_update_his: call write_point_var(his_file%waterlevel_varid,
!      source_array, nthisout) — it gathers nmindobs and writes for you.
!      Use optional scale= for unit conversion (e.g. scale=3600000.0 for mm/hr).
!
! Conditions (subgrid, snapwave, infiltration, etc.) belong in the caller's
! `if (...)` guard around the def + write pair, not inside helpers.
!
! GPU note: when running on GPU (OpenACC), the source array of any new
! time-varying output must be synced back to the host before the write call
! using `!$acc update host(<source>)`.  The helpers themselves run on the
! host; they only read host-resident data.
!
! ----------------------------------------------------------------------------
! HELPER INVENTORY
! ----------------------------------------------------------------------------
!
! Generic NetCDF (bottom of helpers file):
!   ncdef_float_var, ncdef_int_var          one-call var-def + attribute set
!   logical2int(lgc)                        returns 1/.true. or 0/.false.
!   handle_err                              NF90 macro error handler
!
! Module-level static-cell writers:
!   put_static_cell_float(ncid, varid, source, fill, [scale, sw_index, min_value])
!   put_static_cell_mask (ncid, varid, source, [sw_index])
!
! Module-level time-varying cell writers:
!   write_cell_var      (ncid, varid, source, nt, [use_sw_index, check_kcs,
!                                                  scale, min_value])
!   write_cell_var_wet  (ncid, varid, source, zref, nt, [check_wet])
!   write_cell_var_depth(ncid, varid, water_level, bed_level, nt, [check_wet])
!
! Definition wrappers used inside ncoutput_map_init (these reuse the
! enclosing scope's dims_* / coord_str / crsgeo / nc_deflate_level so call
! sites stay one-liners):
!   def_static_cell_float / def_static_cell_int    static cell variables
!   def_time_cell_float                            time-varying cell variables
!   def_maxtime_cell_float                         max-aggregated cell variables
!   add_ugrid_face_attrs      attach mesh2d_face_face_link attrs to a varid
!   def_mesh2d_node_coord     UGRID node coord (float for geographic, double for projected)
!   def_grid_axis_coord       SGRID face/corner coord (always float)
!   put_2d                    nf90_put_var with (/1, 1/) start
!
! Definition wrappers used inside ncoutput_his_init:
!   def_time_point_float                           (points × time)
!   def_his_point_coord                            station coordinate variable
!
! His-update writer (gathers nmindobs internally):
!   write_point_var(varid, source, nt, [scale])
!
! Precompute helpers used inside ncoutput_update_his:
!   compute_uv_at_obs_points   face-averaged, rotated (u,v), magnitude, direction
!   compute_wind_at_obs_points speed + meteorological direction (270 convention)
!
! Precompute helpers used inside ncoutput_update_map (return arrays
! indexed by SFINCS cell number):
!   compute_uv_at_cell_centers
!   compute_pnh_unwrapped
!   compute_subgrid_mean_depth  (also used by ncoutput_update_max)
!
! Finalize helpers (write once at end of simulation):
!   ncoutput_write_timestep_analysis
!   ncoutput_write_tsunami_arrival_time
!
! ----------------------------------------------------------------------------
! GRID-TYPE BRANCH STILL NEEDED WHEN
! ----------------------------------------------------------------------------
!
! - the variable is grid-specific topology (UGRID mesh2d / face_node_connectivity
!   on quadtree vs SGRID face_x / corner_x / sfincsgrid on regular) and the
!   Conventions string that goes with it;
! - the variable itself is grid-agnostic (SnapWave output now reads the
!   snapwave_* node arrays via use_sw_index on both quadtree and regular).
!
! Model-physics asymmetries that propagate into output (not a netCDF concern):
! - dynamic bed level (`store_dynamic_bed_level`) only updates on regular
!   non-subgrid runs, so zb is time-varying there and static everywhere else.
!
! ============================================================================
module sfincs_ncoutput
   !
   use netcdf
   use sfincs_ncoutput_helpers
   !
   implicit none
   !

contains

   subroutine ncoutput_map_init()
   !
   ! Merged init for both regular and quadtree grids.
   ! 1. Initialise dimensions/variables/attributes
   ! 2. Write grid/msk/zb to file
   !
   use sfincs_date
   use sfincs_data
   use sfincs_snapwave
   use quadtree
   !
   implicit none
   !
   integer :: nm, nmq, n, m, ntmx, nn, n_nodes, n_faces, iref
   real*4  :: dxx, dyy
   !
   ! sp_dimids is local scratch; nsd/dims_*/coord_str are module-level in
   ! sfincs_ncoutput_helpers and are set here before any def_* calls.
   integer :: sp_dimids(2)
   !
   ! Regular-grid coordinate arrays (xz/yz at face centres, xg/yg at corners)
   real*4, dimension(:,:), allocatable :: xz, yz, xg, yg
   !
   ! Quadtree topology arrays
   real,      dimension(:),   allocatable :: nodes_x, nodes_y
   integer*4, dimension(:,:), allocatable :: face_nodes
   !
   ! Pre-computed source for subgridslope (crsgeo branch lives in the expression)
   real*4,    dimension(:),   allocatable :: slope_buf
   !
   ! -------------------------------------------------------
   ! Setup dimension abstraction
   ! -------------------------------------------------------
   if (use_quadtree) then
      nsd = 1
      coord_str = ''
   else
      nsd = 2
      coord_str = 'x y'
   endif
   !
   ! -------------------------------------------------------
   ! Quadtree: build node/face arrays before creating file
   ! -------------------------------------------------------
   if (use_quadtree) then
      !
      n_faces = quadtree_nr_points
      n_nodes = quadtree_nr_points * 4
      !
      allocate(nodes_x(n_nodes))
      allocate(nodes_y(n_nodes))
      allocate(face_nodes(4, n_faces))
      !
      nodes_x    = 0.0
      nodes_y    = 0.0
      face_nodes = 0
      nn         = 0
      !
      do nmq = 1, quadtree_nr_points
         n    = quadtree_n(nmq)
         m    = quadtree_m(nmq)
         iref = quadtree_level(nmq)
         dxx  = quadtree_dxr(iref)
         dyy  = quadtree_dyr(iref)
         !
         nn = nn + 1
         nodes_x(nn) = x0 + cosrot*(m-1)*dxx - sinrot*(n-1)*dyy
         nodes_y(nn) = y0 + sinrot*(m-1)*dxx + cosrot*(n-1)*dyy
         face_nodes(1, nmq) = nn
         !
         nn = nn + 1
         nodes_x(nn) = x0 + cosrot*(m  )*dxx - sinrot*(n-1)*dyy
         nodes_y(nn) = y0 + sinrot*(m  )*dxx + cosrot*(n-1)*dyy
         face_nodes(2, nmq) = nn
         !
         nn = nn + 1
         nodes_x(nn) = x0 + cosrot*(m  )*dxx - sinrot*(n  )*dyy
         nodes_y(nn) = y0 + sinrot*(m  )*dxx + cosrot*(n  )*dyy
         face_nodes(3, nmq) = nn
         !
         nn = nn + 1
         nodes_x(nn) = x0 + cosrot*(m-1)*dxx - sinrot*(n  )*dyy
         nodes_y(nn) = y0 + sinrot*(m-1)*dxx + cosrot*(n  )*dyy
         face_nodes(4, nmq) = nn
      enddo
      !
   endif
   !
   ! -------------------------------------------------------
   ! Create NetCDF file
   ! -------------------------------------------------------
   NF90(nf90_create('sfincs_map.nc', ior(NF90_CLOBBER, NF90_NETCDF4), map_file%ncid))
   !
   ! Create dimensions
   !
   if (use_quadtree) then
      !
      NF90(nf90_def_dim(map_file%ncid, 'nmesh2d_node',           n_nodes, map_file%nmesh2d_node_dimid))
      NF90(nf90_def_dim(map_file%ncid, 'nmesh2d_face',           n_faces, map_file%nmesh2d_face_dimid))
      NF90(nf90_def_dim(map_file%ncid, 'max_nmesh2d_face_nodes', 4,       map_file%max_nmesh2d_face_nodes_dimid))
      !
      sp_dimids(1) = map_file%nmesh2d_face_dimid
      !
   else
      !
      NF90(nf90_def_dim(map_file%ncid, 'n',        nmax,     map_file%n_dimid))
      NF90(nf90_def_dim(map_file%ncid, 'm',        mmax,     map_file%m_dimid))
      NF90(nf90_def_dim(map_file%ncid, 'corner_n', nmax + 1, map_file%corner_n_dimid))
      NF90(nf90_def_dim(map_file%ncid, 'corner_m', mmax + 1, map_file%corner_m_dimid))
      !
      sp_dimids(1) = map_file%m_dimid
      sp_dimids(2) = map_file%n_dimid
      !
   endif
   !
   NF90(nf90_def_dim(map_file%ncid, 'time',    NF90_UNLIMITED, map_file%time_dimid))
   ntmx = max(ceiling((t1out - t0out)/dtmaxout), 1)
   NF90(nf90_def_dim(map_file%ncid, 'timemax', ntmx, map_file%timemax_dimid))
   NF90(nf90_def_dim(map_file%ncid, 'runtime', 1,    map_file%runtime_dimid))
   !
   if (store_vegetation) then
      NF90(nf90_def_dim(map_file%ncid, 'nsec', vegetation_vertical_segments, map_file%nsec_dimid)) ! number of vegetation vertical sections
   endif
   !
   ! Build dim arrays
   dims_s (1:nsd)   = sp_dimids(1:nsd)
   dims_st(1:nsd)   = sp_dimids(1:nsd);  dims_st(nsd+1) = map_file%time_dimid
   dims_sm(1:nsd)   = sp_dimids(1:nsd);  dims_sm(nsd+1) = map_file%timemax_dimid
   !
   ! Global metadata
   ! CF version unified at 1.8 across both grid types and the his file.
   if (use_quadtree) then
      NF90(nf90_put_att(map_file%ncid, nf90_global, "Conventions", "CF-1.8 UGRID-1.0 Deltares-0.10"))
   else
      NF90(nf90_put_att(map_file%ncid, nf90_global, "Conventions", "CF-1.8 SGRID-0.3"))
   endif
   NF90(nf90_put_att(map_file%ncid, nf90_global, "Build-Revision-Date-Netcdf-library", trim(nf90_inq_libvers())))
   NF90(nf90_put_att(map_file%ncid, nf90_global, "Producer", "SFINCS model: Super-Fast INundation of CoastS"))
   NF90(nf90_put_att(map_file%ncid, nf90_global, "Build-Revision", trim(build_revision)))
   NF90(nf90_put_att(map_file%ncid, nf90_global, "Build-Date", trim(build_date)))
   NF90(nf90_put_att(map_file%ncid, nf90_global, "title", "SFINCS map netcdf output"))
   !
   call ncoutput_add_params(map_file%ncid, map_file%inp_varid)
   !
   ! -------------------------------------------------------
   ! Grid topology variables (grid-type specific)
   ! -------------------------------------------------------
   if (use_quadtree) then
      !
      NF90(nf90_def_var(map_file%ncid, 'mesh2d', NF90_INT, (/map_file%runtime_dimid/), map_file%mesh2d_varid))
      NF90(nf90_def_var_deflate(map_file%ncid, map_file%mesh2d_varid, 1, 1, nc_deflate_level))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'cf_role',                   'mesh_topology'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'long_name',                 'Topology data of 2D network'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'topology_dimension',        2))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'node_coordinates',          'mesh2d_node_x mesh2d_node_y'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'node_dimension',            'nmesh2d_node'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'max_face_nodes_dimension',  'max_nmesh2d_face_nodes'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'face_node_connectivity',    'mesh2d_face_nodes'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'face_dimension',            'nmesh2d_face'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_varid, 'face_coordinates',          'mesh2d_face_x mesh2d_face_y'))
      !
      call def_mesh2d_node_coord('x', map_file%mesh2d_node_x_varid)
      !
      call def_mesh2d_node_coord('y', map_file%mesh2d_node_y_varid)
      !
      NF90(nf90_def_var(map_file%ncid, 'mesh2d_face_nodes', NF90_INT, (/map_file%max_nmesh2d_face_nodes_dimid, map_file%nmesh2d_face_dimid/), map_file%mesh2d_face_nodes_varid))
      NF90(nf90_def_var_deflate(map_file%ncid, map_file%mesh2d_face_nodes_varid, 1, 1, nc_deflate_level))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'cf_role',   'face_node_connectivity'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'mesh',      'mesh2d'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'location',  'face'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'long_name', 'Mapping from every face to its corner nodes (counterclockwise)'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'start_index', 1))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, '_FillValue', -999))
      !
      ! Face centroid coordinates (required by UGRID and MDAL/QGIS)
      if (crsgeo) then
         NF90(nf90_def_var(map_file%ncid, 'mesh2d_face_x', NF90_FLOAT, (/map_file%nmesh2d_face_dimid/), map_file%mesh2d_face_x_varid))
         NF90(nf90_def_var_deflate(map_file%ncid, map_file%mesh2d_face_x_varid, 1, 1, nc_deflate_level))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_x_varid, 'units',         'degrees_east'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_x_varid, 'standard_name', 'longitude'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_x_varid, 'long_name',     'Characteristic longitude of mesh face'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_x_varid, 'grid_mapping',  'crs'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_x_varid, 'mesh',          'mesh2d'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_x_varid, 'location',      'face'))
         !
         NF90(nf90_def_var(map_file%ncid, 'mesh2d_face_y', NF90_FLOAT, (/map_file%nmesh2d_face_dimid/), map_file%mesh2d_face_y_varid))
         NF90(nf90_def_var_deflate(map_file%ncid, map_file%mesh2d_face_y_varid, 1, 1, nc_deflate_level))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_y_varid, 'units',         'degrees_north'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_y_varid, 'standard_name', 'latitude'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_y_varid, 'long_name',     'Characteristic latitude of mesh face'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_y_varid, 'grid_mapping',  'crs'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_y_varid, 'mesh',          'mesh2d'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_y_varid, 'location',      'face'))
      else
         NF90(nf90_def_var(map_file%ncid, 'mesh2d_face_x', NF90_DOUBLE, (/map_file%nmesh2d_face_dimid/), map_file%mesh2d_face_x_varid))
         NF90(nf90_def_var_deflate(map_file%ncid, map_file%mesh2d_face_x_varid, 1, 1, nc_deflate_level))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_x_varid, 'units',         'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_x_varid, 'standard_name', 'projection_x_coordinate'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_x_varid, 'long_name',     'Characteristic x-coordinate of mesh face'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_x_varid, 'grid_mapping',  'crs'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_x_varid, 'mesh',          'mesh2d'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_x_varid, 'location',      'face'))
         !
         NF90(nf90_def_var(map_file%ncid, 'mesh2d_face_y', NF90_DOUBLE, (/map_file%nmesh2d_face_dimid/), map_file%mesh2d_face_y_varid))
         NF90(nf90_def_var_deflate(map_file%ncid, map_file%mesh2d_face_y_varid, 1, 1, nc_deflate_level))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_y_varid, 'units',         'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_y_varid, 'standard_name', 'projection_y_coordinate'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_y_varid, 'long_name',     'Characteristic y-coordinate of mesh face'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_y_varid, 'grid_mapping',  'crs'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_y_varid, 'mesh',          'mesh2d'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_y_varid, 'location',      'face'))
      endif
      !
      NF90(nf90_def_var(map_file%ncid, 'crs', NF90_INT, map_file%crs_varid))
      NF90(nf90_put_att(map_file%ncid, map_file%crs_varid, 'epsg',      epsg))
      NF90(nf90_put_att(map_file%ncid, map_file%crs_varid, 'epsg_code', 'EPSG:' // trim(epsg_code)))
      !
   else
      !
      ! Regular: face / corner coordinate axes
      call def_grid_axis_coord('x',        'x', (/map_file%m_dimid,        map_file%n_dimid/),        map_file%face_x_varid,   'face_x')
      call def_grid_axis_coord('y',        'y', (/map_file%m_dimid,        map_file%n_dimid/),        map_file%face_y_varid,   'face_y')
      call def_grid_axis_coord('corner_x', 'x', (/map_file%corner_m_dimid, map_file%corner_n_dimid/), map_file%corner_x_varid, 'corner_x')
      call def_grid_axis_coord('corner_y', 'y', (/map_file%corner_m_dimid, map_file%corner_n_dimid/), map_file%corner_y_varid, 'corner_y')
      !
      NF90(nf90_def_var(map_file%ncid, 'crs', NF90_INT, map_file%crs_varid))
      NF90(nf90_put_att(map_file%ncid, map_file%crs_varid, 'epsg',      epsg))
      NF90(nf90_put_att(map_file%ncid, map_file%crs_varid, 'epsg_code', 'EPSG:' // trim(epsg_code)))
      !
      NF90(nf90_def_var(map_file%ncid, 'sfincsgrid', NF90_INT, map_file%grid_varid))
      NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'cf_role',           'grid_topology'))
      NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'topology_dimension', 2))
      NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'node_dimensions',   'n m'))
      NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'face_dimensions',   'n: m:'))
      NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'corner_dimensions', 'corner_n: corner_m:'))
      NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'face_coordinates',  'x y'))
      NF90(nf90_put_att(map_file%ncid, map_file%grid_varid, 'corner_coordinates','corner_x corner_y'))
      !
   endif
   !
   ! -------------------------------------------------------
   ! msk: NF90_INT on both grid types, described via CF flag_values /
   ! flag_meanings. No CF standard_name exists for multi-valued masks
   ! (the previously-used 'land_binary_mask' is strict 0/1 only).
   ! -------------------------------------------------------
   call def_static_cell_int('msk', map_file%msk_varid, 'Active cells mask', &
        units='-', &
        description='inactive=0, active=1, normal_boundary=2, outflow_boundary=3, wavemaker=4, downstream_boundary=5, neumann_boundary=6', &
        flag_values=(/0, 1, 2, 3, 4, 5, 6/), &
        flag_meanings='inactive active normal_boundary outflow_boundary wavemaker downstream_boundary neumann_boundary')
   !
   ! -------------------------------------------------------
   ! Infiltration map (qinf) - both grid types, slightly different
   ! -------------------------------------------------------
   if (infiltration) then
      if (inftype == 'cna') then
         call def_static_cell_float('qinf', map_file%qinf_varid, 'm', 'moisture storage (S) capacity - Curve number', &
              standard_name='S')
      elseif (inftype == 'cnb') then
         call def_static_cell_float('qinf', map_file%qinf_varid, 'm', 'maximum moisture storage (Smax) capacity - Curve number', &
              standard_name='Smax')
      elseif (inftype == 'gai') then
         call def_static_cell_float('qinf', map_file%qinf_varid, 'm', 'suction head at the wetting front - Green and Ampt', &
              standard_name='psi')
      elseif (inftype == 'hor') then
         call def_static_cell_float('qinf', map_file%qinf_varid, 'm', 'initial infiltration rate - Horton', standard_name='f0')
      else
         call def_static_cell_float('qinf', map_file%qinf_varid, 'mm h-1', 'infiltration rate - constant in time', &
              standard_name='qinf')
      endif
   endif
   !
   ! -------------------------------------------------------
   ! zb: time-varying only on regular non-subgrid runs with store_dynamic_bed_level;
   ! static for quadtree, subgrid, or static-bed regular runs.
   ! Def condition matches the write condition in ncoutput_update_map.
   ! -------------------------------------------------------
   if (.not. use_quadtree .and. store_dynamic_bed_level .and. .not. subgrid) then
      call def_time_cell_float('zb', map_file%zb_varid, 'm', 'Bed level above reference level', standard_name='altitude')
   else
      call def_static_cell_float('zb', map_file%zb_varid, 'm', 'Bed level above reference level', standard_name='altitude')
   endif
   !
   ! manning (only meaningful when manning2d, i.e. per-cell field, is set)
   if (.not. subgrid .and. manning2d) then
      call def_static_cell_float('manning', map_file%manning_varid, 's/m^1/3', 'Manning roughness coefficient', standard_name='manning')
   endif
   !
   ! subgrid slope
   if (subgrid .and. store_hsubgrid .and. store_hmean) then
      call def_static_cell_float('subgridslope', map_file%subgridslope_varid, '-', 'Subgrid slope', standard_name='subgrid_slope')
   endif
   !
   ! vegetation stem properties (static, per vertical section)
   if (store_vegetation) then
      call def_static_veg_float('vegetation_stems_cd',       map_file%veg_cd_varid,     '-',   'Bulk drag coefficient per vegetation section',         map_file%nsec_dimid, standard_name='vegetation_stems_cd')
      call def_static_veg_float('vegetation_stems_height',   map_file%veg_ah_varid,     'm',   'Vegetation section thickness',                        map_file%nsec_dimid, standard_name='vegetation_stems_height')
      call def_static_veg_float('vegetation_stems_diameter', map_file%veg_bstems_varid, 'm',   'Diameter of individual vegetation stems per section', map_file%nsec_dimid, standard_name='vegetation_stems_diameter')
      call def_static_veg_float('vegetation_stems_density',  map_file%veg_Nstems_varid, 'm-2', 'Number of stems per unit horizontal area per section', map_file%nsec_dimid, standard_name='vegetation_stems_density')
   endif
   !
   ! -------------------------------------------------------
   ! Time variables
   ! -------------------------------------------------------
   trefstr_iso8601 = date_to_iso8601(trefstr)
   NF90(nf90_def_var(map_file%ncid, 'time', NF90_FLOAT, (/map_file%time_dimid/), map_file%time_varid))
   NF90(nf90_put_att(map_file%ncid, map_file%time_varid, 'units',     'seconds since ' // trim(trefstr_iso8601)))
   NF90(nf90_put_att(map_file%ncid, map_file%time_varid, 'standard_name', 'time'))
   NF90(nf90_put_att(map_file%ncid, map_file%time_varid, 'long_name', 'time_in_seconds_since_' // trim(trefstr_iso8601)))
   !
   ! -------------------------------------------------------
   ! Time-varying water level / depth / velocity
   ! -------------------------------------------------------
   call def_time_cell_float('zs', map_file%zs_varid, 'm', 'Water level', standard_name='sea_surface_height_above_reference_level')
   !
   if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then
      call def_time_cell_float('h', map_file%h_varid, 'm', 'Water depth', standard_name='water_depth')
   endif
   !
   if (store_velocity) then
      call def_time_cell_float('u', map_file%u_varid, 'm s-1', 'Flow velocity x-component', &
           standard_name='eastward_sea_water_velocity')
      !
      call def_time_cell_float('v', map_file%v_varid, 'm s-1', 'Flow velocity y-component', &
           standard_name='northward_sea_water_velocity')
   endif
   !
   ! Subgrid volumes
   if (subgrid) then
      if (store_zvolume) then
         call def_time_cell_float('subgrid_volume', map_file%zvolume_varid, 'm3', 'Subgrid volume in cell', &
              standard_name='subgrid_volume_in_cell')
      endif
      if (store_storagevolume) then
         call def_time_cell_float('storage_volume', map_file%storagevolume_varid, 'm3', 'Storage volume in cell', &
              standard_name='storage_volume_in_cell')
      endif
   endif
   !
   ! Infiltration state vars (Seff / sigma / f). Source arrays scs_Se,
   ! GA_sigma and qinfmap are allocated unconditionally by sfincs_infiltration
   ! whenever the corresponding inftype is active, on both regular and
   ! quadtree grids.
   if (inftype == 'cnb') then
      call def_time_cell_float('Seff', map_file%infstate_varid, 'm', 'current moisture storage (Se) capacity', standard_name='Se')
   elseif (inftype == 'gai') then
      call def_time_cell_float('sigma', map_file%infstate_varid, '-', 'maximum soil moisture deficit', standard_name='sigma')
   elseif (inftype == 'hor') then
      call def_time_cell_float('f', map_file%infstate_varid, 'mm h-1', 'current infiltration capacity', standard_name='f')
   endif
   !
   ! -------------------------------------------------------
   ! Max-output time variable
   ! -------------------------------------------------------
   if (store_maximum_waterlevel) then
      NF90(nf90_def_var(map_file%ncid, 'timemax', NF90_FLOAT, (/map_file%timemax_dimid/), map_file%timemax_varid))
      NF90(nf90_put_att(map_file%ncid, map_file%timemax_varid, 'units',     'seconds since ' // trim(trefstr_iso8601)))
      NF90(nf90_put_att(map_file%ncid, map_file%timemax_varid, 'standard_name', 'time'))
      NF90(nf90_put_att(map_file%ncid, map_file%timemax_varid, 'long_name', 'time_in_seconds_since_' // trim(trefstr_iso8601)))
   endif
   !
   if (store_maximum_waterlevel) then
      call def_maxtime_cell_float('zsmax', map_file%zsmax_varid, 'm', 'Maximum water level', &
           standard_name='maximum_sea_surface_height_above_reference_level')
   endif
   !
   if (store_cumulative_precipitation) then
      call def_maxtime_cell_float('cumprcp', map_file%cumprcp_varid, 'm', 'Cumulative precipitation depth', &
           standard_name='cumulative_precipitation_depth', cell_methods='time: sum')
   endif
   !
   if (store_twet) then
      call def_maxtime_cell_float('tmax', map_file%tmax_varid, 'seconds', 'Time cell was wet', &
           standard_name='duration_wet', cell_methods='time: sum')
   endif
   !
   if (store_t_zsmax) then
      call def_maxtime_cell_float('t_zsmax', map_file%t_zsmax_varid, 'seconds since ' // trim(trefstr_iso8601), &
           'time when zsmax occurs', standard_name='t_zsmax', cell_methods='time: maximum')
   endif
   !
   if (store_maximum_waterlevel) then
      if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then
         call def_maxtime_cell_float('hmax', map_file%hmax_varid, 'm', 'Maximum water depth', &
              standard_name='sea_floor_depth_below_sea_surface', cell_methods='time: maximum')
      endif
   endif
   !
   if (store_maximum_velocity) then
      call def_maxtime_cell_float('vmax', map_file%vmax_varid, 'm s-1', 'Maximum flow velocity', &
           standard_name='maximum_flow_velocity', cell_methods='time: maximum')
   endif
   !
   if (store_maximum_flux) then
      call def_maxtime_cell_float('qmax', map_file%qmax_varid, 'm^2 s-1', 'maximum_flux', standard_name='maximum_flux', &
           cell_methods='time: maximum')
   endif
   !
   if (timestep_analysis) then
      call def_static_cell_float('average_required_timestep', map_file%average_required_timestep_varid, 's', &
           'time-averaged required timestep over the simulation', standard_name='average_required_timestep')
      !
      call def_static_cell_float('percentage_limiting_timestep', map_file%percentage_limiting_varid, '%', &
           'percentage of simulation steps in which this cell was the limiting cell', &
           standard_name='percentage_limiting_timestep')
   endif
   !
   if (store_cumulative_precipitation .and. infiltration) then
      call def_maxtime_cell_float('cuminf', map_file%cuminf_varid, 'm', 'cumulative_infiltration_depth', cell_methods='time: sum')
   endif
   !
   ! -------------------------------------------------------
   ! Meteo variables
   ! -------------------------------------------------------
   if (store_meteo) then
      if (wind) then
         call def_time_cell_float('wind_u', map_file%wind_u_varid, 'm s-1', 'Wind speed u-component', standard_name='eastward_wind')
         !
         call def_time_cell_float('wind_v', map_file%wind_v_varid, 'm s-1', 'Wind speed v-component', standard_name='northward_wind')
         !
         ! windmax (treated like all other max fields: max-time-varying)
         if (store_wind_max .and. meteo3d) then
            call def_maxtime_cell_float('windmax', map_file%windmax_varid, 'm s-1', 'Maximum wind speed', &
                 cell_methods='time: maximum')
         endif
      endif
      !
      if (patmos) then
         call def_time_cell_float('surface_air_pressure', map_file%patm_varid, 'N m-2', 'Surface air pressure', &
              standard_name='surface_air_pressure')
      endif
      !
      ! precipitation_rate (prcp source array is np-shaped on both grids;
      ! sfincs_meteo applies precip identically regardless of grid type)
      if (precip) then
         call def_time_cell_float('precipitation_rate', map_file%precip_varid, 'mm h-1', 'Precipitation rate', &
              standard_name='precipitation_rate')
      endif
   endif
   !
   ! -------------------------------------------------------
   ! SnapWave variables
   ! -------------------------------------------------------
   if (snapwave) then
      !
      ! snapwavemsk: NF90_INT on both grid types, described via CF
      ! flag_values / flag_meanings.
      call def_static_cell_int('snapwavemsk', map_file%snapwavemsk_varid, 'SnapWave active cells mask', &
           units='-', &
           description='inactive=0, active=1, wave_boundary=2, neumann_boundary=3', &
           flag_values=(/0, 1, 2, 3/), &
           flag_meanings='inactive active wave_boundary neumann_boundary')
      !
      call def_time_cell_float('hm0', map_file%hm0_varid, 'm', 'Hm0 wave height', standard_name='hm0_wave_height')
      !
      call def_time_cell_float('hm0ig', map_file%hm0ig_varid, 'm', 'Hm0 infragravity wave height', &
           standard_name='hm0_ig_wave_height')
      !
      call def_time_cell_float('tp', map_file%tp_varid, 's', 'Peak wave period', standard_name='peak_wave_period')
      !
      call def_time_cell_float('tpig', map_file%tpig_varid, 's', 'Peak infragravity wave period', &
           standard_name='peak_ig_wave_period')
      !
      if (store_wave_forces) then
         call def_time_cell_float('fwx', map_file%fwx_varid, 'm', 'Wave force x-component', standard_name='wave_force_x')
         !
         call def_time_cell_float('fwy', map_file%fwy_varid, 'm', 'Wave force y-component', standard_name='wave_force_y')
         !
         call def_time_cell_float('beta', map_file%beta_varid, '-', 'Mean local bed slope', &
              standard_name='directionally_averaged_local_bed_slope')
         !
         call def_time_cell_float('snapwavedepth', map_file%snapwavedepth_varid, 'm', 'Interpolated water depth in Snapwave', &
              standard_name='snapwave_waterdepth')
      endif
      !
      if (store_wave_direction) then
         call def_time_cell_float('wavdir', map_file%wavdir_varid, 'degrees', 'Mean wave angle (deg)', &
              standard_name='mean_wave_direction')
      endif
      !
      if (wavemaker) then
         call def_time_cell_float('zsm', map_file%zsm_varid, 'm', 'Filtered water level', standard_name='mean_water_level')
      endif
      !
   endif
   !
   ! tsunami_arrival_time (single value per cell, written once at finalize)
   if (store_tsunami_arrival_time) then
      call def_static_cell_float('tsunami_arrival_time', map_file%tsunami_arrival_time_varid, &
           '-', 'Tsunami arrival time', standard_name='tsunami_arrival_time')
   endif
   !
   if (nonhydrostatic) then
      call def_time_cell_float('pnonh', map_file%pnonh_varid, 'N m-2', 'Non-hydrostatic pressure')
   endif
   !
   ! Runtime scalars
   call ncdef_float_var(map_file%ncid, 'total_runtime', (/map_file%runtime_dimid/), map_file%total_runtime_varid, &
        's', 'Total model runtime (s)')
   !
   call ncdef_float_var(map_file%ncid, 'average_dt', (/map_file%runtime_dimid/), map_file%average_dt_varid, &
        's', 'Average model time step (s)')
   !
   call ncdef_float_var(map_file%ncid, 'status', (/map_file%runtime_dimid/), map_file%status_varid, &
        '-', 'status of SFINCS simulation - 0 is no error')
   !
   ! -------------------------------------------------------
   ! Finish definitions
   ! -------------------------------------------------------
   NF90(nf90_enddef(map_file%ncid))
   !
   ! -------------------------------------------------------
   ! Write static grid data
   ! -------------------------------------------------------
   !
   ! Topology / coordinate axes — grid-type-specific
   if (use_quadtree) then
      NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_node_x_varid,     nodes_x))
      NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_node_y_varid,     nodes_y))
      NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_face_nodes_varid, face_nodes))
      ! Write face centroid coordinates (cell centres)
      block
         real, allocatable :: face_cx(:), face_cy(:)
         integer :: ifac, iref2
         real    :: dxx2, dyy2
         allocate(face_cx(n_faces), face_cy(n_faces))
         do ifac = 1, n_faces
            n    = quadtree_n(ifac)
            m    = quadtree_m(ifac)
            iref2 = quadtree_level(ifac)
            dxx2  = quadtree_dxr(iref2)
            dyy2  = quadtree_dyr(iref2)
            face_cx(ifac) = x0 + cosrot*(m - 0.5)*dxx2 - sinrot*(n - 0.5)*dyy2
            face_cy(ifac) = y0 + sinrot*(m - 0.5)*dxx2 + cosrot*(n - 0.5)*dyy2
         enddo
         NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_face_x_varid, face_cx))
         NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_face_y_varid, face_cy))
         deallocate(face_cx, face_cy)
      end block
   else
      allocate(xz(mmax,     nmax))
      allocate(yz(mmax,     nmax))
      allocate(xg(mmax + 1, nmax + 1))
      allocate(yg(mmax + 1, nmax + 1))
      do n = 1, nmax
         do m = 1, mmax
            xz(m, n) = x0 + cosrot*(1.0*(m - 0.5))*dx - sinrot*(1.0*(n - 0.5))*dy
            yz(m, n) = y0 + sinrot*(1.0*(m - 0.5))*dx + cosrot*(1.0*(n - 0.5))*dy
         enddo
      enddo
      do n = 1, nmax + 1
         do m = 1, mmax + 1
            xg(m, n) = x0 + cosrot*(1.0*(m - 1))*dx - sinrot*(1.0*(n - 1))*dy
            yg(m, n) = y0 + sinrot*(1.0*(m - 1))*dx + cosrot*(1.0*(n - 1))*dy
         enddo
      enddo
      call put_2d(map_file%face_x_varid,   xz)
      call put_2d(map_file%face_y_varid,   yz)
      call put_2d(map_file%corner_x_varid, xg)
      call put_2d(map_file%corner_y_varid, yg)
   endif
   NF90(nf90_put_var(map_file%ncid, map_file%crs_varid, epsg))
   !
   ! Cell-data writes — uniform via gather/put helpers
   !
   ! zb static write — fires whenever the def chose the static shape.
   ! That is: quadtree, OR regular without dynamic bed level, OR subgrid
   ! (subgrid_z_zmin is a single value per cell so dynamic bed level does
   ! not apply). The complement (regular + dynamic bed level + non-subgrid)
   ! is written each timestep in ncoutput_update_map instead.
   if (use_quadtree .or. .not. store_dynamic_bed_level .or. subgrid) then
      if (subgrid) then
         call put_static_cell_float(map_file%ncid, map_file%zb_varid, subgrid_z_zmin, FILL_VALUE)
      else
         call put_static_cell_float(map_file%ncid, map_file%zb_varid, zb, FILL_VALUE)
      endif
   endif
   !
   ! subgrid slope (precompute the source — denominator depends on crsgeo)
   if (subgrid .and. store_hsubgrid .and. store_hmean) then
      allocate(slope_buf(np))
      if (crsgeo) then
         do nm = 1, np
            slope_buf(nm) = (subgrid_z_zmax(nm) - subgrid_z_zmin(nm)) / sqrt(cell_area_m2(nm))
         enddo
      else
         do nm = 1, np
            slope_buf(nm) = (subgrid_z_zmax(nm) - subgrid_z_zmin(nm)) / sqrt(cell_area(z_flags_iref(nm)))
         enddo
      endif
      call put_static_cell_float(map_file%ncid, map_file%subgridslope_varid, slope_buf, FILL_VALUE)
      deallocate(slope_buf)
   endif
   !
   ! msk (kcs is integer*1 — cast to real*4 for the generic mask helper)
   call put_static_cell_mask(map_file%ncid, map_file%msk_varid, real(kcs, 4))
   !
   ! snapwave msk (real*4 source; uses snapwave's own quadtree index)
   if (snapwave) then
      call put_static_cell_mask(map_file%ncid, map_file%snapwavemsk_varid, snapwave_mask, sw_index=.true.)
   endif
   !
   ! Manning
   if (.not. subgrid .and. manning2d) then
      call put_static_cell_float(map_file%ncid, map_file%manning_varid, rghfield, FILL_VALUE)
   endif
   !
   ! Infiltration map (cna/c2d are stored as mm h-1, others use raw qinffield units)
   if (infiltration) then
      if (inftype == 'con' .or. inftype == 'c2d') then
         call put_static_cell_float(map_file%ncid, map_file%qinf_varid, qinffield, FILL_VALUE, scale=3.6e6)
      else
         call put_static_cell_float(map_file%ncid, map_file%qinf_varid, qinffield, FILL_VALUE)
      endif
   endif
   !
   ! Vegetation stem properties (static, written once at init)
   if (store_vegetation) then
      call put_static_veg_float(map_file%ncid, map_file%veg_cd_varid,     vegetation_stems_cd,       vegetation_vertical_segments, FILL_VALUE)
      call put_static_veg_float(map_file%ncid, map_file%veg_ah_varid,     vegetation_stems_height,   vegetation_vertical_segments, FILL_VALUE)
      call put_static_veg_float(map_file%ncid, map_file%veg_bstems_varid, vegetation_stems_diameter, vegetation_vertical_segments, FILL_VALUE)
      call put_static_veg_float(map_file%ncid, map_file%veg_Nstems_varid, vegetation_stems_density,  vegetation_vertical_segments, FILL_VALUE)
   endif
   !
   NF90(nf90_sync(map_file%ncid))
   !
   end subroutine ncoutput_map_init


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   subroutine ncoutput_his_init()
   !
   ! 1. Initialise dimensions/variables/attributes
   ! 2. write grid/msk/zb to file
   !
   use sfincs_date
   use sfincs_data
   use sfincs_structures
   use sfincs_src_structures, only: nr_src_structures, src_struc_name, src_struc_type, structure_dike_breach
   use sfincs_discharges,     only: src_name, nr_discharge_points
   use sfincs_urban_drainage, only: nr_urban_drainage_zones, urb_zone_name
   !
   implicit none
   !
   integer                      :: istruc
   !
   real*4, dimension(:,:), allocatable :: struc_info
   real*4, dimension(:), allocatable :: struc_x
   real*4, dimension(:), allocatable :: struc_y
   real*4, dimension(:), allocatable :: struc_height
   !
   real*4, dimension(:,:), allocatable :: thindam_info
   real*4, dimension(:), allocatable :: thindam_x
   real*4, dimension(:), allocatable :: thindam_y
   !
   character*256, dimension(:), allocatable :: drain_name_buf
   character*256, dimension(:), allocatable :: river_name_buf
   character*256, dimension(:), allocatable :: urbdrain_name_buf
   !
   if (nobs==0 .and. nrcrosssections==0 .and. nrstructures==0 .and. nr_src_structures==0 .and. .not. (nr_discharge_points>0 .and. store_river_discharge) .and. .not. (nr_urban_drainage_zones>0 .and. store_urban_drainage_discharge) .and. nr_runup_gauges==0) then ! If no observation points, cross-sections, structures, drains, river sources, urban drainage zones or run-up gauges; his file is not created
      return
   endif
   !
   NF90(nf90_create('sfincs_his.nc', ior(NF90_CLOBBER,NF90_NETCDF4), his_file%ncid))
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
   if (nr_src_structures>0) then
      NF90(nf90_def_dim(his_file%ncid, 'drainage', nr_src_structures, his_file%drain_dimid)) ! nr of drainage structures
   endif
   !
   if (nr_discharge_points>0 .and. store_river_discharge) then
      NF90(nf90_def_dim(his_file%ncid, 'rivers', nr_discharge_points, his_file%river_dimid)) ! nr of river point sources
   endif
   !
   if (nr_urban_drainage_zones > 0 .and. store_urban_drainage_discharge) then
      NF90(nf90_def_dim(his_file%ncid, 'urban_drainage_zones', nr_urban_drainage_zones, his_file%urbdrain_dimid)) ! nr of urban drainage zones
   endif
   !
   if (nrstructures>0) then
      NF90(nf90_def_dim(his_file%ncid, 'structures', nrstructures, his_file%structures_dimid)) ! nr of structures (weir)
   endif   
   !
   if (nrthindams>0) then   
      NF90(nf90_def_dim(his_file%ncid, 'thindams', nrthindams, his_file%thindams_dimid)) ! nr of structures (thin dam)
   endif   
   !
   if (nr_runup_gauges > 0) then   
      NF90(nf90_def_dim(his_file%ncid, 'runup_gauges', nr_runup_gauges, his_file%runup_gauges_dimid)) ! runup gauges
   endif
   !
   NF90(nf90_def_dim(his_file%ncid, 'pointnamelength', 256, his_file%pointnamelength_dimid)) ! length of station_name per obs point  
   NF90(nf90_def_dim(his_file%ncid, 'runtime', 1, his_file%runtime_dimid)) ! total_runtime, average_dt    
   !
   ! Some metadata attributes
   NF90(nf90_put_att(his_file%ncid, nf90_global, "Conventions",   "CF-1.8"))
   NF90(nf90_put_att(his_file%ncid, nf90_global, "featureType",   "timeSeries"))
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
   NF90(nf90_put_att(his_file%ncid, his_file%station_name_varid, 'cf_role', 'timeseries_id'))
   !
   if (nrcrosssections>0) then
      NF90(nf90_def_var(his_file%ncid, 'crosssection_name', NF90_CHAR, (/his_file%pointnamelength_dimid, his_file%crosssections_dimid/), his_file%crosssection_name_varid))
   endif      
   !
   if (nr_runup_gauges > 0) then
      NF90(nf90_def_var(his_file%ncid, 'runup_gauge_name', NF90_CHAR, (/his_file%pointnamelength_dimid, his_file%runup_gauges_dimid/), his_file%runup_gauge_name_varid))
   endif
   !
   if (nr_src_structures > 0) then
      NF90(nf90_def_var(his_file%ncid, 'drainage_name', NF90_CHAR, (/his_file%pointnamelength_dimid, his_file%drain_dimid/), his_file%drain_name_varid))
   endif
   !
   if (nr_discharge_points > 0 .and. store_river_discharge) then
      NF90(nf90_def_var(his_file%ncid, 'river_name', NF90_CHAR, (/his_file%pointnamelength_dimid, his_file%river_dimid/), his_file%river_name_varid))
   endif
   !
   if (nr_urban_drainage_zones > 0 .and. store_urban_drainage_discharge) then
      NF90(nf90_def_var(his_file%ncid, 'urban_drainage_zone_name', NF90_CHAR, (/his_file%pointnamelength_dimid, his_file%urbdrain_dimid/), his_file%urbdrain_name_varid))
   endif
   !
   !NF90(nf90_put_att(his_file%ncid, his_file%station_name_varid, 'units', '-')) !not wanted in fews
   !
   ! Domain
   ! Station coordinates: input lat/lon or x/y (CRS-aware via crsgeo)
   call def_his_point_coord('station_x', 'x', his_file%station_x_varid, 'original_x_coordinate_of_station')
   call def_his_point_coord('station_y', 'y', his_file%station_y_varid, 'original_y_coordinate_of_station')
   call def_his_point_coord('point_x',   'x', his_file%point_x_varid,   'point_x')
   call def_his_point_coord('point_y',   'y', his_file%point_y_varid,   'point_y')
   !
   NF90(nf90_def_var(his_file%ncid, 'crs', NF90_INT, his_file%crs_varid)) ! For EPSG code
   NF90(nf90_put_att(his_file%ncid, his_file%crs_varid, 'epsg',      epsg))
   NF90(nf90_put_att(his_file%ncid, his_file%crs_varid, 'epsg_code', 'EPSG:' // trim(epsg_code) ))   !--> add epsg_code like FEWS wants
   !
   call ncdef_float_var(his_file%ncid, 'point_zb', (/his_file%points_dimid/), his_file%zb_varid, &
        'm', 'Bed level above reference level', standard_name='altitude', coordinates=pt_coord)
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
      call ncdef_float_var(his_file%ncid, 'structure_height', (/his_file%structures_dimid/), his_file%structure_height_varid, &
           'm', 'interpolated_structure_height_above_reference_level', standard_name='altitude')
      !
   endif
   !
   if (nrthindams>0) then
      !
      NF90(nf90_def_var(his_file%ncid, 'thindam_x', NF90_FLOAT, (/his_file%thindams_dimid/), his_file%thindam_x_varid))  ! snapped coordinate as used in sfincs
      NF90(nf90_put_att(his_file%ncid, his_file%thindam_x_varid, 'units', 'm'))
      NF90(nf90_put_att(his_file%ncid, his_file%thindam_x_varid, 'standard_name', 'projection_x_coordinate'))
      NF90(nf90_put_att(his_file%ncid, his_file%thindam_x_varid, 'long_name', 'thindam_x'))    
      NF90(nf90_put_att(his_file%ncid, his_file%thindam_x_varid, 'grid_mapping', 'crs'))       
      !
      NF90(nf90_def_var(his_file%ncid, 'thindam_y', NF90_FLOAT, (/his_file%thindams_dimid/), his_file%thindam_y_varid))  ! snapped coordinate as used in sfincs
      NF90(nf90_put_att(his_file%ncid, his_file%thindam_y_varid, 'units', 'm'))
      NF90(nf90_put_att(his_file%ncid, his_file%thindam_y_varid, 'standard_name', 'projection_y_coordinate'))
      NF90(nf90_put_att(his_file%ncid, his_file%thindam_y_varid, 'long_name', 'thindam_y'))    
      NF90(nf90_put_att(his_file%ncid, his_file%thindam_y_varid, 'grid_mapping', 'crs'))   
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
   call def_time_point_float('point_zs', his_file%zs_varid, 'm', 'Water level', &
        standard_name='sea_surface_height_above_reference_level')
   !
   if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then
      call def_time_point_float('point_h', his_file%h_varid, 'm', 'Water depth', standard_name='depth')
   endif
   !
   if (store_velocity) then
      call def_time_point_float('point_u', his_file%u_varid, 'm s-1', 'Flow velocity x-component', &
           standard_name='sea_water_x_velocity')
      !
      call def_time_point_float('point_v', his_file%v_varid, 'm s-1', 'Flow velocity y-component', &
           standard_name='sea_water_y_velocity')
      !
      call def_time_point_float('point_uvmag', his_file%uvmag_varid, 'm s-1', 'Flow velocity magnitude', &
           standard_name='sea_water_velocity')
      !
      call def_time_point_float('point_uvdir', his_file%uvdir_varid, 'degrees', 'Flow velocity bearing (deg)', &
           standard_name='sea_water_velocity_direction')
   endif
   !
   ! Add infiltration
   !
   if (infiltration) then
      call def_time_point_float('point_qinf', his_file%qinf_varid, 'mm hr-1', 'Infiltration rate')
   endif
   !
   if (infiltration) then
      if (inftype == 'cnb') then
         call def_time_point_float('point_S', his_file%S_varid, 'm', 'current moisture storage (Se) capacity')
      elseif (inftype == 'gai') then
         call def_time_point_float('point_S', his_file%S_varid, 'm', 'maximum soil moisture deficit')
      endif
   endif
   !
   if (snapwave) then
      call def_time_point_float('point_hm0', his_file%hm0_varid, 'm', 'Hm0 wave height', standard_name='hm0_wave_height')
      !
      call def_time_point_float('point_hm0ig', his_file%hm0ig_varid, 'm', 'Hm0 infragravity wave height', &
           standard_name='hm0_ig_wave_height')
      !
      call def_time_point_float('point_tp', his_file%tp_varid, 's', 'Peak wave period', standard_name='peak_wave_period')
      !
      call def_time_point_float('point_tpig', his_file%tpig_varid, 's', 'Peak wave period Infragravity wave', &
           standard_name='ig_peak_wave_period')
      !
      if (store_wave_direction) then
         call def_time_point_float('point_wavdir', his_file%wavdir_varid, 'degrees', 'Mean wave angle (deg)', &
              standard_name='mean_wave_direction')
         ! point_dirspr is gathered in ncoutput_update_his but the put_var is
         ! currently commented out — do not define here either, otherwise the
         ! file carries an empty variable.
      endif
      !
      if (wavemaker) then
         call def_time_point_float('point_zsm', his_file%zsm_varid, 'm', 'Filtered water level', &
              standard_name='filtered_water_level')
      endif
      !
      if (store_wave_forces) then
         call def_time_point_float('point_dw', his_file%dw_varid, 'm', 'directionally averaged wave breaking dissipation', &
              standard_name='directionally_averaged_wave_breaking_dissipation')
         !
         call def_time_point_float('point_df', his_file%df_varid, 'm', 'directionally averaged wave friction dissipation', &
              standard_name='directionally_averaged_wave_friction_dissipation')
         !
         call def_time_point_float('point_dwig', his_file%dwig_varid, 'm', 'directionally averaged wave breaking dissipation ig', &
              standard_name='directionally_averaged_wave_breaking_dissipation_ig')
         !
         call def_time_point_float('point_dfig', his_file%dfig_varid, 'm', 'directionally averaged wave friction dissipation ig', &
              standard_name='directionally_averaged_wave_friction_dissipation_ig')
         !
         call def_time_point_float('point_cg', his_file%cg_varid, 'm/s', 'wave group velocity', &
              standard_name='wave_group_velocity')
         !
         call def_time_point_float('point_beta', his_file%beta_varid, '-', 'directionally averaged normalised bed slope', &
              standard_name='directionally_averaged_local_bed_slope')
         !
         call def_time_point_float('point_srcig', his_file%srcig_varid, '-', 'directionally averaged ig energy source', &
              standard_name='directionally_averaged_ig_energy_source')
         !
         call def_time_point_float('point_alphaig', his_file%alphaig_varid, '-', &
              'directionally averaged infragravity waves shoaling factor', &
              standard_name='directionally_averaged_infragravity_waves_shoaling_factor')
      endif
   endif
   !
   if (store_meteo) then
      if (wind) then
         call def_time_point_float('point_wind_speed', his_file%wind_speed_varid, 'm s-1', 'Wind speed')
         !
         call def_time_point_float('point_wind_direction', his_file%wind_dir_varid, 'degrees', 'Wind direction (deg)')
      endif
      if (patmos) then
         call def_time_point_float('point_patm', his_file%patm_varid, 'Pa', 'Surface air pressure')
      endif
      if (precip) then
         call def_time_point_float('point_prcp', his_file%prcp_varid, 'mm hr-1', 'Precipitation rate')
         if (store_cumulative_precipitation) then
            call def_time_point_float('point_cumprcp', his_file%cumprcp_varid, 'm', 'Cumulative precipitation')
         endif
      endif
   endif
   !   
   if (nrcrosssections>0) then
      call ncdef_float_var(his_file%ncid, 'crosssection_discharge', (/his_file%crosssections_dimid, his_file%time_dimid/), his_file%discharge_varid, &
           'm3 s-1', 'discharge', coordinates='crosssection_name')
   endif
   !
   if (ndrn>0) then
      call ncdef_float_var(his_file%ncid, 'drainage_discharge', (/his_file%drain_dimid, his_file%time_dimid/), his_file%drain_varid, &
           'm3 s-1', 'discharge through drainage structure')
   endif
   !
   if (nr_runup_gauges > 0) then
      call ncdef_float_var(his_file%ncid, 'runup_gauge_zs', (/his_file%runup_gauges_dimid, his_file%time_dimid/), his_file%runup_gauge_zs_varid, &
           'm', 'run-up elevation', coordinates='runup_gauge_name')
   endif
   !
   call ncdef_float_var(his_file%ncid, 'total_runtime', (/his_file%runtime_dimid/), his_file%total_runtime_varid, &
        's', 'Total model runtime (s)')
   !
   call ncdef_float_var(his_file%ncid, 'average_dt', (/his_file%runtime_dimid/), his_file%average_dt_varid, &
        's', 'Average model time step (s)')
   !
   call ncdef_float_var(his_file%ncid, 'status', (/his_file%runtime_dimid/), his_file%status_varid, &
        '-', 'status of SFINCS simulation - 0 is no error')
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
   if (nr_runup_gauges > 0) then
      NF90(nf90_put_var(his_file%ncid, his_file%runup_gauge_name_varid, runup_gauge_name))  ! write rug name
   endif
   !
   if (nr_src_structures > 0) then
      !
      ! Copy src_struc_name (length src_struc_name_len = 128) into a length-256 buffer
      ! to match the pointnamelength netCDF dimension used for all his_file name
      ! variables.
      !
      allocate(drain_name_buf(nr_src_structures))
      !
      do istruc = 1, nr_src_structures
         !
         drain_name_buf(istruc) = src_struc_name(istruc)
         !
      enddo
      !
      NF90(nf90_put_var(his_file%ncid, his_file%drain_name_varid, drain_name_buf))  ! write drainage_name
      !
      deallocate(drain_name_buf)
      !
   endif
   !
   if (nr_discharge_points > 0 .and. store_river_discharge) then
      !
      ! Copy src_name (length src_name_len) into a length-256 buffer to match
      ! the pointnamelength netCDF dimension.
      !
      allocate(river_name_buf(nr_discharge_points))
      !
      do istruc = 1, nr_discharge_points
         !
         river_name_buf(istruc) = src_name(istruc)
         !
      enddo
      !
      NF90(nf90_put_var(his_file%ncid, his_file%river_name_varid, river_name_buf))  ! write river_name
      !
      deallocate(river_name_buf)
      !
   endif
   !
   if (nr_urban_drainage_zones > 0 .and. store_urban_drainage_discharge) then
      !
      allocate(urbdrain_name_buf(nr_urban_drainage_zones))
      !
      do istruc = 1, nr_urban_drainage_zones
         !
         urbdrain_name_buf(istruc) = urb_zone_name(istruc)
         !
      enddo
      !
      NF90(nf90_put_var(his_file%ncid, his_file%urbdrain_name_varid, urbdrain_name_buf))  ! write urban_drainage_zone_name
      !
      deallocate(urbdrain_name_buf)
      !
   endif
   !
   if (nrstructures>0) then
      !
      ! Allocate structure
      call give_structure_information(struc_info)
      !
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
   if (nrthindams>0) then
      !
      ! Allocate structure
      call give_thindam_information(thindam_info)
      !
      allocate(thindam_x(nrthindams))
      allocate(thindam_y(nrthindams))
      !
      do istruc = 1, nrthindams
         !
         thindam_x(istruc)        = thindam_info(istruc,1)
         thindam_y(istruc)        = thindam_info(istruc,2)
         !
      enddo      
      NF90(nf90_put_var(his_file%ncid,  his_file%thindam_x_varid, thindam_x)) ! write thindam_x
      !
      NF90(nf90_put_var(his_file%ncid,  his_file%thindam_y_varid, thindam_y)) ! write thindam_y
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
   end subroutine ncoutput_his_init
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   subroutine ncoutput_update_map(t,ntmapout)
   !
   ! Write time-varying output to map file. Single linear flow — the
   ! grid-type branch is hidden inside write_cell_var / write_cell_var_wet
   ! and the two precompute helpers (compute_uv_at_cell_centers,
   ! compute_pnh_unwrapped) below.
   !
   use sfincs_data
   use sfincs_nonhydrostatic
   use sfincs_snapwave
   use quadtree
   !
   implicit none
   !
   real*8  :: t
   integer :: ntmapout
   real*4  :: sq2
   real*4, dimension(:), allocatable :: uxy, vxy, pnh_full
   !
   sq2 = sqrt(2.0)
   !
   NF90(nf90_put_var(map_file%ncid, map_file%time_varid, t, (/ntmapout/)))
   !
   ! -------------------------------------------------------
   ! Water level / depth
   ! -------------------------------------------------------
   if (subgrid) then
      call write_cell_var_wet(map_file%ncid, map_file%zs_varid, real(zs,4), subgrid_z_zmin, ntmapout)
   else
      call write_cell_var_wet(map_file%ncid, map_file%zs_varid, real(zs,4), zb,             ntmapout)
   endif
   !
   ! Optional time-varying zb (regular grid, non-subgrid)
   if (.not. use_quadtree .and. store_dynamic_bed_level .and. .not. subgrid) then
      call write_cell_var(map_file%ncid, map_file%zb_varid, zb, ntmapout)
   endif
   !
   ! h = zs - zref. Quadtree filters wet cells (legacy); regular keeps all.
   if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then
      if (subgrid) then
         call write_cell_var_depth(map_file%ncid, map_file%h_varid, real(zs,4), subgrid_z_zmin, ntmapout, &
              check_wet=use_quadtree)
      else
         call write_cell_var_depth(map_file%ncid, map_file%h_varid, real(zs,4), zb,             ntmapout, &
              check_wet=use_quadtree)
      endif
   endif
   !
   ! -------------------------------------------------------
   ! Velocity (cell-centered, computed once via 4-face average)
   ! -------------------------------------------------------
   if (store_velocity) then
      allocate(uxy(np), vxy(np))
      call compute_uv_at_cell_centers(uxy, vxy)
      call write_cell_var(map_file%ncid, map_file%u_varid, uxy, ntmapout)
      call write_cell_var(map_file%ncid, map_file%v_varid, vxy, ntmapout)
      deallocate(uxy, vxy)
   endif
   !
   ! -------------------------------------------------------
   ! Subgrid volumes
   ! -------------------------------------------------------
   if (subgrid) then
      if (store_zvolume) then
         call write_cell_var(map_file%ncid, map_file%zvolume_varid, real(z_volume,4), ntmapout)
      endif
      if (store_storagevolume) then
         call write_cell_var(map_file%ncid, map_file%storagevolume_varid, storage_volume, ntmapout)
      endif
   endif
   !
   ! -------------------------------------------------------
   ! Infiltration state
   ! -------------------------------------------------------
   if (inftype == 'cnb') then
      call write_cell_var(map_file%ncid, map_file%infstate_varid, scs_Se,   ntmapout)
   elseif (inftype == 'gai') then
      call write_cell_var(map_file%ncid, map_file%infstate_varid, GA_sigma, ntmapout)
   elseif (inftype == 'hor') then
      call write_cell_var(map_file%ncid, map_file%infstate_varid, qinfmap,  ntmapout)
   endif
   !
   ! -------------------------------------------------------
   ! Meteo
   ! -------------------------------------------------------
   if (store_meteo) then
      if (wind) then
         call write_cell_var(map_file%ncid, map_file%wind_u_varid, windu, ntmapout)
         call write_cell_var(map_file%ncid, map_file%wind_v_varid, windv, ntmapout)
      endif
      if (patmos) then
         call write_cell_var(map_file%ncid, map_file%patm_varid, patm, ntmapout)
      endif
      if (precip) then
         call write_cell_var(map_file%ncid, map_file%precip_varid, prcp, ntmapout, scale=3600000.0)
      endif
   endif
   !
   ! -------------------------------------------------------
   ! SnapWave (all fields read the snapwave_* node arrays via use_sw_index,
   ! so quadtree and regular grids share the same output path)
   ! -------------------------------------------------------
   if (snapwave) then
      call write_cell_var(map_file%ncid, map_file%hm0_varid,   snapwave_H,    ntmapout, use_sw_index=.true., scale=sq2)
      call write_cell_var(map_file%ncid, map_file%hm0ig_varid, snapwave_H_ig, ntmapout, use_sw_index=.true., scale=sq2)
      call write_cell_var(map_file%ncid, map_file%tp_varid,    snapwave_Tp,    ntmapout, use_sw_index=.true.)
      call write_cell_var(map_file%ncid, map_file%tpig_varid,  snapwave_Tp_ig, ntmapout, use_sw_index=.true.)
      if (store_wave_forces) then
         call write_cell_var(map_file%ncid, map_file%fwx_varid,           snapwave_Fx,    ntmapout, use_sw_index=.true.)
         call write_cell_var(map_file%ncid, map_file%fwy_varid,           snapwave_Fy,    ntmapout, use_sw_index=.true.)
         call write_cell_var(map_file%ncid, map_file%beta_varid,          snapwave_beta,  ntmapout, use_sw_index=.true.)
         call write_cell_var(map_file%ncid, map_file%snapwavedepth_varid, snapwave_depth, ntmapout, use_sw_index=.true.)
      endif
      if (store_wave_direction) then
         call write_cell_var(map_file%ncid, map_file%wavdir_varid, snapwave_mean_direction, ntmapout, use_sw_index=.true.)
      endif
      if (wavemaker) then
         call write_cell_var(map_file%ncid, map_file%zsm_varid, zsm, ntmapout)
      endif
   endif
   !
   ! -------------------------------------------------------
   ! Non-hydrostatic pressure (precompute then write)
   ! -------------------------------------------------------
   if (nonhydrostatic) then
      allocate(pnh_full(np))
      call compute_pnh_unwrapped(pnh_full)
      call write_cell_var(map_file%ncid, map_file%pnonh_varid, pnh_full, ntmapout)
      deallocate(pnh_full)
   endif
   !
   NF90(nf90_sync(map_file%ncid))
   !
   end subroutine ncoutput_update_map

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   subroutine ncoutput_update_his(t,nthisout)
   !
   ! Write time, zs, u, v, prcp etc. at observation points.
   !
   use sfincs_data
   use sfincs_crosssections
   use sfincs_runup_gauges
   use sfincs_snapwave
   use sfincs_src_structures, only: nr_src_structures, src_struc_q_now, src_struc_breach_width, src_struc_type, structure_dike_breach, src_struc_fraction_open
   use sfincs_discharges,     only: qtsrc, nr_discharge_points
   use sfincs_urban_drainage, only: nr_urban_drainage_zones, urban_drainage_q_total
   !
   implicit none
   !
   integer :: iobs, nm, idrn
   integer :: nthisout
   real*8  :: t
   real*4, dimension(nobs) :: uobs, vobs, uvmag, uvdir
   real*4, dimension(nobs) :: twndmag, twnddir
   real*4, dimension(ndrn) :: q_drain
   real*4, dimension(:), allocatable :: qq, zz
   !
   zobs    = FILL_VALUE
   hobs    = FILL_VALUE
   q_drain = FILL_VALUE
   !
   do iobs = 1, nobs
      nm = nmindobs(iobs)
      if (nm>0) then
         zobs(iobs) = zs(nm)
         if (subgrid) then
            hobs(iobs) = zs(nm) - subgrid_z_zmin(nm)
         else
            hobs(iobs) = zs(nm) - zb(nm)
         endif
      endif
   enddo
   if (store_velocity)          call compute_uv_at_obs_points(uobs, vobs, uvmag, uvdir)
   if (store_meteo .and. wind)  call compute_wind_at_obs_points(twndmag, twnddir)
   !
   NF90(nf90_put_var(his_file%ncid, his_file%time_varid, t, (/nthisout/)))
   !
   NF90(nf90_put_var(his_file%ncid, his_file%zs_varid, zobs, (/1, nthisout/)))
   !
   if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then
      !
      NF90(nf90_put_var(his_file%ncid, his_file%h_varid, hobs, (/1, nthisout/)))
      !
   endif
   !
   if (infiltration) then
      !
      call write_point_var(his_file%qinf_varid, qinfmap, nthisout, scale=3600000.0)
      !
      if (inftype == 'cnb') then
         call write_point_var(his_file%S_varid, scs_Se, nthisout)
      elseif (inftype == 'gai') then
         call write_point_var(his_file%S_varid, GA_sigma, nthisout)
      endif
      !
   endif
   !
   if (snapwave) then
      !
      call write_point_var(his_file%hm0_varid,    snapwave_H,     nthisout, use_sw_index=.true., scale=sqrt(2.0))
      call write_point_var(his_file%hm0ig_varid,  snapwave_H_ig,  nthisout, use_sw_index=.true., scale=sqrt(2.0))
      call write_point_var(his_file%tp_varid,     snapwave_Tp,    nthisout, use_sw_index=.true.)
      call write_point_var(his_file%tpig_varid,   snapwave_Tp_ig, nthisout, use_sw_index=.true.)
      !
      if (store_wave_direction) then
         call write_point_var(his_file%wavdir_varid, snapwave_mean_direction, nthisout, use_sw_index=.true.)
      endif
      !
      if (wavemaker) then
         !
         call write_point_var(his_file%zsm_varid, zsm, nthisout)
         !
      endif
      !
      if (store_wave_forces) then
         !
         call write_point_var(his_file%dw_varid,      snapwave_Dw,      nthisout, use_sw_index=.true.)
         call write_point_var(his_file%df_varid,      snapwave_Df,      nthisout, use_sw_index=.true.)
         call write_point_var(his_file%dwig_varid,    snapwave_Dwig,    nthisout, use_sw_index=.true.)
         call write_point_var(his_file%dfig_varid,    snapwave_Dfig,    nthisout, use_sw_index=.true.)
         call write_point_var(his_file%cg_varid,      snapwave_cg,      nthisout, use_sw_index=.true.)
         call write_point_var(his_file%beta_varid,    snapwave_beta,    nthisout, use_sw_index=.true.)
         call write_point_var(his_file%srcig_varid,   snapwave_srcig,   nthisout, use_sw_index=.true.)
         call write_point_var(his_file%alphaig_varid, snapwave_alphaig, nthisout, use_sw_index=.true.)
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
         call write_point_var(his_file%patm_varid, patm, nthisout)
         !
      endif
      !
      if (precip) then
         !
         call write_point_var(his_file%prcp_varid, prcp, nthisout, scale=3600000.0)
         !
         if (store_cumulative_precipitation) then
            !
            call write_point_var(his_file%cumprcp_varid, cumprcp, nthisout)
            !
         endif
         !
      endif
      !
   endif
   !
   if (nrcrosssections>0) then
      !$acc update host(q)
      ! Get fluxes through cross sections (callee allocates qq)
      call get_discharges_through_crosssections(qq)
      NF90(nf90_put_var(his_file%ncid, his_file%discharge_varid, qq, (/1, nthisout/)))
      if (allocated(qq)) deallocate(qq)
   endif
   !
   if (nr_runup_gauges>0) then
      ! Get run-up elevations (callee allocates zz)
      call get_runup_levels(zz)
      NF90(nf90_put_var(his_file%ncid, his_file%runup_gauge_zs_varid, zz, (/1, nthisout/)))
      if (allocated(zz)) deallocate(zz)
   endif
   !
   if (nr_src_structures>0) then
      !
      !$acc update host(src_struc_q_now, src_struc_fraction_open)
      !
      NF90(nf90_put_var(his_file%ncid, his_file%drain_varid, src_struc_q_now, (/1, nthisout/))) ! write per-structure discharge
      NF90(nf90_put_var(his_file%ncid, his_file%drain_fraction_open_varid, src_struc_fraction_open, (/1, nthisout/))) ! write per-structure gate open fraction
      !
      if (any(src_struc_type == structure_dike_breach)) then
         !$acc update host(src_struc_breach_width)
         NF90(nf90_put_var(his_file%ncid, his_file%breach_width_varid, src_struc_breach_width, (/1, nthisout/))) ! write breach width
      endif
      !
   endif
   !
   if (nr_discharge_points>0 .and. store_river_discharge) then
      !
      !$acc update host(qtsrc)
      ! Get fluxes through drainage structure
      !
      idrn = 0
      do iobs = nsrc + 1, nsrcdrn, 2 !TL: as in sfincs_output.f90
         idrn = idrn + 1
         q_drain(idrn) = qtsrc(iobs)
      enddo
      !
      NF90(nf90_put_var(his_file%ncid, his_file%drain_varid, q_drain, (/1, nthisout/)))
      !
   endif
   !
   if (store_velocity) then
      !
      NF90(nf90_put_var(his_file%ncid, his_file%u_varid,     uobs,  (/1, nthisout/)))
      NF90(nf90_put_var(his_file%ncid, his_file%v_varid,     vobs,  (/1, nthisout/)))
      NF90(nf90_put_var(his_file%ncid, his_file%uvmag_varid, uvmag, (/1, nthisout/)))
      NF90(nf90_put_var(his_file%ncid, his_file%uvdir_varid, uvdir, (/1, nthisout/)))
      !
   endif
   !
   NF90(nf90_sync(his_file%ncid)) !TL: in first test it seems to be faster to let the file update than keep in memory
   !
   end subroutine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   !

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   subroutine ncoutput_update_max(t,ntmaxout)
   !
   ! Write maximum values to map file (handles both regular and quadtree grids)
   !
   use sfincs_data
   use quadtree
   use sfincs_urban_drainage, only: urban_drainage_cumulative_volume
   !
   implicit none
   !
   real*8                            :: t
   integer                           :: ntmaxout, nm
   real*4, dimension(:), allocatable :: hmax_out, hmean
   !
   ! Scalar time of this max-record (defined only when store_maximum_waterlevel)
   if (store_maximum_waterlevel) then
      NF90(nf90_put_var(map_file%ncid, map_file%timemax_varid, t, (/ntmaxout/)))
   endif
   !
   ! Maximum water level
   if (store_maximum_waterlevel) then
      if (subgrid) then
         call write_cell_var_wet(map_file%ncid, map_file%zsmax_varid, zsmax, subgrid_z_zmin, ntmaxout)
      else
         call write_cell_var_wet(map_file%ncid, map_file%zsmax_varid, zsmax, zb,             ntmaxout)
      endif
   endif
   !
   ! Maximum water depth (optional, supports subgrid mean-depth)
   if (store_maximum_waterlevel .and. (.not. subgrid .or. store_hsubgrid)) then
      if (subgrid .and. store_hmean) then
         ! 
         ! Subgrid mean depth needs a per-cell precompute; mask to wet cells and write.
         allocate(hmax_out(np))
         allocate(hmean(np))
         hmax_out = FILL_VALUE
         !
         call compute_subgrid_mean_depth(zsmax, hmean)
         !
         do nm = 1, np
            if ( (zsmax(nm) - subgrid_z_zmin(nm)) > huthresh) then
               hmax_out(nm) = hmean(nm)
            endif
         enddo
         !
         call write_cell_var(map_file%ncid, map_file%hmax_varid, hmax_out, ntmaxout, check_kcs=.true.)
         !
         deallocate(hmax_out)
         deallocate(hmean)
         !
      elseif (subgrid) then
         call write_cell_var_depth(map_file%ncid, map_file%hmax_varid, zsmax, subgrid_z_zmin, ntmaxout)
      else
         call write_cell_var_depth(map_file%ncid, map_file%hmax_varid, zsmax, zb,             ntmaxout)
      endif
   endif
   !
   ! Cumulative rainfall (always when store_cumulative_precipitation) and
   ! cumulative infiltration (only when infiltration is on — same gate as the def)
   if (store_cumulative_precipitation) then
      call write_cell_var(map_file%ncid, map_file%cumprcp_varid, cumprcp, ntmaxout, check_kcs=.true.)
      if (infiltration) then
         call write_cell_var(map_file%ncid, map_file%cuminf_varid, cuminf, ntmaxout, check_kcs=.true.)
      endif
   endif
   !
   ! Maximum flow velocity / flux
   if (store_maximum_velocity) then
      call write_cell_var(map_file%ncid, map_file%vmax_varid, vmax, ntmaxout, check_kcs=.true.)
   endif
   if (store_maximum_flux) then
      call write_cell_var(map_file%ncid, map_file%qmax_varid, qmax, ntmaxout, check_kcs=.true.)
   endif
   !
   ! Duration wet cell
   if (store_twet) then
      call write_cell_var(map_file%ncid, map_file%tmax_varid, twet, ntmaxout, check_kcs=.true.)
   endif
   !
   ! When zsmax occurred (only cells where it actually happened, i.e. t_zsmax > 0)
   if (store_t_zsmax) then
      call write_cell_var(map_file%ncid, map_file%t_zsmax_varid, t_zsmax, ntmaxout, &
           check_kcs=.true., min_value=0.0)
   endif
   !
   ! Maximum wind speed (def is gated on store_meteo .and. wind .and.
   ! store_wind_max .and. meteo3d in ncoutput_map_init; mirror that here)
   if (store_meteo .and. wind .and. store_wind_max .and. meteo3d) then
      call write_cell_var(map_file%ncid, map_file%windmax_varid, windmax, ntmaxout, check_kcs=.true.)
   endif
   !
   end subroutine ncoutput_update_max

   subroutine ncoutput_map_finalize() 
   !
   ! Add total runtime, dtavg to file and close
   !
   use sfincs_data
   !
   implicit none
   !
   if (store_tsunami_arrival_time) then
      !
      call ncoutput_write_tsunami_arrival_time()
      !
   endif
   !
   if (timestep_analysis) then
       !
       call ncoutput_write_timestep_analysis()
       !
   endif
   !
   NF90(nf90_put_var(map_file%ncid, map_file%total_runtime_varid, tfinish_all - tstart_all))
   NF90(nf90_put_var(map_file%ncid, map_file%average_dt_varid,  dtavg))
   NF90(nf90_put_var(map_file%ncid, map_file%status_varid,  error))
   !
   NF90(nf90_close(map_file%ncid))
   !
   end subroutine ncoutput_map_finalize
   !
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
   ! Mirror the early-return condition from ncoutput_his_init exactly: if
   ! none of these are present, no his file was created. (Note: thindams
   ! alone do NOT trigger his-file creation in init, so they're not in
   ! this list either.)
   if (nobs==0 .and. nrcrosssections==0 .and. nrstructures==0 .and. ndrn==0 .and. nr_runup_gauges==0) then
      return
   endif
   !
   implicit none
   !
   if (nobs==0 .and. nrcrosssections==0 .and. nrstructures==0 .and. nrthindams==0 .and. nr_src_structures==0 .and. .not. (nr_discharge_points>0 .and. store_river_discharge) .and. .not. (nr_urban_drainage_zones>0 .and. store_urban_drainage_discharge)) then ! If no observation points, cross-sections, structures (weir or thin dam), drains, river sources or urban drainage zones; hisfile
        return
   endif
   !
   NF90(nf90_put_var(his_file%ncid, his_file%total_runtime_varid, tfinish_all - tstart_all))
   NF90(nf90_put_var(his_file%ncid, his_file%average_dt_varid,  dtavg)) 
   NF90(nf90_put_var(his_file%ncid, his_file%status_varid,  error))       
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
   use sfincs_src_structures, only: drnfile
   use sfincs_discharges,     only: srcfile, disfile, netsrcdisfile
   !
   ! Because of overlapping names, only important specific values from snapwave_data
   use snapwave_data, only: gamma, gammax, alpha, hmin, fw0, fw0_ig, dt, tol, dtheta, crit, nr_sweeps, baldock_exponent, baldock_ratio, &
       igwaves_opt, alpha_ig, gamma_ig, gamma_fac_br, shinc2ig, alphaigfac, baldock_ratio_ig, ig_opt, herbers_opt, tpig_opt, eeinc2ig, tinc2ig, &
       snapwave_jonswapfile, snapwave_encfile, snapwave_bndfile, snapwave_bhsfile, snapwave_btpfile, snapwave_bwdfile, snapwave_bdsfile, upwfile, gridfile, &
       jonswapgam, Tpini, sector, fwratio, fwigratio   
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
        NF90(nf90_put_att(ncid, varid, 'dtmax',dtmax))        
        NF90(nf90_put_att(ncid, varid, 'dtmin',dtmin))  
        NF90(nf90_put_att(ncid, varid, 'hmin_cfl',hmin_cfl))        
        NF90(nf90_put_att(ncid, varid, 'huthresh',huthresh))
        NF90(nf90_put_att(ncid, varid, 'huvmin',huvmin))
        NF90(nf90_put_att(ncid, varid, 'rhoa',rhoa))        
        NF90(nf90_put_att(ncid, varid, 'rhow',rhow))        
        NF90(nf90_put_att(ncid, varid, 'inputformat',inputtype))        
        NF90(nf90_put_att(ncid, varid, 'outputformat',outputtype))   
        NF90(nf90_put_att(ncid, varid, 'outputtype_map',outputtype_map))
        NF90(nf90_put_att(ncid, varid, 'outputtype_his',outputtype_his))
        NF90(nf90_put_att(ncid, varid, 'bndtype',bndtype))
        NF90(nf90_put_att(ncid, varid, 'advection',logical2int(advection)))  
        NF90(nf90_put_att(ncid, varid, 'latitude',latitude))  
        NF90(nf90_put_att(ncid, varid, 'pavbnd',pavbnd))  
        NF90(nf90_put_att(ncid, varid, 'gapres',gapres))  
        NF90(nf90_put_att(ncid, varid, 'baro',baro))  
        NF90(nf90_put_att(ncid, varid, 'utmzone',utmzone))  
        NF90(nf90_put_att(ncid, varid, 'epsg',epsg))  
        NF90(nf90_put_att(ncid, varid, 'epsg_code',epsg_code))  
        NF90(nf90_put_att(ncid, varid, 'advlim',advlim))  
        NF90(nf90_put_att(ncid, varid, 'uvlim',uvlim))  
        NF90(nf90_put_att(ncid, varid, 'uvmax',uvmax))  
        NF90(nf90_put_att(ncid, varid, 'advection_mask',logical2int(advection_mask)))
        NF90(nf90_put_att(ncid, varid, 'friction2d',logical2int(friction2d)))  
        NF90(nf90_put_att(ncid, varid, 'wiggle_suppression',logical2int(wiggle_suppression)))  
        NF90(nf90_put_att(ncid, varid, 'wiggle_factor',wiggle_factor))  
        NF90(nf90_put_att(ncid, varid, 'wiggle_threshold',wiggle_threshold)) 
        NF90(nf90_put_att(ncid, varid, 'slopelim',slopelim))            
        NF90(nf90_put_att(ncid, varid, 'qinf_zmin',qinf_zmin))    
        NF90(nf90_put_att(ncid, varid, 'btfilter', btfilter))                     
        NF90(nf90_put_att(ncid, varid, 'sfacinf', sfacinf))  
        NF90(nf90_put_att(ncid, varid, 'radstr', logical2int(radstr))) 
        NF90(nf90_put_att(ncid, varid, 'crsgeo',logical2int(crsgeo))) 
        NF90(nf90_put_att(ncid, varid, 'coriolis',logical2int(coriolis))) 
        NF90(nf90_put_att(ncid, varid, 'amprblock',logical2int(ampr_block))) 
        NF90(nf90_put_att(ncid, varid, 'spwmergefrac',spw_merge_frac)) 
        NF90(nf90_put_att(ncid, varid, 'usespwprecip',logical2int(use_spw_precip)))         
        NF90(nf90_put_att(ncid, varid, 'global',logical2int(global))) 
        NF90(nf90_put_att(ncid, varid, 'spinup_meteo', logical2int(spinup_meteo))) 
        NF90(nf90_put_att(ncid, varid, 'waveage',waveage)) 
        NF90(nf90_put_att(ncid, varid, 'wavemaker_filter_time', wavemaker_filter_time))         
        NF90(nf90_put_att(ncid, varid, 'wavemaker_filter_fred',wavemaker_filter_fred))          
        NF90(nf90_put_att(ncid, varid, 'wavemaker_hmin',wavemaker_hmin))  
        NF90(nf90_put_att(ncid, varid, 'wavemaker_nfreqs_inc',wavemaker_nfreqs_inc))  
        NF90(nf90_put_att(ncid, varid, 'wavemaker_freqmin_inc',wavemaker_freqmin_inc))  
        NF90(nf90_put_att(ncid, varid, 'wavemaker_freqmax_inc',wavemaker_freqmax_inc))  
        NF90(nf90_put_att(ncid, varid, 'wavemaker_nfreqs_ig',wavemaker_nfreqs_ig))  
        NF90(nf90_put_att(ncid, varid, 'wavemaker_freqmin_ig',wavemaker_freqmin_ig))  
        NF90(nf90_put_att(ncid, varid, 'wavemaker_freqmax_ig',wavemaker_freqmax_ig))  
        NF90(nf90_put_att(ncid, varid, 'wavemaker_tinc2ig',wavemaker_tinc2ig))         
        NF90(nf90_put_att(ncid, varid, 'wavemaker_surfslope',wavemaker_surfslope))         
        NF90(nf90_put_att(ncid, varid, 'wavemaker_hm0_ig_factor',wavemaker_hm0_ig_factor))         
        NF90(nf90_put_att(ncid, varid, 'wavemaker_hm0_inc_factor',wavemaker_hm0_inc_factor))         
        NF90(nf90_put_att(ncid, varid, 'wavemaker_gammax',wavemaker_gammax))         
        NF90(nf90_put_att(ncid, varid, 'wavemaker_tpmin',wavemaker_tpmin))
        NF90(nf90_put_att(ncid, varid, 'horton_kr_kd',horton_kr_kd))
        NF90(nf90_put_att(ncid, varid, 'nuviscdim',nuviscdim))
        NF90(nf90_put_att(ncid, varid, 'nuviscfac',nuviscfac))
        NF90(nf90_put_att(ncid, varid, 'btrelax',btrelax))
        NF90(nf90_put_att(ncid, varid, 'structure_relax',structure_relax))
        NF90(nf90_put_att(ncid, varid, 'wave_enhanced_roughness',logical2int(wave_enhanced_roughness)))
        NF90(nf90_put_att(ncid, varid, 'bathtub',logical2int(bathtub)))
        NF90(nf90_put_att(ncid, varid, 'bathtub_fac_hs',bathtub_fac_hs))
        NF90(nf90_put_att(ncid, varid, 'bathtub_dt',bathtub_dt))
        NF90(nf90_put_att(ncid, varid, 'factor_wind',factor_wind))
        NF90(nf90_put_att(ncid, varid, 'factor_pres',factor_pres))
        NF90(nf90_put_att(ncid, varid, 'factor_prcp',factor_prcp))
        NF90(nf90_put_att(ncid, varid, 'factor_spw_size',factor_spw_size))
        NF90(nf90_put_att(ncid, varid, 'nonh',logical2int(nonhydrostatic)))
        NF90(nf90_put_att(ncid, varid, 'nh_fnudge',nh_fnudge))
        NF90(nf90_put_att(ncid, varid, 'nh_tstop',nh_tstop))
        NF90(nf90_put_att(ncid, varid, 'nh_tol',nh_tol))
        NF90(nf90_put_att(ncid, varid, 'nh_itermax',nh_itermax))
        !
        ! Domain
        !
        NF90(nf90_put_att(ncid, varid, 'qtrfile',qtrfile))        
        NF90(nf90_put_att(ncid, varid, 'depfile',depfile))        
        NF90(nf90_put_att(ncid, varid, 'inifile',zsinifile))        
        NF90(nf90_put_att(ncid, varid, 'rstfile',rstfile))        
        NF90(nf90_put_att(ncid, varid, 'mskfile',mskfile))        
        NF90(nf90_put_att(ncid, varid, 'indexfile',indexfile))        
        NF90(nf90_put_att(ncid, varid, 'sbgfile',sbgfile))        
        NF90(nf90_put_att(ncid, varid, 'thdfile',thdfile))        
        NF90(nf90_put_att(ncid, varid, 'weirfile',weirfile))        
        NF90(nf90_put_att(ncid, varid, 'manningfile',manningfile))    
        NF90(nf90_put_att(ncid, varid, 'drnfile',drnfile))
        NF90(nf90_put_att(ncid, varid, 'rugfile',rugfile))
        NF90(nf90_put_att(ncid, varid, 'cstfile',cstfile))
        NF90(nf90_put_att(ncid, varid, 'volfile',volfile))
        NF90(nf90_put_att(ncid, varid, 'bcafile',bcafile))
        !
        ! Forcing
        !
        NF90(nf90_put_att(ncid, varid, 'bndfile',bndfile))        
        NF90(nf90_put_att(ncid, varid, 'bzsfile',bzsfile))        
        NF90(nf90_put_att(ncid, varid, 'bzifile',bzifile))        
        NF90(nf90_put_att(ncid, varid, 'bdrfile',bdrfile))        
        NF90(nf90_put_att(ncid, varid, 'wavemaker_wfpfile',wavemaker_wfpfile))   
        NF90(nf90_put_att(ncid, varid, 'wavemaker_whifile',wavemaker_whifile))   
        NF90(nf90_put_att(ncid, varid, 'wavemaker_wtifile',wavemaker_wtifile))   
        NF90(nf90_put_att(ncid, varid, 'wavemaker_wstfile',wavemaker_wstfile))   
        NF90(nf90_put_att(ncid, varid, 'srcfile',srcfile))        
        NF90(nf90_put_att(ncid, varid, 'disfile',disfile))        
        NF90(nf90_put_att(ncid, varid, 'spwfile',spwfile))        
        NF90(nf90_put_att(ncid, varid, 'wndfile',wndfile))        
        NF90(nf90_put_att(ncid, varid, 'precipfile',prcpfile))        
        NF90(nf90_put_att(ncid, varid, 'amufile',amufile))        
        NF90(nf90_put_att(ncid, varid, 'amvfile',amvfile))  
        NF90(nf90_put_att(ncid, varid, 'ampfile',ampfile))              
        NF90(nf90_put_att(ncid, varid, 'amprfile',amprfile))  
        NF90(nf90_put_att(ncid, varid, 'infiltrationfile',infiltrationfile))
        NF90(nf90_put_att(ncid, varid, 'infiltrationtype',inftype))
        NF90(nf90_put_att(ncid, varid, 'qinffile',qinffile))
        NF90(nf90_put_att(ncid, varid, 'scsfile',scsfile)) 
        NF90(nf90_put_att(ncid, varid, 'smaxfile',smaxfile)) 
        NF90(nf90_put_att(ncid, varid, 'sefffile',sefffile)) 
        NF90(nf90_put_att(ncid, varid, 'ksfile',ksfile)) 
        NF90(nf90_put_att(ncid, varid, 'psifile',psifile)) 
        NF90(nf90_put_att(ncid, varid, 'sigmafile',sigmafile))
        NF90(nf90_put_att(ncid, varid, 'f0file',f0file))
        NF90(nf90_put_att(ncid, varid, 'fcfile',fcfile))
        NF90(nf90_put_att(ncid, varid, 'kdfile',kdfile))
        NF90(nf90_put_att(ncid, varid, 'z0lfile',z0lfile))
        NF90(nf90_put_att(ncid, varid, 'wavemaker_wvmfile',wavemaker_wvmfile))
        !
        ! Infiltration configuration
        !
        NF90(nf90_put_att(ncid, varid, 'infiltration_file',infiltrationfile))
        NF90(nf90_put_att(ncid, varid, 'infiltration_type',inftype))
        !
        ! Netcdf input
        NF90(nf90_put_att(ncid, varid, 'netbndbzsbzifile',netbndbzsbzifile))
        NF90(nf90_put_att(ncid, varid, 'netsrcdisfile',netsrcdisfile))     
        NF90(nf90_put_att(ncid, varid, 'netamuamvfile',netamuamvfile))     
        NF90(nf90_put_att(ncid, varid, 'netamprfile',netamprfile))  
        NF90(nf90_put_att(ncid, varid, 'netampfile',netampfile))        
        NF90(nf90_put_att(ncid, varid, 'netspwfile',netspwfile))                    
        !
        ! Output
        !
        NF90(nf90_put_att(ncid, varid, 'obsfile',obsfile))   
        NF90(nf90_put_att(ncid, varid, 'nobs',nobs))                    
        NF90(nf90_put_att(ncid, varid, 'crsfile',crsfile))   
        !
        NF90(nf90_put_att(ncid, varid, 'storevelmax',logical2int(store_maximum_velocity)))
        NF90(nf90_put_att(ncid, varid, 'storefluxmax',logical2int(store_maximum_flux)))
        NF90(nf90_put_att(ncid, varid, 'storevel',logical2int(store_velocity)))
        NF90(nf90_put_att(ncid, varid, 'storecumprcp',logical2int(store_cumulative_precipitation)))
        NF90(nf90_put_att(ncid, varid, 'storetwet',logical2int(store_twet)))
        NF90(nf90_put_att(ncid, varid, 'storehsubgrid',logical2int(store_hsubgrid)))
        NF90(nf90_put_att(ncid, varid, 'twet_threshold',twet_threshold))
        NF90(nf90_put_att(ncid, varid, 'store_tsunami_arrival_time',logical2int(store_tsunami_arrival_time)))
        NF90(nf90_put_att(ncid, varid, 'tsunami_arrival_threshold',tsunami_arrival_threshold))
        NF90(nf90_put_att(ncid, varid, 'storeqdrain',logical2int(store_qdrain)))
        NF90(nf90_put_att(ncid, varid, 'storezvolume',logical2int(store_zvolume)))
        NF90(nf90_put_att(ncid, varid, 'writeruntime',logical2int(write_time_output)))
        NF90(nf90_put_att(ncid, varid, 'debug',logical2int(debug)))
        NF90(nf90_put_att(ncid, varid, 'storemeteo',logical2int(store_meteo)))
        NF90(nf90_put_att(ncid, varid, 'storemaxwind',logical2int(store_wind_max))) 
        NF90(nf90_put_att(ncid, varid, 'storefw',logical2int(store_wave_forces)))         
        NF90(nf90_put_att(ncid, varid, 'storewavdir', logical2int(store_wave_direction)))
        NF90(nf90_put_att(ncid, varid, 'storetzsmax',storetzsmax))
        NF90(nf90_put_att(ncid, varid, 'storestoragevolume',storestoragevolume))
        NF90(nf90_put_att(ncid, varid, 'storehmean',logical2int(store_hmean)))
        NF90(nf90_put_att(ncid, varid, 'timestep_analysis',logical2int(timestep_analysis)))
        NF90(nf90_put_att(ncid, varid, 'store_dynamic_bed_level',logical2int(store_dynamic_bed_level)))
        NF90(nf90_put_att(ncid, varid, 'regular_output_on_mesh',logical2int(use_quadtree_output)))
        NF90(nf90_put_att(ncid, varid, 'rugdepth',runup_gauge_depth))
        NF90(nf90_put_att(ncid, varid, 'percentage_done',percdoneval))
        NF90(nf90_put_att(ncid, varid, 'nc_deflate_level',nc_deflate_level))
        !
        NF90(nf90_put_att(ncid, varid, 'cdnrb', cd_nr))   
        NF90(nf90_put_att(ncid, varid, 'cdwnd', cd_wnd))        
        NF90(nf90_put_att(ncid, varid, 'cdval', cd_val))  
        !
        ! Internal code switches - note, you can't store logicals in netcdf, only integers for these type of switches
        NF90(nf90_put_att(ncid, varid, 'manning2d', logical2int(manning2d)))   
        NF90(nf90_put_att(ncid, varid, 'subgrid', logical2int(subgrid)))   
        NF90(nf90_put_att(ncid, varid, 'viscosity', logical2int(viscosity)))   
        NF90(nf90_put_att(ncid, varid, 'wavemaker', logical2int(wavemaker)))         
        NF90(nf90_put_att(ncid, varid, 'wavemaker_spectrum', logical2int(wavemaker_spectrum)))  
        NF90(nf90_put_att(ncid, varid, 'wavemaker_hig', logical2int(wavemaker_hig)))  
        NF90(nf90_put_att(ncid, varid, 'wavemaker_hinc', logical2int(wavemaker_hinc)))          
        !
        ! SnapWave related:
        !
        NF90(nf90_put_att(ncid, varid, 'snapwave', logical2int(snapwave))) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_gamma', gamma)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_gammax', gammax))         
        NF90(nf90_put_att(ncid, varid, 'snapwave_alpha', alpha)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_hmin',hmin)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_fw',fw0)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_fwig',fw0_ig)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_dt',dt)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_tol',tol)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_dtheta',dtheta)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_crit',crit)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_nrsweeps',nr_sweeps)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_baldock_exponent',baldock_exponent))
        NF90(nf90_put_att(ncid, varid, 'snapwave_baldock_ratio',baldock_ratio))
        NF90(nf90_put_att(ncid, varid, 'snapwave_waveforces_ratio',waveforces_ratio))
        NF90(nf90_put_att(ncid, varid, 'snapwave_sector',sector))
        NF90(nf90_put_att(ncid, varid, 'snapwave_Tpini',Tpini))
        NF90(nf90_put_att(ncid, varid, 'snapwave_fw_ratio',fwratio))
        NF90(nf90_put_att(ncid, varid, 'snapwave_fwig_ratio',fwigratio))
        !
        ! SnapWave IG
        !
        NF90(nf90_put_att(ncid, varid, 'snapwave_jonswapgamma',jonswapgam))
        NF90(nf90_put_att(ncid, varid, 'snapwave_igwaves',igwaves_opt))
        NF90(nf90_put_att(ncid, varid, 'snapwave_alpha_ig',alpha_ig)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_gammaig',gamma_ig))
        NF90(nf90_put_att(ncid, varid, 'snapwave_gamma_fac_br',gamma_fac_br))
        NF90(nf90_put_att(ncid, varid, 'snapwave_shinc2ig',shinc2ig))
        NF90(nf90_put_att(ncid, varid, 'snapwave_alphaigfac',alphaigfac)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_baldock_ratio_ig',baldock_ratio_ig)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_ig_opt',ig_opt)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_use_herbers',herbers_opt)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_tpig_opt',tpig_opt)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_eeinc2ig',eeinc2ig)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_Tinc2ig',Tinc2ig)) 
        !
        ! SnapWave input files
        !
        NF90(nf90_put_att(ncid, varid, 'snapwave_ncfile',gridfile))
        NF90(nf90_put_att(ncid, varid, 'snapwave_jonswapfile',snapwave_jonswapfile))
        NF90(nf90_put_att(ncid, varid, 'snapwave_encfile',snapwave_encfile))
        NF90(nf90_put_att(ncid, varid, 'snapwave_upwfile',upwfile))
        NF90(nf90_put_att(ncid, varid, 'snapwave_bndfile',snapwave_bndfile))            
        NF90(nf90_put_att(ncid, varid, 'snapwave_bhsfile',snapwave_bhsfile)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_btpfile',snapwave_btpfile)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_bwdfile',snapwave_bwdfile)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_bdsfile',snapwave_bdsfile))         
        !
   end subroutine
   !
   end module
