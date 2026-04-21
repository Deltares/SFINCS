#define NF90(nf90call) call handle_err(nf90call,__FILE__,__LINE__)
module sfincs_ncoutput
   !
   use netcdf
   use sfincs_nchelper
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
      integer :: zs_varid, zsmax_varid, h_varid, u_varid, v_varid, tmax_varid, Seff_varid, t_zsmax_varid
      integer :: zvolume_varid, storagevolume_varid
      integer :: hmax_varid, vmax_varid, qmax_varid, cumprcp_varid, cuminf_varid, windmax_varid
      integer :: patm_varid, wind_u_varid, wind_v_varid, precip_varid
      integer :: hm0_varid, hm0ig_varid, snapwavemsk_varid, tp_varid, tpig_varid, wavdir_varid, dirspr_varid
      integer :: fwx_varid, fwy_varid, beta_varid, snapwavedepth_varid
      integer :: zsm_varid, tsunami_arrival_time_varid, average_required_timestep_varid, percentage_limiting_varid
      integer :: inp_varid, total_runtime_varid, average_dt_varid, status_varid
      integer :: manning_varid
      integer :: pnonh_varid
      integer :: subgridslope_varid
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
      integer :: crosssections_dimid, structures_dimid, thindams_dimid, drain_dimid, runup_gauges_dimid
      integer :: runtime_dimid
      integer :: point_x_varid, point_y_varid, station_x_varid, station_y_varid, crs_varid, qinf_varid, S_varid
      integer :: station_id_varid, station_name_varid
      integer :: crosssection_name_varid
      integer :: structure_height_varid, structure_x_varid, structure_y_varid
      integer :: thindam_x_varid, thindam_y_varid
      integer :: drain_varid, drain_name_varid
      integer :: zb_varid
      integer :: time_varid
      integer :: zs_varid, h_varid, u_varid, v_varid, prcp_varid, cumprcp_varid, discharge_varid, uvmag_varid, uvdir_varid
      integer :: patm_varid, wind_speed_varid, wind_dir_varid
      integer :: inp_varid, total_runtime_varid, average_dt_varid, status_varid
      integer :: hm0_varid, hm0ig_varid, zsm_varid, tp_varid, tpig_varid, wavdir_varid, dirspr_varid
      integer :: dw_varid, df_varid, dwig_varid, dfig_varid, cg_varid, beta_varid, srcig_varid, alphaig_varid
      integer :: runup_gauge_name_varid, runup_gauge_zs_varid
      !
   end type
   !
   type(map_type) :: map_file
   type(his_type) :: his_file
   !
   real*4, parameter :: FILL_VALUE = -99999.0

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Cell-variable writers
   !
   ! These hide the use_quadtree branch and the scatter/gather pattern
   ! from ncoutput_update_map / ncoutput_update_max. A new time-varying
   ! cell-centered map variable becomes a one-line call.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine write_cell_var(ncid, varid, source, nt, use_sw_index, check_kcs, scale)
      !
      ! Scatter SFINCS cell-centered source(:) into the file-shaped buffer
      ! and put it to the NetCDF variable at timestep nt.
      !   use_sw_index=.true. -> quadtree gather uses index_sw_in_qt (snapwave)
      !   check_kcs=.true.    -> quadtree skips cells with kcs<=0 (max-file style)
      !   scale               -> multiply source by scale on write (unit conversions)
      !
      use sfincs_data
      use sfincs_snapwave, only: index_sw_in_qt
      use quadtree
      !
      integer, intent(in)           :: ncid, varid, nt
      real*4,  intent(in)           :: source(:)
      logical, intent(in), optional :: use_sw_index, check_kcs
      real*4,  intent(in), optional :: scale
      !
      real*4,  allocatable :: buf_q(:), buf_r(:,:)
      logical :: use_sw, filt_kcs
      real*4  :: fac
      integer :: nm, nmq
      !
      use_sw   = .false.; if (present(use_sw_index)) use_sw   = use_sw_index
      filt_kcs = .false.; if (present(check_kcs))    filt_kcs = check_kcs
      fac      = 1.0;     if (present(scale))        fac      = scale
      !
      if (use_quadtree) then
         allocate(buf_q(quadtree_nr_points))
         buf_q = FILL_VALUE
         do nmq = 1, quadtree_nr_points
            if (use_sw) then
               nm = index_sw_in_qt(nmq)
            else
               nm = index_sfincs_in_quadtree(nmq)
            endif
            if (nm > 0) then
               if (filt_kcs) then
                  if (kcs(nm) > 0) buf_q(nmq) = fac * source(nm)
               else
                  buf_q(nmq) = fac * source(nm)
               endif
            endif
         enddo
         NF90(nf90_put_var(ncid, varid, buf_q, (/1, nt/)))
      else
         allocate(buf_r(mmax, nmax))
         buf_r = FILL_VALUE
         do nm = 1, np
            buf_r(z_index_z_m(nm), z_index_z_n(nm)) = fac * source(nm)
         enddo
         NF90(nf90_put_var(ncid, varid, buf_r, (/1, 1, nt/)))
      endif
   end subroutine write_cell_var

   subroutine write_cell_var_wet(ncid, varid, source, zref, nt, output_delta, check_wet)
      !
      ! Wet-cell filtered writer for zs/h/zsmax/hmax.
      ! On quadtree: skips kcs<=0 cells. On both grids: when check_wet (default .true.),
      ! only writes cells with (source-zref) > huthresh. If output_delta, writes the
      ! delta instead of source.
      !
      use sfincs_data
      use quadtree
      !
      integer, intent(in)           :: ncid, varid, nt
      real*4,  intent(in)           :: source(:), zref(:)
      logical, intent(in)           :: output_delta
      logical, intent(in), optional :: check_wet
      !
      real*4,  allocatable :: buf_q(:), buf_r(:,:)
      logical :: filt_wet
      integer :: nm, nmq
      real*4  :: val
      !
      filt_wet = .true.; if (present(check_wet)) filt_wet = check_wet
      !
      if (use_quadtree) then
         allocate(buf_q(quadtree_nr_points))
         buf_q = FILL_VALUE
         do nmq = 1, quadtree_nr_points
            nm = index_sfincs_in_quadtree(nmq)
            if (nm > 0) then
               if (kcs(nm) > 0) then
                  if (.not. filt_wet .or. (source(nm) - zref(nm)) > huthresh) then
                     if (output_delta) then
                        buf_q(nmq) = source(nm) - zref(nm)
                     else
                        buf_q(nmq) = source(nm)
                     endif
                  endif
               endif
            endif
         enddo
         NF90(nf90_put_var(ncid, varid, buf_q, (/1, nt/)))
      else
         allocate(buf_r(mmax, nmax))
         buf_r = FILL_VALUE
         do nm = 1, np
            if (.not. filt_wet .or. (source(nm) - zref(nm)) > huthresh) then
               if (output_delta) then
                  val = source(nm) - zref(nm)
               else
                  val = source(nm)
               endif
               buf_r(z_index_z_m(nm), z_index_z_n(nm)) = val
            endif
         enddo
         NF90(nf90_put_var(ncid, varid, buf_r, (/1, 1, nt/)))
      endif
   end subroutine write_cell_var_wet

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
   ! Dimension abstraction
   integer :: nsd
   integer :: sp_dimids(2), dims_s(3), dims_st(3), dims_sm(3)
   character(4) :: coord_str
   !
   ! Regular-grid temporaries
   real*4, dimension(:,:), allocatable :: zsg, xz, yz, xg, yg
   !
   ! Quadtree temporaries
   real,      dimension(:),   allocatable :: nodes_x, nodes_y
   integer*4, dimension(:,:), allocatable :: face_nodes
   real*4,    dimension(:),   allocatable :: vtmp
   integer*4, dimension(:),   allocatable :: vtmpi
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
      allocate(vtmp(n_faces))
      allocate(vtmpi(n_faces))
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
   ! Build dim arrays
   dims_s (1:nsd)   = sp_dimids(1:nsd)
   dims_st(1:nsd)   = sp_dimids(1:nsd);  dims_st(nsd+1) = map_file%time_dimid
   dims_sm(1:nsd)   = sp_dimids(1:nsd);  dims_sm(nsd+1) = map_file%timemax_dimid
   !
   ! Global metadata
   if (use_quadtree) then
      NF90(nf90_put_att(map_file%ncid, nf90_global, "Conventions", "Conventions = 'CF-1.8 UGRID-1.0 Deltares-0.10'"))
   else
      NF90(nf90_put_att(map_file%ncid, nf90_global, "Conventions", "Conventions = 'CF-1.6, SGRID-0.3"))
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
      !
      if (crsgeo) then
         NF90(nf90_def_var(map_file%ncid, 'mesh2d_node_x', NF90_FLOAT, (/map_file%nmesh2d_node_dimid/), map_file%mesh2d_node_x_varid))
         NF90(nf90_def_var_deflate(map_file%ncid, map_file%mesh2d_node_x_varid, 1, 1, nc_deflate_level))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'units',         'degrees'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'standard_name', 'longitude'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'long_name',     'longitude'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'mesh',          'mesh2d'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'location',      'node'))
         !
         NF90(nf90_def_var(map_file%ncid, 'mesh2d_node_y', NF90_FLOAT, (/map_file%nmesh2d_node_dimid/), map_file%mesh2d_node_y_varid))
         NF90(nf90_def_var_deflate(map_file%ncid, map_file%mesh2d_node_y_varid, 1, 1, nc_deflate_level))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'units',         'degrees'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'standard_name', 'latitude'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'long_name',     'latitude'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'mesh',          'mesh2d'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'location',      'node'))
      else
         NF90(nf90_def_var(map_file%ncid, 'mesh2d_node_x', NF90_DOUBLE, (/map_file%nmesh2d_node_dimid/), map_file%mesh2d_node_x_varid))
         NF90(nf90_def_var_deflate(map_file%ncid, map_file%mesh2d_node_x_varid, 1, 1, nc_deflate_level))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'units',         'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'standard_name', 'projection_x_coordinate'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'long_name',     'x-coordinate of mesh nodes'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'mesh',          'mesh2d'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_x_varid, 'location',      'node'))
         !
         NF90(nf90_def_var(map_file%ncid, 'mesh2d_node_y', NF90_DOUBLE, (/map_file%nmesh2d_node_dimid/), map_file%mesh2d_node_y_varid))
         NF90(nf90_def_var_deflate(map_file%ncid, map_file%mesh2d_node_y_varid, 1, 1, nc_deflate_level))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'units',         'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'standard_name', 'projection_y_coordinate'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'long_name',     'y-coordinate of mesh nodes'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'mesh',          'mesh2d'))
         NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_node_y_varid, 'location',      'node'))
      endif
      !
      NF90(nf90_def_var(map_file%ncid, 'mesh2d_face_nodes', NF90_INT, &
           (/map_file%max_nmesh2d_face_nodes_dimid, map_file%nmesh2d_face_dimid/), map_file%mesh2d_face_nodes_varid))
      NF90(nf90_def_var_deflate(map_file%ncid, map_file%mesh2d_face_nodes_varid, 1, 1, nc_deflate_level))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'cf_role',   'face_node_connectivity'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'mesh',      'mesh2d'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'location',  'face'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'long_name', 'Mapping from every face to its corner nodes (counterclockwise)'))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, 'start_index', 1))
      NF90(nf90_put_att(map_file%ncid, map_file%mesh2d_face_nodes_varid, '_FillValue', -999))
      !
      NF90(nf90_def_var(map_file%ncid, 'crs', NF90_INT, map_file%crs_varid))
      NF90(nf90_put_att(map_file%ncid, map_file%crs_varid, 'EPSG', '-'))
      !
   else
      !
      ! Regular: x, y face coordinates
      NF90(nf90_def_var(map_file%ncid, 'x', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid/), map_file%face_x_varid))
      NF90(nf90_def_var_deflate(map_file%ncid, map_file%face_x_varid, 1, 1, nc_deflate_level))
      NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, '_FillValue', FILL_VALUE))
      if (crsgeo) then
         NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'units',         'degrees'))
         NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'standard_name', 'longitude'))
      else
         NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'units',         'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'standard_name', 'projection_x_coordinate'))
      endif
      NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'long_name',    'face_x'))
      NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'grid_mapping', 'crs'))
      NF90(nf90_put_att(map_file%ncid, map_file%face_x_varid, 'grid',         'sfincsgrid'))
      !
      NF90(nf90_def_var(map_file%ncid, 'y', NF90_FLOAT, (/map_file%m_dimid, map_file%n_dimid/), map_file%face_y_varid))
      NF90(nf90_def_var_deflate(map_file%ncid, map_file%face_y_varid, 1, 1, nc_deflate_level))
      NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, '_FillValue', FILL_VALUE))
      if (crsgeo) then
         NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'units',         'degrees'))
         NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'standard_name', 'latitude'))
      else
         NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'units',         'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'standard_name', 'projection_y_coordinate'))
      endif
      NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'long_name',    'face_y'))
      NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'grid_mapping', 'crs'))
      NF90(nf90_put_att(map_file%ncid, map_file%face_y_varid, 'grid',         'sfincsgrid'))
      !
      NF90(nf90_def_var(map_file%ncid, 'corner_x', NF90_FLOAT, (/map_file%corner_m_dimid, map_file%corner_n_dimid/), map_file%corner_x_varid))
      NF90(nf90_def_var_deflate(map_file%ncid, map_file%corner_x_varid, 1, 1, nc_deflate_level))
      NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, '_FillValue', FILL_VALUE))
      if (crsgeo) then
         NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'units',         'degrees'))
         NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'standard_name', 'longitude'))
      else
         NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'units',         'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'standard_name', 'projection_x_coordinate'))
      endif
      NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'long_name',    'corner_x'))
      NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'grid_mapping', 'crs'))
      NF90(nf90_put_att(map_file%ncid, map_file%corner_x_varid, 'grid',         'sfincsgrid'))
      !
      NF90(nf90_def_var(map_file%ncid, 'corner_y', NF90_FLOAT, (/map_file%corner_m_dimid, map_file%corner_n_dimid/), map_file%corner_y_varid))
      NF90(nf90_def_var_deflate(map_file%ncid, map_file%corner_y_varid, 1, 1, nc_deflate_level))
      NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, '_FillValue', FILL_VALUE))
      if (crsgeo) then
         NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'units',         'degrees'))
         NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'standard_name', 'latitude'))
      else
         NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'units',         'm'))
         NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'standard_name', 'projection_y_coordinate'))
      endif
      NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'long_name',    'corner_y'))
      NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'grid_mapping', 'crs'))
      NF90(nf90_put_att(map_file%ncid, map_file%corner_y_varid, 'grid',         'sfincsgrid'))
      !
      NF90(nf90_def_var(map_file%ncid, 'crs', NF90_INT, map_file%crs_varid))
      NF90(nf90_put_att(map_file%ncid, map_file%crs_varid, 'EPSG',      '-'))
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
   ! msk (INT+no-coord for quadtree, FLOAT+coord for regular)
   ! -------------------------------------------------------
   if (use_quadtree) then
      call ncdef_int_var(map_file%ncid, 'msk', dims_s(1:nsd), map_file%msk_varid, &
           'msk_active_cells', &
           description='inactive=0, active=1, normal_boundary=2, outflow_boundary=3, wavemaker=4')
      NF90(nf90_put_att(map_file%ncid, map_file%msk_varid, 'units',        '-'))
      NF90(nf90_put_att(map_file%ncid, map_file%msk_varid, 'standard_name','mask'))
   else
      call ncdef_float_var(map_file%ncid, 'msk', dims_s(1:nsd), map_file%msk_varid, &
           '-', 'msk_active_cells', standard_name='land_binary_mask', coordinates='x y', &
           description='inactive=0, active=1, normal_boundary=2, outflow_boundary=3')
   endif
   !
   ! -------------------------------------------------------
   ! Infiltration map (qinf) - both grid types, slightly different
   ! -------------------------------------------------------
   if (infiltration) then
      if (inftype == 'cna') then
         call ncdef_float_var(map_file%ncid, 'qinf', dims_s(1:nsd), map_file%qinf_varid, &
              'm', 'moisture storage (S) capacity - Curve number', &
              standard_name='S', coordinates=coord_str)
      elseif (inftype == 'cnb') then
         call ncdef_float_var(map_file%ncid, 'qinf', dims_s(1:nsd), map_file%qinf_varid, &
              'm', 'maximum moisture storage (Smax) capacity - Curve number', &
              standard_name='Smax', coordinates=coord_str)
      elseif (inftype == 'gai') then
         call ncdef_float_var(map_file%ncid, 'qinf', dims_s(1:nsd), map_file%qinf_varid, &
              'm', 'suction head at the wetting front - Green and Ampt', &
              standard_name='psi', coordinates=coord_str)
      elseif (inftype == 'hor') then
         call ncdef_float_var(map_file%ncid, 'qinf', dims_s(1:nsd), map_file%qinf_varid, &
              'm', 'initial infiltration rate - Horton', &
              standard_name='f0', coordinates=coord_str)
      else
         call ncdef_float_var(map_file%ncid, 'qinf', dims_s(1:nsd), map_file%qinf_varid, &
              'mm h-1', 'infiltration rate - constant in time', &
              standard_name='qinf', coordinates=coord_str)
      endif
   endif
   !
   ! -------------------------------------------------------
   ! zb: always static for quadtree; optional time-varying for regular
   ! -------------------------------------------------------
   if (use_quadtree) then
      call ncdef_float_var(map_file%ncid, 'zb', dims_s(1:nsd), map_file%zb_varid, &
           'm', 'bed_level_above_reference_level', standard_name='altitude')
   else
      if (store_dynamic_bed_level) then
         call ncdef_float_var(map_file%ncid, 'zb', dims_st(1:nsd+1), map_file%zb_varid, &
              'm', 'bed_level_above_reference_level', standard_name='altitude', coordinates='x y')
      else
         call ncdef_float_var(map_file%ncid, 'zb', dims_s(1:nsd), map_file%zb_varid, &
              'm', 'bed_level_above_reference_level', standard_name='altitude', coordinates='x y')
      endif
   endif
   !
   ! manning (both grid types, dims_s)
   if (.not. subgrid) then
      call ncdef_float_var(map_file%ncid, 'manning', dims_s(1:nsd), map_file%manning_varid, &
           's/m^1/3', 'manning_roughness', standard_name='manning', coordinates=coord_str)
   endif
   !
   ! subgrid slope
   if (subgrid .and. store_hsubgrid .and. store_hmean) then
      call ncdef_float_var(map_file%ncid, 'subgridslope', dims_s(1:nsd), map_file%subgridslope_varid, &
           '-', 'subgrid_slope', standard_name='subgrid_slope', coordinates=coord_str)
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
   call ncdef_float_var(map_file%ncid, 'zs', dims_st(1:nsd+1), map_file%zs_varid, &
        'm', 'water_level', &
        standard_name='sea_surface_height_above_reference_level', coordinates=coord_str)
   !
   if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then
      call ncdef_float_var(map_file%ncid, 'h', dims_st(1:nsd+1), map_file%h_varid, &
           'm', 'water_depth', &
           standard_name='depth', coordinates=coord_str)
   endif
   !
   if (store_velocity) then
      call ncdef_float_var(map_file%ncid, 'u', dims_st(1:nsd+1), map_file%u_varid, &
           'm s-1', 'flow_velocity_x_direction', &
           standard_name='eastward_sea_water_velocity', coordinates=coord_str)
      call ncdef_float_var(map_file%ncid, 'v', dims_st(1:nsd+1), map_file%v_varid, &
           'm s-1', 'flow_velocity_y_direction', &
           standard_name='northward_sea_water_velocity', coordinates=coord_str)
   endif
   !
   ! Subgrid volumes
   if (subgrid) then
      if (store_zvolume) then
         call ncdef_float_var(map_file%ncid, 'subgrid_volume', dims_st(1:nsd+1), map_file%zvolume_varid, &
              'm3', 'subgrid_volume_in_cell', &
              standard_name='subgrid_volume_in_cell', coordinates=coord_str)
      endif
      if (store_storagevolume) then
         call ncdef_float_var(map_file%ncid, 'storage_volume', dims_st(1:nsd+1), map_file%storagevolume_varid, &
              'm3', 'storage_volume_in_cell', &
              standard_name='storage_volume_in_cell', coordinates=coord_str)
      endif
   endif
   !
   ! Infiltration state vars (regular only: Seff/sigma/f)
   if (.not. use_quadtree) then
      if (inftype == 'cnb') then
         call ncdef_float_var(map_file%ncid, 'Seff', dims_st(1:nsd+1), map_file%Seff_varid, &
              'm', 'current moisture storage (Se) capacity', &
              standard_name='Se', coordinates='corner_x corner_y')
      elseif (inftype == 'gai') then
         call ncdef_float_var(map_file%ncid, 'sigma', dims_st(1:nsd+1), map_file%Seff_varid, &
              '-', 'maximum soil moisture deficit', &
              standard_name='sigma', coordinates='corner_x corner_y')
      elseif (inftype == 'hor') then
         call ncdef_float_var(map_file%ncid, 'f', dims_st(1:nsd+1), map_file%Seff_varid, &
              'mm h-1', 'current infiltration capacity', &
              standard_name='sigma', coordinates='corner_x corner_y')
      endif
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
      call ncdef_float_var(map_file%ncid, 'zsmax', dims_sm(1:nsd+1), map_file%zsmax_varid, &
           'm', 'maximum_water_level', &
           standard_name='maximum of sea_surface_height_above_reference_level', coordinates=coord_str)
   endif
   !
   if (store_cumulative_precipitation) then
      call ncdef_float_var(map_file%ncid, 'cumprcp', dims_sm(1:nsd+1), map_file%cumprcp_varid, &
           'm', 'cumulative_precipitation_depth', &
           standard_name='cumulative_precipitation_depth', &
           cell_methods='time: sum', coordinates=coord_str)
   endif
   !
   if (store_twet) then
      call ncdef_float_var(map_file%ncid, 'tmax', dims_sm(1:nsd+1), map_file%tmax_varid, &
           'seconds', 'duration_wet_cell', &
           standard_name='duration cell is considered wet', &
           cell_methods='time: sum', coordinates=coord_str)
   endif
   !
   if (store_t_zsmax) then
      call ncdef_float_var(map_file%ncid, 't_zsmax', dims_sm(1:nsd+1), map_file%t_zsmax_varid, &
           'seconds since ' // trim(trefstr_iso8601), 'time when zsmax occurs', &
           standard_name='time of max water level', &
           cell_methods='time: max', coordinates=coord_str)
   endif
   !
   if (store_maximum_waterlevel) then
      if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then
         call ncdef_float_var(map_file%ncid, 'hmax', dims_sm(1:nsd+1), map_file%hmax_varid, &
              'm', 'maximum_water_depth', &
              standard_name='sea_floor_depth_below_sea_surface', &
              cell_methods='time: maximum', coordinates=coord_str)
      endif
   endif
   !
   if (store_maximum_velocity) then
      call ncdef_float_var(map_file%ncid, 'vmax', dims_sm(1:nsd+1), map_file%vmax_varid, &
           'm s-1', 'maximum_flow_velocity', &
           standard_name='maximum_flow_velocity', &
           cell_methods='time: maximum', coordinates=coord_str)
   endif
   !
   if (store_maximum_flux) then
      call ncdef_float_var(map_file%ncid, 'qmax', dims_sm(1:nsd+1), map_file%qmax_varid, &
           'm^2 s-1', 'maximum_flux', &
           standard_name='maximum_flux', &
           cell_methods='time: maximum', coordinates=coord_str)
   endif
   !
   if (timestep_analysis) then
      call ncdef_float_var(map_file%ncid, 'average_required_timestep', dims_s(1:nsd), map_file%average_required_timestep_varid, &
           's', 'average_required_timestep', &
           standard_name='average_required_timestep', &
           cell_methods='time: average', coordinates=coord_str)
      call ncdef_float_var(map_file%ncid, 'percentage_limiting_timestep', dims_s(1:nsd), map_file%percentage_limiting_varid, &
           '%', 'percentage_cell_was_limiting_timestep', &
           standard_name='percentage_limiting_timestep', &
           cell_methods='time: maximum', coordinates=coord_str)
   endif
   !
   if (store_cumulative_precipitation) then
      call ncdef_float_var(map_file%ncid, 'cuminf', dims_sm(1:nsd+1), map_file%cuminf_varid, &
           'm', 'cumulative_infiltration_depth', &
           cell_methods='time: sum', coordinates=coord_str)
   endif
   !
   ! -------------------------------------------------------
   ! Meteo variables
   ! -------------------------------------------------------
   if (store_meteo) then
      if (wind) then
         call ncdef_float_var(map_file%ncid, 'wind_u', dims_st(1:nsd+1), map_file%wind_u_varid, &
              'm s-1', 'wind_speed_u', &
              standard_name='eastward_wind', coordinates=coord_str)
         call ncdef_float_var(map_file%ncid, 'wind_v', dims_st(1:nsd+1), map_file%wind_v_varid, &
              'm s-1', 'wind_speed_v', &
              standard_name='northward_wind', coordinates=coord_str)
         !
         ! windmax: different dims and condition for regular vs quadtree
         if (use_quadtree) then
            if (store_wind_max) then
               call ncdef_float_var(map_file%ncid, 'windmax', dims_sm(1:nsd+1), map_file%windmax_varid, &
                    'm s-1', 'maximum_wind_speed', &
                    cell_methods='time: maximum')
            endif
         else
            if (meteo3d .and. store_wind_max) then
               call ncdef_float_var(map_file%ncid, 'windmax', dims_s(1:nsd), map_file%windmax_varid, &
                    'm s-1', 'maximum_wind_speed', &
                    cell_methods='time: sum', coordinates='x y')
            endif
         endif
      endif
      !
      if (patmos) then
         call ncdef_float_var(map_file%ncid, 'surface_air_pressure', dims_st(1:nsd+1), map_file%patm_varid, &
              'N m-2', 'surface_air_pressure', &
              standard_name='surface_air_pressure', coordinates=coord_str)
      endif
      !
      ! precipitation_rate: regular only
      if (.not. use_quadtree) then
         if (precip) then
            call ncdef_float_var(map_file%ncid, 'precipitation_rate', dims_st(1:nsd+1), map_file%precip_varid, &
                 'mm h-1', 'precipitation_rate', &
                 standard_name='precipitation_rate', coordinates=coord_str)
         endif
      endif
   endif
   !
   ! -------------------------------------------------------
   ! SnapWave variables
   ! -------------------------------------------------------
   if (snapwave) then
      !
      ! snapwavemsk: INT no-coord (quadtree), FLOAT+coord (regular)
      if (use_quadtree) then
         call ncdef_int_var(map_file%ncid, 'snapwavemsk', dims_s(1:nsd), map_file%snapwavemsk_varid, &
              'snapwave_msk_active_cells', &
              description='inactive=0, active=1, wave_boundary=2, neumann_boundary=3')
         NF90(nf90_put_att(map_file%ncid, map_file%snapwavemsk_varid, 'units',        '-'))
         NF90(nf90_put_att(map_file%ncid, map_file%snapwavemsk_varid, 'standard_name','snapwavemask'))
      else
         call ncdef_float_var(map_file%ncid, 'snapwavemsk', dims_s(1:nsd), map_file%snapwavemsk_varid, &
              '-', 'snapwave_msk_active_cells', standard_name='snapwavemask', coordinates='x y', &
              description='inactive=0, active=1, wave_boundary=2, neumann_boundary=3')
      endif
      !
      call ncdef_float_var(map_file%ncid, 'hm0', dims_st(1:nsd+1), map_file%hm0_varid, &
           'm', 'Hm0 wave height', standard_name='hm0_wave_height', coordinates=coord_str)
      call ncdef_float_var(map_file%ncid, 'hm0ig', dims_st(1:nsd+1), map_file%hm0ig_varid, &
           'm', 'Hm0 infragravity wave height', standard_name='hm0_ig_wave_height', coordinates=coord_str)
      !
      if (store_wave_forces) then
         call ncdef_float_var(map_file%ncid, 'fwx', dims_st(1:nsd+1), map_file%fwx_varid, &
              'm', 'Wave force in x-direction', standard_name='wave_force_x', coordinates=coord_str)
         call ncdef_float_var(map_file%ncid, 'fwy', dims_st(1:nsd+1), map_file%fwy_varid, &
              'm', 'Wave force in y-direction', standard_name='wave_force_y', coordinates=coord_str)
         call ncdef_float_var(map_file%ncid, 'tp', dims_st(1:nsd+1), map_file%tp_varid, &
              's', 'Peak wave period', standard_name='peak_wave_period', coordinates=coord_str)
         call ncdef_float_var(map_file%ncid, 'tpig', dims_st(1:nsd+1), map_file%tpig_varid, &
              's', 'Peak infragravity wave period', standard_name='peak_ig_wave_period', coordinates=coord_str)
         call ncdef_float_var(map_file%ncid, 'beta', dims_st(1:nsd+1), map_file%beta_varid, &
              '-', 'directionally averaged local bed slope', &
              standard_name='directionally_averaged_local_bed_slope', coordinates=coord_str)
         call ncdef_float_var(map_file%ncid, 'snapwavedepth', dims_st(1:nsd+1), map_file%snapwavedepth_varid, &
              'm', 'Interpolated water depth in Snapwave', &
              standard_name='snapwave_waterdepth', coordinates=coord_str)
      endif
      !
      ! wavdir: quadtree only
      if (use_quadtree .and. store_wave_direction) then
         call ncdef_float_var(map_file%ncid, 'wavdir', dims_st(1:nsd+1), map_file%wavdir_varid, &
              'degrees', 'Mean wave direction', standard_name='mean_wave_direction')
      endif
      !
      if (wavemaker) then
         call ncdef_float_var(map_file%ncid, 'zsm', dims_st(1:nsd+1), map_file%zsm_varid, &
              'm', 'Filtered water level', &
              standard_name='mean_water_level', coordinates=coord_str)
      endif
      !
   endif
   !
   ! tsunami_arrival_time: time-varying for quadtree, static for regular
   if (store_tsunami_arrival_time) then
      if (use_quadtree) then
         call ncdef_float_var(map_file%ncid, 'tsunami_arrival_time', dims_st(1:nsd+1), &
              map_file%tsunami_arrival_time_varid, '-', 'tsunami_arrival_time', &
              standard_name='tsunami_arrival_time')
      else
         call ncdef_float_var(map_file%ncid, 'tsunami_arrival_time', dims_s(1:nsd), &
              map_file%tsunami_arrival_time_varid, '-', 'tsunami_arrival_time', &
              standard_name='tsunami_arrival_time', coordinates='x y')
      endif
   endif
   !
   if (nonhydrostatic) then
      call ncdef_float_var(map_file%ncid, 'pnonh', dims_st(1:nsd+1), map_file%pnonh_varid, &
           'N m-2', 'non_hydrostatic_pressure', coordinates=coord_str)
   endif
   !
   ! Runtime scalars
   call ncdef_float_var(map_file%ncid, 'total_runtime', (/map_file%runtime_dimid/), map_file%total_runtime_varid, &
        's', 'total_model_runtime_in_seconds')
   call ncdef_float_var(map_file%ncid, 'average_dt', (/map_file%runtime_dimid/), map_file%average_dt_varid, &
        's', 'model_average_timestep_in_seconds')
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
   if (use_quadtree) then
      !
      NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_node_x_varid, nodes_x))
      NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_node_y_varid, nodes_y))
      NF90(nf90_put_var(map_file%ncid, map_file%mesh2d_face_nodes_varid, face_nodes))
      NF90(nf90_put_var(map_file%ncid, map_file%crs_varid, epsg))
      !
      ! zb
      vtmp = FILL_VALUE
      if (subgrid) then
         do nmq = 1, quadtree_nr_points
            nm = index_sfincs_in_quadtree(nmq)
            if (nm > 0) vtmp(nmq) = subgrid_z_zmin(nm)
         enddo
      else
         do nmq = 1, quadtree_nr_points
            nm = index_sfincs_in_quadtree(nmq)
            if (nm > 0) vtmp(nmq) = zb(nm)
         enddo
      endif
      NF90(nf90_put_var(map_file%ncid, map_file%zb_varid, vtmp))
      !
      ! subgrid slope
      if (subgrid .and. store_hsubgrid .and. store_hmean) then
         vtmp = FILL_VALUE
         do nmq = 1, quadtree_nr_points
            nm = index_sfincs_in_quadtree(nmq)
            if (nm > 0) then
               if (crsgeo) then
                  vtmp(nmq) = (subgrid_z_zmax(nm) - subgrid_z_zmin(nm)) / sqrt(cell_area_m2(nm))
               else
                  vtmp(nmq) = (subgrid_z_zmax(nm) - subgrid_z_zmin(nm)) / sqrt(cell_area(z_flags_iref(nm)))
               endif
            endif
         enddo
         NF90(nf90_put_var(map_file%ncid, map_file%subgridslope_varid, vtmp))
      endif
      !
      ! msk
      vtmpi = 0
      do nmq = 1, quadtree_nr_points
         nm = index_sfincs_in_quadtree(nmq)
         if (nm > 0) vtmpi(nmq) = kcs(nm)
      enddo
      NF90(nf90_put_var(map_file%ncid, map_file%msk_varid, vtmpi))
      !
      ! snapwave msk
      if (snapwave) then
         vtmpi = 0
         do nmq = 1, quadtree_nr_points
            nm = index_sw_in_qt(nmq)
            if (nm > 0) vtmpi(nmq) = snapwave_mask(nm)
         enddo
         NF90(nf90_put_var(map_file%ncid, map_file%snapwavemsk_varid, vtmpi))
      endif
      !
      ! Manning
      if (.not. subgrid .and. manning2d) then
         vtmp = FILL_VALUE
         do nmq = 1, quadtree_nr_points
            nm = index_sfincs_in_quadtree(nmq)
            if (nm > 0) vtmp(nmq) = rghfield(nm)
         enddo
         NF90(nf90_put_var(map_file%ncid, map_file%manning_varid, vtmp))
      endif
      !
      ! Infiltration map
      if (infiltration) then
         vtmp = FILL_VALUE
         if (inftype == 'con' .or. inftype == 'c2d') then
            do nmq = 1, quadtree_nr_points
               nm = index_sfincs_in_quadtree(nmq)
               if (nm > 0) vtmp(nmq) = qinffield(nm) * 3600 * 1000
            enddo
         else
            do nmq = 1, quadtree_nr_points
               nm = index_sfincs_in_quadtree(nmq)
               if (nm > 0) vtmp(nmq) = qinffield(nm)
            enddo
         endif
         NF90(nf90_put_var(map_file%ncid, map_file%qinf_varid, vtmp))
      endif
      !
   else
      !
      ! Regular grid
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
      do n = 1, nmax + 1
         do m = 1, mmax + 1
            xg(m, n) = x0 + cosrot*(1.0*(m - 1))*dx - sinrot*(1.0*(n - 1))*dy
            yg(m, n) = y0 + sinrot*(1.0*(m - 1))*dx + cosrot*(1.0*(n - 1))*dy
         enddo
      enddo
      !
      NF90(nf90_put_var(map_file%ncid, map_file%face_x_varid,   xz, (/1, 1/)))
      NF90(nf90_put_var(map_file%ncid, map_file%face_y_varid,   yz, (/1, 1/)))
      NF90(nf90_put_var(map_file%ncid, map_file%corner_x_varid, xg, (/1, 1/)))
      NF90(nf90_put_var(map_file%ncid, map_file%corner_y_varid, yg, (/1, 1/)))
      NF90(nf90_put_var(map_file%ncid, map_file%crs_varid, epsg))
      !
      allocate(zsg(mmax, nmax))
      zsg = FILL_VALUE
      !
      ! zb static
      if (.not. store_dynamic_bed_level) then
         do nm = 1, np
            n = z_index_z_n(nm)
            m = z_index_z_m(nm)
            if (subgrid) then
               zsg(m, n) = subgrid_z_zmin(nm)
            else
               zsg(m, n) = zb(nm)
            endif
         enddo
         NF90(nf90_put_var(map_file%ncid, map_file%zb_varid, zsg, (/1, 1/)))
      endif
      !
      zsg = FILL_VALUE
      ! subgrid slope
      if (subgrid .and. store_hsubgrid .and. store_hmean) then
         do nm = 1, np
            n = z_index_z_n(nm)
            m = z_index_z_m(nm)
            if (crsgeo) then
               zsg(m, n) = (subgrid_z_zmax(nm) - subgrid_z_zmin(nm)) / sqrt(cell_area_m2(nm))
            else
               zsg(m, n) = (subgrid_z_zmax(nm) - subgrid_z_zmin(nm)) / sqrt(cell_area(z_flags_iref(nm)))
            endif
         enddo
         NF90(nf90_put_var(map_file%ncid, map_file%subgridslope_varid, zsg, (/1, 1/)))
      endif
      !
      ! msk
      zsg = 0
      do nm = 1, np
         n = z_index_z_n(nm)
         m = z_index_z_m(nm)
         zsg(m, n) = kcs(nm)
      enddo
      NF90(nf90_put_var(map_file%ncid, map_file%msk_varid, zsg, (/1, 1/)))
      !
      ! snapwave msk
      if (snapwave) then
         zsg = 0
         do nm = 1, np
            n = z_index_z_n(nm)
            m = z_index_z_m(nm)
            zsg(m, n) = snapwave_mask(nm)
         enddo
         NF90(nf90_put_var(map_file%ncid, map_file%snapwavemsk_varid, zsg, (/1, 1/)))
      endif
      !
      ! Manning
      if (.not. subgrid .and. manning2d) then
         zsg = FILL_VALUE
         do nm = 1, np
            n = z_index_z_n(nm)
            m = z_index_z_m(nm)
            zsg(m, n) = rghfield(nm)
         enddo
         NF90(nf90_put_var(map_file%ncid, map_file%manning_varid, zsg, (/1, 1/)))
      endif
      !
      ! Infiltration map
      if (infiltration) then
         zsg = FILL_VALUE
         do nm = 1, np
            n = z_index_z_n(nm)
            m = z_index_z_m(nm)
            if (inftype == 'con' .or. inftype == 'c2d') then
               zsg(m, n) = qinffield(nm) * 3600 * 1000
            else
               zsg(m, n) = qinffield(nm)
            endif
         enddo
         NF90(nf90_put_var(map_file%ncid, map_file%qinf_varid, zsg, (/1, 1/)))
      endif
      !
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
   !
   implicit none
   !
   integer                      :: istruc
   character(*), parameter :: pt_coord = 'station_id station_name point_x point_y'
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
   if (nobs==0 .and. nrcrosssections==0 .and. nrstructures==0 .and. ndrn==0 .and. nr_runup_gauges==0) then ! If no observation points, cross-sections, structures, drains or run-up gauges; his file is not created        
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
   if (ndrn>0) then   
      NF90(nf90_def_dim(his_file%ncid, 'drainage', ndrn, his_file%drain_dimid)) ! nr of drainage structures
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
   if (nr_runup_gauges > 0) then
      NF90(nf90_def_var(his_file%ncid, 'runup_gauge_name', NF90_CHAR, (/his_file%pointnamelength_dimid, his_file%runup_gauges_dimid/), his_file%runup_gauge_name_varid))
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
   call ncdef_float_var(his_file%ncid, 'point_zb', (/his_file%points_dimid/), his_file%zb_varid, &
        'm', 'bed_level_above_reference_level', standard_name='altitude', coordinates=pt_coord)
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
   call ncdef_float_var(his_file%ncid, 'point_zs', (/his_file%points_dimid, his_file%time_dimid/), his_file%zs_varid, &
        'm', 'water_level', standard_name='sea_surface_height_above_reference_level', coordinates=pt_coord)
   !
   if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then
      call ncdef_float_var(his_file%ncid, 'point_h', (/his_file%points_dimid, his_file%time_dimid/), his_file%h_varid, &
           'm', 'water_depth', standard_name='depth', coordinates=pt_coord)
   endif
   !
   if (store_velocity) then
      call ncdef_float_var(his_file%ncid, 'point_u', (/his_file%points_dimid, his_file%time_dimid/), his_file%u_varid, &
           'm s-1', 'flow_velocity_x_direction', standard_name='sea_water_x_velocity', coordinates=pt_coord)
      call ncdef_float_var(his_file%ncid, 'point_v', (/his_file%points_dimid, his_file%time_dimid/), his_file%v_varid, &
           'm s-1', 'flow_velocity_y_direction', standard_name='sea_water_y_velocity', coordinates=pt_coord)
      call ncdef_float_var(his_file%ncid, 'point_uvmag', (/his_file%points_dimid, his_file%time_dimid/), his_file%uvmag_varid, &
           'm s-1', 'flow_velocity_magnitude', standard_name='sea_water_velocity', coordinates=pt_coord)
      call ncdef_float_var(his_file%ncid, 'point_uvdir', (/his_file%points_dimid, his_file%time_dimid/), his_file%uvdir_varid, &
           'degrees', 'flow_velocity_direction', standard_name='sea_water_velocity_direction', coordinates=pt_coord)
   endif
   !
   ! Add infiltration
   !
   if (infiltration) then
      call ncdef_float_var(his_file%ncid, 'point_qinf', (/his_file%points_dimid, his_file%time_dimid/), his_file%qinf_varid, &
           'mm hr-1', 'infiltration_rate', coordinates=pt_coord)
   endif
   !
   if (inftype == 'cnb') then
      call ncdef_float_var(his_file%ncid, 'point_S', (/his_file%points_dimid, his_file%time_dimid/), his_file%S_varid, &
           'm', 'current moisture storage (Se) capacity', coordinates=pt_coord)
   endif
   !
   if (inftype == 'gai') then
      call ncdef_float_var(his_file%ncid, 'point_S', (/his_file%points_dimid, his_file%time_dimid/), his_file%S_varid, &
           'm', 'maximum soil moisture deficit', coordinates=pt_coord)
   endif
   !
   if (snapwave) then
      call ncdef_float_var(his_file%ncid, 'point_hm0', (/his_file%points_dimid, his_file%time_dimid/), his_file%hm0_varid, &
           'm', 'Hm0 wave height', standard_name='hm0_wave_height', coordinates=pt_coord)
      call ncdef_float_var(his_file%ncid, 'point_hm0ig', (/his_file%points_dimid, his_file%time_dimid/), his_file%hm0ig_varid, &
           'm', 'Hm0 infragravity wave height', standard_name='hm0_ig_wave_height', coordinates=pt_coord)
      call ncdef_float_var(his_file%ncid, 'point_tp', (/his_file%points_dimid, his_file%time_dimid/), his_file%tp_varid, &
           's', 'Peak wave period', standard_name='peak_wave_period', coordinates=pt_coord)
      call ncdef_float_var(his_file%ncid, 'point_tpig', (/his_file%points_dimid, his_file%time_dimid/), his_file%tpig_varid, &
           's', 'Peak wave period Infragravity wave', standard_name='ig_peak_wave_period', coordinates=pt_coord)
      !
      if (store_wave_direction) then
         call ncdef_float_var(his_file%ncid, 'point_wavdir', (/his_file%points_dimid, his_file%time_dimid/), his_file%wavdir_varid, &
              'degrees', 'Mean wave direction', standard_name='mean_wave_direction', coordinates=pt_coord)
         call ncdef_float_var(his_file%ncid, 'point_dirspr', (/his_file%points_dimid, his_file%time_dimid/), his_file%dirspr_varid, &
              's', 'Wave directional spreading', standard_name='wave_directional_spreading', coordinates=pt_coord)
      endif
      !
      if (wavemaker) then
         call ncdef_float_var(his_file%ncid, 'point_zsm', (/his_file%points_dimid, his_file%time_dimid/), his_file%zsm_varid, &
              'm', 'Filtered water level', standard_name='filtered_water_level', coordinates=pt_coord)
      endif
      !
      if (store_wave_forces) then
         call ncdef_float_var(his_file%ncid, 'point_dw', (/his_file%points_dimid, his_file%time_dimid/), his_file%dw_varid, &
              'm', 'directionally averaged wave breaking dissipation', &
              standard_name='directionally_averaged_wave_breaking_dissipation', coordinates=pt_coord)
         call ncdef_float_var(his_file%ncid, 'point_df', (/his_file%points_dimid, his_file%time_dimid/), his_file%df_varid, &
              'm', 'directionally averaged wave friction dissipation', &
              standard_name='directionally_averaged_wave_friction_dissipation', coordinates=pt_coord)
         call ncdef_float_var(his_file%ncid, 'point_dwig', (/his_file%points_dimid, his_file%time_dimid/), his_file%dwig_varid, &
              'm', 'directionally averaged wave breaking dissipation ig', &
              standard_name='directionally_averaged_wave_breaking_dissipation_ig', coordinates=pt_coord)
         call ncdef_float_var(his_file%ncid, 'point_dfig', (/his_file%points_dimid, his_file%time_dimid/), his_file%dfig_varid, &
              'm', 'directionally averaged wave friction dissipation ig', &
              standard_name='directionally_averaged_wave_friction_dissipation_ig', coordinates=pt_coord)
         call ncdef_float_var(his_file%ncid, 'point_cg', (/his_file%points_dimid, his_file%time_dimid/), his_file%cg_varid, &
              'm/s', 'wave group velocity', standard_name='wave_group_velocity', coordinates=pt_coord)
         call ncdef_float_var(his_file%ncid, 'point_beta', (/his_file%points_dimid, his_file%time_dimid/), his_file%beta_varid, &
              '-', 'directionally averaged normalised bed slope', &
              standard_name='directionally_averaged_local_bed_slope', coordinates=pt_coord)
         call ncdef_float_var(his_file%ncid, 'point_srcig', (/his_file%points_dimid, his_file%time_dimid/), his_file%srcig_varid, &
              '-', 'directionally averaged ig energy source', &
              standard_name='directionally_averaged_ig_energy_source', coordinates=pt_coord)
         call ncdef_float_var(his_file%ncid, 'point_alphaig', (/his_file%points_dimid, his_file%time_dimid/), his_file%alphaig_varid, &
              '-', 'directionally averaged infragravity waves shoaling factor', &
              standard_name='directionally_averaged_infragravity_waves_shoaling_factor', coordinates=pt_coord)
      endif
   endif
   !
   if (store_meteo) then
      if (wind) then
         call ncdef_float_var(his_file%ncid, 'point_wind_speed', (/his_file%points_dimid, his_file%time_dimid/), his_file%wind_speed_varid, &
              'm s-1', 'wind_speed', coordinates=pt_coord)
         call ncdef_float_var(his_file%ncid, 'point_wind_direction', (/his_file%points_dimid, his_file%time_dimid/), his_file%wind_dir_varid, &
              'degrees', 'wind_direction', coordinates=pt_coord)
      endif
      if (patmos) then
         call ncdef_float_var(his_file%ncid, 'point_patm', (/his_file%points_dimid, his_file%time_dimid/), his_file%patm_varid, &
              'Pa', 'surface_air_pressure', coordinates=pt_coord)
      endif
      if (precip) then
         call ncdef_float_var(his_file%ncid, 'point_prcp', (/his_file%points_dimid, his_file%time_dimid/), his_file%prcp_varid, &
              'mm hr-1', 'precipitation_rate', coordinates=pt_coord)
         if (store_cumulative_precipitation) then
            call ncdef_float_var(his_file%ncid, 'point_cumprcp', (/his_file%points_dimid, his_file%time_dimid/), his_file%cumprcp_varid, &
                 'm', 'cumulative_precipitation', coordinates=pt_coord)
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
           'm3 s-1', 'discharge through drainage structure', coordinates='drainage_name')
   endif
   !
   if (nr_runup_gauges > 0) then
      call ncdef_float_var(his_file%ncid, 'runup_gauge_zs', (/his_file%runup_gauges_dimid, his_file%time_dimid/), his_file%runup_gauge_zs_varid, &
           'm', 'run-up elevation', coordinates='runup_gauge_name')
   endif
   !
   call ncdef_float_var(his_file%ncid, 'total_runtime', (/his_file%runtime_dimid/), his_file%total_runtime_varid, &
        's', 'total_model_runtime_in_seconds')
   call ncdef_float_var(his_file%ncid, 'average_dt', (/his_file%runtime_dimid/), his_file%average_dt_varid, &
        's', 'model_average_timestep_in_seconds')
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
   end subroutine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   subroutine ncoutput_update_map(t,ntmapout)
   !
   ! Write time-varying output to map file (handles both regular and quadtree grids)
   !
   use sfincs_data
   use sfincs_nonhydrostatic
   use sfincs_snapwave
   use quadtree
   !
   implicit none
   !
   real*8                       :: t
   real*4                       :: uz, vz, sq2
   !
   integer  :: ntmapout
   !
   integer                      :: nm, nmq, n, m, nmd1, nmu1, ndm1, num1
   real*4, dimension(:,:), allocatable :: zsg, zsgu, zsgv
   real*4, dimension(:),   allocatable :: utmp, vtmp
   !
   NF90(nf90_put_var(map_file%ncid, map_file%time_varid, t, (/ntmapout/))) ! write time
   !
   if (use_quadtree) then
      !
      allocate(utmp(quadtree_nr_points))
      allocate(vtmp(quadtree_nr_points))
      !
      sq2 = sqrt(2.0)
      !
      ! zs: wet-cell filter, write zs itself
      if (subgrid) then
         call write_cell_var_wet(map_file%ncid, map_file%zs_varid, real(zs,4), subgrid_z_zmin, ntmapout, output_delta=.false.)
      else
         call write_cell_var_wet(map_file%ncid, map_file%zs_varid, real(zs,4), zb, ntmapout, output_delta=.false.)
      endif
      !
      ! h: wet-cell filter, write delta (zs - zref)
      if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then
         if (subgrid) then
            call write_cell_var_wet(map_file%ncid, map_file%h_varid, real(zs,4), subgrid_z_zmin, ntmapout, output_delta=.true.)
         else
            call write_cell_var_wet(map_file%ncid, map_file%h_varid, real(zs,4), zb, ntmapout, output_delta=.true.)
         endif
      endif
      !
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
               ! Regular point with four surrounding cells of the same size
               !
               n = z_index_z_n(nm)
               m = z_index_z_m(nm)
               nmd1 = z_index_uv_md(nm)
               nmu1 = z_index_uv_mu(nm)
               ndm1 = z_index_uv_nd(nm)
               num1 = z_index_uv_nu(nm)
               uz = 0.5*(uv(nmd1) + uv(nmu1))
               vz = 0.5*(uv(ndm1) + uv(num1))
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
      if (subgrid) then
         if (store_zvolume) then
            call write_cell_var(map_file%ncid, map_file%zvolume_varid, real(z_volume,4), ntmapout)
         endif
         if (store_storagevolume) then
            call write_cell_var(map_file%ncid, map_file%storagevolume_varid, storage_volume, ntmapout)
         endif
      endif
      !
      if (store_meteo) then
         if (wind) then
            call write_cell_var(map_file%ncid, map_file%wind_u_varid, windu, ntmapout)
            call write_cell_var(map_file%ncid, map_file%wind_v_varid, windv, ntmapout)
         endif
         if (patmos) then
            call write_cell_var(map_file%ncid, map_file%patm_varid, patm, ntmapout)
         endif
      endif
      !
      if (snapwave) then
         call write_cell_var(map_file%ncid, map_file%hm0_varid,   snapwave_H,    ntmapout, use_sw_index=.true., scale=sq2)
         call write_cell_var(map_file%ncid, map_file%hm0ig_varid, snapwave_H_ig, ntmapout, use_sw_index=.true., scale=sq2)
         !
         if (store_wave_forces) then
            call write_cell_var(map_file%ncid, map_file%fwx_varid,           fwx,            ntmapout)
            call write_cell_var(map_file%ncid, map_file%fwy_varid,           fwy,            ntmapout)
            call write_cell_var(map_file%ncid, map_file%tp_varid,            snapwave_Tp,    ntmapout, use_sw_index=.true.)
            call write_cell_var(map_file%ncid, map_file%tpig_varid,          snapwave_Tp_ig, ntmapout, use_sw_index=.true.)
            call write_cell_var(map_file%ncid, map_file%beta_varid,          betamean,       ntmapout)
            call write_cell_var(map_file%ncid, map_file%snapwavedepth_varid, snapwave_depth, ntmapout, use_sw_index=.true.)
         endif
         !
         if (store_wave_direction) then
            call write_cell_var(map_file%ncid, map_file%wavdir_varid, snapwave_mean_direction, ntmapout, use_sw_index=.true.)
         endif
         !
         if (wavemaker) then
            call write_cell_var(map_file%ncid, map_file%zsm_varid, zsm, ntmapout)
         endif
      endif
      !
      if (nonhydrostatic) then
         !
         vtmp = FILL_VALUE     
         !
         do nmq = 1, quadtree_nr_points
            !
            nm = index_sfincs_in_quadtree(nmq)
            !
            if (nm>0) then            
               !
               ! Look up pressure in 'limited' nonh array that only has values where nonh mask is set to active
               !
               if (row_index_of_nm(nm) > 0) then
                  !
                  vtmp(nmq) = pnh(row_index_of_nm(nm))
                  !
               endif   
               !
            endif   
         enddo
         !
         NF90(nf90_put_var(map_file%ncid, map_file%pnonh_varid, vtmp, (/1, ntmapout/)))
         !
      endif

   else
   !
   allocate(zsg(mmax, nmax))    ! still used by velocity & nonhydro blocks below
   !
   ! zs: wet-cell filter
   if (subgrid) then
      call write_cell_var_wet(map_file%ncid, map_file%zs_varid, real(zs,4), subgrid_z_zmin, ntmapout, output_delta=.false.)
   else
      call write_cell_var_wet(map_file%ncid, map_file%zs_varid, real(zs,4), zb, ntmapout, output_delta=.false.)
   endif
   !
   if (store_dynamic_bed_level .and. .not. subgrid) then
      call write_cell_var(map_file%ncid, map_file%zb_varid, zb, ntmapout)
   endif
   !
   ! h: delta (zs - zref), no wet filter (preserves legacy regular-grid behavior)
   if (subgrid .eqv. .false. .or. store_hsubgrid .eqv. .true.) then
      if (subgrid) then
         call write_cell_var_wet(map_file%ncid, map_file%h_varid, real(zs,4), subgrid_z_zmin, ntmapout, output_delta=.true., check_wet=.false.)
      else
         call write_cell_var_wet(map_file%ncid, map_file%h_varid, real(zs,4), zb, ntmapout, output_delta=.true., check_wet=.false.)
      endif
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
         nmd1 = z_index_uv_md(nm)
         nmu1 = z_index_uv_mu(nm)
         uz = 0.0
         if (nmd1>0) then
            uz = uz + 0.5*uv(nmd1)
         endif   
         if (nmu1>0) then
            uz = uz + 0.5*uv(nmu1)
         endif   
         !
         ndm1 = z_index_uv_nd(nm)
         num1 = z_index_uv_nu(nm)
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
   if (subgrid) then
      if (store_zvolume) then
         call write_cell_var(map_file%ncid, map_file%zvolume_varid, real(z_volume,4), ntmapout)
      endif
      if (store_storagevolume) then
         call write_cell_var(map_file%ncid, map_file%storagevolume_varid, storage_volume, ntmapout)
      endif
   endif
   !
   if (inftype == 'cnb') then
      call write_cell_var(map_file%ncid, map_file%Seff_varid, scs_Se,   ntmapout)
   elseif (inftype == 'gai') then
      call write_cell_var(map_file%ncid, map_file%Seff_varid, GA_sigma, ntmapout)
   elseif (inftype == 'hor') then
      call write_cell_var(map_file%ncid, map_file%Seff_varid, qinfmap,  ntmapout)
   endif
   !
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
   if (snapwave) then
      call write_cell_var(map_file%ncid, map_file%hm0_varid,   hm0,    ntmapout)
      call write_cell_var(map_file%ncid, map_file%hm0ig_varid, hm0_ig, ntmapout)
      if (store_wave_forces) then
         call write_cell_var(map_file%ncid, map_file%fwx_varid,           fwx,            ntmapout)
         call write_cell_var(map_file%ncid, map_file%fwy_varid,           fwy,            ntmapout)
         call write_cell_var(map_file%ncid, map_file%tp_varid,            sw_tp,          ntmapout)
         call write_cell_var(map_file%ncid, map_file%tpig_varid,          sw_tp_ig,       ntmapout)
         call write_cell_var(map_file%ncid, map_file%beta_varid,          betamean,       ntmapout)
         call write_cell_var(map_file%ncid, map_file%snapwavedepth_varid, snapwave_depth, ntmapout)
      endif
      if (wavemaker) then
         call write_cell_var(map_file%ncid, map_file%zsm_varid, zsm, ntmapout)
      endif
   endif
   !
   if (nonhydrostatic) then
      !
      zsg = FILL_VALUE       
      !
      do nm = 1, np
         !
         n    = z_index_z_n(nm)
         m    = z_index_z_m(nm)
         !
         ! Look up pressure in 'limited' nonh array that only has values where nonh mask is set to active
         !
         if (row_index_of_nm(nm) > 0) then
            !
            zsg(m, n) = pnh(row_index_of_nm(nm))
            !
         endif   
         !
      enddo
      !
      NF90(nf90_put_var(map_file%ncid, map_file%pnonh_varid, zsg, (/1, 1, ntmapout/))) ! write h
      !
   endif

   endif
   !
   NF90(nf90_sync(map_file%ncid)) !write away intermediate data
   !
   end subroutine ncoutput_update_map

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   subroutine ncoutput_update_his(t,nthisout)
   ! Write time, zs, u, v, prcp of points 
   !
   use sfincs_data   
   use sfincs_crosssections
   use sfincs_runup_gauges
   use sfincs_snapwave
   !
   implicit none   
   !
   integer :: iobs, nm, idrn
   !
   integer :: nthisout      
   integer :: nmd1, nmu1, ndm1, num1
   !
   real*4                  :: uz, vz
   real*8                  :: t
!   real*4, dimension(nobs) :: zobs, hobs
   real*4, dimension(nobs) :: uobs
   real*4, dimension(nobs) :: vobs   
   real*4, dimension(nobs) :: uvmag   
   real*4, dimension(nobs) :: uvdir   
   real*4, dimension(nobs) :: tprcp
   real*4, dimension(nobs) :: tcumprcp
   real*4, dimension(nobs) :: tqinf
   real*4, dimension(nobs) :: tS_effective
   real*4, dimension(nobs) :: tpatm
   real*4, dimension(nobs) :: twndmag
   real*4, dimension(nobs) :: twnddir
   real*4, dimension(nobs) :: hm0obs
   real*4, dimension(nobs) :: hm0igobs
   real*4, dimension(nobs) :: zsmobs
   real*4, dimension(nobs) :: tpobs
   real*4, dimension(nobs) :: tpigobs   
   real*4, dimension(nobs) :: wavdirobs
   real*4, dimension(nobs) :: dirsprobs
   real*4, dimension(ndrn) :: q_drain   
   real*4, dimension(nobs) :: dwobs
   real*4, dimension(nobs) :: dfobs
   real*4, dimension(nobs) :: dwigobs
   real*4, dimension(nobs) :: dfigobs
   real*4, dimension(nobs) :: cgobs
   real*4, dimension(nobs) :: betaobs
   real*4, dimension(nobs) :: srcigobs
   real*4, dimension(nobs) :: alphaigobs
   real*4, dimension(:), allocatable :: qq, zz
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
   tcumprcp     = FILL_VALUE
   tqinf        = FILL_VALUE
   tS_effective = FILL_VALUE
   tpatm        = FILL_VALUE
   twndmag      = FILL_VALUE
   twnddir      = FILL_VALUE
   q_drain      = FILL_VALUE
   dwobs        = FILL_VALUE
   dfobs        = FILL_VALUE
   cgobs        = FILL_VALUE
   betaobs      = FILL_VALUE
   srcigobs     = FILL_VALUE
   alphaigobs   = FILL_VALUE   
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
            ! Regular point with four surrounding cells of the same size
            !
            nmd1 = z_index_uv_md(nm)
            nmu1 = z_index_uv_mu(nm)
            ndm1 = z_index_uv_nd(nm)
            num1 = z_index_uv_mu(nm)
            uz  = 0.5 * (uv(nmd1) + uv(nmu1))
            vz  = 0.5 * (uv(ndm1) + uv(num1))
            !
            uobs(iobs)  = cosrot * uz - sinrot * vz                         
            vobs(iobs)  = sinrot * uz + cosrot * vz
            uvmag(iobs) = sqrt(uobs(iobs)**2 + vobs(iobs)**2)
            uvdir(iobs) = atan2(vobs(iobs), uobs(iobs)) * 180 / pi
            !
         endif
         !
         if (infiltration) then
            !
            tqinf(iobs) = qinfmap(nm) * 3600000 ! show as mm/hr
            ! 
            ! Output for CN and GA method
            !
            if (inftype == 'cnb') then
               !
               tS_effective(iobs) = scs_Se(nm)
               !
            elseif (inftype == 'gai') then
               !
               tS_effective(iobs) = GA_sigma(nm)
               !
            endif
            !
         endif  
         !
         if (store_meteo) then
            !
            if (wind) then
               !
               twndmag(iobs) = sqrt(windu(nm)**2 + windv(nm)**2)
               twnddir(iobs) = 270.0 - atan2(windv(nm), windu(nm)) * 180 / pi
               if (twnddir(iobs) < 0.0) twnddir(iobs) = twnddir(iobs) + 360.0
               if (twnddir(iobs) > 360.0) twnddir(iobs) = twnddir(iobs) - 360.0
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
               tprcp(iobs) = prcp(nm) * 3600000 ! show as mm/hr
               !
               if (store_cumulative_precipitation) then
                  !
                  tcumprcp(iobs) = cumprcp(nm) ! show as m
                  !
               endif   
               !
            endif
            !         
         endif
         !
         if (snapwave) then
            !
            hm0obs(iobs)   = hm0(nm)
            hm0igobs(iobs) = hm0_ig(nm)
            tpobs(iobs)    = sw_tp(nm)
            tpigobs(iobs)  = sw_tp_ig(nm)            
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
            if (store_wave_forces) then
               !
               dwobs(iobs)    = dw(nm)
               dfobs(iobs)    = df(nm)
               dwigobs(iobs)  = dwig(nm)
               dfigobs(iobs)  = dfig(nm)
               cgobs(iobs)    = cg(nm) 
               betaobs(iobs)  = betamean(nm)               
               srcigobs(iobs) = srcig(nm)               
               alphaigobs(iobs) = alphaig(nm)                              
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
      elseif (inftype == 'gai') then
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
      NF90(nf90_put_var(his_file%ncid, his_file%tpig_varid, tpigobs, (/1, nthisout/)))      
      !
      if (store_wave_direction) then
         !
         NF90(nf90_put_var(his_file%ncid, his_file%wavdir_varid, wavdirobs, (/1, nthisout/)))
         !NF90(nf90_put_var(his_file%ncid, his_file%dirspr_varid, dirsprobs, (/1, nthisout/)))
         !
      endif
      !            
      if (wavemaker) then
         !
         NF90(nf90_put_var(his_file%ncid, his_file%zsm_varid, zsmobs, (/1, nthisout/)))
         !
      endif
      !
      if (store_wave_forces) then
         !
         NF90(nf90_put_var(his_file%ncid, his_file%dw_varid, dwobs, (/1, nthisout/)))
         NF90(nf90_put_var(his_file%ncid, his_file%df_varid, dfobs, (/1, nthisout/)))        
         ! 
         NF90(nf90_put_var(his_file%ncid, his_file%dwig_varid, dwigobs, (/1, nthisout/)))
         NF90(nf90_put_var(his_file%ncid, his_file%dfig_varid, dfigobs, (/1, nthisout/)))        
         !
         NF90(nf90_put_var(his_file%ncid, his_file%cg_varid, cgobs, (/1, nthisout/)))
         !
         NF90(nf90_put_var(his_file%ncid, his_file%beta_varid, betaobs, (/1, nthisout/)))
         NF90(nf90_put_var(his_file%ncid, his_file%srcig_varid, srcigobs, (/1, nthisout/)))                  
         NF90(nf90_put_var(his_file%ncid, his_file%alphaig_varid, alphaigobs, (/1, nthisout/)))         
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
         NF90(nf90_put_var(his_file%ncid, his_file%patm_varid, tpatm, (/1, nthisout/))) ! write patmos
         !
      endif   
      !
      if (precip) then
         !
         NF90(nf90_put_var(his_file%ncid, his_file%prcp_varid, tprcp, (/1, nthisout/))) ! write prcp
         !
         if (store_cumulative_precipitation) then
            !
            NF90(nf90_put_var(his_file%ncid, his_file%cumprcp_varid, tcumprcp, (/1, nthisout/))) ! write cumulative prcp
            !
         endif   
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
      NF90(nf90_put_var(his_file%ncid, his_file%discharge_varid, qq, (/1, nthisout/))) ! write discharge
      !
   endif
   !
   if (nr_runup_gauges>0) then
      !
      ! Get run-up elevations
      !
      call get_runup_levels(zz)
      !
      NF90(nf90_put_var(his_file%ncid, his_file%runup_gauge_zs_varid, zz, (/1, nthisout/))) ! write run up level
      !
   endif
   !
   if (ndrn>0) then
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
      NF90(nf90_put_var(his_file%ncid, his_file%drain_varid, q_drain, (/1, nthisout/))) ! write discharge of sink point
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

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   subroutine ncoutput_update_max(t,ntmaxout)
   !
   ! Write maximum values to map file (handles both regular and quadtree grids)
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   real*8                              :: t
   integer                             :: ntmaxout, nm, n, m
   real*4, dimension(:,:), allocatable :: zstmp
   real*4, dimension(:),   allocatable :: hmax_out, hmean
   !
   ! Scalar time of this max-record
   NF90(nf90_put_var(map_file%ncid, map_file%timemax_varid, t, (/ntmaxout/)))
   !
   ! Maximum water level
   if (subgrid) then
      call write_cell_var_wet(map_file%ncid, map_file%zsmax_varid, zsmax, subgrid_z_zmin, ntmaxout, output_delta=.false.)
   else
      call write_cell_var_wet(map_file%ncid, map_file%zsmax_varid, zsmax, zb,             ntmaxout, output_delta=.false.)
   endif
   !
   ! Maximum water depth (optional, supports subgrid mean-depth)
   if (.not. subgrid .or. store_hsubgrid) then
      allocate(hmax_out(np))
      hmax_out = FILL_VALUE
      if (store_hmean .and. subgrid) then
         allocate(hmean(np))
         call compute_subgrid_mean_depth(zsmax, hmean)
      endif
      do nm = 1, np
         if (subgrid) then
            if ( (zsmax(nm) - subgrid_z_zmin(nm)) > huthresh) then
               if (store_hmean) then
                  hmax_out(nm) = hmean(nm)
               else
                  hmax_out(nm) = zsmax(nm) - subgrid_z_zmin(nm)
               endif
            endif
         else
            if ( (zsmax(nm) - zb(nm)) > huthresh) then
               hmax_out(nm) = zsmax(nm) - zb(nm)
            endif
         endif
      enddo
      call write_cell_var(map_file%ncid, map_file%hmax_varid, hmax_out, ntmaxout, check_kcs=.true.)
      deallocate(hmax_out)
      if (allocated(hmean)) deallocate(hmean)
   endif
   !
   ! Cumulative rainfall and infiltration
   if (store_cumulative_precipitation) then
      call write_cell_var(map_file%ncid, map_file%cumprcp_varid, cumprcp, ntmaxout, check_kcs=.true.)
      call write_cell_var(map_file%ncid, map_file%cuminf_varid,  cuminf,  ntmaxout, check_kcs=.true.)
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
   ! When zsmax occurred — quadtree uses helper; regular grid preserves legacy >0 filter
   if (store_t_zsmax) then
      if (use_quadtree) then
         call write_cell_var(map_file%ncid, map_file%t_zsmax_varid, t_zsmax, ntmaxout, check_kcs=.true.)
      else
         allocate(zstmp(mmax, nmax))
         zstmp = FILL_VALUE
         do nm = 1, np
            n = z_index_z_n(nm)
            m = z_index_z_m(nm)
            if (t_zsmax(nm) > 0) zstmp(m, n) = t_zsmax(nm)
         enddo
         NF90(nf90_put_var(map_file%ncid, map_file%t_zsmax_varid, zstmp, (/1, 1, ntmaxout/)))
         deallocate(zstmp)
      endif
   endif
   !
   ! Maximum wind speed
   if (wind .and. store_wind_max .and. meteo3d) then
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
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   subroutine ncoutput_write_timestep_analysis()
   !
   ! Write timestep_analysis_average_required_timestep and timestep_analysis_percentage_limiting_per_cell once at end of simulation (no time dimension)
   !
   use sfincs_data
   use quadtree
   !
   implicit none
   !
   real*4, dimension(:),   allocatable :: vtmp1, vtmp2, vtmp3
   real*4, dimension(:,:), allocatable :: zsg1, zsg2
   integer :: nm, nmq, n, m
   !
   if (use_quadtree) then
      !
      allocate(vtmp1(quadtree_nr_points))
      allocate(vtmp2(quadtree_nr_points))
      !
      vtmp1 = FILL_VALUE
      vtmp2 = FILL_VALUE
      !
      do nm = 1, np
         !
         nmq = index_quadtree_in_sfincs(nm)
         !
         if (timestep_analysis_average_required_timestep_per_cell(nm) > 0.0) then
            !
            vtmp1(nmq) = timestep_analysis_average_required_timestep_per_cell(nm)
            !
         endif
         !
         vtmp2(nmq) = timestep_analysis_percentage_limiting_per_cell(nm)
         ! 
      enddo
      !
      NF90(nf90_put_var(map_file%ncid, map_file%average_required_timestep_varid, vtmp1))
      NF90(nf90_put_var(map_file%ncid, map_file%percentage_limiting_varid, vtmp2))
      !
      deallocate(vtmp1)
      deallocate(vtmp2)
      !
   else
      !
      allocate(zsg1(mmax, nmax))
      allocate(zsg2(mmax, nmax))
      !
      zsg1 = FILL_VALUE
      zsg2 = FILL_VALUE
      !
      do nm = 1, np
         !
         n = z_index_z_n(nm)
         m = z_index_z_m(nm)
         !
         if (timestep_analysis_average_required_timestep_per_cell(nm) > 0.0) then
            !
            zsg1(m, n) = timestep_analysis_average_required_timestep_per_cell(nm)
            !
         endif
         !
         zsg2(m, n) = timestep_analysis_percentage_limiting_per_cell(nm)
         ! 
      enddo
      !
      NF90(nf90_put_var(map_file%ncid, map_file%average_required_timestep_varid, zsg1))
      NF90(nf90_put_var(map_file%ncid, map_file%percentage_limiting_varid, zsg2))
      !
      deallocate(zsg1)
      deallocate(zsg2)
      !
   endif   
   !
   end subroutine ncoutput_write_timestep_analysis
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   subroutine ncoutput_write_tsunami_arrival_time() 
   ! Add tsunami_arrival_time
   use sfincs_data
   !   
   implicit none  
   !
   real*4, dimension(:,:), allocatable :: zsg
   real*4,    dimension(:),   allocatable :: vtmp   
   integer    :: nm, nmq, n, m, n_faces, quadtree_nr_points   
   !
   if (use_quadtree) then
       !
       allocate(vtmp(n_faces))          
       vtmp = FILL_VALUE   
       !
       do nmq = 1, quadtree_nr_points
           nm = index_sfincs_in_quadtree(nmq)
           if (nm>0) then
           vtmp(nmq) = tsunami_arrival_time(nm)
           endif
       enddo
       !
       ! write tsunami_arrival_time_varid              
       NF90(nf90_put_var(map_file%ncid, map_file%tsunami_arrival_time_varid, zsg))         
       !
   else
       ! regular grid
       allocate(zsg(mmax, nmax))   
       zsg = FILL_VALUE
       !   
       do nm = 1, np
           !
           n    = z_index_z_n(nm)
           m    = z_index_z_m(nm)
           !       
           zsg(m, n) = tsunami_arrival_time(nm)
           !
       enddo
       !
       ! write tsunami_arrival_time_varid       
       NF90(nf90_put_var(map_file%ncid, map_file%tsunami_arrival_time_varid, zsg, (/1, 1/)))     
       !
   endif     
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
   if (nobs==0 .and. nrcrosssections==0 .and. nrstructures==0 .and. nrthindams==0 .and. ndrn==0) then ! If no observation points, cross-sections, structures 9weir or thin dam), or drains; hisfile        
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
   !
   ! Because of overlapping names, only important specific values from snapwave_data
   use snapwave_data, only: gamma, gammax, alpha, hmin, fw0, fw0_ig, dt, tol, dtheta, crit, nr_sweeps, baldock_opt, baldock_ratio, &
       igwaves_opt, alpha_ig, gamma_ig, shinc2ig, alphaigfac, baldock_ratio_ig, ig_opt, herbers_opt, tpig_opt, eeinc2ig, tinc2ig, &
       snapwave_jonswapfile, snapwave_encfile, snapwave_bndfile, snapwave_bhsfile, snapwave_btpfile, snapwave_bwdfile, snapwave_bdsfile, upwfile, gridfile
   
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
        NF90(nf90_put_att(ncid, varid, 'rhoa',rhoa))        
        NF90(nf90_put_att(ncid, varid, 'rhow',rhow))        
        NF90(nf90_put_att(ncid, varid, 'inputformat',inputtype))        
        NF90(nf90_put_att(ncid, varid, 'outputformat',outputtype))   
        NF90(nf90_put_att(ncid, varid, 'outputtype_map',outputtype_map))
        NF90(nf90_put_att(ncid, varid, 'outputtype_his',outputtype_his))
        NF90(nf90_put_att(ncid, varid, 'bndtype',bndtype))
        NF90(nf90_put_att(ncid, varid, 'advection',logical2int(advection)))  
        NF90(nf90_put_att(ncid, varid, 'wavemaker_nfreqs_ig', wavemaker_nfreqs_ig))  
        NF90(nf90_put_att(ncid, varid, 'wavemaker_freqmin_ig',wavemaker_freqmin_ig))  
        NF90(nf90_put_att(ncid, varid, 'wavemaker_freqmax_ig',wavemaker_freqmax_ig))  
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
        NF90(nf90_put_att(ncid, varid, 'qinffile',qinffile))   
        NF90(nf90_put_att(ncid, varid, 'scsfile',scsfile)) 
        NF90(nf90_put_att(ncid, varid, 'smaxfile',smaxfile)) 
        NF90(nf90_put_att(ncid, varid, 'sefffile',sefffile)) 
        NF90(nf90_put_att(ncid, varid, 'ksfile',ksfile)) 
        NF90(nf90_put_att(ncid, varid, 'psifile',psifile)) 
        NF90(nf90_put_att(ncid, varid, 'sigmafile',sigmafile)) 
        NF90(nf90_put_att(ncid, varid, 'z0lfile',z0lfile)) 
        NF90(nf90_put_att(ncid, varid, 'wavemaker_wvmfile',wavemaker_wvmfile)) 
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
        NF90(nf90_put_att(ncid, varid, 'storevelmax',storevelmax)) 
        NF90(nf90_put_att(ncid, varid, 'storefluxmax',storefluxmax))        
        NF90(nf90_put_att(ncid, varid, 'storevel',storevel)) 
        NF90(nf90_put_att(ncid, varid, 'storecumprcp',storecumprcp)) 
        NF90(nf90_put_att(ncid, varid, 'storetwet',storetwet)) 
        NF90(nf90_put_att(ncid, varid, 'storehsubgrid',storehsubgrid)) 
        NF90(nf90_put_att(ncid, varid, 'twet_threshold',twet_threshold)) 
        NF90(nf90_put_att(ncid, varid, 'store_tsunami_arrival_time',logical2int(store_tsunami_arrival_time))) 
        NF90(nf90_put_att(ncid, varid, 'tsunami_arrival_threshold',tsunami_arrival_threshold)) 
        NF90(nf90_put_att(ncid, varid, 'storeqdrain',storeqdrain)) 
        NF90(nf90_put_att(ncid, varid, 'storezvolume',storezvolume)) 
        NF90(nf90_put_att(ncid, varid, 'writeruntime',wrttimeoutput)) 
        NF90(nf90_put_att(ncid, varid, 'debug',logical2int(debug))) 
        NF90(nf90_put_att(ncid, varid, 'storemeteo',storemeteo)) 
        NF90(nf90_put_att(ncid, varid, 'storemaxwind',logical2int(store_wind_max))) 
        NF90(nf90_put_att(ncid, varid, 'storefw',logical2int(store_wave_forces)))         
        NF90(nf90_put_att(ncid, varid, 'storewavdir', logical2int(store_wave_direction))) 
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
        NF90(nf90_put_att(ncid, varid, 'snapwave_baldock_opt',baldock_opt)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_baldock_ratio',baldock_ratio)) 
        !
        ! SnapWave IG
        !
        NF90(nf90_put_att(ncid, varid, 'snapwave_igwaves',igwaves_opt))         
        NF90(nf90_put_att(ncid, varid, 'snapwave_alpha_ig',alpha_ig)) 
        NF90(nf90_put_att(ncid, varid, 'snapwave_gammaig',gamma_ig)) 
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
   !
   !
   subroutine compute_subgrid_mean_depth(z, hmean)
   !
   ! This subroutine cannot sit in sfincs_subgrid.f90 because that uses the same netcdf module
   !
   use sfincs_data
   !
   implicit none
   !
   real*4, intent(in)  :: z(np) ! max water level
   real*4, intent(out) :: hmean(np) 
   !
   integer    :: nm, m, n, ivol, ilevel, ip
   real*8     :: volume
   real*4     :: dzvol
   real*4     :: facint
   real*4     :: one_minus_facint 
   !
   ! Compute volumes and mean depths
   !
   do nm = 1, np
      !
      if (z(nm) >= subgrid_z_zmax(nm)) then
         !
         ! Entire cell is wet, no interpolation from table needed
         !
         if (crsgeo) then
            volume = subgrid_z_volmax(nm) + cell_area_m2(nm) * (z(nm) - max(subgrid_z_zmax(nm), -20.0))
         else   
            volume = subgrid_z_volmax(nm) + cell_area(z_flags_iref(nm)) * (z(nm) - max(subgrid_z_zmax(nm), -20.0))
         endif
         !
      else   
         !
         ! Interpolation required
         !
         ivol = 1
         do ilevel = 2, subgrid_nlevels
            if (subgrid_z_dep(ilevel, nm) > z(nm)) then
               ivol = ilevel - 1
               exit
            endif
         enddo
         !
         dzvol  = subgrid_z_volmax(nm) / (subgrid_nlevels - 1)
         facint = (z(nm) - subgrid_z_dep(ivol, nm)) / max(subgrid_z_dep(ivol + 1, nm) - subgrid_z_dep(ivol, nm), 0.001)
         volume = (ivol - 1) * dzvol + facint * dzvol
         !
      endif
      !
      ! Compute mean depth in cell
      !
      if (crsgeo) then
         !
         hmean(nm) = volume / cell_area_m2(nm)
         !
      else   
         !
         hmean(nm) = volume / cell_area(z_flags_iref(nm))
         !
      endif
      !
   enddo
   !
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
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !   
   function logical2int(lgc) result (i)
   !
   implicit none
   !
   logical              :: lgc
   integer              :: i
   !
   if (lgc) then
      i = 1
   else
      i = 0
   endif
   !
   end function   
   !
   end module
