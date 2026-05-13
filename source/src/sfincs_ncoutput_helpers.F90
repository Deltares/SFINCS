#define NF90(nf90call) call handle_err(nf90call,__FILE__,__LINE__)
!
! ============================================================================
! sfincs_ncoutput_helpers — shared types, state and helper subroutines for
! the SFINCS NetCDF output layer.
!
! This module owns:
!   - map_type / his_type derived types and their instances (map_file / his_file)
!   - Fill-value constants FILL_VALUE / FILL_VALUE_INT
!   - Dimension-context module variables set once by ncoutput_map_init before
!     any def_* calls: nsd, dims_s, dims_st, dims_sm, coord_str
!   - pt_coord parameter used by def_time_point_float / def_his_point_coord
!   - All helper subroutines (cell writers, def wrappers, generic ncdef helpers,
!     handle_err)
!
! sfincs_ncoutput uses this module and supplies the public init/update/finalize
! subroutines.
! ============================================================================
module sfincs_ncoutput_helpers
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
      integer :: zs_varid, zsmax_varid, h_varid, u_varid, v_varid, tmax_varid, infstate_varid, t_zsmax_varid
      integer :: zvolume_varid, storagevolume_varid
      integer :: hmax_varid, vmax_varid, qmax_varid, cumprcp_varid, cuminf_varid, windmax_varid
      integer :: patm_varid, wind_u_varid, wind_v_varid, precip_varid
      integer :: hm0_varid, hm0ig_varid, snapwavemsk_varid, tp_varid, tpig_varid, wavdir_varid
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
      integer :: drain_varid
      integer :: zb_varid
      integer :: time_varid
      integer :: zs_varid, h_varid, u_varid, v_varid, prcp_varid, cumprcp_varid, discharge_varid, uvmag_varid, uvdir_varid
      integer :: patm_varid, wind_speed_varid, wind_dir_varid
      integer :: inp_varid, total_runtime_varid, average_dt_varid, status_varid
      integer :: hm0_varid, hm0ig_varid, zsm_varid, tp_varid, tpig_varid, wavdir_varid
      integer :: dw_varid, df_varid, dwig_varid, dfig_varid, cg_varid, beta_varid, srcig_varid, alphaig_varid
      integer :: runup_gauge_name_varid, runup_gauge_zs_varid
      !
   end type
   !
   type(map_type) :: map_file
   type(his_type) :: his_file
   !
   real*4,  parameter :: FILL_VALUE     = -99999.0
   integer, parameter :: FILL_VALUE_INT = -999
   !
   ! pt_coord: coordinate attribute used by his-file point-variable def helpers.
   character(*), parameter :: pt_coord = 'station_id station_name point_x point_y'
   !
   ! Dimension context written once by ncoutput_map_init (before any def_* calls)
   ! and read by the def_* wrappers below.
   integer      :: nsd
   integer      :: dims_s(3), dims_st(3), dims_sm(3)
   character(4) :: coord_str

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

   subroutine write_cell_var(ncid, varid, source, nt, use_sw_index, check_kcs, scale, min_value)
      !
      ! Scatter SFINCS cell-centered source(:) into the file-shaped buffer
      ! and put it to the NetCDF variable at timestep nt.
      !   use_sw_index=.true. -> quadtree gather uses index_sw_in_qt (snapwave)
      !   check_kcs=.true.    -> quadtree skips cells with kcs<=0 (max-file style)
      !   scale               -> multiply source by scale on write (unit conversions)
      !   min_value           -> only write cells where source(nm) > min_value
      !                          (used for t_zsmax to mask never-reached cells)
      !
      use sfincs_data
      use sfincs_snapwave, only: index_sw_in_qt
      use quadtree
      !
      integer, intent(in)           :: ncid, varid, nt
      real*4,  intent(in)           :: source(:)
      logical, intent(in), optional :: use_sw_index, check_kcs
      real*4,  intent(in), optional :: scale, min_value
      !
      real*4,  allocatable :: buf_q(:), buf_r(:,:)
      logical :: use_sw, filt_kcs, filt_min
      real*4  :: fac, vmin
      integer :: nm, nmq
      !
      use_sw   = .false.; if (present(use_sw_index)) use_sw   = use_sw_index
      filt_kcs = .false.; if (present(check_kcs))    filt_kcs = check_kcs
      filt_min = .false.; if (present(min_value))    filt_min = .true.
      fac      = 1.0;     if (present(scale))        fac      = scale
      vmin     = 0.0;     if (present(min_value))    vmin     = min_value
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
            if (nm <= 0) cycle
            if (filt_kcs .and. kcs(nm) <= 0) cycle
            if (filt_min .and. source(nm) <= vmin) cycle
            buf_q(nmq) = fac * source(nm)
         enddo
         NF90(nf90_put_var(ncid, varid, buf_q, (/1, nt/)))
      else
         allocate(buf_r(mmax, nmax))
         buf_r = FILL_VALUE
         do nm = 1, np
            if (filt_kcs .and. kcs(nm) <= 0) cycle
            if (filt_min .and. source(nm) <= vmin) cycle
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

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Static-cell write helpers
   !
   ! Gather a SFINCS-indexed source array into the right per-grid buffer
   ! (1D for quadtree, 2D for regular) and write it via nf90_put_var.
   ! Hides the use_quadtree branch, the index-map choice and the buffer
   ! allocation. Used by ncoutput_map_init (static map fields),
   ! ncoutput_write_tsunami_arrival_time, and ncoutput_write_timestep_analysis.
   !
   !   put_static_cell_float : real*4 source, optional scale,
   !                           optional sw_index for snapwave-indexed sources
   !   put_static_cell_mask  : real*4 source written to NF90_INT (quadtree)
   !                           or NF90_FLOAT (regular); cast int*1 sources
   !                           via real(.,4) at the call site
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine put_static_cell_float(ncid, varid, source, fill, scale, sw_index, min_value)
      use sfincs_data
      use sfincs_snapwave, only: index_sw_in_qt
      use quadtree
      !
      integer, intent(in)           :: ncid, varid
      real*4,  intent(in)           :: source(:)
      real*4,  intent(in)           :: fill
      real*4,  intent(in), optional :: scale     ! default: 1.0
      logical, intent(in), optional :: sw_index  ! default: .false. (use sfincs index)
      real*4,  intent(in), optional :: min_value ! only write cells where source(nm) > min_value
      !
      real*4, allocatable :: buf_q(:), buf_r(:,:)
      real*4  :: sc, vmin
      logical :: usesw, filt_min
      integer :: nmq, nm
      !
      sc = 1.0;       if (present(scale))     sc       = scale
      usesw = .false.; if (present(sw_index)) usesw    = sw_index
      filt_min = .false.; if (present(min_value)) filt_min = .true.
      vmin = 0.0;     if (present(min_value)) vmin     = min_value
      !
      if (use_quadtree) then
         allocate(buf_q(quadtree_nr_points))
         buf_q = fill
         do nmq = 1, quadtree_nr_points
            if (usesw) then
               nm = index_sw_in_qt(nmq)
            else
               nm = index_sfincs_in_quadtree(nmq)
            endif
            if (nm <= 0) cycle
            if (filt_min .and. source(nm) <= vmin) cycle
            buf_q(nmq) = source(nm) * sc
         enddo
         NF90(nf90_put_var(ncid, varid, buf_q))
      else
         allocate(buf_r(mmax, nmax))
         buf_r = fill
         do nm = 1, np
            if (filt_min .and. source(nm) <= vmin) cycle
            buf_r(z_index_z_m(nm), z_index_z_n(nm)) = source(nm) * sc
         enddo
         NF90(nf90_put_var(ncid, varid, buf_r, (/1, 1/)))
      endif
   end subroutine put_static_cell_float

   subroutine put_static_cell_mask(ncid, varid, source, sw_index)
      ! Mask-style write: NF90_INT on both grid types — 1D `face` on quadtree,
      ! 2D `(m, n)` on regular. Source is real*4 — callers pass real(kcs, 4)
      ! for int*1 masks (kcs) or snapwave_mask directly (already real*4).
      use sfincs_data
      use sfincs_snapwave, only: index_sw_in_qt
      use quadtree
      !
      integer, intent(in)           :: ncid, varid
      real*4,  intent(in)           :: source(:)
      logical, intent(in), optional :: sw_index
      !
      integer*4, allocatable :: buf_qi(:)
      integer*4, allocatable :: buf_ri(:,:)
      logical :: usesw
      integer :: nmq, nm
      !
      usesw = .false.
      if (present(sw_index)) usesw = sw_index
      !
      ! Source is real*4 — int() truncates toward zero. Used only for mask
      ! vars (kcs, snapwave_mask) whose values are small non-negative
      ! integers, so truncation is exact.
      if (use_quadtree) then
         allocate(buf_qi(quadtree_nr_points))
         buf_qi = 0
         do nmq = 1, quadtree_nr_points
            if (usesw) then
               nm = index_sw_in_qt(nmq)
            else
               nm = index_sfincs_in_quadtree(nmq)
            endif
            if (nm > 0) buf_qi(nmq) = int(source(nm), 4)
         enddo
         NF90(nf90_put_var(ncid, varid, buf_qi))
      else
         allocate(buf_ri(mmax, nmax))
         ! Initialise to FILL_VALUE_INT so cells outside np (the active SFINCS
         ! domain) read as fill, not as a valid mask=0. Cells with kcs(nm)==0
         ! inside np are still written as 0 (a valid inactive-but-present marker).
         buf_ri = FILL_VALUE_INT
         do nm = 1, np
            buf_ri(z_index_z_m(nm), z_index_z_n(nm)) = int(source(nm), 4)
         enddo
         NF90(nf90_put_var(ncid, varid, buf_ri, (/1, 1/)))
      endif
   end subroutine put_static_cell_mask

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Map-init def wrappers
   !
   ! These read the dimension-context module variables (nsd, dims_*, coord_str)
   ! set at the top of ncoutput_map_init, so call sites carry only var-specific
   ! info (name, varid, units, long_name, optional attrs).
   !
   ! def_static_cell_float   : dims_s(1:nsd)    + coord_str
   ! def_static_cell_int     : dims_s(1:nsd)    (no coord)
   ! def_time_cell_float     : dims_st(1:nsd+1) + coord_str
   ! def_maxtime_cell_float  : dims_sm(1:nsd+1) + coord_str
   ! add_ugrid_face_attrs    : mesh/location on quadtree cell vars
   ! def_mesh2d_node_coord   : UGRID node coord (float geo / double projected)
   ! def_grid_axis_coord     : SGRID face/corner coord (always float)
   ! put_2d                  : nf90_put_var with (/1, 1/) start
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine def_static_cell_float(name, varid, units, long_name, &
                                     standard_name, cell_methods, description)
      use sfincs_data, only: nc_deflate_level
      character(*), intent(in)           :: name
      integer,      intent(out)          :: varid
      character(*), intent(in)           :: units, long_name
      character(*), intent(in), optional :: standard_name, cell_methods, description
      call ncdef_float_var(map_file%ncid, name, dims_s(1:nsd), varid, units, long_name, &
           standard_name=standard_name, coordinates=coord_str,                          &
           cell_methods=cell_methods, description=description,                          &
           deflate_level=nc_deflate_level)
      call add_ugrid_face_attrs(varid)
   end subroutine def_static_cell_float

   subroutine def_static_cell_int(name, varid, long_name, &
                                   units, standard_name, fill_value, description, &
                                   flag_values, flag_meanings)
      use sfincs_data, only: nc_deflate_level
      character(*), intent(in)           :: name
      integer,      intent(out)          :: varid
      character(*), intent(in)           :: long_name
      character(*), intent(in), optional :: units, standard_name
      integer,      intent(in), optional :: fill_value
      character(*), intent(in), optional :: description
      integer,      intent(in), optional :: flag_values(:)
      character(*), intent(in), optional :: flag_meanings
      call ncdef_int_var(map_file%ncid, name, dims_s(1:nsd), varid, long_name, &
           fill_value=fill_value, description=description,                     &
           deflate_level=nc_deflate_level)
      if (present(units)) then
         if (len_trim(units) > 0) NF90(nf90_put_att(map_file%ncid, varid, 'units', trim(units)))
      endif
      if (present(standard_name)) then
         if (len_trim(standard_name) > 0) NF90(nf90_put_att(map_file%ncid, varid, 'standard_name', trim(standard_name)))
      endif
      ! coordinates attribute: empty on quadtree (UGRID uses mesh="mesh2d"
      ! / location="face" instead, added below), 'x y' on SGRID regular grid.
      if (len_trim(coord_str) > 0) then
         NF90(nf90_put_att(map_file%ncid, varid, 'coordinates', trim(coord_str)))
      endif
      ! CF flag-variable convention: flag_values + flag_meanings together
      ! describe the categorical values of the variable.
      if (present(flag_values)) then
         NF90(nf90_put_att(map_file%ncid, varid, 'flag_values', flag_values))
      endif
      if (present(flag_meanings)) then
         if (len_trim(flag_meanings) > 0) then
            NF90(nf90_put_att(map_file%ncid, varid, 'flag_meanings', trim(flag_meanings)))
         endif
      endif
      call add_ugrid_face_attrs(varid)
   end subroutine def_static_cell_int

   subroutine def_time_cell_float(name, varid, units, long_name, &
                                   standard_name, cell_methods, description)
      use sfincs_data, only: nc_deflate_level
      character(*), intent(in)           :: name
      integer,      intent(out)          :: varid
      character(*), intent(in)           :: units, long_name
      character(*), intent(in), optional :: standard_name, cell_methods, description
      call ncdef_float_var(map_file%ncid, name, dims_st(1:nsd+1), varid, units, long_name, &
           standard_name=standard_name, coordinates=coord_str,                             &
           cell_methods=cell_methods, description=description,                             &
           deflate_level=nc_deflate_level)
      call add_ugrid_face_attrs(varid)
   end subroutine def_time_cell_float

   subroutine def_maxtime_cell_float(name, varid, units, long_name, &
                                      standard_name, cell_methods, description)
      use sfincs_data, only: nc_deflate_level
      character(*), intent(in)           :: name
      integer,      intent(out)          :: varid
      character(*), intent(in)           :: units, long_name
      character(*), intent(in), optional :: standard_name, cell_methods, description
      call ncdef_float_var(map_file%ncid, name, dims_sm(1:nsd+1), varid, units, long_name, &
           standard_name=standard_name, coordinates=coord_str,                             &
           cell_methods=cell_methods, description=description,                             &
           deflate_level=nc_deflate_level)
      call add_ugrid_face_attrs(varid)
   end subroutine def_maxtime_cell_float

   subroutine add_ugrid_face_attrs(varid)
      ! On quadtree (UGRID) output, every cell-centred data variable must
      ! carry mesh="mesh2d" and location="face". On regular grid (SGRID) the
      ! topology var declares the grid; cell vars use coordinates='x y' (set
      ! elsewhere) and these attributes are not used.
      use sfincs_data, only: use_quadtree
      integer, intent(in) :: varid
      if (use_quadtree) then
         NF90(nf90_put_att(map_file%ncid, varid, 'mesh',     'mesh2d'))
         NF90(nf90_put_att(map_file%ncid, varid, 'location', 'face'))
      endif
   end subroutine add_ugrid_face_attrs

   subroutine def_mesh2d_node_coord(axis, varid)
      ! Define mesh2d_node_x / mesh2d_node_y for UGRID quadtree topology.
      ! Picks NF90_FLOAT+degrees (geographic) or NF90_DOUBLE+m (projected) via crsgeo.
      use sfincs_data, only: crsgeo, nc_deflate_level
      character(len=*), intent(in)  :: axis   ! 'x' or 'y'
      integer,          intent(out) :: varid
      character(len=64) :: vname
      character(len=32) :: std_name, long_name, units
      integer           :: nctype
      !
      vname = 'mesh2d_node_' // axis
      if (crsgeo) then
         nctype = NF90_FLOAT
         units  = 'degrees'
         if (axis == 'x') then
            std_name  = 'longitude'
            long_name = 'longitude'
         else
            std_name  = 'latitude'
            long_name = 'latitude'
         endif
      else
         nctype = NF90_DOUBLE
         units  = 'm'
         if (axis == 'x') then
            std_name  = 'projection_x_coordinate'
            long_name = 'x-coordinate of mesh nodes'
         else
            std_name  = 'projection_y_coordinate'
            long_name = 'y-coordinate of mesh nodes'
         endif
      endif
      NF90(nf90_def_var(map_file%ncid, trim(vname), nctype, (/map_file%nmesh2d_node_dimid/), varid))
      NF90(nf90_def_var_deflate(map_file%ncid, varid, 1, 1, nc_deflate_level))
      NF90(nf90_put_att(map_file%ncid, varid, 'units',         trim(units)))
      NF90(nf90_put_att(map_file%ncid, varid, 'standard_name', trim(std_name)))
      NF90(nf90_put_att(map_file%ncid, varid, 'long_name',     trim(long_name)))
      NF90(nf90_put_att(map_file%ncid, varid, 'mesh',          'mesh2d'))
      NF90(nf90_put_att(map_file%ncid, varid, 'location',      'node'))
   end subroutine def_mesh2d_node_coord

   subroutine def_grid_axis_coord(name, axis, dimids, varid, long_name)
      ! Define a regular-grid coordinate axis variable (face_x/face_y/corner_x/
      ! corner_y) with the right CF attributes:
      !   crsgeo  : NF90_FLOAT, units='degrees', standard_name=longitude/latitude
      !   else    : NF90_FLOAT, units='m',       standard_name=projection_{x,y}_coordinate
      use sfincs_data, only: crsgeo, nc_deflate_level
      character(len=*), intent(in)  :: name
      character(len=*), intent(in)  :: axis   ! 'x' or 'y'
      integer,          intent(in)  :: dimids(:)
      integer,          intent(out) :: varid
      character(len=*), intent(in)  :: long_name
      character(len=32) :: units, std_name
      !
      if (crsgeo) then
         units = 'degrees'
         if (axis == 'x') then
            std_name = 'longitude'
         else
            std_name = 'latitude'
         endif
      else
         units = 'm'
         if (axis == 'x') then
            std_name = 'projection_x_coordinate'
         else
            std_name = 'projection_y_coordinate'
         endif
      endif
      NF90(nf90_def_var(map_file%ncid, name, NF90_FLOAT, dimids, varid))
      NF90(nf90_def_var_deflate(map_file%ncid, varid, 1, 1, nc_deflate_level))
      NF90(nf90_put_att(map_file%ncid, varid, '_FillValue',    FILL_VALUE))
      NF90(nf90_put_att(map_file%ncid, varid, 'units',         trim(units)))
      NF90(nf90_put_att(map_file%ncid, varid, 'standard_name', trim(std_name)))
      NF90(nf90_put_att(map_file%ncid, varid, 'long_name',     long_name))
      NF90(nf90_put_att(map_file%ncid, varid, 'grid_mapping',  'crs'))
      NF90(nf90_put_att(map_file%ncid, varid, 'grid',          'sfincsgrid'))
   end subroutine def_grid_axis_coord

   subroutine put_2d(varid, arr)
      ! nf90_put_var for a 2D real*4 array (face/corner coords).
      ! Hides the (/1, 1/) start-index boilerplate.
      integer, intent(in) :: varid
      real*4,  intent(in) :: arr(:,:)
      NF90(nf90_put_var(map_file%ncid, varid, arr, (/1, 1/)))
   end subroutine put_2d

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! His-init def wrappers
   !
   ! def_time_point_float : (/points_dimid, time_dimid/) + pt_coord
   ! def_his_point_coord  : CRS-aware station coordinate variable
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine def_time_point_float(name, varid, units, long_name, &
                                    standard_name, cell_methods, description)
      use sfincs_data, only: nc_deflate_level
      character(*), intent(in)           :: name
      integer,      intent(out)          :: varid
      character(*), intent(in)           :: units, long_name
      character(*), intent(in), optional :: standard_name, cell_methods, description
      call ncdef_float_var(his_file%ncid, name,                           &
           (/his_file%points_dimid, his_file%time_dimid/), varid,         &
           units, long_name,                                              &
           standard_name=standard_name, coordinates=pt_coord,             &
           cell_methods=cell_methods, description=description,            &
           deflate_level=nc_deflate_level)
   end subroutine def_time_point_float

   subroutine def_his_point_coord(name, axis, varid, long_name)
      ! Define a his-file station coordinate variable with CRS-aware attrs.
      use sfincs_data, only: crsgeo
      character(len=*), intent(in)  :: name      ! 'station_x', 'point_y', etc.
      character(len=*), intent(in)  :: axis      ! 'x' or 'y'
      integer,          intent(out) :: varid
      character(len=*), intent(in)  :: long_name
      character(len=32) :: units, std_name
      !
      if (crsgeo) then
         units = 'degrees'
         if (axis == 'x') then
            std_name = 'longitude'
         else
            std_name = 'latitude'
         endif
      else
         units = 'm'
         if (axis == 'x') then
            std_name = 'projection_x_coordinate'
         else
            std_name = 'projection_y_coordinate'
         endif
      endif
      NF90(nf90_def_var(his_file%ncid, name, NF90_FLOAT, (/his_file%points_dimid/), varid))
      NF90(nf90_put_att(his_file%ncid, varid, 'units',         trim(units)))
      NF90(nf90_put_att(his_file%ncid, varid, 'standard_name', trim(std_name)))
      NF90(nf90_put_att(his_file%ncid, varid, 'long_name',     long_name))
      NF90(nf90_put_att(his_file%ncid, varid, 'grid_mapping',  'crs'))
   end subroutine def_his_point_coord

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! His-update point writer
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine write_point_var(varid, source, nt, scale)
   ! Gathers source(nmindobs(iobs)) for each obs point and writes to his_file.
   ! Optional scale is applied only to non-fill values.
      use sfincs_data, only: nobs, nmindobs
      integer, intent(in)           :: varid
      real*4,  intent(in)           :: source(:)
      integer, intent(in)           :: nt
      real*4,  intent(in), optional :: scale
      real*4  :: obs(nobs)
      integer :: iobs, nm
      obs = FILL_VALUE
      do iobs = 1, nobs
         nm = nmindobs(iobs)
         if (nm > 0) obs(iobs) = source(nm)
      enddo
      if (present(scale)) then
         where (obs /= FILL_VALUE) obs = obs * scale
      endif
      NF90(nf90_put_var(his_file%ncid, varid, obs, (/1, nt/)))
   end subroutine write_point_var

   subroutine compute_uv_at_obs_points(uobs, vobs, uvmag, uvdir)
   ! Face-averaged, rotated (u,v) and derived magnitude/direction at obs points.
      use sfincs_data, only: nobs, nmindobs, z_index_uv_md, z_index_uv_mu, &
                             z_index_uv_nd, z_index_uv_nu, uv, cosrot, sinrot, pi
      real*4, intent(out) :: uobs(:), vobs(:), uvmag(:), uvdir(:)
      integer :: iobs, nm, nmd1, nmu1, ndm1, num1
      real*4  :: uz, vz
      uobs  = FILL_VALUE
      vobs  = FILL_VALUE
      uvmag = FILL_VALUE
      uvdir = FILL_VALUE
      do iobs = 1, nobs
         nm = nmindobs(iobs)
         if (nm > 0) then
            nmd1 = z_index_uv_md(nm)
            nmu1 = z_index_uv_mu(nm)
            ndm1 = z_index_uv_nd(nm)
            num1 = z_index_uv_nu(nm)
            uz = 0.5 * (uv(nmd1) + uv(nmu1))
            vz = 0.5 * (uv(ndm1) + uv(num1))
            uobs(iobs)  = cosrot * uz - sinrot * vz
            vobs(iobs)  = sinrot * uz + cosrot * vz
            uvmag(iobs) = sqrt(uobs(iobs)**2 + vobs(iobs)**2)
            uvdir(iobs) = atan2(vobs(iobs), uobs(iobs)) * 180 / pi
         endif
      enddo
   end subroutine compute_uv_at_obs_points

   subroutine compute_wind_at_obs_points(wmag, wdir)
   ! Wind speed and meteorological direction (270 - math angle) at obs points.
      use sfincs_data, only: nobs, nmindobs, windu, windv, pi
      real*4, intent(out) :: wmag(:), wdir(:)
      integer :: iobs, nm
      wmag = FILL_VALUE
      wdir = FILL_VALUE
      do iobs = 1, nobs
         nm = nmindobs(iobs)
         if (nm > 0) then
            wmag(iobs) = sqrt(windu(nm)**2 + windv(nm)**2)
            wdir(iobs) = 270.0 - atan2(windv(nm), windu(nm)) * 180 / pi
            if (wdir(iobs) < 0.0)   wdir(iobs) = wdir(iobs) + 360.0
            if (wdir(iobs) > 360.0) wdir(iobs) = wdir(iobs) - 360.0
         endif
      enddo
   end subroutine compute_wind_at_obs_points

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Update-map precompute helpers
   !
   ! Produce SFINCS-indexed (np-shaped) source arrays so that grid-type-aware
   ! write_cell_var handles the gather/put.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine compute_uv_at_cell_centers(uxy, vxy)
      ! Cell-centered (u, v) by averaging the 4 surrounding UV faces.
      ! Boundary cells contribute half (when one of nmd1/nmu1/ndm1/num1 is 0).
      use sfincs_data
      real*4, intent(out) :: uxy(:), vxy(:)
      integer :: nm, nmd1, nmu1, ndm1, num1
      real*4  :: uz, vz
      uxy = FILL_VALUE
      vxy = FILL_VALUE
      do nm = 1, np
         nmd1 = z_index_uv_md(nm)
         nmu1 = z_index_uv_mu(nm)
         ndm1 = z_index_uv_nd(nm)
         num1 = z_index_uv_nu(nm)
         uz = 0.0
         if (nmd1 > 0) uz = uz + 0.5*uv(nmd1)
         if (nmu1 > 0) uz = uz + 0.5*uv(nmu1)
         vz = 0.0
         if (ndm1 > 0) vz = vz + 0.5*uv(ndm1)
         if (num1 > 0) vz = vz + 0.5*uv(num1)
         uxy(nm) = cosrot*uz - sinrot*vz
         vxy(nm) = sinrot*uz + cosrot*vz
      enddo
   end subroutine compute_uv_at_cell_centers

   subroutine compute_pnh_unwrapped(pnh_full)
      ! Map row-indexed nonhydrostatic pressure to SFINCS-indexed array.
      use sfincs_data, only: np
      use sfincs_nonhydrostatic, only: pnh, row_index_of_nm
      real*4, intent(out) :: pnh_full(:)
      integer :: nm
      pnh_full = FILL_VALUE
      do nm = 1, np
         if (row_index_of_nm(nm) > 0) pnh_full(nm) = pnh(row_index_of_nm(nm))
      enddo
   end subroutine compute_pnh_unwrapped

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Finalize-time single-write helpers
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine ncoutput_write_timestep_analysis()
   !
   ! Write per-cell timestep statistics once at the end of simulation
   ! (no time dimension).
   !
   use sfincs_data
   !
   implicit none
   !
   call put_static_cell_float(map_file%ncid, map_file%average_required_timestep_varid, &
        timestep_analysis_average_required_timestep_per_cell, FILL_VALUE, min_value=0.0)
   call put_static_cell_float(map_file%ncid, map_file%percentage_limiting_varid, &
        timestep_analysis_percentage_limiting_per_cell, FILL_VALUE)
   !
   end subroutine ncoutput_write_timestep_analysis

   subroutine ncoutput_write_tsunami_arrival_time()
   ! Write tsunami_arrival_time once at finalize.
   ! min_value=0.0 -> cells the tsunami never reached are written as
   ! _FillValue, not 0 (which would otherwise read as "arrived at t=0").
   use sfincs_data
   !
   implicit none
   !
   call put_static_cell_float(map_file%ncid, map_file%tsunami_arrival_time_varid, tsunami_arrival_time, &
        FILL_VALUE, min_value=0.0)
   !
   end subroutine ncoutput_write_tsunami_arrival_time

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Generic NetCDF variable definition helpers.
   ! ncdef_float_var / ncdef_int_var replace the 6-7 line def_var boilerplate
   ! that repeats across ncoutput_map_init and ncoutput_his_init.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine ncdef_float_var(ncid, varname, dimids, varid, &
                               units, long_name,             &
                               standard_name, coordinates,   &
                               cell_methods, description,    &
                               deflate_level)
      !
      ! Define a NF90_FLOAT variable and set standard attributes in one call.
      ! Optional attributes are only written when present and non-empty.
      !
      integer,      intent(in)           :: ncid
      character(*), intent(in)           :: varname
      integer,      intent(in)           :: dimids(:)
      integer,      intent(out)          :: varid
      character(*), intent(in)           :: units
      character(*), intent(in)           :: long_name
      character(*), intent(in), optional :: standard_name
      character(*), intent(in), optional :: coordinates
      character(*), intent(in), optional :: cell_methods
      character(*), intent(in), optional :: description
      integer,      intent(in), optional :: deflate_level
      !
      integer :: nc_deflate
      !
      nc_deflate = 6
      if (present(deflate_level)) nc_deflate = deflate_level
      !
      NF90(nf90_def_var(ncid, trim(varname), NF90_FLOAT, dimids, varid))
      NF90(nf90_def_var_deflate(ncid, varid, 1, 1, nc_deflate))
      NF90(nf90_put_att(ncid, varid, '_FillValue', FILL_VALUE))
      NF90(nf90_put_att(ncid, varid, 'units', trim(units)))
      if (present(standard_name)) then
         if (len_trim(standard_name) > 0) then
            NF90(nf90_put_att(ncid, varid, 'standard_name', trim(standard_name)))
         endif
      endif
      NF90(nf90_put_att(ncid, varid, 'long_name', trim(long_name)))
      if (present(coordinates)) then
         if (len_trim(coordinates) > 0) then
            NF90(nf90_put_att(ncid, varid, 'coordinates', trim(coordinates)))
         endif
      endif
      if (present(cell_methods)) then
         if (len_trim(cell_methods) > 0) then
            NF90(nf90_put_att(ncid, varid, 'cell_methods', trim(cell_methods)))
         endif
      endif
      if (present(description)) then
         if (len_trim(description) > 0) then
            NF90(nf90_put_att(ncid, varid, 'description', trim(description)))
         endif
      endif
      !
   end subroutine ncdef_float_var

   subroutine ncdef_int_var(ncid, varname, dimids, varid, &
                             long_name,                    &
                             fill_value, description,      &
                             deflate_level)
      !
      ! Define a NF90_INT variable and set standard attributes in one call.
      ! Topology variables with many unique attributes (mesh2d, crs, sfincsgrid,
      ! mesh2d_face_nodes) should still be defined manually.
      !
      integer,      intent(in)           :: ncid
      character(*), intent(in)           :: varname
      integer,      intent(in)           :: dimids(:)
      integer,      intent(out)          :: varid
      character(*), intent(in)           :: long_name
      integer,      intent(in), optional :: fill_value    ! default: FILL_VALUE_INT (-999)
      character(*), intent(in), optional :: description
      integer,      intent(in), optional :: deflate_level ! default: 6
      !
      integer :: nc_deflate, nc_fill
      !
      nc_deflate = 6
      if (present(deflate_level)) nc_deflate = deflate_level
      nc_fill = FILL_VALUE_INT
      if (present(fill_value)) nc_fill = fill_value
      !
      NF90(nf90_def_var(ncid, trim(varname), NF90_INT, dimids, varid))
      NF90(nf90_def_var_deflate(ncid, varid, 1, 1, nc_deflate))
      NF90(nf90_put_att(ncid, varid, '_FillValue', nc_fill))
      NF90(nf90_put_att(ncid, varid, 'long_name', trim(long_name)))
      if (present(description)) then
         if (len_trim(description) > 0) then
            NF90(nf90_put_att(ncid, varid, 'description', trim(description)))
         endif
      endif
      !
   end subroutine ncdef_int_var

   subroutine handle_err(status,file,line)
      ! Reports NetCDF errors to stderr. Does NOT close any file (the failing
      ! call could have been on either map_file%ncid or his_file%ncid; closing
      ! map_file unconditionally was misleading and pre-existing). Both files
      ! are closed by their respective finalize routines on normal shutdown.
      integer,      intent(in) :: status
      character(*), intent(in) :: file
      integer,      intent(in) :: line
      if (status /= nf90_noerr) then
         write(0,'("NETCDF ERROR: ",a,i6,":",a)') file, line, trim(nf90_strerror(status))
      end if
   end subroutine handle_err

end module sfincs_ncoutput_helpers
