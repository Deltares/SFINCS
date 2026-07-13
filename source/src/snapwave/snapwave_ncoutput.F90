#define NF90(nf90call) call handle_err(nf90call,__FILE__,__LINE__)   
module snapwave_ncoutput
   !
   use sfincs_log
   use netcdf       
   !
   implicit none
   !
   contains
   !
   subroutine write_snapwave_mesh(fname, crsgeo)
      !
      use snapwave_data
      !
      implicit none
      !
      character(len=256), intent(in) :: fname
      logical, intent(in)            :: crsgeo
      ! 
      integer :: ncid
      integer :: nmesh2d_node_dimid, nmesh2d_face_dimid, max_nmesh2d_face_nodes_dimid
      integer :: mesh2d_varid
      integer :: mesh2d_node_x_varid, mesh2d_node_y_varid, crs_varid
      integer :: mesh2d_face_nodes_varid
      integer :: zb_varid
      !
      integer, parameter :: nc_deflate_level = 2
      real*4, parameter  :: FILL_VALUE = -99999.0
      !
      ! dimensions
      NF90(nf90_create(trim(fname), ior(NF90_CLOBBER, NF90_NETCDF4), ncid))
      NF90(nf90_def_dim(ncid, 'nmesh2d_node', no_nodes, nmesh2d_node_dimid))
      NF90(nf90_def_dim(ncid, 'nmesh2d_face', no_faces, nmesh2d_face_dimid))
      NF90(nf90_def_dim(ncid, 'max_nmesh2d_face_nodes', 4, max_nmesh2d_face_nodes_dimid))
      !
      ! global attributes
      NF90(nf90_put_att(ncid,nf90_global, "Conventions", "Conventions = 'CF-1.8 UGRID-1.0 Deltares-0.10'"))
      NF90(nf90_put_att(ncid,nf90_global, "Build-Revision-Date-Netcdf-library", trim(nf90_inq_libvers())))
      NF90(nf90_put_att(ncid,nf90_global, "Producer", "SFINCS model: Super-Fast INundation of CoastS"))
      NF90(nf90_put_att(ncid,nf90_global, "Build-Revision", trim(build_revision))) 
      NF90(nf90_put_att(ncid,nf90_global, "Build-Date", trim(build_date)))
      NF90(nf90_put_att(ncid,nf90_global, "title", "Snapwave grid"))
      !
      ! mesh topology
      NF90(nf90_def_var(ncid, 'mesh2d', NF90_INT, mesh2d_varid))
      NF90(nf90_put_att(ncid, mesh2d_varid, 'cf_role', 'mesh_topology'))
      NF90(nf90_put_att(ncid, mesh2d_varid, 'long_name', 'Topology data of 2D network'))
      NF90(nf90_put_att(ncid, mesh2d_varid, 'topology_dimension', 2))
      NF90(nf90_put_att(ncid, mesh2d_varid, 'node_coordinates', 'mesh2d_node_x mesh2d_node_y'))
      NF90(nf90_put_att(ncid, mesh2d_varid, 'node_dimension', 'nmesh2d_node'))
      NF90(nf90_put_att(ncid, mesh2d_varid, 'max_face_nodes_dimension', 'max_nmesh2d_face_nodes'))
      NF90(nf90_put_att(ncid, mesh2d_varid, 'face_node_connectivity', 'mesh2d_face_nodes'))
      NF90(nf90_put_att(ncid, mesh2d_varid, 'face_dimension', 'nmesh2d_face'))
      !
      if (crsgeo) then
         !
         NF90(nf90_def_var(ncid, 'mesh2d_node_x', NF90_FLOAT, (/nmesh2d_node_dimid/), mesh2d_node_x_varid)) ! location of zb, zs etc. in cell centre
         NF90(nf90_def_var_deflate(ncid, mesh2d_node_x_varid, 1, 1, nc_deflate_level))
         NF90(nf90_put_att(ncid, mesh2d_node_x_varid, 'units', 'degrees'))
         NF90(nf90_put_att(ncid, mesh2d_node_x_varid, 'standard_name', 'longitude'))
         NF90(nf90_put_att(ncid, mesh2d_node_x_varid, 'long_name', 'longitude'))
         NF90(nf90_put_att(ncid, mesh2d_node_x_varid, 'mesh', 'mesh2d'))
         NF90(nf90_put_att(ncid, mesh2d_node_x_varid, 'location', 'node'))
         !
         NF90(nf90_def_var(ncid, 'mesh2d_node_y', NF90_FLOAT, (/nmesh2d_node_dimid/), mesh2d_node_y_varid)) ! location of zb, zs etc. in cell centre
         NF90(nf90_def_var_deflate(ncid, mesh2d_node_y_varid, 1, 1, nc_deflate_level))
         NF90(nf90_put_att(ncid, mesh2d_node_y_varid, 'units', 'degrees'))
         NF90(nf90_put_att(ncid, mesh2d_node_y_varid, 'standard_name', 'latitude'))
         NF90(nf90_put_att(ncid, mesh2d_node_y_varid, 'long_name', 'latitude'))
         NF90(nf90_put_att(ncid, mesh2d_node_y_varid, 'mesh', 'mesh2d'))
         NF90(nf90_put_att(ncid, mesh2d_node_y_varid, 'location', 'node'))
         !
      else
         !
         NF90(nf90_def_var(ncid, 'mesh2d_node_x', NF90_DOUBLE, (/nmesh2d_node_dimid/), mesh2d_node_x_varid)) ! location of zb, zs etc. in cell centre
         NF90(nf90_def_var_deflate(ncid, mesh2d_node_x_varid, 1, 1, nc_deflate_level))
         NF90(nf90_put_att(ncid, mesh2d_node_x_varid, 'units', 'm'))
         NF90(nf90_put_att(ncid, mesh2d_node_x_varid, 'standard_name', 'projection_x_coordinate'))
         NF90(nf90_put_att(ncid, mesh2d_node_x_varid, 'long_name', 'x-coordinate of mesh nodes'))
         NF90(nf90_put_att(ncid, mesh2d_node_x_varid, 'mesh', 'mesh2d'))
         NF90(nf90_put_att(ncid, mesh2d_node_x_varid, 'location', 'node'))
         !
         NF90(nf90_def_var(ncid, 'mesh2d_node_y', NF90_DOUBLE, (/nmesh2d_node_dimid/), mesh2d_node_y_varid)) ! location of zb, zs etc. in cell centre
         NF90(nf90_def_var_deflate(ncid, mesh2d_node_y_varid, 1, 1, nc_deflate_level))
         NF90(nf90_put_att(ncid, mesh2d_node_y_varid, 'units', 'm'))
         NF90(nf90_put_att(ncid, mesh2d_node_y_varid, 'standard_name', 'projection_y_coordinate'))
         NF90(nf90_put_att(ncid, mesh2d_node_y_varid, 'long_name', 'y-coordinate of mesh nodes'))
         NF90(nf90_put_att(ncid, mesh2d_node_y_varid, 'mesh', 'mesh2d'))
         NF90(nf90_put_att(ncid, mesh2d_node_y_varid, 'location', 'node'))
         !
      endif
      !
      NF90(nf90_def_var(ncid, 'mesh2d_face_nodes', NF90_INT, (/max_nmesh2d_face_nodes_dimid, nmesh2d_face_dimid/), mesh2d_face_nodes_varid)) ! location of zb, zs etc. in cell centre
      NF90(nf90_def_var_deflate(ncid, mesh2d_face_nodes_varid, 1, 1, nc_deflate_level))
      NF90(nf90_put_att(ncid, mesh2d_face_nodes_varid, 'cf_role', 'face_node_connectivity'))
      NF90(nf90_put_att(ncid, mesh2d_face_nodes_varid, 'mesh', 'mesh2d'))
      NF90(nf90_put_att(ncid, mesh2d_face_nodes_varid, 'location', 'face'))
      NF90(nf90_put_att(ncid, mesh2d_face_nodes_varid, 'long_name', 'Mapping from every face to its corner nodes (counterclockwise)'))
      NF90(nf90_put_att(ncid, mesh2d_face_nodes_varid, 'start_index', 1))
      NF90(nf90_put_att(ncid, mesh2d_face_nodes_varid, '_FillValue', -999))
      !
      NF90(nf90_def_var(ncid, 'crs', NF90_INT, crs_varid)) ! For EPSG code
      NF90(nf90_put_att(ncid, crs_varid, 'EPSG', '-'))
      !
      NF90(nf90_def_var(ncid, 'mesh2d_node_z', NF90_FLOAT, (/nmesh2d_node_dimid/), zb_varid)) ! bed level in cell centre
      NF90(nf90_def_var_deflate(ncid, zb_varid, 1, 1, nc_deflate_level))
      NF90(nf90_put_att(ncid, zb_varid, '_FillValue', FILL_VALUE))   
      NF90(nf90_put_att(ncid, zb_varid, 'units', 'm'))
      NF90(nf90_put_att(ncid, zb_varid, 'standard_name', 'altitude'))
      NF90(nf90_put_att(ncid, zb_varid, 'long_name', 'bed_level_above_reference_level'))
      !
      NF90(nf90_enddef(ncid))
      !
      ! put variables
      NF90(nf90_put_var(ncid, mesh2d_node_x_varid, x)) ! write node x 
      NF90(nf90_put_var(ncid, mesh2d_node_y_varid, y)) ! write node y
      NF90(nf90_put_var(ncid, mesh2d_face_nodes_varid, face_nodes))
      NF90(nf90_put_var(ncid, zb_varid, zb))
      
      ! close file
      NF90(nf90_close(ncid))
   
   end subroutine write_snapwave_mesh
   
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine handle_err(status,file,line)
      !
      integer, intent ( in)    :: status
      character(*), intent(in) :: file
      integer, intent ( in)    :: line
      integer :: status2
   
      if(status /= nf90_noerr) then
      !   !UNIT=6 for stdout and UNIT=0 for stderr.
         write(0,'("NETCDF ERROR: ",a,i6,":",a)') file,line,trim(nf90_strerror(status))
      end if
   end subroutine handle_err
   !
   end module    