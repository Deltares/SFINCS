#define NF90H(nf90call) if (nf90call /= nf90_noerr) call handle_err_h(nf90call, __FILE__, __LINE__)
module sfincs_nchelper
   !
   ! Generic NetCDF helper routines for SFINCS output.
   ! ncdef_float_var / ncdef_int_var replace the 6-7 line def_var boilerplate
   ! that repeats across ncoutput_map_init and ncoutput_his_init.
   !
   use netcdf
   implicit none
   !
   real*4,    parameter :: FILL_VALUE_H   = -99999.0
   integer,   parameter :: FILL_VALUE_INT = -999
   !
contains

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
      NF90H(nf90_def_var(ncid, trim(varname), NF90_FLOAT, dimids, varid))
      NF90H(nf90_def_var_deflate(ncid, varid, 1, 1, nc_deflate))
      NF90H(nf90_put_att(ncid, varid, '_FillValue', FILL_VALUE_H))
      NF90H(nf90_put_att(ncid, varid, 'units', trim(units)))
      if (present(standard_name)) then
         if (len_trim(standard_name) > 0) then
            NF90H(nf90_put_att(ncid, varid, 'standard_name', trim(standard_name)))
         endif
      endif
      NF90H(nf90_put_att(ncid, varid, 'long_name', trim(long_name)))
      if (present(coordinates)) then
         if (len_trim(coordinates) > 0) then
            NF90H(nf90_put_att(ncid, varid, 'coordinates', trim(coordinates)))
         endif
      endif
      if (present(cell_methods)) then
         if (len_trim(cell_methods) > 0) then
            NF90H(nf90_put_att(ncid, varid, 'cell_methods', trim(cell_methods)))
         endif
      endif
      if (present(description)) then
         if (len_trim(description) > 0) then
            NF90H(nf90_put_att(ncid, varid, 'description', trim(description)))
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
      NF90H(nf90_def_var(ncid, trim(varname), NF90_INT, dimids, varid))
      NF90H(nf90_def_var_deflate(ncid, varid, 1, 1, nc_deflate))
      NF90H(nf90_put_att(ncid, varid, '_FillValue', nc_fill))
      NF90H(nf90_put_att(ncid, varid, 'long_name', trim(long_name)))
      if (present(description)) then
         if (len_trim(description) > 0) then
            NF90H(nf90_put_att(ncid, varid, 'description', trim(description)))
         endif
      endif
      !
   end subroutine ncdef_int_var

   subroutine handle_err_h(status, file, line)
      integer,      intent(in) :: status
      character(*), intent(in) :: file
      integer,      intent(in) :: line
      if (status /= nf90_noerr) then
         write(0,'("NETCDF ERROR: ",a,i6,":",a)') file, line, trim(nf90_strerror(status))
      endif
   end subroutine handle_err_h

end module sfincs_nchelper
