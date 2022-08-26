! nf_par_dummy.f90 --
!    Writing NetCDF4 files in parallel is not supported on Windows it seems.
!    The file nf_nc4.f90 defines such routines anyway and we really need it.
!    So provide dummies instead for the routines that are currently missing.
!
function nc_create_par_fortran( string, mode, comm, info, ncid ) bind(C)
    use iso_c_binding
    character(len=1,kind=c_char), dimension(*) :: string
    integer(kind=c_int), value                 :: mode, comm, info, ncid
    integer(kind=c_int)                        :: nc_create_par_fortran

    nc_create_par_fortran = -999

    write(*,*) 'nc_create_par_fortran: parallel reads/writes not supported!'

end function nc_create_par_fortran

function nc_open_par_fortran( string, mode, comm, info, ncid ) bind(C)
    use iso_c_binding
    character(len=1,kind=c_char), dimension(*) :: string
    integer(kind=c_int), value                 :: mode, comm, info, ncid
    integer(kind=c_int)                        :: nc_open_par_fortran

    nc_open_par_fortran = -999

    write(*,*) 'nc_open_par_fortran: parallel reads/writes not supported!'

end function nc_open_par_fortran

function nc_var_par_access( ncid, varid, access ) bind(C)
    use iso_c_binding
    integer(kind=c_int), value                 :: ncid, varid, access
    integer(kind=c_int)                        :: nc_var_par_access

    nc_var_par_access = -999

    write(*,*) 'nc_var_par_access: parallel reads/writes not supported!'

end function nc_var_par_access


