#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
dnl we are using autoconf
AC_INIT([SFINCS], [0.1], [tim.leijnse@deltares.nl])
dnl automake
AM_INIT_AUTOMAKE
dnl and libtool
LT_INIT

dnl this is how to find the source directory
AC_CONFIG_SRCDIR([src/sfincs.f90])
dnl this is th file with header flags, you can use this file in .F90 files
AC_CONFIG_HEADERS([config.h])
dnl m4 macro's in m4 directory.
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
dnl we need a fortran compiler
AC_PROG_FC
dnl also the c compiler for OpenACC/OpenMP checks
AC_PROG_CC
dnl and the preprocessor
AC_PROG_CPP
dnl just to make clear that we have f90 as extension
AC_FC_SRCEXT([f90])
dnl Yes, we want free form Fortran!
AC_FC_FREEFORM
AC_FC_LINE_LENGTH([unlimited])

dnl detect openmp
AC_LANG_PUSH([Fortran])
AC_OPENMP
dnl detect openacc flag
AX_OPENACC
AC_LANG_POP


# NetCDF, the c version
dnl the fortran version is compiled in place so we can link it statically
PKG_CHECK_MODULES(NETCDF, netcdf >= 4)


# Checks for libraries.

# Checks for header files.
dnl we might need to add netcdf here later
AC_CHECK_HEADERS([])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
dnl for later
AC_CHECK_TYPES([])

# Checks for library functions.
dnl for debugging
AC_FUNC_ERROR_AT_LINE
dnl we might need to add netcdf function here to check if it works
AC_CHECK_FUNCS([])

dnl these are our makefiles
AC_CONFIG_FILES([
        Makefile
        src/Makefile
])
dnl we'll try to include netcdf-fortran in the binary so all fortran code is portable
AC_CONFIG_SUBDIRS([third_party_open/netcdf/netcdf-fortran-4.4.4])
AC_OUTPUT
