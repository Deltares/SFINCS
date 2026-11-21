#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
dnl we are using autoconf
AC_INIT([SFINCS], [0.1], [tim.leijnse@deltares.nl])
dnl automake
AM_INIT_AUTOMAKE([subdir-objects])
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

AC_LANG_PUSH([Fortran])
dnl detect openacc flag (only with --enable-openacc)
AX_OPENACC
AS_IF([test "x$enable_openacc" = xyes && test "x$OPENACC_FCFLAGS" = x],
      [AC_MSG_ERROR([--enable-openacc was given, but no OpenACC flag works with $FC])])
AC_MSG_CHECKING([whether to build with OpenACC])
AS_IF([test "x$OPENACC_FCFLAGS" != x],
      [AC_MSG_RESULT([yes ($OPENACC_FCFLAGS)])],
      [AC_MSG_RESULT([no])])
dnl the OpenACC and OpenMP directives in the code cannot be compiled together,
dnl so an OpenACC build is built without OpenMP
AS_IF([test "x$OPENACC_FCFLAGS" != x],
      [AS_IF([test "x$enable_openmp" = xyes],
             [AC_MSG_ERROR([--enable-openacc cannot be combined with --enable-openmp])])
       enable_openmp=no])
dnl detect openmp
AC_OPENMP
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
AC_CONFIG_SUBDIRS([third_party_open/netcdf/netcdf-fortran-4.6.1])
AC_OUTPUT

AS_IF([test "x$OPENACC_FCFLAGS" != x],
      [AC_MSG_NOTICE([OpenACC: enabled ($OPENACC_FCFLAGS), OpenMP disabled; set the GPU target in FCFLAGS, e.g. -gpu=cc86])],
      [AC_MSG_NOTICE([OpenACC: disabled (use --enable-openacc to enable)])])
