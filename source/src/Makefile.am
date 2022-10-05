#-------------------------------------------------------------------------------
#   Makefile for sfincs
#
#   Fedor.Baart@deltares.nl
#   11 nov 2020
#
#   Copyright (C)  Stichting Deltares, 2011-2013.
#-------------------------------------------------------------------------------

bin_PROGRAMS = sfincs
lib_LTLIBRARIES = libsfincs.la


#-------------------------------------------------------------------------------

#all sources for sfincs that go into the library (all but the program)
libsfincs_la_SOURCES = \
	sfincs_date.f90 \
	sfincs_spiderweb.f90 \
	sfincs_data.f90 \
	geometry.f90 \
	quadtree.f90 \
	snapwave/interp.F90 \
	snapwave/snapwave_results.f90 \
	snapwave/snapwave_data.f90 \
	snapwave/snapwave_boundaries.f90 \
	snapwave/snapwave_date.f90 \
	snapwave/snapwave_ncoutput.F90 \
	snapwave/snapwave_domain.f90 \
	snapwave/snapwave_input.f90 \
	snapwave/snapwave_solver.f90 \
	sfincs_input.f90 \
	sfincs_ncinput.F90 \
	sfincs_boundaries.f90 \
	sfincs_continuity.f90 \
	sfincs_crosssections.f90 \
	sfincs_discharges.f90 \
	sfincs_domain.f90 \
	sfincs_structures.f90 \
	sfincs_obspoints.f90 \
	sfincs_snapwave.f90 \
	deg2utm.f90 \
	sfincs_meteo.f90 \
	sfincs_ncoutput.F90 \
	sfincs_output.f90 \
	sfincs_momentum.f90 \
	sfincs_wavemaker.f90 \
	sfincs_lib.f90 \
	sfincs_bmi.f90
        
libsfincs_la_LIBADD = \
	$(top_builddir)/third_party_open/netcdf/netcdf-fortran-4.4.4/fortran/libnetcdff.la
# these files are not included
# sfincs_lib_nonc.f90
# sfincs_output_nonc.f90
# sfincs_boundaries_nonc.f90
# sfincs_meteo_nonc.f90
# sfincs_nonc.f90


# only the sources for the executable (with program)
sfincs_SOURCES = sfincs.f90
# add the sfincs lib as a dependency
# TODO: add netcdf-fortran
sfincs_LDADD = $(builddir)/libsfincs.la
sfincs_LDFLAGS = -static
sfincs_DEPENDENCIES = libsfincs.la

# make sure the mod file is there
sfincs.f90: sfincs_lib.mod
sfincs_lib.mod: libsfincs.la


AM_FCFLAGS=${NETCDF_CFLAGS} ${NETCDF_FORTRAN_CFLAGS} -I$(top_builddir)/third_party_open/netcdf/netcdf-fortran-4.4.4/fortran
AM_FCFLAGS+=${OPENMP_FCFLAGS}
AM_FCFLAGS+=${OPENACC_CFLAGS}
AM_FCFLAGS+=${-O3}


# Manualy add the netcdff (fortran dll), bug in pkg-config file is not fixed until in NetCDF 4.4 HEAD
LIBS+=${NETCDF_LIBS} ${NETCDF_FORTRAN_LIBS}
