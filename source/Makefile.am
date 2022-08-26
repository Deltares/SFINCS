SUBDIRS = \
	third_party_open/netcdf/netcdf-fortran-4.4.4 \
	src

# our m4 stuff is here
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = m4

CFLAGS="$CFLAGS $OPENMP_CFLAGS"

# you might want to add these files to confirm with standard linux software
README: README.md
NEWS:
AUTHORS:
ChangeLog:

.PHONY: test
test:
	cd test/01_noadv && ../../src/sfincs


