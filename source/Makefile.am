SUBDIRS = \
	third_party_open/netcdf/netcdf-fortran-4.6.1 \
	src

# our m4 stuff is here
ACLOCAL_AMFLAGS = -I m4
EXTRA_DIST = m4 README.md

# licences of bundled code, which no Makefile.am lists
EXTRA_DIST += \
	src/snapwave/LICENSE \
	third_party_open/utils/toml-f/LICENSE-Apache \
	third_party_open/utils/toml-f/LICENSE-MIT

CFLAGS="$CFLAGS $OPENMP_CFLAGS"

.PHONY: test
test:
	cd test/01_noadv && ../../src/sfincs
