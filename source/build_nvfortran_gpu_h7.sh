#!/bin/sh

module load nvidia/nvhpc/24.1
module load netcdf

export LDFLAGS="-L${NETCDF_C_LIBRARY}"

find . -name \*.m4|xargs dos2unix && find . -name \*.ac|xargs dos2unix && find . -name \*.am|xargs dos2unix
find . -name \*.f90|xargs dos2unix
find . -name \*.F90|xargs dos2unix
find . -name \*.am|xargs dos2unix

./autogen.sh

./configure FCFLAGS="-acc -Minfo=accel -fast -O3 -gpu=ccall -DSIZEOF_PTRDIFF_T=999" FC=nvfortran --disable-shared --disable-openmp --prefix /u/${USER}/bin/sfincs_nvfortran_gpu

make clean

make

make install
