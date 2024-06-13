#!/bin/sh

module load nvidia/nvhpc/24.1
module load netcdf

export LDFLAGS="-L${NETCDF_C_LIBRARY}"

find . -name \*.m4|xargs dos2unix && find . -name \*.ac|xargs dos2unix && find . -name \*.am|xargs dos2unix
find . -name \*.f90|xargs dos2unix
find . -name \*.F90|xargs dos2unix
find . -name \*.am|xargs dos2unix

#MANPATH=$MANPATH:/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/man; export MANPATH
#PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.5/compilers/bin:$PATH; export PATH

#apt install -y build-essential autoconf automake libtool pkg-config tzdata
#export CONFIG_SHELL=/bin/bash

#autoreconf -vif

./autogen.sh

./configure FCFLAGS="-acc -Minfo=accel -fast -O3 -gpu=ccall -DSIZEOF_PTRDIFF_T=999" FC=nvfortran --disable-shared --disable-openmp --prefix bin/sfincs

make clean

make

make install
