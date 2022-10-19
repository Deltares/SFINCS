#!/bin/sh

# module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc/20.11
# module load nvhpc/20.11

MANPATH=$MANPATH:/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/compilers/man; export MANPATH
PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/compilers/bin:$PATH; export PATH

apt install -y libnetcdf-dev build-essential autoconf automake libtool pkg-config tzdata
export CONFIG_SHELL=/bin/bash

autoreconf -vif

./autogen.sh

./configure FCFLAGS="-acc -Minfo=accel -fast -O3 -gpu=ccall -DSIZEOF_PTRDIFF_T=999" FC=nvfortran --disable-shared --disable-openmp --program-suffix="_async"

make clean

make

make install
