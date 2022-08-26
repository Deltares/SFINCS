#!/bin/sh

MANPATH=$MANPATH:/opt/nvidia/hpc_sdk/Linux_x86_64/21.11/compilers/man; export MANPATH
PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/21.11/compilers/bin:$PATH; export PATH

apt install -y libnetcdf-dev build-essential autoconf automake libtool pkg-config tzdata
export CONFIG_SHELL=/bin/bash

autoreconf -vif

./autogen.sh

./configure FCFLAGS="-acc -Minfo=accel -fast -O3 -gpu=ccall" FC=nvfortran --disable-shared --disable-openmp --program-suffix="_no_async"

make clean

make

make install
