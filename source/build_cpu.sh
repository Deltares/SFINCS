#!/bin/sh

find . -name \*.m4|xargs dos2unix && find . -name \*.ac|xargs dos2unix && find . -name \*.am|xargs dos2unix
find . -name \*.f90|xargs dos2unix
find . -name \*.F90|xargs dos2unix
find . -name \*.am|xargs dos2unix

autoreconf -vif

./autogen.sh

### ./configure FCFLAGS="-fopenmp -O3 -fallow-argument-mismatch -w" FFLAGS="-fopenmp -O3 -fallow-argument-mismatch -w" FC=gfortran --disable-shared --prefix=/usr/local/bin/sfincs/gfortran
./configure FCFLAGS="-fopenmp -O3 -fallow-argument-mismatch -w" FFLAGS="-fopenmp -O3 -fallow-argument-mismatch -w" FC=nvfortran --disable-shared --prefix=/usr/local/bin/sfincs/nvfortran_cpu

make clean

make

make install
