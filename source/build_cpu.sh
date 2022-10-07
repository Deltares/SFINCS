#!/bin/sh

autoreconf -vif

./autogen.sh

./configure FCFLAGS="-fopenmp -O3 -fallow-argument-mismatch -w" FFLAGS="-fopenmp -O3 -fallow-argument-mismatch -w" FC=gfortran --disable-shared

make clean

make

make install
