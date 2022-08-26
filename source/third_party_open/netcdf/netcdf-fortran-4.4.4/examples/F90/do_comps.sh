#!/bin/sh
# This shell script runs the cmp test on the example programs.
# $Id: do_comps.sh 59820 2018-08-15 12:49:57Z markus $

set -e
echo ""
echo "*** Testing that F90 examples produced same files as C examples."
echo "*** checking simple_xy.nc..."
cmp simple_xy.nc ../C/simple_xy.nc

echo "*** checking sfc_pres_temp.nc..."
cmp sfc_pres_temp.nc ../C/sfc_pres_temp.nc

echo "*** checking pres_temp_4D.nc..."
cmp pres_temp_4D.nc ../C/pres_temp_4D.nc

echo "*** All F90 example comparisons worked!"
exit 0
