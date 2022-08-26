#!/bin/sh
# This shell script which creates the fill.nc file from fill.cdl.
# $Id: create_fills.sh 59820 2018-08-15 12:49:57Z markus $

echo
echo "*** Testing creating file with fill values."
set -e
#../ncgen/ncgen -b $srcdir/fills.cdl
cp ${TOPSRCDIR}/nf_test/ref_fills.nc ./fills.nc
echo "*** SUCCESS!"
exit 0
