rm -rf build;mkdir build;cd build
 cmake ..   -DSFINCS_ENABLE_NETCDF=ON   -DNETCDF_PREFIX=/usr -DCMAKE_Fortran_FLAGS="-w";make -j
 ctest -V -R test_sfincs_bmi2
 cd ..
