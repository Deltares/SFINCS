#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "netCDF::netcdff" for configuration "Release"
set_property(TARGET netCDF::netcdff APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(netCDF::netcdff PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C;Fortran"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/netcdff.lib"
  )

list(APPEND _cmake_import_check_targets netCDF::netcdff )
list(APPEND _cmake_import_check_files_for_netCDF::netcdff "${_IMPORT_PREFIX}/lib/netcdff.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
