#----------------------------------------------------------------
# Generated CMake target import file for configuration "RelWithDebInfo".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "netCDF::netcdff" for configuration "RelWithDebInfo"
set_property(TARGET netCDF::netcdff APPEND PROPERTY IMPORTED_CONFIGURATIONS RELWITHDEBINFO)
set_target_properties(netCDF::netcdff PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELWITHDEBINFO "C;Fortran"
  IMPORTED_LOCATION_RELWITHDEBINFO "${_IMPORT_PREFIX}/lib/netcdff.lib"
  )

list(APPEND _cmake_import_check_targets netCDF::netcdff )
list(APPEND _cmake_import_check_files_for_netCDF::netcdff "${_IMPORT_PREFIX}/lib/netcdff.lib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
