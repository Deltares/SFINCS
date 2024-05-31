# NetCDF CXX Configuration Summary

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was netCDF-FortranConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

include(CMakeFindDependencyMacro)

if (0)
  if(EXISTS "")
    set(netCDF_ROOT "")
  endif()
  if(EXISTS "netCDF_DIR-NOTFOUND")
    set(netCDF_DIR "netCDF_DIR-NOTFOUND")
  endif()
  find_dependency(netCDF)
  set(NETCDF_C_LIBRARY ${netCDF_LIBRARIES})
  set(NETCDF_C_INCLUDE_DIR ${netCDF_INCLUDE_DIR})
else()
  set(NETCDF_C_LIBRARY "c:/work/netcdf/netCDF 4.9.2/lib/netcdf.lib")
  set(NETCDF_C_INCLUDE_DIR "c:/work/netcdf/netCDF 4.9.2/include")
endif()

if (NOT TARGET netCDF::netcdf)
  add_library(netCDF::netcdf UNKNOWN IMPORTED)
  set_target_properties(netCDF::netcdf PROPERTIES
    IMPORTED_LOCATION "${NETCDF_C_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_C_INCLUDE_DIR}"
  )
endif()

include("${CMAKE_CURRENT_LIST_DIR}/netcdffTargets.cmake")
