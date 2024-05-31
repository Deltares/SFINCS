# Install script for directory: C:/work/netcdf/netcdf-fortran-main-4.6.1/netcdf-fortran-main/fortran

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files/NC4F")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/objects-Debug/netcdff_c" TYPE FILE FILES
      "nf_lib.obj"
      "nf_v2compat.obj"
      FILES_FROM_DIR "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/fortran/netcdff_c.dir/Debug/")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/objects-Release/netcdff_c" TYPE FILE FILES
      "nf_lib.obj"
      "nf_v2compat.obj"
      FILES_FROM_DIR "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/fortran/netcdff_c.dir/Release/")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/objects-MinSizeRel/netcdff_c" TYPE FILE FILES
      "nf_lib.obj"
      "nf_v2compat.obj"
      FILES_FROM_DIR "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/fortran/netcdff_c.dir/MinSizeRel/")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/objects-RelWithDebInfo/netcdff_c" TYPE FILE FILES
      "nf_lib.obj"
      "nf_v2compat.obj"
      FILES_FROM_DIR "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/fortran/netcdff_c.dir/RelWithDebInfo/")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "shlib" OR NOT CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/fortran/Debug/netcdff.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/fortran/Release/netcdff.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/fortran/MinSizeRel/netcdff.lib")
  elseif(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/fortran/RelWithDebInfo/netcdff.lib")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/fortran/" FILES_MATCHING REGEX "/[^/]*\\.mod$" REGEX "/netcdf\\.inc$")
endif()

