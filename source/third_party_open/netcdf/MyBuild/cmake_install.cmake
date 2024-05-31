# Install script for directory: C:/work/netcdf/netcdf-fortran-main-4.6.1/netcdf-fortran-main

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

if(CMAKE_INSTALL_COMPONENT STREQUAL "utilities" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/netcdf-fortran.pc")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "utilities" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/nf-config")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "libraries" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE FILE FILES "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/libnetcdff.settings")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/netCDF-Fortran/netcdffTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/netCDF-Fortran/netcdffTargets.cmake"
         "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/CMakeFiles/Export/0e98d0e39643b94ae06c9e3e14d8e119/netcdffTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/netCDF-Fortran/netcdffTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/netCDF-Fortran/netcdffTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/netCDF-Fortran" TYPE FILE FILES "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/CMakeFiles/Export/0e98d0e39643b94ae06c9e3e14d8e119/netcdffTargets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/netCDF-Fortran" TYPE FILE FILES "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/CMakeFiles/Export/0e98d0e39643b94ae06c9e3e14d8e119/netcdffTargets-debug.cmake")
  endif()
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Mm][Ii][Nn][Ss][Ii][Zz][Ee][Rr][Ee][Ll])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/netCDF-Fortran" TYPE FILE FILES "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/CMakeFiles/Export/0e98d0e39643b94ae06c9e3e14d8e119/netcdffTargets-minsizerel.cmake")
  endif()
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/netCDF-Fortran" TYPE FILE FILES "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/CMakeFiles/Export/0e98d0e39643b94ae06c9e3e14d8e119/netcdffTargets-relwithdebinfo.cmake")
  endif()
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/netCDF-Fortran" TYPE FILE FILES "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/CMakeFiles/Export/0e98d0e39643b94ae06c9e3e14d8e119/netcdffTargets-release.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "headers" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/netCDF-Fortran" TYPE FILE FILES
    "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/netCDF-FortranConfig.cmake"
    "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/netCDF-FortranConfigVersion.cmake"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/fortran/cmake_install.cmake")
  include("C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/libsrc/cmake_install.cmake")
  include("C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/nf_test/cmake_install.cmake")
  include("C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/nf03_test/cmake_install.cmake")
  include("C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/nf_test4/cmake_install.cmake")
  include("C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/nf03_test4/cmake_install.cmake")
  include("C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/examples/cmake_install.cmake")
  include("C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/docs/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "C:/work/netcdf/netcdf-fortran-main-4.6.1/MyBuild/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
