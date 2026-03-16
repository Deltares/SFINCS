cmake_minimum_required(VERSION 3.19)
project(sfincs LANGUAGES Fortran)

# ---------------- Global PIC for shared libs ----------------
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ---------------- Options ----------------
option(SFINCS_ENABLE_NETCDF "Enable NetCDF-Fortran I/O" ON)
option(SFINCS_ENABLE_OPENMP "Enable OpenMP for Fortran" OFF)
option(SFINCS_BUILD_TESTS "Build SFINCS BMI tests" OFF)

# ---------------- Fortran standard ----------------
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------- Source root ----------------
set(SFINCS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Directory with SFINCS Fortran sources")

function(_append_if_exists outvar)
  foreach(_f IN LISTS ARGN)
    if(IS_ABSOLUTE "${_f}")
      if(EXISTS "${_f}")
        list(APPEND ${outvar} "${_f}")
      else()
        message(WARNING "SFINCS: missing source (skipped): ${_f}")
      endif()
    else()
      if(EXISTS "${SFINCS_SOURCE_DIR}/${_f}")
        list(APPEND ${outvar} "${SFINCS_SOURCE_DIR}/${_f}")
      else()
        message(WARNING "SFINCS: missing source (skipped): ${SFINCS_SOURCE_DIR}/${_f}")
      endif()
    endif()
  endforeach()
  set(${outvar} "${${outvar}}" PARENT_SCOPE)
endfunction()

# Recursive finder for merged Deltares layouts where files may have moved
function(_append_first_match outvar)
  foreach(_name IN LISTS ARGN)
    unset(_hits)
    file(GLOB_RECURSE _hits
      "${SFINCS_SOURCE_DIR}/${_name}"
      "${SFINCS_SOURCE_DIR}/*/${_name}"
      "${SFINCS_SOURCE_DIR}/*/*/${_name}"
      "${SFINCS_SOURCE_DIR}/*/*/*/${_name}"
    )
    list(LENGTH _hits _nhits)
    if(_nhits GREATER 0)
      list(GET _hits 0 _picked)
      list(APPEND ${outvar} "${_picked}")
    else()
      message(WARNING "SFINCS: missing source (skipped): ${_name}")
    endif()
  endforeach()
  list(REMOVE_DUPLICATES ${outvar})
  set(${outvar} "${${outvar}}" PARENT_SCOPE)
endfunction()

# ---------------- Compiler flags ----------------
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(_F_FLAGS
    "-Wall" "-Wextra" "-Wno-unused-parameter" "-Wno-unused-variable"
    "-ffree-line-length-none"
    "-Wno-line-truncation" "-Wno-error=line-truncation"
  )
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(_F_FLAGS "-warn" "all")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
  set(_F_FLAGS "-Mextend")
endif()

# ====================================================================
# NetCDF (manual & simple)
# ====================================================================
set(NETCDF_PREFIX "" CACHE PATH "Prefix for NetCDF (C & Fortran)")
set(NETCDFF_FOUND FALSE)

if(SFINCS_ENABLE_NETCDF)
  if(NOT NETCDF_PREFIX AND DEFINED ENV{NETCDF_PREFIX})
    set(NETCDF_PREFIX "$ENV{NETCDF_PREFIX}")
  endif()

  find_library(NETCDFF_LIB NAMES netcdff HINTS ${NETCDF_PREFIX}/lib ${NETCDF_PREFIX}/lib64)
  find_library(NETCDF_C_LIB NAMES netcdf HINTS ${NETCDF_PREFIX}/lib ${NETCDF_PREFIX}/lib64)
  find_path(NETCDF_INCLUDE_DIR NAMES netcdf.h netcdf.inc HINTS ${NETCDF_PREFIX}/include)
  find_path(NETCDF_MOD_DIR NAMES netcdf.mod typesizes.mod
    HINTS
      ${NETCDF_PREFIX}/include
      ${NETCDF_PREFIX}/include/gfortran
      ${NETCDF_PREFIX}/include/gfortran/modules
      ${NETCDF_PREFIX}/include/ifort
      ${NETCDF_PREFIX}/lib
      ${NETCDF_PREFIX}/lib64)

  if(NETCDFF_LIB AND NETCDF_INCLUDE_DIR)
    set(NETCDFF_FOUND TRUE)
  else()
    message(FATAL_ERROR
      "NetCDF-Fortran not found.\n"
      "Set -DNETCDF_PREFIX=/opt/netcdf... or build with -DSFINCS_ENABLE_NETCDF=OFF.")
  endif()
endif()

# ====================================================================
# iso_c_fortran_bmi (REQUIRED for NGen Fortran BMI ABI)
# ====================================================================
set(_ISO_DEFAULT "${CMAKE_CURRENT_LIST_DIR}/../../../iso_c_fortran_bmi")
set(ISO_C_FORTRAN_BMI_DIR "${_ISO_DEFAULT}" CACHE PATH "Path to iso_c_fortran_bmi source tree")

if(EXISTS "${ISO_C_FORTRAN_BMI_DIR}/CMakeLists.txt")
  add_subdirectory("${ISO_C_FORTRAN_BMI_DIR}" iso_c_fortran_bmi_build)
else()
  message(FATAL_ERROR "ISO_C_FORTRAN_BMI_DIR does not contain a CMakeLists.txt: ${ISO_C_FORTRAN_BMI_DIR}")
endif()

# ====================================================================
# CORE SOURCES
# ====================================================================
set(_CORE_CANDIDATES
  bmi.f90
  deg2utm.f90
  astro.f90
  sfincs_spiderweb.f90
  sfincs_error.f90
  sfincs_log.f90
  sfincs_input.f90
  sfincs_initial_conditions.f90
  sfincs_domain.f90
  sfincs_structures.f90
  sfincs_boundaries.f90
  sfincs_obspoints.f90
  sfincs_crosssections.f90
  sfincs_discharges.f90
  sfincs_meteo.f90
  sfincs_runup_gauges.f90
  geometry.f90
  quadtree.F90
  sfincs_subgrid.F90
  sfincs_infiltration.f90
  sfincs_data.f90
  sfincs_date.f90
  sfincs_output.f90
  sfincs_ncinput.F90
  sfincs_ncoutput.F90
  sfincs_momentum.f90
  sfincs_continuity.f90
  sfincs_nonhydrostatic.f90
  sfincs_openacc.f90
  sfincs_snapwave.f90
  sfincs_wavemaker.f90
)

set(SFINCS_CORE_SOURCES "")
_append_first_match(SFINCS_CORE_SOURCES ${_CORE_CANDIDATES})

set(SNAPWAVE_DIR "${SFINCS_SOURCE_DIR}/snapwave")
if(EXISTS "${SNAPWAVE_DIR}")
  file(GLOB_RECURSE SNAPWAVE_SOURCES
    "${SNAPWAVE_DIR}/*.f90"
    "${SNAPWAVE_DIR}/*.F90"
  )
  if(SNAPWAVE_SOURCES)
    list(APPEND SFINCS_CORE_SOURCES ${SNAPWAVE_SOURCES})
  endif()
endif()

list(REMOVE_DUPLICATES SFINCS_CORE_SOURCES)

set(SFINCS_LIB_SOURCES "")
_append_first_match(SFINCS_LIB_SOURCES sfincs_lib.f90)

if(SFINCS_CORE_SOURCES STREQUAL "" AND SFINCS_LIB_SOURCES STREQUAL "")
  message(FATAL_ERROR "No SFINCS sources found. Check SFINCS_SOURCE_DIR or file names.")
endif()

# ====================================================================
# TARGET: sfincs_core_objs
# ====================================================================
add_library(sfincs_core_objs OBJECT ${SFINCS_CORE_SOURCES})
set_target_properties(sfincs_core_objs PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
  POSITION_INDEPENDENT_CODE ON
)

if(_F_FLAGS)
  target_compile_options(sfincs_core_objs PRIVATE ${_F_FLAGS})
endif()

target_include_directories(sfincs_core_objs PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
  ${SFINCS_SOURCE_DIR}
  ${SNAPWAVE_DIR}
)

if(SFINCS_ENABLE_NETCDF AND NETCDFF_FOUND)
  target_include_directories(sfincs_core_objs PUBLIC ${NETCDF_INCLUDE_DIR})
  if(NETCDF_MOD_DIR)
    target_include_directories(sfincs_core_objs PUBLIC ${NETCDF_MOD_DIR})
  endif()
endif()

find_package(OpenMP COMPONENTS Fortran)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(sfincs_core_objs PUBLIC OpenMP::OpenMP_Fortran)
endif()

# ====================================================================
# TARGET: sfincs_core
# ====================================================================
add_library(sfincs_core STATIC
  $<TARGET_OBJECTS:sfincs_core_objs>
  ${SFINCS_LIB_SOURCES}
)

set_target_properties(sfincs_core PROPERTIES
  OUTPUT_NAME sfincs_core
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
  POSITION_INDEPENDENT_CODE ON
)

if(_F_FLAGS)
  target_compile_options(sfincs_core PRIVATE ${_F_FLAGS})
endif()

target_include_directories(sfincs_core PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
  ${SFINCS_SOURCE_DIR}
  ${SNAPWAVE_DIR}
)

if(SFINCS_ENABLE_NETCDF AND NETCDFF_FOUND)
  if(NETCDF_C_LIB)
    target_link_libraries(sfincs_core PUBLIC ${NETCDFF_LIB} ${NETCDF_C_LIB})
  else()
    target_link_libraries(sfincs_core PUBLIC ${NETCDFF_LIB})
  endif()
endif()

if(OpenMP_Fortran_FOUND)
  target_link_libraries(sfincs_core PUBLIC OpenMP::OpenMP_Fortran)
endif()

if(SFINCS_ENABLE_NETCDF AND NETCDFF_FOUND AND NETCDF_PREFIX)
  target_link_directories(sfincs_core PUBLIC ${NETCDF_PREFIX}/lib ${NETCDF_PREFIX}/lib64)
  if(UNIX AND NOT APPLE)
    target_link_options(sfincs_core PUBLIC
      -Wl,-rpath,${NETCDF_PREFIX}/lib
      -Wl,-rpath,${NETCDF_PREFIX}/lib64)
  endif()
endif()

# ====================================================================
# TARGET: sfincs_bmi (shared)
# ====================================================================
set(SFINCS_BMI_SOURCES "")
_append_first_match(SFINCS_BMI_SOURCES
  sfincs_bmi2.f90
  sfincs_register_bmi.f90
  sfincs_dll.f90
)

if(NOT SFINCS_BMI_SOURCES STREQUAL "")
  add_library(sfincs_bmi SHARED
    ${SFINCS_BMI_SOURCES}
    ${SFINCS_LIB_SOURCES}
    $<TARGET_OBJECTS:sfincs_core_objs>
  )

  set_target_properties(sfincs_bmi PROPERTIES
    OUTPUT_NAME sfincs_bmi
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_bmi
    POSITION_INDEPENDENT_CODE ON
  )

  if(_F_FLAGS)
    target_compile_options(sfincs_bmi PRIVATE ${_F_FLAGS})
  endif()

  target_include_directories(sfincs_bmi PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_bmi
    ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
    ${SFINCS_SOURCE_DIR}
    ${SNAPWAVE_DIR}
  )

  if(SFINCS_ENABLE_NETCDF AND NETCDFF_FOUND)
    if(NETCDF_C_LIB)
      target_link_libraries(sfincs_bmi PUBLIC ${NETCDFF_LIB} ${NETCDF_C_LIB})
    else()
      target_link_libraries(sfincs_bmi PUBLIC ${NETCDFF_LIB})
    endif()
  endif()

  if(OpenMP_Fortran_FOUND)
    target_link_libraries(sfincs_bmi PUBLIC OpenMP::OpenMP_Fortran)
  endif()

  if(TARGET iso_c_bmi)
    get_target_property(_ISO_MOD_DIR iso_c_bmi Fortran_MODULE_DIRECTORY)
    if(_ISO_MOD_DIR)
      target_include_directories(sfincs_bmi PUBLIC ${_ISO_MOD_DIR})
    endif()
    target_link_libraries(sfincs_bmi PUBLIC iso_c_bmi)
  else()
    message(FATAL_ERROR "Expected target 'iso_c_bmi' from iso_c_fortran_bmi, but it was not defined.")
  endif()

  if(SFINCS_ENABLE_NETCDF AND NETCDFF_FOUND AND NETCDF_PREFIX)
    target_link_directories(sfincs_bmi PUBLIC ${NETCDF_PREFIX}/lib ${NETCDF_PREFIX}/lib64)
    if(UNIX AND NOT APPLE)
      target_link_options(sfincs_bmi PUBLIC
        -Wl,-rpath,${NETCDF_PREFIX}/lib
        -Wl,-rpath,${NETCDF_PREFIX}/lib64)
    endif()
  endif()
endif()

# ---------------- Install (optional) ----------------
include(GNUInstallDirs)

install(TARGETS sfincs_core ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
if(TARGET sfincs_bmi)
  install(TARGETS sfincs_bmi LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} OPTIONAL)
if(TARGET sfincs_bmi)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_bmi DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} OPTIONAL)
endif()

# ====================================================================
# Tests (optional; OFF by default while library build is stabilized)
# ====================================================================
if(SFINCS_BUILD_TESTS)
  include(CTest)
  enable_testing()

  set(_TEST_CANDIDATES
    ${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_sfincs_bmi2.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_sfincs_bmi2.f90
  )

  set(_TEST_SRC "")
  foreach(_c ${_TEST_CANDIDATES})
    if(EXISTS "${_c}")
      set(_TEST_SRC "${_c}")
      break()
    endif()
  endforeach()

  if(_TEST_SRC STREQUAL "")
    message(STATUS "SFINCS: test source not found (tests/test_sfincs_bmi2.f90). Skipping test target.")
  else()
    set(SFINCS_BMI_TEST_WORKDIR
        "/home/mohammed.karim/SFINCS/XMI_EXAMPLE"
        CACHE PATH "Directory containing sfincs.inp for BMI test")

    add_executable(test_sfincs_bmi2 ${_TEST_SRC})

    if(_F_FLAGS)
      target_compile_options(test_sfincs_bmi2 PRIVATE ${_F_FLAGS})
    endif()

    set_target_properties(test_sfincs_bmi2 PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/tests
      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    )

    if(TARGET sfincs_bmi)
      target_link_libraries(test_sfincs_bmi2 PRIVATE sfincs_bmi)
      target_include_directories(test_sfincs_bmi2 PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_bmi
        ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
        ${SFINCS_SOURCE_DIR}
      )
    else()
      message(WARNING "Target sfincs_bmi not found; building test with sfincs_bmi2.f90 directly.")
      target_sources(test_sfincs_bmi2 PRIVATE ${SFINCS_SOURCE_DIR}/sfincs_bmi2.f90)
      target_link_libraries(test_sfincs_bmi2 PRIVATE sfincs_core)
      target_include_directories(test_sfincs_bmi2 PRIVATE
        ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
        ${SFINCS_SOURCE_DIR}
      )
    endif()

    add_test(NAME test_sfincs_bmi2 COMMAND test_sfincs_bmi2)
    set_tests_properties(test_sfincs_bmi2 PROPERTIES
      WORKING_DIRECTORY "${SFINCS_BMI_TEST_WORKDIR}"
    )
  endif()

  set(_TEST_REAL "${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_sfincs_bmi2_real.f90")
  if(EXISTS "${_TEST_REAL}" AND TARGET sfincs_bmi)
    add_executable(test_sfincs_bmi2_real ${_TEST_REAL})
    target_link_libraries(test_sfincs_bmi2_real PRIVATE sfincs_bmi)
    target_include_directories(test_sfincs_bmi2_real PRIVATE
      ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_bmi
      ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
      ${SFINCS_SOURCE_DIR}
    )
    add_test(NAME test_sfincs_bmi2_real COMMAND test_sfincs_bmi2_real)
  endif()
endif()
