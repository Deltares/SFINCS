cmake_minimum_required(VERSION 3.19)
project(sfincs LANGUAGES Fortran)

# ---------------- Global PIC for shared libs ----------------
# Ensure all objects are built with -fPIC, so we can link them into libsfincs_bmi.so
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ---------------- Options ----------------
option(SFINCS_ENABLE_NETCDF "Enable NetCDF-Fortran I/O" ON)
option(SFINCS_ENABLE_OPENMP "Enable OpenMP for Fortran" OFF)

# ---------------- Fortran standard ----------------
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
set(CMAKE_Fortran_EXTENSIONS OFF)

# Put executables (if any) in build/bin
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ---------------- Source root ----------------
set(SFINCS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH "Directory with SFINCS Fortran sources")

# Helper: append only files that exist (warn if missing)
function(_append_if_exists outvar)
  foreach(_f IN LISTS ARGN)
    if(EXISTS "${SFINCS_SOURCE_DIR}/${_f}")
      list(APPEND ${outvar} "${SFINCS_SOURCE_DIR}/${_f}")
    else()
      message(WARNING "SFINCS: missing source (skipped): ${SFINCS_SOURCE_DIR}/${_f}")
    endif()
  endforeach()
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
# NetCDF (manual & simple; no generator expressions)
# ====================================================================
set(NETCDF_PREFIX "" CACHE PATH "Prefix for NetCDF (C & Fortran)")
set(NETCDFF_FOUND FALSE)

if(SFINCS_ENABLE_NETCDF)
  if(NOT NETCDF_PREFIX AND DEFINED ENV{NETCDF_PREFIX})
    set(NETCDF_PREFIX "$ENV{NETCDF_PREFIX}")
  endif()

  find_library(NETCDFF_LIB NAMES netcdff
    HINTS ${NETCDF_PREFIX}/lib ${NETCDF_PREFIX}/lib64)
  find_library(NETCDF_C_LIB NAMES netcdf
    HINTS ${NETCDF_PREFIX}/lib ${NETCDF_PREFIX}/lib64)
  find_path(NETCDF_INCLUDE_DIR NAMES netcdf.h netcdf.inc
    HINTS ${NETCDF_PREFIX}/include)
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
      "Set -DNETCDF_PREFIX=/opt/netcdfc_4_7_4_netcdff_4_5_4 or build with -DSFINCS_ENABLE_NETCDF=OFF.")
  endif()
endif()

# ====================================================================
# CORE SOURCES (exact list) + add all snapwave/*.f90,*.F90
# ====================================================================
set(_CORE_CANDIDATES
  bmi.f90
  deg2utm.f90
  sfincs_spiderweb.f90
  sfincs_error.f90
  sfincs_input.f90
  sfincs_domain.f90
  sfincs_structures.f90
  sfincs_boundaries.f90
  sfincs_obspoints.f90
  sfincs_crosssections.f90
  sfincs_discharges.f90
  sfincs_meteo.f90
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
  sfincs_snapwave.f90
  sfincs_wavemaker.f90
)
set(SFINCS_CORE_SOURCES "")
_append_if_exists(SFINCS_CORE_SOURCES ${_CORE_CANDIDATES})

# include all snapwave sources in the same build (module order auto)
set(SNAPWAVE_DIR "${SFINCS_SOURCE_DIR}/snapwave")
if(EXISTS "${SNAPWAVE_DIR}")
  file(GLOB SNAPWAVE_SOURCES
    "${SNAPWAVE_DIR}/*.f90"
    "${SNAPWAVE_DIR}/*.F90")
  if(SNAPWAVE_SOURCES)
    list(APPEND SFINCS_CORE_SOURCES ${SNAPWAVE_SOURCES})
    message(STATUS "SFINCS: added snapwave sources: ${SNAPWAVE_SOURCES}")
  else()
    message(STATUS "SFINCS: snapwave dir present but empty: ${SNAPWAVE_DIR}")
  endif()
else()
  message(STATUS "SFINCS: snapwave directory not found: ${SNAPWAVE_DIR}")
endif()

# ---------------- Glue (sfincs_lib) ----------------
# Keep sfincs_lib.f90 separate so we can compile it *after* objects
set(SFINCS_LIB_SOURCES "")
_append_if_exists(SFINCS_LIB_SOURCES sfincs_lib.f90)

# ---------------- Sanity check ----------------
if(SFINCS_CORE_SOURCES STREQUAL "" AND SFINCS_LIB_SOURCES STREQUAL "")
  message(FATAL_ERROR "No SFINCS sources found. Check SFINCS_SOURCE_DIR or file names.")
endif()

# ====================================================================
# TARGET 1: OBJECT library for all module producers (except sfincs_lib.f90)
# ====================================================================
add_library(sfincs_core_objs OBJECT ${SFINCS_CORE_SOURCES})
set_target_properties(sfincs_core_objs PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
  POSITION_INDEPENDENT_CODE ON          # ensure OBJECT files are built with -fPIC
)
if(_F_FLAGS)
  target_compile_options(sfincs_core_objs PRIVATE ${_F_FLAGS})
endif()
target_include_directories(sfincs_core_objs
  PUBLIC
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
if(SFINCS_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
  target_link_libraries(sfincs_core_objs PUBLIC OpenMP::OpenMP_Fortran)
endif()

# ====================================================================
# TARGET 2: final static lib that adds sfincs_lib.f90 on top of OBJECTS
# ====================================================================
add_library(sfincs_core STATIC
  $<TARGET_OBJECTS:sfincs_core_objs>
  ${SFINCS_LIB_SOURCES}
)
set_target_properties(sfincs_core PROPERTIES
  OUTPUT_NAME sfincs_core
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
  POSITION_INDEPENDENT_CODE ON       # needed so we can link into shared sfincs_bmi
)
if(_F_FLAGS)
  target_compile_options(sfincs_core PRIVATE ${_F_FLAGS})
endif()
target_include_directories(sfincs_core
  PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
    ${SFINCS_SOURCE_DIR}
    ${SNAPWAVE_DIR}
)

# Link NetCDF to the final lib (and set rpath/search paths)
if(SFINCS_ENABLE_NETCDF AND NETCDFF_FOUND)
  if(NETCDF_C_LIB)
    target_link_libraries(sfincs_core PUBLIC ${NETCDFF_LIB} ${NETCDF_C_LIB})
  else()
    target_link_libraries(sfincs_core PUBLIC ${NETCDFF_LIB})
  endif()
  if(NETCDF_PREFIX)
    target_link_directories(sfincs_core PUBLIC
      ${NETCDF_PREFIX}/lib
      ${NETCDF_PREFIX}/lib64)
    if(UNIX AND NOT APPLE)
      target_link_options(sfincs_core PUBLIC
        -Wl,-rpath,${NETCDF_PREFIX}/lib
        -Wl,-rpath,${NETCDF_PREFIX}/lib64)
    endif()
  endif()
endif()

# ====================================================================
# OPTIONAL: BMI wrapper (if sfincs_bmi2.f90 exists)
# ====================================================================
# OPTIONAL: BMI wrapper (if sfincs_bmi2.f90 exists)
set(SFINCS_BMI_SOURCES "")
_append_if_exists(SFINCS_BMI_SOURCES sfincs_bmi2.f90)
# NEW: also add the XMI shim if present
_append_if_exists(SFINCS_BMI_SOURCES sfincs_dll.f90)

if(NOT SFINCS_BMI_SOURCES STREQUAL "")
  add_library(sfincs_bmi SHARED ${SFINCS_BMI_SOURCES})
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
  target_link_libraries(sfincs_bmi PUBLIC sfincs_core)

  if(SFINCS_ENABLE_NETCDF AND NETCDFF_FOUND AND NETCDF_PREFIX)
    target_link_directories(sfincs_bmi PUBLIC
      ${NETCDF_PREFIX}/lib
      ${NETCDF_PREFIX}/lib64)
    if(UNIX AND NOT APPLE)
      target_link_options(sfincs_bmi PUBLIC
        -Wl,-rpath,${NETCDF_PREFIX}/lib
        -Wl,-rpath,${NETCDF_PREFIX}/lib64)
    endif()
  endif()
endif()

# ---------------- Install (optional) ----------------
include(GNUInstallDirs)
install(TARGETS sfincs_core
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
if(TARGET sfincs_bmi)
  install(TARGETS sfincs_bmi
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
endif()
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} OPTIONAL)
if(TARGET sfincs_bmi)
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_bmi
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} OPTIONAL)
endif()

# ====================================================================
# ================================ Tests ==============================
# ====================================================================
include(CTest)
enable_testing()

# Try to locate the test source in two common places:
#   1) <repo_root>/tests/test_sfincs_bmi2.f90        (recommended)
#   2) <this dir>/tests/test_sfincs_bmi2.f90         (fallback)
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
  # Directory that contains a real sfincs.inp for the BMI test.
  # Default is your XMI example, but can be overridden at configure time:
  #   -DSFINCS_BMI_TEST_WORKDIR=/path/to/another/example
  set(SFINCS_BMI_TEST_WORKDIR
      "/home/mohammed.karim/SFINCS/XMI_EXAMPLE"
      CACHE PATH "Directory containing sfincs.inp for BMI test")

  add_executable(test_sfincs_bmi2 ${_TEST_SRC})

  # Use same warning flags for the test
  if(_F_FLAGS)
    target_compile_options(test_sfincs_bmi2 PRIVATE ${_F_FLAGS})
  endif()

  # Keep test .mod files separate
  set_target_properties(test_sfincs_bmi2 PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/tests
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
  )

  # Prefer to link the existing BMI lib if present (Option A)
  if(TARGET sfincs_bmi)
    target_link_libraries(test_sfincs_bmi2 PRIVATE sfincs_bmi)
    target_include_directories(test_sfincs_bmi2 PRIVATE
      ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_bmi
      ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
      ${SFINCS_SOURCE_DIR}
    )
  else()
    message(WARNING "Target sfincs_bmi not found; building test with sfincs_bmi2.f90 directly (Option B).")
    # Compile the BMI wrapper directly into the test; link the core library.
    target_sources(test_sfincs_bmi2 PRIVATE ${SFINCS_SOURCE_DIR}/sfincs_bmi2.f90)
    target_link_libraries(test_sfincs_bmi2 PRIVATE sfincs_core)
    target_include_directories(test_sfincs_bmi2 PRIVATE
      ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
      ${SFINCS_SOURCE_DIR}
    )
  endif()

  # Register with CTest and point it at the real SFINCS example dir
  add_test(NAME test_sfincs_bmi2 COMMAND test_sfincs_bmi2)
  set_tests_properties(test_sfincs_bmi2 PROPERTIES
    WORKING_DIRECTORY "${SFINCS_BMI_TEST_WORKDIR}"
  )
endif()

# ============================================================
# Build XMI shared library: libsfincs_dll.so
# ============================================================
# OPTIONAL: BMI wrapper (if sfincs_bmi2.f90 exists)


set(_TEST_REAL "${CMAKE_CURRENT_SOURCE_DIR}/../tests/test_sfincs_bmi2_real.f90")
if(EXISTS "${_TEST_REAL}")
  add_executable(test_sfincs_bmi2_real ${_TEST_REAL})
  target_link_libraries(test_sfincs_bmi2_real PRIVATE sfincs_bmi)
  target_include_directories(test_sfincs_bmi2_real PRIVATE
    ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_bmi
    ${CMAKE_CURRENT_BINARY_DIR}/mod/sfincs_core
    ${SFINCS_SOURCE_DIR}
  )
  add_test(NAME test_sfincs_bmi2_real COMMAND test_sfincs_bmi2_real)
endif()

