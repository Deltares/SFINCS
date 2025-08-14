cmake_minimum_required(VERSION 3.10)
project(sfincs_bmi_project Fortran)

set(CMAKE_Fortran_STANDARD 2003)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)

# Enable colored output
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  add_compile_options(-fdiagnostics-color=always)
endif()

# Add source files
set(SFincsSources
    bmi.f90
    sfincs_lib2.f90
    sfincs_bmi2.f90
    sfincs_bmi_driver_test.f90
)

# Create executable (for testing)
add_executable(sfincs_bmi_driver ${SFincsSources})

# Create library (for linking in other applications)
add_library(sfincs_bmi STATIC
    bmi.f90
    sfincs_lib2.f90
    sfincs_bmi2.f90
    #sfincs_model_container.f90
)

# Include paths if needed
target_include_directories(sfincs_bmi_driver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(sfincs_bmi PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Link executable with static library
target_link_libraries(sfincs_bmi_driver PRIVATE sfincs_bmi)

