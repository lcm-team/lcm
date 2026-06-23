cmake_minimum_required(VERSION 3.20)

project(lcm LANGUAGES Fortran)

option(LCM_BUILD_COMPAT_STUBS "Build LCM with temporary SAM compatibility stubs" ON)
option(LCM_BUILD_ALL_SAM_SOURCES "Attempt to build all SAM-derived LCM source files" OFF)

find_package(MPI REQUIRED COMPONENTS Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")

add_library(lcm_sam_compat STATIC)

if(LCM_BUILD_COMPAT_STUBS)
  target_sources(lcm_sam_compat
    PRIVATE
      stubs/sam_compat.f90
  )
else()
  message(FATAL_ERROR "Standalone builds currently require LCM_BUILD_COMPAT_STUBS=ON")
endif()

set(LCM_CORE_SOURCES
  micro_user_vars.f90
  micro_vars.f90
  micro_task_util_MPI.f90
  micro_random.f90
  micro_restart.f90
  microphysics.f90
)

set(LCM_SAM_SOURCES
  micro_ds.f90
  micro_boundaries.f90
  micro_coll.f90
  micro_cond.f90
  micro_coupling.f90
  micro_exchange.f90
  micro_init_particles.f90
  micro_sedi.f90
  micro_set.f90
  micro_sgs_mixing.f90
  micro_sgs_uvw.f90
  micro_sort.f90
  micro_user.f90
  micro_uvw.f90
  micro_xyz.f90
)

set(LCM_SOURCES ${LCM_CORE_SOURCES})

if(LCM_BUILD_ALL_SAM_SOURCES)
  list(APPEND LCM_SOURCES ${LCM_SAM_SOURCES})
endif()

add_library(lcm STATIC ${LCM_SOURCES})

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(lcm PRIVATE -cpp -ffree-line-length-none -fallow-argument-mismatch)
  target_compile_options(lcm_sam_compat PRIVATE -cpp -ffree-line-length-none)
endif()

target_link_libraries(lcm
  PRIVATE
    lcm_sam_compat
    MPI::MPI_Fortran
)

target_include_directories(lcm
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>"
)

set_target_properties(lcm lcm_sam_compat
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)
