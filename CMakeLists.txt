cmake_minimum_required(VERSION 3.20)

project(lcm LANGUAGES Fortran)

option(LCM_BUILD_COMPAT_STUBS "Build LCM with temporary SAM compatibility stubs" ON)
option(LCM_BUILD_ALL_SAM_SOURCES "Attempt to build all SAM-derived LCM source files" OFF)

enable_testing()

find_package(MPI REQUIRED COMPONENTS Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod")
set(LCM_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

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
  "${LCM_SOURCE_DIR}/lcm_parameters.f90"
  "${LCM_SOURCE_DIR}/lcm_types.f90"
  "${LCM_SOURCE_DIR}/lcm_grid_state.f90"
  "${LCM_SOURCE_DIR}/lcm_field_state.f90"
  "${LCM_SOURCE_DIR}/lcm_random.f90"
  "${LCM_SOURCE_DIR}/lcm_particle_state.f90"
  "${LCM_SOURCE_DIR}/micro_user_vars.f90"
  "${LCM_SOURCE_DIR}/micro_vars.f90"
  "${LCM_SOURCE_DIR}/lcm_config_state.f90"
  "${LCM_SOURCE_DIR}/lcm_initialization_state.f90"
  "${LCM_SOURCE_DIR}/lcm_advance.f90"
  "${LCM_SOURCE_DIR}/lcm_host_interface.f90"
  "${LCM_SOURCE_DIR}/micro_task_util_MPI.f90"
)

set(LCM_SAM_SOURCES
  "${LCM_SOURCE_DIR}/micro_random.f90"
  "${LCM_SOURCE_DIR}/micro_restart.f90"
  "${LCM_SOURCE_DIR}/microphysics.f90"
  "${LCM_SOURCE_DIR}/micro_ds.f90"
  "${LCM_SOURCE_DIR}/micro_boundaries.f90"
  "${LCM_SOURCE_DIR}/micro_coll.f90"
  "${LCM_SOURCE_DIR}/micro_cond.f90"
  "${LCM_SOURCE_DIR}/micro_coupling.f90"
  "${LCM_SOURCE_DIR}/micro_exchange.f90"
  "${LCM_SOURCE_DIR}/micro_init_particles.f90"
  "${LCM_SOURCE_DIR}/micro_sedi.f90"
  "${LCM_SOURCE_DIR}/micro_set.f90"
  "${LCM_SOURCE_DIR}/micro_sgs_mixing.f90"
  "${LCM_SOURCE_DIR}/micro_sgs_uvw.f90"
  "${LCM_SOURCE_DIR}/micro_sort.f90"
  "${LCM_SOURCE_DIR}/micro_user.f90"
  "${LCM_SOURCE_DIR}/micro_uvw.f90"
  "${LCM_SOURCE_DIR}/micro_xyz.f90"
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

add_executable(test_initial_particle_population
  tests/test_initial_particle_population.f90
)

target_link_libraries(test_initial_particle_population
  PRIVATE
    lcm
)

set_target_properties(test_initial_particle_population
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME initial_particle_population_default_aerosol
  COMMAND test_initial_particle_population
          default_aerosol
)

add_test(
  NAME initial_particle_population_one_mode
  COMMAND test_initial_particle_population
          one_mode
)

add_test(
  NAME initial_particle_population_multi_mode
  COMMAND test_initial_particle_population
          multi_mode
)

add_test(
  NAME initial_particle_population_zero_weight_split
  COMMAND test_initial_particle_population
          zero_weight_split
)

add_test(
  NAME initial_particle_population_usm_aerosol
  COMMAND test_initial_particle_population
          usm_aerosol
)

add_test(
  NAME initial_particle_population_random_aerosol
  COMMAND test_initial_particle_population
          random_aerosol
)

set(SAM_MICRO_INIT_PROBE_MODULE_DIR "${CMAKE_CURRENT_BINARY_DIR}/sam_probe_mod")

add_library(sam_micro_init_particles_probe_objects OBJECT
  stubs/sam_compat.f90
  src/micro_user_vars.f90
  src/micro_vars.f90
  src/micro_random.f90
  tests/sam_micro_init_particles_probe_stubs.f90
  src/micro_init_particles.f90
)

set_target_properties(sam_micro_init_particles_probe_objects
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${SAM_MICRO_INIT_PROBE_MODULE_DIR}"
)

target_include_directories(sam_micro_init_particles_probe_objects
  PRIVATE
    "${SAM_MICRO_INIT_PROBE_MODULE_DIR}"
)

add_executable(sam_micro_init_particles_probe
  tests/sam_micro_init_particles_probe.f90
  $<TARGET_OBJECTS:sam_micro_init_particles_probe_objects>
)

set_target_properties(sam_micro_init_particles_probe
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${SAM_MICRO_INIT_PROBE_MODULE_DIR}"
)

target_include_directories(sam_micro_init_particles_probe
  PRIVATE
    "${SAM_MICRO_INIT_PROBE_MODULE_DIR}"
)

add_executable(test_sam_particle_population_comparison
  tests/test_sam_particle_population_comparison.f90
)

target_link_libraries(test_sam_particle_population_comparison
  PRIVATE
    lcm
)

set_target_properties(test_sam_particle_population_comparison
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_dependencies(test_sam_particle_population_comparison
  sam_micro_init_particles_probe
)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(sam_micro_init_particles_probe_objects
    PRIVATE
      -cpp
      -ffree-line-length-none
      -fallow-argument-mismatch
      -fdefault-real-8
  )
  target_compile_options(sam_micro_init_particles_probe
    PRIVATE
      -cpp
      -ffree-line-length-none
      -fallow-argument-mismatch
      -fdefault-real-8
  )
endif()

add_test(
  NAME initial_particle_population_sam_bin_comparison
  COMMAND test_sam_particle_population_comparison
          $<TARGET_FILE:sam_micro_init_particles_probe>
          bin
)

add_test(
  NAME initial_particle_population_sam_usm_comparison
  COMMAND test_sam_particle_population_comparison
          $<TARGET_FILE:sam_micro_init_particles_probe>
          USM
)
