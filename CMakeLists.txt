cmake_minimum_required(VERSION 3.20)

project(lcm LANGUAGES Fortran)

option(LCM_BUILD_COMPAT_STUBS "Build LCM with temporary SAM compatibility stubs" ON)
option(LCM_BUILD_ALL_SAM_SOURCES "Attempt to build all SAM-derived LCM source files" OFF)
option(LCM_BUILD_EXPERIMENTAL_BENCHMARKS
  "Build opt-in experimental LCM performance benchmarks"
  OFF
)

enable_testing()

add_test(
  NAME mpi_interface_isolation
  COMMAND ${CMAKE_COMMAND}
          -DLCM_SOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/src
          -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/cmake/check_mpi_interface_isolation.cmake
)

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
  "${LCM_SOURCE_DIR}/lcm_particle_mpi_messages.f90"
  "${LCM_SOURCE_DIR}/lcm_particle_mpi_datatype.f90"
  "${LCM_SOURCE_DIR}/lcm_parallel_state.f90"
  "${LCM_SOURCE_DIR}/lcm_boundary_diagnostics.f90"
  "${LCM_SOURCE_DIR}/lcm_particle_properties.f90"
  "${LCM_SOURCE_DIR}/lcm_grid_state.f90"
  "${LCM_SOURCE_DIR}/lcm_exchange_directions.f90"
  "${LCM_SOURCE_DIR}/lcm_particle_exchange_buffers.f90"
  "${LCM_SOURCE_DIR}/lcm_particle_count_exchange.f90"
  "${LCM_SOURCE_DIR}/lcm_particle_mpi_transport.f90"
  "${LCM_SOURCE_DIR}/lcm_field_state.f90"
  "${LCM_SOURCE_DIR}/lcm_random.f90"
  "${LCM_SOURCE_DIR}/lcm_particle_location.f90"
  "${LCM_SOURCE_DIR}/lcm_velocity_interpolation.f90"
  "${LCM_SOURCE_DIR}/lcm_user.f90"
  "${LCM_SOURCE_DIR}/lcm_particle_state.f90"
  "${LCM_SOURCE_DIR}/lcm_particle_exchange.f90"
  "${LCM_SOURCE_DIR}/lcm_particle_boundaries.f90"
  "${LCM_SOURCE_DIR}/lcm_particle_source.f90"
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
  PUBLIC
    MPI::MPI_Fortran
  PRIVATE
    lcm_sam_compat
)

target_include_directories(lcm
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>"
)

set_target_properties(lcm lcm_sam_compat
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

if(LCM_BUILD_EXPERIMENTAL_BENCHMARKS)
  add_subdirectory(benchmarks/experimental/particle_exchange)
endif()

add_executable(test_initial_particle_population
  tests/integration/test_initial_particle_population.f90
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
  NAME initial_particle_population_tracking
  COMMAND test_initial_particle_population
          tracking
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

add_test(
  NAME initial_particle_population_bin_ice
  COMMAND test_initial_particle_population
          bin_ice
)

add_executable(test_minimal_les_host_initialization
  tests/integration/minimal_les_host.f90
  tests/integration/minimal_les_host_references.f90
  tests/integration/test_minimal_les_host_initialization.f90
)

target_link_libraries(test_minimal_les_host_initialization
  PRIVATE
    lcm
)

set_target_properties(test_minimal_les_host_initialization
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME minimal_les_host_initialization_bin
  COMMAND test_minimal_les_host_initialization
          bin
)

add_test(
  NAME minimal_les_host_initialization_usm
  COMMAND test_minimal_les_host_initialization
          USM
)

add_test(
  NAME minimal_les_host_initialization_random
  COMMAND test_minimal_les_host_initialization
          random
)

add_test(
  NAME minimal_les_host_initialization_bin_ice
  COMMAND test_minimal_les_host_initialization
          bin_ice
)

add_executable(test_particle_location
  tests/integration/minimal_les_host.f90
  tests/unit/test_particle_location.f90
)

target_link_libraries(test_particle_location
  PRIVATE
    lcm
)

set_target_properties(test_particle_location
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME particle_location
  COMMAND test_particle_location
)

add_executable(test_particle_cell_ownership
  tests/unit/test_particle_cell_ownership.f90
)

target_link_libraries(test_particle_cell_ownership
  PRIVATE
    lcm
)

set_target_properties(test_particle_cell_ownership
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME particle_cell_ownership
  COMMAND test_particle_cell_ownership
)

add_executable(test_parallel_state
  tests/unit/test_parallel_state.f90
)

target_link_libraries(test_parallel_state
  PRIVATE
    lcm
    MPI::MPI_Fortran
)

set_target_properties(test_parallel_state
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME parallel_state
  COMMAND test_parallel_state
)

add_executable(test_particle_message
  tests/unit/test_particle_message.f90
)

target_link_libraries(test_particle_message
  PRIVATE
    lcm
)

set_target_properties(test_particle_message
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME particle_message
  COMMAND test_particle_message
)

add_executable(test_particle_mpi_datatype
  tests/unit/test_particle_mpi_datatype.f90
)

target_link_libraries(test_particle_mpi_datatype
  PRIVATE
    lcm
    MPI::MPI_Fortran
)

set_target_properties(test_particle_mpi_datatype
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME particle_mpi_datatype_mpi_2
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
          ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_particle_mpi_datatype>
          ${MPIEXEC_POSTFLAGS}
)

add_executable(test_exchange_directions
  tests/unit/test_exchange_directions.f90
)

target_link_libraries(test_exchange_directions
  PRIVATE
    lcm
)

set_target_properties(test_exchange_directions
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME exchange_directions
  COMMAND test_exchange_directions
)

add_executable(test_local_subdomain_bounds
  tests/unit/test_local_subdomain_bounds.f90
)

target_link_libraries(test_local_subdomain_bounds
  PRIVATE
    lcm
)

set_target_properties(test_local_subdomain_bounds
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME local_subdomain_bounds
  COMMAND test_local_subdomain_bounds
)

add_executable(test_particle_exchange_buffers
  tests/unit/test_particle_exchange_buffers.f90
)

target_link_libraries(test_particle_exchange_buffers
  PRIVATE
    lcm
)

set_target_properties(test_particle_exchange_buffers
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME particle_exchange_buffers
  COMMAND test_particle_exchange_buffers
)

add_executable(test_outgoing_particle_message_packing
  tests/unit/test_outgoing_particle_message_packing.f90
)

target_link_libraries(test_outgoing_particle_message_packing
  PRIVATE
    lcm
)

set_target_properties(test_outgoing_particle_message_packing
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME outgoing_particle_message_packing
  COMMAND test_outgoing_particle_message_packing
)

add_executable(test_periodic_particle_destination_wrapping
  tests/unit/test_periodic_particle_destination_wrapping.f90
)

target_link_libraries(test_periodic_particle_destination_wrapping
  PRIVATE
    lcm
)

set_target_properties(test_periodic_particle_destination_wrapping
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME periodic_particle_destination_wrapping
  COMMAND test_periodic_particle_destination_wrapping
)

add_executable(test_particle_count_exchange
  tests/unit/test_particle_count_exchange.f90
)

target_link_libraries(test_particle_count_exchange
  PRIVATE
    lcm
    MPI::MPI_Fortran
)

set_target_properties(test_particle_count_exchange
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME particle_count_exchange
  COMMAND test_particle_count_exchange
)

add_test(
  NAME particle_count_exchange_mpi_2
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
          ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_particle_count_exchange>
          ${MPIEXEC_POSTFLAGS}
)

add_test(
  NAME particle_count_exchange_mpi_2_nonperiodic
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
          ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_particle_count_exchange>
          ${MPIEXEC_POSTFLAGS} nonperiodic
)

add_test(
  NAME particle_count_exchange_mpi_4_cartesian
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
          ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_particle_count_exchange>
          ${MPIEXEC_POSTFLAGS} square
)

add_executable(test_particle_exchange_mpi
  tests/integration/test_particle_exchange_mpi.f90
)

target_link_libraries(test_particle_exchange_mpi
  PRIVATE
    lcm
    MPI::MPI_Fortran
)

set_target_properties(test_particle_exchange_mpi
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME particle_exchange_mpi_4_corner_routes
  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
          ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_particle_exchange_mpi>
          ${MPIEXEC_POSTFLAGS}
)

add_executable(test_minimal_les_host_cell_crossing
  tests/integration/minimal_les_host.f90
  tests/integration/test_minimal_les_host_cell_crossing.f90
)

target_link_libraries(test_minimal_les_host_cell_crossing
  PRIVATE
    lcm
)

set_target_properties(test_minimal_les_host_cell_crossing
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME minimal_les_host_cell_crossing_x
  COMMAND test_minimal_les_host_cell_crossing
          x
)

add_test(
  NAME minimal_les_host_cell_crossing_y
  COMMAND test_minimal_les_host_cell_crossing
          y
)

add_test(
  NAME minimal_les_host_cell_crossing_z
  COMMAND test_minimal_les_host_cell_crossing
          z
)

add_executable(test_minimal_les_host_outgoing_particles
  tests/integration/minimal_les_host.f90
  tests/integration/test_minimal_les_host_outgoing_particles.f90
)

target_link_libraries(test_minimal_les_host_outgoing_particles
  PRIVATE
    lcm
)

set_target_properties(test_minimal_les_host_outgoing_particles
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME minimal_les_host_serial_particle_exchange
  COMMAND test_minimal_les_host_outgoing_particles
)

add_executable(test_field_attachment
  tests/unit/test_field_attachment.f90
)

target_link_libraries(test_field_attachment
  PRIVATE
    lcm
)

set_target_properties(test_field_attachment
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME field_attachment
  COMMAND test_field_attachment
)

add_executable(test_velocity_interpolation
  tests/unit/test_velocity_interpolation.f90
)

target_link_libraries(test_velocity_interpolation
  PRIVATE
    lcm
)

set_target_properties(test_velocity_interpolation
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME velocity_interpolation
  COMMAND test_velocity_interpolation
)

add_executable(test_particle_radius
  tests/unit/test_particle_radius.f90
)

target_link_libraries(test_particle_radius
  PRIVATE
    lcm
)

set_target_properties(test_particle_radius
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME particle_radius
  COMMAND test_particle_radius
)

add_executable(test_boundary_config
  tests/unit/test_boundary_config.f90
)

target_link_libraries(test_boundary_config
  PRIVATE
    lcm
)

set_target_properties(test_boundary_config
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME boundary_config
  COMMAND test_boundary_config
)

add_executable(test_vertical_boundary_bounds
  tests/unit/test_vertical_boundary_bounds.f90
)

target_link_libraries(test_vertical_boundary_bounds
  PRIVATE
    lcm
)

set_target_properties(test_vertical_boundary_bounds
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME vertical_boundary_bounds
  COMMAND test_vertical_boundary_bounds
)

add_executable(test_particle_cell_removal
  tests/unit/test_particle_cell_removal.f90
)

target_link_libraries(test_particle_cell_removal
  PRIVATE
    lcm
)

set_target_properties(test_particle_cell_removal
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME particle_cell_removal
  COMMAND test_particle_cell_removal
)

add_executable(test_minimal_les_host_top_particle_boundary
  tests/integration/minimal_les_host.f90
  tests/integration/test_minimal_les_host_top_particle_boundary.f90
)

target_link_libraries(test_minimal_les_host_top_particle_boundary
  PRIVATE
    lcm
)

set_target_properties(test_minimal_les_host_top_particle_boundary
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME minimal_les_host_top_particle_boundary
  COMMAND test_minimal_les_host_top_particle_boundary
)

add_executable(test_minimal_les_host_bottom_particle_boundary
  tests/integration/minimal_les_host.f90
  tests/integration/test_minimal_les_host_bottom_particle_boundary.f90
)

target_link_libraries(test_minimal_les_host_bottom_particle_boundary
  PRIVATE
    lcm
)

set_target_properties(test_minimal_les_host_bottom_particle_boundary
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME minimal_les_host_bottom_particle_boundary
  COMMAND test_minimal_les_host_bottom_particle_boundary
)

add_executable(test_particle_source
  tests/integration/test_particle_source.f90
)

target_link_libraries(test_particle_source
  PRIVATE
    lcm
)

set_target_properties(test_particle_source
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

add_test(
  NAME particle_source_aerosol
  COMMAND test_particle_source
          aerosol
)

add_test(
  NAME particle_source_in
  COMMAND test_particle_source
          in
)

add_test(
  NAME particle_source_moving_aerosol
  COMMAND test_particle_source
          moving_aerosol
)

add_test(
  NAME particle_source_moving_in
  COMMAND test_particle_source
          moving_in
)
