add_executable(benchmark_particle_exchange
  benchmark_particle_exchange.f90
)

target_link_libraries(benchmark_particle_exchange
  PRIVATE
    lcm
    MPI::MPI_Fortran
)

set_target_properties(benchmark_particle_exchange
  PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(benchmark_particle_exchange
    PRIVATE
      -ffree-line-length-none
  )
endif()
