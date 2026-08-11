if(NOT DEFINED LCM_SOURCE_DIR)
  message(FATAL_ERROR "LCM_SOURCE_DIR is required")
endif()

set(LCM_MPI_BOUNDARY_SOURCES
  "${LCM_SOURCE_DIR}/lcm_parallel_state.f90"
  "${LCM_SOURCE_DIR}/lcm_particle_count_exchange.f90"
)

file(GLOB LCM_STANDALONE_SOURCES "${LCM_SOURCE_DIR}/lcm_*.f90")

foreach(LCM_SOURCE IN LISTS LCM_STANDALONE_SOURCES)
  list(FIND LCM_MPI_BOUNDARY_SOURCES "${LCM_SOURCE}" LCM_MPI_BOUNDARY_INDEX)
  if(NOT LCM_MPI_BOUNDARY_INDEX EQUAL -1)
    continue()
  endif()

  file(STRINGS "${LCM_SOURCE}" LCM_SOURCE_LINES)
  foreach(LCM_SOURCE_LINE IN LISTS LCM_SOURCE_LINES)
    string(TOUPPER "${LCM_SOURCE_LINE}" LCM_SOURCE_LINE_UPPER)

    if(LCM_SOURCE_LINE_UPPER MATCHES "^[ \t]*INCLUDE[ \t]*['\"]MPIF\\.H['\"]")
      message(FATAL_ERROR
        "Direct mpif.h use is outside the LCM MPI boundary: ${LCM_SOURCE}")
    endif()
    if(LCM_SOURCE_LINE_UPPER MATCHES "^[ \t]*USE[ \t]+MPI(_F08)?([ \t]*,|[ \t]*$)")
      message(FATAL_ERROR
        "Direct MPI module use is outside the LCM MPI boundary: ${LCM_SOURCE}")
    endif()
    if(LCM_SOURCE_LINE_UPPER MATCHES "^[ \t]*CALL[ \t]+MPI_")
      message(FATAL_ERROR
        "Direct MPI call is outside the LCM MPI boundary: ${LCM_SOURCE}")
    endif()
  endforeach()
endforeach()
