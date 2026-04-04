function(gwm_add_mpi_test target num_ranks)
  add_executable(${target} ${ARGN})
  target_link_libraries(${target} PRIVATE gwm_core CUDA::cudart MPI::MPI_CXX)
  target_include_directories(${target} PRIVATE tests)
  if(WIN32)
    add_test(
      NAME ${target}
      COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} ${num_ranks}
        $<TARGET_FILE:${target}>
    )
  else()
    add_test(
      NAME ${target}
      COMMAND
        ${MPIEXEC_EXECUTABLE}
        ${MPIEXEC_NUMPROC_FLAG} ${num_ranks}
        ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:${target}>
        ${MPIEXEC_POSTFLAGS}
    )
  endif()
endfunction()
