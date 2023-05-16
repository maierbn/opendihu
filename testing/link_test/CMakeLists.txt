
# Test for MPI installation.
add_executable(test_mpi test_mpi.cpp)
target_link_libraries(test_mpi ${MPI_CXX_LIBRARIES})
add_test(NAME test_mpi COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 test_mpi)