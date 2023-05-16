
# Check for MPI installation.
find_package(MPI COMPONENTS C CXX REQUIRED)

include_directories(${MPI_CXX_INCLUDE_DIRS})