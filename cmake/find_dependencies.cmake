
# Check for MPI installation.
find_package(MPI COMPONENTS C CXX REQUIRED)

# Check for Python3
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

include_directories(${MPI_CXX_INCLUDE_DIRS}
                    ${Python3_INCLUDE_DIRS}
                   )