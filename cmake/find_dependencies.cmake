
# Check for MPI installation.
find_package(MPI COMPONENTS C CXX REQUIRED)

# Check for Python3
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

# Check for Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)
if(${Eigen3_FOUND})
  message(STATUS "Found Eigen3: yes")
endif()

include_directories(${MPI_CXX_INCLUDE_DIRS}
                    ${Python3_INCLUDE_DIRS}
                   )