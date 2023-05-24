
# Check for MPI installation.
find_package(MPI COMPONENTS C CXX REQUIRED)

# Check for zlib installation
find_package(ZLIB REQUIRED)

# Check for Python3
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

# Check for Eigen3
find_package(Eigen3 REQUIRED NO_MODULE)
if(${Eigen3_FOUND})
  message(STATUS "Found Eigen3: yes")
endif()

# Check for libboost
find_package(Boost REQUIRED COMPONENTS log log_setup thread system filesystem program_options unit_test_framework)

# Check for libxml2
find_package(LibXml2 REQUIRED)

# Check for easylogging++
include(cmake/FindEASYLOGGINGPP.cmake)

include_directories(${MPI_CXX_INCLUDE_DIRS}
                    ${Python3_INCLUDE_DIRS}
                    ${LIBXML2_INCLUDE_DIR}
                    ${ZLIB_INCLUDE_DIRS}
                   )