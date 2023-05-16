cmake_minimum_required(VERSION 3.13...3.24)

# package versioning
include(cmake/git.cmake)

# Print package version to console.
message(STATUS "OpenDiHu ${PROJECT_VERSION}")

# Project description
project(OPENDIHU
LANGUAGES C CXX
HOMEPAGE_URL https://opendihu.readthedocs.io/en/latest/introduction.html#overview
DESCRIPTION "OpenDiHu is a scalable Finite Element code for skeletal muscle simulations."
VERSION ${PROJECT_VERSION}
)

# Find the packages of all dependencies.
# The configure is aborted if a package cannot be found.
include(cmake/find_dependencies.cmake)

# Unit testing for dependecies and source code.
include(CTest)
add_subdirectory(testing/link_test)
