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

# Check for MPI installation
find_package(MPI COMPONENTS C CXX REQUIRED)
