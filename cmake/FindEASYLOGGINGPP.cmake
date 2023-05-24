
find_package(PkgConfig QUIET)

pkg_check_modules(EASYLOGGINGPP REQUIRED easyloggingpp)

find_path(ELPP_INCLUDE_DIRS NAMES easylogging++.h)

add_library(
  elpp SHARED
  ${ELPP_INCLUDE_DIRS}/easylogging++.h
  ${ELPP_INCLUDE_DIRS}/easylogging++.cc
)

set(EASYLOGGINGPP_LIBRARIES "elpp")