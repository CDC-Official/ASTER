# ASTER integration only; upstream Eigen sources remain unmodified.
# Use an INTERFACE target without enabling upstream tests, docs or installation.
set(ASTER_EIGEN_SOURCE_DIR "${PROJECT_SOURCE_DIR}/third_party/eigen" CACHE PATH
    "Eigen source root; default is bundled 3.4.0. Use another build directory for version trials.")
set(_aster_eigen_dir "${ASTER_EIGEN_SOURCE_DIR}")
if(NOT EXISTS "${_aster_eigen_dir}/Eigen/Dense")
    message(FATAL_ERROR "Missing Eigen/Dense under ASTER_EIGEN_SOURCE_DIR. Restore the bundled source or specify a valid trial source root.")
endif()
add_library(aster_eigen INTERFACE)
add_library(Eigen3::Eigen ALIAS aster_eigen)
target_include_directories(aster_eigen SYSTEM INTERFACE "${_aster_eigen_dir}")
message(STATUS "ASTER Eigen source: ${ASTER_EIGEN_SOURCE_DIR}")

file(STRINGS "${_aster_eigen_dir}/Eigen/src/Core/util/Macros.h" _aster_eigen_version
    REGEX "^#define EIGEN_(WORLD|MAJOR|MINOR)_VERSION ")
if(NOT _aster_eigen_version STREQUAL "#define EIGEN_WORLD_VERSION 3;#define EIGEN_MAJOR_VERSION 4;#define EIGEN_MINOR_VERSION 0")
    message(FATAL_ERROR "ASTER requires Eigen 3.4.0")
endif()
