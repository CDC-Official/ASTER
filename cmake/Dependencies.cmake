# Dependencies belong to targets, not global compiler flags.
# Pinned header-only dependency, kept verbatim under third_party/eigen.
include("${CMAKE_CURRENT_LIST_DIR}/Eigen.cmake")
find_path(ASTER_BOOST_INCLUDE_DIR NAMES boost/math/distributions/chi_squared.hpp
    HINTS ${ASTER_BOOST_ROOT} /usr/local/include/boost_1_87_0
    PATH_SUFFIXES include REQUIRED)

include("${CMAKE_CURRENT_LIST_DIR}/Sofa.cmake")

# Linux compilers normally discover OpenMP without hints.
find_package(OpenMP QUIET COMPONENTS CXX)
if(NOT OpenMP_CXX_FOUND AND APPLE)
    find_program(ASTER_BREW_EXECUTABLE brew)
    if(ASTER_BREW_EXECUTABLE)
        execute_process(COMMAND "${ASTER_BREW_EXECUTABLE}" --prefix libomp
            OUTPUT_VARIABLE _aster_libomp_prefix OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE _aster_brew_result ERROR_QUIET)
        if(_aster_brew_result EQUAL 0)
            find_path(OpenMP_CXX_INCLUDE_DIR omp.h
                HINTS "${_aster_libomp_prefix}/include")
            find_library(OpenMP_omp_LIBRARY NAMES omp
                HINTS "${_aster_libomp_prefix}/lib")
            if(OpenMP_CXX_INCLUDE_DIR AND OpenMP_omp_LIBRARY)
                set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
                set(OpenMP_CXX_LIB_NAMES omp)
            endif()
        endif()
    endif()
endif()
find_package(OpenMP REQUIRED COMPONENTS CXX)
message(STATUS "ASTER Boost headers: ${ASTER_BOOST_INCLUDE_DIR}")


# FindOpenMP's basic probe may miss compiler/runtime ABI mismatches.
# This tests the dynamic scheduling used by GL-ML, without running an experiment.
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state(RESET)
set(CMAKE_REQUIRED_LIBRARIES OpenMP::OpenMP_CXX)
check_cxx_source_compiles([=[
#include <omp.h>
int main() {
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < 100; ++i) { volatile int value = i; (void)value; }
    return 0;
}
]=] ASTER_OPENMP_DYNAMIC_LINKS)
cmake_pop_check_state()
if(NOT ASTER_OPENMP_DYNAMIC_LINKS)
    message(WARNING
        "OpenMP dynamic scheduling cannot link with this compiler/runtime. "
        "aster_iod_glml need a compatible toolchain. "
        "Other main targets can still be built explicitly. See doc/build.md.")
endif()
