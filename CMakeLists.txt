cmake_minimum_required(VERSION 3.26)
project(ASTER VERSION 1.4.1 LANGUAGES C CXX)
# 实验后缀随源码保存，正式发布使用空字符串。
set(ASTER_VERSION_SUFFIX "")
include(cmake/Version.cmake)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(cmake/Dependencies.cmake)
include(cmake/MathBackend.cmake)
target_link_libraries(aster_sofa PUBLIC ASTER::math_backend)
set_target_properties(aster_sofa_validation PROPERTIES LINKER_LANGUAGE CXX)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
add_library(aster_core STATIC
    src/accel_complete.cpp
    src/accel_drag.cpp
    src/accel_gravity.cpp
    src/accel_solar_radiation.cpp
    src/accel_thirdbody.cpp
    src/analyt_averaging.cpp
    src/analyt_gravity.cpp
    src/analyt_lunisolar.cpp
    src/analyt_solar_radiation.cpp
    src/core_bodn.cpp
    src/core_iod.cpp
    src/core_iod_solver.cpp
    src/core_propagation.cpp
    src/data_gravityfield.cpp
    src/earth_eop.cpp
    src/frame_definition.cpp
    src/frame_earth.cpp
    src/integrator_adams.cpp
    src/integrator_rkf78.cpp
    src/io_common.cpp
    src/io_observation_arc.cpp
    src/io_settings.cpp
    src/math_basic.cpp
    src/math_statistics.cpp
    src/model_atmosphere.cpp
    src/model_ephemeris_analytic.cpp
    src/model_ephemeris_jpl.cpp
    src/model_cislunar_observation.cpp
    src/model_observation.cpp
    src/model_propagation.cpp
    src/model_satellite.cpp
    src/orbit_basic.cpp
    src/orbit_equation.cpp
    src/orbit_expansion.cpp
    src/orbit_geometry.cpp
    src/orbit_statedef.cpp
    src/orbit_stm.cpp
    src/time_transform.cpp
    src/utils_constants.cpp
    src/utils_overwatch.cpp
)
add_library(ASTER::core ALIAS aster_core)
target_include_directories(aster_core PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_include_directories(aster_core SYSTEM PUBLIC "${ASTER_BOOST_INCLUDE_DIR}")
target_link_libraries(aster_core PUBLIC Eigen3::Eigen aster_sofa)
function(aster_add_app name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE ASTER::core)
    target_compile_definitions(${name} PRIVATE ASTER_PROGRAM_VERSION="${ASTER_FULL_VERSION}")
    set_target_properties(${name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
endfunction()
add_subdirectory(apps)
include(CTest)
if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
