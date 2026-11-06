cmake_minimum_required(VERSION 3.26)
project(ASTERJ2PredictionValidation LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 同一预报入口可以链接原版或修正版源码，旧工作树无需增加文件。
get_filename_component(ASTER_DEFAULT_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(ASTER_SOURCE_DIR "${ASTER_DEFAULT_SOURCE}" CACHE PATH "ASTER source tree to validate")
set(BUILD_TESTING OFF CACHE BOOL "Build only the prediction probe" FORCE)
add_subdirectory("${ASTER_SOURCE_DIR}" aster)
add_executable(aster_test_iod_lagrange_prediction main_prediction.cpp)
target_link_libraries(aster_test_iod_lagrange_prediction PRIVATE ASTER::core)

# 与所链接源码的私有OpenLibm放在同一目录，保持两端运行方式一致。
set_target_properties(aster_test_iod_lagrange_prediction PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/aster/bin")
