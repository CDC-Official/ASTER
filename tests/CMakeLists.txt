find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(NAME iod_gl COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_case.py" --case iod_gl --binary $<TARGET_FILE:aster_iod_gl> --output "${PROJECT_BINARY_DIR}/test-output/iod_gl")
set_tests_properties(iod_gl PROPERTIES TIMEOUT 120)
add_test(NAME iod_glml COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_case.py" --case iod_glml --binary $<TARGET_FILE:aster_iod_glml> --output "${PROJECT_BINARY_DIR}/test-output/iod_glml")
set_tests_properties(iod_glml PROPERTIES TIMEOUT 120)
add_test(NAME odn_batch COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_case.py" --case odn_batch --binary $<TARGET_FILE:aster_odn_batch> --output "${PROJECT_BINARY_DIR}/test-output/odn_batch")
set_tests_properties(odn_batch PROPERTIES TIMEOUT 120)
add_test(NAME opn_cartesian COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_case.py" --case opn_cartesian --binary $<TARGET_FILE:aster_opn_cartesian> --output "${PROJECT_BINARY_DIR}/test-output/opn_cartesian")
set_tests_properties(opn_cartesian PROPERTIES TIMEOUT 120)
add_test(NAME opn_nonsingular2 COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_case.py" --case opn_nonsingular2 --binary $<TARGET_FILE:aster_opn_nonsingular2> --output "${PROJECT_BINARY_DIR}/test-output/opn_nonsingular2")
set_tests_properties(opn_nonsingular2 PROPERTIES TIMEOUT 120)
add_test(NAME opn_cislunar COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_case.py" --case opn_cislunar --binary $<TARGET_FILE:aster_opn_cislunar> --output "${PROJECT_BINARY_DIR}/test-output/opn_cislunar")
set_tests_properties(opn_cislunar PROPERTIES TIMEOUT 120)
add_test(NAME convert_osculating_to_mean COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_case.py" --case convert_osculating_to_mean --binary $<TARGET_FILE:aster_convert_osculating_to_mean> --output "${PROJECT_BINARY_DIR}/test-output/convert_osculating_to_mean")
set_tests_properties(convert_osculating_to_mean PROPERTIES TIMEOUT 120)
add_test(NAME convert_keplerian_to_cartesian COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_case.py" --case convert_keplerian_to_cartesian --binary $<TARGET_FILE:aster_convert_keplerian_to_cartesian> --output "${PROJECT_BINARY_DIR}/test-output/convert_keplerian_to_cartesian")
set_tests_properties(convert_keplerian_to_cartesian PROPERTIES TIMEOUT 120)
add_test(NAME convert_cartesian_to_keplerian COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/run_case.py" --case convert_cartesian_to_keplerian --binary $<TARGET_FILE:aster_convert_cartesian_to_keplerian> --output "${PROJECT_BINARY_DIR}/test-output/convert_cartesian_to_keplerian")
set_tests_properties(convert_cartesian_to_keplerian PROPERTIES TIMEOUT 120)
add_executable(aster_test_frame_coordinates core/frame_coordinates.cpp)
target_link_libraries(aster_test_frame_coordinates PRIVATE ASTER::core)
add_test(NAME core_frame_coordinates COMMAND aster_test_frame_coordinates)
set_tests_properties(core_frame_coordinates PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" TIMEOUT 180)
add_executable(aster_test_orbit_averaging core/orbit_averaging.cpp)
target_link_libraries(aster_test_orbit_averaging PRIVATE ASTER::core)
add_test(NAME core_orbit_averaging COMMAND aster_test_orbit_averaging)
set_tests_properties(core_orbit_averaging PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" TIMEOUT 180)
add_executable(aster_test_batch_determination core/batch_determination.cpp)
target_link_libraries(aster_test_batch_determination PRIVATE ASTER::core)
add_test(NAME core_batch_determination COMMAND aster_test_batch_determination)
set_tests_properties(core_batch_determination PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" TIMEOUT 180)
add_executable(aster_test_batch_parameters core/batch_parameters.cpp)
target_link_libraries(aster_test_batch_parameters PRIVATE ASTER::core)
add_test(NAME core_batch_parameters COMMAND aster_test_batch_parameters)
set_tests_properties(core_batch_parameters PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" TIMEOUT 180)
add_executable(aster_test_rkf78_safety core/rkf78_safety.cpp)
target_link_libraries(aster_test_rkf78_safety PRIVATE ASTER::core)
add_test(NAME core_rkf78_safety COMMAND aster_test_rkf78_safety)
set_tests_properties(core_rkf78_safety PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" TIMEOUT 180)
add_executable(aster_test_adams_safety core/adams_safety.cpp)
target_link_libraries(aster_test_adams_safety PRIVATE ASTER::core)
add_test(NAME core_adams_safety COMMAND aster_test_adams_safety)
set_tests_properties(core_adams_safety PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" TIMEOUT 180)
add_executable(aster_test_iod_eccentricity core/iod_eccentricity.cpp)
target_link_libraries(aster_test_iod_eccentricity PRIVATE ASTER::core)
add_test(NAME core_iod_eccentricity COMMAND aster_test_iod_eccentricity "${CMAKE_CURRENT_SOURCE_DIR}/core/iod_eccentricity_reference.tsv" "${PROJECT_BINARY_DIR}/test-output/eccentricity")
set_tests_properties(core_iod_eccentricity PROPERTIES WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}" TIMEOUT 180)
add_test(NAME invalid_inputs COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_failures.py" --bin "${PROJECT_BINARY_DIR}/bin" --output "${PROJECT_BINARY_DIR}/test-output/failures")
set_tests_properties(invalid_inputs PROPERTIES TIMEOUT 120)
add_test(NAME data_integrity COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_data.py")
add_test(NAME iod_batch_handoff COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/check_handoff.py" --binary $<TARGET_FILE:aster_odn_batch> --output "${PROJECT_BINARY_DIR}/test-output/handoff")
set_tests_properties(iod_batch_handoff PROPERTIES TIMEOUT 120)

if(ASTER_MATH_OPENLIBM)
add_test(NAME math_backend_iod COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/math_backend_cases/run_checks.py"
    --bin-dir "${PROJECT_BINARY_DIR}/bin" --output "${PROJECT_BINARY_DIR}/test-output/math_backend")
set_tests_properties(math_backend_iod PROPERTIES TIMEOUT 1200)
endif()
add_test(NAME math_backend_sofa COMMAND aster_sofa_validation)

# 复用原观测生成验收输入与几何、加噪检查。
add_executable(aster_test_observation_geometry_probe core/observation_geometry_probe.cpp)
target_link_libraries(aster_test_observation_geometry_probe PRIVATE ASTER::core)
add_test(NAME observation_runtime COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/observation_runtime.py"
    --binary $<TARGET_FILE:aster_generate_observations> --project "${PROJECT_SOURCE_DIR}"
    --jpl "${PROJECT_SOURCE_DIR}/data/ephemeris/JPLEPH_C.430")
add_test(NAME observation_geometry COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/observation_geometry.py"
    --binary $<TARGET_FILE:aster_generate_observations>
    --probe $<TARGET_FILE:aster_test_observation_geometry_probe> --project "${PROJECT_SOURCE_DIR}"
    --jpl "${PROJECT_SOURCE_DIR}/data/ephemeris/JPLEPH_C.430")
add_test(NAME observation_angular_noise COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/angular_noise_runtime.py"
    --script "${PROJECT_SOURCE_DIR}/tools/add_angular_noise.py")
set_tests_properties(observation_runtime observation_geometry PROPERTIES TIMEOUT 120)
set_tests_properties(observation_angular_noise PROPERTIES TIMEOUT 60)

add_executable(aster_test_cislunar_interpolation core/cislunar_interpolation.cpp)
target_link_libraries(aster_test_cislunar_interpolation PRIVATE ASTER::core)
add_test(NAME cislunar_interpolation COMMAND aster_test_cislunar_interpolation
    "${PROJECT_SOURCE_DIR}/data/ephemeris/JPLEPH_C.430")
set_tests_properties(cislunar_interpolation PROPERTIES TIMEOUT 30)

add_executable(aster_test_cislunar_observation core/cislunar_observation.cpp)
target_link_libraries(aster_test_cislunar_observation PRIVATE ASTER::core)
add_test(NAME core_cislunar_observation COMMAND aster_test_cislunar_observation
    "${PROJECT_SOURCE_DIR}/tests/generate_cislunar_observations/input/JPLEPH_C.440")
add_test(NAME cislunar_observation_runtime COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/check_cislunar_observations.py"
    --binary $<TARGET_FILE:aster_generate_cislunar_observations>
    --case "${CMAKE_CURRENT_SOURCE_DIR}/generate_cislunar_observations"
    --output "${PROJECT_BINARY_DIR}/test-output/cislunar_observation")
set_tests_properties(core_cislunar_observation cislunar_observation_runtime PROPERTIES TIMEOUT 120)

# 有量纲公式验证中心引力和 GGM03C 参数的时标约定。
add_executable(aster_test_gravity_conventions core/gravity_conventions.cpp)
target_link_libraries(aster_test_gravity_conventions PRIVATE ASTER::core)
add_test(NAME gravity_conventions COMMAND aster_test_gravity_conventions
    "${PROJECT_SOURCE_DIR}/data/eop/EOP-All_250516.txt"
    "${PROJECT_SOURCE_DIR}/data/ephemeris/JPLEPH_C.430")
set_tests_properties(gravity_conventions PROPERTIES TIMEOUT 30)

# 固定点 NAO、独立势梯度及两套 GCRS 动力学的固体潮附加修正。
add_executable(aster_test_earth_solid_tide core/earth_solid_tide.cpp)
target_link_libraries(aster_test_earth_solid_tide PRIVATE ASTER::core)
add_test(NAME core_earth_solid_tide COMMAND aster_test_earth_solid_tide
    "${PROJECT_SOURCE_DIR}/data/eop/EOP-All_250516.txt"
    "${PROJECT_SOURCE_DIR}/data/ephemeris/JPLEPH_C.430"
    "${PROJECT_SOURCE_DIR}/tests/earth_solid_tide/nao_acceleration.tsv")
set_tests_properties(core_earth_solid_tide PROPERTIES TIMEOUT 30)

# 六个方向的测站位置转换，与 SOFA 及现有初轨站址计算交叉核对。
add_executable(aster_test_station_coordinates core/station_coordinates.cpp)
target_link_libraries(aster_test_station_coordinates PRIVATE ASTER::core)
add_test(NAME core_station_coordinates COMMAND aster_test_station_coordinates
    "${PROJECT_SOURCE_DIR}/data/eop/EOP-All_250516.txt")
set_tests_properties(core_station_coordinates PROPERTIES TIMEOUT 30)

add_test(NAME station_coordinates_runtime COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/check_station_coordinates.py"
    --binary $<TARGET_FILE:aster_convert_coordinates>
    --probe $<TARGET_FILE:aster_test_station_coordinates>
    --output "${PROJECT_BINARY_DIR}/test-output/station_coordinates")
set_tests_properties(station_coordinates_runtime PROPERTIES TIMEOUT 60)

# 坐标和状态转换入口仅用两套单位的小样本检查往返。
add_test(NAME conversion_roundtrips COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/check_conversion_roundtrips.py"
    --bin-dir "${PROJECT_BINARY_DIR}/bin"
    --output "${PROJECT_BINARY_DIR}/test-output/conversion_roundtrips")
set_tests_properties(conversion_roundtrips PROPERTIES TIMEOUT 60)
