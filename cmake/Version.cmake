# 完整版本由提交中的基础版本和实验后缀组成，不使用每台机器各自的缓存设置。
set(ASTER_FULL_VERSION "${PROJECT_VERSION}${ASTER_VERSION_SUFFIX}")
set(ASTER_SOURCE_COMMIT "unknown")
set(ASTER_SOURCE_STATE "archive")
find_package(Git QUIET)

if(Git_FOUND)
    # 仅识别当前源码自身的仓库，源码归档位于其他仓库内时仍标记为archive。
    execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}" rev-parse --show-toplevel
        OUTPUT_VARIABLE ASTER_GIT_ROOT OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    file(REAL_PATH "${PROJECT_SOURCE_DIR}" ASTER_SOURCE_ROOT)

    if(ASTER_GIT_ROOT STREQUAL ASTER_SOURCE_ROOT)
        execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}" rev-parse HEAD
            OUTPUT_VARIABLE ASTER_SOURCE_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
        execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}" status --porcelain --untracked-files=normal
            OUTPUT_VARIABLE ASTER_SOURCE_CHANGES OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
        set(ASTER_SOURCE_STATE "clean")

        if(ASTER_SOURCE_CHANGES)
            set(ASTER_SOURCE_STATE "modified")
        endif()

        # 提交、暂存区或当前分支引用变化时重新配置，避免旧构建显示过期提交。
        execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}" rev-parse --absolute-git-dir
            OUTPUT_VARIABLE ASTER_GIT_DIR OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
            "${ASTER_GIT_DIR}/HEAD" "${ASTER_GIT_DIR}/index")
        execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}" symbolic-ref -q HEAD
            OUTPUT_VARIABLE ASTER_GIT_REF OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)

        if(ASTER_GIT_REF)
            execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${PROJECT_SOURCE_DIR}" rev-parse --git-path "${ASTER_GIT_REF}"
                OUTPUT_VARIABLE ASTER_GIT_REF_PATH OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
            get_filename_component(ASTER_GIT_REF_PATH "${ASTER_GIT_REF_PATH}" ABSOLUTE BASE_DIR "${PROJECT_SOURCE_DIR}")

            if(EXISTS "${ASTER_GIT_REF_PATH}")
                set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${ASTER_GIT_REF_PATH}")
            endif()
        endif()
    endif()
endif()

configure_file("${PROJECT_SOURCE_DIR}/cmake/version_config.hpp.in" generated/version_config.hpp @ONLY)
message(STATUS "ASTER version=${ASTER_FULL_VERSION}; commit=${ASTER_SOURCE_COMMIT}; source=${ASTER_SOURCE_STATE}")
