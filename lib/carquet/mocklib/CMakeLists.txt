cmake_minimum_required(VERSION 3.16)
project(metricstore VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ---------------------------------------------------------------------------
# Locate the carquet library.
#
# Three resolution strategies, in order:
#   1. A `carquet` target already defined by an enclosing project.
#   2. An installed package (find_package(carquet CONFIG)).
#   3. The in-tree source tree, added as a subdirectory (default: parent dir).
#
# Override the source location with -DCARQUET_SOURCE_DIR=/path/to/carquet.
# ---------------------------------------------------------------------------
set(CARQUET_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." CACHE PATH
    "Path to the carquet source tree used when no installed package is found")

if(TARGET carquet)
    message(STATUS "metricstore: using pre-existing carquet target")
else()
    find_package(carquet CONFIG QUIET)
    if(carquet_FOUND)
        message(STATUS "metricstore: using installed carquet package")
    elseif(EXISTS "${CARQUET_SOURCE_DIR}/CMakeLists.txt")
        message(STATUS "metricstore: building carquet from ${CARQUET_SOURCE_DIR}")
        set(CARQUET_BUILD_CLI OFF CACHE BOOL "" FORCE)
        set(CARQUET_BUILD_DEV OFF CACHE BOOL "" FORCE)
        add_subdirectory("${CARQUET_SOURCE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/carquet")
    else()
        message(FATAL_ERROR
            "Could not find carquet. Install it, or set -DCARQUET_SOURCE_DIR "
            "to its source tree.")
    endif()
endif()

# ---------------------------------------------------------------------------
# The mockup library.
# ---------------------------------------------------------------------------
add_library(metricstore STATIC
    src/metricstore.c
    src/ms_columnar.c
    src/ms_nested.c
    src/ms_arrow.c
    src/ms_diag.c)
target_include_directories(metricstore
    PUBLIC  "${CMAKE_CURRENT_SOURCE_DIR}/include"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(metricstore PUBLIC carquet)

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(metricstore PRIVATE -Wall -Wextra)
endif()

# ---------------------------------------------------------------------------
# Demo application.
# ---------------------------------------------------------------------------
add_executable(telemetry_demo app/telemetry_demo.c)
target_link_libraries(telemetry_demo PRIVATE metricstore)

add_executable(write_sample app/write_sample.c)
target_link_libraries(write_sample PRIVATE metricstore)

# ---------------------------------------------------------------------------
# Round-trip test.
# ---------------------------------------------------------------------------
enable_testing()
add_executable(test_roundtrip tests/test_roundtrip.c)
target_link_libraries(test_roundtrip PRIVATE metricstore)
add_test(NAME roundtrip COMMAND test_roundtrip)
