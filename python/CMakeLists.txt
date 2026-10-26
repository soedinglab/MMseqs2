# Thin CMake entry point for the Python wheel (scikit-build-core).
#
# It builds the C library from the repository root (or from `native/` inside an
# sdist, where the C tree is copied by `sdist.force-include` in pyproject.toml)
# as one self-contained shared library, and installs *only* that library into
# the `carquet` Python package directory. Headers, the CMake package config,
# the pkg-config file and the CLI are deliberately left out of the wheel.
cmake_minimum_required(VERSION 3.16)
project(carquet_python LANGUAGES C)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/native/CMakeLists.txt")
    # Unpacked sdist: the C sources were copied under python/native/.
    set(CARQUET_C_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/native")
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../include/carquet/carquet.h")
    # Git checkout: python/ sits next to the C tree.
    get_filename_component(CARQUET_C_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
else()
    message(FATAL_ERROR "carquet C sources not found (expected ../include/carquet/carquet.h or native/)")
endif()
message(STATUS "carquet Python wheel: building C library from ${CARQUET_C_SOURCE_DIR}")

# Wheel build configuration. These are plain cache overrides of the root
# project's options; pyproject.toml passes the same values via cmake.define so
# a manual `cmake -S python` behaves identically.
set(CARQUET_BUILD_SHARED ON CACHE BOOL "wheel: single shared library" FORCE)
set(CARQUET_BUILD_CLI OFF CACHE BOOL "wheel: no CLI" FORCE)
set(CARQUET_BUILD_DEV OFF CACHE BOOL "wheel: no tests/examples/benchmarks" FORCE)
set(CARQUET_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(CARQUET_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(CARQUET_BUILD_BENCHMARKS OFF CACHE BOOL "" FORCE)
set(CARQUET_BUILD_INTEROP OFF CACHE BOOL "" FORCE)
set(CARQUET_BUNDLE_DEPS ON CACHE BOOL "wheel: static zstd/zlib/lz4/brotli inside libcarquet" FORCE)

# The root CMakeLists.txt only emits install() rules when it is the top-level
# project, so adding it as a subdirectory gives us the `carquet` target with no
# install rules of its own; the fetched codec subprojects' install rules are
# skipped too because only the component below is installed
# (tool.scikit-build.install.components in pyproject.toml).
add_subdirectory("${CARQUET_C_SOURCE_DIR}" carquet-c)

# Land the runtime library directly in the package directory
# (wheel.install-dir = "carquet" makes CMAKE_INSTALL_PREFIX that directory):
#   carquet/libcarquet.so | carquet/libcarquet.dylib | carquet/carquet.dll
install(TARGETS carquet
    COMPONENT carquet_python
    LIBRARY DESTINATION "."
    RUNTIME DESTINATION "."
    ARCHIVE DESTINATION "." EXCLUDE_FROM_ALL   # MSVC import library: not needed by ctypes
)
