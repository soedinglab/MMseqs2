cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(MMseqs CXX C)
message("-- Source Directory: ${CMAKE_CURRENT_SOURCE_DIR}")
message("-- Project Directory: ${PROJECT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(FRAMEWORK_ONLY 0 CACHE BOOL "Framework Mode")
set(HAVE_SANITIZER 0 CACHE BOOL "Have Sanitizers")
set(VERSION_OVERRIDE "" CACHE STRING "Override version string in help and usage messages")

#Sanitizers
if (${HAVE_SANITIZER})
    include(FindUBSan)
    include(FindASan)
    include(FindMSan)
    include(FindTSan)
endif ()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)

# find compiler
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # using Clang
    message("-- Compiler is clang(++)")
    set(CMAKE_COMPILER_IS_CLANG 1)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    message("-- Compiler is GNU ")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.0")
        message(FATAL_ERROR "Insufficient gcc version")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    message("-- Compiler is icc(++)")
    set(CMAKE_COMPILER_IS_ICC 1)
endif ()

# set flags
set(MMSEQS_CXX_FLAGS "-std=c++0x")
if (NOT ${HAVE_NEON})
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -m64")
endif ()

# Compiler-specific features
if (CMAKE_COMPILER_IS_CLANG)
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -stdlib=libc++")
endif ()

# Apple specific features
if (APPLE)
    # macOS SDK started using _Atomic (in ucred.h) which g++ does not support
    # __APPLE_API_STRICT_CONFORMANCE makes sysctl.h not include apis like ucred.h
    # See: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=89864
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -D__APPLE_API_STRICT_CONFORMANCE")
endif ()

# zstd
SET(ZSTD_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/zstd")
SET(CMAKE_INSTALL_LIBDIR bin)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/lib/zstd/build/cmake/CMakeModules")
OPTION(ZSTD_LEGACY_SUPPORT "LEGACY SUPPORT" OFF)
OPTION(ZSTD_BUILD_STATIC "BUILD STATIC LIBRARIES" ON)
OPTION(ZSTD_BUILD_SHARED "BUILD SHARED LIBRARIES" OFF)
OPTION(ZSTD_MULTITHREAD_SUPPORT "MULTITHREADING SUPPORT" OFF)
OPTION(ZSTD_BUILD_PROGRAMS "BUILD PROGRAMS" OFF)
OPTION(ZSTD_BUILD_CONTRIB "BUILD CONTRIB" OFF)
OPTION(ZSTD_BUILD_TESTS "BUILD TESTS" OFF)
add_subdirectory(lib/zstd/build/cmake/lib EXCLUDE_FROM_ALL)

# tinyexpr
OPTION(TE_NAT_LOG "Define the log function as natural logarithm." ON)
include_directories(lib/tinyexpr)
add_subdirectory(lib/tinyexpr EXCLUDE_FROM_ALL)

include_directories(lib)
include_directories(lib/kseq)
include_directories(lib/simd)
include_directories(lib/gzstream)
include_directories(lib/alp)
include_directories(lib/cacode)
include_directories(lib/ksw2)
add_subdirectory(lib/cacode)
add_subdirectory(lib/alp)
add_subdirectory(lib/ksw2)
add_subdirectory(data)
add_subdirectory(src)
if (NOT FRAMEWORK_ONLY)
    add_subdirectory(util)
endif ()
