cmake_minimum_required(VERSION 2.8.9 FATAL_ERROR)
project(MMseqs CXX)
message("-- Source Directory: ${CMAKE_CURRENT_SOURCE_DIR}")
message("-- Project Directory: ${PROJECT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(FRAMEWORK_ONLY 0 CACHE BOOL "Framework Mode")
set(HAVE_SANITIZER 0 CACHE BOOL "Have Sanitizers")

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
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.6.0")
        message(FATAL_ERROR "Insufficient gcc version")
    endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    message("-- Compiler is icc(++)")
    set(CMAKE_COMPILER_IS_ICC 1)
endif ()

# set flags
set(MMSEQS_CXX_FLAGS "-std=c++0x -m64")

# Compiler-specific features
if (CMAKE_COMPILER_IS_CLANG)
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -stdlib=libc++")
endif ()

include_directories(lib)
include_directories(lib/kseq)
include_directories(lib/simd)
include_directories(lib/gzstream)
include_directories(lib/alp)
include_directories(lib/cacode)
add_subdirectory(lib/cacode)
add_subdirectory(lib/alp)
add_subdirectory(data)
add_subdirectory(src)
if (NOT FRAMEWORK_ONLY)
    add_subdirectory(util)
endif ()
