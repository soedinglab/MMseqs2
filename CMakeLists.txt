cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(MMseqs CXX C)
message("-- Source Directory: ${CMAKE_CURRENT_SOURCE_DIR}")
message("-- Project Directory: ${PROJECT_SOURCE_DIR}")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(FRAMEWORK_ONLY 0 CACHE BOOL "Framework Mode")
set(HAVE_SANITIZER 0 CACHE BOOL "Have Sanitizers")
set(INSTALL_UTIL 1 CACHE BOOL "Install util scripts")
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
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
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

# see https://wiki.debian.org/ArchitectureSpecificsMemo for char signedness
set(MMSEQS_CXX_FLAGS "-std=c++0x")
if (CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm.*|ARM.*)")
    set(ARM 1)
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -DNEON=1 -DSSE=1 -fsigned-char -mfpu=neon")
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*)")
    set(ARM 1)
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -DNEON=1 -DSSE=1 -fsigned-char")
elseif (CMAKE_SYSTEM_PROCESSOR MATCHES "PPC64*|ppc64*|powerpc64*")
    set(PPC64 1)
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -DSSE=1 -fsigned-char -std=gnu++0x -mcpu=power8 -mvsx")
elseif (EMSCRIPTEN)
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -DSSE=1 -DWASM=1 -msimd128 -s WASM=1 -s ASSERTIONS=1")
else ()
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -m64")
endif ()

# Compiler-specific features
if (CMAKE_COMPILER_IS_CLANG AND (NOT EMSCRIPTEN))
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

if(CMAKE_SYSTEM_NAME STREQUAL "FreeBSD")
    set(MMSEQS_CXX_FLAGS "${MMSEQS_CXX_FLAGS} -D_WITH_GETLINE")
endif()

# zstd
# We use ZSTD_findDecompressedSize which is only available with ZSTD_STATIC_LINKING_ONLY
# Thus we cannot use a system provided libzstd
set(ZSTD_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib/zstd")
set(CMAKE_INSTALL_LIBDIR bin)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/lib/zstd/build/cmake/CMakeModules")
option(ZSTD_LEGACY_SUPPORT "LEGACY SUPPORT" OFF)
option(ZSTD_BUILD_STATIC "BUILD STATIC LIBRARIES" ON)
option(ZSTD_BUILD_SHARED "BUILD SHARED LIBRARIES" OFF)
option(ZSTD_MULTITHREAD_SUPPORT "MULTITHREADING SUPPORT" OFF)
option(ZSTD_BUILD_PROGRAMS "BUILD PROGRAMS" OFF)
option(ZSTD_BUILD_CONTRIB "BUILD CONTRIB" OFF)
option(ZSTD_BUILD_TESTS "BUILD TESTS" OFF)
include_directories(lib/zstd/lib)
add_subdirectory(lib/zstd/build/cmake/lib EXCLUDE_FROM_ALL)

# tinyexpr
include_directories(lib/tinyexpr)
add_subdirectory(lib/tinyexpr EXCLUDE_FROM_ALL)

# microtar
include_directories(lib/microtar)
add_subdirectory(lib/microtar)

# default cygwin allocator (dlmalloc) locks on every allocation and destroys MT performance
if (CYGWIN)
    add_subdirectory(lib/nedmalloc)
endif()

include_directories(lib)
include_directories(lib/simd)
include_directories(lib/gzstream)
include_directories(lib/alp)
include_directories(lib/cacode)
include_directories(lib/ksw2)
include_directories(lib/xxhash)

add_subdirectory(lib/cacode)
add_subdirectory(lib/alp)
add_subdirectory(lib/ksw2)
add_subdirectory(data)
add_subdirectory(src)
if (NOT FRAMEWORK_ONLY AND INSTALL_UTIL)
    add_subdirectory(util)
endif ()
