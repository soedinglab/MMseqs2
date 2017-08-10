cmake_minimum_required(VERSION 2.8.9 FATAL_ERROR)
project(MMseqs CXX C )
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(HAVE_MPI 0 CACHE BOOL "Have MPI")
set(HAVE_SSE4_1 0 CACHE BOOL "Have SSE4.1")
set(HAVE_AVX2 0 CACHE BOOL "Have AVX2")
set(HAVE_TESTS 1 CACHE BOOL "Have Tests")
set(HAVE_SANITIZER 0 CACHE BOOL "Have Sanitizers")

set(MMSEQS_DIR ${PROJECT_SOURCE_DIR})
set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${PROJECT_SOURCE_DIR}/lib)
message("-- Sourcedir ${CMAKE_CURRENT_SOURCE_DIR}")
message("-- Projectdir ${PROJECT_SOURCE_DIR}")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif(NOT CMAKE_BUILD_TYPE)

# find compiler
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    # using Clang
    message("-- Compiler is clang(++)")
    set(CMAKE_COMPILER_IS_CLANG 1)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # using GCC
    message("-- Compiler is GNU ")
    set(CMAKE_COMPILER_IS_GNUCXX 1)
    set(CMAKE_COMPILER_IS_GNUCCC 1)
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.6.0")
        message(FATAL_ERROR "Insufficient gcc version")
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # using Intel C++
    message("-- Compiler is icc(++)")
    set(CMAKE_COMPILER_IS_ICC 1)
endif()

if(CYGWIN)
    add_definitions("-D_GNU_SOURCE")
endif()

#Sanitizers
if(${HAVE_SANITIZER})
    include(FindUBSan)
    include(FindASan)
    include(FindMSan)
    include(FindTSan)
endif()

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
add_definitions("-DGIT_SHA1=${GIT_SHA1}")

# openmp integration
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOPENMP=1 ${OpenMP_CXX_FLAGS}")
endif()

# needed for concat.h
include(CheckCXXSourceRuns)
check_cxx_source_runs("
        #include <stdlib.h>
        #include <fcntl.h>
        #include <stdio.h>

        int main()
        {
          FILE* tmpf = tmpfile();
          int input_desc = fileno(tmpf);
          int test =  posix_fadvise (input_desc, 0, 0, POSIX_FADV_SEQUENTIAL);
          fclose(tmpf);
        }"
        HAVE_POSIX_FADVISE)
if(HAVE_POSIX_FADVISE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_POSIX_FADVISE=1 ")
endif()


#SSE
if(${HAVE_AVX2})
    add_definitions("-DAVX2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -Wa,-q")
elseif(${HAVE_SSE4_1})
    add_definitions("-DSSE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
else()
    include(CheckSSEFeatures)
    add_definitions(${CXX_DFLAGS})
    if(${HAVE_AVX2_EXTENSIONS})
        add_definitions("-DAVX2")
        # debugging
        #   add_definitions("-DSSE")
    else()
        if(${HAVE_SSE4_1_EXTENSIONS})
            add_definitions("-DSSE")
        else()
            message( FATAL_ERROR "At least SSE4.2 is needed to compile!" )
        endif(${HAVE_SSE4_1_EXTENSIONS})
    endif(${HAVE_AVX2_EXTENSIONS})
endif()

# set flags
set(MMSEQS_WARNINGS "-pedantic -Wall -Wextra -Winline -Wdisabled-optimization -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -std=c++0x -m64 ${MMSEQS_WARNINGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL} -Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDBINFO} -O2 -g")

# Compiler-specific features
if(CMAKE_COMPILER_IS_CLANG)
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

# MPI
if(${HAVE_MPI})
    find_package(MPI QUIET)
    if(MPI_FOUND)
        set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
        set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
        include_directories(${MPI_INCLUDE_PATH})
        ADD_DEFINITIONS("-DHAVE_MPI")
        message("-- Found MPI")
    endif()
endif()

find_package(ZLIB QUIET)
if(ZLIB_FOUND)
    message("-- Found ZLIB")
    include_directories(${ZLIB_INCLUDE_DIRS})
    ADD_DEFINITIONS("-DHAVE_ZLIB")
else()
    message("-- Could not find ZLIB")
endif()

find_package(BZip2 QUIET)
if(BZIP2_FOUND)
    message("-- Found BZLIB")
    include_directories(${BZIP_INCLUDE_DIRS})
    ADD_DEFINITIONS("-DHAVE_BZLIB")
else()
    message("-- Could not find BZLIB")
endif()

include_directories(src/commons)
include_directories(src/prefiltering)
include_directories(src/alignment)
include_directories(src/clustering)
include_directories(src/util)
include_directories(lib)
include_directories(lib/getoptpp)
include_directories(lib/kseq)
include_directories(lib/simd)
include_directories(lib/gzstream)
include_directories(lib/alp)
include_directories(lib/cacode)
add_subdirectory(lib/cacode)
add_subdirectory(lib/alp)
add_subdirectory(data)
add_subdirectory(src)

install(PROGRAMS
	util/bash-completion.sh
	DESTINATION util
)
