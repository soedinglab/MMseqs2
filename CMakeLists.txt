cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(MMseqs CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(HAVE_MPI 0 CACHE BOOL "Have MPI")
set(HAVE_SSE4_1 0 CACHE BOOL "Have SSE4.1")
set(HAVE_AVX2 0 CACHE BOOL "Have AVX2")
set(HAVE_TESTS 1 CACHE BOOL "Have Tests")
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


#include(UseLATEX)
#set(PDFLATEX_COMPILER xelatex)
#add_latex_document(userguide.tex)

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
add_definitions("-DGIT_SHA1=${GIT_SHA1}")

include_directories(src/commons)
include_directories(src/prefiltering)
include_directories(src/alignment)
include_directories(src/clustering)
include_directories(src/util)
include_directories(lib)
include_directories(lib/getoptpp)
include_directories(lib/kseq)
include_directories(lib/simd)
include_directories(lib/seqan)
include_directories(lib/gzstream)
add_subdirectory(src)





install(PROGRAMS
	util/bash-completion.sh
	DESTINATION util
)
