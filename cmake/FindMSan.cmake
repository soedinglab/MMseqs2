#
# The MIT License (MIT)
#
# Copyright (c) 2013 Matthew Arsenault
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

# This module tests if memory sanitizer is supported by the compiler,
# and creates a MSan build type (i.e. set CMAKE_BUILD_TYPE=MSan to use
# it). This sets the following variables:
#
# CMAKE_C_FLAGS_MSAN - Flags to use for C with msan
# CMAKE_CXX_FLAGS_MSAN  - Flags to use for C++ with msan
# HAVE_MEMORY_SANITIZER - True or false if the MSan build type is available

include(CheckCXXCompilerFlag)

# Set -Werror to catch "argument unused during compilation" warnings
set(CMAKE_REQUIRED_FLAGS "-Werror -fsanitize=memory") # Also needs to be a link flag for test to pass
check_cxx_compiler_flag("-fsanitize=memory" HAVE_FLAG_SANITIZE_MEMORY)

unset(CMAKE_REQUIRED_FLAGS)

if(HAVE_FLAG_SANITIZE_MEMORY)
  # Clang 3.2+ use this version
  set(MEMORY_SANITIZER_FLAG "-fsanitize=memory")
endif()

if(NOT MEMORY_SANITIZER_FLAG)
  return()
else(NOT MEMORY_SANITIZER_FLAG)
  set(HAVE_MEMORY_SANITIZER TRUE)
endif()

set(HAVE_MEMORY_SANITIZER TRUE)

set(CMAKE_C_FLAGS_MSAN "-O1 -g ${MEMORY_SANITIZER_FLAG} -fno-omit-frame-pointer -fno-optimize-sibling-calls"
    CACHE STRING "Flags used by the C compiler during MSan builds."
    FORCE)
set(CMAKE_CXX_FLAGS_MSAN "-O1 -g ${MEMORY_SANITIZER_FLAG} -fno-omit-frame-pointer -fno-optimize-sibling-calls"
    CACHE STRING "Flags used by the C++ compiler during MSan builds."
    FORCE)
set(CMAKE_EXE_LINKER_FLAGS_MSAN "${MEMORY_SANITIZER_FLAG}"
    CACHE STRING "Flags used for linking binaries during MSan builds."
    FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_MSAN "${MEMORY_SANITIZER_FLAG}"
    CACHE STRING "Flags used by the shared libraries linker during MSan builds."
    FORCE)
mark_as_advanced(CMAKE_C_FLAGS_MSAN
                 CMAKE_CXX_FLAGS_MSAN
                 CMAKE_EXE_LINKER_FLAGS_MSAN
                 CMAKE_SHARED_LINKER_FLAGS_MSAN)
