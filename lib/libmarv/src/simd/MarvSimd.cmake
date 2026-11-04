# SIMD (CPU) backend: the shared host code compiled as C++ on the shims in src/shims, with striped SIMD gapless kernels
# (marv_simd_aligner.cpp) in place of the GPU kernels.
project(marv LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CheckCXXSourceCompiles)
check_cxx_source_compiles("int main() { _Float16 value = 1; return int(value) - 1; }" MARV_HAVE_FLOAT16)
if(NOT MARV_HAVE_FLOAT16)
    message(FATAL_ERROR "The SIMD host shims require _Float16 support. On x86 Linux use GCC 12 or newer, or a Clang compiler with _Float16 support.")
endif()
option(MARV_SIMD_NATIVE "Optimize SIMD for the build CPU (binary may require this CPU's instruction set)" OFF)

set(MARV_SIMD_DIR ${CMAKE_CURRENT_LIST_DIR})
get_filename_component(MARV_SRC_DIR ${MARV_SIMD_DIR} DIRECTORY)

set(MARV_SIMD_HOST_SOURCES
    ${MARV_SRC_DIR}/marv.cu
    ${MARV_SRC_DIR}/gpudatabase.cu
    ${MARV_SRC_DIR}/dbbatching.cu
    ${MARV_SRC_DIR}/dbdata.cpp
    ${MARV_SRC_DIR}/searches/score_only_search.cu
    ${MARV_SRC_DIR}/searches/score_endpos_search.cu)
set(MARV_SIMD_FORCE_INCLUDE "-include;${MARV_SRC_DIR}/shims/host_backend.h")
set_source_files_properties(${MARV_SIMD_HOST_SOURCES} PROPERTIES
    LANGUAGE CXX
    COMPILE_OPTIONS "-x;c++;${MARV_SIMD_FORCE_INCLUDE}")
set_source_files_properties(${MARV_SIMD_DIR}/marv_simd_aligner.cpp PROPERTIES
    COMPILE_OPTIONS "${MARV_SIMD_FORCE_INCLUDE}")

add_library(marv STATIC
    ${MARV_SIMD_HOST_SOURCES}
    ${MARV_SIMD_DIR}/marv_simd_runtime.cpp
    ${MARV_SIMD_DIR}/marv_simd_aligner.cpp
    ${MARV_SRC_DIR}/marv.h)
target_compile_definitions(marv PRIVATE __SIMD_BACKEND__=1 NO_NVTOOLSEXT AMD_INLINE_ARCH_NAMESPACE_NAME=SIMD_)
target_include_directories(marv PUBLIC ${MARV_SRC_DIR})
target_include_directories(marv PRIVATE ${MARV_SRC_DIR}/shims ${MARV_SRC_DIR}/hpc_helpers)
find_package(OpenMP REQUIRED)
target_link_libraries(marv PUBLIC OpenMP::OpenMP_CXX)
if(MARV_SIMD_NATIVE)
    set_property(SOURCE ${MARV_SIMD_DIR}/marv_simd_aligner.cpp APPEND PROPERTY COMPILE_OPTIONS -march=native)
endif()

set(LIBRARY_ONLY 0 CACHE BOOL "Do not build executables")
if(NOT LIBRARY_ONLY)
    add_executable(test_gapless ${MARV_SRC_DIR}/tests/test_gapless.cpp)
    target_link_libraries(test_gapless PRIVATE marv)
    enable_testing()
    add_test(NAME gapless_exactness_simd128 COMMAND test_gapless)
    set_tests_properties(gapless_exactness_simd128 PROPERTIES ENVIRONMENT MARV_SIMD_WIDTH=16)
    add_test(NAME gapless_exactness_simd256 COMMAND test_gapless)
    set_tests_properties(gapless_exactness_simd256 PROPERTIES ENVIRONMENT MARV_SIMD_WIDTH=32)
    # Match the aligner translation unit's ISA flags when selecting width tests.
    set(_marv_saved_required_flags "${CMAKE_REQUIRED_FLAGS}")
    if(MARV_SIMD_NATIVE)
        string(APPEND CMAKE_REQUIRED_FLAGS " -march=native")
    endif()
    unset(MARV_SIMD_HAS_AVX512BW CACHE)
    check_cxx_source_compiles("#ifndef __AVX512BW__
#error AVX512BW is not enabled
#endif
int main() { return 0; }" MARV_SIMD_HAS_AVX512BW)
    set(CMAKE_REQUIRED_FLAGS "${_marv_saved_required_flags}")
    if(MARV_SIMD_HAS_AVX512BW)
        add_test(NAME gapless_exactness_simd512 COMMAND test_gapless)
        set_tests_properties(gapless_exactness_simd512 PROPERTIES ENVIRONMENT MARV_SIMD_WIDTH=64)
    endif()
endif()
