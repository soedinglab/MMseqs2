# Metal backend: shared host code compiled as C++ on the shims in src/shims, kernels in gapless_metal.metal.
# Needs the Metal compiler from Xcode (DEVELOPER_DIR=/Applications/Xcode.app/Contents/Developer if not found).
if(NOT APPLE)
    message(FATAL_ERROR "ENABLE_METAL requires macOS")
endif()
project(marv LANGUAGES CXX OBJCXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_OBJCXX_STANDARD 20)

set(MARV_METAL_DIR ${CMAKE_CURRENT_LIST_DIR})
get_filename_component(MARV_SRC_DIR ${MARV_METAL_DIR} DIRECTORY)

find_program(XCRUN_EXECUTABLE xcrun REQUIRED)
execute_process(COMMAND ${XCRUN_EXECUTABLE} -sdk macosx --find metal
                RESULT_VARIABLE MARV_METAL_FIND_RESULT OUTPUT_QUIET ERROR_QUIET)
if(NOT MARV_METAL_FIND_RESULT EQUAL 0)
    message(FATAL_ERROR "The Metal shader compiler was not found (xcrun -sdk macosx --find metal failed). "
                        "Install Xcode and, if needed, set DEVELOPER_DIR=/Applications/Xcode.app/Contents/Developer.")
endif()
find_program(XXD_EXECUTABLE xxd REQUIRED)

set(MARV_METAL_AIR ${CMAKE_CURRENT_BINARY_DIR}/gapless_metal.air)
# no extension: xxd derives the symbol names marv_metallib / marv_metallib_len
set(MARV_METALLIB  ${CMAKE_CURRENT_BINARY_DIR}/marv_metallib)
add_custom_command(
    OUTPUT ${MARV_METALLIB}
    COMMAND ${XCRUN_EXECUTABLE} -sdk macosx metal -O2 -I ${MARV_METAL_DIR}
            -c ${MARV_METAL_DIR}/gapless_metal.metal -o ${MARV_METAL_AIR}
    COMMAND ${XCRUN_EXECUTABLE} -sdk macosx metallib ${MARV_METAL_AIR} -o ${MARV_METALLIB}
    DEPENDS ${MARV_METAL_DIR}/gapless_metal.metal ${MARV_METAL_DIR}/marv_metal_shared.h
    COMMENT "Compiling Metal kernels")

set(MARV_METALLIB_EMBED ${CMAKE_CURRENT_BINARY_DIR}/marv_metallib_embed.cpp)
add_custom_command(
    OUTPUT ${MARV_METALLIB_EMBED}
    COMMAND ${XXD_EXECUTABLE} -i marv_metallib > ${MARV_METALLIB_EMBED}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${MARV_METALLIB}
    COMMENT "Embedding the Metal library")
set_source_files_properties(${MARV_METALLIB_EMBED} PROPERTIES GENERATED TRUE)

set(MARV_METAL_HOST_SOURCES
    ${MARV_SRC_DIR}/marv.cu
    ${MARV_SRC_DIR}/gpudatabase.cu
    ${MARV_SRC_DIR}/dbbatching.cu
    ${MARV_SRC_DIR}/dbdata.cpp
    ${MARV_SRC_DIR}/searches/score_only_search.cu
    ${MARV_SRC_DIR}/searches/score_endpos_search.cu)
set(MARV_METAL_FORCE_INCLUDE "-include;${MARV_SRC_DIR}/shims/host_backend.h")
set_source_files_properties(${MARV_METAL_HOST_SOURCES} PROPERTIES
    LANGUAGE CXX
    COMPILE_OPTIONS "-x;c++;${MARV_METAL_FORCE_INCLUDE}")
set_source_files_properties(${MARV_METAL_DIR}/marv_metal_aligner.mm PROPERTIES
    COMPILE_OPTIONS "-fobjc-arc;${MARV_METAL_FORCE_INCLUDE}")
set_source_files_properties(${MARV_METAL_DIR}/cuda_metal_runtime.mm PROPERTIES
    COMPILE_OPTIONS "-fobjc-arc")

add_library(marv STATIC
    ${MARV_METAL_HOST_SOURCES}
    ${MARV_METAL_DIR}/cuda_metal_runtime.mm
    ${MARV_METAL_DIR}/marv_metal_aligner.mm
    ${MARV_METALLIB_EMBED}
    ${MARV_SRC_DIR}/marv.h)
target_compile_definitions(marv PRIVATE __METAL_BACKEND__=1 NO_NVTOOLSEXT AMD_INLINE_ARCH_NAMESPACE_NAME=METAL_)
target_include_directories(marv PUBLIC ${MARV_SRC_DIR})
target_include_directories(marv PRIVATE ${MARV_METAL_DIR} ${MARV_SRC_DIR}/shims ${MARV_SRC_DIR}/hpc_helpers)
target_link_libraries(marv PUBLIC "-framework Metal" "-framework Foundation")

find_package(OpenMP REQUIRED)
target_link_libraries(marv PUBLIC OpenMP::OpenMP_CXX)

set(LIBRARY_ONLY 0 CACHE BOOL "Do not build executables")
if(NOT LIBRARY_ONLY)
    add_executable(test_gapless ${MARV_SRC_DIR}/tests/test_gapless.cpp)
    target_link_libraries(test_gapless PRIVATE marv)
    enable_testing()
    add_test(NAME gapless_exactness COMMAND test_gapless)
endif()
