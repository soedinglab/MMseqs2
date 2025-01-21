cmake_minimum_required(VERSION 3.8)
project(CUDA_INCLUDE_DIRS LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

set(LIBRARY_ONLY 0 CACHE BOOL "Do not build executables")

set(NVCC_FLAGS)
list(APPEND NVCC_FLAGS --extended-lambda)
list(APPEND NVCC_FLAGS --expt-relaxed-constexpr)
list(APPEND NVCC_FLAGS -rdc=true)
string(TOLOWER "${CMAKE_BUILD_TYPE}" BUILD_TYPE_LOWER)
if(BUILD_TYPE_LOWER STREQUAL "debug" OR BUILD_TYPE_LOWER STREQUAL "relwithdebinfo")
    list(APPEND NVCC_FLAGS -lineinfo)
endif()

set(NVCC_FLAGS_OMP ${NVCC_FLAGS})
list(APPEND NVCC_FLAGS_OMP -Xcompiler -fopenmp)


set(ALIGN_SOURCES
    kernels.cuh
    blosum.cu
    pssm.cuh
    pssmkernels_gapless.cuh
    pssmkernels_gapless_instantiation_dpx.cu
    pssmkernels_gapless_instantiation_half2.cu
    pssmkernels_gapless_instantiation_dpx_kernelparamzero.cu
    pssmkernels_gapless_instantiation_half2_kernelparamzero.cu
    pssmkernels_smithwaterman.cuh
    pssmkernels_smithwaterman_instantiation_float.cu
    pssmkernels_smithwaterman_instantiation_dpx.cu
    dbdata.cpp
)

add_library(marv ${ALIGN_SOURCES} marv.cu marv.h)
target_compile_options(marv PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${NVCC_FLAGS}>)
set_target_properties(marv
    PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
)
target_compile_definitions(marv PRIVATE NO_NVTOOLSEXT)

if (NOT LIBRARY_ONLY)
    find_package(CUDAToolkit REQUIRED)
    find_package(ZLIB REQUIRED)
    find_package(OpenMP REQUIRED)

    add_executable(align ${ALIGN_SOURCES} main.cu options.cpp sequence_io.cpp)
    target_compile_options(align PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${NVCC_FLAGS_OMP}>)
    target_link_libraries(align ZLIB::ZLIB CUDA::nvToolsExt OpenMP::OpenMP_CXX)
    set_target_properties(align PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    
    add_executable(makedb makedb.cpp sequence_io.cpp dbdata.cpp)
    target_compile_options(makedb PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${NVCC_FLAGS_OMP}>)
    target_link_libraries(makedb ZLIB::ZLIB OpenMP::OpenMP_CXX)
    set_target_properties(makedb PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
