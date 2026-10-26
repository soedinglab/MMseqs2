cmake_minimum_required(VERSION 3.18)

if(ENABLE_HIP)
    set(MARV_GPU_LANGUAGE HIP)
else()
    set(MARV_GPU_LANGUAGE CUDA)
endif()

project(marv LANGUAGES CXX ${MARV_GPU_LANGUAGE})

if(ENABLE_HIP)
    find_package(HIP REQUIRED)
else()
    find_package(CUDAToolkit REQUIRED)
endif()

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

set(MARV_GPU_SOURCES_COMMON
    blosum.cu
    pssmkernels_gapless_instantiation_half2.cu
    pssmkernels_gapless_instantiation_half2_kernelparamzero.cu
    pssmkernels_smithwaterman_instantiation_float.cu
    marv.cu
    main.cu
)

set(MARV_GPU_SOURCES_DPX
    pssmkernels_gapless_instantiation_dpx.cu
    pssmkernels_gapless_instantiation_dpx_kernelparamzero.cu
    pssmkernels_smithwaterman_instantiation_dpx.cu
)

if(ENABLE_HIP)
    set(MARV_GPU_SOURCES ${MARV_GPU_SOURCES_COMMON})
else()
    set(MARV_GPU_SOURCES ${MARV_GPU_SOURCES_COMMON} ${MARV_GPU_SOURCES_DPX})
endif()

set(INT8_GAPLESS_OBJECTS "")
set(INT8_GAPLESS_DEFINITION "")
if(NOT ENABLE_HIP AND CMAKE_CUDA_COMPILER_VERSION VERSION_GREATER_EQUAL "13.2")
    add_library(pssmkernels_gapless_int8_obj OBJECT pssmkernels_gapless_instantiation_int8.cu)
    set_target_properties(pssmkernels_gapless_int8_obj PROPERTIES
        CUDA_ARCHITECTURES "120f-real;121f-real"
    )
    target_compile_options(pssmkernels_gapless_int8_obj PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:${NVCC_FLAGS}>
    )
    target_compile_definitions(pssmkernels_gapless_int8_obj PRIVATE CUDASW_INT8_GAPLESS)
    set(INT8_GAPLESS_OBJECTS $<TARGET_OBJECTS:pssmkernels_gapless_int8_obj>)
    set(INT8_GAPLESS_DEFINITION CUDASW_INT8_GAPLESS)
endif()

function(configure_marv_gpu_target target)
    target_compile_definitions(${target} PRIVATE
        NO_NVTOOLSEXT
    )

    if(ENABLE_HIP)
        target_compile_options(${target} PRIVATE
            $<$<COMPILE_LANGUAGE:HIP>:-Wno-assume>
            $<$<COMPILE_LANGUAGE:HIP>:-Wno-unused-result>
            $<$<COMPILE_LANGUAGE:HIP>:-Wno-uninitialized>
            $<$<COMPILE_LANGUAGE:HIP>:-Wno-unused-value>
        )

        set_target_properties(${target} PROPERTIES
            HIP_SEPARABLE_COMPILATION ON
        )
    else()
        target_include_directories(${target} PRIVATE ${CUDAToolkit_INCLUDE_DIRS} ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
        target_compile_options(${target} PRIVATE
            $<$<COMPILE_LANGUAGE:CUDA>:${NVCC_FLAGS}>
        )

        set_target_properties(${target} PROPERTIES
            CUDA_SEPARABLE_COMPILATION ON
        )
    endif()
endfunction()

add_library(marv ${ALIGN_SOURCES} marv.cu marv.h ${INT8_GAPLESS_OBJECTS})
set_source_files_properties(${MARV_GPU_SOURCES} PROPERTIES LANGUAGE ${MARV_GPU_LANGUAGE})
configure_marv_gpu_target(marv)
target_compile_definitions(marv PRIVATE ${INT8_GAPLESS_DEFINITION})

if(ENABLE_HIP)
    set_target_properties(marv PROPERTIES
        HIP_SEPARABLE_COMPILATION ON
    )
else()
    set_target_properties(marv PROPERTIES
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
)
endif()

if(NOT LIBRARY_ONLY)
    find_package(ZLIB REQUIRED)
    find_package(OpenMP REQUIRED)

    add_executable(align ${ALIGN_SOURCES} main.cu options.cpp sequence_io.cpp ${INT8_GAPLESS_OBJECTS})
    set_source_files_properties(${MARV_GPU_SOURCES} PROPERTIES LANGUAGE ${MARV_GPU_LANGUAGE})
    configure_marv_gpu_target(align)
    target_compile_definitions(align PRIVATE ${INT8_GAPLESS_DEFINITION})
    if(NOT ENABLE_HIP)
        target_compile_options(align PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-fopenmp>)
    endif()
    target_link_libraries(align PRIVATE ZLIB::ZLIB OpenMP::OpenMP_CXX)
    if(NOT ENABLE_HIP)
        if(TARGET CUDA::nvToolsExt)
            target_link_libraries(align PRIVATE CUDA::nvToolsExt)
        endif()
    endif()
    
    add_executable(makedb makedb.cpp sequence_io.cpp dbdata.cpp)
    target_link_libraries(makedb PRIVATE ZLIB::ZLIB OpenMP::OpenMP_CXX)
endif()
