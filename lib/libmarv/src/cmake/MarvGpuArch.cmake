# libmarv GPU architecture resolution helpers 
#
# CUDA: turns the requested CMAKE_CUDA_ARCHITECTURES into the gencode flags:
#   - per-architecture tuned kernels, SASS only, family "f" variant for 120/121
#   - thrust/CUB core: compiled only for the major-family base. minor-version forward compatibility covers the rest
#   - PTX on the highest base for forward-JIT to newer GPUs
#
# HIP tuning is dispatched at runtime by gcnArchName.

set(MARV_SUPPORTED_ARCHS "75;80;86;87;89;90;100;103;110;120;121"
    CACHE STRING "libmarv-supported CUDA architectures (used to override all/all-major)")

set(MARV_SUPPORTED_HIP_ARCHS "gfx9-generic;gfx9-4-generic;gfx10-1-generic;gfx10-3-generic;gfx11-generic;gfx12-generic"
    CACHE STRING "Default AMD/HIP architectures")

# Resolve CMAKE_CUDA_ARCHITECTURES into a filtered, de-duplicated integer list.
function(marv_resolve_cuda_architectures out_var)
    set(requested "${CMAKE_CUDA_ARCHITECTURES}")
    if(requested STREQUAL "native")
        set(archs "${CMAKE_CUDA_ARCHITECTURES_NATIVE}")
    elseif(requested STREQUAL "all" OR requested STREQUAL "all-major")
        # CMake's ALL / ALL_MAJOR lists are old
        set(archs "${MARV_SUPPORTED_ARCHS}")
    else()
        set(archs "${requested}")
    endif()

    set(result "")
    foreach(arch IN LISTS archs)
        string(REGEX REPLACE "-(real|virtual)$" "" arch "${arch}")
        if(arch MATCHES "^[0-9]+$" AND NOT arch LESS 75)
            list(APPEND result "${arch}")
        endif()
    endforeach()
    list(REMOVE_DUPLICATES result)
    list(SORT result COMPARE NATURAL)

    if(result STREQUAL "")
        message(FATAL_ERROR "marv: no supported CUDA architectures (>= 75) after resolving CMAKE_CUDA_ARCHITECTURES='${requested}'")
    endif()
    set(${out_var} "${result}" PARENT_SCOPE)
endfunction()

# gencode for one arch's tuned kernels: SASS only, family variant for 120/121.
function(marv_arch_object_gencode arch out_var)
    if(arch STREQUAL "120" OR arch STREQUAL "121")
        set(${out_var} "-gencode=arch=compute_${arch}f,code=sm_${arch}" PARENT_SCOPE)
    else()
        set(${out_var} "-gencode=arch=compute_${arch},code=sm_${arch}" PARENT_SCOPE)
    endif()
endfunction()

# Major-family base for an arch; a base cubin runs on all higher minors
function(marv_arch_base arch out_var)
    math(EXPR base "(${arch} / 10) * 10")
    if(base LESS 75)
        set(base 75)
    endif()
    set(${out_var} "${base}" PARENT_SCOPE)
endfunction()

function(marv_core_architectures archs out_var)
    set(bases "")
    foreach(arch IN LISTS archs)
        marv_arch_base("${arch}" base)
        list(APPEND bases "${base}")
    endforeach()
    list(REMOVE_DUPLICATES bases)
    list(SORT bases COMPARE NATURAL)

    set(result "")
    foreach(base IN LISTS bases)
        list(APPEND result "${base}-real")
    endforeach()
    list(GET bases -1 top)
    list(APPEND result "${top}-virtual")
    set(${out_var} "${result}" PARENT_SCOPE)
endfunction()

function(marv_resolve_hip_architectures out_var)
    if(CMAKE_HIP_ARCHITECTURES)
        set(archs "${CMAKE_HIP_ARCHITECTURES}")
    else()
        set(archs "${MARV_SUPPORTED_HIP_ARCHS}")
    endif()
    set(${out_var} "${archs}" PARENT_SCOPE)
endfunction()


function(marv_hip_arch_object_gencode arch out_var)
    set(${out_var} "--offload-arch=${arch}" PARENT_SCOPE)
endfunction()
