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
    if(archs STREQUAL "native")
        get_filename_component(hip_compiler_dir "${CMAKE_HIP_COMPILER}" DIRECTORY)
        find_program(MARV_AMDGPU_ARCH amdgpu-arch HINTS "${hip_compiler_dir}")
        if(MARV_AMDGPU_ARCH)
            execute_process(COMMAND "${MARV_AMDGPU_ARCH}" OUTPUT_VARIABLE archs OUTPUT_STRIP_TRAILING_WHITESPACE)
            string(REPLACE "\n" ";" archs "${archs}")
        endif()
        if(NOT archs OR archs STREQUAL "native")
            message(FATAL_ERROR "marv: could not detect the native HIP architectures, set CMAKE_HIP_ARCHITECTURES explicitly")
        endif()
    endif()
    # one entry per GPU from native detection
    list(REMOVE_DUPLICATES archs)
    set(${out_var} "${archs}" PARENT_SCOPE)
endfunction()


# Tuned AMD config tables an arch or generic target is compiled into; must match amdConfigTable in aligners/amd_arch.cuh
function(marv_hip_tables_for_arch arch out_var)
    string(REGEX REPLACE ":.*$" "" a "${arch}")
    if(a STREQUAL "gfx11-generic")
        set(tables gfx1100 gfx1150)
    elseif(a MATCHES "^gfx115")
        set(tables gfx1150)
    elseif(a MATCHES "^gfx1[012]")
        set(tables gfx1100)
    elseif(a STREQUAL "gfx90a")
        set(tables gfx90a)
    elseif(a STREQUAL "gfx9-4-generic")
        set(tables gfx942 gfx950)
    elseif(a MATCHES "^gfx95")
        set(tables gfx950)
    elseif(a MATCHES "^gfx94")
        set(tables gfx942)
    elseif(a MATCHES "^gfx9")
        set(tables gfx900)
    else()
        message(FATAL_ERROR "marv: no AMD config table for HIP architecture '${arch}'")
    endif()
    set(${out_var} "${tables}" PARENT_SCOPE)
endfunction()

# Sets <prefix>_TABLES and, per table, <prefix>_<table> to the archs it is compiled for
function(marv_hip_group_by_table archs prefix)
    set(all_tables "")
    foreach(arch IN LISTS archs)
        marv_hip_tables_for_arch("${arch}" tables)
        foreach(table IN LISTS tables)
            list(APPEND all_tables "${table}")
            list(APPEND ${prefix}_${table} "${arch}")
        endforeach()
    endforeach()
    list(REMOVE_DUPLICATES all_tables)
    set(${prefix}_TABLES "${all_tables}" PARENT_SCOPE)
    foreach(table IN LISTS all_tables)
        set(${prefix}_${table} "${${prefix}_${table}}" PARENT_SCOPE)
    endforeach()
endfunction()

function(marv_hip_arch_object_gencode archs out_var)
    list(TRANSFORM archs PREPEND "--offload-arch=" OUTPUT_VARIABLE flags)
    set(${out_var} "${flags}" PARENT_SCOPE)
endfunction()
