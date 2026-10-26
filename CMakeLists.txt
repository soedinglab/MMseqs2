cmake_minimum_required(VERSION 3.16)
project(carquet VERSION 0.7.2 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# macOS deployment target - respect user-specified value, only set default if not provided
if(APPLE)
    # Only set a default deployment target if one wasn't already specified
    if(NOT CMAKE_OSX_DEPLOYMENT_TARGET)
        # Get the current macOS version and use it as deployment target
        execute_process(
            COMMAND sw_vers -productVersion
            OUTPUT_VARIABLE MACOS_VERSION
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE SW_VERS_RESULT
        )
        if(SW_VERS_RESULT EQUAL 0 AND MACOS_VERSION)
            # Extract major version (e.g., "15" from "15.2")
            string(REGEX MATCH "^[0-9]+" MACOS_MAJOR_VERSION "${MACOS_VERSION}")
            if(MACOS_MAJOR_VERSION)
                set(CMAKE_OSX_DEPLOYMENT_TARGET "${MACOS_MAJOR_VERSION}.0" CACHE STRING "Minimum macOS version")
            endif()
        endif()
    endif()

    # Suppress linker warnings about library version mismatches from Homebrew LLVM
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -w")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -w")
endif()

# Build options
option(CARQUET_BUILD_DEV "Build all development targets (tests, examples, benchmarks)" OFF)
option(CARQUET_BUILD_TESTS "Build tests" ${CARQUET_BUILD_DEV})
option(CARQUET_BUILD_EXAMPLES "Build example programs" ${CARQUET_BUILD_DEV})
option(CARQUET_BUILD_BENCHMARKS "Build benchmark and profiling programs" ${CARQUET_BUILD_DEV})
option(CARQUET_BUILD_ARROW_CPP_BENCHMARK "Build the optional Arrow C++ comparison benchmark" OFF)
option(CARQUET_BUILD_INTEROP "Build interoperability test programs" ${CARQUET_BUILD_DEV})
option(CARQUET_BUILD_CLI "Build carquet CLI tool" ON)
option(CARQUET_BUILD_SHARED "Build shared library" OFF)
if(CARQUET_BUILD_SHARED)
    # The fetched codec archives (zlib, zstd, lz4, brotli) are linked INTO
    # libcarquet.so, so their objects must be position independent; without
    # this the Linux link fails with "relocation ... can not be used when making
    # a shared object; recompile with -fPIC" (seen on the manylinux wheel build).
    # Setting the variable here initialises POSITION_INDEPENDENT_CODE for every
    # target created below, including the FetchContent subprojects.
    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()
option(CARQUET_ENABLE_SSE "Enable SSE optimizations" ON)
option(CARQUET_ENABLE_AVX "Enable AVX optimizations" ON)
option(CARQUET_ENABLE_AVX2 "Enable AVX2 optimizations" ON)
option(CARQUET_ENABLE_AVX512 "Enable AVX-512 optimizations" ON)
option(CARQUET_ENABLE_NEON "Enable NEON optimizations" ON)
option(CARQUET_ENABLE_SVE "Enable SVE optimizations" OFF)
set(CARQUET_ARROW_CPP_ROOT "" CACHE PATH
    "Optional Arrow/Parquet install prefix used only for benchmark_arrow_cpp")

# External compression libraries
include(FetchContent)

# When ON, ignore any system/conda zstd/zlib/lz4 and always fetch + statically
# bundle them into libcarquet. Produces a self-contained static archive suitable
# for prebuilt release artifacts (see .github/workflows/release.yml).
option(CARQUET_BUNDLE_DEPS "Always fetch and static-bundle zstd/zlib/lz4/brotli (portable release artifacts)" OFF)
# BROTLI (Parquet codec 4). Resolved like the other three (target / system /
# fetch); OFF leaves the slot to carquet_register_codec() and drops the
# dependency (libbrotli's encoder is noticeably larger than lz4 or zstd).
option(CARQUET_ENABLE_BROTLI "Build the BROTLI codec (needs libbrotli; found or fetched)" ON)

# Which zstd / zlib / lz4 carquet links, in order of precedence:
#
#   1. CARQUET_ZSTD_TARGET / CARQUET_ZLIB_TARGET / CARQUET_LZ4_TARGET, when
#      set: the name of a target of the enclosing project. It is linked as is,
#      with no search and no fetch, so this is how a project that builds
#      several variants (a static and a shared zlib, say) picks the one
#      carquet gets.
#   2. A system copy, unless CARQUET_BUNDLE_DEPS is ON.
#   3. A copy fetched and built here. If the enclosing project already builds
#      the library under the target name that fetch would use (libzstd_static,
#      zlibstatic, lz4_static), that target is reused instead: a second copy
#      would only collide with it.
#
# A target taken from 1 or 3 may be imported or built in-tree, static or
# shared; it must carry its include directories as a usage requirement.
foreach(_carquet_dep ZSTD ZLIB LZ4 BROTLI)
    # Not declared when already defined: before CMake 3.21 (CMP0126) creating a
    # cache entry discards a same-named normal variable, which is exactly how
    # an enclosing project sets these before add_subdirectory(carquet).
    if(NOT DEFINED CARQUET_${_carquet_dep}_TARGET)
        set(CARQUET_${_carquet_dep}_TARGET "" CACHE STRING
            "Existing CMake target to link as ${_carquet_dep} instead of searching for or fetching it")
    endif()
    if(CARQUET_${_carquet_dep}_TARGET AND NOT TARGET ${CARQUET_${_carquet_dep}_TARGET})
        message(FATAL_ERROR
            "CARQUET_${_carquet_dep}_TARGET is set to \"${CARQUET_${_carquet_dep}_TARGET}\", "
            "which is not a CMake target here. It must name a target that exists "
            "before carquet is added; leave it empty to let carquet find the library.")
    endif()
endforeach()

# Try to find system zstd first, fetch if not found
if(CARQUET_ZSTD_TARGET)
    message(STATUS "Using zstd target ${CARQUET_ZSTD_TARGET} (CARQUET_ZSTD_TARGET)")
    set(ZSTD_LIBRARIES ${CARQUET_ZSTD_TARGET})
    set(ZSTD_INCLUDE_DIRS "")
elseif(CARQUET_BUNDLE_DEPS)
    message(STATUS "CARQUET_BUNDLE_DEPS=ON — forcing bundled static zstd")
    # Nothing is searched for here, so a ZSTD_FOUND already in scope was set by
    # the enclosing project (or by an earlier configure) and must not count.
    set(ZSTD_FOUND FALSE)
else()
    # Method 1: Try modern CMake config (works best on Windows/vcpkg/conda)
    find_package(zstd CONFIG QUIET)
    if(zstd_FOUND)
        message(STATUS "Found zstd via CMake config")
        # When building a shared library, prefer shared deps to avoid -fPIC issues
        # with static libraries that weren't compiled with -fPIC.
        if(CARQUET_BUILD_SHARED AND TARGET zstd::libzstd_shared)
            set(ZSTD_LIBRARIES zstd::libzstd_shared)
        elseif(TARGET zstd::libzstd_static)
            set(ZSTD_LIBRARIES zstd::libzstd_static)
        elseif(TARGET zstd::libzstd_shared)
            set(ZSTD_LIBRARIES zstd::libzstd_shared)
        endif()
        set(ZSTD_FOUND TRUE)
    else()
        # Method 2: Try pkg-config (Linux/macOS)
        find_package(PkgConfig QUIET)
        if(PkgConfig_FOUND)
            pkg_check_modules(ZSTD QUIET libzstd)
        endif()
    endif()
endif()

# <name>_populated: the enclosing project already fetched this content under
# the same FetchContent name, in which case the fetch below simply reuses it.
FetchContent_GetProperties(zstd POPULATED _carquet_zstd_populated)
if(CARQUET_ZSTD_TARGET)
    # Named explicitly: nothing to resolve.
elseif(ZSTD_FOUND)
    message(STATUS "Found system libzstd")
elseif(TARGET libzstd_static AND NOT _carquet_zstd_populated)
    message(STATUS "Reusing the enclosing project's zstd target: libzstd_static")
    set(ZSTD_LIBRARIES libzstd_static)
    set(ZSTD_INCLUDE_DIRS "")
else()
    # Method 3: Fetch from source
    message(STATUS "Fetching zstd...")
    FetchContent_Declare(
        zstd
        GIT_REPOSITORY https://github.com/facebook/zstd.git
        GIT_TAG        v1.5.6
        GIT_SHALLOW    TRUE
        GIT_PROGRESS   FALSE
        SOURCE_SUBDIR  build/cmake
    )
    set(ZSTD_BUILD_PROGRAMS OFF CACHE BOOL "" FORCE)
    set(ZSTD_BUILD_TESTS OFF CACHE BOOL "" FORCE)
    set(ZSTD_BUILD_SHARED OFF CACHE BOOL "" FORCE)
    set(ZSTD_BUILD_STATIC ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(zstd)
    set(ZSTD_LIBRARIES libzstd_static)
    set(ZSTD_INCLUDE_DIRS ${zstd_SOURCE_DIR}/lib)
endif()

# BROTLI: target / system (pkg-config libbrotlienc+libbrotlidec, or vcpkg's
# unofficial-brotli config) / fetch google/brotli. BROTLI_LIBRARIES empty means
# "built without brotli": the codec then reports UNSUPPORTED_CODEC.
set(BROTLI_LIBRARIES "")
set(BROTLI_INCLUDE_DIRS "")
set(BROTLI_FOUND FALSE)
if(CARQUET_ENABLE_BROTLI)
    if(CARQUET_BROTLI_TARGET)
        message(STATUS "Using brotli target ${CARQUET_BROTLI_TARGET} (CARQUET_BROTLI_TARGET)")
        set(BROTLI_LIBRARIES ${CARQUET_BROTLI_TARGET})
        set(BROTLI_FOUND TRUE)
    elseif(NOT CARQUET_BUNDLE_DEPS)
        find_package(unofficial-brotli CONFIG QUIET)
        if(unofficial-brotli_FOUND AND TARGET unofficial::brotli::brotlienc
           AND TARGET unofficial::brotli::brotlidec)
            message(STATUS "Found brotli via CMake config (unofficial-brotli)")
            set(BROTLI_LIBRARIES unofficial::brotli::brotlienc unofficial::brotli::brotlidec)
            set(BROTLI_FOUND TRUE)
        else()
            find_package(PkgConfig QUIET)
            if(PkgConfig_FOUND)
                pkg_check_modules(BROTLI QUIET libbrotlienc libbrotlidec)
                if(BROTLI_FOUND)
                    message(STATUS "Found system libbrotli via pkg-config")
                    if(BROTLI_LINK_LIBRARIES)
                        set(BROTLI_LIBRARIES ${BROTLI_LINK_LIBRARIES})
                    endif()
                endif()
            endif()
        endif()
    endif()
    FetchContent_GetProperties(brotli POPULATED _carquet_brotli_populated)
    if(BROTLI_FOUND)
        # resolved above
    elseif(TARGET brotlienc AND TARGET brotlidec AND NOT _carquet_brotli_populated)
        message(STATUS "Reusing the enclosing project's brotli targets: brotlienc / brotlidec")
        set(BROTLI_LIBRARIES brotlienc brotlidec)
    else()
        message(STATUS "Fetching brotli...")
        FetchContent_Declare(
            brotli
            GIT_REPOSITORY https://github.com/google/brotli.git
            GIT_TAG        v1.1.0
            GIT_SHALLOW    TRUE
            GIT_PROGRESS   FALSE
        )
        set(BROTLI_DISABLE_TESTS ON CACHE BOOL "" FORCE)
        # brotli follows BUILD_SHARED_LIBS; shadow it with a normal variable for
        # the duration of the add so the enclosing project's setting is untouched.
        # EXCLUDE_FROM_ALL: brotli 1.1.0 always declares its `brotli` CLI target;
        # only the three libraries carquet links get built this way.
        set(_carquet_saved_bsl "${BUILD_SHARED_LIBS}")
        set(BUILD_SHARED_LIBS OFF)
        FetchContent_GetProperties(brotli)
        if(NOT brotli_POPULATED)
            FetchContent_Populate(brotli)
            add_subdirectory(${brotli_SOURCE_DIR} ${brotli_BINARY_DIR} EXCLUDE_FROM_ALL)
        endif()
        set(BUILD_SHARED_LIBS "${_carquet_saved_bsl}")
        set(BROTLI_LIBRARIES brotlienc brotlidec brotlicommon)
        set(BROTLI_INCLUDE_DIRS ${brotli_SOURCE_DIR}/c/include)
        if(NOT MSVC)
            # brotli's CMake does not add -fPIC to its static objects; a shared
            # libcarquet needs it.
            foreach(_bt brotlienc brotlidec brotlicommon)
                set_target_properties(${_bt} PROPERTIES POSITION_INDEPENDENT_CODE ON)
            endforeach()
        endif()
    endif()
else()
    message(STATUS "BROTLI codec disabled (CARQUET_ENABLE_BROTLI=OFF)")
endif()

# Try to find system zlib first, fetch if not found
if(CARQUET_ZLIB_TARGET)
    message(STATUS "Using zlib target ${CARQUET_ZLIB_TARGET} (CARQUET_ZLIB_TARGET)")
    set(ZLIB_LIBRARIES ${CARQUET_ZLIB_TARGET})
    set(ZLIB_INCLUDE_DIRS "")
elseif(NOT CARQUET_BUNDLE_DEPS)
    find_package(ZLIB QUIET)
else()
    set(ZLIB_FOUND FALSE)  # see ZSTD_FOUND above
endif()
FetchContent_GetProperties(zlib POPULATED _carquet_zlib_populated)
if(CARQUET_ZLIB_TARGET)
    # Named explicitly: nothing to resolve.
elseif(ZLIB_FOUND)
    message(STATUS "Found system zlib: ${ZLIB_VERSION_STRING}")
elseif(TARGET zlibstatic AND NOT _carquet_zlib_populated)
    message(STATUS "Reusing the enclosing project's zlib target: zlibstatic")
    set(ZLIB_LIBRARIES zlibstatic)
    set(ZLIB_INCLUDE_DIRS "")
else()
    message(STATUS "Fetching zlib...")
    FetchContent_Declare(
        zlib
        GIT_REPOSITORY https://github.com/madler/zlib.git
        GIT_TAG        v1.3.1
        GIT_SHALLOW    TRUE
        GIT_PROGRESS   FALSE
    )
    set(ZLIB_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(zlib)
    set(ZLIB_LIBRARIES zlibstatic)
    set(ZLIB_INCLUDE_DIRS ${zlib_SOURCE_DIR} ${zlib_BINARY_DIR})
endif()

# Try to find system lz4 first, fetch if not found
# Prefer CMake config (gives proper targets with full paths on all platforms)
if(CARQUET_LZ4_TARGET)
    message(STATUS "Using lz4 target ${CARQUET_LZ4_TARGET} (CARQUET_LZ4_TARGET)")
    set(LZ4_LIBRARIES ${CARQUET_LZ4_TARGET})
    set(LZ4_INCLUDE_DIRS "")
    set(LZ4_FOUND TRUE)
elseif(NOT CARQUET_BUNDLE_DEPS)
    find_package(lz4 CONFIG QUIET)
else()
    set(lz4_FOUND FALSE)  # see ZSTD_FOUND above
    set(LZ4_FOUND FALSE)
endif()
# lz4_FOUND (lower case) is what find_package(lz4 CONFIG) sets: the variable is
# named after the package exactly as written in the call. LZ4_FOUND (upper
# case) is carquet's own "resolved to something linkable" flag, set below and
# by pkg-config. The enclosing project may have run its own find_package(lz4),
# leaving lz4_FOUND set in this scope; a target named explicitly still wins.
if(lz4_FOUND AND NOT CARQUET_LZ4_TARGET)
    message(STATUS "Found lz4 via CMake config")
    if(CARQUET_BUILD_SHARED AND TARGET LZ4::lz4_shared)
        set(LZ4_LIBRARIES LZ4::lz4_shared)
    elseif(TARGET LZ4::lz4_static)
        set(LZ4_LIBRARIES LZ4::lz4_static)
    elseif(TARGET LZ4::lz4_shared)
        set(LZ4_LIBRARIES LZ4::lz4_shared)
    elseif(TARGET lz4::lz4)
        set(LZ4_LIBRARIES lz4::lz4)
    endif()
    # Only mark as found if we actually resolved a target
    if(LZ4_LIBRARIES)
        set(LZ4_FOUND TRUE)
    else()
        message(STATUS "lz4 CMake config found but no known target — trying pkg-config")
        set(LZ4_FOUND FALSE)
    endif()
endif()

if(NOT LZ4_FOUND AND NOT CARQUET_BUNDLE_DEPS)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(LZ4 QUIET liblz4)
        if(LZ4_FOUND)
            message(STATUS "Found system liblz4 via pkg-config")
            # Use LINK_LIBRARIES (full paths) instead of LIBRARIES (bare names)
            # to avoid MSVC linker failing to locate the .lib file
            if(LZ4_LINK_LIBRARIES)
                set(LZ4_LIBRARIES ${LZ4_LINK_LIBRARIES})
            endif()
        endif()
    endif()
endif()

FetchContent_GetProperties(lz4 POPULATED _carquet_lz4_populated)
if(NOT LZ4_FOUND AND TARGET lz4_static AND NOT _carquet_lz4_populated)
    message(STATUS "Reusing the enclosing project's lz4 target: lz4_static")
    set(LZ4_LIBRARIES lz4_static)
    set(LZ4_INCLUDE_DIRS "")
    set(LZ4_FOUND TRUE)
endif()

if(NOT LZ4_FOUND)
    message(STATUS "Fetching lz4...")
    FetchContent_Declare(
        lz4
        GIT_REPOSITORY https://github.com/lz4/lz4.git
        GIT_TAG        v1.10.0
        GIT_SHALLOW    TRUE
        GIT_PROGRESS   FALSE
        SOURCE_SUBDIR  build/cmake
    )
    set(LZ4_BUILD_CLI OFF CACHE BOOL "" FORCE)
    set(LZ4_BUILD_LEGACY_LZ4C OFF CACHE BOOL "" FORCE)
    # Save and restore BUILD_SHARED_LIBS to avoid corrupting cache for rest of build
    set(_saved_BSL ${BUILD_SHARED_LIBS})
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(BUILD_STATIC_LIBS ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(lz4)
    set(BUILD_SHARED_LIBS ${_saved_BSL} CACHE BOOL "" FORCE)
    set(LZ4_LIBRARIES lz4_static)
    set(LZ4_INCLUDE_DIRS ${lz4_SOURCE_DIR}/lib)
endif()

# Compiler detection
if(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
    set(CARQUET_COMPILER_GCC_LIKE TRUE)
elseif(CMAKE_C_COMPILER_ID MATCHES "MSVC")
    set(CARQUET_COMPILER_MSVC TRUE)
endif()

# Architecture detection
include(CheckCSourceCompiles)
check_c_source_compiles("
    #if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    int main() { return 0; }
    #else
    #error Not x86
    #endif
" CARQUET_ARCH_X86)

check_c_source_compiles("
    #if defined(__aarch64__) || defined(__arm64__) || defined(_M_ARM64) || defined(__arm__) || defined(_M_ARM)
    int main() { return 0; }
    #else
    #error Not ARM
    #endif
" CARQUET_ARCH_ARM)

# Compiler flags
if(CARQUET_COMPILER_GCC_LIKE)
    add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-parameter)

    # Debug/Release flags
    set(CMAKE_C_FLAGS_DEBUG "-g -O0 -DCARQUET_DEBUG")
    set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -flto=auto")
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG")

    # When using LTO, ensure ar/ranlib match the compiler toolchain
    # (Apple's ranlib can't parse LLVM bitcode from Homebrew LLVM)
    if(CMAKE_C_COMPILER MATCHES "llvm")
        get_filename_component(_LLVM_BIN_DIR "${CMAKE_C_COMPILER}" DIRECTORY)
        find_program(_LLVM_AR llvm-ar HINTS "${_LLVM_BIN_DIR}" NO_DEFAULT_PATH)
        find_program(_LLVM_RANLIB llvm-ranlib HINTS "${_LLVM_BIN_DIR}" NO_DEFAULT_PATH)
        if(_LLVM_AR)
            set(CMAKE_AR "${_LLVM_AR}" CACHE FILEPATH "Archiver" FORCE)
        endif()
        if(_LLVM_RANLIB)
            set(CMAKE_RANLIB "${_LLVM_RANLIB}" CACHE FILEPATH "Ranlib" FORCE)
        endif()
    endif()
elseif(CARQUET_COMPILER_MSVC)
    # MSVC warning levels and options
    add_compile_options(/W3 /wd4996 /wd4244 /wd4267)

    # Debug/Release flags
    set(CMAKE_C_FLAGS_DEBUG "/Zi /Od /DCARQUET_DEBUG")
    set(CMAKE_C_FLAGS_RELEASE "/O2 /DNDEBUG")
    set(CMAKE_C_FLAGS_RELWITHDEBINFO "/O2 /Zi /DNDEBUG")
endif()

# Optional native architecture optimization
option(CARQUET_NATIVE_ARCH "Build with -march=native for maximum performance (produces a HOST-ONLY binary)" OFF)
if(CARQUET_NATIVE_ARCH AND CARQUET_COMPILER_GCC_LIKE)
    # -march=native applies to the whole library, including the scalar fallbacks
    # and the dispatch layer, so the compiler may auto-vectorize "baseline" code
    # with whatever ISA the build host supports. The resulting binary is NOT
    # portable to older CPUs and the runtime SIMD dispatch no longer guarantees
    # it: run it only on hardware matching (or newer than) the build host. Use
    # it for local benchmarking, not for distributed/packaged builds.
    message(STATUS "CARQUET_NATIVE_ARCH=ON: building a host-only binary (-march=native); not portable to older CPUs")
    add_compile_options(-march=native)
endif()

# Source files
set(CARQUET_CORE_SOURCES
    src/core/arena.c
    src/core/allocator.c
    src/core/buffer.c
    src/core/bitpack.c
    src/core/endian.c
    src/core/error.c
    src/core/geo_wkb.c
)

set(CARQUET_THRIFT_SOURCES
    src/thrift/thrift_decode.c
    src/thrift/thrift_encode.c
    src/thrift/parquet_types.c
)

set(CARQUET_ENCODING_SOURCES
    src/encoding/plain.c
    src/encoding/rle.c
    src/encoding/delta.c
    src/encoding/delta_length.c
    src/encoding/delta_strings.c
    src/encoding/dictionary.c
    src/encoding/byte_stream_split.c
)

set(CARQUET_COMPRESSION_SOURCES
    src/compression/lz4.c
    src/compression/snappy.c
    src/compression/zstd.c
    src/compression/gzip.c
    src/compression/brotli.c
    src/compression/custom.c
)

# Snappy uses SSSE3 (pshufb) for pattern extension on x86.
# SSSE3 is baseline since Core 2 (2006) / x86-64-v2.
if(CARQUET_ARCH_X86 AND CARQUET_ENABLE_SSE)
    if(CARQUET_COMPILER_GCC_LIKE)
        set_source_files_properties(src/compression/snappy.c PROPERTIES
            COMPILE_FLAGS "-mssse3")
    endif()
endif()

set(CARQUET_SIMD_SOURCES
    src/simd/detect.c
    src/simd/dispatch.c
)

# Runtime CPU feature detection must stay correct under mixed LLVM/macOS toolchains.
# Excluding this tiny unit from LTO avoids miscompiled feature flags on Apple Silicon.
if(CARQUET_COMPILER_GCC_LIKE)
    set_source_files_properties(src/simd/detect.c PROPERTIES COMPILE_FLAGS "-fno-lto")
endif()

set(CARQUET_READER_SOURCES
    src/reader/file_reader.c
    src/reader/arrow_schema_read.c
    src/reader/arrow_c_export.c
    src/reader/arrow_c_read.c
    src/reader/row_group_reader.c
    src/reader/column_reader.c
    src/reader/page_reader.c
    src/reader/batch_reader.c
    src/reader/page_filter.c
    src/reader/statistics.c
    src/reader/mmap_reader.c
    src/reader/worker_pool.c
)

set(CARQUET_WRITER_SOURCES
    src/writer/file_writer.c
    src/writer/row_group_writer.c
    src/writer/column_writer.c
    src/writer/page_writer.c
    src/writer/arrow_schema.c
    src/writer/arrow_c_import.c
)

set(CARQUET_METADATA_SOURCES
    src/metadata/schema.c
    src/metadata/bloom_filter.c
    src/metadata/page_index.c
)

set(CARQUET_UTIL_SOURCES
    src/util/crc32.c
    src/util/xxhash.c
)

# x86 SIMD sources — only compile files for enabled instruction sets
if(CARQUET_ARCH_X86)
    # Auto-detect: disable AVX-512 if compiler lacks _mm512_set_epi8 (GCC < 9)
    if(CARQUET_ENABLE_AVX512 AND CARQUET_COMPILER_GCC_LIKE)
        include(CheckCSourceCompiles)
        set(CMAKE_REQUIRED_FLAGS "-mavx512f -mavx512bw -mavx512vl")
        check_c_source_compiles("
            #include <immintrin.h>
            int main(void) {
                __m512i v = _mm512_set_epi8(
                    0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                    16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
                    32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,
                    48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63);
                (void)v; return 0;
            }
        " CARQUET_HAS_AVX512_INTRINSICS)
        unset(CMAKE_REQUIRED_FLAGS)
        if(NOT CARQUET_HAS_AVX512_INTRINSICS)
            message(WARNING "Compiler lacks full AVX-512 intrinsic support — disabling AVX-512")
            set(CARQUET_ENABLE_AVX512 OFF)
        endif()
    endif()

    set(CARQUET_SIMD_X86_SOURCES "")

    # ISA-specific TUs are compiled with per-file -m<isa> flags. Under LTO the
    # code generator can inline these functions into baseline callers, emitting
    # AVX/AVX-512 instructions on paths that run before the runtime CPU-feature
    # dispatch check -> SIGILL on older CPUs. Compile them with -fno-lto (as
    # detect.c already is) so ISA-specific code stays in its own object and is
    # only reached through the dispatch table's function pointers.
    if(CARQUET_ENABLE_SSE)
        list(APPEND CARQUET_SIMD_X86_SOURCES src/simd/x86/sse_ops.c)
        if(CARQUET_COMPILER_GCC_LIKE)
            set_source_files_properties(src/simd/x86/sse_ops.c PROPERTIES
                COMPILE_FLAGS "-msse4.2 -fno-lto")
        endif()
    endif()

    if(CARQUET_ENABLE_AVX)
        list(APPEND CARQUET_SIMD_X86_SOURCES src/simd/x86/avx_ops.c)
        if(CARQUET_COMPILER_GCC_LIKE)
            set_source_files_properties(src/simd/x86/avx_ops.c PROPERTIES
                COMPILE_FLAGS "-mavx -msse4.2 -fno-lto")
        elseif(CARQUET_COMPILER_MSVC)
            set_source_files_properties(src/simd/x86/avx_ops.c PROPERTIES
                COMPILE_FLAGS "/arch:AVX")
        endif()
    endif()

    if(CARQUET_ENABLE_AVX2)
        list(APPEND CARQUET_SIMD_X86_SOURCES src/simd/x86/avx2_ops.c)
        if(CARQUET_COMPILER_GCC_LIKE)
            set_source_files_properties(src/simd/x86/avx2_ops.c PROPERTIES
                COMPILE_FLAGS "-mavx2 -mbmi2 -fno-lto")
        elseif(CARQUET_COMPILER_MSVC)
            set_source_files_properties(src/simd/x86/avx2_ops.c PROPERTIES
                COMPILE_FLAGS "/arch:AVX2")
        endif()
    endif()

    if(CARQUET_ENABLE_AVX512)
        list(APPEND CARQUET_SIMD_X86_SOURCES src/simd/x86/avx512_ops.c)
        if(CARQUET_COMPILER_GCC_LIKE)
            set_source_files_properties(src/simd/x86/avx512_ops.c PROPERTIES
                COMPILE_FLAGS "-mavx512f -mavx512bw -mavx512vl -fno-lto")
        elseif(CARQUET_COMPILER_MSVC)
            set_source_files_properties(src/simd/x86/avx512_ops.c PROPERTIES
                COMPILE_FLAGS "/arch:AVX512")
        endif()
    endif()
endif()

# ARM SIMD sources — only compile files for enabled instruction sets
if(CARQUET_ARCH_ARM)
    set(CARQUET_SIMD_ARM_SOURCES "")

    if(CARQUET_ENABLE_NEON)
        list(APPEND CARQUET_SIMD_ARM_SOURCES src/simd/arm/neon_ops.c)
    endif()

    if(CARQUET_ENABLE_SVE)
        list(APPEND CARQUET_SIMD_ARM_SOURCES src/simd/arm/sve_ops.c)
        if(CARQUET_COMPILER_GCC_LIKE)
            set_source_files_properties(src/simd/arm/sve_ops.c PROPERTIES
                COMPILE_FLAGS "-march=armv8-a+sve")
        endif()
    endif()
endif()

# OpenMP support for parallel column reading
find_package(OpenMP)
if(OpenMP_C_FOUND)
    message(STATUS "OpenMP found - enabling parallel column reading")
endif()

# Combine all sources
set(CARQUET_SOURCES
    ${CARQUET_CORE_SOURCES}
    ${CARQUET_THRIFT_SOURCES}
    ${CARQUET_ENCODING_SOURCES}
    ${CARQUET_COMPRESSION_SOURCES}
    ${CARQUET_SIMD_SOURCES}
    ${CARQUET_READER_SOURCES}
    ${CARQUET_WRITER_SOURCES}
    ${CARQUET_METADATA_SOURCES}
    ${CARQUET_UTIL_SOURCES}
)

if(CARQUET_ARCH_X86)
    list(APPEND CARQUET_SOURCES ${CARQUET_SIMD_X86_SOURCES})
endif()
if(CARQUET_ARCH_ARM)
    list(APPEND CARQUET_SOURCES ${CARQUET_SIMD_ARM_SOURCES})
endif()

# Create library
if(CARQUET_BUILD_SHARED)
    add_library(carquet SHARED ${CARQUET_SOURCES})
    set_target_properties(carquet PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
    target_compile_definitions(carquet
        PUBLIC CARQUET_BUILD_SHARED
        PRIVATE CARQUET_BUILDING_DLL
    )
else()
    add_library(carquet STATIC ${CARQUET_SOURCES})
endif()

target_include_directories(carquet
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${ZSTD_INCLUDE_DIRS}
        ${ZLIB_INCLUDE_DIRS}
        ${LZ4_INCLUDE_DIRS}
        ${BROTLI_INCLUDE_DIRS}
)
if(BROTLI_LIBRARIES)
    target_compile_definitions(carquet PRIVATE CARQUET_HAVE_BROTLI=1)
endif()

# Link external compression libraries
set(CARQUET_COMPRESSION_LIBS ${ZSTD_LIBRARIES} ${LZ4_LIBRARIES} ${BROTLI_LIBRARIES})

function(carquet_private_link_libraries_for_export out_var)
    set(_result)
    foreach(_lib IN LISTS ARGN)
        if(TARGET ${_lib})
            get_target_property(_imported ${_lib} IMPORTED)
            if(NOT _imported AND CARQUET_BUILD_SHARED)
                # Shared carquet bakes the fetched dep into the .so, so consumers
                # of the installed library must NOT link it — keep it build-tree only.
                list(APPEND _result "$<BUILD_INTERFACE:${_lib}>")
            else()
                # Static carquet: the consumer must link the bundled dep archive.
                # These fetched targets are added to carquetTargets below, so a
                # plain reference exports as carquet::<dep> in the install interface.
                list(APPEND _result ${_lib})
            endif()
        else()
            list(APPEND _result ${_lib})
        endif()
    endforeach()
    set(${out_var} ${_result} PARENT_SCOPE)
endfunction()

carquet_private_link_libraries_for_export(
    CARQUET_PRIVATE_COMPRESSION_LIBS
    ${CARQUET_COMPRESSION_LIBS})

# ZLIB_FOUND can also be inherited from a find_package(ZLIB) in the enclosing
# project, so a target named explicitly is checked first.
if(ZLIB_FOUND AND NOT CARQUET_ZLIB_TARGET)
    target_link_libraries(carquet PRIVATE ${CARQUET_PRIVATE_COMPRESSION_LIBS} ZLIB::ZLIB)
else()
    carquet_private_link_libraries_for_export(
        CARQUET_PRIVATE_ZLIB_LIBS
        ${ZLIB_LIBRARIES})
    target_link_libraries(carquet PRIVATE
        ${CARQUET_PRIVATE_COMPRESSION_LIBS}
        ${CARQUET_PRIVATE_ZLIB_LIBS})
endif()

# Link libm for math functions (bloom_filter.c uses log())
# PUBLIC so that consumers of the shared library also link -lm; on some toolchains
# (RHEL with --no-copy-dt-needed-entries) PRIVATE alone is not sufficient.
if(UNIX)
    target_link_libraries(carquet PUBLIC m)
endif()

# Link pthreads for worker pool (POSIX)
if(UNIX)
    find_package(Threads REQUIRED)
    target_link_libraries(carquet PRIVATE Threads::Threads)
endif()

# Link OpenMP for parallel column reading
if(OpenMP_C_FOUND)
    target_link_libraries(carquet PRIVATE OpenMP::OpenMP_C)
endif()

# Architecture definitions
if(CARQUET_ARCH_X86)
    target_compile_definitions(carquet PRIVATE CARQUET_ARCH_X86)
    if(CARQUET_ENABLE_SSE)
        target_compile_definitions(carquet PRIVATE CARQUET_ENABLE_SSE)
    endif()
    if(CARQUET_ENABLE_AVX)
        target_compile_definitions(carquet PRIVATE CARQUET_ENABLE_AVX)
    endif()
    if(CARQUET_ENABLE_AVX2)
        target_compile_definitions(carquet PRIVATE CARQUET_ENABLE_AVX2)
    endif()
    if(CARQUET_ENABLE_AVX512)
        target_compile_definitions(carquet PRIVATE CARQUET_ENABLE_AVX512)
    endif()
endif()

if(CARQUET_ARCH_ARM)
    target_compile_definitions(carquet PRIVATE CARQUET_ARCH_ARM)
    if(CARQUET_ENABLE_NEON)
        target_compile_definitions(carquet PRIVATE CARQUET_ENABLE_NEON)
    endif()
    if(CARQUET_ENABLE_SVE)
        target_compile_definitions(carquet PRIVATE CARQUET_ENABLE_SVE)
    endif()
endif()

# Tests
if(CARQUET_BUILD_TESTS)
    enable_testing()

    # Helper function to configure test targets
    function(add_carquet_test TEST_NAME TEST_SOURCE)
        add_executable(${TEST_NAME} ${TEST_SOURCE})
        target_link_libraries(${TEST_NAME} PRIVATE carquet)
        target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
        # Tests use assert() for validation — keep assertions enabled in all build types
        if(MSVC)
            target_compile_options(${TEST_NAME} PRIVATE /UNDEBUG)
        else()
            target_compile_options(${TEST_NAME} PRIVATE -UNDEBUG)
        endif()
        if(CARQUET_COMPILER_GCC_LIKE)
            target_compile_options(${TEST_NAME} PRIVATE -Wno-unused-result)
        endif()
        add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    endfunction()

    add_carquet_test(test_core tests/test_core.c)
    add_carquet_test(test_fixes_0_7_2 tests/test_fixes_0_7_2.c)
    add_carquet_test(test_thrift tests/test_thrift.c)
    add_carquet_test(test_encodings tests/test_encodings.c)
    add_carquet_test(test_reader tests/test_reader.c)
    add_carquet_test(test_encodings_extended tests/test_encodings_extended.c)
    add_carquet_test(test_compression tests/test_compression.c)
    add_carquet_test(test_utils tests/test_utils.c)
    add_carquet_test(test_production tests/test_production.c)
    add_carquet_test(test_maturity tests/test_maturity.c)
    add_carquet_test(test_large_schema tests/test_large_schema.c)

    # Edge case tests
    add_carquet_test(test_edge_encodings tests/test_edge_encodings.c)
    add_carquet_test(test_edge_compression tests/test_edge_compression.c)
    add_carquet_test(test_edge_malformed tests/test_edge_malformed.c)
    add_carquet_test(test_edge_boundaries tests/test_edge_boundaries.c)
    add_carquet_test(test_edge_io tests/test_edge_io.c)

    # Nested schema tests
    add_carquet_test(test_nested tests/test_nested.c)

    # Bloom filter and page index tests
    add_carquet_test(test_bloom_page_index tests/test_bloom_page_index.c)

    # Memory mapping tests
    add_carquet_test(test_mmap tests/test_mmap.c)

    # Advanced API tests (bloom filter, page index, kv metadata, etc.)
    add_carquet_test(test_advanced_api tests/test_advanced_api.c)

    # carquet-reads-carquet roundtrip for every opt-in encoding
    add_carquet_test(test_encoding_roundtrip tests/test_encoding_roundtrip.c)

    # INT96 writing, opt-in Data Page V2, ARROW:schema metadata
    add_carquet_test(test_writer_extensions tests/test_writer_extensions.c)

    # Arrow C Data Interface export/import bridge
    add_carquet_test(test_arrow_c_data tests/test_arrow_c_data.c)
    add_carquet_test(test_arrow_nested tests/test_arrow_nested.c)
    add_carquet_test(test_p2_conformance tests/test_p2_conformance.c)

    # LIST reconstruction in the batch reader (Arrow offsets + child buffers)
    add_carquet_test(test_batch_nested tests/test_batch_nested.c)

    # Nested write helper: auto-shredding LIST/MAP from offsets + validity
    add_carquet_test(test_nested_write tests/test_nested_write.c)

    # Wide (SIMD) bit-unpack fast path vs independent oracle
    add_carquet_test(test_bitunpack_wide tests/test_bitunpack_wide.c)

    # Page-level filter (predicate pushdown via column index + offset index)
    add_carquet_test(test_page_filter tests/test_page_filter.c)

    # Arbitrary [offset, limit) row-range read (predicate-independent slice)
    add_carquet_test(test_row_range tests/test_row_range.c)

    # Module-level dictionary decode for BYTE_ARRAY / FIXED_LEN_BYTE_ARRAY
    add_carquet_test(test_dict_decode tests/test_dict_decode.c)

    # Signed big-endian ordering of DECIMAL(FLBA/BYTE_ARRAY) min/max statistics
    add_carquet_test(test_decimal_stats tests/test_decimal_stats.c)

    # Per-reason Snappy rejection diagnostics + first-data-page offset resolution
    add_carquet_test(test_snappy_diag tests/test_snappy_diag.c)

    # Page-granular parallel encode + background I/O: byte-identical output
    add_carquet_test(test_parallel_write tests/test_parallel_write.c)

    # Writing dictionary-encoded (categorical) Arrow arrays resolves to values
    add_carquet_test(test_arrow_dictionary tests/test_arrow_dictionary.c)

    # CLI value/type formatting (DATE, TIME, TIMESTAMP rendering)
    add_carquet_test(test_cli_format "tests/test_cli_format.c;src/cli/commands.c")

    # Append row groups to an existing file (metadata carryover, bloom/page-index
    # preservation, nullable + multi-column append, structural-mismatch rejection)
    add_carquet_test(test_append tests/test_append.c)

    # Unit tests for previously-uncovered internal modules
    add_carquet_test(test_float16 tests/test_float16.c)          # binary16 -> binary32
    add_carquet_test(test_geo_wkb tests/test_geo_wkb.c)          # WKB geospatial stats walker
    add_carquet_test(test_worker_pool tests/test_worker_pool.c)  # persistent reader thread pool
    # This test uses C11 <stdatomic.h> directly; MSVC only enables it under an opt-in flag.
    if(MSVC)
        target_compile_options(test_worker_pool PRIVATE /experimental:c11atomics)
    endif()
    add_carquet_test(test_custom_codec tests/test_custom_codec.c) # pluggable compression codec API

    # Full-pipeline "real dataset" regression net: multi-type nullable columns,
    # several row groups, compression, statistics, predicate pushdown, batch read
    add_carquet_test(test_real_world tests/test_real_world.c)
endif()

# Examples
if(CARQUET_BUILD_EXAMPLES)
    # Helper function to configure example targets
    function(add_carquet_example EXAMPLE_NAME EXAMPLE_SOURCE)
        add_executable(${EXAMPLE_NAME} ${EXAMPLE_SOURCE})
        target_link_libraries(${EXAMPLE_NAME} PRIVATE carquet)
        if(CARQUET_COMPILER_GCC_LIKE)
            target_compile_options(${EXAMPLE_NAME} PRIVATE -Wno-unused-result)
        endif()
    endfunction()

    add_carquet_example(example_basic_write_read examples/basic_write_read.c)
    add_carquet_example(example_data_types examples/data_types.c)
    add_carquet_example(example_compression examples/compression_codecs.c)
    add_carquet_example(example_nullable examples/nullable_columns.c)
    add_carquet_example(example_advanced examples/advanced_features.c)
    add_carquet_example(example_append examples/append_rows.c)
    add_carquet_example(example_page_filter examples/page_filter.c)
    add_carquet_example(example_nested examples/nested_data.c)
    add_carquet_example(example_arrow_c_data examples/arrow_c_data.c)
endif()

# Benchmarks
if(CARQUET_BUILD_BENCHMARKS)
    add_executable(benchmark_carquet benchmark/benchmark_carquet.c)
    target_link_libraries(benchmark_carquet PRIVATE carquet)
    if(UNIX AND NOT APPLE)
        target_link_libraries(benchmark_carquet PRIVATE m)  # libm for math functions
    endif()
    if(CARQUET_COMPILER_GCC_LIKE)
        target_compile_options(benchmark_carquet PRIVATE -Wno-unused-result)
    endif()

    # Test file generator (C replacement for generate_test_files.py - no extra deps)
    add_executable(generate_test_files benchmark/generate_test_files.c)
    target_link_libraries(generate_test_files PRIVATE carquet)
    if(CARQUET_COMPILER_GCC_LIKE)
        target_compile_options(generate_test_files PRIVATE -Wno-unused-result)
    endif()

    if(CARQUET_BUILD_ARROW_CPP_BENCHMARK)
        include(CheckLanguage)
        check_language(CXX)
        if(CMAKE_CXX_COMPILER)
            enable_language(CXX)

            if(CARQUET_ARROW_CPP_ROOT)
                find_package(Arrow CONFIG QUIET PATHS ${CARQUET_ARROW_CPP_ROOT})
                find_package(Parquet CONFIG QUIET PATHS ${CARQUET_ARROW_CPP_ROOT})
            endif()

            if(NOT Arrow_FOUND)
                find_package(Arrow CONFIG QUIET)
            endif()
            if(NOT Parquet_FOUND)
                find_package(Parquet CONFIG QUIET)
            endif()

            set(CARQUET_ARROW_TARGET "")
            set(CARQUET_PARQUET_TARGET "")
            if(TARGET Arrow::arrow_shared)
                set(CARQUET_ARROW_TARGET Arrow::arrow_shared)
            elseif(TARGET Arrow::arrow_static)
                set(CARQUET_ARROW_TARGET Arrow::arrow_static)
            endif()
            if(TARGET Parquet::parquet_shared)
                set(CARQUET_PARQUET_TARGET Parquet::parquet_shared)
            elseif(TARGET Parquet::parquet_static)
                set(CARQUET_PARQUET_TARGET Parquet::parquet_static)
            endif()

            if(CARQUET_ARROW_TARGET AND CARQUET_PARQUET_TARGET)
                add_executable(benchmark_arrow_cpp benchmark/benchmark_arrow_cpp.cpp)
                target_link_libraries(benchmark_arrow_cpp PRIVATE
                    ${CARQUET_PARQUET_TARGET}
                    ${CARQUET_ARROW_TARGET}
                )
                set_target_properties(benchmark_arrow_cpp PROPERTIES
                    CXX_STANDARD 20
                    CXX_STANDARD_REQUIRED YES
                )
                if(CARQUET_COMPILER_GCC_LIKE)
                    target_compile_options(benchmark_arrow_cpp PRIVATE -Wno-unused-result)
                endif()
            else()
                message(STATUS
                    "Skipping benchmark_arrow_cpp: Arrow/Parquet C++ libraries were not found. "
                    "Install system packages or configure -DCARQUET_ARROW_CPP_ROOT=/path/to/prefix. "
                    "The benchmark remains fully opt-in to keep normal builds lean.")
            endif()
        else()
            message(STATUS
                "Skipping benchmark_arrow_cpp: no C++ compiler available for the optional Arrow C++ benchmark.")
        endif()
    endif()

    # Profiling benchmark for CPU profiling (sample, Instruments, perf)
    # Only build on Unix systems (uses dirent.h and Unix profiling tools)
    if(NOT WIN32)
        add_executable(profile_core benchmark/profile_core.c)
        target_link_libraries(profile_core PRIVATE carquet)
        if(NOT APPLE)
            target_link_libraries(profile_core PRIVATE m)
        endif()

        # Profiling binaries (used by profiling/run_profiler.py)
        add_executable(profile_read profiling/profile_read.c)
        target_link_libraries(profile_read PRIVATE carquet)
        if(NOT APPLE)
            target_link_libraries(profile_read PRIVATE m)
        endif()

        add_executable(profile_write profiling/profile_write.c)
        target_link_libraries(profile_write PRIVATE carquet)
        if(NOT APPLE)
            target_link_libraries(profile_write PRIVATE m)
        endif()

        add_executable(profile_micro profiling/profile_micro.c)
        target_link_libraries(profile_micro PRIVATE carquet)
        target_include_directories(profile_micro PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
        if(NOT APPLE)
            target_link_libraries(profile_micro PRIVATE m)
        endif()
    endif()
endif()

# Interoperability tests
if(CARQUET_BUILD_INTEROP)
    if(NOT WIN32)
        add_executable(test_interop interop/test_interop.c)
        target_link_libraries(test_interop PRIVATE carquet)
        target_include_directories(test_interop PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
        if(CARQUET_COMPILER_GCC_LIKE)
            target_compile_options(test_interop PRIVATE -Wno-unused-result)
        endif()

        add_executable(roundtrip_writer interop/roundtrip_writer.c)
        target_link_libraries(roundtrip_writer PRIVATE carquet)
        if(CARQUET_COMPILER_GCC_LIKE)
            target_compile_options(roundtrip_writer PRIVATE -Wno-unused-result)
        endif()
    endif()
endif()

# Fuzzing
option(CARQUET_BUILD_FUZZ "Build fuzz targets" OFF)
if(CARQUET_BUILD_FUZZ)
    add_subdirectory(fuzz)
endif()

# Documentation (Doxygen) — optional; the `docs` target only exists when
# Doxygen is actually installed, so a plain build never requires it.
option(CARQUET_BUILD_DOCS "Build API documentation with Doxygen (needs doxygen; graphviz for diagrams)" OFF)
if(CARQUET_BUILD_DOCS AND CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    find_package(Doxygen OPTIONAL_COMPONENTS dot)
    if(DOXYGEN_FOUND)
        set(_carquet_docs_out ${CMAKE_CURRENT_BINARY_DIR}/docs)
        # Reuse the checked-in Doxyfile verbatim, only redirecting its output
        # into the build tree. @INCLUDE pulls in the real config; the trailing
        # override wins because Doxygen applies later settings last.
        set(_carquet_doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile.docs)
        file(WRITE ${_carquet_doxyfile}
            "@INCLUDE_PATH = ${CMAKE_CURRENT_SOURCE_DIR}\n"
            "@INCLUDE      = Doxyfile\n"
            "OUTPUT_DIRECTORY = ${_carquet_docs_out}\n")
        add_custom_target(docs
            COMMAND ${DOXYGEN_EXECUTABLE} ${_carquet_doxyfile}
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}  # so INPUT/IMAGE_PATH resolve
            COMMENT "Generating API docs -> ${_carquet_docs_out}/html/index.html"
            VERBATIM)
        if(NOT TARGET Doxygen::dot)
            message(STATUS "Doxygen found but graphviz 'dot' is missing: docs will build without diagrams")
        endif()
    else()
        message(WARNING "CARQUET_BUILD_DOCS=ON but Doxygen was not found; the 'docs' target is unavailable")
    endif()
endif()

# CLI tool
if(CARQUET_BUILD_CLI)
    add_executable(carquet_cli
        src/cli/main.c
        src/cli/commands.c
        src/cli/codegen.c
        src/cli/codegen_read.c
        src/cli/codegen_write.c
    )
    set_target_properties(carquet_cli PROPERTIES OUTPUT_NAME carquet)
    target_link_libraries(carquet_cli PRIVATE carquet)
    target_include_directories(carquet_cli PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
    if(CARQUET_COMPILER_GCC_LIKE)
        target_compile_options(carquet_cli PRIVATE -Wno-unused-result)
    endif()
endif()

# Installation (skip when used as a subdirectory via FetchContent/add_subdirectory)
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    include(GNUInstallDirs)
    install(TARGETS carquet
        EXPORT carquetTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

    # When carquet is STATIC and we fetched compression deps via FetchContent,
    # CMake propagates those deps to the install export's link interface
    # (consumers must link them transitively). Add them to the same export set
    # so the install(EXPORT) generation succeeds — this is required on Windows
    # static builds where no system zlib/lz4/zstd is available.
    if(NOT CARQUET_BUILD_SHARED)
        set(_carquet_fetched_deps)
        if(NOT ZSTD_FOUND AND TARGET libzstd_static)
            list(APPEND _carquet_fetched_deps libzstd_static)
        endif()
        if(NOT ZLIB_FOUND AND TARGET zlibstatic)
            list(APPEND _carquet_fetched_deps zlibstatic)
        endif()
        if(NOT LZ4_FOUND AND TARGET lz4_static)
            list(APPEND _carquet_fetched_deps lz4_static)
        endif()
        if(BROTLI_LIBRARIES AND NOT BROTLI_FOUND AND TARGET brotlienc)
            list(APPEND _carquet_fetched_deps brotlienc brotlidec brotlicommon)
        endif()
        if(_carquet_fetched_deps)
            # Their INTERFACE_INCLUDE_DIRECTORIES point at build-tree paths
            # which CMake refuses to install. Consumers of carquet only need
            # to link these archives, not include their headers, so wrap any
            # interface include dirs in $<BUILD_INTERFACE:> to suppress them
            # from the installed export.
            foreach(_dep IN LISTS _carquet_fetched_deps)
                get_target_property(_inc ${_dep} INTERFACE_INCLUDE_DIRECTORIES)
                if(_inc)
                    set(_wrapped)
                    foreach(_p IN LISTS _inc)
                        list(APPEND _wrapped "$<BUILD_INTERFACE:${_p}>")
                    endforeach()
                    set_target_properties(${_dep} PROPERTIES
                        INTERFACE_INCLUDE_DIRECTORIES "${_wrapped}")
                endif()
            endforeach()
            install(TARGETS ${_carquet_fetched_deps}
                EXPORT carquetTargets
                ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            )
        endif()
    endif()

    if(CARQUET_BUILD_CLI)
        install(TARGETS carquet_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    endif()

    install(DIRECTORY include/carquet
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    install(EXPORT carquetTargets
        FILE carquetTargets.cmake
        NAMESPACE carquet::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/carquet
    )

    # Package configuration
    include(CMakePackageConfigHelpers)
    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/carquetConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )

    # Record which transitive imported deps the installed target references so
    # the package config can find_dependency() them for consumers.
    set(CARQUET_NEEDS_THREADS 0)
    if(UNIX)
        set(CARQUET_NEEDS_THREADS 1)
    endif()
    set(CARQUET_NEEDS_OPENMP 0)
    if(OpenMP_C_FOUND)
        set(CARQUET_NEEDS_OPENMP 1)
    endif()

    configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/carquetConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/carquetConfig.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/carquet
    )

    # Install the generated package config so find_package(carquet) works from
    # an install tree (without this, only carquetTargets.cmake ships and
    # find_package can't locate the package).
    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/carquetConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/carquetConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/carquet
    )

    # pkg-config file so non-CMake build systems (Meson, Autotools, plain make)
    # can consume carquet via `pkg-config --cflags --libs carquet`.
    # Libs.private lists the transitive static-link deps; it only matters when a
    # consumer links carquet statically (`pkg-config --static`). The standard
    # system link names are assumed — correct for system/distro installs, but
    # not for CARQUET_BUNDLE_DEPS/FetchContent artifacts whose archives carry
    # mangled names (those are meant for CMake export consumption, not .pc).
    set(CARQUET_PC_PRIVATE_LIBS "-lzstd -lz -llz4")
    if(BROTLI_LIBRARIES)
        string(APPEND CARQUET_PC_PRIVATE_LIBS " -lbrotlienc -lbrotlidec -lbrotlicommon")
    endif()
    if(UNIX)
        string(APPEND CARQUET_PC_PRIVATE_LIBS " -lm")
    endif()
    if(CARQUET_NEEDS_THREADS)
        string(APPEND CARQUET_PC_PRIVATE_LIBS " -lpthread")
    endif()
    if(CARQUET_NEEDS_OPENMP AND OpenMP_C_FLAGS)
        string(APPEND CARQUET_PC_PRIVATE_LIBS " ${OpenMP_C_FLAGS}")
    endif()

    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/carquet.pc.in"
        "${CMAKE_CURRENT_BINARY_DIR}/carquet.pc"
        @ONLY
    )
    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/carquet.pc"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
    )
endif()
