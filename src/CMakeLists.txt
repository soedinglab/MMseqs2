set(HAVE_MPI 0 CACHE BOOL "Have MPI")
set(HAVE_AVX2 0 CACHE BOOL "Have AVX2")
set(HAVE_SSE4_1 0 CACHE BOOL "Have SSE4.1")
set(HAVE_NEON 0 CACHE BOOL "Have NEON")
set(HAVE_TESTS 1 CACHE BOOL "Have Tests")
set(HAVE_SHELLCHECK 1 CACHE BOOL "Have ShellCheck")
set(HAVE_GPROF 0 CACHE BOOL "Have GPROF Profiler")
set(ENABLE_WERROR 0 CACHE BOOL "Enable Warnings as Errors")

include(AppendTargetProperty)

add_subdirectory(alignment)
add_subdirectory(clustering)
add_subdirectory(commons)
add_subdirectory(linclust)
add_subdirectory(multihit)
add_subdirectory(prefiltering)
add_subdirectory(taxonomy)
add_subdirectory(util)
add_subdirectory(workflow)

add_library(mmseqs-framework
        $<TARGET_OBJECTS:alp>
        $<TARGET_OBJECTS:ksw2>
        $<TARGET_OBJECTS:cacode>
        ${alignment_header_files}
        ${alignment_source_files}
        ${clustering_header_files}
        ${clustering_source_files}
        ${commons_header_files}
        ${commons_source_files}
        ${prefiltering_header_files}
        ${prefiltering_source_files}
        ${multihit_header_files}
        ${multihit_source_files}
        ${taxonomy_header_files}
        ${taxonomy_source_files}
        ${linclust_source_files}
        ${util_header_files}
        ${util_source_files}
        ${workflow_source_files})

target_include_directories(mmseqs-framework PUBLIC ${CMAKE_BINARY_DIR}/generated)
target_include_directories(mmseqs-framework PUBLIC ${PROJECT_BINARY_DIR}/generated)
target_include_directories(mmseqs-framework PUBLIC alignment)
target_include_directories(mmseqs-framework PUBLIC clustering)
target_include_directories(mmseqs-framework PUBLIC commons)
target_include_directories(mmseqs-framework PUBLIC multihit)
target_include_directories(mmseqs-framework PUBLIC prefiltering)
target_include_directories(mmseqs-framework PUBLIC linclust)
target_include_directories(mmseqs-framework PUBLIC taxonomy)
target_include_directories(mmseqs-framework PUBLIC util)
target_include_directories(mmseqs-framework PUBLIC .)

add_dependencies(mmseqs-framework generated)

if (CMAKE_BUILD_TYPE MATCHES RELEASE OR CMAKE_BUILD_TYPE MATCHES RELWITHDEBINFO)
    if (CMAKE_COMPILER_IS_ICC)
        append_target_property(mmseqs-framework COMPILE_FLAGS -ipo -no-prec-div -xHost)
        append_target_property(mmseqs-framework LINK_FLAGS -ipo -no-prec-div -xHost)
    else ()
        append_target_property(mmseqs-framework COMPILE_FLAGS -ffast-math -ftree-vectorize -fno-strict-aliasing)
        append_target_property(mmseqs-framework LINK_FLAGS -ffast-math -ftree-vectorize -fno-strict-aliasing)
    endif ()
endif ()

append_target_property(mmseqs-framework COMPILE_FLAGS ${MMSEQS_CXX_FLAGS} -fno-exceptions -pedantic -Wall -Wextra -Wdisabled-optimization)
append_target_property(mmseqs-framework LINK_FLAGS ${MMSEQS_CXX_FLAGS} -fno-exceptions -pedantic -Wall -Wextra -Wdisabled-optimization)

if (${ENABLE_WERROR})
    append_target_property(mmseqs-framework COMPILE_FLAGS -Werror)
    append_target_property(mmseqs-framework LINK_FLAGS -Werror)
endif()

if (CYGWIN OR ${HAVE_NEON})
    target_compile_definitions(mmseqs-framework PUBLIC -D_GNU_SOURCE=1)
endif ()

# needed for concat.h
include(CheckCXXSourceRuns)
check_cxx_source_runs("
        #include <stdlib.h>
        #include <fcntl.h>
        #include <stdio.h>

        int main() {
          FILE* tmpf = tmpfile();
          int input_desc = fileno(tmpf);
          int test = posix_fadvise(input_desc, 0, 0, POSIX_FADV_SEQUENTIAL);
          fclose(tmpf);
          return 0;
        }"
        HAVE_POSIX_FADVISE)
if (HAVE_POSIX_FADVISE)
    target_compile_definitions(mmseqs-framework PUBLIC -DHAVE_POSIX_FADVISE=1)
endif ()

check_cxx_source_runs("
        #include <stdlib.h>
        #include <fcntl.h>
        #include <stdio.h>
        #include <sys/mman.h>

        int main() {
          FILE* tmpf = tmpfile();
          void *tmp = mmap(NULL, 32, PROT_READ, MAP_SHARED, fileno(tmpf), 0);
          int test = posix_madvise(tmp, 32, POSIX_MADV_SEQUENTIAL);
          fclose(tmpf);
          return 0;
        }"
        HAVE_POSIX_MADVISE)
if (HAVE_POSIX_MADVISE)
    target_compile_definitions(mmseqs-framework PUBLIC -DHAVE_POSIX_MADVISE=1)
endif ()

#SSE
if (${HAVE_AVX2})
    target_compile_definitions(mmseqs-framework PUBLIC -DAVX2=1)
    append_target_property(mmseqs-framework COMPILE_FLAGS -mavx2 -Wa,-q)
    append_target_property(mmseqs-framework LINK_FLAGS -mavx2 -Wa,-q)
elseif (${HAVE_SSE4_1})
    target_compile_definitions(mmseqs-framework PUBLIC -DSSE=1)
    append_target_property(mmseqs-framework COMPILE_FLAGS -msse4.1)
    append_target_property(mmseqs-framework LINK_FLAGS -msse4.1)
elseif (${HAVE_NEON})
    target_compile_definitions(mmseqs-framework PUBLIC -DSSE=1 -DNEON=1)
else ()
    include(CheckSSEFeatures)
    append_target_property(mmseqs-framework COMPILE_FLAGS ${SSE_FLAGS})
    append_target_property(mmseqs-framework LINK_FLAGS ${SSE_FLAGS})
    if (${HAVE_AVX2_EXTENSIONS})
        target_compile_definitions(mmseqs-framework PUBLIC -DAVX2=1)
        # debugging
        #   list(APPEND MMSEQS_DEFINITIONS -DSSE=1)
    else ()
        if (${HAVE_SSE4_1_EXTENSIONS})
            target_compile_definitions(mmseqs-framework PUBLIC -DSSE=1)
        else ()
            message(FATAL_ERROR "At least SSE4.2 is needed to compile!")
        endif (${HAVE_SSE4_1_EXTENSIONS})
    endif (${HAVE_AVX2_EXTENSIONS})
endif ()

# tinyexpr
target_link_libraries(mmseqs-framework tinyexpr)

find_package(ZLIB QUIET)
if (ZLIB_FOUND)
    message("-- Found ZLIB")
    set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
    set(OLD_CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
    set(OLD_CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
    set(CMAKE_REQUIRED_FLAGS -lz)
    set(CMAKE_REQUIRED_INCLUDES ${ZLIB_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${ZLIB_LIBRARIES})
    check_cxx_source_runs("
        #include <zlib.h>
        int main() { gzFile file; return 0; }"
        HAVE_ZLIB_CHECK)
    set(CMAKE_REQUIRED_FLAGS ${OLD_CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_INCLUDES ${OLD_CMAKE_REQUIRED_INCLUDES})
    set(CMAKE_REQUIRED_LIBRARIES ${OLD_CMAKE_REQUIRED_LIBRARIES})
    if(HAVE_ZLIB_CHECK)
        message("-- ZLIB works")
        target_include_directories(mmseqs-framework PUBLIC ${ZLIB_INCLUDE_DIRS})
        target_compile_definitions(mmseqs-framework PUBLIC -DHAVE_ZLIB=1)
        target_link_libraries(mmseqs-framework ${ZLIB_LIBRARIES})
    else ()
        message("-- ZLIB does not work")
    endif()
else ()
    message("-- Could not find ZLIB")
endif ()

find_package(BZip2 QUIET)
if (BZIP2_FOUND)
    message("-- Found BZLIB")
    set(OLD_CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS})
    set(OLD_CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES})
    set(OLD_CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES})
    set(CMAKE_REQUIRED_FLAGS -lbz2)
    set(CMAKE_REQUIRED_INCLUDES ${BZIP_INCLUDE_DIRS})
    set(CMAKE_REQUIRED_LIBRARIES ${BZIP2_LIBRARIES})
    check_cxx_source_runs("
        #include <bzlib.h>
        int main() { bz_stream stream; return 0; }"
        HAVE_BZLIB_CHECK)
    set(CMAKE_REQUIRED_FLAGS ${OLD_CMAKE_REQUIRED_FLAGS})
    set(CMAKE_REQUIRED_INCLUDES ${OLD_CMAKE_REQUIRED_INCLUDES})
    set(CMAKE_REQUIRED_LIBRARIES ${OLD_CMAKE_REQUIRED_LIBRARIES})
    if(HAVE_BZLIB_CHECK)
        message("-- BZLIB works")
        target_include_directories(mmseqs-framework PUBLIC ${BZIP_INCLUDE_DIRS})
        target_compile_definitions(mmseqs-framework PUBLIC -DHAVE_BZLIB=1)
        target_link_libraries(mmseqs-framework ${BZIP2_LIBRARIES})
    else ()
        message("-- BZLIB does not work")
    endif()
else ()
    message("-- Could not find BZLIB")
endif ()

# MPI
if (${HAVE_MPI})
    find_package(MPI REQUIRED)
    if (MPI_FOUND)
        message("-- Found MPI")
        target_include_directories(mmseqs-framework PUBLIC ${MPI_INCLUDE_PATH})
        #Hack (OMPI_SKIP_MPICXX=1): https://github.com/open-mpi/ompi/issues/5157#issuecomment-388495496
        target_compile_definitions(mmseqs-framework PUBLIC -DHAVE_MPI=1 -DOMPI_SKIP_MPICXX=1)
        target_link_libraries(mmseqs-framework ${MPI_LIBRARIES})
        append_target_property(mmseqs-framework COMPILE_FLAGS ${MPI_COMPILE_FLAGS})
        append_target_property(mmseqs-framework LINK_FLAGS ${MPI_LINK_FLAGS})
    endif ()
endif ()

# openmp integration
find_package(OpenMP)
if (OPENMP_FOUND)
    message("-- Found OpenMP")
    target_compile_definitions(mmseqs-framework PUBLIC -DOPENMP=1)
    target_link_libraries(mmseqs-framework ${OpenMP_CXX_LIBRARIES})
    append_target_property(mmseqs-framework COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
    append_target_property(mmseqs-framework LINK_FLAGS ${OpenMP_CXX_FLAGS})
endif ()

if (${HAVE_GPROF})
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-pg GPROF_FOUND)
    if (GPROF_FOUND)
        append_target_property(mmseqs-framework COMPILE_FLAGS -pg)
        append_target_property(mmseqs-framework LINK_FLAGS -pg)
    else ()
        message(FATAL_ERROR "-- Could not find GPROF")
    endif ()
endif ()

if (NOT FRAMEWORK_ONLY)
    include(MMseqsSetupDerivedTarget)
    add_subdirectory(version)
    set(mmseqs_source_files
            CommandDeclarations.h
            mmseqs.cpp
            )

    add_executable(mmseqs ${mmseqs_source_files})
    mmseqs_setup_derived_target(mmseqs)
    target_link_libraries(mmseqs version)
    install(TARGETS mmseqs DESTINATION bin)

    if (HAVE_TESTS)
        add_subdirectory(test)
    endif ()
endif ()
