set(HAVE_MPI 0 CACHE BOOL "Have MPI")
set(HAVE_TESTS 0 CACHE BOOL "Have Tests")
set(HAVE_SHELLCHECK 1 CACHE BOOL "Have ShellCheck")
set(HAVE_GPROF 0 CACHE BOOL "Have GPROF Profiler")
set(ENABLE_WERROR 0 CACHE BOOL "Enable Warnings as Errors")
#set(DISABLE_LTO 0 CACHE BOOL "Disable link-time optimization in non-debug builds")
set(REQUIRE_OPENMP 1 CACHE BOOL "Require availability of OpenMP")

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
        ${workflow_source_files}
        CommandDeclarations.h
        MMseqsBase.cpp
        )

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

#if (NOT DISABLE_LTO)
#    string(TOUPPER "${CMAKE_BUILD_TYPE}" uppercase_CMAKE_BUILD_TYPE)
#    if (uppercase_CMAKE_BUILD_TYPE MATCHES "^(RELEASE|RELWITHDEBINFO|MINSIZEREL)$")
#        cmake_policy(SET CMP0069 NEW)
#        include(CheckIPOSupported)
#        check_ipo_supported(RESULT IPO_SUPPORTED)
#        if (IPO_SUPPORTED)
#            set_property(TARGET mmseqs-framework PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
#        endif ()
#    endif ()
#endif ()

append_target_property(mmseqs-framework COMPILE_FLAGS ${MMSEQS_CXX_FLAGS} -pedantic -Wall -Wextra -Wdisabled-optimization)
append_target_property(mmseqs-framework LINK_FLAGS ${MMSEQS_CXX_FLAGS} -pedantic -Wall -Wextra -Wdisabled-optimization)
if (NOT EMSCRIPTEN)
    append_target_property(mmseqs-framework COMPILE_FLAGS -fno-exceptions)
    append_target_property(mmseqs-framework LINK_FLAGS -fno-exceptions)
endif()

if (ENABLE_WERROR)
    append_target_property(mmseqs-framework COMPILE_FLAGS -Werror)
    append_target_property(mmseqs-framework LINK_FLAGS -Werror)
endif()

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

if (NOT DISABLE_IPS4O)
    find_package(Atomic)
    if (ATOMIC_FOUND)
        target_link_libraries(mmseqs-framework ${ATOMIC_LIBRARIES})
        target_compile_definitions(mmseqs-framework PUBLIC -DENABLE_IPS4O=1)
        message("-- IPS4O sorting works")
    else ()
        message("-- OMPTL sorting fallback")
    endif ()
else ()
    message("-- OMPTL sorting fallback")
endif ()

target_link_libraries(mmseqs-framework tinyexpr libzstd_static microtar)
if (CYGWIN)
    target_link_libraries(mmseqs-framework nedmalloc)
endif ()

if (EMSCRIPTEN)
    target_compile_definitions(mmseqs-framework PUBLIC -DHAVE_ZLIB=1 -DHAVE_BZLIB=1)
    append_target_property(mmseqs-framework COMPILE_FLAGS -s USE_ZLIB=1 -s USE_BZIP2=1)
    append_target_property(mmseqs-framework LINK_FLAGS -s USE_ZLIB=1 -s USE_BZIP2=1)
else ()
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
endif ()

# MPI
if (HAVE_MPI)
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


find_package(OpenMP QUIET)
if (OPENMP_FOUND)
    message("-- Found OpenMP")
    target_compile_definitions(mmseqs-framework PUBLIC -DOPENMP=1)
    # For GCC we dont want to do this since it breaks macOS static builds
    # It will link libgomp.a internally (through -fopenmp I guess)
    # and also link libgomp.dylib thus breaking static builds
    if (NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        target_link_libraries(mmseqs-framework ${OpenMP_CXX_LIBRARIES})
    endif()
    append_target_property(mmseqs-framework COMPILE_FLAGS ${OpenMP_CXX_FLAGS})
    append_target_property(mmseqs-framework LINK_FLAGS ${OpenMP_CXX_FLAGS})
elseif (REQUIRE_OPENMP)
    message(FATAL_ERROR "-- Could not find OpenMP. Skip check with -DREQUIRE_OPENMP=0.")
endif ()

if (HAVE_GPROF)
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
    set(mmseqs_source_files mmseqs.cpp)

    add_executable(mmseqs${EXE_SUFFIX} ${mmseqs_source_files})
    mmseqs_setup_derived_target(mmseqs${EXE_SUFFIX})
    target_link_libraries(mmseqs${EXE_SUFFIX} version)
    install(TARGETS mmseqs${EXE_SUFFIX} DESTINATION bin)

    if (HAVE_TESTS)
        add_subdirectory(test)
    endif ()
endif ()
