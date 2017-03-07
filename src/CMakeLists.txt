# openmp integration
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOPENMP=1 ${OpenMP_CXX_FLAGS}")
endif()

# needed for concat.h
include(CheckCXXSourceRuns)
check_cxx_source_runs("
        #include <stdlib.h>
        #include <fcntl.h>
        #include <stdio.h>

        int main()
        {
          FILE* tmpf = tmpfile();
          int input_desc = fileno(tmpf);
          int test =  posix_fadvise (input_desc, 0, 0, POSIX_FADV_SEQUENTIAL);
          fclose(tmpf);
        }"
        HAVE_POSIX_FADVISE)
if(HAVE_POSIX_FADVISE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_POSIX_FADVISE=1 ")
endif()


#SSE
if(${HAVE_AVX2})
    add_definitions("-DAVX2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -Wa,-q")
elseif(${HAVE_SSE4_1})
    add_definitions("-DSSE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.1")
else()
    include(CheckSSEFeatures)
    add_definitions(${CXX_DFLAGS})
    if(${HAVE_AVX2_EXTENSIONS})
        add_definitions("-DAVX2")
    # debugging
    #   add_definitions("-DSSE")
    else()
    if(${HAVE_SSE4_1_EXTENSIONS})
        add_definitions("-DSSE")
    else()
        message( FATAL_ERROR "At least SSE4.2 is needed to compile!" )
    endif(${HAVE_SSE4_1_EXTENSIONS})
    endif(${HAVE_AVX2_EXTENSIONS})
endif()

# set flags
set(MMSEQS_WARNINGS "-pedantic -Wall -Wextra -Winline -Wdisabled-optimization -Wno-unused-parameter")
set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -std=c++0x -m64 ${MMSEQS_WARNINGS}")
set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL} -Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDBINFO} -O2 -g")

if(CMAKE_COMPILER_IS_ICC)
    set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -ipo -no-prec-div -xHost -fno-exceptions")
else()
    set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -ffast-math -fno-exceptions -ftree-vectorize -fno-strict-aliasing")
endif()

# Compiler-specific features
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif()

# MPI
if(${HAVE_MPI})
    find_package(MPI QUIET)
    if(MPI_FOUND)
        set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
        set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
        include_directories(${MPI_INCLUDE_PATH})
        ADD_DEFINITIONS("-DHAVE_MPI")
        message("-- Found MPI")

    endif()
endif()

find_package(ZLIB QUIET)
if(ZLIB_FOUND)
    include_directories(${ZLIB_INCLUDE_DIRS})
    ADD_DEFINITIONS("-DHAVE_ZLIB")
endif()
include_directories(generatedfiles)
add_subdirectory(alignment)
add_subdirectory(clustering)
add_subdirectory(commons)
add_subdirectory(prefiltering)
add_subdirectory(util)
add_subdirectory(workflow)
add_subdirectory(test)

set(mmseqs_source_files
    CommandDeclarations.h
    mmseqs.cpp
)


# generate class from data files
FILE(GLOB_RECURSE COMPILED_RESOURCES "${PROJECT_SOURCE_DIR}/data/*")
FOREACH(INPUT_FILE ${COMPILED_RESOURCES})
    get_filename_component(INPUT_FILE_NAME ${INPUT_FILE} NAME)
    SET(OUTPUT_FILE ${PROJECT_SOURCE_DIR}/src/generatedfiles/${INPUT_FILE_NAME}.h)
    execute_process(COMMAND xxd -i ${INPUT_FILE_NAME} ${OUTPUT_FILE}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data/)
    message("-- Compile ${INPUT_FILE_NAME} to ${OUTPUT_FILE}")
    file(READ ${OUTPUT_FILE} CONTENTS)
    #file(APPEND ${PROJECT_SOURCE_DIR}/src/commons/externalfile.h "${CONTENTS}")
    #file(REMOVE ${OUTPUT_FILE})
ENDFOREACH()

add_executable(mmseqs ${mmseqs_source_files} )

target_link_libraries(mmseqs util workflow)
if(MPI_FOUND)
    target_link_libraries(mmseqs ${MPI_LIBRARIES})
    if(MPI_COMPILE_FLAGS)
        set_target_properties(mmseqs PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    endif()

    if(MPI_LINK_FLAGS)
        set_target_properties(mmseqs PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    endif()
endif()

if(ZLIB_FOUND)
    target_link_libraries(mmseqs ${ZLIB_LIBRARIES})
endif()

install(TARGETS mmseqs DESTINATION bin)
