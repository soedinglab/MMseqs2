cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(mmseqs)


# openmp integration
find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DOPENMP=1 ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOPENMP=1 ${OpenMP_CXX_FLAGS}")
endif()

#SSE
INCLUDE(CheckSSEFeatures)
ADD_DEFINITIONS(${CXX_DFLAGS})
if(${HAVE_AVX2_EXTENSIONS})
    ADD_DEFINITIONS("-DAVX2")
#    ADD_DEFINITIONS("-DSSE") #debugging
else()
if(${HAVE_SSE4_2_EXTENSIONS})
    ADD_DEFINITIONS("-DSSE")
else()
    message( FATAL_ERROR "SSE4.2 is needed to run compile! CMake will exit." )
endif(${HAVE_SSE4_2_EXTENSIONS})
endif(${HAVE_AVX2_EXTENSIONS})
# set flags

set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -Wall -Wno-write-strings -std=c++0x -m64")
set(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g ")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "${CMAKE_CXX_FLAGS_MINSIZEREL} -Os -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -ffast-math -ftree-vectorize -fno-strict-aliasing -Wno-write-strings")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDBINFO} -O2 -g")
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")
set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")



# Compiler-specific C++11 activation.
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
	execute_process(
		COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION
	)
	if (NOT (GCC_VERSION VERSION_GREATER 4.4 OR GCC_VERSION VERSION_EQUAL 4.4))
		message(FATAL_ERROR "${PROJECT_NAME} requires g++ 4.4 or greater.")
	endif ()
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
else ()
		message(FATAL_ERROR "Your C++ compiler does not support C++0x.")
endif ()

# MPI
if(${MPI_BUILD})
    find_package(MPI QUIET)
    if(MPI_FOUND)
        set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
        set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
        include_directories(${MPI_INCLUDE_PATH})
        ADD_DEFINITIONS("-DHAVE_MPI")
    endif()
endif()
add_subdirectory(alignment)
add_subdirectory(clustering)
add_subdirectory(commons)
add_subdirectory(prefiltering)
add_subdirectory(util)
add_subdirectory(workflow)
add_subdirectory(orf)
add_subdirectory(test)
set(mmseqs_source_files mmseqs.cpp )
add_executable(mmseqs ${mmseqs_source_files})

target_link_libraries(mmseqs util)
target_link_libraries(mmseqs workflow)
target_link_libraries(mmseqs orf)
target_link_libraries(mmseqs ffindex)
if(MPI_FOUND)
target_link_libraries(mmseqs ${MPI_LIBRARIES})
endif()
install (TARGETS mmseqs
		 DESTINATION bin)

# copy shell scripts to bin
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/workflow/searching.sh
     DESTINATION ${CMAKE_SOURCE_DIR}/bin/
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
     GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/workflow/clustering.sh
     DESTINATION ${CMAKE_SOURCE_DIR}/bin/
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
     GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/workflow/cascaded_clustering.sh
     DESTINATION ${CMAKE_SOURCE_DIR}/bin/
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
     GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/workflow/blastp.sh
     DESTINATION ${CMAKE_SOURCE_DIR}/bin/
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
     GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/workflow/blastpgp.sh
     DESTINATION ${CMAKE_SOURCE_DIR}/bin/
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
     GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
