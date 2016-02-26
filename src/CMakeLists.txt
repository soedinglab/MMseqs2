# openmp integration
find_package(OpenMP)
if(OPENMP_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOPENMP=1 ${OpenMP_CXX_FLAGS}")
endif()

#SSE
if(${HAVE_AVX2})
    add_definitions("-DAVX2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -Wa,-q")
elseif(${HAVE_SSE4_2})
    add_definitions("-DSSE")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse4.2")
else()
    include(CheckSSEFeatures)
    add_definitions(${CXX_DFLAGS})
    if(${HAVE_AVX2_EXTENSIONS})
        add_definitions("-DAVX2")
    # debugging
    #   add_definitions("-DSSE")
    else()
    if(${HAVE_SSE4_2_EXTENSIONS})
        add_definitions("-DSSE")
    else()
        message( FATAL_ERROR "At least SSE4.2 is needed to compile!" )
    endif(${HAVE_SSE4_2_EXTENSIONS})
    endif(${HAVE_AVX2_EXTENSIONS})
endif()

# set flags
set(CMAKE_CXX_FLAGS                "${CMAKE_CXX_FLAGS} -pedantic -Wall -Winline -Wdisabled-optimization -Wno-vla -Wno-sign-compare -Wno-write-strings -std=c++0x -m64")
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
    endif()
endif()
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
add_executable(mmseqs ${mmseqs_source_files})

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

install(TARGETS mmseqs
		DESTINATION bin)

# copy shell scripts to bin
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/workflow/clustering.sh
     DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
     GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/workflow/cascaded_clustering.sh
     DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
     GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/workflow/blastp.sh
     DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
     GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/workflow/blastpgp.sh
     DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
     GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
