if(CMAKE_COMPILER_IS_ICC)
    set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -ipo -no-prec-div -xHost -fno-exceptions")
else()
    set(CMAKE_CXX_FLAGS_RELEASE        "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -ffast-math -fno-exceptions -ftree-vectorize -fno-strict-aliasing")
endif()

include_directories(generatedfiles)
add_subdirectory(alignment)
add_subdirectory(clustering)
add_subdirectory(commons)
add_subdirectory(prefiltering)
add_subdirectory(util)
add_subdirectory(workflow)
add_subdirectory(test)
add_subdirectory(version)

set(mmseqs_source_files
    CommandDeclarations.h
    mmseqs.cpp
)

add_executable(mmseqs ${mmseqs_source_files})

target_link_libraries(mmseqs util workflow version)
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

if(BZIP2_FOUND)
    target_link_libraries(mmseqs ${BZIP2_LIBRARIES})
endif()



install(TARGETS mmseqs DESTINATION bin)
