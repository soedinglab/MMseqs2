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
