set(COMPILED_RESOURCES
        blastp.sh
        blastpgp.sh
        clustering.sh
        cascaded_clustering.sh
        update_clustering.sh
        blosum30.out
        blosum35.out
        blosum40.out
        blosum45.out
        blosum50.out
        blosum55.out
        blosum60.out
        blosum62.out
        blosum65.out
        blosum70.out
        blosum75.out
        blosum80.out
        blosum85.out
        blosum90.out
        blosum95.out
        blosum100.out
        )

set(GENERATED_OUTPUT_HEADERS "")
FOREACH(INPUT_FILE ${COMPILED_RESOURCES})
    get_filename_component(INPUT_FILE_NAME ${INPUT_FILE} NAME)
    set(OUTPUT_FILE ${PROJECT_SOURCE_DIR}/src/generatedfiles/${INPUT_FILE_NAME}.h)
    message("-- Compile ${INPUT_FILE_NAME} to ${OUTPUT_FILE}")
    add_custom_command(OUTPUT ${OUTPUT_FILE}
            COMMAND xxd -i ${INPUT_FILE_NAME} ${OUTPUT_FILE}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data
            MAIN_DEPENDENCY ${INPUT_FILE_NAME})
    set_source_files_properties(${OUTPUT_FILE} PROPERTIES GENERATED TRUE)
    list(APPEND GENERATED_OUTPUT_HEADERS "${OUTPUT_FILE}")
ENDFOREACH()

add_custom_target(generated ALL DEPENDS ${GENERATED_OUTPUT_HEADERS})
