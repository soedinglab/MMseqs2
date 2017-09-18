find_program(XXD_EXECUTABLE xxd)
if(NOT XXD_EXECUTABLE)
    message(FATAL_ERROR "xxd not found in path. It is usually contained in your distributions vim-common package!")
endif()

find_program(SED_EXECUTABLE sed)
if(NOT SED_EXECUTABLE)
    message(FATAL_ERROR "sed not found in path!")
endif()

function(compile_resource INPUT_FILE OUTPUT_FILE)
    get_filename_component(INPUT_FILE_NAME ${INPUT_FILE} NAME)
    set(OUTPUT_FILE ${PROJECT_SOURCE_DIR}/src/generatedfiles/${INPUT_FILE_NAME}.h PARENT_SCOPE)
    add_custom_command(OUTPUT ${OUTPUT_FILE}
            COMMAND ${XXD_EXECUTABLE} -i ${INPUT_FILE_NAME} ${OUTPUT_FILE}
            COMMAND ${SED_EXECUTABLE} 's!unsigned char!static const unsigned char!' < ${OUTPUT_FILE} > ${OUTPUT_FILE}.tmp
            COMMAND mv -f ${OUTPUT_FILE}.tmp ${OUTPUT_FILE}
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data
            MAIN_DEPENDENCY ${INPUT_FILE_NAME})
    set_source_files_properties(${OUTPUT_FILE} PROPERTIES GENERATED TRUE)
endfunction()
