find_program(XXD_EXECUTABLE xxd)
if(NOT XXD_EXECUTABLE)
    message(FATAL_ERROR "xxd not found in path. It is usually contained in your distributions vim-common package!")
endif()

find_program(SED_EXECUTABLE sed)
if(NOT SED_EXECUTABLE)
    message(FATAL_ERROR "sed not found in path!")
endif()

function(compile_resource INPUT_FILE OUTPUT_FILE)
    get_filename_component(INPUT_FILE_NAME ${PROJECT_SOURCE_DIR}/data/${INPUT_FILE} NAME)
    set(OUTPUT_FILE ${PROJECT_BINARY_DIR}/generated/${INPUT_FILE_NAME}.h PARENT_SCOPE)
    add_custom_command(OUTPUT ${PROJECT_BINARY_DIR}/generated/${INPUT_FILE_NAME}.h
            COMMAND mkdir -p ${PROJECT_BINARY_DIR}/generated
            COMMAND ${XXD_EXECUTABLE} -i ${INPUT_FILE} ${PROJECT_BINARY_DIR}/generated/${INPUT_FILE_NAME}.h
            COMMAND ${SED_EXECUTABLE} 's!unsigned char!static const unsigned char!' < ${PROJECT_BINARY_DIR}/generated/${INPUT_FILE_NAME}.h > ${PROJECT_BINARY_DIR}/generated/${INPUT_FILE_NAME}.h.tmp
            COMMAND mv -f ${PROJECT_BINARY_DIR}/generated/${INPUT_FILE_NAME}.h.tmp ${PROJECT_BINARY_DIR}/generated/${INPUT_FILE_NAME}.h
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/data/
            DEPENDS ${PROJECT_SOURCE_DIR}/data/${INPUT_FILE})
    set_source_files_properties(${PROJECT_BINARY_DIR}/generated/${INPUT_FILE_NAME}.h PROPERTIES GENERATED TRUE)
endfunction()
