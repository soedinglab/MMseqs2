set(compile_resource__internal_dir ${CMAKE_CURRENT_LIST_DIR} CACHE INTERNAL "")

find_program(XXD_EXECUTABLE xxd)
if(NOT XXD_EXECUTABLE)
    find_program(PERL_EXECUTABLE perl)
    if(NOT PERL_EXECUTABLE)
        message(FATAL_ERROR "Neither xxd nor perl found in PATH. xxd is usually contained in your distributions vim-common package!")
    else()
        message("-- xxd not found, using xxdi.pl instead")
        set(XXD_EXECUTABLE "${compile_resource__internal_dir}/xxdi.pl")
        set(XXD_PARAMS "")
    endif()
else()
    set(XXD_PARAMS -i)
endif()

find_program(SED_EXECUTABLE sed)
if(NOT SED_EXECUTABLE)
    message(FATAL_ERROR "sed not found in path!")
endif()

if(${HAVE_SHELLCHECK})
    find_program(SHELLCHECK_EXECUTABLE shellcheck)
    if(SHELLCHECK_EXECUTABLE)
        message("-- ShellCheck enabled")
    else()
        message("-- ShellCheck not found")
        set(SHELLCHECK_EXECUTABLE true)
    endif()
else()
    message("-- ShellCheck disabled")
    set(SHELLCHECK_EXECUTABLE true)
endif()

function(compile_resource INPUT_FILE OUTPUT_FILE)
    get_filename_component(INPUT_FILE_NAME ${PROJECT_SOURCE_DIR}/data/${INPUT_FILE} NAME)
    get_filename_component(INPUT_FILE_DIRECTORY ${PROJECT_SOURCE_DIR}/data/${INPUT_FILE} DIRECTORY)
    set(OUTPUT_FILE ${PROJECT_BINARY_DIR}/generated/${INPUT_FILE_NAME}.h)
    set(OUTPUT_FILE ${OUTPUT_FILE} PARENT_SCOPE)
    add_custom_command(OUTPUT ${OUTPUT_FILE}
            COMMAND ${compile_resource__internal_dir}/checkshell.sh ${SHELLCHECK_EXECUTABLE} ${INPUT_FILE_NAME}
            COMMAND mkdir -p ${PROJECT_BINARY_DIR}/generated
            COMMAND ${XXD_EXECUTABLE} ${XXD_PARAMS} ${INPUT_FILE_NAME} > ${OUTPUT_FILE}
            COMMAND ${SED_EXECUTABLE} 's!unsigned char!static const unsigned char!' < ${OUTPUT_FILE} > ${OUTPUT_FILE}.tmp
            COMMAND mv -f ${OUTPUT_FILE}.tmp ${OUTPUT_FILE}
            WORKING_DIRECTORY ${INPUT_FILE_DIRECTORY}
            DEPENDS ${PROJECT_SOURCE_DIR}/data/${INPUT_FILE})
    set_source_files_properties(${OUTPUT_FILE} PROPERTIES GENERATED TRUE)
endfunction()
