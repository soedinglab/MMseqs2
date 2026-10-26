set(compile_resource__internal_dir ${CMAKE_CURRENT_LIST_DIR} CACHE INTERNAL "")

find_program(XXD_EXECUTABLE xxd)
find_program(SED_EXECUTABLE sed)
if(XXD_EXECUTABLE AND SED_EXECUTABLE)
    set(RESOURCE_COMPILER_NATIVE 0 CACHE INTERNAL "")
else()
    message("-- xxd or sed not found, compiling resources with CMake instead")
    set(RESOURCE_COMPILER_NATIVE 1 CACHE INTERNAL "")
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
    if(RESOURCE_COMPILER_NATIVE)
        string(MAKE_C_IDENTIFIER ${INPUT_FILE_NAME} VAR_NAME)
        add_custom_command(OUTPUT ${OUTPUT_FILE}
                COMMAND ${compile_resource__internal_dir}/checkshell.sh ${SHELLCHECK_EXECUTABLE} ${INPUT_FILE_NAME}
                COMMAND ${CMAKE_COMMAND} -E make_directory ${PROJECT_BINARY_DIR}/generated
                COMMAND ${CMAKE_COMMAND} -DIN=${INPUT_FILE_NAME} -DOUT=${OUTPUT_FILE} -DNAME=${VAR_NAME} -P ${compile_resource__internal_dir}/ResourceToString.cmake
                WORKING_DIRECTORY ${INPUT_FILE_DIRECTORY}
                DEPENDS ${PROJECT_SOURCE_DIR}/data/${INPUT_FILE} ${compile_resource__internal_dir}/ResourceToString.cmake)
    else()
        add_custom_command(OUTPUT ${OUTPUT_FILE}
                COMMAND ${compile_resource__internal_dir}/checkshell.sh ${SHELLCHECK_EXECUTABLE} ${INPUT_FILE_NAME}
                COMMAND mkdir -p ${PROJECT_BINARY_DIR}/generated
                COMMAND ${XXD_EXECUTABLE} -i ${INPUT_FILE_NAME} > ${OUTPUT_FILE}
                COMMAND ${SED_EXECUTABLE} -f ${compile_resource__internal_dir}/tostring.sed < ${OUTPUT_FILE} > ${OUTPUT_FILE}.tmp
                COMMAND mv -f ${OUTPUT_FILE}.tmp ${OUTPUT_FILE}
                WORKING_DIRECTORY ${INPUT_FILE_DIRECTORY}
                DEPENDS ${PROJECT_SOURCE_DIR}/data/${INPUT_FILE} ${compile_resource__internal_dir}/tostring.sed)
    endif()
    set_source_files_properties(${OUTPUT_FILE} PROPERTIES GENERATED TRUE)
endfunction()
