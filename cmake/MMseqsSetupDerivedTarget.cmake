include(AppendTargetProperty)

function (mmseqs_setup_derived_target TARGET)
    get_target_property(COMPILE_TMP mmseqs-framework COMPILE_FLAGS)
    get_target_property(LINK_TMP mmseqs-framework LINK_FLAGS)
    get_target_property(DEF_TMP mmseqs-framework COMPILE_DEFINITIONS)
    get_target_property(INCL_TMP mmseqs-framework INCLUDE_DIRECTORIES)

    target_link_libraries(${TARGET} mmseqs-framework)
    append_target_property(${TARGET} COMPILE_FLAGS ${COMPILE_TMP})
    append_target_property(${TARGET} LINK_FLAGS ${LINK_TMP})
    set_property(TARGET ${TARGET} APPEND PROPERTY COMPILE_DEFINITIONS ${DEF_TMP})
    set_property(TARGET ${TARGET} APPEND PROPERTY INCLUDE_DIRECTORIES ${INCL_TMP})
endfunction()

function (restore_exceptions TARGET)
    get_target_property(COMPILE_TMP ${TARGET} COMPILE_FLAGS)
    get_target_property(LINK_TMP ${TARGET} LINK_FLAGS)

    if(COMPILE_TMP MATCHES "-fno-exceptions")
        string(REPLACE "-fno-exceptions" "" COMPILE_TMP "${COMPILE_TMP}")
    endif()

    if(LINK_TMP MATCHES "-fno-exceptions")
        string(REPLACE "-fno-exceptions" "" LINK_TMP "${LINK_TMP}")
    endif()

    set_property(TARGET ${TARGET} PROPERTY COMPILE_FLAGS ${COMPILE_TMP})
    set_property(TARGET ${TARGET} PROPERTY LINK_FLAGS ${LINK_TMP})
endfunction()
