function(mmseqs_setup_test NAME)
    string(TOLOWER ${NAME} BASE_NAME)
    string(REGEX REPLACE "\\.[^.]*$" "" BASE_NAME ${BASE_NAME})
    string(REGEX REPLACE "^test" "test_" BASE_NAME ${BASE_NAME})
    add_executable(${BASE_NAME} ${NAME})

    target_link_libraries(${BASE_NAME} mmseqs-framework)
endfunction()
