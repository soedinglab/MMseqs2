cmake_minimum_required(VERSION 3.15)
project(test_project VERSION 0.1.0)
include(../../test_header.cmake)

corrosion_import_crate(
    MANIFEST_PATH ${CMAKE_CURRENT_SOURCE_DIR}/Cargo.toml
    CRATES member1 member2
    IMPORTED_CRATES imported_crate_list
)

#NOTE: member3 also contains a binary called my_program, but that shouldn't be a problem since it is not imported
add_executable(my_program main.cpp)
target_link_libraries(my_program PUBLIC member1 member2)

# Test that the list of imported crates matches our expectations.
if(NOT DEFINED imported_crate_list)
    message(FATAL_ERROR "Corrosion failed to set the variable passed via IMPORTED_CRATES.")
endif()
set(expected_crates member1 member2)
foreach(crate ${expected_crates})
    if(NOT "${crate}" IN_LIST imported_crate_list)
        message(FATAL_ERROR "Expected ${crate} to be imported, but it wasn't. Imported crate list:\n"
            "${imported_crate_list}"
        )
    endif()
endforeach()
set(additional_crates ${imported_crate_list})
list(REMOVE_ITEM additional_crates ${expected_crates})
if(additional_crates)
    message(FATAL_ERROR "Corrosion unexpectedly imported the following crates: ${additional_crates}")
endif()

