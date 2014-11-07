cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

set(testalignment TestAlignment.cpp )
add_executable(test_alignment ${testalignment})

target_link_libraries(test_alignment util)
target_link_libraries(test_alignment workflow)
target_link_libraries(test_alignment ffindex)
