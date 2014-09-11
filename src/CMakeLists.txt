cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(mmseqs)

add_subdirectory(alignment)
add_subdirectory(clustering)
add_subdirectory(commons)
add_subdirectory(prefiltering)
add_subdirectory(util)
add_subdirectory(workflow)

set(mmseqs_source_files mmseqs.cpp )
add_executable(mmseqs ${mmseqs_source_files})

target_link_libraries(mmseqs util)
target_link_libraries(mmseqs workflow)
target_link_libraries(mmseqs ffindex)
