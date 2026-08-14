cmake_minimum_required(VERSION 3.18)

if (NOT APPLE)
    message(FATAL_ERROR "libmarv-mps (Metal/MPS GPU backend) requires macOS")
endif ()

enable_language(OBJCXX)

add_library(marv STATIC marv_mps.mm)
set_source_files_properties(marv_mps.mm PROPERTIES LANGUAGE OBJCXX)
target_compile_features(marv PUBLIC cxx_std_17)
target_compile_options(marv PRIVATE -fobjc-arc)
set_target_properties(marv PROPERTIES POSITION_INDEPENDENT_CODE ON)

find_library(METAL_FRAMEWORK Metal REQUIRED)
find_library(FOUNDATION_FRAMEWORK Foundation REQUIRED)
target_link_libraries(marv PUBLIC ${METAL_FRAMEWORK} ${FOUNDATION_FRAMEWORK})
