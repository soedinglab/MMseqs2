cmake_minimum_required(VERSION 3.15)
project(test_project VERSION 0.1.0)
include(../../test_header.cmake)

corrosion_import_crate(MANIFEST_PATH rust/Cargo.toml FLAGS --features one)

add_executable(flags-exe main.cpp)
target_link_libraries(flags-exe PUBLIC flags_lib)
corrosion_set_cargo_flags(flags_lib --features two)
corrosion_set_cargo_flags(flags_lib $<TARGET_PROPERTY:flags_lib,more_flags>)

set_property(
    TARGET flags_lib
    APPEND
    PROPERTY more_flags --features three
)
