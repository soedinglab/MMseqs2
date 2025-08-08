cmake_minimum_required(VERSION 3.15)
project(RustupProxy LANGUAGES CXX)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../../../cmake" ${CMAKE_MODULE_PATH})

function(_assert_is_rustup_proxy executable_path)
    execute_process(
        COMMAND
            ${CMAKE_COMMAND} -E env
                RUSTUP_FORCE_ARG0=rustup
            "${executable_path}" --version
        OUTPUT_VARIABLE _VERSION_RAW
        ERROR_VARIABLE _VERSION_STDERR
        RESULT_VARIABLE _VERSION_RESULT
    )

    if(NOT _VERSION_RESULT EQUAL "0")
        message(FATAL_ERROR "`${executable_path} --version` failed with ${_VERSION_RESULT}\n"
            "stderr:\n${_VERSION_STDERR}"
        )
    endif()

    if (NOT _VERSION_RAW MATCHES "rustup [0-9\\.]+")
        message(FATAL_ERROR "`${executable_path} --version` output does not match rustup: ${_VERSION_RAW}\n")
    endif()
endfunction()

set(Rust_RESOLVE_RUSTUP_TOOLCHAINS OFF CACHE BOOL "" FORCE)
find_package(Rust REQUIRED)

if (NOT Rust_FOUND)
    message(FATAL_ERROR "Rustup not found")
endif()

get_property(
    RUSTC_EXECUTABLE
    TARGET Rust::Rustc PROPERTY IMPORTED_LOCATION
)

_assert_is_rustup_proxy(${RUSTC_EXECUTABLE})

get_property(
    CARGO_EXECUTABLE
    TARGET Rust::Cargo PROPERTY IMPORTED_LOCATION
)

_assert_is_rustup_proxy(${CARGO_EXECUTABLE})
