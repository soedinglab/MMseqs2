message(STATUS "Building CMake Generator for Corrosion - This may take a while")

set(generator_src "${CMAKE_CURRENT_BINARY_DIR}/legacy_generator_src")
set(generator_destination "${CMAKE_CURRENT_BINARY_DIR}/legacy_generator")
set(generator_build_quiet "")

file(MAKE_DIRECTORY "${generator_src}")
file(COPY src DESTINATION "${generator_src}")
if(Rust_VERSION VERSION_LESS "1.56")
    message(STATUS "Corrosion Generator: Using Compatibility lock file, due to rust version less than 1.56")
    file(COPY Compat.Cargo.lock Compat.Cargo.toml DESTINATION "${generator_src}")
    file(RENAME "${generator_src}/Compat.Cargo.lock" "${generator_src}/Cargo.lock")
    file(RENAME "${generator_src}/Compat.Cargo.toml" "${generator_src}/Cargo.toml")
else()
    file(COPY Cargo.lock Cargo.toml DESTINATION "${generator_src}")
endif()

# Using cargo install has the advantage of caching the build in the user .cargo directory,
# so likely the rebuild will be very cheap even after deleting the build directory.
execute_process(
        COMMAND ${CMAKE_COMMAND}
        -E env
        # If the Generator is built at configure of a project (instead of being pre-installed)
        # We don't want environment variables like `RUSTFLAGS` affecting the Generator build.
        --unset=RUSTFLAGS
        "CARGO_BUILD_RUSTC=${RUSTC_EXECUTABLE}"
        "${CARGO_EXECUTABLE}" install
        --path "."
        --root "${generator_destination}"
        --locked
        ${_CORROSION_QUIET_OUTPUT_FLAG}
        WORKING_DIRECTORY "${generator_src}"
        RESULT_VARIABLE generator_build_failed
)
if(generator_build_failed)
    message(FATAL_ERROR "Building CMake Generator for Corrosion - failed")
else()
    message(STATUS "Building CMake Generator for Corrosion - done")
endif()
set(host_executable_suffix "")
if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
    set(host_executable_suffix ".exe")
endif()

set(_CORROSION_GENERATOR_EXE
        "${generator_destination}/bin/corrosion-generator${host_executable_suffix}"
)

add_executable(Corrosion::Generator IMPORTED GLOBAL)
set_property(
        TARGET Corrosion::Generator
        PROPERTY IMPORTED_LOCATION "${_CORROSION_GENERATOR_EXE}")

if (CORROSION_DEV_MODE)
    # If you're developing Corrosion, you want to make sure to re-configure whenever the
    # generator changes.
    file(GLOB_RECURSE _RUST_FILES CONFIGURE_DEPENDS generator/src/*.rs)
    file(GLOB _CARGO_FILES CONFIGURE_DEPENDS generator/Cargo.*)
    set_property(
            DIRECTORY APPEND
            PROPERTY CMAKE_CONFIGURE_DEPENDS
            ${_RUST_FILES} ${_CARGO_FILES})
endif()

