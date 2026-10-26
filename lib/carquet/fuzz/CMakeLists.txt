# Fuzzing targets for carquet
#
# Build with:
#   cmake -B build-fuzz -DCMAKE_C_COMPILER=clang -DCARQUET_BUILD_FUZZ=ON \
#         -DCMAKE_BUILD_TYPE=Debug
#   cmake --build build-fuzz
#
# Or for AFL++:
#   CC=afl-clang-fast cmake -B build-afl -DCARQUET_BUILD_FUZZ=ON \
#         -DCARQUET_FUZZ_ENGINE=AFL -DCMAKE_BUILD_TYPE=Debug
#   cmake --build build-afl

cmake_minimum_required(VERSION 3.16)

# Fuzzing options
option(CARQUET_FUZZ_SANITIZERS "Enable sanitizers for fuzzing" ON)
set(CARQUET_FUZZ_ENGINE "libFuzzer" CACHE STRING "Fuzzing engine: libFuzzer or AFL")

# Sanitizer flags
if(CARQUET_FUZZ_SANITIZERS)
    set(FUZZ_SANITIZER_FLAGS "-fsanitize=address" "-fsanitize=undefined" "-fno-omit-frame-pointer")
endif()

# Function to add a fuzz target
function(add_fuzz_target TARGET_NAME SOURCE_FILE)
    add_executable(${TARGET_NAME} ${SOURCE_FILE})
    target_link_libraries(${TARGET_NAME} PRIVATE carquet)
    target_include_directories(${TARGET_NAME} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
    )

    if(CARQUET_FUZZ_ENGINE STREQUAL "libFuzzer")
        target_compile_options(${TARGET_NAME} PRIVATE -fsanitize=fuzzer ${FUZZ_SANITIZER_FLAGS})
        target_link_options(${TARGET_NAME} PRIVATE -fsanitize=fuzzer ${FUZZ_SANITIZER_FLAGS})
    elseif(CARQUET_FUZZ_ENGINE STREQUAL "AFL")
        target_compile_definitions(${TARGET_NAME} PRIVATE AFL_MAIN)
        if(CARQUET_FUZZ_SANITIZERS)
            target_compile_options(${TARGET_NAME} PRIVATE ${FUZZ_SANITIZER_FLAGS})
            target_link_options(${TARGET_NAME} PRIVATE ${FUZZ_SANITIZER_FLAGS})
        endif()
    endif()
endfunction()

# Add fuzz targets
# Note: fuzz_reader covers all code paths through the public API
# fuzz_compression and fuzz_thrift test internal decoders directly
add_fuzz_target(fuzz_reader fuzz_reader.c)
add_fuzz_target(fuzz_compression fuzz_compression.c)
# Encodings fuzzer tests RLE, delta, and dictionary decoders
add_fuzz_target(fuzz_encodings fuzz_encodings.c)

# Thrift fuzzer tests the Thrift compact protocol decoder
add_fuzz_target(fuzz_thrift fuzz_thrift.c)

# Roundtrip fuzzer tests encode->decode consistency
add_fuzz_target(fuzz_roundtrip fuzz_roundtrip.c)

# Writer fuzzer tests schema creation and data writing
add_fuzz_target(fuzz_writer fuzz_writer.c)

# Page-filter fuzzer (v0.6.0): exercises the page-level predicate evaluator,
# which compares attacker-controlled column-index/offset-index stats against
# caller predicate values.
add_fuzz_target(fuzz_page_filter fuzz_page_filter.c)

# Append fuzzer (v0.6.0): exercises carquet_writer_open_append parsing the
# trailing footer of an untrusted existing file, schema-match validation, and
# the close-time footer rewrite over restored row-group metadata.
add_fuzz_target(fuzz_append fuzz_append.c)

# Arrow-schema fuzzer: exercises the FlatBuffer reader that recovers per-field
# custom_metadata (variable labels) from an untrusted "ARROW:schema" blob.
add_fuzz_target(fuzz_arrow_schema fuzz_arrow_schema.c)

# Arrow C Data Interface fuzzer: exercises the format-string parser reached via
# carquet_arrow_import_schema with attacker-controlled Arrow `format` strings.
add_fuzz_target(fuzz_arrow_c_data fuzz_arrow_c_data.c)

# Nested write fuzzer: drives carquet_writer_write_list_column (auto-shredding)
# with arbitrary offsets/validity/child values over LIST and MAP schemas.
add_fuzz_target(fuzz_nested_write fuzz_nested_write.c)

# Nested Arrow reassembler fuzzer: reads arbitrary bytes as a Parquet file and
# reassembles each row group into a nested ArrowArray via carquet_reader_read_arrow.
add_fuzz_target(fuzz_arrow_read fuzz_arrow_read.c)

# Create corpus directories
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/corpus_reader)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/corpus_compression)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/corpus_encodings)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/corpus_thrift)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/corpus_roundtrip)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/corpus_writer)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/corpus_page_filter)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/corpus_append)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/corpus_arrow_c_data)
