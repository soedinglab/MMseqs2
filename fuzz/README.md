# Carquet Fuzzing

Fuzz targets for testing carquet with random and malformed inputs using libFuzzer or AFL++.

## Quick Start

```bash
# Build and run the reader fuzzer for 5 minutes
python3 fuzz/run_fuzzer.py reader

# Run all fuzzers (5 min each)
python3 fuzz/run_fuzzer.py all

# Run for 1 hour with 4 parallel jobs
python3 fuzz/run_fuzzer.py all --time 3600 --jobs 4

# Just build
python3 fuzz/run_fuzzer.py build

# List targets and corpus status
python3 fuzz/run_fuzzer.py list

# Drop old crash artifacts before a fresh run
python3 fuzz/run_fuzzer.py reader --clean-artifacts
```

## Requirements

- **Clang** with libFuzzer support (Clang 6.0+)
- On macOS: `brew install llvm`
- On Linux: `apt install clang`

## Fuzz Targets

| Target | Modes | Coverage |
|--------|-------|----------|
| `fuzz_reader` | 1 | Full file reader: buffer open, schema, batch reader, row-group filters, dictionary-preserving batches, column reader, type names, bloom filter query, column/offset index parsing, key-value metadata, column chunk metadata |
| `fuzz_writer` | 1 | Flat-schema creation, all primitive physical types, valid logical annotations, all codecs, opt-in encodings, Data Page V2, Arrow schema metadata, sorting columns, per-column bloom config, buffer/file outputs, prebuffered readback |
| `fuzz_nested_write` | 1 | Nested write helper (`carquet_writer_write_list_column`): arbitrary offsets/validity/child values over LIST and MAP schemas, then readback through `carquet_row_batch_column_list` |
| `fuzz_arrow_read` | 1 | Nested Arrow reassembler (`carquet_reader_read_arrow`): reassembles every row group of arbitrary file bytes into a nested `ArrowArray` tree (struct/list/map, any depth) and releases it |
| `fuzz_reader_file` | 1 | Path / `FILE*` / mmap entry points on a temp file: the fread page-reader path (unreachable from `open_buffer`), `carquet_validate_file`, `carquet_get_file_info`, column skip/has_next/remaining, row-group statistics + `row_group_matches` / `filter_row_groups`, geospatial stats, Arrow field metadata / type refinements, bloom float probes, `list_offsets` / `count_rows`, `read_arrow` on the fread path, mmap batch reader + `carquet_batch_reader_last_error` |
| `fuzz_arrow_write` | 1 | Differential Arrow write roundtrip: valid Arrow arrays built from fuzz bytes (primitives incl. coerced int8/uint16/uint64, date32, timestamp, decimal128, fixed-size binary, large_utf8, dictionary, plus list / large_list / struct / map, validity bitmaps and slice offsets at every level) → `carquet_writer_write_arrow` with fuzz-chosen codec, page version, page size, per-leaf encodings, bloom filters, page index → `carquet_reader_read_arrow`; every value compared, bloom filters checked for false negatives |
| `fuzz_arrow_stream` | 1 | Arrow stream scan (`carquet_reader_export_arrow_stream`) over arbitrary file bytes, every batch validated against the stream schema (layout per format, buffer extents, offsets, dictionary codes): `ARROW:schema` restore (dictionary encoding, offset widening, views, fixed-size lists), a predicate on the first leaf on every other input (exact filter + row compaction), fuzz-chosen batch size, every batch re-exported through `carquet_arrow_ref_*` views and slices |
| `fuzz_compression` | 3 | Snappy/LZ4_RAW/GZIP/ZSTD/BROTLI/Hadoop-framed LZ4: decompress malformed, compress-decompress roundtrip, undersized buffer, Snappy per-reason diagnostics (`carquet_snappy_decompress_diag` + formatter) |
| `fuzz_encodings` | 18 | RLE, Delta INT32/INT64, Plain (bool/int32/int64/float/double), Dictionary (int32/int64/float/double), BSS (float/double/generic), Delta Length Byte Array, Delta Byte Array, RLE levels |
| `fuzz_thrift` | 8 | Primitives, struct parsing, containers (list/map), file metadata, page headers, encoder-decoder roundtrip, ColumnIndex parsing, OffsetIndex parsing |
| `fuzz_roundtrip` | 16 | Delta INT32/INT64, LZ4/Snappy/GZIP/ZSTD, BSS float/double/generic fixed-width payloads, RLE, Dictionary INT32, Plain INT32, Delta Length Byte Array, Delta Byte Array |

## Commands

```bash
# Fuzz a specific target
python3 fuzz/run_fuzzer.py <target> [--time N] [--jobs N]

# Run all targets
python3 fuzz/run_fuzzer.py all [--time N] [--jobs N]

# Build only
python3 fuzz/run_fuzzer.py build [--no-sanitizers]

# Clean old crash artifacts before fuzzing
python3 fuzz/run_fuzzer.py <target> --clean-artifacts

# List targets and corpus info
python3 fuzz/run_fuzzer.py list

# Minimize a corpus
python3 fuzz/run_fuzzer.py minimize <target>

# Generate coverage report
python3 fuzz/run_fuzzer.py coverage <target>
```

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--time N` | 300 | Seconds per target |
| `--jobs N` | 1 | Parallel fuzzing jobs |
| `--no-sanitizers` | off | Disable ASan/UBSan |
| `--build-dir DIR` | `build-fuzz` | Build directory |
| `--clean-artifacts` | off | Delete saved `artifacts_<target>/` directories before fuzzing |

Extra libFuzzer flags can be passed after `--`:
```bash
python3 fuzz/run_fuzzer.py reader --time 600 -- -max_len=4096 -only_ascii=1
```

## Building Manually

### libFuzzer (recommended)

```bash
cmake -B build-fuzz \
    -DCMAKE_C_COMPILER=clang \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_C_FLAGS="-fsanitize=address,undefined -fno-omit-frame-pointer" \
    -DCARQUET_BUILD_FUZZ=ON \
    -DCARQUET_BUILD_TESTS=OFF

cmake --build build-fuzz
```

### AFL++

```bash
CC=afl-clang-fast cmake -B build-afl \
    -DCARQUET_BUILD_FUZZ=ON \
    -DCARQUET_FUZZ_ENGINE=AFL \
    -DCMAKE_BUILD_TYPE=Debug

cmake --build build-afl
afl-fuzz -i corpus_reader -o findings -- ./build-afl/fuzz/fuzz_reader_afl @@
```

## Notes

- **Coverage instrumentation**: `fuzz/CMakeLists.txt` compiles libcarquet itself with `-fsanitize=fuzzer-no-link`. Before Oct 2026 only the harness was instrumented, so libFuzzer saw a few dozen edges and fuzzed blind (e.g. `fuzz_reader` stalled at cov 70; instrumented it starts above 2000).
- **LeakSanitizer suppressions**: ZSTD / zlib contexts cached per thread in pthread TSD look leaked to LSan on macOS (it does not scan the TSD of live OpenMP threads). `run_fuzzer.py` sets `LSAN_OPTIONS=suppressions=fuzz/lsan.supp`; pass it yourself when running a target by hand.
- Checked-in seeds live in `fuzz/corpus/<target>/` or `fuzz/corpus/fuzz_<target>/`; both are picked up.
- File-level regressions for fuzzer findings go in `tests/test_fuzz_regressions.c`.

## Crash Analysis

When a crash is found:

1. Artifacts are saved to `build-fuzz/fuzz/artifacts_<target>/`
2. Reproduce: `./build-fuzz/fuzz/fuzz_<target> <artifact>`
3. Stack trace: `ASAN_OPTIONS=symbolize=1 ./build-fuzz/fuzz/fuzz_<target> <artifact>`
