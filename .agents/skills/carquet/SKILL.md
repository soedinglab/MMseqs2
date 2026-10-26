---
name: carquet
description: Integrate the carquet C11 Parquet library into a C or C++ codebase. Use when a project needs to read or write Apache Parquet files from C/C++ without Arrow, when wiring carquet into CMake / FetchContent / add_subdirectory / pkg-config / xmake, when choosing between carquet's reader and writer APIs (writer, batch reader, column reader, nested lists/maps, Arrow C Data bridge), when writing nullable or nested columns (definition/repetition levels, null bitmaps), when mapping C types to Parquet physical/logical types, when debugging carquet error codes, buffer lifetimes or link errors, or when checking a carquet-produced file with the `carquet` CLI, PyArrow or DuckDB.
---

# Integrating carquet

carquet is a pure C11 library (`#include <carquet/carquet.h>`, CMake target `carquet`) that reads and writes Apache Parquet with Snappy / GZIP / LZ4 / ZSTD / Brotli, statistics, bloom filters, page indexes, nested schemas and an Arrow C Data Interface bridge, with no Arrow dependency. This skill is the shortest path from "this project needs Parquet" to a verified integration. Read a reference file only when its step needs it:

| File | Read it when |
|---|---|
| [references/integration.md](references/integration.md) | wiring the build (CMake variants, reusing the project's zstd/zlib/lz4/brotli, pkg-config, xmake, C++, MSVC, LTO) |
| [references/api-map.md](references/api-map.md) | you need the exact name, arguments or purpose of any of the 131 public functions, or the fields of an options struct |
| [references/recipes.md](references/recipes.md) | writing the code: 13 complete patterns (flat/nullable/string tables, temporal and decimal types, lists and maps, projected scans, predicate and page filters, column reader with levels, append, buffer output, Arrow bridge, dictionary-preserving reads, per-column tuning, custom codecs, error reporting) |
| [references/pitfalls.md](references/pitfalls.md) | something misbehaves: symptom → cause → fix tables for writing, reading, build/link, Arrow and concurrency |
| [references/types.md](references/types.md) | choosing physical/logical types, C buffer types, encodings, codecs, Arrow format mapping |

## Procedure

### 1. Establish the project's shape (2 minutes, no code yet)

Answer these from the repository before writing anything; each one changes the integration:

- Build system: CMake (installed package, `FetchContent`, vendored subdirectory), plain Makefile / pkg-config, xmake, or Bazel/Meson (then link the static library and `include/` by hand).
- Does the project already build zstd, zlib, lz4 or brotli? If yes, pass `CARQUET_ZSTD_TARGET` / `CARQUET_ZLIB_TARGET` / `CARQUET_LZ4_TARGET` / `CARQUET_BROTLI_TARGET` so carquet links those instead of fetching copies. `CARQUET_ENABLE_BROTLI=OFF` drops the brotli dependency entirely (the codec then reports `UNSUPPORTED_CODEC`).
- C or C++ consumer, compiler family, and whether LTO objects from a different compiler would be linked (then build carquet with the same compiler or turn `CMAKE_INTERPROCEDURAL_OPTIMIZATION` off).
- Data shape: flat columns only, nullable columns, lists/maps, deeper nesting, or data already living in Arrow arrays. Approximate row count and whether files are read back by this project, by PyArrow/DuckDB/Spark, or both.
- Read pattern: full scans, column projection, point lookups / predicates, row windows, or metadata inspection only.

### 2. Wire the dependency

Use the matching recipe in [references/integration.md](references/integration.md). The invariants: target `carquet` (static by default; `CARQUET_BUILD_SHARED=ON` for shared), set `CARQUET_BUILD_CLI OFF` and `CARQUET_BUILD_DEV OFF` in embedded builds, CMake ≥ 3.16, compilers GCC / Clang / MSVC / MinGW, OpenMP optional. Prove the link before touching application code:

```c
#include <carquet/carquet.h>
#include <stdio.h>
int main(void) { printf("carquet %s\n", carquet_version()); return 0; }
```

### 3. Pick the API surface

| Need | Use | Recipe |
|---|---|---|
| Write flat columns (required or nullable) | `carquet_schema_add_column()` + `carquet_writer_create()` + `carquet_writer_write_batch()` | 1, 2 |
| Write `LIST<T>` / `MAP<K,V>` one level deep | `carquet_schema_add_list()` / `add_map()` + `carquet_writer_write_list_column()` (Arrow offsets + validity, no level maths) | 3 |
| Write deeper nesting, or data already in Arrow | `carquet_arrow_import_schema()` + `carquet_writer_write_arrow()` | 9 |
| Append to an existing file | `carquet_writer_open_append()` | 8 |
| Produce bytes in memory | `carquet_writer_create_buffer()` + `carquet_writer_get_buffer()` after `close()` | 8 |
| Scan rows, project columns, parallel decode | `carquet_batch_reader_create()` + `carquet_batch_reader_next()` + `carquet_row_batch_column()` | 4 |
| Skip data with predicates | automatic row-group pruning + `carquet_batch_reader_set_page_filter()` (file needs a page index) + `set_row_range()` | 5 |
| One column of one row group, raw levels, manual skipping | `carquet_reader_get_column()` + `carquet_column_read_batch_ex()` + `carquet_list_offsets()` | 6 |
| Single-level lists through the batch reader | `carquet_row_batch_column_list()` | 7 |
| Hand data to Arrow / DuckDB / PyArrow in-process | `carquet_arrow_export_batch()` (flat) or `carquet_reader_read_arrow()` (any depth) | 9 |
| Low-cardinality strings at speed | `preserve_dictionaries` + `carquet_row_batch_column_dictionary()` | 10 |
| Inspect without reading data | `carquet_get_file_info()`, `carquet_validate_file()`, `carquet_reader_column_statistics()`, `carquet_reader_column_chunk_metadata()` | 5, 13 |

Defaults worth knowing before tuning: PLAIN encoding (BYTE_STREAM_SPLIT for floats under a codec) and no dictionary, 128 MB row groups, 1 MB pages, statistics and CRC on, page index and bloom filters off, `async_io` on, ZSTD level 1 when ZSTD is chosen. For files meant to be queried later, turn on `write_page_index` and consider `write_bloom_filters` and `RLE_DICTIONARY` for low-cardinality columns.

### 4. Write the code against the three conventions that cause nearly every bug

1. **Nullable values are dense on write and in the column reader, row-aligned in the batch reader.** `carquet_writer_write_batch(w, col, values, num_values, def_levels, rep_levels)`: `num_values` is the logical row count, `def_levels` has one entry per row (`== max_def` present, lower = null), `values` holds only the present values packed contiguously; `carquet_column_read_batch()` returns the same dense layout plus levels. The batch reader instead returns one slot per logical row: `carquet_row_batch_column()` gives `num_values == rows`, null slots zero-filled (`{NULL, 0}` for BYTE_ARRAY), and a null bitmap indexed by row (bit set = present, LSB first, NULL = all present). Get `max_def` / `max_rep` from `carquet_schema_max_def_level()` / `max_rep_level()`, never hard-code them.
2. **Row counts must agree.** Every leaf column must have received the same number of logical rows before `carquet_writer_new_row_group()` and `carquet_writer_close()`. A failed write path ends with `carquet_writer_abort()`, never `close()`.
3. **Ownership is explicit.** The writer copies the schema (free it after `create`); reader-returned pointers (schema, metadata strings, statistics bytes) live until `carquet_reader_close()`; batch pointers live until `carquet_row_batch_free()`; `carquet_reader_open_buffer()` does not copy the buffer; `carquet_writer_get_buffer()` hands you a heap block to `free()`; the Arrow bridge **consumes** the `ArrowArray` / `ArrowSchema` you pass to `write_arrow` / `import_schema` (never release them again) and makes you the owner of what `read_arrow` / `export_*` return (release exactly once).

Everything returns `carquet_status_t` (`CARQUET_OK` is 0) or a pointer plus a filled `carquet_error_t` (`CARQUET_ERROR_INIT`; C++: `carquet_error_t err{};`). `carquet_batch_reader_next()` ends a scan with `CARQUET_ERROR_END_OF_DATA`; any other non-OK status is a failure whose detail is in `carquet_batch_reader_last_error()`. Format errors with `carquet_error_format()`.

### 5. Verify before declaring success

1. Build with warnings on; run the smallest round trip (write a handful of rows including a null and an empty string, read them back, compare values and null positions).
2. Inspect the file with the CLI that ships with the library: `carquet schema f.parquet`, `carquet info f.parquet`, `carquet validate f.parquet`, `carquet head -n 5 f.parquet`, `carquet stat f.parquet`.
3. Cross-check with another reader when the files leave the project: `python3 -c "import pyarrow.parquet as pq; print(pq.read_table('f.parquet'))"` or DuckDB `SELECT * FROM read_parquet('f.parquet')`.
4. For multi-threaded consumers confirm one handle per thread (handles are not shared; parallelism is internal via `num_threads` / OpenMP), and that `carquet_register_codec()` / `carquet_set_allocator()` run at startup.
5. If performance matters, measure first; then apply in this order: `use_mmap = true`, project columns, keep `batch_size` at 65536, `preserve_dictionaries` for string-heavy scans, page index + filters for selective reads, ZSTD level 1 or LZ4_RAW for the size/speed balance, larger row groups (tens of MB) so parallel encode and async I/O pay off. Details: `docs/performance.md` in the carquet tree.

### 6. When something fails

Match the symptom in [references/pitfalls.md](references/pitfalls.md) first; the top entries are sparse-vs-dense nullable values, mismatched row counts, kept pointers after free, missing page index for filters, double release of Arrow structs, and private codec libraries not linked outside CMake. Status codes that point at the caller: `INVALID_ARGUMENT` (shape, levels, options), `TYPE_MISMATCH` (buffer type), `INVALID_STATE` (row-count or lifecycle), `NOT_IMPLEMENTED` (unsupported nesting shape or Arrow type), `PAGE_INDEX_REQUIRED`, `INVALID_SCHEMA` (append). Codes that point at the file: `DECODE`, `INVALID_COMPRESSED_DATA`, `CRC_MISMATCH`, `INVALID_FOOTER`, `FILE_TRUNCATED`.

## Where things are in the carquet tree

`include/carquet/carquet.h` (whole API, Doxygen comments), `include/carquet/types.h` and `error.h`; `docs/` (manual: `reading.md`, `writing.md`, `nested-data.md`, `performance.md`, `error-handling.md`); `examples/` (one runnable program per topic); `mocklib/` (a complete downstream library using 117 public functions, with the three-way CMake resolution integrators should copy); `interop/` (PyArrow / DuckDB round-trip scripts); `CHANGELOG.md` (behaviour changes per version).
