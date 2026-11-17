# Pitfalls: symptom, cause, fix

Ordered by how often they bite integrators. When a user reports one of the symptoms in the left column, check the cause before reading any further.

## Writing

| Symptom | Cause | Fix |
|---|---|---|
| Nullable column reads back with garbage or shifted values; PyArrow shows wrong numbers next to nulls | `values` passed as a sparse array (one slot per row, placeholders for nulls) | Pack only the present values contiguously; `num_values` stays the logical row count; one `def_level` per row (`max_def` = present) |
| `carquet_writer_close()` or `new_row_group()` returns `CARQUET_ERROR_INVALID_STATE` / `INVALID_ARGUMENT` | Columns have different logical row counts | Every leaf column must receive the same number of rows (nulls included) before a row group ends; interleave `write_batch` calls per column |
| `write_batch` returns `CARQUET_ERROR_TYPE_MISMATCH` or values look byte-swapped | Buffer element type does not match the physical type (e.g. `int` for INT64, `bool` for BOOLEAN packed as bits, `float` for DOUBLE) | Use the C type of the physical type: `int32_t`, `int64_t`, `float`, `double`, `uint8_t` 0/1 per boolean, `carquet_byte_array_t`, packed `uint8_t[type_length]` |
| Strings written as empty / crash on write | `carquet_byte_array_t.data == NULL` with `length > 0`, or `length` negative | `data` may be NULL only when `length == 0`; `length` is `int32_t` |
| Per-column encoding/compression/bloom setting has no effect | Override called after the first `write_batch` | Call `carquet_writer_set_column_*()` after `create` and before any write |
| File much larger than expected | Default encoding is PLAIN (fast reads), no dictionary | Opt in: `carquet_writer_set_column_encoding(w, i, CARQUET_ENCODING_RLE_DICTIONARY)` for low-cardinality columns, `DELTA_BINARY_PACKED` for sorted ints, ZSTD level 1 |
| Readers (Spark / DuckDB) cannot prune row groups or pages | Statistics or page index disabled, or BYTE_ARRAY min/max truncated at 32 bytes by design | Keep `write_statistics = true` (default), set `write_page_index = true`, consider `write_bloom_filters` for point lookups |
| Output file invalid after an error | `carquet_writer_close()` called after a failed write, or never called | On failure call `carquet_writer_abort()`; `close()` only on the success path. After `abort()` the file is not a valid Parquet file |
| Writing from a worker thread crashes with `carquet_writer_create_file(FILE*)` | `async_io` writes from a background thread while the caller also uses the `FILE*` | Set `opts.async_io = false`, or do not touch the `FILE*` until `close()` |
| `carquet_writer_open_append()` returns `CARQUET_ERROR_INVALID_SCHEMA` | Supplied schema differs from the file's leaves (name, physical type or repetition) | Rebuild the schema from the file (`carquet_reader_schema()` introspection) and pass the same shape |
| Timestamps off by 1000× in PyArrow | Unit in `params.timestamp.unit` does not match the integers written | MICROS means microseconds since epoch in an `int64_t`; or let `coerce_timestamps` rescale |
| DECIMAL rejected by `add_column` | Precision too large for the physical type | INT32 ≤ 9 digits, INT64 ≤ 18, FLBA(N) ≤ floor((8N−1)·log10 2), BYTE_ARRAY unlimited |
| `write_list_column` returns `INVALID_ARGUMENT` | Null list on a REQUIRED list, null element on a REQUIRED leaf, or offsets not starting at 0 / decreasing | Declare the list / element OPTIONAL where nulls occur; offsets must be non-decreasing from 0 |
| Peak memory higher than one row group | `async_io` keeps one encoded row group in flight; the pipelined encode holds the pending row group's input while the next one fills | Lower `row_group_size`, or set `async_io = false` |

## Reading

| Symptom | Cause | Fix |
|---|---|---|
| Loop ends early with no visible error | The `while (next() == OK && batch)` idiom swallows the status | Keep the status; anything other than `CARQUET_OK` / `CARQUET_ERROR_END_OF_DATA` is an error: read `carquet_batch_reader_last_error()` |
| Values of an OPTIONAL column misaligned with rows | Mixing up the two layouts | Batch reader (`carquet_row_batch_column`): one slot per row, null slots zero-filled, bitmap indexed by row. Column reader (`carquet_column_read_batch`) and writer: dense present-only values plus one `def_level` per row |
| Crash or stale data after `carquet_row_batch_free()` | Pointers from `carquet_row_batch_column*()` kept past the free | Copy out before freeing the batch; the next `next()` also recycles buffers |
| Crash after `carquet_reader_close()` | Schema node pointers or metadata strings kept past close | Everything returned by reader metadata APIs is owned by the reader |
| `carquet_reader_open_buffer()` reader returns corrupt data | Source buffer freed or modified while the reader was alive | The reader does not copy the buffer; keep it alive and immutable until close |
| `next()` returns `CARQUET_ERROR_PAGE_INDEX_REQUIRED` | Page filter on a column whose chunk has no offset index | Write with `write_page_index = true`, or filter manually per row. A file with an offset index but no column index (statistics disabled) is read in full, not an error |
| `carquet_row_batch_column()` returns `INVALID_ARGUMENT` for a list column | Repeated columns use the list accessor | `carquet_row_batch_column_list()`; `max_rep_level > 1` needs `carquet_reader_read_arrow()` or the column reader |
| `next()` returns `NOT_IMPLEMENTED` | Projected a list-of-list column, or combined a page filter / row range with repeated columns | Use `carquet_reader_read_arrow()` for deep nesting; filters and row ranges are for flat projections |
| `carquet_row_batch_column_dictionary()` returns `INVALID_ARGUMENT` | Column not dictionary-encoded in this row group, or `preserve_dictionaries` not set | Fall back to `carquet_row_batch_column()`; dictionary preservation is per row group |
| `carquet_column_read_batch()` returns fewer values than expected, no error | Count is in level entries; a mid-batch error is indistinguishable from end of column | Use `carquet_column_read_batch_ex()` with an error out-parameter; the values array is dense |
| INT96 column cannot be filtered | INT96 has no defined sort order | Filter after reading, or write TIMESTAMP(INT64) |
| Statistics `min_value` for BYTE_ARRAY look like a `carquet_byte_array_t` | They are raw bytes plus `min_value_size` | Print with `%.*s`; for numeric types `memcpy` into the C type (unaligned) |
| `carquet_reader_column_statistics()` has `has_min_max == false` | INT96, GEOMETRY/GEOGRAPHY, statistics disabled, or an all-NaN FLOAT/DOUBLE chunk | Treat as "unknown": read the data |
| Zero-copy expected but `carquet_reader_can_zero_copy()` is false | Needs mmap + uncompressed + PLAIN + fixed width + no definition levels | Write that column REQUIRED, PLAIN, UNCOMPRESSED; open with `use_mmap = true` |
| Reading is slow in small batches on an OPTIONAL column | Fixed in v0.7.2 (dense cursor); older versions recounted levels per call | Upgrade, or read with larger `max_values` |

## Build and link

| Symptom | Cause | Fix |
|---|---|---|
| `undefined reference to ZSTD_*` / `inflate` / `LZ4_*` when linking a static carquet | Private dependencies not propagated outside CMake | Link with CMake's `carquet` target (propagates them), or `pkg-config --static --libs carquet`, or add `-lzstd -lz -llz4` (+ OpenMP runtime) by hand |
| `undefined reference to omp_*` / `GOMP_*` | carquet built with OpenMP, consumer not linking the runtime | Use the CMake target (adds `OpenMP::OpenMP_C`), or add `-fopenmp` / `-lomp` |
| LTO link error ("bitcode", "unknown file format", `ld: archive member not an object`) | carquet built in Release with LTO by one compiler, consumer linked with another | Build carquet with the consumer's compiler, or `-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF` |
| `find_package(carquet)` fails | Not installed, or non-default prefix | `cmake --install build --prefix <p>` then `-DCMAKE_PREFIX_PATH=<p>`; or `FetchContent` / `add_subdirectory` |
| Configure fetches zstd/zlib/lz4 even though the project builds them | carquet searches system libs first, then fetches | Set `CARQUET_ZSTD_TARGET` / `CARQUET_ZLIB_TARGET` / `CARQUET_LZ4_TARGET` to the project's targets before adding carquet |
| MSVC: C11 features rejected | Old toolset | MSVC 2015+ builds carquet; consumers compile their own code with their own standard, the header needs only C99-level syntax |
| `carquet` CLI binary built and installed although only the library is wanted | `CARQUET_BUILD_CLI` defaults to ON | `-DCARQUET_BUILD_CLI=OFF` |

## Arrow bridge

| Symptom | Cause | Fix |
|---|---|---|
| Double free / crash after `carquet_writer_write_arrow()` | Caller released the `ArrowArray` / `ArrowSchema` after the call | Both are consumed (released) by carquet, on success and failure; never release them again |
| Leak reported after `carquet_reader_read_arrow()` / `export_batch()` | Caller did not call `release` | Call `out_array.release(&out_array)` and `out_schema.release(&out_schema)` exactly once |
| `NOT_IMPLEMENTED` on import | Union or run-end-encoded Arrow array, or unsupported dictionary value type | Convert on the Arrow side first; dictionaries with primitive values, sliced arrays, fixed-size lists, decimals, date64 / seconds timestamps are handled |
| `INVALID_ARGUMENT`: leaf count mismatch | Writer schema and Arrow schema differ in leaves (struct children count as leaves) | Build the writer schema with `carquet_arrow_import_schema()` from the same Arrow schema |
| `export_batch` returns `NOT_IMPLEMENTED` | Projection used, STRUCT column, or nesting deeper than one level | Export from an unprojected batch, or use `carquet_reader_read_arrow()` |

## Concurrency and lifetime summary

- One reader / writer / batch reader handle per thread at a time; different handles on different threads are independent. Internal parallelism is controlled by `num_threads` and OpenMP.
- `carquet_register_codec()` and `carquet_set_allocator()` are process-wide: call them at startup, before any reader or writer exists.
- `carquet_init()` is optional (lazy on first use); `carquet_cleanup()` frees per-thread codec contexts and is only needed for leak-checker-clean shutdown.
