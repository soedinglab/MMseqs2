# Carquet Manual

This manual tracks the public API in `include/carquet/carquet.h`. It is meant to answer "which API should I use?" and "what lifetime or level rules do I need to get right?" without repeating every declaration from the header.

For runnable code, start with:

- [`examples/basic_write_read.c`](../examples/basic_write_read.c)
- [`examples/nullable_columns.c`](../examples/nullable_columns.c)
- [`examples/data_types.c`](../examples/data_types.c)
- [`examples/advanced_features.c`](../examples/advanced_features.c)

## Pick the Right Surface

| Task | Use |
| --- | --- |
| Scan rows for analytics, project a few columns, or parallelize reads | [`reading.md`](./reading.md) and `carquet_batch_reader_t` |
| Read one column from one row group, or work with raw definition/repetition levels | [`reading.md`](./reading.md) and `carquet_column_reader_t` |
| Write flat data to a file | [`writing.md`](./writing.md) and `carquet_writer_t` |
| Add row groups to an existing file | [`writing.md`](./writing.md#append-to-an-existing-file) and `carquet_writer_open_append()` |
| Write nested schemas or nullable/repeated data | [`nested-data.md`](./nested-data.md) |
| Exchange data with Arrow (PyArrow/DuckDB/nanoarrow) via the C Data Interface | [`reading.md`](./reading.md#export-to-arrow-c-data-interface) / [`writing.md`](./writing.md#import-from-arrow-c-data-interface) and `carquet_arrow_*` |
| Push performance on local files or string-heavy datasets | [`performance.md`](./performance.md) |
| Inspect footer metadata only | `carquet_get_file_info()` and `carquet_validate_file()` |
| Understand error codes, type mappings, or level conventions | [`error-handling.md`](./error-handling.md) |

## Rules Worth Remembering

- `carquet_init()` is optional. Call it explicitly only if you want deterministic startup or early CPU feature detection.
- `carquet_cleanup()` is mainly useful at process shutdown if you want a clean leak checker run.
- `carquet_writer_create()`, `carquet_writer_create_file()`, and `carquet_writer_create_buffer()` copy the schema. You can free the schema after writer creation.
- `carquet_reader_schema()` returns a pointer owned by the reader. It stays valid until `carquet_reader_close()`.
- `carquet_reader_open_buffer()` does not copy the source buffer. The buffer must stay alive and unchanged for the full reader lifetime.
- `carquet_writer_get_buffer()` transfers a heap buffer to the caller after `carquet_writer_close()`. Free it with `free()`.
- For nullable writes, `num_values` is the logical row count, but the `values` array contains only present values packed contiguously.
- Batch column pointers returned by `carquet_row_batch_column()` and `carquet_row_batch_column_dictionary()` should be treated as valid only until `carquet_row_batch_free()`.
- All columns in a writer must advance by the same logical row count before `carquet_writer_new_row_group()` or `carquet_writer_close()`.

## Manual Pages

- [`reading.md`](./reading.md): opening files, batch reads, low-level column reads, statistics, metadata inspection
- [`writing.md`](./writing.md): schema creation, required and nullable writes, row groups, metadata, buffer-backed output
- [`nested-data.md`](./nested-data.md): lists, maps, groups, definition levels, repetition levels, schema introspection
- [`performance.md`](./performance.md): mmap, zero-copy, dictionary-preserving reads, prebuffering, tuning knobs
- [`error-handling.md`](./error-handling.md): status codes, rich error context, type mapping, null bitmap and level conventions
