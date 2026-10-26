# Performance and Tuning

## Start With the Cheap Wins

- Read with `carquet_batch_reader_t` unless you need raw level streams.
- Project only the columns you need.
- Enable `use_mmap` in `carquet_reader_options_t` for local file workloads.
- Leave `num_threads = 0` unless you have a reason to pin a specific thread count.

## mmap and Zero-Copy

`mmap` is the main read-side performance lever for local files. Check the actual mode with `carquet_reader_is_mmap()`.

Carquet can expose page data without copying, but only for a narrow fast path:

- mmap reader
- uncompressed column chunk
- PLAIN encoding
- fixed-width physical type
- no definition levels

Check eligibility with `carquet_reader_can_zero_copy()`.

In practice, this is most valuable for large required numeric columns and fixed-size binary data.

## Batch Size and Threads

Relevant knobs:

- `carquet_batch_reader_config_t.batch_size`
- `carquet_batch_reader_config_t.num_threads`
- `carquet_reader_options_t.num_threads`

Practical defaults:

- leave batch size at `65536` rows until you measure otherwise
- use larger batches for simple scans
- use smaller batches when downstream processing is heavy or latency-sensitive
- keep auto thread selection for wide or compressed datasets

## Preserve Dictionaries for String-Heavy Workloads

If a projected column is dictionary-encoded, `preserve_dictionaries = true` avoids materializing the final values and returns `uint32_t` indices instead.

That is often a large win for string-heavy scans because it skips per-row scatter/gather work.

Workflow:

1. Set `cfg.preserve_dictionaries = true`.
2. Read a batch.
3. Call `carquet_row_batch_column_dictionary()`.
4. If the column was not dictionary-preserved, the call returns `CARQUET_ERROR_INVALID_ARGUMENT`.

For `BYTE_ARRAY` dictionaries, use the returned `dictionary_offsets` table for O(1) lookup into `dictionary_data`.

## Prebuffer on Non-mmap Storage

`carquet_reader_prebuffer()` is useful when all of these are true:

- you are not using mmap
- you will read many columns from the same row group
- the storage is high-latency or seek-heavy

It coalesces file reads so later column readers can serve data from the prebuffered cache. On mmap readers it is a no-op.

## Writer Parallelism

With OpenMP, a row group with more than one column is encoded and compressed in parallel when it is finalized (the next `carquet_writer_new_row_group()` call, an automatic flush, or `carquet_writer_close()`): `carquet_writer_write_batch()` only stashes the input. Flat fixed-width columns (REQUIRED INT32/INT64/FLOAT/DOUBLE/INT96/FIXED_LEN_BYTE_ARRAY with PLAIN or BYTE_STREAM_SPLIT, no bloom filter or page index) are split into one task per page, so a few wide numeric columns still keep every core busy; other columns are one task each. The thread count follows `OMP_NUM_THREADS` / `omp_set_num_threads()`.

With `async_io` (on by default) the finished row group is written by a background thread while the next one is encoded. The kernel copy into the page cache is often the single largest cost of an uncompressed write, so overlapping it matters even without a codec. Keep row groups a few MB or larger: the overlap and the page-level split both need more than one page per column to pay off.

## Write Files for Future Reads

Writer settings change how much pruning future readers can do:

- `write_statistics = true`: enables row-group min/max pruning
- `write_page_index = true`: enables page-level pruning
- `write_bloom_filters = true`: enables fast membership checks

You can override settings per column before writing:

- `carquet_writer_set_column_encoding()`
- `carquet_writer_set_column_compression()`
- `carquet_writer_set_column_statistics()`
- `carquet_writer_set_column_bloom_filter()`

## Checksums, Initialization, and Diagnostics

- Keep `verify_checksums = true` when correctness matters more than peak throughput.
- Call `carquet_init()` explicitly if you want deterministic SIMD setup timing.
- Use `carquet_get_cpu_info()` for diagnostics or benchmark logs.
- Call `carquet_cleanup()` before process exit if you want clean valgrind-style reports.
