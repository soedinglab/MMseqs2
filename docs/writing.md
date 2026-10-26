# Writing Files

## Build the Schema First

The writer is column-oriented. Define the schema, then write one leaf column at a time.

```c
carquet_error_t err = CARQUET_ERROR_INIT;
carquet_schema_t* schema = carquet_schema_create(&err);
if (!schema) {
    fprintf(stderr, "%s\n", err.message);
    return 1;
}

carquet_logical_type_t string_type = { .id = CARQUET_LOGICAL_STRING };

carquet_schema_add_column(schema, "id",
    CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
carquet_schema_add_column(schema, "name",
    CARQUET_PHYSICAL_BYTE_ARRAY, &string_type, CARQUET_REPETITION_REQUIRED, 0, 0);
carquet_schema_add_column(schema, "score",
    CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);
```

For nested schemas, see [`nested-data.md`](./nested-data.md).

### Newer Logical Types

Modern logical annotations are passed through the same schema APIs as older types. `GEOMETRY` and `GEOGRAPHY` use `BYTE_ARRAY` storage, and `VARIANT` has a helper that creates the standard unshredded group layout:

```c
carquet_logical_type_t geometry_type = { .id = CARQUET_LOGICAL_GEOMETRY };
carquet_logical_type_t geography_type = {
    .id = CARQUET_LOGICAL_GEOGRAPHY,
    .params.geography = {
        .algorithm = CARQUET_GEOSPATIAL_EDGE_SPHERICAL,
        .has_algorithm = true,
    },
};

carquet_schema_add_column(schema, "shape",
    CARQUET_PHYSICAL_BYTE_ARRAY, &geometry_type,
    CARQUET_REPETITION_OPTIONAL, 0, 0);
carquet_schema_add_column(schema, "region",
    CARQUET_PHYSICAL_BYTE_ARRAY, &geography_type,
    CARQUET_REPETITION_OPTIONAL, 0, 0);

int32_t payload_group = carquet_schema_add_variant(
    schema, "payload", CARQUET_REPETITION_OPTIONAL, 0);
```

When a modern logical type has a legacy equivalent, the writer also emits the corresponding `ConvertedType` compatibility annotation automatically. Footer `column_orders` metadata is emitted automatically for written leaf columns.

## Create the Writer

```c
carquet_writer_options_t opts;
carquet_writer_options_init(&opts);
opts.compression = CARQUET_COMPRESSION_ZSTD;
opts.compression_level = 1;
opts.row_group_size = 128 * 1024 * 1024;
opts.page_size = 1 * 1024 * 1024;
opts.write_statistics = true;
opts.write_crc = true;
opts.write_page_index = true;
opts.write_bloom_filters = true;
opts.write_arrow_schema = false;  /* embed Arrow IPC schema as ARROW:schema (opt-in) */
opts.data_page_version = 1;       /* 1 = DATA_PAGE (default); 2 = DATA_PAGE_V2 */
opts.file_format_version = 2;     /* FileMetaData.version; 1 for very old readers */
opts.coerce_timestamps = false;   /* rescale all TIMESTAMP columns to one unit */
opts.write_batch_size = 0;        /* 0 = automatic internal batching */
opts.async_io = true;             /* write row groups on a background thread */

carquet_writer_t* writer = carquet_writer_create("out.parquet", schema, &opts, &err);
if (!writer) {
    fprintf(stderr, "%s\n", err.message);
    carquet_schema_free(schema);
    return 1;
}

carquet_schema_free(schema);  /* Safe after writer creation */
```

Two opt-in options, both off by default so default output bytes are unchanged:

- `write_arrow_schema`: when `true`, embeds the original Arrow schema as a base64-encoded Arrow IPC Schema message under the `ARROW:schema` footer key, so Arrow/PyArrow recover Arrow-specific type information losslessly. Both flat and nested schemas are emitted — a Parquet LIST becomes Arrow `List<element>`, a MAP becomes `Map<entries: struct<key, value>>`, and a plain group becomes a `Struct`. Emitted only when you have not already set that key yourself; any element that cannot be faithfully mapped aborts emission so the blob never disagrees with the Parquet schema.
- `data_page_version`: `2` writes `DATA_PAGE_V2` (repetition/definition levels stored uncompressed and outside the compressed value region, matching parquet-cpp); any other value keeps the default V1 path.
- `coerce_timestamps` / `coerce_timestamp_unit` / `allow_timestamp_truncation`: when `coerce_timestamps` is true, every `TIMESTAMP` column is rescaled to `coerce_timestamp_unit` (and its metadata emitted at that unit) regardless of the unit declared in the schema — the equivalent of PyArrow's `coerce_timestamps`. A coarser target loses precision and is rejected unless `allow_timestamp_truncation` is true (PyArrow's `allow_truncated_timestamps`).
- `write_batch_size`: caps how many values are processed per internal chunk before a page flush is considered (PyArrow's `write_batch_size`); `0` keeps the automatic page-size-derived heuristic. This tunes streaming/memory behavior, not the output format.
- `file_format_version`: the value written into `FileMetaData.version` in the footer. Default `2`; set to `1` for the very small set of historical readers that reject version-2 files. Independent of `data_page_version` (which controls page format) and of carquet's always-compatible metadata (modern `LogicalType` + legacy `ConvertedType`); any value other than `1` is treated as `2`.
- `async_io`: when `true` (the default from `carquet_writer_options_init()`), each finalized row group is written by a background thread while the next one is buffered, encoded and compressed, so file writes overlap with encoding. At most one row group is in flight (peak memory grows by about one compressed row group). A failed background write surfaces from the next row-group flush or from `carquet_writer_close()`. Output bytes are identical either way; set it to `false` to keep every write on the calling thread (for example when the `FILE*` handed to `carquet_writer_create_file()` must not be touched from another thread).

Other writer entry points:

- `carquet_writer_create_file(FILE*, ...)`
- `carquet_writer_create_buffer(...)`

## Per-Field Metadata (Variable Labels / Descriptions)

To attach a human-readable label or description to a column — the kind of
annotation survey/statistical tools carry — use `carquet_schema_set_field_metadata()`.
It records Arrow's `Field.custom_metadata`, the standard, reader-agnostic home
for schema-level annotations (it is *not* the per-row-group
`ColumnMetaData.key_value_metadata`). With `write_arrow_schema` enabled, the
pairs are emitted into the file-level `ARROW:schema` footer blob, so PyArrow,
Parquet viewers, and any Arrow-compatible reader surface them automatically.

```c
carquet_schema_add_column(schema, "Sex", CARQUET_PHYSICAL_INT32,
                          NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
/* Address the column just added by its schema element index. */
int32_t idx = carquet_schema_num_elements(schema) - 1;
carquet_schema_set_field_metadata(schema, idx, "Label", "Sex of Respondent");

carquet_writer_options_t opts;
carquet_writer_options_init(&opts);
opts.write_arrow_schema = true;   /* required for the labels to be emitted */
```

Setting the same key twice replaces the value; distinct keys accumulate. Only
flat top-level fields are emitted, matching the `ARROW:schema` writer. Read the
labels back with the [reader's column-metadata accessors](reading.md#metadata-bloom-filters-and-page-indexes).

## Write Required Columns

Required columns are straightforward: one value per logical row, no level arrays.

```c
int64_t ids[] = {1, 2, 3};
carquet_byte_array_t names[] = {
    {(uint8_t*)"alice", 5},
    {(uint8_t*)"bob", 3},
    {(uint8_t*)"carol", 5},
};

carquet_writer_write_batch(writer, 0, ids, 3, NULL, NULL);
carquet_writer_write_batch(writer, 1, names, 3, NULL, NULL);
```

## Write Nullable Columns Correctly

This is the rule most users get wrong:

- `num_values` is the logical row count.
- `def_levels` has one entry per logical row.
- `values` contains only the present values, packed contiguously.

Example for logical rows `[1.5, NULL, 3.5, NULL, 5.5]`:

```c
double values[] = {1.5, 3.5, 5.5};
int16_t def_levels[] = {1, 0, 1, 0, 1};

carquet_writer_write_batch(writer, 2, values, 5, def_levels, NULL);
```

You can query the required definition level from the schema with `carquet_schema_max_def_level()` when you are generating levels programmatically.

## Write Nested Columns Without Computing Levels

For a single-level `LIST<T>` or `MAP<K,V>` (schemas built with `carquet_schema_add_list()` / `carquet_schema_add_map()`), you do not have to precompute definition/repetition levels. `carquet_writer_write_list_column()` takes the Arrow columnar layout — `int32` offsets plus optional list- and element-level validity bitmaps (Arrow's LSB-first, bit set ⇒ present) — and shreds it internally. It is the write-side inverse of `carquet_row_batch_column_list()`.

```c
/* Lists: [[10, 20], [], NULL, [30]] on the leaf element column (index 1). */
int32_t offsets[]      = {0, 2, 2, 2, 3};      /* num_lists + 1 entries       */
int32_t values[]       = {10, 20, 30};          /* flattened child elements     */
uint8_t list_validity  = 0x0B;                  /* bits 0,1,3 set; row 2 null   */
carquet_writer_write_list_column(writer, 1, 4, offsets, &list_validity,
                                 values, /*value_validity=*/NULL, &err);
```

For a `MAP<K,V>` call it twice with the *same* `offsets` and `list_validity` — once for the key leaf (`value_validity` NULL; keys are REQUIRED) and once for the value leaf. Offsets must start at 0 and be non-decreasing; a null list requires an OPTIONAL container and a null element an OPTIONAL leaf, else the call returns `CARQUET_ERROR_INVALID_ARGUMENT`. Deeper nesting returns `CARQUET_ERROR_NOT_IMPLEMENTED`. See [`nested-data.md`](./nested-data.md).

## Row Groups, Metadata, and Per-Column Overrides

Important writer invariants:

- Every column must advance by the same logical row count.
- Call `carquet_writer_new_row_group()` only when all columns are aligned.
- `carquet_writer_close()` has the same requirement.

Useful APIs before the first write:

- `carquet_writer_add_metadata()`: add footer key-value pairs
- `carquet_writer_set_column_encoding()`
- `carquet_writer_set_column_compression()`
- `carquet_writer_set_column_statistics()`
- `carquet_writer_set_column_bloom_filter()`
- `carquet_writer_set_sorting_columns()`: declare the row-group sort order (written to every row group's metadata; the writer records the declaration only and does not sort or verify the data)

Per-column overrides must be set after writer creation and before writing data.

### Custom Compression Codecs

`carquet_register_codec(codec, impl)` installs a user-supplied `carquet_custom_codec_t` (compress / decompress / compress_bound + opaque `user_data`) against any `carquet_compression_t` slot. A registered codec wins over the built-in for that slot, so this both fills the slot carquet does not ship a built-in for (`LZO`, and `BROTLI` when built with `CARQUET_ENABLE_BROTLI=OFF`) and lets you swap a built-in for an alternative implementation (e.g. a hardware-accelerated `GZIP`). Pass `impl = NULL` to unregister. Registering against `UNCOMPRESSED` is rejected — that path has a no-copy fast lane that must not be intercepted. Registrations are process-wide and not safe to mutate while reader or writer threads are mid-compress / mid-decompress; install codecs at startup before opening files.

### Encoding defaults

By default columns use `PLAIN`, with automatic `BYTE_STREAM_SPLIT` for `FLOAT`/`DOUBLE` columns when a compression codec is set. This favors carquet's fast (near zero-copy) read path; it does not dictionary-encode by default.

Opt into another encoding per column with `carquet_writer_set_column_encoding()`: `RLE_DICTIONARY` (emits a `PLAIN` dictionary page plus `RLE_DICTIONARY` data pages, with automatic fallback to `PLAIN` when the dictionary would exceed `dictionary_page_size` or the data is effectively all-unique; eligible for `INT32`/`INT64`/`FLOAT`/`DOUBLE`/`BYTE_ARRAY`/`FIXED_LEN_BYTE_ARRAY`), `BYTE_STREAM_SPLIT` (FLOAT/DOUBLE/INT32/INT64/FLBA), `DELTA_BINARY_PACKED` (INT32/INT64), `DELTA_LENGTH_BYTE_ARRAY` / `DELTA_BYTE_ARRAY` (BYTE_ARRAY), or `RLE` for `BOOLEAN` (run-length/bit-packed hybrid; efficient for long boolean runs). Dictionary encoding produces smaller files but is slower to read, so it is a deliberate choice rather than the default.

## Column Statistics

With `opts.write_statistics = true` (the default), the writer records per-row-group min/max and null counts for every primitive type — `INT32`, `INT64`, `FLOAT`, `DOUBLE`, `BOOLEAN`, `BYTE_ARRAY`, and `FIXED_LEN_BYTE_ARRAY`. Stats are aggregated across pages and used by `carquet_reader_filter_row_groups()` for predicate pushdown. Dictionary-encoded columns additionally carry an exact `distinct_count` (the dictionary size).

The writer also emits **Parquet 2.9 SizeStatistics** (`ColumnMetaData.size_statistics`) where they carry information: `unencoded_byte_array_data_bytes` for `BYTE_ARRAY` columns, and definition/repetition level histograms for nullable / nested / repeated columns. The page index carries the matching per-page level histograms and the spec's `unencoded_byte_array_data_bytes` (OffsetIndex field 2). Flat required non-`BYTE_ARRAY` columns are unaffected.

`BYTE_ARRAY` min/max are stored using lexicographic order; values longer than 32 bytes are truncated per the Parquet spec (min is truncated to a prefix; max is truncated and incremented so the stored bound is still a valid upper bound), and `is_min/max_value_exact` is set to `false` on truncation.

Type-specific statistics behavior, all automatic:

- `FLOAT16` (`FIXED_LEN_BYTE_ARRAY(2)`) min/max are ordered by the represented floating-point value with NaNs skipped, and a zero bound is normalized to `-0.0` (min) / `+0.0` (max), per the spec — not byte-lexicographically.
- `INT96` has no defined sort order, so no min/max statistics are written for it (matching parquet-cpp).
- `GEOMETRY` / `GEOGRAPHY` columns instead get `GeospatialStatistics` written into the column metadata: a coordinate bounding box (`xmin/xmax/ymin/ymax`, plus `z`/`m` when present) and the set of ISO-WKB geometry type codes, computed by parsing the column's WKB values. Regular min/max remain suppressed for these types.

Use `carquet_writer_set_column_statistics(writer, idx, false)` to disable stats for a single column while keeping them on globally.

```c
carquet_writer_options_t opts;
carquet_writer_options_init(&opts);
opts.write_statistics = true;  /* default — shown for clarity */

carquet_writer_t* writer = carquet_writer_create("out.parquet", schema, &opts, &err);

/* Skip stats for one column (e.g. a large opaque blob). */
carquet_writer_set_column_statistics(writer, 2, false);
```

See [`reading.md`](./reading.md#column-statistics) for how to read these stats back and use them for predicate pushdown.

## Append to an Existing File

`carquet_writer_open_append()` opens an existing Parquet file with read+write access and returns a writer positioned to add new row groups. The closing footer is rewritten to list the existing row groups followed by the new ones. Existing bloom filters, page indexes, and key-value metadata are preserved.

```c
carquet_schema_t* schema = /* ...same shape as the existing file's columns... */;

carquet_writer_t* writer = carquet_writer_open_append(
    "events.parquet", schema, /* options= */ NULL, &err);
if (!writer) {
    fprintf(stderr, "%s\n", err.message);
    carquet_schema_free(schema);
    return 1;
}
carquet_schema_free(schema);

int32_t batch[] = { /* new rows */ };
carquet_writer_write_batch(writer, 0, batch, sizeof(batch) / sizeof(*batch), NULL, NULL);
carquet_writer_close(writer);
```

What the writer validates and preserves:

- The supplied schema must match the existing file's leaf columns by count, name, physical type, and repetition. Logical-type metadata for already-written row groups stays exactly as recorded in the existing footer; only new row groups follow the current `schema` + `options`.
- Existing key-value footer metadata is carried over. Further `carquet_writer_add_metadata()` calls add to it.
- Bloom filters and page indexes for the existing row groups stay valid — they live between the row group data and the old footer, which is the region the writer overwrites.
- A schema mismatch surfaces as `CARQUET_ERROR_INVALID_SCHEMA` and the file is left untouched.

When this is the right tool:

- Streaming ingestion, log accumulation, periodic snapshots — anywhere you would otherwise read-then-rewrite.
- The footer is rewritten on every close, so each append pays the cost of one footer serialization, not a full file rewrite.

Limitations to be aware of:

- The writer needs read+write access to the file; `r+b` open mode is used internally.
- Existing bloom filters and page indexes are preserved as-is; they are not merged with bloom filters / page indexes you write for the new row groups (each row group's auxiliary structures stay scoped to that row group, which is how Parquet readers expect them).

## Write to Memory Instead of a File

Using an existing schema, you can write to memory and retrieve the final Parquet bytes after close:

```c
carquet_writer_t* writer = carquet_writer_create_buffer(schema, NULL, &err);
if (!writer) {
    fprintf(stderr, "%s\n", err.message);
    return 1;
}

int32_t xs[] = {10, 20, 30};
carquet_writer_write_batch(writer, 0, xs, 3, NULL, NULL);
carquet_writer_close(writer);

void* buffer = NULL;
size_t size = 0;
if (carquet_writer_get_buffer(writer, &buffer, &size) == CARQUET_OK) {
    /* buffer now belongs to the caller */
    free(buffer);
}
```

This is useful for RPC payloads, tests, and embedding Parquet output inside a larger container format.

## Import from Arrow (C Data Interface)

To write data that already lives in an Arrow `ArrowArray` (from PyArrow, DuckDB, nanoarrow, any C Data Interface producer) without a manual copy, hand the struct array straight to the writer:

```c
/* writer schema and Arrow array must have matching columns, in order */
carquet_status_t st = carquet_writer_write_arrow(writer, &array, &schema, &err);
/* write_arrow consumes (releases) both `array` and `schema` — do not release again */
```

- `carquet_writer_write_arrow()` runs a generic Dremel record-shredder over the Arrow array tree, converting validity bitmaps to Parquet definition/repetition levels and compacting values. It handles **arbitrary-depth nesting** — `struct` (`+s`), `list` (`+l`), `large_list` (`+L`) and `map` (`+m`) composed to any depth — mapping each leaf to the corresponding writer column. It follows Arrow "move" semantics: the passed `array` and `schema` are consumed (their `release` callbacks are invoked) on both success and failure. See [`nested-data.md`](./nested-data.md#arbitrary-depth-nesting-via-the-arrow-c-data-interface).
- `carquet_arrow_import_schema()` turns an Arrow `"+s"` struct schema into a `carquet_schema_t*` you can build a writer from — mapping each child's Arrow `format` string back to a physical + logical type.
- Nested arrays are supported (see `carquet_writer_write_arrow()` above). Sliced arrays (`offset != 0`, at any level, including a struct offset that applies to its children) are written as the slice. Dictionary-encoded arrays (pandas categoricals, `pyarrow.DictionaryArray` with primitive values) are resolved to their values; a null dictionary slot yields a null row. Fixed-size lists (`+w:n`) are written as `LIST`. Unions and run-end-encoded arrays are rejected with a clear error.
- Arrow types without a direct Parquet physical equivalent are coerced the way pyarrow does: `decimal32/64/128/256` (`d:p,s[,bits]`) become `FIXED_LEN_BYTE_ARRAY` `DECIMAL(p, s)` with parquet-cpp's byte length for the precision and big-endian two's-complement bytes; `date64` (`tdm`) becomes `DATE` in days; `time32[s]` (`tts`) becomes `TIME(MILLIS)`; `timestamp[s]` (`tss:tz`) becomes `TIMESTAMP(MILLIS)`; `int8/16` and `uint8/16` become `INT32` with an `INTEGER(8|16, signed)` annotation.

## Failure Path

If any write step fails, call `carquet_writer_abort()` unless `carquet_writer_close()` has already succeeded. `abort()` frees writer resources but does not produce a valid Parquet file.
