# Reading Files

## Choose the Reader API

Use the batch reader unless you specifically need raw Parquet level streams.

| Need | API |
| --- | --- |
| Project a subset of columns, iterate in row batches, parallelize work | `carquet_batch_reader_t` |
| Read one column from one row group, inspect `def_levels` / `rep_levels`, or skip values manually | `carquet_column_reader_t` |
| Only inspect footer metadata | `carquet_get_file_info()`, `carquet_validate_file()` |

## Open a Reader

```c
#include <carquet/carquet.h>

carquet_error_t err = CARQUET_ERROR_INIT;

carquet_reader_options_t opts;
carquet_reader_options_init(&opts);
opts.use_mmap = true;          /* Good default for local files */
opts.verify_checksums = true;  /* Keep enabled unless you trust the source */
opts.num_threads = 0;          /* Auto */

carquet_reader_t* reader = carquet_reader_open("data.parquet", &opts, &err);
if (!reader) {
    char buf[512];
    carquet_error_format(&err, buf, sizeof(buf));
    fprintf(stderr, "%s\n", buf);
    return 1;
}
```

Other entry points:

- `carquet_reader_open_file(FILE*, ...)`: caller keeps ownership of the `FILE*`
- `carquet_reader_open_buffer(const void* buffer, size_t size, ...)`: buffer must outlive the reader

Useful metadata calls right after open:

- `carquet_reader_schema()`
- `carquet_reader_num_rows()`
- `carquet_reader_num_row_groups()`
- `carquet_reader_num_columns()`
- `carquet_reader_is_mmap()`

## Batch Reader Workflow

This is the default path for scans and analytics.

```c
carquet_batch_reader_config_t cfg;
carquet_batch_reader_config_init(&cfg);

const char* cols[] = {"id", "price"};
cfg.column_names = cols;
cfg.num_column_names = 2;
cfg.batch_size = 65536;
cfg.num_threads = 0;

carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &err);
if (!br) {
    fprintf(stderr, "%s\n", err.message);
    carquet_reader_close(reader);
    return 1;
}

carquet_row_batch_t* batch = NULL;
while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
    const void* id_data;
    const uint8_t* id_nulls;
    int64_t id_count;

    carquet_row_batch_column(batch, 0, &id_data, &id_nulls, &id_count);
    const int64_t* ids = id_data;

    carquet_row_batch_free(batch);
    batch = NULL;
}

carquet_batch_reader_free(br);
carquet_reader_close(reader);
```

Notes:

- The loop above exits on any non-`CARQUET_OK` status. To tell end-of-file from a genuine failure, keep the status and consult `carquet_batch_reader_last_error()` — `carquet_batch_reader_next()` has no error out-parameter, so that accessor is the only way to get the message (see [error-handling.md](error-handling.md#per-read-error-reporting-on-the-batch-reader)).
- `column_indices` takes precedence over `column_names`.
- `row_group_filter` lets you skip entire row groups before any data pages are read. Use it with `carquet_reader_row_group_matches()` or `carquet_reader_column_statistics()`.
- `carquet_row_batch_num_columns()` is the number of projected columns, not the file-wide total.
- For nullable columns, the null bitmap uses `1` for present and `0` for null:

```c
bool is_null = null_bitmap && !(null_bitmap[i / 8] & (1u << (i % 8)));
```

- `BYTE_ARRAY` data comes back as `const carquet_byte_array_t*`.
- `FIXED_LEN_BYTE_ARRAY` data comes back as tightly packed bytes. Use the schema's type length to stride the buffer correctly.

## Column Reader Workflow

Use a column reader when you need explicit control over row groups or Parquet levels.

```c
carquet_column_reader_t* col = carquet_reader_get_column(reader, 0, 2, &err);
if (!col) {
    fprintf(stderr, "%s\n", err.message);
    carquet_reader_close(reader);
    return 1;
}

int32_t values[1024];
int16_t def_levels[1024];
int64_t n;

while ((n = carquet_column_read_batch(col, values, 1024, def_levels, NULL)) > 0) {
    /* values contains only materialized values for this physical type */
}

carquet_column_reader_free(col);
```

`carquet_column_read_batch()` returns `-1` for any failure and cannot distinguish a page-read error that truncates a batch after some values were already read from a clean short read at end-of-column — the `n > 0` loop above treats both as "keep going / done". When you need to tell those apart (e.g. so incomplete data is never processed as complete), use `carquet_column_read_batch_ex()`, which reports a distinct status code and message through a `carquet_error_t` out-parameter:

```c
carquet_error_t err = CARQUET_ERROR_INIT;
int64_t n = carquet_column_read_batch_ex(col, values, 1024, def_levels, NULL, &err);
if (n < 0) {
    /* hard failure, nothing read — err.code / err.message say why */
} else if (carquet_error_is_set(&err)) {
    /* n values are valid, but a read error truncated the batch:
       do NOT treat n as a clean end-of-column result */
} else {
    /* clean read of n values (n < requested ⇒ end of column) */
}
```

The error out-parameter is cleared on entry, so `err.code == CARQUET_OK` reliably means "no failure", independent of the count. Passing `NULL` for the error makes it behave exactly like `carquet_column_read_batch()`.

Also available:

- `carquet_column_skip()`
- `carquet_column_has_next()`
- `carquet_column_remaining()`

## Reading Nested Data (lists, structs, maps)

Carquet reads nested schemas — structs (groups), `LIST`, and `MAP` — at the **physical leaf-column level**, which is how Parquet stores nested data on disk. It does not currently reassemble leaves into materialized nested value objects for you (that is, there is no `List<Struct<…>>`-shaped result type); instead you read each leaf column plus its repetition and definition levels and reconstruct the structure yourself. This is a deliberate boundary: the format-level support is complete (nested files written by carquet are read back as native nested types by Arrow/PyArrow and DuckDB), but high-level nested materialization is tracked as future work in `TODO.md`.

What this means in practice:

- A nested field expands to one column chunk **per leaf**. For example a PyArrow `struct<a, b>` is two leaf columns (`a`, `b`); a `map<string,int32>` is two leaves (`key`, `value`); a `list<int32>` is one leaf (`element`). `carquet_reader_num_columns()` therefore counts leaves, not top-level fields.
- Read each leaf with the column-reader workflow above, requesting **definition and repetition levels**. Definition levels tell you which entries are null / which optional ancestors are present; repetition levels tell you where each list/row begins.
- Use the reconstruction helpers to recover structure from the levels:
  - `carquet_count_rows(rep_levels, n, ...)` — count logical (top-level) rows from repetition levels.
  - `carquet_list_offsets(rep_levels, n, list_level, offsets, max)` — compute list boundaries (offsets) for a given repetition level.

```c
/* list<int32> "tags": one leaf, max_def = 3, max_rep = 1 */
carquet_column_reader_t* col = carquet_reader_get_column(reader, 0, leaf_index, &err);

int32_t values[1024];
int16_t def_levels[1024];
int16_t rep_levels[1024];
int64_t n = carquet_column_read_batch(col, values, 1024, def_levels, rep_levels);

int64_t offsets[256];
int64_t num_lists = carquet_list_offsets(rep_levels, n, 1, offsets, 256);
/* offsets[i]..offsets[i+1] delimits the i-th list; def_level < max_def marks
 * null lists / null elements per the Parquet level rules. */

carquet_column_reader_free(col);
```

If you need values presented as native nested types rather than leaves + levels, that reassembly currently lives in the calling application (or a higher-level binding). See `tests/test_nested.c` for end-to-end reconstruction examples.

### Single-Level Lists Through the Batch Reader

For **single-level** repeated columns (`max_rep_level == 1` — i.e. a `LIST<T>` element or a `MAP` key/value leaf), the batch reader reconstructs Arrow's `List<T>` layout for you. When any projected column is repeated, `carquet_batch_reader_next()` reads a whole row group per batch and each repeated column is read via `carquet_row_batch_column_list()` instead of `carquet_row_batch_column()` (which rejects list columns):

```c
const int32_t* offsets;      /* num_lists + 1 entries */
int64_t num_lists, num_values;
const void* values;          /* flattened child (element) array */
const uint8_t* value_validity;  /* child validity bitmap, or NULL */
const uint8_t* list_validity;   /* list-level validity bitmap, or NULL */

carquet_row_batch_column_list(batch, col, &offsets, &num_lists,
                              &values, &value_validity, &num_values, &list_validity);

/* Row i's elements are values[offsets[i] .. offsets[i+1]-1].
 * Element k is null if value_validity && bit k is clear;
 * the whole list (row i) is null if list_validity && bit i is clear. */
```

This is the standard Arrow list encoding (offsets + child array + two validity bitmaps), so it maps directly onto Arrow / a downstream consumer. Deeper nesting (`max_rep_level > 1`, e.g. list-of-list) is not yet materialized this way — `carquet_batch_reader_next()` returns `CARQUET_ERROR_NOT_IMPLEMENTED` — and still needs the manual leaf + levels approach above. Page filters are not combined with repeated-column projections in this release.

## Predicate Pushdown and Cheap Inspection

Use row-group statistics before you start reading payload pages:

- `carquet_reader_column_statistics()`: get min/max/null counts for one row group + column
- `carquet_reader_row_group_matches()`: ask whether one row group might satisfy a predicate
- `carquet_reader_filter_row_groups()`: collect all candidate row groups at once

Use footer-only helpers when you do not need to build a full reader:

- `carquet_get_file_info()`
- `carquet_validate_file()`

### Column Statistics

`carquet_reader_column_statistics()` returns aggregated stats for a single (row group, column). The `min_value` / `max_value` fields are raw bytes — interpret them according to the column's physical type. For `FLOAT16` columns the two bytes are the little-endian IEEE half representation of the numeric min/max (not a lexicographic bound). `INT96` and `GEOMETRY`/`GEOGRAPHY` columns report no min/max (`has_min_max == false`). Legacy files whose V1 pages use the deprecated `BIT_PACKED` level encoding are decoded transparently.

For `GEOMETRY`/`GEOGRAPHY` columns, call `carquet_reader_geospatial_statistics()` to get the coordinate bounding box (`xmin/xmax/ymin/ymax`, plus `z`/`m` when present) and the set of ISO-WKB geometry type codes. It returns `CARQUET_ERROR_INVALID_METADATA` for columns that carry no geospatial statistics (not an error). The `carquet stat` CLI prints this automatically.

```c
carquet_column_statistics_t s;
if (carquet_reader_column_statistics(reader, /*row_group=*/0, /*column=*/0, &s)
        == CARQUET_OK && s.has_min_max) {
    /* Fixed-width numeric types: cast directly. */
    int64_t min_id = *(const int64_t*)s.min_value;
    int64_t max_id = *(const int64_t*)s.max_value;
    printf("id [%lld, %lld] nulls=%lld\n",
           (long long)min_id, (long long)max_id, (long long)s.null_count);
}
```

For `BYTE_ARRAY` columns the payload is the bytes themselves with length in `min_value_size` / `max_value_size` — **not** a `carquet_byte_array_t` struct:

```c
if (s.has_min_max) {
    printf("min=%.*s max=%.*s\n",
           s.min_value_size, (const char*)s.min_value,
           s.max_value_size, (const char*)s.max_value);
}
```

Long byte-array max values written by carquet are truncated at 32 bytes and incremented, so the stored bound is an upper bound but not necessarily an exact value present in the column.

To prune row groups in bulk, pass a typed value to `carquet_reader_filter_row_groups()`:

```c
int32_t threshold = 5000;
int32_t matches[64];
int32_t n = carquet_reader_filter_row_groups(
    reader, /*column=*/0, CARQUET_COMPARE_GT,
    &threshold, sizeof(threshold), matches, 64);
```

## Filtering Rows With a Page Filter

The batch reader can attach a page-level filter that skips data pages whose column-index statistics prove no value can satisfy the predicate. Only the matching pages are decompressed and decoded; pruned pages are never read from disk past their offset-index entry. This applies to both the sequential reader and the parallel pipeline.

```c
carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &err);

int32_t lo = 1000;
int32_t hi = 2000;
carquet_filter_clause_t clause = {
    .column_index = 3,                 /* file column to filter on */
    .op = CARQUET_FILTER_RANGE,
    .has_lo = true, .lo = &lo, .lo_size = sizeof(lo),
    .has_hi = true, .hi = &hi, .hi_size = sizeof(hi),
};
carquet_batch_reader_set_page_filter(br, &clause, 1);
```

Each clause is a `(column_index, op, value)` triple; multiple clauses are AND'd together and may reference different columns (whether or not those columns are part of the projection — a predicate column that is not projected is read **only** via its column + offset index, never decompressed). Supported ops cover `EQ`, `NE`, `LT`, `LE`, `GT`, `GE`, `RANGE`, `IN`, `IS_NULL`, and `IS_NOT_NULL`.

```c
int64_t in_values[] = {10, 20, 30};
carquet_filter_clause_t in_clause = {
    .column_index = 0,
    .op = CARQUET_FILTER_IN,
    .values = in_values,
    .value_count = 3,
};
```

Pass `clauses = NULL, count = 0` to clear an active filter.

After reading, `carquet_batch_reader_rows_skipped()` reports the running count of rows the filter pruned — useful for confirming selectivity:

```c
int64_t skipped = carquet_batch_reader_rows_skipped(br);
```

**Automatic row-group pruning (statistics + bloom filters).** Before any page index is consulted, an installed filter automatically drops whole row groups whose ColumnChunk statistics (min/max/null-count) or — for `EQ`/`IN` clauses — bloom filter prove they cannot match. This needs no user callback and no page index: a predicate that falls outside every row group's range skips the file cleanly (rather than erroring), and a bloom filter drops a row group whose value sits inside the min/max range but is provably absent. It is purely additive to the page-level filtering below; it only ever skips provably-empty row groups. (This supersedes the older manual pattern of calling `carquet_reader_filter_row_groups()` and translating the result into a `row_group_filter` callback.)

**Requirements and semantics**

- Page-*level* filtering (narrowing to row ranges *within* a surviving row group) still requires `write_page_index = true` for every column the filter references. When a row group survives the automatic pruning above but a referenced column has no page index, the next `carquet_batch_reader_next()` returns `CARQUET_ERROR_PAGE_INDEX_REQUIRED`.
- `INT96` columns have no defined sort order per the Parquet spec and are rejected with `CARQUET_ERROR_INVALID_ARGUMENT` at filter-set time.
- For BYTE_ARRAY columns with truncated min/max stats, the filter is conservative: a page may be kept even when no row in it actually matches. Filtering is page-granular by design — rows inside a matching page that fail the predicate are still returned and the caller must apply the exact predicate to the batch if needed.
- Floats whose predicate value is NaN match nothing under ordered ops (`EQ`, `LT`, …), mirroring Arrow semantics.
- The clauses array and any pointers it references must remain valid until the next `set_page_filter()` call or until the batch reader is freed; no copy is made.

## Metadata, Bloom Filters, and Page Indexes

Carquet exposes the optional metadata structures that many readers hide:

- Key-value footer metadata: `carquet_reader_num_metadata()`, `carquet_reader_get_metadata()`, `carquet_reader_find_metadata()`
- Per-field metadata (variable labels/descriptions): `carquet_reader_column_num_metadata()`, `carquet_reader_column_get_metadata()`, and `carquet_reader_column_find_metadata(reader, column, "Label")`. These recover Arrow `Field.custom_metadata` from the file's `ARROW:schema` blob (written by carquet — see [Per-Field Metadata](writing.md#per-field-metadata-variable-labels--descriptions) — or by PyArrow / Arrow C++). Returns nothing when the file has no such blob.
- Arrow type refinements: `carquet_reader_column_arrow_type_refinement(reader, column)` recovers per-leaf Arrow types the Parquet type system alone cannot express — 64-bit-offset `LargeUtf8` / `LargeBinary` — from the `ARROW:schema` blob, returning a `carquet_arrow_type_refinement_t`. The column's values still read exactly as their Parquet physical type; the refinement is informational (e.g. so a caller re-exporting to Arrow can restore the original large-offset type). Returns `CARQUET_ARROW_REFINE_NONE` when absent.
- Bloom filters: `carquet_reader_get_bloom_filter()` plus `carquet_bloom_filter_check_*()`
- Page indexes: `carquet_reader_get_column_index()`, `carquet_reader_get_offset_index()`
- Column chunk metadata: `carquet_reader_column_chunk_metadata()`

Use these when you need diagnostics, custom pruning, or interoperability checks. For a compact end-to-end example, see [`examples/advanced_features.c`](../examples/advanced_features.c).

## Export to Arrow (C Data Interface)

To hand read data to the Arrow ecosystem (PyArrow, DuckDB, nanoarrow, any C Data Interface consumer) without a bespoke copy, export a row batch to a standard `ArrowArray` + `ArrowSchema`:

```c
struct ArrowSchema aschema;
struct ArrowArray  aarray;
const carquet_schema_t* schema = carquet_reader_schema(reader);

if (carquet_arrow_export_batch(batch, schema, &aschema, &aarray, &err) == CARQUET_OK) {
    /* Hand &aschema / &aarray to any Arrow consumer, then release exactly once: */
    aarray.release(&aarray);
    aschema.release(&aschema);
}
```

- `carquet_arrow_export_batch()` covers **flat** columns plus **single-level `LIST<T>` and `MAP<K,V>`** from a row batch — a list becomes an Arrow `List` child (offsets + flattened element) and a map an Arrow `Map` child (`entries` struct of key/value). STRUCT and deeper-than-single-level nesting return `CARQUET_ERROR_NOT_IMPLEMENTED` here; use `carquet_reader_read_arrow()` (below) for those. The batch must have been read without column projection (leaf count must match the schema).
- Every exported buffer is an **independent copy owned by the struct**, so the export stays valid after the source batch is freed or the batch reader advances. The consumer owns it and must call the `release` callback once.
- `carquet_arrow_export_schema()` exports just the schema (a `"+s"` struct), now **including nested types** — a Parquet LIST becomes Arrow `List<element>`, a MAP becomes `Map<entries: struct<key, value>>`, and a group becomes a `Struct`. BOOLEAN is bit-packed, UTF8/binary get `int32` offset buffers, and validity bitmaps are already Arrow-native (LSB-first, present bit set).
- The struct definitions are declared inline in `carquet.h` under the standard `ARROW_C_DATA_INTERFACE` guard, so including a real Arrow `abi.h`/`nanoarrow.h` alongside is safe.

### Nested data: `carquet_reader_read_arrow()`

For **arbitrary-depth nesting** (`LIST<LIST<T>>`, `LIST<STRUCT<...>>`, `MAP<K, LIST<V>>`, …), read a whole row group straight into a nested `ArrowArray` tree instead of a flat row batch:

```c
struct ArrowSchema aschema;
struct ArrowArray  aarray;
if (carquet_reader_read_arrow(reader, /*row_group=*/0, &aschema, &aarray, &err) == CARQUET_OK) {
    /* aarray is a "+s" struct array; its children mirror the file's top-level
       fields with full struct/list/map nesting. Release once when done. */
    aarray.release(&aarray);
    aschema.release(&aschema);
}
```

This reassembles the nested structure from the columns' repetition/definition levels (the read-side inverse of `carquet_writer_write_arrow()`). See [`nested-data.md`](./nested-data.md#arbitrary-depth-nesting-via-the-arrow-c-data-interface).

## Lifetime and Ownership

- Close column readers with `carquet_column_reader_free()`.
- Free each row batch with `carquet_row_batch_free()` before asking for the next one.
- Pointers returned from batch APIs belong to the batch. Do not keep them after freeing the batch.
- Pointers returned from reader metadata APIs belong to the reader. Do not keep them after closing the reader.
