# Recipes

Complete, compilable patterns for the integrations users ask for most. Every snippet assumes `#include <carquet/carquet.h>` and `carquet_error_t err = CARQUET_ERROR_INIT;`. Error handling is shown once per recipe and elided afterwards; in real code check every call. Runnable versions live in `examples/`.

## 1. Flat table with required, nullable and string columns

```c
carquet_schema_t* schema = carquet_schema_create(&err);
carquet_logical_type_t str = { .id = CARQUET_LOGICAL_STRING };
carquet_schema_add_column(schema, "id",    CARQUET_PHYSICAL_INT64,      NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
carquet_schema_add_column(schema, "name",  CARQUET_PHYSICAL_BYTE_ARRAY, &str, CARQUET_REPETITION_REQUIRED, 0, 0);
carquet_schema_add_column(schema, "score", CARQUET_PHYSICAL_DOUBLE,     NULL, CARQUET_REPETITION_OPTIONAL, 0, 0);

carquet_writer_options_t o; carquet_writer_options_init(&o);
o.compression = CARQUET_COMPRESSION_ZSTD;      /* level 1 by default: best size/speed trade-off */
o.write_page_index = true;                     /* lets readers prune pages and use page filters */
carquet_writer_t* w = carquet_writer_create("people.parquet", schema, &o, &err);
if (!w) { fprintf(stderr, "%s\n", err.message); carquet_schema_free(schema); return 1; }
carquet_schema_free(schema);                   /* writer keeps its own copy */

int64_t ids[] = {1, 2, 3, 4, 5};
carquet_byte_array_t names[] = { {(uint8_t*)"ann", 3}, {(uint8_t*)"bob", 3}, {(uint8_t*)"", 0},
                                 {(uint8_t*)"dee", 3}, {(uint8_t*)"eve", 3} };
/* logical scores: 1.5, NULL, 3.5, NULL, 5.5  ->  dense values + one def level per row */
double scores[]    = {1.5, 3.5, 5.5};
int16_t score_def[] = {1, 0, 1, 0, 1};        /* 1 == carquet_schema_max_def_level(schema, 2) */

carquet_writer_write_batch(w, 0, ids,    5, NULL, NULL);
carquet_writer_write_batch(w, 1, names,  5, NULL, NULL);
carquet_writer_write_batch(w, 2, scores, 5, score_def, NULL);
if (carquet_writer_close(w) != CARQUET_OK) { carquet_writer_abort(w); return 1; }
```

Streaming: call `write_batch` repeatedly per column in any interleaving; every column must have received the same number of logical rows before `carquet_writer_new_row_group()` (optional, start a new row group) and before `close()`. Row groups also flush automatically at `row_group_size` bytes.

## 2. Timestamps, dates, decimals, UUID, fixed-size binary

```c
carquet_logical_type_t ts  = { .id = CARQUET_LOGICAL_TIMESTAMP,
                               .params.timestamp = { .unit = CARQUET_TIME_UNIT_MICROS, .is_adjusted_to_utc = true } };
carquet_logical_type_t day = { .id = CARQUET_LOGICAL_DATE };                        /* INT32 days since epoch */
carquet_logical_type_t dec = { .id = CARQUET_LOGICAL_DECIMAL, .params.decimal = { .precision = 18, .scale = 2 } };
carquet_logical_type_t u8  = { .id = CARQUET_LOGICAL_INTEGER, .params.integer = { .bit_width = 8, .is_signed = false } };
carquet_logical_type_t uuid = { .id = CARQUET_LOGICAL_UUID };

carquet_schema_add_column(schema, "ts",    CARQUET_PHYSICAL_INT64, &ts,  CARQUET_REPETITION_REQUIRED, 0, 0);   /* int64 micros */
carquet_schema_add_column(schema, "day",   CARQUET_PHYSICAL_INT32, &day, CARQUET_REPETITION_REQUIRED, 0, 0);   /* int32 days   */
carquet_schema_add_column(schema, "price", CARQUET_PHYSICAL_INT64, &dec, CARQUET_REPETITION_REQUIRED, 0, 0);   /* int64 unscaled: 12.34 -> 1234 */
carquet_schema_add_column(schema, "flags", CARQUET_PHYSICAL_INT32, &u8,  CARQUET_REPETITION_REQUIRED, 0, 0);   /* still int32 in memory */
carquet_schema_add_column(schema, "uuid",  CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY, &uuid, CARQUET_REPETITION_REQUIRED, 16, 0); /* uint8_t[16] packed */
```

Rules: DECIMAL on INT32 allows precision ≤ 9, on INT64 ≤ 18, wider precisions use `FIXED_LEN_BYTE_ARRAY` with big-endian two's-complement bytes (length N holds floor((8N−1)·log10 2) digits; 16 bytes = 38 digits) or `BYTE_ARRAY`. TIME(MILLIS) is INT32, TIME(MICROS/NANOS) and every TIMESTAMP are INT64. INT96 timestamps are readable and writable but deprecated: prefer TIMESTAMP. The in-memory type never changes with the logical annotation (INTEGER(8) is still an `int32_t` buffer).

## 3. Lists and maps without computing levels

```c
/* schema: id INT32, tags LIST<STRING> (optional list, required elements) */
carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
carquet_schema_add_list(schema, "tags", CARQUET_PHYSICAL_BYTE_ARRAY, &str, CARQUET_REPETITION_OPTIONAL, 0, 0);
/* leaf columns are now: 0 = id, 1 = tags.list.element */

/* rows: ["a","b"], [], NULL, ["c"]  in Arrow layout */
int32_t offsets[] = {0, 2, 2, 2, 3};                        /* num_lists + 1 */
carquet_byte_array_t elems[] = { {(uint8_t*)"a",1}, {(uint8_t*)"b",1}, {(uint8_t*)"c",1} };
uint8_t list_valid = 0x0B;                                  /* bits 0,1,3 set: row 2 is a null list */
carquet_writer_write_list_column(w, 1, 4, offsets, &list_valid, elems, /*value_validity=*/NULL, &err);
```

A `MAP<K,V>` built with `carquet_schema_add_map()` has two leaves (key, value): call `write_list_column` twice with the same `offsets` and `list_validity`, keys with `value_validity = NULL` (keys are REQUIRED). Reading back through the batch reader gives the same layout via `carquet_row_batch_column_list()` (recipe 7). Deeper nesting: use the Arrow bridge (recipe 9) or write levels by hand (`docs/nested-data.md`).

## 4. Scan with projection, batch by batch

```c
carquet_reader_options_t ro; carquet_reader_options_init(&ro);
ro.use_mmap = true;                                   /* local files: biggest read-side win */
carquet_reader_t* r = carquet_reader_open("people.parquet", &ro, &err);

carquet_batch_reader_config_t cfg; carquet_batch_reader_config_init(&cfg);
const char* cols[] = {"id", "score"};
cfg.column_names = cols; cfg.num_column_names = 2;    /* or column_indices + num_columns */
cfg.batch_size = 65536;
carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);

carquet_row_batch_t* batch = NULL;
carquet_status_t st;
while ((st = carquet_batch_reader_next(br, &batch)) == CARQUET_OK && batch) {
    int64_t n = carquet_row_batch_num_rows(batch);
    const void* d; const uint8_t* valid; int64_t cnt;
    carquet_row_batch_column(batch, 0, &d, &valid, &cnt);           /* index within the projection */
    const int64_t* id = (const int64_t*)d;
    carquet_row_batch_column(batch, 1, &d, &valid, &cnt);
    const double* score = (const double*)d;      /* ROW-ALIGNED: cnt == n, null slots hold 0 */
    for (int64_t i = 0; i < n; i++) {
        bool present = !valid || (valid[i >> 3] & (1u << (i & 7)));
        if (present) { /* use id[i], score[i] */ }
    }
    carquet_row_batch_free(batch); batch = NULL;
}
if (st != CARQUET_OK && st != CARQUET_ERROR_END_OF_DATA) {
    const carquet_error_t* e = carquet_batch_reader_last_error(br);
    fprintf(stderr, "%s\n", e ? e->message : carquet_status_string(st));
}
carquet_batch_reader_free(br);
carquet_reader_close(r);
```

The batch reader is row-aligned: for an OPTIONAL column `cnt` equals the batch row count, null rows are zero-filled slots (`{NULL, 0}` for BYTE_ARRAY) and `valid` bit `i` describes slot `i`. Only the column reader (recipe 6) and the writer use the dense layout. BYTE_ARRAY columns return `const carquet_byte_array_t*`; FIXED_LEN_BYTE_ARRAY returns packed bytes, stride = `carquet_schema_node_type_length()`.

## 5. Predicate pushdown and page filters

```c
/* 1) cheap: which row groups can contain id > 5000? */
int64_t threshold = 5000;
int32_t rgs[64];
int32_t n = carquet_reader_filter_row_groups(r, /*column=*/0, CARQUET_COMPARE_GT, &threshold, sizeof threshold, rgs, 64);

/* 2) page-level: only decompress pages whose min/max can match (file needs write_page_index) */
int64_t lo = 1000, hi = 2000;
carquet_filter_clause_t clause = { .column_index = 0, .op = CARQUET_FILTER_RANGE,
                                   .has_lo = true, .lo = &lo, .lo_size = sizeof lo,
                                   .has_hi = true, .hi = &hi, .hi_size = sizeof hi };
carquet_batch_reader_set_page_filter(br, &clause, 1);   /* clauses are AND-ed; NULL,0 clears; pointers must outlive use */
/* ... read as in recipe 4; rows inside a kept page may still fail the predicate: re-check per row ... */
int64_t skipped = carquet_batch_reader_rows_skipped(br);
```

Equality on strings: `.op = CARQUET_FILTER_EQ, .value = "ann", .value_size = 3` (raw bytes, no `carquet_byte_array_t`). Row groups are also pruned automatically by statistics and bloom filters as soon as a filter is installed. A row window independent of predicates: `carquet_batch_reader_set_row_range(br, offset, limit)`.

## 6. Column reader with levels (one row group, one column)

```c
carquet_column_reader_t* c = carquet_reader_get_column(r, /*row_group=*/0, /*leaf=*/1, &err);
int16_t max_def = carquet_schema_max_def_level(carquet_reader_schema(r), 1);
int32_t  vals[1024]; int16_t def[1024]; int16_t rep[1024];
for (;;) {
    carquet_error_t rerr = CARQUET_ERROR_INIT;
    int64_t n = carquet_column_read_batch_ex(c, vals, 1024, def, rep, &rerr);
    if (n < 0) { /* hard failure: rerr.message */ break; }
    if (carquet_error_is_set(&rerr)) { /* n values valid, batch truncated by an error */ }
    if (n == 0) break;
    int64_t lists = carquet_count_rows(rep, n);           /* top-level rows in this batch */
    int64_t offs[1025];
    carquet_list_offsets(rep, n, /*list_rep_level=*/1, offs, 1025);
    /* def[i] == max_def -> value present; smaller -> null at some level */
}
carquet_column_reader_free(c);
```

`n` counts level entries (list slots), not rows and not present values; unlike the batch reader, the column reader's values array is dense (present values only, in order). `carquet_column_skip()` advances without decoding.

## 7. Lists through the batch reader (single level)

```c
const int32_t* offsets; int64_t num_lists; const void* vals; const uint8_t* elem_valid; int64_t num_vals; const uint8_t* list_valid;
carquet_row_batch_column_list(batch, /*projected col*/ 1, &offsets, &num_lists, &vals, &elem_valid, &num_vals, &list_valid);
for (int64_t i = 0; i < num_lists; i++) {
    bool null_list = list_valid && !(list_valid[i >> 3] & (1u << (i & 7)));
    for (int32_t k = offsets[i]; k < offsets[i + 1]; k++) { /* element k of ((const T*)vals) */ }
}
```

When any projected column is repeated the batch reader returns one row group per batch, and `carquet_row_batch_column()` rejects that column. `max_rep_level > 1` (list of list) is not materialised here: use recipe 9.

## 8. Append, buffer output, metadata

```c
/* append row groups to an existing file (schema must match leaf names/types/repetitions) */
carquet_writer_t* a = carquet_writer_open_append("events.parquet", schema, NULL, &err);
carquet_writer_write_batch(a, 0, more_ids, n, NULL, NULL);
carquet_writer_close(a);

/* write to memory */
carquet_writer_t* m = carquet_writer_create_buffer(schema, NULL, &err);
/* ... write_batch ... */
carquet_writer_close(m);
void* bytes; size_t size;
carquet_writer_get_buffer(m, &bytes, &size);          /* only after close(); caller free()s */
carquet_reader_t* rb = carquet_reader_open_buffer(bytes, size, NULL, &err);   /* bytes must outlive rb */

/* footer key/value metadata and per-field labels */
carquet_writer_add_metadata(w, "producer", "myapp 1.2");        /* before the first write */
carquet_schema_set_field_metadata(schema, carquet_schema_num_elements(schema) - 1, "Label", "Age in years");
o.write_arrow_schema = true;                                    /* needed for per-field metadata to be emitted */
const char* v = carquet_reader_find_metadata(r, "producer");
const char* label = carquet_reader_column_find_metadata(r, /*leaf*/ 3, "Label");
```

## 9. Arrow C Data Interface (in-process exchange with Arrow, DuckDB, PyArrow, nanoarrow)

```c
/* write anything nested from an Arrow struct array; both structs are CONSUMED (release called) */
struct ArrowArray  arr;  struct ArrowSchema sch;   /* filled by any Arrow producer */
carquet_schema_t* cs = NULL;
carquet_arrow_import_schema(&sch_copy, &cs, &err);             /* consumes sch_copy; build the writer from cs */
carquet_writer_t* w = carquet_writer_create("nested.parquet", cs, NULL, &err);
carquet_writer_write_arrow(w, &arr, &sch, &err);               /* consumes arr and sch, success or failure */
carquet_writer_close(w);

/* read one row group back as a nested Arrow tree (owned copies, release exactly once) */
struct ArrowSchema out_s; struct ArrowArray out_a;
if (carquet_reader_read_arrow(r, 0, &out_s, &out_a, &err) == CARQUET_OK) {
    /* out_a.children[i] mirror the file's top-level fields */
    out_a.release(&out_a); out_s.release(&out_s);
}

/* flat fast path from a row batch (no projection): */
carquet_arrow_export_batch(batch, carquet_reader_schema(r), &out_s, &out_a, &err);
```

Import handles sliced arrays, dictionary arrays (resolved to values), fixed-size lists, decimals, date64/second timestamps; unions and run-end encoding are rejected. The `ArrowSchema` / `ArrowArray` struct definitions are in `carquet.h` under the standard `ARROW_C_DATA_INTERFACE` guard, so including Arrow's `abi.h` or nanoarrow alongside is safe.

## 10. Dictionary-preserving string scans

```c
cfg.preserve_dictionaries = true;
/* ... next() ... */
const uint32_t* idx; const uint8_t* valid; int64_t n; const uint8_t* dict; int32_t dict_count; const uint32_t* dict_offsets;
if (carquet_row_batch_column_dictionary(batch, 0, &idx, &valid, &n, &dict, &dict_count, &dict_offsets) == CARQUET_OK) {
    /* row i (if present) is dictionary entry idx[i]; BYTE_ARRAY entry k starts at dict + dict_offsets[k]:
       4-byte little-endian length, then the bytes. Fixed-width types: dict_offsets is NULL, entry k at dict + k * width. */
} else {
    /* column was not dictionary-encoded in this row group: fall back to carquet_row_batch_column() */
}
```

Skips materialising every string; large win for low-cardinality text. Only columns written with `RLE_DICTIONARY` (opt in per column on the write side) benefit.

## 11. Per-column tuning before the first write

```c
carquet_writer_set_column_encoding(w, 1, CARQUET_ENCODING_RLE_DICTIONARY);   /* low-cardinality strings */
carquet_writer_set_column_encoding(w, 2, CARQUET_ENCODING_DELTA_BINARY_PACKED); /* sorted ints */
carquet_writer_set_column_compression(w, 3, CARQUET_COMPRESSION_UNCOMPRESSED, 0); /* already-compressed blobs */
carquet_writer_set_column_statistics(w, 3, false);
carquet_writer_set_column_bloom_filter_options(w, 0, true, /*ndv=*/1000000, /*fpp=*/0.01);
carquet_sorting_column_t sc = { .column_index = 0, .descending = false, .nulls_first = false };
carquet_writer_set_sorting_columns(w, &sc, 1);                               /* metadata only, data is not sorted */
```

Defaults are PLAIN (BYTE_STREAM_SPLIT for FLOAT/DOUBLE under a codec), no dictionary: fastest reads. Dictionary encoding is a size trade-off you opt into.

## 12. Custom or hardware codec

```c
static carquet_status_t my_compress(const uint8_t* s, size_t n, uint8_t* d, size_t cap, size_t* out, int32_t level, void* ud);
static carquet_status_t my_decompress(const uint8_t* s, size_t n, uint8_t* d, size_t cap, size_t* out, void* ud);
static size_t my_bound(size_t n, void* ud);
carquet_custom_codec_t impl = { my_compress, my_decompress, my_bound, /*user_data=*/NULL };
carquet_register_codec(CARQUET_COMPRESSION_BROTLI, &impl);   /* process-wide; install at startup; NULL unregisters */
```

Fills the BROTLI / LZO slots carquet has no built-in for, or replaces a built-in (e.g. an accelerated GZIP). `UNCOMPRESSED` cannot be intercepted.

## 13. Robust error reporting

```c
char buf[512];
if (!reader) { carquet_error_format(&err, buf, sizeof buf); log_error("%s", buf); }   /* code, message, offset, column, row group */
const char* hint = carquet_error_recovery_hint(err.code);                             /* may be NULL */
```

Codes: `CARQUET_ERROR_INVALID_ARGUMENT` (shape / type / level mistakes in your call), `TYPE_MISMATCH`, `COLUMN_NOT_FOUND`, `END_OF_DATA` (normal terminator), `PAGE_INDEX_REQUIRED` (page filter on a file without offset index), `INVALID_SCHEMA` (append with a different schema), `DECODE` / `INVALID_COMPRESSED_DATA` / `CRC_MISMATCH` (corrupt input), `NOT_IMPLEMENTED` (unsupported nesting shape or Arrow type), `OUT_OF_MEMORY`.
