# Error Handling and Type Reference

## Error Model

Every fallible function returns `carquet_status_t` (an enum) or uses an `carquet_error_t*` out-parameter. Programming errors (NULL where nonnull is required) trigger assertions in debug builds; runtime errors return codes.

```c
carquet_error_t err = CARQUET_ERROR_INIT;
carquet_reader_t* r = carquet_reader_open("data.parquet", NULL, &err);
if (!r) {
    fprintf(stderr, "[%s] %s\n",
            carquet_status_string(err.code), err.message);
    return 1;
}
```

## Status Codes

| Range | Category | Examples |
|---|---|---|
| 0 | Success | `CARQUET_OK` |
| 1-4 | General | `INVALID_ARGUMENT`, `OUT_OF_MEMORY`, `NOT_IMPLEMENTED`, `INTERNAL` |
| 10-15 | File I/O | `FILE_NOT_FOUND`, `FILE_OPEN`, `FILE_READ`, `FILE_WRITE`, `FILE_TRUNCATED` |
| 20-26 | Format | `INVALID_MAGIC`, `INVALID_FOOTER`, `INVALID_SCHEMA`, `INVALID_PAGE` |
| 30-33 | Thrift | `THRIFT_DECODE`, `THRIFT_ENCODE`, `THRIFT_TRUNCATED` |
| 40-44 | Encoding | `DECODE`, `ENCODE`, `DICTIONARY_NOT_FOUND`, `INVALID_RLE` |
| 50-53 | Compression | `COMPRESSION`, `DECOMPRESSION`, `UNSUPPORTED_CODEC` |
| 60-63 | Data | `TYPE_MISMATCH`, `COLUMN_NOT_FOUND`, `ROW_GROUP_NOT_FOUND`, `END_OF_DATA` |
| 70-71 | Checksum | `CHECKSUM`, `CRC_MISMATCH` |
| 80-82 | State | `INVALID_STATE`, `ALREADY_CLOSED`, `NOT_OPEN` |

## Rich Error Context

`carquet_error_t` carries more than a code and message:

```c
typedef struct carquet_error {
    carquet_status_t code;
    char message[256];
    const char* file;       /* source file where error was set */
    int line;
    const char* function;
    int64_t offset;         /* file offset, -1 if not applicable */
    int32_t column_index;   /* -1 if not applicable */
    int32_t row_group_index;
} carquet_error_t;
```

Format everything into a single string:

```c
char buf[512];
carquet_error_format(&err, buf, sizeof(buf));
fprintf(stderr, "%s\n", buf);
```

Get a recovery hint (returns NULL if none available):

```c
const char* hint = carquet_error_recovery_hint(err.code);
if (hint) fprintf(stderr, "Hint: %s\n", hint);
```

Check recoverability (e.g., transient I/O vs permanent corruption):

```c
if (carquet_error_is_recoverable(err.code)) {
    /* retry may succeed */
}
```

## Macros for Library Internals

These are useful if you extend carquet or write custom encodings:

```c
CARQUET_RETURN_IF_ERROR(status);          /* early return on failure */
CARQUET_SET_ERROR(err, code, fmt, ...);   /* set error with source location */
CARQUET_CHECK(cond, err, code, fmt, ...); /* assert + set error + return */
CARQUET_SUCCEEDED(status);                /* status == CARQUET_OK */
CARQUET_FAILED(status);                   /* status != CARQUET_OK */
```

## Per-Read Error Reporting on the Column Reader

`carquet_column_read_batch()` returns the number of values read, or `-1` on any error — it cannot say *which* error, and it reports a mid-batch read failure that already produced some values by returning that partial count (indistinguishable from a clean end-of-column short read). `carquet_column_read_batch_ex()` takes the same arguments plus a trailing `carquet_error_t*` and follows the out-parameter convention used everywhere else in the API:

| Return | `error` | Meaning |
|---|---|---|
| `>= 0` | unset (`CARQUET_OK`) | Clean read; a value below the requested count means end-of-column. |
| `> 0` | set | Partial read — values are valid but a read error truncated the batch; the rest was not read. |
| `-1` | set | Hard failure, nothing read (`INVALID_ARGUMENT`, `TYPE_MISMATCH`, `OUT_OF_MEMORY`, or the propagated page/decode/I-O code). |

The error is cleared on entry, so `error->code == CARQUET_OK` is a reliable "no failure" signal regardless of the returned count. Passing `NULL` for `error` reproduces `carquet_column_read_batch()` exactly.

## Per-Read Error Reporting on the Batch Reader

`carquet_batch_reader_next()` has no error out-parameter — it can only return a status code, and on a decode failure that status is the coarse `CARQUET_ERROR_DECODE`. `carquet_batch_reader_last_error()` recovers the detail that was built deeper in the reader:

```c
carquet_row_batch_t* batch = NULL;
carquet_status_t st;
while ((st = carquet_batch_reader_next(br, &batch)) == CARQUET_OK && batch) {
    /* ... */
    carquet_row_batch_free(batch);
    batch = NULL;
}
if (st != CARQUET_OK && st != CARQUET_ERROR_END_OF_DATA) {
    const carquet_error_t* e = carquet_batch_reader_last_error(br);
    fprintf(stderr, "read failed: %s\n", e ? e->message : carquet_status_string(st));
}
```

- Returns `NULL` when the last `carquet_batch_reader_next()` did not fail. `CARQUET_ERROR_END_OF_DATA` is the normal terminator and leaves the error unset.
- When the failure came from a column read, the returned error is that column's verbatim error and its `code` is the *underlying, more specific* status. A corrupt Snappy page makes `next()` return `CARQUET_ERROR_DECODE` while `e->code` is `CARQUET_ERROR_INVALID_COMPRESSED_DATA` and `e->message` names the page's absolute file offset and the exact codec check that fired:

  ```
  Failed to decompress data page at file offset 2214592512: snappy[10]: copy offset
  points before start of output [src=25 dst=20 declared=9242833 tag=30 op=100]
  (codec=1, compressed=1599289, uncompressed=9242833)
  ```

- For failures outside a column read the message falls back to the status name and `code` equals what `next()` returned.
- With several projected columns, the reported error is the first failing column in projection order.
- The pointer is owned by the batch reader and is invalidated by the next `carquet_batch_reader_next()` call or by `carquet_batch_reader_free()`; copy the message out if you need to keep it.

## Physical Type ↔ C Type Mapping

When calling `carquet_writer_write_batch()` or reading from `carquet_column_read_batch()` / `carquet_column_read_batch_ex()` / `carquet_row_batch_column()`, the `values` buffer must match the column's physical type:

| Physical Type | C Type | Size | Notes |
|---|---|---|---|
| `CARQUET_PHYSICAL_BOOLEAN` | `uint8_t` | 1 | One byte per value (not bit-packed in the API) |
| `CARQUET_PHYSICAL_INT32` | `int32_t` | 4 | |
| `CARQUET_PHYSICAL_INT64` | `int64_t` | 8 | |
| `CARQUET_PHYSICAL_INT96` | `carquet_int96_t` | 12 | Deprecated legacy timestamps; readable and writable (PLAIN only, no statistics) |
| `CARQUET_PHYSICAL_FLOAT` | `float` | 4 | IEEE 754 |
| `CARQUET_PHYSICAL_DOUBLE` | `double` | 8 | IEEE 754 |
| `CARQUET_PHYSICAL_BYTE_ARRAY` | `carquet_byte_array_t` | ptr+len | Variable-length; data pointer valid until batch/reader freed |
| `CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY` | `uint8_t[N]` | N | N = `type_length` from schema; tightly packed |

## Logical Type Quick Reference

Logical types annotate physical types with semantic meaning. Set via `carquet_logical_type_t` when adding columns.

| Logical Type | Physical Type | Params |
|---|---|---|
| `STRING` | `BYTE_ARRAY` | — |
| `DATE` | `INT32` | Days since 1970-01-01 |
| `TIME` | `INT32` or `INT64` | `.time.unit`, `.time.is_adjusted_to_utc` |
| `TIMESTAMP` | `INT64` | `.timestamp.unit`, `.timestamp.is_adjusted_to_utc` |
| `DECIMAL` | `INT32`, `INT64`, `FIXED_LEN_BYTE_ARRAY`, or `BYTE_ARRAY` | `.decimal.precision`, `.decimal.scale` |
| `INTEGER` | `INT32` or `INT64` | `.integer.bit_width`, `.integer.is_signed` |
| `UUID` | `FIXED_LEN_BYTE_ARRAY` (16) | — |
| `JSON` | `BYTE_ARRAY` | — |
| `BSON` | `BYTE_ARRAY` | — |
| `FLOAT16` | `FIXED_LEN_BYTE_ARRAY` (2) | — |
| `LIST` | group | Use `carquet_schema_add_list()` |
| `MAP` | group | Use `carquet_schema_add_map()` |

## Null Bitmap Convention

The batch reader returns a null bitmap where bit `i` is **set** (1) when value `i` is **present** (not null):

```c
bool is_null = null_bitmap && !(null_bitmap[i / 8] & (1u << (i % 8)));
```

If `null_bitmap` is NULL, all values are present (required column).

## Definition / Repetition Level Convention

When writing nullable columns:

- `def_levels` has one entry per **logical row**
- `values` contains only **non-null** values, packed contiguously
- `def_level == max_def_level` means the value is present
- `def_level < max_def_level` means null (at the level indicated)

When writing repeated columns, `rep_levels` additionally encodes list boundaries:

- `rep_level == 0` starts a new top-level record
- `rep_level > 0` continues the current repeated group

Query max levels from the schema rather than hard-coding them:

```c
int16_t max_def = carquet_schema_max_def_level(schema, leaf_index);
int16_t max_rep = carquet_schema_max_rep_level(schema, leaf_index);
```
