# MetricStore — a carquet real-world testing mockup

MetricStore is a small **mockup downstream library** that consumes [carquet](../) the way a real application would. It models a realistic domain — time-series telemetry / observability events — and drives a broad slice of the public carquet API end to end. Its purpose is *real-world integration testing*: build against carquet as an external dependency, round-trip data, and assert the results.

It re-implements no Parquet logic of its own. The public API (`include/metricstore/metricstore.h`) is deliberately carquet-agnostic (no carquet types leak through), so it reads like an independent product.

It drives **117 of carquet's 143 public functions (~82%)**. The remaining 26 are niche paths (geospatial stats, custom-codec registration, dictionary-preserving reads, FILE*-based create/open, group/map/variant schema builders, the `[offset, limit)` row window, the batch reader's last-error accessor, the low-level column reader, the Arrow import-side read, and the `carquet_error_*` helpers).

## Layout

```
mocklib/
  include/metricstore/metricstore.h   Public API (writer, query, introspection, ...)
  src/metricstore.c                    Schema, writer, query core, introspection
  src/ms_columnar.c                    Buffer I/O, append, low-level column reads, int64 bloom
  src/ms_nested.c                      LIST<STRING> round-trip + rep-level helpers
  src/ms_arrow.c                       Arrow C Data Interface bridge
  src/ms_diag.c                        CPU/version diagnostics, custom allocator
  app/telemetry_demo.c                 End-to-end scenario runner (CLI)
  app/write_sample.c                   Emit a file for external inspection
  tests/test_roundtrip.c               Self-checking round-trip test (37 assertions)
  CMakeLists.txt
```

## What it exercises

The store models this schema (8 leaf columns, mixed physical + logical types):

| # | column       | physical / logical                     | tuning applied                     |
|---|--------------|----------------------------------------|------------------------------------|
| 0 | `event_id`   | INT64                                  | `DELTA_BINARY_PACKED`, bloom       |
| 1 | `ts`         | INT64 / `TIMESTAMP(micros, utc)`       | statistics drive pushdown          |
| 2 | `host`       | BYTE_ARRAY / `STRING`                  | `RLE_DICTIONARY`, bloom            |
| 3 | `region`     | BYTE_ARRAY / `STRING`                  | `RLE_DICTIONARY`                    |
| 4 | `metric`     | BYTE_ARRAY / `STRING`                  | `RLE_DICTIONARY`, bloom            |
| 5 | `value`      | DOUBLE                                 | `BYTE_STREAM_SPLIT`                 |
| 6 | `error_code` | INT32 / `INTEGER(32, signed)` OPTIONAL | nullable (definition levels)       |
| 7 | `session_id` | FIXED_LEN_BYTE_ARRAY(16) / `UUID`      | fixed-length                        |

Across the write and read paths it touches, among others:

- **Write:** `carquet_schema_create` / `add_column` (with logical types) / `set_field_metadata`; `writer_options` (compression, statistics, bloom filters, page index, Arrow schema, CRC); `set_column_encoding` / `set_column_compression` / `set_column_bloom_filter`; `add_metadata`; `write_batch` for required **and** nullable columns; `new_row_group`; `close` / `abort`.
- **Read:** `reader_open` (incl. mmap option); `batch_reader` with **column projection** and a **row-group filter for predicate pushdown**; per-column **statistics**; **bloom-filter** membership (`get_bloom_filter` / `check_bytes`); chunk metadata, key-value metadata, `get_file_info`, `validate_file`.

## Building

By default it builds carquet from the parent source tree via `add_subdirectory`:

```bash
cd mocklib
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

Alternative carquet sources are auto-detected in this order:

1. A `carquet` target already defined by an enclosing project.
2. An installed package — `find_package(carquet CONFIG)`.
3. The source tree at `-DCARQUET_SOURCE_DIR=/path/to/carquet` (default: `..`).

## Running

```bash
# Self-checking round-trip test (also registered with ctest)
./build/test_roundtrip
ctest --test-dir build --output-on-failure

# End-to-end demo: ingest N synthetic events, then describe + query
./build/telemetry_demo /tmp/telemetry.parquet 200000
```

The demo ingests events in streaming chunks (multiple row groups), prints a full
file description (schema, encodings, compression ratios, statistics, embedded
Arrow field labels), validates the file, runs filtered/aggregated queries with
statistics-based row-group pruning, and probes the host bloom filter for present
and absent values.
