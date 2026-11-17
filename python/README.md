# carquet for Python

Python bindings for [carquet](https://github.com/Vitruves/carquet), a pure C11 Apache Parquet reader and writer. The wheel bundles the C library (zstd, zlib, LZ4, Snappy and Brotli compiled in) and talks to it through `ctypes`; decoding, decompression, filtering, encoding and I/O all run in C, without the GIL. **It has no dependencies**: no PyArrow, no NumPy needed.

```bash
pip install carquet
pip install "carquet[pandas]"   # with pandas + NumPy for read_pandas() / to_pandas() / DataFrame input
```

Requires Python 3.9+.

Reads return a `carquet.Table`: its columns stay in the Arrow buffers the C library decoded into, and it implements the [Arrow PyCapsule interface](https://arrow.apache.org/docs/format/CDataInterface/PyCapsuleInterface.html), so polars, DuckDB, PyArrow, nanoarrow... take it without a copy, and it converts to Python values, NumPy and pandas on its own.

## Contents

- [Quick start](#quick-start) and [common tasks](#common-tasks)
- Reading: [`read_table`](#read_table), [`carquet.Table`](#carquettable), [Arrow streams](#arrow-streams), [`ParquetFile`](#parquetfile), [filters](#filters)
- Writing: [`write_table`](#write_table), [`ParquetWriter`](#parquetwriter)
- [Metadata](#metadata), [schema](#schema), [type fidelity](#type-fidelity)
- [Performance](#performance), [errors](#errors), [other helpers](#other-helpers), [interoperability](#interoperability)
- [Installing from source](#installing-from-source), [supported platforms](#supported-platforms), [limitations](#limitations)

Release notes are in the [CHANGELOG](https://github.com/Vitruves/carquet/blob/main/CHANGELOG.md).

## Quick start

```python
import carquet as cq

t = cq.read_table("data.parquet", columns=["id", "price"], filters=[("price", ">", 100)])   # carquet.Table
t.num_rows, t.column_names, t.schema
t.column("price").to_pylist()              # Python values
t.column("price").to_numpy()               # NumPy (zero-copy for null-free fixed-width columns)
t.to_pandas(); t.to_pydict(); t.to_pylist()
polars.DataFrame(t); duckdb.sql("SELECT * FROM t"); pyarrow.table(t)    # whatever you have installed

df = cq.read_pandas("data.parquet")        # index and dtypes restored from the pandas metadata

cq.write_table({"id": [1, 2, 3], "name": ["a", None, "c"]}, "out.parquet")   # dicts, DataFrames, NumPy...
cq.write_table(polars_df, "out.parquet")   # or any Arrow table: polars, PyArrow, DuckDB, carquet

with cq.ParquetFile("data.parquet") as pf:
    pf.schema                              # carquet.Schema (Arrow types)
    pf.parquet_schema                      # physical / logical types, levels
    pf.metadata.row_group(0).column("price").statistics   # Statistics(min=..., max=..., null_count=...)
    for batch in pf.iter_batches(batch_size=50_000, columns=["id"]):   # carquet.RecordBatch
        ...
```

Everything below is also in the interactive help: `help(carquet)` gives this overview, and `help(carquet.read_table)`, `help(carquet.write_table)`, `help(carquet.Table)`, `help(carquet.ParquetFile)` or `help(carquet.ParquetWriter)` the details (every public function, method and property has a docstring).

## Common tasks

| I want to... | Do this |
|---|---|
| read a whole file | `t = cq.read_table("data.parquet")` |
| read only some columns | `cq.read_table(path, columns=["id", "price"])` |
| read only some rows | `cq.read_table(path, filters=[("price", ">", 100), ("country", "in", ["FR", "DE"])])` |
| get a pandas DataFrame | `cq.read_pandas(path)` or `t.to_pandas()` |
| get NumPy arrays | `t.column("price").to_numpy()` |
| get plain Python values | `t.to_pydict()` (per column) or `t.to_pylist()` (per row) |
| hand the data to polars / DuckDB / PyArrow | `polars.DataFrame(t)`, `duckdb.sql("SELECT * FROM t")`, `pyarrow.table(t)`: no copy |
| read a big file piece by piece | `for batch in cq.iter_batches(path, batch_size=50_000): ...` |
| look at the first rows | `cq.ParquetFile(path).head(10)` |
| see column names and types | `cq.read_schema(path)` |
| see row counts, sizes, statistics | `cq.info(path)`, `cq.read_metadata(path)` |
| write Python lists / NumPy arrays | `cq.write_table({"id": [1, 2], "name": ["a", None]}, "out.parquet")` |
| write a DataFrame (pandas, polars) or any Arrow table | `cq.write_table(df, "out.parquet")` |
| choose the compression | `cq.write_table(data, path, compression="snappy")` (`"zstd"` by default) |
| add rows to an existing file | `cq.write_table(more, path, append=True)` |
| write a big file piece by piece | `with cq.ParquetWriter(path, schema) as w: w.write(batch)` |

## Reading

### `read_table`

```python
cq.read_table(source, *, columns=None, filters=None, row_groups=None, use_threads=True,
              memory_map=True, verify_checksums=True, restore_arrow_schema=True) -> carquet.Table
```

`source` is a path, a bytes-like object (`bytes`, `bytearray`, `memoryview`, any buffer), a binary file object (`io.BytesIO`, an open file), a directory (every `*.parquet` file under it, sorted; hidden and `_`-prefixed entries skipped) or a list of any of these (concatenated). Buffers are read in place, without a copy.

- `columns`: top-level fields to read, in the order wanted. Only the projected columns are read and decoded. A dotted path (`"address.city"`) selects its top-level field.
- `filters`: a row predicate, see [Filters](#filters). The result is exact; row groups and pages that statistics rule out are not read.
- `row_groups`: indices of the row groups to read, in that order.
- `use_threads`: decode columns in parallel.
- `memory_map`: map the file instead of buffered reads.
- `restore_arrow_schema`: apply the Arrow types stored in the file's `ARROW:schema` entry (time zones, dictionary / categorical, large and view types, durations, the null type, fixed-size lists, field metadata), as `pyarrow.parquet` does. The file's key/value metadata is always the schema metadata.

`cq.read_pandas(source, **kwargs)` is `read_table(...).to_pandas()`.

### `carquet.Table`

| Member | |
|---|---|
| `num_rows`, `num_columns`, `shape`, `column_names`, `schema` | sizes and the `carquet.Schema` (`Field`s with `name`, `type` (a `DataType`: `str(t)` is PyArrow's spelling, `t.unit`, `t.tz`, `t.precision`, ...), `nullable`, `metadata`) |
| `column(name or i)`, `t[name]` | a `ChunkedArray`: `to_pylist()`, `to_numpy()`, `t[i]`, `null_count`, `chunks` |
| `select(columns)`, `slice(offset, length)`, `to_batches()` | views, no copy |
| `to_pydict()`, `to_pylist()` | Python values (`datetime` with time zones, `Decimal`, `date`, nested `list` / `dict`) |
| `to_pandas()` | NumPy-backed DataFrame with PyArrow's default dtype mapping; the `pandas` metadata restores the index, categoricals, nullable dtypes and time zones |
| `to_polars()`, `to_pyarrow()` | through the PyCapsule interface (needs polars / PyArrow) |
| `equals(other)` | same schema and values; `other` may be any Arrow table |
| `__arrow_c_stream__`, `__arrow_c_schema__` | Arrow PyCapsule interface |

`cq.table(data)` builds a `Table` from any PyCapsule producer (zero-copy), a `{name: column}` dict or a pandas DataFrame.

### Arrow streams

`pf.arrow_stream(columns=None, *, filters=None, row_groups=None, batch_size=65536, use_threads=True)` (or the `ParquetFile` itself, through `__arrow_c_stream__`) hands batches to a consumer as it pulls them, without materialising the file:

```python
with cq.ParquetFile("data.parquet") as pf:
    s = pf.arrow_stream(["id", "amount"], filters=[("amount", ">", 0)])
duckdb.sql("SELECT sum(amount) FROM s")      # the stream outlives the file; each read starts afresh
```

Batches hold up to `batch_size` rows when no selected column is a list, map or struct, one row group otherwise (`iter_batches()` cuts those down to `batch_size`).

### `ParquetFile`

```python
cq.ParquetFile(source, *, memory_map=True, verify_checksums=True, num_threads=0,
               buffer_size=None, restore_arrow_schema=True)
```

An open file; use it as a context manager or call `close()`. Not thread-safe: open one per thread (separate `ParquetFile` objects can be used concurrently, since the GIL is released).

| Member | Returns |
|---|---|
| `schema` / `schema_arrow` | `carquet.Schema` of what `read()` returns |
| `parquet_schema` | `ParquetSchema` (see [Schema](#schema)) |
| `metadata` | `FileMetaData` (row groups, column chunks, statistics) |
| `key_value_metadata` | footer key/value metadata as a `dict` |
| `num_rows`, `num_row_groups`, `num_columns` | counts (`num_columns` counts leaf columns) |
| `read(columns=None, *, filters=None, row_groups=None, use_threads=True)` | `carquet.Table` |
| `read_row_group(i, columns=None)`, `read_row_groups(indices, columns=None)` | `carquet.Table` |
| `iter_batches(batch_size=65536, columns=None, *, filters=None, row_groups=None)` | iterator of `carquet.RecordBatch` (never spanning row groups) |
| `iter_row_groups(columns=None)` | iterator of `carquet.Table`, one per row group |
| `slice(offset, length=None, columns=None)`, `head(n=5, columns=None)` | `carquet.Table` of a row window; row groups outside it are not read, leading pages are skipped through the offset index when present |
| `arrow_stream(columns=None, *, filters=None, row_groups=None, batch_size=65536)` | an [Arrow stream](#arrow-streams) |
| `count_rows(filters=None)` | row count, or the number of matching rows (only the filter columns are read) |
| `filter_row_groups(filters)` | indices of the row groups that statistics cannot rule out |
| `row_group(i)` | `RowGroupMetaData` |
| `bloom_filter(column, row_group=0)` | `BloomFilter` or `None` |
| `column_index(column, row_group=0)`, `offset_index(column, row_group=0)` | `ColumnIndex` / `OffsetIndex` or `None` |

`cq.iter_batches(source, batch_size, **kwargs)` streams a file (or several) and closes it at the end. Closing a `ParquetFile` stops its iterators.

### Filters

`filters` uses the `pyarrow.parquet` DNF form: a list of `(column, op, value)` tuples is an AND, a list of such lists is an OR of ANDs.

```python
cq.read_table(path, filters=[("year", ">=", 2020), ("country", "in", ["FR", "DE"])])
cq.read_table(path, filters=[[("id", "<", 10)], [("status", "is null")]])        # OR
cq.read_table(path, filters=[("address.zip", "==", "75001")])                    # a struct field
```

Ops: `==` (`=`), `!=`, `<`, `<=`, `>`, `>=`, `in`, `not in`, `is null`, `is not null` (the last two take no value; `("x", "==", None)` means `is null`). Values are converted with the column's logical type: `datetime` / `date` / `numpy.datetime64` / `pandas.Timestamp` for timestamps and dates, `Decimal` for decimals, `str` for strings, `uuid.UUID` for UUIDs. A column inside a list or map cannot be filtered on.

The predicate runs in the C library, in three stages: row groups whose statistics rule it out are skipped; for a single AND of flat columns in a file written with a page index, pages are pruned the same way; every batch is then filtered exactly. Filter columns that are not in `columns` are read and dropped. The semantics are PyArrow's: a comparison involving a null is false, `not in` keeps nulls unless `None` is in the set, NaN equals nothing (and `!=` holds for it). A value the column cannot hold exactly is compared exactly too: `("i", "<", 2.5)` on an integer column, a float64 value against a float32 column, a datetime against a date column or a timestamp finer than the column's unit are rewritten on the neighbouring representable values. PyArrow `Expression` filters are not accepted.

## Writing

### `write_table`

```python
cq.write_table(data, where, *, row_group_size=None, **options) -> None
```

`data` is any Arrow table or stream (a `carquet.Table`, a PyArrow / polars table, a DuckDB result, anything implementing the Arrow PyCapsule interface: written without a Python-side copy), a pandas DataFrame (with the `pandas` metadata that restores its index and dtypes, PyArrow's layout), a `{name: column}` dict (Python lists, NumPy arrays, pandas Series), or a list / iterator of these. Python values map to Arrow types as PyArrow's do: `int` -> int64, `float` -> double, `str` -> string, `bytes` -> binary, `bool`, `datetime` (time zones kept), `date`, `time`, `timedelta`, `Decimal` (precision and scale inferred), `list` -> list, `dict` -> struct, all-`None` -> null. `where` is a path or a writable binary file object. If writing fails, a file created by the call is removed.

| Option | Meaning |
|---|---|
| `compression` | `"zstd"` (default), `"snappy"`, `"gzip"`, `"brotli"`, `"lz4"` (LZ4_RAW, what PyArrow/DuckDB/Spark write), `"lz4_hadoop"` (legacy framed LZ4), `"none"`; or `{column: codec}` (unlisted columns are uncompressed) |
| `compression_level` | an int (zstd 1-22, gzip 1-9, brotli 0-11) or `{column: level}`; `None` = codec default |
| `row_group_size` | rows per row group (default 1 048 576) |
| `data_page_size` | target page size in bytes (default 1 MiB), or `{column: bytes}` |
| `max_rows_per_page` | rows per data page |
| `use_dictionary` | `True` or a list of columns to dictionary-encode; the writer falls back to PLAIN for a chunk whose dictionary does not pay off. Default: carquet's own choice (PLAIN, BYTE_STREAM_SPLIT for floats) |
| `column_encoding` | one encoding for every column or `{column: encoding}`: `"plain"`, `"dictionary"`, `"delta_binary_packed"` (`"delta"`), `"delta_length_byte_array"`, `"delta_byte_array"`, `"byte_stream_split"`, `"rle"` |
| `use_byte_stream_split` | `True` (float columns) or a list of columns |
| `write_statistics` | `True` (default), `False`, or a list of columns |
| `max_statistics_size` | maximum stored size of binary min/max values (bytes, default 32) |
| `write_page_index` | write the column / offset indexes (page-level pruning) |
| `bloom_filters` | `True`, a list of columns, or `{column: {"ndv": n, "fpp": p}}` |
| `sorting_columns` | declared sort order: names or `(name, "descending")` tuples (recorded, not verified) |
| `coerce_timestamps` / `allow_truncated_timestamps` | rescale every timestamp to `"ms"` / `"us"` / `"ns"`; truncation raises unless allowed |
| `data_page_version` | `"1.0"` (default) or `"2.0"` |
| `version` | footer format version: `"2.6"` (default) or `"1.0"` for very old readers |
| `write_page_checksum` (alias `write_crc`) | page CRC32 checksums (default on) |
| `write_batch_size`, `dictionary_pagesize_limit` | internal value-batch size, dictionary page size limit |
| `store_schema` | embed the Arrow schema (`ARROW:schema`, built by the C library) and the schema metadata (e.g. pandas) so readers restore exact Arrow types (default on) |
| `metadata` | extra `{key: value}` footer metadata |
| `append` | add the rows as new row groups to an existing file |

Brotli is compiled into the published wheels. A library built without it raises `CarquetError` with `status == StatusCode.UNSUPPORTED_CODEC`.

### `ParquetWriter`

```python
with cq.ParquetWriter("out.parquet", schema, compression="zstd", row_group_size=500_000) as w:
    for batch in source:
        w.write(batch)                 # Table, RecordBatch, reader, DataFrame, dict...
    w.add_metadata("origin", "etl-v2")
```

`schema` is a `carquet.Schema` or anything implementing `__arrow_c_schema__` (a PyArrow schema...). Takes the same options as `write_table`. Rows are packed into row groups of `row_group_size` rows however they arrive (`flush()` ends the current row group early). The data's columns must be the writer's (`ValueError` otherwise); narrower integers go into int64 columns and float32 into float64 ones, a nullable field without nulls into a non-nullable column, and a null into a non-nullable column raises `CarquetError`. On an exception inside the `with` block the file is discarded (`abort()`); `num_rows` counts the rows written so far.

`cq.ParquetWriter(path, append=True)` opens an existing file and adds row groups after the existing ones; the schema defaults to the file's.

## Metadata

`cq.read_metadata(source)` (or `ParquetFile.metadata`) returns a `FileMetaData`; everything is copied out of the footer, so these objects stay valid after the file is closed.

```python
md = cq.read_metadata("data.parquet")
md.num_rows, md.num_row_groups, md.num_columns, md.created_by, md.format_version, md.file_size
md.key_value_metadata                      # {"ARROW:schema": ..., ...}
rg = md.row_group(0)                       # RowGroupMetaData: num_rows, total_byte_size, total_compressed_size
col = rg.column("price")                   # by leaf index, dotted path or unique leaf name
col.compression, col.encodings, col.total_compressed_size, col.has_dictionary_page, col.has_bloom_filter
col.statistics.min, col.statistics.max, col.statistics.null_count, col.statistics.distinct_count
md.to_dict()
```

Statistics values are decoded with the column's logical type (`str`, `datetime`, `date`, `time`, `Decimal`, `uuid.UUID`, ...); the stored bytes are in `min_raw` / `max_raw`. GEOMETRY / GEOGRAPHY chunks also carry `geospatial_statistics` (bounding box, Z / M ranges, geometry types).

Page indexes and bloom filters:

```python
with cq.ParquetFile("data.parquet") as pf:
    ci = pf.column_index("price", row_group=0)      # per-page min / max / null_count
    oi = pf.offset_index("price", row_group=0)      # per-page offset / size / first row
    bf = pf.bloom_filter("user_id", row_group=0)
    if bf is not None and 12345 not in bf:          # False = definitely absent
        ...
```

## Schema

`ParquetFile.parquet_schema` / `cq.read_parquet_schema(source)` describe the Parquet side of the schema. A `ParquetSchema` is a sequence of `ColumnSchema`, one per leaf column:

```python
s = cq.read_parquet_schema("data.parquet")
len(s), s.names, s.fields                  # leaf count, dotted leaf paths, top-level field names
c = s["amount"]                            # by index, dotted path or unique leaf name
c.physical_type, c.logical_type, c.repetition, c.max_definition_level, c.max_repetition_level
c.logical_type.precision, c.logical_type.scale, str(c.logical_type)    # 'DECIMAL(12, 2)'
s.elements                                 # every node (groups and leaves) in pre-order
print(s)                                   # indented tree
```

`cq.read_schema(source)` returns the Arrow schema (`carquet.Schema`).

## Type fidelity

What `write_table` writes, `read_table` reads back identically (`table.equals`), and `pyarrow.parquet.read_table` agrees. Files written by PyArrow, DuckDB, Spark, polars and others read with the types `pyarrow.parquet` produces:

- DECIMAL columns (fixed-size binary, INT32, INT64 or binary storage) come back as `decimal128` / `decimal256`, INT96 timestamps as `timestamp[ns]` (converted in C).
- The `ARROW:schema` entry restores time zones, `dictionary` (pandas categorical), `large_string` / `large_binary` / `large_list`, string / binary views, durations, the null type, fixed-size lists, decimal widths and field metadata; like PyArrow, `date64`, `time32[s]` and `timestamp[s]` come back in the Parquet unit, and list children keep their Parquet name. The file's other key/value metadata becomes the schema metadata (so `to_pandas()` restores the index).
- `write_table` stores the Arrow schema of the data (built by the C library from the Arrow types), so nested and dictionary schemas round-trip too.
- Types without a direct Parquet mapping are written as their storage type and restored from `ARROW:schema`, as PyArrow does: `duration` as INT64, views as strings / binary, `null` as an all-null INT32 annotated Null. `fixed_size_list` round-trips (PyArrow itself cannot read a `fixed_size_list` column that contains null lists).

Dictionary columns come back with one dictionary per batch, in first-appearance order.

## Performance

- Every C call releases the GIL: separate `ParquetFile` objects read in parallel from Python threads.
- Within a read, columns are decoded in parallel (OpenMP on Linux wheels; a persistent worker pool on macOS / Windows wheels), and batches are exported to Arrow in parallel per column.
- Column projection, row-group selection, statistics / page-index pruning and row windows (`slice`, `head`) are pushed down to C: unselected data is never decompressed.
- Memory-mapped files and in-memory buffers are read in place.
- Reads hand Arrow buffers out in place: `carquet.Table` exports views (reference-counted in C), never copies.
- Writing hands each batch of Arrow input to C without a Python-side copy; pages are encoded in parallel and row groups written on a background thread.

`cq.cpu_info()` lists the SIMD features the library dispatches to.

## Errors

Every failing C call raises `carquet.CarquetError`. It carries `code` (the numeric `carquet_status_t`), `status` (a `carquet.StatusCode` `IntEnum`, e.g. `StatusCode.FILE_NOT_FOUND`), `message` (the library's detailed message, never empty) and, when the library sets them, `offset`, `column_index` and `row_group_index`. Invalid arguments raise the usual `TypeError` / `ValueError` / `KeyError` (unknown column) / `IndexError` (row group out of range).

```python
try:
    cq.info("missing.parquet")
except cq.CarquetError as e:
    print(e.status.name, e.code, e.message)
```

## Other helpers

- `cq.info(path)`: `FileInfo` dataclass (`num_rows`, `num_row_groups`, `num_columns`, `file_size`, `version`, `created_by`) from the footer only.
- `cq.validate(path)`: structural validation (magic bytes, footer, metadata); raises `CarquetError`.
- `cq.__version__`, `cq.lib_version()`: package version and the loaded library's `carquet_version()`, both from `CARQUET_VERSION_STRING` in `include/carquet/carquet.h`.
- `cq.cpu_info()`: runtime SIMD detection (`has_neon`, `has_avx2`, `has_avx512f`, ...).

## Interoperability

Files written by carquet are standard Parquet (format version 2 footer, optional DATA_PAGE_V2, page index, bloom filters), readable by every mainstream implementation:

```python
import pyarrow.parquet as pq
pq.read_table("out.parquet")

import duckdb
duckdb.sql("SELECT count(*) FROM 'out.parquet'")

import polars as pl
pl.read_parquet("out.parquet")
```

## Installing from source

`pip install carquet` takes a prebuilt wheel on the platforms below. Elsewhere, or from a checkout, pip compiles the C library (CMake >= 3.16, a C11 compiler and git are needed):

```bash
pip install ./python          # from the repository root
```

How the wheel is built, the development loop and the release process are in [BUILDING.md](https://github.com/Vitruves/carquet/blob/main/python/BUILDING.md).

## Supported platforms

| Platform | Wheels | Notes |
|---|---|---|
| Linux x86_64 / aarch64 (glibc 2.28+) | yes | SSE/AVX2/AVX-512 or NEON dispatched at runtime; OpenMP |
| macOS 11+ arm64, macOS 10.13+ x86_64 | yes | worker-pool parallelism, no OpenMP |
| Windows x64 | yes | MSVC build, worker-pool parallelism; Windows loader path is covered by CI only |
| musl, 32-bit, PyPy | source build only | not covered by the release tests |

## Limitations

- No Parquet modular encryption (not implemented in the C library).
- Structs and lists nested more than one level deep are assembled one row group at a time (`carquet_reader_read_arrow_fields`); `iter_batches()` slices them afterwards, so peak memory is one row group of those columns. Flat columns and single-level `list` / `map` columns stream at `batch_size`.
- Page-level pruning applies to a single AND of flat columns; OR-filters and nested columns use row-group statistics only (results are exact either way).
- A dotted path in `columns` selects its whole top-level field.
- Unions and run-end-encoded Arrow arrays are rejected on write (no Parquet equivalent); dictionary-encoded Arrow columns are written as their values (and restored as dictionaries through `ARROW:schema`).
- `created_by` is always the library's own string (`Carquet`); the page index records `boundary_order` as `UNORDERED`.
