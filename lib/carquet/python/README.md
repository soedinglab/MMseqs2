# carquet for Python

Python bindings for [carquet](https://github.com/Vitruves/carquet), a pure C11 Apache Parquet reader and writer. The wheel bundles the C library (with zstd, zlib, LZ4, Snappy and Brotli compiled in) and talks to it through `ctypes`; tables cross the boundary via the [Arrow C Data Interface](https://arrow.apache.org/docs/format/CDataInterface.html), so a `pyarrow.Table` is produced or consumed without any Python-side decoding or copying of the Parquet logic.

```bash
pip install carquet            # library only: info(), validate(), lib_version()
pip install "carquet[arrow]"   # + pyarrow, needed for read_table / write_table / schema / open_reader
```

Requires Python 3.9+. PyArrow (>= 14) is an optional runtime dependency: it is imported lazily by the table-level functions only.

## Quick start

```python
import pyarrow as pa
import carquet

table = pa.table({"id": pa.array(range(5), pa.int64()), "name": ["a", "b", None, "d", "e"]})

carquet.write_table(table, "out.parquet", compression="zstd")
back = carquet.read_table("out.parquet")          # pyarrow.Table
assert back.equals(table)

fi = carquet.info("out.parquet")                  # footer only, no pyarrow needed
print(fi.num_rows, fi.num_row_groups, fi.num_columns, fi.file_size, fi.created_by)

carquet.validate("out.parquet")                   # raises carquet.CarquetError on a damaged file
print(carquet.__version__, carquet.lib_version())
```

## API

### `carquet.read_table(path, *, columns=None, use_mmap=True) -> pyarrow.Table`

Reads every row group with `carquet_reader_read_arrow()`; each row group becomes one chunk of the returned table, nested `struct` / `list` / `map` columns included. `columns` selects (and orders) top-level fields by name after the read. `use_mmap=False` reads through buffered I/O instead of a memory map.

### `carquet.write_table(table, path, *, compression="zstd", compression_level=None, row_group_size=None, write_page_index=False, write_statistics=True, data_page_version=1, store_schema=True, write_crc=True) -> None`

Writes a `pyarrow.Table` (or a single `RecordBatch`). The Arrow schema is converted with `carquet_arrow_import_schema()` and every row group is handed over as a struct `ArrowArray` to `carquet_writer_write_arrow()`.

| Argument | Meaning |
|---|---|
| `compression` | `"none"`/`"uncompressed"`, `"snappy"`, `"gzip"`, `"lz4"` (LZ4_RAW, what PyArrow/DuckDB/Spark write), `"lz4_hadoop"` (legacy framed LZ4), `"zstd"` (default), `"brotli"` |
| `compression_level` | codec level (zstd 1-22, gzip 1-9, brotli 0-11); `None` uses the codec default |
| `row_group_size` | maximum rows per row group; default 1 048 576 like `pyarrow.parquet.write_table` |
| `write_page_index` | emit the column index / offset index (page-level pruning for readers that support it) |
| `write_statistics` | min/max statistics per column chunk |
| `data_page_version` | `1` (DATA_PAGE) or `2` (DATA_PAGE_V2) |
| `store_schema` | embed the Arrow schema as the `ARROW:schema` footer entry (flat schemas) so PyArrow restores exact Arrow types |
| `write_crc` | page CRC32 checksums |

Brotli is compiled into the published wheels. A build of the library without it makes `write_table(..., compression="brotli")` raise `CarquetError` with `status == StatusCode.UNSUPPORTED_CODEC`.

### `carquet.schema(path) -> pyarrow.Schema`

The file schema via `carquet_arrow_export_schema()`, without reading data.

### `carquet.info(path) -> FileInfo`

Dataclass with `path`, `num_rows`, `num_row_groups`, `num_columns` (leaf columns), `file_size`, `version` (Parquet format version in the footer) and `created_by`. Only the footer is read (`carquet_get_file_info`). Works without PyArrow.

### `carquet.validate(path) -> None`

Structural validation (`carquet_validate_file`: magic bytes, footer, metadata). Raises `CarquetError` on failure. Works without PyArrow.

### `carquet.open_reader(path, *, use_mmap=True, verify_checksums=True, num_threads=0) -> Reader`

Keeps the file open for row-group-level access:

```python
with carquet.open_reader("out.parquet") as r:
    print(r.num_rows, r.num_row_groups, r.num_columns, r.schema)
    first = r.read_row_group(0)            # pyarrow.RecordBatch
    for batch in r.iter_batches():          # all row groups, in order
        ...
    whole = r.read_all(columns=["id"])      # pyarrow.Table
```

### Errors

Every failing C call raises `carquet.CarquetError`. It carries `code` (the numeric `carquet_status_t`), `status` (a `carquet.StatusCode` `IntEnum`, e.g. `StatusCode.FILE_NOT_FOUND`), `message` (the library's detailed message, never empty) and, when the library sets them, `offset`, `column_index` and `row_group_index`.

```python
try:
    carquet.info("missing.parquet")
except carquet.CarquetError as e:
    print(e.status.name, e.code, e.message)
```

### Versions

`carquet.__version__` is the package version; `carquet.lib_version()` is `carquet_version()` from the loaded shared library. Both come from `CARQUET_VERSION_STRING` in `include/carquet/carquet.h` at build time.

## Interoperability

Files written by `carquet.write_table` are standard Parquet (format version 2 footer, PLAIN / RLE_DICTIONARY encodings, optional DATA_PAGE_V2 and page index), readable by every mainstream implementation:

```python
import pyarrow.parquet as pq
pq.read_table("out.parquet")

import duckdb
duckdb.sql("SELECT count(*) FROM 'out.parquet'")

import polars as pl
pl.read_parquet("out.parquet")
```

The reverse direction works the same: `carquet.read_table` opens files produced by PyArrow, DuckDB, Spark, polars and others (the test suite round-trips through PyArrow in both directions).

## How the library is bundled

The wheel contains the pure-Python package and one shared library placed next to `__init__.py`: `carquet/libcarquet.so` (Linux), `carquet/libcarquet.dylib` (macOS) or `carquet/carquet.dll` (Windows). It is built from the repository's own `CMakeLists.txt` with `CARQUET_BUILD_SHARED=ON`, `CARQUET_BUILD_CLI=OFF`, `CARQUET_BUILD_DEV=OFF` and `CARQUET_BUNDLE_DEPS=ON`, so zstd, zlib, LZ4 and Brotli are fetched at build time and linked statically into that single file (Snappy and LZ4_RAW are implemented inside carquet). The result depends only on the C runtime and system libraries. Linux wheels keep OpenMP (parallel column reads, page-parallel encoding); `auditwheel` grafts `libgomp` into the wheel. The macOS and Windows wheels are built without OpenMP so they carry no `libomp`/`vcomp` dependency.

The loader (`carquet/_ffi.py`) looks for the library next to the package, then in `lib/`, `lib64/` and `bin/` under it (layouts a plain `cmake --install` produces), then on the system path. Set the `CARQUET_LIBRARY` environment variable to an explicit file to override.

## Building the wheel locally

From the repository root, with CMake >= 3.16, a C11 compiler and git (for the codec sources) available:

```bash
pip wheel ./python -w dist/          # or: python -m build python/ --outdir dist/
pip install dist/carquet-*.whl
pytest python/tests                  # needs pytest + pyarrow
```

`pip install -e ./python` is not supported (the shared library must be installed into the package); use a regular install. If your shell exports `CC`, `CFLAGS`, `LDFLAGS` or `SDKROOT` for another toolchain, clear them for the build (`env -u CC -u CFLAGS -u CPPFLAGS -u LDFLAGS -u SDKROOT pip wheel ./python -w dist/`).

The build is driven by `python/CMakeLists.txt`, a thin wrapper that adds the repository root as a subdirectory and installs only the runtime library (CMake component `carquet_python`) into the package directory. `pyproject.toml` (scikit-build-core) does the rest.

### Source distribution

The sdist root is `python/`, so the C tree one level up cannot be referenced as `..` once unpacked. `tool.scikit-build.sdist.force-include` therefore copies `CMakeLists.txt`, `LICENSE`, `include/`, `src/` and `cmake/` into `native/` inside the sdist, `python/CMakeLists.txt` prefers `native/` when it exists, and a `[[tool.scikit-build.overrides]]` entry (`if.from-sdist = true`) points the version regex at `native/include/carquet/carquet.h`. `python -m build --sdist python/` produces it; `pip install carquet-<version>.tar.gz` compiles it like the git checkout.

### Release wheels

`.github/workflows/wheels.yml` builds wheels with [cibuildwheel](https://cibuildwheel.pypa.io/) for CPython 3.9 - 3.13 on Linux (x86_64, aarch64; manylinux_2_28), macOS (arm64, x86_64) and Windows (AMD64), runs this test suite against each wheel, builds the sdist and publishes everything to PyPI via trusted publishing when a `v*` tag is pushed.

## Supported platforms

| Platform | Wheels | Notes |
|---|---|---|
| Linux x86_64 / aarch64 (glibc 2.28+) | yes | SSE/AVX2/AVX-512 or NEON dispatched at runtime; OpenMP on |
| macOS 11+ arm64, macOS 10.13+ x86_64 | yes | no OpenMP |
| Windows x64 | yes | MSVC build, no OpenMP; Windows loader path is covered by CI only |
| musl, 32-bit, PyPy | source build only | no PyArrow wheels exist there to test against |

## Limitations

- No Parquet modular encryption (not implemented in the C library).
- Nested data (struct / list / map) goes through the Arrow bridge (`carquet_reader_read_arrow` / `carquet_writer_write_arrow`), which assembles each row group as a fresh copy; the zero-copy batch export in the C API is not used by the Python package.
- `columns=` is applied after the row group is read (no column projection in the Arrow read path yet).
- Row groups are materialised one at a time; there is no streaming/iterator API below row-group granularity.
- Unions and run-end-encoded Arrow arrays are rejected on write (no Parquet equivalent); dictionary-encoded Arrow columns are written as their values.
- `created_by` is always the library's own string (`Carquet <version>`).
