"""carquet -- Apache Parquet reading and writing backed by the carquet C library.

Tables cross the C boundary through the Arrow C Data Interface, so
:func:`read_table` / :func:`write_table` move PyArrow data without any
Python-side decoding. :func:`info`, :func:`validate`, :func:`lib_version` and
:data:`__version__` do not need PyArrow at all.
"""

from __future__ import annotations

import ctypes
import dataclasses
import os
from ctypes import addressof, byref
from typing import TYPE_CHECKING, Any, Iterable, Iterator, Optional, Sequence, Union

from . import _ffi
from ._arrow import ArrowArray, ArrowSchema, release_array, release_schema
from .errors import CarquetError, StatusCode

if TYPE_CHECKING:  # pragma: no cover
    import pyarrow as pa

__all__ = [
    "COMPRESSION_CODECS",
    "CarquetError",
    "FileInfo",
    "Reader",
    "StatusCode",
    "__version__",
    "info",
    "lib_version",
    "open_reader",
    "read_table",
    "schema",
    "validate",
    "write_table",
]

PathLike = Union[str, "os.PathLike[str]", bytes]


def lib_version() -> str:
    """Version string reported by the loaded C library (``carquet_version()``)."""
    return _ffi.get_lib().carquet_version().decode("ascii")


def _package_version() -> str:
    try:
        from importlib.metadata import version

        return version("carquet")
    except Exception:  # not installed as a distribution (e.g. run from a checkout)
        try:
            return lib_version()
        except Exception:
            return "0.0.0"


__version__ = _package_version()

# Parquet compression codec ids (carquet_compression_t). "lz4" is the modern
# LZ4_RAW codec (7) that PyArrow/DuckDB/Spark write; "lz4_hadoop" is the legacy
# Hadoop-framed codec (5).
COMPRESSION_CODECS = {
    "none": 0,
    "uncompressed": 0,
    "snappy": 1,
    "gzip": 2,
    "brotli": 4,
    "lz4_hadoop": 5,
    "zstd": 6,
    "lz4": 7,
    "lz4_raw": 7,
}

_DEFAULT_ROW_GROUP_ROWS = 1 << 20  # same default as pyarrow.parquet.write_table


# --------------------------------------------------------------------------
# PyArrow helpers
# --------------------------------------------------------------------------
def _pa():
    try:
        import pyarrow as pa
    except ImportError as exc:  # pragma: no cover - exercised only without pyarrow
        raise ImportError(
            "pyarrow is required for this carquet function; install it with "
            "`pip install pyarrow` or `pip install 'carquet[arrow]'`"
        ) from exc
    return pa


def _import_batch(pa, array: ArrowArray, c_schema: ArrowSchema) -> "pa.RecordBatch":
    """Move a top-level struct ArrowArray/ArrowSchema pair into a RecordBatch.

    PyArrow takes ownership of both structs (it marks them released).
    """
    try:
        return pa.RecordBatch._import_from_c(addressof(array), addressof(c_schema))
    except BaseException:
        release_array(array)
        release_schema(c_schema)
        raise


def _import_schema(pa, c_schema: ArrowSchema) -> "pa.Schema":
    try:
        return pa.Schema._import_from_c(addressof(c_schema))
    except BaseException:
        release_schema(c_schema)
        raise


# --------------------------------------------------------------------------
# Metadata helpers (no PyArrow needed)
# --------------------------------------------------------------------------
@dataclasses.dataclass(frozen=True)
class FileInfo:
    """Footer-level facts about a Parquet file (``carquet_file_info_t``)."""

    path: str
    num_rows: int
    num_row_groups: int
    num_columns: int
    file_size: int
    version: int
    created_by: str


def info(path: PathLike) -> FileInfo:
    """Read only the footer of ``path`` and return a :class:`FileInfo`.

    Raises :class:`CarquetError` if the file is missing or not a Parquet file.
    """
    lib = _ffi.get_lib()
    out = _ffi.FileInfoStruct()
    err = _ffi.CarquetErrorStruct()
    st = lib.carquet_get_file_info(_ffi.fsencode(path), byref(out), byref(err))
    _ffi.check(st, err, str(os.fsdecode(path)))
    return FileInfo(
        path=os.fsdecode(path),
        num_rows=int(out.num_rows),
        num_row_groups=int(out.num_row_groups),
        num_columns=int(out.num_columns),
        file_size=int(out.file_size),
        version=int(out.version),
        created_by=out.created_by.decode("utf-8", "replace"),
    )


def validate(path: PathLike) -> None:
    """Structurally validate ``path`` (magic, footer, metadata); raise :class:`CarquetError` on failure."""
    lib = _ffi.get_lib()
    err = _ffi.CarquetErrorStruct()
    st = lib.carquet_validate_file(_ffi.fsencode(path), byref(err))
    _ffi.check(st, err, str(os.fsdecode(path)))


# --------------------------------------------------------------------------
# Reader
# --------------------------------------------------------------------------
class Reader:
    """An open Parquet file. Use as a context manager or call :meth:`close`.

    Each row group is materialised as one :class:`pyarrow.RecordBatch` by
    ``carquet_reader_read_arrow()`` (nested struct/list/map columns included).
    """

    def __init__(
        self,
        path: PathLike,
        *,
        use_mmap: bool = True,
        verify_checksums: bool = True,
        num_threads: int = 0,
    ) -> None:
        lib = _ffi.get_lib()
        opts = _ffi.ReaderOptions()
        lib.carquet_reader_options_init(byref(opts))
        opts.use_mmap = bool(use_mmap)
        opts.verify_checksums = bool(verify_checksums)
        opts.num_threads = int(num_threads)
        err = _ffi.CarquetErrorStruct()
        self.path = os.fsdecode(path)
        self._handle = lib.carquet_reader_open(_ffi.fsencode(path), byref(opts), byref(err))
        if not self._handle:
            _ffi.check(err.code or StatusCode.FILE_OPEN, err, self.path)
        self.num_rows: int = int(lib.carquet_reader_num_rows(self._handle))
        self.num_row_groups: int = int(lib.carquet_reader_num_row_groups(self._handle))
        self.num_columns: int = int(lib.carquet_reader_num_columns(self._handle))
        self._schema: Optional[Any] = None

    # -- lifecycle ---------------------------------------------------------
    @property
    def closed(self) -> bool:
        return not self._handle

    def _require_open(self) -> None:
        if not self._handle:
            raise CarquetError(StatusCode.ALREADY_CLOSED, f"{self.path}: reader is closed")

    def close(self) -> None:
        h, self._handle = self._handle, None
        if h:
            _ffi.get_lib().carquet_reader_close(h)

    def __enter__(self) -> "Reader":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def __del__(self) -> None:  # pragma: no cover - best effort
        try:
            self.close()
        except Exception:
            pass

    def __repr__(self) -> str:
        state = "closed" if self.closed else f"rows={self.num_rows} row_groups={self.num_row_groups} columns={self.num_columns}"
        return f"<carquet.Reader {self.path!r} {state}>"

    # -- data --------------------------------------------------------------
    @property
    def schema(self) -> "pa.Schema":
        """The file schema as a :class:`pyarrow.Schema` (``carquet_arrow_export_schema``)."""
        if self._schema is None:
            self._require_open()
            pa = _pa()
            lib = _ffi.get_lib()
            c_schema = ArrowSchema()
            err = _ffi.CarquetErrorStruct()
            st = lib.carquet_arrow_export_schema(lib.carquet_reader_schema(self._handle), byref(c_schema), byref(err))
            _ffi.check(st, err, self.path)
            self._schema = _import_schema(pa, c_schema)
        return self._schema

    def read_row_group(self, index: int) -> "pa.RecordBatch":
        """Read row group ``index`` (0-based) as a :class:`pyarrow.RecordBatch`."""
        self._require_open()
        pa = _pa()
        if not 0 <= index < self.num_row_groups:
            raise IndexError(f"row group {index} out of range (file has {self.num_row_groups})")
        lib = _ffi.get_lib()
        array = ArrowArray()
        c_schema = ArrowSchema()
        err = _ffi.CarquetErrorStruct()
        st = lib.carquet_reader_read_arrow(self._handle, int(index), byref(c_schema), byref(array), byref(err))
        _ffi.check(st, err, self.path)
        return _import_batch(pa, array, c_schema)

    def iter_batches(self) -> Iterator["pa.RecordBatch"]:
        """Yield every row group in order."""
        for i in range(self.num_row_groups):
            yield self.read_row_group(i)

    def read_all(self, columns: Optional[Sequence[str]] = None) -> "pa.Table":
        """Read the whole file into a :class:`pyarrow.Table` (one chunk per row group)."""
        pa = _pa()
        batches = list(self.iter_batches())
        table = pa.Table.from_batches(batches, schema=self.schema if not batches else None)
        if columns is not None:
            table = table.select(list(columns))
        return table


def open_reader(path: PathLike, *, use_mmap: bool = True, verify_checksums: bool = True, num_threads: int = 0) -> Reader:
    """Open ``path`` for reading; the returned :class:`Reader` is a context manager."""
    return Reader(path, use_mmap=use_mmap, verify_checksums=verify_checksums, num_threads=num_threads)


def read_table(path: PathLike, *, columns: Optional[Sequence[str]] = None, use_mmap: bool = True) -> "pa.Table":
    """Read a Parquet file into a :class:`pyarrow.Table`.

    ``columns`` selects (and orders) top-level fields by name after the read;
    every row group becomes one chunk of the resulting table.
    """
    with Reader(path, use_mmap=use_mmap) as reader:
        return reader.read_all(columns)


def schema(path: PathLike) -> "pa.Schema":
    """Return the Arrow schema of ``path`` without reading any data."""
    with Reader(path, use_mmap=False) as reader:
        return reader.schema


# --------------------------------------------------------------------------
# Writer
# --------------------------------------------------------------------------
def _codec_id(compression: Union[str, int, None]) -> int:
    if compression is None:
        return 0
    if isinstance(compression, int):
        return compression
    key = compression.lower()
    try:
        return COMPRESSION_CODECS[key]
    except KeyError:
        raise ValueError(f"unknown compression {compression!r}; expected one of {sorted(set(COMPRESSION_CODECS))}") from None


def _iter_write_batches(pa, table: "pa.Table", row_group_size: Optional[int]) -> Iterator["pa.RecordBatch"]:
    """Slice ``table`` into row groups of at most ``row_group_size`` rows."""
    step = int(row_group_size) if row_group_size else _DEFAULT_ROW_GROUP_ROWS
    if step <= 0:
        raise ValueError("row_group_size must be a positive integer")
    n = table.num_rows
    for start in range(0, n, step):
        chunk = table.slice(start, step).combine_chunks()
        yield from chunk.to_batches()


def write_table(
    table: "Union[pa.Table, pa.RecordBatch]",
    path: PathLike,
    *,
    compression: Union[str, int, None] = "zstd",
    compression_level: Optional[int] = None,
    row_group_size: Optional[int] = None,
    write_page_index: bool = False,
    write_statistics: bool = True,
    data_page_version: int = 1,
    store_schema: bool = True,
    write_crc: bool = True,
) -> None:
    """Write a :class:`pyarrow.Table` (or RecordBatch) to ``path``.

    Args:
        compression: ``"none"``/``"uncompressed"``, ``"snappy"``, ``"gzip"``, ``"lz4"`` (LZ4_RAW),
            ``"lz4_hadoop"``, ``"zstd"`` (default) or ``"brotli"``. A codec the library was built
            without raises :class:`CarquetError` with ``StatusCode.UNSUPPORTED_CODEC``.
        compression_level: codec level (zstd 1-22, gzip 1-9, brotli 0-11); ``None`` = codec default.
        row_group_size: maximum rows per row group (default 1 048 576, like PyArrow).
        write_page_index: emit the column/offset page index.
        write_statistics: emit min/max statistics.
        data_page_version: 1 (DATA_PAGE) or 2 (DATA_PAGE_V2).
        store_schema: embed the Arrow schema (``ARROW:schema`` footer entry) so PyArrow restores exact Arrow types; flat schemas only.
        write_crc: write page CRC32 checksums.
    """
    pa = _pa()
    if isinstance(table, pa.RecordBatch):
        table = pa.Table.from_batches([table])
    if not isinstance(table, pa.Table):
        raise TypeError(f"expected pyarrow.Table or pyarrow.RecordBatch, got {type(table).__name__}")
    codec = _codec_id(compression)
    if data_page_version not in (1, 2):
        raise ValueError("data_page_version must be 1 or 2")

    lib = _ffi.get_lib()
    err = _ffi.CarquetErrorStruct()

    # 1. Arrow schema -> carquet schema (import consumes the ArrowSchema).
    c_schema = ArrowSchema()
    table.schema._export_to_c(addressof(c_schema))
    cq_schema = ctypes.c_void_p()
    st = lib.carquet_arrow_import_schema(byref(c_schema), byref(cq_schema), byref(err))
    _ffi.check(st, err, "schema import")

    # 2. Writer options + writer (the writer copies the schema).
    writer = None
    try:
        opts = _ffi.WriterOptions()
        lib.carquet_writer_options_init(byref(opts))
        opts.compression = codec
        opts.compression_level = int(compression_level) if compression_level is not None else 0
        opts.write_statistics = bool(write_statistics)
        opts.write_page_index = bool(write_page_index)
        opts.write_crc = bool(write_crc)
        opts.data_page_version = int(data_page_version)
        opts.write_arrow_schema = bool(store_schema)
        err = _ffi.CarquetErrorStruct()
        writer = lib.carquet_writer_create(_ffi.fsencode(path), cq_schema, byref(opts), byref(err))
        if not writer:
            _ffi.check(err.code or StatusCode.FILE_OPEN, err, os.fsdecode(path))
    finally:
        lib.carquet_schema_free(cq_schema)

    # 3. One row group per batch. write_arrow consumes both structs.
    try:
        for batch in _iter_write_batches(pa, table, row_group_size):
            array = ArrowArray()
            b_schema = ArrowSchema()
            batch._export_to_c(addressof(array), addressof(b_schema))
            err = _ffi.CarquetErrorStruct()
            st = lib.carquet_writer_write_arrow(writer, byref(array), byref(b_schema), byref(err))
            _ffi.check(st, err, "write")
            _ffi.check(lib.carquet_writer_new_row_group(writer), None, "row group flush")
        w, writer = writer, None  # close() invalidates the handle even on failure
        _ffi.check(lib.carquet_writer_close(w), None, "close")
    except BaseException:
        if writer is not None:
            lib.carquet_writer_abort(writer)
        raise
