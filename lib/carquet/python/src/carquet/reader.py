"""Reading Parquet files: :class:`ParquetFile` and the ``read_*`` functions.

Every read is one C Arrow stream (``carquet_reader_export_arrow_stream``):
column projection, row-group selection, row windows, filters (statistics,
page index and the exact predicate), parallel decoding and the assembly of
nested types all run in C, without the GIL. The batches it returns stay in
the buffers C allocated; :class:`carquet.Table` reads them in place and hands
them to any Arrow consumer through the Arrow PyCapsule interface. PyArrow is
not needed.
"""

from __future__ import annotations

import ctypes
import dataclasses
import os
import threading
import weakref
from ctypes import byref, c_int32
from typing import Any, Iterable, Iterator, Optional, Sequence, Union

from . import _ffi
from ._arrow import ArrowArray, ArrowArrayStream, ArrowSchema, release_schema, stream_capsule
from .arrow import RecordBatch, Schema, Table, _stream_error
from .errors import CarquetError, StatusCode
from .filters import _ALL, Filters, compile_filters, normalize_filters
from .metadata import (
    BloomFilter,
    ColumnIndex,
    FileMetaData,
    OffsetIndex,
    RowGroupMetaData,
    read_key_value_metadata,
)
from .schema import ColumnSchema, ParquetSchema

__all__ = [
    "ArrowStream",
    "FileInfo",
    "ParquetFile",
    "info",
    "iter_batches",
    "read_metadata",
    "read_pandas",
    "read_schema",
    "read_table",
    "validate",
]

PathLike = Union[str, "os.PathLike[str]"]
Source = Union[PathLike, bytes, bytearray, memoryview, Any]

DEFAULT_BATCH_SIZE = 65536
# Batch size used when the whole result is materialised (read / read_table):
# large enough that a row group usually comes back as a single chunk.
_READ_BATCH_SIZE = 1 << 20


def _is_path(source: Any) -> bool:
    return isinstance(source, (str, os.PathLike))


# --------------------------------------------------------------------------
# Footer-only helpers (no PyArrow needed)
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
# Top-level fields
# --------------------------------------------------------------------------
class _Field:
    """A top-level field and its leaf columns."""

    __slots__ = ("ordinal", "name", "leaves")

    def __init__(self, ordinal: int, name: str, leaves: list) -> None:
        self.ordinal = ordinal
        self.name = name
        self.leaves = leaves


def _plan_fields(schema: ParquetSchema) -> list:
    fields: list = []
    by_name: dict = {}
    for col in schema:
        top = col.parts[0]
        f = by_name.get(top)
        if f is None:
            f = _Field(len(fields), top, [])
            by_name[top] = f
            fields.append(f)
        f.leaves.append(col)
    return fields


class _Scan:
    """One C Arrow stream over a :class:`ParquetFile`, iterated as RecordBatches."""

    def __init__(self, pf: "ParquetFile", ordinals: list, *, row_groups: Optional[list] = None,
                 window: Optional[tuple] = None, predicate: Any = None, batch_size: int = _READ_BATCH_SIZE,
                 use_threads: bool = True, split: bool = False) -> None:
        self._pf = pf
        self._stream: Optional[ArrowArrayStream] = None
        self._interrupted = False
        # The C stream returns a row group per batch for list / map / struct
        # columns: with `split`, larger batches are cut into views of at most
        # batch_size rows.
        self._split = int(batch_size) if split else 0
        self._pending: Optional[RecordBatch] = None
        lib = _ffi.get_lib()
        opts = _ffi.ArrowStreamOptions()
        lib.carquet_arrow_stream_options_init(byref(opts))
        keep: list = []
        if ordinals:
            idx = (c_int32 * len(ordinals))(*ordinals)
            keep.append(idx)
            opts.field_indices = idx
            opts.num_fields = len(ordinals)
        opts.batch_size = max(1, min(int(batch_size), 2**31 - 1))
        opts.num_threads = pf._num_threads if use_threads else 1
        if row_groups is not None:
            rgs = (c_int32 * max(1, len(row_groups)))(*row_groups)
            keep.append(rgs)
            opts.row_groups = rgs
            opts.num_row_groups = len(row_groups)
        if window is not None:
            opts.row_offset, opts.row_limit = window
        if predicate is not None:
            opts.predicate = ctypes.pointer(predicate.struct)
        opts.restore_arrow_schema = pf._restore
        stream = ArrowArrayStream()
        err = _ffi.CarquetErrorStruct()
        with pf._lock:
            st = lib.carquet_reader_export_arrow_stream(pf._require_open(), byref(opts), byref(stream), byref(err))
        _ffi.check(st, err, pf._name)
        self._stream = stream
        try:
            c_schema = ArrowSchema()
            if stream.get_schema(byref(stream), byref(c_schema)) != 0:
                raise CarquetError(StatusCode.INTERNAL, f"{pf._name}: {_stream_error(stream)}")
            self.schema = Schema(c_schema)
        except BaseException:
            self.close()
            raise
        pf._streams.add(self)

    def __iter__(self) -> "_Scan":
        return self

    def __next__(self) -> RecordBatch:
        if self._pending is not None:
            b = self._pending
            self._pending = None
            return self._cut(b)
        return self._cut(self._next_raw())

    def _cut(self, b: RecordBatch) -> RecordBatch:
        if self._split and b.num_rows > self._split:
            self._pending = b.slice(self._split)
            return b.slice(0, self._split)
        return b

    def _next_raw(self) -> RecordBatch:
        pf = self._pf
        with pf._lock:
            if self._stream is None:
                if self._interrupted:
                    raise CarquetError(StatusCode.ALREADY_CLOSED,
                                       f"{pf._name}: file closed while batches were being read")
                raise StopIteration
            arr = ArrowArray()
            rc = self._stream.get_next(byref(self._stream), byref(arr))
            if rc != 0:
                msg = _stream_error(self._stream)
                self.close()
                raise CarquetError(StatusCode.FILE_READ if rc == 5 else StatusCode.INTERNAL, f"{pf._name}: {msg}")
            if not arr.release:
                self.close()
                raise StopIteration
        return RecordBatch._from_c(self.schema, arr)

    def close(self, interrupted: bool = False) -> None:
        if interrupted and self._stream is not None:
            self._interrupted = True
        s, self._stream = self._stream, None
        if s is not None and s.release:
            s.release(byref(s))

    def __del__(self) -> None:  # pragma: no cover - best effort
        try:
            self.close()
        except Exception:
            pass


# --------------------------------------------------------------------------
# ParquetFile
# --------------------------------------------------------------------------
class ParquetFile:
    """An open Parquet file.

    ``source`` is a path, a bytes-like object (``bytes``, ``bytearray``,
    ``memoryview``, any buffer) or a binary file object; buffers are read in
    place, without copying. Use it as a context manager or call
    :meth:`close`::

        with carquet.ParquetFile("data.parquet") as pf:
            pf.schema                     # carquet.Schema (Arrow types)
            pf.parquet_schema             # carquet.ParquetSchema (physical/logical types)
            pf.metadata.row_group(0).column("id").statistics
            t = pf.read(columns=["id", "name"], filters=[("id", ">", 100)])   # carquet.Table
            for batch in pf.iter_batches(batch_size=10_000):
                ...

    Args:
        memory_map: map the file instead of buffered reads (paths only).
        verify_checksums: verify page CRC32 checksums when present.
        num_threads: decoder threads (0 = one per core).
        buffer_size: read buffer size in bytes for buffered (non-mmap) I/O.
        restore_arrow_schema: apply the Arrow types stored in the file's
            ``ARROW:schema`` entry (time zones, large/dictionary types, field
            metadata...), as ``pyarrow.parquet`` does.

    A ``ParquetFile`` is not thread-safe; open one per thread.
    """

    def __init__(
        self,
        source: Source,
        *,
        memory_map: bool = True,
        verify_checksums: bool = True,
        num_threads: int = 0,
        buffer_size: Optional[int] = None,
        restore_arrow_schema: bool = True,
        use_mmap: Optional[bool] = None,
    ) -> None:
        if use_mmap is not None:  # 0.7.2 spelling
            memory_map = bool(use_mmap)
        self._handle: Optional[int] = None
        self._keepalive: Any = None
        self._buffer_size: Optional[int] = None
        self._buffer_addr: Optional[int] = None
        self._streams: "weakref.WeakSet" = weakref.WeakSet()
        self._lock = threading.RLock()
        self._num_threads = int(num_threads)
        self._restore = bool(restore_arrow_schema)

        lib = _ffi.get_lib()
        opts = _ffi.ReaderOptions()
        lib.carquet_reader_options_init(byref(opts))
        opts.verify_checksums = bool(verify_checksums)
        opts.num_threads = self._num_threads
        if buffer_size is not None:
            opts.buffer_size = int(buffer_size)
        err = _ffi.CarquetErrorStruct()

        if _is_path(source):
            self.path: Optional[str] = os.fsdecode(source)
            opts.use_mmap = bool(memory_map)
            handle = lib.carquet_reader_open(_ffi.fsencode(source), byref(opts), byref(err))
        else:
            self.path = None
            data = source
            if hasattr(source, "getbuffer"):  # io.BytesIO: zero-copy view
                data = source.getbuffer()
            elif hasattr(source, "read") and not isinstance(source, (bytes, bytearray, memoryview)):
                data = source.read()
            try:
                addr, size, keep = _ffi.buffer_address(data)
            except TypeError:
                raise TypeError(
                    "source must be a path, a bytes-like object or a binary file object, "
                    f"got {type(source).__name__}"
                ) from None
            self._keepalive = keep
            self._buffer_size = size
            self._buffer_addr = addr
            opts.use_mmap = False
            handle = lib.carquet_reader_open_buffer(addr or None, size, byref(opts), byref(err))
        if not handle:
            _ffi.check(err.code or StatusCode.FILE_OPEN, err, self.path or "<buffer>")
        self._handle = handle
        self._reader_opts = opts  # to open the independent readers of arrow_stream()
        self._memory_map = bool(lib.carquet_reader_is_mmap(handle))
        self._schema_handle = lib.carquet_reader_schema(handle)

        self.parquet_schema: ParquetSchema = ParquetSchema._from_handle(self._schema_handle)
        self.num_rows: int = int(lib.carquet_reader_num_rows(handle))
        self.num_row_groups: int = int(lib.carquet_reader_num_row_groups(handle))
        self.num_columns: int = int(lib.carquet_reader_num_columns(handle))
        self._fields = _plan_fields(self.parquet_schema)
        self._field_by_name = {f.name: f for f in self._fields}
        self._schema: Optional[Schema] = None
        self._metadata: Optional[FileMetaData] = None
        self._kv: Optional[dict] = None

    # -- lifecycle ---------------------------------------------------------
    @property
    def closed(self) -> bool:
        """True once the file has been closed."""
        return not self._handle

    def _require_open(self) -> int:
        if not self._handle:
            raise CarquetError(StatusCode.ALREADY_CLOSED, f"{self._name}: file is closed")
        return self._handle

    @property
    def _name(self) -> str:
        return self.path or "<buffer>"

    def close(self) -> None:
        """Release the file (and stop any batch iterator still reading it)."""
        with self._lock:
            for s in list(self._streams):
                s.close(interrupted=True)
            h, self._handle = self._handle, None
            if h:
                _ffi.get_lib().carquet_reader_close(h)
            self._keepalive = None

    def __enter__(self) -> "ParquetFile":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def __del__(self) -> None:  # pragma: no cover - best effort
        try:
            self.close()
        except Exception:
            pass

    def __repr__(self) -> str:
        if self.closed:
            return f"<carquet.ParquetFile {self._name!r} closed>"
        return (
            f"<carquet.ParquetFile {self._name!r} rows={self.num_rows} "
            f"row_groups={self.num_row_groups} columns={self.num_columns}>"
        )

    # -- schema & metadata ---------------------------------------------------
    @property
    def schema(self) -> Schema:
        """The Arrow schema of what :meth:`read` returns (a :class:`carquet.Schema`)."""
        if self._schema is None:
            self._schema = _Scan(self, [f.ordinal for f in self._fields], window=(0, 0)).schema
        return self._schema

    @property
    def schema_arrow(self) -> Schema:
        """Alias of :attr:`schema` (``pyarrow.parquet.ParquetFile`` spelling)."""
        return self.schema

    @property
    def key_value_metadata(self) -> dict:
        """The footer's key/value metadata (``str`` -> ``str`` or ``None``)."""
        if self._kv is None:
            self._kv = read_key_value_metadata(self._require_open())
        return self._kv

    @property
    def metadata(self) -> FileMetaData:
        """Footer metadata: row groups, column chunks, statistics (a :class:`FileMetaData`)."""
        if self._metadata is None:
            handle = self._require_open()
            created_by = version = file_size = None
            if self.path is not None:
                try:
                    fi = info(self.path)
                    created_by, version, file_size = fi.created_by, fi.version, fi.file_size
                except CarquetError:
                    pass
            else:
                file_size = self._buffer_size
            self._metadata = FileMetaData(
                num_rows=self.num_rows,
                schema=self.parquet_schema,
                row_groups=[RowGroupMetaData._read(handle, i, self.parquet_schema) for i in range(self.num_row_groups)],
                key_value_metadata=dict(self.key_value_metadata),
                created_by=created_by,
                format_version=version,
                file_size=file_size,
            )
        return self._metadata

    def row_group(self, i: int) -> RowGroupMetaData:
        """Metadata of row group ``i`` (same as ``metadata.row_group(i)``)."""
        return self.metadata.row_group(i)

    def _column(self, column: Union[int, str]) -> ColumnSchema:
        return self.parquet_schema.column(column)

    def _check_rg(self, row_group: int) -> int:
        if not -self.num_row_groups <= row_group < self.num_row_groups:
            raise IndexError(f"row group {row_group} out of range (file has {self.num_row_groups})")
        return row_group % self.num_row_groups if self.num_row_groups else row_group

    def bloom_filter(self, column: Union[int, str], row_group: int = 0) -> Optional[BloomFilter]:
        """The bloom filter of ``column`` in ``row_group``, or ``None`` if none was written."""
        handle = self._require_open()
        col = self._column(column)
        rg = self._check_rg(row_group)
        err = _ffi.CarquetErrorStruct()
        bf = _ffi.get_lib().carquet_reader_get_bloom_filter(handle, rg, col.index, byref(err))
        return BloomFilter(bf, col) if bf else None

    def column_index(self, column: Union[int, str], row_group: int = 0) -> Optional[ColumnIndex]:
        """Per-page statistics (Parquet ColumnIndex) of a column chunk, or ``None``."""
        return ColumnIndex._read(self._require_open(), self._check_rg(row_group), self._column(column))

    def offset_index(self, column: Union[int, str], row_group: int = 0) -> Optional[OffsetIndex]:
        """Per-page file locations (Parquet OffsetIndex) of a column chunk, or ``None``."""
        return OffsetIndex._read(self._require_open(), self._check_rg(row_group), self._column(column))

    # -- projection --------------------------------------------------------
    def _resolve(self, columns: Optional[Sequence[str]]) -> list:
        """Top-level field ordinals for ``columns`` (all fields when ``None``)."""
        if columns is None:
            return [f.ordinal for f in self._fields]
        if isinstance(columns, str):
            raise TypeError("columns must be a list of names, not a str")
        out: list = []
        for name in columns:
            if not isinstance(name, str):
                raise TypeError(f"column names must be str, got {name!r}")
            f = self._field_by_name.get(name)
            if f is None:
                try:
                    top = self.parquet_schema.column(name).parts[0]
                except KeyError:
                    top = name.split(".", 1)[0] if self.parquet_schema.leaves_of(name) else None
                f = self._field_by_name.get(top) if top else None
            if f is None:
                raise KeyError(f"column {name!r} not found in {self._name} (fields: {[x.name for x in self._fields]})")
            if f.ordinal not in out:
                out.append(f.ordinal)
        return out

    def _check_row_groups(self, row_groups: Optional[Iterable[int]]) -> Optional[list]:
        if row_groups is None:
            return None
        out = []
        for i in row_groups:
            rg = self._check_rg(int(i))
            if rg in out:
                raise ValueError(f"row group {rg} listed twice")
            out.append(rg)
        return out

    def _predicate(self, flt: Optional[Filters]):
        """``(predicate or None, no_rows)`` for a normalised filter."""
        if flt is None:
            return None, False
        p = compile_filters(self.parquet_schema, flt)
        if p is _ALL:
            return None, False
        return p, p is None

    def _filter_fields(self, flt: Optional[Filters]) -> list:
        out: list = []
        for name in (flt.columns if flt is not None else []):
            for o in self._resolve([name]):
                if o not in out:
                    out.append(o)
        return out

    def filter_row_groups(self, filters: Any) -> list:
        """Indices of the row groups that may contain rows matching ``filters``.

        Uses only column statistics (no data is read); a returned row group may
        still hold no matching row.
        """
        flt = normalize_filters(filters)
        if flt is None:
            return list(range(self.num_row_groups))
        p = compile_filters(self.parquet_schema, flt, skip_repeated=True)
        if p is _ALL:
            return list(range(self.num_row_groups))
        if p is None:
            return []
        out = (c_int32 * max(1, self.num_row_groups))()
        err = _ffi.CarquetErrorStruct()
        with self._lock:
            n = _ffi.get_lib().carquet_reader_predicate_row_groups(self._require_open(), byref(p.struct), out,
                                                                   self.num_row_groups, byref(err))
        if n < 0:
            _ffi.check(-n, err, self._name)
        return list(out[:n])

    # -- reading -------------------------------------------------------------
    def _empty(self, ordinals: list) -> Table:
        return Table(_Scan(self, ordinals, window=(0, 0)).schema, [])

    def _no_columns(self, rgs: Optional[list], flt: Optional[Filters], window: Optional[tuple] = None) -> Table:
        """A table without columns: only the row count."""
        pred, none = self._predicate(flt)
        schema = self.schema._select([])
        if none:
            return Table(schema, [], 0)
        if rgs == []:
            return Table(schema, [], 0)
        if pred is None:
            md = self.metadata
            if window is not None:
                off, lim = window
                n = max(0, min(self.num_rows - off, lim if lim >= 0 else self.num_rows))
            else:
                n = sum(md.row_group(rg).num_rows for rg in (rgs if rgs is not None else range(self.num_row_groups)))
            return Table(schema, [], n)
        fields = self._filter_fields(flt)
        n = sum(b.num_rows for b in _Scan(self, fields, row_groups=rgs, window=window, predicate=pred))
        return Table(schema, [], n)

    def read(
        self,
        columns: Optional[Sequence[str]] = None,
        *,
        filters: Any = None,
        row_groups: Optional[Iterable[int]] = None,
        use_threads: bool = True,
    ) -> Table:
        """Read the file (or some row groups) into a :class:`carquet.Table`.

        Args:
            columns: top-level fields to read, in the order wanted (a dotted
                path selects its top-level field). ``None`` reads every field.
            filters: row filter in DNF form, e.g. ``[("year", ">=", 2020),
                ("country", "in", ["FR", "DE"])]`` (see :mod:`carquet.filters`).
                Row groups and pages ruled out by statistics are not read; the
                result is filtered exactly.
            row_groups: only read these row groups (in this order).
            use_threads: allow multi-threaded decoding.
        """
        flt = normalize_filters(filters)
        ordinals = self._resolve(columns)
        rgs = self._check_row_groups(row_groups)
        if not ordinals and self._fields:
            return self._no_columns(rgs, flt)
        pred, none = self._predicate(flt)
        if none or rgs == []:
            return self._empty(ordinals)
        scan = _Scan(self, ordinals, row_groups=rgs, predicate=pred, use_threads=use_threads)
        return Table(scan.schema, list(scan))

    def read_row_group(self, i: int, columns: Optional[Sequence[str]] = None, *, use_threads: bool = True) -> Table:
        """Read row group ``i``."""
        return self.read(columns, row_groups=[i], use_threads=use_threads)

    def read_row_groups(self, row_groups: Iterable[int], columns: Optional[Sequence[str]] = None, *,
                        use_threads: bool = True) -> Table:
        """Read the given row groups (in that order) into one table."""
        return self.read(columns, row_groups=list(row_groups), use_threads=use_threads)

    def read_all(self, columns: Optional[Sequence[str]] = None, **kwargs) -> Table:
        """Alias of :meth:`read`."""
        return self.read(columns, **kwargs)

    def iter_batches(
        self,
        batch_size: int = DEFAULT_BATCH_SIZE,
        columns: Optional[Sequence[str]] = None,
        *,
        filters: Any = None,
        row_groups: Optional[Iterable[int]] = None,
        use_threads: bool = True,
    ) -> Iterator[RecordBatch]:
        """Stream the file as :class:`carquet.RecordBatch` es of at most ``batch_size`` rows.

        Batches never span row groups (a selection with list, map or struct
        columns comes one row group per batch). With ``filters``, each batch
        is filtered exactly (so batches may be shorter, and empty ones are
        skipped).
        """
        if int(batch_size) <= 0:
            raise ValueError("batch_size must be a positive integer")
        flt = normalize_filters(filters)
        ordinals = self._resolve(columns)
        rgs = self._check_row_groups(row_groups)
        pred, none = self._predicate(flt)
        if none or rgs == []:
            return iter(())
        if not ordinals and self._fields:
            return self._row_count_batches(rgs, flt, int(batch_size))
        return _Scan(self, ordinals, row_groups=rgs, predicate=pred, batch_size=int(batch_size),
                     use_threads=use_threads, split=True)

    def _row_count_batches(self, rgs: Optional[list], flt: Optional[Filters], batch_size: int) -> Iterator:
        """Batches without columns, carrying only row counts (``columns=[]``)."""
        schema = self.schema._select([])
        pred, _none = self._predicate(flt)
        if pred is None:
            md = self.metadata
            for rg in (rgs if rgs is not None else range(self.num_row_groups)):
                n = md.row_group(rg).num_rows
                for start in range(0, n, batch_size):
                    yield RecordBatch(schema, None, 0, min(batch_size, n - start))
            return
        for b in _Scan(self, self._filter_fields(flt), row_groups=rgs, predicate=pred, batch_size=batch_size,
                       split=True):
            yield RecordBatch(schema, None, 0, b.num_rows)

    def iter_row_groups(self, columns: Optional[Sequence[str]] = None, **kwargs) -> Iterator[Table]:
        """Yield one table per row group."""
        for i in range(self.num_row_groups):
            yield self.read(columns, row_groups=[i], **kwargs)

    def slice(self, offset: int, length: Optional[int] = None, columns: Optional[Sequence[str]] = None, *,
              use_threads: bool = True) -> Table:
        """Read rows ``[offset, offset + length)`` (to the end when ``length`` is None).

        Row groups outside the window are not read; for flat columns, pages
        before ``offset`` are skipped through the offset index when present.
        """
        offset = int(offset)
        if offset < 0:
            offset = max(0, self.num_rows + offset)
        if length is not None and int(length) < 0:
            raise ValueError("length must be >= 0")
        limit = -1 if length is None else int(length)
        ordinals = self._resolve(columns)
        if not ordinals and self._fields:
            return self._no_columns(None, None, (offset, limit))
        scan = _Scan(self, ordinals, window=(offset, limit), use_threads=use_threads)
        return Table(scan.schema, list(scan))

    def head(self, n: int = 5, columns: Optional[Sequence[str]] = None) -> Table:
        """The first ``n`` rows."""
        return self.slice(0, n, columns)

    def count_rows(self, filters: Any = None) -> int:
        """Number of rows, or of rows matching ``filters`` (only the filter columns are read)."""
        flt = normalize_filters(filters)
        if flt is None:
            return self.num_rows
        return self._no_columns(None, flt).num_rows

    # -- Arrow streams -----------------------------------------------------
    def arrow_stream(
        self,
        columns: Optional[Sequence[str]] = None,
        *,
        filters: Any = None,
        row_groups: Optional[Iterable[int]] = None,
        batch_size: int = DEFAULT_BATCH_SIZE,
        use_threads: bool = True,
    ) -> "ArrowStream":
        """The file as an Arrow C stream, for any Arrow consumer.

        The result implements the Arrow PyCapsule interface
        (``__arrow_c_stream__``): ``polars.DataFrame(stream)``,
        ``duckdb.sql("SELECT * FROM stream")``, ``pyarrow.table(stream)``,
        nanoarrow... Unlike :meth:`read`, nothing is materialised: the consumer
        pulls the batches. Every ``__arrow_c_stream__`` call reads the file
        again through its own handle, so the stream stays usable after this
        file is closed.
        """
        self._require_open()
        if int(batch_size) <= 0:
            raise ValueError("batch_size must be a positive integer")
        ordinals = self._resolve(columns)
        if not ordinals and self._fields:
            raise ValueError("arrow_stream needs at least one column")
        flt = normalize_filters(filters)
        pred, none = self._predicate(flt)
        return ArrowStream(self.path, self._buffer_addr, self._buffer_size, self._keepalive, self._reader_opts,
                           ordinals, self._check_row_groups(row_groups), pred, none,
                           min(int(batch_size), 2**31 - 1), self._num_threads if use_threads else 1, self._name,
                           self._restore)

    def __arrow_c_stream__(self, requested_schema: Any = None) -> object:
        """Arrow PyCapsule interface: every column, as :meth:`arrow_stream`."""
        return self.arrow_stream().__arrow_c_stream__(requested_schema)


class ArrowStream:
    """An Arrow C stream over a Parquet file (see :meth:`ParquetFile.arrow_stream`).

    It does not depend on the :class:`ParquetFile` it came from, and each
    ``__arrow_c_stream__`` call opens a fresh read (DuckDB, for one, asks twice).
    """

    def __init__(self, path, buffer_addr, buffer_size, keepalive, reader_opts, ordinals: list,
                 row_groups: Optional[list], predicate: Any, no_rows: bool, batch_size: int, num_threads: int,
                 name: str, restore: bool = True) -> None:
        self._path = path
        self._buffer = (buffer_addr, buffer_size, keepalive)
        self._opts = _ffi.ReaderOptions.from_buffer_copy(reader_opts)
        self._ordinals = list(ordinals)
        self._row_groups = row_groups
        self._predicate = predicate
        self._no_rows = no_rows
        self._batch_size = batch_size
        self._num_threads = num_threads
        self._name = name
        self._restore = restore

    def __arrow_c_stream__(self, requested_schema: Any = None) -> object:
        """Arrow PyCapsule interface. ``requested_schema`` is ignored (the
        protocol lets a producer keep its own types)."""
        lib = _ffi.get_lib()
        err = _ffi.CarquetErrorStruct()
        if self._path is not None:
            handle = lib.carquet_reader_open(_ffi.fsencode(self._path), byref(self._opts), byref(err))
            keep = None
        else:
            addr, size, keep = self._buffer
            handle = lib.carquet_reader_open_buffer(addr or None, size, byref(self._opts), byref(err))
        if not handle:
            _ffi.check(err.code or StatusCode.FILE_OPEN, err, self._name)
        opts = _ffi.ArrowStreamOptions()
        lib.carquet_arrow_stream_options_init(byref(opts))
        idx = (c_int32 * max(1, len(self._ordinals)))(*self._ordinals)
        opts.field_indices = idx
        opts.num_fields = len(self._ordinals)
        opts.batch_size = self._batch_size
        opts.num_threads = self._num_threads
        opts.close_reader = True
        opts.restore_arrow_schema = self._restore
        if self._row_groups is not None:
            rgs = (c_int32 * max(1, len(self._row_groups)))(*self._row_groups)
            opts.row_groups = rgs
            opts.num_row_groups = len(self._row_groups)
        if self._no_rows:
            opts.row_offset, opts.row_limit, opts.row_groups, opts.num_row_groups = 0, 0, None, 0
        elif self._predicate is not None:
            opts.predicate = ctypes.pointer(self._predicate.struct)
        stream = ArrowArrayStream()
        st = lib.carquet_reader_export_arrow_stream(handle, byref(opts), byref(stream), byref(err))
        if st != _ffi.OK:
            lib.carquet_reader_close(handle)
            _ffi.check(st, err, self._name)
        return stream_capsule(stream, keep)

    def read_all(self) -> Table:
        """Read the stream into a :class:`carquet.Table`."""
        from .arrow import table

        return table(self)


# --------------------------------------------------------------------------
# Module-level functions
# --------------------------------------------------------------------------
def _sources(source: Any) -> Optional[list]:
    """Expand a directory or a list of paths into files; ``None`` for a single source."""
    if isinstance(source, (list, tuple)):
        return list(source)
    if _is_path(source) and os.path.isdir(source):
        out = []
        for root, dirs, files in os.walk(source):
            dirs[:] = sorted(d for d in dirs if not d.startswith((".", "_")))
            for name in sorted(files):
                if not name.startswith((".", "_")) and name.endswith((".parquet", ".parq", ".pq")):
                    out.append(os.path.join(root, name))
        return out
    return None


def _concat(tables: list) -> Table:
    first = tables[0]
    for t in tables[1:]:
        if t.schema != first.schema:
            raise ValueError("the files have different schemas and cannot be concatenated")
    batches = [b for t in tables for b in t.to_batches()]
    return Table(first.schema, batches, sum(t.num_rows for t in tables))


def read_table(
    source: Source,
    *,
    columns: Optional[Sequence[str]] = None,
    filters: Any = None,
    row_groups: Optional[Iterable[int]] = None,
    use_threads: bool = True,
    memory_map: bool = True,
    verify_checksums: bool = True,
    restore_arrow_schema: bool = True,
    use_mmap: Optional[bool] = None,
) -> Table:
    """Read a Parquet file into a :class:`carquet.Table`.

    ::

        t = carquet.read_table("data.parquet")
        t = carquet.read_table("data.parquet", columns=["id", "price"],
                               filters=[("price", ">", 100)])
        t.to_pandas()        # or t.to_pydict(), t.column("price").to_numpy(),
                             # polars.DataFrame(t), pyarrow.table(t)

    Args:
        source: a path, a bytes-like object (read in place), a binary file
            object, a directory (every ``*.parquet`` file under it, sorted)
            or a list of any of these. Several files are concatenated; their
            schemas must match.
        columns: names of the columns to read, in the order wanted (default:
            all). Only these columns are read from the file.
        filters: keep only the rows that match, e.g. ``[("year", ">=", 2020),
            ("country", "in", ["FR", "DE"])]``. A list of tuples is an AND, a
            list of lists an OR of ANDs. Operators: ``==``, ``!=``, ``<``,
            ``<=``, ``>``, ``>=``, ``in``, ``not in``, ``is null``,
            ``is not null``. The result is exact, and data that statistics
            rule out is not read at all.
        row_groups: indices of the row groups to read (default: all).
        use_threads: decode columns in parallel (default True).
        memory_map: map the file in memory instead of buffered reads
            (default True; ``use_mmap`` is the old name).
        verify_checksums: check page checksums when the file has them.
        restore_arrow_schema: give columns the exact types they were written
            with (time zones, categoricals...) when the file recorded them
            (default True).

    Raises :class:`CarquetError` when the file cannot be read, ``KeyError``
    for an unknown column.
    """
    many = _sources(source)
    opts = dict(memory_map=memory_map, verify_checksums=verify_checksums,
                restore_arrow_schema=restore_arrow_schema, use_mmap=use_mmap)
    if many is not None:
        if row_groups is not None:
            raise ValueError("row_groups cannot be combined with several source files")
        if not many:
            raise ValueError(f"no Parquet files found in {source!r}")
        return _concat([read_table(s, columns=columns, filters=filters, use_threads=use_threads, **opts)
                        for s in many])
    with ParquetFile(source, **opts) as pf:
        return pf.read(columns, filters=filters, row_groups=row_groups, use_threads=use_threads)


def read_pandas(source: Source, **kwargs):
    """Read a Parquet file into a pandas DataFrame (needs pandas).

    Takes the arguments of :func:`read_table` (``columns``, ``filters``...).
    A file written from a DataFrame gets its index, categoricals and nullable
    dtypes back.
    """
    return read_table(source, **kwargs).to_pandas()


def iter_batches(
    source: Source,
    batch_size: int = DEFAULT_BATCH_SIZE,
    *,
    columns: Optional[Sequence[str]] = None,
    filters: Any = None,
    row_groups: Optional[Iterable[int]] = None,
    use_threads: bool = True,
    memory_map: bool = True,
    verify_checksums: bool = True,
) -> Iterator[RecordBatch]:
    """Stream a Parquet file (or several) as RecordBatches; the file is closed at the end."""
    many = _sources(source)
    for s in many if many is not None else [source]:
        with ParquetFile(s, memory_map=memory_map, verify_checksums=verify_checksums) as pf:
            yield from pf.iter_batches(batch_size, columns, filters=filters, row_groups=row_groups,
                                       use_threads=use_threads)


def read_schema(source: Source) -> Schema:
    """The Arrow schema of a Parquet file, without reading any data."""
    with ParquetFile(source, memory_map=False) as pf:
        return pf.schema


def read_metadata(source: Source) -> FileMetaData:
    """The footer metadata (:class:`FileMetaData`) of a Parquet file, without reading data."""
    with ParquetFile(source, memory_map=False) as pf:
        return pf.metadata


def read_parquet_schema(source: Source) -> ParquetSchema:
    """The Parquet-level schema (:class:`ParquetSchema`) of a file."""
    with ParquetFile(source, memory_map=False) as pf:
        return pf.parquet_schema
