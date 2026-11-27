"""Writing Parquet files: :class:`ParquetWriter` and :func:`write_table`.

Data reaches the C writer through the Arrow C Data Interface
(``carquet_writer_write_arrow``), one batch at a time. Anything that speaks
the Arrow PyCapsule interface (carquet, PyArrow, polars, DuckDB, nanoarrow...)
is written without a Python-side copy; dicts, lists, NumPy arrays and pandas
DataFrames are converted by :mod:`carquet._build`. Encoding, compression and
I/O run in C without the GIL (pages are encoded in parallel and row groups are
written on a background thread). PyArrow is not needed.
"""

from __future__ import annotations

import ctypes
import os
from ctypes import byref
from typing import Any, Iterable, Iterator, Mapping, Optional, Sequence, Union

from . import _ffi
from ._arrow import ArrowArray, ArrowSchema, take_from_capsule, SCHEMA_CAPSULE
from ._build import from_mapping, from_pandas, is_pandas, table_from_columns
from .arrow import RecordBatch, Schema, StreamReader, Table
from .arrow import table as _to_table
from ._types import (
    ENCODING_NAMES,
    RLE_DICTIONARY,
    codec_id,
    encoding_id,
    time_unit_id,
)
from .errors import CarquetError, StatusCode
from .schema import ParquetSchema

__all__ = ["ParquetWriter", "write_table"]

DEFAULT_ROW_GROUP_SIZE = 1 << 20  # rows, like pyarrow.parquet.write_table


# --------------------------------------------------------------------------
# Input normalisation
# --------------------------------------------------------------------------
def to_batches(data: Any) -> tuple:
    """``(Schema, iterator of RecordBatch)`` for any supported tabular input."""
    if isinstance(data, Table):
        return data.schema, iter(data.to_batches())
    if isinstance(data, RecordBatch):
        return data.schema, iter([data])
    if is_pandas(data):
        cols, md = from_pandas(data)
        t = table_from_columns(cols, md)
        return t.schema, iter(t.to_batches())
    if hasattr(data, "__arrow_c_stream__"):  # PyArrow, polars, DuckDB, nanoarrow, ArrowStream...
        r = StreamReader(data)
        return r.schema, r
    if hasattr(data, "__arrow_c_array__"):
        t = _to_table(data)
        return t.schema, iter(t.to_batches())
    if isinstance(data, Mapping):
        t = from_mapping(data)
        return t.schema, iter(t.to_batches())
    if isinstance(data, (list, tuple)) or hasattr(data, "__next__"):
        it = iter(data)
        try:
            first = next(it)
        except StopIteration:
            raise ValueError("cannot infer a schema from an empty sequence of batches") from None
        schema, first_batches = to_batches(first)
        rest = (b for item in it for b in to_batches(item)[1])
        return schema, _chain(first_batches, rest)
    raise TypeError(
        "expected a carquet Table / RecordBatch, an object implementing the Arrow PyCapsule interface "
        "(PyArrow, polars, DuckDB...), a pandas DataFrame, a dict of columns or a list of batches; "
        f"got {type(data).__name__}"
    )


def _chain(first: Iterable, rest: Iterable) -> Iterator:
    yield from first
    yield from rest


def _schema_of(schema: Any) -> Schema:
    if isinstance(schema, Schema):
        return schema
    if hasattr(schema, "__arrow_c_schema__"):
        return Schema(take_from_capsule(schema.__arrow_c_schema__(), SCHEMA_CAPSULE, ArrowSchema))
    raise TypeError(f"schema must be a carquet.Schema or implement __arrow_c_schema__, got {type(schema).__name__}")


def _same_names(a: Schema, b: Schema) -> bool:
    """The data's columns line up with the writer's (types are checked by the
    C import, which also widens int32 -> int64 and float -> double)."""
    return a.names == b.names


def _conform_nullability(src: ArrowSchema, ref: ArrowSchema, depth: int = 0) -> None:
    """Give the data's schema the writer's nullability flags: a nullable field
    without nulls goes into a REQUIRED column (the C import refuses an actual
    null there), a non-nullable one into an OPTIONAL column."""
    src.flags = (src.flags & ~2) | (ref.flags & 2)
    if depth < 64:
        for i in range(min(src.n_children, ref.n_children)):
            _conform_nullability(src.children[i].contents, ref.children[i].contents, depth + 1)


def _format_version(version: Any) -> int:
    if version is None:
        return 2
    s = str(version)
    if s in ("1", "1.0"):
        return 1
    if s.startswith("2"):
        return 2
    raise ValueError(f"version must be '1.0' or '2.x', got {version!r}")


def _page_version(v: Any) -> int:
    s = str(v)
    if s in ("1", "1.0"):
        return 1
    if s in ("2", "2.0"):
        return 2
    raise ValueError(f"data_page_version must be 1 or 2 ('1.0' / '2.0'), got {v!r}")


class _Columns:
    """Resolves user column names to the writer's leaf column ordinals."""

    def __init__(self, schema: ParquetSchema) -> None:
        self.schema = schema

    def leaves(self, name: str) -> list:
        try:
            return [self.schema.column(name).index]
        except KeyError:
            pass
        hits = [c.index for c in self.schema.leaves_of(name)]
        if not hits:
            raise KeyError(f"column {name!r} not found (leaf columns: {self.schema.names})")
        return hits

    def select(self, spec: Union[bool, Sequence[str], None]) -> Optional[set]:
        """Leaf set for a ``bool | list[str]`` option (``None`` = option not used)."""
        if spec is None or spec is False:
            return None
        if spec is True:
            return {c.index for c in self.schema}
        if isinstance(spec, str):
            spec = [spec]
        out: set = set()
        for name in spec:
            out.update(self.leaves(name))
        return out


def _per_column(value: Any, columns: _Columns) -> tuple:
    """``(global_value, {leaf: value})`` for options that take a value or a ``{column: value}`` dict."""
    if isinstance(value, Mapping):
        per: dict = {}
        for name, v in value.items():
            for leaf in columns.leaves(name):
                per[leaf] = v
        return None, per
    return value, {}


def _sorting(spec: Any, columns: _Columns) -> list:
    out = []
    for item in spec:
        nulls_first = False
        if hasattr(item, "column_index") and hasattr(item, "descending"):  # pyarrow.parquet.SortingColumn
            idx, desc, nulls_first = int(item.column_index), bool(item.descending), bool(getattr(item, "nulls_first", False))
        else:
            if isinstance(item, str):
                name, order = item, "ascending"
            elif isinstance(item, (tuple, list)) and 2 <= len(item) <= 3:
                name, order = item[0], item[1]
                if len(item) == 3:
                    np = item[2]
                    nulls_first = np in ("at_start", "first", True)
            else:
                raise ValueError(f"sorting column must be a name or (name, 'ascending'|'descending'[, nulls]); got {item!r}")
            if isinstance(order, bool):
                desc = order
            elif str(order).lower() in ("ascending", "asc"):
                desc = False
            elif str(order).lower() in ("descending", "desc"):
                desc = True
            else:
                raise ValueError(f"sort order must be 'ascending' or 'descending', got {order!r}")
            leaves = columns.leaves(name)
            if len(leaves) != 1:
                raise ValueError(f"sorting column {name!r} must name a single leaf column")
            idx = leaves[0]
        out.append((idx, desc, nulls_first))
    return out


# --------------------------------------------------------------------------
# ParquetWriter
# --------------------------------------------------------------------------
class ParquetWriter:
    """Write a Parquet file piece by piece (:func:`write_table` writes one in a single call).

    ::

        with carquet.ParquetWriter("out.parquet", schema, compression="zstd") as w:
            for batch in batches:
                w.write(batch)          # Table, RecordBatch, DataFrame, dict...

    ``where`` is a path or a writable binary file object (the file is then
    built in memory and written to it on :meth:`close`). ``schema`` is a
    :class:`carquet.Schema` (``table.schema``, ``carquet.read_schema(path)``)
    or any Arrow schema, PyArrow's included. Every piece written must have the
    schema's columns. Rows are packed into row groups of ``row_group_size``
    rows however they arrive. Leaving the ``with`` block closes the file; if
    the block raises, the file is discarded.

    With ``append=True`` the file at ``where`` must exist; new row groups are
    added after the existing ones (``schema`` defaults to the file's schema).

    Options (all keyword-only; names follow ``pyarrow.parquet.write_table``):

    * ``compression``: ``"zstd"`` (default), ``"snappy"``, ``"gzip"``, ``"brotli"``,
      ``"lz4"`` (LZ4_RAW), ``"lz4_hadoop"``, ``"none"``; or ``{column: codec}``
      (unlisted columns are uncompressed). ``compression_level``: an int or a
      ``{column: level}`` dict.
    * ``row_group_size``: rows per row group (default 1 048 576).
    * ``data_page_size``: target page size in bytes, or ``{column: bytes}``.
    * ``max_rows_per_page``, ``write_batch_size``, ``dictionary_pagesize_limit``.
    * ``use_dictionary``: ``True`` / list of columns to dictionary-encode
      (default: carquet's PLAIN encoding). ``column_encoding``: an encoding for
      every column or ``{column: encoding}`` -- ``"plain"``, ``"dictionary"``,
      ``"delta_binary_packed"``, ``"delta_length_byte_array"``,
      ``"delta_byte_array"``, ``"byte_stream_split"``, ``"rle"``.
      ``use_byte_stream_split``: ``True`` (float columns) or a list of columns.
    * ``write_statistics``: ``True`` / ``False`` / list of columns;
      ``max_statistics_size`` caps stored min/max of binary columns (bytes).
    * ``write_page_index``: write the column / offset indexes (page pruning).
    * ``bloom_filters``: ``True``, a list of columns, or ``{column: {"ndv": n,
      "fpp": p}}``.
    * ``sorting_columns``: declared sort order -- names, ``(name, "descending")``
      tuples or ``pyarrow.parquet.SortingColumn`` objects (not verified).
    * ``coerce_timestamps`` (``"ms"`` / ``"us"`` / ``"ns"``) and
      ``allow_truncated_timestamps``.
    * ``data_page_version``: ``"1.0"`` (default) or ``"2.0"``; ``version``: the
      footer format version (``"2.6"`` default, ``"1.0"`` for very old readers).
    * ``write_page_checksum`` (alias ``write_crc``): page CRC32 (default on).
    * ``store_schema``: embed the Arrow schema (``ARROW:schema``) and the
      schema metadata (e.g. pandas) so readers restore exact Arrow types.
    * ``metadata``: extra ``{key: value}`` footer metadata.
    """

    def __init__(
        self,
        where: Any,
        schema: Any = None,
        *,
        append: bool = False,
        compression: Any = "zstd",
        compression_level: Any = None,
        row_group_size: Optional[int] = None,
        data_page_size: Any = None,
        data_page_version: Any = 1,
        version: Any = "2.6",
        use_dictionary: Any = None,
        column_encoding: Any = None,
        use_byte_stream_split: Any = None,
        write_statistics: Any = True,
        max_statistics_size: Optional[int] = None,
        write_page_index: bool = False,
        write_page_checksum: Optional[bool] = None,
        write_crc: Optional[bool] = None,
        bloom_filters: Any = None,
        sorting_columns: Optional[Sequence[Any]] = None,
        coerce_timestamps: Optional[str] = None,
        allow_truncated_timestamps: bool = False,
        write_batch_size: Optional[int] = None,
        max_rows_per_page: Optional[int] = None,
        dictionary_pagesize_limit: Optional[int] = None,
        store_schema: Optional[bool] = None,
        metadata: Optional[Mapping[str, Any]] = None,
        page_size: Any = None,
    ) -> None:
        lib = _ffi.get_lib()
        self._writer: Optional[int] = None
        self._sink = None
        self.path: Optional[str] = None
        self._created_file = False
        self.num_rows = 0
        self._rows_in_group = 0

        if _is_path_like(where):
            self.path = os.fsdecode(where)
        elif hasattr(where, "write"):
            if append:
                raise ValueError("append=True needs a path")
            self._sink = where
        else:
            raise TypeError(f"where must be a path or a writable binary file object, got {type(where).__name__}")

        if schema is None:
            if not append:
                raise TypeError("schema is required (or use write_table / append=True)")
            from .reader import ParquetFile

            with ParquetFile(self.path, memory_map=False) as pf:
                schema = pf.schema
        schema = _schema_of(schema)
        self.schema = schema
        self.row_group_size = int(row_group_size) if row_group_size else DEFAULT_ROW_GROUP_SIZE
        if self.row_group_size <= 0:
            raise ValueError("row_group_size must be a positive integer")
        if store_schema is None:
            store_schema = not append

        # Arrow schema -> carquet schema (the import consumes the ArrowSchema).
        c_schema = schema._copy()
        cq_schema = ctypes.c_void_p()
        err = _ffi.CarquetErrorStruct()
        _ffi.check(lib.carquet_arrow_import_schema(byref(c_schema), byref(cq_schema), byref(err)), err, "schema")
        try:
            self.parquet_schema = ParquetSchema._from_handle(cq_schema)
            cols = _Columns(self.parquet_schema)

            opts = _ffi.WriterOptions()
            lib.carquet_writer_options_init(byref(opts))
            codec, per_codec = _per_column(compression, cols)
            opts.compression = codec_id(codec) if codec is not None or not per_codec else 0
            level, per_level = _per_column(compression_level, cols)
            opts.compression_level = int(level) if level is not None else 0
            if page_size is not None and data_page_size is None:
                data_page_size = page_size
            psize, per_psize = _per_column(data_page_size, cols)
            if psize is not None:
                opts.page_size = int(psize)
            opts.data_page_version = _page_version(data_page_version)
            opts.file_format_version = _format_version(version)
            opts.write_statistics = write_statistics is not False
            opts.write_page_index = bool(write_page_index)
            crc = write_page_checksum if write_page_checksum is not None else write_crc
            if crc is not None:
                opts.write_crc = bool(crc)
            # ARROW:schema is built from the Arrow schema itself below (exact
            # types: time zones, dictionaries, large / view types, durations,
            # child names, metadata), not from the Parquet schema.
            opts.write_arrow_schema = False
            if coerce_timestamps is not None:
                opts.coerce_timestamps = True
                opts.coerce_timestamp_unit = time_unit_id(coerce_timestamps)
            opts.allow_timestamp_truncation = bool(allow_truncated_timestamps)
            if write_batch_size is not None:
                opts.write_batch_size = int(write_batch_size)
            if max_rows_per_page is not None:
                opts.max_rows_per_page = int(max_rows_per_page)
                # the C writer checks the row cap once per write batch: make
                # the batches no larger than the cap so pages hold exactly it
                if write_batch_size is None and int(max_rows_per_page) > 0:
                    opts.write_batch_size = int(max_rows_per_page)
            if dictionary_pagesize_limit is not None:
                opts.dictionary_page_size = int(dictionary_pagesize_limit)

            err = _ffi.CarquetErrorStruct()
            if self._sink is not None:
                writer = lib.carquet_writer_create_buffer(cq_schema, byref(opts), byref(err))
            elif append:
                writer = lib.carquet_writer_open_append(_ffi.fsencode(self.path), cq_schema, byref(opts), byref(err))
            else:
                writer = lib.carquet_writer_create(_ffi.fsencode(self.path), cq_schema, byref(opts), byref(err))
                self._created_file = bool(writer)
            if not writer:
                _ffi.check(err.code or StatusCode.FILE_OPEN, err, self.path or "<buffer>")
        finally:
            lib.carquet_schema_free(cq_schema)
        self._writer = writer

        try:
            self._configure_columns(cols, per_codec, codec, per_level, level, per_psize, use_dictionary,
                                    column_encoding, use_byte_stream_split, write_statistics,
                                    max_statistics_size, bloom_filters, sorting_columns)
            if store_schema:
                st = lib.carquet_writer_set_arrow_schema(self._writer, byref(schema._c))
                if st == StatusCode.NOT_IMPLEMENTED:  # a type with no IPC form: the Parquet types only
                    pass
                else:
                    _ffi.check(st, None, "ARROW:schema")
            for k, v in (metadata or {}).items():
                self.add_metadata(k, v)
        except BaseException:
            self._abort()
            raise

    def _configure_columns(self, cols, per_codec, codec, per_level, level, per_psize, use_dictionary,
                           column_encoding, use_byte_stream_split, write_statistics, max_statistics_size,
                           bloom_filters, sorting_columns) -> None:
        lib = _ffi.get_lib()
        w = self._writer

        def check(st, what):
            _ffi.check(st, None, what)

        # encodings: use_dictionary < use_byte_stream_split < column_encoding
        enc: dict = {}
        dict_cols = cols.select(use_dictionary)
        if dict_cols:
            for leaf in dict_cols:
                if use_dictionary is not True or cols.schema[leaf].physical_type != "BOOLEAN":
                    enc[leaf] = RLE_DICTIONARY
        if use_byte_stream_split:
            if use_byte_stream_split is True:
                bss = {c.index for c in cols.schema if c.physical_type in ("FLOAT", "DOUBLE")}
            else:
                bss = cols.select(use_byte_stream_split) or set()
            for leaf in bss:
                enc[leaf] = ENCODING_NAMES["byte_stream_split"]
        if column_encoding is not None:
            g, per = _per_column(column_encoding, cols)
            if g is not None:
                for c in cols.schema:
                    enc[c.index] = encoding_id(g)
            for leaf, e in per.items():
                enc[leaf] = encoding_id(e)
        for leaf, e in sorted(enc.items()):
            check(lib.carquet_writer_set_column_encoding(w, leaf, e), f"encoding of {cols.schema[leaf].path}")

        for leaf in sorted(set(per_codec) | set(per_level)):
            c = codec_id(per_codec[leaf]) if leaf in per_codec else codec_id(codec)
            lvl = per_level.get(leaf, level)
            check(lib.carquet_writer_set_column_compression(w, leaf, c, int(lvl) if lvl is not None else 0),
                  f"compression of {cols.schema[leaf].path}")
        for leaf, size in per_psize.items():
            check(lib.carquet_writer_set_column_page_size(w, leaf, int(size)), f"page size of {cols.schema[leaf].path}")

        if write_statistics is not True and write_statistics is not False:
            keep = cols.select(write_statistics) or set()
            for c in cols.schema:
                if c.index not in keep:
                    check(lib.carquet_writer_set_column_statistics(w, c.index, False), "statistics")
        if max_statistics_size is not None:
            check(lib.carquet_writer_set_max_statistics_size(w, int(max_statistics_size)), "max_statistics_size")

        if bloom_filters:
            if isinstance(bloom_filters, Mapping):
                items = []
                for name, spec in bloom_filters.items():
                    for leaf in cols.leaves(name):
                        items.append((leaf, spec))
            else:
                items = [(leaf, True) for leaf in sorted(cols.select(bloom_filters) or ())]
            for leaf, spec in items:
                if spec is False or spec is None:
                    continue
                ndv, fpp = 0, 0.0
                if isinstance(spec, Mapping):
                    ndv, fpp = int(spec.get("ndv", 0) or 0), float(spec.get("fpp", 0.0) or 0.0)
                elif isinstance(spec, (tuple, list)) and len(spec) == 2:
                    ndv, fpp = int(spec[0] or 0), float(spec[1] or 0.0)
                check(lib.carquet_writer_set_column_bloom_filter_options(w, leaf, True, ndv, fpp), "bloom filter")

        if sorting_columns:
            entries = _sorting(sorting_columns, cols)
            arr = (_ffi.SortingColumnStruct * len(entries))()
            for i, (idx, desc, nf) in enumerate(entries):
                arr[i].column_index, arr[i].descending, arr[i].nulls_first = idx, desc, nf
            check(lib.carquet_writer_set_sorting_columns(w, arr, len(entries)), "sorting_columns")

    # -- writing -------------------------------------------------------------
    @property
    def closed(self) -> bool:
        """True once the writer has been closed or aborted."""
        return not self._writer

    def _require_open(self) -> int:
        if not self._writer:
            raise CarquetError(StatusCode.ALREADY_CLOSED, "writer is closed")
        return self._writer

    def add_metadata(self, key: str, value: Optional[str]) -> None:
        """Add a footer key/value entry (before :meth:`close`)."""
        w = self._require_open()
        k = key.encode("utf-8") if isinstance(key, str) else bytes(key)
        if value is None:
            v = None
        elif isinstance(value, bytes):
            v = value
        else:
            v = str(value).encode("utf-8")
        _ffi.check(_ffi.get_lib().carquet_writer_add_metadata(w, k, v), None, f"metadata {key!r}")

    def write(self, data: Any, row_group_size: Optional[int] = None) -> None:
        """Write a table, record batch, Arrow stream (PyArrow, polars, DuckDB...),
        DataFrame or dict of columns."""
        self._require_open()
        schema, batches = to_batches(data)
        if not _same_names(schema, self.schema):
            raise ValueError(f"data columns {schema.names} do not match the writer schema {self.schema.names}")
        rg_size = int(row_group_size) if row_group_size else self.row_group_size
        for batch in batches:
            self._write_batch(batch, rg_size)

    write_table = write
    write_batch = write

    def _write_batch(self, batch: RecordBatch, rg_size: int) -> None:
        if batch is None or batch.num_rows == 0:
            return
        start, n = 0, batch.num_rows
        while start < n:
            room = rg_size - self._rows_in_group
            take = min(room, n - start)
            part = batch if (start == 0 and take == n) else batch.slice(start, take)
            self._write_arrow(part)
            start += take
            self._rows_in_group += take
            self.num_rows += take
            if self._rows_in_group >= rg_size:
                self.flush()

    def _write_arrow(self, batch: RecordBatch) -> None:
        lib = _ffi.get_lib()
        array = batch._export()
        c_schema = batch.schema._copy()
        _conform_nullability(c_schema, self.schema._c)
        err = _ffi.CarquetErrorStruct()
        st = lib.carquet_writer_write_arrow(self._writer, byref(array), byref(c_schema), byref(err))
        _ffi.check(st, err, "write")

    def flush(self) -> None:
        """End the current row group now (no-op when it is empty)."""
        w = self._require_open()
        if self._rows_in_group:
            _ffi.check(_ffi.get_lib().carquet_writer_new_row_group(w), None, "row group flush")
            self._rows_in_group = 0

    def close(self) -> None:
        """Finish the file (footer) and release the writer. Idempotent."""
        w, self._writer = self._writer, None
        if not w:
            return
        lib = _ffi.get_lib()
        st = lib.carquet_writer_close(w)
        if self._sink is None:
            _ffi.check(st, None, f"close {self.path}")
            return
        if st != _ffi.OK:
            lib.carquet_writer_abort(w)
            _ffi.check(st, None, "close")
        ptr = ctypes.c_void_p()
        size = ctypes.c_size_t()
        st = lib.carquet_writer_get_buffer(w, byref(ptr), byref(size))
        if st != _ffi.OK:
            lib.carquet_writer_abort(w)
            _ffi.check(st, None, "get buffer")
        try:
            view = (ctypes.c_char * size.value).from_address(ptr.value)
            self._sink.write(memoryview(view).cast("B"))
        finally:
            _ffi.free_buffer(ptr.value)

    def _abort(self) -> None:
        w, self._writer = self._writer, None
        if w:
            _ffi.get_lib().carquet_writer_abort(w)
        if self._created_file and self.path:
            try:
                os.remove(self.path)
            except OSError:
                pass

    def abort(self) -> None:
        """Discard the file being written (a new file is removed; nothing is finalised)."""
        self._abort()

    def __enter__(self) -> "ParquetWriter":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if exc_type is None:
            self.close()
        else:
            self._abort()

    def __del__(self) -> None:  # pragma: no cover - best effort
        # like pyarrow.parquet.ParquetWriter: an unclosed writer is finished, not discarded
        try:
            if self._writer:
                self.close()
        except Exception:
            pass

    def __repr__(self) -> str:
        target = self.path or type(self._sink).__name__
        state = "closed" if self.closed else f"rows={self.num_rows}"
        return f"<carquet.ParquetWriter {target!r} {state}>"


def _is_path_like(obj: Any) -> bool:
    return isinstance(obj, (str, os.PathLike))


def write_table(table: Any, where: Any, *, row_group_size: Optional[int] = None, **options) -> None:
    """Write data to a Parquet file.

    ::

        carquet.write_table({"id": [1, 2, 3], "name": ["a", None, "c"]}, "out.parquet")
        carquet.write_table(df, "out.parquet", compression="snappy")     # pandas / polars
        carquet.write_table(more_rows, "out.parquet", append=True)

    Args:
        table: what to write. A ``{column name: values}`` dict (Python lists,
            NumPy arrays, pandas Series), a pandas DataFrame, a
            :class:`carquet.Table`, or any Arrow table or stream (PyArrow,
            polars, DuckDB...). In lists, None is a null; ``int``, ``float``,
            ``str``, ``bytes``, ``bool``, ``datetime``, ``date``, ``time``,
            ``timedelta``, ``Decimal``, ``list`` and ``dict`` values are
            understood.
        where: a path, or a writable binary file object (``io.BytesIO``...).
        row_group_size: rows per row group (default 1 048 576).
        **options: the most useful ones are ``compression`` (``"zstd"`` by
            default, ``"snappy"``, ``"gzip"``, ``"brotli"``, ``"lz4"``,
            ``"none"``), ``compression_level``, ``use_dictionary`` (True or a
            list of columns with few distinct values), ``write_page_index``
            and ``bloom_filters`` (faster filtered reads), ``metadata`` (a
            ``{key: value}`` dict stored in the file) and ``append`` (add the
            rows to an existing file). The full list is in
            ``help(carquet.ParquetWriter)``.

    If writing fails, a file created by this call is removed.
    """
    schema, batches = to_batches(table)
    w = ParquetWriter(where, None if options.get("append") else schema, row_group_size=row_group_size, **options)
    try:
        if not _same_names(schema, w.schema):
            raise ValueError(f"data columns {schema.names} do not match the file schema {w.schema.names}")
        for b in batches:
            w._write_batch(b, w.row_group_size)
        w.close()
    except BaseException:
        w._abort()
        raise
