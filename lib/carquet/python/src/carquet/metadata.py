"""File / row-group / column-chunk metadata, statistics, page indexes and bloom filters.

Everything here is a snapshot copied out of the footer while the file is open,
so these objects stay valid after the file is closed and need no PyArrow.
"""

from __future__ import annotations

import ctypes
import dataclasses
import struct
from ctypes import byref
from typing import Any, Iterator, Optional, Union

from . import _ffi
from ._types import BOUNDARY_ORDERS, COMPRESSIONS, ENCODINGS, PHYSICAL_TYPES
from ._values import decode_value, encode_value
from .schema import ColumnSchema, ParquetSchema

__all__ = [
    "BloomFilter",
    "ColumnChunkMetaData",
    "ColumnIndex",
    "FileMetaData",
    "GeospatialStatistics",
    "OffsetIndex",
    "PageLocation",
    "PageStatistics",
    "RowGroupMetaData",
    "Statistics",
]


def _copy(ptr: Optional[int], size: int) -> Optional[bytes]:
    # A zero-length value (an empty string minimum) may come with a NULL pointer.
    if size == 0:
        return b""
    if not ptr or size < 0:
        return None
    return ctypes.string_at(ptr, size)


# --------------------------------------------------------------------------
# Statistics
# --------------------------------------------------------------------------
@dataclasses.dataclass(frozen=True)
class Statistics:
    """Column-chunk statistics (``carquet_column_statistics_t``).

    ``min`` / ``max`` are decoded according to the column's logical type
    (``str`` for strings, ``datetime`` for timestamps, ``Decimal`` for
    decimals, ...); ``min_raw`` / ``max_raw`` are the stored physical bytes.
    Counts the writer did not record are ``None``.
    """

    column: ColumnSchema = dataclasses.field(repr=False, compare=False)
    has_min_max: bool
    min_raw: Optional[bytes]
    max_raw: Optional[bytes]
    null_count: Optional[int]
    distinct_count: Optional[int]
    num_values: int

    @property
    def min(self) -> Any:
        """Smallest value of the chunk, decoded with the column's type (``str``, ``datetime``, ``Decimal``...); None when the file stores no min/max."""
        return decode_value(self.min_raw, self.column) if self.has_min_max else None

    @property
    def max(self) -> Any:
        """Largest value of the chunk, decoded with the column's type; None when the file stores no min/max."""
        return decode_value(self.max_raw, self.column) if self.has_min_max else None

    @property
    def has_null_count(self) -> bool:
        """True when the file stores the number of nulls (``null_count`` is then not None)."""
        return self.null_count is not None

    @property
    def has_distinct_count(self) -> bool:
        """True when the file stores the number of distinct values (``distinct_count`` is then not None)."""
        return self.distinct_count is not None

    @property
    def physical_type(self) -> str:
        """Parquet physical type of the column (``"INT64"``, ``"BYTE_ARRAY"``...)."""
        return self.column.physical_type

    def to_dict(self) -> dict:
        """The statistics as a plain dict."""
        return {
            "has_min_max": self.has_min_max,
            "min": self.min,
            "max": self.max,
            "null_count": self.null_count,
            "distinct_count": self.distinct_count,
            "num_values": self.num_values,
        }

    def __repr__(self) -> str:
        mm = f"min={self.min!r}, max={self.max!r}, " if self.has_min_max else ""
        return (
            f"Statistics({mm}null_count={self.null_count}, distinct_count={self.distinct_count}, "
            f"num_values={self.num_values})"
        )


@dataclasses.dataclass(frozen=True)
class GeospatialStatistics:
    """GeospatialStatistics of a GEOMETRY / GEOGRAPHY column chunk.

    ``bbox`` is ``(xmin, ymin, xmax, ymax)`` or ``None``; ``z_range`` /
    ``m_range`` are ``(min, max)`` when the Z / M dimension is present;
    ``geometry_types`` lists the ISO-WKB type codes seen (empty = unknown).
    """

    bbox: Optional[tuple]
    z_range: Optional[tuple]
    m_range: Optional[tuple]
    geometry_types: tuple


# --------------------------------------------------------------------------
# Column chunk / row group / file
# --------------------------------------------------------------------------
@dataclasses.dataclass(frozen=True)
class ColumnChunkMetaData:
    """Metadata of one column chunk (one leaf column within one row group)."""

    column: ColumnSchema = dataclasses.field(repr=False)
    row_group: int
    physical_type: str
    compression: str
    encodings: tuple
    num_values: int
    total_compressed_size: int
    total_uncompressed_size: int
    data_page_offset: int
    dictionary_page_offset: Optional[int]
    has_bloom_filter: bool
    has_column_index: bool
    has_offset_index: bool
    statistics: Optional[Statistics]
    geospatial_statistics: Optional[GeospatialStatistics] = None

    @property
    def path_in_schema(self) -> str:
        """Dotted path of the column (``"address.city"``)."""
        return self.column.path

    @property
    def column_index(self) -> int:
        """Leaf ordinal of this column (``ParquetSchema.column(i)``)."""
        return self.column.index

    @property
    def has_dictionary_page(self) -> bool:
        """True when the chunk is dictionary-encoded (it starts with a dictionary page)."""
        return self.dictionary_page_offset is not None

    @property
    def compression_ratio(self) -> float:
        """``total_uncompressed_size / total_compressed_size`` (1.0 for an empty chunk)."""
        if self.total_compressed_size <= 0:
            return 1.0
        return self.total_uncompressed_size / self.total_compressed_size

    def to_dict(self) -> dict:
        """The chunk metadata as a plain dict (statistics included)."""
        d = {
            "path_in_schema": self.path_in_schema,
            "physical_type": self.physical_type,
            "compression": self.compression,
            "encodings": self.encodings,
            "num_values": self.num_values,
            "total_compressed_size": self.total_compressed_size,
            "total_uncompressed_size": self.total_uncompressed_size,
            "data_page_offset": self.data_page_offset,
            "dictionary_page_offset": self.dictionary_page_offset,
            "has_bloom_filter": self.has_bloom_filter,
            "has_column_index": self.has_column_index,
            "has_offset_index": self.has_offset_index,
            "statistics": self.statistics.to_dict() if self.statistics is not None else None,
        }
        if self.geospatial_statistics is not None:
            d["geospatial_statistics"] = dataclasses.asdict(self.geospatial_statistics)
        return d

    @classmethod
    def _read(cls, handle: int, rg: int, col: ColumnSchema) -> "ColumnChunkMetaData":
        lib = _ffi.get_lib()
        m = _ffi.ColumnChunkMetadataStruct()
        _ffi.check(lib.carquet_reader_column_chunk_metadata(handle, rg, col.index, byref(m)), None,
                   f"row group {rg} column {col.path}")
        n_enc = max(0, min(int(m.num_encodings), 4))
        encodings = tuple(ENCODINGS.get(int(m.encodings[i]), str(int(m.encodings[i]))) for i in range(n_enc))

        stats = None
        s = _ffi.ColumnStatisticsStruct()
        if lib.carquet_reader_column_statistics(handle, rg, col.index, byref(s)) == _ffi.OK:
            if s.has_min_max or s.has_null_count or s.has_distinct_count:
                stats = Statistics(
                    column=col,
                    has_min_max=bool(s.has_min_max),
                    min_raw=_copy(s.min_value, int(s.min_value_size)) if s.has_min_max else None,
                    max_raw=_copy(s.max_value, int(s.max_value_size)) if s.has_min_max else None,
                    null_count=int(s.null_count) if s.has_null_count else None,
                    distinct_count=int(s.distinct_count) if s.has_distinct_count else None,
                    num_values=int(s.num_values),
                )

        geo = None
        if col.logical_type is not None and col.logical_type.name in ("GEOMETRY", "GEOGRAPHY"):
            g = _ffi.GeospatialStatisticsStruct()
            if lib.carquet_reader_geospatial_statistics(handle, rg, col.index, byref(g)) == _ffi.OK:
                n_types = max(0, min(int(g.num_geometry_types), _ffi.MAX_GEOSPATIAL_TYPES))
                geo = GeospatialStatistics(
                    bbox=(g.xmin, g.ymin, g.xmax, g.ymax) if g.has_bbox else None,
                    z_range=(g.zmin, g.zmax) if g.has_z else None,
                    m_range=(g.mmin, g.mmax) if g.has_m else None,
                    geometry_types=tuple(int(g.geometry_types[i]) for i in range(n_types)),
                )

        return cls(
            column=col,
            row_group=rg,
            physical_type=PHYSICAL_TYPES.get(int(m.type), col.physical_type),
            compression=COMPRESSIONS.get(int(m.codec), str(int(m.codec))),
            encodings=encodings,
            num_values=int(m.num_values),
            total_compressed_size=int(m.total_compressed_size),
            total_uncompressed_size=int(m.total_uncompressed_size),
            data_page_offset=int(m.data_page_offset),
            dictionary_page_offset=int(m.dictionary_page_offset) if m.has_dictionary_page else None,
            has_bloom_filter=bool(m.has_bloom_filter),
            has_column_index=bool(m.has_column_index),
            has_offset_index=bool(m.has_offset_index),
            statistics=stats,
            geospatial_statistics=geo,
        )


class RowGroupMetaData:
    """Metadata of one row group: sizes plus one :class:`ColumnChunkMetaData` per leaf."""

    def __init__(self, index: int, num_rows: int, total_byte_size: int, total_compressed_size: int,
                 columns: list, schema: ParquetSchema) -> None:
        self.index = index
        self.num_rows = num_rows
        self.total_byte_size = total_byte_size
        self.total_compressed_size = total_compressed_size
        self._columns = columns
        self._schema = schema

    @property
    def num_columns(self) -> int:
        """Number of leaf columns."""
        return len(self._columns)

    @property
    def columns(self) -> list:
        """One :class:`ColumnChunkMetaData` per leaf column, in schema order."""
        return list(self._columns)

    def column(self, key: Union[int, str]) -> ColumnChunkMetaData:
        """Column chunk by leaf ordinal, dotted path or unique leaf name."""
        if isinstance(key, int) and not isinstance(key, bool):
            return self._columns[key]
        return self._columns[self._schema.index(key)]

    def __getitem__(self, key: Union[int, str]) -> ColumnChunkMetaData:
        return self.column(key)

    def __iter__(self) -> Iterator[ColumnChunkMetaData]:
        return iter(self._columns)

    def __len__(self) -> int:
        return len(self._columns)

    def to_dict(self) -> dict:
        """The row group metadata as a plain dict (its column chunks included)."""
        return {
            "num_rows": self.num_rows,
            "num_columns": self.num_columns,
            "total_byte_size": self.total_byte_size,
            "total_compressed_size": self.total_compressed_size,
            "columns": [c.to_dict() for c in self._columns],
        }

    def __repr__(self) -> str:
        return (
            f"<carquet.RowGroupMetaData index={self.index} num_rows={self.num_rows} "
            f"num_columns={self.num_columns} total_byte_size={self.total_byte_size} "
            f"total_compressed_size={self.total_compressed_size}>"
        )

    @classmethod
    def _read(cls, handle: int, rg: int, schema: ParquetSchema) -> "RowGroupMetaData":
        lib = _ffi.get_lib()
        m = _ffi.RowGroupMetadataStruct()
        _ffi.check(lib.carquet_reader_row_group_metadata(handle, rg, byref(m)), None, f"row group {rg}")
        cols = [ColumnChunkMetaData._read(handle, rg, c) for c in schema]
        return cls(rg, int(m.num_rows), int(m.total_byte_size), int(m.total_compressed_size), cols, schema)


class FileMetaData:
    """Footer metadata of a Parquet file.

    ``row_groups`` / ``row_group(i)`` give per-row-group (and per-column-chunk)
    detail; ``key_value_metadata`` is the footer's key/value map (``bytes``
    keys/values are decoded as UTF-8).
    """

    def __init__(self, *, num_rows: int, schema: ParquetSchema, row_groups: list,
                 key_value_metadata: dict, created_by: Optional[str], format_version: Optional[int],
                 file_size: Optional[int]) -> None:
        self.num_rows = num_rows
        self.schema = schema
        self._row_groups = row_groups
        self.key_value_metadata = key_value_metadata
        self.created_by = created_by
        self.format_version = format_version
        self.file_size = file_size

    @property
    def num_row_groups(self) -> int:
        """Number of row groups in the file."""
        return len(self._row_groups)

    @property
    def num_columns(self) -> int:
        """Number of leaf columns."""
        return len(self.schema)

    @property
    def row_groups(self) -> list:
        """Every :class:`RowGroupMetaData`, in file order."""
        return list(self._row_groups)

    def row_group(self, i: int) -> RowGroupMetaData:
        """The :class:`RowGroupMetaData` of row group ``i`` (negative indices count from the end). Raises ``IndexError`` when out of range."""
        if not -len(self._row_groups) <= i < len(self._row_groups):
            raise IndexError(f"row group {i} out of range (file has {len(self._row_groups)})")
        return self._row_groups[i]

    @property
    def metadata(self) -> dict:
        """Alias of :attr:`key_value_metadata`."""
        return self.key_value_metadata

    def to_dict(self) -> dict:
        """The whole footer as a plain dict (handy with ``json.dumps(..., default=str)``)."""
        return {
            "num_rows": self.num_rows,
            "num_row_groups": self.num_row_groups,
            "num_columns": self.num_columns,
            "created_by": self.created_by,
            "format_version": self.format_version,
            "file_size": self.file_size,
            "key_value_metadata": dict(self.key_value_metadata),
            "row_groups": [rg.to_dict() for rg in self._row_groups],
        }

    def __repr__(self) -> str:
        return (
            f"<carquet.FileMetaData num_rows={self.num_rows} num_row_groups={self.num_row_groups} "
            f"num_columns={self.num_columns} created_by={self.created_by!r}>"
        )


def read_key_value_metadata(handle: int) -> dict:
    lib = _ffi.get_lib()
    out: dict = {}
    for i in range(int(lib.carquet_reader_num_metadata(handle))):
        k = ctypes.c_char_p()
        v = ctypes.c_char_p()
        if lib.carquet_reader_get_metadata(handle, i, byref(k), byref(v)) == _ffi.OK and k.value is not None:
            out[_ffi.decode(k.value)] = _ffi.decode(v.value) if v.value is not None else None
    return out


# --------------------------------------------------------------------------
# Page index
# --------------------------------------------------------------------------
@dataclasses.dataclass(frozen=True)
class PageStatistics:
    """Statistics of one data page, from the column index."""

    column: ColumnSchema = dataclasses.field(repr=False, compare=False)
    min_raw: Optional[bytes]
    max_raw: Optional[bytes]
    null_count: Optional[int]
    is_null_page: bool

    @property
    def min(self) -> Any:
        """Smallest value of the page, decoded with the column's type; None for a page of nulls only."""
        return None if self.is_null_page else decode_value(self.min_raw, self.column)

    @property
    def max(self) -> Any:
        """Largest value of the page, decoded with the column's type; None for a page of nulls only."""
        return None if self.is_null_page else decode_value(self.max_raw, self.column)

    def __repr__(self) -> str:
        return f"PageStatistics(min={self.min!r}, max={self.max!r}, null_count={self.null_count}, is_null_page={self.is_null_page})"


class ColumnIndex:
    """Per-page min/max statistics of a column chunk (Parquet ColumnIndex)."""

    def __init__(self, pages: list, boundary_order: str) -> None:
        self.pages = pages
        self.boundary_order = boundary_order

    def __len__(self) -> int:
        return len(self.pages)

    def __iter__(self) -> Iterator[PageStatistics]:
        return iter(self.pages)

    def __getitem__(self, i: int) -> PageStatistics:
        return self.pages[i]

    @property
    def num_pages(self) -> int:
        """Number of data pages described."""
        return len(self.pages)

    def __repr__(self) -> str:
        return f"<carquet.ColumnIndex num_pages={len(self.pages)} boundary_order={self.boundary_order}>"

    @classmethod
    def _read(cls, handle: int, rg: int, col: ColumnSchema) -> Optional["ColumnIndex"]:
        lib = _ffi.get_lib()
        err = _ffi.CarquetErrorStruct()
        ci = lib.carquet_reader_get_column_index(handle, rg, col.index, byref(err))
        if not ci:
            return None
        try:
            pages = []
            st = _ffi.PageStatsStruct()
            for p in range(int(lib.carquet_column_index_num_pages(ci))):
                _ffi.check(lib.carquet_column_index_get_page_stats(ci, p, byref(st)), None, "column index")
                pages.append(
                    PageStatistics(
                        column=col,
                        min_raw=_copy(st.min_value, int(st.min_value_size)),
                        max_raw=_copy(st.max_value, int(st.max_value_size)),
                        null_count=int(st.null_count) if st.has_null_count else None,
                        is_null_page=bool(st.is_null_page),
                    )
                )
            order = BOUNDARY_ORDERS.get(int(lib.carquet_column_index_boundary_order(ci)), "UNORDERED")
            return cls(pages, order)
        finally:
            lib.carquet_column_index_free(ci)


@dataclasses.dataclass(frozen=True)
class PageLocation:
    """Location of one data page, from the offset index."""

    offset: int
    compressed_size: int
    first_row_index: int


class OffsetIndex:
    """Per-page file locations of a column chunk (Parquet OffsetIndex)."""

    def __init__(self, pages: list) -> None:
        self.pages = pages

    def __len__(self) -> int:
        return len(self.pages)

    def __iter__(self) -> Iterator[PageLocation]:
        return iter(self.pages)

    def __getitem__(self, i: int) -> PageLocation:
        return self.pages[i]

    @property
    def num_pages(self) -> int:
        """Number of data pages described."""
        return len(self.pages)

    def __repr__(self) -> str:
        return f"<carquet.OffsetIndex num_pages={len(self.pages)}>"

    @classmethod
    def _read(cls, handle: int, rg: int, col: ColumnSchema) -> Optional["OffsetIndex"]:
        lib = _ffi.get_lib()
        err = _ffi.CarquetErrorStruct()
        oi = lib.carquet_reader_get_offset_index(handle, rg, col.index, byref(err))
        if not oi:
            return None
        try:
            pages = []
            loc = _ffi.PageLocationStruct()
            for p in range(int(lib.carquet_offset_index_num_pages(oi))):
                _ffi.check(lib.carquet_offset_index_get_page_location(oi, p, byref(loc)), None, "offset index")
                pages.append(PageLocation(int(loc.offset), int(loc.compressed_size), int(loc.first_row_index)))
            return cls(pages)
        finally:
            lib.carquet_offset_index_free(oi)


# --------------------------------------------------------------------------
# Bloom filter
# --------------------------------------------------------------------------
class BloomFilter:
    """Split-block bloom filter of one column chunk.

    ``value in bf`` (or :meth:`might_contain`) is ``False`` only when the value
    is definitely absent from the chunk; ``True`` means "possibly present".
    Values are converted with the column's logical type (``str`` for strings,
    ``datetime`` for timestamps, ...).
    """

    def __init__(self, handle: int, column: ColumnSchema) -> None:
        self._handle = handle
        self.column = column

    def might_contain(self, value: Any) -> bool:
        """False when ``value`` is certainly not in the column chunk, True when it may be (``value in bloom_filter`` does the same)."""
        if not self._handle:
            raise ValueError("bloom filter is closed")
        lib = _ffi.get_lib()
        try:
            raw = encode_value(value, self.column)
        except OverflowError:  # outside the column type's range: cannot be in the chunk
            return False
        pt = self.column.physical_type
        if pt == "INT32":
            return bool(lib.carquet_bloom_filter_check_i32(self._handle, struct.unpack("<i", raw)[0]))
        if pt == "INT64":
            return bool(lib.carquet_bloom_filter_check_i64(self._handle, struct.unpack("<q", raw)[0]))
        if pt == "FLOAT":
            return bool(lib.carquet_bloom_filter_check_float(self._handle, struct.unpack("<f", raw)[0]))
        if pt == "DOUBLE":
            return bool(lib.carquet_bloom_filter_check_double(self._handle, struct.unpack("<d", raw)[0]))
        addr, n, keep = _ffi.buffer_address(raw)
        return bool(lib.carquet_bloom_filter_check_bytes(self._handle, addr or None, n))

    def __contains__(self, value: Any) -> bool:
        return self.might_contain(value)

    @property
    def size(self) -> int:
        """Size of the filter bitset in bytes."""
        return int(_ffi.get_lib().carquet_bloom_filter_size(self._handle)) if self._handle else 0

    def close(self) -> None:
        """Free the filter (also done when the object is garbage-collected)."""
        h, self._handle = self._handle, None
        if h:
            _ffi.get_lib().carquet_bloom_filter_destroy(h)

    def __del__(self) -> None:  # pragma: no cover - best effort
        try:
            self.close()
        except Exception:
            pass

    def __repr__(self) -> str:
        return f"<carquet.BloomFilter column={self.column.path!r} size={self.size}>"
