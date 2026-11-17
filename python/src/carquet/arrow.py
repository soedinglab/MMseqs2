"""carquet's in-memory tables: Arrow data, without PyArrow.

Every array lives in buffers the C library allocated (Arrow C Data
Interface); Python holds a reference and reads them in place. The objects
implement the Arrow PyCapsule interface, so any Arrow library takes them
without a copy::

    t = carquet.read_table("data.parquet")      # carquet.Table
    polars.DataFrame(t)                          # __arrow_c_stream__
    duckdb.sql("SELECT * FROM t")
    pyarrow.table(t)                             # if PyArrow is installed
    t.to_pandas(), t.to_pydict(), t.column("x").to_numpy()
"""

from __future__ import annotations

import ctypes
from ctypes import byref, c_int32
from typing import Any, Iterator, Optional, Sequence, Union

from . import _ffi
from ._arrow import (
    ARRAY_CAPSULE,
    SCHEMA_CAPSULE,
    STREAM_CAPSULE,
    ArrowArray,
    ArrowArrayStream,
    ArrowSchema,
    make_capsule,
    release_schema,
    stream_capsule,
    take_from_capsule,
)
from .errors import CarquetError, StatusCode

__all__ = ["DataType", "Field", "Schema", "Array", "ChunkedArray", "RecordBatch", "Table", "table"]


# --------------------------------------------------------------------------
# Types
# --------------------------------------------------------------------------
_SIMPLE = {
    "n": "null", "b": "bool", "c": "int8", "C": "uint8", "s": "int16", "S": "uint16",
    "i": "int32", "I": "uint32", "l": "int64", "L": "uint64", "e": "halffloat",
    "f": "float", "g": "double", "z": "binary", "Z": "large_binary", "u": "string",
    "U": "large_string", "vz": "binary_view", "vu": "string_view",
    "tdD": "date32[day]", "tdm": "date64[ms]", "tts": "time32[s]", "ttm": "time32[ms]",
    "ttu": "time64[us]", "ttn": "time64[ns]", "tiM": "month_interval",
    "tiD": "day_time_interval", "tin": "month_day_nano_interval",
}
_UNITS = {"s": "s", "m": "ms", "u": "us", "n": "ns"}


class DataType:
    """The type of a column: ``int64``, ``string``, ``timestamp[us, tz=UTC]``...

    ``str(t)`` gives the usual Arrow spelling and ``t.id`` the type family
    (``"timestamp"``, ``"decimal128"``, ``"list"``...). Details are in
    properties that are None when they do not apply: ``unit`` and ``tz``
    (temporal types), ``precision`` and ``scale`` (decimals), ``value_type``
    (lists and dictionaries), ``fields`` (structs), ``key_type`` and
    ``item_type`` (maps). ``t.format`` is the Arrow C Data Interface format
    string. Types compare equal to PyArrow's (``t == pyarrow.int64()``).
    """

    __slots__ = ("format", "children", "dictionary", "ordered", "keys_sorted")

    def __init__(self, format: str, children: Sequence["Field"] = (), dictionary: Optional["DataType"] = None,
                 ordered: bool = False, keys_sorted: bool = False) -> None:
        self.format = format
        self.children = list(children)
        self.dictionary = dictionary
        self.ordered = ordered
        self.keys_sorted = keys_sorted

    # identification ---------------------------------------------------------
    @property
    def id(self) -> str:
        """Type family: ``int64``, ``string``, ``timestamp``, ``decimal128``, ``list``, ``struct``..."""
        if self.dictionary is not None:
            return "dictionary"
        f = self.format
        if f in _SIMPLE:
            return _SIMPLE[f].split("[")[0]
        if f.startswith("d:"):
            return f"decimal{self.bit_width}"
        if f.startswith("w:"):
            return "fixed_size_binary"
        if f.startswith("ts"):
            return "timestamp"
        if f.startswith("tD"):
            return "duration"
        return {"+l": "list", "+L": "large_list", "+s": "struct", "+m": "map", "+vl": "list_view",
                "+vL": "large_list_view", "+r": "run_end_encoded"}.get(
            f, "fixed_size_list" if f.startswith("+w:") else ("union" if f.startswith("+u") else f))

    @property
    def unit(self) -> Optional[str]:
        """Time unit of a timestamp, time or duration (``"s"``, ``"ms"``, ``"us"``, ``"ns"``), ``"day"`` / ``"ms"`` for dates; None for other types."""
        f = self.format
        if f.startswith(("ts", "tD")) or f in ("tts", "ttm", "ttu", "ttn"):
            return _UNITS[f[2]]
        if f in ("tdD", "tdm"):
            return "day" if f == "tdD" else "ms"
        return None

    @property
    def tz(self) -> Optional[str]:
        """Time zone of a timestamp (``"UTC"``, ``"Europe/Paris"``, ``"+02:00"``), or None when it has none."""
        if self.format.startswith("ts") and len(self.format) > 4:
            return self.format[4:]
        return None

    @property
    def precision(self) -> Optional[int]:
        """Total number of digits of a decimal type; None for other types."""
        return int(self.format[2:].split(",")[0]) if self.format.startswith("d:") else None

    @property
    def scale(self) -> Optional[int]:
        """Digits after the decimal point of a decimal type; None for other types."""
        return int(self.format[2:].split(",")[1]) if self.format.startswith("d:") else None

    @property
    def bit_width(self) -> Optional[int]:
        """Size of one value in bits for fixed-width types (128 or 256 for decimals); None for strings, binary and nested types."""
        f = self.format
        if f.startswith("d:"):
            parts = f[2:].split(",")
            return int(parts[2]) if len(parts) > 2 else 128
        widths = {"b": 1, "c": 8, "C": 8, "s": 16, "S": 16, "i": 32, "I": 32, "l": 64, "L": 64,
                  "e": 16, "f": 32, "g": 64, "tdD": 32, "tdm": 64, "tts": 32, "ttm": 32, "ttu": 64, "ttn": 64}
        if f in widths:
            return widths[f]
        if f.startswith(("ts", "tD")):
            return 64
        if f.startswith("w:"):
            return 8 * int(f[2:])
        return None

    @property
    def byte_width(self) -> Optional[int]:
        """Size of one value in bytes for fixed-width types; None for booleans and variable-width types."""
        bw = self.bit_width
        return None if bw is None or self.format == "b" else bw // 8

    @property
    def list_size(self) -> Optional[int]:
        """Number of items in every list of a fixed-size list type; None for other types."""
        return int(self.format[3:]) if self.format.startswith("+w:") else None

    @property
    def value_field(self) -> Optional["Field"]:
        """The item :class:`Field` of a list type (name, type, nullability); None for other types."""
        return self.children[0] if self.format in ("+l", "+L") or self.format.startswith("+w:") else None

    @property
    def value_type(self) -> Optional["DataType"]:
        """Type of the items of a list, or of the values of a dictionary (categorical) type; None for other types."""
        if self.dictionary is not None:
            return self.dictionary
        vf = self.value_field
        return vf.type if vf is not None else None

    @property
    def index_type(self) -> Optional["DataType"]:
        """Integer type of the codes of a dictionary (categorical) type; None for other types."""
        return DataType(self.format) if self.dictionary is not None else None

    @property
    def fields(self) -> list:
        """The child :class:`Field` objects of a struct type (empty list for other types)."""
        return list(self.children) if self.format == "+s" else []

    @property
    def key_type(self) -> Optional["DataType"]:
        """Type of the keys of a map type; None for other types."""
        return self.children[0].type.children[0].type if self.format == "+m" else None

    @property
    def item_type(self) -> Optional["DataType"]:
        """Type of the values of a map type; None for other types."""
        return self.children[0].type.children[1].type if self.format == "+m" else None

    @property
    def is_nested(self) -> bool:
        """True for list, struct, map and other types made of child types."""
        return self.format.startswith("+")

    def __str__(self) -> str:
        if self.dictionary is not None:
            return (f"dictionary<values={self.dictionary}, indices={DataType(self.format)}, "
                    f"ordered={int(self.ordered)}>")
        f = self.format
        if f in _SIMPLE:
            return _SIMPLE[f]
        if f.startswith("d:"):
            return f"decimal{self.bit_width}({self.precision}, {self.scale})"
        if f.startswith("w:"):
            return f"fixed_size_binary[{f[2:]}]"
        if f.startswith("ts"):
            return f"timestamp[{self.unit}" + (f", tz={self.tz}]" if self.tz else "]")
        if f.startswith("tD"):
            return f"duration[{self.unit}]"
        if f in ("+l", "+L"):
            c = self.children[0]
            return f"{'large_list' if f == '+L' else 'list'}<{c.name}: {c.type}>"
        if f.startswith("+w:"):
            c = self.children[0]
            return f"fixed_size_list<{c.name}: {c.type}>[{self.list_size}]"
        if f == "+s":
            return "struct<" + ", ".join(f"{c.name}: {c.type}" for c in self.children) + ">"
        if f == "+m":
            return f"map<{self.key_type}, {self.item_type}>"
        return f

    __repr__ = __str__

    def __eq__(self, other: object) -> bool:
        other = _as_type(other)
        if not (isinstance(other, DataType) and self.format == other.format and self.dictionary == other.dictionary
                and self.ordered == other.ordered and len(self.children) == len(other.children)):
            return False
        if self.format == "+s":
            return self.children == other.children
        # list / map children: only their types count (Arrow's own equality)
        return all(a.type == b.type for a, b in zip(self.children, other.children))

    def __hash__(self) -> int:
        return hash((self.format, tuple(self.children)))


class Field:
    """One column of a :class:`Schema` (or one child of a nested type).

    Attributes: ``name``, ``type`` (a :class:`DataType`), ``nullable`` and
    ``metadata`` (a ``{bytes: bytes}`` dict).
    """

    __slots__ = ("name", "type", "nullable", "metadata")

    def __init__(self, name: str, type: DataType, nullable: bool = True, metadata: Optional[dict] = None) -> None:
        self.name = name
        self.type = type
        self.nullable = nullable
        self.metadata = metadata or {}

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Field) and hasattr(other, "__arrow_c_schema__"):
            other = _parse_capsule_field(other)
        return (isinstance(other, Field) and self.name == other.name and self.type == other.type
                and self.nullable == other.nullable)

    def __hash__(self) -> int:
        return hash((self.name, self.type))

    def __repr__(self) -> str:
        return f"{self.name}: {self.type}" + ("" if self.nullable else " not null")


def _parse_capsule_field(obj: Any) -> Field:
    c = take_from_capsule(obj.__arrow_c_schema__(), SCHEMA_CAPSULE, ArrowSchema)
    try:
        return _parse_node(c)
    finally:
        release_schema(c)


def _as_type(other: Any) -> Any:
    """A :class:`DataType` for any object exposing ``__arrow_c_schema__``
    (a PyArrow ``DataType`` for instance)."""
    if isinstance(other, DataType) or not hasattr(other, "__arrow_c_schema__"):
        return other
    return _parse_capsule_field(other).type


def _metadata(ptr: int) -> dict:
    """Decode an Arrow metadata blob (int32 count, then length-prefixed key/value pairs)."""
    if not ptr:
        return {}
    n = ctypes.c_int32.from_address(ptr).value
    pos = ptr + 4
    out = {}
    for _ in range(max(n, 0)):
        kl = ctypes.c_int32.from_address(pos).value
        key = ctypes.string_at(pos + 4, kl)
        pos += 4 + kl
        vl = ctypes.c_int32.from_address(pos).value
        val = ctypes.string_at(pos + 4, vl)
        pos += 4 + vl
        out[key] = val
    return out


def _parse_node(s: ArrowSchema) -> Field:
    fmt = s.format.decode("utf-8")
    children = [_parse_node(s.children[i].contents) for i in range(s.n_children)]
    dictionary = None
    if s.dictionary:
        dictionary = _parse_node(s.dictionary.contents).type
    t = DataType(fmt, children, dictionary, ordered=bool(s.flags & 1), keys_sorted=bool(s.flags & 4))
    name = s.name.decode("utf-8", "replace") if s.name else ""
    return Field(name, t, bool(s.flags & 2), _metadata(s.metadata))


class Schema:
    """The columns of a table: their names and types, plus key/value metadata.

    ::

        s = carquet.read_schema("data.parquet")       # or table.schema
        s.names, s.types, len(s)
        s.field("price").type, s["price"].nullable
        s.metadata                                    # {bytes: bytes}

    Iterating gives the :class:`Field` objects. Arrow libraries accept it
    where they expect a schema (``pyarrow.schema(s)``).
    """

    def __init__(self, c_schema: ArrowSchema) -> None:
        """Take ownership of a populated top-level struct ``ArrowSchema``."""
        self._c = c_schema
        root = _parse_node(c_schema)
        self.fields: list = root.type.children
        self.metadata: dict = root.metadata
        self._index = {f.name: i for i, f in reversed(list(enumerate(self.fields)))}

    def __del__(self) -> None:  # pragma: no cover - best effort
        try:
            release_schema(self._c)
        except Exception:
            pass

    @property
    def names(self) -> list:
        """Column names, in order."""
        return [f.name for f in self.fields]

    @property
    def types(self) -> list:
        """The :class:`DataType` of every column, in order."""
        return [f.type for f in self.fields]

    def __len__(self) -> int:
        return len(self.fields)

    def __iter__(self) -> Iterator[Field]:
        return iter(self.fields)

    def __getitem__(self, key: Union[int, str]) -> Field:
        return self.field(key)

    def field(self, key: Union[int, str]) -> Field:
        """The :class:`Field` at a position or with a name (``schema["price"]`` does the same). Raises ``KeyError`` for an unknown name."""
        if isinstance(key, str):
            i = self.get_field_index(key)
            if i < 0:
                raise KeyError(f"no field named {key!r} (fields: {self.names})")
            return self.fields[i]
        return self.fields[key]

    def get_field_index(self, name: str) -> int:
        """Position of the column called ``name``, or -1 when there is none."""
        return self._index.get(name, -1)

    def _copy(self) -> ArrowSchema:
        out = ArrowSchema()
        _check(_lib().carquet_arrow_schema_copy(byref(self._c), byref(out)), "schema copy")
        return out

    def _select(self, indices: Sequence[int]) -> "Schema":
        c = self._copy()
        idx = (c_int32 * max(1, len(indices)))(*indices)
        st = _lib().carquet_arrow_schema_select(byref(c), idx, len(indices))
        if st != _ffi.OK:
            release_schema(c)
            _check(st, "schema select")
        return Schema(c)

    def __arrow_c_schema__(self) -> object:
        return make_capsule(self._copy(), SCHEMA_CAPSULE)

    def __eq__(self, other: object) -> bool:
        other = _as_schema(other)
        return isinstance(other, Schema) and self.fields == other.fields

    def equals(self, other: Any, check_metadata: bool = False) -> bool:
        """Same field names, types and nullability as ``other`` (a carquet schema or any Arrow schema, PyArrow's too). ``check_metadata=True`` also compares the key/value metadata."""
        other = _as_schema(other)
        return self == other and (not check_metadata or self.metadata == other.metadata)

    def __repr__(self) -> str:
        lines = [repr(f) for f in self.fields]
        if self.metadata:
            lines.append(f"-- metadata: {sorted(k.decode('utf-8', 'replace') for k in self.metadata)}")
        return "\n".join(lines)


def _as_schema(other: Any) -> Any:
    if isinstance(other, Schema) or not hasattr(other, "__arrow_c_schema__"):
        return other
    try:
        return Schema(take_from_capsule(other.__arrow_c_schema__(), SCHEMA_CAPSULE, ArrowSchema))
    except Exception:
        return other


def _lib():
    lib = _ffi.get_lib()
    if not _ffi.has_symbol("carquet_arrow_ref_create"):
        raise CarquetError(StatusCode.NOT_IMPLEMENTED, "this libcarquet predates carquet 0.8 (no Arrow refs)")
    return lib


def _check(st: int, what: str) -> None:
    if st != _ffi.OK:
        raise CarquetError(st, f"{what} failed")


class _Ref:
    """Owner of one C ``carquet_arrow_ref_t`` (an exported array, shared)."""

    __slots__ = ("ptr", "__weakref__")

    def __init__(self, array: ArrowArray) -> None:
        ptr = _lib().carquet_arrow_ref_create(byref(array))
        if not ptr:
            if array.release:
                array.release(byref(array))
            raise MemoryError("carquet_arrow_ref_create failed")
        self.ptr = ptr

    def array(self) -> ArrowArray:
        return _lib().carquet_arrow_ref_array(self.ptr).contents

    def __del__(self) -> None:  # pragma: no cover - best effort
        try:
            _ffi.get_lib().carquet_arrow_ref_release(self.ptr)
        except Exception:
            pass


# --------------------------------------------------------------------------
# Arrays
# --------------------------------------------------------------------------
class Array:
    """The values of one column in one batch (a chunk of a :class:`ChunkedArray`).

    ::

        a = batch.column("price")
        len(a), a.type, a.null_count
        a[0], a.to_pylist(), a.to_numpy(), a.slice(10, 5)

    It is a view of buffers held by the C library: slicing copies nothing,
    and Arrow libraries take it as is (``pyarrow.array(a)``).
    """

    __slots__ = ("_ref", "_child", "_offset", "_length", "_field", "_schema")

    def __init__(self, ref: _Ref, child: Optional[int], offset: int, length: int, field: Field,
                 schema: Schema) -> None:
        self._ref = ref
        self._child = child
        self._offset = offset
        self._length = length
        self._field = field
        self._schema = schema  # the batch schema (C copy source for __arrow_c_array__)

    @property
    def type(self) -> DataType:
        """The :class:`DataType` of the values."""
        return self._field.type

    def __len__(self) -> int:
        return self._length

    def _node(self) -> tuple:
        """``(ArrowArray, absolute index of row 0 in its buffers)``."""
        root = self._ref.array()
        if self._child is None:
            return root, root.offset + self._offset
        node = root.children[self._child].contents
        return node, root.offset + self._offset + node.offset

    @property
    def null_count(self) -> int:
        """Number of null values."""
        from ._decode import null_count

        node, start = self._node()
        return null_count(node, start, self._length)

    def slice(self, offset: int = 0, length: Optional[int] = None) -> "Array":
        """Rows ``offset`` to ``offset + length`` (to the end when ``length`` is None) as a new Array. A view: nothing is copied."""
        offset = max(0, min(int(offset), self._length))
        n = self._length - offset if length is None else max(0, min(int(length), self._length - offset))
        return Array(self._ref, self._child, self._offset + offset, n, self._field, self._schema)

    def to_pylist(self) -> list:
        """The values as a Python list, None for nulls.

        Strings come back as ``str``, timestamps as ``datetime`` (time-zone aware
        when the column has a zone), dates as ``date``, decimals as ``Decimal``,
        lists as ``list``, structs as ``dict`` and maps as lists of ``(key, value)``.
        """
        from ._decode import to_pylist

        node, start = self._node()
        return to_pylist(node, self.type, start, self._length)

    def _raw(self) -> list:
        from ._decode import to_pylist

        node, start = self._node()
        return to_pylist(node, self.type, start, self._length, raw=True)

    def to_numpy(self, zero_copy_only: bool = False):
        """The values as a NumPy array (needs NumPy).

        Numbers, dates and timestamps without nulls are a view of the Arrow buffer:
        no copy, and the array keeps the data alive. With nulls the result is a
        copy where integers become float64 with NaN and dates / timestamps get NaT.
        Strings, booleans with nulls, decimals and nested types give an object
        array. ``zero_copy_only=True`` raises ``ValueError`` instead of copying.
        """
        from ._decode import to_numpy

        node, start = self._node()
        return to_numpy(node, self.type, start, self._length, keep=self._ref, zero_copy_only=zero_copy_only)

    def __getitem__(self, i: int) -> Any:
        if i < 0:
            i += self._length
        if not 0 <= i < self._length:
            raise IndexError("index out of range")
        return self.slice(i, 1).to_pylist()[0]

    def __iter__(self) -> Iterator[Any]:
        return iter(self.to_pylist())

    def _export(self) -> ArrowArray:
        out = ArrowArray()
        lib = _lib()
        if self._child is None:
            st = lib.carquet_arrow_ref_export(self._ref.ptr, self._offset, self._length, byref(out))
        else:
            st = lib.carquet_arrow_ref_export_child(self._ref.ptr, self._child, self._offset, self._length,
                                                    byref(out))
        _check(st, "array export")
        return out

    def _c_schema(self) -> ArrowSchema:
        out = ArrowSchema()
        if self._child is None:
            _check(_lib().carquet_arrow_schema_copy(byref(self._schema._c), byref(out)), "schema copy")
        else:
            _check(_lib().carquet_arrow_schema_copy_child(byref(self._schema._c), self._child, byref(out)),
                   "schema copy")
        return out

    def __arrow_c_array__(self, requested_schema: Any = None) -> tuple:
        return make_capsule(self._c_schema(), SCHEMA_CAPSULE), make_capsule(self._export(), ARRAY_CAPSULE)

    def equals(self, other: Any) -> bool:
        """Same type and values (``other``: any Arrow array, PyArrow's too)."""
        other = chunked_array(other)
        return self.type == other.type and self._raw() == other._raw()

    def __repr__(self) -> str:
        vals = self.slice(0, 10).to_pylist()
        more = ", ..." if self._length > 10 else ""
        return f"<carquet.Array {self.type} len={self._length}> [{', '.join(map(repr, vals))}{more}]"


class ChunkedArray:
    """A column of a :class:`Table`.

    ::

        col = table.column("price")
        len(col), col.type, col.null_count
        col[0], col.to_pylist(), col.to_numpy()

    It is stored as one or more :class:`Array` chunks (``col.chunks``),
    usually one per batch read from the file; the methods above treat them as
    one sequence.
    """

    def __init__(self, chunks: list, field: Field) -> None:
        self.chunks = chunks
        self._field = field

    @property
    def type(self) -> DataType:
        """The :class:`DataType` of the values."""
        return self._field.type

    @property
    def num_chunks(self) -> int:
        """Number of chunks (usually one per batch read from the file)."""
        return len(self.chunks)

    def chunk(self, i: int) -> Array:
        """The ``i``-th chunk, an :class:`Array`."""
        return self.chunks[i]

    def __len__(self) -> int:
        return sum(len(c) for c in self.chunks)

    @property
    def null_count(self) -> int:
        """Number of null values in the whole column."""
        return sum(c.null_count for c in self.chunks)

    def to_pylist(self) -> list:
        """The whole column as a Python list, None for nulls (see :meth:`Array.to_pylist` for the value types)."""
        out: list = []
        for c in self.chunks:
            out.extend(c.to_pylist())
        return out

    def to_numpy(self, zero_copy_only: bool = False):
        """The whole column as one NumPy array (needs NumPy).

        A single chunk of numbers without nulls is returned without a copy;
        several chunks are concatenated. See :meth:`Array.to_numpy` for how nulls
        and non-numeric types come out. ``zero_copy_only=True`` raises
        ``ValueError`` when a copy would be needed.
        """
        from ._decode import concat_numpy

        parts = [c.to_numpy(zero_copy_only=zero_copy_only) for c in self.chunks]
        if zero_copy_only and len(parts) > 1:
            raise ValueError("a column of several chunks cannot be exported to NumPy without a copy")
        return concat_numpy(parts, self.type)

    def __getitem__(self, i: int) -> Any:
        n = len(self)
        if i < 0:
            i += n
        for c in self.chunks:
            if i < len(c):
                return c[i]
            i -= len(c)
        raise IndexError("index out of range")

    def __iter__(self) -> Iterator[Any]:
        for c in self.chunks:
            yield from c

    def __arrow_c_stream__(self, requested_schema: Any = None) -> object:
        if not self.chunks:
            raise CarquetError(StatusCode.INVALID_STATE, "an empty column has no Arrow schema to export")
        return _list_stream(self.chunks[0]._c_schema(), [c._export() for c in self.chunks])

    def _raw(self) -> list:
        out: list = []
        for c in self.chunks:
            out.extend(c._raw())
        return out

    def equals(self, other: Any) -> bool:
        """Same type and values (``other``: any Arrow array or chunked array, PyArrow's too)."""
        other = chunked_array(other)
        return self.type == other.type and self._raw() == other._raw()

    def __repr__(self) -> str:
        return f"<carquet.ChunkedArray {self.type} len={len(self)} chunks={len(self.chunks)}>"


def _list_stream(c_schema: ArrowSchema, arrays: list) -> object:
    arr = (ArrowArray * max(1, len(arrays)))()
    for i, a in enumerate(arrays):
        ctypes.pointer(arr[i])[0] = a
        a.release = type(a.release)()
    stream = ArrowArrayStream()
    st = _lib().carquet_arrow_stream_from_arrays(byref(c_schema), arr, len(arrays), byref(stream))
    if st != _ffi.OK:
        release_schema(c_schema)
        for i in range(len(arrays)):
            if arr[i].release:
                arr[i].release(byref(arr[i]))
        _check(st, "stream export")
    return stream_capsule(stream, None)


# --------------------------------------------------------------------------
# Batches and tables
# --------------------------------------------------------------------------
class RecordBatch:
    """A batch of rows: what :meth:`ParquetFile.iter_batches` yields.

    Works like a small :class:`Table` (``num_rows``, ``column_names``,
    ``column()``, ``slice()``, ``to_pydict()``, ``to_pylist()``), with each
    column an :class:`Array`. Arrow libraries take it without a copy.
    """

    def __init__(self, schema: Schema, ref: Optional[_Ref], offset: int, num_rows: int,
                 cols: Optional[list] = None, base_schema: Optional[Schema] = None) -> None:
        self.schema = schema
        self._ref = ref
        self._offset = offset
        self.num_rows = num_rows
        self._cols = cols              # children of the held array, when a selection
        self._base_schema = base_schema or schema

    @classmethod
    def _from_c(cls, schema: Schema, array: ArrowArray) -> "RecordBatch":
        n = int(array.length)
        return cls(schema, _Ref(array), 0, n)

    @property
    def num_columns(self) -> int:
        """Number of columns."""
        return len(self.schema)

    @property
    def column_names(self) -> list:
        """Column names, in order."""
        return self.schema.names

    def column(self, key: Union[int, str]) -> Array:
        """One column, by position or name, as an :class:`Array`. Raises ``KeyError`` for an unknown name."""
        i = self.schema.get_field_index(key) if isinstance(key, str) else int(key)
        if isinstance(key, str) and i < 0:
            raise KeyError(f"no column named {key!r}")
        child = self._cols[i] if self._cols is not None else i
        return Array(self._ref, child, self._offset, self.num_rows, self.schema.field(i), self._base_schema)

    @property
    def columns(self) -> list:
        """Every column, as a list of :class:`Array`."""
        return [self.column(i) for i in range(self.num_columns)]

    def __len__(self) -> int:
        return self.num_rows

    def slice(self, offset: int = 0, length: Optional[int] = None) -> "RecordBatch":
        """Rows ``offset`` to ``offset + length`` (to the end when ``length`` is None) as a new batch. A view: nothing is copied."""
        offset = max(0, min(int(offset), self.num_rows))
        n = self.num_rows - offset if length is None else max(0, min(int(length), self.num_rows - offset))
        return RecordBatch(self.schema, self._ref, self._offset + offset, n, self._cols, self._base_schema)

    def select(self, idx: list, schema: Schema) -> "RecordBatch":
        """The columns at positions ``idx`` under ``schema`` (the schema of that selection). Used by :meth:`Table.select`, which is the easier way to pick columns."""
        cols = [self._cols[i] for i in idx] if self._cols is not None else list(idx)
        return RecordBatch(schema, self._ref, self._offset, self.num_rows, cols, self._base_schema)

    def _export(self) -> ArrowArray:
        out = ArrowArray()
        lib = _lib()
        if self._cols is None:
            st = lib.carquet_arrow_ref_export(self._ref.ptr, self._offset, self.num_rows, byref(out))
        else:
            idx = (c_int32 * max(1, len(self._cols)))(*self._cols)
            st = lib.carquet_arrow_ref_export_select(self._ref.ptr, idx, len(self._cols), self._offset,
                                                     self.num_rows, byref(out))
        _check(st, "batch export")
        return out

    def __arrow_c_array__(self, requested_schema: Any = None) -> tuple:
        return make_capsule(self.schema._copy(), SCHEMA_CAPSULE), make_capsule(self._export(), ARRAY_CAPSULE)

    def __arrow_c_stream__(self, requested_schema: Any = None) -> object:
        return _list_stream(self.schema._copy(), [self._export()])

    def to_pydict(self) -> dict:
        """``{column name: list of values}`` (see :meth:`Array.to_pylist` for the value types)."""
        return {name: self.column(i).to_pylist() for i, name in enumerate(self.column_names)}

    def to_pylist(self) -> list:
        """One ``{column name: value}`` dict per row."""
        cols = self.to_pydict()
        names = list(cols)
        return [dict(zip(names, row)) for row in zip(*cols.values())] if names else [{}] * self.num_rows

    def __repr__(self) -> str:
        return f"<carquet.RecordBatch rows={self.num_rows}>\n{self.schema!r}"


class Table:
    """What reads return: named columns of equal length.

    ::

        t = carquet.read_table("data.parquet")
        t.num_rows, t.column_names, t.schema
        t.column("price").to_pylist()        # one column as Python values
        t.column("price").to_numpy()         # ... or as a NumPy array
        t.to_pydict()                        # {name: [values]}
        t.to_pylist()                        # [{name: value}, ...] one dict per row
        t.to_pandas()                        # pandas DataFrame
        t.select(["id", "price"]).slice(0, 10)

    The data stays in the buffers the C library decoded into: ``select`` and
    ``slice`` are views, and Arrow libraries take the table without a copy
    (``polars.DataFrame(t)``, ``pyarrow.table(t)``, ``duckdb.sql("SELECT *
    FROM t")``). Build one from your own data with :func:`carquet.table`.
    """

    def __init__(self, schema: Schema, batches: list, num_rows: Optional[int] = None) -> None:
        self.schema = schema
        self._batches = [b for b in batches if b.num_rows > 0] if len(schema) else []
        self.num_rows = sum(b.num_rows for b in batches) if num_rows is None else int(num_rows)

    # shape ------------------------------------------------------------------
    @property
    def num_columns(self) -> int:
        """Number of columns."""
        return len(self.schema)

    @property
    def column_names(self) -> list:
        """Column names, in order."""
        return self.schema.names

    @property
    def shape(self) -> tuple:
        """``(num_rows, num_columns)``."""
        return (self.num_rows, self.num_columns)

    def __len__(self) -> int:
        return self.num_rows

    # access -----------------------------------------------------------------
    def column(self, key: Union[int, str]) -> ChunkedArray:
        """One column, by position or name, as a :class:`ChunkedArray` (``table["price"]`` does the same). Raises ``KeyError`` for an unknown name."""
        i = self.schema.get_field_index(key) if isinstance(key, str) else int(key)
        if isinstance(key, str) and i < 0:
            raise KeyError(f"no column named {key!r} (columns: {self.column_names})")
        return ChunkedArray([b.column(i) for b in self._batches], self.schema.field(i))

    __getitem__ = column

    @property
    def columns(self) -> list:
        """Every column, as a list of :class:`ChunkedArray`."""
        return [self.column(i) for i in range(self.num_columns)]

    def to_batches(self) -> list:
        """The :class:`RecordBatch` objects the table is made of (no copy)."""
        return list(self._batches)

    def select(self, columns: Sequence[Union[int, str]]) -> "Table":
        """A table with only these columns (names or positions), in the order given. A view: nothing is copied."""
        idx = [self.schema.get_field_index(c) if isinstance(c, str) else int(c) for c in columns]
        if any(i < 0 for i in idx):
            raise KeyError(f"unknown column in {list(columns)}")
        schema = self.schema._select(idx)
        return Table(schema, [b.select(idx, schema) for b in self._batches], self.num_rows)

    def slice(self, offset: int = 0, length: Optional[int] = None) -> "Table":
        """Rows ``offset`` to ``offset + length`` (to the end when ``length`` is None) as a new table. A view: nothing is copied."""
        offset = max(0, min(int(offset), self.num_rows))
        end = self.num_rows if length is None else min(self.num_rows, offset + max(0, int(length)))
        out, row = [], 0
        for b in self._batches:
            lo, hi = max(offset, row), min(end, row + b.num_rows)
            if hi > lo:
                out.append(b.slice(lo - row, hi - lo))
            row += b.num_rows
        return Table(self.schema, out, end - offset)

    # conversions ------------------------------------------------------------
    def to_pydict(self) -> dict:
        """``{column name: list of values}``, None for nulls (see :meth:`Array.to_pylist` for the value types)."""
        return {name: self.column(i).to_pylist() for i, name in enumerate(self.column_names)}

    def to_pylist(self) -> list:
        """One ``{column name: value}`` dict per row."""
        cols = self.to_pydict()
        names = list(cols)
        return [dict(zip(names, row)) for row in zip(*cols.values())] if names else [{}] * self.num_rows

    def to_pandas(self, **kwargs):
        """The table as a pandas DataFrame (needs pandas and NumPy).

        Numbers, dates and timestamps become NumPy-backed columns, strings
        pandas strings, dictionary columns categoricals; an integer column
        with nulls becomes float64 with NaN, as with PyArrow. A file written
        from a DataFrame gets its index, nullable dtypes and time zones back
        (from the ``pandas`` metadata stored in the file).
        """
        from ._decode import table_to_pandas

        return table_to_pandas(self, **kwargs)

    def to_polars(self):
        """A ``polars.DataFrame`` sharing the table's buffers (needs polars)."""
        import polars

        return polars.DataFrame(self)

    def to_pyarrow(self):
        """A ``pyarrow.Table`` (needs PyArrow; zero-copy)."""
        import pyarrow

        return pyarrow.table(self)

    # Arrow PyCapsule interface ------------------------------------------------
    def __arrow_c_stream__(self, requested_schema: Any = None) -> object:
        return _list_stream(self.schema._copy(), [b._export() for b in self._batches])

    def __arrow_c_schema__(self) -> object:
        return self.schema.__arrow_c_schema__()

    def equals(self, other: Any, check_metadata: bool = False) -> bool:
        """Same schema and values; ``other`` may be any Arrow table (PyArrow,
        polars...). Dictionary columns compare by value."""
        if not isinstance(other, Table):
            other = table(other)
        return (self.schema.equals(other.schema, check_metadata) and self.num_rows == other.num_rows
                and all(self.column(i)._raw() == other.column(i)._raw() for i in range(self.num_columns)))

    def __eq__(self, other: object) -> bool:
        return isinstance(other, Table) and self.equals(other)

    __hash__ = None  # type: ignore[assignment]

    def __repr__(self) -> str:
        return f"<carquet.Table rows={self.num_rows} columns={self.num_columns}>\n{self.schema!r}"


# --------------------------------------------------------------------------
# Import from any Arrow producer
# --------------------------------------------------------------------------
def table(data: Any) -> Table:
    """A :class:`Table` from

    * any object implementing the Arrow PyCapsule interface
      (``__arrow_c_stream__`` / ``__arrow_c_array__``): PyArrow and polars
      tables, DuckDB results, nanoarrow arrays... (not copied);
    * a ``{name: column}`` dict (Python lists, NumPy arrays, pandas Series,
      Arrow arrays) or a pandas DataFrame (with its ``pandas`` metadata).
    """
    if isinstance(data, Table):
        return data
    from . import _build

    if _build.is_pandas(data):
        cols, md = _build.from_pandas(data)
        return _build.table_from_columns(cols, md)
    if isinstance(data, dict):
        return _build.from_mapping(data)
    if hasattr(data, "__arrow_c_stream__"):
        stream = take_from_capsule(data.__arrow_c_stream__(), STREAM_CAPSULE, ArrowArrayStream)
        try:
            c_schema = ArrowSchema()
            if stream.get_schema(byref(stream), byref(c_schema)) != 0:
                raise CarquetError(StatusCode.INVALID_STATE, _stream_error(stream))
            schema = Schema(c_schema)
            batches = []
            while True:
                arr = ArrowArray()
                if stream.get_next(byref(stream), byref(arr)) != 0:
                    raise CarquetError(StatusCode.INVALID_STATE, _stream_error(stream))
                if not arr.release:
                    break
                batches.append(RecordBatch._from_c(schema, arr))
            return Table(schema, batches)
        finally:
            if stream.release:
                stream.release(byref(stream))
    if hasattr(data, "__arrow_c_array__"):
        sc, ac = data.__arrow_c_array__()
        c_schema = take_from_capsule(sc, SCHEMA_CAPSULE, ArrowSchema)
        arr = take_from_capsule(ac, ARRAY_CAPSULE, ArrowArray)
        schema = Schema(c_schema)
        if schema._c.format != b"+s":
            raise TypeError("an Arrow array must be a struct array (a record batch) to become a table")
        return Table(schema, [RecordBatch._from_c(schema, arr)])
    raise TypeError(f"expected an object implementing the Arrow PyCapsule interface, got {type(data).__name__}")


class StreamReader:
    """Batches pulled one at a time from any Arrow stream producer."""

    def __init__(self, data: Any) -> None:
        self._stream: Optional[ArrowArrayStream] = take_from_capsule(data.__arrow_c_stream__(), STREAM_CAPSULE,
                                                                     ArrowArrayStream)
        c_schema = ArrowSchema()
        if self._stream.get_schema(byref(self._stream), byref(c_schema)) != 0:
            msg = _stream_error(self._stream)
            self.close()
            raise CarquetError(StatusCode.INVALID_STATE, msg)
        self.schema = Schema(c_schema)

    def __iter__(self) -> "StreamReader":
        return self

    def __next__(self) -> RecordBatch:
        if self._stream is None:
            raise StopIteration
        arr = ArrowArray()
        if self._stream.get_next(byref(self._stream), byref(arr)) != 0:
            msg = _stream_error(self._stream)
            self.close()
            raise CarquetError(StatusCode.INVALID_STATE, msg)
        if not arr.release:
            self.close()
            raise StopIteration
        return RecordBatch._from_c(self.schema, arr)

    def close(self) -> None:
        s, self._stream = self._stream, None
        if s is not None and s.release:
            s.release(byref(s))

    def __del__(self) -> None:  # pragma: no cover - best effort
        try:
            self.close()
        except Exception:
            pass


def chunked_array(data: Any) -> ChunkedArray:
    """A :class:`ChunkedArray` from a carquet array or any Arrow array /
    chunked array producer (``__arrow_c_array__`` / ``__arrow_c_stream__``)."""
    if isinstance(data, ChunkedArray):
        return data
    if isinstance(data, Array):
        return ChunkedArray([data], data._field)
    if hasattr(data, "__arrow_c_stream__"):
        stream = take_from_capsule(data.__arrow_c_stream__(), STREAM_CAPSULE, ArrowArrayStream)
        try:
            c_schema = ArrowSchema()
            if stream.get_schema(byref(stream), byref(c_schema)) != 0:
                raise CarquetError(StatusCode.INVALID_STATE, _stream_error(stream))
            holder = _SchemaHolder(c_schema)
            chunks = []
            while True:
                arr = ArrowArray()
                if stream.get_next(byref(stream), byref(arr)) != 0:
                    raise CarquetError(StatusCode.INVALID_STATE, _stream_error(stream))
                if not arr.release:
                    break
                n = int(arr.length)
                chunks.append(Array(_Ref(arr), None, 0, n, holder.field, holder))
            return ChunkedArray(chunks, holder.field)
        finally:
            if stream.release:
                stream.release(byref(stream))
    if hasattr(data, "__arrow_c_array__"):
        sc, ac = data.__arrow_c_array__()
        holder = _SchemaHolder(take_from_capsule(sc, SCHEMA_CAPSULE, ArrowSchema))
        arr = take_from_capsule(ac, ARRAY_CAPSULE, ArrowArray)
        n = int(arr.length)
        return ChunkedArray([Array(_Ref(arr), None, 0, n, holder.field, holder)], holder.field)
    raise TypeError(f"expected an Arrow array, got {type(data).__name__}")


class _SchemaHolder:
    """Owner of the C schema of a single (non-struct) array."""

    def __init__(self, c_schema: ArrowSchema) -> None:
        self._c = c_schema
        self.field = _parse_node(c_schema)

    def __del__(self) -> None:  # pragma: no cover - best effort
        try:
            release_schema(self._c)
        except Exception:
            pass


def _stream_error(stream: ArrowArrayStream) -> str:
    msg = stream.get_last_error(byref(stream)) if stream.get_last_error else None
    return msg.decode("utf-8", "replace") if msg else "Arrow stream error"
