"""Build Arrow data from Python objects, NumPy arrays and pandas DataFrames.

The C library assembles the arrays (``carquet_arrow_array_build`` copies the
buffers given here, so nothing Python owns is referenced afterwards); this
module decides the Arrow type of each column and lays its buffers out:

* lists of Python values: ``None`` is null; ``bool``, ``int``, ``float``,
  ``str``, ``bytes``, ``datetime`` (naive or aware), ``date``, ``time``,
  ``timedelta``, ``Decimal``, ``list`` / ``tuple`` (list type) and ``dict``
  (struct type, keys in first-seen order);
* NumPy arrays by dtype (``datetime64`` / ``timedelta64`` included, ``NaT``
  is null, object arrays are inferred as lists);
* pandas columns, including nullable extension dtypes, categoricals
  (dictionary type) and time zones, plus the ``pandas`` schema metadata that
  restores the index and dtypes on read (the layout PyArrow writes).
"""

from __future__ import annotations

import ctypes
import datetime as dt
import decimal
import json
import struct
from array import array
from ctypes import byref, c_int64, c_void_p
from typing import Any, Optional, Sequence

from . import _ffi
from ._arrow import ArrowArray, ArrowSchema, release_array, release_schema
from .errors import CarquetError, StatusCode

_EPOCH = dt.datetime(1970, 1, 1)
_EPOCH_UTC = dt.datetime(1970, 1, 1, tzinfo=dt.timezone.utc)
_EPOCH_DATE = dt.date(1970, 1, 1)
NULLABLE = 2


class Built:
    """A built column: its C schema and array (owned until moved into a parent)."""

    __slots__ = ("schema", "array")

    def __init__(self, schema: ArrowSchema, array: ArrowArray) -> None:
        self.schema = schema
        self.array = array

    def release(self) -> None:
        release_array(self.array)
        release_schema(self.schema)


def _lib():
    lib = _ffi.get_lib()
    if not _ffi.has_symbol("carquet_arrow_array_build"):
        raise CarquetError(StatusCode.NOT_IMPLEMENTED, "this libcarquet predates carquet 0.8")
    return lib


def metadata_blob(md: Optional[dict]) -> Optional[bytes]:
    if not md:
        return None
    out = [struct.pack("<i", len(md))]
    for k, v in md.items():
        k = k.encode("utf-8") if isinstance(k, str) else bytes(k)
        v = v.encode("utf-8") if isinstance(v, str) else bytes(v)
        out += [struct.pack("<i", len(k)), k, struct.pack("<i", len(v)), v]
    return b"".join(out)


def make(fmt: str, name: str, length: int, null_count: int, buffers: Sequence[Optional[bytes]],
         children: Sequence[Built] = (), dictionary: Optional[Built] = None, nullable: bool = True,
         flags: int = 0, metadata: Optional[dict] = None) -> Built:
    """Assemble one column from its buffers (copied by C) and built children (moved)."""
    lib = _lib()
    n = len(buffers)
    keep = [b if b is None or isinstance(b, bytes) else bytes(b) for b in buffers]
    ptrs = (c_void_p * max(1, n))(*[ctypes.cast(ctypes.c_char_p(b), c_void_p).value if b is not None else None
                                    for b in keep])
    sizes = (c_int64 * max(1, n))(*[len(b) if b is not None else 0 for b in keep])
    kids_a = (ArrowArray * max(1, len(children)))()
    kids_s = (ArrowSchema * max(1, len(children)))()
    for i, c in enumerate(children):
        ctypes.pointer(kids_a[i])[0] = c.array
        ctypes.pointer(kids_s[i])[0] = c.schema
        c.array.release = type(c.array.release)()
        c.schema.release = type(c.schema.release)()
    arr = ArrowArray()
    st = lib.carquet_arrow_array_build(length, null_count, n, ptrs, sizes, len(children), kids_a,
                                       byref(dictionary.array) if dictionary else None, byref(arr))
    if st != _ffi.OK:
        for i in range(len(children)):
            release_array(kids_a[i])
            release_schema(kids_s[i])
        raise MemoryError("carquet_arrow_array_build failed")
    sch = ArrowSchema()
    md = metadata_blob(metadata)
    st = lib.carquet_arrow_schema_build(fmt.encode(), name.encode("utf-8"), flags | (NULLABLE if nullable else 0),
                                        md, len(children), kids_s,
                                        byref(dictionary.schema) if dictionary else None, byref(sch))
    if st != _ffi.OK:
        release_array(arr)
        raise MemoryError("carquet_arrow_schema_build failed")
    return Built(sch, arr)


def _bitmap(valid: Sequence[bool]) -> bytes:
    out = bytearray((len(valid) + 7) // 8)
    for i, v in enumerate(valid):
        if v:
            out[i >> 3] |= 1 << (i & 7)
    return bytes(out)


def _validity(values: Sequence[Any]) -> tuple:
    valid = [v is not None for v in values]
    nulls = valid.count(False)
    return (_bitmap(valid) if nulls else None), nulls


# --------------------------------------------------------------------------
# Python values
# --------------------------------------------------------------------------
def _is_na(v: Any) -> bool:
    if v is None:
        return True
    if isinstance(v, float) and v != v:
        return False  # NaN is a float value, not a null (as in PyArrow from lists)
    t = type(v)
    return t.__name__ in ("NAType", "NaTType") and t.__module__.startswith("pandas")


def _kind(values: Sequence[Any]) -> str:
    kinds = set()
    for v in values:
        if v is None:
            continue
        if isinstance(v, bool):
            kinds.add("bool")
        elif isinstance(v, int):
            kinds.add("int")
        elif isinstance(v, float):
            kinds.add("float")
        elif isinstance(v, str):
            kinds.add("str")
        elif isinstance(v, (bytes, bytearray, memoryview)):
            kinds.add("bytes")
        elif isinstance(v, dt.datetime):
            kinds.add("datetime")
        elif isinstance(v, dt.date):
            kinds.add("date")
        elif isinstance(v, dt.time):
            kinds.add("time")
        elif isinstance(v, dt.timedelta):
            kinds.add("timedelta")
        elif isinstance(v, decimal.Decimal):
            kinds.add("decimal")
        elif isinstance(v, (list, tuple)):
            kinds.add("list")
        elif isinstance(v, dict):
            kinds.add("dict")
        elif hasattr(v, "dtype") and hasattr(v, "item"):  # NumPy scalar
            kinds.add(_kind([v.item()]))
        elif hasattr(v, "tolist") and hasattr(v, "dtype"):  # NumPy array as a list value
            kinds.add("list")
        else:
            raise TypeError(f"cannot convert a {type(v).__name__} value to Arrow")
    if not kinds:
        return "null"
    if kinds == {"int", "float"}:
        return "float"
    if len(kinds) > 1:
        raise TypeError(f"a column mixes values of types {sorted(kinds)}")
    return kinds.pop()


def _tz_name(tz: Any) -> str:
    if tz is dt.timezone.utc:
        return "UTC"
    key = getattr(tz, "key", None) or getattr(tz, "zone", None)
    if key:
        return str(key)
    off = tz.utcoffset(None)
    if off is None:
        return "UTC"
    mins = int(off.total_seconds() // 60)
    sign = "-" if mins < 0 else "+"
    mins = abs(mins)
    return f"{sign}{mins // 60:02d}:{mins % 60:02d}"


def _norm(v: Any) -> Any:
    if _is_na(v):
        return None
    if hasattr(v, "dtype") and getattr(v, "ndim", None) == 0 and hasattr(v, "item"):  # NumPy scalar
        return v.item()
    return v


def from_pylist(values: Sequence[Any], name: str = "", nullable: bool = True) -> Built:
    values = [_norm(v) for v in values]
    kind = _kind(values)
    n = len(values)
    valid, nulls = _validity(values)
    if kind == "null":
        return make("n", name, n, n, [], nullable=True)
    if kind == "bool":
        return make("b", name, n, nulls, [valid, _bitmap([bool(v) for v in values])], nullable=nullable)
    if kind == "int":
        big = any(v is not None and v > 2**63 - 1 for v in values)
        try:
            data = array("Q" if big else "q", [0 if v is None else v for v in values])
        except OverflowError:
            raise OverflowError(f"column {name!r}: integer out of the 64-bit range") from None
        return make("L" if big else "l", name, n, nulls, [valid, data.tobytes()], nullable=nullable)
    if kind == "float":
        data = array("d", [0.0 if v is None else float(v) for v in values])
        return make("g", name, n, nulls, [valid, data.tobytes()], nullable=nullable)
    if kind in ("str", "bytes"):
        enc = [b"" if v is None else (v.encode("utf-8") if kind == "str" else bytes(v)) for v in values]
        total = sum(len(b) for b in enc)
        large = total > 2**31 - 1
        offs = array("q" if large else "i", [0])
        acc = 0
        for b in enc:
            acc += len(b)
            offs.append(acc)
        fmt = ("U" if large else "u") if kind == "str" else ("Z" if large else "z")
        return make(fmt, name, n, nulls, [valid, offs.tobytes(), b"".join(enc)], nullable=nullable)
    if kind == "datetime":
        aware = {v.tzinfo is not None for v in values if v is not None}
        if len(aware) > 1:
            raise TypeError(f"column {name!r} mixes naive and time-zone aware datetimes")
        tz = ""
        if aware == {True}:
            tz = _tz_name(next(v.tzinfo for v in values if v is not None))
            ticks = [0 if v is None else _us(v - _EPOCH_UTC) for v in values]
        else:
            ticks = [0 if v is None else _us(v - _EPOCH) for v in values]
        return make(f"tsu:{tz}", name, n, nulls, [valid, array("q", ticks).tobytes()], nullable=nullable)
    if kind == "date":
        days = [0 if v is None else (v - _EPOCH_DATE).days for v in values]
        return make("tdD", name, n, nulls, [valid, array("i", days).tobytes()], nullable=nullable)
    if kind == "time":
        us = [0 if v is None else ((v.hour * 60 + v.minute) * 60 + v.second) * 10**6 + v.microsecond for v in values]
        return make("ttu", name, n, nulls, [valid, array("q", us).tobytes()], nullable=nullable)
    if kind == "timedelta":
        return make("tDu", name, n, nulls, [valid, array("q", [0 if v is None else _us(v) for v in values]).tobytes()],
                    nullable=nullable)
    if kind == "decimal":
        return _decimals(values, name, valid, nulls, nullable)
    if kind == "list":
        offs = array("i", [0])
        flat: list = []
        for v in values:
            if v is not None:
                flat.extend(v.tolist() if hasattr(v, "tolist") else v)
            offs.append(len(flat))
        child = from_pylist(flat, "item", True)
        return make("+l", name, n, nulls, [valid, offs.tobytes()], [child], nullable=nullable)
    if kind == "dict":
        keys: list = []
        for v in values:
            if v is not None:
                for k in v:
                    if k not in keys:
                        keys.append(k)
        kids = [from_pylist([None if v is None else v.get(k) for v in values], str(k), True) for k in keys]
        return make("+s", name, n, nulls, [valid], kids, nullable=nullable)
    raise TypeError(f"cannot convert column {name!r} to Arrow")  # pragma: no cover


def _us(td: dt.timedelta) -> int:
    return (td.days * 86400 + td.seconds) * 10**6 + td.microseconds


def _decimals(values: list, name: str, valid, nulls: int, nullable: bool) -> Built:
    scale = 0
    for v in values:
        if v is not None:
            if not v.is_finite():
                raise ValueError(f"column {name!r}: {v} cannot be stored as a decimal")
            scale = max(scale, -v.as_tuple().exponent)
    ints = [0 if v is None else int(v.scaleb(scale)) for v in values]
    precision = max([len(str(abs(i))) for i in ints] + [scale + 1, 1])
    if precision > 76:
        raise ValueError(f"column {name!r}: decimals need {precision} digits (at most 76)")
    width = 16 if precision <= 38 else 32
    data = b"".join(i.to_bytes(width, "little", signed=True) for i in ints)
    fmt = f"d:{precision},{scale}" if width == 16 else f"d:{precision},{scale},256"
    return make(fmt, name, len(values), nulls, [valid, data], nullable=nullable)


# --------------------------------------------------------------------------
# NumPy
# --------------------------------------------------------------------------
_NP_INT = {("i", 1): "c", ("u", 1): "C", ("i", 2): "s", ("u", 2): "S", ("i", 4): "i", ("u", 4): "I",
           ("i", 8): "l", ("u", 8): "L"}
_NP_UNIT = {"s": "s", "ms": "m", "us": "u", "ns": "n"}


def from_numpy(arr, name: str = "", mask=None, nullable: bool = True, nan_is_null: bool = False) -> Built:
    """``mask``: True where the value is null."""
    import numpy as np

    arr = np.asarray(arr)
    if arr.ndim != 1:
        raise ValueError(f"column {name!r} must be one-dimensional, got shape {arr.shape}")
    kind = arr.dtype.kind
    n = len(arr)
    if kind in "fc" and nan_is_null:
        nan = np.isnan(arr) if kind == "f" else np.zeros(n, bool)
        mask = nan if mask is None else (np.asarray(mask, bool) | nan)
    if kind in "Mm":
        nat = np.isnat(arr)
        mask = nat if mask is None else (np.asarray(mask, bool) | nat)
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        nulls = int(mask.sum())
        valid = np.packbits(~mask, bitorder="little").tobytes() if nulls else None
    else:
        nulls, valid = 0, None
    if kind == "b":
        return make("b", name, n, nulls, [valid, np.packbits(arr, bitorder="little").tobytes()], nullable=nullable)
    if kind in "iu":
        fmt = _NP_INT[(kind, arr.dtype.itemsize)]
        return make(fmt, name, n, nulls, [valid, np.ascontiguousarray(arr, arr.dtype.newbyteorder("=")).tobytes()],
                    nullable=nullable)
    if kind == "f":
        fmt = {2: "e", 4: "f", 8: "g"}.get(arr.dtype.itemsize)
        if fmt is None:
            arr, fmt = arr.astype("f8"), "g"
        return make(fmt, name, n, nulls, [valid, np.ascontiguousarray(arr, arr.dtype.newbyteorder("=")).tobytes()],
                    nullable=nullable)
    if kind == "M":
        unit = np.datetime_data(arr.dtype)[0]
        if unit == "D":
            return make("tdD", name, n, nulls, [valid, arr.astype("i4").tobytes()], nullable=nullable)
        if unit not in _NP_UNIT:
            arr, unit = arr.astype("datetime64[ns]"), "ns"
        return make(f"ts{_NP_UNIT[unit]}:", name, n, nulls, [valid, arr.view("i8").tobytes()], nullable=nullable)
    if kind == "m":
        unit = np.datetime_data(arr.dtype)[0]
        if unit not in _NP_UNIT:
            arr, unit = arr.astype("timedelta64[ns]"), "ns"
        return make(f"tD{_NP_UNIT[unit]}", name, n, nulls, [valid, arr.view("i8").tobytes()], nullable=nullable)
    values = arr.tolist()
    if mask is not None:
        values = [None if m else v for v, m in zip(values, mask)]
    if kind == "O":
        values = [None if (isinstance(v, float) and v != v and nan_is_null) else v for v in values]
    return from_pylist(values, name, nullable)


# --------------------------------------------------------------------------
# pandas
# --------------------------------------------------------------------------
def _pandas_column(s, name: str):
    """``(Built, pandas metadata entry)`` for one pandas Series."""
    import numpy as np
    import pandas as pd

    dtype = s.dtype
    meta_extra = None
    if isinstance(dtype, pd.CategoricalDtype):
        cats = from_pandas_values(dtype.categories, "")
        codes = s.cat.codes.to_numpy()
        fmt = {1: "c", 2: "s", 4: "i", 8: "l"}[codes.dtype.itemsize]
        mask = codes < 0
        nulls = int(mask.sum())
        valid = np.packbits(~mask, bitorder="little").tobytes() if nulls else None
        b = make(fmt, name, len(codes), nulls, [valid, np.where(mask, 0, codes).astype(codes.dtype).tobytes()],
                 dictionary=cats, flags=1 if dtype.ordered else 0)
        return b, ("categorical", str(codes.dtype), {"num_categories": len(dtype.categories),
                                                     "ordered": bool(dtype.ordered)})
    if isinstance(dtype, pd.DatetimeTZDtype):
        # the column's own unit (pandas >= 2: s / ms / us / ns), as PyArrow writes it
        values = s.dt.tz_convert("UTC").dt.tz_localize(None).to_numpy()
        if values.dtype.kind != "M" or np.datetime_data(values.dtype)[0] not in _NP_UNIT:
            values = values.astype("datetime64[ns]")
        b = from_numpy(values, name)
        unit_fmt = b.schema.format.decode()
        release_schema(b.schema)
        b.schema = make_schema(unit_fmt[:4] + str(dtype.tz), name, True)
        return b, ("datetimetz", str(values.dtype), {"timezone": str(dtype.tz)})
    if isinstance(dtype, np.dtype) and dtype.kind != "O":
        b = from_numpy(s.to_numpy(), name, nan_is_null=dtype.kind == "f")
        return b, (_pandas_type(dtype), str(dtype), meta_extra)
    if hasattr(dtype, "numpy_dtype") and getattr(dtype, "na_value", None) is pd.NA and dtype.name not in ("string",):
        # masked extension arrays: Int64, UInt8, Float64, boolean...
        mask = s.isna().to_numpy()
        values = s.to_numpy(dtype=dtype.numpy_dtype, na_value=0 if dtype.kind != "b" else False)
        b = from_numpy(values, name, mask=mask)
        return b, (_pandas_type(dtype.numpy_dtype), dtype.name, None)
    values = [None if (v is None or v is pd.NA or (isinstance(v, float) and v != v)) else v for v in s.tolist()]
    b = from_pylist(values, name)
    kind = _kind(values)
    ptype = {"str": "unicode", "bytes": "bytes", "decimal": "decimal", "date": "date", "datetime": "datetime",
             "list": "list[unknown]", "dict": "object", "null": "empty", "bool": "bool"}.get(kind, "object")
    numpy_type = str(dtype) if str(dtype) in ("string", "str") else "object"
    extra = None
    fmt = b.schema.format.decode()
    if kind == "decimal" and fmt.startswith("d:"):
        precision, scale = fmt[2:].split(",")[:2]
        extra = {"precision": int(precision), "scale": int(scale)}
    return b, (ptype, numpy_type, extra)


def make_schema(fmt: str, name: str, nullable: bool) -> ArrowSchema:
    sch = ArrowSchema()
    st = _lib().carquet_arrow_schema_build(fmt.encode(), name.encode("utf-8"), NULLABLE if nullable else 0,
                                           None, 0, None, None, byref(sch))
    if st != _ffi.OK:
        raise MemoryError("carquet_arrow_schema_build failed")
    return sch


def from_pandas_values(values, name: str) -> Built:
    import numpy as np

    arr = values.to_numpy() if hasattr(values, "to_numpy") else np.asarray(values)
    if arr.dtype.kind == "O":
        return from_pylist(arr.tolist(), name, False)
    return from_numpy(arr, name, nullable=False)


def _pandas_type(dtype) -> str:
    k = dtype.kind
    if k == "b":
        return "bool"
    if k in "iuf":
        return str(dtype)
    if k == "M":
        return "datetime"
    if k == "m":
        return "timedelta"
    return "object"


def from_pandas(df, preserve_index: Optional[bool] = None) -> tuple:
    """``(list of (name, Built), metadata dict)`` for a DataFrame, as
    ``pyarrow.Table.from_pandas`` lays it out."""
    import pandas as pd

    import carquet

    cols = []
    entries = []
    for name in df.columns:
        b, (ptype, ntype, extra) = _pandas_column(df[name], str(name))
        cols.append((str(name), b))
        entries.append({"name": str(name), "field_name": str(name), "pandas_type": ptype, "numpy_type": ntype,
                        "metadata": extra})
    index_columns: list = []
    idx = df.index
    if isinstance(idx, pd.RangeIndex) and preserve_index is not True:
        if preserve_index is None:
            index_columns.append({"kind": "range", "name": idx.name, "start": int(idx.start),
                                  "stop": int(idx.stop), "step": int(idx.step)})
    elif preserve_index is not False:
        levels = idx.levels if isinstance(idx, pd.MultiIndex) else [idx]
        for i, level in enumerate(levels):
            values = idx.get_level_values(i) if isinstance(idx, pd.MultiIndex) else idx
            field = str(level.name) if level.name is not None else f"__index_level_{i}__"
            b, (ptype, ntype, extra) = _pandas_column(pd.Series(values, index=None), field)
            cols.append((field, b))
            index_columns.append(field)
            entries.append({"name": level.name, "field_name": field, "pandas_type": ptype, "numpy_type": ntype,
                            "metadata": extra})
    meta = {
        "index_columns": index_columns,
        "column_indexes": [{"name": df.columns.name, "field_name": df.columns.name, "pandas_type": "unicode",
                            "numpy_type": "str" if str(df.columns.dtype) == "str" else "object",
                            "metadata": {"encoding": "UTF-8"}}],
        "columns": entries,
        "creator": {"library": "carquet", "version": carquet.__version__},
        "pandas_version": pd.__version__,
    }
    return cols, {"pandas": json.dumps(meta)}


def _from_pandas_module(obj: Any) -> bool:
    mod = type(obj).__module__
    return mod == "pandas" or mod.startswith("pandas.")


def is_pandas(obj: Any) -> bool:
    return _from_pandas_module(obj) and hasattr(obj, "columns") and hasattr(obj, "dtypes")


def table_from_columns(cols: list, metadata: Optional[dict] = None):
    """A :class:`carquet.Table` from ``[(name, Built)]`` of equal lengths."""
    from .arrow import RecordBatch, Schema, Table

    lengths = {b.array.length for _, b in cols}
    if len(lengths) > 1:
        for _, b in cols:
            b.release()
        raise ValueError(f"columns have different lengths: {sorted(lengths)}")
    n = lengths.pop() if lengths else 0
    top = make("+s", "", n, 0, [None], [b for _, b in cols], nullable=False, metadata=metadata)
    schema = Schema(top.schema)
    return Table(schema, [RecordBatch._from_c(schema, top.array)], n)


def from_mapping(data, schema: Any = None):
    """A table from ``{name: column}`` (lists, NumPy arrays, pandas Series,
    carquet / Arrow arrays)."""
    cols = []
    try:
        for name, col in data.items():
            cols.append((str(name), column(col, str(name))))
    except BaseException:
        for _, b in cols:
            b.release()
        raise
    return table_from_columns(cols)


def column(col: Any, name: str) -> Built:
    if is_pandas_series(col):
        return _pandas_column(col, name)[0]
    if hasattr(col, "__arrow_c_array__"):
        sc, ac = col.__arrow_c_array__()
        from ._arrow import ARRAY_CAPSULE, SCHEMA_CAPSULE, take_from_capsule

        c_schema = take_from_capsule(sc, SCHEMA_CAPSULE, ArrowSchema)
        arr = take_from_capsule(ac, ARRAY_CAPSULE, ArrowArray)
        try:
            return Built(_renamed(c_schema, name), arr)
        except BaseException:
            release_array(arr)
            raise
    if hasattr(col, "dtype") and hasattr(col, "ndim"):
        return from_numpy(col, name)
    if hasattr(col, "to_pylist"):
        return from_pylist(col.to_pylist(), name)
    return from_pylist(list(col), name)


def _renamed(c_schema: ArrowSchema, name: str) -> ArrowSchema:
    """The producer's schema node under another name: its children and
    dictionary are moved into a node built by C, the rest is released."""
    n = c_schema.n_children
    kids = (ArrowSchema * max(1, n))()
    for i in range(n):
        child = c_schema.children[i].contents
        ctypes.pointer(kids[i])[0] = child
        child.release = type(child.release)()
    dictionary = None
    if c_schema.dictionary:
        src = c_schema.dictionary.contents
        dictionary = ArrowSchema()
        ctypes.pointer(dictionary)[0] = src
        src.release = type(src.release)()
    # The metadata block holds NUL bytes: hand its address over, not a Python copy.
    md = ctypes.cast(c_schema.metadata, ctypes.c_char_p) if c_schema.metadata else None
    out = ArrowSchema()
    st = _lib().carquet_arrow_schema_build(c_schema.format, name.encode("utf-8"), c_schema.flags, md, n, kids,
                                           byref(dictionary) if dictionary is not None else None, byref(out))
    release_schema(c_schema)
    if st != _ffi.OK:
        for i in range(n):
            release_schema(kids[i])
        if dictionary is not None:
            release_schema(dictionary)
        raise MemoryError("carquet_arrow_schema_build failed")
    return out


def is_pandas_series(obj: Any) -> bool:
    return _from_pandas_module(obj) and hasattr(obj, "dtype") and hasattr(obj, "index") \
        and not hasattr(obj, "columns")
