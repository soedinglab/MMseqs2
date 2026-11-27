"""Read Arrow C Data arrays into Python values, NumPy arrays and pandas.

Every function takes a C ``ArrowArray`` node, its :class:`~carquet.arrow.DataType`,
``start`` (the buffer index of the first row: logical row + every offset on
the way down) and a row count. Values are decoded in bulk where Python
allows it (``memoryview.cast`` for fixed-width data); NumPy is only imported
by the NumPy and pandas conversions.
"""

from __future__ import annotations

import ctypes
import datetime as dt
import decimal
import json
import struct
from typing import Any

from .errors import CarquetError, StatusCode

_EPOCH_DATE = dt.date(1970, 1, 1)
_EPOCH = dt.datetime(1970, 1, 1)
_UTC = dt.timezone.utc
_FIXED = {"c": ("b", 1), "C": ("B", 1), "s": ("h", 2), "S": ("H", 2), "i": ("i", 4), "I": ("I", 4),
          "l": ("q", 8), "L": ("Q", 8), "e": ("e", 2), "f": ("f", 4), "g": ("d", 8)}
_PER_SECOND = {"s": 1, "m": 1000, "u": 1_000_000, "n": 1_000_000_000}


def _buf(node, i: int) -> int:
    if i >= node.n_buffers:
        return 0
    return node.buffers[i] or 0


def _bits(addr: int, start: int, n: int) -> list:
    if n <= 0:
        return []
    raw = ctypes.string_at(addr + (start >> 3), ((start & 7) + n + 7) >> 3)
    s = start & 7
    return [(raw[(j + s) >> 3] >> ((j + s) & 7)) & 1 for j in range(n)]


def _validity(node, start: int, n: int):
    """List of 0/1 per row, or None when every row is valid."""
    addr = _buf(node, 0)
    if not addr or node.null_count == 0:
        return None
    return _bits(addr, start, n)


def null_count(node, start: int, n: int) -> int:
    if node.n_buffers == 0 and n:  # null type
        return n
    v = _validity(node, start, n)
    return 0 if v is None else n - sum(v)


def _fixed(node, start: int, n: int, code: str, width: int) -> list:
    if n <= 0:
        return []
    raw = ctypes.string_at(_buf(node, 1) + start * width, n * width)
    if code == "e":  # memoryview.cast() only takes half floats from Python 3.12
        return list(struct.unpack(f"={n}e", raw))
    return memoryview(raw).cast(code).tolist()


def _offsets(node, start: int, n: int, large: bool) -> list:
    return _fixed(node, start, n + 1, "q" if large else "i", 8 if large else 4)


def _tz(name: str):
    if name in ("UTC", "utc", "Z", "+00:00"):
        return _UTC
    if name[:1] in "+-" and ":" in name:
        h, m = name[1:].split(":")
        delta = dt.timedelta(hours=int(h), minutes=int(m))
        return dt.timezone(-delta if name[0] == "-" else delta)
    try:
        from zoneinfo import ZoneInfo

        return ZoneInfo(name)
    except Exception:  # pragma: no cover - unknown zone name: keep UTC
        return _UTC


def to_pylist(node, t, start: int, n: int, raw: bool = False) -> list:  # noqa: C901 - one branch per Arrow type
    """Python values; with ``raw``, temporal values stay integer ticks (exact
    comparisons: ``datetime`` holds microseconds only)."""
    f = t.format
    if n <= 0:
        return []
    if raw and f[0] == "t" and t.dictionary is None:
        if f in ("tdD", "tts", "ttm"):
            out = _fixed(node, start, n, "i", 4)
        elif f.startswith("ti"):
            w = {"tiM": 4, "tiD": 8, "tin": 16}[f]
            data = ctypes.string_at(_buf(node, 1) + start * w, n * w)
            out = [data[i * w:(i + 1) * w] for i in range(n)]
        else:
            out = _fixed(node, start, n, "q", 8)
        valid = _validity(node, start, n)
        return out if valid is None else [v if ok else None for v, ok in zip(out, valid)]
    if t.dictionary is not None:
        idx = _fixed(node, start, n, *_FIXED[f])
        d = node.dictionary.contents
        values = to_pylist(d, t.dictionary, d.offset, d.length, raw)
        valid = _validity(node, start, n)
        if valid is not None:  # the index of a null slot is unspecified
            return [values[i] if ok else None for i, ok in zip(idx, valid)]
        out = [values[i] for i in idx]
    elif f == "n":
        return [None] * n
    elif f in _FIXED:
        out = _fixed(node, start, n, *_FIXED[f])
    elif f == "b":
        out = [bool(b) for b in _bits(_buf(node, 1), start, n)]
    elif f in ("u", "z", "U", "Z"):
        off = _offsets(node, start, n, f in ("U", "Z"))
        base = off[0]
        data = ctypes.string_at(_buf(node, 2) + base, off[-1] - base) if off[-1] > base else b""
        out = [data[off[i] - base:off[i + 1] - base] for i in range(n)]
        if f in ("u", "U"):
            out = [b.decode("utf-8") for b in out]
    elif f in ("vu", "vz"):
        out = []
        views = ctypes.string_at(_buf(node, 1) + 16 * start, 16 * n)
        valid = _validity(node, start, n)
        for i in range(n):
            if valid is not None and not valid[i]:  # the view of a null slot is unspecified
                out.append(None)
                continue
            ln = int.from_bytes(views[16 * i:16 * i + 4], "little", signed=True)
            if ln <= 12:
                b = views[16 * i + 4:16 * i + 4 + ln]
            else:
                bi = int.from_bytes(views[16 * i + 8:16 * i + 12], "little", signed=True)
                bo = int.from_bytes(views[16 * i + 12:16 * i + 16], "little", signed=True)
                b = ctypes.string_at(_buf(node, 2 + bi) + bo, ln)
            out.append(b.decode("utf-8") if f == "vu" else b)
    elif f.startswith("d:"):
        w = t.bit_width // 8
        raw = ctypes.string_at(_buf(node, 1) + start * w, n * w)
        scale = t.scale
        out = [decimal.Decimal(int.from_bytes(raw[i * w:(i + 1) * w], "little", signed=True)).scaleb(-scale)
               for i in range(n)]
    elif f.startswith("w:"):
        w = int(f[2:])
        raw = ctypes.string_at(_buf(node, 1) + start * w, n * w)
        out = [raw[i * w:(i + 1) * w] for i in range(n)]
    elif f == "tdD":
        out = [_EPOCH_DATE + dt.timedelta(days=v) for v in _fixed(node, start, n, "i", 4)]
    elif f == "tdm":
        out = [_EPOCH_DATE + dt.timedelta(milliseconds=v) for v in _fixed(node, start, n, "q", 8)]
    elif f in ("tts", "ttm", "ttu", "ttn"):
        per = _PER_SECOND[f[2]]
        vals = _fixed(node, start, n, "i" if f in ("tts", "ttm") else "q", 4 if f in ("tts", "ttm") else 8)
        out = [(dt.datetime.min + dt.timedelta(microseconds=v * 1_000_000 // per)).time() for v in vals]
    elif f.startswith("ts"):
        per = _PER_SECOND[f[2]]
        vals = _fixed(node, start, n, "q", 8)
        tz = f[4:]
        out = [_EPOCH + dt.timedelta(microseconds=v * 1_000_000 // per) if per <= 1_000_000
               else _EPOCH + dt.timedelta(microseconds=v // 1000) for v in vals]
        if tz:
            zone = _tz(tz)
            out = [v.replace(tzinfo=_UTC).astimezone(zone) for v in out]
    elif f.startswith("tD"):
        per = _PER_SECOND[f[2]]
        out = [dt.timedelta(microseconds=v * 1_000_000 // per) if per <= 1_000_000
               else dt.timedelta(microseconds=v // 1000) for v in _fixed(node, start, n, "q", 8)]
    elif f in ("+l", "+L"):
        off = _offsets(node, start, n, f == "+L")
        child = node.children[0].contents
        vals = to_pylist(child, t.children[0].type, off[0] + child.offset, off[-1] - off[0], raw)
        base = off[0]
        out = [vals[off[i] - base:off[i + 1] - base] for i in range(n)]
    elif f.startswith("+w:"):
        size = t.list_size
        child = node.children[0].contents
        vals = to_pylist(child, t.children[0].type, start * size + child.offset, n * size, raw)
        out = [vals[i * size:(i + 1) * size] for i in range(n)]
    elif f == "+s":
        cols = []
        for k, fld in enumerate(t.children):
            c = node.children[k].contents
            cols.append(to_pylist(c, fld.type, start + c.offset, n, raw))
        names = [fld.name for fld in t.children]
        out = [dict(zip(names, row)) for row in zip(*cols)] if names else [{} for _ in range(n)]
    elif f == "+m":
        off = _offsets(node, start, n, False)
        entries = node.children[0].contents
        et = t.children[0].type
        base = off[0]
        m = off[-1] - base
        ks = to_pylist(entries.children[0].contents, et.children[0].type,
                       base + entries.offset + entries.children[0].contents.offset, m, raw)
        vs = to_pylist(entries.children[1].contents, et.children[1].type,
                       base + entries.offset + entries.children[1].contents.offset, m, raw)
        pairs = list(zip(ks, vs))
        out = [pairs[off[i] - base:off[i + 1] - base] for i in range(n)]
    else:
        raise CarquetError(StatusCode.NOT_IMPLEMENTED, f"reading Arrow type {t} into Python is not supported")
    valid = _validity(node, start, n)
    if valid is not None:
        out = [v if ok else None for v, ok in zip(out, valid)]
    return out


# --------------------------------------------------------------------------
# NumPy
# --------------------------------------------------------------------------
_NP_FIXED = {"c": "i1", "C": "u1", "s": "i2", "S": "u2", "i": "i4", "I": "u4", "l": "i8", "L": "u8",
             "e": "f2", "f": "f4", "g": "f8"}


def _np_dtype(f: str):
    if f in _NP_FIXED:
        return _NP_FIXED[f]
    if f == "tdD":
        return "datetime64[D]"
    if f == "tdm":
        return "datetime64[ms]"
    if f.startswith("ts"):
        return f"datetime64[{ {'s': 's', 'm': 'ms', 'u': 'us', 'n': 'ns'}[f[2]] }]"
    if f.startswith("tD"):
        return f"timedelta64[{ {'s': 's', 'm': 'ms', 'u': 'us', 'n': 'ns'}[f[2]] }]"
    return None


class _Holder:
    """Exposes C memory to NumPy while keeping its owner alive (``array.base``)."""

    def __init__(self, addr: int, n: int, typestr: str, keep: Any) -> None:
        self.__array_interface__ = {"data": (addr, True), "shape": (n,), "typestr": typestr, "version": 3}
        self._keep = keep


def to_numpy(node, t, start: int, n: int, keep: Any = None, zero_copy_only: bool = False):
    import numpy as np

    f = t.format
    dtype = _np_dtype(f) if t.dictionary is None else None
    valid = _validity(node, start, n)
    if dtype == "datetime64[D]":
        # date32 is 32-bit days and datetime64[D] 64-bit: widen (a copy), never view
        if zero_copy_only:
            raise ValueError(f"a {t} column cannot be exported to NumPy without a copy")
        if n == 0:
            return np.empty(0, dtype=dtype)
        days = np.asarray(_Holder(_buf(node, 1) + start * 4, n, np.dtype("i4").str, keep))
        out = days.astype("i8").view(dtype)
        if valid is not None:
            out[~np.asarray(valid, dtype=bool)] = np.datetime64("NaT", "D")
        return out
    if dtype is not None:
        width = np.dtype(dtype).itemsize
        if n == 0:
            return np.empty(0, dtype=dtype)
        typestr = np.dtype(dtype).str
        view = np.asarray(_Holder(_buf(node, 1) + start * width, n, typestr, keep))
        if dtype.startswith(("datetime64", "timedelta64")):
            view = view.view(dtype)
        if valid is None:
            return view
        if zero_copy_only:
            raise ValueError(f"a {t} column with nulls cannot be exported to NumPy without a copy")
        mask = ~np.asarray(valid, dtype=bool)
        if dtype.startswith(("datetime64", "timedelta64")):
            out = view.copy()
            out[mask] = np.array("NaT", dtype=dtype)
            return out
        out = view.astype("f8" if f not in ("e", "f") else view.dtype, copy=True)
        out[mask] = np.nan
        return out
    if zero_copy_only:
        raise ValueError(f"a {t} column cannot be exported to NumPy without a copy")
    if f == "b" and valid is None and t.dictionary is None:
        return np.asarray(to_pylist(node, t, start, n), dtype=bool)
    out = np.empty(n, dtype=object)
    out[:] = to_pylist(node, t, start, n)
    return out


def concat_numpy(parts: list, t):
    import numpy as np

    if not parts:
        dtype = _np_dtype(t.format) if t.dictionary is None else None
        return np.empty(0, dtype=dtype or object)
    if len(parts) == 1:
        return parts[0]
    kinds = {p.dtype for p in parts}
    if len(kinds) > 1:  # e.g. int64 chunks next to float64 chunks with nulls
        if any(k == object for k in kinds):
            return np.concatenate([p.astype(object) for p in parts])
        return np.concatenate([p.astype("f8") for p in parts])
    return np.concatenate(parts)


# --------------------------------------------------------------------------
# pandas
# --------------------------------------------------------------------------
_PANDAS_EXT = {"Int8", "Int16", "Int32", "Int64", "UInt8", "UInt16", "UInt32", "UInt64",
               "Float32", "Float64", "boolean", "string"}


def _pandas_infers_str(pd) -> bool:
    """Whether string columns are the ``str`` dtype: always from pandas 3 (where
    reading the option is deprecated), opt-in through ``future.infer_string`` before."""
    try:
        if int(pd.__version__.split(".", 1)[0]) >= 3:
            return True
        return bool(pd.get_option("future.infer_string"))
    except Exception:  # pragma: no cover - pandas without the option
        return False


def _list_cells(np, values: list, t):
    """A list column the way PyArrow's ``to_pandas`` lays it out: one ndarray
    per row (numeric children typed, nulls as NaN), ``None`` for a null list."""
    child = t.children[0].type
    cf = child.format
    nested = cf.startswith(("+l", "+L", "+w"))
    dtype = None
    if child.dictionary is None:
        dtype = "bool" if cf == "b" else _NP_FIXED.get(cf)
    out = np.empty(len(values), dtype=object)
    for i, v in enumerate(values):
        if v is None:
            continue
        if nested:
            out[i] = _list_cells(np, v, child)
        elif dtype is not None and None not in v:
            out[i] = np.array(v, dtype=dtype)
        elif dtype is not None and dtype != "bool":
            out[i] = np.array([np.nan if x is None else x for x in v], dtype="f8")
        else:
            cell = np.empty(len(v), dtype=object)
            cell[:] = v
            out[i] = cell
    return out


def _pandas_major(pd) -> int:
    try:
        return int(str(pd.__version__).split(".", 1)[0])
    except ValueError:  # pragma: no cover - unusual version string
        return 2


def _column_to_pandas(pd, np, col, numpy_type: Any):
    t = col.type
    if t.dictionary is not None:
        cats = None
        codes_parts = []
        for chunk in col.chunks:
            node, start = chunk._node()
            d = node.dictionary.contents
            values = to_pylist(d, t.dictionary, d.offset, d.length)
            codes = _fixed(node, start, len(chunk), *_FIXED[t.format])
            valid = _validity(node, start, len(chunk))
            if cats is None:
                cats = values
            remap = {i: (cats.index(v) if v in cats else None) for i, v in enumerate(values)}
            for i, v in enumerate(values):
                if remap[i] is None:
                    cats.append(v)
                    remap[i] = len(cats) - 1
            if valid is not None:  # the code of a null slot is unspecified
                codes = [remap[c] if ok else -1 for c, ok in zip(codes, valid)]
            else:
                codes = [remap[c] for c in codes]
            codes_parts.extend(codes)
        return pd.Categorical.from_codes(codes_parts, categories=cats or [], ordered=t.ordered)
    if isinstance(numpy_type, str) and numpy_type in _PANDAS_EXT:
        return pd.array(col.to_pylist(), dtype=numpy_type)
    if t.format == "tdD":
        # PyArrow's default (date_as_object): datetime.date objects, None for nulls
        out = np.empty(len(col), dtype=object)
        out[:] = col.to_pylist()
        return out
    if t.format.startswith(("+l", "+L", "+w")):
        return _list_cells(np, col.to_pylist(), t)
    arr = col.to_numpy()
    if t.format in ("u", "U", "vu") and arr.dtype == object:
        # pandas 3 (or `future.infer_string`): string columns are `str` dtype. The
        # DataFrame constructor infers it from the values, but not for an
        # all-null column, which PyArrow also returns as `str`.
        if _pandas_infers_str(pd):
            return pd.array(arr, dtype=pd.StringDtype(na_value=np.nan))
    # pandas >= 2 holds s / ms / us / ns timestamps and PyArrow keeps the file's
    # unit there; older pandas only has nanoseconds.
    keep_unit = _pandas_major(pd) >= 2
    if t.format.startswith("ts") and t.tz:
        s = pd.Series(arr.astype("datetime64[ns]") if arr.dtype != object and not keep_unit else arr)
        return s.dt.tz_localize("UTC").dt.tz_convert(t.tz).array
    if arr.dtype.kind == "M" and arr.dtype != np.dtype("datetime64[ns]") and not keep_unit:
        try:
            return arr.astype("datetime64[ns]")
        except (OverflowError, ValueError):  # pragma: no cover - out of the ns range
            return arr
    return arr


def table_to_pandas(table, **_kwargs):
    import numpy as np
    import pandas as pd

    meta = None
    raw = table.schema.metadata.get(b"pandas")
    if raw:
        try:
            meta = json.loads(raw)
        except ValueError:
            meta = None
    numpy_types = {}
    field_names = {}
    if meta:
        for c in meta.get("columns", []):
            numpy_types[c.get("field_name", c.get("name"))] = c.get("numpy_type")
            field_names[c.get("field_name", c.get("name"))] = c.get("name")
    data = {}
    for i, name in enumerate(table.column_names):
        data[name] = _column_to_pandas(pd, np, table.column(i), numpy_types.get(name))
    df = pd.DataFrame(data, columns=table.column_names) if data else pd.DataFrame(index=range(table.num_rows))
    if meta:
        index_cols = meta.get("index_columns", [])
        named = [c for c in index_cols if isinstance(c, str) and c in df.columns]
        if named and len(named) == len(index_cols):
            df = df.set_index(named)
            df.index.names = [field_names.get(c) if not str(c).startswith("__index_level_") else None
                              for c in named]
        elif len(index_cols) == 1 and isinstance(index_cols[0], dict) and index_cols[0].get("kind") == "range":
            r = index_cols[0]
            idx = pd.RangeIndex(r.get("start", 0), r.get("stop", 0), r.get("step", 1), name=r.get("name"))
            if len(idx) == len(df):
                df.index = idx
        col_idx = meta.get("column_indexes") or []
        if col_idx and isinstance(col_idx[0], dict) and col_idx[0].get("name") is not None:
            df.columns.name = col_idx[0]["name"]
    return df
