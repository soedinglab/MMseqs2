"""Row filters: normalisation and translation to the C predicate.

``filters`` follows the ``pyarrow.parquet`` DNF convention:

* a list of ``(column, op, value)`` tuples is a conjunction (AND);
* a list of such lists is a disjunction (OR) of conjunctions.

Supported ops: ``==`` (``=``), ``!=``, ``<``, ``<=``, ``>``, ``>=``, ``in``,
``not in``, ``is null``, ``is not null`` (the last two take no value; a
``("x", "==", None)`` / ``("x", "!=", None)`` clause means the same).

The predicate runs in the C library (``carquet_predicate_t``): row groups
whose statistics rule it out are skipped, pages are pruned through the page
index for a single conjunction, and every batch is filtered exactly. The
semantics are those of PyArrow's: a comparison involving a null is false,
``not in`` keeps nulls unless ``None`` is in the set, ``in`` matches them only
when it is. A value the column cannot hold exactly (``x < 2.5`` on an integer
column, a float32 column against a float64 value, a datetime against a date
column, a timestamp finer than the column's unit...) is compared exactly:
the clause is rewritten on the neighbouring representable values.
"""

from __future__ import annotations

import ctypes
import math
import struct
from fractions import Fraction
from typing import Any, Optional

from . import _ffi
from ._values import InexactValue, encode_value, is_nan
from .schema import ColumnSchema, ParquetSchema

__all__ = ["normalize_filters"]

_OPS = {
    "==": "==", "=": "==", "eq": "==",
    "!=": "!=", "<>": "!=", "ne": "!=",
    "<": "<", "lt": "<", "<=": "<=", "le": "<=",
    ">": ">", "gt": ">", ">=": ">=", "ge": ">=",
    "in": "in", "not in": "not in", "not_in": "not in",
    "is null": "is null", "is_null": "is null",
    "is not null": "is not null", "is_not_null": "is not null",
}

# carquet_filter_op_t
_C_OP = {"==": 0, "!=": 1, "<": 2, "<=": 3, ">": 4, ">=": 5, "in": 7, "is null": 8, "is not null": 9, "not in": 10}
_PACKED = ("BOOLEAN", "INT32", "INT64", "INT96", "FLOAT", "DOUBLE", "FIXED_LEN_BYTE_ARRAY")


class Filters:
    """A normalised filter: a list of conjunctions of ``(column, op, value)``."""

    def __init__(self, dnf: list) -> None:
        self.dnf = dnf

    @property
    def columns(self) -> list:
        """Column names the predicate references, in first-use order."""
        out: list = []
        for conj in self.dnf:
            for col, _op, _v in conj:
                if col not in out:
                    out.append(col)
        return out


def normalize_filters(filters: Any) -> Optional[Filters]:
    if filters is None:
        return None
    t = type(filters)
    if t.__name__ == "Expression" and t.__module__.startswith("pyarrow"):
        raise TypeError("pyarrow.compute.Expression filters are not supported; use the DNF form "
                        "[(column, op, value), ...]")
    if isinstance(filters, tuple) and filters and isinstance(filters[0], str):
        filters = [filters]  # a single clause
    if not isinstance(filters, (list, tuple)):
        raise TypeError(f"filters must be a list of (column, op, value) tuples or a list of such lists; "
                        f"got {type(filters).__name__}")
    if not filters:
        return None
    if all(isinstance(c, tuple) for c in filters):
        groups = [list(filters)]
    elif all(isinstance(c, (list, tuple)) and (not c or isinstance(c[0], tuple)) for c in filters):
        groups = [list(c) for c in filters]
    else:
        raise TypeError("filters mixes clauses and clause lists; use [(...), ...] or [[(...)], [(...)]]")
    dnf = []
    for conj in groups:
        if not conj:
            raise ValueError("an empty conjunction in filters would match every row; remove it")
        dnf.append([_clause(c) for c in conj])
    return Filters(dnf)


def _clause(c: tuple) -> tuple:
    if len(c) == 2:
        col, op = c
        val = None
    elif len(c) == 3:
        col, op, val = c
    else:
        raise ValueError(f"filter clause must be (column, op, value) or (column, op); got {c!r}")
    if not isinstance(col, str):
        raise TypeError(f"filter column must be a str, got {col!r}")
    key = str(op).strip().lower()
    if key not in _OPS:
        raise ValueError(f"unsupported filter op {op!r}; expected one of {sorted(set(_OPS.values()))}")
    op = _OPS[key]
    if val is None and op == "==":
        op = "is null"
    elif val is None and op == "!=":
        op = "is not null"
    if op in ("in", "not in"):
        if isinstance(val, (str, bytes)) or not hasattr(val, "__iter__"):
            raise TypeError(f"'{op}' needs a collection of values, got {val!r}")
        val = list(val)
    elif op not in ("is null", "is not null") and len(c) == 2:
        raise ValueError(f"filter op {op!r} needs a value")
    return (col, op, val)


# --------------------------------------------------------------------------
# Exact rewriting of values the column cannot hold
# --------------------------------------------------------------------------
_FALSE = object()


def _rational(value: Any) -> Optional[Fraction]:
    import decimal

    if isinstance(value, bool):
        return Fraction(int(value))
    if isinstance(value, int):
        return Fraction(value)
    if isinstance(value, float) and math.isfinite(value):
        return Fraction(value)
    if isinstance(value, decimal.Decimal) and value.is_finite():
        return Fraction(value)
    if hasattr(value, "__index__"):
        return Fraction(int(value.__index__()))
    return None


def _int_range(col: ColumnSchema) -> tuple:
    width = 32 if col.physical_type == "INT32" else 64
    lt = col.logical_type
    if lt is not None and lt.name == "INTEGER" and lt.is_signed is False:
        return 0, 2**width - 1
    return -(2 ** (width - 1)), 2 ** (width - 1) - 1


def _pack_int(col: ColumnSchema, v: int) -> bytes:
    lo, _ = _int_range(col)
    fmt = ("<I" if lo == 0 else "<i") if col.physical_type == "INT32" else ("<Q" if lo == 0 else "<q")
    return struct.pack(fmt, v)


def _ticks(value: Any, col: ColumnSchema) -> Optional[Fraction]:
    """Exact position of a temporal value on the column's integer scale."""
    import datetime as dt

    lt = col.logical_type.name if col.logical_type else None
    if lt == "DATE":
        if isinstance(value, dt.datetime):
            if value.tzinfo is not None:
                value = value.astimezone(dt.timezone.utc).replace(tzinfo=None)
            d = value - dt.datetime(1970, 1, 1)
            return Fraction(d.days) + Fraction(d.seconds * 10**6 + d.microseconds, 86400 * 10**6)
        return None
    per_s = {"MILLIS": 10**3, "MICROS": 10**6, "NANOS": 10**9}.get(
        (col.logical_type.unit if col.logical_type else None) or ("NANOS" if col.physical_type == "INT96" else "MILLIS"))
    if lt == "TIME" and isinstance(value, dt.time):
        us = ((value.hour * 60 + value.minute) * 60 + value.second) * 10**6 + value.microsecond
        return Fraction(us * per_s, 10**6)
    if lt == "TIMESTAMP" or col.physical_type == "INT96":
        ns = getattr(value, "value", None)
        if isinstance(ns, int) and hasattr(value, "nanosecond"):  # pandas.Timestamp
            return Fraction(ns * per_s, 10**9)
        t = type(value)
        if t.__module__ == "numpy" and t.__name__ == "datetime64":
            import numpy as np

            return Fraction(int(value.astype("datetime64[ns]").astype(np.int64)) * per_s, 10**9)
        if isinstance(value, dt.datetime):
            if value.tzinfo is not None:
                value = value.astimezone(dt.timezone.utc).replace(tzinfo=None)
            d = value - dt.datetime(1970, 1, 1)
            return Fraction(((d.days * 86400 + d.seconds) * 10**6 + d.microseconds) * per_s, 10**6)
    return None


def _decimal_width(col: ColumnSchema) -> int:
    if col.physical_type == "INT32":
        return 4
    if col.physical_type == "INT64":
        return 8
    if col.physical_type == "FIXED_LEN_BYTE_ARRAY":
        return col.type_length or 16
    return 16 if (col.precision or 0) <= 38 else 32


def _encode_decimal(col: ColumnSchema, v: int) -> bytes:
    if col.physical_type in ("INT32", "INT64"):
        return struct.pack("<i" if col.physical_type == "INT32" else "<q", v)
    if col.physical_type == "FIXED_LEN_BYTE_ARRAY":
        return v.to_bytes(col.type_length or 16, "big", signed=True)
    return v.to_bytes(max(1, (v.bit_length() + 8) // 8), "big", signed=True)


_FLOAT_FMT = {"FLOAT": ("<f", "<I", 32), "DOUBLE": ("<d", "<Q", 64), "FLOAT16": ("<e", "<H", 16)}


def _float_neighbours(r: Fraction, fmt: str, ifmt: str, bits: int) -> tuple:
    """Largest float of the format below ``r`` and smallest above (``r`` is not representable)."""
    sign_bit = 1 << (bits - 1)

    def bits_of(x: float) -> int:
        return struct.unpack(ifmt, struct.pack(fmt, x))[0]

    def from_bits(u: int) -> float:
        return struct.unpack(fmt, struct.pack(ifmt, u))[0]

    def next_up(x: float) -> float:
        if x == math.inf:
            return x
        if x == 0.0:
            return from_bits(1)
        u = bits_of(x)
        return from_bits(u + 1 if not u & sign_bit else u - 1)

    def next_down(x: float) -> float:
        return -next_up(-x)

    try:
        f = struct.unpack(fmt, struct.pack(fmt, float(r)))[0]
    except (OverflowError, struct.error):
        f = math.inf if r > 0 else -math.inf
    if math.isinf(f):
        return (next_down(f), f) if f > 0 else (f, next_up(f))
    return (f, next_up(f)) if Fraction(f) < r else (next_down(f), f)


def _bracket(col: ColumnSchema, value: Any):
    """``(lo, hi)``: encodings of the closest column values below and above
    ``value`` (``None`` past the column's range), for a value the column
    cannot hold exactly. Raises ``TypeError`` when the types do not compare."""
    lt = col.logical_type.name if col.logical_type else None
    pt = col.physical_type
    if lt == "DECIMAL":
        r = _rational(value)
        if r is None and isinstance(value, str):
            import decimal

            r = _rational(decimal.Decimal(value))
        if r is None:
            raise TypeError(f"cannot compare a decimal column with {type(value).__name__}")
        r *= 10 ** (col.scale or 0)
        w = _decimal_width(col)
        mn, mx = -(2 ** (8 * w - 1)), 2 ** (8 * w - 1) - 1
        lo_i, hi_i = math.floor(r), math.ceil(r)
        lo = _encode_decimal(col, min(lo_i, mx)) if lo_i >= mn else None
        hi = _encode_decimal(col, max(hi_i, mn)) if hi_i <= mx else None
        return lo, hi
    if lt in ("DATE", "TIME", "TIMESTAMP") or pt == "INT96":
        r = _ticks(value, col)
        if r is None:
            raise TypeError(f"cannot compare a {lt or pt} column with {type(value).__name__}")
        if pt == "INT96":
            return encode_value(math.floor(r), col), encode_value(math.ceil(r), col)
        mn, mx = _int_range(col)
        lo_i, hi_i = math.floor(r), math.ceil(r)
        return (_pack_int(col, min(lo_i, mx)) if lo_i >= mn else None,
                _pack_int(col, max(hi_i, mn)) if hi_i <= mx else None)
    if pt in ("INT32", "INT64"):
        if isinstance(value, float) and math.isinf(value):
            mn, mx = _int_range(col)
            return (_pack_int(col, mx), None) if value > 0 else (None, _pack_int(col, mn))
        r = _rational(value)
        if r is None:
            raise TypeError(f"cannot compare an integer column with {type(value).__name__}")
        mn, mx = _int_range(col)
        lo_i, hi_i = math.floor(r), math.ceil(r)
        return (_pack_int(col, min(lo_i, mx)) if lo_i >= mn else None,
                _pack_int(col, max(hi_i, mn)) if hi_i <= mx else None)
    key = "FLOAT16" if lt == "FLOAT16" else pt
    if key in _FLOAT_FMT:
        r = _rational(value)
        if r is None:
            raise TypeError(f"cannot compare a floating-point column with {type(value).__name__}")
        fmt, ifmt, bits = _FLOAT_FMT[key]
        lo, hi = _float_neighbours(r, fmt, ifmt, bits)
        return struct.pack(fmt, lo), struct.pack(fmt, hi)
    raise TypeError(f"cannot compare {col.path} with {value!r}")


def _encode_exact(col: ColumnSchema, value: Any):
    """The encoded value, or ``(lo, hi)`` from :func:`_bracket` when inexact."""
    try:
        return encode_value(value, col, exact=True), None
    except (InexactValue, OverflowError):
        return None, _bracket(col, value)
    except ValueError:
        if col.physical_type in ("BYTE_ARRAY", "FIXED_LEN_BYTE_ARRAY", "BOOLEAN") and \
                not (col.logical_type and col.logical_type.name == "DECIMAL"):
            raise
        return None, _bracket(col, value)


def _variants(col: ColumnSchema, op: str, val: Any) -> list:
    """The clause as alternatives (OR) of exact C clauses ``(leaf, op, payload)``."""
    leaf = col.index
    if op in ("is null", "is not null"):
        return [[(leaf, op, None)]]
    if op in ("in", "not in"):
        has_null = any(v is None for v in val)
        encoded = []
        for v in val:
            if v is None:
                continue
            if is_nan(v) and col.physical_type not in ("FLOAT", "DOUBLE") and not (
                    col.logical_type and col.logical_type.name == "FLOAT16"):
                continue  # NaN matches nothing but a float NaN
            enc, br = _encode_exact(col, v)
            if enc is not None:
                encoded.append(enc)
        alts = []
        if op == "in":
            if encoded:
                alts.append([(leaf, "in", encoded)])
            if has_null:
                alts.append([(leaf, "is null", None)])
        else:
            base = [(leaf, "not in", encoded)] if encoded else [(leaf, "is not null", None)]
            alts.append(base)
            if not has_null:
                alts.append([(leaf, "is null", None)])
        return alts
    if is_nan(val):
        # IEEE: NaN equals nothing; every non-null value differs from it.
        return [[(leaf, "is not null", None)]] if op == "!=" else []
    enc, br = _encode_exact(col, val)
    if enc is not None:
        return [[(leaf, op, enc)]]
    lo, hi = br
    if op == "==":
        return []
    if op == "!=":
        return [[(leaf, "is not null", None)]]
    if op in ("<", "<="):
        return [[(leaf, "<=", lo)]] if lo is not None else []
    return [[(leaf, ">=", hi)]] if hi is not None else []


class Predicate:
    """A ``carquet_predicate_t`` and the buffers it points to."""

    def __init__(self, conjunctions: list, byte_array_leaves: set) -> None:
        clauses = [c for conj in conjunctions for c in conj]
        self._keep: list = []
        self.columns = sorted({c[0] for c in clauses})
        arr = (_ffi.FilterClauseStruct * max(1, len(clauses)))()
        self._keep.append(arr)
        for i, (leaf, op, payload) in enumerate(clauses):
            c = arr[i]
            c.column_index = leaf
            c.op = _C_OP[op]
            if op in ("in", "not in"):
                if leaf in byte_array_leaves:
                    vals = (_ffi.ByteArrayStruct * len(payload))()
                    for k, e in enumerate(payload):
                        b = ctypes.create_string_buffer(e, max(1, len(e)))
                        self._keep.append(b)
                        vals[k].data = ctypes.addressof(b)
                        vals[k].length = len(e)
                    self._keep.append(vals)
                    c.values = ctypes.addressof(vals)
                else:
                    blob = ctypes.create_string_buffer(b"".join(payload), max(1, sum(len(e) for e in payload)))
                    self._keep.append(blob)
                    c.values = ctypes.addressof(blob)
                c.value_count = len(payload)
            elif payload is not None:
                b = ctypes.create_string_buffer(payload, max(1, len(payload)))
                self._keep.append(b)
                c.value = ctypes.addressof(b)
                c.value_size = len(payload)
        sizes = (ctypes.c_int32 * max(1, len(conjunctions)))(*[len(c) for c in conjunctions])
        self._keep.append(sizes)
        self.struct = _ffi.PredicateStruct()
        self.struct.clauses = arr
        self.struct.conjunction_sizes = sizes
        self.struct.num_conjunctions = len(conjunctions)


def compile_filters(schema: ParquetSchema, flt: Filters, skip_repeated: bool = False):
    """``Predicate`` for ``flt``, or ``None`` when no row can match.

    With ``skip_repeated``, clauses on columns inside a list or map are dropped
    (treated as true), for statistics-only uses."""
    conjunctions: list = []
    byte_array = set()
    for conj in flt.dnf:
        alts: list = [[]]
        dropped_all = False
        for name, op, val in conj:
            try:
                col = _column(schema, name)
            except ValueError:
                if skip_repeated:  # a group: statistics cannot tell
                    continue
                raise
            if col.max_repetition_level > 0:
                if skip_repeated:
                    continue
                raise ValueError(f"cannot filter on {name!r}: it is inside a list or map")
            if col.physical_type == "BYTE_ARRAY":
                byte_array.add(col.index)
            v = _variants(col, op, val)
            if not v:
                dropped_all = True
                break
            alts = [a + b for a in alts for b in v]
        if dropped_all:
            continue
        for a in alts:
            if a:
                conjunctions.append(a)
            else:  # every clause was skipped: this conjunction matches everything
                return _ALL
    if not conjunctions:
        return None
    return Predicate(conjunctions, byte_array)


_ALL = object()  # a predicate that keeps every row


def _column(schema: ParquetSchema, name: str) -> ColumnSchema:
    try:
        return schema.column(name)
    except KeyError:
        leaves = schema.leaves_of(name) if hasattr(schema, "leaves_of") else []
        if leaves:
            raise ValueError(f"cannot filter on {name!r}: it is a group, name one of its leaf columns") from None
        raise KeyError(f"filter column {name!r} not found") from None
