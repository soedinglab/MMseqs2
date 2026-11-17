"""Conversions between Python scalars and Parquet physical value bytes.

Used for statistics / page-index decoding (bytes -> Python value) and for
predicate pushdown and bloom-filter probes (Python value -> the column's
physical little-endian representation, as the C API expects it).
"""

from __future__ import annotations

import datetime
import decimal
import math
import struct
import uuid
from typing import Any, Optional

from .schema import ColumnSchema

_EPOCH = datetime.datetime(1970, 1, 1)
_EPOCH_DATE = datetime.date(1970, 1, 1)
_UTC = datetime.timezone.utc
_UNIT_PER_SECOND = {"MILLIS": 10**3, "MICROS": 10**6, "NANOS": 10**9}
_JULIAN_EPOCH_DAY = 2440588  # Julian day number of 1970-01-01

_STRING_LIKE = ("STRING", "ENUM", "JSON")


def _lt_name(col: ColumnSchema) -> Optional[str]:
    return col.logical_type.name if col.logical_type is not None else None


def _unsigned(col: ColumnSchema) -> bool:
    lt = col.logical_type
    return lt is not None and lt.name == "INTEGER" and lt.is_signed is False


# --------------------------------------------------------------------------
# bytes -> Python
# --------------------------------------------------------------------------
def decode_value(raw: Optional[bytes], col: ColumnSchema) -> Any:
    """Decode one physical value (statistics / page index) to a Python object.

    Falls back to the plain number or the raw bytes when the logical value
    cannot be represented (e.g. a timestamp outside ``datetime``'s range).
    """
    if raw is None:
        return None
    pt = col.physical_type
    lt = _lt_name(col)
    try:
        if pt == "BOOLEAN":
            return bool(raw[0]) if raw else None
        if pt in ("INT32", "INT64"):
            width = 4 if pt == "INT32" else 8
            if len(raw) < width:
                return raw
            fmt = ("<I" if width == 4 else "<Q") if _unsigned(col) else ("<i" if width == 4 else "<q")
            v = struct.unpack(fmt, raw[:width])[0]
            return _decode_int_logical(v, col, lt)
        if pt == "FLOAT":
            return struct.unpack("<f", raw[:4])[0] if len(raw) >= 4 else raw
        if pt == "DOUBLE":
            return struct.unpack("<d", raw[:8])[0] if len(raw) >= 8 else raw
        if pt == "INT96":
            if len(raw) < 12:
                return raw
            nanos, jd = struct.unpack("<qI", raw[:12])
            return _EPOCH + datetime.timedelta(days=jd - _JULIAN_EPOCH_DAY, microseconds=nanos // 1000)
        if pt in ("BYTE_ARRAY", "FIXED_LEN_BYTE_ARRAY"):
            if lt in _STRING_LIKE:
                return raw.decode("utf-8", "replace")
            if lt == "DECIMAL":
                unscaled = int.from_bytes(raw, "big", signed=True)
                return decimal.Decimal(unscaled).scaleb(-(col.scale or 0))
            if lt == "UUID" and len(raw) == 16:
                return uuid.UUID(bytes=bytes(raw))
            if lt == "FLOAT16" and len(raw) == 2:
                return struct.unpack("<e", raw)[0]
            return bytes(raw)
    except (OverflowError, ValueError, OSError, struct.error):
        return bytes(raw)
    return bytes(raw)


def _decode_int_logical(v: int, col: ColumnSchema, lt: Optional[str]) -> Any:
    if lt == "DECIMAL":
        return decimal.Decimal(v).scaleb(-(col.scale or 0))
    if lt == "DATE":
        try:
            return _EPOCH_DATE + datetime.timedelta(days=v)
        except OverflowError:
            return v
    if lt == "TIMESTAMP":
        unit = col.logical_type.unit or "MICROS"
        per_s = _UNIT_PER_SECOND[unit]
        try:
            ts = _EPOCH + datetime.timedelta(microseconds=v * 10**6 // per_s)
        except OverflowError:
            return v
        return ts.replace(tzinfo=_UTC) if col.logical_type.is_adjusted_to_utc else ts
    if lt == "TIME":
        unit = col.logical_type.unit or "MILLIS"
        per_s = _UNIT_PER_SECOND[unit]
        us = v * 10**6 // per_s
        if not 0 <= us < 86400 * 10**6:
            return v
        return (datetime.datetime.min + datetime.timedelta(microseconds=us)).time()
    return v


# --------------------------------------------------------------------------
# Python -> bytes
# --------------------------------------------------------------------------
def _as_int(value: Any) -> int:
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    if hasattr(value, "__index__"):
        return int(value.__index__())
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, (float, decimal.Decimal)):
        if isinstance(value, float) and math.isinf(value) or isinstance(value, decimal.Decimal) and value.is_infinite():
            raise OverflowError(f"{value!r} is out of an integer column's range")
        if isinstance(value, decimal.Decimal) and value == value.to_integral_value():
            return int(value)
        raise InexactValue(f"{value!r} is not an integer")
    raise TypeError(f"expected an integer, got {type(value).__name__}")


class InexactValue(ValueError):
    """The value has no exact representation in the column's type (pushdown would round)."""


def _rescale(amount: int, from_per_s: int, to_per_s: int, exact: bool) -> int:
    q, r = divmod(amount * to_per_s, from_per_s)
    if r and exact:
        raise InexactValue("value is finer than the column's time unit")
    return q


def _timestamp_units(value: Any, unit: str, exact: bool) -> int:
    per_s = _UNIT_PER_SECOND[unit]
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    ns = getattr(value, "value", None)  # pandas.Timestamp: exact nanoseconds since epoch (UTC)
    if isinstance(ns, int) and hasattr(value, "nanosecond"):
        return _rescale(ns, 10**9, per_s, exact)
    if _is_datetime64(value):
        import numpy as np

        return _rescale(int(value.astype("datetime64[ns]").astype(np.int64)), 10**9, per_s, exact)
    if isinstance(value, datetime.datetime):
        if value.tzinfo is not None:
            value = value.astimezone(_UTC).replace(tzinfo=None)
        d = value - _EPOCH
        us = (d.days * 86400 + d.seconds) * 10**6 + d.microseconds
        return _rescale(us, 10**6, per_s, exact)
    if isinstance(value, datetime.date):
        return (value - _EPOCH_DATE).days * 86400 * per_s
    raise TypeError(f"cannot compare a timestamp column with {type(value).__name__}")


def _is_datetime64(value: Any) -> bool:
    t = type(value)
    return t.__module__ == "numpy" and t.__name__ == "datetime64"


def _date_days(value: Any, exact: bool) -> int:
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if _is_datetime64(value):
        import numpy as np

        days = value.astype("datetime64[D]")
        if exact and days != value:
            raise InexactValue("datetime64 value is not a whole day")
        return int(days.astype(np.int64))
    if isinstance(value, datetime.datetime):
        if exact and (value.time() != datetime.time(0) or value.tzinfo is not None):
            raise InexactValue("datetime value is not a whole day")
        value = value.date()
    if isinstance(value, datetime.date):
        return (value - _EPOCH_DATE).days
    raise TypeError(f"cannot compare a date column with {type(value).__name__}")


def _time_units(value: Any, unit: str, exact: bool) -> int:
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    if isinstance(value, datetime.time):
        us = ((value.hour * 60 + value.minute) * 60 + value.second) * 10**6 + value.microsecond
        return _rescale(us, 10**6, _UNIT_PER_SECOND[unit], exact)
    raise TypeError(f"cannot compare a time column with {type(value).__name__}")


def _unscaled(value: Any, scale: int) -> int:
    if isinstance(value, decimal.Decimal):
        d = value
    elif isinstance(value, (int, str)) and not isinstance(value, bool):
        d = decimal.Decimal(value)
    elif isinstance(value, float):
        d = decimal.Decimal(repr(value))
    else:
        raise TypeError(f"cannot compare a decimal column with {type(value).__name__}")
    if not d.is_finite():
        raise InexactValue(f"{value!r} is not a finite decimal")
    scaled = d.scaleb(scale)
    if scaled != scaled.to_integral_value():
        raise InexactValue(f"{value!r} has more than {scale} fractional digits")
    return int(scaled)


def _float(value: Any, fmt: str, exact: bool) -> bytes:
    f = float(value)
    raw = struct.pack(fmt, f)
    if exact:
        back = struct.unpack(fmt, raw)[0]
        if back != f and f == f:
            raise InexactValue(f"{value!r} is not exactly representable in the column's float type")
    return raw


def encode_value(value: Any, col: ColumnSchema, exact: bool = False) -> bytes:
    """Encode ``value`` in ``col``'s physical representation.

    Raises ``TypeError`` / ``ValueError`` / ``OverflowError`` when the value
    cannot be expressed in the column's type; callers treat that as "this
    predicate cannot be pushed down". With ``exact=True`` a value that would
    have to be rounded (a float32 approximation, a timestamp finer than the
    column unit, a datetime compared with a date) raises :class:`InexactValue`
    too, since pruning on the rounded value could drop matching rows.
    """
    try:
        return _encode(value, col, exact)
    except struct.error as exc:  # out of the physical type's range
        raise OverflowError(str(exc)) from None


def _encode(value: Any, col: ColumnSchema, exact: bool) -> bytes:
    pt = col.physical_type
    lt = _lt_name(col)
    if value is None:
        raise TypeError("None cannot be compared; use 'is null' / 'is not null'")
    if pt == "BOOLEAN":
        if not isinstance(value, (bool, int)):
            raise TypeError(f"cannot compare a boolean column with {type(value).__name__}")
        return b"\x01" if value else b"\x00"
    if pt in ("INT32", "INT64"):
        if lt == "DATE":
            v = _date_days(value, exact)
        elif lt == "TIMESTAMP":
            v = _timestamp_units(value, col.logical_type.unit or "MICROS", exact)
        elif lt == "TIME":
            v = _time_units(value, col.logical_type.unit or "MILLIS", exact)
        elif lt == "DECIMAL":
            v = _unscaled(value, col.scale or 0)
        else:
            v = _as_int(value)
        width = 4 if pt == "INT32" else 8
        if _unsigned(col):
            return struct.pack("<I" if width == 4 else "<Q", v)
        return struct.pack("<i" if width == 4 else "<q", v)
    if pt == "FLOAT":
        return _float(value, "<f", exact)
    if pt == "DOUBLE":
        return _float(value, "<d", exact)
    if pt == "BYTE_ARRAY":
        if lt == "DECIMAL":
            v = _unscaled(value, col.scale or 0)
            n = max(1, (v.bit_length() + 8) // 8)
            return v.to_bytes(n, "big", signed=True)
        if isinstance(value, str):
            return value.encode("utf-8")
        if isinstance(value, (bytes, bytearray, memoryview)):
            return bytes(value)
        raise TypeError(f"cannot compare a {lt or 'binary'} column with {type(value).__name__}")
    if pt == "FIXED_LEN_BYTE_ARRAY":
        n = col.type_length or 0
        if lt == "DECIMAL":
            return _unscaled(value, col.scale or 0).to_bytes(n, "big", signed=True)
        if lt == "FLOAT16":
            return _float(value, "<e", exact)
        if lt == "UUID":
            if isinstance(value, uuid.UUID):
                return value.bytes
            if isinstance(value, str):
                return uuid.UUID(value).bytes
        if isinstance(value, str):
            value = value.encode("utf-8")
        if isinstance(value, (bytes, bytearray, memoryview)):
            b = bytes(value)
            if len(b) != n:
                raise ValueError(f"expected {n} bytes for {col.path}, got {len(b)}")
            return b
        raise TypeError(f"cannot compare a fixed-length binary column with {type(value).__name__}")
    if pt == "INT96":
        ns = _timestamp_units(value, "NANOS", exact)
        day, nanos = divmod(ns, 86400 * 10**9)
        return struct.pack("<qI", nanos, day + _JULIAN_EPOCH_DAY)
    raise TypeError(f"predicates on {pt} columns are not supported")


def is_nan(value: Any) -> bool:
    try:
        return isinstance(value, float) and value != value or (
            isinstance(value, decimal.Decimal) and value.is_nan())
    except Exception:
        return False


def is_float_column(col: ColumnSchema) -> bool:
    return col.physical_type in ("FLOAT", "DOUBLE") or (
        col.logical_type is not None and col.logical_type.name == "FLOAT16")
