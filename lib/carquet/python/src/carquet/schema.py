"""Parquet-level schema description (physical / logical types, levels).

:class:`ParquetSchema` is a snapshot taken from the C schema when a file is
opened: it stays valid after the file is closed and needs no PyArrow. The Arrow
view of the same schema is :attr:`carquet.ParquetFile.schema`.
"""

from __future__ import annotations

import dataclasses
from ctypes import c_char_p
from typing import Iterator, Optional, Union

from . import _ffi
from ._types import (
    EDGE_ALGORITHMS,
    LOGICAL_TYPES,
    PHYSICAL_TYPES,
    REPETITIONS,
    TIME_UNITS,
)

__all__ = ["ColumnSchema", "LogicalType", "ParquetSchema", "SchemaElement"]

_MAX_PATH_DEPTH = 256


@dataclasses.dataclass(frozen=True)
class LogicalType:
    """A Parquet logical type annotation (``carquet_logical_type_t``).

    ``name`` is the annotation (``"STRING"``, ``"TIMESTAMP"``, ``"DECIMAL"``, ...);
    the remaining attributes are set only for the annotations that carry them.
    ``str()`` gives a compact form such as ``TIMESTAMP(MICROS, UTC)``.
    """

    name: str
    precision: Optional[int] = None  # DECIMAL
    scale: Optional[int] = None  # DECIMAL
    bit_width: Optional[int] = None  # INTEGER
    is_signed: Optional[bool] = None  # INTEGER
    unit: Optional[str] = None  # TIME / TIMESTAMP: "MILLIS" | "MICROS" | "NANOS"
    is_adjusted_to_utc: Optional[bool] = None  # TIME / TIMESTAMP
    crs: Optional[str] = None  # GEOMETRY / GEOGRAPHY
    algorithm: Optional[str] = None  # GEOGRAPHY edge interpolation
    specification_version: Optional[int] = None  # VARIANT

    def __str__(self) -> str:
        n = self.name
        if n == "DECIMAL":
            return f"DECIMAL({self.precision}, {self.scale})"
        if n == "INTEGER":
            return f"INTEGER({self.bit_width}, {'signed' if self.is_signed else 'unsigned'})"
        if n in ("TIME", "TIMESTAMP"):
            return f"{n}({self.unit}{', UTC' if self.is_adjusted_to_utc else ''})"
        if n in ("GEOMETRY", "GEOGRAPHY"):
            parts = [p for p in (self.crs, self.algorithm if n == "GEOGRAPHY" else None) if p]
            return f"{n}({', '.join(parts)})" if parts else n
        return n

    @classmethod
    def _from_c(cls, lt) -> Optional["LogicalType"]:
        name = LOGICAL_TYPES.get(int(lt.id), f"UNKNOWN_{int(lt.id)}")
        p = lt.params
        if name == "UNKNOWN":
            return None
        if name == "DECIMAL":
            return cls(name, precision=int(p.decimal.precision), scale=int(p.decimal.scale))
        if name == "INTEGER":
            return cls(name, bit_width=int(p.integer.bit_width), is_signed=bool(p.integer.is_signed))
        if name in ("TIME", "TIMESTAMP"):
            t = p.time if name == "TIME" else p.timestamp
            return cls(name, unit=TIME_UNITS.get(int(t.unit), "MICROS"), is_adjusted_to_utc=bool(t.is_adjusted_to_utc))
        if name == "GEOMETRY":
            return cls(name, crs=_ffi.decode(p.geometry.crs) or None)
        if name == "GEOGRAPHY":
            g = p.geography
            algo = EDGE_ALGORITHMS.get(int(g.algorithm)) if g.has_algorithm else None
            return cls(name, crs=_ffi.decode(g.crs) or None, algorithm=algo)
        if name == "VARIANT":
            return cls(name, specification_version=int(p.variant.specification_version) or 1)
        return cls(name)


@dataclasses.dataclass(frozen=True)
class SchemaElement:
    """One node of the Parquet schema tree (a group or a leaf column)."""

    name: str
    path: tuple  # names from the root (excluded) down to this node
    is_leaf: bool
    repetition: str  # "REQUIRED" | "OPTIONAL" | "REPEATED"
    logical_type: Optional[LogicalType]
    physical_type: Optional[str]  # leaves only
    type_length: Optional[int]  # FIXED_LEN_BYTE_ARRAY leaves only

    @property
    def depth(self) -> int:
        """Nesting depth: 1 for a top-level field."""
        return len(self.path)


@dataclasses.dataclass(frozen=True)
class ColumnSchema:
    """A leaf column of a Parquet schema.

    ``index`` is the leaf column ordinal used throughout the file metadata
    (``RowGroupMetaData.column(index)``); ``path`` is the dotted path from the
    root (``"address.city"``, ``"tags.list.element"``).
    """

    index: int
    name: str
    parts: tuple
    physical_type: str
    logical_type: Optional[LogicalType]
    repetition: str
    max_definition_level: int
    max_repetition_level: int
    type_length: Optional[int]

    @property
    def path(self) -> str:
        """Dotted path of the column (``"address.city"``)."""
        return ".".join(self.parts)

    @property
    def nullable(self) -> bool:
        """True when a value of this leaf can be null (some level is optional)."""
        return self.max_definition_level > self.max_repetition_level

    @property
    def is_nested(self) -> bool:
        """True when the column sits inside a struct, list or map."""
        return len(self.parts) > 1

    @property
    def precision(self) -> Optional[int]:
        """Total digits of a DECIMAL column; None otherwise."""
        return self.logical_type.precision if self.logical_type is not None else None

    @property
    def scale(self) -> Optional[int]:
        """Digits after the decimal point of a DECIMAL column; None otherwise."""
        return self.logical_type.scale if self.logical_type is not None else None

    def __str__(self) -> str:
        lt = f" {self.logical_type}" if self.logical_type is not None else ""
        tl = f"({self.type_length})" if self.type_length else ""
        return f"{self.path}: {self.repetition} {self.physical_type}{tl}{lt}"


class ParquetSchema:
    """The Parquet schema of a file: leaf columns plus the group structure.

    Behaves like a read-only sequence of :class:`ColumnSchema` (one per leaf
    column); index it by position, by dotted path or by unique leaf name::

        s = pf.parquet_schema
        len(s), s.names            # leaf count, dotted leaf paths
        s["id"].physical_type      # 'INT64'
        s.column("address.city")
        s.fields                   # top-level field names
    """

    def __init__(self, columns: list, elements: list) -> None:
        self._columns = columns
        self._elements = elements
        self._by_path = {c.path: c for c in columns}
        by_name: dict = {}
        for c in columns:
            by_name.setdefault(c.name, []).append(c)
        self._by_name = by_name

    # -- construction from the C schema ------------------------------------
    @classmethod
    def _from_handle(cls, handle: int) -> "ParquetSchema":
        lib = _ffi.get_lib()
        ncols = int(lib.carquet_schema_num_columns(handle))
        paths = []
        buf = (c_char_p * _MAX_PATH_DEPTH)()
        for i in range(ncols):
            n = int(lib.carquet_schema_column_path(handle, i, buf, _MAX_PATH_DEPTH))
            parts = tuple(_ffi.decode(buf[k]) for k in range(n))
            if not parts:
                parts = (_ffi.decode(lib.carquet_schema_column_name(handle, i)),)
            paths.append(parts)

        nelem = int(lib.carquet_schema_num_elements(handle))
        raw = []  # (name, is_leaf, repetition, logical, physical, type_length)
        for e in range(1, nelem):
            node = lib.carquet_schema_get_element(handle, e)
            if not node:
                continue
            is_leaf = bool(lib.carquet_schema_node_is_leaf(node))
            lt_ptr = lib.carquet_schema_node_logical_type(node)
            lt = LogicalType._from_c(lt_ptr.contents) if lt_ptr else None
            phys = tl = None
            if is_leaf:
                phys = PHYSICAL_TYPES.get(int(lib.carquet_schema_node_physical_type(node)), "UNKNOWN")
                tl = int(lib.carquet_schema_node_type_length(node)) or None
            raw.append(
                (
                    _ffi.decode(lib.carquet_schema_node_name(node)),
                    is_leaf,
                    REPETITIONS.get(int(lib.carquet_schema_node_repetition(node)), "OPTIONAL"),
                    lt,
                    phys,
                    tl,
                )
            )

        # Elements are stored in pre-order; leaves appear in leaf-column order.
        # The groups seen since the previous leaf are the innermost ancestors of
        # the next leaf, which fixes their depth (and so their path).
        elements: list = []
        columns: list = []
        pending: list = []
        leaf = 0
        for name, is_leaf, rep, lt, phys, tl in raw:
            if not is_leaf:
                pending.append((name, rep, lt))
                continue
            parts = paths[leaf] if leaf < len(paths) else (name,)
            first = len(parts) - 1 - len(pending)
            for j, (gname, grep, glt) in enumerate(pending):
                d = first + j + 1
                gpath = parts[:d] if d >= 1 and parts[d - 1] == gname else (gname,)
                elements.append(SchemaElement(gname, gpath, False, grep, glt, None, None))
            pending = []
            elements.append(SchemaElement(name, parts, True, rep, lt, phys, tl))
            columns.append(
                ColumnSchema(
                    index=leaf,
                    name=name,
                    parts=parts,
                    physical_type=phys or "UNKNOWN",
                    logical_type=lt,
                    repetition=rep,
                    max_definition_level=int(lib.carquet_schema_max_def_level(handle, leaf)),
                    max_repetition_level=int(lib.carquet_schema_max_rep_level(handle, leaf)),
                    type_length=tl,
                )
            )
            leaf += 1
        for gname, grep, glt in pending:  # trailing groups without leaves
            elements.append(SchemaElement(gname, (gname,), False, grep, glt, None, None))
        return cls(columns, elements)

    # -- sequence protocol --------------------------------------------------
    def __len__(self) -> int:
        return len(self._columns)

    def __iter__(self) -> Iterator[ColumnSchema]:
        return iter(self._columns)

    def __getitem__(self, key: Union[int, str]) -> ColumnSchema:
        return self.column(key)

    def __contains__(self, key: object) -> bool:
        if not isinstance(key, str):
            return False
        return key in self._by_path or len(self._by_name.get(key, ())) == 1

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ParquetSchema):
            return NotImplemented
        return self._columns == other._columns and self._elements == other._elements

    def __hash__(self) -> int:
        return hash(tuple(self._columns))

    def column(self, key: Union[int, str]) -> ColumnSchema:
        """Leaf column by ordinal, dotted path, or (unique) leaf name."""
        if isinstance(key, int) and not isinstance(key, bool):
            if not -len(self._columns) <= key < len(self._columns):
                raise IndexError(f"column {key} out of range (schema has {len(self._columns)} leaf columns)")
            return self._columns[key]
        if not isinstance(key, str):
            raise TypeError(f"column key must be int or str, got {type(key).__name__}")
        c = self._by_path.get(key)
        if c is not None:
            return c
        hits = self._by_name.get(key, [])
        if len(hits) == 1:
            return hits[0]
        if hits:
            raise KeyError(f"leaf name {key!r} is ambiguous; use one of {[h.path for h in hits]}")
        raise KeyError(f"no column {key!r} in schema (columns: {self.names})")

    def index(self, key: str) -> int:
        """Leaf ordinal of the column at dotted path / unique leaf name ``key``."""
        return self.column(key).index

    # -- views --------------------------------------------------------------
    @property
    def columns(self) -> list:
        """Every leaf column, in file order."""
        return list(self._columns)

    @property
    def names(self) -> list:
        """Dotted paths of the leaf columns, in file order."""
        return [c.path for c in self._columns]

    @property
    def elements(self) -> list:
        """Every schema node (groups and leaves) in pre-order, root excluded."""
        return list(self._elements)

    @property
    def fields(self) -> list:
        """Names of the top-level fields (the root's children), in order."""
        return [e.name for e in self._elements if e.depth == 1]

    def leaves_of(self, field: str) -> list:
        """Leaf columns under top-level field (or group path) ``field``."""
        parts = tuple(field.split("."))
        return [c for c in self._columns if c.parts[: len(parts)] == parts]

    def __repr__(self) -> str:
        lines = ["<carquet.ParquetSchema>"]
        for e in self._elements:
            pad = "  " * e.depth
            lt = f" {e.logical_type}" if e.logical_type is not None else ""
            if e.is_leaf:
                tl = f"({e.type_length})" if e.type_length else ""
                lines.append(f"{pad}{e.name}: {e.repetition} {e.physical_type}{tl}{lt}")
            else:
                lines.append(f"{pad}{e.name}: {e.repetition} group{lt}")
        return "\n".join(lines)

