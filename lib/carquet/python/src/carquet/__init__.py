"""carquet -- read and write Apache Parquet files, fast, with no dependencies.

The package wraps the carquet C library. It needs nothing else: no PyArrow,
no NumPy (pandas, NumPy, polars and PyArrow are used only when you ask for
their objects).

Read a file::

    import carquet as cq

    t = cq.read_table("data.parquet")                      # a carquet.Table
    t = cq.read_table("data.parquet", columns=["id", "price"],
                      filters=[("price", ">", 100)])       # only what you need
    t.num_rows, t.column_names, t.schema

Use the data::

    t.to_pydict()                    # {"id": [...], "price": [...]}
    t.to_pylist()                    # [{"id": 1, "price": 120.0}, ...]
    t.column("price").to_numpy()     # NumPy array
    t.to_pandas()                    # pandas DataFrame (or cq.read_pandas(path))
    polars.DataFrame(t), pyarrow.table(t), duckdb.sql("SELECT * FROM t")

Write a file::

    cq.write_table({"id": [1, 2, 3], "name": ["a", None, "c"]}, "out.parquet")
    cq.write_table(df, "out.parquet", compression="snappy")   # pandas, polars, PyArrow...

Large files, a batch at a time::

    with cq.ParquetFile("data.parquet") as pf:
        print(pf.num_rows, pf.schema)
        for batch in pf.iter_batches(batch_size=50_000, columns=["id"]):
            ...

    with cq.ParquetWriter("out.parquet", t.schema) as w:
        for batch in batches:
            w.write(batch)

Look inside a file without reading its data::

    cq.info("data.parquet")                  # rows, row groups, size
    cq.read_schema("data.parquet")           # column names and types
    cq.read_metadata("data.parquet")         # row groups, statistics, encodings

Where to look next: ``help(cq.read_table)``, ``help(cq.write_table)``,
``help(cq.Table)``, ``help(cq.ParquetFile)``, ``help(cq.ParquetWriter)``.
Errors from the C library are raised as :class:`CarquetError`.
"""

from __future__ import annotations

from . import _ffi
from ._types import COMPRESSION_CODECS
from .arrow import Array, ChunkedArray, DataType, Field, RecordBatch, Schema, Table, table
from .errors import CarquetError, StatusCode
from .metadata import (
    BloomFilter,
    ColumnChunkMetaData,
    ColumnIndex,
    FileMetaData,
    GeospatialStatistics,
    OffsetIndex,
    PageLocation,
    PageStatistics,
    RowGroupMetaData,
    Statistics,
)
from .reader import (
    ArrowStream,
    FileInfo,
    ParquetFile,
    info,
    iter_batches,
    read_metadata,
    read_pandas,
    read_parquet_schema,
    read_schema,
    read_table,
    validate,
)
from .schema import ColumnSchema, LogicalType, ParquetSchema, SchemaElement
from .writer import ParquetWriter, write_table

__all__ = [
    "COMPRESSION_CODECS",
    "Array",
    "ArrowStream",
    "ChunkedArray",
    "DataType",
    "Field",
    "RecordBatch",
    "Schema",
    "Table",
    "table",
    "BloomFilter",
    "CarquetError",
    "ColumnChunkMetaData",
    "ColumnIndex",
    "ColumnSchema",
    "FileInfo",
    "FileMetaData",
    "GeospatialStatistics",
    "LogicalType",
    "OffsetIndex",
    "PageLocation",
    "PageStatistics",
    "ParquetFile",
    "ParquetSchema",
    "ParquetWriter",
    "Reader",
    "RowGroupMetaData",
    "SchemaElement",
    "Statistics",
    "StatusCode",
    "__version__",
    "cpu_info",
    "info",
    "iter_batches",
    "lib_version",
    "open_reader",
    "read_metadata",
    "read_pandas",
    "read_parquet_schema",
    "read_schema",
    "read_table",
    "schema",
    "validate",
    "write_table",
]


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


def cpu_info() -> dict:
    """SIMD features detected at runtime (which optimised kernels the library uses)."""
    c = _ffi.get_lib().carquet_get_cpu_info().contents
    return {name: (int(getattr(c, name)) if name == "sve_vector_length" else bool(getattr(c, name)))
            for name, _ in c._fields_}


# Names from carquet 0.7.2 and earlier, kept as aliases.
Reader = ParquetFile
schema = read_schema


def open_reader(path, **kwargs) -> ParquetFile:
    """Alias of :class:`ParquetFile` (``use_mmap=`` is accepted for ``memory_map=``)."""
    return ParquetFile(path, **kwargs)
