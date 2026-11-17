"""Name tables for the C enums (types.h), and parsing of user-facing names.

Metadata objects report enum values as upper-case strings, the way
``pyarrow.parquet`` does (``"ZSTD"``, ``"RLE_DICTIONARY"``, ``"INT64"``);
options accept any case and a few common aliases.
"""

from __future__ import annotations

from typing import Union

PHYSICAL_TYPES = {
    0: "BOOLEAN",
    1: "INT32",
    2: "INT64",
    3: "INT96",
    4: "FLOAT",
    5: "DOUBLE",
    6: "BYTE_ARRAY",
    7: "FIXED_LEN_BYTE_ARRAY",
}
BOOLEAN, INT32, INT64, INT96, FLOAT, DOUBLE, BYTE_ARRAY, FIXED_LEN_BYTE_ARRAY = range(8)

REPETITIONS = {0: "REQUIRED", 1: "OPTIONAL", 2: "REPEATED"}

COMPRESSIONS = {
    0: "UNCOMPRESSED",
    1: "SNAPPY",
    2: "GZIP",
    3: "LZO",
    4: "BROTLI",
    5: "LZ4",
    6: "ZSTD",
    7: "LZ4_RAW",
}

# Accepted spellings for the ``compression`` option. "lz4" is the modern
# LZ4_RAW codec (7) that PyArrow/DuckDB/Spark write; "lz4_hadoop" is the legacy
# Hadoop-framed codec (5).
COMPRESSION_CODECS = {
    "none": 0,
    "uncompressed": 0,
    "snappy": 1,
    "gzip": 2,
    "lzo": 3,
    "brotli": 4,
    "lz4_hadoop": 5,
    "zstd": 6,
    "lz4": 7,
    "lz4_raw": 7,
}

ENCODINGS = {
    0: "PLAIN",
    2: "PLAIN_DICTIONARY",
    3: "RLE",
    4: "BIT_PACKED",
    5: "DELTA_BINARY_PACKED",
    6: "DELTA_LENGTH_BYTE_ARRAY",
    7: "DELTA_BYTE_ARRAY",
    8: "RLE_DICTIONARY",
    9: "BYTE_STREAM_SPLIT",
}
RLE_DICTIONARY = 8

ENCODING_NAMES = {
    "plain": 0,
    "rle": 3,
    "delta": 5,
    "delta_binary_packed": 5,
    "delta_length_byte_array": 6,
    "delta_byte_array": 7,
    "dictionary": 8,
    "rle_dictionary": 8,
    "plain_dictionary": 8,
    "byte_stream_split": 9,
    "bss": 9,
}

LOGICAL_TYPES = {
    0: "UNKNOWN",
    1: "STRING",
    2: "MAP",
    3: "LIST",
    4: "ENUM",
    5: "DECIMAL",
    6: "DATE",
    7: "TIME",
    8: "TIMESTAMP",
    9: "INTEGER",
    10: "NULL",
    11: "JSON",
    12: "BSON",
    13: "UUID",
    14: "FLOAT16",
    15: "VARIANT",
    16: "GEOMETRY",
    17: "GEOGRAPHY",
    18: "INTERVAL",
}

TIME_UNITS = {0: "MILLIS", 1: "MICROS", 2: "NANOS"}
TIME_UNIT_NAMES = {"ms": 0, "millis": 0, "us": 1, "micros": 1, "ns": 2, "nanos": 2}

EDGE_ALGORITHMS = {0: "SPHERICAL", 1: "VINCENTY", 2: "THOMAS", 3: "ANDOYER", 4: "KARNEY"}

BOUNDARY_ORDERS = {0: "UNORDERED", 1: "ASCENDING", 2: "DESCENDING"}


def codec_id(compression: Union[str, int, None]) -> int:
    if compression is None:
        return 0
    if isinstance(compression, bool):
        raise TypeError("compression must be a codec name, not a bool")
    if isinstance(compression, int):
        if compression not in COMPRESSIONS:
            raise ValueError(f"unknown compression codec id {compression}")
        return compression
    key = str(compression).lower()
    try:
        return COMPRESSION_CODECS[key]
    except KeyError:
        raise ValueError(
            f"unknown compression {compression!r}; expected one of {sorted(set(COMPRESSION_CODECS))}"
        ) from None


def encoding_id(encoding: Union[str, int]) -> int:
    if isinstance(encoding, int) and not isinstance(encoding, bool):
        if encoding not in ENCODINGS:
            raise ValueError(f"unknown encoding id {encoding}")
        return encoding
    key = str(encoding).lower()
    try:
        return ENCODING_NAMES[key]
    except KeyError:
        raise ValueError(f"unknown encoding {encoding!r}; expected one of {sorted(ENCODING_NAMES)}") from None


def time_unit_id(unit: str) -> int:
    try:
        return TIME_UNIT_NAMES[str(unit).lower()]
    except KeyError:
        raise ValueError(f"unknown time unit {unit!r}; expected 'ms', 'us' or 'ns'") from None
