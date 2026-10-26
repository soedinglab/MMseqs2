"""Error types for the carquet Python package."""

from __future__ import annotations

import enum

__all__ = ["CarquetError", "StatusCode"]


class StatusCode(enum.IntEnum):
    """Mirror of ``carquet_status_t`` (include/carquet/error.h)."""

    OK = 0
    INVALID_ARGUMENT = 1
    OUT_OF_MEMORY = 2
    NOT_IMPLEMENTED = 3
    INTERNAL = 4
    FILE_NOT_FOUND = 10
    FILE_OPEN = 11
    FILE_READ = 12
    FILE_WRITE = 13
    FILE_SEEK = 14
    FILE_TRUNCATED = 15
    INVALID_MAGIC = 20
    INVALID_FOOTER = 21
    INVALID_SCHEMA = 22
    INVALID_METADATA = 23
    INVALID_PAGE = 24
    INVALID_ENCODING = 25
    VERSION_NOT_SUPPORTED = 26
    THRIFT_DECODE = 30
    THRIFT_ENCODE = 31
    THRIFT_INVALID_TYPE = 32
    THRIFT_TRUNCATED = 33
    DECODE = 40
    ENCODE = 41
    DICTIONARY_NOT_FOUND = 42
    INVALID_RLE = 43
    INVALID_DELTA = 44
    COMPRESSION = 50
    DECOMPRESSION = 51
    UNSUPPORTED_CODEC = 52
    INVALID_COMPRESSED_DATA = 53
    TYPE_MISMATCH = 60
    COLUMN_NOT_FOUND = 61
    ROW_GROUP_NOT_FOUND = 62
    END_OF_DATA = 63
    CHECKSUM = 70
    CRC_MISMATCH = 71
    INVALID_STATE = 80
    ALREADY_CLOSED = 81
    NOT_OPEN = 82
    PAGE_INDEX_REQUIRED = 90

    @classmethod
    def _missing_(cls, value):  # unknown future codes still map to something printable
        obj = int.__new__(cls, value)
        obj._name_ = f"UNKNOWN_{value}"
        obj._value_ = value
        return obj


class CarquetError(Exception):
    """Raised when a carquet C call fails.

    Attributes:
        code: numeric ``carquet_status_t`` value.
        status: the same code as a :class:`StatusCode` (``status.name`` gives the symbolic name).
        message: the detailed message produced by the library (never empty: falls back to ``carquet_status_string``).
        offset, column_index, row_group_index: optional context (``-1`` or ``None`` when not applicable).
    """

    def __init__(
        self,
        code: int,
        message: str,
        *,
        offset: int | None = None,
        column_index: int | None = None,
        row_group_index: int | None = None,
    ) -> None:
        self.code = int(code)
        self.status = StatusCode(self.code)
        self.message = message
        self.offset = offset
        self.column_index = column_index
        self.row_group_index = row_group_index
        super().__init__(self.__str__())

    def __str__(self) -> str:
        text = f"{self.message} [{self.status.name} ({self.code})]"
        ctx = []
        if self.row_group_index is not None and self.row_group_index >= 0:
            ctx.append(f"row group {self.row_group_index}")
        if self.column_index is not None and self.column_index >= 0:
            ctx.append(f"column {self.column_index}")
        if self.offset is not None and self.offset >= 0:
            ctx.append(f"offset {self.offset}")
        if ctx:
            text += " (" + ", ".join(ctx) + ")"
        return text
