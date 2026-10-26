"""ctypes binding to libcarquet: loader, struct layouts and prototypes.

Struct layouts mirror include/carquet/{carquet,error}.h field by field; the
tests in python/tests check the writer-option defaults against the header's
documented values, which catches any drift in field order or type.
"""

from __future__ import annotations

import ctypes
import os
import sys
from ctypes import (
    POINTER,
    Structure,
    c_bool,
    c_char,
    c_char_p,
    c_int,
    c_int32,
    c_int64,
    c_size_t,
    c_void_p,
)
from pathlib import Path
from typing import Optional

from ._arrow import ArrowArray, ArrowSchema
from .errors import CarquetError

__all__ = [
    "ERROR_MESSAGE_MAX",
    "CarquetErrorStruct",
    "FileInfoStruct",
    "ReaderOptions",
    "WriterOptions",
    "check",
    "find_library_path",
    "fsencode",
    "get_lib",
]

ERROR_MESSAGE_MAX = 256  # CARQUET_ERROR_MESSAGE_MAX
CREATED_BY_MAX = 256  # CARQUET_CREATED_BY_MAX


# --------------------------------------------------------------------------
# Struct layouts
# --------------------------------------------------------------------------
class CarquetErrorStruct(Structure):
    """carquet_error_t (error.h)."""

    _fields_ = [
        ("code", c_int),  # carquet_status_t (enum)
        ("message", c_char * ERROR_MESSAGE_MAX),
        ("file", c_char_p),  # const char*
        ("line", c_int),
        ("function", c_char_p),  # const char*
        ("offset", c_int64),
        ("column_index", c_int32),
        ("row_group_index", c_int32),
    ]

    def __init__(self) -> None:  # matches CARQUET_ERROR_INIT
        super().__init__()
        self.code = 0
        self.offset = -1
        self.column_index = -1
        self.row_group_index = -1


class FileInfoStruct(Structure):
    """carquet_file_info_t (carquet.h)."""

    _fields_ = [
        ("file_size", c_int64),
        ("num_rows", c_int64),
        ("num_row_groups", c_int32),
        ("num_columns", c_int32),
        ("version", c_int32),
        ("created_by", c_char * CREATED_BY_MAX),
    ]


class ReaderOptions(Structure):
    """carquet_reader_options_t (carquet.h)."""

    _fields_ = [
        ("use_mmap", c_bool),
        ("verify_checksums", c_bool),
        ("buffer_size", c_size_t),
        ("num_threads", c_int32),
    ]


class WriterOptions(Structure):
    """carquet_writer_options_t (carquet.h) -- field order and types matter."""

    _fields_ = [
        ("compression", c_int),  # carquet_compression_t
        ("compression_level", c_int32),
        ("row_group_size", c_int64),
        ("page_size", c_int64),
        ("write_statistics", c_bool),
        ("write_crc", c_bool),
        ("write_page_index", c_bool),
        ("write_bloom_filters", c_bool),
        ("dictionary_encoding", c_int),  # carquet_encoding_t
        ("dictionary_page_size", c_int64),
        ("created_by", c_char_p),  # const char*
        ("max_rows_per_page", c_int64),
        ("write_arrow_schema", c_bool),
        ("data_page_version", c_int32),
        ("coerce_timestamps", c_bool),
        ("coerce_timestamp_unit", c_int),  # carquet_time_unit_t
        ("allow_timestamp_truncation", c_bool),
        ("write_batch_size", c_int64),
        ("file_format_version", c_int32),
        ("async_io", c_bool),
    ]


# --------------------------------------------------------------------------
# Library loading
# --------------------------------------------------------------------------
def _library_names() -> list[str]:
    if sys.platform == "win32":
        return ["carquet.dll", "libcarquet.dll"]
    if sys.platform == "darwin":
        return ["libcarquet.dylib", "libcarquet.so"]
    return ["libcarquet.so", "libcarquet.dylib"]


def find_library_path() -> Optional[Path]:
    """Locate the bundled shared library.

    Search order: ``$CARQUET_LIBRARY`` (explicit file), the package directory
    (where the wheel puts it), then ``lib/``, ``lib64/`` and ``bin/`` under the
    package (layouts produced by a plain ``cmake --install`` with the default
    GNUInstallDirs destinations), and finally the system search path.
    """
    env = os.environ.get("CARQUET_LIBRARY")
    if env:
        p = Path(env)
        if p.is_file():
            return p
    pkg = Path(__file__).resolve().parent
    names = _library_names()
    for sub in ("", "lib", "lib64", "bin", "carquet.libs", ".libs"):
        base = pkg / sub if sub else pkg
        for name in names:
            cand = base / name
            if cand.is_file():
                return cand
    try:
        from ctypes.util import find_library

        found = find_library("carquet")
    except Exception:  # pragma: no cover - defensive
        found = None
    return Path(found) if found else None


_LIB: Optional[ctypes.CDLL] = None


def _load() -> ctypes.CDLL:
    path = find_library_path()
    if path is None:
        raise ImportError(
            "carquet: shared library not found (looked for "
            + ", ".join(_library_names())
            + f" next to {Path(__file__).resolve().parent}). "
            "Reinstall the wheel, or point CARQUET_LIBRARY at libcarquet."
        )
    if sys.platform == "win32" and hasattr(os, "add_dll_directory"):
        try:
            os.add_dll_directory(str(path.parent))
        except OSError:
            pass
    lib = ctypes.CDLL(str(path))
    _declare(lib)
    return lib


def get_lib() -> ctypes.CDLL:
    """Return the loaded ``libcarquet`` (loaded on first use)."""
    global _LIB
    if _LIB is None:
        _LIB = _load()
    return _LIB


# --------------------------------------------------------------------------
# Prototypes
# --------------------------------------------------------------------------
_ErrP = POINTER(CarquetErrorStruct)


def _declare(lib: ctypes.CDLL) -> None:
    def proto(name, restype, *argtypes):
        fn = getattr(lib, name)
        fn.restype = restype
        fn.argtypes = list(argtypes)
        return fn

    # version / status
    proto("carquet_version", c_char_p)
    proto("carquet_status_string", c_char_p, c_int)
    proto("carquet_compression_name", c_char_p, c_int)

    # utilities
    proto("carquet_get_file_info", c_int, c_char_p, POINTER(FileInfoStruct), _ErrP)
    proto("carquet_validate_file", c_int, c_char_p, _ErrP)

    # reader
    proto("carquet_reader_options_init", None, POINTER(ReaderOptions))
    proto("carquet_reader_open", c_void_p, c_char_p, POINTER(ReaderOptions), _ErrP)
    proto("carquet_reader_close", None, c_void_p)
    proto("carquet_reader_schema", c_void_p, c_void_p)
    proto("carquet_reader_num_rows", c_int64, c_void_p)
    proto("carquet_reader_num_row_groups", c_int32, c_void_p)
    proto("carquet_reader_num_columns", c_int32, c_void_p)
    proto(
        "carquet_reader_read_arrow",
        c_int,
        c_void_p,
        c_int32,
        POINTER(ArrowSchema),
        POINTER(ArrowArray),
        _ErrP,
    )

    # schema <-> Arrow
    proto("carquet_arrow_export_schema", c_int, c_void_p, POINTER(ArrowSchema), _ErrP)
    proto("carquet_arrow_import_schema", c_int, POINTER(ArrowSchema), POINTER(c_void_p), _ErrP)
    proto("carquet_schema_free", None, c_void_p)

    # writer
    proto("carquet_writer_options_init", None, POINTER(WriterOptions))
    proto("carquet_writer_create", c_void_p, c_char_p, c_void_p, POINTER(WriterOptions), _ErrP)
    proto(
        "carquet_writer_write_arrow",
        c_int,
        c_void_p,
        POINTER(ArrowArray),
        POINTER(ArrowSchema),
        _ErrP,
    )
    proto("carquet_writer_new_row_group", c_int, c_void_p)
    proto("carquet_writer_close", c_int, c_void_p)
    proto("carquet_writer_abort", None, c_void_p)


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def fsencode(path) -> bytes:
    return os.fsencode(os.fspath(path))


def status_string(code: int) -> str:
    s = get_lib().carquet_status_string(int(code))
    return s.decode("utf-8", "replace") if s else f"status {code}"


def check(status: int, err: Optional[CarquetErrorStruct] = None, context: str = "") -> None:
    """Raise :class:`CarquetError` when ``status != CARQUET_OK``."""
    if status == 0:
        return
    code = int(status)
    message = ""
    offset = column = row_group = None
    if err is not None:
        if err.code != 0:
            code = int(err.code)
        message = err.message.decode("utf-8", "replace").strip()
        offset, column, row_group = int(err.offset), int(err.column_index), int(err.row_group_index)
    if not message:
        message = status_string(code)
    if context:
        message = f"{context}: {message}"
    raise CarquetError(code, message, offset=offset, column_index=column, row_group_index=row_group)
