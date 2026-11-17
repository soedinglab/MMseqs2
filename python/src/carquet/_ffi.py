"""ctypes binding to libcarquet: loader, struct layouts and prototypes.

Struct layouts mirror include/carquet/{carquet,types,error}.h field by field.
``python/tests/test_ffi.py`` checks every ``sizeof`` against the values the C
compiler reports for the same structs (and the writer-option defaults against
the header's documented values), which catches any drift in field order or
type.

Every prototype is declared on a plain :class:`ctypes.CDLL`, so the GIL is
released for the duration of each C call: decoding, decompression and encoding
run without holding it.
"""

from __future__ import annotations

import ctypes
import os
import sys
from ctypes import (
    CFUNCTYPE,
    POINTER,
    Structure,
    Union,
    c_bool,
    c_char,
    c_char_p,
    c_double,
    c_float,
    c_int,
    c_int8,
    c_int16,
    c_int32,
    c_int64,
    c_size_t,
    c_uint8,
    c_void_p,
)
from pathlib import Path
from typing import Optional

from ._arrow import ArrowArray, ArrowArrayStream, ArrowSchema
from .errors import CarquetError

ERROR_MESSAGE_MAX = 256  # CARQUET_ERROR_MESSAGE_MAX
CREATED_BY_MAX = 256  # CARQUET_CREATED_BY_MAX
GEOSPATIAL_CRS_MAX = 128  # CARQUET_GEOSPATIAL_CRS_MAX
MAX_GEOSPATIAL_TYPES = 64  # CARQUET_MAX_GEOSPATIAL_TYPES

# carquet_status_t values used by the bindings themselves
OK = 0
END_OF_DATA = 63


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


class _DecimalParams(Structure):
    _fields_ = [("precision", c_int32), ("scale", c_int32)]


class _IntegerParams(Structure):
    _fields_ = [("bit_width", c_int8), ("is_signed", c_bool)]


class _TimeParams(Structure):
    _fields_ = [("unit", c_int), ("is_adjusted_to_utc", c_bool)]


class _VariantParams(Structure):
    _fields_ = [("specification_version", c_int8)]


class _GeometryParams(Structure):
    _fields_ = [("crs", c_char * GEOSPATIAL_CRS_MAX)]


class _GeographyParams(Structure):
    _fields_ = [
        ("crs", c_char * GEOSPATIAL_CRS_MAX),
        ("algorithm", c_int),
        ("has_algorithm", c_bool),
    ]


class _LogicalParams(Union):
    _fields_ = [
        ("decimal", _DecimalParams),
        ("integer", _IntegerParams),
        ("time", _TimeParams),
        ("timestamp", _TimeParams),
        ("variant", _VariantParams),
        ("geometry", _GeometryParams),
        ("geography", _GeographyParams),
    ]


class LogicalTypeStruct(Structure):
    """carquet_logical_type_t (types.h)."""

    _fields_ = [("id", c_int), ("params", _LogicalParams)]


class RowGroupMetadataStruct(Structure):
    """carquet_row_group_metadata_t."""

    _fields_ = [
        ("num_rows", c_int64),
        ("total_byte_size", c_int64),
        ("total_compressed_size", c_int64),
    ]


class ColumnStatisticsStruct(Structure):
    """carquet_column_statistics_t."""

    _fields_ = [
        ("has_min_max", c_bool),
        ("has_null_count", c_bool),
        ("has_distinct_count", c_bool),
        ("null_count", c_int64),
        ("distinct_count", c_int64),
        ("num_values", c_int64),
        ("min_value", c_void_p),
        ("max_value", c_void_p),
        ("min_value_size", c_int32),
        ("max_value_size", c_int32),
    ]


class ColumnChunkMetadataStruct(Structure):
    """carquet_column_chunk_metadata_t."""

    _fields_ = [
        ("type", c_int),
        ("codec", c_int),
        ("num_values", c_int64),
        ("total_compressed_size", c_int64),
        ("total_uncompressed_size", c_int64),
        ("data_page_offset", c_int64),
        ("has_dictionary_page", c_bool),
        ("dictionary_page_offset", c_int64),
        ("num_encodings", c_int32),
        ("encodings", c_int * 4),
        ("has_bloom_filter", c_bool),
        ("has_column_index", c_bool),
        ("has_offset_index", c_bool),
    ]


class GeospatialStatisticsStruct(Structure):
    """carquet_geospatial_statistics_t."""

    _fields_ = [
        ("has_bbox", c_bool),
        ("xmin", c_double),
        ("xmax", c_double),
        ("ymin", c_double),
        ("ymax", c_double),
        ("has_z", c_bool),
        ("zmin", c_double),
        ("zmax", c_double),
        ("has_m", c_bool),
        ("mmin", c_double),
        ("mmax", c_double),
        ("num_geometry_types", c_int32),
        ("geometry_types", c_int32 * MAX_GEOSPATIAL_TYPES),
    ]


class PageStatsStruct(Structure):
    """carquet_page_stats_t."""

    _fields_ = [
        ("null_count", c_int64),
        ("min_value", c_void_p),
        ("min_value_size", c_int32),
        ("max_value", c_void_p),
        ("max_value_size", c_int32),
        ("is_null_page", c_bool),
        ("has_null_count", c_bool),
    ]


class PageLocationStruct(Structure):
    """carquet_page_location_t."""

    _fields_ = [
        ("offset", c_int64),
        ("compressed_size", c_int32),
        ("first_row_index", c_int64),
    ]


class FilterClauseStruct(Structure):
    """carquet_filter_clause_t."""

    _fields_ = [
        ("column_index", c_int32),
        ("op", c_int),
        ("value", c_void_p),
        ("value_size", c_int32),
        ("lo", c_void_p),
        ("lo_size", c_int32),
        ("hi", c_void_p),
        ("hi_size", c_int32),
        ("has_lo", c_bool),
        ("has_hi", c_bool),
        ("values", c_void_p),
        ("value_count", c_int32),
    ]


class SortingColumnStruct(Structure):
    """carquet_sorting_column_t."""

    _fields_ = [("column_index", c_int32), ("descending", c_bool), ("nulls_first", c_bool)]


class ByteArrayStruct(Structure):
    """carquet_byte_array_t."""

    _fields_ = [("data", c_void_p), ("length", c_int32)]


class CpuInfoStruct(Structure):
    """carquet_cpu_info_t."""

    _fields_ = [
        ("has_sse2", c_bool),
        ("has_sse41", c_bool),
        ("has_sse42", c_bool),
        ("has_avx", c_bool),
        ("has_avx2", c_bool),
        ("has_avx512f", c_bool),
        ("has_avx512bw", c_bool),
        ("has_avx512vl", c_bool),
        ("has_avx512vbmi", c_bool),
        ("has_neon", c_bool),
        ("has_sve", c_bool),
        ("sve_vector_length", c_int),
    ]


# bool (*)(const carquet_reader_t*, int32_t row_group_index, void* user_data)
RowGroupFilterFn = CFUNCTYPE(c_bool, c_void_p, c_int32, c_void_p)


class BatchReaderConfig(Structure):
    """carquet_batch_reader_config_t."""

    _fields_ = [
        ("batch_size", c_int32),
        ("num_threads", c_int32),
        ("use_mmap", c_bool),
        ("column_indices", POINTER(c_int32)),
        ("num_columns", c_int32),
        ("column_names", POINTER(c_char_p)),
        ("num_column_names", c_int32),
        ("row_group_filter", RowGroupFilterFn),
        ("row_group_filter_ctx", c_void_p),
        ("preserve_dictionaries", c_bool),
        ("thread_pool", c_void_p),
    ]


class PredicateStruct(Structure):
    """carquet_predicate_t."""

    _fields_ = [
        ("clauses", POINTER(FilterClauseStruct)),
        ("conjunction_sizes", POINTER(c_int32)),
        ("num_conjunctions", c_int32),
    ]


class ArrowStreamOptions(Structure):
    """carquet_arrow_stream_options_t."""

    _fields_ = [
        ("field_indices", POINTER(c_int32)),
        ("num_fields", c_int32),
        ("batch_size", c_int32),
        ("num_threads", c_int32),
        ("close_reader", c_bool),
        ("row_groups", POINTER(c_int32)),
        ("num_row_groups", c_int32),
        ("row_offset", c_int64),
        ("row_limit", c_int64),
        ("predicate", POINTER(PredicateStruct)),
        ("restore_arrow_schema", c_bool),
    ]


class AllocatorStruct(Structure):
    """carquet_allocator_t (only ``free`` is used, to release writer buffers)."""

    _fields_ = [
        ("malloc", CFUNCTYPE(c_void_p, c_size_t, c_void_p)),
        ("realloc", CFUNCTYPE(c_void_p, c_void_p, c_size_t, c_void_p)),
        ("free", CFUNCTYPE(None, c_void_p, c_void_p)),
        ("ctx", c_void_p),
    ]


# sizeof() of each struct as reported by the C compiler for the supported
# 64-bit targets; asserted by the test suite.
EXPECTED_SIZES = {
    "LogicalTypeStruct": 140,
    "BatchReaderConfig": 80,
    "ColumnStatisticsStruct": 56,
    "ColumnChunkMetadataStruct": 80,
    "GeospatialStatisticsStruct": 352,
    "PageStatsStruct": 40,
    "PageLocationStruct": 24,
    "FilterClauseStruct": 72,
    "SortingColumnStruct": 8,
    "CpuInfoStruct": 16,
    "WriterOptions": 96,
    "ReaderOptions": 24,
    "RowGroupMetadataStruct": 24,
    "ByteArrayStruct": 16,
    "AllocatorStruct": 32,
    "ArrowStreamOptions": 72,
    "PredicateStruct": 24,
}


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
    def proto(name, restype, *argtypes, optional=False):
        try:
            fn = getattr(lib, name)
        except AttributeError:
            if optional:  # symbol added after the library this binding was built for
                return None
            raise
        fn.restype = restype
        fn.argtypes = list(argtypes)
        return fn

    vp = c_void_p

    # version / status / runtime
    proto("carquet_version", c_char_p)
    proto("carquet_status_string", c_char_p, c_int)
    proto("carquet_compression_name", c_char_p, c_int)
    proto("carquet_get_cpu_info", POINTER(CpuInfoStruct))
    proto("carquet_get_allocator", POINTER(AllocatorStruct))

    # utilities
    proto("carquet_get_file_info", c_int, c_char_p, POINTER(FileInfoStruct), _ErrP)
    proto("carquet_validate_file", c_int, c_char_p, _ErrP)

    # schema
    proto("carquet_schema_free", None, vp)
    proto("carquet_schema_num_columns", c_int32, vp)
    proto("carquet_schema_num_elements", c_int32, vp)
    proto("carquet_schema_get_element", vp, vp, c_int32)
    proto("carquet_schema_find_column", c_int32, vp, c_char_p)
    proto("carquet_schema_max_def_level", c_int16, vp, c_int32)
    proto("carquet_schema_max_rep_level", c_int16, vp, c_int32)
    proto("carquet_schema_column_name", c_char_p, vp, c_int32)
    proto("carquet_schema_column_type", c_int, vp, c_int32)
    proto("carquet_schema_column_path", c_int32, vp, c_int32, POINTER(c_char_p), c_int32)
    proto("carquet_schema_node_name", c_char_p, vp)
    proto("carquet_schema_node_is_leaf", c_bool, vp)
    proto("carquet_schema_node_physical_type", c_int, vp)
    proto("carquet_schema_node_logical_type", POINTER(LogicalTypeStruct), vp)
    proto("carquet_schema_node_repetition", c_int, vp)
    proto("carquet_schema_node_type_length", c_int32, vp)

    # reader
    proto("carquet_reader_options_init", None, POINTER(ReaderOptions))
    proto("carquet_reader_open", vp, c_char_p, POINTER(ReaderOptions), _ErrP)
    proto("carquet_reader_open_buffer", vp, vp, c_size_t, POINTER(ReaderOptions), _ErrP)
    proto("carquet_reader_close", None, vp)
    proto("carquet_reader_schema", vp, vp)
    proto("carquet_reader_num_rows", c_int64, vp)
    proto("carquet_reader_num_row_groups", c_int32, vp)
    proto("carquet_reader_num_columns", c_int32, vp)
    proto("carquet_reader_is_mmap", c_bool, vp)
    proto("carquet_reader_row_group_metadata", c_int, vp, c_int32, POINTER(RowGroupMetadataStruct))
    proto("carquet_reader_column_chunk_metadata", c_int, vp, c_int32, c_int32, POINTER(ColumnChunkMetadataStruct))
    proto("carquet_reader_column_statistics", c_int, vp, c_int32, c_int32, POINTER(ColumnStatisticsStruct))
    proto("carquet_reader_geospatial_statistics", c_int, vp, c_int32, c_int32, POINTER(GeospatialStatisticsStruct))
    proto("carquet_reader_row_group_matches", c_int, vp, c_int32, c_int32, c_int, vp, c_int32, POINTER(c_bool))
    proto("carquet_reader_num_metadata", c_int32, vp)
    proto("carquet_reader_get_metadata", c_int, vp, c_int32, POINTER(c_char_p), POINTER(c_char_p))
    proto("carquet_reader_find_metadata", c_char_p, vp, c_char_p)
    proto("carquet_reader_column_num_metadata", c_int32, vp, c_int32)
    proto("carquet_reader_column_get_metadata", c_int, vp, c_int32, c_int32, POINTER(c_char_p), POINTER(c_char_p))
    proto("carquet_reader_column_arrow_type_refinement", c_int, vp, c_int32)
    proto("carquet_reader_read_arrow", c_int, vp, c_int32, POINTER(ArrowSchema), POINTER(ArrowArray), _ErrP)
    proto(
        "carquet_reader_read_arrow_fields",
        c_int,
        vp,
        c_int32,
        POINTER(c_int32),
        c_int32,
        POINTER(ArrowSchema),
        POINTER(ArrowArray),
        _ErrP,
        optional=True,
    )

    # bloom filters
    proto("carquet_reader_get_bloom_filter", vp, vp, c_int32, c_int32, _ErrP)
    proto("carquet_bloom_filter_check_i32", c_bool, vp, c_int32)
    proto("carquet_bloom_filter_check_i64", c_bool, vp, c_int64)
    proto("carquet_bloom_filter_check_float", c_bool, vp, c_float)
    proto("carquet_bloom_filter_check_double", c_bool, vp, c_double)
    proto("carquet_bloom_filter_check_bytes", c_bool, vp, vp, c_size_t)
    proto("carquet_bloom_filter_size", c_size_t, vp)
    proto("carquet_bloom_filter_destroy", None, vp)

    # page index
    proto("carquet_reader_get_column_index", vp, vp, c_int32, c_int32, _ErrP)
    proto("carquet_column_index_num_pages", c_int32, vp)
    proto("carquet_column_index_get_page_stats", c_int, vp, c_int32, POINTER(PageStatsStruct))
    proto("carquet_column_index_boundary_order", c_int32, vp)
    proto("carquet_column_index_free", None, vp)
    proto("carquet_reader_get_offset_index", vp, vp, c_int32, c_int32, _ErrP)
    proto("carquet_offset_index_num_pages", c_int32, vp)
    proto("carquet_offset_index_get_page_location", c_int, vp, c_int32, POINTER(PageLocationStruct))
    proto("carquet_offset_index_free", None, vp)

    # batch reader
    proto("carquet_batch_reader_config_init", None, POINTER(BatchReaderConfig))
    proto("carquet_batch_reader_create", vp, vp, POINTER(BatchReaderConfig), _ErrP)
    proto("carquet_batch_reader_next", c_int, vp, POINTER(c_void_p))
    proto("carquet_batch_reader_last_error", POINTER(CarquetErrorStruct), vp)
    proto("carquet_batch_reader_free", None, vp)
    proto("carquet_batch_reader_set_page_filter", c_int, vp, POINTER(FilterClauseStruct), c_int32)
    proto("carquet_batch_reader_rows_skipped", c_int64, vp)
    proto("carquet_batch_reader_set_row_range", c_int, vp, c_int64, c_int64)
    proto("carquet_row_batch_num_rows", c_int64, vp)
    proto("carquet_row_batch_num_columns", c_int32, vp)
    proto("carquet_row_batch_free", None, vp)
    proto("carquet_thread_pool_create", vp, c_int32)
    proto("carquet_thread_pool_destroy", None, vp)

    # schema <-> Arrow
    proto("carquet_arrow_export_schema", c_int, vp, POINTER(ArrowSchema), _ErrP)
    proto("carquet_arrow_export_batch", c_int, vp, vp, POINTER(ArrowSchema), POINTER(ArrowArray), _ErrP)
    proto("carquet_arrow_import_schema", c_int, POINTER(ArrowSchema), POINTER(c_void_p), _ErrP)
    proto("carquet_arrow_stream_options_init", None, POINTER(ArrowStreamOptions), optional=True)
    proto("carquet_reader_predicate_row_groups", c_int32, vp, POINTER(PredicateStruct), POINTER(c_int32), c_int32,
          _ErrP, optional=True)
    proto("carquet_arrow_ref_create", vp, POINTER(ArrowArray), optional=True)
    proto("carquet_arrow_ref_array", POINTER(ArrowArray), vp, optional=True)
    proto("carquet_arrow_ref_export", c_int, vp, c_int64, c_int64, POINTER(ArrowArray), optional=True)
    proto("carquet_arrow_ref_export_child", c_int, vp, c_int32, c_int64, c_int64, POINTER(ArrowArray), optional=True)
    proto("carquet_arrow_ref_export_select", c_int, vp, POINTER(c_int32), c_int32, c_int64, c_int64,
          POINTER(ArrowArray), optional=True)
    proto("carquet_arrow_ref_release", None, vp, optional=True)
    proto("carquet_arrow_schema_copy", c_int, POINTER(ArrowSchema), POINTER(ArrowSchema), optional=True)
    proto("carquet_arrow_schema_copy_child", c_int, POINTER(ArrowSchema), c_int32, POINTER(ArrowSchema), optional=True)
    proto("carquet_arrow_schema_select", c_int, POINTER(ArrowSchema), POINTER(c_int32), c_int32, optional=True)
    proto("carquet_arrow_stream_from_arrays", c_int, POINTER(ArrowSchema), POINTER(ArrowArray), c_int64,
          POINTER(ArrowArrayStream), optional=True)
    proto("carquet_arrow_array_build", c_int, c_int64, c_int64, c_int32, POINTER(c_void_p), POINTER(c_int64),
          c_int32, POINTER(ArrowArray), POINTER(ArrowArray), POINTER(ArrowArray), optional=True)
    proto("carquet_arrow_schema_build", c_int, c_char_p, c_char_p, c_int64, c_char_p, c_int32,
          POINTER(ArrowSchema), POINTER(ArrowSchema), POINTER(ArrowSchema), optional=True)
    proto("carquet_writer_set_arrow_schema", c_int, vp, POINTER(ArrowSchema), optional=True)
    proto("carquet_arrow_capsule_init", None, vp, optional=True)
    proto("carquet_arrow_struct_alloc", vp, c_size_t, optional=True)
    proto(
        "carquet_reader_export_arrow_stream",
        c_int,
        vp,
        POINTER(ArrowStreamOptions),
        POINTER(ArrowArrayStream),
        _ErrP,
        optional=True,
    )

    # writer
    proto("carquet_writer_options_init", None, POINTER(WriterOptions))
    proto("carquet_writer_create", vp, c_char_p, vp, POINTER(WriterOptions), _ErrP)
    proto("carquet_writer_create_buffer", vp, vp, POINTER(WriterOptions), _ErrP)
    proto("carquet_writer_open_append", vp, c_char_p, vp, POINTER(WriterOptions), _ErrP)
    proto("carquet_writer_get_buffer", c_int, vp, POINTER(c_void_p), POINTER(c_size_t))
    proto("carquet_writer_write_arrow", c_int, vp, POINTER(ArrowArray), POINTER(ArrowSchema), _ErrP)
    proto("carquet_writer_new_row_group", c_int, vp)
    proto("carquet_writer_num_columns", c_int32, vp)
    proto("carquet_writer_add_metadata", c_int, vp, c_char_p, c_char_p)
    proto("carquet_writer_set_column_encoding", c_int, vp, c_int32, c_int)
    proto("carquet_writer_set_column_compression", c_int, vp, c_int32, c_int, c_int32)
    proto("carquet_writer_set_column_page_size", c_int, vp, c_int32, c_int64)
    proto("carquet_writer_set_max_statistics_size", c_int, vp, c_int64)
    proto("carquet_writer_set_column_statistics", c_int, vp, c_int32, c_bool)
    proto("carquet_writer_set_column_bloom_filter", c_int, vp, c_int32, c_bool)
    proto("carquet_writer_set_column_bloom_filter_options", c_int, vp, c_int32, c_bool, c_int64, c_double)
    proto("carquet_writer_set_sorting_columns", c_int, vp, POINTER(SortingColumnStruct), c_int32)
    proto("carquet_writer_close", c_int, vp)
    proto("carquet_writer_abort", None, vp)


def has_symbol(name: str) -> bool:
    """True when the loaded library exports ``name`` (newer optional entry points)."""
    fn = getattr(get_lib(), name, None)
    return fn is not None and fn.argtypes is not None


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def fsencode(path) -> bytes:
    return os.fsencode(os.fspath(path))


def status_string(code: int) -> str:
    s = get_lib().carquet_status_string(int(code))
    return s.decode("utf-8", "replace") if s else f"status {code}"


def decode(raw: Optional[bytes]) -> str:
    return raw.decode("utf-8", "replace") if raw else ""


def check(status: int, err: Optional[CarquetErrorStruct] = None, context: str = "") -> None:
    """Raise :class:`CarquetError` when ``status != CARQUET_OK``."""
    if status == 0:
        return
    raise make_error(int(status), err, context)


def make_error(code: int, err: Optional[CarquetErrorStruct] = None, context: str = "") -> CarquetError:
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
    return CarquetError(code, message, offset=offset, column_index=column, row_group_index=row_group)


def free_buffer(ptr: int) -> None:
    """Release memory the library allocated and handed over (e.g. a writer buffer)."""
    alloc = get_lib().carquet_get_allocator().contents
    alloc.free(ptr, alloc.ctx)


def buffer_address(obj) -> tuple[int, int, object]:
    """Return ``(address, size, keepalive)`` for a bytes-like object without copying.

    ``keepalive`` must be referenced for as long as the address is used. A
    read-only buffer that ctypes cannot address directly is copied once.
    """
    if isinstance(obj, bytes):
        return ctypes.cast(c_char_p(obj), c_void_p).value or 0, len(obj), obj
    try:
        import pyarrow as pa  # noqa: F401 - optional

        if isinstance(obj, pa.Buffer):
            return obj.address, obj.size, obj
    except ImportError:  # pragma: no cover
        pass
    mv = memoryview(obj).cast("B")
    n = mv.nbytes
    if n == 0:
        return 0, 0, obj
    if not mv.readonly:
        holder = (c_uint8 * n).from_buffer(mv)
        return ctypes.addressof(holder), n, (holder, mv, obj)
    data = mv.tobytes()
    return ctypes.cast(c_char_p(data), c_void_p).value or 0, n, data

