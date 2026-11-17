"""ctypes layer: struct layouts, library defaults, version and runtime helpers."""

from __future__ import annotations

import ctypes
import re
from pathlib import Path

import carquet
from carquet import _ffi

HEADER = Path(__file__).resolve().parents[2] / "include" / "carquet" / "carquet.h"


def test_struct_sizes_match_c():
    if ctypes.sizeof(ctypes.c_void_p) != 8:  # pragma: no cover - 64-bit wheels only
        return
    for name, size in _ffi.EXPECTED_SIZES.items():
        assert ctypes.sizeof(getattr(_ffi, name)) == size, name
    assert ctypes.sizeof(_ffi.CarquetErrorStruct) == 4 + 256 + 4 + 8 + 4 + 4 + 8 + 8 + 4 + 4
    assert ctypes.sizeof(_ffi.FileInfoStruct) == 8 + 8 + 4 + 4 + 4 + 256 + 4


def test_struct_offsets():
    """Field offsets reported by the C compiler for the trickier layouts."""
    assert _ffi.LogicalTypeStruct.params.offset == 4
    assert _ffi.BatchReaderConfig.preserve_dictionaries.offset == 64
    assert _ffi.BatchReaderConfig.thread_pool.offset == 72
    assert _ffi.ColumnStatisticsStruct.min_value_size.offset == 48
    assert _ffi.ColumnChunkMetadataStruct.encodings.offset == 60
    assert _ffi.ColumnChunkMetadataStruct.has_offset_index.offset == 78
    assert _ffi.GeospatialStatisticsStruct.geometry_types.offset == 92
    assert _ffi.PageStatsStruct.has_null_count.offset == 37
    assert _ffi.FilterClauseStruct.values.offset == 56
    assert _ffi.FilterClauseStruct.value_count.offset == 64
    assert _ffi.CpuInfoStruct.sve_vector_length.offset == 12


def test_writer_and_reader_defaults():
    """Option defaults documented in carquet.h, read back through the ctypes layout."""
    lib = _ffi.get_lib()
    opts = _ffi.WriterOptions()
    lib.carquet_writer_options_init(ctypes.byref(opts))
    assert opts.compression == 0
    assert opts.compression_level == 0
    assert opts.row_group_size == 128 * 1024 * 1024
    assert opts.page_size == 1024 * 1024
    assert opts.write_statistics is True
    assert opts.write_crc is True
    assert opts.write_page_index is False
    assert opts.dictionary_encoding == 8
    assert opts.dictionary_page_size == 1024 * 1024
    assert opts.created_by == b"Carquet"
    assert opts.max_rows_per_page == 0
    assert opts.data_page_version == 1
    assert opts.coerce_timestamp_unit == 1
    assert opts.write_batch_size == 0
    assert opts.file_format_version == 2
    assert opts.async_io is True

    ropts = _ffi.ReaderOptions()
    lib.carquet_reader_options_init(ctypes.byref(ropts))
    assert ropts.use_mmap is False
    assert ropts.verify_checksums is True
    assert ropts.buffer_size == 65536
    assert ropts.num_threads == 0

    cfg = _ffi.BatchReaderConfig()
    lib.carquet_batch_reader_config_init(ctypes.byref(cfg))
    assert cfg.batch_size == 65536
    assert cfg.num_threads == 0
    assert not cfg.column_indices
    assert cfg.preserve_dictionaries is False


def test_version_matches_header():
    assert re.fullmatch(r"\d+\.\d+\.\d+", carquet.lib_version())
    assert carquet.__version__.startswith(carquet.lib_version())
    if HEADER.is_file():
        m = re.search(r'#define\s+CARQUET_VERSION_STRING\s+"([^"]+)"', HEADER.read_text())
        assert m and m.group(1) == carquet.lib_version()


def test_cpu_info():
    info = carquet.cpu_info()
    assert set(info) >= {"has_neon", "has_avx2", "has_sse42", "sve_vector_length"}
    assert all(isinstance(v, bool) for k, v in info.items() if k != "sve_vector_length")


def test_projected_read_entry_point():
    assert _ffi.has_symbol("carquet_reader_read_arrow_fields")


def test_status_strings():
    assert _ffi.status_string(0)
    err = _ffi.make_error(int(carquet.StatusCode.FILE_NOT_FOUND), None, "ctx")
    assert err.status == carquet.StatusCode.FILE_NOT_FOUND
    assert str(err).startswith("ctx: ")
