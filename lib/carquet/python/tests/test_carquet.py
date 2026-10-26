"""pytest suite for the carquet Python package (run against an installed wheel)."""

from __future__ import annotations

import ctypes
import re
from pathlib import Path

import pytest

import carquet
from carquet import CarquetError, StatusCode

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")

HEADER = Path(__file__).resolve().parents[2] / "include" / "carquet" / "carquet.h"
CODECS = ["none", "snappy", "zstd", "lz4", "brotli"]


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------
def make_table(n: int = 1000) -> "pa.Table":
    import datetime as dt

    ids = pa.array(range(n), pa.int64())
    values = pa.array([i * 0.5 for i in range(n)], pa.float64())
    names = pa.array([None if i % 7 == 3 else f"name-{i}" for i in range(n)], pa.string())
    epoch = dt.datetime(2024, 1, 1)
    ts = pa.array([epoch + dt.timedelta(microseconds=1_000_003 * i) for i in range(n)], pa.timestamp("us"))
    lists = pa.array(
        [None if i % 11 == 5 else [j for j in range(i % 4)] for i in range(n)],
        pa.list_(pa.int32()),
    )
    struct = pa.StructArray.from_arrays(
        [
            pa.array([i % 100 for i in range(n)], pa.int32()),
            pa.array([None if i % 13 == 0 else f"s{i}" for i in range(n)], pa.string()),
        ],
        names=["a", "b"],
    )
    return pa.table(
        {"id": ids, "value": values, "name": names, "ts": ts, "tags": lists, "point": struct}
    )


@pytest.fixture(scope="module")
def table() -> "pa.Table":
    return make_table()


@pytest.fixture
def written(tmp_path, table):
    path = tmp_path / "sample.parquet"
    carquet.write_table(table, path, compression="zstd")
    return path


def assert_same_data(actual: "pa.Table", expected: "pa.Table") -> None:
    """Logical equality: same columns, same values; schemas may differ in list child names / nullability."""
    assert actual.num_rows == expected.num_rows
    assert actual.column_names == expected.column_names
    assert actual.schema.equals(expected.schema, check_metadata=False) or actual.cast(expected.schema).schema.equals(
        expected.schema, check_metadata=False
    )
    casted = actual.cast(expected.schema)
    for name in expected.column_names:
        assert casted.column(name).equals(expected.column(name)), name


# --------------------------------------------------------------------------
# version / metadata (no pyarrow needed)
# --------------------------------------------------------------------------
def test_version_matches_header():
    assert re.fullmatch(r"\d+\.\d+\.\d+", carquet.lib_version())
    assert carquet.__version__.startswith(carquet.lib_version())
    if HEADER.is_file():
        m = re.search(r'#define\s+CARQUET_VERSION_STRING\s+"([^"]+)"', HEADER.read_text())
        assert m and m.group(1) == carquet.lib_version()


def test_struct_layouts():
    """Writer-option defaults documented in carquet.h, read back through the ctypes layout."""
    from carquet import _ffi

    assert ctypes.sizeof(_ffi.WriterOptions) == 96
    assert ctypes.sizeof(_ffi.ReaderOptions) == 24
    assert ctypes.sizeof(_ffi.CarquetErrorStruct) == 4 + 256 + 4 + 8 + 4 + 4 + 8 + 8 + 4 + 4  # incl. padding before `file`
    assert ctypes.sizeof(_ffi.FileInfoStruct) == 8 + 8 + 4 + 4 + 4 + 256 + 4

    lib = _ffi.get_lib()
    opts = _ffi.WriterOptions()
    lib.carquet_writer_options_init(ctypes.byref(opts))
    assert opts.compression == 0  # UNCOMPRESSED (file_writer.c; the header comment still says SNAPPY)
    assert opts.compression_level == 0
    assert opts.row_group_size == 128 * 1024 * 1024
    assert opts.page_size == 1024 * 1024
    assert opts.write_statistics is True
    assert opts.write_crc is True
    assert opts.write_page_index is False
    assert opts.dictionary_encoding == 8  # RLE_DICTIONARY (file_writer.c; header comment says PLAIN_DICTIONARY)
    assert opts.dictionary_page_size == 1024 * 1024
    assert opts.created_by == b"Carquet"
    assert opts.max_rows_per_page == 0
    assert opts.data_page_version == 1
    assert opts.coerce_timestamp_unit == 1  # MICROS
    assert opts.write_batch_size == 0
    assert opts.file_format_version == 2
    assert opts.async_io is True

    ropts = _ffi.ReaderOptions()
    lib.carquet_reader_options_init(ctypes.byref(ropts))
    assert ropts.use_mmap is False
    assert ropts.verify_checksums is True
    assert ropts.buffer_size == 65536
    assert ropts.num_threads == 0


def test_info_and_validate(written, table):
    fi = carquet.info(written)
    assert fi.num_rows == table.num_rows
    assert fi.num_columns == 7  # leaf columns: id, value, name, ts, tags.list.element, point.a, point.b
    assert fi.num_row_groups == 1
    assert fi.file_size == written.stat().st_size
    assert fi.version == 2
    assert "Carquet" in fi.created_by
    carquet.validate(written)  # must not raise


def test_missing_file_raises():
    with pytest.raises(CarquetError) as ei:
        carquet.info("/definitely/not/here.parquet")
    assert ei.value.message
    assert ei.value.code != 0
    assert ei.value.status in (StatusCode.FILE_NOT_FOUND, StatusCode.FILE_OPEN)
    with pytest.raises(CarquetError):
        carquet.validate("/definitely/not/here.parquet")
    with pytest.raises(CarquetError):
        carquet.read_table("/definitely/not/here.parquet")


def test_not_parquet_raises(tmp_path):
    bad = tmp_path / "bad.parquet"
    bad.write_bytes(b"this is not a parquet file at all" * 10)
    with pytest.raises(CarquetError) as ei:
        carquet.validate(bad)
    assert ei.value.message
    with pytest.raises(CarquetError):
        carquet.info(bad)


# --------------------------------------------------------------------------
# round trips
# --------------------------------------------------------------------------
@pytest.mark.parametrize("compression", CODECS)
def test_roundtrip(tmp_path, table, compression):
    path = tmp_path / f"rt_{compression}.parquet"
    try:
        carquet.write_table(table, path, compression=compression)
    except CarquetError as exc:
        if compression == "brotli" and exc.status == StatusCode.UNSUPPORTED_CODEC:
            pytest.skip("library built without brotli")
        raise
    out = carquet.read_table(path)
    assert_same_data(out, table)
    # and PyArrow agrees with what carquet wrote
    assert_same_data(pq.read_table(path), table)
    meta = pq.ParquetFile(path).metadata
    expected = {"none": {"UNCOMPRESSED"}, "lz4": {"LZ4", "LZ4_RAW"}}.get(compression, {compression.upper()})
    assert meta.row_group(0).column(0).compression in expected


def test_row_group_size(tmp_path, table):
    path = tmp_path / "rg.parquet"
    carquet.write_table(table, path, row_group_size=300)
    fi = carquet.info(path)
    assert fi.num_row_groups == 4  # 1000 rows / 300
    out = carquet.read_table(path)
    assert out.column("id").num_chunks == 4
    assert_same_data(out, table)
    assert pq.ParquetFile(path).metadata.num_row_groups == 4


def test_schema(written, table):
    s = carquet.schema(written)
    assert isinstance(s, pa.Schema)
    assert s.names == table.schema.names
    for got, want in zip(s, table.schema):
        assert got.name == want.name
        # equal up to nullability (and the Parquet "element" list-child name): a zero-copy cast must exist
        assert got.type == want.type or pa.array([], got.type).cast(want.type).type == want.type


def test_columns_selection(written, table):
    out = carquet.read_table(written, columns=["name", "id"])
    assert out.column_names == ["name", "id"]
    assert out.column("id").equals(table.column("id"))


def test_open_reader(written, table):
    with carquet.open_reader(written) as r:
        assert r.num_rows == table.num_rows
        assert r.num_row_groups == 1
        assert r.num_columns >= table.num_columns
        batch = r.read_row_group(0)
        assert isinstance(batch, pa.RecordBatch)
        assert batch.num_rows == table.num_rows
        assert r.schema.names == table.schema.names
        with pytest.raises(IndexError):
            r.read_row_group(5)
        assert "rows=" in repr(r)
    assert r.closed
    with pytest.raises(CarquetError):
        r.read_row_group(0)


def test_pyarrow_reads_carquet_output_and_back(tmp_path, table):
    path = tmp_path / "x.parquet"
    carquet.write_table(table, path, compression="snappy", write_page_index=True, data_page_version=2)
    assert_same_data(pq.read_table(path), table)
    # the other direction: a PyArrow-written file read by carquet
    path2 = tmp_path / "y.parquet"
    pq.write_table(table, path2, compression="zstd")
    assert_same_data(carquet.read_table(path2), table)
    assert_same_data(carquet.read_table(path2, use_mmap=False), table)


def test_empty_table(tmp_path, table):
    path = tmp_path / "empty.parquet"
    carquet.write_table(table.schema.empty_table(), path)
    out = carquet.read_table(path)
    assert out.num_rows == 0
    assert out.schema.names == table.schema.names
    assert carquet.info(path).num_rows == 0


def test_record_batch_input(tmp_path, table):
    path = tmp_path / "batch.parquet"
    carquet.write_table(table.to_batches()[0], path, compression="none")
    assert carquet.info(path).num_rows == table.num_rows


def test_bad_compression_name(tmp_path, table):
    with pytest.raises(ValueError):
        carquet.write_table(table, tmp_path / "z.parquet", compression="lzma")


def test_compression_level(tmp_path, table):
    path = tmp_path / "lvl.parquet"
    carquet.write_table(table, path, compression="zstd", compression_level=10)
    assert_same_data(carquet.read_table(path), table)
