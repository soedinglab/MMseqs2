"""Metadata: FileMetaData, statistics, schema, page indexes (most need no PyArrow)."""

from __future__ import annotations

import datetime as dt
import decimal
import uuid

import pytest

import carquet

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")


def test_file_metadata(tmp_path, table):
    path = tmp_path / "m.parquet"
    carquet.write_table(table, path, row_group_size=400, metadata={"who": "me"})
    md = carquet.read_metadata(path)
    assert md.num_rows == 1000
    assert md.num_row_groups == 3
    assert md.num_columns == 7
    assert md.created_by.startswith("Carquet")
    assert md.format_version == 2
    assert md.file_size == path.stat().st_size
    assert md.key_value_metadata["who"] == "me"
    assert md.metadata is md.key_value_metadata
    rg = md.row_group(0)
    assert rg.num_rows == 400 and rg.num_columns == 7 and len(rg) == 7
    assert rg.total_compressed_size > 0
    col = rg.column("name")
    assert col is rg["name"] and col is rg[2]
    assert col.path_in_schema == "name" and col.column_index == 2
    assert col.physical_type == "BYTE_ARRAY"
    assert col.compression == "ZSTD"
    assert col.num_values == 400
    assert col.total_uncompressed_size > 0 and col.compression_ratio > 0
    assert rg.column("point.b").path_in_schema == "point.b"
    d = md.to_dict()
    assert d["num_rows"] == 1000 and len(d["row_groups"]) == 3
    assert "Statistics(" in repr(col.statistics)
    with pytest.raises(IndexError):
        md.row_group(3)


def test_statistics_decoding(tmp_path):
    utc = dt.timezone.utc
    u = uuid.UUID("12345678-1234-5678-1234-567812345678")
    t = pa.table({
        "i": pa.array([5, -3, None, 9], pa.int32()),
        "u": pa.array([1, 2**32 - 1, 0, 7], pa.uint32()),
        "f": pa.array([1.5, -2.25, 3.0, None], pa.float64()),
        "s": pa.array(["pear", "apple", None, "zoo"]),
        "b": pa.array([True, False, True, None]),
        "d": pa.array([dt.date(2020, 1, 2), dt.date(2019, 5, 6), None, dt.date(2021, 1, 1)]),
        "ts": pa.array([dt.datetime(2024, 1, 1, tzinfo=utc), dt.datetime(2023, 1, 1, tzinfo=utc), None, None],
                       pa.timestamp("us", tz="UTC")),
        "dec": pa.array([decimal.Decimal("1.50"), decimal.Decimal("-7.25"), None, decimal.Decimal("3")],
                        pa.decimal128(6, 2)),
        "uuid": pa.array([u.bytes] * 4, pa.binary(16)),
    })
    path = tmp_path / "stats.parquet"
    carquet.write_table(t, path)
    rg = carquet.read_metadata(path).row_group(0)
    st = {c.path_in_schema: c.statistics for c in rg}
    assert (st["i"].min, st["i"].max, st["i"].null_count) == (-3, 9, 1)
    assert (st["u"].min, st["u"].max) == (0, 2**32 - 1)
    assert (st["f"].min, st["f"].max) == (-2.25, 3.0)
    assert (st["s"].min, st["s"].max) == ("apple", "zoo")
    assert (st["b"].min, st["b"].max) == (False, True)
    assert (st["d"].min, st["d"].max) == (dt.date(2019, 5, 6), dt.date(2021, 1, 1))
    assert (st["ts"].min, st["ts"].max) == (dt.datetime(2023, 1, 1, tzinfo=utc), dt.datetime(2024, 1, 1, tzinfo=utc))
    assert (st["dec"].min, st["dec"].max) == (decimal.Decimal("-7.25"), decimal.Decimal("3.00"))
    assert st["uuid"].min == u.bytes
    assert st["i"].has_null_count and st["i"].physical_type == "INT32"
    assert st["i"].to_dict()["min"] == -3


def test_parquet_schema(tmp_path, typed_table):
    path = tmp_path / "s.parquet"
    carquet.write_table(typed_table, path)
    ps = carquet.read_parquet_schema(path)
    assert len(ps) == len(ps.columns) == carquet.info(path).num_columns
    assert ps["dec"].logical_type.name == "DECIMAL"
    assert (ps["dec"].precision, ps["dec"].scale) == (9, 2)
    assert str(ps["dec"].logical_type) == "DECIMAL(9, 2)"
    assert ps["ts_tz"].logical_type.unit == "MILLIS" and ps["ts_tz"].logical_type.is_adjusted_to_utc
    assert ps["i8"].logical_type.bit_width == 8 and ps["i8"].logical_type.is_signed
    assert ps["u64"].logical_type.is_signed is False
    assert ps["f16"].logical_type.name == "FLOAT16" and ps["f16"].type_length == 2
    assert ps["bin"].logical_type is None and ps["bin"].physical_type == "BYTE_ARRAY"
    m_leaves = ps.leaves_of("m")
    assert [c.path for c in m_leaves] == ["m.key_value.key", "m.key_value.value"]
    assert all(c.max_repetition_level == 1 for c in m_leaves)
    nested = ps.column("nested.list.element.y.list.element")
    assert nested.max_repetition_level == 2 and nested.is_nested
    assert ps.fields == typed_table.column_names
    kinds = {e.name: e for e in ps.elements if e.depth == 1}
    assert kinds["m"].logical_type.name == "MAP" and not kinds["m"].is_leaf
    assert "x" in ps and "nope" not in ps
    assert ps["element"].path == "nested.list.element.y.list.element"  # unique leaf name
    with pytest.raises(KeyError):
        ps["nope"]
    with pytest.raises(IndexError):
        ps[999]
    assert "dec: OPTIONAL FIXED_LEN_BYTE_ARRAY" in repr(ps)
    assert ps == carquet.read_parquet_schema(path)


def test_page_indexes(tmp_path, table):
    path = tmp_path / "pi.parquet"
    carquet.write_table(table, path, write_page_index=True, max_rows_per_page=100)
    with carquet.ParquetFile(path) as pf:
        ci = pf.column_index("id")
        oi = pf.offset_index("id")
        assert len(ci) == len(oi) == 10
        assert ci.boundary_order in ("UNORDERED", "ASCENDING")  # the C writer records UNORDERED
        assert [(p.min, p.max) for p in ci][:2] == [(0, 99), (100, 199)]
        assert ci[0].null_count == 0 and not ci[0].is_null_page
        assert [p.first_row_index for p in oi] == list(range(0, 1000, 100))
        assert oi[1].offset > oi[0].offset and oi[0].compressed_size > 0
        assert pf.column_index("name")[0].min.startswith("name-")
        assert "num_pages=10" in repr(ci)
    carquet.write_table(table, path)
    with carquet.ParquetFile(path) as pf:
        assert pf.column_index("id") is None and pf.offset_index("id") is None
        assert pf.metadata.row_group(0).column("id").has_column_index is False


def test_metadata_survives_close(written):
    pf = carquet.ParquetFile(written)
    md = pf.metadata
    ps = pf.parquet_schema
    pf.close()
    assert md.row_group(0).column("id").statistics.max == 999
    assert ps["id"].physical_type == "INT64"


def test_metadata_matches_pyarrow(tmp_path, table):
    path = tmp_path / "x.parquet"
    pq.write_table(table, path, row_group_size=250, compression="snappy")
    ours = carquet.read_metadata(path)
    theirs = pq.ParquetFile(path).metadata
    assert ours.num_rows == theirs.num_rows
    assert ours.num_row_groups == theirs.num_row_groups
    for i in range(theirs.num_row_groups):
        for j in range(theirs.num_columns):
            a, b = ours.row_group(i).column(j), theirs.row_group(i).column(j)
            assert a.path_in_schema == b.path_in_schema
            assert a.compression == b.compression
            assert a.total_compressed_size == b.total_compressed_size
            assert set(a.encodings) == set(b.encodings)
            if b.is_stats_set and b.statistics.has_min_max:
                assert a.statistics.min == b.statistics.min
                assert a.statistics.max == b.statistics.max
                assert a.statistics.null_count == b.statistics.null_count
