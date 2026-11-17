"""Cross-library compatibility: what DuckDB / polars / PyArrow / pandas need from
carquet's files, and what carquet must give back for theirs."""

from __future__ import annotations

import decimal
import uuid

import pytest

import carquet

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")


def test_decimal_scale_zero_is_readable_elsewhere(tmp_path):
    """DECIMAL(p, 0) must carry `scale` next to `precision` in the schema element:
    DuckDB, polars and fastparquet refuse the file otherwise."""
    path = tmp_path / "dec0.parquet"
    t = pa.table({"d": pa.array([decimal.Decimal(12345), None, decimal.Decimal(-7)], pa.decimal128(5, 0))})
    carquet.write_table(t, path)
    assert pq.read_table(path).equals(t)
    try:
        import duckdb
    except ImportError:
        duckdb = None
    if duckdb is not None:
        rows = duckdb.sql(f"select d from read_parquet('{path}')").fetchall()
        assert rows == [(decimal.Decimal(12345),), (None,), (decimal.Decimal(-7),)]
    try:
        import polars
    except ImportError:
        polars = None
    if polars is not None:
        assert polars.read_parquet(path)["d"].to_list() == [decimal.Decimal(12345), None, decimal.Decimal(-7)]
    if duckdb is None and polars is None:
        # the footer bytes themselves: converted_type DECIMAL, scale 0, precision 5
        assert b"\x25\x0a\x15\x00\x15\x0a" in path.read_bytes()


def test_uuid_and_json_extensions_get_logical_types(tmp_path):
    if not hasattr(pa, "uuid"):
        pytest.skip("pyarrow without the uuid extension type")
    path = tmp_path / "ext.parquet"
    ids = [uuid.UUID(int=i * 7919 + 1).bytes for i in range(5)]
    cols = {"u": pa.array(ids, pa.uuid()), "plain": pa.array(ids, pa.binary(16))}
    if hasattr(pa, "json_"):
        cols["j"] = pa.array(['{"a": 1}'] * 5, pa.json_())
    t = pa.table(cols)
    carquet.write_table(t, path)
    schema = carquet.read_parquet_schema(path)
    assert schema.column("u").logical_type.name == "UUID"
    assert schema.column("plain").logical_type is None
    if "j" in cols:
        assert schema.column("j").logical_type.name == "JSON"
    assert pq.read_table(path).equals(t)   # ARROW:schema restores the extension types
    assert pq.ParquetFile(path).schema.column(0).logical_type.type == "UUID"


def test_empty_string_minimum(tmp_path):
    path = tmp_path / "empty_min.parquet"
    pq.write_table(pa.table({"s": ["", "b", None], "b": [b"", b"\x01", None]}), path)
    rg = carquet.read_metadata(path).row_group(0)
    s = rg.column("s").statistics
    assert (s.min, s.max, s.min_raw) == ("", "b", b"")
    assert rg.column("b").statistics.min == b""


def test_window_of_nested_file_imports_everywhere(tmp_path, table):
    """A row window on the row-group path (nested columns) is a batch with an
    offset; PyArrow refuses a record batch whose struct carries one."""
    path = tmp_path / "nested.parquet"
    pq.write_table(table, path, row_group_size=400)
    with carquet.ParquetFile(path) as pf:
        for offset, length in ((333, 100), (0, 50), (390, 30), (999, 10)):
            got = pf.slice(offset, length)
            assert pa.table(got).to_pylist() == table.slice(offset, length).to_pylist()
            for b in got.to_batches():
                assert pa.record_batch(b).num_rows == b.num_rows


def test_required_field_in_null_struct_has_no_nulls(tmp_path):
    """A non-nullable child of a null struct must not come back with nulls:
    pyarrow.parquet refuses to write such a table."""
    ty = pa.struct([pa.field("m", pa.binary(), nullable=False), pa.field("v", pa.int32()),
                    pa.field("inner", pa.struct([pa.field("q", pa.int64(), nullable=False)]), nullable=False)])
    t = pa.table({"var": pa.array([None, {"m": b"x", "v": None, "inner": {"q": 1}},
                                   {"m": b"y", "v": 2, "inner": {"q": 2}}], ty)})
    path = tmp_path / "req_in_opt.parquet"
    pq.write_table(t, path)
    got = pa.table(carquet.read_table(path))
    assert got.equals(t)
    var = got.column("var").chunk(0)
    assert var.null_count == 1
    assert var.field("m").null_count == 0 and var.field("inner").field("q").null_count == 0
    pq.write_table(got, tmp_path / "again.parquet")
    assert pq.read_table(tmp_path / "again.parquet").equals(t)


def test_all_null_string_column_to_pandas(tmp_path):
    pd = pytest.importorskip("pandas")
    path = tmp_path / "allnull.parquet"
    t = pa.table({"s": pa.array([None, None], pa.string()), "x": pa.array(["a", None], pa.string())})
    pq.write_table(t, path)
    got = carquet.read_table(path).to_pandas()
    expected = pq.read_table(path).to_pandas()
    assert got["s"].dtype == expected["s"].dtype
    assert got["x"].dtype == expected["x"].dtype
    pd.testing.assert_frame_equal(got, expected)


def test_date32_to_numpy_and_pandas(tmp_path):
    """date32 is 32-bit: viewing its buffer as datetime64[D] (64-bit) read past
    the buffer and returned garbage dates."""
    import datetime as dt

    np = pytest.importorskip("numpy")
    days = [dt.date(2020, 1, 1) + dt.timedelta(days=i) for i in range(500)]
    days[7] = None
    path = tmp_path / "dates.parquet"
    pq.write_table(pa.table({"d": pa.array(days, pa.date32())}), path)
    t = carquet.read_table(path)
    arr = t.column("d").to_numpy()
    assert arr.dtype == np.dtype("datetime64[D]") and len(arr) == 500
    assert np.isnat(arr[7]) and arr[0] == np.datetime64("2020-01-01") and arr[499] == np.datetime64(days[499])
    pd = pytest.importorskip("pandas")
    pd.testing.assert_frame_equal(t.to_pandas(), pq.read_table(path).to_pandas())


def test_pandas_dtypes_follow_pyarrow(tmp_path):
    """Timestamp units, tz-aware columns, list cells and the pandas metadata:
    carquet's read of a file and PyArrow's give the same DataFrame, and a
    DataFrame written by carquet reads back the same through PyArrow."""
    import datetime as dt

    pd = pytest.importorskip("pandas")
    np = pytest.importorskip("numpy")
    n = 60
    df = pd.DataFrame({
        "ts": pd.date_range("2024-01-01", periods=n, freq="37min"),
        "ts_ns": pd.date_range("2024-01-01", periods=n, freq="37min").astype("datetime64[ns]"),
        "ts_tz": pd.date_range("2024-01-01", periods=n, freq="h", tz="America/New_York"),
        "td": pd.to_timedelta(np.arange(n), unit="s"),
        "date": [dt.date(2020, 1, 1) + dt.timedelta(days=i) for i in range(n)],
        "dec": [decimal.Decimal(i) / 100 for i in range(n)],
        "lst": [[i, i + 1] if i % 3 else None for i in range(n)],
        "s": [None if i % 9 == 0 else f"s{i}" for i in range(n)],
    })
    ref = tmp_path / "pyarrow.parquet"
    df.to_parquet(ref, engine="pyarrow")
    expected = pd.read_parquet(ref, engine="pyarrow")
    pd.testing.assert_frame_equal(carquet.read_pandas(ref), expected)

    ours = tmp_path / "carquet.parquet"
    carquet.write_table(df, ours)
    pd.testing.assert_frame_equal(pd.read_parquet(ours, engine="pyarrow"), expected)
    pd.testing.assert_frame_equal(carquet.read_pandas(ours), expected)
    a, b = pq.read_schema(ours), pq.read_schema(ref)
    for name in ("ts", "ts_ns", "ts_tz", "td", "date", "dec"):
        assert a.field(name).type == b.field(name).type, name


def test_truncated_string_statistics_stay_utf8(tmp_path):
    """A truncated STRING max used to be cut inside a character (then byte-incremented),
    so PyArrow raised UnicodeDecodeError on `statistics.max`."""
    path = tmp_path / "utf8stats.parquet"
    values = ["a", "ÿ" * 40, "￿" * 20, "z" * 31 + "é" * 10]
    carquet.write_table(pa.table({"s": values}), path)
    st = pq.read_metadata(path).row_group(0).column(0).statistics
    assert st.has_min_max
    assert st.min <= min(values) and st.max >= max(values)
    assert len(st.max.encode()) <= 32
    assert carquet.read_table(path, filters=[("s", "==", "￿" * 20)]).num_rows == 1
