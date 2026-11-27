"""Reading: ParquetFile, read_table, projections, sources, streaming, fidelity."""

from __future__ import annotations

import io

import pytest

import carquet
from carquet import CarquetError, StatusCode
from conftest import assert_tables_equal

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")


def test_roundtrip_exact(tmp_path, table):
    path = tmp_path / "t.parquet"
    carquet.write_table(table, path)
    out = carquet.read_table(path)
    assert out.equals(table)
    assert pq.read_table(path).equals(table)


def test_typed_roundtrip_exact(tmp_path, typed_table):
    """Every type written by carquet comes back identical through carquet and PyArrow."""
    path = tmp_path / "typed.parquet"
    carquet.write_table(typed_table, path, row_group_size=64)
    out = carquet.read_table(path)
    for name in typed_table.column_names:
        assert out.schema.field(name).type == typed_table.schema.field(name).type, name
    assert_tables_equal(out, typed_table)
    assert_tables_equal(pq.read_table(path), typed_table)


@pytest.mark.parametrize("compression", ["NONE", "SNAPPY", "ZSTD", "GZIP", "LZ4", "BROTLI"])
def test_reads_pyarrow_files(tmp_path, typed_table, compression):
    path = tmp_path / f"pa_{compression}.parquet"
    pq.write_table(typed_table, path, compression=compression, row_group_size=50)
    expected = pq.read_table(path)
    out = carquet.read_table(path)
    assert out.schema.equals(expected.schema)
    assert_tables_equal(out, expected)


def test_reads_pyarrow_v2_pages_and_page_index(tmp_path, table):
    path = tmp_path / "v2.parquet"
    pq.write_table(table, path, data_page_version="2.0", write_page_index=True, use_dictionary=False)
    assert carquet.read_table(path).equals(pq.read_table(path))


def test_int96_timestamps(tmp_path):
    import datetime as dt

    t = pa.table({"ts": pa.array([dt.datetime(2020, 5, 17, 1, 2, 3, 456789), None, dt.datetime(1999, 12, 31)],
                                 pa.timestamp("ns"))})
    path = tmp_path / "int96.parquet"
    pq.write_table(t, path, use_deprecated_int96_timestamps=True, store_schema=False)
    assert carquet.read_parquet_schema(path)["ts"].physical_type == "INT96"
    out = carquet.read_table(path)
    assert out.schema.field("ts").type == pa.timestamp("ns")
    assert out.equals(pq.read_table(path))


def test_integer_backed_decimals(tmp_path):
    import decimal

    t = pa.table({"d32": pa.array([decimal.Decimal("1.23"), None, decimal.Decimal("-9.99")], pa.decimal128(5, 2)),
                  "d64": pa.array([decimal.Decimal("123456789.1234"), decimal.Decimal("-1"), None],
                                  pa.decimal128(15, 4))})
    path = tmp_path / "decint.parquet"
    try:
        pq.write_table(t, path, store_decimal_as_integer=True)
    except TypeError:
        pytest.skip("pyarrow without store_decimal_as_integer")
    phys = {c.name: c.physical_type for c in carquet.read_parquet_schema(path)}
    out = carquet.read_table(path)
    assert out.equals(t), phys


def test_columns_projection(written, table):
    out = carquet.read_table(written, columns=["name", "id"])
    assert out.column_names == ["name", "id"]
    assert out.column("id").equals(table.column("id"))
    # nested fields go through the nested assembler with field projection
    out = carquet.read_table(written, columns=["point", "id"])
    assert out.column_names == ["point", "id"]
    assert out.column("point").equals(table.column("point"))
    # a dotted path selects its top-level field
    assert carquet.read_table(written, columns=["point.b"]).column_names == ["point"]
    assert carquet.read_table(written, columns=["tags"]).column("tags").equals(table.column("tags"))


def test_unknown_column(written):
    with pytest.raises(KeyError):
        carquet.read_table(written, columns=["nope"])
    with pytest.raises(TypeError):
        carquet.read_table(written, columns="id")


def test_row_groups(tmp_path, table):
    path = tmp_path / "rg.parquet"
    carquet.write_table(table, path, row_group_size=300)
    with carquet.ParquetFile(path) as pf:
        assert pf.num_row_groups == 4
        assert pf.read_row_group(1).column("id").to_pylist() == list(range(300, 600))
        both = pf.read_row_groups([3, 0], columns=["id", "point"])
        assert both.column("id").to_pylist() == list(range(900, 1000)) + list(range(0, 300))
        assert [t.num_rows for t in pf.iter_row_groups(["id"])] == [300, 300, 300, 100]
        assert pf.read(row_groups=[2]).num_rows == 300
        with pytest.raises(IndexError):
            pf.read_row_group(7)
        assert pf.read(row_groups=[]).num_rows == 0


def test_iter_batches(tmp_path, table):
    path = tmp_path / "it.parquet"
    carquet.write_table(table, path, row_group_size=400)
    with carquet.ParquetFile(path) as pf:
        sizes = [b.num_rows for b in pf.iter_batches(batch_size=150)]
        assert sum(sizes) == table.num_rows and max(sizes) <= 150
        flat = pa.Table.from_batches([pa.record_batch(b) for b in pf.iter_batches(batch_size=150, columns=["id", "name"])])
        assert flat.equals(table.select(["id", "name"]))
        nested = pa.Table.from_batches([pa.record_batch(b) for b in pf.iter_batches(batch_size=128, columns=["point"])])
        assert nested.equals(table.select(["point"]))
        with pytest.raises(ValueError):
            next(pf.iter_batches(batch_size=0))
    # module-level helper closes the file itself
    assert sum(b.num_rows for b in carquet.iter_batches(path, 333)) == table.num_rows


def test_slice_and_head(tmp_path, table):
    path = tmp_path / "s.parquet"
    carquet.write_table(table, path, row_group_size=256)
    with carquet.ParquetFile(path) as pf:
        assert pf.slice(250, 20).equals(table.slice(250, 20))  # spans a row-group boundary, nested columns
        assert pf.slice(250, 20, columns=["id", "name"]).equals(table.select(["id", "name"]).slice(250, 20))
        assert pf.slice(990).num_rows == 10
        assert pf.slice(-5).column("id").to_pylist() == list(range(995, 1000))
        assert pf.slice(5000).num_rows == 0
        assert pf.head(3).equals(table.slice(0, 3))
        assert pf.head(3, columns=["id"]).column("id").to_pylist() == [0, 1, 2]


@pytest.mark.parametrize("kind", ["bytes", "bytearray", "memoryview", "buffer", "bytesio", "file"])
def test_buffer_sources(tmp_path, table, kind):
    path = tmp_path / "b.parquet"
    carquet.write_table(table, path)
    data = path.read_bytes()
    if kind == "bytes":
        src = data
    elif kind == "bytearray":
        src = bytearray(data)
    elif kind == "memoryview":
        src = memoryview(data)
    elif kind == "buffer":
        src = pa.py_buffer(data)
    elif kind == "bytesio":
        src = io.BytesIO(data)
    else:
        src = open(path, "rb")
    try:
        assert carquet.read_table(src).equals(table)
    finally:
        if kind == "file":
            src.close()
    with carquet.ParquetFile(data) as pf:
        assert pf.path is None
        assert pf.metadata.file_size == len(data)
        assert "<buffer>" in repr(pf)


def test_directory_and_list_sources(tmp_path, table):
    d = tmp_path / "ds"
    (d / "sub").mkdir(parents=True)
    carquet.write_table(table.slice(0, 400), d / "a.parquet")
    carquet.write_table(table.slice(400), d / "sub" / "b.parquet")
    (d / "_SUCCESS").write_text("")
    assert carquet.read_table(d).equals(table)
    assert carquet.read_table([d / "a.parquet", d / "sub" / "b.parquet"], columns=["id"]).num_rows == 1000
    with pytest.raises(ValueError):
        carquet.read_table(tmp_path / "ds" / "sub" / "..", row_groups=[0])


def test_buffered_io_and_threads(written, table):
    assert carquet.read_table(written, memory_map=False).equals(table)
    assert carquet.read_table(written, use_threads=False).equals(table)
    assert carquet.read_table(written, use_mmap=False).equals(table)  # 0.7.2 spelling
    with carquet.ParquetFile(written, num_threads=2, verify_checksums=False) as pf:
        assert pf.read().equals(table)


def test_empty_and_zero_column_selection(tmp_path, table):
    path = tmp_path / "empty.parquet"
    carquet.write_table(table.schema.empty_table(), path)
    out = carquet.read_table(path)
    assert out.num_rows == 0 and out.schema.names == table.schema.names
    assert carquet.read_table(path, columns=["id"]).num_rows == 0


def test_closed_file(written):
    pf = carquet.ParquetFile(written)
    it = pf.iter_batches(batch_size=100, columns=["id", "name"])
    next(it)
    pf.close()
    assert pf.closed
    with pytest.raises(CarquetError) as ei:  # cut short, not silently finished
        next(it)
    assert ei.value.status == StatusCode.ALREADY_CLOSED
    with pytest.raises(CarquetError) as ei:
        pf.read()
    assert ei.value.status == StatusCode.ALREADY_CLOSED
    pf.close()  # idempotent
    assert "closed" in repr(pf)


def test_restore_arrow_schema_toggle(tmp_path, typed_table):
    path = tmp_path / "r.parquet"
    carquet.write_table(typed_table.select(["ts_tz", "cat", "lstr"]), path)
    raw = carquet.read_table(path, restore_arrow_schema=False)
    assert raw.schema.field("ts_tz").type == pa.timestamp("ms", tz="UTC")
    assert raw.schema.field("cat").type == pa.string()
    assert raw.schema.field("lstr").type == pa.string()
    restored = carquet.read_table(path)
    assert restored.schema.field("cat").type == pa.dictionary(pa.int32(), pa.string())
    assert restored.schema.field("lstr").type == pa.large_string()


def test_schema_metadata_and_pandas(tmp_path):
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame({"a": [1, 2, 3], "b": ["x", None, "z"]}, index=pd.Index([10, 20, 30], name="idx"))
    path = tmp_path / "pd.parquet"
    carquet.write_table(df, path)
    back = carquet.read_pandas(path)
    pd.testing.assert_frame_equal(back, df)
    pd.testing.assert_frame_equal(pq.read_pandas(path).to_pandas(), df)
    assert b"pandas" in carquet.read_table(path).schema.metadata


def test_read_schema_helpers(written, table):
    s = carquet.read_schema(written)
    assert s.equals(table.schema, check_metadata=False)
    assert carquet.schema(written).names == table.schema.names  # 0.7.2 alias
    ps = carquet.read_parquet_schema(written)
    assert ps.names == ["id", "value", "name", "ts", "tags.list.element", "point.a", "point.b"]


def test_reader_alias_and_open_reader(written, table):
    with carquet.open_reader(written, use_mmap=False) as r:
        assert isinstance(r, carquet.Reader)
        assert r.read_all(columns=["id"]).column("id").equals(table.column("id"))
        assert r.schema_arrow.names == table.schema.names


def test_count_rows(tmp_path, table):
    path = tmp_path / "c.parquet"
    carquet.write_table(table, path, row_group_size=100)
    with carquet.ParquetFile(path) as pf:
        assert pf.count_rows() == 1000
        assert pf.count_rows([("id", "<", 250)]) == 250
        assert pf.count_rows([("name", "is null")]) == table.column("name").null_count


def test_threaded_readers(tmp_path, table):
    """Separate ParquetFile objects can be used from several threads at once (the GIL is released)."""
    from concurrent.futures import ThreadPoolExecutor

    path = tmp_path / "mt.parquet"
    carquet.write_table(table, path, row_group_size=200)

    def work(i):
        with carquet.ParquetFile(path) as pf:
            return pf.read(row_groups=[i % 5]).num_rows

    with ThreadPoolExecutor(4) as ex:
        assert sum(ex.map(work, range(20))) == 20 * 200


def test_slice_with_list_and_map(tmp_path):
    t = pa.table({"id": list(range(50)), "l": pa.array([[i] * (i % 3) for i in range(50)], pa.list_(pa.int32())),
                  "m": pa.array([[("k", i)] for i in range(50)], pa.map_(pa.string(), pa.int64()))})
    path = tmp_path / "lm.parquet"
    carquet.write_table(t, path, row_group_size=16)
    with carquet.ParquetFile(path) as pf:
        assert pf.head(3).equals(t.slice(0, 3))
        assert pf.slice(10, 20).equals(t.slice(10, 20))
        assert pf.slice(40).equals(t.slice(40))


def test_sliced_nested_conversions(tmp_path):
    import decimal

    t = pa.table({
        "s": pa.array([{"d": decimal.Decimal(i), "k": i} if i % 5 else None for i in range(20)],
                      pa.struct([("d", pa.decimal128(9, 0)), ("k", pa.int64())])),
        "m": pa.array([[("a", decimal.Decimal(i).scaleb(-1))] for i in range(20)],
                      pa.map_(pa.string(), pa.decimal128(6, 1))),
        "l": pa.array([[decimal.Decimal(i)] for i in range(20)], pa.list_(pa.decimal128(4, 0))),
    })
    path = tmp_path / "sn.parquet"
    carquet.write_table(t, path)
    with carquet.ParquetFile(path) as pf:
        assert pa.Table.from_batches([pa.record_batch(b) for b in pf.iter_batches(batch_size=7)]).equals(t)
        assert pf.slice(5, 5).equals(t.slice(5, 5))


def test_zero_columns(tmp_path, table):
    path = tmp_path / "z.parquet"
    carquet.write_table(table, path, row_group_size=300)
    out = carquet.read_table(path, columns=[])
    assert (out.num_rows, out.num_columns) == (1000, 0)
    assert carquet.read_table(path, columns=[], filters=[("id", "<", 10)]).num_rows == 10
    with carquet.ParquetFile(path) as pf:
        assert sum(b.num_rows for b in pf.iter_batches(256, columns=[])) == 1000
        assert pf.slice(10, 5, columns=[]).num_rows == 5
