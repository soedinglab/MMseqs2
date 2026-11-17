"""The package without PyArrow: writing Python / NumPy / pandas data, reading
carquet Tables, filters, streams. ``test_no_pyarrow`` runs this module with
PyArrow made unimportable."""

from __future__ import annotations

import datetime as dt
import decimal
import math

import pytest

import carquet as cq

ROWS = 1000


def _data(n: int = ROWS) -> dict:
    utc = dt.timezone.utc
    return {
        "id": list(range(n)),
        "value": [i * 0.5 for i in range(n)],
        "name": [None if i % 7 == 3 else f"name-{i}" for i in range(n)],
        "flag": [i % 3 == 0 for i in range(n)],
        "ts": [dt.datetime(2024, 1, 1) + dt.timedelta(seconds=i) for i in range(n)],
        "ts_utc": [dt.datetime(2024, 1, 1, tzinfo=utc) + dt.timedelta(minutes=i) for i in range(n)],
        "day": [dt.date(2020, 1, 1) + dt.timedelta(days=i) for i in range(n)],
        "amount": [None if i % 11 == 0 else decimal.Decimal(i * 25 - 3000).scaleb(-2) for i in range(n)],
        "blob": [bytes([i % 256]) * (i % 4) for i in range(n)],
        "tags": [None if i % 13 == 5 else [j for j in range(i % 4)] for i in range(n)],
        "point": [{"a": i % 100, "b": None if i % 9 == 0 else f"s{i}"} for i in range(n)],
    }


@pytest.fixture(scope="module")
def data() -> dict:
    return _data()


@pytest.fixture(scope="module")
def path(tmp_path_factory, data):
    p = tmp_path_factory.mktemp("core") / "core.parquet"
    cq.write_table(data, p, row_group_size=300, write_page_index=True, max_rows_per_page=50)
    return p


def test_roundtrip_python_values(path, data):
    t = cq.read_table(path)
    assert isinstance(t, cq.Table)
    assert t.num_rows == ROWS and t.column_names == list(data)
    assert t.to_pydict() == data
    types = {f.name: str(f.type) for f in t.schema}
    assert types["id"] == "int64" and types["name"] == "string" and types["ts"] == "timestamp[us]"
    assert types["ts_utc"] == "timestamp[us, tz=UTC]" and types["amount"] == "decimal128(5, 2)"
    assert types["tags"].startswith("list<") and types["point"] == "struct<a: int64, b: string>"
    assert [b.num_rows for b in t.to_batches()] == [300, 300, 300, 100]


def test_projection_and_values(path):
    t = cq.read_table(path, columns=["name", "id"])
    assert t.column_names == ["name", "id"]
    assert t.column("id")[5] == 5 and t["name"][3] is None
    assert t.slice(10, 3).column("id").to_pylist() == [10, 11, 12]
    assert t.select(["id"]).column_names == ["id"]
    assert t.to_pylist()[2] == {"name": "name-2", "id": 2}


def test_filters_exact(path, data):
    def ids(filters, **kw):
        return cq.read_table(path, columns=["id"], filters=filters, **kw).column("id").to_pylist()

    assert ids([("id", ">=", 995)]) == [995, 996, 997, 998, 999]
    assert ids([("name", "is null"), ("id", "<", 30)]) == [3, 10, 17, 24]
    assert ids([[("id", "<", 2)], [("id", ">", 997)]]) == [0, 1, 998, 999]
    assert ids([("point.a", "==", 42), ("id", "<", 300)]) == [42, 142, 242]
    assert ids([("amount", "<", decimal.Decimal("-29.5"))]) == [1]
    assert ids([("id", "<", 2.5)]) == [0, 1, 2]                    # inexact value, exact result
    assert ids([("value", "==", 0.25)]) == []
    assert ids([("day", ">=", dt.date(2022, 9, 25))]) == [998, 999]
    assert ids([("ts_utc", "<", dt.datetime(2024, 1, 1, 0, 2, tzinfo=dt.timezone.utc))]) == [0, 1]
    assert ids([("name", "in", ["name-1", "name-2", None]), ("id", "<", 4)]) == [1, 2, 3]
    assert ids([("id", "in", [])]) == []
    with cq.ParquetFile(path) as pf:
        assert pf.filter_row_groups([("id", ">=", 950)]) == [3]
        assert pf.count_rows([("flag", "==", True)]) == sum(data["flag"])
        assert pf.count_rows() == ROWS
    with pytest.raises(ValueError):
        ids([("tags", "==", 1)])


def test_batches_slices_and_streams(path, data):
    with cq.ParquetFile(path) as pf:
        sizes = [b.num_rows for b in pf.iter_batches(128)]
        assert sum(sizes) == ROWS and max(sizes) <= 128
        assert pf.slice(298, 4, ["id"]).column("id").to_pylist() == [298, 299, 300, 301]
        assert pf.head(2).to_pydict()["name"] == ["name-0", "name-1"]
        stream = pf.arrow_stream(["id", "point"], filters=[("id", "<", 5)])
    again = cq.table(stream)                      # a carquet stream read back by carquet
    assert again.to_pydict() == {"id": [0, 1, 2, 3, 4], "point": data["point"][:5]}
    copy = cq.table(cq.read_table(path))          # PyCapsule round trip
    assert copy.equals(cq.read_table(path))


def test_numpy(tmp_path):
    np = pytest.importorskip("numpy")
    a = np.arange(10, dtype=np.int32)
    f = np.linspace(0, 1, 10)
    t64 = np.array(["2020-01-01T00:00:00.000000001", "NaT"] * 5, dtype="datetime64[ns]")
    p = tmp_path / "np.parquet"
    cq.write_table({"a": a, "f": f, "t": t64}, p)
    t = cq.read_table(p)
    col = t.column("a").to_numpy()
    assert col.dtype == np.int32 and col.tolist() == list(range(10))
    one = t.column("f").chunk(0).to_numpy(zero_copy_only=True)   # no copy
    assert one.base is not None and math.isclose(one[-1], 1.0)
    ts = t.column("t").to_numpy()
    assert ts.dtype == np.dtype("datetime64[ns]") and np.isnat(ts[1]) and ts[0] == t64[0]


def test_pandas_roundtrip(tmp_path):
    pd = pytest.importorskip("pandas")
    df = pd.DataFrame({
        "a": [1, 2, 3],
        "b": ["x", None, "z"],
        "c": pd.Categorical(["u", "v", "u"]),
        "d": pd.array([1, None, 3], dtype="Int64"),
        "e": pd.date_range("2020-01-01", periods=3, tz="Europe/Paris"),
        "f": [1.0, float("nan"), 3.0],
    }, index=pd.Index([10, 20, 30], name="k"))
    p = tmp_path / "pd.parquet"
    cq.write_table(df, p)
    back = cq.read_pandas(p)
    assert list(back.index) == [10, 20, 30] and back.index.name == "k"
    assert back["c"].dtype.name == "category" and str(back["d"].dtype) == "Int64"
    assert str(back["e"].dt.tz) == "Europe/Paris" and math.isnan(back["f"][20])
    assert back["b"].tolist()[0] == "x"
    assert isinstance(cq.read_pandas(p, columns=["a"]).index, pd.Index)


def test_writer_options_and_metadata(tmp_path, data):
    p = tmp_path / "opts.parquet"
    with cq.ParquetWriter(p, cq.table(data).schema, compression="snappy", row_group_size=400,
                          use_dictionary=["name"], bloom_filters=["id"], metadata={"k": "v"}) as w:
        w.write(data)
    md = cq.read_metadata(p)
    assert md.num_rows == ROWS and md.num_row_groups == 3 and md.metadata["k"] == "v"
    assert md.row_group(0).column("id").compression == "SNAPPY"
    with cq.ParquetFile(p) as pf:
        assert 5 in pf.bloom_filter("id") and pf.read().to_pydict() == data
    with pytest.raises(ValueError):
        with cq.ParquetWriter(tmp_path / "bad.parquet", cq.table(data).schema) as w:
            w.write({"other": [1]})


def test_widening_and_nullability(tmp_path):
    p = tmp_path / "w.parquet"
    schema = cq.table({"x": [1]}).schema                       # int64
    with cq.ParquetWriter(p, schema) as w:
        w.write({"x": __import__("array").array("i", [1, 2]).tolist()})
    np = pytest.importorskip("numpy")
    with cq.ParquetWriter(p, schema) as w:
        w.write({"x": np.array([3, 4], dtype=np.int32)})       # int32 widened into int64
    assert cq.read_table(p).column("x").to_pylist() == [3, 4]


def test_table_from_arrow_arrays_keeps_names_and_metadata():
    """Nested / dictionary Arrow arrays in a mapping take the mapping's names
    (they kept the producer's, usually empty), and a schema that carries
    metadata is accepted (it failed with a ctypes ArgumentError)."""
    pa = pytest.importorskip("pyarrow")
    t = cq.table({
        "lst": pa.array([[1, 2], None]),
        "dic": pa.array(["a", "b"]).dictionary_encode(),
        "num": pa.array([1, 2]),
    })
    assert t.column_names == ["lst", "dic", "num"]
    assert t.to_pylist() == [{"lst": [1, 2], "dic": "a", "num": 1}, {"lst": None, "dic": "b", "num": 2}]

    schema = pa.schema([pa.field("a", pa.int64())], metadata={"k": "v"})
    rb = pa.RecordBatch.from_pylist([{"a": 1}, {"a": 2}], schema=schema)
    t = cq.table({"s": rb})
    assert t.column_names == ["s"]
    assert t.to_pylist() == [{"s": {"a": 1}}, {"s": {"a": 2}}]


def test_dictionary_null_slot_with_out_of_range_index():
    """The index under a null slot is unspecified: it must not be looked up."""
    pa = pytest.importorskip("pyarrow")
    indices = pa.array([0, 99, 1], pa.int8())
    d = pa.DictionaryArray.from_buffers(
        pa.dictionary(pa.int8(), pa.string()), 3,
        [pa.py_buffer(bytes([0b101])), indices.buffers()[1]], pa.array(["a", "b"]))
    t = cq.table({"d": d})
    assert t.column("d").to_pylist() == ["a", None, "b"]
    pd = pytest.importorskip("pandas")
    assert t.to_pandas()["d"].astype(object).where(lambda s: s.notna(), None).tolist() == ["a", None, "b"]


def test_help_covers_the_public_api():
    """``help(carquet)`` stays complete: every public function, class, method
    and property has a docstring."""
    import inspect

    missing = []
    assert inspect.getdoc(cq)
    for name in cq.__all__:
        obj = getattr(cq, name)
        if inspect.isclass(obj):
            if not inspect.getdoc(obj):
                missing.append(name)
            for attr, value in vars(obj).items():
                if attr.startswith("_"):
                    continue
                if isinstance(value, property):
                    value = value.fget
                elif isinstance(value, (classmethod, staticmethod)):
                    value = value.__func__
                elif not inspect.isfunction(value):
                    continue
                if not inspect.getdoc(value):
                    missing.append(f"{name}.{attr}")
        elif callable(obj) and not inspect.getdoc(obj):
            missing.append(name)
    assert not missing, f"no docstring: {sorted(set(missing))}"
