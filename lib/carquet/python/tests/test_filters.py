"""Row filters: DNF semantics, statistics / page-index pushdown, exactness."""

from __future__ import annotations

import datetime as dt
import decimal

import pytest

import carquet

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")
pc = pytest.importorskip("pyarrow.compute")


def to_expression(filters, top_level):
    """The DNF filter as a pyarrow.compute expression: the oracle carquet's
    C predicate is checked against."""
    dnf = carquet.filters.normalize_filters(filters).dnf
    tops = set(top_level)

    def ref(name):
        return pc.field(name) if name in tops or "." not in name else pc.field(*name.split("."))

    disj = None
    for conj in dnf:
        e = None
        for col, op, val in conj:
            f = ref(col)
            c = {"==": lambda: f == val, "!=": lambda: f != val, "<": lambda: f < val, "<=": lambda: f <= val,
                 ">": lambda: f > val, ">=": lambda: f >= val, "in": lambda: f.isin(list(val)),
                 "not in": lambda: ~f.isin(list(val)), "is null": lambda: f.is_null(),
                 "is not null": lambda: f.is_valid()}[op]()
            e = c if e is None else (e & c)
        disj = e if disj is None else (disj | e)
    return disj


@pytest.fixture(scope="module")
def sorted_file(tmp_path_factory, table):
    path = tmp_path_factory.mktemp("f") / "sorted.parquet"
    carquet.write_table(table, path, row_group_size=100, write_page_index=True, max_rows_per_page=25)
    return path


FILTERS = [
    [("id", ">", 950)],
    [("id", ">=", 100), ("id", "<", 120)],
    [("id", "==", 5)],
    [("id", "=", 5)],
    [("id", "!=", 5)],
    [("id", "in", [1, 500, 999, 5000])],
    [("id", "not in", list(range(10, 1000)))],
    [("name", "is null")],
    [("name", "is not null"), ("id", "<", 20)],
    [("name", "==", None)],
    [("value", "<=", 10.0)],
    [("name", "==", "name-42")],
    [("name", "in", ["name-1", "name-2"])],
    [("ts", ">", dt.datetime(2024, 1, 1, 0, 15))],
    [[("id", "<", 3)], [("id", ">", 997)]],
    [("point.a", "==", 7)],
    [("id", ">", 10_000)],
]


@pytest.mark.parametrize("filters", FILTERS, ids=[str(f) for f in FILTERS])
def test_filters_match_pyarrow_compute(sorted_file, table, filters):
    expected = table.filter(to_expression(filters, table.column_names))
    out = carquet.read_table(sorted_file, filters=filters)
    assert out.equals(expected)
    # projection that does not include the filter columns
    out = carquet.read_table(sorted_file, columns=["value"], filters=filters)
    assert out.column_names == ["value"]
    assert out.column("value").equals(expected.column("value"))
    # streaming
    ids = [i for b in carquet.iter_batches(sorted_file, 64, columns=["id"], filters=filters)
           for i in b.column("id").to_pylist()]
    assert ids == expected.column("id").to_pylist()


def test_row_group_pruning(sorted_file):
    with carquet.ParquetFile(sorted_file) as pf:
        assert pf.filter_row_groups([("id", ">=", 950)]) == [9]
        assert pf.filter_row_groups([("id", "in", [5, 905])]) == [0, 9]
        assert pf.filter_row_groups([[("id", "<", 50)], [("id", ">", 850)]]) == [0, 8, 9]
        assert pf.filter_row_groups([("id", ">", 10_000)]) == []
        assert pf.filter_row_groups([("id", "in", [])]) == []
        assert pf.filter_row_groups([("name", "==", "zzz")]) == []
        # not pushable -> conservative
        assert len(pf.filter_row_groups([("id", "not in", [1])])) == 10
        assert len(pf.filter_row_groups([("tags", "==", 1)])) == 10


def test_typed_filter_values(tmp_path):
    t = pa.table({
        "d": pa.array([dt.date(2020, 1, 1) + dt.timedelta(days=i) for i in range(100)]),
        "ts": pa.array([dt.datetime(2020, 1, 1, tzinfo=dt.timezone.utc) + dt.timedelta(hours=i) for i in range(100)],
                       pa.timestamp("ms", tz="UTC")),
        "dec": pa.array([decimal.Decimal(i).scaleb(-1) for i in range(100)], pa.decimal128(5, 1)),
        "u": pa.array(range(100), pa.uint8()),
    })
    path = tmp_path / "typed.parquet"
    carquet.write_table(t, path, row_group_size=10)
    with carquet.ParquetFile(path) as pf:
        assert pf.filter_row_groups([("d", ">=", dt.date(2020, 3, 31))]) == [9]
        assert pf.filter_row_groups([("ts", "<", dt.datetime(2020, 1, 1, 5, tzinfo=dt.timezone.utc))]) == [0]
        assert pf.filter_row_groups([("dec", ">", decimal.Decimal("9.5"))]) == [9]
        assert pf.filter_row_groups([("u", "==", 42)]) == [4]
        out = pf.read(filters=[("dec", "<", decimal.Decimal("0.3"))])
        assert out.column("dec").to_pylist() == [decimal.Decimal("0.0"), decimal.Decimal("0.1"), decimal.Decimal("0.2")]
        # a value of the wrong type is an error
        with pytest.raises(TypeError):
            pf.read(filters=[("d", "==", "not a date")])


def test_expression_filter_rejected(sorted_file):
    """pyarrow.compute expressions need PyArrow: carquet takes the DNF form only."""
    with pytest.raises(TypeError):
        carquet.read_table(sorted_file, filters=pc.field("id") < 205)


def test_invalid_filters(sorted_file):
    with pytest.raises(ValueError):
        carquet.read_table(sorted_file, filters=[("id", "~", 1)])
    with pytest.raises(TypeError):
        carquet.read_table(sorted_file, filters=[("id", "in", 3)])
    with pytest.raises(TypeError):
        carquet.read_table(sorted_file, filters="id > 3")
    with pytest.raises(ValueError):
        carquet.read_table(sorted_file, filters=[[]])
    with pytest.raises(ValueError):
        carquet.read_table(sorted_file, filters=[("id", ">")])
    assert carquet.read_table(sorted_file, filters=[]).num_rows == 1000
    assert carquet.read_table(sorted_file, filters=("id", "<", 2)).num_rows == 2  # a bare clause


def test_page_filter_without_page_index(tmp_path, table):
    path = tmp_path / "nopi.parquet"
    carquet.write_table(table, path, row_group_size=500)
    assert carquet.read_table(path, filters=[("id", "<", 7)]).num_rows == 7
    assert carquet.read_table(path, filters=[("id", "in", [3, 600])]).column("id").to_pylist() == [3, 600]


def test_pyarrow_file_filters(tmp_path, table):
    path = tmp_path / "pa.parquet"
    pq.write_table(table, path, row_group_size=128, write_page_index=True, data_page_size=512)
    expected_table = pq.read_table(path)
    for f in FILTERS:
        expr = to_expression(f, table.column_names)
        assert carquet.read_table(path, filters=f).equals(expected_table.filter(expr)), f


# ---- pushdown must never drop a matching row (regressions) ----------------
def _check(tmp_path, table, filters, **write_opts):
    path = tmp_path / "p.parquet"
    carquet.write_table(table, path, write_page_index=True, max_rows_per_page=2, row_group_size=4, **write_opts)
    expected = table.filter(to_expression(filters, table.column_names))
    got = carquet.read_table(path, filters=filters)
    assert repr(got.to_pylist()) == repr(expected.to_pylist()), filters  # repr: NaN-aware


def test_in_on_fixed_len_columns(tmp_path):
    import uuid

    u = [uuid.UUID(int=i) for i in range(6)]
    t = pa.table({
        "fx": pa.array([i.to_bytes(4, "big") for i in range(6)], pa.binary(4)),
        "dec": pa.array([decimal.Decimal(i).scaleb(-2) for i in range(6)], pa.decimal128(5, 2)),
        "u": pa.array([x.bytes for x in u], pa.binary(16)),
    })
    _check(tmp_path, t, [("fx", "in", [b"\0\0\0\1", b"\0\0\0\5"])])
    _check(tmp_path, t, [("dec", "in", [decimal.Decimal("0.03")])])
    _check(tmp_path, t, [("u", "in", [u[4].bytes])])


def test_inexact_values_are_not_pruned(tmp_path):
    t = pa.table({
        "f": pa.array([0.1, 0.1, 0.5, 0.5, 0.7, 0.1], pa.float32()),
        "ts": pa.array([1, 1, 2, 2, 3, 3], pa.timestamp("ms")),
        "d": pa.array([dt.date(2024, 1, 1)] * 2 + [dt.date(2024, 1, 2)] * 4),
    })
    _check(tmp_path, t, [("f", ">", 0.1)])
    _check(tmp_path, t, [("f", "<=", 0.1)])
    _check(tmp_path, t, [("ts", "<", dt.datetime(1970, 1, 1, 0, 0, 0, 1500))])
    _check(tmp_path, t, [("ts", ">=", dt.datetime(1970, 1, 1, 0, 0, 0, 2500))])
    _check(tmp_path, t, [("d", "<", dt.datetime(2024, 1, 1, 12))])


def test_nan_rows_kept(tmp_path):
    nan = float("nan")
    t = pa.table({"x": pa.array([1.0, nan, 3.0, 4.0, nan, nan, 1.0, 2.0])})
    _check(tmp_path, t, [("x", "!=", 1.0)])
    _check(tmp_path, t, [("x", "in", [nan, 3.0])])
    _check(tmp_path, t, [("x", ">", 2.0)])


def test_out_of_range_values(tmp_path):
    t = pa.table({"i": pa.array(range(8), pa.int32()), "u": pa.array(range(8), pa.uint32())})
    _check(tmp_path, t, [("i", "<", 2**40)])
    _check(tmp_path, t, [("u", ">", -1)])
    path = tmp_path / "bf.parquet"
    carquet.write_table(t, path, bloom_filters=["u"])
    with carquet.ParquetFile(path) as pf:
        bf = pf.bloom_filter("u")
        assert -1 not in bf and 2**40 not in bf and 3 in bf
