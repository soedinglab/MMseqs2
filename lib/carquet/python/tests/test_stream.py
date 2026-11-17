"""ParquetFile.arrow_stream() / __arrow_c_stream__ (Arrow PyCapsule interface)."""

from __future__ import annotations

import datetime as dt
import decimal
import gc
import subprocess
import sys
import textwrap

import pytest

import carquet
from carquet import _ffi

from conftest import make_table

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")

pytestmark = pytest.mark.skipif(
    not _ffi.has_symbol("carquet_reader_export_arrow_stream"), reason="libcarquet without the Arrow stream"
)


def _typed_table(n: int = 300) -> "pa.Table":
    return pa.table(
        {
            "id": pa.array(range(n), pa.int64()),
            "d9": pa.array([None if i % 5 == 0 else decimal.Decimal(i * 37 - 4000).scaleb(-2) for i in range(n)],
                           pa.decimal128(9, 2)),
            "d38": pa.array([decimal.Decimal(-(10**30) + i).scaleb(-5) for i in range(n)], pa.decimal128(38, 5)),
            "d50": pa.array([decimal.Decimal(10**45 * (-1) ** i + i) for i in range(n)], pa.decimal256(50, 0)),
            "name": pa.array([None if i % 7 == 3 else f"n{i}" for i in range(n)], pa.string()),
        }
    )


def test_stream_matches_pyarrow_types_and_values(tmp_path):
    t = _typed_table()
    path = tmp_path / "typed.parquet"
    pq.write_table(t, path, row_group_size=100, store_schema=False)
    with carquet.ParquetFile(path) as pf:
        got = pa.table(pf)  # __arrow_c_stream__
    expected = pq.read_table(path)
    assert got.schema.types == expected.schema.types
    assert got.equals(expected)


def test_int96_and_integer_decimals(tmp_path):
    epoch = dt.datetime(1999, 12, 31, 23, 59, 59, 123456)
    t = pa.table({
        "ts": pa.array([epoch + dt.timedelta(days=i, microseconds=i) for i in range(50)], pa.timestamp("ns")),
        "d": pa.array([decimal.Decimal(i - 25).scaleb(-1) for i in range(50)], pa.decimal128(7, 1)),
    })
    path = tmp_path / "int96.parquet"
    try:
        pq.write_table(t, path, use_deprecated_int96_timestamps=True, store_schema=False,
                       store_decimal_as_integer=True)
    except TypeError:
        pq.write_table(t, path, use_deprecated_int96_timestamps=True, store_schema=False)
    with carquet.ParquetFile(path) as pf:
        assert pf.parquet_schema.column("ts").physical_type == "INT96"
        got = pf.arrow_stream().read_all()
    assert got.equals(pq.read_table(path))


def test_projection_batches_and_nested(tmp_path):
    t = make_table(1000)
    path = tmp_path / "t.parquet"
    pq.write_table(t, path, row_group_size=400)
    with carquet.ParquetFile(path) as pf:
        # flat: batch_size honoured, no batch spans row groups
        reader = pa.RecordBatchReader.from_stream(pf.arrow_stream(["name", "id"], batch_size=150))
        assert reader.schema.names == ["name", "id"]
        assert [b.num_rows for b in reader] == [150, 150, 100, 150, 150, 100, 150, 50]
        # a list: the batch reader reads repeated columns a row group at a time
        reader = pa.RecordBatchReader.from_stream(pf.arrow_stream(["tags", "id"], batch_size=150))
        assert reader.schema.names == ["tags", "id"]
        assert [b.num_rows for b in reader] == [400, 400, 200]
        # a struct: one batch per row group
        reader = pa.RecordBatchReader.from_stream(pf.arrow_stream(["point", "id"], batch_size=150))
        assert [b.num_rows for b in reader] == [400, 400, 200]
        full = pf.arrow_stream().read_all()
    assert full.select(["id", "value", "name", "tags"]).equals(t.select(["id", "value", "name", "tags"]))
    assert full.column("point").to_pylist() == t.column("point").to_pylist()


def test_stream_outlives_file_and_buffer(tmp_path):
    t = make_table(500)
    sink = pa.BufferOutputStream()
    pq.write_table(t, sink)
    data = bytearray(sink.getvalue().to_pybytes())
    pf = carquet.ParquetFile(data)
    capsule = pf.arrow_stream(["id", "name"]).__arrow_c_stream__()
    pf.close()
    del pf, data
    gc.collect()
    got = pa.RecordBatchReader._import_from_c_capsule(capsule).read_all()
    assert got.equals(t.select(["id", "name"]))

    path = tmp_path / "p.parquet"
    pq.write_table(t, path)
    with carquet.ParquetFile(path) as pf:
        stream = pf.arrow_stream(["id"])
    assert pa.table(stream).column("id").equals(t.column("id"))


def test_stream_reusable_and_errors(tmp_path):
    path = tmp_path / "t.parquet"
    t = make_table(10)
    pq.write_table(t, path)
    with carquet.ParquetFile(path) as pf:
        s = pf.arrow_stream(["id"])
        assert pa.table(s).equals(pa.table(s))  # each call reads afresh
        with pytest.raises(KeyError):
            pf.arrow_stream(["nope"])
        with pytest.raises(ValueError):
            pf.arrow_stream(batch_size=0)
        # a capsule nobody consumes is released with it
        cap = pf.arrow_stream().__arrow_c_stream__()
        del cap
        gc.collect()
    with pytest.raises(carquet.CarquetError):
        pf.arrow_stream()


def test_duckdb_reads_the_stream(tmp_path):
    duckdb = pytest.importorskip("duckdb")
    path = tmp_path / "t.parquet"
    t = _typed_table(1000)
    pq.write_table(t, path, row_group_size=300)
    with carquet.ParquetFile(path) as pf:
        stream = pf.arrow_stream(["id", "d9"])
        ((n, s_id, s_d9),) = duckdb.sql("SELECT count(*), sum(id), sum(d9) FROM stream").fetchall()
    assert n == 1000 and s_id == sum(range(1000))
    assert s_d9 == sum(v for v in t.column("d9").to_pylist() if v is not None)


def test_stream_without_pyarrow(tmp_path):
    """The capsule is produced and readable with pyarrow unimportable."""
    path = tmp_path / "t.parquet"
    pq.write_table(_typed_table(250), path, row_group_size=100)
    code = textwrap.dedent(
        f"""
        import ctypes, sys
        sys.modules["pyarrow"] = None  # any import of pyarrow now fails
        import carquet
        from carquet._arrow import ArrowArray, ArrowArrayStream, ArrowSchema
        pf = carquet.ParquetFile({str(path)!r})
        cap = pf.__arrow_c_stream__()
        get = ctypes.pythonapi.PyCapsule_GetPointer
        get.restype = ctypes.c_void_p
        get.argtypes = [ctypes.py_object, ctypes.c_char_p]
        s = ArrowArrayStream.from_address(get(cap, b"arrow_array_stream"))
        sc = ArrowSchema()
        assert s.get_schema(ctypes.byref(s), ctypes.byref(sc)) == 0
        fmts = [sc.children[i].contents.format.decode() for i in range(sc.n_children)]
        sc.release(ctypes.byref(sc))
        rows = 0
        while True:
            a = ArrowArray()
            assert s.get_next(ctypes.byref(s), ctypes.byref(a)) == 0
            if not a.release:
                break
            rows += a.length
            a.release(ctypes.byref(a))
        print(fmts, rows)
        assert "pyarrow" not in [m for m in sys.modules if sys.modules[m] is not None]
        """
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=False)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "['l', 'd:9,2', 'd:38,5', 'd:50,0,256', 'u'] 250"
