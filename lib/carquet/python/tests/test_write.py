"""Writing: write_table options, ParquetWriter, append, sinks and inputs."""

from __future__ import annotations

import io

import pytest

import carquet
from carquet import CarquetError, StatusCode

pa = pytest.importorskip("pyarrow")
pq = pytest.importorskip("pyarrow.parquet")

CODECS = ["none", "snappy", "gzip", "zstd", "lz4", "lz4_hadoop", "brotli"]


@pytest.mark.parametrize("compression", CODECS)
def test_codecs(tmp_path, table, compression):
    path = tmp_path / f"c_{compression}.parquet"
    try:
        carquet.write_table(table, path, compression=compression)
    except CarquetError as exc:
        if compression == "brotli" and exc.status == StatusCode.UNSUPPORTED_CODEC:
            pytest.skip("library built without brotli")
        raise
    assert carquet.read_table(path).equals(table)
    if compression != "lz4_hadoop":
        assert pq.read_table(path).equals(table)
    expected = {"none": "UNCOMPRESSED", "lz4": "LZ4_RAW", "lz4_hadoop": "LZ4"}.get(compression, compression.upper())
    assert carquet.read_metadata(path).row_group(0).column(0).compression == expected


def test_compression_options(tmp_path, table):
    path = tmp_path / "lvl.parquet"
    carquet.write_table(table, path, compression="zstd", compression_level=10)
    assert carquet.read_table(path).equals(table)
    carquet.write_table(table, path, compression={"id": "snappy", "name": "gzip"},
                        compression_level={"name": 9})
    rg = carquet.read_metadata(path).row_group(0)
    assert rg.column("id").compression == "SNAPPY"
    assert rg.column("name").compression == "GZIP"
    assert rg.column("value").compression == "UNCOMPRESSED"
    assert carquet.read_table(path).equals(table)
    with pytest.raises(ValueError):
        carquet.write_table(table, tmp_path / "z.parquet", compression="lzma")
    with pytest.raises(KeyError):
        carquet.write_table(table, tmp_path / "z.parquet", compression={"nope": "zstd"})
    assert not (tmp_path / "z.parquet").exists()


def test_encodings(tmp_path, table):
    path = tmp_path / "enc.parquet"
    # low-cardinality strings, so the writer keeps the dictionary (it falls back
    # to PLAIN when a dictionary does not pay off)
    table = table.set_column(2, "name", pa.array([None if i % 7 == 3 else f"n{i % 20}" for i in range(1000)]))
    carquet.write_table(
        table, path,
        use_dictionary=["name"],
        column_encoding={"id": "delta_binary_packed", "tags": "plain"},
        use_byte_stream_split=["value"],
    )
    rg = carquet.read_metadata(path).row_group(0)
    assert "RLE_DICTIONARY" in rg.column("name").encodings
    assert rg.column("name").has_dictionary_page
    assert "DELTA_BINARY_PACKED" in rg.column("id").encodings
    assert "BYTE_STREAM_SPLIT" in rg.column("value").encodings
    assert carquet.read_table(path).equals(table)
    assert pq.read_table(path).equals(table)

    carquet.write_table(table, path, use_dictionary=True)
    rg = carquet.read_metadata(path).row_group(0)
    assert "RLE_DICTIONARY" in rg.column("name").encodings and "RLE_DICTIONARY" in rg.column("point.a").encodings
    carquet.write_table(table.select(["id"]), path, column_encoding="delta")
    assert "DELTA_BINARY_PACKED" in carquet.read_metadata(path).row_group(0).column(0).encodings
    with pytest.raises(ValueError):
        carquet.write_table(table, path, column_encoding="bogus")


def test_statistics_options(tmp_path, table):
    path = tmp_path / "st.parquet"
    carquet.write_table(table, path, write_statistics=["id"])
    rg = carquet.read_metadata(path).row_group(0)
    assert rg.column("id").statistics is not None and rg.column("id").statistics.has_min_max
    assert rg.column("value").statistics is None or not rg.column("value").statistics.has_min_max
    carquet.write_table(table, path, write_statistics=False)
    assert all(c.statistics is None or not c.statistics.has_min_max for c in carquet.read_metadata(path).row_group(0))
    long = pa.table({"s": ["x" * 100, "y" * 100]})
    carquet.write_table(long, path, max_statistics_size=8)
    st = carquet.read_metadata(path).row_group(0).column("s").statistics
    assert len(st.min_raw) <= 8 and len(st.max_raw) <= 8


def test_page_options(tmp_path, table):
    path = tmp_path / "pg.parquet"
    carquet.write_table(table, path, write_page_index=True, max_rows_per_page=100, data_page_version="2.0")
    with carquet.ParquetFile(path) as pf:
        oi = pf.offset_index("id")
        assert [p.first_row_index for p in oi] == list(range(0, 1000, 100))
        assert pf.read().equals(table)
    carquet.write_table(table, path, data_page_size=1024, write_page_index=True, write_batch_size=64)
    with carquet.ParquetFile(path) as pf:
        assert len(pf.offset_index("id")) >= 7
        assert pf.read().equals(table)
    carquet.write_table(table, path, data_page_size={"name": 512}, version="1.0", write_page_checksum=False)
    assert carquet.info(path).version == 1
    assert pq.read_table(path).equals(table)


def test_bloom_filters(tmp_path, table):
    path = tmp_path / "bf.parquet"
    carquet.write_table(table, path, bloom_filters={"id": {"ndv": 1000, "fpp": 0.01}, "name": True})
    with carquet.ParquetFile(path) as pf:
        rg = pf.metadata.row_group(0)
        assert rg.column("id").has_bloom_filter and rg.column("name").has_bloom_filter
        assert not rg.column("value").has_bloom_filter
        bf = pf.bloom_filter("id")
        assert all(i in bf for i in range(0, 1000, 37))
        assert sum(i in bf for i in range(10_000, 11_000)) < 100
        bn = pf.bloom_filter("name")
        assert "name-1" in bn
        assert bn.size > 0
        assert pf.bloom_filter("value") is None
        bf.close()
        with pytest.raises(ValueError):
            bf.might_contain(1)
    carquet.write_table(table, path, bloom_filters=["value"])
    assert carquet.read_metadata(path).row_group(0).column("value").has_bloom_filter


def test_sorting_columns_and_metadata(tmp_path, table):
    path = tmp_path / "sort.parquet"
    carquet.write_table(table, path, sorting_columns=["id", ("value", "descending")],
                        metadata={"origin": "test", "n": 3})
    md = pq.ParquetFile(path).metadata
    sc = md.row_group(0).sorting_columns
    assert [(c.column_index, c.descending) for c in sc] == [(0, False), (1, True)]
    kv = carquet.read_metadata(path).key_value_metadata
    assert kv["origin"] == "test" and kv["n"] == "3"
    assert "ARROW:schema" in kv
    if hasattr(pq, "SortingColumn"):
        carquet.write_table(table, path, sorting_columns=[pq.SortingColumn(2, descending=True, nulls_first=True)])
        c = pq.ParquetFile(path).metadata.row_group(0).sorting_columns[0]
        assert (c.column_index, c.descending, c.nulls_first) == (2, True, True)


def test_coerce_timestamps(tmp_path):
    import datetime as dt

    t = pa.table({"ts": pa.array([dt.datetime(2024, 1, 1, 0, 0, 0, 1500)], pa.timestamp("us"))})
    path = tmp_path / "ts.parquet"
    with pytest.raises(CarquetError):
        carquet.write_table(t, path, coerce_timestamps="ms")
    carquet.write_table(t, path, coerce_timestamps="ms", allow_truncated_timestamps=True, store_schema=False)
    out = carquet.read_table(path)
    assert out.schema.field("ts").type == pa.timestamp("ms")
    assert out.column("ts")[0] == dt.datetime(2024, 1, 1, 0, 0, 0, 1000)


def test_parquet_writer_packs_row_groups(tmp_path, table):
    path = tmp_path / "w.parquet"
    with carquet.ParquetWriter(path, table.schema, row_group_size=300) as w:
        for b in table.to_batches(max_chunksize=70):
            w.write(b)
        assert w.num_rows == 1000
    md = carquet.read_metadata(path)
    assert [rg.num_rows for rg in md.row_groups] == [300, 300, 300, 100]
    assert carquet.read_table(path).equals(table)


def test_parquet_writer_flush_and_inputs(tmp_path, table):
    pd = pytest.importorskip("pandas")
    path = tmp_path / "w2.parquet"
    small = table.select(["id", "name"]).slice(0, 10)
    with carquet.ParquetWriter(path, small.schema) as w:
        w.write_table(small)
        w.flush()
        w.write_batch(small.to_batches()[0])
        w.write({"id": [100, 101], "name": ["a", None]})
        w.write(pa.RecordBatchReader.from_batches(small.schema, small.to_batches()))
        w.write(small.to_pandas())
        w.add_metadata("k", "v")
        w.flush()
        w.flush()  # empty row group: no-op
    md = carquet.read_metadata(path)
    assert md.num_rows == 42
    assert md.row_groups[0].num_rows == 10
    assert md.key_value_metadata["k"] == "v"


def test_writer_casts_compatible_data(tmp_path):
    schema = pa.schema([("a", pa.int64()), ("b", pa.string())])
    path = tmp_path / "cast.parquet"
    with carquet.ParquetWriter(path, schema) as w:
        w.write(pa.table({"a": pa.array([1, 2], pa.int32()), "b": ["x", "y"]}))
        with pytest.raises(ValueError):
            w.write(pa.table({"zz": [1]}))
    assert carquet.read_table(path).schema.field("a").type == pa.int64()


def test_writer_errors_remove_partial_file(tmp_path, table):
    path = tmp_path / "bad.parquet"
    with pytest.raises(RuntimeError):
        with carquet.ParquetWriter(path, table.schema) as w:
            w.write(table)
            raise RuntimeError("boom")
    assert not path.exists()
    w = carquet.ParquetWriter(path, table.schema)
    w.abort()
    assert not path.exists() and w.closed
    with pytest.raises(CarquetError):
        w.write(table)
    with pytest.raises(TypeError):
        carquet.ParquetWriter(path)
    with pytest.raises(TypeError):
        carquet.ParquetWriter(12345, table.schema)


def test_append(tmp_path, table):
    path = tmp_path / "app.parquet"
    carquet.write_table(table.slice(0, 600), path, row_group_size=300)
    carquet.write_table(table.slice(600), path, append=True)
    with carquet.ParquetWriter(path, append=True) as w:  # schema taken from the file
        w.write(table.slice(0, 10))
    md = carquet.read_metadata(path)
    assert [rg.num_rows for rg in md.row_groups] == [300, 300, 400, 10]
    out = carquet.read_table(path)
    assert out.slice(0, 1000).equals(table)
    assert pq.read_table(path).num_rows == 1010
    assert list(md.key_value_metadata).count("ARROW:schema") == 1
    with pytest.raises(ValueError):
        carquet.ParquetWriter(io.BytesIO(), table.schema, append=True)


def test_file_like_sinks(tmp_path, table):
    buf = io.BytesIO()
    carquet.write_table(table, buf, compression="snappy")
    data = buf.getvalue()
    assert data[:4] == b"PAR1" and data[-4:] == b"PAR1"
    assert carquet.read_table(data).equals(table)
    assert pq.read_table(io.BytesIO(data)).equals(table)
    with open(tmp_path / "f.parquet", "wb") as f:
        with carquet.ParquetWriter(f, table.schema) as w:
            w.write(table)
    assert carquet.read_table(tmp_path / "f.parquet").equals(table)


def test_inputs(tmp_path, table):
    path = tmp_path / "in.parquet"
    carquet.write_table(table.to_batches()[0], path)
    assert carquet.read_table(path).equals(table)
    carquet.write_table(table.to_batches(max_chunksize=100), path)
    assert carquet.read_table(path).equals(table)
    carquet.write_table({"x": [1, 2, 3]}, path)
    assert carquet.read_table(path).column("x").to_pylist() == [1, 2, 3]
    carquet.write_table(pa.RecordBatchReader.from_batches(table.schema, table.to_batches()), path)
    assert carquet.read_table(path).equals(table)
    with pytest.raises(TypeError):
        carquet.write_table(42, path)
    with pytest.raises(ValueError):
        carquet.write_table([], path)


def test_polars_input(tmp_path):
    pl = pytest.importorskip("polars")
    df = pl.DataFrame({"a": [1, 2, 3], "b": ["x", "y", None]})
    path = tmp_path / "pl.parquet"
    carquet.write_table(df, path)
    assert pl.DataFrame(carquet.read_table(path)).equals(df)
    assert carquet.read_table(path).to_polars().equals(df)


def test_empty_and_wide(tmp_path):
    path = tmp_path / "wide.parquet"
    wide = pa.table({f"c{i}": pa.array([i, i + 1], pa.int32()) for i in range(300)})
    carquet.write_table(wide, path)
    assert carquet.read_table(path).equals(wide)
    assert carquet.read_table(path, columns=["c299", "c7"]).column_names == ["c299", "c7"]


def test_row_group_size_validation(tmp_path, table):
    with pytest.raises(ValueError):
        carquet.write_table(table, tmp_path / "x.parquet", row_group_size=-1)
    with pytest.raises(ValueError):
        carquet.write_table(table, tmp_path / "x.parquet", data_page_version=3)
    with pytest.raises(ValueError):
        carquet.write_table(table, tmp_path / "x.parquet", coerce_timestamps="hours")


def test_storage_types_roundtrip(tmp_path):
    """Types without a Parquet mapping in the C bridge go through their storage type and come back."""
    cols = {
        "dur": pa.array([1, 2, None], pa.duration("ms")),
        "null": pa.array([None] * 3, pa.null()),
        "fsl": pa.array([[1, 2], [3, 4], None], pa.list_(pa.int32(), 2)),
        "vec": pa.array([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]], pa.list_(pa.float32(), 2)),
        "st": pa.array([{"d": 5, "x": [1]}] * 3, pa.struct([("d", pa.duration("s")), ("x", pa.list_(pa.int8()))])),
    }
    if hasattr(pa, "string_view"):
        cols["sv"] = pa.array(["a", "b", None], pa.string_view())
    t = pa.table(cols)
    path = tmp_path / "storage.parquet"
    carquet.write_table(t, path, row_group_size=2)
    assert carquet.read_table(path).equals(t)
    assert carquet.read_parquet_schema(path)["dur"].physical_type == "INT64"


def test_nullability_mismatch_is_conformed(tmp_path):
    path = tmp_path / "nn.parquet"
    schema = pa.schema([pa.field("x", pa.int64(), nullable=False)])
    with carquet.ParquetWriter(path, schema) as w:
        w.write(pa.table({"x": [1, 2, 3]}))  # nullable field, no nulls
    assert carquet.read_parquet_schema(path)["x"].repetition == "REQUIRED"
    carquet.write_table(pa.table({"x": [4]}), path, append=True)  # append to a REQUIRED column
    assert carquet.read_table(path).column("x").to_pylist() == [1, 2, 3, 4]
    with carquet.ParquetWriter(tmp_path / "n2.parquet", pa.schema([("x", pa.int64())])) as w:
        w.write(pa.table({"x": pa.array([1])}, schema=schema))
    # a null where the file's column is REQUIRED is refused, not written as garbage
    with pytest.raises(CarquetError):
        carquet.write_table(pa.table({"x": [None]}), path, append=True)


def test_unclosed_writer_is_finished_on_gc(tmp_path, table):
    import gc

    path = tmp_path / "gc.parquet"
    w = carquet.ParquetWriter(path, table.schema)
    w.write(table)
    del w
    gc.collect()
    assert carquet.read_table(path).equals(table)
