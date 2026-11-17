#!/usr/bin/env python3
"""
Generate Parquet test files using various libraries for interoperability testing.

Requirements:
    pip install pyarrow pandas duckdb fastparquet

Usage:
    python generate_test_files.py [output_dir]
"""

import os
import sys
import json
from pathlib import Path
from decimal import Decimal

# Try to import each library, track availability
AVAILABLE_LIBS = {}

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    AVAILABLE_LIBS['pyarrow'] = pa.__version__
except ImportError:
    print("Warning: pyarrow not available")

try:
    import pandas as pd
    AVAILABLE_LIBS['pandas'] = pd.__version__
except ImportError:
    print("Warning: pandas not available")

try:
    import duckdb
    AVAILABLE_LIBS['duckdb'] = duckdb.__version__
except ImportError:
    print("Warning: duckdb not available")

try:
    import fastparquet
    AVAILABLE_LIBS['fastparquet'] = fastparquet.__version__
except ImportError:
    print("Warning: fastparquet not available")


def generate_pyarrow_files(output_dir: Path):
    """Generate test files using PyArrow."""
    if 'pyarrow' not in AVAILABLE_LIBS:
        return []

    files = []
    pa_dir = output_dir / "pyarrow"
    pa_dir.mkdir(exist_ok=True)

    # 1. Basic types - uncompressed
    print("  Generating: basic_types_uncompressed.parquet")
    table = pa.table({
        'int32_col': pa.array([1, 2, 3, None, 5], type=pa.int32()),
        'int64_col': pa.array([100, 200, 300, 400, None], type=pa.int64()),
        'float_col': pa.array([1.1, 2.2, None, 4.4, 5.5], type=pa.float32()),
        'double_col': pa.array([1.11, 2.22, 3.33, 4.44, 5.55], type=pa.float64()),
        'bool_col': pa.array([True, False, True, None, False], type=pa.bool_()),
        'string_col': pa.array(['hello', 'world', None, 'test', 'data'], type=pa.string()),
    })
    path = pa_dir / "basic_types_uncompressed.parquet"
    pq.write_table(table, path, compression='NONE')
    files.append(('basic_types_uncompressed', str(path), 'pyarrow', 'NONE'))

    # 2. Basic types - SNAPPY
    print("  Generating: basic_types_snappy.parquet")
    path = pa_dir / "basic_types_snappy.parquet"
    pq.write_table(table, path, compression='SNAPPY')
    files.append(('basic_types_snappy', str(path), 'pyarrow', 'SNAPPY'))

    # 3. Basic types - GZIP
    print("  Generating: basic_types_gzip.parquet")
    path = pa_dir / "basic_types_gzip.parquet"
    pq.write_table(table, path, compression='GZIP')
    files.append(('basic_types_gzip', str(path), 'pyarrow', 'GZIP'))

    # 4. Basic types - ZSTD
    print("  Generating: basic_types_zstd.parquet")
    path = pa_dir / "basic_types_zstd.parquet"
    pq.write_table(table, path, compression='ZSTD')
    files.append(('basic_types_zstd', str(path), 'pyarrow', 'ZSTD'))

    # 5. Basic types - LZ4
    print("  Generating: basic_types_lz4.parquet")
    path = pa_dir / "basic_types_lz4.parquet"
    pq.write_table(table, path, compression='LZ4')
    files.append(('basic_types_lz4', str(path), 'pyarrow', 'LZ4'))

    # 6. Large dataset (1M rows)
    print("  Generating: large_dataset.parquet")
    import numpy as np
    n = 1_000_000
    large_table = pa.table({
        'id': pa.array(range(n), type=pa.int64()),
        'value': pa.array(np.random.randn(n), type=pa.float64()),
        'category': pa.array([f'cat_{i % 100}' for i in range(n)], type=pa.string()),
    })
    path = pa_dir / "large_dataset.parquet"
    pq.write_table(large_table, path, compression='ZSTD', row_group_size=100000)
    files.append(('large_dataset', str(path), 'pyarrow', 'ZSTD'))

    # 7. Dictionary encoding
    print("  Generating: dictionary_encoded.parquet")
    dict_table = pa.table({
        'category': pa.array(['A', 'B', 'A', 'C', 'B', 'A'] * 1000, type=pa.string()),
        'value': pa.array(range(6000), type=pa.int32()),
    })
    path = pa_dir / "dictionary_encoded.parquet"
    pq.write_table(dict_table, path, use_dictionary=True)
    files.append(('dictionary_encoded', str(path), 'pyarrow', 'NONE'))

    # 8. Multiple row groups
    print("  Generating: multiple_row_groups.parquet")
    path = pa_dir / "multiple_row_groups.parquet"
    pq.write_table(large_table[:100000], path, row_group_size=10000)
    files.append(('multiple_row_groups', str(path), 'pyarrow', 'NONE'))

    # 9. All null column
    print("  Generating: all_nulls.parquet")
    null_table = pa.table({
        'all_null_int': pa.array([None] * 100, type=pa.int32()),
        'all_null_str': pa.array([None] * 100, type=pa.string()),
        'has_values': pa.array(range(100), type=pa.int32()),
    })
    path = pa_dir / "all_nulls.parquet"
    pq.write_table(null_table, path)
    files.append(('all_nulls', str(path), 'pyarrow', 'NONE'))

    # 10. Empty table
    print("  Generating: empty_table.parquet")
    empty_table = pa.table({
        'col1': pa.array([], type=pa.int32()),
        'col2': pa.array([], type=pa.string()),
    })
    path = pa_dir / "empty_table.parquet"
    pq.write_table(empty_table, path)
    files.append(('empty_table', str(path), 'pyarrow', 'NONE'))

    # 11. Binary data
    print("  Generating: binary_data.parquet")
    binary_table = pa.table({
        'binary_col': pa.array([b'\x00\x01\x02', b'\xff\xfe\xfd', None, b'hello'], type=pa.binary()),
        'fixed_binary': pa.array([b'AAAA', b'BBBB', b'CCCC', b'DDDD'], type=pa.binary(4)),
    })
    path = pa_dir / "binary_data.parquet"
    pq.write_table(binary_table, path)
    files.append(('binary_data', str(path), 'pyarrow', 'NONE'))

    # 12. Date and time types
    print("  Generating: datetime_types.parquet")
    from datetime import date, datetime
    dt_table = pa.table({
        'date_col': pa.array([date(2024, 1, 1), date(2024, 6, 15), None], type=pa.date32()),
        'timestamp_col': pa.array([
            datetime(2024, 1, 1, 12, 0, 0),
            datetime(2024, 6, 15, 18, 30, 0),
            None
        ], type=pa.timestamp('us')),
    })
    path = pa_dir / "datetime_types.parquet"
    pq.write_table(dt_table, path)
    files.append(('datetime_types', str(path), 'pyarrow', 'NONE'))

    # 13. Decimal types
    print("  Generating: decimal_types.parquet")
    dec_table = pa.table({
        'decimal_col': pa.array([Decimal('123.45'), Decimal('678.90'), Decimal('-111.11')], type=pa.decimal128(10, 2)),
    })
    path = pa_dir / "decimal_types.parquet"
    pq.write_table(dec_table, path)
    files.append(('decimal_types', str(path), 'pyarrow', 'NONE'))

    # 14. Delta encoding (integers with small deltas)
    print("  Generating: delta_encoded.parquet")
    delta_table = pa.table({
        'sequential': pa.array(range(10000), type=pa.int64()),
        'timestamps': pa.array([1704067200 + i for i in range(10000)], type=pa.int64()),
    })
    path = pa_dir / "delta_encoded.parquet"
    pq.write_table(delta_table, path, use_dictionary=False)
    files.append(('delta_encoded', str(path), 'pyarrow', 'NONE'))

    files.extend(generate_pyarrow_feature_files(pa_dir))
    return files


def _pa_write_optional(name, pa_dir, table, files, label, **kwargs):
    """Write one feature file; skip (with a note) if this PyArrow version
    lacks a writer option rather than failing the whole generation."""
    path = pa_dir / f"{name}.parquet"
    try:
        pq.write_table(table, path, **kwargs)
    except (TypeError, ValueError, NotImplementedError, pa.ArrowException) as e:
        print(f"  Skipping: {name}.parquet ({e})")
        return
    print(f"  Generating: {name}.parquet")
    files.append((name, str(path), 'pyarrow', label))


def generate_pyarrow_feature_files(pa_dir: Path):
    """Encodings / page formats / codecs carquet must decode from PyArrow.

    Values are deterministic and every file spans several pages (small
    data_page_size) and row groups, so run_interop.py's full-value oracle
    exercises page and row-group boundaries, not only the first page."""
    files = []
    n = 20000

    def nullable(values, every):
        return [None if i % every == 0 else v for i, v in enumerate(values)]

    ints32 = nullable([(i * 7919) % 100003 - 50000 for i in range(n)], 11)
    ints64 = nullable([i * 1000003 - 7 * 10**9 for i in range(n)], 13)
    floats = nullable([i * 0.25 - 999.5 for i in range(n)], 7)
    doubles = nullable([i * 1e-3 - 3.5 for i in range(n)], 9)
    words = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta']
    strs = nullable([f"{words[i % 7]}_{i // 3}" for i in range(n)], 5)
    bools = nullable([i % 3 == 0 for i in range(n)], 17)

    flat = pa.table({
        'i32': pa.array(ints32, type=pa.int32()),
        'i64': pa.array(ints64, type=pa.int64()),
        'f32': pa.array(floats, type=pa.float32()),
        'f64': pa.array(doubles, type=pa.float64()),
        'str': pa.array(strs, type=pa.string()),
        'flag': pa.array(bools, type=pa.bool_()),
    })
    common = dict(row_group_size=7000, data_page_size=2048)

    # DELTA_BINARY_PACKED / DELTA_LENGTH_BYTE_ARRAY / DELTA_BYTE_ARRAY
    delta = flat.append_column('str_dba', pa.array(strs, type=pa.string()))
    _pa_write_optional('enc_delta_all', pa_dir, delta, files, 'ZSTD',
                       compression='ZSTD', use_dictionary=False,
                       column_encoding={'i32': 'DELTA_BINARY_PACKED',
                                        'i64': 'DELTA_BINARY_PACKED',
                                        'str': 'DELTA_LENGTH_BYTE_ARRAY',
                                        'str_dba': 'DELTA_BYTE_ARRAY'},
                       **common)
    # BYTE_STREAM_SPLIT on FLOAT / DOUBLE / INT32 / INT64
    _pa_write_optional('enc_bss_all', pa_dir, flat.select(['i32', 'i64', 'f32', 'f64']),
                       files, 'SNAPPY', compression='SNAPPY', use_dictionary=False,
                       column_encoding={'i32': 'BYTE_STREAM_SPLIT',
                                        'i64': 'BYTE_STREAM_SPLIT',
                                        'f32': 'BYTE_STREAM_SPLIT',
                                        'f64': 'BYTE_STREAM_SPLIT'},
                       **common)
    # Data Page V2 with dictionary + nulls, and V2 without dictionary
    _pa_write_optional('pagev2_dict', pa_dir, flat, files, 'ZSTD',
                       compression='ZSTD', data_page_version='2.0',
                       use_dictionary=True, **common)
    _pa_write_optional('pagev2_plain_gzip', pa_dir, flat, files, 'GZIP',
                       compression='GZIP', data_page_version='2.0',
                       use_dictionary=False, **common)
    # Page index + page CRCs (carquet verifies CRCs when asked)
    _pa_write_optional('page_index_crc_lz4raw', pa_dir, flat, files, 'LZ4_RAW',
                       compression='LZ4', write_page_index=True,
                       write_page_checksum=True, **common)
    # BROTLI
    _pa_write_optional('brotli', pa_dir, flat, files, 'BROTLI',
                       compression='BROTLI', **common)
    # Wider types: float16, large_string / large_binary, uint, fixed binary
    wide = pa.table({
        'f16': pa.array([float((i % 2048) - 1024) for i in range(n)], type=pa.float16()),
        'lstr': pa.array(strs, type=pa.large_string()),
        'lbin': pa.array([None if s is None else s.encode() for s in strs], type=pa.large_binary()),
        'u16': pa.array([(i * 13) % 65536 for i in range(n)], type=pa.uint16()),
        'u64': pa.array([(i * 11400714819323198485) % (1 << 64) for i in range(n)], type=pa.uint64()),
        'i8': pa.array([(i % 256) - 128 for i in range(n)], type=pa.int8()),
    })
    _pa_write_optional('wide_types', pa_dir, wide, files, 'ZSTD',
                       compression='ZSTD', **common)
    # Nested (decode-checked by --dir; not value-comparable by the flat oracle)
    nested = pa.table({
        'id': pa.array(range(2000), type=pa.int32()),
        'tags': pa.array([None if i % 9 == 0 else [i, i + 1][: i % 3]
                          for i in range(2000)], type=pa.list_(pa.int32())),
        'kv': pa.array([[('k%d' % (i % 5), i)] if i % 4 else []
                        for i in range(2000)], type=pa.map_(pa.string(), pa.int64())),
        'st': pa.array([None if i % 6 == 0 else {'a': i, 'b': str(i)}
                        for i in range(2000)],
                       type=pa.struct([('a', pa.int32()), ('b', pa.string())])),
    })
    _pa_write_optional('nested_list_map_struct', pa_dir, nested, files, 'SNAPPY',
                       compression='SNAPPY', row_group_size=700, data_page_size=1024)
    return files


def generate_duckdb_files(output_dir: Path):
    """Generate test files using DuckDB."""
    if 'duckdb' not in AVAILABLE_LIBS:
        return []

    files = []
    duck_dir = output_dir / "duckdb"
    duck_dir.mkdir(exist_ok=True)

    conn = duckdb.connect()

    # 1. Basic types
    print("  Generating: duckdb_basic.parquet")
    path = duck_dir / "duckdb_basic.parquet"
    conn.execute(f"""
        COPY (
            SELECT
                i::INT32 as int_col,
                i::INT64 as bigint_col,
                i::DOUBLE as double_col,
                'row_' || i as string_col
            FROM range(1000) t(i)
        ) TO '{path}' (FORMAT PARQUET)
    """)
    files.append(('duckdb_basic', str(path), 'duckdb', 'NONE'))

    # 2. With compression
    print("  Generating: duckdb_zstd.parquet")
    path = duck_dir / "duckdb_zstd.parquet"
    conn.execute(f"""
        COPY (
            SELECT
                i::INT64 as id,
                random() as value,
                'category_' || (i % 10) as category
            FROM range(100000) t(i)
        ) TO '{path}' (FORMAT PARQUET, COMPRESSION ZSTD)
    """)
    files.append(('duckdb_zstd', str(path), 'duckdb', 'ZSTD'))

    # 3. With nulls
    print("  Generating: duckdb_nulls.parquet")
    path = duck_dir / "duckdb_nulls.parquet"
    conn.execute(f"""
        COPY (
            SELECT
                CASE WHEN i % 3 = 0 THEN NULL ELSE i END as nullable_int,
                CASE WHEN i % 5 = 0 THEN NULL ELSE 'val_' || i END as nullable_str
            FROM range(1000) t(i)
        ) TO '{path}' (FORMAT PARQUET)
    """)
    files.append(('duckdb_nulls', str(path), 'duckdb', 'NONE'))

    conn.close()
    return files


def generate_fastparquet_files(output_dir: Path):
    """Generate test files using fastparquet."""
    if 'fastparquet' not in AVAILABLE_LIBS or 'pandas' not in AVAILABLE_LIBS:
        return []

    files = []
    fp_dir = output_dir / "fastparquet"
    fp_dir.mkdir(exist_ok=True)

    # 1. Basic types
    print("  Generating: fastparquet_basic.parquet")
    df = pd.DataFrame({
        'int_col': [1, 2, 3, 4, 5],
        'float_col': [1.1, 2.2, 3.3, 4.4, 5.5],
        'string_col': ['a', 'b', 'c', 'd', 'e'],
    })
    path = fp_dir / "fastparquet_basic.parquet"
    fastparquet.write(str(path), df)
    files.append(('fastparquet_basic', str(path), 'fastparquet', 'NONE'))

    # 2. With compression
    print("  Generating: fastparquet_snappy.parquet")
    path = fp_dir / "fastparquet_snappy.parquet"
    fastparquet.write(str(path), df, compression='SNAPPY')
    files.append(('fastparquet_snappy', str(path), 'fastparquet', 'SNAPPY'))

    return files


def main():
    output_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "test_files")
    output_dir.mkdir(exist_ok=True)

    print(f"Output directory: {output_dir}")
    print(f"Available libraries: {AVAILABLE_LIBS}")
    print()

    all_files = []

    # Generate files with each library
    print("Generating PyArrow files...")
    all_files.extend(generate_pyarrow_files(output_dir))

    print("\nGenerating DuckDB files...")
    all_files.extend(generate_duckdb_files(output_dir))

    print("\nGenerating fastparquet files...")
    all_files.extend(generate_fastparquet_files(output_dir))

    # Write manifest
    manifest = {
        'libraries': AVAILABLE_LIBS,
        'files': [
            {
                'name': name,
                'path': path,
                'producer': producer,
                'compression': compression
            }
            for name, path, producer, compression in all_files
        ]
    }

    manifest_path = output_dir / "manifest.json"
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    print(f"\nGenerated {len(all_files)} test files")
    print(f"Manifest written to: {manifest_path}")

    return 0


if __name__ == '__main__':
    sys.exit(main())
