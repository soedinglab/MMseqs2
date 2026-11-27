#!/usr/bin/env python3
"""
Benchmark for fastparquet - compare with Carquet.

Same data, iteration counts, trimmed median and cache purge as
benchmark_carquet.c / benchmark_pyarrow.py. Settings are aligned where
fastparquet allows it (no dictionary, same row group size, zstd level, LZ4_RAW
for "lz4"); fastparquet cannot write BYTE_STREAM_SPLIT or page checksums and
reads through pandas, so its files and read path differ from the others.

"<size>_nullable" datasets use pandas' masked dtypes (Int64 / Float64) with the
shared null pattern; `category` is a pandas Categorical, which is the only way
fastparquet writes a dictionary-encoded column (its codes are computed when the
data is generated, outside the timed write). There is no mmap read option, so
no read_nommap figure is reported.
"""

import os
import time

import fastparquet
import pandas as pd

from bench_common import (
    DICTIONARY_COLUMNS,
    WARMUP_ITERATIONS,
    generate_columns,
    get_benchmark_zstd_level,
    get_iterations,
    get_temp_dir,
    is_nullable_dataset,
    null_masks,
    purge_file_cache,
    trimmed_median,
)

# Carquet's "lz4" is LZ4_RAW; fastparquet's legacy "LZ4" id is not.
CODECS = {"none": None, "snappy": "SNAPPY", "zstd": "ZSTD", "lz4": "LZ4_RAW"}


def generate_data(num_rows, nullable=False):
    """Pre-generate test data using the same RNG sequence as benchmark_carquet.c."""
    columns = generate_columns(num_rows)
    if not nullable:
        return pd.DataFrame(columns)
    masks = null_masks(num_rows)
    data = {}
    for name, values in columns.items():
        if values.dtype.kind == "f":
            arr = pd.arrays.FloatingArray(values, masks[name])
        else:
            arr = pd.arrays.IntegerArray(values, masks[name])
        data[name] = pd.Categorical(arr) if name in DICTIONARY_COLUMNS else arr
    return pd.DataFrame(data)


def expected_null_counts(num_rows):
    return {name: int(mask.sum()) for name, mask in null_masks(num_rows).items()}


def write_compression(compression_name, nullable=False):
    codec = CODECS[compression_name]
    if codec == "ZSTD":
        comp = {"_default": {"type": "ZSTD", "args": {"level": get_benchmark_zstd_level()}}}
        if nullable:
            # fastparquet calls .upper() on the codec when it writes a
            # Categorical's dictionary page, so those columns cannot take the
            # {"type", "args"} form and use cramjam's default zstd level.
            for name in DICTIONARY_COLUMNS:
                comp[name] = "ZSTD"
        return comp
    return codec


def benchmark_write(filename, df, compression, nullable=False):
    """Write a Parquet file and return time_ms."""
    start = time.perf_counter()
    num_rows = len(df)
    rg_size = num_rows // 10 if num_rows > 1_000_000 else 100_000
    fastparquet.write(
        filename,
        df,
        row_group_offsets=rg_size,
        compression=compression,
        write_index=False,
        has_nulls=nullable,
    )
    return (time.perf_counter() - start) * 1000


def benchmark_read(filename, expected_rows, expected_nulls=None):
    """Read a Parquet file and return time_ms."""
    start = time.perf_counter()
    df = fastparquet.ParquetFile(filename).to_pandas()

    nulls = None
    checksum = 0
    if expected_nulls is None:
        for col_name in df.columns:
            checksum += df[col_name].to_numpy()[0]
    else:
        # Masked / categorical columns: count nulls from the masks instead of
        # to_numpy(), which would box them into object arrays.
        nulls = {col_name: int(df[col_name].isna().sum()) for col_name in df.columns}

    elapsed_ms = (time.perf_counter() - start) * 1000

    if len(df) != expected_rows:
        raise ValueError(f"Row count mismatch: {len(df)} vs {expected_rows}")
    if expected_nulls is not None and nulls != expected_nulls:
        raise ValueError(f"Null count mismatch: {nulls} vs {expected_nulls}")

    return elapsed_ms


def run_benchmark(name, num_rows, compression_name):
    """Run a single benchmark configuration."""
    filename = os.path.join(
        get_temp_dir(),
        f"benchmark_{name}_{compression_name}_fastparquet.parquet",
    )
    nullable = is_nullable_dataset(name)
    compression = write_compression(compression_name, nullable)

    print(f"\n=== {name} ({num_rows:,} rows, {compression_name}) ===")

    df = generate_data(num_rows, nullable)
    expected_nulls = expected_null_counts(num_rows) if nullable else None
    iters = get_iterations(name)

    for _ in range(WARMUP_ITERATIONS):
        benchmark_write(filename, df, compression, nullable)

    write_times = []
    for _ in range(iters):
        write_times.append(benchmark_write(filename, df, compression, nullable))

    file_size = os.path.getsize(filename)

    # Purge cache once, then warmup reads (first read = cold, rest = warm)
    purge_file_cache(filename)

    for _ in range(WARMUP_ITERATIONS):
        benchmark_read(filename, num_rows, expected_nulls)

    read_times = []
    for _ in range(iters):
        read_times.append(benchmark_read(filename, num_rows, expected_nulls))

    write_med = trimmed_median(write_times)
    read_med = trimmed_median(read_times)
    rows_per_sec_write = (num_rows / write_med) * 1000
    rows_per_sec_read = (num_rows / read_med) * 1000

    print(f"  Write: {write_med:.2f} ms ({rows_per_sec_write / 1e6:.2f} M rows/sec)")
    print(f"  Read:  {read_med:.2f} ms ({rows_per_sec_read / 1e6:.2f} M rows/sec)")
    print(f"  File:  {file_size / (1024 * 1024):.2f} MB ({file_size / num_rows:.2f} bytes/row)")

    print(f"CSV:fastparquet,{name},{compression_name},{num_rows},"
          f"{write_med:.2f},{read_med:.2f},{file_size}")

    os.remove(filename)


def main():
    print("fastparquet Benchmark")
    print(f"fastparquet version: {fastparquet.__version__}")

    for name, rows in [("small", 100_000), ("medium", 1_000_000), ("large", 10_000_000),
                       ("xlarge", 100_000_000)]:
        for comp_name in CODECS:
            run_benchmark(name, rows, comp_name)

    print("\nBenchmark complete.")


if __name__ == "__main__":
    main()
