#!/usr/bin/env python3
"""
Benchmark for polars' native Parquet reader / writer - compare with Carquet.

polars reads and writes Parquet with its own Rust implementation
(polars-parquet), not through Arrow C++ / PyArrow: `use_pyarrow` stays False.

Same data, iteration counts, trimmed median and cache purge as
benchmark_carquet.c / benchmark_pyarrow.py. Settings are aligned where polars
allows it (same row group size, zstd level, page size and LZ4_RAW for "lz4",
statistics on); polars cannot write BYTE_STREAM_SPLIT or page checksums and
picks dictionary encoding itself (no per-column switch), so its files differ
from the others. Its native reader has no memory-map switch (`memory_map` only
applies to `use_pyarrow=True`), so no read_nommap figure is reported.

"<size>_nullable" datasets carry the shared null pattern (bench_common.null_masks).
"""

import os
import time

import numpy as np
import polars as pl

from bench_common import (
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

# polars' "lz4" writes LZ4_RAW, like carquet's.
CODECS = {"none": "uncompressed", "snappy": "snappy", "zstd": "zstd", "lz4": "lz4"}


def generate_data(num_rows, nullable=False):
    """Pre-generate test data using the same RNG sequence as benchmark_carquet.c."""
    columns = generate_columns(num_rows)
    if not nullable:
        return pl.DataFrame(columns)
    masks = null_masks(num_rows)
    return pl.DataFrame([
        pl.Series(name, values).scatter(np.flatnonzero(masks[name]), None)
        for name, values in columns.items()
    ])


def expected_null_counts(num_rows):
    return {name: int(mask.sum()) for name, mask in null_masks(num_rows).items()}


def benchmark_write(filename, df, compression_name):
    """Write a Parquet file and return time_ms."""
    start = time.perf_counter()
    num_rows = df.height
    rg_size = num_rows // 10 if num_rows > 1_000_000 else 100_000
    kwargs = {"compression": CODECS[compression_name], "row_group_size": rg_size}
    if compression_name == "zstd":
        kwargs["compression_level"] = get_benchmark_zstd_level()
    df.write_parquet(filename, **kwargs)
    return (time.perf_counter() - start) * 1000


def benchmark_read(filename, expected_rows, expected_nulls=None):
    """Read a Parquet file and return time_ms."""
    start = time.perf_counter()
    df = pl.read_parquet(filename)

    nulls = None
    checksum = 0
    if expected_nulls is None:
        for col_name in df.columns:
            checksum += df[col_name].to_numpy()[0]
    else:
        nulls = {col_name: df[col_name].null_count() for col_name in df.columns}

    elapsed_ms = (time.perf_counter() - start) * 1000

    if df.height != expected_rows:
        raise ValueError(f"Row count mismatch: {df.height} vs {expected_rows}")
    if expected_nulls is not None and nulls != expected_nulls:
        raise ValueError(f"Null count mismatch: {nulls} vs {expected_nulls}")

    return elapsed_ms


def run_benchmark(name, num_rows, compression_name):
    """Run a single benchmark configuration."""
    filename = os.path.join(
        get_temp_dir(),
        f"benchmark_{name}_{compression_name}_polars.parquet",
    )
    nullable = is_nullable_dataset(name)

    print(f"\n=== {name} ({num_rows:,} rows, {compression_name}) ===")

    df = generate_data(num_rows, nullable)
    expected_nulls = expected_null_counts(num_rows) if nullable else None
    iters = get_iterations(name)

    for _ in range(WARMUP_ITERATIONS):
        benchmark_write(filename, df, compression_name)

    write_times = []
    for _ in range(iters):
        write_times.append(benchmark_write(filename, df, compression_name))

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

    print(f"CSV:polars,{name},{compression_name},{num_rows},"
          f"{write_med:.2f},{read_med:.2f},{file_size}")

    os.remove(filename)


def main():
    print("polars Benchmark")
    print(f"polars version: {pl.__version__}")

    for name, rows in [("small", 100_000), ("medium", 1_000_000), ("large", 10_000_000),
                       ("xlarge", 100_000_000)]:
        for comp_name in CODECS:
            run_benchmark(name, rows, comp_name)

    print("\nBenchmark complete.")


if __name__ == "__main__":
    main()
