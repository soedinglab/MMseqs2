#!/usr/bin/env python3
"""
Benchmark for PyArrow Parquet - compare with Carquet.

Methodology mirrors benchmark_carquet.c:
  - Write and read benchmarked separately
  - Page cache purged between write and read phases
  - 3 warmup + 11 measured iterations, trimmed median (drop min + max)
  - Data pre-generated outside timing loops
  - Arrow write/read options aligned to the Carquet benchmark where supported
  - Reads timed with memory_map=True ("read") and memory_map=False ("read_nommap")
  - "<size>_nullable" datasets: nullable columns with the shared null pattern
    (bench_common.null_masks) and dictionary encoding on `category` only
"""

import os
import time

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

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


def generate_data(num_rows, nullable=False):
    """Pre-generate test data using the same RNG sequence as benchmark_carquet.c."""
    columns = generate_columns(num_rows)
    if not nullable:
        return pa.table(columns)
    masks = null_masks(num_rows)
    return pa.table({name: pa.array(values, mask=masks[name])
                     for name, values in columns.items()})


def expected_null_counts(num_rows):
    return {name: int(mask.sum()) for name, mask in null_masks(num_rows).items()}


def benchmark_write(filename, table, compression, nullable=False):
    """Write a Parquet file and return time_ms."""
    start = time.perf_counter()
    num_rows = len(table)
    rg_size = num_rows // 10 if num_rows > 1_000_000 else 100_000
    kwargs = {
        "compression": compression,
        "row_group_size": rg_size,
        "use_dictionary": list(DICTIONARY_COLUMNS) if nullable else False,
        "store_schema": False,
        "write_page_checksum": True,
    }
    if compression is not None:
        kwargs["use_byte_stream_split"] = ["value"]
    if compression == "zstd":
        kwargs["compression_level"] = get_benchmark_zstd_level()
    pq.write_table(table, filename, **kwargs)
    return (time.perf_counter() - start) * 1000


def touch_nullable(table):
    """Per-column null counts from the validity bitmaps, plus first/last
    value of each chunk (like the C loops). to_numpy() would turn nullable
    ints into float64 + NaN copies, so the value buffers are viewed instead."""
    checksum = 0
    nulls = {}
    for col_name in table.column_names:
        col_nulls = 0
        for chunk in table[col_name].chunks:
            n = len(chunk)
            if n == 0:
                continue
            col_nulls += n - pc.is_valid(chunk).true_count
            values = np.frombuffer(chunk.buffers()[1], dtype=chunk.type.to_pandas_dtype())
            checksum += values[chunk.offset] + values[chunk.offset + n - 1]
        nulls[col_name] = col_nulls
    return checksum, nulls


def benchmark_read(filename, expected_rows, memory_map=True, expected_nulls=None):
    """Read a Parquet file and return time_ms."""
    start = time.perf_counter()
    table = pq.read_table(
        filename,
        memory_map=memory_map,
        use_threads=True,
        page_checksum_verification=False,
    )

    nulls = None
    checksum = 0
    if expected_nulls is None:
        for col_name in table.column_names:
            arr = table[col_name].to_numpy()
            checksum += arr[0]
    else:
        checksum, nulls = touch_nullable(table)

    elapsed_ms = (time.perf_counter() - start) * 1000

    if len(table) != expected_rows:
        raise ValueError(f"Row count mismatch: {len(table)} vs {expected_rows}")
    if expected_nulls is not None and nulls != expected_nulls:
        raise ValueError(f"Null count mismatch: {nulls} vs {expected_nulls}")

    return elapsed_ms


def run_benchmark(name, num_rows, compression, compression_name):
    """Run a single benchmark configuration."""
    filename = os.path.join(
        get_temp_dir(),
        f"benchmark_{name}_{compression_name}_pyarrow.parquet",
    )

    print(f"\n=== {name} ({num_rows:,} rows, {compression_name}) ===")

    nullable = is_nullable_dataset(name)
    table = generate_data(num_rows, nullable)
    expected_nulls = expected_null_counts(num_rows) if nullable else None

    iters = get_iterations(name)

    # Warmup writes
    for _ in range(WARMUP_ITERATIONS):
        benchmark_write(filename, table, compression, nullable)

    # Benchmark writes
    write_times = []
    for _ in range(iters):
        write_times.append(benchmark_write(filename, table, compression, nullable))

    file_size = os.path.getsize(filename)

    # Purge cache once, then warmup reads (first read = cold, rest = warm)
    purge_file_cache(filename)

    for _ in range(WARMUP_ITERATIONS):
        benchmark_read(filename, num_rows, True, expected_nulls)

    # Benchmark reads (warm cache — realistic for most workloads)
    read_times = []
    for _ in range(iters):
        read_times.append(benchmark_read(filename, num_rows, True, expected_nulls))

    # Same protocol without memory mapping (pyarrow's default), after the
    # mmap pass so that one is unchanged.
    purge_file_cache(filename)

    for _ in range(WARMUP_ITERATIONS):
        benchmark_read(filename, num_rows, False, expected_nulls)

    read_nommap_times = []
    for _ in range(iters):
        read_nommap_times.append(benchmark_read(filename, num_rows, False, expected_nulls))

    write_med = trimmed_median(write_times)
    read_med = trimmed_median(read_times)
    read_nommap_med = trimmed_median(read_nommap_times)
    rows_per_sec_write = (num_rows / write_med) * 1000
    rows_per_sec_read = (num_rows / read_med) * 1000

    print(f"  Write: {write_med:.2f} ms ({rows_per_sec_write / 1e6:.2f} M rows/sec)")
    print(f"  Read:  {read_med:.2f} ms ({rows_per_sec_read / 1e6:.2f} M rows/sec)")
    print(f"  Read (no mmap): {read_nommap_med:.2f} ms "
          f"({(num_rows / read_nommap_med) * 1000 / 1e6:.2f} M rows/sec)")
    print(f"  File:  {file_size / (1024 * 1024):.2f} MB ({file_size / num_rows:.2f} bytes/row)")

    # 8th field (read_nommap) is optional for the orchestrator's parser.
    print(f"CSV:pyarrow,{name},{compression_name},{num_rows},"
          f"{write_med:.2f},{read_med:.2f},{file_size},{read_nommap_med:.2f}")

    os.remove(filename)


def main():
    print("PyArrow Benchmark")
    print(f"PyArrow version: {pa.__version__}")

    for name, rows in [("small", 100_000), ("medium", 1_000_000), ("large", 10_000_000),
                        ("xlarge", 100_000_000)]:
        for comp_name, comp in [("none", None), ("snappy", "snappy"), ("zstd", "zstd"), ("lz4", "lz4")]:
            run_benchmark(name, rows, comp, comp_name)

    print("\nBenchmark complete.")


if __name__ == "__main__":
    main()
