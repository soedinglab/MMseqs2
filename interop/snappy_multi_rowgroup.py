#!/usr/bin/env python3
"""
Roundtrip a large, many-row-group, Snappy-compressed PyArrow file through carquet.

This reproduces the shape of a field report where a 3.6 GB / 92-row-group /
36-column Snappy file failed part-way through with
CARQUET_ERROR_INVALID_COMPRESSED_DATA while PyArrow read the same file fine.
The characteristics that matter, and that this script reproduces, are:

  * many row groups, so per-chunk state (dictionary page, first-data-page
    offset) is re-derived dozens of times;
  * absolute page offsets past INT32_MAX, when --rows is large enough;
  * wide BYTE_ARRAY columns, whose per-page statistics push page headers well
    past the reader's initial header read window;
  * high-cardinality string columns, which make the writer fall back from
    dictionary to plain encoding mid-chunk;
  * a mix of dictionary-encoded, plain, and numeric columns in one file.

Usage:
    python snappy_multi_rowgroup.py [--rows N] [--row-groups N] [--keep] [--out PATH]

Verification runs `carquet export` over the whole file (a full decode of every
page of every column) and compares the row count against PyArrow's.
Set CARQUET_BIN to point at the carquet CLI; it defaults to ../build/carquet.
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

NUM_COLUMNS = 36


def build_table(rows: int, rng: np.random.Generator, wide_len: int) -> pa.Table:
    """36 columns: one wide high-entropy string, one long-statistics string,
    then a repeating mix of int64 / double / low-cardinality string."""
    names, cols = [], []

    raw = rng.integers(33, 126, size=rows * wide_len, dtype=np.uint8).tobytes()
    offsets = np.arange(0, (rows + 1) * wide_len, wide_len, dtype=np.int32)
    cols.append(pa.StringArray.from_buffers(rows, pa.py_buffer(offsets), pa.py_buffer(raw)))
    names.append("payload")

    # Long values so the page-header statistics (untruncated min/max) exceed
    # any small fixed header read window.
    long_prefix = "L" * 300
    cols.append(pa.array([long_prefix + str(v) for v in rng.integers(0, 5, rows)],
                         type=pa.string()))
    names.append("longstr")

    for c in range(NUM_COLUMNS - 2):
        if c % 3 == 0:
            cols.append(pa.array(rng.integers(-2**40, 2**40, rows), type=pa.int64()))
        elif c % 3 == 1:
            cols.append(pa.array(rng.normal(size=rows)))
        else:
            cols.append(pa.array([f"cat_{v}" for v in rng.integers(0, 50, rows)],
                                 type=pa.string()))
        names.append(f"c{c}")

    return pa.table(cols, names=names)


def write_file(path: Path, row_groups: int, rows: int, wide_len: int) -> None:
    rng = np.random.default_rng(20260731)
    writer = None
    try:
        for i in range(row_groups):
            table = build_table(rows, rng, wide_len)
            if writer is None:
                writer = pq.ParquetWriter(
                    path, table.schema,
                    compression="snappy",
                    data_page_size=16 * 1024 * 1024,
                    dictionary_pagesize_limit=1024 * 1024,
                    write_statistics=True,
                    write_page_index=False,
                    version="1.0",
                )
            writer.write_table(table, row_group_size=rows)
            print(f"  row group {i + 1}/{row_groups} "
                  f"({os.path.getsize(path) / 1e9:.3f} GB)", flush=True)
    finally:
        if writer is not None:
            writer.close()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", type=int, default=50_000,
                    help="rows per row group (default 50000)")
    ap.add_argument("--row-groups", type=int, default=8,
                    help="number of row groups (default 8)")
    ap.add_argument("--wide-len", type=int, default=40,
                    help="byte length of the wide string column (default 40)")
    ap.add_argument("--out", type=Path, default=None,
                    help="output file (default: a temporary file)")
    ap.add_argument("--keep", action="store_true",
                    help="keep the generated file")
    args = ap.parse_args()

    carquet_bin = Path(os.environ.get(
        "CARQUET_BIN", Path(__file__).resolve().parent.parent / "build" / "carquet"))
    if not carquet_bin.exists():
        print(f"carquet CLI not found at {carquet_bin}; set CARQUET_BIN", file=sys.stderr)
        return 2

    tmpdir = None
    if args.out is None:
        tmpdir = tempfile.TemporaryDirectory()
        out = Path(tmpdir.name) / "snappy_multi_rowgroup.parquet"
    else:
        out = args.out

    try:
        print(f"Writing {args.row_groups} row groups x {args.rows} rows "
              f"x {NUM_COLUMNS} columns (snappy) to {out}")
        write_file(out, args.row_groups, args.rows, args.wide_len)

        meta = pq.ParquetFile(out).metadata
        print(f"PyArrow: {meta.num_rows} rows, {meta.num_columns} columns, "
              f"{meta.num_row_groups} row groups, {os.path.getsize(out)} bytes")

        # Full decode of every page of every column.
        print(f"carquet export {out} ...")
        proc = subprocess.run([str(carquet_bin), "export", str(out)],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if proc.returncode != 0:
            print("FAIL: carquet export exited "
                  f"{proc.returncode}: {proc.stderr.decode(errors='replace')[-4000:]}",
                  file=sys.stderr)
            return 1

        # Header row plus one line per row (values contain no newlines here).
        exported = proc.stdout.count(b"\n") - 1
        if exported != meta.num_rows:
            print(f"FAIL: carquet exported {exported} rows, expected {meta.num_rows}",
                  file=sys.stderr)
            return 1

        print(f"OK: carquet decoded all {exported} rows of {out.name}")
        return 0
    finally:
        if tmpdir is not None and not args.keep:
            tmpdir.cleanup()


if __name__ == "__main__":
    sys.exit(main())
