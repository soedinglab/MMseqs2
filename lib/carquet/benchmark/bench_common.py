"""Shared helpers for the Python Parquet benchmarks (PyArrow, fastparquet, polars).

Mirrors benchmark_carquet.c: same LCG data, iteration counts, trimmed median,
page-cache purge and temp directory, and the null pattern of the
"<size>_nullable" datasets.
"""

import os
import platform
from functools import lru_cache

import numpy as np

WARMUP_ITERATIONS = 3
BENCH_ITERATIONS = {"small": 51, "medium": 21, "large": 11, "xlarge": 11}
LCG_A = np.uint64(1103515245)
LCG_C = np.uint64(12345)
LCG_MASK = np.uint64(0xFFFFFFFF)
LCG_CHUNK_ROWS = 1_000_000
BOX_MULLER_PI = 3.14159265358979
COLUMN_NAMES = ("id", "value", "category")

# "<size>_nullable" datasets: same values, OPTIONAL columns with NULL_PERCENT %
# nulls each, dictionary encoding on `category`. Must match benchmark_carquet.c
# and benchmark_arrow_cpp.cpp.
NULLABLE_SUFFIX = "_nullable"
NULL_PERCENT = 7
DICTIONARY_COLUMNS = ("category",)


def get_benchmark_zstd_level():
    raw = os.getenv("CARQUET_BENCH_ZSTD_LEVEL", "1").strip()
    try:
        level = int(raw)
    except ValueError:
        return 1
    return max(1, min(level, 22))


def trimmed_median(values):
    """Sort, drop min and max, return median of the rest."""
    s = sorted(values)
    n = len(s)
    if n <= 2:
        return s[n // 2]
    trimmed = s[1:-1]
    m = len(trimmed)
    if m % 2 == 1:
        return trimmed[m // 2]
    return (trimmed[m // 2 - 1] + trimmed[m // 2]) / 2


def purge_file_cache(filename):
    """Evict file from OS page cache so reads hit storage."""
    if platform.system() == "Darwin":
        import fcntl, shutil
        tmp = filename + ".nocache"
        # Copy through F_NOCACHE fd to avoid populating cache
        src = os.open(filename, os.O_RDONLY)
        dst = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_TRUNC, 0o644)
        fcntl.fcntl(dst, fcntl.F_NOCACHE, 1)
        fcntl.fcntl(src, fcntl.F_NOCACHE, 1)
        while True:
            buf = os.read(src, 262144)
            if not buf:
                break
            os.write(dst, buf)
        os.close(src)
        os.close(dst)
        os.unlink(filename)
        os.rename(tmp, filename)
    elif platform.system() == "Linux":
        fd = os.open(filename, os.O_RDONLY)
        try:
            size = os.fstat(fd).st_size
            os.posix_fadvise(fd, 0, size, os.POSIX_FADV_DONTNEED)
        finally:
            os.close(fd)


def get_temp_dir():
    override = os.getenv("CARQUET_BENCH_TMPDIR", "").strip()
    if override:
        return override
    return "/tmp"


@lru_cache(maxsize=8)
def _lcg_advance_params(count):
    """Return affine params for advancing the C benchmark LCG by count steps."""
    idx = np.arange(count + 1, dtype=np.uint64)
    mul = np.ones(count + 1, dtype=np.uint64)
    add = np.zeros(count + 1, dtype=np.uint64)

    cur_mul = LCG_A
    cur_add = LCG_C
    bit = 1
    while bit <= count:
        mask = (idx & bit) != 0
        mul_mask = mul[mask]
        add_mask = add[mask]
        mul[mask] = (cur_mul * mul_mask) & LCG_MASK
        add[mask] = (cur_mul * add_mask + cur_add) & LCG_MASK
        cur_add = (cur_add * (cur_mul + np.uint64(1))) & LCG_MASK
        cur_mul = (cur_mul * cur_mul) & LCG_MASK
        bit <<= 1

    return mul[1:], add[1:]


def _lcg_rand_chunk(state, count):
    """Generate count outputs from the benchmark LCG, matching the C code."""
    mul, add = _lcg_advance_params(count)
    states = (mul * np.uint64(state) + add) & LCG_MASK
    values = ((states >> np.uint64(16)) & np.uint64(0x7FFF)).astype(np.uint32)
    return values, int(states[-1]) if count else state


def generate_columns(num_rows):
    """Pre-generate test data using the same RNG sequence as benchmark_carquet.c."""
    ids = np.empty(num_rows, dtype=np.int64)
    values = np.empty(num_rows, dtype=np.float64)
    categories = np.empty(num_rows, dtype=np.int32)

    state = 42
    offset = 0
    while offset < num_rows:
        chunk_rows = min(LCG_CHUNK_ROWS, num_rows - offset)
        rand, state = _lcg_rand_chunk(state, chunk_rows * 4)
        rand = rand.reshape(chunk_rows, 4)

        ids[offset:offset + chunk_rows] = 1_000_000 + (rand[:, 0].astype(np.int64) % 9_000_000)
        u1 = (rand[:, 1].astype(np.float64) + 1.0) / 32768.0
        u2 = (rand[:, 2].astype(np.float64) + 1.0) / 32768.0
        normal = np.sqrt(-2.0 * np.log(u1)) * np.cos(2.0 * BOX_MULLER_PI * u2)
        values[offset:offset + chunk_rows] = np.abs(100.0 + 50.0 * normal)
        categories[offset:offset + chunk_rows] = (rand[:, 3] % 100).astype(np.int32)
        offset += chunk_rows

    return {"id": ids, "value": values, "category": categories}


def is_nullable_dataset(name):
    return name.endswith(NULLABLE_SUFFIX)


def base_dataset(name):
    """Size name of a dataset ("large_nullable" -> "large")."""
    return name[:-len(NULLABLE_SUFFIX)] if is_nullable_dataset(name) else name


def null_masks(num_rows):
    """Per-column null masks (True = null) of the nullable datasets.

    lowbias32 hash of (row ^ per-column seed), identical to bench_is_null() in
    benchmark_carquet.c / benchmark_arrow_cpp.cpp.
    """
    rows = np.arange(num_rows, dtype=np.uint32)
    masks = {}
    for col, name in enumerate(COLUMN_NAMES):
        x = rows ^ np.uint32((0x9E3779B9 * (col + 1)) & 0xFFFFFFFF)
        x ^= x >> np.uint32(16)
        x *= np.uint32(0x7FEB352D)
        x ^= x >> np.uint32(15)
        x *= np.uint32(0x846CA68B)
        x ^= x >> np.uint32(16)
        masks[name] = (x % np.uint32(100)) < np.uint32(NULL_PERCENT)
    return masks


def get_iterations(name):
    iter_override = os.environ.get("CARQUET_BENCH_ITERATIONS")
    if iter_override and int(iter_override) > 0:
        return min(int(iter_override), 201)
    return BENCH_ITERATIONS.get(base_dataset(name), 21)
