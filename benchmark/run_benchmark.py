#!/usr/bin/env python3
"""
Carquet benchmark orchestrator.

Runs each configuration independently with cooldown pauses between them
to avoid thermal throttling on fanless laptops (e.g. MacBook Air).

Usage:
    python3 benchmark/run_benchmark.py                # full benchmark
    python3 benchmark/run_benchmark.py --quick         # large only (10M rows)
    python3 benchmark/run_benchmark.py --skip-xlarge   # skip 100M row configs
    python3 benchmark/run_benchmark.py --skip-small --skip-medium  # large+ only
    python3 benchmark/run_benchmark.py --iterations 51  # 51 iters for stable numbers
    python3 benchmark/run_benchmark.py --cooldown 5    # 5s cooldown between configs
    python3 benchmark/run_benchmark.py --no-pyarrow    # skip PyArrow
    python3 benchmark/run_benchmark.py --no-arrow-cpp  # skip Arrow C++
"""

import argparse
import datetime
import json
import os
import platform
import subprocess
import sys
import time

# ── ANSI helpers ─────────────────────────────────────────────────────────────

BOLD = "\033[1m"
DIM = "\033[2m"
RST = "\033[0m"
RED = "\033[31m"
GRN = "\033[32m"
YLW = "\033[33m"
BLU = "\033[34m"
CYN = "\033[36m"
WHT = "\033[37m"

# Disable colors if not a terminal
if not sys.stdout.isatty():
    BOLD = DIM = RST = RED = GRN = YLW = BLU = CYN = WHT = ""

def _bar(char="─", width=62):
    return DIM + char * width + RST


# ── System info ──────────────────────────────────────────────────────────────

def get_cpu_name():
    """Get CPU model string."""
    if platform.system() == "Darwin":
        try:
            out = subprocess.check_output(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                stderr=subprocess.DEVNULL, text=True
            ).strip()
            if out:
                return out
        except Exception:
            pass
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return platform.processor() or "unknown"


def get_mem_gb():
    """Get total RAM in GB."""
    if platform.system() == "Darwin":
        try:
            out = subprocess.check_output(
                ["sysctl", "-n", "hw.memsize"],
                stderr=subprocess.DEVNULL, text=True
            ).strip()
            return int(out) / (1024 ** 3)
        except Exception:
            pass
    return 0


def get_system_info():
    """Collect machine specs as a dict."""
    return {
        "cpu": get_cpu_name(),
        "ram_gb": round(get_mem_gb()),
        "os": platform.platform(),
        "arch": platform.machine(),
        "python": platform.python_version(),
    }


def get_benchmark_version(binary):
    """Get version string from a benchmark binary."""
    try:
        out = subprocess.check_output(
            [binary, "--version"], stderr=subprocess.DEVNULL, text=True
        ).strip()
        return out
    except Exception:
        return "unknown"


def print_system_info(info):
    print(f"  {DIM}CPU:{RST}  {info['cpu']}")
    if info["ram_gb"] > 0:
        print(f"  {DIM}RAM:{RST}  {info['ram_gb']} GB")
    print(f"  {DIM}OS:{RST}   {info['os']}")


# ── Benchmark runner ─────────────────────────────────────────────────────────

CONFIGS = [
    # (dataset, compression, rows) — largest first for best thermal conditions
    ("xlarge", "none",   100_000_000),
    ("xlarge", "snappy", 100_000_000),
    ("xlarge", "zstd",   100_000_000),
    ("xlarge", "lz4",    100_000_000),
    ("large",  "none",    10_000_000),
    ("large",  "snappy",  10_000_000),
    ("large",  "zstd",    10_000_000),
    ("large",  "lz4",     10_000_000),
    ("medium", "none",     1_000_000),
    ("medium", "snappy",   1_000_000),
    ("medium", "zstd",     1_000_000),
    ("medium", "lz4",      1_000_000),
    ("small",  "none",       100_000),
    ("small",  "snappy",     100_000),
    ("small",  "zstd",       100_000),
    ("small",  "lz4",        100_000),
]
DATASET_ORDER = ["xlarge", "large", "medium", "small"]
COMPRESSION_ORDER = []
for _dataset, _compression, _rows in CONFIGS:
    if _compression not in COMPRESSION_ORDER:
        COMPRESSION_ORDER.append(_compression)


def get_benchmark_zstd_level():
    raw = os.getenv("CARQUET_BENCH_ZSTD_LEVEL", "1").strip()
    try:
        level = int(raw)
    except ValueError:
        return 1
    return max(1, min(level, 22))


def parse_csv_line(line):
    """Parse a CSV: prefixed output line from benchmarks."""
    if not line.startswith("CSV:"):
        return None
    parts = line[4:].strip().split(",")
    return {
        "library":     parts[0],
        "dataset":     parts[1],
        "compression": parts[2],
        "rows":        int(parts[3]),
        "write_ms":    float(parts[4]),
        "read_ms":     float(parts[5]),
        "file_bytes":  int(parts[6]),
    }


def result_is_valid(result):
    if not result:
        return False
    return (
        result["rows"] > 0 and
        result["write_ms"] > 0 and
        result["read_ms"] > 0 and
        result["file_bytes"] > 0
    )


def stderr_tail(stderr):
    if not stderr:
        return None
    lines = [line.strip() for line in stderr.splitlines() if line.strip()]
    if not lines:
        return None
    return lines[-1]


def run_benchmark_process(argv, timeout_seconds):
    env = os.environ.copy()
    try:
        proc = subprocess.run(
            argv,
            capture_output=True, text=True, timeout=timeout_seconds, env=env
        )
    except subprocess.TimeoutExpired:
        return None, f"timed out after {timeout_seconds}s"
    parsed = None
    for line in proc.stdout.splitlines():
        r = parse_csv_line(line)
        if r:
            parsed = r
            break

    if proc.returncode != 0:
        return None, stderr_tail(proc.stderr) or f"exit code {proc.returncode}"
    if not result_is_valid(parsed):
        return None, stderr_tail(proc.stderr) or "benchmark returned invalid or zero timings"
    return parsed, None


def run_carquet_single(binary, dataset, compression, timeout_seconds):
    """Run carquet benchmark for a single config."""
    return run_benchmark_process([binary, dataset, compression], timeout_seconds)


def run_arrow_cpp_single(binary, dataset, compression, timeout_seconds):
    """Run Arrow C++ benchmark for a single config."""
    return run_benchmark_process([binary, dataset, compression], timeout_seconds)


def run_pyarrow_single(python, script, dataset, compression, rows, timeout_seconds):
    """Run pyarrow benchmark for a single config by invoking a helper snippet."""
    code = f"""
import sys, os
sys.path.insert(0, os.path.dirname({script!r}))
from benchmark_pyarrow import run_benchmark
comp_map = {{"none": None, "snappy": "snappy", "zstd": "zstd", "lz4": "lz4"}}
run_benchmark({dataset!r}, {rows}, comp_map[{compression!r}], {compression!r})
"""
    return run_benchmark_process([python, "-c", code], timeout_seconds)


def crossread_result_is_valid(result):
    """Validate a cross-read result (write_ms is 0, only read matters)."""
    if not result:
        return False
    return result["rows"] > 0 and result["read_ms"] > 0


def run_crossread_process(argv, timeout_seconds):
    """Run a cross-read benchmark subprocess (write_ms=0 is expected)."""
    env = os.environ.copy()
    try:
        proc = subprocess.run(
            argv,
            capture_output=True, text=True, timeout=timeout_seconds, env=env
        )
    except subprocess.TimeoutExpired:
        return None, f"timed out after {timeout_seconds}s"
    parsed = None
    for line in proc.stdout.splitlines():
        r = parse_csv_line(line)
        if r:
            parsed = r
            break

    if proc.returncode != 0:
        return None, stderr_tail(proc.stderr) or f"exit code {proc.returncode}"
    if not crossread_result_is_valid(parsed):
        return None, stderr_tail(proc.stderr) or "cross-read returned invalid timings"
    return parsed, None


def run_carquet_read_file(binary, filename, timeout_seconds):
    """Run carquet benchmark in read-only mode on an external file."""
    return run_crossread_process([binary, "read", filename], timeout_seconds)


def run_arrow_cpp_read_file(binary, filename, timeout_seconds):
    """Run Arrow C++ benchmark in read-only mode on an external file."""
    return run_crossread_process([binary, "read", filename], timeout_seconds)


def cooldown(seconds):
    """Visual cooldown timer."""
    if seconds <= 0:
        return
    for i in range(seconds, 0, -1):
        print(f"\r  {DIM}cooldown {i}s ...{RST}", end="", flush=True)
        time.sleep(1)
    print(f"\r{' ':40}\r", end="", flush=True)


# ── Result formatting ────────────────────────────────────────────────────────

def fmt_ms(ms):
    if ms is None:
        return "    -   "
    if ms >= 100:
        return f"{ms:>7.1f}ms"
    return f"{ms:>7.2f}ms"


def fmt_throughput(rows, ms):
    if ms is None or ms == 0:
        return ""
    mrows = (rows / ms) * 1000 / 1e6
    return f"{DIM}({mrows:>.0f}M r/s){RST}"


def fmt_speedup(carquet_ms, pyarrow_ms):
    if carquet_ms is None or pyarrow_ms is None or carquet_ms == 0:
        return "       "
    ratio = pyarrow_ms / carquet_ms
    if ratio >= 1.5:
        color = GRN
    elif ratio >= 1.0:
        color = WHT
    elif ratio >= 0.8:
        color = YLW
    else:
        color = RED
    return f"{color}{BOLD}{ratio:>5.2f}x{RST}"


def fmt_size(b):
    mb = b / (1024 * 1024)
    if mb >= 10:
        return f"{mb:>6.1f}MB"
    return f"{mb:>6.2f}MB"


def have_results(results, library):
    return any(lib == library for lib, _ in results.keys())


def print_library_table(results, library, label):
    """Print a single-library summary table."""
    print()
    print(f"  {BOLD}{label}{RST}")
    print(f"  {BOLD}{'':15} {'Write':^18}{'':2}{'Read':^18}{'':2}{'Size':^10}{RST}")
    print(f"  {'':15} {'time':>8} {'M r/s':>8}"
          f"  {'time':>8} {'M r/s':>8}"
          f"  {'bytes':>8}")
    print(f"  {_bar()}")

    for dataset_idx, dataset in enumerate(DATASET_ORDER):
        printed_dataset = False
        for compression in COMPRESSION_ORDER:
            key = (dataset, compression)
            result = results.get((library, key))
            if not result:
                continue

            printed_dataset = True
            label = f"{dataset}/{compression}"
            rows = result["rows"]
            line = f"  {BOLD}{label:<15}{RST}"
            def throughput(r, ms):
                if ms and ms > 0:
                    return f"{(r/ms)*1000/1e6:>7.1f}M"
                return "      -"
            line += f" {fmt_ms(result['write_ms'])} {throughput(rows, result['write_ms'])}"
            line += f"  {fmt_ms(result['read_ms'])} {throughput(rows, result['read_ms'])}"
            line += f"  {fmt_size(result['file_bytes']):>8}"
            print(line)
        has_later = any(results.get((library, (later, comp)))
                        for later in DATASET_ORDER[dataset_idx + 1:]
                        for comp in COMPRESSION_ORDER)
        if printed_dataset and has_later:
            print()

def print_comparison_table(results, competitor, competitor_label):
    """Print a Carquet-vs-competitor comparison table."""
    print()
    print(f"  {BOLD}Carquet vs {competitor_label}{RST}")
    print(f"  {BOLD}{'':15} {'Write':^24}{'':2}{'Read':^24}{'':2}{'Size':^16}{RST}")
    print(f"  {'':15} {'Carquet':>8} {competitor_label:>8} {'ratio':>6}"
          f"  {'Carquet':>8} {competitor_label:>8} {'ratio':>6}"
          f"  {'Carquet':>7} {competitor_label:>7}")
    print(f"  {_bar()}")

    for dataset_idx, dataset in enumerate(DATASET_ORDER):
        printed_dataset = False
        for compression in COMPRESSION_ORDER:
            key = (dataset, compression)
            carquet = results.get(("carquet", key))
            other = results.get((competitor, key))
            if not carquet and not other:
                continue

            printed_dataset = True
            label = f"{dataset}/{compression}"
            c_w = carquet["write_ms"] if carquet else None
            c_r = carquet["read_ms"] if carquet else None
            c_sz = carquet["file_bytes"] if carquet else None
            o_w = other["write_ms"] if other else None
            o_r = other["read_ms"] if other else None
            o_sz = other["file_bytes"] if other else None

            line = f"  {BOLD}{label:<15}{RST}"
            line += f" {fmt_ms(c_w)} {fmt_ms(o_w)} {fmt_speedup(c_w, o_w)}"
            line += f"  {fmt_ms(c_r)} {fmt_ms(o_r)} {fmt_speedup(c_r, o_r)}"
            line += f"  {(fmt_size(c_sz) if c_sz else ''):>7}"
            line += f" {(fmt_size(o_sz) if o_sz else ''):>7}"
            print(line)
        has_later = any(results.get(("carquet", (later, comp))) or
                        results.get((competitor, (later, comp)))
                        for later in DATASET_ORDER[dataset_idx + 1:]
                        for comp in COMPRESSION_ORDER)
        if printed_dataset and has_later:
            print()

    print()
    print(f"  {DIM}ratio = {competitor_label} time / Carquet time (higher = Carquet faster){RST}")


def print_crossread_table(crossread_results):
    """Print a same-file cross-read comparison table."""
    if not crossread_results:
        return

    print()
    print(f"  {BOLD}Same-file cross-read comparison{RST}")
    print(f"  {DIM}(both libraries read the SAME Parquet file){RST}")
    print(f"  {BOLD}{'':15} {'Writer':>10} {'Carquet R':>10} {'Arrow R':>10} {'ratio':>7}{RST}")
    print(f"  {_bar()}")

    for dataset_idx, dataset in enumerate(DATASET_ORDER):
        printed_dataset = False
        for compression in COMPRESSION_ORDER:
            key = (dataset, compression)

            # Carquet-written file results
            cw = crossread_results.get(("carquet_file", key))
            # Arrow-written file results
            aw = crossread_results.get(("arrow_file", key))

            if cw:
                printed_dataset = True
                label = f"{dataset}/{compression}"
                c_r = cw.get("carquet_read_ms")
                a_r = cw.get("arrow_read_ms")
                line = f"  {BOLD}{label:<15}{RST}"
                line += f" {'Carquet':>10}"
                line += f" {fmt_ms(c_r)}"
                line += f" {fmt_ms(a_r)}"
                line += f" {fmt_speedup(c_r, a_r)}"
                print(line)

            if aw:
                printed_dataset = True
                label = f"{dataset}/{compression}"
                c_r = aw.get("carquet_read_ms")
                a_r = aw.get("arrow_read_ms")
                line = f"  {BOLD}{label:<15}{RST}"
                line += f" {'Arrow':>10}"
                line += f" {fmt_ms(c_r)}"
                line += f" {fmt_ms(a_r)}"
                line += f" {fmt_speedup(c_r, a_r)}"
                print(line)

            has_later = any(
                crossread_results.get(("carquet_file", (later, comp))) or
                crossread_results.get(("arrow_file", (later, comp)))
                for later in DATASET_ORDER[dataset_idx + 1:]
                for comp in COMPRESSION_ORDER
            )
            if printed_dataset and has_later and (cw or aw):
                # Only add spacing between dataset groups
                pass
        if printed_dataset:
            has_later = any(
                crossread_results.get(("carquet_file", (later, comp))) or
                crossread_results.get(("arrow_file", (later, comp)))
                for later in DATASET_ORDER[dataset_idx + 1:]
                for comp in COMPRESSION_ORDER
            )
            if has_later:
                print()

    print()
    print(f"  {DIM}ratio = Arrow read time / Carquet read time (higher = Carquet faster){RST}")


def print_final_tables(results, crossread_results=None):
    """Print final summary tables for the libraries that ran."""
    has_carquet = have_results(results, "carquet")
    has_arrow_cpp = have_results(results, "arrow_cpp")
    has_pyarrow = have_results(results, "pyarrow")

    if has_carquet:
        if has_arrow_cpp:
            print_comparison_table(results, "arrow_cpp", "Arrow C++")
        if has_pyarrow:
            print_comparison_table(results, "pyarrow", "PyArrow")
        if not has_arrow_cpp and not has_pyarrow:
            print_library_table(results, "carquet", "Carquet")
    else:
        if has_arrow_cpp:
            print_library_table(results, "arrow_cpp", "Arrow C++")
        if has_pyarrow:
            print_library_table(results, "pyarrow", "PyArrow")

    if crossread_results:
        print_crossread_table(crossread_results)


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Carquet benchmark runner with optional Arrow C++ and PyArrow comparisons"
    )
    parser.add_argument("--quick", action="store_true",
                        help="Run only large configs (10M rows, all codecs)")
    parser.add_argument("--skip-xlarge", action="store_true",
                        help="Skip xlarge (100M rows) configs")
    parser.add_argument("--skip-large", action="store_true",
                        help="Skip large (10M rows) configs")
    parser.add_argument("--skip-medium", action="store_true",
                        help="Skip medium (1M rows) configs")
    parser.add_argument("--skip-small", action="store_true",
                        help="Skip small (100K rows) configs")
    parser.add_argument("--cooldown", type=int, default=3,
                        help="Seconds between configs (default: 3)")
    parser.add_argument("--timeout", type=int, default=300,
                        help="Per-library subprocess timeout in seconds (default: 300)")
    parser.add_argument("--temp-dir",
                        help="Directory for temporary benchmark parquet files")
    parser.add_argument("--no-arrow-cpp", action="store_true",
                        help="Skip Arrow C++ benchmarks")
    parser.add_argument("--no-pyarrow", action="store_true",
                        help="Skip PyArrow benchmarks")
    parser.add_argument("--no-carquet", action="store_true",
                        help="Skip Carquet benchmarks")
    parser.add_argument("--iterations", type=int, default=0,
                        help="Override iteration count for all sizes (default: auto)")
    parser.add_argument("--no-crossread", action="store_true",
                        help="Skip same-file cross-read benchmarks")
    args = parser.parse_args()

    if args.iterations > 0:
        os.environ["CARQUET_BENCH_ITERATIONS"] = str(args.iterations)

    if args.temp_dir:
        temp_dir = os.path.abspath(args.temp_dir)
        if not os.path.isdir(temp_dir):
            print(f"{RED}Error:{RST} temp dir does not exist: {temp_dir}")
            sys.exit(1)
        os.environ["CARQUET_BENCH_TMPDIR"] = temp_dir

    # Locate binaries
    script_dir = os.path.dirname(os.path.abspath(__file__))
    build_dir = os.path.join(script_dir, "..", "build")
    carquet_bin = os.path.join(build_dir, "benchmark_carquet")
    arrow_cpp_bin = os.path.join(build_dir, "benchmark_arrow_cpp")

    if not args.no_carquet and not os.path.isfile(carquet_bin):
        print(f"{RED}Error:{RST} {carquet_bin} not found. Build the benchmark first:")
        print(f"  cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_BUILD_BENCHMARKS=ON -DCARQUET_NATIVE_ARCH=ON && cmake --build build")
        print(f"  (CARQUET_BUILD_BENCHMARKS is OFF by default; CARQUET_NATIVE_ARCH=ON gives the fastest, most representative numbers)")
        sys.exit(1)

    # Collect metadata
    sys_info = get_system_info()
    carquet_version = get_benchmark_version(carquet_bin) if not args.no_carquet else "unknown"
    timestamp = datetime.date.today().isoformat()

    has_arrow_cpp = not args.no_arrow_cpp and os.path.isfile(arrow_cpp_bin)
    arrow_cpp_version = get_benchmark_version(arrow_cpp_bin) if has_arrow_cpp else None
    if not args.no_arrow_cpp and not has_arrow_cpp:
        print(f"{YLW}Warning:{RST} Arrow C++ benchmark not built, skipping")
        print("         Build with -DCARQUET_BUILD_BENCHMARKS=ON "
              "-DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON")
        print("         (requires Arrow/Parquet C++ libs; on macOS: brew install apache-arrow,")
        print("          or point CMake at a prefix with -DCARQUET_ARROW_CPP_ROOT=/path)")

    # Find Python with pyarrow
    python = sys.executable
    has_pyarrow = not args.no_pyarrow
    pyarrow_script = os.path.join(script_dir, "benchmark_pyarrow.py")
    pyarrow_version = None

    if has_pyarrow:
        try:
            out = subprocess.check_output(
                [python, "-c", "import pyarrow; print(pyarrow.__version__)"],
                stderr=subprocess.DEVNULL, text=True
            ).strip()
            pyarrow_version = out
        except Exception:
            print(f"{YLW}Warning:{RST} PyArrow not available, skipping")
            has_pyarrow = False

    benchmark_runners = []
    if not args.no_carquet:
        benchmark_runners.append((
            "carquet",
            "Carquet",
            lambda dataset, compression, rows: run_carquet_single(
                carquet_bin, dataset, compression, args.timeout)
        ))
    if has_arrow_cpp:
        benchmark_runners.append((
            "arrow_cpp",
            "Arrow C++",
            lambda dataset, compression, rows: run_arrow_cpp_single(
                arrow_cpp_bin, dataset, compression, args.timeout)
        ))
    if has_pyarrow:
        benchmark_runners.append((
            "pyarrow",
            "PyArrow",
            lambda dataset, compression, rows: run_pyarrow_single(
                python, pyarrow_script, dataset, compression, rows, args.timeout)
        ))

    if not benchmark_runners:
        print(f"{RED}Error:{RST} no benchmarks are available to run")
        sys.exit(1)

    # Select configs
    configs = CONFIGS
    if args.quick:
        configs = [c for c in configs if c[0] == "large"]
    else:
        skip = set()
        if args.skip_xlarge:
            skip.add("xlarge")
        if args.skip_large:
            skip.add("large")
        if args.skip_medium:
            skip.add("medium")
        if args.skip_small:
            skip.add("small")
        if skip:
            configs = [c for c in configs if c[0] not in skip]

    # ── Header ──
    print()
    print(f"  {BOLD}Carquet {carquet_version} Benchmark{RST}")
    print(f"  {_bar()}")
    print_system_info(sys_info)
    if has_arrow_cpp:
        print(f"  {DIM}vs:{RST}   Arrow C++ {arrow_cpp_version}")
    if has_pyarrow:
        print(f"  {DIM}vs:{RST}   PyArrow {pyarrow_version}")
    zstd_level = get_benchmark_zstd_level()
    iter_desc = f"{args.iterations} iters (override)" if args.iterations > 0 else "11-51 iters by size"
    print(f"  {DIM}cfg:{RST}  {len(configs)} configs, {args.cooldown}s cooldown, "
          f"trimmed median, cache purged, {iter_desc}")
    print(f"  {DIM}zstd:{RST} level {zstd_level}")
    print(f"  {DIM}timeout:{RST} {args.timeout}s per library/config")
    if args.temp_dir:
        print(f"  {DIM}tmp:{RST}  {os.environ['CARQUET_BENCH_TMPDIR']}")
    print(f"  {_bar()}")
    print()

    # ── Initial cooldown ──
    print(f"  {DIM}Letting system settle...{RST}")
    cooldown(args.cooldown)
    print()

    # ── Determine temp directory for file path construction ──
    bench_tmpdir = os.environ.get("CARQUET_BENCH_TMPDIR", "/tmp")

    # Cross-read is only possible when both carquet and arrow_cpp are enabled
    do_crossread = (not args.no_crossread and not args.no_carquet
                    and has_arrow_cpp)

    # ── Run benchmarks ──
    results = {}
    crossread_results = {}
    total = len(configs)
    inter_library_cooldown = 0 if args.cooldown <= 0 else max(1, args.cooldown // 2)

    for i, (dataset, compression, rows) in enumerate(configs):
        key = (dataset, compression)
        tag = f"[{i+1}/{total}]"
        label = f"{dataset}/{compression}"
        print(f"  {CYN}{tag}{RST} {BOLD}{label}{RST} {DIM}({rows:,} rows){RST}")

        config_results = {}
        for runner_idx, (lib_key, lib_label, runner) in enumerate(benchmark_runners):
            t0 = time.time()
            result, error_detail = runner(dataset, compression, rows)
            elapsed = time.time() - t0
            if result:
                results[(lib_key, key)] = result
                config_results[lib_key] = result
                mr = (rows / result["read_ms"]) * 1000 / 1e6
                sz = fmt_size(result["file_bytes"])
                print(f"       {lib_label:<9} W {fmt_ms(result['write_ms'])}  "
                      f"R {fmt_ms(result['read_ms'])} "
                      f"{DIM}({mr:.0f}M r/s, {sz}, {elapsed:.0f}s){RST}")
            else:
                print(f"       {lib_label:<9} {RED}failed{RST}")
                if error_detail:
                    print(f"       {DIM}{error_detail}{RST}")

            if runner_idx < len(benchmark_runners) - 1:
                cooldown(inter_library_cooldown)

        carquet_r = config_results.get("carquet")
        arrow_cpp_r = config_results.get("arrow_cpp")
        pyarrow_r = config_results.get("pyarrow")
        if carquet_r and arrow_cpp_r:
            ws = fmt_speedup(carquet_r["write_ms"], arrow_cpp_r["write_ms"])
            rs = fmt_speedup(carquet_r["read_ms"], arrow_cpp_r["read_ms"])
            print(f"       {DIM}vs Arrow C++  W {ws}  R {rs}{RST}")
        if carquet_r and pyarrow_r:
            ws = fmt_speedup(carquet_r["write_ms"], pyarrow_r["write_ms"])
            rs = fmt_speedup(carquet_r["read_ms"], pyarrow_r["read_ms"])
            print(f"       {DIM}vs PyArrow    W {ws}  R {rs}{RST}")

        # ── Cross-read: each library reads the other's file ──
        carquet_file = os.path.join(
            bench_tmpdir, f"benchmark_{dataset}_{compression}_carquet.parquet")
        arrow_file = os.path.join(
            bench_tmpdir, f"benchmark_{dataset}_{compression}_arrow_cpp.parquet")

        if do_crossread and carquet_r and arrow_cpp_r:
            print(f"       {DIM}cross-read...{RST}")

            # Arrow C++ reads Carquet-written file
            if os.path.isfile(carquet_file):
                cooldown(inter_library_cooldown)
                xr_arrow, xr_err = run_arrow_cpp_read_file(
                    arrow_cpp_bin, carquet_file, args.timeout)
                xr_carquet, xr_err2 = run_carquet_read_file(
                    carquet_bin, carquet_file, args.timeout)

                cr_entry = {}
                if crossread_result_is_valid(xr_arrow):
                    cr_entry["arrow_read_ms"] = xr_arrow["read_ms"]
                if crossread_result_is_valid(xr_carquet):
                    cr_entry["carquet_read_ms"] = xr_carquet["read_ms"]
                if cr_entry:
                    crossread_results[("carquet_file", key)] = cr_entry
                    c_ms = cr_entry.get("carquet_read_ms")
                    a_ms = cr_entry.get("arrow_read_ms")
                    print(f"       {DIM}Carquet file: "
                          f"Carquet R {fmt_ms(c_ms)}  Arrow R {fmt_ms(a_ms)}  "
                          f"{fmt_speedup(c_ms, a_ms)}{RST}")

            # Carquet reads Arrow-written file
            if os.path.isfile(arrow_file):
                cooldown(inter_library_cooldown)
                xr_carquet2, xr_err3 = run_carquet_read_file(
                    carquet_bin, arrow_file, args.timeout)
                xr_arrow2, xr_err4 = run_arrow_cpp_read_file(
                    arrow_cpp_bin, arrow_file, args.timeout)

                ar_entry = {}
                if crossread_result_is_valid(xr_carquet2):
                    ar_entry["carquet_read_ms"] = xr_carquet2["read_ms"]
                if crossread_result_is_valid(xr_arrow2):
                    ar_entry["arrow_read_ms"] = xr_arrow2["read_ms"]
                if ar_entry:
                    crossread_results[("arrow_file", key)] = ar_entry
                    c_ms = ar_entry.get("carquet_read_ms")
                    a_ms = ar_entry.get("arrow_read_ms")
                    print(f"       {DIM}Arrow file:   "
                          f"Carquet R {fmt_ms(c_ms)}  Arrow R {fmt_ms(a_ms)}  "
                          f"{fmt_speedup(c_ms, a_ms)}{RST}")

        # ── Clean up benchmark files ──
        for f in [carquet_file, arrow_file]:
            try:
                if os.path.isfile(f):
                    os.remove(f)
            except OSError:
                pass

        # Cooldown between configs
        if i < total - 1:
            cooldown(args.cooldown)
        print()

    # ── Final table ──
    print(f"  {_bar('═')}")
    print(f"  {BOLD}Results{RST}")
    print(f"  {_bar('═')}")
    print_final_tables(results, crossread_results)
    print()

    # ── Save JSON report ──
    report = {
        "carquet_version": carquet_version,
        "timestamp": timestamp,
        "system": sys_info,
        "config": {
            "warmup_iterations": 3,
            "bench_iterations": args.iterations if args.iterations > 0 else {"small": 51, "medium": 21, "large": 11, "xlarge": 11},
            "statistic": "trimmed_median",
            "cache_purge": True,
            "cooldown_seconds": args.cooldown,
            "subprocess_timeout_seconds": args.timeout,
        },
        "results": [],
    }
    if args.temp_dir:
        report["config"]["temp_dir"] = os.environ["CARQUET_BENCH_TMPDIR"]
    if arrow_cpp_version:
        report["arrow_cpp_version"] = arrow_cpp_version
    if pyarrow_version:
        report["pyarrow_version"] = pyarrow_version

    for (lib, (ds, comp)), r in sorted(results.items()):
        report["results"].append({
            "library": lib,
            "dataset": ds,
            "compression": comp,
            "rows": r["rows"],
            "write_ms": round(r["write_ms"], 2),
            "read_ms": round(r["read_ms"], 2),
            "file_bytes": r["file_bytes"],
        })

    if crossread_results:
        report["crossread_results"] = []
        for (writer, (ds, comp)), cr in sorted(crossread_results.items()):
            entry = {
                "writer": writer,
                "dataset": ds,
                "compression": comp,
            }
            if "carquet_read_ms" in cr:
                entry["carquet_read_ms"] = round(cr["carquet_read_ms"], 2)
            if "arrow_read_ms" in cr:
                entry["arrow_read_ms"] = round(cr["arrow_read_ms"], 2)
            report["crossread_results"].append(entry)

    date_str = datetime.date.today().strftime("%Y%m%d")
    os_tag = f"{platform.system().lower()}_{platform.machine()}"
    json_name = f"bench_{carquet_version}_{os_tag}_{date_str}.json"
    json_path = os.path.join(script_dir, json_name)
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2)

    print(f"  {DIM}Saved:{RST} {json_name}")
    print()


if __name__ == "__main__":
    main()
