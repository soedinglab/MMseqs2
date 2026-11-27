#!/usr/bin/env python3
"""
Carquet profiling orchestrator.

Cross-platform profiling with macOS (sample/Instruments) and Linux (perf)
support, flamegraph generation, and micro-benchmark analysis.

Usage:
    python3 profiling/run_profiler.py full                      # full profiling suite
    python3 profiling/run_profiler.py stat                      # quick timing stats
    python3 profiling/run_profiler.py stat --component rle      # specific component
    python3 profiling/run_profiler.py record                    # CPU profile recording
    python3 profiling/run_profiler.py flamegraph                # generate flamegraph SVG
    python3 profiling/run_profiler.py micro                     # micro-benchmarks
    python3 profiling/run_profiler.py compare                   # scalar vs SIMD
    python3 profiling/run_profiler.py build                     # just build
"""

import argparse
import datetime
import os
import platform
import re
import shutil
import subprocess
import sys
import time

# -- ANSI helpers -------------------------------------------------------------

BOLD = "\033[1m"
DIM = "\033[2m"
RST = "\033[0m"
RED = "\033[31m"
GRN = "\033[32m"
YLW = "\033[33m"
BLU = "\033[34m"
CYN = "\033[36m"
WHT = "\033[37m"

if not sys.stdout.isatty():
    BOLD = DIM = RST = RED = GRN = YLW = BLU = CYN = WHT = ""


def _bar(char="\u2500", width=62):
    return DIM + char * width + RST


def info(msg):
    print(f"  {BLU}[INFO]{RST} {msg}")


def ok(msg):
    print(f"  {GRN}[OK]{RST}   {msg}")


def warn(msg):
    print(f"  {YLW}[WARN]{RST} {msg}")


def error(msg):
    print(f"  {RED}[ERR]{RST}  {msg}")


# -- System info --------------------------------------------------------------

def get_cpu_name():
    if platform.system() == "Darwin":
        try:
            out = subprocess.check_output(
                ["sysctl", "-n", "machdep.cpu.brand_string"],
                stderr=subprocess.DEVNULL, text=True,
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
    if platform.system() == "Darwin":
        try:
            out = subprocess.check_output(
                ["sysctl", "-n", "hw.memsize"],
                stderr=subprocess.DEVNULL, text=True,
            ).strip()
            return int(out) / (1024 ** 3)
        except Exception:
            pass
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal"):
                    kb = int(line.split()[1])
                    return kb / (1024 ** 2)
    except Exception:
        pass
    return 0


def is_macos():
    return platform.system() == "Darwin"


def is_linux():
    return platform.system() == "Linux"


def core_count():
    try:
        return os.cpu_count() or 4
    except Exception:
        return 4


# -- Tool detection -----------------------------------------------------------

def has_tool(name):
    return shutil.which(name) is not None


def find_profiling_tools():
    """Detect available profiling tools on the current platform."""
    tools = {}
    if is_macos():
        if has_tool("sample"):
            tools["sample"] = shutil.which("sample")
        if has_tool("xctrace"):
            tools["xctrace"] = shutil.which("xctrace")
        if has_tool("instruments"):
            tools["instruments"] = shutil.which("instruments")
    if is_linux():
        if has_tool("perf"):
            tools["perf"] = shutil.which("perf")
    # FlameGraph tools (cross-platform)
    flamegraph_dir = os.environ.get("FLAMEGRAPH_DIR", os.path.expanduser("~/FlameGraph"))
    if os.path.isdir(flamegraph_dir):
        tools["flamegraph_dir"] = flamegraph_dir
    return tools


# -- Build --------------------------------------------------------------------

def build_binaries(project_dir, build_dir):
    """Build profiling binaries with debug symbols and optimization."""
    print()
    print(f"  {BOLD}Building profiling binaries{RST}")
    print(f"  {_bar()}")
    print(f"  {DIM}build:{RST}  {build_dir}")

    cmake_flags = "-g -fno-omit-frame-pointer"
    print(f"  {DIM}flags:{RST}  {cmake_flags}")
    print(f"  {_bar()}")
    print()

    # Configure
    cmake_args = [
        "cmake", "-B", build_dir,
        "-DCMAKE_BUILD_TYPE=RelWithDebInfo",
        f"-DCMAKE_C_FLAGS={cmake_flags}",
        "-DCARQUET_BUILD_BENCHMARKS=ON",
        "-DCARQUET_BUILD_TESTS=OFF",
        "-DCARQUET_BUILD_EXAMPLES=OFF",
        project_dir,
    ]

    print(f"  {DIM}Configuring...{RST}", flush=True)
    result = subprocess.run(cmake_args, capture_output=True, text=True)
    if result.returncode != 0:
        error("CMake configure failed")
        if result.stderr:
            for line in result.stderr.strip().splitlines()[-5:]:
                print(f"    {DIM}{line}{RST}")
        return False

    # Build
    jobs = core_count()
    print(f"  {DIM}Building ({jobs} jobs)...{RST}", flush=True)

    build_cmd = ["cmake", "--build", build_dir, f"-j{jobs}",
                 "--target", "profile_read", "--target", "profile_write",
                 "--target", "profile_micro"]
    result = subprocess.run(build_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        # Fallback: build all (some CMake versions don't support multiple --target)
        result = subprocess.run(
            ["cmake", "--build", build_dir, f"-j{jobs}"],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            error("Build failed")
            if result.stderr:
                for line in result.stderr.strip().splitlines()[-10:]:
                    print(f"    {DIM}{line}{RST}")
            return False

    # Verify
    built = []
    missing = []
    for name in ["profile_read", "profile_write", "profile_micro"]:
        path = os.path.join(build_dir, name)
        if os.path.isfile(path):
            built.append(name)
        else:
            missing.append(name)

    print()
    for name in built:
        print(f"  {GRN}built{RST}  {name}")
    for name in missing:
        print(f"  {RED}missing{RST}  {name}")
    print()

    if missing:
        error(f"{len(missing)} target(s) not built")
        return False

    ok("All profiling binaries built")
    return True


# -- Run profiling binary -----------------------------------------------------

def run_binary(binary, args, capture=False, timeout=600):
    """Run a profiling binary and optionally capture output."""
    argv = [binary] + args
    if capture:
        result = subprocess.run(
            argv, capture_output=True, text=True, timeout=timeout,
        )
        return result.stdout, result.stderr, result.returncode
    else:
        result = subprocess.run(argv, timeout=timeout)
        return None, None, result.returncode


def parse_profile_read_output(stdout):
    """Parse profile_read output for timing results."""
    results = {}
    for line in stdout.splitlines():
        line = line.strip()
        m = re.match(r"Avg read time:\s+([\d.]+)\s+ms", line)
        if m:
            results["avg_read_ms"] = float(m.group(1))
        m = re.match(r"Min read time:\s+([\d.]+)\s+ms", line)
        if m:
            results["min_read_ms"] = float(m.group(1))
        m = re.match(r"Max read time:\s+([\d.]+)\s+ms", line)
        if m:
            results["max_read_ms"] = float(m.group(1))
        m = re.match(r"Throughput:\s+([\d.]+)\s+M rows/sec", line)
        if m:
            results["throughput_mrows"] = float(m.group(1))
        m = re.match(r"Write time:\s+([\d.]+)\s+ms", line)
        if m:
            results["write_ms"] = float(m.group(1))
        m = re.match(r"File size:\s+([\d.]+)\s+MB", line)
        if m:
            results["file_size_mb"] = float(m.group(1))
    return results


def parse_micro_output(stdout):
    """Parse profile_micro output for benchmark results."""
    results = []
    for line in stdout.splitlines():
        m = re.match(r"\s+(.+?):\s+([\d.]+)\s+ns/value,\s+([\d.]+)\s+M values/sec", line)
        if m:
            results.append({
                "name": m.group(1).strip(),
                "ns_per_value": float(m.group(2)),
                "mvalues_per_sec": float(m.group(3)),
            })
            continue
        m = re.match(r"\s+(.+?):\s+([\d.]+)\s+MB/sec", line)
        if m:
            results.append({
                "name": m.group(1).strip(),
                "mb_per_sec": float(m.group(2)),
            })
            continue
        m = re.match(r"\s+(.+?):\s+([\d.]+)\s+ns/call", line)
        if m:
            results.append({
                "name": m.group(1).strip(),
                "ns_per_call": float(m.group(2)),
            })
    return results


# -- macOS: sample-based profiling --------------------------------------------

def macos_sample_profile(binary, args, duration, output_path):
    """Use macOS `sample` command for CPU profiling."""
    info(f"Starting profiling target...")

    # Launch the process
    proc = subprocess.Popen(
        [binary] + args,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    )

    pid = proc.pid
    info(f"Sampling PID {pid} for {duration}s...")

    sample_file = output_path + "_sample.txt"
    try:
        sample_result = subprocess.run(
            ["sample", str(pid), str(duration), "-f", sample_file],
            capture_output=True, text=True, timeout=duration + 30,
        )
        if sample_result.returncode != 0:
            warn(f"sample returned {sample_result.returncode}")
            if sample_result.stderr:
                for line in sample_result.stderr.strip().splitlines()[-3:]:
                    print(f"    {DIM}{line}{RST}")
    except subprocess.TimeoutExpired:
        warn("sample timed out")
    except Exception as e:
        warn(f"sample failed: {e}")

    # Wait for the process to finish
    try:
        proc.wait(timeout=duration + 60)
    except subprocess.TimeoutExpired:
        proc.terminate()
        proc.wait()

    if os.path.isfile(sample_file) and os.path.getsize(sample_file) > 0:
        ok(f"Sample data saved to {os.path.basename(sample_file)}")
        return sample_file
    else:
        warn("No sample data captured")
        return None


def macos_xctrace_profile(binary, args, output_path):
    """Use Instruments via xctrace for detailed profiling."""
    trace_file = output_path + ".trace"
    info(f"Recording with Instruments (Time Profiler)...")

    argv = [
        "xctrace", "record",
        "--template", "Time Profiler",
        "--output", trace_file,
        "--launch", "--",
    ] + [binary] + args

    result = subprocess.run(argv, capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        warn(f"xctrace returned {result.returncode}")
        if result.stderr:
            for line in result.stderr.strip().splitlines()[-3:]:
                print(f"    {DIM}{line}{RST}")
        return None

    if os.path.exists(trace_file):
        ok(f"Instruments trace saved to {os.path.basename(trace_file)}")
        info(f"Open with: open {trace_file}")
        return trace_file

    warn("No trace file generated")
    return None


def macos_sample_to_collapsed(sample_file, output_path):
    """Convert macOS sample output to collapsed stack format for flamegraph."""
    collapsed_file = output_path + "_collapsed.txt"

    try:
        with open(sample_file, "r") as f:
            content = f.read()
    except Exception as e:
        warn(f"Cannot read sample file: {e}")
        return None

    # Parse the sample output format into collapsed stacks
    # macOS sample output has indented call trees like:
    #   2701 Thread_8749954   DispatchQueue_1: ...
    #     2701 start  (in dyld) ...
    #       2701 main  (in profile_read) ...
    #         2688 read_batch_loop  (in profile_read) ...
    stacks = []
    current_stack = []
    current_count = 0

    for line in content.splitlines():
        stripped = line.rstrip()
        if not stripped:
            continue

        # Match sample lines: leading spaces + count + function_name
        m = re.match(r"^(\s+)(\d+)\s+(.+?)(?:\s+\(in\s+(.+?)\))?", stripped)
        if not m:
            # Flush if we had a stack
            if current_stack and current_count > 0:
                stacks.append((list(current_stack), current_count))
                current_stack = []
                current_count = 0
            continue

        indent = len(m.group(1))
        count = int(m.group(2))
        func = m.group(3).strip()
        # module = m.group(4) or ""

        # Clean up function name
        func = re.sub(r'\s+\+\s+\d+.*$', '', func)
        func = re.sub(r'\s+\[.*$', '', func)
        func = func.strip()

        if not func or func.startswith("?"):
            continue

        # Determine depth from indentation (2 spaces per level typically)
        depth = indent // 2

        # Trim stack to current depth
        while len(current_stack) > depth:
            if current_stack and current_count > 0:
                stacks.append((list(current_stack), current_count))
            current_stack.pop()

        if len(current_stack) < depth:
            # Padding if we jumped levels
            while len(current_stack) < depth:
                current_stack.append("???")

        if len(current_stack) == depth:
            current_stack.append(func)
        else:
            current_stack[depth] = func

        current_count = count

    # Flush last
    if current_stack and current_count > 0:
        stacks.append((list(current_stack), current_count))

    # Write collapsed format
    # A simpler approach: just extract leaf-most stacks
    collapsed = {}
    for stack, count in stacks:
        # Only output stacks at their deepest level
        key = ";".join(stack)
        if key in collapsed:
            collapsed[key] = max(collapsed[key], count)
        else:
            collapsed[key] = count

    if not collapsed:
        warn("Could not parse sample data into collapsed stacks")
        return None

    with open(collapsed_file, "w") as f:
        for stack_str, count in sorted(collapsed.items()):
            f.write(f"{stack_str} {count}\n")

    ok(f"Collapsed stacks: {len(collapsed)} unique ({os.path.basename(collapsed_file)})")
    return collapsed_file


# -- Linux: perf-based profiling ----------------------------------------------

def linux_perf_stat(binary, args, output_path):
    """Run perf stat for hardware counter statistics."""
    stat_file = output_path + "_stat.txt"
    info("Running perf stat...")

    events = "cycles,instructions,cache-references,cache-misses,branches,branch-misses"
    argv = ["perf", "stat", "-e", events] + [binary] + args

    result = subprocess.run(
        argv, capture_output=True, text=True, timeout=600,
    )

    output = result.stdout + "\n" + result.stderr
    with open(stat_file, "w") as f:
        f.write(output)

    # Print key stats
    for line in output.splitlines():
        line = line.strip()
        if any(k in line for k in ["cycles", "instructions", "cache", "branches", "insn per cycle"]):
            print(f"    {DIM}{line}{RST}")

    ok(f"Stats saved to {os.path.basename(stat_file)}")

    # Memory stats
    mem_events = "L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses"
    argv = ["perf", "stat", "-e", mem_events] + [binary] + args
    result = subprocess.run(argv, capture_output=True, text=True, timeout=600)
    mem_output = result.stdout + "\n" + result.stderr

    with open(stat_file, "a") as f:
        f.write("\n--- Memory Stats ---\n")
        f.write(mem_output)

    for line in mem_output.splitlines():
        line = line.strip()
        if any(k in line for k in ["L1-dcache", "LLC"]):
            print(f"    {DIM}{line}{RST}")

    return stat_file


def linux_perf_record(binary, args, output_path):
    """Record with perf for detailed profiling."""
    data_file = output_path + "_perf.data"
    info("Recording with perf (dwarf call graph)...")

    argv = [
        "perf", "record", "-g", "--call-graph", "dwarf", "-F", "999",
        "-o", data_file,
    ] + [binary] + args

    result = subprocess.run(argv, capture_output=True, text=True, timeout=600)
    if result.returncode != 0:
        warn(f"perf record returned {result.returncode}")

    if os.path.isfile(data_file) and os.path.getsize(data_file) > 0:
        ok(f"Recording saved to {os.path.basename(data_file)}")
        return data_file

    warn("No perf data captured")
    return None


def linux_perf_report(data_file, output_path):
    """Generate perf report from recorded data."""
    if not data_file or not os.path.isfile(data_file):
        error("No perf data found. Run 'record' first.")
        return None

    report_file = output_path + "_report.txt"
    top_file = output_path + "_top_functions.txt"

    info("Generating perf report...")

    # Hierarchical report
    subprocess.run(
        ["perf", "report", "-i", data_file, "--hierarchy", "--stdio",
         "--no-children", "--percent-limit", "0.5"],
        stdout=open(report_file, "w"), stderr=subprocess.DEVNULL,
    )

    # Top functions
    subprocess.run(
        ["perf", "report", "-i", data_file, "--stdio",
         "--no-children", "--percent-limit", "0.5", "-n"],
        stdout=open(top_file, "w"), stderr=subprocess.DEVNULL,
    )

    ok(f"Reports saved")

    # Print top functions
    print()
    print(f"  {BOLD}Top Functions{RST}")
    print(f"  {_bar()}")
    try:
        with open(top_file, "r") as f:
            count = 0
            for line in f:
                if re.match(r"\s+[\d.]+%", line) and count < 15:
                    print(f"  {line.rstrip()}")
                    count += 1
    except Exception:
        pass

    return report_file


def linux_perf_annotate(data_file, output_path):
    """Generate source annotations for hot functions."""
    if not data_file or not os.path.isfile(data_file):
        return

    functions = [
        "carquet_rle_decoder_get",
        "carquet_rle_decoder_get_batch",
        "carquet_rle_decode_levels",
        "read_batch_loop",
        "carquet_dispatch_gather_i32",
        "carquet_lz4_decompress",
    ]

    ann_dir = output_path + "_annotations"
    os.makedirs(ann_dir, exist_ok=True)

    info("Generating source annotations...")
    annotated = 0
    for func in functions:
        out_file = os.path.join(ann_dir, f"{func}.txt")
        result = subprocess.run(
            ["perf", "annotate", "-i", data_file, "-s", func, "--stdio"],
            stdout=open(out_file, "w"), stderr=subprocess.DEVNULL,
        )
        if os.path.isfile(out_file) and os.path.getsize(out_file) > 100:
            annotated += 1
        else:
            try:
                os.remove(out_file)
            except OSError:
                pass

    if annotated > 0:
        ok(f"Annotated {annotated} function(s) in {os.path.basename(ann_dir)}/")


def linux_perf_to_collapsed(data_file, output_path):
    """Convert perf data to collapsed stack format."""
    collapsed_file = output_path + "_collapsed.txt"

    # Try perf script -> stackcollapse-perf.pl
    # First check if we have FlameGraph tools
    flamegraph_dir = os.environ.get("FLAMEGRAPH_DIR", os.path.expanduser("~/FlameGraph"))
    collapse_script = os.path.join(flamegraph_dir, "stackcollapse-perf.pl")

    if os.path.isfile(collapse_script):
        perf_script = subprocess.Popen(
            ["perf", "script", "-i", data_file],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
        )
        with open(collapsed_file, "w") as f:
            collapse = subprocess.run(
                [collapse_script, "--all"],
                stdin=perf_script.stdout, stdout=f, stderr=subprocess.DEVNULL,
            )
        perf_script.wait()
        if collapse.returncode == 0 and os.path.getsize(collapsed_file) > 0:
            return collapsed_file

    # Fallback: try inferno or flamegraph-rs if installed
    if has_tool("inferno-collapse-perf"):
        perf_script = subprocess.Popen(
            ["perf", "script", "-i", data_file],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
        )
        with open(collapsed_file, "w") as f:
            subprocess.run(
                ["inferno-collapse-perf"],
                stdin=perf_script.stdout, stdout=f, stderr=subprocess.DEVNULL,
            )
        perf_script.wait()
        if os.path.getsize(collapsed_file) > 0:
            return collapsed_file

    warn("No stack collapse tool found (FlameGraph or inferno)")
    return None


# -- Flamegraph generation (cross-platform) -----------------------------------

def generate_flamegraph(collapsed_file, output_path, title, tools):
    """Generate flamegraph SVG from collapsed stack data."""
    if not collapsed_file or not os.path.isfile(collapsed_file):
        error("No collapsed stack data available")
        return None

    flamegraph_svg = output_path + "_flamegraph.svg"
    icicle_svg = output_path + "_icicle.svg"

    flamegraph_dir = tools.get("flamegraph_dir")
    fg_script = os.path.join(flamegraph_dir, "flamegraph.pl") if flamegraph_dir else None

    if fg_script and os.path.isfile(fg_script):
        info("Generating flamegraph with FlameGraph toolkit...")

        # Normal flamegraph
        with open(collapsed_file, "r") as fin, open(flamegraph_svg, "w") as fout:
            subprocess.run(
                [fg_script, "--title", title, "--width", "1800", "--colors", "hot"],
                stdin=fin, stdout=fout, stderr=subprocess.DEVNULL,
            )

        if os.path.isfile(flamegraph_svg) and os.path.getsize(flamegraph_svg) > 0:
            ok(f"Flamegraph: {os.path.basename(flamegraph_svg)}")
        else:
            warn("Flamegraph generation failed")
            flamegraph_svg = None

        # Icicle graph (reversed)
        with open(collapsed_file, "r") as fin, open(icicle_svg, "w") as fout:
            subprocess.run(
                [fg_script, "--title", f"{title} (Reversed)",
                 "--reverse", "--inverted", "--width", "1800", "--colors", "hot"],
                stdin=fin, stdout=fout, stderr=subprocess.DEVNULL,
            )

        if os.path.isfile(icicle_svg) and os.path.getsize(icicle_svg) > 0:
            ok(f"Icicle graph: {os.path.basename(icicle_svg)}")

        return flamegraph_svg

    # Fallback: try inferno
    if has_tool("inferno-flamegraph"):
        info("Generating flamegraph with inferno...")
        with open(collapsed_file, "r") as fin, open(flamegraph_svg, "w") as fout:
            subprocess.run(
                ["inferno-flamegraph", "--title", title],
                stdin=fin, stdout=fout, stderr=subprocess.DEVNULL,
            )
        if os.path.isfile(flamegraph_svg) and os.path.getsize(flamegraph_svg) > 0:
            ok(f"Flamegraph: {os.path.basename(flamegraph_svg)}")
            return flamegraph_svg

    error("No flamegraph tool found")
    print(f"    {DIM}Install FlameGraph: git clone https://github.com/brendangregg/FlameGraph ~/FlameGraph{RST}")
    print(f"    {DIM}Or install inferno: cargo install inferno{RST}")
    return None


# -- Stat mode (timing only, cross-platform) ----------------------------------

def run_stat(build_dir, config, tools):
    """Run timing statistics by executing the profiling binaries."""
    binary = os.path.join(build_dir, "profile_read")
    micro_binary = os.path.join(build_dir, "profile_micro")

    if config.component != "all" and config.component != "read":
        # Run micro-benchmark for specific component
        if not os.path.isfile(micro_binary):
            error(f"profile_micro not found at {micro_binary}")
            return

        args = ["--component", config.component,
                "--count", "1000000", "--iterations", "100"]

        if is_linux() and "perf" in tools:
            linux_perf_stat(micro_binary, args, config.output_prefix)
        else:
            info(f"Running micro-benchmark: {config.component}")
            stdout, _, rc = run_binary(micro_binary, args, capture=True)
            if rc != 0:
                error("Micro-benchmark failed")
                return
            results = parse_micro_output(stdout)
            print()
            print_micro_results(results)
        return

    # Full read path
    if not os.path.isfile(binary):
        error(f"profile_read not found at {binary}")
        return

    args = ["-r", str(config.rows), "-i", str(config.iterations), "-v"]
    if config.dictionary:
        args.append("-d")
    if config.compression > 0:
        args.extend(["-c", str(config.compression)])
    if config.nulls > 0:
        args.extend(["-n", str(config.nulls)])

    if is_linux() and "perf" in tools:
        linux_perf_stat(binary, args, config.output_prefix)
    else:
        info(f"Running profile_read ({config.rows:,} rows, {config.iterations} iterations)")
        stdout, _, rc = run_binary(binary, args, capture=True)
        if rc != 0:
            error("profile_read failed")
            return
        results = parse_profile_read_output(stdout)
        print()
        print_read_results(results, config)


def print_read_results(results, config):
    """Pretty-print profile_read results."""
    print(f"  {BOLD}Read Path Results{RST}")
    print(f"  {_bar()}")

    if "write_ms" in results:
        print(f"  {DIM}write:{RST}      {results['write_ms']:>8.2f} ms")
    if "file_size_mb" in results:
        print(f"  {DIM}file size:{RST}  {results['file_size_mb']:>8.2f} MB")

    if "avg_read_ms" in results:
        avg = results["avg_read_ms"]
        low = results.get("min_read_ms", avg)
        high = results.get("max_read_ms", avg)
        tp = results.get("throughput_mrows", 0)

        print(f"  {DIM}avg read:{RST}   {avg:>8.2f} ms")
        print(f"  {DIM}min read:{RST}   {low:>8.2f} ms")
        print(f"  {DIM}max read:{RST}   {high:>8.2f} ms")

        if tp > 0:
            color = GRN if tp > 75 else (YLW if tp > 30 else RED)
            print(f"  {DIM}throughput:{RST} {color}{BOLD}{tp:>7.2f} M rows/sec{RST}")

        # Variance indicator
        if high > 0 and low > 0:
            variance_pct = ((high - low) / avg) * 100
            if variance_pct > 20:
                warn(f"High variance: {variance_pct:.0f}% (consider more iterations or cooldown)")

    print(f"  {_bar()}")


def print_micro_results(results):
    """Pretty-print micro-benchmark results."""
    print(f"  {BOLD}Micro-benchmark Results{RST}")
    print(f"  {_bar()}")

    for r in results:
        name = r["name"]
        if "ns_per_value" in r:
            ns = r["ns_per_value"]
            mv = r["mvalues_per_sec"]
            color = GRN if mv > 500 else (YLW if mv > 100 else WHT)
            print(f"  {name:<35} {ns:>7.2f} ns/val  {color}{mv:>7.2f} M val/s{RST}")
        elif "mb_per_sec" in r:
            mbs = r["mb_per_sec"]
            color = GRN if mbs > 1000 else (YLW if mbs > 200 else WHT)
            print(f"  {name:<35} {color}{mbs:>7.2f} MB/s{RST}")
        elif "ns_per_call" in r:
            ns = r["ns_per_call"]
            print(f"  {name:<35} {ns:>7.2f} ns/call")

    print(f"  {_bar()}")


# -- Record mode --------------------------------------------------------------

def run_record(build_dir, config, tools):
    """Record CPU profile data."""
    binary = os.path.join(build_dir, "profile_read")

    if config.component != "all" and config.component != "read":
        binary = os.path.join(build_dir, "profile_micro")
        args = ["--component", config.component,
                "--count", "2000000", "--iterations", "200"]
    else:
        args = ["-r", str(config.rows), "-i", str(config.iterations), "-v"]
        if config.dictionary:
            args.append("-d")
        if config.compression > 0:
            args.extend(["-c", str(config.compression)])
        if config.nulls > 0:
            args.extend(["-n", str(config.nulls)])

    if not os.path.isfile(binary):
        error(f"{os.path.basename(binary)} not found")
        return None

    if is_linux() and "perf" in tools:
        return linux_perf_record(binary, args, config.output_prefix)

    if is_macos():
        # Prefer xctrace if available, fall back to sample
        if "xctrace" in tools:
            return macos_xctrace_profile(binary, args, config.output_prefix)
        elif "sample" in tools:
            # Estimate duration from rows/iterations
            duration = max(10, min(60, config.rows // 1000000 * config.iterations))
            return macos_sample_profile(binary, args, duration, config.output_prefix)
        else:
            error("No profiling tool available (need sample or xctrace)")
            return None

    error("No supported profiler for this platform")
    return None


# -- Report mode --------------------------------------------------------------

def run_report(build_dir, config, tools):
    """Generate and display profiling report."""
    if is_linux() and "perf" in tools:
        data_file = config.output_prefix + "_perf.data"
        if not os.path.isfile(data_file):
            error("No perf data found. Run 'record' first.")
            return
        linux_perf_report(data_file, config.output_prefix)
        linux_perf_annotate(data_file, config.output_prefix)
        return

    if is_macos():
        sample_file = config.output_prefix + "_sample.txt"
        trace_file = config.output_prefix + ".trace"

        if os.path.isfile(sample_file):
            info(f"Sample data available: {os.path.basename(sample_file)}")
            # Parse and display top functions from sample
            print_sample_top_functions(sample_file)
        elif os.path.exists(trace_file):
            info(f"Instruments trace available: {os.path.basename(trace_file)}")
            info(f"Open with: open {trace_file}")
        else:
            warn("No profiling data found. Run 'record' first.")


def print_sample_top_functions(sample_file):
    """Extract and print top functions from macOS sample output."""
    try:
        with open(sample_file, "r") as f:
            content = f.read()
    except Exception:
        return

    # Collect function sample counts
    func_counts = {}
    for line in content.splitlines():
        m = re.match(r"\s+(\d+)\s+(.+?)(?:\s+\(in\s+(.+?)\))?", line)
        if m:
            count = int(m.group(1))
            func = m.group(2).strip()
            module = m.group(3) or ""
            func = re.sub(r'\s+\+\s+\d+.*$', '', func).strip()
            if func and not func.startswith("?") and "Thread_" not in func:
                key = f"{func} ({module})" if module else func
                func_counts[key] = max(func_counts.get(key, 0), count)

    if not func_counts:
        return

    # Sort by count, show top 20
    sorted_funcs = sorted(func_counts.items(), key=lambda x: -x[1])
    total = sorted_funcs[0][1] if sorted_funcs else 1

    print()
    print(f"  {BOLD}Top Functions (from sample){RST}")
    print(f"  {_bar()}")

    for func, count in sorted_funcs[:20]:
        pct = (count / total) * 100
        bar_w = 20
        filled = int(pct / 100 * bar_w)
        bar = GRN + "\u2588" * filled + DIM + "\u2591" * (bar_w - filled) + RST
        print(f"  {bar} {pct:>5.1f}%  {func}")

    print(f"  {_bar()}")


# -- Flamegraph mode ----------------------------------------------------------

def run_flamegraph(build_dir, config, tools):
    """Generate flamegraph from profiling data."""
    collapsed_file = None

    if is_linux() and "perf" in tools:
        data_file = config.output_prefix + "_perf.data"
        if not os.path.isfile(data_file):
            info("No perf data found, recording first...")
            data_file = run_record(build_dir, config, tools)
        if data_file:
            collapsed_file = linux_perf_to_collapsed(data_file, config.output_prefix)

    elif is_macos():
        sample_file = config.output_prefix + "_sample.txt"
        if not os.path.isfile(sample_file):
            info("No sample data found, recording first...")
            binary = os.path.join(build_dir, "profile_read")
            args = ["-r", str(config.rows), "-i", str(config.iterations)]
            if config.dictionary:
                args.append("-d")
            if config.compression > 0:
                args.extend(["-c", str(config.compression)])
            duration = max(10, min(60, config.rows // 1000000 * config.iterations))
            sample_file = macos_sample_profile(binary, args, duration, config.output_prefix)

        if sample_file:
            collapsed_file = macos_sample_to_collapsed(sample_file, config.output_prefix)

    title = f"Carquet Read Path ({config.rows:,} rows)"
    generate_flamegraph(collapsed_file, config.output_prefix, title, tools)


# -- Micro mode ---------------------------------------------------------------

def run_micro(build_dir, config, tools):
    """Run micro-benchmarks."""
    binary = os.path.join(build_dir, "profile_micro")
    if not os.path.isfile(binary):
        error(f"profile_micro not found at {binary}")
        return

    args = ["--component", config.component,
            "--count", "2000000", "--iterations", "200"]

    info(f"Running micro-benchmarks (component: {config.component})")
    stdout, _, rc = run_binary(binary, args, capture=True)
    if rc != 0:
        error("Micro-benchmark failed")
        return

    results = parse_micro_output(stdout)
    print()
    print_micro_results(results)

    # Save raw output
    out_file = config.output_prefix + "_micro.txt"
    with open(out_file, "w") as f:
        f.write(stdout)
    ok(f"Raw output saved to {os.path.basename(out_file)}")


# -- Compare mode -------------------------------------------------------------

def run_compare(build_dir, config, tools):
    """Compare scalar vs SIMD implementations."""
    binary = os.path.join(build_dir, "profile_micro")
    if not os.path.isfile(binary):
        error(f"profile_micro not found at {binary}")
        return

    components = ["gather", "null", "dispatch"]
    all_results = {}

    for comp in components:
        args = ["--component", comp, "--count", "2000000", "--iterations", "200"]
        if comp == "dispatch":
            args = ["--component", comp, "--iterations", "500000"]

        tag = f"[{components.index(comp)+1}/{len(components)}]"
        print(f"  {CYN}{tag}{RST} {BOLD}{comp}{RST}")

        stdout, _, rc = run_binary(binary, args, capture=True)
        if rc != 0:
            error(f"{comp} benchmark failed")
            continue

        results = parse_micro_output(stdout)
        all_results[comp] = results

        for r in results:
            name = r["name"]
            if "ns_per_value" in r:
                print(f"       {name:<30} {r['ns_per_value']:>7.2f} ns/val  {r['mvalues_per_sec']:>7.2f} M val/s")
            elif "mb_per_sec" in r:
                print(f"       {name:<30} {r['mb_per_sec']:>7.2f} MB/s")
            elif "ns_per_call" in r:
                print(f"       {name:<30} {r['ns_per_call']:>7.2f} ns/call")
        print()

    # Print speedup summary
    print(f"  {BOLD}Speedup Summary (dispatch/SIMD vs scalar){RST}")
    print(f"  {_bar()}")

    for comp, results in all_results.items():
        scalar_val = None
        dispatch_val = None
        for r in results:
            name = r["name"].lower()
            if "scalar" in name:
                scalar_val = r.get("ns_per_value") or r.get("ns_per_call")
            elif "dispatch" in name:
                dispatch_val = r.get("ns_per_value") or r.get("ns_per_call")

        if scalar_val and dispatch_val and dispatch_val > 0:
            speedup = scalar_val / dispatch_val
            color = GRN if speedup > 1.5 else (YLW if speedup > 1.0 else RED)
            print(f"  {comp:<20} {color}{BOLD}{speedup:>5.2f}x{RST}")

    print(f"  {_bar()}")

    # Save
    out_file = config.output_prefix + "_compare.txt"
    with open(out_file, "w") as f:
        for comp, results in all_results.items():
            f.write(f"=== {comp} ===\n")
            for r in results:
                f.write(f"  {r}\n")
            f.write("\n")
    ok(f"Comparison saved to {os.path.basename(out_file)}")


# -- Write mode ---------------------------------------------------------------

def run_write(build_dir, config, args):
    """Run write path profiling."""
    binary = os.path.join(build_dir, "profile_write")
    if not os.path.isfile(binary):
        error(f"profile_write not found at {binary}")
        return

    bin_args = ["-r", str(config.rows), "-i", str(config.iterations), "-v"]
    if args.dictionary:
        bin_args.append("-d")
    if args.compression > 0:
        bin_args.extend(["-c", str(args.compression)])
    if args.nulls > 0:
        bin_args.extend(["-n", str(args.nulls)])
    if args.columns != 5:
        bin_args.extend(["-C", str(args.columns)])
    if args.bloom:
        bin_args.append("-B")
    if args.page_index:
        bin_args.append("-P")
    if args.overhead:
        bin_args.append("--overhead")
    if args.codecs:
        bin_args.append("--codecs")
    if args.scaling:
        bin_args.append("--scaling")
    if args.encoding_cmp:
        bin_args.append("--encoding")

    info(f"Running write profiler ({config.rows:,} rows, {args.columns} cols)")
    stdout, _, rc = run_binary(binary, bin_args, capture=True, timeout=1200)
    if rc != 0:
        error("profile_write failed")
        return

    # Print the raw output (it's already well-formatted)
    print()
    if stdout:
        for line in stdout.splitlines():
            print(f"  {line}")

    # Save
    out_file = config.output_prefix + "_write.txt"
    with open(out_file, "w") as f:
        f.write(stdout or "")
    ok(f"Output saved to {os.path.basename(out_file)}")


# -- Full mode ----------------------------------------------------------------

def run_full(build_dir, config, tools):
    """Run full profiling suite."""
    steps = ["stat", "record", "report", "flamegraph"]
    total = len(steps)

    for i, step in enumerate(steps):
        tag = f"[{i+1}/{total}]"
        print()
        print(f"  {CYN}{tag}{RST} {BOLD}{step}{RST}")
        print(f"  {_bar()}")

        if step == "stat":
            run_stat(build_dir, config, tools)
        elif step == "record":
            run_record(build_dir, config, tools)
        elif step == "report":
            run_report(build_dir, config, tools)
        elif step == "flamegraph":
            run_flamegraph(build_dir, config, tools)

    # Summary of generated files
    print()
    dbar = _bar("\u2550")
    print(f"  {dbar}")
    print(f"  {BOLD}Output Files{RST}")
    print(f"  {dbar}")

    output_dir = os.path.dirname(config.output_prefix)
    prefix_base = os.path.basename(config.output_prefix)
    if os.path.isdir(output_dir):
        for entry in sorted(os.scandir(output_dir), key=lambda e: e.name):
            if entry.name.startswith(prefix_base) and entry.is_file():
                size_kb = entry.stat().st_size / 1024
                print(f"  {DIM}{size_kb:>8.1f} KB{RST}  {entry.name}")

    print()


# -- Config -------------------------------------------------------------------

class ProfileConfig:
    def __init__(self, args, output_dir):
        self.rows = args.rows
        self.iterations = args.iterations
        self.component = args.component
        self.dictionary = args.dictionary
        self.compression = args.compression
        self.nulls = args.nulls
        self.output_prefix = os.path.join(
            output_dir,
            f"carquet_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}",
        )


# -- Main ---------------------------------------------------------------------

COMPRESSION_NAMES = {0: "none", 1: "snappy", 2: "zstd", 3: "lz4"}
NULL_NAMES = {0: "none", 1: "10%", 2: "30%", 3: "50%"}


def main():
    parser = argparse.ArgumentParser(
        description="Carquet profiling orchestrator (macOS + Linux)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""modes:
  full        Full profiling suite (stat + record + report + flamegraph)
  stat        Quick timing statistics
  record      Record CPU profile (perf on Linux, sample/Instruments on macOS)
  report      Display profiling report from recorded data
  flamegraph  Generate flamegraph SVG
  micro       Run isolated micro-benchmarks
  compare     Compare scalar vs SIMD implementations
  write       Write path profiling (encoding, compression, features)
  build       Just build profiling binaries
""",
    )
    parser.add_argument("mode", nargs="?", default="full",
                        choices=["full", "stat", "record", "report",
                                 "flamegraph", "micro", "compare",
                                 "write", "build"],
                        help="Profiling mode (default: full)")
    parser.add_argument("--rows", type=int, default=5_000_000,
                        help="Number of rows (default: 5,000,000)")
    parser.add_argument("--iterations", type=int, default=5,
                        help="Number of iterations (default: 5)")
    parser.add_argument("--component", default="all",
                        choices=["all", "read", "rle", "gather", "null",
                                 "compression", "dispatch", "minmax", "bss",
                                 "prefix_sum", "hash", "bloom", "encoding"],
                        help="Component to profile (default: all)")
    parser.add_argument("--dictionary", action="store_true",
                        help="Enable dictionary encoding")
    parser.add_argument("--compression", type=int, default=0,
                        choices=[0, 1, 2, 3],
                        help="Compression: 0=none, 1=snappy, 2=zstd, 3=lz4")
    parser.add_argument("--nulls", type=int, default=0,
                        choices=[0, 1, 2, 3],
                        help="Null ratio: 0=none, 1=10%%, 2=30%%, 3=50%%")
    parser.add_argument("--columns", type=int, default=5,
                        help="Number of columns for write profiling (default: 5)")
    parser.add_argument("--overhead", action="store_true",
                        help="Measure feature overhead (write mode)")
    parser.add_argument("--codecs", action="store_true",
                        help="Compare compression codecs (write mode)")
    parser.add_argument("--scaling", action="store_true",
                        help="Measure column count scaling (write mode)")
    parser.add_argument("--encoding-cmp", action="store_true",
                        help="Compare plain vs dictionary encoding (write mode)")
    parser.add_argument("--bloom", action="store_true",
                        help="Enable bloom filters (write mode)")
    parser.add_argument("--page-index", action="store_true",
                        help="Enable page index (write mode)")
    parser.add_argument("--output", help="Output directory (default: profiling/output)")
    parser.add_argument("--build-dir", help="Build directory (default: <project>/build)")

    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    build_dir = args.build_dir or os.path.join(project_dir, "build")
    output_dir = args.output or os.path.join(script_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    config = ProfileConfig(args, output_dir)

    # Detect tools
    tools = find_profiling_tools()

    # -- Header --
    print()
    print(f"  {BOLD}Carquet Profiling Suite{RST}")
    print(f"  {_bar()}")
    print(f"  {DIM}CPU:{RST}         {get_cpu_name()}")
    mem = get_mem_gb()
    if mem > 0:
        print(f"  {DIM}RAM:{RST}         {mem:.0f} GB")
    print(f"  {DIM}OS:{RST}          {platform.platform()}")
    print(f"  {DIM}mode:{RST}        {args.mode}")
    print(f"  {DIM}rows:{RST}        {args.rows:,}")
    print(f"  {DIM}iterations:{RST}  {args.iterations}")
    print(f"  {DIM}component:{RST}   {args.component}")
    print(f"  {DIM}compression:{RST} {COMPRESSION_NAMES[args.compression]}")
    print(f"  {DIM}nulls:{RST}       {NULL_NAMES[args.nulls]}")
    print(f"  {DIM}dictionary:{RST}  {'yes' if args.dictionary else 'no'}")

    # Tool availability
    tool_names = []
    if "perf" in tools:
        tool_names.append("perf")
    if "sample" in tools:
        tool_names.append("sample")
    if "xctrace" in tools:
        tool_names.append("xctrace")
    if "flamegraph_dir" in tools:
        tool_names.append("FlameGraph")
    if has_tool("inferno-flamegraph"):
        tool_names.append("inferno")
    print(f"  {DIM}tools:{RST}       {', '.join(tool_names) if tool_names else 'none detected'}")

    if not tool_names and args.mode not in ("build", "stat", "micro", "compare"):
        warn("No platform profiling tools detected")
        if is_macos():
            print(f"    {DIM}sample and xctrace should be available with Xcode Command Line Tools{RST}")
        elif is_linux():
            print(f"    {DIM}Install perf: sudo apt install linux-tools-generic{RST}")

    print(f"  {_bar()}")

    # -- Build --
    if args.mode == "build":
        build_binaries(project_dir, build_dir)
        return

    # Auto-build if needed
    profile_read = os.path.join(build_dir, "profile_read")
    profile_write = os.path.join(build_dir, "profile_write")
    profile_micro = os.path.join(build_dir, "profile_micro")
    if (not os.path.isfile(profile_read) or not os.path.isfile(profile_write)
            or not os.path.isfile(profile_micro)):
        if not build_binaries(project_dir, build_dir):
            sys.exit(1)

    # -- Dispatch --
    if args.mode == "full":
        run_full(build_dir, config, tools)
    elif args.mode == "stat":
        run_stat(build_dir, config, tools)
    elif args.mode == "record":
        run_record(build_dir, config, tools)
    elif args.mode == "report":
        run_report(build_dir, config, tools)
    elif args.mode == "flamegraph":
        run_flamegraph(build_dir, config, tools)
    elif args.mode == "micro":
        run_micro(build_dir, config, tools)
    elif args.mode == "compare":
        run_compare(build_dir, config, tools)
    elif args.mode == "write":
        run_write(build_dir, config, args)

    print()


if __name__ == "__main__":
    main()
