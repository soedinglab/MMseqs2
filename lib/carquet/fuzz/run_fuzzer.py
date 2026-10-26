#!/usr/bin/env python3
"""
Carquet fuzzing orchestrator.

Builds and runs libFuzzer / AFL++ fuzz targets with sanitizers,
manages corpus directories, and reports results.

Usage:
    python3 fuzz/run_fuzzer.py reader                  # fuzz reader for 5 min
    python3 fuzz/run_fuzzer.py all --time 3600          # all targets, 1 hour each
    python3 fuzz/run_fuzzer.py compression --jobs 4     # 4 parallel jobs
    python3 fuzz/run_fuzzer.py build                    # just build
    python3 fuzz/run_fuzzer.py list                     # list available targets
    python3 fuzz/run_fuzzer.py minimize reader          # minimize reader corpus
    python3 fuzz/run_fuzzer.py coverage reader          # coverage report
"""

import argparse
import datetime
import glob
import os
import platform
import re
import shutil
import subprocess
import sys
import threading
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

if not sys.stdout.isatty():
    BOLD = DIM = RST = RED = GRN = YLW = BLU = CYN = WHT = ""


def _bar(char="\u2500", width=62):
    return DIM + char * width + RST


# ── Target definitions ───────────────────────────────────────────────────────

TARGETS = [
    ("reader",       "Full Parquet file reader (file format, metadata, all decoders)"),
    ("writer",       "Schema creation, data writing, write-then-read roundtrip"),
    ("compression",  "Compression codecs: Snappy, LZ4, GZIP, ZSTD (encode + decode)"),
    ("encodings",    "Encoding decoders: RLE, Delta, Plain, Dictionary, BSS, DeltaLen"),
    ("thrift",       "Thrift compact protocol: primitives, structs, file/page metadata"),
    ("roundtrip",    "Encode-decode consistency for delta, LZ4, BSS, RLE, Snappy, dict"),
    ("page_filter",  "Page-level predicate filter over column-index/offset-index stats"),
    ("append",       "open_append footer parse, schema-match, restore + footer rewrite"),
    ("arrow_schema", "ARROW:schema FlatBuffer reader recovering per-field custom_metadata"),
    ("arrow_c_data", "Arrow C Data Interface format-string parser (import_schema)"),
    ("nested_write", "Nested write helper: auto-shredding LIST/MAP offsets+validity"),
    ("arrow_read",   "Nested Arrow reassembler (carquet_reader_read_arrow) over file bytes"),
]
TARGET_NAMES = [t[0] for t in TARGETS]

# Targets that consume whole Parquet files and therefore benefit from the
# real-world apache/parquet-testing seeds fetched by fuzz/fetch_corpus.sh.
EXTERNAL_SEED_TARGETS = {"reader", "page_filter", "arrow_read"}


# ── Seed corpus ──────────────────────────────────────────────────────────────

SEEDS = {
    "reader": {
        "seed_magic": b"PAR1",
        "seed_minimal": b"PAR1\x00\x00\x00\x00PAR1",
        "seed_empty_footer": (
            b"PAR1"
            b"\x15\x00\x15\x00\x15\x00\x00"  # minimal thrift
            b"\x07\x00\x00\x00"               # footer length
            b"PAR1"
        ),
    },
    "writer": {
        "seed_simple": b"\x01\x0a\x00\x00\x00\x01",
        "seed_multi": b"\x03\x05\x00\x04\x55\x01\x02\x03",
        "seed_snappy": b"\x02\x64\x00\x01\x00\x01\x04",
        "seed_nullable": b"\x02\x08\x00\x00\xff\x01\x03",
        "seed_all_types": b"\x06\x04\x00\x00\x00\x00\x01\x02\x03\x04\x05",
    },
    "compression": {
        "seed_snappy": b"\x00\x00",
        "seed_lz4": b"\x01\x00",
        "seed_gzip": b"\x02\x00",
        "seed_zstd": b"\x03\x00",
        "seed_snappy_rt": b"\x04hello world test data",
        "seed_lz4_rt": b"\x05the quick brown fox jumps",
        "seed_gzip_rt": b"\x06repeated repeated repeated",
        "seed_zstd_rt": b"\x07\x00\x01\x02\x03\x04\x05\x06\x07",
    },
    "encodings": {
        "seed_rle": b"\x00\x08\x00",
        "seed_delta32": b"\x01\x10\x00",
        "seed_delta64": b"\x02\x10\x00",
        "seed_plain32": b"\x03\x10\x00\x00\x00\x00",
        "seed_plain64": b"\x04\x10\x00\x00\x00\x00\x00\x00\x00\x00",
        "seed_plain_double": b"\x05\x10\x00\x00\x00\x00\x00\x00\x00\x00",
        "seed_dict32": b"\x06\x10\x01\x00\x00\x00\x02\x00",
        "seed_dict64": b"\x07\x10\x01\x00\x00\x00\x00\x00\x00\x00\x02\x00",
        "seed_dict_float": b"\x08\x10\x00\x00\x80\x3f\x02\x00",
        "seed_dict_double": b"\x09\x10\x00\x00\x00\x00\x00\x00\xf0\x3f\x02\x00",
        "seed_bss_float": b"\x0a\x10\x00\x00\x00\x00\x00\x00\x00\x00",
        "seed_bss_double": (
            b"\x0b\x10\x00\x00\x00\x00\x00\x00\x00\x00"
            b"\x00\x00\x00\x00\x00\x00\x00\x00"
        ),
        "seed_plain_bool": b"\x0c\x08\xff",
        "seed_plain_float": b"\x0d\x04\x00\x00\x80\x3f",
        "seed_delta_len": b"\x0e\x04\x03\x00\x00\x00abc",
        "seed_rle_levels": b"\x0f\x08\x00",
        "seed_dict_bytearray": b"\x10\x02\x03\x00\x00\x00foo\x02\x00",
    },
    "thrift": {
        "seed_empty": b"\x00\x00",
        "seed_struct": b"\x15\x00\x00",
        "seed_list": b"\x19\x00\x00\x00",
        "seed_nested": b"\x15\x15\x00\x00\x00",
        "seed_map": b"\x1b\x55\x00\x00",
    },
    "roundtrip": {
        "seed_delta32": b"\x00\x01\x00\x00\x00\x02\x00\x00\x00\x03\x00\x00\x00",
        "seed_delta64": (
            b"\x01\x01\x00\x00\x00\x00\x00\x00\x00"
            b"\x02\x00\x00\x00\x00\x00\x00\x00"
        ),
        "seed_lz4": b"\x02hello world test data",
        "seed_bss_float": b"\x03\x00\x00\x80\x3f\x00\x00\x00\x40\x00\x00\x40\x40",
        "seed_bss_double": (
            b"\x04\x00\x00\x00\x00\x00\x00\xf0\x3f"
            b"\x00\x00\x00\x00\x00\x00\x00\x40"
        ),
        "seed_rle": b"\x05\x04\x01\x01\x01\x01\x02\x02\x02\x02",
        "seed_snappy": b"\x06the quick brown fox",
        "seed_dict32": b"\x07\x01\x00\x00\x00\x02\x00\x00\x00\x01\x00\x00\x00",
        "seed_plain32": b"\x08\x01\x00\x00\x00\x02\x00\x00\x00\x03\x00\x00\x00",
        "seed_gzip": b"\x09abcdefghijklmnopqrstuvwxyz",
        "seed_zstd": b"\x0arepeat repeat repeat repeat",
    },
    "arrow_schema": {
        "seed_pyarrow": (
            b'\xff\xff\xff\xff\x88\x01\x00\x00\x10\x00\x00\x00\x00\x00\n\x00\x0c\x00\n\x00\t\x00\x04\x00\n\x00\x00\x00\x10\x00\x00\x00\x00\x01\x04\x00\x08\x00\x08\x00\x00\x00\x04\x00\x08\x00\x00\x00'
            b'\x04\x00\x00\x00\x03\x00\x00\x00\xd4\x00\x00\x00L\x00\x00\x00\x10\x00\x00\x00\x0c\x00\x10\x00\x0c\x00\x0b\x00\n\x00\x04\x00\x0c\x00\x00\x00\x1c\x00\x00\x00\x00\x00\x04\x01\x04\x00\x00\x00'
            b'\x04\x00\x00\x00Name\x00\x00\x00\x00\x04\x00\x06\x00\x04\x00\x00\x00\x00\x00\x12\x00\x14\x00\x10\x00\x00\x00\x0f\x00\x08\x00\x00\x00\x00\x00\x04\x00\x12\x00\x00\x00\x10\x00\x00\x00'
            b'\\\x00\x00\x00\x00\x00\x00\x02D\x00\x00\x00\x01\x00\x00\x00\x0c\x00\x00\x00\x08\x00\x0c\x00\x08\x00\x04\x00\x08\x00\x00\x00\x08\x00\x00\x00\x1c\x00\x00\x00\x11\x00\x00\x00Sex '
            b'of Respondent\x00\x00\x00\x05\x00\x00\x00Label\x00\x00\x00\x03\x00\x00\x00Sex\x00\x08\x00\x0e\x00\x08\x00\x07\x00\x08\x00\x00\x00'
            b'\x00\x00\x00\x01 \x00\x00\x00\x00\x00\x12\x00\x14\x00\x10\x00\x0f\x00\x0e\x00\x08\x00\x00\x00\x00\x00\x04\x00\x12\x00\x00\x00\x10\x00\x00\x00p\x00\x00\x00\x00\x00\x02\x01P\x00\x00\x00'
            b'\x01\x00\x00\x00\x0c\x00\x00\x00\x08\x00\x0c\x00\x08\x00\x04\x00\x08\x00\x00\x00\x08\x00\x00\x00(\x00\x00\x00\x1f\x00\x00\x00Numeric With Lab'
            b'els and Missing\x00\x05\x00\x00\x00Label\x00\x00\x00\t\x00\x00\x00SurveyVar\x00\x00\x00\x08\x00\x0c\x00'
            b'\x08\x00\x07\x00\x08\x00\x00\x00\x00\x00\x00\x01 \x00\x00\x00'
        ),
    },
}


# ── System info ──────────────────────────────────────────────────────────────

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


def find_clang():
    """Find a clang with libFuzzer support."""
    candidates = [
        "/opt/homebrew/opt/llvm/bin/clang",
        "/usr/local/opt/llvm/bin/clang",
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c

    clang = shutil.which("clang")
    if clang:
        try:
            out = subprocess.check_output(
                [clang, "--print-runtime-dir"],
                stderr=subprocess.DEVNULL, text=True,
            ).strip()
            if "lib/clang" in out:
                return clang
        except Exception:
            pass
    return None


def core_count():
    try:
        return os.cpu_count() or 4
    except Exception:
        return 4


# ── Build ────────────────────────────────────────────────────────────────────

def build_fuzzers(project_dir, build_dir, clang, sanitizers=True):
    print()
    print(f"  {BOLD}Building fuzz targets{RST}")
    print(f"  {_bar()}")
    print(f"  {DIM}compiler:{RST}  {clang}")
    print(f"  {DIM}build:{RST}     {build_dir}")

    san_flags = ""
    if sanitizers:
        san_flags = "-fsanitize=address,undefined -fno-omit-frame-pointer -fno-sanitize-recover=all"
        print(f"  {DIM}sanitizers:{RST} ASan + UBSan")
    else:
        print(f"  {DIM}sanitizers:{RST} disabled")

    print(f"  {_bar()}")
    print()

    # Configure
    cmake_args = [
        "cmake", "-B", build_dir,
        f"-DCMAKE_C_COMPILER={clang}",
        "-DCMAKE_BUILD_TYPE=Debug",
        f"-DCMAKE_C_FLAGS={san_flags}",
        "-DCARQUET_BUILD_FUZZ=ON",
        "-DCARQUET_BUILD_TESTS=OFF",
        "-DCARQUET_BUILD_EXAMPLES=OFF",
        "-DCARQUET_BUILD_BENCHMARKS=OFF",
        project_dir,
    ]

    print(f"  {DIM}Configuring...{RST}", flush=True)
    result = subprocess.run(cmake_args, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  {RED}CMake configure failed{RST}")
        if result.stderr:
            for line in result.stderr.strip().splitlines()[-5:]:
                print(f"    {DIM}{line}{RST}")
        return False

    # Build
    jobs = core_count()
    print(f"  {DIM}Building ({jobs} jobs)...{RST}", flush=True)
    result = subprocess.run(
        ["cmake", "--build", build_dir, f"-j{jobs}"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"  {RED}Build failed{RST}")
        if result.stderr:
            for line in result.stderr.strip().splitlines()[-10:]:
                print(f"    {DIM}{line}{RST}")
        return False

    # Verify targets
    built = []
    missing = []
    for name, _ in TARGETS:
        path = os.path.join(build_dir, "fuzz", f"fuzz_{name}")
        if os.path.isfile(path):
            built.append(name)
        else:
            missing.append(name)

    print()
    for name in built:
        print(f"  {GRN}built{RST}  fuzz_{name}")
    for name in missing:
        print(f"  {RED}missing{RST}  fuzz_{name}")
    print()

    if missing:
        print(f"  {YLW}Warning: {len(missing)} target(s) not built{RST}")
        return len(built) > 0

    print(f"  {GRN}All {len(built)} targets built successfully{RST}")
    return True


# ── Corpus management ────────────────────────────────────────────────────────

def ensure_corpus(build_dir, target, script_dir):
    """Create and seed corpus directory if needed."""
    corpus_dir = os.path.join(build_dir, "fuzz", f"corpus_{target}")
    os.makedirs(corpus_dir, exist_ok=True)

    if not os.listdir(corpus_dir):
        if target in SEEDS:
            for name, data in SEEDS[target].items():
                with open(os.path.join(corpus_dir, name), "wb") as f:
                    f.write(data)
        # Some targets (e.g. page_filter) need real binary seeds too large to
        # inline as byte literals. Copy any checked-in seeds from
        # fuzz/corpus/<target>/ if present.
        checked_in = os.path.join(script_dir, "corpus", target)
        if os.path.isdir(checked_in):
            for entry in os.scandir(checked_in):
                if entry.is_file():
                    shutil.copy2(entry.path, corpus_dir)

        # Real-world seeds from apache/parquet-testing (fetched on demand by
        # fuzz/fetch_corpus.sh into fuzz/external/). These full Parquet files
        # massively widen "valid file" coverage for the targets that consume
        # whole files. Names are prefixed to avoid collisions.
        if target in EXTERNAL_SEED_TARGETS:
            ext_root = os.path.join(script_dir, "external", "parquet-testing")
            for sub in ("data", os.path.join("data", "geospatial"), "bad_data"):
                src_dir = os.path.join(ext_root, sub)
                if not os.path.isdir(src_dir):
                    continue
                tag = sub.replace(os.sep, "_")
                for entry in os.scandir(src_dir):
                    if entry.is_file() and entry.name.endswith(".parquet"):
                        dst = os.path.join(corpus_dir, f"pt_{tag}_{entry.name}")
                        shutil.copy2(entry.path, dst)

    return corpus_dir


def corpus_stats(corpus_dir):
    """Return (file_count, total_bytes) for a corpus directory."""
    if not os.path.isdir(corpus_dir):
        return 0, 0
    count = 0
    total = 0
    for entry in os.scandir(corpus_dir):
        if entry.is_file():
            count += 1
            total += entry.stat().st_size
    return count, total


def clean_artifacts(build_dir, targets):
    """Remove saved crash artifacts for the selected targets."""
    removed = 0
    for target in targets:
        artifact_dir = os.path.join(build_dir, "fuzz", f"artifacts_{target}")
        if os.path.isdir(artifact_dir):
            shutil.rmtree(artifact_dir)
            removed += 1
    return removed


# ── Fuzzer runner ────────────────────────────────────────────────────────────

def setup_lib_paths():
    """Set library paths for Homebrew/Linuxbrew OpenMP."""
    if platform.system() == "Darwin":
        for lib_dir in ["/opt/homebrew/lib", "/usr/local/opt/llvm/lib"]:
            if os.path.isdir(lib_dir):
                current = os.environ.get("DYLD_LIBRARY_PATH", "")
                os.environ["DYLD_LIBRARY_PATH"] = f"{lib_dir}:{current}" if current else lib_dir
                return
    else:
        # Linux: look for libomp
        for candidate in ["/home/linuxbrew/.linuxbrew/lib"]:
            if os.path.isdir(candidate):
                current = os.environ.get("LD_LIBRARY_PATH", "")
                os.environ["LD_LIBRARY_PATH"] = f"{candidate}:{current}" if current else candidate
                return


def parse_crash_info(stderr_text):
    """Parse sanitizer/fuzzer stderr to extract structured crash info.

    Returns list of dicts with keys:
        bug_type: 'heap-buffer-overflow', 'undefined-behavior', etc.
        summary: one-line SUMMARY from sanitizer
        location: source file:line where the bug was detected
        frames: list of (func, file_line) tuples for the call stack
        artifact: path to the crash artifact
        input_hex: hex bytes of the crashing input
        input_b64: base64 of the crashing input
    """
    import re

    crashes = []
    current = {}

    for line in stderr_text.splitlines():
        line_s = line.strip()

        # UBSan runtime error
        m = re.match(r"(.+:\d+:\d+): runtime error: (.+)", line_s)
        if m:
            current["location"] = m.group(1)
            current["bug_type"] = "undefined-behavior"
            current["detail"] = m.group(2)

        # ASan error header
        m = re.match(r"==\d+==ERROR: AddressSanitizer: (.+?) on address", line_s)
        if m:
            current["bug_type"] = m.group(1)

        # SUMMARY line
        m = re.match(r"SUMMARY: (\w+): (.+)", line_s)
        if m:
            sanitizer = m.group(1)
            summary = m.group(2)
            if sanitizer != "libFuzzer":
                current["summary"] = summary

        # Stack frame: #N 0xADDR in func file:line
        m = re.match(r"#(\d+)\s+0x[0-9a-f]+\s+in\s+(\S+)\s+(\S+)", line_s)
        if m:
            frame_num = int(m.group(1))
            func = m.group(2)
            file_loc = m.group(3)
            if "frames" not in current:
                current["frames"] = []
            current["frames"].append((func, file_loc))
            # Use first non-sanitizer frame as location if not set
            if "location" not in current and not any(
                s in func for s in ["sanitizer", "fuzzer", "Fuzzer", "sigtramp",
                                    "pthread_kill", "abort", "__ubsan", "__asan"]
            ):
                current["location"] = f"{func} {file_loc}"

        # Artifact path
        m = re.search(r"Test unit written to (.+)", line_s)
        if m:
            current["artifact"] = m.group(1)

        # Base64 of crashing input
        m = re.match(r"Base64: (.+)", line_s)
        if m:
            current["input_b64"] = m.group(1)

        # Hex dump of crashing input
        m = re.match(r"(0x[0-9a-f]+(?:,0x[0-9a-f]+)+),?$", line_s)
        if m:
            current["input_hex"] = line_s

        # libFuzzer deadly signal = end of a crash report
        if "ERROR: libFuzzer: deadly signal" in line_s:
            if current.get("bug_type") or current.get("summary"):
                crashes.append(current)
            current = {}

    # Flush last crash if not terminated by deadly signal
    if current.get("bug_type") or current.get("summary"):
        crashes.append(current)

    return crashes


def first_interesting_frame(frames):
    """Find the first stack frame in user code (not sanitizer/fuzzer internals)."""
    skip = {"sanitizer", "fuzzer", "Fuzzer", "sigtramp", "pthread_kill",
            "abort", "__ubsan", "__asan", "Die", "Abort", "start+"}
    for func, loc in frames:
        if not any(s in func for s in skip):
            return func, loc
    return None, None


def print_crash_report(crashes, target, build_dir):
    """Print a structured crash report for one target."""
    if not crashes:
        return

    seen = set()
    unique = []
    for c in crashes:
        key = (c.get("bug_type", ""), c.get("detail", c.get("summary", "")))
        if key not in seen:
            seen.add(key)
            unique.append(c)

    for i, c in enumerate(unique):
        bug = c.get("bug_type", "unknown")
        detail = c.get("detail", c.get("summary", ""))
        loc = c.get("location", "")

        # Clean up location path to be relative
        if loc:
            for prefix in ["/Users/", "/home/"]:
                idx = loc.find("/carquet/")
                if idx >= 0:
                    loc = loc[idx + len("/carquet/"):]
                    break

        print(f"       {RED}{bug}{RST}: {detail}")
        if loc:
            print(f"       {DIM}at {loc}{RST}")

        # Show relevant stack frames (skip sanitizer internals)
        frames = c.get("frames", [])
        func, frame_loc = first_interesting_frame(frames)
        if func and func not in loc:
            if frame_loc:
                for prefix in ["/Users/", "/home/"]:
                    idx = frame_loc.find("/carquet/")
                    if idx >= 0:
                        frame_loc = frame_loc[idx + len("/carquet/"):]
                        break
            print(f"       {DIM}in {func} ({frame_loc}){RST}")

        # Show artifact path
        artifact = c.get("artifact", "")
        if artifact:
            # Make relative if possible
            rel = artifact
            idx = rel.find("build-fuzz/")
            if idx >= 0:
                rel = rel[idx:]
            print(f"       {DIM}artifact: {rel}{RST}")

        # Show crashing input (compact)
        b64 = c.get("input_b64", "")
        if b64:
            display = b64 if len(b64) <= 60 else b64[:57] + "..."
            print(f"       {DIM}input(b64): {display}{RST}")


def _parse_libfuzzer_line(line):
    """Extract live stats from a single libFuzzer stderr line.

    libFuzzer output format:
      #1234  NEW    cov: 42 ft: 156 corp: 23/1234b lim: 4096 exec/s: 567 rss: 45Mb
      #5678  REDUCE cov: 42 ft: 156 corp: 23/1200b lim: 4096 exec/s: 890 rss: 45Mb

    Returns dict with extracted values, or None if not a stats line.
    """
    line = line.strip()
    if not line.startswith("#"):
        return None

    info = {}

    # Event type: INITED, NEW, REDUCE, pulse, DONE, RELOAD
    for tag in ("INITED", "NEW", "REDUCE", "pulse", "DONE", "RELOAD"):
        if tag in line:
            info["event"] = tag
            break

    m = re.search(r"#(\d+)", line)
    if m:
        info["runs"] = int(m.group(1))

    m = re.search(r"cov:\s*(\d+)", line)
    if m:
        info["cov"] = int(m.group(1))

    m = re.search(r"ft:\s*(\d+)", line)
    if m:
        info["ft"] = int(m.group(1))

    m = re.search(r"corp:\s*(\d+)/(\d+\w?)", line)
    if m:
        info["corp_n"] = int(m.group(1))
        info["corp_size"] = m.group(2)

    m = re.search(r"exec/s:\s*(\d+)", line)
    if m:
        info["exec_s"] = int(m.group(1))

    m = re.search(r"rss:\s*(\d+)", line)
    if m:
        info["rss_mb"] = int(m.group(1))

    return info if info else None


def _format_status_line(target, stats, elapsed, max_time, crashes):
    """Build the live status line string."""
    # Progress bar
    pct = min(elapsed / max_time, 1.0) if max_time > 0 else 1.0
    bar_w = 20
    filled = int(pct * bar_w)
    bar = GRN + "\u2588" * filled + DIM + "\u2591" * (bar_w - filled) + RST

    # Time
    t_str = f"{int(elapsed):>3}s/{max_time}s"

    # Stats
    parts = []
    cov = stats.get("cov")
    if cov is not None:
        parts.append(f"cov:{cov}")
    ft = stats.get("ft")
    if ft is not None:
        parts.append(f"ft:{ft}")
    corp_n = stats.get("corp_n")
    if corp_n is not None:
        corp_s = stats.get("corp_size", "")
        parts.append(f"corp:{corp_n}/{corp_s}")
    exec_s = stats.get("exec_s")
    if exec_s is not None:
        parts.append(f"exec/s:{exec_s:,}")

    stats_str = "  ".join(parts)

    crash_str = ""
    if crashes > 0:
        crash_str = f"  {RED}\u2716 {crashes} crash{'es' if crashes > 1 else ''}{RST}"

    workers = stats.get("workers")
    workers_str = f"  {CYN}[{workers}w]{RST}" if workers else ""

    event = stats.get("event", "")
    event_str = ""
    if event == "NEW":
        event_str = f" {GRN}NEW{RST}"
    elif event == "REDUCE":
        event_str = f" {BLU}REDUCE{RST}"

    return f"\r  {bar} {DIM}{t_str}{RST}  {DIM}{stats_str}{RST}{workers_str}{event_str}{crash_str}"


def _poll_worker_logs(log_dir):
    """Read last stats from fuzz-*.log worker files and aggregate.

    With -jobs=N, libFuzzer writes each worker's output to fuzz-{JOBID}.log
    in the working directory. We read the tail of each file, parse the last
    stats line, and aggregate across all active workers.
    """
    total_runs = 0
    total_exec_s = 0
    max_cov = 0
    max_ft = 0
    max_corp_n = 0
    max_corp_size = ""
    active = 0

    for logfile in glob.glob(os.path.join(log_dir, "fuzz-*.log")):
        try:
            with open(logfile, "rb") as f:
                f.seek(0, 2)
                size = f.tell()
                f.seek(max(0, size - 4096))
                tail = f.read().decode("utf-8", errors="replace")

            for line in reversed(tail.splitlines()):
                info = _parse_libfuzzer_line(line)
                if info and "runs" in info:
                    active += 1
                    total_runs += info.get("runs", 0)
                    total_exec_s += info.get("exec_s", 0)
                    cov = info.get("cov", 0)
                    if cov > max_cov:
                        max_cov = cov
                    ft = info.get("ft", 0)
                    if ft > max_ft:
                        max_ft = ft
                    cn = info.get("corp_n", 0)
                    if cn > max_corp_n:
                        max_corp_n = cn
                        max_corp_size = info.get("corp_size", "")
                    break
        except OSError:
            pass

    if active == 0:
        return None

    return {
        "runs": total_runs,
        "exec_s": total_exec_s,
        "cov": max_cov,
        "ft": max_ft,
        "corp_n": max_corp_n,
        "corp_size": max_corp_size,
        "workers": active,
    }


def run_single_fuzzer(build_dir, target, script_dir, max_time, jobs, extra_args):
    """Run a single fuzz target with live monitoring. Returns (info_dict, corpus_count, elapsed)."""
    fuzzer = os.path.join(build_dir, "fuzz", f"fuzz_{target}")
    if not os.path.isfile(fuzzer):
        print(f"  {RED}not found:{RST} fuzz_{target}")
        return None, 0, 0

    corpus_dir = ensure_corpus(build_dir, target, script_dir)
    c_files, c_bytes = corpus_stats(corpus_dir)

    # Dictionary
    dict_file = os.path.join(script_dir, "parquet.dict")
    dict_opt = [f"-dict={dict_file}"] if os.path.isfile(dict_file) else []

    # Crash/artifact output
    artifact_dir = os.path.join(build_dir, "fuzz", f"artifacts_{target}")
    os.makedirs(artifact_dir, exist_ok=True)

    print(f"  {DIM}corpus:{RST}    {c_files} files ({c_bytes:,} bytes)")
    if dict_opt:
        print(f"  {DIM}dict:{RST}      parquet.dict")
    if jobs > 1:
        print(f"  {DIM}jobs:{RST}      {jobs}")

    # Ensure all paths are absolute (needed when we set cwd for multi-job)
    fuzzer = os.path.abspath(fuzzer)
    corpus_dir = os.path.abspath(corpus_dir)
    artifact_dir = os.path.abspath(artifact_dir)
    if dict_opt:
        dict_file = os.path.abspath(dict_file)
        dict_opt = [f"-dict={dict_file}"]

    argv = [
        fuzzer,
        corpus_dir,
        f"-max_total_time={max_time}",
        f"-artifact_prefix={artifact_dir}/",
        "-print_final_stats=1",
    ]
    argv.extend(dict_opt)
    if jobs > 1:
        argv.extend([f"-jobs={jobs}", f"-workers={jobs}"])
    argv.extend(extra_args)

    is_tty = sys.stdout.isatty()
    t0 = time.time()
    stderr_lines = []
    live_stats = {}
    live_crashes = 0
    proc = None

    try:
        # Multi-job: run in artifact_dir so fuzz-*.log files go there
        cwd = artifact_dir if jobs > 1 else None

        proc = subprocess.Popen(
            argv,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
            cwd=cwd,
        )

        if jobs <= 1:
            # ── Single-job: stream stderr line by line ──
            for line in proc.stderr:
                stderr_lines.append(line)
                elapsed = time.time() - t0

                if "ERROR: AddressSanitizer:" in line or "runtime error:" in line:
                    live_crashes += 1
                if "Test unit written to" in line:
                    live_crashes = max(live_crashes, 1)

                info = _parse_libfuzzer_line(line)
                if info:
                    live_stats.update(info)
                    if is_tty:
                        status = _format_status_line(
                            target, live_stats, elapsed, max_time, live_crashes)
                        print(status, end="", flush=True)

            proc.wait()
        else:
            # ── Multi-job: drain stderr in thread, poll worker logs ──
            def drain_stderr():
                nonlocal live_crashes
                for line in proc.stderr:
                    stderr_lines.append(line)
                    if "Test unit written to" in line:
                        live_crashes += 1

            reader = threading.Thread(target=drain_stderr, daemon=True)
            reader.start()

            while proc.poll() is None:
                time.sleep(0.5)
                elapsed = time.time() - t0

                # Aggregate stats from worker log files
                stats = _poll_worker_logs(artifact_dir)
                if stats:
                    live_stats.update(stats)

                # Count crash artifacts
                try:
                    crash_count = sum(
                        1 for f in os.scandir(artifact_dir)
                        if f.is_file() and f.name.startswith("crash-")
                    )
                    live_crashes = max(live_crashes, crash_count)
                except OSError:
                    pass

                if is_tty:
                    status = _format_status_line(
                        target, live_stats, elapsed, max_time, live_crashes)
                    print(status, end="", flush=True)

            reader.join(timeout=5)

            # Clean up worker log files
            for f in glob.glob(os.path.join(artifact_dir, "fuzz-*.log")):
                try:
                    os.remove(f)
                except OSError:
                    pass

    except KeyboardInterrupt:
        if proc and proc.poll() is None:
            proc.terminate()
            proc.wait()
        if is_tty:
            print()  # clear partial line
        print(f"  {YLW}interrupted{RST}")

    elapsed = time.time() - t0

    # Clear the live status line
    if is_tty and live_stats:
        print("\r" + " " * 100 + "\r", end="")

    stderr_text = "".join(stderr_lines)

    # Parse crash details from stderr
    crash_infos = parse_crash_info(stderr_text)

    # Count artifact files
    artifact_count = 0
    if os.path.isdir(artifact_dir):
        artifact_count = sum(
            1 for f in os.scandir(artifact_dir)
            if f.is_file() and (f.name.startswith("crash-") or f.name.startswith("leak-")
                                or f.name.startswith("timeout-"))
        )

    # Parse final stats from stderr
    cov = execs = None
    for line_s in stderr_text.splitlines():
        line_s = line_s.strip()
        if line_s.startswith("stat::number_of_executed_units:"):
            try:
                execs = int(line_s.split(":")[2].strip())
            except (IndexError, ValueError):
                pass

    # Coverage from last stats update (prefer live, fallback to parsed)
    cov = live_stats.get("cov")
    if cov is None:
        for line in reversed(stderr_text.splitlines()):
            m = re.search(r"cov: (\d+)", line)
            if m:
                cov = int(m.group(1))
                break

    c_files_after, c_bytes_after = corpus_stats(corpus_dir)

    return {
        "crashes": crash_infos,
        "artifact_count": artifact_count,
        "corpus_count": c_files_after,
        "corpus_bytes": c_bytes_after,
        "elapsed": elapsed,
        "coverage": cov,
        "execs": execs,
    }, c_files_after, elapsed


# ── Corpus minimization ─────────────────────────────────────────────────────

def minimize_corpus(build_dir, target, script_dir):
    """Minimize a corpus directory using libFuzzer's -merge=1."""
    fuzzer = os.path.join(build_dir, "fuzz", f"fuzz_{target}")
    if not os.path.isfile(fuzzer):
        print(f"  {RED}not found:{RST} fuzz_{target}")
        return False

    corpus_dir = os.path.join(build_dir, "fuzz", f"corpus_{target}")
    if not os.path.isdir(corpus_dir):
        print(f"  {YLW}no corpus to minimize{RST}")
        return False

    before_files, before_bytes = corpus_stats(corpus_dir)
    if before_files == 0:
        print(f"  {YLW}empty corpus{RST}")
        return False

    print(f"  {DIM}before:{RST}  {before_files} files ({before_bytes:,} bytes)")

    # Merge into a new directory
    merged = corpus_dir + "_merged"
    os.makedirs(merged, exist_ok=True)

    result = subprocess.run(
        [fuzzer, "-merge=1", merged, corpus_dir],
        capture_output=True, text=True, timeout=300,
    )

    if result.returncode != 0:
        print(f"  {RED}merge failed{RST}")
        shutil.rmtree(merged, ignore_errors=True)
        return False

    # Swap
    backup = corpus_dir + "_backup"
    os.rename(corpus_dir, backup)
    os.rename(merged, corpus_dir)
    shutil.rmtree(backup)

    after_files, after_bytes = corpus_stats(corpus_dir)
    reduction = (1 - after_files / before_files) * 100 if before_files else 0
    print(f"  {GRN}after:{RST}   {after_files} files ({after_bytes:,} bytes)")
    print(f"  {DIM}reduced:{RST} {reduction:.0f}%")
    return True


# ── Coverage report ──────────────────────────────────────────────────────────

def run_coverage(build_dir, project_dir, target, script_dir):
    """Build with coverage, run fuzzer briefly, generate report."""
    cov_dir = os.path.join(project_dir, "build-fuzz-cov")
    clang = find_clang()
    if not clang:
        print(f"  {RED}no clang found{RST}")
        return False

    print(f"  {DIM}Building with coverage instrumentation...{RST}")
    san_flags = "-fprofile-instr-generate -fcoverage-mapping"

    result = subprocess.run([
        "cmake", "-B", cov_dir,
        f"-DCMAKE_C_COMPILER={clang}",
        "-DCMAKE_BUILD_TYPE=Debug",
        f"-DCMAKE_C_FLAGS={san_flags}",
        "-DCARQUET_BUILD_FUZZ=ON",
        "-DCARQUET_BUILD_TESTS=OFF",
        "-DCARQUET_BUILD_EXAMPLES=OFF",
        "-DCARQUET_BUILD_BENCHMARKS=OFF",
        project_dir,
    ], capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  {RED}configure failed{RST}")
        return False

    result = subprocess.run(
        ["cmake", "--build", cov_dir, f"-j{core_count()}"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"  {RED}build failed{RST}")
        return False

    fuzzer = os.path.join(cov_dir, "fuzz", f"fuzz_{target}")
    corpus_dir = os.path.join(build_dir, "fuzz", f"corpus_{target}")
    if not os.path.isdir(corpus_dir):
        corpus_dir = ensure_corpus(cov_dir, target, script_dir)

    print(f"  {DIM}Running corpus for coverage...{RST}")
    env = os.environ.copy()
    profraw = os.path.join(cov_dir, "default.profraw")
    env["LLVM_PROFILE_FILE"] = profraw

    subprocess.run(
        [fuzzer, corpus_dir, "-max_total_time=30", "-runs=0"],
        env=env, capture_output=True, timeout=60,
    )

    if not os.path.isfile(profraw):
        print(f"  {YLW}no profile data generated{RST}")
        return False

    # Convert profile
    llvm_dir = os.path.dirname(clang)
    profdata = os.path.join(cov_dir, "default.profdata")
    profdata_bin = os.path.join(llvm_dir, "llvm-profdata")
    if not os.path.isfile(profdata_bin):
        profdata_bin = shutil.which("llvm-profdata") or "llvm-profdata"

    subprocess.run(
        [profdata_bin, "merge", "-sparse", profraw, "-o", profdata],
        capture_output=True,
    )

    # Show report
    cov_bin = os.path.join(llvm_dir, "llvm-cov")
    if not os.path.isfile(cov_bin):
        cov_bin = shutil.which("llvm-cov") or "llvm-cov"

    subprocess.run([
        cov_bin, "report", fuzzer,
        f"-instr-profile={profdata}",
        f"--ignore-filename-regex=.*/(zstd|zlib|lz4_raw)/.*",
    ])
    return True


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Carquet fuzzing orchestrator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""targets:
  reader        Full Parquet file reader
  writer        Schema creation and data writing
  compression   Compression codecs (Snappy, LZ4, GZIP, ZSTD)
  encodings     Encoding decoders (RLE, Delta, Plain, Dictionary, BSS)
  thrift        Thrift compact protocol decoder
  roundtrip     Encode-decode roundtrip consistency
  all           Run all targets sequentially

commands:
  build         Build fuzz targets only
  list          List available targets
  minimize      Minimize corpus (pass target name after)
  coverage      Generate coverage report (pass target name after)
""",
    )
    parser.add_argument("command", nargs="?", default="help",
                        help="Target name, 'all', 'build', 'list', 'minimize', or 'coverage'")
    parser.add_argument("target", nargs="?",
                        help="Target for minimize/coverage commands")
    parser.add_argument("--time", type=int, default=300,
                        help="Seconds per target (default: 300)")
    parser.add_argument("--jobs", type=int, default=1,
                        help="Parallel fuzzing jobs (default: 1)")
    parser.add_argument("--no-sanitizers", action="store_true",
                        help="Disable ASan/UBSan")
    parser.add_argument("--build-dir",
                        help="Build directory (default: <project>/build-fuzz)")
    parser.add_argument("--clean-artifacts", action="store_true",
                        help="Delete saved crash artifacts for selected targets before fuzzing")
    parser.add_argument("extra", nargs="*",
                        help="Extra libFuzzer arguments (e.g. -max_len=4096)")

    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(script_dir)
    build_dir = args.build_dir or os.path.join(project_dir, "build-fuzz")

    # ── Help ──
    if args.command in ("help", "--help", "-h"):
        parser.print_help()
        sys.exit(0)

    # ── List ──
    if args.command == "list":
        print()
        print(f"  {BOLD}Fuzz Targets{RST}")
        print(f"  {_bar()}")
        for name, desc in TARGETS:
            fuzzer = os.path.join(build_dir, "fuzz", f"fuzz_{name}")
            status = f"{GRN}built{RST}" if os.path.isfile(fuzzer) else f"{DIM}not built{RST}"
            print(f"  {status}  {BOLD}{name:<14}{RST} {DIM}{desc}{RST}")

            corpus_dir = os.path.join(build_dir, "fuzz", f"corpus_{name}")
            c_files, c_bytes = corpus_stats(corpus_dir)
            if c_files > 0:
                print(f"            {DIM}corpus: {c_files} files ({c_bytes:,} bytes){RST}")

            artifacts_dir = os.path.join(build_dir, "fuzz", f"artifacts_{name}")
            if os.path.isdir(artifacts_dir):
                crashes = sum(1 for f in os.scandir(artifacts_dir) if f.is_file())
                if crashes:
                    print(f"            {RED}artifacts: {crashes} file(s){RST}")
        print()
        sys.exit(0)

    # ── Build ──
    if args.command == "build":
        clang = find_clang()
        if not clang:
            print(f"\n  {RED}No clang with libFuzzer found.{RST}")
            print(f"  {DIM}On macOS: brew install llvm{RST}")
            print(f"  {DIM}On Linux: apt install clang{RST}\n")
            sys.exit(1)
        ok = build_fuzzers(project_dir, build_dir, clang, not args.no_sanitizers)
        sys.exit(0 if ok else 1)

    # ── Minimize ──
    if args.command == "minimize":
        target = args.target
        if not target or target not in TARGET_NAMES:
            print(f"  {RED}Usage: {sys.argv[0]} minimize <target>{RST}")
            sys.exit(1)
        print()
        print(f"  {BOLD}Minimizing corpus: {target}{RST}")
        print(f"  {_bar()}")
        minimize_corpus(build_dir, target, script_dir)
        print()
        sys.exit(0)

    # ── Coverage ──
    if args.command == "coverage":
        target = args.target
        if not target or target not in TARGET_NAMES:
            print(f"  {RED}Usage: {sys.argv[0]} coverage <target>{RST}")
            sys.exit(1)
        print()
        print(f"  {BOLD}Coverage report: {target}{RST}")
        print(f"  {_bar()}")
        run_coverage(build_dir, project_dir, target, script_dir)
        print()
        sys.exit(0)

    # ── Fuzz targets ──
    if args.command == "all":
        targets_to_run = TARGET_NAMES[:]
    elif args.command in TARGET_NAMES:
        targets_to_run = [args.command]
    else:
        print(f"  {RED}Unknown command or target: {args.command}{RST}")
        parser.print_help()
        sys.exit(1)

    setup_lib_paths()

    # Auto-build if needed
    any_missing = any(
        not os.path.isfile(os.path.join(build_dir, "fuzz", f"fuzz_{t}"))
        for t in targets_to_run
    )
    if any_missing:
        clang = find_clang()
        if not clang:
            print(f"\n  {RED}No clang with libFuzzer found.{RST}")
            print(f"  {DIM}On macOS: brew install llvm{RST}")
            print(f"  {DIM}On Linux: apt install clang{RST}\n")
            sys.exit(1)
        if not build_fuzzers(project_dir, build_dir, clang, not args.no_sanitizers):
            sys.exit(1)

    # ── Header ──
    print()
    print(f"  {BOLD}Carquet Fuzzing{RST}")
    print(f"  {_bar()}")
    print(f"  {DIM}CPU:{RST}       {get_cpu_name()}")
    print(f"  {DIM}targets:{RST}   {', '.join(targets_to_run)}")
    print(f"  {DIM}time:{RST}      {args.time}s per target")
    if args.jobs > 1:
        print(f"  {DIM}jobs:{RST}      {args.jobs}")
    if args.clean_artifacts:
        print(f"  {DIM}artifacts:{RST} clean before run")
    san_label = "disabled" if args.no_sanitizers else "ASan + UBSan"
    print(f"  {DIM}sanitizers:{RST} {san_label}")
    print(f"  {_bar()}")
    print()

    if args.clean_artifacts:
        removed = clean_artifacts(build_dir, targets_to_run)
        if removed > 0:
            noun = "directory" if removed == 1 else "directories"
            print(f"  {DIM}removed:{RST}   {removed} artifact {noun}")
            print()

    # ── Run ──
    results = {}
    total = len(targets_to_run)

    for i, target in enumerate(targets_to_run):
        tag = f"[{i+1}/{total}]"
        print(f"  {CYN}{tag}{RST} {BOLD}fuzz_{target}{RST}")

        info, corpus_count, elapsed = run_single_fuzzer(
            build_dir, target, script_dir, args.time, args.jobs, args.extra,
        )
        results[target] = info  # dict or None if not found

        if info is None:
            print(f"       {RED}skipped (not built){RST}")
        else:
            crash_list = info.get("crashes", [])
            n_crashes = len(crash_list) + info.get("artifact_count", 0)
            cov = info.get("coverage")
            execs = info.get("execs")

            if crash_list:
                print_crash_report(crash_list, target, build_dir)

            if n_crashes > 0 and not crash_list:
                print(f"       {RED}{n_crashes} crash artifact(s) found{RST}")
            elif not crash_list:
                stats_parts = [f"corpus: {info['corpus_count']} files"]
                if cov is not None:
                    stats_parts.append(f"cov: {cov}")
                if execs is not None:
                    stats_parts.append(f"execs: {execs:,}")
                stats_parts.append(f"{elapsed:.0f}s")
                stats_str = ", ".join(stats_parts)
                print(f"       {GRN}no crashes{RST}  {DIM}{stats_str}{RST}")

        if i < total - 1:
            print()

    # ── Summary ──
    print()
    dbar = _bar("\u2550")
    print(f"  {dbar}")
    print(f"  {BOLD}Summary{RST}")
    print(f"  {dbar}")
    print()

    total_crashes = 0
    total_time = 0.0
    for target in targets_to_run:
        info = results[target]

        if info is None:
            status = f"{DIM}skipped{RST}"
            corpus_n = 0
            elapsed = 0.0
            cov_str = ""
        else:
            crash_list = info.get("crashes", [])
            n_crashes = len(crash_list) + info.get("artifact_count", 0)
            elapsed = info.get("elapsed", 0.0)
            corpus_n = info.get("corpus_count", 0)
            cov = info.get("coverage")
            cov_str = f"  cov: {cov:>5}" if cov is not None else ""

            if n_crashes > 0:
                status = f"{RED}{n_crashes} CRASH(ES){RST}"
                total_crashes += n_crashes
            else:
                status = f"{GRN}clean{RST}"

        total_time += elapsed
        print(f"  {BOLD}{target:<14}{RST} {status:<28} "
              f"{DIM}corpus: {corpus_n:>5}{cov_str}  time: {elapsed:>5.0f}s{RST}")

    print(f"  {_bar()}")
    print(f"  {DIM}total time:{RST} {total_time:.0f}s  "
          f"{DIM}date:{RST} {datetime.date.today().isoformat()}")

    if total_crashes > 0:
        print(f"\n  {RED}{BOLD}{total_crashes} total crash(es) found!{RST}")
        print(f"  {DIM}Artifacts in: {build_dir}/fuzz/artifacts_<target>/{RST}")
        print(f"  {DIM}Reproduce:    ./build-fuzz/fuzz/fuzz_<target> <artifact>{RST}")
        print(f"  {DIM}Stack trace:  ASAN_OPTIONS=symbolize=1 ./fuzz_<target> <artifact>{RST}")
    else:
        print(f"\n  {GRN}All targets clean{RST}")
    print()


if __name__ == "__main__":
    main()
