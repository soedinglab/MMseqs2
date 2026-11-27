#!/usr/bin/env python3
"""
Carquet interoperability test suite.

Tests bidirectional compatibility with PyArrow, DuckDB, and fastparquet:
  - Read:  Carquet reads files written by other libraries
  - Write: Other libraries read files written by Carquet

Usage:
    python3 interop/run_interop.py                    # full suite
    python3 interop/run_interop.py --read-only         # carquet reads only
    python3 interop/run_interop.py --write-only        # roundtrip only
    python3 interop/run_interop.py -v                  # verbose output
    python3 interop/run_interop.py --keep-files        # keep temp files
    python3 interop/run_interop.py --build-dir ../build
"""

import argparse
import datetime
import json
import os
import platform
import subprocess
import sys
import tempfile
import shutil
from pathlib import Path

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


# ── System info ──────────────────────────────────────────────────────────────

def get_cpu_name():
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


def get_system_info():
    return {
        "cpu": get_cpu_name(),
        "os": platform.platform(),
        "arch": platform.machine(),
    }


def get_carquet_version(build_dir):
    """Try to get carquet version from a built binary."""
    for binary in ["benchmark_carquet", "test_interop"]:
        path = os.path.join(build_dir, binary)
        if os.path.isfile(path):
            try:
                out = subprocess.check_output(
                    [path, "--version"], stderr=subprocess.DEVNULL, text=True,
                    timeout=5
                ).strip()
                if out:
                    return out
            except Exception:
                pass
    return "unknown"


# ── Library detection ────────────────────────────────────────────────────────

def detect_libraries():
    """Detect available Python libraries and their versions."""
    libs = {}
    for name in ["pyarrow", "duckdb", "fastparquet", "pandas", "numpy"]:
        try:
            mod = __import__(name)
            libs[name] = getattr(mod, "__version__", "?")
        except ImportError:
            pass
    return libs


# ── Value comparison (shared by PyArrow and DuckDB verifiers) ─────────────────

_EPOCH_DATE = datetime.date(1970, 1, 1)
_EPOCH_DT = datetime.datetime(1970, 1, 1, tzinfo=datetime.timezone.utc)


def _to_days(v):
    """Normalize a date/datetime to days since the Unix epoch."""
    if isinstance(v, datetime.datetime):
        v = v.date()
    if isinstance(v, datetime.date):
        return (v - _EPOCH_DATE).days
    return v


def _to_micros(v):
    """Normalize a datetime to microseconds since the Unix epoch (UTC)."""
    if isinstance(v, datetime.datetime):
        dt = v if v.tzinfo is not None else v.replace(tzinfo=datetime.timezone.utc)
        return round((dt - _EPOCH_DT).total_seconds() * 1_000_000)
    return v


def _to_uuid_hex(v):
    """Normalize bytes / uuid.UUID / str to lowercase dashless hex."""
    if isinstance(v, bytes):
        return v.hex()
    return str(v).replace("-", "").lower()


def compare_first_values(col_name, actual, exp, col_type, errors):
    """Compare the first values of a column against expected, by logical type."""
    if col_type == "float":
        for i, (a, b) in enumerate(zip(actual, exp)):
            if a is None and b is None:
                continue
            if a is None or b is None or abs(a - b) > 1e-4:
                errors.append(f"{col_name}[{i}]: {a} != {b}")
                break
    elif col_type == "double":
        for i, (a, b) in enumerate(zip(actual, exp)):
            if a is None and b is None:
                continue
            if a is None or b is None or abs(a - b) > 1e-10:
                errors.append(f"{col_name}[{i}]: {a} != {b}")
                break
    elif col_type == "string":
        decoded = [
            s.decode("utf-8") if isinstance(s, bytes) else s for s in actual
        ]
        if decoded != exp:
            errors.append(f"{col_name}: {decoded} != {exp}")
    elif col_type == "date":
        norm = [None if a is None else _to_days(a) for a in actual]
        if norm != exp:
            errors.append(f"{col_name}: {norm} != {exp}")
    elif col_type == "timestamp_us":
        norm = [None if a is None else _to_micros(a) for a in actual]
        if norm != exp:
            errors.append(f"{col_name}: {norm} != {exp}")
    elif col_type == "decimal":
        norm = [None if a is None else str(a) for a in actual]
        if norm != exp:
            errors.append(f"{col_name}: {norm} != {exp}")
    elif col_type == "uuid":
        norm = [None if a is None else _to_uuid_hex(a) for a in actual]
        if norm != exp:
            errors.append(f"{col_name}: {norm} != {exp}")
    else:
        if actual != exp:
            errors.append(f"{col_name}: {actual} != {exp}")


# ── Decode oracle: carquet vs PyArrow on the same file ───────────────────────
#
# carquet decodes a file (via `test_interop --json`) and PyArrow decodes the
# same file; we cross-check row counts and the first values of each column.
# PyArrow values are normalized to carquet's *physical* representation; column
# types we can't map unambiguously are skipped rather than failed, so the check
# never produces false positives on features outside the comparable set.

CORPUS_REL = Path("fuzz") / "external" / "parquet-testing"
_SKIP = object()  # sentinel: this column/value is not oracle-comparable


def find_corpus(project_dir):
    """Return the parquet-testing corpus dir if it has been fetched, else None."""
    d = project_dir / CORPUS_REL
    return d if (d / "data").is_dir() else None


def carquet_json(test_interop_bin, path, first_n=None):
    """Decode a file with carquet. Returns (summary_dict_or_None, crashed).
    first_n: emit that many values per column (across all pages / row
    groups) instead of test_interop's default first 5."""
    extra = ["--first", str(first_n)] if first_n else []
    try:
        proc = subprocess.run(
            [str(test_interop_bin), "--json", str(path)] + extra,
            capture_output=True, text=True, timeout=60
        )
    except subprocess.TimeoutExpired:
        return None, True
    if proc.returncode < 0:          # killed by a signal == crash
        return None, True
    try:
        obj = json.loads(proc.stdout.strip().splitlines()[-1])
    except (json.JSONDecodeError, IndexError):
        return None, False
    return (obj if obj.get("ok") else None), False


def _pa_normalize(col_type, value):
    """Map a PyArrow value to carquet's physical representation, or _SKIP."""
    import pyarrow as pa
    if value is None:
        return None
    if pa.types.is_boolean(col_type):
        return bool(value)
    if pa.types.is_integer(col_type):
        v = int(value)
        # carquet reports the physical (signed) INT32 / INT64 value.
        if pa.types.is_unsigned_integer(col_type) and col_type.bit_width >= 32:
            bits = col_type.bit_width
            if v >= 1 << (bits - 1):
                v -= 1 << bits
        return v
    if pa.types.is_floating(col_type):
        f = float(value)
        # carquet emits null for NaN/Inf; not meaningfully comparable.
        return _SKIP if (f != f or f in (float("inf"), float("-inf"))) else f
    if pa.types.is_string(col_type) or pa.types.is_large_string(col_type):
        return value.encode("utf-8").hex()
    if pa.types.is_binary(col_type) or pa.types.is_large_binary(col_type):
        return value.hex()
    if pa.types.is_date32(col_type):
        return (value - _EPOCH_DATE).days
    return _SKIP


def pyarrow_summary(path, first_n=5):
    """Decode with PyArrow. Returns dict {num_rows, columns, nested} or None.

    `nested` flags a schema with list/struct/map columns: carquet exposes leaf
    columns while PyArrow exposes top-level fields, so the two column models
    don't line up and the file is not oracle-comparable (it's skipped, not
    failed) — flat-file value comparison is where decode bugs actually surface.
    """
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
        table = pq.read_table(path)
    except Exception:
        return None
    nested = any(pa.types.is_nested(f.type) for f in table.schema)
    cols = []
    for field in table.schema:
        try:
            raw = table.column(field.name).to_pylist()[:first_n]
            first = [_pa_normalize(field.type, v) for v in raw]
        except Exception:
            first = [_SKIP]
        cols.append({"name": field.name, "first": first})
    return {"num_rows": table.num_rows, "columns": cols, "nested": nested}


def compare_decode(cq, pa_sum):
    """Compare a carquet summary against a PyArrow summary. Returns error list."""
    errors = []
    if cq["num_rows"] != pa_sum["num_rows"]:
        errors.append(f"rows {cq['num_rows']} != {pa_sum['num_rows']}")
    cq_cols = cq.get("columns", [])
    pa_cols = pa_sum["columns"]
    if len(cq_cols) != len(pa_cols):
        errors.append(f"columns {len(cq_cols)} != {len(pa_cols)}")
        return errors
    for cc, pc in zip(cq_cols, pa_cols):
        for i, (a, b) in enumerate(zip(cc.get("first", []), pc["first"])):
            if b is _SKIP or a == "<bin>":
                continue
            if a is None and b is None:
                continue
            if a is None or b is None:
                errors.append(f"{cc['name']}[{i}]: {a} != {b}")
                break
            if isinstance(b, float):
                if not isinstance(a, (int, float)) or abs(a - b) > 1e-6 * (1 + abs(b)):
                    errors.append(f"{cc['name']}[{i}]: {a} != {b}")
                    break
            elif a != b:
                errors.append(f"{cc['name']}[{i}]: {a} != {b}")
                break
    return errors


def run_corpus_tests(test_interop_bin, corpus_dir, verbose):
    """Read the Apache parquet-testing corpus. data/ files are oracle-checked
    against PyArrow; bad_data/ files must be rejected without crashing."""
    results = {"data": {"checked": 0, "passed": 0, "skipped": 0, "failed": 0},
               "bad": {"checked": 0, "passed": 0, "failed": 0}}

    data_dir = corpus_dir / "data"
    files = sorted(p for p in data_dir.rglob("*.parquet"))
    import pyarrow.parquet as pq
    for p in files:
        results["data"]["checked"] += 1
        # Compare every value (up to 20000 rows), not just the first five:
        # an empty carquet decode used to pass the zip()-based comparison.
        try:
            n = max(1, min(pq.read_metadata(p).num_rows, 20000))
        except Exception:
            n = 5
        cq, crashed = carquet_json(test_interop_bin, p, first_n=n)
        if crashed:
            results["data"]["failed"] += 1
            print(f"  {RED}CRASH{RST} {p.name}")
            continue
        pa_sum = pyarrow_summary(p, first_n=n)
        if cq is None or pa_sum is None:
            # One side can't read it (unsupported feature / encrypted / nested).
            results["data"]["skipped"] += 1
            if verbose:
                who = "carquet" if cq is None else "pyarrow"
                print(f"    {DIM}skip {p.name} ({who} cannot read){RST}")
            continue
        if pa_sum["nested"] or len(cq.get("columns", [])) != len(pa_sum["columns"]):
            # Nested / non-aligned column models: not oracle-comparable.
            results["data"]["skipped"] += 1
            if verbose:
                print(f"    {DIM}skip {p.name} (nested / column model differs){RST}")
            continue
        errs = compare_decode(cq, pa_sum)
        if not errs:
            for cc, pc in zip(cq["columns"], pa_sum["columns"]):
                if len(cc.get("first", [])) != len(pc["first"]):
                    errs.append(f"{cc['name']}: {len(cc.get('first', []))} values "
                                f"!= {len(pc['first'])}")
                    break
        if errs and ("corrupt" in p.name or "malformed" in p.name):
            # Intentionally broken files PyArrow tolerates (it skips CRC
            # checks by default) but carquet rejects: expected.
            results["data"]["skipped"] += 1
            if verbose:
                print(f"    {DIM}skip {p.name} (carquet rejects: {errs[0]}){RST}")
        elif errs:
            results["data"]["failed"] += 1
            print(f"  {RED}FAIL{RST} {p.name}: {errs[0]}")
        else:
            results["data"]["passed"] += 1

    bad_dir = corpus_dir / "bad_data"
    if bad_dir.is_dir():
        for p in sorted(bad_dir.rglob("*.parquet")):
            results["bad"]["checked"] += 1
            _, crashed = carquet_json(test_interop_bin, p)
            if crashed:
                results["bad"]["failed"] += 1
                print(f"  {RED}CRASH on malformed{RST} {p.name}")
            else:
                results["bad"]["passed"] += 1

    return results


def run_generated_oracle(test_interop_bin, test_dir, verbose, cap=20000):
    """Full-value oracle over the generated files: every value of every flat
    column (up to `cap` rows) decoded by carquet must equal PyArrow's decode.
    run_read_tests() only proves the files decode without error."""
    import pyarrow.parquet as pq
    res = {"checked": 0, "passed": 0, "skipped": 0, "failed": 0}
    for p in sorted(Path(test_dir).rglob("*.parquet")):
        res["checked"] += 1
        try:
            rows = pq.read_metadata(p).num_rows
        except Exception:
            res["skipped"] += 1
            continue
        n = max(1, min(rows, cap))
        cq, crashed = carquet_json(test_interop_bin, p, first_n=n)
        if crashed:
            res["failed"] += 1
            print(f"  {RED}CRASH{RST} {p.name}")
            continue
        pa_sum = pyarrow_summary(p, first_n=n)
        if (cq is None or pa_sum is None or pa_sum["nested"] or
                len(cq.get("columns", [])) != len(pa_sum["columns"])):
            res["skipped"] += 1
            if verbose:
                print(f"    {DIM}oracle skip {p.name}{RST}")
            continue
        errs = compare_decode(cq, pa_sum)
        if not errs:
            # compare_decode zips: make sure carquet emitted every value.
            for cc, pc in zip(cq["columns"], pa_sum["columns"]):
                if len(cc.get("first", [])) != len(pc["first"]):
                    errs.append(f"{cc['name']}: {len(cc.get('first', []))} values "
                                f"!= {len(pc['first'])}")
                    break
        if errs:
            res["failed"] += 1
            print(f"  {RED}FAIL{RST} {p.parent.name}/{p.name}: {errs[0]}")
        else:
            res["passed"] += 1
    return res


# ── Read tests: Carquet reads files from other libraries ─────────────────────

def run_generate(script_dir, output_dir, verbose):
    """Generate test files using generate_test_files.py."""
    gen_script = script_dir / "generate_test_files.py"
    if not gen_script.exists():
        print(f"  {RED}Error:{RST} {gen_script} not found")
        return False

    proc = subprocess.run(
        [sys.executable, str(gen_script), str(output_dir)],
        capture_output=True, text=True, timeout=120
    )
    if proc.returncode != 0:
        print(f"  {RED}Error generating test files:{RST}")
        print(proc.stderr)
        return False

    if verbose:
        for line in proc.stdout.splitlines():
            print(f"    {DIM}{line}{RST}")

    # Count generated files
    count = sum(1 for f in output_dir.rglob("*.parquet"))
    print(f"  Generated {BOLD}{count}{RST} test files")
    return True


def run_read_tests(test_interop_bin, test_dir, verbose):
    """Run test_interop to read files from other libraries."""
    cmd = [str(test_interop_bin), "--dir", str(test_dir)]
    if verbose:
        cmd.append("-v")

    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

    results = {
        "tested": 0,
        "passed": 0,
        "failed": 0,
        "output": proc.stdout,
    }

    # Parse summary from test_interop output
    for line in proc.stdout.splitlines():
        if line.startswith("Files tested:"):
            try:
                results["tested"] = int(line.split(":")[1].strip())
            except (ValueError, IndexError):
                pass
        elif line.startswith("Passed:"):
            try:
                results["passed"] = int(line.split(":")[1].strip())
            except (ValueError, IndexError):
                pass
        elif line.startswith("Failed:"):
            try:
                results["failed"] = int(line.split(":")[1].strip())
            except (ValueError, IndexError):
                pass

    if verbose:
        for line in proc.stdout.splitlines():
            print(f"    {line}")

    return results


# ── Write tests: Other libraries read Carquet output ─────────────────────────

def run_roundtrip_writer(roundtrip_bin, output_dir):
    """Run roundtrip_writer to generate carquet files, return expected JSON."""
    proc = subprocess.run(
        [str(roundtrip_bin), str(output_dir)],
        capture_output=True, text=True, timeout=60
    )
    if proc.returncode != 0:
        return None
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError:
        return None


def verify_pyarrow(path, expected, file_info):
    """Verify a carquet-written file with PyArrow."""
    try:
        import pyarrow.parquet as pq
    except ImportError:
        return ["pyarrow not available"]

    errors = []
    try:
        table = pq.read_table(path)
    except Exception as e:
        return [f"Failed to read: {e}"]

    # Per-file num_rows / verification override the global manifest values
    # (the append file carries doubled totals, for example).
    exp_rows = file_info.get("num_rows", expected["num_rows"])
    if table.num_rows != exp_rows:
        errors.append(f"Row count: {table.num_rows} != {exp_rows}")

    cols = file_info["columns"]

    # Check first values for each column
    for col_name, col_info in cols.items():
        try:
            actual = table.column(col_name).to_pylist()[:5]
            exp = col_info["first"]
            col_type = col_info.get("type", "")
            compare_first_values(col_name, actual, exp, col_type, errors)
        except Exception as e:
            errors.append(f"{col_name}: {e}")

    # Verify null counts
    verification = file_info.get("verification", expected.get("verification", {}))
    for key in ["null_count_string_col", "null_count_nullable_int"]:
        if key not in verification:
            continue
        col_name = key.replace("null_count_", "")
        try:
            actual_nulls = sum(
                1 for v in table.column(col_name).to_pylist() if v is None
            )
            if actual_nulls != verification[key]:
                errors.append(f"{col_name} nulls: {actual_nulls} != {verification[key]}")
        except Exception:
            pass

    # Verify aggregates
    if "int32_sum" in verification:
        try:
            actual_sum = sum(table.column("int32_col").to_pylist())
            if actual_sum != verification["int32_sum"]:
                errors.append(f"int32 sum: {actual_sum} != {verification['int32_sum']}")
        except Exception:
            pass

    # Verify Arrow per-field custom_metadata (variable labels) recovered from
    # the ARROW:schema footer blob.
    for fname, want in file_info.get("field_metadata", {}).items():
        try:
            field = table.schema.field(fname)
        except Exception as e:
            errors.append(f"field_metadata {fname}: field missing ({e})")
            continue
        md = {k.decode(): v.decode() for k, v in (field.metadata or {}).items()}
        for k, v in want.items():
            if md.get(k) != v:
                errors.append(
                    f"field_metadata {fname}[{k}]: {md.get(k)!r} != {v!r}")

    if not errors:
        errors.extend(verify_full_pyarrow(path, expected, file_info))
    return errors


def verify_duckdb(path, expected, file_info):
    """Verify a carquet-written file with DuckDB."""
    try:
        import duckdb
    except ImportError:
        return ["duckdb not available"]

    errors = []
    conn = None
    try:
        conn = duckdb.connect()
        result = conn.execute("SELECT * FROM read_parquet(?)", [path])
        rows = result.fetchall()
        col_names = [desc[0] for desc in result.description]
        errors.extend(verify_full_duckdb(path, expected, file_info, conn))
    except Exception as e:
        return [f"Failed to read: {e}"]
    finally:
        if conn is not None:
            conn.close()

    exp_rows = file_info.get("num_rows", expected["num_rows"])
    if len(rows) != exp_rows:
        errors.append(f"Row count: {len(rows)} != {exp_rows}")

    cols = file_info["columns"]

    # Check first values for each column without requiring pandas.
    for col_name, col_info in cols.items():
        try:
            idx = col_names.index(col_name)
            actual = [row[idx] for row in rows[:5]]
            exp = col_info["first"]
            col_type = col_info.get("type", "")
            compare_first_values(col_name, actual, exp, col_type, errors)
        except Exception as e:
            errors.append(f"{col_name}: {e}")

    verification = file_info.get("verification", expected.get("verification", {}))
    for key in ["null_count_string_col", "null_count_nullable_int"]:
        if key not in verification:
            continue
        col_name = key.replace("null_count_", "")
        try:
            idx = col_names.index(col_name)
            actual_nulls = sum(1 for row in rows if row[idx] is None)
            if actual_nulls != verification[key]:
                errors.append(f"{col_name} nulls: {actual_nulls} != {verification[key]}")
        except Exception:
            pass

    if "int32_sum" in verification:
        try:
            idx = col_names.index("int32_col")
            actual_sum = sum(row[idx] for row in rows if row[idx] is not None)
            if actual_sum != verification["int32_sum"]:
                errors.append(f"int32 sum: {actual_sum} != {verification['int32_sum']}")
        except Exception:
            pass

    if "last_int32" in verification:
        try:
            idx = col_names.index("int32_col")
            actual_last = rows[-1][idx] if rows else None
            if actual_last != verification["last_int32"]:
                errors.append(f"last int32: {actual_last} != {verification['last_int32']}")
        except Exception:
            pass

    return errors


# ── Full-value + footer verification (carquet → other libs) ──────────────────
#
# roundtrip_writer tags each flat file with a "dataset" (values are pure
# functions of the row index, regenerated here) and an "expect" block (codec,
# per-column encodings, page index, bloom columns) checked against the footer.
# The first-5 comparison above only proves the start of the first page; these
# cover every page and row group, and prove the requested encodings were
# actually applied rather than silently falling back to PLAIN.

_SAMPLE_STRINGS = ["hello", "world", "carquet", "parquet", "test",
                   "alpha", "beta", "gamma", "delta", "epsilon"]


def _basic_expected(n, copies=1):
    one = {
        "bool_col": [i % 2 == 0 for i in range(n)],
        "int32_col": [i * 10 - 5000 for i in range(n)],
        "int64_col": [i * 1000000 - 2500000000 for i in range(n)],
        "float_col": [i * 0.5 - 1250.0 for i in range(n)],
        "double_col": [i * 0.125 - 312.5 for i in range(n)],
        "string_col": [None if i % 7 == 0 else _SAMPLE_STRINGS[i % 10] for i in range(n)],
        "nullable_int": [None if i % 5 == 0 else i * 100 for i in range(n)],
    }
    return {k: v * copies for k, v in one.items()}


def _logical2_expected(n):
    """Mirror of write_logical2_file() in roundtrip_writer.c, normalized to
    raw integers (times / timestamps), unscaled ints (decimals), str."""
    m64 = (1 << 64) - 1
    return {
        "i8": [(i % 256) - 128 for i in range(n)],
        "u8": [i % 256 for i in range(n)],
        "i16_nullable": [None if i % 3 == 0 else ((i * 7) % 65536) - 32768 for i in range(n)],
        "u16": [(i * 13) % 65536 for i in range(n)],
        "u32": [(i * 2654435761) & 0xFFFFFFFF for i in range(n)],
        "u64": [(i * 11400714819323198485) & m64 for i in range(n)],
        "time_ms": [(i * 17321) % 86400000 for i in range(n)],
        "time_us": [(i * 1000003) % 86400000000 for i in range(n)],
        "time_ns": [(i * 1000000007) % 86400000000000 for i in range(n)],
        "ts_ms_local": [1609459200000 + i * 61001 for i in range(n)],
        "ts_ns_utc": [1609459200000000000 + i * 1000000123 for i in range(n)],
        "dec64": [i * 12345 - 30000000 for i in range(n)],
        "dec_flba": [(i - 2500) * 1000000007 for i in range(n)],
        "f16": [float((i % 2048) - 1024) for i in range(n)],
        "json": ['{"i": %d}' % i for i in range(n)],
        "enum": [_SAMPLE_STRINGS[i % 10] for i in range(n)],
    }


def _dataset_expected(file_info, n):
    ds = file_info.get("dataset")
    if ds == "basic":
        return _basic_expected(n)
    if ds == "basic_x2":
        return _basic_expected(n, copies=2)
    if ds == "logical2":
        return _logical2_expected(n)
    return None


def _pa_column_values(col):
    """PyArrow column -> list normalized like _logical2_expected."""
    import pyarrow as pa
    import pyarrow.compute as pc
    t = col.type
    if isinstance(t, pa.ExtensionType):
        col = pa.chunked_array([c.storage for c in col.chunks], type=t.storage_type)
        t = t.storage_type
    if pa.types.is_time32(t):
        col = pc.cast(col, pa.int32())
    elif pa.types.is_time64(t) or pa.types.is_timestamp(t):
        col = pc.cast(col, pa.int64())
    elif pa.types.is_float16(t):
        col = pc.cast(col, pa.float32())
    vals = col.to_pylist()
    if pa.types.is_decimal(t):
        vals = [None if v is None else int(v.scaleb(t.scale)) for v in vals]
    elif pa.types.is_binary(t) or pa.types.is_large_binary(t):
        vals = [None if v is None else v.decode("utf-8") for v in vals]
    return vals


def _first_diff(name, got, want, errors):
    if len(got) != len(want):
        errors.append(f"{name}: {len(got)} values != {len(want)}")
        return
    for i, (a, b) in enumerate(zip(got, want)):
        if a != b:
            errors.append(f"{name}[{i}]: {a!r} != {b!r}")
            return


_PA_CODEC_ALIASES = {"LZ4_RAW": {"LZ4_RAW", "LZ4"}, "LZ4": {"LZ4_HADOOP", "UNKNOWN"}}


def verify_full_pyarrow(path, expected, file_info):
    """Every value of every column, plus the footer "expect" block."""
    import pyarrow.parquet as pq
    errors = []
    n = expected["num_rows"]
    want = _dataset_expected(file_info, n)
    try:
        table = pq.read_table(path)
    except Exception as e:
        return [f"Failed to read: {e}"]
    if want is not None:
        for name, vals in want.items():
            try:
                _first_diff(name, _pa_column_values(table.column(name)), vals, errors)
            except Exception as e:
                errors.append(f"{name}: {e}")

    exp = file_info.get("expect")
    if not exp:
        return errors
    md = pq.ParquetFile(path).metadata
    names = [md.schema.column(i).path for i in range(md.num_columns)]
    for rg in range(md.num_row_groups):
        r = md.row_group(rg)
        for ci, cname in enumerate(names):
            cc = r.column(ci)
            # PyArrow names LZ4_RAW "LZ4" and has no name for the Hadoop
            # LZ4 codec; DuckDB's exact Parquet codec name is checked too.
            if cc.compression not in _PA_CODEC_ALIASES.get(exp["codec"], {exp["codec"]}):
                errors.append(f"rg{rg} {cname}: codec {cc.compression} != {exp['codec']}")
            want_enc = exp.get("encodings", {}).get(cname)
            if want_enc and want_enc not in cc.encodings:
                errors.append(f"rg{rg} {cname}: encodings {cc.encodings} lack {want_enc}")
            if exp.get("page_index") is not None:
                for attr in ("has_column_index", "has_offset_index"):
                    if bool(getattr(cc, attr)) != bool(exp["page_index"]):
                        errors.append(f"rg{rg} {cname}: {attr}={getattr(cc, attr)} "
                                      f"!= {exp['page_index']}")
    return errors


def verify_full_duckdb(path, expected, file_info, conn):
    """Every value of the flat columns DuckDB maps 1:1, plus bloom filters."""
    errors = []
    n = expected["num_rows"]
    want = _dataset_expected(file_info, n)
    if want is not None:
        if file_info.get("dataset") == "logical2":
            # DuckDB applies session time zones to TIME/TIMESTAMP WITH TZ;
            # those are covered by PyArrow. Compare the rest.
            skip = {"time_ms", "time_us", "time_ns", "ts_ms_local", "ts_ns_utc"}
            want = {k: v for k, v in want.items() if k not in skip}
            scale = {"dec64": 4, "dec_flba": 2}
        else:
            scale = {}
        cols = ", ".join(f'"{c}"' for c in want)
        rows = conn.execute(f"SELECT {cols} FROM read_parquet(?)", [path]).fetchall()
        for j, (name, vals) in enumerate(want.items()):
            # Un-annotated BYTE_ARRAY comes back as BLOB (bytes).
            got = [r[j].decode("utf-8") if isinstance(r[j], (bytes, bytearray)) else r[j]
                   for r in rows]
            if name in scale:
                got = [None if v is None else int(v.scaleb(scale[name])) for v in got]
            _first_diff(name, got, vals, errors)

    exp = file_info.get("expect") or {}
    bloom_cols = exp.get("bloom")
    if bloom_cols is None:
        return errors
    meta = conn.execute(
        "SELECT row_group_id, path_in_schema, bloom_filter_offset, row_group_num_rows, "
        "compression FROM parquet_metadata(?)", [path]).fetchall()
    rg_rows = {}
    for rg, col, off, nrows, codec in meta:
        rg_rows[rg] = nrows
        if exp.get("codec") and codec != exp["codec"]:
            errors.append(f"rg{rg} {col}: DuckDB codec {codec} != {exp['codec']}")
        if (off is not None) != (col in bloom_cols):
            errors.append(f"rg{rg} {col}: bloom_filter_offset={off} but bloom expected="
                          f"{col in bloom_cols}")
    if not bloom_cols or want is None:
        return errors
    starts, acc = {}, 0
    for rg in sorted(rg_rows):
        starts[rg] = acc
        acc += rg_rows[rg]

    def rg_of(i):
        return max(rg for rg in starts if starts[rg] <= i % n)

    absent = {"int32_col": [7 + 10 * k for k in range(16)],
              "int64_col": [123 + 1000000 * k for k in range(16)],
              "string_col": [f"absent-{k}" for k in range(16)]}
    excluded = probes = 0
    for col in bloom_cols:
        vals = want[col]
        for i in range(0, len(vals), max(1, len(vals) // 25)):
            v = vals[i]
            if v is None:
                continue
            res = dict(conn.execute("SELECT row_group_id, bloom_filter_excludes FROM "
                                    "parquet_bloom_probe(?, ?, ?)", [path, col, v]).fetchall())
            if res.get(rg_of(i), False):
                errors.append(f"bloom {col}: present value {v!r} (row {i}) excluded")
                break
        for v in absent[col]:
            for _, ex in conn.execute("SELECT row_group_id, bloom_filter_excludes FROM "
                                      "parquet_bloom_probe(?, ?, ?)", [path, col, v]).fetchall():
                probes += 1
                excluded += bool(ex)
    if probes and excluded < probes // 2:
        errors.append(f"bloom: only {excluded}/{probes} absent-value probes excluded")
    return errors


# ── Nested write verification (carquet → other libs read nested) ─────────────
#
# The nested fixture (LIST + STRUCT) has a fixed expected structure; the column
# model is nested rather than flat, so it gets a dedicated verifier instead of
# the column-by-column "first values" path.

_NESTED_ID = [1, 2, 3, 4]
_NESTED_TAGS = [[100, 200], None, [300], [400, 500, 600]]
_NESTED_INFO = [{"name": "a", "age": 10}, {"name": "b", "age": 20},
                {"name": "c", "age": 30}, {"name": "d", "age": 40}]


def _norm_struct(v):
    """Normalize a struct value (dict / pyarrow scalar) for comparison."""
    if isinstance(v, dict):
        name = v.get("name")
        if isinstance(name, bytes):
            name = name.decode("utf-8")
        return {"name": name, "age": v.get("age")}
    return v


_NESTED_DEEP_MATRIX = [[[1, 2], [3]], [], [[4]]]  # list<list<int32>>


def verify_nested_pyarrow(path):
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
        t = pq.read_table(path)
    except ImportError:
        return ["pyarrow not available"]
    except Exception as e:
        return [f"Failed to read: {e}"]
    errors = []
    # Deep-nested file (list<list<int32>>) written via carquet_writer_write_arrow.
    if "matrix" in t.schema.names:
        try:
            got = t.column("matrix").to_pylist()
            if got != _NESTED_DEEP_MATRIX:
                errors.append(f"matrix {got} != {_NESTED_DEEP_MATRIX}")
            mt = t.schema.field("matrix").type
            if not (pa.types.is_list(mt) and pa.types.is_list(mt.value_type)):
                errors.append(f"matrix not list<list>: {mt}")
        except Exception as e:
            errors.append(str(e))
        return errors
    try:
        if t.column("id").to_pylist() != _NESTED_ID:
            errors.append("id mismatch")
        if t.column("tags").to_pylist() != _NESTED_TAGS:
            errors.append(f"tags {t.column('tags').to_pylist()} != {_NESTED_TAGS}")
        info = [_norm_struct(v) for v in t.column("info").to_pylist()]
        if info != _NESTED_INFO:
            errors.append(f"info {info} != {_NESTED_INFO}")
        # carquet now emits ARROW:schema for nested schemas: the raw footer must
        # carry the blob (PyArrow consumes it into schema types, so check the
        # raw file metadata) and PyArrow must read the correct nested types.
        raw_meta = pq.read_metadata(path).metadata or {}
        if b"ARROW:schema" not in raw_meta:
            errors.append("ARROW:schema blob missing for nested file")
        import pyarrow as pa
        if not pa.types.is_list(t.schema.field("tags").type):
            errors.append(f"tags not list: {t.schema.field('tags').type}")
        if not pa.types.is_struct(t.schema.field("info").type):
            errors.append(f"info not struct: {t.schema.field('info').type}")
    except Exception as e:
        errors.append(str(e))
    return errors


def verify_nested_duckdb(path):
    try:
        import duckdb
    except ImportError:
        return ["duckdb not available"]
    conn = None
    try:
        conn = duckdb.connect()
        cols = [d[0] for d in conn.execute(
            "SELECT * FROM read_parquet(?) LIMIT 0", [path]).description]
        if "matrix" in cols:
            rows = conn.execute(
                "SELECT matrix FROM read_parquet(?)", [path]).fetchall()
            got = [r[0] for r in rows]
            if got != _NESTED_DEEP_MATRIX:
                return [f"matrix {got} != {_NESTED_DEEP_MATRIX}"]
            return []
        rows = conn.execute(
            "SELECT id, tags, info FROM read_parquet(?) ORDER BY id", [path]
        ).fetchall()
    except Exception as e:
        return [f"Failed to read: {e}"]
    finally:
        if conn is not None:
            conn.close()
    errors = []
    if [r[0] for r in rows] != _NESTED_ID:
        errors.append("id mismatch")
    if [r[1] for r in rows] != _NESTED_TAGS:
        errors.append(f"tags {[r[1] for r in rows]} != {_NESTED_TAGS}")
    info = [_norm_struct(r[2]) for r in rows]
    if info != _NESTED_INFO:
        errors.append(f"info {info} != {_NESTED_INFO}")
    return errors


def run_write_tests(roundtrip_bin, verbose):
    """Run roundtrip tests: carquet writes, other libraries verify."""
    results = []

    with tempfile.TemporaryDirectory() as tmpdir:
        expected = run_roundtrip_writer(roundtrip_bin, tmpdir)
        if not expected:
            print(f"  {RED}Error:{RST} roundtrip_writer failed")
            return results

        num_files = len(expected.get("files", []))
        print(f"  Generated {BOLD}{num_files}{RST} files "
              f"({expected.get('num_rows', '?')} rows each)")

        for i, file_info in enumerate(expected.get("files", [])):
            path = file_info["path"]
            compression = file_info["compression"]
            tag = f"[{i + 1}/{num_files}]"

            entry = {
                "compression": compression,
                "pyarrow": None,
                "duckdb": None,
            }

            is_nested = file_info.get("nested", False)

            # PyArrow verification
            pa_errors = (verify_nested_pyarrow(path) if is_nested
                         else verify_pyarrow(path, expected, file_info))
            if pa_errors and pa_errors != ["pyarrow not available"]:
                entry["pyarrow"] = "FAIL"
                pa_status = f"{RED}FAIL{RST}"
                if verbose:
                    for e in pa_errors:
                        print(f"      {RED}{e}{RST}")
            elif pa_errors:
                entry["pyarrow"] = "SKIP"
                pa_status = f"{DIM}skip{RST}"
            else:
                entry["pyarrow"] = "PASS"
                pa_status = f"{GRN}OK{RST}"

            # DuckDB verification. DuckDB does not support the deprecated
            # Hadoop-framed LZ4 (Parquet codec 5); PyArrow does, which is the
            # interop proof for that wire format.
            if compression == "lz4" or not file_info.get("expect", {}).get("duckdb", True):
                # Also BYTE_STREAM_SPLIT on integers ("duckdb": false), which
                # DuckDB 1.4 cannot decode.
                db_errors = ["duckdb not available"]  # treated as SKIP below
            else:
                db_errors = (verify_nested_duckdb(path) if is_nested
                             else verify_duckdb(path, expected, file_info))
            if db_errors and db_errors != ["duckdb not available"]:
                entry["duckdb"] = "FAIL"
                db_status = f"{RED}FAIL{RST}"
                if verbose:
                    for e in db_errors:
                        print(f"      {RED}{e}{RST}")
            elif db_errors:
                entry["duckdb"] = "SKIP"
                db_status = f"{DIM}skip{RST}"
            else:
                entry["duckdb"] = "PASS"
                db_status = f"{GRN}OK{RST}"

            print(f"  {CYN}{tag}{RST} {compression:<28} "
                  f"PyArrow {pa_status}  DuckDB {db_status}")

            results.append(entry)

    return results


# ── Summary ──────────────────────────────────────────────────────────────────

def print_read_summary(read_results):
    """Print read test summary."""
    tested = read_results["tested"]
    passed = read_results["passed"]
    failed = read_results["failed"]

    if failed > 0:
        color = RED
    elif passed == tested and tested > 0:
        color = GRN
    else:
        color = YLW

    print(f"  {BOLD}Read:{RST}   {color}{passed}/{tested} passed{RST}", end="")
    if failed > 0:
        print(f"  ({RED}{failed} failed{RST})", end="")
    oracle = read_results.get("oracle")
    if oracle:
        ocolor = RED if oracle["failed"] else GRN
        print(f"; full-value oracle {ocolor}{oracle['passed']} match{RST}"
              f" ({oracle['skipped']} nested skipped, {oracle['failed']} failed)", end="")
    print()


def print_write_summary(write_results):
    """Print write test summary table."""
    if not write_results:
        print(f"  {BOLD}Write:{RST}  {DIM}skipped{RST}")
        return

    pa_pass = sum(1 for r in write_results if r["pyarrow"] == "PASS")
    pa_total = sum(1 for r in write_results if r["pyarrow"] != "SKIP")
    db_pass = sum(1 for r in write_results if r["duckdb"] == "PASS")
    db_total = sum(1 for r in write_results if r["duckdb"] != "SKIP")

    pa_color = GRN if pa_pass == pa_total and pa_total > 0 else RED
    db_color = GRN if db_pass == db_total and db_total > 0 else RED

    parts = []
    if pa_total > 0:
        parts.append(f"PyArrow {pa_color}{pa_pass}/{pa_total}{RST}")
    if db_total > 0:
        parts.append(f"DuckDB {db_color}{db_pass}/{db_total}{RST}")

    print(f"  {BOLD}Write:{RST}  {', '.join(parts)}")

    # Detail table
    print()
    width = max([14] + [len(r["compression"]) for r in write_results])
    print(f"  {'Variant':<{width}} {'PyArrow':>8} {'DuckDB':>8}")
    print(f"  {_bar(width=width + 18)}")

    for r in write_results:
        pa = r["pyarrow"] or "-"
        db = r["duckdb"] or "-"

        def color_status(s):
            # Pad before colouring: ANSI codes would otherwise count toward
            # the field width (and they are empty when stdout is not a TTY).
            color = GRN if s == "PASS" else RED if s == "FAIL" else DIM
            return f"{color}{s:>8}{RST}"

        print(f"  {r['compression']:<{width}} {color_status(pa)} {color_status(db)}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Carquet interoperability test suite"
    )
    parser.add_argument("--read-only", action="store_true",
                        help="Only test reading (carquet reads other libs' files)")
    parser.add_argument("--write-only", action="store_true",
                        help="Only test writing (other libs read carquet files)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose output")
    parser.add_argument("--keep-files", action="store_true",
                        help="Keep generated test files after run")
    parser.add_argument("--build-dir", default=None,
                        help="Path to CMake build directory (default: auto-detect)")
    args = parser.parse_args()

    script_dir = Path(__file__).parent
    project_dir = script_dir.parent

    # Find build directory
    if args.build_dir:
        build_dir = Path(args.build_dir)
    else:
        build_dir = project_dir / "build"

    # Locate binaries
    test_interop_bin = build_dir / "test_interop"
    roundtrip_bin = build_dir / "roundtrip_writer"

    need_read = not args.write_only
    need_write = not args.read_only

    if need_read and not test_interop_bin.is_file():
        print(f"{RED}Error:{RST} {test_interop_bin} not found. Build first:")
        print(f"  cmake -B build -DCARQUET_BUILD_INTEROP=ON && cmake --build build")
        sys.exit(1)

    if need_write and not roundtrip_bin.is_file():
        print(f"{RED}Error:{RST} {roundtrip_bin} not found. Build first:")
        print(f"  cmake -B build -DCARQUET_BUILD_INTEROP=ON && cmake --build build")
        sys.exit(1)

    # Detect libraries
    libs = detect_libraries()
    sys_info = get_system_info()
    carquet_ver = get_carquet_version(str(build_dir))

    # The full read suite cross-checks against the Apache parquet-testing
    # corpus; without it we run a degraded suite over the synthetic files only.
    corpus_dir = find_corpus(project_dir)
    have_pa = "pyarrow" in libs

    # ── Header ──
    print()
    print(f"  {BOLD}Carquet {carquet_ver} Interoperability Tests{RST}")
    print(f"  {_bar()}")
    print(f"  {DIM}CPU:{RST}  {sys_info['cpu']}")
    print(f"  {DIM}OS:{RST}   {sys_info['os']}")
    lib_strs = [f"{k} {v}" for k, v in libs.items()
                if k in ("pyarrow", "duckdb", "fastparquet")]
    if lib_strs:
        print(f"  {DIM}Libs:{RST} {', '.join(lib_strs)}")
    if need_read:
        if corpus_dir and have_pa:
            n = sum(1 for _ in (corpus_dir / "data").rglob("*.parquet"))
            print(f"  {DIM}Mode:{RST} {GRN}full{RST} "
                  f"{DIM}(Apache parquet-testing corpus: {n} files){RST}")
        else:
            reason = "corpus not fetched" if not corpus_dir else "pyarrow missing"
            print(f"  {DIM}Mode:{RST} {YLW}degraded{RST} {DIM}({reason}){RST}")
            if not corpus_dir:
                print(f"  {YLW}Run {BOLD}fuzz/fetch_corpus.sh{RST}{YLW} for the "
                      f"full interop suite (real-world Parquet files).{RST}")
    print(f"  {_bar()}")
    print()

    read_results = None
    write_results = None
    corpus_results = None
    has_failures = False

    # ── Read phase ──
    if need_read:
        print(f"  {BOLD}Phase 1: Read Tests{RST} {DIM}(carquet reads other libs){RST}")
        print(f"  {_bar(width=42)}")

        # Generate test files
        if args.keep_files:
            test_dir = project_dir / "interop" / "test_files"
            test_dir.mkdir(exist_ok=True)
            cleanup_read = False
        else:
            test_dir = Path(tempfile.mkdtemp(prefix="carquet_interop_read_"))
            cleanup_read = True

        try:
            if not run_generate(script_dir, test_dir, args.verbose):
                has_failures = True
            else:
                read_results = run_read_tests(
                    test_interop_bin, test_dir, args.verbose
                )
                if read_results["failed"] > 0:
                    has_failures = True
                if have_pa:
                    oracle = run_generated_oracle(test_interop_bin, test_dir,
                                                  args.verbose)
                    read_results["oracle"] = oracle
                    print(f"  Full-value oracle: {oracle['passed']} files match "
                          f"PyArrow, {oracle['skipped']} skipped (nested), "
                          f"{oracle['failed']} failed")
                    if oracle["failed"] > 0:
                        has_failures = True
        finally:
            if cleanup_read and test_dir.exists():
                shutil.rmtree(test_dir, ignore_errors=True)

        print()

        # ── Corpus phase (full mode only) ──
        if corpus_dir and have_pa:
            print(f"  {BOLD}Phase 1b: Corpus{RST} "
                  f"{DIM}(Apache parquet-testing, oracle = PyArrow){RST}")
            print(f"  {_bar(width=42)}")
            corpus_results = run_corpus_tests(test_interop_bin, corpus_dir,
                                              args.verbose)
            if (corpus_results["data"]["failed"] > 0 or
                    corpus_results["bad"]["failed"] > 0):
                has_failures = True
            d = corpus_results["data"]
            b = corpus_results["bad"]
            print(f"  data: {GRN}{d['passed']} ok{RST}, "
                  f"{DIM}{d['skipped']} skipped{RST}"
                  + (f", {RED}{d['failed']} failed{RST}" if d['failed'] else "")
                  + f"  ({d['checked']} files)")
            print(f"  bad_data: {GRN}{b['passed']} rejected cleanly{RST}"
                  + (f", {RED}{b['failed']} crashed{RST}" if b['failed'] else "")
                  + f"  ({b['checked']} files)")
            print()

    # ── Write phase ──
    if need_write:
        phase_num = "2" if need_read else "1"
        print(f"  {BOLD}Phase {phase_num}: Write Tests{RST} "
              f"{DIM}(other libs read carquet){RST}")
        print(f"  {_bar(width=42)}")

        write_results = run_write_tests(roundtrip_bin, args.verbose)
        write_failures = sum(
            1 for r in write_results
            if r["pyarrow"] == "FAIL" or r["duckdb"] == "FAIL"
        )
        if write_failures > 0:
            has_failures = True

        print()

    # ── Summary ──
    print("  " + _bar("\u2550"))
    print(f"  {BOLD}Summary{RST}")
    print("  " + _bar("\u2550"))
    print()

    if read_results:
        print_read_summary(read_results)

    if corpus_results is not None:
        d = corpus_results["data"]
        b = corpus_results["bad"]
        color = RED if (d["failed"] or b["failed"]) else GRN
        print(f"  {BOLD}Corpus:{RST} {color}{d['passed']} verified{RST}, "
              f"{d['skipped']} skipped, {d['failed']} failed; "
              f"{b['passed']}/{b['checked']} malformed handled")

    if write_results is not None:
        print_write_summary(write_results)
    elif need_write:
        print(f"  {BOLD}Write:{RST}  {RED}failed to generate{RST}")

    print()

    if has_failures:
        print(f"  {RED}{BOLD}SOME TESTS FAILED{RST}")
    else:
        print(f"  {GRN}{BOLD}ALL TESTS PASSED{RST}")

    # ── JSON report ──
    report = {
        "carquet_version": carquet_ver,
        "timestamp": datetime.date.today().isoformat(),
        "system": sys_info,
        "libraries": libs,
    }
    if read_results:
        report["read"] = {
            "tested": read_results["tested"],
            "passed": read_results["passed"],
            "failed": read_results["failed"],
        }
        if "oracle" in read_results:
            report["read"]["oracle"] = read_results["oracle"]
    if corpus_results is not None:
        report["corpus"] = corpus_results
    if write_results is not None:
        report["write"] = [
            {
                "compression": r["compression"],
                "pyarrow": r["pyarrow"],
                "duckdb": r["duckdb"],
            }
            for r in write_results
        ]

    date_str = datetime.date.today().strftime("%Y%m%d")
    os_tag = f"{platform.system().lower()}_{platform.machine()}"
    json_name = f"interop_{carquet_ver}_{os_tag}_{date_str}.json"
    json_path = script_dir / json_name
    with open(json_path, "w") as f:
        json.dump(report, f, indent=2)

    print(f"  {DIM}Saved:{RST} {json_name}")
    print()

    sys.exit(1 if has_failures else 0)


if __name__ == "__main__":
    main()
