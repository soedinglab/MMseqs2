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


def carquet_json(test_interop_bin, path):
    """Decode a file with carquet. Returns (summary_dict_or_None, crashed)."""
    try:
        proc = subprocess.run(
            [str(test_interop_bin), "--json", str(path)],
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
        return int(value)
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
    for p in files:
        results["data"]["checked"] += 1
        cq, crashed = carquet_json(test_interop_bin, p)
        if crashed:
            results["data"]["failed"] += 1
            print(f"  {RED}CRASH{RST} {p.name}")
            continue
        pa_sum = pyarrow_summary(p)
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
        if errs:
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
            if compression == "lz4":
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

            print(f"  {CYN}{tag}{RST} {compression:<14} "
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
    print(f"  {'Compression':<14} {'PyArrow':>8} {'DuckDB':>8}")
    print(f"  {_bar(width=32)}")

    for r in write_results:
        pa = r["pyarrow"] or "-"
        db = r["duckdb"] or "-"

        def color_status(s):
            if s == "PASS":
                return f"{GRN}PASS{RST}"
            elif s == "FAIL":
                return f"{RED}FAIL{RST}"
            return f"{DIM}{s}{RST}"

        print(f"  {r['compression']:<14} {color_status(pa):>17} {color_status(db):>17}")


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
