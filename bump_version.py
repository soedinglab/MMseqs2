#!/usr/bin/env python3
"""Bump the project version.

Detects the current version from CMakeLists.txt, finds every occurrence in
git-tracked files, and rewrites them to the new version.

Usage:
    python bump_version.py --version 0.4.1
    python bump_version.py --version 0.4.1 --dry-run
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent

# Files whose historical version mentions must NOT be rewritten.
SKIP_PATHS = {
    "CHANGELOG.md",
    "bump_version.py",
}
SKIP_PREFIXES = ("interop/interop_",)  # dated snapshot files

VERSION_HEADER = REPO / "include" / "carquet" / "carquet.h"
MACRO_RE = {
    "MAJOR": re.compile(r"(#define\s+CARQUET_VERSION_MAJOR\s+)\d+"),
    "MINOR": re.compile(r"(#define\s+CARQUET_VERSION_MINOR\s+)\d+"),
    "PATCH": re.compile(r"(#define\s+CARQUET_VERSION_PATCH\s+)\d+"),
}

# xmake.lua carries its own `set_version("X.Y.Z")`. It is handled explicitly
# (not via the generic git-tracked pass) so the bump works even when the file
# is untracked, and so only the set_version line is touched — never some other
# version-looking string that may live in the build script.
XMAKE_FILE = REPO / "xmake.lua"
XMAKE_RE = re.compile(r'(set_version\(")\d+\.\d+\.\d+("\))')


def current_version() -> str:
    text = (REPO / "CMakeLists.txt").read_text()
    m = re.search(r"project\([^)]*VERSION\s+(\d+\.\d+\.\d+)", text)
    if not m:
        sys.exit("error: could not find project(... VERSION ...) in CMakeLists.txt")
    return m.group(1)


def tracked_files() -> list[Path]:
    out = subprocess.check_output(["git", "ls-files"], cwd=REPO, text=True)
    return [REPO / p for p in out.splitlines() if p]


def is_skipped(rel: str) -> bool:
    return rel in SKIP_PATHS or rel.startswith(SKIP_PREFIXES)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", required=True, help="new version (X.Y.Z)")
    ap.add_argument("--dry-run", action="store_true", help="show changes, do not write")
    args = ap.parse_args()

    if not re.fullmatch(r"\d+\.\d+\.\d+", args.version):
        sys.exit(f"error: --version must be X.Y.Z, got {args.version!r}")

    old = current_version()
    new = args.version
    if old == new:
        print(f"already at {old}, nothing to do")
        return 0

    pattern = re.compile(r"\b" + re.escape(old) + r"\b")
    hits: list[tuple[Path, int, str]] = []
    skipped: list[tuple[Path, int, str]] = []

    for path in tracked_files():
        rel = path.relative_to(REPO).as_posix()
        try:
            text = path.read_text()
        except (UnicodeDecodeError, FileNotFoundError):
            continue
        if not pattern.search(text):
            continue
        bucket = skipped if is_skipped(rel) else hits
        for i, line in enumerate(text.splitlines(), 1):
            if pattern.search(line):
                bucket.append((path, i, line))

    print(f"current version: {old}")
    print(f"new version:     {new}")
    print()

    files_to_change = sorted({p for p, _, _ in hits})
    print(f"will rewrite {len(hits)} occurrence(s) in {len(files_to_change)} file(s):")
    for p, i, line in hits:
        print(f"  {p.relative_to(REPO)}:{i}: {line.strip()}")

    major, minor, patch = new.split(".")
    header_text = VERSION_HEADER.read_text()
    macro_updates = []
    for name, value in (("MAJOR", major), ("MINOR", minor), ("PATCH", patch)):
        m = MACRO_RE[name].search(header_text)
        if m and m.group(0) != f"{m.group(1)}{value}":
            macro_updates.append((name, m.group(0), f"{m.group(1)}{value}"))
    if macro_updates:
        print()
        print(f"will update {len(macro_updates)} macro(s) in {VERSION_HEADER.relative_to(REPO)}:")
        for name, old_line, new_line in macro_updates:
            print(f"  {name}: {old_line!r} -> {new_line!r}")

    xmake_text = XMAKE_FILE.read_text() if XMAKE_FILE.exists() else None
    xmake_m = XMAKE_RE.search(xmake_text) if xmake_text is not None else None
    if xmake_m and xmake_m.group(0) != f"{xmake_m.group(1)}{new}{xmake_m.group(2)}":
        print()
        print(f"will update set_version in {XMAKE_FILE.relative_to(REPO)}:")
        print(f"  {xmake_m.group(0)!r} -> {xmake_m.group(1) + new + xmake_m.group(2)!r}")

    if skipped:
        print()
        print(f"skipped {len(skipped)} occurrence(s) in historical files:")
        for p, i, line in skipped:
            print(f"  {p.relative_to(REPO)}:{i}: {line.strip()}")

    if args.dry_run or not (hits or macro_updates or xmake_m):
        return 0

    for path in files_to_change:
        text = path.read_text()
        path.write_text(pattern.sub(new, text))

    if xmake_m:
        new_xmake = XMAKE_RE.sub(rf"\g<1>{new}\g<2>", xmake_text)
        if new_xmake != xmake_text:
            XMAKE_FILE.write_text(new_xmake)
            print(f"  updated set_version in {XMAKE_FILE.relative_to(REPO)}")

    major, minor, patch = new.split(".")
    header_text = VERSION_HEADER.read_text()
    new_header = MACRO_RE["MAJOR"].sub(rf"\g<1>{major}", header_text)
    new_header = MACRO_RE["MINOR"].sub(rf"\g<1>{minor}", new_header)
    new_header = MACRO_RE["PATCH"].sub(rf"\g<1>{patch}", new_header)
    if new_header != header_text:
        VERSION_HEADER.write_text(new_header)
        print(f"  updated MAJOR/MINOR/PATCH macros in {VERSION_HEADER.relative_to(REPO)}")

    print()
    print(f"bumped {old} -> {new}. review with: git diff")
    return 0


if __name__ == "__main__":
    sys.exit(main())
