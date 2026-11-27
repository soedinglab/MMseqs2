#!/usr/bin/env python3
"""Bump the project version.

Every place that carries the *current* version is listed explicitly in
TARGETS below and rewritten in place. Other mentions of the old version
(benchmark captions, report file names, changelog entries) are history and
are only listed for review, never rewritten.

Derived automatically, nothing to edit: the Python wheel version (read from
CARQUET_VERSION_STRING by scikit-build-core), carquet.pc and the CMake
package config (@PROJECT_VERSION@), CARQUET_VERSION_NUMBER, the release
workflows (version taken from the pushed tag).

Usage:
    python3 bump_version.py --version 0.7.4            # rewrite
    python3 bump_version.py --version 0.7.4 --dry-run  # show what would change
    python3 bump_version.py --check                    # all targets agree?
"""

import argparse
import datetime
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent
V = r"\d+\.\d+\.\d+"

# FetchContent snippets that pin carquet to the latest release tag.
PIN_RE = r"(GIT_REPOSITORY\s+https://github\.com/Vitruves/carquet(?:\.git)?\s+GIT_TAG\s+v)(" + V + r")"

# (path, regex with the version in group 2, replacement builder, required)
# `rep(new, parts, today)` returns the text for group 2.
TARGETS = [
    ("CMakeLists.txt", r"(project\(carquet\s+VERSION\s+)(" + V + r")", lambda n, p, d: n, True),
    ("xmake.lua", r'(set_version\(")(' + V + r')', lambda n, p, d: n, True),
    ("include/carquet/carquet.h", r'(#define\s+CARQUET_VERSION_STRING\s+")(' + V + r')', lambda n, p, d: n, True),
    ("include/carquet/carquet.h", r"(#define\s+CARQUET_VERSION_MAJOR\s+)(\d+)", lambda n, p, d: p[0], True),
    ("include/carquet/carquet.h", r"(#define\s+CARQUET_VERSION_MINOR\s+)(\d+)", lambda n, p, d: p[1], True),
    ("include/carquet/carquet.h", r"(#define\s+CARQUET_VERSION_PATCH\s+)(\d+)", lambda n, p, d: p[2], True),
    ("Doxyfile", r"(PROJECT_NUMBER\s*=\s*)(" + V + r")", lambda n, p, d: n, True),
    ("CITATION.cff", r"(^version:\s*)(" + V + r")", lambda n, p, d: n, True),
    ("CITATION.cff", r'(^date-released:\s*")(\d{4}-\d{2}-\d{2})', lambda n, p, d: d, True),
    ("site/index.html", PIN_RE, lambda n, p, d: n, False),
    (".agents/skills/carquet/references/integration.md", PIN_RE, lambda n, p, d: n, False),
    (".claude/skills/carquet/references/integration.md", PIN_RE, lambda n, p, d: n, False),
]

# Version-only targets used by --check (the date and pins are not versions of the tree).
CHECK_SKIP = {"date-released", "GIT_TAG"}

# Historical files whose old-version mentions are expected and not worth listing.
REVIEW_SKIP = ("CHANGELOG.md", "bump_version.py", "benchmark/bench_", "interop/interop_")


def current_version() -> str:
    m = re.search(r"project\(carquet\s+VERSION\s+(" + V + r")", (REPO / "CMakeLists.txt").read_text())
    if not m:
        sys.exit("error: no project(carquet VERSION X.Y.Z) in CMakeLists.txt")
    return m.group(1)


def check() -> int:
    seen = []
    for rel, pat, _, required in TARGETS:
        if any(k in pat for k in CHECK_SKIP):
            continue
        path = REPO / rel
        m = re.search(pat, path.read_text(), re.M) if path.exists() else None
        if not m:
            if required:
                seen.append((rel, pat, None))
            continue
        seen.append((rel, pat, m.group(2)))
    head = REPO / "include/carquet/carquet.h"
    ht = head.read_text()
    parts = [re.search(r"CARQUET_VERSION_%s\s+(\d+)" % k, ht).group(1) for k in ("MAJOR", "MINOR", "PATCH")]
    want = current_version()
    ok = True
    for rel, pat, val in seen:
        if "MAJOR" in pat or "MINOR" in pat or "PATCH" in pat:
            continue
        flag = "ok" if val == want else "MISMATCH"
        ok &= val == want
        print(f"  {flag:8} {rel}: {val}")
    macro = ".".join(parts)
    print(f"  {'ok' if macro == want else 'MISMATCH':8} include/carquet/carquet.h MAJOR.MINOR.PATCH: {macro}")
    ok &= macro == want
    print(f"\n{'all targets at ' + want if ok else 'version mismatch'}")
    return 0 if ok else 1


def tracked_files() -> list[str]:
    out = subprocess.check_output(["git", "ls-files", "--cached", "--others", "--exclude-standard"], cwd=REPO, text=True)
    return [p for p in out.splitlines() if p]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--version", help="new version (X.Y.Z)")
    ap.add_argument("--date", default=datetime.date.today().isoformat(), help="CITATION.cff date-released (default: today)")
    ap.add_argument("--dry-run", action="store_true", help="show changes, do not write")
    ap.add_argument("--check", action="store_true", help="verify every target carries the CMakeLists.txt version")
    args = ap.parse_args()

    if args.check:
        return check()
    if not args.version or not re.fullmatch(V, args.version):
        sys.exit("error: --version X.Y.Z is required (or use --check)")

    old, new = current_version(), args.version
    parts = new.split(".")
    print(f"current version: {old}\nnew version:     {new}\n")

    texts: dict[str, str] = {}
    changes = []
    for rel, pat, rep, required in TARGETS:
        path = REPO / rel
        if not path.exists():
            if required:
                sys.exit(f"error: {rel} not found")
            continue
        text = texts.setdefault(rel, path.read_text())
        rx = re.compile(pat, re.M)
        if not rx.search(text):
            if required:
                sys.exit(f"error: version pattern not found in {rel}: {pat}")
            continue
        value = rep(new, parts, args.date)
        for m in rx.finditer(text):
            if m.group(2) != value:
                line = text.count("\n", 0, m.start(2)) + 1
                changes.append(f"  {rel}:{line}: {m.group(2)} -> {value}")
        texts[rel] = rx.sub(lambda m: m.group(1) + value, text)

    print(f"will update {len(changes)} value(s):" if changes else "nothing to update")
    print("\n".join(changes))

    # Remaining mentions of the old version: listed for a human, never rewritten.
    word = re.compile(r"(?<![\w.])" + re.escape(old) + r"(?![\w.]*\d)")
    leftovers = []
    for rel in tracked_files():
        if rel.startswith(REVIEW_SKIP):
            continue
        try:
            text = texts.get(rel) or (REPO / rel).read_text()
        except (UnicodeDecodeError, FileNotFoundError, IsADirectoryError):
            continue
        for i, line in enumerate(text.splitlines(), 1):
            if word.search(line):
                leftovers.append(f"  {rel}:{i}: {line.strip()[:140]}")
    if leftovers and old != new:
        print(f"\nleft as is ({len(leftovers)} mention(s) of {old}, usually history; review):")
        print("\n".join(leftovers))

    if args.dry_run or not changes:
        return 0
    for rel, text in texts.items():
        if (REPO / rel).read_text() != text:
            (REPO / rel).write_text(text)
    if old == new:
        print(f"\nsynced every target to {new}. Review with git diff.")
    else:
        print(f"\nbumped {old} -> {new}. Next: add a '## v{new}' section to CHANGELOG.md, then review with git diff.")
    return check()


if __name__ == "__main__":
    sys.exit(main())
