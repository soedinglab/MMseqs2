# Contributing to Carquet

Thanks for contributing. Carquet is a pure-C Parquet library that other tools trust to produce and read correct bytes, so the bar for merging is deliberately high: a change is not "done" until it is tested, fuzzed, sanitizer-clean, documented, and green on every CI platform.

This file is the **routine gate — every pull request must pass all of it**, no matter how small the change. Cutting an actual release has extra steps (version bump, fuzz soak, tagging); those live in **[RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md)** and are the maintainer's job.

If in doubt, do more, not less. A skipped step is a regression waiting to be shipped. When a change genuinely cannot reach one of these (e.g. a fuzz target it does not touch), say so explicitly in the PR description with the reason — silent omission is not allowed.

## Getting set up

```bash
# Debug build with tests (what you develop against)
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCARQUET_BUILD_TESTS=ON
cmake --build build
cd build && ctest --output-on-failure
```

Compilers supported: GCC, Clang, MSVC. See [CLAUDE.md](CLAUDE.md) for the architecture overview, directory layout, and the full set of build options.

## Style

Write code that reads like the code around it — match the surrounding naming, error-handling idioms (`carquet_status_t` returns, `carquet_error_t*` out-params, `CARQUET_RETURN_IF_ERROR`), and comment density. Physical types are the storage format; logical types are the semantic interpretation. New test files must be wired into **both** `CMakeLists.txt` and `xmake.lua`.

---

## The gate — every pull request

### 1. Tests pass

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCARQUET_BUILD_TESTS=ON
cmake --build build
cd build && ctest --output-on-failure
```

All targets must pass (currently 30 test files under `tests/`). A change that touches encoding, reading, writing, metadata, or SIMD must run the full suite, not just the one test you edited.

### 2. Unit tests for your change

Any new public function, type, error code, or behavior needs its own unit test **in the same PR** — not "later". The test lives in the matching `tests/test_*.c` (or a new `tests/test_<feature>.c` wired into both `CMakeLists.txt` and `xmake.lua`). It must cover:

- the happy path,
- at least one error / rejection path (bad input, malformed file, allocation-failure where reachable),
- boundary values (empty, one element, max-values, null-heavy columns).

**If you wrote a wire-format encoder** (Thrift compact-protocol metadata, or a Parquet page/encoding byte stream), the test must include a **byte-level assertion** on the produced bytes — assert the exact serialized output against a known-good hand-verified buffer, not just that it round-trips through your own decoder. Round-trip-only tests hide a mismatched-but-self-consistent encoder that no other reader (Arrow, parquet-mr, DuckDB) will accept. Cross-check the bytes against the Parquet/Thrift spec and, where practical, against what pyarrow writes for the same logical input.

### 3. ASan + UBSan clean on the test suite

The suite must run clean under AddressSanitizer **and** UndefinedBehaviorSanitizer — zero errors, zero leaks.

```bash
cmake -B build-san -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_FLAGS="-fsanitize=address,undefined -fno-sanitize-recover=all -fno-omit-frame-pointer" \
  -DCARQUET_BUILD_TESTS=ON
cmake --build build-san
cd build-san && ctest --output-on-failure
```

`-fno-sanitize-recover=all` makes UBSan abort (non-zero exit) on the first violation instead of merely printing — so ctest actually fails. A leak that is "only in the test harness" is still a leak: fix it or document why it is unreachable in library code.

### 4. Fuzzing — at least 60s per reachable target

Run every fuzz target that your change can reach, **at least 1 minute each**, with sanitizers on (the default). Zero crashes, zero leaks, zero timeouts.

```bash
# Build once, then run every target for 60s each:
python3 fuzz/run_fuzzer.py all --time 60

# Or a single target while iterating:
python3 fuzz/run_fuzzer.py reader --time 60
```

Current targets (`run_fuzzer.py list`): `reader`, `writer`, `compression`, `encodings`, `thrift`, `roundtrip`, `page_filter`, `append`.

**Update the fuzzers when the feature benefits from it.** If your change adds:

- a new encoding / codec / logical type → extend `fuzz_encodings.c` / `fuzz_compression.c` / `fuzz_roundtrip.c` so the new path is actually exercised;
- a new reader entry point or metadata parser → extend `fuzz_reader.c` / `fuzz_thrift.c`;
- a new writer surface (options, append, filters) → extend `fuzz_writer.c` / `fuzz_append.c` / `fuzz_page_filter.c`;
- a genuinely new attack surface with no matching target → add a new `fuzz_<feature>.c` and register it in `fuzz/CMakeLists.txt`, the `fuzzers` list in `xmake.lua`, and `run_fuzzer.py`.

A new decoder that no fuzzer reaches is not "fuzzed clean" — it is unfuzzed. Growing the corpus (`fuzz/corpus/`, `fetch_corpus.sh`) for the new path is encouraged before the timed run. The much longer soak run happens at release time — see [RELEASE_CHECKLIST.md](RELEASE_CHECKLIST.md).

### 5. Docs updated

If behavior, API, or invariants changed, update the docs in the same PR:

- `docs/` — the relevant guide (`reading.md`, `error-handling.md`, etc.);
- `include/carquet/*.h` — the doc comment on every new/changed public symbol;
- **the component README that owns the area you changed** — e.g. `fuzz/README.md` for a new fuzz target, `docs/README.md` for a new guide, `profiling/README.md` for profiling changes;
- **the main `README.md`** only for user-visible changes worth surfacing on the front page (a new feature, a new build option, a changed default). Keep it **concise** — the main README is a landing page, not a manual. Add a short line or a table row and **link out** to the relevant `docs/` guide, component README, or header; do not paste large explanations or API dumps into it;
- `CLAUDE.md` if the change alters the architecture, build options, or directory layout;
- examples under `examples/` if a new capability deserves a runnable demo.

Undocumented public API is treated as incomplete.

### 6. CHANGELOG entry

Add a bullet to the **top unreleased section** of `CHANGELOG.md`, under the right heading (`API`, `New Features`, `Security / Bug Fixes`, `Interoperability`, `Performance`, `Robustness`). Describe what changed and, for a fix, *what was wrong and how it manifested* — match the existing entries' level of detail. Note explicitly whether the change is ABI-compatible and whether default output bytes are unchanged. (The maintainer rolls these bullets into a dated release section at release time — you just add the bullet.)

### 7. Interop still holds

If your change touches the reader, writer, or any encoding/metadata on the wire, run the interop suite and confirm carquet still round-trips against pyarrow (and the other libs the suite drives):

```bash
python3 interop/run_interop.py -v
```

**Add coverage when you add wire surface.** A new encoding, logical type, codec, or writer option must get a case in `interop/roundtrip_writer.c` (and the generator/driver in `run_interop.py` / `generate_test_files.py`) so the new path is actually proven to interoperate. Shipping a wire feature with no interop case is shipping it untested. Any known, deliberate gap (e.g. DuckDB not reading BSS-int) must be an explicitly asserted, documented exception in the harness — not a silent skip.

### 8. CI is green on all platforms

The PR cannot merge until **all** GitHub Actions jobs pass (`.github/workflows/cpp.yml`):

- **`build`** matrix — macOS, Windows (MSVC), Ubuntu x86-64, Ubuntu ARM — each with `BUILD_SHARED_LIBS` both `ON` and `OFF`;
- **`redhat-build`** — Rocky Linux 9 (GCC).

"Green on my machine" is not green. Cross-platform failures (MSVC warnings-as-errors, big-endian assumptions, ARM SIMD, shared-vs-static symbol visibility) are common and must be fixed, not retried. If CI is flaky, fix the flake — do not merge over a red run.

### 9. Benchmark — if the change may affect performance

If your change touches anything on a hot path — an encoding/decoding routine, compression, SIMD dispatch or a kernel, the bit-packer, the reader/writer inner loops, memory layout, or allocation behavior — you must benchmark it and show it is not a regression. "It's probably fine" is not evidence.

Build the benchmarks in **Release** with native arch (the numbers are meaningless in Debug or without `-march=native`), and build the **Arrow C++ comparison** so carquet is measured against a real-world baseline. The benchmark runner looks for the binaries in `./build`, so build into `build`:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release \
  -DCARQUET_BUILD_BENCHMARKS=ON \
  -DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON \
  -DCARQUET_NATIVE_ARCH=ON \
  -DCARQUET_ARROW_CPP_ROOT=/path/to/arrow   # only if Arrow/Parquet isn't found automatically
cmake --build build
```

`CARQUET_BUILD_ARROW_CPP_BENCHMARK` needs a C++ compiler and an Arrow + Parquet install (Homebrew `apache-arrow`, a distro package, or a source build — point `CARQUET_ARROW_CPP_ROOT` at its prefix if CMake can't find it). Then run the comparison — Arrow C++ is used as the baseline automatically once `build/benchmark_arrow_cpp` exists:

```bash
python3 benchmark/run_benchmark.py --quick   # carquet vs Arrow C++ (baseline) vs PyArrow
python3 benchmark/run_benchmark.py           # full sweep (small→xlarge, all codecs) before a release-worthy perf claim
```

Requirements:

- **No regression** versus the pre-change numbers on the paths you touched. Compare against the most recent committed `benchmark/bench_<version>_<platform>_<date>.json` snapshot; if none is comparable, capture a baseline by stashing your change and running first.
- **Report the numbers in the PR** — the relevant rows of the comparison table (carquet vs Arrow C++), not just "no regression". If the change is a speedup, quantify it the way the CHANGELOG entries do (e.g. "~1.8× on Apple Silicon, output byte-for-byte identical").
- A **deliberate** trade-off (slower but more correct/safer) must be called out and justified in the PR and the CHANGELOG.

---

## Quick reference — the per-PR gate

| # | Requirement |
|---|-------------|
| 1 | Tests pass (`ctest`, all targets) |
| 2 | Unit tests for the change (+ byte-level assertion for new encoders) |
| 3 | ASan + UBSan clean on the suite |
| 4 | Fuzz ≥60s per reachable target; fuzzers updated for new features |
| 5 | Docs updated (`docs/`, headers, README/CLAUDE, examples) |
| 6 | CHANGELOG bullet in the unreleased section |
| 7 | Interop validated via `run_interop.py` (+ new cases) if the wire was touched |
| 8 | CI green on every platform in the matrix |
| 9 | Benchmark vs Arrow C++ (Release + native), no regression, numbers in the PR — if a hot path was touched |

## Reporting security issues

If you find a memory-safety bug or a crafted-file crash, please treat it as security-sensitive: include a minimal reproducer (ideally a fuzz corpus input) and describe the out-of-bounds / overflow condition rather than posting it as a routine issue.
