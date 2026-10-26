<h1 align="center">
  <img src="res/img/carquet_banner.png" alt="Carquet" width="840" />
</h1>

<p align="center">
  A fast, pure C library for reading and writing Apache Parquet files.
</p>

<p align="center">
  <a href="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml"><img src="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml/badge.svg" alt="Build" /></a>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue" alt="Platform" />
  <img src="https://img.shields.io/badge/C-C11-blue" alt="C Standard" />
  <a href="https://github.com/Vitruves/carquet/releases"><img src="https://img.shields.io/github/v/release/Vitruves/carquet?color=blue" alt="Release" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="License" /></a>
  <a href="https://github.com/Vitruves/carquet/stargazers"><img src="https://img.shields.io/github/stars/Vitruves/carquet?style=flat&color=blue" alt="Stars" /></a>
  <a href="https://github.com/Vitruves/carquet/graphs/traffic"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/Vitruves/carquet/traffic/clones.json" alt="Clones" /></a>
  <a href="https://pepy.tech/projects/carquet"><img src="https://img.shields.io/pepy/dt/carquet?color=blue&label=pypi%20downloads" alt="PyPI downloads" /></a>
  <br/>
  <img src="https://img.shields.io/badge/SIMD-SSE4.2%20%C2%B7%20AVX2%20%C2%B7%20AVX--512%20%C2%B7%20NEON%20%C2%B7%20SVE-555" alt="SIMD" />
  <img src="https://img.shields.io/badge/deps-zstd%20%C2%B7%20zlib%20%C2%B7%20lz4%20%C2%B7%20brotli-555" alt="Dependencies" />
</p>

<p align="center">
  <a href="#highlights">Highlights</a> &nbsp;·&nbsp;
  <a href="#installation">Installation</a> &nbsp;·&nbsp;
  <a href="#c-api">C API</a> &nbsp;·&nbsp;
  <a href="#cli-tool">CLI Tool</a> &nbsp;·&nbsp;
  <a href="#benchmarks">Benchmarks</a> &nbsp;·&nbsp;
  <a href="#parquet-feature-support">Feature Support</a> &nbsp;·&nbsp;
  <a href="#interoperability">Interoperability</a> &nbsp;·&nbsp;
  <a href="#ready-for-ai-coding-agents">AI agents</a> &nbsp;·&nbsp;
  <a href="docs/README.md">Manual</a> &nbsp;·&nbsp;
  <a href="CHANGELOG.md">Changelog</a>
</p>

<p align="center">
  <b>English</b> &nbsp;·&nbsp;
  <a href="README.zh-CN.md">简体中文</a> &nbsp;·&nbsp;
  <a href="README.ja.md">日本語</a> &nbsp;·&nbsp;
  <a href="README.ko.md">한국어</a> &nbsp;·&nbsp;
  <a href="README.de.md">Deutsch</a> &nbsp;·&nbsp;
  <a href="README.fr.md">Français</a> &nbsp;·&nbsp;
  <a href="README.pt-BR.md">Português (Brasil)</a>
</p>

<br/>

---

<!-- ────────────────────────────  HIGHLIGHTS  ──────────────────────────── -->

## Highlights

**Small and self-contained**

- **Pure C11** with four external dependencies (zstd, zlib, lz4, brotli), all auto-fetched by CMake when missing; Brotli can be switched off
- **~500 KB library** (stripped Release build, codecs linked separately) vs ~50 MB+ for Arrow
- **Built-in CLI** for file inspection (`schema`, `info`, `head`, `tail`, `stat`, …) and C code generation (`codegen`)

**Fast**

- **39x faster reads** than Arrow C++ on uncompressed data (mmap zero-copy), **100x faster** than PyArrow
- **1.3–2.1x faster compressed reads** than Arrow C++ on the same file (cross-read benchmark, Apple M3; 1.6–2.6x on x86)
- **Writes 1.8–5.7x faster** than Arrow C++ on Apple M3 (1.0–1.6x on the older x86 run)
- Reads 10M uncompressed rows in **0.36 ms** (mmap zero-copy on Apple M3)
- SIMD-optimized (SSE4.2, AVX2, AVX-512, NEON, SVE) with runtime detection and scalar fallbacks

**Complete and compatible**

- Broad Parquet coverage: all physical types and encodings, Snappy / GZIP / LZ4 / ZSTD / Brotli, nested schemas, bloom filters, page indexes ([full matrix](#parquet-feature-support))
- PyArrow, DuckDB, Spark compatible out of the box
- [Arrow C Data Interface](docs/reading.md#export-to-arrow-c-data-interface) bridge (`carquet_arrow_export_*` / `carquet_arrow_import_*` / `carquet_writer_write_arrow` / `carquet_reader_read_arrow`) for zero-dependency, copy-light interchange with the Arrow ecosystem, with **arbitrary-depth nested read/write** (struct / list / map at any depth)

<a id="ready-for-ai-coding-agents"></a>

**Ready for AI coding agents**

- [`.agents/skills/carquet/`](.agents/skills/carquet/) is an [Agent Skill](https://agentskills.io) (Claude Code, Codex CLI) that teaches an agent to integrate carquet into your project: build wiring, which API for which task, complete recipes, the level / null-bitmap / lifetime rules, type mappings and the usual pitfalls. Copy it into your repository and the agent uses it on its own:

  ```bash
  cp -r path/to/carquet/.agents/skills/carquet .claude/skills/   # Claude Code
  cp -r path/to/carquet/.agents/skills/carquet .agents/skills/   # Codex CLI
  ```

<br/>

### Thirty seconds of carquet

```c
#include <carquet/carquet.h>

carquet_error_t err = CARQUET_ERROR_INIT;
carquet_schema_t* schema = carquet_schema_create(&err);
carquet_schema_add_column(schema, "id", CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

carquet_writer_t* w = carquet_writer_create("ids.parquet", schema, NULL, &err);
int64_t ids[] = {1, 2, 3};
carquet_writer_write_batch(w, 0, ids, 3, NULL, NULL);
carquet_writer_close(w);

carquet_reader_t* r = carquet_reader_open("ids.parquet", NULL, &err);
printf("%lld rows\n", (long long)carquet_reader_num_rows(r));   /* 3 rows */
carquet_reader_close(r);
carquet_schema_free(schema);
```

Full write and read loops, with batches and nullable columns, are in [C API](#c-api).

<br/>

### Where carquet fits

carquet is the small, embeddable option: a ~500 KB C11 library with four auto-fetched dependencies that builds in under a minute, where Arrow C++ weighs 50 MB+ and a PyArrow wheel about 100 MB. On the 10M-row benchmark below it writes compressed files 1.8–5.7x faster than Arrow C++ and reads them 1.3–2.1x faster from the very same file; uncompressed data is served straight from the mmap in 0.36 ms. Arrow stays the right choice when you need its ecosystem (compute kernels, datasets, Flight, encryption); carquet is for the C and C++ programs, embedded systems and language runtimes that only need to read and write Parquet quickly and correctly.

<br/>

---

<!-- ────────────────────────────  INSTALLATION  ──────────────────────────── -->

## Installation

Carquet builds with CMake (primary) or [xmake](https://xmake.io), on Linux, macOS, and Windows.

### Requirements

- C11 compiler (GCC 4.9+, Clang 3.4+, MSVC 2015+)
- CMake 3.16+ (or [xmake](https://xmake.io), see [Building with xmake](#building-with-xmake))
- zstd, zlib, lz4, brotli (auto-fetched if missing; Brotli optional, `-DCARQUET_ENABLE_BROTLI=OFF`)
- OpenMP (optional, for parallel column reading and writing)

<br/>

### Quick Start (make)

The `make` wrapper drives an optimized CMake build and a `/usr/local` install:

```bash
git clone https://github.com/Vitruves/carquet.git
cd carquet
make                              # optimized build (run `make help` to list all targets)
sudo make install                 # install to /usr/local (override with PREFIX=/opt/carquet)
```

<br/>

### Build & Install (CMake)

Invoke CMake directly when you need specific build options or a custom prefix:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release   # add options, e.g. -DCARQUET_BUILD_SHARED=ON
cmake --build build -j$(nproc)
sudo cmake --install build --prefix /usr/local
```

Either path installs:

- `libcarquet.a` (or `.so` / `.dylib` with `-DCARQUET_BUILD_SHARED=ON`)
- `include/carquet/` headers
- `carquet` CLI binary
- `carquet.pc` (pkg-config) and CMake package config for `find_package(carquet)`

After installation, link with `-lcarquet`, or resolve flags via `pkg-config --cflags --libs carquet`.

<br/>

### Python (pip)

The `carquet` wheel bundles the library and exposes it through the Arrow C Data Interface, so tables flow to and from PyArrow without copies of the Parquet decoding logic:

```bash
pip install carquet
```

```python
import carquet
table = carquet.read_table("data.parquet")                 # pyarrow.Table
carquet.write_table(table, "out.parquet", compression="zstd")
print(carquet.info("out.parquet").num_rows)
```

The package lives in [`python/`](python/) (scikit-build-core, ctypes, PyArrow optional at runtime); see its [README](python/README.md) for the API and the wheel build.

<br/>

### Build Options

Pass these to the CMake configure step (`-DOPTION=ON|OFF`):

<div align="center">

| Option | Default | Description |
|--------|---------|-------------|
| `CARQUET_BUILD_DEV` | OFF | Build everything (tests, examples, benchmarks, interop) |
| `CARQUET_BUILD_TESTS` | OFF | Build the test suite only |
| `CARQUET_BUILD_CLI` | ON | Build the `carquet` CLI tool |
| `CARQUET_BUILD_SHARED` | OFF | Build a shared library instead of static |
| `CARQUET_NATIVE_ARCH` | OFF | `-march=native` for maximum performance (host-only binary) |
| `CARQUET_ENABLE_SVE` | OFF | ARM SVE (experimental) |
| `CARQUET_ENABLE_BROTLI` | ON | BROTLI codec (libbrotli found or fetched; OFF drops the dependency) |

</div>

All x86 SIMD (SSE, AVX, AVX2, AVX-512) and ARM NEON are auto-detected and enabled by default.

<details>
<summary><b>All build options</b></summary>

<br/>

<div align="center">

| Option | Default | Description |
|--------|---------|-------------|
| `CARQUET_BUILD_EXAMPLES` | OFF | Build example programs |
| `CARQUET_BUILD_BENCHMARKS` | OFF | Build benchmark and profiling programs |
| `CARQUET_BUILD_ARROW_CPP_BENCHMARK` | OFF | Optional Arrow C++ comparison benchmark |
| `CARQUET_BUILD_INTEROP` | OFF | Build interoperability tests |
| `CARQUET_BUILD_FUZZ` | OFF | Build fuzz targets |
| `CARQUET_BUILD_DOCS` | OFF | Add a `docs` target (needs Doxygen) |
| `CARQUET_BUNDLE_DEPS` | OFF | Always fetch and static-bundle zstd/zlib/lz4/brotli (portable release artifacts) |
| `CARQUET_ZSTD_TARGET` / `CARQUET_ZLIB_TARGET` / `CARQUET_LZ4_TARGET` / `CARQUET_BROTLI_TARGET` | *(empty)* | Existing CMake target(s) to link for that library instead of searching for or fetching it |
| `CARQUET_ARROW_CPP_ROOT` | *(empty)* | Path to a custom Arrow C++ install for the comparison benchmark |
| `CARQUET_ENABLE_SSE` | ON | SSE optimizations (x86, auto-detected) |
| `CARQUET_ENABLE_AVX` | ON | AVX optimizations (x86, auto-detected) |
| `CARQUET_ENABLE_AVX2` | ON | AVX2 optimizations (x86, auto-detected) |
| `CARQUET_ENABLE_AVX512` | ON | AVX-512 optimizations (x86, auto-detected) |
| `CARQUET_ENABLE_NEON` | ON | NEON optimizations (ARM, auto-detected) |

</div>

</details>

<details>
<summary><b>Using your own zstd / zlib / lz4</b> (embedding carquet in a larger CMake project)</summary>

<br/>

When carquet is pulled in with `add_subdirectory()` or `FetchContent`, it resolves each compression library in this order:

1. The target named by `CARQUET_ZSTD_TARGET`, `CARQUET_ZLIB_TARGET` or `CARQUET_LZ4_TARGET`, if set.
2. A system copy (skipped when `CARQUET_BUNDLE_DEPS` is ON).
3. A copy it fetches and builds itself. If your project already builds the library under its usual target name (`libzstd_static`, `zlibstatic`, `lz4_static`), carquet reuses that target instead of fetching a second copy.

A project that vendors these libraries and has no system copy therefore needs no configuration. To choose explicitly, for example your own zlib over the system one, or your shared build over your static one, name the target:

```cmake
set(CARQUET_ZLIB_TARGET zlibstatic)
add_subdirectory(carquet)
```

The same can be done without editing anything: `cmake -DCARQUET_ZLIB_TARGET=zlib ...`. The three variables are independent. A named target is linked as is, with no search and no fetch, and a name that is not a target is a configure error rather than a silent fallback. The target can be imported or built in-tree, static or shared, and must carry its include directories as a usage requirement.

</details>

<details>
<summary><b>Development build</b> (tests, examples, benchmarks)</summary>

<br/>

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCARQUET_BUILD_DEV=ON
cmake --build build -j$(nproc)
cd build && ctest --output-on-failure
```

[CONTRIBUTING.md](CONTRIBUTING.md) describes the per-pull-request gate (tests, sanitizers, fuzzing, interop, benchmarks).

</details>

<details>
<summary><b>API documentation</b> (Doxygen)</summary>

<br/>

The public headers are annotated with [Doxygen](https://www.doxygen.nl) comments. Generating the HTML reference needs `doxygen` installed (plus optional `graphviz` for include/dependency diagrams); neither is required for a normal build.

```bash
# CMake — enable the target at configure time, then build it
cmake -B build -DCARQUET_BUILD_DOCS=ON
cmake --build build --target docs

# xmake — a standalone task, no configure flag needed
xmake docs
```

Both write to `build/docs/html/index.html`. If Doxygen is not installed the target/task simply reports that and does nothing.

</details>

<details>
<summary><b>AI coding agents</b> (Claude Code, Codex: a ready-made skill)</summary>

<br/>

[`.agents/skills/carquet/`](.agents/skills/carquet/) is an [Agent Skill](https://agentskills.io) that teaches a coding agent how to integrate carquet into your project: build wiring (CMake, FetchContent, pkg-config, xmake), which API to use for what, complete recipes, the level / null-bitmap / lifetime rules, type mappings and the usual pitfalls with their fixes. Copy or link it into your own repository and the agent picks it up automatically:

```bash
# Claude Code
mkdir -p .claude/skills && cp -r path/to/carquet/.agents/skills/carquet .claude/skills/
# Codex CLI
mkdir -p .agents/skills && cp -r path/to/carquet/.agents/skills/carquet .agents/skills/
# vendored carquet: a symlink keeps it in sync
ln -s ../../third_party/carquet/.agents/skills/carquet .claude/skills/carquet
```

The same files also work as a personal skill (`~/.claude/skills/`, `~/.codex/skills/`).

</details>

<br/>

### Building with xmake

An [xmake](https://xmake.io) build (`xmake.lua`) is provided as an alternative to CMake, with the same options and defaults; zstd/zlib/lz4 are linked statically so binaries are self-contained.

```bash
xmake                    # build the static library + `carquet` CLI (release)
xmake f --dev=y && xmake # add tests, examples, benchmarks, interop
xmake test               # run the test suite
```

<details>
<summary><b>Configure options</b> (<code>xmake f --option=y|n</code>)</summary>

<br/>

Options mirror the CMake ones, without the `CARQUET_` prefix and lower-case:

<div align="center">

| xmake option | Default | Description |
|--------------|---------|-------------|
| `--dev` | n | Build tests, examples, benchmarks and interop |
| `--tests` / `--examples` / `--benchmarks` / `--interop` | n | Build one group individually |
| `--cli` | y | Build the `carquet` CLI tool |
| `--shared` | n | Build a shared library instead of static |
| `--openmp` | y | OpenMP parallel column reading (auto-disabled if unavailable) |
| `--native_arch` | n | `-march=native` for maximum performance (host-only binary) |
| `--sse` / `--avx` / `--avx2` / `--avx512` / `--neon` | y | SIMD instruction sets (auto-detected) |
| `--sve` | n | ARM SVE (experimental) |
| `--brotli` | y | BROTLI codec (libbrotli from the xmake repo) |
| `--fuzz` | n | Build fuzz targets (use `--toolchain=clang`) |

</div>

```bash
xmake f -m release --shared=y            # shared library
xmake f --dev=y --avx512=n && xmake      # dev build, AVX-512 disabled
```

</details>

<br/>

---

<!-- ────────────────────────────  C API  ──────────────────────────── -->

## C API

A Write and a Read example below; the [manual in `docs/`](docs/README.md) covers everything else.

### Write a Parquet File

```c
#include <carquet/carquet.h>

int main(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    // Define schema
    carquet_schema_t* schema = carquet_schema_create(&err);
    carquet_schema_add_column(schema, "id",    CARQUET_PHYSICAL_INT64,  NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    // Configure writer
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;

    // Write
    carquet_writer_t* w = carquet_writer_create("output.parquet", schema, &opts, &err);

    int64_t ids[]    = {1, 2, 3, 4, 5};
    double values[]  = {1.1, 2.2, 3.3, 4.4, 5.5};
    carquet_writer_write_batch(w, 0, ids, 5, NULL, NULL);
    carquet_writer_write_batch(w, 1, values, 5, NULL, NULL);
    carquet_writer_close(w);

    carquet_schema_free(schema);
    return 0;
}
```

<br/>

### Read a Parquet File

```c
#include <carquet/carquet.h>
#include <stdio.h>

int main(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    // Open with mmap for best read performance
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = true;

    carquet_reader_t* r = carquet_reader_open("output.parquet", &opts, &err);
    if (!r) { printf("Error: %s\n", err.message); return 1; }

    printf("Rows: %lld, Columns: %d\n",
           (long long)carquet_reader_num_rows(r),
           carquet_reader_num_columns(r));

    // Batch reader for efficient iteration
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 65536;

    carquet_batch_reader_t* br = carquet_batch_reader_create(r, &cfg, &err);
    carquet_row_batch_t* batch = NULL;

    while (carquet_batch_reader_next(br, &batch) == CARQUET_OK && batch) {
        const void* data;
        const uint8_t* nulls;
        int64_t n;
        carquet_row_batch_column(batch, 0, &data, &nulls, &n);
        const int64_t* ids = (const int64_t*)data;
        // process ids[0..n-1] ...
        carquet_row_batch_free(batch);
        batch = NULL;
    }

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    return 0;
}
```

<br/>

### More Recipes

Everything beyond flat read/write lives in the manual; each page links to a runnable example:

<div align="center">

| You want to… | See |
|---|---|
| Nullable columns, row groups, buffer output | [`docs/writing.md`](docs/writing.md) |
| Lists, maps, groups, definition/repetition levels | [`docs/nested-data.md`](docs/nested-data.md) |
| Column projection, statistics, metadata inspection | [`docs/reading.md`](docs/reading.md) |
| Predicate pushdown, page-level filtering | [`docs/reading.md`](docs/reading.md) |
| Append row groups to an existing file | [`docs/writing.md`](docs/writing.md#append-to-an-existing-file) |
| Compression, custom codecs, writer tuning | [`docs/writing.md`](docs/writing.md), [`docs/performance.md`](docs/performance.md) |
| mmap, zero-copy, prebuffering, I/O coalescing | [`docs/performance.md`](docs/performance.md) |
| Error codes and recovery hints | [`docs/error-handling.md`](docs/error-handling.md) |

</div>

<br/>

### API Reference

The full API is in [`include/carquet/carquet.h`](include/carquet/carquet.h). Key types:

<div align="center">

| Type | Purpose |
|------|---------|
| `carquet_reader_t` | File reader (open from path, `FILE*`, or memory buffer) |
| `carquet_writer_t` | File writer |
| `carquet_batch_reader_t` | High-level batch iteration |
| `carquet_schema_t` | Schema definition and introspection |
| `carquet_error_t` | Rich error info (code, message, source location, recovery hint) |

</div>

The [manual](docs/README.md) explains which surface to use when; the source layout and architecture are documented in [CONTRIBUTING.md](CONTRIBUTING.md). A browsable HTML reference can be generated from the headers with Doxygen, see the *API documentation* fold under [Build Options](#build-options).

<details>
<summary><b>Example application: MetricStore</b> (<code>mocklib/</code>)</summary>

<br/>

[`mocklib/`](mocklib/) is **MetricStore**, a complete, self-contained example application built on top of carquet, useful both as a reference for real-world usage and as an end-to-end integration test of the public API.

It models a time-series telemetry store (ingest events, then introspect and query them) and, in doing so, exercises **117 of carquet's 143 public functions (~82 %)**:

- **Write side**: schema construction with logical types, per-column encoding / compression / bloom tuning, page indexes and statistics
- **Read side**: column projection, predicate pushdown, page-level filtering, bloom membership, the Arrow C Data Interface bridge, nested `LIST` reconstruction, buffer/append I/O and metadata introspection

It links carquet the way any downstream project would (`find_package(carquet)` / `add_subdirectory`), ships a self-checking round-trip test, and its `write_sample` binary emits a file for external inspection. See [`mocklib/README.md`](mocklib/README.md).

</details>

<br/>

---

<!-- ────────────────────────────  CLI TOOL  ──────────────────────────── -->

## CLI Tool

Carquet ships with a command-line tool for inspecting Parquet files and generating C reader code. It is built and installed by default alongside the library.

```
Commands:
  schema     Print file schema
  info       Print detailed file metadata
  head       Print first N rows
  tail       Print last N rows
  cat        Print rows with slicing/column/row filtering
  count      Print total row count
  columns    List column names (one per line)
  stat       Print column statistics
  validate   Verify file integrity
  sample     Print N random rows
  export     Write rows to stdout as CSV
  codegen    Generate C reader code
```

Typical usage:

```bash
carquet schema data.parquet
carquet head -n 20 data.parquet
carquet stat data.parquet
carquet validate data.parquet
```

<br/>

### Row Filtering

`cat`, `count`, `head`, `tail`, `sample`, and `export` accept `-p / --filter EXPR` to push a row predicate down to the page level: only pages whose column-index min/max can match the predicate are decompressed.

```bash
carquet cat -p "price > 100 AND status = 'active'" data.parquet
carquet count --filter "id >= 1000" data.parquet
carquet export --filter "ts IS NOT NULL" -c id,ts data.parquet
```

The grammar is `column OP value [AND column OP value]...` with `OP` ∈ {`=`, `==`, `!=`, `<>`, `<`, `<=`, `>`, `>=`}, plus `column IS NULL` / `column IS NOT NULL`.

> Filtering requires the file to have a page index (`write_page_index = true`).

<br/>

### Code Generation

Generate a complete, compilable C reader from any Parquet file's schema:

```bash
carquet codegen -f data.parquet -o reader.c
# Generated: reader.c
# Compile:   clang -o reader reader.c -I.../include -L.../build -lcarquet ...

./reader                    # reads data.parquet (embedded as default)
./reader other.parquet      # override with a different file
```

<details>
<summary><b>codegen options</b></summary>

<br/>

<div align="center">

| Flag | Description |
|------|-------------|
| `-f`, `--file FILE` | Parquet file to inspect the schema from |
| `-o`, `--output FILE` | Output source file (default: stdout) |
| `--mmap` | Use memory-mapped I/O in the generated code |
| `--skeleton` | Generate an empty `process_batch` for custom logic |
| `-c`, `--columns COLS` | Comma-separated column filter |
| `-b`, `--batch-size N` | Batch size (default: 1024) |

</div>

</details>

<br/>

---

<!-- ────────────────────────────  BENCHMARKS  ──────────────────────────── -->

## Benchmarks

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="res/img/bench_10m_m3_dark.png">
    <img src="res/img/bench_10m_m3_light.png" alt="Write and read time for 10 million rows on Apple M3: carquet vs Arrow C++ vs PyArrow, per codec" width="900">
  </picture>
</p>

At 10M rows (the most representative size); higher ratio = Carquet faster.

ARM (Apple M3): Carquet 0.7.2 vs Arrow C++ 24.0.0, measured 2026-10-20. x86 (Xeon D-1531): Carquet 0.4.4 vs Arrow C++ 23.0.1; that run predates the writer parallelism of 0.7.1 and has not been re-measured since.

<div align="center">

| | x86 (Xeon D-1531) | | ARM (Apple M3) | |
|---|---|---|---|---|
| **Codec** | **Write** | **Read** | **Write** | **Read** |
| snappy | **1.55x** | **1.25x** | **4.63x** | **1.66x** |
| zstd | **1.31x** | **1.04x** | **5.66x** | **1.36x** |
| lz4 | **1.02x** | 0.83x | **4.47x** | **1.76x** |
| none | **1.13x** | **40.6x**\* | **1.83x** | **38.9x**\* |

</div>

\* Uncompressed reads use mmap zero-copy, see the methodology note.

Compressed reads involve full decompression and decoding of every value, no shortcuts, and both libraries use the same system lz4/zstd shared libraries, so the raw codec speed is identical. The most meaningful comparison is the **same-file cross-read** table (in the full results), where both libraries read the exact same Parquet file: Carquet reads compressed data **1.3–2.1x faster** than Arrow C++ on that apples-to-apples test (x86 run: 1.6–2.6x).

<details>
<summary><b>Benchmark methodology</b></summary>

<br/>

All benchmarks use identical data (deterministic LCG PRNG), identical Parquet settings (no dictionary, BYTE_STREAM_SPLIT for floats, page checksums, mmap reads), trimmed median of 11–51 iterations, with the OS page cache purged between write and read phases and a cooldown between configurations. Schema: 3 columns (INT64, DOUBLE, INT32). Compared against the Arrow C++ low-level Parquet reader (bypassing Arrow Table materialization) and PyArrow.

The **same-file cross-read** benchmark is the fairest comparison: both libraries read the exact same Parquet file (written by one, read by both). This eliminates differences in page sizes, encoding choices, and row group layout.

**Uncompressed reads** marked with \* use Carquet's **mmap zero-copy path**: for PLAIN-encoded, uncompressed, fixed-size, required columns, the batch reader returns pointers directly into the memory-mapped file with no memcpy. Arrow always materializes into its own buffers. **The compressed read numbers are the most representative measure of end-to-end read throughput.**

</details>

<details>
<summary><b>Full x86 results</b> (Intel Xeon D-1531, Linux)</summary>

<br/>

*12 threads @ 2.2GHz, 32GB RAM, Ubuntu 24.04, ZSTD level 1*

#### 10M rows vs Arrow C++

<div align="center">

| Codec | Carquet Write | Arrow C++ Write | W ratio | Carquet Read | Arrow C++ Read | R ratio | Size |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **1557ms** | 1766ms | **1.13x** | **1.25ms** | 50.8ms | **40.6x**\* | 190.7MB |
| snappy | **1002ms** | 1549ms | **1.55x** | **78ms** | 97.8ms | **1.25x** | 125.1MB |
| zstd | **1311ms** | 1714ms | **1.31x** | **76.8ms** | 80.2ms | **1.04x** | 95.3MB |
| lz4 | **1521ms** | 1554ms | **1.02x** | 59.1ms | **49.0ms** | 0.83x | 122.9MB |

</div>

<br/>

#### 1M rows vs Arrow C++

<div align="center">

| Codec | Carquet Write | Arrow C++ Write | W ratio | Carquet Read | Arrow C++ Read | R ratio |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **180ms** | 196ms | **1.09x** | **0.22ms** | 6.2ms | **28x**\* |
| snappy | **141ms** | 148ms | **1.05x** | **8.1ms** | 11.6ms | **1.44x** |
| zstd | **131ms** | 185ms | **1.41x** | 10.3ms | **9.1ms** | 0.88x |
| lz4 | **143ms** | 149ms | **1.04x** | 8.5ms | **6.1ms** | 0.72x |

</div>

<br/>

#### 100K rows vs Arrow C++

<div align="center">

| Codec | Carquet Write | Arrow C++ Write | W ratio | Carquet Read | Arrow C++ Read | R ratio |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **14.1ms** | 18.4ms | **1.30x** | **0.11ms** | 2.18ms | **19.8x**\* |
| snappy | **10.1ms** | 10.6ms | **1.05x** | **1.27ms** | 5.97ms | **4.70x** |
| zstd | **8.7ms** | 14.1ms | **1.62x** | **1.58ms** | 3.88ms | **2.46x** |
| lz4 | **9.6ms** | 11.0ms | **1.14x** | **0.77ms** | 2.78ms | **3.61x** |

</div>

<br/>

#### Same-file cross-read (10M rows)

Both libraries read the **same** Parquet file, the fairest apples-to-apples comparison.

<div align="center">

| Codec | Writer | Carquet Read | Arrow C++ Read | Ratio |
|-------|--------|-------------|----------------|-------|
| none | Carquet | **0.99ms** | 73.6ms | **74x**\* |
| none | Arrow | **7.6ms** | 51.2ms | **6.8x**\* |
| snappy | Carquet | **41.0ms** | 107ms | **2.61x** |
| snappy | Arrow | **43.4ms** | 101ms | **2.33x** |
| zstd | Carquet | **46.1ms** | 88.4ms | **1.92x** |
| zstd | Arrow | **49.1ms** | 79.5ms | **1.62x** |
| lz4 | Carquet | **34.8ms** | 74.8ms | **2.15x** |
| lz4 | Arrow | **27.4ms** | 52.0ms | **1.90x** |

</div>

<br/>

#### 10M rows vs PyArrow

<div align="center">

| Codec | Carquet Write | PyArrow Write | W ratio | Carquet Read | PyArrow Read | R ratio |
|-------|--------------|---------------|---------|-------------|--------------|---------|
| none | **1557ms** | 1806ms | **1.16x** | **1.25ms** | 213ms | **170x**\* |
| snappy | **1002ms** | 1649ms | **1.65x** | **78ms** | 384ms | **4.91x** |
| zstd | **1311ms** | 1796ms | **1.37x** | **76.8ms** | 369ms | **4.81x** |
| lz4 | **1521ms** | 1676ms | **1.10x** | **59.1ms** | 281ms | **4.76x** |

</div>

\* Zero-copy mmap path

</details>

<details>
<summary><b>Full ARM results</b> (Apple M3, macOS)</summary>

<br/>

*Carquet 0.7.2, MacBook Air M3, 16GB RAM, macOS 26.5, Arrow C++ 24.0.0, PyArrow 23.0.1, ZSTD level 1. The 10M-row tables were measured on 2026-10-20; the 1M and 100K tables are the 0.7.1 run of 2026-10-08 (not re-measured).*

#### 10M rows vs Arrow C++

<div align="center">

| Codec | Carquet Write | Arrow C++ Write | W ratio | Carquet Read | Arrow C++ Read | R ratio | Size |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **74.7ms** | 137ms | **1.83x** | **0.36ms** | 14.0ms | **38.9x**\* | 200.0MB |
| snappy | **52.7ms** | 244ms | **4.63x** | **12.4ms** | 20.6ms | **1.66x** | 131.2MB |
| zstd | **60.4ms** | 342ms | **5.66x** | **19.8ms** | 26.9ms | **1.36x** | 99.9MB |
| lz4 | **57.2ms** | 256ms | **4.47x** | **9.22ms** | 16.2ms | **1.76x** | 128.9MB |

</div>

<br/>

#### 1M rows vs Arrow C++

<div align="center">

| Codec | Carquet Write | Arrow C++ Write | W ratio | Carquet Read | Arrow C++ Read | R ratio |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **6.48ms** | 14.82ms | **2.29x** | **0.06ms** | 1.63ms | **27.2x**\* |
| snappy | **14.28ms** | 26.44ms | **1.85x** | **2.14ms** | 2.85ms | **1.33x** |
| zstd | **14.76ms** | 36.17ms | **2.45x** | **2.70ms** | 3.46ms | **1.28x** |
| lz4 | **12.30ms** | 26.10ms | **2.12x** | **1.04ms** | 1.87ms | **1.80x** |

</div>

<br/>

#### 100K rows vs Arrow C++

<div align="center">

| Codec | Carquet Write | Arrow C++ Write | W ratio | Carquet Read | Arrow C++ Read | R ratio |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **1.07ms** | 1.79ms | **1.67x** | **0.02ms** | 0.25ms | **12.5x**\* |
| snappy | **1.60ms** | 2.63ms | **1.64x** | **0.34ms** | 0.91ms | **2.68x** |
| zstd | **1.73ms** | 3.63ms | **2.10x** | **0.63ms** | 1.26ms | **2.00x** |
| lz4 | **1.52ms** | 2.62ms | **1.72x** | **0.25ms** | 0.55ms | **2.20x** |

</div>

<br/>

#### Same-file cross-read (10M rows)

Both libraries read the **same** Parquet file, the fairest apples-to-apples comparison.

<div align="center">

| Codec | Writer | Carquet Read | Arrow C++ Read | Ratio |
|-------|--------|-------------|----------------|-------|
| none | Carquet | **0.31ms** | 15.7ms | **50.7x**\* |
| none | Arrow | **0.93ms** | 13.7ms | **14.8x**\* |
| snappy | Carquet | **12.6ms** | 22.0ms | **1.74x** |
| snappy | Arrow | **14.3ms** | 21.8ms | **1.52x** |
| zstd | Carquet | **21.3ms** | 28.0ms | **1.31x** |
| zstd | Arrow | **21.1ms** | 30.4ms | **1.44x** |
| lz4 | Carquet | **8.48ms** | 17.7ms | **2.09x** |
| lz4 | Arrow | **7.90ms** | 16.3ms | **2.06x** |

</div>

<br/>

#### 10M rows vs PyArrow

<div align="center">

| Codec | Carquet Write | PyArrow Write | W ratio | Carquet Read | PyArrow Read | R ratio |
|-------|--------------|---------------|---------|-------------|--------------|---------|
| none | **74.7ms** | 202ms | **2.71x** | **0.36ms** | 36.2ms | **100.5x**\* |
| snappy | **52.7ms** | 312ms | **5.93x** | **12.4ms** | 46.6ms | **3.75x** |
| zstd | **60.4ms** | 410ms | **6.78x** | **19.8ms** | 75.6ms | **3.82x** |
| lz4 | **57.2ms** | 321ms | **5.62x** | **9.22ms** | 40.0ms | **4.33x** |

</div>

<br/>

#### 1M rows vs PyArrow

<div align="center">

| Codec | Carquet Write | PyArrow Write | W ratio | Carquet Read | PyArrow Read | R ratio |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **6.48ms** | 20.17ms | **3.11x** | **0.06ms** | 2.79ms | **46.5x**\* |
| snappy | **14.28ms** | 31.54ms | **2.21x** | **2.14ms** | 4.01ms | **1.87x** |
| zstd | **14.76ms** | 42.51ms | **2.88x** | **2.70ms** | 4.42ms | **1.64x** |
| lz4 | **12.30ms** | 32.03ms | **2.60x** | **1.04ms** | 3.22ms | **3.10x** |

</div>

<br/>

#### 100K rows vs PyArrow

<div align="center">

| Codec | Carquet Write | PyArrow Write | W ratio | Carquet Read | PyArrow Read | R ratio |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **1.07ms** | 2.12ms | **1.98x** | **0.02ms** | 0.24ms | **12.0x**\* |
| snappy | **1.60ms** | 3.17ms | **1.98x** | **0.34ms** | 0.60ms | **1.76x** |
| zstd | **1.73ms** | 4.16ms | **2.40x** | **0.63ms** | 0.81ms | **1.29x** |
| lz4 | **1.52ms** | 3.21ms | **2.11x** | **0.25ms** | 0.42ms | **1.68x** |

</div>

\* Zero-copy mmap path

</details>

<details>
<summary><b>Running the benchmarks yourself</b></summary>

<br/>

```bash
# Build with max optimizations
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON -DCARQUET_BUILD_DEV=ON
cmake --build build -j$(nproc)

cd build
./benchmark_carquet                     # Carquet standalone
python3 ../benchmark/run_benchmark.py   # Full comparison (+ PyArrow, + Arrow C++)

# Skip 100M-row (xlarge) configs — they write ~2GB files per codec
# and can take 30+ minutes depending on hardware
python3 ../benchmark/run_benchmark.py --skip-xlarge

# Override ZSTD level (default: 1)
CARQUET_BENCH_ZSTD_LEVEL=3 python3 ../benchmark/run_benchmark.py
```

To include the Arrow C++ comparison, enable its target at configure time:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON \
  -DCARQUET_BUILD_BENCHMARKS=ON \
  -DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON
cmake --build build -j$(nproc)

# Or point at a custom Arrow install
cmake -B build ... -DCARQUET_ARROW_CPP_ROOT=/path/to/arrow-prefix
```

The Arrow C++ benchmark uses the low-level `parquet::ParquetFileReader` API (bypassing Arrow Table materialization overhead) with parallel row group readers. The **same-file cross-read** mode has both libraries read the exact same Parquet file, eliminating differences in page sizes, encoding, and row group layout. Both benchmarks use identical data, row group sizing, no dictionary, page checksums, mmap reads, BYTE_STREAM_SPLIT for floats.

</details>

<br/>

---

<!-- ────────────────────────────  FEATURE SUPPORT  ──────────────────────────── -->

## Parquet Feature Support

What the library implements today, by area.

### Format

<div align="center">

| Feature | Status |
|---------|--------|
| Physical types | ✓  All 8 (BOOLEAN through FIXED_LEN_BYTE_ARRAY) |
| Logical types | ✓  STRING, DATE, TIME, TIMESTAMP, DECIMAL, UUID, JSON, INTERVAL, FLOAT16, VARIANT, GEOMETRY, GEOGRAPHY |
| Encodings | ✓  PLAIN, RLE, DICTIONARY, DELTA_BINARY_PACKED, DELTA_LENGTH_BYTE_ARRAY, DELTA_BYTE_ARRAY, BYTE_STREAM_SPLIT (read + write) |
| Data Page versions | ✓  V1 (default) and V2 (read + write) |
| Compression | ✓  UNCOMPRESSED, SNAPPY, GZIP, LZ4 (Hadoop-framed, codec 5), LZ4_RAW (codec 7), ZSTD, BROTLI (`CARQUET_ENABLE_BROTLI`, default ON); LZO via `carquet_register_codec()` |
| Nested schemas | ✓  Groups, lists, maps with definition/repetition levels; single-level LIST/MAP auto-shredding on write (`carquet_writer_write_list_column`) and List reconstruction on read; nested `ARROW:schema` emission |
| Encryption | ✗ Not supported; an encrypted-footer file (`PARE`) is reported as `CARQUET_ERROR_NOT_IMPLEMENTED` with an explicit message |

</div>

<details>
<summary><b>Indexing &amp; pruning</b></summary>

<br/>

<div align="center">

| Feature | Status |
|---------|--------|
| Bloom filters | ✓  Read, write, and query (`carquet_bloom_filter_check_*`); consulted automatically for row-group pruning |
| Page indexes | ✓  Column index + offset index (read + write + per-page stats access) |
| Statistics | ✓  Min/max/null count per column chunk; exact `distinct_count` for dictionary columns; Parquet 2.9 SizeStatistics (unencoded byte-array bytes + level histograms) |
| Predicate pushdown | ✓  Automatic row-group pruning via statistics + bloom filters (no callback needed); page-level filtering via column index (`carquet_batch_reader_set_page_filter`) |
| Geospatial statistics | ✓  Per-column bounding box and geometry types for GEOMETRY/GEOGRAPHY (`carquet_reader_geospatial_statistics`) |
| Column projection | ✓  Read only selected columns |
| Row window | ✓  Read an `[offset, limit)` slice of rows, skipping pages entirely (`carquet_batch_reader_set_row_range`) |

</div>

</details>

<details>
<summary><b>I/O &amp; extensibility</b></summary>

<br/>

<div align="center">

| Feature | Status |
|---------|--------|
| Append | ✓  Add row groups to an existing file (`carquet_writer_open_append`) |
| Custom codecs | ✓  Register a custom compress/decompress impl per codec slot (`carquet_register_codec`) |
| Key-value metadata | ✓  Read and write arbitrary footer metadata |
| Per-field metadata | ✓  Arrow `Field.custom_metadata` (variable labels/descriptions) via `ARROW:schema` (read + write) |
| Per-column options | ✓  Per-column encoding, compression, statistics, bloom filter |
| Buffer writer | ✓  Write Parquet to an in-memory buffer |
| CRC32 | ✓  Page-level verification (HW-accelerated on ARM) |
| Memory-mapped I/O | ✓  Zero-copy reads for uncompressed PLAIN data |
| Dictionary-preserving reads | ✓  Return RLE indices + dictionary instead of materialized values (`preserve_dictionaries`, `carquet_row_batch_column_dictionary`) |
| I/O coalescing | ✓  Pre-buffer multi-column reads in a single I/O |
| Speculative footer | ✓  Single-I/O file open for most files |
| OpenMP parallel reads | ✓  When available |
| Parallel writes | ✓  Page-level parallel encode/compress across columns (OpenMP) + background row-group I/O (`async_io`) |

</div>

</details>

<br/>

---

<!-- ────────────────────────────  INTEROPERABILITY  ──────────────────────────── -->

## Interoperability

Carquet files are fully compatible with PyArrow, DuckDB, Spark, and any Parquet reader.

```python
import pyarrow.parquet as pq
table = pq.read_table("carquet_output.parquet")  # just works
```

```sql
-- DuckDB
SELECT * FROM read_parquet('carquet_output.parquet');
```

<details>
<summary><b>Interop test suite</b></summary>

<br/>

Bidirectional round-trips (C writes / Python reads, and back):

```bash
cmake -B build -DCARQUET_BUILD_INTEROP=ON && cmake --build build
python3 interop/run_interop.py
```

Large multi-row-group Snappy files (many row groups, wide `BYTE_ARRAY` columns, page offsets past 2 GiB) have their own generator, which builds the file with PyArrow and checks that carquet decodes every page:

```bash
CARQUET_BIN=build/carquet python3 interop/snappy_multi_rowgroup.py --rows 150000 --row-groups 92
```

</details>

<br/>

---

<!-- ────────────────────────────  CONTRIBUTING & LICENSE  ──────────────────────────── -->

## Contributing, Citing, License

- **Contributing**: [CONTRIBUTING.md](CONTRIBUTING.md) covers setup, code style and the per-pull-request gate (tests, sanitizers, fuzzing, interop, benchmarks). Release notes are in [CHANGELOG.md](CHANGELOG.md).
- **Citing**: a [`CITATION.cff`](CITATION.cff) is included; GitHub's "Cite this repository" button formats it as BibTeX or APA.
- **License**: MIT, see [LICENSE](LICENSE).
