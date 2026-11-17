# Building carquet

Every way to configure and build the library, the `carquet` CLI, the tests and the benchmarks. For the short version (clone, `make`, install), see [Installation](../README.md#installation) in the README.

carquet needs a C11 compiler (GCC 4.9+, Clang 3.4+, MSVC 2015+) and CMake 3.16+ or [xmake](https://xmake.io). zstd, zlib, lz4 and brotli are used from the system when present and fetched otherwise; OpenMP is optional.

## `make` targets

Every target configures `build/` and builds it in parallel; `make help` prints the list. Each one has an [xmake equivalent](#xmake).

| Command | Build type | What you get |
|---|---|---|
| `make` / `make release` | Release | Static library + `carquet` CLI (the default) |
| `make native` | Release | Same, tuned for this CPU with `-march=native` (binary does not run on older CPUs) |
| `make shared` | Release | Shared library (`.so` / `.dylib` / `.dll`) instead of static |
| `make benchmark` | Release | Benchmark programs (`benchmark_carquet`, profiling targets) |
| `make debug` | Debug | Library + CLI with assertions, no optimization |
| `make test` | Debug | Test suite, then runs it with `ctest` |
| `make dev` | Debug | Everything: tests, examples, benchmarks, interop programs |
| `make fuzz` | Debug | libFuzzer targets (use a Clang toolchain) |
| `make install` | Release | Builds `release`, then installs to `PREFIX` |
| `make clean` | | Deletes the build directory |

| Variable | Default | Example |
|---|---|---|
| `PREFIX` | CMake default (`/usr/local` on Unix) | `make install PREFIX=$HOME/.local` |
| `BUILD` | `build` | `make test BUILD=build-debug` (keep several configurations side by side) |
| `JOBS` | number of cores | `make JOBS=4` |
| `CMAKE` | `cmake` | `make CMAKE=/opt/cmake/bin/cmake` |

All targets share `build/` by default, so `make test` after `make benchmark` reconfigures it as Debug. Use `BUILD=` to keep a Release and a Debug tree apart.

## CMake options

For anything the `make` targets don't cover, call CMake directly and pass options with `-DOPTION=ON|OFF`:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_BUILD_SHARED=ON -DCARQUET_ENABLE_BROTLI=OFF
cmake --build build -j
sudo cmake --install build --prefix /usr/local
```

| What to build | Default | |
|---|:---:|---|
| `CARQUET_BUILD_CLI` | ON | `carquet` command-line tool |
| `CARQUET_BUILD_SHARED` | OFF | Shared instead of static library |
| `CARQUET_BUILD_DEV` | OFF | Turns on tests, examples, benchmarks and interop at once |
| `CARQUET_BUILD_TESTS` | = DEV | Test suite (`ctest`) |
| `CARQUET_BUILD_EXAMPLES` | = DEV | Programs in [`examples/`](../examples/) |
| `CARQUET_BUILD_BENCHMARKS` | = DEV | Benchmark and profiling programs |
| `CARQUET_BUILD_ARROW_CPP_BENCHMARK` | OFF | Arrow C++ comparison benchmark (needs Arrow C++) |
| `CARQUET_BUILD_INTEROP` | = DEV | C / Python round-trip programs |
| `CARQUET_BUILD_FUZZ` | OFF | Fuzz targets |
| `CARQUET_BUILD_DOCS` | OFF | `docs` target, Doxygen HTML into `build/docs/html/` |

| Dependencies | Default | |
|---|:---:|---|
| `CARQUET_ENABLE_BROTLI` | ON | BROTLI codec; OFF drops the libbrotli dependency |
| `CARQUET_BUNDLE_DEPS` | OFF | Always fetch and statically bundle zstd / zlib / lz4 / brotli (portable artifacts, Python wheel) |
| `CARQUET_ZSTD_TARGET`, `_ZLIB_`, `_LZ4_`, `_BROTLI_TARGET` | empty | Link an existing CMake target instead of searching or fetching (see below) |
| `CARQUET_ARROW_CPP_ROOT` | empty | Arrow C++ install prefix for the comparison benchmark |
| `CMAKE_DISABLE_FIND_PACKAGE_OpenMP` | OFF | ON builds without OpenMP (no `libgomp` / `vcomp140.dll` at runtime) |

| CPU | Default | |
|---|:---:|---|
| `CARQUET_NATIVE_ARCH` | OFF | `-march=native`; fastest, host-only binary |
| `CARQUET_ENABLE_SSE`, `_AVX`, `_AVX2`, `_AVX512` | ON | x86 SIMD kernels, chosen at runtime from what the CPU supports |
| `CARQUET_ENABLE_NEON` | ON | ARM NEON kernels |
| `CARQUET_ENABLE_SVE` | OFF | ARM SVE kernels (experimental) |
| `CARQUET_MSVC_HARDENING` | ON | MSVC: `/guard:cf` and `/CETCOMPAT` |

SIMD kernels are compiled in and selected at runtime with scalar fallbacks, so a default build runs on any CPU of its architecture. Only `CARQUET_NATIVE_ARCH` ties the binary to the build machine.

## xmake

`xmake.lua` is a full alternative to CMake with the same options and defaults. Configure once with `xmake f`, then `xmake` builds; the configuration is remembered until the next `xmake f`.

| `make` target | xmake equivalent |
|---|---|
| `make` | `xmake` |
| `make native` | `xmake f --native_arch=y && xmake` |
| `make shared` | `xmake f --shared=y && xmake` |
| `make benchmark` | `xmake f --benchmarks=y && xmake` |
| `make debug` | `xmake f -m debug && xmake` |
| `make test` | `xmake f -m debug --tests=y && xmake && xmake test` |
| `make dev` | `xmake f -m debug --dev=y && xmake` |
| `make fuzz` | `xmake f -m debug --fuzz=y --toolchain=clang && xmake` |
| `make clean` | `xmake clean` (or `xmake f -c` to also reset the configuration) |
| `cmake --build build --target docs` | `xmake docs` (no configure flag needed) |

Options mirror the CMake ones without the `CARQUET_` prefix, in lower case, as `--option=y|n`. Combine them freely: `xmake f -m release --shared=y --brotli=n --avx512=n`.

| xmake option | Default | Description |
|--------------|---------|-------------|
| `--dev` | n | Build tests, examples, benchmarks and interop |
| `--tests` / `--examples` / `--benchmarks` / `--interop` | n | Build one group individually |
| `--arrow_cpp_benchmark` | n | Arrow C++ comparison benchmark (uses Arrow from the xmake repo if available) |
| `--cli` | y | Build the `carquet` CLI tool |
| `--shared` | n | Build a shared library instead of static |
| `--openmp` | y | OpenMP parallel column reading (auto-disabled if unavailable) |
| `--native_arch` | n | `-march=native` for maximum performance (host-only binary) |
| `--sse` / `--avx` / `--avx2` / `--avx512` / `--neon` | y | SIMD instruction sets (auto-detected) |
| `--sve` | n | ARM SVE (experimental) |
| `--brotli` | y | BROTLI codec (libbrotli from the xmake repo) |
| `--fuzz` | n | Build fuzz targets (use `--toolchain=clang`); then `--fuzz_engine=libFuzzer\|AFL`, `--fuzz_sanitizers=y\|n` |

## Windows (MSVC, clang-cl, MinGW-w64)

MSVC (x64 and Win32), clang-cl and MinGW-w64 GCC are all supported, static or shared. `CMakePresets.json` has `msvc`, `msvc-shared`, `clang-cl` and `mingw` presets:

```powershell
# from a Developer PowerShell for VS (or with MinGW-w64 gcc + ninja on PATH for "mingw")
cmake --preset msvc
cmake --build --preset msvc
ctest --preset msvc
```

- When configuring by hand, set the C++ compiler as well as the C one (for example `CC=clang-cl CXX=clang-cl`): the fetched zstd project enables C++, and CMake otherwise picks a mismatched driver.
- With `-DCARQUET_BUILD_SHARED=ON` you get `carquet.dll` plus the `carquet.lib` import library (`libcarquet.dll` / `libcarquet.dll.a` with MinGW). Consumers must define `CARQUET_BUILD_SHARED`; `find_package(carquet)` and `carquet.pc` add it for you. Put the DLL next to your executable or on `PATH`.
- Release memory that carquet hands you, such as `carquet_writer_get_buffer()` output, with `carquet_free()` and not `free()`. The DLL and your program may use different C runtime heaps (Debug vs Release CRT, `/MT`, MinGW vs MSVC).
- File paths may be UTF-8 (non-ASCII names work). The `carquet` CLI accepts Unicode arguments and prints UTF-8 to the console.
- MinGW builds with OpenMP need `libgomp-1.dll` and `libwinpthread-1.dll` at runtime. MSVC builds need `vcomp140.dll` (from the Visual C++ Redistributable). Pass `-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=ON` for a build without that dependency.
- Windows Defender's machine-learning heuristic may quarantine `test_encoding_roundtrip.exe` and `test_fuzz_regressions.exe` as `Trojan:Win64/CobaltStrike.RDH!MTB` (a false positive on their embedded binary test data; ctest then reports them as "Not Run"). Exclude the build directory from scanning to run the full suite.

## Using your own zstd / zlib / lz4 / brotli

For embedding carquet in a larger CMake project.

When carquet is pulled in with `add_subdirectory()` or `FetchContent`, it resolves each compression library in this order:

1. The target named by `CARQUET_ZSTD_TARGET`, `CARQUET_ZLIB_TARGET`, `CARQUET_LZ4_TARGET` or `CARQUET_BROTLI_TARGET`, if set.
2. A system copy (skipped when `CARQUET_BUNDLE_DEPS` is ON).
3. A copy it fetches and builds itself. If your project already builds the library under its usual target name (`libzstd_static`, `zlibstatic`, `lz4_static`), carquet reuses that target instead of fetching a second copy.

A project that vendors these libraries and has no system copy therefore needs no configuration. To choose explicitly, name the target:

```cmake
set(CARQUET_ZLIB_TARGET zlibstatic)
add_subdirectory(carquet)
```

A named target is linked as is, with no search and no fetch, and a name that is not a target is a configure error rather than a silent fallback. The target can be imported or built in-tree, static or shared, and must carry its include directories as a usage requirement.

## Tests, sanitizers and API docs

```bash
make test                                    # Debug build + ctest
make dev && ctest --test-dir build           # everything, including examples and interop programs
cmake -B build -DCARQUET_BUILD_DOCS=ON && cmake --build build --target docs   # Doxygen → build/docs/html/
```

[CONTRIBUTING.md](../CONTRIBUTING.md) describes the per-pull-request gate (tests, sanitizers, fuzzing, interop, benchmarks).
