# Integration recipes

Copy-paste build wiring for the common setups. The CMake target is always `carquet`; headers are under `include/carquet/`, the umbrella header is `<carquet/carquet.h>`.

## CMake: installed package

```cmake
find_package(carquet CONFIG REQUIRED)          # carquetConfig.cmake from `cmake --install`
add_executable(app main.c)
target_link_libraries(app PRIVATE carquet)     # headers, zstd/zlib/lz4, Threads/OpenMP propagate
```

Install carquet once with `cmake -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build && sudo cmake --install build` (or `make && sudo make install` at the repo root). A non-default prefix needs `-DCMAKE_PREFIX_PATH=/opt/carquet` on the consumer side.

## CMake: FetchContent (no install step, pinned version)

```cmake
include(FetchContent)
set(CARQUET_BUILD_CLI OFF CACHE BOOL "" FORCE)   # library only
set(CARQUET_BUILD_DEV OFF CACHE BOOL "" FORCE)
FetchContent_Declare(carquet
    GIT_REPOSITORY https://github.com/Vitruves/carquet.git
    GIT_TAG        v0.7.2)                      # pin a release tag
FetchContent_MakeAvailable(carquet)
target_link_libraries(app PRIVATE carquet)
```

carquet's own `FetchContent` fallback fetches zstd, zlib, lz4 and brotli when they are not on the system; nothing else is needed offline-wise beyond network access at configure time. `-DCARQUET_ENABLE_BROTLI=OFF` removes the brotli dependency (the codec then returns `CARQUET_ERROR_UNSUPPORTED_CODEC`).

## CMake: vendored subdirectory

```cmake
set(CARQUET_BUILD_CLI OFF CACHE BOOL "" FORCE)
set(CARQUET_BUILD_DEV OFF CACHE BOOL "" FORCE)
add_subdirectory(third_party/carquet)
target_link_libraries(app PRIVATE carquet)
```

The three-way resolution used by the shipped example consumer (`mocklib/CMakeLists.txt`) is a good template: reuse a `carquet` target if the enclosing project already defines one, else `find_package`, else `add_subdirectory`.

## Reusing the project's own zstd / zlib / lz4

If the project already builds these (static or shared, imported or in-tree), name the targets before carquet is added so it links them instead of searching or fetching:

```cmake
set(CARQUET_ZSTD_TARGET   libzstd_static)
set(CARQUET_ZLIB_TARGET   zlibstatic)
set(CARQUET_LZ4_TARGET    lz4_static)
set(CARQUET_BROTLI_TARGET "brotlienc;brotlidec")   # encoder + decoder targets
add_subdirectory(third_party/carquet)
```

A name that is not a target is a configure error (no silent fallback). The target must carry its include directories as a usage requirement. `CARQUET_BUNDLE_DEPS=ON` does the opposite: always fetch and static-bundle the three codecs for portable artifacts.

## Options that matter to consumers

| Option | Default | Set when |
|---|---|---|
| `CARQUET_BUILD_SHARED` | OFF | you want `libcarquet.so` / `.dylib` / `.dll` |
| `CARQUET_BUILD_CLI` | ON | set OFF in embedded builds (no `carquet` binary) |
| `CARQUET_BUILD_DEV` | OFF | tests, examples, benchmarks, interop (never for consumers) |
| `CARQUET_NATIVE_ARCH` | OFF | host-only binaries; otherwise SIMD is dispatched at runtime |
| `CARQUET_ENABLE_SSE/AVX/AVX2/AVX512/NEON` | ON | disable to shrink or to work around a toolchain |
| `CARQUET_ENABLE_BROTLI` | ON | OFF drops libbrotli (encoder is ~0.5–1 MB static); BROTLI then needs `carquet_register_codec()` |
| `CMAKE_INTERPROCEDURAL_OPTIMIZATION` | ON in Release | set OFF if the consumer links with a different compiler than the one that built carquet (LTO bitcode is not portable across compilers) |

## pkg-config / plain Makefiles

```sh
CFLAGS  += $(shell pkg-config --cflags carquet)
LDFLAGS += $(shell pkg-config --libs carquet)
# static carquet: add the private libs too
LDFLAGS += $(shell pkg-config --libs --static carquet)      # -lzstd -lz -llz4 -lbrotlienc -lbrotlidec -lbrotlicommon (+ OpenMP runtime when built with it)
```

Manual link line on macOS with Homebrew codecs: `cc main.c -I<prefix>/include <prefix>/lib/libcarquet.a -lzstd -lz -llz4 $(brew --prefix libomp)/lib/libomp.dylib`.

## xmake

carquet ships an `xmake.lua` mirroring the CMake options (`xmake f --shared=y`, `--dev=y`, `--openmp=n`, ...). Build it with `xmake`; codecs are linked statically, so a consumer only needs `add_includedirs("carquet/include")` plus `add_linkdirs`/`add_links("carquet")` on the build output, or the CMake route above.

## C++ projects

The header is `extern "C"`-wrapped; include it directly. Wrap handles in RAII (`std::unique_ptr<carquet_reader_t, decltype(&carquet_reader_close)>`) and translate `carquet_status_t` to exceptions at the boundary if that is the project's convention. Designated initialisers in the examples (`CARQUET_ERROR_INIT`) are C11; in C++ write `carquet_error_t err{};`.

## Windows / MSVC

Builds with MSVC 2015+ and MinGW. No `%zu`-style portability traps in the public API. The library uses the Win32 TLS API for per-thread codec contexts; nothing to configure. Static linking is the default; for a DLL set `CARQUET_BUILD_SHARED=ON`.

## Sanity check after wiring

```c
#include <carquet/carquet.h>
#include <stdio.h>
int main(void) { printf("carquet %s\n", carquet_version()); return 0; }
```

Then run the round trip from `examples/basic_write_read.c` and open the output with `carquet schema out.parquet` and, if available, PyArrow or DuckDB.
