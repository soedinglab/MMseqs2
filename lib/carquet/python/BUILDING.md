# Building and packaging carquet for Python

How the wheel is put together, for contributors and for people who build from source. Using the package is documented in [README.md](README.md).

## How the library is bundled

The wheel contains the pure-Python package and one shared library placed next to `__init__.py`: `carquet/libcarquet.so` (Linux), `carquet/libcarquet.dylib` (macOS) or `carquet/carquet.dll` (Windows). It is built from the repository's own `CMakeLists.txt` with `CARQUET_BUILD_SHARED=ON`, `CARQUET_BUILD_CLI=OFF`, `CARQUET_BUILD_DEV=OFF` and `CARQUET_BUNDLE_DEPS=ON`, so zstd, zlib, LZ4 and Brotli are fetched at build time and linked statically into that single file (Snappy and LZ4_RAW are implemented inside carquet). The result depends only on the C runtime and system libraries. Linux wheels keep OpenMP; `auditwheel` grafts `libgomp` into the wheel. The macOS and Windows wheels are built without OpenMP (no `libomp` / `vcomp` dependency) and parallelise reads with the library's own worker pool.

The loader (`carquet/_ffi.py`) looks for the library next to the package, then in `lib/`, `lib64/` and `bin/` under it (layouts a plain `cmake --install` produces), then on the system path. Set the `CARQUET_LIBRARY` environment variable to an explicit file to override.

## Building the wheel locally

From the repository root, with CMake >= 3.16, a C11 compiler and git (for the codec sources) available:

```bash
pip wheel ./python -w dist/          # or: python -m build python/ --outdir dist/
pip install dist/carquet-*.whl
pytest python/tests                  # PyArrow, pandas, polars and DuckDB serve as oracles: their tests skip when absent
```

`pip install -e ./python` is not supported (the shared library must be installed into the package); use a regular install. For development without reinstalling, build the shared library once (`cmake -B build-py -DCMAKE_BUILD_TYPE=Release -DCARQUET_BUILD_SHARED=ON -DCARQUET_BUNDLE_DEPS=ON && cmake --build build-py`) and run `CARQUET_LIBRARY=build-py/libcarquet.dylib PYTHONPATH=python/src pytest python/tests`. If your shell exports `CC`, `CFLAGS`, `LDFLAGS` or `SDKROOT` for another toolchain, clear them for the build (`env -u CC -u CFLAGS -u CPPFLAGS -u LDFLAGS -u SDKROOT pip wheel ./python -w dist/`).

The build is driven by `python/CMakeLists.txt`, a thin wrapper that adds the repository root as a subdirectory and installs only the runtime library (CMake component `carquet_python`) into the package directory. `pyproject.toml` (scikit-build-core) does the rest.

### Source distribution

The sdist root is `python/`, so the C tree one level up cannot be referenced as `..` once unpacked. `tool.scikit-build.sdist.force-include` therefore copies `CMakeLists.txt`, `LICENSE`, `include/`, `src/` and `cmake/` into `native/` inside the sdist, `python/CMakeLists.txt` prefers `native/` when it exists, and a `[[tool.scikit-build.overrides]]` entry (`if.from-sdist = true`) points the version regex at `native/include/carquet/carquet.h`. `python -m build --sdist python/` produces it; `pip install carquet-<version>.tar.gz` compiles it like the git checkout.

### Release wheels

`.github/workflows/wheels.yml` builds wheels with [cibuildwheel](https://cibuildwheel.pypa.io/) for CPython 3.9 - 3.13 on Linux (x86_64, aarch64; manylinux_2_28), macOS (arm64, x86_64) and Windows (AMD64), runs this test suite against each wheel, builds the sdist and publishes everything to PyPI via trusted publishing when a `v*` tag is pushed.
