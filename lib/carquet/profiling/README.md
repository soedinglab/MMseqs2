# Carquet Profiling Suite

Cross-platform profiling tools for identifying performance bottlenecks in Carquet read/write paths. Supports macOS (sample, Instruments) and Linux (perf).

## Quick Start

```bash
# Build and run full profiling
python3 profiling/run_profiler.py full

# Quick timing statistics
python3 profiling/run_profiler.py stat

# Write path profiling
python3 profiling/run_profiler.py write

# Generate flamegraph
python3 profiling/run_profiler.py flamegraph

# Run micro-benchmarks
python3 profiling/run_profiler.py micro

# Compare scalar vs SIMD
python3 profiling/run_profiler.py compare

# Just build profiling binaries
python3 profiling/run_profiler.py build
```

## Prerequisites

### macOS
- Xcode Command Line Tools (provides `sample` and `xctrace`)
- No special setup required — `sample` works out of the box

### Linux
- `perf` — Linux performance profiler
  ```bash
  sudo apt install linux-tools-generic linux-tools-$(uname -r)
  ```

### Optional (for flamegraphs)
- [FlameGraph](https://github.com/brendangregg/FlameGraph)
  ```bash
  git clone https://github.com/brendangregg/FlameGraph ~/FlameGraph
  ```
- Or [inferno](https://github.com/jonhoo/inferno) (Rust alternative)
  ```bash
  cargo install inferno
  ```

## Profiling Binaries

### `profile_read`

Full read path profiler that exercises:
- Dictionary encoding with gather operations
- RLE level decoding
- Null bitmap construction
- Various compression codecs
- SIMD dispatch paths

```bash
# Basic usage
./build/profile_read -r 5000000 -i 10

# With dictionary encoding and 10% nulls
./build/profile_read -r 10000000 -d -n 1

# With ZSTD compression
./build/profile_read -r 5000000 -c 2
```

Options:
| Flag | Description | Default |
|------|-------------|---------|
| `-r, --rows N` | Number of rows | 10000000 |
| `-b, --batch N` | Batch size | 262144 |
| `-i, --iterations N` | Test iterations | 10 |
| `-d, --dictionary` | Enable dictionary encoding | off |
| `-n, --nulls MODE` | Null ratio: 0=none, 1=10%, 2=30%, 3=50% | 0 |
| `-c, --compression N` | 0=none, 1=snappy, 2=zstd, 3=lz4 | 0 |
| `-v, --verbose` | Verbose output | off |

### `profile_write`

Full write path profiler that exercises:
- Plain, dictionary, delta, byte-stream-split encoding
- RLE level encoding
- Statistics computation (SIMD minmax)
- Compression codecs (including GZIP)
- Bloom filter insertion
- CRC32 computation
- Page index building

```bash
# Basic usage
./build/profile_write -r 5000000 -i 5

# With dictionary encoding, bloom filters, and ZSTD
./build/profile_write -r 5000000 -d -B -c 2

# Measure feature overhead (statistics, CRC, bloom, page index)
./build/profile_write -r 5000000 --overhead

# Compare compression codecs
./build/profile_write -r 5000000 --codecs

# Compare plain vs dictionary encoding
./build/profile_write -r 5000000 --encoding

# Measure column count scaling
./build/profile_write -r 5000000 --scaling

# Many columns with all features
./build/profile_write -r 1000000 -C 50 -B -P -d -c 2
```

Options:
| Flag | Description | Default |
|------|-------------|---------|
| `-r, --rows N` | Number of rows | 10000000 |
| `-g, --rowgroup N` | Row group size | 1000000 |
| `-i, --iterations N` | Test iterations | 5 |
| `-d, --dictionary` | Enable dictionary encoding | off |
| `-n, --nulls MODE` | Null ratio: 0=none, 1=10%, 2=30%, 3=50% | 0 |
| `-c, --compression N` | 0=none, 1=snappy, 2=zstd, 3=lz4, 4=gzip | 0 |
| `-C, --columns N` | Number of columns | 5 |
| `-S, --no-stats` | Disable statistics | on |
| `-B, --bloom` | Enable bloom filters | off |
| `-P, --page-index` | Enable page index | off |
| `--no-crc` | Disable CRC32 | on |
| `--overhead` | Measure feature overhead | off |
| `--codecs` | Compare compression codecs | off |
| `--encoding` | Compare plain vs dictionary | off |
| `--scaling` | Column count scaling | off |

### `profile_micro`

Micro-benchmarks for isolated component profiling:

```bash
# All components
./build/profile_micro --component all

# Specific components
./build/profile_micro --component rle --count 2000000 --iterations 500
./build/profile_micro --component encoding
./build/profile_micro --component compression
./build/profile_micro --component minmax
./build/profile_micro --component bloom
```

Components:
| Component | What it benchmarks |
|-----------|-------------------|
| `rle` | RLE encoding & decoding (single, batch, levels API) |
| `gather` | Dictionary gather (scalar vs SIMD, L1/L2/L3/memory) |
| `null` | Null bitmap (count_non_nulls, build_bitmap, scalar vs dispatch) |
| `compression` | LZ4, Snappy, ZSTD (L1/L3), GZIP — compress & decompress |
| `dispatch` | SIMD dispatch function call overhead |
| `minmax` | Min/max statistics (scalar vs dispatch, copy+minmax, all types) |
| `bss` | Byte-stream-split encode/decode (float & double, dispatch) |
| `prefix_sum` | Prefix sum / delta decoding (scalar vs dispatch) |
| `hash` | CRC32 (software), CRC32C (dispatch/HW), XXHash64 |
| `bloom` | Bloom filter insert & check throughput |
| `encoding` | Plain, delta, dictionary, BSS encode/decode |

## Profiling Modes (Python orchestrator)

### `full` — Complete Analysis
Runs stat → record → report → flamegraph in sequence.

### `stat` — Quick Timing
Runs the profiling binaries and reports timing results. On Linux with perf, also reports hardware counters.

```bash
python3 profiling/run_profiler.py stat --rows 5000000
python3 profiling/run_profiler.py stat --component rle
```

### `write` — Write Path
Runs `profile_write` with configurable features and analysis modes.

```bash
python3 profiling/run_profiler.py write --rows 5000000
python3 profiling/run_profiler.py write --overhead --codecs --encoding-cmp
```

### `record` — CPU Profile
Records a CPU profile using the best available tool:
- **macOS**: `xctrace` (Instruments) or `sample`
- **Linux**: `perf record`

```bash
python3 profiling/run_profiler.py record --rows 10000000
```

### `flamegraph` — Visual Call Graph
Generates an interactive SVG flamegraph. Records first if no data exists.

```bash
python3 profiling/run_profiler.py flamegraph --rows 10000000
```

### `compare` — Scalar vs SIMD
Benchmarks gather, null bitmap, and dispatch overhead, showing SIMD speedup ratios.

```bash
python3 profiling/run_profiler.py compare
```

### `micro` — Micro-benchmarks
Runs isolated component benchmarks via `profile_micro`.

```bash
python3 profiling/run_profiler.py micro
python3 profiling/run_profiler.py micro --component encoding
```

## Manual Profiling

### macOS

```bash
# CPU sampling (no special permissions needed)
sample ./build/profile_read -r 5000000 10 -f output.txt

# Instruments (Time Profiler)
xctrace record --template 'Time Profiler' --launch -- ./build/profile_read -r 5000000
open *.trace

# Allocations profiling
xctrace record --template 'Allocations' --launch -- ./build/profile_write -r 1000000
```

### Linux

```bash
# Record with call graph
perf record -g --call-graph dwarf ./build/profile_read -r 5000000

# View report
perf report --hierarchy

# Source annotation
perf annotate carquet_rle_decoder_get

# Hardware counters
perf stat -e cycles,instructions,cache-misses ./build/profile_read -r 1000000

# Cache analysis
perf stat -e L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses \
    ./build/profile_read -r 5000000
```

## Interpreting Results

### Key Metrics to Watch

1. **Instructions per Cycle (IPC)**
   - Good: > 1.5
   - Bad: < 0.5 (memory-bound)

2. **Cache Miss Rate**
   - L1 miss: < 5% is good
   - LLC miss: < 1% is good

3. **Branch Misprediction**
   - < 2% is good
   - > 5% indicates branchy code

### Common Bottlenecks

| Symptom | Cause | Solution |
|---------|-------|----------|
| Low IPC, high LLC misses | Random memory access | Add prefetching |
| High branch mispredictions | Unpredictable conditionals | Branchless algorithms |
| High function call overhead | Small functions called often | Inline or batch |
| Visible in flamegraph | Hot function | Optimize that function |

## Output Files

Profiling output is saved to `profiling/output/`:

| File | Description |
|------|-------------|
| `*_stat.txt` | CPU statistics |
| `*_perf.data` | Raw perf recording (Linux) |
| `*_sample.txt` | CPU sample data (macOS) |
| `*.trace` | Instruments trace (macOS) |
| `*_report.txt` | Function report |
| `*_flamegraph.svg` | Interactive flamegraph |
| `*_icicle.svg` | Reversed (icicle) flamegraph |
| `*_micro.txt` | Micro-benchmark results |
| `*_write.txt` | Write profiling results |
| `*_compare.txt` | Scalar vs SIMD comparison |
| `*_annotations/` | Source annotations (Linux) |

## Architecture-Specific Notes

### x86-64
- SSE4.2, AVX2, AVX-512 paths available
- Prefetching critical for dictionary gather
- Check dispatch overhead vs inline SIMD

### ARM64 (Apple Silicon)
- NEON paths typically faster than x86 SSE
- Hardware CRC32C acceleration
- SVE available but experimental
- Less dispatch overhead due to simpler calling convention
