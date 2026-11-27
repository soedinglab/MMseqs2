<h1 align="center">
  <img src="res/img/carquet_banner.png" alt="Carquet" width="840" />
</h1>

<p align="center">
  Eine schnelle, reine C-Bibliothek zum Lesen und Schreiben von Apache-Parquet-Dateien.
</p>

<p align="center">
  <a href="https://vitruves.github.io/carquet/"><img src="https://img.shields.io/badge/Website-vitruves.github.io%2Fcarquet-2a78d6?style=for-the-badge" alt="Webseite" /></a>
</p>

<p align="center">
  <a href="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml"><img src="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml/badge.svg" alt="Build" /></a>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue" alt="Plattform" />
  <img src="https://img.shields.io/badge/C-C11-blue" alt="C-Standard" />
  <a href="https://github.com/Vitruves/carquet/releases"><img src="https://img.shields.io/github/v/release/Vitruves/carquet?color=blue" alt="Release" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="Lizenz" /></a>
  <a href="https://github.com/Vitruves/carquet/graphs/traffic"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/Vitruves/carquet/traffic/clones.json" alt="Clones" /></a>
  <br/>
  <a href="https://pypi.org/project/carquet/"><img src="https://img.shields.io/pypi/v/carquet?color=blue&label=pypi" alt="PyPI" /></a>
  <a href="https://pypi.org/project/carquet/"><img src="https://img.shields.io/pypi/pyversions/carquet" alt="Python-Versionen" /></a>
  <a href="https://pepy.tech/projects/carquet"><img src="https://img.shields.io/pepy/dt/carquet?color=blue&label=pypi%20downloads" alt="PyPI-Downloads" /></a>
  <img src="https://img.shields.io/badge/SIMD-SSE4.2%20%C2%B7%20AVX2%20%C2%B7%20AVX--512%20%C2%B7%20NEON%20%C2%B7%20SVE-eb6834" alt="SIMD" />
</p>

<p align="center">
  <a href="#highlights">Highlights</a> &nbsp;·&nbsp;
  <a href="#benchmarks">Benchmarks</a> &nbsp;·&nbsp;
  <a href="#installation">Installation</a> &nbsp;·&nbsp;
  <a href="#build-konfigurationen">Build-Konfigurationen</a> &nbsp;·&nbsp;
  <a href="#c-api">C-API</a> &nbsp;·&nbsp;
  <a href="#python-api">Python-API</a> &nbsp;·&nbsp;
  <a href="#cli-tool">CLI</a> &nbsp;·&nbsp;
  <a href="#unterstützte-parquet-funktionen">Funktionsumfang</a> &nbsp;·&nbsp;
  <a href="docs/README.md">Handbuch</a> &nbsp;·&nbsp;
  <a href="CHANGELOG.md">Changelog</a>
</p>

<p align="center">
  <a href="README.md">English</a> &nbsp;·&nbsp;
  <a href="README.zh-CN.md">简体中文</a> &nbsp;·&nbsp;
  <a href="README.ja.md">日本語</a> &nbsp;·&nbsp;
  <a href="README.ko.md">한국어</a> &nbsp;·&nbsp;
  <b>Deutsch</b> &nbsp;·&nbsp;
  <a href="README.fr.md">Français</a> &nbsp;·&nbsp;
  <a href="README.pt-BR.md">Português (Brasil)</a>
</p>

<br/>

---

<!-- ────────────────────────────  HIGHLIGHTS  ──────────────────────────── -->

## Highlights

<div align="center">

| | |
|---|---|
| ⚡ **Schnell** | Schreibt **2.0–5.6x schneller** als Arrow C++ und liest komprimierte Daten aus ein und derselben Datei **1.5–2.1x schneller** (Apple M3, 10 Mio. Zeilen). Unkomprimiertes Lesen wird in **0.28 ms** direkt aus dem mmap bedient. SIMD für SSE4.2 / AVX2 / AVX-512 / NEON / SVE mit Laufzeit-Dispatch. |
| 🪶 **Klein** | Reines C11, eine **~500 KB** große Bibliothek (Arrow C++: 50 MB+). Vier Kompressionsabhängigkeiten (zstd, zlib, lz4, brotli), die bei Bedarf automatisch nachgeladen werden. Baut in unter einer Minute mit einem einzigen `make` oder `xmake`. |
| 📦 **Vollständig** | Alle physischen Typen und Encodings, Snappy / GZIP / LZ4 / ZSTD / Brotli, V1- + V2-Pages, verschachtelte Schemata in beliebiger Tiefe, Bloom-Filter, Page-Indexe, statistikbasiertes Pruning, Append. [Vollständige Matrix ↓](#unterstützte-parquet-funktionen) |
| 🔗 **Interoperabel** | Dateien werden out of the box von PyArrow, DuckDB und Spark gelesen. Abhängigkeitsfreie [Arrow C Data Interface](docs/reading.md#export-to-arrow-c-data-interface)-Brücke (struct / list / map in beliebiger Tiefe) und eine vollständige [Python-API](#python-api) ohne Abhängigkeiten. |
| 🛠 **Werkzeuge** | `carquet`-CLI zum Inspizieren, Filtern, Exportieren und Validieren von Dateien, dazu `codegen`, das einen kompilierfertigen C-Reader erzeugt. Ein fertiger [Agent Skill](#ai-coding-agents) für Claude Code und Codex. |

</div>

<br/>

---

<!-- ────────────────────────────  BENCHMARKS  ──────────────────────────── -->

## Benchmarks

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="res/img/bench_10m_m3_dark.png">
    <source media="(prefers-color-scheme: light)" srcset="res/img/bench_10m_m3_light.png">
    <img src="res/img/bench_10m_m3_light.png" alt="Schreib- und Lesezeit für 10 Millionen Zeilen auf Apple M3: carquet vs. Arrow C++ und Python-Parquet-Bibliotheken, pro Codec" width="900">
  </picture>
</p>

<p align="center"><b><a href="https://vitruves.github.io/carquet/">Den interaktiven Benchmark auf der Webseite erkunden →</a></b><br/><sub>zwischen Lesen und Schreiben sowie zwischen 100K, 1M und 10M Zeilen umschalten</sub></p>

Beschleunigung von carquet gegenüber Arrow C++ bei 10 Mio. Zeilen (höher = carquet schneller):

<div align="center">

| | ARM · Apple M3 | | x86 · Windows (Haswell, MSVC) | | x86 · Linux (Xeon D-1531) | |
|---|---:|---:|---:|---:|---:|---:|
| **Codec** | **Schreiben** | **Lesen** | **Schreiben** | **Lesen** | **Schreiben** | **Lesen** |
| snappy | **4.77x** | **1.63x** | **3.39x** | **1.17x** | **1.55x** | **1.25x** |
| zstd | **5.56x** | **1.40x** | **3.83x** | **1.13x** | **1.31x** | **1.04x** |
| lz4 | **4.62x** | **1.85x** | **2.96x** | **1.05x** | **1.02x** | 0.83x |
| none | **2.04x** | **52.0x**\* | **2.10x** | **62.0x**\* | **1.13x** | **40.6x**\* |

</div>

<sub>\* mmap-Zero-Copy-Pfad. Apple M3: carquet 0.8.0 vs. Arrow C++ 26.0.1 (2026-11-09). Windows: carquet 0.7.3 vs. Arrow C++ 24.0.0 (2026-11-02). Linux: carquet 0.4.4 vs. Arrow C++ 23.0.1, vor der Writer-Parallelisierung von 0.7.1, seither nicht neu gemessen.</sub>

Komprimiertes Lesen dekodiert jeden einzelnen Wert, und unter Linux und macOS linken beide Bibliotheken dieselben lz4 / zstd des Systems, die Codec-Geschwindigkeit ist also identisch. Der fairste Vergleich ist der **Cross-Read derselben Datei** (beide Bibliotheken lesen exakt dieselbe Datei): carquet ist auf M3 **1.5–2.1x schneller** als Arrow C++ (2.0–4.1x schneller als polars), unter Linux 1.6–2.6x und auf der Windows-VM 1.04–1.30x, deren kalte Lesezugriffe durch die virtuelle Festplatte begrenzt sind.

<details>
<summary><b>Benchmark-Methodik</b></summary>

<br/>

Alle Benchmarks verwenden identische Daten (deterministischer LCG-PRNG), identische Parquet-Einstellungen (kein Dictionary, BYTE_STREAM_SPLIT für Gleitkommazahlen, Page-Prüfsummen, mmap-Lesen) und den getrimmten Median aus 11–51 Iterationen, wobei der Page-Cache des Betriebssystems zwischen Schreib- und Lesephase geleert wird und zwischen den Konfigurationen eine Abkühlpause liegt. Schema: 3 Spalten (INT64, DOUBLE, INT32). Verglichen wird mit dem Low-Level-Parquet-Reader von Arrow C++ (unter Umgehung der Arrow-Table-Materialisierung), mit PyArrow, mit fastparquet sowie mit dem nativen (Rust-)Reader und -Writer von polars.

Der Page-Cache wird pro Datei geleert: unter Linux mit `posix_fadvise(DONTNEED)`; unter macOS und Windows wird die Datei ohne Caching kopiert (`F_NOCACHE`, `FILE_FLAG_NO_BUFFERING`) und durch die Kopie ersetzt (`benchmark/bench_purge_win32.h`). Unter Windows linkt carquet die selbst gebauten zstd / lz4 und Arrow C++ die von conda-forge, sodass sich die Codec-Geschwindigkeit zwischen beiden leicht unterscheiden kann.

Der Benchmark **Cross-Read derselben Datei** ist der fairste Vergleich: Beide Bibliotheken lesen exakt dieselbe Parquet-Datei (von einer geschrieben, von beiden gelesen). Das eliminiert Unterschiede bei Page-Größen, Encoding-Entscheidungen und Row-Group-Layout.

Mit \* markiertes **unkomprimiertes Lesen** nutzt Carquets **mmap-Zero-Copy-Pfad**: Für PLAIN-kodierte, unkomprimierte, required Spalten fester Größe liefert der Batch-Reader Zeiger direkt in die speicherabgebildete Datei, ohne memcpy. Arrow materialisiert immer in eigene Puffer. **Die Werte für komprimiertes Lesen sind das repräsentativste Maß für den End-to-End-Lesedurchsatz.**

</details>

<details>
<summary><b>Vollständige x86-Ergebnisse</b> (Intel Xeon D-1531, Linux)</summary>

<br/>

*12 Threads @ 2.2GHz, 32GB RAM, Ubuntu 24.04, ZSTD-Level 1*

#### 10 Mio. Zeilen vs. Arrow C++

<div align="center">

| Codec | Carquet Schreiben | Arrow C++ Schreiben | Verh. Schreiben | Carquet Lesen | Arrow C++ Lesen | Verh. Lesen | Größe |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **1557ms** | 1766ms | **1.13x** | **1.25ms** | 50.8ms | **40.6x**\* | 190.7MB |
| snappy | **1002ms** | 1549ms | **1.55x** | **78ms** | 97.8ms | **1.25x** | 125.1MB |
| zstd | **1311ms** | 1714ms | **1.31x** | **76.8ms** | 80.2ms | **1.04x** | 95.3MB |
| lz4 | **1521ms** | 1554ms | **1.02x** | 59.1ms | **49.0ms** | 0.83x | 122.9MB |

</div>

<br/>

#### 1 Mio. Zeilen vs. Arrow C++

<div align="center">

| Codec | Carquet Schreiben | Arrow C++ Schreiben | Verh. Schreiben | Carquet Lesen | Arrow C++ Lesen | Verh. Lesen |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **180ms** | 196ms | **1.09x** | **0.22ms** | 6.2ms | **28x**\* |
| snappy | **141ms** | 148ms | **1.05x** | **8.1ms** | 11.6ms | **1.44x** |
| zstd | **131ms** | 185ms | **1.41x** | 10.3ms | **9.1ms** | 0.88x |
| lz4 | **143ms** | 149ms | **1.04x** | 8.5ms | **6.1ms** | 0.72x |

</div>

<br/>

#### 100K Zeilen vs. Arrow C++

<div align="center">

| Codec | Carquet Schreiben | Arrow C++ Schreiben | Verh. Schreiben | Carquet Lesen | Arrow C++ Lesen | Verh. Lesen |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **14.1ms** | 18.4ms | **1.30x** | **0.11ms** | 2.18ms | **19.8x**\* |
| snappy | **10.1ms** | 10.6ms | **1.05x** | **1.27ms** | 5.97ms | **4.70x** |
| zstd | **8.7ms** | 14.1ms | **1.62x** | **1.58ms** | 3.88ms | **2.46x** |
| lz4 | **9.6ms** | 11.0ms | **1.14x** | **0.77ms** | 2.78ms | **3.61x** |

</div>

<br/>

#### Cross-Read derselben Datei (10 Mio. Zeilen)

Beide Bibliotheken lesen **dieselbe** Parquet-Datei, der fairste direkte Vergleich.

<div align="center">

| Codec | Writer | Carquet Lesen | Arrow C++ Lesen | Verhältnis |
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

#### 10 Mio. Zeilen vs. PyArrow

<div align="center">

| Codec | Carquet Schreiben | PyArrow Schreiben | Verh. Schreiben | Carquet Lesen | PyArrow Lesen | Verh. Lesen |
|-------|--------------|---------------|---------|-------------|--------------|---------|
| none | **1557ms** | 1806ms | **1.16x** | **1.25ms** | 213ms | **170x**\* |
| snappy | **1002ms** | 1649ms | **1.65x** | **78ms** | 384ms | **4.91x** |
| zstd | **1311ms** | 1796ms | **1.37x** | **76.8ms** | 369ms | **4.81x** |
| lz4 | **1521ms** | 1676ms | **1.10x** | **59.1ms** | 281ms | **4.76x** |

</div>

\* Zero-Copy-mmap-Pfad

</details>

<details>
<summary><b>Vollständige ARM-Ergebnisse</b> (Apple M3, macOS)</summary>

<br/>

*Carquet 0.8.0, MacBook Air M3, 16GB RAM, macOS 26.5, Arrow C++ 26.0.1, PyArrow 26.0.1, fastparquet 2026.9.0, polars 2.1.0, ZSTD-Level 1, gemessen am 2026-11-09 (10M-Zeilen-Konfigurationen aus einem zweiten Lauf am selben Tag: die erste Konfiguration des ersten Laufs war gestört). Rohbericht: [`benchmark/bench_0.8.0_darwin_arm64_20261109.json`](benchmark/bench_0.8.0_darwin_arm64_20261109.json).*

#### 10 Mio. Zeilen vs. Arrow C++

<div align="center">

| Codec | Carquet Schreiben | Arrow C++ Schreiben | Verh. Schreiben | Carquet Lesen | Arrow C++ Lesen | Verh. Lesen | Größe |
|-------|--------------|---------|---------|-------------|---------|---------|------|
| none | **57.7ms** | 117ms | **2.04x** | **0.28ms** | 14.6ms | **52.0x**\* | 190.7MB |
| snappy | **49.3ms** | 235ms | **4.77x** | **14.0ms** | 22.7ms | **1.63x** | 125.1MB |
| zstd | **61.0ms** | 339ms | **5.56x** | **20.6ms** | 28.8ms | **1.40x** | 95.3MB |
| lz4 | **52.9ms** | 244ms | **4.62x** | **9.83ms** | 18.2ms | **1.85x** | 122.7MB |

</div>

<br/>

#### 1 Mio. Zeilen vs. Arrow C++

<div align="center">

| Codec | Carquet Schreiben | Arrow C++ Schreiben | Verh. Schreiben | Carquet Lesen | Arrow C++ Lesen | Verh. Lesen |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **5.87ms** | 12.3ms | **2.10x** | **0.06ms** | 1.59ms | **26.5x**\* |
| snappy | **7.40ms** | 25.9ms | **3.49x** | **1.30ms** | 2.53ms | **1.95x** |
| zstd | **8.32ms** | 34.1ms | **4.10x** | **2.36ms** | 3.17ms | **1.34x** |
| lz4 | **7.70ms** | 25.7ms | **3.34x** | **1.06ms** | 2.15ms | **2.03x** |

</div>

<br/>

#### 100K Zeilen vs. Arrow C++

<div align="center">

| Codec | Carquet Schreiben | Arrow C++ Schreiben | Verh. Schreiben | Carquet Lesen | Arrow C++ Lesen | Verh. Lesen |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **0.87ms** | 1.65ms | **1.90x** | **0.02ms** | 0.23ms | **11.5x**\* |
| snappy | **1.52ms** | 2.62ms | **1.72x** | **0.38ms** | 0.91ms | **2.39x** |
| zstd | **1.90ms** | 3.69ms | **1.94x** | **0.64ms** | 1.32ms | **2.06x** |
| lz4 | **1.66ms** | 2.63ms | **1.58x** | **0.30ms** | 0.63ms | **2.10x** |

</div>

<br/>

#### 10 Mio. Zeilen vs. PyArrow

<div align="center">

| Codec | Carquet Schreiben | PyArrow Schreiben | Verh. Schreiben | Carquet Lesen | PyArrow Lesen | Verh. Lesen |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **57.7ms** | 165ms | **2.87x** | **0.28ms** | 26.1ms | **93.1x**\* |
| snappy | **49.3ms** | 276ms | **5.60x** | **14.0ms** | 39.3ms | **2.82x** |
| zstd | **61.0ms** | 387ms | **6.34x** | **20.6ms** | 41.3ms | **2.00x** |
| lz4 | **52.9ms** | 287ms | **5.43x** | **9.83ms** | 29.0ms | **2.95x** |

</div>

<br/>

#### 1 Mio. Zeilen vs. PyArrow

<div align="center">

| Codec | Carquet Schreiben | PyArrow Schreiben | Verh. Schreiben | Carquet Lesen | PyArrow Lesen | Verh. Lesen |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **5.87ms** | 17.5ms | **2.98x** | **0.06ms** | 2.79ms | **46.5x**\* |
| snappy | **7.40ms** | 28.9ms | **3.90x** | **1.30ms** | 3.61ms | **2.78x** |
| zstd | **8.32ms** | 39.4ms | **4.74x** | **2.36ms** | 4.00ms | **1.69x** |
| lz4 | **7.70ms** | 30.5ms | **3.97x** | **1.06ms** | 3.22ms | **3.04x** |

</div>

<br/>

#### 100K Zeilen vs. PyArrow

<div align="center">

| Codec | Carquet Schreiben | PyArrow Schreiben | Verh. Schreiben | Carquet Lesen | PyArrow Lesen | Verh. Lesen |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **0.87ms** | 2.10ms | **2.41x** | **0.02ms** | 0.26ms | **13.0x**\* |
| snappy | **1.52ms** | 3.00ms | **1.97x** | **0.38ms** | 0.67ms | **1.76x** |
| zstd | **1.90ms** | 4.13ms | **2.17x** | **0.64ms** | 0.73ms | **1.14x** |
| lz4 | **1.66ms** | 3.04ms | **1.83x** | **0.30ms** | 0.44ms | **1.47x** |

</div>

<br/>

#### 10 Mio. Zeilen vs. fastparquet

<div align="center">

| Codec | Carquet Schreiben | fastparquet Schreiben | Verh. Schreiben | Carquet Lesen | fastparquet Lesen | Verh. Lesen |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **57.7ms** | 64.6ms | **1.12x** | **0.28ms** | 55.7ms | **198.8x**\* |
| snappy | **49.3ms** | 208ms | **4.22x** | **14.0ms** | 103ms | **7.40x** |
| zstd | **61.0ms** | 325ms | **5.33x** | **20.6ms** | 185ms | **8.95x** |
| lz4 | **52.9ms** | 244ms | **4.61x** | **9.83ms** | 55.4ms | **5.64x** |

</div>

<br/>

#### 1 Mio. Zeilen vs. fastparquet

<div align="center">

| Codec | Carquet Schreiben | fastparquet Schreiben | Verh. Schreiben | Carquet Lesen | fastparquet Lesen | Verh. Lesen |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **5.87ms** | 8.01ms | **1.36x** | **0.06ms** | 3.08ms | **51.3x**\* |
| snappy | **7.40ms** | 22.8ms | **3.08x** | **1.30ms** | 10.2ms | **7.85x** |
| zstd | **8.32ms** | 33.5ms | **4.03x** | **2.36ms** | 18.0ms | **7.61x** |
| lz4 | **7.70ms** | 22.4ms | **2.92x** | **1.06ms** | 5.37ms | **5.07x** |

</div>

<br/>

#### 100K Zeilen vs. fastparquet

<div align="center">

| Codec | Carquet Schreiben | fastparquet Schreiben | Verh. Schreiben | Carquet Lesen | fastparquet Lesen | Verh. Lesen |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **0.87ms** | 1.08ms | **1.24x** | **0.02ms** | 0.58ms | **29.0x**\* |
| snappy | **1.52ms** | 2.43ms | **1.60x** | **0.38ms** | 1.20ms | **3.16x** |
| zstd | **1.90ms** | 3.47ms | **1.83x** | **0.64ms** | 2.06ms | **3.22x** |
| lz4 | **1.66ms** | 2.31ms | **1.39x** | **0.30ms** | 0.73ms | **2.43x** |

</div>

<br/>

#### 10 Mio. Zeilen vs. polars

<div align="center">

| Codec | Carquet Schreiben | polars Schreiben | Verh. Schreiben | Carquet Lesen | polars Lesen | Verh. Lesen | Carquet Größe | polars Größe |
|-------|--------------|---------|---------|-------------|---------|---------|------|------|
| none | 57.7ms | **37.5ms** | 0.65x | **0.28ms** | 28.7ms | **102.4x**\* | 190.7MB | 105.1MB |
| snappy | 49.3ms | **38.1ms** | 0.77x | **14.0ms** | 28.3ms | **2.03x** | 125.1MB | 103.8MB |
| zstd | 61.0ms | **48.6ms** | 0.80x | **20.6ms** | 42.1ms | **2.04x** | 95.3MB | 98.1MB |
| lz4 | 52.9ms | **40.0ms** | 0.76x | **9.83ms** | 29.4ms | **2.99x** | 122.7MB | 104.2MB |

</div>

<br/>

#### 1 Mio. Zeilen vs. polars

<div align="center">

| Codec | Carquet Schreiben | polars Schreiben | Verh. Schreiben | Carquet Lesen | polars Lesen | Verh. Lesen | Carquet Größe | polars Größe |
|-------|--------------|---------|---------|-------------|---------|---------|------|------|
| none | 5.87ms | **4.53ms** | 0.77x | **0.06ms** | 3.14ms | **52.3x**\* | 19.1MB | 12.6MB |
| snappy | 7.40ms | **4.82ms** | 0.65x | **1.30ms** | 3.25ms | **2.50x** | 12.5MB | 11.5MB |
| zstd | 8.32ms | **5.83ms** | 0.70x | **2.36ms** | 4.12ms | **1.75x** | 9.6MB | 10.1MB |
| lz4 | 7.70ms | **4.77ms** | 0.62x | **1.06ms** | 3.26ms | **3.08x** | 12.3MB | 11.5MB |

</div>

<br/>

#### 100K Zeilen vs. polars

<div align="center">

| Codec | Carquet Schreiben | polars Schreiben | Verh. Schreiben | Carquet Lesen | polars Lesen | Verh. Lesen | Carquet Größe | polars Größe |
|-------|--------------|---------|---------|-------------|---------|---------|------|------|
| none | 0.87ms | **0.72ms** | 0.83x | **0.02ms** | 0.33ms | **16.5x**\* | 1.9MB | 1.3MB |
| snappy | 1.52ms | **0.91ms** | 0.60x | **0.38ms** | 0.46ms | **1.21x** | 1.3MB | 1.1MB |
| zstd | 1.90ms | **1.05ms** | 0.55x | **0.64ms** | 1.08ms | **1.69x** | 1.0MB | 1.0MB |
| lz4 | 1.66ms | **0.96ms** | 0.58x | **0.30ms** | 0.44ms | **1.47x** | 1.2MB | 1.1MB |

</div>

<br/>

#### Cross-Read derselben Datei (10 Mio. Zeilen)

Beide Bibliotheken lesen **dieselbe** Parquet-Datei, der fairste direkte Vergleich.

<div align="center">

| Codec | Writer | Carquet Lesen | Arrow C++ Lesen | Verhältnis |
|-------|--------|-------------|----------------|-------|
| none | Carquet | **0.30ms** | 17.1ms | **56.8x**\* |
| none | Arrow | **0.94ms** | 17.4ms | **18.5x**\* |
| snappy | Carquet | **12.7ms** | 23.3ms | **1.83x** |
| snappy | Arrow | **11.5ms** | 24.4ms | **2.13x** |
| zstd | Carquet | **20.2ms** | 30.3ms | **1.50x** |
| zstd | Arrow | **19.8ms** | 30.7ms | **1.55x** |
| lz4 | Carquet | **9.96ms** | 19.2ms | **1.93x** |
| lz4 | Arrow | **8.90ms** | 17.8ms | **2.00x** |

</div>

<br/>

#### Cross-Read derselben Datei vs. polars (10 Mio. Zeilen)

Jede Bibliothek liest die Datei, die der carquet-Benchmark schreibt, und die Datei, die polars mit seinen eigenen Einstellungen schreibt (Dictionary-Encoding auf `id` und `category`).

<div align="center">

| Codec | Writer | Carquet Lesen | polars Lesen | Verhältnis | Datei |
|-------|--------|-------------|-------------|-------|------|
| none | Carquet | **0.28ms** | 36.4ms | **130.0x**\* | 190.7MB |
| none | polars | **10.3ms** | 27.8ms | **2.71x** | 105.1MB |
| snappy | Carquet | **13.2ms** | 53.9ms | **4.09x** | 125.1MB |
| snappy | polars | **12.4ms** | 30.5ms | **2.45x** | 103.8MB |
| zstd | Carquet | **21.7ms** | 51.8ms | **2.39x** | 95.3MB |
| zstd | polars | **17.4ms** | 34.7ms | **1.99x** | 98.1MB |
| lz4 | Carquet | **10.7ms** | 37.7ms | **3.52x** | 122.7MB |
| lz4 | polars | **11.5ms** | 26.7ms | **2.31x** | 104.2MB |

</div>

\* Zero-Copy-mmap-Pfad. fastparquet liest über pandas und kann weder BYTE_STREAM_SPLIT noch Page-Prüfsummen schreiben, daher sind seine Zeilen weniger direkt vergleichbar.

polars verwendet seinen eigenen Rust-Parquet-Reader und -Writer (`use_pyarrow=False`, sein Standard), mit denselben Daten, derselben Row-Group-Größe, demselben zstd-Level, denselben Statistiken und LZ4_RAW für `lz4`. Es kann weder BYTE_STREAM_SPLIT noch Page-Prüfsummen schreiben und hat keinen Schalter für Dictionary-Encoding: Es kodiert `id` und `category` von sich aus per Dictionary, daher sind seine Dateien auf diesen Daten unkomprimiert 45% kleiner als die dictionary-freien von carquet und mit snappy / lz4 15–17% kleiner (mit zstd etwa gleich groß), und seine Schreibvorgänge schicken weniger Bytes durch Codec und Festplatte. Mit Dictionary-Encoding auf denselben beiden Spalten schreibt carquet die 10 Mio. Zeilen in 45 ms (none), 51 ms (snappy), 51 ms (zstd) und 48 ms (lz4), gegenüber 38, 38, 49 und 40 ms bei polars. Lesezugriffe vergleicht die Tabelle zum Cross-Read derselben Datei fair: carquet liest die Datei jeder der beiden Bibliotheken mit Codec 2.0–4.1x schneller.

</details>

<details>
<summary><b>Vollständige Windows-Ergebnisse</b> (Intel-Haswell-VM, MSVC)</summary>

<br/>

*Carquet 0.7.3, gebaut mit MSVC 19.44 (x64 Release, OpenMP), OpenStack-VM mit 4 vCPU (Intel Haswell, 3.1 GHz), 14.6GB RAM, virtuelle SAS-Festplatte, Windows Server 2025, Arrow C++ 24.0.0 (conda-forge), ZSTD-Level 1, gemessen am 2026-11-02. Kalte Lesezugriffe gehen auf die virtuelle Festplatte, was die Leseverhältnisse verringert; Schreiben und der Zero-Copy-Pfad sind CPU-gebunden.*

#### 10 Mio. Zeilen vs. Arrow C++

<div align="center">

| Codec | Carquet Schreiben | Arrow C++ Schreiben | Verh. Schreiben | Carquet Lesen | Arrow C++ Lesen | Verh. Lesen | Größe |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **182ms** | 382ms | **2.10x** | **2.04ms** | 126.5ms | **62.0x**\* | 190.7MB |
| snappy | **175ms** | 594ms | **3.39x** | **115.6ms** | 134.9ms | **1.17x** | 125.1MB |
| zstd | **242ms** | 927ms | **3.83x** | **140.7ms** | 159.2ms | **1.13x** | 95.3MB |
| lz4 | **178ms** | 527ms | **2.96x** | **95.5ms** | 99.9ms | **1.05x** | 122.9MB |

</div>

<br/>

#### 1 Mio. Zeilen vs. Arrow C++

<div align="center">

| Codec | Carquet Schreiben | Arrow C++ Schreiben | Verh. Schreiben | Carquet Lesen | Arrow C++ Lesen | Verh. Lesen |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **15.6ms** | 37.6ms | **2.41x** | **0.44ms** | 13.72ms | **31.2x**\* |
| snappy | **25.4ms** | 61.2ms | **2.41x** | **14.01ms** | 16.31ms | **1.16x** |
| zstd | **30.4ms** | 94.3ms | **3.10x** | **17.09ms** | 18.38ms | **1.08x** |
| lz4 | **26.4ms** | 53.8ms | **2.04x** | 13.07ms | **12.89ms** | 0.99x |

</div>

<br/>

#### 100K Zeilen vs. Arrow C++

<div align="center">

| Codec | Carquet Schreiben | Arrow C++ Schreiben | Verh. Schreiben | Carquet Lesen | Arrow C++ Lesen | Verh. Lesen |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **3.53ms** | 3.95ms | **1.12x** | **0.16ms** | 2.13ms | **13.3x**\* |
| snappy | **5.60ms** | 6.23ms | **1.11x** | **2.79ms** | 3.92ms | **1.41x** |
| zstd | **7.11ms** | 9.79ms | **1.38x** | **2.86ms** | 4.36ms | **1.52x** |
| lz4 | **5.75ms** | 5.82ms | **1.01x** | 2.22ms | **2.15ms** | 0.97x |

</div>

<br/>

#### Cross-Read derselben Datei (10 Mio. Zeilen)

<div align="center">

| Codec | Writer | Carquet Lesen | Arrow C++ Lesen | Verhältnis |
|-------|--------|-------------|----------------|-------|
| none | Carquet | **2.65ms** | 116.9ms | **44.1x**\* |
| none | Arrow | **5.29ms** | 126.6ms | **23.9x**\* |
| snappy | Carquet | **115.2ms** | 135.0ms | **1.17x** |
| snappy | Arrow | **112.7ms** | 146.2ms | **1.30x** |
| zstd | Carquet | **137.7ms** | 150.5ms | **1.09x** |
| zstd | Arrow | **134.9ms** | 155.3ms | **1.15x** |
| lz4 | Carquet | **92.1ms** | 95.6ms | **1.04x** |
| lz4 | Arrow | **89.7ms** | 98.1ms | **1.09x** |

</div>

\* Zero-Copy-mmap-Pfad

</details>

<details>
<summary><b>Benchmarks selbst ausführen</b></summary>

<br/>

```bash
make benchmark                                    # Release build of the benchmark programs into build/
./build/benchmark_carquet                         # carquet alone
python3 benchmark/run_benchmark.py --skip-xlarge  # full comparison with PyArrow (+ Arrow C++ if built)
```

<div align="center">

| `run_benchmark.py`-Flag / -Variable | Wirkung |
|---|---|
| `--quick` | Nur die 10M-Zeilen-Konfigurationen (das, was das Diagramm zeigt) |
| `--skip-xlarge` | Die 100M-Zeilen-Konfigurationen überspringen (~2 GB pro Codec, 30+ Minuten) |
| `--fastparquet` | Zusätzlich fastparquet vergleichen (`pip install fastparquet pandas`) |
| `--polars` | Zusätzlich den eigenen Rust-Parquet-Reader / -Writer von polars vergleichen (`pip install polars`) |
| `--no-pyarrow` / `--no-arrow-cpp` | PyArrow / Arrow C++ überspringen |
| `--iterations N` | Feste Anzahl an Iterationen (Standard: automatisch, 11–51) |
| `--build-dir DIR` / `--temp-dir DIR` | Wo die Binaries liegen / wohin die Testdateien geschrieben werden |
| `CARQUET_BENCH_ZSTD_LEVEL=3` | ZSTD-Level (Standard 1) |

</div>

Benchmarken Sie immer einen **Release**-Build: `make benchmark` konfiguriert `build/` als Release neu; ein übrig gebliebenes Debug-`build/` lässt carquet um ein Mehrfaches langsamer erscheinen. `make native` + `-DCARQUET_BUILD_BENCHMARKS=ON` liefert die auf den Host abgestimmten Werte.

fastparquet verwendet dieselben Daten, dieselbe Row-Group-Größe, kein Dictionary und dasselbe zstd-Level sowie LZ4_RAW für `lz4`. Es kann weder BYTE_STREAM_SPLIT noch Page-Prüfsummen schreiben und liest über pandas, daher sind seine Werte weniger direkt vergleichbar als die von Arrow.

polars läuft mit seinem nativen Reader und Writer (nicht über PyArrow), denselben Daten, derselben Row-Group-Größe, demselben zstd-Level und LZ4_RAW für `lz4`. Es kann weder BYTE_STREAM_SPLIT noch Page-Prüfsummen schreiben und wählt das Dictionary-Encoding selbst (kein Schalter), wodurch seine Dateien auf diesen Daten kleiner werden.

Der Arrow-C++-Vergleich benötigt ein eigenes Target (Arrow C++ wird auf dem System gefunden oder über `CARQUET_ARROW_CPP_ROOT` angegeben):

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON \
  -DCARQUET_BUILD_BENCHMARKS=ON -DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON \
  [-DCARQUET_ARROW_CPP_ROOT=/path/to/arrow-prefix]
cmake --build build -j
```

Er verwendet die Low-Level-API `parquet::ParquetFileReader` (ohne Arrow-Table-Materialisierung) mit parallelen Row-Group-Readern.

Unter Windows (Developer PowerShell for VS) kann Arrow C++ von conda-forge kommen (`micromamba create -p C:\arrow -c conda-forge libarrow libparquet`); seine DLLs müssen beim Ausführen des Benchmarks im `PATH` liegen:

```powershell
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCARQUET_BUILD_BENCHMARKS=ON `
  -DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON -DCARQUET_ARROW_CPP_ROOT=C:/arrow/Library
cmake --build build
$env:PATH = "C:\arrow\Library\bin;$env:PATH"
python benchmark\run_benchmark.py --build-dir build --skip-xlarge --no-pyarrow
```

Halten Sie den Build- und den `--temp-dir`-Ordner aus dem Echtzeit-Virenscan heraus: Windows Defender scannt jede Datei, die der Benchmark schreibt.

Das README-Diagramm wird aus dem JSON erzeugt, das ein Lauf in `benchmark/` hinterlässt, in einer hellen und einer dunklen Variante, zwischen denen GitHub je nach Theme des Lesers wählt:

```bash
python3 benchmark/plot_readme_chart.py benchmark/bench_<version>_<os>_<arch>_<date>.json   # writes res/img/bench_10m_m3_{light,dark}.png
```

</details>

<br/>

---

<!-- ────────────────────────────  KURZER BLICK  ──────────────────────────── -->

## Carquet in dreißig Sekunden

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

> **Wo carquet hineinpasst.** Arrow bleibt die richtige Wahl, wenn Sie sein Ökosystem brauchen (Compute-Kernel, Datasets, Flight, Verschlüsselung). carquet ist für C- und C++-Programme, eingebettete Systeme und Sprachlaufzeiten gedacht, die Parquet lediglich schnell und korrekt lesen und schreiben müssen.

<br/>

---

<!-- ────────────────────────────  INSTALLATION  ──────────────────────────── -->

## Installation

Dieser Abschnitt installiert die **C-Bibliothek, ihre Header und die `carquet`-CLI**. Python-Nutzer brauchen die separate [Python-API](#python-api) (`pip install carquet`).

**Voraussetzungen:** ein C11-Compiler (GCC 4.9+, Clang 3.4+, MSVC 2015+) und CMake 3.16+ (oder [xmake](#aus-dem-quellcode-mit-xmake)). zstd, zlib, lz4 und brotli werden vom System verwendet, falls vorhanden, und andernfalls nachgeladen. OpenMP ist optional (paralleles Lesen und Schreiben).

<br/>

### Aus dem Quellcode, mit `make`

Der einfachste Weg. `make` ist ein dünner Wrapper, der CMake mit den richtigen Flags aufruft; es gibt nichts weiter zu konfigurieren:

```bash
git clone https://github.com/Vitruves/carquet.git && cd carquet
make                  # optimized Release build into build/
sudo make install     # library, headers, CLI, pkg-config and CMake package files
```

Danach linken Sie mit `-lcarquet`, `pkg-config --cflags --libs carquet` oder `find_package(carquet)` in CMake.

<br/>

### Aus dem Quellcode, mit `xmake`

[xmake](https://xmake.io) benötigt sonst nichts Installiertes: Es lädt fest versionierte zstd 1.5.6, zlib 1.4, lz4 1.10.0 und brotli 1.1.0 herunter und linkt sie statisch, sodass die `carquet`-CLI keine Laufzeitabhängigkeit hat.

```bash
git clone https://github.com/Vitruves/carquet.git && cd carquet
xmake                         # Release static library + CLI into build/<plat>/<arch>/release/
xmake install -o /usr/local   # headers, library and CLI (sudo if the prefix needs it)
```

Eine statische Installation enthält nur `libcarquet.a`, linken Sie die Codecs also selbst (`-lcarquet -lzstd -lz -llz4`, dazu `-lbrotlienc -lbrotlidec -lbrotlicommon`, sofern nicht mit `--brotli=n` gebaut). Konfigurieren Sie stattdessen mit `xmake f --shared=y`, um eine eigenständige Shared Library mit eingebauten Codecs zu erhalten.

<br/>

### In Ihrem CMake-Projekt

```cmake
include(FetchContent)
FetchContent_Declare(carquet GIT_REPOSITORY https://github.com/Vitruves/carquet.git GIT_TAG main)
FetchContent_MakeAvailable(carquet)
target_link_libraries(myapp PRIVATE carquet)
```

Oder, nach `make install`: `find_package(carquet REQUIRED)` und `carquet::carquet` linken.

<br/>

---

<!-- ────────────────────────────  BUILD-KONFIGURATIONEN  ──────────────────────────── -->

## Build-Konfigurationen

### `make`-Targets

Jedes Target konfiguriert `build/` und baut es parallel. `make help` gibt diese Liste aus. Zu jedem gibt es ein [xmake-Äquivalent](#xmake).

<div align="center">

| Befehl | Build-Typ | Was Sie erhalten |
|---|---|---|
| `make` / `make release` | Release | Statische Bibliothek + `carquet`-CLI (der Standard) |
| `make native` | Release | Dasselbe, mit `-march=native` auf diese CPU abgestimmt (das Binary läuft nicht auf älteren CPUs) |
| `make shared` | Release | Shared Library (`.so` / `.dylib` / `.dll`) statt einer statischen |
| `make benchmark` | Release | Benchmark-Programme (`benchmark_carquet`, Profiling-Targets) |
| `make debug` | Debug | Bibliothek + CLI mit Assertions, ohne Optimierung |
| `make test` | Debug | Testsuite, die anschließend mit `ctest` ausgeführt wird |
| `make dev` | Debug | Alles: Tests, Beispiele, Benchmarks, Interop-Programme |
| `make fuzz` | Debug | libFuzzer-Targets (mit einer Clang-Toolchain) |
| `make install` | Release | Baut `release` und installiert dann nach `PREFIX` |
| `make clean` | | Löscht das Build-Verzeichnis |

</div>

<div align="center">

| Variable | Standard | Beispiel |
|---|---|---|
| `PREFIX` | CMake-Standard (`/usr/local` unter Unix) | `make install PREFIX=$HOME/.local` |
| `BUILD` | `build` | `make test BUILD=build-debug` (mehrere Konfigurationen nebeneinander halten) |
| `JOBS` | Anzahl der Kerne | `make JOBS=4` |
| `CMAKE` | `cmake` | `make CMAKE=/opt/cmake/bin/cmake` |

</div>

Alle Targets teilen sich standardmäßig `build/`, sodass `make test` nach `make benchmark` es als Debug neu konfiguriert. Verwenden Sie `BUILD=`, um einen Release- und einen Debug-Baum getrennt zu halten.

<br/>

### CMake-Optionen

Für alles, was die `make`-Targets nicht abdecken, rufen Sie CMake direkt auf und übergeben Optionen mit `-DOPTION=ON|OFF`:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_BUILD_SHARED=ON -DCARQUET_ENABLE_BROTLI=OFF
cmake --build build -j
sudo cmake --install build --prefix /usr/local
```

<div align="center">

| Was gebaut wird | Standard | |
|---|:---:|---|
| `CARQUET_BUILD_CLI` | ON | `carquet`-Kommandozeilenwerkzeug |
| `CARQUET_BUILD_SHARED` | OFF | Shared statt statischer Bibliothek |
| `CARQUET_BUILD_DEV` | OFF | Schaltet Tests, Beispiele, Benchmarks und Interop auf einmal ein |
| `CARQUET_BUILD_TESTS` | = DEV | Testsuite (`ctest`) |
| `CARQUET_BUILD_EXAMPLES` | = DEV | Programme in [`examples/`](examples/) |
| `CARQUET_BUILD_BENCHMARKS` | = DEV | Benchmark- und Profiling-Programme |
| `CARQUET_BUILD_ARROW_CPP_BENCHMARK` | OFF | Arrow-C++-Vergleichsbenchmark (benötigt Arrow C++) |
| `CARQUET_BUILD_INTEROP` | = DEV | C- / Python-Round-Trip-Programme |
| `CARQUET_BUILD_FUZZ` | OFF | Fuzz-Targets |
| `CARQUET_BUILD_DOCS` | OFF | `docs`-Target, Doxygen-HTML nach `build/docs/html/` |

| Abhängigkeiten | Standard | |
|---|:---:|---|
| `CARQUET_ENABLE_BROTLI` | ON | BROTLI-Codec; OFF entfernt die libbrotli-Abhängigkeit |
| `CARQUET_BUNDLE_DEPS` | OFF | zstd / zlib / lz4 / brotli immer nachladen und statisch einbinden (portable Artefakte, Python-Wheel) |
| `CARQUET_ZSTD_TARGET`, `_ZLIB_`, `_LZ4_`, `_BROTLI_TARGET` | leer | Ein vorhandenes CMake-Target linken, statt zu suchen oder nachzuladen (siehe unten) |
| `CARQUET_ARROW_CPP_ROOT` | leer | Installationspräfix von Arrow C++ für den Vergleichsbenchmark |
| `CMAKE_DISABLE_FIND_PACKAGE_OpenMP` | OFF | ON baut ohne OpenMP (kein `libgomp` / `vcomp140.dll` zur Laufzeit) |

| CPU | Standard | |
|---|:---:|---|
| `CARQUET_NATIVE_ARCH` | OFF | `-march=native`; am schnellsten, Binary nur für den Host |
| `CARQUET_ENABLE_SSE`, `_AVX`, `_AVX2`, `_AVX512` | ON | x86-SIMD-Kernel, zur Laufzeit danach ausgewählt, was die CPU unterstützt |
| `CARQUET_ENABLE_NEON` | ON | ARM-NEON-Kernel |
| `CARQUET_ENABLE_SVE` | OFF | ARM-SVE-Kernel (experimentell) |
| `CARQUET_MSVC_HARDENING` | ON | MSVC: `/guard:cf` und `/CETCOMPAT` |

</div>

SIMD-Kernel werden einkompiliert und zur Laufzeit mit skalaren Fallbacks ausgewählt, sodass ein Standard-Build auf jeder CPU seiner Architektur läuft. Nur `CARQUET_NATIVE_ARCH` bindet das Binary an die Build-Maschine.

<br/>

### xmake

`xmake.lua` ist eine vollwertige Alternative zu CMake mit denselben Optionen und Standardwerten. Einmal mit `xmake f` konfigurieren, dann baut `xmake`; die Konfiguration bleibt bis zum nächsten `xmake f` erhalten.

<div align="center">

| `make`-Target | xmake-Äquivalent |
|---|---|
| `make` | `xmake` |
| `make native` | `xmake f --native_arch=y && xmake` |
| `make shared` | `xmake f --shared=y && xmake` |
| `make benchmark` | `xmake f --benchmarks=y && xmake` |
| `make debug` | `xmake f -m debug && xmake` |
| `make test` | `xmake f -m debug --tests=y && xmake && xmake test` |
| `make dev` | `xmake f -m debug --dev=y && xmake` |
| `make fuzz` | `xmake f -m debug --fuzz=y --toolchain=clang && xmake` |
| `make clean` | `xmake clean` (oder `xmake f -c`, um auch die Konfiguration zurückzusetzen) |
| `cmake --build build --target docs` | `xmake docs` (kein Konfigurations-Flag nötig) |

</div>

Die Optionen entsprechen denen von CMake, ohne das Präfix `CARQUET_`, in Kleinbuchstaben, als `--option=y|n`. Sie lassen sich frei kombinieren: `xmake f -m release --shared=y --brotli=n --avx512=n`.

<div align="center">

| xmake-Option | Standard | Beschreibung |
|--------------|---------|-------------|
| `--dev` | n | Tests, Beispiele, Benchmarks und Interop bauen |
| `--tests` / `--examples` / `--benchmarks` / `--interop` | n | Eine Gruppe einzeln bauen |
| `--arrow_cpp_benchmark` | n | Arrow-C++-Vergleichsbenchmark (verwendet Arrow aus dem xmake-Repo, falls verfügbar) |
| `--cli` | y | Das `carquet`-CLI-Tool bauen |
| `--shared` | n | Eine Shared Library statt einer statischen bauen |
| `--openmp` | y | Paralleles Lesen von Spalten mit OpenMP (automatisch deaktiviert, falls nicht verfügbar) |
| `--native_arch` | n | `-march=native` für maximale Leistung (Binary nur für den Host) |
| `--sse` / `--avx` / `--avx2` / `--avx512` / `--neon` | y | SIMD-Befehlssätze (automatisch erkannt) |
| `--sve` | n | ARM SVE (experimentell) |
| `--brotli` | y | BROTLI-Codec (libbrotli aus dem xmake-Repo) |
| `--fuzz` | n | Fuzz-Targets bauen (mit `--toolchain=clang`); dann `--fuzz_engine=libFuzzer\|AFL`, `--fuzz_sanitizers=y\|n` |

</div>

<details>
<summary><b>Windows</b> (MSVC, clang-cl, MinGW-w64)</summary>

<br/>

MSVC (x64 und Win32), clang-cl und MinGW-w64 GCC werden alle unterstützt, statisch oder shared. `CMakePresets.json` enthält die Presets `msvc`, `msvc-shared`, `clang-cl` und `mingw`:

```powershell
# from a Developer PowerShell for VS (or with MinGW-w64 gcc + ninja on PATH for "mingw")
cmake --preset msvc
cmake --build --preset msvc
ctest --preset msvc
```

- Wenn Sie von Hand konfigurieren, setzen Sie den C++-Compiler ebenso wie den C-Compiler (zum Beispiel `CC=clang-cl CXX=clang-cl`): Das nachgeladene zstd-Projekt aktiviert C++, und CMake wählt sonst einen nicht passenden Treiber.
- Mit `-DCARQUET_BUILD_SHARED=ON` erhalten Sie `carquet.dll` plus die Importbibliothek `carquet.lib` (`libcarquet.dll` / `libcarquet.dll.a` mit MinGW). Nutzer müssen `CARQUET_BUILD_SHARED` definieren; `find_package(carquet)` und `carquet.pc` erledigen das für Sie. Legen Sie die DLL neben Ihre ausführbare Datei oder in den `PATH`.
- Geben Sie Speicher, den carquet Ihnen übergibt, etwa die Ausgabe von `carquet_writer_get_buffer()`, mit `carquet_free()` frei und nicht mit `free()`. Die DLL und Ihr Programm können unterschiedliche Heaps der C-Laufzeit verwenden (Debug- vs. Release-CRT, `/MT`, MinGW vs. MSVC).
- Dateipfade dürfen UTF-8 sein (Nicht-ASCII-Namen funktionieren). Die `carquet`-CLI akzeptiert Unicode-Argumente und gibt UTF-8 auf der Konsole aus.
- MinGW-Builds mit OpenMP benötigen zur Laufzeit `libgomp-1.dll` und `libwinpthread-1.dll`. MSVC-Builds benötigen `vcomp140.dll` (aus dem Visual C++ Redistributable). Übergeben Sie `-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=ON` für einen Build ohne diese Abhängigkeit.
- Die Machine-Learning-Heuristik von Windows Defender stellt `test_encoding_roundtrip.exe` und `test_fuzz_regressions.exe` möglicherweise als `Trojan:Win64/CobaltStrike.RDH!MTB` unter Quarantäne (ein Fehlalarm wegen ihrer eingebetteten binären Testdaten; ctest meldet sie dann als „Not Run"). Schließen Sie das Build-Verzeichnis vom Scan aus, um die vollständige Suite auszuführen.

</details>

<details>
<summary><b>Eigene zstd / zlib / lz4 / brotli verwenden</b> (carquet in ein größeres CMake-Projekt einbetten)</summary>

<br/>

Wird carquet per `add_subdirectory()` oder `FetchContent` eingebunden, löst es jede Kompressionsbibliothek in dieser Reihenfolge auf:

1. Das durch `CARQUET_ZSTD_TARGET`, `CARQUET_ZLIB_TARGET`, `CARQUET_LZ4_TARGET` oder `CARQUET_BROTLI_TARGET` benannte Target, falls gesetzt.
2. Eine Systemkopie (wird übersprungen, wenn `CARQUET_BUNDLE_DEPS` ON ist).
3. Eine Kopie, die es selbst nachlädt und baut. Baut Ihr Projekt die Bibliothek bereits unter ihrem üblichen Target-Namen (`libzstd_static`, `zlibstatic`, `lz4_static`), verwendet carquet dieses Target weiter, statt eine zweite Kopie nachzuladen.

Ein Projekt, das diese Bibliotheken selbst mitliefert und keine Systemkopie hat, braucht daher keine Konfiguration. Um explizit zu wählen, benennen Sie das Target:

```cmake
set(CARQUET_ZLIB_TARGET zlibstatic)
add_subdirectory(carquet)
```

Ein benanntes Target wird unverändert gelinkt, ohne Suche und ohne Nachladen, und ein Name, der kein Target ist, ist ein Konfigurationsfehler und kein stiller Fallback. Das Target kann importiert oder im Baum gebaut, statisch oder shared sein und muss seine Include-Verzeichnisse als Usage Requirement mitführen.

</details>

<details>
<summary><b>Tests, Sanitizer und API-Dokumentation</b></summary>

<br/>

```bash
make test                                    # Debug build + ctest
make dev && ctest --test-dir build           # everything, including examples and interop programs
cmake -B build -DCARQUET_BUILD_DOCS=ON && cmake --build build --target docs   # Doxygen → build/docs/html/
```

[CONTRIBUTING.md](CONTRIBUTING.md) beschreibt das Gate pro Pull Request (Tests, Sanitizer, Fuzzing, Interop, Benchmarks).

</details>

<a id="ai-coding-agents"></a>

<details>
<summary><b>KI-Coding-Agenten</b> (ein fertiger Skill für Claude Code und Codex)</summary>

<br/>

[`.agents/skills/carquet/`](.agents/skills/carquet/) ist ein [Agent Skill](https://agentskills.io), der einem Coding-Agenten beibringt, wie carquet in Ihr Projekt integriert wird: Build-Anbindung (CMake, FetchContent, pkg-config, xmake), welche API wofür, vollständige Rezepte, die Regeln zu Levels / Null-Bitmaps / Lebensdauern, Typzuordnungen und die üblichen Fallstricke samt ihren Lösungen. Kopieren oder verlinken Sie ihn in Ihr Repository, und der Agent greift automatisch darauf zu:

```bash
mkdir -p .claude/skills && cp -r path/to/carquet/.agents/skills/carquet .claude/skills/   # Claude Code
mkdir -p .agents/skills && cp -r path/to/carquet/.agents/skills/carquet .agents/skills/   # Codex CLI
ln -s ../../third_party/carquet/.agents/skills/carquet .claude/skills/carquet               # vendored: stays in sync
```

Dieselben Dateien funktionieren auch als persönlicher Skill (`~/.claude/skills/`, `~/.codex/skills/`).

</details>

<br/>

---

<!-- ────────────────────────────  C-API  ──────────────────────────── -->

## C-API

Unten ein Schreib- und ein Lesebeispiel; alles Weitere behandelt das [Handbuch in `docs/`](docs/README.md).

<br/>

### Eine Parquet-Datei schreiben

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

### Eine Parquet-Datei lesen

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

### Weitere Rezepte

Alles jenseits von flachem Lesen/Schreiben steht im Handbuch; jede Seite verlinkt ein lauffähiges Beispiel:

<div align="center">

| Sie möchten… | Siehe |
|---|---|
| Nullable Spalten, Row Groups, Ausgabe in einen Puffer | [`docs/writing.md`](docs/writing.md) |
| Listen, Maps, Gruppen, Definition-/Repetition-Levels | [`docs/nested-data.md`](docs/nested-data.md) |
| Spaltenprojektion, Statistiken, Metadaten-Inspektion | [`docs/reading.md`](docs/reading.md) |
| Predicate Pushdown, Filterung auf Page-Ebene | [`docs/reading.md`](docs/reading.md) |
| Row Groups an eine bestehende Datei anhängen | [`docs/writing.md`](docs/writing.md#append-to-an-existing-file) |
| Kompression, eigene Codecs, Writer-Tuning | [`docs/writing.md`](docs/writing.md), [`docs/performance.md`](docs/performance.md) |
| mmap, Zero-Copy, Prebuffering, I/O-Coalescing | [`docs/performance.md`](docs/performance.md) |
| Fehlercodes und Hinweise zur Fehlerbehebung | [`docs/error-handling.md`](docs/error-handling.md) |

</div>

<br/>

### API-Referenz

Die vollständige API steht in [`include/carquet/carquet.h`](include/carquet/carquet.h). Wichtige Typen:

<div align="center">

| Typ | Zweck |
|------|---------|
| `carquet_reader_t` | Datei-Reader (öffnen aus Pfad, `FILE*` oder Speicherpuffer) |
| `carquet_writer_t` | Datei-Writer |
| `carquet_batch_reader_t` | High-Level-Batch-Iteration |
| `carquet_schema_t` | Schemadefinition und -introspektion |
| `carquet_error_t` | Ausführliche Fehlerinformation (Code, Meldung, Quellposition, Hinweis zur Behebung) |

</div>

Das [Handbuch](docs/README.md) erklärt, wann welche Schnittstelle zu verwenden ist; Quelltextaufbau und Architektur sind in [CONTRIBUTING.md](CONTRIBUTING.md) dokumentiert. Eine durchsuchbare HTML-Referenz lässt sich mit Doxygen aus den Headern erzeugen, siehe *Tests, Sanitizer und API-Dokumentation* unter [Build-Konfigurationen](#build-konfigurationen).

<details>
<summary><b>Beispielanwendung: MetricStore</b> (<code>mocklib/</code>)</summary>

<br/>

[`mocklib/`](mocklib/) ist **MetricStore**, eine vollständige, eigenständige Beispielanwendung auf Basis von carquet, nützlich sowohl als Referenz für den praktischen Einsatz als auch als End-to-End-Integrationstest der öffentlichen API.

Sie modelliert einen Zeitreihen-Telemetriespeicher (Ereignisse aufnehmen, dann inspizieren und abfragen) und nutzt dabei **117 der 143 öffentlichen Funktionen von carquet (~82 %)**:

- **Schreibseite**: Schemaaufbau mit logischen Typen, Encoding / Kompression / Bloom-Tuning pro Spalte, Page-Indexe und Statistiken
- **Leseseite**: Spaltenprojektion, Predicate Pushdown, Filterung auf Page-Ebene, Bloom-Mitgliedschaftstests, die Arrow-C-Data-Interface-Brücke, Rekonstruktion verschachtelter `LIST`-Spalten, Puffer-/Append-I/O und Metadaten-Introspektion

Sie linkt carquet so, wie es jedes nachgelagerte Projekt täte (`find_package(carquet)` / `add_subdirectory`), enthält einen selbstprüfenden Round-Trip-Test, und ihr `write_sample`-Binary erzeugt eine Datei zur externen Inspektion. Siehe [`mocklib/README.md`](mocklib/README.md).

</details>

<br/>

---

<!-- ────────────────────────────  PYTHON-API  ──────────────────────────── -->

## Python-API

`carquet` auf PyPI ist kein Wrapper um die obigen C-Aufrufe: Es ist eine **pythonische API im Stil von `pyarrow.parquet`** (`read_table`, `ParquetFile`, `ParquetWriter`, DNF-`filters`) auf Basis derselben Engine. Das Wheel bündelt die C-Bibliothek mit allen einkompilierten Codecs; die Daten gelangen über das Arrow C Data Interface nach Python, sodass Dekodierung, Dekompression und I/O in C ohne GIL laufen. **Sie hat keine Abhängigkeiten**: Lesevorgänge liefern eine `carquet.Table`, die sich selbst in Python-Werte, NumPy und pandas umwandelt und die polars, DuckDB und PyArrow ohne Kopie übernehmen.

```bash
pip install carquet              # no dependencies: no PyArrow, no NumPy needed
pip install "carquet[pandas]"    # with pandas + NumPy, for read_pandas() / to_pandas()
```

```python
import carquet as cq

table = cq.read_table("data.parquet", columns=["id", "price"], filters=[("price", ">", 100)])
df    = cq.read_pandas("data.parquet")                          # pandas.DataFrame, index restored
cq.write_table(table, "out.parquet", compression="zstd")        # also DataFrame, polars, dict, batches

with cq.ParquetFile("data.parquet") as pf:
    for batch in pf.iter_batches(batch_size=50_000, columns=["id"]):
        ...
    pf.metadata.row_group(0).column("price").statistics       # min / max / null_count
```

<div align="center">

| Sie möchten… | Verwenden |
|---|---|
| Eine Datei, ein Verzeichnis, Bytes oder ein Dateiobjekt lesen | `cq.read_table`, `cq.read_pandas` |
| Große Dateien streamen | `ParquetFile.iter_batches`, `cq.iter_batches` |
| Mit Statistiken, Page-Indexen, Bloom-Filtern prunen | `filters=` (DNF, exakte Ergebnisse), `ParquetFile.bloom_filter` |
| Tabellen, DataFrames, polars, DuckDB-Ergebnisse schreiben | `cq.write_table`, `cq.ParquetWriter` (inkrementell, `append=True`) |
| Daten an pandas, NumPy, polars, DuckDB, PyArrow übergeben | `t.to_pandas()`, `t.column("x").to_numpy()`, `polars.DataFrame(t)`, `pyarrow.table(t)` |
| Eine Datei inspizieren, ohne ihre Daten zu lesen | `cq.read_metadata`, `cq.read_parquet_schema`, `cq.info`, `cq.validate` |

</div>

Typen durchlaufen den Round-Trip mit `pyarrow.parquet` exakt (Zeitzonen, Kategorien, Dezimalzahlen, Large-Typen, verschachtelte Daten). Die vollständige Referenz steht in [`python/README.md`](python/README.md).

<br/>

---

<!-- ────────────────────────────  CLI-TOOL  ──────────────────────────── -->

## CLI-Tool

Carquet wird mit einem Kommandozeilenwerkzeug zum Inspizieren von Parquet-Dateien und zum Generieren von C-Reader-Code ausgeliefert. Es wird standardmäßig zusammen mit der Bibliothek gebaut und installiert.

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

Typische Verwendung:

```bash
carquet schema data.parquet
carquet head -n 20 data.parquet
carquet stat data.parquet
carquet validate data.parquet
```

<br/>

### Zeilenfilterung

`cat`, `count`, `head`, `tail`, `sample` und `export` akzeptieren `-p / --filter EXPR`, um ein Zeilenprädikat bis auf die Page-Ebene durchzureichen: Nur Pages, deren Min/Max im Column Index zum Prädikat passen können, werden dekomprimiert.

```bash
carquet cat -p "price > 100 AND status = 'active'" data.parquet
carquet count --filter "id >= 1000" data.parquet
carquet export --filter "ts IS NOT NULL" -c id,ts data.parquet
```

Die Grammatik lautet `column OP value [AND column OP value]...` mit `OP` ∈ {`=`, `==`, `!=`, `<>`, `<`, `<=`, `>`, `>=`}, dazu `column IS NULL` / `column IS NOT NULL`.

> Die Filterung setzt voraus, dass die Datei einen Page-Index hat (`write_page_index = true`).

<br/>

### Code-Generierung

Erzeugen Sie aus dem Schema einer beliebigen Parquet-Datei einen vollständigen, kompilierbaren C-Reader:

```bash
carquet codegen -f data.parquet -o reader.c
# Generated: reader.c
# Compile:   clang -o reader reader.c -I.../include -L.../build -lcarquet ...

./reader                    # reads data.parquet (embedded as default)
./reader other.parquet      # override with a different file
```

<details>
<summary><b>codegen-Optionen</b></summary>

<br/>

<div align="center">

| Flag | Beschreibung |
|------|-------------|
| `-f`, `--file FILE` | Parquet-Datei, deren Schema inspiziert wird |
| `-o`, `--output FILE` | Ausgabe-Quelldatei (Standard: stdout) |
| `--mmap` | Memory-Mapped I/O im generierten Code verwenden |
| `--skeleton` | Ein leeres `process_batch` für eigene Logik generieren |
| `-c`, `--columns COLS` | Kommagetrennter Spaltenfilter |
| `-b`, `--batch-size N` | Batch-Größe (Standard: 1024) |

</div>

</details>

<br/>

---

<!-- ────────────────────────────  UNTERSTÜTZTE FUNKTIONEN  ──────────────────────────── -->

## Unterstützte Parquet-Funktionen

Was die Bibliothek heute implementiert, nach Bereichen.

<br/>

### Format

<div align="center">

| Funktion | Status |
|---------|--------|
| Physische Typen | ✓  Alle 8 (BOOLEAN bis FIXED_LEN_BYTE_ARRAY) |
| Logische Typen | ✓  STRING, DATE, TIME, TIMESTAMP, DECIMAL, UUID, JSON, INTERVAL, FLOAT16, VARIANT, GEOMETRY, GEOGRAPHY |
| Encodings | ✓  PLAIN, RLE, DICTIONARY, DELTA_BINARY_PACKED, DELTA_LENGTH_BYTE_ARRAY, DELTA_BYTE_ARRAY, BYTE_STREAM_SPLIT (Lesen + Schreiben) |
| Data-Page-Versionen | ✓  V1 (Standard) und V2 (Lesen + Schreiben) |
| Kompression | ✓  UNCOMPRESSED, SNAPPY, GZIP, LZ4 (Hadoop-Framing, Codec 5), LZ4_RAW (Codec 7), ZSTD, BROTLI (`CARQUET_ENABLE_BROTLI`, Standard ON); LZO über `carquet_register_codec()` |
| Verschachtelte Schemata | ✓  Gruppen, Listen, Maps mit Definition-/Repetition-Levels; automatisches Shredding einstufiger LIST/MAP beim Schreiben (`carquet_writer_write_list_column`) und List-Rekonstruktion beim Lesen; Emission eines verschachtelten `ARROW:schema` |
| Verschlüsselung | ✗ Nicht unterstützt; eine Datei mit verschlüsseltem Footer (`PARE`) wird als `CARQUET_ERROR_NOT_IMPLEMENTED` mit einer expliziten Meldung gemeldet |

</div>

<details>
<summary><b>Indizierung &amp; Pruning</b></summary>

<br/>

<div align="center">

| Funktion | Status |
|---------|--------|
| Bloom-Filter | ✓  Lesen, Schreiben und Abfragen (`carquet_bloom_filter_check_*`); werden automatisch für das Row-Group-Pruning herangezogen |
| Page-Indexe | ✓  Column Index + Offset Index (Lesen + Schreiben + Zugriff auf Statistiken pro Page) |
| Statistiken | ✓  Min/Max/Null-Anzahl pro Column Chunk; exakter `distinct_count` für Dictionary-Spalten; Parquet-2.9-SizeStatistics (unkodierte Byte-Array-Bytes + Level-Histogramme) |
| Predicate Pushdown | ✓  Automatisches Row-Group-Pruning über Statistiken + Bloom-Filter (kein Callback nötig); Filterung auf Page-Ebene über den Column Index (`carquet_batch_reader_set_page_filter`) |
| Geostatistiken | ✓  Bounding Box und Geometrietypen pro Spalte für GEOMETRY/GEOGRAPHY (`carquet_reader_geospatial_statistics`) |
| Spaltenprojektion | ✓  Nur ausgewählte Spalten lesen |
| Zeilenfenster | ✓  Einen `[offset, limit)`-Ausschnitt der Zeilen lesen und Pages dabei vollständig überspringen (`carquet_batch_reader_set_row_range`) |

</div>

</details>

<details>
<summary><b>I/O &amp; Erweiterbarkeit</b></summary>

<br/>

<div align="center">

| Funktion | Status |
|---------|--------|
| Append | ✓  Row Groups an eine bestehende Datei anhängen (`carquet_writer_open_append`) |
| Eigene Codecs | ✓  Eigene Compress-/Decompress-Implementierung pro Codec-Slot registrieren (`carquet_register_codec`) |
| Key-Value-Metadaten | ✓  Beliebige Footer-Metadaten lesen und schreiben |
| Metadaten pro Feld | ✓  Arrow `Field.custom_metadata` (Variablenbezeichnungen/-beschreibungen) über `ARROW:schema` (Lesen + Schreiben) |
| Optionen pro Spalte | ✓  Encoding, Kompression, Statistiken, Bloom-Filter pro Spalte |
| Puffer-Writer | ✓  Parquet in einen Speicherpuffer schreiben |
| CRC32 | ✓  Verifikation auf Page-Ebene (hardwarebeschleunigt auf ARM) |
| Memory-Mapped I/O | ✓  Zero-Copy-Lesen für unkomprimierte PLAIN-Daten |
| Dictionary-erhaltendes Lesen | ✓  RLE-Indizes + Dictionary statt materialisierter Werte zurückgeben (`preserve_dictionaries`, `carquet_row_batch_column_dictionary`) |
| I/O-Coalescing | ✓  Mehrspaltiges Lesen in einem einzigen I/O vorpuffern |
| Spekulativer Footer | ✓  Öffnen der Datei mit einem einzigen I/O bei den meisten Dateien |
| Paralleles Lesen mit OpenMP | ✓  Sofern verfügbar |
| Paralleles Schreiben | ✓  Paralleles Kodieren/Komprimieren auf Page-Ebene über Spalten hinweg (OpenMP) + Row-Group-I/O im Hintergrund (`async_io`) |

</div>

</details>

<br/>

---

<!-- ────────────────────────────  INTEROPERABILITÄT  ──────────────────────────── -->

## Interoperabilität

Carquet-Dateien sind vollständig kompatibel mit PyArrow, DuckDB, Spark und jedem anderen Parquet-Reader.

```python
import pyarrow.parquet as pq
table = pq.read_table("carquet_output.parquet")  # just works
```

```sql
-- DuckDB
SELECT * FROM read_parquet('carquet_output.parquet');
```

<details>
<summary><b>Interop-Testsuite</b></summary>

<br/>

Bidirektionale Round-Trips (C schreibt / Python liest, und umgekehrt):

```bash
cmake -B build -DCARQUET_BUILD_INTEROP=ON && cmake --build build
python3 interop/run_interop.py
```

Große Snappy-Dateien mit vielen Row Groups (viele Row Groups, breite `BYTE_ARRAY`-Spalten, Page-Offsets jenseits von 2 GiB) haben einen eigenen Generator, der die Datei mit PyArrow erzeugt und prüft, dass carquet jede Page dekodiert:

```bash
CARQUET_BIN=build/carquet python3 interop/snappy_multi_rowgroup.py --rows 150000 --row-groups 92
```

</details>

<br/>

---

<!-- ────────────────────────────  MITWIRKEN & LIZENZ  ──────────────────────────── -->

## Mitwirken, Zitieren, Lizenz

- **Mitwirken**: [CONTRIBUTING.md](CONTRIBUTING.md) behandelt Einrichtung, Code-Stil und das Gate pro Pull Request (Tests, Sanitizer, Fuzzing, Interop, Benchmarks). Release Notes stehen in [CHANGELOG.md](CHANGELOG.md).
- **Zitieren**: eine [`CITATION.cff`](CITATION.cff) liegt bei; GitHubs Schaltfläche „Cite this repository" formatiert sie als BibTeX oder APA.
- **Lizenz**: MIT, siehe [LICENSE](LICENSE).
