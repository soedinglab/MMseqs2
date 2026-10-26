<h1 align="center">
  <img src="res/img/carquet_banner.png" alt="Carquet" width="840" />
</h1>

<p align="center">
  Eine schnelle, reine C-Bibliothek zum Lesen und Schreiben von Apache-Parquet-Dateien.
</p>

<p align="center">
  <a href="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml"><img src="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml/badge.svg" alt="Build" /></a>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue" alt="Plattform" />
  <img src="https://img.shields.io/badge/C-C11-blue" alt="C-Standard" />
  <a href="https://github.com/Vitruves/carquet/releases"><img src="https://img.shields.io/github/v/release/Vitruves/carquet?color=blue" alt="Release" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="Lizenz" /></a>
  <a href="https://github.com/Vitruves/carquet/stargazers"><img src="https://img.shields.io/github/stars/Vitruves/carquet?style=flat&color=blue" alt="Stars" /></a>
  <a href="https://github.com/Vitruves/carquet/graphs/traffic"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/Vitruves/carquet/traffic/clones.json" alt="Clones" /></a>
  <br/>
  <img src="https://img.shields.io/badge/SIMD-SSE4.2%20%C2%B7%20AVX2%20%C2%B7%20AVX--512%20%C2%B7%20NEON%20%C2%B7%20SVE-555" alt="SIMD" />
  <img src="https://img.shields.io/badge/deps-zstd%20%C2%B7%20zlib%20%C2%B7%20lz4%20%C2%B7%20brotli-555" alt="Abhängigkeiten" />
</p>

<p align="center">
  <a href="#highlights">Highlights</a> &nbsp;·&nbsp;
  <a href="#installation">Installation</a> &nbsp;·&nbsp;
  <a href="#c-api">C-API</a> &nbsp;·&nbsp;
  <a href="#cli-tool">CLI-Tool</a> &nbsp;·&nbsp;
  <a href="#benchmarks">Benchmarks</a> &nbsp;·&nbsp;
  <a href="#unterstützte-parquet-funktionen">Funktionsumfang</a> &nbsp;·&nbsp;
  <a href="#interoperabilität">Interoperabilität</a> &nbsp;·&nbsp;
  <a href="#ready-for-ai-coding-agents">KI-Agenten</a> &nbsp;·&nbsp;
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

**Klein und eigenständig**

- **Reines C11** mit vier externen Abhängigkeiten (zstd, zlib, lz4, brotli), die CMake bei Bedarf automatisch nachlädt; Brotli lässt sich abschalten
- **~500 KB Bibliothek** (gestrippter Release-Build, Codecs separat gelinkt) gegenüber ~50 MB+ bei Arrow
- **Integrierte CLI** zur Dateiinspektion (`schema`, `info`, `head`, `tail`, `stat`, …) und zur C-Code-Generierung (`codegen`)

**Schnell**

- **39x schnelleres Lesen** als Arrow C++ bei unkomprimierten Daten (mmap Zero-Copy), **100x schneller** als PyArrow
- **1.3–2.1x schnelleres komprimiertes Lesen** als Arrow C++ bei derselben Datei (Cross-Read-Benchmark, Apple M3; 1.6–2.6x auf x86)
- **Schreibt 1.8–5.7x schneller** als Arrow C++ auf Apple M3 (1.0–1.6x im älteren x86-Lauf)
- Liest 10 Mio. unkomprimierte Zeilen in **0.36 ms** (mmap Zero-Copy auf Apple M3)
- SIMD-optimiert (SSE4.2, AVX2, AVX-512, NEON, SVE) mit Laufzeiterkennung und skalaren Fallbacks

**Vollständig und kompatibel**

- Breite Parquet-Abdeckung: alle physischen Typen und Encodings, Snappy / GZIP / LZ4 / ZSTD / Brotli, verschachtelte Schemata, Bloom-Filter, Page-Indexe ([vollständige Matrix](#unterstützte-parquet-funktionen))
- Out of the box kompatibel mit PyArrow, DuckDB und Spark
- [Arrow C Data Interface](docs/reading.md#export-to-arrow-c-data-interface)-Brücke (`carquet_arrow_export_*` / `carquet_arrow_import_*` / `carquet_writer_write_arrow` / `carquet_reader_read_arrow`) für abhängigkeitsfreien, kopierarmen Austausch mit dem Arrow-Ökosystem, einschließlich **verschachteltem Lesen/Schreiben in beliebiger Tiefe** (struct / list / map auf jeder Ebene)

<a id="ready-for-ai-coding-agents"></a>

**Bereit für KI-Coding-Agenten**

- [`.agents/skills/carquet/`](.agents/skills/carquet/) ist ein [Agent Skill](https://agentskills.io) (Claude Code, Codex CLI), der einem Agenten beibringt, carquet in Ihr Projekt zu integrieren: Build-Anbindung, welche API für welche Aufgabe, vollständige Rezepte, die Regeln zu Levels / Null-Bitmaps / Lebensdauern, Typzuordnungen und die üblichen Fallstricke. Kopieren Sie ihn in Ihr Repository, und der Agent nutzt ihn von selbst:

  ```bash
  cp -r path/to/carquet/.agents/skills/carquet .claude/skills/   # Claude Code
  cp -r path/to/carquet/.agents/skills/carquet .agents/skills/   # Codex CLI
  ```

<br/>

### Carquet in dreißig Sekunden

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
printf("%lld rows\n", (long long)carquet_reader_num_rows(r));   /* 3 Zeilen */
carquet_reader_close(r);
carquet_schema_free(schema);
```

Vollständige Schreib- und Leseschleifen, mit Batches und nullable Spalten, finden Sie unter [C-API](#c-api).

<br/>

### Wo carquet hineinpasst

carquet ist die kleine, einbettbare Option: eine ~500 KB große C11-Bibliothek mit vier automatisch nachgeladenen Abhängigkeiten, die in unter einer Minute baut, während Arrow C++ 50 MB+ und ein PyArrow-Wheel rund 100 MB wiegen. Im Benchmark mit 10 Mio. Zeilen unten schreibt es komprimierte Dateien 1.8–5.7x schneller als Arrow C++ und liest sie aus ein und derselben Datei 1.3–2.1x schneller; unkomprimierte Daten werden in 0.36 ms direkt aus dem mmap bedient. Arrow bleibt die richtige Wahl, wenn Sie sein Ökosystem brauchen (Compute-Kernel, Datasets, Flight, Verschlüsselung); carquet ist für die C- und C++-Programme, eingebetteten Systeme und Sprachlaufzeiten gedacht, die Parquet lediglich schnell und korrekt lesen und schreiben müssen.

<br/>

---

<!-- ────────────────────────────  INSTALLATION  ──────────────────────────── -->

## Installation

Carquet wird mit CMake (primär) oder [xmake](https://xmake.io) gebaut, unter Linux, macOS und Windows.

### Voraussetzungen

- C11-Compiler (GCC 4.9+, Clang 3.4+, MSVC 2015+)
- CMake 3.16+ (oder [xmake](https://xmake.io), siehe [Bauen mit xmake](#bauen-mit-xmake))
- zstd, zlib, lz4, brotli (werden bei Bedarf automatisch nachgeladen; Brotli optional, `-DCARQUET_ENABLE_BROTLI=OFF`)
- OpenMP (optional, für paralleles Lesen und Schreiben von Spalten)

<br/>

### Schnellstart (make)

Der `make`-Wrapper steuert einen optimierten CMake-Build und eine Installation nach `/usr/local`:

```bash
git clone https://github.com/Vitruves/carquet.git
cd carquet
make                              # optimierter Build (`make help` listet alle Ziele auf)
sudo make install                 # Installation nach /usr/local (überschreibbar mit PREFIX=/opt/carquet)
```

<br/>

### Bauen & Installieren (CMake)

Rufen Sie CMake direkt auf, wenn Sie bestimmte Build-Optionen oder ein eigenes Präfix benötigen:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release   # Optionen ergänzen, z. B. -DCARQUET_BUILD_SHARED=ON
cmake --build build -j$(nproc)
sudo cmake --install build --prefix /usr/local
```

Beide Wege installieren:

- `libcarquet.a` (oder `.so` / `.dylib` mit `-DCARQUET_BUILD_SHARED=ON`)
- die Header unter `include/carquet/`
- das `carquet`-CLI-Binary
- `carquet.pc` (pkg-config) und die CMake-Package-Konfiguration für `find_package(carquet)`

Nach der Installation linken Sie mit `-lcarquet` oder ermitteln die Flags über `pkg-config --cflags --libs carquet`.

<br/>

### Python (pip)

Das `carquet`-Wheel bündelt die Bibliothek und stellt sie über das Arrow C Data Interface bereit, sodass Tabellen ohne Kopien der Parquet-Dekodierlogik von und nach PyArrow fließen:

```bash
pip install carquet
```

```python
import carquet
table = carquet.read_table("data.parquet")                 # pyarrow.Table
carquet.write_table(table, "out.parquet", compression="zstd")
print(carquet.info("out.parquet").num_rows)
```

Das Paket liegt in [`python/`](python/) (scikit-build-core, ctypes, PyArrow zur Laufzeit optional); API und Wheel-Build sind in seiner [README](python/README.md) beschrieben.

<br/>

### Build-Optionen

Übergeben Sie diese im CMake-Konfigurationsschritt (`-DOPTION=ON|OFF`):

<div align="center">

| Option | Standard | Beschreibung |
|--------|---------|-------------|
| `CARQUET_BUILD_DEV` | OFF | Alles bauen (Tests, Beispiele, Benchmarks, Interop) |
| `CARQUET_BUILD_TESTS` | OFF | Nur die Testsuite bauen |
| `CARQUET_BUILD_CLI` | ON | Das `carquet`-CLI-Tool bauen |
| `CARQUET_BUILD_SHARED` | OFF | Eine Shared Library statt einer statischen bauen |
| `CARQUET_NATIVE_ARCH` | OFF | `-march=native` für maximale Leistung (Binary nur für den Host) |
| `CARQUET_ENABLE_SVE` | OFF | ARM SVE (experimentell) |
| `CARQUET_ENABLE_BROTLI` | ON | BROTLI-Codec (libbrotli gefunden oder nachgeladen; OFF entfernt die Abhängigkeit) |

</div>

Alle x86-SIMD-Erweiterungen (SSE, AVX, AVX2, AVX-512) und ARM NEON werden automatisch erkannt und sind standardmäßig aktiviert.

<details>
<summary><b>Alle Build-Optionen</b></summary>

<br/>

<div align="center">

| Option | Standard | Beschreibung |
|--------|---------|-------------|
| `CARQUET_BUILD_EXAMPLES` | OFF | Beispielprogramme bauen |
| `CARQUET_BUILD_BENCHMARKS` | OFF | Benchmark- und Profiling-Programme bauen |
| `CARQUET_BUILD_ARROW_CPP_BENCHMARK` | OFF | Optionaler Arrow-C++-Vergleichsbenchmark |
| `CARQUET_BUILD_INTEROP` | OFF | Interoperabilitätstests bauen |
| `CARQUET_BUILD_FUZZ` | OFF | Fuzz-Targets bauen |
| `CARQUET_BUILD_DOCS` | OFF | Ein `docs`-Target hinzufügen (benötigt Doxygen) |
| `CARQUET_BUNDLE_DEPS` | OFF | zstd/zlib/lz4/brotli immer nachladen und statisch einbinden (portable Release-Artefakte) |
| `CARQUET_ZSTD_TARGET` / `CARQUET_ZLIB_TARGET` / `CARQUET_LZ4_TARGET` / `CARQUET_BROTLI_TARGET` | *(leer)* | Vorhandene(s) CMake-Target(s), die für diese Bibliothek gelinkt werden, statt sie zu suchen oder nachzuladen |
| `CARQUET_ARROW_CPP_ROOT` | *(leer)* | Pfad zu einer eigenen Arrow-C++-Installation für den Vergleichsbenchmark |
| `CARQUET_ENABLE_SSE` | ON | SSE-Optimierungen (x86, automatisch erkannt) |
| `CARQUET_ENABLE_AVX` | ON | AVX-Optimierungen (x86, automatisch erkannt) |
| `CARQUET_ENABLE_AVX2` | ON | AVX2-Optimierungen (x86, automatisch erkannt) |
| `CARQUET_ENABLE_AVX512` | ON | AVX-512-Optimierungen (x86, automatisch erkannt) |
| `CARQUET_ENABLE_NEON` | ON | NEON-Optimierungen (ARM, automatisch erkannt) |

</div>

</details>

<details>
<summary><b>Eigene zstd / zlib / lz4 verwenden</b> (carquet in ein größeres CMake-Projekt einbetten)</summary>

<br/>

Wird carquet per `add_subdirectory()` oder `FetchContent` eingebunden, löst es jede Kompressionsbibliothek in dieser Reihenfolge auf:

1. Das durch `CARQUET_ZSTD_TARGET`, `CARQUET_ZLIB_TARGET` oder `CARQUET_LZ4_TARGET` benannte Target, falls gesetzt.
2. Eine Systemkopie (wird übersprungen, wenn `CARQUET_BUNDLE_DEPS` ON ist).
3. Eine Kopie, die es selbst nachlädt und baut. Baut Ihr Projekt die Bibliothek bereits unter ihrem üblichen Target-Namen (`libzstd_static`, `zlibstatic`, `lz4_static`), verwendet carquet dieses Target weiter, statt eine zweite Kopie nachzuladen.

Ein Projekt, das diese Bibliotheken selbst mitliefert und keine Systemkopie hat, braucht daher keine Konfiguration. Um explizit zu wählen, etwa Ihre eigene zlib statt der des Systems oder Ihren Shared-Build statt des statischen, benennen Sie das Target:

```cmake
set(CARQUET_ZLIB_TARGET zlibstatic)
add_subdirectory(carquet)
```

Dasselbe geht auch, ohne etwas zu bearbeiten: `cmake -DCARQUET_ZLIB_TARGET=zlib ...`. Die drei Variablen sind voneinander unabhängig. Ein benanntes Target wird unverändert gelinkt, ohne Suche und ohne Nachladen, und ein Name, der kein Target ist, ist ein Konfigurationsfehler und kein stiller Fallback. Das Target kann importiert oder im Baum gebaut, statisch oder shared sein und muss seine Include-Verzeichnisse als Usage Requirement mitführen.

</details>

<details>
<summary><b>Entwicklungs-Build</b> (Tests, Beispiele, Benchmarks)</summary>

<br/>

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCARQUET_BUILD_DEV=ON
cmake --build build -j$(nproc)
cd build && ctest --output-on-failure
```

[CONTRIBUTING.md](CONTRIBUTING.md) beschreibt das Gate pro Pull Request (Tests, Sanitizer, Fuzzing, Interop, Benchmarks).

</details>

<details>
<summary><b>API-Dokumentation</b> (Doxygen)</summary>

<br/>

Die öffentlichen Header sind mit [Doxygen](https://www.doxygen.nl)-Kommentaren annotiert. Zum Erzeugen der HTML-Referenz muss `doxygen` installiert sein (optional dazu `graphviz` für Include-/Abhängigkeitsdiagramme); für einen normalen Build wird keines von beiden benötigt.

```bash
# CMake — das Target beim Konfigurieren aktivieren, dann bauen
cmake -B build -DCARQUET_BUILD_DOCS=ON
cmake --build build --target docs

# xmake — eigenständiger Task, kein Konfigurations-Flag nötig
xmake docs
```

Beide schreiben nach `build/docs/html/index.html`. Ist Doxygen nicht installiert, meldet das Target bzw. der Task das lediglich und tut nichts.

</details>

<details>
<summary><b>KI-Coding-Agenten</b> (Claude Code, Codex: ein fertiger Skill)</summary>

<br/>

[`.agents/skills/carquet/`](.agents/skills/carquet/) ist ein [Agent Skill](https://agentskills.io), der einem Coding-Agenten beibringt, wie carquet in Ihr Projekt integriert wird: Build-Anbindung (CMake, FetchContent, pkg-config, xmake), welche API wofür, vollständige Rezepte, die Regeln zu Levels / Null-Bitmaps / Lebensdauern, Typzuordnungen und die üblichen Fallstricke samt ihren Lösungen. Kopieren oder verlinken Sie ihn in Ihr eigenes Repository, und der Agent greift automatisch darauf zu:

```bash
# Claude Code
mkdir -p .claude/skills && cp -r path/to/carquet/.agents/skills/carquet .claude/skills/
# Codex CLI
mkdir -p .agents/skills && cp -r path/to/carquet/.agents/skills/carquet .agents/skills/
# mitgeliefertes (vendored) carquet: ein Symlink hält ihn synchron
ln -s ../../third_party/carquet/.agents/skills/carquet .claude/skills/carquet
```

Dieselben Dateien funktionieren auch als persönlicher Skill (`~/.claude/skills/`, `~/.codex/skills/`).

</details>

<br/>

### Bauen mit xmake

Als Alternative zu CMake liegt ein [xmake](https://xmake.io)-Build (`xmake.lua`) bei, mit denselben Optionen und Standardwerten; zstd/zlib/lz4 werden statisch gelinkt, sodass die Binaries eigenständig sind.

```bash
xmake                    # statische Bibliothek + `carquet`-CLI bauen (Release)
xmake f --dev=y && xmake # Tests, Beispiele, Benchmarks, Interop hinzufügen
xmake test               # Testsuite ausführen
```

<details>
<summary><b>Konfigurationsoptionen</b> (<code>xmake f --option=y|n</code>)</summary>

<br/>

Die Optionen entsprechen denen von CMake, ohne das Präfix `CARQUET_` und in Kleinbuchstaben:

<div align="center">

| xmake-Option | Standard | Beschreibung |
|--------------|---------|-------------|
| `--dev` | n | Tests, Beispiele, Benchmarks und Interop bauen |
| `--tests` / `--examples` / `--benchmarks` / `--interop` | n | Eine Gruppe einzeln bauen |
| `--cli` | y | Das `carquet`-CLI-Tool bauen |
| `--shared` | n | Eine Shared Library statt einer statischen bauen |
| `--openmp` | y | Paralleles Lesen von Spalten mit OpenMP (automatisch deaktiviert, falls nicht verfügbar) |
| `--native_arch` | n | `-march=native` für maximale Leistung (Binary nur für den Host) |
| `--sse` / `--avx` / `--avx2` / `--avx512` / `--neon` | y | SIMD-Befehlssätze (automatisch erkannt) |
| `--sve` | n | ARM SVE (experimentell) |
| `--brotli` | y | BROTLI-Codec (libbrotli aus dem xmake-Repo) |
| `--fuzz` | n | Fuzz-Targets bauen (mit `--toolchain=clang`) |

</div>

```bash
xmake f -m release --shared=y            # Shared Library
xmake f --dev=y --avx512=n && xmake      # Dev-Build, AVX-512 deaktiviert
```

</details>

<br/>

---

<!-- ────────────────────────────  C API  ──────────────────────────── -->

## C-API

Unten ein Schreib- und ein Lesebeispiel; alles Weitere behandelt das [Handbuch in `docs/`](docs/README.md).

### Eine Parquet-Datei schreiben

```c
#include <carquet/carquet.h>

int main(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    // Schema definieren
    carquet_schema_t* schema = carquet_schema_create(&err);
    carquet_schema_add_column(schema, "id",    CARQUET_PHYSICAL_INT64,  NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    // Writer konfigurieren
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;

    // Schreiben
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

    // Mit mmap öffnen für die beste Leseleistung
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = true;

    carquet_reader_t* r = carquet_reader_open("output.parquet", &opts, &err);
    if (!r) { printf("Error: %s\n", err.message); return 1; }

    printf("Rows: %lld, Columns: %d\n",
           (long long)carquet_reader_num_rows(r),
           carquet_reader_num_columns(r));

    // Batch-Reader für effiziente Iteration
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
        // ids[0..n-1] verarbeiten ...
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

Das [Handbuch](docs/README.md) erklärt, wann welche Schnittstelle zu verwenden ist; Quelltextaufbau und Architektur sind in [CONTRIBUTING.md](CONTRIBUTING.md) dokumentiert. Eine durchsuchbare HTML-Referenz lässt sich mit Doxygen aus den Headern erzeugen, siehe den ausklappbaren Abschnitt *API-Dokumentation* unter [Build-Optionen](#build-optionen).

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

<!-- ────────────────────────────  CLI TOOL  ──────────────────────────── -->

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

./reader                    # liest data.parquet (als Standard eingebettet)
./reader other.parquet      # mit einer anderen Datei überschreiben
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

<!-- ────────────────────────────  BENCHMARKS  ──────────────────────────── -->

## Benchmarks

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="res/img/bench_10m_m3_dark.png">
    <img src="res/img/bench_10m_m3_light.png" alt="Schreib- und Lesezeit für 10 Millionen Zeilen auf Apple M3: carquet vs. Arrow C++ vs. PyArrow, pro Codec" width="900">
  </picture>
</p>

Bei 10 Mio. Zeilen (der repräsentativsten Größe); höheres Verhältnis = Carquet schneller.

ARM (Apple M3): Carquet 0.7.2 vs. Arrow C++ 24.0.0, gemessen am 2026-10-20. x86 (Xeon D-1531): Carquet 0.4.4 vs. Arrow C++ 23.0.1; dieser Lauf liegt vor der Writer-Parallelisierung von 0.7.1 und wurde seither nicht neu gemessen.

<div align="center">

| | x86 (Xeon D-1531) | | ARM (Apple M3) | |
|---|---|---|---|---|
| **Codec** | **Schreiben** | **Lesen** | **Schreiben** | **Lesen** |
| snappy | **1.55x** | **1.25x** | **4.63x** | **1.66x** |
| zstd | **1.31x** | **1.04x** | **5.66x** | **1.36x** |
| lz4 | **1.02x** | 0.83x | **4.47x** | **1.76x** |
| none | **1.13x** | **40.6x**\* | **1.83x** | **38.9x**\* |

</div>

\* Unkomprimiertes Lesen nutzt mmap Zero-Copy, siehe den Hinweis zur Methodik.

Komprimiertes Lesen umfasst die vollständige Dekompression und Dekodierung jedes einzelnen Werts, ohne Abkürzungen, und beide Bibliotheken verwenden dieselben lz4-/zstd-Shared-Libraries des Systems, die reine Codec-Geschwindigkeit ist also identisch. Der aussagekräftigste Vergleich ist die Tabelle **Cross-Read derselben Datei** (in den vollständigen Ergebnissen), bei der beide Bibliotheken exakt dieselbe Parquet-Datei lesen: In diesem direkten Vergleich liest Carquet komprimierte Daten **1.3–2.1x schneller** als Arrow C++ (x86-Lauf: 1.6–2.6x).

<details>
<summary><b>Benchmark-Methodik</b></summary>

<br/>

Alle Benchmarks verwenden identische Daten (deterministischer LCG-PRNG), identische Parquet-Einstellungen (kein Dictionary, BYTE_STREAM_SPLIT für Gleitkommazahlen, Page-Prüfsummen, mmap-Lesen) und den getrimmten Median aus 11–51 Iterationen, wobei der Page-Cache des Betriebssystems zwischen Schreib- und Lesephase geleert wird und zwischen den Konfigurationen eine Abkühlpause liegt. Schema: 3 Spalten (INT64, DOUBLE, INT32). Verglichen wird mit dem Low-Level-Parquet-Reader von Arrow C++ (unter Umgehung der Arrow-Table-Materialisierung) und mit PyArrow.

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

*Carquet 0.7.2, MacBook Air M3, 16GB RAM, macOS 26.5, Arrow C++ 24.0.0, PyArrow 23.0.1, ZSTD-Level 1. Die Tabellen mit 10 Mio. Zeilen wurden am 2026-10-20 gemessen; die Tabellen mit 1 Mio. und 100K Zeilen stammen aus dem 0.7.1-Lauf vom 2026-10-08 (nicht neu gemessen).*

#### 10 Mio. Zeilen vs. Arrow C++

<div align="center">

| Codec | Carquet Schreiben | Arrow C++ Schreiben | Verh. Schreiben | Carquet Lesen | Arrow C++ Lesen | Verh. Lesen | Größe |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **74.7ms** | 137ms | **1.83x** | **0.36ms** | 14.0ms | **38.9x**\* | 200.0MB |
| snappy | **52.7ms** | 244ms | **4.63x** | **12.4ms** | 20.6ms | **1.66x** | 131.2MB |
| zstd | **60.4ms** | 342ms | **5.66x** | **19.8ms** | 26.9ms | **1.36x** | 99.9MB |
| lz4 | **57.2ms** | 256ms | **4.47x** | **9.22ms** | 16.2ms | **1.76x** | 128.9MB |

</div>

<br/>

#### 1 Mio. Zeilen vs. Arrow C++

<div align="center">

| Codec | Carquet Schreiben | Arrow C++ Schreiben | Verh. Schreiben | Carquet Lesen | Arrow C++ Lesen | Verh. Lesen |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **6.48ms** | 14.82ms | **2.29x** | **0.06ms** | 1.63ms | **27.2x**\* |
| snappy | **14.28ms** | 26.44ms | **1.85x** | **2.14ms** | 2.85ms | **1.33x** |
| zstd | **14.76ms** | 36.17ms | **2.45x** | **2.70ms** | 3.46ms | **1.28x** |
| lz4 | **12.30ms** | 26.10ms | **2.12x** | **1.04ms** | 1.87ms | **1.80x** |

</div>

<br/>

#### 100K Zeilen vs. Arrow C++

<div align="center">

| Codec | Carquet Schreiben | Arrow C++ Schreiben | Verh. Schreiben | Carquet Lesen | Arrow C++ Lesen | Verh. Lesen |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **1.07ms** | 1.79ms | **1.67x** | **0.02ms** | 0.25ms | **12.5x**\* |
| snappy | **1.60ms** | 2.63ms | **1.64x** | **0.34ms** | 0.91ms | **2.68x** |
| zstd | **1.73ms** | 3.63ms | **2.10x** | **0.63ms** | 1.26ms | **2.00x** |
| lz4 | **1.52ms** | 2.62ms | **1.72x** | **0.25ms** | 0.55ms | **2.20x** |

</div>

<br/>

#### Cross-Read derselben Datei (10 Mio. Zeilen)

Beide Bibliotheken lesen **dieselbe** Parquet-Datei, der fairste direkte Vergleich.

<div align="center">

| Codec | Writer | Carquet Lesen | Arrow C++ Lesen | Verhältnis |
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

#### 10 Mio. Zeilen vs. PyArrow

<div align="center">

| Codec | Carquet Schreiben | PyArrow Schreiben | Verh. Schreiben | Carquet Lesen | PyArrow Lesen | Verh. Lesen |
|-------|--------------|---------------|---------|-------------|--------------|---------|
| none | **74.7ms** | 202ms | **2.71x** | **0.36ms** | 36.2ms | **100.5x**\* |
| snappy | **52.7ms** | 312ms | **5.93x** | **12.4ms** | 46.6ms | **3.75x** |
| zstd | **60.4ms** | 410ms | **6.78x** | **19.8ms** | 75.6ms | **3.82x** |
| lz4 | **57.2ms** | 321ms | **5.62x** | **9.22ms** | 40.0ms | **4.33x** |

</div>

<br/>

#### 1 Mio. Zeilen vs. PyArrow

<div align="center">

| Codec | Carquet Schreiben | PyArrow Schreiben | Verh. Schreiben | Carquet Lesen | PyArrow Lesen | Verh. Lesen |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **6.48ms** | 20.17ms | **3.11x** | **0.06ms** | 2.79ms | **46.5x**\* |
| snappy | **14.28ms** | 31.54ms | **2.21x** | **2.14ms** | 4.01ms | **1.87x** |
| zstd | **14.76ms** | 42.51ms | **2.88x** | **2.70ms** | 4.42ms | **1.64x** |
| lz4 | **12.30ms** | 32.03ms | **2.60x** | **1.04ms** | 3.22ms | **3.10x** |

</div>

<br/>

#### 100K Zeilen vs. PyArrow

<div align="center">

| Codec | Carquet Schreiben | PyArrow Schreiben | Verh. Schreiben | Carquet Lesen | PyArrow Lesen | Verh. Lesen |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **1.07ms** | 2.12ms | **1.98x** | **0.02ms** | 0.24ms | **12.0x**\* |
| snappy | **1.60ms** | 3.17ms | **1.98x** | **0.34ms** | 0.60ms | **1.76x** |
| zstd | **1.73ms** | 4.16ms | **2.40x** | **0.63ms** | 0.81ms | **1.29x** |
| lz4 | **1.52ms** | 3.21ms | **2.11x** | **0.25ms** | 0.42ms | **1.68x** |

</div>

\* Zero-Copy-mmap-Pfad

</details>

<details>
<summary><b>Benchmarks selbst ausführen</b></summary>

<br/>

```bash
# Build mit maximalen Optimierungen
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON -DCARQUET_BUILD_DEV=ON
cmake --build build -j$(nproc)

cd build
./benchmark_carquet                     # Carquet allein
python3 ../benchmark/run_benchmark.py   # Vollständiger Vergleich (+ PyArrow, + Arrow C++)

# Die 100M-Zeilen-Konfigurationen (xlarge) überspringen — sie schreiben ~2GB Dateien pro Codec
# und können je nach Hardware 30+ Minuten dauern
python3 ../benchmark/run_benchmark.py --skip-xlarge

# ZSTD-Level überschreiben (Standard: 1)
CARQUET_BENCH_ZSTD_LEVEL=3 python3 ../benchmark/run_benchmark.py
```

Um den Arrow-C++-Vergleich einzubeziehen, aktivieren Sie sein Target beim Konfigurieren:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON \
  -DCARQUET_BUILD_BENCHMARKS=ON \
  -DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON
cmake --build build -j$(nproc)

# Oder auf eine eigene Arrow-Installation zeigen
cmake -B build ... -DCARQUET_ARROW_CPP_ROOT=/path/to/arrow-prefix
```

Der Arrow-C++-Benchmark verwendet die Low-Level-API `parquet::ParquetFileReader` (unter Umgehung des Overheads der Arrow-Table-Materialisierung) mit parallelen Row-Group-Readern. Im Modus **Cross-Read derselben Datei** lesen beide Bibliotheken exakt dieselbe Parquet-Datei, was Unterschiede bei Page-Größen, Encoding und Row-Group-Layout eliminiert. Beide Benchmarks verwenden identische Daten, identische Row-Group-Größen, kein Dictionary, Page-Prüfsummen, mmap-Lesen und BYTE_STREAM_SPLIT für Gleitkommazahlen.

</details>

<br/>

---

<!-- ────────────────────────────  FEATURE SUPPORT  ──────────────────────────── -->

## Unterstützte Parquet-Funktionen

Was die Bibliothek heute implementiert, nach Bereichen.

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

<!-- ────────────────────────────  INTEROPERABILITY  ──────────────────────────── -->

## Interoperabilität

Carquet-Dateien sind vollständig kompatibel mit PyArrow, DuckDB, Spark und jedem anderen Parquet-Reader.

```python
import pyarrow.parquet as pq
table = pq.read_table("carquet_output.parquet")  # funktioniert einfach
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

<!-- ────────────────────────────  CONTRIBUTING & LICENSE  ──────────────────────────── -->

## Mitwirken, Zitieren, Lizenz

- **Mitwirken**: [CONTRIBUTING.md](CONTRIBUTING.md) behandelt Einrichtung, Code-Stil und das Gate pro Pull Request (Tests, Sanitizer, Fuzzing, Interop, Benchmarks). Release Notes stehen in [CHANGELOG.md](CHANGELOG.md).
- **Zitieren**: eine [`CITATION.cff`](CITATION.cff) liegt bei; GitHubs Schaltfläche „Cite this repository" formatiert sie als BibTeX oder APA.
- **Lizenz**: MIT, siehe [LICENSE](LICENSE).
