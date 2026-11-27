<h1 align="center">
  <img src="res/img/carquet_banner.png" alt="Carquet" width="840" />
</h1>

<p align="center">
  Uma biblioteca rápida, em C puro, para ler e escrever arquivos Apache Parquet.
</p>

<p align="center">
  <a href="https://vitruves.github.io/carquet/"><img src="https://img.shields.io/badge/Website-vitruves.github.io%2Fcarquet-2a78d6?style=for-the-badge" alt="Site" /></a>
</p>

<p align="center">
  <a href="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml"><img src="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml/badge.svg" alt="Build" /></a>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue" alt="Plataforma" />
  <img src="https://img.shields.io/badge/C-C11-blue" alt="Padrão C" />
  <a href="https://github.com/Vitruves/carquet/releases"><img src="https://img.shields.io/github/v/release/Vitruves/carquet?color=blue" alt="Release" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="Licença" /></a>
  <a href="https://github.com/Vitruves/carquet/graphs/traffic"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/Vitruves/carquet/traffic/clones.json" alt="Clones" /></a>
  <br/>
  <a href="https://pypi.org/project/carquet/"><img src="https://img.shields.io/pypi/v/carquet?color=blue&label=pypi" alt="PyPI" /></a>
  <a href="https://pypi.org/project/carquet/"><img src="https://img.shields.io/pypi/pyversions/carquet" alt="Versões do Python" /></a>
  <a href="https://pepy.tech/projects/carquet"><img src="https://img.shields.io/pepy/dt/carquet?color=blue&label=pypi%20downloads" alt="Downloads no PyPI" /></a>
  <img src="https://img.shields.io/badge/SIMD-SSE4.2%20%C2%B7%20AVX2%20%C2%B7%20AVX--512%20%C2%B7%20NEON%20%C2%B7%20SVE-eb6834" alt="SIMD" />
</p>

<p align="center">
  <a href="#destaques">Destaques</a> &nbsp;·&nbsp;
  <a href="#benchmarks">Benchmarks</a> &nbsp;·&nbsp;
  <a href="#instalação">Instalação</a> &nbsp;·&nbsp;
  <a href="#configurações-de-build">Configurações de build</a> &nbsp;·&nbsp;
  <a href="#api-c">API C</a> &nbsp;·&nbsp;
  <a href="#api-python">API Python</a> &nbsp;·&nbsp;
  <a href="#ferramenta-de-linha-de-comando">CLI</a> &nbsp;·&nbsp;
  <a href="#suporte-a-recursos-do-parquet">Recursos</a> &nbsp;·&nbsp;
  <a href="docs/README.md">Manual</a> &nbsp;·&nbsp;
  <a href="CHANGELOG.md">Changelog</a>
</p>

<p align="center">
  <a href="README.md">English</a> &nbsp;·&nbsp;
  <a href="README.zh-CN.md">简体中文</a> &nbsp;·&nbsp;
  <a href="README.ja.md">日本語</a> &nbsp;·&nbsp;
  <a href="README.ko.md">한국어</a> &nbsp;·&nbsp;
  <a href="README.de.md">Deutsch</a> &nbsp;·&nbsp;
  <a href="README.fr.md">Français</a> &nbsp;·&nbsp;
  <b>Português (Brasil)</b>
</p>

<br/>

---

<!-- ────────────────────────────  DESTAQUES  ──────────────────────────── -->

## Destaques

<div align="center">

| | |
|---|---|
| ⚡ **Rápida** | Escreve **2.0–5.6x mais rápido** que o Arrow C++ e lê dados comprimidos **1.5–2.1x mais rápido** a partir do mesmíssimo arquivo (Apple M3, 10M de linhas). Leituras não comprimidas são servidas diretamente do mmap em **0.28 ms**. SIMD para SSE4.2 / AVX2 / AVX-512 / NEON / SVE com seleção em tempo de execução. |
| 🪶 **Pequena** | C11 puro, uma biblioteca de **~500 KB** (Arrow C++: 50 MB+). Quatro dependências de compressão (zstd, zlib, lz4, brotli), baixadas automaticamente quando ausentes. Compila em menos de um minuto com um único `make` ou `xmake`. |
| 📦 **Completa** | Todos os tipos físicos e encodings, Snappy / GZIP / LZ4 / ZSTD / Brotli, páginas V1 + V2, schemas aninhados em qualquer profundidade, bloom filters, índices de página, poda baseada em estatísticas, append. [Matriz completa ↓](#suporte-a-recursos-do-parquet) |
| 🔗 **Interoperável** | Arquivos lidos por PyArrow, DuckDB e Spark sem configuração adicional. Ponte sem dependências com a [Arrow C Data Interface](docs/reading.md#export-to-arrow-c-data-interface) (struct / list / map em qualquer profundidade) e uma [API Python](#api-python) completa, sem dependências. |
| 🛠 **Ferramentas** | CLI `carquet` para inspecionar, filtrar, exportar e validar arquivos, mais `codegen` para gerar um leitor C pronto para compilar. Uma [skill de agente](#ai-coding-agents) pronta para Claude Code e Codex. |

</div>

<br/>

---

<!-- ────────────────────────────  BENCHMARKS  ──────────────────────────── -->

## Benchmarks

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="res/img/bench_10m_m3_dark.png">
    <source media="(prefers-color-scheme: light)" srcset="res/img/bench_10m_m3_light.png">
    <img src="res/img/bench_10m_m3_light.png" alt="Tempo de escrita e leitura de 10 milhões de linhas em Apple M3: carquet vs Arrow C++ e bibliotecas Parquet para Python, por codec" width="900">
  </picture>
</p>

<p align="center"><b><a href="https://vitruves.github.io/carquet/">Explore o benchmark interativo no site →</a></b><br/><sub>alterne entre leitura e escrita, e entre 100K, 1M e 10M de linhas</sub></p>

Ganho de velocidade do carquet sobre o Arrow C++ com 10M de linhas (maior = carquet mais rápido):

<div align="center">

| | ARM · Apple M3 | | x86 · Windows (Haswell, MSVC) | | x86 · Linux (Xeon D-1531) | |
|---|---:|---:|---:|---:|---:|---:|
| **Codec** | **Escrita** | **Leitura** | **Escrita** | **Leitura** | **Escrita** | **Leitura** |
| snappy | **4.77x** | **1.63x** | **3.39x** | **1.17x** | **1.55x** | **1.25x** |
| zstd | **5.56x** | **1.40x** | **3.83x** | **1.13x** | **1.31x** | **1.04x** |
| lz4 | **4.62x** | **1.85x** | **2.96x** | **1.05x** | **1.02x** | 0.83x |
| none | **2.04x** | **52.0x**\* | **2.10x** | **62.0x**\* | **1.13x** | **40.6x**\* |

</div>

<sub>\* Caminho mmap zero-copy. Apple M3: carquet 0.8.0 vs Arrow C++ 26.0.1 (2026-11-09). Windows: carquet 0.7.3 vs Arrow C++ 24.0.0 (2026-11-02). Linux: carquet 0.4.4 vs Arrow C++ 23.0.1, antes do paralelismo do writer da 0.7.1, não medido novamente desde então.</sub>

As leituras comprimidas decodificam cada valor, e em Linux e macOS ambas as bibliotecas linkam as mesmas lz4 / zstd do sistema, de modo que a velocidade do codec é idêntica. A comparação mais justa é a **leitura cruzada do mesmo arquivo** (ambas as bibliotecas leem exatamente o mesmo arquivo): o carquet é **1.5–2.1x mais rápido** que o Arrow C++ no M3 (2.0–4.1x mais rápido que o polars), 1.6–2.6x em Linux e 1.04–1.30x na VM Windows, cujas leituras a frio são limitadas pelo disco virtual.

<details>
<summary><b>Metodologia dos benchmarks</b></summary>

<br/>

Todos os benchmarks usam dados idênticos (PRNG LCG determinístico), configurações Parquet idênticas (sem dicionário, BYTE_STREAM_SPLIT para floats, checksums de página, leituras via mmap), mediana aparada de 11–51 iterações, com o cache de páginas do SO limpo entre as fases de escrita e leitura e um intervalo de resfriamento entre as configurações. Schema: 3 colunas (INT64, DOUBLE, INT32). Comparado com o leitor Parquet de baixo nível do Arrow C++ (sem materialização em Arrow Table), com o PyArrow, com o fastparquet e com o leitor e o escritor nativos (Rust) do polars.

O cache de páginas é limpo por arquivo: `posix_fadvise(DONTNEED)` em Linux; em macOS e Windows o arquivo é copiado sem cache (`F_NOCACHE`, `FILE_FLAG_NO_BUFFERING`) e a cópia o substitui (`benchmark/bench_purge_win32.h`). No Windows, o carquet linka as zstd / lz4 que ele mesmo compila e o Arrow C++ as do conda-forge, de modo que a velocidade do codec pode diferir ligeiramente entre os dois.

O benchmark de **leitura cruzada do mesmo arquivo** é a comparação mais justa: ambas as bibliotecas leem exatamente o mesmo arquivo Parquet (escrito por uma, lido por ambas). Isso elimina diferenças de tamanho de página, escolhas de encoding e layout de row groups.

As **leituras não comprimidas** marcadas com \* usam o **caminho mmap zero-copy** do carquet: para colunas obrigatórias, de tamanho fixo, não comprimidas e com encoding PLAIN, o batch reader retorna ponteiros diretamente para o arquivo mapeado em memória, sem memcpy. O Arrow sempre materializa em seus próprios buffers. **Os números de leitura comprimida são a medida mais representativa da vazão de leitura de ponta a ponta.**

</details>

<details>
<summary><b>Resultados completos em x86</b> (Intel Xeon D-1531, Linux)</summary>

<br/>

*12 threads @ 2.2GHz, 32GB de RAM, Ubuntu 24.04, ZSTD nível 1*

#### 10M de linhas vs Arrow C++

<div align="center">

| Codec | Escrita Carquet | Escrita Arrow C++ | Razão E | Leitura Carquet | Leitura Arrow C++ | Razão L | Tamanho |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **1557ms** | 1766ms | **1.13x** | **1.25ms** | 50.8ms | **40.6x**\* | 190.7MB |
| snappy | **1002ms** | 1549ms | **1.55x** | **78ms** | 97.8ms | **1.25x** | 125.1MB |
| zstd | **1311ms** | 1714ms | **1.31x** | **76.8ms** | 80.2ms | **1.04x** | 95.3MB |
| lz4 | **1521ms** | 1554ms | **1.02x** | 59.1ms | **49.0ms** | 0.83x | 122.9MB |

</div>

<br/>

#### 1M de linhas vs Arrow C++

<div align="center">

| Codec | Escrita Carquet | Escrita Arrow C++ | Razão E | Leitura Carquet | Leitura Arrow C++ | Razão L |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **180ms** | 196ms | **1.09x** | **0.22ms** | 6.2ms | **28x**\* |
| snappy | **141ms** | 148ms | **1.05x** | **8.1ms** | 11.6ms | **1.44x** |
| zstd | **131ms** | 185ms | **1.41x** | 10.3ms | **9.1ms** | 0.88x |
| lz4 | **143ms** | 149ms | **1.04x** | 8.5ms | **6.1ms** | 0.72x |

</div>

<br/>

#### 100K linhas vs Arrow C++

<div align="center">

| Codec | Escrita Carquet | Escrita Arrow C++ | Razão E | Leitura Carquet | Leitura Arrow C++ | Razão L |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **14.1ms** | 18.4ms | **1.30x** | **0.11ms** | 2.18ms | **19.8x**\* |
| snappy | **10.1ms** | 10.6ms | **1.05x** | **1.27ms** | 5.97ms | **4.70x** |
| zstd | **8.7ms** | 14.1ms | **1.62x** | **1.58ms** | 3.88ms | **2.46x** |
| lz4 | **9.6ms** | 11.0ms | **1.14x** | **0.77ms** | 2.78ms | **3.61x** |

</div>

<br/>

#### Leitura cruzada do mesmo arquivo (10M de linhas)

Ambas as bibliotecas leem o **mesmo** arquivo Parquet, a comparação mais justa em condições idênticas.

<div align="center">

| Codec | Escritor | Leitura Carquet | Leitura Arrow C++ | Razão |
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

#### 10M de linhas vs PyArrow

<div align="center">

| Codec | Escrita Carquet | Escrita PyArrow | Razão E | Leitura Carquet | Leitura PyArrow | Razão L |
|-------|--------------|---------------|---------|-------------|--------------|---------|
| none | **1557ms** | 1806ms | **1.16x** | **1.25ms** | 213ms | **170x**\* |
| snappy | **1002ms** | 1649ms | **1.65x** | **78ms** | 384ms | **4.91x** |
| zstd | **1311ms** | 1796ms | **1.37x** | **76.8ms** | 369ms | **4.81x** |
| lz4 | **1521ms** | 1676ms | **1.10x** | **59.1ms** | 281ms | **4.76x** |

</div>

\* Caminho mmap zero-copy

</details>

<details>
<summary><b>Resultados completos em ARM</b> (Apple M3, macOS)</summary>

<br/>

*Carquet 0.8.0, MacBook Air M3, 16GB de RAM, macOS 26.5, Arrow C++ 26.0.1, PyArrow 26.0.1, fastparquet 2026.9.0, polars 2.1.0, ZSTD nível 1, medido em 2026-11-09 (configurações de 10M de linhas de uma segunda execução no mesmo dia: a primeira configuração da primeira execução foi perturbada). Relatório bruto: [`benchmark/bench_0.8.0_darwin_arm64_20261109.json`](benchmark/bench_0.8.0_darwin_arm64_20261109.json).*

#### 10M de linhas vs Arrow C++

<div align="center">

| Codec | Escrita Carquet | Escrita Arrow C++ | Razão E | Leitura Carquet | Leitura Arrow C++ | Razão L | Tamanho |
|-------|--------------|---------|---------|-------------|---------|---------|------|
| none | **57.7ms** | 117ms | **2.04x** | **0.28ms** | 14.6ms | **52.0x**\* | 190.7MB |
| snappy | **49.3ms** | 235ms | **4.77x** | **14.0ms** | 22.7ms | **1.63x** | 125.1MB |
| zstd | **61.0ms** | 339ms | **5.56x** | **20.6ms** | 28.8ms | **1.40x** | 95.3MB |
| lz4 | **52.9ms** | 244ms | **4.62x** | **9.83ms** | 18.2ms | **1.85x** | 122.7MB |

</div>

<br/>

#### 1M de linhas vs Arrow C++

<div align="center">

| Codec | Escrita Carquet | Escrita Arrow C++ | Razão E | Leitura Carquet | Leitura Arrow C++ | Razão L |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **5.87ms** | 12.3ms | **2.10x** | **0.06ms** | 1.59ms | **26.5x**\* |
| snappy | **7.40ms** | 25.9ms | **3.49x** | **1.30ms** | 2.53ms | **1.95x** |
| zstd | **8.32ms** | 34.1ms | **4.10x** | **2.36ms** | 3.17ms | **1.34x** |
| lz4 | **7.70ms** | 25.7ms | **3.34x** | **1.06ms** | 2.15ms | **2.03x** |

</div>

<br/>

#### 100K linhas vs Arrow C++

<div align="center">

| Codec | Escrita Carquet | Escrita Arrow C++ | Razão E | Leitura Carquet | Leitura Arrow C++ | Razão L |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **0.87ms** | 1.65ms | **1.90x** | **0.02ms** | 0.23ms | **11.5x**\* |
| snappy | **1.52ms** | 2.62ms | **1.72x** | **0.38ms** | 0.91ms | **2.39x** |
| zstd | **1.90ms** | 3.69ms | **1.94x** | **0.64ms** | 1.32ms | **2.06x** |
| lz4 | **1.66ms** | 2.63ms | **1.58x** | **0.30ms** | 0.63ms | **2.10x** |

</div>

<br/>

#### 10M de linhas vs PyArrow

<div align="center">

| Codec | Escrita Carquet | Escrita PyArrow | Razão E | Leitura Carquet | Leitura PyArrow | Razão L |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **57.7ms** | 165ms | **2.87x** | **0.28ms** | 26.1ms | **93.1x**\* |
| snappy | **49.3ms** | 276ms | **5.60x** | **14.0ms** | 39.3ms | **2.82x** |
| zstd | **61.0ms** | 387ms | **6.34x** | **20.6ms** | 41.3ms | **2.00x** |
| lz4 | **52.9ms** | 287ms | **5.43x** | **9.83ms** | 29.0ms | **2.95x** |

</div>

<br/>

#### 1M de linhas vs PyArrow

<div align="center">

| Codec | Escrita Carquet | Escrita PyArrow | Razão E | Leitura Carquet | Leitura PyArrow | Razão L |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **5.87ms** | 17.5ms | **2.98x** | **0.06ms** | 2.79ms | **46.5x**\* |
| snappy | **7.40ms** | 28.9ms | **3.90x** | **1.30ms** | 3.61ms | **2.78x** |
| zstd | **8.32ms** | 39.4ms | **4.74x** | **2.36ms** | 4.00ms | **1.69x** |
| lz4 | **7.70ms** | 30.5ms | **3.97x** | **1.06ms** | 3.22ms | **3.04x** |

</div>

<br/>

#### 100K linhas vs PyArrow

<div align="center">

| Codec | Escrita Carquet | Escrita PyArrow | Razão E | Leitura Carquet | Leitura PyArrow | Razão L |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **0.87ms** | 2.10ms | **2.41x** | **0.02ms** | 0.26ms | **13.0x**\* |
| snappy | **1.52ms** | 3.00ms | **1.97x** | **0.38ms** | 0.67ms | **1.76x** |
| zstd | **1.90ms** | 4.13ms | **2.17x** | **0.64ms** | 0.73ms | **1.14x** |
| lz4 | **1.66ms** | 3.04ms | **1.83x** | **0.30ms** | 0.44ms | **1.47x** |

</div>

<br/>

#### 10M de linhas vs fastparquet

<div align="center">

| Codec | Escrita Carquet | Escrita fastparquet | Razão E | Leitura Carquet | Leitura fastparquet | Razão L |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **57.7ms** | 64.6ms | **1.12x** | **0.28ms** | 55.7ms | **198.8x**\* |
| snappy | **49.3ms** | 208ms | **4.22x** | **14.0ms** | 103ms | **7.40x** |
| zstd | **61.0ms** | 325ms | **5.33x** | **20.6ms** | 185ms | **8.95x** |
| lz4 | **52.9ms** | 244ms | **4.61x** | **9.83ms** | 55.4ms | **5.64x** |

</div>

<br/>

#### 1M de linhas vs fastparquet

<div align="center">

| Codec | Escrita Carquet | Escrita fastparquet | Razão E | Leitura Carquet | Leitura fastparquet | Razão L |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **5.87ms** | 8.01ms | **1.36x** | **0.06ms** | 3.08ms | **51.3x**\* |
| snappy | **7.40ms** | 22.8ms | **3.08x** | **1.30ms** | 10.2ms | **7.85x** |
| zstd | **8.32ms** | 33.5ms | **4.03x** | **2.36ms** | 18.0ms | **7.61x** |
| lz4 | **7.70ms** | 22.4ms | **2.92x** | **1.06ms** | 5.37ms | **5.07x** |

</div>

<br/>

#### 100K linhas vs fastparquet

<div align="center">

| Codec | Escrita Carquet | Escrita fastparquet | Razão E | Leitura Carquet | Leitura fastparquet | Razão L |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **0.87ms** | 1.08ms | **1.24x** | **0.02ms** | 0.58ms | **29.0x**\* |
| snappy | **1.52ms** | 2.43ms | **1.60x** | **0.38ms** | 1.20ms | **3.16x** |
| zstd | **1.90ms** | 3.47ms | **1.83x** | **0.64ms** | 2.06ms | **3.22x** |
| lz4 | **1.66ms** | 2.31ms | **1.39x** | **0.30ms** | 0.73ms | **2.43x** |

</div>

<br/>

#### 10M de linhas vs polars

<div align="center">

| Codec | Escrita Carquet | Escrita polars | Razão E | Leitura Carquet | Leitura polars | Razão L | Tamanho Carquet | Tamanho polars |
|-------|--------------|---------|---------|-------------|---------|---------|------|------|
| none | 57.7ms | **37.5ms** | 0.65x | **0.28ms** | 28.7ms | **102.4x**\* | 190.7MB | 105.1MB |
| snappy | 49.3ms | **38.1ms** | 0.77x | **14.0ms** | 28.3ms | **2.03x** | 125.1MB | 103.8MB |
| zstd | 61.0ms | **48.6ms** | 0.80x | **20.6ms** | 42.1ms | **2.04x** | 95.3MB | 98.1MB |
| lz4 | 52.9ms | **40.0ms** | 0.76x | **9.83ms** | 29.4ms | **2.99x** | 122.7MB | 104.2MB |

</div>

<br/>

#### 1M de linhas vs polars

<div align="center">

| Codec | Escrita Carquet | Escrita polars | Razão E | Leitura Carquet | Leitura polars | Razão L | Tamanho Carquet | Tamanho polars |
|-------|--------------|---------|---------|-------------|---------|---------|------|------|
| none | 5.87ms | **4.53ms** | 0.77x | **0.06ms** | 3.14ms | **52.3x**\* | 19.1MB | 12.6MB |
| snappy | 7.40ms | **4.82ms** | 0.65x | **1.30ms** | 3.25ms | **2.50x** | 12.5MB | 11.5MB |
| zstd | 8.32ms | **5.83ms** | 0.70x | **2.36ms** | 4.12ms | **1.75x** | 9.6MB | 10.1MB |
| lz4 | 7.70ms | **4.77ms** | 0.62x | **1.06ms** | 3.26ms | **3.08x** | 12.3MB | 11.5MB |

</div>

<br/>

#### 100K linhas vs polars

<div align="center">

| Codec | Escrita Carquet | Escrita polars | Razão E | Leitura Carquet | Leitura polars | Razão L | Tamanho Carquet | Tamanho polars |
|-------|--------------|---------|---------|-------------|---------|---------|------|------|
| none | 0.87ms | **0.72ms** | 0.83x | **0.02ms** | 0.33ms | **16.5x**\* | 1.9MB | 1.3MB |
| snappy | 1.52ms | **0.91ms** | 0.60x | **0.38ms** | 0.46ms | **1.21x** | 1.3MB | 1.1MB |
| zstd | 1.90ms | **1.05ms** | 0.55x | **0.64ms** | 1.08ms | **1.69x** | 1.0MB | 1.0MB |
| lz4 | 1.66ms | **0.96ms** | 0.58x | **0.30ms** | 0.44ms | **1.47x** | 1.2MB | 1.1MB |

</div>

<br/>

#### Leitura cruzada do mesmo arquivo (10M de linhas)

Ambas as bibliotecas leem o **mesmo** arquivo Parquet, a comparação mais justa em condições idênticas.

<div align="center">

| Codec | Escritor | Leitura Carquet | Leitura Arrow C++ | Razão |
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

#### Leitura cruzada do mesmo arquivo vs polars (10M de linhas)

Cada biblioteca lê o arquivo que o benchmark do carquet escreve e o arquivo que o polars escreve com suas próprias configurações (codificação por dicionário em `id` e `category`).

<div align="center">

| Codec | Escritor | Leitura Carquet | Leitura polars | Razão | Arquivo |
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

\* Caminho mmap zero-copy. O fastparquet lê através do pandas e não consegue escrever BYTE_STREAM_SPLIT nem checksums de página, de modo que suas linhas são menos diretamente comparáveis.

O polars usa seu próprio leitor e escritor Parquet em Rust (`use_pyarrow=False`, seu padrão), com os mesmos dados, o mesmo tamanho de row group, o mesmo nível de zstd, as mesmas estatísticas e LZ4_RAW para `lz4`. Ele não consegue escrever BYTE_STREAM_SPLIT nem checksums de página e não tem opção para a codificação por dicionário: ele codifica `id` e `category` por dicionário por conta própria, de modo que, com estes dados, seus arquivos são 45% menores que os do carquet (sem dicionário) sem compressão e 15–17% menores com snappy / lz4 (praticamente do mesmo tamanho com zstd), e suas escritas passam menos bytes pelo codec e pelo disco. Com codificação por dicionário nas mesmas duas colunas, o carquet escreve as 10M de linhas em 45 ms (none), 51 ms (snappy), 51 ms (zstd) e 48 ms (lz4), contra 38, 38, 49 e 40 ms do polars. As leituras são comparadas de forma justa pela tabela de leitura cruzada do mesmo arquivo: o carquet lê o arquivo de qualquer uma das bibliotecas 2.0–4.1x mais rápido com um codec.

</details>

<details>
<summary><b>Resultados completos em Windows</b> (VM Intel Haswell, MSVC)</summary>

<br/>

*Carquet 0.7.3 compilado com MSVC 19.44 (x64 Release, OpenMP), VM OpenStack com 4 vCPU (Intel Haswell, 3.1 GHz), 14.6GB de RAM, disco SAS virtual, Windows Server 2025, Arrow C++ 24.0.0 (conda-forge), ZSTD nível 1, medido em 2026-11-02. As leituras a frio vão ao disco virtual, o que estreita as razões de leitura; as escritas e o caminho zero-copy são limitados pela CPU.*

#### 10M de linhas vs Arrow C++

<div align="center">

| Codec | Escrita Carquet | Escrita Arrow C++ | Razão E | Leitura Carquet | Leitura Arrow C++ | Razão L | Tamanho |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **182ms** | 382ms | **2.10x** | **2.04ms** | 126.5ms | **62.0x**\* | 190.7MB |
| snappy | **175ms** | 594ms | **3.39x** | **115.6ms** | 134.9ms | **1.17x** | 125.1MB |
| zstd | **242ms** | 927ms | **3.83x** | **140.7ms** | 159.2ms | **1.13x** | 95.3MB |
| lz4 | **178ms** | 527ms | **2.96x** | **95.5ms** | 99.9ms | **1.05x** | 122.9MB |

</div>

<br/>

#### 1M de linhas vs Arrow C++

<div align="center">

| Codec | Escrita Carquet | Escrita Arrow C++ | Razão E | Leitura Carquet | Leitura Arrow C++ | Razão L |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **15.6ms** | 37.6ms | **2.41x** | **0.44ms** | 13.72ms | **31.2x**\* |
| snappy | **25.4ms** | 61.2ms | **2.41x** | **14.01ms** | 16.31ms | **1.16x** |
| zstd | **30.4ms** | 94.3ms | **3.10x** | **17.09ms** | 18.38ms | **1.08x** |
| lz4 | **26.4ms** | 53.8ms | **2.04x** | 13.07ms | **12.89ms** | 0.99x |

</div>

<br/>

#### 100K linhas vs Arrow C++

<div align="center">

| Codec | Escrita Carquet | Escrita Arrow C++ | Razão E | Leitura Carquet | Leitura Arrow C++ | Razão L |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **3.53ms** | 3.95ms | **1.12x** | **0.16ms** | 2.13ms | **13.3x**\* |
| snappy | **5.60ms** | 6.23ms | **1.11x** | **2.79ms** | 3.92ms | **1.41x** |
| zstd | **7.11ms** | 9.79ms | **1.38x** | **2.86ms** | 4.36ms | **1.52x** |
| lz4 | **5.75ms** | 5.82ms | **1.01x** | 2.22ms | **2.15ms** | 0.97x |

</div>

<br/>

#### Leitura cruzada do mesmo arquivo (10M de linhas)

<div align="center">

| Codec | Escritor | Leitura Carquet | Leitura Arrow C++ | Razão |
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

\* Caminho mmap zero-copy

</details>

<details>
<summary><b>Executando os benchmarks você mesmo</b></summary>

<br/>

```bash
make benchmark                                    # Release build of the benchmark programs into build/
./build/benchmark_carquet                         # carquet alone
python3 benchmark/run_benchmark.py --skip-xlarge  # full comparison with PyArrow (+ Arrow C++ if built)
```

<div align="center">

| Flag / variável do `run_benchmark.py` | Efeito |
|---|---|
| `--quick` | Apenas as configurações de 10M de linhas (o que o gráfico mostra) |
| `--skip-xlarge` | Pula as configurações de 100M de linhas (~2 GB por codec, mais de 30 minutos) |
| `--fastparquet` | Também compara com o fastparquet (`pip install fastparquet pandas`) |
| `--polars` | Também compara com o leitor / escritor Parquet próprio do polars, em Rust (`pip install polars`) |
| `--no-pyarrow` / `--no-arrow-cpp` | Pula o PyArrow / o Arrow C++ |
| `--iterations N` | Número fixo de iterações (padrão: automático, 11–51) |
| `--build-dir DIR` / `--temp-dir DIR` | Onde ficam os binários / onde os arquivos de teste são escritos |
| `CARQUET_BENCH_ZSTD_LEVEL=3` | Nível do ZSTD (padrão 1) |

</div>

Sempre faça benchmark de um build **Release**: `make benchmark` reconfigura `build/` como Release; um `build/` Debug esquecido faz o carquet parecer várias vezes mais lento. `make native` + `-DCARQUET_BUILD_BENCHMARKS=ON` fornece os números ajustados para a máquina.

O fastparquet usa os mesmos dados, o mesmo tamanho de row group, nenhum dicionário e o mesmo nível de zstd, e LZ4_RAW para `lz4`. Ele não consegue escrever BYTE_STREAM_SPLIT nem checksums de página, e lê através do pandas, de modo que seus números são menos diretamente comparáveis que os do Arrow.

O polars roda com seu leitor e escritor nativos (não através do PyArrow), com os mesmos dados, o mesmo tamanho de row group, o mesmo nível de zstd e LZ4_RAW para `lz4`. Ele não consegue escrever BYTE_STREAM_SPLIT nem checksums de página, e escolhe a codificação por dicionário por conta própria (sem opção), o que torna seus arquivos menores com estes dados.

A comparação com o Arrow C++ precisa de um alvo próprio (Arrow C++ encontrado no sistema, ou indicado com `CARQUET_ARROW_CPP_ROOT`):

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON \
  -DCARQUET_BUILD_BENCHMARKS=ON -DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON \
  [-DCARQUET_ARROW_CPP_ROOT=/path/to/arrow-prefix]
cmake --build build -j
```

Ele usa a API de baixo nível `parquet::ParquetFileReader` (sem materialização em Arrow Table) com leitores de row group paralelos.

No Windows (Developer PowerShell for VS), o Arrow C++ pode vir do conda-forge (`micromamba create -p C:\arrow -c conda-forge libarrow libparquet`); suas DLLs precisam estar no `PATH` quando o benchmark é executado:

```powershell
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCARQUET_BUILD_BENCHMARKS=ON `
  -DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON -DCARQUET_ARROW_CPP_ROOT=C:/arrow/Library
cmake --build build
$env:PATH = "C:\arrow\Library\bin;$env:PATH"
python benchmark\run_benchmark.py --build-dir build --skip-xlarge --no-pyarrow
```

Mantenha as pastas de build e de `--temp-dir` fora da varredura de antivírus em tempo real: o Windows Defender examina cada arquivo que o benchmark escreve.

O gráfico do README é desenhado a partir do JSON que uma execução deixa em `benchmark/`, em uma variante clara e uma escura que o GitHub escolhe conforme o tema do leitor:

```bash
python3 benchmark/plot_readme_chart.py benchmark/bench_<version>_<os>_<arch>_<date>.json   # writes res/img/bench_10m_m3_{light,dark}.png
```

</details>

<br/>

---

<!-- ────────────────────────────  VISÃO RÁPIDA  ──────────────────────────── -->

## Trinta segundos de carquet

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

> **Onde o carquet se encaixa.** O Arrow continua sendo a escolha certa quando você precisa do ecossistema dele (kernels de computação, datasets, Flight, criptografia). O carquet é para programas em C e C++, sistemas embarcados e runtimes de linguagens que só precisam ler e escrever Parquet com rapidez e correção.

<br/>

---

<!-- ────────────────────────────  INSTALAÇÃO  ──────────────────────────── -->

## Instalação

Esta seção instala a **biblioteca C, seus headers e a CLI `carquet`**. Usuários de Python devem usar a [API Python](#api-python) separada (`pip install carquet`).

**Requisitos:** um compilador C11 (GCC 4.9+, Clang 3.4+, MSVC 2015+) e CMake 3.16+ (ou [xmake](#a-partir-do-código-fonte-com-xmake)). zstd, zlib, lz4 e brotli são usados a partir do sistema quando presentes e baixados caso contrário. O OpenMP é opcional (leituras e escritas paralelas).

<br/>

### A partir do código-fonte, com `make`

O caminho mais simples. `make` é um wrapper fino que chama o CMake com as flags certas; não há mais nada a configurar:

```bash
git clone https://github.com/Vitruves/carquet.git && cd carquet
make                  # optimized Release build into build/
sudo make install     # library, headers, CLI, pkg-config and CMake package files
```

Em seguida, linke com `-lcarquet`, `pkg-config --cflags --libs carquet` ou `find_package(carquet)` no CMake.

<br/>

### A partir do código-fonte, com `xmake`

O [xmake](https://xmake.io) não precisa de mais nada instalado: ele baixa versões fixadas de zstd 1.5.6, zlib 1.4, lz4 1.10.0 e brotli 1.1.0 e as linka estaticamente, de modo que a CLI `carquet` não tem nenhuma dependência em tempo de execução.

```bash
git clone https://github.com/Vitruves/carquet.git && cd carquet
xmake                         # Release static library + CLI into build/<plat>/<arch>/release/
xmake install -o /usr/local   # headers, library and CLI (sudo if the prefix needs it)
```

Uma instalação estática contém apenas `libcarquet.a`, então linke os codecs você mesmo (`-lcarquet -lzstd -lz -llz4`, mais `-lbrotlienc -lbrotlidec -lbrotlicommon`, a menos que tenha compilado com `--brotli=n`). Configure com `xmake f --shared=y` para obter uma biblioteca compartilhada autocontida com os codecs embutidos.

<br/>

### No seu projeto CMake

```cmake
include(FetchContent)
FetchContent_Declare(carquet GIT_REPOSITORY https://github.com/Vitruves/carquet.git GIT_TAG main)
FetchContent_MakeAvailable(carquet)
target_link_libraries(myapp PRIVATE carquet)
```

Ou, após `make install`: `find_package(carquet REQUIRED)` e linke `carquet::carquet`.

<br/>

---

<!-- ────────────────────────────  CONFIGURAÇÕES DE BUILD  ──────────────────────────── -->

## Configurações de build

### Alvos do `make`

Cada alvo configura `build/` e o compila em paralelo. `make help` imprime esta lista. Cada um tem um [equivalente no xmake](#xmake).

<div align="center">

| Comando | Tipo de build | O que você obtém |
|---|---|---|
| `make` / `make release` | Release | Biblioteca estática + CLI `carquet` (o padrão) |
| `make native` | Release | O mesmo, ajustado para esta CPU com `-march=native` (o binário não roda em CPUs mais antigas) |
| `make shared` | Release | Biblioteca compartilhada (`.so` / `.dylib` / `.dll`) em vez de estática |
| `make benchmark` | Release | Programas de benchmark (`benchmark_carquet`, alvos de profiling) |
| `make debug` | Debug | Biblioteca + CLI com asserções, sem otimização |
| `make test` | Debug | Suíte de testes, depois a executa com `ctest` |
| `make dev` | Debug | Tudo: testes, exemplos, benchmarks, programas de interop |
| `make fuzz` | Debug | Alvos libFuzzer (use uma toolchain Clang) |
| `make install` | Release | Compila `release` e depois instala em `PREFIX` |
| `make clean` | | Apaga o diretório de build |

</div>

<div align="center">

| Variável | Padrão | Exemplo |
|---|---|---|
| `PREFIX` | Padrão do CMake (`/usr/local` em Unix) | `make install PREFIX=$HOME/.local` |
| `BUILD` | `build` | `make test BUILD=build-debug` (mantém várias configurações lado a lado) |
| `JOBS` | número de núcleos | `make JOBS=4` |
| `CMAKE` | `cmake` | `make CMAKE=/opt/cmake/bin/cmake` |

</div>

Todos os alvos compartilham `build/` por padrão, então `make test` depois de `make benchmark` o reconfigura como Debug. Use `BUILD=` para manter uma árvore Release e uma Debug separadas.

<br/>

### Opções do CMake

Para qualquer coisa que os alvos do `make` não cubram, chame o CMake diretamente e passe opções com `-DOPTION=ON|OFF`:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_BUILD_SHARED=ON -DCARQUET_ENABLE_BROTLI=OFF
cmake --build build -j
sudo cmake --install build --prefix /usr/local
```

<div align="center">

| O que compilar | Padrão | |
|---|:---:|---|
| `CARQUET_BUILD_CLI` | ON | Ferramenta de linha de comando `carquet` |
| `CARQUET_BUILD_SHARED` | OFF | Biblioteca compartilhada em vez de estática |
| `CARQUET_BUILD_DEV` | OFF | Ativa testes, exemplos, benchmarks e interop de uma só vez |
| `CARQUET_BUILD_TESTS` | = DEV | Suíte de testes (`ctest`) |
| `CARQUET_BUILD_EXAMPLES` | = DEV | Programas em [`examples/`](examples/) |
| `CARQUET_BUILD_BENCHMARKS` | = DEV | Programas de benchmark e profiling |
| `CARQUET_BUILD_ARROW_CPP_BENCHMARK` | OFF | Benchmark de comparação com Arrow C++ (requer Arrow C++) |
| `CARQUET_BUILD_INTEROP` | = DEV | Programas de ida e volta C / Python |
| `CARQUET_BUILD_FUZZ` | OFF | Alvos de fuzzing |
| `CARQUET_BUILD_DOCS` | OFF | Alvo `docs`, HTML Doxygen em `build/docs/html/` |

| Dependências | Padrão | |
|---|:---:|---|
| `CARQUET_ENABLE_BROTLI` | ON | Codec BROTLI; OFF remove a dependência da libbrotli |
| `CARQUET_BUNDLE_DEPS` | OFF | Sempre baixa e embute estaticamente zstd / zlib / lz4 / brotli (artefatos portáveis, wheel Python) |
| `CARQUET_ZSTD_TARGET`, `_ZLIB_`, `_LZ4_`, `_BROTLI_TARGET` | vazio | Linka um alvo CMake existente em vez de procurar ou baixar (veja abaixo) |
| `CARQUET_ARROW_CPP_ROOT` | vazio | Prefixo de instalação do Arrow C++ para o benchmark de comparação |
| `CMAKE_DISABLE_FIND_PACKAGE_OpenMP` | OFF | ON compila sem OpenMP (sem `libgomp` / `vcomp140.dll` em tempo de execução) |

| CPU | Padrão | |
|---|:---:|---|
| `CARQUET_NATIVE_ARCH` | OFF | `-march=native`; o mais rápido, binário só para a máquina atual |
| `CARQUET_ENABLE_SSE`, `_AVX`, `_AVX2`, `_AVX512` | ON | Kernels SIMD x86, escolhidos em tempo de execução conforme o que a CPU suporta |
| `CARQUET_ENABLE_NEON` | ON | Kernels ARM NEON |
| `CARQUET_ENABLE_SVE` | OFF | Kernels ARM SVE (experimental) |
| `CARQUET_MSVC_HARDENING` | ON | MSVC: `/guard:cf` e `/CETCOMPAT` |

</div>

Os kernels SIMD são compilados no binário e selecionados em tempo de execução, com fallbacks escalares, de modo que um build padrão roda em qualquer CPU da sua arquitetura. Apenas `CARQUET_NATIVE_ARCH` amarra o binário à máquina de build.

<br/>

### xmake

`xmake.lua` é uma alternativa completa ao CMake, com as mesmas opções e padrões. Configure uma vez com `xmake f` e depois `xmake` compila; a configuração é lembrada até o próximo `xmake f`.

<div align="center">

| Alvo do `make` | Equivalente no xmake |
|---|---|
| `make` | `xmake` |
| `make native` | `xmake f --native_arch=y && xmake` |
| `make shared` | `xmake f --shared=y && xmake` |
| `make benchmark` | `xmake f --benchmarks=y && xmake` |
| `make debug` | `xmake f -m debug && xmake` |
| `make test` | `xmake f -m debug --tests=y && xmake && xmake test` |
| `make dev` | `xmake f -m debug --dev=y && xmake` |
| `make fuzz` | `xmake f -m debug --fuzz=y --toolchain=clang && xmake` |
| `make clean` | `xmake clean` (ou `xmake f -c` para também redefinir a configuração) |
| `cmake --build build --target docs` | `xmake docs` (nenhuma flag de configuração necessária) |

</div>

As opções espelham as do CMake, sem o prefixo `CARQUET_`, em minúsculas, como `--option=y|n`. Combine-as livremente: `xmake f -m release --shared=y --brotli=n --avx512=n`.

<div align="center">

| Opção xmake | Padrão | Descrição |
|--------------|---------|-------------|
| `--dev` | n | Compila testes, exemplos, benchmarks e interop |
| `--tests` / `--examples` / `--benchmarks` / `--interop` | n | Compila um grupo individualmente |
| `--arrow_cpp_benchmark` | n | Benchmark de comparação com Arrow C++ (usa o Arrow do repositório do xmake, se disponível) |
| `--cli` | y | Compila a ferramenta CLI `carquet` |
| `--shared` | n | Compila uma biblioteca compartilhada em vez de estática |
| `--openmp` | y | Leitura paralela de colunas com OpenMP (desativada automaticamente se indisponível) |
| `--native_arch` | n | `-march=native` para desempenho máximo (binário só para a máquina atual) |
| `--sse` / `--avx` / `--avx2` / `--avx512` / `--neon` | y | Conjuntos de instruções SIMD (detectados automaticamente) |
| `--sve` | n | ARM SVE (experimental) |
| `--brotli` | y | Codec BROTLI (libbrotli do repositório do xmake) |
| `--fuzz` | n | Compila os alvos de fuzzing (use `--toolchain=clang`); depois `--fuzz_engine=libFuzzer\|AFL`, `--fuzz_sanitizers=y\|n` |

</div>

<details>
<summary><b>Windows</b> (MSVC, clang-cl, MinGW-w64)</summary>

<br/>

MSVC (x64 e Win32), clang-cl e MinGW-w64 GCC são todos suportados, em build estático ou compartilhado. `CMakePresets.json` tem os presets `msvc`, `msvc-shared`, `clang-cl` e `mingw`:

```powershell
# from a Developer PowerShell for VS (or with MinGW-w64 gcc + ninja on PATH for "mingw")
cmake --preset msvc
cmake --build --preset msvc
ctest --preset msvc
```

- Ao configurar manualmente, defina o compilador C++ além do C (por exemplo `CC=clang-cl CXX=clang-cl`): o projeto zstd baixado habilita C++, e o CMake, caso contrário, escolhe um driver incompatível.
- Com `-DCARQUET_BUILD_SHARED=ON` você obtém `carquet.dll` mais a biblioteca de importação `carquet.lib` (`libcarquet.dll` / `libcarquet.dll.a` com MinGW). Os consumidores precisam definir `CARQUET_BUILD_SHARED`; `find_package(carquet)` e `carquet.pc` fazem isso por você. Coloque a DLL ao lado do seu executável ou no `PATH`.
- Libere a memória que o carquet entrega a você, como a saída de `carquet_writer_get_buffer()`, com `carquet_free()` e não com `free()`. A DLL e o seu programa podem usar heaps de runtime C diferentes (CRT Debug vs Release, `/MT`, MinGW vs MSVC).
- Os caminhos de arquivo podem ser UTF-8 (nomes não ASCII funcionam). A CLI `carquet` aceita argumentos Unicode e imprime UTF-8 no console.
- Builds MinGW com OpenMP precisam de `libgomp-1.dll` e `libwinpthread-1.dll` em tempo de execução. Builds MSVC precisam de `vcomp140.dll` (do Visual C++ Redistributable). Passe `-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=ON` para um build sem essa dependência.
- A heurística de aprendizado de máquina do Windows Defender pode colocar em quarentena `test_encoding_roundtrip.exe` e `test_fuzz_regressions.exe` como `Trojan:Win64/CobaltStrike.RDH!MTB` (um falso positivo causado pelos dados binários de teste embutidos neles; o ctest então os reporta como "Not Run"). Exclua o diretório de build da varredura para executar a suíte completa.

</details>

<details>
<summary><b>Usando seus próprios zstd / zlib / lz4 / brotli</b> (embutindo o carquet em um projeto CMake maior)</summary>

<br/>

Quando o carquet é incluído com `add_subdirectory()` ou `FetchContent`, ele resolve cada biblioteca de compressão nesta ordem:

1. O alvo nomeado por `CARQUET_ZSTD_TARGET`, `CARQUET_ZLIB_TARGET`, `CARQUET_LZ4_TARGET` ou `CARQUET_BROTLI_TARGET`, se definido.
2. Uma cópia do sistema (ignorada quando `CARQUET_BUNDLE_DEPS` está ON).
3. Uma cópia que ele mesmo baixa e compila. Se o seu projeto já compila a biblioteca sob o nome de alvo habitual (`libzstd_static`, `zlibstatic`, `lz4_static`), o carquet reutiliza esse alvo em vez de baixar uma segunda cópia.

Um projeto que vendoriza essas bibliotecas e não tem cópia do sistema, portanto, não precisa de configuração alguma. Para escolher explicitamente, nomeie o alvo:

```cmake
set(CARQUET_ZLIB_TARGET zlibstatic)
add_subdirectory(carquet)
```

Um alvo nomeado é linkado como está, sem busca e sem download, e um nome que não é um alvo gera um erro de configuração em vez de um fallback silencioso. O alvo pode ser importado ou compilado na árvore, estático ou compartilhado, e deve carregar seus diretórios de include como requisito de uso.

</details>

<details>
<summary><b>Testes, sanitizers e documentação da API</b></summary>

<br/>

```bash
make test                                    # Debug build + ctest
make dev && ctest --test-dir build           # everything, including examples and interop programs
cmake -B build -DCARQUET_BUILD_DOCS=ON && cmake --build build --target docs   # Doxygen → build/docs/html/
```

[CONTRIBUTING.md](CONTRIBUTING.md) descreve a verificação exigida por pull request (testes, sanitizers, fuzzing, interop, benchmarks).

</details>

<a id="ai-coding-agents"></a>

<details>
<summary><b>Agentes de IA de programação</b> (uma skill pronta para Claude Code e Codex)</summary>

<br/>

[`.agents/skills/carquet/`](.agents/skills/carquet/) é uma [Agent Skill](https://agentskills.io) que ensina um agente de programação a integrar o carquet ao seu projeto: configuração do build (CMake, FetchContent, pkg-config, xmake), qual API usar para quê, receitas completas, as regras de níveis / null-bitmap / tempo de vida, mapeamentos de tipos e as armadilhas habituais com suas correções. Copie-a ou crie um link para ela no seu repositório e o agente a detecta automaticamente:

```bash
mkdir -p .claude/skills && cp -r path/to/carquet/.agents/skills/carquet .claude/skills/   # Claude Code
mkdir -p .agents/skills && cp -r path/to/carquet/.agents/skills/carquet .agents/skills/   # Codex CLI
ln -s ../../third_party/carquet/.agents/skills/carquet .claude/skills/carquet               # vendored: stays in sync
```

Os mesmos arquivos também funcionam como skill pessoal (`~/.claude/skills/`, `~/.codex/skills/`).

</details>

<br/>

---

<!-- ────────────────────────────  API C  ──────────────────────────── -->

## API C

Um exemplo de escrita e um de leitura abaixo; o [manual em `docs/`](docs/README.md) cobre todo o resto.

<br/>

### Escrever um arquivo Parquet

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

### Ler um arquivo Parquet

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

### Mais receitas

Tudo além da leitura/escrita plana está no manual; cada página aponta para um exemplo executável:

<div align="center">

| Você quer… | Veja |
|---|---|
| Colunas anuláveis, row groups, saída em buffer | [`docs/writing.md`](docs/writing.md) |
| Listas, maps, grupos, níveis de definição/repetição | [`docs/nested-data.md`](docs/nested-data.md) |
| Projeção de colunas, estatísticas, inspeção de metadados | [`docs/reading.md`](docs/reading.md) |
| Predicate pushdown, filtragem no nível de página | [`docs/reading.md`](docs/reading.md) |
| Anexar row groups a um arquivo existente | [`docs/writing.md`](docs/writing.md#append-to-an-existing-file) |
| Compressão, codecs personalizados, ajuste do writer | [`docs/writing.md`](docs/writing.md), [`docs/performance.md`](docs/performance.md) |
| mmap, zero-copy, prebuffering, coalescência de I/O | [`docs/performance.md`](docs/performance.md) |
| Códigos de erro e dicas de recuperação | [`docs/error-handling.md`](docs/error-handling.md) |

</div>

<br/>

### Referência da API

A API completa está em [`include/carquet/carquet.h`](include/carquet/carquet.h). Tipos principais:

<div align="center">

| Tipo | Finalidade |
|------|---------|
| `carquet_reader_t` | Leitor de arquivos (abre a partir de caminho, `FILE*` ou buffer em memória) |
| `carquet_writer_t` | Escritor de arquivos |
| `carquet_batch_reader_t` | Iteração em batches de alto nível |
| `carquet_schema_t` | Definição e introspecção de schema |
| `carquet_error_t` | Informação de erro detalhada (código, mensagem, local no código-fonte, dica de recuperação) |

</div>

O [manual](docs/README.md) explica qual superfície usar em cada situação; o layout do código-fonte e a arquitetura estão documentados em [CONTRIBUTING.md](CONTRIBUTING.md). Uma referência HTML navegável pode ser gerada a partir dos headers com Doxygen, veja *Testes, sanitizers e documentação da API* em [Configurações de build](#configurações-de-build).

<details>
<summary><b>Aplicação de exemplo: MetricStore</b> (<code>mocklib/</code>)</summary>

<br/>

[`mocklib/`](mocklib/) é o **MetricStore**, uma aplicação de exemplo completa e autocontida construída sobre o carquet, útil tanto como referência de uso real quanto como teste de integração de ponta a ponta da API pública.

Ela modela um armazenamento de telemetria em séries temporais (ingere eventos, depois os inspeciona e consulta) e, ao fazê-lo, exercita **117 das 143 funções públicas do carquet (~82 %)**:

- **Lado da escrita**: construção de schema com tipos lógicos, ajuste de encoding / compressão / bloom por coluna, índices de página e estatísticas
- **Lado da leitura**: projeção de colunas, predicate pushdown, filtragem no nível de página, teste de pertinência em bloom filter, a ponte com a Arrow C Data Interface, reconstrução de `LIST` aninhadas, I/O em buffer/append e introspecção de metadados

Ela linka o carquet como qualquer projeto downstream faria (`find_package(carquet)` / `add_subdirectory`), traz um teste de ida e volta autoverificável, e seu binário `write_sample` emite um arquivo para inspeção externa. Veja [`mocklib/README.md`](mocklib/README.md).

</details>

<br/>

---

<!-- ────────────────────────────  API PYTHON  ──────────────────────────── -->

## API Python

O `carquet` no PyPI não é um wrapper das chamadas C acima: é uma **API pythônica, no estilo do `pyarrow.parquet`** (`read_table`, `ParquetFile`, `ParquetWriter`, `filters` em DNF) construída sobre o mesmo motor. O wheel embute a biblioteca C com todos os codecs compilados; os dados chegam ao Python pela Arrow C Data Interface, de modo que decodificação, descompressão e I/O rodam em C sem o GIL. **Não tem dependências**: as leituras devolvem uma `carquet.Table` que se converte sozinha em valores Python, NumPy e pandas, e que polars, DuckDB e PyArrow recebem sem cópia.

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

| Você quer… | Use |
|---|---|
| Ler um arquivo, um diretório, bytes ou um objeto de arquivo | `cq.read_table`, `cq.read_pandas` |
| Processar arquivos grandes em streaming | `ParquetFile.iter_batches`, `cq.iter_batches` |
| Podar com estatísticas, índices de página, bloom filters | `filters=` (DNF, resultados exatos), `ParquetFile.bloom_filter` |
| Escrever tabelas, DataFrames, polars, resultados do DuckDB | `cq.write_table`, `cq.ParquetWriter` (incremental, `append=True`) |
| Passar os dados para pandas, NumPy, polars, DuckDB, PyArrow | `t.to_pandas()`, `t.column("x").to_numpy()`, `polars.DataFrame(t)`, `pyarrow.table(t)` |
| Inspecionar um arquivo sem ler seus dados | `cq.read_metadata`, `cq.read_parquet_schema`, `cq.info`, `cq.validate` |

</div>

Os tipos fazem a ida e volta exatamente como com o `pyarrow.parquet` (fusos horários, categóricos, decimais, tipos large, dados aninhados). A referência completa está em [`python/README.md`](python/README.md).

<br/>

---

<!-- ────────────────────────────  FERRAMENTA CLI  ──────────────────────────── -->

## Ferramenta de linha de comando

O carquet inclui uma ferramenta de linha de comando para inspecionar arquivos Parquet e gerar código C de leitura. Ela é compilada e instalada por padrão junto com a biblioteca.

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

Uso típico:

```bash
carquet schema data.parquet
carquet head -n 20 data.parquet
carquet stat data.parquet
carquet validate data.parquet
```

<br/>

### Filtragem de linhas

`cat`, `count`, `head`, `tail`, `sample` e `export` aceitam `-p / --filter EXPR` para empurrar um predicado de linha até o nível de página: só são descomprimidas as páginas cujo min/max do índice de colunas pode satisfazer o predicado.

```bash
carquet cat -p "price > 100 AND status = 'active'" data.parquet
carquet count --filter "id >= 1000" data.parquet
carquet export --filter "ts IS NOT NULL" -c id,ts data.parquet
```

A gramática é `column OP value [AND column OP value]...` com `OP` ∈ {`=`, `==`, `!=`, `<>`, `<`, `<=`, `>`, `>=`}, mais `column IS NULL` / `column IS NOT NULL`.

> A filtragem exige que o arquivo tenha um índice de página (`write_page_index = true`).

<br/>

### Geração de código

Gere um leitor C completo e compilável a partir do schema de qualquer arquivo Parquet:

```bash
carquet codegen -f data.parquet -o reader.c
# Generated: reader.c
# Compile:   clang -o reader reader.c -I.../include -L.../build -lcarquet ...

./reader                    # reads data.parquet (embedded as default)
./reader other.parquet      # override with a different file
```

<details>
<summary><b>Opções do codegen</b></summary>

<br/>

<div align="center">

| Flag | Descrição |
|------|-------------|
| `-f`, `--file FILE` | Arquivo Parquet do qual inspecionar o schema |
| `-o`, `--output FILE` | Arquivo-fonte de saída (padrão: stdout) |
| `--mmap` | Usa I/O mapeado em memória no código gerado |
| `--skeleton` | Gera um `process_batch` vazio para lógica personalizada |
| `-c`, `--columns COLS` | Filtro de colunas separadas por vírgula |
| `-b`, `--batch-size N` | Tamanho do batch (padrão: 1024) |

</div>

</details>

<br/>

---

<!-- ────────────────────────────  RECURSOS SUPORTADOS  ──────────────────────────── -->

## Suporte a recursos do Parquet

O que a biblioteca implementa hoje, por área.

<br/>

### Formato

<div align="center">

| Recurso | Status |
|---------|--------|
| Tipos físicos | ✓  Todos os 8 (de BOOLEAN a FIXED_LEN_BYTE_ARRAY) |
| Tipos lógicos | ✓  STRING, DATE, TIME, TIMESTAMP, DECIMAL, UUID, JSON, INTERVAL, FLOAT16, VARIANT, GEOMETRY, GEOGRAPHY |
| Encodings | ✓  PLAIN, RLE, DICTIONARY, DELTA_BINARY_PACKED, DELTA_LENGTH_BYTE_ARRAY, DELTA_BYTE_ARRAY, BYTE_STREAM_SPLIT (leitura + escrita) |
| Versões de Data Page | ✓  V1 (padrão) e V2 (leitura + escrita) |
| Compressão | ✓  UNCOMPRESSED, SNAPPY, GZIP, LZ4 (com framing Hadoop, codec 5), LZ4_RAW (codec 7), ZSTD, BROTLI (`CARQUET_ENABLE_BROTLI`, padrão ON); LZO via `carquet_register_codec()` |
| Schemas aninhados | ✓  Grupos, listas, maps com níveis de definição/repetição; auto-shredding de LIST/MAP de um nível na escrita (`carquet_writer_write_list_column`) e reconstrução de List na leitura; emissão de `ARROW:schema` aninhado |
| Criptografia | ✗ Não suportada; um arquivo com footer criptografado (`PARE`) é reportado como `CARQUET_ERROR_NOT_IMPLEMENTED` com uma mensagem explícita |

</div>

<details>
<summary><b>Indexação &amp; poda</b></summary>

<br/>

<div align="center">

| Recurso | Status |
|---------|--------|
| Bloom filters | ✓  Leitura, escrita e consulta (`carquet_bloom_filter_check_*`); consultados automaticamente para poda de row groups |
| Índices de página | ✓  Índice de colunas + índice de offsets (leitura + escrita + acesso a estatísticas por página) |
| Estatísticas | ✓  Min/max/contagem de nulos por column chunk; `distinct_count` exato para colunas com dicionário; SizeStatistics do Parquet 2.9 (bytes não codificados de byte-array + histogramas de níveis) |
| Predicate pushdown | ✓  Poda automática de row groups via estatísticas + bloom filters (sem necessidade de callback); filtragem no nível de página via índice de colunas (`carquet_batch_reader_set_page_filter`) |
| Estatísticas geoespaciais | ✓  Bounding box e tipos de geometria por coluna para GEOMETRY/GEOGRAPHY (`carquet_reader_geospatial_statistics`) |
| Projeção de colunas | ✓  Lê apenas as colunas selecionadas |
| Janela de linhas | ✓  Lê uma fatia `[offset, limit)` de linhas, pulando páginas inteiras (`carquet_batch_reader_set_row_range`) |

</div>

</details>

<details>
<summary><b>I/O &amp; extensibilidade</b></summary>

<br/>

<div align="center">

| Recurso | Status |
|---------|--------|
| Append | ✓  Adiciona row groups a um arquivo existente (`carquet_writer_open_append`) |
| Codecs personalizados | ✓  Registra uma implementação própria de compressão/descompressão por slot de codec (`carquet_register_codec`) |
| Metadados chave-valor | ✓  Lê e escreve metadados arbitrários no footer |
| Metadados por campo | ✓  `Field.custom_metadata` do Arrow (rótulos/descrições de variáveis) via `ARROW:schema` (leitura + escrita) |
| Opções por coluna | ✓  Encoding, compressão, estatísticas e bloom filter por coluna |
| Escrita em buffer | ✓  Escreve Parquet em um buffer na memória |
| CRC32 | ✓  Verificação no nível de página (acelerada por hardware em ARM) |
| I/O mapeado em memória | ✓  Leituras zero-copy para dados PLAIN não comprimidos |
| Leituras preservando dicionário | ✓  Retorna índices RLE + dicionário em vez de valores materializados (`preserve_dictionaries`, `carquet_row_batch_column_dictionary`) |
| Coalescência de I/O | ✓  Pré-carrega leituras de várias colunas em um único I/O |
| Footer especulativo | ✓  Abertura de arquivo com um único I/O para a maioria dos arquivos |
| Leituras paralelas com OpenMP | ✓  Quando disponível |
| Escritas paralelas | ✓  Encode/compressão paralelos no nível de página entre colunas (OpenMP) + I/O de row groups em segundo plano (`async_io`) |

</div>

</details>

<br/>

---

<!-- ────────────────────────────  INTEROPERABILIDADE  ──────────────────────────── -->

## Interoperabilidade

Os arquivos do carquet são totalmente compatíveis com PyArrow, DuckDB, Spark e qualquer leitor Parquet.

```python
import pyarrow.parquet as pq
table = pq.read_table("carquet_output.parquet")  # just works
```

```sql
-- DuckDB
SELECT * FROM read_parquet('carquet_output.parquet');
```

<details>
<summary><b>Suíte de testes de interop</b></summary>

<br/>

Idas e voltas bidirecionais (C escreve / Python lê, e vice-versa):

```bash
cmake -B build -DCARQUET_BUILD_INTEROP=ON && cmake --build build
python3 interop/run_interop.py
```

Arquivos Snappy grandes com vários row groups (muitos row groups, colunas `BYTE_ARRAY` largas, offsets de página além de 2 GiB) têm seu próprio gerador, que constrói o arquivo com PyArrow e verifica que o carquet decodifica cada página:

```bash
CARQUET_BIN=build/carquet python3 interop/snappy_multi_rowgroup.py --rows 150000 --row-groups 92
```

</details>

<br/>

---

<!-- ────────────────────────────  CONTRIBUIÇÃO E LICENÇA  ──────────────────────────── -->

## Contribuição, citação, licença

- **Contribuição**: [CONTRIBUTING.md](CONTRIBUTING.md) cobre a configuração, o estilo de código e a verificação exigida por pull request (testes, sanitizers, fuzzing, interop, benchmarks). As notas de release estão em [CHANGELOG.md](CHANGELOG.md).
- **Citação**: um [`CITATION.cff`](CITATION.cff) está incluído; o botão "Cite this repository" do GitHub o formata como BibTeX ou APA.
- **Licença**: MIT, veja [LICENSE](LICENSE).
