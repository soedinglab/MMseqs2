<h1 align="center">
  <img src="res/img/carquet_banner.png" alt="Carquet" width="840" />
</h1>

<p align="center">
  Apache Parquet ファイルを読み書きするための、高速なピュア C ライブラリです。
</p>

<p align="center">
  <a href="https://vitruves.github.io/carquet/"><img src="https://img.shields.io/badge/Website-vitruves.github.io%2Fcarquet-2a78d6?style=for-the-badge" alt="ウェブサイト" /></a>
</p>

<p align="center">
  <a href="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml"><img src="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml/badge.svg" alt="ビルド" /></a>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue" alt="プラットフォーム" />
  <img src="https://img.shields.io/badge/C-C11-blue" alt="C 規格" />
  <a href="https://github.com/Vitruves/carquet/releases"><img src="https://img.shields.io/github/v/release/Vitruves/carquet?color=blue" alt="リリース" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="ライセンス" /></a>
  <a href="https://github.com/Vitruves/carquet/graphs/traffic"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/Vitruves/carquet/traffic/clones.json" alt="クローン数" /></a>
  <br/>
  <a href="https://pypi.org/project/carquet/"><img src="https://img.shields.io/pypi/v/carquet?color=blue&label=pypi" alt="PyPI" /></a>
  <a href="https://pypi.org/project/carquet/"><img src="https://img.shields.io/pypi/pyversions/carquet" alt="Python バージョン" /></a>
  <a href="https://pepy.tech/projects/carquet"><img src="https://img.shields.io/pepy/dt/carquet?color=blue&label=pypi%20downloads" alt="PyPI ダウンロード数" /></a>
  <img src="https://img.shields.io/badge/SIMD-SSE4.2%20%C2%B7%20AVX2%20%C2%B7%20AVX--512%20%C2%B7%20NEON%20%C2%B7%20SVE-eb6834" alt="SIMD" />
</p>

<p align="center">
  <a href="#ハイライト">ハイライト</a> &nbsp;·&nbsp;
  <a href="#ベンチマーク">ベンチマーク</a> &nbsp;·&nbsp;
  <a href="#インストール">インストール</a> &nbsp;·&nbsp;
  <a href="#ビルド構成">ビルド構成</a> &nbsp;·&nbsp;
  <a href="#c-api">C API</a> &nbsp;·&nbsp;
  <a href="#python-api">Python API</a> &nbsp;·&nbsp;
  <a href="#cli-ツール">CLI</a> &nbsp;·&nbsp;
  <a href="#parquet-機能サポート">機能</a> &nbsp;·&nbsp;
  <a href="docs/README.md">マニュアル</a> &nbsp;·&nbsp;
  <a href="CHANGELOG.md">変更履歴</a>
</p>

<p align="center">
  <a href="README.md">English</a> &nbsp;·&nbsp;
  <a href="README.zh-CN.md">简体中文</a> &nbsp;·&nbsp;
  <b>日本語</b> &nbsp;·&nbsp;
  <a href="README.ko.md">한국어</a> &nbsp;·&nbsp;
  <a href="README.de.md">Deutsch</a> &nbsp;·&nbsp;
  <a href="README.fr.md">Français</a> &nbsp;·&nbsp;
  <a href="README.pt-BR.md">Português (Brasil)</a>
</p>

<br/>

---

<!-- ────────────────────────────  ハイライト  ──────────────────────────── -->

## ハイライト

<div align="center">

| | |
|---|---|
| ⚡ **高速** | 書き込みは Arrow C++ より **2.0–5.6x** 速く、まったく同じファイルの圧縮データ読み込みは **1.5–2.1x** 高速です（Apple M3、1000 万行）。非圧縮データの読み込みは mmap から直接 **0.28 ms** で提供されます。SSE4.2 / AVX2 / AVX-512 / NEON / SVE 向け SIMD を実行時ディスパッチで利用します。 |
| 🪶 **小さい** | ピュア C11、ライブラリは **約 500 KB**（Arrow C++：50 MB 以上）。圧縮の依存関係は 4 つ（zstd、zlib、lz4、brotli）で、見つからない場合は自動取得されます。`make` または `xmake` 一発で 1 分以内にビルドできます。 |
| 📦 **完全** | すべての物理型とエンコーディング、Snappy / GZIP / LZ4 / ZSTD / Brotli、V1 + V2 ページ、任意の深さのネストされたスキーマ、ブルームフィルタ、ページインデックス、統計に基づくプルーニング、追記。[全対応表 ↓](#parquet-機能サポート) |
| 🔗 **相互運用** | 生成ファイルは PyArrow、DuckDB、Spark でそのまま読めます。依存ゼロの [Arrow C Data Interface](docs/reading.md#export-to-arrow-c-data-interface) ブリッジ（任意の深さの struct / list / map）と、依存パッケージのない完全な [Python API](#python-api) を備えます。 |
| 🛠 **ツール** | ファイルの検査、フィルタリング、エクスポート、検証を行う `carquet` CLI と、すぐにコンパイルできる C リーダーを出力する `codegen`。Claude Code と Codex 向けのすぐに使える [エージェントスキル](#ai-coding-agents) もあります。 |

</div>

<br/>

---

<!-- ────────────────────────────  ベンチマーク  ──────────────────────────── -->

## ベンチマーク

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="res/img/bench_10m_m3_dark.png">
    <source media="(prefers-color-scheme: light)" srcset="res/img/bench_10m_m3_light.png">
    <img src="res/img/bench_10m_m3_light.png" alt="Apple M3 における 1000 万行の書き込み・読み込み時間：carquet vs Arrow C++ および Python の Parquet ライブラリ、コーデック別" width="900">
  </picture>
</p>

<p align="center"><b><a href="https://vitruves.github.io/carquet/">ウェブサイトでインタラクティブなベンチマークを見る →</a></b><br/><sub>読み込みと書き込み、10 万行・100 万行・1000 万行を切り替えられます</sub></p>

1000 万行での Arrow C++ に対する carquet の高速化率（高いほど carquet が高速）：

<div align="center">

| | ARM · Apple M3 | | x86 · Windows (Haswell, MSVC) | | x86 · Linux (Xeon D-1531) | |
|---|---:|---:|---:|---:|---:|---:|
| **コーデック** | **書き込み** | **読み込み** | **書き込み** | **読み込み** | **書き込み** | **読み込み** |
| snappy | **4.77x** | **1.63x** | **3.39x** | **1.17x** | **1.55x** | **1.25x** |
| zstd | **5.56x** | **1.40x** | **3.83x** | **1.13x** | **1.31x** | **1.04x** |
| lz4 | **4.62x** | **1.85x** | **2.96x** | **1.05x** | **1.02x** | 0.83x |
| none | **2.04x** | **52.0x**\* | **2.10x** | **62.0x**\* | **1.13x** | **40.6x**\* |

</div>

<sub>\* mmap ゼロコピーパス。Apple M3：carquet 0.8.0 vs Arrow C++ 26.0.1（2026-11-09）。Windows：carquet 0.7.3 vs Arrow C++ 24.0.0（2026-11-02）。Linux：carquet 0.4.4 vs Arrow C++ 23.0.1。0.7.1 のライター並列化より前の計測で、その後再計測されていません。</sub>

圧縮読み込みはすべての値をデコードします。また Linux と macOS では両ライブラリが同じシステムの lz4 / zstd をリンクするため、コーデック速度は同一です。最も公平な比較は**同一ファイルのクロスリード**（両ライブラリがまったく同じファイルを読む）で、carquet は M3 で Arrow C++ より **1.5–2.1x**（polars より 2.0–4.1x）、Linux で 1.6–2.6x、Windows VM で 1.04–1.30x 高速です。Windows VM のコールド読み込みは仮想ディスクに律速されます。

<details>
<summary><b>ベンチマークの方法論</b></summary>

<br/>

すべてのベンチマークは、同一のデータ（決定論的 LCG PRNG）、同一の Parquet 設定（辞書なし、浮動小数点に BYTE_STREAM_SPLIT、ページチェックサム、mmap 読み込み）を使用し、11〜51 回の反復のトリム中央値を取り、書き込みフェーズと読み込みフェーズの間に OS ページキャッシュをパージし、構成間にクールダウンを挟んでいます。スキーマ：3 カラム（INT64、DOUBLE、INT32）。比較対象は Arrow C++ の低レベル Parquet リーダー（Arrow Table のマテリアライズを回避）、PyArrow、fastparquet、polars のネイティブ（Rust）リーダーとライターです。

ページキャッシュはファイルごとにパージされます。Linux では `posix_fadvise(DONTNEED)`、macOS と Windows ではファイルをキャッシュなしでコピーし（`F_NOCACHE`、`FILE_FLAG_NO_BUFFERING`）、そのコピーで置き換えます（`benchmark/bench_purge_win32.h`）。Windows では carquet は自身でビルドした zstd / lz4 を、Arrow C++ は conda-forge のものをリンクするため、両者のコーデック速度はわずかに異なる場合があります。

**同一ファイルのクロスリード**ベンチマークが最も公平な比較です。両ライブラリがまったく同じ Parquet ファイル（一方が書き、両方が読む）を読みます。これにより、ページサイズ、エンコーディングの選択、行グループのレイアウトの違いが排除されます。

\* 付きの**非圧縮読み込み**は、Carquet の **mmap ゼロコピーパス**を使用します。PLAIN エンコード、非圧縮、固定サイズ、required なカラムに対して、バッチリーダーは memcpy なしでメモリマップされたファイルへのポインタを直接返します。Arrow は常に自身のバッファへマテリアライズします。**圧縮読み込みの数値が、エンドツーエンドの読み込みスループットの最も代表的な指標です。**

</details>

<details>
<summary><b>x86 の完全な結果</b>（Intel Xeon D-1531、Linux）</summary>

<br/>

*12 スレッド @ 2.2GHz、32GB RAM、Ubuntu 24.04、ZSTD レベル 1*

#### 1000 万行 vs Arrow C++

<div align="center">

| コーデック | Carquet 書き込み | Arrow C++ 書き込み | 書き込み比 | Carquet 読み込み | Arrow C++ 読み込み | 読み込み比 | サイズ |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **1557ms** | 1766ms | **1.13x** | **1.25ms** | 50.8ms | **40.6x**\* | 190.7MB |
| snappy | **1002ms** | 1549ms | **1.55x** | **78ms** | 97.8ms | **1.25x** | 125.1MB |
| zstd | **1311ms** | 1714ms | **1.31x** | **76.8ms** | 80.2ms | **1.04x** | 95.3MB |
| lz4 | **1521ms** | 1554ms | **1.02x** | 59.1ms | **49.0ms** | 0.83x | 122.9MB |

</div>

<br/>

#### 100 万行 vs Arrow C++

<div align="center">

| コーデック | Carquet 書き込み | Arrow C++ 書き込み | 書き込み比 | Carquet 読み込み | Arrow C++ 読み込み | 読み込み比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **180ms** | 196ms | **1.09x** | **0.22ms** | 6.2ms | **28x**\* |
| snappy | **141ms** | 148ms | **1.05x** | **8.1ms** | 11.6ms | **1.44x** |
| zstd | **131ms** | 185ms | **1.41x** | 10.3ms | **9.1ms** | 0.88x |
| lz4 | **143ms** | 149ms | **1.04x** | 8.5ms | **6.1ms** | 0.72x |

</div>

<br/>

#### 10 万行 vs Arrow C++

<div align="center">

| コーデック | Carquet 書き込み | Arrow C++ 書き込み | 書き込み比 | Carquet 読み込み | Arrow C++ 読み込み | 読み込み比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **14.1ms** | 18.4ms | **1.30x** | **0.11ms** | 2.18ms | **19.8x**\* |
| snappy | **10.1ms** | 10.6ms | **1.05x** | **1.27ms** | 5.97ms | **4.70x** |
| zstd | **8.7ms** | 14.1ms | **1.62x** | **1.58ms** | 3.88ms | **2.46x** |
| lz4 | **9.6ms** | 11.0ms | **1.14x** | **0.77ms** | 2.78ms | **3.61x** |

</div>

<br/>

#### 同一ファイルのクロスリード（1000 万行）

両ライブラリが**同じ** Parquet ファイルを読みます。最も公平な同条件比較です。

<div align="center">

| コーデック | ライター | Carquet 読み込み | Arrow C++ 読み込み | 比率 |
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

#### 1000 万行 vs PyArrow

<div align="center">

| コーデック | Carquet 書き込み | PyArrow 書き込み | 書き込み比 | Carquet 読み込み | PyArrow 読み込み | 読み込み比 |
|-------|--------------|---------------|---------|-------------|--------------|---------|
| none | **1557ms** | 1806ms | **1.16x** | **1.25ms** | 213ms | **170x**\* |
| snappy | **1002ms** | 1649ms | **1.65x** | **78ms** | 384ms | **4.91x** |
| zstd | **1311ms** | 1796ms | **1.37x** | **76.8ms** | 369ms | **4.81x** |
| lz4 | **1521ms** | 1676ms | **1.10x** | **59.1ms** | 281ms | **4.76x** |

</div>

\* ゼロコピー mmap パス

</details>

<details>
<summary><b>ARM の完全な結果</b>（Apple M3、macOS）</summary>

<br/>

*Carquet 0.8.0、MacBook Air M3、16GB RAM、macOS 26.5、Arrow C++ 26.0.1、PyArrow 26.0.1、fastparquet 2026.9.0、polars 2.1.0、ZSTD レベル 1、2026-11-09 計測（1000 万行の構成は同日の 2 回目の実行から。1 回目の実行の最初の構成が乱れたため）。生のレポート：[`benchmark/bench_0.8.0_darwin_arm64_20261109.json`](benchmark/bench_0.8.0_darwin_arm64_20261109.json)。*

#### 1000 万行 vs Arrow C++

<div align="center">

| コーデック | Carquet 書き込み | Arrow C++ 書き込み | 書き込み比 | Carquet 読み込み | Arrow C++ 読み込み | 読み込み比 | サイズ |
|-------|--------------|---------|---------|-------------|---------|---------|------|
| none | **57.7ms** | 117ms | **2.04x** | **0.28ms** | 14.6ms | **52.0x**\* | 190.7MB |
| snappy | **49.3ms** | 235ms | **4.77x** | **14.0ms** | 22.7ms | **1.63x** | 125.1MB |
| zstd | **61.0ms** | 339ms | **5.56x** | **20.6ms** | 28.8ms | **1.40x** | 95.3MB |
| lz4 | **52.9ms** | 244ms | **4.62x** | **9.83ms** | 18.2ms | **1.85x** | 122.7MB |

</div>

<br/>

#### 100 万行 vs Arrow C++

<div align="center">

| コーデック | Carquet 書き込み | Arrow C++ 書き込み | 書き込み比 | Carquet 読み込み | Arrow C++ 読み込み | 読み込み比 |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **5.87ms** | 12.3ms | **2.10x** | **0.06ms** | 1.59ms | **26.5x**\* |
| snappy | **7.40ms** | 25.9ms | **3.49x** | **1.30ms** | 2.53ms | **1.95x** |
| zstd | **8.32ms** | 34.1ms | **4.10x** | **2.36ms** | 3.17ms | **1.34x** |
| lz4 | **7.70ms** | 25.7ms | **3.34x** | **1.06ms** | 2.15ms | **2.03x** |

</div>

<br/>

#### 10 万行 vs Arrow C++

<div align="center">

| コーデック | Carquet 書き込み | Arrow C++ 書き込み | 書き込み比 | Carquet 読み込み | Arrow C++ 読み込み | 読み込み比 |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **0.87ms** | 1.65ms | **1.90x** | **0.02ms** | 0.23ms | **11.5x**\* |
| snappy | **1.52ms** | 2.62ms | **1.72x** | **0.38ms** | 0.91ms | **2.39x** |
| zstd | **1.90ms** | 3.69ms | **1.94x** | **0.64ms** | 1.32ms | **2.06x** |
| lz4 | **1.66ms** | 2.63ms | **1.58x** | **0.30ms** | 0.63ms | **2.10x** |

</div>

<br/>

#### 1000 万行 vs PyArrow

<div align="center">

| コーデック | Carquet 書き込み | PyArrow 書き込み | 書き込み比 | Carquet 読み込み | PyArrow 読み込み | 読み込み比 |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **57.7ms** | 165ms | **2.87x** | **0.28ms** | 26.1ms | **93.1x**\* |
| snappy | **49.3ms** | 276ms | **5.60x** | **14.0ms** | 39.3ms | **2.82x** |
| zstd | **61.0ms** | 387ms | **6.34x** | **20.6ms** | 41.3ms | **2.00x** |
| lz4 | **52.9ms** | 287ms | **5.43x** | **9.83ms** | 29.0ms | **2.95x** |

</div>

<br/>

#### 100 万行 vs PyArrow

<div align="center">

| コーデック | Carquet 書き込み | PyArrow 書き込み | 書き込み比 | Carquet 読み込み | PyArrow 読み込み | 読み込み比 |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **5.87ms** | 17.5ms | **2.98x** | **0.06ms** | 2.79ms | **46.5x**\* |
| snappy | **7.40ms** | 28.9ms | **3.90x** | **1.30ms** | 3.61ms | **2.78x** |
| zstd | **8.32ms** | 39.4ms | **4.74x** | **2.36ms** | 4.00ms | **1.69x** |
| lz4 | **7.70ms** | 30.5ms | **3.97x** | **1.06ms** | 3.22ms | **3.04x** |

</div>

<br/>

#### 10 万行 vs PyArrow

<div align="center">

| コーデック | Carquet 書き込み | PyArrow 書き込み | 書き込み比 | Carquet 読み込み | PyArrow 読み込み | 読み込み比 |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **0.87ms** | 2.10ms | **2.41x** | **0.02ms** | 0.26ms | **13.0x**\* |
| snappy | **1.52ms** | 3.00ms | **1.97x** | **0.38ms** | 0.67ms | **1.76x** |
| zstd | **1.90ms** | 4.13ms | **2.17x** | **0.64ms** | 0.73ms | **1.14x** |
| lz4 | **1.66ms** | 3.04ms | **1.83x** | **0.30ms** | 0.44ms | **1.47x** |

</div>

<br/>

#### 1000 万行 vs fastparquet

<div align="center">

| コーデック | Carquet 書き込み | fastparquet 書き込み | 書き込み比 | Carquet 読み込み | fastparquet 読み込み | 読み込み比 |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **57.7ms** | 64.6ms | **1.12x** | **0.28ms** | 55.7ms | **198.8x**\* |
| snappy | **49.3ms** | 208ms | **4.22x** | **14.0ms** | 103ms | **7.40x** |
| zstd | **61.0ms** | 325ms | **5.33x** | **20.6ms** | 185ms | **8.95x** |
| lz4 | **52.9ms** | 244ms | **4.61x** | **9.83ms** | 55.4ms | **5.64x** |

</div>

<br/>

#### 100 万行 vs fastparquet

<div align="center">

| コーデック | Carquet 書き込み | fastparquet 書き込み | 書き込み比 | Carquet 読み込み | fastparquet 読み込み | 読み込み比 |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **5.87ms** | 8.01ms | **1.36x** | **0.06ms** | 3.08ms | **51.3x**\* |
| snappy | **7.40ms** | 22.8ms | **3.08x** | **1.30ms** | 10.2ms | **7.85x** |
| zstd | **8.32ms** | 33.5ms | **4.03x** | **2.36ms** | 18.0ms | **7.61x** |
| lz4 | **7.70ms** | 22.4ms | **2.92x** | **1.06ms** | 5.37ms | **5.07x** |

</div>

<br/>

#### 10 万行 vs fastparquet

<div align="center">

| コーデック | Carquet 書き込み | fastparquet 書き込み | 書き込み比 | Carquet 読み込み | fastparquet 読み込み | 読み込み比 |
|-------|--------------|---------|---------|-------------|---------|---------|
| none | **0.87ms** | 1.08ms | **1.24x** | **0.02ms** | 0.58ms | **29.0x**\* |
| snappy | **1.52ms** | 2.43ms | **1.60x** | **0.38ms** | 1.20ms | **3.16x** |
| zstd | **1.90ms** | 3.47ms | **1.83x** | **0.64ms** | 2.06ms | **3.22x** |
| lz4 | **1.66ms** | 2.31ms | **1.39x** | **0.30ms** | 0.73ms | **2.43x** |

</div>

<br/>

#### 1000 万行 vs polars

<div align="center">

| コーデック | Carquet 書き込み | polars 書き込み | 書き込み比 | Carquet 読み込み | polars 読み込み | 読み込み比 | Carquet サイズ | polars サイズ |
|-------|--------------|---------|---------|-------------|---------|---------|------|------|
| none | 57.7ms | **37.5ms** | 0.65x | **0.28ms** | 28.7ms | **102.4x**\* | 190.7MB | 105.1MB |
| snappy | 49.3ms | **38.1ms** | 0.77x | **14.0ms** | 28.3ms | **2.03x** | 125.1MB | 103.8MB |
| zstd | 61.0ms | **48.6ms** | 0.80x | **20.6ms** | 42.1ms | **2.04x** | 95.3MB | 98.1MB |
| lz4 | 52.9ms | **40.0ms** | 0.76x | **9.83ms** | 29.4ms | **2.99x** | 122.7MB | 104.2MB |

</div>

<br/>

#### 100 万行 vs polars

<div align="center">

| コーデック | Carquet 書き込み | polars 書き込み | 書き込み比 | Carquet 読み込み | polars 読み込み | 読み込み比 | Carquet サイズ | polars サイズ |
|-------|--------------|---------|---------|-------------|---------|---------|------|------|
| none | 5.87ms | **4.53ms** | 0.77x | **0.06ms** | 3.14ms | **52.3x**\* | 19.1MB | 12.6MB |
| snappy | 7.40ms | **4.82ms** | 0.65x | **1.30ms** | 3.25ms | **2.50x** | 12.5MB | 11.5MB |
| zstd | 8.32ms | **5.83ms** | 0.70x | **2.36ms** | 4.12ms | **1.75x** | 9.6MB | 10.1MB |
| lz4 | 7.70ms | **4.77ms** | 0.62x | **1.06ms** | 3.26ms | **3.08x** | 12.3MB | 11.5MB |

</div>

<br/>

#### 10 万行 vs polars

<div align="center">

| コーデック | Carquet 書き込み | polars 書き込み | 書き込み比 | Carquet 読み込み | polars 読み込み | 読み込み比 | Carquet サイズ | polars サイズ |
|-------|--------------|---------|---------|-------------|---------|---------|------|------|
| none | 0.87ms | **0.72ms** | 0.83x | **0.02ms** | 0.33ms | **16.5x**\* | 1.9MB | 1.3MB |
| snappy | 1.52ms | **0.91ms** | 0.60x | **0.38ms** | 0.46ms | **1.21x** | 1.3MB | 1.1MB |
| zstd | 1.90ms | **1.05ms** | 0.55x | **0.64ms** | 1.08ms | **1.69x** | 1.0MB | 1.0MB |
| lz4 | 1.66ms | **0.96ms** | 0.58x | **0.30ms** | 0.44ms | **1.47x** | 1.2MB | 1.1MB |

</div>

<br/>

#### 同一ファイルのクロスリード（1000 万行）

両ライブラリが**同じ** Parquet ファイルを読みます。最も公平な同条件比較です。

<div align="center">

| コーデック | ライター | Carquet 読み込み | Arrow C++ 読み込み | 比率 |
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

#### 同一ファイルのクロスリード vs polars（1000 万行）

各ライブラリが、carquet のベンチマークが書くファイルと、polars が独自の設定（`id` と `category` に辞書エンコーディング）で書くファイルを読みます。

<div align="center">

| コーデック | ライター | Carquet 読み込み | polars 読み込み | 比率 | ファイル |
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

\* ゼロコピー mmap パス。fastparquet は pandas を介して読み込み、BYTE_STREAM_SPLIT やページチェックサムを書き込めないため、その行は直接比較しにくくなっています。

polars は独自の Rust 製 Parquet リーダーとライター（`use_pyarrow=False`、デフォルト）を使い、同じデータ、行グループサイズ、zstd レベル、統計情報を使い、`lz4` には LZ4_RAW を使います。BYTE_STREAM_SPLIT やページチェックサムは書き込めず、辞書エンコーディングを切り替える設定もありません。`id` と `category` を自動的に辞書エンコードするため、このデータでは非圧縮で carquet の辞書なしファイルより 45% 小さく、snappy / lz4 では 15–17% 小さくなり（zstd ではほぼ同じサイズ）、書き込みでコーデックとディスクを通るバイト数も少なくなります。同じ 2 カラムに辞書エンコーディングを有効にすると、carquet は 1000 万行を 45 ms（none）、51 ms（snappy）、51 ms（zstd）、48 ms（lz4）で書き込み、polars はそれぞれ 38、38、49、40 ms です。読み込みは同一ファイルの表で公平に比較できます。コーデックありの場合、carquet はどちらのライブラリのファイルも 2.0–4.1x 速く読みます。

</details>

<details>
<summary><b>Windows の完全な結果</b>（Intel Haswell VM、MSVC）</summary>

<br/>

*Carquet 0.7.3 を MSVC 19.44（x64 Release、OpenMP）でビルド、4 vCPU（Intel Haswell、3.1 GHz）の OpenStack VM、14.6GB RAM、仮想 SAS ディスク、Windows Server 2025、Arrow C++ 24.0.0（conda-forge）、ZSTD レベル 1、2026-11-02 計測。コールド読み込みは仮想ディスクへ行くため読み込み比は縮まります。書き込みとゼロコピーパスは CPU 律速です。*

#### 1000 万行 vs Arrow C++

<div align="center">

| コーデック | Carquet 書き込み | Arrow C++ 書き込み | 書き込み比 | Carquet 読み込み | Arrow C++ 読み込み | 読み込み比 | サイズ |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **182ms** | 382ms | **2.10x** | **2.04ms** | 126.5ms | **62.0x**\* | 190.7MB |
| snappy | **175ms** | 594ms | **3.39x** | **115.6ms** | 134.9ms | **1.17x** | 125.1MB |
| zstd | **242ms** | 927ms | **3.83x** | **140.7ms** | 159.2ms | **1.13x** | 95.3MB |
| lz4 | **178ms** | 527ms | **2.96x** | **95.5ms** | 99.9ms | **1.05x** | 122.9MB |

</div>

<br/>

#### 100 万行 vs Arrow C++

<div align="center">

| コーデック | Carquet 書き込み | Arrow C++ 書き込み | 書き込み比 | Carquet 読み込み | Arrow C++ 読み込み | 読み込み比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **15.6ms** | 37.6ms | **2.41x** | **0.44ms** | 13.72ms | **31.2x**\* |
| snappy | **25.4ms** | 61.2ms | **2.41x** | **14.01ms** | 16.31ms | **1.16x** |
| zstd | **30.4ms** | 94.3ms | **3.10x** | **17.09ms** | 18.38ms | **1.08x** |
| lz4 | **26.4ms** | 53.8ms | **2.04x** | 13.07ms | **12.89ms** | 0.99x |

</div>

<br/>

#### 10 万行 vs Arrow C++

<div align="center">

| コーデック | Carquet 書き込み | Arrow C++ 書き込み | 書き込み比 | Carquet 読み込み | Arrow C++ 読み込み | 読み込み比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **3.53ms** | 3.95ms | **1.12x** | **0.16ms** | 2.13ms | **13.3x**\* |
| snappy | **5.60ms** | 6.23ms | **1.11x** | **2.79ms** | 3.92ms | **1.41x** |
| zstd | **7.11ms** | 9.79ms | **1.38x** | **2.86ms** | 4.36ms | **1.52x** |
| lz4 | **5.75ms** | 5.82ms | **1.01x** | 2.22ms | **2.15ms** | 0.97x |

</div>

<br/>

#### 同一ファイルのクロスリード（1000 万行）

<div align="center">

| コーデック | ライター | Carquet 読み込み | Arrow C++ 読み込み | 比率 |
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

\* ゼロコピー mmap パス

</details>

<details>
<summary><b>ベンチマークを自分で実行する</b></summary>

<br/>

```bash
make benchmark                                    # Release build of the benchmark programs into build/
./build/benchmark_carquet                         # carquet alone
python3 benchmark/run_benchmark.py --skip-xlarge  # full comparison with PyArrow (+ Arrow C++ if built)
```

<div align="center">

| `run_benchmark.py` のフラグ / 変数 | 効果 |
|---|---|
| `--quick` | 1000 万行の構成のみ（グラフに表示されるもの） |
| `--skip-xlarge` | 1 億行の構成をスキップ（コーデックごとに約 2 GB、30 分以上） |
| `--fastparquet` | fastparquet とも比較（`pip install fastparquet pandas`） |
| `--polars` | polars 独自の Rust 製 Parquet リーダー / ライターとも比較（`pip install polars`） |
| `--no-pyarrow` / `--no-arrow-cpp` | PyArrow / Arrow C++ をスキップ |
| `--iterations N` | 反復回数を固定（デフォルト：自動、11〜51） |
| `--build-dir DIR` / `--temp-dir DIR` | バイナリの場所 / テストファイルの書き込み先 |
| `CARQUET_BENCH_ZSTD_LEVEL=3` | ZSTD レベル（デフォルト 1） |

</div>

ベンチマークは必ず **Release** ビルドで行ってください。`make benchmark` は `build/` を Release として再構成します。Debug の `build/` が残っていると、carquet が数倍遅く見えます。`make native` + `-DCARQUET_BUILD_BENCHMARKS=ON` でホストに最適化した数値が得られます。

fastparquet には同じデータ、行グループサイズ、辞書なし、zstd レベルを使い、`lz4` には LZ4_RAW を使います。fastparquet は BYTE_STREAM_SPLIT やページチェックサムを書き込めず、pandas を介して読み込むため、その数値は Arrow のものほど直接比較できません。

polars はネイティブのリーダーとライター（PyArrow を介さない）で実行し、同じデータ、行グループサイズ、zstd レベルを使い、`lz4` には LZ4_RAW を使います。BYTE_STREAM_SPLIT やページチェックサムは書き込めず、辞書エンコーディングを自分で選ぶ（切り替え設定なし）ため、このデータではファイルが小さくなります。

Arrow C++ との比較には専用のターゲットが必要です（システム上の Arrow C++ を検出するか、`CARQUET_ARROW_CPP_ROOT` で指定）：

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON \
  -DCARQUET_BUILD_BENCHMARKS=ON -DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON \
  [-DCARQUET_ARROW_CPP_ROOT=/path/to/arrow-prefix]
cmake --build build -j
```

低レベルの `parquet::ParquetFileReader` API（Arrow Table のマテリアライズなし）を並列行グループリーダーとともに使用します。

Windows（Developer PowerShell for VS）では、Arrow C++ を conda-forge から入手できます（`micromamba create -p C:\arrow -c conda-forge libarrow libparquet`）。ベンチマーク実行時にはその DLL が `PATH` 上にある必要があります：

```powershell
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCARQUET_BUILD_BENCHMARKS=ON `
  -DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON -DCARQUET_ARROW_CPP_ROOT=C:/arrow/Library
cmake --build build
$env:PATH = "C:\arrow\Library\bin;$env:PATH"
python benchmark\run_benchmark.py --build-dir build --skip-xlarge --no-pyarrow
```

ビルドフォルダと `--temp-dir` フォルダはリアルタイムのウイルス対策スキャンから除外してください。Windows Defender はベンチマークが書き込むすべてのファイルをスキャンします。

README のグラフは、実行後に `benchmark/` に残る JSON から、GitHub が閲覧者のテーマに応じて選ぶライト版とダーク版の 2 種類で描画されます：

```bash
python3 benchmark/plot_readme_chart.py benchmark/bench_<version>_<os>_<arch>_<date>.json   # writes res/img/bench_10m_m3_{light,dark}.png
```

</details>

<br/>

---

<!-- ────────────────────────────  クイックルック  ──────────────────────────── -->

## 30 秒でわかる carquet

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

> **carquet の位置づけ。** Arrow のエコシステム（計算カーネル、データセット、Flight、暗号化）が必要な場合は引き続き Arrow が適切な選択です。carquet は、Parquet を高速かつ正確に読み書きすることだけを必要とする C / C++ プログラム、組み込みシステム、言語ランタイムのためのものです。

<br/>

---

<!-- ────────────────────────────  インストール  ──────────────────────────── -->

## インストール

このセクションでは **C ライブラリ、そのヘッダー、`carquet` CLI** をインストールします。Python ユーザーは別の [Python API](#python-api)（`pip install carquet`）を利用してください。

**要件：** C11 コンパイラ（GCC 4.9+、Clang 3.4+、MSVC 2015+）と CMake 3.16+（または [xmake](#ソースから-xmake-で)）。zstd、zlib、lz4、brotli はシステムにあればそれを使い、なければ取得します。OpenMP は省略可能です（並列読み書き用）。

<br/>

### ソースから `make` で

最も簡単な方法です。`make` は適切なフラグで CMake を呼び出す薄いラッパーで、ほかに設定するものはありません：

```bash
git clone https://github.com/Vitruves/carquet.git && cd carquet
make                  # optimized Release build into build/
sudo make install     # library, headers, CLI, pkg-config and CMake package files
```

その後は `-lcarquet`、`pkg-config --cflags --libs carquet`、または CMake の `find_package(carquet)` でリンクします。

<br/>

### ソースから `xmake` で

[xmake](https://xmake.io) はほかに何もインストールする必要がありません。固定バージョンの zstd 1.5.6、zlib 1.4、lz4 1.10.0、brotli 1.1.0 をダウンロードして静的にリンクするため、`carquet` CLI には実行時の依存がありません。

```bash
git clone https://github.com/Vitruves/carquet.git && cd carquet
xmake                         # Release static library + CLI into build/<plat>/<arch>/release/
xmake install -o /usr/local   # headers, library and CLI (sudo if the prefix needs it)
```

静的インストールには `libcarquet.a` しか含まれないため、コーデックは自分でリンクしてください（`-lcarquet -lzstd -lz -llz4`、`--brotli=n` でビルドしていなければ `-lbrotlienc -lbrotlidec -lbrotlicommon` も）。代わりに `xmake f --shared=y` で構成すると、コーデックを内蔵した自己完結型の共有ライブラリになります。

<br/>

### CMake プロジェクトで使う

```cmake
include(FetchContent)
FetchContent_Declare(carquet GIT_REPOSITORY https://github.com/Vitruves/carquet.git GIT_TAG main)
FetchContent_MakeAvailable(carquet)
target_link_libraries(myapp PRIVATE carquet)
```

または `make install` の後に `find_package(carquet REQUIRED)` を使い、`carquet::carquet` をリンクします。

<br/>

---

<!-- ────────────────────────────  ビルド構成  ──────────────────────────── -->

## ビルド構成

### `make` ターゲット

どのターゲットも `build/` を構成し、並列にビルドします。`make help` でこの一覧を表示できます。いずれにも [xmake の同等コマンド](#xmake)があります。

<div align="center">

| コマンド | ビルドタイプ | 得られるもの |
|---|---|---|
| `make` / `make release` | Release | 静的ライブラリ + `carquet` CLI（デフォルト） |
| `make native` | Release | 同上、`-march=native` でこの CPU 向けに最適化（古い CPU では動作しないバイナリ） |
| `make shared` | Release | 静的ではなく共有ライブラリ（`.so` / `.dylib` / `.dll`） |
| `make benchmark` | Release | ベンチマークプログラム（`benchmark_carquet`、プロファイリングターゲット） |
| `make debug` | Debug | アサーション有効、最適化なしのライブラリ + CLI |
| `make test` | Debug | テストスイートをビルドし、`ctest` で実行 |
| `make dev` | Debug | すべて：テスト、サンプル、ベンチマーク、相互運用プログラム |
| `make fuzz` | Debug | libFuzzer ターゲット（Clang ツールチェーンを使用） |
| `make install` | Release | `release` をビルドしてから `PREFIX` へインストール |
| `make clean` | | ビルドディレクトリを削除 |

</div>

<div align="center">

| 変数 | デフォルト | 例 |
|---|---|---|
| `PREFIX` | CMake のデフォルト（Unix では `/usr/local`） | `make install PREFIX=$HOME/.local` |
| `BUILD` | `build` | `make test BUILD=build-debug`（複数の構成を並べて保持） |
| `JOBS` | コア数 | `make JOBS=4` |
| `CMAKE` | `cmake` | `make CMAKE=/opt/cmake/bin/cmake` |

</div>

デフォルトではすべてのターゲットが `build/` を共有するため、`make benchmark` の後に `make test` を実行すると Debug として再構成されます。Release と Debug のツリーを分けておくには `BUILD=` を使ってください。

<br/>

### CMake オプション

`make` ターゲットでカバーされないことは、CMake を直接呼び出し、`-DOPTION=ON|OFF` でオプションを渡します：

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_BUILD_SHARED=ON -DCARQUET_ENABLE_BROTLI=OFF
cmake --build build -j
sudo cmake --install build --prefix /usr/local
```

<div align="center">

| ビルド対象 | デフォルト | |
|---|:---:|---|
| `CARQUET_BUILD_CLI` | ON | `carquet` コマンドラインツール |
| `CARQUET_BUILD_SHARED` | OFF | 静的ではなく共有ライブラリ |
| `CARQUET_BUILD_DEV` | OFF | テスト、サンプル、ベンチマーク、相互運用をまとめて有効化 |
| `CARQUET_BUILD_TESTS` | = DEV | テストスイート（`ctest`） |
| `CARQUET_BUILD_EXAMPLES` | = DEV | [`examples/`](examples/) のプログラム |
| `CARQUET_BUILD_BENCHMARKS` | = DEV | ベンチマークとプロファイリングのプログラム |
| `CARQUET_BUILD_ARROW_CPP_BENCHMARK` | OFF | Arrow C++ 比較ベンチマーク（Arrow C++ が必要） |
| `CARQUET_BUILD_INTEROP` | = DEV | C / Python のラウンドトリッププログラム |
| `CARQUET_BUILD_FUZZ` | OFF | ファズターゲット |
| `CARQUET_BUILD_DOCS` | OFF | `docs` ターゲット、Doxygen HTML を `build/docs/html/` に出力 |

| 依存関係 | デフォルト | |
|---|:---:|---|
| `CARQUET_ENABLE_BROTLI` | ON | BROTLI コーデック。OFF で libbrotli への依存を外す |
| `CARQUET_BUNDLE_DEPS` | OFF | zstd / zlib / lz4 / brotli を常に取得して静的に同梱（ポータブルな成果物、Python ホイール） |
| `CARQUET_ZSTD_TARGET`, `_ZLIB_`, `_LZ4_`, `_BROTLI_TARGET` | 空 | 検索や取得の代わりに既存の CMake ターゲットをリンク（下記参照） |
| `CARQUET_ARROW_CPP_ROOT` | 空 | 比較ベンチマーク用の Arrow C++ インストールプレフィックス |
| `CMAKE_DISABLE_FIND_PACKAGE_OpenMP` | OFF | ON で OpenMP なしでビルド（実行時に `libgomp` / `vcomp140.dll` 不要） |

| CPU | デフォルト | |
|---|:---:|---|
| `CARQUET_NATIVE_ARCH` | OFF | `-march=native`。最速だがホスト専用のバイナリ |
| `CARQUET_ENABLE_SSE`, `_AVX`, `_AVX2`, `_AVX512` | ON | x86 SIMD カーネル。CPU の対応状況から実行時に選択 |
| `CARQUET_ENABLE_NEON` | ON | ARM NEON カーネル |
| `CARQUET_ENABLE_SVE` | OFF | ARM SVE カーネル（実験的） |
| `CARQUET_MSVC_HARDENING` | ON | MSVC：`/guard:cf` と `/CETCOMPAT` |

</div>

SIMD カーネルは組み込まれたうえで実行時に選択され、スカラーフォールバックもあるため、デフォルトのビルドはそのアーキテクチャのどの CPU でも動作します。バイナリをビルドマシンに縛るのは `CARQUET_NATIVE_ARCH` だけです。

<br/>

### xmake

`xmake.lua` は CMake の完全な代替で、オプションとデフォルトは同じです。`xmake f` で一度構成すれば、あとは `xmake` でビルドできます。構成は次の `xmake f` まで記憶されます。

<div align="center">

| `make` ターゲット | xmake の同等コマンド |
|---|---|
| `make` | `xmake` |
| `make native` | `xmake f --native_arch=y && xmake` |
| `make shared` | `xmake f --shared=y && xmake` |
| `make benchmark` | `xmake f --benchmarks=y && xmake` |
| `make debug` | `xmake f -m debug && xmake` |
| `make test` | `xmake f -m debug --tests=y && xmake && xmake test` |
| `make dev` | `xmake f -m debug --dev=y && xmake` |
| `make fuzz` | `xmake f -m debug --fuzz=y --toolchain=clang && xmake` |
| `make clean` | `xmake clean`（構成もリセットするには `xmake f -c`） |
| `cmake --build build --target docs` | `xmake docs`（configure フラグは不要） |

</div>

オプションは CMake のものと対応しており、`CARQUET_` プレフィックスを除いた小文字で `--option=y|n` の形式です。自由に組み合わせられます：`xmake f -m release --shared=y --brotli=n --avx512=n`。

<div align="center">

| xmake オプション | デフォルト | 説明 |
|--------------|---------|-------------|
| `--dev` | n | テスト、サンプル、ベンチマーク、相互運用をビルド |
| `--tests` / `--examples` / `--benchmarks` / `--interop` | n | グループを個別にビルド |
| `--arrow_cpp_benchmark` | n | Arrow C++ 比較ベンチマーク（利用可能なら xmake リポジトリの Arrow を使用） |
| `--cli` | y | `carquet` CLI ツールをビルド |
| `--shared` | n | 静的ライブラリの代わりに共有ライブラリをビルド |
| `--openmp` | y | OpenMP によるカラムの並列読み込み（利用できない場合は自動的に無効化） |
| `--native_arch` | n | 最大性能のための `-march=native`（ホスト専用バイナリ） |
| `--sse` / `--avx` / `--avx2` / `--avx512` / `--neon` | y | SIMD 命令セット（自動検出） |
| `--sve` | n | ARM SVE（実験的） |
| `--brotli` | y | BROTLI コーデック（xmake リポジトリの libbrotli） |
| `--fuzz` | n | ファズターゲットをビルド（`--toolchain=clang` を使用）。続けて `--fuzz_engine=libFuzzer\|AFL`、`--fuzz_sanitizers=y\|n` |

</div>

<details>
<summary><b>Windows</b>（MSVC、clang-cl、MinGW-w64）</summary>

<br/>

MSVC（x64 と Win32）、clang-cl、MinGW-w64 GCC のいずれにも、静的・共有の両方で対応しています。`CMakePresets.json` には `msvc`、`msvc-shared`、`clang-cl`、`mingw` のプリセットがあります：

```powershell
# from a Developer PowerShell for VS (or with MinGW-w64 gcc + ninja on PATH for "mingw")
cmake --preset msvc
cmake --build --preset msvc
ctest --preset msvc
```

- 手動で構成する場合は、C コンパイラだけでなく C++ コンパイラも指定してください（例：`CC=clang-cl CXX=clang-cl`）。取得される zstd プロジェクトは C++ を有効にするため、指定しないと CMake が不一致のドライバを選んでしまいます。
- `-DCARQUET_BUILD_SHARED=ON` では `carquet.dll` とインポートライブラリ `carquet.lib`（MinGW では `libcarquet.dll` / `libcarquet.dll.a`）が得られます。利用側は `CARQUET_BUILD_SHARED` を定義する必要があります。`find_package(carquet)` と `carquet.pc` は自動で追加します。DLL は実行ファイルの隣か `PATH` 上に置いてください。
- `carquet_writer_get_buffer()` の出力など、carquet から受け取ったメモリは `free()` ではなく `carquet_free()` で解放してください。DLL とプログラムが異なる C ランタイムヒープ（Debug と Release の CRT、`/MT`、MinGW と MSVC）を使っている場合があります。
- ファイルパスには UTF-8 を使えます（非 ASCII の名前も動作します）。`carquet` CLI は Unicode の引数を受け付け、コンソールに UTF-8 で出力します。
- OpenMP 付きの MinGW ビルドには実行時に `libgomp-1.dll` と `libwinpthread-1.dll` が必要です。MSVC ビルドには `vcomp140.dll`（Visual C++ 再頒布可能パッケージに含まれる）が必要です。その依存のないビルドには `-DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=ON` を渡してください。
- Windows Defender の機械学習ヒューリスティックが、`test_encoding_roundtrip.exe` と `test_fuzz_regressions.exe` を `Trojan:Win64/CobaltStrike.RDH!MTB` として隔離することがあります（埋め込まれたバイナリのテストデータによる誤検知で、ctest はこれらを「Not Run」と報告します）。スイート全体を実行するには、ビルドディレクトリをスキャン対象から除外してください。

</details>

<details>
<summary><b>独自の zstd / zlib / lz4 / brotli を使う</b>（より大きな CMake プロジェクトへ carquet を組み込む場合）</summary>

<br/>

carquet を `add_subdirectory()` または `FetchContent` で取り込むと、各圧縮ライブラリは次の順序で解決されます：

1. `CARQUET_ZSTD_TARGET`、`CARQUET_ZLIB_TARGET`、`CARQUET_LZ4_TARGET`、`CARQUET_BROTLI_TARGET` で指定されたターゲット（設定されている場合）。
2. システムにインストールされたコピー（`CARQUET_BUNDLE_DEPS` が ON の場合はスキップ）。
3. carquet 自身が取得してビルドするコピー。プロジェクトが既にそのライブラリを通常のターゲット名（`libzstd_static`、`zlibstatic`、`lz4_static`）でビルドしている場合、carquet は 2 つ目のコピーを取得する代わりにそのターゲットを再利用します。

したがって、これらのライブラリをベンダリングしておりシステムにコピーがないプロジェクトでは、設定は不要です。明示的に選択するには、ターゲットを名指しします：

```cmake
set(CARQUET_ZLIB_TARGET zlibstatic)
add_subdirectory(carquet)
```

名指しされたターゲットは検索も取得も行わずそのままリンクされ、ターゲットではない名前は暗黙のフォールバックではなく configure エラーになります。ターゲットはインポート済みでもツリー内ビルドでも、静的でも共有でもかまいませんが、インクルードディレクトリを usage requirement として持っている必要があります。

</details>

<details>
<summary><b>テスト、サニタイザ、API ドキュメント</b></summary>

<br/>

```bash
make test                                    # Debug build + ctest
make dev && ctest --test-dir build           # everything, including examples and interop programs
cmake -B build -DCARQUET_BUILD_DOCS=ON && cmake --build build --target docs   # Doxygen → build/docs/html/
```

[CONTRIBUTING.md](CONTRIBUTING.md) に、プルリクエストごとのゲート（テスト、サニタイザ、ファジング、相互運用テスト、ベンチマーク）が記載されています。

</details>

<a id="ai-coding-agents"></a>

<details>
<summary><b>AI コーディングエージェント</b>（Claude Code と Codex 向けのすぐに使えるスキル）</summary>

<br/>

[`.agents/skills/carquet/`](.agents/skills/carquet/) は、コーディングエージェントに carquet をプロジェクトへ統合する方法を教える [Agent Skill](https://agentskills.io) です：ビルドの配線（CMake、FetchContent、pkg-config、xmake）、何にどの API を使うか、完全なレシピ、レベル / null ビットマップ / ライフタイムのルール、型マッピング、よくある落とし穴とその対処法。自分のリポジトリにコピーまたはリンクすれば、エージェントが自動的に認識します：

```bash
mkdir -p .claude/skills && cp -r path/to/carquet/.agents/skills/carquet .claude/skills/   # Claude Code
mkdir -p .agents/skills && cp -r path/to/carquet/.agents/skills/carquet .agents/skills/   # Codex CLI
ln -s ../../third_party/carquet/.agents/skills/carquet .claude/skills/carquet               # vendored: stays in sync
```

同じファイルは個人用スキル（`~/.claude/skills/`、`~/.codex/skills/`）としても機能します。

</details>

<br/>

---

<!-- ────────────────────────────  C API  ──────────────────────────── -->

## C API

以下に書き込みと読み込みの例を示します。それ以外はすべて [`docs/` のマニュアル](docs/README.md)で扱っています。

<br/>

### Parquet ファイルの書き込み

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

### Parquet ファイルの読み込み

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

### その他のレシピ

フラットな読み書き以外はすべてマニュアルにあります。各ページは実行可能なサンプルにリンクしています：

<div align="center">

| やりたいこと | 参照先 |
|---|---|
| nullable カラム、行グループ、バッファ出力 | [`docs/writing.md`](docs/writing.md) |
| リスト、マップ、グループ、定義/反復レベル | [`docs/nested-data.md`](docs/nested-data.md) |
| カラム射影、統計、メタデータの検査 | [`docs/reading.md`](docs/reading.md) |
| 述語プッシュダウン、ページレベルのフィルタリング | [`docs/reading.md`](docs/reading.md) |
| 既存ファイルへの行グループの追記 | [`docs/writing.md`](docs/writing.md#append-to-an-existing-file) |
| 圧縮、カスタムコーデック、ライターのチューニング | [`docs/writing.md`](docs/writing.md), [`docs/performance.md`](docs/performance.md) |
| mmap、ゼロコピー、プリバッファリング、I/O 結合 | [`docs/performance.md`](docs/performance.md) |
| エラーコードと復旧のヒント | [`docs/error-handling.md`](docs/error-handling.md) |

</div>

<br/>

### API リファレンス

完全な API は [`include/carquet/carquet.h`](include/carquet/carquet.h) にあります。主な型：

<div align="center">

| 型 | 用途 |
|------|---------|
| `carquet_reader_t` | ファイルリーダー（パス、`FILE*`、またはメモリバッファから開く） |
| `carquet_writer_t` | ファイルライター |
| `carquet_batch_reader_t` | 高レベルのバッチ反復 |
| `carquet_schema_t` | スキーマの定義とイントロスペクション |
| `carquet_error_t` | 詳細なエラー情報（コード、メッセージ、ソース位置、復旧のヒント） |

</div>

[マニュアル](docs/README.md)は、どの API をいつ使うべきかを説明しています。ソースのレイアウトとアーキテクチャは [CONTRIBUTING.md](CONTRIBUTING.md) に記載されています。ブラウズ可能な HTML リファレンスは Doxygen でヘッダーから生成できます。[ビルド構成](#ビルド構成)の *テスト、サニタイザ、API ドキュメント* を参照してください。

<details>
<summary><b>サンプルアプリケーション：MetricStore</b>（<code>mocklib/</code>）</summary>

<br/>

[`mocklib/`](mocklib/) は **MetricStore** という、carquet 上に構築された完全で自己完結したサンプルアプリケーションです。実際の使用例のリファレンスとしても、公開 API のエンドツーエンド統合テストとしても役立ちます。

時系列テレメトリストア（イベントを取り込み、検査してクエリする）をモデル化しており、その過程で **carquet の公開関数 143 個のうち 117 個（約 82 %）** を使用します：

- **書き込み側**：論理型を使ったスキーマ構築、カラムごとのエンコーディング / 圧縮 / ブルームフィルタ調整、ページインデックスと統計
- **読み込み側**：カラム射影、述語プッシュダウン、ページレベルのフィルタリング、ブルームフィルタによるメンバーシップ判定、Arrow C Data Interface ブリッジ、ネストした `LIST` の再構築、バッファ/追記 I/O、メタデータのイントロスペクション

下流のプロジェクトと同じ方法で carquet をリンクし（`find_package(carquet)` / `add_subdirectory`）、自己検証するラウンドトリップテストを同梱し、`write_sample` バイナリは外部検査用のファイルを出力します。[`mocklib/README.md`](mocklib/README.md) を参照してください。

</details>

<br/>

---

<!-- ────────────────────────────  PYTHON API  ──────────────────────────── -->

## Python API

PyPI の `carquet` は上記の C 呼び出しのラッパーではありません。同じエンジン上に構築された、**Python らしい `pyarrow.parquet` 風の API**（`read_table`、`ParquetFile`、`ParquetWriter`、DNF `filters`）です。ホイールにはすべてのコーデックを組み込んだ C ライブラリが同梱されています。データは Arrow C Data Interface を通じて Python に渡されるため、デコード、展開、I/O は GIL なしで C の中で実行されます。**依存パッケージはありません**。読み取り結果は `carquet.Table` で、それ自体で Python の値、NumPy、pandas に変換でき、polars、DuckDB、PyArrow にはコピーなしで渡せます。

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

| やりたいこと | 使うもの |
|---|---|
| ファイル、ディレクトリ、バイト列、ファイルオブジェクトを読む | `cq.read_table`, `cq.read_pandas` |
| 大きなファイルをストリーミングする | `ParquetFile.iter_batches`, `cq.iter_batches` |
| 統計、ページインデックス、ブルームフィルタでプルーニングする | `filters=`（DNF、厳密な結果）、`ParquetFile.bloom_filter` |
| テーブル、DataFrame、polars、DuckDB の結果を書き込む | `cq.write_table`, `cq.ParquetWriter`（インクリメンタル、`append=True`） |
| pandas、NumPy、polars、DuckDB、PyArrow にデータを渡す | `t.to_pandas()`, `t.column("x").to_numpy()`, `polars.DataFrame(t)`, `pyarrow.table(t)` |
| データを読まずにファイルを検査する | `cq.read_metadata`, `cq.read_parquet_schema`, `cq.info`, `cq.validate` |

</div>

型は `pyarrow.parquet` と正確にラウンドトリップします（タイムゾーン、カテゴリカル、デシマル、large 型、ネストされたデータ）。完全なリファレンスは [`python/README.md`](python/README.md) にあります。

<br/>

---

<!-- ────────────────────────────  CLI ツール  ──────────────────────────── -->

## CLI ツール

Carquet には、Parquet ファイルの検査と C リーダーコードの生成を行うコマンドラインツールが付属しています。ライブラリと一緒にデフォルトでビルド・インストールされます。

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

典型的な使い方：

```bash
carquet schema data.parquet
carquet head -n 20 data.parquet
carquet stat data.parquet
carquet validate data.parquet
```

<br/>

### 行フィルタリング

`cat`、`count`、`head`、`tail`、`sample`、`export` は `-p / --filter EXPR` を受け付け、行述語をページレベルまでプッシュダウンします。カラムインデックスの min/max が述語に一致しうるページだけが展開されます。

```bash
carquet cat -p "price > 100 AND status = 'active'" data.parquet
carquet count --filter "id >= 1000" data.parquet
carquet export --filter "ts IS NOT NULL" -c id,ts data.parquet
```

文法は `column OP value [AND column OP value]...` で、`OP` ∈ {`=`, `==`, `!=`, `<>`, `<`, `<=`, `>`, `>=`} です。加えて `column IS NULL` / `column IS NOT NULL` が使えます。

> フィルタリングにはファイルにページインデックスが必要です（`write_page_index = true`）。

<br/>

### コード生成

任意の Parquet ファイルのスキーマから、完全でコンパイル可能な C リーダーを生成します：

```bash
carquet codegen -f data.parquet -o reader.c
# Generated: reader.c
# Compile:   clang -o reader reader.c -I.../include -L.../build -lcarquet ...

./reader                    # reads data.parquet (embedded as default)
./reader other.parquet      # override with a different file
```

<details>
<summary><b>codegen オプション</b></summary>

<br/>

<div align="center">

| フラグ | 説明 |
|------|-------------|
| `-f`, `--file FILE` | スキーマを調べる Parquet ファイル |
| `-o`, `--output FILE` | 出力ソースファイル（デフォルト：stdout） |
| `--mmap` | 生成コードでメモリマップド I/O を使用 |
| `--skeleton` | カスタムロジック用に空の `process_batch` を生成 |
| `-c`, `--columns COLS` | カンマ区切りのカラムフィルタ |
| `-b`, `--batch-size N` | バッチサイズ（デフォルト：1024） |

</div>

</details>

<br/>

---

<!-- ────────────────────────────  機能サポート  ──────────────────────────── -->

## Parquet 機能サポート

ライブラリが現在実装している内容を、領域別に示します。

<br/>

### フォーマット

<div align="center">

| 機能 | 状況 |
|---------|--------|
| 物理型 | ✓  全 8 種（BOOLEAN から FIXED_LEN_BYTE_ARRAY まで） |
| 論理型 | ✓  STRING, DATE, TIME, TIMESTAMP, DECIMAL, UUID, JSON, INTERVAL, FLOAT16, VARIANT, GEOMETRY, GEOGRAPHY |
| エンコーディング | ✓  PLAIN, RLE, DICTIONARY, DELTA_BINARY_PACKED, DELTA_LENGTH_BYTE_ARRAY, DELTA_BYTE_ARRAY, BYTE_STREAM_SPLIT（読み書き） |
| データページバージョン | ✓  V1（デフォルト）と V2（読み書き） |
| 圧縮 | ✓  UNCOMPRESSED, SNAPPY, GZIP, LZ4（Hadoop フレーム、コーデック 5）, LZ4_RAW（コーデック 7）, ZSTD, BROTLI（`CARQUET_ENABLE_BROTLI`、デフォルト ON）。LZO は `carquet_register_codec()` 経由 |
| ネストされたスキーマ | ✓  定義/反復レベルを持つグループ、リスト、マップ。書き込み時の単一レベル LIST/MAP 自動シュレッディング（`carquet_writer_write_list_column`）と読み込み時の List 再構築。ネストされた `ARROW:schema` の出力 |
| 暗号化 | ✗ 未対応。暗号化フッターのファイル（`PARE`）は、明示的なメッセージ付きで `CARQUET_ERROR_NOT_IMPLEMENTED` として報告されます |

</div>

<details>
<summary><b>インデックスとプルーニング</b></summary>

<br/>

<div align="center">

| 機能 | 状況 |
|---------|--------|
| ブルームフィルタ | ✓  読み込み、書き込み、クエリ（`carquet_bloom_filter_check_*`）。行グループのプルーニングに自動的に利用 |
| ページインデックス | ✓  カラムインデックス + オフセットインデックス（読み書き + ページごとの統計アクセス） |
| 統計 | ✓  カラムチャンクごとの min/max/null 数。辞書カラムの正確な `distinct_count`。Parquet 2.9 SizeStatistics（非エンコードのバイト配列バイト数 + レベルヒストグラム） |
| 述語プッシュダウン | ✓  統計 + ブルームフィルタによる自動行グループプルーニング（コールバック不要）。カラムインデックスによるページレベルのフィルタリング（`carquet_batch_reader_set_page_filter`） |
| 地理空間統計 | ✓  GEOMETRY/GEOGRAPHY のカラムごとのバウンディングボックスとジオメトリ型（`carquet_reader_geospatial_statistics`） |
| カラム射影 | ✓  選択したカラムのみ読み込み |
| 行ウィンドウ | ✓  行の `[offset, limit)` スライスを読み込み、ページを丸ごとスキップ（`carquet_batch_reader_set_row_range`） |

</div>

</details>

<details>
<summary><b>I/O と拡張性</b></summary>

<br/>

<div align="center">

| 機能 | 状況 |
|---------|--------|
| 追記 | ✓  既存ファイルへ行グループを追加（`carquet_writer_open_append`） |
| カスタムコーデック | ✓  コーデックスロットごとにカスタムの圧縮/展開実装を登録（`carquet_register_codec`） |
| キー・バリューメタデータ | ✓  任意のフッターメタデータの読み書き |
| フィールドごとのメタデータ | ✓  `ARROW:schema` 経由の Arrow `Field.custom_metadata`（変数ラベル/説明）（読み書き） |
| カラムごとのオプション | ✓  カラムごとのエンコーディング、圧縮、統計、ブルームフィルタ |
| バッファライター | ✓  Parquet をメモリ内バッファへ書き込み |
| CRC32 | ✓  ページレベルの検証（ARM ではハードウェアアクセラレーション） |
| メモリマップド I/O | ✓  非圧縮 PLAIN データのゼロコピー読み込み |
| 辞書保持読み込み | ✓  マテリアライズされた値の代わりに RLE インデックス + 辞書を返す（`preserve_dictionaries`、`carquet_row_batch_column_dictionary`） |
| I/O 結合 | ✓  複数カラムの読み込みを単一の I/O でプリバッファ |
| 投機的フッター読み込み | ✓  ほとんどのファイルで単一 I/O でのオープン |
| OpenMP 並列読み込み | ✓  利用可能な場合 |
| 並列書き込み | ✓  カラム横断のページレベル並列エンコード/圧縮（OpenMP）+ バックグラウンドの行グループ I/O（`async_io`） |

</div>

</details>

<br/>

---

<!-- ────────────────────────────  相互運用性  ──────────────────────────── -->

## 相互運用性

Carquet のファイルは、PyArrow、DuckDB、Spark、およびあらゆる Parquet リーダーと完全に互換です。

```python
import pyarrow.parquet as pq
table = pq.read_table("carquet_output.parquet")  # just works
```

```sql
-- DuckDB
SELECT * FROM read_parquet('carquet_output.parquet');
```

<details>
<summary><b>相互運用テストスイート</b></summary>

<br/>

双方向のラウンドトリップ（C が書き Python が読む、およびその逆）：

```bash
cmake -B build -DCARQUET_BUILD_INTEROP=ON && cmake --build build
python3 interop/run_interop.py
```

大きなマルチ行グループ Snappy ファイル（多数の行グループ、幅広い `BYTE_ARRAY` カラム、2 GiB を超えるページオフセット）には専用のジェネレータがあり、PyArrow でファイルを構築して carquet がすべてのページをデコードできることを確認します：

```bash
CARQUET_BIN=build/carquet python3 interop/snappy_multi_rowgroup.py --rows 150000 --row-groups 92
```

</details>

<br/>

---

<!-- ────────────────────────────  コントリビューションとライセンス  ──────────────────────────── -->

## コントリビューション、引用、ライセンス

- **コントリビューション**：[CONTRIBUTING.md](CONTRIBUTING.md) に、セットアップ、コードスタイル、プルリクエストごとのゲート（テスト、サニタイザ、ファジング、相互運用テスト、ベンチマーク）が記載されています。リリースノートは [CHANGELOG.md](CHANGELOG.md) にあります。
- **引用**：[`CITATION.cff`](CITATION.cff) を同梱しています。GitHub の「Cite this repository」ボタンで BibTeX または APA 形式に整形できます。
- **ライセンス**：MIT、[LICENSE](LICENSE) を参照してください。
