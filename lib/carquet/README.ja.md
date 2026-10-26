<h1 align="center">
  <img src="res/img/carquet_banner.png" alt="Carquet" width="840" />
</h1>

<p align="center">
  Apache Parquet ファイルを読み書きするための、高速なピュア C ライブラリです。
</p>

<p align="center">
  <a href="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml"><img src="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml/badge.svg" alt="ビルド" /></a>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue" alt="プラットフォーム" />
  <img src="https://img.shields.io/badge/C-C11-blue" alt="C 規格" />
  <a href="https://github.com/Vitruves/carquet/releases"><img src="https://img.shields.io/github/v/release/Vitruves/carquet?color=blue" alt="リリース" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="ライセンス" /></a>
  <a href="https://github.com/Vitruves/carquet/stargazers"><img src="https://img.shields.io/github/stars/Vitruves/carquet?style=flat&color=blue" alt="スター" /></a>
  <a href="https://github.com/Vitruves/carquet/graphs/traffic"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/Vitruves/carquet/traffic/clones.json" alt="クローン数" /></a>
  <br/>
  <img src="https://img.shields.io/badge/SIMD-SSE4.2%20%C2%B7%20AVX2%20%C2%B7%20AVX--512%20%C2%B7%20NEON%20%C2%B7%20SVE-555" alt="SIMD" />
  <img src="https://img.shields.io/badge/deps-zstd%20%C2%B7%20zlib%20%C2%B7%20lz4%20%C2%B7%20brotli-555" alt="依存関係" />
</p>

<p align="center">
  <a href="#ハイライト">ハイライト</a> &nbsp;·&nbsp;
  <a href="#インストール">インストール</a> &nbsp;·&nbsp;
  <a href="#c-api">C API</a> &nbsp;·&nbsp;
  <a href="#cli-ツール">CLI ツール</a> &nbsp;·&nbsp;
  <a href="#ベンチマーク">ベンチマーク</a> &nbsp;·&nbsp;
  <a href="#parquet-機能サポート">機能サポート</a> &nbsp;·&nbsp;
  <a href="#相互運用性">相互運用性</a> &nbsp;·&nbsp;
  <a href="#ready-for-ai-coding-agents">AI エージェント</a> &nbsp;·&nbsp;
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

<!-- ────────────────────────────  HIGHLIGHTS  ──────────────────────────── -->

## ハイライト

**小さく自己完結**

- **ピュア C11**。外部依存は 4 つ（zstd、zlib、lz4、brotli）のみで、見つからない場合は CMake がすべて自動取得します。Brotli は無効化できます
- **ライブラリ約 500 KB**（strip 済み Release ビルド、コーデックは別途リンク）。Arrow は約 50 MB 以上
- **CLI 内蔵**。ファイルの検査（`schema`、`info`、`head`、`tail`、`stat`、…）と C コード生成（`codegen`）

**高速**

- 非圧縮データの読み込みは Arrow C++ の **39 倍**（mmap ゼロコピー）、PyArrow の **100 倍**
- 同一ファイルの圧縮データ読み込みは Arrow C++ の **1.3〜2.1 倍**（クロスリードベンチマーク、Apple M3。x86 では 1.6〜2.6 倍）
- 書き込みは Apple M3 で Arrow C++ の **1.8〜5.7 倍**（旧い x86 の計測では 1.0〜1.6 倍）
- 1000 万行の非圧縮データを **0.36 ms** で読み込み（Apple M3 での mmap ゼロコピー）
- SIMD 最適化（SSE4.2、AVX2、AVX-512、NEON、SVE）。実行時検出とスカラーフォールバック付き

**完全で互換性が高い**

- 幅広い Parquet カバレッジ：すべての物理型とエンコーディング、Snappy / GZIP / LZ4 / ZSTD / Brotli、ネストされたスキーマ、ブルームフィルタ、ページインデックス（[全対応表](#parquet-機能サポート)）
- PyArrow、DuckDB、Spark とそのまま互換
- [Arrow C Data Interface](docs/reading.md#export-to-arrow-c-data-interface) ブリッジ（`carquet_arrow_export_*` / `carquet_arrow_import_*` / `carquet_writer_write_arrow` / `carquet_reader_read_arrow`）により、依存ゼロかつコピー最小で Arrow エコシステムとデータ交換でき、**任意の深さのネスト読み書き**（任意の深さの struct / list / map）に対応

<a id="ready-for-ai-coding-agents"></a>

**AI コーディングエージェント対応**

- [`.agents/skills/carquet/`](.agents/skills/carquet/) は [Agent Skill](https://agentskills.io)（Claude Code、Codex CLI）で、carquet をプロジェクトに統合する方法をエージェントに教えます：ビルドの配線、タスクごとに使う API、完全なレシピ、レベル / null ビットマップ / ライフタイムのルール、型マッピング、よくある落とし穴。リポジトリにコピーすれば、エージェントが自動的に利用します：

  ```bash
  cp -r path/to/carquet/.agents/skills/carquet .claude/skills/   # Claude Code
  cp -r path/to/carquet/.agents/skills/carquet .agents/skills/   # Codex CLI
  ```

<br/>

### 30 秒でわかる carquet

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
printf("%lld rows\n", (long long)carquet_reader_num_rows(r));   /* 3 行 */
carquet_reader_close(r);
carquet_schema_free(schema);
```

バッチと nullable カラムを含む完全な書き込み・読み込みループは [C API](#c-api) にあります。

<br/>

### carquet の位置づけ

carquet は小さく組み込み可能な選択肢です。自動取得される 4 つの依存関係を持つ約 500 KB の C11 ライブラリで、1 分以内にビルドできます。対して Arrow C++ は 50 MB 以上、PyArrow のホイールは約 100 MB あります。下記の 1000 万行ベンチマークでは、圧縮ファイルの書き込みは Arrow C++ より 1.8〜5.7 倍速く、まったく同じファイルの読み込みは 1.3〜2.1 倍速く、非圧縮データは mmap から直接 0.36 ms で提供されます。Arrow のエコシステム（計算カーネル、データセット、Flight、暗号化）が必要な場合は引き続き Arrow が適切な選択です。carquet は、Parquet を高速かつ正確に読み書きすることだけを必要とする C / C++ プログラム、組み込みシステム、言語ランタイムのためのものです。

<br/>

---

<!-- ────────────────────────────  INSTALLATION  ──────────────────────────── -->

## インストール

Carquet は CMake（主）または [xmake](https://xmake.io) で、Linux、macOS、Windows 上でビルドできます。

### 要件

- C11 コンパイラ（GCC 4.9+、Clang 3.4+、MSVC 2015+）
- CMake 3.16+（または [xmake](https://xmake.io)。[xmake でのビルド](#xmake-でのビルド)を参照）
- zstd、zlib、lz4、brotli（見つからない場合は自動取得。Brotli は省略可、`-DCARQUET_ENABLE_BROTLI=OFF`）
- OpenMP（省略可、カラムの並列読み書き用）

<br/>

### クイックスタート（make）

`make` ラッパーは、最適化された CMake ビルドと `/usr/local` へのインストールを実行します：

```bash
git clone https://github.com/Vitruves/carquet.git
cd carquet
make                              # 最適化ビルド（`make help` で全ターゲットを一覧表示）
sudo make install                 # /usr/local へインストール（PREFIX=/opt/carquet で変更可）
```

<br/>

### ビルドとインストール（CMake）

特定のビルドオプションやカスタムプレフィックスが必要な場合は、CMake を直接呼び出します：

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release   # オプションを追加、例：-DCARQUET_BUILD_SHARED=ON
cmake --build build -j$(nproc)
sudo cmake --install build --prefix /usr/local
```

どちらの方法でも以下がインストールされます：

- `libcarquet.a`（`-DCARQUET_BUILD_SHARED=ON` の場合は `.so` / `.dylib`）
- `include/carquet/` ヘッダー
- `carquet` CLI バイナリ
- `carquet.pc`（pkg-config）と `find_package(carquet)` 用の CMake パッケージ設定

インストール後は `-lcarquet` でリンクするか、`pkg-config --cflags --libs carquet` でフラグを解決してください。

<br/>

### Python（pip）

`carquet` ホイールはライブラリを同梱し、Arrow C Data Interface を通じて公開します。そのため、Parquet デコードロジックを重複させることなく、テーブルを PyArrow と相互にやり取りできます：

```bash
pip install carquet
```

```python
import carquet
table = carquet.read_table("data.parquet")                 # pyarrow.Table
carquet.write_table(table, "out.parquet", compression="zstd")
print(carquet.info("out.parquet").num_rows)
```

パッケージは [`python/`](python/) にあります（scikit-build-core、ctypes、PyArrow は実行時に省略可）。API とホイールのビルドについては、その [README](python/README.md) を参照してください。

<br/>

### ビルドオプション

CMake の configure ステップに渡します（`-DOPTION=ON|OFF`）：

<div align="center">

| オプション | デフォルト | 説明 |
|--------|---------|-------------|
| `CARQUET_BUILD_DEV` | OFF | すべてをビルド（テスト、サンプル、ベンチマーク、相互運用テスト） |
| `CARQUET_BUILD_TESTS` | OFF | テストスイートのみビルド |
| `CARQUET_BUILD_CLI` | ON | `carquet` CLI ツールをビルド |
| `CARQUET_BUILD_SHARED` | OFF | 静的ライブラリの代わりに共有ライブラリをビルド |
| `CARQUET_NATIVE_ARCH` | OFF | 最大性能のための `-march=native`（ホスト専用バイナリ） |
| `CARQUET_ENABLE_SVE` | OFF | ARM SVE（実験的） |
| `CARQUET_ENABLE_BROTLI` | ON | BROTLI コーデック（libbrotli を検出または取得。OFF で依存を外す） |

</div>

x86 SIMD（SSE、AVX、AVX2、AVX-512）と ARM NEON はすべて自動検出され、デフォルトで有効です。

<details>
<summary><b>すべてのビルドオプション</b></summary>

<br/>

<div align="center">

| オプション | デフォルト | 説明 |
|--------|---------|-------------|
| `CARQUET_BUILD_EXAMPLES` | OFF | サンプルプログラムをビルド |
| `CARQUET_BUILD_BENCHMARKS` | OFF | ベンチマークとプロファイリングプログラムをビルド |
| `CARQUET_BUILD_ARROW_CPP_BENCHMARK` | OFF | オプションの Arrow C++ 比較ベンチマーク |
| `CARQUET_BUILD_INTEROP` | OFF | 相互運用テストをビルド |
| `CARQUET_BUILD_FUZZ` | OFF | ファズターゲットをビルド |
| `CARQUET_BUILD_DOCS` | OFF | `docs` ターゲットを追加（Doxygen が必要） |
| `CARQUET_BUNDLE_DEPS` | OFF | zstd/zlib/lz4/brotli を常に取得して静的に同梱（ポータブルなリリース成果物向け） |
| `CARQUET_ZSTD_TARGET` / `CARQUET_ZLIB_TARGET` / `CARQUET_LZ4_TARGET` / `CARQUET_BROTLI_TARGET` | *（空）* | 検索や取得の代わりにリンクする、そのライブラリの既存 CMake ターゲット |
| `CARQUET_ARROW_CPP_ROOT` | *（空）* | 比較ベンチマーク用のカスタム Arrow C++ インストールへのパス |
| `CARQUET_ENABLE_SSE` | ON | SSE 最適化（x86、自動検出） |
| `CARQUET_ENABLE_AVX` | ON | AVX 最適化（x86、自動検出） |
| `CARQUET_ENABLE_AVX2` | ON | AVX2 最適化（x86、自動検出） |
| `CARQUET_ENABLE_AVX512` | ON | AVX-512 最適化（x86、自動検出） |
| `CARQUET_ENABLE_NEON` | ON | NEON 最適化（ARM、自動検出） |

</div>

</details>

<details>
<summary><b>独自の zstd / zlib / lz4 を使う</b>（より大きな CMake プロジェクトへ carquet を組み込む場合）</summary>

<br/>

carquet を `add_subdirectory()` または `FetchContent` で取り込むと、各圧縮ライブラリは次の順序で解決されます：

1. `CARQUET_ZSTD_TARGET`、`CARQUET_ZLIB_TARGET`、`CARQUET_LZ4_TARGET` で指定されたターゲット（設定されている場合）。
2. システムにインストールされたコピー（`CARQUET_BUNDLE_DEPS` が ON の場合はスキップ）。
3. carquet 自身が取得してビルドするコピー。プロジェクトが既にそのライブラリを通常のターゲット名（`libzstd_static`、`zlibstatic`、`lz4_static`）でビルドしている場合、carquet は 2 つ目のコピーを取得する代わりにそのターゲットを再利用します。

したがって、これらのライブラリをベンダリングしておりシステムにコピーがないプロジェクトでは、設定は不要です。明示的に選択したい場合、たとえばシステムの zlib ではなく自前の zlib を、あるいは静的ビルドではなく共有ビルドを使いたい場合は、ターゲットを名指しします：

```cmake
set(CARQUET_ZLIB_TARGET zlibstatic)
add_subdirectory(carquet)
```

何も編集せずに同じことができます：`cmake -DCARQUET_ZLIB_TARGET=zlib ...`。3 つの変数は互いに独立です。名指しされたターゲットは検索も取得も行わずそのままリンクされ、ターゲットではない名前は暗黙のフォールバックではなく configure エラーになります。ターゲットはインポート済みでもツリー内ビルドでも、静的でも共有でもかまいませんが、インクルードディレクトリを usage requirement として持っている必要があります。

</details>

<details>
<summary><b>開発ビルド</b>（テスト、サンプル、ベンチマーク）</summary>

<br/>

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCARQUET_BUILD_DEV=ON
cmake --build build -j$(nproc)
cd build && ctest --output-on-failure
```

[CONTRIBUTING.md](CONTRIBUTING.md) に、プルリクエストごとのゲート（テスト、サニタイザ、ファジング、相互運用テスト、ベンチマーク）が記載されています。

</details>

<details>
<summary><b>API ドキュメント</b>（Doxygen）</summary>

<br/>

公開ヘッダーには [Doxygen](https://www.doxygen.nl) コメントが付いています。HTML リファレンスを生成するには `doxygen` のインストールが必要です（インクルード/依存関係図のために `graphviz` も任意で）。通常のビルドにはどちらも不要です。

```bash
# CMake — configure 時にターゲットを有効にしてからビルドする
cmake -B build -DCARQUET_BUILD_DOCS=ON
cmake --build build --target docs

# xmake — 独立したタスク、configure フラグは不要
xmake docs
```

どちらも `build/docs/html/index.html` に出力します。Doxygen がインストールされていない場合、ターゲット/タスクはその旨を報告するだけで何もしません。

</details>

<details>
<summary><b>AI コーディングエージェント</b>（Claude Code、Codex：すぐに使えるスキル）</summary>

<br/>

[`.agents/skills/carquet/`](.agents/skills/carquet/) は、コーディングエージェントに carquet をプロジェクトへ統合する方法を教える [Agent Skill](https://agentskills.io) です：ビルドの配線（CMake、FetchContent、pkg-config、xmake）、何にどの API を使うか、完全なレシピ、レベル / null ビットマップ / ライフタイムのルール、型マッピング、よくある落とし穴とその対処法。自分のリポジトリにコピーまたはリンクすれば、エージェントが自動的に認識します：

```bash
# Claude Code
mkdir -p .claude/skills && cp -r path/to/carquet/.agents/skills/carquet .claude/skills/
# Codex CLI
mkdir -p .agents/skills && cp -r path/to/carquet/.agents/skills/carquet .agents/skills/
# ベンダリングした carquet：シンボリックリンクで同期を保つ
ln -s ../../third_party/carquet/.agents/skills/carquet .claude/skills/carquet
```

同じファイルは個人用スキル（`~/.claude/skills/`、`~/.codex/skills/`）としても機能します。

</details>

<br/>

### xmake でのビルド

CMake の代替として [xmake](https://xmake.io) ビルド（`xmake.lua`）が用意されており、オプションとデフォルトは同じです。zstd/zlib/lz4 は静的にリンクされるため、バイナリは自己完結しています。

```bash
xmake                    # 静的ライブラリ + `carquet` CLI をビルド（release）
xmake f --dev=y && xmake # テスト、サンプル、ベンチマーク、相互運用テストを追加
xmake test               # テストスイートを実行
```

<details>
<summary><b>configure オプション</b>（<code>xmake f --option=y|n</code>）</summary>

<br/>

オプションは CMake のものと対応しており、`CARQUET_` プレフィックスを除いた小文字です：

<div align="center">

| xmake オプション | デフォルト | 説明 |
|--------------|---------|-------------|
| `--dev` | n | テスト、サンプル、ベンチマーク、相互運用テストをビルド |
| `--tests` / `--examples` / `--benchmarks` / `--interop` | n | グループを個別にビルド |
| `--cli` | y | `carquet` CLI ツールをビルド |
| `--shared` | n | 静的ライブラリの代わりに共有ライブラリをビルド |
| `--openmp` | y | OpenMP によるカラムの並列読み込み（利用できない場合は自動的に無効化） |
| `--native_arch` | n | 最大性能のための `-march=native`（ホスト専用バイナリ） |
| `--sse` / `--avx` / `--avx2` / `--avx512` / `--neon` | y | SIMD 命令セット（自動検出） |
| `--sve` | n | ARM SVE（実験的） |
| `--brotli` | y | BROTLI コーデック（xmake リポジトリの libbrotli） |
| `--fuzz` | n | ファズターゲットをビルド（`--toolchain=clang` を使用） |

</div>

```bash
xmake f -m release --shared=y            # 共有ライブラリ
xmake f --dev=y --avx512=n && xmake      # 開発ビルド、AVX-512 無効
```

</details>

<br/>

---

<!-- ────────────────────────────  C API  ──────────────────────────── -->

## C API

以下に書き込みと読み込みの例を示します。それ以外はすべて [`docs/` のマニュアル](docs/README.md)で扱っています。

### Parquet ファイルの書き込み

```c
#include <carquet/carquet.h>

int main(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    // スキーマを定義
    carquet_schema_t* schema = carquet_schema_create(&err);
    carquet_schema_add_column(schema, "id",    CARQUET_PHYSICAL_INT64,  NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    // ライターを設定
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;

    // 書き込み
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

    // 最高の読み込み性能のために mmap で開く
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = true;

    carquet_reader_t* r = carquet_reader_open("output.parquet", &opts, &err);
    if (!r) { printf("Error: %s\n", err.message); return 1; }

    printf("Rows: %lld, Columns: %d\n",
           (long long)carquet_reader_num_rows(r),
           carquet_reader_num_columns(r));

    // 効率的な反復のためのバッチリーダー
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
        // ids[0..n-1] を処理 ...
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

[マニュアル](docs/README.md)は、どの API をいつ使うべきかを説明しています。ソースのレイアウトとアーキテクチャは [CONTRIBUTING.md](CONTRIBUTING.md) に記載されています。ブラウズ可能な HTML リファレンスは Doxygen でヘッダーから生成できます。[ビルドオプション](#ビルドオプション)の *API ドキュメント* の折りたたみを参照してください。

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

<!-- ────────────────────────────  CLI TOOL  ──────────────────────────── -->

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

./reader                    # data.parquet を読む（デフォルトとして埋め込み済み）
./reader other.parquet      # 別のファイルで上書き
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

<!-- ────────────────────────────  BENCHMARKS  ──────────────────────────── -->

## ベンチマーク

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="res/img/bench_10m_m3_dark.png">
    <img src="res/img/bench_10m_m3_light.png" alt="Apple M3 における 1000 万行の書き込み・読み込み時間：carquet vs Arrow C++ vs PyArrow、コーデック別" width="900">
  </picture>
</p>

1000 万行（最も代表的なサイズ）での結果です。比率が高いほど Carquet が高速です。

ARM（Apple M3）：Carquet 0.7.2 vs Arrow C++ 24.0.0、2026-10-20 計測。x86（Xeon D-1531）：Carquet 0.4.4 vs Arrow C++ 23.0.1。この計測は 0.7.1 のライター並列化より前のもので、その後再計測されていません。

<div align="center">

| | x86 (Xeon D-1531) | | ARM (Apple M3) | |
|---|---|---|---|---|
| **コーデック** | **書き込み** | **読み込み** | **書き込み** | **読み込み** |
| snappy | **1.55x** | **1.25x** | **4.63x** | **1.66x** |
| zstd | **1.31x** | **1.04x** | **5.66x** | **1.36x** |
| lz4 | **1.02x** | 0.83x | **4.47x** | **1.76x** |
| none | **1.13x** | **40.6x**\* | **1.83x** | **38.9x**\* |

</div>

\* 非圧縮読み込みは mmap ゼロコピーを使用します。方法論の注記を参照してください。

圧縮読み込みはすべての値の完全な展開とデコードを伴い、近道はありません。また、両ライブラリは同じシステムの lz4/zstd 共有ライブラリを使用しているため、生のコーデック速度は同一です。最も意味のある比較は（完全な結果の中の）**同一ファイルのクロスリード**表で、両ライブラリがまったく同じ Parquet ファイルを読みます。この同条件のテストで、Carquet は圧縮データを Arrow C++ より **1.3〜2.1 倍**速く読みます（x86 の計測では 1.6〜2.6 倍）。

<details>
<summary><b>ベンチマークの方法論</b></summary>

<br/>

すべてのベンチマークは、同一のデータ（決定論的 LCG PRNG）、同一の Parquet 設定（辞書なし、浮動小数点に BYTE_STREAM_SPLIT、ページチェックサム、mmap 読み込み）を使用し、11〜51 回の反復のトリム中央値を取り、書き込みフェーズと読み込みフェーズの間に OS ページキャッシュをパージし、構成間にクールダウンを挟んでいます。スキーマ：3 カラム（INT64、DOUBLE、INT32）。比較対象は Arrow C++ の低レベル Parquet リーダー（Arrow Table のマテリアライズを回避）と PyArrow です。

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

*Carquet 0.7.2、MacBook Air M3、16GB RAM、macOS 26.5、Arrow C++ 24.0.0、PyArrow 23.0.1、ZSTD レベル 1。1000 万行の表は 2026-10-20 に計測。100 万行と 10 万行の表は 2026-10-08 の 0.7.1 の計測です（再計測していません）。*

#### 1000 万行 vs Arrow C++

<div align="center">

| コーデック | Carquet 書き込み | Arrow C++ 書き込み | 書き込み比 | Carquet 読み込み | Arrow C++ 読み込み | 読み込み比 | サイズ |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **74.7ms** | 137ms | **1.83x** | **0.36ms** | 14.0ms | **38.9x**\* | 200.0MB |
| snappy | **52.7ms** | 244ms | **4.63x** | **12.4ms** | 20.6ms | **1.66x** | 131.2MB |
| zstd | **60.4ms** | 342ms | **5.66x** | **19.8ms** | 26.9ms | **1.36x** | 99.9MB |
| lz4 | **57.2ms** | 256ms | **4.47x** | **9.22ms** | 16.2ms | **1.76x** | 128.9MB |

</div>

<br/>

#### 100 万行 vs Arrow C++

<div align="center">

| コーデック | Carquet 書き込み | Arrow C++ 書き込み | 書き込み比 | Carquet 読み込み | Arrow C++ 読み込み | 読み込み比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **6.48ms** | 14.82ms | **2.29x** | **0.06ms** | 1.63ms | **27.2x**\* |
| snappy | **14.28ms** | 26.44ms | **1.85x** | **2.14ms** | 2.85ms | **1.33x** |
| zstd | **14.76ms** | 36.17ms | **2.45x** | **2.70ms** | 3.46ms | **1.28x** |
| lz4 | **12.30ms** | 26.10ms | **2.12x** | **1.04ms** | 1.87ms | **1.80x** |

</div>

<br/>

#### 10 万行 vs Arrow C++

<div align="center">

| コーデック | Carquet 書き込み | Arrow C++ 書き込み | 書き込み比 | Carquet 読み込み | Arrow C++ 読み込み | 読み込み比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **1.07ms** | 1.79ms | **1.67x** | **0.02ms** | 0.25ms | **12.5x**\* |
| snappy | **1.60ms** | 2.63ms | **1.64x** | **0.34ms** | 0.91ms | **2.68x** |
| zstd | **1.73ms** | 3.63ms | **2.10x** | **0.63ms** | 1.26ms | **2.00x** |
| lz4 | **1.52ms** | 2.62ms | **1.72x** | **0.25ms** | 0.55ms | **2.20x** |

</div>

<br/>

#### 同一ファイルのクロスリード（1000 万行）

両ライブラリが**同じ** Parquet ファイルを読みます。最も公平な同条件比較です。

<div align="center">

| コーデック | ライター | Carquet 読み込み | Arrow C++ 読み込み | 比率 |
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

#### 1000 万行 vs PyArrow

<div align="center">

| コーデック | Carquet 書き込み | PyArrow 書き込み | 書き込み比 | Carquet 読み込み | PyArrow 読み込み | 読み込み比 |
|-------|--------------|---------------|---------|-------------|--------------|---------|
| none | **74.7ms** | 202ms | **2.71x** | **0.36ms** | 36.2ms | **100.5x**\* |
| snappy | **52.7ms** | 312ms | **5.93x** | **12.4ms** | 46.6ms | **3.75x** |
| zstd | **60.4ms** | 410ms | **6.78x** | **19.8ms** | 75.6ms | **3.82x** |
| lz4 | **57.2ms** | 321ms | **5.62x** | **9.22ms** | 40.0ms | **4.33x** |

</div>

<br/>

#### 100 万行 vs PyArrow

<div align="center">

| コーデック | Carquet 書き込み | PyArrow 書き込み | 書き込み比 | Carquet 読み込み | PyArrow 読み込み | 読み込み比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **6.48ms** | 20.17ms | **3.11x** | **0.06ms** | 2.79ms | **46.5x**\* |
| snappy | **14.28ms** | 31.54ms | **2.21x** | **2.14ms** | 4.01ms | **1.87x** |
| zstd | **14.76ms** | 42.51ms | **2.88x** | **2.70ms** | 4.42ms | **1.64x** |
| lz4 | **12.30ms** | 32.03ms | **2.60x** | **1.04ms** | 3.22ms | **3.10x** |

</div>

<br/>

#### 10 万行 vs PyArrow

<div align="center">

| コーデック | Carquet 書き込み | PyArrow 書き込み | 書き込み比 | Carquet 読み込み | PyArrow 読み込み | 読み込み比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **1.07ms** | 2.12ms | **1.98x** | **0.02ms** | 0.24ms | **12.0x**\* |
| snappy | **1.60ms** | 3.17ms | **1.98x** | **0.34ms** | 0.60ms | **1.76x** |
| zstd | **1.73ms** | 4.16ms | **2.40x** | **0.63ms** | 0.81ms | **1.29x** |
| lz4 | **1.52ms** | 3.21ms | **2.11x** | **0.25ms** | 0.42ms | **1.68x** |

</div>

\* ゼロコピー mmap パス

</details>

<details>
<summary><b>ベンチマークを自分で実行する</b></summary>

<br/>

```bash
# 最大限の最適化でビルド
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON -DCARQUET_BUILD_DEV=ON
cmake --build build -j$(nproc)

cd build
./benchmark_carquet                     # Carquet 単体
python3 ../benchmark/run_benchmark.py   # 完全な比較（+ PyArrow、+ Arrow C++）

# 1 億行（xlarge）の構成をスキップ — コーデックごとに約 2GB のファイルを書き、
# ハードウェアによっては 30 分以上かかる
python3 ../benchmark/run_benchmark.py --skip-xlarge

# ZSTD レベルを上書き（デフォルト：1）
CARQUET_BENCH_ZSTD_LEVEL=3 python3 ../benchmark/run_benchmark.py
```

Arrow C++ との比較を含めるには、configure 時にそのターゲットを有効にします：

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON \
  -DCARQUET_BUILD_BENCHMARKS=ON \
  -DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON
cmake --build build -j$(nproc)

# またはカスタムの Arrow インストールを指定
cmake -B build ... -DCARQUET_ARROW_CPP_ROOT=/path/to/arrow-prefix
```

Arrow C++ ベンチマークは、低レベルの `parquet::ParquetFileReader` API（Arrow Table マテリアライズのオーバーヘッドを回避）を並列行グループリーダーとともに使用します。**同一ファイルのクロスリード**モードでは両ライブラリがまったく同じ Parquet ファイルを読み、ページサイズ、エンコーディング、行グループレイアウトの違いを排除します。両ベンチマークは、同一のデータ、行グループサイズ、辞書なし、ページチェックサム、mmap 読み込み、浮動小数点に BYTE_STREAM_SPLIT を使用します。

</details>

<br/>

---

<!-- ────────────────────────────  FEATURE SUPPORT  ──────────────────────────── -->

## Parquet 機能サポート

ライブラリが現在実装している内容を、領域別に示します。

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

<!-- ────────────────────────────  INTEROPERABILITY  ──────────────────────────── -->

## 相互運用性

Carquet のファイルは、PyArrow、DuckDB、Spark、およびあらゆる Parquet リーダーと完全に互換です。

```python
import pyarrow.parquet as pq
table = pq.read_table("carquet_output.parquet")  # そのまま動きます
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

<!-- ────────────────────────────  CONTRIBUTING & LICENSE  ──────────────────────────── -->

## コントリビューション、引用、ライセンス

- **コントリビューション**：[CONTRIBUTING.md](CONTRIBUTING.md) に、セットアップ、コードスタイル、プルリクエストごとのゲート（テスト、サニタイザ、ファジング、相互運用テスト、ベンチマーク）が記載されています。リリースノートは [CHANGELOG.md](CHANGELOG.md) にあります。
- **引用**：[`CITATION.cff`](CITATION.cff) を同梱しています。GitHub の「Cite this repository」ボタンで BibTeX または APA 形式に整形できます。
- **ライセンス**：MIT、[LICENSE](LICENSE) を参照してください。
