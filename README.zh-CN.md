<h1 align="center">
  <img src="res/img/carquet_banner.png" alt="Carquet" width="840" />
</h1>

<p align="center">
  一个快速、纯 C 实现的 Apache Parquet 文件读写库。
</p>

<p align="center">
  <a href="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml"><img src="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml/badge.svg" alt="构建" /></a>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue" alt="平台" />
  <img src="https://img.shields.io/badge/C-C11-blue" alt="C 标准" />
  <a href="https://github.com/Vitruves/carquet/releases"><img src="https://img.shields.io/github/v/release/Vitruves/carquet?color=blue" alt="发布版本" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="许可证" /></a>
  <a href="https://github.com/Vitruves/carquet/stargazers"><img src="https://img.shields.io/github/stars/Vitruves/carquet?style=flat&color=blue" alt="Star 数" /></a>
  <a href="https://github.com/Vitruves/carquet/graphs/traffic"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/Vitruves/carquet/traffic/clones.json" alt="克隆数" /></a>
  <br/>
  <img src="https://img.shields.io/badge/SIMD-SSE4.2%20%C2%B7%20AVX2%20%C2%B7%20AVX--512%20%C2%B7%20NEON%20%C2%B7%20SVE-555" alt="SIMD" />
  <img src="https://img.shields.io/badge/deps-zstd%20%C2%B7%20zlib%20%C2%B7%20lz4%20%C2%B7%20brotli-555" alt="依赖" />
</p>

<p align="center">
  <a href="#亮点">亮点</a> &nbsp;·&nbsp;
  <a href="#安装">安装</a> &nbsp;·&nbsp;
  <a href="#c-api">C API</a> &nbsp;·&nbsp;
  <a href="#cli-工具">CLI 工具</a> &nbsp;·&nbsp;
  <a href="#基准测试">基准测试</a> &nbsp;·&nbsp;
  <a href="#parquet-功能支持">功能支持</a> &nbsp;·&nbsp;
  <a href="#互操作性">互操作性</a> &nbsp;·&nbsp;
  <a href="#ready-for-ai-coding-agents">AI 编程助手</a> &nbsp;·&nbsp;
  <a href="docs/README.md">手册</a> &nbsp;·&nbsp;
  <a href="CHANGELOG.md">更新日志</a>
</p>

<p align="center">
  <a href="README.md">English</a> &nbsp;·&nbsp;
  <b>简体中文</b> &nbsp;·&nbsp;
  <a href="README.ja.md">日本語</a> &nbsp;·&nbsp;
  <a href="README.ko.md">한국어</a> &nbsp;·&nbsp;
  <a href="README.de.md">Deutsch</a> &nbsp;·&nbsp;
  <a href="README.fr.md">Français</a> &nbsp;·&nbsp;
  <a href="README.pt-BR.md">Português (Brasil)</a>
</p>

<br/>

---

<!-- ────────────────────────────  HIGHLIGHTS  ──────────────────────────── -->

## 亮点

**小巧、自包含**

- **纯 C11**，仅四个外部依赖（zstd、zlib、lz4、brotli），缺失时由 CMake 自动获取；Brotli 可以关闭
- **约 500 KB 的库**（Release 构建去除符号后，编解码器单独链接），而 Arrow 约 50 MB 以上
- **内置 CLI**，用于文件检视（`schema`、`info`、`head`、`tail`、`stat` 等）和 C 代码生成（`codegen`）

**快速**

- 未压缩数据的读取比 Arrow C++ **快 39 倍**（mmap 零拷贝），比 PyArrow **快 100 倍**
- 对同一文件的压缩读取比 Arrow C++ **快 1.3–2.1 倍**（交叉读取基准测试，Apple M3；x86 上为 1.6–2.6 倍）
- 在 Apple M3 上**写入比 Arrow C++ 快 1.8–5.7 倍**（在较早的 x86 测试中为 1.0–1.6 倍）
- 读取 1000 万行未压缩数据仅需 **0.36 ms**（Apple M3 上的 mmap 零拷贝）
- SIMD 优化（SSE4.2、AVX2、AVX-512、NEON、SVE），运行时检测并提供标量回退实现

**完整且兼容**

- 广泛的 Parquet 覆盖：所有物理类型和编码，Snappy / GZIP / LZ4 / ZSTD / Brotli，嵌套模式，布隆过滤器，页索引（[完整矩阵](#parquet-功能支持)）
- 开箱即用地兼容 PyArrow、DuckDB、Spark
- [Arrow C Data Interface](docs/reading.md#export-to-arrow-c-data-interface) 桥接（`carquet_arrow_export_*` / `carquet_arrow_import_*` / `carquet_writer_write_arrow` / `carquet_reader_read_arrow`），与 Arrow 生态系统进行零依赖、少拷贝的数据交换，支持**任意深度的嵌套读写**（任意深度的 struct / list / map）

<a id="ready-for-ai-coding-agents"></a>

**面向 AI 编程助手，开箱即用**

- [`.agents/skills/carquet/`](.agents/skills/carquet/) 是一个 [Agent Skill](https://agentskills.io)（Claude Code、Codex CLI），用于教会编程助手将 carquet 集成到你的项目中：构建接线、何种任务用何种 API、完整的示例配方、层级 / 空值位图 / 生命周期规则、类型映射以及常见陷阱。把它复制到你的仓库中，助手就会自行使用：

  ```bash
  cp -r path/to/carquet/.agents/skills/carquet .claude/skills/   # Claude Code
  cp -r path/to/carquet/.agents/skills/carquet .agents/skills/   # Codex CLI
  ```

<br/>

### 三十秒上手 carquet

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

完整的写入与读取循环（含批处理和可空列）见 [C API](#c-api)。

<br/>

### carquet 的定位

carquet 是小巧、可嵌入的选择：一个约 500 KB 的 C11 库，仅四个自动获取的依赖，一分钟内即可构建完成；相比之下 Arrow C++ 重达 50 MB 以上，PyArrow wheel 约 100 MB。在下文的 1000 万行基准测试中，它写入压缩文件比 Arrow C++ 快 1.8–5.7 倍，读取同一个文件快 1.3–2.1 倍；未压缩数据直接从 mmap 提供，仅需 0.36 ms。当你需要 Arrow 的生态系统（计算内核、数据集、Flight、加密）时，Arrow 仍然是正确的选择；carquet 则面向那些只需要快速、正确地读写 Parquet 的 C 和 C++ 程序、嵌入式系统和语言运行时。

<br/>

---

<!-- ────────────────────────────  INSTALLATION  ──────────────────────────── -->

## 安装

Carquet 可通过 CMake（主要方式）或 [xmake](https://xmake.io) 构建，支持 Linux、macOS 和 Windows。

### 环境要求

- C11 编译器（GCC 4.9+、Clang 3.4+、MSVC 2015+）
- CMake 3.16+（或 [xmake](https://xmake.io)，见[使用 xmake 构建](#使用-xmake-构建)）
- zstd、zlib、lz4、brotli（缺失时自动获取；Brotli 可选，`-DCARQUET_ENABLE_BROTLI=OFF`）
- OpenMP（可选，用于并行列读取和写入）

<br/>

### 快速开始（make）

`make` 封装脚本会驱动一次优化的 CMake 构建并安装到 `/usr/local`：

```bash
git clone https://github.com/Vitruves/carquet.git
cd carquet
make                              # 优化构建（运行 `make help` 列出所有目标）
sudo make install                 # 安装到 /usr/local（用 PREFIX=/opt/carquet 覆盖）
```

<br/>

### 构建与安装（CMake）

当你需要特定的构建选项或自定义安装前缀时，直接调用 CMake：

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release   # 添加选项，例如 -DCARQUET_BUILD_SHARED=ON
cmake --build build -j$(nproc)
sudo cmake --install build --prefix /usr/local
```

两种方式都会安装：

- `libcarquet.a`（使用 `-DCARQUET_BUILD_SHARED=ON` 时为 `.so` / `.dylib`）
- `include/carquet/` 头文件
- `carquet` CLI 可执行文件
- `carquet.pc`（pkg-config）以及供 `find_package(carquet)` 使用的 CMake 包配置

安装完成后，使用 `-lcarquet` 链接，或通过 `pkg-config --cflags --libs carquet` 获取编译标志。

<br/>

### Python（pip）

`carquet` wheel 内置了该库，并通过 Arrow C Data Interface 暴露出来，因此表数据可以在 PyArrow 之间来回流动，而无需重复一份 Parquet 解码逻辑：

```bash
pip install carquet
```

```python
import carquet
table = carquet.read_table("data.parquet")                 # pyarrow.Table
carquet.write_table(table, "out.parquet", compression="zstd")
print(carquet.info("out.parquet").num_rows)
```

该包位于 [`python/`](python/)（scikit-build-core、ctypes，PyArrow 在运行时可选）；API 和 wheel 构建方法见其 [README](python/README.md)。

<br/>

### 构建选项

在 CMake 配置步骤中传入这些选项（`-DOPTION=ON|OFF`）：

<div align="center">

| 选项 | 默认值 | 说明 |
|--------|---------|-------------|
| `CARQUET_BUILD_DEV` | OFF | 构建所有内容（测试、示例、基准测试、互操作） |
| `CARQUET_BUILD_TESTS` | OFF | 仅构建测试套件 |
| `CARQUET_BUILD_CLI` | ON | 构建 `carquet` CLI 工具 |
| `CARQUET_BUILD_SHARED` | OFF | 构建共享库而非静态库 |
| `CARQUET_NATIVE_ARCH` | OFF | `-march=native` 以获得最高性能（仅限本机运行的二进制） |
| `CARQUET_ENABLE_SVE` | OFF | ARM SVE（实验性） |
| `CARQUET_ENABLE_BROTLI` | ON | BROTLI 编解码器（查找或获取 libbrotli；OFF 则去掉该依赖） |

</div>

所有 x86 SIMD（SSE、AVX、AVX2、AVX-512）和 ARM NEON 均自动检测并默认启用。

<details>
<summary><b>全部构建选项</b></summary>

<br/>

<div align="center">

| 选项 | 默认值 | 说明 |
|--------|---------|-------------|
| `CARQUET_BUILD_EXAMPLES` | OFF | 构建示例程序 |
| `CARQUET_BUILD_BENCHMARKS` | OFF | 构建基准测试和性能分析程序 |
| `CARQUET_BUILD_ARROW_CPP_BENCHMARK` | OFF | 可选的 Arrow C++ 对比基准测试 |
| `CARQUET_BUILD_INTEROP` | OFF | 构建互操作性测试 |
| `CARQUET_BUILD_FUZZ` | OFF | 构建模糊测试目标 |
| `CARQUET_BUILD_DOCS` | OFF | 添加 `docs` 目标（需要 Doxygen） |
| `CARQUET_BUNDLE_DEPS` | OFF | 始终获取并静态打包 zstd/zlib/lz4/brotli（可移植的发布产物） |
| `CARQUET_ZSTD_TARGET` / `CARQUET_ZLIB_TARGET` / `CARQUET_LZ4_TARGET` / `CARQUET_BROTLI_TARGET` | *（空）* | 为该库链接已有的 CMake 目标，而不是查找或获取它 |
| `CARQUET_ARROW_CPP_ROOT` | *（空）* | 用于对比基准测试的自定义 Arrow C++ 安装路径 |
| `CARQUET_ENABLE_SSE` | ON | SSE 优化（x86，自动检测） |
| `CARQUET_ENABLE_AVX` | ON | AVX 优化（x86，自动检测） |
| `CARQUET_ENABLE_AVX2` | ON | AVX2 优化（x86，自动检测） |
| `CARQUET_ENABLE_AVX512` | ON | AVX-512 优化（x86，自动检测） |
| `CARQUET_ENABLE_NEON` | ON | NEON 优化（ARM，自动检测） |

</div>

</details>

<details>
<summary><b>使用你自己的 zstd / zlib / lz4</b>（将 carquet 嵌入更大的 CMake 项目）</summary>

<br/>

当 carquet 通过 `add_subdirectory()` 或 `FetchContent` 引入时，它按以下顺序解析每个压缩库：

1. 由 `CARQUET_ZSTD_TARGET`、`CARQUET_ZLIB_TARGET` 或 `CARQUET_LZ4_TARGET` 指定的目标（若已设置）。
2. 系统已安装的副本（`CARQUET_BUNDLE_DEPS` 为 ON 时跳过）。
3. 它自行获取并构建的副本。如果你的项目已经以常用目标名（`libzstd_static`、`zlibstatic`、`lz4_static`）构建了该库，carquet 会复用该目标，而不是再获取一份。

因此，自带这些库且系统中没有副本的项目无需任何配置。若要显式选择，例如用你自己的 zlib 而非系统的，或用你的共享构建而非静态构建，请指定目标名：

```cmake
set(CARQUET_ZLIB_TARGET zlibstatic)
add_subdirectory(carquet)
```

不修改任何文件也能做到同样的事：`cmake -DCARQUET_ZLIB_TARGET=zlib ...`。这三个变量相互独立。指定的目标会原样链接，不做查找也不获取；如果给出的名字不是一个目标，则会产生配置错误，而不是静默回退。该目标可以是导入的或在源码树内构建的、静态的或共享的，并且必须以使用要求（usage requirement）的形式携带其头文件目录。

</details>

<details>
<summary><b>开发构建</b>（测试、示例、基准测试）</summary>

<br/>

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCARQUET_BUILD_DEV=ON
cmake --build build -j$(nproc)
cd build && ctest --output-on-failure
```

[CONTRIBUTING.md](CONTRIBUTING.md) 描述了每个 Pull Request 必须通过的检查关卡（测试、Sanitizer、模糊测试、互操作、基准测试）。

</details>

<details>
<summary><b>API 文档</b>（Doxygen）</summary>

<br/>

公共头文件带有 [Doxygen](https://www.doxygen.nl) 注释。生成 HTML 参考文档需要安装 `doxygen`（可选安装 `graphviz` 以生成包含/依赖关系图）；常规构建不需要二者。

```bash
# CMake —— 在配置时启用该目标，然后构建它
cmake -B build -DCARQUET_BUILD_DOCS=ON
cmake --build build --target docs

# xmake —— 独立任务，无需配置标志
xmake docs
```

两者都输出到 `build/docs/html/index.html`。如果未安装 Doxygen，该目标/任务只会报告这一情况，不做任何事。

</details>

<details>
<summary><b>AI 编程助手</b>（Claude Code、Codex：现成的技能）</summary>

<br/>

[`.agents/skills/carquet/`](.agents/skills/carquet/) 是一个 [Agent Skill](https://agentskills.io)，用于教会编程助手如何将 carquet 集成到你的项目中：构建接线（CMake、FetchContent、pkg-config、xmake）、何种任务用何种 API、完整的示例配方、层级 / 空值位图 / 生命周期规则、类型映射以及常见陷阱及其修复方法。把它复制或链接到你自己的仓库中，助手会自动识别并使用：

```bash
# Claude Code
mkdir -p .claude/skills && cp -r path/to/carquet/.agents/skills/carquet .claude/skills/
# Codex CLI
mkdir -p .agents/skills && cp -r path/to/carquet/.agents/skills/carquet .agents/skills/
# 以 vendored 方式内置 carquet 时：用符号链接保持同步
ln -s ../../third_party/carquet/.agents/skills/carquet .claude/skills/carquet
```

这些文件同样可以作为个人技能使用（`~/.claude/skills/`、`~/.codex/skills/`）。

</details>

<br/>

### 使用 xmake 构建

提供了 [xmake](https://xmake.io) 构建（`xmake.lua`）作为 CMake 的替代方案，选项和默认值相同；zstd/zlib/lz4 静态链接，因此生成的二进制文件是自包含的。

```bash
xmake                    # 构建静态库 + `carquet` CLI（release）
xmake f --dev=y && xmake # 加上测试、示例、基准测试、互操作
xmake test               # 运行测试套件
```

<details>
<summary><b>配置选项</b>（<code>xmake f --option=y|n</code>）</summary>

<br/>

选项与 CMake 的一一对应，去掉 `CARQUET_` 前缀并改为小写：

<div align="center">

| xmake 选项 | 默认值 | 说明 |
|--------------|---------|-------------|
| `--dev` | n | 构建测试、示例、基准测试和互操作 |
| `--tests` / `--examples` / `--benchmarks` / `--interop` | n | 单独构建某一组 |
| `--cli` | y | 构建 `carquet` CLI 工具 |
| `--shared` | n | 构建共享库而非静态库 |
| `--openmp` | y | OpenMP 并行列读取（不可用时自动禁用） |
| `--native_arch` | n | `-march=native` 以获得最高性能（仅限本机运行的二进制） |
| `--sse` / `--avx` / `--avx2` / `--avx512` / `--neon` | y | SIMD 指令集（自动检测） |
| `--sve` | n | ARM SVE（实验性） |
| `--brotli` | y | BROTLI 编解码器（来自 xmake 仓库的 libbrotli） |
| `--fuzz` | n | 构建模糊测试目标（使用 `--toolchain=clang`） |

</div>

```bash
xmake f -m release --shared=y            # 共享库
xmake f --dev=y --avx512=n && xmake      # 开发构建，禁用 AVX-512
```

</details>

<br/>

---

<!-- ────────────────────────────  C API  ──────────────────────────── -->

## C API

下面是一个写入示例和一个读取示例；其余内容见 [`docs/` 中的手册](docs/README.md)。

### 写入 Parquet 文件

```c
#include <carquet/carquet.h>

int main(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    // 定义模式
    carquet_schema_t* schema = carquet_schema_create(&err);
    carquet_schema_add_column(schema, "id",    CARQUET_PHYSICAL_INT64,  NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    // 配置写入器
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;

    // 写入
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

### 读取 Parquet 文件

```c
#include <carquet/carquet.h>
#include <stdio.h>

int main(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    // 使用 mmap 打开以获得最佳读取性能
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = true;

    carquet_reader_t* r = carquet_reader_open("output.parquet", &opts, &err);
    if (!r) { printf("Error: %s\n", err.message); return 1; }

    printf("Rows: %lld, Columns: %d\n",
           (long long)carquet_reader_num_rows(r),
           carquet_reader_num_columns(r));

    // 使用批量读取器高效迭代
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
        // 处理 ids[0..n-1] ...
        carquet_row_batch_free(batch);
        batch = NULL;
    }

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    return 0;
}
```

<br/>

### 更多示例

平面读写之外的所有内容都在手册中；每一页都链接到一个可运行的示例：

<div align="center">

| 你想要… | 参见 |
|---|---|
| 可空列、行组、缓冲区输出 | [`docs/writing.md`](docs/writing.md) |
| 列表、映射、组、定义/重复层级 | [`docs/nested-data.md`](docs/nested-data.md) |
| 列投影、统计信息、元数据检视 | [`docs/reading.md`](docs/reading.md) |
| 谓词下推、页级过滤 | [`docs/reading.md`](docs/reading.md) |
| 向现有文件追加行组 | [`docs/writing.md`](docs/writing.md#append-to-an-existing-file) |
| 压缩、自定义编解码器、写入器调优 | [`docs/writing.md`](docs/writing.md)、[`docs/performance.md`](docs/performance.md) |
| mmap、零拷贝、预缓冲、I/O 合并 | [`docs/performance.md`](docs/performance.md) |
| 错误码与恢复提示 | [`docs/error-handling.md`](docs/error-handling.md) |

</div>

<br/>

### API 参考

完整 API 位于 [`include/carquet/carquet.h`](include/carquet/carquet.h)。主要类型：

<div align="center">

| 类型 | 用途 |
|------|---------|
| `carquet_reader_t` | 文件读取器（可从路径、`FILE*` 或内存缓冲区打开） |
| `carquet_writer_t` | 文件写入器 |
| `carquet_batch_reader_t` | 高层批量迭代 |
| `carquet_schema_t` | 模式定义与自省 |
| `carquet_error_t` | 丰富的错误信息（错误码、消息、源码位置、恢复提示） |

</div>

[手册](docs/README.md)说明了何时使用哪一层接口；源码布局和架构记录在 [CONTRIBUTING.md](CONTRIBUTING.md) 中。可以用 Doxygen 从头文件生成可浏览的 HTML 参考文档，见[构建选项](#构建选项)下的 *API 文档* 折叠区。

<details>
<summary><b>示例应用：MetricStore</b>（<code>mocklib/</code>）</summary>

<br/>

[`mocklib/`](mocklib/) 是 **MetricStore**，一个构建在 carquet 之上的完整、自包含的示例应用，既可作为真实使用场景的参考，也可作为公共 API 的端到端集成测试。

它模拟了一个时间序列遥测存储（摄入事件，然后进行检视和查询），在此过程中调用了 **carquet 的 143 个公共函数中的 117 个（约 82 %）**：

- **写入侧**：带逻辑类型的模式构建、每列的编码 / 压缩 / 布隆过滤器调优、页索引和统计信息
- **读取侧**：列投影、谓词下推、页级过滤、布隆过滤器成员查询、Arrow C Data Interface 桥接、嵌套 `LIST` 重建、缓冲区/追加 I/O 以及元数据自省

它像任何下游项目一样链接 carquet（`find_package(carquet)` / `add_subdirectory`），附带一个自检的往返测试，其 `write_sample` 可执行文件会生成一个文件供外部检视。见 [`mocklib/README.md`](mocklib/README.md)。

</details>

<br/>

---

<!-- ────────────────────────────  CLI TOOL  ──────────────────────────── -->

## CLI 工具

Carquet 附带一个命令行工具，用于检视 Parquet 文件和生成 C 读取器代码。它默认与库一起构建和安装。

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

典型用法：

```bash
carquet schema data.parquet
carquet head -n 20 data.parquet
carquet stat data.parquet
carquet validate data.parquet
```

<br/>

### 行过滤

`cat`、`count`、`head`、`tail`、`sample` 和 `export` 接受 `-p / --filter EXPR`，将行谓词下推到页级别：只有列索引 min/max 可能匹配该谓词的页才会被解压。

```bash
carquet cat -p "price > 100 AND status = 'active'" data.parquet
carquet count --filter "id >= 1000" data.parquet
carquet export --filter "ts IS NOT NULL" -c id,ts data.parquet
```

语法为 `column OP value [AND column OP value]...`，其中 `OP` ∈ {`=`, `==`, `!=`, `<>`, `<`, `<=`, `>`, `>=`}，另有 `column IS NULL` / `column IS NOT NULL`。

> 过滤要求文件带有页索引（`write_page_index = true`）。

<br/>

### 代码生成

从任意 Parquet 文件的模式生成一个完整、可编译的 C 读取器：

```bash
carquet codegen -f data.parquet -o reader.c
# Generated: reader.c
# Compile:   clang -o reader reader.c -I.../include -L.../build -lcarquet ...

./reader                    # 读取 data.parquet（作为默认值内嵌）
./reader other.parquet      # 改用另一个文件
```

<details>
<summary><b>codegen 选项</b></summary>

<br/>

<div align="center">

| 标志 | 说明 |
|------|-------------|
| `-f`, `--file FILE` | 用于读取模式的 Parquet 文件 |
| `-o`, `--output FILE` | 输出源文件（默认：stdout） |
| `--mmap` | 在生成的代码中使用内存映射 I/O |
| `--skeleton` | 生成空的 `process_batch` 以便填入自定义逻辑 |
| `-c`, `--columns COLS` | 逗号分隔的列过滤器 |
| `-b`, `--batch-size N` | 批大小（默认：1024） |

</div>

</details>

<br/>

---

<!-- ────────────────────────────  BENCHMARKS  ──────────────────────────── -->

## 基准测试

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="res/img/bench_10m_m3_dark.png">
    <img src="res/img/bench_10m_m3_light.png" alt="Apple M3 上 1000 万行的写入和读取耗时：carquet、Arrow C++ 与 PyArrow 按编解码器对比" width="900">
  </picture>
</p>

在 1000 万行（最具代表性的规模）下；比值越高 = Carquet 越快。

ARM（Apple M3）：Carquet 0.7.2 对比 Arrow C++ 24.0.0，测量于 2026-10-20。x86（Xeon D-1531）：Carquet 0.4.4 对比 Arrow C++ 23.0.1；该次测试早于 0.7.1 引入的写入器并行化，此后未重新测量。

<div align="center">

| | x86 (Xeon D-1531) | | ARM (Apple M3) | |
|---|---|---|---|---|
| **编解码器** | **写入** | **读取** | **写入** | **读取** |
| snappy | **1.55x** | **1.25x** | **4.63x** | **1.66x** |
| zstd | **1.31x** | **1.04x** | **5.66x** | **1.36x** |
| lz4 | **1.02x** | 0.83x | **4.47x** | **1.76x** |
| none | **1.13x** | **40.6x**\* | **1.83x** | **38.9x**\* |

</div>

\* 未压缩读取使用 mmap 零拷贝，见方法说明。

压缩读取涉及对每个值的完整解压和解码，没有任何捷径，且两个库使用相同的系统 lz4/zstd 共享库，因此原始编解码器速度完全相同。最有意义的对比是**同文件交叉读取**表（见完整结果），其中两个库读取完全相同的 Parquet 文件：在这一同等条件的测试中，Carquet 读取压缩数据比 Arrow C++ **快 1.3–2.1 倍**（x86 测试：1.6–2.6 倍）。

<details>
<summary><b>基准测试方法</b></summary>

<br/>

所有基准测试使用相同的数据（确定性 LCG 伪随机数生成器）、相同的 Parquet 设置（无字典、浮点数使用 BYTE_STREAM_SPLIT、页校验和、mmap 读取），取 11–51 次迭代的截尾中位数，在写入和读取阶段之间清空操作系统页缓存，并在各配置之间留有冷却时间。模式：3 列（INT64、DOUBLE、INT32）。对比对象为 Arrow C++ 底层 Parquet 读取器（绕过 Arrow Table 物化）和 PyArrow。

**同文件交叉读取**基准测试是最公平的对比：两个库读取完全相同的 Parquet 文件（由其中一个写入，两者都读取）。这消除了页大小、编码选择和行组布局上的差异。

标有 \* 的**未压缩读取**使用 Carquet 的 **mmap 零拷贝路径**：对于 PLAIN 编码、未压缩、定长、required 的列，批量读取器直接返回指向内存映射文件的指针，没有 memcpy。Arrow 总是物化到自己的缓冲区中。**压缩读取数据是端到端读取吞吐量最具代表性的衡量指标。**

</details>

<details>
<summary><b>完整 x86 结果</b>（Intel Xeon D-1531，Linux）</summary>

<br/>

*12 线程 @ 2.2GHz，32GB 内存，Ubuntu 24.04，ZSTD 级别 1*

#### 1000 万行 vs Arrow C++

<div align="center">

| 编解码器 | Carquet 写入 | Arrow C++ 写入 | 写入比 | Carquet 读取 | Arrow C++ 读取 | 读取比 | 大小 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **1557ms** | 1766ms | **1.13x** | **1.25ms** | 50.8ms | **40.6x**\* | 190.7MB |
| snappy | **1002ms** | 1549ms | **1.55x** | **78ms** | 97.8ms | **1.25x** | 125.1MB |
| zstd | **1311ms** | 1714ms | **1.31x** | **76.8ms** | 80.2ms | **1.04x** | 95.3MB |
| lz4 | **1521ms** | 1554ms | **1.02x** | 59.1ms | **49.0ms** | 0.83x | 122.9MB |

</div>

<br/>

#### 100 万行 vs Arrow C++

<div align="center">

| 编解码器 | Carquet 写入 | Arrow C++ 写入 | 写入比 | Carquet 读取 | Arrow C++ 读取 | 读取比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **180ms** | 196ms | **1.09x** | **0.22ms** | 6.2ms | **28x**\* |
| snappy | **141ms** | 148ms | **1.05x** | **8.1ms** | 11.6ms | **1.44x** |
| zstd | **131ms** | 185ms | **1.41x** | 10.3ms | **9.1ms** | 0.88x |
| lz4 | **143ms** | 149ms | **1.04x** | 8.5ms | **6.1ms** | 0.72x |

</div>

<br/>

#### 10 万行 vs Arrow C++

<div align="center">

| 编解码器 | Carquet 写入 | Arrow C++ 写入 | 写入比 | Carquet 读取 | Arrow C++ 读取 | 读取比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **14.1ms** | 18.4ms | **1.30x** | **0.11ms** | 2.18ms | **19.8x**\* |
| snappy | **10.1ms** | 10.6ms | **1.05x** | **1.27ms** | 5.97ms | **4.70x** |
| zstd | **8.7ms** | 14.1ms | **1.62x** | **1.58ms** | 3.88ms | **2.46x** |
| lz4 | **9.6ms** | 11.0ms | **1.14x** | **0.77ms** | 2.78ms | **3.61x** |

</div>

<br/>

#### 同文件交叉读取（1000 万行）

两个库读取**同一个** Parquet 文件，这是最公平的同等条件对比。

<div align="center">

| 编解码器 | 写入方 | Carquet 读取 | Arrow C++ 读取 | 比值 |
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

| 编解码器 | Carquet 写入 | PyArrow 写入 | 写入比 | Carquet 读取 | PyArrow 读取 | 读取比 |
|-------|--------------|---------------|---------|-------------|--------------|---------|
| none | **1557ms** | 1806ms | **1.16x** | **1.25ms** | 213ms | **170x**\* |
| snappy | **1002ms** | 1649ms | **1.65x** | **78ms** | 384ms | **4.91x** |
| zstd | **1311ms** | 1796ms | **1.37x** | **76.8ms** | 369ms | **4.81x** |
| lz4 | **1521ms** | 1676ms | **1.10x** | **59.1ms** | 281ms | **4.76x** |

</div>

\* 零拷贝 mmap 路径

</details>

<details>
<summary><b>完整 ARM 结果</b>（Apple M3，macOS）</summary>

<br/>

*Carquet 0.7.2，MacBook Air M3，16GB 内存，macOS 26.5，Arrow C++ 24.0.0，PyArrow 23.0.1，ZSTD 级别 1。1000 万行的表格测量于 2026-10-20；100 万行和 10 万行的表格为 2026-10-08 的 0.7.1 测试结果（未重新测量）。*

#### 1000 万行 vs Arrow C++

<div align="center">

| 编解码器 | Carquet 写入 | Arrow C++ 写入 | 写入比 | Carquet 读取 | Arrow C++ 读取 | 读取比 | 大小 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **74.7ms** | 137ms | **1.83x** | **0.36ms** | 14.0ms | **38.9x**\* | 200.0MB |
| snappy | **52.7ms** | 244ms | **4.63x** | **12.4ms** | 20.6ms | **1.66x** | 131.2MB |
| zstd | **60.4ms** | 342ms | **5.66x** | **19.8ms** | 26.9ms | **1.36x** | 99.9MB |
| lz4 | **57.2ms** | 256ms | **4.47x** | **9.22ms** | 16.2ms | **1.76x** | 128.9MB |

</div>

<br/>

#### 100 万行 vs Arrow C++

<div align="center">

| 编解码器 | Carquet 写入 | Arrow C++ 写入 | 写入比 | Carquet 读取 | Arrow C++ 读取 | 读取比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **6.48ms** | 14.82ms | **2.29x** | **0.06ms** | 1.63ms | **27.2x**\* |
| snappy | **14.28ms** | 26.44ms | **1.85x** | **2.14ms** | 2.85ms | **1.33x** |
| zstd | **14.76ms** | 36.17ms | **2.45x** | **2.70ms** | 3.46ms | **1.28x** |
| lz4 | **12.30ms** | 26.10ms | **2.12x** | **1.04ms** | 1.87ms | **1.80x** |

</div>

<br/>

#### 10 万行 vs Arrow C++

<div align="center">

| 编解码器 | Carquet 写入 | Arrow C++ 写入 | 写入比 | Carquet 读取 | Arrow C++ 读取 | 读取比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **1.07ms** | 1.79ms | **1.67x** | **0.02ms** | 0.25ms | **12.5x**\* |
| snappy | **1.60ms** | 2.63ms | **1.64x** | **0.34ms** | 0.91ms | **2.68x** |
| zstd | **1.73ms** | 3.63ms | **2.10x** | **0.63ms** | 1.26ms | **2.00x** |
| lz4 | **1.52ms** | 2.62ms | **1.72x** | **0.25ms** | 0.55ms | **2.20x** |

</div>

<br/>

#### 同文件交叉读取（1000 万行）

两个库读取**同一个** Parquet 文件，这是最公平的同等条件对比。

<div align="center">

| 编解码器 | 写入方 | Carquet 读取 | Arrow C++ 读取 | 比值 |
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

| 编解码器 | Carquet 写入 | PyArrow 写入 | 写入比 | Carquet 读取 | PyArrow 读取 | 读取比 |
|-------|--------------|---------------|---------|-------------|--------------|---------|
| none | **74.7ms** | 202ms | **2.71x** | **0.36ms** | 36.2ms | **100.5x**\* |
| snappy | **52.7ms** | 312ms | **5.93x** | **12.4ms** | 46.6ms | **3.75x** |
| zstd | **60.4ms** | 410ms | **6.78x** | **19.8ms** | 75.6ms | **3.82x** |
| lz4 | **57.2ms** | 321ms | **5.62x** | **9.22ms** | 40.0ms | **4.33x** |

</div>

<br/>

#### 100 万行 vs PyArrow

<div align="center">

| 编解码器 | Carquet 写入 | PyArrow 写入 | 写入比 | Carquet 读取 | PyArrow 读取 | 读取比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **6.48ms** | 20.17ms | **3.11x** | **0.06ms** | 2.79ms | **46.5x**\* |
| snappy | **14.28ms** | 31.54ms | **2.21x** | **2.14ms** | 4.01ms | **1.87x** |
| zstd | **14.76ms** | 42.51ms | **2.88x** | **2.70ms** | 4.42ms | **1.64x** |
| lz4 | **12.30ms** | 32.03ms | **2.60x** | **1.04ms** | 3.22ms | **3.10x** |

</div>

<br/>

#### 10 万行 vs PyArrow

<div align="center">

| 编解码器 | Carquet 写入 | PyArrow 写入 | 写入比 | Carquet 读取 | PyArrow 读取 | 读取比 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **1.07ms** | 2.12ms | **1.98x** | **0.02ms** | 0.24ms | **12.0x**\* |
| snappy | **1.60ms** | 3.17ms | **1.98x** | **0.34ms** | 0.60ms | **1.76x** |
| zstd | **1.73ms** | 4.16ms | **2.40x** | **0.63ms** | 0.81ms | **1.29x** |
| lz4 | **1.52ms** | 3.21ms | **2.11x** | **0.25ms** | 0.42ms | **1.68x** |

</div>

\* 零拷贝 mmap 路径

</details>

<details>
<summary><b>自行运行基准测试</b></summary>

<br/>

```bash
# 以最高优化级别构建
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON -DCARQUET_BUILD_DEV=ON
cmake --build build -j$(nproc)

cd build
./benchmark_carquet                     # Carquet 单独运行
python3 ../benchmark/run_benchmark.py   # 完整对比（+ PyArrow，+ Arrow C++）

# 跳过 1 亿行（xlarge）配置 —— 每种编解码器会写入约 2GB 文件，
# 视硬件情况可能耗时 30 分钟以上
python3 ../benchmark/run_benchmark.py --skip-xlarge

# 覆盖 ZSTD 级别（默认：1）
CARQUET_BENCH_ZSTD_LEVEL=3 python3 ../benchmark/run_benchmark.py
```

要包含 Arrow C++ 对比，请在配置时启用其目标：

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON \
  -DCARQUET_BUILD_BENCHMARKS=ON \
  -DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON
cmake --build build -j$(nproc)

# 或指向自定义的 Arrow 安装
cmake -B build ... -DCARQUET_ARROW_CPP_ROOT=/path/to/arrow-prefix
```

Arrow C++ 基准测试使用底层 `parquet::ParquetFileReader` API（绕过 Arrow Table 物化开销）并采用并行行组读取器。**同文件交叉读取**模式让两个库读取完全相同的 Parquet 文件，消除页大小、编码和行组布局上的差异。两个基准测试使用相同的数据、行组大小、无字典、页校验和、mmap 读取、浮点数使用 BYTE_STREAM_SPLIT。

</details>

<br/>

---

<!-- ────────────────────────────  FEATURE SUPPORT  ──────────────────────────── -->

## Parquet 功能支持

本库目前实现的功能，按领域列出。

### 格式

<div align="center">

| 功能 | 状态 |
|---------|--------|
| 物理类型 | ✓  全部 8 种（从 BOOLEAN 到 FIXED_LEN_BYTE_ARRAY） |
| 逻辑类型 | ✓  STRING, DATE, TIME, TIMESTAMP, DECIMAL, UUID, JSON, INTERVAL, FLOAT16, VARIANT, GEOMETRY, GEOGRAPHY |
| 编码 | ✓  PLAIN, RLE, DICTIONARY, DELTA_BINARY_PACKED, DELTA_LENGTH_BYTE_ARRAY, DELTA_BYTE_ARRAY, BYTE_STREAM_SPLIT（读 + 写） |
| 数据页版本 | ✓  V1（默认）和 V2（读 + 写） |
| 压缩 | ✓  UNCOMPRESSED、SNAPPY、GZIP、LZ4（Hadoop 帧格式，编解码器 5）、LZ4_RAW（编解码器 7）、ZSTD、BROTLI（`CARQUET_ENABLE_BROTLI`，默认 ON）；LZO 可通过 `carquet_register_codec()` 注册 |
| 嵌套模式 | ✓  带定义/重复层级的组、列表、映射；写入时单层 LIST/MAP 自动拆分（`carquet_writer_write_list_column`），读取时重建 List；生成嵌套的 `ARROW:schema` |
| 加密 | ✗ 不支持；加密页脚文件（`PARE`）会以 `CARQUET_ERROR_NOT_IMPLEMENTED` 报告并附带明确的错误消息 |

</div>

<details>
<summary><b>索引与剪枝</b></summary>

<br/>

<div align="center">

| 功能 | 状态 |
|---------|--------|
| 布隆过滤器 | ✓  读取、写入和查询（`carquet_bloom_filter_check_*`）；自动用于行组剪枝 |
| 页索引 | ✓  列索引 + 偏移索引（读 + 写 + 按页访问统计信息） |
| 统计信息 | ✓  每个列块的 min/max/空值计数；字典列的精确 `distinct_count`；Parquet 2.9 SizeStatistics（未编码的字节数组字节数 + 层级直方图） |
| 谓词下推 | ✓  通过统计信息 + 布隆过滤器自动进行行组剪枝（无需回调）；通过列索引进行页级过滤（`carquet_batch_reader_set_page_filter`） |
| 地理空间统计信息 | ✓  GEOMETRY/GEOGRAPHY 的每列边界框和几何类型（`carquet_reader_geospatial_statistics`） |
| 列投影 | ✓  只读取选定的列 |
| 行窗口 | ✓  读取 `[offset, limit)` 范围的行，整页跳过（`carquet_batch_reader_set_row_range`） |

</div>

</details>

<details>
<summary><b>I/O 与可扩展性</b></summary>

<br/>

<div align="center">

| 功能 | 状态 |
|---------|--------|
| 追加 | ✓  向现有文件添加行组（`carquet_writer_open_append`） |
| 自定义编解码器 | ✓  为每个编解码器槽位注册自定义的压缩/解压实现（`carquet_register_codec`） |
| 键值元数据 | ✓  读写任意页脚元数据 |
| 字段级元数据 | ✓  通过 `ARROW:schema` 支持 Arrow `Field.custom_metadata`（变量标签/描述）（读 + 写） |
| 列级选项 | ✓  每列独立的编码、压缩、统计信息、布隆过滤器 |
| 缓冲区写入器 | ✓  将 Parquet 写入内存缓冲区 |
| CRC32 | ✓  页级校验（ARM 上硬件加速） |
| 内存映射 I/O | ✓  未压缩 PLAIN 数据的零拷贝读取 |
| 保留字典的读取 | ✓  返回 RLE 索引 + 字典而非物化后的值（`preserve_dictionaries`、`carquet_row_batch_column_dictionary`） |
| I/O 合并 | ✓  单次 I/O 预缓冲多列读取 |
| 推测性页脚读取 | ✓  大多数文件只需一次 I/O 即可打开 |
| OpenMP 并行读取 | ✓  可用时启用 |
| 并行写入 | ✓  跨列的页级并行编码/压缩（OpenMP）+ 后台行组 I/O（`async_io`） |

</div>

</details>

<br/>

---

<!-- ────────────────────────────  INTEROPERABILITY  ──────────────────────────── -->

## 互操作性

Carquet 生成的文件与 PyArrow、DuckDB、Spark 以及任何 Parquet 读取器完全兼容。

```python
import pyarrow.parquet as pq
table = pq.read_table("carquet_output.parquet")  # 直接可用
```

```sql
-- DuckDB
SELECT * FROM read_parquet('carquet_output.parquet');
```

<details>
<summary><b>互操作测试套件</b></summary>

<br/>

双向往返测试（C 写入 / Python 读取，以及反向）：

```bash
cmake -B build -DCARQUET_BUILD_INTEROP=ON && cmake --build build
python3 interop/run_interop.py
```

大型多行组 Snappy 文件（许多行组、宽 `BYTE_ARRAY` 列、超过 2 GiB 的页偏移）有自己的生成器，它用 PyArrow 构建文件并检查 carquet 能解码每一页：

```bash
CARQUET_BIN=build/carquet python3 interop/snappy_multi_rowgroup.py --rows 150000 --row-groups 92
```

</details>

<br/>

---

<!-- ────────────────────────────  CONTRIBUTING & LICENSE  ──────────────────────────── -->

## 贡献、引用与许可证

- **贡献**：[CONTRIBUTING.md](CONTRIBUTING.md) 涵盖环境搭建、代码风格以及每个 Pull Request 必须通过的检查关卡（测试、Sanitizer、模糊测试、互操作、基准测试）。发布说明见 [CHANGELOG.md](CHANGELOG.md)。
- **引用**：仓库包含 [`CITATION.cff`](CITATION.cff)；GitHub 的 "Cite this repository" 按钮可将其格式化为 BibTeX 或 APA。
- **许可证**：MIT，见 [LICENSE](LICENSE)。
