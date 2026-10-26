<h1 align="center">
  <img src="res/img/carquet_banner.png" alt="Carquet" width="840" />
</h1>

<p align="center">
  Apache Parquet 파일을 읽고 쓰기 위한 빠른 순수 C 라이브러리입니다.
</p>

<p align="center">
  <a href="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml"><img src="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml/badge.svg" alt="빌드" /></a>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue" alt="플랫폼" />
  <img src="https://img.shields.io/badge/C-C11-blue" alt="C 표준" />
  <a href="https://github.com/Vitruves/carquet/releases"><img src="https://img.shields.io/github/v/release/Vitruves/carquet?color=blue" alt="릴리스" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="라이선스" /></a>
  <a href="https://github.com/Vitruves/carquet/stargazers"><img src="https://img.shields.io/github/stars/Vitruves/carquet?style=flat&color=blue" alt="스타" /></a>
  <a href="https://github.com/Vitruves/carquet/graphs/traffic"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/Vitruves/carquet/traffic/clones.json" alt="클론" /></a>
  <br/>
  <img src="https://img.shields.io/badge/SIMD-SSE4.2%20%C2%B7%20AVX2%20%C2%B7%20AVX--512%20%C2%B7%20NEON%20%C2%B7%20SVE-555" alt="SIMD" />
  <img src="https://img.shields.io/badge/deps-zstd%20%C2%B7%20zlib%20%C2%B7%20lz4%20%C2%B7%20brotli-555" alt="의존성" />
</p>

<p align="center">
  <a href="#주요-특징">주요 특징</a> &nbsp;·&nbsp;
  <a href="#설치">설치</a> &nbsp;·&nbsp;
  <a href="#c-api">C API</a> &nbsp;·&nbsp;
  <a href="#cli-도구">CLI 도구</a> &nbsp;·&nbsp;
  <a href="#벤치마크">벤치마크</a> &nbsp;·&nbsp;
  <a href="#parquet-기능-지원">기능 지원</a> &nbsp;·&nbsp;
  <a href="#상호-운용성">상호 운용성</a> &nbsp;·&nbsp;
  <a href="#ready-for-ai-coding-agents">AI 에이전트</a> &nbsp;·&nbsp;
  <a href="docs/README.md">매뉴얼</a> &nbsp;·&nbsp;
  <a href="CHANGELOG.md">변경 이력</a>
</p>

<p align="center">
  <a href="README.md">English</a> &nbsp;·&nbsp;
  <a href="README.zh-CN.md">简体中文</a> &nbsp;·&nbsp;
  <a href="README.ja.md">日本語</a> &nbsp;·&nbsp;
  <b>한국어</b> &nbsp;·&nbsp;
  <a href="README.de.md">Deutsch</a> &nbsp;·&nbsp;
  <a href="README.fr.md">Français</a> &nbsp;·&nbsp;
  <a href="README.pt-BR.md">Português (Brasil)</a>
</p>

<br/>

---

<!-- ────────────────────────────  HIGHLIGHTS  ──────────────────────────── -->

## 주요 특징

**작고 자체 완결적**

- **순수 C11**로 작성되었으며 외부 의존성은 네 가지(zstd, zlib, lz4, brotli)뿐이고, 없으면 CMake가 모두 자동으로 가져옵니다. Brotli는 끌 수 있습니다
- **약 500 KB 라이브러리**(스트립된 Release 빌드, 코덱은 별도 링크) vs Arrow의 약 50 MB 이상
- 파일 검사(`schema`, `info`, `head`, `tail`, `stat`, …)와 C 코드 생성(`codegen`)을 위한 **내장 CLI**

**빠릅니다**

- 비압축 데이터에서 Arrow C++보다 **39배 빠른 읽기**(mmap 제로 카피), PyArrow보다 **100배 빠름**
- 동일한 파일에 대해 Arrow C++보다 **1.3–2.1배 빠른 압축 읽기**(교차 읽기 벤치마크, Apple M3; x86에서는 1.6–2.6배)
- Apple M3에서 Arrow C++보다 **1.8–5.7배 빠른 쓰기**(이전 x86 실행에서는 1.0–1.6배)
- 비압축 1천만 행을 **0.36 ms**에 읽음(Apple M3, mmap 제로 카피)
- 런타임 감지와 스칼라 폴백을 갖춘 SIMD 최적화(SSE4.2, AVX2, AVX-512, NEON, SVE)

**완전하고 호환됩니다**

- 폭넓은 Parquet 지원: 모든 물리 타입과 인코딩, Snappy / GZIP / LZ4 / ZSTD / Brotli, 중첩 스키마, 블룸 필터, 페이지 인덱스([전체 매트릭스](#parquet-기능-지원))
- PyArrow, DuckDB, Spark와 별도 설정 없이 호환
- Arrow 생태계와 의존성 없이 복사를 최소화하며 데이터를 주고받기 위한 [Arrow C Data Interface](docs/reading.md#export-to-arrow-c-data-interface) 브리지(`carquet_arrow_export_*` / `carquet_arrow_import_*` / `carquet_writer_write_arrow` / `carquet_reader_read_arrow`), **임의 깊이의 중첩 읽기/쓰기** 지원(struct / list / map, 어떤 깊이든)

<a id="ready-for-ai-coding-agents"></a>

**AI 코딩 에이전트 준비 완료**

- [`.agents/skills/carquet/`](.agents/skills/carquet/)는 에이전트에게 carquet을 프로젝트에 통합하는 방법을 가르치는 [Agent Skill](https://agentskills.io)(Claude Code, Codex CLI)입니다: 빌드 연결, 작업별로 어떤 API를 쓸지, 완전한 레시피, 레벨 / 널 비트맵 / 수명 규칙, 타입 매핑과 흔한 함정을 다룹니다. 리포지토리에 복사해 두면 에이전트가 스스로 사용합니다:

  ```bash
  cp -r path/to/carquet/.agents/skills/carquet .claude/skills/   # Claude Code
  cp -r path/to/carquet/.agents/skills/carquet .agents/skills/   # Codex CLI
  ```

<br/>

### 30초로 보는 carquet

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

배치와 널 허용 컬럼을 포함한 전체 쓰기 및 읽기 루프는 [C API](#c-api)에 있습니다.

<br/>

### carquet이 적합한 자리

carquet은 작고 임베드 가능한 선택지입니다. 자동으로 가져오는 네 가지 의존성을 가진 약 500 KB의 C11 라이브러리로 1분 안에 빌드되는 반면, Arrow C++는 50 MB 이상, PyArrow 휠은 약 100 MB에 달합니다. 아래 1천만 행 벤치마크에서 압축 파일을 Arrow C++보다 1.8–5.7배 빠르게 쓰고, 바로 그 동일한 파일을 1.3–2.1배 빠르게 읽습니다. 비압축 데이터는 mmap에서 곧바로 0.36 ms에 제공됩니다. Arrow의 생태계(컴퓨트 커널, 데이터셋, Flight, 암호화)가 필요하다면 여전히 Arrow가 올바른 선택입니다. carquet은 Parquet을 빠르고 정확하게 읽고 쓰기만 하면 되는 C 및 C++ 프로그램, 임베디드 시스템, 언어 런타임을 위한 것입니다.

<br/>

---

<!-- ────────────────────────────  INSTALLATION  ──────────────────────────── -->

## 설치

Carquet은 Linux, macOS, Windows에서 CMake(기본) 또는 [xmake](https://xmake.io)로 빌드합니다.

### 요구 사항

- C11 컴파일러(GCC 4.9+, Clang 3.4+, MSVC 2015+)
- CMake 3.16+(또는 [xmake](https://xmake.io), [xmake로 빌드하기](#xmake로-빌드하기) 참조)
- zstd, zlib, lz4, brotli(없으면 자동으로 가져옴; Brotli는 선택 사항, `-DCARQUET_ENABLE_BROTLI=OFF`)
- OpenMP(선택 사항, 병렬 컬럼 읽기 및 쓰기용)

<br/>

### 빠른 시작 (make)

`make` 래퍼는 최적화된 CMake 빌드와 `/usr/local` 설치를 수행합니다:

```bash
git clone https://github.com/Vitruves/carquet.git
cd carquet
make                              # 최적화 빌드 (`make help`로 전체 타깃 목록 확인)
sudo make install                 # /usr/local에 설치 (PREFIX=/opt/carquet으로 변경 가능)
```

<br/>

### 빌드 및 설치 (CMake)

특정 빌드 옵션이나 사용자 지정 접두사가 필요하면 CMake를 직접 호출합니다:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release   # 옵션 추가, 예: -DCARQUET_BUILD_SHARED=ON
cmake --build build -j$(nproc)
sudo cmake --install build --prefix /usr/local
```

어느 방법이든 다음을 설치합니다:

- `libcarquet.a`(`-DCARQUET_BUILD_SHARED=ON`이면 `.so` / `.dylib`)
- `include/carquet/` 헤더
- `carquet` CLI 바이너리
- `carquet.pc`(pkg-config) 및 `find_package(carquet)`용 CMake 패키지 설정

설치 후에는 `-lcarquet`으로 링크하거나 `pkg-config --cflags --libs carquet`으로 플래그를 가져옵니다.

<br/>

### Python (pip)

`carquet` 휠은 라이브러리를 번들하고 Arrow C Data Interface를 통해 노출하므로, Parquet 디코딩 로직을 중복하지 않고도 테이블이 PyArrow와 양방향으로 오갑니다:

```bash
pip install carquet
```

```python
import carquet
table = carquet.read_table("data.parquet")                 # pyarrow.Table
carquet.write_table(table, "out.parquet", compression="zstd")
print(carquet.info("out.parquet").num_rows)
```

패키지는 [`python/`](python/)에 있습니다(scikit-build-core, ctypes, PyArrow는 런타임 선택 사항). API와 휠 빌드는 해당 [README](python/README.md)를 참조하세요.

<br/>

### 빌드 옵션

CMake 구성 단계에 전달합니다(`-DOPTION=ON|OFF`):

<div align="center">

| 옵션 | 기본값 | 설명 |
|--------|---------|-------------|
| `CARQUET_BUILD_DEV` | OFF | 모든 것을 빌드(테스트, 예제, 벤치마크, 상호 운용성) |
| `CARQUET_BUILD_TESTS` | OFF | 테스트 스위트만 빌드 |
| `CARQUET_BUILD_CLI` | ON | `carquet` CLI 도구 빌드 |
| `CARQUET_BUILD_SHARED` | OFF | 정적 대신 공유 라이브러리 빌드 |
| `CARQUET_NATIVE_ARCH` | OFF | 최대 성능을 위한 `-march=native`(호스트 전용 바이너리) |
| `CARQUET_ENABLE_SVE` | OFF | ARM SVE(실험적) |
| `CARQUET_ENABLE_BROTLI` | ON | BROTLI 코덱(libbrotli를 찾거나 가져옴; OFF이면 의존성 제거) |

</div>

모든 x86 SIMD(SSE, AVX, AVX2, AVX-512)와 ARM NEON은 자동 감지되어 기본적으로 활성화됩니다.

<details>
<summary><b>전체 빌드 옵션</b></summary>

<br/>

<div align="center">

| 옵션 | 기본값 | 설명 |
|--------|---------|-------------|
| `CARQUET_BUILD_EXAMPLES` | OFF | 예제 프로그램 빌드 |
| `CARQUET_BUILD_BENCHMARKS` | OFF | 벤치마크 및 프로파일링 프로그램 빌드 |
| `CARQUET_BUILD_ARROW_CPP_BENCHMARK` | OFF | 선택적 Arrow C++ 비교 벤치마크 |
| `CARQUET_BUILD_INTEROP` | OFF | 상호 운용성 테스트 빌드 |
| `CARQUET_BUILD_FUZZ` | OFF | 퍼즈 타깃 빌드 |
| `CARQUET_BUILD_DOCS` | OFF | `docs` 타깃 추가(Doxygen 필요) |
| `CARQUET_BUNDLE_DEPS` | OFF | zstd/zlib/lz4/brotli를 항상 가져와 정적으로 번들(이식 가능한 릴리스 아티팩트) |
| `CARQUET_ZSTD_TARGET` / `CARQUET_ZLIB_TARGET` / `CARQUET_LZ4_TARGET` / `CARQUET_BROTLI_TARGET` | *(비어 있음)* | 해당 라이브러리를 검색하거나 가져오는 대신 링크할 기존 CMake 타깃 |
| `CARQUET_ARROW_CPP_ROOT` | *(비어 있음)* | 비교 벤치마크용 사용자 지정 Arrow C++ 설치 경로 |
| `CARQUET_ENABLE_SSE` | ON | SSE 최적화(x86, 자동 감지) |
| `CARQUET_ENABLE_AVX` | ON | AVX 최적화(x86, 자동 감지) |
| `CARQUET_ENABLE_AVX2` | ON | AVX2 최적화(x86, 자동 감지) |
| `CARQUET_ENABLE_AVX512` | ON | AVX-512 최적화(x86, 자동 감지) |
| `CARQUET_ENABLE_NEON` | ON | NEON 최적화(ARM, 자동 감지) |

</div>

</details>

<details>
<summary><b>직접 준비한 zstd / zlib / lz4 사용하기</b> (더 큰 CMake 프로젝트에 carquet 임베드)</summary>

<br/>

carquet을 `add_subdirectory()`나 `FetchContent`로 가져오면 각 압축 라이브러리를 다음 순서로 결정합니다:

1. `CARQUET_ZSTD_TARGET`, `CARQUET_ZLIB_TARGET` 또는 `CARQUET_LZ4_TARGET`에 지정된 타깃(설정된 경우).
2. 시스템 사본(`CARQUET_BUNDLE_DEPS`가 ON이면 건너뜀).
3. 직접 가져와 빌드하는 사본. 프로젝트가 이미 해당 라이브러리를 일반적인 타깃 이름(`libzstd_static`, `zlibstatic`, `lz4_static`)으로 빌드하고 있다면, carquet은 두 번째 사본을 가져오는 대신 그 타깃을 재사용합니다.

따라서 이 라이브러리들을 벤더링하고 시스템 사본이 없는 프로젝트는 아무 설정도 필요 없습니다. 예를 들어 시스템 zlib 대신 자체 zlib을, 또는 정적 빌드 대신 공유 빌드를 명시적으로 선택하려면 타깃 이름을 지정합니다:

```cmake
set(CARQUET_ZLIB_TARGET zlibstatic)
add_subdirectory(carquet)
```

아무것도 편집하지 않고도 같은 일을 할 수 있습니다: `cmake -DCARQUET_ZLIB_TARGET=zlib ...`. 세 변수는 서로 독립적입니다. 지정된 타깃은 검색이나 가져오기 없이 그대로 링크되며, 타깃이 아닌 이름은 조용한 폴백이 아니라 구성 오류가 됩니다. 타깃은 임포트된 것이든 트리 내에서 빌드된 것이든, 정적이든 공유든 상관없지만, 인클루드 디렉터리를 사용 요구 사항(usage requirement)으로 가지고 있어야 합니다.

</details>

<details>
<summary><b>개발 빌드</b> (테스트, 예제, 벤치마크)</summary>

<br/>

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCARQUET_BUILD_DEV=ON
cmake --build build -j$(nproc)
cd build && ctest --output-on-failure
```

[CONTRIBUTING.md](CONTRIBUTING.md)에 풀 리퀘스트별 게이트(테스트, 새니타이저, 퍼징, 상호 운용성, 벤치마크)가 설명되어 있습니다.

</details>

<details>
<summary><b>API 문서</b> (Doxygen)</summary>

<br/>

공개 헤더에는 [Doxygen](https://www.doxygen.nl) 주석이 달려 있습니다. HTML 레퍼런스를 생성하려면 `doxygen`이 설치되어 있어야 하며(인클루드/의존성 다이어그램을 위한 `graphviz`는 선택 사항), 일반 빌드에는 둘 다 필요하지 않습니다.

```bash
# CMake — 구성 시 타깃을 활성화한 뒤 빌드
cmake -B build -DCARQUET_BUILD_DOCS=ON
cmake --build build --target docs

# xmake — 독립 태스크, 구성 플래그 불필요
xmake docs
```

둘 다 `build/docs/html/index.html`에 출력합니다. Doxygen이 설치되어 있지 않으면 타깃/태스크는 그 사실만 보고하고 아무것도 하지 않습니다.

</details>

<details>
<summary><b>AI 코딩 에이전트</b> (Claude Code, Codex: 바로 쓸 수 있는 스킬)</summary>

<br/>

[`.agents/skills/carquet/`](.agents/skills/carquet/)는 코딩 에이전트에게 carquet을 프로젝트에 통합하는 방법을 가르치는 [Agent Skill](https://agentskills.io)입니다: 빌드 연결(CMake, FetchContent, pkg-config, xmake), 어떤 API를 어디에 쓸지, 완전한 레시피, 레벨 / 널 비트맵 / 수명 규칙, 타입 매핑, 그리고 흔한 함정과 해결책을 다룹니다. 자신의 리포지토리에 복사하거나 링크해 두면 에이전트가 자동으로 인식합니다:

```bash
# Claude Code
mkdir -p .claude/skills && cp -r path/to/carquet/.agents/skills/carquet .claude/skills/
# Codex CLI
mkdir -p .agents/skills && cp -r path/to/carquet/.agents/skills/carquet .agents/skills/
# 벤더링한 carquet: 심볼릭 링크로 동기화 유지
ln -s ../../third_party/carquet/.agents/skills/carquet .claude/skills/carquet
```

같은 파일을 개인 스킬(`~/.claude/skills/`, `~/.codex/skills/`)로도 사용할 수 있습니다.

</details>

<br/>

### xmake로 빌드하기

CMake의 대안으로 동일한 옵션과 기본값을 가진 [xmake](https://xmake.io) 빌드(`xmake.lua`)가 제공됩니다. zstd/zlib/lz4는 정적으로 링크되므로 바이너리가 자체 완결적입니다.

```bash
xmake                    # 정적 라이브러리 + `carquet` CLI 빌드 (release)
xmake f --dev=y && xmake # 테스트, 예제, 벤치마크, 상호 운용성 추가
xmake test               # 테스트 스위트 실행
```

<details>
<summary><b>구성 옵션</b> (<code>xmake f --option=y|n</code>)</summary>

<br/>

옵션은 CMake 옵션과 동일하되 `CARQUET_` 접두사가 없고 소문자입니다:

<div align="center">

| xmake 옵션 | 기본값 | 설명 |
|--------------|---------|-------------|
| `--dev` | n | 테스트, 예제, 벤치마크, 상호 운용성 빌드 |
| `--tests` / `--examples` / `--benchmarks` / `--interop` | n | 그룹 하나만 개별 빌드 |
| `--cli` | y | `carquet` CLI 도구 빌드 |
| `--shared` | n | 정적 대신 공유 라이브러리 빌드 |
| `--openmp` | y | OpenMP 병렬 컬럼 읽기(사용 불가 시 자동 비활성화) |
| `--native_arch` | n | 최대 성능을 위한 `-march=native`(호스트 전용 바이너리) |
| `--sse` / `--avx` / `--avx2` / `--avx512` / `--neon` | y | SIMD 명령어 집합(자동 감지) |
| `--sve` | n | ARM SVE(실험적) |
| `--brotli` | y | BROTLI 코덱(xmake 저장소의 libbrotli) |
| `--fuzz` | n | 퍼즈 타깃 빌드(`--toolchain=clang` 사용) |

</div>

```bash
xmake f -m release --shared=y            # 공유 라이브러리
xmake f --dev=y --avx512=n && xmake      # 개발 빌드, AVX-512 비활성화
```

</details>

<br/>

---

<!-- ────────────────────────────  C API  ──────────────────────────── -->

## C API

아래에 쓰기와 읽기 예제가 하나씩 있습니다. 나머지는 모두 [`docs/`의 매뉴얼](docs/README.md)에서 다룹니다.

### Parquet 파일 쓰기

```c
#include <carquet/carquet.h>

int main(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    // 스키마 정의
    carquet_schema_t* schema = carquet_schema_create(&err);
    carquet_schema_add_column(schema, "id",    CARQUET_PHYSICAL_INT64,  NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    // 라이터 구성
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;

    // 쓰기
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

### Parquet 파일 읽기

```c
#include <carquet/carquet.h>
#include <stdio.h>

int main(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    // 최상의 읽기 성능을 위해 mmap으로 열기
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = true;

    carquet_reader_t* r = carquet_reader_open("output.parquet", &opts, &err);
    if (!r) { printf("Error: %s\n", err.message); return 1; }

    printf("Rows: %lld, Columns: %d\n",
           (long long)carquet_reader_num_rows(r),
           carquet_reader_num_columns(r));

    // 효율적인 순회를 위한 배치 리더
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
        // ids[0..n-1] 처리 ...
        carquet_row_batch_free(batch);
        batch = NULL;
    }

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    return 0;
}
```

<br/>

### 더 많은 레시피

플랫 읽기/쓰기를 넘어서는 모든 내용은 매뉴얼에 있으며, 각 페이지는 실행 가능한 예제로 연결됩니다:

<div align="center">

| 하려는 작업 | 참조 |
|---|---|
| 널 허용 컬럼, 행 그룹, 버퍼 출력 | [`docs/writing.md`](docs/writing.md) |
| 리스트, 맵, 그룹, definition/repetition 레벨 | [`docs/nested-data.md`](docs/nested-data.md) |
| 컬럼 프로젝션, 통계, 메타데이터 검사 | [`docs/reading.md`](docs/reading.md) |
| 프레디킷 푸시다운, 페이지 수준 필터링 | [`docs/reading.md`](docs/reading.md) |
| 기존 파일에 행 그룹 추가 | [`docs/writing.md`](docs/writing.md#append-to-an-existing-file) |
| 압축, 사용자 지정 코덱, 라이터 튜닝 | [`docs/writing.md`](docs/writing.md), [`docs/performance.md`](docs/performance.md) |
| mmap, 제로 카피, 프리버퍼링, I/O 병합 | [`docs/performance.md`](docs/performance.md) |
| 오류 코드와 복구 힌트 | [`docs/error-handling.md`](docs/error-handling.md) |

</div>

<br/>

### API 레퍼런스

전체 API는 [`include/carquet/carquet.h`](include/carquet/carquet.h)에 있습니다. 주요 타입:

<div align="center">

| 타입 | 용도 |
|------|---------|
| `carquet_reader_t` | 파일 리더(경로, `FILE*` 또는 메모리 버퍼에서 열기) |
| `carquet_writer_t` | 파일 라이터 |
| `carquet_batch_reader_t` | 고수준 배치 순회 |
| `carquet_schema_t` | 스키마 정의 및 인트로스펙션 |
| `carquet_error_t` | 풍부한 오류 정보(코드, 메시지, 소스 위치, 복구 힌트) |

</div>

[매뉴얼](docs/README.md)은 언제 어떤 인터페이스를 써야 하는지 설명하고, 소스 레이아웃과 아키텍처는 [CONTRIBUTING.md](CONTRIBUTING.md)에 문서화되어 있습니다. 탐색 가능한 HTML 레퍼런스는 Doxygen으로 헤더에서 생성할 수 있습니다. [빌드 옵션](#빌드-옵션) 아래의 *API 문서* 접힘 항목을 참조하세요.

<details>
<summary><b>예제 애플리케이션: MetricStore</b> (<code>mocklib/</code>)</summary>

<br/>

[`mocklib/`](mocklib/)은 carquet 위에 구축된 완전하고 자체 완결적인 예제 애플리케이션 **MetricStore**로, 실제 사용 사례의 참고 자료이자 공개 API의 엔드투엔드 통합 테스트로 유용합니다.

시계열 텔레메트리 저장소(이벤트를 수집한 뒤 검사하고 질의)를 모델링하며, 그 과정에서 **carquet의 공개 함수 143개 중 117개(약 82 %)**를 사용합니다:

- **쓰기 측**: 논리 타입을 포함한 스키마 구성, 컬럼별 인코딩 / 압축 / 블룸 튜닝, 페이지 인덱스와 통계
- **읽기 측**: 컬럼 프로젝션, 프레디킷 푸시다운, 페이지 수준 필터링, 블룸 멤버십, Arrow C Data Interface 브리지, 중첩 `LIST` 재구성, 버퍼/추가 I/O 및 메타데이터 인트로스펙션

여느 다운스트림 프로젝트와 같은 방식으로 carquet을 링크하고(`find_package(carquet)` / `add_subdirectory`), 자체 검증 라운드트립 테스트를 포함하며, `write_sample` 바이너리는 외부 검사를 위한 파일을 생성합니다. [`mocklib/README.md`](mocklib/README.md)를 참조하세요.

</details>

<br/>

---

<!-- ────────────────────────────  CLI TOOL  ──────────────────────────── -->

## CLI 도구

Carquet에는 Parquet 파일을 검사하고 C 리더 코드를 생성하는 명령줄 도구가 포함되어 있습니다. 라이브러리와 함께 기본으로 빌드 및 설치됩니다.

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

일반적인 사용법:

```bash
carquet schema data.parquet
carquet head -n 20 data.parquet
carquet stat data.parquet
carquet validate data.parquet
```

<br/>

### 행 필터링

`cat`, `count`, `head`, `tail`, `sample`, `export`는 `-p / --filter EXPR`를 받아 행 프레디킷을 페이지 수준까지 푸시다운합니다. 컬럼 인덱스의 min/max가 프레디킷과 일치할 수 있는 페이지만 압축 해제됩니다.

```bash
carquet cat -p "price > 100 AND status = 'active'" data.parquet
carquet count --filter "id >= 1000" data.parquet
carquet export --filter "ts IS NOT NULL" -c id,ts data.parquet
```

문법은 `column OP value [AND column OP value]...`이며 `OP` ∈ {`=`, `==`, `!=`, `<>`, `<`, `<=`, `>`, `>=`}이고, 추가로 `column IS NULL` / `column IS NOT NULL`을 지원합니다.

> 필터링하려면 파일에 페이지 인덱스가 있어야 합니다(`write_page_index = true`).

<br/>

### 코드 생성

어떤 Parquet 파일의 스키마에서든 완전하고 컴파일 가능한 C 리더를 생성합니다:

```bash
carquet codegen -f data.parquet -o reader.c
# Generated: reader.c
# Compile:   clang -o reader reader.c -I.../include -L.../build -lcarquet ...

./reader                    # data.parquet 읽기 (기본값으로 내장됨)
./reader other.parquet      # 다른 파일로 대체
```

<details>
<summary><b>codegen 옵션</b></summary>

<br/>

<div align="center">

| 플래그 | 설명 |
|------|-------------|
| `-f`, `--file FILE` | 스키마를 검사할 Parquet 파일 |
| `-o`, `--output FILE` | 출력 소스 파일(기본값: stdout) |
| `--mmap` | 생성된 코드에서 메모리 매핑 I/O 사용 |
| `--skeleton` | 사용자 지정 로직을 위한 빈 `process_batch` 생성 |
| `-c`, `--columns COLS` | 쉼표로 구분된 컬럼 필터 |
| `-b`, `--batch-size N` | 배치 크기(기본값: 1024) |

</div>

</details>

<br/>

---

<!-- ────────────────────────────  BENCHMARKS  ──────────────────────────── -->

## 벤치마크

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="res/img/bench_10m_m3_dark.png">
    <img src="res/img/bench_10m_m3_light.png" alt="Apple M3에서 1천만 행의 쓰기 및 읽기 시간: 코덱별 carquet vs Arrow C++ vs PyArrow" width="900">
  </picture>
</p>

1천만 행 기준(가장 대표적인 크기); 비율이 높을수록 Carquet이 더 빠릅니다.

ARM(Apple M3): Carquet 0.7.2 vs Arrow C++ 24.0.0, 2026-10-20 측정. x86(Xeon D-1531): Carquet 0.4.4 vs Arrow C++ 23.0.1; 이 실행은 0.7.1의 라이터 병렬화 이전에 측정된 것으로, 이후 다시 측정되지 않았습니다.

<div align="center">

| | x86 (Xeon D-1531) | | ARM (Apple M3) | |
|---|---|---|---|---|
| **코덱** | **쓰기** | **읽기** | **쓰기** | **읽기** |
| snappy | **1.55x** | **1.25x** | **4.63x** | **1.66x** |
| zstd | **1.31x** | **1.04x** | **5.66x** | **1.36x** |
| lz4 | **1.02x** | 0.83x | **4.47x** | **1.76x** |
| none | **1.13x** | **40.6x**\* | **1.83x** | **38.9x**\* |

</div>

\* 비압축 읽기는 mmap 제로 카피를 사용합니다. 방법론 노트를 참조하세요.

압축 읽기는 모든 값을 완전히 압축 해제하고 디코딩하며 지름길은 없습니다. 두 라이브러리 모두 동일한 시스템 lz4/zstd 공유 라이브러리를 사용하므로 코덱 자체의 속도는 동일합니다. 가장 의미 있는 비교는 (전체 결과에 있는) **동일 파일 교차 읽기** 표로, 두 라이브러리가 정확히 같은 Parquet 파일을 읽습니다. 이 동일 조건 테스트에서 Carquet은 압축 데이터를 Arrow C++보다 **1.3–2.1배 빠르게** 읽습니다(x86 실행: 1.6–2.6배).

<details>
<summary><b>벤치마크 방법론</b></summary>

<br/>

모든 벤치마크는 동일한 데이터(결정론적 LCG PRNG), 동일한 Parquet 설정(딕셔너리 없음, 부동소수점에 BYTE_STREAM_SPLIT, 페이지 체크섬, mmap 읽기)을 사용하며, 11–51회 반복의 절사 중앙값을 취하고, 쓰기와 읽기 단계 사이에 OS 페이지 캐시를 비우며 구성 간에 쿨다운을 둡니다. 스키마: 3개 컬럼(INT64, DOUBLE, INT32). Arrow C++ 저수준 Parquet 리더(Arrow Table 구체화를 우회) 및 PyArrow와 비교했습니다.

**동일 파일 교차 읽기** 벤치마크가 가장 공정한 비교입니다. 두 라이브러리가 정확히 같은 Parquet 파일(한쪽이 쓰고 양쪽이 읽음)을 읽으므로 페이지 크기, 인코딩 선택, 행 그룹 레이아웃의 차이가 제거됩니다.

\*로 표시된 **비압축 읽기**는 Carquet의 **mmap 제로 카피 경로**를 사용합니다. PLAIN 인코딩, 비압축, 고정 크기, required 컬럼의 경우 배치 리더가 memcpy 없이 메모리 매핑된 파일을 직접 가리키는 포인터를 반환합니다. Arrow는 항상 자체 버퍼로 구체화합니다. **압축 읽기 수치가 엔드투엔드 읽기 처리량을 가장 대표적으로 보여주는 척도입니다.**

</details>

<details>
<summary><b>전체 x86 결과</b> (Intel Xeon D-1531, Linux)</summary>

<br/>

*12 스레드 @ 2.2GHz, 32GB RAM, Ubuntu 24.04, ZSTD 레벨 1*

#### 1천만 행 vs Arrow C++

<div align="center">

| 코덱 | Carquet 쓰기 | Arrow C++ 쓰기 | 쓰기 비율 | Carquet 읽기 | Arrow C++ 읽기 | 읽기 비율 | 크기 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **1557ms** | 1766ms | **1.13x** | **1.25ms** | 50.8ms | **40.6x**\* | 190.7MB |
| snappy | **1002ms** | 1549ms | **1.55x** | **78ms** | 97.8ms | **1.25x** | 125.1MB |
| zstd | **1311ms** | 1714ms | **1.31x** | **76.8ms** | 80.2ms | **1.04x** | 95.3MB |
| lz4 | **1521ms** | 1554ms | **1.02x** | 59.1ms | **49.0ms** | 0.83x | 122.9MB |

</div>

<br/>

#### 100만 행 vs Arrow C++

<div align="center">

| 코덱 | Carquet 쓰기 | Arrow C++ 쓰기 | 쓰기 비율 | Carquet 읽기 | Arrow C++ 읽기 | 읽기 비율 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **180ms** | 196ms | **1.09x** | **0.22ms** | 6.2ms | **28x**\* |
| snappy | **141ms** | 148ms | **1.05x** | **8.1ms** | 11.6ms | **1.44x** |
| zstd | **131ms** | 185ms | **1.41x** | 10.3ms | **9.1ms** | 0.88x |
| lz4 | **143ms** | 149ms | **1.04x** | 8.5ms | **6.1ms** | 0.72x |

</div>

<br/>

#### 10만 행 vs Arrow C++

<div align="center">

| 코덱 | Carquet 쓰기 | Arrow C++ 쓰기 | 쓰기 비율 | Carquet 읽기 | Arrow C++ 읽기 | 읽기 비율 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **14.1ms** | 18.4ms | **1.30x** | **0.11ms** | 2.18ms | **19.8x**\* |
| snappy | **10.1ms** | 10.6ms | **1.05x** | **1.27ms** | 5.97ms | **4.70x** |
| zstd | **8.7ms** | 14.1ms | **1.62x** | **1.58ms** | 3.88ms | **2.46x** |
| lz4 | **9.6ms** | 11.0ms | **1.14x** | **0.77ms** | 2.78ms | **3.61x** |

</div>

<br/>

#### 동일 파일 교차 읽기 (1천만 행)

두 라이브러리가 **동일한** Parquet 파일을 읽는, 가장 공정한 동일 조건 비교입니다.

<div align="center">

| 코덱 | 라이터 | Carquet 읽기 | Arrow C++ 읽기 | 비율 |
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

#### 1천만 행 vs PyArrow

<div align="center">

| 코덱 | Carquet 쓰기 | PyArrow 쓰기 | 쓰기 비율 | Carquet 읽기 | PyArrow 읽기 | 읽기 비율 |
|-------|--------------|---------------|---------|-------------|--------------|---------|
| none | **1557ms** | 1806ms | **1.16x** | **1.25ms** | 213ms | **170x**\* |
| snappy | **1002ms** | 1649ms | **1.65x** | **78ms** | 384ms | **4.91x** |
| zstd | **1311ms** | 1796ms | **1.37x** | **76.8ms** | 369ms | **4.81x** |
| lz4 | **1521ms** | 1676ms | **1.10x** | **59.1ms** | 281ms | **4.76x** |

</div>

\* 제로 카피 mmap 경로

</details>

<details>
<summary><b>전체 ARM 결과</b> (Apple M3, macOS)</summary>

<br/>

*Carquet 0.7.2, MacBook Air M3, 16GB RAM, macOS 26.5, Arrow C++ 24.0.0, PyArrow 23.0.1, ZSTD 레벨 1. 1천만 행 표는 2026-10-20에 측정했고, 100만 및 10만 행 표는 2026-10-08의 0.7.1 실행 결과입니다(재측정하지 않음).*

#### 1천만 행 vs Arrow C++

<div align="center">

| 코덱 | Carquet 쓰기 | Arrow C++ 쓰기 | 쓰기 비율 | Carquet 읽기 | Arrow C++ 읽기 | 읽기 비율 | 크기 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **74.7ms** | 137ms | **1.83x** | **0.36ms** | 14.0ms | **38.9x**\* | 200.0MB |
| snappy | **52.7ms** | 244ms | **4.63x** | **12.4ms** | 20.6ms | **1.66x** | 131.2MB |
| zstd | **60.4ms** | 342ms | **5.66x** | **19.8ms** | 26.9ms | **1.36x** | 99.9MB |
| lz4 | **57.2ms** | 256ms | **4.47x** | **9.22ms** | 16.2ms | **1.76x** | 128.9MB |

</div>

<br/>

#### 100만 행 vs Arrow C++

<div align="center">

| 코덱 | Carquet 쓰기 | Arrow C++ 쓰기 | 쓰기 비율 | Carquet 읽기 | Arrow C++ 읽기 | 읽기 비율 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **6.48ms** | 14.82ms | **2.29x** | **0.06ms** | 1.63ms | **27.2x**\* |
| snappy | **14.28ms** | 26.44ms | **1.85x** | **2.14ms** | 2.85ms | **1.33x** |
| zstd | **14.76ms** | 36.17ms | **2.45x** | **2.70ms** | 3.46ms | **1.28x** |
| lz4 | **12.30ms** | 26.10ms | **2.12x** | **1.04ms** | 1.87ms | **1.80x** |

</div>

<br/>

#### 10만 행 vs Arrow C++

<div align="center">

| 코덱 | Carquet 쓰기 | Arrow C++ 쓰기 | 쓰기 비율 | Carquet 읽기 | Arrow C++ 읽기 | 읽기 비율 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **1.07ms** | 1.79ms | **1.67x** | **0.02ms** | 0.25ms | **12.5x**\* |
| snappy | **1.60ms** | 2.63ms | **1.64x** | **0.34ms** | 0.91ms | **2.68x** |
| zstd | **1.73ms** | 3.63ms | **2.10x** | **0.63ms** | 1.26ms | **2.00x** |
| lz4 | **1.52ms** | 2.62ms | **1.72x** | **0.25ms** | 0.55ms | **2.20x** |

</div>

<br/>

#### 동일 파일 교차 읽기 (1천만 행)

두 라이브러리가 **동일한** Parquet 파일을 읽는, 가장 공정한 동일 조건 비교입니다.

<div align="center">

| 코덱 | 라이터 | Carquet 읽기 | Arrow C++ 읽기 | 비율 |
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

#### 1천만 행 vs PyArrow

<div align="center">

| 코덱 | Carquet 쓰기 | PyArrow 쓰기 | 쓰기 비율 | Carquet 읽기 | PyArrow 읽기 | 읽기 비율 |
|-------|--------------|---------------|---------|-------------|--------------|---------|
| none | **74.7ms** | 202ms | **2.71x** | **0.36ms** | 36.2ms | **100.5x**\* |
| snappy | **52.7ms** | 312ms | **5.93x** | **12.4ms** | 46.6ms | **3.75x** |
| zstd | **60.4ms** | 410ms | **6.78x** | **19.8ms** | 75.6ms | **3.82x** |
| lz4 | **57.2ms** | 321ms | **5.62x** | **9.22ms** | 40.0ms | **4.33x** |

</div>

<br/>

#### 100만 행 vs PyArrow

<div align="center">

| 코덱 | Carquet 쓰기 | PyArrow 쓰기 | 쓰기 비율 | Carquet 읽기 | PyArrow 읽기 | 읽기 비율 |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **6.48ms** | 20.17ms | **3.11x** | **0.06ms** | 2.79ms | **46.5x**\* |
| snappy | **14.28ms** | 31.54ms | **2.21x** | **2.14ms** | 4.01ms | **1.87x** |
| zstd | **14.76ms** | 42.51ms | **2.88x** | **2.70ms** | 4.42ms | **1.64x** |
| lz4 | **12.30ms** | 32.03ms | **2.60x** | **1.04ms** | 3.22ms | **3.10x** |

</div>

<br/>

#### 10만 행 vs PyArrow

<div align="center">

| 코덱 | Carquet 쓰기 | PyArrow 쓰기 | 쓰기 비율 | Carquet 읽기 | PyArrow 읽기 | 읽기 비율 |
|-------|--------------|-----------------|---------|-------------|--------------|---------|
| none | **1.07ms** | 2.12ms | **1.98x** | **0.02ms** | 0.24ms | **12.0x**\* |
| snappy | **1.60ms** | 3.17ms | **1.98x** | **0.34ms** | 0.60ms | **1.76x** |
| zstd | **1.73ms** | 4.16ms | **2.40x** | **0.63ms** | 0.81ms | **1.29x** |
| lz4 | **1.52ms** | 3.21ms | **2.11x** | **0.25ms** | 0.42ms | **1.68x** |

</div>

\* 제로 카피 mmap 경로

</details>

<details>
<summary><b>직접 벤치마크 실행하기</b></summary>

<br/>

```bash
# 최대 최적화로 빌드
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON -DCARQUET_BUILD_DEV=ON
cmake --build build -j$(nproc)

cd build
./benchmark_carquet                     # Carquet 단독
python3 ../benchmark/run_benchmark.py   # 전체 비교 (+ PyArrow, + Arrow C++)

# 1억 행(xlarge) 구성 건너뛰기 — 코덱당 약 2GB 파일을 쓰며
# 하드웨어에 따라 30분 이상 걸릴 수 있음
python3 ../benchmark/run_benchmark.py --skip-xlarge

# ZSTD 레벨 변경 (기본값: 1)
CARQUET_BENCH_ZSTD_LEVEL=3 python3 ../benchmark/run_benchmark.py
```

Arrow C++ 비교를 포함하려면 구성 시 해당 타깃을 활성화합니다:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON \
  -DCARQUET_BUILD_BENCHMARKS=ON \
  -DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON
cmake --build build -j$(nproc)

# 또는 사용자 지정 Arrow 설치 경로 지정
cmake -B build ... -DCARQUET_ARROW_CPP_ROOT=/path/to/arrow-prefix
```

Arrow C++ 벤치마크는 병렬 행 그룹 리더와 함께 저수준 `parquet::ParquetFileReader` API(Arrow Table 구체화 오버헤드 우회)를 사용합니다. **동일 파일 교차 읽기** 모드에서는 두 라이브러리가 정확히 같은 Parquet 파일을 읽어 페이지 크기, 인코딩, 행 그룹 레이아웃의 차이를 제거합니다. 두 벤치마크 모두 동일한 데이터, 행 그룹 크기, 딕셔너리 없음, 페이지 체크섬, mmap 읽기, 부동소수점에 BYTE_STREAM_SPLIT을 사용합니다.

</details>

<br/>

---

<!-- ────────────────────────────  FEATURE SUPPORT  ──────────────────────────── -->

## Parquet 기능 지원

현재 라이브러리가 구현하는 내용을 영역별로 정리했습니다.

### 포맷

<div align="center">

| 기능 | 상태 |
|---------|--------|
| 물리 타입 | ✓  8개 모두(BOOLEAN부터 FIXED_LEN_BYTE_ARRAY까지) |
| 논리 타입 | ✓  STRING, DATE, TIME, TIMESTAMP, DECIMAL, UUID, JSON, INTERVAL, FLOAT16, VARIANT, GEOMETRY, GEOGRAPHY |
| 인코딩 | ✓  PLAIN, RLE, DICTIONARY, DELTA_BINARY_PACKED, DELTA_LENGTH_BYTE_ARRAY, DELTA_BYTE_ARRAY, BYTE_STREAM_SPLIT(읽기 + 쓰기) |
| 데이터 페이지 버전 | ✓  V1(기본값) 및 V2(읽기 + 쓰기) |
| 압축 | ✓  UNCOMPRESSED, SNAPPY, GZIP, LZ4(Hadoop 프레임, 코덱 5), LZ4_RAW(코덱 7), ZSTD, BROTLI(`CARQUET_ENABLE_BROTLI`, 기본 ON); LZO는 `carquet_register_codec()`으로 |
| 중첩 스키마 | ✓  definition/repetition 레벨을 사용하는 그룹, 리스트, 맵; 쓰기 시 단일 레벨 LIST/MAP 자동 셰레딩(`carquet_writer_write_list_column`)과 읽기 시 List 재구성; 중첩 `ARROW:schema` 생성 |
| 암호화 | ✗ 지원하지 않음; 암호화된 푸터 파일(`PARE`)은 명시적 메시지와 함께 `CARQUET_ERROR_NOT_IMPLEMENTED`로 보고됨 |

</div>

<details>
<summary><b>인덱싱 &amp; 프루닝</b></summary>

<br/>

<div align="center">

| 기능 | 상태 |
|---------|--------|
| 블룸 필터 | ✓  읽기, 쓰기, 질의(`carquet_bloom_filter_check_*`); 행 그룹 프루닝에 자동으로 활용 |
| 페이지 인덱스 | ✓  컬럼 인덱스 + 오프셋 인덱스(읽기 + 쓰기 + 페이지별 통계 접근) |
| 통계 | ✓  컬럼 청크별 min/max/null count; 딕셔너리 컬럼의 정확한 `distinct_count`; Parquet 2.9 SizeStatistics(비인코딩 바이트 배열 바이트 수 + 레벨 히스토그램) |
| 프레디킷 푸시다운 | ✓  통계 + 블룸 필터를 통한 자동 행 그룹 프루닝(콜백 불필요); 컬럼 인덱스를 통한 페이지 수준 필터링(`carquet_batch_reader_set_page_filter`) |
| 지리공간 통계 | ✓  GEOMETRY/GEOGRAPHY의 컬럼별 경계 상자 및 지오메트리 타입(`carquet_reader_geospatial_statistics`) |
| 컬럼 프로젝션 | ✓  선택한 컬럼만 읽기 |
| 행 범위 | ✓  행의 `[offset, limit)` 구간을 읽고 페이지를 통째로 건너뜀(`carquet_batch_reader_set_row_range`) |

</div>

</details>

<details>
<summary><b>I/O &amp; 확장성</b></summary>

<br/>

<div align="center">

| 기능 | 상태 |
|---------|--------|
| 추가(Append) | ✓  기존 파일에 행 그룹 추가(`carquet_writer_open_append`) |
| 사용자 지정 코덱 | ✓  코덱 슬롯별로 사용자 지정 압축/압축 해제 구현 등록(`carquet_register_codec`) |
| 키-값 메타데이터 | ✓  임의의 푸터 메타데이터 읽기 및 쓰기 |
| 필드별 메타데이터 | ✓  `ARROW:schema`를 통한 Arrow `Field.custom_metadata`(변수 레이블/설명)(읽기 + 쓰기) |
| 컬럼별 옵션 | ✓  컬럼별 인코딩, 압축, 통계, 블룸 필터 |
| 버퍼 라이터 | ✓  메모리 내 버퍼에 Parquet 쓰기 |
| CRC32 | ✓  페이지 수준 검증(ARM에서 하드웨어 가속) |
| 메모리 매핑 I/O | ✓  비압축 PLAIN 데이터의 제로 카피 읽기 |
| 딕셔너리 보존 읽기 | ✓  구체화된 값 대신 RLE 인덱스 + 딕셔너리 반환(`preserve_dictionaries`, `carquet_row_batch_column_dictionary`) |
| I/O 병합 | ✓  다중 컬럼 읽기를 단일 I/O로 프리버퍼링 |
| 추측적 푸터 읽기 | ✓  대부분의 파일을 단일 I/O로 열기 |
| OpenMP 병렬 읽기 | ✓  사용 가능할 때 |
| 병렬 쓰기 | ✓  컬럼 전반의 페이지 수준 병렬 인코딩/압축(OpenMP) + 백그라운드 행 그룹 I/O(`async_io`) |

</div>

</details>

<br/>

---

<!-- ────────────────────────────  INTEROPERABILITY  ──────────────────────────── -->

## 상호 운용성

Carquet 파일은 PyArrow, DuckDB, Spark 및 모든 Parquet 리더와 완전히 호환됩니다.

```python
import pyarrow.parquet as pq
table = pq.read_table("carquet_output.parquet")  # 그냥 동작합니다
```

```sql
-- DuckDB
SELECT * FROM read_parquet('carquet_output.parquet');
```

<details>
<summary><b>상호 운용성 테스트 스위트</b></summary>

<br/>

양방향 라운드트립(C가 쓰고 Python이 읽기, 그리고 그 반대):

```bash
cmake -B build -DCARQUET_BUILD_INTEROP=ON && cmake --build build
python3 interop/run_interop.py
```

대용량 다중 행 그룹 Snappy 파일(많은 행 그룹, 넓은 `BYTE_ARRAY` 컬럼, 2 GiB를 넘는 페이지 오프셋)에는 별도의 생성기가 있어, PyArrow로 파일을 만들고 carquet이 모든 페이지를 디코딩하는지 확인합니다:

```bash
CARQUET_BIN=build/carquet python3 interop/snappy_multi_rowgroup.py --rows 150000 --row-groups 92
```

</details>

<br/>

---

<!-- ────────────────────────────  CONTRIBUTING & LICENSE  ──────────────────────────── -->

## 기여, 인용, 라이선스

- **기여**: [CONTRIBUTING.md](CONTRIBUTING.md)는 설정, 코드 스타일, 풀 리퀘스트별 게이트(테스트, 새니타이저, 퍼징, 상호 운용성, 벤치마크)를 다룹니다. 릴리스 노트는 [CHANGELOG.md](CHANGELOG.md)에 있습니다.
- **인용**: [`CITATION.cff`](CITATION.cff)가 포함되어 있습니다. GitHub의 "Cite this repository" 버튼이 이를 BibTeX 또는 APA 형식으로 만들어 줍니다.
- **라이선스**: MIT, [LICENSE](LICENSE)를 참조하세요.
