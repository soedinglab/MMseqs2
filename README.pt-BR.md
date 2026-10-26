<h1 align="center">
  <img src="res/img/carquet_banner.png" alt="Carquet" width="840" />
</h1>

<p align="center">
  Uma biblioteca rápida, em C puro, para ler e escrever arquivos Apache Parquet.
</p>

<p align="center">
  <a href="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml"><img src="https://github.com/Vitruves/carquet/actions/workflows/cpp.yml/badge.svg" alt="Build" /></a>
  <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue" alt="Plataforma" />
  <img src="https://img.shields.io/badge/C-C11-blue" alt="Padrão C" />
  <a href="https://github.com/Vitruves/carquet/releases"><img src="https://img.shields.io/github/v/release/Vitruves/carquet?color=blue" alt="Release" /></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="Licença" /></a>
  <a href="https://github.com/Vitruves/carquet/stargazers"><img src="https://img.shields.io/github/stars/Vitruves/carquet?style=flat&color=blue" alt="Estrelas" /></a>
  <a href="https://github.com/Vitruves/carquet/graphs/traffic"><img src="https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/Vitruves/carquet/traffic/clones.json" alt="Clones" /></a>
  <br/>
  <img src="https://img.shields.io/badge/SIMD-SSE4.2%20%C2%B7%20AVX2%20%C2%B7%20AVX--512%20%C2%B7%20NEON%20%C2%B7%20SVE-555" alt="SIMD" />
  <img src="https://img.shields.io/badge/deps-zstd%20%C2%B7%20zlib%20%C2%B7%20lz4%20%C2%B7%20brotli-555" alt="Dependências" />
</p>

<p align="center">
  <a href="#destaques">Destaques</a> &nbsp;·&nbsp;
  <a href="#instalação">Instalação</a> &nbsp;·&nbsp;
  <a href="#api-c">API C</a> &nbsp;·&nbsp;
  <a href="#ferramenta-de-linha-de-comando">Ferramenta CLI</a> &nbsp;·&nbsp;
  <a href="#benchmarks">Benchmarks</a> &nbsp;·&nbsp;
  <a href="#suporte-a-recursos-do-parquet">Recursos suportados</a> &nbsp;·&nbsp;
  <a href="#interoperabilidade">Interoperabilidade</a> &nbsp;·&nbsp;
  <a href="#ready-for-ai-coding-agents">Agentes de IA</a> &nbsp;·&nbsp;
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

**Pequena e autocontida**

- **C11 puro** com quatro dependências externas (zstd, zlib, lz4, brotli), todas baixadas automaticamente pelo CMake quando ausentes; o Brotli pode ser desativado
- **Biblioteca de ~500 KB** (build Release com símbolos removidos, codecs linkados separadamente) contra 50 MB+ do Arrow
- **CLI integrada** para inspeção de arquivos (`schema`, `info`, `head`, `tail`, `stat`, …) e geração de código C (`codegen`)

**Rápida**

- **Leituras 39x mais rápidas** que o Arrow C++ em dados não comprimidos (mmap zero-copy), **100x mais rápidas** que o PyArrow
- **Leituras comprimidas 1,3–2,1x mais rápidas** que o Arrow C++ sobre o mesmo arquivo (benchmark de leitura cruzada, Apple M3; 1,6–2,6x em x86)
- **Escritas 1,8–5,7x mais rápidas** que o Arrow C++ em Apple M3 (1,0–1,6x na execução x86 mais antiga)
- Lê 10M de linhas não comprimidas em **0,36 ms** (mmap zero-copy em Apple M3)
- Otimizada com SIMD (SSE4.2, AVX2, AVX-512, NEON, SVE) com detecção em tempo de execução e fallbacks escalares

**Completa e compatível**

- Ampla cobertura do Parquet: todos os tipos físicos e encodings, Snappy / GZIP / LZ4 / ZSTD / Brotli, schemas aninhados, bloom filters, índices de página ([matriz completa](#suporte-a-recursos-do-parquet))
- Compatível com PyArrow, DuckDB e Spark sem configuração adicional
- Ponte com a [Arrow C Data Interface](docs/reading.md#export-to-arrow-c-data-interface) (`carquet_arrow_export_*` / `carquet_arrow_import_*` / `carquet_writer_write_arrow` / `carquet_reader_read_arrow`) para intercâmbio sem dependências e com poucas cópias com o ecossistema Arrow, com **leitura/escrita aninhada em profundidade arbitrária** (struct / list / map em qualquer profundidade)

<a id="ready-for-ai-coding-agents"></a>

**Pronta para agentes de IA de programação**

- [`.agents/skills/carquet/`](.agents/skills/carquet/) é uma [Agent Skill](https://agentskills.io) (Claude Code, Codex CLI) que ensina um agente a integrar o carquet ao seu projeto: configuração do build, qual API usar para cada tarefa, receitas completas, as regras de níveis / null-bitmap / tempo de vida, mapeamentos de tipos e as armadilhas habituais. Copie-a para o seu repositório e o agente a usa por conta própria:

  ```bash
  cp -r path/to/carquet/.agents/skills/carquet .claude/skills/   # Claude Code
  cp -r path/to/carquet/.agents/skills/carquet .agents/skills/   # Codex CLI
  ```

<br/>

### Trinta segundos de carquet

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

Os loops completos de escrita e leitura, com batches e colunas anuláveis, estão em [API C](#api-c).

<br/>

### Onde o carquet se encaixa

O carquet é a opção pequena e embarcável: uma biblioteca C11 de ~500 KB com quatro dependências baixadas automaticamente que compila em menos de um minuto, enquanto o Arrow C++ pesa 50 MB+ e um wheel do PyArrow cerca de 100 MB. No benchmark de 10M de linhas abaixo, ele escreve arquivos comprimidos 1,8–5,7x mais rápido que o Arrow C++ e os lê 1,3–2,1x mais rápido a partir do mesmíssimo arquivo; dados não comprimidos são servidos diretamente do mmap em 0,36 ms. O Arrow continua sendo a escolha certa quando você precisa do ecossistema dele (kernels de computação, datasets, Flight, criptografia); o carquet é para os programas em C e C++, sistemas embarcados e runtimes de linguagens que só precisam ler e escrever Parquet com rapidez e correção.

<br/>

---

<!-- ────────────────────────────  INSTALAÇÃO  ──────────────────────────── -->

## Instalação

O carquet compila com CMake (principal) ou [xmake](https://xmake.io), em Linux, macOS e Windows.

### Requisitos

- Compilador C11 (GCC 4.9+, Clang 3.4+, MSVC 2015+)
- CMake 3.16+ (ou [xmake](https://xmake.io), veja [Compilando com xmake](#compilando-com-xmake))
- zstd, zlib, lz4, brotli (baixados automaticamente se ausentes; Brotli opcional, `-DCARQUET_ENABLE_BROTLI=OFF`)
- OpenMP (opcional, para leitura e escrita paralela de colunas)

<br/>

### Início rápido (make)

O wrapper `make` conduz um build CMake otimizado e uma instalação em `/usr/local`:

```bash
git clone https://github.com/Vitruves/carquet.git
cd carquet
make                              # build otimizado (execute `make help` para listar todos os alvos)
sudo make install                 # instala em /usr/local (sobrescreva com PREFIX=/opt/carquet)
```

<br/>

### Build e instalação (CMake)

Invoque o CMake diretamente quando precisar de opções de build específicas ou de um prefixo personalizado:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release   # adicione opções, p. ex. -DCARQUET_BUILD_SHARED=ON
cmake --build build -j$(nproc)
sudo cmake --install build --prefix /usr/local
```

Qualquer um dos caminhos instala:

- `libcarquet.a` (ou `.so` / `.dylib` com `-DCARQUET_BUILD_SHARED=ON`)
- headers em `include/carquet/`
- binário da CLI `carquet`
- `carquet.pc` (pkg-config) e configuração de pacote CMake para `find_package(carquet)`

Após a instalação, linke com `-lcarquet` ou resolva as flags via `pkg-config --cflags --libs carquet`.

<br/>

### Python (pip)

O wheel `carquet` embute a biblioteca e a expõe pela Arrow C Data Interface, de modo que tabelas fluem de e para o PyArrow sem cópias da lógica de decodificação do Parquet:

```bash
pip install carquet
```

```python
import carquet
table = carquet.read_table("data.parquet")                 # pyarrow.Table
carquet.write_table(table, "out.parquet", compression="zstd")
print(carquet.info("out.parquet").num_rows)
```

O pacote fica em [`python/`](python/) (scikit-build-core, ctypes, PyArrow opcional em tempo de execução); veja o [README](python/README.md) dele para a API e o build do wheel.

<br/>

### Opções de build

Passe estas opções na etapa de configuração do CMake (`-DOPTION=ON|OFF`):

<div align="center">

| Opção | Padrão | Descrição |
|--------|---------|-------------|
| `CARQUET_BUILD_DEV` | OFF | Compila tudo (testes, exemplos, benchmarks, interop) |
| `CARQUET_BUILD_TESTS` | OFF | Compila apenas a suíte de testes |
| `CARQUET_BUILD_CLI` | ON | Compila a ferramenta CLI `carquet` |
| `CARQUET_BUILD_SHARED` | OFF | Compila uma biblioteca compartilhada em vez de estática |
| `CARQUET_NATIVE_ARCH` | OFF | `-march=native` para desempenho máximo (binário só para a máquina atual) |
| `CARQUET_ENABLE_SVE` | OFF | ARM SVE (experimental) |
| `CARQUET_ENABLE_BROTLI` | ON | Codec BROTLI (libbrotli encontrada ou baixada; OFF remove a dependência) |

</div>

Todo o SIMD x86 (SSE, AVX, AVX2, AVX-512) e o ARM NEON são detectados automaticamente e habilitados por padrão.

<details>
<summary><b>Todas as opções de build</b></summary>

<br/>

<div align="center">

| Opção | Padrão | Descrição |
|--------|---------|-------------|
| `CARQUET_BUILD_EXAMPLES` | OFF | Compila os programas de exemplo |
| `CARQUET_BUILD_BENCHMARKS` | OFF | Compila os programas de benchmark e profiling |
| `CARQUET_BUILD_ARROW_CPP_BENCHMARK` | OFF | Benchmark opcional de comparação com Arrow C++ |
| `CARQUET_BUILD_INTEROP` | OFF | Compila os testes de interoperabilidade |
| `CARQUET_BUILD_FUZZ` | OFF | Compila os alvos de fuzzing |
| `CARQUET_BUILD_DOCS` | OFF | Adiciona um alvo `docs` (requer Doxygen) |
| `CARQUET_BUNDLE_DEPS` | OFF | Sempre baixa e embute estaticamente zstd/zlib/lz4/brotli (artefatos de release portáveis) |
| `CARQUET_ZSTD_TARGET` / `CARQUET_ZLIB_TARGET` / `CARQUET_LZ4_TARGET` / `CARQUET_BROTLI_TARGET` | *(vazio)* | Alvo(s) CMake existente(s) a linkar para essa biblioteca, em vez de procurá-la ou baixá-la |
| `CARQUET_ARROW_CPP_ROOT` | *(vazio)* | Caminho para uma instalação personalizada do Arrow C++ para o benchmark de comparação |
| `CARQUET_ENABLE_SSE` | ON | Otimizações SSE (x86, detectadas automaticamente) |
| `CARQUET_ENABLE_AVX` | ON | Otimizações AVX (x86, detectadas automaticamente) |
| `CARQUET_ENABLE_AVX2` | ON | Otimizações AVX2 (x86, detectadas automaticamente) |
| `CARQUET_ENABLE_AVX512` | ON | Otimizações AVX-512 (x86, detectadas automaticamente) |
| `CARQUET_ENABLE_NEON` | ON | Otimizações NEON (ARM, detectadas automaticamente) |

</div>

</details>

<details>
<summary><b>Usando seus próprios zstd / zlib / lz4</b> (embutindo o carquet em um projeto CMake maior)</summary>

<br/>

Quando o carquet é incluído com `add_subdirectory()` ou `FetchContent`, ele resolve cada biblioteca de compressão nesta ordem:

1. O alvo nomeado por `CARQUET_ZSTD_TARGET`, `CARQUET_ZLIB_TARGET` ou `CARQUET_LZ4_TARGET`, se definido.
2. Uma cópia do sistema (ignorada quando `CARQUET_BUNDLE_DEPS` está ON).
3. Uma cópia que ele mesmo baixa e compila. Se o seu projeto já compila a biblioteca sob o nome de alvo habitual (`libzstd_static`, `zlibstatic`, `lz4_static`), o carquet reutiliza esse alvo em vez de baixar uma segunda cópia.

Um projeto que vendoriza essas bibliotecas e não tem cópia do sistema, portanto, não precisa de configuração alguma. Para escolher explicitamente, por exemplo sua própria zlib em vez da do sistema, ou seu build compartilhado em vez do estático, nomeie o alvo:

```cmake
set(CARQUET_ZLIB_TARGET zlibstatic)
add_subdirectory(carquet)
```

O mesmo pode ser feito sem editar nada: `cmake -DCARQUET_ZLIB_TARGET=zlib ...`. As três variáveis são independentes. Um alvo nomeado é linkado como está, sem busca e sem download, e um nome que não é um alvo gera um erro de configuração em vez de um fallback silencioso. O alvo pode ser importado ou compilado na árvore, estático ou compartilhado, e deve carregar seus diretórios de include como requisito de uso.

</details>

<details>
<summary><b>Build de desenvolvimento</b> (testes, exemplos, benchmarks)</summary>

<br/>

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug -DCARQUET_BUILD_DEV=ON
cmake --build build -j$(nproc)
cd build && ctest --output-on-failure
```

[CONTRIBUTING.md](CONTRIBUTING.md) descreve a verificação exigida por pull request (testes, sanitizers, fuzzing, interop, benchmarks).

</details>

<details>
<summary><b>Documentação da API</b> (Doxygen)</summary>

<br/>

Os headers públicos são anotados com comentários [Doxygen](https://www.doxygen.nl). Gerar a referência HTML requer o `doxygen` instalado (mais o `graphviz`, opcional, para diagramas de include/dependências); nenhum dos dois é necessário para um build normal.

```bash
# CMake — habilite o alvo na configuração e depois compile-o
cmake -B build -DCARQUET_BUILD_DOCS=ON
cmake --build build --target docs

# xmake — uma task independente, sem flag de configuração
xmake docs
```

Ambos escrevem em `build/docs/html/index.html`. Se o Doxygen não estiver instalado, o alvo/task apenas informa isso e não faz nada.

</details>

<details>
<summary><b>Agentes de IA de programação</b> (Claude Code, Codex: uma skill pronta para uso)</summary>

<br/>

[`.agents/skills/carquet/`](.agents/skills/carquet/) é uma [Agent Skill](https://agentskills.io) que ensina um agente de programação a integrar o carquet ao seu projeto: configuração do build (CMake, FetchContent, pkg-config, xmake), qual API usar para quê, receitas completas, as regras de níveis / null-bitmap / tempo de vida, mapeamentos de tipos e as armadilhas habituais com suas correções. Copie-a ou crie um link para ela no seu próprio repositório e o agente a detecta automaticamente:

```bash
# Claude Code
mkdir -p .claude/skills && cp -r path/to/carquet/.agents/skills/carquet .claude/skills/
# Codex CLI
mkdir -p .agents/skills && cp -r path/to/carquet/.agents/skills/carquet .agents/skills/
# carquet vendorizado: um symlink mantém tudo sincronizado
ln -s ../../third_party/carquet/.agents/skills/carquet .claude/skills/carquet
```

Os mesmos arquivos também funcionam como skill pessoal (`~/.claude/skills/`, `~/.codex/skills/`).

</details>

<br/>

### Compilando com xmake

Um build [xmake](https://xmake.io) (`xmake.lua`) é fornecido como alternativa ao CMake, com as mesmas opções e padrões; zstd/zlib/lz4 são linkados estaticamente, de modo que os binários são autocontidos.

```bash
xmake                    # compila a biblioteca estática + CLI `carquet` (release)
xmake f --dev=y && xmake # adiciona testes, exemplos, benchmarks, interop
xmake test               # executa a suíte de testes
```

<details>
<summary><b>Opções de configuração</b> (<code>xmake f --option=y|n</code>)</summary>

<br/>

As opções espelham as do CMake, sem o prefixo `CARQUET_` e em minúsculas:

<div align="center">

| Opção xmake | Padrão | Descrição |
|--------------|---------|-------------|
| `--dev` | n | Compila testes, exemplos, benchmarks e interop |
| `--tests` / `--examples` / `--benchmarks` / `--interop` | n | Compila um grupo individualmente |
| `--cli` | y | Compila a ferramenta CLI `carquet` |
| `--shared` | n | Compila uma biblioteca compartilhada em vez de estática |
| `--openmp` | y | Leitura paralela de colunas com OpenMP (desativada automaticamente se indisponível) |
| `--native_arch` | n | `-march=native` para desempenho máximo (binário só para a máquina atual) |
| `--sse` / `--avx` / `--avx2` / `--avx512` / `--neon` | y | Conjuntos de instruções SIMD (detectados automaticamente) |
| `--sve` | n | ARM SVE (experimental) |
| `--brotli` | y | Codec BROTLI (libbrotli do repositório do xmake) |
| `--fuzz` | n | Compila os alvos de fuzzing (use `--toolchain=clang`) |

</div>

```bash
xmake f -m release --shared=y            # biblioteca compartilhada
xmake f --dev=y --avx512=n && xmake      # build de desenvolvimento, AVX-512 desativado
```

</details>

<br/>

---

<!-- ────────────────────────────  API C  ──────────────────────────── -->

## API C

Um exemplo de escrita e um de leitura abaixo; o [manual em `docs/`](docs/README.md) cobre todo o resto.

### Escrever um arquivo Parquet

```c
#include <carquet/carquet.h>

int main(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    // Define o schema
    carquet_schema_t* schema = carquet_schema_create(&err);
    carquet_schema_add_column(schema, "id",    CARQUET_PHYSICAL_INT64,  NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(schema, "value", CARQUET_PHYSICAL_DOUBLE, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);

    // Configura o writer
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    opts.compression = CARQUET_COMPRESSION_ZSTD;

    // Escreve
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

    // Abre com mmap para o melhor desempenho de leitura
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.use_mmap = true;

    carquet_reader_t* r = carquet_reader_open("output.parquet", &opts, &err);
    if (!r) { printf("Error: %s\n", err.message); return 1; }

    printf("Rows: %lld, Columns: %d\n",
           (long long)carquet_reader_num_rows(r),
           carquet_reader_num_columns(r));

    // Batch reader para iteração eficiente
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
        // processa ids[0..n-1] ...
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

O [manual](docs/README.md) explica qual superfície usar em cada situação; o layout do código-fonte e a arquitetura estão documentados em [CONTRIBUTING.md](CONTRIBUTING.md). Uma referência HTML navegável pode ser gerada a partir dos headers com Doxygen, veja a seção recolhível *Documentação da API* em [Opções de build](#opções-de-build).

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

./reader                    # lê data.parquet (embutido como padrão)
./reader other.parquet      # sobrescreve com outro arquivo
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

<!-- ────────────────────────────  BENCHMARKS  ──────────────────────────── -->

## Benchmarks

<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="res/img/bench_10m_m3_dark.png">
    <img src="res/img/bench_10m_m3_light.png" alt="Tempo de escrita e leitura de 10 milhões de linhas em Apple M3: carquet vs Arrow C++ vs PyArrow, por codec" width="900">
  </picture>
</p>

Com 10M de linhas (o tamanho mais representativo); razão maior = carquet mais rápido.

ARM (Apple M3): Carquet 0.7.2 vs Arrow C++ 24.0.0, medido em 2026-10-20. x86 (Xeon D-1531): Carquet 0.4.4 vs Arrow C++ 23.0.1; essa execução é anterior ao paralelismo do writer da 0.7.1 e não foi medida novamente desde então.

<div align="center">

| | x86 (Xeon D-1531) | | ARM (Apple M3) | |
|---|---|---|---|---|
| **Codec** | **Escrita** | **Leitura** | **Escrita** | **Leitura** |
| snappy | **1.55x** | **1.25x** | **4.63x** | **1.66x** |
| zstd | **1.31x** | **1.04x** | **5.66x** | **1.36x** |
| lz4 | **1.02x** | 0.83x | **4.47x** | **1.76x** |
| none | **1.13x** | **40.6x**\* | **1.83x** | **38.9x**\* |

</div>

\* As leituras não comprimidas usam mmap zero-copy, veja a nota de metodologia.

As leituras comprimidas envolvem descompressão e decodificação completas de cada valor, sem atalhos, e ambas as bibliotecas usam as mesmas bibliotecas compartilhadas lz4/zstd do sistema, de modo que a velocidade bruta do codec é idêntica. A comparação mais significativa é a tabela de **leitura cruzada do mesmo arquivo** (nos resultados completos), em que ambas as bibliotecas leem exatamente o mesmo arquivo Parquet: o carquet lê dados comprimidos **1,3–2,1x mais rápido** que o Arrow C++ nesse teste em condições idênticas (execução x86: 1,6–2,6x).

<details>
<summary><b>Metodologia dos benchmarks</b></summary>

<br/>

Todos os benchmarks usam dados idênticos (PRNG LCG determinístico), configurações Parquet idênticas (sem dicionário, BYTE_STREAM_SPLIT para floats, checksums de página, leituras via mmap), mediana aparada de 11–51 iterações, com o cache de páginas do SO limpo entre as fases de escrita e leitura e um intervalo de resfriamento entre as configurações. Schema: 3 colunas (INT64, DOUBLE, INT32). Comparado com o leitor Parquet de baixo nível do Arrow C++ (sem materialização em Arrow Table) e com o PyArrow.

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

*Carquet 0.7.2, MacBook Air M3, 16GB de RAM, macOS 26.5, Arrow C++ 24.0.0, PyArrow 23.0.1, ZSTD nível 1. As tabelas de 10M de linhas foram medidas em 2026-10-20; as tabelas de 1M e 100K são da execução da 0.7.1 de 2026-10-08 (não medidas novamente).*

#### 10M de linhas vs Arrow C++

<div align="center">

| Codec | Escrita Carquet | Escrita Arrow C++ | Razão E | Leitura Carquet | Leitura Arrow C++ | Razão L | Tamanho |
|-------|--------------|-----------------|---------|-------------|----------------|---------|------|
| none | **74.7ms** | 137ms | **1.83x** | **0.36ms** | 14.0ms | **38.9x**\* | 200.0MB |
| snappy | **52.7ms** | 244ms | **4.63x** | **12.4ms** | 20.6ms | **1.66x** | 131.2MB |
| zstd | **60.4ms** | 342ms | **5.66x** | **19.8ms** | 26.9ms | **1.36x** | 99.9MB |
| lz4 | **57.2ms** | 256ms | **4.47x** | **9.22ms** | 16.2ms | **1.76x** | 128.9MB |

</div>

<br/>

#### 1M de linhas vs Arrow C++

<div align="center">

| Codec | Escrita Carquet | Escrita Arrow C++ | Razão E | Leitura Carquet | Leitura Arrow C++ | Razão L |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **6.48ms** | 14.82ms | **2.29x** | **0.06ms** | 1.63ms | **27.2x**\* |
| snappy | **14.28ms** | 26.44ms | **1.85x** | **2.14ms** | 2.85ms | **1.33x** |
| zstd | **14.76ms** | 36.17ms | **2.45x** | **2.70ms** | 3.46ms | **1.28x** |
| lz4 | **12.30ms** | 26.10ms | **2.12x** | **1.04ms** | 1.87ms | **1.80x** |

</div>

<br/>

#### 100K linhas vs Arrow C++

<div align="center">

| Codec | Escrita Carquet | Escrita Arrow C++ | Razão E | Leitura Carquet | Leitura Arrow C++ | Razão L |
|-------|--------------|-----------------|---------|-------------|----------------|---------|
| none | **1.07ms** | 1.79ms | **1.67x** | **0.02ms** | 0.25ms | **12.5x**\* |
| snappy | **1.60ms** | 2.63ms | **1.64x** | **0.34ms** | 0.91ms | **2.68x** |
| zstd | **1.73ms** | 3.63ms | **2.10x** | **0.63ms** | 1.26ms | **2.00x** |
| lz4 | **1.52ms** | 2.62ms | **1.72x** | **0.25ms** | 0.55ms | **2.20x** |

</div>

<br/>

#### Leitura cruzada do mesmo arquivo (10M de linhas)

Ambas as bibliotecas leem o **mesmo** arquivo Parquet, a comparação mais justa em condições idênticas.

<div align="center">

| Codec | Escritor | Leitura Carquet | Leitura Arrow C++ | Razão |
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

#### 10M de linhas vs PyArrow

<div align="center">

| Codec | Escrita Carquet | Escrita PyArrow | Razão E | Leitura Carquet | Leitura PyArrow | Razão L |
|-------|--------------|---------------|---------|-------------|--------------|---------|
| none | **74.7ms** | 202ms | **2.71x** | **0.36ms** | 36.2ms | **100.5x**\* |
| snappy | **52.7ms** | 312ms | **5.93x** | **12.4ms** | 46.6ms | **3.75x** |
| zstd | **60.4ms** | 410ms | **6.78x** | **19.8ms** | 75.6ms | **3.82x** |
| lz4 | **57.2ms** | 321ms | **5.62x** | **9.22ms** | 40.0ms | **4.33x** |

</div>

<br/>

#### 1M de linhas vs PyArrow

<div align="center">

| Codec | Escrita Carquet | Escrita PyArrow | Razão E | Leitura Carquet | Leitura PyArrow | Razão L |
|-------|--------------|-----------------|---------|-------------|--------------|---------|
| none | **6.48ms** | 20.17ms | **3.11x** | **0.06ms** | 2.79ms | **46.5x**\* |
| snappy | **14.28ms** | 31.54ms | **2.21x** | **2.14ms** | 4.01ms | **1.87x** |
| zstd | **14.76ms** | 42.51ms | **2.88x** | **2.70ms** | 4.42ms | **1.64x** |
| lz4 | **12.30ms** | 32.03ms | **2.60x** | **1.04ms** | 3.22ms | **3.10x** |

</div>

<br/>

#### 100K linhas vs PyArrow

<div align="center">

| Codec | Escrita Carquet | Escrita PyArrow | Razão E | Leitura Carquet | Leitura PyArrow | Razão L |
|-------|--------------|---------------|---------|-------------|--------------|---------|
| none | **1.07ms** | 2.12ms | **1.98x** | **0.02ms** | 0.24ms | **12.0x**\* |
| snappy | **1.60ms** | 3.17ms | **1.98x** | **0.34ms** | 0.60ms | **1.76x** |
| zstd | **1.73ms** | 4.16ms | **2.40x** | **0.63ms** | 0.81ms | **1.29x** |
| lz4 | **1.52ms** | 3.21ms | **2.11x** | **0.25ms** | 0.42ms | **1.68x** |

</div>

\* Caminho mmap zero-copy

</details>

<details>
<summary><b>Executando os benchmarks você mesmo</b></summary>

<br/>

```bash
# Build com otimizações máximas
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON -DCARQUET_BUILD_DEV=ON
cmake --build build -j$(nproc)

cd build
./benchmark_carquet                     # Carquet isoladamente
python3 ../benchmark/run_benchmark.py   # Comparação completa (+ PyArrow, + Arrow C++)

# Pula as configurações de 100M de linhas (xlarge) — elas escrevem arquivos de ~2GB por codec
# e podem levar mais de 30 minutos dependendo do hardware
python3 ../benchmark/run_benchmark.py --skip-xlarge

# Sobrescreve o nível do ZSTD (padrão: 1)
CARQUET_BENCH_ZSTD_LEVEL=3 python3 ../benchmark/run_benchmark.py
```

Para incluir a comparação com Arrow C++, habilite o alvo correspondente na configuração:

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release -DCARQUET_NATIVE_ARCH=ON \
  -DCARQUET_BUILD_BENCHMARKS=ON \
  -DCARQUET_BUILD_ARROW_CPP_BENCHMARK=ON
cmake --build build -j$(nproc)

# Ou aponte para uma instalação personalizada do Arrow
cmake -B build ... -DCARQUET_ARROW_CPP_ROOT=/path/to/arrow-prefix
```

O benchmark do Arrow C++ usa a API de baixo nível `parquet::ParquetFileReader` (evitando o custo de materialização em Arrow Table) com leitores de row group paralelos. O modo de **leitura cruzada do mesmo arquivo** faz ambas as bibliotecas lerem exatamente o mesmo arquivo Parquet, eliminando diferenças de tamanho de página, encoding e layout de row groups. Ambos os benchmarks usam dados idênticos, mesmo dimensionamento de row groups, sem dicionário, checksums de página, leituras via mmap, BYTE_STREAM_SPLIT para floats.

</details>

<br/>

---

<!-- ────────────────────────────  RECURSOS SUPORTADOS  ──────────────────────────── -->

## Suporte a recursos do Parquet

O que a biblioteca implementa hoje, por área.

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
table = pq.read_table("carquet_output.parquet")  # simplesmente funciona
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
