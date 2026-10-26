/**
 * @file codegen_read.c
 * @brief Code generation: reads a parquet file's schema and generates
 *        type-correct C source code for reading files with that schema.
 */

#include "cli.h"
#include "reader/reader_internal.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <inttypes.h>
#ifdef _WIN32
#include <io.h>
#include <direct.h>
#include <stdlib.h>      /* _MAX_PATH, _fullpath */
#define codegen_getcwd _getcwd
#define CODEGEN_PATH_MAX _MAX_PATH
#else
#include <unistd.h>
#include <limits.h>      /* PATH_MAX */
#define codegen_getcwd getcwd
#ifdef PATH_MAX
#define CODEGEN_PATH_MAX PATH_MAX
#else
#define CODEGEN_PATH_MAX 4096
#endif
#endif

/* ── Helpers ──────────────────────────────────────────────────────────── */

static void sanitize_ident(const char* name, char* out, size_t out_size) {
    size_t j = 0;
    for (size_t i = 0; name[i] && j < out_size - 1; i++) {
        char ch = name[i];
        if (isalnum((unsigned char)ch) || ch == '_')
            out[j++] = ch;
        else if (ch == '.' || ch == '-' || ch == ' ')
            out[j++] = '_';
    }
    if (j == 0 && out_size > 1) out[j++] = '_';
    out[j] = '\0';
    if (isdigit((unsigned char)out[0]) && j + 1 < out_size) {
        memmove(out + 1, out, j + 1);
        out[0] = '_';
    }
}

static const char* c_type_for(carquet_physical_type_t phys) {
    switch (phys) {
        case CARQUET_PHYSICAL_BOOLEAN:    return "uint8_t";
        case CARQUET_PHYSICAL_INT32:      return "int32_t";
        case CARQUET_PHYSICAL_INT64:      return "int64_t";
        case CARQUET_PHYSICAL_FLOAT:      return "float";
        case CARQUET_PHYSICAL_DOUBLE:     return "double";
        case CARQUET_PHYSICAL_BYTE_ARRAY: return "carquet_byte_array_t";
        case CARQUET_PHYSICAL_INT96:      return "carquet_int96_t";
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY: return "uint8_t";
        default:                          return "uint8_t";
    }
}

static void type_comment(carquet_physical_type_t phys,
                         const carquet_logical_type_t* lt,
                         carquet_field_repetition_t rep,
                         char* buf, size_t buf_size) {
    char type_str[64];
    cli_format_type(phys, lt, type_str, sizeof(type_str));
    snprintf(buf, buf_size, "%s, %s", type_str, cli_repetition_name(rep));
}

static bool column_matches_filter(const char* name, const char* filter) {
    if (!filter) return true;
    const char* p = filter;
    size_t name_len = strlen(name);
    while (*p) {
        const char* comma = strchr(p, ',');
        size_t tok_len = comma ? (size_t)(comma - p) : strlen(p);
        while (tok_len > 0 && p[tok_len - 1] == ' ') tok_len--;
        const char* start = p;
        while (*start == ' ' && tok_len > 0) { start++; tok_len--; }
        if (tok_len == name_len && strncmp(start, name, tok_len) == 0)
            return true;
        p = comma ? comma + 1 : p + strlen(p);
    }
    return false;
}

/* ── Line-counting fprintf wrapper ────────────────────────────────────── */

static int g_line;  /* current output line number (1-based) */

static void emit(FILE* out, const char* fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    char buf[4096];
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    /* Count newlines */
    for (const char* p = buf; *p; p++) {
        if (*p == '\n') g_line++;
    }
    fputs(buf, out);
}

/* ── Build instruction detection ──────────────────────────────────────── */

static int file_exists(const char* path) {
#ifdef _WIN32
    return _access(path, 0) == 0;
#else
    return access(path, F_OK) == 0;
#endif
}

static const char* detect_compiler(void) {
    const char* cc_env = getenv("CC");
    if (cc_env && cc_env[0]) return cc_env;
    static const char* candidates[] = {
#ifdef _WIN32
        "cl",
#elif defined(__APPLE__)
        "/opt/homebrew/opt/llvm/bin/clang",
        "/usr/local/opt/llvm/bin/clang",
#endif
        NULL
    };
    for (int i = 0; candidates[i]; i++) {
        if (file_exists(candidates[i])) return candidates[i];
    }
#ifdef _WIN32
    return "cl";
#else
    return "cc";
#endif
}

static void derive_binary_name(const char* output_path, char* buf, size_t buf_size) {
    if (!output_path) { snprintf(buf, buf_size, "reader"); return; }
    const char* slash = strrchr(output_path, '/');
#ifdef _WIN32
    const char* bs = strrchr(output_path, '\\');
    if (bs && (!slash || bs > slash)) slash = bs;
#endif
    const char* base = slash ? slash + 1 : output_path;
    snprintf(buf, buf_size, "%s", base);
    size_t len = strlen(buf);
    if (len > 2 && buf[len - 2] == '.' && buf[len - 1] == 'c')
        buf[len - 2] = '\0';
}

static void detect_link_deps(char* buf, size_t buf_size) {
    buf[0] = '\0';
    size_t off = 0;
    static const char* flags[] = { "-lzstd", "-lz", "-llz4", "-lm", NULL };
    for (int i = 0; flags[i]; i++) {
        int n = snprintf(buf + off, buf_size - off, " %s", flags[i]);
        if (n > 0) off += (size_t)n;
    }
#ifndef _WIN32
    { int n = snprintf(buf + off, buf_size - off, " -lpthread"); if (n > 0) off += (size_t)n; }
#endif
#ifdef __APPLE__
    { int n = snprintf(buf + off, buf_size - off, " -Wl,-w"); if (n > 0) off += (size_t)n; }
#endif
}

static void detect_build_line(const codegen_opts_t* opts, char* out, size_t out_size) {
    char binary_name[256];
    derive_binary_name(opts->output_path, binary_name, sizeof(binary_name));
    const char* source_name = opts->output_path ? opts->output_path : "reader.c";
    const char* compiler = detect_compiler();
    char deps[256];
    detect_link_deps(deps, sizeof(deps));
    char probe[1024], cwd[512];

    /* 1. Local repo */
    if (codegen_getcwd(cwd, sizeof(cwd))) {
        snprintf(probe, sizeof(probe), "%s/include/carquet/carquet.h", cwd);
        if (file_exists(probe)) {
            char lib_probe[1024];
            snprintf(lib_probe, sizeof(lib_probe), "%s/build/libcarquet.a", cwd);
            if (file_exists(lib_probe)) {
                snprintf(out, out_size, "%s -o %s %s -I%s/include -L%s/build -lcarquet%s",
                         compiler, binary_name, source_name, cwd, cwd, deps);
                return;
            }
            snprintf(out, out_size, "%s -o %s %s -I%s/include -L/path/to/lib -lcarquet%s",
                     compiler, binary_name, source_name, cwd, deps);
            return;
        }
        char parent[512];
        snprintf(parent, sizeof(parent), "%s/..", cwd);
        snprintf(probe, sizeof(probe), "%s/include/carquet/carquet.h", parent);
        if (file_exists(probe)) {
            snprintf(out, out_size, "%s -o %s %s -I%s/include -L%s -lcarquet%s",
                     compiler, binary_name, source_name, parent, cwd, deps);
            return;
        }
    }

    /* 2. System-wide */
    static const char* sys[] = {
        "/usr/local/include/carquet/carquet.h",
        "/usr/include/carquet/carquet.h",
        "/opt/homebrew/include/carquet/carquet.h",
        NULL
    };
    for (int i = 0; sys[i]; i++) {
        if (file_exists(sys[i])) {
            snprintf(out, out_size, "%s -o %s %s -lcarquet", compiler, binary_name, source_name);
            return;
        }
    }

    /* 3. Fallback */
    snprintf(out, out_size, "%s -o %s %s -I/path/to/include -L/path/to/lib -lcarquet%s",
             compiler, binary_name, source_name, deps);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Code generation: --read
 * ══════════════════════════════════════════════════════════════════════════ */

int cmd_codegen_read(FILE* out, carquet_reader_t* reader,
                     const codegen_opts_t* opts,
                     codegen_hints_t* hints) {
    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    int64_t total_rows = carquet_reader_num_rows(reader);
    int32_t batch_size = opts->batch_size;

    bool* include = calloc((size_t)num_cols, sizeof(bool));
    int32_t included_count = 0;
    for (int32_t c = 0; c < num_cols; c++) {
        if (column_matches_filter(carquet_schema_column_name(schema, c), opts->columns)) {
            include[c] = true;
            included_count++;
        }
    }
    if (included_count == 0) {
        fprintf(stderr, "error: no columns match the filter\n");
        free(include);
        return 1;
    }

    detect_build_line(opts, hints->build_line, sizeof(hints->build_line));
    g_line = 1;

    /* ── Header ───────────────────────────────────────────────────── */
    emit(out,
        "/*\n"
        " * Auto-generated by: carquet codegen --read\n"
        " * Source file:        %s\n"
        " * Schema:             %" PRId64 " rows, %d columns\n"
        " *\n"
        " * Build:\n"
        " *   %s\n"
        " */\n\n"
        "#include <stdio.h>\n"
        "#include <stdlib.h>\n"
        "#include <inttypes.h>\n"
        "#include <carquet/carquet.h>\n\n",
        opts->input_path ? opts->input_path : "<file.parquet>",
        total_rows, num_cols, hints->build_line);

    /* ── Schema documentation ─────────────────────────────────────── */
    emit(out, "/*\n * Schema:\n");
    for (int32_t c = 0; c < num_cols; c++) {
        if (!include[c]) continue;
        const carquet_schema_node_t* node = carquet_schema_get_element(schema,
            schema->leaf_indices[c]);
        char desc[128];
        type_comment(carquet_schema_column_type(schema, c),
                    carquet_schema_node_logical_type(node),
                    carquet_schema_node_repetition(node), desc, sizeof(desc));
        emit(out, " *   [%d] %-30s %s\n", c,
             carquet_schema_column_name(schema, c), desc);
    }
    emit(out, " */\n\n");

    /* ── Process callback ──────────────────────────────────────────── */
    emit(out,
        "/* Called once per batch of rows read from each row group. */\n"
        "static void process_batch(\n"
        "    int32_t row_group,\n"
        "    int64_t batch_offset,\n"
        "    int64_t count,\n");

    for (int32_t c = 0; c < num_cols; c++) {
        if (!include[c]) continue;
        const carquet_schema_node_t* node = carquet_schema_get_element(schema,
            schema->leaf_indices[c]);
        carquet_physical_type_t phys = carquet_schema_column_type(schema, c);
        bool nullable = carquet_schema_node_repetition(node) != CARQUET_REPETITION_REQUIRED;
        int32_t tl = carquet_schema_node_type_length(node);
        const char* ctype = c_type_for(phys);
        char ident[128];
        sanitize_ident(carquet_schema_column_name(schema, c), ident, sizeof(ident));

        if (phys == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY)
            emit(out, "    const uint8_t* %s,  /* [count * %d] */\n", ident, tl);
        else
            emit(out, "    const %s* %s,\n", ctype, ident);
        if (nullable)
            emit(out, "    const int16_t* %s_def,\n", ident);
    }

    emit(out, "    int dummy)\n{\n");

    /* Record the line where process_batch body starts */
    hints->process_batch_line = g_line;

    if (opts->skeleton) {
        /* Empty body for user to fill in */
        emit(out,
            "    /* TODO: implement your processing logic here */\n"
            "    (void)row_group; (void)batch_offset; (void)count; (void)dummy;\n");
        /* Suppress unused-parameter warnings for all column args */
        for (int32_t c = 0; c < num_cols; c++) {
            if (!include[c]) continue;
            const carquet_schema_node_t* node = carquet_schema_get_element(schema,
                schema->leaf_indices[c]);
            bool nullable = carquet_schema_node_repetition(node) != CARQUET_REPETITION_REQUIRED;
            char ident[128];
            sanitize_ident(carquet_schema_column_name(schema, c), ident, sizeof(ident));
            if (nullable)
                emit(out, "    (void)%s; (void)%s_def;\n", ident, ident);
            else
                emit(out, "    (void)%s;\n", ident);
        }
    } else {
        /* Default: print each row as tab-separated values */
        emit(out,
            "    (void)row_group; (void)batch_offset; (void)dummy;\n");

        for (int32_t c = 0; c < num_cols; c++) {
            if (!include[c]) continue;
            const carquet_schema_node_t* node = carquet_schema_get_element(schema,
                schema->leaf_indices[c]);
            bool nullable = carquet_schema_node_repetition(node) != CARQUET_REPETITION_REQUIRED;
            if (nullable) {
                char ident[128];
                sanitize_ident(carquet_schema_column_name(schema, c), ident, sizeof(ident));
                emit(out, "    int64_t %s_value_index = 0;\n", ident);
            }
        }

        emit(out, "    for (int64_t i = 0; i < count; i++) {\n");

        int col_printed = 0;
        for (int32_t c = 0; c < num_cols; c++) {
            if (!include[c]) continue;
            const carquet_schema_node_t* node = carquet_schema_get_element(schema,
                schema->leaf_indices[c]);
            carquet_physical_type_t phys = carquet_schema_column_type(schema, c);
            bool nullable = carquet_schema_node_repetition(node) != CARQUET_REPETITION_REQUIRED;
            int32_t tl = carquet_schema_node_type_length(node);
            char ident[128];
            sanitize_ident(carquet_schema_column_name(schema, c), ident, sizeof(ident));
            const char* sep = col_printed > 0 ? "\\t" : "";

            if (nullable) {
                int16_t max_def = carquet_schema_node_max_def_level(node);
                emit(out,
                    "        if (%s_def[i] < %d) {\n"
                    "            printf(\"%s\");\n"
                    "        } else {\n", ident, max_def, sep);
                /* Print value (indented inside else) */
                switch (phys) {
                    case CARQUET_PHYSICAL_BOOLEAN:
                        emit(out, "            printf(\"%s%%s\", %s[%s_value_index] ? \"true\" : \"false\");\n", sep, ident, ident); break;
                    case CARQUET_PHYSICAL_INT32:
                        emit(out, "            printf(\"%s%%\" PRId32, %s[%s_value_index]);\n", sep, ident, ident); break;
                    case CARQUET_PHYSICAL_INT64:
                        emit(out, "            printf(\"%s%%\" PRId64, %s[%s_value_index]);\n", sep, ident, ident); break;
                    case CARQUET_PHYSICAL_FLOAT:
                        emit(out, "            printf(\"%s%%g\", (double)%s[%s_value_index]);\n", sep, ident, ident); break;
                    case CARQUET_PHYSICAL_DOUBLE:
                        emit(out, "            printf(\"%s%%g\", %s[%s_value_index]);\n", sep, ident, ident); break;
                    case CARQUET_PHYSICAL_BYTE_ARRAY:
                        emit(out, "            printf(\"%s%%.*s\", %s[%s_value_index].length, (const char*)%s[%s_value_index].data);\n", sep, ident, ident, ident, ident); break;
                    case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
                        emit(out, "            for (int32_t b = 0; b < %d; b++) printf(\"%s%%02x\", %s[%s_value_index * %d + b]);\n", tl, sep, ident, ident, tl); break;
                    default:
                        emit(out, "            printf(\"%s?\");\n", sep); break;
                }
                emit(out, "            %s_value_index++;\n", ident);
                emit(out, "        }\n");
            } else {
                switch (phys) {
                    case CARQUET_PHYSICAL_BOOLEAN:
                        emit(out, "        printf(\"%s%%s\", %s[i] ? \"true\" : \"false\");\n", sep, ident); break;
                    case CARQUET_PHYSICAL_INT32:
                        emit(out, "        printf(\"%s%%\" PRId32, %s[i]);\n", sep, ident); break;
                    case CARQUET_PHYSICAL_INT64:
                        emit(out, "        printf(\"%s%%\" PRId64, %s[i]);\n", sep, ident); break;
                    case CARQUET_PHYSICAL_FLOAT:
                        emit(out, "        printf(\"%s%%g\", (double)%s[i]);\n", sep, ident); break;
                    case CARQUET_PHYSICAL_DOUBLE:
                        emit(out, "        printf(\"%s%%g\", %s[i]);\n", sep, ident); break;
                    case CARQUET_PHYSICAL_BYTE_ARRAY:
                        emit(out, "        printf(\"%s%%.*s\", %s[i].length, (const char*)%s[i].data);\n", sep, ident, ident); break;
                    case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
                        emit(out, "        for (int32_t b = 0; b < %d; b++) printf(\"%s%%02x\", %s[i * %d + b]);\n", tl, sep, ident, tl); break;
                    default:
                        emit(out, "        printf(\"%s?\");\n", sep); break;
                }
            }
            col_printed++;
        }
        emit(out,
            "        printf(\"\\n\");\n"
            "    }\n");
    }
    emit(out, "}\n\n");

    /* ── main() ───────────────────────────────────────────────────── */

    /* DEFAULT_FILE: resolve to absolute, or use placeholder.
     * Sized for the escaped form of a full PATH_MAX path (each byte may double). */
    char escaped_path[CODEGEN_PATH_MAX * 2];
    if (opts->input_path) {
        /* realpath() requires a buffer of at least PATH_MAX bytes regardless of
         * the actual path length; an undersized buffer trips glibc's
         * _FORTIFY_SOURCE check (__realpath_chk -> "buffer overflow detected")
         * in optimized builds where fortification is active. */
        char abs_input[CODEGEN_PATH_MAX];
#ifdef _WIN32
        if (!_fullpath(abs_input, opts->input_path, sizeof(abs_input)))
            snprintf(abs_input, sizeof(abs_input), "%s", opts->input_path);
#else
        /* Pass NULL so realpath() allocates a buffer of the system's actual
         * PATH_MAX itself; writing into a fixed CODEGEN_PATH_MAX stack buffer
         * overflows when the runtime PATH_MAX exceeds the compile-time fallback. */
        char* resolved = realpath(opts->input_path, NULL);
        if (resolved) {
            snprintf(abs_input, sizeof(abs_input), "%s", resolved);
            free(resolved);
        } else {
            snprintf(abs_input, sizeof(abs_input), "%s", opts->input_path);
        }
#endif
        size_t j = 0;
        for (size_t i = 0; abs_input[i] && j < sizeof(escaped_path) - 2; i++) {
            if (abs_input[i] == '\\') escaped_path[j++] = '\\';
            escaped_path[j++] = abs_input[i];
        }
        escaped_path[j] = '\0';
    } else {
        snprintf(escaped_path, sizeof(escaped_path), "/path/to/file.parquet");
    }

    hints->default_file_line = g_line + 1; /* next line emitted is #define */
    emit(out, "#define DEFAULT_FILE \"%s\"\n\n", escaped_path);

    /* reader_open — with mmap option if requested */
    if (opts->use_mmap) {
        emit(out,
            "int main(int argc, char** argv) {\n"
            "    const char* path = (argc >= 2) ? argv[1] : DEFAULT_FILE;\n\n"
            "    carquet_error_t err = CARQUET_ERROR_INIT;\n"
            "    carquet_reader_options_t ropts;\n"
            "    carquet_reader_options_init(&ropts);\n"
            "    ropts.use_mmap = true;\n\n"
            "    carquet_reader_t* reader = carquet_reader_open(path, &ropts, &err);\n"
            "    if (!reader) {\n"
            "        fprintf(stderr, \"Error opening file: %%s\\n\", err.message);\n"
            "        return 1;\n"
            "    }\n\n");
    } else {
        emit(out,
            "int main(int argc, char** argv) {\n"
            "    const char* path = (argc >= 2) ? argv[1] : DEFAULT_FILE;\n\n"
            "    carquet_error_t err = CARQUET_ERROR_INIT;\n"
            "    carquet_reader_t* reader = carquet_reader_open(path, NULL, &err);\n"
            "    if (!reader) {\n"
            "        fprintf(stderr, \"Error opening file: %%s\\n\", err.message);\n"
            "        return 1;\n"
            "    }\n\n");
    }

    emit(out,
        "    int64_t total_rows = carquet_reader_num_rows(reader);\n"
        "    int32_t num_row_groups = carquet_reader_num_row_groups(reader);\n"
        "    printf(\"Reading %%\" PRId64 \" rows from %%d row groups\\n\",\n"
        "           total_rows, num_row_groups);\n\n");

    emit(out,
        "    if (carquet_reader_num_columns(reader) != %d) {\n"
        "        fprintf(stderr, \"Error: expected %d columns, got %%d\\n\",\n"
        "                carquet_reader_num_columns(reader));\n"
        "        carquet_reader_close(reader);\n"
        "        return 1;\n"
        "    }\n\n", num_cols, num_cols);

    /* ── Row group loop ───────────────────────────────────────────── */
    emit(out, "    for (int32_t rg = 0; rg < num_row_groups; rg++) {\n");

    for (int32_t c = 0; c < num_cols; c++) {
        if (!include[c]) continue;
        const carquet_schema_node_t* node = carquet_schema_get_element(schema, schema->leaf_indices[c]);
        carquet_physical_type_t phys = carquet_schema_column_type(schema, c);
        const carquet_logical_type_t* lt = carquet_schema_node_logical_type(node);
        bool nullable = carquet_schema_node_repetition(node) != CARQUET_REPETITION_REQUIRED;
        int32_t tl = carquet_schema_node_type_length(node);
        const char* ctype = c_type_for(phys);
        char ident[128], desc[128];
        sanitize_ident(carquet_schema_column_name(schema, c), ident, sizeof(ident));
        type_comment(phys, lt, carquet_schema_node_repetition(node), desc, sizeof(desc));

        emit(out, "\n        /* Column %d: %s (%s) */\n", c,
             carquet_schema_column_name(schema, c), desc);
        if (phys == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY)
            emit(out, "        uint8_t %s_buf[%d * %d];\n", ident, batch_size, tl);
        else
            emit(out, "        %s %s_buf[%d];\n", ctype, ident, batch_size);
        if (nullable)
            emit(out, "        int16_t %s_def[%d];\n", ident, batch_size);
    }
    emit(out, "\n");

    /* Open column readers */
    for (int32_t c = 0; c < num_cols; c++) {
        if (!include[c]) continue;
        char ident[128];
        sanitize_ident(carquet_schema_column_name(schema, c), ident, sizeof(ident));
        emit(out,
            "        carquet_column_reader_t* col_%s =\n"
            "            carquet_reader_get_column(reader, rg, %d, &err);\n"
            "        if (!col_%s) {\n"
            "            fprintf(stderr, \"Error reading column %d: %%s\\n\", err.message);\n"
            "            carquet_reader_close(reader);\n"
            "            return 1;\n"
            "        }\n\n", ident, c, ident, c);
    }

    /* Batch read loop */
    emit(out,
        "        int64_t batch_offset = 0;\n"
        "        int done = 0;\n"
        "        while (!done) {\n"
        "            int64_t count = 0;\n\n");

    bool first_col = true;
    for (int32_t c = 0; c < num_cols; c++) {
        if (!include[c]) continue;
        const carquet_schema_node_t* node = carquet_schema_get_element(schema, schema->leaf_indices[c]);
        bool nullable = carquet_schema_node_repetition(node) != CARQUET_REPETITION_REQUIRED;
        char ident[128], def_arg[140];
        sanitize_ident(carquet_schema_column_name(schema, c), ident, sizeof(ident));
        snprintf(def_arg, sizeof(def_arg), nullable ? "%s_def" : "NULL", ident);

        if (first_col) {
            emit(out,
                "            count = carquet_column_read_batch(\n"
                "                col_%s, %s_buf, %d, %s, NULL);\n"
                "            if (count <= 0) { done = 1; break; }\n\n",
                ident, ident, batch_size, def_arg);
            first_col = false;
        } else {
            emit(out,
                "            (void)carquet_column_read_batch(\n"
                "                col_%s, %s_buf, count, %s, NULL);\n\n",
                ident, ident, def_arg);
        }
    }

    /* Call process_batch */
    emit(out, "            process_batch(rg, batch_offset, count,\n");
    for (int32_t c = 0; c < num_cols; c++) {
        if (!include[c]) continue;
        const carquet_schema_node_t* node = carquet_schema_get_element(schema, schema->leaf_indices[c]);
        bool nullable = carquet_schema_node_repetition(node) != CARQUET_REPETITION_REQUIRED;
        char ident[128];
        sanitize_ident(carquet_schema_column_name(schema, c), ident, sizeof(ident));
        emit(out, "                %s_buf,\n", ident);
        if (nullable) emit(out, "                %s_def,\n", ident);
    }
    emit(out,
        "                0);\n"
        "            batch_offset += count;\n"
        "        }\n\n");

    /* Free column readers */
    for (int32_t c = 0; c < num_cols; c++) {
        if (!include[c]) continue;
        char ident[128];
        sanitize_ident(carquet_schema_column_name(schema, c), ident, sizeof(ident));
        emit(out, "        carquet_column_reader_free(col_%s);\n", ident);
    }

    emit(out,
        "    }\n\n"
        "    carquet_reader_close(reader);\n"
        "    printf(\"Done.\\n\");\n"
        "    return 0;\n"
        "}\n");

    free(include);
    return 0;
}
