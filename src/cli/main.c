/**
 * @file main.c
 * @brief Entry point for the carquet CLI tool
 */

#include "cli.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── Help text ────────────────────────────────────────────────────────── */

static void print_usage(void) {
    fprintf(stderr,
        "carquet %s - Parquet file inspector and code generator\n"
        "\n"
        "Usage: carquet <command> [options] <file.parquet>\n"
        "\n"
        "Commands:\n"
        "  schema     Print file schema\n"
        "  info       Print detailed file metadata\n"
        "  head       Print first N rows\n"
        "  tail       Print last N rows\n"
        "  cat        Print rows with optional slicing and column filter\n"
        "  count      Print total row count\n"
        "  columns    List column names (one per line)\n"
        "  stat       Print column statistics\n"
        "  validate   Verify file integrity\n"
        "  sample     Print N random rows\n"
        "  export     Write rows in another format (csv)\n"
        "  codegen    Generate C reader code\n"
        "\n"
        "Run 'carquet <command> -h' for command-specific help.\n",
        CARQUET_VERSION_STRING);
}

static void print_help_schema(void) {
    fprintf(stderr,
        "Usage: carquet schema <file.parquet>\n"
        "\n"
        "Print the Parquet schema in a human-readable tree format.\n"
        "Shows physical types, logical types, and repetition levels.\n"
        "\n"
        "Arguments:\n"
        "  <file.parquet>    Input Parquet file\n");
}

static void print_help_info(void) {
    fprintf(stderr,
        "Usage: carquet info <file.parquet>\n"
        "\n"
        "Print detailed file metadata.\n"
        "\n"
        "Arguments:\n"
        "  <file.parquet>    Input Parquet file\n"
        "\n"
        "Output includes:\n"
        "  - File path, creator, row/column/row-group counts\n"
        "  - Key-value metadata\n"
        "  - Per-column type information and nullability\n"
        "  - Per-row-group size and compression ratio\n");
}

static void print_help_head(void) {
    fprintf(stderr,
        "Usage: carquet head [-n NUM] [-p EXPR] <file.parquet>\n"
        "\n"
        "Print the first N rows in a tabular format. With --filter, prints the\n"
        "first N rows matching the predicate; pages that cannot match are\n"
        "skipped without decompression (requires a file written with\n"
        "write_page_index = true).\n"
        "\n"
        "Arguments:\n"
        "  <file.parquet>    Input Parquet file\n"
        "\n"
        "Options:\n"
        "  -n NUM            Number of rows to display (default: %d)\n"
        "  -p, --filter EXPR Filter expression (see `carquet cat -h`)\n",
        CLI_DEFAULT_NUM_ROWS);
}

static void print_help_tail(void) {
    fprintf(stderr,
        "Usage: carquet tail [-n NUM] <file.parquet>\n"
        "\n"
        "Print the last N rows in a tabular format.\n"
        "\n"
        "Arguments:\n"
        "  <file.parquet>    Input Parquet file\n"
        "\n"
        "Options:\n"
        "  -n NUM            Number of rows to display (default: %d)\n",
        CLI_DEFAULT_NUM_ROWS);
}

static void print_help_count(void) {
    fprintf(stderr,
        "Usage: carquet count [-p EXPR] <file.parquet>\n"
        "\n"
        "Print the total number of rows. Output is a single integer,\n"
        "suitable for use in shell scripts. With --filter, prints the number\n"
        "of rows that match the predicate (page-level pruning skips work).\n"
        "\n"
        "Arguments:\n"
        "  <file.parquet>    Input Parquet file\n"
        "\n"
        "Options:\n"
        "  -p, --filter EXPR Filter expression (see `carquet cat -h`)\n");
}

static void print_help_columns(void) {
    fprintf(stderr,
        "Usage: carquet columns <file.parquet>\n"
        "\n"
        "List column names, one per line. Useful for scripting:\n"
        "  carquet columns data.parquet | grep timestamp\n"
        "\n"
        "Arguments:\n"
        "  <file.parquet>    Input Parquet file\n");
}

static void print_help_stat(void) {
    fprintf(stderr,
        "Usage: carquet stat <file.parquet>\n"
        "\n"
        "Print column statistics (min, max, null count) per row group.\n"
        "Shows '-' when statistics are not available.\n"
        "\n"
        "Arguments:\n"
        "  <file.parquet>    Input Parquet file\n");
}

static void print_help_validate(void) {
    fprintf(stderr,
        "Usage: carquet validate <file.parquet>\n"
        "\n"
        "Verify file integrity by reading all pages with CRC32 checksum\n"
        "verification. Reports OK or lists page read errors.\n"
        "\n"
        "Arguments:\n"
        "  <file.parquet>    Input Parquet file\n");
}

static void print_help_sample(void) {
    fprintf(stderr,
        "Usage: carquet sample [-n NUM] <file.parquet>\n"
        "\n"
        "Print N random rows in a tabular format.\n"
        "\n"
        "Arguments:\n"
        "  <file.parquet>    Input Parquet file\n"
        "\n"
        "Options:\n"
        "  -n NUM            Number of rows to sample (default: %d)\n",
        CLI_DEFAULT_NUM_ROWS);
}

static void print_help_cat(void) {
    fprintf(stderr,
        "Usage: carquet cat [options] <file.parquet>\n"
        "\n"
        "Print rows in a tabular format with optional slicing, column\n"
        "projection, and row-predicate filtering. Unlike head/tail, supports\n"
        "arbitrary row offsets.\n"
        "\n"
        "Arguments:\n"
        "  <file.parquet>    Input Parquet file\n"
        "\n"
        "Options:\n"
        "  -n, --limit N     Number of rows to print (default: all)\n"
        "  -s, --offset N    Skip the first N rows (default: 0)\n"
        "  -c, --columns L   Comma-separated column names (default: all)\n"
        "  -p, --filter EXPR Filter expression (see below)\n"
        "\n"
        "Filter expression grammar (case-insensitive keywords):\n"
        "  expr   := clause (AND clause)*\n"
        "  clause := column OP value\n"
        "          | column IS NULL\n"
        "          | column IS NOT NULL\n"
        "  OP     := = | == | != | <> | < | <= | > | >=\n"
        "  value  := signed_number | 'single-quoted string' | TRUE | FALSE\n"
        "\n"
        "Page-level filtering: pages whose column-index min/max prove no\n"
        "value can match the predicate are skipped without decompression.\n"
        "The file must have been written with write_page_index = true for\n"
        "every column referenced by the filter. INT96 columns are rejected\n"
        "(no defined sort order). Page granularity means rows inside a\n"
        "matching page that fail the predicate are still returned; pipe\n"
        "through awk/grep for exact post-filtering.\n"
        "\n"
        "Examples:\n"
        "  carquet cat -n 1000 data.parquet\n"
        "  carquet cat -c id,name --offset 5000 -n 100 data.parquet\n"
        "  carquet cat -p 'age >= 30 AND status = \\'active\\'' data.parquet\n"
        "  carquet cat -p 'ts >= 1700000000 AND ts < 1700001000' -c ts,event log.parquet\n");
}

static void print_help_export(void) {
    fprintf(stderr,
        "Usage: carquet export [options] <file.parquet>\n"
        "\n"
        "Write rows to stdout in another format. Currently supports CSV.\n"
        "Output is RFC 4180 quoted (header row + comma-separated values).\n"
        "\n"
        "Arguments:\n"
        "  <file.parquet>    Input Parquet file\n"
        "\n"
        "Options:\n"
        "  --format FMT      Output format: csv (default)\n"
        "  -n, --limit N     Number of rows to export (default: all)\n"
        "  -s, --offset N    Skip the first N rows (default: 0)\n"
        "  -c, --columns L   Comma-separated column names (default: all)\n"
        "  -p, --filter EXPR Filter expression (see `carquet cat -h`)\n"
        "\n"
        "Examples:\n"
        "  carquet export data.parquet > data.csv\n"
        "  carquet export -c id,name -n 1000 data.parquet | head\n"
        "  carquet export -p 'price > 100' data.parquet > expensive.csv\n");
}

static void print_help_codegen(void) {
    fprintf(stderr,
        "Usage: carquet codegen [options]\n"
        "\n"
        "Generate type-correct C source code for reading a Parquet file.\n"
        "Inspects the schema of a real file and emits a complete, compilable\n"
        "C program tailored to that schema.\n"
        "\n"
        "Mode:\n"
        "  -r, --read           Generate reader code (default)\n"
        "  -w, --write          Generate writer code (not yet implemented)\n"
        "\n"
        "Options:\n"
        "  -f, --file FILE      Parquet file to inspect schema from\n"
        "  -o, --output FILE    Output source file (default: stdout)\n"
        "  -b, --batch-size N   Batch size in generated code (default: %d)\n"
        "  -c, --columns COLS   Comma-separated column filter\n"
        "  --mmap               Use memory-mapped I/O in generated code\n"
        "  --skeleton           Generate empty process_batch for custom logic\n"
        "\n"
        "Examples:\n"
        "  carquet codegen -r -f data.parquet -o reader.c\n"
        "  carquet codegen -f data.parquet --mmap --skeleton -o reader.c\n"
        "  carquet codegen -f data.parquet -c id,name -o reader.c\n",
        CLI_DEFAULT_BATCH_SIZE);
}

/* ── Argument helpers ─────────────────────────────────────────────────── */

static int is_help_flag(const char* arg) {
    return strcmp(arg, "-h") == 0 || strcmp(arg, "--help") == 0;
}

static int parse_int64(const char* str, int64_t* out) {
    char* end;
    long long val = strtoll(str, &end, 10);
    if (*end != '\0' || end == str || val < 0) return -1;
    *out = (int64_t)val;
    return 0;
}

/* ── main ─────────────────────────────────────────────────────────────── */

int main(int argc, char** argv) {
    if (argc < 2) {
        print_usage();
        return 1;
    }

    const char* cmd = argv[1];

    /* Top-level help */
    if (is_help_flag(cmd) || strcmp(cmd, "help") == 0) {
        print_usage();
        return 0;
    }

    /* ── cat / export ───────────────────────────────────────────────── */
    if (strcmp(cmd, "cat") == 0 || strcmp(cmd, "export") == 0) {
        bool is_export = strcmp(cmd, "export") == 0;
        for (int i = 2; i < argc; i++) {
            if (is_help_flag(argv[i])) {
                if (is_export) print_help_export(); else print_help_cat();
                return 0;
            }
        }

        row_select_opts_t opts;
        opts.offset = 0;
        opts.limit = -1;
        opts.columns = NULL;
        opts.filter = NULL;
        const char* file_path = NULL;
        const char* format = "csv";  /* only used by export */

        for (int i = 2; i < argc; i++) {
            if ((strcmp(argv[i], "-n") == 0 || strcmp(argv[i], "--limit") == 0) && i + 1 < argc) {
                int64_t v;
                if (parse_int64(argv[++i], &v) != 0) {
                    fprintf(stderr, "error: invalid --limit '%s'\n", argv[i]);
                    return 1;
                }
                opts.limit = v;
            } else if ((strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--offset") == 0) && i + 1 < argc) {
                int64_t v;
                if (parse_int64(argv[++i], &v) != 0) {
                    fprintf(stderr, "error: invalid --offset '%s'\n", argv[i]);
                    return 1;
                }
                opts.offset = v;
            } else if ((strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--columns") == 0) && i + 1 < argc) {
                opts.columns = argv[++i];
            } else if ((strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--filter") == 0) && i + 1 < argc) {
                opts.filter = argv[++i];
            } else if (is_export && strcmp(argv[i], "--format") == 0 && i + 1 < argc) {
                format = argv[++i];
            } else if (argv[i][0] != '-') {
                file_path = argv[i];
            } else {
                fprintf(stderr, "error: unknown option '%s' for '%s'\n\n", argv[i], cmd);
                if (is_export) print_help_export(); else print_help_cat();
                return 1;
            }
        }

        if (!file_path) {
            fprintf(stderr, "error: no input file specified\n\n");
            if (is_export) print_help_export(); else print_help_cat();
            return 1;
        }

        if (is_export) {
            export_format_t fmt;
            if (strcmp(format, "csv") == 0) {
                fmt = CLI_EXPORT_CSV;
            } else {
                fprintf(stderr, "error: unsupported --format '%s' (expected: csv)\n", format);
                return 1;
            }
            return cmd_export(file_path, &opts, fmt);
        }
        return cmd_cat(file_path, &opts);
    }

    /* ── codegen ────────────────────────────────────────────────────── */
    if (strcmp(cmd, "codegen") == 0) {
        /* Check for help */
        for (int i = 2; i < argc; i++) {
            if (is_help_flag(argv[i])) {
                print_help_codegen();
                return 0;
            }
        }

        codegen_opts_t opts = {0};
        opts.batch_size = CLI_DEFAULT_BATCH_SIZE;
        opts.mode = 0; /* read */

        for (int i = 2; i < argc; i++) {
            if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--read") == 0) {
                opts.mode = 0;
            } else if (strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "--write") == 0) {
                opts.mode = 1;
            } else if ((strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--file") == 0) && i + 1 < argc) {
                opts.input_path = argv[++i];
            } else if ((strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) && i + 1 < argc) {
                opts.output_path = argv[++i];
            } else if ((strcmp(argv[i], "-b") == 0 || strcmp(argv[i], "--batch-size") == 0) && i + 1 < argc) {
                int64_t val;
                if (parse_int64(argv[++i], &val) != 0 || val <= 0) {
                    fprintf(stderr, "error: invalid batch size '%s'\n", argv[i]);
                    return 1;
                }
                opts.batch_size = (int32_t)val;
            } else if ((strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--columns") == 0) && i + 1 < argc) {
                opts.columns = argv[++i];
            } else if (strcmp(argv[i], "--mmap") == 0) {
                opts.use_mmap = true;
            } else if (strcmp(argv[i], "--skeleton") == 0) {
                opts.skeleton = true;
            } else {
                fprintf(stderr, "error: unknown codegen option '%s'\n\n", argv[i]);
                print_help_codegen();
                return 1;
            }
        }

        return cmd_codegen(&opts);
    }

    /* ── All other commands: parse [-n NUM] [-h] <file> ────────────── */

    /* Dispatch help based on command name */
    typedef void (*help_fn)(void);
    struct { const char* name; help_fn help; } help_table[] = {
        {"schema",   print_help_schema},
        {"info",     print_help_info},
        {"head",     print_help_head},
        {"tail",     print_help_tail},
        {"count",    print_help_count},
        {"columns",  print_help_columns},
        {"stat",     print_help_stat},
        {"validate", print_help_validate},
        {"sample",   print_help_sample},
    };
    int num_cmds = (int)(sizeof(help_table) / sizeof(help_table[0]));

    /* Check for -h in any position */
    for (int i = 2; i < argc; i++) {
        if (is_help_flag(argv[i])) {
            for (int j = 0; j < num_cmds; j++) {
                if (strcmp(cmd, help_table[j].name) == 0) {
                    help_table[j].help();
                    return 0;
                }
            }
            /* Unknown command with -h */
            fprintf(stderr, "error: unknown command '%s'\n\n", cmd);
            print_usage();
            return 1;
        }
    }

    int64_t num_rows = CLI_DEFAULT_NUM_ROWS;
    const char* file_path = NULL;
    const char* filter = NULL;

    for (int i = 2; i < argc; i++) {
        if ((strcmp(argv[i], "-n") == 0) && i + 1 < argc) {
            if (parse_int64(argv[++i], &num_rows) != 0) {
                fprintf(stderr, "error: invalid number '%s'\n", argv[i]);
                return 1;
            }
        } else if ((strcmp(argv[i], "-p") == 0 ||
                    strcmp(argv[i], "--filter") == 0) && i + 1 < argc) {
            filter = argv[++i];
        } else if (argv[i][0] != '-') {
            file_path = argv[i];
        } else {
            fprintf(stderr, "error: unknown option '%s' for '%s'\n\n", argv[i], cmd);
            /* Try to show command-specific help */
            for (int j = 0; j < num_cmds; j++) {
                if (strcmp(cmd, help_table[j].name) == 0) {
                    help_table[j].help();
                    return 1;
                }
            }
            print_usage();
            return 1;
        }
    }

    if (!file_path) {
        fprintf(stderr, "error: no input file specified\n\n");
        /* Show command-specific help if valid command */
        for (int j = 0; j < num_cmds; j++) {
            if (strcmp(cmd, help_table[j].name) == 0) {
                help_table[j].help();
                return 1;
            }
        }
        print_usage();
        return 1;
    }

    if (strcmp(cmd, "schema") == 0)        return cmd_schema(file_path);
    if (strcmp(cmd, "info") == 0)          return cmd_info(file_path);
    if (strcmp(cmd, "head") == 0)          return cmd_head(file_path, num_rows, filter);
    if (strcmp(cmd, "tail") == 0)          return cmd_tail(file_path, num_rows, filter);
    if (strcmp(cmd, "count") == 0)         return cmd_count(file_path, filter);
    if (strcmp(cmd, "columns") == 0)       return cmd_columns(file_path);
    if (strcmp(cmd, "stat") == 0)          return cmd_stat(file_path);
    if (strcmp(cmd, "validate") == 0)      return cmd_validate(file_path);
    if (strcmp(cmd, "sample") == 0)        return cmd_sample(file_path, num_rows, filter);

    fprintf(stderr, "error: unknown command '%s'\n\n", cmd);
    print_usage();
    return 1;
}
