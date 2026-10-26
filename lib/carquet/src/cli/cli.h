/**
 * @file cli.h
 * @brief Shared declarations for carquet CLI commands
 */

#ifndef CARQUET_CLI_H
#define CARQUET_CLI_H

#include <carquet/carquet.h>
#include <stdio.h>
#include <stdint.h>

/* Maximum columns we support displaying in head/tail */
#define CLI_MAX_DISPLAY_COLS 256

/* Default number of rows for head/tail */
#define CLI_DEFAULT_NUM_ROWS 10

/* Default batch size for codegen */
#define CLI_DEFAULT_BATCH_SIZE 1024

/* ── Command handlers ─────────────────────────────────────────────────── */

int cmd_schema(const char* path);
int cmd_info(const char* path);
int cmd_head(const char* path, int64_t n, const char* filter);
int cmd_tail(const char* path, int64_t n, const char* filter);
int cmd_count(const char* path, const char* filter);
int cmd_columns(const char* path);
int cmd_stat(const char* path);
int cmd_validate(const char* path);
int cmd_sample(const char* path, int64_t n, const char* filter);

/* Options for `cat` and `export` (subset selection + slicing). */
typedef struct row_select_opts {
    int64_t     offset;    /* rows to skip from the start */
    int64_t     limit;     /* -1 = all remaining */
    const char* columns;   /* comma-separated names; NULL = all columns */
    const char* filter;    /* page-filter expression, NULL = no filter */
} row_select_opts_t;

typedef enum export_format {
    CLI_EXPORT_CSV = 0,
} export_format_t;

int cmd_cat(const char* path, const row_select_opts_t* opts);
int cmd_export(const char* path, const row_select_opts_t* opts, export_format_t fmt);

/* ── Codegen ──────────────────────────────────────────────────────────── */

typedef struct codegen_opts {
    const char* input_path;     /* NULL = no file (generate placeholder) */
    const char* output_path;    /* NULL = stdout */
    int32_t     batch_size;
    const char* columns;        /* comma-separated column filter, NULL = all */
    int         mode;           /* 0 = read, 1 = write */
    bool        use_mmap;       /* generate mmap-based reader */
    bool        skeleton;       /* empty process_batch body */
} codegen_opts_t;

/** Hints returned by codegen for user messages */
typedef struct codegen_hints {
    char build_line[2048];
    int  default_file_line;     /* line number of DEFAULT_FILE, 0 = not emitted */
    int  process_batch_line;    /* line number of process_batch body, 0 = N/A */
} codegen_hints_t;

int cmd_codegen(const codegen_opts_t* opts);

/* Internal: called by cmd_codegen */
int cmd_codegen_read(FILE* out, carquet_reader_t* reader,
                     const codegen_opts_t* opts,
                     codegen_hints_t* hints);
int cmd_codegen_write(const codegen_opts_t* opts);

/* ── Helpers ──────────────────────────────────────────────────────────── */

/** Format a physical+logical type into a human-readable string */
void cli_format_type(carquet_physical_type_t phys,
                     const carquet_logical_type_t* logical,
                     char* buf, size_t buf_size);

/** Format a byte count as human-readable (e.g. "1.2 MB") */
void cli_format_bytes(int64_t bytes, char* buf, size_t buf_size);

/** Format a value to string based on physical type. Returns buf. */
const char* cli_format_value(carquet_physical_type_t type,
                             const void* value, int32_t type_len,
                             const carquet_logical_type_t* logical,
                             char* buf, size_t buf_size);

/** Repetition name */
const char* cli_repetition_name(carquet_field_repetition_t rep);

#endif /* CARQUET_CLI_H */
