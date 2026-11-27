/**
 * @file commands.c
 * @brief Implementation of carquet CLI commands
 */

#include "cli.h"
#include "core/compat.h"
#include "core/float16.h"
#include "reader/reader_internal.h"
#include "thrift/parquet_types.h"
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>

/* Shared layout constants (defined early so commands above the table
 * helpers below can reference them). */
#define MAX_COL_WIDTH 40
#define MAX_VALUE_BUF 256

/* Forward declaration: implementation lives further down with the other
 * tabular-output helpers. */
static void print_dyn_table(const char* const* headers, int32_t num_cols,
                             const char* const* cells, int64_t num_rows);

/* Types and forward declarations for the filter / batch-reader path
 * used by cmd_count, cmd_head, cmd_cat, and cmd_export. The
 * implementations live further down. */

typedef struct str_matrix {
    char** cells;     /* [num_rows * num_cols], heap-strdup'd; may be NULL */
    int64_t num_rows;
    int32_t num_cols;
} str_matrix_t;

typedef struct cli_filter_storage {
    carquet_filter_clause_t* clauses;
    int32_t count;
    int32_t capacity;
    /* Backing storage: one heap buffer per clause's value (or NULL for
     * IS_NULL / IS_NOT_NULL). The clause's value pointer aliases into
     * this buffer; freeing storage invalidates every clause. */
    uint8_t** blobs;
    int32_t num_blobs;
} cli_filter_storage_t;

static int cli_parse_filter(const char* expr,
                            const carquet_schema_t* schema,
                            int32_t num_cols,
                            cli_filter_storage_t* out,
                            char* err, size_t errsz);
static void cli_filter_free(cli_filter_storage_t* s);
static int read_rows_filtered(carquet_reader_t* reader,
                              const carquet_schema_t* schema,
                              const int32_t* col_indices, int32_t num_sel_cols,
                              int64_t offset, int64_t limit,
                              const cli_filter_storage_t* filter,
                              const char* null_text, str_matrix_t* out);
static int read_rows(carquet_reader_t* reader,
                     const carquet_schema_t* schema,
                     const int32_t* col_indices, int32_t num_sel_cols,
                     int64_t offset, int64_t limit, const int64_t* picks,
                     const char* null_text, str_matrix_t* out);
static int64_t count_rows_filtered(carquet_reader_t* reader,
                                   const cli_filter_storage_t* filter);
static void matrix_free(str_matrix_t* m);

/* ══════════════════════════════════════════════════════════════════════════
 * Helpers
 * ══════════════════════════════════════════════════════════════════════════ */

const char* cli_repetition_name(carquet_field_repetition_t rep) {
    switch (rep) {
        case CARQUET_REPETITION_REQUIRED: return "REQUIRED";
        case CARQUET_REPETITION_OPTIONAL: return "OPTIONAL";
        case CARQUET_REPETITION_REPEATED: return "REPEATED";
        default: return "?";
    }
}

void cli_format_type(carquet_physical_type_t phys,
                     const carquet_logical_type_t* logical,
                     char* buf, size_t buf_size)
{
    const char* base = carquet_physical_type_name(phys);
    if (!logical || logical->id == CARQUET_LOGICAL_UNKNOWN) {
        snprintf(buf, buf_size, "%s", base);
        return;
    }
    switch (logical->id) {
        case CARQUET_LOGICAL_STRING:    snprintf(buf, buf_size, "STRING"); break;
        case CARQUET_LOGICAL_DATE:      snprintf(buf, buf_size, "DATE"); break;
        case CARQUET_LOGICAL_UUID:      snprintf(buf, buf_size, "UUID"); break;
        case CARQUET_LOGICAL_JSON:      snprintf(buf, buf_size, "JSON"); break;
        case CARQUET_LOGICAL_ENUM:      snprintf(buf, buf_size, "ENUM"); break;
        case CARQUET_LOGICAL_LIST:      snprintf(buf, buf_size, "LIST"); break;
        case CARQUET_LOGICAL_MAP:       snprintf(buf, buf_size, "MAP"); break;
        case CARQUET_LOGICAL_FLOAT16:   snprintf(buf, buf_size, "FLOAT16"); break;
        case CARQUET_LOGICAL_VARIANT:   snprintf(buf, buf_size, "VARIANT"); break;
        case CARQUET_LOGICAL_GEOMETRY:  snprintf(buf, buf_size, "GEOMETRY"); break;
        case CARQUET_LOGICAL_GEOGRAPHY: snprintf(buf, buf_size, "GEOGRAPHY"); break;
        case CARQUET_LOGICAL_NULL:      snprintf(buf, buf_size, "NULL"); break;
        case CARQUET_LOGICAL_BSON:      snprintf(buf, buf_size, "BSON"); break;
        case CARQUET_LOGICAL_INTERVAL:  snprintf(buf, buf_size, "INTERVAL"); break;
        case CARQUET_LOGICAL_DECIMAL:
            snprintf(buf, buf_size, "DECIMAL(%d,%d)",
                     logical->params.decimal.precision,
                     logical->params.decimal.scale);
            break;
        case CARQUET_LOGICAL_INTEGER:
            snprintf(buf, buf_size, "%sINT%d",
                     logical->params.integer.is_signed ? "" : "U",
                     logical->params.integer.bit_width);
            break;
        case CARQUET_LOGICAL_TIME: {
            const char* unit = "?";
            switch (logical->params.time.unit) {
                case CARQUET_TIME_UNIT_MILLIS: unit = "ms"; break;
                case CARQUET_TIME_UNIT_MICROS: unit = "us"; break;
                case CARQUET_TIME_UNIT_NANOS:  unit = "ns"; break;
            }
            snprintf(buf, buf_size, "TIME(%s%s)", unit,
                     logical->params.time.is_adjusted_to_utc ? ",UTC" : "");
            break;
        }
        case CARQUET_LOGICAL_TIMESTAMP: {
            const char* unit = "?";
            switch (logical->params.timestamp.unit) {
                case CARQUET_TIME_UNIT_MILLIS: unit = "ms"; break;
                case CARQUET_TIME_UNIT_MICROS: unit = "us"; break;
                case CARQUET_TIME_UNIT_NANOS:  unit = "ns"; break;
            }
            snprintf(buf, buf_size, "TIMESTAMP(%s%s)", unit,
                     logical->params.timestamp.is_adjusted_to_utc ? ",UTC" : "");
            break;
        }
        default:
            snprintf(buf, buf_size, "%s", base);
            break;
    }
}

void cli_format_bytes(int64_t bytes, char* buf, size_t buf_size) {
    if (bytes < 1024)
        snprintf(buf, buf_size, "%" PRId64 " B", bytes);
    else if (bytes < 1024 * 1024)
        snprintf(buf, buf_size, "%.1f KB", bytes / 1024.0);
    else if (bytes < 1024LL * 1024 * 1024)
        snprintf(buf, buf_size, "%.1f MB", bytes / (1024.0 * 1024));
    else
        snprintf(buf, buf_size, "%.2f GB", bytes / (1024.0 * 1024 * 1024));
}

/* Days since 1970-01-01 -> proleptic Gregorian y/m/d (Howard Hinnant's
 * civil_from_days). Deliberately not gmtime_r/gmtime_s: the MSVC CRT rejects
 * any time outside [1970-01-01, 3000-12-31] and leaves `struct tm` untouched,
 * so Parquet's full DATE/TIMESTAMP range (year 1 through year 9999, and
 * pre-epoch values everywhere) printed uninitialized garbage on Windows. This
 * is total over the whole int64 domain and identical on every platform. */
static void civil_from_days(int64_t z, int64_t* year, int* month, int* day) {
    z += 719468; /* shift the epoch to 0000-03-01 so leap day lands last */
    int64_t era = (z >= 0 ? z : z - 146096) / 146097;
    int64_t doe = z - era * 146097;                                   /* [0, 146096] */
    int64_t yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365; /* [0, 399] */
    int64_t doy = doe - (365 * yoe + yoe / 4 - yoe / 100);            /* [0, 365] */
    int64_t mp  = (5 * doy + 2) / 153;                                /* [0, 11], March-based */
    *day   = (int)(doy - (153 * mp + 2) / 5 + 1);                     /* [1, 31] */
    *month = (int)(mp < 10 ? mp + 3 : mp - 9);                        /* [1, 12] */
    *year  = yoe + era * 400 + (*month <= 2);
}

/* Split a signed count into a floored quotient and a non-negative remainder,
 * so pre-epoch values decompose correctly (truncating division would round
 * toward zero and yield a negative remainder). */
static int64_t floor_div_mod(int64_t val, int64_t divisor, int64_t* rem) {
    int64_t q = val / divisor;
    int64_t r = val % divisor;
    if (r < 0) { r += divisor; q -= 1; }
    *rem = r;
    return q;
}

const char* cli_format_value(carquet_physical_type_t type,
                             const void* value, int32_t type_len,
                             const carquet_logical_type_t* logical,
                             char* buf, size_t buf_size)
{
    if (!value) { snprintf(buf, buf_size, "null"); return buf; }

    /* Handle logical type formatting */
    if (logical && logical->id == CARQUET_LOGICAL_DATE && type == CARQUET_PHYSICAL_INT32) {
        int64_t y; int m, d;
        civil_from_days(*(const int32_t*)value, &y, &m, &d);
        snprintf(buf, buf_size, "%04" PRId64 "-%02d-%02d", y, m, d);
        return buf;
    }

    /* TIME is a time-of-day offset from midnight: MILLIS is INT32-backed,
     * MICROS/NANOS are INT64-backed. Rendered as HH:MM:SS.frac with the
     * fractional part always at the unit's full precision. */
    if (logical && logical->id == CARQUET_LOGICAL_TIME &&
        (type == CARQUET_PHYSICAL_INT32 || type == CARQUET_PHYSICAL_INT64)) {
        int64_t val = (type == CARQUET_PHYSICAL_INT32)
                          ? (int64_t)*(const int32_t*)value
                          : *(const int64_t*)value;
        uint64_t divisor = 1000;
        int frac_digits = 3;
        switch (logical->params.time.unit) {
            case CARQUET_TIME_UNIT_MILLIS: divisor = 1000ULL;       frac_digits = 3; break;
            case CARQUET_TIME_UNIT_MICROS: divisor = 1000000ULL;    frac_digits = 6; break;
            case CARQUET_TIME_UNIT_NANOS:  divisor = 1000000000ULL; frac_digits = 9; break;
        }
        /* Spec-valid values are in [0, 24h); out-of-range ones still print
         * readably (a sign and/or an hour field past 23) instead of wrapping.
         * Negation goes through the magnitude to stay INT64_MIN-safe. */
        bool neg = val < 0;
        uint64_t mag = neg ? (uint64_t)(-(val + 1)) + 1ULL : (uint64_t)val;
        uint64_t secs = mag / divisor;
        uint32_t frac = (uint32_t)(mag % divisor);
        snprintf(buf, buf_size, "%s%02" PRIu64 ":%02" PRIu64 ":%02" PRIu64 ".%0*" PRIu32,
                 neg ? "-" : "", secs / 3600, (secs / 60) % 60, secs % 60,
                 frac_digits, frac);
        return buf;
    }

    if (logical && logical->id == CARQUET_LOGICAL_TIMESTAMP &&
        type == CARQUET_PHYSICAL_INT64) {
        int64_t val = *(const int64_t*)value;
        const char* frac_fmt = "";
        int64_t divisor = 1;
        switch (logical->params.timestamp.unit) {
            case CARQUET_TIME_UNIT_MILLIS:
                divisor = 1000;
                frac_fmt = ".%03d";
                break;
            case CARQUET_TIME_UNIT_MICROS:
                divisor = 1000000;
                frac_fmt = ".%06d";
                break;
            case CARQUET_TIME_UNIT_NANOS:
                divisor = 1000000000LL;
                frac_fmt = ".%09d";
                break;
        }
        /* Floor division at each step so pre-epoch (negative) values split
         * correctly: truncating division would push the quotient up by one and
         * make the remainder negative (e.g. -999 ms -> 0 s and -999 ms instead
         * of -1 s and 1 ms). */
        int64_t frac_val, secs_of_day;
        int64_t secs = floor_div_mod(val, divisor, &frac_val);
        int64_t days = floor_div_mod(secs, 86400, &secs_of_day);
        int64_t y; int m, d;
        civil_from_days(days, &y, &m, &d);
        int n = snprintf(buf, buf_size,
                         "%04" PRId64 "-%02d-%02dT%02" PRId64 ":%02" PRId64 ":%02" PRId64,
                         y, m, d, secs_of_day / 3600, (secs_of_day / 60) % 60,
                         secs_of_day % 60);
        if (frac_val != 0 && n > 0 && (size_t)n < buf_size)
            snprintf(buf + n, buf_size - (size_t)n, frac_fmt, (int)frac_val);
        return buf;
    }

    switch (type) {
        case CARQUET_PHYSICAL_BOOLEAN:
            snprintf(buf, buf_size, "%s", *(const uint8_t*)value ? "true" : "false");
            break;
        case CARQUET_PHYSICAL_INT32:
            snprintf(buf, buf_size, "%" PRId32, *(const int32_t*)value);
            break;
        case CARQUET_PHYSICAL_INT64:
            snprintf(buf, buf_size, "%" PRId64, *(const int64_t*)value);
            break;
        case CARQUET_PHYSICAL_FLOAT:
            snprintf(buf, buf_size, "%g", (double)*(const float*)value);
            break;
        case CARQUET_PHYSICAL_DOUBLE:
            snprintf(buf, buf_size, "%g", *(const double*)value);
            break;
        case CARQUET_PHYSICAL_BYTE_ARRAY: {
            const carquet_byte_array_t* ba = (const carquet_byte_array_t*)value;
            /* Check if it looks like a string (logical STRING or UTF8) */
            bool is_string = logical && (logical->id == CARQUET_LOGICAL_STRING ||
                                          logical->id == CARQUET_LOGICAL_JSON ||
                                          logical->id == CARQUET_LOGICAL_ENUM);
            if (is_string || 1) {
                /* Try to print as string, truncate if long */
                int32_t len = ba->length;
                int32_t max_len = (int32_t)(buf_size - 1);
                if (len > max_len) len = max_len;
                memcpy(buf, ba->data, (size_t)len);
                buf[len] = '\0';
            }
            break;
        }
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY: {
            /* FLOAT16 (FLBA length 2) prints as its float value, not hex. */
            if (logical && logical->id == CARQUET_LOGICAL_FLOAT16 &&
                type_len == 2) {
                const uint8_t* b = (const uint8_t*)value;
                snprintf(buf, buf_size, "%g",
                         (double)carquet_half_to_float(
                             (uint16_t)(b[0] | (b[1] << 8))));
                break;
            }
            /* Print as hex */
            const uint8_t* bytes = (const uint8_t*)value;
            int32_t len = type_len;
            if (len > (int32_t)(buf_size / 2 - 1)) len = (int32_t)(buf_size / 2 - 1);
            for (int32_t i = 0; i < len; i++)
                snprintf(buf + i * 2, buf_size - (size_t)(i * 2), "%02x", bytes[i]);
            break;
        }
        case CARQUET_PHYSICAL_INT96: {
            const uint32_t* v96 = (const uint32_t*)value;
            snprintf(buf, buf_size, "0x%08x%08x%08x", v96[2], v96[1], v96[0]);
            break;
        }
        default:
            snprintf(buf, buf_size, "?");
            break;
    }
    return buf;
}

const char* cli_format_stat_value(carquet_physical_type_t type,
                                  const void* raw, int32_t raw_size,
                                  int32_t type_len,
                                  const carquet_logical_type_t* logical,
                                  char* buf, size_t buf_size)
{
    /* Raw stat bytes come straight from the file: check the width before
     * cli_format_value() dereferences a fixed-size value, and copy numeric
     * values into an aligned buffer (arena bytes carry no alignment). */
    size_t need;
    switch (type) {
        case CARQUET_PHYSICAL_BOOLEAN: need = 1; break;
        case CARQUET_PHYSICAL_INT32:
        case CARQUET_PHYSICAL_FLOAT:   need = 4; break;
        case CARQUET_PHYSICAL_INT64:
        case CARQUET_PHYSICAL_DOUBLE:  need = 8; break;
        case CARQUET_PHYSICAL_INT96:   need = 12; break;
        case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
            need = type_len > 0 ? (size_t)type_len : SIZE_MAX; break;
        default: need = SIZE_MAX; break;  /* BYTE_ARRAY: caller wraps it */
    }
    if (!raw || raw_size < 0 || (size_t)raw_size < need) {
        snprintf(buf, buf_size, "(invalid)");
        return buf;
    }
    if (type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY) {
        return cli_format_value(type, raw, type_len, logical, buf, buf_size);
    }
    uint64_t aligned[2] = {0, 0};
    memcpy(aligned, raw, need);
    return cli_format_value(type, aligned, type_len, logical, buf, buf_size);
}

static carquet_reader_t* open_or_die(const char* path, carquet_error_t* err) {
    carquet_reader_t* reader = carquet_reader_open(path, NULL, err);
    if (!reader) {
        fprintf(stderr, "error: %s\n", err->message);
    }
    return reader;
}

/* ══════════════════════════════════════════════════════════════════════════
 * cmd_schema
 * ══════════════════════════════════════════════════════════════════════════ */

int cmd_schema(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = open_or_die(path, &err);
    if (!reader) return 1;

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t n = carquet_schema_num_elements(schema);

    printf("message schema {\n");

    /*
     * Parquet schema is stored as a flat list with num_children to define tree
     * structure (Thrift-style DFS pre-order). We track depth with a stack of
     * remaining children counts.
     */
    int depth_stack[64] = {0};
    int depth = 0;
    /* Element 0 is root "schema", its children count tells us how many
     * top-level elements follow */
    const carquet_schema_node_t* root = carquet_schema_get_element(schema, 0);
    /* Get num_children from internal struct */
    const parquet_schema_element_t* root_elem = (const parquet_schema_element_t*)root;
    depth_stack[0] = root_elem->num_children;
    depth = 1;

    for (int32_t i = 1; i < n; i++) {
        const carquet_schema_node_t* node = carquet_schema_get_element(schema, i);
        const parquet_schema_element_t* elem = (const parquet_schema_element_t*)node;

        /* Indent */
        for (int d = 0; d < depth; d++) printf("  ");

        if (carquet_schema_node_is_leaf(node)) {
            char type_buf[64];
            cli_format_type(carquet_schema_node_physical_type(node),
                           carquet_schema_node_logical_type(node),
                           type_buf, sizeof(type_buf));
            printf("%s %s %s",
                   cli_repetition_name(carquet_schema_node_repetition(node)),
                   type_buf,
                   carquet_schema_node_name(node));

            int32_t tl = carquet_schema_node_type_length(node);
            if (tl > 0) printf(" (length=%d)", tl);
            printf(";\n");
        } else {
            /* Group node */
            const carquet_logical_type_t* lt = carquet_schema_node_logical_type(node);
            const char* annotation = "";
            if (lt) {
                switch (lt->id) {
                    case CARQUET_LOGICAL_LIST: annotation = " (LIST)"; break;
                    case CARQUET_LOGICAL_MAP:  annotation = " (MAP)"; break;
                    case CARQUET_LOGICAL_VARIANT: annotation = " (VARIANT)"; break;
                    case CARQUET_LOGICAL_GEOMETRY: annotation = " (GEOMETRY)"; break;
                    case CARQUET_LOGICAL_GEOGRAPHY: annotation = " (GEOGRAPHY)"; break;
                    default: break;
                }
            }
            printf("%s group %s%s {\n",
                   cli_repetition_name(carquet_schema_node_repetition(node)),
                   carquet_schema_node_name(node),
                   annotation);

            /* Push children count */
            if (depth < 63) {
                depth++;
                depth_stack[depth - 1] = elem->num_children;
            }
            continue; /* Don't decrement children count yet */
        }

        /* Decrement parent's children count and close groups */
        depth_stack[depth - 1]--;
        while (depth > 1 && depth_stack[depth - 1] == 0) {
            depth--;
            for (int d = 0; d < depth; d++) printf("  ");
            printf("}\n");
            if (depth > 0) depth_stack[depth - 1]--;
        }
    }

    printf("}\n");
    carquet_reader_close(reader);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * cmd_info
 * ══════════════════════════════════════════════════════════════════════════ */

int cmd_info(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = open_or_die(path, &err);
    if (!reader) return 1;

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int64_t total_rows = carquet_reader_num_rows(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    int32_t num_rgs = carquet_reader_num_row_groups(reader);

    /* Access internal metadata for created_by and key-value metadata */
    const parquet_file_metadata_t* meta = &reader->metadata;

    printf("File:         %s\n", path);
    if (meta && meta->created_by)
        printf("Created by:   %s\n", meta->created_by);
    printf("Rows:         %" PRId64 "\n", total_rows);
    printf("Columns:      %d\n", num_cols);
    printf("Row groups:   %d\n", num_rgs);

    /* Key-value metadata */
    if (meta && meta->num_key_value > 0) {
        printf("\nKey-value metadata:\n");
        for (int32_t i = 0; i < meta->num_key_value; i++) {
            const char* val = meta->key_value_metadata[i].value;
            if (val && strlen(val) > 60) {
                printf("  %-20s %.57s...\n", meta->key_value_metadata[i].key, val);
            } else {
                printf("  %-20s %s\n", meta->key_value_metadata[i].key,
                       val ? val : "(null)");
            }
        }
    }

    /* Column details */
    printf("\nColumns:\n");
    printf("  %-4s %-30s %-20s %-10s\n", "#", "Name", "Type", "Nullable");
    printf("  %-4s %-30s %-20s %-10s\n", "---", "---", "---", "---");
    for (int32_t c = 0; c < num_cols; c++) {
        char type_buf[64];
        cli_format_type(carquet_schema_column_type(schema, c),
                       carquet_schema_node_logical_type(
                           carquet_schema_get_element(schema,
                               schema->leaf_indices[c])),
                       type_buf, sizeof(type_buf));

        const carquet_schema_node_t* node = carquet_schema_get_element(schema,
            schema->leaf_indices[c]);
        bool nullable = carquet_schema_node_repetition(node) != CARQUET_REPETITION_REQUIRED;

        char idx[8];
        snprintf(idx, sizeof(idx), "%d", c);
        printf("  %-4s %-30s %-20s %-10s\n", idx,
               carquet_schema_column_name(schema, c),
               type_buf, nullable ? "yes" : "no");
    }

    /* Row group details */
    printf("\nRow groups:\n");
    printf("  %-4s %-15s %-15s %-15s %-10s\n",
           "#", "Rows", "Uncompressed", "Compressed", "Ratio");
    printf("  %-4s %-15s %-15s %-15s %-10s\n",
           "---", "---", "---", "---", "---");
    for (int32_t rg = 0; rg < num_rgs; rg++) {
        carquet_row_group_metadata_t rgm;
        if (carquet_reader_row_group_metadata(reader, rg, &rgm) != CARQUET_OK)
            continue;
        char uncomp[32], comp[32], ratio[16], idx[8];
        cli_format_bytes(rgm.total_byte_size, uncomp, sizeof(uncomp));
        cli_format_bytes(rgm.total_compressed_size, comp, sizeof(comp));
        if (rgm.total_byte_size > 0)
            snprintf(ratio, sizeof(ratio), "%.1fx",
                     (double)rgm.total_byte_size / (double)rgm.total_compressed_size);
        else
            snprintf(ratio, sizeof(ratio), "-");
        snprintf(idx, sizeof(idx), "%d", rg);
        printf("  %-4s %-15" PRId64 " %-15s %-15s %-10s\n",
               idx, rgm.num_rows, uncomp, comp, ratio);
    }

    /* Sort order: dump the first row group's sorting_columns if present.
     * Carquet writers record the same list on every row group, so the
     * first is representative. */
    if (num_rgs > 0 && meta && meta->row_groups[0].num_sorting_columns > 0) {
        const parquet_row_group_t* rg0 = &meta->row_groups[0];
        printf("\nSort order:\n");
        for (int32_t i = 0; i < rg0->num_sorting_columns; i++) {
            const parquet_sorting_column_t* sc = &rg0->sorting_columns[i];
            const char* nm = (sc->column_idx >= 0 && sc->column_idx < num_cols)
                ? carquet_schema_column_name(schema, sc->column_idx) : "?";
            printf("  %s %s NULLS %s\n", nm,
                   sc->descending ? "DESC" : "ASC",
                   sc->nulls_first ? "FIRST" : "LAST");
        }
    }

    /* Page index summary: per-column, sampled from the first row group
     * (all row groups produced by carquet have the same coverage). Shown
     * only when at least one column has a column index, so files that
     * were written without page index don't add a section. */
    if (num_rgs > 0) {
        bool any = false;
        for (int32_t c = 0; c < num_cols; c++) {
            carquet_error_t ie = CARQUET_ERROR_INIT;
            carquet_column_index_t* ci =
                carquet_reader_get_column_index(reader, 0, c, &ie);
            if (ci) {
                any = true;
                carquet_column_index_free(ci);
                break;
            }
        }
        if (any) {
            printf("\nPage index:\n");
            printf("  %-4s %-30s %-8s %-12s\n",
                   "#", "Name", "Pages", "Boundary");
            printf("  %-4s %-30s %-8s %-12s\n",
                   "---", "---", "---", "---");
            for (int32_t c = 0; c < num_cols; c++) {
                carquet_error_t ie = CARQUET_ERROR_INIT;
                carquet_column_index_t* ci =
                    carquet_reader_get_column_index(reader, 0, c, &ie);
                const char* nm = carquet_schema_column_name(schema, c);
                if (!ci) {
                    char idx[8];
                    snprintf(idx, sizeof(idx), "%d", c);
                    printf("  %-4s %-30s %-8s %-12s\n",
                           idx, nm ? nm : "?", "-", "-");
                    continue;
                }
                int32_t np = carquet_column_index_num_pages(ci);
                int32_t bo = carquet_column_index_boundary_order(ci);
                const char* bo_name = "UNORDERED";
                if (bo == 1) bo_name = "ASCENDING";
                else if (bo == 2) bo_name = "DESCENDING";
                char idx[8], pages[16];
                snprintf(idx, sizeof(idx), "%d", c);
                snprintf(pages, sizeof(pages), "%d", np);
                printf("  %-4s %-30s %-8s %-12s\n",
                       idx, nm ? nm : "?", pages, bo_name);
                carquet_column_index_free(ci);
            }
        }
    }

    carquet_reader_close(reader);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * cmd_count
 * ══════════════════════════════════════════════════════════════════════════ */

int cmd_count(const char* path, const char* filter) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = open_or_die(path, &err);
    if (!reader) return 1;

    if (!filter) {
        printf("%" PRId64 "\n", carquet_reader_num_rows(reader));
        carquet_reader_close(reader);
        return 0;
    }

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    cli_filter_storage_t fs;
    char ferr[256];
    if (cli_parse_filter(filter, schema, num_cols, &fs, ferr, sizeof(ferr)) != 0) {
        fprintf(stderr, "error: %s\n", ferr);
        carquet_reader_close(reader);
        return 1;
    }
    int64_t total = count_rows_filtered(reader, &fs);
    cli_filter_free(&fs);
    carquet_reader_close(reader);
    if (total < 0) return 1;
    printf("%" PRId64 "\n", total);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * cmd_columns
 * ══════════════════════════════════════════════════════════════════════════ */

int cmd_columns(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = open_or_die(path, &err);
    if (!reader) return 1;

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    for (int32_t c = 0; c < num_cols; c++) {
        printf("%s\n", carquet_schema_column_name(schema, c));
    }
    carquet_reader_close(reader);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * cmd_stat
 * ══════════════════════════════════════════════════════════════════════════ */

int cmd_stat(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = open_or_die(path, &err);
    if (!reader) return 1;

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    int32_t num_rgs = carquet_reader_num_row_groups(reader);

    static const char* const HEADERS[] = {"Column", "Type", "Nulls", "Min", "Max"};
    const int32_t NCOLS = 5;

    char** cells = calloc((size_t)num_cols * NCOLS, sizeof(char*));
    if (!cells) {
        carquet_reader_close(reader);
        return 1;
    }

    for (int32_t rg = 0; rg < num_rgs; rg++) {
        if (num_rgs > 1)
            printf("Row group %d:\n", rg);

        for (int32_t c = 0; c < num_cols; c++) {
            carquet_column_statistics_t stats;
            carquet_physical_type_t phys = carquet_schema_column_type(schema, c);
            const carquet_schema_node_t* node = carquet_schema_get_element(schema,
                schema->leaf_indices[c]);
            const carquet_logical_type_t* lt = carquet_schema_node_logical_type(node);
            int32_t tl = carquet_schema_node_type_length(node);

            char type_buf[64];
            cli_format_type(phys, lt, type_buf, sizeof(type_buf));

            char nulls[32] = "-";
            char min_buf[MAX_VALUE_BUF] = "-";
            char max_buf[MAX_VALUE_BUF] = "-";

            if (carquet_reader_column_statistics(reader, rg, c, &stats) == CARQUET_OK) {
                if (stats.has_null_count)
                    snprintf(nulls, sizeof(nulls), "%" PRId64, stats.null_count);
                if (stats.has_min_max) {
                    /* stats.min_value / max_value are raw bytes for BYTE_ARRAY.
                     * cli_format_value expects a carquet_byte_array_t* for that
                     * physical type, so wrap the raw bytes here. */
                    if (phys == CARQUET_PHYSICAL_BYTE_ARRAY) {
                        carquet_byte_array_t min_ba = {
                            .data = (uint8_t*)(uintptr_t)stats.min_value,
                            .length = stats.min_value_size
                        };
                        carquet_byte_array_t max_ba = {
                            .data = (uint8_t*)(uintptr_t)stats.max_value,
                            .length = stats.max_value_size
                        };
                        cli_format_value(phys, &min_ba, tl, lt,
                                         min_buf, sizeof(min_buf));
                        cli_format_value(phys, &max_ba, tl, lt,
                                         max_buf, sizeof(max_buf));
                    } else {
                        cli_format_stat_value(phys, stats.min_value,
                                              stats.min_value_size, tl, lt,
                                              min_buf, sizeof(min_buf));
                        cli_format_stat_value(phys, stats.max_value,
                                              stats.max_value_size, tl, lt,
                                              max_buf, sizeof(max_buf));
                    }
                }
            }

            /* GEOMETRY/GEOGRAPHY have no min/max; surface the bounding box
             * and ISO-WKB type codes from GeospatialStatistics instead. */
            if (lt && (lt->id == CARQUET_LOGICAL_GEOMETRY ||
                       lt->id == CARQUET_LOGICAL_GEOGRAPHY)) {
                carquet_geospatial_statistics_t gs;
                if (carquet_reader_geospatial_statistics(reader, rg, c, &gs)
                        == CARQUET_OK) {
                    if (gs.has_bbox) {
                        char zb[48] = "";
                        if (gs.has_z)
                            snprintf(zb, sizeof(zb), " z[%g,%g]",
                                     gs.zmin, gs.zmax);
                        snprintf(min_buf, sizeof(min_buf),
                                 "bbox x[%g,%g] y[%g,%g]%s",
                                 gs.xmin, gs.xmax, gs.ymin, gs.ymax, zb);
                    }
                    int off = snprintf(max_buf, sizeof(max_buf), "types[");
                    for (int32_t t = 0; t < gs.num_geometry_types &&
                         off < (int)sizeof(max_buf) - 8; t++) {
                        off += snprintf(max_buf + off, sizeof(max_buf) - off,
                                        "%s%d", t ? "," : "",
                                        gs.geometry_types[t]);
                    }
                    snprintf(max_buf + off, sizeof(max_buf) - off, "]");
                }
            }

            cells[c * NCOLS + 0] = carquet_heap_strdup(carquet_schema_column_name(schema, c));
            cells[c * NCOLS + 1] = carquet_heap_strdup(type_buf);
            cells[c * NCOLS + 2] = carquet_heap_strdup(nulls);
            cells[c * NCOLS + 3] = carquet_heap_strdup(min_buf);
            cells[c * NCOLS + 4] = carquet_heap_strdup(max_buf);
        }

        print_dyn_table(HEADERS, NCOLS, (const char* const*)cells, num_cols);

        /* Free this row group's cells before reusing the buffer. */
        for (int32_t i = 0; i < num_cols * NCOLS; i++) {
            free(cells[i]);
            cells[i] = NULL;
        }
        if (rg < num_rgs - 1) printf("\n");
    }

    free(cells);
    carquet_reader_close(reader);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * cmd_validate
 * ══════════════════════════════════════════════════════════════════════════ */

int cmd_validate(const char* path) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* Open with checksum verification enabled */
    carquet_reader_options_t opts;
    carquet_reader_options_init(&opts);
    opts.verify_checksums = true;

    carquet_reader_t* reader = carquet_reader_open(path, &opts, &err);
    if (!reader) {
        fprintf(stderr, "INVALID: %s\n", err.message);
        return 1;
    }

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    int32_t num_rgs = carquet_reader_num_row_groups(reader);
    int64_t total_rows = carquet_reader_num_rows(reader);
    int errors = 0;

    /* Try to read every column in every row group */
    for (int32_t rg = 0; rg < num_rgs; rg++) {
        for (int32_t c = 0; c < num_cols; c++) {
            carquet_column_reader_t* col = carquet_reader_get_column(reader, rg, c, &err);
            if (!col) {
                fprintf(stderr, "  ERROR: rg=%d col=%d (%s): %s\n",
                        rg, c, carquet_schema_column_name(schema, c), err.message);
                errors++;
                continue;
            }

            /* Read through all pages to trigger CRC checks */
            carquet_physical_type_t phys = carquet_schema_column_type(schema, c);
            int32_t elem_size = carquet_physical_type_size(phys);

            if (elem_size > 0) {
                /* Fixed-size type */
                uint8_t buf[8192];
                int64_t batch = (int64_t)(sizeof(buf) / (size_t)elem_size);
                while (carquet_column_read_batch(col, buf, batch, NULL, NULL) > 0)
                    ;
            } else {
                /* Variable-length type */
                carquet_byte_array_t buf[256];
                while (carquet_column_read_batch(col, buf, 256, NULL, NULL) > 0)
                    ;
            }

            carquet_column_reader_free(col);
        }
    }

    if (errors == 0) {
        printf("OK: %" PRId64 " rows, %d columns, %d row groups - all pages valid\n",
               total_rows, num_cols, num_rgs);
    } else {
        printf("ERRORS: %d page read failures\n", errors);
    }

    carquet_reader_close(reader);
    return errors > 0 ? 1 : 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * cmd_head / cmd_tail / cmd_sample
 *
 * All three read whole rows through read_rows() (shared with cat / export),
 * which walks each column by rows, not by level slots, so LIST / MAP columns
 * print one list per row and stay aligned with the flat columns.
 * ══════════════════════════════════════════════════════════════════════════ */

/* Print rows [offset, offset + n) — or the sorted row numbers in `picks` —
 * of every column as a table, nulls shown as "null". */
static int print_rows(carquet_reader_t* reader, int64_t offset, int64_t n,
                      const int64_t* picks) {
    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    if (n <= 0 || num_cols <= 0) return 0;

    int32_t* sel = malloc((size_t)num_cols * sizeof(int32_t));
    const char** headers = malloc((size_t)num_cols * sizeof(const char*));
    if (!sel || !headers) {
        free(sel);
        free((void*)headers);
        fprintf(stderr, "error: out of memory\n");
        return 1;
    }
    for (int32_t c = 0; c < num_cols; c++) {
        sel[c] = c;
        headers[c] = carquet_schema_column_name(schema, c);
    }
    str_matrix_t mat = {0};
    int rc = read_rows(reader, schema, sel, num_cols, offset, n, picks, "null", &mat);
    if (rc == 0) {
        print_dyn_table(headers, num_cols, (const char* const*)mat.cells, mat.num_rows);
    }
    matrix_free(&mat);
    free(sel);
    free((void*)headers);
    return rc == 0 ? 0 : 1;
}

int cmd_head(const char* path, int64_t n, const char* filter) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = open_or_die(path, &err);
    if (!reader) return 1;

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    int64_t total = carquet_reader_num_rows(reader);
    /* Cap before either path: the filtered path sizes its output from n. */
    if (n > total) n = total;
    if (n <= 0 || num_cols <= 0) {
        carquet_reader_close(reader);
        return 0;
    }
    if (filter) {
        cli_filter_storage_t fs;
        char ferr[256];
        if (cli_parse_filter(filter, schema, num_cols, &fs, ferr, sizeof(ferr)) != 0) {
            fprintf(stderr, "error: %s\n", ferr);
            carquet_reader_close(reader);
            return 1;
        }
        int32_t* sel = malloc((size_t)num_cols * sizeof(int32_t));
        const char** headers = malloc((size_t)num_cols * sizeof(const char*));
        int rc = -1;
        str_matrix_t mat = {0};
        if (sel && headers) {
            for (int32_t c = 0; c < num_cols; c++) {
                sel[c] = c;
                headers[c] = carquet_schema_column_name(schema, c);
            }
            rc = read_rows_filtered(reader, schema, sel, num_cols, 0, n,
                                    &fs, "null", &mat);
            if (rc == 0) {
                print_dyn_table(headers, num_cols,
                                (const char* const*)mat.cells, mat.num_rows);
            }
        }
        matrix_free(&mat);
        free(sel);
        free((void*)headers);
        cli_filter_free(&fs);
        carquet_reader_close(reader);
        return rc == 0 ? 0 : 1;
    }

    int rc = print_rows(reader, 0, n, NULL);
    carquet_reader_close(reader);
    return rc;
}

int cmd_tail(const char* path, int64_t n, const char* filter) {
    if (filter) {
        fprintf(stderr,
            "error: --filter is not supported with `tail` (would require\n"
            "materializing every matching row to find the last N). Use\n"
            "`cat --filter` and pipe through `tail` instead.\n");
        return 1;
    }
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = open_or_die(path, &err);
    if (!reader) return 1;

    int64_t total = carquet_reader_num_rows(reader);
    if (n > total) n = total;
    int rc = n > 0 ? print_rows(reader, total - n, n, NULL) : 0;
    carquet_reader_close(reader);
    return rc;
}

static int compare_int64(const void* a, const void* b) {
    int64_t va = *(const int64_t*)a;
    int64_t vb = *(const int64_t*)b;
    return (va > vb) - (va < vb);
}

int cmd_sample(const char* path, int64_t n, const char* filter) {
    if (filter) {
        fprintf(stderr,
            "error: --filter is not supported with `sample` (would need a\n"
            "two-pass scan to count matching rows before picking random\n"
            "indices). Use `cat --filter` and pipe through `shuf | head`.\n");
        return 1;
    }
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = open_or_die(path, &err);
    if (!reader) return 1;

    int64_t total = carquet_reader_num_rows(reader);
    if (n > total) n = total;
    if (n <= 0) {
        carquet_reader_close(reader);
        return 0;
    }

    /* n random row numbers, sorted and de-duplicated. */
    srand((unsigned)time(NULL));
    int64_t* indices = calloc((size_t)n, sizeof(int64_t));
    if (!indices) {
        carquet_reader_close(reader);
        fprintf(stderr, "error: out of memory\n");
        return 1;
    }
    for (int64_t i = 0; i < n; i++) {
        indices[i] = (int64_t)(((uint64_t)rand() * (uint64_t)rand()) % (uint64_t)total);
    }
    qsort(indices, (size_t)n, sizeof(int64_t), compare_int64);
    int64_t unique = 1;
    for (int64_t i = 1; i < n; i++) {
        if (indices[i] != indices[unique - 1]) indices[unique++] = indices[i];
    }

    int rc = print_rows(reader, 0, unique, indices);
    free(indices);
    carquet_reader_close(reader);
    return rc;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Dynamic-width tabular printer
 *
 * Generic header + cells output that auto-sizes each column to the widest
 * value (capped at MAX_COL_WIDTH). Used by `cat` and `stat` so both commands
 * produce the same clean two-space-separated layout regardless of content
 * width. `cells` is a flat array indexed as cells[row * num_cols + col];
 * a NULL entry prints empty.
 * ══════════════════════════════════════════════════════════════════════════ */

static void print_dyn_table(const char* const* headers, int32_t num_cols,
                             const char* const* cells, int64_t num_rows) {
    if (num_cols <= 0) return;

    int* widths = calloc((size_t)num_cols, sizeof(int));
    if (!widths) return;

    for (int32_t c = 0; c < num_cols; c++) {
        int len = (int)strlen(headers[c]);
        widths[c] = len < MAX_COL_WIDTH ? len : MAX_COL_WIDTH;
    }
    for (int64_t r = 0; r < num_rows; r++) {
        for (int32_t c = 0; c < num_cols; c++) {
            const char* v = cells[r * num_cols + c];
            if (!v) continue;
            int len = (int)strlen(v);
            if (len > MAX_COL_WIDTH) len = MAX_COL_WIDTH;
            if (len > widths[c]) widths[c] = len;
        }
    }

    printf("  ");
    for (int32_t c = 0; c < num_cols; c++) {
        if (c > 0) printf("  ");
        printf("%-*.*s", widths[c], widths[c], headers[c]);
    }
    printf("\n  ");
    for (int32_t c = 0; c < num_cols; c++) {
        if (c > 0) printf("  ");
        for (int w = 0; w < widths[c]; w++) putchar('-');
    }
    printf("\n");

    for (int64_t r = 0; r < num_rows; r++) {
        printf("  ");
        for (int32_t c = 0; c < num_cols; c++) {
            const char* v = cells[r * num_cols + c];
            if (!v) v = "";
            if (c > 0) printf("  ");
            printf("%-*.*s", widths[c], widths[c], v);
        }
        printf("\n");
    }
    free(widths);
}

/* ══════════════════════════════════════════════════════════════════════════
 * Shared row-extraction for cmd_cat and cmd_export
 *
 * Both commands need to read N rows starting at an offset, optionally
 * restricted to a column subset, and turn each value into a string. The
 * heavy lifting (per-column read + skip across row groups) lives here.
 * ══════════════════════════════════════════════════════════════════════════ */

/* Match `name` against the comma-separated list in `filter`. NULL filter
 * matches everything. Leading/trailing whitespace per token is tolerated. */
static bool name_in_filter(const char* name, const char* filter) {
    if (!filter) return true;
    const char* p = filter;
    size_t name_len = strlen(name);
    while (*p) {
        while (*p == ' ' || *p == '\t') p++;
        const char* comma = strchr(p, ',');
        size_t tok_len = comma ? (size_t)(comma - p) : strlen(p);
        while (tok_len > 0 && (p[tok_len - 1] == ' ' || p[tok_len - 1] == '\t'))
            tok_len--;
        if (tok_len == name_len && strncmp(p, name, tok_len) == 0)
            return true;
        p = comma ? comma + 1 : p + strlen(p);
    }
    return false;
}

/* Resolve the column filter into a list of column indices. Returns the
 * number of selected columns, or -1 if a name didn't match the schema.
 * On success, *out is a malloc'd array the caller must free. */
static int32_t resolve_columns(const carquet_schema_t* schema,
                                int32_t num_cols, const char* filter,
                                int32_t** out) {
    int32_t* sel = malloc((size_t)num_cols * sizeof(int32_t));
    if (!sel) return -1;
    int32_t n = 0;
    for (int32_t c = 0; c < num_cols; c++) {
        const char* nm = carquet_schema_column_name(schema, c);
        if (name_in_filter(nm, filter)) {
            sel[n++] = c;
        }
    }
    if (filter && n == 0) {
        free(sel);
        return -1;
    }
    *out = sel;
    return n;
}

/* String matrix used to hold formatted values before display/export. */
static void matrix_free(str_matrix_t* m) {
    if (m->cells) {
        int64_t total = m->num_rows * m->num_cols;
        for (int64_t i = 0; i < total; i++) free(m->cells[i]);
        free(m->cells);
    }
    m->cells = NULL;
    m->num_rows = 0;
}

/* Allocate an empty rows x cols matrix; -1 on overflow or out of memory. */
static int matrix_alloc(str_matrix_t* m, int64_t rows, int32_t cols) {
    m->num_cols = cols;
    m->num_rows = 0;
    m->cells = NULL;
    if (rows <= 0 || cols <= 0) return 0;
    if ((uint64_t)rows > SIZE_MAX / sizeof(char*) / (uint64_t)cols) return -1;
    m->cells = calloc((size_t)rows * (size_t)cols, sizeof(char*));
    return m->cells ? 0 : -1;
}

/* Bounded string builder: appends until full, then ends the text in "...". */
typedef struct {
    char* buf;
    size_t cap;
    size_t len;
    bool full;
} cell_text_t;

static void cell_puts(cell_text_t* t, const char* s) {
    if (t->full || t->cap < 4) return;
    size_t n = strlen(s);
    if (t->len + n < t->cap) {
        memcpy(t->buf + t->len, s, n + 1);
        t->len += n;
        return;
    }
    if (t->len > t->cap - 4) t->len = t->cap - 4;
    memcpy(t->buf + t->len, "...", 4);
    t->len += 3;
    t->full = true;
}

/*
 * Row cursor over one leaf column, across row groups.
 *
 * Column readers return level slots: a REPEATED leaf has several slots per
 * row (a row starts at rep == 0) and carquet_column_skip() skips slots. The
 * cursor reads slots in chunks and groups them into rows, so every column
 * advances by the same rows. Row groups are delimited by their row count;
 * a row group whose column cannot be read prints "?" for its rows.
 */
typedef struct {
    carquet_reader_t* reader;
    int32_t col;
    carquet_physical_type_t phys;
    const carquet_logical_type_t* lt;
    int32_t tl;
    size_t stride;              /* bytes per decoded value; 0 = unprintable type */
    int16_t max_def, max_rep;
    int16_t* list_def;          /* [1..max_rep]: def level once list j has an entry */
    int32_t num_rgs, rg;        /* current row group (-1 before the first) */
    int64_t rg_start, rg_end;   /* its rows: [rg_start, rg_end) */
    int64_t pos_row;            /* row the column reader is positioned at */
    carquet_column_reader_t* cr;
    bool broken;                /* current row group unreadable */
    int64_t chunk;              /* slots per refill */
    uint8_t* vals;
    int16_t* def;
    int16_t* rep;
    int64_t n, pos, dense;      /* buffered slots, next slot, its dense index */
} col_cursor_t;

static void cursor_free(col_cursor_t* c) {
    if (c->cr) carquet_column_reader_free(c->cr);
    free(c->list_def);
    free(c->vals);
    free(c->def);
    free(c->rep);
    memset(c, 0, sizeof(*c));
}

static int cursor_init(col_cursor_t* c, carquet_reader_t* reader,
                       const carquet_schema_t* schema, int32_t col) {
    memset(c, 0, sizeof(*c));
    c->reader = reader;
    c->col = col;
    c->rg = -1;
    c->num_rgs = carquet_reader_num_row_groups(reader);
    c->phys = carquet_schema_column_type(schema, col);
    int32_t leaf = schema->leaf_indices[col];
    const carquet_schema_node_t* node = carquet_schema_get_element(schema, leaf);
    c->lt = carquet_schema_node_logical_type(node);
    c->tl = carquet_schema_node_type_length(node);
    c->max_def = carquet_schema_max_def_level(schema, col);
    c->max_rep = carquet_schema_max_rep_level(schema, col);

    if (c->phys == CARQUET_PHYSICAL_BYTE_ARRAY) {
        c->stride = sizeof(carquet_byte_array_t);
    } else if (c->phys == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY) {
        c->stride = c->tl > 0 ? (size_t)c->tl : 0;
    } else {
        int32_t es = carquet_physical_type_size(c->phys);
        c->stride = es > 0 ? (size_t)es : 0;
    }
    if (c->stride == 0) return 0;   /* malformed type: every row prints "?" */

    /* Definition level reached by each repeated ancestor, root to leaf. */
    if (c->max_rep > 0) {
        c->list_def = calloc((size_t)c->max_rep + 1, sizeof(int16_t));
        int32_t depth = 0;
        for (int32_t e = leaf; e > 0; e = schema->parent_indices[e]) depth++;
        int32_t* path = malloc((size_t)depth * sizeof(int32_t));
        if (!c->list_def || !path) { free(path); return -1; }
        int32_t k = depth;
        for (int32_t e = leaf; e > 0; e = schema->parent_indices[e]) path[--k] = e;
        int16_t d = 0, r = 0;
        for (int32_t i = 0; i < depth; i++) {
            const parquet_schema_element_t* el = &schema->elements[path[i]];
            if (el->repetition_type == CARQUET_REPETITION_OPTIONAL) d++;
            if (el->repetition_type == CARQUET_REPETITION_REPEATED) {
                d++;
                if (r < c->max_rep) c->list_def[++r] = d;
            }
        }
        free(path);
    }

    /* Keep the value buffer small for wide FLBA (type_length is file-set). */
    c->chunk = (int64_t)((4u << 20) / c->stride);
    if (c->chunk > 1024) c->chunk = 1024;
    if (c->chunk < 1) c->chunk = 1;
    c->vals = malloc((size_t)c->chunk * c->stride);
    if (c->max_def > 0) c->def = malloc((size_t)c->chunk * sizeof(int16_t));
    if (c->max_rep > 0) c->rep = malloc((size_t)c->chunk * sizeof(int16_t));
    if (!c->vals || (c->max_def > 0 && !c->def) || (c->max_rep > 0 && !c->rep)) return -1;
    return 0;
}

/* Open the row group holding `row` (forward only: row >= rg_start). */
static void cursor_enter_row_group(col_cursor_t* c, int64_t row) {
    if (c->rg >= 0 && row < c->rg_end) return;
    if (c->cr) { carquet_column_reader_free(c->cr); c->cr = NULL; }
    c->n = c->pos = c->dense = 0;
    int64_t start = c->rg >= 0 ? c->rg_end : 0;
    while (++c->rg < c->num_rgs) {
        carquet_row_group_metadata_t rgm;
        int64_t rows = 0;
        if (carquet_reader_row_group_metadata(c->reader, c->rg, &rgm) == CARQUET_OK &&
            rgm.num_rows > 0) {
            rows = rgm.num_rows;
        }
        c->rg_start = start;
        c->rg_end = start + rows;
        if (row < c->rg_end) break;
        start = c->rg_end;
    }
    c->pos_row = c->rg_start;
    c->broken = c->rg >= c->num_rgs || c->stride == 0;
    if (!c->broken) {
        carquet_error_t err = CARQUET_ERROR_INIT;
        c->cr = carquet_reader_get_column(c->reader, c->rg, c->col, &err);
        c->broken = c->cr == NULL;
    }
}

/* Make at least one slot available; false at the end of the chunk. */
static bool cursor_fill(col_cursor_t* c) {
    if (c->pos < c->n) return true;
    if (c->broken || !c->cr) return false;
    int64_t got = carquet_column_read_batch(c->cr, c->vals, c->chunk, c->def, c->rep);
    c->n = got > 0 ? got : 0;
    c->pos = c->dense = 0;
    return c->n > 0;
}

/* Consume one slot; returns its value when present, NULL for a null slot. */
static const void* cursor_take_slot(col_cursor_t* c, int16_t* d, int16_t* r) {
    *d = c->def ? c->def[c->pos] : c->max_def;
    *r = c->rep ? c->rep[c->pos] : 0;
    c->pos++;
    if (*d < c->max_def) return NULL;
    return c->vals + (size_t)(c->dense++) * c->stride;
}

/* Consume the next row of the current row group and, when out != NULL,
 * format it. A repeated leaf prints as nested lists ("[1, null, 3]", "[]");
 * null_text is used for a null row. */
static void cursor_take_row(col_cursor_t* c, cell_text_t* out, const char* null_text) {
    if (!cursor_fill(c)) {
        c->broken = true;   /* fewer values than rows: stop trusting this chunk */
        if (out) cell_puts(out, "?");
        return;
    }
    int16_t open = 0;
    for (int64_t t = 0; ; t++) {
        int16_t d, r;
        const void* vp = cursor_take_slot(c, &d, &r);
        if (out) {
            char vbuf[MAX_VALUE_BUF];
            if (c->max_rep == 0) {
                cell_puts(out, vp ? cli_format_value(c->phys, vp, c->tl, c->lt,
                                                     vbuf, sizeof(vbuf))
                                  : null_text);
            } else {
                if (t == 0) r = 0;
                if (r > c->max_rep) r = c->max_rep;
                while (open > r) { cell_puts(out, "]"); open--; }
                while (open < r) { cell_puts(out, "["); open++; }
                if (t > 0) cell_puts(out, ", ");
                /* Lists r+1..m have an entry in this slot. */
                int16_t m = r;
                while (m < c->max_rep && d >= c->list_def[m + 1]) m++;
                while (open < m) { cell_puts(out, "["); open++; }
                if (m < c->max_rep) {
                    /* List m+1 has no entry: empty if its parent is defined. */
                    if (d == c->list_def[m + 1] - 1) cell_puts(out, "[]");
                    else cell_puts(out, m == 0 ? null_text : "null");
                } else {
                    cell_puts(out, vp ? cli_format_value(c->phys, vp, c->tl, c->lt,
                                                         vbuf, sizeof(vbuf))
                                      : "null");
                }
            }
        }
        if (c->max_rep == 0) break;
        if (!cursor_fill(c) || c->rep[c->pos] == 0) break;
    }
    while (out && open > 0) { cell_puts(out, "]"); open--; }
}

/* Format row `row` (rows are requested in increasing order) into out. */
static void cursor_read_row(col_cursor_t* c, int64_t row, cell_text_t* out,
                            const char* null_text) {
    cursor_enter_row_group(c, row);
    if (!c->broken && c->pos_row < row) {
        if (c->max_rep == 0) {
            /* Flat: one slot per row; drop buffered slots, skip the rest. */
            while (c->pos_row < row && c->pos < c->n) {
                int16_t d, r;
                (void)cursor_take_slot(c, &d, &r);
                c->pos_row++;
            }
            if (c->pos_row < row) {
                int64_t want = row - c->pos_row;
                if (carquet_column_skip(c->cr, want) != want) c->broken = true;
            }
        } else {
            while (c->pos_row < row && !c->broken) {
                cursor_take_row(c, NULL, NULL);
                c->pos_row++;
            }
        }
        c->pos_row = row;
    }
    if (c->broken) cell_puts(out, "?");
    else cursor_take_row(c, out, null_text);
    c->pos_row++;
}

/* Fill `out` with `limit` rows of the selected columns: rows offset.. or,
 * when `picks` is non-NULL, the sorted row numbers picks[0..limit). Nulls
 * print as null_text. Returns 0, or -1 on allocation failure. */
static int read_rows(carquet_reader_t* reader,
                     const carquet_schema_t* schema,
                     const int32_t* col_indices, int32_t num_sel_cols,
                     int64_t offset, int64_t limit, const int64_t* picks,
                     const char* null_text, str_matrix_t* out) {
    int64_t total = carquet_reader_num_rows(reader);
    if (!picks) {
        if (offset < 0) offset = 0;
        if (limit > total - offset) limit = total - offset;
    }
    if (matrix_alloc(out, limit, num_sel_cols) != 0) {
        fprintf(stderr, "error: out of memory\n");
        return -1;
    }
    if (!out->cells) return 0;

    for (int32_t i = 0; i < num_sel_cols; i++) {
        col_cursor_t cur;
        if (cursor_init(&cur, reader, schema, col_indices[i]) != 0) {
            cursor_free(&cur);
            fprintf(stderr, "error: out of memory\n");
            out->num_rows = limit;
            return -1;
        }
        for (int64_t r = 0; r < limit; r++) {
            int64_t row = picks ? picks[r] : offset + r;
            if (row < 0 || row >= total) break;
            char text[MAX_VALUE_BUF];
            cell_text_t ct = {text, sizeof(text), 0, false};
            text[0] = '\0';
            cursor_read_row(&cur, row, &ct, null_text);
            out->cells[r * num_sel_cols + i] = carquet_heap_strdup(text);
        }
        cursor_free(&cur);
    }
    out->num_rows = limit;
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * Filter expression parser
 *
 * Grammar (case-insensitive keywords):
 *   filter   := clause (AND clause)*
 *   clause   := name op value
 *             | name 'IS' 'NULL'
 *             | name 'IS' 'NOT' 'NULL'
 *   op       := '=' | '==' | '!=' | '<>' | '<' | '<=' | '>' | '>='
 *   value    := signed_number | quoted_string | TRUE | FALSE
 *
 * Strings use single quotes; embedded quotes are doubled ('it''s'). Each
 * value is converted to the column's physical type (INT32/64, FLOAT/DOUBLE,
 * BOOLEAN, BYTE_ARRAY). FIXED_LEN_BYTE_ARRAY / FLOAT16 / INT96 columns are
 * not supported via the CLI grammar — they need raw bytes the parser would
 * have to encode, which is out of scope; use the library API for those.
 * ══════════════════════════════════════════════════════════════════════════ */

static void cli_filter_free(cli_filter_storage_t* s) {
    if (!s) return;
    if (s->blobs) {
        for (int32_t i = 0; i < s->num_blobs; i++) free(s->blobs[i]);
        free(s->blobs);
    }
    free(s->clauses);
    memset(s, 0, sizeof(*s));
}

static int cli_filter_grow(cli_filter_storage_t* s) {
    int32_t new_cap = s->capacity > 0 ? s->capacity * 2 : 4;
    carquet_filter_clause_t* nc = realloc(s->clauses,
        (size_t)new_cap * sizeof(carquet_filter_clause_t));
    if (!nc) return -1;
    s->clauses = nc;   /* realloc may have moved it: publish before the next one can fail */
    uint8_t** nb = realloc(s->blobs, (size_t)new_cap * sizeof(uint8_t*));
    if (!nb) return -1;
    s->blobs = nb;
    s->capacity = new_cap;
    return 0;
}

static void filter_skip_ws(const char** p) {
    while (**p == ' ' || **p == '\t' || **p == '\n' || **p == '\r') (*p)++;
}

/* Compare a literal keyword case-insensitively; on match, advance *p
 * past the keyword (caller still needs to require trailing whitespace
 * or end-of-input). */
static bool filter_match_kw(const char** p, const char* kw) {
    const char* s = *p;
    size_t n = strlen(kw);
    for (size_t i = 0; i < n; i++) {
        char a = s[i];
        char b = kw[i];
        if (a >= 'a' && a <= 'z') a = (char)(a - 'a' + 'A');
        if (b >= 'a' && b <= 'z') b = (char)(b - 'a' + 'A');
        if (a != b) return false;
    }
    /* Must be followed by a delimiter (not a continuing identifier). */
    char c = s[n];
    if (c && (c >= 'a' && c <= 'z')) return false;
    if (c && (c >= 'A' && c <= 'Z')) return false;
    if (c && c == '_') return false;
    if (c && (c >= '0' && c <= '9')) return false;
    *p = s + n;
    return true;
}

static bool filter_is_ident_char(char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z')
        || (c >= '0' && c <= '9') || c == '_' || c == '-' || c == '.';
}

static int filter_parse_ident(const char** p, char* out, size_t cap) {
    filter_skip_ws(p);
    size_t n = 0;
    while (filter_is_ident_char(**p)) {
        if (n + 1 >= cap) return -1;
        out[n++] = **p;
        (*p)++;
    }
    if (n == 0) return -1;
    out[n] = 0;
    return 0;
}

static int filter_lookup_column(const carquet_schema_t* schema,
                                int32_t num_cols, const char* name) {
    for (int32_t c = 0; c < num_cols; c++) {
        const char* cn = carquet_schema_column_name(schema, c);
        if (cn && strcmp(cn, name) == 0) return c;
    }
    return -1;
}

static int filter_parse_op(const char** p, carquet_filter_op_t* out) {
    filter_skip_ws(p);
    const char* s = *p;
    if (s[0] == '!' && s[1] == '=') { *out = CARQUET_FILTER_NE; *p = s + 2; return 0; }
    if (s[0] == '<' && s[1] == '>') { *out = CARQUET_FILTER_NE; *p = s + 2; return 0; }
    if (s[0] == '<' && s[1] == '=') { *out = CARQUET_FILTER_LE; *p = s + 2; return 0; }
    if (s[0] == '>' && s[1] == '=') { *out = CARQUET_FILTER_GE; *p = s + 2; return 0; }
    if (s[0] == '=' && s[1] == '=') { *out = CARQUET_FILTER_EQ; *p = s + 2; return 0; }
    if (s[0] == '=') { *out = CARQUET_FILTER_EQ; *p = s + 1; return 0; }
    if (s[0] == '<') { *out = CARQUET_FILTER_LT; *p = s + 1; return 0; }
    if (s[0] == '>') { *out = CARQUET_FILTER_GT; *p = s + 1; return 0; }
    return -1;
}

/* Decode a single-quoted string literal. Returns a freshly malloc'd
 * buffer holding the unescaped bytes; sets *len_out. */
static uint8_t* filter_parse_string(const char** p, int32_t* len_out,
                                    char* err, size_t errsz) {
    if (**p != '\'') {
        snprintf(err, errsz, "expected string literal at: %.20s", *p);
        return NULL;
    }
    (*p)++;
    size_t cap = 16;
    uint8_t* buf = malloc(cap);
    if (!buf) return NULL;
    size_t n = 0;
    while (**p) {
        if (**p == '\'') {
            if ((*p)[1] == '\'') {
                if (n + 1 > cap) {
                    cap *= 2;
                    uint8_t* nb = realloc(buf, cap);
                    if (!nb) { free(buf); return NULL; }
                    buf = nb;
                }
                buf[n++] = '\'';
                *p += 2;
                continue;
            }
            (*p)++;
            *len_out = (int32_t)n;
            return buf;
        }
        if (n + 1 > cap) {
            cap *= 2;
            uint8_t* nb = realloc(buf, cap);
            if (!nb) { free(buf); return NULL; }
            buf = nb;
        }
        buf[n++] = (uint8_t)**p;
        (*p)++;
    }
    snprintf(err, errsz, "unterminated string literal");
    free(buf);
    return NULL;
}

/* Convert a parsed literal value into the column's native binary format
 * and stash it in a freshly malloc'd buffer of the right size. Stores
 * the resulting (ptr, size) on `clause`. */
static int filter_encode_value(const carquet_schema_t* schema,
                               int32_t col_idx, const char* val_start,
                               const char* val_end,
                               carquet_filter_clause_t* clause,
                               uint8_t** out_blob,
                               char* err, size_t errsz) {
    carquet_physical_type_t phys = carquet_schema_column_type(schema, col_idx);
    char buf[128];
    size_t vlen = (size_t)(val_end - val_start);
    if (vlen >= sizeof(buf)) {
        snprintf(err, errsz, "numeric literal too long");
        return -1;
    }
    memcpy(buf, val_start, vlen);
    buf[vlen] = 0;

    switch (phys) {
        case CARQUET_PHYSICAL_BOOLEAN: {
            uint8_t* p = malloc(1);
            if (!p) return -1;
            if (strcmp(buf, "true") == 0 || strcmp(buf, "TRUE") == 0 ||
                strcmp(buf, "1") == 0) {
                p[0] = 1;
            } else if (strcmp(buf, "false") == 0 || strcmp(buf, "FALSE") == 0 ||
                       strcmp(buf, "0") == 0) {
                p[0] = 0;
            } else {
                free(p);
                snprintf(err, errsz, "boolean expects true/false/0/1, got '%s'", buf);
                return -1;
            }
            clause->value = p;
            clause->value_size = 1;
            *out_blob = p;
            return 0;
        }
        case CARQUET_PHYSICAL_INT32: {
            char* end;
            long long v = strtoll(buf, &end, 10);
            if (*end != 0 || end == buf) {
                snprintf(err, errsz, "expected INT32, got '%s'", buf);
                return -1;
            }
            int32_t v32 = (int32_t)v;
            uint8_t* p = malloc(4);
            if (!p) return -1;
            memcpy(p, &v32, 4);
            clause->value = p;
            clause->value_size = 4;
            *out_blob = p;
            return 0;
        }
        case CARQUET_PHYSICAL_INT64: {
            char* end;
            long long v = strtoll(buf, &end, 10);
            if (*end != 0 || end == buf) {
                snprintf(err, errsz, "expected INT64, got '%s'", buf);
                return -1;
            }
            int64_t v64 = (int64_t)v;
            uint8_t* p = malloc(8);
            if (!p) return -1;
            memcpy(p, &v64, 8);
            clause->value = p;
            clause->value_size = 8;
            *out_blob = p;
            return 0;
        }
        case CARQUET_PHYSICAL_FLOAT: {
            char* end;
            double v = strtod(buf, &end);
            if (*end != 0 || end == buf) {
                snprintf(err, errsz, "expected FLOAT, got '%s'", buf);
                return -1;
            }
            float vf = (float)v;
            uint8_t* p = malloc(4);
            if (!p) return -1;
            memcpy(p, &vf, 4);
            clause->value = p;
            clause->value_size = 4;
            *out_blob = p;
            return 0;
        }
        case CARQUET_PHYSICAL_DOUBLE: {
            char* end;
            double v = strtod(buf, &end);
            if (*end != 0 || end == buf) {
                snprintf(err, errsz, "expected DOUBLE, got '%s'", buf);
                return -1;
            }
            uint8_t* p = malloc(8);
            if (!p) return -1;
            memcpy(p, &v, 8);
            clause->value = p;
            clause->value_size = 8;
            *out_blob = p;
            return 0;
        }
        default:
            snprintf(err, errsz,
                "filter literal type unsupported for column physical type %d "
                "(use a string literal for BYTE_ARRAY)", (int)phys);
            return -1;
    }
}

/* Parse a single clause and append it to the storage. */
static int filter_parse_clause(const char** p, const carquet_schema_t* schema,
                               int32_t num_cols, cli_filter_storage_t* s,
                               char* err, size_t errsz) {
    if (s->count >= s->capacity && cli_filter_grow(s) != 0) {
        snprintf(err, errsz, "out of memory");
        return -1;
    }
    carquet_filter_clause_t* c = &s->clauses[s->count];
    memset(c, 0, sizeof(*c));
    s->blobs[s->count] = NULL;
    /* The slot is owned from here on, even if this clause fails to parse
     * after its value was allocated. */
    if (s->num_blobs <= s->count) s->num_blobs = s->count + 1;

    char name[128];
    if (filter_parse_ident(p, name, sizeof(name)) != 0) {
        snprintf(err, errsz, "expected column name at: %.20s", *p);
        return -1;
    }
    int32_t col = filter_lookup_column(schema, num_cols, name);
    if (col < 0) {
        snprintf(err, errsz, "unknown column '%s'", name);
        return -1;
    }
    c->column_index = col;

    /* IS [NOT] NULL */
    filter_skip_ws(p);
    const char* save = *p;
    if (filter_match_kw(p, "IS")) {
        filter_skip_ws(p);
        if (filter_match_kw(p, "NOT")) {
            filter_skip_ws(p);
            if (!filter_match_kw(p, "NULL")) {
                snprintf(err, errsz, "expected NULL after IS NOT");
                return -1;
            }
            c->op = CARQUET_FILTER_IS_NOT_NULL;
        } else if (filter_match_kw(p, "NULL")) {
            c->op = CARQUET_FILTER_IS_NULL;
        } else {
            snprintf(err, errsz, "expected NULL after IS");
            return -1;
        }
        s->count++;
        return 0;
    }
    *p = save;

    /* op value */
    if (filter_parse_op(p, &c->op) != 0) {
        snprintf(err, errsz, "expected comparison operator at: %.20s", *p);
        return -1;
    }

    filter_skip_ws(p);
    carquet_physical_type_t phys = carquet_schema_column_type(schema, col);
    if (phys == CARQUET_PHYSICAL_BYTE_ARRAY) {
        if (**p != '\'') {
            snprintf(err, errsz,
                "expected single-quoted string for BYTE_ARRAY column '%s'",
                name);
            return -1;
        }
        int32_t slen = 0;
        uint8_t* sval = filter_parse_string(p, &slen, err, errsz);
        if (!sval) return -1;
        c->value = sval;
        c->value_size = slen;
        s->blobs[s->count] = sval;
    } else {
        const char* start = *p;
        /* Allow leading sign + digits, dot, exponent. */
        if (**p == '+' || **p == '-') (*p)++;
        while ((**p >= '0' && **p <= '9') || **p == '.' ||
               **p == 'e' || **p == 'E' || **p == '+' || **p == '-' ||
               (**p >= 'a' && **p <= 'z') || (**p >= 'A' && **p <= 'Z')) {
            (*p)++;
        }
        if (*p == start) {
            snprintf(err, errsz, "expected literal at: %.20s", start);
            return -1;
        }
        uint8_t* blob = NULL;
        if (filter_encode_value(schema, col, start, *p, c, &blob,
                                err, errsz) != 0) {
            return -1;
        }
        s->blobs[s->count] = blob;
    }
    s->count++;
    return 0;
}

/* Parse the full expression and populate storage. Returns 0 on success. */
static int cli_parse_filter(const char* expr,
                            const carquet_schema_t* schema,
                            int32_t num_cols,
                            cli_filter_storage_t* out,
                            char* err, size_t errsz) {
    memset(out, 0, sizeof(*out));
    const char* p = expr;
    for (;;) {
        if (filter_parse_clause(&p, schema, num_cols, out, err, errsz) != 0) {
            cli_filter_free(out);
            return -1;
        }
        filter_skip_ws(&p);
        if (*p == 0) return 0;
        if (!filter_match_kw(&p, "AND")) {
            snprintf(err, errsz, "expected AND or end of expression at: %.20s", p);
            cli_filter_free(out);
            return -1;
        }
    }
}

/* ══════════════════════════════════════════════════════════════════════════
 * Filtered read path — uses the batch reader (so set_page_filter works).
 *
 * Returns -1 on hard error; the file's batch-reader error code is mapped
 * to a stderr message. Skips `offset` matching rows and emits up to
 * `limit` matching rows into the matrix.
 * ══════════════════════════════════════════════════════════════════════════ */

/* Format row r of a single-level list column of a batch as "[a, null, b]". */
static void format_batch_list_cell(const carquet_row_batch_t* batch, int32_t c,
                                   int64_t r, carquet_physical_type_t phys,
                                   int32_t tl, const carquet_logical_type_t* lt,
                                   const char* null_text, char* text, size_t cap) {
    cell_text_t ct = {text, cap, 0, false};
    text[0] = '\0';
    const int32_t* offsets;
    int64_t num_lists, num_values;
    const void* values;
    const uint8_t* vvalid;
    const uint8_t* lvalid;
    if (carquet_row_batch_column_list(batch, c, &offsets, &num_lists, &values,
                                      &vvalid, &num_values, &lvalid) != CARQUET_OK ||
        r >= num_lists) {
        cell_puts(&ct, "?");
        return;
    }
    if (lvalid && !(lvalid[r / 8] & (1u << (r % 8)))) {
        cell_puts(&ct, null_text);
        return;
    }
    size_t stride = phys == CARQUET_PHYSICAL_BYTE_ARRAY ? sizeof(carquet_byte_array_t)
                  : phys == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY ? (size_t)(tl > 0 ? tl : 0)
                  : (size_t)carquet_physical_type_size(phys);
    cell_puts(&ct, "[");
    for (int32_t k = offsets[r]; k < offsets[r + 1] && k < num_values; k++) {
        if (k > offsets[r]) cell_puts(&ct, ", ");
        if ((vvalid && !(vvalid[k / 8] & (1u << (k % 8)))) || stride == 0) {
            cell_puts(&ct, "null");
        } else {
            char vbuf[MAX_VALUE_BUF];
            cell_puts(&ct, cli_format_value(phys, (const uint8_t*)values + (size_t)k * stride,
                                            tl, lt, vbuf, sizeof(vbuf)));
        }
    }
    cell_puts(&ct, "]");
}

static int read_rows_filtered(carquet_reader_t* reader,
                              const carquet_schema_t* schema,
                              const int32_t* col_indices, int32_t num_sel_cols,
                              int64_t offset, int64_t limit,
                              const cli_filter_storage_t* filter,
                              const char* null_text, str_matrix_t* out) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    /* limit * num_sel_cols cells: overflow-checked (head -n is user-set). */
    if (matrix_alloc(out, limit, num_sel_cols) != 0) {
        fprintf(stderr, "error: out of memory\n");
        return -1;
    }
    if (!out->cells) return 0;

    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 4096;
    cfg.column_indices = col_indices;
    cfg.num_columns = num_sel_cols;

    carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &err);
    if (!br) {
        fprintf(stderr, "error: %s\n", err.message);
        return -1;
    }
    if (filter && filter->count > 0) {
        carquet_status_t st = carquet_batch_reader_set_page_filter(
            br, filter->clauses, filter->count);
        if (st != CARQUET_OK) {
            fprintf(stderr,
                "error: invalid filter (status %d)\n", (int)st);
            carquet_batch_reader_free(br);
            return -1;
        }
    }

    /* Cache type metadata for cell formatting. */
    carquet_physical_type_t* phys = malloc((size_t)num_sel_cols * sizeof(*phys));
    const carquet_logical_type_t** lts = malloc((size_t)num_sel_cols * sizeof(*lts));
    int32_t* tls = malloc((size_t)num_sel_cols * sizeof(*tls));
    int16_t* max_defs = malloc((size_t)num_sel_cols * sizeof(*max_defs));
    int16_t* max_reps = malloc((size_t)num_sel_cols * sizeof(*max_reps));
    if (!phys || !lts || !tls || !max_defs || !max_reps) {
        free(phys); free(lts); free(tls); free(max_defs); free(max_reps);
        carquet_batch_reader_free(br);
        return -1;
    }
    for (int32_t c = 0; c < num_sel_cols; c++) {
        int32_t file_col = col_indices[c];
        phys[c] = carquet_schema_column_type(schema, file_col);
        const carquet_schema_node_t* node = carquet_schema_get_element(schema,
            schema->leaf_indices[file_col]);
        lts[c] = carquet_schema_node_logical_type(node);
        tls[c] = carquet_schema_node_type_length(node);
        max_defs[c] = carquet_schema_max_def_level(schema, file_col);
        max_reps[c] = carquet_schema_max_rep_level(schema, file_col);
    }

    int64_t skipped = 0;
    int64_t produced = 0;
    int rc = 0;
    carquet_row_batch_t* batch = NULL;
    while (produced < limit) {
        carquet_status_t st = carquet_batch_reader_next(br, &batch);
        if (st != CARQUET_OK || !batch) {
            if (st != CARQUET_OK && st != CARQUET_ERROR_END_OF_DATA) {
                if (st == CARQUET_ERROR_PAGE_INDEX_REQUIRED) {
                    fprintf(stderr,
                        "error: filter requires a page index but the file\n"
                        "has none for at least one referenced column.\n"
                        "Re-write the file with write_page_index = true.\n");
                } else {
                    fprintf(stderr,
                        "error: filtered read failed: %s\n",
                        carquet_status_string(st));
                }
                rc = -1;
            }
            break;
        }
        int64_t batch_rows = carquet_row_batch_num_rows(batch);
        for (int64_t r = 0; r < batch_rows && produced < limit; r++) {
            if (skipped < offset) { skipped++; continue; }
            for (int32_t c = 0; c < num_sel_cols; c++) {
                if (max_reps[c] > 0) {
                    /* Repeated leaf: the batch holds it as an Arrow list. */
                    char text[MAX_VALUE_BUF];
                    format_batch_list_cell(batch, c, r, phys[c], tls[c], lts[c],
                                           null_text, text, sizeof(text));
                    out->cells[produced * num_sel_cols + c] = carquet_heap_strdup(text);
                    continue;
                }
                const void* data;
                const uint8_t* nb;
                int64_t n;
                if (carquet_row_batch_column(batch, c, &data, &nb, &n)
                    != CARQUET_OK) continue;
                char vbuf[MAX_VALUE_BUF];
                const char* cell;
                bool is_null = false;
                if (nb && max_defs[c] > 0) {
                    is_null = (nb[r / 8] & (1u << (r % 8))) == 0;
                }
                if (is_null) {
                    cell = null_text;
                } else {
                    const void* vp;
                    int32_t tl = tls[c];
                    if (phys[c] == CARQUET_PHYSICAL_BYTE_ARRAY) {
                        vp = &((const carquet_byte_array_t*)data)[r];
                    } else if (phys[c] == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY) {
                        vp = (const uint8_t*)data + (size_t)r * (size_t)tl;
                    } else {
                        int32_t es = carquet_physical_type_size(phys[c]);
                        vp = (const uint8_t*)data + (size_t)r * (size_t)es;
                    }
                    cli_format_value(phys[c], vp, tl, lts[c],
                                     vbuf, sizeof(vbuf));
                    cell = vbuf;
                }
                out->cells[produced * num_sel_cols + c] =
                    carquet_heap_strdup(cell);
            }
            produced++;
        }
        carquet_row_batch_free(batch);
        batch = NULL;
    }

    out->num_rows = produced;
    free(phys); free(lts); free(tls); free(max_defs); free(max_reps);
    carquet_batch_reader_free(br);
    return rc;
}

/* Count matching rows under a filter via the batch reader. */
static int64_t count_rows_filtered(carquet_reader_t* reader,
                                   const cli_filter_storage_t* filter) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_batch_reader_config_t cfg;
    carquet_batch_reader_config_init(&cfg);
    cfg.batch_size = 65536;
    /* Project just column 0 to minimize materialization cost — we only
     * care about row counts. */
    int32_t one = 0;
    cfg.column_indices = &one;
    cfg.num_columns = 1;

    carquet_batch_reader_t* br = carquet_batch_reader_create(reader, &cfg, &err);
    if (!br) {
        fprintf(stderr, "error: %s\n", err.message);
        return -1;
    }
    if (filter && filter->count > 0) {
        carquet_status_t st = carquet_batch_reader_set_page_filter(
            br, filter->clauses, filter->count);
        if (st != CARQUET_OK) {
            fprintf(stderr, "error: invalid filter (status %d)\n", (int)st);
            carquet_batch_reader_free(br);
            return -1;
        }
    }
    int64_t total = 0;
    carquet_row_batch_t* batch = NULL;
    carquet_status_t st;
    while ((st = carquet_batch_reader_next(br, &batch)) == CARQUET_OK && batch) {
        total += carquet_row_batch_num_rows(batch);
        carquet_row_batch_free(batch);
        batch = NULL;
    }
    if (st != CARQUET_OK && st != CARQUET_ERROR_END_OF_DATA) {
        if (st == CARQUET_ERROR_PAGE_INDEX_REQUIRED) {
            fprintf(stderr,
                "error: filter requires a page index but the file has\n"
                "none for at least one referenced column. Re-write the\n"
                "file with write_page_index = true.\n");
        } else {
            fprintf(stderr,
                "error: filtered read failed: %s\n",
                carquet_status_string(st));
        }
        carquet_batch_reader_free(br);
        return -1;
    }
    carquet_batch_reader_free(br);
    return total;
}

/* ══════════════════════════════════════════════════════════════════════════
 * cmd_cat — print rows with optional slicing and column filter
 * ══════════════════════════════════════════════════════════════════════════ */

int cmd_cat(const char* path, const row_select_opts_t* opts) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = open_or_die(path, &err);
    if (!reader) return 1;

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    int64_t total = carquet_reader_num_rows(reader);

    int64_t offset = opts->offset < 0 ? 0 : opts->offset;
    if (offset > total) offset = total;
    int64_t limit = opts->limit < 0 ? (total - offset) : opts->limit;
    if (limit > total - offset) limit = total - offset;

    int32_t* sel = NULL;
    int32_t num_sel = resolve_columns(schema, num_cols, opts->columns, &sel);
    if (num_sel < 0) {
        fprintf(stderr, "error: no columns matched filter '%s'\n",
                opts->columns ? opts->columns : "");
        free(sel);
        carquet_reader_close(reader);
        return 1;
    }
    if (limit <= 0 || num_sel == 0) {
        free(sel);
        carquet_reader_close(reader);
        return 0;
    }

    cli_filter_storage_t fs;
    bool has_filter = false;
    if (opts->filter) {
        char ferr[256];
        if (cli_parse_filter(opts->filter, schema, num_cols, &fs,
                             ferr, sizeof(ferr)) != 0) {
            fprintf(stderr, "error: %s\n", ferr);
            free(sel);
            carquet_reader_close(reader);
            return 1;
        }
        has_filter = true;
    }

    str_matrix_t mat = {0};
    int rc;
    if (has_filter) {
        rc = read_rows_filtered(reader, schema, sel, num_sel, offset, limit,
                                &fs, "", &mat);
    } else {
        rc = read_rows(reader, schema, sel, num_sel, offset, limit, NULL, "", &mat);
    }
    if (rc != 0) {
        if (has_filter) cli_filter_free(&fs);
        matrix_free(&mat);
        free(sel);
        carquet_reader_close(reader);
        return 1;
    }

    const char** headers = malloc((size_t)num_sel * sizeof(const char*));
    for (int32_t c = 0; c < num_sel; c++) {
        headers[c] = carquet_schema_column_name(schema, sel[c]);
    }
    print_dyn_table(headers, num_sel, (const char* const*)mat.cells, mat.num_rows);
    free(headers);

    matrix_free(&mat);
    if (has_filter) cli_filter_free(&fs);
    free(sel);
    carquet_reader_close(reader);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * cmd_export --format csv — write rows to stdout as CSV
 *
 * RFC 4180 quoting: fields containing comma, quote, CR, or LF are wrapped
 * in double quotes; embedded quotes are doubled. Header row first.
 * ══════════════════════════════════════════════════════════════════════════ */

static void emit_csv_field(const char* v) {
    if (!v) v = "";
    bool needs_quote = false;
    for (const char* p = v; *p; p++) {
        if (*p == ',' || *p == '"' || *p == '\n' || *p == '\r') {
            needs_quote = true;
            break;
        }
    }
    if (!needs_quote) {
        fputs(v, stdout);
        return;
    }
    fputc('"', stdout);
    for (const char* p = v; *p; p++) {
        if (*p == '"') fputc('"', stdout);
        fputc(*p, stdout);
    }
    fputc('"', stdout);
}

int cmd_export(const char* path, const row_select_opts_t* opts, export_format_t fmt) {
    if (fmt != CLI_EXPORT_CSV) {
        fprintf(stderr, "error: unsupported export format\n");
        return 1;
    }

    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = open_or_die(path, &err);
    if (!reader) return 1;

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    int64_t total = carquet_reader_num_rows(reader);

    int64_t offset = opts->offset < 0 ? 0 : opts->offset;
    if (offset > total) offset = total;
    int64_t limit = opts->limit < 0 ? (total - offset) : opts->limit;
    if (limit > total - offset) limit = total - offset;

    int32_t* sel = NULL;
    int32_t num_sel = resolve_columns(schema, num_cols, opts->columns, &sel);
    if (num_sel < 0) {
        fprintf(stderr, "error: no columns matched filter '%s'\n",
                opts->columns ? opts->columns : "");
        free(sel);
        carquet_reader_close(reader);
        return 1;
    }

    /* Header row (always emitted, even when limit==0). */
    for (int32_t c = 0; c < num_sel; c++) {
        if (c > 0) fputc(',', stdout);
        emit_csv_field(carquet_schema_column_name(schema, sel[c]));
    }
    fputc('\n', stdout);

    if (limit <= 0 || num_sel == 0) {
        free(sel);
        carquet_reader_close(reader);
        return 0;
    }

    cli_filter_storage_t fs;
    bool has_filter = false;
    if (opts->filter) {
        char ferr[256];
        if (cli_parse_filter(opts->filter, schema, num_cols, &fs,
                             ferr, sizeof(ferr)) != 0) {
            fprintf(stderr, "error: %s\n", ferr);
            free(sel);
            carquet_reader_close(reader);
            return 1;
        }
        has_filter = true;
    }

    str_matrix_t mat = {0};
    int rc;
    if (has_filter) {
        rc = read_rows_filtered(reader, schema, sel, num_sel, offset, limit,
                                &fs, "", &mat);
    } else {
        rc = read_rows(reader, schema, sel, num_sel, offset, limit, NULL, "", &mat);
    }
    if (rc != 0) {
        if (has_filter) cli_filter_free(&fs);
        matrix_free(&mat);
        free(sel);
        carquet_reader_close(reader);
        return 1;
    }

    for (int64_t r = 0; r < mat.num_rows; r++) {
        for (int32_t c = 0; c < num_sel; c++) {
            if (c > 0) fputc(',', stdout);
            emit_csv_field(mat.cells[r * num_sel + c]);
        }
        fputc('\n', stdout);
    }

    matrix_free(&mat);
    if (has_filter) cli_filter_free(&fs);
    free(sel);
    carquet_reader_close(reader);
    return 0;
}
