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
                              str_matrix_t* out);
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
                        cli_format_value(phys, stats.min_value, tl, lt,
                                         min_buf, sizeof(min_buf));
                        cli_format_value(phys, stats.max_value, tl, lt,
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
 * Table display helpers for head/tail/sample
 * ══════════════════════════════════════════════════════════════════════════ */

typedef struct {
    char** cells;       /* [row * num_cols + col] */
    int*   widths;      /* per column */
    int32_t num_cols;
    int64_t num_rows;
    int64_t capacity;
    const carquet_schema_t* schema;
} table_t;

static void table_init(table_t* t, const carquet_schema_t* schema, int32_t num_cols, int64_t cap) {
    t->schema = schema;
    t->num_cols = num_cols;
    t->num_rows = 0;
    t->capacity = cap;
    t->cells = calloc((size_t)(cap * num_cols), sizeof(char*));
    t->widths = calloc((size_t)num_cols, sizeof(int));

    /* Initialize widths from column names */
    for (int32_t c = 0; c < num_cols; c++) {
        const char* name = carquet_schema_column_name(schema, c);
        int len = (int)strlen(name);
        t->widths[c] = len < MAX_COL_WIDTH ? len : MAX_COL_WIDTH;
    }
}

static void table_add_cell(table_t* t, int64_t row, int32_t col, const char* value) {
    if (row >= t->capacity || col >= t->num_cols) return;
    t->cells[row * t->num_cols + col] = carquet_heap_strdup(value);
    int len = (int)strlen(value);
    if (len > MAX_COL_WIDTH) len = MAX_COL_WIDTH;
    if (len > t->widths[col]) t->widths[col] = len;
    if (row >= t->num_rows) t->num_rows = row + 1;
}

static void table_print(const table_t* t) {
    /* Header */
    printf("  ");
    for (int32_t c = 0; c < t->num_cols; c++) {
        if (c > 0) printf("  ");
        printf("%-*.*s", t->widths[c], t->widths[c],
               carquet_schema_column_name(t->schema, c));
    }
    printf("\n  ");
    for (int32_t c = 0; c < t->num_cols; c++) {
        if (c > 0) printf("  ");
        for (int w = 0; w < t->widths[c]; w++) putchar('-');
    }
    printf("\n");

    /* Rows */
    for (int64_t r = 0; r < t->num_rows; r++) {
        printf("  ");
        for (int32_t c = 0; c < t->num_cols; c++) {
            if (c > 0) printf("  ");
            const char* val = t->cells[r * t->num_cols + c];
            if (!val) val = "";
            printf("%-*.*s", t->widths[c], t->widths[c], val);
        }
        printf("\n");
    }
}

static void table_free(table_t* t) {
    if (t->cells) {
        for (int64_t i = 0; i < t->capacity * t->num_cols; i++)
            free(t->cells[i]);
        free(t->cells);
    }
    free(t->widths);
}

/* ══════════════════════════════════════════════════════════════════════════
 * cmd_head
 * ══════════════════════════════════════════════════════════════════════════ */

int cmd_head(const char* path, int64_t n, const char* filter) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_t* reader = open_or_die(path, &err);
    if (!reader) return 1;

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    int64_t total = carquet_reader_num_rows(reader);
    if (filter) {
        cli_filter_storage_t fs;
        char ferr[256];
        if (cli_parse_filter(filter, schema, num_cols, &fs, ferr, sizeof(ferr)) != 0) {
            fprintf(stderr, "error: %s\n", ferr);
            carquet_reader_close(reader);
            return 1;
        }
        int32_t* sel = malloc((size_t)num_cols * sizeof(int32_t));
        for (int32_t c = 0; c < num_cols; c++) sel[c] = c;
        str_matrix_t mat = {0};
        int rc = read_rows_filtered(reader, schema, sel, num_cols, 0, n,
                                    &fs, &mat);
        if (rc == 0) {
            const char** headers = malloc((size_t)num_cols * sizeof(const char*));
            for (int32_t c = 0; c < num_cols; c++)
                headers[c] = carquet_schema_column_name(schema, c);
            print_dyn_table(headers, num_cols,
                            (const char* const*)mat.cells, mat.num_rows);
            free(headers);
        }
        matrix_free(&mat);
        free(sel);
        cli_filter_free(&fs);
        carquet_reader_close(reader);
        return rc == 0 ? 0 : 1;
    }
    if (n > total) n = total;
    if (n <= 0 || num_cols <= 0) {
        carquet_reader_close(reader);
        return 0;
    }

    table_t tbl;
    table_init(&tbl, schema, num_cols, n);

    /* Read n rows from first row group(s) */
    for (int32_t c = 0; c < num_cols; c++) {
        carquet_physical_type_t phys = carquet_schema_column_type(schema, c);
        const carquet_schema_node_t* node = carquet_schema_get_element(schema,
            schema->leaf_indices[c]);
        const carquet_logical_type_t* lt = carquet_schema_node_logical_type(node);
        int32_t tl = carquet_schema_node_type_length(node);
        bool nullable = carquet_schema_node_repetition(node) != CARQUET_REPETITION_REQUIRED;
        int16_t max_def = carquet_schema_node_max_def_level(node);

        int64_t rows_read = 0;
        for (int32_t rg = 0; rg < carquet_reader_num_row_groups(reader) && rows_read < n; rg++) {
            carquet_column_reader_t* col = carquet_reader_get_column(reader, rg, c, &err);
            if (!col) continue;

            int64_t want = n - rows_read;

            /* Allocate buffer based on type */
            int32_t elem_size = carquet_physical_type_size(phys);
            void* buf;
            int16_t* def = NULL;
            if (phys == CARQUET_PHYSICAL_BYTE_ARRAY) {
                buf = calloc((size_t)want, sizeof(carquet_byte_array_t));
            } else if (phys == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY) {
                buf = calloc((size_t)want, (size_t)tl);
            } else {
                buf = calloc((size_t)want, (size_t)elem_size);
            }
            if (nullable)
                def = calloc((size_t)want, sizeof(int16_t));

            int64_t got = carquet_column_read_batch(col, buf, want, def, NULL);

            /* read_batch packs non-null values densely (no slot for nulls),
             * so buffer addressing advances only on present rows. */
            int64_t dense = 0;
            for (int64_t i = 0; i < got && rows_read + i < n; i++) {
                char vbuf[MAX_VALUE_BUF];
                if (nullable && def && def[i] < max_def) {
                    table_add_cell(&tbl, rows_read + i, c, "null");
                } else {
                    const void* vp = NULL;
                    if (phys == CARQUET_PHYSICAL_BYTE_ARRAY)
                        vp = &((carquet_byte_array_t*)buf)[dense];
                    else if (phys == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY)
                        vp = (uint8_t*)buf + dense * tl;
                    else
                        vp = (uint8_t*)buf + dense * elem_size;
                    dense++;

                    cli_format_value(phys, vp, tl, lt, vbuf, sizeof(vbuf));
                    table_add_cell(&tbl, rows_read + i, c, vbuf);
                }
            }

            rows_read += got;
            free(buf);
            free(def);
            carquet_column_reader_free(col);
        }
    }

    table_print(&tbl);
    table_free(&tbl);
    carquet_reader_close(reader);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * cmd_tail
 * ══════════════════════════════════════════════════════════════════════════ */

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

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    int32_t num_rgs = carquet_reader_num_row_groups(reader);
    int64_t total = carquet_reader_num_rows(reader);
    if (n > total) n = total;
    if (n <= 0 || num_cols <= 0) {
        carquet_reader_close(reader);
        return 0;
    }

    /* Figure out where to start reading:
     * skip_rows = total - n
     * Find the row group containing the start offset */
    int64_t skip_rows = total - n;

    table_t tbl;
    table_init(&tbl, schema, num_cols, n);

    for (int32_t c = 0; c < num_cols; c++) {
        carquet_physical_type_t phys = carquet_schema_column_type(schema, c);
        const carquet_schema_node_t* node = carquet_schema_get_element(schema,
            schema->leaf_indices[c]);
        const carquet_logical_type_t* lt = carquet_schema_node_logical_type(node);
        int32_t tl = carquet_schema_node_type_length(node);
        bool nullable = carquet_schema_node_repetition(node) != CARQUET_REPETITION_REQUIRED;
        int16_t max_def = carquet_schema_node_max_def_level(node);

        int64_t rows_seen = 0;
        int64_t rows_output = 0;

        for (int32_t rg = 0; rg < num_rgs && rows_output < n; rg++) {
            carquet_row_group_metadata_t rgm;
            (void)carquet_reader_row_group_metadata(reader, rg, &rgm);

            /* Skip entire row groups before the start */
            if (rows_seen + rgm.num_rows <= skip_rows) {
                rows_seen += rgm.num_rows;
                continue;
            }

            carquet_column_reader_t* col = carquet_reader_get_column(reader, rg, c, &err);
            if (!col) continue;

            /* Skip rows within this row group */
            int64_t skip_in_rg = skip_rows - rows_seen;
            if (skip_in_rg < 0) skip_in_rg = 0;
            if (skip_in_rg > 0)
                carquet_column_skip(col, skip_in_rg);

            int64_t want = n - rows_output;
            int32_t elem_size = carquet_physical_type_size(phys);
            void* buf;
            int16_t* def = NULL;
            if (phys == CARQUET_PHYSICAL_BYTE_ARRAY) {
                buf = calloc((size_t)want, sizeof(carquet_byte_array_t));
            } else if (phys == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY) {
                buf = calloc((size_t)want, (size_t)tl);
            } else {
                buf = calloc((size_t)want, (size_t)elem_size);
            }
            if (nullable)
                def = calloc((size_t)want, sizeof(int16_t));

            int64_t got = carquet_column_read_batch(col, buf, want, def, NULL);

            /* read_batch packs non-null values densely (no slot for nulls),
             * so buffer addressing advances only on present rows. */
            int64_t dense = 0;
            for (int64_t i = 0; i < got && rows_output < n; i++) {
                char vbuf[MAX_VALUE_BUF];
                if (nullable && def && def[i] < max_def) {
                    table_add_cell(&tbl, rows_output, c, "null");
                } else {
                    const void* vp = NULL;
                    if (phys == CARQUET_PHYSICAL_BYTE_ARRAY)
                        vp = &((carquet_byte_array_t*)buf)[dense];
                    else if (phys == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY)
                        vp = (uint8_t*)buf + dense * tl;
                    else
                        vp = (uint8_t*)buf + dense * elem_size;
                    dense++;

                    cli_format_value(phys, vp, tl, lt, vbuf, sizeof(vbuf));
                    table_add_cell(&tbl, rows_output, c, vbuf);
                }
                rows_output++;
            }

            rows_seen += rgm.num_rows;
            free(buf);
            free(def);
            carquet_column_reader_free(col);
        }
    }

    table_print(&tbl);
    table_free(&tbl);
    carquet_reader_close(reader);
    return 0;
}

/* ══════════════════════════════════════════════════════════════════════════
 * cmd_sample
 * ══════════════════════════════════════════════════════════════════════════ */

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

    const carquet_schema_t* schema = carquet_reader_schema(reader);
    int32_t num_cols = carquet_reader_num_columns(reader);
    int64_t total = carquet_reader_num_rows(reader);
    if (n > total) n = total;
    if (n <= 0 || num_cols <= 0) {
        carquet_reader_close(reader);
        return 0;
    }

    /* Generate n sorted random row indices using reservoir sampling.
     * For simplicity, just pick n random indices. */
    srand((unsigned)time(NULL));
    int64_t* indices = calloc((size_t)n, sizeof(int64_t));
    for (int64_t i = 0; i < n; i++) {
        indices[i] = ((int64_t)rand() * rand()) % total;
    }
    qsort(indices, (size_t)n, sizeof(int64_t), compare_int64);

    /* Remove duplicates */
    int64_t unique = 1;
    for (int64_t i = 1; i < n; i++) {
        if (indices[i] != indices[unique - 1])
            indices[unique++] = indices[i];
    }
    n = unique;

    /* Read sampled rows. For each column, we use head-style reading
     * with skip to jump to each sampled row. */
    table_t tbl;
    table_init(&tbl, schema, num_cols, n);
    void* valbuf = NULL;        /* one-value read buffer, sized per column */
    size_t valbuf_cap = 0;

    for (int32_t c = 0; c < num_cols; c++) {
        carquet_physical_type_t phys = carquet_schema_column_type(schema, c);
        const carquet_schema_node_t* node = carquet_schema_get_element(schema,
            schema->leaf_indices[c]);
        const carquet_logical_type_t* lt = carquet_schema_node_logical_type(node);
        int32_t tl = carquet_schema_node_type_length(node);
        bool nullable = carquet_schema_node_repetition(node) != CARQUET_REPETITION_REQUIRED;
        int16_t max_def = carquet_schema_node_max_def_level(node);
        int32_t num_rgs = carquet_reader_num_row_groups(reader);

        int64_t sample_idx = 0;     /* Index into sorted indices[] */
        int64_t rg_row_start = 0;   /* Absolute row offset of current row group */

        for (int32_t rg = 0; rg < num_rgs && sample_idx < n; rg++) {
            carquet_row_group_metadata_t rgm;
            (void)carquet_reader_row_group_metadata(reader, rg, &rgm);
            int64_t rg_row_end = rg_row_start + rgm.num_rows;

            /* Skip row groups with no sampled rows */
            if (sample_idx < n && indices[sample_idx] >= rg_row_end) {
                rg_row_start = rg_row_end;
                continue;
            }

            carquet_column_reader_t* col = carquet_reader_get_column(reader, rg, c, &err);
            if (!col) { rg_row_start = rg_row_end; continue; }

            int64_t pos_in_rg = 0; /* Current position within the row group */

            while (sample_idx < n && indices[sample_idx] < rg_row_end) {
                int64_t target_in_rg = indices[sample_idx] - rg_row_start;
                int64_t skip = target_in_rg - pos_in_rg;
                if (skip > 0) {
                    carquet_column_skip(col, skip);
                    pos_in_rg += skip;
                }

                /* Read one value. A FIXED_LEN_BYTE_ARRAY value is
                 * `type_length` bytes (file-controlled, up to 16 MiB), so
                 * size the buffer from the column, not a fixed union. */
                size_t val_cap = sizeof(carquet_byte_array_t) > 16
                    ? sizeof(carquet_byte_array_t) : 16;
                if (phys == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY && tl > 0 &&
                    (size_t)tl > val_cap) {
                    val_cap = (size_t)tl;
                }
                if (val_cap > valbuf_cap) {
                    void* nb = realloc(valbuf, val_cap);
                    if (!nb) { free(valbuf); valbuf = NULL; valbuf_cap = 0; break; }
                    valbuf = nb; valbuf_cap = val_cap;
                }
                int16_t def_level = 0;

                int64_t got = carquet_column_read_batch(col, valbuf, 1,
                    nullable ? &def_level : NULL, NULL);
                pos_in_rg++;

                char vbuf[MAX_VALUE_BUF];
                if (got <= 0) {
                    table_add_cell(&tbl, sample_idx, c, "?");
                } else if (nullable && def_level < max_def) {
                    table_add_cell(&tbl, sample_idx, c, "null");
                } else {
                    cli_format_value(phys, valbuf, tl, lt, vbuf, sizeof(vbuf));
                    table_add_cell(&tbl, sample_idx, c, vbuf);
                }

                sample_idx++;
            }

            carquet_column_reader_free(col);
            rg_row_start = rg_row_end;
        }
    }

    table_print(&tbl);
    table_free(&tbl);
    free(valbuf);
    free(indices);
    carquet_reader_close(reader);
    return 0;
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
}

/* Read the requested column at `col_index`, skipping `offset` rows and
 * filling at most `limit` formatted strings into `matrix` at column slot
 * `dst_col`. Returns the number of rows actually filled. */
static int64_t read_column_strings(carquet_reader_t* reader,
                                    const carquet_schema_t* schema,
                                    int32_t col_index,
                                    int64_t offset, int64_t limit,
                                    str_matrix_t* matrix, int32_t dst_col) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_physical_type_t phys = carquet_schema_column_type(schema, col_index);
    const carquet_schema_node_t* node = carquet_schema_get_element(schema,
        schema->leaf_indices[col_index]);
    const carquet_logical_type_t* lt = carquet_schema_node_logical_type(node);
    int32_t tl = carquet_schema_node_type_length(node);
    bool nullable = carquet_schema_node_repetition(node) != CARQUET_REPETITION_REQUIRED;
    int16_t max_def = carquet_schema_node_max_def_level(node);

    int32_t num_rgs = carquet_reader_num_row_groups(reader);
    int64_t rows_seen = 0;
    int64_t rows_output = 0;

    for (int32_t rg = 0; rg < num_rgs && rows_output < limit; rg++) {
        carquet_row_group_metadata_t rgm;
        (void)carquet_reader_row_group_metadata(reader, rg, &rgm);

        if (rows_seen + rgm.num_rows <= offset) {
            rows_seen += rgm.num_rows;
            continue;
        }

        carquet_column_reader_t* col = carquet_reader_get_column(reader, rg,
            col_index, &err);
        if (!col) { rows_seen += rgm.num_rows; continue; }

        int64_t skip_in_rg = offset - rows_seen;
        if (skip_in_rg < 0) skip_in_rg = 0;
        if (skip_in_rg > 0) carquet_column_skip(col, skip_in_rg);

        int64_t want = limit - rows_output;
        int64_t rg_remaining = rgm.num_rows - skip_in_rg;
        if (want > rg_remaining) want = rg_remaining;

        int32_t elem_size = carquet_physical_type_size(phys);
        void* buf;
        int16_t* def = NULL;
        if (phys == CARQUET_PHYSICAL_BYTE_ARRAY)
            buf = calloc((size_t)want, sizeof(carquet_byte_array_t));
        else if (phys == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY)
            buf = calloc((size_t)want, (size_t)tl);
        else
            buf = calloc((size_t)want, (size_t)elem_size);
        if (nullable) def = calloc((size_t)want, sizeof(int16_t));

        int64_t got = carquet_column_read_batch(col, buf, want, def, NULL);

        /* read_batch packs non-null values densely (no slot for nulls),
         * so buffer addressing advances only on present rows. */
        int64_t dense = 0;
        for (int64_t i = 0; i < got && rows_output < limit; i++) {
            char vbuf[MAX_VALUE_BUF];
            const char* cell;
            if (nullable && def && def[i] < max_def) {
                cell = "";
            } else {
                const void* vp;
                if (phys == CARQUET_PHYSICAL_BYTE_ARRAY)
                    vp = &((carquet_byte_array_t*)buf)[dense];
                else if (phys == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY)
                    vp = (uint8_t*)buf + dense * tl;
                else
                    vp = (uint8_t*)buf + dense * elem_size;
                dense++;
                cli_format_value(phys, vp, tl, lt, vbuf, sizeof(vbuf));
                cell = vbuf;
            }
            matrix->cells[rows_output * matrix->num_cols + dst_col] =
                carquet_heap_strdup(cell);
            rows_output++;
        }

        rows_seen += rgm.num_rows;
        free(buf);
        free(def);
        carquet_column_reader_free(col);
    }

    return rows_output;
}

static int read_rows(carquet_reader_t* reader,
                     const carquet_schema_t* schema,
                     const int32_t* col_indices, int32_t num_sel_cols,
                     int64_t offset, int64_t limit,
                     str_matrix_t* out) {
    out->num_cols = num_sel_cols;
    out->num_rows = limit;
    out->cells = calloc((size_t)(limit * num_sel_cols), sizeof(char*));
    if (!out->cells) return -1;

    int64_t produced = 0;
    for (int32_t i = 0; i < num_sel_cols; i++) {
        int64_t n = read_column_strings(reader, schema, col_indices[i],
                                         offset, limit, out, i);
        if (n > produced) produced = n;
    }
    out->num_rows = produced;
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

static int read_rows_filtered(carquet_reader_t* reader,
                              const carquet_schema_t* schema,
                              const int32_t* col_indices, int32_t num_sel_cols,
                              int64_t offset, int64_t limit,
                              const cli_filter_storage_t* filter,
                              str_matrix_t* out) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    out->num_cols = num_sel_cols;
    out->num_rows = 0;
    out->cells = NULL;
    if (limit <= 0 || num_sel_cols <= 0) return 0;

    out->cells = calloc((size_t)(limit * num_sel_cols), sizeof(char*));
    if (!out->cells) return -1;

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
    if (!phys || !lts || !tls || !max_defs) {
        free(phys); free(lts); free(tls); free(max_defs);
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
        max_defs[c] = carquet_schema_node_max_def_level(node);
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
                    cell = "";
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
    free(phys); free(lts); free(tls); free(max_defs);
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
                                &fs, &mat);
    } else {
        rc = read_rows(reader, schema, sel, num_sel, offset, limit, &mat);
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
                                &fs, &mat);
    } else {
        rc = read_rows(reader, schema, sel, num_sel, offset, limit, &mat);
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
