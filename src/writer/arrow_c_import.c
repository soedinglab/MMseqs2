/**
 * @file arrow_c_import.c
 * @brief Import Arrow C Data Interface structs into Carquet.
 *
 * Accepts standard `ArrowSchema` / `ArrowArray` structs (see carquet.h) and
 * either builds a Carquet schema (@ref carquet_arrow_import_schema) or writes a
 * struct array to a writer (@ref carquet_writer_write_arrow).
 *
 * Nesting: both entry points handle arbitrary depth. `carquet_writer_write_arrow`
 * runs a generic Dremel record-shredding pass that walks the Arrow (schema,
 * array) tree and the writer's leaf columns in lockstep, producing, for every
 * leaf, its repetition levels, definition levels, and dense (present-only)
 * values, then hands each to @ref carquet_writer_write_batch. This covers
 * struct, list, large-list and map nodes composed to any depth.
 *
 * Following Arrow "move" semantics, both entry points consume the structs they
 * are given: the `release` callback is invoked before returning, on success and
 * on failure alike.
 *
 * Uses only the public Carquet API (plus the standard C allocator for the
 * transient shredding buffers).
 */

#include <carquet/carquet.h>

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Present-bit accessor for an Arrow validity bitmap (LSB-first). */
#define ARROW_VALID(buf, i) (((buf)[(size_t)(i) >> 3] >> ((i) & 7)) & 1u)

/* Bound recursion depth so a hostile / cyclic-looking schema cannot blow the
 * stack. Real-world nesting is a handful of levels. */
#define ARROW_MAX_NEST_DEPTH 64

static void release_schema(struct ArrowSchema* s) {
    if (s && s->release) s->release(s);
}
static void release_array(struct ArrowArray* a) {
    if (a && a->release) a->release(a);
}

/* ============================================================================
 * Arrow format string classification
 * ============================================================================
 */

typedef enum {
    ANODE_PRIMITIVE,
    ANODE_STRUCT,   /* "+s" */
    ANODE_LIST,     /* "+l" (int32 offsets) or "+L" (int64 offsets) */
    ANODE_MAP,      /* "+m" */
    ANODE_UNSUPPORTED
} anode_kind_t;

/* "+w:<n>" (fixed-size list): the element count per slot, 0 if `fmt` is not a
 * well-formed fixed-size-list format. */
static int64_t fixed_list_size(const char* fmt) {
    if (!fmt || fmt[0] != '+' || fmt[1] != 'w' || fmt[2] != ':' || fmt[3] == '\0') return 0;
    char* end = NULL;
    long n = strtol(fmt + 3, &end, 10);
    if (n <= 0 || n > INT32_MAX || !end || *end != '\0') return 0;
    return (int64_t)n;
}

static anode_kind_t classify(const char* fmt) {
    if (!fmt || !fmt[0]) return ANODE_UNSUPPORTED;
    if (fmt[0] != '+') return ANODE_PRIMITIVE;
    switch (fmt[1]) {
    case 's': return ANODE_STRUCT;
    case 'l': case 'L': return ANODE_LIST;
    case 'm': return ANODE_MAP;
    case 'w': return fixed_list_size(fmt) > 0 ? ANODE_LIST : ANODE_UNSUPPORTED;
    default:  return ANODE_UNSUPPORTED;  /* +ud/+us unions, +r run-end, ... */
    }
}

/* ============================================================================
 * Arrow format string -> Carquet physical/logical type
 * ============================================================================
 */
/* Value conversion applied while shredding, when the Arrow in-memory
 * representation is not the Parquet physical one. */
typedef enum {
    CONV_NONE = 0,
    CONV_I8, CONV_U8, CONV_I16, CONV_U16,  /* Arrow narrow ints -> INT32 */
    CONV_DATE64,    /* int64 milliseconds -> INT32 days (DATE) */
    CONV_TIME32_S,  /* int32 seconds -> INT32 milliseconds (TIME(MILLIS)) */
    CONV_TS_S,      /* int64 seconds -> INT64 milliseconds (TIMESTAMP(MILLIS)) */
    CONV_DECIMAL    /* little-endian two's complement -> big-endian FLBA (DECIMAL) */
} arrow_conv_t;

/* Byte length of the FIXED_LEN_BYTE_ARRAY used to store a DECIMAL of the given
 * precision. Matches parquet-cpp's DecimalSize() up to 38 digits so that files
 * written from an Arrow decimal128 are byte-identical to pyarrow's; wider
 * precisions (decimal256) take the smallest length that holds them. */
static int32_t decimal_flba_length(int32_t precision) {
    static const int8_t table[39] = {
        0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 6, 6, 6, 7, 7, 8, 8, 9, 9, 9,
        10, 10, 11, 11, 11, 12, 12, 13, 13, 14, 14, 14, 15, 15, 16, 16, 16 };
    if (precision <= 0) return 0;
    if (precision <= 38) return table[precision];
    if (precision > 76) return 0;
    /* smallest n with floor((8n - 1) * log10(2)) >= precision */
    for (int32_t n = 17; n <= 32; n++) {
        if ((int32_t)((8.0 * n - 1.0) * 0.30102999566398120) >= precision) return n;
    }
    return 0;
}

static carquet_status_t parse_arrow_format(
    const char* fmt,
    carquet_physical_type_t* pt,
    carquet_logical_type_t* lt,
    bool* has_lt,
    int32_t* type_length,
    arrow_conv_t* conv,
    size_t* src_stride) {

    *has_lt = false;
    *type_length = 0;
    *conv = CONV_NONE;
    *src_stride = 0;
    memset(lt, 0, sizeof(*lt));

    if (!fmt || !fmt[0]) return CARQUET_ERROR_INVALID_ARGUMENT;

    /* Single-character primitives */
    if (fmt[1] == '\0') {
        switch (fmt[0]) {
        case 'b': *pt = CARQUET_PHYSICAL_BOOLEAN; return CARQUET_OK;
        case 'c': *pt = CARQUET_PHYSICAL_INT32; *has_lt = true;
                  lt->id = CARQUET_LOGICAL_INTEGER;
                  lt->params.integer.bit_width = 8; lt->params.integer.is_signed = true;
                  *conv = CONV_I8; *src_stride = 1;
                  return CARQUET_OK;
        case 'C': *pt = CARQUET_PHYSICAL_INT32; *has_lt = true;
                  lt->id = CARQUET_LOGICAL_INTEGER;
                  lt->params.integer.bit_width = 8; lt->params.integer.is_signed = false;
                  *conv = CONV_U8; *src_stride = 1;
                  return CARQUET_OK;
        case 's': *pt = CARQUET_PHYSICAL_INT32; *has_lt = true;
                  lt->id = CARQUET_LOGICAL_INTEGER;
                  lt->params.integer.bit_width = 16; lt->params.integer.is_signed = true;
                  *conv = CONV_I16; *src_stride = 2;
                  return CARQUET_OK;
        case 'S': *pt = CARQUET_PHYSICAL_INT32; *has_lt = true;
                  lt->id = CARQUET_LOGICAL_INTEGER;
                  lt->params.integer.bit_width = 16; lt->params.integer.is_signed = false;
                  *conv = CONV_U16; *src_stride = 2;
                  return CARQUET_OK;
        case 'i': *pt = CARQUET_PHYSICAL_INT32; return CARQUET_OK;
        case 'I': *pt = CARQUET_PHYSICAL_INT32; *has_lt = true;
                  lt->id = CARQUET_LOGICAL_INTEGER;
                  lt->params.integer.bit_width = 32; lt->params.integer.is_signed = false;
                  return CARQUET_OK;
        case 'l': *pt = CARQUET_PHYSICAL_INT64; return CARQUET_OK;
        case 'L': *pt = CARQUET_PHYSICAL_INT64; *has_lt = true;
                  lt->id = CARQUET_LOGICAL_INTEGER;
                  lt->params.integer.bit_width = 64; lt->params.integer.is_signed = false;
                  return CARQUET_OK;
        case 'f': *pt = CARQUET_PHYSICAL_FLOAT; return CARQUET_OK;
        case 'g': *pt = CARQUET_PHYSICAL_DOUBLE; return CARQUET_OK;
        case 'e': *pt = CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY; *type_length = 2;
                  *has_lt = true; lt->id = CARQUET_LOGICAL_FLOAT16; return CARQUET_OK;
        case 'u': case 'U': *pt = CARQUET_PHYSICAL_BYTE_ARRAY; *has_lt = true;
                  lt->id = CARQUET_LOGICAL_STRING; return CARQUET_OK;   /* utf8 / large utf8 */
        case 'z': case 'Z': *pt = CARQUET_PHYSICAL_BYTE_ARRAY; return CARQUET_OK; /* binary / large binary */
        default: return CARQUET_ERROR_NOT_IMPLEMENTED;
        }
    }

    /* fixed_size_binary "w:<n>" */
    if (fmt[0] == 'w' && fmt[1] == ':') {
        long n = strtol(fmt + 2, NULL, 10);
        if (n <= 0 || n > (16 * 1024 * 1024)) return CARQUET_ERROR_INVALID_ARGUMENT;
        *pt = CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY;
        *type_length = (int32_t)n;
        return CARQUET_OK;
    }

    /* decimal "d:<precision>,<scale>[,<bitwidth>]" -> FIXED_LEN_BYTE_ARRAY DECIMAL.
     * Arrow stores little-endian two's complement of 32/64/128/256 bits; Parquet
     * wants the big-endian minimal-length representation (see CONV_DECIMAL). */
    if (fmt[0] == 'd' && fmt[1] == ':') {
        char* end = NULL;
        long precision = strtol(fmt + 2, &end, 10);
        if (!end || *end != ',' || precision <= 0 || precision > 76) {
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        long scale = strtol(end + 1, &end, 10);
        if (!end || scale < 0 || scale > precision) return CARQUET_ERROR_INVALID_ARGUMENT;
        long bits = 128;
        if (*end == ',') {
            bits = strtol(end + 1, &end, 10);
            if (!end || *end != '\0') return CARQUET_ERROR_INVALID_ARGUMENT;
            if (bits != 32 && bits != 64 && bits != 128 && bits != 256) {
                return CARQUET_ERROR_INVALID_ARGUMENT;
            }
        } else if (*end != '\0') {
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        int32_t len = decimal_flba_length((int32_t)precision);
        if (len <= 0 || (size_t)len > (size_t)bits / 8) return CARQUET_ERROR_INVALID_ARGUMENT;
        *pt = CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY;
        *type_length = len;
        *has_lt = true;
        lt->id = CARQUET_LOGICAL_DECIMAL;
        lt->params.decimal.precision = (int32_t)precision;
        lt->params.decimal.scale = (int32_t)scale;
        *conv = CONV_DECIMAL;
        *src_stride = (size_t)bits / 8;
        return CARQUET_OK;
    }

    /* Temporal: "td*", "tt*", "ts*" */
    if (fmt[0] == 't') {
        if (fmt[1] == 'd') {              /* date */
            if (fmt[2] == 'D' && fmt[3] == '\0') {
                *pt = CARQUET_PHYSICAL_INT32; *has_lt = true;
                lt->id = CARQUET_LOGICAL_DATE; return CARQUET_OK;
            }
            if (fmt[2] == 'm' && fmt[3] == '\0') {   /* date64 (milliseconds) */
                *pt = CARQUET_PHYSICAL_INT32; *has_lt = true;
                lt->id = CARQUET_LOGICAL_DATE;
                *conv = CONV_DATE64; *src_stride = 8;
                return CARQUET_OK;
            }
            return CARQUET_ERROR_NOT_IMPLEMENTED;
        }
        if (fmt[1] == 't') {              /* time */
            *has_lt = true; lt->id = CARQUET_LOGICAL_TIME;
            lt->params.time.is_adjusted_to_utc = false;
            if (fmt[2] == 's' && fmt[3] == '\0') {   /* time32[s] -> TIME(MILLIS) */
                *pt = CARQUET_PHYSICAL_INT32;
                lt->params.time.unit = CARQUET_TIME_UNIT_MILLIS;
                *conv = CONV_TIME32_S; *src_stride = 4;
                return CARQUET_OK;
            }
            if (fmt[2] == 'm' && fmt[3] == '\0') {
                *pt = CARQUET_PHYSICAL_INT32;
                lt->params.time.unit = CARQUET_TIME_UNIT_MILLIS; return CARQUET_OK;
            }
            if (fmt[2] == 'u' && fmt[3] == '\0') {
                *pt = CARQUET_PHYSICAL_INT64;
                lt->params.time.unit = CARQUET_TIME_UNIT_MICROS; return CARQUET_OK;
            }
            if (fmt[2] == 'n' && fmt[3] == '\0') {
                *pt = CARQUET_PHYSICAL_INT64;
                lt->params.time.unit = CARQUET_TIME_UNIT_NANOS; return CARQUET_OK;
            }
            return CARQUET_ERROR_NOT_IMPLEMENTED;
        }
        if (fmt[1] == 's') {              /* timestamp "ts<unit>:<tz>" */
            if (fmt[2] == '\0' || fmt[3] != ':') return CARQUET_ERROR_INVALID_ARGUMENT;
            *pt = CARQUET_PHYSICAL_INT64; *has_lt = true;
            lt->id = CARQUET_LOGICAL_TIMESTAMP;
            if (fmt[2] == 'm') lt->params.timestamp.unit = CARQUET_TIME_UNIT_MILLIS;
            else if (fmt[2] == 'u') lt->params.timestamp.unit = CARQUET_TIME_UNIT_MICROS;
            else if (fmt[2] == 'n') lt->params.timestamp.unit = CARQUET_TIME_UNIT_NANOS;
            else if (fmt[2] == 's') {     /* timestamp[s] -> TIMESTAMP(MILLIS) */
                lt->params.timestamp.unit = CARQUET_TIME_UNIT_MILLIS;
                *conv = CONV_TS_S; *src_stride = 8;
            }
            else return CARQUET_ERROR_NOT_IMPLEMENTED;
            lt->params.timestamp.is_adjusted_to_utc = (fmt[4] != '\0');
            return CARQUET_OK;
        }
    }

    return CARQUET_ERROR_NOT_IMPLEMENTED;
}

/* In-memory byte width of a fixed-size physical type. 0 for BYTE_ARRAY. */
static size_t fixed_stride(carquet_physical_type_t pt, int32_t type_length) {
    switch (pt) {
    case CARQUET_PHYSICAL_BOOLEAN:  return 1;
    case CARQUET_PHYSICAL_INT32:    return 4;
    case CARQUET_PHYSICAL_INT64:    return 8;
    case CARQUET_PHYSICAL_INT96:    return 12;
    case CARQUET_PHYSICAL_FLOAT:    return 4;
    case CARQUET_PHYSICAL_DOUBLE:   return 8;
    case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY: return (size_t)type_length;
    default: return 0;
    }
}

/* Number of leaf columns under an Arrow schema node. */
static int64_t arrow_leaf_count(const struct ArrowSchema* s) {
    if (!s || !s->format) return 0;
    switch (classify(s->format)) {
    case ANODE_PRIMITIVE: return 1;
    case ANODE_STRUCT: {
        int64_t n = 0;
        for (int64_t i = 0; i < s->n_children; i++) n += arrow_leaf_count(s->children[i]);
        return n;
    }
    case ANODE_LIST:
    case ANODE_MAP:
        return (s->n_children == 1) ? arrow_leaf_count(s->children[0]) : 0;
    default:
        return 0;
    }
}

/* ============================================================================
 * Schema import (arbitrary depth)
 * ============================================================================
 */

/* Recursively add one Arrow field (and its subtree) under `parent_elem`. */
static carquet_status_t import_field(
    carquet_schema_t* cs,
    const struct ArrowSchema* s,
    int32_t parent_elem,
    const char* name_override,
    bool force_required,
    int depth,
    carquet_error_t* error) {

    if (depth > ARROW_MAX_NEST_DEPTH) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "Arrow import: schema nesting too deep");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    if (!s || !s->format) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "Arrow import: null child schema");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    const char* name = name_override ? name_override : (s->name ? s->name : "");
    bool nullable = (s->flags & ARROW_FLAG_NULLABLE) != 0;
    carquet_field_repetition_t rep = (nullable && !force_required)
        ? CARQUET_REPETITION_OPTIONAL : CARQUET_REPETITION_REQUIRED;

    switch (classify(s->format)) {
    case ANODE_STRUCT: {
        int32_t g = carquet_schema_add_group(cs, name, rep, parent_elem);
        if (g < 0) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INTERNAL, "add_group failed");
            return CARQUET_ERROR_INTERNAL;
        }
        for (int64_t i = 0; i < s->n_children; i++) {
            carquet_status_t rc = import_field(cs, s->children[i], g, NULL, false,
                                               depth + 1, error);
            if (rc != CARQUET_OK) return rc;
        }
        return CARQUET_OK;
    }
    case ANODE_LIST: {
        if (s->n_children != 1) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "Arrow import: list \"%s\" must have exactly one child", name);
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        int32_t inner = carquet_schema_add_list_group(cs, name, rep, parent_elem);
        if (inner < 0) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INTERNAL, "add_list_group failed");
            return CARQUET_ERROR_INTERNAL;
        }
        return import_field(cs, s->children[0], inner, "element", false,
                            depth + 1, error);
    }
    case ANODE_MAP: {
        if (s->n_children != 1 || !s->children[0] ||
            classify(s->children[0]->format) != ANODE_STRUCT ||
            s->children[0]->n_children != 2) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "Arrow import: map \"%s\" must have a 2-field struct child", name);
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        int32_t inner = carquet_schema_add_map_group(cs, name, rep, parent_elem);
        if (inner < 0) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INTERNAL, "add_map_group failed");
            return CARQUET_ERROR_INTERNAL;
        }
        const struct ArrowSchema* entries = s->children[0];
        carquet_status_t rc = import_field(cs, entries->children[0], inner, "key",
                                           /*force_required=*/true, depth + 1, error);
        if (rc != CARQUET_OK) return rc;
        return import_field(cs, entries->children[1], inner, "value", false,
                            depth + 1, error);
    }
    case ANODE_PRIMITIVE: {
        /* Dictionary-encoded array (e.g. a pandas categorical): `s->format` is
         * only the integer index type; the column's real type is the dictionary
         * value type in `s->dictionary`. Resolve to the value type so the
         * indices are materialized to values on write instead of being written
         * verbatim as an integer column. Only primitive dictionary values are
         * supported (nested dictionaries are rare and would be silent data loss
         * otherwise). */
        const char* value_format = s->format;
        if (s->dictionary) {
            if (!s->dictionary->format ||
                classify(s->dictionary->format) != ANODE_PRIMITIVE) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
                    "Arrow import: only primitive dictionary values are supported "
                    "for field \"%s\"", name);
                return CARQUET_ERROR_NOT_IMPLEMENTED;
            }
            value_format = s->dictionary->format;
        }
        carquet_physical_type_t pt;
        carquet_logical_type_t lt;
        bool has_lt;
        int32_t tlen;
        arrow_conv_t conv;
        size_t src_stride;
        carquet_status_t rc = parse_arrow_format(value_format, &pt, &lt, &has_lt, &tlen,
                                                 &conv, &src_stride);
        if (rc != CARQUET_OK) {
            CARQUET_SET_ERROR(error, rc,
                "Arrow import: unsupported format \"%s\" for field \"%s\"",
                value_format, name);
            return rc;
        }
        rc = carquet_schema_add_column(cs, name, pt, has_lt ? &lt : NULL, rep, tlen,
                                       parent_elem);
        if (rc != CARQUET_OK) {
            CARQUET_SET_ERROR(error, rc, "Arrow import: add_column failed for \"%s\"", name);
            return rc;
        }
        return CARQUET_OK;
    }
    default:
        CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
            "Arrow import: unsupported nested type \"%s\" for field \"%s\"",
            s->format, name);
        return CARQUET_ERROR_NOT_IMPLEMENTED;
    }
}

carquet_status_t carquet_arrow_import_schema(
    struct ArrowSchema* schema,
    carquet_schema_t** out,
    carquet_error_t* error) {

    if (!schema || !out) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "NULL schema or out");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    *out = NULL;

    carquet_status_t rc = CARQUET_OK;
    carquet_schema_t* cs = NULL;

    if (!schema->format || strcmp(schema->format, "+s") != 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "Arrow import: top-level schema must be a struct (\"+s\")");
        rc = CARQUET_ERROR_INVALID_ARGUMENT;
        goto done;
    }

    cs = carquet_schema_create(error);
    if (!cs) { rc = CARQUET_ERROR_OUT_OF_MEMORY; goto done; }

    for (int64_t i = 0; i < schema->n_children; i++) {
        rc = import_field(cs, schema->children[i], /*parent=*/0, NULL, false, 0, error);
        if (rc != CARQUET_OK) goto done;
    }

done:
    if (rc == CARQUET_OK) {
        *out = cs;
    } else if (cs) {
        carquet_schema_free(cs);
    }
    release_schema(schema);   /* consume, per Arrow move semantics */
    return rc;
}

/* ============================================================================
 * Generic Dremel shredder (write path)
 * ============================================================================
 */

/* Per-leaf accumulator: repetition/definition levels plus dense present-only
 * values, grown as the record recursion visits the leaf. */
typedef struct {
    carquet_physical_type_t pt;
    bool is_bool;
    bool is_bytearray;
    bool off64;             /* byte-array with 64-bit offsets (large utf8/binary) */
    size_t stride;          /* fixed-width byte stride written (0 for byte array) */
    arrow_conv_t conv;      /* value conversion from the Arrow representation */
    size_t src_stride;      /* Arrow in-memory stride when conv != CONV_NONE */
    int16_t max_def;        /* definition level of a present value */
    int16_t max_rep;        /* repetition level of the leaf */

    int16_t* def;
    int16_t* rep;
    int64_t  nlevels;
    int64_t  cap_levels;

    uint8_t* fixed;         /* fixed-width / boolean dense buffer */
    int64_t  fixed_bytes;
    int64_t  cap_fixed;

    carquet_byte_array_t* ba;  /* byte-array dense buffer */
    int64_t  nba;
    int64_t  cap_ba;

    bool oom;
} leaf_acc_t;

static bool grow_levels(leaf_acc_t* a, int64_t need) {
    if (a->nlevels + need <= a->cap_levels) return true;
    int64_t cap = a->cap_levels ? a->cap_levels : 64;
    while (cap < a->nlevels + need) cap *= 2;
    int16_t* nd = (int16_t*)realloc(a->def, (size_t)cap * sizeof(int16_t));
    int16_t* nr = (int16_t*)realloc(a->rep, (size_t)cap * sizeof(int16_t));
    if (nd) a->def = nd;
    if (nr) a->rep = nr;
    if (!nd || !nr) { a->oom = true; return false; }
    a->cap_levels = cap;
    return true;
}

static void append_level(leaf_acc_t* a, int16_t rep, int16_t def) {
    if (!grow_levels(a, 1)) return;
    a->rep[a->nlevels] = rep;
    a->def[a->nlevels] = def;
    a->nlevels++;
}

static bool grow_fixed(leaf_acc_t* a, int64_t need) {
    if (a->fixed_bytes + need <= a->cap_fixed) return true;
    int64_t cap = a->cap_fixed ? a->cap_fixed : 256;
    while (cap < a->fixed_bytes + need) cap *= 2;
    uint8_t* n = (uint8_t*)realloc(a->fixed, (size_t)cap);
    if (!n) { a->oom = true; return false; }
    a->fixed = n; a->cap_fixed = cap;
    return true;
}

static bool grow_ba(leaf_acc_t* a, int64_t need) {
    if (a->nba + need <= a->cap_ba) return true;
    int64_t cap = a->cap_ba ? a->cap_ba : 64;
    while (cap < a->nba + need) cap *= 2;
    carquet_byte_array_t* n =
        (carquet_byte_array_t*)realloc(a->ba, (size_t)cap * sizeof(carquet_byte_array_t));
    if (!n) { a->oom = true; return false; }
    a->ba = n; a->cap_ba = cap;
    return true;
}

/* Context threaded through the record recursion. */
typedef struct {
    leaf_acc_t* accs;
    int32_t num_leaves;
    carquet_error_t* error;
    carquet_status_t rc;
} shred_ctx_t;

static bool arr_is_null(const struct ArrowArray* a, int64_t i) {
    const uint8_t* v = (a->n_buffers >= 1) ? (const uint8_t*)a->buffers[0] : NULL;
    if (!v) return false;
    return !ARROW_VALID(v, i);
}

/* Append the present value at source index i of the leaf array `a`. */
static void emit_value(shred_ctx_t* ctx, leaf_acc_t* acc,
                       const struct ArrowArray* a, int64_t i) {
    if (acc->is_bool) {
        const uint8_t* bits = (a->n_buffers >= 2) ? (const uint8_t*)a->buffers[1] : NULL;
        if (!bits) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }
        if (!grow_fixed(acc, 1)) return;
        acc->fixed[acc->fixed_bytes] = (uint8_t)ARROW_VALID(bits, i);
        acc->fixed_bytes += 1;
        return;
    }
    if (acc->is_bytearray) {
        if (a->n_buffers < 3) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }
        const uint8_t* bytes = (const uint8_t*)a->buffers[2];
        int64_t s, e;
        if (acc->off64) {
            const int64_t* off = (const int64_t*)a->buffers[1];
            if (!off) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }
            s = off[i]; e = off[i + 1];
        } else {
            const int32_t* off = (const int32_t*)a->buffers[1];
            if (!off) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }
            s = off[i]; e = off[i + 1];
        }
        if (s < 0 || e < s) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }
        if (!grow_ba(acc, 1)) return;
        acc->ba[acc->nba].data = (uint8_t*)(bytes ? bytes + s : NULL);
        acc->ba[acc->nba].length = (uint32_t)(e - s);
        acc->nba++;
        return;
    }
    /* fixed-width */
    const uint8_t* d = (a->n_buffers >= 2) ? (const uint8_t*)a->buffers[1] : NULL;
    if (!d || acc->stride == 0) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }
    if (!grow_fixed(acc, (int64_t)acc->stride)) return;
    uint8_t* dst = acc->fixed + acc->fixed_bytes;
    if (acc->conv != CONV_NONE) {
        const uint8_t* src = d + (size_t)i * acc->src_stride;
        switch (acc->conv) {
        case CONV_I8:  { int32_t v = *(const int8_t*)src; memcpy(dst, &v, 4); break; }
        case CONV_U8:  { int32_t v = *src; memcpy(dst, &v, 4); break; }
        case CONV_I16: { int16_t t; memcpy(&t, src, 2); int32_t v = t; memcpy(dst, &v, 4); break; }
        case CONV_U16: { uint16_t t; memcpy(&t, src, 2); int32_t v = t; memcpy(dst, &v, 4); break; }
        case CONV_DATE64: {
            /* milliseconds since epoch -> days since epoch, floored so that
             * negative timestamps map to the calendar day they fall in. */
            int64_t ms; memcpy(&ms, src, 8);
            int64_t days = ms / 86400000;
            if (ms % 86400000 < 0) days--;
            if (days < INT32_MIN || days > INT32_MAX) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }
            int32_t v = (int32_t)days; memcpy(dst, &v, 4); break;
        }
        case CONV_TIME32_S: {
            int32_t sec; memcpy(&sec, src, 4);
            if (sec < 0 || sec > 86399) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }
            int32_t v = sec * 1000; memcpy(dst, &v, 4); break;
        }
        case CONV_TS_S: {
            int64_t sec; memcpy(&sec, src, 8);
            if (sec > INT64_MAX / 1000 || sec < INT64_MIN / 1000) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }
            int64_t v = sec * 1000; memcpy(dst, &v, 8); break;
        }
        case CONV_DECIMAL: {
            /* Keep the low `stride` bytes of the little-endian two's complement
             * value (the precision guarantees the value fits) and reverse them
             * into Parquet's big-endian order. */
            for (size_t b = 0; b < acc->stride; b++) dst[acc->stride - 1 - b] = src[b];
            break;
        }
        default: ctx->rc = CARQUET_ERROR_INTERNAL; return;
        }
        acc->fixed_bytes += (int64_t)acc->stride;
        return;
    }
    memcpy(dst, d + (size_t)i * acc->stride, acc->stride);
    acc->fixed_bytes += (int64_t)acc->stride;
}

/* Read the dictionary index at position `i` from a dictionary-encoded array's
 * index buffer (buffers[1]), widening it to int64. `index_fmt` is the Arrow
 * format of the index type (a single-character signed/unsigned integer code).
 * Returns false on an unsupported index type or a missing buffer. */
static bool read_dict_index(const char* index_fmt, const struct ArrowArray* a,
                            int64_t i, int64_t* out) {
    const void* buf = (a->n_buffers >= 2) ? a->buffers[1] : NULL;
    if (!buf || !index_fmt || !index_fmt[0] || index_fmt[1] != '\0') return false;
    switch (index_fmt[0]) {
    case 'c': *out = ((const int8_t*)buf)[i];             return true;  /* int8  */
    case 'C': *out = ((const uint8_t*)buf)[i];            return true;  /* uint8 */
    case 's': *out = ((const int16_t*)buf)[i];            return true;  /* int16 */
    case 'S': *out = ((const uint16_t*)buf)[i];           return true;  /* uint16 */
    case 'i': *out = ((const int32_t*)buf)[i];            return true;  /* int32 */
    case 'I': *out = ((const uint32_t*)buf)[i];           return true;  /* uint32 */
    case 'l': *out = ((const int64_t*)buf)[i];            return true;  /* int64 */
    case 'L': *out = (int64_t)((const uint64_t*)buf)[i];  return true;  /* uint64 */
    default:  return false;
    }
}

/* Emit one absent/empty placeholder entry (rep, def) for every leaf covered by
 * a subtree spanning columns [base, base+count). */
static void emit_placeholder_range(shred_ctx_t* ctx, int32_t base, int64_t count,
                                   int16_t rep, int16_t def) {
    for (int64_t c = 0; c < count; c++) {
        append_level(&ctx->accs[base + c], rep, def);
        if (ctx->accs[base + c].oom) ctx->rc = CARQUET_ERROR_OUT_OF_MEMORY;
    }
}

/*
 * Fast path for the overwhelmingly common case: a top-level, non-nested,
 * non-dictionary, fixed-width primitive leaf (a flat column of a flat table).
 * Instead of recursing once per row through visit()/emit_value() — a function
 * call plus a capacity check plus a 4/8-byte memcpy per element — bulk-copy the
 * contiguous Arrow values buffer and fill the definition levels in one pass.
 *
 * Returns false if (s, a) is not eligible, leaving the accumulator untouched so
 * the caller falls back to the generic row-by-row shredder. When it returns
 * true it has fully handled the column; the caller must then check ctx->rc for
 * an out-of-memory condition.
 */
static bool shred_flat_fast(shred_ctx_t* ctx, leaf_acc_t* acc,
                            const struct ArrowSchema* s,
                            const struct ArrowArray* a) {
    if (classify(s->format) != ANODE_PRIMITIVE) return false;
    if (s->dictionary) return false;                 /* indices need resolving  */
    if (acc->max_rep != 0) return false;             /* top-level (flat) only   */
    if (acc->is_bool || acc->is_bytearray) return false;  /* fixed-width only   */
    if (acc->stride == 0) return false;
    if (acc->conv != CONV_NONE) return false;        /* needs conversion, not memcpy */
    if (a->offset < 0) return false;                 /* visit() reports the error */
    const uint8_t* d = (a->n_buffers >= 2) ? (const uint8_t*)a->buffers[1] : NULL;
    if (!d) return false;

    int64_t n = a->length;
    if (n <= 0) return true;                         /* eligible, nothing to do */
    size_t stride = acc->stride;
    /* Sliced array: its buffers start `offset` logical elements before the
     * first element of this slice. */
    int64_t off = a->offset;
    d += (size_t)off * stride;
    bool nullable = (s->flags & ARROW_FLAG_NULLABLE) != 0;
    const uint8_t* valid = (nullable && a->n_buffers >= 1)
                               ? (const uint8_t*)a->buffers[0] : NULL;

    if (!grow_levels(acc, n)) { ctx->rc = CARQUET_ERROR_OUT_OF_MEMORY; return true; }
    if (!grow_fixed(acc, n * (int64_t)stride)) { ctx->rc = CARQUET_ERROR_OUT_OF_MEMORY; return true; }

    int16_t* def = acc->def;
    int64_t lv = acc->nlevels;
    uint8_t* dst = acc->fixed + acc->fixed_bytes;
    /* rep is left untouched: max_rep == 0 means write_batch is passed rep = NULL. */

    if (!valid) {
        /* No validity bitmap -> every row present: one bulk copy. */
        memcpy(dst, d, (size_t)n * stride);
        acc->fixed_bytes += n * (int64_t)stride;
        int16_t defp = acc->max_def;                 /* 0 (REQUIRED) or 1 (OPTIONAL) */
        if (defp != 0) {
            for (int64_t r = 0; r < n; r++) def[lv + r] = defp;
        }
    } else {
        /* Mixed nulls: copy present values densely, def = present ? 1 : 0.
         * `valid` non-NULL implies nullable, hence max_def >= 1 and a real def
         * array is written. */
        int64_t nb = 0;
        for (int64_t r = 0; r < n; r++) {
            int present = (int)ARROW_VALID(valid, off + r);
            def[lv + r] = (int16_t)present;
            if (present) {
                memcpy(dst + (size_t)nb * stride, d + (size_t)r * stride, stride);
                nb++;
            }
        }
        acc->fixed_bytes += nb * (int64_t)stride;
    }
    acc->nlevels = lv + n;
    return true;
}

/*
 * Visit element `idx` of the Arrow node (s, a), appending shredded levels and
 * values to the leaf accumulators.
 *
 *   def_in    definition level already accounted for by present ancestors
 *   rep_in    repetition level to stamp on the first leaf value produced here
 *   rep_depth number of repeated groups entered so far (rep level of the
 *             innermost repeated ancestor)
 *   base      column index of the leftmost leaf under this node
 */
static void visit(shred_ctx_t* ctx, const struct ArrowSchema* s,
                  const struct ArrowArray* a, int64_t idx,
                  int16_t def_in, int16_t rep_in, int16_t rep_depth,
                  int32_t base, int depth) {
    if (ctx->rc != CARQUET_OK) return;
    if (depth > ARROW_MAX_NEST_DEPTH) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }
    if (a->offset < 0 || idx < 0 || idx >= a->length) {
        ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return;
    }
    /* Sliced array: `idx` is logical (relative to the slice); every buffer of
     * this node is indexed by the physical position `offset + idx`. Children
     * of a struct are indexed by the same physical position (the Arrow C Data
     * Interface applies a struct's offset to its children), while a list's
     * offsets buffer already yields logical child positions. */
    int64_t pidx = a->offset + idx;

    anode_kind_t kind = classify(s->format);
    bool nullable = (s->flags & ARROW_FLAG_NULLABLE) != 0;
    int64_t subtree_leaves = arrow_leaf_count(s);

    /* Absent: record one null placeholder at def_in for every leaf below. */
    if (nullable && arr_is_null(a, pidx)) {
        emit_placeholder_range(ctx, base, subtree_leaves, rep_in, def_in);
        return;
    }
    int16_t def_present = (int16_t)(def_in + (nullable ? 1 : 0));

    switch (kind) {
    case ANODE_PRIMITIVE: {
        leaf_acc_t* acc = &ctx->accs[base];
        append_level(acc, rep_in, def_present);
        if (acc->oom) { ctx->rc = CARQUET_ERROR_OUT_OF_MEMORY; return; }
        if (s->dictionary) {
            /* Dictionary-encoded: `a` holds indices, `a->dictionary` the values.
             * `acc` is configured for the (resolved) value type, so emit the
             * looked-up value. Null entries were already handled above via the
             * index array's validity bitmap. */
            const struct ArrowArray* dict = a->dictionary;
            int64_t k;
            if (!dict || dict->offset < 0 || !read_dict_index(s->format, a, pidx, &k) ||
                k < 0 || k >= dict->length) {
                ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT;
                return;
            }
            int64_t pk = dict->offset + k;
            if (arr_is_null(dict, pk)) {
                /* The dictionary slot itself is null: the row is null. Replace
                 * the level just appended with the absent placeholder. */
                acc->def[acc->nlevels - 1] = def_in;
                return;
            }
            emit_value(ctx, acc, dict, pk);
        } else {
            emit_value(ctx, acc, a, pidx);
        }
        return;
    }
    case ANODE_STRUCT: {
        int32_t child_base = base;
        for (int64_t c = 0; c < s->n_children; c++) {
            if (!a->children || !a->children[c] || !s->children[c]) {
                ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return;
            }
            visit(ctx, s->children[c], a->children[c], pidx,
                  def_present, rep_in, rep_depth, child_base, depth + 1);
            child_base += (int32_t)arrow_leaf_count(s->children[c]);
        }
        return;
    }
    case ANODE_LIST:
    case ANODE_MAP: {
        if (s->n_children != 1 || !a->children || !a->children[0] || !s->children[0]) {
            ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return;
        }
        int64_t lo, hi;
        int64_t fixed_n = fixed_list_size(s->format);
        if (fixed_n > 0) {
            /* "+w:n": no offsets buffer, slot `pidx` owns [pidx*n, pidx*n + n). */
            if (pidx > INT64_MAX / fixed_n) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }
            lo = pidx * fixed_n; hi = lo + fixed_n;
        } else {
            bool large = (s->format[1] == 'L');
            if (a->n_buffers < 2 || !a->buffers[1]) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }
            if (large) {
                const int64_t* off = (const int64_t*)a->buffers[1];
                lo = off[pidx]; hi = off[pidx + 1];
            } else {
                const int32_t* off = (const int32_t*)a->buffers[1];
                lo = off[pidx]; hi = off[pidx + 1];
            }
        }
        if (lo < 0 || hi < lo) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }

        if (hi == lo) {
            /* Present but empty: one placeholder at def_present per leaf below. */
            emit_placeholder_range(ctx, base, subtree_leaves, rep_in, def_present);
            return;
        }
        int16_t rep_child = (int16_t)(rep_depth + 1);
        const struct ArrowSchema* cs = s->children[0];
        const struct ArrowArray* ca = a->children[0];
        if (ca->length < hi) { ctx->rc = CARQUET_ERROR_INVALID_ARGUMENT; return; }
        for (int64_t j = lo; j < hi; j++) {
            int16_t r = (j == lo) ? rep_in : rep_child;
            visit(ctx, cs, ca, j, (int16_t)(def_present + 1), r, rep_child,
                  base, depth + 1);
            if (ctx->rc != CARQUET_OK) return;
        }
        return;
    }
    default:
        ctx->rc = CARQUET_ERROR_NOT_IMPLEMENTED;
        return;
    }
}

/* Pre-walk: assign leaf column indices in DFS order and populate each leaf
 * accumulator's type descriptor. Mirrors the DFS order the writer's schema was
 * built in, so column N corresponds to accs[N]. */
static carquet_status_t assign_leaves(const struct ArrowSchema* s, int32_t* next_col,
                                      leaf_acc_t* accs, int32_t cap,
                                      int16_t cur_def, int16_t cur_rep,
                                      int depth, carquet_error_t* error) {
    if (depth > ARROW_MAX_NEST_DEPTH || !s || !s->format) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "Arrow import: bad schema tree");
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    bool nullable = (s->flags & ARROW_FLAG_NULLABLE) != 0;
    switch (classify(s->format)) {
    case ANODE_PRIMITIVE: {
        if (*next_col >= cap) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "Arrow import: more leaves than writer columns");
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        /* Dictionary-encoded leaf: configure the accumulator for the dictionary
         * VALUE type (indices are resolved to values during shredding), matching
         * the column type import_field() derived from the same dictionary. */
        const char* value_format = s->dictionary ? s->dictionary->format : s->format;
        if (s->dictionary && (!value_format || classify(value_format) != ANODE_PRIMITIVE)) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
                "Arrow import: only primitive dictionary values are supported");
            return CARQUET_ERROR_NOT_IMPLEMENTED;
        }
        carquet_physical_type_t pt;
        carquet_logical_type_t lt;
        bool has_lt;
        int32_t tlen;
        arrow_conv_t conv;
        size_t src_stride;
        carquet_status_t rc = parse_arrow_format(value_format, &pt, &lt, &has_lt, &tlen,
                                                 &conv, &src_stride);
        if (rc != CARQUET_OK) {
            CARQUET_SET_ERROR(error, rc, "Arrow import: unsupported format \"%s\"", value_format);
            return rc;
        }
        leaf_acc_t* acc = &accs[*next_col];
        acc->pt = pt;
        acc->is_bool = (pt == CARQUET_PHYSICAL_BOOLEAN);
        acc->is_bytearray = (pt == CARQUET_PHYSICAL_BYTE_ARRAY);
        acc->off64 = (value_format[0] == 'U' || value_format[0] == 'Z');
        acc->stride = fixed_stride(pt, tlen);
        acc->conv = conv;
        acc->src_stride = src_stride;
        acc->max_def = (int16_t)(cur_def + (nullable ? 1 : 0));
        acc->max_rep = cur_rep;
        (*next_col)++;
        return CARQUET_OK;
    }
    case ANODE_STRUCT: {
        int16_t d = (int16_t)(cur_def + (nullable ? 1 : 0));
        for (int64_t i = 0; i < s->n_children; i++) {
            carquet_status_t rc = assign_leaves(s->children[i], next_col, accs, cap,
                                                d, cur_rep, depth + 1, error);
            if (rc != CARQUET_OK) return rc;
        }
        return CARQUET_OK;
    }
    case ANODE_LIST:
    case ANODE_MAP:
        if (s->n_children != 1) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "Arrow import: malformed nested schema");
            return CARQUET_ERROR_INVALID_ARGUMENT;
        }
        /* outer nullability (+def) plus the implicit REPEATED group (+def,+rep) */
        return assign_leaves(s->children[0], next_col, accs, cap,
                             (int16_t)(cur_def + (nullable ? 1 : 0) + 1),
                             (int16_t)(cur_rep + 1), depth + 1, error);
    default:
        CARQUET_SET_ERROR(error, CARQUET_ERROR_NOT_IMPLEMENTED,
            "Arrow import: unsupported type \"%s\"", s->format);
        return CARQUET_ERROR_NOT_IMPLEMENTED;
    }
}

carquet_status_t carquet_writer_write_arrow(
    carquet_writer_t* writer,
    struct ArrowArray* array,
    struct ArrowSchema* schema,
    carquet_error_t* error) {

    if (!writer || !array || !schema) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT, "NULL argument");
        release_array(array);
        release_schema(schema);
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    carquet_status_t rc = CARQUET_OK;
    leaf_acc_t* accs = NULL;
    int32_t num_cols = carquet_writer_num_columns(writer);

    if (!schema->format || strcmp(schema->format, "+s") != 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "Arrow import: top-level schema must be a struct (\"+s\")");
        rc = CARQUET_ERROR_INVALID_ARGUMENT; goto done;
    }
    if (array->offset < 0 || array->length < 0) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "Arrow import: negative array offset or length");
        rc = CARQUET_ERROR_INVALID_ARGUMENT; goto done;
    }
    if (array->n_children != schema->n_children) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
            "Arrow import: array/schema child count mismatch (%lld vs %lld)",
            (long long)array->n_children, (long long)schema->n_children);
        rc = CARQUET_ERROR_INVALID_ARGUMENT; goto done;
    }

    /* Assign leaf columns and validate the leaf count equals the writer's. */
    accs = (leaf_acc_t*)calloc((size_t)(num_cols > 0 ? num_cols : 1), sizeof(leaf_acc_t));
    if (!accs) { rc = CARQUET_ERROR_OUT_OF_MEMORY; goto done; }
    {
        int32_t next_col = 0;
        for (int64_t i = 0; i < schema->n_children; i++) {
            rc = assign_leaves(schema->children[i], &next_col, accs, num_cols, 0, 0, 0, error);
            if (rc != CARQUET_OK) goto done;
        }
        if (next_col != num_cols) {
            CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                "Arrow import: schema has %d leaf columns, writer expects %d",
                (int)next_col, (int)num_cols);
            rc = CARQUET_ERROR_INVALID_ARGUMENT; goto done;
        }
    }

    /* Shred: for each top-level field, walk every record. */
    {
        shred_ctx_t ctx = { accs, num_cols, error, CARQUET_OK };
        int32_t base = 0;
        for (int64_t i = 0; i < schema->n_children; i++) {
            struct ArrowArray* carray = array->children[i];
            struct ArrowSchema* cschema = schema->children[i];
            if (!carray || !cschema) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                    "Arrow import: null child at %lld", (long long)i);
                rc = CARQUET_ERROR_INVALID_ARGUMENT; goto done;
            }
            /* A struct's offset applies to its children: child row r lives at
             * logical child position array->offset + r. */
            if (carray->offset < 0 || carray->length < array->offset + array->length) {
                CARQUET_SET_ERROR(error, CARQUET_ERROR_INVALID_ARGUMENT,
                    "Arrow import: child %lld length %lld < struct offset %lld + length %lld",
                    (long long)i, (long long)carray->length,
                    (long long)array->offset, (long long)array->length);
                rc = CARQUET_ERROR_INVALID_ARGUMENT; goto done;
            }
            struct ArrowArray view = *carray;
            view.offset = carray->offset + array->offset;
            view.length = array->length;
            if (!shred_flat_fast(&ctx, &accs[base], cschema, &view)) {
                for (int64_t r = 0; r < array->length; r++) {
                    visit(&ctx, cschema, &view, r, 0, 0, 0, base, 0);
                    if (ctx.rc != CARQUET_OK) break;
                }
            }
            if (ctx.rc != CARQUET_OK) {
                rc = ctx.rc;
                CARQUET_SET_ERROR(error, rc,
                    "Arrow import: shredding failed for field %lld", (long long)i);
                goto done;
            }
            base += (int32_t)arrow_leaf_count(cschema);
        }
    }

    /* Write each leaf column once, in column order. */
    for (int32_t c = 0; c < num_cols; c++) {
        leaf_acc_t* acc = &accs[c];
        const void* values = acc->is_bytearray ? (const void*)acc->ba
                                               : (const void*)acc->fixed;
        if (!values) values = "";  /* non-NULL sentinel; present count is 0 */
        /* Pass level arrays only where the schema carries them: a REQUIRED flat
         * leaf must be written with def == rep == NULL. */
        const int16_t* def = (acc->max_def > 0) ? acc->def : NULL;
        const int16_t* rep = (acc->max_rep > 0) ? acc->rep : NULL;
        rc = carquet_writer_write_batch(writer, c, values, acc->nlevels, def, rep);
        if (rc != CARQUET_OK) {
            CARQUET_SET_ERROR(error, rc,
                "Arrow import: write_batch failed for column %d", (int)c);
            goto done;
        }
    }

done:
    if (accs) {
        for (int32_t c = 0; c < num_cols; c++) {
            free(accs[c].def);
            free(accs[c].rep);
            free(accs[c].fixed);
            free(accs[c].ba);
        }
        free(accs);
    }
    release_array(array);
    release_schema(schema);
    return rc;
}
