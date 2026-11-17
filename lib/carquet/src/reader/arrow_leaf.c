/**
 * @file arrow_leaf.c
 * @brief Parquet leaf -> Arrow value mapping (see arrow_leaf.h).
 */
#include "arrow_leaf.h"

#include <stdio.h>
#include <string.h>

#include "core/endian.h"

enum {
    CV_NARROW,          /* INT32 INTEGER(8|16) -> 1 or 2 bytes */
    CV_DECIMAL_INT32,   /* sign-extended into dst_width */
    CV_DECIMAL_INT64,
    CV_DECIMAL_FLBA,    /* big-endian two's complement, src_width <= dst_width */
    CV_DECIMAL_BYTES,   /* big-endian two's complement of any length */
    CV_INT96_NANOS      /* nanos of day (LE int64) + Julian day (LE uint32) */
};

#define JULIAN_UNIX_EPOCH 2440588
#define NS_PER_DAY 86400000000000ULL

static size_t physical_width(carquet_physical_type_t pt, int32_t type_length) {
    switch (pt) {
    case CARQUET_PHYSICAL_BOOLEAN: return 1;
    case CARQUET_PHYSICAL_INT32: case CARQUET_PHYSICAL_FLOAT: return 4;
    case CARQUET_PHYSICAL_INT64: case CARQUET_PHYSICAL_DOUBLE: return 8;
    case CARQUET_PHYSICAL_INT96: return 12;
    case CARQUET_PHYSICAL_BYTE_ARRAY: return sizeof(carquet_byte_array_t);
    case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY: return type_length > 0 ? (size_t)type_length : 0;
    default: return 0;
    }
}

/* decimal128 up to 38 digits, decimal256 up to 76; 0 when Arrow has none. */
static size_t decimal_width(const carquet_logical_type_t* lt) {
    int32_t p = lt->params.decimal.precision;
    if (p < 1) return 0;
    if (p <= 38) return 16;
    if (p <= 76) return 32;
    return 0;
}

static void set_convert(carquet_arrow_leaf_t* out, uint8_t cv, size_t dst_width) {
    out->kind = CARQUET_ARROW_LEAF_CONVERT;
    out->convert = cv;
    out->dst_width = dst_width;
}

carquet_status_t carquet_arrow_leaf_plan(carquet_physical_type_t pt, int32_t type_length,
                                         const carquet_logical_type_t* lt,
                                         carquet_arrow_leaf_t* out) {
    memset(out, 0, sizeof(*out));
    out->physical = pt;
    out->src_width = physical_width(pt, type_length);
    out->dst_width = out->src_width;
    out->kind = CARQUET_ARROW_LEAF_COPY;
    char* f = out->format;
    const size_t fs = sizeof(out->format);
    carquet_logical_type_id_t lid = lt ? lt->id : CARQUET_LOGICAL_UNKNOWN;

    if (lid == CARQUET_LOGICAL_DECIMAL &&
        (pt == CARQUET_PHYSICAL_INT32 || pt == CARQUET_PHYSICAL_INT64 ||
         pt == CARQUET_PHYSICAL_BYTE_ARRAY ||
         (pt == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY && type_length > 0))) {
        size_t w = decimal_width(lt);
        /* An FLBA wider than the Arrow decimal keeps its storage type. */
        if (w && !(pt == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY && (size_t)type_length > w)) {
            snprintf(f, fs, w == 16 ? "d:%d,%d" : "d:%d,%d,256",
                     (int)lt->params.decimal.precision, (int)lt->params.decimal.scale);
            set_convert(out, pt == CARQUET_PHYSICAL_INT32 ? CV_DECIMAL_INT32
                           : pt == CARQUET_PHYSICAL_INT64 ? CV_DECIMAL_INT64
                           : pt == CARQUET_PHYSICAL_BYTE_ARRAY ? CV_DECIMAL_BYTES
                           : CV_DECIMAL_FLBA, w);
            return CARQUET_OK;
        }
    }

    if (lid == CARQUET_LOGICAL_NULL) {
        snprintf(f, fs, "n");
        out->kind = CARQUET_ARROW_LEAF_NULL;
        out->dst_width = 0;
        return CARQUET_OK;
    }

    switch (pt) {
    case CARQUET_PHYSICAL_BOOLEAN:
        snprintf(f, fs, "b");
        out->kind = CARQUET_ARROW_LEAF_BOOL;
        out->dst_width = 0;
        return CARQUET_OK;

    case CARQUET_PHYSICAL_INT32:
        if (lid == CARQUET_LOGICAL_DATE) {
            snprintf(f, fs, "tdD");   /* date32[days] */
        } else if (lid == CARQUET_LOGICAL_TIME) {
            snprintf(f, fs, "ttm");   /* time32[ms] */
        } else if (lid == CARQUET_LOGICAL_INTEGER) {
            int bw = lt->params.integer.bit_width;
            bool s = lt->params.integer.is_signed;
            if (bw == 8 || bw == 16) {
                snprintf(f, fs, bw == 8 ? (s ? "c" : "C") : (s ? "s" : "S"));
                set_convert(out, CV_NARROW, (size_t)bw / 8);
            } else {
                snprintf(f, fs, s ? "i" : "I");
            }
        } else {
            snprintf(f, fs, "i");
        }
        return CARQUET_OK;

    case CARQUET_PHYSICAL_INT64:
        if (lid == CARQUET_LOGICAL_TIMESTAMP) {
            char u = lt->params.timestamp.unit == CARQUET_TIME_UNIT_MILLIS ? 'm'
                   : lt->params.timestamp.unit == CARQUET_TIME_UNIT_MICROS ? 'u' : 'n';
            snprintf(f, fs, "ts%c:%s", u, lt->params.timestamp.is_adjusted_to_utc ? "UTC" : "");
        } else if (lid == CARQUET_LOGICAL_TIME) {
            snprintf(f, fs, "tt%c", lt->params.time.unit == CARQUET_TIME_UNIT_MICROS ? 'u' : 'n');
        } else if (lid == CARQUET_LOGICAL_INTEGER) {
            snprintf(f, fs, lt->params.integer.is_signed ? "l" : "L");
        } else {
            snprintf(f, fs, "l");
        }
        return CARQUET_OK;

    case CARQUET_PHYSICAL_INT96:
        /* As pyarrow: timestamp[ns], no time zone. */
        snprintf(f, fs, "tsn:");
        set_convert(out, CV_INT96_NANOS, 8);
        return CARQUET_OK;

    case CARQUET_PHYSICAL_FLOAT:
        snprintf(f, fs, "f");
        return CARQUET_OK;

    case CARQUET_PHYSICAL_DOUBLE:
        snprintf(f, fs, "g");
        return CARQUET_OK;

    case CARQUET_PHYSICAL_BYTE_ARRAY:
        snprintf(f, fs, (lid == CARQUET_LOGICAL_STRING || lid == CARQUET_LOGICAL_ENUM ||
                         lid == CARQUET_LOGICAL_JSON) ? "u" : "z");
        out->kind = CARQUET_ARROW_LEAF_BINARY;
        out->dst_width = 0;
        return CARQUET_OK;

    case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
        if (type_length <= 0) return CARQUET_ERROR_INVALID_ARGUMENT;
        if (lid == CARQUET_LOGICAL_FLOAT16 && type_length == 2) {
            snprintf(f, fs, "e");     /* halffloat */
        } else {
            snprintf(f, fs, "w:%d", (int)type_length);
        }
        return CARQUET_OK;

    default:
        return CARQUET_ERROR_NOT_IMPLEMENTED;
    }
}

/* Write the little-endian two's complement `le` (dst_width bytes) to dst in
 * the host's byte order, which is how Arrow lays out decimals. */
static void store_native(uint8_t* dst, const uint8_t* le, size_t width) {
#if CARQUET_LITTLE_ENDIAN
    memcpy(dst, le, width);
#else
    for (size_t i = 0; i < width; i++) dst[i] = le[width - 1 - i];
#endif
}

/* Big-endian two's complement `be` (len bytes) -> dst (width bytes, native).
 * Leading bytes that only repeat the sign are dropped; false if the value
 * still does not fit. */
static bool store_be_decimal(uint8_t* dst, const uint8_t* be, size_t len, size_t width) {
    uint8_t le[32];
    uint8_t sign = (len > 0 && (be[0] & 0x80)) ? 0xFF : 0x00;
    while (len > width && be[0] == sign && ((be[1] & 0x80) == (sign & 0x80))) {
        be++;
        len--;
    }
    if (len > width) {
        memset(dst, 0, width);
        return false;
    }
    memset(le, sign, width);
    for (size_t i = 0; i < len; i++) le[i] = be[len - 1 - i];
    store_native(dst, le, width);
    return true;
}

static void store_int_decimal(uint8_t* dst, int64_t v, size_t width) {
    uint8_t le[32];
    uint64_t u = (uint64_t)v;
    for (size_t i = 0; i < 8; i++) le[i] = (uint8_t)(u >> (8 * i));
    memset(le + 8, v < 0 ? 0xFF : 0x00, width - 8);
    store_native(dst, le, width);
}

bool carquet_arrow_leaf_store(const carquet_arrow_leaf_t* plan, uint8_t* dst, const void* src) {
    if (plan->kind == CARQUET_ARROW_LEAF_COPY) {
        memcpy(dst, src, plan->dst_width);
        return true;
    }
    switch (plan->convert) {
    case CV_NARROW: {
        int32_t v;
        memcpy(&v, src, sizeof(v));
        if (plan->dst_width == 1) {
            dst[0] = (uint8_t)v;
        } else {
            uint16_t t = (uint16_t)v;   /* two's complement: serves int16 and uint16 */
            memcpy(dst, &t, sizeof(t));
        }
        return true;
    }
    case CV_DECIMAL_INT32: {
        int32_t v;
        memcpy(&v, src, sizeof(v));
        store_int_decimal(dst, v, plan->dst_width);
        return true;
    }
    case CV_DECIMAL_INT64: {
        int64_t v;
        memcpy(&v, src, sizeof(v));
        store_int_decimal(dst, v, plan->dst_width);
        return true;
    }
    case CV_DECIMAL_FLBA:
        return store_be_decimal(dst, (const uint8_t*)src, plan->src_width, plan->dst_width);
    case CV_DECIMAL_BYTES: {
        carquet_byte_array_t ba;
        memcpy(&ba, src, sizeof(ba));
        if (ba.length <= 0 || !ba.data) {
            memset(dst, 0, plan->dst_width);
            return true;
        }
        return store_be_decimal(dst, ba.data, (size_t)ba.length, plan->dst_width);
    }
    case CV_INT96_NANOS: {
        const uint8_t* b = (const uint8_t*)src;
        uint64_t nanos = 0;
        uint32_t day = 0;
        for (int i = 7; i >= 0; i--) nanos = (nanos << 8) | b[i];
        for (int i = 11; i >= 8; i--) day = (day << 8) | b[i];
        /* Unsigned arithmetic wraps like pyarrow's int64 for out-of-range days. */
        uint64_t ns = (uint64_t)((int64_t)day - JULIAN_UNIX_EPOCH) * NS_PER_DAY + nanos;
        int64_t out;
        memcpy(&out, &ns, sizeof(out));
        memcpy(dst, &out, sizeof(out));
        return true;
    }
    default:
        memset(dst, 0, plan->dst_width);
        return false;
    }
}
