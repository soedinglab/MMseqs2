/**
 * @file thrift_decode.c
 * @brief Thrift Compact Protocol decoder implementation
 */

#include "core/allocator.h"
#include "thrift_decode.h"
#include "core/endian.h"
#include <stdlib.h>
#include <string.h>

/* ============================================================================
 * Internal Helpers
 * ============================================================================
 */

static void set_error(thrift_decoder_t* dec, carquet_status_t status, const char* msg) {
    if (dec->status == CARQUET_OK) {
        dec->status = status;
        if (msg) {
            strncpy(dec->error_message, msg, sizeof(dec->error_message) - 1);
            dec->error_message[sizeof(dec->error_message) - 1] = '\0';
        }
    }
}

static inline bool has_bytes(thrift_decoder_t* dec, size_t n) {
    return carquet_buffer_reader_has(&dec->reader, n);
}

static inline uint8_t read_byte_raw(thrift_decoder_t* dec) {
    if (!has_bytes(dec, 1)) {
        set_error(dec, CARQUET_ERROR_THRIFT_TRUNCATED, "Unexpected end of data");
        return 0;
    }
    uint8_t b;
    carquet_buffer_reader_read_byte(&dec->reader, &b);
    return b;
}

/* ============================================================================
 * Decoder Lifecycle
 * ============================================================================
 */

void thrift_decoder_init(thrift_decoder_t* dec, const uint8_t* data, size_t size) {
    memset(dec, 0, sizeof(*dec));
    carquet_buffer_reader_init_data(&dec->reader, data, size);
    dec->nesting_level = 0;
    dec->bool_pending = false;
    dec->status = CARQUET_OK;
}

void thrift_decoder_init_reader(thrift_decoder_t* dec,
                                 const carquet_buffer_reader_t* reader) {
    memset(dec, 0, sizeof(*dec));
    dec->reader = *reader;
    dec->nesting_level = 0;
    dec->bool_pending = false;
    dec->status = CARQUET_OK;
}

/* ============================================================================
 * Varint Reading
 * ============================================================================
 */

uint64_t thrift_read_varint(thrift_decoder_t* dec) {
    uint64_t result = 0;
    int shift = 0;

    while (shift < 64) {
        if (!has_bytes(dec, 1)) {
            set_error(dec, CARQUET_ERROR_THRIFT_TRUNCATED, "Truncated varint");
            return 0;
        }

        uint8_t byte = read_byte_raw(dec);
        result |= (uint64_t)(byte & 0x7F) << shift;

        if ((byte & 0x80) == 0) {
            return result;
        }

        shift += 7;
    }

    set_error(dec, CARQUET_ERROR_THRIFT_DECODE, "Varint overflow");
    return 0;
}

int64_t thrift_read_zigzag(thrift_decoder_t* dec) {
    uint64_t n = thrift_read_varint(dec);
    return carquet_zigzag_decode64(n);
}

/* ============================================================================
 * Primitive Reading
 * ============================================================================
 */

int8_t thrift_read_byte(thrift_decoder_t* dec) {
    return (int8_t)read_byte_raw(dec);
}

int16_t thrift_read_i16(thrift_decoder_t* dec) {
    return (int16_t)thrift_read_zigzag(dec);
}

int32_t thrift_read_i32(thrift_decoder_t* dec) {
    return (int32_t)thrift_read_zigzag(dec);
}

int64_t thrift_read_i64(thrift_decoder_t* dec) {
    return thrift_read_zigzag(dec);
}

double thrift_read_double(thrift_decoder_t* dec) {
    if (!has_bytes(dec, 8)) {
        set_error(dec, CARQUET_ERROR_THRIFT_TRUNCATED, "Truncated double");
        return 0.0;
    }

    /* Doubles are stored as 8 bytes, little-endian */
    const uint8_t* p = carquet_buffer_reader_peek(&dec->reader);
    double result = carquet_read_f64_le(p);
    carquet_buffer_reader_skip(&dec->reader, 8);
    return result;
}

bool thrift_read_bool(thrift_decoder_t* dec) {
    /* If we have a pending boolean from field header, use it */
    if (dec->bool_pending) {
        dec->bool_pending = false;
        return dec->bool_value;
    }

    /* Otherwise read a byte */
    return read_byte_raw(dec) == 1;
}

const uint8_t* thrift_read_binary(thrift_decoder_t* dec, int32_t* length) {
    /* Length is a varint */
    int32_t len = (int32_t)thrift_read_varint(dec);

    if (len < 0) {
        set_error(dec, CARQUET_ERROR_THRIFT_DECODE, "Negative binary length");
        *length = 0;
        return NULL;
    }

    if (!has_bytes(dec, (size_t)len)) {
        set_error(dec, CARQUET_ERROR_THRIFT_TRUNCATED, "Truncated binary data");
        *length = 0;
        return NULL;
    }

    *length = len;
    const uint8_t* data = carquet_buffer_reader_peek(&dec->reader);
    carquet_buffer_reader_skip(&dec->reader, (size_t)len);
    return data;
}

char* thrift_read_string_alloc(thrift_decoder_t* dec) {
    int32_t length;
    const uint8_t* data = thrift_read_binary(dec, &length);

    if (!data && length == 0 && dec->status == CARQUET_OK) {
        /* Empty string */
        char* str = (char*)carquet_mem_malloc(1);
        if (str) str[0] = '\0';
        return str;
    }

    if (!data) {
        return NULL;
    }

    char* str = (char*)carquet_mem_malloc((size_t)length + 1);
    if (!str) {
        set_error(dec, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate string");
        return NULL;
    }

    memcpy(str, data, (size_t)length);
    str[length] = '\0';
    return str;
}

void thrift_read_uuid(thrift_decoder_t* dec, uint8_t uuid[16]) {
    if (!has_bytes(dec, 16)) {
        set_error(dec, CARQUET_ERROR_THRIFT_TRUNCATED, "Truncated UUID");
        memset(uuid, 0, 16);
        return;
    }

    carquet_buffer_reader_read(&dec->reader, uuid, 16);
}

/* ============================================================================
 * Struct Reading
 * ============================================================================
 */

void thrift_read_struct_begin(thrift_decoder_t* dec) {
    if (dec->nesting_level >= THRIFT_MAX_NESTING) {
        set_error(dec, CARQUET_ERROR_THRIFT_DECODE, "Struct nesting too deep");
        return;
    }

    dec->last_field_id[dec->nesting_level] = 0;
    dec->nesting_level++;
}

void thrift_read_struct_end(thrift_decoder_t* dec) {
    if (dec->nesting_level > 0) {
        dec->nesting_level--;
    }
}

bool thrift_read_field_begin(thrift_decoder_t* dec,
                              thrift_type_t* type,
                              int16_t* field_id) {
    if (dec->status != CARQUET_OK) {
        *type = THRIFT_TYPE_STOP;
        *field_id = 0;
        return false;
    }

    uint8_t header = read_byte_raw(dec);

    if (header == 0) {
        /* STOP field */
        *type = THRIFT_TYPE_STOP;
        *field_id = 0;
        return false;
    }

    /* Lower 4 bits are the type */
    *type = (thrift_type_t)(header & 0x0F);

    /* Upper 4 bits are the field ID delta (if non-zero) */
    int16_t delta = (header >> 4) & 0x0F;

    int16_t prev_field_id = 0;
    if (dec->nesting_level > 0) {
        prev_field_id = dec->last_field_id[dec->nesting_level - 1];
    }

    if (delta == 0) {
        /* Field ID is encoded as zigzag varint */
        *field_id = thrift_read_i16(dec);
    } else {
        /* Field ID is delta from previous */
        *field_id = prev_field_id + delta;
    }

    /* Update last field ID */
    if (dec->nesting_level > 0) {
        dec->last_field_id[dec->nesting_level - 1] = *field_id;
    }

    /* Handle embedded boolean values */
    if (*type == THRIFT_TYPE_TRUE) {
        dec->bool_pending = true;
        dec->bool_value = true;
        *type = THRIFT_TYPE_TRUE;
    } else if (*type == THRIFT_TYPE_FALSE) {
        dec->bool_pending = true;
        dec->bool_value = false;
        *type = THRIFT_TYPE_FALSE;
    }

    return true;
}

void thrift_skip_field(thrift_decoder_t* dec, thrift_type_t type) {
    thrift_skip(dec, type);
}

/* ============================================================================
 * Container Reading
 * ============================================================================
 */

void thrift_read_list_begin(thrift_decoder_t* dec,
                             thrift_type_t* elem_type,
                             int32_t* count) {
    uint8_t header = read_byte_raw(dec);

    /* Lower 4 bits are element type */
    *elem_type = (thrift_type_t)(header & 0x0F);

    /* Upper 4 bits are size if <= 14 */
    int32_t size_and_type = (header >> 4) & 0x0F;

    if (size_and_type == 0x0F) {
        /* Size is encoded as a separate varint */
        *count = (int32_t)thrift_read_varint(dec);
    } else {
        *count = size_and_type;
    }

    if (*count < 0) {
        set_error(dec, CARQUET_ERROR_THRIFT_DECODE, "Negative list size");
        *count = 0;
        return;
    }

    /* Each list element consumes at least 1 byte, so count cannot exceed
     * remaining data.  This prevents billion-iteration busy loops from
     * malicious varints in tiny payloads.
     *
     * Report this as TRUNCATED, not DECODE: callers that parse out of a
     * growable window (see read_and_parse_page_header_fread) only enlarge the
     * window on TRUNCATED, so misclassifying it turns "need more bytes" into a
     * hard parse failure on an otherwise valid structure. */
    size_t remaining = carquet_buffer_reader_remaining(&dec->reader);
    if ((size_t)*count > remaining) {
        set_error(dec, CARQUET_ERROR_THRIFT_TRUNCATED, "List count exceeds remaining data");
        *count = 0;
    }
}

void thrift_read_set_begin(thrift_decoder_t* dec,
                            thrift_type_t* elem_type,
                            int32_t* count) {
    /* Set has the same encoding as list */
    thrift_read_list_begin(dec, elem_type, count);
}

void thrift_read_map_begin(thrift_decoder_t* dec,
                            thrift_type_t* key_type,
                            thrift_type_t* value_type,
                            int32_t* count) {
    /* Size first */
    *count = (int32_t)thrift_read_varint(dec);

    if (*count < 0) {
        set_error(dec, CARQUET_ERROR_THRIFT_DECODE, "Negative map size");
        *count = 0;
        *key_type = THRIFT_TYPE_STOP;
        *value_type = THRIFT_TYPE_STOP;
        return;
    }

    if (*count == 0) {
        *key_type = THRIFT_TYPE_STOP;
        *value_type = THRIFT_TYPE_STOP;
        return;
    }

    /* Each map entry consumes at least 1 byte (types byte already read
     * separately), so count cannot exceed remaining data.  This prevents
     * billion-iteration busy loops from malicious varints. */
    size_t remaining = carquet_buffer_reader_remaining(&dec->reader);
    if ((size_t)*count > remaining) {
        set_error(dec, CARQUET_ERROR_THRIFT_TRUNCATED, "Map count exceeds remaining data");
        *count = 0;
        *key_type = THRIFT_TYPE_STOP;
        *value_type = THRIFT_TYPE_STOP;
        return;
    }

    /* Key and value types in one byte */
    uint8_t types = read_byte_raw(dec);
    *key_type = (thrift_type_t)((types >> 4) & 0x0F);
    *value_type = (thrift_type_t)(types & 0x0F);
}

/* ============================================================================
 * Skip Functions
 * ============================================================================
 */

void thrift_skip(thrift_decoder_t* dec, thrift_type_t type);

/* Skip one container element. Only a struct field carries a boolean inside
 * its header; in a list, set or map each boolean is a byte of its own. */
static void thrift_skip_element(thrift_decoder_t* dec, thrift_type_t type) {
    if (type == THRIFT_TYPE_TRUE || type == THRIFT_TYPE_FALSE) {
        if (dec->status == CARQUET_OK &&
            carquet_buffer_reader_skip(&dec->reader, 1) != CARQUET_OK) {
            set_error(dec, CARQUET_ERROR_THRIFT_DECODE, "Truncated boolean element");
        }
        return;
    }
    thrift_skip(dec, type);
}

void thrift_skip(thrift_decoder_t* dec, thrift_type_t type) {
    if (dec->status != CARQUET_OK) {
        return;
    }

    switch (type) {
        case THRIFT_TYPE_STOP:
            /* STOP is a struct terminator, never a value type to skip.
             * Treating it as a no-op would cause infinite loops when it
             * appears as a container element type from malformed data. */
            set_error(dec, CARQUET_ERROR_THRIFT_DECODE, "Cannot skip STOP type");
            break;

        case THRIFT_TYPE_TRUE:
        case THRIFT_TYPE_FALSE:
            /* Boolean value is embedded in type, nothing to skip */
            dec->bool_pending = false;
            break;

        case THRIFT_TYPE_BYTE:
            carquet_buffer_reader_skip(&dec->reader, 1);
            break;

        case THRIFT_TYPE_I16:
        case THRIFT_TYPE_I32:
        case THRIFT_TYPE_I64:
            thrift_read_varint(dec);  /* Skip the varint */
            break;

        case THRIFT_TYPE_DOUBLE:
            carquet_buffer_reader_skip(&dec->reader, 8);
            break;

        case THRIFT_TYPE_BINARY: {
            int32_t len;
            thrift_read_binary(dec, &len);  /* Advances past the data */
            break;
        }

        case THRIFT_TYPE_LIST:
        case THRIFT_TYPE_SET: {
            thrift_type_t elem_type;
            int32_t count;
            thrift_read_list_begin(dec, &elem_type, &count);
            /* Bound recursion depth: a malformed footer can nest LIST/SET/MAP
             * containers arbitrarily (one C stack frame per level), which is not
             * covered by the STRUCT-only nesting cap. Charge each container
             * level against nesting_level so a hostile file cannot exhaust the
             * stack. read_struct_begin's own >= cap check keeps last_field_id[]
             * indexing in range even when the depth was spent on containers. */
            if (dec->nesting_level >= THRIFT_MAX_NESTING) {
                set_error(dec, CARQUET_ERROR_THRIFT_DECODE, "Thrift nesting too deep");
                break;
            }
            dec->nesting_level++;
            for (int32_t i = 0; i < count && dec->status == CARQUET_OK; i++) {
                thrift_skip_element(dec, elem_type);
            }
            dec->nesting_level--;
            break;
        }

        case THRIFT_TYPE_MAP: {
            thrift_type_t key_type, value_type;
            int32_t count;
            thrift_read_map_begin(dec, &key_type, &value_type, &count);
            if (dec->nesting_level >= THRIFT_MAX_NESTING) {
                set_error(dec, CARQUET_ERROR_THRIFT_DECODE, "Thrift nesting too deep");
                break;
            }
            dec->nesting_level++;
            for (int32_t i = 0; i < count && dec->status == CARQUET_OK; i++) {
                thrift_skip_element(dec, key_type);
                thrift_skip_element(dec, value_type);
            }
            dec->nesting_level--;
            break;
        }

        case THRIFT_TYPE_STRUCT: {
            thrift_read_struct_begin(dec);
            thrift_type_t field_type;
            int16_t field_id;
            while (thrift_read_field_begin(dec, &field_type, &field_id)) {
                thrift_skip(dec, field_type);
            }
            thrift_read_struct_end(dec);
            break;
        }

        case THRIFT_TYPE_UUID:
            carquet_buffer_reader_skip(&dec->reader, 16);
            break;

        default:
            set_error(dec, CARQUET_ERROR_THRIFT_INVALID_TYPE, "Unknown type to skip");
            break;
    }
}

/* ============================================================================
 * Utility Functions
 * ============================================================================
 */

const char* thrift_type_name(thrift_type_t type) {
    switch (type) {
        case THRIFT_TYPE_STOP: return "STOP";
        case THRIFT_TYPE_TRUE: return "TRUE";
        case THRIFT_TYPE_FALSE: return "FALSE";
        case THRIFT_TYPE_BYTE: return "BYTE";
        case THRIFT_TYPE_I16: return "I16";
        case THRIFT_TYPE_I32: return "I32";
        case THRIFT_TYPE_I64: return "I64";
        case THRIFT_TYPE_DOUBLE: return "DOUBLE";
        case THRIFT_TYPE_BINARY: return "BINARY";
        case THRIFT_TYPE_LIST: return "LIST";
        case THRIFT_TYPE_SET: return "SET";
        case THRIFT_TYPE_MAP: return "MAP";
        case THRIFT_TYPE_STRUCT: return "STRUCT";
        case THRIFT_TYPE_UUID: return "UUID";
        default: return "UNKNOWN";
    }
}
