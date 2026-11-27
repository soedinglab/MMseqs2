/**
 * @file thrift_encode.c
 * @brief Thrift Compact Protocol encoder implementation
 */

#include "thrift_encode.h"
#include "thrift_decode.h"  /* For type constants */
#include "core/endian.h"
#include <string.h>

/* ============================================================================
 * Internal Helpers
 * ============================================================================
 */

static void set_error(thrift_encoder_t* enc, carquet_status_t status) {
    if (enc->status == CARQUET_OK) {
        enc->status = status;
    }
}

/* ============================================================================
 * Encoder Lifecycle
 * ============================================================================
 */

void thrift_encoder_init(thrift_encoder_t* enc, carquet_buffer_t* buffer) {
    memset(enc, 0, sizeof(*enc));
    enc->buffer = buffer;
    enc->nesting_level = 0;
    enc->status = CARQUET_OK;
}

/* ============================================================================
 * Varint Writing
 * ============================================================================
 */

void thrift_write_varint(thrift_encoder_t* enc, uint64_t value) {
    uint8_t buf[10];
    int len = 0;

    while (value >= 0x80) {
        buf[len++] = (uint8_t)((value & 0x7F) | 0x80);
        value >>= 7;
    }
    buf[len++] = (uint8_t)value;

    if (carquet_buffer_append(enc->buffer, buf, len) != CARQUET_OK) {
        set_error(enc, CARQUET_ERROR_OUT_OF_MEMORY);
    }
}

void thrift_write_zigzag(thrift_encoder_t* enc, int64_t value) {
    uint64_t encoded = carquet_zigzag_encode64(value);
    thrift_write_varint(enc, encoded);
}

/* ============================================================================
 * Primitive Writing
 * ============================================================================
 */

void thrift_write_byte(thrift_encoder_t* enc, int8_t value) {
    if (carquet_buffer_append_byte(enc->buffer, (uint8_t)value) != CARQUET_OK) {
        set_error(enc, CARQUET_ERROR_OUT_OF_MEMORY);
    }
}

void thrift_write_i16(thrift_encoder_t* enc, int16_t value) {
    thrift_write_zigzag(enc, value);
}

void thrift_write_i32(thrift_encoder_t* enc, int32_t value) {
    thrift_write_zigzag(enc, value);
}

void thrift_write_i64(thrift_encoder_t* enc, int64_t value) {
    thrift_write_zigzag(enc, value);
}

void thrift_write_double(thrift_encoder_t* enc, double value) {
    if (carquet_buffer_append_f64_le(enc->buffer, value) != CARQUET_OK) {
        set_error(enc, CARQUET_ERROR_OUT_OF_MEMORY);
    }
}

void thrift_write_bool(thrift_encoder_t* enc, bool value) {
    /* When writing a standalone bool (not in a field header), use a byte */
    thrift_write_byte(enc, value ? 1 : 0);
}

void thrift_write_binary(thrift_encoder_t* enc, const uint8_t* data, int32_t length) {
    thrift_write_varint(enc, (uint64_t)length);

    if (length > 0 && data) {
        if (carquet_buffer_append(enc->buffer, data, (size_t)length) != CARQUET_OK) {
            set_error(enc, CARQUET_ERROR_OUT_OF_MEMORY);
        }
    }
}

void thrift_write_string(thrift_encoder_t* enc, const char* str) {
    if (!str) {
        thrift_write_binary(enc, NULL, 0);
        return;
    }
    thrift_write_binary(enc, (const uint8_t*)str, (int32_t)strlen(str));
}

void thrift_write_uuid(thrift_encoder_t* enc, const uint8_t uuid[16]) {
    if (carquet_buffer_append(enc->buffer, uuid, 16) != CARQUET_OK) {
        set_error(enc, CARQUET_ERROR_OUT_OF_MEMORY);
    }
}

/* ============================================================================
 * Struct Writing
 * ============================================================================
 */

void thrift_write_struct_begin(thrift_encoder_t* enc) {
    if (enc->nesting_level >= THRIFT_ENCODER_MAX_NESTING) {
        set_error(enc, CARQUET_ERROR_THRIFT_ENCODE);
        return;
    }

    enc->last_field_id[enc->nesting_level] = 0;
    enc->nesting_level++;
}

void thrift_write_struct_end(thrift_encoder_t* enc) {
    thrift_write_field_stop(enc);

    if (enc->nesting_level > 0) {
        enc->nesting_level--;
    }
}

void thrift_write_field_header(thrift_encoder_t* enc, int type, int16_t field_id) {
    int16_t last_id = 0;
    if (enc->nesting_level > 0) {
        last_id = enc->last_field_id[enc->nesting_level - 1];
    }

    int16_t delta = field_id - last_id;

    if (delta > 0 && delta <= 15) {
        /* Use compact form: delta in upper nibble, type in lower */
        uint8_t header = (uint8_t)(((delta & 0x0F) << 4) | (type & 0x0F));
        thrift_write_byte(enc, (int8_t)header);
    } else {
        /* Use extended form: type byte followed by field ID */
        thrift_write_byte(enc, (int8_t)(type & 0x0F));
        thrift_write_i16(enc, field_id);
    }

    /* Update last field ID */
    if (enc->nesting_level > 0) {
        enc->last_field_id[enc->nesting_level - 1] = field_id;
    }
}

void thrift_write_field_stop(thrift_encoder_t* enc) {
    thrift_write_byte(enc, 0);
}

/* ============================================================================
 * Container Writing
 * ============================================================================
 */

void thrift_write_list_begin(thrift_encoder_t* enc, int elem_type, int32_t count) {
    if (count < 15) {
        /* Compact form: count in upper nibble */
        uint8_t header = (uint8_t)(((count & 0x0F) << 4) | (elem_type & 0x0F));
        thrift_write_byte(enc, (int8_t)header);
    } else {
        /* Extended form: 0xF in upper nibble, followed by varint count */
        uint8_t header = (uint8_t)((0x0F << 4) | (elem_type & 0x0F));
        thrift_write_byte(enc, (int8_t)header);
        thrift_write_varint(enc, (uint64_t)count);
    }
}

void thrift_write_set_begin(thrift_encoder_t* enc, int elem_type, int32_t count) {
    /* Set has the same encoding as list */
    thrift_write_list_begin(enc, elem_type, count);
}

void thrift_write_map_begin(thrift_encoder_t* enc,
                             int key_type, int value_type, int32_t count) {
    if (count == 0) {
        thrift_write_byte(enc, 0);
        return;
    }

    thrift_write_varint(enc, (uint64_t)count);

    uint8_t types = (uint8_t)(((key_type & 0x0F) << 4) | (value_type & 0x0F));
    thrift_write_byte(enc, (int8_t)types);
}
