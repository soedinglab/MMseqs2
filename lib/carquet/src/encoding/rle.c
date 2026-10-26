/**
 * @file rle.c
 * @brief RLE/Bit-packing hybrid encoding implementation
 */

#include "rle.h"
#include "core/endian.h"
#include "core/bitpack.h"
#include <string.h>

/* SIMD dispatch for optimized level fill (AVX2/AVX-512/NEON/SVE) */
extern void carquet_dispatch_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value);

/* ============================================================================
 * Internal Helpers
 * ============================================================================
 */

static int read_varint(const uint8_t* data, size_t size, size_t* pos, uint32_t* out) {
    uint32_t result = 0;
    int shift = 0;
    size_t p = *pos;

    while (p < size && shift < 32) {
        uint8_t byte = data[p++];
        result |= (uint32_t)(byte & 0x7F) << shift;
        if ((byte & 0x80) == 0) {
            *pos = p;
            *out = result;
            return 0;
        }
        shift += 7;
    }

    return -1;  /* Truncated or overflow */
}

static bool start_new_run(carquet_rle_decoder_t* dec) {
    /* Loop rather than tail-recurse over empty runs: a hostile page can be a
     * long run of zero-length-run headers, which would otherwise recurse one
     * C stack frame per header (only elided under optimization/TCO). */
    for (;;) {
        if (dec->pos >= dec->size) {
            return false;
        }

        /* Read header */
        uint32_t header;
        if (read_varint(dec->data, dec->size, &dec->pos, &header) < 0) {
            dec->status = CARQUET_ERROR_INVALID_RLE;
            return false;
        }

        if ((header & 1) == 0) {
            /* RLE run */
            dec->in_rle_run = true;
            dec->run_remaining = (int64_t)(header >> 1);

            if (dec->run_remaining == 0) {
                /* Empty run, try next */
                continue;
            }

            /* Read the repeated value (ceil(bit_width/8) bytes) */
            int value_bytes = (dec->bit_width + 7) / 8;
            if (dec->pos + (size_t)value_bytes > dec->size) {
                dec->status = CARQUET_ERROR_INVALID_RLE;
                return false;
            }

            dec->rle_value = 0;
            for (int i = 0; i < value_bytes; i++) {
                dec->rle_value |= (uint32_t)dec->data[dec->pos++] << (i * 8);
            }
            dec->rle_value &= dec->value_mask;

        } else {
            /* Bit-packed run */
            dec->in_rle_run = false;
            int num_groups = (int)(header >> 1);  /* Number of 8-value groups */
            dec->run_remaining = (int64_t)num_groups * 8;

            if (dec->run_remaining == 0) {
                continue;
            }

            /* We'll decode 8 values at a time into the buffer */
            dec->bitpack_pos = 0;
            dec->bitpack_count = 0;
        }

        return true;
    }
}

static bool fill_bitpack_buffer(carquet_rle_decoder_t* dec) {
    if (dec->run_remaining <= 0) {
        return false;
    }

    /* Read 8 packed values */
    size_t bytes_needed = (size_t)dec->bit_width;  /* 8 values * bit_width bits = bit_width bytes */
    if (dec->pos + bytes_needed > dec->size) {
        dec->status = CARQUET_ERROR_INVALID_RLE;
        return false;
    }

    carquet_bitunpack8_32(dec->data + dec->pos, dec->bit_width, dec->bitpack_buffer);
    dec->pos += bytes_needed;
    dec->bitpack_pos = 0;
    dec->bitpack_count = 8;

    return true;
}

/* ============================================================================
 * RLE Decoder
 * ============================================================================
 */

void carquet_rle_decoder_init(
    carquet_rle_decoder_t* dec,
    const uint8_t* data,
    size_t size,
    int bit_width) {

    memset(dec, 0, sizeof(*dec));
    dec->data = data;
    dec->size = size;
    /* RLE values are uint32_t — bit_width must be 0..32 */
    if (bit_width < 0 || bit_width > 32) {
        dec->bit_width = 0;
        dec->value_mask = 0;
        dec->status = CARQUET_ERROR_INVALID_RLE;
        return;
    }
    dec->bit_width = bit_width;
    dec->value_mask = bit_width >= 32 ? ~0U : (1U << bit_width) - 1;
    dec->status = CARQUET_OK;
}

bool carquet_rle_decoder_has_next(const carquet_rle_decoder_t* dec) {
    if (dec->status != CARQUET_OK) {
        return false;
    }
    if (dec->run_remaining > 0) {
        return true;
    }
    return dec->pos < dec->size;
}

uint32_t carquet_rle_decoder_get(carquet_rle_decoder_t* dec) {
    if (dec->status != CARQUET_OK) {
        return 0;
    }

    /* Need new run? */
    if (dec->run_remaining <= 0) {
        if (!start_new_run(dec)) {
            return 0;
        }
    }

    if (dec->in_rle_run) {
        dec->run_remaining--;
        return dec->rle_value;
    } else {
        /* Bit-packed run */
        if (dec->bitpack_pos >= dec->bitpack_count) {
            if (!fill_bitpack_buffer(dec)) {
                return 0;
            }
        }

        dec->run_remaining--;
        return dec->bitpack_buffer[dec->bitpack_pos++];
    }
}

int64_t carquet_rle_decoder_get_batch(
    carquet_rle_decoder_t* dec,
    uint32_t* output,
    int64_t count) {

    int64_t read = 0;

    while (read < count && carquet_rle_decoder_has_next(dec)) {
        /* Need new run? */
        if (dec->run_remaining <= 0) {
            if (!start_new_run(dec)) {
                break;
            }
        }

        if (dec->in_rle_run) {
            /* Fill with repeated value */
            int64_t to_fill = count - read;
            if (to_fill > dec->run_remaining) {
                to_fill = dec->run_remaining;
            }

            for (int64_t i = 0; i < to_fill; i++) {
                output[read++] = dec->rle_value;
            }
            dec->run_remaining -= to_fill;

        } else {
            /* Bit-packed run */
            while (read < count && dec->run_remaining > 0) {
                if (dec->bitpack_pos >= dec->bitpack_count) {
                    if (!fill_bitpack_buffer(dec)) {
                        break;
                    }
                }

                /* Copy from buffer */
                while (read < count && dec->bitpack_pos < dec->bitpack_count &&
                       dec->run_remaining > 0) {
                    output[read++] = dec->bitpack_buffer[dec->bitpack_pos++];
                    dec->run_remaining--;
                }
            }
        }
    }

    return read;
}

int64_t carquet_rle_decoder_skip(
    carquet_rle_decoder_t* dec,
    int64_t count) {

    int64_t skipped = 0;

    while (skipped < count && carquet_rle_decoder_has_next(dec)) {
        if (dec->run_remaining <= 0) {
            if (!start_new_run(dec)) {
                break;
            }
        }

        int64_t to_skip = count - skipped;
        if (to_skip > dec->run_remaining) {
            to_skip = dec->run_remaining;
        }

        if (dec->in_rle_run) {
            /* Easy - just reduce count */
            skipped += to_skip;
            dec->run_remaining -= to_skip;
        } else {
            /* Need to actually advance through bit-packed data */
            while (to_skip > 0 && dec->run_remaining > 0) {
                if (dec->bitpack_pos >= dec->bitpack_count) {
                    if (!fill_bitpack_buffer(dec)) {
                        break;
                    }
                }

                int64_t can_skip = dec->bitpack_count - dec->bitpack_pos;
                if (can_skip > to_skip) can_skip = to_skip;
                if (can_skip > dec->run_remaining) can_skip = dec->run_remaining;

                dec->bitpack_pos += (int)can_skip;
                dec->run_remaining -= can_skip;
                skipped += can_skip;
                to_skip -= can_skip;
            }
        }
    }

    return skipped;
}

/* ============================================================================
 * RLE Encoder
 * ============================================================================
 */

static void write_varint(carquet_buffer_t* buf, uint32_t value) {
    uint8_t bytes[5];
    int len = 0;

    while (value >= 0x80) {
        bytes[len++] = (uint8_t)((value & 0x7F) | 0x80);
        value >>= 7;
    }
    bytes[len++] = (uint8_t)value;

    carquet_buffer_append(buf, bytes, (size_t)len);
}

static inline int rle_value_bytes(int bit_width) {
    return (bit_width + 7) / 8;
}

static inline void write_value_bytes(uint8_t* dst, uint32_t value, int value_bytes) {
    /* Unrolled for common cases (1, 2, 4 bytes) */
    dst[0] = (uint8_t)value;
    if (value_bytes > 1) dst[1] = (uint8_t)(value >> 8);
    if (value_bytes > 2) dst[2] = (uint8_t)(value >> 16);
    if (value_bytes > 3) dst[3] = (uint8_t)(value >> 24);
}

static void flush_rle(carquet_rle_encoder_t* enc) {
    /* The run header is `count << 1` in a 32-bit varint, so one RLE run holds
     * at most 2^31 - 1 values; a longer run (int64 repeat_count, reachable via
     * put_repeat) is emitted as several back-to-back runs. */
    while (enc->repeat_count > 0) {
        int64_t run = enc->repeat_count;
        if (run > 0x7FFFFFFF) run = 0x7FFFFFFF;

        /* Write RLE header + value in a single buffer append where possible */
        int vb = rle_value_bytes(enc->bit_width);
        uint8_t buf[9]; /* max 5 (varint) + 4 (value) */
        int len = 0;

        /* Inline varint encoding */
        uint32_t header = (uint32_t)run << 1;
        while (header >= 0x80) {
            buf[len++] = (uint8_t)((header & 0x7F) | 0x80);
            header >>= 7;
        }
        buf[len++] = (uint8_t)header;

        /* Append value bytes */
        write_value_bytes(buf + len, enc->prev_value, vb);
        len += vb;

        carquet_buffer_append(enc->buffer, buf, (size_t)len);
        enc->repeat_count -= run;
    }
}

/* Emit the staged complete groups as one bit-packed run. */
static void emit_staged_groups(carquet_rle_encoder_t* enc) {
    if (enc->staged_groups == 0) return;
    write_varint(enc->buffer, (uint32_t)((enc->staged_groups << 1) | 1));
    carquet_buffer_append(enc->buffer, enc->staged,
                          (size_t)enc->staged_groups * (size_t)enc->bit_width);
    enc->staged_groups = 0;
}

/* Pack the 8 buffered values into the staging area; emit when full. */
static void stage_full_group(carquet_rle_encoder_t* enc) {
    carquet_bitpack8_32(enc->bitpack_buffer, enc->bit_width,
                        enc->staged + (size_t)enc->staged_groups * (size_t)enc->bit_width);
    enc->staged_groups++;
    enc->bitpack_count = 0;
    if (enc->staged_groups == 63) {
        emit_staged_groups(enc);
    }
}

static void flush_bitpack_as_rle(carquet_rle_encoder_t* enc) {
    /* Emit remaining buffered values as individual RLE runs.
     * Used when we have a partial group (< 8 values) in the middle of the
     * stream (an RLE run follows): a bit-packed group always decodes to 8
     * values, so it cannot be padded there. */
    int vb = rle_value_bytes(enc->bit_width);

    int i = 0;
    while (i < enc->bitpack_count) {
        uint32_t val = enc->bitpack_buffer[i];
        int64_t run = 1;
        while (i + run < enc->bitpack_count && enc->bitpack_buffer[i + run] == val) {
            run++;
        }
        /* Write RLE header + value in single append */
        uint8_t buf[9];
        int len = 0;
        uint32_t header = (uint32_t)(run << 1);
        while (header >= 0x80) {
            buf[len++] = (uint8_t)((header & 0x7F) | 0x80);
            header >>= 7;
        }
        buf[len++] = (uint8_t)header;
        write_value_bytes(buf + len, val, vb);
        len += vb;
        carquet_buffer_append(enc->buffer, buf, (size_t)len);
        i += (int)run;
    }

    enc->bitpack_count = 0;
    enc->bitpack_total = 0;
}

/* Flush the bit-packed state: staged complete groups as one run, then any
 * partial group as RLE runs. The partial group is never zero-padded, even at
 * the end of the stream: the page writer concatenates one encoded stream per
 * add_values() call into a single levels block, so every stream must decode
 * to exactly its own value count. */
static void flush_bitpack(carquet_rle_encoder_t* enc) {
    if (enc->bitpack_count == 8) {
        stage_full_group(enc);
    }
    emit_staged_groups(enc);
    if (enc->bitpack_count > 0) {
        flush_bitpack_as_rle(enc);
    }
    enc->bitpack_count = 0;
    enc->bitpack_total = 0;
}

void carquet_rle_encoder_init(
    carquet_rle_encoder_t* enc,
    carquet_buffer_t* buffer,
    int bit_width) {

    memset(enc, 0, sizeof(*enc));
    enc->buffer = buffer;
    enc->bit_width = bit_width;
    enc->status = CARQUET_OK;
}

carquet_status_t carquet_rle_encoder_put(
    carquet_rle_encoder_t* enc,
    uint32_t value) {

    if (enc->status != CARQUET_OK) {
        return enc->status;
    }

    if (!enc->has_prev) {
        enc->prev_value = value;
        enc->repeat_count = 1;
        enc->has_prev = true;
        return CARQUET_OK;
    }

    if (value == enc->prev_value) {
        enc->repeat_count++;
        return CARQUET_OK;
    }

    /* Value changed */
    if (enc->repeat_count >= 8) {
        /* Flush as RLE */
        flush_bitpack(enc);  /* Flush any pending bit-pack */
        flush_rle(enc);
    } else {
        /* Add to bit-pack buffer */
        for (int64_t i = 0; i < enc->repeat_count; i++) {
            enc->bitpack_buffer[enc->bitpack_count++] = enc->prev_value;
            enc->bitpack_total++;

            if (enc->bitpack_count == 8) {
                stage_full_group(enc);
            }
        }
        enc->repeat_count = 0;
    }

    enc->prev_value = value;
    enc->repeat_count = 1;
    return CARQUET_OK;
}

carquet_status_t carquet_rle_encoder_put_repeat(
    carquet_rle_encoder_t* enc,
    uint32_t value,
    int64_t count) {

    if (count <= 0 || enc->status != CARQUET_OK) {
        return enc->status;
    }

    /* If same value as current run, just extend the count */
    if (enc->has_prev && value == enc->prev_value) {
        enc->repeat_count += count;
        return CARQUET_OK;
    }

    /* First value or value changed: flush then set up new run */
    if (enc->has_prev) {
        carquet_status_t status = carquet_rle_encoder_put(enc, value);
        if (status != CARQUET_OK) return status;
        /* put() set repeat_count=1 for the new value, add remaining */
        enc->repeat_count += count - 1;
    } else {
        enc->prev_value = value;
        enc->repeat_count = count;
        enc->has_prev = true;
    }
    return CARQUET_OK;
}

carquet_status_t carquet_rle_encoder_flush(carquet_rle_encoder_t* enc) {
    if (enc->status != CARQUET_OK) {
        return enc->status;
    }

    if (enc->repeat_count >= 8) {
        flush_bitpack(enc);
        flush_rle(enc);
    } else if (enc->repeat_count > 0) {
        for (int64_t i = 0; i < enc->repeat_count; i++) {
            enc->bitpack_buffer[enc->bitpack_count++] = enc->prev_value;
            enc->bitpack_total++;

            if (enc->bitpack_count == 8) {
                stage_full_group(enc);
            }
        }
        enc->repeat_count = 0;
    }
    /* End of stream: emit staged groups, then the partial tail as RLE runs. */
    flush_bitpack(enc);

    return CARQUET_OK;
}

/* ============================================================================
 * Convenience Functions
 * ============================================================================
 */

int64_t carquet_rle_decode_all(
    const uint8_t* input,
    size_t input_size,
    int bit_width,
    uint32_t* output,
    int64_t max_values) {

    carquet_rle_decoder_t dec;
    carquet_rle_decoder_init(&dec, input, input_size, bit_width);
    return carquet_rle_decoder_get_batch(&dec, output, max_values);
}

int64_t carquet_rle_decode_levels(
    const uint8_t* input,
    size_t input_size,
    int bit_width,
    int16_t* output,
    int64_t max_values) {

    if (max_values <= 0 || input_size == 0 || bit_width < 0 || bit_width > 32) {
        return 0;
    }

    /* Fast path: decode directly without per-value function calls */
    size_t pos = 0;
    int64_t count = 0;
    uint32_t value_mask = bit_width >= 32 ? ~0U : (1U << bit_width) - 1;
    int value_bytes = (bit_width + 7) / 8;

    while (count < max_values && pos < input_size) {
        /* Read varint header inline */
        uint32_t header = 0;
        int shift = 0;
        while (pos < input_size && shift < 32) {
            uint8_t byte = input[pos++];
            header |= (uint32_t)(byte & 0x7F) << shift;
            if ((byte & 0x80) == 0) break;
            shift += 7;
        }

        if ((header & 1) == 0) {
            /* RLE run: fill output with repeated value */
            int64_t run_length = (int64_t)(header >> 1);
            if (run_length == 0) continue;

            if (pos + (size_t)value_bytes > input_size) break;

            /* Read the repeated value */
            uint32_t rle_value = 0;
            for (int i = 0; i < value_bytes; i++) {
                rle_value |= (uint32_t)input[pos++] << (i * 8);
            }
            rle_value &= value_mask;
            int16_t val16 = (int16_t)rle_value;

            /* Fill output in bulk */
            int64_t to_fill = run_length;
            if (count + to_fill > max_values) {
                to_fill = max_values - count;
            }

            /* Use SIMD-dispatched fill (AVX2/AVX-512/NEON/SVE at runtime) */
            carquet_dispatch_fill_def_levels(output + count, to_fill, val16);
            count += to_fill;

        } else {
            /* Bit-packed run: decode 8 values at a time */
            int num_groups = (int)(header >> 1);
            int64_t run_length = (int64_t)num_groups * 8;
            if (run_length == 0) continue;

            size_t bytes_per_group = (size_t)bit_width;

            /* Fast path: the whole bit-packed run is present in the input
             * and fits in the output. Bulk-decode through the wide-SIMD
             * unpacker (carquet_bitunpack_32) instead of one group at a
             * time. The division form of the bound avoids size_t overflow.
             * Truncated input or an output cap mid-run (rare) falls through
             * to the exact per-8 loop below, whose semantics are unchanged. */
            if (bit_width >= 1 &&
                (size_t)num_groups <= (input_size - pos) / bytes_per_group &&
                count + run_length <= max_values) {
                int64_t done = 0;
                uint32_t tmp[256];  /* 256 is a multiple of 8 */
                while (done < run_length) {
                    int64_t chunk = run_length - done;
                    if (chunk > 256) chunk = 256;
                    size_t used = carquet_bitunpack_32(input + pos,
                                        (size_t)chunk, bit_width, tmp);
                    pos += used;
                    for (int64_t k = 0; k < chunk; k++) {
                        output[count++] = (int16_t)tmp[k];
                    }
                    done += chunk;
                }
                continue;
            }

            for (int g = 0; g < num_groups && count < max_values; g++) {
                if (pos + bytes_per_group > input_size) break;

                /* Unpack 8 values */
                uint32_t temp[8];
                carquet_bitunpack8_32(input + pos, bit_width, temp);
                pos += bytes_per_group;

                /* Convert to int16_t and store */
                int64_t to_store = 8;
                if (count + to_store > max_values) {
                    to_store = max_values - count;
                }

                for (int64_t i = 0; i < to_store; i++) {
                    output[count++] = (int16_t)temp[i];
                }
            }
        }
    }

    return count;
}

int64_t carquet_rle_decode_to_bitmap(
    const uint8_t* input,
    size_t input_size,
    uint8_t* bitmap,
    int64_t max_values,
    int64_t* non_null_count) {

    if (max_values <= 0 || input_size == 0) {
        if (non_null_count) *non_null_count = 0;
        return 0;
    }

    /* Pre-clear the bitmap */
    size_t bitmap_bytes = ((size_t)max_values + 7) / 8;
    memset(bitmap, 0, bitmap_bytes);

    size_t pos = 0;
    int64_t count = 0;
    int64_t nn_count = 0;

    /* For max_def_level == 1: bit_width is 1, value is 1 byte */
    const int value_bytes = 1;

    while (count < max_values && pos < input_size) {
        /* Read varint header inline */
        uint32_t header = 0;
        int shift = 0;
        while (pos < input_size && shift < 32) {
            uint8_t byte = input[pos++];
            header |= (uint32_t)(byte & 0x7F) << shift;
            if ((byte & 0x80) == 0) break;
            shift += 7;
        }

        if ((header & 1) == 0) {
            /* RLE run */
            int64_t run_length = (int64_t)(header >> 1);
            if (run_length == 0) continue;

            if (pos + (size_t)value_bytes > input_size) break;

            uint32_t rle_value = input[pos++] & 1;

            int64_t to_fill = run_length;
            if (count + to_fill > max_values) {
                to_fill = max_values - count;
            }

            if (rle_value == 1) {
                /* Present: set bits in bitmap (convention: bit set = present) */
                for (int64_t i = 0; i < to_fill; i++) {
                    int64_t idx = count + i;
                    bitmap[idx / 8] |= (uint8_t)(1 << (idx % 8));
                }
                nn_count += to_fill;
            } else {
                /* Null: bits stay clear (already memset to 0) */
            }

            count += to_fill;

        } else {
            /* Bit-packed run: 8 values per group, 1 bit each = 1 byte per group */
            int num_groups = (int)(header >> 1);
            if (num_groups == 0) continue;

            for (int g = 0; g < num_groups && count < max_values; g++) {
                if (pos >= input_size) break;

                uint8_t packed_byte = input[pos++];

                /* Each bit in the byte is one def_level value (0 or 1) */
                int64_t to_store = 8;
                if (count + to_store > max_values) {
                    to_store = max_values - count;
                }

                for (int64_t i = 0; i < to_store; i++) {
                    uint8_t bit = (packed_byte >> i) & 1;
                    if (bit == 1) {
                        /* Present: set bit (convention: bit set = present) */
                        bitmap[(count + i) / 8] |= (uint8_t)(1 << ((count + i) % 8));
                        nn_count++;
                    }
                }
                count += to_store;
            }
        }
    }

    if (non_null_count) *non_null_count = nn_count;
    return count;
}

carquet_status_t carquet_rle_encode_all(
    const uint32_t* input,
    int64_t count,
    int bit_width,
    carquet_buffer_t* output) {

    carquet_rle_encoder_t enc;
    carquet_rle_encoder_init(&enc, output, bit_width);

    for (int64_t i = 0; i < count; i++) {
        carquet_status_t status = carquet_rle_encoder_put(&enc, input[i]);
        if (status != CARQUET_OK) return status;
    }

    return carquet_rle_encoder_flush(&enc);
}

carquet_status_t carquet_rle_encode_levels(
    const int16_t* input,
    int64_t count,
    int bit_width,
    carquet_buffer_t* output) {

    carquet_rle_encoder_t enc;
    carquet_rle_encoder_init(&enc, output, bit_width);

    for (int64_t i = 0; i < count; i++) {
        carquet_status_t status = carquet_rle_encoder_put(&enc, (uint32_t)input[i]);
        if (status != CARQUET_OK) return status;
    }

    return carquet_rle_encoder_flush(&enc);
}

int64_t carquet_rle_decode_levels_prefixed(
    const uint8_t* input,
    size_t input_size,
    int bit_width,
    int16_t* output,
    int64_t max_values,
    size_t* bytes_consumed) {

    if (input_size < 4) {
        if (bytes_consumed) *bytes_consumed = 0;
        return -1;
    }

    /* Read 4-byte length prefix (little-endian). Widen to 64-bit before adding
     * the 4-byte prefix: `4 + rle_length` in 32-bit wraps when rle_length is
     * near UINT32_MAX (e.g. 0xFFFFFFFF -> 3), which would pass this check and
     * hand ~4 GB to the decoder over a tiny buffer. */
    uint32_t rle_length = carquet_read_u32_le(input);
    if ((uint64_t)rle_length + 4 > (uint64_t)input_size) {
        if (bytes_consumed) *bytes_consumed = 0;
        return -1;
    }

    int64_t count = carquet_rle_decode_levels(
        input + 4, rle_length, bit_width, output, max_values);

    if (bytes_consumed) {
        *bytes_consumed = 4 + rle_length;
    }

    return count;
}
