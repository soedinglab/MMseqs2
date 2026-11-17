/**
 * @file rle.c
 * @brief RLE/Bit-packing hybrid encoding implementation
 */

#include "rle.h"
#include "core/endian.h"
#include "core/bitpack.h"
#include "core/allocator.h"
#include <stdint.h>
#include <string.h>

/* SIMD dispatch for optimized level fill (AVX2/AVX-512/NEON/SVE) */
extern void carquet_dispatch_fill_def_levels(int16_t* def_levels, int64_t count, int16_t value);
/* SIMD dispatch: pack 0/1 levels into a bit array (1-bit level encoder) */
extern bool carquet_dispatch_levels_to_bits(const int16_t* levels, int64_t count, uint64_t* bits);

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
            /* Bit-packed run. Bulk fast path: while the 8-value staging
             * buffer is empty, unpack whole groups straight into the output
             * through carquet_bitunpack_32 (wide SIMD kernels, one dispatch
             * lookup per call) instead of one dispatched group + copy loop
             * per 8 values. Only groups fully present in the input are taken;
             * a truncated tail falls through to the per-group loop below,
             * which reports the error exactly as before. */
            if (dec->bitpack_pos >= dec->bitpack_count && dec->bit_width > 0 &&
                dec->pos < dec->size) {
                int64_t want = count - read;
                if (want > dec->run_remaining) {
                    want = dec->run_remaining;
                }
                size_t groups = (size_t)(want / 8);
                size_t avail = (dec->size - dec->pos) / (size_t)dec->bit_width;
                if (groups > avail) {
                    groups = avail;
                }
                if (groups > 0) {
                    carquet_bitunpack_32(dec->data + dec->pos, groups * 8,
                                         dec->bit_width, output + read);
                    dec->pos += groups * (size_t)dec->bit_width;
                    read += (int64_t)(groups * 8);
                    dec->run_remaining -= (int64_t)(groups * 8);
                }
            }
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

/* Every encoder write goes through here. A failed append (out of memory)
 * is recorded in enc->status, later writes are dropped, and put / flush
 * return the error instead of CARQUET_OK with a truncated stream. */
static void enc_append(carquet_rle_encoder_t* enc, const void* data, size_t len) {
    if (enc->status != CARQUET_OK) return;
    carquet_status_t st = carquet_buffer_append(enc->buffer, data, len);
    if (st != CARQUET_OK) enc->status = st;
}

static void write_varint(carquet_rle_encoder_t* enc, uint32_t value) {
    uint8_t bytes[5];
    int len = 0;

    while (value >= 0x80) {
        bytes[len++] = (uint8_t)((value & 0x7F) | 0x80);
        value >>= 7;
    }
    bytes[len++] = (uint8_t)value;

    enc_append(enc, bytes, (size_t)len);
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

        enc_append(enc, buf, (size_t)len);
        enc->repeat_count -= run;
    }
}

/* Emit the staged complete groups as one bit-packed run. */
static void emit_staged_groups(carquet_rle_encoder_t* enc) {
    if (enc->staged_groups == 0) return;
    write_varint(enc, (uint32_t)((enc->staged_groups << 1) | 1));
    enc_append(enc, enc->staged,
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
        enc_append(enc, buf, (size_t)len);
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

/* The encoder's step for one value; inlined into the bulk encoders below,
 * which call it for every value that does not simply extend the current
 * run. */
static inline carquet_status_t rle_put(
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
    return enc->status;
}

carquet_status_t carquet_rle_encoder_put(
    carquet_rle_encoder_t* enc,
    uint32_t value) {
    return rle_put(enc, value);
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

    return enc->status;
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

/* 1-bit levels of `nbytes` packed bytes (LSB first) as int16: out[8i + j] =
 * bit j of src[i]. One multiply spreads a nibble x over four int16 lanes:
 * x * (1 + 2^15 + 2^30 + 2^45) holds copies of x at bits 0, 15, 30 and 45
 * (disjoint, so no carries), which puts bit k of x at bit 16k. */
static void expand_bits_to_levels(const uint8_t* src, size_t nbytes, int16_t* out) {
#if CARQUET_LITTLE_ENDIAN
    const uint64_t spread = 0x0000200040008001ULL, lanes = 0x0001000100010001ULL;
    for (size_t i = 0; i < nbytes; i++) {
        uint64_t lo = ((uint64_t)(src[i] & 0x0Fu) * spread) & lanes;
        uint64_t hi = ((uint64_t)(src[i] >> 4) * spread) & lanes;
        memcpy(out + 8 * i, &lo, sizeof(lo));
        memcpy(out + 8 * i + 4, &hi, sizeof(hi));
    }
#else
    for (size_t i = 0; i < nbytes; i++) {
        for (int j = 0; j < 8; j++) out[8 * i + (size_t)j] = (int16_t)((src[i] >> j) & 1);
    }
#endif
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

            /* Short runs (the common case for scattered nulls) are cheaper
             * inline than through the dispatched SIMD fill. */
            if (to_fill <= 16) {
                int16_t* dst = output + count;
                for (int64_t k = 0; k < to_fill; k++) dst[k] = val16;
            } else {
                carquet_dispatch_fill_def_levels(output + count, to_fill, val16);
            }
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
                if (bit_width == 1) {
                    /* One byte per group: expand the bits directly. */
                    expand_bits_to_levels(input + pos, (size_t)num_groups, output + count);
                    pos += (size_t)num_groups;
                    count += run_length;
                    continue;
                }
                int64_t done = 0;
                uint32_t tmp[256];  /* 256 is a multiple of 8 */
                while (done < run_length) {
                    int64_t chunk = run_length - done;
                    if (chunk > 256) chunk = 256;
                    size_t used = carquet_bitunpack_32(input + pos,
                                        (size_t)chunk, bit_width, tmp);
                    pos += used;
                    /* Base-pointer form: MSVC refuses to vectorize the
                     * output[count++] shape (C5002 reason 1203). */
                    int16_t* dst = output + count;
                    for (int64_t k = 0; k < chunk; k++) {
                        dst[k] = (int16_t)tmp[k];
                    }
                    count += chunk;
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

/* Set bits [at, at + n) of `bm` (n > 0). Touches only the bytes covering
 * the range. */
static void bitmap_set_range(uint8_t* bm, int64_t at, int64_t n) {
    int64_t end = at + n;
    int64_t first = at >> 3;
    int64_t last = (end - 1) >> 3;
    uint8_t head = (uint8_t)(0xFFu << (at & 7));
    uint8_t tail = (uint8_t)(0xFFu >> (7 - ((end - 1) & 7)));
    if (first == last) {
        bm[first] |= (uint8_t)(head & tail);
        return;
    }
    bm[first] |= head;
    if (last - first > 1) {
        memset(bm + first + 1, 0xFF, (size_t)(last - first - 1));
    }
    bm[last] |= tail;
}

/* Append the first `n` bits of `src` (LSB first) to `bm` at bit `at`, whose
 * bits from `at` on are still zero. Reads ceil(n/8) bytes of `src` and
 * writes only the bytes of `bm` covering [at, at + n). */
static void bitmap_append(uint8_t* bm, int64_t at, const uint8_t* src, int64_t n) {
    uint8_t* d = bm + (at >> 3);
    unsigned sh = (unsigned)(at & 7);
    if (sh == 0) {
        size_t full = (size_t)(n >> 3);
        memcpy(d, src, full);
        if (n & 7) d[full] = (uint8_t)(src[full] & ((1u << (n & 7)) - 1));
        return;
    }
    int64_t i = 0;
    /* d[0] keeps its low sh bits; a full word's top sh bits land in d[8],
     * which is inside the range because all 64 source bits are. */
    for (; n - i >= 64; i += 64, src += 8, d += 8) {
        uint64_t w = carquet_read_u64_le(src);
        carquet_write_u64_le(d, carquet_read_u64_le(d) | (w << sh));
        d[8] = (uint8_t)(w >> (64 - sh));
    }
    for (; n - i >= 8; i += 8, src++, d++) {
        unsigned b = *src;
        d[0] |= (uint8_t)(b << sh);
        d[1] = (uint8_t)(b >> (8 - sh));
    }
    if (n - i > 0) {
        unsigned r = (unsigned)(n - i);
        unsigned b = *src & ((1u << r) - 1);
        d[0] |= (uint8_t)(b << sh);
        if (sh + r > 8) d[1] = (uint8_t)(b >> (8 - sh));
    }
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

    /* Pre-clear the bitmap and its 8 slack bytes: every run is then OR-ed
     * in, and a short one as a single 64-bit read-modify-write at the byte
     * holding its first bit (which the slack keeps inside the buffer). */
    size_t bitmap_bytes = ((size_t)max_values + 7) / 8;
    memset(bitmap, 0, bitmap_bytes + 8);

    size_t pos = 0;
    int64_t count = 0;

    /* Same stream walk as carquet_rle_decode_levels at bit_width 1 (so the
     * two agree on every input, truncated ones included): a bit-packed
     * group of 8 one-bit values is one byte whose bits, LSB first, are the
     * levels, so it is the bitmap byte itself; an RLE run is a bit range. */
    while (count < max_values && pos < input_size) {
        uint32_t header = input[pos++];
        if (header & 0x80) {
            header &= 0x7F;
            int shift = 7;
            while (pos < input_size && shift < 32) {
                uint8_t byte = input[pos++];
                header |= (uint32_t)(byte & 0x7F) << shift;
                if ((byte & 0x80) == 0) break;
                shift += 7;
            }
        }
        int64_t room = max_values - count;
        uint8_t* d = bitmap + (count >> 3);
        unsigned sh = (unsigned)(count & 7);

        if ((header & 1) == 0) {
            int64_t run = (int64_t)(header >> 1);
            if (run == 0) continue;
            if (pos >= input_size) break;
            uint64_t value = input[pos++] & 1u;
            if (run > room) run = room;
            if (run + (int64_t)sh <= 64) {
                uint64_t m = (run == 64 ? ~(uint64_t)0 : (((uint64_t)1 << run) - 1)) << sh;
                carquet_write_u64_le(d, carquet_read_u64_le(d) | (m & ((uint64_t)0 - value)));
            } else if (value) {
                bitmap_set_range(bitmap, count, run);
            }
            count += run;
        } else {
            size_t num_groups = (size_t)(header >> 1);
            if (num_groups == 0) continue;
            /* Whole groups only: a truncated tail stops the decode, as the
             * per-group loop of carquet_rle_decode_levels does. */
            size_t avail = input_size - pos;
            if (num_groups > avail) num_groups = avail;
            int64_t vals = (int64_t)num_groups * 8;
            if (vals > room) vals = room;
            if (vals <= 56 && avail >= 8) {
                uint64_t v = carquet_read_u64_le(input + pos) & (((uint64_t)1 << vals) - 1);
                carquet_write_u64_le(d, carquet_read_u64_le(d) | (v << sh));
            } else {
                bitmap_append(bitmap, count, input + pos, vals);
            }
            pos += num_groups;
            count += vals;
        }
    }

    if (non_null_count) {
        /* Bits past `count` are still zero. */
        size_t used = (size_t)((count + 7) / 8);
        size_t b = 0;
        int64_t nn = 0;
        for (; b + 8 <= used; b += 8) nn += carquet_popcount64(carquet_read_u64_le(bitmap + b));
        for (; b < used; b++) nn += carquet_popcount32(bitmap[b]);
        *non_null_count = nn;
    }
    return count;
}

/* Stepwise bulk encoders, the fallback of the batch encoders below. Output
 * is exactly that of one carquet_rle_encoder_put() per value: a value equal
 * to the previous one only extends the current run (put's whole effect is
 * repeat_count++), so whole stretches of equal values are measured with word
 * compares and added at once, and every other value goes through the same
 * step as put. */
static carquet_status_t encode_all_stepwise(
    const uint32_t* input,
    int64_t count,
    int bit_width,
    carquet_buffer_t* output) {

    carquet_rle_encoder_t enc;
    carquet_rle_encoder_init(&enc, output, bit_width);

    int64_t i = 0;
    while (i < count) {
        uint32_t v = input[i];
        if (enc.has_prev && v == enc.prev_value) {
            int64_t j = i + 1;
            uint64_t pattern = ((uint64_t)v << 32) | v;
            while (j + 2 <= count) {
                uint64_t w;
                memcpy(&w, input + j, sizeof(w));
                if (w != pattern) break;
                j += 2;
            }
            while (j < count && input[j] == v) j++;
            enc.repeat_count += j - i;
            i = j;
            continue;
        }
        carquet_status_t status = rle_put(&enc, v);
        if (status != CARQUET_OK) return status;
        i++;
    }

    return carquet_rle_encoder_flush(&enc);
}

static carquet_status_t encode_levels_stepwise(
    const int16_t* input,
    int64_t count,
    int bit_width,
    carquet_buffer_t* output) {

    carquet_rle_encoder_t enc;
    carquet_rle_encoder_init(&enc, output, bit_width);

    int64_t i = 0;
    while (i < count) {
        uint32_t v = (uint32_t)input[i];
        if (enc.has_prev && v == enc.prev_value) {
            int64_t j = i + 1;
            uint64_t lane = (uint16_t)input[i];
            uint64_t pattern = lane | (lane << 16) | (lane << 32) | (lane << 48);
            while (j + 4 <= count) {
                uint64_t w;
                memcpy(&w, input + j, sizeof(w));
                if (w != pattern) break;
                j += 4;
            }
            while (j < count && input[j] == input[i]) j++;
            enc.repeat_count += j - i;
            i = j;
            continue;
        }
        carquet_status_t status = rle_put(&enc, v);
        if (status != CARQUET_OK) return status;
        i++;
    }

    return carquet_rle_encoder_flush(&enc);
}

/* ============================================================================
 * Batch encoders
 * ============================================================================
 *
 * carquet_rle_encode_all() and carquet_rle_encode_levels() produce exactly
 * the bytes of one carquet_rle_encoder_put() per value followed by
 * carquet_rle_encoder_flush() (the serial and page-parallel writer paths and
 * tests/test_encodings_extended.c rely on it). That output depends only on
 * the maximal runs of equal values:
 *   - every maximal run of at least 8 values is an RLE run (flush_rle, split
 *     into runs of at most 2^31 - 1 values);
 *   - the values between two such runs (a "stretch": runs shorter than 8)
 *     are bit-packed from the stretch start: its complete groups of 8 as
 *     bit-packed runs of at most 63 groups (stage_full_group), then its last
 *     (length % 8) values as one RLE run per maximal sub-run
 *     (flush_bitpack_as_rle: a partial group is never zero-padded).
 * So these encoders locate the long runs directly and emit each stretch in
 * bulk. The find_long_* helpers return the first q >= i that starts 8 equal
 * values; when i starts a run (i == 0, or i ends the previous long run), q
 * starts a maximal run too, because an earlier start in [i, q) of the same
 * run would also start 8 equal values.
 *
 * Output is written straight into the buffer's spare capacity (rle_sink_t)
 * instead of one carquet_buffer_append per run.
 */

typedef struct {
    carquet_buffer_t* buf;
    uint8_t* p;
    uint8_t* end;
    carquet_status_t status;
} rle_sink_t;

/* The batch encoders write into spare capacity, which a wrapped
 * (non-owning) buffer cannot grow; those keep the stepwise path. */
static inline bool rle_batch_buffer_ok(const carquet_buffer_t* buf) {
    return buf->owns_data || buf->data == NULL;
}

static bool rle_sink_grow(rle_sink_t* s, size_t n) {
    size_t used = (size_t)(s->p - s->buf->data);
    s->buf->size = used;
    if (s->status != CARQUET_OK) return false;
    if (n > SIZE_MAX - used) {
        s->status = CARQUET_ERROR_OUT_OF_MEMORY;
        return false;
    }
    carquet_status_t st = carquet_buffer_reserve(s->buf, used + n);
    if (st != CARQUET_OK) {
        s->status = st;
        return false;
    }
    s->p = s->buf->data + used;
    s->end = s->buf->data + s->buf->capacity;
    return true;
}

static bool rle_sink_init(rle_sink_t* s, carquet_buffer_t* buf) {
    s->buf = buf;
    s->status = CARQUET_OK;
    size_t size = buf->size;
    if (size > SIZE_MAX - 64) {
        s->status = CARQUET_ERROR_OUT_OF_MEMORY;
        return false;
    }
    carquet_status_t st = carquet_buffer_reserve(buf, size + 64);
    if (st != CARQUET_OK) {
        s->status = st;
        return false;
    }
    s->p = buf->data + size;
    s->end = buf->data + buf->capacity;
    return true;
}

/* Make room for n more bytes at s->p. */
static inline bool rle_sink_need(rle_sink_t* s, size_t n) {
    return (size_t)(s->end - s->p) >= n || rle_sink_grow(s, n);
}

static carquet_status_t rle_sink_finish(rle_sink_t* s) {
    s->buf->size = (size_t)(s->p - s->buf->data);
    return s->status;
}

static inline int rle_ctz64(uint64_t v) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_ctzll(v);
#else
    int n = 0;
    while (!(v & 1)) { v >>= 1; n++; }
    return n;
#endif
}

/* One RLE run of count copies of value (flush_rle). */
static inline void rle_sink_run(rle_sink_t* s, uint32_t value, int64_t count, int vb) {
    while (count > 0) {
        int64_t run = count > 0x7FFFFFFF ? 0x7FFFFFFF : count;
        if (!rle_sink_need(s, 10)) return;
        uint8_t* p = s->p;
        uint32_t header = (uint32_t)run << 1;
        while (header >= 0x80) {
            *p++ = (uint8_t)((header & 0x7F) | 0x80);
            header >>= 7;
        }
        *p++ = (uint8_t)header;
        write_value_bytes(p, value, vb);
        s->p = p + vb;
        count -= run;
    }
}

/* ---- uint32 values (dictionary indices) ---- */

/* First q >= i starting 8 equal values (see above), or n; *run_end gets the
 * end of that maximal run. A run of 8 or more covers a whole aligned block
 * of 4, so only all-equal aligned blocks are examined further; blocks are
 * visited in order, so the first long run found is the first one. */
static int64_t find_long_u32(const uint32_t* x, int64_t i, int64_t n, int64_t* run_end) {
    int64_t b = (i + 3) & ~(int64_t)3;
    for (; b + 4 <= n; b += 4) {
        uint32_t v = x[b];
        if (((x[b + 1] ^ v) | (x[b + 2] ^ v) | (x[b + 3] ^ v)) != 0) continue;
        int64_t p = b;
        while (p > i && x[p - 1] == v) p--;
        int64_t e = b + 4;
        uint64_t pattern = ((uint64_t)v << 32) | v;
        while (e + 2 <= n) {
            uint64_t w;
            memcpy(&w, x + e, sizeof(w));
            if (w != pattern) break;
            e += 2;
        }
        while (e < n && x[e] == v) e++;
        if (e - p >= 8) {
            *run_end = e;
            return p;
        }
    }
    return n;
}

/* Pack one group of 8 values: the bytes of carquet_bitpack8_32. For widths
 * up to 8 the group is assembled in one 64-bit word and stored whole (8
 * bytes, so the caller leaves 8 bytes of slack). */
static inline void rle_pack8_u32(const uint32_t* v, int bw, uint8_t* out) {
    if (bw <= 8) {
        uint32_t mask = (1u << bw) - 1u;
        uint64_t acc = 0;
        for (int k = 0; k < 8; k++) {
            acc |= (uint64_t)(v[k] & mask) << (k * bw);
        }
        carquet_write_u64_le(out, acc);
    } else {
        carquet_bitpack8_32(v, bw, out);
    }
}

/* The last t (< 8) values of a stretch: one RLE run per maximal sub-run. */
static void emit_tail_u32(rle_sink_t* s, const uint32_t* x, int t, int vb) {
    if (t == 0 || !rle_sink_need(s, 64)) return;
    uint8_t* p = s->p;
    int j = 0;
    while (j < t) {
        uint32_t v = x[j];
        int r = 1;
        while (j + r < t && x[j + r] == v) r++;
        *p++ = (uint8_t)(r << 1);
        write_value_bytes(p, v, vb);
        p += vb;
        j += r;
    }
    s->p = p;
}

static void emit_stretch_u32(rle_sink_t* s, const uint32_t* x, int64_t len, int bw, int vb) {
    int64_t groups = len >> 3;
    while (groups > 0) {
        int g = groups > 63 ? 63 : (int)groups;
        if (!rle_sink_need(s, 1 + (size_t)g * (size_t)bw + 8)) return;
        uint8_t* p = s->p;
        *p++ = (uint8_t)((g << 1) | 1);
        if (bw > 0) {
            for (int k = 0; k < g; k++) {
                rle_pack8_u32(x, bw, p);
                x += 8;
                p += bw;
            }
        } else {
            x += 8 * g;
        }
        s->p = p;
        groups -= g;
    }
    emit_tail_u32(s, x, (int)(len & 7), vb);
}

static carquet_status_t encode_all_batch(
    const uint32_t* x, int64_t n, int bw, carquet_buffer_t* output) {
    rle_sink_t s;
    if (!rle_sink_init(&s, output)) return s.status;
    int vb = rle_value_bytes(bw);
    int64_t start = 0;
    while (s.status == CARQUET_OK) {
        int64_t e = n;
        int64_t q = find_long_u32(x, start, n, &e);
        if (q >= n) break;
        emit_stretch_u32(&s, x + start, q - start, bw, vb);
        rle_sink_run(&s, x[q], e - q, vb);
        start = e;
    }
    emit_stretch_u32(&s, x + start, n - start, bw, vb);
    return rle_sink_finish(&s);
}

/* ---- 1-bit levels as a bit array ----
 * Level i is bit i % 64 of w[i / 64]; bits past n are zero and the array
 * has one zero word past the last one holding levels, so rle_bits_at(p)
 * may be called for any p < n. A bit-packed group of eight 1-bit values is
 * exactly the byte of the bit array that holds them, so groups are copied
 * out 64 bits at a time. */

static inline uint64_t rle_bits_at(const uint64_t* w, int64_t pos) {
    int64_t k = pos >> 6;
    unsigned sh = (unsigned)(pos & 63);
    uint64_t v = w[k] >> sh;
    if (sh) v |= w[k + 1] << (64 - sh);
    return v;
}

/* First q >= i starting 8 equal bits, or n. e has bit k set when bits
 * i + k and i + k + 1 are equal; 8 equal bits from k are 7 such bits in a
 * row, which m keeps (m = e[k..k+6] all set). */
static int64_t find_long_bits(const uint64_t* w, int64_t i, int64_t n) {
    while (i + 8 <= n) {
        uint64_t e = ~(rle_bits_at(w, i) ^ rle_bits_at(w, i + 1));
        int64_t lim = n - 1 - i;            /* comparable pairs, >= 7 */
        if (lim < 64) {
            e &= (1ULL << lim) - 1;
        } else {
            lim = 64;
        }
        uint64_t m = e & (e >> 1);
        m &= m >> 2;
        m &= m >> 3;
        if (m) return i + rle_ctz64(m);
        i += lim - 6;
    }
    return n;
}

/* First p' >= p whose bit differs from v, or n. */
static int64_t run_end_bits(const uint64_t* w, int64_t p, int64_t n, uint32_t v) {
    uint64_t flip = v ? ~0ULL : 0;
    while (p < n) {
        uint64_t d = rle_bits_at(w, p) ^ flip;
        if (d) {
            p += rle_ctz64(d);
            return p < n ? p : n;
        }
        p += 64;
    }
    return n;
}

static void emit_stretch_bits(rle_sink_t* s, const uint64_t* w, int64_t pos, int64_t end) {
    int64_t groups = (end - pos) >> 3;
    while (groups > 0) {
        int g = groups > 63 ? 63 : (int)groups;
        if (!rle_sink_need(s, (size_t)g + 9)) return;
        uint8_t* p = s->p;
        *p++ = (uint8_t)((g << 1) | 1);
        for (int k = 0; k < g; k += 8) {
            carquet_write_u64_le(p + k, rle_bits_at(w, pos + 8 * (int64_t)k));
        }
        s->p = p + g;
        pos += 8 * (int64_t)g;
        groups -= g;
    }
    int t = (int)(end - pos);
    if (t == 0 || !rle_sink_need(s, 16)) return;
    /* Tail runs: each run ends at the first bit differing from its value,
     * or at the sentinel bit t. */
    uint32_t bits = (uint32_t)rle_bits_at(w, pos) & ((1u << t) - 1u);
    uint32_t sentinel = 1u << t;
    uint8_t* p = s->p;
    int j = 0;
    while (j < t) {
        uint32_t v = (bits >> j) & 1u;
        uint32_t diff = ((bits ^ (0u - v)) | sentinel) >> j;
        int r = carquet_ctz32(diff);
        p[0] = (uint8_t)(r << 1);
        p[1] = (uint8_t)v;
        p += 2;
        j += r;
    }
    s->p = p;
}

static carquet_status_t encode_bits_batch(
    const uint64_t* w, int64_t n, carquet_buffer_t* output) {
    rle_sink_t s;
    if (!rle_sink_init(&s, output)) return s.status;
    int64_t start = 0;
    while (s.status == CARQUET_OK) {
        int64_t q = find_long_bits(w, start, n);
        if (q >= n) break;
        uint32_t v = (uint32_t)(w[q >> 6] >> (q & 63)) & 1u;
        int64_t e = run_end_bits(w, q + 8, n, v);
        emit_stretch_bits(&s, w, start, q);
        rle_sink_run(&s, v, e - q, 1);
        start = e;
    }
    emit_stretch_bits(&s, w, start, n);
    return rle_sink_finish(&s);
}

carquet_status_t carquet_rle_encode_all(
    const uint32_t* input,
    int64_t count,
    int bit_width,
    carquet_buffer_t* output) {

    if (count > 0 && bit_width >= 0 && bit_width <= 32 &&
        rle_batch_buffer_ok(output)) {
        return encode_all_batch(input, count, bit_width, output);
    }
    return encode_all_stepwise(input, count, bit_width, output);
}

/* Levels of a 1-bit stream (max level 1: flat nullable definition levels,
 * single-level list repetition levels) go through the bit array; a stream
 * holding any other value, and wider levels, take the stepwise path. */
enum { RLE_STACK_BIT_WORDS = 66 };

carquet_status_t carquet_rle_encode_levels(
    const int16_t* input,
    int64_t count,
    int bit_width,
    carquet_buffer_t* output) {
    return carquet_rle_encode_levels_counted(input, count, bit_width, output, NULL);
}

carquet_status_t carquet_rle_encode_levels_counted(
    const int16_t* input,
    int64_t count,
    int bit_width,
    carquet_buffer_t* output,
    int64_t* num_ones) {

    if (num_ones) *num_ones = -1;
    if (bit_width == 1 && count > 0 && rle_batch_buffer_ok(output)) {
        uint64_t stack_words[RLE_STACK_BIT_WORDS];
        size_t nwords = (size_t)((count + 63) / 64) + 1;
        uint64_t* w = stack_words;
        if (nwords > RLE_STACK_BIT_WORDS) {
            w = (nwords <= SIZE_MAX / sizeof(uint64_t))
                ? (uint64_t*)carquet_mem_malloc(nwords * sizeof(uint64_t))
                : NULL;
        }
        if (w) {
            bool valid = carquet_dispatch_levels_to_bits(input, count, w);
            w[nwords - 1] = 0;
            carquet_status_t status = CARQUET_OK;
            if (valid) {
                status = encode_bits_batch(w, count, output);
                if (num_ones && status == CARQUET_OK) {
                    int64_t ones = 0;
                    for (size_t k = 0; k + 1 < nwords; k++) ones += carquet_popcount64(w[k]);
                    *num_ones = ones;
                }
            }
            if (w != stack_words) carquet_mem_free(w);
            if (valid) return status;
        }
    }
    return encode_levels_stepwise(input, count, bit_width, output);
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
