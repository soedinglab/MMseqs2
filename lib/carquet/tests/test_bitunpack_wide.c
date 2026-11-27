/**
 * @file test_bitunpack_wide.c
 * @brief Correctness oracle for the wide (SIMD) bit-unpack fast path.
 *
 * Validates carquet_bitunpack_32() (which feeds DELTA decode) and the RLE
 * decoders (which feed dictionary indices and definition/repetition levels)
 * against an INDEPENDENT, trivially-correct scalar bit reader. The whole
 * point is to not reuse any library bit logic in the oracle, so a bug in a
 * wide SIMD kernel cannot hide behind a shared helper.
 */

#include "core/bitpack.h"
#include "encoding/rle.h"
#include "core/buffer.h"
#include <carquet/carquet.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TEST_PASS(name) printf("[PASS] %s\n", name)
#define TEST_FAIL(name, msg) do { printf("[FAIL] %s: %s\n", name, msg); return 1; } while(0)

/* Independent oracle: read `count` values of `bit_width` bits, LSB-first,
 * little-endian byte order — the documented packing contract. No library
 * code involved. */
static void oracle_unpack(const uint8_t* in, size_t count, int bit_width,
                          uint32_t* out) {
    size_t bitpos = 0;
    for (size_t i = 0; i < count; i++) {
        uint32_t v = 0;
        for (int b = 0; b < bit_width; b++) {
            size_t bp = bitpos + (size_t)b;
            uint32_t bit = (in[bp >> 3] >> (bp & 7)) & 1u;
            v |= bit << b;
        }
        out[i] = v;
        bitpos += (size_t)bit_width;
    }
}

static uint32_t rng_state = 0x12345678u;
static uint32_t rng(void) {
    /* xorshift32 — deterministic, defined unsigned arithmetic */
    uint32_t x = rng_state;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    rng_state = x;
    return x;
}

/* Exhaustive-ish: every bit width, many counts that straddle wide-chunk
 * (16/32-value) and 8-value and tail boundaries, random packed input. */
static int test_bitunpack_32_vs_oracle(void) {
    const char* name = "bitunpack_32_vs_oracle";
    static const size_t counts[] = {
        1, 2, 7, 8, 9, 15, 16, 17, 31, 32, 33, 40, 63, 64, 65,
        96, 97, 127, 128, 129, 200, 255, 256, 257, 512, 1000, 4096
    };
    uint32_t got[4200];
    uint32_t exp[4200];
    uint8_t packed[4200 * 4 + 64];

    for (int bw = 1; bw <= 32; bw++) {
        for (size_t ci = 0; ci < sizeof(counts) / sizeof(counts[0]); ci++) {
            size_t count = counts[ci];
            size_t nbytes = carquet_packed_size(count, bw);
            for (size_t k = 0; k < nbytes + 16; k++) packed[k] = (uint8_t)rng();

            memset(got, 0xCC, sizeof(uint32_t) * count);
            size_t consumed = carquet_bitunpack_32(packed, count, bw, got);
            oracle_unpack(packed, count, bw, exp);

            if (consumed != nbytes) {
                printf("  bw=%d count=%zu consumed=%zu want=%zu\n",
                       bw, count, consumed, nbytes);
                TEST_FAIL(name, "wrong consumed byte count");
            }
            for (size_t i = 0; i < count; i++) {
                if (got[i] != exp[i]) {
                    printf("  bw=%d count=%zu i=%zu got=%u exp=%u\n",
                           bw, count, i, got[i], exp[i]);
                    TEST_FAIL(name, "value mismatch vs oracle");
                }
            }
        }
    }
    TEST_PASS(name);
    return 0;
}

/* RLE bit-packed path (levels + dictionary indices) vs oracle-built values. */
static int test_rle_bitpacked_vs_oracle(void) {
    const char* name = "rle_bitpacked_vs_oracle";

    for (int bw = 1; bw <= 24; bw++) {
        uint32_t mask = (bw >= 32) ? ~0u : ((1u << bw) - 1u);
        for (int trial = 0; trial < 12; trial++) {
            int64_t n = 1 + (int64_t)(rng() % 5000);
            uint32_t* vals = malloc(sizeof(uint32_t) * (size_t)n);
            for (int64_t i = 0; i < n; i++) {
                /* mix long repeats (RLE runs) and noise (bit-packed runs) */
                if ((rng() & 3) == 0) {
                    uint32_t r = rng() & mask;
                    int64_t runlen = 1 + (rng() % 40);
                    for (int64_t j = 0; j < runlen && i < n; j++, i++) vals[i] = r;
                    i--;
                } else {
                    vals[i] = rng() & mask;
                }
            }

            carquet_buffer_t buf;
            carquet_buffer_init(&buf);
            carquet_rle_encoder_t enc;
            carquet_rle_encoder_init(&enc, &buf, bw);
            for (int64_t i = 0; i < n; i++) carquet_rle_encoder_put(&enc, vals[i]);
            carquet_rle_encoder_flush(&enc);

            const uint8_t* data = carquet_buffer_data(&buf);
            size_t size = carquet_buffer_size(&buf);

            /* Path A: dictionary-index decode (carquet_rle_decode_all). */
            uint32_t* out = malloc(sizeof(uint32_t) * (size_t)n);
            int64_t dec = carquet_rle_decode_all(data, size, bw, out, n);
            if (dec != n) { free(vals); free(out); carquet_buffer_destroy(&buf);
                TEST_FAIL(name, "decode_all wrong count"); }
            for (int64_t i = 0; i < n; i++) if (out[i] != vals[i]) {
                printf("  all bw=%d i=%lld got=%u exp=%u\n",
                       bw, (long long)i, out[i], vals[i]);
                free(vals); free(out); carquet_buffer_destroy(&buf);
                TEST_FAIL(name, "decode_all value mismatch");
            }

            /* Path B: level decode (carquet_rle_decode_levels), bw<=16. */
            if (bw <= 16) {
                int16_t* lv = malloc(sizeof(int16_t) * (size_t)n);
                int64_t ld = carquet_rle_decode_levels(data, size, bw, lv, n);
                if (ld != n) { free(vals); free(out); free(lv);
                    carquet_buffer_destroy(&buf);
                    TEST_FAIL(name, "decode_levels wrong count"); }
                for (int64_t i = 0; i < n; i++) {
                    if ((uint16_t)lv[i] != (uint16_t)vals[i]) {
                        printf("  lvl bw=%d i=%lld got=%d exp=%u\n",
                               bw, (long long)i, lv[i], vals[i]);
                        free(vals); free(out); free(lv);
                        carquet_buffer_destroy(&buf);
                        TEST_FAIL(name, "decode_levels value mismatch");
                    }
                }
                free(lv);
            }

            /* Path C: partial cap (exercise max_values truncation). */
            int64_t cap = n / 3 + 1;
            int64_t dc = carquet_rle_decode_all(data, size, bw, out, cap);
            if (dc != cap) { free(vals); free(out); carquet_buffer_destroy(&buf);
                TEST_FAIL(name, "decode_all cap wrong count"); }
            for (int64_t i = 0; i < cap; i++) if (out[i] != vals[i]) {
                free(vals); free(out); carquet_buffer_destroy(&buf);
                TEST_FAIL(name, "decode_all cap value mismatch");
            }

            free(vals); free(out);
            carquet_buffer_destroy(&buf);
        }
    }
    TEST_PASS(name);
    return 0;
}

int main(void) {
    int rc = 0;
    printf("=== Wide bit-unpack correctness ===\n\n");
    rc |= test_bitunpack_32_vs_oracle();
    rc |= test_rle_bitpacked_vs_oracle();
    if (rc == 0) printf("\nAll wide bit-unpack tests passed.\n");
    return rc;
}
